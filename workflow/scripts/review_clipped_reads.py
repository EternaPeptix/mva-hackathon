#!/usr/bin/env python3
"""Fit full seed queries; remap canonical seed and mate ends to the whole reference.

Distances and remapped alignments are descriptive sensitivity results, not calls.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time
import pysam

from clipped_read_helpers import fit_full_query
from joint_haplotype_helpers import enumerate_haplotypes
from review_joint_haplotypes import collect_native_context, window_edits


def require(value, message):
    if not value:
        raise ValueError(message)


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def pin(path, expected=None):
    path = Path(path).resolve()
    before = path.stat()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    after = path.stat()
    require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), 'input changed while hashing')
    result = dict(path=str(path), bytes=after.st_size, mtime_ns=after.st_mtime_ns, sha256=digest.hexdigest())
    if expected:
        require(result['sha256'] == expected['sha256'], 'input SHA differs: ' + str(path))
        require('bytes' not in expected or result['bytes'] == expected['bytes'], 'input bytes differ')
    return result


def identity(read):
    require(read.has_tag('RG') and read.get_tag('RG') and read.query_name, 'missing original fragment identity')
    require(not (read.is_read1 and read.is_read2), 'read has both end-number flags')
    number = 1 if read.is_read1 else 2 if read.is_read2 else 0
    require(not read.is_paired or number in (1, 2), 'paired read lacks end number')
    return read.get_tag('RG'), read.query_name, number


def describe(read, identifier, roles):
    seq, qual = read.query_sequence, read.query_qualities
    forward, forward_quality = read.get_forward_sequence(), read.get_forward_qualities()
    clips = []
    query_position = 0
    for index, (operation, length) in enumerate(read.cigartuples or []):
        if operation == 4:
            clips.append(dict(side='left' if index == 0 or all(op == 5 for op, _ in read.cigartuples[:index]) else 'right',
                query_start0=query_position, query_end0=query_position + length,
                sequence=seq[query_position:query_position + length] if seq else None,
                raw_base_qualities=list(qual[query_position:query_position + length]) if qual is not None else None))
        if operation in (0, 1, 4, 7, 8):
            query_position += length
    key = identity(read)
    reasons = []
    if read.flag & 0xF04:
        reasons.append('original_flag_0xF04')
    if read.mapping_quality == 255:
        reasons.append('original_mapping_quality_unavailable_255')
    if not seq or re.fullmatch('[ACGT]+', seq) is None:
        reasons.append('query_missing_or_nonconcrete_ACGT')
    fastq_reason = ('query_missing' if not forward else 'quality_missing_or_not_FASTQ_representable'
        if forward_quality is None or len(forward_quality) != len(forward)
        or any(q < 0 or q > 93 for q in forward_quality) else None)
    return dict(id=identifier, RG=key[0], QNAME=key[1], end_number=key[2], roles=sorted(roles),
        original_SAM=read.to_string(), flag=read.flag, original_MQ=read.mapping_quality,
        chrom=read.reference_name, start0=read.reference_start, end0=read.reference_end,
        CIGAR=read.cigarstring, strand='reverse' if read.is_reverse else 'forward',
        mate_chrom=read.next_reference_name, mate_start0=read.next_reference_start, TLEN=read.template_length,
        original_alignment_query=seq, raw_base_qualities=list(qual) if qual is not None else None,
        forward_query=forward, forward_base_qualities=list(forward_quality) if forward_quality is not None else None,
        terminal_clips=clips, local_fit_ineligibility=reasons,
        remap_ineligibility=fastq_reason, raw_quality_summary=dict(
            minimum=min(qual) if qual is not None and len(qual) else None,
            below20=sum(q < 20 for q in qual) if qual is not None else None,
            below30=sum(q < 30 for q in qual) if qual is not None else None))


def recover_ends(config, pins):
    seed = config['seed_window']
    with pysam.AlignmentFile(pins['alignment']['path'], reference_filename=pins['reference']['path'],
            index_filename=pins['alignment_index']['path']) as original, \
            pysam.AlignmentFile(pins['seed_bam']['path'], 'rb', index_filename=pins['seed_bam_index']['path']) as subset:
        full = list(original.fetch(seed['chrom'], seed['start0'], seed['end0']))
        retained = list(subset.fetch(seed['chrom'], seed['start0'], seed['end0']))
        require([r.to_string() for r in full] == [r.to_string() for r in retained], 'original and pinned seed records differ')
        require(original.header.to_dict()['SQ'] == subset.header.to_dict()['SQ'] and
                original.header.to_dict().get('RG') == subset.header.to_dict().get('RG'), 'seed headers differ')
        groups = {r['ID']: r for r in original.header.to_dict().get('RG', [])}
        require(len(groups) == config['expected_read_groups'] and
                {r.get('SM') for r in groups.values()} == {config['sample']}, 'unexpected sample/read groups')
        require(len(full) <= config['maximum_seed_records'], 'seed record bound exceeded')
        canonical, roles, seed_keys = {}, defaultdict(set), []
        raw = []
        for read in full:
            require(read.has_tag('RG') and read.get_tag('RG') in groups, 'undeclared seed read group')
            raw.append(read.to_string())
            if read.is_secondary or read.is_supplementary:
                continue
            key = identity(read)
            require(key not in canonical, 'repeated original primary end identity')
            canonical[key] = read
            roles[key].add('seed')
            seed_keys.append(key)
        retrieval = []
        for key in seed_keys:
            read = canonical[key]
            entry = dict(seed_identity=key)
            retrieval.append(entry)
            if read.flag & 0xF04:
                entry['status'] = 'seed_excluded_original_flags'
                continue
            if not read.is_paired or read.mate_is_unmapped:
                entry['status'] = 'unpaired_or_mate_reported_unmapped'
                continue
            wanted = key[:2] + (3 - key[2],)
            matches = [mate for mate in original.fetch(read.next_reference_name,
                read.next_reference_start, read.next_reference_start + 1)
                if not (mate.is_secondary or mate.is_supplementary) and mate.has_tag('RG')
                and (mate.get_tag('RG'), mate.query_name, 1 if mate.is_read1 else 2 if mate.is_read2 else 0) == wanted
                and mate.reference_start == read.next_reference_start]
            if len(matches) != 1:
                entry.update(status='mate_lookup_missing_or_ambiguous', candidate_records=len(matches))
                continue
            mate = matches[0]
            reciprocal = (mate.reference_id == read.next_reference_id and
                mate.reference_start == read.next_reference_start and mate.next_reference_id == read.reference_id and
                mate.next_reference_start == read.reference_start and not mate.mate_is_unmapped)
            if not reciprocal:
                entry['status'] = 'mate_coordinates_not_reciprocal'
                continue
            require(identity(mate) == wanted, 'mate identity changed')
            require(wanted not in canonical or canonical[wanted].to_string() == mate.to_string(), 'mate/seed payload disagrees')
            canonical[wanted] = mate
            roles[wanted].add('mate')
            entry.update(status='retrieved_and_reciprocal', mate_identity=wanted)
        require(len(canonical) <= config['maximum_canonical_ends'], 'canonical end bound exceeded')
    ids = {key: 'Q%05d' % index for index, key in enumerate(sorted(canonical), 1)}
    rows = {ids[key]: describe(read, ids[key], roles[key]) for key, read in canonical.items()}
    for entry in retrieval:
        entry['seed_id'] = ids[tuple(entry.pop('seed_identity'))]
        if 'mate_identity' in entry:
            entry['mate_id'] = ids[tuple(entry.pop('mate_identity'))]
    return rows, retrieval, dict(raw_seed_records=len(raw), primary_seed_ends=len(seed_keys),
        canonical_ends=len(rows), unique_fragments=len({key[:2] for key in canonical}),
        decoded_seed_record_payload_sha256=hashlib.sha256(('\n'.join(raw) + '\n').encode()).hexdigest(),
        raw_seed_SAM=raw)


def fastq(row, identifier):
    return '@' + identifier + '\n' + row['forward_query'] + '\n+\n' + ''.join(chr(q + 33) for q in row['forward_base_qualities']) + '\n'


def remap_arm(pins, config, rows, output, command, paired=False):
    mapping, pairs = {}, []
    if paired:
        fragments = defaultdict(dict)
        for identifier, row in rows.items():
            fragments[(row['RG'], row['QNAME'])][row['end_number']] = identifier
        for ends in fragments.values():
            if set(ends) != {1, 2}:
                continue
            first, second = (rows[ends[n]] for n in (1, 2))
            if any(row['remap_ineligibility'] or row['flag'] & 0xF04 for row in (first, second)):
                continue
            if any(not row['flag'] & 1 or row['flag'] & 8 for row in (first, second)):
                continue
            if not (first['mate_chrom'] == second['chrom'] and second['mate_chrom'] == first['chrom'] and
                    first['mate_start0'] == second['start0'] and second['mate_start0'] == first['start0']):
                continue
            pair_id = 'F%05d' % (len(pairs) + 1)
            pairs.append((pair_id, first, second))
            mapping[(pair_id, 1)], mapping[(pair_id, 2)] = first['id'], second['id']
        label = 'paired-whole-reference'
        paths = [output / (label + '.R%d.fastq' % n) for n in (1, 2)]
        for number, path in enumerate(paths, 1):
            path.write_text(''.join(fastq(item[number], item[0]) for item in pairs))
    else:
        label = 'single-whole-reference'
        eligible = [row for row in rows.values() if not row['remap_ineligibility']]
        mapping = {(row['id'], 0): row['id'] for row in eligible}
        paths = [output / (label + '.fastq')]
        paths[0].write_text(''.join(fastq(row, row['id']) for row in eligible))
    if not mapping:
        return dict(status='no_eligible_queries', paired=paired, reported_alignments=[])
    sam = output / (label + '.sam')
    command([pins['bwa']['path'], 'mem', '-a', '-Y', '-t', str(config['bwa_threads']),
             pins['reference']['path'], *map(str, paths)], label, stdout_path=sam)
    observations, primary = [], Counter()
    with pysam.AlignmentFile(str(sam), 'r') as reader:
        for read in reader:
            number = 1 if read.is_read1 else 2 if read.is_read2 else 0
            require((read.query_name, number) in mapping, 'unexpected remapped query identity')
            identifier = mapping[(read.query_name, number)]
            row = rows[identifier]
            if read.query_sequence is not None:
                require(read.get_forward_sequence() == row['forward_query'], 'remap query sequence changed')
                require(read.get_forward_qualities() is not None and
                    list(read.get_forward_qualities()) == row['forward_base_qualities'], 'remap query qualities changed')
            is_primary = not (read.is_secondary or read.is_supplementary)
            if is_primary:
                require(read.query_sequence is not None, 'primary remap omitted query')
                primary[identifier] += 1
            observations.append(dict(canonical_end_id=identifier, export_query_name=read.query_name,
                flag=read.flag, is_primary=is_primary, chrom=read.reference_name,
                start0=read.reference_start, end0=read.reference_end, MQ=read.mapping_quality,
                CIGAR=read.cigarstring, strand='reverse' if read.is_reverse else 'forward',
                aligned_query_bases=read.query_alignment_length, query_length=read.query_length,
                AS=read.get_tag('AS') if read.has_tag('AS') else None,
                NM=read.get_tag('NM') if read.has_tag('NM') else None,
                mate_chrom=read.next_reference_name, mate_start0=read.next_reference_start,
                TLEN=read.template_length, original_SAM=read.to_string(),
                query_payload_verified=read.query_sequence is not None))
    require(dict(primary) == {identifier: 1 for identifier in mapping.values()}, 'missing/repeated primary remap')
    return dict(status='complete', paired=paired, exported_ends=len(mapping), exported_pairs=len(pairs) if paired else None,
        original_flag_counts=dict(Counter(rows[identifier]['flag'] for identifier in mapping.values())),
        original_flag_excluded_ends_exported=sum(bool(rows[identifier]['flag'] & 0xF04) for identifier in mapping.values()),
        selection_scope='Reciprocal flag-eligible original pairs with valid FASTQ payloads' if paired else 'All canonical ends with valid FASTQ payloads, including explicitly flagged contextual observations',
        export_identity=[dict(query=name, end_number=n, canonical_end_id=identifier) for (name, n), identifier in mapping.items()],
        reported_alignments=observations, primary_query_payloads_verified=True,
        limitation='BWA reported alternatives are not exhaustive. Paired remapping estimates insert-size behavior from this selected subset; original mapping and mate context remain authoritative observations.')


def run(args):
    root, output = args.root.resolve(), args.output.resolve()
    require(not output.exists() and not output.is_relative_to(Path(__file__).resolve().parents[2]), 'new output outside source required')
    output.mkdir(parents=True)
    os.nice(10)
    start = time.monotonic()
    state = dict(status='running', pid=os.getpid(), started_at=stamp(), stage='verify_inputs')
    save(output / 'status.json', state)
    try:
        config = json.loads(args.config.read_text())
        pins = {'config': pin(args.config)}
        for name, expected in config['inputs'].items():
            pins[name] = pin(root / expected['path'], expected)
        for name, digest in config['source_sha256'].items():
            pins['source/' + name] = pin(Path(__file__).with_name(name), {'sha256': digest})
        require(pins['reference_index']['path'] == pins['reference']['path'] + '.fai', 'reference index path differs')
        index = json.loads(Path(pins['bwa_index_record']['path']).read_text())
        require(index['reference_sha256'] == pins['reference']['sha256'], 'BWA reference index binding differs')
        for suffix, digest in index['index_sha256'].items():
            require(pins['bwa_index' + suffix]['sha256'] == digest and
                pins['bwa_index' + suffix]['path'] == pins['reference']['path'] + suffix, 'BWA index differs')
        require(config['bwa_threads'] in (1, 2) and config['maximum_hypothesis_states'] <= 256, 'unbounded compute configuration')
        rows, mates, recovery = recover_ends(config, pins)
        save(output / 'original-ends.json', dict(ends=rows, mate_retrieval=mates, recovery=recovery))
        contexts, census = collect_native_context(config, pins, config['targets'])
        save(output / 'native-context.json', dict(targets=contexts, source_census=census))
        commands = []

        def command(argv, label, stdout_path=None):
            state['stage'] = label
            save(output / 'status.json', state)
            entry = dict(argv=list(map(str, argv)), label=label, started_at=stamp())
            commands.append(entry)
            save(output / 'commands.json', commands)
            with (output / (label + '.log')).open('wb') as log:
                if stdout_path:
                    with stdout_path.open('wb') as stdout:
                        process = subprocess.run(entry['argv'], stdout=stdout, stderr=log)
                else:
                    process = subprocess.run(entry['argv'], stdout=log, stderr=subprocess.STDOUT)
            entry.update(exit_code=process.returncode, completed_at=stamp())
            save(output / 'commands.json', commands)
            require(process.returncode == 0, 'native command failed: ' + label)

        panels, fits, total_cells = [], {}, 0
        callers = sorted({s['caller'] for s in config['caller_sources']})
        with pysam.FastaFile(pins['reference']['path']) as fasta:
            for target in config['targets']:
                for flank in config['window_flanks_bp']:
                    left, right = target['core_start0'] - flank, target['core_end0'] + flank
                    require(0 <= left < right <= fasta.get_reference_length(target['chrom']), 'fit window outside reference')
                    ref = fasta.fetch(target['chrom'], left, right).upper()
                    edits, omitted = window_edits(contexts[target['id']], left, right)
                    family = enumerate_haplotypes(target['chrom'], left, ref, edits, callers=callers,
                        max_states=config['maximum_hypothesis_states'], max_haplotype_bases=config['maximum_haplotype_bases'])
                    panel_id = '%s_F%d' % (target['id'], flank)
                    panel = dict(id=panel_id, chrom=target['chrom'], start0=left, end0=right,
                        reference_sequence=ref, native_edits=edits, omitted_native_alternatives=omitted,
                        enumeration={k: v for k, v in family.items() if k != 'haplotypes'}, hypotheses={}, end_comparisons=[])
                    reference_id = None
                    for i, hypothesis in enumerate(family['haplotypes']):
                        hid = panel_id + '_H%04d' % i
                        panel['hypotheses'][hid] = hypothesis
                        if hypothesis['sequence'] == ref:
                            reference_id = hid
                    require(reference_id is not None, 'missing reference hypothesis')
                    panel['reference_hypothesis_id'] = reference_id
                    for identifier, row in rows.items():
                        if 'seed' not in row['roles']:
                            continue
                        comparison = dict(canonical_end_id=identifier, status='ineligible',
                            reasons=row['local_fit_ineligibility'], fits={})
                        panel['end_comparisons'].append(comparison)
                        if row['local_fit_ineligibility']:
                            continue
                        query = row['original_alignment_query']
                        state['stage'] = 'full_query_fit/' + panel_id + '/' + identifier
                        save(output / 'status.json', state)
                        for hid, hypothesis in panel['hypotheses'].items():
                            sequence = hypothesis['sequence']
                            cache_key = hashlib.sha256((query + '\0' + sequence).encode()).hexdigest()
                            if cache_key not in fits:
                                cells = (len(query) + 1) * (len(sequence) + 1)
                                total_cells += cells
                                require(total_cells <= config['maximum_total_fit_cells'], 'total fitting work cap exceeded')
                                result = fit_full_query(query, sequence, max_cells=config['maximum_fit_cells'], max_paths=config['maximum_tracebacks'])
                                fits[cache_key] = dict(query=query, target_sequence=sequence, result=result)
                            comparison['fits'][hid] = cache_key
                        costs = {hid: fits[key]['result']['best_cost'] for hid, key in comparison['fits'].items()}
                        best = min(costs.values())
                        comparison.update(status='scored_descriptively', costs=costs, minimum_cost=best,
                            reference_cost=costs[reference_id], best_hypothesis_ids=[hid for hid, cost in costs.items() if cost == best],
                            quality_used_in_score=False, query_orientation='original_SAM_alignment_orientation')
                    panels.append(panel)
        evaluated_cells = sum(f['result']['dp_cells_evaluated'] for f in fits.values())
        save(output / 'local-fits.json', dict(panels=panels, fits=fits, budgeted_fit_cells=total_cells,
            evaluated_DP_cells=evaluated_cells,
            scope='Exact full-query unit edit distances for seed ends only against bounded per-caller sequence panels; recovered mate-only ends receive genome remapping and original-context review, not local scores. No allele-support or genotype labels. Nonconcrete/flag-excluded/MQ255 seeds remain explicit. Raw qualities are retained but do not weight distances.'))
        remaps = dict(single=remap_arm(pins, config, rows, output, command),
            paired=remap_arm(pins, config, rows, output, command, paired=True))
        save(output / 'remapping.json', remaps)
        for fingerprint in pins.values():
            now = Path(fingerprint['path']).stat()
            require((now.st_size, now.st_mtime_ns) == (fingerprint['bytes'], fingerprint['mtime_ns']), 'input changed during review')
        summary = dict(status='complete', completed_at=stamp(), elapsed_seconds=time.monotonic() - start,
            inputs=pins, recovery={k: v for k, v in recovery.items() if k != 'raw_seed_SAM'},
            local_fit_panels=len(panels), unique_fits=len(fits), budgeted_fit_cells=total_cells,
            evaluated_DP_cells=evaluated_cells,
            local_fit_scope='Seed ends only; recovered mate-only ends are retained for original context and genome remapping',
            remap_counts={k: {name: v.get(name) for name in ('status','exported_ends','exported_pairs','primary_query_payloads_verified')} for k, v in remaps.items()},
            limitations=['Post hoc local and remapping comparisons reuse original molecules and cannot be added as independent evidence.',
                'Small-variant panels do not exhaust longer insertions, complex assemblies or alternative loci. Poor fit remains unexplained.',
                'Unit edit costs ignore base qualities and are not likelihoods, genotype scores or EPCR. Full original raw qualities are retained.',
                'Traceback caps preserve exact minimum cost and optimal endpoints but may leave start positions, CIGARs and allele coverage incomplete.',
                'Fitting coordinates are offsets in each constructed sequence, not reference coordinates. Best fits can avoid the candidate bases; cost ties or advantages alone are not allele evidence.',
                'Full-query fitting is in original SAM orientation. Genome remapping considers BWA-reported orientations and split alternatives, not exhaustive locus origin.',
                'Paired remapping uses a selected subset and re-estimates insert-size statistics; it does not reproduce the full original alignment process.',
                'No new native call, parental phase, functional mechanism, diagnosis or treatment benefit is inferred.'],
            outputs={p.name: pin(p) for p in sorted(output.iterdir()) if p.is_file() and p.name not in ('summary.json','status.json')},
            runtime=dict(python=os.sys.version, pysam=pysam.__version__))
        save(output / 'summary.json', summary)
        state.update(status='complete', stage='completed', completed_at=summary['completed_at'], summary_sha256=pin(output / 'summary.json')['sha256'])
        save(output / 'status.json', state)
        return summary
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__, error=str(error))
        save(output / 'status.json', state)
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root','config','output'):
        parser.add_argument('--' + name, required=True, type=Path)
    run(parser.parse_args())
