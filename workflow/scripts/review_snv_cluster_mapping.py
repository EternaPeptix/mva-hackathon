#!/usr/bin/env python3
"""Bounded SNV-cluster read, mate, caller and full-reference mapping review.

This preserves observed sequence and missingness. Remapping the same reads is a
sensitivity analysis, not independent replication or proof of their true origin.
Detailed reads, sequences and case outputs must remain outside the source tree.
"""
import argparse
from collections import Counter, defaultdict
import datetime
import itertools
import json
import os
from pathlib import Path
import subprocess
import time

import pysam

from review_cnv_context import fingerprint, require, save, sha
from review_paralog_candidate_pair import (record_details, remap_details,
                                          compare_reported_mappings, variant_details)


def snv_states(read, sites, minimum_base_quality):
    """Read-orientation bases from actual aligned positions; missing is never REF."""
    mapping = {r: q for q, r in read.get_aligned_pairs() if r is not None}
    states = []
    for position, reference, alternate in sites:
        q = mapping.get(position - 1)
        if q is None or read.query_qualities is None or read.query_qualities[q] < minimum_base_quality:
            states.append(None)
        else:
            base = read.query_sequence[q].upper()
            states.append('0' if base == reference else '1' if base == alternate else '?')
    return states


def fragment_configuration(rows, indices, minimum_mapping_quality):
    """One read must cover all selected sites; collapse by RG/template, retain conflict."""
    grouped = defaultdict(list)
    excluded = Counter()
    for row in rows:
        if row['flag'] & 0xF04:
            excluded['flag'] += 1
        elif row['mapping_quality'] == 255 or row['mapping_quality'] < minimum_mapping_quality:
            excluded['mapping_quality'] += 1
        elif any(row['states'][i] is None for i in indices):
            excluded['missing_or_low_quality_target'] += 1
        else:
            states = [row['states'][i] for i in indices]
            label = 'other' if '?' in states else '|'.join(states)
            grouped[(row['RG'], row['query_name'])].append((label, row['id']))
    fragments = []
    for (rg, name), observations in sorted(grouped.items()):
        labels = {label for label, identifier in observations}
        fragments.append({'RG': rg, 'query_name': name,
                          'label': next(iter(labels)) if len(labels) == 1 else 'mate_conflict',
                          'observation_ids': [identifier for label, identifier in observations]})
    counts = Counter(f['label'] for f in fragments)
    labels = ['|'.join(bits) for bits in itertools.product('01', repeat=len(indices))]
    return {'site_indices': indices, 'minimum_mapping_quality': minimum_mapping_quality,
            'counts': {label: counts[label] for label in labels + ['other', 'mate_conflict']},
            'excluded_read_counts': dict(excluded), 'fragments': fragments}


def location(row, loci):
    if row['flag'] & 4:
        return 'unmapped'
    names = [name for name, locus in loci.items() if row['chrom'] == locus['chrom']
             and row['start_0based'] < locus['end_exclusive']
             and row['end_exclusive'] is not None and row['end_exclusive'] > locus['start_0based']]
    return '+'.join(sorted(names)) if names else 'outside_annotated_loci'


def allele_context(configuration, observations, mates, mappings, loci):
    result = {}
    for label in configuration['counts']:
        fragments = [f for f in configuration['fragments'] if f['label'] == label]
        ids = {i for f in fragments for i in f['observation_ids']}
        rows = [observations[i] for i in ids]
        primary = [mappings[i]['primary'] for i in ids if i in mappings]
        mate_groups = []
        for fragment in fragments:
            selected = {}
            for identifier in fragment['observation_ids']:
                row = observations[identifier]
                if row.get('mate_id'):
                    mate = mates[row['mate_id']]
                    selected[(mate['RG'], mate['query_name'], mate['read_number'])] = mate
            mate_groups.append(list(selected.values()))
        result[label] = {
            'fragments': len(fragments), 'read_observations': len(rows),
            'by_read_group': dict(Counter(f['RG'] for f in fragments)),
            'original_read_strand': dict(Counter('reverse' if r['reverse'] else 'forward' for r in rows)),
            'original_read_MQ': dict(Counter(r['mapping_quality'] for r in rows)),
            'original_read_NM': dict(Counter(r['tags'].get('NM') for r in rows)),
            'original_read_CIGAR': dict(Counter(r['cigar'] for r in rows)),
            'single_end_primary_MQ': dict(Counter(r['mapping_quality'] for r in primary)),
            'single_end_primary_location': dict(Counter(location(r, loci) for r in primary)),
            'single_end_primary_NM': dict(Counter(r['tags'].get('NM') for r in primary)),
            'single_end_gene_minus_other_AS': dict(Counter(str(mappings[i]['comparison']['gene_minus_other_AS']) for i in ids if i in mappings)),
            'fragments_with_retrieved_mate': sum(bool(ms) for ms in mate_groups),
            'mate_primary_locations_nonexclusive_fragment_counts': {
                name: sum(any(m['mapping_quality'] != 255 and m['mapping_quality'] >= 20
                              and not m['flag'] & 0xF04 and location(m, loci) == name for m in ms)
                          for ms in mate_groups)
                for name in list(loci) + ['outside_annotated_loci']},
            'single_end_mate_locations_nonexclusive_fragment_counts': {
                name: sum(any(mappings[m['id']]['primary']['mapping_quality'] != 255
                              and mappings[m['id']]['primary']['mapping_quality'] >= 20
                              and not m['flag'] & 0xF04
                              and location(mappings[m['id']]['primary'], loci) == name for m in ms)
                          for ms in mate_groups)
                for name in list(loci) + ['outside_annotated_loci']},
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'config', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    root, output = args.root.resolve(), args.output.resolve()
    require(not output.exists() and not output.is_relative_to(Path(__file__).resolve().parents[2]),
            'Use a new output directory outside the source tree')
    config = json.loads(args.config.read_text())
    output.mkdir(parents=True)
    os.nice(10)
    tick = time.monotonic()
    state = {'status': 'running', 'pid': os.getpid(), 'started_at': datetime.datetime.now(datetime.timezone.utc).isoformat()}
    save(output/'status.json', state)
    try:
        target = root/config['target_directory']
        require(sha(target/'status.json') == config['target_status_sha256'], 'Target receipt differs')
        target_state = json.loads((target/'status.json').read_text())
        require(target_state['status'] == 'complete', 'Target analysis incomplete')
        for name, expected in target_state['output_sha256'].items():
            require(sha(target/name) == expected, 'Target output differs: ' + name)
        alignment_state = root/'state/raw-reads.json'
        require(sha(alignment_state) == target_state['source_state_sha256'], 'Alignment receipt differs')
        original = json.loads(alignment_state.read_text())
        source = Path(original['cram_path'])
        reference = root/config['bwa_reference']
        require(sha(source) == target_state['source_cram_sha256'], 'Source CRAM differs')
        require(sha(reference) == target_state['reference_sha256'], 'Reference differs')
        index_record = root/config['bwa_index_record']
        index = json.loads(index_record.read_text())
        require(index['reference_sha256'] == target_state['reference_sha256'], 'BWA reference differs')
        annotation = root/config['annotation']['path']
        require(sha(annotation) == config['annotation']['sha256'], 'Locus annotation differs')
        inputs = {name: fingerprint(path) for name, path in {
            'target_status': target/'status.json', 'alignment_state': alignment_state,
            'source': source, 'source_index': Path(str(source)+'.crai'), 'reference': reference,
            'reference_index': Path(str(reference)+'.fai'), 'bwa_index_record': index_record,
            'annotation': annotation,
            'bwa': root/config['bwa_binary'], 'configuration': args.config, 'script': Path(__file__),
            'mapping_helper': Path(__file__).with_name('review_paralog_candidate_pair.py'),
            'io_helper': Path(__file__).with_name('review_cnv_context.py')}.items()}
        for suffix, expected in index['index_sha256'].items():
            p = Path(str(reference)+suffix)
            require(sha(p) == expected, 'BWA index differs')
            inputs['bwa_index'+suffix] = fingerprint(p)
        monitored = {v['path']:(Path(v['path']).stat().st_size,Path(v['path']).stat().st_mtime_ns) for v in inputs.values()}
        chrom, (start, end) = config['chrom'], config['read_window_0based']
        sites = config['sites']
        with pysam.FastaFile(str(reference)) as fa:
            for pos, ref, alt in sites:
                require(len(ref) == len(alt) == 1 and ref != alt and fa.fetch(chrom,pos-1,pos).upper() == ref,
                        'Invalid target/reference allele')
        observations, mates, queries = {}, {}, {}
        with pysam.AlignmentFile(str(target/'targets.bam'),'rb') as bam, pysam.AlignmentFile(str(source),'rc',reference_filename=str(reference)) as original_bam:
            for number, read in enumerate(bam.fetch(chrom,start,end)):
                identifier = 'read_'+str(number)
                row = record_details(read,identifier)
                row['states'] = snv_states(read,sites,config['minimum_base_quality'])
                if not read.is_secondary and not read.is_supplementary:
                    require(read.query_sequence is not None, 'Primary read lacks sequence')
                    queries[identifier] = read.get_forward_sequence()
                    if not read.is_paired or read.mate_is_unmapped:
                        row['mate_status'] = 'unpaired_or_mate_unmapped'
                    else:
                        try:
                            mate = original_bam.mate(read)
                        except ValueError as error:
                            row['mate_status'], row['mate_error'] = 'lookup_failed', str(error)
                        else:
                            require(mate.has_tag('RG') and mate.get_tag('RG') == row['RG'] and mate.query_name == read.query_name
                                    and mate.is_read1 != read.is_read1 and mate.is_read2 != read.is_read2
                                    and not mate.is_secondary and not mate.is_supplementary, 'Mate identity differs')
                            require(mate.reference_id == read.next_reference_id and mate.reference_start == read.next_reference_start
                                    and mate.next_reference_id == read.reference_id and mate.next_reference_start == read.reference_start,
                                    'Mate coordinates are not reciprocal')
                            mid = identifier+'_mate'
                            row['mate_status'], row['mate_id'] = 'retrieved_reciprocal_identity_checked', mid
                            mates[mid] = record_details(mate,mid)
                            queries[mid] = mate.get_forward_sequence()
                observations[identifier] = row
        require(all(s and set(s.upper()) <= set('ACGTN') for s in queries.values()), 'Invalid query sequence')
        query_path = output/'queries.fa'
        query_path.write_text(''.join('>'+name+'\n'+seq+'\n' for name,seq in queries.items()))
        command = [str(root/config['bwa_binary']),'mem','-a','-Y','-t','2',str(reference),str(query_path)]
        save(output/'commands.json',[command])
        with (output/'remapped.sam').open('w') as sam, (output/'bwa.log').open('w') as log:
            subprocess.run(command,stdout=sam,stderr=log,check=True)
        remaps = defaultdict(list)
        with pysam.AlignmentFile(str(output/'remapped.sam'),'r') as sam:
            for read in sam:
                require(read.query_name in queries,'Unexpected remap query')
                row = remap_details(read)
                row['query_length'] = len(queries[read.query_name])
                remaps[read.query_name].append(row)
        require(set(remaps) == set(queries),'Missing remap query')
        mappings = {}
        for name, rows in remaps.items():
            comparison = compare_reported_mappings(rows,config['annotated_loci'][config['target_gene']])
            primary = next(r for r in rows if not r['flag'] & (256|2048))
            mappings[name] = {'comparison': comparison, 'primary': primary, 'reported_alignments': rows}
        combinations = [[i] for i in range(len(sites))] + [list(p) for p in itertools.combinations(range(len(sites)),2)] + [list(range(len(sites)))]
        phases, contexts = {}, {}
        for mq in config['mapping_quality_thresholds']:
            phases[str(mq)] = [fragment_configuration(list(observations.values()),indices,mq) for indices in combinations]
            contexts[str(mq)] = [allele_context(c,observations,mates,mappings,config['annotated_loci']) for c in phases[str(mq)][:len(sites)]]
        variants = {}
        for name in ['gatk.allsites.vcf.gz','bcftools.allsites.vcf.gz']:
            with pysam.VariantFile(str(target/name)) as vcf:
                variants[name] = [variant_details(r) for r in vcf.fetch(chrom,start,end)]
        save(output/'read-review.json',{'observations':observations,'mates':mates,'configurations':phases})
        save(output/'remapping.json',mappings)
        save(output/'regional-calls.json',variants)
        require(monitored == {p:(Path(p).stat().st_size,Path(p).stat().st_mtime_ns) for p in monitored},'Input changed')
        summary = {'completed_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'elapsed_seconds':time.monotonic()-tick,
                   'configuration':config,'inputs':inputs,'pysam':pysam.__version__,
                   'observations':len(observations),'mate_statuses':dict(Counter(r.get('mate_status','secondary_or_supplementary') for r in observations.values())),
                   'unique_original_fragments':len({(r['RG'],r['query_name']) for r in observations.values()}),
                   'configurations':{mq:[{k:v for k,v in c.items() if k!='fragments'} for c in cs] for mq,cs in phases.items()},
                   'per_site_mapping_context':contexts,
                   'outputs':{p.name:fingerprint(p) for p in output.iterdir() if p.name not in ['summary.json','status.json']},
                   'limits':['Same reads and reference; remapping is not independent locus validation. BWA-reported alternatives are not exhaustive.',
                             'Direct aligned-base states use raw base quality, no BAQ, exclude flags 0xF04 and MQ255; one read must span all sites in each configuration.',
                             'Fragments are RG/QNAME units, not independent biological samples. Per-site and mate categories overlap and must not be summed.',
                             'Missing bases are never reference. No complete assembled haplotype, unique locus origin, trans pair, clinical diagnosis or gene exclusion is established.',
                             'Annotated locus boundaries describe the pinned reference and cannot assign the origin of a sequence absent from that reference.']}
        save(output/'summary.json',summary)
        state.update(status='complete',completed_at=summary['completed_at'],summary_sha256=sha(output/'summary.json'))
        save(output/'status.json',state)
        print(json.dumps({'status':'complete','observations':len(observations),'elapsed_seconds':summary['elapsed_seconds']}))
    except BaseException as error:
        state.update(status='failed',error_type=type(error).__name__,error=str(error))
        save(output/'status.json',state)
        raise


if __name__ == '__main__':
    main()
