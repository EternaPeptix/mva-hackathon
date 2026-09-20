#!/usr/bin/env python3
"""Describe mapped CDS depth for a frozen gene panel, not clinical callability."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

from coding_depth_targets import prepare_targets, require, summarize_depth


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def fingerprint(path):
    path = Path(path).resolve()
    before = path.stat()
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    after = path.stat()
    require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), 'input changed while hashing')
    return dict(path=str(path), bytes=after.st_size, mtime_ns=after.st_mtime_ns, sha256=h.hexdigest())


def save(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')
    os.replace(temporary, path)


def check_quality_availability(read):
    if read.flag & 0xF04:
        return
    require(read.mapping_quality != 255, 'unknown mapping quality needs separate handling')
    qualities = read.query_qualities
    require(qualities is not None and len(qualities) == read.query_length
            and not any(value == 255 for value in qualities),
            'unknown base qualities need separate handling; cannot treat them as passing BQ')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--expected-config-sha256', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    require(not args.output.exists(), 'existing output preserved')
    config_pin = fingerprint(args.config)
    require(config_pin['sha256'] == args.expected_config_sha256, 'configuration identity differs')
    config = json.loads(args.config.read_text())
    require(config['schema_version'] == 1 and config['assembly'] == 'GRCh38', 'configuration schema differs')
    require(config['depth_modes'] == [dict(name='MQ0_BQ1', minimum_MQ=0, minimum_BQ=1),
        dict(name='MQ20_BQ20', minimum_MQ=20, minimum_BQ=20),
        dict(name='MQ60_BQ20', minimum_MQ=60, minimum_BQ=20)], 'depth modes differ')
    require(type(config['expected_read_groups']) is int and config['expected_read_groups'] > 0,
            'invalid read-group count')
    require(set(config['inputs']) == {'alignment', 'alignment_index', 'reference', 'reference_index', 'gff', 'samtools'},
            'input identities incomplete')
    expected_dependencies = {'measure_coding_depth.py', 'coding_depth_targets.py',
                             'survey_chromosome_depth.py', 'prepare_gatk_compatible_cram.py'}
    require(set(config['source_sha256']) == expected_dependencies, 'source identities incomplete')
    pins = {'config': config_pin}
    for name, expected in config['inputs'].items():
        pins[name] = fingerprint(expected['path'])
        require(pins[name]['sha256'] == expected['sha256'], 'input identity differs: ' + name)
    for name, expected in config['source_sha256'].items():
        pins['source/' + name] = fingerprint(Path(__file__).with_name(name))
        require(pins['source/' + name]['sha256'] == expected, 'source identity differs: ' + name)
    alignment = Path(pins['alignment']['path'])
    reference = Path(pins['reference']['path'])
    require(pins['reference_index']['path'] == str(reference) + '.fai', 'reference index path differs')
    require(pins['alignment_index']['path'] == str(alignment) + ('.crai' if alignment.suffix == '.cram' else '.bai'),
            'alignment index path differs')
    import pysam
    from survey_chromosome_depth import collect_depth
    os.nice(10)
    with pysam.FastaFile(str(reference)) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
        design = prepare_targets(pins['gff']['path'], config['genes'], lengths)
        sequences = [fasta.fetch(w['chrom'], w['start'], w['end']).upper() for w in design['windows']]
    samtools = Path(pins['samtools']['path'])
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    state = dict(status='running', pid=os.getpid(), started_at=stamp(), inputs=pins,
                 scope='Descriptive mapped-read depth on a fixed reference-only CDS union; no genotype or clinical negative inference.')
    save(out / 'status.json', state)
    commands = []
    started = time.monotonic()
    def run(command, label):
        record = dict(label=label, command=list(map(str, command)), started_at=stamp())
        commands.append(record)
        save(out / 'commands.json', commands)
        with (out / (label + '.log')).open('w') as log:
            result = subprocess.run(record['command'], stdout=log, stderr=subprocess.STDOUT)
        record.update(completed_at=stamp(), exit_code=result.returncode)
        save(out / 'commands.json', commands)
        require(result.returncode == 0, 'command failed: ' + label)
    try:
        design.update(status='prepared', prepared_at=stamp(), config=config,
                      coordinates='0-based half-open; all CDS phases retained, not translated')
        save(out / 'selection.json', design)
        bed = out / 'regions.bed'
        bed.write_text(''.join(f"{w['chrom']}\t{w['start']}\t{w['end']}\n" for w in design['windows']))
        bam = out / 'selected.bam'
        run([samtools, 'view', '-@', '1', '--no-PG', '-T', reference, '-b', '-M', '-L', bed,
             '--save-counts', out / 'extraction-counts.json', '-o', bam, alignment], 'extract-preserving-records')
        run([samtools, 'quickcheck', '-v', bam], 'quickcheck')
        run([samtools, 'index', bam], 'index')
        counts, quality, flags = Counter(), Counter(), Counter()
        with pysam.AlignmentFile(str(bam)) as reads:
            groups = reads.header.to_dict().get('RG', [])
            require(len(groups) == config['expected_read_groups'] and len({g['ID'] for g in groups}) == len(groups)
                    and {g.get('SM') for g in groups} == {config['sample']}, 'sample/read-group identity differs')
            require(list(zip(reads.references, reads.lengths)) == list(lengths.items()), 'alignment/reference dictionaries differ')
            allowed_groups = {g['ID'] for g in groups}
            for read in reads:
                require(read.has_tag('RG') and read.get_tag('RG') in allowed_groups, 'unregistered read group')
                counts[read.get_tag('RG')] += 1
                flags[read.flag] += 1
                check_quality_availability(read)
                if not read.flag & 0xF04:
                    quality[read.mapping_quality] += 1
        require(not quality[255], 'unknown mapping quality needs separate handling')
        results = {}
        for mode in config['depth_modes']:
            name = mode['name']
            command = [str(x) for x in [samtools, 'depth', '-aa', '-b', bed, '-q', mode['minimum_BQ'],
                '-Q', mode['minimum_MQ'], '-G', '0xF04', '-s', bam]]
            row = dict(label=name, command=command, started_at=stamp())
            commands.append(row)
            save(out / 'commands.json', commands)
            with (out / (name + '.log')).open('w') as log, gzip.open(out / (name + '.depth.tsv.gz'), 'wt') as raw:
                process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log, text=True)
                try:
                    windows, n = collect_depth(process.stdout, design['windows'], sequences, raw)
                except BaseException:
                    process.terminate()
                    process.wait()
                    raise
                finally:
                    process.stdout.close()
                code = process.wait()
            row.update(completed_at=stamp(), exit_code=code)
            save(out / 'commands.json', commands)
            require(code == 0 and n == design['panel_union_bases'], 'native depth failed or denominator differs')
            result = dict(mode=mode, position_count=n, groups=summarize_depth(windows, design), windows=windows)
            save(out / (name + '.json'), result)
            results[name] = dict(mode=mode, position_count=n, groups=result['groups'])
        for pin in pins.values():
            stat = Path(pin['path']).stat()
            require((stat.st_size, stat.st_mtime_ns) == (pin['bytes'], pin['mtime_ns']), 'input changed during analysis')
        summary = dict(status='complete', completed_at=stamp(), inputs=pins,
            genes=[{k: g[k] for k in ('id', 'name', 'chrom', 'strand', 'CDS_union_bases', 'transcript_ids')} for g in design['genes']],
            panel_union_bases=design['panel_union_bases'], gene_union_bases_sum=design['gene_union_bases_sum'],
            read_records_by_group={g: counts[g] for g in sorted(allowed_groups)}, flags=dict(flags),
            eligible_flag_mapping_qualities=dict(quality), depth_modes=results,
            limitations=['CDS-only union across all annotated transcripts of the declared genes; no splice, regulatory, repeat or genome-wide coverage claim.',
                'Explicit zero depth and non-ACGT reference bases remain separate; no uncovered base is filled with a reference genotype.',
                'Thresholds 1/10/20/30 are descriptive depth bins, not validated genotype sensitivity or clinical thresholds.',
                'Native depth excludes deletions and reference skips, applies no BAQ, excludes flags0xF04 and suppresses overlapping mates using pinned samtools rules.',
                'Depth is retained aligned-read evidence after mate suppression, not independently established molecule count, locus origin or allele balance.',
                'Reference selection is fixed before reads. No GC adjustment, external mappability mask, contamination estimate, or diploid/mosaic sensitivity is inferred.',
                'Quality modes reuse the same reads. Genes can share CDS bases; their totals must not be summed as independent panel bases.'])
        save(out / 'summary.json', summary)
        state.update(status='complete', completed_at=stamp(), elapsed_seconds=time.monotonic() - started,
            summary_sha256=fingerprint(out / 'summary.json')['sha256'],
            outputs={p.name: fingerprint(p) for p in out.iterdir() if p.is_file() and p.name != 'status.json'})
        save(out / 'status.json', state)
        print(json.dumps(dict(status='complete', genes=len(design['genes']), panel_union_bases=design['panel_union_bases'])))
    except BaseException as error:
        state.update(status='failed', failed_at=stamp(), error=repr(error))
        save(out / 'status.json', state)
        raise


if __name__ == '__main__':
    main()
