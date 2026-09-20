#!/usr/bin/env python3
"""Run a traceable Exomiser analysis without changing the installed application."""
import argparse
import csv
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

import yaml


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def save_json(path, record):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(record, indent=2) + '\n')
    os.replace(temporary, path)


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def samples(vcf):
    opener = gzip.open if vcf.name.endswith('.gz') else open
    with opener(vcf, 'rt') as source:
        for line in source:
            if line.startswith('#CHROM\t'):
                return line.rstrip().split('\t')[9:]
            if not line.startswith('#'):
                break
    raise ValueError('VCF sample header is missing')


def run(args):
    vcf, output, exomiser, data = (Path(p).resolve() for p in (args.vcf, args.output, args.exomiser, args.data))
    config = yaml.safe_load(Path(args.analysis).read_text())
    analysis = config['analysis']
    if analysis.get('genomeAssembly') != 'hg38':
        raise ValueError('this runner requires the pinned hg38 analysis')
    if analysis['proband'] not in samples(vcf):
        raise ValueError('configured proband is absent from the supplied VCF header')
    for name in ('%s_hg38' % args.data_version, '%s_phenotype' % args.data_version):
        if not (data / name).is_dir():
            raise ValueError('required extracted database is missing: ' + name)
    jars = list(exomiser.glob('exomiser-cli-*.jar'))
    if len(jars) != 1:
        raise ValueError('expected exactly one Exomiser CLI jar')
    if args.threads < 1:
        raise ValueError('threads must be positive')
    analysis['vcf'] = str(vcf)
    for term in args.add_hpo:
        if term not in analysis['hpoIds']:
            analysis['hpoIds'].append(term)
    config['outputOptions']['outputDirectory'] = str(output)
    output.mkdir(parents=True, exist_ok=True)
    with (output / '.analysis.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if any(path.name != '.analysis.lock' for path in output.iterdir()):
            raise ValueError('existing run preserved; choose a new output directory')
        generated = output / 'analysis.yml'
        generated.write_text(yaml.safe_dump(config, sort_keys=False))
        # A blank hg19 data-version still satisfies Spring's ConditionalOnProperty.
        # Remove the active setting from a per-run copy instead of passing an empty value.
        properties = output / 'application.properties'
        properties.write_text('\n'.join(line for line in (exomiser / 'application.properties').read_text().splitlines()
                                       if not re.match(r'\s*exomiser\.hg19\.data-version\s*[:=]', line)) + '\n')
        command = [args.java, '-XX:ActiveProcessorCount=%d' % args.threads,
                   '-Xmx' + args.heap,
                   '-Dexomiser.data-directory=' + str(data),
                   '-Dexomiser.hg38.data-version=' + args.data_version,
                   '-Dexomiser.phenotype.data-version=' + args.data_version,
                   '-Dspring.config.location=' + properties.as_uri(),
                   '-jar', str(jars[0]), 'analyse', '--analysis', str(generated)]
        status = {'status': 'preparing', 'started_at': stamp(), 'pid': os.getpid(),
                  'input_vcf': str(vcf), 'analysis_sha256': sha256(generated),
                  'application_properties_sha256': sha256(properties),
                  'runner_sha256': sha256(Path(__file__)),
                  'hpo_ids': analysis['hpoIds'], 'threads': args.threads, 'heap': args.heap,
                  'data_version': args.data_version, 'command': command,
                  'interpretation': 'Genome-wide ranking replication; prior candidate was already known.'}
        status_path = output / 'analysis-status.json'
        save_json(status_path, status)
        try:
            status['input_vcf_sha256'] = sha256(vcf)
            status['jar_sha256'] = sha256(jars[0])
            status['status'] = 'running'
            save_json(status_path, status)
            print('[exomiser] running with %d CPUs and %s heap' % (args.threads, args.heap), flush=True)
            started = time.monotonic()
            with (output / 'exomiser.log').open('w') as log:
                subprocess.run(command, cwd=exomiser, stdout=log, stderr=subprocess.STDOUT, check=True)
            tables = list(output.glob('*.genes.tsv'))
            if len(tables) != 1:
                raise ValueError('Exomiser exited without the expected gene ranking table')
            with tables[0].open() as source:
                reader = csv.DictReader(source, delimiter='\t')
                fields = reader.fieldnames or []
                rows = list(reader)
            if 'GENE_SYMBOL' not in fields:
                raise ValueError('gene ranking table has an unexpected schema')
            rank_key = '#RANK' if '#RANK' in fields else 'RANK'
            if rank_key not in fields or any(not row[rank_key].isdigit() or int(row[rank_key]) < 1 for row in rows):
                raise ValueError('gene ranking table has missing or invalid ranks')
            summary_keys = (rank_key, 'GENE_SYMBOL', 'MOI', 'EXOMISER_GENE_COMBINED_SCORE')
            gene_summary = [{key: row[key] for key in summary_keys if key in row}
                            for row in rows[:10]]
            summary = {'analysis_finished_at': stamp(), 'hpo_ids': analysis['hpoIds'],
                       'reported_gene_rows': len(rows), 'top_genes': gene_summary,
                       'distinct_reported_genes': len({row['GENE_SYMBOL'] for row in rows}),
                       'rank_unit': 'Gene and inheritance-model combination; one gene may occupy several ranks.',
                       'prior_candidate_ranks': [row.get(rank_key) for row in rows if row['GENE_SYMBOL'] == 'BUB1B'],
                       'prior_candidate_entries': [{key: row[key] for key in summary_keys if key in row}
                                                   for row in rows if row['GENE_SYMBOL'] == 'BUB1B'],
                       'scope': 'Supplied VCF ranking; independent raw-read calling and phase validation remain separate.'}
            save_json(output / 'gene-summary.json', summary)
            status.update(status='complete', completed_at=stamp(), elapsed_seconds=round(time.monotonic() - started, 2),
                          gene_table_sha256=sha256(tables[0]), gene_summary=summary)
            save_json(status_path, status)
            print('[done] verified gene ranking table and saved summary', flush=True)
            return summary
        except BaseException as error:
            status.update(status='failed', failed_at=stamp(), error_type=type(error).__name__)
            save_json(status_path, status)
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--vcf', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-e', '--exomiser', required=True)
    parser.add_argument('-d', '--data', required=True)
    parser.add_argument('-y', '--analysis', required=True)
    parser.add_argument('-t', '--threads', type=int, default=8)
    parser.add_argument('-m', '--heap', default='16g')
    parser.add_argument('--data-version', default='2512')
    parser.add_argument('--java', default=os.environ.get('JAVA', 'java'))
    parser.add_argument('--add-hpo', action='append', default=[])
    run(parser.parse_args())


if __name__ == '__main__':
    main()
