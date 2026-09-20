#!/usr/bin/env python3
"""Run pinned LIRICAL with traceable inputs, strict HPO checks and preserved outputs."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

import pysam

DATA_FILES = ['hp.json', 'phenotype.hpoa', 'hgnc_complete_set.txt', 'mim2gene_medgen',
              'en_product6.xml'] + [assembly + '_' + source + '.ser'
                                    for assembly in ('hg19', 'hg38')
                                    for source in ('ensembl', 'refseq', 'refseq_curated')]


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def stamp():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def save(path, value):
    temporary = Path(str(path) + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    os.replace(temporary, path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--data', required=True, type=Path)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--add-hpo', action='append', default=[])
    parser.add_argument('--omit-hpo', action='append', default=[])
    parser.add_argument('--vcf', type=Path)
    parser.add_argument('--exomiser-data', type=Path)
    args = parser.parse_args()
    root, data, output = (p.resolve() for p in (args.root, args.data, args.output))
    config = json.loads(args.config.read_text())
    observed = [x for x in config['observed'] if x not in args.omit_hpo] + args.add_hpo
    negated = config['negated']
    if not observed or len(set(observed + negated)) != len(observed + negated):
        raise ValueError('empty, duplicate or contradictory HPO input')
    if any(not re.fullmatch(r'HP:\d{7}', term) for term in observed + negated):
        raise ValueError('invalid HPO identifier')
    if set(args.omit_hpo) - set(config['observed']):
        raise ValueError('cannot omit an HPO term absent from the input')
    if bool(args.vcf) != bool(args.exomiser_data):
        raise ValueError('VCF and Exomiser data directory must be supplied together')
    java = root / 'tools/bio/lib/jvm/bin/java'
    jar = root / 'lirical/lirical-cli-2.4.1/lirical-cli-2.4.1.jar'
    resources = {name: digest(data / name) for name in DATA_FILES}
    command = [str(java), '-XX:ActiveProcessorCount=2', '-Xmx8g', '-jar', str(jar),
               'prioritize', '-d', str(data), '-p', ','.join(observed),
               '--sample-id', config['sample_id'], '--sex', config['sex'],
               '--validation-policy', 'STRICT', '--parallelism', '2',
               '--transcript-db', 'ENSEMBL', '-f', 'tsv', '-f', 'json', '-f', 'html',
               '-o', str(output), '-x', 'lirical']
    if negated:
        command += ['-n', ','.join(negated)]
    extra = {}
    if args.vcf:
        with pysam.VariantFile(str(args.vcf)) as vcf:
            if list(vcf.header.samples) != [config['sample_id']]:
                raise ValueError('unexpected VCF sample header')
        extra = {'input_vcf_sha256': digest(args.vcf), 'exomiser_files': {
            p.name: digest(p) for p in sorted(args.exomiser_data.glob('*.mv.db'))}}
        command += ['--vcf', str(args.vcf.resolve()), '--assembly', 'hg38',
                    '-ed38', str(args.exomiser_data.resolve())]
    output.mkdir(parents=True, exist_ok=False)
    os.nice(10)
    with (output / '.lirical.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        state = {'status': 'running', 'pid': os.getpid(), 'started_at': stamp(),
                 'command': command, 'observed': observed, 'negated': negated,
                 'config_sha256': digest(args.config), 'runner_sha256': digest(Path(__file__)),
                 'jar_sha256': digest(jar), 'resource_hashes': resources,
                 'mode': 'genotype_aware' if args.vcf else 'phenotype_only',
                 'scope': 'Algorithm comparison with previously known candidates, not blind discovery. Model probabilities are not EPCR or clinical certainty.', **extra}
        save(output / 'run-status.json', state)
        started = time.monotonic()
        try:
            with (output / 'lirical.log').open('w') as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True,
                               pass_fds=(lock.fileno(),))
            expected = [output / ('lirical.' + ext) for ext in ('tsv', 'json', 'html')]
            if any(not p.is_file() or p.stat().st_size == 0 for p in expected):
                raise ValueError('missing nonempty LIRICAL outputs')
            json.loads((output / 'lirical.json').read_text())
            if not (output / 'lirical.tsv').read_text().strip():
                raise ValueError('empty ranking table')
            state.update(status='complete', completed_at=stamp(),
                         elapsed_seconds=round(time.monotonic() - started, 3),
                         output_hashes={p.name: digest(p) for p in expected})
            save(output / 'run-status.json', state)
            print('[done] ' + str(output), flush=True)
        except BaseException as error:
            state.update(status='failed', failed_at=stamp(), error_type=type(error).__name__)
            save(output / 'run-status.json', state)
            raise


if __name__ == '__main__':
    main()
