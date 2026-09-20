#!/usr/bin/env python3
"""Build an explicit source/public-control snapshot without internal Git history.

This is a release preparation helper, not a privacy certification or a publisher.
It retains all workflow files selected by the source manifest and a bounded list
of public/synthetic/runtime receipts. Patient-result files are not selected.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import re


RECEIPTS = [
    'germline-public-inputs-20260919.json',
    'germline-public-reference-20260919-v2.json',
    'germline-public-pilot-complete-20260919.json',
    'germline-public-comparator-controls-20260919.json',
    'germline-public-comparison-20260919.json',
    'germline-public-comparison-audit-20260919.json',
    'germline-public-audit-controls-20260919.json',
    'germline-multiallelic-controls-20260919.json',
    'germline-retained-multiallelic-failed-20260919.json',
    'germline-retained-multiallelic-20260919-v2.json',
    'germline-retained-multiallelic-comparison-20260919.json',
    'germline-retained-multiallelic-comparison-audit-20260919.json',
    'germline-representation-reconciliation-20260919.json',
    'germline-representation-loss-review-20260919.json',
    'germline-representation-validation-20260919.json',
    'rtg-install-20260919.json',
    'rtg-runtime-20260919.json',
    'runtime-conda-rebuild-20260919.json',
    'runtime-python-rebuild-20260919.json',
    'runtime-structural-comparator-20260919.json',
    'runtime-wheel-preparation-20260919.json',
    'runtime-wheel-sources-20260919.json',
    'runtime-published-source-build-20260920.json',
    'runtime-published-source-validation-20260920.json',
    'runtime-published-source-public-sv-20260920.json',
    'gatk-cram-compatibility-native-control-20260920.json',
    'structural-cram30-control-20260920.json',
    'cnv-context-public-sources-20260920.json',
    'cnv-context-method-controls-20260920.json',
    'structural-cnv-context-controls-20260920.json',
]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def release_alignment_comments(data):
    text = data.decode()
    left_marker = '  # Illumina header is authoritative:\n'
    right_marker = '  # Read the first line WITHOUT'
    if text.count(left_marker) != 1 or text.count(right_marker) != 1:
        raise ValueError('Expected alignment comment block not found exactly once')
    start, stop = text.index(left_marker), text.index(right_marker)
    if stop <= start or any(not line.lstrip().startswith('#') for line in text[start:stop].splitlines()):
        raise ValueError('Release comment range includes executable text')
    replacement = (
        '  # Illumina header fields are instrument:run:flowcell:lane:tile:x:y.\n'
        '  # The acquisition-specific example is omitted from this release copy.\n'
        '  # Use flowcell.lane as the platform unit. A sample barcode can carry\n'
        '  # sequencing error and is not used to define the platform unit.\n'
    )
    result = text[:start] + replacement + text[stop:]
    def executable_lines(value):
        return [line for line in value.splitlines() if not line.lstrip().startswith('#')]
    if executable_lines(text) != executable_lines(result):
        raise ValueError('Non-comment lines changed')
    return result.encode(), {'kind': 'comment_only_acquisition_metadata_replacement',
                             'non_comment_lines_identical': True,
                             'source_comment_block_sha256': digest(text[start:stop].encode())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'manifest', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    if output.exists() or output.is_relative_to(source):
        parser.error('Use a new external snapshot directory')
    manifest_bytes = args.manifest.read_bytes()
    manifest = json.loads(manifest_bytes)
    pins = {entry['path']: entry for entry in manifest['files']}
    selected = {name: name for name in pins if name.startswith('workflow/')}
    selected['workflow/scripts/export_repository_snapshot.py'] = 'workflow/scripts/export_repository_snapshot.py'
    selected.update({'.gitignore': '.gitignore', 'LICENSE': 'LICENSE',
                     'README.md': 'publication/README.md',
                     'REPRODUCIBILITY.md': 'publication/REPRODUCIBILITY.md'})
    for name in RECEIPTS:
        selected['results/' + name] = 'results/' + name
    payloads, rows, transformations = {}, [], []
    for target, origin in sorted(selected.items()):
        if Path(target).is_absolute() or '..' in Path(target).parts or '.git' in Path(target).parts:
            raise ValueError('Unsafe release path')
        original = source / origin
        if original.is_symlink() or not original.is_file():
            raise ValueError('Release sources must be regular files')
        data = original.read_bytes()
        if origin in pins and origin != '.gitignore' and digest(data) != pins[origin]['sha256']:
            raise ValueError('Selected source differs from its sealed manifest: ' + origin)
        transformed = data
        if target == 'workflow/scripts/01_align.sh':
            transformed, note = release_alignment_comments(data)
            transformations.append(dict(note, path=target))
        payloads[target] = transformed
        rows.append({'path': target, 'source_path': origin,
                     'source_mode': original.stat().st_mode & 0o777,
                     'source_bytes': len(data),
                     'source_sha256': digest(data), 'release_bytes': len(transformed),
                     'release_sha256': digest(transformed)})
    findings = []
    patterns = {
        'credential_like_token': r'(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}|sk-[A-Za-z0-9_-]{24,}|hf_[A-Za-z0-9]{25,}|xox[baprs]-[A-Za-z0-9-]{20,}|AKIA[0-9A-Z]{16})',
        'private_key': r'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----',
        'acquisition_read_identifier': r'[A-Z0-9_-]{3,}:\d+:[A-Z0-9_-]+:\d+:\d+:\d+:\d+',
    }
    forbidden = re.compile(r'\.(?:fastq|fq|bam|cram|sam|vcf|bcf|gvcf|fa|fasta|fna)(?:\.|$)|\.(?:bai|crai|tbi|csi|pytor)$', re.I)
    for target, data in payloads.items():
        if forbidden.search(target):
            findings.append({'path': target, 'kind': 'genomic_file_extension'})
        text = data.decode()
        for label, pattern in patterns.items():
            if re.search(pattern, text):
                findings.append({'path': target, 'kind': label})
        if target.startswith('results/') and ('hc_reason' in text or 'ecdsafail' in text):
            findings.append({'path': target, 'kind': 'unrelated_project_content'})
    if findings:
        raise ValueError('Release patterns require review: ' + json.dumps(findings))
    for target in ('README.md', 'REPRODUCIBILITY.md'):
        for link in re.findall(r'\]\(([^)]+)\)', payloads[target].decode()):
            if re.match(r'[a-z]+://|#|mailto:', link):
                continue
            if link.split('#')[0] not in payloads:
                raise ValueError('Release document link not selected: ' + link)
    output.mkdir(parents=True)
    modes = {row['path']: row['source_mode'] for row in rows}
    for target, data in payloads.items():
        path = output / target
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        path.chmod(modes[target])
        if digest(path.read_bytes()) != digest(data):
            raise ValueError('Release copy differs')
        if (path.stat().st_mode & 0o777) != modes[target]:
            raise ValueError('Release permissions differ')
    summary = {'created_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'scope': 'Complete selected workflow source/configuration inventory plus public/synthetic/runtime evidence; no patient-result files or internal Git history.',
               'source_manifest_sha256': digest(manifest_bytes),
               'exporter_sha256': digest(Path(__file__).read_bytes()),
               'files': rows, 'file_count': len(rows),
               'total_release_bytes': sum(row['release_bytes'] for row in rows),
               'transformations': transformations, 'release_pattern_findings': findings,
               'limits': ['Pattern checks are not a privacy certification; this explicit scope was selected after source/history review.',
                          'Named target/HPO configurations are included; patient-derived output tables, internal history and source data are absent.',
                          'Python wheels can be rebuilt from pinned publisher inputs with the tested source recipe; standalone binaries and reference/genomic inputs are acquired separately.',
                          'A source snapshot is not completion of the patient analysis, final submissions or public release.']}
    (output / 'SOURCE-SNAPSHOT.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({'output': str(output), 'files': len(rows),
                      'bytes': summary['total_release_bytes'], 'transformations': transformations}))


if __name__ == '__main__':
    main()
