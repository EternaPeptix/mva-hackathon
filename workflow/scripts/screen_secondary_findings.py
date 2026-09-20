#!/usr/bin/env python3
"""Reproducible ClinVar/ACMG-SF exact-allele triage, not clinical reporting.

Run on a compute host. Public annotations and individual genotypes never leave
the host. Both sides are reference checked and left normalized by bcftools.
Disease-specific eligibility and inheritance require review of the matches.
"""
import argparse
from collections import Counter, defaultdict
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

import pysam

PRIMARY = tuple(str(x) for x in range(1, 23)) + ('X', 'Y')
PLP = {'Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic'}
REVIEW_STARS = {
    'practice_guideline': 4,
    'reviewed_by_expert_panel': 3,
    'criteria_provided,multiple_submitters,no_conflicts': 2,
    'criteria_provided,single_submitter': 1,
    'criteria_provided,conflicting_classifications': 1,
    'criteria_provided,conflicting_interpretations': 1,
    'no_assertion_criteria_provided': 0,
    'no_classification_provided': 0,
    'no_assertion_provided': 0,
    'no_classification_for_the_individual_variant': 0,
    # Current VCF spelling of the individual-variant/haplotype-only status.
    'no_classification_for_the_single_variant': 0,
    'no_classifications_from_unflagged_records': 0,
    'no_assertion_for_the_individual_variant': 0,
    '': 0,
}


def timestamp():
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())


def fingerprint(path):
    path = Path(path)
    before = path.stat()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError('input changed during hashing: ' + path.name)
    return {'path': str(path.resolve()), 'bytes': after.st_size,
            'mtime_ns': after.st_mtime_ns, 'sha256': digest.hexdigest()}


def values(info, field):
    if field not in info:
        return ''
    value = info.get(field)
    if value is None:
        return ''
    return ','.join(str(x) for x in value) if isinstance(value, tuple) else str(value)


def has_plp_signal(text):
    for term in re.split(r'[|,]', text):
        # VCF underscores encode spaces, including the space before the
        # submission-count suffix: Likely_pathogenic_(1). Haplotype entries
        # prefix the term with a Variation ID, e.g. 440849:Pathogenic.
        term = re.sub(r'^\d+:', '', term.strip('_ '))
        term = re.sub(r'_?\(\d+\)$', '', term).strip('_ ')
        if term in PLP:
            return True
    return False


def germline_stratum(info):
    """Only a clean aggregate P/LP classification enters the P/LP strata.

    Conflicting and haplotype-only P/LP assertions remain review records. Somatic
    oncogenicity/clinical-impact fields never count as germline pathogenicity.
    Unknown review strings fail closed, rather than silently receiving stars.
    """
    classification = values(info, 'CLNSIG')
    status = values(info, 'CLNREVSTAT').replace(',_', ',')
    parts = set(re.split(r'[|,]', classification))
    clean = bool(parts) and parts <= PLP
    conflict = values(info, 'CLNSIGCONF')
    included = values(info, 'CLNSIGINCL')
    if not (clean or has_plp_signal(classification) or has_plp_signal(conflict) or has_plp_signal(included)):
        return None
    if status not in REVIEW_STARS:
        raise ValueError('unrecognized ClinVar review status: ' + status)
    stars = REVIEW_STARS[status]
    if clean and stars >= 2:
        return 'aggregate_plp_2plus_stars', stars
    if clean and stars == 1 and 'conflict' not in status:
        return 'aggregate_plp_1star', stars
    if clean and stars == 0:
        return 'aggregate_plp_0star_review', stars
    if has_plp_signal(classification) or has_plp_signal(conflict):
        return 'conflicting_or_mixed_classification_review', stars
    return 'haplotype_or_genotype_only_review', stars


def gene_symbols(info):
    return {part.rsplit(':', 1)[0] for part in values(info, 'GENEINFO').split('|') if ':' in part}


def concrete_alleles(record):
    return bool(record.alts) and all(re.fullmatch('[ACGT]+', allele) for allele in record.alleles)


def allele_key(record):
    if len(record.alts or ()) != 1:
        raise ValueError('expected a normalized biallelic record')
    return (record.contig, record.pos, record.ref, record.alts[0])


@contextmanager
def normalized_stream(bcftools, reference, source, log_path, regions=None):
    """Bounded streaming normalization; retain the original allele identity.

    No atomization of complex alleles. Splitting a multi-ALT record does not
    establish phase; source GT/PGT/PID are recovered from the original VCF.
    """
    commands = []
    processes = []
    with open(log_path, 'wb') as log:
        try:
            command = [str(bcftools), 'norm', '--fasta-ref', str(reference),
                       '--check-ref', 'e', '--multiallelics', '-any',
                       '--old-rec-tag', 'SF_ORIG', '-Ou']
            if regions:
                view = [str(bcftools), 'view', '--regions', ','.join(regions), '-Ou', str(source)]
                commands.append(view)
                processes.append(subprocess.Popen(view, stdout=subprocess.PIPE, stderr=log))
                process = subprocess.Popen(command, stdin=processes[0].stdout,
                                           stdout=subprocess.PIPE, stderr=log)
                processes[0].stdout.close()
            else:
                command.append(str(source))
                process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log)
            commands.append(command)
            processes.append(process)
            with pysam.VariantFile(process.stdout) as reader:
                yield reader, commands
            process.stdout.close()
            codes = [p.wait() for p in processes]
            if any(codes):
                raise RuntimeError('normalization failed; see ' + str(log_path))
        finally:
            for p in processes:
                if p.poll() is None:
                    p.terminate()
                    try:
                        p.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        p.kill()
                        p.wait()


def original_call(record, source, sample_name):
    marker = values(record.info, 'SF_ORIG')
    if marker:
        fields = marker.split('|')
        # bcftools omits USED_ALT_IDX when only left alignment changed a
        # biallelic record; it includes it when splitting multiple ALTs.
        if len(fields) == 4 and ',' not in fields[3]:
            fields.append('1')
        if len(fields) != 5:
            raise ValueError('unsupported original-allele tag')
        chrom, pos, ref, alternate, index = fields
        pos, index = int(pos), int(index)
        alternatives = tuple(alternate.split(','))
    else:
        chrom, pos, ref, alternatives, index = record.contig, record.pos, record.ref, record.alts, 1
    matches = [r for r in source.fetch(chrom, pos - 1, pos)
               if r.pos == pos and r.ref == ref and r.alts == alternatives]
    if len(matches) != 1 or not 1 <= index <= len(alternatives):
        raise ValueError('original VCF allele cannot be uniquely recovered')
    original = matches[0]
    sample = original.samples[sample_name]
    gt = sample.get('GT') or ()
    if any(x is not None and not 0 <= x <= len(alternatives) for x in gt):
        raise ValueError('invalid source genotype allele index')
    ad = sample.get('AD')
    delimiter = '|' if sample.phased else '/'
    return {
        'source_allele': {'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alternatives[index - 1]},
        'source_all_alts': list(alternatives), 'source_alt_index': index,
        'GT': delimiter.join('.' if x is None else str(x) for x in gt),
        'GT_values': list(gt), 'GT_phased': sample.phased,
        'alt_copies': gt.count(index), 'called_ploidy': len(gt),
        'fully_called': bool(gt) and None not in gt,
        'DP': sample.get('DP'), 'GQ': sample.get('GQ'),
        'AD': list(ad) if ad is not None else None,
        'alt_depth': ad[index] if ad is not None and len(ad) > index else None,
        'PGT': sample.get('PGT'), 'PID': sample.get('PID'), 'PS': sample.get('PS'),
        'FILTER': list(original.filter), 'QUAL': original.qual,
    }


def quality_failures(call, min_dp=10, min_gq=20):
    reasons = []
    if call['alt_copies'] <= 0:
        reasons.append('no_called_alternate')
    if not call['fully_called']:
        reasons.append('partial_or_missing_genotype')
    if call['FILTER'] != ['PASS']:
        reasons.append('record_not_explicit_PASS')
    if call['DP'] is None or call['DP'] < min_dp:
        reasons.append('depth_missing_or_below_threshold')
    if call['GQ'] is None or call['GQ'] < min_gq:
        reasons.append('genotype_quality_missing_or_below_threshold')
    return reasons


def load_panel(path, expected_count):
    data = json.loads(Path(path).read_text())
    genes = {row['symbol']: row for row in data['rows']}
    if len(genes) != expected_count or len(data['rows']) != expected_count or data['total'] != expected_count:
        raise ValueError('unexpected panel membership/count; review the source before screening')
    return genes


def run(args):
    os.nice(10)
    started = time.monotonic()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    state = {'status': 'running', 'started_at': timestamp(), 'pid': os.getpid(),
             'scope': 'Exact normalized ClinVar allele triage in ACMG SF genes; NOT final reportability.',
             'sample': args.sample, 'thresholds': {'min_DP': args.min_dp, 'min_GQ': args.min_gq},
             'commands': []}
    status_path = output / 'run-status.json'
    status_path.write_text(json.dumps(state, indent=2) + '\n')
    try:
        genes = load_panel(args.panel, args.expected_genes)
        state['inputs'] = {name: fingerprint(path) for name, path in {
            'vcf': args.vcf, 'vcf_index': Path(str(args.vcf) + '.tbi'),
            'clinvar': args.clinvar, 'clinvar_index': Path(str(args.clinvar) + '.tbi'),
            'panel': args.panel, 'reference': args.reference,
            'reference_index': Path(str(args.reference) + '.fai'),
            'script': Path(__file__), 'bcftools': args.bcftools}.items()}
        state['tools'] = {'pysam': pysam.__version__, 'bcftools': subprocess.check_output(
            [str(args.bcftools), '--version'], text=True).splitlines()[:2]}
        state['panel_gene_count'] = len(genes)
        public_counts, sample_counts = Counter(), Counter()
        selected = output / 'public-clinvar-selected.vcf.gz'
        with pysam.VariantFile(args.clinvar) as reader:
            if '##reference=GRCh38\n' not in str(reader.header):
                raise ValueError('ClinVar input is not declared GRCh38')
            if 'SF_ORIG' in reader.header.info:
                raise ValueError('reserved normalization tag already present')
            state['clinvar_file_date'] = [line for line in str(reader.header).splitlines()
                                          if line.startswith('##fileDate=')]
            header = reader.header.copy()
            with pysam.FastaFile(args.reference) as reference:
                for chrom in PRIMARY:
                    if chrom not in header.contigs:
                        header.contigs.add(chrom, length=reference.get_reference_length(chrom))
            with pysam.VariantFile(selected, 'wz', header=header) as writer:
                for record in reader:
                    public_counts['all_records'] += 1
                    if not gene_symbols(record.info) & genes.keys():
                        continue
                    public_counts['panel_gene_records'] += 1
                    stratum = germline_stratum(record.info)
                    if stratum is None:
                        continue
                    public_counts['panel_records_with_plp_signal'] += 1
                    if record.contig not in PRIMARY or not concrete_alleles(record):
                        public_counts['non_primary_or_non_concrete_excluded'] += 1
                        continue
                    public_counts['selected_' + stratum[0]] += 1
                    copy = record.copy()
                    copy.translate(header)
                    writer.write(copy)
        lookup = defaultdict(list)
        with normalized_stream(args.bcftools, args.reference, selected,
                               output / 'clinvar-normalization.log') as (reader, commands):
            state['commands'].extend(commands)
            for record in reader:
                key = allele_key(record)
                stratum, stars = germline_stratum(record.info)
                lookup[key].append({'variation_id': record.id,
                    'genes': sorted(gene_symbols(record.info) & genes.keys()),
                    'classification': values(record.info, 'CLNSIG'),
                    'conflicting_classifications': values(record.info, 'CLNSIGCONF'),
                    'included_haplotype_classifications': values(record.info, 'CLNSIGINCL'),
                    'review_status': values(record.info, 'CLNREVSTAT'), 'stars': stars,
                    'stratum': stratum, 'conditions': values(record.info, 'CLNDN'),
                    'condition_ids': values(record.info, 'CLNDISDB'),
                    'molecular_consequences': values(record.info, 'MC'),
                    'clinvar_source_allele': values(record.info, 'SF_ORIG')})
        public_counts['normalized_unique_alleles'] = len(lookup)
        hits = []
        with pysam.VariantFile(args.vcf) as source:
            if list(source.header.samples) != [args.sample]:
                raise ValueError('expected exactly the specified sample')
            if 'SF_ORIG' in source.header.info:
                raise ValueError('reserved normalization tag already present')
            regions = [chrom for chrom in PRIMARY if chrom in source.header.contigs]
            if not set(str(i) for i in range(1, 23)) <= set(regions):
                raise ValueError('input must have matching, unprefixed nuclear contig names')
            state['contigs_in_scope'] = regions
            state['contigs_outside_scope'] = [c for c in source.header.contigs if c not in regions]
            with normalized_stream(args.bcftools, args.reference, args.vcf,
                                   output / 'sample-normalization.log', regions) as (reader, commands):
                state['commands'].extend(commands)
                for record in reader:
                    sample_counts['normalized_nuclear_allele_records'] += 1
                    if not concrete_alleles(record):
                        sample_counts['non_concrete_alleles_excluded'] += 1
                        continue
                    key = allele_key(record)
                    if key not in lookup:
                        continue
                    sample_counts['exact_matching_allele_records'] += 1
                    call = original_call(record, source, args.sample)
                    failures = quality_failures(call, args.min_dp, args.min_gq)
                    for failure in failures:
                        sample_counts['matching_records_' + failure] += 1
                    if call['alt_copies'] <= 0:
                        continue
                    sample_counts['called_alt_matching_records'] += 1
                    hits.append({'normalized_allele': dict(zip(('chrom', 'pos', 'ref', 'alt'), key)),
                                 'source_call': call, 'quality_failures': failures,
                                 'clinvar': lookup[key],
                                 'reportability': 'pending_disease_variant_and_inheritance_review'})
        # Do not count duplicated normalized records as independent variants.
        identities = [tuple(hit['normalized_allele'].values()) for hit in hits]
        if len(identities) != len(set(identities)):
            raise ValueError('duplicate normalized subject alleles require explicit reconciliation')
        pass_by_stratum = Counter()
        for hit in hits:
            if not hit['quality_failures']:
                pass_by_stratum.update({x['stratum'] for x in hit['clinvar']})
        result = {
            'status': 'complete', 'sample': args.sample, 'assembly': 'GRCh38',
            'panel_version': 'ACMG SF v3.3, live ClinGen snapshot',
            'panel_genes': sorted(genes), 'public_counts': dict(public_counts),
            'sample_counts': dict(sample_counts),
            'quality_pass_unique_alleles_by_stratum': dict(pass_by_stratum),
            'matches': hits,
            'limitations': [
                'This output is triage, not clinical classification, ACMG reportability, or absence of secondary findings.',
                'Only exact reference-normalized concrete alleles with ClinVar GENEINFO membership in the panel are screened.',
                'No de novo classification, read validation, copy-number/structural/repeat/mosaic analysis, or callability guarantee.',
                'Complex alleles are not atomized; equivalent multi-record haplotypes may be missed.',
                'No frequency filter was applied. VUS-only and unclassified variants are not classified by this screen.',
                'Aggregate ClinVar P/LP may concern an ineligible disease. All positives require condition-specific review.',
                'AR pairs, sex-related eligibility, phase, special variant restrictions and flagged evidence need review.',
                'The 2-star/1-star split is a review stratum chosen here, not an ACMG eligibility requirement.',
                'PASS, DP and GQ thresholds are analytical choices, not a clinical accuracy guarantee.',
            ]}
        (output / 'triage.json').write_text(json.dumps(result, indent=2) + '\n')
        # Preserve current public guidance only for matched genes, outside the repo.
        matched_genes = {g for h in hits for c in h['clinvar'] for g in c['genes']}
        (output / 'matched-gene-guidance.json').write_text(json.dumps(
            {g: genes[g] for g in sorted(matched_genes)}, indent=2) + '\n')
        for name, record in state['inputs'].items():
            now = Path(record['path']).stat()
            if (now.st_size, now.st_mtime_ns) != (record['bytes'], record['mtime_ns']):
                raise ValueError('input changed during screen: ' + name)
        state.update(status='complete', completed_at=timestamp(), elapsed_seconds=time.monotonic() - started,
                     output_hashes={name: fingerprint(output / name)['sha256'] for name in
                                    ('triage.json', 'matched-gene-guidance.json', 'public-clinvar-selected.vcf.gz')})
        print(json.dumps({'status': state['status'], 'public_counts': dict(public_counts),
                          'sample_counts': dict(sample_counts), 'quality_pass_by_stratum': dict(pass_by_stratum),
                          'elapsed_seconds': state['elapsed_seconds']}), flush=True)
    except BaseException as error:
        state.update(status='failed', failed_at=timestamp(), error=str(error),
                     elapsed_seconds=time.monotonic() - started)
        raise
    finally:
        status_path.write_text(json.dumps(state, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('vcf', 'clinvar', 'panel', 'reference', 'bcftools', 'output'):
        parser.add_argument('--' + name, required=True, type=Path)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--expected-genes', type=int, default=84)
    parser.add_argument('--min-dp', type=int, default=10)
    parser.add_argument('--min-gq', type=int, default=20)
    args = parser.parse_args()
    if min(args.min_dp, args.min_gq) < 0:
        parser.error('quality thresholds must be nonnegative')
    run(args)


if __name__ == '__main__':
    main()
