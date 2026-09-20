"""Reference-only CDS union and descriptive depth denominators; no genotyping."""
from collections import Counter, defaultdict
import gzip
from urllib.parse import unquote


def require(value, message):
    if not value:
        raise ValueError(message)


def prepare_targets(gff, names, reference_lengths):
    require(isinstance(names, list) and names and len(set(names)) == len(names),
            'explicit unique gene names required')
    require(all(isinstance(n, str) and n for n in names), 'invalid gene name')
    primary = set(map(str, range(1, 23))) | {'X', 'Y'}
    genes, transcripts, coding = {}, {}, []
    with gzip.open(gff, 'rt') as stream:
        for ordinal, line in enumerate(stream, 1):
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            require(len(f) == 9, 'annotation columns differ')
            if f[0] not in primary:
                continue
            attributes = {}
            for item in f[8].split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    require(key not in attributes, 'duplicate annotation attribute')
                    attributes[key] = unquote(value)
            ident = attributes.get('ID', '')
            if not (ident.startswith(('gene:', 'transcript:')) or f[2] == 'CDS'):
                continue
            start, end = int(f[3]) - 1, int(f[4])
            require(0 <= start < end and f[0] in reference_lengths
                    and end <= reference_lengths[f[0]], 'annotation outside reference')
            require(f[6] in ('+', '-'), 'unknown feature strand')
            row = dict(chrom=f[0], start=start, end=end, strand=f[6], ordinal=ordinal)
            if ident.startswith('gene:'):
                require(ident not in genes, 'duplicate gene identity')
                genes[ident] = dict(row, id=ident, name=attributes.get('Name'),
                                    biotype=attributes.get('biotype'))
            elif ident.startswith('transcript:'):
                parent = attributes.get('Parent', '').split(',')
                require(len(parent) == 1 and parent[0].startswith('gene:'), 'ambiguous transcript parent')
                require(ident not in transcripts, 'duplicate transcript identity')
                transcripts[ident] = dict(row, parent=parent[0])
            else:
                parents = attributes.get('Parent', '').split(',')
                require(parents and len(set(parents)) == len(parents)
                        and all(x.startswith('transcript:') for x in parents), 'invalid CDS parent')
                require(f[7] in ('0', '1', '2'), 'missing CDS phase')
                coding.append(dict(row, parents=parents, phase=int(f[7])))
    selected = {}
    for name in names:
        found = [g for g in genes.values() if g['name'] == name]
        require(len(found) == 1, 'gene name missing or ambiguous: ' + name)
        g = found[0]
        require(g['biotype'] == 'protein_coding', 'selected gene is not protein coding')
        selected[g['id']] = dict(g, transcript_ids=set(), cds_rows=[])
    for t, row in transcripts.items():
        require(row['parent'] in genes, 'unresolved transcript gene')
        gene = genes[row['parent']]
        require((row['chrom'], row['strand']) == (gene['chrom'], gene['strand'])
                and gene['start'] <= row['start'] < row['end'] <= gene['end'],
                'transcript geometry disagrees with gene')
    for row in coding:
        for parent in row['parents']:
            require(parent in transcripts, 'unresolved CDS transcript')
            transcript = transcripts[parent]
            require((row['chrom'], row['strand']) == (transcript['chrom'], transcript['strand'])
                    and transcript['start'] <= row['start'] < row['end'] <= transcript['end'],
                    'CDS geometry disagrees with transcript')
            gene = selected.get(transcript['parent'])
            if gene is not None:
                gene['transcript_ids'].add(parent)
                gene['cds_rows'].append(dict(row, transcript_id=parent))
    events = defaultdict(lambda: defaultdict(Counter))
    for ident, gene in selected.items():
        require(gene['cds_rows'], 'selected gene has no annotated CDS: ' + gene['name'])
        intervals = sorted({(row['start'], row['end']) for row in gene['cds_rows']})
        merged = []
        for start, end in intervals:
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        gene['cds_union'] = merged
        gene['CDS_union_bases'] = sum(b - a for a, b in merged)
        gene['transcript_ids'] = sorted(gene['transcript_ids'])
        for start, end in merged:
            events[gene['chrom']][start][ident] += 1
            events[gene['chrom']][end][ident] -= 1
    windows = []
    for chrom in reference_lengths:
        active = Counter()
        previous = None
        for position, changes in sorted(events[chrom].items()):
            members = sorted(k for k, count in active.items() if count)
            if previous is not None and previous < position and members:
                windows.append(dict(chrom=chrom, start=previous, end=position, gene_ids=members))
            active.update(changes)
            require(all(x in (0, 1) for x in active.values()), 'invalid interval membership')
            previous = position
        require(not any(active.values()), 'unclosed target interval')
    require(windows, 'no target windows')
    for ident, gene in selected.items():
        require(sum(w['end'] - w['start'] for w in windows if ident in w['gene_ids'])
                == gene['CDS_union_bases'], 'gene denominator lost during segmentation')
    return dict(genes=sorted(selected.values(), key=lambda g: names.index(g['name'])), windows=windows,
                panel_union_bases=sum(w['end'] - w['start'] for w in windows),
                gene_union_bases_sum=sum(g['CDS_union_bases'] for g in selected.values()),
                scope='Union of all annotated CDS per selected gene, across transcripts; primary nuclear reference only. No splice flank, promoter, UTR or other noncoding target.')


def summarize_depth(rows, design):
    require(len(rows) == len(design['windows']), 'target/depth row count differs')
    groups = {'panel_union': []}
    groups.update({g['id']: [] for g in design['genes']})
    for row, expected in zip(rows, design['windows']):
        require(all(row[k] == expected[k] for k in ('chrom', 'start', 'end', 'gene_ids')),
                'depth target identity differs')
        groups['panel_union'].append(row)
        for ident in expected['gene_ids']:
            groups[ident].append(row)
    summaries = {}
    for name, selected in groups.items():
        result = {}
        for label in ('span', 'ACGT'):
            hist = Counter()
            for row in selected:
                h = row[label + '_depth']['histogram']
                require(all(type(k) is int and k >= 0 and type(v) is int and v >= 0
                            for k, v in h.items()), 'invalid depth histogram')
                hist.update(h)
            total = sum(hist.values())
            result[label] = dict(bases=total, zero_depth_bases=hist[0],
                sum_depth=sum(d * n for d, n in hist.items()),
                mean=sum(d * n for d, n in hist.items()) / total if total else None,
                bases_at_least={str(t): sum(n for d, n in hist.items() if d >= t) for t in (1, 10, 20, 30)},
                histogram=dict(sorted(hist.items())))
        result['non_ACGT_reference_bases'] = result['span']['bases'] - result['ACGT']['bases']
        expected_bases = design['panel_union_bases'] if name == 'panel_union' else next(
            g['CDS_union_bases'] for g in design['genes'] if g['id'] == name)
        require(result['span']['bases'] == expected_bases and result['non_ACGT_reference_bases'] >= 0,
                'base denominator differs')
        summaries[name] = result
    return summaries
