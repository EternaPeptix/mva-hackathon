#!/usr/bin/env python3
"""Reconstructed seven-example replay of need()/merged(); not a full audit test."""
import argparse
import ast
import datetime
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw = args.source.read_bytes()
    tree = ast.parse(raw)
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in {'need', 'merged'}]
    if [n.name for n in nodes] != ['need', 'merged']:
        raise ValueError('Expected exactly need and merged in that order')
    module = ast.fix_missing_locations(ast.Module(body=nodes, type_ignores=[]))
    scope = {}
    exec(compile(module, '<isolated need/merged definitions>', 'exec'), scope)
    cases = [
        ('disjoint_list', [[1, 3], [5, 7]], [[1, 3], [5, 7]]),
        ('disjoint_tuple', [(1, 3), (5, 7)], [[1, 3], [5, 7]]),
        ('overlap_list', [[1, 4], [3, 6]], [[1, 6]]),
        ('overlap_tuple', [(1, 4), (3, 6)], [[1, 6]]),
        ('adjacent_list', [[1, 4], [4, 6]], [[1, 6]]),
        ('adjacent_tuple', [(1, 4), (4, 6)], [[1, 6]]),
        ('duplicate_mixed', [[1, 4], (1, 4)], [[1, 4]]),
    ]
    results = []
    for name, values, expected in cases:
        actual = scope['merged'](values)
        if actual != expected:
            raise ValueError(name)
        results.append(dict(name=name, input=values,
                            input_item_types=[type(item).__name__ for item in values],
                            expected=expected, actual=actual, passed=True))
    receipt = dict(
        status='complete', at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source_sha256=hashlib.sha256(raw).hexdigest(),
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        isolated_function_ast_sha256=hashlib.sha256(ast.dump(module, include_attributes=False).encode()).hexdigest(),
        functions_tested=['need', 'merged'], examples_run=len(results), examples_passed=len(results),
        cases=results,
        provenance='Reconstructed standalone replay of retained seven interval-helper examples; historical inline runner was not retained. Tuple types restored from author confirmation.',
        scope='Only isolated need()/merged() definitions execute. No full-checker qualification, saved case data, GFF, alignment, external command, or network access.',
    )
    with args.output.open('x') as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({'status': receipt['status'], 'examples_passed': len(results)}))


if __name__ == '__main__':
    main()
