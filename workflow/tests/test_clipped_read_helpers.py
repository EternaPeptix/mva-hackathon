"""Fictional controls for complete-query fits, independent of native read I/O."""
import itertools
from pathlib import Path
import random
import re
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from clipped_read_helpers import fit_full_query


def global_distance(query, target):
    """Independent ordinary global distance used on every target substring."""
    row = list(range(len(target) + 1))
    for i, base in enumerate(query, 1):
        next_row = [i]
        for j, other in enumerate(target, 1):
            next_row.append(min(row[j] + 1, next_row[-1] + 1,
                                row[j - 1] + (base != other)))
        row = next_row
    return row[-1]


def substring_oracle(query, target):
    outcomes = [(global_distance(query, target[start:end]), start, end)
                for start in range(len(target) + 1)
                for end in range(start, len(target) + 1)]
    best = min(row[0] for row in outcomes)
    return best, {end for cost, _, end in outcomes if cost == best}


def exhaustive_short_oracle(query, target):
    """Enumerate edit scripts directly; no production DP/traceback recurrence."""
    rows = []
    def walk(start, i, j, operations, cost):
        if i == len(query):
            # Any trailing deletion is more costly than ending at this target
            # coordinate, so it cannot improve a query-global/target-local fit.
            rows.append((cost, start, j, operations))
            return
        walk(start, i + 1, j, operations + 'I', cost + 1)
        if j < len(target):
            operation = '=' if query[i] == target[j] else 'X'
            walk(start, i + 1, j + 1, operations + operation, cost + (operation == 'X'))
            # A leading target deletion is a free omitted target prefix instead.
            if i > 0:
                walk(start, i, j + 1, operations + 'D', cost + 1)
    for start in range(len(target) + 1):
        walk(start, 0, start, '', 0)
    best = min(row[0] for row in rows)
    return {(start, end, operations) for cost, start, end, operations in rows if cost == best}


def expand_cigar(cigar):
    parts = re.findall(r'([1-9][0-9]*)([=XID])', cigar)
    if ''.join(number + operation for number, operation in parts) != cigar:
        raise AssertionError('Malformed CIGAR: ' + cigar)
    return ''.join(operation * int(number) for number, operation in parts)


class FullQueryFitTests(unittest.TestCase):
    def assert_replays(self, query, target, result):
        identities = set()
        for alignment in result['alignments']:
            i, j, score = 0, alignment['start0'], 0
            self.assertLessEqual(j, alignment['end0'])
            self.assertGreaterEqual(j, 0)
            self.assertLessEqual(alignment['end0'], len(target))
            operations = expand_cigar(alignment['cigar'])
            identity = (j, alignment['end0'], operations)
            self.assertNotIn(identity, identities)
            identities.add(identity)
            for operation in operations:
                if operation in '=X':
                    self.assertEqual(query[i] == target[j], operation == '=')
                    i += 1
                    j += 1
                    score += operation == 'X'
                elif operation == 'I':
                    i += 1
                    score += 1
                else:
                    j += 1
                    score += 1
            self.assertEqual(i, len(query))
            self.assertEqual(j, alignment['end0'])
            self.assertEqual((alignment['query_start0'], alignment['query_end0']), (0, len(query)))
            self.assertEqual(score, result['best_cost'])
            self.assertEqual(alignment['edit_cost'], score)
            self.assertIn(j, result['optimal_end0'])
        return identities

    def test_short_strings_against_independent_substring_and_script_oracles(self):
        strings = [''.join(bases) for length in range(1, 4)
                   for bases in itertools.product('AC', repeat=length)]
        for query in strings:
            for target in strings:
                with self.subTest(query=query, target=target):
                    result = fit_full_query(query, target)
                    best, ends = substring_oracle(query, target)
                    self.assertEqual(result['best_cost'], best)
                    self.assertEqual(set(result['optimal_end0']), ends)
                    expected = exhaustive_short_oracle(query, target)
                    observed = self.assert_replays(query, target, result)
                    self.assertEqual(observed, expected)
                    self.assertEqual(result['optimal_path_count_capped'], len(expected))
                    self.assertTrue(result['traceback_complete'])

    def test_exact_fit_and_query_ends_are_not_freely_clipped(self):
        exact = fit_full_query('ACGT', 'TTACGTAA')
        self.assertEqual(exact['best_cost'], 0)
        self.assertEqual(exact['alignments'][0]['cigar'], '4=')
        self.assertEqual((exact['alignments'][0]['start0'], exact['alignments'][0]['end0']), (2, 6))
        self.assertTrue(exact['unique_optimal_alignment'])
        for query in ('TACGT', 'ACGTA', 'TACGTA'):
            result = fit_full_query(query, 'ACGT')
            self.assertEqual(result['best_cost'], len(query) - 4)
            self.assert_replays(query, 'ACGT', result)

    def test_149_base_true_two_base_insertion_and_novel_123_base_clip(self):
        randomizer = random.Random(12092026)
        reference = ''.join(randomizer.choice('ACGT') for _ in range(400))
        alternate = reference[:220] + 'GC' + reference[220:]
        query = reference[150:220] + 'GC' + reference[220:297]
        self.assertEqual(len(query), 149)
        reference_fit = fit_full_query(query, reference)
        alternate_fit = fit_full_query(query, alternate)
        self.assertEqual(reference_fit['best_cost'], 2)
        self.assertEqual(alternate_fit['best_cost'], 0)
        self.assert_replays(query, reference, reference_fit)
        self.assert_replays(query, alternate, alternate_fit)
        long_clip = 'A' * 26 + 'C' * 123
        result = fit_full_query(long_clip, 'A' * 150 + 'CC' + 'A' * 150)
        self.assertEqual(result['best_cost'], 121)
        self.assert_replays(long_clip, 'A' * 150 + 'CC' + 'A' * 150, result)
        self.assertEqual(fit_full_query(long_clip, 'A' * 300)['best_cost'], 123)

    def test_repeat_endpoints_traceback_ties_and_cap_are_explicit(self):
        repeat = fit_full_query('AA', 'AAAA', max_paths=2)
        self.assertEqual(repeat['best_cost'], 0)
        self.assertEqual(repeat['optimal_end0'], [2, 3, 4])
        self.assertEqual(repeat['optimal_path_count_capped'], 3)
        self.assertEqual(len(repeat['alignments']), 2)
        self.assertFalse(repeat['path_count_is_exact'])
        self.assertFalse(repeat['traceback_complete'])
        self.assertFalse(repeat['unique_optimal_alignment'])
        self.assertTrue(repeat['saturated'])
        self.assertEqual(repeat['status'], 'capped')
        self.assertEqual(repeat['optimal_path_count_exact'], 3)
        self.assert_replays('AA', 'AAAA', repeat)
        tied = fit_full_query('AAA', 'A', max_paths=3)
        self.assertEqual(tied['optimal_end0'], [1])
        self.assertEqual(tied['optimal_path_count_capped'], 3)
        self.assertTrue(tied['path_count_is_exact'])
        self.assertTrue(tied['traceback_complete'])
        self.assertFalse(tied['unique_optimal_alignment'])
        self.assertEqual({row['cigar'] for row in tied['alignments']}, {'2I1=', '1I1=1I', '1=2I'})
        overflow = fit_full_query('A' * 80, 'A' * 40, max_paths=2)
        self.assertTrue(overflow['saturated'])
        self.assertEqual(overflow['optimal_path_count_capped'], 3)
        self.assertEqual(len(overflow['alignments']), 2)
        self.assert_replays('A' * 80, 'A' * 40, overflow)

    def test_exact_substring_fast_path_preserves_overlaps_and_cell_contract(self):
        query, target = 'ACAC', 'ACACACAC'
        result = fit_full_query(query, target, max_paths=1)
        self.assertTrue(result['exact_substring_fast_path'])
        self.assertEqual(result['optimal_end0'], [4, 6, 8])
        self.assertEqual(result['optimal_path_count_exact'], 3)
        self.assertEqual(result['optimal_path_count_capped'], 2)
        self.assertFalse(result['path_count_is_exact'])
        self.assertFalse(result['traceback_complete'])
        self.assertEqual(result['dp_cells'], 45)
        self.assertEqual(result['dp_cells_evaluated'], 0)
        self.assertEqual(result['traceback_steps'], 0)
        self.assert_replays(query, target, result)
        with self.assertRaisesRegex(ValueError, 'cell cap'):
            fit_full_query(query, target, max_cells=44)
        nonexact = fit_full_query('AAC', 'AAT')
        self.assertFalse(nonexact['exact_substring_fast_path'])
        self.assertEqual(nonexact['dp_cells_evaluated'], nonexact['dp_cells'])

    def test_target_empty_segment_is_a_costly_placement_not_missing_query(self):
        result = fit_full_query('A', 'C')
        self.assertEqual(result['best_cost'], 1)
        self.assertEqual(result['optimal_end0'], [0, 1])
        self.assertEqual(self.assert_replays('A', 'C', result),
                         {(0, 0, 'I'), (0, 1, 'X'), (1, 1, 'I')})

    def test_clear_alphabet_bounds_and_no_recursive_traceback(self):
        for query, target in [('', 'A'), ('A', ''), ('N', 'N'), ('a', 'A'), ('A-', 'AA'), (None, 'A')]:
            with self.subTest(query=query, target=target), self.assertRaises(ValueError):
                fit_full_query(query, target)
        for kwargs in ({'max_cells': 3}, {'max_cells': 0}, {'max_cells': True},
                       {'max_cells': 10000001}, {'max_paths': 0}, {'max_paths': False},
                       {'max_paths': 4097}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                fit_full_query('A', 'A', **kwargs)
        self.assertEqual(fit_full_query('A', 'A', max_cells=4)['dp_cells'], 4)
        long = fit_full_query('A' * 1500, 'A', max_paths=1)
        self.assertTrue(long['saturated'])
        self.assertEqual(long['best_cost'], 1499)
        self.assert_replays('A' * 1500, 'A', long)


if __name__ == '__main__':
    unittest.main()
