"""Exact bounded full-query fitting; sequence compatibility is not locus identity.

The query is global and target flanks are free. Every query base is consumed;
query-end insertions cost one per base. CIGAR uses =/X/I/D in query-to-target
orientation. Coordinates are zero-based, half-open on the supplied target.
Only nonempty uppercase ACGT strings are accepted. Qualities, orientation,
physical-end selection and fragment aggregation belong to the caller.
"""
from array import array


def _positive_integer(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(name + ' must be a positive integer')


def _cigar(reverse_operations):
    pieces = []
    previous, count = None, 0
    for operation in reversed(reverse_operations):
        if operation == previous:
            count += 1
        else:
            if previous is not None:
                pieces.append(str(count) + previous)
            previous, count = operation, 1
    if previous is not None:
        pieces.append(str(count) + previous)
    return ''.join(pieces)


def fit_full_query(query, target, max_cells=1000000, max_paths=256):
    """Return exact unit edit cost/endpoints and bounded optimal CIGAR witnesses.

    All optimal target endpoints are exact, even when traceback is capped. The
    optimal path count is saturated at max_paths+1; path_count_is_exact is false
    when that sentinel is reached. Each witness is a distinct (start0,end0,CIGAR)
    alignment. A capped result does not imply that the listed witnesses exhaust
    the possible starts or allele coverage. Traversal is iterative and bounded
    by max_paths times alignment length, rather than an unbounded path search.

    Exact substrings take a proven zero-cost fast path: every overlapping exact
    occurrence is counted, because a zero-cost fit can contain only matches.
    optimal_path_count_exact reports that total even when witnesses are capped;
    it is otherwise null if the DP count saturated. The capped-count fields keep
    their original semantics. All cell/path bounds are enforced before this path.

    Storage is O((len(query)+1)*(len(target)+1)) for unsigned DP costs plus two
    count rows and one traceback stack. The caller-supplied cell cap is checked
    before allocation. A hard 10-million-cell ceiling prevents accidental huge
    allocations; max_paths has a hard 4096 ceiling. Empty target substrings are
    legal placements and consume the complete nonempty query as insertions.
    """
    for sequence, name in ((query, 'query'), (target, 'target')):
        if not isinstance(sequence, str) or not sequence or set(sequence) - set('ACGT'):
            raise ValueError(name + ' must be nonempty uppercase ACGT')
    _positive_integer(max_cells, 'max_cells')
    _positive_integer(max_paths, 'max_paths')
    if max_cells > 10000000:
        raise ValueError('max_cells exceeds hard ceiling 10000000')
    if max_paths > 4096:
        raise ValueError('max_paths exceeds hard ceiling 4096')
    m, n = len(query), len(target)
    width = n + 1
    cells = (m + 1) * width
    if cells > max_cells:
        raise ValueError('full-query DP cell cap exceeded')
    exact_starts, search_start = [], 0
    while True:
        occurrence = target.find(query, search_start)
        if occurrence < 0:
            break
        exact_starts.append(occurrence)
        search_start = occurrence + 1  # Preserve overlapping repeat placements.
    if exact_starts:
        exact_count = len(exact_starts)
        saturated = exact_count > max_paths
        return dict(best_cost=0, optimal_end0=[start + m for start in exact_starts],
                    optimal_path_count_capped=min(exact_count, max_paths + 1),
                    optimal_path_count_exact=exact_count,
                    path_count_is_exact=not saturated, traceback_complete=not saturated,
                    saturated=saturated, status='capped' if saturated else 'complete',
                    unique_optimal_alignment=exact_count == 1,
                    alignments=[dict(start0=start, end0=start + m, cigar=str(m) + '=',
                                     edit_cost=0, query_start0=0, query_end0=m)
                                for start in exact_starts[:max_paths]],
                    query_length=m, target_length=n, dp_cells=cells,
                    dp_cells_evaluated=0, max_paths=max_paths, traceback_steps=0,
                    exact_substring_fast_path=True,
                    score_model='unit_edit_query_global_target_local',
                    alphabet_policy='nonempty_uppercase_ACGT_only')
    costs = array('I', [0]) * cells
    count_limit = max_paths + 1
    previous_counts = [1] * width
    for i in range(1, m + 1):
        base, above = i * width, (i - 1) * width
        costs[base] = i
        counts = [1] + [0] * n
        for j in range(1, n + 1):
            diagonal = costs[above + j - 1] + (query[i - 1] != target[j - 1])
            insertion = costs[above + j] + 1
            deletion = costs[base + j - 1] + 1
            best = min(diagonal, insertion, deletion)
            costs[base + j] = best
            count = 0
            if diagonal == best:
                count += previous_counts[j - 1]
            if insertion == best:
                count += previous_counts[j]
            if deletion == best:
                count += counts[j - 1]
            counts[j] = min(count, count_limit)
        previous_counts = counts
    final = m * width
    best_cost = min(costs[final:final + width])
    ends = [j for j in range(width) if costs[final + j] == best_cost]
    path_count = min(count_limit, sum(previous_counts[j] for j in ends))
    saturated = path_count == count_limit

    alignments, seen = [], set()
    traceback_steps = 0
    for end0 in ends:
        # Frames hold query index, target index and next predecessor choice.
        # The shared operation list avoids quadratic copying along a long path.
        stack, reverse_operations = [(m, end0, 0)], []
        while stack and len(alignments) < max_paths:
            traceback_steps += 1
            i, j, choice = stack[-1]
            if i == 0:
                cigar = _cigar(reverse_operations)
                identity = (j, end0, cigar)
                if identity in seen:
                    raise RuntimeError('duplicate optimal alignment witness')
                seen.add(identity)
                alignments.append(dict(start0=j, end0=end0, cigar=cigar,
                                       edit_cost=best_cost, query_start0=0, query_end0=m))
                stack.pop()
                if reverse_operations:
                    reverse_operations.pop()
                continue
            if choice == 3:
                stack.pop()
                if reverse_operations:
                    reverse_operations.pop()
                continue
            stack[-1] = (i, j, choice + 1)
            value = costs[i * width + j]
            predecessor = None
            if choice == 0 and j:
                mismatch = query[i - 1] != target[j - 1]
                if costs[(i - 1) * width + j - 1] + mismatch == value:
                    predecessor = (i - 1, j - 1, 'X' if mismatch else '=')
            elif choice == 1:
                if costs[(i - 1) * width + j] + 1 == value:
                    predecessor = (i - 1, j, 'I')
            elif choice == 2 and j:
                if costs[i * width + j - 1] + 1 == value:
                    predecessor = (i, j - 1, 'D')
            if predecessor is not None:
                next_i, next_j, operation = predecessor
                reverse_operations.append(operation)
                stack.append((next_i, next_j, 0))
        if len(alignments) == max_paths:
            break
    if len(alignments) != min(path_count, max_paths):
        raise RuntimeError('optimal path count and witnesses disagree')
    return dict(best_cost=best_cost, optimal_end0=ends,
                optimal_path_count_capped=path_count, path_count_is_exact=not saturated,
                optimal_path_count_exact=None if saturated else path_count,
                traceback_complete=not saturated, saturated=saturated,
                status='capped' if saturated else 'complete',
                unique_optimal_alignment=path_count == 1,
                alignments=alignments, query_length=m, target_length=n,
                dp_cells=cells, dp_cells_evaluated=cells, max_paths=max_paths,
                traceback_steps=traceback_steps, exact_substring_fast_path=False,
                score_model='unit_edit_query_global_target_local',
                alphabet_policy='nonempty_uppercase_ACGT_only')
