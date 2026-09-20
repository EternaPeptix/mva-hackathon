#!/usr/bin/env python3
"""
Characterise the Track 1 scoring function empirically.

Imports the organisers' evaluation.py UNMODIFIED and probes its behaviour. Every table in
report/scoring-analysis.md is produced here, so the claims in that document are reproducible
rather than asserted.

Usage:
    python3 score_probe.py /path/to/space/checkout

where the checkout is the challenge Space source (containing evaluation.py). Fetch it with:
    for f in evaluation.py groundtruth.py config.py; do
      curl -sLO "https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/raw/main/$f"
    done
"""
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()))
import evaluation as E  # noqa: E402

# The two reported alleles, in the naming convention the scorer expects.
A = ("chr15", 40209701, "T", "G")
B = ("chr15", 40220612, "T", "G")
TRUE = frozenset([A, B])
JUNK = [("chr%d" % i, 1000 * i, "A", "T") for i in range(1, 12)]


def score(rows, truth=TRUE):
    """rows: [(variant_iterable, epcr)] already in submission order."""
    srs = [
        E.SubmissionRow(variants=frozenset(v), epcr=e, rank=i + 1)
        for i, (v, e) in enumerate(rows)
    ]
    return E.score_proband("P", srs, truth)


def hr(title):
    print("\n" + title)
    print("-" * len(title))


hr("1. Non-monotonicity: does adding the CORRECT row ever LOWER rank points?")
print("   Rank-1 row holds a partial match; the correct pair is inserted at rank f.")
print("   f | without correct row | with correct row | delta")
for f in range(2, 11):
    base = [(frozenset([A, JUNK[0]]), 0.99)]
    base += [(frozenset([JUNK[k]]), 0.9 - 0.01 * k) for k in range(1, f - 1)]
    without = score(base).rank_points
    with_ = score(base + [(TRUE, 0.9 - 0.01 * (f - 1))]).rank_points
    delta = with_ - without
    flag = "  <-- PENALISED" if delta < 0 else ("  (tie)" if delta == 0 else "")
    print(f"  {f:2d} | {without:18.1f} | {with_:16.1f} | {delta:+6.1f}{flag}")

hr("2. Duplicate EPCR values cap F-max")
for tie in (False, True):
    e2 = 0.95 if tie else 0.80
    r = score([(TRUE, 0.95), (frozenset([JUNK[0]]), e2)])
    print(f"  junk EPCR={e2}  tie={tie!s:5}  pts={r.rank_points:5.1f}  f_max={r.f_max:.4f}")

hr("3. With a rank-1 correct pair and distinct EPCRs, extra rows are free")
for n in (0, 1, 3, 5, 9):
    rows = [(TRUE, 0.99)] + [(frozenset([JUNK[k]]), 0.90 - 0.05 * k) for k in range(n)]
    r = score(rows)
    print(f"  {n:2d} extra rows -> pts={r.rank_points:5.1f}  f_max={r.f_max:.4f}")

hr("4. Chromosome-naming hedge payoff matrix")
A2, B2 = ("15", 40209701, "T", "G"), ("15", 40220612, "T", "G")
TRUE_BARE = frozenset([A2, B2])
for gt_label, gt in (("chr-prefixed", TRUE), ("bare", TRUE_BARE)):
    for label, rows in (
        ("chr only", [(TRUE, 0.99)]),
        ("chr@1 + bare@2", [(TRUE, 0.99), (TRUE_BARE, 0.80)]),
    ):
        r = score(rows, gt)
        print(f"  GT {gt_label:12} | {label:15} -> pts={r.rank_points:5.1f} f_max={r.f_max:.4f}")

hr("5. Tier tables (the source of finding 1)")
ranks = range(1, 11)
print("  rank:       " + "".join(f"{i:6d}" for i in ranks))
print("  full:       " + "".join(f"{E._rank_to_points(i):6.0f}" for i in ranks))
print("  partial(1/2):" + "".join(f"{0.5 * E._rank_to_points(i):6.1f}" for i in ranks))
