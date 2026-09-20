#!/usr/bin/env python3
"""
Audit the public Track 1 leaderboard: is its visible ranking carrying information?

Motivation. On 2026-08-26, every entry on the Track 1 leaderboard held the maximum
score (rank_points 100.0, f_max 1.000). When a metric saturates, the *order* the
board displays is decided entirely by its tie-break. This script fetches the live
board and tests, against real data, which of three candidate tie-breaks explains the
observed order: submission time, alphabetical username, or nothing (random).

Reading the Space source predicts the answer. utils.best_per_team ends with

    sorted_rows = sorted(best.values(), key=lambda r: (-r["rank_points"], -r["f_max"]))

Python's sorted is stable, so equal-scoring rows keep the order they had in
best.values(). That is dict insertion order, which is the order rows arrive from
load_leaderboard(), which iterates api.list_repo_files(...) — repository paths of
the form  track1_submissions/<hf_username>/sub_<timestamp>.json  in lexicographic
order. Lexicographic by path means grouped by username, in *byte* order: digits
before uppercase before lowercase.

So the prediction is: the board is sorted by ASCII username, and submission time is
irrelevant. This script checks that prediction rather than asserting it.

Usage:  python3 leaderboard_audit.py [--json out.json]
Requires no credentials; the leaderboard endpoint is public.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from typing import Any

SPACE = "https://sagebio-rare-disease-real-kid-mva-hackathon-2026.hf.space"


def fetch_leaderboard() -> list[list[Any]]:
    """Call the Space's Gradio endpoint and return its raw table rows."""
    start = subprocess.run(
        ["curl", "-sS", "-m", "60", "-X", "POST", f"{SPACE}/gradio_api/call/leaderboard_df",
         "-H", "Content-Type: application/json", "-d", '{"data":[]}'],
        capture_output=True, text=True, check=True,
    )
    event_id = json.loads(start.stdout)["event_id"]
    stream = subprocess.run(
        ["curl", "-sS", "-m", "90", "-N", f"{SPACE}/gradio_api/call/leaderboard_df/{event_id}"],
        capture_output=True, text=True, check=True,
    )
    payloads = [ln[5:] for ln in stream.stdout.splitlines() if ln.startswith("data:")]
    if not payloads:
        sys.exit("no data frame returned by the Space")
    return json.loads(payloads[-1])[0]["data"]


def username_of(participant: str) -> str:
    """Rows render as 'Display (@user)' or '@user'; the sort key is the username."""
    m = re.search(r"@([^)\s]+)", participant)
    return m.group(1) if m else participant


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", help="write the parsed board and verdicts here")
    args = ap.parse_args()

    rows = fetch_leaderboard()
    board = [
        {"place": r[0], "user": username_of(r[1]), "rank_points": float(r[2]),
         "f_max": float(r[3]), "file": r[4], "submitted_at": r[5]}
        for r in rows
    ]

    top = max(e["rank_points"] for e in board)
    tied = [e for e in board if e["rank_points"] == top
            and e["f_max"] == max(x["f_max"] for x in board if x["rank_points"] == top)]

    users = [e["user"] for e in tied]
    times = [e["submitted_at"] for e in tied]

    # Three candidate explanations for the displayed order of the tied block.
    verdicts = {
        "ascii_username":      users == sorted(users),                 # byte order
        "casefold_username":   [u.lower() for u in users] == sorted(u.lower() for u in users),
        "oldest_first":        times == sorted(times),
        "newest_first":        times == sorted(times, reverse=True),
    }

    print(f"leaderboard rows            : {len(board)}")
    print(f"tied at the maximum score   : {len(tied)}  "
          f"(rank_points={top}, f_max={tied[0]['f_max'] if tied else float('nan')})")
    print(f"saturated                   : {len(tied) == len(board)}")
    print()
    for name, ok in verdicts.items():
        print(f"  tie order explained by {name:<20} {ok}")
    print()

    if verdicts["ascii_username"] and not (verdicts["oldest_first"] or verdicts["newest_first"]):
        print("CONCLUSION: the displayed ranking of the tied block is alphabetical by")
        print("username in byte order, and is unrelated to submission time or to any")
        print("property of the submission. With the metric saturated, the medal")
        print("positions carry no information about submission quality.")

    print("\n  #   username              submitted            file")
    for e in board:
        print(f"  {str(e['place']):<3} {e['user']:<21} {e['submitted_at']:<20} {e['file'][:44]}")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"board": board, "verdicts": verdicts,
                       "n_tied": len(tied), "saturated": len(tied) == len(board)}, fh, indent=2)
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
