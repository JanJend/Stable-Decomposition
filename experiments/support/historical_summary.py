#!/usr/bin/env python3
"""Calculate a scaling table from measured runs, preserving censored results."""
import argparse
import json
import math
from pathlib import Path

from . import HERE
from .scaling import SWITCHES, speedup, summarize


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=HERE / "optimization_scaling/torus_gaussian_h1_eps1.jsonl")
    parser.add_argument("--output", type=Path, default=HERE / "optimization_scaling/torus_gaussian_h1_eps1_summary.json")
    args = parser.parse_args()
    records = [json.loads(line) for line in args.raw.read_text().splitlines()]
    hashes = {row["worker_sha256"] for row in records}
    if len(hashes) != 1:
        raise ValueError("Choose observations from a single worker binary")
    rows = []
    for points in sorted({record["n_points"] for record in records}):
        current = [record for record in records if record["n_points"] == points]
        baseline = summarize([r for r in current if r["mask"] == 0])
        variants = []
        signatures = set()
        if baseline["status"] == "ok":
            signatures.add(tuple(baseline["signature"]))
        for bit, name in enumerate(SWITCHES):
            result = summarize([r for r in current if r["mask"] == 1 << bit])
            if result["status"] == "ok":
                signatures.add(tuple(result["signature"]))
            variants.append(dict(name=name, **result, speedup=speedup(baseline, result)))
        if len(signatures) > 1:
            raise ValueError(f"Optimization output/iteration metadata differs at {points} points")
        meta = current[0]
        rows.append(dict(points=points, epsilon=meta["epsilon"], generators=meta["generators"],
                         relations=meta["relations"], input_sha256=meta["input_sha256"],
                         baseline=baseline, variants=variants))
    result = dict(family=records[0]["family"], homology=records[0]["homology"], seed=records[0]["seed"],
                  epsilon_percent=records[0]["epsilon_percent"], worker_sha256=next(iter(hashes)),
                  speedup_definition="100*(baseline_seconds/variant_seconds-1)", rows=rows)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print("Points | " + " | ".join(SWITCHES))
    for row in rows:
        values = []
        for variant in row["variants"]:
            gain = variant["speedup"]
            if gain["kind"] == "unavailable":
                values.append("n/a")
            else:
                prefix = {"exact": "", "lower_bound": ">", "upper_bound": "<"}[gain["kind"]]
                percent = gain["percent"]
                if gain["kind"] == "lower_bound":
                    percent = math.floor(percent * 10) / 10
                elif gain["kind"] == "upper_bound":
                    percent = math.ceil(percent * 10) / 10
                values.append(f"{prefix}{percent:,.1f}%" if prefix else f"{percent:+,.1f}%")
        print(str(row["points"]) + " | " + " | ".join(values))


if __name__ == "__main__":
    main()
