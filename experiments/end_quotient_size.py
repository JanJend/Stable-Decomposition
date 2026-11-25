#!/usr/bin/env python3
"""Measure only End_0, End_(2 epsilon), and quotient dimensions in three size tables."""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import shutil
import subprocess

from benchmark_end_sizes import EPSILON_PERCENTAGES, MEMORY_MB, TIMEOUT_SECONDS, measure, parse_output
from support import DATA, HERE, ROOT, input_lock, select_samples
from support.paper_tables import ROW_END, compact_cell, latex
from support.runtime import digest
from support.scaling import presentation_info

POINT_COUNTS = [320, 640, 1280]
FAMILIES = [("torus_gauss", "Gaussian torus", "1"),
            ("gauss_clust", "Clusters", "0"), ("immuno", "Immuno", "1")]
DIMENSIONS = [("dim End_0", r"$\dim\mathrm{End}_0(X)$", "End_0_dimension"),
              ("dim End_(2 epsilon)", r"$\dim\mathrm{End}_{2\varepsilon}(X)$", "End_d_basis_I_size"),
              ("dim Q", r"$\dim Q$", "quotient_basis_I_size")]
REPORT_FILES = ["results.tex", *(f"results_{n}.{suffix}" for n in POINT_COUNTS for suffix in ("csv", "tex"))]


def select_inputs(manifest, *, point_counts=POINT_COUNTS):
    with manifest.open(newline="") as stream:
        selected = select_samples(list(csv.DictReader(stream)), point_counts)
    inputs = []
    for n in point_counts:
        for family, _, homology in FAMILIES:
            matches = [row for row in selected if row["family"] == family
                       and int(row["target_points"]) == n and row["homology"] == homology
                       and (family == "immuno" or row["seed"] == "1")]
            if len(matches) != 1:
                raise ValueError(f"Expected one {family} H{homology} sample at target size {n} (synthetic sample 1)")
            inputs.append(dict(matches[0], input=str(manifest.parent / matches[0]["scc"])))
    return inputs


def write_reports(directory, observations, inputs):
    measured = {(int(row["target_points"]), row["family"], row["epsilon_percent"]): row
                for row in observations}
    tables = []
    for n in POINT_COUNTS:
        columns = [(family, percentage) for family, _, _ in FAMILIES for percentage in EPSILON_PERCENTAGES]
        rows = []
        for label, tex_label, key in DIMENSIONS:
            values = [measured.get((n, family, percentage), {}) for family, percentage in columns]
            rows.append((label, tex_label, [row.get(key, row.get("status", "pending")) for row in values]))
        with (directory / f"results_{n}.csv").open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["Dimension", *(f"{family} epsilon={percentage:g}%" for family, percentage in columns)])
            writer.writerows([label, *values] for label, _, values in rows)
        immune = next(row for row in inputs if int(row["target_points"]) == n and row["family"] == "immuno")
        caption = (f"End and quotient dimensions at target size $n={n}$. "
                   r"$Q=\mathrm{End}_{2\varepsilon}(X)/(\eta\mathrm{End}_0(X))$. "
                   r"Epsilon scales each bounding-box side; synthetic sample 1. "
                   f"The immune cloud has {immune['n_points']} cells. "
                   r"-- means pending; timeout/mem./error mark incomplete measurements.")
        lines = [r"\begin{table}[tbp]", r"\centering", r"\begingroup\scriptsize",
                 r"\setlength{\tabcolsep}{3pt}", r"\renewcommand{\arraystretch}{1.05}",
                 r"\caption{" + caption + "}", r"\begin{tabular}{@{}l*{9}{r}@{}}", r"\hline",
                 " & " + " & ".join(r"\multicolumn{3}{c}{" + label + "}" for _, label, _ in FAMILIES) + ROW_END,
                 r"$\varepsilon$ (\%) & " + " & ".join(f"{percentage:g}" for _, percentage in columns) + ROW_END,
                 r"\hline"]
        for _, label, values in rows:
            cells = [str(value) if isinstance(value, int) else compact_cell(value) for value in values]
            lines.append(label + " & " + " & ".join(latex(cell) for cell in cells) + ROW_END)
        lines += [r"\hline", r"\end{tabular}", r"\endgroup", r"\end{table}"]
        table = "\n".join(lines) + "\n"
        (directory / f"results_{n}.tex").write_text(table)
        tables.append(table)
    (directory / "results.tex").write_text("% Generated automatically; no seed averaging.\n" + "\n".join(tables))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DATA / "manifest.csv")
    parser.add_argument("--worker", type=Path, default=ROOT / "build/local-release/hom_reduction_benchmark")
    parser.add_argument("--output-dir", type=Path, default=HERE / "results/end_quotient_size")
    parser.add_argument("--timeout", type=float, default=TIMEOUT_SECONDS, help="Seconds per input/epsilon (default: 600)")
    parser.add_argument("--memory-limit-mb", type=float, default=MEMORY_MB)
    args = parser.parse_args()
    if (not math.isfinite(args.timeout) or args.timeout <= 0
            or not math.isfinite(args.memory_limit_mb) or args.memory_limit_mb < 0):
        parser.error("Timeout must be finite and positive; memory must be finite and nonnegative")
    args.worker, args.manifest = args.worker.resolve(), args.manifest.resolve()
    if not args.worker.is_file():
        parser.error("Build the worker first: cmake --build build/local-release --target hom_reduction_benchmark")
    help_text = subprocess.run([str(args.worker), "--help"], check=True, capture_output=True, text=True, timeout=10).stdout
    if "hom_algorithm=optimised" not in help_text or "--dimensions-only" not in help_text:
        parser.error("Rebuild hom_reduction_benchmark: this experiment requires A-1/2 and --dimensions-only")
    if not args.manifest.is_file():
        parser.error(f"Prepare the shared inputs first: missing {args.manifest}")
    with input_lock(args.manifest.parent):
        try:
            inputs = select_inputs(args.manifest)
        except ValueError as error:
            parser.error(str(error))
        for sample in inputs:
            path = Path(sample["input"])
            if not path.is_file() or sample.get("generation_status", "ok") != "ok":
                parser.error(f"Missing or ungenerated SCC input: {path}")
            sample["input_sha256"] = digest(path)
            if sample.get("scc_sha256", sample["input_sha256"]) != sample["input_sha256"]:
                parser.error(f"SCC differs from the shared manifest: {path}")
            sample["generators"], sample["relations"], sample["bounding_box_sides"] = presentation_info(path, per_axis=True)
            sample["bounding_box_extent"] = max(sample["bounding_box_sides"])
        worker_hash = digest(args.worker)
        run = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        directory = args.output_dir.resolve()
        logs = directory / run
        logs.mkdir(parents=True)
        (logs / "settings.json").write_text(json.dumps(dict(
            inputs=inputs, manifest=str(args.manifest), manifest_sha256=digest(args.manifest),
            worker=str(args.worker), worker_sha256=worker_hash, hom_algorithm="optimised", dimensions_only=True,
            points=POINT_COUNTS, families=[family for family, _, _ in FAMILIES], synthetic_seed=1,
            percentages=EPSILON_PERCENTAGES, timeout_seconds=args.timeout, memory_limit_mb=args.memory_limit_mb), indent=2) + "\n")
        for name in REPORT_FILES:
            if (directory / name).exists():
                shutil.copyfile(directory / name, logs / ("previous_" + name))
        observations = []

        def report():
            write_reports(logs, observations, inputs)
            for name in REPORT_FILES:
                shutil.copyfile(logs / name, directory / name)

        report()
        failed = False
        print(f"End quotient sizes: {len(inputs)} inputs x {len(EPSILON_PERCENTAGES)} epsilons; results: {logs}", flush=True)
        try:
            with (logs / "observations.jsonl").open("a") as records:
                for sample in inputs:
                    for percentage in EPSILON_PERCENTAGES:
                        path = Path(sample["input"])
                        if digest(path) != sample["input_sha256"]:
                            raise SystemExit(f"Input changed during the experiment: {path}; saved results are in {logs}")
                        prefix = logs / f"{path.stem}_epsilon_{percentage:g}pct"
                        stdout_path, stderr_path = Path(str(prefix) + ".stdout.txt"), Path(str(prefix) + ".stderr.txt")
                        command = [str(args.worker), str(path), str(percentage), "--percent", "--dimensions-only",
                                   "--memory-limit-mb", str(args.memory_limit_mb)]
                        epsilon = tuple(side * (percentage / 100) for side in sample["bounding_box_sides"])
                        row = dict(sample, run=run, worker_sha256=worker_hash, epsilon_percent=percentage,
                                   epsilon=epsilon, shift=tuple(2 * value for value in epsilon),
                                   timeout_seconds=args.timeout,
                                   stdout_log=str(stdout_path.relative_to(directory)), stderr_log=str(stderr_path.relative_to(directory)))
                        print(f"RUN {path.stem}: epsilon={percentage:g}% ({row['epsilon']})", flush=True)
                        status, seconds, error = measure(command, stdout_path, stderr_path, args.timeout)
                        try:
                            row.update(parse_output(stdout_path.read_text()))
                        except (ValueError, KeyError) as parse_error:
                            status, error = "error", f"Cannot parse worker output: {parse_error}"
                        if status == "ok" and (row.get("hom_algorithm") != "optimised"
                                                or any(key not in row for _, _, key in DIMENSIONS)):
                            status, error = "error", "Worker output is incomplete or uses the wrong Hom algorithm"
                        row.update(status=status, process_seconds=seconds, error=error, command=command)
                        Path(str(prefix) + ".json").write_text(json.dumps(row, indent=2) + "\n")
                        records.write(json.dumps(row, separators=(",", ":")) + "\n")
                        records.flush()
                        observations.append(row)
                        report()
                        failed |= status != "ok"
                        print(f"  {status}: {seconds:.2f}s; " + ", ".join(
                            f"{label}={row.get(key, '?')}" for label, _, key in DIMENSIONS), flush=True)
        except KeyboardInterrupt:
            print(f"Stopped. Completed measurements and partial tables are saved in {logs}", flush=True)
            return 130
        print(f"Results: {directory / 'results.tex'} and results_{{320,640,1280}}.csv", flush=True)
        return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
