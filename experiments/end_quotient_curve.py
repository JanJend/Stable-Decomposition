#!/usr/bin/env python3
"""Sweep epsilon from 0% to 5% in 0.1% steps on three 320- or 640-point examples."""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import shutil
import subprocess

from end_quotient_size import DIMENSIONS, FAMILIES, MEMORY_MB, TIMEOUT_SECONDS, measure, parse_output, select_inputs
from support import DATA, HERE, ROOT, input_lock
from support.runtime import digest
from support.scaling import presentation_info

EPSILON_PERCENTAGES = [step / 10 for step in range(51)]
CSV_COLUMNS = ("family", "dataset", "n_points", "homology", "epsilon_percent", "epsilon_x", "epsilon_y",
               "End_0_dimension", "End_2epsilon_dimension", "quotient_dimension", "status", "error")


def write_csv(path, inputs, observations):
    measured = {(row["family"], row["epsilon_percent"]): row for row in observations}
    temporary = path.with_suffix(".csv.tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for sample in inputs:
            for percentage in EPSILON_PERCENTAGES:
                row = measured.get((sample["family"], percentage), {})
                epsilon = [side * percentage / 100 for side in sample["bounding_box_sides"]]
                writer.writerow(dict(
                    **{key: sample[key] for key in CSV_COLUMNS[:4]}, epsilon_percent=f"{percentage:.1f}",
                    epsilon_x=epsilon[0], epsilon_y=epsilon[1],
                    End_0_dimension=row.get("End_0_dimension", ""),
                    End_2epsilon_dimension=row.get("End_d_basis_I_size", ""),
                    quotient_dimension=row.get("quotient_basis_I_size", ""),
                    status=row.get("status", "pending"), error=row.get("error", "")))
    temporary.replace(path)


def plot_csv(path):
    """Open one figure with two measured curves per dataset; leave saving to the user."""
    import matplotlib.pyplot as plt

    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    figure, axes = plt.subplots(1, 3, figsize=(11, 3.6))
    titles = ("Noisy Torus", "Noisy Clusters", "Immune Cells")
    for axis, (family, _, homology), title in zip(axes, FAMILIES, titles):
        values = sorted((row for row in rows if row["family"] == family),
                        key=lambda row: float(row["epsilon_percent"]))
        x = [float(row["epsilon_percent"]) for row in values]
        for column, legend, style in (
            ("End_2epsilon_dimension", r"$\dim\,\mathrm{End}_{2\varepsilon}(X)$", "-"),
            ("quotient_dimension", r"$\dim\,\mathrm{End}_{2\varepsilon}(X)/(\eta\mathrm{End}_0(X))$", "--"),
        ):
            y = [float(row[column]) if row[column] else math.nan for row in values]
            axis.plot(x, y, style, label=legend, linewidth=1.6)
        axis.set(title=fr"{title} $H_{homology}$", xlabel=r"$\varepsilon$ (%)", xlim=(0, 5))
        axis.set_ylim(bottom=0)
        axis.spines[["top", "right"]].set_visible(False)
        axis.ticklabel_format(axis="y", style="plain", useOffset=False)
    figure.legend(*axes[0].get_legend_handles_labels(), loc="upper center", ncol=2, frameon=False)
    figure.tight_layout(rect=(0, 0, 1, 0.85))
    plt.show()
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=int, choices=(320, 640), default=320,
                        help="Target point count for all three inputs (default: 320)")
    parser.add_argument("--manifest", type=Path, default=DATA / "manifest.csv")
    parser.add_argument("--worker", type=Path, default=ROOT / "build/local-release/hom_reduction_benchmark")
    parser.add_argument("--output-dir", type=Path,
                        help="Results directory (default: results/end_quotient_curve or end_quotient_curve_640)")
    parser.add_argument("--timeout", type=float, default=TIMEOUT_SECONDS, help="Seconds per input/epsilon (default: 600)")
    parser.add_argument("--memory-limit-mb", type=float, default=MEMORY_MB)
    plotting = parser.add_mutually_exclusive_group()
    plotting.add_argument("--no-plot", action="store_true", help="Compute and save CSV without opening a figure")
    plotting.add_argument("--plot-only", action="store_true", help="Open the existing CSV without recomputing")
    args = parser.parse_args()
    directory = (args.output_dir or HERE / "results" / (
        "end_quotient_curve" if args.points == 320 else f"end_quotient_curve_{args.points}")).resolve()
    csv_path = directory / "results.csv"
    if args.plot_only:
        if not csv_path.is_file():
            parser.error(f"Missing results: {csv_path}")
        plot_csv(csv_path)
        return 0
    if (not math.isfinite(args.timeout) or args.timeout <= 0
            or not math.isfinite(args.memory_limit_mb) or args.memory_limit_mb < 0):
        parser.error("Timeout must be finite and positive; memory must be finite and nonnegative")
    if not args.no_plot:
        try:
            import matplotlib.pyplot  # Check before starting the computation.
        except ImportError:
            parser.error("Plotting requires Matplotlib; use a Python environment containing it or pass --no-plot")
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
            inputs = select_inputs(args.manifest, point_counts=[args.points])
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
        run = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
        logs = directory / run
        logs.mkdir(parents=True)
        (logs / "settings.json").write_text(json.dumps(dict(
            inputs=inputs, manifest=str(args.manifest), manifest_sha256=digest(args.manifest),
            worker=str(args.worker), worker_sha256=digest(args.worker), hom_algorithm="optimised", dimensions_only=True,
            points=[args.points], percentages=EPSILON_PERCENTAGES, synthetic_seed=1,
            epsilon_convention="Each generator/relation bounding-box side times epsilon_percent/100; shift = 2*epsilon",
            quotient="End_(2 epsilon)(X) / (eta End_0(X))",
            timeout_seconds=args.timeout, memory_limit_mb=args.memory_limit_mb), indent=2) + "\n")
        if csv_path.exists():
            shutil.copyfile(csv_path, logs / "previous_results.csv")
        observations = []

        def report():
            write_csv(logs / "results.csv", inputs, observations)
            shutil.copyfile(logs / "results.csv", csv_path)

        report()
        failed = False
        print(f"End quotient curves: 3 inputs x 51 epsilons; CSV: {csv_path}", flush=True)
        try:
            with (logs / "observations.jsonl").open("a") as records:
                for sample in inputs:
                    path = Path(sample["input"])
                    for percentage in EPSILON_PERCENTAGES:
                        if digest(path) != sample["input_sha256"]:
                            raise SystemExit(f"Input changed during the experiment: {path}; saved results are in {logs}")
                        prefix = logs / f"{path.stem}_epsilon_{percentage:.1f}pct"
                        stdout, stderr = Path(str(prefix) + ".stdout.txt"), Path(str(prefix) + ".stderr.txt")
                        command = [str(args.worker), str(path), f"{percentage:.1f}", "--percent", "--dimensions-only",
                                   "--memory-limit-mb", str(args.memory_limit_mb)]
                        print(f"RUN {sample['family']}: epsilon={percentage:.1f}%", flush=True)
                        status, seconds, error = measure(command, stdout, stderr, args.timeout)
                        row = dict(family=sample["family"], dataset=sample["dataset"], epsilon_percent=percentage)
                        try:
                            row.update(parse_output(stdout.read_text()))
                        except (ValueError, KeyError) as parse_error:
                            status, error = "error", f"Cannot parse worker output: {parse_error}"
                        if status == "ok" and (row.get("hom_algorithm") != "optimised"
                                                or any(key not in row for _, _, key in DIMENSIONS)):
                            status, error = "error", "Worker output is incomplete or uses the wrong Hom algorithm"
                        row.update(status=status, process_seconds=seconds, error=error, command=command)
                        records.write(json.dumps(row) + "\n")
                        records.flush()
                        observations.append(row)
                        report()
                        failed |= status != "ok"
                        print(f"  {status}: {seconds:.2f}s; End_2eps={row.get('End_d_basis_I_size', '?')}, "
                              f"quotient={row.get('quotient_basis_I_size', '?')}", flush=True)
        except KeyboardInterrupt:
            print(f"Stopped. Partial CSV: {csv_path}", flush=True)
            return 130
    print(f"Results: {csv_path}", flush=True)
    if not args.no_plot:
        plot_csv(csv_path)
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
