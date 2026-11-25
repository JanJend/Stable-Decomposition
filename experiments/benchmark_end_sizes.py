#!/usr/bin/env python3
"""Measure the End-space reduction funnel on fixed point sets."""
import argparse
import csv
from contextlib import nullcontext
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import time

from support import DATA, HERE, ROOT, input_lock, select_samples
from support.runtime import digest
from support.scaling import presentation_info
from support.paper_tables import ROW_END, compact_cell, latex
from support.generator_options import add_generator_arguments, generator_settings, generator_flags

# Default, independent of the optimization experiment: six sample-1 clouds
# at 640/1280 target points, plus the complete LargeHypoxicRegion H1 cloud.
POINT_COUNTS = [640, 1280]
FAMILIES = ["torus_gauss", "immuno", "gauss_clust"]
EPSILON_PERCENTAGES = [0.5, 1, 2]
TIMEOUT_SECONDS = 10 * 60
LARGE_TIMEOUT_SECONDS = 60 * 60
MEMORY_MB = 12288
STAGES = ("End_d_basis", "quotient_basis", "quotient_coordinates", "composition_generators")
COLUMNS = (
    ("dim End(X)_0", "End_0_dimension"),
    ("dim End(X)_(2eps)", "End_d_basis_I_size"),
    ("dim Q", "quotient_basis_I_size"),
    ("|B_l| (I)", "composition_generators_I_size"),
    ("|B_r| (K)", "composition_generators_K_size"),
    ("I restricted max/pass", "I_restriction_max"),
    ("K restricted max/pass", "K_restriction_max"),
)
IDENTIFIERS = ("Point set", "n", "epsilon (%)")
TIMINGS = (
    ("End_0 seconds", "End_0_seconds"), ("End_2eps seconds", "End_d_basis_seconds"),
    ("Quotient seconds", "quotient_basis_seconds"), ("Coordinates seconds", "quotient_coordinates_seconds"),
    ("Composition seconds", "composition_generators_seconds"),
    ("I restriction seconds", "I_restriction_seconds"), ("K restriction seconds", "K_restriction_seconds"),
    ("I phase seconds", "I_phase_seconds"), ("K phase seconds", "K_phase_seconds"),
    ("Process seconds", "process_seconds"),
)
REPORT_FILES = ("results.csv", "results.tex", "timings.csv", "restrictions.csv")


def parse_output(output):
    """Keep completed stages even when a later stage fails or times out."""
    result = {"restrictions": []}
    for line in output.splitlines():
        values = dict(re.findall(r"(\w+)=([^\s]+)", line))
        if "End_0_dimension" in values:
            result.update(End_0_dimension=int(values["End_0_dimension"]), End_0_seconds=float(values["seconds"]))
        if "epsilon" in values:
            for key in ("epsilon", "shift"):
                value = values[key]
                result[key] = tuple(map(float, value.strip("()").split(","))) if value.startswith("(") else float(value)
        if "hom_algorithm" in values:
            result["hom_algorithm"] = values["hom_algorithm"]
        if line.startswith("restriction_summary "):
            phase = values["phase"]
            result.update({f"{phase}_restriction_max": int(values["max_retained"]),
                           f"{phase}_restriction_total": int(values["total_retained"]),
                           f"{phase}_restriction_passes": int(values["passes"]),
                           f"{phase}_phase_seconds": float(values["phase_seconds"]),
                           f"{phase}_restriction_seconds": float(values["selection_seconds"])})
        columns = line.split("\t")
        if len(columns) in (5, 6) and columns[0] == "restriction":
            restriction = dict(phase=columns[1], iteration=int(columns[2]),
                               candidates=int(columns[3]), retained=int(columns[4]))
            if len(columns) == 6:
                restriction["seconds"] = float(columns[5])
            result["restrictions"].append(restriction)
        if len(columns) == 4 and columns[0] in STAGES:
            result.update({f"{columns[0]}_I_size": int(columns[1]),
                           f"{columns[0]}_K_size": int(columns[2]),
                           f"{columns[0]}_seconds": float(columns[3])})
    return result


def measure(command, stdout_path, stderr_path, timeout):
    start = time.monotonic()
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr)
        try:
            while True:
                remaining = timeout - (time.monotonic() - start)
                if remaining <= 0:
                    process.kill()
                    process.wait()
                    return "timeout", time.monotonic() - start, f"Exceeded {timeout:g} seconds"
                try:
                    code = process.wait(timeout=min(30, remaining))
                    status = "ok" if code == 0 else "memory_limit" if code == 75 else "error"
                    return status, time.monotonic() - start, "" if code == 0 else f"Worker exited with status {code}"
                except subprocess.TimeoutExpired:
                    lines = stdout_path.read_text().splitlines()
                    print(f"  {time.monotonic()-start:.0f}s: {lines[-1] if lines else 'preparing input'}", flush=True)
        finally:
            if process.poll() is None:
                process.kill()
            process.wait()


def select_inputs(manifest, *, immuno_points=None, only_640=False):
    with manifest.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if immuno_points is not None:
        samples = [row for row in rows if row["family"] == "immuno"
                   and int(row["target_points"]) == immuno_points and row["homology"] == "1"]
        if len(samples) != 1:
            raise ValueError(f"Expected exactly one immune H1 presentation at target size {immuno_points}")
        return [dict(row, input=str(manifest.parent / row["scc"])) for row in samples]
    points = [640] if only_640 else POINT_COUNTS
    selected = select_samples(rows, points)
    samples = []
    homology = {"torus_gauss": "1", "immuno": "1", "gauss_clust": "0"}
    for family in FAMILIES:
        for target in points:
            matches = [row for row in selected if row["family"] == family
                       and int(row["target_points"]) == target and row["homology"] == homology[family]
                       and (family == "immuno" or row["seed"] == "1")]
            if len(matches) != 1:
                raise ValueError(f"Expected one {family} sample at target size {target} (synthetic sample 1)")
            samples.extend(matches)
    if not only_640:
        large = [row for row in rows if row["dataset"] == "LargeHypoxicRegion" and row["homology"] == "1"]
        if len(large) != 1:
            raise ValueError("Expected exactly one complete LargeHypoxicRegion H1 presentation")
        samples.extend(large)
    return [dict(row, input=str(manifest.parent / row["scc"])) for row in samples]


def report_rows(observations, inputs, percentages):
    measured = {(row["input"], row["epsilon_percent"]): row for row in observations}
    rows = []
    for sample in inputs:
        for percentage in percentages:
            observation = measured.get((sample["input"], percentage), {})
            row = dict(pointset="LHR" if sample["dataset"] == "LargeHypoxicRegion" else sample["family"],
                       n_points=sample["n_points"], epsilon_percent=percentage, observation=observation)
            for _, key in (*COLUMNS, *TIMINGS):
                row[key] = observation.get(key, observation.get("status", "pending"))
            rows.append(row)
    return rows


def write_reports(directory, rows, *, generators_only=False, time_steps=False):
    def identifiers(row):
        return dict(zip(IDENTIFIERS, (row["pointset"], row["n_points"], row["epsilon_percent"])))

    columns = COLUMNS[:5] if generators_only else COLUMNS
    timings = (*TIMINGS[:5], TIMINGS[-1]) if generators_only else TIMINGS
    for name, columns in (("results.csv", columns), ("timings.csv", timings)):
        with (directory / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=(*IDENTIFIERS, *(label for label, _ in columns)))
            writer.writeheader()
            for row in rows:
                writer.writerow(dict(identifiers(row), **{label: row[key] for label, key in columns}))
    with (directory / "restrictions.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=(*IDENTIFIERS, "phase", "iteration", "candidates", "retained",
                                                   *(("seconds",) if time_steps else ())))
        writer.writeheader()
        for row in rows:
            for restriction in row["observation"].get("restrictions", []):
                writer.writerow(dict(identifiers(row), **restriction))
    if generators_only:
        tables = []
        for name, columns, headings, caption in (
            ("dimensions.tex", COLUMNS[:3], [r"$\dim E$", r"$\dim H$", r"$\dim Q$"],
             r"Dimensions over $\mathbf{F}_2$: $E=\mathrm{End}_0(X)$, "
             r"$H=\mathrm{End}_{2\varepsilon}(X)$, and $Q=H/(\eta E)$, "
             r"where $\eta:X\to X(2\varepsilon)$ is the canonical shift. "
             r"Tori and immune samples use $H_1$; clusters use $H_0$. "
             r"Synthetic sample 1; epsilon scales each presentation bounding-box side."),
            ("generators.tex", COLUMNS[3:5], [r"$|B_l|$", r"$|B_r|$"],
             r"Left and right $E$-generating families of $Q$ for the same inputs and epsilons. "
             r"Nilpotent reduction is followed by Complement-orbit with the recorded search budgets. "
             r"These are certified generating families; their sizes need not be minimal."),
        ):
            lines = [r"\begin{table}[tbp]", r"\centering", r"\small",
                     r"\caption{" + caption + "}",
                     r"\label{tab:end-" + name.removesuffix(".tex") + "}",
                     r"\begin{tabular}{@{}l" + "r" * (2 + len(columns)) + "@{}}", r"\hline",
                     " & ".join(["Point set", "$n$", r"$\varepsilon$ (\%)", *headings]) + ROW_END, r"\hline"]
            for row in rows:
                pointset = {"torus_gauss": "Gaussian torus", "gauss_clust": "Gaussian clusters",
                            "immuno": "Immune cells"}.get(row["pointset"], row["pointset"])
                values = [pointset, row["n_points"] or "--", f"{row['epsilon_percent']:g}"]
                values += [str(row[key]) if isinstance(row[key], int) else compact_cell(row[key]) for _, key in columns]
                lines.append(" & ".join(latex(str(value)) for value in values) + ROW_END)
            lines += [r"\hline", r"\end{tabular}", r"\end{table}"]
            table = "\n".join(lines) + "\n"
            (directory / name).write_text(table)
            tables.append(table)
        (directory / "results.tex").write_text("% Generated automatically; no seed averaging.\n" + "\n".join(tables))
        return
    caption = (r"End-space reduction funnel. $E=\mathrm{End}_0(X)$, "
               r"$H=\mathrm{End}_{2\varepsilon}(X)$, $Q=H/(\eta E)$. "
               r"$B_l$ and $B_r$ are the retained left/right generating families after "
               r"Nilpotent reduction followed by Complement-orbit; their columns count maps. "
               r"$r_I^{\max},r_K^{\max}$ are maximum retained restriction counts over all actual passes, "
               r"not just the terminal pass. Per-pass counts: \texttt{restrictions.csv}. "
               r"LHR = LargeHypoxicRegion. $n$ is the actual point count; synthetic sample 1. "
               r"Epsilon scales each bounding-box side by its percentage. -- means pending.")
    headings = ["Point set", "$n$", r"$\varepsilon$ (\%)", r"$\dim E$", r"$\dim H$", r"$\dim Q$",
                r"$|B_l|$", r"$|B_r|$", r"$r_I^{\max}$", r"$r_K^{\max}$"]
    lines = ["% Generated automatically. End-size experiment; no seed averaging.",
             "% Restriction columns are maxima across actual pass starts; see restrictions.csv for every pass.",
             r"\begin{table}[tbp]", r"\centering", r"\begingroup\scriptsize",
             r"\setlength{\tabcolsep}{1.5pt}", r"\renewcommand{\arraystretch}{1.05}",
             r"\caption{" + caption + "}", r"\begin{tabular}{@{}lrrrrrrrrr@{}}", r"\hline",
             " & ".join(headings) + ROW_END, r"\hline"]
    for row in rows:
        values = [row["pointset"], row["n_points"] if row["n_points"] != "" else "--", f"{row['epsilon_percent']:g}"]
        values += [str(row[key]) if isinstance(row[key], int) else compact_cell(row[key]) for _, key in COLUMNS]
        lines.append(" & ".join(latex(value) for value in values) + ROW_END)
    lines += [r"\hline", r"\end{tabular}", r"\endgroup", r"\end{table}"]
    (directory / "results.tex").write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--input", type=Path, help="Override the shared-manifest selection with one SCC file")
    selection.add_argument("--640-only", dest="only_640", action="store_true",
                           help="Measure only the three predetermined 640-point sample-1 inputs")
    selection.add_argument("--immuno-only", type=int, choices=POINT_COUNTS,
                           help="Measure only the immune H1 sample at this target size")
    selection.add_argument("--immuno-640-only", dest="immuno_only", action="store_const", const=640,
                           help="Alias for --immuno-only 640")
    parser.add_argument("--manifest", type=Path, default=DATA / "manifest.csv")
    parser.add_argument("--families", nargs="+", choices=FAMILIES,
                        help="Keep only these families from the manifest selection")
    parser.add_argument("--generators-only", action="store_true",
                        help="Stop after left/right E-generator selection; write separate dimension and generator tables")
    parser.add_argument("--time-steps", action="store_true",
                        help="Include reduction seconds for every I/K restriction pass in restrictions.csv (default: off)")
    parser.add_argument("--percentages", type=float, nargs="+", default=EPSILON_PERCENTAGES,
                        help="Percentages of each generator/relation bounding-box side (default: 0.5 1 2)")
    parser.add_argument("--worker", type=Path, default=ROOT / "build/local-release/hom_reduction_benchmark")
    parser.add_argument("--output-dir", type=Path,
                        help="Results directory (default: results/End_sizes, End_sizes_640, or End_sizes_immuno_N for a focused run)")
    parser.add_argument("--timeout", type=float, default=TIMEOUT_SECONDS,
                        help="Seconds per input and epsilon, except LargeHypoxicRegion H1 (default: 600)")
    parser.add_argument("--large-timeout", type=float, default=LARGE_TIMEOUT_SECONDS,
                        help="Seconds per epsilon for LargeHypoxicRegion H1 (default: 3600)")
    parser.add_argument("--memory-limit-mb", type=float, default=MEMORY_MB)
    parser.add_argument("--nilpotent-reduction", action=argparse.BooleanOptionalAction, default=True,
                        help="Measure nilpotent reduction compression (default: enabled)")
    parser.add_argument("--complement-orbit-generators", action=argparse.BooleanOptionalAction, default=True,
                        help="Measure complement-orbit compression (default: enabled)")
    add_generator_arguments(parser)
    args = parser.parse_args()
    if args.input is not None and args.families:
        parser.error("--families filters manifest inputs and cannot be combined with --input")
    search_settings = generator_settings(args)
    enabled = {"nilpotent_reduction": args.nilpotent_reduction,
               "complement_orbit": args.complement_orbit_generators,
               "combined": args.nilpotent_reduction and args.complement_orbit_generators}
    method_flags = [f"--{'' if enabled[method] else 'no-'}{flag}"
                    for method, flag in (("nilpotent_reduction", "nilpotent-reduction"),
                                         ("complement_orbit", "complement-orbit-generators"))] + generator_flags(search_settings)
    if (any(not math.isfinite(p) or p < 0 for p in args.percentages)
            or not math.isfinite(args.timeout) or args.timeout <= 0
            or not math.isfinite(args.large_timeout) or args.large_timeout <= 0
            or not math.isfinite(args.memory_limit_mb) or args.memory_limit_mb < 0):
        parser.error("Percentages/memory must be finite and nonnegative; timeout must be finite and positive")
    args.worker, args.manifest = args.worker.resolve(), args.manifest.resolve()
    if not args.worker.is_file():
        parser.error("Build the worker first: cmake --build build/local-release --target hom_reduction_benchmark")
    help_text = subprocess.run([str(args.worker), "--help"], check=True, capture_output=True, text=True, timeout=10).stdout
    if "hom_algorithm=optimised" not in help_text or "end_size_pipeline=1" not in help_text:
        parser.error("Rebuild hom_reduction_benchmark: the experiment requires A-1/2 and the I/K restriction pipeline")
    if args.generators_only and "--generators-only" not in help_text:
        parser.error("Rebuild hom_reduction_benchmark: --generators-only is required")
    if args.time_steps and "--time-steps" not in help_text:
        parser.error("Rebuild hom_reduction_benchmark: --time-steps is required")
    if args.input is None and not args.manifest.is_file():
        parser.error(f"Prepare the shared inputs first: missing {args.manifest}")
    percentages = list(dict.fromkeys(args.percentages))
    lock = input_lock(args.manifest.parent) if args.input is None else nullcontext()
    with lock:
        if args.input is None:
            try:
                inputs = select_inputs(args.manifest, immuno_points=args.immuno_only, only_640=args.only_640)
            except ValueError as error:
                parser.error(str(error))
            if args.families:
                inputs = [sample for sample in inputs if sample["family"] in args.families]
        else:
            path = args.input.resolve()
            name = re.sub(r"_H\d+$", "", path.stem)
            homology = re.search(r"_H(\d+)$", path.stem)
            inputs = [dict(dataset=name, family=name, n_points="", seed="", input=str(path),
                           homology=homology[1] if homology else "")]
        if not inputs:
            parser.error("No End-size inputs match the requested selection")
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
        directory = (args.output_dir or HERE / "results" / (
            "End_sizes_640" if args.only_640 else
            f"End_sizes_immuno_{args.immuno_only}" if args.immuno_only is not None else "End_sizes")).resolve()
        logs = directory / run
        logs.mkdir(parents=True)
        (logs / "settings.json").write_text(json.dumps(dict(
            inputs=inputs, manifest=str(args.manifest) if args.input is None else None,
            manifest_sha256=digest(args.manifest) if args.input is None else None,
            worker=str(args.worker), worker_sha256=worker_hash, hom_algorithm="optimised",
            percentages=percentages, points=[640] if args.only_640 else [args.immuno_only] if args.immuno_only is not None else POINT_COUNTS,
            families=list(dict.fromkeys(sample["family"] for sample in inputs)), synthetic_seed=1,
            adaptive=not args.generators_only, generators_only=args.generators_only, time_steps=args.time_steps,
            methods=enabled, generator_settings=search_settings,
            timeout_seconds=args.timeout, large_timeout_seconds=args.large_timeout,
            memory_limit_mb=args.memory_limit_mb), indent=2) + "\n")
        report_files = (*REPORT_FILES, "dimensions.tex", "generators.tex") if args.generators_only else REPORT_FILES
        for name in report_files:
            if (directory / name).exists():
                shutil.copyfile(directory / name, logs / ("previous_" + name))
        observations = []

        def report():
            rows = report_rows(observations, inputs, percentages)
            write_reports(logs, rows, generators_only=args.generators_only, time_steps=args.time_steps)
            for name in report_files:
                shutil.copyfile(logs / name, directory / name)

        report()
        failed = False
        print(f"End sizes: {len(inputs)} inputs x {len(percentages)} epsilons; results: {logs}", flush=True)
        try:
            with (logs / "observations.jsonl").open("a") as records:
                for sample in inputs:
                    for percentage in percentages:
                        path = Path(sample["input"])
                        if digest(path) != sample["input_sha256"]:
                            raise SystemExit(f"Input changed during the experiment: {path}; saved results are in {logs}")
                        prefix = logs / f"{path.stem}_epsilon_{percentage:g}pct"
                        stdout_path, stderr_path = Path(str(prefix) + ".stdout.txt"), Path(str(prefix) + ".stderr.txt")
                        command = [str(args.worker), str(path), str(percentage), "--percent",
                                   "--memory-limit-mb", str(args.memory_limit_mb), *method_flags]
                        if args.generators_only:
                            command.append("--generators-only")
                        if args.time_steps:
                            command.append("--time-steps")
                        epsilon = tuple(side * (percentage / 100) for side in sample["bounding_box_sides"])
                        timeout = args.large_timeout if (sample["dataset"] == "LargeHypoxicRegion"
                                                         and sample.get("homology") == "1") else args.timeout
                        row = dict(sample, run=run, worker_sha256=worker_hash, epsilon_percent=percentage,
                                   epsilon=epsilon, shift=tuple(2 * value for value in epsilon),
                                   timeout_seconds=timeout,
                                   methods=enabled, generator_settings=search_settings,
                                   stdout_log=str(stdout_path.relative_to(directory)), stderr_log=str(stderr_path.relative_to(directory)))
                        print(f"RUN {path.stem}: epsilon={percentage:g}% ({row['epsilon']})", flush=True)
                        status, seconds, error = measure(command, stdout_path, stderr_path, timeout)
                        try:
                            row.update(parse_output(stdout_path.read_text()))
                        except (ValueError, KeyError) as parse_error:
                            status, error = "error", f"Cannot parse worker output: {parse_error}"
                        required = [key for _, key in (
                            (*COLUMNS[:5], *TIMINGS[:5]) if args.generators_only else (*COLUMNS, *TIMINGS)
                        ) if key != "process_seconds"]
                        if status == "ok" and (row.get("hom_algorithm") != "optimised" or any(key not in row for key in required)):
                            status, error = "error", "Worker output is incomplete or uses the wrong Hom algorithm"
                        if status == "ok" and not args.generators_only:
                            for phase in ("I", "K"):
                                trace = [item for item in row["restrictions"] if item["phase"] == phase]
                                if (len(trace) != row[f"{phase}_restriction_passes"] or
                                    max((item["retained"] for item in trace), default=0) != row[f"{phase}_restriction_max"] or
                                    sum(item["retained"] for item in trace) != row[f"{phase}_restriction_total"] or
                                    any(item["retained"] > item["candidates"] for item in trace) or
                                    (args.time_steps and any(not math.isfinite(item.get("seconds", math.nan))
                                                            or item["seconds"] < 0 for item in trace))):
                                    status, error = "error", f"Inconsistent {phase} per-pass restriction counts"
                        row.update(status=status, process_seconds=seconds, error=error, command=command)
                        Path(str(prefix) + ".json").write_text(json.dumps(row, indent=2) + "\n")
                        records.write(json.dumps(row, separators=(",", ":")) + "\n")
                        records.flush()
                        observations.append(row)
                        report()
                        failed |= status != "ok"
                        print(f"  {status}: {seconds:.2f}s; full Hom={row.get('End_d_basis_I_size', '?')}, "
                              f"quotient={row.get('quotient_basis_I_size', '?')}", flush=True)
        except KeyboardInterrupt:
            print(f"Stopped. Completed measurements and partial tables are saved in {logs}", flush=True)
            return 130
        print(f"Results: {directory / 'results.csv'}; {directory / 'results.tex'}", flush=True)
        return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
