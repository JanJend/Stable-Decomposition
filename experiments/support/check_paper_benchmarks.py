"""Small end-to-end check: python3 -m experiments.support.check_paper_benchmarks."""
import csv
import json
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import tempfile
from unittest.mock import patch

from . import HERE, ROOT, input_lock
from .paper_benchmarks import column_storage, execute, main as benchmark
from .paper_tables import VARIANTS, SPEEDUP_NAMES, aggregate, cell, speedup_cell, write_latex, write_reports
from .runtime import digest
from .scaling import SWITCHES
from .generator_options import COUNT_DEFAULTS, add_generator_arguments, generator_sweep
import argparse


def check_configurations(directory):
    """Assert the switches actually passed, and keep epsilon measurements separate."""
    fake = directory / "fake-pruning"
    fake.write_text(f"#!{sys.executable}\n" + "import sys\n" +
                    "print('Column storage: compact')\n" +
                    f"switches = {SWITCHES!r}\n" +
                    "algorithm = next((arg.split('=', 1)[1] for arg in sys.argv if arg.startswith('--hom-algorithm=')), 'optimised')\n" +
                    "print('Module pruning options (0=off, 1=on): ' + ' '.join(\n"
                    "    name.replace('-', '_') + '=' + str(int('--full-end-space' in sys.argv) if name == 'full-end-space' else int('--no-' + name not in sys.argv))\n"
                    "    for name in switches) + ' hom_algorithm=' + algorithm)\n"
                    "print('Composition generator search: ' + ' '.join(arg[2:] for arg in sys.argv if arg.startswith('--composition-')))\n"
                    "print('Initial Hom quotient: 8; I maps: 3; K maps: 2')\n"
                    "print('I iterations: 1; K iterations: 1')\n"
                    "epsilon = float(sys.argv[sys.argv.index('--epsilon') + 1]) if '--epsilon' in sys.argv else 0\n"
                    "print(f'Total 1 {epsilon:.6f}')\n")
    fake.chmod(0o755)
    assert column_storage(fake) == "compact"
    assert VARIANTS == ["baseline", *SWITCHES]
    observations = []
    for percentage in (1, 8):
        for variant in [*VARIANTS, "all-optimizations"]:
            result = execute(fake, directory / "unused.scc", percentage / 100, variant, 10, 0, "compact")
            assert result["status"] == "ok", result
            flags = {argument[2:] for argument in result["command"] if argument.startswith("--")}
            optimizations = set(SWITCHES) - {"full-end-space"}
            enabled = {name for name in optimizations if 'no-' + name not in flags}
            if "full-end-space" in flags:
                enabled.add("full-end-space")
            expected = optimizations if variant == "all-optimizations" else ({variant} if variant in SWITCHES else set())
            assert enabled == expected, (variant, enabled)
            assert ("--no-single-solve-quotient" not in result["command"]) == (
                variant in ("single-solve-quotient", "all-optimizations"))
            assert ("--update-i-target" in result["command"]) == ("update-i-target" in expected)
            assert not flags & {"fuse-i-preimage", "defer-k-image-minimization", "throttle-progress", "dense-kernel-pivots", "no-dense-kernel-pivots"}
            assert "hom-algorithm=optimised" in flags
            observations.append(dict(dataset="fake", family="fake", seed=42, n_points=20,
                                     target_points=20, homology=1, generators=1, relations=1,
                                     epsilon_percent=percentage, epsilon=percentage / 100,
                                     variant=variant, repetition=1, timeout_seconds=10, **result))
    vector_result = execute(fake, directory / "unused.scc", (0.01, 0.001), "baseline", 10, 0, "compact")
    command = vector_result["command"]
    position = command.index("--epsilon-vector")
    assert command[position + 1:position + 3] == ["0.01", "0.001"]
    assert "--epsilon" not in command
    parser = argparse.ArgumentParser()
    add_generator_arguments(parser, rounds=[3, 6, 12], modes=["hybrid"])
    configs = generator_sweep(parser.parse_args(["--composition-orbit-seed-mode", "sum", "random"]))
    assert len(configs) == 6
    for name, settings in configs.items():
        result = execute(fake, directory / "unused.scc", 0.01, name, 10, 0, "compact", settings)
        assert result["status"] == "ok", result
        assert result["generator_settings"]["composition-orbit-rounds"] == str(settings["composition-orbit-rounds"])
        assert result["i_maps"] == 3 and result["k_maps"] == 2
        assert "--no-nilpotent-reduction" not in result["command"]
        assert "--no-complement-orbit-generators" not in result["command"]
    assert execute(fake, directory / "unused.scc", 0.01, "baseline", 10, 0, "vector")["status"] == "error"
    write_reports(observations, directory, "optimizations")
    with (directory / "table.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 2 and [float(row["epsilon_percent"]) for row in rows] == [1, 8]
    assert [float(row["baseline_seconds"]) for row in rows] == [0.01, 0.08]
    assert {name[:-8] for name in rows[0] if name.endswith("_seconds")} == set(VARIANTS)
    assert "optional switches off" in (directory / "table.tex").read_text()
    assert "Fast I" in (directory / "table.tex").read_text()
    assert "update-i-target_seconds" in rows[0]
    with (directory / "speedups.csv").open() as stream:
        percentages = list(csv.DictReader(stream))
    assert len(percentages) == 2
    assert {name.removesuffix("_speedup_percent") for name in percentages[0]
            if name.endswith("_speedup_percent")} == set(SPEEDUP_NAMES)
    assert all(row[name + "_speedup_percent"] == "+0.0%" for row in percentages for name in SPEEDUP_NAMES)
    # The Quot column measures the optimized baseline against the full-End
    # control, so it must reverse the ratio used by the other comparisons.
    altered = [dict(row, pruning_seconds=row["pruning_seconds"] * 4)
               if row["variant"] == "full-end-space" else row for row in observations]
    write_reports(altered, directory, "optimizations")
    with (directory / "speedups.csv").open() as stream:
        assert all(row["hom-quotient_speedup_percent"] == "+300.0%" for row in csv.DictReader(stream))
    speedup_latex = (directory / "speedups.tex").read_text()
    assert " & Base & " not in speedup_latex
    assert speedup_latex.index(" & Nilpotent & ") < speedup_latex.index(" & Adaptive & ")
    assert "Orbit" in speedup_latex
    assert "Direct Q" in speedup_latex
    assert "End0 is still computed separately when needed for module actions" in speedup_latex
    altered = [dict(row, pruning_seconds=row["pruning_seconds"] / 4)
               if row["variant"] == "single-solve-quotient" else row for row in observations]
    write_reports(altered, directory, "optimizations")
    with (directory / "speedups.csv").open() as stream:
        assert all(row["single-solve-quotient_speedup_percent"] == "+300.0%" for row in csv.DictReader(stream))
    baseline = dict(completed=1, status="ok", pruning_seconds=2)
    assert speedup_cell(baseline, dict(baseline, pruning_seconds=1)) == "+100.0%"
    assert speedup_cell(baseline, dict(baseline, pruning_seconds=4)) == "-50.0%"
    assert speedup_cell(baseline, dict(baseline, status="partial")) == "+0.0% [partial]"
    assert speedup_cell(baseline, dict(baseline, completed=0, status="timeout")) == "n/a"
    assert speedup_cell(dict(baseline, completed=0, status="memory_limit"), baseline) == "n/a"
    assert speedup_cell(baseline, dict(baseline, pruning_seconds=0)) == "n/a"
    assert speedup_cell(baseline, None) == "pending"
    mixed = aggregate([observations[0], dict(observations[0], status="timeout")])
    assert mixed["status"] == "partial" and "[1/2]" in cell(mixed)


def check_latex_layout(directory):
    """All 27 averaged rows fit a standard portrait table without page changes."""
    entry = dict(column_storage="compact")
    entry.update((name + "_speedup_percent", "+850.3% [partial]") for name in SPEEDUP_NAMES)
    entries = [dict(entry, family=family, n_points=n, epsilon_percent=epsilon)
               for family in ("gauss_clust", "immuno", "torus_gauss")
               for n in (160, 320, 640) for epsilon in (0.5, 1, 2)]
    write_latex(directory / "wide.tex", entries, "optimizations", percentages=True)
    source = (directory / "wide.tex").read_text()
    assert source.count(r"\begin{table}") == source.count(r"\begin{tabular}") == source.count(r"\caption{") == 1
    assert not any(token in source for token in (r"\begin{landscape}", r"\resizebox", r"\geometry", "Seed", "$H$"))
    assert "Point set" in source and r"$\epsilon$ (\%)" in source
    assert r"gauss\_clust" in source and "immuno" in source
    assert sum(" & " in line for line in source.splitlines()) == len(entries) + 1
    if not shutil.which("pdflatex"):
        return
    document = directory / "layout.tex"
    document.write_text("\\documentclass{article}\n"
                        "\\begin{document}\\input{wide.tex}\\end{document}\n")
    for _ in range(2):
        result = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", document.name],
                                cwd=directory, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout[-3000:]
    assert "Overfull" not in result.stdout, result.stdout


def check_entry_points():
    # Execute each public script as __main__, without launching the full experiment.
    sys.path.insert(0, str(HERE))
    try:
        from support import paper_benchmarks
        for kind in ("runtime", "optimizations"):
            with patch.object(paper_benchmarks, "main") as run, patch("subprocess.run") as generate:
                with patch.object(sys, "argv", [f"benchmark_{kind}.py"]):
                    with patch.object(Path, "is_file", return_value=False):
                        settings = runpy.run_path(str(HERE / f"benchmark_{kind}.py"), run_name="__main__")
                run.assert_called_once()
                assert run.call_args.args == (kind,)
                assert run.call_args.kwargs["points"] == settings["POINT_COUNTS"]
                assert run.call_args.kwargs["epsilon_percents"] == settings["EPSILON_PERCENTAGES"]
                if kind == "optimizations":
                    assert set(run.call_args.kwargs["families"]) == {"torus_gauss", "immuno", "gauss_clust"}
                    configs = run.call_args.kwargs["generator_configs"]
                    assert list(configs) == [f"generators-{mode}-{rounds}"
                                             for mode in settings["ORBIT_SEED_MODES"]
                                             for rounds in settings["ORBIT_ROUNDS"]]
                else:
                    assert "families" not in run.call_args.kwargs
                assert "manifest" not in run.call_args.kwargs, "Both experiments use the shared default manifest"
                generate.assert_not_called()
    finally:
        sys.path.pop(0)


def check_seed_averages(directory):
    base = dict(family="gauss_clust", n_points=80, target_points=80, homology=0,
                generators=1, relations=1, epsilon_percent=1, epsilon=0.01,
                column_storage="compact", timeout_seconds=10, status="ok")
    observations = []
    # Unequal repeat counts must not give one seed extra weight. Seed medians
    # are 1, 3, 5; their mean is 3, distinct from the pooled timing median.
    for seed, times in enumerate(([1, 1, 100], [2, 4], [5]), 1):
        for repetition, seconds in enumerate(times, 1):
            observations.append(dict(base, dataset=f"cluster_80_{seed}", seed=seed, variant="baseline",
                                     repetition=repetition, pruning_seconds=seconds, process_seconds=seconds))
        for variant, seconds in (("projected-preimages", (1, 1, 4)[seed-1]), ("full-end-space", (2, 6, 10)[seed-1])):
            observations.append(dict(base, dataset=f"cluster_80_{seed}", seed=seed, variant=variant,
                                     repetition=1, pruning_seconds=seconds, process_seconds=seconds))
    counts = {("gauss_clust", 80, 0): 3}
    write_reports(observations, directory, "optimizations", seed_counts=counts)
    with (directory / "table.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 1 and "seed" not in rows[0] and "homology" not in rows[0]
    assert float(rows[0]["baseline_seconds"]) == 3
    assert float(rows[0]["projected-preimages_seconds"]) == 2
    with (directory / "speedups.csv").open() as stream:
        percentages = list(csv.DictReader(stream))
    assert percentages[0]["projected-preimages_speedup_percent"] == "+50.0%", percentages
    assert percentages[0]["hom-quotient_speedup_percent"] == "+100.0%", percentages
    with (directory / "summary.csv").open() as stream:
        assert len(list(csv.DictReader(stream))) == 9, "Keep individual seed summaries"
    write_reports([row for row in observations if row["seed"] != 3], directory, "optimizations", seed_counts=counts)
    with (directory / "table.csv").open() as stream:
        assert next(csv.DictReader(stream))["baseline_seconds"] == "pending"
    failed = [dict(row, status="timeout") if row["seed"] == 3 and row["variant"] == "baseline" else row
              for row in observations]
    write_reports(failed, directory, "optimizations", seed_counts=counts)
    with (directory / "table.csv").open() as stream:
        assert next(csv.DictReader(stream))["baseline_seconds"].startswith("timeout")
    with (directory / "speedups.csv").open() as stream:
        assert next(csv.DictReader(stream))["projected-preimages_speedup_percent"] == "n/a"


def check_optimization_selection(directory):
    settings = dict(points=[80, 160, 320], families=["torus_gauss", "immuno", "gauss_clust"],
                    epsilon_percents=[1, 8], repeats=1, memory_mb=0,
                    pruning=directory / "fake-pruning", output_dir=directory / "selection")
    result = dict(status="ok", pruning_seconds=0.001, process_seconds=0.002, column_storage="compact")
    with patch("experiments.support.paper_benchmarks.execute", return_value=result):
        output = benchmark("optimizations", **settings)
    rows = [json.loads(line) for line in (output / "observations.jsonl").read_text().splitlines()]
    assert len(rows) == 21 * 2 * len(VARIANTS)
    assert {row["family"] for row in rows} == set(settings["families"])
    assert {row["n_points"] for row in rows} == set(settings["points"])
    assert {row["epsilon_percent"] for row in rows} == {1, 8}
    with (output / "table.csv").open() as stream:
        assert len(list(csv.DictReader(stream))) == 18
    assert json.loads((output / "run.json").read_text())["configuration"]["families"] == sorted(settings["families"])


def check_input_safety(directory):
    bundle = directory / "inputs"
    bundle.mkdir()
    child = ("from pathlib import Path\nimport sys\n"
             "from experiments.support import input_lock\n"
             "with input_lock(Path(sys.argv[1]), writing=sys.argv[2] == 'write'):\n    pass\n")
    for writing in (False, True):
        with input_lock(bundle, writing=writing):
            blocked = subprocess.run([sys.executable, "-c", child, str(bundle), "read" if writing else "write"],
                                     cwd=ROOT, capture_output=True, text=True)
            assert blocked.returncode != 0 and "Wait for it to finish" in blocked.stderr, blocked
            if not writing:
                with input_lock(bundle):
                    pass  # Independent readers can share the inputs.
    released = subprocess.run([sys.executable, "-c", child, str(bundle), "write"], cwd=ROOT)
    assert released.returncode == 0

    presentation = bundle / "valid.scc"
    presentation.write_text("scc2020\n2\n0 0 0\n")
    first = dict(dataset="first", family="fake", seed=1, n_points=20, target_points=20,
                 homology=0, generation_status="ok", scc=presentation.name, scc_sha256=digest(presentation))
    manifest = bundle / "manifest.csv"
    with manifest.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(first))
        writer.writeheader()
        writer.writerows([first, dict(first, dataset="last", seed=2, scc_sha256="stale")])
    output = directory / "invalid-input-results"
    with patch("experiments.support.paper_benchmarks.execute") as timing:
        try:
            benchmark("optimizations", points=[20], epsilon_percents=[1], manifest=manifest,
                      pruning=directory / "fake-pruning", output_dir=output)
        except SystemExit as error:
            assert "different from the manifest" in str(error) and str(presentation) in str(error)
        else:
            raise AssertionError("Changed input was not rejected before timing")
        timing.assert_not_called()
        assert not output.exists(), "Invalid inputs must not create a partial experiment"


def main():
    check_entry_points()
    pruning = ROOT / "build/pruning"
    # Missing builds must produce an instruction, never trigger compilation.
    with patch.object(Path, "is_file", return_value=False), patch("subprocess.run") as processes:
        for kind in ("runtime", "optimizations"):
            try:
                benchmark(kind, points=[80], epsilon_percents=[1])
            except SystemExit as error:
                assert str(pruning) in str(error) and "Build it with CMake" in str(error)
            else:
                raise AssertionError("Missing executable was not rejected")
        processes.assert_not_called()
    with tempfile.TemporaryDirectory(prefix="paper-benchmark-check-") as temporary:
        temporary = Path(temporary).resolve()
        check_configurations(temporary)
        check_seed_averages(temporary)
        check_optimization_selection(temporary)
        check_input_safety(temporary)
        check_latex_layout(temporary)
        with (HERE / "data/manifest.csv").open() as stream:
            samples = [row for row in csv.DictReader(stream) if row["target_points"] == "20" and
                       (row["family"] == "immuno" or (row["family"] == "torus_gauss" and row["seed"] == "1" and row["homology"] == "1"))]
        assert len(samples) == 2
        for row in samples:
            row["scc"] = str(HERE / "data" / row["scc"])
        manifest = temporary / "manifest.csv"
        with manifest.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(samples[0]))
            writer.writeheader()
            writer.writerows(samples)
        for kind in ("runtime", "optimizations"):
            output = temporary / kind
            settings = dict(points=[20], epsilon_percents=[1, 8], repeats=2, repeat_below=60,
                            memory_mb=0, output_dir=output, manifest=manifest)
            if kind == "optimizations":
                settings["generator_configs"] = {
                    f"generators-{mode}-{rounds}": dict(COUNT_DEFAULTS, **{
                        "composition-orbit-rounds": rounds, "composition-orbit-seed-mode": mode})
                    for mode, rounds in (("sum", 3), ("random", 6), ("hybrid", 12))}
            with patch("subprocess.run", wraps=subprocess.run) as processes:
                benchmark(kind, **settings)
                commands = [call.args[0] for call in processes.call_args_list]
                assert [str(pruning), "--help"] in commands
                assert all(Path(command[0]).name != "cmake" for command in commands)
            raw = output / "observations.jsonl"
            before = raw.read_bytes()
            observations = [json.loads(line) for line in before.splitlines()]
            assert len(observations) == (8 if kind == "runtime" else 8 * (len(VARIANTS) + 3))
            assert all(row["status"] == "ok" and row["pruning_seconds"] > 0 for row in observations), observations
            assert all(row["process_seconds"] >= row["pruning_seconds"] for row in observations)
            for row in observations:
                assert len(row["epsilon"]) == 2
                position = row["command"].index("--epsilon-vector")
                assert list(map(float, row["command"][position + 1:position + 3])) == row["epsilon"]
            config = json.loads((output / "run.json").read_text())["configuration"]
            assert all(row["pruning_sha256"] == config["pruning_sha256"] for row in observations)
            assert {row["family"] for row in observations} == {"immuno", "torus_gauss"}
            assert {row["epsilon_percent"] for row in observations} == {1, 8}
            if kind == "optimizations":
                with (output / "generators.csv").open() as stream:
                    generator_rows = list(csv.DictReader(stream))
                assert len(generator_rows) == 12
                assert {int(row["composition-orbit-rounds"]) for row in generator_rows} == {3, 6, 12}
                assert all(row["i_maps"] != "" and row["k_maps"] != "" for row in generator_rows)
                assert config["generator_configs"] == settings["generator_configs"]
            # A partial run must finish the missing epsilon/repetition, without duplicates.
            raw.write_bytes(b"\n".join(before.splitlines()[:-1]) + b"\n")
            benchmark(kind, **settings)
            after = raw.read_bytes()
            restored = [json.loads(line) for line in after.splitlines()]
            assert len(restored) == len(observations)
            assert restored[:-1] == observations[:-1]
            for key in ("dataset", "epsilon_percent", "variant", "repetition"):
                assert restored[-1][key] == observations[-1][key]
            table = (output / "table.csv").read_bytes()
            speedups = (output / "speedups.csv").read_bytes() if kind == "optimizations" else None
            (output / "table.csv").unlink()
            if speedups is not None:
                (output / "speedups.csv").unlink()
            assert benchmark(kind, **settings) == output
            assert raw.read_bytes() == after and (output / "table.csv").read_bytes() == table
            if speedups is not None:
                assert (output / "speedups.csv").read_bytes() == speedups
                assert (output / "speedups.tex").is_file()
            else:
                assert not (output / "speedups.csv").exists()
            with (output / "table.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            assert len(rows) == 4
            # A changed protocol gets fresh results; subsequent runs resume that directory.
            changed_settings = dict(settings, epsilon_percents=[2], repeats=1)
            changed = benchmark(kind, **changed_settings)
            changed_raw = (changed / "observations.jsonl").read_bytes()
            assert changed != output and raw.read_bytes() == after
            assert benchmark(kind, **changed_settings) == changed
            assert (changed / "observations.jsonl").read_bytes() == changed_raw
        small = HERE / "data/presentations/torus_gauss/torus_gauss_160_1_H1.scc"
        timeout = execute(pruning, small, 0.005, "baseline", 0.01, 0)
        memory = execute(pruning, small, 0.005, "baseline", 10, 1)
        assert timeout["status"] == "timeout" and "pruning_seconds" not in timeout
        assert memory["status"] == "memory_limit" and "pruning_seconds" not in memory
        if shutil.which("pdflatex"):
            document = temporary / "tables.tex"
            document.write_text("\\documentclass{article}\n\\usepackage{longtable,pdflscape}\n"
                                "\\begin{document}\n\\input{runtime/table.tex}\n"
                                "\\input{optimizations/table.tex}\n"
                                "\\input{optimizations/speedups.tex}\n\\input{optimizations/generators.tex}\n\\end{document}\n")
            result = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", document.name],
                                    cwd=temporary, capture_output=True, text=True)
            assert result.returncode == 0, result.stdout[-3000:]
    print("PASS: input locks and preflight, entry points, individual switches, real/synthetic inputs, epsilon sweep, repeats, resume, tables and limits")


if __name__ == "__main__":
    main()
