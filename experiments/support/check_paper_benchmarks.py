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

from . import HERE, ROOT
from .paper_benchmarks import column_storage, execute, main as benchmark
from .paper_tables import VARIANTS, aggregate, cell, write_reports
from .scaling import SWITCHES


def check_configurations(directory):
    """Assert the switches actually passed, and keep epsilon measurements separate."""
    fake = directory / "fake-pruning"
    fake.write_text(f"#!{sys.executable}\n" + "import sys\n" +
                    "print('Column storage: compact')\n" +
                    f"switches = {SWITCHES!r}\n" +
                    "print('Module pruning options (0=off, 1=on): ' + ' '.join(\n"
                    "    name.replace('-', '_') + '=' + str(int('--all-optimizations' in sys.argv or\n"
                    "    '--' + name in sys.argv)) for name in switches))\n"
                    "print('I iterations: 1; K iterations: 1')\n"
                    "epsilon = float(sys.argv[sys.argv.index('--epsilon') + 1]) if '--epsilon' in sys.argv else 0\n"
                    "print(f'Total 1 {epsilon:.6f}')\n")
    fake.chmod(0o755)
    assert column_storage(fake) == "compact"
    assert VARIANTS == ["baseline", *SWITCHES]
    observations = []
    for percentage in (1, 10):
        for variant in [*VARIANTS, "all-optimizations"]:
            result = execute(fake, directory / "unused.scc", percentage / 100, variant, 10, 0, "compact")
            assert result["status"] == "ok", result
            flags = {argument[2:] for argument in result["command"] if argument.startswith("--")}
            enabled = set(SWITCHES) if "all-optimizations" in flags else flags & set(SWITCHES)
            expected = set(SWITCHES) if variant == "all-optimizations" else ({variant} if variant != "baseline" else set())
            assert enabled == expected, (variant, enabled)
            observations.append(dict(dataset="fake", family="fake", seed=42, n_points=20,
                                     target_points=20, homology=1, generators=1, relations=1,
                                     epsilon_percent=percentage, epsilon=percentage / 100,
                                     variant=variant, repetition=1, timeout_seconds=10, **result))
    assert execute(fake, directory / "unused.scc", 0.01, "baseline", 10, 0, "vector")["status"] == "error"
    write_reports(observations, directory, "optimizations")
    with (directory / "table.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 2 and [float(row["epsilon_percent"]) for row in rows] == [1, 10]
    assert [float(row["baseline_seconds"]) for row in rows] == [0.01, 0.1]
    assert {name[:-8] for name in rows[0] if name.endswith("_seconds")} == set(VARIANTS)
    assert "all switches off" in (directory / "table.tex").read_text()
    mixed = aggregate([observations[0], dict(observations[0], status="timeout")])
    assert mixed["status"] == "partial" and "[1/2]" in cell(mixed)


def check_entry_points():
    # Execute each public script as __main__, without launching the full experiment.
    sys.path.insert(0, str(HERE))
    try:
        from support import paper_benchmarks
        for kind in ("runtime", "optimizations"):
            with patch.object(paper_benchmarks, "main") as run, patch("subprocess.run") as generate:
                with patch.object(sys, "argv", [f"benchmark_{kind}.py"]):
                    with patch.object(Path, "is_file", return_value=False):
                        runpy.run_path(str(HERE / f"benchmark_{kind}.py"), run_name="__main__")
                run.assert_called_once()
                assert run.call_args.args == (kind,)
                assert run.call_args.kwargs["points"] == (
                    [80, 160, 320, 640] if kind == "optimizations" else [20, 32, 64, 128, 256, 512, 1024, 2048])
                assert run.call_args.kwargs["epsilon_percents"] == [1, 2, 5, 10]
                if kind == "optimizations":
                    assert run.call_args.kwargs["manifest"] == HERE / "optimization_data/manifest.csv"
                    generate.assert_called_once_with(
                        [sys.executable, str(HERE / "generate_data.py"), "--sizes", "80", "160", "320", "640",
                         "--output-dir", str(HERE / "optimization_data")], check=True)
                    with patch.object(sys, "argv", [f"benchmark_{kind}.py"]), patch.object(Path, "is_file", return_value=True):
                        runpy.run_path(str(HERE / f"benchmark_{kind}.py"), run_name="__main__")
                    assert generate.call_count == 1, "Prepared inputs should be reused"
                else:
                    generate.assert_not_called()
    finally:
        sys.path.pop(0)


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
        with (HERE / "manifest.csv").open() as stream:
            samples = [row for row in csv.DictReader(stream) if row["target_points"] == "20" and
                       (row["family"] == "immuno" or (row["family"] == "torus_gaussian" and row["seed"] == "42"))]
        assert len(samples) == 2
        for row in samples:
            row["scc"] = str(HERE / row["scc"])
        manifest = temporary / "manifest.csv"
        with manifest.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(samples[0]))
            writer.writeheader()
            writer.writerows(samples)
        for kind in ("runtime", "optimizations"):
            output = temporary / kind
            settings = dict(points=[20], epsilon_percents=[1, 10], repeats=2, repeat_below=60,
                            memory_mb=0, output_dir=output, manifest=manifest)
            with patch("subprocess.run", wraps=subprocess.run) as processes:
                benchmark(kind, **settings)
                commands = [call.args[0] for call in processes.call_args_list]
                assert [str(pruning), "--help"] in commands
                assert all(Path(command[0]).name != "cmake" for command in commands)
            raw = output / "observations.jsonl"
            before = raw.read_bytes()
            observations = [json.loads(line) for line in before.splitlines()]
            assert len(observations) == (8 if kind == "runtime" else 88)
            assert all(row["status"] == "ok" and row["pruning_seconds"] > 0 for row in observations), observations
            assert all(row["process_seconds"] >= row["pruning_seconds"] for row in observations)
            config = json.loads((output / "run.json").read_text())["configuration"]
            assert all(row["pruning_sha256"] == config["pruning_sha256"] for row in observations)
            assert {row["family"] for row in observations} == {"immuno", "torus_gaussian"}
            assert {row["epsilon_percent"] for row in observations} == {1, 10}
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
            (output / "table.csv").unlink()
            assert benchmark(kind, **settings) == output
            assert raw.read_bytes() == after and (output / "table.csv").read_bytes() == table
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
        small = HERE / "test_presentations/torus_gaussian/torus_gaussian_n128_seed42_h1.scc"
        timeout = execute(pruning, small, 0.005, "baseline", 0.01, 0)
        memory = execute(pruning, small, 0.005, "baseline", 10, 1)
        assert timeout["status"] == "timeout" and "pruning_seconds" not in timeout
        assert memory["status"] == "memory_limit" and "pruning_seconds" not in memory
        if shutil.which("pdflatex"):
            document = temporary / "tables.tex"
            document.write_text("\\documentclass{article}\n\\usepackage{longtable,pdflscape}\n"
                                "\\begin{document}\n\\input{runtime/table.tex}\n"
                                "\\input{optimizations/table.tex}\n\\end{document}\n")
            result = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", document.name],
                                    cwd=temporary, capture_output=True, text=True)
            assert result.returncode == 0, result.stdout[-3000:]
    print("PASS: no-argument entry points, individual switches, real/synthetic inputs, epsilon sweep, repeats, resume, tables and limits")


if __name__ == "__main__":
    main()
