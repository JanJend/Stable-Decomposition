"""End-size pipeline check: python3 -m experiments.support.check_end_sizes."""
import csv
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from . import DATA, HERE, ROOT
from .runtime import digest


def read_rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def check_seven_inputs(experiment, temporary):
    selected = experiment.select_inputs(DATA / "manifest.csv")
    assert [row["dataset"] for row in selected] == [
        "torus_gauss_640_1", "torus_gauss_1280_1",
        "T_E_ROI_40_locations_CD8", "T_A_ROI_36_locations_CD8",
        "gauss_clust_640_1", "gauss_clust_1280_1", "LargeHypoxicRegion"]
    assert [int(row["n_points"]) for row in selected] == [640, 1280, 640, 1281, 640, 1280, 15724]
    focused = experiment.select_inputs(DATA / "manifest.csv", only_640=True)
    assert [row["dataset"] for row in focused] == [
        "torus_gauss_640_1", "T_E_ROI_40_locations_CD8", "gauss_clust_640_1"]
    assert all(int(row["n_points"]) == 640 for row in focused)
    focused = experiment.select_inputs(DATA / "manifest.csv", immuno_points=640)
    assert [(row["dataset"], int(row["n_points"])) for row in focused] == [("T_E_ROI_40_locations_CD8", 640)]
    focused = experiment.select_inputs(DATA / "manifest.csv", immuno_points=1280)
    assert [(row["dataset"], int(row["n_points"])) for row in focused] == [("T_A_ROI_36_locations_CD8", 1281)]
    assert experiment.EPSILON_PERCENTAGES == [0.5, 1, 2]
    fixture = temporary / "inputs"
    fixture.mkdir()
    # Include distracting sizes, seeds and the second hypoxic region: selection
    # must still contain exactly the six requested clouds plus Region (not 2).
    records = []
    for family in ("torus_gauss", "immuno", "gauss_clust", "S2"):
        for n in (80, 320, 640, 1280):
            for seed in ((1,) if family == "immuno" else (1, 2, 3)):
                name = f"{family}_{n}_{seed}"
                path = fixture / (name + ".scc")
                path.write_text("scc2020\n2\n0 2 0\n0 2 ;\n2 0 ;\n")
                records.append(dict(dataset=name, family=family, target_points=n,
                                    n_points=1281 if family == "immuno" and n == 1280 else n,
                                    seed="" if family == "immuno" else seed,
                                    homology=0 if family == "gauss_clust" else 1,
                                    scc=path.name, scc_sha256=digest(path), generation_status="ok"))
    for name in ("LargeHypoxicRegion", "LargeHypoxicRegion2"):
        path = fixture / (name + ".scc")
        path.write_text("scc2020\n2\n0 2 0\n0 2 ;\n2 0 ;\n")
        records.append(dict(dataset=name, family="immuno", target_points=15724, n_points=15724,
                            seed="", homology=1, scc=path.name, scc_sha256=digest(path), generation_status="ok"))
    manifest = fixture / "manifest.csv"
    with manifest.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    worker = fixture / "fake-worker"
    worker.write_text(f"#!{sys.executable}\nimport sys\n"
                      "print('hom_algorithm=optimised')\nprint('end_size_pipeline=1')\n"
                      "print('--generators-only')\n"
                      "if '--help' in sys.argv: sys.exit(0)\n"
                      "print('End_0_dimension=3 seconds=0.1')\n"
                      "print('End_d_basis\\t8\\t8\\t0.1')\n"
                      "print('quotient_basis\\t5\\t5\\t0.1')\n"
                      "print('quotient_coordinates\\t5\\t5\\t0.1')\n"
                      "print('composition_generators\\t3\\t2\\t0.1')\n"
                      "if '--generators-only' in sys.argv: sys.exit(0)\n"
                      "print('restriction\\tI\\t1\\t3\\t2')\n"
                      "print('restriction\\tI\\t2\\t3\\t0')\n"
                      "print('restriction\\tK\\t1\\t2\\t0')\n"
                      "print('restriction\\tK\\t2\\t2\\t1')\n"
                      "print('restriction\\tK\\t3\\t2\\t0')\n"
                      "print('restriction_summary phase=I passes=2 max_retained=2 total_retained=2 phase_seconds=0.2 selection_seconds=0.1')\n"
                      "print('restriction_summary phase=K passes=3 max_retained=1 total_retained=1 phase_seconds=0.2 selection_seconds=0.1')\n")
    worker.chmod(0o755)
    output = temporary / "batch"
    subprocess.run([sys.executable, str(HERE / "benchmark_end_sizes.py"), "--manifest", str(manifest),
                    "--worker", str(worker), "--output-dir", str(output)], check=True)
    rows = read_rows(output / "results.csv")
    assert len(rows) == 21 and len(rows[0]) == 10
    for row in rows:
        assert [row[label] for label, _ in experiment.COLUMNS] == ["3", "8", "5", "3", "2", "2", "1"]
    assert len(read_rows(output / "restrictions.csv")) == 21 * 5
    logs = next(output.glob("*/settings.json")).parent
    settings = json.loads((logs / "settings.json").read_text())
    assert settings["points"] == [640, 1280] and settings["adaptive"] is True
    observations = [json.loads(line) for line in (logs / "observations.jsonl").read_text().splitlines()]
    assert len(observations) == 21 and all(row["status"] == "ok" for row in observations)
    assert all(row["timeout_seconds"] == (3600 if row["dataset"] == "LargeHypoxicRegion" else 600) for row in observations)
    for target, selection in ((640, ["--immuno-640-only"]), (1280, ["--immuno-only", "1280"])):
        focused_output = temporary / f"immune-{target}"
        subprocess.run([sys.executable, str(HERE / "benchmark_end_sizes.py"), "--manifest", str(manifest),
                        "--worker", str(worker), "--output-dir", str(focused_output), *selection], check=True)
        focused_rows = read_rows(focused_output / "results.csv")
        assert len(focused_rows) == 3 and {row["n"] for row in focused_rows} == {str(1281 if target == 1280 else 640)}
        assert {row["Point set"] for row in focused_rows} == {"immuno"}
        assert [float(row["epsilon (%)"]) for row in focused_rows] == [0.5, 1, 2]
        focused_settings = json.loads(next(focused_output.glob("*/settings.json")).read_text())
        assert focused_settings["points"] == [target] and focused_settings["families"] == ["immuno"]
        assert [row["dataset"] for row in focused_settings["inputs"]] == [f"immuno_{target}_1"]
    focused_output = temporary / "focused-640"
    subprocess.run([sys.executable, str(HERE / "benchmark_end_sizes.py"), "--manifest", str(manifest),
                    "--worker", str(worker), "--output-dir", str(focused_output), "--640-only",
                    "--families", "torus_gauss", "gauss_clust", "--generators-only"], check=True)
    focused_rows = read_rows(focused_output / "results.csv")
    assert len(focused_rows) == 6 and len(focused_rows[0]) == 8
    assert {row["Point set"] for row in focused_rows} == {"torus_gauss", "gauss_clust"}
    assert all(row["n"] == "640" and row["|B_l| (I)"] == "3" and row["|B_r| (K)"] == "2" for row in focused_rows)
    focused_settings = json.loads(next(focused_output.glob("*/settings.json")).read_text())
    assert focused_settings["points"] == [640] and focused_settings["generators_only"] and not focused_settings["adaptive"]
    assert focused_settings["families"] == ["torus_gauss", "gauss_clust"]
    assert not read_rows(focused_output / "restrictions.csv")
    assert (focused_output / "results.tex").read_text().count(r"\begin{table}") == 2
    assert (focused_output / "dimensions.tex").is_file() and (focused_output / "generators.tex").is_file()
    # An unfinished K phase must retain E,H,Q,B_l,B_r and completed I selection.
    partial = experiment.parse_output("End_0_dimension=3 seconds=0.1\nEnd_d_basis\t8\t8\t0.1\n"
                                      "quotient_basis\t5\t5\t0.1\ncomposition_generators\t3\t2\t0.1\n"
                                      "restriction\tI\t1\t3\t2\n"
                                      "restriction_summary phase=I passes=1 max_retained=2 total_retained=2 phase_seconds=0.2 selection_seconds=0.1\n")
    partial.update(input=settings["inputs"][0]["input"], epsilon_percent=1, status="timeout")
    partial_rows = experiment.report_rows([partial], settings["inputs"], [1, 8])
    assert partial_rows[0]["I_restriction_max"] == 2 and partial_rows[0]["K_restriction_max"] == "timeout"
    assert partial_rows[0]["End_0_dimension"] == 3 and partial_rows[1]["End_0_dimension"] == "pending"
    tex = (output / "results.tex").read_text()
    assert tex.count(r"\begin{table}") == tex.count(r"\begin{tabular}") == 1
    assert sum(" & " in line for line in tex.splitlines()) == 22
    assert not any(token in tex for token in (r"\begin{landscape}", r"\resizebox", r"\geometry", "Seed"))
    assert "max" in tex and "not just the terminal pass" in tex
    if shutil.which("pdflatex"):
        (output / "tables.tex").write_text("\\documentclass{article}\n"
                                           "\\begin{document}\\input{results.tex}\\end{document}\n")
        result = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "tables.tex"],
                                cwd=output, capture_output=True, text=True)
        assert result.returncode == 0 and "Overfull" not in result.stdout, result.stdout[-4000:]


def main():
    sys.path.insert(0, str(HERE))
    try:
        import benchmark_end_sizes as experiment
    finally:
        sys.path.pop(0)
    with tempfile.TemporaryDirectory(prefix="end-funnel-check-") as name:
        temporary = Path(name)
        check_seven_inputs(experiment, temporary)
        source = temporary / "incomparable generators.scc"
        source.write_text("scc2020\n2\n0 2 0\n0 2 ;\n2 0 ;\n")
        output = temporary / "known-dimensions"
        output.mkdir()
        (output / "results.csv").write_text("old_column\nold_result\n")
        command = [sys.executable, str(HERE / "benchmark_end_sizes.py"), "--input", str(source), "--output-dir", str(output)]
        subprocess.run([*command, "--percentages", "0", "0.5", "50"], check=True)
        rows = read_rows(output / "results.csv")
        assert [row["dim End(X)_0"] for row in rows] == ["2", "2", "2"]
        assert [row["dim End(X)_(2eps)"] for row in rows] == ["2", "2", "4"]
        assert [row["dim Q"] for row in rows] == ["0", "0", "2"]
        assert [row["|B_l| (I)"] for row in rows] == ["0", "0", "1"]
        assert [row["|B_r| (K)"] for row in rows] == ["0", "0", "1"]
        assert [row["I restricted max/pass"] for row in rows] == ["0", "0", "1"]
        assert [row["K restricted max/pass"] for row in rows] == ["0", "0", "0"]
        passes = read_rows(output / "restrictions.csv")
        assert [int(row["retained"]) for row in passes if row["phase"] == "I"] == [1, 0]
        assert [int(row["retained"]) for row in passes if row["phase"] == "K"] == [0]
        assert all("seconds" not in row for row in passes), "Per-pass timing must default to off"
        assert not json.loads(next(output.glob("*/settings.json")).read_text())["time_steps"]
        timed_output = temporary / "timed-restrictions"
        subprocess.run([sys.executable, str(HERE / "benchmark_end_sizes.py"), "--input", str(source),
                        "--output-dir", str(timed_output), "--percentages", "50", "--time-steps"], check=True)
        assert read_rows(timed_output / "results.csv") == [rows[-1]], "Timing must preserve the reduction counts"
        timed_passes = read_rows(timed_output / "restrictions.csv")
        assert [{key: value for key, value in row.items() if key != "seconds"} for row in timed_passes] == passes
        assert all(math.isfinite(float(row["seconds"])) and float(row["seconds"]) >= 0 for row in timed_passes)
        assert any(float(row["seconds"]) > 0 for row in timed_passes)
        timed_settings = json.loads(next(timed_output.glob("*/settings.json")).read_text())
        assert timed_settings["time_steps"] is True
        timed_totals = read_rows(timed_output / "timings.csv")[0]
        for phase in ("I", "K"):
            total = sum(float(row["seconds"]) for row in timed_passes if row["phase"] == phase)
            assert total <= float(timed_totals[f"{phase} restriction seconds"]) + 1e-6
        generators_output = temporary / "known-generators"
        subprocess.run([sys.executable, str(HERE / "benchmark_end_sizes.py"), "--input", str(source),
                        "--output-dir", str(generators_output), "--percentages", "50", "--generators-only"], check=True)
        generator_row = read_rows(generators_output / "results.csv")[0]
        assert [generator_row[label] for label, _ in experiment.COLUMNS[:5]] == ["2", "4", "2", "1", "1"]
        assert not read_rows(generators_output / "restrictions.csv")
        generator_logs = next(generators_output.glob("*/settings.json")).parent
        stdout = next(generator_logs.glob("*.stdout.txt")).read_text()
        assert "composition_generators\t1\t1\t" in stdout and "restriction" not in stdout
        logs = next(output.glob("*/settings.json")).parent
        assert (logs / "previous_results.csv").read_text() == "old_column\nold_result\n"
        for flags, count in ((["--no-nilpotent-reduction"], "1"),
                             (["--no-complement-orbit-generators"], "1"),
                             (["--no-nilpotent-reduction", "--no-complement-orbit-generators"], "2")):
            subprocess.run([*command, "--percentages", "50", *flags], check=True)
            row = read_rows(output / "results.csv")[0]
            assert row["|B_l| (I)"] == row["|B_r| (K)"] == count
        subprocess.run([*command, "--percentages", "50", "--composition-orbit-rounds", "6",
                        "--composition-orbit-seed-mode", "random"], check=True)
        latest = max(output.glob("*/settings.json"))
        assert json.loads(latest.read_text())["generator_settings"]["composition-orbit-rounds"] == 6
        timed_out = subprocess.run([*command, "--timeout", "0.000001"])
        assert timed_out.returncode == 1
        assert all(row["dim End(X)_0"] == "timeout" for row in read_rows(output / "results.csv"))
        assert len(read_rows(logs / "results.csv")) == 3, "Preserve the first run"
    print("PASS: default and focused 640-point inputs, family filters, generator-only tables, E/H/Q dimensions, left/right generators, actual I/K pass restrictions, optional per-pass timing, partial output, history and controls")


if __name__ == "__main__":
    main()
