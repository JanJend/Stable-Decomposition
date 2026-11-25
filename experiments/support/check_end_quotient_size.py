"""Run with python3 -m experiments.support.check_end_quotient_size."""
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from . import DATA, HERE, ROOT


def main():
    sys.path.insert(0, str(HERE))
    import end_quotient_size as experiment

    selected = experiment.select_inputs(DATA / "manifest.csv")
    assert [row["dataset"] for row in selected] == [
        "torus_gauss_320_1", "gauss_clust_320_1", "T_D_ROI_50_locations_FoxP3",
        "torus_gauss_640_1", "gauss_clust_640_1", "T_E_ROI_40_locations_CD8",
        "torus_gauss_1280_1", "gauss_clust_1280_1", "T_A_ROI_36_locations_CD8"]
    with tempfile.TemporaryDirectory(prefix="end-quotient-check-") as name:
        temporary = Path(name)
        records = []
        for n in (80, 320, 640, 1280):
            for family, _, homology in experiment.FAMILIES:
                for seed in (("",) if family == "immuno" else ("1", "2")):
                    path = temporary / f"{family}_{n}_{seed}.scc"
                    path.write_text("scc2020\n2\n0 2 0\n0 2 ;\n2 0 ;\n")
                    records.append(dict(dataset=path.stem, family=family, target_points=n,
                                        n_points=1281 if n == 1280 and family == "immuno" else n,
                                        homology=homology, seed=seed, scc=path.name, generation_status="ok"))
        manifest = temporary / "manifest.csv"
        with manifest.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)
        worker = temporary / "fake-worker"
        worker.write_text(f"#!{sys.executable}\n"
                          "import pathlib, sys\n"
                          "print('hom_algorithm=optimised')\n"
                          "if '--help' in sys.argv:\n"
                          "    print('--dimensions-only'); sys.exit(0)\n"
                          "assert '--dimensions-only' in sys.argv and '--percent' in sys.argv\n"
                          "n = int(pathlib.Path(sys.argv[1]).stem.split('_')[-2])\n"
                          "family = pathlib.Path(sys.argv[1]).stem.rsplit('_', 2)[0]\n"
                          "base = n + ['torus_gauss', 'gauss_clust', 'immuno'].index(family)\n"
                          "delta = int(float(sys.argv[2]) * 10)\n"
                          "print(f'End_0_dimension={base} seconds=0.1')\n"
                          "print(f'End_d_basis\\t{base+delta}\\t{base+delta}\\t0.1')\n"
                          "print(f'quotient_basis\\t{delta}\\t{delta}\\t0.1')\n")
        worker.chmod(0o755)
        output = temporary / "results"
        output.mkdir()
        (output / "results_320.csv").write_text("previous result\n")
        command = [sys.executable, str(HERE / "end_quotient_size.py"), "--manifest", str(manifest),
                   "--worker", str(worker), "--output-dir", str(output)]
        subprocess.run(command, check=True, capture_output=True, text=True)
        logs = next(output.glob("*/settings.json")).parent
        assert (logs / "previous_results_320.csv").read_text() == "previous result\n"
        settings = json.loads((logs / "settings.json").read_text())
        assert settings["points"] == [320, 640, 1280] and settings["percentages"] == [0.5, 1, 2]
        assert settings["dimensions_only"] and len(settings["inputs"]) == 9
        observations = [json.loads(line) for line in (logs / "observations.jsonl").read_text().splitlines()]
        assert len(observations) == 27 and all(row["status"] == "ok" for row in observations)
        for n in experiment.POINT_COUNTS:
            with (output / f"results_{n}.csv").open(newline="") as stream:
                rows = list(csv.reader(stream))
            assert len(rows) == 4 and all(len(row) == 10 for row in rows)
            assert [row[0] for row in rows[1:]] == [label for label, _, _ in experiment.DIMENSIONS]
            bases = [n + family for family in range(3) for _ in range(3)]
            deltas = [5, 10, 20] * 3
            assert rows[1][1:] == list(map(str, bases))
            assert rows[2][1:] == [str(base + delta) for base, delta in zip(bases, deltas)]
            assert rows[3][1:] == list(map(str, deltas))
        tex = (output / "results.tex").read_text()
        assert tex.count(r"\begin{table}") == tex.count(r"\begin{tabular}") == 3
        assert "1281 cells" in tex
        if shutil.which("pdflatex"):
            (output / "tables.tex").write_text("\\documentclass{article}\n"
                                               "\\begin{document}\\input{results.tex}\\end{document}\n")
            result = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "tables.tex"],
                                    cwd=output, capture_output=True, text=True)
            assert result.returncode == 0 and "Overfull" not in result.stdout, result.stdout[-4000:]
        partial = dict(observations[0], status="timeout")
        del partial["End_d_basis_I_size"], partial["quotient_basis_I_size"]
        experiment.write_reports(output, [partial], settings["inputs"])
        with (output / "results_320.csv").open(newline="") as stream:
            rows = list(csv.reader(stream))
        assert rows[1][1:3] == ["320", "pending"] and rows[2][1:3] == ["timeout", "pending"]
        # Real worker: two incomparable free generators give a nontrivial quotient
        # at 50%; even comparison/adaptive flags must not run any later stages.
        real_worker = ROOT / "build/local-release/hom_reduction_benchmark"
        result = subprocess.run([str(real_worker), str(temporary / "torus_gauss_320_1.scc"), "0", "50",
                                 "--percent", "--dimensions-only", "--compare-methods", "--adaptive"],
                                check=True, capture_output=True, text=True, timeout=30)
        assert "End_0_dimension=2" in result.stdout
        assert "End_d_basis\t2\t2\t" in result.stdout and "End_d_basis\t4\t4\t" in result.stdout
        assert "quotient_basis\t0\t0\t" in result.stdout and "quotient_basis\t2\t2\t" in result.stdout
        assert not any(stage in result.stdout for stage in
                       ("quotient_coordinates", "composition_generators", "restriction", "adaptive_phase", "full_"))
        timed_out = subprocess.run([*command, "--timeout", "0.000001"], capture_output=True, text=True)
        assert timed_out.returncode == 1
        with (output / "results_320.csv").open(newline="") as stream:
            assert all(value == "timeout" for row in list(csv.reader(stream))[1:] for value in row[1:])
    print("PASS: nine inputs, 27 measurements, three 3x9 tables, archived logs, partial results, timeouts and dimension-only worker")


if __name__ == "__main__":
    main()
