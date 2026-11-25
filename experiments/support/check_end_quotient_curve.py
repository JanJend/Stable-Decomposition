"""Run with MPLBACKEND=Agg python -m experiments.support.check_end_quotient_curve."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

from . import DATA, HERE


def main():
    sys.path.insert(0, str(HERE))
    import end_quotient_curve as experiment

    selected = experiment.select_inputs(DATA / "manifest.csv", point_counts=[320])
    assert [row["dataset"] for row in selected] == [
        "torus_gauss_320_1", "gauss_clust_320_1", "T_D_ROI_50_locations_FoxP3"]
    selected = experiment.select_inputs(DATA / "manifest.csv", point_counts=[640])
    assert [row["dataset"] for row in selected] == [
        "torus_gauss_640_1", "gauss_clust_640_1", "T_E_ROI_40_locations_CD8"]
    with tempfile.TemporaryDirectory(prefix="end-curve-check-") as name:
        temporary = Path(name)
        samples = []
        for family, _, homology in experiment.FAMILIES:
            for n in (320, 640):
                path = temporary / f"{family}_{n}.scc"
                path.write_text("scc2020\n2\n0 2 0\n0 2 ;\n2 0 ;\n")
                samples.append(dict(dataset=path.stem, family=family, target_points=n, n_points=n,
                                    homology=homology, seed="1" if family != "immuno" else "",
                                    scc=path.name, generation_status="ok"))
        manifest = temporary / "manifest.csv"
        with manifest.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(samples[0]))
            writer.writeheader()
            writer.writerows(samples)
        worker = temporary / "fake-worker"
        worker.write_text(f"#!{sys.executable}\n"
                          "import sys\n"
                          "print('hom_algorithm=optimised')\n"
                          "if '--help' in sys.argv:\n"
                          "    print('--dimensions-only'); sys.exit(0)\n"
                          "assert '--percent' in sys.argv and '--dimensions-only' in sys.argv\n"
                          "delta = round(float(sys.argv[2]) * 10)\n"
                          "print('End_0_dimension=100 seconds=0.1')\n"
                          "print(f'End_d_basis\\t{100+delta}\\t{100+delta}\\t0.1')\n"
                          "if 'gauss_clust' in sys.argv[1] and delta == 25: sys.exit(1)\n"
                          "print(f'quotient_basis\\t{delta}\\t{delta}\\t0.1')\n")
        worker.chmod(0o755)
        output = temporary / "results"
        output.mkdir()
        (output / "results.csv").write_text("previous results\n")
        command = [sys.executable, str(HERE / "end_quotient_curve.py"), "--manifest", str(manifest),
                   "--worker", str(worker), "--output-dir", str(output), "--no-plot"]
        result = subprocess.run(command, capture_output=True, text=True, timeout=60)
        assert result.returncode == 1, result.stdout + result.stderr
        with (output / "results.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == 153
        for family, _, _ in experiment.FAMILIES:
            values = [row for row in rows if row["family"] == family]
            assert [row["epsilon_percent"] for row in values] == [f"{i / 10:.1f}" for i in range(51)]
            for i, row in enumerate(values):
                assert all(math.isclose(float(row[key]), i / 500) for key in ("epsilon_x", "epsilon_y"))
                assert row["End_0_dimension"] == "100" and row["End_2epsilon_dimension"] == str(100 + i)
                if family == "gauss_clust" and i == 25:
                    assert row["status"] == "error" and not row["quotient_dimension"]
                else:
                    assert row["status"] == "ok" and row["quotient_dimension"] == str(i)
        logs = next(output.glob("*/settings.json")).parent
        assert (logs / "results.csv").read_bytes() == (output / "results.csv").read_bytes()
        assert (logs / "previous_results.csv").read_text() == "previous results\n"
        settings = json.loads((logs / "settings.json").read_text())
        assert settings["points"] == [320] and settings["percentages"] == experiment.EPSILON_PERCENTAGES
        assert len((logs / "observations.jsonl").read_text().splitlines()) == 153
        import matplotlib.pyplot as plt
        from unittest.mock import patch

        with patch.object(plt, "show") as show:
            figure = experiment.plot_csv(output / "results.csv")
            show.assert_called_once()
        assert len(figure.axes) == 3 and len(figure.legends) == 1
        assert all(len(axis.lines) == 2 and len(axis.lines[0].get_xdata()) == 51 for axis in figure.axes)
        assert math.isnan(figure.axes[1].lines[1].get_ydata()[25])
        figure.canvas.draw()
        plt.close(figure)
        assert not any(output.glob("*.png")) and not any(output.glob("*.pdf"))
        result = subprocess.run([*command, "--points", "640"], capture_output=True, text=True, timeout=60)
        assert result.returncode == 1, result.stdout + result.stderr
        with (output / "results.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == 153 and {row["n_points"] for row in rows} == {"640"}
        assert {row["dataset"] for row in rows} == {f"{family}_640" for family, _, _ in experiment.FAMILIES}
        settings = json.loads(sorted(output.glob("*/settings.json"))[-1].read_text())
        assert settings["points"] == [640]
    print("PASS: 320/640-point selections, 153 rows per sweep, endpoints, partial failures, archives and two curves per panel")


if __name__ == "__main__":
    main()
