from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from mbpi_core import command_string, default_paths, ensure_out_dir, write_algorithm_spec, write_run_log


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the complete MBPI reproducible pipeline.")
    parser.add_argument("--root", default=".", help="Repository root.")
    parser.add_argument("--out-dir", default="results/mbpi", help="Output directory.")
    parser.add_argument("--bootstrap-B", type=int, default=500, help="Bootstrap replicate count.")
    parser.add_argument(
        "--bootstrap-mode",
        choices=["fast", "full"],
        default="fast",
        help="Bootstrap recomputation mode.",
    )
    parser.add_argument("--simulation-reps", type=int, default=50, help="Replicates per simulation setting.")
    parser.add_argument("--seed", type=int, default=20260608, help="Random seed.")
    return parser.parse_args()


def run_step(cmd: list[str], cwd: Path) -> None:
    print("Running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()
    paths = default_paths(root, Path(args.out_dir))
    ensure_out_dir(paths.out_dir)
    pipeline_command = command_string(sys.argv)
    py = sys.executable
    steps = [
        [
            py,
            "scripts/mbpi/compute_mbpi.py",
            "--root",
            str(root),
            "--out-dir",
            str(paths.out_dir.relative_to(root)),
            "--write-spec",
            "--pipeline-command",
            pipeline_command,
        ],
        [
            py,
            "scripts/mbpi/bootstrap_mbpi.py",
            "--root",
            str(root),
            "--out-dir",
            str(paths.out_dir.relative_to(root)),
            "--B",
            str(args.bootstrap_B),
            "--mode",
            args.bootstrap_mode,
            "--seed",
            str(args.seed),
        ],
        [
            py,
            "scripts/mbpi/simulate_mbpi_calibration.py",
            "--root",
            str(root),
            "--out-dir",
            str(paths.out_dir.relative_to(root)),
            "--reps",
            str(args.simulation_reps),
            "--seed",
            str(args.seed),
        ],
        [
            py,
            "scripts/mbpi/ablate_mbpi.py",
            "--root",
            str(root),
            "--out-dir",
            str(paths.out_dir.relative_to(root)),
        ],
        [
            py,
            "scripts/mbpi/plot_mbpi_outputs.py",
            "--root",
            str(root),
            "--out-dir",
            str(paths.out_dir.relative_to(root)),
        ],
    ]
    warnings = []
    if args.bootstrap_mode == "fast":
        warnings.append(
            "Fast bootstrap recomputed CDDR and disease-axis angle from resampled CLR matrices; "
            "prediction, feature-stability, and genus-level heterogeneity components were held at observed summary values."
        )
    for step in steps:
        run_step(step, root)
    write_algorithm_spec(paths, pipeline_command)
    write_run_log(paths, pipeline_command, warnings=warnings)
    print(f"MBPI pipeline complete. Outputs: {paths.out_dir}", flush=True)


if __name__ == "__main__":
    main()
