from __future__ import annotations

import argparse
from pathlib import Path

from mbpi_core import (
    aggregate_scores,
    command_string,
    default_paths,
    ensure_out_dir,
    write_algorithm_spec,
    write_observed_outputs,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute observed MBPI components and audit outputs.")
    parser.add_argument("--root", default=".", help="Repository root.")
    parser.add_argument("--out-dir", default="results/mbpi", help="Output directory.")
    parser.add_argument(
        "--write-spec",
        action="store_true",
        help="Write MBPI_algorithm_specification.md using the supplied command string.",
    )
    parser.add_argument(
        "--pipeline-command",
        default="python scripts/mbpi/run_mbpi_pipeline.py",
        help="Exact reproducibility command to place in the technical specification.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    paths = default_paths(Path(args.root), Path(args.out_dir))
    ensure_out_dir(paths.out_dir)
    components, classifier_components = write_observed_outputs(paths)
    if args.write_spec:
        write_algorithm_spec(paths, args.pipeline_command)
    mbpi = aggregate_scores(components, include_domain=False)
    mbpi_plus = aggregate_scores(components, include_domain=True)
    print(f"MBPI observed: {mbpi:.6f}")
    print(f"MBPI_plus observed: {mbpi_plus:.6f}")
    print(f"Components written to {paths.out_dir / 'mbpi_components_observed.csv'}")
    if not classifier_components.empty:
        print(f"Classifier diagnostics written to {paths.out_dir / 'mbpi_classifier_components.csv'}")


if __name__ == "__main__":
    main()
