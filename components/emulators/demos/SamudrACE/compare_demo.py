"""Compare the SamudrACE demo's outputs with fme's own rollout.

    python compare_demo.py run

Differences are reported relative to each field's standard deviation. The first
step should match exactly, so this exits with an error if it doesn't. Later steps
differ slightly because the demo holds global_mean_co2 constant.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

FIRST_STEP_TOLERANCE = 1e-6


def compare_step(step_dir: Path, reference_dir: Path) -> tuple[int, int, list]:
    exact = 0
    differences = []
    outputs = sorted(step_dir.glob("*.bin"))
    if not outputs:
        raise FileNotFoundError(f"No outputs in {step_dir}")
    for path in outputs:
        actual = np.fromfile(path)
        reference = np.fromfile(reference_dir / path.name)
        difference = np.abs(actual - reference).max()
        exact += difference == 0
        scale = reference.std() or 1.0
        differences.append((difference / scale, path.stem))
    return exact, len(outputs), sorted(differences, reverse=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()

    step_dirs = sorted(
        (args.run_dir / "outputs").iterdir(),
        key=lambda path: int(path.name.removeprefix("step_")),
    )
    first_step_worst = None
    for step_dir in step_dirs:
        exact, total, differences = compare_step(
            step_dir, args.run_dir / "reference" / step_dir.name
        )
        worst, name = differences[0]
        print(
            f"{step_dir.name}: {exact}/{total} fields identical; "
            f"largest difference {worst:.2e} of a standard deviation ({name})"
        )
        if first_step_worst is None:
            first_step_worst = worst

    if first_step_worst > FIRST_STEP_TOLERANCE:
        sys.exit(
            f"ERROR: the first step differs from fme by {first_step_worst:.2e}, "
            f"more than the tolerance of {FIRST_STEP_TOLERANCE:.0e}"
        )


if __name__ == "__main__":
    main()
