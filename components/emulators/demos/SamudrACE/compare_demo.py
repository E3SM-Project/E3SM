"""Compare the SamudrACE demo's outputs with fme's own coupled rollout.

    python compare_demo.py run

Differences are reported relative to each field's standard deviation. NaNs in the
same places count as equal. This exits with an error if either of these fails:

- The first window of atmosphere steps should match exactly, because it doesn't
  depend on the ocean's output yet.
- The first ocean step and the first atmosphere step after it should be very close.
  They can't match exactly: this package averages the fluxes in double precision,
  and fme averages them, and the constant LANDFRAC, in float32. Differences from that
  rounding then grow over time, so later steps aren't checked.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

FIRST_WINDOW_TOLERANCE = 1e-6
FIRST_COUPLING_TOLERANCE = 1e-3


def compare_step(step_dir: Path, reference_dir: Path) -> tuple[int, int, float, str]:
    exact = 0
    worst, worst_name = 0.0, ""
    outputs = sorted(step_dir.glob("*.bin"))
    if not outputs:
        raise FileNotFoundError(f"No outputs in {step_dir}")
    for path in outputs:
        actual = np.fromfile(path)
        reference = np.fromfile(reference_dir / path.name)
        if np.array_equal(actual, reference, equal_nan=True):
            exact += 1
            continue
        same_nans = np.array_equal(np.isnan(actual), np.isnan(reference))
        difference = (
            np.nanmax(np.abs(actual - reference)) / (np.nanstd(reference) or 1.0)
            if same_nans
            else np.inf
        )
        if difference > worst:
            worst, worst_name = difference, path.stem
    return exact, len(outputs), worst, worst_name


def step_directories(directory: Path) -> list[Path]:
    return sorted(
        directory.iterdir(), key=lambda path: int(path.name.removeprefix("step_"))
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()

    n_inner_steps = int((args.run_dir / "n_inner_steps.txt").read_text())
    first_window_worst = 0.0
    first_coupling_worst = 0.0
    for component in ("atmosphere", "ocean"):
        steps = step_directories(args.run_dir / "outputs" / component)
        print(f"{component}: {len(steps)} step(s)")
        for number, step_dir in enumerate(steps):
            exact, total, worst, name = compare_step(
                step_dir, args.run_dir / "reference" / component / step_dir.name
            )
            detail = f"; largest difference {worst:.2e} ({name})" if name else ""
            print(f"  {step_dir.name}: {exact}/{total} identical{detail}")
            if component == "atmosphere" and number < n_inner_steps:
                first_window_worst = max(first_window_worst, worst)
            if number == (n_inner_steps if component == "atmosphere" else 0):
                first_coupling_worst = max(first_coupling_worst, worst)

    if first_window_worst > FIRST_WINDOW_TOLERANCE:
        sys.exit(
            f"ERROR: the first window differs from fme by {first_window_worst:.2e}, "
            f"more than the tolerance of {FIRST_WINDOW_TOLERANCE:.0e}"
        )
    if first_coupling_worst > FIRST_COUPLING_TOLERANCE:
        sys.exit(
            f"ERROR: the first steps after coupling differ from fme by "
            f"{first_coupling_worst:.2e}, more than the tolerance of "
            f"{FIRST_COUPLING_TOLERANCE:.0e}"
        )


if __name__ == "__main__":
    main()
