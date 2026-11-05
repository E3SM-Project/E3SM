"""Write the inputs, configs, and reference results for the SamudrACE demo.

Run in an environment with fme installed, from any directory:

    python prepare_demo.py --data-dir e3smv3 --fields fields.yaml --run-dir run

The reference results come from fme's own rollout of the atmosphere, with the
forcing built by fme's coupled stepper, so they show what the demo should produce.
"""

from __future__ import annotations

import argparse
import site
import sys
from pathlib import Path

import numpy as np
import torch
import xarray as xr
from fme.core.optimization import NullOptimization
from fme.core.random_state import RandomState
from fme.core.stepper_state import StepperState
from fme.coupled.inference.evaluator import load_stepper

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "config"))
from create_emulator_config_files import (
    CoupledField,
    CouplingPlan,
    EmulatorConfig,
    Field,
    Model,
    PrognosticField,
    emulator_config_from_field_list,
    prognostic_fields,
)

EMULATOR_NAME = "atm-emulator"
DIMENSIONS = ["lat", "lon"]
SURFACE_TEMPERATURE_SOURCE = "sst"


def write_field(path: Path, values: torch.Tensor) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    values.cpu().numpy().astype(np.float64).tofile(path)


def load_initial_state(path: Path, names: list[str]) -> dict[str, torch.Tensor]:
    with xr.open_dataset(path) as data:
        first = data.isel(time=0)
        return {
            name: torch.as_tensor(first[name].values, dtype=torch.float32)[None, None]
            for name in names
        }


def load_forcing(path: Path, names: list[str], n_times: int) -> dict[str, torch.Tensor]:
    """Return tensors of shape [sample, time, lat, lon] from the start of the file."""
    with xr.open_dataset(path) as data:
        window = data.isel(time=slice(0, n_times))
        n_lat, n_lon = data.sizes["lat"], data.sizes["lon"]
        forcing = {}
        for name in names:
            values = window[name]
            tensor = torch.as_tensor(values.values, dtype=torch.float32)
            if values.dims == ("time",):
                tensor = tensor[:, None, None].expand(n_times, n_lat, n_lon)
            elif "time" not in values.dims:
                tensor = tensor[None].expand(n_times, n_lat, n_lon)
            forcing[name] = tensor[None]
        return forcing


def write_configs(
    run_dir: Path, fields: Path, checkpoint: Path, options: dict[str, str]
) -> EmulatorConfig:
    atmosphere = emulator_config_from_field_list(
        fields,
        "atmosphere",
        name=EMULATOR_NAME,
        model=Model(backend="python", path=checkpoint.resolve(), options=options),
        dimensions=DIMENSIONS,
    )
    source = Field(name=SURFACE_TEMPERATURE_SOURCE, dimensions=DIMENSIONS)
    atmosphere.inputs = [
        *(item for item in atmosphere.inputs ),
        source,
    ]

    prognostic = prognostic_fields(atmosphere)
    prognostic_names = {item.name for item in prognostic}
    inputs_only = [i for i in atmosphere.inputs if i.name not in prognostic_names]
    outputs_only = [i for i in atmosphere.outputs if i.name not in prognostic_names]
    cpl7 = EmulatorConfig(
        name="cpl7",
        model=Model(backend="stub"),
        dimensions=[],
        inputs=[*prognostic, *outputs_only],
        outputs=[*prognostic, *inputs_only],
    )
    atmosphere.write(run_dir / "emulator_config.yaml")
    CouplingPlan(
        fields=[
            *(
                PrognosticField(field=item, component=atmosphere, external=cpl7)
                for item in prognostic
            ),
            *(
                CoupledField(
                    field=item, kind="state", sources=[cpl7], destinations=[atmosphere]
                )
                for item in inputs_only
            ),
            *(
                CoupledField(
                    field=item, kind="flux", sources=[atmosphere], destinations=[cpl7]
                )
                for item in outputs_only
            ),
        ]
    ).write(run_dir / "coupling_plan.yaml")
    return atmosphere


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--fields", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--flavor", default="C1A0")
    parser.add_argument("--year", type=int, default=1985)
    parser.add_argument("--n-steps", type=int, default=2)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    args.run_dir.mkdir(parents=True, exist_ok=True)

    checkpoint = args.data_dir / "checkpoints" / f"{args.flavor}.tar"
    coupled = load_stepper(checkpoint)
    atmosphere = coupled.atmosphere
    atmosphere.set_eval()
    if not 1 <= args.n_steps <= coupled.n_inner_steps:
        raise ValueError(f"--n-steps must be in [1, {coupled.n_inner_steps}]")

    ic_dir = args.data_dir / "initial_conditions"
    atmosphere_ic = load_initial_state(
        ic_dir / "restart-ic-1985-2015_atmosphere.nc", atmosphere.prognostic_names
    )
    ocean_ic = load_initial_state(
        ic_dir / "restart-ic-1985-2015_ocean.nc",
        coupled._config.ocean_to_atmosphere_forcing_names,
    )
    window = coupled._config.atmosphere_forcing_window_names
    forcing_file = args.data_dir / "forcing_data" / f"atmosphere-forcing-{args.year}.nc"
    with xr.open_dataset(forcing_file) as data:
        available = [name for name in window if name in data]
    forcing = load_forcing(forcing_file, available, coupled.n_inner_steps + 1)
    forcing = coupled._get_atmosphere_forcings(forcing, ocean_ic)

    step_obj = atmosphere._step_obj
    next_step_names = step_obj.next_step_forcing_names
    surface_temperature = step_obj.surface_temperature_name

    # Each step's inputs: forcing at the start, except fields from the end of the step
    for name in set(atmosphere.config.input_only_names):
        offset = 1 if name in next_step_names else 0
        for step in range(args.n_steps):
            write_field(
                args.run_dir / "inputs" / f"step_{step}" / f"{name}.bin",
                forcing[name][0, step + offset],
            )
    # The caller's surface temperature source is masked the way fme masks it
    write_field(
        args.run_dir / "inputs" / f"{SURFACE_TEMPERATURE_SOURCE}.bin",
        forcing[surface_temperature][0, 0],
    )
    for name, values in atmosphere_ic.items():
        write_field(args.run_dir / "initial_state" / f"{name}.bin", values[0, 0])

    options = {
        "python_module": "fme_component",
        "python_path": ":".join(
            [str(Path(__file__).resolve().parent), *site.getsitepackages()]
        ),
        "fme_component": "atmosphere",
        "seed": str(args.seed),
        "surface_temperature_source": SURFACE_TEMPERATURE_SOURCE,
    }
    write_configs(args.run_dir, args.fields, checkpoint, options)

    # Reference: fme's own rollout, using the true forcing at the end of each step
    with torch.no_grad():
        rollout = atmosphere.predict_generator(
            atmosphere_ic,
            forcing,
            args.n_steps,
            NullOptimization(),
            labels=None,
            stepper_state=StepperState(random_state=RandomState.from_seed(args.seed)),
        )
        for step, result in enumerate(rollout):
            for name, values in result.output.items():
                write_field(
                    args.run_dir / "reference" / f"step_{step}" / f"{name}.bin",
                    values[0],
                )
    (args.run_dir / "n_steps.txt").write_text(f"{args.n_steps}\n", encoding="utf-8")
    print(f"Wrote {args.run_dir}")


if __name__ == "__main__":
    main()
