"""Write the inputs, configs, and reference results for the SamudrACE demo.

Run in an environment with fme installed, from any directory:

    python prepare_demo.py --data-dir e3smv3 --fields fields.yaml --run-dir run

The reference results come from fme's own coupled rollout, so they show what the
demo should produce.
"""

from __future__ import annotations

import argparse
import site
import sys
from pathlib import Path

import numpy as np
import torch
import xarray as xr
from fme.ace.data_loading.batch_data import BatchData, PrognosticState
from fme.core.optimization import NullOptimization
from fme.core.rand import alternate_seed
from fme.coupled.data_loading.batch_data import CoupledBatchData, CoupledPrognosticState
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

DIMENSIONS = ["lat", "lon"]
ATMOSPHERE_NAME = "atm-emulator"
OCEAN_NAME = "ocn-emulator"
# The driver computes this one each step, instead of reading it from a file
COMPUTED_BY_DRIVER = {"global_mean_co2"}


def write_field(path: Path, values: torch.Tensor) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    values.cpu().numpy().astype(np.float64).tofile(path)


def load_initial_state(path: Path, names: list[str]) -> tuple[dict, xr.DataArray]:
    """Return tensors of shape [sample, time, lat, lon], and their time."""
    with xr.open_dataset(path) as data:
        first = data.isel(time=[0])
        state = {
            name: torch.as_tensor(first[name].values, dtype=torch.float32)[None]
            for name in names
        }
        return state, first["time"].load().expand_dims("sample", 0)


def load_forcing(path: Path, names: list[str], n_times: int) -> tuple[dict, xr.DataArray]:
    """Return tensors of shape [sample, time, lat, lon], and their time."""
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
        return forcing, window["time"].load().expand_dims("sample", 0)


def write_configs(
    run_dir: Path, fields: Path, checkpoint: Path, options: dict, seeds: dict, coupled
) -> None:
    config = coupled._config
    extra_field = lambda name: Field(name=name, dimensions=DIMENSIONS)

    def component_config(component, name, seed):
        return emulator_config_from_field_list(
            fields,
            component,
            name=name,
            model=Model(
                backend="python",
                path=checkpoint.resolve(),
                options={**options, "fme_component": component, "seed": str(seed)},
            ),
            dimensions=DIMENSIONS,
        )

    atmosphere = component_config("atmosphere", ATMOSPHERE_NAME, seeds["atmosphere"])
    ocean = component_config("ocean", OCEAN_NAME, seeds["ocean"])

    # The atmosphere gets the ocean's fields, and derives its fractions from them
    fractions = config.ocean_fraction_prediction
    derived = {config.ocean_fraction_name, fractions.sea_ice_fraction_name_in_atmosphere}
    ocean_to_atmosphere = [config.sst_name, fractions.sea_ice_fraction_name]
    atmosphere.inputs = [
        *(item for item in atmosphere.inputs if item.name not in derived),
        *(extra_field(name) for name in ocean_to_atmosphere),
    ]
    # The ocean also needs the fractions that scale the atmosphere's fluxes
    ocean.inputs = [
        *ocean.inputs,
        extra_field(config._atmosphere_flux_fraction_name),
        extra_field(config._static_flux_fraction_name),
    ]

    fluxes = set(config.atmosphere_to_ocean_forcing_names)
    atmosphere_prognostic = prognostic_fields(atmosphere)
    ocean_prognostic = prognostic_fields(ocean)
    prognostic_names = {i.name for i in (*atmosphere_prognostic, *ocean_prognostic)}
    routed_from_elsewhere = fluxes | set(ocean_to_atmosphere)

    def external_inputs(component):
        return [
            item
            for item in component.inputs
            if item.name not in prognostic_names | routed_from_elsewhere
        ]

    atmosphere_forcing = external_inputs(atmosphere)
    ocean_forcing = external_inputs(ocean)
    atmosphere_only = [i for i in atmosphere_forcing if i.name not in {f.name for f in ocean_forcing}]
    ocean_only = [i for i in ocean_forcing if i.name not in {f.name for f in atmosphere_forcing}]
    shared = [i for i in atmosphere_forcing if i.name in {f.name for f in ocean_forcing}]
    other_outputs = [
        item
        for item in atmosphere.outputs
        if item.name not in prognostic_names | fluxes
    ]

    cpl7 = EmulatorConfig(
        name="cpl7",
        model=Model(backend="stub"),
        dimensions=[],
        inputs=[*atmosphere_prognostic, *ocean_prognostic, *other_outputs],
        outputs=[
            *atmosphere_prognostic,
            *ocean_prognostic,
            *atmosphere_only,
            *ocean_only,
            *shared,
        ],
    )

    def field_named(items, name):
        return next(item for item in items if item.name == name)

    atmosphere.write(run_dir / "atmosphere_config.yaml")
    ocean.write(run_dir / "ocean_config.yaml")
    CouplingPlan(
        fields=[
            *(
                PrognosticField(field=item, component=atmosphere, external=cpl7)
                for item in atmosphere_prognostic
            ),
            *(
                PrognosticField(field=item, component=ocean, external=cpl7)
                for item in ocean_prognostic
            ),
            # Averaged, because the ocean runs less often than the atmosphere
            *(
                CoupledField(
                    field=field_named(atmosphere.outputs, name),
                    kind="flux",
                    sources=[atmosphere],
                    destinations=[ocean],
                )
                for name in sorted(fluxes)
            ),
            # Held between ocean steps. Separate from the fields' prognostic routes
            *(
                CoupledField(
                    field=field_named(ocean.outputs, name),
                    kind="state",
                    sources=[ocean],
                    destinations=[atmosphere],
                    route_id=f"{name}_to_atmosphere",
                )
                for name in ocean_to_atmosphere
            ),
            *(
                CoupledField(
                    field=item, kind="state", sources=[cpl7], destinations=[atmosphere]
                )
                for item in atmosphere_only
            ),
            *(
                CoupledField(
                    field=item, kind="state", sources=[cpl7], destinations=[ocean]
                )
                for item in ocean_only
            ),
            *(
                CoupledField(
                    field=item,
                    kind="state",
                    sources=[cpl7],
                    destinations=[atmosphere, ocean],
                )
                for item in shared
            ),
            *(
                CoupledField(
                    field=item, kind="flux", sources=[atmosphere], destinations=[cpl7]
                )
                for item in other_outputs
            ),
        ]
    ).write(run_dir / "coupling_plan.yaml")


def write_forcing(run_dir: Path, forcing: dict, directory: str, n_steps: int, offset):
    """Write fields that change once per step, and the others once."""
    for name, values in forcing.items():
        if name in COMPUTED_BY_DRIVER:
            continue
        if torch.equal(values[:, 0], values[:, -1]) and (values == values[:, :1]).all():
            write_field(run_dir / "inputs" / "fixed" / f"{name}.bin", values[0, 0])
            continue
        for step in range(n_steps):
            write_field(
                run_dir / "inputs" / f"{directory}_{step}" / f"{name}.bin",
                values[0, step + offset(name)],
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--fields", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--flavor", default="C1A0")
    parser.add_argument("--year", type=int, default=1985)
    parser.add_argument("--n-coupled-steps", type=int, default=2)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    args.run_dir.mkdir(parents=True, exist_ok=True)

    checkpoint = args.data_dir / "checkpoints" / f"{args.flavor}.tar"
    coupled = load_stepper(checkpoint)
    coupled.set_eval()
    config = coupled._config
    atmosphere = coupled.atmosphere
    ocean = coupled.ocean
    n_inner = coupled.n_inner_steps
    n_coupled = args.n_coupled_steps
    n_atmosphere = n_coupled * n_inner

    ic_dir = args.data_dir / "initial_conditions"
    atmosphere_ic, atmosphere_time = load_initial_state(
        ic_dir / "restart-ic-1985-2015_atmosphere.nc", atmosphere.prognostic_names
    )
    ocean_ic, ocean_time = load_initial_state(
        ic_dir / "restart-ic-1985-2015_ocean.nc", ocean.prognostic_names
    )
    forcing_dir = args.data_dir / "forcing_data"
    atmosphere_forcing, atmosphere_forcing_time = load_forcing(
        forcing_dir / f"atmosphere-forcing-{args.year}.nc",
        config.atmosphere_forcing_window_names,
        n_atmosphere + 1,
    )
    ocean_forcing, ocean_forcing_time = load_forcing(
        forcing_dir / f"ocean-forcing-{args.year}.nc",
        config.ocean_forcing_window_names,
        n_coupled + 1,
    )

    # Per-step inputs. Forcing is from the start of the step, except for fields
    # the model wants from the end of it.
    atmosphere_next = atmosphere._step_obj.next_step_forcing_names
    write_forcing(
        args.run_dir,
        atmosphere_forcing,
        "step",
        n_atmosphere,
        lambda name: 1 if name in atmosphere_next else 0,
    )
    write_forcing(
        args.run_dir, ocean_forcing, "coupled_step", n_coupled, lambda name: 0
    )
    for name, values in {**atmosphere_ic, **ocean_ic}.items():
        write_field(args.run_dir / "initial_state" / f"{name}.bin", values[0, 0])
    # The atmosphere runs first, so its routes from the ocean start with the
    # ocean's initial state
    fractions = config.ocean_fraction_prediction
    for name in (config.sst_name, fractions.sea_ice_fraction_name):
        write_field(
            args.run_dir / "route_initial_state" / f"{name}_to_atmosphere.bin",
            ocean_ic[name][0, 0],
        )

    options = {
        "python_module": "fme_component",
        # An embedded interpreter doesn't see a virtual environment's packages
        "python_path": ":".join(
            [str(Path(__file__).resolve().parent), *site.getsitepackages()]
        ),
    }
    # fme seeds the atmosphere with the seed, and the ocean with another one
    seeds = {"atmosphere": args.seed, "ocean": alternate_seed(args.seed)}
    write_configs(args.run_dir, args.fields, checkpoint, options, seeds, coupled)
    (args.run_dir / "n_coupled_steps.txt").write_text(f"{n_coupled}\n", encoding="utf-8")
    (args.run_dir / "n_inner_steps.txt").write_text(f"{n_inner}\n", encoding="utf-8")

    # Reference: fme's own coupled rollout
    initial_condition = CoupledPrognosticState(
        ocean_data=PrognosticState(BatchData(data=ocean_ic, time=ocean_time)),
        atmosphere_data=PrognosticState(BatchData(data=atmosphere_ic, time=atmosphere_time)),
    ).apply_config_seed(args.seed)
    forcing = CoupledBatchData(
        ocean_data=BatchData(data=ocean_forcing, time=ocean_forcing_time),
        atmosphere_data=BatchData(data=atmosphere_forcing, time=atmosphere_forcing_time),
    )
    with torch.no_grad():
        for prediction in coupled.get_prediction_generator(
            initial_condition, forcing, NullOptimization()
        ):
            for name, values in prediction.data.items():
                write_field(
                    args.run_dir
                    / "reference"
                    / prediction.realm
                    / f"step_{prediction.step}"
                    / f"{name}.bin",
                    values[0],
                )
    print(f"Wrote {args.run_dir}")


if __name__ == "__main__":
    main()
