"""Write the fields of an fme coupled checkpoint to a YAML field list.

Run this in an environment with fme installed:

    python extract_fme_fields.py checkpoint.tar -o fields.yaml

Pass the field list to `emulator_config_from_field_list` in
create_emulator_config_files.py, which needs neither fme nor torch.

Inputs and outputs keep the order of the model's `in_names` and `out_names`.
Don't use the stepper's `input_names` and `output_names` instead: they are built
from sets, so their order changes between Python processes.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import torch
import yaml
from fme.coupled.stepper import CoupledStepperConfig


def component_fields(state: dict, forcing: list[str]) -> dict:
    step = state["stepper"]["step"]
    if step.get("type") != "single_module":
        raise ValueError(f"Unsupported fme step type: {step.get('type')!r}")
    return {
        "timestep": state["timedelta"],
        "inputs": list(step["config"]["in_names"]),
        "outputs": list(step["config"]["out_names"]),
        "forcing": sorted(forcing),
    }


def extract_fields(checkpoint: str | Path) -> dict:
    # weights_only=False unpickles arbitrary objects: only use trusted checkpoints.
    # fme loads its checkpoints the same way.
    state = torch.load(checkpoint, map_location="cpu", weights_only=False)
    config_state = state["stepper"]["config"]
    if "atmosphere" not in config_state or "ocean" not in config_state:
        raise ValueError(f"{checkpoint} is not a coupled fme checkpoint")

    atmosphere = config_state["atmosphere"]
    ocean = config_state["ocean"]
    config = CoupledStepperConfig.from_state(config_state)
    return {
        "source": Path(checkpoint).name,
        "components": {
            "atmosphere": component_fields(
                atmosphere, config.atmosphere_forcing_exogenous_names
            ),
            "ocean": component_fields(ocean, config.ocean_forcing_exogenous_names),
        },
        "exchanges": {
            "atmosphere_to_ocean": sorted(config.atmosphere_to_ocean_forcing_names),
            "ocean_to_atmosphere": sorted(config.ocean_to_atmosphere_forcing_names),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("checkpoint", type=Path, help="fme coupled checkpoint")
    parser.add_argument("-o", "--output", type=Path, required=True, help="YAML output")
    args = parser.parse_args()

    fields = extract_fields(args.checkpoint)
    args.output.write_text(yaml.safe_dump(fields, sort_keys=False), encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
