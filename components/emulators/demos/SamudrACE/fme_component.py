"""Run one component of an fme coupled checkpoint through the Python backend.

Each `infer` call advances the component by one of its own timesteps. The caller
supplies every input and takes every output, so this module holds only what has
no place in the coupling plan: the corrector and random state, and the steps that
fme's coupled stepper applies between its components.

Options read from the emulator config, all strings:
- `fme_component`: `atmosphere` or `ocean`.
- `seed`: seeds the component's random state. Unset uses the global torch RNG.

Atmosphere inputs are the model's own, except that `OCNFRAC` and `ICEFRAC` are
derived from the ocean's `sst` and `ocean_sea_ice_fraction`, which are inputs here.
`SOLIN` must hold its value at the end of the step.

Ocean inputs are the model's own, with the atmosphere's fluxes averaged over the
ocean's timestep. They also include `atmosphere_flux_fraction` and
`momentum_flux_fraction`, which are used to scale those fluxes.

Other forcing is reused for the end of the step.
"""

from functools import lru_cache

import torch
from fme.core import get_device
from fme.core.random_state import RandomState
from fme.core.step.args import StepArgs
from fme.core.stepper_state import StepperState
from fme.coupled.inference.evaluator import load_stepper

TIME_DIM = 1


@lru_cache(maxsize=None)
def _load_coupled(path):
    """Both components share one checkpoint, which is too large to load twice."""
    coupled = load_stepper(path)
    coupled.set_eval()
    return coupled


def create_emulator(config):
    return FmeComponent(config)


class FmeComponent:
    def __init__(self, config):
        self._coupled = _load_coupled(config["model_path"])
        self._component = config["fme_component"]
        self._stepper = getattr(self._coupled, self._component)
        self._device = get_device()
        self._n_calls = 0
        self._state = None
        if "seed" in config:
            self._state = StepperState(
                random_state=RandomState.from_seed(int(config["seed"]))
            )

    def _to_tensor(self, values):
        # Copies: inputs are read-only views, and the model runs in float32
        tensor = torch.tensor(values, dtype=torch.float32, device=self._device)
        return tensor.unsqueeze(0)  # sample dimension

    def _prepare_atmosphere(self, step_input):
        """Return the inputs and the surface temperature to prescribe."""
        config = self._coupled._config
        fractions = config.ocean_fraction_prediction
        surface_temperature = self._stepper.surface_temperature_name
        from_ocean = self._coupled._forcings_from_ocean_with_ocean_fraction(
            {
                surface_temperature: step_input.pop(config.sst_name),
                fractions.sea_ice_fraction_name: step_input.pop(
                    fractions.sea_ice_fraction_name
                ),
            },
            step_input,
        )
        step_input[config.ocean_fraction_name] = from_ocean[config.ocean_fraction_name]
        sea_ice = fractions.sea_ice_fraction_name_in_atmosphere
        step_input[sea_ice] = from_ocean[sea_ice]

        # The ocean's state only changes between windows. Blending with the
        # ocean fraction isn't idempotent, so this mustn't repeat within a window.
        if self._n_calls % self._coupled.n_inner_steps == 0:
            step_input = self._stepper.prescribe_sst(
                mask_data=from_ocean, gen_data=step_input, target_data=from_ocean
            )
        return step_input, from_ocean[surface_temperature]

    def _prepare_ocean(self, step_input):
        """Apply the flux scaling from CoupledStepper._get_ocean_forcings."""
        config = self._coupled._config

        def pop_fraction(name):
            return torch.nan_to_num(step_input.pop(name), nan=1.0)

        flux_fraction = pop_fraction(config._atmosphere_flux_fraction_name)
        static_fraction = pop_fraction(config._static_flux_fraction_name)

        names = config.ocean_fraction_prediction
        ocean_state = {
            name: step_input[name].unsqueeze(TIME_DIM)
            for name in (config.sst_name, names.sea_ice_fraction_name)
        }
        land = {
            names.land_fraction_name: step_input[
                names.land_fraction_name
            ].unsqueeze(TIME_DIM)
        }
        open_water = self._coupled._get_open_water_fraction(ocean_state, land)
        step_input = config.open_water_flux_scaling.scale(
            step_input, open_water.squeeze(TIME_DIM) * flux_fraction
        )
        step_input = config.static_flux_scaling.scale(step_input, static_fraction)
        return step_input, None

    def _build_args(self, inputs):
        step_input = {name: self._to_tensor(values) for name, values in inputs.items()}
        prepare = (
            self._prepare_atmosphere
            if self._component == "atmosphere"
            else self._prepare_ocean
        )
        step_input, surface_temperature = prepare(step_input)

        step = self._stepper._step_obj
        next_step_input = {}
        for name in step.next_step_input_names:
            if name == step.surface_temperature_name:
                next_step_input[name] = surface_temperature
            else:
                next_step_input[name] = step_input[name]
        return StepArgs(
            input=step_input,
            next_step_input_data=next_step_input,
            stepper_state=self._state,
        )

    @torch.no_grad()
    def infer(self, inputs, outputs):
        result = self._stepper.step(self._build_args(inputs))
        self._state = result.stepper_state
        self._n_calls += 1
        for name, values in outputs.items():
            values[:] = result.output[name][0].double().cpu().numpy()
