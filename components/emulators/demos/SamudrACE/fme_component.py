"""Run one component of an fme coupled checkpoint through the Python backend.

Each `infer` call advances the component by one of its own timesteps. The caller
supplies every input and takes every output, so this module holds only what has
no place in the coupling plan: the corrector and random state.

Options read from the emulator config, all strings:
- `fme_component`: `atmosphere` or `ocean`.
- `seed`: seeds the component's random state. Unset uses the global torch RNG.
- `surface_temperature_source`: input field holding the surface temperature to
  prescribe at the end of the step, such as `sst`. It isn't passed to the network.

Fields in the step's `next_step_forcing_names`, such as `SOLIN`, must hold their
values at the end of the step. Other forcing is reused for the end of the step.
"""

import torch
from fme.core import get_device
from fme.core.random_state import RandomState
from fme.core.step.args import StepArgs
from fme.core.stepper_state import StepperState
from fme.coupled.inference.evaluator import load_stepper


def create_emulator(config):
    return FmeComponent(config)


class FmeComponent:
    def __init__(self, config):
        coupled = load_stepper(config["model_path"])
        self._stepper = getattr(coupled, config["fme_component"])
        self._stepper.set_eval()
        self._device = get_device()
        self._surface_temperature_source = config.get("surface_temperature_source")
        self._state = None
        if "seed" in config:
            self._state = StepperState(
                random_state=RandomState.from_seed(int(config["seed"]))
            )

    def _to_tensor(self, values):
        tensor = torch.tensor(values, dtype=torch.float32, device=self._device)
        return tensor.unsqueeze(0)

    def _build_args(self, inputs):
        step_input = {name: self._to_tensor(values) for name, values in inputs.items()}
        shape = next(iter(step_input.values())).shape

        step = self._stepper._step_obj
        surface_temperature = None
        if self._surface_temperature_source is not None:
            surface_temperature = step_input.pop(self._surface_temperature_source)

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
        for name, values in outputs.items():
            values[:] = result.output[name][0].double().cpu().numpy()
