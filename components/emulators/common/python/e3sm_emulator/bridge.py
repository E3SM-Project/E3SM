"""The entry point the C++ backend calls.

``PythonBackend`` imports this module, calls :func:`create_emulator` once with
everything from the namelist plus the context the coupler supplied, and then
calls ``infer(inputs, outputs)`` on whatever comes back, once per timestep.

``inputs`` and ``outputs`` are dicts of numpy arrays that *view E3SM memory*.
Reading an input reads the component's own field; writing an output writes it.
Neither is copied, so an emulator writes its results with a slice
assignment::

    def infer(self, inputs, outputs):
        outputs["dT"][:] = self.model(inputs["T"])

Input arrays are read-only, so a bug in a model cannot corrupt component
state.  Everything is float64, because every E3SM field is ``real(r8)``;
casting to whatever precision a model was trained in is one call and belongs
on this side of the boundary.
"""

from __future__ import annotations

import importlib

from .context import Context

#: Emulator name -> the module that provides ``build(config, context)``.
_EMULATORS = {
    "example": "e3sm_emulator.example",
}


def create_emulator(config: dict):
    """Build the emulator named by ``config['emulator']``.

    Args:
        config: Every ``inference.*`` setting from the component namelist, as
            strings, plus ``model_path``, ``inputs``, ``outputs``, ``verbose``
            and a ``context`` dict.

    Returns:
        An object with ``infer(inputs, outputs)`` and, optionally,
        ``finalize()``.
    """
    context = Context.from_dict(config.get("context", {}))
    name = str(config.get("emulator", "example")).lower()

    if name not in _EMULATORS:
        raise ValueError(
            f"Unknown emulator '{name}'. Set `inference.emulator` to one of: "
            f"{', '.join(sorted(_EMULATORS))}. To run your own, point "
            "`inference.python_module` at a module of yours that provides "
            "create_emulator(config)."
        )

    emulator = importlib.import_module(_EMULATORS[name]).build(config, context)
    if config.get("verbose") and context.is_root:
        print(f"[e3sm_emulator] {name}: {context.describe()}", flush=True)
    return emulator
