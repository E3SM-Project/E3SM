"""A Python emulator with no dependencies beyond numpy, for the C++ tests.

Exercises the whole bridge -- the factory contract, the context payload, the
zero-copy views -- without needing torch or a checkpoint.
"""

import numpy as np


def create_emulator(config):
    return FixtureEmulator(config)


class FixtureEmulator:
    def __init__(self, config):
        self.config = config
        self.context = config.get("context", {})
        self.steps = 0
        # Written out so the C++ side can check that the factory ran and that
        # what it received survived the trip.
        with open(config.get("report_path", "emulator_fixture_report.txt"), "w") as f:
            f.write(f"model_path={config.get('model_path', '')}\n")
            f.write(f"scale={config.get('scale', '')}\n")
            f.write(f"rank={self.context.get('rank')}\n")
            f.write(f"world_size={self.context.get('world_size')}\n")
            f.write(f"nx={self.context.get('nx')} ny={self.context.get('ny')}\n")
            gids = np.asarray(self.context.get("col_gids", [])).tolist()
            f.write(f"gids={gids}\n")

    def infer(self, inputs, outputs):
        self.steps += 1
        scale = float(self.config.get("scale", 2.0))

        # Inputs must arrive as read-only views, so that a bug in a model
        # cannot corrupt component state.  Fail loudly if they do not: a
        # passing infer() is then proof, not an assumption.
        x = inputs["x"]
        if x.flags.writeable:
            raise AssertionError("input 'x' arrived writeable")
        if not outputs["y"].flags.writeable:
            raise AssertionError("output 'y' arrived read-only")

        # Slice assignment writes straight into the component's own array.
        outputs["y"][:] = scale * x + self.steps

    def finalize(self):
        pass
