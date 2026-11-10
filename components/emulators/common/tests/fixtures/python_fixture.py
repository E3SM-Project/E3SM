"""Python models the Python backend test runs. They need only numpy."""


def create_emulator(config):
    return Model(config)


def create_broken(config):
    raise RuntimeError("deliberate failure in the factory")


class Model:
    """y = scale * x + step, and the sum of x over its last dimension."""

    def __init__(self, config):
        if config["model_path"] != "fixture.ckpt":
            raise ValueError(f"unexpected model_path {config['model_path']!r}")
        # Options arrive as strings
        self.scale = float(config["scale"])
        self.steps = 0
        self.finalized = False

    def infer(self, inputs, outputs):
        x = inputs["x"]
        if x.flags.writeable:
            raise AssertionError("input 'x' is writeable")
        if x.dtype != "float64" or not x.flags.c_contiguous:
            raise AssertionError("input 'x' is not a C-ordered float64 array")
        if self.scale < 0:
            raise RuntimeError("deliberate failure in infer")

        self.steps += 1
        outputs["y"][:] = self.scale * x + self.steps
        outputs["row_sum"][:] = x.sum(axis=-1)

    def finalize(self):
        if self.finalized:
            raise AssertionError("finalize() was called twice")
        self.finalized = True
