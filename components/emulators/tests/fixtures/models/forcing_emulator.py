# forcing_emulator.py
import numpy as np


def create_emulator(config):
    return ForcingEmulator(config)


class ForcingEmulator:
    def __init__(self, config):
        pass

    def infer(self, inputs, outputs):
        temperature = outputs["temperature"]
        soil_moisture = outputs["soil_moisture"]

        if temperature.size == 0:
            raise ValueError("ForcingEmulator has no temperature samples")

        ncols = temperature.size
        if soil_moisture.size % ncols != 0:
            raise ValueError(
                "soil_moisture size is not divisible by the number of columns"
            )

        nlevels = soil_moisture.size // ncols

        columns = np.arange(ncols, dtype=np.float64)
        levels = np.arange(nlevels, dtype=np.float64)

        temperature[...] = (280.0 + columns).reshape(temperature.shape)

        moisture = (
            0.10
            + 0.01 * columns[:, None]
            + 0.001 * levels[None, :]
        )
        soil_moisture[...] = moisture.reshape(soil_moisture.shape)
