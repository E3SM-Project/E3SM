import numpy as np


class LandEmulator:
    def __init__(self, config):
        pass

    def infer(self, inputs, outputs):

        temperature = inputs["temperature"]
        soil_moisture = inputs["soil_moisture"]
        surface_flux = outputs["surface_flux"]

        surface_flux[...] = 2.0 * temperature + 100.0 * soil_moisture.mean(axis=1)
