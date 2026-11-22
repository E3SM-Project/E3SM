import numpy as np


MASK_VAL = -9.99
class LandEmulator:
    def __init__(self, config):
        pass

    def infer(self, inputs, outputs):

        temperature = inputs["temperature"]
        soil_moisture = inputs["soil_moisture"]
        filter_column = inputs["filter_column"]
        surface_flux = outputs["surface_flux"]
        columns = filter_column[filter_column >= 0]

        surface_flux[:] = MASK_VAL

        surface_flux[columns] = 2.0 * temperature + 100.0 * soil_moisture.mean(axis=1)
