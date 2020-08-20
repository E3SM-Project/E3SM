from enum import Enum
import numpy as np
import winds.parameters as Parameters
import hurricane_model as Hurricane
import structures as Geogrid
import matplotlib.pyplot as plt
import math


class PROFILE_TYPE(Enum):
    HOLLAND = 'holland'
    WILLOUGHBY = 'willoughby'

class WindModel:
    def __init__(self, params: Parameters, curr_hurricane: Hurricane, grid: Geogrid):
        self.profile_type = params.wind_profile_type
        if self.profile_type == PROFILE_TYPE.HOLLAND:

            # Distance between the hurricane eye and the grid points
            # Great circle distance in km
            r = np.power(np.sin((grid.lat - curr_hurricane.center[1]) * 0.5), 2) + \
                np.cos(grid.lat) * np.cos(curr_hurricane.center[1]) * \
                np.power(np.sin((grid.lon - curr_hurricane.center[0]) * 0.5), 2)
            r = 2.0 * params.earth_radius * np.arcsin(np.sqrt(r))

            # Compute pressure
            self.pressure_profile = holland_pressure_profile(curr_hurricane, r)

            # Compute wind speed
            self.wind_speed_profile = holland_windspeed_profile(params, curr_hurricane, r)

            # plt.scatter(grid.lon, grid.lat, s=10., c=self.wind_speed_profile, alpha=1.)
            # plt.show()

            # Compute wind components
            self.u, self.v = compute_components(self.wind_speed_profile, curr_hurricane, grid)

        else:
            raise 'Profile models other than Holland are not currently supported.'


def holland_pressure_profile(hurricane: Hurricane, r: np.ndarray):
        """
        :param hurricane: class type Hurricane
        :param r: distance between the eye of the hurricane and the grid points in km
        """
        return hurricane.pcentral + hurricane.deltap * np.exp(-np.power(hurricane.extent / r ,hurricane.b))


def holland_windspeed_profile(params: Parameters, hurricane: Hurricane, r: np.ndarray, coriolis=False):
        """
        :param params: class parameters
        :param hurricane: class Hurricane
        :param r: distance between the eye of the hurricane and the grid points in km
        :param coriolis: coriolis factor in rad/hrs
        """

        # Holland equation assumes:
        # deltap in Pa
        # density in kg/m3
        # and returns m/s
        units_factor = 100.  # To convert the deltap from mbar to Pascals


        y = np.power(hurricane.extent / r, hurricane.b)
        exp_term = units_factor*(hurricane.deltap / params.rho) * hurricane.b * y * np.exp(-y)
        if coriolis is True:
            v = np.sqrt(exp_term + 0.25 * np.power(r * params.f, 2)) + 0.5 * r * params.f
        else:
            v = np.sqrt(exp_term)

        v *= 3.6 # Conversion from m/s to km/h

        return v

def compute_components(wind_speed_profile, curr_hurricane: Hurricane, grid: Geogrid) -> (np.ndarray, np.ndarray):
    # Compute components of vg
    theta = np.arctan2(grid.lat - curr_hurricane.center[1], grid.lon - curr_hurricane.center[0])
    theta += math.pi * 0.5
    vg_x = wind_speed_profile * np.cos(theta)
    vg_y = wind_speed_profile * np.sin(theta)

    # Compute total velocity
    ratio = wind_speed_profile / curr_hurricane.vmax
    u = vg_x + curr_hurricane.vforward[0] * ratio
    v = vg_y + curr_hurricane.vforward[1] * ratio

    return u, v