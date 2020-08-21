from winds.wind_model import PROFILE_TYPE
from winds.parameters import Parameters
import math

def test_parameters():
    gridsize = [10, 10]
    nr = 100
    wind_profile_type = PROFILE_TYPE.HOLLAND
    grid_position = [-106.0,35.0]
    cellsize = 2.0
    siderealDay = 23.934  # A sidereal day in hrs.
    omega = 2.0 * math.pi / siderealDay  # The Earth's rotation rate in rad/hr.
    rho = 1.225e9  # Air density at sea level in kg/m^3.
    distance_unit = 'kilometers'
    time_unit = 'hours'
    pressure_unit = 'millibars'
    # The Coriolis parameter should be 2*omega*sin(pi*|phi|/360), for phi in degrees latitude [-90,90].

    params = Parameters(gridsize,nr,wind_profile_type)



def eval_coriolis(lat,omega):
    return 2*omega * math.sin(math.pi*math.fabs(lat)/360)
