import math
from winds.wind_model import PROFILE_TYPE

class Parameters:
    def __init__(self, mean_lat: float, wind_profile_type=PROFILE_TYPE.HOLLAND):
        """
        Constructor.
        :param mean_lat: mean latitude of the hurricane trajectory to compute the Coroilis factor in radians
        Units are km, hr, and millibars for distance, wind, and pressure respectively, and lat in decimal degrees.
        """
        self.siderealDay = 23.934 # A sidereal day in hrs.
        self.omega = 2.0 * math.pi / self.siderealDay # The Earth's rotation rate in rad/hr.

        self.rho = 1.15 # Air density at sea level in kg/m^3.
        self.wind_profile_type = wind_profile_type # The particular wind profile model being used.

        # Earth radius in km
        self.earth_radius = 6371.1


    def get_coriolis(self, lat: float) -> float:
        """
        Returns the Coriolis parameter for a given latitude.
        :param lat: in radians
        :return: coriolis factor in rad/s to be consistent with Holland's model units
        """
        # The Coriolis parameter = 2*omega*sin(|phi|)
        # 3600 to convert omega in rad/s
        return 2.0 * self.omega / 3600. * math.sin(math.fabs(lat))


    def get_pressure_unit(self):
        return 'millibars'


    def get_distance_unit(self):
        return 'kilometers'


    def get_time_unit(self):
        return 'hours'
