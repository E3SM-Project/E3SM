import datetime

class Hurricane:
    def __init__(self, center: tuple, extent: float, pcentral: float, deltap: float,
                 vmax: float, b: float, time: float, initial_datetime: datetime.datetime):
        self.center = center            # Position of the eye (lon,lat) in radians as tuple.
        self.extent = extent            # The maximum extent of the hurricane in kilometers.
        self.vforward = []              # Forward velocity [ve, vn] in km/hr.
        self.pcentral = pcentral        # Central pressure in millibars.
        self.deltap = deltap            # Pressure difference in millibars.
        self.vmax = vmax                # The maximum gradient wind [ve, vn] in km/hr.
        self.b = b                      # The Holland parameter, conventionally in the range [0.5,2.5]
        self.time = time                # Time of this trajectory point in hours.
        self.ref_time = initial_datetime


    def set_vf(self, vf: tuple):
        self.vforward = vf
