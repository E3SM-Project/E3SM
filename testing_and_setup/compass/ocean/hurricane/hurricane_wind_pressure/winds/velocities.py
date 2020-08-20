import math


class Velocities:

    def __init__(self, vfe, vfn, vmax):
        """
        Initialize with the forward velocity components.
        :param vfe: Eastward forward velocity (x-component in the Earth frame) in km/hr.
        :param vfn: Northward forward velocity component (y-component in the Earth frame) in km/hr.
        """
        self.vf = []
        self.vfmagn = []
        self.xunitv = []
        self.yunitv = []

        self.set_vforward(vfe, vfn)
        self.vmax = vmax

    def set_vforward(self, vfe, vfn):
        self.vf = [vfe, vfn]
        self.vfmagn = math.sqrt(pow(vfe, 2) + pow(vfn, 2))
        self.xunitv = [vfn/self.vfmagn, -vfe/self.vfmagn]
        self.yunitv = [vfe/self.vfmagn, vfn/self.vfmagn]


    def compute_wind_vector(self, vg, xe, yn):
        """
        Returns the velocity components [ve,vn] given the tangential gradient wind speed.
        :param vg: The tangential (theta) gradient wind speed in the hurricane frame in km/hr.
        :param xe: The eastern component of position relative to the local origin (the hurricane eye) in km.
        :param yn: The northern component of position relative to the local origin (the hurricane eye) in km.
        :return: [ve,vn] the eastward and nortward components of the wind velocity in the Earth frame in km/hr.
        """
        rmagn = math.sqrt(xe*xe + yn*yn)

        costheta = (xe*self.xunitv[0] + yn*self.xunitv[1])/rmagn
        sintheta = -(xe*self.xunitv[1] - yn*self.xunitv[0])/rmagn
        theta_unitv = [-sintheta*self.xunitv[0]+costheta*self.yunitv[0],
                       -sintheta*self.xunitv[1]+costheta*self.yunitv[1]]
        vgtheta = [theta_unitv[0]*vg, theta_unitv[1]*vg]
        vfcorr = vg/self.vmax
        ve = self.vf[0]*vfcorr + vgtheta[0]
        vn = self.vf[1]*vfcorr + vgtheta[1]

        return [ve, vn]

