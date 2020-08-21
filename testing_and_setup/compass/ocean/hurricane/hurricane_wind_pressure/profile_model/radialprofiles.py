import numpy as np
import math

class RadialProfile():

    def __init__(self,n,extent):
        self.profile = np.zeros(n,dtype=np.float64)
        self.rvals = np.zeros(n,dtype=np.float64)
        self.n = n
        self.extent = extent
        self.dr = extent/(n-1)
        for i in range(0,n):
            self.rvals[i] = i*self.dr

    def getValue(self,r):
        if r<0 or r>self.extent:
            return 0.0
        else:
            k = int(r/self.dr)
            return self.rvals[k]

class PressureProfile(RadialProfile):
    def __init__(self,n,extent,pcentral,deltap,rmax):
        super().__init__(n,extent)
        self.pcentral = pcentral
        self.deltap = deltap
        self.rmax = rmax

class HollandPressureProfile(PressureProfile):
    def __init__(self,n,extent,pcentral,deltap,rmax,b):
        super().__init__(n,extent,pcentral,deltap,rmax)
        self.b = b
        for i in range(0,self.n):
            r = self.rvals[i]
            if r>0:
                p = self.pcentral + self.deltap*math.exp(-pow(self.rmax/r,b))
            else:
                p = pcentral
            self.profile[i] = p

class WindSpeedProfile(RadialProfile):
    def __init__(self,n,extent,rmax):
        super().__init__(n,extent)
        self.rmax = rmax
        self.vmax = 0

    def getVmax(self):
        if self.vmax==0:
            for i in range(0,self.n):
                self.vmax = max(self.vmax,self.profile[i])
        return self.vmax

class HollandWindSpeedProfile(WindSpeedProfile):
    def __init__(self,n,extent,rmax,deltap,rho,f,b,coriolis=False):
        super().__init__(n,extent,rmax)
        self.units_factor = 100  # To convert the leading term to m/s
        # This factor comes from adopting millibars instead of Pascals, and km/hr instead of m/s.
        self.deltap = deltap
        self.rho = rho
        self.f = f
        self.b = b
        for i in range(0,self.n):
            r = self.rvals[i]
            if r>0:
                y = pow(rmax/r,b)
                exp_term = self.units_factor*(deltap/rho)*b*y*math.exp(-y)
                if coriolis == True:
                    v = math.sqrt(exp_term + 0.25*pow(r,2)*pow(f,2))+0.5*r*f
                else:
                    v = math.sqrt(exp_term)
            else:
                v = 0.0
            self.profile[i] = v * 3.6 # to convert to km/h
