"""Implementation of the exact solutions for verifying isothermal ice sheet
models. This implementation is based on matlab code by Bueler et al 2005.

Ed Bueler, Craig S. Lingle, Jed A. Kallen-Brown, David N. Covey, and Latrice
N. Bowman, Exact solutions and the verification of numerical models for
isothermal ice sheets, to appear in J. Glaciology"""

class Model:
    """Base class for exact solution ice sheet models."""

    def __init__(self,n,L):
        """Initialise exact solution.

        n: flow law exponent
        L: radius of margin. """

        self.Lkm = L
        self.n   = n

        # some parameters
        self.SperA = 31556926. # seconds per year
        self.Ayr = 1.e-16      # Pre-exponential flow law [Pa^3 a]
        self.rho = 910.        # density of ice [kg/m^3]
        self.g = 9.81          # gravity [m/s^2]
        
        self.calc_constants()    

    def calc_constants(self):
        """Derive some constants from parameters."""

        # flow law parameter
        self.A = self.Ayr/self.SperA
        # overall constant in dormation discharge q_f
        self.Gam = 2.*(self.rho*self.g)**self.n*self.A/(self.n+2)
        # convert to meters
        self.L = self.Lkm*1000.

class ModelAE(Model):

    def __init__(self,M0,n,L):
        """Initialise exact solution.

        M0: accumulation rate in m/a
        n: flow law exponent
        L:  margin radius in km"""

        self.M0yr = M0
        self.bsflags = False

        # some parameters
        self.nut = 8.0e9       # visocity of till [Pa s]

        Model.__init__(self,n,L)

        self.calc_constants()

    def calc_constants(self):
        """Derive some constants from parameters."""

        Model.calc_constants(self)
        self.M0 = self.M0yr/self.SperA
        # constant in H_s(r)
        self.Cs = (2**(self.n-1)*self.M0/self.Gam)**(1./(2*self.n+2))

    def getH(self,r):
        """Get ice thickness.

        r: radius [m]
        return: ice thickness [m]"""

        H = 0.
        if r<self.L:
            H = self.Cs * (self.L**(1+1./self.n) - r**(1.+1./self.n))**(self.n/(2*self.n+2))
        return H

    def getddr(self, r):
        """Get derivatives of ice thickness.

        r: radius [m]
        return: (dHs, ddHs)"""

        if r>=self.L:
            raise ValueError, "radius r=%f is outside region L=%f"%(r,self.L)

        chi = self.L**(1+1./self.n) - r**(1.+1./self.n)
        dHs  = -(self.Cs/2.)*r**(1./self.n)*chi**((self.n-2)/(2*self.n+2))
        ddHs = -(self.Cs/(2*self.n))*(chi**((-3.*self.n-4)/(2*self.n+2))) * \
               (r**((1-self.n)/self.n)*chi + ((self.n+2)/2.)*r**(2/self.n))

        return (dHs, ddHs)

    def getMb(self, r):
        """Get mass balance.

        r: radius [m]
        return: mass balance."""

        mb = 0.
        if r<self.L:
            mb = self.M0

        return mb

class ModelBC(Model):

    def __init__(self,n,lam, H0, R0, tDel):
        """Initialise exact solution.

        n:    flow law exponent
        lam:  accumulation paramter (M_lam = lam t^-1 H)
        H0:   central thickness at t0 (m)
        R0:   radius of magin at t0 (km)
        tDel: time to advance from t0 (tf = t0+tDel) (years)
              if tDel<0 run from t=0 to t0"""
        
        self.lam = lam
        self.H0 = H0
        self.tDel = tDel
        self.R0km = R0
        Model.__init__(self,n,R0)

        self.calc_constants()

    def calc_constants(self):
        """Derive some constants from parameters."""

        Model.calc_constants(self)
        self.R0 = self.R0km*1000.
        self.alf = (2.-(self.n+1)*self.lam)/(5.*self.n+3.)
        self.bet = (1.+(2*self.n+1.)*self.lam)/(5.*self.n+3.)

        # time since creation 
        self.t0 = (self.bet/self.Gam) * ( (2.*n+1.)/(n+1.) )**n * (self.L**(n+1.)/self.H0**(2.*n+1.))
        self.tDel = self.tDel*self.SperA
        if self.tDel < 0:
            self.tf = self.t0
        else:
            self.tf = self.t0 + self.tDel

        # some more constants
        s0 = self.t0**(-self.bet)*self.L
        rmax = self.tf**self.bet * s0
        self.L = 1.1*rmax

    def getH(self,r,t):
        """Get ice thickness.

        r: radius [m]
        t: time (a)
        return: ice thickness [m]"""

        tt = t*self.SperA
        rscl = ( tt**(-self.bet)*r ) / ( self.t0**(-self.bet)*self.R0 )
        temp = max(0., 1-rscl**((self.n+1.)/n))

        H = self.H0*(tt/self.t0)**(-self.alf)*temp**(self.n/(2.*self.n+1.))

    def getMb(self, r, t):
        """Get mass balance.

        r: radius [m]
        t: time (a)
        return: mass balance."""

        tt = t*self.SperA
        if self.tDel<0:
            mb = 0.
        else:
            mb = self.lam/tt*self.getH(r,t)

        return mb
