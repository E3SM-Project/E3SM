#! /usr/bin/env python
# plot basal sediment models

# Copyright (C) 2005, 2010
# Glimmer-CISM contributors - see AUTHORS file for list of contributors
#
# This file is part of Glimmer-CISM.
#
# Glimmer-CISM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or (at
# your option) any later version.
#
# Glimmer-CISM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.
#
# Glimmer-CISM is hosted on BerliOS.de:
# https://developer.berlios.de/projects/glimmer-cism/

import matplotlib.pyplot
import matplotlib
import numpy, math, sys
from  pygsl import integrate

def calc_alpha(Nz,phi):
    return (1/(Nz*math.tan(math.radians(phi))))
def calc_beta(N0,Nz,c,phi):
    return - (N0+c/math.tan(math.radians(phi)))/Nz

def calc_sigma(z,N0,Nz,c,phi):
    """Calculate yield strength.

    z: depth
    N0: effective pressure
    Nz: effective pressure gradient
    c: cohesion
    phi: angle of internal friction"""

    return (N0+Nz*z)*math.tan(math.radians(phi))+c

def calc_za(tau,N0,Nz,c,phi):
    """Calculate depth of deforming layer.

    tau: basal shear stress
    N0: effective pressure
    Nz: effective pressure gradient
    c: cohesion
    phi: angle of internal friction"""

    za = numpy.zeros(len(tau),numpy.float)
    beta  = calc_beta(N0,Nz,c,phi)
    alpha = calc_alpha(Nz,phi)
    
    za = numpy.minimum(0., alpha*tau+beta)
    
    return za

def flow1(z,params):
    """flow law 1.

    tau: basal shear stress
    z: depth
    A: flow law factor
    m: exponent of effective pressure
    n: exponent of basal shear stress
    N0: effective pressure
    Nz: effective pressure gradient"""

    tau = params[0]
    A   = params[1]
    m   = params[2]
    n   = params[3]
    N0  = params[4]
    Nz  = params[5]
    c   = params[6]
    phi = params[7]
    upper  = params[8]

    if za == None:
        return A*tau**n/(N0+Nz*z)**m
    else:
        return (upper-z)*A*tau**n/(N0+Nz*z)**m

def flow2(z,params):
    """flow law 1.

    tau: basal shear stress
    z: depth
    A: flow law factor
    m: exponent of effective pressure
    n: exponent of basal shear stress
    N0: effective pressure
    Nz: effective pressure gradient
    c: cohesion
    phi: angle of internal friction"""

    tau = params[0]
    A   = params[1]
    m   = params[2]
    n   = params[3]
    N0  = params[4]
    Nz  = params[5]
    c   = params[6]
    phi = params[7]
    upper  = params[8]

    if upper == None:
        return A*abs(tau-calc_sigma(z,N0,Nz,c,phi))**n/(N0+Nz*z)**m
    else:
        return (upper-z)*A*abs(tau-calc_sigma(z,N0,Nz,c,phi))**n/(N0+Nz*z)**m
     
def calc_vsld(tau,A,m,n,N0,Nz,c,phi,sigma=False):

    za = calc_za(tau,N0,Nz,c,phi)
    v = numpy.zeros(len(tau),numpy.float)

    if sigma:
        f = flow2
    else:
        f = flow1

    for i in range(0,len(v)):
        sys = integrate.gsl_function(f,[tau[i],A,m,n,N0,Nz,c,phi,None])
        flag,result,error,num = integrate.qng(sys,za[i],0.,1e-8, 1e-8)
        if flag != 0:
            print flag
        v[i] = result
    return v

def calc_vavg(tau,A,m,n,N0,Nz,c,phi,sigma=False):
    """Calculating average sediment velocities.
    
    tau: basal shear stress
    A: flow law factor
    m: exponent of effective pressure
    n: exponent of basal shear stress
    N0: effective pressure
    Nz: effective pressure gradient
    c: cohesion
    phi: angle of internal friction"""

    za = calc_za(tau,N0,Nz,c,phi)
    v = numpy.zeros(len(tau),numpy.float)
    
    if sigma:
        f = flow2
    else:
        f = flow1

    for i in range(0,len(v)):
        sys = integrate.gsl_function(f,[tau[i],A,m,n,N0,Nz,c,phi,0.])
        flag,result,error,num = integrate.qng(sys,za[i],0.,1e-6, 1e-6)
        if flag != 0:
            print flag
        if za[i]!=0:
            v[i] = -result/za[i]
    return v    

if __name__ == '__main__':

    c   = [3.75,15,30,70,70,150,150,250]
    phi = [32,30,27,32,30,32,32,35]
    tau = numpy.arange(0.,100.,1.,numpy.float)
    
    fig = matplotlib.pyplot.figure(figsize=(10,10))

    # create axes
    axprops = {}
    za = fig.add_subplot(3,1,3,**axprops)
    za.set_ylabel("depth of deforming\nsediment layer [m]")
    za.set_xlabel("shear stress [kPa]")
    axprops['sharex'] = za
    
    vsld = fig.add_subplot(3,1,2,**axprops)
    vsld.set_ylabel("sliding velocity [m/a]")

    vavg = fig.add_subplot(3,1,1,**axprops)
    vavg.set_ylabel("average sediment velocity [m/a]")

    def p(tau,N0,Nz,c,phi, colour, name):
        za.plot(tau,calc_za(tau, N0,Nz,c,phi),color=colour,label=name)
        vsld.plot(tau,calc_vsld(tau,107.11,1.35,0.77 ,N0,Nz,c,phi,True),color=colour)
        vsld.plot(tau,calc_vsld(tau, 380.86,2,1,N0,Nz,c,phi,True),color=colour,ls=':')

        vavg.plot(tau,calc_vavg(tau, 34.8,1.8,1.33, N0,Nz,c,phi,False),color=colour)
        vavg.plot(tau,calc_vavg(tau, 380.86,2,1,N0,Nz,c,phi,True),color=colour,ls=':')
        
    
    # Breidamerkurkoekull
    p(tau,50., -10., 3.75, 32., 'red', "Breidamerkurkoekull, N=50kPa")

    # typical till
    p(tau,50., -10., 15, 30, 'green', "typical till, N=50kPa")

    # Breidamerkurkoekull
    p(tau,20., -10., 3.75, 32., 'blue', "Breidamerkurkoekull, N=20kPa")

    # typical till
    p(tau,20., -10., 15, 30, 'cyan', "typical till, N=20kPa")
    
    za.legend(loc=3)

    if len(sys.argv)<2:
        matplotlib.pyplot.show()
    else:
        matplotlib.pyplot.savefig(sys.argv[1])
