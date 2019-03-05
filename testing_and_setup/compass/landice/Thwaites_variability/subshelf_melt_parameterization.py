#!/usr/bin/env python

import numpy
import matplotlib.pyplot as plt

# constants

# Made up (typical for MISOMIP1, for example)
SRef = 34.4

# Constant for the non-local, quadratic melt parameterizaitotn
Kappa = 8.5

# slope of the ice draft
alpha = 1e-2

# thickness of the plume (m)
D = 30.


D=15.0
Kappa=11.0

print "Kappa={}, D={}, alpha={}".format(Kappa, D, alpha)
# Entrainment parameter (Jenkins 1991)
E0 = 0.036

# Stanton number (Jenkins et al. 2010)
St = 5.9e-4


iDepth=-1

zUppers =  numpy.linspace(-0.,-600., 25)
#zUppers =  numpy.array([-200., -400.])
#zUppers =  numpy.array([-300.])
#zUppers =  numpy.array([-500.])
meanMeltRates = numpy.zeros((len(zUppers),))
GLMeltRates = numpy.zeros((len(zUppers),))
for zUpper in zUppers:
    iDepth = iDepth+1
    print "iDepth={}, zUpper={}".format(iDepth, zUpper)

    zGL = -1500.
##    zGL =-688.
#    zGL =-662.
##    zGL =-644.

    z = numpy.linspace(0., zGL-100, 1001)

    zLower = zUpper - 400.

    zSill = -663.

    TUpper = -1.
    TLower = 1.2

    TRegional = numpy.zeros(z.shape)

    TRegional[z >= zUpper] = TUpper
    TRegional[z <= zLower] = TLower
    mask = numpy.logical_and(z < zUpper, z > zLower)
    TRegional[mask] = \
        (TUpper - TLower)*(z[mask] - zLower)/(zUpper - zLower) + TLower

    TSill = numpy.interp(zSill, z[::-1], TRegional[::-1])

    TCavity = TRegional.copy()
    TCavity[z <= zSill] = TSill

    # Jenkins (1991)
    TFreeze = 0.0901 - 0.0575*SRef + 7.61e-4*z

    # length scale of decay of plume temperature
    zeta = alpha*D/(E0*alpha + St)

    # limit of exponental decay of T
    TInfinity = (E0*alpha*TCavity + St*TFreeze)/(E0*alpha + St)

    TGL = numpy.interp(zGL, z[::-1], TCavity[::-1])

    exponent = numpy.minimum(-(z-zGL)/zeta, 0.)

    TPlume = TInfinity + (TGL - TInfinity)*numpy.exp(exponent)

    plt.figure(iDepth, figsize=(9,4))
    plt.subplot(1,2,1)
    plt.title("kappa={}, D={}, alpha={}, zGL={}".format(Kappa, D, alpha, zGL))
    plt.plot(TRegional, z, label='TRegional')
    plt.plot(TCavity, z, label='TCavity')
    plt.plot(TFreeze, z, label='TFreeze')
    plt.plot(TInfinity, z, label='TInfinity')
    plt.plot(TPlume, z, label='TPlume')
    plt.legend()

    ThermalForcing = TPlume - TFreeze

    TFMean = numpy.mean(ThermalForcing)

    print "TFmean=", TFMean

    meltRate = Kappa*TFMean*ThermalForcing

    plt.subplot(1,2,2)
    plt.plot(meltRate, z, label='meltRate')
    plt.legend()
    plt.grid()

    #plt.draw()

    print "melt_GL={}, melt_200m={}".format(meltRate[numpy.nonzero(z<=zGL)[0][0]], meltRate[numpy.nonzero(z<=-200.0)[0][0]])
    meanMeltRates[iDepth] = meltRate[z>zGL].mean()
    print "integrated melt rate=", meanMeltRates[iDepth]

    GLMeltRates[iDepth] = meltRate[numpy.argmin(numpy.absolute((z-zGL))).min()]
if len(zUppers) > 1:
  # plot summarizing everything for mean melt rate
  plt.figure(200)
  plt.plot(meanMeltRates, zUppers, '.')
  plt.plot(meanMeltRates, zUppers-400.0, '.')
  plotMeltMax=meanMeltRates.max()
  plt.plot([0.0, plotMeltMax], [zSill, zSill],'r--')
  plt.plot([0.0, plotMeltMax], [-700, -700],'k-')
  plt.plot([0.0, plotMeltMax], [-400, -400],'k:')
  plt.plot([0.0, plotMeltMax], [-1000, -1000],'k:')
#  plt.plot([0.0, plotMeltMax], [-800, -800],'k-')
  plt.xlabel('mean melt rate (m/yr)')
  plt.ylabel('depth of thermocline (m)')
  plt.title("kappa={}, D={}, alpha={}, zGL={}, zSill={}".format(Kappa, D, alpha, zGL, zSill))
  plt.grid()

  # repeat for melt at GL
  plt.figure(201)
  plt.plot(GLMeltRates, zUppers, '.')
  plt.plot(GLMeltRates, zUppers-400.0, '.')
  plotMeltMax=GLMeltRates.max()
  plt.plot([0.0, plotMeltMax], [zSill, zSill],'r--')
  plt.plot([0.0, plotMeltMax], [-700, -700],'k-')
  plt.plot([0.0, plotMeltMax], [-400, -400],'k:')
  plt.plot([0.0, plotMeltMax], [-1000, -1000],'k:')
#  plt.plot([0.0, plotMeltMax], [-800, -800],'k-')
  plt.xlabel('GL melt rate (m/yr)')
  plt.ylabel('depth of thermocline (m)')
  plt.title("kappa={}, D={}, alpha={}, zGL={}, zSill={}".format(Kappa, D, alpha, zGL, zSill))
  plt.grid()

plt.show()
