#!/usr/bin/env python

import numpy
import matplotlib.pyplot as plt
import os.path
import netCDF4
from matplotlib import cm

# constants

# Made up (typical for MISOMIP1, for example)
SRef = 34.4

# Constant for the non-local, quadratic melt parameterizaitotn
Kappa = 8.5

# slope of the ice draft
alpha = 1e-2

# thickness of the plume (m)
D = 30.

print "Kappa={}, D={}, alpha={}".format(Kappa, D, alpha)
# Entrainment parameter (Jenkins 1991)
E0 = 0.036

# Stanton number (Jenkins et al. 2010)
St = 5.9e-4


iDepth=-1

# Set depth of grounding line here
#    zGL = -1500.
##    zGL =-688.
##    zGL =-644.
zGL =-693. # value MALI gets for initial condition of Thwaites domain
zCF = -173.  # calving front depth from Thwaites initial condition

#zGL=-1453.0; zCF=-202.0  # TG sim at yr 325
#zGL=-1107.0; zCF=-101.0  # TG sim at yr 225
#zGL=-1235.0; zCF=-125.0  # TG sim approx at both 285 & 430


# use a list to see how the melt rate changes with different profiles.
# use a list with one entry to just plot up a single case.
zLow = -1001.
zUppers =  numpy.linspace(-0.,-600., 13); zLow = zGL-400.  # for variability histogram/profile plot (7)
#zUppers =  numpy.linspace(-0.,-600., 5); zLow = zGL-400.  # For variability range plot in paper (Fig. 3)
#zUppers =  numpy.array([-200., -400.])
#zUppers =  numpy.array([-500.])
#zUppers =  numpy.array([-300.]); zLow = zGL-100.0 # "standard" value for Pine Island Bay

zSill = -663. # estimate used for Thwaites

# idealized Pine Island Bay temperatures
TUpper = -1.
TLower = 1.2


# =============

meanMeltRates = numpy.zeros((len(zUppers),))
GLMeltRates = numpy.zeros((len(zUppers),))

figVar = plt.figure(150, facecolor='w', figsize=(14, 4))
nrow = 1; ncol = 3;

# temperature plot
axTvar = figVar.add_subplot(nrow, ncol, 1)
plt.xlabel('ocean temperature ($^{\circ}$C)')
plt.ylabel('depth (m)')
plt.grid()
axTvar.text(-0.15, 0.95, 'a', transform=axTvar.transAxes, fontsize=14, fontweight='bold')


# melt plot
axMeltvar = figVar.add_subplot(nrow, ncol, 2)
plt.xlabel('melt rate (m yr$^{-1}$)')
plt.ylabel('depth (m)')
plt.grid()
axMeltvar.text(-0.15, 0.95, 'b', transform=axMeltvar.transAxes, fontsize=14, fontweight='bold')

# histogram
axHist = figVar.add_subplot(nrow, ncol, 3)
#plt.xlabel('ice shelf area (km$^2$)')
plt.xlabel('fraction of ice shelf area')
plt.ylabel('depth (m)')
axHist.text(-0.15, 0.95, 'c', transform=axHist.transAxes, fontsize=14, fontweight='bold')

if len(zUppers)>6:
   nOneSide = (len(zUppers)-1)/2
   colorsOneSide = numpy.array( [ cm.cool(x) for x in numpy.linspace(0.0, 1.0, nOneSide) ])
   #print colorsOneSide
else:
   nOneSide = 2
   colorsOneSide = numpy.flip(numpy.array( [ cm.bwr(x) for x in numpy.linspace(0.0, 1.0, nOneSide) ]),0)
   colors = numpy.vstack( (colorsOneSide, [0,0,0,1], numpy.flip(colorsOneSide,0)))
colors = numpy.vstack( (colorsOneSide, [0,0,0,1], numpy.flip(colorsOneSide,0)))

for zUpper in zUppers:
    iDepth = iDepth+1
    print "iDepth={}, zUpper={}".format(iDepth, zUpper)

    z = numpy.linspace(0., zLow, 1001)

    zLower = zUpper - 400.

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

    ind = numpy.where(numpy.logical_and(z>zGL, z<zCF))[0]
    TFMean = numpy.mean(ThermalForcing[ind])

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

    # plot onto temp/melt var plot
    axTvar.plot(TRegional, z, label='$T_{regional}$', color=colors[iDepth,:])
    axMeltvar.plot(meltRate, z, color=colors[iDepth,:])

    if len(zUppers) == 1:
        # If only one plot was requested, also do a 'clean' version for the paper
        #plt.rc('text', usetex=True)
        fig = plt.figure(101, facecolor='w', figsize=(12, 5))
        nrow = 1; ncol = 3;

        # temperature plot
        axT = fig.add_subplot(nrow, ncol, 1)
        plt.xlabel('ocean temperature ($^{\circ}$C)')
        plt.ylabel('depth (m)')
        plt.grid()
        axT.text(-0.15, 0.95, 'a', transform=axT.transAxes, fontsize=14, fontweight='bold')

        plt.plot(TRegional, z, label='$T_{regional}$')
        plt.plot(TCavity, z, label='$T_{cavity}$')
        plt.plot(TFreeze, z, label='$T_{freeze}$')
        plt.plot(TInfinity, z, label='$T_{infinity}$')
        plt.plot(TPlume, z, label='$T_{plume}$')
        plt.plot([TFreeze.min(), TRegional.max()], zGL*numpy.array([1.0, 1.0]), '--k', label='$z_{GL}$')
        plt.plot([TFreeze.min(), TRegional.max()], zCF*numpy.array([1.0, 1.0]), ':k', label='$z_{CF}$')
        plt.plot([TFreeze.min(), TRegional.max()], zSill*numpy.array([1.0, 1.0]), '-.k', label='$z_{sill}$')
        plt.legend()

        # TF plot
        axTF = fig.add_subplot(nrow, ncol, 2)
        plt.xlabel('thermal forcing ($^{\circ}$C)')
        plt.ylabel('depth (m)')
        plt.grid()
        plt.plot(ThermalForcing, z, label='$TF_{local}$')
        plt.plot(TFMean*numpy.array([1.0, 1.0]), [zGL, zCF], label='$TF_{mean}$', linewidth=3)
        plt.plot([0.0, ThermalForcing.max()], zGL*numpy.array([1.0, 1.0]), '--k', label='$z_{GL}$')
        plt.plot([0.0, ThermalForcing.max()], zCF*numpy.array([1.0, 1.0]), ':k', label='$z_{CF}$')
        plt.legend()
        axTF.text(-0.15, 0.95, 'b', transform=axTF.transAxes, fontsize=14, fontweight='bold')


        # melt plot
        axMelt = fig.add_subplot(nrow, ncol, 3)
        plt.xlabel('melt rate (m yr$^{-1}$)')
        plt.ylabel('depth (m)')
        plt.grid()
        axMelt.text(-0.15, 0.95, 'c', transform=axMelt.transAxes, fontsize=14, fontweight='bold')
        # plot Rignot obs if available
        obsFile = '/Users/mhoffman/Documents/PAPERS_PRESENTATIONS/2017_Thwaites_variability/melt_param_testing/iceshelf_melt_param_test/thwaites_1-8km_resolution.cleaned.withRignotMelt.nc'
        if os.path.isfile(obsFile):
            f = netCDF4.Dataset(obsFile, 'r')
            meltObs = f.variables['floatingBasalMassBal'][0,:] / 910.0 * 3600.0 * 24.0 * 365.0
            #ind = numpy.nonzero(meltObs != 0.0)
            lowerSurface = f.variables['lowerSurface'][0,:]
            xCell = f.variables['xCell'][:]
            yCell = f.variables['yCell'][:]
            f.close()
            # divide into east and west shelf regions
            x1=-1590948.400363433; y1=-459735.6052551331;
            x2=-1531877.338559778; y2=-440731.18578141753;
            m = (y2-y1)/(x2-x1); b = y1-m*x1
            ind1 = numpy.nonzero(numpy.logical_and(meltObs!=0.0, yCell>=m*xCell+b))
            ind2 = numpy.nonzero(numpy.logical_and(meltObs!=0.0, yCell< m*xCell+b))

            plt.plot(meltObs[ind1], lowerSurface[ind1], '.', label='$obs_{east}$', markersize = 1)
            plt.plot(meltObs[ind2], lowerSurface[ind2], '.', label='$obs_{west}$', markersize = 1)
        # now plot param.
        plt.plot(meltRate, z, label='model')
        plt.plot([0.0, meltRate.max()], zGL*numpy.array([1.0, 1.0]), '--k', label='$z_{GL}$')
        plt.plot([0.0, meltRate.max()], zCF*numpy.array([1.0, 1.0]), ':k', label='$z_{CF}$')
        plt.legend()
        axMelt.set_ylim(axTF.get_ylim())

        plt.tight_layout()


# clean up temp/melt var plot
xmin=-1.0; xmax=1.2;
axTvar.plot([xmin, xmax], zGL*numpy.array([1.0, 1.0]), '--k', label='$z_{GL}$')
axTvar.plot([xmin, xmax], zCF*numpy.array([1.0, 1.0]), ':k', label='$z_{CF}$')
axTvar.plot([xmin, xmax], zSill*numpy.array([1.0, 1.0]), '-.k', label='$z_{sill}$')

axMeltvar.set_ylim(axTvar.get_ylim())
xmax=72.0
axMeltvar.plot([0.0, xmax], zGL*numpy.array([1.0, 1.0]), '--k', label='$z_{GL}$')
axMeltvar.plot([0.0, xmax], zCF*numpy.array([1.0, 1.0]), ':k', label='$z_{CF}$')
axMeltvar.plot([0.0, xmax], zSill*numpy.array([1.0, 1.0]), '-.k', label='$z_{sill}$')

# plot histogram data
outFile = '/Users/mhoffman/Documents/PAPERS_PRESENTATIONS/2017_Thwaites_variability/2018_OUTPUT/D30_control/output.nc'
if os.path.isfile(outFile):
    f = netCDF4.Dataset(outFile, 'r')
    areaCell = f.variables['areaCell'][:] / 1000.0/1000.0
    lsrf285 = f.variables['lowerSurface'][285,:]
    mask285 = f.variables['cellMask'][285,:]
    ind285 = numpy.nonzero(mask285&4>0)
    bins = numpy.linspace(-1400.0, 0.0, 31)
    axHist.hist(lsrf285[ind285], bins=bins, weights=areaCell[ind285]/areaCell[ind285].sum(), orientation="horizontal", alpha=0.5, label='year 285')

    lsrf430 = f.variables['lowerSurface'][430,:]
    mask430 = f.variables['cellMask'][430,:]
    ind430 = numpy.nonzero(mask430&4>0)
    axHist.hist(lsrf430[ind430], bins=bins, weights=areaCell[ind430]/areaCell[ind430].sum(), orientation="horizontal", alpha=0.5, label='year 430')
    axHist.legend()



figVar.tight_layout()



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
  #plt.plot([0.0, plotMeltMax], [-800, -800],'k-')
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
  #plt.plot([0.0, plotMeltMax], [-800, -800],'k-')
  plt.xlabel('GL melt rate (m/yr)')
  plt.ylabel('depth of thermocline (m)')
  plt.title("kappa={}, D={}, alpha={}, zGL={}, zSill={}".format(Kappa, D, alpha, zGL, zSill))
  plt.grid()

plt.show()
