#!/usr/bin/env python
'''
Plots profiles for hydro-ship test case
'''
import sys
import numpy as np
import netCDF4
#import datetime
# import math
# from pylab import *
from optparse import OptionParser
import matplotlib.pyplot as plt
from matplotlib import cm
# from matplotlib.contour import QuadContourSet
# import time

secInYr = 3600.0 * 24.0 * 365.0  # Note: this may be slightly wrong for some calendar types!

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="file to visualize", metavar="FILE")
parser.add_option("-t", "--time", dest="time", help="time step to visualize (0 based)", metavar="TIME")
parser.add_option("-s", "--save", action="store_true", dest="saveimages", help="include this flag to save plots as files")
parser.add_option("-n", "--nodisp", action="store_true", dest="hidefigs", help="include this flag to not display plots (usually used with -s)")
parser.add_option("-3", dest="A3", action="store_true", help="plot GLADS results for experiment 3")
parser.add_option("-5", dest="A5", action="store_true", help="plot GLADS results for experiment 5")

options, args = parser.parse_args()

if not options.filename:
	print "No filename provided. Using output.nc."
        options.filename = "output.nc"

if not options.time:
	print "No time provided. Using time -1."
        time_slice = -1
else:
        time_slice = int(options.time)


if options.A3 and options.A5:
   sys.exit("Only one of -3 and -5 can be specified.")

f = netCDF4.Dataset(options.filename,'r')

#xtime = f.variables['xtime'][:]
xCell = f.variables['xCell'][:]
yCell = f.variables['yCell'][:]
yVertex = f.variables['yVertex'][:]
xEdge = f.variables['xEdge'][:]
#yEdge = f.variables['yEdge'][:]
h = f.variables['waterThickness'][time_slice,:]
u = f.variables['waterVelocityCellX'][time_slice,:]
N = f.variables['effectivePressure'][time_slice,:]
days = f.variables['daysSinceStart'][:]
#basalmelt = f.variables['basalMeltInput'][time_slice,:]
#surfmelt = f.variables['externalWaterInput'][time_slice,:]
xtime= f.variables['xtime']
areaCell = f.variables['areaCell'][:]

q = u*h


#print "attempting to get input data from landice_grid.nc!"
#fin = netCDF4.Dataset('landice_grid.nc','r')
#H = fin.variables['thickness'][0,:]

print "Using time level ", time_slice, ", which is xtime=",''.join( xtime[time_slice,:])

# Find center row  - currently files are set up to have central row at y=0
unique_ys=np.unique(yCell[:])
centerY=unique_ys[len(unique_ys)/2]
print "number of ys, center y index, center Y value", len(unique_ys), len(unique_ys)/2, centerY
ind = np.nonzero(yCell[:] == centerY)
x = xCell[ind]/1000.0

# calculate mean,min,max for all x values for needed variables
allx=np.unique(xCell[:])
N_mean = np.zeros(allx.shape)
N_min = np.zeros(allx.shape)
N_max = np.zeros(allx.shape)
q_mean = np.zeros(allx.shape)
q_min = np.zeros(allx.shape)
q_max = np.zeros(allx.shape)
for i in range(len(allx)):
    N_mean[i] = N[ xCell == allx[i] ].mean()
    N_min[i] = N[ xCell == allx[i] ].min()
    N_max[i] = N[ xCell == allx[i] ].max()
    q_mean[i] = q[ xCell == allx[i] ].mean()
    q_min[i] = q[ xCell == allx[i] ].min()
    q_max[i] = q[ xCell == allx[i] ].max()

print "start plotting."

fig = plt.figure(1, facecolor='w')
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312, sharex=ax1)
ax3 = fig.add_subplot(313, sharex=ax1)

if options.A3:
   #x,N_mean,N_min,N_max,q_mean,q_min,q_max,Q_max
   data = np.array([
   [1250,3.6938e+05,8918,8.3148e+05,0.00057512,0.00054189,0.00060403,0.0005687],
   [3750,1.1954e+06,8.6108e+05,1.5055e+06,0.00055975,0.00051241,0.00061171,1.385e-05],
   [6250,1.7328e+06,1.5285e+06,1.8977e+06,0.00054136,0.00049619,0.00058063,4.908e-05],
   [8750,1.9756e+06,1.9151e+06,2.0184e+06,0.00052717,0.00051715,0.00054088,6.225e-05],
   [11250,2.0292e+06,2.0191e+06,2.0342e+06,0.00051337,0.00050593,0.00052117,6.208e-05],
   [13750,2.0343e+06,2.0337e+06,2.0349e+06,0.00049898,0.00049205,0.00050688,5.7544e-05],
   [16250,2.0333e+06,2.0328e+06,2.034e+06,0.00048492,0.00047704,0.00049275,5.0174e-05],
   [18750,2.0334e+06,2.0328e+06,2.0341e+06,0.00047067,0.00046278,0.00047807,4.3253e-05],
   [21250,2.0356e+06,2.0342e+06,2.0376e+06,0.00045614,0.00044857,0.00046326,3.7132e-05],
   [23750,2.0401e+06,2.0376e+06,2.0431e+06,0.00044145,0.00043417,0.00044871,3.3286e-05],
   [26250,2.0466e+06,2.0432e+06,2.0505e+06,0.000427,0.00041967,0.00043462,2.8756e-05],
   [28750,2.055e+06,2.0508e+06,2.0597e+06,0.00041219,0.00040515,0.00041988,2.4631e-05],
   [31250,2.0648e+06,2.0598e+06,2.0705e+06,0.00039822,0.0003912,0.00040547,2.2199e-05],
   [33750,2.0769e+06,2.0707e+06,2.083e+06,0.00038385,0.00037639,0.00039104,1.9196e-05],
   [36250,2.0898e+06,2.0832e+06,2.0968e+06,0.0003691,0.00036202,0.00037614,1.6603e-05],
   [38750,2.1048e+06,2.0973e+06,2.1128e+06,0.00035445,0.00034747,0.00036198,1.4817e-05],
   [41250,2.1219e+06,2.1131e+06,2.1303e+06,0.00034021,0.00033275,0.00034758,1.263e-05],
   [43750,2.1398e+06,2.1305e+06,2.1489e+06,0.00032575,0.00031846,0.00033287,1.0957e-05],
   [46250,2.1591e+06,2.1499e+06,2.1697e+06,0.00031096,0.0003043,0.00031857,9.5511e-06],
   [48750,2.1805e+06,2.1704e+06,2.1925e+06,0.00029685,0.00028949,0.00030424,8.2628e-06],
   [51250,2.2057e+06,2.193e+06,2.2175e+06,0.00028231,0.00027498,0.00028954,7.0915e-06],
   [53750,2.2314e+06,2.2186e+06,2.2442e+06,0.0002678,0.0002608,0.00027504,5.953e-06],
   [56250,2.259e+06,2.2444e+06,2.2727e+06,0.00025315,0.00024601,0.00026062,4.9826e-06],
   [58750,2.2887e+06,2.2735e+06,2.3033e+06,0.00023867,0.00023181,0.00024572,4.185e-06],
   [61250,2.3212e+06,2.3056e+06,2.338e+06,0.0002241,0.00021699,0.00023169,3.5131e-06],
   [63750,2.3571e+06,2.3385e+06,2.3749e+06,0.00020952,0.00020255,0.00021702,2.9005e-06],
   [66250,2.3949e+06,2.3759e+06,2.415e+06,0.0001953,0.00018817,0.00020254,2.3651e-06],
   [68750,2.438e+06,2.4168e+06,2.4592e+06,0.00018099,0.00017366,0.00018823,1.8838e-06],
   [71250,2.4837e+06,2.4621e+06,2.5086e+06,0.00016647,0.00015922,0.00017381,1.4542e-06],
   [73750,2.5348e+06,2.51e+06,2.562e+06,0.0001516,0.00014474,0.0001594,1.1529e-06],
   [76250,2.5908e+06,2.5622e+06,2.6215e+06,0.0001375,0.00013076,0.00014454,8.6251e-07],
   [78750,2.6556e+06,2.6219e+06,2.6876e+06,0.00012298,0.0001159,0.00013037,6.2556e-07],
   [81250,2.7256e+06,2.6895e+06,2.7603e+06,0.00010867,0.0001013,0.00011553,4.1966e-07],
   [83750,2.8049e+06,2.7631e+06,2.8462e+06,9.4164e-05,8.6859e-05,0.00010129,2.9364e-07],
   [86250,2.8947e+06,2.8512e+06,2.9432e+06,7.9323e-05,7.2377e-05,8.7046e-05,1.7608e-07],
   [88750,2.9973e+06,2.9436e+06,3.0541e+06,6.4644e-05,5.7607e-05,7.2224e-05,9.9334e-08],
   [91250,3.1202e+06,3.0561e+06,3.184e+06,5.0362e-05,4.3453e-05,5.765e-05,4.8309e-08],
   [93750,3.2618e+06,3.1878e+06,3.3358e+06,3.6105e-05,2.8436e-05,4.3565e-05,1.5298e-08],
   [96250,3.4222e+06,3.3441e+06,3.5114e+06,2.2029e-05,1.4291e-05,2.9408e-05,3.7029e-09],
   [98750,3.6054e+06,3.5227e+06,3.6783e+06,7.1355e-06,1.5613e-06,1.5243e-05,2.3962e-10],])

if options.A5:
   #x,N_mean,N_min,N_max,q_mean,q_min,q_max,Q_max
   data = np.array([
   [1250,3.619e+05,8918,7.9712e+05,0.0010021,0.00017955,0.001585,39.744],
   [3750,1.0909e+06,8.1947e+05,1.3215e+06,0.0013428,0.0012158,0.0015387,38.601],
   [6250,1.4278e+06,1.3165e+06,1.5309e+06,0.0012151,0.00097281,0.0014001,38.361],
   [8750,1.53e+06,1.4328e+06,1.6546e+06,0.001147,0.00089554,0.0013547,37.825],
   [11250,1.5416e+06,1.429e+06,1.6859e+06,0.0011332,0.0008849,0.0013364,36.779],
   [13750,1.5158e+06,1.4051e+06,1.6878e+06,0.0011489,0.00091201,0.0013329,35.373],
   [16250,1.4954e+06,1.384e+06,1.6386e+06,0.0011667,0.00093384,0.0013389,33.595],
   [18750,1.4606e+06,1.3643e+06,1.5949e+06,0.0011937,0.00094267,0.0013531,31.617],
   [21250,1.4328e+06,1.3424e+06,1.5689e+06,0.0012208,0.001004,0.0013729,29.621],
   [23750,1.4017e+06,1.3192e+06,1.5327e+06,0.0012504,0.0010469,0.0014016,27.568],
   [26250,1.3721e+06,1.2872e+06,1.5025e+06,0.0012757,0.0010689,0.0014265,25.957],
   [28750,1.351e+06,1.2621e+06,1.4725e+06,0.0013092,0.0011136,0.0014584,23.621],
   [31250,1.3192e+06,1.2396e+06,1.4348e+06,0.0013414,0.0011444,0.0014804,21.843],
   [33750,1.2893e+06,1.2179e+06,1.406e+06,0.0013761,0.0011585,0.0014989,20.081],
   [36250,1.2651e+06,1.2029e+06,1.3636e+06,0.0014143,0.0012063,0.0015224,17.839],
   [38750,1.232e+06,1.1845e+06,1.3467e+06,0.00148,0.0012409,0.0017218,15.268],
   [41250,1.2015e+06,1.1576e+06,1.2575e+06,0.0014802,0.0013333,0.0015753,13.513],
   [43750,1.1874e+06,1.1379e+06,1.2476e+06,0.0014928,0.0013516,0.0016143,12.013],
   [46250,1.169e+06,1.1104e+06,1.2367e+06,0.0015339,0.0013754,0.001671,10.481],
   [48750,1.1426e+06,1.0807e+06,1.2262e+06,0.0016014,0.0014019,0.0017374,8.3152],
   [51250,1.1058e+06,1.0392e+06,1.1898e+06,0.0016687,0.001452,0.0018008,6.1913],
   [53750,1.0677e+06,1.0072e+06,1.1592e+06,0.0017757,0.0015356,0.0020891,4.4044],
   [56250,1.014e+06,9.6914e+05,1.1225e+06,0.001886,0.0016167,0.0022956,2.2276],
   [58750,9.7705e+05,9.617e+05,1.035e+06,0.0018537,0.0017904,0.0021797,0.0036676],
   [61250,9.748e+05,9.63e+05,9.8723e+05,0.0017411,0.0016804,0.0018283,0.0015423],
   [63750,9.9444e+05,9.7833e+05,1.0079e+06,0.0016276,0.0015713,0.001692,0.00091574],
   [66250,1.0255e+06,1.0078e+06,1.0438e+06,0.0015172,0.0014609,0.0015759,0.00062522],
   [68750,1.0656e+06,1.0448e+06,1.0858e+06,0.0014061,0.001349,0.0014638,0.00042488],
   [71250,1.11e+06,1.0888e+06,1.1342e+06,0.0012935,0.0012371,0.0013503,0.00028731],
   [73750,1.1597e+06,1.1356e+06,1.186e+06,0.001178,0.0011249,0.0012386,0.00021198],
   [76250,1.2134e+06,1.186e+06,1.2423e+06,0.0010685,0.0010165,0.0011229,0.00015782],
   [78750,1.274e+06,1.2427e+06,1.3035e+06,0.00095566,0.00090092,0.0010127,0.00011363],
   [81250,1.3379e+06,1.3052e+06,1.3691e+06,0.00084453,0.00078734,0.00089811,7.546e-05],
   [83750,1.4086e+06,1.3716e+06,1.445e+06,0.00073179,0.00067526,0.00078665,5.2317e-05],
   [86250,1.4873e+06,1.4494e+06,1.5294e+06,0.00061648,0.00056304,0.00067593,3.1047e-05],
   [88750,1.5761e+06,1.5297e+06,1.6251e+06,0.00050237,0.0004486,0.0005615,1.7476e-05],
   [91250,1.6824e+06,1.6269e+06,1.7379e+06,0.0003914,0.00033855,0.00044755,8.5425e-06],
   [93750,1.8065e+06,1.7412e+06,1.8724e+06,0.00028062,0.00022185,0.00033801,2.7783e-06],
   [96250,1.9517e+06,1.88e+06,2.0349e+06,0.0001711,0.00011134,0.00022772,7.0808e-07],
   [98750,2.126e+06,2.0454e+06,2.1977e+06,5.5348e-05,1.2068e-05,0.00011773,5.022e-08], ])


if options.A3 or options.A5:
   G_x = data[:,0]
   G_Nmean = data[:,1]
   G_Nmin = data[:,2]
   G_Nmax = data[:,3]
   G_qmean = data[:,4]
   G_qmin = data[:,5]
   G_qmax = data[:,6]
   G_Qmax = data[:,7]

   # plot GLADS data
   lw = 3  # lineweight to use
   ax1.plot(G_x/1000.0, G_Nmin/1.0e6, 'g--', linewidth=lw)
   ax1.plot(G_x/1000.0, G_Nmean/1.0e6, 'g-', linewidth=lw, label='GLADS mean/range')
   ax1.plot(G_x/1000.0, G_Nmax/1.0e6, 'g--', linewidth=lw)

   ax2.plot(G_x/1000.0, G_qmin, 'g--', linewidth=lw)
   ax2.plot(G_x/1000.0, G_qmean, 'g-', linewidth=lw, label='GLADS mean/range')
   ax2.plot(G_x/1000.0, G_qmax, 'g--', linewidth=lw)

   ax3.plot(G_x/1000.0, G_Qmax, 'g--', linewidth=lw, label='GLADS max')

# panel 1: effective pressure
plt.sca(ax1)
#plt.plot(x, N[ind] / 1.0e6, '.-g')  # this just plots the centerline profile
plt.plot(allx/1000.0, N_mean / 1.0e6, '-b', label='MPAS mean/range')
plt.plot(allx/1000.0, N_min / 1.0e6, '--b')
plt.plot(allx/1000.0, N_max / 1.0e6, '--b')

plt.xlabel('X-position (km)')
plt.ylabel('effecive pressure (MPa)')
plt.xlim( (0, 100.0) )
plt.grid(True)
plt.legend(loc='best')

# panel 2: sheet flux
plt.sca(ax2)
#plt.plot(x, np.absolute(h[ind] * u[ind]), '.-g') # this plots centerline profile
plt.plot(allx/1000.0, np.absolute(q_mean), '-b', label='MPAS mean/range')
plt.plot(allx/1000.0, np.absolute(q_min), '--b')
plt.plot(allx/1000.0, np.absolute(q_max), '--b')

plt.xlabel('X-position (km)')
plt.ylabel('sheet water flux (m^2/s)')
plt.grid(True)
plt.legend(loc='best')

# panel 3: channel flux
plt.sca(ax3)
try:
   channelDischarge = f.variables['channelDischarge'][time_slice,:]
   allxEdge=np.unique(xEdge[:])
   allxEdge2=100000.0 - (allxEdge - 100000.0)
   Q_max = np.zeros(allxEdge.shape)
   Q_sum = np.zeros(allxEdge.shape)
   for i in range(len(allxEdge)):
     Q_max[i] = np.absolute(channelDischarge[ xEdge == allxEdge[i] ]).max()
     Q_sum[i] = np.absolute(channelDischarge[ xEdge == allxEdge[i] ]).sum()

   plt.plot(allxEdge/1000.0, np.absolute(Q_max), 'bo', label='MPAS max')
   plt.plot(allxEdge/1000.0, np.absolute(Q_sum), 'bx', label='MPAS sum')

except:
   print "Skipping plotting of channel output."

plt.xlabel('X-position (km)')
plt.ylabel('channel water flux (m^3/s)')
plt.grid(True)
plt.legend(loc='best')





# plot how close to SS we are
fig = plt.figure(2, facecolor='w')
ax1 = fig.add_subplot(211)
for i in ind:
    plt.plot(days/365.0, f.variables['waterThickness'][:,i])
plt.xlabel('Years since start')
plt.ylabel('water thickness (m)')
plt.grid(True)

ax = fig.add_subplot(212, sharex=ax1)
for i in ind:
    plt.plot(days/365.0, f.variables['effectivePressure'][:,i]/1.0e6)
plt.xlabel('Years since start')
plt.ylabel('effective pressure (MPa)')
plt.grid(True)


# plot time steps for various
try:
   dtA=f.variables['deltatSGHadvec'][:]
   dtD=f.variables['deltatSGHdiffu'][:]
   dtP=f.variables['deltatSGHpressure'][:]
   dtC=f.variables['deltatSGHchannel'][:]
   fig = plt.figure(3, facecolor='w')
   plt.plot(days/365.0, dtA, label='A')
   plt.plot(days/365.0, dtD, label='D')
   plt.plot(days/365.0, dtP, label='P')
   plt.plot(days/365.0, dtC, label='C')
   plt.legend()
   plt.xlabel('Time (yr)')
   plt.ylabel('Allowable time step (s)')
except:
   pass

print "plotting complete"

plt.draw()
if options.saveimages:
        print "Saving figures to files."
        plt.savefig('GL-position.png')




if options.hidefigs:
     print "Plot display disabled with -n argument."
else:
     plt.show()

