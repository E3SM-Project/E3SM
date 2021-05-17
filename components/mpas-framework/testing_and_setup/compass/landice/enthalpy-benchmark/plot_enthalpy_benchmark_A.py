#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import loadmat

anaData = loadmat('./enthA_analy_result.mat')
basalMelt = anaData['basalMelt']

SPY = 31556926
G = 0.042
kc = 2.1
rhow = 1000.0
Lw = 3.34e5

dz = 2.5

data1 = Dataset('./output1.nc','r')
data2 = Dataset('./output2.nc','r')
data3 = Dataset('./output3.nc','r')

yr1 = data1.variables['daysSinceStart'][:]/365.0
yr2 = data2.variables['daysSinceStart'][:]/365.0
yr3 = data3.variables['daysSinceStart'][:]/365.0

basalT1 = data1.variables['basalTemperature'][:,:]
basalT2 = data2.variables['basalTemperature'][:,:]
basalT3 = data3.variables['basalTemperature'][:,:]
basalMeanT1 = np.mean(basalT1, axis=1)
basalMeanT2 = np.mean(basalT2, axis=1)
basalMeanT3 = np.mean(basalT3, axis=1)

basalBmb1 = data1.variables['groundedBasalMassBal'][:,:]
basalBmb2 = data2.variables['groundedBasalMassBal'][:,:]
basalBmb3 = data3.variables['groundedBasalMassBal'][:,:]
basalMeanBmb1 = np.mean(basalBmb1, axis=1)
basalMeanBmb2 = np.mean(basalBmb2, axis=1)
basalMeanBmb3 = np.mean(basalBmb3, axis=1)

basalWaterThickness1 = data1.variables['basalWaterThickness'][:,:]
basalWaterThickness2 = data2.variables['basalWaterThickness'][:,:]
basalWaterThickness3 = data3.variables['basalWaterThickness'][:,:]
basalMeanWaterThickness1 = np.mean(basalWaterThickness1, axis=1)
basalMeanWaterThickness2 = np.mean(basalWaterThickness2, axis=1)
basalMeanWaterThickness3 = np.mean(basalWaterThickness3, axis=1)

T1 = data1.variables['temperature'][:,:,:]
T2 = data2.variables['temperature'][:,:,:]
T3 = data3.variables['temperature'][:,:,:]
TMean1 = np.mean(T1, axis=1)
TMean2 = np.mean(T2, axis=1)
TMean3 = np.mean(T3, axis=1)
TMean1_nvert = TMean1[:,-1]
TMean2_nvert = TMean2[:,-1]
TMean3_nvert = TMean3[:,-1]

basalbmb1 = SPY*(G + kc*(TMean1_nvert - basalMeanT1)/dz)/(rhow*Lw)
basalbmb2 = SPY*(G + kc*(TMean2_nvert - basalMeanT2)/dz)/(rhow*Lw)
basalbmb3 = SPY*(G + kc*(TMean3_nvert - basalMeanT3)/dz)/(rhow*Lw)

Hw1 = np.copy(basalbmb1)
Hw2 = np.copy(basalbmb2)
Hw3 = np.copy(basalbmb3)

for i in range(len(basalbmb1)):
    Hw1[i] = sum(basalbmb1[0:i])*10
for i in range(len(basalbmb2)):
    Hw2[i] = sum(basalbmb2[0:i])*10
for i in range(len(basalbmb3)):
    Hw3[i] = sum(basalbmb3[0:i])*10

year = np.concatenate([yr1[1::], yr2, yr3])/1000.0
basalMeanT = np.concatenate([basalMeanT1[1::],basalMeanT2,basalMeanT3])
basalMeanBmb = np.concatenate([basalMeanBmb1[1::],basalMeanBmb2,basalMeanBmb3])
basalMeanWaterThickness = np.concatenate([basalMeanWaterThickness1[1::],basalMeanWaterThickness2,basalMeanWaterThickness3])
#basalbmb = np.concatenate([basalbmb1[1::],basalbmb2,basalbmb3])
#Hw = np.concatenate([Hw1[1::],Hw2,Hw3])


plt.figure (1)
plt.subplot(311)
plt.plot(year,basalMeanT-273.15)
plt.ylabel('$T_{\\rm b}$ ($^\circ \\rm C$)')
plt.text(10, -28, '(a)', fontsize=20)
plt.grid(True)

plt.subplot(312)
plt.plot(year,-basalMeanBmb*SPY)
plt.plot(basalMelt[1,:]/1000.0, basalMelt[0,:],linewidth=2)
plt.ylabel('$a_{\\rm b}$ (mm a$^{-1}$ w.e.)')
plt.text(10, -1.6, '(b)', fontsize=20)
plt.grid(True)

plt.subplot(313)
plt.plot(year,basalMeanWaterThickness*910.0/1000.0)
plt.ylabel('$H_{\\rm w}$ (m)')
plt.xlabel('Year (ka)')
plt.text(10, 8, '(c)', fontsize=20)
plt.grid(True)


# Create image plot
plotname = 'enthalpy_A_results.png'
plt.savefig(plotname, dpi=150)
print('Saved plot as ' + plotname)

displayImage = False
if displayImage:
   # Note: To get interactive plot, need to comment line at beginning of script "matplotlib.use('Agg')"
   plt.show()
