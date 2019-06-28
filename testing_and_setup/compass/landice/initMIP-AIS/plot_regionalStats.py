#!/usr/bin/env python
'''
Script to plot regionalStats figure used in MALI description paper.
'''

from __future__ import absolute_import, division, print_function, unicode_literals

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import copy

##########################################################################
## Change settings here !!
##########################################################################

#regionalData = Dataset('./ctrl_forcing/regionalStats_ctrl/regionalStats.nc')
#regionalData = Dataset('BMB/regionalStats.nc')
regionalData = Dataset('regionalStats.BMB.annual.nc')
regionalDataC = Dataset('regionalStats.control.annual.nc')
#regionalData = Dataset('regionalStats.BMB.nc')
#regionalDataC = Dataset('regionalStats.control.nc')
#regionalDataC = Dataset('/scratch1/scratchdirs/hoffman2/ais-2km/initMIP/BMB/continue_after_flood_fill_yr65/regionalStats.nc')

targetRegions = [["FR_West","FR_East","FR_Central"],
                 ["Amery_West","Amery_East","Amery_Central"],
                 ["Ross_EAIS","Ross_WAIS_N","Ross_WAIS_S"],
               ["Totten"],
                 ["Thwaites","Pine_Island"]]
# regions that we want to see; refer to regionImbieDic below; 
# the code will plot the sum of regions in each sub-list which is thought as a combined drainage basin

targetRegionVars = ['regionalSumFloatingBasalMassBal',
                    'regionalVolumeAboveFloatation',
                    'regionalSumGroundingLineFlux',
                    'regionalGroundedIceArea']
#                    'regionalSumCalvingFlux',
#                    'regionalFloatingIceArea']
# variables that we want to see

targetRegionVars_scaleFunc = [lambda x: x*(-1)*1.0e-12, 
                              #lambda x: x/x[0]*100-100, 
                              lambda x: (x-x[0])*910.0/1.0e12, 
                              lambda x: x*1.0e-12, 
                              lambda x: (x-x[0])/1.0e6,   #x/x[0]*100-100, 
                              lambda x: x*1.0e-12, 
                              lambda x: x/x[0]*100-100]
# scaling for each variable; change yLabels for different scalings

yLabels = ['Ice shelf basal melt (Gt yr$^{-1}$)', 
   #'VAF (% change)', 
   'Volume above\nfloatation change (Gt)', 
   'Grounding line flux (Gt yr$^{-1}$)', 
   'Grounded area change (km$^2$)', 
   'calving flux (Gt yr$^{-1}$)', 
   'Floating area (% change)']
xLabels = 'Year'

grpLegendNames = ['Filchner-Ronne','Lambert-Amery','Ross','Thwaites-PIG']
grpLegendNames = ['Filchner-Ronne','Lambert-Amery','Ross','Totten','Thwaites-Pine Island']
#grpColorsLines = ['-k','-m','-b','-r']
#grpColorsLines = ['--k','--m','--b','--r']

regionImbieDic = {"FR_West":0, "FR_Central":11, "FR_East":20,
                  "Amery_West":26, "Amery_Central":1, "Amery_East":2,
                  "Totten":4,
                  "Ross_WAIS_S":9, "Ross_WAIS_N":10, "Ross_EAIS":8,
                  "Thwaites":13,
                  "Pine_Island":14,
                  "Totten":4}
# a complete dictionary for Antarctica regions and their IMBIE numbers
# TODO need further completion

deltaT = 1.
# output timestep, unit: year; for example, deltaT = 1/12. if timestep is one month

subplot_col = 2
#subplot_row = 3
subplot_row = 2

###########################################################################



subplot_num = subplot_row*subplot_col


regionGrpNum = len(targetRegions)
regionNumEachGrp = copy.copy(targetRegions)

for i in range(regionGrpNum):
    regionNumEachGrp[i] = len(targetRegions[i])

regionNum = sum(regionNumEachGrp)


regionVarNum = len(targetRegionVars)
assert(subplot_num >= regionVarNum)


regionId = []
for i in range(regionGrpNum):
    for j in range(regionNumEachGrp[i]):
        regionId.append(regionImbieDic[targetRegions[i][j]])


fig = plt.figure(1)
cmap = plt.get_cmap("tab10")

for i in range(regionVarNum):

    ax = fig.add_subplot(subplot_row,subplot_col, i+1)

    targetRegionVarData = regionalData[targetRegionVars[i]][:,regionId]
    dataArray = np.array(targetRegionVarData)
    years = regionalData['daysSinceStart'][:] / 365.0

    targetRegionVarDataC = regionalDataC[targetRegionVars[i]][:,regionId]
    dataArrayC = np.array(targetRegionVarDataC)
    yearsC = regionalDataC['daysSinceStart'][:] / 365.0

    for j in range(regionGrpNum):
        #targetGrpRegionVarData = np.zeros((len(years),))
        #targetGrpRegionVarData = dataArray[:,currentGrpLen]
        if j == 0:
            currentGrpLen = regionNumEachGrp[j]
            targetGrpRegionVarData = np.sum(dataArray[:,0:currentGrpLen],axis=1)

            targetGrpRegionVarDataC = np.sum(dataArrayC[:,0:currentGrpLen],axis=1)
        else:
            previousGrpLen = sum(regionNumEachGrp[0:j])
            currentGrpLen = regionNumEachGrp[j]
            targetGrpRegionVarData = np.sum(dataArray[:,previousGrpLen:previousGrpLen+currentGrpLen],axis=1)
            targetGrpRegionVarDataC = np.sum(dataArrayC[:,previousGrpLen:previousGrpLen+currentGrpLen],axis=1)

        yearNum = len(targetGrpRegionVarData)*deltaT
        annualMeanData = np.mean(np.reshape(targetGrpRegionVarData,(int(yearNum),int(1.0/deltaT))),axis=1)

        scaleFunc = targetRegionVars_scaleFunc[i]
        #annualMeanDataScaled = scaleFunc(annualMeanData)
        annualMeanDataScaled = scaleFunc(targetGrpRegionVarData)
        annualMeanDataScaledC = scaleFunc(targetGrpRegionVarDataC)

        #plt.plot(years,targetGrpRegionVarData,grpColorsLines[j])
        plt.plot(years,annualMeanDataScaled,color=cmap(j), linewidth=2, label=grpLegendNames[j])
        plt.plot(yearsC,annualMeanDataScaledC,color=cmap(j),linewidth=0.5, label=None)

    plt.ylabel(yLabels[i])
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
#        plt.xlim(0,yearNum)

    if i == regionVarNum-1 or i == regionVarNum-2:
        plt.xlabel(xLabels)

    if i == 3:
        #plt.legend((grpLegendNames))
        plt.legend()
plt.tight_layout()
plt.show()
