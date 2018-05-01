from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import copy

##########################################################################
## Change settings here !!
##########################################################################

regionalData = Dataset('./ctrl_forcing/regionalStats_ctrl/regionalStats.nc')
#regionalData = Dataset('./bmb_forcing/regionalStats_bmb/regionalStats.nc')

targetRegions = [["FR_West","FR_East","FR_Central"],
                 ["Amery_West","Amery_East","Amery_Central"],
                 ["Ross_EAIS","Ross_WAIS_N","Ross_WAIS_S"],
                 ["Thwaites","Pine_Island"]]
# regions that we want to see; refer to regionImbieDic below; 
# the code will plot the sum of regions in each sub-list which is thought as a combined drainage basin

targetRegionVars = ['regionalSumFloatingBasalMassBal',
                    'regionalVolumeAboveFloatation',
                    'regionalSumGroundingLineFlux',
                    'regionalGroundedIceArea',
                    'regionalSumCalvingFlux',
                    'regionalFloatingIceArea']
# variables that we want to see

targetRegionVars_scaleFunc = [lambda x: x*(-1)*1.0e-12, 
                              lambda x: x/x[0]*100, 
                              lambda x: x*1.0e-12, 
                              lambda x: x/x[0]*100, 
                              lambda x: x*1.0e-12, 
                              lambda x: x/x[0]*100]
# scaling for each variable; change yLabels for different scalings

yLabels = ['Basal Melt (Gt yr$^{-1}$)', 'VAF (% change)', 'GL flux (Gt yr$^{-1}$)', 'Grounded area (% change)', 'calving flux (Gt yr$^{-1}$)', 'Floating area (% change)']
xLabels = 'Year'

grpLegendNames = ['Filchner-Ronne','Lambert-Amery','Ross','Thwaites-PIG']
grpColorsLines = ['-k','-m','-b','-r']
#grpColorsLines = ['--k','--m','--b','--r']

regionImbieDic = {"FR_West":0, "FR_Central":11, "FR_East":20,
                  "Amery_West":26, "Amery_Central":1, "Amery_East":2,
                  "Totten":4,
                  "Ross_WAIS_S":9, "Ross_WAIS_N":10, "Ross_EAIS":8,
                  "Thwaites":13,
                  "Pine_Island":14}
# a complete dictionary for Antarctica regions and their IMBIE numbers
# TODO need further completion

deltaT = 1.
# output timestep, unit: year; for example, deltaT = 1/12. if timestep is one month

subplot_col = 2
subplot_row = 3

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

for i in range(regionVarNum):

    fig.add_subplot(subplot_row,subplot_col, i+1)

    targetRegionVarData = regionalData[targetRegionVars[i]][1::,regionId]
    dataArray = np.array(targetRegionVarData)

    for j in range(regionGrpNum):
        if j == 0:
            currentGrpLen = regionNumEachGrp[j]
            targetGrpRegionVarData = np.sum(dataArray[:,0:currentGrpLen],axis=1)
        else:
            previousGrpLen = sum(regionNumEachGrp[0:j])
            currentGrpLen = regionNumEachGrp[j]
            targetGrpRegionVarData = np.sum(dataArray[:,previousGrpLen:previousGrpLen+currentGrpLen],axis=1)

        yearNum = len(targetGrpRegionVarData)*deltaT
        annualMeanData = np.mean(np.reshape(targetGrpRegionVarData,(int(yearNum),int(1.0/deltaT))),axis=1)

        scaleFunc = targetRegionVars_scaleFunc[i]
        annualMeanDataScaled = scaleFunc(annualMeanData)

        years = np.array(range(int(yearNum)))+1
        plt.plot(years,annualMeanDataScaled,grpColorsLines[j])
        plt.ylabel(yLabels[i])
        plt.xlim(0,yearNum)

    if i == regionVarNum-1 or i == regionVarNum-2:
        plt.xlabel(xLabels)

    if i == 0:
        plt.legend((grpLegendNames))

plt.show()
