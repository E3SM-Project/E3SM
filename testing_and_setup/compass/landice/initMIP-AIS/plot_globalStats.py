from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import copy

data = Dataset('./ctrl_forcing/globalStats_ctrl.nc')

targetVars = ['totalFloatingBasalMassBal',
                    'volumeAboveFloatation',
                    'groundingLineFlux',
                    'groundedIceArea',
                    'totalCalvingFlux',
                    'floatingIceArea']
# variables that we want to see

targetVars_scaleFunc = [lambda x: x*(-1)*1.0e-12, 
                        lambda x: x/x[0]*100, 
                        lambda x: x*1.0e-12, 
                        lambda x: x/x[0]*100, 
                        lambda x: x*1.0e-12, 
                        lambda x: x/x[0]*100]
# scaling for each variable; change yLabels for different scalings

yLabels = ['Basal Melt (Gt yr$^{-1}$)', 'VAF (% change)', 'GL flux (Gt yr$^{-1}$)', 'Grounded area (% change)', 'calving flux (Gt yr$^{-1}$)', 'Floating area (% change)']
xLabels = 'Year'

grpColorsLines = '-k'
deltaT = 1.

subplot_col = 2
subplot_row = 3

subplot_num = subplot_row*subplot_col

varNum = len(targetVars)

fig = plt.figure(1)

for i in range(varNum):

    fig.add_subplot(subplot_row,subplot_col, i+1)
    targetVarData = data[targetVars[i]][1::]
    dataArray = np.array(targetVarData)

    yearNum = len(dataArray)*deltaT
    annualMeanData = np.mean(np.reshape(dataArray,(int(yearNum),int(1.0/deltaT))),axis=1)

    scaleFunc = targetVars_scaleFunc[i]
    annualMeanDataScaled = scaleFunc(annualMeanData)

    years = np.array(range(int(yearNum)))+1
    #plt.plot(targetGrpRegionVarData)
    plt.plot(years,annualMeanDataScaled,grpColorsLines)
    plt.ylabel(yLabels[i])
    plt.xlim(0,yearNum)

    if i == varNum-1 or i == varNum-2:
        plt.xlabel(xLabels)

    #if i == 0:
    #    plt.legend((grpLegendNames))

plt.show()
