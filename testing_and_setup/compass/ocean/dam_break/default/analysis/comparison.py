#!/usr/bin/env python
"""

Dam break comparison betewen MPAS-O, data, and ROMS results.

Zhendong Cao and Phillip J. Wolfram
04/05/2019

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import xarray
import pandas as pd
from scipy import spatial
from PIL import Image
from collections import OrderedDict

# render statically by default
plt.switch_backend('agg')

## read output.nc
data = xarray.open_dataset('output.nc')
xCell = data.xCell.values
yCell = data.yCell.values
ssh = data.ssh.values
ssh = ssh+0.6

## read station coordinates
StaData = pd.read_csv('stations/stationCoords.csv',header=None)
StaName = StaData.iloc[:,0].values
StaCoord = StaData.iloc[:,1:]

## find the nearest cell of each station
# coordinate shift from MPAS-O grid to dam break case
xCell = 13 - xCell
yCell = 13 - yCell

matrix = np.array([xCell,yCell])
tree = spatial.KDTree(list(zip(*matrix)))
StaCell = tree.query(StaCoord)[1]

## cells representing station locations
Station = OrderedDict(list(zip(StaName,StaCell)))

ii=0

## plot
fig = plt.figure()
for cell in Station:
    ii+=1
    plt.subplot(3,2,ii+1)

    plt.tight_layout()

  # MPAS-O simulation results
    mpaso=plt.plot(np.linspace(0,100,ssh.shape[0]),ssh[:,Station[cell]],color='#228B22',linewidth=2,alpha=0.6)

  # Measured data
    data = pd.read_csv('stations/'+'Station'+cell+'.csv',header=None)
    measured=plt.scatter(data[0]*10,data[1],4,marker='o',color='k')

  # ROMS simulation results (Warner et al., 2013)
    sim = pd.read_csv('stations/'+cell+'-sim.csv',header=None)
    roms=plt.scatter(sim[0]*10,sim[1],4,marker='v',color='b')

    plt.xlim(0,100)
    plt.xticks(np.arange(0,101,20), np.arange(0,11,2))
    plt.ylim(0,0.7)
    plt.yticks(np.arange(0,0.7,0.2))
    plt.text(35,0.5,'Station '+cell)

    if ii%2==0:
        plt.ylabel('h (m)')
    if ii>=4:
        plt.xlabel('time (s)')

    plt.legend([mpaso[0],measured,roms],['MPAS-O','Measured','ROMS'],
            fontsize='xx-small',frameon=False)

#station location map
im = Image.open('dam_break.png')
im2 = im.resize((650, 300))
plt.subplot(3,2,1)
plt.imshow(im2,interpolation='bicubic')
#plt.axis('off')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.locator_params(axis='x',nbins=3)
plt.xticks([0,325,650],[4,2,0])
plt.locator_params(axis='y',nbins=5)
plt.yticks([0,75,150,220,300],reversed([2,1.5,1,0.5,0]))
#plt.show()
plt.savefig('dam_break_comparison.pdf',dpi=600)
