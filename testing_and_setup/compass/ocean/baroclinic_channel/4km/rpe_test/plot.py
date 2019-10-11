import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

fig = plt.gcf()
fig.set_size_inches(8.0,10.0)
nRow=1 #2
nCol=5
nu=['1','5','10','100','200']
iTime=[0,1]
time=['0 days','20 days']
for iCol in range(nCol):
    ncfile = Dataset('output_'+str(iCol+1)+'.nc','r')
    #ncfile = Dataset('initial_state.nc','r')
    var = ncfile.variables['temperature']
    #xtime = ncfile.variables['xtime']
    for iRow in range(nRow):
        print(var)
        plt.subplot(nRow, nCol, iRow*nCol+iCol+1) 
        ax = plt.imshow(np.reshape(var[iTime[iRow],:,0],[126,40]),extent=[0,160,0,500])
        plt.gca().invert_yaxis()
        plt.clim([11.8,13])
        plt.jet()
        if iRow==nRow-1:
            plt.xlabel('x, km')
        if iCol==0:
            plt.ylabel('y, km')
        plt.colorbar()
        #print(xtime[iTime[iCol],11:13])
        plt.title('time='+time[iRow]+', nu='+nu[iCol])
    ncfile.close()
plt.savefig('sections_baroclinic_channel_4km.png')
