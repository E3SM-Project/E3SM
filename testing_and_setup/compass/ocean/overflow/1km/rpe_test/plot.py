import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy

fig = plt.gcf()
fig.set_size_inches(8.0,10.0)
nRow=1 #6
nCol=2
nu=['0.01','0.1','1','10','100','1000']
iTime=[3,6]
time=['3 hrs','6 hrs']
for iRow in range(nRow):
    ncfile = Dataset('output_'+str(iRow+1)+'.nc','r')
    var = ncfile.variables['temperature']
    xtime = ncfile.variables['xtime']
    for iCol in range(nCol):
        plt.subplot(nRow, nCol, iRow*nCol+iCol+1) 
        ax = plt.imshow(var[iTime[iCol],0::4,:].T,extent=[0,200,2000,0],aspect=2)
        plt.clim([10,20])
        plt.jet()
        if iRow==nRow-1:
            plt.xlabel('x, km')
        if iCol==0:
            plt.ylabel('depth, m')
        plt.colorbar()
        #print(xtime[iTime[iCol],11:13])
        plt.title('time='+time[iCol]+', nu='+nu[iRow])
    ncfile.close()
plt.savefig('sections_overflow.png')
