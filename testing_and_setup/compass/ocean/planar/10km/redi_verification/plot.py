import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np

fig = plt.gcf()
fig.set_size_inches(8.0,10.0)
nRow=1 #2
nCol=3
iTime=[0,1]

ncfileMesh = Dataset('../base_mesh/planar_hex_mesh.nc','r')
nx=ncfileMesh.getncattr('nx')
ny=ncfileMesh.getncattr('ny')
ncfileMesh.close()

ncfile = Dataset('output.nc','r')
ncfileIC = Dataset('../initial_state/initial_state.nc','r')
#var = ncfile.variables['temperature']
var = ncfile.variables['tracer1']
varSim = ncfile.variables['tracer1Tend']
term1 = ncfileIC.variables['term1']
term2 = ncfileIC.variables['term2']
term3 = ncfileIC.variables['term3']
varExact = term1[:] + term2[:] + term3[:]
xtime = ncfile.variables['xtime']
layer = [0,2,4] # vertical layers
#print(var.shape)

iCol=1
plt.subplot(nRow, nCol, 1)
tmp = np.reshape(varSim[1,:,layer[iCol]],[ny,nx])
print(np.shape(tmp))
ax = plt.imshow(tmp[:,1:nx-1])
plt.title('sim:tend')
#ax.set_clim(-8,8)

plt.subplot(nRow, nCol, 2)
ax = plt.imshow(np.reshape(varExact[0,:,layer[iCol]],[ny,nx])[:,1:nx-1])
plt.title('exact')
#ax.set_clim(-8,8)

plt.subplot(nRow, nCol, 3)
#ax = plt.imshow(np.reshape(var[0,:,layer[iCol]],[ny,nx]))
#plt.title('tracer1')
tmp = np.reshape(varSim[1,:,layer[iCol]],[ny,nx])
diff = tmp[1:ny-1,1:nx-1]
ax = plt.imshow(diff)
plt.title('diff')

for iCol in range(nCol):
    for iRow in range(nRow):
        plt.subplot(nRow, nCol, iRow*nCol+iCol+1) 
#        ax = plt.imshow(np.reshape(varSim[1,:,layer[iCol]],[ny,nx])) #,extent=[0,160,0,500])
        plt.axis('off')
#        #diff = varSim[1,:,layer[iCol]] - 600.0*varEx[0,:,layer[iCol]]
#        #print(max(abs(diff)))
#        #ax = plt.imshow(np.reshape(diff,[ny,nx])) #,extent=[0,160,0,500])
        plt.jet()
#        #plt.clim([11.8,13])
#        #if iRow==nRow-1:
#        #    plt.xlabel('x, km')
#        #if iCol==0:
#        #    plt.ylabel('y, km')
        plt.colorbar()
#        #print(xtime[iTime[iCol],11:13])
#        plt.title('title')

ncfile.close()
ncfileIC.close()
#plt.savefig('redi_output.png')
plt.savefig('term2.png')
