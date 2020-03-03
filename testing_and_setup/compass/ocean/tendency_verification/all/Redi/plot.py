#!/usr/bin/env python
'''
This script plots results from MPAS-Ocean convergence test.
'''
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

fig = plt.gcf()
fig.set_size_inches(20.0, 15.0)

# mrp: read from file mesh after nx,ny attributes are added:
nx = 10  # ncfileMesh.getncattr('nx')
ny = 10  # ncfileMesh.getncattr('ny')
iz = 5

tests = ['test1', 'test2', 'test3']
ntests = len(tests)
grids = ['res1', 'res2', 'res3', 'res4']
dx = [10, 20, 50, 100]
ngrids = len(grids)
tracers = ['tracer1', 'tracer2', 'tracer3']

difL2 = np.zeros([ntests, ngrids, 3])
errL2 = np.zeros([ntests, ngrids, 3])
for i in range(ntests):
    test = tests[i]
    for j in range(len(grids)):
        grid = grids[j]
        ncfileIC = Dataset('../' + test + '_' + grid + '/init.nc', 'r')
        ncfile = Dataset('../' + test + '_' + grid + '/output.nc', 'r')
        for k in range(len(tracers)):
            tracer = tracers[k]
            sol = np.reshape(
                ncfileIC.variables[tracer + 'Tend'][0, :, iz], [ny, nx])
            var = np.reshape(
                ncfile.variables[tracer + 'Tend'][0, :, iz], [ny, nx])
            dif = abs(var[1:ny - 1, 1:nx - 1] - sol[1:ny - 1, 1:nx - 1])
            err = abs((var[1:ny - 1, 1:nx - 1] - sol[1:ny - \
                      1, 1:nx - 1]) / sol[1:ny - 1, 1:nx - 1])
            difL2[i, j, k] = np.sqrt(np.mean(dif[:]**2))
            errL2[i, j, k] = np.sqrt(np.mean(err[:]**2))
            #errL2[i,j,k] = np.max(err[:])
        ncfileIC.close()
        ncfile.close()

for i in range(ntests):
    test = tests[i]
    plt.subplot(ntests, 3, 3 * i + 1)
    for k in range(len(tracers)):
        tracer = tracers[k]
        plt.loglog(dx, difL2[i, :, k], '-x', label=tracer)

    plt.ylabel('diff: rms(exact[:] - calc[:])')
    plt.legend()
    plt.grid()
plt.xlabel('cell width, km')

for i in range(ntests):
    test = tests[i]
    plt.subplot(ntests, 3, 3 * i + 2)
    for k in range(len(tracers)):
        tracer = tracers[k]
        plt.loglog(dx, errL2[i, :, k], '-x', label=tracer)

    plt.title('Error in Redi tendancy term, ' + test)
    plt.ylabel('error: rms((exact[:] - calc[:])/exact[:])')
    plt.legend()
    plt.grid()
plt.xlabel('cell width, km')

for i in range(ntests):
    test = tests[i]
    plt.subplot(ntests, 3, 3 * i + 3)
    for k in range(len(tracers)):
        tracer = tracers[k]
        plt.loglog(dx, difL2[i, :, k] /
                   np.max(difL2[i, :, k]), '-x', label=tracer)

    plt.ylabel('normalized diff')
    plt.grid()
plt.xlabel('cell width, km')

plt.savefig('output.png')
