#!/usr/bin/env python
'''
This script plots results from MPAS-Ocean convergence test.
'''
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')


# mrp: read from file mesh after nx,ny attributes are added:
nx = 10  # ncfileMesh.getncattr('nx')
ny = 10  # ncfileMesh.getncattr('ny')
iz = 5

tests = ['test1', 'test2', 'test3']
nTests = len(tests)
grids = ['res1', 'res2', 'res3', 'res4']
dx = [10, 20, 50, 100]
nGrids = len(grids)
tracers = ['tracer1', 'tracer2', 'tracer3']
nTracers = len(tracers)

difL2 = np.zeros([nTests, nGrids, nTracers])
errL2 = np.zeros([nTests, nGrids, nTracers])
for i in range(nTests):
    test = tests[i]

    for k in range(nTracers):
        fig = plt.gcf()
        plt.clf()
        fig.set_size_inches(20.0, 20.0)

        for j in range(nGrids):
            grid = grids[j]
            ncfileIC = Dataset('../' + test + '_' + grid + '/init.nc', 'r')
            ncfile = Dataset('../' + test + '_' + grid + '/output.nc', 'r')
            tracer = tracers[k]
            sol = np.reshape(
                ncfileIC.variables[tracer + 'Tend'][0, :, iz], [ny, nx])
            var = np.reshape(
                ncfile.variables[tracer + 'Tend'][0, :, iz], [ny, nx])
            dif = abs(var[1:ny - 1, 1:nx - 2] - sol[1:ny - 1, 1:nx - 2])
            err = abs((var[1:ny - 1, 1:nx - 2] - sol[1:ny - \
                      1, 1:nx - 2]) / sol[1:ny - 1, 1:nx - 2])
            difL2[i, j, k] = np.sqrt(np.mean(dif[:]**2))
            errL2[i, j, k] = np.sqrt(np.mean(err[:]**2))

            # --- Every other row in y needs to average two neighbors in x on planar hex mesh
            var_avg = var
            sol_avg = sol
            for iy in range(0, ny, 2):
                for ix in range(1, nx - 2):
                    var_avg[iy, ix] = (var[iy, ix + 1] + var[iy, ix]) / 2.0
                    sol_avg[iy, ix] = (sol[iy, ix + 1] + sol[iy, ix]) / 2.0

            ax1 = plt.subplot(nGrids, 3, 3 * j + 1)
            plt.imshow(var_avg[1:ny - 1, 1:nx - 2])
            plt.colorbar()
            plt.title(test + ' ' + tracer + ' ' + grid + ' computed')
            
            ax2 = plt.subplot(nGrids, 3, 3 * j + 2)
            plt.imshow(sol_avg[1:ny - 1, 1:nx - 2])
            plt.colorbar()
            plt.title(test + ' ' + tracer + ' ' + grid + ' solution')

            ax3 = plt.subplot(nGrids, 3, 3 * j + 3)
            plt.imshow(var_avg[1:ny - 1, 1:nx - 2] - sol_avg[1:ny - 1, 1:nx - 2])
            plt.colorbar()
            plt.title(test + ' ' + tracer + ' ' + grid + ' error')

            ncfileIC.close()
            ncfile.close()

        plt.savefig(test + '_' + tracer + '_sections.png')


fig = plt.gcf()
plt.clf()
fig.set_size_inches(20.0, 15.0)

for i in range(nTests):
    test = tests[i]
    plt.subplot(nTests, 3, 3 * i + 1)
    for k in range(len(tracers)):
        tracer = tracers[k]
        plt.loglog(dx, difL2[i, :, k], '-x', label=tracer)

    plt.ylabel('diff: rms(exact[:] - calc[:])')
    plt.legend()
    plt.grid()
plt.xlabel('cell width, km')

for i in range(nTests):
    test = tests[i]
    plt.subplot(nTests, 3, 3 * i + 2)
    for k in range(len(tracers)):
        tracer = tracers[k]
        plt.loglog(dx, errL2[i, :, k], '-x', label=tracer)

    plt.title('Error in Redi tendancy term, ' + test)
    plt.ylabel('error: rms((exact[:] - calc[:])/exact[:])')
    plt.legend()
    plt.grid()
plt.xlabel('cell width, km')

for i in range(nTests):
    test = tests[i]
    plt.subplot(nTests, 3, 3 * i + 3)
    for k in range(len(tracers)):
        tracer = tracers[k]
        plt.loglog(dx, difL2[i, :, k] /
                   np.max(difL2[i, :, k]), '-x', label=tracer)

    plt.ylabel('normalized diff')
    plt.grid()
plt.xlabel('cell width, km')

plt.savefig('convergence.png')
