#!/usr/bin/env python
'''
This script creates an initial condition file for MPAS-Ocean.
'''
import os
import shutil
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset

Ly = 3000.0e3
Lz = 4800.0

def main():
    # {{{

    shutil.copy2('base_mesh.nc', 'initial_state.nc')
    ds = Dataset('initial_state.nc', 'a', format='NETCDF3_64BIT_OFFSET')

    vertical_init(ds)
    tracer_init(ds)
    velocity_init(ds)
    coriolis_init(ds)
    others_init(ds)

    ds.close()
# }}}

def vertical_init(ds):
    thicknessAllLayers = 100.0  # [m] for evenly spaced layers
    nVertLevels = int(Lz / thicknessAllLayers)
    minLayers = 3
# {{{
    # create new variables # {{{
    ds.createDimension('nVertLevels', nVertLevels)
    refLayerThickness = ds.createVariable(
        'refLayerThickness', np.float64, ('nVertLevels',))
    maxLevelCell = ds.createVariable('maxLevelCell', np.int32, ('nCells',))
    refBottomDepth = ds.createVariable(
        'refBottomDepth', np.float64, ('nVertLevels',))
    refZMid = ds.createVariable('refZMid', np.float64, ('nVertLevels',))
    bottomDepth = ds.createVariable('bottomDepth', np.float64, ('nCells',))
    bottomDepthObserved = ds.createVariable(
        'bottomDepthObserved', np.float64, ('nCells',))
    layerThickness = ds.createVariable(
        'layerThickness', np.float64, ('Time', 'nCells', 'nVertLevels',))
    restingThickness = ds.createVariable(
        'restingThickness', np.float64, ('nCells', 'nVertLevels',))
    vertCoordMovementWeights = ds.createVariable(
        'vertCoordMovementWeights', np.float64, ('nVertLevels',))
    # }}}

    # obtain dimensions and mesh variables # {{{
    nCells = len(ds.dimensions['nCells'])
    xCell = ds.variables['xCell']
    yCell = ds.variables['yCell']
    # }}}

    # evenly spaced vertical grid
    refLayerThickness[:] = thicknessAllLayers

    # Create other variables from refLayerThickness
    refBottomDepth[0] = refLayerThickness[0]
    refZMid[0] = -0.5 * refLayerThickness[0]
    for k in range(1, nVertLevels):
        refBottomDepth[k] = refBottomDepth[k - 1] + refLayerThickness[k]
        refZMid[k] = -refBottomDepth[k - 1] - 0.5 * refLayerThickness[k]
    vertCoordMovementWeights[:] = 1.0

    # flat bottom, no bathymetry
    #maxLevelCell[:] = nVertLevels
    #bottomDepth[:] = refBottomDepth[nVertLevels-1]
    #bottomDepthObserved[:] = refBottomDepth[nVertLevels-1]
    # for k in range(nVertLevels):
    #    layerThickness[0,:,k] = refLayerThickness[k]
    #    restingThickness[:,k] = refLayerThickness[k]

    # Define bottom depth: parabola
    for iCell in range(0, nCells):
        x = xCell[iCell]
        y = yCell[iCell]
        bottomDepthObserved[iCell] \
            = 1.1 * Lz * (1.0 - ((y - Ly / 2) / (Ly / 2))**2) - 100

    # full cells, not partial
    # initialize to very bottom:
    maxLevelCell[:] = nVertLevels
    bottomDepth[:] = refBottomDepth[nVertLevels - 1]
    for k in range(nVertLevels):
        layerThickness[0, :, k] = refLayerThickness[k]
        restingThickness[:, k] = refLayerThickness[k]
    for iCell in range(0, nCells):
        x = xCell[iCell]
        y = yCell[iCell]
        for k in range(nVertLevels):
            if bottomDepthObserved[iCell] < refBottomDepth[k]:
                maxLevelCell[iCell] = max(k, minLayers)
                bottomDepth[iCell] = refBottomDepth[maxLevelCell[iCell] - 1]
                break
# }}}

def tracer_init(ds):
    slope = 0.001
    # temperature: linear, match slope
    Tmin = 5.0
    Tx = 0.0
    Ty = 15.0 / Ly
    Tz = Ty / slope

    # salinity: linear, match slope
    Smin = 15.0
    Sx = 0.0
    Sy = 15.0 / Ly
    Sz = Sy / slope

    # tracer1: Gaussian
    y0 = Ly / 2
    yr = Ly / 4
    z0 = -Lz / 3
    zr = Lz / 4
# {{{

    # create new variables # {{{
    tracer1 = ds.createVariable(
        'tracer1', np.float64, ('Time', 'nCells', 'nVertLevels',))
    tracer2 = ds.createVariable(
        'tracer2', np.float64, ('Time', 'nCells', 'nVertLevels',))
    tracer3 = ds.createVariable(
        'tracer3', np.float64, ('Time', 'nCells', 'nVertLevels',))
    temperature = ds.createVariable(
        'temperature', np.float64, ('Time', 'nCells', 'nVertLevels',))
    salinity = ds.createVariable(
        'salinity', np.float64, ('Time', 'nCells', 'nVertLevels',))
    layerThickness = ds.variables['layerThickness']
    # }}}

    # obtain dimensions and mesh variables # {{{
    nVertLevels = len(ds.dimensions['nVertLevels'])
    nCells = len(ds.dimensions['nCells'])
    xCell = ds.variables['xCell']
    yCell = ds.variables['yCell']
    refZMid = ds.variables['refZMid']
    refBottomDepth = ds.variables['refBottomDepth']
    # }}}
    for iCell in range(0, nCells):
        x = xCell[iCell]
        y = yCell[iCell]
        for k in range(0, nVertLevels):
            z = refZMid[k]

            temperature[0, iCell, k] \
                = Tx * x + Ty * y + Tz * z
            salinity[0, iCell, k] \
                = Sx * x + Sy * y + Sz * z
            tracer1[0, iCell, k] \
                = 1.0 + np.exp(
                    -((y - y0) / yr)**2
                    - ((z - z0) / zr)**2)
            tracer2[0, iCell, k] = 1.0
            tracer3[0, iCell, k] = 1.0

        tracer2[0, iCell, 10:20] = int(2 + np.cos(y * 4 * 2 * np.pi / Ly))

        if ((y > Ly / 4) & (y < Ly / 2)) | (y > 3 * Ly / 4):
            tracer3[0, iCell, 0:20] = 2.0

    # Normalize T&S:
    temperature[:] += Tmin - np.min(temperature[:])
    salinity[:] += Smin - np.min(salinity[:])
    print(
        'Temperature ranges from ', np.min(
            temperature[:]), ' to ', np.max(
            temperature[:]))
    print(
        'Salinity ranges from ', np.min(
            salinity[:]), ' to ', np.max(
            salinity[:]))
# }}}

def velocity_init(ds):
    # {{{
    normalVelocity = ds.createVariable(
        'normalVelocity', np.float64, ('Time', 'nEdges', 'nVertLevels',))
    normalVelocity[:] = 0.0
# }}}

def coriolis_init(ds):
    # {{{
    fEdge = ds.createVariable('fEdge', np.float64, ('nEdges',))
    fEdge[:] = 0.0
    fVertex = ds.createVariable('fVertex', np.float64, ('nVertices',))
    fVertex[:] = 0.0
    fCell = ds.createVariable('fCell', np.float64, ('nCells',))
    fCell[:] = 0.0
# }}}

def others_init(ds):
    # {{{
    surfaceStress = ds.createVariable(
        'surfaceStress', np.float64, ('Time', 'nEdges',))
    surfaceStress[:] = 0.0
    atmosphericPressure = ds.createVariable(
        'atmosphericPressure', np.float64, ('Time', 'nCells',))
    atmosphericPressure[:] = 0.0
    boundaryLayerDepth = ds.createVariable(
        'boundaryLayerDepth', np.float64, ('Time', 'nCells',))
    boundaryLayerDepth[:] = 0.0
# }}}

if __name__ == '__main__':
    # If called as a primary module, run main
    main()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
