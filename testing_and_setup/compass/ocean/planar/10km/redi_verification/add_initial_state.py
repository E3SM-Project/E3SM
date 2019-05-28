#!/usr/bin/env python
'''
This script creates an initial condition file for MPAS-Ocean.
'''
import os
import shutil
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset

def main(): 
# {{{

    shutil.copy2('base_mesh.nc','initial_state.nc')
    ds = Dataset('initial_state.nc', 'a', format='NETCDF3_64BIT_OFFSET')

    vertical_init(ds)
    tracer_init(ds)
    velocity_init(ds)
    coriolis_init(ds)
    others_init(ds)

    ds.close()
# }}}

def vertical_init(ds):
# {{{

    # config settings
    nVertLevels = 10
    thicknessAllLayers = 2 # [m] for evenly spaced layers

    # create new variables # {{{
    ds.createDimension('nVertLevels', nVertLevels)
    refLayerThickness = ds.createVariable('refLayerThickness', np.float64, ('nVertLevels',))
    maxLevelCell = ds.createVariable('maxLevelCell', np.int32, ('nCells',))
    refBottomDepth = ds.createVariable('refBottomDepth', np.float64, ('nVertLevels',))
    refZMid = ds.createVariable('refZMid', np.float64, ('nVertLevels',))
    bottomDepth = ds.createVariable('bottomDepth', np.float64, ('nCells',))
    bottomDepthObserved = ds.createVariable('bottomDepthObserved', np.float64, ('nCells',))
    layerThickness = ds.createVariable('layerThickness', np.float64, ('Time','nCells','nVertLevels',))
    restingThickness = ds.createVariable('restingThickness', np.float64, ('nCells','nVertLevels',))
    vertCoordMovementWeights = ds.createVariable('vertCoordMovementWeights', np.float64, ('nVertLevels',))
    # }}}

    # evenly spaced vertical grid
    refLayerThickness[:] = thicknessAllLayers

    # Create other variables from refLayerThickness
    refBottomDepth[0] = refLayerThickness[0]
    refZMid[0] = -0.5*refLayerThickness[0]
    for k in range(1,nVertLevels):
        refBottomDepth[k] = refBottomDepth[k-1] + refLayerThickness[k]
        refZMid[k] = -refBottomDepth[k-1]-0.5*refLayerThickness[k]
    vertCoordMovementWeights[:] = 1.0

    # flat bottom, no bathymetry
    maxLevelCell[:] = nVertLevels
    bottomDepth[:] = refBottomDepth[nVertLevels-1]
    bottomDepthObserved[:] = refBottomDepth[nVertLevels-1]
    for k in range(nVertLevels):
        layerThickness[0,:,k] = refLayerThickness[k]
        restingThickness[:,k] = refLayerThickness[k]
# }}}
# define the passive debug tracer
def PDef(x,y,z):
    return P0*(x**px * z**pz)
    #return P0*(y**px * z**pz)

# define the passive debug tracer
def TDef(x,y,z):
    return Tx*x**pTx + Tz*z**pTz
    #return Tx*y**pTx + Tz*z**pTz

def tracer_init(ds):
# {{{

    px=2;py=0;pz=2;pTx=2;pTy=0;pTz=2;
    # config settings
    h = 2 # thickness of all layers
    Kappa = 600.0 # kappa for Redi
    # tracer1
    P0 = 10000

    # temperature
    T0 = 20
    Tx = 5/1000
    Ty = 0
    Tz = 10.0/1000

    S0 = 35

    # create new variables # {{{
    tracer1 = ds.createVariable('tracer1', np.float64, ('Time','nCells','nVertLevels',))
    term1 = ds.createVariable('term1', np.float64, ('Time','nCells','nVertLevels',))
    term2 = ds.createVariable('term2', np.float64, ('Time','nCells','nVertLevels',))
    term3 = ds.createVariable('term3', np.float64, ('Time','nCells','nVertLevels',))
    temperature = ds.createVariable('temperature', np.float64, ('Time','nCells','nVertLevels',))
    salinity = ds.createVariable('salinity', np.float64, ('Time','nCells','nVertLevels',))
    layerThickness = ds.variables['layerThickness']
    # }}}

    # obtain dimensions and mesh variables # {{{
    nVertLevels = len(ds.dimensions['nVertLevels'])
    nCells = len(ds.dimensions['nCells'])
    xCell = ds.variables['xCell']
    yCell = ds.variables['yCell']
    # For periodic domains, the max cell coordinate is also the domain width
    Lx = max(xCell) 
    Ly = max(yCell)
    refZMid = ds.variables['refZMid']
    refBottomDepth = ds.variables['refBottomDepth']
    H = max(refBottomDepth)
    # }}}

    for iCell in range(0,nCells):
        x = xCell[iCell]
        y = yCell[iCell]
        for k in range(0,nVertLevels):
            z = refZMid[k]

            tracer1[0,iCell,k] \
                = P0*(x**px * y**py * z**pz)
                #= P0*(y**px * z**pz)
            temperature[0,iCell,k] \
                = Tx*x**pTx + Ty*y**pTy + Tz*z**pTz
                #= Tx*y**pTx  + Tz*z**pTz

            salinity[0,iCell,k] = S0
            #2*P0*z     6*P0*Tx*x**2/Tz      4*P0*Tx*x**2/Tz
            term1[0,iCell,k] = h*Kappa* \
                0
            term2[0,iCell,k] = h*Kappa* \
                6*P0*Tx/Tz
            term3[0,iCell,k] = h*Kappa* \
                2*P0*Tx/Tz
    print('term1  min/max: %12.5E %12.5E'%(np.min(term1[0,:,:]),np.max(term1[0,:,:])))
    print('term2  min/max: %12.5E %12.5E'%(np.min(term2[0,:,:]),np.max(term2[0,:,:])))
    print('term3  min/max: %12.5E %12.5E'%(np.min(term3[0,:,:]),np.max(term3[0,:,:])))
    print('sum123 min/max: %12.5E %12.5E'%(np.min(term1[0,:,:]+term2[0,:,:]+term3[0,:,:]),np.max(term1[0,:,:]+term2[0,:,:]+term3[0,:,:])))
# }}}

def velocity_init(ds):
# {{{
    normalVelocity = ds.createVariable('normalVelocity', np.float64, ('Time','nEdges','nVertLevels',))
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
    surfaceStress = ds.createVariable('surfaceStress', np.float64, ('Time','nEdges',))
    surfaceStress[:] = 0.0
    atmosphericPressure = ds.createVariable('atmosphericPressure', np.float64, ('Time','nCells',))
    atmosphericPressure[:] = 0.0
    boundaryLayerDepth = ds.createVariable('boundaryLayerDepth', np.float64, ('Time','nCells',))
    boundaryLayerDepth[:] = 0.0
# }}}

if __name__ == '__main__':
    # If called as a primary module, run main
    main()

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
