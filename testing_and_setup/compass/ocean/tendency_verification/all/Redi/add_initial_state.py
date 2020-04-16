#!/usr/bin/env python
'''
This script creates an initial condition file for MPAS-Ocean.
'''
# import packages {{{
import os
import shutil
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import argparse
# }}}

# output from polynomial_validation.py:
# {{{
#    Temperature T                  tracer P             Slope Sx                   term 1     term 2                                                      term 3                         powers
# 1   Tx*(x + x0) + Ty + Tz*z        P0*z*(x + x0)        -Tx/Tz                     0          -P0*Tx/Tz                                                   -P0*Tx/Tz                      px=1;py=0;pz=1;pTx=1;pTy=0;pTz=1;
# 2   Tx*(x + x0) + Ty + Tz*z        P0*z*(x + x0)**2     -Tx/Tz                     2*P0*z     -P0*Tx*(2*x + 2*x0)/Tz                                      -P0*Tx*(2*x + 2*x0)/Tz         px=2;py=0;pz=1;pTx=1;pTy=0;pTz=1;
# 3   Tx*(x + x0) + Ty + Tz*z        P0*z**2*(x + x0)     -Tx/Tz                     0          -2*P0*Tx*z/Tz                                               -2*P0*Tx*z/Tz                  px=1;py=0;pz=2;pTx=1;pTy=0;pTz=1;
# 4   Tx*(x + x0) + Ty + Tz*z        P0*z**2*(x + x0)**2  -Tx/Tz                     2*P0*z**2  -2*P0*Tx*z*(2*x + 2*x0)/Tz                                  -2*P0*Tx*z*(2*x + 2*x0)/Tz     px=2;py=0;pz=2;pTx=1;pTy=0;pTz=1;
# 5   Tx*(x + x0)**2 + Ty + Tz*z     P0*z*(x + x0)        -Tx*(2*x + 2*x0)/Tz        0          -2*P0*Tx*(x + x0)/Tz - P0*Tx*(2*x + 2*x0)/Tz                -P0*Tx*(2*x + 2*x0)/Tz         px=1;py=0;pz=1;pTx=2;pTy=0;pTz=1;
# 6   Tx*(x + x0)**2 + Ty + Tz*z     P0*z*(x + x0)**2     -Tx*(2*x + 2*x0)/Tz        2*P0*z     -2*P0*Tx*(x + x0)**2/Tz - P0*Tx*(2*x + 2*x0)**2/Tz          -P0*Tx*(2*x + 2*x0)**2/Tz      px=2;py=0;pz=1;pTx=2;pTy=0;pTz=1;
# 7   Tx*(x + x0)**2 + Ty + Tz*z     P0*z**2*(x + x0)     -Tx*(2*x + 2*x0)/Tz        0          -4*P0*Tx*z*(x + x0)/Tz - 2*P0*Tx*z*(2*x + 2*x0)/Tz          -2*P0*Tx*z*(2*x + 2*x0)/Tz     px=1;py=0;pz=2;pTx=2;pTy=0;pTz=1;
# 8   Tx*(x + x0)**2 + Ty + Tz*z     P0*z**2*(x + x0)**2  -Tx*(2*x + 2*x0)/Tz        2*P0*z**2  -4*P0*Tx*z*(x + x0)**2/Tz - 2*P0*Tx*z*(2*x + 2*x0)**2/Tz    -2*P0*Tx*z*(2*x + 2*x0)**2/Tz  px=2;py=0;pz=2;pTx=2;pTy=0;pTz=1;
# 9   Tx*(x + x0) + Ty + Tz*z**2     P0*z*(x + x0)        -Tx/(2*Tz*z)               0          -P0*Tx/(2*Tz*z)                                             0                              px=1;py=0;pz=1;pTx=1;pTy=0;pTz=2;
# 10  Tx*(x + x0) + Ty + Tz*z**2     P0*z*(x + x0)**2     -Tx/(2*Tz*z)               2*P0*z     -P0*Tx*(2*x + 2*x0)/(2*Tz*z)                                0                              px=2;py=0;pz=1;pTx=1;pTy=0;pTz=2;
# 11  Tx*(x + x0) + Ty + Tz*z**2     P0*z**2*(x + x0)     -Tx/(2*Tz*z)               0          -P0*Tx/Tz                                                   -P0*Tx/(2*Tz)                  px=1;py=0;pz=2;pTx=1;pTy=0;pTz=2;
# 12  Tx*(x + x0) + Ty + Tz*z**2     P0*z**2*(x + x0)**2  -Tx/(2*Tz*z)               2*P0*z**2  -P0*Tx*(2*x + 2*x0)/Tz                                      -P0*Tx*(2*x + 2*x0)/(2*Tz)     px=2;py=0;pz=2;pTx=1;pTy=0;pTz=2;
# 13  Tx*(x + x0)**2 + Ty + Tz*z**2  P0*z*(x + x0)        -Tx*(2*x + 2*x0)/(2*Tz*z)  0          -P0*Tx*(x + x0)/(Tz*z) - P0*Tx*(2*x + 2*x0)/(2*Tz*z)        0                              px=1;py=0;pz=1;pTx=2;pTy=0;pTz=2;
# 14  Tx*(x + x0)**2 + Ty + Tz*z**2  P0*z*(x + x0)**2     -Tx*(2*x + 2*x0)/(2*Tz*z)  2*P0*z     -P0*Tx*(x + x0)**2/(Tz*z) - P0*Tx*(2*x + 2*x0)**2/(2*Tz*z)  0                              px=2;py=0;pz=1;pTx=2;pTy=0;pTz=2;
# 15  Tx*(x + x0)**2 + Ty + Tz*z**2  P0*z**2*(x + x0)     -Tx*(2*x + 2*x0)/(2*Tz*z)  0          -2*P0*Tx*(x + x0)/Tz - P0*Tx*(2*x + 2*x0)/Tz                -P0*Tx*(2*x + 2*x0)/(2*Tz)     px=1;py=0;pz=2;pTx=2;pTy=0;pTz=2;
# 16  Tx*(x + x0)**2 + Ty + Tz*z**2  P0*z**2*(x + x0)**2  -Tx*(2*x + 2*x0)/(2*Tz*z)  2*P0*z**2  -2*P0*Tx*(x + x0)**2/Tz - P0*Tx*(2*x + 2*x0)**2/Tz          -P0*Tx*(2*x + 2*x0)**2/(2*Tz)  px=2;py=0;pz=2;pTx=2;pTy=0;pTz=2;
# }}}


def main():
    # {{{
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', dest='input_file',
                        default='base_mesh.nc',
                        help='Input file, containing base mesh'
                        )
    parser.add_argument('-o', '--output_file', dest='output_file',
                        default='initial_state.nc',
                        help='Output file, containing initial variables'
                        )
    parser.add_argument('-t', '--test', dest='test',
                        default=1,
                        type=int,
                        help='Test number 1, 2, or 3'
                        )
    parser.add_argument('-L', '--nVertLevels', dest='nVertLevels',
                        default=10,
                        help='Number of vertical levels'
                        )
    parser.add_argument(
        '-H',
        '--thicknessAllLayers',
        dest='thicknessAllLayers',
        default=10,
        type=float,
        help='thickness of each layer, [m]')
    thicknessAllLayers = parser.parse_args().thicknessAllLayers
    nVertLevels = parser.parse_args().nVertLevels
    test = parser.parse_args().test
    if (test != 1) & (test != 2) & (test != 3):
        print('Invalid choice of test. Setting to test=1.')
        test = 1

    input_file = parser.parse_args().input_file
    output_file = parser.parse_args().output_file
    shutil.copy2(input_file, output_file)
    ds = Dataset(output_file, 'a', format='NETCDF3_64BIT_OFFSET')

    vertical_init(ds, thicknessAllLayers, nVertLevels)
    tracer_init(ds, thicknessAllLayers, parser.parse_args().test)
    velocity_init(ds)
    coriolis_init(ds)
    others_init(ds)

    ds.close()
# }}}


def vertical_init(ds, thicknessAllLayers, nVertLevels):
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

    # evenly spaced vertical grid
    refLayerThickness[:] = thicknessAllLayers
    # make first layer deep to avoid z^2 derivative problems near zero.
    refLayerThickness[0] = 100

    # Create other variables from refLayerThickness
    refBottomDepth[0] = refLayerThickness[0]
    refZMid[0] = -0.5 * refLayerThickness[0]
    for k in range(1, nVertLevels):
        refBottomDepth[k] = refBottomDepth[k - 1] + refLayerThickness[k]
        refZMid[k] = -refBottomDepth[k - 1] - 0.5 * refLayerThickness[k]
    vertCoordMovementWeights[:] = 1.0

    # flat bottom, no bathymetry
    maxLevelCell[:] = nVertLevels
    bottomDepth[:] = refBottomDepth[nVertLevels - 1]
    bottomDepthObserved[:] = refBottomDepth[nVertLevels - 1]
    for k in range(nVertLevels):
        layerThickness[0, :, k] = refLayerThickness[k]
        restingThickness[:, k] = refLayerThickness[k]
# }}}


def tracer_init(ds, thicknessAllLayers, test):
    K = 600.0  # kappa for Redi
    # tracers
    x0 = 1e2
    P0 = 3.0e-3
    Tx = 5 / 1.0e10
    Ty = 0
    Tz = 10.0 / 1000
    S0 = 35
# {{{
    h = thicknessAllLayers
    # create new variables # {{{
    temperature = ds.createVariable(
        'temperature', np.float64, ('Time', 'nCells', 'nVertLevels',))
    salinity = ds.createVariable(
        'salinity', np.float64, ('Time', 'nCells', 'nVertLevels',))
    tracer1 = ds.createVariable(
        'tracer1', np.float64, ('Time', 'nCells', 'nVertLevels',))
    tracer2 = ds.createVariable(
        'tracer2', np.float64, ('Time', 'nCells', 'nVertLevels',))
    tracer3 = ds.createVariable(
        'tracer3', np.float64, ('Time', 'nCells', 'nVertLevels',))
    tracer1Tend = ds.createVariable(
        'tracer1Tend', np.float64, ('Time', 'nCells', 'nVertLevels',))
    tracer2Tend = ds.createVariable(
        'tracer2Tend', np.float64, ('Time', 'nCells', 'nVertLevels',))
    tracer3Tend = ds.createVariable(
        'tracer3Tend', np.float64, ('Time', 'nCells', 'nVertLevels',))
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

    for iCell in range(0, nCells):
        x = xCell[iCell]
        y = yCell[iCell]
        for k in range(0, nVertLevels):
            z = refZMid[k]

            # }}}
            salinity[0, iCell, k] = S0
# TEST 1:
#    Temperature T                  tracer P             Slope Sx                   term 1     term 2                                                      term 3
#1   Tx*(x + x0) + Ty + Tz*z        P0*z*(x + x0)        -Tx/Tz                     0          -P0*Tx/Tz                                                   -P0*Tx/Tz
#2   Tx*(x + x0) + Ty + Tz*z        P0*z*(x + x0)**2     -Tx/Tz                     2*P0*z     -P0*Tx*(2*x + 2*x0)/Tz                                      -P0*Tx*(2*x + 2*x0)/Tz
#3   Tx*(x + x0) + Ty + Tz*z        P0*z**2*(x + x0)     -Tx/Tz                     0          -2*P0*Tx*z/Tz                                               -2*P0*Tx*z/Tz
            if test==1:
                temperature[0,iCell,k] = Tx*(x + x0) + Ty + Tz*z
                tracer1[0,iCell,k] = P0*z*(x + x0)   
                tracer2[0,iCell,k] = P0*z*(x + x0)**2
                tracer3[0,iCell,k] = P0*z**2*(x + x0)
                                           #  term 1     term 2                      term 3
                tracer1Tend[0,iCell,k] = h*K*(0          -P0*Tx/Tz                   -P0*Tx/Tz)
                tracer2Tend[0,iCell,k] = h*K*(2*P0*z     -P0*Tx*(2*x + 2*x0)/Tz      -P0*Tx*(2*x + 2*x0)/Tz)
                tracer3Tend[0,iCell,k] = h*K*(0          -2*P0*Tx*z/Tz               -2*P0*Tx*z/Tz)
# TEST 2:
#    Temperature T                  tracer P             Slope Sx                   term 1     term 2                                                      term 3
#5   Tx*(x + x0)**2 + Ty + Tz*z     P0*z*(x + x0)        -Tx*(2*x + 2*x0)/Tz        0          -2*P0*Tx*(x + x0)/Tz - P0*Tx*(2*x + 2*x0)/Tz                -P0*Tx*(2*x + 2*x0)/Tz
#6   Tx*(x + x0)**2 + Ty + Tz*z     P0*z*(x + x0)**2     -Tx*(2*x + 2*x0)/Tz        2*P0*z     -2*P0*Tx*(x + x0)**2/Tz - P0*Tx*(2*x + 2*x0)**2/Tz          -P0*Tx*(2*x + 2*x0)**2/Tz
#7   Tx*(x + x0)**2 + Ty + Tz*z     P0*z**2*(x + x0)     -Tx*(2*x + 2*x0)/Tz        0          -4*P0*Tx*z*(x + x0)/Tz - 2*P0*Tx*z*(2*x + 2*x0)/Tz          -2*P0*Tx*z*(2*x + 2*x0)/Tz
            elif test==2:
                temperature[0,iCell,k] = Tx*(x + x0)**2 + Ty + Tz*z
                tracer1[0,iCell,k] = P0*z*(x + x0)   
                tracer2[0,iCell,k] = P0*z*(x + x0)**2
                tracer3[0,iCell,k] = P0*z**2*(x + x0)
                                           #  term 1  term 2                                               term 3
                tracer1Tend[0,iCell,k] = h*K*(0       -2*P0*Tx*(x + x0)/Tz - P0*Tx*(2*x + 2*x0)/Tz         -P0*Tx*(2*x + 2*x0)/Tz)
                tracer2Tend[0,iCell,k] = h*K*(2*P0*z  -2*P0*Tx*(x + x0)**2/Tz - P0*Tx*(2*x + 2*x0)**2/Tz   -P0*Tx*(2*x + 2*x0)**2/Tz)
                tracer3Tend[0,iCell,k] = h*K*(0       -4*P0*Tx*z*(x + x0)/Tz - 2*P0*Tx*z*(2*x + 2*x0)/Tz   -2*P0*Tx*z*(2*x + 2*x0)/Tz)
# TEST 3:
#    Temperature T                  tracer P             Slope Sx                   term 1     term 2                                                      term 3
#14  Tx*(x + x0)**2 + Ty + Tz*z**2  P0*z*(x + x0)**2     -Tx*(2*x + 2*x0)/(2*Tz*z)  2*P0*z     -P0*Tx*(x + x0)**2/(Tz*z) - P0*Tx*(2*x + 2*x0)**2/(2*Tz*z)  0
#15  Tx*(x + x0)**2 + Ty + Tz*z**2  P0*z**2*(x + x0)     -Tx*(2*x + 2*x0)/(2*Tz*z)  0          -2*P0*Tx*(x + x0)/Tz - P0*Tx*(2*x + 2*x0)/Tz                -P0*Tx*(2*x + 2*x0)/(2*Tz)
#16  Tx*(x + x0)**2 + Ty + Tz*z**2  P0*z**2*(x + x0)**2  -Tx*(2*x + 2*x0)/(2*Tz*z)  2*P0*z**2  -2*P0*Tx*(x + x0)**2/Tz - P0*Tx*(2*x + 2*x0)**2/Tz          -P0*Tx*(2*x + 2*x0)**2/(2*Tz)
            elif test==3:
                temperature[0,iCell,k] = Tx*(x + x0)**2 + Ty + Tz*z**2
                tracer1[0,iCell,k] = P0*z*(x + x0)**2
                tracer2[0,iCell,k] = P0*z**2*(x + x0)
                tracer3[0,iCell,k] = P0*z**2*(x + x0)**2
                                           #  term 1     term 2                                                      term 3
                tracer1Tend[0,iCell,k] = h*K*(2*P0*z     -P0*Tx*(x + x0)**2/(Tz*z) - P0*Tx*(2*x + 2*x0)**2/(2*Tz*z)  +0)
                tracer2Tend[0,iCell,k] = h*K*(0          -2*P0*Tx*(x + x0)/Tz - P0*Tx*(2*x + 2*x0)/Tz                -P0*Tx*(2*x + 2*x0)/(2*Tz))
                tracer3Tend[0,iCell,k] = h*K*(2*P0*z**2  -2*P0*Tx*(x + x0)**2/Tz - P0*Tx*(2*x + 2*x0)**2/Tz          -P0*Tx*(2*x + 2*x0)**2/(2*Tz))

    # set min T to -1.8
    temperature[:] += -1.8 - np.min(temperature[:])
    # set all tracers to be positive, because of Redi limiter for negative tracers.
    tracer1[:] += - np.min(tracer1[:])
    tracer2[:] += - np.min(tracer2[:])
    tracer3[:] += - np.min(tracer3[:])


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
