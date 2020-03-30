def dz_z(z,Hmax,epsilon,dzmax):
    import numpy as np
    return dzmax*np.tanh(-z*np.pi/Hmax) + epsilon

from netCDF4 import Dataset
import numpy as np
"""
Write vertical grid to a netcdf file
"""
def create_vertical_grid(bottom_depth, nz, dz1, dz2, maxit=1000,
    epsilon=1e-2, outFile='MPAS-Ocean_vertical_grid.nc'):

  # open a new netCDF file for writing.
  ncfile = Dataset(outFile,'w')
  # create the depth_t dimension.
  ncfile.createDimension('dim_layer',nz)

  refBottomDepth = ncfile.createVariable('refBottomDepth',np.dtype('float64').char,('dim_layer'))
  refMidDepth = ncfile.createVariable('refMidDepth',np.dtype('float64').char,('dim_layer'))
  refLayerThickness = ncfile.createVariable('refLayerThickness',np.dtype('float64').char,('dim_layer'))

  Hmax = bottom_depth
  nLayers = 0

  layerThickness = np.zeros(nz)
  dz = [epsilon]
  z = [0]
  count=0

  while nLayers != nz and count < maxit:
    zval = -epsilon
    dz = [epsilon]
    z = [0]
    nLayers = 0
    while zval > -Hmax:
        difference = dz_z(zval,Hmax,epsilon,dz2) - zval
        while abs(difference) > 0.3:
            zval -= epsilon
            difference = dz_z(zval,Hmax,epsilon,dz2) -z[nLayers] + zval
        z.append(zval)
        dz.append(dz_z(zval,Hmax,epsilon,dz2))
        nLayers += 1
        zval -= epsilon

    dz_arr = np.asarray(dz)
    ind = abs(dz_arr - dz1).argmin()

    dztemp = dz_arr[ind:]
    nLayers = len(dztemp)
    change = nz-nLayers
    dz2 -= float(change)
    count += 1

  if count >= maxit:
    print("Error: grid did not converge, adjust parameters")
  else:
    layerThickness[:nLayers] = dztemp

    nVertLevels = nz 
    botDepth = np.zeros(nVertLevels)
    midDepth = np.zeros(nVertLevels)
    botDepth[0] = layerThickness[0]
    midDepth[0] = 0.5*layerThickness[0]

    for i in range(1,nVertLevels):
      botDepth[i] = botDepth[i-1] + layerThickness[i]
      midDepth[i] = midDepth[i-1] + 0.5*(layerThickness[i] + layerThickness[i-1])

    refBottomDepth[:] = botDepth
    refMidDepth[:] = midDepth
    refLayerThickness[:] = layerThickness[:nVertLevels]

    ncfile.close()

if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument("-o", "--output_file_pattern", dest="output_filename_pattern",
      default='MPAS-Ocean_vertical_grid.nc',
      help="MPAS Filename pattern for file output of vertical grid.", metavar="NAME")
  parser.add_argument("-bd", "--bottom_depth", dest="bottom_depth",
      help="bottom depth for the chosen vertical coordinate", 
      type=float,required=True)
  parser.add_argument("-nz", "--num_vert_levels", dest="nz",
      help="Number of vertical levels for the grid", type=int,
      required=True)
  parser.add_argument("-dz1", "--layer1_thickness", dest="dz1",
      help="Target thickness of the first layer", 
      type=float, required=True)
  parser.add_argument("-dz2", "--maxLayer_thickness", dest="dz2",
      help="Target maximum thickness in column", 
      type=float, required=True)
  parser.add_argument("-eps", "--error_tolerance", dest="epsilon",
      help="Threshold for iterations", type=float, default=1e-2)
  parser.add_argument("-maxit", "--max_iterations", dest="maxit",
      help="maximum number of iterations for grid convergences",
      type=int, default=1000)
  args = parser.parse_args()

  create_vertical_grid(args.bottom_depth,args.nz,args.dz1,
      args.dz2,maxit=args.maxit,epsilon=args.epsilon,
      outFile=args.output_filename_pattern)
