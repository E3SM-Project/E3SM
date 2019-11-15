#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import sys
import subprocess
import argparse
import numpy
from netCDF4 import Dataset


parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "--iteration_count",
    dest="iteration_count",
    default=1,
    type=int,
    help="The number of iterations between init and forward mode for computing"
         " a balanced land-ice pressure.")
parser.add_argument(
    "--first_iteration",
    dest="first_iteration",
    default=0,
    type=int,
    help="The iteration to start from (for continuing iteration if iterrupted "
         "or insufficient)")
parser.add_argument(
    "--plot_globalStats",
    dest="plot_globalStats",
    action='store_true',
    help="If present, plot mean and max KE, min layer thickness and mean "
         "temperature for debugging.")
parser.add_argument(
    "--variable_to_modify",
    dest="variable_to_modify",
    default='ssh',
    help="Which variable, either ssh or landIcePressure, to modify at each "
         "iteration.")

args = parser.parse_args()
dev_null = open('/dev/null', 'w')
error = False

if args.variable_to_modify not in ['ssh', 'landIcePressure']:
    print("Error: unknown variable to modify", args.variable_to_modify)

if args.plot_globalStats:
    subprocess.check_call(['mkdir',
                           '-p',
                           'statsPlots'],
                          stdout=dev_null,
                          stderr=dev_null)

for iterIndex in range(args.first_iteration, args.iteration_count):
    print(" * Iteration %i/%i" % (iterIndex + 1, args.iteration_count))

    subprocess.check_call(['ln',
                           '-sfn',
                           'init%i.nc' % iterIndex,
                           'init.nc'],
                          stdout=dev_null,
                          stderr=dev_null)

    print("   * Running forward model")
    subprocess.check_call(
        ['./run_model.py'],
        stdout=dev_null,
        stderr=dev_null)
    print("   - Complete")

    if args.plot_globalStats:
        print("   * Plotting stats")
        subprocess.check_call(
            ['./plot_globalStats.py', '--out_dir=statsPlots',
             '--iteration={}'.format(iterIndex), 'kineticEnergyCellMax',
             'kineticEnergyCellAvg', 'layerThicknessMin'], stdout=dev_null,
            stderr=dev_null)
        print("   - Complete")

    print("   * Updating SSH or land-ice pressure")

    # copy the init file first

    subprocess.check_call(['cp', 'init{}.nc'.format(iterIndex),
                           'init{}.nc'.format(iterIndex + 1)],
                          stdout=dev_null, stderr=dev_null)
    subprocess.check_call(['ln',
                           '-sfn',
                           'init%i.nc' % (iterIndex + 1),
                           'init.nc'],
                          stdout=dev_null,
                          stderr=dev_null)
    initFile = Dataset('init.nc', 'r+')

    nVertLevels = len(initFile.dimensions['nVertLevels'])
    initSSH = initFile.variables['ssh'][0, :]
    bottomDepth = initFile.variables['bottomDepth'][:]
    modifySSHMask = initFile.variables['modifySSHMask'][0, :]
    landIcePressure = initFile.variables['landIcePressure'][0, :]
    lonCell = initFile.variables['lonCell'][:]
    latCell = initFile.variables['latCell'][:]
    maxLevelCell = initFile.variables['maxLevelCell'][:]

    inSSHFile = Dataset('output_ssh.nc', 'r')
    nTime = len(inSSHFile.dimensions['Time'])
    finalSSH = inSSHFile.variables['ssh'][nTime - 1, :]
    topDensity = inSSHFile.variables['density'][nTime - 1, :, 0]
    inSSHFile.close()

    mask = numpy.logical_and(maxLevelCell > 0, modifySSHMask == 1)

    deltaSSH = mask * (finalSSH - initSSH)

    # then, modifty the SSH or land-ice pressure
    if args.variable_to_modify == 'ssh':
        initFile.variables['ssh'][0, :] = finalSSH
        # also update the landIceDraft variable, which will be used to
        # compensate for the SSH due to land-ice pressure when computing
        # sea-surface tilt
        initFile.variables['landIceDraft'][0, :] = finalSSH
        # we also need to stretch layerThickness to be compatible with the new
        # SSH
        stretch = (finalSSH + bottomDepth) / (initSSH + bottomDepth)
        layerThickness = initFile.variables['layerThickness']
        for k in range(nVertLevels):
            layerThickness[0, :, k] *= stretch
    else:
        # Moving the SSH up or down by deltaSSH would change the land-ice
        # pressure by density(SSH)*g*deltaSSH. If deltaSSH is positive (moving
        # up), it means the land-ice pressure is too small and if deltaSSH is
        # negative (moving down), it means land-ice pressure is too large, the
        # sign of the second term makes sense.
        gravity = 9.80616
        deltaLandIcePressure = topDensity * gravity * deltaSSH

        landIcePressure = numpy.maximum(
            0.0, landIcePressure + deltaLandIcePressure)

        initFile.variables['landIcePressure'][0, :] = landIcePressure

        finalSSH = initSSH

    initFile.close()

    # Write the largest change in SSH and its lon/lat to a file
    logFile = open('maxDeltaSSH_%03i.log' % iterIndex, 'w')

    indices = numpy.nonzero(landIcePressure)[0]
    index = numpy.argmax(numpy.abs(deltaSSH[indices]))
    iCell = indices[index]
    logFile.write('deltaSSHMax: {:g}, lon/lat: {:f} {:f}, ssh: {:g}, '
                  'landIcePressure: {:g}\n'.format(
                      deltaSSH[iCell], numpy.rad2deg(lonCell[iCell]),
                      numpy.rad2deg(latCell[iCell]), finalSSH[iCell],
                      landIcePressure[iCell]))
    logFile.close()

    print("   - Complete")

if error:
    sys.exit(1)
else:
    sys.exit(0)
