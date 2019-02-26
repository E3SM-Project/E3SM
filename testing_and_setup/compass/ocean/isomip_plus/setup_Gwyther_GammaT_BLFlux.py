#!/usr/bin/env python
"""
Sets up a paramter study that varies GammaT and GammaS (the
heat and salt transfer coefficients) according to the
specification of the ISOMIP+ Ocean0 experiment
"""

import numpy
import subprocess

gammaTs = numpy.array(
    [0.01, 0.02, 0.03, 0.04, 0.05, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11,
     0.12, 0.15, 0.2])
#gammaTs = numpy.array([0.01])
gammaSs = gammaTs/35.

blFlux = numpy.array([2., 5., 10., 15., 20., 25., 30., 35., 40.])
#blFlux = numpy.array([2.])

GammaTs, BLFluxs = numpy.meshgrid(gammaTs, blFlux)
GammaSs, _ = numpy.meshgrid(gammaSs, blFlux)

GammaTStrings = ['{:g}'.format(GammaT) for GammaT in GammaTs.ravel()]
GammaSStrings = ['{:g}'.format(GammaS) for GammaS in GammaSs.ravel()]
BLFluxStrings = ['{:g}'.format(depth) for depth in BLFluxs.ravel()]

args = ['../../utility_scripts/make_parameter_study_configs.py', '-t', 
        'template_Gwyther_GammaT_BLFlux.xml', '-o', '2km/Ocean0/config_GammaT_BLFlux',
        '-p', 'GammaT={}'.format(','.join(GammaTStrings)),
        'GammaS={}'.format(','.join(GammaSStrings)), 
         'blFlux={}'.format(','.join(BLFluxStrings))]
        
print 'running {}'.format(' '.join(args))
subprocess.check_call(args)
