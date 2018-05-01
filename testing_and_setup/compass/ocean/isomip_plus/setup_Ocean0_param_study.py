#!/usr/bin/env python
"""
Sets up a paramter study that varies GammaT and GammaS (the
heat and salt transfer coefficients) according to the
specification of the ISOMIP+ Ocean0 experiment
"""

import numpy
import subprocess

GammaTs = numpy.linspace(0.002,0.02,10)
GammaSs = GammaTs/35.

GammaTStrings = []
GammaSStrings = []
for index in range(len(GammaTs)):
    GammaTStrings.append('%g'%GammaTs[index])
    GammaSStrings.append('%g'%GammaSs[index])

args = ['../../utility_scripts/make_parameter_study_configs.py', '-t', 
        'template_Ocean0_param_study.xml', '-o', '2km/Ocean0/config_GammaT',
        '-p', 'GammaT=%s'%','.join(GammaTStrings),
        'GammaS=%s'%','.join(GammaSStrings)]
        
print 'running %s'%' '.join(args)
subprocess.check_call(args)
