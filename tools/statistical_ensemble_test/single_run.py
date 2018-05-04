#! /usr/bin/env python 
import sys, getopt, os 
import numpy as np 
import Nio 
import utility
import time


#
# Command options 
#
def disp_usage(callType):
    print '\n Sets up CESM cases for either an ensemble of runs or a small (size 3)\n'
    print 'test set (default). Then Use pyCECT utilities to create an ensemble \n'
    print 'summary file or to evaluate the small test set of runs against the ensemble.\n'
    print '  '
    print '  ----------------------------'
    print '   Args for ensemble.py :'
    print '  ----------------------------'
    print '   ensemble.py'
    print '   -h                      : prints out this usage message'


    print 'Required flags:'
    print '  --case <name>    Case name passed on to create_newcase'
    print '  --mach <name>    Machine name passed on to create_newcase'
    print ' ' 
    print 'Optional flags (+ all "--" options to create_newcase): '
    print '  --project <num>  Project number to charge in job scripts'
    if callType == 'single_run.py':
       print '  --pertlim <num>     Run CAM with non-zero pertlim'
   print '  --walltime <hr:mn> Amount of walltime requested (default = 4:30, or 0:10 with --uf enabled)'
   print '  --compiler <name>  Compiler to use (default = same as Machine default) '
   print '  --compset <name>   Compset to use (default = F2000)'
   print '  --res <name>       Resolution to run (default = f19_f19)'
   print '  --uf               Enable ninth time step runs (ultra-fast mode) - otherwise the default is 12-month runs'
   if callType == 'ensemble.py': 
       print '  --nb               Disables auto building the root case of the ensemble.'
       print '  --ns               Disables auto submitting any members of the ensemble.'
       print '  --ensemble <size>  Build the ensemble (instead of building 3 cases with random pertlim values for verification),'
       print '                     and specify the number of ensemble members to generate (e.g.: 151 for annual averages or 350 for ultra-fast mode)'
   else:
       print '  --nb               Disables building (and submitting) the single case.'
       print '  --ns               Disables submitting the single case.'

