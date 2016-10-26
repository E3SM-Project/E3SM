#!/usr/bin/env python

import argparse
import os
import subprocess
import glob

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-f", "--folder", dest="folder", help="Folder for plots", required=True)
parser.add_argument("-e", "--expt", dest="expt", help="Experiment number (0, 1 or 2)", required=True)
parser.add_argument("-a", "--avconv", dest="avconv", help="full path to avconv")

args = parser.parse_args()


folder = args.folder
expt = args.expt

if(args.avconv is None):
    avconv = 'avconv'
else:
    avconv = args.avconv

initFile = 'init.nc'
print folder
print 'expt ', expt

try:
    os.makedirs('%s/plotLogs'%folder)
except OSError:
    pass


logFile = open('%s/plotLogs/barotropic.log'%folder, 'w')
args = ['./viz/computeBarotropicStreamfunction.py', folder]

print 'running %s'%' '.join(args)
barotropicProcess = subprocess.Popen(args, stdout=logFile, stderr=logFile)


logFile = open('%s/plotLogs/overturning.log'%folder, 'w')
args = ['./viz/computeOverturningStreamfunction.py', folder]

print 'running %s'%' '.join(args)
overturningProcess = subprocess.Popen(args, stdout=logFile, stderr=logFile)

status = barotropicProcess.wait()

if(status != 0):
    print 'Error running computeBarotropicStreamfunction.py'
    exit(1)

status = overturningProcess.wait()

if(status != 0):
    print 'Error running computeOverturningStreamfunction.py'
    exit(1)

logFile = open('%s/plotLogs/plot.log'%folder, 'w')
args = ['./viz/plotResults.py', '--inFolder=%s'%folder, '--outImageFolder=%s/plots'%folder,
        '--expt=%s'%expt]

print 'running %s'%' '.join(args)
subprocess.check_call(args, stdout=logFile, stderr=logFile)


framesPerSecond = '30'
try:
    os.makedirs('%s/movies'%folder)
except OSError:
    pass


processes = []

for fileName in glob.glob('%s/plots/*0001.png'%folder):
    prefix = fileName.split('/')[-1][:-9]
    print prefix
    logFile = open('%s/plotLogs/%s.log'%(folder,prefix), 'w')
    args = [avconv, '-y', '-r', framesPerSecond, '-i', '%s/plots/%s_%%04d.png'%(folder, prefix),
            '-b', '32000k', '-r', framesPerSecond, '%s/movies/%s.mp4'%(folder, prefix)]
 
    print 'running %s'%' '.join(args)
    process = subprocess.Popen(args, stdout=logFile, stderr=logFile)
    
    processes.append(process)

for process in processes:
    status = process.wait()
    
    if(status != 0):
        print 'Error running avconv'

