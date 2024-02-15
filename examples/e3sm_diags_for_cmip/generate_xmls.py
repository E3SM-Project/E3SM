# Script to generate xml files

import glob
import os.path
import shlex
import subprocess

from datetime import date

def run_command(command):
    p1 = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout, stderr) = p1.communicate()
    status = p1.returncode
    if status != 0:
        print('stdout:\n%s', stdout)
        print('stderr:\n%s', stderr)
        raise Exception


# Experiment ro process
experiment = "amip"
#experiment = "historical"

# Dry-run? (only scan available files)
dry_run = False

# Output directory for xml files
destination = '/home/zhang40/e3sm_diags_for_CMIP6/CMIP6_20240111'

# Input search paths
paths = ('/p/user_pub/work/CMIP6/CMIP','')  # for E3SM CMIP archive only
#paths = ('/p/css03/esgf_publish/CMIP6/CMIP', '/p/user_pub/work/CMIP6/CMIP')

# Search patterm
patterns = ('*/*/%s/r1i1p1f1/Amon/' % (experiment),)
#patterns = ('UCSB/*/%s/r*i*p*f*/Amon/' % (experiment),)

# Output file to log included input netCDF files
name = "%s_%s.log" % (experiment,date.today().strftime("%y%m%d"))
flog = open(name ,"w")
print(name)

# Find available simulations
simulations = []
for path in paths:
    for pattern in patterns:
        simulations.extend(glob.glob(os.path.join(path, pattern)))

# Loop over simulations
#for simulation in simulations[0:3]:
for simulation in simulations:
    print('\n=== %s ===' % (simulation))
    # Split source path, remove empty strings
    p = list(filter(None, simulation.split('/')))
    # Destination directory
    dest = os.path.join(destination, *p[-6:-1])
    if not os.path.exists(dest):
        os.makedirs(dest)
    # Create xml files for each available variable
    vars = sorted(os.listdir(simulation))
    # Skip is certain specific variables are not found
    if not set(['pr',]).issubset(vars):
        print("Skipping...")
        continue
    for var in vars:
        versions = glob.glob(os.path.join(simulation,var,'g*/v[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]'))
        versions = sorted(versions)
        mostRecent = versions[-1]
        # Now, extract first and last date
        print('most recent version', mostRecent)
        files = sorted(os.listdir(mostRecent))
        # First time stamp
        first = files[0].split('_')[-1]
        first = first[0:6]
        # Last time stamp
        last = files[-1].split('_')[-1]
        last = last[7:13]
        # cdscan command
        inputFiles = glob.glob(mostRecent+'/*.nc')
        # Output file names
        for f in inputFiles:
            flog.write(f+'\n')
        inputFiles = ' '.join(inputFiles)
        command = 'cdscan -x %s/%s_%s_%s.xml %s' % (dest,var,first,last,inputFiles)
        if not dry_run:
            print(command)
            run_command(command)

# Close log output file
flog.close()
