#!/usr/bin/env python

# Set up subdirectories for the various MISMIP+ experiments:
# Ice0, Ice1r, Ice1ra, Ice1rr, Ice1rax, Ice1rrx, Ice2r, Ice2ra, Ice2r, Ice2rax, Ice2rrx
# Note: Ice1rax is the optional extension of Ice1ra from year 200 to 1000,
#       and similarly for the other Ice*x experiments.

# The namelist and streams files for each experiment should already
# have been created in the directory from which this script is launched;
# the script simply moves them to the subdirectories.

import sys, os
import shutil

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-x", "--expt", dest="experiment", type='string', help="MISMIP+ experiment(s) to set up", metavar="EXPT")
options, args = parser.parse_args()

if options.experiment:
    if options.experiment == 'all':
        # Set up subdirectories for all experiments
        experiments = ['Spinup', 'Ice0', 'Ice1r', 'Ice1ra', 'Ice1rr', 'Ice1rax', 'Ice1rrx', 'Ice2r', 'Ice2ra', 'Ice2rr', 'Ice2rax', 'Ice2rrx']
    else:
        experiments = [options.experiment]
else:
    sys.exit('Error: No experiment specified.  Please specify experiment(s) with the -x option')


print 'Experiments:', experiments

# Loop through experiments
for expt in experiments:
    print 'Setting up directory for experiment', expt

    # Make the subdirectory if it does not exist already
    try: 
        os.mkdir(expt)
    except:
        pass

    # Go to the subdirectory
    try:
        os.chdir(expt)
    except:
        sys.exit('Error, could not change to subdirectory')

    # Move the appropriate namelist and stream files from the parent directory.
    # Note: In the subdirectory, the expt prefix (e.g., 'Ice0') is not included.
    #       So there is no need for the -n and -s specifiers when launching a run.
    namelistFile = '../namelist.landice.' + expt
    shutil.move(namelistFile, './namelist.landice')

    streamsFile = '../streams.landice.' + expt
    shutil.move(streamsFile, './streams.landice')

    # Link to the executable in the parent directory
    executableName = 'landice_model'
    os.symlink('../landice_model', 'landice_model')

    # Link to any and all graph partition files in the parent directory
    # Note: If a new file is needed, it can be created using metis.
    #       For example to run on 128 cores on LANL IC:
    #       > module load metis
    #       > gpmetis graph.info 128
    #       This creates a file called graph.info.part.128

    for file in os.listdir('..'):
        if file.startswith('graph.info.part'):
            os.symlink('../' + file, file)

    # Link to the albany input file in the parent directory
    os.symlink('../' + 'albany_input.xml', 'albany_input.xml')

    # Link to the appropriate restart file and timestamp
    # No restart file needed for the Spinup experiment
    # Note: The symlinks will initially be empty (except for landice_grid.nc).
    #       Ice0, Ice1r and Ice2r must follow the Spinup.
    #       Ice1ra and Ice1rr must follow Ice1r; Ice2ra and Ice2rr must follow Ice2r.
    #       Ice1rax must follow Ice1ra, and similiary for the other Ice*x.

    if expt == 'Spinup':
        # Start from landice_grid.nc
        gridfile = 'landice_grid.nc'
        griddir = '../'
        os.symlink(griddir + gridfile, gridfile)
    elif expt =='Ice0' or expt=='Ice1r' or expt=='Ice2r':
        # Start from restart file at the end of Spinup, but call it landice_grid.nc,
        #  so the run is treated as a cold start
        # Note: This requires a one-line NCO command to rename the Spinup restart file
        #       while removing the xtime variable
        gridfile = 'landice_grid.nc'
        griddir = '../Spinup/'
        os.symlink(griddir + gridfile, gridfile)    
    else:
        # Start from the appropriate restart file
        if expt=='Ice1ra' or expt=='Ice1rr':
            restartYear = 100
            restartdir = '../Ice1r/'
        elif expt=='Ice1rax':
            restartYear = 200
            restartdir = '../Ice1ra/'
        elif expt=='Ice1rrx':
            restartYear = 200
            restartdir = '../Ice1rr/'
        elif expt=='Ice2ra' or expt=='Ice2rr':
            restartYear = 100
            restartdir = '../Ice2r/'
        elif expt=='Ice2rax':
            restartYear = 200
            restartdir = '../Ice2ra/'
        elif expt=='Ice2rrx':
            restartYear = 200
            restartdir = '../Ice2rr/'

        # Link to the restart file
        restartfile = 'restart_00' + str(restartYear) + '.nc'
        os.symlink(restartdir + restartfile, restartfile)

        # Create the restart_timestamp file
        # Not using symbolic links because these allow files to be rewritten
        #  from other directories
        timestampFile = open('restart_timestamp', 'w')
        restartTimestamp = ' ' + str(restartYear) + '-01-01_00:00:00'
        timestampFile.write(restartTimestamp + '\n')
        timestampFile.close()

    # Go back to the main directory and continue
    os.chdir('..')
