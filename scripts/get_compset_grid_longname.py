#!/usr/bin/env python2.7

from Tools.standard_script_setup import *
from CIME.utils        import expect, run_cmd
from CIME.case         import Case

import argparse, doctest, tempfile

logger = logging.getLogger(__name__)
                
###############################################################################
#def _main_func():
###############################################################################

# hard-wired for now - case, compset and grid input 
compset_input = "A"
grid_input    = "f19_g16_rx1"
casename      = "foop"

# Determine 
if os.path.exists(casename):
    caseroot = os.path.dirname(casename)
    casename = os.path.basename(casename)
else:
    caseroot = os.path.join(os.getcwd(),casename)

# Set the case object
case = Case(caseroot)

# Set values for env_case.xml 
case.set_value("CASE",casename)
case.set_value("CASEROOT",caseroot)

cimeroot  = CIME.utils.get_cime_root()
srcroot  = os.path.join(cimeroot,"/../")
case.set_value("SRCROOT",srcroot);

# Configure the Case
case.configure(compset_input, grid_input)

# Make the case directory
os.mkdir(caseroot)
os.chdir(caseroot)

# Write out the case files
for file in case._env_entryid_files:
    print "file is ",file
    file.write()







