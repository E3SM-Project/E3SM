#!/usr/bin/env python2.7

from Tools.standard_script_setup import *
from CIME.utils        import expect, run_cmd
from CIME.case         import Case

from CIME.XML.files    import Files
from CIME.XML.compsets import Compsets
from CIME.XML.grids    import Grids

import argparse, doctest, tempfile

logger = logging.getLogger(__name__)

###############################################################################
def get_compset_longname(compset_name):
###############################################################################

    files = Files()
    components = files.get_components("COMPSETS_SPEC_FILE")
    logger.info("components are %s" % components)
    
    # Loop through all of the files listed in COMPSETS_SPEC_FILE and find the file
    # that has a match for either the alias or the longname in that order
    for component in components:
        
        # Determine the compsets file for the target component
        file = files.get_value("COMPSETS_SPEC_FILE", {"component":component})
        
        # If the file exists, read it and see if there is a match for the compset alias or longname
        if (os.path.isfile(file)):
            compsets = Compsets(file)
            match = compsets.get_compset_match(name=compset_name)
            if match is not None:
                compset_lname = match
                logger.debug("Successful compset match %s found in file %s " %(match, file))
                return match

        logger.debug("Could not find a compset match for either alias or longname in %s" %file)
                
###############################################################################
def get_grid_longname(grid_name, compset=None):
###############################################################################

    files = Files()
    gridfile = files.get_value("GRIDS_SPEC_FILE")
    logger.info("Grid specification file is %s" % gridfile)

    grids = Grids(gridfile)
    match = grids.get_grid_match(name=grid_name, compset=compset)

    expect ((match is not None),
            "No valid grid found for input grid %s and compset %s" %(grid_name, compset))

    return match

###############################################################################
def get_compset_components(compset):
###############################################################################

    compsetRE = re.compile(r"_")
    components = compsetRE.split(compset)
    print "DEBUG: components are ",components

###############################################################################
#def _main_func():
###############################################################################

cime_root  = CIME.utils.get_cime_root()
cime_model = CIME.utils.get_model()

compset_input = "A"
grid_input    = "f19_g16_rx1"

casedir = "./foop"

# make case directory
cwd = os.getcwd()
casedir = cwd + "/foop"
os.mkdir(casedir)
os.chdir(casedir)

# write xml files in case directory
case = Case(case_root=casedir)
case.configure(compset_input, grid_input)
for file in case._env_entryid_files:
    file.write()




