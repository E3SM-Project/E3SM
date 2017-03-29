#!/usr/bin/env python

"""Generator of html file for compsets
"""

# Typically ignore this.
# pylint: disable=invalid-name

# Disable these because this is our standard setup
# pylint: disable=wildcard-import,unused-wildcard-import,wrong-import-position

import os, sys, re, glob
import datetime

CIMEROOT = os.environ.get("CIMEROOT")
if CIMEROOT is None:
    raise SystemExit("ERROR: must set CIMEROOT environment variable")
sys.path.append(os.path.join(CIMEROOT, "scripts", "Tools"))

SRCROOT = os.environ.get("SRTROOT")
if SRCROOT is None:
    SRCROOT = os.path.dirname(CIMEROOT)
os.environ['SRCROOT'] = SRCROOT

os.environ['MODEL'] = "cesm"

from standard_script_setup import *
from CIME.utils import expect
from CIME.XML.entry_id import GenericXML
from CIME.XML.files    import Files
from CIME.XML.compsets import Compsets

# check for  dependency module
try:
    import jinja2
except:
    raise SystemExit("ERROR: compdef2html.py depends on the jinja2 template module. " /
                     "Install using 'pip --user install jinja2'")

# global variables
_now = datetime.datetime.now().strftime('%Y-%m-%d')
_comps = ['CAM', 'CLM', 'CISM', 'POP2', 'CICE', 'RTM', 'MOSART', 'WW3', 
          'Driver', 'DATM', 'DESP', 'DICE', 'DLND', 'DOCN', 'DROF', 'DWAV']
logger = logging.getLogger(__name__)


###############################################################################
def commandline_options():
###############################################################################

    """ Process the command line arguments.                                                                                                                                    
    """
    parser = argparse.ArgumentParser(
        description='Read all the config_compset.xml files and generate a corresponding HTML file.')

    CIME.utils.setup_standard_logging_options(parser)

    parser.add_argument('--configfile', nargs=1, required=True,
                        help='Fully quailfied path to config_files.xml.')

    parser.add_argument('--htmlfile', nargs=1, required=True,
                        help='Fully quailfied path to output HTML file.')

    options = parser.parse_args()

    CIME.utils.handle_standard_logging_options(options)

    return options

###############################################################################
def _main_func(options, work_dir):
###############################################################################

    """Construct compsets from an XML file."""

    # Create a definition object from the xml file
    filename = options.configfile[0]
    expect(os.path.isfile(filename), "File %s does not exist"%filename)
    try:
        definition = GenericXML(infile=filename)
    except:
        sys.exit("Error: unable to parse file %s" %filename)
        
    # Initialize a variables for the html template
    all_compsets = dict()
    html_dict = dict()
    cesm_version = 'CESM2.0'

    # read in all the component namelist files
    files = Files()
    components = files.get_components("COMPSETS_SPEC_FILE")   

    # loop through the components and read in the config_compset.xml file
    for component in components:
        # Determine the compsets file for this component
        compsets_filename = files.get_value("COMPSETS_SPEC_FILE", {"component":component})

        # check to make sure compsets_filename exists before opening and reading
        if (os.path.isfile(compsets_filename)):
            compsets = Compsets(compsets_filename)
            all_compsets = compsets.return_all_values()
## START HERE 
## need to also get info from compsets.get_compset_match and compsets.get_compset_var_settings





    # load up jinja template
    templateLoader = jinja2.FileSystemLoader( searchpath='{0}/templates'.format(work_dir) )
    templateEnv = jinja2.Environment( loader=templateLoader )

    # TODO - get the cesm_version for the CIME root
    tmplFile = 'nmldef2html.tmpl'
    template = templateEnv.get_template( tmplFile )
    templateVars = { 'html_dict'    : html_dict,
                     'today'        : _now,
                     'cesm_version' : cesm_version }
        
    # render the template
    nml_tmpl = template.render( templateVars )

    # write the output file
    with open( options.htmlfile[0], 'w') as html:
        html.write(nml_tmpl)

    return 0

###############################################################################

if __name__ == "__main__":

    options = commandline_options()
    work_dir = os.path.join(CIMEROOT,"scripts","Tools","xml2html")
    try:
        status = _main_func(options, work_dir)
        sys.exit(status)
    except Exception as error:
        print(str(error))
        sys.exit(1)




