#!/usr/bin/env python

"""Generator of html file for machines
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

SRCROOT = os.environ.get("SRCROOT")
if SRCROOT is None:
    SRCROOT = os.path.dirname(CIMEROOT)
os.environ['SRCROOT'] = SRCROOT

from standard_script_setup import *
from CIME.utils import expect
from CIME.XML.entry_id import GenericXML
from CIME.XML.files    import Files
from CIME.XML.machines import Machines

# check for  dependency module
try:
    import jinja2
except:
    raise SystemExit("ERROR: compdef2html.py depends on the jinja2 template module. " /
                     "Install using 'pip --user install jinja2'")

# global variables
_now = datetime.datetime.now().strftime('%Y-%m-%d')
logger = logging.getLogger(__name__)


###############################################################################
def commandline_options():
###############################################################################

    """ Process the command line arguments.                                                                                                                                    
    """
    parser = argparse.ArgumentParser(
        description='Read all the config_machines.xml files and generate a corresponding HTML file.')

    CIME.utils.setup_standard_logging_options(parser)

    parser.add_argument('--htmlfile', nargs=1, required=True,
                        help='Fully quailfied path to output HTML file.')

    parser.add_argument('--version', nargs=1, required=True,
                        help='Model version (e.g. CESM2.0)')

    parser.add_argument('--supported', nargs=1, required=True,
                        help='Comma separated list of supported' \
                            'machines that have passed statistical validation and' \
                            'may be used for large coupled experiments.')

    parser.add_argument('--tested', nargs=1, required=True,
                        help='Comma separated list of tested machines' \
                            'but not necessarily appropriate for large coupled experiments.')

    options = parser.parse_args()

    CIME.utils.parse_args_and_handle_standard_logging_options(options)

    return options

###############################################################################
def _main_func(options, work_dir):
###############################################################################

    """Construct machines html from an XML file."""
        
    # Initialize a variables for the html template
    mach_dict = dict()
    model_version = options.version[0]

    # get the machine config file
    files = Files()
    config_file = files.get_value("MACHINES_SPEC_FILE")
    expect(os.path.isfile(config_file),
           "Cannot find config_file {} on disk".format(config_file))

    # instantiate a machines object and read XML values into a dictionary
    machines = Machines(config_file, machine="Query")
    mach_dict = machines.return_all_values()

    # intialize the supported and tested keys
    for machine in mach_dict.keys():
        mach_dict[machine]['supported'] = False
        mach_dict[machine]['tested'] = False

    # loop through the list of supported machines and flag in the dictionary
    supported = options.supported[0].split(',')
    for machine in supported:
        mach_dict[machine]['supported'] = True

    # loop through the list of tested machines and flag in the dictionary
    tested = options.tested[0].split(',')
    for machine in tested:
        mach_dict[machine]['tested'] = True

    # load up jinja template
    templateLoader = jinja2.FileSystemLoader( searchpath='{0}/templates'.format(work_dir) )
    templateEnv = jinja2.Environment( loader=templateLoader )

    # TODO - get the cesm_version for the CIME root
    tmplFile = 'machdef2html.tmpl'
    template = templateEnv.get_template( tmplFile )
    templateVars = { 'mach_dict'     : mach_dict,
                     'today'         : _now,
                     'model_version' : model_version }
        
    # render the template
    mach_tmpl = template.render( templateVars )

    # write the output file
    with open( options.htmlfile[0], 'w') as html:
        html.write(mach_tmpl)

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




