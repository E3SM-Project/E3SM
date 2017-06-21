#!/usr/bin/env python

"""Generator of html file for a given caseroot env_*.xml input file 
"""

# Typically ignore this.
# pylint: disable=invalid-name

# Disable these because this is our standard setup
# pylint: disable=wildcard-import,unused-wildcard-import,wrong-import-position

import os, sys
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

# check dependency modules
try:
    import jinja2
except:
    raise SystemExit("ERROR: casedef2html.py depends on the python jinja2 template module. " /
                     "Install using 'pip --user install jinja2'")
try:
    import lxml.etree as ET
except:
    raise SystemExit("ERROR: casedef2html.py depends on the python lxml module. " /
                     "Install using 'pip --user install lxml'")

# global variables
_now = datetime.datetime.now().strftime('%Y-%m-%d')
logger = logging.getLogger(__name__)

###############################################################################
def commandline_options():
###############################################################################

    """ Process the command line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Read a CASEROOT env_*.xml file and generate a corresponding HTML file ' \
            'using a corresponding XSL style sheet transformation.')

    CIME.utils.setup_standard_logging_options(parser)

    parser.add_argument('--htmlfile', nargs=1, required=True,
                        help='Fully quailfied path to output HTML file.')

    parser.add_argument('--xmlfile', nargs=1, required=True,
                        help='Fully qualified input XML file to be transformed to html.')

    parser.add_argument('--xslfile', nargs=1, required=True,
                        help='Fully qualified path to input XSL file corresponding to the input xmlfile.')

    options = parser.parse_args()

    CIME.utils.parse_args_and_handle_standard_logging_options(options)

    return options

###############################################################################
def _main_func(options, work_dir):
###############################################################################

    """Construct html from an XML file."""

    # check the xml file
    xmlfile = options.xmlfile[0]
    expect(os.path.isfile(xmlfile),
           "Cannot find xml file {} on disk".format(xmlfile))

    # check the corresponding xsl stylesheet transform file
    xslfile = options.xslfile[0]
    expect(os.path.isfile(xslfile),
           "Cannot find xsl file {} on disk".format(xslfile))

    # transform
    dom = ET.parse(xmlfile)
    xslt = ET.parse(xslfile)
    transform = ET.XSLT(xslt)
    newdom = transform(dom)

    # output transformed html
    outfile = open(options.htmlfile[0], 'w')
    outfile.write(ET.tostring(newdom, pretty_print=True))
    outfile.close()

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




