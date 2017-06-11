#!/usr/bin/env python

"""Generator of html file for compsets
"""

# Typically ignore this.
# pylint: disable=invalid-name

# Disable these because this is our standard setup
# pylint: disable=wildcard-import,unused-wildcard-import,wrong-import-position

import os, sys, re, glob
import datetime
import collections

CIMEROOT = os.environ.get("CIMEROOT")
if CIMEROOT is None:
    raise SystemExit("ERROR: must set CIMEROOT environment variable")
sys.path.append(os.path.join(CIMEROOT, "scripts", "Tools"))

SRCROOT = os.environ.get("SRCROOT")
if SRCROOT is None:
    SRCROOT = os.path.dirname(CIMEROOT)
os.environ['SRCROOT'] = SRCROOT

from standard_script_setup import *
from CIME.utils            import expect, get_model, get_cime_root
from CIME.XML.entry_id     import GenericXML
from CIME.XML.files        import Files
from CIME.XML.component    import Component
from CIME.XML.compsets     import Compsets
from CIME.case             import Case
import xml.etree.ElementTree as etree

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
        description='Read all the config_compset.xml files and generate a corresponding HTML file.')

    CIME.utils.setup_standard_logging_options(parser)

    parser.add_argument('--htmlfile', nargs=1, required=True,
                        help='Fully quailfied path to output HTML file.')

    parser.add_argument('--version', nargs=1, required=True,
                        help='Model version (e.g. CESM2.0)')

    options = parser.parse_args()

    CIME.utils.parse_args_and_handle_standard_logging_options(options)

    return options

###############################################################################
def get_descs(desc_file):
###############################################################################

    descs = dict()
    ordered_descs = collections.OrderedDict()

    # check that desc_file exists
    expect(os.path.isfile(desc_file), "File %s does not exist"%desc_file)

    xml_tree = etree.ElementTree()
    xml_tree.parse(desc_file)
    descrips = xml_tree.findall("description/desc")

    for descrip in descrips:
        compset_attr = descrip.get("compset")
        descs[compset_attr] = descrip.text
        
        # parse out the warning texts
        descs[compset_attr] = descs[compset_attr].replace('-----------------------------WARNING ------------------------------------------------',
                                                                       ' -- WARNING -- ')
        descs[compset_attr] = descs[compset_attr].replace('-------------------------------------------------------------------------------------',
                                                                       '')

    ordered_descs = collections.OrderedDict(sorted(descs.items(), key=lambda t: t[0]))

    return ordered_descs

###############################################################################
def _main_func(options, work_dir):
###############################################################################

    """Construct compsets html from an XML file."""
        
    # Initialize variables
    ordered_compsets = collections.OrderedDict()
    html_dict = dict()
    model_version = options.version[0]

    # read in all the component config_compsets.xml files
    files = Files()
    components = files.get_components("COMPSETS_SPEC_FILE")   

    # loop through the components and read in the config_compset.xml file to get all the compsets
    for component in components:
        # Determine the compsets file for this component
        config_file = files.get_value("COMPSETS_SPEC_FILE", {"component":component})

        expect((config_file),
               "Cannot find any config_compsets.xml file for %s" %component)

        # Check that file exists on disk
        if (os.path.isfile(config_file)):
            compsets = Compsets(config_file)
            # TODO finish getting the test_compsets from the compsets module
            help_text, all_compsets, science_compsets, test_compsets = compsets.return_all_values()

            # get the all_compsets sorted
            ordered_compsets = collections.OrderedDict(sorted(all_compsets.items(), key=lambda t: t[0]))

            # load up the html_dict
            html_dict[component] = { 'compsets' : ordered_compsets }

    # loop through the compsets and get the matching the descriptions
##    for compset in ordered_compsets:

    # load up jinja template
    templateLoader = jinja2.FileSystemLoader( searchpath='{0}/templates'.format(work_dir) )
    templateEnv = jinja2.Environment( loader=templateLoader )
    tmplFile = 'compdef2html.tmpl'
    template = templateEnv.get_template( tmplFile )
    templateVars = { 'html_dict'     : html_dict,
                     'all_descs'     : all_descs,
                     'today'         : _now,
                     'model_version' : model_version,
                     'help'          : help_text }
        
    # render the template
    comp_tmpl = template.render( templateVars )

    # write the output file
    with open( options.htmlfile[0], 'w') as html:
        html.write(comp_tmpl)

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




