#!/usr/bin/env python

"""Generator of html file for a given caseroot env_*.xml input file 
"""

# Typically ignore this.
# pylint: disable=invalid-name

# Disable these because this is our standard setup
# pylint: disable=wildcard-import,unused-wildcard-import,wrong-import-position

import os, sys, re, glob
import datetime
import xml.etree.ElementTree as etree

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
from CIME.XML.env_test              import EnvTest
from CIME.XML.env_mach_specific     import EnvMachSpecific
from CIME.XML.env_case              import EnvCase
from CIME.XML.env_mach_pes          import EnvMachPes
from CIME.XML.env_build             import EnvBuild
from CIME.XML.env_run               import EnvRun
from CIME.XML.env_archive           import EnvArchive
from CIME.XML.env_batch             import EnvBatch

# check for  dependency module
try:
    import jinja2
except:
    raise SystemExit("ERROR: compdef2html.py depends on the jinja2 template module. " /
                     "Install using 'pip --user install jinja2'")

# global variables
_now = datetime.datetime.now().strftime('%Y-%m-%d')
_xmlfiles = ['env_archive.xml','env_batch.xml','env_build.xml','env_case.xml',
             'env_mach_pes.xml','env_mach_specific.xml','env_run.xml','env_test.xml']
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

    parser.add_argument('--xmlfile', nargs=1, required=True, choices=_xmlfiles,
                        help='Must be one of {0}.'.format(_xmlfiles))

    parser.add_argument('--caseroot', nargs=1, required=True,
                        help='caseroot location of input xmlfile')

    options = parser.parse_args()

    CIME.utils.parse_args_and_handle_standard_logging_options(options)

    return options


###############################################################################
def parse_archive(env_file, xml_dict):
###############################################################################

    # get an archive object
    archive = EnvArchive(infile=env_file)

    # traverse the archive object and load the xml_dict
    for archive_entry in archive.get_entries():
        compname, compclass = archive.get_entry_info(archive_entry)
        rest_file_extensions = archive.get_rest_file_extensions(archive_entry)
        hist_file_extensions = archive.get_hist_file_extensions(archive_entry)
        rpointer_contents = archive.get_rpointer_contents(archive_entry)
        rest_hist_varname = archive.get_entry_value('rest_history_varname', archive_entry)

        xml_dict[compclass] = { 'compname'             : compname,
                               'rest_file_extensions' : rest_file_extensions,
                               'hist_file_extensions' : hist_file_extensions,
                               'rpointers'            : rpointer_contents,
                               'rest_history_varname' : rest_hist_varname }
    return xml_dict

###############################################################################
def parse_batch(env_file, xml_dict):
###############################################################################

    batch = etree.parse(env_file)
    root = batch.getroot()
    
    # get the header text
    header = root.find('header').text

    # loop throught the batch_system elements
    batch_systems = root.findall('batchsystem')


    batch = EnvBatch(infile=env_file)
    header = batch.get_nodes('header')
    batch_systems = batch.get_nodes('batch_system')
    groups = batch.get_nodes('group')
    for batch_system in batch_systems.

    # load the xml_dict
    return xml_dict

###############################################################################
def parse_build(env_file, xml_dict):
###############################################################################

    return xml_dict


###############################################################################
def parse_case(env_file, xml_dict):
###############################################################################

    return xml_dict


###############################################################################
def parse_mach_pes(env_file, xml_dict):
###############################################################################

    return xml_dict


###############################################################################
def parse_mach_specific(env_file, xml_dict):
###############################################################################

    return xml_dict


###############################################################################
def parse_run(env_file, xml_dict):
###############################################################################

    return xml_dict


###############################################################################
def parse_test(env_file, xml_dict):
###############################################################################

    return xml_dict

###############################################################################
def _main_func(options, work_dir):
###############################################################################

    """Construct html from an XML file."""
        
    # Initialize a variables for the html template
    xml_dict = dict()
    model_version = options.version[0]

    xmlfile = options.xmlfile[0]
    caseroot = options.caseroot[0]
    env_file = os.path.join(caseroot, xmlfile)
    expect(os.path.isfile(env_file),
           "Cannot find env_file {} on disk".format(env_file))

    # call the appropriate file handler routine
    if 'env_archive' in env_file:
        xml_dict = parse_archive(env_file, xml_dict)
        tmplFile = 'env_archive.tmpl'

    elif 'env_batch' in env_file:
        xml_dict = parse_batch(env_file, xml_dict)
        tmplFile = 'env_batch.tmpl'

    elif 'env_build' in env_file:
        xml_dict = parse_build(env_file, xml_dict)
        tmplFile = 'env_build.tmpl'

    elif 'env_case' in env_file:
        xml_dict = parse_case(env_file, xml_dict)
        tmplFile = 'env_case.tmpl'

    elif 'env_mach_pes' in env_file:
        xml_dict = parse_mach_pes(env_file, xml_dict)
        tmplFile = 'env_mach_pes.tmpl'

    elif 'env_mach_specific' in env_file:
        xml_dict = parse_mach_specific(env_file, xml_dict)
        tmplFile = 'env_mach_specific.tmpl'

    elif 'env_run' in env_file:
        xml_dict = parse_run(env_file, xml_dict)
        tmplFile = 'env_run.tmpl'

    elif 'env_test' in env_file:
        xml_dict = parse_test(env_file, xml_dict)
        tmplFile = 'env_test.tmpl'

    else:
        print("Invalid caseroot file {0}. exiting...".format(env_file))
        sys.exit(1)

    # load up jinja template
    templateLoader = jinja2.FileSystemLoader( searchpath='{0}/templates'.format(work_dir) )
    templateEnv = jinja2.Environment( loader=templateLoader )

    # TODO - get the cesm_version for the CIME root
    template = templateEnv.get_template( tmplFile )
    templateVars = { 'xml_dict'     : xml_dict,
                     'today'         : _now,
                     'model_version' : model_version }
        
    # render the template
    xml_tmpl = template.render( templateVars )

    # write the output file
    with open( options.htmlfile[0], 'w') as html:
        html.write(xml_tmpl)

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




