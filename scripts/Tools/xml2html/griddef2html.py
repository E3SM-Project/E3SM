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

from standard_script_setup import *
from CIME.utils import expect
from CIME.XML.entry_id import GenericXML

# check for  dependency module
try:
    import jinja2
except:
    raise SystemExit("ERROR: nmldef2html.py depends on the jinja2 template module. " /
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
        description='Read the component namelist file and generate a corresponding HTML file.')

    CIME.utils.setup_standard_logging_options(parser)

    parser.add_argument('--nmlfile', nargs=1, required=True,
                        help='Fully nquailfied path to input namelist XML file.')

    parser.add_argument('--comp', nargs=1, required=True, choices=_comps, 
                        help='Component name.')

    parser.add_argument('--htmlfile', nargs=1, required=True,
                        help='Fully quailfied path to output HTML file.')

    options = parser.parse_args()

    CIME.utils.handle_standard_logging_options(options)

    return options

###############################################################################
def _main_func(options, nmldoc_dir):
###############################################################################

    """Construct a `NamelistDefinition` from an XML file."""

    # Create a definition object from the xml file
    filename = options.nmlfile[0]
    expect(os.path.isfile(filename), "File %s does not exist"%filename)
    try:
        definition = GenericXML(infile=filename)
    except:
        sys.exit("Error: unable to parse file %s" %filename)
        
    # Determine if have new or old schema
    basepath = os.path.dirname(filename)
    default_files = glob.glob(os.path.join(basepath,"namelist_defaults*.xml"))
    defaults = []
    if len(default_files) > 0:
        schema = "old"
        for default_file in default_files:
            defaults.append(GenericXML(infile=default_file))
    else:
        schema = "new"

    # Initialize a variables for the html template
    html_dict = dict()
    cesm_version = 'CESM2.0'
    comp = ''
    if options.comp:
        comp = options.comp[0]

    # Create a dictionary with a category key and a list of all entry nodes for each key
    category_dict = dict()
    for node in definition.get_nodes("entry"):
        if schema == "new":
            category = definition.get_element_text("category", root=node)
        else:
            category = node.get("category")
        if category in category_dict:
            category_dict[category].append(node)
        else:
            category_dict[category] = [ node ]

    # Loop over each category and load up the html_dict
    for category in category_dict:

        # Create a dictionary of groups with a group key and an array of group nodes for each key
        groups_dict = dict()
        for node in category_dict[category]:
            if schema == "new":
                group = definition.get_element_text("group", root=node)
            else:
                group = node.get("group") 
            if group in groups_dict:
                groups_dict[group].append(node) 
            else:
                groups_dict[group] = [ node ]

        # Loop over the keys
        group_list = list()
        for group_name in groups_dict:

            # Loop over the nodes in each group
            for node in groups_dict[group_name]:

                # Determine the name
                # @ is used in a namelist to put the same namelist variable in multiple groups
                # in the write phase, all characters in the namelist variable name after 
                # the @ and including the @ should be removed
                name = node.get("id")
                if "@" in name:
                    name = re.sub('@.+$', "", name)

                # Create the information for this node - start with the description
                if schema == "new":
                    raw_desc = definition.get_element_text("desc", root=node)
                else:
                    raw_desc = node.text 

                if raw_desc is not None: 
                    desc = ' '.join(raw_desc.split())
                    short_desc = desc[:75] + (desc[75:] and '...')
                else:
                    desc = ''

                # add type
                if schema == "new":
                    entry_type = definition.get_element_text("type", root=node)
                else:
                    entry_type = node.get("type")

                # add valid_values
                if schema == "new":
                    valid_values = definition.get_element_text("valid_values", root=node)
                else:
                    valid_values = node.get("valid_values")
                    
                if entry_type == "logical":
                    valid_values = ".true.,.false"
                else:
                    if not valid_values:
                        valid_values = "any " + entry_type
                        if "char" in valid_values:
                            valid_values = "any char"

                if valid_values is not None:
                    valid_values = valid_values.split(',')
                    
                # add default values
                values = ""
                if schema == "new":
                    value_nodes = definition.get_nodes('value', root=node)
                    if value_nodes is not None and len(value_nodes) > 0:
                        for value_node in value_nodes:
                            try:
                                value = value_node.text.strip()
                            except:
                                value = 'undefined'
                            if value_node.attrib:
                                values += "value is %s for: %s <br/>" %(value, value_node.attrib)
                            else:
                                values += "value: %s <br/>" %(value)
                else:
                    for default in defaults:
                        value_nodes = default.get_nodes(name)
                        if len(value_nodes) > 0:
                            for value_node in value_nodes:
                                if value_node.attrib:
                                    values += "value is %s for: %s <br/>" %(value_node.text, value_node.attrib)
                                else:
                                    values += "value: %s <br/>" %(value_node.text)

                # create the node dictionary
                node_dict = { 'name'        : name,
                              'desc'        : desc,
                              'short_desc'  : short_desc,
                              'entry_type'  : entry_type,
                              'valid_values': valid_values,
                              'value'       : values,
                              'group_name'  : group_name }

                # append this node_dict to the group_list
                group_list.append(node_dict)

            # update the group_list for this category in the html_dict
            category_group = category
            html_dict[category_group] = group_list

    # load up jinja template
    templateLoader = jinja2.FileSystemLoader( searchpath='{0}/templates'.format(nmldoc_dir) )
    templateEnv = jinja2.Environment( loader=templateLoader )

    # TODO - get the cesm_version for the CIME root
    tmplFile = 'nmldef2html.tmpl'
    template = templateEnv.get_template( tmplFile )
    templateVars = { 'html_dict'    : html_dict,
                     'today'        : _now,
                     'cesm_version' : cesm_version,
                     'comp'         : comp }
        
    # render the template
    nml_tmpl = template.render( templateVars )

    # write the output file
    with open( options.htmlfile[0], 'w') as html:
        html.write(nml_tmpl)

    return 0

###############################################################################

if __name__ == "__main__":

    options = commandline_options()
    nmldoc_dir = os.path.join(CIMEROOT,"scripts","Tools","nmldoc")
    try:
        status = _main_func(options, nmldoc_dir)
        sys.exit(status)
    except Exception as error:
        print(str(error))
        sys.exit(1)




