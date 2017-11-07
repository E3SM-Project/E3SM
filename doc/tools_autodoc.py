#!/usr/bin/env python
"""script to auto generate rst documentation for cime/scripts/Tools 
   user facing utilities.

"""

from __future__ import print_function

import sys

if sys.hexversion < 0x02070000:
    print(70 * "*")
    print("ERROR: {0} requires python >= 2.7.x. ".format(sys.argv[0]))
    print("It appears that you are running python {0}".format(
        ".".join(str(x) for x in sys.version_info[0:3])))
    print(70 * "*")
    sys.exit(1)

import argparse
import os
import traceback
import shutil
from string import Template

if sys.version_info[0] == 2:
    from ConfigParser import SafeConfigParser as config_parser
else:
    from configparser import ConfigParser as config_parser

# -------------------------------------------------------------------------------
#
# User input
#
# -------------------------------------------------------------------------------

def commandline_options():
    """Process the command line arguments.

    """
    parser = argparse.ArgumentParser(
        description='Auto generate rst documentation for cime/scripts/Tools.')

    parser.add_argument('--backtrace', action='store_true',
                        help='show exception backtraces as extra debugging '
                        'output')

    parser.add_argument('--debug', action='store_true',
                        help='extra debugging output')

    parser.add_argument('--config', nargs=1, default=['tools_autodoc.cfg'],
                        help='path to config file')

    options = parser.parse_args()
    return options



def read_config_file(filename):
    """Read the configuration file and process

    """
    print("tools_autodoc.py - Reading configuration file : {0}".format(filename))

    cfg_file = os.path.abspath(filename)
    if not os.path.isfile(cfg_file):
        raise RuntimeError("Could not find config file: {0}".format(cfg_file))

    config = config_parser()
    config.read(cfg_file)

    return config


# -------------------------------------------------------------------------------
#
# main
#
# -------------------------------------------------------------------------------

def main(options):
    config = read_config_file(options.config[0])

    tools_dir = config.get('tools','tools_dir')
    tools_dir = os.path.abspath(tools_dir)

    # get list of files to exclude
    exclude_files = config.get('tools','exclude_files').split()

    # get list of files to exclude
    exclude_ext = config.get('tools','exclude_ext').split()

    # get list of files to exclude
    exclude_prefix = config.get('tools','exclude_prefix').split()

    # get a list of all files in the tools_dir
    all_files = next(os.walk(tools_dir))[2]

    tools_files = list()
    # exclude files 
    for f in all_files:
        f = f.strip()
        include = True
        for e in exclude_files:
            if f == e.strip():
                include = False
        for e in exclude_ext:
            if f.endswith(e.strip()):
                include = False
        for e in exclude_prefix:
            if f.startswith(e.strip()):
                include = False
        if include:
            tools_files.append(f)

    # copy the index.rst.template to index.rst
    doc_dir = config.get('doc','doc_dir')
    doc_dir = os.path.abspath(doc_dir)

    index_template = config.get('doc','index_template')
    index_rst_file = index_template.split('.')[0:-1]
    index_template = os.path.join(doc_dir,index_template)
    index_rst_file = '.'.join(index_rst_file)
    index_rst_file = os.path.join(doc_dir,index_rst_file)

    shutil.copy2(index_template, index_rst_file)

    # open index_rst_file in append mode
    with open(index_rst_file,'a') as index_rst:
        for f in tools_files:
            index_rst.write('   {0}\n'.format(f))

    tool_template = Template('''
.. _$tool_name:

####################################################
$tool_name 
####################################################

**$tool_name** is a script in CIMEROOT/scripts/Tools.

.. toctree::
   :maxdepth: 1

.. command-output:: ./$tool_name --help
   :cwd: ../../$tools_dir
''')

    tools_dir = config.get('tools','tools_dir')
    for f in tools_files:
        tool_file = os.path.join(doc_dir, '{0}.rst'.format(f))
        with open(tool_file,'w') as tf:
            contents = tool_template.substitute(tool_name=f, tools_dir=tools_dir)
            tf.write(contents)

    return 0


if __name__ == "__main__":
    options = commandline_options()
    try:
        status = main(options)
        sys.exit(status)
    except Exception as error:
        print(str(error))
        if options.backtrace:
            traceback.print_exc()
        sys.exit(1)
