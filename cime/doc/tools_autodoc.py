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
import re
import stat
import shutil
from string import Template
import traceback

if sys.version_info[0] == 2:
    from ConfigParser import SafeConfigParser as config_parser
else:
    from configparser import ConfigParser as config_parser

# define rst templates
_tool_template = Template('''
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

_script_template = Template('''
.. _$script_name:

####################################################
$script_name 
####################################################

**$script_name** is a script in CIMEROOT/scripts.

.. toctree::
   :maxdepth: 1

.. command-output:: ./$script_name --help
   :cwd: ../../$scripts_dir
''')

_tmpl_template = Template('''
.. _$tmpl_name:

####################################################
$tmpl_name 
####################################################

**$tmpl_name** is a script template in CIMEROOT/config/cesm/machines
that is dynamically created in the CASEROOT. Note that scripts
**.case.run** and **.case.test** are hidden files and should never
be run manually.

.. toctree::
   :maxdepth: 1

.. command-output:: ./$tmpl_name --help
   :cwd: ./temp_files
''')


# -------------------------------------------------------------------------------
# User input
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


# -------------------------------------------------------------------------------
# read the tools_autodoc.cfg configuration file
# -------------------------------------------------------------------------------
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
# create the rst files for the Tools configuration settings
# -------------------------------------------------------------------------------
def get_tools(config, doc_dir):

    # get the input tools dir
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

    tools_dir = config.get('tools','tools_dir')
    for f in tools_files:
        tool_file = os.path.join(doc_dir, '{0}.rst'.format(f))
        with open(tool_file,'w') as tf:
            contents = _tool_template.substitute(tool_name=f, tools_dir=tools_dir)
            tf.write(contents)

    return tools_files


# -------------------------------------------------------------------------------
# create the rst files for the scripts configuration settings
# -------------------------------------------------------------------------------
def get_scripts(config, doc_dir):

    # get the input scripts dir
    scripts_dir = config.get('scripts','scripts_dir')
    scripts_dir = os.path.abspath(scripts_dir)

    # get list of files to exclude
    exclude_files = config.get('scripts','exclude_files').split()

    # get list of files to exclude
    exclude_ext = config.get('scripts','exclude_ext').split()

    # get list of files to exclude
    exclude_prefix = config.get('scripts','exclude_prefix').split()

    # get a list of all files in the scripts_dir
    all_files = next(os.walk(scripts_dir))[2]

    scripts_files = list()
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
            scripts_files.append(f)

    scripts_dir = config.get('scripts','scripts_dir')
    for f in scripts_files:
        script_file = os.path.join(doc_dir, '{0}.rst'.format(f))
        with open(script_file,'w') as tf:
            contents = _script_template.substitute(script_name=f, scripts_dir=scripts_dir)
            tf.write(contents)

    return scripts_files

# -------------------------------------------------------------------------------
# get the template files and substitute the {{...}} patterns so they can be
# run with the --help command
# -------------------------------------------------------------------------------
def get_templates(config, doc_dir):

    # get the input template dir 
    templates_dir = config.get('templates','templates_dir')
    templates_dir = os.path.abspath(templates_dir)

    # get list of files to exclude
    exclude_files = config.get('templates','exclude_files').split()

    # get list of files to exclude
    exclude_ext = config.get('templates','exclude_ext').split()

    # get list of files to exclude
    exclude_prefix = config.get('templates','exclude_prefix').split()

    # get a list of all files in the templates_dir
    all_files = next(os.walk(templates_dir))[2]

    template_files = list()
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
            template_files.append(f)

    # create temporary files with the {{..}} stripped out
    temp_files = list()
    temp_dir = '{0}/temp_files'.format(doc_dir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    for fname in template_files:
        with open(os.path.join(templates_dir,fname),'r') as f:
            content = f.read()
        content = content.replace("{{ batchdirectives }}", "# {{ batchdirective }}", 1)
        content = content.replace("os.chdir( '{{ caseroot }}')", "# os.chdir( '{{ caseroot }}')", 1)
        content = content.replace('os.path.join("{{ cimeroot }}", "scripts", "Tools")', 'os.path.join("../../../..","scripts", "Tools")',1)
        # create a temporary file
        tf = fname.split('.')
        tfname = '.'.join(tf[1:])
        if tfname == 'st_archive':
            tfname = 'case.st_archive'
        tfile = os.path.join(temp_dir, tfname)
        with open(tfile, 'w') as tf:
            tf.write(content)
        
        temp_files.append(tfname)

    for f in temp_files:
        tmpl_file = os.path.join(doc_dir, '{0}.rst'.format(f))
        with open(tmpl_file,'w') as tf:
            contents = _tmpl_template.substitute(tmpl_name=f, temp_dir=temp_dir)
            tf.write(contents)
            tf.close()
        exefile = os.path.join(doc_dir, 'temp_files/{0}'.format(f))
        st = os.stat(exefile)
        os.chmod(exefile, st.st_mode | stat.S_IEXEC)

    return temp_files

# -------------------------------------------------------------------------------
# main
# -------------------------------------------------------------------------------

def main(options):
    all_files = list()
    config = read_config_file(options.config[0])

    # get the output doc dir
    doc_dir = config.get('doc','doc_dir')
    doc_dir = os.path.abspath(doc_dir)

    # gather the files from different locations in the CIMEROOT
    tools_files = get_tools(config, doc_dir)
    scripts_files = get_scripts(config, doc_dir)
    template_files = get_templates(config, doc_dir)

    all_files = tools_files + scripts_files + template_files
    all_files.sort()

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
        for f in all_files:
            index_rst.write('   {0}\n'.format(f))

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
