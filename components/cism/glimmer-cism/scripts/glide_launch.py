#! /usr/bin/env python

# Copyright (C) 2005, 2006, 2007, 2009
# Glimmer-CISM contributors - see AUTHORS file for list of contributors
#
# This file is part of Glimmer-CISM.
#
# Glimmer-CISM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or (at
# your option) any later version.
#
# Glimmer-CISM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.
#
# Glimmer-CISM is hosted on BerliOS.de:
# https://developer.berlios.de/projects/glimmer-cism/

# a handy python script which will launch glide and record the execution time of the model

import os,sys,ConfigParser,os.path,optparse

usage = """%prog [options] config_file
launch glide and record execution time of the model. The model binary can be
either set using the environment variable GLIDE_MODEL or is automatically
determined from the config file.

config_file is the name of the model configuration to be used."""


def get_runtype(config):
    """Determine which model to run given configuration.

    First we check for the presence of the environment variable GLIDE_MODEL, if
    that fails we try to determine which binary to run given the config file.

    config: ConfigParser object"""

    try:
        model = os.environ['GLIDE_MODEL']
    except KeyError:
        # simple_glide
        if ('EISMINT-1 fixed margin' in config.sections() 
            or 'EISMINT-1 moving margin' in config.sections()
            or 'EISMINT-2' in config.sections()):
            model = 'simple_glide'
            # eis_glide
        elif ('EIS ELA' in config.sections()):
            model = 'eis_glide'
            # no idea
        else:
            raise KeyError, 'no idea what model I should start'
    return model

def find_prefix(bname):
    """figure out full path to binary bname"""

    print os.path.dirname(bname)

    #if os.path.dirname(bname)!="":
        

if __name__ == '__main__':


    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-m", "--model", help="name of model binary to be launched", metavar="BINARY")
    parser.add_option("-r", "--results", help="name of file where timing info is stored (default: results)", metavar="RESULTS", default="results")
    parser.add_option("-s", "--submit-sge",action="store_true",default=False,help="submit job to Sun Grid Engine")
    parser.add_option("-o", "--submit-options",default="",help="set additional options for cluster submission")
    (options, args) = parser.parse_args()


    if len(args) == 1:
        configname = args[0]
    else:
        parser.error("no configuration file specified")

    config = ConfigParser.ConfigParser()
    config.readfp(open(configname))

    prefix = os.path.abspath(os.path.dirname(sys.argv[0]))
    p = prefix.split(os.sep)[-1]
    sge_script = prefix
    if p == 'bin':
        sge_script = os.path.abspath(os.path.join(prefix,'..','share','glimmer'))
    sge_script = os.path.join(sge_script,'qsub_glide.sh')
    if options.model == None:
        model = os.path.join(prefix,get_runtype(config))
    else:
        model = options.model

    if not os.path.isfile(model):
        sys.stderr.write("Cannot find model executable %s\n"%model)
        sys.exit(0)
    prog = "%s -r %s %s"%(model,options.results,configname)

    if options.submit_sge:
        if not os.path.isfile(sge_script):
            sys.stderr.write("Cannot find model submission script %s"%sge_script)
            sys.exit(0)
        prog = "qsub %s %s %s"%(options.submit_options,sge_script,prog)

#    try:
#        retcode = subprocess.call(prog,shell=True)
#        if retcode < 0:
#            sys.stderr.write("glide model %s was terminated by signal %d\n"%(model,-retcode))
#    except OSError, e:
#        sys.stderr.write("Execution failed: %s\n"%e)
    print os.popen(prog,'r').read()
