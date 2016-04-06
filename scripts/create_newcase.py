#!/usr/bin/env python2.7

from Tools.standard_script_setup import *

from CIME.utils         import expect, get_model, run_cmd
from CIME.XML.files     import Files
from CIME.XML.component import Component
from CIME.XML.compsets  import Compsets
from CIME.XML.grids     import Grids
from CIME.XML.machines  import Machines
from CIME.case          import Case

from CIME.XML.env_case  import EnvCase  # temporary
from CIME.XML.env_build import EnvBuild # temporary
from CIME.XML.env_run   import EnvRun   # temporary
from CIME.XML.env_batch import EnvBatch # temporary
from CIME.XML.env_test  import EnvTest  # temporary

from shutil import copyfile

import os
import stat
import argparse, doctest
import pdb

logger = logging.getLogger(__name__)

###############################################################################
def parse_command_line(args):
###############################################################################

    cime_model = CIME.utils.get_model()

    parser = argparse.ArgumentParser()

    CIME.utils.setup_standard_logging_options(parser)

    parser.add_argument("-case", "--case", required=True,
                        help="(required) Specify the case name. "
                        "If not a full pathname, then the current working directory is assumed")

    parser.add_argument("-compset", "--compset", required=True,
                        help="(required) Specify a compset. "
                        "To see list of current compsets, use the utility manage_case in this directory")

    parser.add_argument("-res", "--resolution", required=True,
                        help="(required) Specify a model grid resolution. "
                        "To see list of current compsets, use the utility manage_case in this directory")

    parser.add_argument("-mach", "--machine", required=True,
                        help="(required) Specify a machine. "
                        "To see list of current compsets, use the utility manage_case in this directory")

    parser.add_argument("-compiler", "--compiler", 
                        help="Specify a compiler. "
                        "To see list of supported compilers for each machine, use the utility manage_case in this directory")

    parser.add_argument("-mpilib", "--mpilib", 
                        help="Specify the mpilib. "
                        "To see list of supported mpilibs for each machine, use the utility manage_case in this directory. "
                        "The default is mpi-serial, but will be replaced by default mpi library for the target machine.")
                        
    parser.add_argument("-project", "--project", type=int, default="00000",
                        help="Specify a project id")

    parser.add_argument("-pes_file", "--pes_file", 
                        help="Specify full pathname of pes file to use."
                        "This file will contain settings that will overwrite the default settings")

    parser.add_argument("-pecount", "--pecount", choices=('S','M','L','X1','X2'), default="M",
                        help="Specify a target size description for the number of cores")

    parser.add_argument("-petype", "--petype", choices=('threaded','mpionly','unset'), default="unset",
                        help="Force pes to be all threaded or all mpi")

    parser.add_argument("-mach_dir", "--mach-dir", 
                        help="Specify the locations of the Machines directory, other than the default"
                        "The default is CIMEROOT/machines")

    parser.add_argument("-user_mods_dir", "--user-mods-dir", 
                        help="Path to directory with user_nl_* files and xmlchange "
                        "commands to utilize. This can also include SourceMods")  

    parser.add_argument("-user_compset", "--user-compset", 
                        help="Long name for new user compset to use" 
                        "This assumes that all of the compset settings in the"
                        "long name have been defined by the target model components")

    parser.add_argument("-user_grid_file", "--user-grid-file", 
                        help="Full pathname of grid file to use"
                        "This should be a copy of cime_config/config_grids.xml"
                        "with the new user grid changes added to it")

    parser.add_argument("-user_pes_setby", "--user-pes-setby", 
                        help="Only used and required for --user-compset argument."
                        "Component that sets the pe-layouts and pio settings"
                        "For CESM this is [allactive, cam, clm, cice, cism, pop")

    parser.add_argument("-srcroot", "--srcroot", 
                        help="Alternative path for source root directory. By defualt this is set to"
                        "cimeroot/../")

    # -sharedlibroot       Used for re-using build components when building multiple cases, default is \$EXEROOT

    args = parser.parse_args()

    CIME.utils.handle_standard_logging_options(args)

    if args.srcroot is not None:
        expect(os.path.isdir(args.srcroot),
               "Input non-default directory srcroot %s does not exist " %args.srcroot)
        args.srcroot = os.path.abspath(args.srcroot)

    if args.user_grid_file is not None:
        expect(os.path.isfile(args.user_grid_file),
               "User_grid_file %s does not exist " %args.user_grid_file)


    return args.case, args.compset, args.resolution, args.machine, args.compiler,\
        args.mpilib, args.project, args.pes_file, args.pecount, args.petype, \
        args.mach_dir, args.user_mods_dir, args.user_compset, args.user_pes_setby, \
        args.user_grid_file, args.srcroot
    
    

###############################################################################
def _create_caseroot_tools(case, caseroot, cimeroot, machine):
###############################################################################

    cime_model = CIME.utils.get_model()
    machobj = Machines(machine=machine)
    machines_dir = machobj.get_machines_dir()

    # setup executable files in caseroot/
    exefiles = (cimeroot + "/scripts/Tools/case.setup",  
                cimeroot + "/scripts/Tools/case.build",
                cimeroot + "/scripts/Tools/case.submit",
                cimeroot + "/scripts/Tools/preview_namelists", 
                cimeroot + "/scripts/Tools/testcase.setup", 
                cimeroot + "/scripts/Tools/check_input_data",
                cimeroot + "/scripts/Tools/check_case", 
                cimeroot + "/scripts/Tools/archive_metadata.sh", 
                cimeroot + "/scripts/Tools/create_production_test",
                cimeroot + "/scripts/Tools/xmlchange", 
                cimeroot + "/scripts/Tools/xmlquery")
    try:
        for exefile in exefiles:
            destfile = caseroot + "/" + os.path.basename(exefile)
            os.symlink(exefile, destfile)
    except Exception as e:
        logger.warning("FAILED to set up exefiles: %s" % str(e))

    # set up utility files in caseroot/Tools/
    toolfiles = (cimeroot + "/cime_config/" + cime_model + "/archive.xml",
                 cimeroot + "/scripts/Tools/check_lockedfiles", 
                 cimeroot + "/scripts/Tools/lt_archive.sh", 
                 cimeroot + "/scripts/Tools/st_archive", 
                 cimeroot + "/scripts/Tools/getTiming", 
                 cimeroot + "/scripts/Tools/compare_namelists.pl",
                 machines_dir + "/taskmaker.pl", 
                 machines_dir + "/Makefile",
                 machines_dir + "/mkSrcfiles", 
                 machines_dir + "/mkDepends") 
    try:
        for toolfile in toolfiles:
            destfile = caseroot + "/Tools/" + os.path.basename(toolfile)
            os.symlink(toolfile, destfile)
    except Exception as e:
        logger.warning("FAILED to set up toolfiles: %s %s %s" % (str(e), toolfile, destfile))

    # set up infon files
    infofiles = (cimeroot + "/scripts/Tools/README.post_process")
    #FIXME - the following does not work
    # print "DEBUG: infofiles are ",infofiles
    #    try:
    #        for infofile in infofiles:
    #            print "DEBUG: infofile is %s, %s"  %(infofile, os.path.basename(infofile))
    #            dst_file = caseroot + "/" + os.path.basename(infofile)
    #            copyfile(infofile, dst_file)
    #            os.chmod(dst_file, os.stat(dst_file).st_mode | stat.S_IXUSR | stat.S_IXGRP)
    #    except Exception as e:
    #        logger.warning("FAILED to set up infofiles: %s" % str(e))

            
###############################################################################
def _create_caseroot_sourcemods(case, caseroot, cimeroot, components):
###############################################################################

    for component in components:
        directory = caseroot + "/SourceMods/src." + component
        if not os.path.exists(directory):
            os.makedirs(directory)

    directory = caseroot + "/SourceMods/src.share"
    if not os.path.exists(directory):
        os.makedirs(directory)

    directory = caseroot + "/SourceMods/src.drv"
    if not os.path.exists(directory):
        os.makedirs(directory)

    cime_model = CIME.utils.get_model()
    if cime_model is "cesm":
        # Note: this is CESM specific, given that we are referencing cism explitly
        if "cism" in components:
            directory = caseroot + "/SourceMods/src.cism/glimmer-cism"
            if not os.path.exists(directory):
                os.makedirs(directory)
                readme_file = os.path.join(directory, "README")
                
                str_to_write = """
                Put source mods for the glimmer-cism library in the glimmer-cism subdirectory
                This includes any files that are in the glimmer-cism subdirectory of $cimeroot/../components/cism
                Anything else (e.g., mods to source_glc or drivers) goes in this directory, NOT in glimmer-cism/"""
                
                with open(readme_file, "w") as fd:
                    fd.write(str_to_write)


###############################################################################
def _create_caseroot(case, caseroot, cimeroot, components, machine):
###############################################################################

    #logger.info(" Creating Caseroot %s" %caseroot)
    logger.info(" Creating Caseroot")

    # Make the case directory
    os.mkdir(caseroot)
    os.chdir(caseroot)

    # Create relevant directories in $caseroot
    newdirs = ("SourceMods", "LockedFiles", "Buildconf", "Tools")
    for newdir in newdirs:
        os.mkdir(newdir)

    # Open a new README.case file in $caseroot
    #file = "$caseroot/README.case";
    #fh = IO::File->new($file, '>' ) or $logger->logdie( "can't open file: $file\n");
    #print $fh "$commandline\n\n\n";
    #fh->close;

    _create_caseroot_sourcemods(case, caseroot, cimeroot, components);  
    _create_caseroot_tools(case, caseroot, cimeroot, machine);

###############################################################################
def _main_func():
###############################################################################

    #pdb.set_trace()
    
    case, compset, grid, machine, \
        compiler, mpilib, project, pes_file, pecount, petype, \
        machine_dir, usermods_dir, user_compset, user_pes_setby, user_grid_file, srcroot \
        = parse_command_line(sys.argv)

    cimeroot  = CIME.utils.get_cime_root()
    if srcroot is None:
        srcroot  = os.path.join(cimeroot,"/../")

    # Determine 
    if os.path.exists(case):
        caseroot = os.path.dirname(case)
        case = os.path.basename(case)
    else:
        caseroot = os.path.join(os.getcwd(),case)

    # Set the case object
    caseobj = Case(caseroot)
    
    # Set values for env_case.xml 
    caseobj.set_value("CASE", case)
    caseobj.set_value("CASEROOT", caseroot)
    caseobj.set_value("SRCROOT", srcroot)
    caseobj.set_value("MACH", machine)
    if user_grid_file is not None:
        caseobj.set_value("GRIDS_SPEC_FILE", user_grid_file);
    
    # Configure the Case
    caseobj.configure(compset, grid, machine)
    
    components = caseobj.get_compset_components()
    _create_caseroot(case, caseroot, cimeroot, components, machine)
    
    # Write out the case files
    for file in caseobj._env_entryid_files:
        file.write()
    caseobj._env_generic_files[0].write()


###############################################################################

if __name__ == "__main__":
    _main_func()
