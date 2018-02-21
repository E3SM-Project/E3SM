"""
API for checking input for testcase
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import get_model, SharedArea
from CIME.XML.inputdata import Inputdata

import fnmatch, glob, shutil

logger = logging.getLogger(__name__)

def find_files(rootdir, pattern):
    """
    recursively find all files matching a pattern
    """
    result = []
    for root, _, files in os.walk(rootdir):
        for filename in files:
            if (fnmatch.fnmatch(filename, pattern)):
                result.append(os.path.join(root, filename))

    return result

def download_if_in_repo(svn_loc, input_data_root, rel_path):
    """
    Return True if successfully downloaded
    """
    rel_path = rel_path.strip('/')
    full_url = os.path.join(svn_loc, rel_path)
    logger.info("full url {}".format(full_url))
    err = run_cmd("svn --non-interactive --trust-server-cert ls {}".format(svn_loc))[0]
    if err != 0:
        logging.warning(
"""
Could not connect to svn repo '{0}'
This is most likely either a credential, proxy, or network issue .
To check connection and store your credential run 'svn ls {0}' and permanently store your password""".format(svn_loc))
        return False

    full_path = os.path.join(input_data_root, rel_path)
    logging.info("Trying to download file: '{}' to path '{}'".format(full_url, full_path))
    # Make sure local path exists, create if it does not
    if(not os.path.exists(os.path.dirname(full_path))):
        os.makedirs(os.path.dirname(full_path))

    stat, out, err = run_cmd("svn --non-interactive --trust-server-cert ls {}".format(full_url))
    if (stat != 0):
        logging.warning("FAIL: SVN repo '{}' does not have file '{}'\nReason:{}\n{}\n".format(svn_loc, full_url, out, err))
        return False
    else:
        # Use umask to make sure files are group read/writable. As long as parent directories
        # have +s, then everything should work.
        with SharedArea():
            stat, output, errput = \
                run_cmd("svn --non-interactive --trust-server-cert export {} {}".format(full_url, full_path))
            if (stat != 0):
                logging.warning("svn export failed with output: {} and errput {}\n".format(output, errput))
                return False
            else:
                logging.info("SUCCESS\n")
                return True

###############################################################################
def check_all_input_data(case, input_data_root=None, data_list_dir="Buildconf", download=True):
###############################################################################
    success = False
    inputdata = Inputdata()

    while not success:
        protocal, address = inputdata.get_next_server()
        expect(protocal is not None, "Failed to find input data")
        logger.info("Checking server {} with protocal {}".format(address, protocal))
        if protocal == 'svn':
            success = check_input_data(case=case, svn_loc=address, download=download,
                                       input_data_root=input_data_root, data_list_dir=data_list_dir)
        else:
            expect(False,"Unknown inputdata server protocal {}".format(protocal))

    stage_refcase(case)

def stage_refcase(case):
    get_refcase  = case.get_value("GET_REFCASE")
    run_type     = case.get_value("RUN_TYPE")
    continue_run = case.get_value("CONTINUE_RUN")

    # We do not fully populate the inputdata directory on every
    # machine and do not expect every user to download the 3TB+ of
    # data in our inputdata repository. This code checks for the
    # existence of inputdata in the local inputdata directory and
    # attempts to download data from the server if it's needed and
    # missing.
    if get_refcase and run_type != "startup" and not continue_run:
        din_loc_root = case.get_value("DIN_LOC_ROOT")
        run_refdate  = case.get_value("RUN_REFDATE")
        run_refcase  = case.get_value("RUN_REFCASE")
        run_refdir   = case.get_value("RUN_REFDIR")
        rundir       = case.get_value("RUNDIR")

        refdir = os.path.join(din_loc_root, run_refdir, run_refcase, run_refdate)
        expect(os.path.isdir(refdir),
"""
*****************************************************************
prestage ERROR: $refdir is not on local disk
obtain this data from the svn input data repository
> mkdir -p {}
> cd {}
> cd ..
> svn export --force https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/{}
or set GET_REFCASE to FALSE in env_run.xml
and prestage the restart data to $RUNDIR manually
*****************************************************************
""".format(refdir, refdir, refdir))

        logger.info(" - Prestaging REFCASE ({}) to {}".format(refdir, rundir))

        # prestage the reference case's files.

        if (not os.path.exists(rundir)):
            logger.debug("Creating run directory: {}".format(rundir))
            os.makedirs(rundir)

        # copy the refcases' rpointer files to the run directory
        for rpointerfile in glob.iglob(os.path.join("{}","*rpointer*").format(refdir)):
            logger.info("Copy rpointer {}".format(rpointerfile))
            shutil.copy(rpointerfile, rundir)

        # link everything else

        for rcfile in glob.iglob(os.path.join(refdir,"*")):
            rcbaseline = os.path.basename(rcfile)
            if not os.path.exists("{}/{}".format(rundir, rcbaseline)):
                logger.info("Staging file {}".format(rcfile))
                os.symlink(rcfile, "{}/{}".format(rundir, rcbaseline))

        for cam2file in  glob.iglob(os.path.join("{}","*.cam2.*").format(rundir)):
            camfile = cam2file.replace("cam2", "cam")
            os.symlink(cam2file, camfile)

def check_input_data(case, svn_loc=None, input_data_root=None, data_list_dir="Buildconf", download=False):
    """
    Return True if no files missing
    """
    # Fill in defaults as needed

    input_data_root = case.get_value("DIN_LOC_ROOT") if input_data_root is None else input_data_root

    expect(os.path.isdir(input_data_root), "Invalid input_data_root directory: '{}'".format(input_data_root))
    expect(os.path.isdir(data_list_dir), "Invalid data_list_dir directory: '{}'".format(data_list_dir))

    data_list_files = find_files(data_list_dir, "*.input_data_list")
    expect(data_list_files, "No .input_data_list files found in dir '{}'".format(data_list_dir))

    no_files_missing = True
    for data_list_file in data_list_files:
        logging.info("Loading input file list: '{}'".format(data_list_file))
        with open(data_list_file, "r") as fd:
            lines = fd.readlines()

        for line in lines:
            line = line.strip()
            if (line and not line.startswith("#")):
                tokens = line.split('=')
                description, full_path = tokens[0].strip(), tokens[1].strip()
                if(full_path):
                    # expand xml variables
                    full_path = case.get_resolved_value(full_path)
                    rel_path  = full_path.replace(input_data_root, "")
                    model = os.path.basename(data_list_file).split('.')[0]

                    if ("/" in rel_path and rel_path == full_path):
                        # User pointing to a file outside of input_data_root, we cannot determine
                        # rel_path, and so cannot download the file. If it already exists, we can
                        # proceed
                        if not os.path.exists(full_path):
                            logging.warning("  Model {} missing file {} = '{}'".format(model, description, full_path))
                            if download:
                                logging.warning("    Cannot download file since it lives outside of the input_data_root '{}'".format(input_data_root))
                            no_files_missing = False
                        else:
                            logging.debug("  Found input file: '{}'".format(full_path))

                    else:
                        # There are some special values of rel_path that
                        # we need to ignore - some of the component models
                        # set things like 'NULL' or 'same_as_TS' -
                        # basically if rel_path does not contain '/' (a
                        # directory tree) you can assume it's a special
                        # value and ignore it (perhaps with a warning)
                        if ("/" in rel_path and not os.path.exists(full_path)):
                            logging.warning("  Model {} missing file {} = '{}'".format(model, description, full_path))

                            if (download):
                                success = download_if_in_repo(svn_loc, input_data_root, rel_path)
                            # if not download
                            else:
                                no_files_missing = False
                        else:
                            logging.debug("  Already had input file: '{}'".format(full_path))

                else:
                    model = os.path.basename(data_list_file).split('.')[0]
                    logging.warning("Model {} no file specified for {}".format(model, description))

    return no_files_missing
