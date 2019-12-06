"""
API for checking input for testcase
"""
from CIME.XML.standard_module_setup import *
from CIME.utils import SharedArea, find_files, safe_copy, expect
from CIME.XML.inputdata import Inputdata
import CIME.Servers

import glob, hashlib, shutil

logger = logging.getLogger(__name__)
# The inputdata_checksum.dat file will be read into this hash if it's available
chksum_hash = dict()
local_chksum_file = 'inputdata_checksum.dat'

def _download_checksum_file(rundir):
    """
    Download the checksum files from each server and merge them into rundir.
    """
    inputdata = Inputdata()
    protocol = "svn"
    # download and merge all available chksum files.
    while protocol is not None:
        protocol, address, user, passwd, chksum_file = inputdata.get_next_server()
        if protocol not in vars(CIME.Servers):
            logger.warning("Client protocol {} not enabled".format(protocol))
            continue
        logger.info("Using protocol {} with user {} and passwd {}".format(protocol, user, passwd))
        if protocol == "svn":
            server = CIME.Servers.SVN(address, user, passwd)
        elif protocol == "gftp":
            server = CIME.Servers.GridFTP(address, user, passwd)
        elif protocol == "ftp":
            server = CIME.Servers.FTP(address, user, passwd)
        elif protocol == "wget":
            server = CIME.Servers.WGET(address, user, passwd)
        else:
            expect(False, "Unsupported inputdata protocol: {}".format(protocol))

        if not chksum_file:
            continue

        success = False
        rel_path = chksum_file
        full_path = os.path.join(rundir, local_chksum_file)
        new_file = full_path + '.raw'
        protocol = type(server).__name__
        logging.info("Trying to download file: '{}' to path '{}' using {} protocol.".format(rel_path, new_file, protocol))
        tmpfile = None
        if os.path.isfile(full_path):
            tmpfile = full_path+".tmp"
            os.rename(full_path, tmpfile)
        # Use umask to make sure files are group read/writable. As long as parent directories
        # have +s, then everything should work.
        with SharedArea():
            success = server.getfile(rel_path, new_file)
            if success:
                _reformat_chksum_file(full_path, new_file)
                if tmpfile:
                    _merge_chksum_files(full_path, tmpfile)
                chksum_hash.clear()
            else:
                if tmpfile and os.path.isfile(tmpfile):
                    os.rename(tmpfile, full_path)
                    logger.warning("Could not automatically download file "+full_path+
                                   " Restoring existing version.")
                else:
                    logger.warning("Could not automatically download file {}".
                                   format(full_path))


def _reformat_chksum_file(chksum_file, server_file):
    """
    The checksum file on the server has 8 space seperated columns, I need the first and last ones.
    This function gets the first and last column of server_file and saves it to chksum_file
    """
    with open(server_file) as fd, open(chksum_file,"w") as fout:
        lines = fd.readlines()
        for line in lines:
            lsplit = line.split()
            if len(lsplit) < 8 or ' DIR ' in line:
                continue

            # remove the first directory ('inputdata/') from the filename
            chksum = lsplit[0]
            fname = (lsplit[7]).split('/',1)[1]
            fout.write(" ".join((chksum, fname))+"\n")
    os.remove(server_file)

def _merge_chksum_files(new_file, old_file):
    """
    If more than one server checksum file is available, this merges the files and removes
    any duplicate lines
    """
    with open(old_file) as fin:
        lines = fin.readlines()
    with open(new_file) as fin:
        lines += fin.readlines()
    lines = set(lines)
    with open(new_file, "w") as fout:
        fout.write("".join(lines))
    os.remove(old_file)



def _download_if_in_repo(server, input_data_root, rel_path, isdirectory=False):
    """
    Return True if successfully downloaded
    server is an object handle of type CIME.Servers
    input_data_root is the local path to inputdata (DIN_LOC_ROOT)
    rel_path is the path to the file or directory relative to input_data_root
    user is the user name of the person running the script
    isdirectory indicates that this is a directory download rather than a single file
    """
    if not (rel_path or server.fileexists(rel_path)):
        return False

    full_path = os.path.join(input_data_root, rel_path)
    logging.info("Trying to download file: '{}' to path '{}' using {} protocol.".format(rel_path, full_path, type(server).__name__))
    # Make sure local path exists, create if it does not
    if isdirectory or full_path.endswith(os.sep):
        if not os.path.exists(full_path):
            logger.info("Creating directory {}".format(full_path))
            os.makedirs(full_path+".tmp")
        isdirectory = True
    elif not os.path.exists(os.path.dirname(full_path)):
        os.makedirs(os.path.dirname(full_path))

    # Use umask to make sure files are group read/writable. As long as parent directories
    # have +s, then everything should work.
    with SharedArea():
        if isdirectory:
            success = server.getdirectory(rel_path, full_path+".tmp")
            # this is intended to prevent a race condition in which
            # one case attempts to use a refdir before another one has
            # completed the download
            if success:
                os.rename(full_path+".tmp",full_path)
            else:
                shutil.rmtree(full_path+".tmp")
        else:
            success = server.getfile(rel_path, full_path)
    return success

def check_all_input_data(self, protocol=None, address=None, input_data_root=None, data_list_dir="Buildconf",
                         download=True, chksum=False):
    """
    Read through all files of the form *.input_data_list in the data_list_dir directory.  These files
    contain a list of input and boundary files needed by each model component.  For each file in the
    list confirm that it is available in input_data_root and if not (optionally download it from a
    server at address using protocol.  Perform a chksum of the downloaded file.
    """
    success = False
    if protocol is not None and address is not None:
        success = self.check_input_data(protocol=protocol, address=address, download=download,
                                        input_data_root=input_data_root, data_list_dir=data_list_dir, chksum=chksum)
    else:
        if chksum:
            _download_checksum_file(self.get_value("RUNDIR"))

        success = self.check_input_data(protocol=protocol, address=address, download=False,
                                        input_data_root=input_data_root, data_list_dir=data_list_dir, chksum=chksum)
        if download and not success:
            if not chksum:
                _download_checksum_file(self.get_value("RUNDIR"))
            success = _downloadfromserver(self, input_data_root, data_list_dir)

    expect(not download or (download and success), "Could not find all inputdata on any server")
    self.stage_refcase(input_data_root=input_data_root, data_list_dir=data_list_dir)
    return success

def _downloadfromserver(case, input_data_root, data_list_dir):
    """
    Download files
    """
    success = False
    protocol = 'svn'
    inputdata = Inputdata()
    if not input_data_root:
        input_data_root = case.get_value('DIN_LOC_ROOT')

    while not success and protocol is not None:
        protocol, address, user, passwd, _ = inputdata.get_next_server()
        logger.info("Checking server {} with protocol {}".format(address, protocol))
        success = case.check_input_data(protocol=protocol, address=address, download=True,
                                        input_data_root=input_data_root,
                                        data_list_dir=data_list_dir,
                                        user=user, passwd=passwd)
    return success

def stage_refcase(self, input_data_root=None, data_list_dir=None):
    """
    Get a REFCASE for a hybrid or branch run
    This is the only case in which we are downloading an entire directory instead of
    a single file at a time.
    """
    get_refcase  = self.get_value("GET_REFCASE")
    run_type     = self.get_value("RUN_TYPE")
    continue_run = self.get_value("CONTINUE_RUN")

    # We do not fully populate the inputdata directory on every
    # machine and do not expect every user to download the 3TB+ of
    # data in our inputdata repository. This code checks for the
    # existence of inputdata in the local inputdata directory and
    # attempts to download data from the server if it's needed and
    # missing.
    if get_refcase and run_type != "startup" and not continue_run:
        din_loc_root = self.get_value("DIN_LOC_ROOT")
        run_refdate  = self.get_value("RUN_REFDATE")
        run_refcase  = self.get_value("RUN_REFCASE")
        run_refdir   = self.get_value("RUN_REFDIR")
        rundir       = self.get_value("RUNDIR")

        if os.path.isabs(run_refdir):
            refdir = run_refdir
            expect(os.path.isdir(refdir), "Reference case directory {} does not exist or is not readable".format(refdir))

        else:
            refdir = os.path.join(din_loc_root, run_refdir, run_refcase, run_refdate)
            if not os.path.isdir(refdir):
                logger.warning("Refcase not found in {}, will attempt to download from inputdata".format(refdir))
                with open(os.path.join("Buildconf","refcase.input_data_list"),"w") as fd:
                    fd.write("refdir = {}{}".format(refdir, os.sep))
                if input_data_root is None:
                    input_data_root = din_loc_root
                if data_list_dir is None:
                    data_list_dir = "Buildconf"
                success = _downloadfromserver(self, input_data_root=input_data_root, data_list_dir=data_list_dir)
                expect(success, "Could not download refcase from any server")

        logger.info(" - Prestaging REFCASE ({}) to {}".format(refdir, rundir))

        # prestage the reference case's files.

        if (not os.path.exists(rundir)):
            logger.debug("Creating run directory: {}".format(rundir))
            os.makedirs(rundir)
        rpointerfile = None
        # copy the refcases' rpointer files to the run directory
        for rpointerfile in glob.iglob(os.path.join("{}","*rpointer*").format(refdir)):
            logger.info("Copy rpointer {}".format(rpointerfile))
            safe_copy(rpointerfile, rundir)
        expect(rpointerfile,"Reference case directory {} does not contain any rpointer files".format(refdir))
        # link everything else

        for rcfile in glob.iglob(os.path.join(refdir,"*")):
            rcbaseline = os.path.basename(rcfile)
            if not os.path.exists("{}/{}".format(rundir, rcbaseline)):
                logger.info("Staging file {}".format(rcfile))
                os.symlink(rcfile, "{}/{}".format(rundir, rcbaseline))
        # Backward compatibility, some old refcases have cam2 in the name
        # link to local cam file.
        for cam2file in  glob.iglob(os.path.join("{}","*.cam2.*").format(rundir)):
            camfile = cam2file.replace("cam2", "cam")
            os.symlink(cam2file, camfile)
    elif not get_refcase and run_type != "startup":
        logger.info("GET_REFCASE is false, the user is expected to stage the refcase to the run directory.")
        if os.path.exists(os.path.join("Buildconf","refcase.input_data_list")):
            os.remove(os.path.join("Buildconf","refcase.input_data_list"))
    return True

def check_input_data(case, protocol="svn", address=None, input_data_root=None, data_list_dir="Buildconf",
                     download=False, user=None, passwd=None, chksum=False):
    """
    For a given case check for the relevant input data as specified in data_list_dir/*.input_data_list
    in the directory input_data_root, if not found optionally download it using the servers specified
    in config_inputdata.xml.  If a chksum file is available compute the chksum and compare it to that
    in the file.
    Return True if no files missing
    """
    case.load_env(reset=True)
    rundir = case.get_value("RUNDIR")
    # Fill in defaults as needed
    input_data_root = case.get_value("DIN_LOC_ROOT") if input_data_root is None else input_data_root

    expect(os.path.isdir(data_list_dir), "Invalid data_list_dir directory: '{}'".format(data_list_dir))

    data_list_files = find_files(data_list_dir, "*.input_data_list")
    expect(data_list_files, "No .input_data_list files found in dir '{}'".format(data_list_dir))

    no_files_missing = True
    if download:
        if protocol not in vars(CIME.Servers):
            logger.warning("Client protocol {} not enabled".format(protocol))
            return False
        logger.info("Using protocol {} with user {} and passwd {}".format(protocol, user, passwd))
        if protocol == "svn":
            server = CIME.Servers.SVN(address, user, passwd)
        elif protocol == "gftp":
            server = CIME.Servers.GridFTP(address, user, passwd)
        elif protocol == "ftp":
            server = CIME.Servers.FTP(address, user, passwd)
        elif protocol == "wget":
            server = CIME.Servers.WGET(address, user, passwd)
        else:
            expect(False, "Unsupported inputdata protocol: {}".format(protocol))

    for data_list_file in data_list_files:
        logging.info("Loading input file list: '{}'".format(data_list_file))
        with open(data_list_file, "r") as fd:
            lines = fd.readlines()

        for line in lines:
            line = line.strip()
            if (line and not line.startswith("#")):
                tokens = line.split('=')
                description, full_path = tokens[0].strip(), tokens[1].strip()
                if description.endswith('datapath'):
                    continue
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
                            logging.warning("Model {} missing file {} = '{}'".format(model, description, full_path))
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
                        isdirectory=rel_path.endswith(os.sep)

                        if ("/" in rel_path and not os.path.exists(full_path)):
                            logger.warning("  Model {} missing file {} = '{}'".format(model, description, full_path))
                            no_files_missing = False

                            if (download):
                                no_files_missing = _download_if_in_repo(server,
                                                                        input_data_root, rel_path.strip(os.sep),
                                                                        isdirectory=isdirectory)
                                if no_files_missing:
                                    verify_chksum(input_data_root, rundir, rel_path.strip(os.sep), isdirectory)
                        else:
                            if chksum:
                                verify_chksum(input_data_root, rundir, rel_path.strip(os.sep), isdirectory)
                                logger.info("Chksum passed for file {}".format(os.path.join(input_data_root,rel_path)))
                            logging.debug("  Already had input file: '{}'".format(full_path))
                else:
                    model = os.path.basename(data_list_file).split('.')[0]
                    logging.warning("Model {} no file specified for {}".format(model, description))

    return no_files_missing

def verify_chksum(input_data_root, rundir, filename, isdirectory):
    """
    For file in filename perform a chksum and compare the result to that stored in
    the local checksumfile, if isdirectory chksum all files in the directory of form *.*
    """
    hashfile = os.path.join(rundir, local_chksum_file)
    if not chksum_hash:
        if not os.path.isfile(hashfile):
            logger.warning("Failed to find or download file {}".format(hashfile))
            return

        with open(hashfile) as fd:
            lines = fd.readlines()
            for line in lines:
                fchksum, fname = line.split()
                if fname in chksum_hash:
                    expect(chksum_hash[fname] == fchksum, " Inconsistent hashes in chksum for file {}".format(fname))
                else:
                    chksum_hash[fname] = fchksum

    if isdirectory:
        filenames = glob.glob(os.path.join(filename,"*.*"))
    else:
        filenames = [filename]
    for fname in filenames:
        if not os.sep in fname:
            continue
        chksum = md5(os.path.join(input_data_root, fname))
        if chksum_hash:
            if not fname in chksum_hash:
                logger.warning("Did not find hash for file {} in chksum file {}".format(filename, hashfile))
            else:
                expect(chksum == chksum_hash[fname],
                       "chksum mismatch for file {} expected {} found {}".
                       format(os.path.join(input_data_root,fname),chksum, chksum_hash[fname]))

def md5(fname):
    """
    performs an md5 sum one chunk at a time to avoid memory issues with large files.
    """
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
