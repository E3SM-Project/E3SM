"""
API for checking input for testcase
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import get_model

import fnmatch

# Should probably be in XML somewhere
SVN_LOCS = {
    "acme" : "https://acme-svn2.ornl.gov/acme-repo/acme/inputdata",
    "cesm" : "https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata"
}

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

    full_path = os.path.join(input_data_root, rel_path)
    logging.info("Trying to download file: '%s' to path '%s'" % (full_url, full_path))
    # Make sure local path exists, create if it does not
    if(not os.path.exists(os.path.dirname(full_path))):
        os.makedirs(os.path.dirname(full_path))

    stat, out, err = run_cmd("svn --non-interactive --trust-server-cert ls %s" % full_url)
    if (stat != 0):
        logging.warning("FAIL: SVN repo '%s' does not have file '%s'\nReason:%s\n%s\n" % (svn_loc, full_url, out, err))
        return False
    else:
        stat, output, errput = \
            run_cmd("svn --non-interactive --trust-server-cert export %s %s" % (full_url, full_path))
        if (stat != 0):
            logging.warning("svn export failed with output: %s and errput %s\n" % (output, errput))
            return False
        else:
            # Make sure it is group r/w
            os.chmod(full_path, 0664)
            logging.info("SUCCESS\n")
            return True

def check_input_data(case, svn_loc=None, input_data_root=None, data_list_dir="Buildconf", download=False):
    """
    Return True if no files missing
    """
    # Fill in defaults as needed
    svn_loc = SVN_LOCS[get_model()] if svn_loc is None else svn_loc
    input_data_root = case.get_value("DIN_LOC_ROOT") if input_data_root is None else input_data_root

    expect(os.path.isdir(input_data_root), "Invalid input_data_root directory: '%s'" % input_data_root)
    expect(os.path.isdir(data_list_dir), "Invalid data_list_dir directory: '%s'" % data_list_dir)

    data_list_files = find_files(data_list_dir, "*.input_data_list")
    expect(data_list_files, "No .input_data_list files found in dir '%s'" % data_list_dir)

    no_files_missing = True
    for data_list_file in data_list_files:
        logging.info("Loading input file: '%s'" % data_list_file)
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

                    # There are some special values of rel_path that
                    # we need to ignore - some of the component models
                    # set things like 'NULL' or 'same_as_TS' -
                    # basically if rel_path does not contain '/' (a
                    # directory tree) you can assume it's a special
                    # value and ignore it (perhaps with a warning)
                    if ("/" in rel_path and not os.path.exists(full_path)):
                        model = os.path.basename(data_list_file).split('.')[0]
                        logging.warning("Model %s missing file %s = '%s'" % (model,description,full_path))

                        if (download):
                            success = download_if_in_repo(svn_loc, input_data_root, rel_path)
                            if (not success):
                                # If ACME, try CESM repo as backup
                                if (get_model() == "acme" and svn_loc != SVN_LOCS["cesm"]):
                                    success = download_if_in_repo(SVN_LOCS["cesm"], input_data_root, rel_path)
                                    if (not success):
                                        no_files_missing = False
                                else:
                                    no_files_missing = False
                        # if not download
                        else:
                            no_files_missing = False
                    else:
                        logging.debug("Already had input file: '%s'" % full_path)

                else:
                    model = os.path.basename(data_list_file).split('.')[0]
                    logging.warning("Model %s no file specified for %s"%(model,description))

    return no_files_missing
