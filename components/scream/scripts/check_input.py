
import sys, os

# Add CIME libs to sys path
_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..","..","cime")
_LIB_DIR = os.path.join(_CIMEROOT, "scripts", "lib")
sys.path.append(_LIB_DIR)

from CIME.case.check_input_data import _download_if_in_repo
from CIME.utils import safe_copy, expect
from CIME.Servers import WGET, SVN, FTP
from CIME.XML.inputdata import Inputdata

###############################################################################
def download_file(input_root, the_file):
###############################################################################
    inputdata = Inputdata()
    protocol = "wget"
    success = False
    while not success and protocol is not None:
        protocol, address, user, passwd, _, ic_filepath, _ = inputdata.get_next_server()
        if protocol is not None:
            if protocol == "svn":
                server = SVN(address, user, passwd)
            elif protocol == "gftp":
                server = GridFTP(address, user, passwd)
            elif protocol == "ftp":
                server = FTP.ftp_login(address, user, passwd)
            elif protocol == "wget":
                server = WGET.wget_login(address, user, passwd)
            else:
                expect(False, "Unsupported inputdata protocol: {}".format(protocol))

            print("  Attempting to download {} from {} with protocol {}".format(the_file, address, protocol))
            success = _download_if_in_repo(server,
                                           input_root, the_file.strip(os.sep))
            print("  {}".format("SUCCESS" if success else "FAILED"))

    return success

###############################################################################
def check_input(input_root, tgt_dir, files):
###############################################################################
    any_fails = False

    for the_file in files:
        print("Checking for file {} within {} needed in {}.".format(the_file, input_root, tgt_dir))
        basename = os.path.basename(the_file)
        full_tgt = os.path.join(tgt_dir, basename)
        if not os.path.exists(full_tgt):
            full_path = os.path.join(input_root, the_file)
            if not os.path.exists(full_path):
                print("  Input file {} needs to be downloaded.".format(full_path))
                success = download_file(input_root, the_file)
                if success:
                    print("  Copying {} to {}.".format(full_path, tgt_dir))
                    safe_copy(full_path, tgt_dir)
                else:
                    print("  Could not download file {}".format(the_file))
                    any_fails = True
            else:
                print("  Input file {} already downloaded, copying to {}.".format(full_path, tgt_dir))
                safe_copy(full_path, tgt_dir)

        else:
            print("  Input file {} already exists.".format(full_tgt))

    return not any_fails
