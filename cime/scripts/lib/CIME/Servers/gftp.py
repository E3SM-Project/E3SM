"""
GridFTP Server class.  Interact with a server using GridFTP protocol
"""
# pylint: disable=super-init-not-called
from CIME.XML.standard_module_setup import *
from CIME.Servers.generic_server import GenericServer
from CIME.utils import run_cmd

logger = logging.getLogger(__name__)

class GridFTP(GenericServer):
    def __init__(self, address, user='', passwd=''):
        self._root_address = address

    def fileexists(self, rel_path):
        stat,out,err = run_cmd("globus-url-copy -list {}".format(os.path.join(self._root_address, os.path.dirname(rel_path))+os.sep))
        if stat or os.path.basename(rel_path) not in out:
            logging.warning("FAIL: File {} not found.\nstat={} error={}".format(rel_path, stat, err))
            return False
        return True

    def getfile(self, rel_path, full_path):
        stat, _,err = run_cmd("globus-url-copy -v {} file://{}".format(os.path.join(self._root_address, rel_path), full_path))

        if (stat != 0):
            logging.warning("FAIL: GridFTP repo '{}' does not have file '{}' error={}\n".
                            format(self._root_address,rel_path, err))
            return False
        return True

    def getdirectory(self, rel_path, full_path):
        stat, _,err = run_cmd("globus-url-copy -v -r {}{} file://{}{}".format(os.path.join(self._root_address, rel_path), os.sep, full_path, os.sep))

        if (stat != 0):
            logging.warning("FAIL: GridFTP repo '{}' does not have directory '{}' error={}\n".
                            format(self._root_address,rel_path, err))
            return False
        return True
