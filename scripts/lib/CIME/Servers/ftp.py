"""
FTP Server class.  Interact with a server using FTP protocol
"""
# pylint: disable=super-init-not-called
from CIME.XML.standard_module_setup import *
from CIME.Servers.generic_server import GenericServer
from ftplib import FTP as FTPpy

logger = logging.getLogger(__name__)
# I think that multiple inheritence would be useful here, but I couldnt make it work
# in a py2/3 compatible way.
class FTP(GenericServer):
    def __init__(self, address, user='', passwd=''):
        ftp_server, root_address = address.split('/', 1)
        logger.info("server address {} root path {}".format(ftp_server, root_address))
        self.ftp = FTPpy(ftp_server)

        self._ftp_server = address

        stat = self.ftp.login(user, passwd)
        logger.debug("login stat {}".format(stat))
        if "Login successful" not in stat:
            logging.warning("FAIL: Could not login to ftp server {}\n error {}".format(ftp_server, stat))
            return None
        stat = self.ftp.cwd(root_address)
        logger.debug("cwd {} stat {}".format(root_address,stat))
        if "Directory successfully changed" not in stat:
            logging.warning("FAIL: Could not cd to server root directory {}\n error {}".format(root_address, stat))
            return None

    def fileexists(self, rel_path):
        try:
            stat = self.ftp.nlst(rel_path)
        except:
            logger.warning("ERROR from ftp server, trying next server")
            return False
        if rel_path not in stat:
            if not stat or not stat[0].startswith(rel_path):
                logging.warning("FAIL: File {} not found.\nerror {}".format(rel_path, stat))
                return False
        return True

    def getfile(self, rel_path, full_path):
        try:
            stat = self.ftp.retrbinary('RETR {}'.format(rel_path), open(full_path, "wb").write)
        except:
            logger.warning("ERROR from ftp server, trying next server")
            return False

        if (stat != '226 Transfer complete.'):
            logging.warning("FAIL: Failed to retreve file '{}' from FTP repo '{}' stat={}\n".
                            format(rel_path, self._ftp_server, stat))
            return False
        return True

    def getdirectory(self, rel_path, full_path):
        try:
            stat = self.ftp.nlst(rel_path)
        except:
            logger.warning("ERROR from ftp server, trying next server")
            return False

        for _file in stat:
            self.getfile(_file, full_path+os.sep+os.path.basename(_file))
