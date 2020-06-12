"""
FTP Server class.  Interact with a server using FTP protocol
"""
# pylint: disable=super-init-not-called
from CIME.XML.standard_module_setup import *
from CIME.Servers.generic_server import GenericServer
from CIME.utils import Timeout
from ftplib import FTP as FTPpy
from ftplib import all_errors as all_ftp_errors
import socket

logger = logging.getLogger(__name__)
# I think that multiple inheritence would be useful here, but I couldnt make it work
# in a py2/3 compatible way.
class FTP(GenericServer):
    def __init__(self, address, user='', passwd='', server=None):
        if not user:
            user = ''
        if not passwd:
            passwd = ''
        expect(server," Must call via ftp_login function")
        root_address = address.split('/', 1)[1]
        self.ftp = server
        self._ftp_server = address
        stat = self.ftp.login(user, passwd)
        logger.debug("login stat {}".format(stat))
        if "Login successful" not in stat:
            logging.warning("FAIL: Could not login to ftp server {}\n error {}".format(address, stat))
            return None

        stat = self.ftp.cwd(root_address)
        logger.debug("cwd {} stat {}".format(root_address,stat))
        if "Directory successfully changed" not in stat:
            logging.warning("FAIL: Could not cd to server root directory {}\n error {}".format(root_address, stat))
            return None

    @classmethod
    def ftp_login(cls, address, user='', passwd=''):
        ftp_server, root_address = address.split('/', 1)
        logger.info("server address {} root path {}".format(ftp_server, root_address))
        try:
            with Timeout(60):
                ftp = FTPpy(ftp_server)

        except socket.error as e:
            logger.warning("ftp login timeout! {} ".format(e))
            return None
        except RuntimeError:
            logger.warning("ftp login timeout!")
            return None
        return cls(address, user=user, passwd=passwd, server=ftp)

    def fileexists(self, rel_path):
        try:
            stat = self.ftp.nlst(rel_path)
        except all_ftp_errors:
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
        except all_ftp_errors:
            if os.path.isfile(full_path):
                os.remove(full_path)
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
        except all_ftp_errors:
            logger.warning("ERROR from ftp server, trying next server")
            return False

        for _file in stat:
            self.getfile(_file, full_path+os.sep+os.path.basename(_file))
