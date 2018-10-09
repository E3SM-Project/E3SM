"""
SVN Server class.  Interact with a server using SVN protocol
"""
# pylint: disable=super-init-not-called
from CIME.XML.standard_module_setup import *
from CIME.Servers.generic_server import GenericServer

logger = logging.getLogger(__name__)

class SVN(GenericServer):
    def __init__(self, address, user='', passwd=''):
        self._args = ''
        if user:
            self._args += "--username {}".format(user)
        if passwd:
            self._args += "--password {}".format(passwd)

        self._svn_loc = address

        err = run_cmd("svn --non-interactive --trust-server-cert {} ls {}".format(self._args, address))[0]
        if err != 0:
            logging.warning(
"""
Could not connect to svn repo '{0}'
This is most likely either a credential, proxy, or network issue .
To check connection and store your credential run 'svn ls {0}' and permanently store your password""".format(address))
            return None

    def fileexists(self, rel_path):
        full_url = os.path.join(self._svn_loc, rel_path)
        stat, out, err = run_cmd("svn --non-interactive --trust-server-cert {} ls {}".format(self._args, full_url))
        if (stat != 0):
            logging.warning("FAIL: SVN repo '{}' does not have file '{}'\nReason:{}\n{}\n".format(self._svn_loc, full_url, out, err))
            return False
        return True

    def getfile(self, rel_path, full_path):
        full_url = os.path.join(self._svn_loc, rel_path)
        stat, output, errput = \
            run_cmd("svn --non-interactive --trust-server-cert {} export {} {}".format(self._args, full_url, full_path))
        if (stat != 0):
            logging.warning("svn export failed with output: {} and errput {}\n".format(output, errput))
            return False
        else:
            logging.info("SUCCESS\n")
            return True

    def getdirectory(self, rel_path, full_path):
        full_url = os.path.join(self._svn_loc, rel_path)
        stat, output, errput = \
            run_cmd("svn --non-interactive --trust-server-cert {} export --force {} {}".format(self._args, full_url, full_path))
        if (stat != 0):
            logging.warning("svn export failed with output: {} and errput {}\n".format(output, errput))
            return False
        else:
            logging.info("SUCCESS\n")
            return True
