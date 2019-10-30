"""
WGET Server class.  Interact with a server using WGET protocol
"""
# pylint: disable=super-init-not-called
from CIME.XML.standard_module_setup import *
from CIME.Servers.generic_server import GenericServer
logger = logging.getLogger(__name__)

class WGET(GenericServer):
    def __init__(self, address, user='', passwd=''):
        self._args = ''
        if user:
            self._args += "--user {} ".format(user)
        if passwd:
            self._args += "--password {} ".format(passwd)
        self._server_loc = address

        cmd = "wget {} --no-check-certificate --spider {}".format(self._args, address)
        err, output, _ = run_cmd(cmd, combine_output=True)
        expect(err == 0,"Could not connect to repo via '{}'\nThis is most likely either a proxy, or network issue.\nOutput:\n{}".format(cmd, output))

    def fileexists(self, rel_path):
        full_url = os.path.join(self._server_loc, rel_path)
        stat, out, err = run_cmd("wget {} --no-check-certificate --spider {}".format(self._args, full_url))
        if (stat != 0):
            logging.warning("FAIL: Repo '{}' does not have file '{}'\nReason:{}\n{}\n".format(self._server_loc, full_url, out, err))
            return False
        return True

    def getfile(self, rel_path, full_path):
        full_url = os.path.join(self._server_loc, rel_path)
        stat, output, errput = \
                run_cmd("wget {} {} -nc --no-check-certificate --output-document {}".format(self._args, full_url, full_path))
        if (stat != 0):
            logging.warning("wget failed with output: {} and errput {}\n".format(output, errput))
            # wget puts an empty file if it fails.
            try:
                os.remove(full_path)
            except OSError:
                pass
            return False
        else:
            logging.info("SUCCESS\n")
            return True

    def getdirectory(self, rel_path, full_path):
        full_url = os.path.join(self._server_loc, rel_path)
        stat, output, errput = \
            run_cmd("wget  {} {} -r -N --no-check-certificate --no-directories ".format(self._args, full_url+os.sep), from_dir=full_path)
        logger.debug(output)
        logger.debug(errput)
        if (stat != 0):
            logging.warning("wget failed with output: {} and errput {}\n".format(output, errput))
            # wget puts an empty file if it fails.
            try:
                os.remove(full_path)
            except OSError:
                pass
            return False
        else:
            logging.info("SUCCESS\n")
            return True
