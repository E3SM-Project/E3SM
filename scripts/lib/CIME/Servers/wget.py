"""
WGET Server class.  Interact with a server using WGET protocol
"""
# pylint: disable=super-init-not-called
from CIME.XML.standard_module_setup import *
from CIME.Servers.generic_server import GenericServer
logger = logging.getLogger(__name__)

import signal
from contextlib import contextmanager

class TimeoutException(Exception): pass

@contextmanager
def time_limit(seconds):
    #pylint: disable=unused-argument
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


class WGET(GenericServer):
    def __init__(self, address, user='', passwd=''):
        self._args = ''
        if user:
            self._args += "--user {} ".format(user)
        if passwd:
            self._args += "--password {} ".format(passwd)
        self._server_loc = address

    @classmethod
    def wget_login(cls, address, user='', passwd=''):
        args = ''
        if user:
            args += "--user {} ".format(user)
        if passwd:
            args += "--password {} ".format(passwd)

        try:
            with time_limit(30):
                err,_,_ = run_cmd("wget {} --spider {}".format(args, address), verbose=True)
                expect(err == 0,"Could not connect to repo '{0}'\nThis is most likely either a proxy, or network issue .")
            return cls(address,user=user,passwd=passwd)
        except TimeoutException as e:
            logger.warning("wget login timeout! {}".format(e))
            return None


    def fileexists(self, rel_path):
        full_url = os.path.join(self._server_loc, rel_path)
        stat, out, err = run_cmd("wget {} --spider {}".format(self._args, full_url))
        if (stat != 0):
            logging.warning("FAIL: Repo '{}' does not have file '{}'\nReason:{}\n{}\n".format(self._server_loc, full_url, out.encode('utf-8'), err.encode('utf-8')))
            return False
        return True

    def getfile(self, rel_path, full_path):
        full_url = os.path.join(self._server_loc, rel_path)
        stat, output, errput = \
                run_cmd("wget {} {} -nc --output-document {}".format(self._args, full_url, full_path))
        if (stat != 0):
            logging.warning("wget failed with output: {} and errput {}\n".format(output.encode('utf-8'), errput.encode('utf-8')))
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
            run_cmd("wget  {} {} -r -N --no-directories ".format(self._args, full_url+os.sep), from_dir=full_path)
        logger.debug(output)
        logger.debug(errput)
        if (stat != 0):
            logging.warning("wget failed with output: {} and errput {}\n".format(output.encode('utf-8'), errput.encode('utf-8')))
            # wget puts an empty file if it fails.
            try:
                os.remove(full_path)
            except OSError:
                pass
            return False
        else:
            logging.info("SUCCESS\n")
            return True
