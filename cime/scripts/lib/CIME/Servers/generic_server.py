"""
Generic Server class.  There should be little or no functionality in this class, it serves only
to make sure that specific server classes maintain a consistant argument list and functionality
so that they are interchangable objects
"""
# pylint: disable=unused-argument

from CIME.XML.standard_module_setup import *
from socket import _GLOBAL_DEFAULT_TIMEOUT
logger = logging.getLogger(__name__)

class GenericServer(object):
    def __init__(self, host=' ',user=' ', passwd=' ', acct=' ', timeout=_GLOBAL_DEFAULT_TIMEOUT):
        raise NotImplementedError

    def fileexists(self, rel_path):
        '''  Returns True if rel_path exists on server '''
        raise NotImplementedError

    def getfile(self, rel_path, full_path):
        ''' Get file from rel_path on server and place in location full_path on client
        fail if full_path already exists on client, return True if successful '''
        raise NotImplementedError
