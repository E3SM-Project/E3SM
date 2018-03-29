has_gftp = True
try:
    from globus_sdk import TransferClient
except ImportError:
    has_gftp = False

from CIME.Servers.ftp import FTP
from CIME.Servers.svn import SVN
from CIME.Servers.wget import WGET
if has_gftp:
    from CIME.Servers.gftp import GridFTP
