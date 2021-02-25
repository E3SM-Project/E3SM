#pylint: disable=import-error
from distutils.spawn import find_executable
has_gftp = find_executable("globus-url-copy")
has_svn = find_executable("svn")
has_wget = find_executable("wget")
has_ftp = True
try:
    from ftplib import FTP
except ImportError:
    has_ftp = False
if has_ftp:
    from CIME.Servers.ftp import FTP
if has_svn:
    from CIME.Servers.svn import SVN
if has_wget:
    from CIME.Servers.wget import WGET
if has_gftp:
    from CIME.Servers.gftp import GridFTP
