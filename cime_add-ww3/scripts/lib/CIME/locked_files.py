from CIME.XML.standard_module_setup import *
from CIME.utils import safe_copy
from CIME.XML.generic_xml import GenericXML

logger = logging.getLogger(__name__)

LOCKED_DIR = "LockedFiles"

def lock_file(filename, caseroot=None, newname=None):
    expect("/" not in filename, "Please just provide basename of locked file")
    caseroot = os.getcwd() if caseroot is None else caseroot
    newname = filename if newname is None else newname
    fulllockdir = os.path.join(caseroot, LOCKED_DIR)
    if not os.path.exists(fulllockdir):
        os.mkdir(fulllockdir)

    logging.debug("Locking file {}".format(filename))

    # JGF: It is extremely dangerous to alter our database (xml files) without
    # going through the standard API. The copy below invalidates all existing
    # GenericXML instances that represent this file and all caching that may
    # have involved this file. We should probably seek a safer way of locking
    # files.
    safe_copy(os.path.join(caseroot, filename), os.path.join(fulllockdir, newname))
    GenericXML.invalidate(os.path.join(fulllockdir, newname))

def unlock_file(filename, caseroot=None):
    expect("/" not in filename, "Please just provide basename of locked file")
    caseroot = os.getcwd() if caseroot is None else caseroot
    locked_path = os.path.join(caseroot, LOCKED_DIR, filename)
    if os.path.exists(locked_path):
        os.remove(locked_path)

    logging.debug("Unlocking file {}".format(filename))

def is_locked(filename, caseroot=None):
    expect("/" not in filename, "Please just provide basename of locked file")
    caseroot = os.getcwd() if caseroot is None else caseroot
    return os.path.exists(os.path.join(caseroot, LOCKED_DIR, filename))
