"""
API for checking locked files
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.env_build import EnvBuild
from CIME.XML.env_case import EnvCase
from CIME.XML.env_mach_pes import EnvMachPes

import glob

def check_lockedfiles(caseroot=None):
    """
    Check that all lockedfiles match what's in case

    If caseroot is not specified, it is set to the current working directory
    """
    caseroot = os.getcwd() if caseroot is None else caseroot
    lockedfiles = glob.glob(os.path.join(caseroot, "LockedFiles", "*.xml"))
    for lfile in lockedfiles:
        fpart = os.path.basename(lfile)
        cfile = os.path.join(caseroot, fpart)
        if os.path.isfile(cfile):
            objname = fpart.split('.')[0]
            logging.info("Checking file %s"%objname)
            if objname == "env_build":
                f1obj = EnvBuild(caseroot, cfile)
                f2obj = EnvBuild(caseroot, lfile)
            elif objname == "env_mach_pes":
                f1obj = EnvMachPes(caseroot, cfile)
                f2obj = EnvMachPes(caseroot, lfile)
            elif objname == "env_case":
                f1obj = EnvCase(caseroot, cfile)
                f2obj = EnvCase(caseroot, lfile)

            diffs = f1obj.compare_xml(f2obj)
            if diffs:
                logging.warn("File %s has been modified"%lfile)
                for key in diffs.keys():
                    print("  found difference in %s : case %s locked %s" %
                          (key, repr(diffs[key][0]), repr(diffs[key][1])))
                if objname == "env_mach_pes":
                    expect(False, "Invoke case.setup --reset ")
                elif objname == "env_case":
                    expect(False, "Cannot change file env_case.xml, please"
                           " recover the original copy from LockedFiles")
                elif objname == "env_build":
                    logging.warn("Setting build complete to False")
                    f1obj.set_value("BUILD_COMPLETE", False)
                    if "PIO_VERSION" in diffs.keys():
                        f1obj.set_value("BUILD_STATUS", 2)
                        f1obj.write()
                        logging.critical("Changing PIO_VERSION requires running "
                                         "case.build --clean-all and rebuilding")
                    else:
                        f1obj.set_value("BUILD_STATUS", 1)
                        f1obj.write()
                else:
                    # JGF : Error?
                    pass
