"""
common implementation for building components

These are used by components/<model_type>/<component>/cime_config/buildlib
"""
from CIME.XML.standard_module_setup import *
from CIME.case import Case

logger = logging.getLogger(__name__)

def build_x_lib(caseroot, model, comp):
    case = Case(caseroot)
    casetools = case.get_value("CASETOOLS")
    cimeroot = case.get_value("CIMEROOT")
    objroot = case.get_value("OBJROOT")
    mach = case.get_value("MACH")
    libroot = case.get_value("LIBROOT")
    gmake_j = case.get_value("GMAKE_J")
    gmake = case.get_value("GMAKE")

    # Write directory list
    with open('Filepath', 'w') as out:
        out.write(os.path.join(caseroot, "SourceMods", "src.%s", comp)
                  + "\n")
        out.write(os.path.join(cimeroot, "components", "xcpl_comps",
                               "xshare") + "\n")
        out.write(os.path.join(cimeroot, "components", "xcpl_comps",
                               comp, "cpl") + "\n")

    # generate command line argument values
    complib = os.path.join(libroot, "lib%s.a" % model)
    makefile = os.path.join(casetools, "Makefile")
    macfile = os.path.join(caseroot, "Macros.%s" % mach)

    # build
    cmd = "%s complib -j %d MODEL=%s COMPLIB=%s -f %s MACFILE=%s" % \
           (gmake, gmake_j, model, complib, makefile, macfile)
    rc, out, err = run_cmd(cmd, from_dir=os.path.join(objroot,model,'obj'), ok_to_fail=True)
    expect(rc == 0, "Command %s failed rc=%d\nout=%s\nerr=%s" % (cmd, rc,
                                                                 out, err))
    logger.info("Command %s completed with output %s\nerr %s", cmd, out, err)
