"""
functions for building CIME models
"""
from CIME.XML.standard_module_setup  import *
from CIME.utils                 import get_model, append_status
from CIME.preview_namelists     import preview_namelists
from CIME.check_input_data      import check_input_data
from CIME.provenance            import save_build_provenance

import glob, shutil, time, threading, gzip, subprocess

logger = logging.getLogger(__name__)

def stringify_bool(val):
    val = False if val is None else val
    expect(type(val) is bool, "Wrong type for val '%s'" % repr(val))
    return "TRUE" if val else "FALSE"

###############################################################################
def build_model(build_threaded, exeroot, clm_config_opts, incroot, complist,
                lid, caseroot, cimeroot):
###############################################################################

    logs = []

    thread_bad_results = []
    for model, comp, nthrds, _, config_dir in complist:

        # aquap has a dependency on atm so we will build it after the threaded loop
        if comp == "aquap":
            logger.debug("Skip aquap ocn build here")
            continue

        # coupler handled seperately
        if model == "cpl":
            continue

        # special case for clm
        # clm 4_0 is not a shared library and must be built here
        # clm 4_5 and newer is a shared library (but not in ACME) and should be built in build_libraries
        if get_model() != "acme":
            if comp == "clm":
                if "clm4_0" in clm_config_opts:
                    logger.info("         - Building clm4_0 Library ")
                else:
                    continue
        else:
            logger.info("         - Building clm Library ")

        smp = nthrds > 1 or build_threaded

        bldroot = os.path.join(exeroot, model, "obj")
        libroot = os.path.join(exeroot, "lib")
        file_build = os.path.join(exeroot, "%s.bldlog.%s" % (model, lid))
        logger.debug("bldroot is %s" % bldroot)
        logger.debug("libroot is %s" % libroot)

        # make sure bldroot and libroot exist
        for build_dir in [bldroot, libroot]:
            if not os.path.exists(build_dir):
                os.makedirs(build_dir)

        # build the component library
        t = threading.Thread(target=_build_model_thread,
            args=(config_dir, model, caseroot, bldroot, libroot, incroot, file_build,
                  thread_bad_results, smp))
        t.start()

        for mod_file in glob.glob(os.path.join(bldroot, "*_[Cc][Oo][Mm][Pp]_*.mod")):
            shutil.copy(mod_file, incroot)

        logs.append(file_build)

    # Wait for threads to finish
    while(threading.active_count() > 1):
        time.sleep(1)

    # aquap has a dependancy on atm so we build it after the threaded loop

    for model, comp, nthrds, _, config_dir in complist:
        smp = nthrds > 1 or build_threaded
        if comp == "aquap":
            logger.debug("Now build aquap ocn component")
            _build_model_thread(config_dir, comp, caseroot, bldroot, libroot, incroot, file_build,
                                thread_bad_results, smp)

    expect(not thread_bad_results, "\n".join(thread_bad_results))

    #
    # Now build the executable
    #

    cime_model = get_model()
    file_build = os.path.join(exeroot, "%s.bldlog.%s" % (cime_model, lid))

    config_dir = os.path.join(cimeroot, "driver_cpl", "cime_config")
    f = open(file_build, "w")
    stat = run_cmd("%s/buildexe %s %s %s" %
                   (config_dir, caseroot, bldroot, libroot),
                   from_dir=bldroot, verbose=True, arg_stdout=f,
                   arg_stderr=subprocess.STDOUT)[0]
    f.close()

    expect(stat == 0, "ERROR: buildexe failed, cat %s" % file_build)

    # Copy the just-built ${MODEL}.exe to ${MODEL}.exe.$LID
    shutil.copy("%s/%s.exe" % (exeroot, cime_model), "%s/%s.exe.%s" % (exeroot, cime_model, lid))

    logs.append(file_build)

    return logs

###############################################################################
def post_build(case, logs):
###############################################################################

    logdir = case.get_value("LOGDIR")

    #zip build logs to CASEROOT/logs
    if logdir:
        bldlogdir = os.path.join(logdir, "bld")
        if not os.path.exists(bldlogdir):
            os.makedirs(bldlogdir)

    for log in logs:
        logger.debug("Copying build log %s to %s"%(log,bldlogdir))
        with open(log, 'rb') as f_in:
            with gzip.open("%s.gz"%log, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        if "sharedlibroot" not in log:
            shutil.copy("%s.gz"%log,os.path.join(bldlogdir,"%s.gz"%os.path.basename(log)))

    # Set XML to indicate build complete
    case.set_value("BUILD_COMPLETE", True)
    case.set_value("BUILD_STATUS", 0)
    if "SMP_VALUE" in os.environ:
        case.set_value("SMP_BUILD", os.environ["SMP_VALUE"])
    case.flush()

    if os.path.exists("LockedFiles/env_build.xml"):
        os.remove("LockedFiles/env_build.xml")

    shutil.copy("env_build.xml", "LockedFiles")

    # must ensure there's an lid
    lid = os.environ["LID"] if "LID" in os.environ else run_cmd_no_fail("date +%y%m%d-%H%M%S")
    save_build_provenance(case, lid=lid)

###############################################################################
def case_build(caseroot, case, sharedlib_only=False, model_only=False):
###############################################################################

    t1 = time.time()

    expect(not (sharedlib_only and model_only),
           "Contradiction: both sharedlib_only and model_only")

    logger.info("sharedlib_only is %s" % sharedlib_only)
    logger.info("model_only is %s" % model_only)

    expect(os.path.isdir(caseroot), "'%s' is not a valid directory" % caseroot)
    os.chdir(caseroot)

    expect(os.path.exists("case.run"),
           "ERROR: must invoke case.setup script before calling build script ")

    cimeroot = case.get_value("CIMEROOT")

    comp_classes = case.get_value("COMP_CLASSES").split(',')

    if not sharedlib_only:
        check_all_input_data(case)

    run_cmd_no_fail("./Tools/check_lockedfiles --caseroot %s" % caseroot)

    # Retrieve relevant case data
    # This environment variable gets set for cesm Make and
    # needs to be unset before building again.
    if "MODEL" in os.environ.keys():
        del os.environ["MODEL"]
    build_threaded      = case.get_value("BUILD_THREADED")
    casetools           = case.get_value("CASETOOLS")
    exeroot             = case.get_value("EXEROOT")
    incroot             = case.get_value("INCROOT")
    libroot             = case.get_value("LIBROOT")
    sharedlibroot       = case.get_value("SHAREDLIBROOT")

    complist = []
    for comp_class in comp_classes:
        if comp_class == "DRV":
            comp_class = "CPL"
            ninst = 1
            config_dir = None
        else:
            ninst = case.get_value("NINST_%s"%comp_class)
            config_dir = os.path.dirname(case.get_value("CONFIG_%s_FILE"%comp_class))

        comp = case.get_value("COMP_%s"%comp_class)
        thrds =  case.get_value("NTHRDS_%s"%comp_class)
        complist.append((comp_class.lower(), comp, thrds, ninst, config_dir ))
        os.environ["COMP_%s"%comp_class] = comp
    machines_file       = case.get_value("MACHINES_SPEC_FILE")
    ocn_submodel        = case.get_value("OCN_SUBMODEL")
    profile_papi_enable = case.get_value("PROFILE_PAPI_ENABLE")
    compiler            = case.get_value("COMPILER")
    comp_interface      = case.get_value("COMP_INTERFACE")
    mpilib              = case.get_value("MPILIB")
    use_esmf_lib        = case.get_value("USE_ESMF_LIB")
    debug               = case.get_value("DEBUG")
    ninst_build         = case.get_value("NINST_BUILD")
    smp_value           = case.get_value("SMP_VALUE")
    clm_use_petsc       = case.get_value("CLM_USE_PETSC")
    cism_use_trilinos   = case.get_value("CISM_USE_TRILINOS")
    mpasli_use_albany   = case.get_value("MPASLI_USE_ALBANY")
    clm_config_opts     = case.get_value("CLM_CONFIG_OPTS")
    cam_config_opts     = case.get_value("CAM_CONFIG_OPTS")
    pio_config_opts     = case.get_value("PIO_CONFIG_OPTS")
    ninst_value         = case.get_value("NINST_VALUE")
    mach                = case.get_value("MACH")
    os_                 = case.get_value("OS")
    # Load some params into env
    os.environ["CIMEROOT"]             = cimeroot
    os.environ["CASETOOLS"]            = casetools
    os.environ["EXEROOT"]              = exeroot
    os.environ["INCROOT"]              = incroot
    os.environ["LIBROOT"]              = libroot
    os.environ["SHAREDLIBROOT"]        = sharedlibroot
    os.environ["CASEROOT"]             = caseroot
    os.environ["COMPILER"]             = compiler
    os.environ["COMP_INTERFACE"]       = comp_interface
    os.environ["NINST_VALUE"]          = str(ninst_value)
    os.environ["BUILD_THREADED"]       = stringify_bool(build_threaded)
    os.environ["MACH"]                 = mach
    os.environ["USE_ESMF_LIB"]         = stringify_bool(use_esmf_lib)
    os.environ["MPILIB"]               = mpilib
    os.environ["DEBUG"]                = stringify_bool(debug)
    os.environ["OS"]                   = os_
    os.environ["CLM_CONFIG_OPTS"]      = clm_config_opts     if clm_config_opts     is not None else ""
    os.environ["CAM_CONFIG_OPTS"]      = cam_config_opts     if cam_config_opts     is not None else ""
    os.environ["PIO_CONFIG_OPTS"]      = pio_config_opts     if pio_config_opts     is not None else ""
    os.environ["OCN_SUBMODEL"]         = ocn_submodel        if ocn_submodel        is not None else ""
    os.environ["PROFILE_PAPI_ENABLE"]  = stringify_bool(profile_papi_enable)
    os.environ["CLM_USE_PETSC"]        = stringify_bool(clm_use_petsc)
    os.environ["CISM_USE_TRILINOS"]    = stringify_bool(cism_use_trilinos)
    os.environ["MPASLI_USE_ALBANY"]    = stringify_bool(mpasli_use_albany)

    # This is a timestamp for the build , not the same as the testid,
    # and this case may not be a test anyway. For a production
    # experiment there may be many builds of the same case.
    lid               = run_cmd_no_fail("date +%y%m%d-%H%M%S")
    os.environ["LID"] = lid

    # Set the overall USE_PETSC variable to TRUE if any of the
    # *_USE_PETSC variables are TRUE.
    # For now, there is just the one CLM_USE_PETSC variable, but in
    # the future there may be others -- so USE_PETSC will be true if
    # ANY of those are true.

    use_petsc = clm_use_petsc
    case.set_value("USE_PETSC", use_petsc)
    os.environ["USE_PETSC"] = stringify_bool(use_petsc)

    # Set the overall USE_TRILINOS variable to TRUE if any of the
    # *_USE_TRILINOS variables are TRUE.
    # For now, there is just the one CISM_USE_TRILINOS variable, but in
    # the future there may be others -- so USE_TRILINOS will be true if
    # ANY of those are true.

    use_trilinos = False if cism_use_trilinos is None else cism_use_trilinos
    case.set_value("USE_TRILINOS", use_trilinos)
    os.environ["USE_TRILINOS"] = stringify_bool(use_trilinos)

    # Set the overall USE_ALBANY variable to TRUE if any of the
    # *_USE_ALBANY variables are TRUE.
    # For now, there is just the one MPASLI_USE_ALBANY variable, but in
    # the future there may be others -- so USE_ALBANY will be true if
    # ANY of those are true.

    use_albany = stringify_bool(mpasli_use_albany)
    case.set_value("USE_ALBANY", use_albany)
    os.environ["USE_ALBANY"] = use_albany

    # Load modules
    env_module = case.get_env("mach_specific")
    env_module.load_env_for_case(compiler=case.get_value("COMPILER"),
                                 debug=case.get_value("DEBUG"),
                                 mpilib=case.get_value("MPILIB"))


    # Need to flush case xml to disk before calling preview_namelists
    case.flush()

    if not sharedlib_only:
        preview_namelists(case)

    build_checks(case, build_threaded, comp_interface, use_esmf_lib, debug, compiler, mpilib,
                 sharedlibroot, complist, ninst_build, smp_value)

    t2 = time.time()
    logs = []

    if not model_only:
        logs = build_libraries(case, exeroot, caseroot, cimeroot, libroot, mpilib, lid,
                               machines_file)

    if not sharedlib_only:
        logs.extend(build_model(build_threaded, exeroot, clm_config_opts, incroot, complist,
                                lid, caseroot, cimeroot))

    if not sharedlib_only:
        post_build(case, logs)

    t3 = time.time()

    logger.info("Time spent not building: %f sec" % (t2 - t1))
    logger.info("Time spent building: %f sec" % (t3 - t2))

    return True

###############################################################################
def check_all_input_data(case):
###############################################################################

    success = check_input_data(case=case, download=True)
    expect(success, "Failed to download input data")

    get_refcase  = case.get_value("GET_REFCASE")
    run_type     = case.get_value("RUN_TYPE")
    continue_run = case.get_value("CONTINUE_RUN")

    # We do not fully populate the inputdata directory on every
    # machine and do not expect every user to download the 3TB+ of
    # data in our inputdata repository. This code checks for the
    # existence of inputdata in the local inputdata directory and
    # attempts to download data from the server if it's needed and
    # missing.
    if get_refcase and run_type != "startup" and not continue_run:
        din_loc_root = case.get_value("DIN_LOC_ROOT")
        run_refdate  = case.get_value("RUN_REFDATE")
        run_refcase  = case.get_value("RUN_REFCASE")
        run_refdir   = case.get_value("RUN_REFDIR")
        rundir       = case.get_value("RUNDIR")

        refdir = os.path.join(din_loc_root, run_refdir, run_refcase, run_refdate)
        expect(os.path.isdir(refdir),
"""
*****************************************************************
prestage ERROR: $refdir is not on local disk
obtain this data from the svn input data repository
> mkdir -p %s
> cd %s
> cd ..
> svn export --force https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/%s
or set GET_REFCASE to FALSE in env_run.xml
and prestage the restart data to $RUNDIR manually
*****************************************************************""" % (refdir, refdir, refdir))

        logger.info(" - Prestaging REFCASE (%s) to %s" % (refdir, rundir))

        # prestage the reference case's files.

        if (not os.path.exists(rundir)):
            logger.debug("Creating run directory: %s"%rundir)
            os.makedirs(rundir)

        for rcfile in glob.iglob(os.path.join(refdir,"*%s*"%run_refcase)):
            logger.debug("Staging file %s"%rcfile)
            rcbaseline = os.path.basename(rcfile)
            if not os.path.exists("%s/%s" % (rundir, rcbaseline)):
                os.symlink(rcfile, "%s/%s" % ((rundir, rcbaseline)))

        # copy the refcases' rpointer files to the run directory
        for rpointerfile in  glob.iglob(os.path.join("%s","*rpointer*") % (refdir)):
            logger.debug("Copy rpointer %s"%rpointerfile)
            shutil.copy(rpointerfile, rundir)


        for cam2file in  glob.iglob(os.path.join("%s","*.cam2.*") % rundir):
            camfile = cam2file.replace("cam2", "cam")
            os.symlink(cam2file, camfile)


###############################################################################
def build_checks(case, build_threaded, comp_interface, use_esmf_lib, debug, compiler, mpilib,
                 sharedlibroot, complist, ninst_build, smp_value):
###############################################################################

    ninst_value  = case.get_value("NINST_VALUE")
    smp_build    = case.get_value("SMP_BUILD")
    build_status = case.get_value("BUILD_STATUS")

    smpstr = ""
    inststr = ""
    for model, _, nthrds, ninst, _ in complist:
        if nthrds > 1:
            build_threaded = True
        if build_threaded:
            smpstr += "%s1"%model[0]
        else:
            smpstr += "%s0"%model[0]
        inststr += "%s%d"%(model[0],ninst)

    if build_threaded:
        os.environ["SMP"] = "TRUE"
    else:
        os.environ["SMP"] = "FALSE"
    case.set_value("SMP_VALUE", smpstr)
    os.environ["SMP_VALUE"] = smpstr
    case.set_value("NINST_VALUE", inststr)
    os.environ["NINST_VALUE"] = inststr


    debugdir = "debug" if debug else "nodebug"
    threaddir = "threads" if (os.environ["SMP"] == "TRUE" or build_threaded) else "nothreads"
    sharedpath = os.path.join(sharedlibroot, compiler, mpilib, debugdir, threaddir)
    logger.debug("compiler=%s mpilib=%s debugdir=%s threaddir=%s"%
                 (compiler,mpilib,debugdir,threaddir))
    os.environ["SHAREDPATH"] = sharedpath

    expect( ninst_build == ninst_value or ninst_build == "0",
            """
ERROR, NINST VALUES HAVE CHANGED
  NINST_BUILD = %s
  NINST_VALUE = %s
  A manual clean of your obj directories is strongly recommended
  You should execute the following:
    ./case.build --clean
  Then rerun the build script interactively
  ---- OR ----
  You can override this error message at your own risk by executing:
    ./xmlchange -file env_build.xml -id NINST_BUILD -val 0
  Then rerun the build script interactively
""" % (ninst_build, ninst_value))

    expect( smp_build == smpstr or smp_build == "0",
            """
ERROR, SMP VALUES HAVE CHANGED
  SMP_BUILD = %s
  SMP_VALUE = %s
  smpstr = %s
  A manual clean of your obj directories is strongly recommended
  You should execute the following:
    ./case.build --clean
  Then rerun the build script interactively
  ---- OR ----
  You can override this error message at your own risk by executing:
    ./xmlchange -file env_build.xml -id SMP_BUILD -val 0
  Then rerun the build script interactively
""" % (smp_build, smp_value, smpstr))

    expect(build_status == 0,
           """
ERROR env_build HAS CHANGED
  A manual clean of your obj directories is required
  You should execute the following:
    ./case.build --clean-all
""")

    expect(comp_interface != "ESMF" or use_esmf_lib,
           """
ERROR COMP_INTERFACE IS ESMF BUT USE_ESMF_LIB IS NOT TRUE
  SET USE_ESMF_LIB to TRUE with:
    ./xmlchange -file env_build.xml -id USE_ESMF_LIB -value TRUE
""")

    expect(mpilib != "mpi-serial" or not use_esmf_lib,
           """
ERROR MPILIB is mpi-serial and USE_ESMF_LIB IS TRUE
  MPILIB can only be used with an ESMF library built with mpiuni on
  Set USE_ESMF_LIB to FALSE with
    ./xmlchange -file env_build.xml -id USE_ESMF_LIB -val FALSE
  ---- OR ----
  Make sure the ESMF_LIBDIR used was built with mipuni (or change it to one that was)
  And comment out this if block in Tools/models_buildexe
""")

    case.set_value("BUILD_COMPLETE", False)

    case.flush()

###############################################################################
def build_libraries(case, exeroot, caseroot, cimeroot, libroot, mpilib, lid, machines_file):
###############################################################################

    if (mpilib == "mpi-serial"):
        for header_to_copy in glob.glob(os.path.join(cimeroot, "externals/mct/mpi-serial/*.h")):
            shutil.copy(header_to_copy, os.path.join(libroot, "include"))

    sharedpath = os.environ["SHAREDPATH"]
    shared_lib = os.path.join(sharedpath, "lib")
    shared_inc = os.path.join(sharedpath, "include")
    for shared_item in [shared_lib, shared_inc]:
        if (not os.path.exists(shared_item)):
            os.makedirs(shared_item)

    libs = ["mct", "gptl", "pio", "csm_share"]
    logs = []

    for lib in libs:
        full_lib_path = os.path.join(sharedpath, lib)
        if (not os.path.exists(full_lib_path)):
            os.makedirs(full_lib_path)

        file_build = os.path.join(sharedpath, "%s.bldlog.%s" % (lib, lid))
        my_file = os.path.join(os.path.dirname(machines_file), "buildlib.%s" % lib)
        with open(file_build, "w") as fd:
            fd.write("Current env:\n%s" % "\n".join(["  %s = %s" % (env, os.environ[env]) for env in sorted(os.environ)]))
            fd.flush()
            stat = run_cmd("%s %s %s 2>&1" %
                           (my_file, sharedpath, caseroot),
                           from_dir=exeroot,
                           verbose=True, arg_stdout=fd)[0]
        expect(stat == 0, "ERROR: buildlib.%s failed, cat %s" % (lib, file_build))
        logs.append(file_build)

    # clm not a shared lib for ACME
    if get_model() != "acme":
        comp_lnd = case.get_value("COMP_LND")
        clm_config_opts = case.get_value("CLM_CONFIG_OPTS")
        if comp_lnd == "clm" and not "clm4_0" in clm_config_opts:
            logging.info("         - Building clm4_5/clm5_0 Library ")
            esmfdir = "esmf" if case.get_value("USE_ESMF_LIB") else "noesmf"
            sharedpath = os.environ["SHAREDPATH"]
            bldroot = os.path.join(sharedpath, case.get_value("COMP_INTERFACE"), esmfdir, "clm","obj" )
            libroot = os.path.join(sharedpath, case.get_value("COMP_INTERFACE"), esmfdir, "lib")
            incroot = os.path.join(sharedpath, case.get_value("COMP_INTERFACE"), esmfdir, "include")
            file_build = os.path.join(exeroot, "lnd.bldlog.%s" %  lid)
            config_lnd_dir = os.path.dirname(case.get_value("CONFIG_LND_FILE"))

            for ndir in [bldroot, libroot]:
                if (not os.path.isdir(ndir)):
                    os.makedirs(ndir)

            smp = "SMP" in os.environ and os.environ["SMP"] == "TRUE"
            _build_model_thread(config_lnd_dir, "lnd", caseroot, bldroot, libroot, incroot, file_build, logs, smp)

    return logs

###############################################################################
def _build_model_thread(config_dir, compclass, caseroot, bldroot, libroot, incroot, file_build,
                        thread_bad_results, smp):
###############################################################################
    with open(file_build, "w") as fd:
        stat = run_cmd("MODEL=%s SMP=%s %s/buildlib %s %s %s " %
                       (compclass, stringify_bool(smp), config_dir, caseroot, bldroot, libroot),
                       from_dir=bldroot, verbose=True, arg_stdout=fd,
                       arg_stderr=subprocess.STDOUT)[0]
    if (stat != 0):
        thread_bad_results.append("ERROR: %s.buildlib failed, see %s" % (compclass, file_build))

    for mod_file in glob.glob(os.path.join(bldroot, "*_[Cc][Oo][Mm][Pp]_*.mod")):
        shutil.copy(mod_file, incroot)

###############################################################################
def clean(case, cleanlist=None):
###############################################################################

    clm_config_opts = case.get_value("CLM_CONFIG_OPTS")
    comp_lnd = case.get_value("COMP_LND")
    if cleanlist is None:
        cleanlist = case.get_value("COMP_CLASSES").split(',')
        cleanlist = [x.lower().replace('drv','cpl') for x in cleanlist]
        testcase        = case.get_value("TESTCASE")
        # we only want to clean clm here if it is clm4_0 otherwise remove
        # it from the cleanlist
        if testcase is not None and comp_lnd == "clm" and \
                clm_config_opts is not None and "lnd" in cleanlist and \
                "clm4_0" not in clm_config_opts:
            cleanlist.remove('lnd')

    debug           = case.get_value("DEBUG")
    use_esmf_lib    = case.get_value("USE_ESMF_LIB")
    build_threaded  = case.get_value("BUILD_THREADED")
    gmake           = case.get_value("GMAKE")
    caseroot        = case.get_value("CASEROOT")
    casetools       = case.get_value("CASETOOLS")

    os.environ["DEBUG"]           = stringify_bool(debug)
    os.environ["USE_ESMF_LIB"]    = stringify_bool(use_esmf_lib)
    os.environ["BUILD_THREADED"]  = stringify_bool(build_threaded)
    os.environ["CASEROOT"]        = case.get_value("CASEROOT")
    os.environ["COMP_INTERFACE"]  = case.get_value("COMP_INTERFACE")
    os.environ["PIO_VERSION"]     = str(case.get_value("PIO_VERSION"))
    os.environ["CLM_CONFIG_OPTS"] = clm_config_opts  if clm_config_opts is not None else ""

    cmd = gmake + " -f " + casetools + "/Makefile"
    for item in cleanlist:
        cmd = cmd + " clean" + item
    logger.info("calling %s "%(cmd))
    run_cmd_no_fail(cmd)

    # unlink Locked files directory
    locked_env_build = os.path.join(caseroot,"LockedFiles/env_build.xml")
    if os.path.isfile(locked_env_build):
        os.unlink(locked_env_build)

    # reset following values in xml files
    case.set_value("SMP_BUILD",str(0))
    case.set_value("NINST_BUILD",str(0))
    case.set_value("BUILD_STATUS",str(0))
    case.set_value("BUILD_COMPLETE","FALSE")
    case.flush()

    # append call of to CaseStatus
    msg = "cleanbuild %s "%" ".join(cleanlist)
    append_status(msg, caseroot=caseroot, sfile="CaseStatus")
