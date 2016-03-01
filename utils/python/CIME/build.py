"""
functions for building CIME models
"""
from XML.standard_module_setup import *
from case import Case
from utils import expect, run_cmd, get_model
from env_module import EnvModule
from CIME.preview_namelists import preview_namelists
from CIME.check_input_data import check_input_data

import glob, shutil, time, threading

###############################################################################
def build_model(case, build_threaded, exeroot, clm_config_opts, incroot,
                comp_atm, comp_lnd, comp_ice, comp_ocn, comp_glc, comp_wav, comp_rof,
                nthrds_atm, nthrds_lnd, nthrds_ice, nthrds_ocn, nthrds_glc, nthrds_wav,
                nthrds_rof, lid, caseroot, cimeroot, use_esmf_lib, comp_interface):
###############################################################################
    config_atm_dir = os.path.dirname(case.get_value("CONFIG_ATM_FILE"))
    config_lnd_dir = os.path.dirname(case.get_value("CONFIG_LND_FILE"))
    config_ice_dir = os.path.dirname(case.get_value("CONFIG_ICE_FILE"))
    config_ocn_dir = os.path.dirname(case.get_value("CONFIG_OCN_FILE"))
    config_glc_dir = os.path.dirname(case.get_value("CONFIG_GLC_FILE"))
    config_wav_dir = os.path.dirname(case.get_value("CONFIG_WAV_FILE"))
    config_rof_dir = os.path.dirname(case.get_value("CONFIG_ROF_FILE"))

    # Also defines build order
    models_build_data = [("atm", comp_atm, config_atm_dir, nthrds_atm),
                         ("lnd", comp_lnd, config_lnd_dir, nthrds_lnd),
                         ("ice", comp_ice, config_ice_dir, nthrds_ice),
                         ("ocn", comp_ocn, config_ocn_dir, nthrds_ocn),
                         ("glc", comp_glc, config_glc_dir, nthrds_glc),
                         ("wav", comp_wav, config_wav_dir, nthrds_wav),
                         ("rof", comp_rof, config_rof_dir, nthrds_rof)]

    logs = []
    overall_smp = os.environ["SMP"]
    sharedpath = os.environ["SHAREDPATH"]

    thread_bad_results = []
    for model, comp, config_dir, nthrds in models_build_data:
        os.environ["MODEL"] = model
        if nthrds > 1 or build_threaded == "TRUE":
            os.environ["SMP"] = "TRUE"
        else:
            os.environ["SMP"] = "FALSE"

        bldroot = exeroot # What is this for?

        objdir = os.path.join(exeroot, model, "obj")
        libdir = os.path.join(exeroot, model)
        compspec = comp

        # Special case for clm
        if comp == "clm":
            esmfdir = "esmf" if use_esmf_lib == "TRUE" else "noesmf"
            if "clm4_0" in clm_config_opts:
                logging.info("         - Building clm4_0 Library ")
                compspec = "lnd"
            else:
                logging.info("         - Building clm4_5/clm5_0 Library ")
                bldroot = os.path.join(sharedpath, comp_interface, esmfdir)
                objdir = os.path.join(bldroot, comp, "obj")
                libdir = os.path.join(bldroot, "lib")
                compspec = "clm"

        logging.debug("bldroot is %s" % bldroot)
        logging.debug("objdir is %s" % objdir)
        logging.debug("libdir is %s" % libdir)

        # Make sure obj, lib dirs exist
        for build_dir in [objdir, libdir]:
            if not os.path.exists(build_dir):
                os.makedirs(build_dir)

        file_build = os.path.join(exeroot, "%s.bldlog.%s" % (model, lid))

        # build the component library
        t = threading.Thread(target=_build_model_thread,
            args=(config_dir, caseroot, bldroot, compspec, file_build,
                  exeroot, model, comp, objdir, incroot, thread_bad_results))
        t.start()

        for mod_file in glob.glob(os.path.join(objdir, "*_[Cc][Oo][Mm][Pp]_*.mod")):
            shutil.copy(mod_file, incroot)

        logs.append(file_build)

    # Wait for threads to finish
    while(threading.active_count() > 1):
        time.sleep(1)

    expect(not thread_bad_results, "\n".join(thread_bad_results))

    #
    # Now build the executable
    #

    os.environ["SMP"] = overall_smp

    cime_model = get_model()
    file_build = os.path.join(exeroot, "%s.bldlog.%s" % (cime_model, lid))

    stat = run_cmd("%s/driver_cpl/cime_config/buildexe %s >> %s 2>&1" %
                   (cimeroot, caseroot, file_build),
                   ok_to_fail=True,
                   verbose=True)[0]
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
            run_cmd("gzip %s" % log, from_dir=bldlogdir)

    # Set XML to indicate build complete
    case.set_value("BUILD_COMPLETE", "TRUE")
    case.set_value("BUILD_STATUS", "0")
    case.set_value("SMP_BUILD", os.environ["SMP_VALUE"])
    case.flush()

    if os.path.exists("LockedFiles/env_build.xml"):
        os.remove("LockedFiles/env_build.xml")

    shutil.copy("env_build.xml", "LockedFiles")

###############################################################################
def case_build(caseroot, case=None, testmode=False, sharedlib_only=False, model_only=False):
###############################################################################
    t1 = time.time()

    expect(not (sharedlib_only and model_only),
           "Contradiction: both sharedlib_only and model_only")

    logging.info("sharedlib_only is %s" % sharedlib_only)
    logging.info("model_only is %s" % model_only)

    expect(os.path.isdir(caseroot), "'%s' is not a valid directory" % caseroot)
    os.chdir(caseroot)

    expect(os.path.exists("case.run"),
           "ERROR: must invoke case.setup script before calling build script ")

    case = Case() if case is None else case
    testcase = case.get_value("TESTCASE")
    cimeroot = case.get_value("CIMEROOT")
    expect(not (testcase is not None and
                os.path.exists("%s/scripts/Testing/Testcases/%s_build.csh" %
                               (cimeroot, testcase)) and not testmode),
           "%s build must be invoked via case.testbuild script" % testcase)

    if not sharedlib_only:
        check_all_input_data(case)

    run_cmd("./Tools/check_lockedfiles --caseroot %s" % caseroot)

    # Retrieve relevant case data
    build_threaded      = case.get_value("BUILD_THREADED")
    casetools           = case.get_value("CASETOOLS")
    exeroot             = case.get_value("EXEROOT")
    incroot             = case.get_value("INCROOT")
    libroot             = case.get_value("LIBROOT")
    sharedlibroot       = case.get_value("SHAREDLIBROOT")
    comp_atm            = case.get_value("COMP_ATM")
    comp_lnd            = case.get_value("COMP_LND")
    comp_ice            = case.get_value("COMP_ICE")
    comp_ocn            = case.get_value("COMP_OCN")
    comp_glc            = case.get_value("COMP_GLC")
    comp_wav            = case.get_value("COMP_WAV")
    comp_rof            = case.get_value("COMP_ROF")
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
    comp_cpl            = case.get_value("COMP_CPL")
    machines_file       = case.get_value("MACHINES_SPEC_FILE")
    ocn_submodel        = case.get_value("OCN_SUBMODEL")
    profile_papi_enable = case.get_value("PROFILE_PAPI_ENABLE")
    nthrds_cpl          = int(case.get_value("NTHRDS_CPL"))
    nthrds_atm          = int(case.get_value("NTHRDS_ATM"))
    nthrds_lnd          = int(case.get_value("NTHRDS_LND"))
    nthrds_ice          = int(case.get_value("NTHRDS_ICE"))
    nthrds_ocn          = int(case.get_value("NTHRDS_OCN"))
    nthrds_glc          = int(case.get_value("NTHRDS_GLC"))
    nthrds_wav          = int(case.get_value("NTHRDS_WAV"))
    nthrds_rof          = int(case.get_value("NTHRDS_ROF"))

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
    os.environ["NINST_VALUE"]          = ninst_value
    os.environ["BUILD_THREADED"]       = build_threaded
    os.environ["MACH"]                 = mach
    os.environ["USE_ESMF_LIB"]         = use_esmf_lib
    os.environ["MPILIB"]               = mpilib
    os.environ["DEBUG"]                = debug
    os.environ["OS"]                   = os_
    os.environ["COMP_CPL"]             = comp_cpl
    os.environ["COMP_ATM"]             = comp_atm
    os.environ["COMP_LND"]             = comp_lnd
    os.environ["COMP_ICE"]             = comp_ice
    os.environ["COMP_OCN"]             = comp_ocn
    os.environ["COMP_GLC"]             = comp_glc
    os.environ["COMP_WAV"]             = comp_wav
    os.environ["COMP_ROF"]             = comp_rof
    os.environ["CLM_CONFIG_OPTS"]      = clm_config_opts     if clm_config_opts     is not None else ""
    os.environ["CAM_CONFIG_OPTS"]      = cam_config_opts     if cam_config_opts     is not None else ""
    os.environ["PIO_CONFIG_OPTS"]      = pio_config_opts     if pio_config_opts     is not None else ""
    os.environ["OCN_SUBMODEL"]         = ocn_submodel        if ocn_submodel        is not None else ""
    os.environ["PROFILE_PAPI_ENABLE"]  = profile_papi_enable if profile_papi_enable is not None else ""
    os.environ["CLM_USE_PETSC"]        = clm_use_petsc       if clm_use_petsc       is not None else ""
    os.environ["CISM_USE_TRILINOS"]    = cism_use_trilinos   if cism_use_trilinos   is not None else ""
    os.environ["MPASLI_USE_ALBANY"]    = mpasli_use_albany   if mpasli_use_albany   is not None else ""

    # This is a timestamp for the build , not the same as the testid,
    # and this case may not be a test anyway. For a production
    # experiment there may be many builds of the same case.
    lid               = run_cmd("date +%y%m%d-%H%M%S")
    os.environ["LID"] = lid

    # Set the overall USE_PETSC variable to TRUE if any of the
    # XXX_USE_PETSC variables are TRUE.
    # For now, there is just the one CLM_USE_PETSC variable, but in
    # the future there may be others -- so USE_PETSC will be true if
    # ANY of those are true.

    use_petsc = "TRUE" if clm_use_petsc == "TRUE" else "FALSE"
    case.set_value("USE_PETSC", use_petsc)
    os.environ["USE_PETSC"] = use_petsc

    # Set the overall USE_TRILINOS variable to TRUE if any of the
    # XXX_USE_TRILINOS variables are TRUE.
    # For now, there is just the one CISM_USE_TRILINOS variable, but in
    # the future there may be others -- so USE_TRILINOS will be true if
    # ANY of those are true.

    use_trilinos = "TRUE" if cism_use_trilinos == "TRUE" else "FALSE"
    case.set_value("USE_TRILINOS", use_trilinos)
    os.environ["USE_TRILINOS"] = use_trilinos

    # Set the overall USE_ALBANY variable to TRUE if any of the
    # XXX_USE_ALBANY variables are TRUE.
    # For now, there is just the one MPASLI_USE_ALBANY variable, but in
    # the future there may be others -- so USE_ALBANY will be true if
    # ANY of those are true.

    use_albany = "TRUE" if mpasli_use_albany == "TRUE" else "FALSE"
    case.set_value("USE_ALBANY", use_albany)
    os.environ["USE_ALBANY"] = use_albany

    # Load modules
    env_module = EnvModule(mach, compiler, cimeroot, caseroot, mpilib, debug)
    env_module.load_env_for_case()

    # Need to flush case xml to disk before calling preview_namelists
    case.flush()

    if not sharedlib_only:
        run_cmd("./preview_namelists")

    build_checks(case, build_threaded, comp_interface, use_esmf_lib, compiler, mpilib, debug, sharedlibroot,
                 nthrds_cpl, nthrds_atm, nthrds_lnd, nthrds_ice, nthrds_ocn, nthrds_glc, nthrds_wav,
                 nthrds_rof, ninst_build, smp_value)

    t2 = time.time()
    logs = []

    if not model_only:
        logs = build_libraries(exeroot, caseroot, cimeroot, libroot, mpilib, lid, machines_file)
    if not sharedlib_only:
        logs.extend(build_model(case, build_threaded, exeroot, clm_config_opts, incroot,
                                comp_atm,   comp_lnd,   comp_ice,   comp_ocn,   comp_glc,   comp_wav,   comp_rof,
                                nthrds_atm, nthrds_lnd, nthrds_ice, nthrds_ocn, nthrds_glc, nthrds_wav, nthrds_rof,
                                lid, caseroot, cimeroot, use_esmf_lib, comp_interface))

    if not sharedlib_only:
        post_build(case, logs)

    t3 = time.time()

    logging.info("Time spent not building: %f sec" % (t2 - t1))
    logging.info("Time spent building: %f sec" % (t3 - t2))

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
    if get_refcase == "TRUE" and run_type != "startup" and continue_run == "FALSE":
        din_loc_root = case.get_value("DIN_LOC_ROOT")
        run_refdate  = case.get_value("RUN_REFDATE")
        run_refcase  = case.get_value("RUN_REFCASE")
        run_refdir   = case.get_value("RUN_REFDIR")
        rundir       = case.get_value("RUNDIR")

        refdir = os.path.join(run_refdir, run_refcase, run_refdate)
        expect(os.path.isdir(refdir),
"""
*****************************************************************
ccsm_prestage ERROR: $refdir is not on local disk
obtain this data from the svn input data repository
> mkdir -p %s
> cd %s
> cd ..
> svn export --force https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/%s
or set GET_REFCASE to FALSE in env_run.xml
and prestage the restart data to $RUNDIR manually
*****************************************************************""" % (refdir, refdir, refdir))

        logging.info(" - Prestaging REFCASE (%s) to %s" % (refdir, rundir))

        # prestage the reference case's files.

        if (not os.path.exists(rundir)):
            os.makedirs(rundir)

        refcasefiles = glob.glob("%s/%s/*%s*" % (din_loc_root, refdir, run_refcase))
        for rcfile in refcasefiles:
            rcbaseline = os.path.basename(rcfile)
            if not os.path.exists("%s/%s" % (rundir, rcbaseline)):
                os.symlink(rcfile, "%s/%s" % ((rundir, rcbaseline)))

            # copy the refcases' rpointer files to the run directory
            rpointerfiles = glob.glob("%s/%s/*rpointer*" % (din_loc_root, refdir))
            for rpointerfile in rpointerfiles:
                shutil.copy(rpointerfile, rundir)

            cam2_list = glob.glob("%s/*.cam2.*" % rundir)
            for cam2file in cam2_list:
                camfile = cam2file.replace("cam2", "cam")
                os.symlink(cam2file, camfile)

            allrundirfiles = glob.glob("%s/*" % rundir)
            for runfile in allrundirfiles:
                os.chmod(runfile, 0755)

###############################################################################
def build_checks(case, build_threaded, comp_interface, use_esmf_lib, debug, compiler, mpilib,
                 sharedlibroot, nthrds_cpl, nthrds_atm, nthrds_lnd, nthrds_ice, nthrds_ocn,
                 nthrds_glc, nthrds_wav, nthrds_rof, ninst_build, smp_value):
###############################################################################
    ninst_atm    = int(case.get_value("NINST_ATM"))
    ninst_lnd    = int(case.get_value("NINST_LND"))
    ninst_ice    = int(case.get_value("NINST_ICE"))
    ninst_ocn    = int(case.get_value("NINST_OCN"))
    ninst_glc    = int(case.get_value("NINST_GLC"))
    ninst_wav    = int(case.get_value("NINST_WAV"))
    ninst_rof    = int(case.get_value("NINST_ROF"))
    ninst_value  = case.get_value("NINST_VALUE")
    smp_build    = case.get_value("SMP_BUILD")
    build_status = int(case.get_value("BUILD_STATUS"))

    atmstr = 1 if (build_threaded == "TRUE" or nthrds_atm > 1) else 0
    lndstr = 1 if (build_threaded == "TRUE" or nthrds_lnd > 1) else 0
    icestr = 1 if (build_threaded == "TRUE" or nthrds_ice > 1) else 0 # not in perl
    ocnstr = 1 if (build_threaded == "TRUE" or nthrds_ocn > 1) else 0
    rofstr = 1 if (build_threaded == "TRUE" or nthrds_rof > 1) else 0
    glcstr = 1 if (build_threaded == "TRUE" or nthrds_glc > 1) else 0
    wavstr = 1 if (build_threaded == "TRUE" or nthrds_wav > 1) else 0
    cplstr = 1 if (build_threaded == "TRUE" or nthrds_cpl > 1) else 0

    if (nthrds_atm > 1 or nthrds_lnd > 1 or nthrds_ice > 1 or nthrds_ocn > 1 or
        nthrds_rof > 1 or nthrds_glc > 1 or nthrds_wav > 1 or nthrds_cpl > 1):
        os.environ["SMP"] = "TRUE"
    else:
        os.environ["SMP"] = "FALSE"

    debugdir = "debug" if debug else "nodebug"
    threaddir = "threads" if (os.environ["SMP"] == "TRUE" or build_threaded == "TRUE") else "nothreads"
    sharedpath = os.path.join(sharedlibroot, compiler, mpilib, debugdir, threaddir)
    os.environ["SHAREDPATH"] = sharedpath

    smpstr = "a%dl%dr%di%do%dg%dw%dc%d" % \
        (atmstr, lndstr, rofstr, icestr, ocnstr,  glcstr, wavstr, cplstr)

    case.set_value("SMP_VALUE", smpstr)
    os.environ["SMP_VALUE"] = smpstr

    inststr = "a%dl%dr%di%do%dg%dw%d" % \
        (ninst_atm, ninst_lnd, ninst_rof, ninst_ice, ninst_ocn, ninst_glc, ninst_wav)

    case.set_value("NINST_VALUE", inststr)
    os.environ["NINST_VALUE"] = inststr

    expect( ninst_build == ninst_value or ninst_build == "0",
            """
ERROR, NINST VALUES HAVE CHANGED
  NINST_BUILD = %s
  NINST_VALUE = %s
  A manual clean of your obj directories is strongly recommended
  You should execute the following:
    ./case.clean_build
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
    ./case.clean_build
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
    ./case.clean_build all
""")

    expect(comp_interface != "ESMF" or use_esmf_lib == "TRUE",
           """
ERROR COMP_INTERFACE IS ESMF BUT USE_ESMF_LIB IS NOT TRUE
  SET USE_ESMF_LIB to TRUE with:
    ./xmlchange -file env_build.xml -id USE_ESMF_LIB -value TRUE
""")

    expect(mpilib != "mpi-serial" or use_esmf_lib != "TRUE",
           """
ERROR MPILIB is mpi-serial and USE_ESMF_LIB IS TRUE
  MPILIB can only be used with an ESMF library built with mpiuni on
  Set USE_ESMF_LIB to FALSE with
    ./xmlchange -file env_build.xml -id USE_ESMF_LIB -val FALSE
  ---- OR ----
  Make sure the ESMF_LIBDIR used was built with mipuni (or change it to one that was)
  And comment out this if block in Tools/models_buildexe
""")

    case.set_value("BUILD_COMPLETE", "FALSE")

    case.flush()

###############################################################################
def build_libraries(exeroot, caseroot, cimeroot, libroot, mpilib, lid, machines_file):
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
        with open(file_build, "w") as fd:
            fd.write("Current env:\n%s" % "\n".join(["  %s = %s" % (env, os.environ[env]) for env in sorted(os.environ)]))

        my_file = os.path.join(os.path.dirname(machines_file), "buildlib.%s" % lib)
        stat = run_cmd("%s %s %s >> %s 2>&1" %
                       (my_file, sharedpath, caseroot, file_build),
                       from_dir=exeroot,
                       ok_to_fail=True, verbose=True)[0]
        expect(stat == 0, "ERROR: buildlib.%s failed, cat %s" % (lib, file_build))
        logs.append(file_build)

    return logs

###############################################################################
def _build_model_thread(config_dir, caseroot, bldroot, compspec, file_build,
                        exeroot, model, comp, objdir, incroot, thread_bad_results):
###############################################################################
    stat = run_cmd("%s/buildlib %s %s %s >> %s 2>&1" %
                   (config_dir, caseroot, bldroot, compspec, file_build),
                   from_dir=os.path.join(exeroot, model), ok_to_fail=True)[0]
    if (stat != 0):
        thread_bad_results.append("ERROR: %s.buildlib failed, see %s" % (comp, file_build))

    for mod_file in glob.glob(os.path.join(objdir, "*_[Cc][Oo][Mm][Pp]_*.mod")):
        shutil.copy(mod_file, incroot)
