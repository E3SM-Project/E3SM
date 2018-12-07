"""
functions for building CIME models
"""
import glob, shutil, time, threading, subprocess
from CIME.XML.standard_module_setup  import *
from CIME.utils                 import get_model, analyze_build_log, stringify_bool, run_and_log_case_status, get_timestamp, run_sub_or_cmd, run_cmd, get_batch_script_for_job, gzip_existing_file, safe_copy
from CIME.provenance            import save_build_provenance as save_build_provenance_sub
from CIME.locked_files          import lock_file, unlock_file

logger = logging.getLogger(__name__)

###############################################################################
def _build_model(build_threaded, exeroot, clm_config_opts, incroot, complist,
                 lid, caseroot, cimeroot, compiler, buildlist, comp_interface):
###############################################################################
    logs = []

    thread_bad_results = []
    for model, comp, nthrds, _, config_dir in complist:
        if buildlist is not None and model.lower() not in buildlist:
            continue

        # aquap has a dependency on atm so we will build it after the threaded loop
        if comp == "aquap":
            logger.debug("Skip aquap ocn build here")
            continue

        # coupler handled seperately
        if model == "cpl":
            continue

        # special case for clm
        # clm 4_0 is not a shared (as in sharedlibs, shared by all tests) library and must be built here
        # clm 4_5 and newer is a shared library (but not in E3SM) and should be built in build_libraries
        if get_model() != "e3sm":
            if comp == "clm":
                if "clm4_0" in clm_config_opts:
                    logger.info("         - Building clm4_0 Library ")
                else:
                    continue
        else:
            logger.info("         - Building {} Library ".format(model))

        smp = nthrds > 1 or build_threaded

        bldroot = os.path.join(exeroot, model, "obj")
        libroot = os.path.join(exeroot, "lib")
        file_build = os.path.join(exeroot, "{}.bldlog.{}".format(model, lid))
        logger.debug("bldroot is {}".format(bldroot))
        logger.debug("libroot is {}".format(libroot))

        # make sure bldroot and libroot exist
        for build_dir in [bldroot, libroot]:
            if not os.path.exists(build_dir):
                os.makedirs(build_dir)

        # build the component library
        # thread_bad_results captures error output from thread (expected to be empty)
        # logs is a list of log files to be compressed and added to the case logs/bld directory
        t = threading.Thread(target=_build_model_thread,
            args=(config_dir, model, comp, caseroot, libroot, bldroot, incroot, file_build,
                  thread_bad_results, smp, compiler))
        t.start()

        logs.append(file_build)

    # Wait for threads to finish
    while(threading.active_count() > 1):
        time.sleep(1)

    expect(not thread_bad_results, "\n".join(thread_bad_results))

    #
    # Now build the executable
    #

    if not buildlist:
        cime_model = get_model()
        file_build = os.path.join(exeroot, "{}.bldlog.{}".format(cime_model, lid))

        config_dir = os.path.join(cimeroot, "src", "drivers", comp_interface, "cime_config")
        bldroot = os.path.join(exeroot, "cpl", "obj")
        if not os.path.isdir(bldroot):
            os.makedirs(bldroot)
        logger.info("Building {} with output to {} ".format(cime_model, file_build))

        with open(file_build, "w") as fd:
            stat = run_cmd("{}/buildexe {} {} {} "
                       .format(config_dir, caseroot, libroot, bldroot),
                       from_dir=bldroot,  arg_stdout=fd,
                       arg_stderr=subprocess.STDOUT)[0]

        analyze_build_log("{} exe".format(cime_model), file_build, compiler)
        expect(stat == 0, "BUILD FAIL: buildexe failed, cat {}".format(file_build))

        # Copy the just-built ${MODEL}.exe to ${MODEL}.exe.$LID
        safe_copy("{}/{}.exe".format(exeroot, cime_model), "{}/{}.exe.{}".format(exeroot, cime_model, lid))

        logs.append(file_build)

    return logs

###############################################################################
def _build_checks(case, build_threaded, comp_interface, use_esmf_lib,
                  debug, compiler, mpilib, complist, ninst_build, smp_value,
                  model_only, buildlist):
###############################################################################
    """
    check if a build needs to be done and warn if a clean is warrented first
    returns the relative sharedpath directory for sharedlibraries
    """
    ninst_value  = case.get_value("NINST_VALUE")
    smp_build    = case.get_value("SMP_BUILD")
    build_status = case.get_value("BUILD_STATUS")
    expect(comp_interface in ("mct", "moab", "nuopc"),
           "Only supporting mct nuopc, or moab comp_interfaces at this time, found {}".format(comp_interface))
    smpstr = ""
    inststr = ""
    for model, _, nthrds, ninst, _ in complist:
        if nthrds > 1:
            build_threaded = True
        if build_threaded:
            smpstr += "{}1".format(model[0])
        else:
            smpstr += "{}0".format(model[0])
        inststr += "{}{:d}".format((model[0]),ninst)

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
    sharedpath = os.path.join(compiler, mpilib, debugdir, threaddir, comp_interface)

    logger.debug("compiler={} mpilib={} debugdir={} threaddir={}"
                 .format(compiler,mpilib,debugdir,threaddir))

    expect(ninst_build == ninst_value or ninst_build == "0",
            """
ERROR, NINST VALUES HAVE CHANGED
  NINST_BUILD = {}
  NINST_VALUE = {}
  A manual clean of your obj directories is strongly recommended
  You should execute the following:
    ./case.build --clean
  Then rerun the build script interactively
  ---- OR ----
  You can override this error message at your own risk by executing:
    ./xmlchange -file env_build.xml -id NINST_BUILD -val 0
  Then rerun the build script interactively
""".format(ninst_build, ninst_value))

    expect(smp_build == smpstr or smp_build == "0",
            """
ERROR, SMP VALUES HAVE CHANGED
  SMP_BUILD = {}
  SMP_VALUE = {}
  smpstr = {}
  A manual clean of your obj directories is strongly recommended
  You should execute the following:
    ./case.build --clean
  Then rerun the build script interactively
  ---- OR ----
  You can override this error message at your own risk by executing:
    ./xmlchange -file env_build.xml -id SMP_BUILD -val 0
  Then rerun the build script interactively
""".format(smp_build, smp_value, smpstr))

    expect(build_status == 0,
           """
ERROR env_build HAS CHANGED
  A manual clean of your obj directories is required
  You should execute the following:
    ./case.build --clean-all
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

    # User may have rm -rf their build directory
    case.create_dirs()

    case.flush()
    if not model_only and not buildlist:
        logger.info("Generating component namelists as part of build")
        case.create_namelists()

    return sharedpath

###############################################################################
def _build_libraries(case, exeroot, sharedpath, caseroot, cimeroot, libroot, lid, compiler, buildlist, comp_interface):
###############################################################################

    shared_lib = os.path.join(exeroot, sharedpath, "lib")
    shared_inc = os.path.join(exeroot, sharedpath, "include")
    for shared_item in [shared_lib, shared_inc]:
        if (not os.path.exists(shared_item)):
            os.makedirs(shared_item)
    mpilib = case.get_value("MPILIB")
    libs = ["gptl", "mct", "pio", "csm_share"]
    if mpilib == "mpi-serial":
        libs.insert(0, mpilib)
    logs = []
    sharedlibroot = os.path.abspath(case.get_value("SHAREDLIBROOT"))
    for lib in libs:
        if buildlist is not None and lib not in buildlist:
            continue

        if lib == "csm_share":
            # csm_share adds its own dir name
            full_lib_path = os.path.join(sharedlibroot, sharedpath)
        elif lib == "mpi-serial":
            full_lib_path = os.path.join(sharedlibroot, sharedpath, "mct", lib)
        else:
            full_lib_path = os.path.join(sharedlibroot, sharedpath, lib)
        # pio build creates its own directory
        if (lib != "pio" and not os.path.exists(full_lib_path)):
            os.makedirs(full_lib_path)

        file_build = os.path.join(exeroot, "{}.bldlog.{}".format(lib, lid))
        my_file = os.path.join(cimeroot, "src", "build_scripts", "buildlib.{}".format(lib))
        logger.info("Building {} with output to file {}".format(lib,file_build))

        run_sub_or_cmd(my_file, [full_lib_path, os.path.join(exeroot, sharedpath), caseroot], 'buildlib',
                       [full_lib_path, os.path.join(exeroot, sharedpath), case], logfile=file_build)

        analyze_build_log(lib, file_build, compiler)
        logs.append(file_build)
        if lib == "pio":
            bldlog = open(file_build, "r")
            for line in bldlog:
                if re.search("Current setting for", line):
                    logger.warning(line)

    # clm not a shared lib for E3SM
    if get_model() != "e3sm" and (buildlist is None or "lnd" in buildlist):
        comp_lnd = case.get_value("COMP_LND")
        clm_config_opts = case.get_value("CLM_CONFIG_OPTS")
        if comp_lnd == "clm" and "clm4_0" not in clm_config_opts:
            logging.info("         - Building clm4_5/clm5_0 Library ")
            esmfdir = "esmf" if case.get_value("USE_ESMF_LIB") else "noesmf"
            bldroot = os.path.join(sharedlibroot, sharedpath, comp_interface, esmfdir, "clm","obj" )
            libroot = os.path.join(exeroot, sharedpath, comp_interface, esmfdir, "lib")
            incroot = os.path.join(exeroot, sharedpath, comp_interface, esmfdir, "include")
            file_build = os.path.join(exeroot, "lnd.bldlog.{}".format( lid))
            config_lnd_dir = os.path.dirname(case.get_value("CONFIG_LND_FILE"))

            for ndir in [bldroot, libroot, incroot]:
                if (not os.path.isdir(ndir)):
                    os.makedirs(ndir)

            smp = "SMP" in os.environ and os.environ["SMP"] == "TRUE"
            # thread_bad_results captures error output from thread (expected to be empty)
            # logs is a list of log files to be compressed and added to the case logs/bld directory
            thread_bad_results = []
            _build_model_thread(config_lnd_dir, "lnd", comp_lnd, caseroot, libroot, bldroot, incroot,
                                file_build, thread_bad_results, smp, compiler)
            logs.append(file_build)
            expect(not thread_bad_results, "\n".join(thread_bad_results))

    return logs

###############################################################################
def _build_model_thread(config_dir, compclass, compname, caseroot, libroot, bldroot, incroot, file_build,
                        thread_bad_results, smp, compiler):
###############################################################################
    logger.info("Building {} with output to {}".format(compclass, file_build))
    t1 = time.time()
    cmd = os.path.join(caseroot, "SourceMods", "src." + compname, "buildlib")
    if os.path.isfile(cmd):
        logger.warning("WARNING: using local buildlib script for {}".format(compname))
    else:
        cmd = os.path.join(config_dir, "buildlib")
        expect(os.path.isfile(cmd), "Could not find buildlib for {}".format(compname))

    with open(file_build, "w") as fd:
        stat = run_cmd("MODEL={} SMP={} {} {} {} {} "
                       .format(compclass, stringify_bool(smp), cmd, caseroot, libroot, bldroot),
                       from_dir=bldroot,  arg_stdout=fd,
                       arg_stderr=subprocess.STDOUT)[0]
    analyze_build_log(compclass, file_build, compiler)
    if (stat != 0):
        thread_bad_results.append("BUILD FAIL: {}.buildlib failed, cat {}".format(compname, file_build))

    for mod_file in glob.glob(os.path.join(bldroot, "*_[Cc][Oo][Mm][Pp]_*.mod")):
        safe_copy(mod_file, incroot)

    t2 = time.time()
    logger.info("{} built in {:f} seconds".format(compname, (t2 - t1)))

###############################################################################
def _clean_impl(case, cleanlist, clean_all, clean_depends):
###############################################################################
    exeroot = os.path.abspath(case.get_value("EXEROOT"))
    if clean_all:
        # If cleanlist is empty just remove the bld directory
        expect(exeroot is not None,"No EXEROOT defined in case")
        if os.path.isdir(exeroot):
            logging.info("cleaning directory {}".format(exeroot))
            shutil.rmtree(exeroot)
        # if clean_all is True also remove the sharedlibpath
        sharedlibroot = os.path.abspath(case.get_value("SHAREDLIBROOT"))
        expect(sharedlibroot is not None,"No SHAREDLIBROOT defined in case")
        if sharedlibroot != exeroot and os.path.isdir(sharedlibroot):
            logging.warning("cleaning directory {}".format(sharedlibroot))
            shutil.rmtree(sharedlibroot)
    else:
        expect((cleanlist is not None and len(cleanlist) > 0) or
               (clean_depends is not None and len(clean_depends)),"Empty cleanlist not expected")
        debug           = case.get_value("DEBUG")
        use_esmf_lib    = case.get_value("USE_ESMF_LIB")
        build_threaded  = case.get_build_threaded()
        gmake           = case.get_value("GMAKE")
        caseroot        = os.path.abspath(case.get_value("CASEROOT"))
        casetools       = case.get_value("CASETOOLS")
        clm_config_opts = case.get_value("CLM_CONFIG_OPTS")

        os.environ["DEBUG"]           = stringify_bool(debug)
        os.environ["USE_ESMF_LIB"]    = stringify_bool(use_esmf_lib)
        os.environ["BUILD_THREADED"]  = stringify_bool(build_threaded)
        os.environ["CASEROOT"]        = caseroot
        os.environ["COMP_INTERFACE"]  = case.get_value("COMP_INTERFACE")
        os.environ["PIO_VERSION"]     = str(case.get_value("PIO_VERSION"))
        os.environ["CLM_CONFIG_OPTS"] = clm_config_opts  if clm_config_opts is not None else ""

        cmd = gmake + " -f " + os.path.join(casetools, "Makefile")
        if cleanlist is not None:
            for item in cleanlist:
                tcmd = cmd + " clean" + item
                logger.info("calling {} ".format(tcmd))
                run_cmd_no_fail(tcmd)
        else:
            for item in clean_depends:
                tcmd = cmd + " clean_depends" + item
                logger.info("calling {} ".format(tcmd))
                run_cmd_no_fail(tcmd)

    # unlink Locked files directory
    unlock_file("env_build.xml")

    # reset following values in xml files
    case.set_value("SMP_BUILD",str(0))
    case.set_value("NINST_BUILD",str(0))
    case.set_value("BUILD_STATUS",str(0))
    case.set_value("BUILD_COMPLETE","FALSE")
    case.flush()

###############################################################################
def _case_build_impl(caseroot, case, sharedlib_only, model_only, buildlist,
                     save_build_provenance):
###############################################################################

    t1 = time.time()

    expect(not (sharedlib_only and model_only),
           "Contradiction: both sharedlib_only and model_only")
    logger.info("Building case in directory {}".format(caseroot))
    logger.info("sharedlib_only is {}".format(sharedlib_only))
    logger.info("model_only is {}".format(model_only))

    expect(os.path.isdir(caseroot), "'{}' is not a valid directory".format(caseroot))
    os.chdir(caseroot)

    expect(os.path.exists(get_batch_script_for_job(case.get_primary_job())),
           "ERROR: must invoke case.setup script before calling build script ")

    cimeroot = case.get_value("CIMEROOT")

    comp_classes = case.get_values("COMP_CLASSES")

    case.check_lockedfiles(skip="env_batch")

    # Retrieve relevant case data
    # This environment variable gets set for cesm Make and
    # needs to be unset before building again.
    if "MODEL" in os.environ:
        del os.environ["MODEL"]
    build_threaded      = case.get_build_threaded()
    casetools           = case.get_value("CASETOOLS")
    exeroot             = os.path.abspath(case.get_value("EXEROOT"))
    incroot             = os.path.abspath(case.get_value("INCROOT"))
    libroot             = os.path.abspath(case.get_value("LIBROOT"))
    sharedlibroot       = os.path.abspath(case.get_value("SHAREDLIBROOT"))
    multi_driver = case.get_value("MULTI_DRIVER")
    complist = []
    ninst = 1
    for comp_class in comp_classes:
        if comp_class == "CPL":
            config_dir = None
            if multi_driver:
                ninst = case.get_value("NINST_MAX")
        else:
            config_dir = os.path.dirname(case.get_value("CONFIG_{}_FILE".format(comp_class)))
            if multi_driver:
                ninst = 1
            else:
                ninst = case.get_value("NINST_{}".format(comp_class))

        comp = case.get_value("COMP_{}".format(comp_class))
        thrds =  case.get_value("NTHRDS_{}".format(comp_class))
        expect(ninst is not None,"Failed to get ninst for comp_class {}".format(comp_class))
        complist.append((comp_class.lower(), comp, thrds, ninst, config_dir ))
        os.environ["COMP_{}".format(comp_class)] = comp

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
    mali_use_albany     = case.get_value("MALI_USE_ALBANY")
    use_moab            = case.get_value("USE_MOAB")
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
    os.environ["MALI_USE_ALBANY"]      = stringify_bool(mali_use_albany)
    os.environ["USE_MOAB"]             = stringify_bool(use_moab)

    if get_model() == "e3sm" and mach == "titan" and compiler == "pgiacc":
        case.set_value("CAM_TARGET", "preqx_acc")

    # This is a timestamp for the build , not the same as the testid,
    # and this case may not be a test anyway. For a production
    # experiment there may be many builds of the same case.
    lid               = get_timestamp("%y%m%d-%H%M%S")
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
    # For now, there is just the one MALI_USE_ALBANY variable, but in
    # the future there may be others -- so USE_ALBANY will be true if
    # ANY of those are true.

    use_albany = stringify_bool(mali_use_albany)
    case.set_value("USE_ALBANY", use_albany)
    os.environ["USE_ALBANY"] = use_albany

    # Load modules
    case.load_env()

    sharedpath = _build_checks(case, build_threaded, comp_interface,
                               use_esmf_lib, debug, compiler, mpilib,
                               complist, ninst_build, smp_value, model_only, buildlist)

    t2 = time.time()
    logs = []

    if not model_only:
        logs = _build_libraries(case, exeroot, sharedpath, caseroot,
                                cimeroot, libroot, lid, compiler, buildlist, comp_interface)

    if not sharedlib_only:
        os.environ["INSTALL_SHAREDPATH"] = os.path.join(exeroot, sharedpath) # for MPAS makefile generators
        logs.extend(_build_model(build_threaded, exeroot, clm_config_opts, incroot, complist,
                                lid, caseroot, cimeroot, compiler, buildlist, comp_interface))

        if not buildlist:
            # in case component build scripts updated the xml files, update the case object
            case.read_xml()
            # Note, doing buildlists will never result in the system thinking the build is complete

    post_build(case, logs, build_complete=not (buildlist or sharedlib_only),
               save_build_provenance=save_build_provenance)

    t3 = time.time()

    if not sharedlib_only:
        logger.info("Time spent not building: {:f} sec".format(t2 - t1))
        logger.info("Time spent building: {:f} sec".format(t3 - t2))
        logger.info("MODEL BUILD HAS FINISHED SUCCESSFULLY")

    return True

###############################################################################
def post_build(case, logs, build_complete=False, save_build_provenance=True):
###############################################################################
    for log in logs:
        gzip_existing_file(log)

    if build_complete:
        # must ensure there's an lid
        lid = os.environ["LID"] if "LID" in os.environ else get_timestamp("%y%m%d-%H%M%S")
        if save_build_provenance:
            save_build_provenance_sub(case, lid=lid)
        # Set XML to indicate build complete
        case.set_value("BUILD_COMPLETE", True)
        case.set_value("BUILD_STATUS", 0)
        if "SMP_VALUE" in os.environ:
            case.set_value("SMP_BUILD", os.environ["SMP_VALUE"])
            case.flush()

        lock_file("env_build.xml")


###############################################################################
def case_build(caseroot, case, sharedlib_only=False, model_only=False, buildlist=None, save_build_provenance=True):
###############################################################################
    functor = lambda: _case_build_impl(caseroot, case, sharedlib_only, model_only, buildlist,
                                       save_build_provenance)
    return run_and_log_case_status(functor, "case.build", caseroot=caseroot)

###############################################################################
def clean(case, cleanlist=None, clean_all=False, clean_depends=None):
###############################################################################
    functor = lambda: _clean_impl(case, cleanlist, clean_all, clean_depends)
    return run_and_log_case_status(functor, "build.clean", caseroot=case.get_value("CASEROOT"))
