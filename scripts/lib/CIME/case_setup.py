"""
Library for case.setup.
"""

from CIME.XML.standard_module_setup import *

from CIME.check_lockedfiles import *
from CIME.preview_namelists import create_dirs, create_namelists
from CIME.XML.env_mach_pes  import EnvMachPes
from CIME.XML.machines      import Machines
from CIME.BuildTools.configure import configure
from CIME.utils             import append_status, get_cime_root
from CIME.test_status       import *

import shutil

logger = logging.getLogger(__name__)

###############################################################################
def _build_usernl_files(case, model, comp):
###############################################################################
    """
    Create user_nl_xxx files, expects cwd is caseroot
    """
    model = model.upper()
    if model == "DRV":
        model_file = case.get_value("CONFIG_CPL_FILE")
    else:
        model_file = case.get_value("CONFIG_%s_FILE" % model)
    expect(model_file is not None,
           "Could not locate CONFIG_%s_FILE in config_files.xml"%model)
    model_dir = os.path.dirname(model_file)

    expect(os.path.isdir(model_dir),
           "cannot find cime_config directory %s for component %s" % (model_dir, comp))

    if comp == "cpl":
        if not os.path.exists("user_nl_cpl"):
            shutil.copy(os.path.join(model_dir, "user_nl_cpl"), ".")
    else:
        ninst = case.get_value("NINST_%s" % model)
        nlfile = "user_nl_%s" % comp
        model_nl = os.path.join(model_dir, nlfile)
        if ninst > 1:
            for inst_counter in xrange(1, ninst+1):
                inst_nlfile = "%s_%04d" % (nlfile, inst_counter)
                if not os.path.exists(inst_nlfile):
                    # If there is a user_nl_foo in the case directory, copy it
                    # to user_nl_foo_INST; otherwise, copy the original
                    # user_nl_foo from model_dir
                    if os.path.exists(nlfile):
                        shutil.copy(nlfile, inst_nlfile)
                    elif os.path.exists(model_nl):
                        shutil.copy(model_nl, inst_nlfile)
        else:
            # ninst = 1
            if not os.path.exists(nlfile):
                if os.path.exists(model_nl):
                    shutil.copy(model_nl, nlfile)

###############################################################################
def _case_setup_impl(case, caseroot, clean=False, test_mode=False, reset=False, adjust_pio=True):
###############################################################################
    os.chdir(caseroot)
    msg = "case.setup starting"
    append_status(msg, caseroot=caseroot, sfile="CaseStatus")

    cimeroot = get_cime_root(case)

    # Check that $DIN_LOC_ROOT exists - and abort if not a namelist compare tests
    din_loc_root = case.get_value("DIN_LOC_ROOT")
    testcase     = case.get_value("TESTCASE")
    expect(not (not os.path.isdir(din_loc_root) and testcase != "SBN"),
           "inputdata root is not a directory: \"$din_loc_root\" ")

    # Check that userdefine settings are specified before expanding variable
    for vid, value in case:
        expect(not (type(value) is str and "USERDEFINED_required_build" in value),
               "Parameter '%s' must be defined" % vid)

    # Create batch script
    if reset or clean:
        # back up relevant files
        if os.path.exists("case.run"):
            os.remove("case.run")

        # only do the following if are NOT in testmode
        if not test_mode:
            # rebuild the models (even on restart)
            case.set_value("BUILD_COMPLETE", False)

            # backup and then clean test script
            if os.path.exists("case.test"):
                os.remove("case.test")
                logger.info("Successfully cleaned test script case.test")

            if os.path.exists("case.testdriver"):
                os.remove("case.testdriver")
                logger.info("Successfully cleaned test script case.testdriver")

        logger.info("Successfully cleaned batch script case.run")

        msg = "case.setup clean complete"
        append_status(msg, caseroot=caseroot, sfile="CaseStatus")

    if not clean:
        case.load_env()

        models = case.get_values("COMP_CLASSES")
        mach = case.get_value("MACH")
        compiler = case.get_value("COMPILER")
        debug = case.get_value("DEBUG")
        mpilib = case.get_value("MPILIB")
        sysos = case.get_value("OS")
        expect(mach is not None, "xml variable MACH is not set")

        # creates the Macros.make, Depends.compiler, Depends.machine, Depends.machine.compiler
        # and env_mach_specific.xml if they don't already exist.
        if not os.path.isfile("Macros.make") or not os.path.isfile("env_mach_specific.xml"):
            configure(Machines(machine=mach), caseroot, ["Makefile"], compiler, mpilib, debug, sysos)

        # Set tasks to 1 if mpi-serial library
        if mpilib == "mpi-serial":
            for vid, value in case:
                if vid.startswith("NTASKS_") and value != 1:
                    case.set_value(vid, 1)

        # Check ninst.
        # In CIME there can be multiple instances of each component model (an ensemble) NINST is the instance of that component.
        for comp in models:
            if comp == "CPL":
                continue
            ninst  = case.get_value("NINST_%s" % comp)
            ntasks = case.get_value("NTASKS_%s" % comp)
            if ninst > ntasks:
                if ntasks == 1:
                    case.set_value("NTASKS_%s" % comp, ninst)
                else:
                    expect(False, "NINST_%s value %d greater than NTASKS_%s %d" % (comp, ninst, comp, ntasks))

        if os.path.exists("case.run"):
            logger.info("Machine/Decomp/Pes configuration has already been done ...skipping")
        else:
            check_pelayouts_require_rebuild(case, models)

            unlock_file("env_build.xml")

            case.flush()
            check_lockedfiles()
            env_mach_pes = case.get_env("mach_pes")
            pestot = env_mach_pes.get_total_tasks(models)
            logger.debug("at update TOTALPES = %s"%pestot)
            case.set_value("TOTALPES", pestot)
            thread_count = env_mach_pes.get_max_thread_count(models)
            if thread_count > 1:
                case.set_value("BUILD_THREADED", True)

            expect(not (case.get_value("BUILD_THREADED")  and compiler == "nag"),
                   "it is not possible to run with OpenMP if using the NAG Fortran compiler")
            cost_pes = env_mach_pes.get_cost_pes(pestot, thread_count, machine=case.get_value("MACH"))
            case.set_value("COST_PES", cost_pes)

            # create batch files
            logger.info("Creating batch script case.run")
            env_batch = case.get_env("batch")
            num_nodes = env_mach_pes.get_total_nodes(pestot, thread_count)
            tasks_per_node = env_mach_pes.get_tasks_per_node(pestot, thread_count)
            for job in env_batch.get_jobs():
                input_batch_script  = os.path.join(case.get_value("MACHDIR"), env_batch.get_value('template', subgroup=job))
                if job == "case.test" and testcase is not None and not test_mode:
                    logger.info("Writing %s script" % job)
                    testscript = os.path.join(cimeroot, "scripts", "Testing", "Testcases", "%s_script" % testcase)
                    # Short term fix to be removed when csh tests are removed
                    if not os.path.exists(testscript):
                        env_batch.make_batch_script(input_batch_script, job, case, pestot, tasks_per_node, num_nodes, thread_count)
                elif job != "case.test":
                    logger.info("Writing %s script from input template %s" % (job, input_batch_script))
                    env_batch.make_batch_script(input_batch_script, job, case, pestot, tasks_per_node, num_nodes, thread_count)

            # Make sure pio settings are consistant
            if adjust_pio:
                adjust_pio_layout(case, tasks_per_node)

            # Make a copy of env_mach_pes.xml in order to be able
            # to check that it does not change once case.setup is invoked
            logger.info("Locking file env_mach_pes.xml")
            case.flush()
            logger.debug("at copy TOTALPES = %s"%case.get_value("TOTALPES"))
            lock_file("env_mach_pes.xml")

        # Create user_nl files for the required number of instances
        if not os.path.exists("user_nl_cpl"):
            logger.info("Creating user_nl_xxx files for components and cpl")
        # loop over models
        for model in models:
            comp = case.get_value("COMP_%s" % model)
            logger.debug("Building %s usernl files"%model)
            _build_usernl_files(case, model, comp)
            if comp == "cism":
                run_cmd_no_fail("%s/../components/cism/cime_config/cism.template %s" % (cimeroot, caseroot))

        _build_usernl_files(case, "drv", "cpl")

        # Create needed directories for case
        create_dirs(case)

        logger.info("If an old case build already exists, might want to run \'case.build --clean\' before building")

        # Create test script if appropriate
        # Short term fix to be removed when csh tests are removed
        if os.path.exists("env_test.xml"):
            if not os.path.exists("case.test"):
                logger.info("Starting testcase.setup")
                run_cmd_no_fail("./testcase.setup -caseroot %s" % caseroot)
                logger.info("Finished testcase.setup")

        # Some tests need namelists created here (ERP) - so do this if are in test mode
        if test_mode:
            logger.info("Generating component namelists as part of setup")
            create_namelists(case)

        msg = "case.setup complete"
        append_status(msg, caseroot=caseroot, sfile="CaseStatus")

        # Record env information
        env_module = case.get_env("mach_specific")
        env_module.make_env_mach_specific_file(compiler, debug, mpilib, "sh")
        env_module.make_env_mach_specific_file(compiler, debug, mpilib, "csh")
        env_module.save_all_env_info("software_environment.txt")


def adjust_pio_layout(case, new_pio_stride):

    models = case.get_values("COMP_CLASSES")
    for comp in models:
        pio_stride = case.get_value("PIO_STRIDE_%s"%comp)
        pio_numtasks = case.get_value("PIO_NUMTASKS_%s"%comp)
        ntasks = case.get_value("NTASKS_%s"%comp)
        new_stride = min(ntasks, new_pio_stride)
        new_numtasks = max(1, ntasks//new_stride)
        if pio_stride != new_stride:
            logger.info("Resetting  PIO_STRIDE_%s to %s"%(comp, new_stride))
            case.set_value("PIO_STRIDE_%s"%comp, new_stride)
            if pio_numtasks != new_numtasks:
                logger.info("Resetting  PIO_NUMTASKS_%s to %s"%(comp, new_numtasks))
                case.set_value("PIO_NUMTASKS_%s"%comp, new_numtasks)


###############################################################################
def case_setup(case, clean=False, test_mode=False, reset=False, no_status=False, adjust_pio=True):
###############################################################################
    caseroot, casebaseid = case.get_value("CASEROOT"), case.get_value("CASEBASEID")
    if case.get_value("TEST") and not no_status:
        test_name = casebaseid if casebaseid is not None else case.get_value("CASE")
        with TestStatus(test_dir=caseroot, test_name=test_name) as ts:
            try:
                _case_setup_impl(case, caseroot, clean=clean, test_mode=test_mode, reset=reset, adjust_pio=adjust_pio)
            except:
                ts.set_status(SETUP_PHASE, TEST_FAIL_STATUS)
                raise
            else:
                ts.set_status(SETUP_PHASE, TEST_PASS_STATUS)
    else:
        _case_setup_impl(case, caseroot, clean=clean, test_mode=test_mode, reset=reset)
