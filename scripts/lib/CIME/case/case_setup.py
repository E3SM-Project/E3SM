"""
Library for case.setup.
case_setup is a member of class Case from file case.py
"""

from CIME.XML.standard_module_setup import *

from CIME.XML.machines      import Machines
from CIME.BuildTools.configure import configure
from CIME.utils             import get_cime_root, run_and_log_case_status, get_model, get_batch_script_for_job, safe_copy
from CIME.test_status       import *
from CIME.locked_files      import unlock_file, lock_file

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
        model_file = case.get_value("CONFIG_{}_FILE".format(model))
    expect(model_file is not None,
           "Could not locate CONFIG_{}_FILE in config_files.xml".format(model))
    model_dir = os.path.dirname(model_file)

    expect(os.path.isdir(model_dir),
           "cannot find cime_config directory {} for component {}".format(model_dir, comp))
    comp_interface = case.get_value("COMP_INTERFACE")
    multi_driver = case.get_value("MULTI_DRIVER")
    ninst = 1

    if multi_driver:
        ninst_max = case.get_value("NINST_MAX")
        if comp_interface != "nuopc" and model not in ("DRV","CPL","ESP"):
            ninst_model = case.get_value("NINST_{}".format(model))
            expect(ninst_model==ninst_max,"MULTI_DRIVER mode, all components must have same NINST value.  NINST_{} != {}".format(model,ninst_max))
    if comp == "cpl":
        if not os.path.exists("user_nl_cpl"):
            safe_copy(os.path.join(model_dir, "user_nl_cpl"), ".")
    else:
        if comp_interface == "nuopc":
            ninst = case.get_value("NINST")
        elif ninst == 1:
            ninst = case.get_value("NINST_{}".format(model))
        nlfile = "user_nl_{}".format(comp)
        model_nl = os.path.join(model_dir, nlfile)
        if ninst > 1:
            for inst_counter in range(1, ninst+1):
                inst_nlfile = "{}_{:04d}".format(nlfile, inst_counter)
                if not os.path.exists(inst_nlfile):
                    # If there is a user_nl_foo in the case directory, copy it
                    # to user_nl_foo_INST; otherwise, copy the original
                    # user_nl_foo from model_dir
                    if os.path.exists(nlfile):
                        safe_copy(nlfile, inst_nlfile)
                    elif os.path.exists(model_nl):
                        safe_copy(model_nl, inst_nlfile)
        else:
            # ninst = 1
            if not os.path.exists(nlfile):
                if os.path.exists(model_nl):
                    safe_copy(model_nl, nlfile)

###############################################################################
def _case_setup_impl(case, caseroot, clean=False, test_mode=False, reset=False):
###############################################################################
    os.chdir(caseroot)

    non_local = case.get_value("NONLOCAL")

    # Check that $DIN_LOC_ROOT exists - and abort if not a namelist compare tests
    if not non_local:
        din_loc_root = case.get_value("DIN_LOC_ROOT")
        testcase     = case.get_value("TESTCASE")
        expect(not (not os.path.isdir(din_loc_root) and testcase != "SBN"),
               "inputdata root is not a directory: {}".format(din_loc_root))

    # Remove batch scripts
    if reset or clean:
        # clean batch script
        batch_script = get_batch_script_for_job(case.get_primary_job())
        if os.path.exists(batch_script):
            os.remove(batch_script)
            logger.info("Successfully cleaned batch script {}".format(batch_script))

        if not test_mode:
            # rebuild the models (even on restart)
            case.set_value("BUILD_COMPLETE", False)

    if not clean:
        if not non_local:
            case.load_env()

        models = case.get_values("COMP_CLASSES")
        mach = case.get_value("MACH")
        compiler = case.get_value("COMPILER")
        debug = case.get_value("DEBUG")
        mpilib = case.get_value("MPILIB")
        sysos = case.get_value("OS")
        comp_interface = case.get_value("COMP_INTERFACE")
        expect(mach is not None, "xml variable MACH is not set")

        # creates the Macros.make, Depends.compiler, Depends.machine, Depends.machine.compiler
        # and env_mach_specific.xml if they don't already exist.
        if not os.path.isfile("Macros.make") or not os.path.isfile("env_mach_specific.xml"):
            configure(Machines(machine=mach), caseroot, ["Makefile"], compiler, mpilib, debug, comp_interface, sysos)

        # Set tasks to 1 if mpi-serial library
        if mpilib == "mpi-serial":
            for vid, value in case:
                if vid.startswith("NTASKS") and value != 1:
                    case.set_value(vid, 1)

        # Check ninst.
        # In CIME there can be multiple instances of each component model (an ensemble) NINST is the instance of that component.
        comp_interface = case.get_value("COMP_INTERFACE")
        if comp_interface == "nuopc":
            ninst  = case.get_value("NINST")

        multi_driver = case.get_value("MULTI_DRIVER")

        for comp in models:
            ntasks = case.get_value("NTASKS_{}".format(comp))
            if comp == "CPL":
                continue
            if comp_interface != "nuopc":
                ninst  = case.get_value("NINST_{}".format(comp))
            if multi_driver:
                if comp_interface != "nuopc":
                    expect(case.get_value("NINST_LAYOUT_{}".format(comp)) == "concurrent",
                           "If multi_driver is TRUE, NINST_LAYOUT_{} must be concurrent".format(comp))
                case.set_value("NTASKS_PER_INST_{}".format(comp), ntasks)
            else:
                if ninst > ntasks:
                    if ntasks == 1:
                        case.set_value("NTASKS_{}".format(comp), ninst)
                        ntasks = ninst
                    else:
                        expect(False, "NINST_{comp} value {ninst} greater than NTASKS_{comp} {ntasks}".format(comp=comp, ninst=ninst, ntasks=ntasks))

                case.set_value("NTASKS_PER_INST_{}".format(comp), max(1,int(ntasks / ninst)))

        if os.path.exists(get_batch_script_for_job(case.get_primary_job())):
            logger.info("Machine/Decomp/Pes configuration has already been done ...skipping")

            case.initialize_derived_attributes()

            case.set_value("SMP_PRESENT", case.get_build_threaded())

        else:
            case.check_pelayouts_require_rebuild(models)

            unlock_file("env_build.xml")
            unlock_file("env_batch.xml")

            case.flush()
            case.check_lockedfiles()

            case.initialize_derived_attributes()

            cost_per_node = case.get_value("COSTPES_PER_NODE")
            case.set_value("COST_PES", case.num_nodes * cost_per_node)
            threaded = case.get_build_threaded()
            case.set_value("SMP_PRESENT", threaded)
            if threaded and case.total_tasks * case.thread_count > cost_per_node:
                smt_factor = max(1.0,int(case.get_value("MAX_TASKS_PER_NODE") / cost_per_node))
                case.set_value("TOTALPES", int(case.total_tasks * max(1.0,float(case.thread_count) / smt_factor)))
            else:
                case.set_value("TOTALPES", case.total_tasks*case.thread_count)

            # May need to select new batch settings if pelayout changed (e.g. problem is now too big for prev-selected queue)
            env_batch = case.get_env("batch")
            env_batch.set_job_defaults([(case.get_primary_job(), {})], case)

            # create batch files
            env_batch.make_all_batch_files(case)
            if get_model() == "e3sm" and not case.get_value("TEST"):
                input_batch_script = os.path.join(case.get_value("MACHDIR"), "template.case.run.sh")
                env_batch.make_batch_script(input_batch_script, "case.run", case, outfile=get_batch_script_for_job("case.run.sh"))

            # Make a copy of env_mach_pes.xml in order to be able
            # to check that it does not change once case.setup is invoked
            case.flush()
            logger.debug("at copy TOTALPES = {}".format(case.get_value("TOTALPES")))
            lock_file("env_mach_pes.xml")
            lock_file("env_batch.xml")

        # Create user_nl files for the required number of instances
        if not os.path.exists("user_nl_cpl"):
            logger.info("Creating user_nl_xxx files for components and cpl")

        # loop over models
        for model in models:
            comp = case.get_value("COMP_{}".format(model))
            logger.debug("Building {} usernl files".format(model))
            _build_usernl_files(case, model, comp)
            if comp == "cism":
                glcroot = case.get_value("COMP_ROOT_DIR_GLC")
                run_cmd_no_fail("{}/cime_config/cism.template {}".format(glcroot, caseroot))

        _build_usernl_files(case, "drv", "cpl")

        # Create needed directories for case
        case.create_dirs()

        logger.info("If an old case build already exists, might want to run \'case.build --clean\' before building")

        # Some tests need namelists created here (ERP) - so do this if we are in test mode
        if (test_mode or get_model() == "e3sm") and not non_local:
            logger.info("Generating component namelists as part of setup")
            case.create_namelists()

        # Record env information
        env_module = case.get_env("mach_specific")
        env_module.make_env_mach_specific_file("sh", case)
        env_module.make_env_mach_specific_file("csh", case)
        if not non_local:
            env_module.save_all_env_info("software_environment.txt")

        logger.info("You can now run './preview_run' to get more info on how your case will be run")

###############################################################################
def case_setup(self, clean=False, test_mode=False, reset=False):
###############################################################################
    caseroot, casebaseid = self.get_value("CASEROOT"), self.get_value("CASEBASEID")
    phase = "setup.clean" if clean else "case.setup"
    functor = lambda: _case_setup_impl(self, caseroot, clean, test_mode, reset)

    if self.get_value("TEST") and not test_mode:
        test_name = casebaseid if casebaseid is not None else self.get_value("CASE")
        with TestStatus(test_dir=caseroot, test_name=test_name) as ts:
            try:
                run_and_log_case_status(functor, phase, caseroot=caseroot)
            except:
                ts.set_status(SETUP_PHASE, TEST_FAIL_STATUS)
                raise
            else:
                if clean:
                    ts.set_status(SETUP_PHASE, TEST_PEND_STATUS)
                else:
                    ts.set_status(SETUP_PHASE, TEST_PASS_STATUS)
    else:
        run_and_log_case_status(functor, phase, caseroot=caseroot)
