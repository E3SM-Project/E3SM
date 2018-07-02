import CIME.wait_for_tests
from CIME.utils import expect
from CIME.case import Case

import os, shutil, glob, signal, logging

###############################################################################
def cleanup_queue(test_root, test_id):
###############################################################################
    """
    Delete all jobs left in the queue
    """
    for teststatus_file in glob.iglob("{}/*{}*/TestStatus".format(test_root, test_id)):
        case_dir = os.path.dirname(teststatus_file)
        with Case(case_dir, read_only=True) as case:
            jobmap = case.get_job_info()
            jobkills = []
            for jobname, jobid in jobmap.items():
                logging.warning("Found leftover batch job {} ({}) that need to be deleted".format(jobid, jobname))
                jobkills.append(jobid)

            case.cancel_batch_jobs(jobkills)

###############################################################################
def jenkins_generic_job(generate_baselines, submit_to_cdash, no_batch,
                        baseline_name,
                        arg_cdash_build_name, cdash_project,
                        arg_test_suite,
                        cdash_build_group, baseline_compare,
                        scratch_root, parallel_jobs, walltime,
                        machine, compiler, real_baseline_name):
###############################################################################
    """
    Return True if all tests passed
    """
    use_batch = machine.has_batch_system() and not no_batch
    test_suite = machine.get_value("TESTS")
    proxy = machine.get_value("PROXY")
    test_suite = test_suite if arg_test_suite is None else arg_test_suite
    test_root = os.path.join(scratch_root, "J")
    run_area = os.path.dirname(os.path.dirname(machine.get_value("RUNDIR")))

    if (use_batch):
        batch_system = machine.get_value("BATCH_SYSTEM")
        expect(batch_system is not None, "Bad XML. Batch machine has no batch_system configuration.")

    #
    # Env changes
    #

    if (submit_to_cdash and proxy is not None):
        os.environ["http_proxy"] = proxy

    if (not os.path.isdir(scratch_root)):
        os.makedirs(scratch_root)

    # Important, need to set up signal handlers before we officially
    # kick off tests. We don't want this process getting killed outright
    # since it's critical that the cleanup in the finally block gets run
    CIME.wait_for_tests.set_up_signal_handlers()

    #
    # Clean up leftovers from previous run of jenkins_generic_job. This will
    # break the previous run of jenkins_generic_job if it's still running. Set up
    # the Jenkins jobs with timeouts to avoid this.
    #

    mach_comp = "{}_{}".format(machine.get_machine_name(), compiler)

    # Remove the old CTest XML
    if (os.path.isdir("Testing")):
        shutil.rmtree("Testing")

    # Remove old dirs
    test_id_root = "J{}{}".format(baseline_name.capitalize(), test_suite.replace("e3sm_", "").capitalize())
    for clutter_area in [scratch_root, test_root, run_area]:
        for old_file in glob.glob("{}/*{}*{}*".format(clutter_area, mach_comp, test_id_root)):
            if (os.path.isdir(old_file)):
                shutil.rmtree(old_file)
            else:
                os.remove(old_file)

    #
    # Set up create_test command and run it
    #

    test_id = "%s%s" % (test_id_root, CIME.utils.get_timestamp())
    create_test_args = [test_suite, "--test-root %s" % test_root, "-t %s" % test_id, "--machine %s" % machine.get_machine_name(), "--compiler %s" % compiler]
    if (generate_baselines):
        create_test_args.append("-g -b " + real_baseline_name)
    elif (baseline_compare):
        create_test_args.append("-c -b " + real_baseline_name)

    if scratch_root != machine.get_value("CIME_OUTPUT_ROOT"):
        create_test_args.append("--output-root=" + scratch_root)

    if no_batch:
        create_test_args.append("--no-batch")

    if parallel_jobs is not None:
        create_test_args.append("-j {:d}".format(parallel_jobs))

    if walltime is not None:
        create_test_args.append(" --walltime " + walltime)

    create_test_cmd = "./create_test " + " ".join(create_test_args)

    if (not CIME.wait_for_tests.SIGNAL_RECEIVED):
        create_test_stat = CIME.utils.run_cmd(create_test_cmd, from_dir=CIME.utils.get_scripts_root(),
                                             verbose=True, arg_stdout=None, arg_stderr=None)[0]
        # Create_test should have either passed, detected failing tests, or timed out
        expect(create_test_stat in [0, CIME.utils.TESTS_FAILED_ERR_CODE, -signal.SIGTERM],
               "Create_test script FAILED with error code '{:d}'!".format(create_test_stat))

    #
    # Wait for tests
    #

    if (submit_to_cdash):
        cdash_build_name = "_".join([test_suite, baseline_name, compiler]) if arg_cdash_build_name is None else arg_cdash_build_name
    else:
        cdash_build_name = None

    os.environ["CIME_MACHINE"] = machine.get_machine_name()

    if submit_to_cdash:
        logging.info("To resubmit to dashboard: wait_for_tests {}/*{}/TestStatus --no-wait -b {}".format(test_root, test_id, cdash_build_name))

    tests_passed = CIME.wait_for_tests.wait_for_tests(glob.glob("{}/*{}/TestStatus".format(test_root, test_id)),
                                                 no_wait=not use_batch, # wait if using queue
                                                 check_throughput=False, # don't check throughput
                                                 check_memory=False, # don't check memory
                                                 ignore_namelists=False, # don't ignore namelist diffs
                                                 cdash_build_name=cdash_build_name,
                                                 cdash_project=cdash_project,
                                                 cdash_build_group=cdash_build_group)

    if use_batch and CIME.wait_for_tests.SIGNAL_RECEIVED:
        # Cleanup
        cleanup_queue(test_root, test_id)

    return tests_passed
