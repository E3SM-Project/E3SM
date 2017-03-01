import CIME.wait_for_tests
from CIME.utils import expect
from CIME.XML.machines import Machines

import os, shutil, glob, signal, logging

_MACHINE = Machines()

###############################################################################
def cleanup_queue(set_of_jobs_we_created):
###############################################################################
    """
    Delete all jobs left in the queue
    """
    current_jobs = set(CIME.utils.get_my_queued_jobs())
    jobs_to_delete = set_of_jobs_we_created & current_jobs

    if (jobs_to_delete):
        logging.warning("Found leftover batch jobs that need to be deleted: %s" % ", ".join(jobs_to_delete))
        success = CIME.utils.delete_jobs(jobs_to_delete)
        if not success:
            logging.warning("FAILED to clean up leftover jobs!")

###############################################################################
def jenkins_generic_job(generate_baselines, submit_to_cdash, no_batch,
                        baseline_name,
                        arg_cdash_build_name, cdash_project,
                        arg_test_suite,
                        cdash_build_group, baseline_compare,
                        scratch_root, parallel_jobs):
###############################################################################
    """
    Return True if all tests passed
    """
    use_batch = _MACHINE.has_batch_system() and not no_batch
    compiler = _MACHINE.get_default_compiler()
    test_suite = _MACHINE.get_value("TESTS")
    proxy = _MACHINE.get_value("PROXY")
    test_suite = test_suite if arg_test_suite is None else arg_test_suite
    test_root = os.path.join(scratch_root, "jenkins")

    if (use_batch):
        batch_system = _MACHINE.get_value("BATCH_SYSTEM")
        expect(batch_system is not None, "Bad XML. Batch machine has no batch_system configuration.")

    #
    # Env changes
    #

    if (submit_to_cdash and proxy is not None):
        os.environ["http_proxy"] = proxy

    #
    # Update submodules (Jenkins is struggling with this at the moment)
    #

    CIME.utils.run_cmd_no_fail("git submodule update --init", from_dir=CIME.utils.get_cime_root())

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

    # Remove the old CTest XML
    if (os.path.isdir("Testing")):
        shutil.rmtree("Testing")

    # Remove the old build/run dirs
    test_id_root = "jenkins_%s" % baseline_name
    for old_dir in glob.glob("%s/*%s*" % (scratch_root, test_id_root)):
        shutil.rmtree(old_dir)

    # Remove the old cases
    for old_file in glob.glob("%s/*%s*" % (test_root, test_id_root)):
        if (os.path.isdir(old_file)):
            shutil.rmtree(old_file)
        else:
            os.remove(old_file)

    #
    # Make note of things already in the queue so we know not to delete
    # them if we timeout
    #
    preexisting_queued_jobs = []
    if (use_batch):
        preexisting_queued_jobs = CIME.utils.get_my_queued_jobs()

    #
    # Set up create_test command and run it
    #

    baseline_args = ""
    if (generate_baselines):
        baseline_args = "-g -b %s" % baseline_name
    elif (baseline_compare == "yes"):
        baseline_args = "-c -b %s" % baseline_name

    batch_args = "--no-batch" if no_batch else ""
    pjob_arg = "" if parallel_jobs is None else "-j %d" % parallel_jobs

    test_id = "%s_%s" % (test_id_root, CIME.utils.get_timestamp())
    create_test_cmd = "./create_test %s --test-root %s -t %s %s %s %s" % \
                      (test_suite, test_root, test_id, baseline_args, batch_args, pjob_arg)

    if (not CIME.wait_for_tests.SIGNAL_RECEIVED):
        create_test_stat = CIME.utils.run_cmd(create_test_cmd, from_dir=CIME.utils.get_scripts_root(),
                                             verbose=True, arg_stdout=None, arg_stderr=None)[0]
        # Create_test should have either passed, detected failing tests, or timed out
        expect(create_test_stat in [0, CIME.utils.TESTS_FAILED_ERR_CODE, -signal.SIGTERM],
               "Create_test script FAILED with error code '%d'!" % create_test_stat)

    if (use_batch):
        # This is not fullproof. Any jobs that happened to be
        # submitted by this user while create_test was running will be
        # potentially deleted. This is still a big improvement over the
        # previous implementation which just assumed all queued jobs for this
        # user came from create_test.
        # TODO: change this to probe test_root for jobs ids
        #
        our_jobs = set(CIME.utils.get_my_queued_jobs()) - set(preexisting_queued_jobs)

    #
    # Wait for tests
    #

    if (submit_to_cdash):
        cdash_build_name = "_".join([test_suite, baseline_name, compiler]) if arg_cdash_build_name is None else arg_cdash_build_name
    else:
        cdash_build_name = None

    tests_passed = CIME.wait_for_tests.wait_for_tests(glob.glob("%s/*%s/TestStatus" % (test_root, test_id)),
                                                 no_wait=not use_batch, # wait if using queue
                                                 check_throughput=False, # don't check throughput
                                                 check_memory=False, # don't check memory
                                                 ignore_namelists=False, # don't ignore namelist diffs
                                                 cdash_build_name=cdash_build_name,
                                                 cdash_project=cdash_project,
                                                 cdash_build_group=cdash_build_group)
    if (not tests_passed and use_batch and CIME.wait_for_tests.SIGNAL_RECEIVED):
        # Cleanup
        cleanup_queue(our_jobs)

    return tests_passed
