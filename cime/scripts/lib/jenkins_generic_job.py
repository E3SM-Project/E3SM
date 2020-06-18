import CIME.wait_for_tests
from CIME.utils import expect, run_cmd_no_fail
from CIME.case import Case

import os, shutil, glob, signal, logging, threading, sys, re, tarfile, time

##############################################################################
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
def delete_old_test_data(mach_comp, test_id_root, scratch_root, test_root, run_area, build_area, archive_area, avoid_test_id):
###############################################################################
    # Remove old dirs
    for clutter_area in [scratch_root, test_root, run_area, build_area, archive_area]:
        for old_file in glob.glob("{}/*{}*{}*".format(clutter_area, mach_comp, test_id_root)):
            if avoid_test_id not in old_file:
                logging.info("TEST ARCHIVER: Removing {}".format(old_file))
                if (os.path.isdir(old_file)):
                    shutil.rmtree(old_file)
                else:
                    os.remove(old_file)

###############################################################################
def scan_for_test_ids(old_test_archive, mach_comp, test_id_root):
###############################################################################
    results = set([])
    test_id_re = re.compile(".+[.]([^.]+)")
    for item in glob.glob("{}/{}/*{}*{}*".format(old_test_archive, "old_cases", mach_comp, test_id_root)):
        filename = os.path.basename(item)
        the_match = test_id_re.match(filename)
        if the_match:
            test_id = the_match.groups()[0]
            results.add(test_id)

    return list(results)

###############################################################################
def archive_old_test_data(machine, mach_comp, test_id_root, scratch_root, test_root, old_test_archive, avoid_test_id):
###############################################################################

    gb_allowed    = machine.get_value("MAX_GB_OLD_TEST_DATA")
    gb_allowed    = 500 if gb_allowed is None else gb_allowed
    bytes_allowed = gb_allowed * 1000000000
    expect(bytes_allowed > 0, "Machine {} does not support test archiving".format(machine.get_machine_name()))

    # Remove old cs.status, cs.submit. I don't think there's any value to leaving these around
    # or archiving them
    for old_cs_file in glob.glob("{}/cs.*".format(scratch_root)):
        if avoid_test_id not in old_cs_file:
            logging.info("TEST ARCHIVER: Removing {}".format(old_cs_file))
            os.remove(old_cs_file)

    # Remove the old CTest XML, same reason as above
    if (os.path.isdir("Testing")):
        logging.info("TEST ARCHIVER: Removing {}".format(os.path.join(os.getcwd(), "Testing")))
        shutil.rmtree("Testing")

    if not os.path.exists(old_test_archive):
        os.mkdir(old_test_archive)

    # Archive old data by looking at old test cases
    for old_case in glob.glob("{}/*{}*{}*".format(test_root, mach_comp, test_id_root)):
        if avoid_test_id not in old_case:
            logging.info("TEST ARCHIVER: archiving case {}".format(old_case))
            exeroot, rundir, archdir = run_cmd_no_fail("./xmlquery EXEROOT RUNDIR DOUT_S_ROOT --value", from_dir=old_case).split(",")

            for the_dir, target_area in [(exeroot, "old_builds"), (rundir, "old_runs"), (archdir, "old_archives"), (old_case, "old_cases")]:
                if os.path.exists(the_dir):
                    start_time = time.time()
                    logging.info("TEST ARCHIVER:   archiving {} to {}".format(the_dir, os.path.join(old_test_archive, target_area)))
                    if not os.path.exists(os.path.join(old_test_archive, target_area)):
                        os.mkdir(os.path.join(old_test_archive, target_area))

                    old_case_name = os.path.basename(old_case)
                    with tarfile.open(os.path.join(old_test_archive, target_area, "{}.tar.gz".format(old_case_name)), "w:gz") as tfd:
                        tfd.add(the_dir, arcname=old_case_name)

                    shutil.rmtree(the_dir)

                    # Remove parent dir if it's empty
                    parent_dir = os.path.dirname(the_dir)
                    if not os.listdir(parent_dir) or os.listdir(parent_dir) == ["case2_output_root"]:
                        shutil.rmtree(parent_dir)

                    end_time = time.time()
                    logging.info("TEST ARCHIVER:   archiving {} took {} seconds".format(the_dir, int(end_time - start_time)))

    # Check size of archive
    bytes_of_old_test_data = int(run_cmd_no_fail("du -sb {}".format(old_test_archive)).split()[0])
    if bytes_of_old_test_data > bytes_allowed:
        logging.info("TEST ARCHIVER: Too much test data, {}GB (actual) > {}GB (limit)".format(bytes_of_old_test_data / 1000000000, bytes_allowed / 1000000000))
        old_test_ids = scan_for_test_ids(old_test_archive, mach_comp, test_id_root)
        for old_test_id in sorted(old_test_ids):
            logging.info("TEST ARCHIVER:   Removing old data for test {}".format(old_test_id))
            for item in ["old_cases", "old_builds", "old_runs", "old_archives"]:
                for dir_to_rm in glob.glob("{}/{}/*{}*{}*".format(old_test_archive, item, mach_comp, old_test_id)):
                    logging.info("TEST ARCHIVER:     Removing {}".format(dir_to_rm))
                    if (os.path.isdir(dir_to_rm)):
                        shutil.rmtree(dir_to_rm)
                    else:
                        os.remove(dir_to_rm)

            bytes_of_old_test_data = int(run_cmd_no_fail("du -sb {}".format(old_test_archive)).split()[0])
            if bytes_of_old_test_data < bytes_allowed:
                break

    else:
        logging.info("TEST ARCHIVER: Test data is within accepted bounds, {}GB (actual) < {}GB (limit)".format(bytes_of_old_test_data / 1000000000, bytes_allowed / 1000000000))

###############################################################################
def handle_old_test_data(machine, compiler, test_id_root, scratch_root, test_root, avoid_test_id):
###############################################################################
    run_area = os.path.dirname(os.path.dirname(machine.get_value("RUNDIR"))) # Assumes XXX/$CASE/run
    build_area = os.path.dirname(os.path.dirname(machine.get_value("EXEROOT"))) # Assumes XXX/$CASE/build
    archive_area = os.path.dirname(machine.get_value("DOUT_S_ROOT")) # Assumes XXX/archive/$CASE
    old_test_archive = os.path.join(scratch_root, "old_test_archive")

    mach_comp = "{}_{}".format(machine.get_machine_name(), compiler)

    try:
        archive_old_test_data(machine, mach_comp, test_id_root, scratch_root, test_root, old_test_archive, avoid_test_id)
    except Exception:
        logging.warning("TEST ARCHIVER: Archiving of old test data FAILED: {}\nDeleting data instead".format(sys.exc_info()[1]))
        delete_old_test_data(mach_comp, test_id_root, scratch_root, test_root, run_area, build_area, archive_area, avoid_test_id)

###############################################################################
def jenkins_generic_job(generate_baselines, submit_to_cdash, no_batch,
                        baseline_name,
                        arg_cdash_build_name, cdash_project,
                        arg_test_suite,
                        cdash_build_group, baseline_compare,
                        scratch_root, parallel_jobs, walltime,
                        machine, compiler, real_baseline_name, baseline_root, update_success):
###############################################################################
    """
    Return True if all tests passed
    """
    use_batch = machine.has_batch_system() and not no_batch
    test_suite = machine.get_value("TESTS")
    proxy = machine.get_value("PROXY")
    test_suite = test_suite if arg_test_suite is None else arg_test_suite
    test_root = os.path.join(scratch_root, "J")

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

    test_id_root = "J{}{}".format(baseline_name.capitalize(), test_suite.replace("e3sm_", "").capitalize())
    test_id = "%s%s" % (test_id_root, CIME.utils.get_timestamp())
    archiver_thread = threading.Thread(target=handle_old_test_data, args=(machine, compiler, test_id_root, scratch_root, test_root, test_id))
    archiver_thread.start()

    #
    # Set up create_test command and run it
    #

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

    if baseline_root is not None:
        create_test_args.append(" --baseline-root " + baseline_root)

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
                                                      cdash_build_group=cdash_build_group,
                                                      update_success=update_success)

    logging.info("TEST ARCHIVER: Waiting for archiver thread")
    archiver_thread.join()
    logging.info("TEST ARCHIVER: Waiting for archiver finished")

    if use_batch and CIME.wait_for_tests.SIGNAL_RECEIVED:
        # Cleanup
        cleanup_queue(test_root, test_id)

    return tests_passed
