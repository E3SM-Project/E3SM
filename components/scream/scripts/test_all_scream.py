from utils import run_cmd, run_cmd_no_fail, check_minimum_python_version, get_current_head,     \
    get_current_commit, get_current_branch, expect, is_repo_clean, cleanup_repo,  \
    get_common_ancestor, merge_git_ref, checkout_git_ref, print_last_commit

from machines_specs import get_mach_compilation_resources, get_mach_testing_resources, \
                           get_mach_baseline_root_dir, setup_mach_env

check_minimum_python_version(3, 4)

import os, shutil, pathlib
import concurrent.futures as threading3
import itertools

###############################################################################
class TestAllScream(object):
###############################################################################

    ###########################################################################
    def __init__(self, cxx, kokkos=None, submit=False, parallel=False, fast_fail=False, baseline_ref=None,
                 baseline_dir=None, machine=None, no_tests=False, keep_tree=False,
                 custom_cmake_opts=(), custom_env_vars=(), tests=(),
                 integration_test="JENKINS_HOME" in os.environ, root_dir=None, dry_run=False,
                 make_parallel_level=0, ctest_parallel_level=0):
    ###########################################################################

        self._cxx                     = cxx
        self._kokkos                  = kokkos
        self._submit                  = submit
        self._parallel                = parallel
        self._fast_fail               = fast_fail
        self._baseline_ref            = baseline_ref
        self._machine                 = machine
        self._perform_tests           = not no_tests
        self._keep_tree               = keep_tree
        self._baseline_dir            = baseline_dir
        self._custom_cmake_opts       = custom_cmake_opts
        self._custom_env_vars         = custom_env_vars
        self._tests                   = tests
        self._root_dir                = root_dir
        self._integration_test        = integration_test
        self._dry_run                 = dry_run
        self._must_generate_baselines = False
        self._testing_dir             = "ctest-build"

        ############################################
        #  Sanity checks and helper structs setup  #
        ############################################

        expect (not self._baseline_dir or self._testing_dir != self._baseline_dir,
                "Error! For your safety, do NOT use 'ctest-build' to store baselines. Move them to a different directory.")

        expect(not (self._baseline_ref and self._baseline_dir),
               "Makes no sense to specify a baseline generation commit if using pre-existing baselines ")

        self._tests_cmake_args = {"dbg" : [("CMAKE_BUILD_TYPE", "Debug")],
                                  "sp"  : [("CMAKE_BUILD_TYPE", "Debug"),
                                           ("SCREAM_DOUBLE_PRECISION", "False")],
                                  "fpe" : [("CMAKE_BUILD_TYPE", "Debug"),
                                           ("SCREAM_PACK_SIZE", "1"),
                                           ("SCREAM_SMALL_PACK_SIZE", "1")]}

        self._test_full_names = { "dbg" : "full_debug",
                                  "sp"  : "full_sp_debug",
                                  "fpe" : "debug_nopack_fpe"}

        if not self._tests:
            self._tests = ["dbg", "sp", "fpe"]
        else:
            for t in self._tests:
                expect(t in self._test_full_names,
                       "Requested test '{}' is not supported by test-all-scream, please choose from: {}".\
                           format(t, ", ".join(self._test_full_names.keys())))

        # Compute root dir
        if not self._root_dir:
            self._root_dir = pathlib.Path(__file__).resolve().parent.parent
        else:
            self._root_dir = pathlib.Path(self._root_dir).resolve()
            expect(self._root_dir.is_dir() and self._root_dir.parts()[-2:] == ('scream', 'components'),
                   "Bad root-dir '{}', should be: $scream_repo/components/scream".format(self._root_dir))

        os.chdir(str(self._root_dir)) # needed, or else every git command will need repo=root_dir
        expect(get_current_commit(), "Root dir: {}, does not appear to be a git repo".format(self._root_dir))

        self._original_branch = get_current_branch()
        self._original_commit = get_current_commit()

        if not self._kokkos:
            expect(self._machine, "If no kokkos provided, must provide machine name for internal kokkos build")
        if self._submit:
            expect(self._machine, "If dashboard submit request, must provide machine name")

        print_last_commit(git_ref=self._original_branch)

        ###################################
        #      Compute baseline info      #
        ###################################

        default_baselines_root_dir = pathlib.Path(self._testing_dir,"baselines")
        if self._baseline_dir is None:
            if self._baseline_ref is None:
                # Compute baseline ref
                if self._keep_tree:
                    self._baseline_ref = "HEAD"
                elif self._integration_test:
                    self._baseline_ref = "origin/master"
                    merge_git_ref(git_ref="origin/master",verbose=True)
                else:
                    self._baseline_ref = get_common_ancestor("origin/master")
                    # Prefer a symbolic ref if possible
                    if self._baseline_ref is None or self._baseline_ref == get_current_commit(commit="origin/master"):
                        self._baseline_ref = "origin/master"
            self._must_generate_baselines = True

            self._baseline_dir = pathlib.Path(default_baselines_root_dir).resolve()

        else:
            # We treat the "AUTO" string as a request for automatic baseline dir.
            if self._baseline_dir == "AUTO":
                self._baseline_dir = get_mach_baseline_root_dir(self._machine,default_baselines_root_dir)

            self._baseline_dir = pathlib.Path(self._baseline_dir).resolve();

            # Make sure the baseline root directory exists
            expect(self._baseline_dir.is_dir(), "Baseline_dir {} is not a dir".format(self._baseline_dir))

            if self._integration_test:
                self._baseline_ref = "origin/master"
                merge_git_ref(git_ref=self._baseline_ref,verbose=True)
            else:
                for test in self._tests:
                    test_baseline_dir = self.get_preexisting_baseline(test)
                    expect(test_baseline_dir.is_dir(), "Missing baseline {}".format(test_baseline_dir))

        # Name of the file used to store/check the git sha of the repo used to generate baselines,
        # and name of the file used to store/check the builds for which baselines are available
        # Store it once to avoid typos-like bugs
        self._baseline_sha_file = pathlib.Path(self._baseline_dir,"baseline_git_sha").resolve()
        self._baseline_names_file = pathlib.Path(self._baseline_dir,"baseline_names").resolve()

        if self._integration_test:
            master_sha = get_current_commit(commit=self._baseline_ref)
            if not self.baselines_are_present():
                print ("Some baselines were not found. Rebuilding them.")
                self._must_generate_baselines = True
            elif self.baselines_are_expired(expected_baseline_sha=master_sha):
                print ("Baselines expired. Rebuilding them.")
                self._must_generate_baselines = True
            else:
                print ("Baselines found and not expired. Skipping baselines generation.")

        if self._must_generate_baselines:
            print("Using commit {} to generate baselines".format(self._baseline_ref))

        ##################################################
        #   Deduce how many testing resources per test   #
        ##################################################

        if ctest_parallel_level > 0:
            ctest_max_jobs = ctest_parallel_level
            print("Note: honoring requested value for ctest parallel level: {}".format(ctest_max_jobs))
        else:
            ctest_max_jobs = get_mach_testing_resources(self._machine)
            print("Note: no value passed for --ctest-parallel-level. Using the default for this machine: {}".format(ctest_max_jobs))

        self._testing_res_count = {"dbg" : ctest_max_jobs,
                                   "sp"  : ctest_max_jobs,
                                   "fpe" : ctest_max_jobs}

        # Deduce how many compilation resources per test
        if make_parallel_level > 0:
            make_max_jobs = make_parallel_level
            print("Note: honoring requested value for make parallel level: {}".format(make_max_jobs))
        else:
            make_max_jobs = get_mach_compilation_resources(self._machine)
            print("Note: no value passed for --make-parallel-level. Using the default for this machine: {}".format(make_max_jobs))

        self._compile_res_count = {"dbg" : make_max_jobs,
                                   "sp"  : make_max_jobs,
                                   "fpe" : make_max_jobs}

        if self._parallel:
            # We need to be aware that other builds may be running too.
            # (Do not oversubscribe the machine)
            make_remainder = make_max_jobs % len(self._tests)
            make_count     = make_max_jobs // len(self._tests)
            ctest_remainder = ctest_max_jobs % len(self._tests)
            ctest_count     = ctest_max_jobs // len(self._tests)

            # In case we have more tests than cores (unlikely)
            if make_count == 0:
                make_count = 1
            if ctest_count == 0:
                ctest_count = 1

            for test in self._tests:
                self._compile_res_count[test] = make_count
                if self._tests.index(test)<make_remainder:
                    self._compile_res_count[test] = make_count + 1

                self._testing_res_count[test] = ctest_count
                if self._tests.index(test)<ctest_remainder:
                    self._testing_res_count[test] = ctest_count + 1

                print("test {} can use {} jobs to compile, and {} jobs for testing".format(test,self._compile_res_count[test],self._testing_res_count[test]))

        if self._keep_tree:
            expect(not is_repo_clean(silent=True), "Makes no sense to use --keep-tree when repo is clean")
            expect(not self._integration_test, "Should not be doing keep-tree with integration testing")
            print("WARNING! You have uncommitted changes in your repo.",
                  "         The PASS/FAIL status may depend on these changes",
                  "         so if you want to keep them, don't forget to create a commit.",sep="\n")
            if self._baseline_dir is None:
                # Make sure the baseline ref is HEAD
                expect(self._baseline_ref == "HEAD",
                       "The option --keep-tree is only available when testing against pre-built baselines "
                       "(--baseline-dir) or HEAD (-b HEAD)")
        else:
            expect(is_repo_clean(),
                   "Repo must be clean before running. If testing against HEAD or pre-built baselines, "
                   "you can pass `--keep-tree` to allow non-clean repo.")

    ###############################################################################
    def get_baseline_dir(self, test):
    ###############################################################################
        return pathlib.Path(self._baseline_dir, self._test_full_names[test])

    ###############################################################################
    def get_test_dir(self, test):
    ###############################################################################
        return pathlib.Path(self._root_dir, self._testing_dir, self._test_full_names[test])

    ###############################################################################
    def get_preexisting_baseline(self, test):
    ###############################################################################
        expect(self._baseline_dir is not None, "Cannot supply preexisting baseline without baseline dir")
        return pathlib.Path(self._baseline_dir, self._test_full_names[test], "data")

    ###############################################################################
    def baselines_are_present (self):
    ###############################################################################
        # Check that all baselines are present (one subdir for all values of self._tests)

        # Sanity check
        expect(self._baseline_dir is not None, "Error! This routine should only be called when testing against pre-existing baselines.")

        # Even if a single baseline is missing, we consider all the baselines not present
        for test in self._tests:
            test_baseline_dir = pathlib.Path(self._baseline_dir, self._test_full_names[test], "data")
            if not test_baseline_dir.is_dir():
                return False

        # Note: inside this script we don't know what kind of file should be in the baseline dirs.
        #       If the actual files are missing, some other part of the testing will crash.
        return True

    ###############################################################################
    def baselines_are_expired (self, expected_baseline_sha):
    ###############################################################################
        # Baselines are expired if either:
        #  2) there is no file in baseline_dir containing the sha of the baselines
        #  3) the baselines sha does not match the one passed to this function

        # Sanity check
        expect(self._baseline_dir is not None, "Error! This routine should only be called when testing against pre-existing baselines.")

        # The file specifying what baselines were built during last baselines generation msut be there
        if not self._baseline_names_file.exists():
            return True

        # It might happen that we generate baselines for all build types, then later on
        # for some reason we manually generate baselines for only one build type. The other
        # baselines will still be there, but may be expired. Therefore, we check the
        # baselines_names file, to see what baselines were built last time. If all the
        # baselines we need are there, then we're good
        valid_baselines = run_cmd_no_fail("cat {}".format(self._baseline_names_file.resolve()))
        for test in self._tests:
           if not test in valid_baselines:
                return True

        # No sha file => baselines expired
        if not self._baseline_sha_file.exists():
            return True

        # Different sha => baselines expired
        baseline_sha = run_cmd_no_fail("cat {}".format(self._baseline_sha_file))
        return expected_baseline_sha != baseline_sha

    ###############################################################################
    def get_machine_file(self):
    ###############################################################################
        expect(self._machine is not None, "Cannot get machine file without machine")
        return pathlib.Path(self._root_dir, "cmake", "machine-files", "{}.cmake".format(self._machine))

    ###############################################################################
    def generate_cmake_config(self, extra_configs, for_ctest=False):
    ###############################################################################
        if self._kokkos:
            kokkos_cmake = "-DKokkos_DIR={}".format(self._kokkos)
        else:
            kokkos_cmake = "-C {}".format(self.get_machine_file())

        result = "{}-DCMAKE_CXX_COMPILER={} {}".format("" if for_ctest else "cmake ", self._cxx, kokkos_cmake)
        for key, value in extra_configs:
            result += " -D{}={}".format(key, value)

        for custom_opt in self._custom_cmake_opts:
            if "=" in custom_opt:
                name, value = custom_opt.split("=", 1)
                # Some effort is needed to ensure quotes are perserved
                result += " -D{}='{}'".format(name, value)
            else:
                result += " -D{}".format(custom_opt)

        return result

    ###############################################################################
    def get_taskset_id(self, test):
    ###############################################################################
        # Note: we need to loop through the whole list, since the compile_res_count
        #       might not be the same for all test.

        it = itertools.takewhile(lambda name: name!=test, self._tests)
        offset = sum(self._compile_res_count[prevs] for prevs in it)

        start = offset
        end   = offset + self._compile_res_count[test] - 1

        return start, end

    ###############################################################################
    def generate_ctest_config(self, cmake_config, extra_configs, test):
    ###############################################################################
        name = self._test_full_names[test]
        result = ""
        if self._submit:
            result += "CIME_MACHINE={} ".format(self._machine)

        result += "CTEST_PARALLEL_LEVEL={} ctest -V --output-on-failure ".format(self._testing_res_count[test])

        if self._baseline_dir is not None:
            cmake_config += " -DSCREAM_TEST_DATA_DIR={}".format(self.get_preexisting_baseline(test))

        if not self._submit:
            result += "-DNO_SUBMIT=True "

        for key, value in extra_configs:
            result += "-D{}={} ".format(key, value)

        result += "-DBUILD_NAME_MOD={} ".format(name)
        result += '-S {}/cmake/ctest_script.cmake -DCMAKE_COMMAND="{}" '.format(self._root_dir, cmake_config)

        # Ctest can only competently manage test pinning across a single instance of ctest. For
        # multiple concurrent instances of ctest, we have to help it.
        if self._parallel:
            start, end = self.get_taskset_id(test)
            result = "taskset -c {}-{} sh -c '{}'".format(start,end,result)

        return result

    ###############################################################################
    def generate_baselines(self, test):
    ###############################################################################
        test_dir = self.get_baseline_dir(test)

        cmake_config = self.generate_cmake_config(self._tests_cmake_args[test])
        cmake_config += " -DSCREAM_BASELINES_ONLY=ON"

        print("===============================================================================")
        print("Generating baseline for test {} with config '{}'".format(self._test_full_names[test], cmake_config))
        print("===============================================================================")

        # We cannot just crash if we fail to generate baselines, since we would
        # not get a dashboard report if we did that. Instead, just ensure there is
        # no baseline file to compare against if there's a problem.
        stat, _, err = run_cmd("{} {}".format(cmake_config, self._root_dir), from_dir=test_dir, verbose=True, dry_run=self._dry_run)
        if stat!= 0:
            print ("WARNING: Failed to configure baselines:\n{}".format(err))
            return False

        cmd = "make -j{} && make -j{} baseline".format(self._compile_res_count[test],self._testing_res_count[test])
        if self._parallel:
            start, end = self.get_taskset_id(test)
            cmd = "taskset -c {}-{} sh -c '{}'".format(start,end,cmd)

        stat, _, err = run_cmd(cmd, from_dir=test_dir, verbose=True, dry_run=self._dry_run)

        if stat != 0:
            print("WARNING: Failed to create baselines:\n{}".format(err))
            return False

        return True

    ###############################################################################
    def generate_all_baselines(self):
    ###############################################################################
        git_head_ref = get_current_head()

        print("###############################################################################")
        print("Generating baselines for ref {}".format(self._baseline_ref))
        print("###############################################################################")

        # First, create build directories (one per test)
        for test in self._tests:
            test_dir = self.get_baseline_dir(test)

            # Create this test's build dir
            if test_dir.exists():
                shutil.rmtree(str(test_dir))

            test_dir.mkdir(parents=True)

        checkout_git_ref(self._baseline_ref, verbose=True)

        success = True
        num_workers = len(self._tests) if self._parallel else 1
        with threading3.ProcessPoolExecutor(max_workers=num_workers) as executor:

            future_to_test = {
                executor.submit(self.generate_baselines, test) : test
                for test in self._tests}

            for future in threading3.as_completed(future_to_test):
                test = future_to_test[future]
                success &= future.result()

                if not success and self._fast_fail:
                    print('Generation of baselines for build {} failed'.format(self._test_full_names[test]))
                    return False

        if success:
            # Store the sha used for baselines generation
            run_cmd_no_fail("echo '{}' > {}".format(get_current_commit(commit=self._baseline_ref),self._baseline_sha_file))
            # Store the name of the builds for which we created a baseline
            tmp_string = ""
            for test in self._tests:
               tmp_string += " {}".format(test) 
            run_cmd_no_fail("echo '{}' > {}".format(tmp_string,self._baseline_names_file))

        checkout_git_ref(git_head_ref, verbose=True)

        return success

    ###############################################################################
    def run_test(self, test):
    ###############################################################################
        git_head = get_current_head()

        print("===============================================================================")
        print("Testing '{}' for test '{}'".format(git_head, self._test_full_names[test]))
        print("===============================================================================")

        test_dir = self.get_test_dir(test)
        cmake_config = self.generate_cmake_config(self._tests_cmake_args[test], for_ctest=True)
        ctest_config = self.generate_ctest_config(cmake_config, [], test)

        # This directory might have been used also to build the model to generate baselines.
        # Although it's ok to build in the same dir, we MUST make sure to erase cmake's cache
        # and internal files from the previous build (CMakeCache.txt and CMakeFiles folder)
        run_cmd_no_fail("rm -rf CMake*", from_dir=test_dir)

        success = run_cmd(ctest_config, from_dir=test_dir, arg_stdout=None, arg_stderr=None, verbose=True, dry_run=self._dry_run)[0] == 0

        return success

    ###############################################################################
    def run_all_tests(self):
    ###############################################################################
        print("###############################################################################")
        print("Running tests!")
        print("###############################################################################")

        # First, create build directories (one per test)
        for test in self._tests:
            test_dir = self.get_test_dir(test)

            # Create this test's build dir
            if test_dir.exists():
                shutil.rmtree(str(test_dir))

            test_dir.mkdir(parents=True)

        success = True
        tests_success = {
            test : False
            for test in self._tests}

        num_workers = len(self._tests) if self._parallel else 1
        with threading3.ProcessPoolExecutor(max_workers=num_workers) as executor:
            future_to_test = {
                executor.submit(self.run_test,test) : test
                for test in self._tests}

            for future in threading3.as_completed(future_to_test):
                test = future_to_test[future]
                tests_success[test] = future.result()
                success &= tests_success[test]
                # If failed, and fast fail is requested, return immediately
                # Note: this is effective only if num_worksers=1
                if not success and self._fast_fail:
                    break

        for t,s in tests_success.items():
            if not s:
                print("Build type {} failed. Here's a list of failed tests:".format(self._test_full_names[t]))
                test_dir = self.get_test_dir(t)
                test_results_dir = pathlib.Path(test_dir, "Testing", "Temporary")
                for result in test_results_dir.glob("LastTestsFailed*"):
                    print(result.read_text())

        return success

    ###############################################################################
    def test_all_scream(self):
    ###############################################################################

        # Setup the env on this machine
        setup_mach_env(self._machine)

        # Add any override the user may have requested
        for env_var in self._custom_env_vars:
            key,val = env_var.split("=",2)
            os.environ.update( { key : val } )

        success = True
        try:
            # If needed, generate baselines first
            if self._must_generate_baselines:
                expect(self._baseline_ref is not None, "Missing baseline ref")

                success = self.generate_all_baselines()
                if not success:
                    print ("Error(s) occurred during baselines generation phase")
                    return False

            # If requested, run tests
            if self._perform_tests:
                success &= self.run_all_tests()
                if not success:
                    print ("Error(s) occurred during test phase")

        finally:
            if not self._keep_tree:
                # Cleanup the repo if needed
                cleanup_repo(self._original_branch, self._original_commit)

        return success
