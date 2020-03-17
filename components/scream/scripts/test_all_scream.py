from utils import run_cmd, check_minimum_python_version, run_cmd_no_fail, get_current_head,     \
                  get_current_commit, get_current_branch, expect, is_repo_clean, cleanup_repo,  \
                  get_common_ancestor, merge_git_ref, checkout_git_ref, print_last_commit

check_minimum_python_version(3, 4)

import os, shutil
import concurrent.futures as threading3

###############################################################################
class TestAllScream(object):
###############################################################################

    ###########################################################################
    def __init__(self, cxx, kokkos=None, submit=False, parallel=False, fast_fail=False, baseline_ref=None,
                 baseline_dir=None, machine=None, no_tests=False, keep_tree=False, custom_cmake_opts=(), tests=(),
                 integration_test="JENKINS_HOME" in os.environ):
    ###########################################################################

        self._cxx               = cxx
        self._kokkos            = kokkos
        self._submit            = submit
        self._parallel          = parallel
        self._fast_fail         = fast_fail
        self._baseline_ref      = baseline_ref
        self._machine           = machine
        self._perform_tests     = not no_tests
        self._keep_tree         = keep_tree
        self._baseline_dir      = baseline_dir
        self._custom_cmake_opts = custom_cmake_opts
        self._tests             = tests
        self._src_dir           = os.getcwd()
        self._integration_test  = integration_test
        self._original_branch   = get_current_branch()
        self._original_commit   = get_current_commit()

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

        expect(self._src_dir.endswith("components/scream"), "Run from $scream_repo/components/scream")
        if not self._kokkos:
            expect(self._machine, "If no kokkos provided, must provide machine name for internal kokkos build")
        if self._submit:
            expect(self._machine, "If dashboard submit request, must provide machine name")

        print_last_commit()

        # Compute baseline info
        expect(not (self._baseline_ref and self._baseline_dir),
               "Makes no sense to specify a baseline generation commit if using pre-existing baselines ")
        if self._baseline_dir is None:
            if self._baseline_ref is None:
                # Compute baseline ref
                if self._keep_tree:
                    self._baseline_ref = "HEAD"
                elif self._integration_test:
                    self._baseline_ref = "origin/master"
                    merge_git_ref(git_ref="origin/master")
                else:
                    self._baseline_ref = get_common_ancestor("origin/master")
                    # Prefer a symbolic ref if possible
                    if self._baseline_ref is None or self._baseline_ref == get_current_commit(commit="origin/master"):
                        self._baseline_ref = "origin/master"

                print("Using baseline commit {}".format(self._baseline_ref))
        else:
            if self._integration_test:
                if get_current_commit() != get_current_commit(commit="origin/master"):
                    merge_git_ref(git_ref="origin/master")

            print("NOTE: baselines for each build type BT must be in '{}/BT/data'. We don't check this, "
                  "but there will be errors if the baselines are not found.".format(self._baseline_dir))

        # Deduce how many resources per test
        self._proc_count = 4 # default
        proc_set = False
        if "CTEST_PARALLEL_LEVEL" in os.environ:
            try:
                self._proc_count = int(os.environ["CTEST_PARALLEL_LEVEL"])
                proc_set = True
            except ValueError:
                pass

        if not proc_set:
            print("WARNING: env var CTEST_PARALLEL_LEVEL unset, defaulting to {} which probably underutilizes your machine".\
                      format(self._proc_count))

        if self._parallel:
            # We need to be aware that other builds may be running too.
            # (Do not oversubscribe the machine)
            self._proc_count = self._proc_count // len(self._tests)

            # In case we have more tests than cores (unlikely)
            if self._proc_count == 0:
                self._proc_count = 1

        if self._keep_tree:
            expect(not is_repo_clean(), "Makes no sense to use --keep-tree when repo is clean")
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
    def generate_cmake_config(self, extra_configs, for_ctest=False):
    ###############################################################################
        if self._kokkos:
            kokkos_cmake = "-DKokkos_DIR={}".format(self._kokkos)
        else:
            kokkos_cmake = "-C {}/cmake/machine-files/{}.cmake".format(self._src_dir, self._machine)

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
        myid = self._tests.index(test)
        start = myid * self._proc_count
        end   = (myid+1) * self._proc_count - 1

        return start, end

    ###############################################################################
    def generate_ctest_config(self, cmake_config, extra_configs, test):
    ###############################################################################
        name = self._test_full_names[test]
        result = ""
        if self._submit:
            result += "CIME_MACHINE={} ".format(self._machine)

        result += "CTEST_PARALLEL_LEVEL={} ctest -V --output-on-failure ".format(self._proc_count)

        if self._baseline_dir is not None:
            cmake_config += " -DSCREAM_TEST_DATA_DIR={}/{}/data".format(self._baseline_dir,name)

        if not self._submit:
            result += "-DNO_SUBMIT=True "

        for key, value in extra_configs:
            result += "-D{}={} ".format(key, value)

        result += "-DBUILD_NAME_MOD={} ".format(name)
        result += '-S {}/cmake/ctest_script.cmake -DCMAKE_COMMAND="{}" '.format(self._src_dir, cmake_config)

        # Ctest can only competently manage test pinning across a single instance of ctest. For
        # multiple concurrent instances of ctest, we have to help it.
        if self._parallel:
            start, end = self.get_taskset_id(test)
            result = "taskset -c {}-{} sh -c '{}'".format(start,end,result)

        return result

    ###############################################################################
    def generate_baselines(self, test):
    ###############################################################################
        name = self._test_full_names[test]
        test_dir = "ctest-build/{}".format(name)

        cmake_config = self.generate_cmake_config(self._tests_cmake_args[test])

        print("Generating baseline for build type {} with config '{}'".format(name, cmake_config))

        # We cannot just crash if we fail to generate baselines, since we would
        # not get a dashboard report if we did that. Instead, just ensure there is
        # no baseline file to compare against if there's a problem.
        stat, _, err = run_cmd("{} {}".format(cmake_config, self._src_dir), from_dir=test_dir, verbose=True)
        if stat!= 0:
            print ("WARNING: Failed to configure baselines:\n{}".format(err))
            return False

        cmd = "make -j{} && make baseline".format(self._proc_count)
        if self._parallel:
            start, end = self.get_taskset_id(test)
            cmd = "taskset -c {}-{} sh -c '{}'".format(start,end,cmd)

        stat, _, err = run_cmd(cmd, from_dir=test_dir, verbose=True)

        if stat != 0:
            print("WARNING: Failed to create baselines:\n{}".format(err))
            return False

        return True

    ###############################################################################
    def generate_all_baselines(self):
    ###############################################################################
        git_head_commit     = get_current_commit()
        git_head_ref        = get_current_head()
        git_baseline_commit = get_current_commit(commit=self._baseline_ref)

        print("Generating baselines for ref {}".format(self._baseline_ref))

        checkout_git_ref(git_ref=self._baseline_ref,verbose=True)

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

        checkout_git_ref(git_ref=git_head_ref,verbose=True)

        return success

    ###############################################################################
    def run_test(self, test):
    ###############################################################################
        name = self._test_full_names[test]
        git_head = get_current_head()
        print("Testing '{}' for build type '{}'".format(git_head,name))

        test_dir = "ctest-build/{}".format(name)
        cmake_config = self.generate_cmake_config(self._tests_cmake_args[test], for_ctest=True)
        ctest_config = self.generate_ctest_config(cmake_config, [], test)

        success = run_cmd(ctest_config, from_dir=test_dir, arg_stdout=None, arg_stderr=None, verbose=True)[0] == 0

        return success

    ###############################################################################
    def run_all_tests(self):
    ###############################################################################
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
                name = self._test_full_names[t]
                print("Build type {} failed. Here's a list of failed tests:".format(name))
                out = run_cmd("cat ctest-build/{}/Testing/Temporary/LastTestsFailed*".format(name))[1]
                print(out)

        return success

    ###############################################################################
    def test_all_scream(self):
    ###############################################################################
        success = True
        try:
            # First, create build directories (one per test)
            for test in self._tests:
                # Get this test's build dir name and cmake args
                full_name = self._test_full_names[test]
                test_dir = "./ctest-build/{}".format(full_name)

                # Create this test's build dir
                if os.path.exists(test_dir):
                    shutil.rmtree(test_dir)

                os.makedirs(test_dir)

            if self._baseline_dir is None:
                # Second, generate baselines
                expect(self._baseline_ref is not None, "Missing baseline ref")

                success = self.generate_all_baselines()
                if not success:
                    print ("Error(s) occurred during baselines generation phase")

            if self._perform_tests:
                # Finally, run the tests
                success &= self.run_all_tests()
                if not success:
                    print ("Error(s) occurred during test phase")

            if not self._keep_tree:
                # Cleanup the repo if needed
                cleanup_repo(self._original_head, self._original_commit)

        except:
            if not self._keep_tree:
                # Cleanup the repo if needed, then re-throw whatever the exception was
                cleanup_repo(self._original_head, self._original_commit)
            raise

        return success
