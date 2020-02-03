from utils import run_cmd, check_minimum_python_version, get_current_head, run_cmd_no_fail, get_current_commit, expect, is_repo_clean
check_minimum_python_version(3, 4)

import os, shutil
import concurrent.futures as threading3

###############################################################################
class TestAllScream(object):
###############################################################################

    ###########################################################################
    def __init__(self, cxx, kokkos, submit, parallel, fast_fail, baseline_ref, baseline_dir, machine, no_tests, keep_tree, custom_cmake_opts, tests):
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

        if not self._baseline_dir == "NONE":
            print ("Ignoring baseline ref {}, and using baselines in directory {} instead".format(self._baseline_ref, self._baseline_dir))
            print ("NOTE: baselines for each build type BT must be in '{}/BT/data'. We don't check this, but there will be errors if the baselines are not found.".format(self._baseline_dir))

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
            if self._baseline_dir=="NONE":
                # Make sure the baseline ref is HEAD
                expect(self._baseline_ref=="HEAD","The option --keep-tree is only available when testing against pre-built baselines (--baseline-dir) or HEAD (-b HEAD)")
        else:
            expect(is_repo_clean(),"Repo must be clean before running. If testing against HEAD or pre-built baselines, you can pass `--keep-tree` to allow non-clean repo.")
                

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

        if not self._baseline_dir == "NONE":
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
    def generate_baselines(self, test, cleanup):
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

        cmd = "make -j{} && make baseline".format(self._proc_count);
        if self._parallel:
            start, end = self.get_taskset_id(test)
            cmd = "taskset -c {}-{} sh -c '{}'".format(start,end,cmd)

        stat, _, err = run_cmd(cmd, from_dir=test_dir, verbose=True)

        if stat != 0:
            print("WARNING: Failed to create baselines:\n{}".format(err))
            return False

        if cleanup:        
            run_cmd_no_fail("ls | grep -v data | xargs rm -rf ", from_dir=test_dir)

        return True

    ###############################################################################
    def generate_all_baselines(self, git_baseline_head, git_head):
    ###############################################################################
        print("Generating baselines for ref {}".format(git_baseline_head))

        if git_baseline_head != "HEAD":
            expect(is_repo_clean(), "If baseline commit is not HEAD, then the repo must be clean before running")
            run_cmd_no_fail("git checkout {}".format(git_baseline_head))
            print("  Switched to {} ({})".format(git_baseline_head, get_current_commit()))


        cleanup = git_baseline_head != "HEAD"
        success = True
        num_workers = len(self._tests) if self._parallel else 1
        with threading3.ProcessPoolExecutor(max_workers=num_workers) as executor:

            future_to_test = {
                executor.submit(self.generate_baselines, test, cleanup) : test
                for test in self._tests}

            for future in threading3.as_completed(future_to_test):
                test = future_to_test[future]
                success &= future.result()

                if not success and self._fast_fail:
                    print('Generation of baselines for build {} failed'.format(self._test_full_names[test]))
                    return False

        if git_baseline_head != "HEAD":
            run_cmd_no_fail("git checkout {}".format(git_head))
            print("  Switched back to {} ({})".format(git_head, get_current_commit()))

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
        git_head_commit = get_current_commit()
        git_head = get_current_head()

        print("Testing git ref {} ({})".format(git_head, git_head_commit))

        success = True
        # First, create build directories (one per test)
        for test in self._tests:
            # Get this test's build dir name and cmake args
            full_name = self._test_full_names[test]
            test_dir = "./ctest-build/{}".format(full_name)

            # Create this test's build dir
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)

            os.makedirs(test_dir)

        if self._baseline_dir=="NONE":
            # Second, generate baselines
            git_baseline_commit = get_current_commit(commit=self._baseline_ref)
            if git_baseline_commit == git_head_commit:
                self._baseline_ref = None
                print("WARNING: baseline commit is same as current HEAD")

            git_baseline_head = "HEAD" if self._baseline_ref is None else self._baseline_ref
            success = self.generate_all_baselines(git_baseline_head,git_head)
            if not success:
                print ("Error(s) occurred during baselines generation phase")
                return success

        if self._perform_tests:
            # Finally, run the tests
            success &= self.run_all_tests()
            if not success:
                print ("Error(s) occurred during test phase")

        return success
