from utils import run_cmd, check_minimum_python_version, get_current_head, run_cmd_no_fail, get_current_commit
check_minimum_python_version(3, 4)

import os, shutil
import concurrent.futures as threading3

###############################################################################
class TestAllScream(object):
###############################################################################

    ###########################################################################
    def __init__(self, cxx, kokkos, submit, parallel, fast_fail, baseline, machine, custom_cmake_opts, tests):
    ###########################################################################

        self._cxx               = cxx
        self._kokkos            = kokkos
        self._submit            = submit
        self._parallel          = parallel 
        self._fast_fail         = fast_fail
        self._baseline          = baseline
        self._machine           = machine
        self._custom_cmake_opts = custom_cmake_opts
        self._tests             = tests
        self._src_dir           = os.getcwd()
        if not self._tests:
            self._tests = ["dbg", "sp", "fpe"]

        self._tests_cmake_args = {"dbg" : [("CMAKE_BUILD_TYPE", "Debug")],
                                  "sp"  : [("CMAKE_BUILD_TYPE", "Debug"),
                                           ("SCREAM_DOUBLE_PRECISION", "False")],
                                  "fpe" : [("CMAKE_BUILD_TYPE", "Debug"),
                                           ("SCREAM_PACK_SIZE", "1"),
                                           ("SCREAM_SMALL_PACK_SIZE", "1")]}

        self._test_full_names = { "dbg" : "full_debug",
                                  "sp"  : "full_sp_debug",
                                  "fpe" : "debug_nopack_fpe"}

        # Deduce how many resources per test
        self._proc_count = 4 # default
        if "CTEST_PARALLEL_LEVEL" in os.environ:
            try:
                self._proc_count = int(os.environ["CTEST_PARALLEL_LEVEL"])
            except ValueError:
                pass

        if self._parallel:
            # We need to be aware that other builds may be running too.
            # (Do not oversubscribe the machine)
            self._proc_count = self._proc_count // len(self._tests)

            # In case we have more tests than cores (unlikely)
            if self._proc_count == 0:
                self._proc_count = 1

    ###############################################################################
    def generate_cmake_config(self, extra_configs):
    ###############################################################################
        if self._kokkos:
            kokkos_cmake = "-DKokkos_DIR={}".format(self._kokkos)
        else:
            kokkos_cmake = "-C {}/cmake/machine-files/{}.cmake".format(self._src_dir,self._machine)

        result = "cmake -DCMAKE_CXX_COMPILER={} {}".format(self._cxx, kokkos_cmake)
        for key, value in extra_configs:
            result += " -D{}={}".format(key, value)

        if self._custom_cmake_opts:
            result += " {}".format(self._custom_cmake_opts)

        return result

    ###############################################################################
    def generate_ctest_config(self, cmake_config, extra_configs, test):
    ###############################################################################
        name = self._test_full_names[test]
        result = ""
        if self._submit:
            result += "CIME_MACHINE={} ".format(self._machine)

        result += "ctest -V --output-on-failure "

        if not self._submit:
            result += "-DNO_SUBMIT=True "

        for key, value in extra_configs:
            result += "-D{}={} ".format(key, value)

        result += "-DBUILD_NAME_MOD={} ".format(name)
        result += "-DBUILD_PARALLEL_RESOURCES={} ".format(self._proc_count)

        result += '-S {}/cmake/ctest_script.cmake -DCMAKE_COMMAND="{}" '.format(self._src_dir,cmake_config)

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
        stat, _, err = run_cmd("{} {}".format(cmake_config,self._src_dir),from_dir=test_dir, verbose=True)
        if stat == 0:

            stat, _, err = run_cmd("make -j{} && make baseline".format(self._proc_count), from_dir=test_dir, verbose=True)

        if stat != 0:
            print("WARNING: Failed to create baselines:\n{}".format(err))

        baseline_files = run_cmd_no_fail("find . -name ""*baseline*"" -type f").split()
        datas = []
        for baseline_file in baseline_files:
            if stat == 0:
                with open(baseline_file, "rb") as fd:
                    datas.append(fd.read())
            else:
                os.remove(baseline_file)

        # Clean baseline build files
        run_cmd_no_fail("/bin/rm -rf *",from_dir=test_dir) # Clean out baseline build

        return baseline_files, datas

    ###############################################################################
    def generate_all_baselines(self,git_baseline_head,git_head):
    ###############################################################################
        print("Generating baselines for ref {}".format(git_baseline_head))

        run_cmd_no_fail("git checkout {}".format(git_baseline_head))
        print("  Switched to {} ({})".format(git_baseline_head,get_current_commit()))

        success = True
        num_workers = len(self._tests) if self._parallel else 1
        with threading3.ProcessPoolExecutor(max_workers=num_workers) as executor:

            future_to_test = {
                executor.submit(self.generate_baselines,test) : test
                for test in self._tests}

            for future in threading3.as_completed(future_to_test):
                test = future_to_test[future]
                try:
                    baseline_files, datas = future.result()
                    for filepath, data in zip(baseline_files, datas):
                        if not os.path.isdir(os.path.dirname(filepath)):
                            os.makedirs(os.path.dirname(filepath))
                        with open(filepath, "wb") as fd: 
                            fd.write(data)
                except RuntimeError:
                    print('Generation of baselines for build {} failed'.format(self._test_full_names[test]))
                    success &= False

        print("  Switched back to {} ({})".format(git_head,get_current_commit()))

        return success

    ###############################################################################
    def run_test(self, test):
    ###############################################################################
        name = self._test_full_names[test]
        git_head = get_current_head()
        print("Testing '{}' for build type '{}'".format(git_head,name))

        test_dir = "ctest-build/{}".format(name)
        cmake_config = self.generate_cmake_config(self._tests_cmake_args[test])
        ctest_config = self.generate_ctest_config(cmake_config, [], test)

        success = run_cmd(ctest_config, from_dir=test_dir, arg_stderr=None, verbose=True)[0] == 0

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
                    name = self._test_full_names[test]
                    print('Build type {} failed. Here''s a list of failed tests:'.format(name))
                    stat,out,err = run_cmd("cat ctest-build/{}/Testing/Temporary/LastTestsFailed*".format(name))
                    print(out.strip())

                    return success

        for t,s in tests_success.items():
            if not s:
                name = self._test_full_names[t]
                print('Build type {} failed. Here''s a list of failed tests:'.format(name))
                stat,out,err = run_cmd("cat ctest-build/{}/Testing/Temporary/LastTestsFailed*".format(name))
                print(out.strip())
        return success

    ###############################################################################
    def test_all_scream(self):
    ###############################################################################
        git_head_commit = get_current_commit()
        git_head = get_current_head()

        print("Testing git ref {} ({})".format(git_head, git_head_commit))

        # First, create build directories (one per test)
        for test in self._tests:
            # Get this test's build dir name and cmake args
            full_name = self._test_full_names[test]
            test_dir = "ctest-build/{}".format(full_name)

            # Create this test's build dir
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
            os.mkdir(test_dir)

        # Second, generate baselines
        git_baseline_commit = get_current_commit(commit=self._baseline)
        if git_baseline_commit == git_head_commit:
            self._baseline = None
            print("WARNING: baseline commit is same as current HEAD")

        git_baseline_head = "HEAD" if self._baseline is None else self._baseline
        success = self.generate_all_baselines(git_baseline_head,git_head)
        if not success:
            print ("Error(s) occurred during baselines generation phase")
            return success

        # Finally, run the tests
        success &= self.run_all_tests()
        if not success:
            print ("Error(s) occurred during test phase")

        return success
