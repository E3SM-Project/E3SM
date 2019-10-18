from utils import run_cmd, check_minimum_python_version, get_current_head, run_cmd_no_fail, get_current_commit
check_minimum_python_version(3, 4)

import os, shutil

###############################################################################
class TestAllScream(object):
###############################################################################

    ###########################################################################
    def __init__(self, cxx, kokkos, submit, baseline, machine, custom_cmake_opts, tests):
    ###########################################################################
        self._cxx               = cxx
        self._kokkos            = kokkos
        self._submit            = submit
        self._baseline          = baseline
        self._machine           = machine
        self._custom_cmake_opts = custom_cmake_opts
        self._tests             = tests

    ###############################################################################
    def generate_cmake_config(self, extra_configs):
    ###############################################################################
        if self._kokkos:
            kokkos_cmake = "-DKokkos_DIR={}".format(self._kokkos)
        else:
            kokkos_cmake = "-C ../cmake/machine-files/{}.cmake".format(self._machine)

        result = "cmake -DCMAKE_CXX_COMPILER={} {}".format(self._cxx, kokkos_cmake)
        for key, value in extra_configs:
            result += " -D{}={}".format(key, value)

        if self._custom_cmake_opts:
            result += " {}".format(self._custom_cmake_opts)

        return result

    ###############################################################################
    def generate_ctest_config(self, cmake_config, extra_configs, build_name):
    ###############################################################################
        result = ""
        if self._submit:
            result += "CIME_MACHINE={} ".format(self._machine)

        result += "ctest -V --output-on-failure "

        if not self._submit:
            result += "-DNO_SUBMIT=True "

        for key, value in extra_configs:
            result += "-D{}={} ".format(key, value)

        result += "-DBUILD_NAME_MOD=_{} ".format(build_name)

        result += '-S ../cmake/ctest_script.cmake -DCMAKE_COMMAND="{}" '.format(cmake_config)

        return result

    ###############################################################################
    def generate_baselines(self, cmake_config, git_head):
    ###############################################################################
        print("Generating baseline for {} with config '{}'".format("HEAD" if self._baseline is None else self._baseline, cmake_config))

        if self._baseline is not None:
            run_cmd_no_fail("git checkout {}".format(self._baseline))
            print("  Switched to {}".format(get_current_commit()))

        # We cannot just crash if we fail to generate baselines, since we would
        # not get a dashboard report if we did that. Instead, just ensure there is
        # no baseline file to compare against if there's a problem.
        stat, _, err = run_cmd("{} ..".format(cmake_config), verbose=True)
        if stat == 0:
            proc_count = None
            if "CTEST_PARALLEL_LEVEL" in os.environ:
                try:
                    proc_count = int(os.environ["CTEST_PARALLEL_LEVEL"])
                except ValueError:
                    pass

            if not proc_count:
                proc_count = 8 # Default

            stat, _, err = run_cmd("make -j{} && make baseline".format(proc_count), verbose=True)

        if stat != 0:
            print("WARNING: Failed to create baselines:\n{}".format(err))

        baseline_files = run_cmd_no_fail('find . -name "*baseline*" -type f').split()
        datas = []
        for baseline_file in baseline_files:
            if stat == 0:
                with open(baseline_file, "rb") as fd:
                    datas.append(fd.read())
            else:
                os.remove(baseline_file)

        if self._baseline is not None:
            run_cmd_no_fail("git checkout {}".format(git_head))
            print("  Switched back to {}".format(get_current_commit()))

        return baseline_files, datas

    ###############################################################################
    def run_test(self, extra_cmake_configs, extra_ctest_configs, build_name, git_head):
    ###############################################################################
        cmake_config = self.generate_cmake_config(extra_cmake_configs)

        # Clean out whatever might have been left in the build area from
        # previous tests
        run_cmd_no_fail("/bin/rm -rf *")

        if ("BUILD_ONLY", "True") not in extra_ctest_configs:
            filepaths, datas = self.generate_baselines(cmake_config, git_head)
            run_cmd_no_fail("/bin/rm -rf *") # Clean out baseline build
            for filepath, data in zip(filepaths, datas):
                if not os.path.isdir(os.path.dirname(filepath)):
                    os.makedirs(os.path.dirname(filepath))
                with open(filepath, "wb") as fd:
                    fd.write(data)

        ctest_config = self.generate_ctest_config(cmake_config, extra_ctest_configs, build_name)

        return run_cmd(ctest_config, arg_stdout=None, arg_stderr=None, verbose=True)[0] == 0

    ###############################################################################
    def test_all_scream(self):
    ###############################################################################
        git_head_commit = get_current_commit()
        git_baseline_commit = get_current_commit(commit=self._baseline)
        git_head = get_current_head()
        print("Testing {} ({})".format(git_head, git_head_commit))
        if git_baseline_commit == git_head_commit:
            self._baseline = None
            print("WARNING: baseline commit is same as current HEAD")

        success = True
        try:
            if os.path.exists("ctest-build"):
                shutil.rmtree("ctest-build")

            os.mkdir("ctest-build")
            os.chdir("ctest-build")

            # A full debug test
            if not self._tests or "dbg" in self._tests:
                success &= self.run_test([("CMAKE_BUILD_TYPE", "Debug")],
                                         [], "full_debug", git_head)

            # A full debug single precision
            if not self._tests or "sp" in self._tests:
                success &= self.run_test([("CMAKE_BUILD_TYPE", "Debug"), ("SCREAM_DOUBLE_PRECISION", "False")],
                                         [], "full_sp_debug", git_head)

            # A full debug test with packsize=1 and FPE
            if not self._tests or "fpe" in self._tests:
                if self._machine not in ["waterman", "white"]:
                    success &= self.run_test([("CMAKE_BUILD_TYPE", "Debug"), ("SCREAM_PACK_SIZE", "1"), ("SCREAM_SMALL_PACK_SIZE", "1")],
                                             [], "debug_nopack_fpe", git_head)

        finally:
            if self._baseline is not None:
                run_cmd_no_fail("git checkout {}".format(git_head))

        return success
