from utils import run_cmd, run_cmd_no_fail, expect, check_minimum_python_version, ensure_psutil, \
    SharedArea
from git_utils import get_current_head, get_current_commit

from test_factory import create_tests, COV, CSR

from machines_specs import get_machine, setup_mach_env, is_machine_supported, \
                           logical_cores_per_physical_core, get_cpu_ids_from_slurm_env_var

check_minimum_python_version(3, 4)

import os, shutil
import concurrent.futures as threading3
import itertools
import json

ensure_psutil()
import psutil
import re

from pathlib import Path

###############################################################################
class TestAllScream(object):
###############################################################################

    ###########################################################################
    def __init__(self, cxx_compiler=None, f90_compiler=None, c_compiler=None,
                 submit=False, parallel=False, generate=False,
                 baseline_dir=None, machine=None, config_only=False,
                 custom_cmake_opts=(), custom_env_vars=(), preserve_env=False, tests=(),
                 local=False, root_dir=None, work_dir=None,
                 quick_rerun=False,quick_rerun_failed=False,
                 make_parallel_level=0, ctest_parallel_level=0,
                 extra_verbose=False, limit_test_regex=None, test_level="at", test_size=None):
    ###########################################################################

        # When using scripts-tests, we can't pass "-l" to test-all-scream,
        # but we can pass "-m local". So if machine="local", reset things
        # as if local=True and machine=None
        if machine=="local":
            local = True
            machine = None

        self._cxx_compiler            = cxx_compiler
        self._f90_compiler            = f90_compiler
        self._c_compiler              = c_compiler
        self._submit                  = submit
        self._parallel                = parallel
        self._config_only             = config_only
        self._baseline_dir            = baseline_dir
        self._custom_cmake_opts       = custom_cmake_opts
        self._custom_env_vars         = custom_env_vars
        self._preserve_env            = preserve_env
        self._root_dir                = root_dir
        self._work_dir                = None if work_dir is None else Path(work_dir)
        self._quick_rerun             = quick_rerun
        self._quick_rerun_failed      = quick_rerun_failed
        self._extra_verbose           = extra_verbose
        self._limit_test_regex        = limit_test_regex
        self._test_level              = test_level
        self._test_size               = test_size
        self._generate                = generate

        if self._quick_rerun_failed:
            self._quick_rerun = True

        ############################################
        #  Sanity checks and helper structs setup  #
        ############################################

        # Quick rerun skips config phase, and config-only runs only config. You can't ask for both...
        expect (not (self._quick_rerun and self._config_only),
                "Makes no sense to ask for --quick-rerun and --config-only at the same time")

        # Probe machine if none was specified
        if machine is not None:
            self._machine = get_machine(machine)
            expect (not local, "Specifying a machine while passing '-l,--local' is ambiguous.")
        else:
            # We could potentially integrate more with CIME here to do actual
            # nodename probing.
            if "SCREAM_MACHINE" in os.environ and is_machine_supported(os.environ["SCREAM_MACHINE"]):
                self._machine = get_machine(os.environ["SCREAM_MACHINE"])
            else:
                expect(local,
                       "test-all-scream requires either the machine arg (-m $machine) or the -l flag,"
                       "which makes it look for machine specs in '~/.cime/scream_mach_specs.py'.")
                self._machine = get_machine("local")

        # Compute root dir (where repo is) and work dir (where build/test will happen)
        if not self._root_dir:
            self._root_dir = Path(__file__).resolve().parent.parent
        else:
            self._root_dir = Path(self._root_dir).resolve()
            expect(self._root_dir.is_dir() and list(self._root_dir.parts)[-2:] == ["components","eamxx"],
                   f"Bad root-dir '{self._root_dir}', should end with: /components/eamxx")

        # Make our test objects! Change mem to default mem-check test for current platform
        if "mem" in tests:
            tests[tests.index("mem")] = "csm" if self._machine.gpu_arch=="cuda" else "valg"
        self._tests = create_tests(tests, self)

        if self._work_dir is not None:
            self._work_dir = Path(self._work_dir).absolute()
        else:
            self._work_dir = self._root_dir.absolute().joinpath("ctest-build")

        self._work_dir.mkdir(parents=True, exist_ok=True)

        os.chdir(str(self._root_dir)) # needed, or else every git command will need repo=root_dir

        ###################################
        #  Compilation/testing resources  #
        ###################################

        # Deduce how many compilation resources per test
        make_max_jobs = self._machine.num_bld_res
        if make_parallel_level > 0:
            expect(make_parallel_level <= make_max_jobs,
                   f"Requested make_parallel_level {make_parallel_level} is more than max available {make_max_jobs}")
            make_max_jobs = make_parallel_level
            print(f"Note: honoring requested value for make parallel level: {make_max_jobs}")
        else:
            print(f"Note: no value passed for --make-parallel-level. Using the default for this machine: {make_max_jobs}")

        ctest_max_jobs = self._machine.num_run_res
        if ctest_parallel_level > 0:
            expect(ctest_parallel_level <= ctest_max_jobs,
                   f"Requested ctest_parallel_level {ctest_parallel_level} is more than max available {ctest_max_jobs}")
            ctest_max_jobs = ctest_parallel_level
            print(f"Note: honoring requested value for ctest parallel level: {ctest_max_jobs}")
        elif "CTEST_PARALLEL_LEVEL" in os.environ:
            env_val = int(os.environ["CTEST_PARALLEL_LEVEL"])
            expect(env_val <= ctest_max_jobs,
                   f"CTEST_PARALLEL_LEVEL env {env_val} is more than max available {ctest_max_jobs}")
            ctest_max_jobs = env_val
            print(f"Note: honoring environment value for ctest parallel level: {ctest_max_jobs}")
        else:
            print(f"Note: no value passed for --ctest-parallel-level. Using the default for this machine: {ctest_max_jobs}")

        self._ctest_max_jobs = ctest_max_jobs

        for test in self._tests:
            test.testing_res_count = ctest_max_jobs
            test.compile_res_count = make_max_jobs

        if self._parallel:
            # We need to be aware that other builds may be running too.
            # (Do not oversubscribe the machine)
            log_per_phys = logical_cores_per_physical_core()

            # Avoid splitting physical cores across test types
            make_jobs_per_test  = ((make_max_jobs  // len(self._tests)) // log_per_phys) * log_per_phys
            if self._machine.uses_gpu():
                ctest_jobs_per_test = ctest_max_jobs // len(self._tests)
            else:
                ctest_jobs_per_test = ((ctest_max_jobs // len(self._tests)) // log_per_phys) * log_per_phys

            # The current system of selecting cores explicitly with taskset will not work
            # if we try to oversubscribe. We would need to implement some kind of wrap-around
            # mechanism
            if make_jobs_per_test == 0 or ctest_jobs_per_test == 0:
                expect(False, "test-all-scream does not currently support oversubscription. "
                              "Either run fewer test types or turn off parallel testing")

            for test in self._tests:
                test.testing_res_count = ctest_jobs_per_test
                test.compile_res_count = make_jobs_per_test
                print(f"Test {test} can use {test.compile_res_count} jobs to compile, and {test.testing_res_count} jobs for test")

        # Unless the user claims to know what he/she is doing, we setup the env.
        # Need to happen before compiler probing
        if not self._preserve_env:
            # Setup the env on this machine
            setup_mach_env(self._machine.name, ctest_j=ctest_max_jobs)

        ############################################
        #           Check repo status              #
        ############################################

        expect(get_current_commit(), f"Root dir: {self._root_dir}, does not appear to be a git repo")

        ###################################
        #      Compute baseline info      #
        ###################################

        expect (not self._baseline_dir or self._work_dir != self._baseline_dir,
                f"Error! For your safety, do NOT use '{self._work_dir}' (the work_dir) to store baselines. Move them to a different directory (even a subdirectory if that works).")

        # These two dir are special dir for "on-the-fly baselines" and "machine's official baselines"
        local_baseline_dir = self._work_dir/"baselines"
        auto_dir = Path(self._machine.baselines_dir).expanduser().absolute()
        # Handle the "fake" auto case, used in scripts tests
        if "SCREAM_FAKE_AUTO" in os.environ:
            auto_dir = auto_dir / "fake"

        if self._baseline_dir == "LOCAL":
            self._baseline_dir = local_baseline_dir
        elif self._baseline_dir == "AUTO":
            self._baseline_dir = auto_dir
        elif self._baseline_dir is not None:
            self._baseline_dir = Path(self._baseline_dir).absolute()

        # Make the baseline dir, if not already existing.
        if self._generate:
            expect(self._baseline_dir is not None, "Cannot generate without --baseline-dir")
            self.create_tests_dirs(self._baseline_dir, clean=False)

        # Check baselines status
        if self._baseline_dir is not None:
            print (f"Checking baselines directory: {self._baseline_dir}")
            missing_baselines = self.check_baselines_are_present()
            expect (len(missing_baselines)==0 or self._generate,
                    f"Missing baselines for builds {missing_baselines}. Re-run with -g to generate them")
        else:
            print("baseline_dir is None! Skipping all baseline tests")

        ############################################
        #    Deduce compilers if needed/possible   #
        ############################################

        if self._cxx_compiler is None:
            self._cxx_compiler = self._machine.cxx_compiler
        if self._f90_compiler is None:
            self._f90_compiler = self._machine.ftn_compiler
        if self._c_compiler is None:
            self._c_compiler = self._machine.c_compiler

    ###############################################################################
    def create_tests_dirs(self, root, clean):
    ###############################################################################

        # Make sure the tests root directory exists
        root.mkdir(parents=True,exist_ok=True)

        # Create build directories (one per test)
        for test in self._tests:

            test_dir = self.get_test_dir(root, test)

            if test_dir.exists() and clean:
                # LB: without '._str', I was getting the error
                # TypeError: lstat: illegal type for path parameter
                shutil.rmtree(str(test_dir))

            # Create this built type's build dir (if not already existing)
            test_dir.mkdir(parents=True,exist_ok=True)

            # Create the 'data' subdir (if not already existing)
            (test_dir / "data").mkdir(parents=False,exist_ok=True)

    ###############################################################################
    def set_baseline_file_sha(self, test):
    ###############################################################################
        sha = get_current_commit()
        baseline_file = (self.get_preexisting_baseline(test).parent)/"baseline_git_sha"
        with baseline_file.open("w", encoding="utf-8") as fd:
            return fd.write(sha)

    ###############################################################################
    def get_test_dir(self, root, test):
    ###############################################################################
        return root/str(test)

    ###############################################################################
    def get_root_dir(self):
    ###############################################################################
        return self._root_dir

    ###############################################################################
    def get_machine(self):
    ###############################################################################
        return self._machine

    ###############################################################################
    def get_preexisting_baseline(self, test):
    ###############################################################################
        expect(self._baseline_dir is not None, "Cannot supply preexisting baseline without baseline dir")
        return self._baseline_dir/str(test)/"data"

    ###############################################################################
    def check_baselines_are_present(self):
    ###############################################################################
        """
        Check that all baselines are present (one subdir for all values of self._tests)
        Note: if test.uses_baselines=False, skip the check
        """
        # Sanity check
        expect(self._baseline_dir is not None,
                "Error! Baseline directory not correctly set.")

        missing = []
        for test in self._tests:
            if test.uses_baselines:
                data_dir = self.get_preexisting_baseline(test)
                if not data_dir.is_dir():
                    test.baselines_missing = True
                    missing.append(test.longname)
                    print(f" -> Test {test} is missing baselines")
                else:
                    print(f" -> Test {test} appears to have baselines")
            else:
                print(f" -> Test {test} does not use baselines")

        return missing

    ###############################################################################
    def generate_cmake_config(self, test, for_ctest=False):
    ###############################################################################

        # Ctest only needs config options, and doesn't need the leading 'cmake '
        result  = f"{'' if for_ctest else 'cmake '}-C {self._machine.mach_file}"

        # Netcdf should be available. But if the user is doing a testing session
        # where all netcdf-related code is disabled, he/she should be able to run
        # even if no netcdf is available
        stat, f_path, _ = run_cmd("nf-config --prefix")
        if stat == 0:
            result += f" -DNetCDF_Fortran_PATH={f_path}"
        stat, c_path, _ = run_cmd("nc-config --prefix")
        if stat == 0:
            result += f" -DNetCDF_C_PATH={c_path}"
        stat, pc_path, _ = run_cmd("pnetcdf-config --prefix")
        if stat == 0:
            result += f" -DPnetCDF_C_PATH={pc_path}"

        # Test-specific cmake options
        for key, value in test.cmake_args:
            result += f" -D{key}={value}"

        # The output coming from all tests at the same time will be a mixed-up mess
        # unless we tell test-launcher to buffer all output
        if self._extra_verbose:
            result += " -DEKAT_TEST_LAUNCHER_BUFFER=True"

        # Define test level
        if self._test_level == "at":
            pass # Nothing to do, this is the default in the cmake system
        elif self._test_level == "nightly":
            result += " -DSCREAM_TEST_LEVEL=NIGHTLY"
        elif self._test_level == "experimental":
            result += " -DSCREAM_TEST_LEVEL=EXPERIMENTAL"

        # Define test size. Default is to let it be undefined and CMake will pick
        # unless test has a special default test length
        if self._test_size:
            result += f" -DSCREAM_TEST_SIZE={self._test_size.upper()}"
        elif test.default_test_len:
            result += f" -DSCREAM_TEST_SIZE={test.default_test_len.upper()}"

        # User-requested config options
        custom_opts_keys = []
        for custom_opt in self._custom_cmake_opts:
            expect ("=" in custom_opt, "Error! Syntax error in custom cmake options. Should be `VAR_NAME=VALUE`.")
            if "=" in custom_opt:
                name, value = custom_opt.split("=", 1)
                # Some effort is needed to ensure quotes are perserved
                result += f" -D{name}='{value}'"
                custom_opts_keys.append(name)

        # Common config options (unless already specified by the user)
        if "CMAKE_CXX_COMPILER" not in custom_opts_keys:
            result += f" -DCMAKE_CXX_COMPILER={self._cxx_compiler}"
        if "CMAKE_C_COMPILER" not in custom_opts_keys:
            result += f" -DCMAKE_C_COMPILER={self._c_compiler}"
        if "CMAKE_Fortran_COMPILER" not in custom_opts_keys:
            result += f" -DCMAKE_Fortran_COMPILER={self._f90_compiler}"

        if "SCREAM_DYNAMICS_DYCORE" not in custom_opts_keys:
            result += " -DSCREAM_DYNAMICS_DYCORE=HOMME"

        # For the compute-sanitizer tool 'racecheck', if no option --racecheck-num-workers
        # is provided, it will attempt to use all threads available on node. This can cause
        # issues if other test cases are being run in parallel. If the option was not specified,
        # limit the number of threads availible to racecheck to the number of compile resources.
        if self._parallel and isinstance(test, CSR) and not '--racecheck-num-workers' in result:
            new_option = ' --racecheck-num-workers=' + str(test.compile_res_count)
            index = result.index('--tool=racecheck') + len('--tool=racecheck')
            result = result[:index] + new_option + result[index:]

        return result

    ###############################################################################
    def get_taskset_resources(self, test, for_compile=True):
    ###############################################################################
        res_name = "compile_res_count" if for_compile else "testing_res_count"

        if not for_compile and self._machine.uses_gpu():
            # For GPUs, the cpu affinity is irrelevant. Just assume all GPUS are open
            affinity_cp = list(range(self._ctest_max_jobs))
        elif "SLURM_CPU_BIND_LIST" in os.environ:
            affinity_cp = get_cpu_ids_from_slurm_env_var()
        else:
            this_process = psutil.Process()
            affinity_cp = list(this_process.cpu_affinity())

        affinity_cp.sort()

        if self._parallel:
            it = itertools.takewhile(lambda item: item != test, self._tests)
            offset = sum(getattr(prevs, res_name) for prevs in it)
        else:
            offset = 0

        expect(offset < len(affinity_cp),
               f"Offset {offset} out of bounds (max={len(affinity_cp)}) for test {test}\naffinity_cp: {affinity_cp}")
        resources = []
        for i in range(0, getattr(test, res_name)):
            resources.append(affinity_cp[offset+i])

        return resources

    ###############################################################################
    def create_ctest_resource_file(self, test, build_dir):
    ###############################################################################
        # Create a json file in the test build dir, which ctest will then use
        # to schedule tests in parallel.
        # In the resource file, we have N res groups with 1 slot, with N being
        # what's in test.testing_res_count. On CPU machines, res groups
        # are cores, on GPU machines, res groups are GPUs. In other words, a
        # res group is where we usually bind an individual MPI rank.
        # The id of the res groups on is offset so that it is unique across all builds

        resources = self.get_taskset_resources(test, for_compile=False)

        data = {}

        # This is the only version numbering supported by ctest, so far
        data["version"] = {"major":1,"minor":0}

        # We add leading zeroes to ensure that ids will sort correctly
        # both alphabetically and numerically
        devices = []
        for res_id in resources:
            devices.append({"id":f"{res_id:05d}"})

        # Add resource groups
        data["local"] = [{"devices":devices}]

        with (build_dir/"ctest_resource_file.json").open("w", encoding="utf-8") as outfile:
            json.dump(data,outfile,indent=2)

        return len(resources)

    ###############################################################################
    def generate_ctest_config(self, cmake_config, extra_configs, test):
    ###############################################################################
        result = ""

        test_dir = self.get_test_dir(self._work_dir,test)
        num_test_res = self.create_ctest_resource_file(test,test_dir)
        cmake_config += f" -DSCREAM_TEST_MAX_TOTAL_THREADS={num_test_res}"
        verbosity = "-V --output-on-failure" if not self._extra_verbose else "-VV"

        result += f"SCREAM_BUILD_PARALLEL_LEVEL={test.compile_res_count} CTEST_PARALLEL_LEVEL={test.testing_res_count} ctest {verbosity} "
        result += f"--resource-spec-file {test_dir}/ctest_resource_file.json "

        if self._baseline_dir is not None and test.uses_baselines:
            cmake_config += f" -DSCREAM_BASELINES_DIR={self.get_preexisting_baseline(test).parent}"
        else:
            cmake_config += " -DSCREAM_ENABLE_BASELINE_TESTS=Off"

        if not self._submit:
            result += "-DNO_SUBMIT=True "

        if isinstance(test, COV):
            result += "-DDO_COVERAGE=True "

        for key, value in extra_configs:
            result += f"-D{key}={value} "

        work_dir = self._work_dir/str(test)
        result += f"-DBUILD_WORK_DIR={work_dir} "

        build_name_mod = str(test)
        result += f"-DBUILD_NAME_MOD={build_name_mod} "

        if self._limit_test_regex:
            result += f"-DINCLUDE_REGEX={self._limit_test_regex} "
        result += f'-S {self._root_dir}/cmake/ctest_script.cmake -DCTEST_SITE={self._machine.name} -DCMAKE_COMMAND="{cmake_config}" '

        # Ctest can only competently manage test pinning across a single instance of ctest. For
        # multiple concurrent instances of ctest, we have to help it. It's OK to use the compile_res_count
        # taskset resources even though the ctest script is also running the tests
        if self._parallel:
            resources = self.get_taskset_resources(test)
            result = result.replace("'", r"'\''") # handle nested quoting
            result = f"taskset -c {','.join([str(r) for r in resources])} sh -c '{result}'"

        return result

    ###############################################################################
    def generate_baselines(self, test):
    ###############################################################################
        expect(test.uses_baselines,
               f"Something is off. generate_baseline should have not be called for test {test}")

        baseline_dir = self.get_test_dir(self._baseline_dir, test)
        test_dir = self.get_test_dir(self._work_dir, test)
        if test_dir.exists():
            shutil.rmtree(test_dir)
        test_dir.mkdir()

        num_test_res = self.create_ctest_resource_file(test,test_dir)
        cmake_config = self.generate_cmake_config(test)
        cmake_config +=  " -DSCREAM_ONLY_GENERATE_BASELINES=ON"
        cmake_config += f" -DSCREAM_BASELINES_DIR={baseline_dir}"
        cmake_config += f" -DSCREAM_TEST_MAX_TOTAL_THREADS={num_test_res}"

        print("===============================================================================")
        print(f"Generating baseline for test {test} with config '{cmake_config}'")
        print("===============================================================================")

        # We cannot just crash if we fail to generate baselines, since we would
        # not get a dashboard report if we did that. Instead, just ensure there is
        # no baseline file to compare against if there's a problem.
        stat, _, err = run_cmd(f"{cmake_config} {self._root_dir}",
                               from_dir=test_dir, verbose=True)
        if stat != 0:
            print (f"WARNING: Failed to create baselines (config phase):\n{err}")
            return False

        cmd = f"make -j{test.compile_res_count}"
        if self._parallel:
            resources = self.get_taskset_resources(test)
            cmd = f"taskset -c {','.join([str(r) for r in resources])} sh -c '{cmd}'"

        stat, _, err = run_cmd(cmd, from_dir=test_dir, verbose=True)

        if stat != 0:
            print (f"WARNING: Failed to create baselines (build phase):\n{err}")
            return False

        cmd  = f"ctest -j{test.testing_res_count}"
        cmd +=  " -L baseline_gen"
        cmd += f" --resource-spec-file {test_dir}/ctest_resource_file.json"
        stat, _, err = run_cmd(cmd, from_dir=test_dir, verbose=True)

        if stat != 0:
            print (f"WARNING: Failed to create baselines (run phase):\n{err}")
            return False

        # Read list of nc files to copy to baseline dir
        with open(test_dir/"data/baseline_list","r",encoding="utf-8") as fd:
            files = fd.read().splitlines()

            with SharedArea():
                for fn in files:
                    # In case appending to the file leaves an empty line at the end
                    if fn != "":
                        src = Path(fn)
                        dst = baseline_dir / "data" / src.name
                        shutil.copyfile(src, dst)

        # Store the sha used for baselines generation. This is only for record
        # keeping.
        self.set_baseline_file_sha(test)
        test.baselines_missing = False

        # Clean up the directory by removing everything
        shutil.rmtree(test_dir)

        return True

    ###############################################################################
    def generate_all_baselines(self):
    ###############################################################################
        git_head = get_current_head()

        tests_needing_baselines = self.baselines_to_be_generated()
        if len(tests_needing_baselines)==0:
            return True

        print("###############################################################################")
        print(f"Generating baselines using '{git_head}'")
        print("###############################################################################")

        success = True
        num_workers = len(tests_needing_baselines) if self._parallel else 1
        with threading3.ProcessPoolExecutor(max_workers=num_workers) as executor:

            future_to_test = {
                executor.submit(self.generate_baselines, test) : test
                for test in tests_needing_baselines}

            for future in threading3.as_completed(future_to_test):
                test = future_to_test[future]
                success &= future.result()

        return success

    ###############################################################################
    def run_test(self, test):
    ###############################################################################
        git_head = get_current_head()

        print("===============================================================================")
        print(f"Testing '{git_head}' for test '{test}'")
        print("===============================================================================")

        test_dir = self.get_test_dir(self._work_dir,test)
        cmake_config = self.generate_cmake_config(test, for_ctest=True)
        ctest_config = self.generate_ctest_config(cmake_config, [], test)

        if self._config_only:
            ctest_config += "-DCONFIG_ONLY=TRUE"

        if self._quick_rerun and (test_dir/"CMakeCache.txt").is_file():
            # Do not purge bld dir, and do not rerun config step.
            # Note: make will still rerun cmake if some cmake file has changed
            if self._quick_rerun_failed:
                ctest_config += "--rerun-failed "
        else:
            # This directory might have been used before during another test-all-scream run.
            # Although it's ok to build in the same dir, we MUST make sure to erase cmake's cache
            # and internal files from the previous build (CMakeCache.txt and CMakeFiles folder),
            # Otherwise, we may not pick up changes in certain cmake vars that are already cached.
            run_cmd_no_fail("rm -rf CMake*", from_dir=test_dir)

        success = run_cmd(ctest_config, from_dir=test_dir, arg_stdout=None, arg_stderr=None, verbose=True)[0] == 0

        return success

    ###############################################################################
    def run_all_tests(self):
    ###############################################################################
        print("###############################################################################")
        print("Running tests!")
        print("###############################################################################")

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

        for t,s in tests_success.items():
            if not s:
                last_test = self.get_last_ctest_file(t,"TestsFailed")
                last_build  = self.get_last_ctest_file(t,"Build")
                last_config = self.get_last_ctest_file(t,"Configure")
                if last_test is not None:
                    print(f"Build type {t} failed at testing time. Here's a list of failed tests:")
                    print (last_test.read_text())
                elif last_build is not None:
                    print(f"Build type {t} failed at build time. Here's the build log:")
                    print (last_build.read_text())
                elif last_config is not None:
                    print(f"Build type {t} failed at config time. Here's the config log:\n\n")
                    print (last_config.read_text())
                else:
                    print(f"Build type {t} failed before configure step.")

        return success

    ###############################################################################
    def get_last_ctest_file(self,test,phase):
    ###############################################################################
        test_dir = self.get_test_dir(self._work_dir,test)
        test_results_dir = Path(test_dir,"Testing","Temporary")
        files = list(test_results_dir.glob(f"Last{phase}*"))
        if len(files)>0:
            curr_tag=0
            curr_idx=0
            latest = None
            # ctest creates files of the form Last{phase}_$TIMESTAMP-$IDX.log
            # we split the name, and sue $TIMESTAMP to pick the newest, and, in case
            # of tie, $IDX as tiebreaker
            for file in files:
                file_no_path = file.name
                tokens = re.split(r"_|-|\.",str(file_no_path))
                if latest is None:
                    latest = file
                    curr_tag = int(tokens[1])
                    curr_idx = int(tokens[2])
                else:
                    if int(tokens[1])>curr_tag:
                        latest = file
                        curr_tag = int(tokens[1])
                        curr_idx = int(tokens[2])
                    elif int(tokens[1])==curr_tag and int(tokens[2])>curr_idx:
                        latest = file
                        curr_tag = int(tokens[1])
                        curr_idx = int(tokens[2])

            return latest
        else:
            return None

    ###############################################################################
    def baselines_to_be_generated(self):
    ###############################################################################
        """
        Return list of baselines to generate. Baselines will always be generated
        for tests that use baselines.
        """
        ret = []
        for test in self._tests:
            if test.uses_baselines:
                ret.append(test)

        return ret

    ###############################################################################
    def test_all_scream(self):
    ###############################################################################

        # Add any override the user may have requested
        for env_var in self._custom_env_vars:
            key,val = env_var.split("=",2)
            os.environ.update( { key : val } )

        success = True

        if self._generate:
            success = self.generate_all_baselines()

        else:
            # First, create build directories (one per test). If existing, nuke the content
            self.create_tests_dirs(self._work_dir, not self._quick_rerun)

            success &= self.run_all_tests()
            if not success:
                print ("Error(s) occurred during test phase")

        return success
