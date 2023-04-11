from utils import run_cmd, run_cmd_no_fail, expect, check_minimum_python_version, ensure_psutil
from git_utils import get_current_head, get_current_commit, get_current_branch, is_repo_clean, \
    cleanup_repo, merge_git_ref, checkout_git_ref, git_refs_difference, print_last_commit

from machines_specs import get_mach_compilation_resources, get_mach_testing_resources, \
    get_mach_baseline_root_dir, setup_mach_env, is_cuda_machine, \
    get_mach_cxx_compiler, get_mach_f90_compiler, get_mach_c_compiler, is_machine_supported, \
    logical_cores_per_physical_core

check_minimum_python_version(3, 4)

import os, shutil
import concurrent.futures as threading3
import itertools
import json

ensure_psutil()
import psutil
import re

from collections import OrderedDict
from pathlib import Path

###############################################################################
class TestProperty(object):
###############################################################################

    def __init__(self, longname, description, cmake_args, uses_baselines=True, on_by_default=True):
        # What the user uses to select tests via test-all-scream CLI.
        # Should also match the class name when converted to caps
        self.shortname      = type(self).__name__.lower()

        # A longer name used to name baseline and test directories for a test.
        # Also used in output/error messages to refer to the test
        self.longname       = longname

        # A longer decription of the test
        self.description    = description

        # Cmake config args for this test
        self.cmake_args     = cmake_args

        # Does the test do baseline testing
        self.uses_baselines = uses_baselines

        # Should this test be run if the user did not specify tests at all?
        self.on_by_default  = on_by_default

        #
        # Properties not set by constructor (Set by the main TestAllScream object)
        #

        # Resources used by this test.
        self.compile_res_count = None
        self.testing_res_count = None

        # Does this test need baselines
        self.missing_baselines = False

        #
        # Common
        #

        if not self.uses_baselines:
            self.cmake_args += [("SCREAM_ENABLE_BASELINE_TESTS", "False")]

    # Tests will generally be referred to via their longname
    def disable_baselines(self):
        if self.uses_baselines:
            self.uses_baselines = False
            self.cmake_args += [("SCREAM_ENABLE_BASELINE_TESTS", "False")]

    # Tests will generally be referred to via their longname
    def __str__(self):
        return self.longname

###############################################################################
class DBG(TestProperty):
###############################################################################

    CMAKE_ARGS = [("CMAKE_BUILD_TYPE", "Debug"), ("EKAT_DEFAULT_BFB", "True")]

    def __init__(self):
        TestProperty.__init__(
            self,
            "full_debug",
            "debug",
            self.CMAKE_ARGS,
        )

###############################################################################
class SP(TestProperty):
###############################################################################

    def __init__(self):
        TestProperty.__init__(
            self,
            "full_sp_debug",
            "debug single precision",
            DBG.CMAKE_ARGS + [("SCREAM_DOUBLE_PRECISION", "False")],
        )

###############################################################################
class FPE(TestProperty):
###############################################################################

    def __init__(self):
        TestProperty.__init__(
            self,
            "debug_nopack_fpe",
            "debug pksize=1 floating point exceptions on",
            DBG.CMAKE_ARGS + [("SCREAM_PACK_SIZE", "1"), ("SCREAM_FPE","True")],
            uses_baselines=False
        )

###############################################################################
class OPT(TestProperty):
###############################################################################

    def __init__(self):
        TestProperty.__init__(
            self,
            "release",
            "release",
            [("CMAKE_BUILD_TYPE", "Release")],
        )

###############################################################################
class COV(TestProperty):
###############################################################################

    def __init__(self):
        TestProperty.__init__(
            self,
            "coverage",
            "debug coverage",
            [("CMAKE_BUILD_TYPE", "Debug"), ("EKAT_ENABLE_COVERAGE", "True")],
            uses_baselines=False,
            on_by_default=False
        )

###############################################################################
def test_factory(user_req_tests, machine, mem_check):
###############################################################################
    testclasses = TestProperty.__subclasses__()
    if not user_req_tests:
        result = [testclass() for testclass in testclasses
            if (testclass().on_by_default and not (is_cuda_machine(machine) and testclass().shortname == "fpe"))]
    else:
        valid_names = [testclass().shortname for testclass in testclasses]
        for user_req_test in user_req_tests:
            expect(user_req_test in valid_names, f"'{user_req_test}' is not a known test")

        result = [testclass() for testclass in testclasses if testclass().shortname in user_req_tests]

    # Mangle test full name if mem_check is on
    if mem_check:
        for test in result:
            test.longname += "_cuda_mem_check" if is_cuda_machine(machine) else "_valgrind"

    return result

###############################################################################
class TestAllScream(object):
###############################################################################

    ###########################################################################
    @classmethod
    def get_test_name_desc(cls):
    ###########################################################################
        """
        Returns a dict mapping short test names to full names
        """
        testclasses = TestProperty.__subclasses__()
        return OrderedDict([(testc().shortname, testc().description) for testc in testclasses])

    ###########################################################################
    def __init__(self, cxx_compiler=None, f90_compiler=None, c_compiler=None,
                 submit=False, parallel=False, fast_fail=False,
                 baseline_ref=None, baseline_dir=None, machine=None, no_tests=False, config_only=False, keep_tree=False,
                 custom_cmake_opts=(), custom_env_vars=(), preserve_env=False, tests=(),
                 integration_test=False, local=False, root_dir=None, work_dir=None,
                 quick_rerun=False,quick_rerun_failed=False,dry_run=False,
                 make_parallel_level=0, ctest_parallel_level=0, update_expired_baselines=False,
                 extra_verbose=False, limit_test_regex=None, test_level="at",
                 mem_check=False, force_baseline_regen=False):
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
        self._fast_fail               = fast_fail
        self._baseline_ref            = baseline_ref
        self._machine                 = machine
        self._local                   = local
        self._perform_tests           = not no_tests
        self._config_only             = config_only
        self._keep_tree               = keep_tree
        self._baseline_dir            = baseline_dir
        self._custom_cmake_opts       = custom_cmake_opts
        self._custom_env_vars         = custom_env_vars
        self._preserve_env            = preserve_env
        self._root_dir                = root_dir
        self._work_dir                = None if work_dir is None else Path(work_dir)
        self._integration_test        = integration_test
        self._quick_rerun             = quick_rerun
        self._quick_rerun_failed      = quick_rerun_failed
        self._dry_run                 = dry_run
        self._update_expired_baselines= update_expired_baselines
        self._extra_verbose           = extra_verbose
        self._limit_test_regex        = limit_test_regex
        self._test_level              = test_level
        self._mem_check               = mem_check
        self._force_baseline_regen    = force_baseline_regen

        # Not all builds are ment to perform comparisons against pre-built baselines

        if self._quick_rerun_failed:
            self._quick_rerun = True

        ############################################
        #  Sanity checks and helper structs setup  #
        ############################################

        # Quick rerun skips config phase, and config-only runs only config. You can't ask for both...
        expect (not (self._quick_rerun and self._config_only),
                "Makes no sense to ask for --quick-rerun and --config-only at the same time")

        # Probe machine if none was specified
        if self._machine is None:
            # We could potentially integrate more with CIME here to do actual
            # nodename probing.
            if "SCREAM_MACHINE" in os.environ and is_machine_supported(os.environ["SCREAM_MACHINE"]):
                self._machine = os.environ["SCREAM_MACHINE"]
            else:
                expect(self._local,
                       "test-all-scream requires either the machine arg (-m $machine) or the -l flag,"
                       "which makes it look for machine specs in '~/.cime/scream_mach_specs.py'.")
                self._machine = "local"
        else:
            expect (not self._local, "Specifying a machine while passing '-l,--local' is ambiguous.")

        # Make our test objects!
        self._tests = test_factory(tests, self._machine, self._mem_check)

        # Compute root dir (where repo is) and work dir (where build/test will happen)
        if not self._root_dir:
            self._root_dir = Path(__file__).resolve().parent.parent
        else:
            self._root_dir = Path(self._root_dir).resolve()
            expect(self._root_dir.is_dir() and self._root_dir.parts()[-2:] == ('scream', 'components'),
                   f"Bad root-dir '{self._root_dir}', should be: $scream_repo/components/eamxx")

        if self._work_dir is not None:
            self._work_dir = Path(self._work_dir).absolute()
            expect(self._work_dir.is_dir(),
                   f"Error! Work directory '{self._work_dir}' does not exist.")
        else:
            self._work_dir = self._root_dir.absolute().joinpath("ctest-build")
            self._work_dir.mkdir(exist_ok=True)

        os.chdir(str(self._root_dir)) # needed, or else every git command will need repo=root_dir
        expect(get_current_commit(), f"Root dir: {self._root_dir}, does not appear to be a git repo")

        # Print some info on the branch
        self._original_branch = get_current_branch()
        self._original_commit = get_current_commit()

        print_last_commit(git_ref=self._original_branch, dry_run=self._dry_run)

        ###################################
        #  Compilation/testing resources  #
        ###################################

        # Deduce how many compilation resources per test
        make_max_jobs = get_mach_compilation_resources()
        if make_parallel_level > 0:
            expect(make_parallel_level <= make_max_jobs,
                   f"Requested make_parallel_level {make_parallel_level} is more than max available {make_max_jobs}")
            make_max_jobs = make_parallel_level
            print(f"Note: honoring requested value for make parallel level: {make_max_jobs}")
        else:
            print(f"Note: no value passed for --make-parallel-level. Using the default for this machine: {make_max_jobs}")

        ctest_max_jobs = get_mach_testing_resources(self._machine)
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
            if is_cuda_machine(self._machine):
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
            setup_mach_env(self._machine, ctest_j=ctest_max_jobs)

        ###################################
        #      Compute baseline info      #
        ###################################

        expect (not self._baseline_dir or self._work_dir != self._baseline_dir,
                f"Error! For your safety, do NOT use '{self._work_dir}' to store baselines. Move them to a different directory (even a subdirectory if that works).")

        # If no baseline ref/dir was provided, use default master baseline dir for this machine
        # NOTE: if user specifies baseline ref, baseline dir will be set later to a path within work dir
        if self._baseline_dir is None and self._baseline_ref is None:
            self._baseline_dir = "AUTO"
            print ("No '--baseline-dir XYZ' nor '-b XYZ' provided. Testing against default baselines dir for this machine.")

        # If -k was used, make sure it's allowed
        if self._keep_tree:
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
                # Make sure the baseline ref is unset (or HEAD)
                expect(self._baseline_ref is None or self._baseline_ref == "HEAD",
                       "The option --keep-tree is only available when testing against pre-built baselines "
                       "(--baseline-dir) or HEAD (-b HEAD)")
        else:
            expect(self._dry_run or is_repo_clean(),
                   "Repo must be clean before running. If testing against HEAD or pre-built baselines, "
                   "you can pass `--keep-tree` to allow non-clean repo.")

        # For integration test, enforce baseline_ref==origin/master, and proceed to merge origin/master
        if self._integration_test:
            expect (self._baseline_ref is None or self._baseline_ref=="origin/master",
                    "Error! Integration tests cannot be done against an arbitrary baseline ref.")

            # Set baseline ref and merge it
            self._baseline_ref = "origin/master"
            merge_git_ref(git_ref=self._baseline_ref, verbose=True, dry_run=self._dry_run)

            # Always update expired baselines if this is an integration test
            self._update_expired_baselines = True

        # By now, we should have at least one between baseline_dir and baseline_ref set (possibly both)
        default_baselines_root_dir = self._work_dir/"baselines"
        if self._baseline_dir is None:
            # Use default baseline dir, and create it if necessary
            self._baseline_dir = Path(default_baselines_root_dir).absolute()
            self.create_tests_dirs(self._baseline_dir, True) # Wipe out previous baselines

        else:
            if self._baseline_dir == "AUTO":
                expect (self._baseline_ref is None or self._baseline_ref == 'origin/master',
                        "Do not specify `-b XYZ` when using `--baseline-dir AUTO`. The AUTO baseline dir should be used for the master baselines only.\n"
                        "       `-b XYZ` needs to probably build baselines for ref XYZ. However, no baselines will be built if the dir already contains baselines.\n")
                # We treat the "AUTO" string as a request for automatic baseline dir.
                auto_dir = get_mach_baseline_root_dir(self._machine)
                self._baseline_dir = Path(auto_dir) if auto_dir else default_baselines_root_dir
                if "SCREAM_FAKE_AUTO" in os.environ:
                    self._baseline_dir = self._baseline_dir/"fake"
            else:
                self._baseline_dir = Path(self._baseline_dir).absolute()

            # Make sure the baseline folders exist (but do not purge content if they exist)
            self.create_tests_dirs(self._baseline_dir, False)

        # Do not do baseline operations if mem checking is on
        print (f"Checking baselines directory: {self._baseline_dir}")
        if self._mem_check:
            for test in self._tests:
                test.disable_baselines()
        else:
            self.baselines_are_present()
            if self._update_expired_baselines:
                self.baselines_are_expired()

        ############################################
        #    Deduce compilers if needed/possible   #
        ############################################

        if self._cxx_compiler is None:
            self._cxx_compiler = get_mach_cxx_compiler(self._machine)
        if self._f90_compiler is None:
            self._f90_compiler = get_mach_f90_compiler(self._machine)
        if self._c_compiler is None:
            self._c_compiler = get_mach_c_compiler(self._machine)

        if not self._dry_run:
            self._f90_compiler = run_cmd_no_fail(f"which {self._f90_compiler}")
            self._cxx_compiler = run_cmd_no_fail(f"which {self._cxx_compiler}")
            self._c_compiler   = run_cmd_no_fail(f"which {self._c_compiler}")

    ###############################################################################
    def create_tests_dirs(self, root, clean):
    ###############################################################################

        # Make sure the baseline root directory exists
        root.mkdir(parents=True,exist_ok=True)

        # Create build directories (one per test)
        for test in self._tests:

            test_dir = self.get_test_dir(root, test)

            if test_dir.exists() and clean:
                # LB: without '._str', I was getting the error
                # TypeError: lstat: illegal type for path parameter
                shutil.rmtree(str(test_dir))

            # Create this baseline's build dir
            if not test_dir.exists():
                test_dir.mkdir(parents=True)

    ###############################################################################
    def get_baseline_file_sha(self, test):
    ###############################################################################
        baseline_file = (self.get_preexisting_baseline(test).parent)/"baseline_git_sha"
        if baseline_file.exists():
            with baseline_file.open("r", encoding="utf-8") as fd:
                return fd.read().strip()

        return None

    ###############################################################################
    def set_baseline_file_sha(self, test, sha):
    ###############################################################################
        baseline_file = (self.get_preexisting_baseline(test).parent)/"baseline_git_sha"
        with baseline_file.open("w", encoding="utf-8") as fd:
            return fd.write(sha)

    ###############################################################################
    def get_test_dir(self, root, test):
    ###############################################################################
        return root/str(test)

    ###############################################################################
    def get_preexisting_baseline(self, test):
    ###############################################################################
        expect(self._baseline_dir is not None, "Cannot supply preexisting baseline without baseline dir")
        return self._baseline_dir/str(test)/"data"

    ###############################################################################
    def baselines_are_present(self):
    ###############################################################################
        """
        Check that all baselines are present (one subdir for all values of self._tests)
        Note: if test.uses_baselines=False, skip the check
        """
        # Sanity check
        expect(self._baseline_dir is not None,
                "Error! Baseline directory not correctly set.")

        for test in self._tests:
            if test.uses_baselines:
                data_dir = self.get_preexisting_baseline(test)
                if not data_dir.is_dir():
                    test.missing_baselines = True
                    print(f" -> Test {test} is missing baselines")
                else:
                    print(f" -> Test {test} appears to have baselines")
            else:
                print(f" -> Test {test} does not use baselines")

    ###############################################################################
    def baselines_are_expired(self):
    ###############################################################################
        """
        Baselines are expired if either:
          1) there is no file in baseline_dir containing the sha of the baselines
          2) the baselines sha does not match baseline_ref
        """
        baseline_ref_sha = get_current_commit(commit=self._baseline_ref)

        # Sanity check
        expect(self._baseline_dir is not None, "Error! This routine should only be called when testing against pre-existing baselines.")

        for test in self._tests:
            if test.uses_baselines and not test.missing_baselines:
                # this test is not missing a baseline, but it may be expired.

                baseline_file_sha = self.get_baseline_file_sha(test)
                if baseline_file_sha is None:
                    test.missing_baselines = True
                    print(f" -> Test {test} has no stored sha so must be considered expired")
                else:
                    num_ref_is_behind_file, num_ref_is_ahead_file = git_refs_difference(baseline_file_sha, baseline_ref_sha)

                    # If the copy in our repo is behind, then we need to update the repo
                    expect (num_ref_is_behind_file==0 or not self._integration_test,
f"""Error! Your repo seems stale, since the baseline sha in your repo is behind
the one last used to generated them. We do *not* allow an integration
test to replace baselines with older ones, for security reasons.
If this is a legitimate case where baselines need to be 'rewound',
e.g. b/c of a (hopefully VERY RARE) force push to master, then
remove existing baselines first. Otherwise, please run 'git fetch $remote'.
 - baseline_ref: {self._baseline_ref}
 - repo baseline sha: {baseline_ref_sha}
 - last used baseline sha: {baseline_file_sha}""")

                    # If the copy in our repo is not ahead, then baselines are not expired
                    if num_ref_is_ahead_file > 0 or self._force_baseline_regen:
                        test.missing_baselines = True
                        reason = "forcing baseline regen" if self._force_baseline_regen \
                                 else "f{self._baseline_ref} is ahead of the baseline commit by {num_ref_is_ahead_file}"
                        print(f" -> Test {test} baselines are expired because {reason}")
                    else:
                        print(f" -> Test {test} baselines are valid and do not need to be regenerated")

    ###############################################################################
    def get_machine_file(self):
    ###############################################################################
        if self._local:
            return Path("~/.cime/scream_mach_file.cmake").expanduser()
        else:
            return self._root_dir/"cmake"/"machine-files"/f"{self._machine}.cmake"

    ###############################################################################
    def generate_cmake_config(self, extra_configs, for_ctest=False):
    ###############################################################################

        # Ctest only needs config options, and doesn't need the leading 'cmake '
        result  = f"{'' if for_ctest else 'cmake '}-C {self.get_machine_file()}"

        # Netcdf should be available. But if the user is doing a testing session
        # where all netcdf-related code is disabled, he/she should be able to run
        # even if no netcdf is available
        stat, f_path, _ = run_cmd("nf-config --prefix")
        if stat == 0:
            result += f" -DNetCDF_Fortran_PATH={f_path}"
        stat, c_path, _ = run_cmd("nc-config --prefix")
        if stat == 0:
            result += f" -DNetCDF_C_PATH={c_path}"

        # Test-specific cmake options
        for key, value in extra_configs:
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

        # Add configs for mem checking
        if self._mem_check:
            # valgrind/cmc slow down things a lot, we need to run tests in short mode
            result += " -DSCREAM_TEST_SIZE=SHORT"
            result += f" -DEKAT_ENABLE_{'CUDA_MEMCHECK' if is_cuda_machine(self._machine) else 'VALGRIND'}=True"

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

        return result

    ###############################################################################
    def get_taskset_range(self, test, for_compile=True):
    ###############################################################################
        res_name = "compile_res_count" if for_compile else "testing_res_count"

        if not for_compile and is_cuda_machine(self._machine):
            # For GPUs, the cpu affinity is irrelevant. Just assume all GPUS are open
            affinity_cp = list(range(self._ctest_max_jobs))
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
        start = affinity_cp[offset]
        end = start
        for i in range(1, getattr(test, res_name)):
            expect(affinity_cp[offset+i] == start+i, f"Could not get contiguous range for test {test}")
            end = affinity_cp[offset+i]

        return start, end

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

        start, end = self.get_taskset_range(test, for_compile=False)

        data = {}

        # This is the only version numbering supported by ctest, so far
        data['version'] = {"major":1,"minor":0}

        # We add leading zeroes to ensure that ids will sort correctly
        # both alphabetically and numerically
        devices = []
        for res_id in range(start,end+1):
            devices.append({"id":f"{res_id:05d}"})

        # Add resource groups
        data['local'] = [{"devices":devices}]

        with (build_dir/"ctest_resource_file.json").open('w', encoding="utf-8") as outfile:
            json.dump(data,outfile,indent=2)

        return (end-start)+1

    ###############################################################################
    def generate_ctest_config(self, cmake_config, extra_configs, test):
    ###############################################################################
        result = ""
        if self._submit:
            result += f"SCREAM_MACHINE={self._machine} "

        test_dir = self.get_test_dir(self._work_dir,test)
        num_test_res = self.create_ctest_resource_file(test,test_dir)
        cmake_config += f" -DSCREAM_TEST_MAX_TOTAL_THREADS={num_test_res}"
        verbosity = "-V --output-on-failure" if not self._extra_verbose else "-VV"

        result += f"SCREAM_BUILD_PARALLEL_LEVEL={test.compile_res_count} CTEST_PARALLEL_LEVEL={test.testing_res_count} ctest {verbosity} "
        result += f"--resource-spec-file {test_dir}/ctest_resource_file.json "

        if self._baseline_dir is not None and test.uses_baselines:
            cmake_config += f" -DSCREAM_TEST_DATA_DIR={self.get_preexisting_baseline(test)}"

        if not self._submit:
            result += "-DNO_SUBMIT=True "

        if isinstance(test, COV):
            result += "-DDO_COVERAGE=True "

        for key, value in extra_configs:
            result += f"-D{key}={value} "

        work_dir = self._work_dir/str(test)
        result += f"-DBUILD_WORK_DIR={work_dir} "

        build_name_mod = str(test)
        if self._mem_check:
            if self._submit:
                # Need backwards compatible build name for db submission
                expect(len(self._tests) == 1, "Expected only one mem-check test being submitted")
                build_name_mod = "cuda_mem_check" if is_cuda_machine(self._machine) else "valgrind"

        result += f"-DBUILD_NAME_MOD={build_name_mod} "

        if self._limit_test_regex:
            result += f"-DINCLUDE_REGEX={self._limit_test_regex} "
        result += f'-S {self._root_dir}/cmake/ctest_script.cmake -DCMAKE_COMMAND="{cmake_config}" '

        # Ctest can only competently manage test pinning across a single instance of ctest. For
        # multiple concurrent instances of ctest, we have to help it. It's OK to use the compile_res_count
        # taskset range even though the ctest script is also running the tests
        if self._parallel:
            start, end = self.get_taskset_range(test)
            result = f"taskset -c {start}-{end} sh -c '{result}'"

        return result

    ###############################################################################
    def generate_baselines(self, test, commit):
    ###############################################################################
        expect(test.uses_baselines,
               f"Something is off. generate_baseline should have not be called for test {test}")

        test_dir = self.get_test_dir(self._baseline_dir, test)

        cmake_config = self.generate_cmake_config(test.cmake_args)
        cmake_config += " -DSCREAM_BASELINES_ONLY=ON"
        cmake_config += f" -DSCREAM_TEST_DATA_DIR={test_dir}/data"

        print("===============================================================================")
        print(f"Generating baseline for test {test} with config '{cmake_config}'")
        print("===============================================================================")

        success = True

        try:
            # We cannot just crash if we fail to generate baselines, since we would
            # not get a dashboard report if we did that. Instead, just ensure there is
            # no baseline file to compare against if there's a problem.
            stat, _, err = run_cmd(f"{cmake_config} {self._root_dir}",
                                   from_dir=test_dir, verbose=True, dry_run=self._dry_run)
            if stat != 0:
                print (f"WARNING: Failed to configure baselines:\n{err}")
                success = False

            else:
                cmd = f"make -j{test.compile_res_count} && make -j{test.testing_res_count} baseline"
                if self._parallel:
                    start, end = self.get_taskset_range(test)
                    cmd = f"taskset -c {start}-{end} sh -c '{cmd}'"

                stat, _, err = run_cmd(cmd, from_dir=test_dir, verbose=True, dry_run=self._dry_run)

                if stat != 0:
                    print(f"WARNING: Failed to create baselines:\n{err}")
                    success = False

        finally:
            # Clean up the directory, by removing everything but the 'data' subfolder. This must
            # happen unconditionally or else subsequent runs could be corrupted
            run_cmd_no_fail(r"find -maxdepth 1 -not -name data ! -path . -exec rm -rf {} \;",
                            from_dir=test_dir, verbose=True, dry_run=self._dry_run)

        if success:
            # Store the sha used for baselines generation
            self.set_baseline_file_sha(test, commit)
            test.missing_baselines = False

        return success

    ###############################################################################
    def generate_all_baselines(self):
    ###############################################################################
        git_head_ref = get_current_head()

        print("###############################################################################")
        print(f"Generating baselines for ref {self._baseline_ref}")
        print("###############################################################################")

        commit = get_current_commit(commit=self._baseline_ref)

        # Switch to the baseline commit
        checkout_git_ref(self._baseline_ref, verbose=True, dry_run=self._dry_run)

        success = True
        tests_needing_baselines = [test for test in self._tests if test.missing_baselines]
        num_workers = len(tests_needing_baselines) if self._parallel else 1
        with threading3.ProcessPoolExecutor(max_workers=num_workers) as executor:

            future_to_test = {
                executor.submit(self.generate_baselines, test, commit) : test
                for test in tests_needing_baselines}

            for future in threading3.as_completed(future_to_test):
                test = future_to_test[future]
                success &= future.result()

                if not success and self._fast_fail:
                    print(f"Generation of baselines for test {test} failed")
                    return False

        # Switch back to the branch commit
        checkout_git_ref(git_head_ref, verbose=True, dry_run=self._dry_run)

        return success

    ###############################################################################
    def run_test(self, test):
    ###############################################################################
        git_head = get_current_head()

        print("===============================================================================")
        print(f"Testing '{git_head}' for test '{test}'")
        print("===============================================================================")

        test_dir = self.get_test_dir(self._work_dir,test)
        cmake_config = self.generate_cmake_config(test.cmake_args, for_ctest=True)
        ctest_config = self.generate_ctest_config(cmake_config, [], test)

        if self._config_only:
            ctest_config += "-DCONFIG_ONLY=TRUE"

        if self._quick_rerun and (test_dir/"CMakeCache.txt").is_file():
            # Do not purge bld dir, and do not rerun config step.
            # Note: make will still rerun cmake if some cmake file has changed
            if self._quick_rerun_failed:
                ctest_config += "--rerun-failed "
        else:
            # This directory might have been used also to build the model to generate baselines.
            # Although it's ok to build in the same dir, we MUST make sure to erase cmake's cache
            # and internal files from the previous build (CMakeCache.txt and CMakeFiles folder)
            run_cmd_no_fail("rm -rf CMake*", from_dir=test_dir, dry_run=self._dry_run)

        success = run_cmd(ctest_config, from_dir=test_dir, arg_stdout=None, arg_stderr=None, verbose=True, dry_run=self._dry_run)[0] == 0

        return success

    ###############################################################################
    def run_all_tests(self):
    ###############################################################################
        print("###############################################################################")
        print("Running tests!")
        print("###############################################################################")

        # First, create build directories (one per test). If existing, nuke the content
        self.create_tests_dirs(self._work_dir, not self._quick_rerun)

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
                last_test = self.get_last_ctest_file(t,"Tests")
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
                tokens = re.split(r'_|-|\.',str(file_no_path))
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
    def test_all_scream(self):
    ###############################################################################

        # Add any override the user may have requested
        for env_var in self._custom_env_vars:
            key,val = env_var.split("=",2)
            os.environ.update( { key : val } )

        success = True
        try:
            # If needed, generate baselines first
            tests_needing_baselines = [test for test in self._tests if test.missing_baselines]
            if tests_needing_baselines:
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
                cleanup_repo(self._original_branch, self._original_commit, dry_run=self._dry_run)

        return success
