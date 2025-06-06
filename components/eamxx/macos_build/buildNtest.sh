#! /bin/zsh

cleanCache() {
  rm -f $1/cmakeCache.txt
}

build_dir="ctest-build/full_debug"

cleanCache "${build_dir}"

# ==============================================================================
outfile='build-test.log'
outfile_old='build-test.old.log'
# ==============================================================================
makeNtest_ncores="6"
make_ncores="${makeNtest_ncores}"
test_ncores="4"
pbuild="-p --make-parallel-level ${make_ncores}"
ptest="--ctest-parallel-level ${test_ncores}"
# ==============================================================================
local='-l'

test_config='dbg'
t_config="-t ${test_config}"

vflag="--extra-verbose"
notest="--no-tests"
quick_flag="--quick-rerun"
misc_flag="--quick-rerun-failed"

bl_loc="-b ${HOME}/scream-data"
bl_gen="-g"
# ==============================================================================
# test_regex="[^\s]*wetscav[^\s]*"
# tchoice="--limit-test-regex ${test_regex}"
# ==============================================================================
base_flags="${local} ${pbuild} ${ptest}"
# base_flags="${local}"
# ==============================================================================
bl_test_flags="${base_flags} ${bl_loc}"
bl_gen_flags="${base_flags} ${bl_loc} ${bl_gen}"
# ==============================================================================
fullDbg_flags="${base_flags} ${t_config}"
verbose_fullDbg_flags="${fullDbg_flags} ${vflag}"
testRegex_verbose_fullDbg_flags="${verbose_fullDbg_flags} ${tchoice}"
quickRR_testRegex_verbose_fullDbg_flags="${testRegex_verbose_fullDbg_flags} ${quick_flag}"
# ==============================================================================
test_flags=${verbose_fullDbg_flags}
# ==============================================================================
# ==============================================================================
echo "TEST FLAGS = ${test_flags}"
eval "./scripts/test-all-eamxx ${test_flags}" |& tee "${outfile}"
# testing_status=${PIPESTATUS[0]}
# cp ${outfile} ${outfile_old}
# if (( testing_status > 0 )); then
#   echo "****ERRORS IN MAKE--STATUS = ${testing_status}****"
#   exit ${testing_status}
# fi
# ==============================================================================

# usage:
# test-all-scream <ARGS> [--verbose]
# OR
# test-all-scream --help

# EXAMPLES (assumes user is on machine mappy):
#     # Run all tests on current machine using the SCREAM-approved env for this machine
#     > cd $scream_repo/components/eamxx
#     > ./scripts/test-all-scream -m mappy

#     # Run all tests on current machine with defaut behavior except using your current shell env
#     > cd $scream_repo/components/eamxx
#     > ./scripts/test-all-scream --preserve-env -m mappy

#     # Run all tests on current machine with default behavior except using non-default baselines
#     > cd $scream_repo/components/eamxx
#     > ./scripts/test-all-scream -m mappy --baseline-dir=PATH_TO_BASELINES

#     # Run all tests on current machine with default behavior except using local baselines
#     > cd $scream_repo/components/eamxx
#     > ./scripts/test-all-scream -m mappy --baseline-dir=LOCAL

#     # Run only the dbg test on current machine with default behavior otherwise
#     > cd $scream_repo/components/eamxx
#     > ./scripts/test-all-scream -m mappy -t dbg

# Drive ctest testing of scream for a complete set of tests. This will be our
# gold standard to determine if the code is working or not on the current platform.
# For batch machines, this script expects to already be on a compute node; it does not
# do batch submissions.

# IMPORTANT: the default behavior of this script *changes your environment*,
#            by loading machine-specific modules and setting machine-specific
#            env vars. To prevent this behavior, use --preserve-env flag.

# Baselines: By default, test-all-scream will not run baseline tests. If you set
# -b AUTO, baseline tests will be done with the pre-existing public baselines in
# the location specified by the machine spec. You can regenerate
# baselines any time by using the -g flag, but be aware this will impact everyone
# if you regenerate the public baselines. You can change the target baseline area
# using -b $path. You can also use the magic word "LOCAL" to have test-all-scream
# pick a local directory for you if you want to manage your own baselines within
# the current repo. If -g is provided, no tests will be run; -g means generate only.

# The general workflow for baseline-changing PRs is:
# 1) Issue PR
# 2) AT will fail with differences in the baseline tests
# 3) Review and merge.
# 4) Ask JimF or LucaB for a bless to baselines on mappy and weaver.

# If you are developing a branch and baselines tests are failing unexpectedly, it is
# likely that your branch has fallen out of date. You should upstream merge or rebase
# your branch.

# To sum up, the baseline handling for test-all-scream should basically match what we
# do for create_test tests, the only difference is that test-all-scream does baseline
# comparison tests by default.

# optional arguments:
#   -h, --help            show this help message and exit
#   -cxx CXX_COMPILER, --cxx-compiler CXX_COMPILER
#                         C++ compiler (default: None)
#   -f90 F90_COMPILER, --f90-compiler F90_COMPILER
#                         F90 compiler (default: None)
#   --c-compiler C_COMPILER
#                         C compiler (default: None)
#   -s, --submit          Submit results to dashboad (default: False)
#   -p, --parallel        Launch the different build types stacks in parallel
#                         (default: False)
#   -g, --generate        Instruct test-all-scream to generate baselines from
#                         current commit. Skips tests (default: False)
#   -b BASELINE_DIR, --baseline-dir BASELINE_DIR
#                         Directory where baselines should be read from (or
#                         written to, if -g is used). Default is None which
#                         skips all baseline tests. AUTO means use public
#                         baselines. You can also use LOCAL to manage baselines
#                         in your local work dir. Lastly, you can provide a path
#                         here as well. (default: None)
#   -m MACHINE, --machine MACHINE
#                         Provide machine name. This is *always* required. It
#                         can, but does nothave to, match SCREAM_MACHINE. You
#                         can decorate this with compilerinfo if a machine
#                         supports multiple compiler types. This value will
#                         beused as the CTEST_SITE for cdash if the tests are
#                         submitted. It isexpected that a scream machine file
#                         exists for this value. (default: None)
#   --config-only         In the testing phase, only run config step, skip build
#                         and tests (default: False)
#   -c CUSTOM_CMAKE_OPTS, --custom-cmake-opts CUSTOM_CMAKE_OPTS
#                         Extra custom options to pass to cmake. Can use
#                         multiple times for multiple cmake options. The -D is
#                         added for you (default: [])
#   -e CUSTOM_ENV_VARS, --custom-env-vars CUSTOM_ENV_VARS
#                         Extra custom environment variables to be used. These
#                         will override(if applicable) whatever was found in
#                         machine_specs. Each -e flagsupports a single env var,
#                         so to pass multiple env var, do -e 'KEY1=VALUE1' -e
#                         'KEY2=VALUE2' (default: [])
#   --preserve-env        Whether to skip machine env setup, and preserve the
#                         current user env (useful to manually test new modules)
#                         (default: False)
#   -t TESTS, --test TESTS
#                         Only run specific test configurations, choices='dbg'
#                         (debug), 'sp' (debug single precision), 'fpe' (debug
#                         pksize=1 floating point exceptions on), 'opt'
#                         (release), 'cov' (debug coverage), 'valg' (debug with
#                         valgrind), 'csm' (debug with compute sanitizer
#                         memcheck), 'csr' (debug with compute sanitizer
#                         racecheck), 'csi' (debug with compute sanitizer
#                         initcheck), 'css' (debug with compute sanitizer
#                         synccheck) (default: [])
#   -l, --local           Allow to not specify a machine name, and have test-
#                         all-scream to lookfor '~/.cime/scream_mach_specs.py'
#                         for machine specifications. (default: False)
#   -r ROOT_DIR, --root-dir ROOT_DIR
#                         The root directory of the scream src you want to test.
#                         Default will be the scream src containing this script.
#                         (default: None)
#   -w WORK_DIR, --work-dir WORK_DIR
#                         The work directory where all the building/testing will
#                         happen. Defaults to ${root_dir}/ctest-build (default:
#                         None)
#   --quick-rerun         Do not clean the build dir, and do not reconfigure.
#                         Just (incremental) build and test. (default: False)
#   --quick-rerun-failed  Do not clean the build dir, and do not reconfigure.
#                         Just (incremental) build and retest failed tests only.
#                         (default: False)
#   --make-parallel-level MAKE_PARALLEL_LEVEL
#                         Max number of jobs to be created during compilation.
#                         If not provided, use default for given machine.
#                         (default: 0)
#   --ctest-parallel-level CTEST_PARALLEL_LEVEL
#                         Force to use this value for CTEST_PARALLEL_LEVEL. If
#                         not provided, use default for given machine. (default:
#                         0)
#   -x, --extra-verbose   Have ctest run in extra-verbose mode, which should
#                         print full output from all tests (default: False)
#   --limit-test-regex LIMIT_TEST_REGEX
#                         Limit ctest to running only tests that match this
#                         regex (default: None)
#   --test-level {at,nightly,experimental}
#                         Set the testing level. (default: at)
#   --test-size {short,medium,long}
#                         Set the testing level. Defaults to medium unless the
#                         test is cov or mem_check(short is default in those
#                         cases). (default: None)
