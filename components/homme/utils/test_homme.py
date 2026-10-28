#!/usr/bin/env python3
"""
Simple script to run ctest (in dashboard mode) to test Homme
"""

import argparse
import pathlib
import subprocess
import os
import sys
import shutil
import shlex
import getpass

DEFAULT_ROOT = pathlib.Path(__file__).resolve().parent.parent
E3SM_ROOT = DEFAULT_ROOT.parent.parent

# pylint: disable=wrong-import-position, import-error, no-name-in-module
sys.path.append( str(E3SM_ROOT / "cime") )
import CIME.XML.machines
# pylint: enable=wrong-import-position, import-error, no-name-in-module

def _bool_str(value):
    return "TRUE" if value else "FALSE"


def _normalize_cmake_opt(opt):
    return opt if opt.startswith("-D") else f"-D{opt}"


def _parse_args(argv):

    parser = argparse.ArgumentParser(
        description="Drive standalone HOMME testing through ctest -S."
    )

    parser.add_argument(
        "-m",
        "--machine",
        required=True,
        help="Name of the machine. MUST be a valid CIME machine. Also,"
             "components/homme/cmake/machineFiles/{machine}.cmake MUST exist.",
    )
    parser.add_argument(
        "--compiler",
        default=None,
        help="In case the machine supports 2+ compilers, the user can ask for a specific one",
    )
    parser.add_argument(
        "--machine-file",
        default=None,
        help="Name (without -bfb suffix and .cmake extension) of the file in components/homme/cmake/machineFiles "
             "to use, instead of the one named like the machine. Useful if the machine supports 2+ compilers, "
             "and each needs its own file (e.g., ghci-gnu).",
    )
    parser.add_argument(
        "-r",
        "--root-dir",
        default=DEFAULT_ROOT,
        help="Path to components/homme. Defaults to script-relative path.",
    )
    parser.add_argument(
        "-w",
        "--work-dir",
        help="Build/test work directory. Defaults to <root-dir>/ctest-build.",
    )

    parser.add_argument(
        "-g",
        "--generate",
        action="store_true",
        help="If true, generate baselines, otherwise compare against existing ones"
    )
    parser.add_argument(
        "--no-baseline",
        action="store_true",
        help="Do not compare against baselines (e.g., for --bfb, where we only care about f90 vs cxx). "
             "Cannot be used with --generate.",
    )
    parser.add_argument(
        "--baseline-root",
        default=None,
        help="Baseline root directory.",
    )
    parser.add_argument(
        "-b",
        "--baseline-name",
        default=getpass.getuser(),
        help="Baseline case name under baseline-root (e.g., master).",
    )

    parser.add_argument(
        "--build-parallel-level",
        type=int,
        default=-1,
        help="Parallel level for building test-execs.",
    )
    parser.add_argument(
        "--run-parallel-level",
        type=int,
        default=-1,
        help="Parallel level for baseline/check target invocations.",
    )

    parser.add_argument(
        "--bfb",
        action="store_true",
        help="Use the '<machine>-bfb.cmake' machine file and enable cxx-f90 BFB testing.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Configure with CMAKE_BUILD_TYPE=Debug.",
    )
    parser.add_argument(
        "--cprnc-dir",
        help="Directory containing cprnc executable (passed as CPRNC_DIR).",
    )
    parser.add_argument(
        "--cmake-opt",
        action="append",
        default=[],
        help="Extra cmake option, with or without -D prefix. May repeat.",
    )

    parser.add_argument(
        "-s",
        "--submit",
        action="store_true",
        help="Submit results to CDash after testing.",
    )
    parser.add_argument(
        "--ctest-arg",
        action="append",
        default=[],
        help="Extra argument forwarded to ctest (e.g., --ctest-args=-V).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print command and exit without executing.",
    )

    return parser.parse_args(argv[1:])


def expect (cond, msg):
    """
    Wrap condition check. Helps lower pylint complain about "too many branches"
    """
    if not cond:
        raise SystemExit(msg)

def ensure_paths(args,machine):
    """
    Ensures the root, work, and baseline dirs exist 
    """
    expect(not (args.no_baseline and args.generate),
           "ERROR: --no-baseline and --generate cannot be used together")

    root_dir = pathlib.Path(args.root_dir).resolve()

    expect(root_dir.is_dir() and root_dir.name == "homme",
            f"ERROR: root-dir must point to /some/path/components/homme. Got: {root_dir}")

    if args.work_dir:
        work_dir = pathlib.Path(args.work_dir).resolve()
    else:
        work_dir = root_dir / "ctest-build"

    if args.no_baseline:
        return root_dir, work_dir, None

    # We may need some mach defaults here
    if args.compiler:
        compiler = args.compiler
    else:
        compiler = machine.get_value('COMPILER')

    if args.baseline_root:
        baseline_dir = pathlib.Path(args.baseline_root).resolve()
    else:
        machine.set_value('COMPILER', compiler)
        raw_name = machine.get_value('BASELINE_ROOT', resolved=True)
        baseline_dir = pathlib.Path(raw_name).resolve() / args.baseline_name / "homme"

    if args.generate:
        # Ok if dir does not exist, but we MUST be able to create it
        try:
            baseline_dir.mkdir(parents=True, exist_ok=True)
        except (PermissionError, FileNotFoundError):
            print(f"ERROR: Cannot create '{baseline_dir}'")
            raise

        # Double check write access just to be completely sure
        expect(os.access(baseline_dir, os.W_OK),
               f"ERROR: Cannot create write to '{baseline_dir}'")
    else:
        expect(baseline_dir.is_dir(),
               f"ERROR: baseline dir does not exist: '{baseline_dir}'")
        expect(os.access(baseline_dir, os.R_OK),
               f"ERROR: baseline dir '{baseline_dir}' is not readable.")

    return root_dir, work_dir, baseline_dir

def _build_ctest_cmd(args,root_dir,work_dir,baseline_dir,machine):
    cmake_opts = [_normalize_cmake_opt(opt) for opt in args.cmake_opt]

    ctest_script = root_dir / "cmake" / "ctest_script.cmake"

    make_j = args.build_parallel_level
    if args.build_parallel_level <= 0:
        make_j = machine.get_value('GMAKE_J')
    ctest_j = args.run_parallel_level
    if args.run_parallel_level <= 0:
        ctest_j = machine.get_value('MAX_TASKS_PER_NODE')
    # The machine file is named after the machine, unless the user asks for a different one
    machine_file = args.machine_file or args.machine
    if args.bfb:
        machine_file += "-bfb"

    cmd  =  "ctest -V --output-on-failure"
    cmd +=  " " + " ".join(args.ctest_arg)
    cmd += f" -S {ctest_script}"
    cmd += f" -D HOMME_ROOT={root_dir}"
    cmd += f" -D BUILD_WORK_DIR={work_dir}"
    cmd += f" -D MACHINE_FILE={machine_file}"
    cmd += f" -D BFB_TESTING={_bool_str(args.bfb)}"
    if args.no_baseline:
        cmd += " -D NO_BASELINE=TRUE"
    else:
        cmd += f" -D BASELINE_DIR={baseline_dir}"
    cmd += f" -D DEBUG={_bool_str(args.debug)}"
    cmd += f" -D GENERATE={args.generate}"
    cmd += f" -D BUILD_PARALLEL_LEVEL={make_j}"
    cmd += f" -D RUN_PARALLEL_LEVEL={ctest_j}"
    cmd += f" -D SUBMIT={_bool_str(args.submit)}"

    if args.cprnc_dir:
        cmd += f" -D CPRNC_DIR={pathlib.Path(args.cprnc_dir).resolve()}"
    if cmake_opts:
        # The opts are joined with ';' (a cmake list), so quote it to keep the shell from splitting it
        cmd += " -D " + shlex.quote(f"EXTRA_CMAKE_OPTIONS={';'.join(cmake_opts)}")

    return cmd

def clean_work_dir(work_dir,root_dir):
    """
    Remove the build left in work_dir by a previous run (if any), so that stale state cannot affect this run.
    E.g., an old CMakeCache.txt makes cmake silently ignore a changed --cmake-opt or machine file, and old test
    outputs in work_dir/tests/baseline would end up in the baselines when generating them.

    We only remove a dir if it is a cmake build dir of this homme source tree. If the dir is not empty and we
    cannot tell that, we abort rather than risk removing something else (e.g., if the user passes a source dir).
    """
    if not work_dir.exists():
        return

    expect(work_dir.is_dir() and not work_dir.is_symlink(),
           f"ERROR: work dir '{work_dir}' is not a directory")

    if not any(work_dir.iterdir()):
        return

    # CMake records the source dir that a build dir was configured from
    src_dir = None
    cache = work_dir / "CMakeCache.txt"
    if cache.is_file():
        for line in cache.read_text(errors="replace").splitlines():
            if line.startswith("CMAKE_HOME_DIRECTORY:"):
                src_dir = pathlib.Path(line.split("=",1)[1]).resolve()
                break

    expect(src_dir == root_dir,
           f"ERROR: work dir '{work_dir}' is not empty, and it is not a cmake build dir of '{root_dir}'.\n"
           f"       Refusing to remove it. Remove it manually, or choose a different --work-dir.")

    print(f"Removing previous build dir '{work_dir}'")
    shutil.rmtree(work_dir)

def main(argv):
    """
    Main function
    """
    # 0. Parse args and validate paths
    args = _parse_args(argv)

    config_machines = str(E3SM_ROOT / "cime_config/machines/config_machines.xml")
    machine = CIME.XML.machines.Machines(config_machines,machine=args.machine)

    root_dir, work_dir, baseline_dir = ensure_paths(args,machine)

    cmd = _build_ctest_cmd(args,root_dir,work_dir,baseline_dir,machine)

    print(f"RUN: {cmd}")
    if args.dry_run:
        return 0

    # Always start from a clean build dir, so no stale state (e.g., an old CMakeCache.txt) is silently reused
    clean_work_dir(work_dir,root_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    # 1. Run the command and capture the process result object
    res = subprocess.run(cmd, cwd=str(work_dir), check=False, shell=True)

    # 2. Check if the CTest process failed
    if res.returncode != 0:
        print(f"ERROR: ctest failed (exit code {res.returncode}). The logs are in '{work_dir / 'Testing' / 'Temporary'}'")
        return res.returncode

    # 3. Process succeeded! If needed, copy to baseline dir
    if args.generate:
        print(f"Baseline generation successful. Copying baselines/artifacts to {baseline_dir}...")
        try:
            tgt = baseline_dir
            src = work_dir
            shutil.copytree(src / 'tests/baseline', tgt / 'tests/baseline', dirs_exist_ok=True)
            shutil.copytree(src / 'Testing', tgt / 'Testing', dirs_exist_ok=True)
            shutil.copy2(src / 'CMakeCache.txt', tgt / 'CMakeCache.txt')
            shutil.copy2(src / 'homme_git_sha.h', tgt / 'homme_git_sha.h')
            print("Successfully copied baseline and artifacts.")
        except (OSError, ValueError) as e:
            print(f"ERROR: Failed to copy to baseline dir: {e}")
            return 1

    # 4. Return the successful 0 code
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
