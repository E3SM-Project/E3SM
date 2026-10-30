"""
A script to wrap all e3sm compilations. This script should be enabled
by prefixing it to all compiler and link calls. This can be done easily in
CMake by using the RULE_LAUNCH_COMPILE and RULE_LAUNCH_LINK global properties.

This script will allow us to do whatever bookeeping, timing, logging, etc that
we want in our build system.

We want this script to be super-lean, so we do not load any of the standard CIME
stuff.
"""

import sys, subprocess, time, psutil, os

###############################################################################
def run_cmd(args):
###############################################################################
    t1 = time.time()

    # Start the subprocess
    proc = subprocess.Popen(args)

    # Monitor memory usage throughout the process
    peak_memory = 0  # in bytes
    try:
        process = psutil.Process(proc.pid)
        while proc.poll() is None:
            try:
                mem_info = process.memory_info()
                peak_memory = max(peak_memory, mem_info.rss)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
            time.sleep(0.1)

        # Get final memory info
        try:
            mem_info = process.memory_info()
            peak_memory = max(peak_memory, mem_info.rss)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        pass

    result = proc.returncode
    t2 = time.time()

    arglen = len(args)
    target = None
    for idx, arg in enumerate(args):
        if arg == "-o" and idx + 1 < arglen:
            target = args[idx + 1]
            break

        if arg.startswith("lib") and arg.endswith(".a"):
            target = arg

    if target is None:
        target = " ".join(args)

    # Convert peak_memory to MB for readability
    peak_memory_mb = peak_memory / (1024 * 1024)

    # Printing the time this way ensures no mixing of output with other
    # concurrent processes
    sys.stderr.write("Target {} built in {:f} seconds, peak memory: {:.1f} MB\n".format(
        target, (t2 - t1), peak_memory_mb))

    return result

###############################################################################
def parse_command_line(args, _):
###############################################################################
    return args[1:]

###############################################################################
def _main_func(description):
###############################################################################
    cmd_args = parse_command_line(sys.argv, description)

    result = run_cmd(cmd_args)

    sys.exit(result)

###############################################################################

if __name__ == "__main__":
    _main_func(__doc__)
