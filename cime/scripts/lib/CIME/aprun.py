"""
Aprun is far too complex to handle purely through XML. We need python
code to compute and assemble aprun commands.
"""

from CIME.XML.standard_module_setup import *

import math

logger = logging.getLogger(__name__)

###############################################################################
def _get_aprun_cmd_for_case_impl(ntasks, nthreads, rootpes, pstrids,
                                 max_tasks_per_node, max_mpitasks_per_node,
                                 pio_numtasks, pio_async_interface,
                                 compiler, machine, run_exe):
###############################################################################
    """
    No one really understands this code, but we can at least test it.

    >>> ntasks = [512, 675, 168, 512, 128, 168, 168, 512, 1]
    >>> nthreads = [2, 2, 2, 2, 4, 2, 2, 2, 1]
    >>> rootpes = [0, 0, 512, 0, 680, 512, 512, 0, 0]
    >>> pstrids = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    >>> max_tasks_per_node = 16
    >>> max_mpitasks_per_node = 16
    >>> pio_numtasks = -1
    >>> pio_async_interface = False
    >>> compiler = "pgi"
    >>> machine = "titan"
    >>> run_exe = "e3sm.exe"
    >>> _get_aprun_cmd_for_case_impl(ntasks, nthreads, rootpes, pstrids, max_tasks_per_node, max_mpitasks_per_node, pio_numtasks, pio_async_interface, compiler, machine, run_exe)
    (' -S 4 -n 680 -N 8 -d 2 e3sm.exe : -S 2 -n 128 -N 4 -d 4 e3sm.exe ', 117, 808, 4, 4)
    >>> compiler = "intel"
    >>> _get_aprun_cmd_for_case_impl(ntasks, nthreads, rootpes, pstrids, max_tasks_per_node, max_mpitasks_per_node, pio_numtasks, pio_async_interface, compiler, machine, run_exe)
    (' -S 4 -cc numa_node -n 680 -N 8 -d 2 e3sm.exe : -S 2 -cc numa_node -n 128 -N 4 -d 4 e3sm.exe ', 117, 808, 4, 4)

    >>> ntasks = [64, 64, 64, 64, 64, 64, 64, 64, 1]
    >>> nthreads = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    >>> rootpes = [0, 0, 0, 0, 0, 0, 0, 0, 0]
    >>> pstrids = [1, 1, 1, 1, 1, 1, 1, 1, 1]
    >>> _get_aprun_cmd_for_case_impl(ntasks, nthreads, rootpes, pstrids, max_tasks_per_node, max_mpitasks_per_node, pio_numtasks, pio_async_interface, compiler, machine, run_exe)
    (' -S 8 -cc numa_node -n 64 -N 16 -d 1 e3sm.exe ', 4, 64, 16, 1)
    """
    max_tasks_per_node = 1 if max_tasks_per_node < 1 else max_tasks_per_node

    total_tasks = 0
    for ntask, rootpe, pstrid in zip(ntasks, rootpes, pstrids):
        tt = rootpe + (ntask - 1) * pstrid + 1
        total_tasks = max(tt, total_tasks)

    # Check if we need to add pio's tasks to the total task count
    if pio_async_interface:
        total_tasks += pio_numtasks if pio_numtasks > 0 else max_mpitasks_per_node

    # Compute max threads for each mpi task
    maxt = [0] * total_tasks
    for ntask, nthrd, rootpe, pstrid in zip(ntasks, nthreads, rootpes, pstrids):
        c2 = 0
        while c2 < ntask:
            s = rootpe + c2 * pstrid
            if nthrd > maxt[s]:
                maxt[s] = nthrd

            c2 += 1

    # make sure all maxt values at least 1
    for c1 in range(0, total_tasks):
        if maxt[c1] < 1:
            maxt[c1] = 1

    # Compute task and thread settings for batch commands
    tasks_per_node, min_tasks_per_node, task_count, thread_count, max_thread_count, total_node_count, total_task_count, aprun_args = \
        0, max_mpitasks_per_node, 1, maxt[0], maxt[0], 0, 0, ""
    c1list = list(range(1, total_tasks))
    c1list.append(None)
    for c1 in c1list:
        if c1 is None or maxt[c1] != thread_count:
            tasks_per_node = min(max_mpitasks_per_node, int(max_tasks_per_node / thread_count))

            tasks_per_node = min(task_count, tasks_per_node)

            # Compute for every subset
            task_per_numa = int(math.ceil(tasks_per_node / 2.0))
            # Option for Titan
            if machine == "titan" and tasks_per_node > 1:
                aprun_args += " -S {:d}".format(task_per_numa)
                if compiler == "intel":
                    aprun_args += " -cc numa_node"

            aprun_args += " -n {:d} -N {:d} -d {:d} {} {}".format(task_count, tasks_per_node, thread_count, run_exe, "" if c1 is None else ":")

            node_count = int(math.ceil(float(task_count) / tasks_per_node))
            total_node_count += node_count
            total_task_count += task_count

            if tasks_per_node < min_tasks_per_node:
                min_tasks_per_node = tasks_per_node

            if c1 is not None:
                thread_count = maxt[c1]
                max_thread_count = max(max_thread_count, maxt[c1])
                task_count = 1

        else:
            task_count += 1

    return aprun_args, total_node_count, total_task_count, min_tasks_per_node, max_thread_count

###############################################################################
def get_aprun_cmd_for_case(case, run_exe, overrides=None):
###############################################################################
    """
    Given a case, construct and return the aprun command and optimized node count
    """
    models = case.get_values("COMP_CLASSES")
    ntasks, nthreads, rootpes, pstrids = [], [], [], []
    for model in models:
        model = "CPL" if model == "DRV" else model
        for the_list, item_name in zip([ntasks, nthreads, rootpes, pstrids],
                                       ["NTASKS", "NTHRDS", "ROOTPE", "PSTRID"]):
            the_list.append(case.get_value("_".join([item_name, model])))
    max_tasks_per_node = case.get_value("MAX_TASKS_PER_NODE")
    if overrides:
        if 'max_tasks_per_node' in overrides:
            max_tasks_per_node = overrides['max_tasks_per_node'] 
        if 'total_tasks' in overrides:
            ntasks = [overrides['total_tasks'] if x > 1 else x for x in ntasks]
        if 'thread_count' in overrides:
            nthreads = [overrides['thread_count'] if x > 1 else x for x in nthreads]
    


    return _get_aprun_cmd_for_case_impl(ntasks, nthreads, rootpes, pstrids,
                                        max_tasks_per_node,
                                        case.get_value("MAX_MPITASKS_PER_NODE"),
                                        case.get_value("PIO_NUMTASKS"),
                                        case.get_value("PIO_ASYNC_INTERFACE"),
                                        case.get_value("COMPILER"),
                                        case.get_value("MACH"),
                                        run_exe)
