"""
Aprun is far too complex to handle purely through XML. We need python
code to compute and assemble aprun commands.
"""

from CIME.XML.standard_module_setup import *

import math

logger = logging.getLogger(__name__)

###############################################################################
def _get_aprun_cmd_for_case_impl(ntasks, nthreads, rootpes, pstrids,
                                 max_tasks_per_node, pes_per_node,
                                 pio_numtasks, pio_async_interface,
                                 compiler, machine, run_exe):
###############################################################################
    """
    No one really understands this code, but we can at least test it.

    >>> ntasks = []
    >>> nthreads = []
    >>> rootpes = []
    >>> pstrids = []
    >>> max_tasks_per_node = 1
    >>> pes_per_node = 1
    >>> pio_num_tasks = 1
    >>> pio_async_interface = True
    >>> compiler = "pgi"
    >>> machine = "titan"
    >>> run_exe = "acme.exe"
    >>> _get_aprun_cmd_for_case_impl(ntasks, nthreads, rootpes, pstrids, max_tasks_per_node, pes_per_node, pio_numtasks, pio_async_interface, compiler, machine, run_exe)
    ''
    """
    max_tasks_per_node = 1 if max_tasks_per_node < 1 else max_tasks_per_node

    total_tasks = 0
    for ntask, rootpe, pstrid in zip(ntasks, rootpes, pstrids):
        tt = rootpe + (ntask - 1) * pstrid + 1
        total_tasks = max(tt, total_tasks)

    # Check if we need to add pio's tasks to the total task count
    if pio_async_interface:
        total_tasks += pio_numtasks if pio_numtasks > 0 else pes_per_node

    # Compute max threads for each mpi task
    maxt = [0] * total_tasks
    for ntask, nthrd, rootpe, pstrid in zip(ntasks, nthreads, rootpes, pstrids):
        c2 = 0
        while c2 < ntask:
            s = rootpe + c2 * pstrid
            if nthrd > maxt[s]:
                maxt[s] = nthrd

            c2 += 1

    logger.info("total tasks is: %s" % total_tasks)

    # make sure all maxt values at least 1, don't know why we start at index 1
    for c1 in xrange(1, total_tasks):
        if maxt[c1] < 1:
            maxt[c1] = 1

    # Compute task and thread settings for batch commands
    tasks_per_node, task_count, thread_count, max_thread_count, aprun = \
        0, 1, maxt[0], maxt[0], ""
    for c1 in xrange(1, total_tasks):
        if maxt[c1] != thread_count:
            tasks_per_node = min(pes_per_node, max_tasks_per_node / thread_count)

            tasks_per_node = min(task_count, tasks_per_node)

            # Compute for every subset
            task_per_numa = int(math.ceil(tasks_per_node / 2.0))
            # Option for Titan
            if machine == "titan" and tasks_per_node > 1:
                if compiler == "intel":
                    aprun += " -S %d -cc numa_node " % task_per_numa
                else:
                    aprun += " -S %d " % task_per_numa

            aprun += " -n %d -N %d -d %d %s :" % (task_count, tasks_per_node, thread_count, run_exe)

            thread_count = maxt[c1]
            max_thread_count = max(max_thread_count, maxt[c1])
            task_count = 1

        else:
            task_count += 1

    if pes_per_node > 0:
        tasks_per_node = min(pes_per_node, max_tasks_per_node / thread_count)
    else:
        tasks_per_node = max_tasks_per_node / thread_count

    tasks_per_node = min(task_count, tasks_per_node)

    task_per_numa = int(math.ceil(tasks_per_node / 2.0))

    # Special option for Titan or intel compiler
    if compiler == "intel" or machine == "titan" and tasks_per_node > 1:
        aprun += " -S %d -cc numa_node " % task_per_numa

    aprun += " -n %d -N %d -d %d %s " % (task_count, tasks_per_node, thread_count, run_exe)

    return aprun

###############################################################################
def get_aprun_cmd_for_case(case, run_exe):
###############################################################################
    """
    Given a case, construct and return the aprun command
    """
    models = case.get_values("COMP_CLASSES")
    ntasks, nthreads, rootpes, pstrids = [], [], [], []
    for model in models:
        model = "CPL" if model == "DRV" else model
        for the_list, item_name in zip([ntasks, nthreads, rootpes, pstrids],
                                       ["NTASKS", "NTHREADS", "ROOTPE", "PSTRID"]):
            the_list.append(case.get_value("_".join([item_name, model])))

    return _get_aprun_cmd_for_case_impl(ntasks, nthreads, rootpes, pstrids,
                                        case.get_value("MAX_TASKS_PER_NODE"),
                                        case.get_value("PES_PER_NODE"),
                                        case.get_value("PIO_NUMTASKS"),
                                        case.get_value("PIO_ASYNC_INTERFACE"),
                                        case.get_value("COMPILER"),
                                        case.get_value("MACH"),
                                        run_exe)
