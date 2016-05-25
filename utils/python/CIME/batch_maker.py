"""
API for interacting with batch system

Provide a class hierarchy to ease the job of making batch scripts for
each CIME-ported machine.  We have class hierarchy and a factory class
to facilitate getting the right BatchMaker class for the appropriate
machine.

BatchMaker: This is the base class which contains functionality
common to making batch scripts for every machine. This class should
never be instantiated directly, use get_batch_maker for this.

We can have subclasses based on the batch system type, then subclasses
of those based on the machine.
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, run_cmd, get_model, convert_to_seconds, get_project, transform_vars
from CIME.case import Case
from CIME.XML.env_batch import EnvBatch
from CIME.XML.batch import Batch
from CIME.XML.machines import Machines
from CIME.task_maker import TaskMaker
from CIME.XML.lt_archive import LTArchive

import re, stat

logger = logging.getLogger(__name__)

class BatchMaker(object):

    def __init__(self, job, case=None):
        """
        Class constructor.  We need to know where in the filesystem we are,
        so caseroot, case, machroot, machine, cimeroot
        """
        self.case = case if case is not None else Case()
        self.job = job

        self.override_node_count = None

        # set up paths to the template files, this could and should be
        # extracted out somehow??
        self.job_id = self.case.get_value("CASE")

        if "pleiades" in self.case.get_value("MACH"):
            # pleiades jobname needs to be limited to 15 chars
            self.job_id = self.job_id[:15]

        self.output_error_path = self.case.get_value("CASE")

        self.env_batch = EnvBatch()

        self.machine = Machines(machine=self.case.get_value("MACH"))
        self.batch_system = self.machine.get_batch_system_type()
        self.batch_parser = Batch(self.batch_system, machine=self.machine.get_machine_name())

        self._initialize()

    def _initialize(self):
        self._set_task_info()
        self._set_queue()
        self._set_wall_time()
        self._set_project()
        self._set_lt_archive_options()
        self._set_model_run()
        self._set_batch_directives()
        self.env_batch.write()

    def _set_task_info(self):
        """
        Uses TaskMaker.pm to get the appropriate pe layout values for the run.
        Set the value as instance variables in the object.
        This can also be called from overrideNodeCount, in which case values can be
        manually overridden.
        """
        self.task_maker = TaskMaker(case=self.case)
        self.sumpes = self.task_maker.full_sum
        self.tasks_per_node = self.task_maker.task_per_node
        self.max_tasks_per_node = self.task_maker.MAX_TASKS_PER_NODE
        self.tasks_per_numa = self.task_maker.task_per_numa
        self.fullsum = self.task_maker.full_sum
        self.task_count = self.task_maker.full_sum
        self.sumtasks = self.task_maker.total_tasks
        self.num_tasks = self.task_maker.total_tasks
        self.totaltasks = self.task_maker.total_tasks
        self.maxthreads = self.task_maker.max_threads
        self.taskgeometry = self.task_maker.task_geom
        self.threadgeometry = self.task_maker.thread_geom
        self.taskcount = self.task_maker.task_count
        self.num_nodes = self.task_maker.node_count
        self.thread_count = self.task_maker.thread_count
        self.pedocumentation = self.task_maker.document()
        self.ptile = self.task_maker.ptile

        if self.override_node_count is not None:
            self.sumpes = self.override_node_count
            self.totaltasks = self.override_node_count
            self.fullsum = self.override_node_count
            self.sumtasks = self.override_node_count
            self.task_count = self.override_node_count
            self.num_nodes = self.override_node_count
            self.pedocumentation = ""

    def _set_queue(self):
        self.queue = self.env_batch.get_value("JOB_QUEUE", subgroup=self.job)
        self.ccsm_estcost = self.case.get_value("CCSM_ESTCOST")
        self.wall_time_max = None
        if self.queue:
            return

        # Make sure to check default queue first.
        all_queues = []
        all_queues.append( self.machine.get_default_queue())
        all_queues = all_queues + self.machine.get_all_queues()
        for queue in all_queues:
            if queue is not None:
                jobmin = queue.get("jobmin")
                jobmax = queue.get("jobmax")
                # if the fullsum is between the min and max # jobs, then use this queue.
                if jobmin is not None and jobmax is not None and self.fullsum >= int(jobmin) and self.fullsum <= int(jobmax):
                    self.queue = queue.text
                    self.wall_time_max = queue.get("walltimemax")
                    break

        if self.queue:
            self.case.set_value("JOB_QUEUE", self.queue, subgroup=self.job)
            logger.info("Using queue %s for job %s" % (self.queue,self.job))

    def _set_wall_time(self):
        # Get the wallclock time from env_batch.xml if its defined there
        # otherwise get the default from config_machines.xml
        # and set it in env_batch.xml
        self.wall_time = self.env_batch.get_value("JOB_WALLCLOCK_TIME", subgroup=self.job)

        # go through the walltime elements, and if our estimated cost is greater than the element's estimated cost,
        # then set the walltime.
        if not self.wall_time:
            for wall_time in self.machine.get_walltimes():
                estcost = wall_time.get("ccsm_estcost")
                if estcost is not None and self.ccsm_estcost > int(estcost):
                    self.wall_time = wall_time.text
        else:
            return

        # if we didn't find a walltime previously, use the default.
        if not self.wall_time:
            wall_time_node = self.machine.get_default_walltime()
            if wall_time_node is not None:
                self.wall_time = wall_time_node.text
            else:
                self.wall_time = "0"

        if self.wall_time_max and \
                convert_to_seconds(self.wall_time_max) < convert_to_seconds(self.wall_time):
            self.wall_time = self.wall_time_max

        self.env_batch.set_value("JOB_WALLCLOCK_TIME", self.wall_time, subgroup=self.job)

    def _set_project(self):
        if self.env_batch.get_value("PROJECT_REQUIRED", subgroup=self.job):
            self.project = self.env_batch.get_value("PROJECT", subgroup=self.job)
            if not self.project:
                project = get_project()
                self.env_batch.set_value("PROJECT", project, subgroup=self.job)
        else:
            self.project = None

    def _set_lt_archive_options(self):
        """
        Get the long-term archiver options from $CIMEROOT/cime_config/cesm/machines
        These options will be used when creating the lt_archive run scrip
        """
        lt_archive = LTArchive(self.machine.get_machine_name())
        self.lt_archive_args = lt_archive.get_lt_archive_args()

    def _set_model_run(self):
        """
        set the model run command per machine.
        """
        expect(self.batch_system,
"""
No batch system type configured for this machine!  Please see config_batch.xml
within model's Machines directory, and add a batch system type for this machine
""")

        self.mpirun = self.machine.get_full_mpirun(self, self.case, self.job)

    def _set_batch_directives(self):
        """
        Set the batch directives for this machine from config_batch.xml
        """
        batch_directives = self.batch_parser.get_batch_directives(self)

        self.batchdirectives = "\n".join(batch_directives)

    def transform_vars_bm(self, text, default=None):
        """
        Do the variable substitution for any variables that need transforms
        recursively. Use the self BatchMaker for substitutions
        """
        return transform_vars(text, case=self.case, subgroup=self.job, check_members=self, default=default)

    def make_batch_script(self, input_filename, output_filename):
        expect(os.path.exists(input_filename), "input file '%s' does not exist" % input_filename)

        template_text = self.transform_vars_bm(open(input_filename, "r").read())

        with open(output_filename, "w") as fd:
            fd.write(template_text)

        os.chmod(output_filename, os.stat(output_filename).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    def get_batch_system_type_for_machine(self):
        return self.batch_system

    def set_job(self, job):
        self.job = job
        self._initialize()

# TODO: Machine-specific overrides

def get_batch_maker(job, case=None):
    """
    Simple factory class to get the right BatchMaker class for each machine.
    The only downside to this strategy is that we have to have a BatchMaker_${machine}
    class for every machine we port.
    TODO: REFACTOR this so that if no machine or batch class is found, then the base
    class is returned.
    """
    case = case if case is not None else Case()

    machine = case.get_value("MACH")
    batch_maker = BatchMaker(job, case)
    subclassname = "BatchMaker_%s" % machine
    if "pleiades" in machine:
        new_machine = machine.replace("-", "_")
        subclassname = "BatchMaker_%s" % new_machine

    if subclassname in globals():
        return globals()[subclassname](job, case)
    else:
        return batch_maker


