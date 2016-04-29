"""
BatchUtils class for submitting jobs and managing dependancies.
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, run_cmd, get_model, convert_to_seconds, get_project
from CIME.case import Case
from CIME.batch_maker import get_batch_maker
#from CIME.XML.env_batch import EnvBatch
from CIME.XML.batch import Batch
#from CIME.XML.machines import Machines
#from CIME.task_maker import TaskMaker
#from CIME.XML.lt_archive import LTArchive
#import re, stat

logger = logging.getLogger(__name__)

class BatchUtils(object):

    def __init__(self, job, case=None):
        """
        Class constructor.  We need to know where in the filesystem we are,
        so caseroot, case, machroot, machine, cimeroot
        """
        self.case = case if case is not None else Case()
        self.job = job
        self.batchobj = None
        self.override_node_count = None
        self.batchmaker = None
        self.dependancies = list()

    def submit_jobs(self, no_batch=False):
        for job in self.dependancies:
            self.submit_single_job(job, no_batch=no_batch)

    def submit_single_job(self, job, no_batch=False):
        if no_batch:
            batchtype = "none"
        else:
            batchtype = self.case.get_value("batch_system",subgroup=None)
        if batchtype == "none":
            logger.info("Starting job script %s"%job)
            run_cmd(os.path.join(".", job))
            return
        if self.batchobj is None:
            self.batchobj = Batch(batch_system=batchtype,
                                  machine=self.case.get_value("MACH"))
        submitargs = self.get_submit_args(job, batchtype)
        batchsubmit = self.case.get_value("BATCHSUBMIT",subgroup=None)
        batchredirect = self.case.get_value("BATCHREDIRECT",subgroup=None)
        submitcmd = batchsubmit + " " + submitargs + " " + batchredirect + \
            " " + job
        logger.info("Submitting job script %s"%submitcmd)
        output = run_cmd(submitcmd)
        jobid = self.get_job_id(batchtype, output)
        logger.debug("Submitted job id is %s"%jobid)
        return jobid

    def get_job_id(self, batchtype, output):

        jobid_pattern = self.batchobj.get_value("jobid_pattern")
        expect(jobid_pattern is not None,"Could not find jobid_pattern in config_batch")
        jobid = re.search(jobid_pattern, output).group(1)
        return jobid

    def get_submit_args(self, job, batchtype):
        if self.batchmaker is None:
            self.batchmaker = get_batch_maker(job, case=self.case)
        task_count = self.case.get_value("task_count", subgroup=job)
        if task_count == "default":
            self.batchmaker.override_node_count = None
        else:
            self.batchmaker.override_node_count = int(task_count)

        self.batchmaker.set_job(job)

        submit_args = self.batchobj.get_submit_args()
        print submit_args
        submitargs=""
        for flag, name in submit_args:
            if name is None:
                submitargs+=" %s"%flag
            else:
                val = self.case.get_value(name,subgroup=job)
                if val is not None and len(val) > 0:
                    if flag.rfind("=", len(flag)-1, len(flag)) >= 0:
                        submitargs+=" %s%s"%(flag,val)
                    else:
                        submitargs+=" %s %s"%(flag,val)
        return submitargs
