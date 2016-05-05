"""
BatchUtils class for submitting jobs and managing dependancies.
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import expect, run_cmd, get_model, convert_to_seconds, get_project
from CIME.case import Case
from CIME.batch_maker import get_batch_maker
from CIME.XML.batch import Batch

logger = logging.getLogger(__name__)

class BatchUtils(object):

    def __init__(self, job, case=None, prereq_jobid=None):
        """
        Class constructor.  We need to know where in the filesystem we are,
        so caseroot, case, machroot, machine, cimeroot
        prereq_jobid allows the user to set a prereq to the first job in this workflow
        """
        self.case = case if case is not None else Case()
        self.job = job
        self.batchobj = None
        self.override_node_count = None
        self.batchmaker = None
        self.jobs = []
        self.prereq_jobid = prereq_jobid

    def submit_jobs(self, no_batch=False):
        if no_batch:
            batchtype = "none"
        else:
            batchtype = self.case.get_value("batch_system",subgroup=None)

        if self.batchobj is None:
            self.batchobj = Batch(batch_system=batchtype,
                                  machine=self.case.get_value("MACH"))
        env_batch = self.case._get_env("batch")
        potential_jobs = env_batch.get_jobs()

        for job, jobdict in potential_jobs:
            if job == self.job:
                self.jobs.append((job,None))
            else:
                try:
                    prereq = eval(self.case.get_resolved_value(jobdict['prereq']))
                except:
                    expect(False,"Unable to evaluate prereq expression '%s' for job '%s'"%(jobdict['prereq'],job))
                if prereq:
                    self.jobs.append((job,jobdict['dependancy']))
        depid = {}
        for job, dependancy in self.jobs:
            if dependancy is not None:
                deps = dependancy.split()
            else:
                deps = []
            jobid = ""
            if job == self.job and self.prereq_jobid is not None:
                jobid = self.prereq_jobid
            for dep in deps:
                if dep in depid.keys():
                    jobid += " "+str(depid[dep])
#TODO: doubt these will be used
#               elif dep == "and":
#                   jobid += " && "
#               elif dep == "or":
#                   jobid += " || "


            slen = len(jobid)
            if slen == 0:
                jobid = None

            depid[job] = self.submit_single_job(job, batchtype, jobid)

    def submit_single_job(self, job, batchtype, depid=None):
        caseroot = self.case.get_value("CASEROOT")
        if batchtype == "none":
            logger.info("Starting job script %s"%job)
            run_cmd("%s --caseroot %s"%(os.path.join(".", job),caseroot))
            return

        submitargs = self.get_submit_args(job, batchtype)

        if depid is not None:
            dep_string = self.batchobj.get_value("depend_string")
            dep_string = dep_string.replace("jobid",depid)

            submitargs += " "+dep_string

        batchsubmit = self.case.get_value("BATCHSUBMIT",subgroup=None)
        batchredirect = self.case.get_value("BATCHREDIRECT",subgroup=None)
        submitcmd = batchsubmit + " " + submitargs + " " + batchredirect + \
            " " + job + " --caseroot %s"%caseroot
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
        logger.debug("Submit args are %s"%" ".join(submitargs))
        return submitargs
