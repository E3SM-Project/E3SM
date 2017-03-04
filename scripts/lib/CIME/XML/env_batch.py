"""
Interface to the env_batch.xml file.  This class inherits from EnvBase
"""
import stat
import time
import re
import math
from CIME.XML.standard_module_setup import *
from CIME.XML.env_base import EnvBase
from CIME.utils import transform_vars, get_cime_root
from copy import deepcopy

logger = logging.getLogger(__name__)

# pragma pylint: disable=attribute-defined-outside-init

class EnvBatch(EnvBase):

    def __init__(self, case_root=None, infile="env_batch.xml"):
        """
        initialize an object interface to file env_batch.xml in the case directory
        """
        self.prereq_jobid = None
        self.batchtype = None
        # This arbitrary setting should always be overwritten
        self._default_walltime = "00:20:00"
        EnvBase.__init__(self, case_root, infile)

    def set_value(self, item, value, subgroup=None, ignore_type=False):
        """
        Override the entry_id set_value function with some special cases for this class
        """
        val = None
        if item == "JOB_WALLCLOCK_TIME":
            #Most systems use %H:%M:%S format for wallclock but LSF
            #uses %H:%M this code corrects the value passed in to be
            #the correct format - if we find we have more exceptions
            #than this we may need to generalize this further
            walltime_format = self.get_value("walltime_format", subgroup=None)
            if walltime_format is not None and walltime_format.count(":") != value.count(":"): # pylint: disable=maybe-no-member
                if value.count(":") == 1:
                    t = time.strptime(value,"%H:%M")
                elif value.count(":") == 2:
                    t = time.strptime(value,"%H:%M:%S")
                else:
                    expect(False, "could not interpret format for wallclock time %s"%value)

                value = time.strftime(walltime_format, t)

        # allow the user to set item for all jobs if subgroup is not provided
        if subgroup is None:
            nodes = self.get_nodes("entry", {"id":item})
            for node in nodes:
                self._set_value(node, value, vid=item, ignore_type=ignore_type)
                val = value
        else:
            val = self.set_value(item, value, subgroup=subgroup, ignore_type=ignore_type)

        return val

    def get_value(self, item, attribute=None, resolved=True, subgroup="case.run"):
        """
        Must default subgroup to something in order to provide single return value
        """

        value = None
        if subgroup is None:
            nodes = self.get_nodes(item, attribute)
            if len(nodes) == 1:
                node = nodes[0]
                value = node.text
                if resolved:
                    value = self.get_resolved_value(value)
            elif not nodes:
                value = EnvBase.get_value(self,item,attribute,resolved)
        else:
            value = EnvBase.get_value(self, item, attribute=attribute, resolved=resolved, subgroup=subgroup)

        return value

    def get_type_info(self, vid):
        nodes = self.get_nodes("entry",{"id":vid})
        type_info = None
        for node in nodes:
            new_type_info = self._get_type_info(node)
            if type_info is None:
                type_info = new_type_info
            else:
                expect( type_info == new_type_info,
                        "Inconsistent type_info for entry id=%s %s %s" % (vid, new_type_info, type_info))
        return type_info

    def get_jobs(self):
        groups = self.get_nodes("group")
        results = []
        for group in groups:
            if group.get("id") not in ["job_submission", "config_batch"]:
                results.append(group.get("id"))

        return results

    def create_job_groups(self, batch_jobs):
        # Subtle: in order to support dynamic batch jobs, we need to remove the
        # job_submission group and replace with job-based groups

        orig_group = self.get_optional_node("group", {"id":"job_submission"})
        expect(orig_group, "Looks like job groups have already been created")

        childnodes = []
        for child in reversed(orig_group):
            childnodes.append(deepcopy(child))
            orig_group.remove(child)

        self.root.remove(orig_group)

        for name, jdict in batch_jobs:
            new_job_group = ET.Element("group")
            new_job_group.set("id", name)
            for field in jdict.keys():
                val = jdict[field]
                node = ET.SubElement(new_job_group, "entry", {"id":field,"value":val})
                tnode = ET.SubElement(node, "type")
                tnode.text = "char"

            for child in childnodes:
                new_job_group.append(deepcopy(child))

            self.root.append(new_job_group)

    def cleanupnode(self, node):
        if node.get("id") == "batch_system":
            fnode = node.find(".//file")
            node.remove(fnode)
            gnode = node.find(".//group")
            node.remove(gnode)
            vnode = node.find(".//values")
            if vnode is not None:
                node.remove(vnode)
        else:
            node = EnvBase.cleanupnode(self, node)
        return node

    def set_batch_system(self, batchobj, batch_system_type=None):
        if batch_system_type is not None:
            self.set_batch_system_type(batch_system_type)
        if batchobj.batch_system_node is not None:
            self.root.append(deepcopy(batchobj.batch_system_node))
        if batchobj.machine_node is not None:
            self.root.append(deepcopy(batchobj.machine_node))

    def make_batch_script(self, input_template, job, case, total_tasks, tasks_per_node, num_nodes, thread_count):
        expect(os.path.exists(input_template), "input file '%s' does not exist" % input_template)

        self.tasks_per_node = tasks_per_node
        self.num_tasks = total_tasks
        self.tasks_per_numa = tasks_per_node / 2
        self.thread_count = thread_count
        task_count = self.get_value("task_count", subgroup=job)

        if task_count == "default":
            self.total_tasks = total_tasks
            self.num_nodes = num_nodes
        else:
            self.total_tasks = task_count
            self.num_nodes = int(math.ceil(float(task_count)/float(tasks_per_node)))

        self.pedocumentation = ""
        self.job_id = case.get_value("CASE") + os.path.splitext(job)[1]
        if "pleiades" in case.get_value("MACH"):
            # pleiades jobname needs to be limited to 15 chars
            self.job_id = self.job_id[:15]
        self.output_error_path = self.job_id

        self.batchdirectives = self.get_batch_directives(case, job)

        output_text = transform_vars(open(input_template,"r").read(), case=case, subgroup=job, check_members=self)
        with open(job, "w") as fd:
            fd.write(output_text)
        os.chmod(job, os.stat(job).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    def set_job_defaults(self, batch_jobs, pesize=None, walltime=None, force_queue=None):
        if self.batchtype is None:
            self.batchtype = self.get_batch_system_type()

        if self.batchtype == 'none':
            return

        for job, jsect in batch_jobs:
            task_count = jsect["task_count"]
            if task_count is None or task_count == "default":
                task_count = pesize
            else:
                task_count = int(task_count)

            queue = force_queue if force_queue is not None else self.select_best_queue(task_count, job)
            self.set_value("JOB_QUEUE", queue, subgroup=job)

            walltime = self.get_max_walltime(queue) if walltime is None else walltime
            if walltime is None:
                logger.warn("Could not find a queue matching task count %d, falling back to deprecated default walltime parameter"%task_count)
                #if the user names a queue which is not defined in config_batch.xml and does not set a
                #walltime, fall back to the max walltime in the default queue
                if force_queue:
                    self.get_default_queue()
                walltime = self._default_walltime

            self.set_value("JOB_WALLCLOCK_TIME", walltime, subgroup=job)
            logger.debug("Job %s queue %s walltime %s" % (job, queue, walltime))

    def get_batch_directives(self, case, job, raw=False):
        """
        """
        result = []
        directive_prefix = self.get_node("batch_directive").text
        directive_prefix = "" if directive_prefix is None else directive_prefix

        roots = self.get_nodes("batch_system")
        for root in roots:
            if root is not None:
                nodes = self.get_nodes("directive", root=root)
                for node in nodes:
                    directive = self.get_resolved_value("" if node.text is None else node.text)
                    default = node.get("default")
                    if not raw:
                        directive = transform_vars(directive, case=case, subgroup=job, default=default, check_members=self)
                    elif default is not None:
                        directive = transform_vars(directive, default=default)
                    result.append("%s %s" % (directive_prefix, directive))

        return "\n".join(result)

    def get_submit_args(self, case, job):
        '''
        return a list of touples (flag, name)
        '''
        submitargs = " "
        bs_nodes = self.get_nodes("batch_system")
        submit_arg_nodes = []
        for node in bs_nodes:
            submit_arg_nodes += self.get_nodes("arg",root=node)
        for arg in submit_arg_nodes:
            flag = arg.get("flag")
            name = arg.get("name")
            if self.batchtype == "cobalt" and job == "case.st_archive":
                if flag == "-n":
                    name = 'task_count'
                if flag == "--mode":
                    continue

            if name is None:
                submitargs+=" %s"%flag
            else:
                if name.startswith("$"):
                    name = name[1:]
                val = case.get_value(name, subgroup=job)
                if val is None:
                    val = case.get_resolved_value(name)

                if val is not None and len(str(val)) > 0 and val != "None":
                    # Try to evaluate val
                    try:
                        rval = eval(val)
                    except:
                        rval = val
                    # need a correction for tasks per node
                    if flag == "-n" and rval<= 0:
                        rval = 1

                    if flag.rfind("=", len(flag)-1, len(flag)) >= 0 or\
                       flag.rfind(":", len(flag)-1, len(flag)) >= 0:
                        submitargs+=" %s%s"%(flag,str(rval).strip())
                    else:
                        submitargs+=" %s %s"%(flag,str(rval).strip())

        return submitargs

    def submit_jobs(self, case, no_batch=False, job=None, batch_args=None):
        alljobs = self.get_jobs()
        startindex = 0
        jobs = []
        if job is not None:
            expect(job in alljobs, "Do not know about batch job %s"%job)
            startindex = alljobs.index(job)

        for index, job in enumerate(alljobs):
            logger.debug( "Index %d job %s startindex %d" % (index, job, startindex))
            if index < startindex:
                continue
            try:
                prereq = self.get_value('prereq', subgroup=job, resolved=False)
                if prereq is None:
                    prereq = True
                else:
                    prereq = case.get_resolved_value(prereq)
                    prereq = eval(prereq)
            except:
                expect(False,"Unable to evaluate prereq expression '%s' for job '%s'"%(self.get_value('prereq',subgroup=job), job))
            if prereq:
                jobs.append((job,self.get_value('dependency', subgroup=job)))
            if self.batchtype == "cobalt":
                break

        depid = {}
        for job, dependency in jobs:
            if dependency is not None:
                deps = dependency.split()
            else:
                deps = []
            jobid = ""
            if self.prereq_jobid is not None:
                jobid = self.prereq_jobid
            for dep in deps:
                if dep in depid.keys() and depid[dep] is not None:
                    jobid += " "+str(depid[dep])
#TODO: doubt these will be used
#               elif dep == "and":
#                   jobid += " && "
#               elif dep == "or":
#                   jobid += " || "


            slen = len(jobid)
            if slen == 0:
                jobid = None

            logger.warn("job is %s"%job)
            depid[job] = self.submit_single_job(case, job, jobid, no_batch=no_batch, batch_args=batch_args)
            if self.batchtype == "cobalt":
                break

        return sorted(list(depid.values()))

    def submit_single_job(self, case, job, depid=None, no_batch=False, batch_args=None):
        logger.warn("Submit job %s"%job)
        caseroot = case.get_value("CASEROOT")
        batch_system = self.get_value("BATCH_SYSTEM", subgroup=None)
        if batch_system is None or batch_system == "none" or no_batch:
            # Import here to avoid circular include
            from CIME.case_test       import case_test # pylint: disable=unused-variable
            from CIME.case_run        import case_run # pylint: disable=unused-variable
            from CIME.case_st_archive import case_st_archive # pylint: disable=unused-variable
            from CIME.case_lt_archive import case_lt_archive # pylint: disable=unused-variable

            logger.info("Starting job script %s" % job)

            # Hack until all testcases are ported to python
            testcase = case.get_value("TESTCASE")
            cimeroot = get_cime_root()
            testscript = os.path.join(cimeroot, "scripts", "Testing", "Testcases", "%s_script" % testcase)
            if job == "case.test" and testcase is not None and os.path.exists(testscript):
                run_cmd_no_fail("%s --caseroot %s" % (os.path.join(".", job), caseroot))
            else:
                # This is what we want longterm
                function_name = job.replace(".", "_")
                success = locals()[function_name](case)
                expect(success, "%s failed" % function_name)

            return

        submitargs = self.get_submit_args(case, job)

        if depid is not None:
            dep_string = self.get_value("depend_string", subgroup=None)
            dep_string = dep_string.replace("jobid",depid.strip()) # pylint: disable=maybe-no-member
            submitargs += " " + dep_string

        if batch_args is not None:
            submitargs += " " + batch_args

        batchsubmit = self.get_value("batch_submit", subgroup=None)
        expect(batchsubmit is not None,
               "Unable to determine the correct command for batch submission.")
        batchredirect = self.get_value("batch_redirect", subgroup=None)
        submitcmd = ''
        for string in (batchsubmit, submitargs, batchredirect, job):
            if  string is not None:
                submitcmd += string + " "

        logger.info("Submitting job script %s"%submitcmd)
        output = run_cmd_no_fail(submitcmd)
        jobid = self.get_job_id(output)
        logger.info("Submitted job id is %s"%jobid)
        return jobid

    def get_batch_system_type(self):
        nodes = self.get_nodes("batch_system")
        for node in nodes:
            type_ = node.get("type")
            if type_ is not None:
                self.batchtype = type_
        return self.batchtype

    def set_batch_system_type(self, batchtype):
        self.batchtype = batchtype

    def get_job_id(self, output):
        jobid_pattern = self.get_value("jobid_pattern", subgroup=None)
        expect(jobid_pattern is not None, "Could not find jobid_pattern in env_batch.xml")
        search_match = re.search(jobid_pattern, output)
        expect(search_match is not None,
               "Couldn't match jobid_pattern '%s' within submit output:\n '%s'" % (jobid_pattern, output))
        jobid = search_match.group(1)
        return jobid

    def select_best_queue(self, num_pes, job=None):
        # Make sure to check default queue first.
        all_queues = []
        all_queues.append( self.get_default_queue())
        all_queues = all_queues + self.get_all_queues()
        for queue in all_queues:
            if queue is not None:
                jobmin = queue.get("jobmin")
                jobmax = queue.get("jobmax")
                jobname = queue.get("jobname")
                if jobname is not None:
                    if job == jobname:
                        return queue.text
                # if the fullsum is between the min and max # jobs, then use this queue.
                elif jobmin is not None and jobmax is not None and num_pes >= int(jobmin) and num_pes <= int(jobmax):
                    return queue.text
        return None

    def get_max_walltime(self, queue):
        for queue_node in self.get_all_queues():
            if queue_node.text == queue:
                return queue_node.get("walltimemax")

    def get_default_queue(self):
        node = self.get_optional_node("queue", attributes={"default" : "true"})
        if node is None:
            node = self.get_optional_node("queue")
        expect(node is not None, "No queues found")
        self._default_walltime = node.get("walltimemax")
        return(node)

    def get_all_queues(self):
        return self.get_nodes("queue")

    def get_nodes(self, nodename, attributes=None, root=None, xpath=None):
        if nodename in ("JOB_WALLCLOCK_TIME", "PROJECT", "PROJECT_REQUIRED",
                        "JOB_QUEUE"):
            nodes = EnvBase.get_nodes(self, "entry", attributes={"id":nodename},
                                        root=root, xpath=xpath)
        else:
            nodes =  EnvBase.get_nodes(self, nodename, attributes, root, xpath)
        return nodes
