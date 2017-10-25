"""
Interface to the env_batch.xml file.  This class inherits from EnvBase
"""

from CIME.XML.standard_module_setup import *
from CIME.utils import format_time
from CIME.XML.env_base import EnvBase
from CIME.utils import transform_vars, get_cime_root, convert_to_seconds

from copy import deepcopy
from collections import OrderedDict
import stat, re, math

logger = logging.getLogger(__name__)

# pragma pylint: disable=attribute-defined-outside-init

class EnvBatch(EnvBase):

    def __init__(self, case_root=None, infile="env_batch.xml"):
        """
        initialize an object interface to file env_batch.xml in the case directory
        """
        self._batchtype = None
        # This arbitrary setting should always be overwritten
        self._default_walltime = "00:20:00"
        schema = os.path.join(get_cime_root(), "config", "xml_schemas", "env_batch.xsd")
        EnvBase.__init__(self, case_root, infile, schema=schema)

    # pylint: disable=arguments-differ
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
                    t_spec = "%H:%M"
                elif value.count(":") == 2:
                    t_spec = "%H:%M:%S"
                else:
                    expect(False, "could not interpret format for wallclock time {}".format(value))
                value = format_time(walltime_format, t_spec, value)

        # allow the user to set item for all jobs if subgroup is not provided
        if subgroup is None:
            nodes = self.get_nodes("entry", {"id":item})
            for node in nodes:
                self._set_value(node, value, vid=item, ignore_type=ignore_type)
                val = value
        else:
            group = self.get_optional_node("group", {"id":subgroup})
            if group is not None:
                node = self.get_optional_node("entry", {"id":item}, root=group)
                if node is not None:
                    val = self._set_value(node, value, vid=item, ignore_type=ignore_type)

        return val

    # pylint: disable=arguments-differ
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
                        "Inconsistent type_info for entry id={} {} {}".format(vid, new_type_info, type_info))
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

    def make_batch_script(self, input_template, job, case):
        expect(os.path.exists(input_template), "input file '{}' does not exist".format(input_template))

        task_count = self.get_value("task_count", subgroup=job)
        overrides = {}
        if task_count is not None:
            overrides["total_tasks"] = int(task_count)
            overrides["num_nodes"]   = int(math.ceil(float(task_count)/float(case.tasks_per_node)))

        overrides["pedocumentation"] = "" # TODO?
        overrides["job_id"] = case.get_value("CASE") + os.path.splitext(job)[1]
        if "pleiades" in case.get_value("MACH"):
            # pleiades jobname needs to be limited to 15 chars
            overrides["job_id"] = overrides["job_id"][:15]

        overrides["batchdirectives"] = self.get_batch_directives(case, job, overrides=overrides)

        output_text = transform_vars(open(input_template,"r").read(), case=case, subgroup=job, overrides=overrides)
        with open(job, "w") as fd:
            fd.write(output_text)
        os.chmod(job, os.stat(job).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    def set_job_defaults(self, batch_jobs, case):
        walltime    = case.get_value("USER_REQUESTED_WALLTIME") if case.get_value("USER_REQUESTED_WALLTIME") else None
        force_queue = case.get_value("USER_REQUESTED_QUEUE") if case.get_value("USER_REQUESTED_QUEUE") else None

        if self._batchtype is None:
            self._batchtype = self.get_batch_system_type()

        if self._batchtype == 'none':
            return

        for job, jsect in batch_jobs:
            task_count = jsect["task_count"] if "task_count" in jsect else None
            if task_count is None:
                node_count = case.num_nodes
            else:
                node_count = int(math.ceil(float(task_count)/float(case.tasks_per_node)))

            if force_queue:
                if not self.queue_meets_spec(force_queue, node_count, walltime=walltime, job=job):
                    logger.warning("WARNING: User-requested queue '{}' does not meet requirements for job '{}'".format(force_queue, job))
                queue = force_queue
            else:
                queue = self.select_best_queue(node_count, walltime=walltime, job=job)
                if queue is None and walltime is not None:
                    # Try to see if walltime was the holdup
                    queue = self.select_best_queue(node_count, walltime=None, job=job)
                    if queue is not None:
                        # It was, override the walltime if a test, otherwise just warn the user
                        new_walltime = self.get_queue_specs(queue)[3]
                        expect(new_walltime is not None, "Should never make it here")
                        logger.warning("WARNING: Requested walltime '{}' could not be matched by any queue".format(walltime))
                        if case.get_value("TEST"):
                            logger.warning("  Using walltime '{}' instead".format(new_walltime))
                            walltime = new_walltime
                        else:
                            logger.warning("  Continuing with suspect walltime, batch submission may fail")

                if queue is None:
                    logger.warning("WARNING: No queue on this system met the requirements for this job. Falling back to defaults")
                    default_queue_node = self.get_default_queue()
                    queue = default_queue_node.text
                    walltime = self.get_queue_specs(queue)[3]

            if walltime is None:
                # Figure out walltime
                specs = self.get_queue_specs(queue)
                if specs is None:
                    # Queue is unknown, use specs from default queue
                    walltime = self.get_default_queue().get("walltimemax")
                else:
                    walltime = specs[3]

                walltime = self._default_walltime if walltime is None else walltime # last-chance fallback

            self.set_value("JOB_QUEUE", queue, subgroup=job)
            self.set_value("JOB_WALLCLOCK_TIME", walltime, subgroup=job)
            logger.debug("Job {} queue {} walltime {}".format(job, queue, walltime))

    def get_batch_directives(self, case, job, overrides=None):
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
                    if default is None:
                        directive = transform_vars(directive, case=case, subgroup=job, default=default, overrides=overrides)
                    else:
                        directive = transform_vars(directive, default=default)

                    result.append("{} {}".format(directive_prefix, directive))

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
            if self._batchtype == "cobalt" and job == "case.st_archive":
                if flag == "-n":
                    name = 'task_count'
                if flag == "--mode":
                    continue

            if name is None:
                submitargs+=" {}".format(flag)
            else:
                if name.startswith("$"):
                    name = name[1:]

                if '$' in name:
                    # We have a complex expression and must rely on get_resolved_value.
                    # Hopefully, none of the values require subgroup
                    val = case.get_resolved_value(name)
                else:
                    val = case.get_value(name, subgroup=job)

                if val is not None and len(str(val)) > 0 and val != "None":
                    # Try to evaluate val if it contains any whitespace
                    if " " in val:
                        try:
                            rval = eval(val)
                        except:
                            rval = val
                    else:
                        rval = val
                    # need a correction for tasks per node
                    if flag == "-n" and rval<= 0:
                        rval = 1

                    if flag == "-q" and rval == "batch" and case.get_value("MACH") == "blues":
                        # Special case. Do not provide '-q batch' for blues
                        continue

                    if flag.rfind("=", len(flag)-1, len(flag)) >= 0 or\
                       flag.rfind(":", len(flag)-1, len(flag)) >= 0:
                        submitargs+=" {}{}".format(flag,str(rval).strip())
                    else:
                        submitargs+=" {} {}".format(flag,str(rval).strip())

        return submitargs

    def submit_jobs(self, case, no_batch=False, job=None, user_prereq=None,
                    skip_pnl=False, mail_user=None, mail_type='never',
                    batch_args=None, dry_run=False):
        alljobs = self.get_jobs()
        startindex = 0
        jobs = []
        firstjob = job
        if job is not None:
            expect(job in alljobs, "Do not know about batch job {}".format(job))
            startindex = alljobs.index(job)

        for index, job in enumerate(alljobs):
            logger.debug( "Index {:d} job {} startindex {:d}".format(index, job, startindex))
            if index < startindex:
                continue
            try:
                prereq = self.get_value('prereq', subgroup=job, resolved=False)
                if prereq is None or job == firstjob or (dry_run and prereq == "$BUILD_COMPLETE"):
                    prereq = True
                else:
                    prereq = case.get_resolved_value(prereq)
                    prereq = eval(prereq)
            except:
                expect(False,"Unable to evaluate prereq expression '{}' for job '{}'".format(self.get_value('prereq',subgroup=job), job))

            if prereq:
                jobs.append((job, self.get_value('dependency', subgroup=job)))

            if self._batchtype == "cobalt":
                break

        depid = OrderedDict()
        jobcmds = []

        for job, dependency in jobs:
            if dependency is not None:
                deps = dependency.split()
            else:
                deps = []
            dep_jobs = []
            if user_prereq is not None:
                dep_jobs.append(user_prereq)
            for dep in deps:
                if dep in depid.keys() and depid[dep] is not None:
                    dep_jobs.append(str(depid[dep]))

            logger.warning("job {} depends on {}".format(job, dep_jobs))
            result = self._submit_single_job(case, job,
                                             dep_jobs=dep_jobs,
                                             no_batch=no_batch,
                                             skip_pnl=skip_pnl,
                                             mail_user=mail_user,
                                             mail_type=mail_type,
                                             batch_args=batch_args,
                                             dry_run=dry_run)
            batch_job_id = str(alljobs.index(job)) if dry_run else result
            depid[job] = batch_job_id
            jobcmds.append( (job, result) )
            if self._batchtype == "cobalt":
                break

        if dry_run:
            return jobcmds
        else:
            return depid

    def _submit_single_job(self, case, job, dep_jobs=None, no_batch=False,
                           skip_pnl=False, mail_user=None, mail_type='never',
                           batch_args=None, dry_run=False):
        logger.warning("Submit job {}".format(job))
        batch_system = self.get_value("BATCH_SYSTEM", subgroup=None)
        if batch_system is None or batch_system == "none" or no_batch:
            # Import here to avoid circular include
            from CIME.case_test       import case_test # pylint: disable=unused-variable
            from CIME.case_run        import case_run # pylint: disable=unused-variable
            from CIME.case_st_archive import case_st_archive # pylint: disable=unused-variable

            logger.info("Starting job script {}".format(job))

            function_name = job.replace(".", "_")
            if not dry_run:
                locals()[function_name](case)

            return

        submitargs = self.get_submit_args(case, job)
        args_override = self.get_value("BATCH_COMMAND_FLAGS", subgroup=job)
        if args_override:
            submitargs = args_override

        if dep_jobs is not None and len(dep_jobs) > 0:
            logger.info("dependencies: {}".format(dep_jobs))
            dep_string = self.get_value("depend_string", subgroup=None)
            separator_string = self.get_value("depend_separator", subgroup=None)
            expect("jobid" in dep_string, "depend_string is missing jobid for prerequisite jobs")
            dep_ids_str = str(dep_jobs[0])
            for dep_id in dep_jobs[1:]:
                dep_ids_str += separator_string + str(dep_id)
            dep_string = dep_string.replace("jobid",dep_ids_str.strip()) # pylint: disable=maybe-no-member
            submitargs += " " + dep_string

        if batch_args is not None:
            submitargs += " " + batch_args

        if mail_user is not None:
            mail_user_flag = self.get_value('batch_mail_flag', subgroup=None)
            if mail_user_flag is not None:
                submitargs += " " + mail_user_flag + " " + mail_user
        if 'never' not in mail_type:
            mail_type_flag, mail_type = self.get_batch_mail_type(mail_type)
            if mail_type_flag is not None:
                submitargs += " " + mail_type_flag + " " + mail_type

        batchsubmit = self.get_value("batch_submit", subgroup=None)
        expect(batchsubmit is not None,
               "Unable to determine the correct command for batch submission.")
        batchredirect = self.get_value("batch_redirect", subgroup=None)
        submitcmd = ''
        for string in (batchsubmit, submitargs, batchredirect, job):
            if  string is not None:
                submitcmd += string + " "

        if job == 'case.run' and skip_pnl:
            batch_env_flag = self.get_value("batch_env", subgroup=None)
            if not batch_env_flag:
                submitcmd += " --skip-preview-namelist"
            else:
                submitcmd += " {} ARGS_FOR_SCRIPT='--skip-preview-namelist'".format(batch_env_flag)

        if dry_run:
            return submitcmd
        else:
            logger.info("Submitting job script {}".format(submitcmd))
            output = run_cmd_no_fail(submitcmd, combine_output=True)
            jobid = self.get_job_id(output)
            logger.info("Submitted job id is {}".format(jobid))
            return jobid

    def get_batch_mail_type(self, mail_type='never'):
        mail_type_flag = self.get_value("batch_mail_type_flag", subgroup=None)
        raw =  self.get_value("batch_mail_type", subgroup=None)
        mail_types = [item.strip() for item in raw.split(",")] # pylint: disable=no-member
        idx = ["never", "all", "begin", "end", "fail"].index(mail_type)

        return mail_type_flag, mail_types[idx]

    def get_batch_system_type(self):
        nodes = self.get_nodes("batch_system")
        for node in nodes:
            type_ = node.get("type")
            if type_ is not None:
                self._batchtype = type_
        return self._batchtype

    def set_batch_system_type(self, batchtype):
        self._batchtype = batchtype

    def get_job_id(self, output):
        jobid_pattern = self.get_value("jobid_pattern", subgroup=None)
        expect(jobid_pattern is not None, "Could not find jobid_pattern in env_batch.xml")
        search_match = re.search(jobid_pattern, output)
        expect(search_match is not None,
               "Couldn't match jobid_pattern '{}' within submit output:\n '{}'".format(jobid_pattern, output))
        jobid = search_match.group(1)
        return jobid

    def queue_meets_spec(self, queue, num_nodes, walltime=None, job=None):
        specs = self.get_queue_specs(queue)
        if specs is None:
            logger.warning("WARNING: queue '{}' is unknown to this system".format(queue))
            return True

        nodemin, nodemax, jobname, walltimemax, strict = specs

        # A job name match automatically meets spec
        if job is not None and jobname is not None:
            return jobname == job

        if nodemin is not None and num_nodes < int(nodemin) or \
           nodemax is not None and num_nodes > int(nodemax):
            return False

        if walltime is not None and walltimemax is not None and strict:
            walltime_s = convert_to_seconds(walltime)
            walltimemax_s = convert_to_seconds(walltimemax)
            if walltime_s > walltimemax_s:
                return False

        return True

    def select_best_queue(self, num_nodes, walltime=None, job=None):
        # Make sure to check default queue first.
        all_queues = []
        all_queues.append( self.get_default_queue())
        all_queues = all_queues + self.get_all_queues()
        for queue in all_queues:
            if queue is not None:
                qname = queue.text
                if self.queue_meets_spec(qname, num_nodes, walltime=walltime, job=job):
                    return qname

        return None

    def get_queue_specs(self, queue):
        """
        Get queue specifications by name.

        Returns (nodemin, nodemax, jobname, walltimemax, is_strict)
        """
        for queue_node in self.get_all_queues():
            if queue_node.text == queue:
                nodemin = queue_node.get("nodemin")
                nodemax = queue_node.get("nodemax")
                jobname = queue_node.get("jobname")
                walltimemax = queue_node.get("walltimemax")
                strict = queue_node.get("strict") == "true"

                return nodemin, nodemax, jobname, walltimemax, strict

        return None

    def get_default_queue(self):
        node = self.get_optional_node("queue", attributes={"default" : "true"})
        if node is None:
            node = self.get_optional_node("queue")
        expect(node is not None, "No queues found")
        return node

    def get_all_queues(self):
        return self.get_nodes("queue")

    def get_nodes(self, nodename, attributes=None, root=None, xpath=None):
        if nodename in ("JOB_WALLCLOCK_TIME", "PROJECT", "CHARGE_ACCOUNT", 
                        "PROJECT_REQUIRED", "JOB_QUEUE", "BATCH_COMMAND_FLAGS"):
            nodes = EnvBase.get_nodes(self, "entry", attributes={"id":nodename},
                                        root=root, xpath=xpath)
        else:
            nodes =  EnvBase.get_nodes(self, nodename, attributes, root, xpath)
        return nodes

    def get_status(self, jobid):
        batch_query = self.get_optional_node("batch_query")
        if batch_query is None:
            logger.warning("Batch queries not supported on this platform")
        else:
            cmd = batch_query.text + " "
            if "per_job_arg" in batch_query.attrib:
                cmd += batch_query.get("per_job_arg") + " "

            cmd += jobid

            status, out, err = run_cmd(cmd)
            if status != 0:
                logger.warning("Batch query command '{}' failed with error '{}'".format(cmd, err))
            else:
                return out.strip()

    def cancel_job(self, jobid):
        batch_cancel = self.get_optional_node("batch_cancel")
        if batch_cancel is None:
            logger.warning("Batch cancellation not supported on this platform")
            return False
        else:
            cmd = batch_cancel.text + " "  + str(jobid)

            status, out, err = run_cmd(cmd)
            if status != 0:
                logger.warning("Batch cancel command '{}' failed with error '{}'".format(cmd, out + "\n" + err))
            else:
                return True
