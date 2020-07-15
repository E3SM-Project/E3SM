"""
Interface to the env_workflow.xml file.  This class inherits from EnvBase
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.env_base import EnvBase
from CIME.utils import get_cime_root
import re, math

logger = logging.getLogger(__name__)

# pragma pylint: disable=attribute-defined-outside-init

class EnvWorkflow(EnvBase):

    def __init__(self, case_root=None, infile="env_workflow.xml", read_only=False):
        """
        initialize an object interface to file env_workflow.xml in the case directory
        """
        # This arbitrary setting should always be overwritten
        #        schema = os.path.join(get_cime_root(), "config", "xml_schemas", "env_workflow.xsd")
        # TODO: define schema for this file
        schema = None
        super(EnvWorkflow,self).__init__(case_root, infile, schema=schema, read_only=read_only)

    def create_job_groups(self, batch_jobs, is_test):
        # Subtle: in order to support dynamic batch jobs, we need to remove the
        # job_submission group and replace with job-based groups
        orig_group = self.get_optional_child("group", {"id":"job_submission"},
                                    err_msg="Looks like job groups have already been created")
        expect(orig_group, "No workflow groups found")
        orig_group_children = super(EnvWorkflow, self).get_children(root=orig_group)

        childnodes = []
        for child in reversed(orig_group_children):
            childnodes.append(child)

        self.remove_child(orig_group)

        for name, jdict in batch_jobs:
            if name == "case.run" and is_test:
                pass # skip
            elif name == "case.test" and not is_test:
                pass # skip
            elif name == "case.run.sh":
                pass # skip
            else:
                new_job_group = self.make_child("group", {"id":name})
                for field in jdict.keys():
                    if field == "runtime_parameters":
                        continue
                    val = jdict[field]
                    node = self.make_child("entry", {"id":field,"value":val}, root=new_job_group)
                    self.make_child("type", root=node, text="char")

                for child in childnodes:
                    self.add_child(self.copy(child), root=new_job_group)

    def get_jobs(self):
        groups = self.get_children("group")
        results = []
        for group in groups:
            results.append(self.get(group, "id"))
        return results

    def get_type_info(self, vid):
        gnodes = self.get_children("group")
        type_info = None
        for gnode in gnodes:
            nodes = self.get_children("entry",{"id":vid}, root=gnode)
            type_info = None
            for node in nodes:
                new_type_info = self._get_type_info(node)
                if type_info is None:
                    type_info = new_type_info
                else:
                    expect( type_info == new_type_info,
                            "Inconsistent type_info for entry id={} {} {}".format(vid, new_type_info, type_info))
        return type_info

    def get_job_specs(self, case, job):
        task_count = case.get_resolved_value(self.get_value("task_count", subgroup=job))
        tasks_per_node = case.get_resolved_value(self.get_value("tasks_per_node", subgroup=job))
        thread_count = case.get_resolved_value(self.get_value("thread_count", subgroup=job))
        num_nodes = None
        if task_count is not None and tasks_per_node is not None:
            task_count = int(task_count)
            num_nodes   = int(math.ceil(float(task_count)/float(tasks_per_node)))
            tasks_per_node =  task_count//num_nodes
        if not thread_count:
            thread_count = 1

        return task_count, num_nodes, tasks_per_node, thread_count

    # pylint: disable=arguments-differ
    def get_value(self, item, attribute=None, resolved=True, subgroup="PRIMARY"):
        """
        Must default subgroup to something in order to provide single return value
        """
        value = None
        if subgroup == "PRIMARY":
            subgroup = "case.test" if "case.test" in self.get_jobs() else "case.run"

        #pylint: disable=assignment-from-none
        if value is None:
            value = super(EnvWorkflow, self).get_value(item, attribute=attribute, resolved=resolved, subgroup=subgroup)

        return value

    # pylint: disable=arguments-differ
    def set_value(self, item, value, subgroup=None, ignore_type=False):
        """
        Override the entry_id set_value function with some special cases for this class
        """
        val = None

        # allow the user to set item for all jobs if subgroup is not provided
        if subgroup is None:
            gnodes = self.get_children("group")
            for gnode in gnodes:
                node = self.get_optional_child("entry", {"id":item}, root=gnode)
                if node is not None:
                    self._set_value(node, value, vid=item, ignore_type=ignore_type)
                    val = value
        else:
            group = self.get_optional_child("group", {"id":subgroup})
            if group is not None:
                node = self.get_optional_child("entry", {"id":item}, root=group)
                if node is not None:
                    val = self._set_value(node, value, vid=item, ignore_type=ignore_type)

        return val

    def get_children(self, name=None, attributes=None, root=None):
        if name in ("JOB_WALLCLOCK_TIME", "PROJECT", "CHARGE_ACCOUNT",
                    "JOB_QUEUE", "BATCH_COMMAND_FLAGS"):
            nodes = super(EnvWorkflow, self).get_children("entry", attributes={"id":name}, root=root)
        else:
            nodes = super(EnvWorkflow, self).get_children(name, attributes=attributes, root=root)

        return nodes
