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

from XML.standard_module_setup import *
from CIME.utils import expect, run_cmd, get_model
from CIME.case import Case
from CIME.XML.env_batch import EnvBatch
from CIME.XML.batch import Batch
from CIME.XML.machines import Machines

import re

logger = logging.getLogger(__name__)

class BatchMaker(object):

    def __init__(self, case=None):
        """
        Class constructor.  We need to know where in the filesystem we are,
        so caseroot, case, machroot, machine, cimeroot
        """
        self.case = case if case is not None else Case()

        # set up paths to the template files, this could and should be
        # extracted out somehow??
        self.job_id = case.get_value("CASE")
        if "pleiades" in case.get_value("MACH"):
            # pleiades jobname needs to be limited to 15 chars
            self.job_id = self.job_id[:15]

        self.output_error_path = case.get_value("CASE")

        self.env_batch = EnvBatch()

        self.batch_parser = Batch()
        self.config_machines_parser = Machines(machine=case.get_value("MACH"))
        self.batch_system = self.config_machines_parser.get_batch_system_type()

        self._set_task_info()
        self._set_queue()
        self._set_wall_time()
        self._set_project()
        self._set_lt_archive_options()
        self._set_model_run()
        self._set_batch_directives()

    def _set_task_info(self):
        

    def _transform_vars(self, text):
        """
        Do the variable substitution for any variables that need transforms
        recursively.
        """
        directive_re = re.compile(r"{{ (\w+) }}")
        lines = text.splitlines()
        for line in lines:
            # loop through directive line, replacing each string enclosed with
            # template characters with the necessary values.
            while directive_re.search(line):
                m = directive_re.search(line)
                variable = m.groups()[0]
                whole_match = m.group()

                repl = case.get_value(variable.upper())
                if repl is not None:
                    line = line.replace(whole_match, repl)
                elif hasattr(self, variable.lower()):
                    line = line.replace(whole_match, getattr(self, variable.lower()))
                else:
                    logger.warn("Could not replace variable '%s'" % variable)

        return "\n".join(lines)

    def make_batch_script(self, input_filename, output_filename):
        expect(os.path.exists(input_filename), "input file '%s' does not exist" % input_filename)

        template_text = transform_vars(open(input_filename, "r").read())

        with open(output_filename, "w") as fd:
            fd.write(template_text)

        os.chmod(output_filename, os.stat(output_filename) | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    def get_batch_system_type_for_machine(self):
        return self.batch_system

    def get_batch_directives(self):
        # TODO
        pass

    def get_field(self, field_name):


def get_batch_maker(case=None):
    """
    Simple factory class to get the right BatchMaker class for each machine.
    The only downside to this strategy is that we have to have a BatchMaker_${machine}
    class for every machine we port.
    TODO: REFACTOR this so that if no machine or batch class is found, then the base
    class is returned.
    """
    case = case if case is not None else Case()

    machine = case.get_value("MACH")
    batch_maker = BatchMaker(case)
    subclassname = "BatchMaker_%s" % machine
    if "pleiades" in machine:
        new_machine = machine.replace("-", "_")
        subclassname = "BatchMaker_%s" % new_machine

    if subclassname in globals():
        return globals[subclassname](case)
    else:
        return batch_maker
