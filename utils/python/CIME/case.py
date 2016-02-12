"""
Wrapper around all env XML for a case.

All interaction with and between the module files in XML/ takes place
through the Case module.
"""

from XML.standard_module_setup import *
import CIME.utils
from CIME.utils import expect, run_cmd
from CIME.XML.machines import Machines

from CIME.XML.env_test          import EnvTest
from CIME.XML.env_mach_specific import EnvMachSpecific
from CIME.XML.env_case          import EnvCase
from CIME.XML.env_mach_pes      import EnvMachPes
from CIME.XML.env_build         import EnvBuild
from CIME.XML.env_run           import EnvRun
from CIME.XML.env_archive       import EnvArchive
from CIME.XML.env_batch         import EnvBatch

class Case(object):

    def __init__(self, case_root=os.getcwd()):
        expect(os.path.isdir(case_root),
               "Case root directory '%s' does not exist" % case_root)

        self._env_files = []
        self._env_files_that_need_rewrite = set()

        self._env_files.append(EnvTest(case_root))
        self._env_files.append(EnvRun(case_root))
        self._env_files.append(EnvMachSpecific(case_root))
        self._env_files.append(EnvCase(case_root))
        self._env_files.append(EnvMachPes(case_root))
        self._env_files.append(EnvBuild(case_root))
        self._env_files.append(EnvArchive(case_root))
        self._env_files.append(EnvBatch(case_root))

    def __del__(self):
        self.flush()

    def flush(self):
        for env_file in self._env_files_that_need_rewrite:
            env_file.write()

        self._env_files_that_need_rewrite = set()

    def get_value(self, item):
        for env_file in self._env_files:
            result = env_file.get_value(item)
            if (result is not None):
                return self.get_resolved_value(result)

        logging.info("Not able to retreive value for item '%s'" % item)

    def get_resolved_value(self, item):
        # TODO HACK - surely there is a better way?
        num_unresolved = item.count("$")
        if (num_unresolved > 0):
            for env_file in self._env_files:
                result = env_file.get_resolved_value(item)
                if (result.count("$") < num_unresolved):
                    num_unresolved = result.count("$")
                    item = result
                    if ("$" not in item):
                        return item

            logging.warning("Not able to fully resolve item '%s'" % item)

        return item

    def set_value(self, item, value):
        for env_file in self._env_files:
            result = env_file.set_value(item, value)
            if (result is not None):
                self._env_files_that_need_rewrite.add(env_file)
                return result

        logging.warning("Not able to set value for item '%s'" % item)

