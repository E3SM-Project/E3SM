"""
Interface to the env_mach_pes.xml file.  This class inherits from EnvBase
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.env_base import EnvBase
import math

logger = logging.getLogger(__name__)

class EnvMachPes(EnvBase):

    def __init__(self, case_root=None, infile="env_mach_pes.xml"):
        """
        initialize an object interface to file env_mach_pes.xml in the case directory
        """
        EnvBase.__init__(self, case_root, infile)

    def get_value(self, vid, attribute=None, resolved=True, subgroup=None, pes_per_node=None): # pylint: disable=arguments-differ
        value = EnvBase.get_value(self, vid, attribute, resolved, subgroup)
        if "NTASKS" in vid or "ROOTPE" in vid and pes_per_node is None:
            pes_per_node = self.get_value("PES_PER_NODE")

            if "NTASKS" in vid and value < 0:
                value = -1*value*pes_per_node
            if "ROOTPE" in vid and value < 0:
                value = -1*value*pes_per_node
        return value

    def set_value(self, vid, value, subgroup=None, ignore_type=False, pes_per_node=None): # pylint: disable=arguments-differ
        """
        Set the value of an entry-id field to value
        Returns the value or None if not found
        subgroup is ignored in the general routine and applied in specific methods
        """
        val = None
        node = self.get_optional_node("entry", {"id":vid})
        if node is not None:
            val = self._set_value(node, value, vid, subgroup, ignore_type, pes_per_node=pes_per_node)
        return val



    def _set_value(self, node, value, vid=None, subgroup=None, ignore_type=False, pes_per_node=None): # pylint: disable=arguments-differ
        if vid is None:
            vid = node.get("id")

        if "NTASKS" in vid or "ROOTPE" in vid and pes_per_node is None:
            pes_per_node = self.get_value("PES_PER_NODE")

        if "NTASKS" in vid and value < 0:
            value = -1*value*pes_per_node
        if "ROOTPE" in vid and value < 0:
            value = -1*value*pes_per_node
        val = EnvBase._set_value(self, node, value, vid, subgroup, ignore_type)
        return val

    def get_max_thread_count(self, comp_classes):
        ''' Find the maximum number of openmp threads for any component in the case '''
        max_threads = 1
        for comp in comp_classes:
            if comp == "DRV":
                comp = "CPL"
            threads = self.get_value("NTHRDS_%s"%comp)
            expect(threads is not None, "Error no thread count found for component class %s"%comp)
            if threads > max_threads:
                max_threads = threads
        return max_threads

    def get_cost_pes(self, totaltasks, max_thread_count, machine=None):
        """
        figure out the value of COST_PES which is the pe value used to estimate model cost
        """
        mtpn = self.get_value("MAX_TASKS_PER_NODE")
        pespn = self.get_value("PES_PER_NODE")
        # This is hardcoded because on yellowstone by default we
        # run with 15 pes per node
        # but pay for 16 pes per node.  See github issue #518
        if machine is not None and machine == "yellowstone":
            pespn = 16
        pestot = totaltasks
        if mtpn > pespn and pestot > pespn:
            pestot = pestot * (mtpn // pespn)
            return pespn * self.get_total_nodes(pestot, max_thread_count)
        else:
            # reset cost_pes to totalpes
            return 0

    def get_total_tasks(self, comp_classes):
        total_tasks = 0
        for comp in comp_classes:
            if comp == "DRV":
                comp = "CPL"
            ntasks = self.get_value("NTASKS_%s"%comp)
            rootpe = self.get_value("ROOTPE_%s"%comp)
            pstrid = self.get_value("PSTRID_%s"%comp)
            tt = rootpe + (ntasks - 1) * pstrid + 1
            total_tasks = max(tt, total_tasks)
        if total_tasks < self.get_value("PES_PER_NODE"): # all PEs fit in 1 node
            self.set_value("PES_PER_NODE", total_tasks)
        return total_tasks

    def get_tasks_per_node(self, total_tasks, max_thread_count):
        tasks_per_node = min(self.get_value("MAX_TASKS_PER_NODE")/ max_thread_count,
                             self.get_value("PES_PER_NODE"), total_tasks)
        return tasks_per_node

    def get_total_nodes(self, total_tasks, max_thread_count):
        tasks_per_node = self.get_tasks_per_node(total_tasks, max_thread_count)
        num_nodes = int(math.ceil(float(total_tasks) / tasks_per_node))
        return num_nodes
