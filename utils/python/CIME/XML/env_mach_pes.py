"""
Interface to the env_mach_pes.xml file.  This class inherits from EntryID
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.env_base import EnvBase
from CIME.utils import convert_to_type, convert_to_string
import math

logger = logging.getLogger(__name__)

class EnvMachPes(EnvBase):

    def __init__(self, case_root=None, infile="env_mach_pes.xml"):
        """
        initialize an object interface to file env_mach_pes.xml in the case directory
        """
        EnvBase.__init__(self, case_root, infile)
        self._component_value_list = ["NTASKS"]

    def check_if_comp_var(self, vid, attribute=None):
        '''
        Returns the varid and component for
        '''
        comp = None
        if vid in self._component_value_list:
            if attribute is not None:
                if "component" in attribute:
                    comp = attribute["component"]
            return vid, comp
        parts = vid.split("_")
        if len(parts) == 2 and parts[0] in self._component_value_list:
            return parts[0], parts[1]
        return vid, comp



    def get_value(self, vid, attribute=None, resolved=True, subgroup=None, pes_per_node=None): # pylint: disable=arguments-differ
        value = None
        vid, comp = self.check_if_comp_var(vid, attribute)

        if vid in self._component_value_list and comp is None:
            logger.debug("Not enough info to get value for %s"%vid)
            return value
        elif comp is not None:
            if attribute is None:
                attribute = {"component" : comp}
            else:
                attribute["component"] = comp
            node = self.get_optional_node("entry", {"id":vid})
            if node is not None:
                type_str = self._get_type_info(node)
                value = convert_to_type(self.get_element_text("value", attribute, root=node), type_str, vid)

        if value is None:
            value = EnvBase.get_value(self, vid, attribute, resolved, subgroup)

        if "NTASKS" in vid or "ROOTPE" in vid and pes_per_node is None:
            pes_per_node = self.get_value("PES_PER_NODE")

            if "NTASKS" in vid and value < 0:
                value = -1*value*pes_per_node
            if "ROOTPE" in vid and value < 0:
                value = -1*value*pes_per_node
        return value

    def set_value(self, vid, value, subgroup=None, ignore_type=False):
        """
        Set the value of an entry-id field to value
        Returns the value or None if not found
        subgroup is ignored in the general routine and applied in specific methods
        """
        vid, comp = self.check_if_comp_var(vid, None)
        val = None
        node = self.get_optional_node("entry", {"id":vid})
        if node is not None:
            val = self._set_value(node, value, vid, subgroup, ignore_type, component=comp)
        return val

    def _set_value(self, node, value, vid=None, subgroup=None, ignore_type=False, component=None): # pylint: disable=arguments-differ
        if vid is None:
            vid = node.get("id")

        if vid in self._component_value_list:
            attribute = {"component":component}
            type_str = self._get_type_info(node)
            val = self.set_element_text("value", convert_to_string(value, type_str, vid), attribute, root=node)
            return val
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
        return total_tasks

    def get_tasks_per_node(self, total_tasks, max_thread_count):
        tasks_per_node = min(self.get_value("MAX_TASKS_PER_NODE")/ max_thread_count,
                             self.get_value("PES_PER_NODE"), total_tasks)
        return tasks_per_node

    def get_total_nodes(self, total_tasks, max_thread_count):
        tasks_per_node = self.get_tasks_per_node(total_tasks, max_thread_count)
        num_nodes = int(math.ceil(float(total_tasks) / tasks_per_node))
        return num_nodes

    def cleanupnode(self, node):
        """
        Remove the <group>, <file>, <values> and <value> childnodes from node
        """
        fnode = node.find(".//file")
        node.remove(fnode)
        gnode = node.find(".//group")
        node.remove(gnode)
        dnode = node.find(".//default_value")
        if dnode is not None:
            node.remove(dnode)
        if node.get("id") not in self._component_value_list:
            vnode = node.find(".//values")
            if vnode is not None:
                node.remove(vnode)
        return node


    def get_nodes_by_id(self, varid):
        varid, _ = self.check_if_comp_var(varid, attribute)
        return EnvBase.get_nodes_by_id(self, varid)


