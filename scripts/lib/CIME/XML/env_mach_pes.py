"""
Interface to the env_mach_pes.xml file.  This class inherits from EntryID
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.env_base import EnvBase
import math

logger = logging.getLogger(__name__)

class EnvMachPes(EnvBase):

    def __init__(self, case_root=None, infile="env_mach_pes.xml", components=None):
        """
        initialize an object interface to file env_mach_pes.xml in the case directory
        """
        self._components = components
        self._component_value_list = ["NTASKS", "NTHRDS", "NINST",
                                      "ROOTPE", "PSTRID", "NINST_LAYOUT"]
        schema = os.path.join(get_cime_root(), "config", "xml_schemas", "env_mach_pes.xsd")
        EnvBase.__init__(self, case_root, infile, schema=schema)

    def set_value(self, vid, value, subgroup=None, ignore_type=False):
        """
        Set the value of an entry-id field to value
        Returns the value or None if not found
        subgroup is ignored in the general routine and applied in specific methods
        """
        oldcomps = self._components[:]
        nvid, comp, _ = self.check_if_comp_var(vid, None)
        if nvid == "NINST":
            if self.get_value("MULTI_COUPLER"):
                self._components = ["CPL"]
                if comp is None:
                    comp = "CPL"
                expect(comp == "CPL" or value == 1,"Cannot change {} when MULTI_COUPLER flag is TRUE".format(vid))
            elif value != 1:
                if 'CPL' in self._components:
                    self._components.remove('CPL')
                if 'ESP' in self._components and value != 1:
                    self._components.remove('ESP')
                expect(comp is None or comp != "CPL","Cannot change NINST_CPL if MULTI_COUPLER flag is FALSE")
        # Toggling the
        if vid == "MULTI_COUPLER":
            newval = EnvBase.set_value(self, vid, value,subgroup=subgroup, ignore_type=ignore_type)
            if value:
                maxinst = 1
                for tcomp in self._components:
                    if tcomp in ["CPL"]:
                        continue
                    tcomp_inst = "NINST_{}".format(tcomp)
                    maxinst = max(maxinst, self.get_value(tcomp_inst))
                    EnvBase.set_value(self,tcomp_inst, 1)
                    self.set_value("NINST_CPL", maxinst)
            else:
                ninst = self.get_value("NINST_CPL")
                newval = EnvBase.set_value(self, vid, value,subgroup=subgroup, ignore_type=ignore_type)
                cpl_index = None
                esp_index = None
                if 'CPL' in self._components:
                    cpl_index = self._components.index('CPL')
                    del self._components[cpl_index]
                if 'ESP' in self._components:
                    esp_index = self._components.index('ESP')
                    del self._components[esp_index]

                EnvBase.set_value(self, "NINST", ninst)
                if cpl_index is not None:
                    self._components.insert(cpl_index,'CPL')
                if esp_index is not None:
                    self._components.insert(esp_index,'ESP')
                EnvBase.set_value(self, "NINST_CPL", 1)
        else:
            newval = EnvBase.set_value(self, vid, value,subgroup=subgroup, ignore_type=ignore_type)
        self._components = oldcomps
        return newval

    def get_value(self, vid, attribute=None, resolved=True, subgroup=None, pes_per_node=None): # pylint: disable=arguments-differ
        value = EnvBase.get_value(self, vid, attribute, resolved, subgroup)

        if "NTASKS" in vid or "ROOTPE" in vid:
            if pes_per_node is None:
                pes_per_node = self.get_value("PES_PER_NODE")
            if value is not None and value < 0:
                value = -1*value*pes_per_node

        return value

    def get_max_thread_count(self, comp_classes):
        ''' Find the maximum number of openmp threads for any component in the case '''
        max_threads = 1
        for comp in comp_classes:
            threads = self.get_value("NTHRDS",attribute={"component":comp})
            expect(threads is not None, "Error no thread count found for component class {}".format(comp))
            if threads > max_threads:
                max_threads = threads
        return max_threads

    def get_cost_pes(self, totaltasks, max_thread_count, machine=None):
        """
        figure out the value of COST_PES which is the pe value used to estimate model cost
        """
        expect(totaltasks > 0,"totaltasks > 0 expected totaltasks = {}".format(totaltasks))
        pespn = self.get_value("PES_PER_NODE")
        num_nodes, spare_nodes = self.get_total_nodes(totaltasks, max_thread_count)
        num_nodes += spare_nodes
        # This is hardcoded because on yellowstone by default we
        # run with 15 pes per node
        # but pay for 16 pes per node.  See github issue #518
        if machine is not None and machine == "yellowstone":
            pespn = 16
        return num_nodes * pespn

    def get_total_tasks(self, comp_classes):
        total_tasks = 0
        for comp in comp_classes:
            ntasks = self.get_value("NTASKS", attribute={"component":comp})
            rootpe = self.get_value("ROOTPE", attribute={"component":comp})
            pstrid = self.get_value("PSTRID", attribute={"component":comp})
            tt = rootpe + (ntasks - 1) * pstrid + 1
            total_tasks = max(tt, total_tasks)
        if self.get_value("MULTI_COUPLER"):
            total_tasks *= self.get_value("NINST_CPL")
        return total_tasks

    def get_tasks_per_node(self, total_tasks, max_thread_count):
        expect(total_tasks > 0,"totaltasks > 0 expected totaltasks = {}".format(total_tasks))
        tasks_per_node = min(self.get_value("MAX_TASKS_PER_NODE")/ max_thread_count,
                             self.get_value("PES_PER_NODE"), total_tasks)
        return tasks_per_node if tasks_per_node > 0 else 1

    def get_total_nodes(self, total_tasks, max_thread_count):
        """
        Return (num_active_nodes, num_spare_nodes)
        """
        tasks_per_node = self.get_tasks_per_node(total_tasks, max_thread_count)
        num_nodes = int(math.ceil(float(total_tasks) / tasks_per_node))
        return num_nodes, self.get_spare_nodes(num_nodes)

    def get_spare_nodes(self, num_nodes):
        force_spare_nodes = self.get_value("FORCE_SPARE_NODES")
        if force_spare_nodes != -999:
            return force_spare_nodes

        if self.get_value("ALLOCATE_SPARE_NODES"):
            ten_pct = int(math.ceil(float(num_nodes) * 0.1))
            if ten_pct < 1:
                return 1 # Always provide at lease one spare node
            elif ten_pct > 10:
                return 10 # Never provide more than 10 spare nodes
            else:
                return ten_pct
        else:
            return 0
