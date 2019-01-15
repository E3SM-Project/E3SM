"""
Interface to the env_mach_pes.xml file.  This class inherits from EntryID
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.env_base import EnvBase
import math

logger = logging.getLogger(__name__)

class EnvMachPes(EnvBase):

    def __init__(self, case_root=None, infile="env_mach_pes.xml", components=None, read_only=False):
        """
        initialize an object interface to file env_mach_pes.xml in the case directory
        """
        self._components = components
        schema = os.path.join(get_cime_root(), "config", "xml_schemas", "env_mach_pes.xsd")
        EnvBase.__init__(self, case_root, infile, schema=schema, read_only=read_only)

    def add_comment(self, comment):
        if comment is not None:
            node = self.make_child("comment", text=comment)
            # make_child adds to the end of the file but we want it to follow the header
            # so we need to remove it and add it in the correct position
            self.remove_child(node)
            self.add_child(node, position=1)

    def get_value(self, vid, attribute=None, resolved=True, subgroup=None, max_mpitasks_per_node=None): # pylint: disable=arguments-differ
        # Special variable NINST_MAX is used to determine the number of
        # drivers in multi-driver mode.
        if vid == "NINST_MAX":
            # in the nuopc driver there is only a single NINST value
            value = 1
            for comp in self._components:
                if comp != "CPL":
                    value = max(value, self.get_value("NINST_{}".format(comp)))
            return value

        value = EnvBase.get_value(self, vid, attribute, resolved, subgroup)

        if "NTASKS" in vid or "ROOTPE" in vid:
            if max_mpitasks_per_node is None:
                max_mpitasks_per_node = self.get_value("MAX_MPITASKS_PER_NODE")
            if value is not None and value < 0:
                value = -1*value*max_mpitasks_per_node
        # in the nuopc driver there is only one NINST value
        # so that NINST_{comp} = NINST
        if "NINST_" in vid and value is None:
            value = self.get_value("NINST")
        return value

    def set_value(self, vid, value, subgroup=None, ignore_type=False):
        """
        Set the value of an entry-id field to value
        Returns the value or None if not found
        subgroup is ignored in the general routine and applied in specific methods
        """
        if vid == "MULTI_DRIVER" and value:
            ninst_max = self.get_value("NINST_MAX")
            for comp in self._components:
                if comp == "CPL":
                    continue
                ninst = self.get_value("NINST_{}".format(comp))
                expect(ninst == ninst_max,
                       "All components must have the same NINST value in multi_driver mode.  NINST_{}={} shoud be {}".format(comp,ninst,ninst_max))
        if "NTASKS" in vid or "NTHRDS" in vid:
            expect(value != 0, "Cannot set NTASKS or NTHRDS to 0")


        return EnvBase.set_value(self, vid, value, subgroup=subgroup, ignore_type=ignore_type)


    def get_max_thread_count(self, comp_classes):
        ''' Find the maximum number of openmp threads for any component in the case '''
        max_threads = 1
        for comp in comp_classes:
            threads = self.get_value("NTHRDS",attribute={"compclass":comp})
            expect(threads is not None, "Error no thread count found for component class {}".format(comp))
            if threads > max_threads:
                max_threads = threads
        return max_threads

    def get_total_tasks(self, comp_classes):
        total_tasks = 0
        maxinst = self.get_value("NINST")
        if maxinst:
            comp_interface = "nuopc"
        else:
            comp_interface = 'unknown'
            maxinst = 1
        for comp in comp_classes:
            ntasks = self.get_value("NTASKS", attribute={"compclass":comp})
            rootpe = self.get_value("ROOTPE", attribute={"compclass":comp})
            pstrid = self.get_value("PSTRID", attribute={"compclass":comp})
            if comp != "CPL" and comp_interface!="nuopc":
                ninst = self.get_value("NINST", attribute={"compclass":comp})
                maxinst = max(maxinst, ninst)
            tt = rootpe + (ntasks - 1) * pstrid + 1
            total_tasks = max(tt, total_tasks)
        if self.get_value("MULTI_DRIVER"):
            total_tasks *= maxinst
        return total_tasks

    def get_tasks_per_node(self, total_tasks, max_thread_count):
        expect(total_tasks > 0,"totaltasks > 0 expected, totaltasks = {}".format(total_tasks))
        tasks_per_node = min(self.get_value("MAX_TASKS_PER_NODE")// max_thread_count,
                             self.get_value("MAX_MPITASKS_PER_NODE"), total_tasks)
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
