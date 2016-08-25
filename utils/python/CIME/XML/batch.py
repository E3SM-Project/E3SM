"""
Interface to the config_batch.xml file.  This class inherits from GenericXML.py

The batch_system type="foo" blocks define most things. Machine-specific overrides
can be defined by providing a batch_system MACH="mach" block.
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.utils import expect, get_cime_root, get_model

logger = logging.getLogger(__name__)

class Batch(GenericXML):

    def __init__(self, batch_system=None, machine=None, infile=None):
        """
        initialize an object
        """
        if infile is None:
            infile = os.path.join(get_cime_root(), "cime_config", get_model(), "machines", "config_batch.xml")

        GenericXML.__init__(self, infile)

        self.batch_system_node = None
        self.machine_node      = None
        self.batch_system      = batch_system
        self.machine           = machine

        #Append the contents of $HOME/.cime/config_batch.xml if it exists
        #This could cause problems if node matchs are repeated when only one is expected
        infile = os.path.join(os.environ.get("HOME"),".cime","config_batch.xml")
        if os.path.exists(infile):
            GenericXML.read(self, infile)

        if self.batch_system is not None:
            self.set_batch_system(self.batch_system, machine=machine)

    def get_batch_system(self):
        """
        Return the name of the batch system
        """
        return self.batch_system

    def get_optional_batch_node(self, nodename, attributes=None):
        """
        Return data on a node for a batch system
        """
        expect(self.batch_system_node is not None, "Batch system not set, use parent get_node?")

        if self.machine_node is not None:
            result = self.get_optional_node(nodename, attributes, root=self.machine_node)
            if result is None:
                return self.get_optional_node(nodename, attributes, root=self.batch_system_node)
            else:
                return result
        else:
            return self.get_optional_node(nodename, attributes, root=self.batch_system_node)

    def set_batch_system(self, batch_system, machine=None):
        """
        Sets the batch system block in the Batch object
        """
        machine = machine if machine is not None else self.machine
        if self.batch_system != batch_system or self.batch_system_node is None:
            nodes = self.get_nodes("batch_system",{"type" : batch_system})
            for node in nodes:
                mach = node.get("MACH")
                if mach is None:
                    self.batch_system_node = node
                elif mach == machine:
                    self.machine = machine
                    self.machine_node = node

            expect(self.batch_system_node is not None, "No batch system '%s' found" % batch_system)

        return batch_system

    def get_value(self, name, attribute=None, resolved=True, subgroup=None):
        """
        Get Value of fields in the config_batch.xml file
        """
        expect(self.batch_system_node is not None, "Batch object has no batch system defined")
        expect(subgroup is None, "This class does not support subgroups")
        value = None

        node = self.get_optional_batch_node(name)
        if node is not None:
            value = node.text

        if value is None:
            # if all else fails
            value = GenericXML.get_value(self, name, attribute, resolved, subgroup)

        if resolved:
            if value is not None:
                value = self.get_resolved_value(value)
            elif name in os.environ:
                value = os.environ[name]

        return value

    def get_batch_jobs(self):
        """
        Return a list of jobs with the first element the name of the case script
        and the second a dict of qualifiers for the job
        """
        jobs = []
        bnode = self.get_node("batch_jobs")
        for jnode in bnode:
            if jnode.tag == "job":
                name = jnode.get("name")
                jdict = {}
                for child in jnode:
                    jdict[child.tag] = child.text

            jobs.append((name, jdict))

        return jobs
