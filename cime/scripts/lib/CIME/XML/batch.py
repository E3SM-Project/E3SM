"""
Interface to the config_batch.xml file.  This class inherits from GenericXML.py

The batch_system type="foo" blocks define most things. Machine-specific overrides
can be defined by providing a batch_system MACH="mach" block.
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.XML.files import Files
from CIME.utils import expect

logger = logging.getLogger(__name__)

class Batch(GenericXML):

    def __init__(self, batch_system=None, machine=None, infile=None, files=None, extra_machines_dir=None):
        """
        initialize an object

        If extra_machines_dir is provided, it should be a string giving a path to an
        additional directory that will be searched for a config_batch.xml file; if
        found, the contents of this file will be appended to the standard
        config_batch.xml. An empty string is treated the same as None.
        """
        if files is None:
            files = Files()
        if infile is None:
            infile = files.get_value("BATCH_SPEC_FILE")

        schema = files.get_schema("BATCH_SPEC_FILE")

        GenericXML.__init__(self, infile, schema=schema)

        self.batch_system_node = None
        self.machine_node      = None
        self.batch_system      = batch_system
        self.machine           = machine

        # Append the contents of $HOME/.cime/config_batch.xml if it exists.
        #
        # Also append the contents of a config_batch.xml file in the directory given by
        # extra_machines_dir, if present.
        #
        # This could cause problems if node matches are repeated when only one is expected.
        infile = os.path.join(os.environ.get("HOME"),".cime","config_batch.xml")
        if os.path.exists(infile):
            GenericXML.read(self, infile)
        if extra_machines_dir:
            infile = os.path.join(extra_machines_dir, "config_batch.xml")
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
            result = self.get_optional_child(nodename, attributes, root=self.machine_node)
            if result is None:
                return self.get_optional_child(nodename, attributes, root=self.batch_system_node)
            else:
                return result
        else:
            return self.get_optional_child(nodename, attributes, root=self.batch_system_node)

    def set_batch_system(self, batch_system, machine=None):
        """
        Sets the batch system block in the Batch object
        """
        machine = machine if machine is not None else self.machine
        if self.batch_system != batch_system or self.batch_system_node is None:
            nodes = self.get_children("batch_system",{"type" : batch_system})
            for node in nodes:
                mach = self.get(node, "MACH")
                if mach is None:
                    self.batch_system_node = node
                elif mach == machine:
                    self.machine = machine
                    self.machine_node = node

            expect(self.batch_system_node is not None, "No batch system '{}' found".format(batch_system))

        return batch_system

    #pylint: disable=arguments-differ
    def get_value(self, name, attribute=None, resolved=True, subgroup=None):
        """
        Get Value of fields in the config_batch.xml file
        """
        expect(self.batch_system_node is not None, "Batch object has no batch system defined")
        expect(subgroup is None, "This class does not support subgroups")
        value = None

        node = self.get_optional_batch_node(name)
        if node is not None:
            value = self.text(node)

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
        bnode = self.get_optional_child("batch_jobs")
        if bnode:
            for jnode in self.get_children(root=bnode):
                if self.name(jnode) == "job":
                    name = self.get(jnode, "name")
                    jdict = {}
                    for child in self.get_children(root=jnode):
                        jdict[self.name(child)] = self.text(child)

                    jobs.append((name, jdict))

        return jobs
