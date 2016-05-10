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
            self.batch_system_node = self.get_optional_node("batch_system", {"type" : batch_system})
            expect(self.batch_system_node is not None, "No batch system '%s' found" % batch_system)

            if machine is not None:
                self.machine_node = self.get_optional_node("batch_system", {"MACH" : machine})
                self.machine = machine

        return batch_system

    def get_value(self, name, resolved=True):
        """
        Get Value of fields in the config_batch.xml file
        """
        expect(self.batch_system_node is not None, "Batch object has no batch system defined")
        value = None

        node = self.get_optional_batch_node(name)
        if node is not None:
            value = node.text

        if value is None:
            # if all else fails
            value = GenericXML.get_value(self, name)

        if resolved:
            if value is not None:
                value = self.get_resolved_value(value)
            elif name in os.environ:
                value = os.environ[name]

        return value

    def get_batch_directives(self, batch_maker=None):
        """
        """
        result = []
        directive_prefix = self.get_node("batch_directive", root=self.batch_system_node).text
        directive_prefix = "" if directive_prefix is None else directive_prefix

        roots = [self.machine_node, self.batch_system_node]
        for root in roots:
            if root is not None:
                nodes = self.get_nodes("directive", root=root)
                for node in nodes:
                    directive = self.get_resolved_value("" if node.text is None else node.text)
                    if batch_maker is None:
                        default = node.get("default")
                        if default is not None:
                            directive_re = re.compile(r"{{ (\w+) }}", flags=re.M)
                            m = directive_re.search(directive)
                            if m is not None:
                                whole_match = m.group()
                                directive = directive.replace(whole_match, default)
                    else:
                        directive = batch_maker.transform_vars(directive, default=node.get("default"))

                    result.append("%s %s" % (directive_prefix, directive))

        return result

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

    def get_submit_args(self):
        '''
        return a list of touples (flag, name)
        '''
        arg_nodes = self.get_nodes("arg", root=self.batch_system_node)
        if self.machine_node is not None:
            arg_nodes.extend(self.get_nodes("arg", root=self.machine_node))

        values = []
        for arg in arg_nodes:
            values.append((arg.get("flag"), arg.get("name")))

        return values
