"""
Interface to the config_compilers.xml file.  This class inherits from GenericXML.py
"""

from standard_module_setup import *
from generic_xml import GenericXML
from files import Files
from CIME.utils import expect

logger = logging.getLogger(__name__)

class Compilers(GenericXML):

    def __init__(self, machine=None, os_= None, mpilib=None, infile=None):
        """
        initialize an object
        """
        if (infile is None):
            infile = os.path.join(get_cime_root(), "cime_config", get_model(), "machines", "config_compilers.xml")

        GenericXML.__init__(self, infile)

        self.machine  = machine
        self.os       = os_
        self.mpilib   = mpilib
        self.compiler = None
        self.parent   = None
        self.name     = None

    def get_compiler(self):
        """
        Return the name of the compiler
        """
        return self.name

    def get_node(self, nodename, attributes=None):
        """
        Return data on a node for a compiler
        """
        expect(self.compiler is not None, "Compiler not set, use parent get_node?")
        result = GenericXML.get_optional_node(self, nodename, attributes, root=self.compiler)
        if result is None:
            return GenericXML.get_node(self, nodename, attributes, root=self.parent)

    def get_optional_node(self, nodename, attributes=None):
        """
        Return data on a node for a compiler
        """
        expect(self.compiler is not None, "Compiler not set, use parent get_node?")
        result = GenericXML.get_optional_node(self, nodename, attributes, root=self.compiler)
        if result is None:
            return GenericXML.get_optional_node(self, nodename, attributes, root=self.parent)

    def _is_compatible(self, compiler_node, machine, os_, mpilib):
        for xmlid, value in [ ("MACH", machine), ("OS", os_), ("MPILIB", mpilib) ]:
            if value is not None and xmlid in compiler_node.attrib and value != compiler_node.attrib[xmlid]:
                return False

        return True

    def list_available_compilers(self, machine=None, os_=None, mpilib=None):
        """
        Return a list of compilers defined for a given CIME_MODEL
        """
        machine = machine if machine else self.machine
        os_     = os_ if os_ else self.os
        mpilib  = mpilib if mpilib else self.mpilib

        compilers = set()
        nodes  = self.get_nodes("compiler")
        for node in nodes:
            if self._is_compatible(node, machine, os_, mpilib):
                compiler = node.attrib["COMPILER"]
                compilers.add(compiler)

        return list(compilers)

    def set_compiler(self, compiler, machine=None, os_=None, mpilib=None):
        """
        Sets the compiler block in the Compilers object

        >>> machobj = Compilers(machine="melvin")
        >>> machobj.set_compiler("gnu")
        >>> machobj.get_compiler()
        'gnu'
        >>> machobj.set_compiler("trump")
        Traceback (most recent call last):
        ...
        SystemExit: ERROR: No compiler trump found
        """
        machine = machine if machine else self.machine
        os_     = os_ if os_ else self.os
        mpilib  = mpilib if mpilib else self.mpilib

        if self.name != compiler or self.machine != machine or self.os != os_ or self.mpilib != mpilib:
            nodes = GenericXML.get_nodes(self, "compiler", {"COMPILER" : compiler})
            most_general_match = (None, 99)
            most_specific_match = (None, 0)
            for node in nodes:
                if self._is_compatible(node, machine, os_, mpilib):
                    num_attrib = len(node.attrib)

                    expect(num_attrib != most_general_match[1], "Ambiguous compiler match")
                    expect(num_attrib != most_specific_match[1], "Ambiguous compiler match")

                    if num_attrib < most_general_match[1]:
                        most_general_match = (node, num_attrib)

                    if num_attrib > most_specific_match[1]:
                        most_specific_match = (node, num_attrib)


            expect(most_specific_match[0] is not None, "No compiler %s found" % compiler)

            self.compiler = most_specific_match[0]
            self.parent   = most_general_match[0]

            self.name    = compiler
            self.machine = machine
            self.os      = os_
            self.mpilib  = mpilib

        return compiler

    def get_value(self, name, resolved=True):
        """
        Get Value of fields in the config_compilers.xml file
        """
        expect(self.compiler is not None, "Compiler object has no compiler defined")
        value = None

        node = self.get_optional_node(name)
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
