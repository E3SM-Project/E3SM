"""
Interface to the config_compilers.xml file.  This class inherits from GenericXML.py
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.utils import expect, get_cime_root, get_model

logger = logging.getLogger(__name__)

class Compilers(GenericXML):

    def __init__(self, compiler=None, machine=None, os_= None, mpilib=None, infile=None):
        """
        initialize an object
        """
        if (infile is None):
            infile = os.path.join(get_cime_root(), "cime_config", get_model(), "machines", "config_compilers.xml")

        GenericXML.__init__(self, infile)

        self.machine       = machine
        self.os            = os_
        self.mpilib        = mpilib
        self.compiler_node = None
        self.parent_node   = None
        self.compiler      = compiler

        if self.compiler is not None:
            self.set_compiler(compiler)

    def get_compiler(self):
        """
        Return the name of the compiler
        """
        return self.compiler

    def get_compiler_node(self, nodename, attributes=None):
        """
        Return data on a node for a compiler
        """
        expect(self.compiler_node is not None, "Compiler not set, use parent get_node?")
        result = self.get_optional_node(nodename, attributes, root=self.compiler_node)
        if result is None:
            return self.get_node(nodename, attributes, root=self.parent_node)

    def get_optional_compiler_node(self, nodename, attributes=None):
        """
        Return data on a node for a compiler
        """
        expect(self.compiler_node is not None, "Compiler not set, use parent get_node?")
        result = self.get_optional_node(nodename, attributes, root=self.compiler_node)
        if result is None:
            return self.get_optional_node(nodename, attributes, root=self.parent_node)

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

        if self.compiler != compiler or self.machine != machine or self.os != os_ or self.mpilib != mpilib or self.compiler_node is None:
            nodes = self.get_nodes("compiler", {"COMPILER" : compiler})
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

            self.compiler_node = most_specific_match[0]
            self.parent_node   = most_general_match[0]

            self.compiler = compiler
            self.machine  = machine
            self.os       = os_
            self.mpilib   = mpilib

    def get_value(self, name, resolved=True):
        """
        Get Value of fields in the config_compilers.xml file
        """
        expect(self.compiler_node is not None, "Compiler object has no compiler defined")
        value = None

        node = self.get_optional_compiler_node(name)
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
