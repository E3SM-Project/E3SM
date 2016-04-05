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

        self.machine        = machine
        self.os             = os_
        self.mpilib         = mpilib
        self.compiler_nodes = None # Listed from last to first
        self.compiler       = compiler

        if self.compiler is not None:
            self.set_compiler(compiler)

    def get_compiler(self):
        """
        Return the name of the compiler
        """
        return self.compiler

    def get_optional_compiler_node(self, nodename, attributes=None):
        """
        Return data on a node for a compiler
        """
        expect(self.compiler_nodes is not None, "Compiler not set, use parent get_node?")
        for compiler_node in self.compiler_nodes:
            result = self.get_optional_node(nodename, attributes, root=compiler_node)
            if result is not None:
                return result

        return None

    def _is_compatible(self, compiler_node, compiler, machine, os_, mpilib):
        for xmlid, value in [ ("COMPILER", compiler), ("MACH", machine), ("OS", os_), ("MPILIB", mpilib) ]:
            if value is not None and xmlid in compiler_node.attrib and value != compiler_node.get(xmlid):
                return False

        return True

    def set_compiler(self, compiler, machine=None, os_=None, mpilib=None):
        """
        Sets the compiler block in the Compilers object

        >>> machobj = Compilers(machine="melvin")
        >>> machobj.set_compiler("gnu")
        >>> machobj.get_compiler()
        'gnu'
        """
        machine = machine if machine else self.machine
        os_     = os_ if os_ else self.os
        mpilib  = mpilib if mpilib else self.mpilib

        if self.compiler != compiler or self.machine != machine or self.os != os_ or self.mpilib != mpilib or self.compiler_nodes is None:
            self.compiler_nodes = []
            nodes = self.get_nodes("compiler")
            for node in nodes:
                if self._is_compatible(node, compiler, machine, os_, mpilib):
                    self.compiler_nodes.append(node)

            self.compiler_nodes.reverse()

            self.compiler = compiler
            self.machine  = machine
            self.os       = os_
            self.mpilib   = mpilib

    def get_value(self, name, resolved=True):
        """
        Get Value of fields in the config_compilers.xml file
        """
        expect(self.compiler_nodes is not None, "Compiler object has no compiler defined")
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
