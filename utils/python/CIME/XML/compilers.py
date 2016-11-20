"""
Interface to the config_compilers.xml file.  This class inherits from GenericXML.py
"""

from CIME.XML.standard_module_setup import *
from CIME.XML.generic_xml import GenericXML
from CIME.XML.files import Files
from CIME.BuildTools.macrowriterbase import write_macros_file_v1 #write_macros_file_v2

logger = logging.getLogger(__name__)

class Compilers(GenericXML):

    def __init__(self, compiler=None, machine=None, os_= None, mpilib=None, infile=None, files=None):
        """
        initialize an object
        """
        if infile is None:
            if files is None:
                files = Files()
            infile = files.get_value("COMPILERS_SPEC_FILE")
            schema = files.get_schema("COMPILERS_SPEC_FILE")

        GenericXML.__init__(self, infile, schema)
        self._version = self.get_version()

        self.machine        = machine
        self.os             = os_
        self.mpilib         = mpilib
        self.compiler_nodes = None # Listed from last to first
        self.compiler       = compiler
        #Append the contents of $HOME/.cime/config_compilers.xml if it exists
        #This could cause problems if node matchs are repeated when only one is expected
        infile = os.path.join(os.environ.get("HOME"),".cime","config_compilers.xml")
        if os.path.exists(infile):
            GenericXML.read(self, infile)

        if self.compiler is not None:
            self.set_compiler(compiler)

        if self._version != "1.0":
            # Run an XPath query to extract the list of flag variable names.
            ns = {"xs": "http://www.w3.org/2001/XMLSchema"}
            flag_xpath = ".//xs:group[@name='compilerVars']/xs:choice/xs:element[@type='flagsVar']"
            flag_elems = ET.parse(schema).getroot().findall(flag_xpath, ns)
            self.flag_vars = set(elem.get('name') for elem in flag_elems)



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

    def get_value(self, name, attribute=None, resolved=True, subgroup=None):
        """
        Get Value of fields in the config_compilers.xml file
        """
        expect(self.compiler_nodes is not None, "Compiler object has no compiler defined")
        expect(subgroup is None, "This class does not support subgroups")
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

    def write_macros_file(self, macros_file="Macros", output_format="make"):
        if self._version == "1.0":
            # Parse the xml settings into the $macros hash structure
            # put conditional settings in the _COND_ portion of the hash
            # and handle them seperately
            macros = {"_COND_" : {}}

            # Do worst matches first
            for compiler_node in reversed(self.compiler_nodes):
                _add_to_macros(compiler_node, macros)
            return write_macros_file_v1(macros, self.compiler, self.os,
                                        self.machine, macros_file,
                                        output_format)
        else:
            expect(False, "Need to fix this")
#            return write_macros_file_v2(macros_file, output_format)

def _add_to_macros(node, macros):
    for child in node:
        name = child.tag
        attrib = child.attrib
        value = child.text

        if not attrib:
            if name.startswith("ADD_"):
                basename = name[4:]
                if basename in macros:
                    macros[basename] = "%s %s" % (macros[basename], value)
                elif name in macros:
                    macros[name] = "%s %s" % (macros[name], value)
                else:
                    macros[name] = value
            else:
                macros[name] = value

        else:
            cond_macros = macros["_COND_"]
            for key, value2 in attrib.iteritems():
                if key not in cond_macros:
                    cond_macros[key] = {}
                if value2 not in cond_macros[key]:
                    cond_macros[key][value2] = {}
                cond_macros = cond_macros[key][value2]

            cond_macros[name] = value
