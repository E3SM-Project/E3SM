"""
Interface to the config_component.xml files.  This class inherits from EntryID.py
"""
from CIME.XML.standard_module_setup import *

from CIME.XML.entry_id import EntryID
from CIME.XML.files import Files
from CIME.utils import get_cime_root

logger = logging.getLogger(__name__)

class Component(EntryID):

    def __init__(self, infile):
        """
        initialize an object
        """
        files = Files()
        schema = None
        # not checking schema on external components yet
        cimeroot = get_cime_root()
        if  cimeroot in os.path.abspath(infile):
            schema = files.get_schema("CONFIG_CPL_FILE")

        EntryID.__init__(self, infile, schema=schema)

    def get_value(self, name, attribute=None, resolved=False, subgroup=None):
        expect(subgroup is None, "This class does not support subgroups")
        return EntryID.get_value(self, name, attribute, resolved)

    def get_valid_model_components(self):
        """
        return a list of all possible valid generic (e.g. atm, clm, ...) model components
        from the entries in the model CONFIG_CPL_FILE
        """
        components = []
        comps_node = self.get_node("entry", {"id":"COMP_CLASSES"})
        comps = self.get_default_value(comps_node)
        components = comps.split(',')
        return components

    def _get_value_match(self, node, attributes=None, exact_match=False):
        match_value = None
        match_max = 0
        match_count = 0
        match_values = []
        expect(not exact_match, " exact_match not implemented in this method")
        expect(node is not None," Empty node in _get_value_match")
        values = self.get_optional_node("values", root=node)
        if values is None:
            return
        # use the default_value if present
        val_node = self.get_optional_node("default_value", root=node)
        if val_node is None:
            logger.debug("No default_value for %s"%node.get("id"))
            return val_node
        value = val_node.text
        if value is not None and len(value) > 0 and value != "UNSET":
            match_values.append(value)
        for valnode in self.get_nodes("value", root=node):
            # loop through all the keys in valnode (value nodes) attributes
            for key,value in valnode.attrib.iteritems():
                # determine if key is in attributes dictionary
                match_count = 0
                if attributes is not None and key in attributes:
                    if re.search(value, attributes[key]):
                        logger.debug("Value %s and key %s match with value %s"%(value, key, attributes[key]))
                        match_count += 1
                    else:
                        match_count = 0
                        break
            if match_count > 0:
                # append the current result
                if values.get("modifier") == "additive":
                    match_values.append(valnode.text)

                # replace the current result if it already contains the new value
                # otherwise append the current result
                elif values.get("modifier") == "merge":
                    if valnode.text in match_values:
                        del match_values[:]
                    match_values.append(valnode.text)

                # take the *last* best match
                elif match_count >= match_max:
                    del match_values[:]
                    match_max = match_count
                    match_value = valnode.text

        if len(match_values) > 0:
            match_value = " ".join(match_values)

        return match_value

    def print_values(self):
        """
        print values for help and description in target config_component.xml file
        """
        rootnode = self.get_node("help")
        helptext = rootnode.text

        rootnode = self.get_node("description")
        compsets = {}
        descs = self.get_nodes("desc", root=rootnode)
        for desc in descs:
            attrib = desc.get("compset")
            text = desc.text
            compsets[attrib] = text

        logger.info(" %s" %helptext)
        for v in compsets.iteritems():
            label, definition = v
            logger.info("   %20s : %s" %(label, definition))
