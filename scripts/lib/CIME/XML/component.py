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

    #pylint: disable=arguments-differ
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
        """
        return the best match for the node <values> entries
        Note that a component object uses a different matching algorithm than an entryid object
        For a component object the _get_value_match used is below  and is not the one in entry_id.py
        """
        match_value = None
        match_max = 0
        match_count = 0
        match_values = []
        expect(not exact_match, " exact_match not implemented in this method")
        expect(node is not None," Empty node in _get_value_match")
        values = self.get_optional_node("values", root=node)
        if values is None:
            return

        # determine match_type if there is a tie 
        # ASSUME a default of "last" if "match" attribute is not there
        match_type = values.get("match", default="last")

        # use the default_value if present
        val_node = self.get_optional_node("default_value", root=node)
        if val_node is None:
            logger.debug("No default_value for {}".format(node.get("id")))
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
                        logger.debug("Value {} and key {} match with value {}".format(value, key, attributes[key]))
                        match_count += 1
                    else:
                        match_count = 0
                        break

            # a match is found
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

                else: 
                    if match_type == "last":
                        # take the *last* best match
                        if match_count >= match_max:
                            del match_values[:]
                            match_max = match_count
                            match_value = valnode.text
                    elif match_type == "first":
                        # take the *first* best match
                        if match_count > match_max:
                            del match_values[:]
                            match_max = match_count
                            match_value = valnode.text
                    else:
                        expect(False, "match attribute can only have a value of 'last' or 'first'")

        if len(match_values) > 0:
            match_value = " ".join(match_values)

        return match_value

    #pylint: disable=arguments-differ
    def get_description(self, compsetname):
        rootnode = self.get_node("description")
        desc_nodes = self.get_nodes("desc", root=rootnode)
        desc = ""
        for node in desc_nodes:
            compsetmatch = node.get("compset")
            if re.search(compsetmatch, compsetname):
                desc += node.text

        return desc

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

        logger.info(" {}".format(helptext))
        for v in sorted(compsets.iteritems()):
            label, definition = v
            logger.info("   {:20s} : {}".format(label, definition))
