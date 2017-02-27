"""
Base class for env files .  This class inherits from EntryID.py
"""
import string
from CIME.XML.standard_module_setup import *
from CIME.XML.entry_id import EntryID
from CIME.XML.headers import Headers
from CIME.utils import convert_to_type, convert_to_string
logger = logging.getLogger(__name__)

class EnvBase(EntryID):

    def __init__(self, case_root, infile, schema=None):
        if case_root is None:
            case_root = os.getcwd()

        if os.path.isabs(infile):
            fullpath = infile
        else:
            fullpath = os.path.join(case_root, infile)

        EntryID.__init__(self, fullpath, schema=schema)
        if not os.path.isfile(fullpath):
            headerobj = Headers()
            headernode = headerobj.get_header_node(os.path.basename(fullpath))
            self.root.append(headernode)

    def set_components(self, components):
        if hasattr(self, '_components'):
            # pylint: disable=attribute-defined-outside-init
            self._components = components

    def check_if_comp_var(self, vid, attribute=None):
        # pylint: disable=no-member
        if not hasattr(self, "_component_value_list") or\
                (self.get_nodes("entry", {"id" : vid}) and \
                     not vid in self._component_value_list):
            return vid, None, False

        comp = None
        if vid in self._component_value_list:
            if attribute is not None:
                if "component" in attribute:
                    comp = attribute["component"]
            return vid, comp, True

        for comp in self._components:
            if "_"+comp in vid:
                vid = string.replace(vid, '_'+comp, '', 1)
                break
            elif comp+"_" in vid:
                vid = string.replace(vid, comp+'_', '', 1)
                break
        if vid in self._component_value_list:
            return vid, comp, True
        return vid, None, False

    def get_value(self, vid, attribute=None, resolved=True, subgroup=None):
        """
        Get a value for entry with id attribute vid.
        or from the values field if the attribute argument is provided
        and matches
        """
        value = None
        vid, comp, iscompvar = self.check_if_comp_var(vid, attribute)
        logger.debug( "vid %s comp %s iscompvar %s"%(vid, comp, iscompvar))
        if iscompvar:
            if comp is None:
                if subgroup is not None:
                    comp = subgroup
                else:
                    logger.debug("Not enough info to get value for %s"%vid)
                    return value
            if attribute is None:
                attribute = {"component" : comp}
            else:
                attribute["component"] = comp
            node = self.get_optional_node("entry", {"id":vid})
            if node is not None:
                type_str = self._get_type_info(node)
                val = self.get_element_text("value", attribute, root=node)
                if val is not None:
                    if val.startswith("$"):
                        value = val
                    else:
                        value = convert_to_type(val,type_str, vid)
                return value
        return EntryID.get_value(self, vid, attribute=attribute, resolved=resolved, subgroup=subgroup)

    def set_value(self, vid, value, subgroup=None, ignore_type=False):
        """
        Set the value of an entry-id field to value
        Returns the value or None if not found
        subgroup is ignored in the general routine and applied in specific methods
        """
        vid, comp, iscompvar = self.check_if_comp_var(vid, None)
        val = None
        node = self.get_optional_node("entry", {"id":vid})
        if node is not None:
            if iscompvar and comp is None:
                # pylint: disable=no-member
                for comp in self._components:
                    val = self._set_value(node, value, vid, subgroup, ignore_type, component=comp)
            else:
                val = self._set_value(node, value, vid, subgroup, ignore_type, component=comp)
        return val
    # pylint: disable=arguments-differ
    def _set_value(self, node, value, vid=None, subgroup=None, ignore_type=False, component=None):
        if vid is None:
            vid = node.get("id")
        vid, _, iscompvar = self.check_if_comp_var(vid, None)

        if iscompvar:
            attribute = {"component":component}
            type_str = self._get_type_info(node)
            val = self.set_element_text("value", convert_to_string(value, type_str, vid), attribute, root=node)
            return val
        val = EntryID._set_value(self, node, value, vid, subgroup, ignore_type)
        return val

    def get_nodes_by_id(self, varid):
        varid, _, _ = self.check_if_comp_var(varid, None)
        return EntryID.get_nodes_by_id(self, varid)

    def cleanupnode(self, node):
        """
        Remove the <group>, <file>, <values> and <value> childnodes from node
        """
        fnode = node.find(".//file")
        node.remove(fnode)
        gnode = node.find(".//group")
        node.remove(gnode)
        dnode = node.find(".//default_value")
        if dnode is not None:
            node.remove(dnode)
        vnode = node.find(".//values")
        vid = node.get("id")
        _, _, iscompvar = self.check_if_comp_var(vid)
        if vnode is not None and not iscompvar:
            node.remove(vnode)
        return node
