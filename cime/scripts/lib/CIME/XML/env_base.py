"""
Base class for env files.  This class inherits from EntryID.py
"""
from CIME.XML.standard_module_setup import *
from CIME.XML.entry_id import EntryID
from CIME.XML.headers import Headers
from CIME.utils import convert_to_type
logger = logging.getLogger(__name__)

class EnvBase(EntryID):

    def __init__(self, case_root, infile, schema=None, read_only=False):
        if case_root is None:
            case_root = os.getcwd()

        if os.path.isabs(infile):
            fullpath = infile
        else:
            fullpath = os.path.join(case_root, infile)

        EntryID.__init__(self, fullpath, schema=schema, read_only=read_only)

        self._id_map = None
        self._group_map = None

        if not os.path.isfile(fullpath):
            headerobj = Headers()
            headernode = headerobj.get_header_node(os.path.basename(fullpath))
            self.add_child(headernode)
        else:
            self._setup_cache()

    def _setup_cache(self):
        self._id_map = {}    # map id directly to nodes
        self._group_map = {} # map group name to entry id dict

        group_elems = self.get_children("group")
        for group_elem in group_elems:
            group_name = self.get(group_elem, "id")
            expect(group_name not in self._group_map, "Repeat group '{}'".format(group_name))
            group_map = {}
            self._group_map[group_name] = group_map
            entry_elems = self.get_children("entry", root=group_elem)
            for entry_elem in entry_elems:
                entry_id = self.get(entry_elem, "id")
                expect(entry_id not in group_map, "Repeat entry '{}' in group '{}'".format(entry_id, group_name))
                group_map[entry_id] = entry_elem
                if entry_id in self._id_map:
                    self._id_map[entry_id].append(entry_elem)
                else:
                    self._id_map[entry_id] = [entry_elem]

        self.lock()

    def change_file(self, newfile, copy=False):
        self.unlock()
        EntryID.change_file(self, newfile, copy=copy)
        self._setup_cache()

    def get_children(self, name=None, attributes=None, root=None):
        if self.locked and name == "entry" and attributes is not None and attributes.keys() == ["id"]:
            entry_id = attributes["id"]
            if root is None or self.name(root) == "file":
                if entry_id in self._id_map:
                    return self._id_map[entry_id]
                else:
                    return []
            else:
                expect(self.name(root) == "group", "Unexpected elem '{}' for {}, attrs {}".format(self.name(root), self.filename, self.attrib(root)))
                group_id = self.get(root, "id")
                if group_id in self._group_map and entry_id in self._group_map[group_id]:
                    return [self._group_map[group_id][entry_id]]
                else:
                    return []

        else:
            # Non-compliant look up
            return EntryID.get_children(self, name=name, attributes=attributes, root=root)

    def scan_children(self, nodename, attributes=None, root=None):
        if self.locked and nodename == "entry" and attributes is not None and attributes.keys() == ["id"]:
            return EnvBase.get_children(self, name=nodename, attributes=attributes, root=root)
        else:
            return EntryID.scan_children(self, nodename, attributes=attributes, root=root)

    def set_components(self, components):
        if hasattr(self, '_components'):
            # pylint: disable=attribute-defined-outside-init
            self._components = components

    def check_if_comp_var(self, vid, attribute=None, node=None):
        comp = None
        if node is None:
            nodes = self.scan_children("entry", {"id" : vid})
            if len(nodes):
                node = nodes[0]

        if node:
            valnodes = self.scan_children("value", attributes={"compclass":None}, root=node)
            if len(valnodes) == 0:
                logger.debug("vid {} is not a compvar".format(vid))
                return vid, None, False
            else:
                logger.debug("vid {} is a compvar".format(vid))
                if attribute is not None:
                    comp = attribute["compclass"]
                return vid, comp, True
        else:
            if hasattr(self, "_components") and self._components:
                new_vid = None
                for comp in self._components:
                    if vid.endswith('_'+comp):
                        new_vid = vid.replace('_'+comp, '', 1)
                    elif vid.startswith(comp+'_'):
                        new_vid = vid.replace(comp+'_', '', 1)
                    elif '_' + comp + '_' in vid:
                        new_vid = vid.replace(comp+'_','', 1)
                    if new_vid is not None:
                        break
                if new_vid is not None:
                    logger.debug("vid {} is a compvar with comp {}".format(vid, comp))
                    return new_vid, comp, True

        return vid, None, False

    def get_value(self, vid, attribute=None, resolved=True, subgroup=None):
        """
        Get a value for entry with id attribute vid.
        or from the values field if the attribute argument is provided
        and matches
        """
        value = None
        vid, comp, iscompvar = self.check_if_comp_var(vid, attribute)
        logger.debug("vid {} comp {} iscompvar {}".format(vid, comp, iscompvar))
        if iscompvar:
            if comp is None:
                if subgroup is not None:
                    comp = subgroup
                else:
                    logger.debug("Not enough info to get value for {}".format(vid))
                    return value
            if attribute is None:
                attribute = {"compclass" : comp}
            else:
                attribute["compclass"] = comp
            node = self.scan_optional_child("entry", {"id":vid})
            if node is not None:
                type_str = self._get_type_info(node)
                values = self.get_optional_child("values", root=node)
                node = values if values is not None else node
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
        root = self.root if subgroup is None else self.get_optional_child("group", {"id":subgroup})
        node = self.scan_optional_child("entry", {"id":vid}, root=root)
        if node is not None:
            if iscompvar and comp is None:
                # pylint: disable=no-member
                for comp in self._components:
                    val = self._set_value(node, value, vid, subgroup, ignore_type, compclass=comp)
            else:
                val = self._set_value(node, value, vid, subgroup, ignore_type, compclass=comp)
        return val

    # pylint: disable=arguments-differ
    def _set_value(self, node, value, vid=None, subgroup=None, ignore_type=False, compclass=None):
        if vid is None:
            vid = self.get(node, "id")
        vid, _, iscompvar = self.check_if_comp_var(vid, node=node)

        if iscompvar:
            expect(compclass is not None, "compclass must be specified if is comp var")
            attribute = {"compclass":compclass}
            str_value = self.get_valid_value_string(node, value, vid, ignore_type)
            values = self.get_optional_child("values", root=node)
            node = values if values is not None else node
            val = self.set_element_text("value", str_value, attribute, root=node)
        else:
            val = EntryID._set_value(self, node, value, vid, subgroup, ignore_type)
        return val

    def get_nodes_by_id(self, varid):
        varid, _, _ = self.check_if_comp_var(varid, None)
        return EntryID.get_nodes_by_id(self, varid)

    def cleanupnode(self, node):
        """
        Remove the <group>, <file>, <values> and <value> childnodes from node
        """
        fnode = self.get_child("file", root=node)
        self.remove_child(fnode, node)
        gnode = self.get_child("group", root=node)
        self.remove_child(gnode, node)
        dnode = self.get_optional_child("default_value", root=node)
        if dnode is not None:
            self.remove_child(dnode, node)

        vnode = self.get_optional_child("values", root=node)
        if vnode is not None:
            componentatt = self.get_children("value", attributes={"component":"ATM"}, root=vnode)
            # backward compatibility (compclasses and component were mixed
            # now we seperated into component and compclass)
            if len(componentatt) > 0:
                for ccnode in self.get_children("value", attributes={"component":None}, root=vnode):
                    val = self.get(ccnode, "component")
                    self.pop(ccnode, "component")
                    self.set(ccnode, "compclass", val)

            compclassatt = self.get_children("value", attributes={"compclass":None}, root=vnode)
            if len(compclassatt) == 0:
                self.remove_child(vnode, root=node)

        return node
