"""
Interface to the config_files.xml file.  This class inherits from EntryID.py
"""
import re
from CIME.XML.standard_module_setup import *

from CIME.XML.entry_id import EntryID
from CIME.utils import expect, get_cime_root, get_model

logger = logging.getLogger(__name__)

class Files(EntryID):

    def __init__(self, comp_interface="mct"):
        """
        initialize an object

        >>> files = Files()
        >>> files.get_value('CASEFILE_HEADERS',resolved=False)
        '$CIMEROOT/config/config_headers.xml'
        """
        cimeroot = get_cime_root()
        infile = os.path.join(cimeroot, "config", get_model(), "config_files.xml")
        expect(os.path.isfile(infile), "Could not find or open file {}".format(infile))
        schema = os.path.join(cimeroot, "config", "xml_schemas", "entry_id.xsd")
        EntryID.__init__(self, infile, schema=schema)
        config_files_override = os.path.join(os.path.dirname(cimeroot),".config_files.xml")
        # variables COMP_ROOT_DIR_{} are mutable, all other variables are read only
        self.COMP_ROOT_DIR = {}
        self._comp_interface = comp_interface
        self._cpl_comp = {}
        # .config_file.xml at the top level may overwrite COMP_ROOT_DIR_ nodes in config_files

        if os.path.isfile(config_files_override):
            self.read(config_files_override)
            self.overwrite_existing_entries()

    def get_value(self, vid, attribute=None, resolved=True, subgroup=None):
        if vid == "COMP_ROOT_DIR_CPL":
            if self._cpl_comp:
                attribute = self._cpl_comp
            elif attribute:
                self._cpl_comp = attribute
            else:
                self._cpl_comp['component'] = 'cpl'
        if "COMP_ROOT_DIR" in vid:
            if vid in self.COMP_ROOT_DIR:
                if attribute is not None:
                    if vid+attribute["component"] in self.COMP_ROOT_DIR:
                        return self.COMP_ROOT_DIR[vid+attribute["component"]]
                else:
                    return self.COMP_ROOT_DIR[vid]

        newatt = {"comp_interface":self._comp_interface}
        if attribute:
            newatt.update(attribute)
        value = super(Files, self).get_value(vid, attribute=newatt, resolved=False, subgroup=subgroup)
        if value is None and attribute is not None:
            value = super(Files, self).get_value(vid, attribute=attribute, resolved=False, subgroup=subgroup)
        if value is None:
            value = super(Files, self).get_value(vid, attribute=None, resolved=False, subgroup=subgroup)

        if "COMP_ROOT_DIR" not in vid and value is not None and "COMP_ROOT_DIR" in value:
            m = re.search("(COMP_ROOT_DIR_[^/]+)/", value)
            comp_root_dir_var_name = m.group(1)
            newatt = {"comp_interface":self._comp_interface}
            if attribute:
                newatt.update(attribute)

            crd_node = self.scan_optional_child(comp_root_dir_var_name, attributes=newatt)
            if crd_node:
                comp_root_dir = self.get_value(comp_root_dir_var_name, attribute=newatt, resolved=False, subgroup=subgroup)
            else:
                comp_root_dir = self.get_value(comp_root_dir_var_name, attribute=attribute, resolved=False, subgroup=subgroup)
            self.set_value(comp_root_dir_var_name, comp_root_dir,subgroup=attribute)
            if resolved:
                value = value.replace("$"+comp_root_dir_var_name, comp_root_dir)

        if resolved and value is not None:
            value = value.replace("$COMP_INTERFACE", self._comp_interface)
            value = self.get_resolved_value(value)
        return value

    def set_value(self, vid, value,subgroup=None,ignore_type=False):
        if "COMP_ROOT_DIR" in vid:
            if subgroup is not None:
                self.COMP_ROOT_DIR[vid+subgroup["component"]] = value
            else:
                self.COMP_ROOT_DIR[vid] = value

        else:
            expect(False, "Attempt to set a nonmutable variable {}".format(vid))
        return value


    def get_schema(self, nodename, attributes=None):
        node = self.get_optional_child("entry", {"id":nodename})
        schemanode = self.get_optional_child("schema", root=node, attributes=attributes)
        if schemanode is not None:
            logger.debug("Found schema for {}".format(nodename))
            return self.get_resolved_value(self.text(schemanode))
        return None

    def get_components(self, nodename):
        node = self.get_optional_child("entry", {"id":nodename})
        if node is not None:
            valnodes = self.get_children("value", root=self.get_child("values", root=node))
            values = []
            for valnode in valnodes:
                value = self.get(valnode, "component")
                values.append(value)
            return values

        return None
