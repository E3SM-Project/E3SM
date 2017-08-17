"""
Common interface to XML files which follow the entry id format,
this is an abstract class and is expected to
be used by other XML interface modules and not directly.
"""
from CIME.XML.standard_module_setup import *
from CIME.utils import expect, convert_to_string, convert_to_type
from CIME.XML.generic_xml import GenericXML

from copy import deepcopy

logger = logging.getLogger(__name__)

class EntryID(GenericXML):

    def __init__(self, infile=None, schema=None):
        GenericXML.__init__(self, infile, schema)
        self.groups={}

    def get_default_value(self, node, attributes=None):
        """
        Set the value of an entry to the default value for that entry
        """
        value = self._get_value_match(node, attributes)
        if value is None:
            # Fall back to default value
            value = self.get_element_text("default_value", root=node)
        else:
            logger.debug("node is {} value is {}".format(node.get("id"), value))

        if value is None:
            logger.debug("For vid {} value is none".format(node.get("id")))
            value = ""

        return value

    def set_default_value(self, vid, val):
        node = self.get_optional_node("entry", {"id":vid})
        if node is not None:
            val = self.set_element_text("default_value", val, root=node)
            if val is None:
                logger.warn("Called set_default_value on a node without default_value field")

        return val

    def get_value_match(self, vid, attributes=None, exact_match=False, entry_node=None):
        # Handle this case:
        # <entry id ...>
        #  <values>
        #   <value A="a1">X</value>
        #   <value A="a2">Y</value>
        #   <value A="a3" B="b1">Z</value>
        #  </values>
        # </entry>

        if entry_node is not None:
            value = self._get_value_match(entry_node, attributes, exact_match)
        else:
            node = self.get_optional_node("entry", {"id":vid})
            value = None
            if node is not None:
                value = self._get_value_match(node, attributes, exact_match)
        logger.debug("(get_value_match) vid {} value {}".format(vid, value))
        return value

    def _get_value_match(self, node, attributes=None, exact_match=False):
        '''
        Note that the component class has a specific version of this function
        '''
        # if there is a <values> element - check to see if there is a match attribute
        # if there is NOT a match attribute, then set the default to "first"
        # this is different than the component class _get_value_match where the default is "last"
        values_node = self.get_optional_node("values", root=node)
        if values_node is not None:
            match_type = values_node.get("match", default="first")
        else:
            match_type = "first"

        # Store nodes that match the attributes and their scores.
        matches = []
        nodes = self.get_nodes("value", root=node)
        for vnode in nodes:
            # For each node in the list start a score.
            score = 0
            if attributes:
                for attribute in vnode.keys():
                    # For each attribute, add to the score.
                    score += 1
                    # If some attribute is specified that we don't know about,
                    # or the values don't match, it's not a match we want.
                    if exact_match:
                        if attribute not in attributes or \
                                attributes[attribute] != vnode.get(attribute):
                            score = -1
                            break
                    else:
                        if attribute not in attributes or not \
                                re.search(vnode.get(attribute),attributes[attribute]):
                            score = -1
                            break

            # Add valid matches to the list.
            if score >= 0:
                matches.append((score, vnode))

        if not matches:
            return None

        # Get maximum score using either a "last" or "first" match in case of a tie
        max_score = -1
        mnode = None
        for score,node in matches:
            if match_type == "last":
                # take the *last* best match
                if score >= max_score:
                    max_score = score
                    mnode = node
            elif match_type == "first":
                # take the *first* best match
                if score > max_score:
                    max_score = score
                    mnode = node
            else:
                expect(False, 
                       "match attribute can only have a value of 'last' or 'first', value is %s" %match_type)

        return mnode.text

    def get_node_element_info(self, vid, element_name):
        node = self.get_optional_node("entry", {"id":vid})
        if node is None:
            return None
        else:
            return self._get_node_element_info(node, element_name)

    def _get_node_element_info(self, node, element_name):
        return self.get_element_text(element_name, root=node)

    def _get_type_info(self, node):
        if node is None:
            return None
        val = self._get_node_element_info(node, "type")
        if val is None:
            return "char"
        return val

    def get_type_info(self, vid):
        vid, _, _ = self.check_if_comp_var(vid, None)
        node = self.get_optional_node("entry", {"id":vid})
        return self._get_type_info(node)

    # pylint: disable=unused-argument
    def check_if_comp_var(self, vid, attribute=None):
        # handled in classes
        return vid, None, False

    def _get_default(self, node):
        return self._get_node_element_info(node, "default_value")

    # Get description , expect child with tag "description" for parent node
    def get_description (self, node):
        return self._get_node_element_info(node, "desc")

    # Get group , expect node with tag "group"
    # entry id nodes are children of group nodes
    def get_groups(self, node):
        groups = self.get_nodes("group")
        result = []
        nodes = []
        vid = node.get("id")
        for group in groups:
            nodes = self.get_nodes("entry", attributes={"id":vid}, root=group)
            if nodes:
                result.append(group.get("id"))

        return result

    def get_valid_values(self, vid):
        node = self.get_optional_node("entry", {"id":vid})
        if node is None:
            return None
        return self._get_valid_values(node)

    def _get_valid_values(self, node):
        valid_values = self.get_element_text("valid_values", root=node)
        valid_values_list = None
        if valid_values is not None:
            valid_values_list = [item.lstrip() for item in valid_values.split(',')]
        return valid_values_list

    def set_valid_values(self, vid, new_valid_values):
        node = self.get_optional_node("entry", {"id":vid})
        if node is None:
            return None
        return self._set_valid_values(node, new_valid_values)

    def get_nodes_by_id(self, vid):
        return self.get_nodes("entry", {"id":vid})

    def _set_valid_values(self, node, new_valid_values):
        old_vv = self._get_valid_values(node)
        if old_vv is None:
            vv_node = ET.Element("valid_values")
            vv_node.text = new_valid_values
            self.add_child(vv_node)
            logger.debug("Adding valid_values {} for {}".format(new_valid_values, node.get("id")))
        else:
            vv_text = self.set_element_text("valid_values", new_valid_values, root=node)
            logger.debug("Replacing valid_values {} with {} for {}".format(old_vv, vv_text, node.get("id")))

        current_value = node.get("value")
        valid_values_list = self._get_valid_values(node)
        if current_value is not None and current_value not in valid_values_list:
            logger.warn("WARNING: Current setting for {} not in new valid values. Updating setting to \"{}\"".format(node.get("id"), valid_values_list[0]))
            self._set_value(node, valid_values_list[0])
        return new_valid_values

    def _set_value(self, node, value, vid=None, subgroup=None, ignore_type=False):
        """
        Set the value of an entry-id field to value
        Returns the value or None if not found
        subgroup is ignored in the general routine and applied in specific methods
        """
        expect(subgroup is None, "Subgroup not supported")
        str_value = self.get_valid_value_string(node, value, vid, ignore_type)
        node.set("value", str_value)
        return value

    def get_valid_value_string(self, node, value,vid=None,  ignore_type=False):
        valid_values = self._get_valid_values(node)
        if ignore_type:
            expect(type(value) is str, "Value must be type string if ignore_type is true")
            str_value = value
            return str_value
        type_str = self._get_type_info(node)
        str_value = convert_to_string(value, type_str, vid)

        if valid_values is not None and not str_value.startswith('$'):
            expect(str_value in valid_values, "Did not find {} in valid values for {}: {}".format(value, vid, valid_values))
        return str_value

    def set_value(self, vid, value, subgroup=None, ignore_type=False):
        """
        Set the value of an entry-id field to value
        Returns the value or None if not found
        subgroup is ignored in the general routine and applied in specific methods
        """
        val = None
        root = self.root if subgroup is None else self.get_optional_node("group", {"id":subgroup})
        node = self.get_optional_node("entry", {"id":vid}, root=root)
        if node is not None:
            val = self._set_value(node, value, vid, subgroup, ignore_type)
        return val

    def get_values(self, vid, attribute=None, resolved=True, subgroup=None):
        """
        Same functionality as get_value but it returns a list, if the
        value in xml contains commas the list have multiple elements split on
        commas
        """
        results = []
        node = self.get_optional_node("entry", {"id":vid})
        if node is None:
            return results
        str_result = self._get_value(node, attribute=attribute, resolved=resolved, subgroup=subgroup)
        str_results = str_result.split(',')
        for result in str_results:
            # Return value as right type if we were able to fully resolve
            # otherwise, we have to leave as string.
            if "$" in result:
                results.append(result)
            else:
                type_str = self._get_type_info(node)
                results.append( convert_to_type(result, type_str, vid))
        return results
    #pylint: disable=arguments-differ
    def get_value(self, vid, attribute=None, resolved=True, subgroup=None):
        """
        Get a value for entry with id attribute vid.
        or from the values field if the attribute argument is provided
        and matches
        """
        root = self.root if subgroup is None else self.get_optional_node("group", {"id":subgroup})
        node = self.get_optional_node("entry", {"id":vid}, root=root)
        if node is None:
            return
        val = self._get_value(node, attribute=attribute, resolved=resolved, subgroup=subgroup)
        # Return value as right type if we were able to fully resolve
        # otherwise, we have to leave as string.
        if val is None:
            return val
        elif "$" in val:
            return val
        else:
            type_str = self._get_type_info(node)
            return convert_to_type(val, type_str, vid)

    def _get_value(self, node, attribute=None, resolved=True, subgroup=None):
        """
        internal get_value, does not convert to type
        """
        logger.debug("(_get_value) ({}, {}, {})".format(attribute, resolved, subgroup))
        val = None
        if node is None:
            logger.debug("No node")
            return val

        logger.debug("Found node {} with attributes {}".format(node.tag , node.attrib))
        if attribute:
            val = self.get_element_text("value", attributes=attribute, root=node)
        elif node.get("value") is not None:
            val = node.get("value")
        else:
            val = self.get_default_value(node)

        if resolved:
            val = self.get_resolved_value(val)

        return val

    def get_child_content(self, vid, childname):
        val = None
        node = self.get_optional_node("entry", {"id" : vid})
        if node is not None:
            val = self.get_element_text(childname, root=node)
        return val

    def get_elements_from_child_content(self, childname, childcontent):
        nodes = self.get_nodes("entry")
        elements = []
        for node in nodes:
            content = self.get_element_text(childname, root=node)
            expect(content is not None,"No childname {} for id {}".format(childname, node.get("id")))
            if content == childcontent:
                elements.append(deepcopy(node))

        return elements

    def add_elements_by_group(self, srcobj, attributes=None, infile=None):
        """
        Add elements from srcobj to self under the appropriate
        group element, entries to be added must have a child element
        <file> with value "infile"
        """
        if infile is None:
            infile = os.path.basename(self.filename)
        # First get the list of entries in srcobj with matching file children

        nodelist = srcobj.get_elements_from_child_content('file', infile)

        # For matchs found: Remove {<group>, <file>, <values>}
        # children from each entry and set the default value for the
        # new entries in self - putting the entries as children of
        # group elements in file $file
        for src_node in nodelist:
            node  = deepcopy(src_node)
            gnode = src_node.find(".//group")
            gname = gnode.text
            if gname is None:
                gname = "group_not_set"
            # If group with id=$gname does not exist in self.groups
            # then create the group node and add it to infile file
            if gname not in self.groups.keys():
                newgroup = ET.Element("group")
                newgroup.set("id",gname)
                # initialize an empty list
                self.groups[gname] = newgroup
                self.add_child(newgroup)

            # Remove {<group>, <file>, <values>} from the entry element
            node = self.cleanupnode(node)

            # Add the entry element to the group
            self.groups[gname].append(node)

            # Set the default value, it may be determined by a regular
            # expression match to a dictionary value in attributes matching a
            # value attribute in node
            value = srcobj.get_default_value(src_node, attributes)
            if value is not None and len(value):
                self._set_value(node,value)
            logger.debug ("Adding to group " + gname)

        return nodelist

    def cleanupnode(self, node):
        """
        in env_base.py, not expected to get here
        """
        expect(False, " Not expected to be here {}".format(node.get("id")))

    def compare_xml(self, other):
        xmldiffs = {}
        f1nodes = self.get_nodes("entry")
        for node in f1nodes:
            vid = node.get("id")
            f2match = other.get_optional_node("entry", attributes={"id":vid})
            if f2match is not None:
                f1val = self.get_value(vid, resolved=False)
                if f1val is not None:
                    f2val = other.get_value(vid, resolved=False)
                    if f1val != f2val:
                        xmldiffs[vid] = [f1val, f2val]
                else:
                    f1val = ET.tostring(node, method="text")
                    f2val = ET.tostring(f2match, method="text")
                    if f2val != f1val:
                        f1value_nodes = self.get_nodes("value", root=node)
                        for valnode in f1value_nodes:
                            f2valnodes = other.get_nodes("value", root=f2match, attributes=valnode.attrib)
                            for f2valnode in f2valnodes:
                                if valnode.attrib is None and f2valnode.attrib is None or \
                                   f2valnode.attrib == valnode.attrib:
                                    if other.get_resolved_value(f2valnode.text) != self.get_resolved_value(valnode.text):
                                        xmldiffs["{}:{}".format(vid, valnode.attrib)] = [valnode.text, f2valnode.text]

        return xmldiffs

    def __iter__(self):
        for node in self.get_nodes("entry"):
            vid = node.get("id")
            yield vid, self.get_value(vid)
