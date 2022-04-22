import os, sys, copy
import xml.etree.ElementTree as ET

_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..","..","cime")
sys.path.append(_CIMEROOT)

from CIME.utils import expect

###############################################################################
def parse_string_as_list (string):
###############################################################################
    """
    Takes a string representation of nested list and creates
    a nested list of stirng. For instance, with
        s = "(a,b,(c,d),e)
        l = parse_string_as_list
    we would have l = ['a', 'b', ['c','d'], 'e']
    """

    sub_open = string.find('(')
    sub_close = string.rfind(')')
    last = len(string)-1

    if string[0]==',' or string[-1]==',':
        print ("malformed: starts or ends with ','")
        raise

    if sub_open*sub_close<0:
        print ("malformed: unmatched open/close parentheses")
        raise

    if sub_open<0:
        return string.split(',')

    if sub_open==0 and sub_close==last:
        return parse_string_as_list(string[1:-1]);

    # Parse stuff before '(' and after ')'
    if sub_open==0:
        pre = []
    else:
        if string[sub_open-1]!=',':
            return []
        pre = string[0:sub_open-1].split(',')

    if sub_close==last:
        post = []
    else:
        if string[sub_close+1]!=',':
            return []
        post = string[sub_close+2:].split(',')

    # Parse middle, then glue together
    middle = parse_string_as_list(string[sub_open+1:sub_close])

    result = pre
    result.append(middle)
    result.extend(post)
    return result

###############################################################################
def find_node (root,name):
###############################################################################
    """
    Finds node with given name inside the root element, with a depth-search
    strategy (i.e., follow children before siblings).
    WARNING: this function does not check for uniqueness. If there are
             multiple matches, the first match is returned.
    """

    if root.tag==name:
        return root

    for elem in root:
        found = find_node(elem,name)
        if found is not None:
            return found

    return None

###############################################################################
def get_child (root,name,remove=False):
###############################################################################
    """
    Get children with given name. If not found, throws an exception.
    Optionally, the child can be removed from the parent.
    """

    expect (len(root.findall(name))==1,
            "There must be exactly one {} entry inside {}".format(name,root.tag))
    child = root.find(name)
    if remove:
        root.remove(child)

    return child

###############################################################################
def has_child (root,name):
###############################################################################
    """
    Check if root element has child with given name
    """

    return False if root.find(name) is None else True

###############################################################################
def refine_type(entry, force_type=None):
###############################################################################
    """
    Try to convert the text entry to the appropriate type based on its contents.
    """
    # We want to preserve strings representing lists
    if (entry[0]=="(" and entry[-1]==")") or \
       (entry[0]=="[" and entry[-1]=="]") :
        return entry

    if "," in entry:
        result = [refine_type(item.strip(), force_type=force_type) for item in entry.split(",") if item.strip() != ""]
        expected_type = type(result[0])
        for item in result[1:]:
            expect(isinstance(item, expected_type), "List '{}' has inconsistent types inside".format(entry))

        return result

    if force_type:
        try:
            if force_type == "logical":
                if entry.upper() == "TRUE":
                    return True
                elif entry.upper() == "FALSE":
                    return False
                else:
                    return bool(int(entry))

            elif force_type == "integer":
                return int(entry)
            elif force_type == "real":
                return float(entry)
            elif force_type == "string":
                return str(entry)
            else:
                expect(False, "Bad force_type '{}'".format(force_type))
                return None

        except ValueError:
            expect(False, "Could not use '{}' as type '{}'".format(entry, force_type))
            return None

    if entry.upper() == "TRUE":
        return True
    elif entry.upper() == "FALSE":
        return False

    try:
        v = int(entry)
        return v
    except ValueError:
        pass

    try:
        v = float(entry)
        return v
    except ValueError:
        return entry

###############################################################################
def derive_type(entry):
###############################################################################
    refined_value = refine_type(entry)
    if isinstance(refined_value, list):
        refined_value = refined_value[0]

    if isinstance(refined_value, bool):
        return "logical"
    elif isinstance(refined_value, int):
        return "integer"
    elif isinstance(refined_value, float):
        return "real"
    elif isinstance(refined_value, str):
        return "string"
    else:
        expect(False, "Couldn't derive type of '{}'".format(entry))
        return None

###############################################################################
def check_value(elem, value):
###############################################################################
    """
    Check that a parameter's value is in the valid list
    """
    v = value
    vtype = None
    if "type" in elem.attrib.keys():
        vtype = elem.attrib["type"]
        v = refine_type(v,force_type=vtype)
        expect (v is not None,
                "Error! Value '{}' for element '{}' does not satisfy the constraint type={}"
                .format(value,elem.tag,vtype))

    if "valid_values" in elem.attrib.keys():
        valids_str = elem.attrib["valid_values"]
        valids = [refine_type(item.strip(), force_type=vtype) for item in valids_str.split(",")]
        expect(v in valids,
                "Invalid value '{}' for element '{}'. Value not in the valid list ('{}')".format(value, elem.tag, valids))

###############################################################################
def check_all_values(root):
###############################################################################
    """
    Check that all values in the xml tree do not violate their metadata
    """

    has_children = len(root)>0
    if has_children:
        for c in root:
            check_all_values(c)
    else:
        if "type" not in root.attrib.keys():
            root.attrib["type"] = derive_type(root.text)
        check_value(root,root.text)

###############################################################################
def resolve_inheritance (root,elem):
###############################################################################
    """
    If elem inherits from another node within $root, this function adds all
    children of its "parent" to elem. If parent also inherits, first
    resolve parent recursively. If parent is not found, throw an exception
    """

    if "inherit" in elem.attrib.keys():
        parent_name = elem.attrib["inherit"]

        parent = find_node(root,parent_name)
        expect (elem is not None,
                "Error! Parent {} of {} not found within root {}"
                .format(parent_name,elem.tag,root.tag))

        # Make sure the parent is fully resolved
        resolve_inheritance(root,parent)

        del elem.attrib["inherit"]

        for entry in parent:
            # Add the parent's default only if this element does not
            # have a more specialized version
            if not has_child(elem,entry.tag):
                new_entry = copy.deepcopy(entry)
                elem.append(new_entry)

    for child in elem:
        resolve_inheritance(root,child)

###############################################################################
def resolve_all_inheritances (root):
###############################################################################
    """
    Resolve all inheritances in the root tree
    """

    for elem in root:
        resolve_inheritance(root,elem)

###############################################################################
def get_valid_selectors(xml_root):
###############################################################################
    """
    Extract the <selector> node from the xml root, verifying
    its integrity, and returning selectors as a dict.
    """

    # Get the right XML element, and iterate over its children
    selectors_elem = get_child(xml_root,"selectors",remove=True)
    selectors = {}
    for selector in selectors_elem:
        expect(selector.tag == "selector",
               "Expected selector tag, not {}".format(selector.tag))

        selector_name  = selector.attrib["name"]
        selector_env   = selector.attrib["case_env"]
        if "regex" in selector.attrib:
            selector_regex = selector.attrib["regex"]
        else:
            selector_regex = "(.*)" # Just grab the whole thing

        selectors[selector_name] = (selector_env, selector_regex)

    return selectors 

###############################################################################
def gen_atm_proc_group (ap_names_list, atm_procs_defaults, name=None):
###############################################################################
    """
    Given a (possibly nested) list of atm procs names, and the defaults
    section for each atm proc, builds an XML node containing the tree
    of atm procs to be included in the input file, along with their parameters.
    """

    valid_blocks = []
    for child in atm_procs_defaults:
        valid_blocks.append(child.tag)

    if name is None:
        # This is a group defined on the fly, via (a,b,c) syntax.
        # Create empty group, make it inherit from default atm proc group type
        ap_group = ET.Element("__APG__")
        ap_group.attrib["inherit"] = "atm_proc_group"
        resolve_inheritance(atm_procs_defaults,ap_group)
    else:
        # This is a named group, so we *expect* it to be in the atm_procs_default tree
        ap_group = copy.deepcopy(get_child(atm_procs_defaults,name))

    # The 'atm_procs_list' child will always be formatted as a nested list,
    # like (a,((b,c,d),e),f). However, the name of the process needs to be
    # without parentheses to be a valid XML tag, so we form the tag of
    # a group as "group.ap1_ap2_..._apN." (with apX possibly itself a group tag)
    names = []
    nested_list = []
    for ap in ap_names_list:
        # There are two ways to declare a group:
        #  - have its string be "(ap1,ap2,...,apXYZ)", with each ap possibly
        #    itself a group string. This group is built on the fly based
        #    on the building blocks specs.
        #  - declare an element block (e.g., "my_awesome_group") in XML defaults,
        #    and make it inherit from the atm_proc_group default section
        #    In this case, the element must have the "atm_procs_list" child,
        #    storing a string "(ap1,ap2,...,apXYZ)", so that we can generate the group
        if isinstance(ap,list):
            # This entry of the group is itself a group, declared via string list syntax.
            # Recurse, and append all entries
            sub_group = gen_atm_proc_group(ap,atm_procs_defaults)
            inner_str = get_child(sub_group,"atm_procs_list")
            nested_list.append(inner_str.text)
            names.append(sub_group.tag)

            ap_group.insert(len(ap_group),sub_group)
        else:
            expect (ap in valid_blocks, "Unrecognized atm proc name: {}".format(ap))

            default = get_child(atm_procs_defaults,ap)

            #  ap_group.append(proc)
            if has_child(default,"Type") and get_child(default,"Type").text=="Group":
                # This entry of the group is itself a group, with pre-defined
                # defaults. Let's add its entries to it
                sub_group_procs_list_str = get_child(default,"atm_procs_list").text
                sub_group_procs_list = parse_string_as_list(sub_group_procs_list_str)
                proc = gen_atm_proc_group (sub_group_procs_list, atm_procs_defaults, name=ap)
            else:
                # Just use the defaults
                proc = copy.deepcopy(default)

            nested_list.append(proc.tag)
            names.append(proc.tag)
            ap_group.insert(len(ap_group),proc)

    # Set the atm proc list in here
    ap_group_procs_list = get_child(ap_group,"atm_procs_list")
    ap_group_procs_list.text = "(" + ",".join(nested_list) + ")"
    ap_group_procs_list.attrib["locked"] = "true"
    if ap_group.tag=="__APG__":
        ap_group.tag = "group." + "_".join(names) + "."

    return ap_group

###############################################################################
def gen_atm_procs(atm_procs_list_str, atm_procs_defaults):
###############################################################################
    """
    Generate the atm proc group representing the full atmosphere, given
    the string representation and the defaults section.
    """

    ap_names_list = parse_string_as_list(atm_procs_list_str)

    return gen_atm_proc_group (ap_names_list, atm_procs_defaults)
