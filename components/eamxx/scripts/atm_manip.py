"""
Retrieve nodes from EAMxx XML config file.
"""

import sys, os, re, pathlib
from collections import namedtuple

# Used for doctests
import xml.etree.ElementTree as ET # pylint: disable=unused-import

# Add path to cime_config folder
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "cime_config"))
from eamxx_buildnml_impl import check_value, is_array_type, get_child, find_node, \
        is_open_node, is_leaf, get_leaf_attribs
from utils import expect, run_cmd_no_fail

ATMCHANGE_SEP = "-ATMCHANGE_SEP-"
ATMCHANGE_BUFF_XML_NAME = "SCREAM_ATMCHANGE_BUFFER"

###############################################################################
def apply_atm_procs_list_changes_from_buffer(case, xml):
###############################################################################
    atmchg_buffer = case.get_value(ATMCHANGE_BUFF_XML_NAME)
    any_change = False
    if atmchg_buffer:
        atmchgs = unbuffer_changes(case)

        for chg in atmchgs:
            if "atm_procs_list" in chg:
                atm_config_chg_impl(xml, chg)
                any_change = True
    return any_change

###############################################################################
def apply_non_atm_procs_list_changes_from_buffer(case, xml):
###############################################################################
    atmchg_buffer = case.get_value(ATMCHANGE_BUFF_XML_NAME)
    if atmchg_buffer:
        atmchgs = unbuffer_changes(case)

        for chg in atmchgs:
            if "atm_procs_list" not in chg:
                atm_config_chg_impl(xml, chg)

###############################################################################
def buffer_changes(changes):
###############################################################################
    """
    Take a list of raw changes and buffer them in the XML case settings. Raw changes
    are what goes to atm_config_chg_impl.
    """
    # Commas confuse xmlchange and so need to be escaped.
    changes_str = ATMCHANGE_SEP.join(changes).replace(",",r"\,")

    run_cmd_no_fail(f"./xmlchange --append {ATMCHANGE_BUFF_XML_NAME}='{changes_str}{ATMCHANGE_SEP}'")

###############################################################################
def unbuffer_changes(case):
###############################################################################
    """
    From a case, get and return a list of raw changes
    """
    atmchg_buffer = case.get_value(ATMCHANGE_BUFF_XML_NAME)
    atmchgs = []
    for item in atmchg_buffer.split(ATMCHANGE_SEP):
        if item.strip():
            atmchgs.append(item.replace(r"\,", ",").strip())

    return atmchgs

###############################################################################
def reset_buffer():
###############################################################################
    run_cmd_no_fail(f"./xmlchange {ATMCHANGE_BUFF_XML_NAME}=''")

###############################################################################
def get_changes_for_node(xml_root, node_name, changes):
###############################################################################
    """
    Return the subset of changes from `changes` that do NOT target the given
    XML node or any of its descendants.  This is the filtering kernel used by
    reset_node_changes and can be called directly for testing.

    NOTE: `changes` is a list of already-unescaped change strings
          (commas are NOT escaped with backslash).

    >>> xml = '''
    ... <root>
    ...     <foo type="integer">0</foo>
    ...     <bar type="integer">10</bar>
    ...     <sub>
    ...         <child1 type="integer">1</child1>
    ...         <child2 type="integer">2</child2>
    ...     </sub>
    ... </root>
    ... '''
    >>> import xml.etree.ElementTree as ET
    >>> tree = ET.fromstring(xml)
    >>> ################ ERROR: node not found #######################
    >>> get_changes_for_node(tree, 'nonexistent', ['foo=5'])
    Traceback (most recent call last):
    SystemExit: ERROR: 'nonexistent' did not match any node in the XML file
    >>> ################ Reset leaf node (foo), keep bar change #####
    >>> get_changes_for_node(tree, 'foo', ['foo=5', 'bar=2'])
    ['bar=2']
    >>> ################ Reset leaf node (foo), no matching change ##
    >>> get_changes_for_node(tree, 'foo', ['bar=2', 'sub::child1=99'])
    ['bar=2', 'sub::child1=99']
    >>> ################ Reset non-leaf node (sub), removes children
    >>> get_changes_for_node(tree, 'sub', ['sub::child1=99', 'bar=3'])
    ['bar=3']
    >>> ################ Reset child only, keep sibling change #######
    >>> get_changes_for_node(tree, 'sub::child1', ['sub::child1=99', 'sub::child2=5', 'bar=3'])
    ['sub::child2=5', 'bar=3']
    >>> ################ Empty changes list ##########################
    >>> get_changes_for_node(tree, 'foo', [])
    []
    >>> ################ Add/rm changes target the parent node ########
    >>> get_changes_for_node(tree, 'sub', ['sub::new:=1', 'sub::child1~=', 'bar=2'])
    ['bar=2']
    >>> get_changes_for_node(tree, 'sub::child1', ['sub::child1~=', 'sub::new:=1'])
    ['sub::new:=1']
    >>> get_changes_for_node(tree, 'foo', ['sub::new:=1'])
    ['sub::new:=1']
    """
    reset_targets = get_xml_nodes(xml_root, node_name)
    expect(len(reset_targets) > 0,
           f"'{node_name}' did not match any node in the XML file")

    if not changes:
        return []

    parent_map = create_parent_map(xml_root)

    filtered_changes = []
    for chg in changes:
        chg_node_name, _, chg_op = parse_change(chg)
        if chg_op in ["add","rm"]:
            # The leaf may not exist (yet/anymore), so identify the change by its
            # parent node and leaf name
            parent_name, leaf_name = split_parent_name(chg_node_name)
            affects_reset = any(
                is_anchestor_of(reset_target, parent, parent_map) or
                (reset_target.tag==leaf_name and parent_map[reset_target] is parent)
                for parent in get_xml_nodes(xml_root, parent_name)
                for reset_target in reset_targets
            )
            if not affects_reset:
                filtered_changes.append(chg)
            continue

        chg_nodes = get_xml_nodes(xml_root, chg_node_name)
        # is_anchestor_of(A, B, ...) returns True when A == B too, so
        # this covers both direct matches and descendant matches.
        affects_reset = any(
            is_anchestor_of(reset_target, chg_node, parent_map)
            for chg_node in chg_nodes
            for reset_target in reset_targets
        )
        if not affects_reset:
            filtered_changes.append(chg)

    return filtered_changes

###############################################################################
def reset_node_changes(xml_root, node_name):
###############################################################################
    """
    Remove from the SCREAM_ATMCHANGE_BUFFER all changes that target
    the given node or any of its descendants.
    """
    # Get current buffered changes
    buff_str = run_cmd_no_fail(f"./xmlquery {ATMCHANGE_BUFF_XML_NAME} --value")
    changes = []
    for item in buff_str.split(ATMCHANGE_SEP):
        if item.strip():
            changes.append(item.replace(r"\,", ",").strip())

    # Filter changes (also validates node_name exists in xml_root)
    filtered_changes = get_changes_for_node(xml_root, node_name, changes)

    # Reset buffer and write back the filtered changes
    run_cmd_no_fail(f"./xmlchange {ATMCHANGE_BUFF_XML_NAME}=''")
    if filtered_changes:
        buffer_changes(filtered_changes)

###############################################################################
def get_xml_nodes(xml_root, name):
###############################################################################
    """
    Find all elements matching a name where name uses '::' syntax

    >>> xml = '''
    ... <root>
    ...     <prop1>one</prop1>
    ...     <sub>
    ...         <prop1>two</prop1>
    ...         <prop2 type="integer" valid_values="1,2">2</prop2>
    ...     </sub>
    ... </root>
    ... '''
    >>> import xml.etree.ElementTree as ET
    >>> tree = ET.fromstring(xml)
    >>> ################ INVALID SYNTAX #######################
    >>> get_xml_nodes(tree,'sub::::prop1')
    Traceback (most recent call last):
    SystemExit: ERROR: Invalid xml node name format, 'sub::::prop1' contains '::::'
    >>> ################ VALID USAGE #######################
    >>> get_xml_nodes(tree,'invalid::prop1')
    []
    >>> [item.text for item in get_xml_nodes(tree,'ANY::prop1')]
    ['one', 'two']
    >>> [item.text for item in get_xml_nodes(tree,'::prop1')]
    ['one']
    >>> [item.text for item in get_xml_nodes(tree,'prop2')]
    ['2']
    >>> item = get_xml_nodes(tree,'prop2')[0]
    >>> parent_map = create_parent_map(tree)
    >>> [p.tag for p in get_parents(item, parent_map)]
    ['sub']
    """
    expect('::::' not in name,
            f"Invalid xml node name format, '{name}' contains '::::'")

    tokens = name.split("::")
    expect (tokens[-1] != '', "Input query string ends with '::'. It should end with an actual node name")
    if 'ANY' in tokens:
        multiple_hits_ok = True

        # Check there's only ONE 'ANY' token
        expect(tokens.count('ANY') == 1, "Invalid xml node name format, multiple 'ANY' tokens found.")

        # Split tokens list into two parts: before and after 'ANY'
        before_any = tokens[:tokens.index('ANY')]
        after_any = tokens[tokens.index('ANY') + 1:]

        # The case where name starts with ::ANY is delicate, since before_any=[''], and this
        # trips the call to get_xml_nodes. Since ANY and ::ANY are conceptually the same,
        # we set before_any=[] if before_any==['']
        if before_any == ['']:
            before_any = []

        # Call get_xml_nodes with '::'.join(before_any) to get the new root
        new_root = get_xml_nodes(xml_root, '::'.join(before_any)) if before_any else [xml_root]

        # Reset xml_root to new_root for the next search
        xml_root = new_root[0]

        # Use new_root to find all matches for whatever comes after 'ANY::'
        name = '*' if not after_any else '::'.join(after_any)
    else:
        multiple_hits_ok = False

    if name.startswith("::"):
        prefix = "./"  # search immediate children only
        name = name[2:]
    else:
        prefix = ".//"  # search entire tree


    # Handle case without ANY
    try:
        xpath_str = prefix + name.replace("::", "/")
        result = xml_root.findall(xpath_str)
    except SyntaxError as e:
        expect(False, f"Invalid syntax '{name}' -> {e}")

    # Note: don't check that len(result)>0, since user may be ok with 0 matches
    if not multiple_hits_ok and len(result)>1:
        parent_map = create_parent_map(xml_root)
        error_str = f"{name} is ambiguous. Use ANY in the node path to allow multiple matches. Matches:\n"
        for node in result:
            parents = get_parents(node, parent_map)
            name = "::".join(e.tag for e in parents)
            name = node.tag if name=="" else name + "::" + node.tag
            error_str += "  " + name + "\n"

        expect(False, error_str)

    return result

###############################################################################
def modify_ap_list(group, ap_list_str, append_this, remove_this=False, defaults_xml=None):
###############################################################################
    """
    Modify the atm_procs_list entry of this XML node (which is an atm proc group).
    This routine can only be used to add an atm proc group OR to remove some
    atm procs.
    NOTE: defaults_xml is not None only in doctest tests. In regular use cases,
    we load the defaults from the defaults xml file in eamxx/cime_config
    >>> xml = '''
    ... <dummy_defaults>
    ...     <atmosphere_processes_defaults>
    ...         <atm_proc_group>
    ...             <atm_procs_list type="array(string)"/>
    ...         </atm_proc_group>
    ...         <p1>
    ...             <my_param>1</my_param>
    ...         </p1>
    ...         <p2>
    ...             <my_param>2</my_param>
    ...         </p2>
    ...     </atmosphere_processes_defaults>
    ... </dummy_defaults>
    ... '''
    >>> from eamxx_buildnml_impl import has_child
    >>> import xml.etree.ElementTree as ET
    >>> defaults = ET.fromstring(xml)
    >>> node = ET.Element("my_group")
    >>> node.append(ET.Element("atm_procs_list"))
    >>> get_child(node,"atm_procs_list").text = ""
    >>> modify_ap_list(node,"p1,p2",False,False,defaults)
    True
    >>> get_child(node,"atm_procs_list").text
    'p1,p2'
    >>> modify_ap_list(node,"p1",True,False,defaults)
    True
    >>> get_child(node,"atm_procs_list").text
    'p1,p2,p1'
    >>> modify_ap_list(node,"p1",False,True,defaults)
    True
    >>> get_child(node,"atm_procs_list").text
    'p2,p1'
    >>> modify_ap_list(node,"p1,p3",False,False,defaults)
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot modify ap list for group 'my_group'
    Process 'p3' not found in XML tree 'dummy_defaults'
    >>> modify_ap_list(node,"p3",False,True,defaults)
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot remove 'p3' from atm_procs_list of group 'my_group': not found in current list 'p2,p1'
    """
    curr_apl = get_child(group,"atm_procs_list")

    ap_list = ap_list_str.split(",")

    if remove_this:
        curr_list = curr_apl.text.split(",") if curr_apl.text else []
        for ap in ap_list:
            expect(ap in curr_list,
                   f"Cannot remove '{ap}' from atm_procs_list of group '{group.tag}': "
                   f"not found in current list '{curr_apl.text}'")
            curr_list.remove(ap)
        new_text = ','.join(curr_list)
        if curr_apl.text == new_text:
            return False
        curr_apl.text = new_text
        return True

    if curr_apl.text==ap_list_str:
        return False

    expect (len(ap_list)==len(set(ap_list)),
            "Input list of atm procs contains repetitions")

    # To avoid buffering a change that has an invalid atm proc name,
    # we load the defaults here, and check that the added procs exist
    if defaults_xml is None:
        defaults_xml_file = pathlib.Path(__file__).parent.parent.resolve() / "cime_config/namelist_defaults_eamxx.xml"
        with open(defaults_xml_file, "r") as fd:
            defaults_xml = ET.parse(fd).getroot()
    procs_defaults = get_child(defaults_xml,"atmosphere_processes_defaults")
    for ap in ap_list:
        expect (find_node(procs_defaults,ap) is not None,
                f"Cannot modify ap list for group '{group.tag}'\n"
                f"Process '{ap}' not found in XML tree '{defaults_xml.tag}'")

    # Update the 'atm_procs_list' in this node
    if append_this:
        curr_apl.text = ','.join(curr_apl.text.split(",")+ap_list)
    else:
        curr_apl.text = ','.join(ap_list)

    return True

###############################################################################
def is_locked_impl(node):
###############################################################################
    return "locked" in node.attrib.keys() and str(node.attrib["locked"]).upper() == "TRUE"

###############################################################################
def is_locked(xml_root, node):
###############################################################################
    if is_locked_impl(node):
        return True
    else:
        parent_map = create_parent_map(xml_root)
        parents = get_parents(node, parent_map)
        for parent in parents:
            if is_locked_impl(parent):
                return True

    return False

###############################################################################
def split_parent_name(name):
###############################################################################
    """
    Split a node name into the name of its parent node and its own name.

    >>> split_parent_name('a::b::c')
    ('a::b', 'c')
    >>> split_parent_name('ANY::b')
    ('ANY', 'b')
    >>> split_parent_name('b')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot add/remove leaf 'b': the name must include the name of the parent open node (e.g., 'parent::b')
    """
    parent_name, sep, leaf_name = name.rpartition("::")
    expect (sep!="" and parent_name!="",
            f"Cannot add/remove leaf '{name}': the name must include the name of the parent open node (e.g., 'parent::{leaf_name or name}')")
    return parent_name, leaf_name

###############################################################################
def check_can_edit_leaves(xml_root, parent, leaf_name, action):
###############################################################################
    expect (is_open_node(parent),
            f"Cannot {action} leaf '{leaf_name}': '{parent.tag}' is not an open node.\n"
            "Only nodes marked as open=\"true\" in the defaults file allow adding/removing leaves.")
    expect (not is_locked(xml_root, parent),
            f"Cannot change {parent.tag}, it is locked")

###############################################################################
def add_leaf(xml_root, name, value):
###############################################################################
    """
    Add a new leaf to an open node. The leaf inherits its metadata from the
    leaf_* attributes of the open node. Returns True (a leaf was added).
    """
    parent_name, leaf_name = split_parent_name(name)
    parents = get_xml_nodes(xml_root, parent_name)
    expect (len(parents)>0, f"Cannot add leaf '{name}': '{parent_name}' did not match any node")
    expect (len(parents)==1,
            f"Cannot add leaf '{name}': '{parent_name}' matches multiple nodes. Please, be more specific.")
    parent = parents[0]

    check_can_edit_leaves(xml_root, parent, leaf_name, "add")
    expect (re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", leaf_name) is not None and leaf_name!="ANY",
            f"Invalid leaf name '{leaf_name}'. Names must be alphanumeric (underscores allowed) and cannot start with a digit.")
    expect (parent.find(leaf_name) is None,
            f"Cannot add leaf '{leaf_name}' to '{parent.tag}': a leaf with that name already exists.\n"
            f"  Use '{name}=<value>' to modify it.")

    leaf = ET.SubElement(parent, leaf_name)
    leaf.attrib.update(get_leaf_attribs(parent))
    expect ("type" in leaf.attrib,
            f"Open node '{parent.tag}' does not specify 'leaf_type'. Please, contact developers.")

    try:
        check_value(leaf, value)
    except BaseException:
        # Don't leave a half-baked leaf in the tree
        parent.remove(leaf)
        raise
    leaf.text = value

    return True

###############################################################################
def rm_leaf(xml_root, name):
###############################################################################
    """
    Remove leaves from open nodes. Returns True (a leaf was removed).
    """
    matches = get_xml_nodes(xml_root, name)
    expect (len(matches)>0, f"{name} did not match any items")

    parent_map = create_parent_map(xml_root)
    for node in matches:
        expect (is_leaf(node), f"Cannot remove '{node.tag}': it is not a leaf")
        parent = parent_map[node]
        check_can_edit_leaves(xml_root, parent, node.tag, "remove")
        parent.remove(node)

    return True

###############################################################################
def apply_change(xml_root, node, new_value, append_this, remove_this=False):
###############################################################################
    any_change = False

    # User can change the list of atm procs in a group doing ./atmchange group_name=a,b,c
    # If we detect that this node is an atm proc group, don't modify the text, but do something els
    if node.tag=="atm_procs_list":
        parent_map = create_parent_map(xml_root)
        group = get_parents(node,parent_map)[-1]
        return modify_ap_list(group, new_value, append_this, remove_this)

    if append_this:

        expect (not is_locked(xml_root, node), f"Cannot change {node.tag}, it is locked")
        expect ("type" in node.attrib.keys(),
                f"Error! Missing type information for {node.tag}")
        type_ = node.attrib["type"]
        expect (is_array_type(type_) or type_=="string",
                "Error! Can only append with array and string types.\n"
                f"    - name: {node.tag}\n"
                f"    - type: {type_}")

        if node.text is None:
            node.text = ""

        if is_array_type(type_) and node.text!="":
            node.text += ", " + new_value
        else:
            node.text += new_value

        any_change = True

    elif remove_this:

        expect (not is_locked(xml_root, node), f"Cannot change {node.tag}, it is locked")
        expect ("type" in node.attrib.keys(),
                f"Error! Missing type information for {node.tag}")
        type_ = node.attrib["type"]
        expect (is_array_type(type_),
                "Error! Can only remove with array types.\n"
                f"    - name: {node.tag}\n"
                f"    - type: {type_}")

        curr_list = [v.strip() for v in node.text.split(",")] if node.text else []
        remove_list = [v.strip() for v in new_value.split(",")]
        for v in remove_list:
            expect (v in curr_list,
                    f"Error! Value '{v}' not found in {node.tag}. "
                    f"Current value: {node.text}")
            curr_list.remove(v)
        node.text = ",".join(curr_list)
        any_change = True

    elif node.text != new_value:
        expect (not is_locked(xml_root, node), f"Cannot change {node.tag}, it is locked")
        check_value(node,new_value)
        node.text = new_value
        any_change = True

    return any_change

###############################################################################
# A change request. The 'op' field can be one of
#  - set    : set the value of an existing leaf
#  - append : append to an existing array (or string) leaf
#  - remove : remove entries from an existing array leaf
#  - add    : add a new leaf to an open node
#  - rm     : remove a leaf from an open node (value must be empty)
# The 'add' and 'rm' operations are internal: they are how the --add and --rm
# options of atmchange are stored in the buffer (as 'A::B:=value' and 'A::B~='),
# so they can be replayed. They are not accepted directly from the command line.
Change = namedtuple("Change", ["name", "value", "op"])
CHANGE_OPS = {"": "set", "+": "append", "-": "remove", ":": "add", "~": "rm"}
CHANGE_OPS_INV = {v:k for k,v in CHANGE_OPS.items()}
CHANGE_RE = re.compile(r"^([^=]*?)([+\-:~]?)=(.*)$", re.DOTALL)

###############################################################################
def parse_change(change):
###############################################################################
    """
    >>> parse_change("a+=2")
    Change(name='a', value='2', op='append')
    >>> parse_change("a-=2")
    Change(name='a', value='2', op='remove')
    >>> parse_change("a=hello")
    Change(name='a', value='hello', op='set')
    >>> parse_change("a::b:=1,2")
    Change(name='a::b', value='1,2', op='add')
    >>> parse_change("a::b~=")
    Change(name='a::b', value='', op='rm')
    >>> parse_change("a::b=name=2")
    Change(name='a::b', value='name=2', op='set')
    >>> parse_change("a::b~=2")
    Traceback (most recent call last):
    SystemExit: ERROR: Invalid change request 'a::b~=2'. The 'rm' operation does not accept a value
    >>> parse_change("a::=2")
    Traceback (most recent call last):
    SystemExit: ERROR: Invalid change request 'a::=2'. The node name cannot end with '::'
    """
    m = CHANGE_RE.match(change)
    expect (m is not None and m.group(1)!="",
        f"Invalid change request '{change}'. Valid formats are:\n"
        f"  - A[::B[...]=value\n"
        f"  - A[::B[...]+=value  (implies append for this change)\n"
        f"  - A[::B[...]-=value  (implies removal for this change, arrays only)")

    name, op, value = m.group(1), CHANGE_OPS[m.group(2)], m.group(3)
    expect (not name.endswith(":"),
            f"Invalid change request '{change}'. The node name cannot end with '::'")
    expect (op!="rm" or value=="",
            f"Invalid change request '{change}'. The 'rm' operation does not accept a value")

    return Change(name, value, op)

###############################################################################
def format_change(name, value, op):
###############################################################################
    """
    Inverse of parse_change

    >>> format_change('a::b','1','add')
    'a::b:=1'
    >>> parse_change(format_change('a::b','','rm'))
    Change(name='a::b', value='', op='rm')
    """
    return f"{name}{CHANGE_OPS_INV[op]}={value}"

###############################################################################
def atm_config_chg_impl(xml_root, change):
###############################################################################
    """
    >>> xml = '''
    ... <root>
    ...   <a type="array(int)">1,2,3</a>
    ...   <b type="array(int)">1</b>
    ...   <c type="int">1</c>
    ...   <d type="string">one</d>
    ...   <e type="array(string)">one</e>
    ...   <prop1>one</prop1>
    ...   <sub>
    ...     <prop1>two</prop1>
    ...     <prop2 type="integer" valid_values="1,2">2</prop2>
    ...   </sub>
    ...   <sub2 locked="true">
    ...     <subsub2>
    ...       <subsubsub2>
    ...         <lprop2>hi</lprop2>
    ...       </subsubsub2>
    ...     </subsub2>
    ...   </sub2>
    ...   <sub3>
    ...     <subsub3>
    ...       <subsubsub3 locked="true">
    ...         <lprop3>hi</lprop3>
    ...       </subsubsub3>
    ...     </subsub3>
    ...   </sub3>
    ...   <sub4>
    ...     <subsub4>
    ...       <subsubsub4>
    ...         <lprop4 locked="true">hi</lprop4>
    ...       </subsubsub4>
    ...     </subsub4>
    ...   </sub4>
    ... </root>
    ... '''
    >>> import xml.etree.ElementTree as ET
    >>> tree = ET.fromstring(xml)
    >>> ################ INVALID SYNTAX #######################
    >>> atm_config_chg_impl(tree,'prop1->2')
    Traceback (most recent call last):
    SystemExit: ERROR: Invalid change request 'prop1->2'. Valid formats are:
      - A[::B[...]=value
      - A[::B[...]+=value  (implies append for this change)
      - A[::B[...]-=value  (implies removal for this change, arrays only)
    >>> ################ INVALID TYPE #######################
    >>> atm_config_chg_impl(tree,'prop2=two')
    Traceback (most recent call last):
    CIME.core.exceptions.CIMEError: ERROR: Could not refine 'two' as type 'integer':
    could not convert string to float: 'two'
    >>> ################ INVALID VALUE #######################
    >>> atm_config_chg_impl(tree,'prop2=3')
    Traceback (most recent call last):
    CIME.core.exceptions.CIMEError: ERROR: Invalid value '3' for element 'prop2'. Value not in the valid list ('[1, 2]')
    >>> ################ AMBIGUOUS CHANGE #######################
    >>> atm_config_chg_impl(tree,'prop1=three')
    Traceback (most recent call last):
    SystemExit: ERROR: prop1 is ambiguous. Use ANY in the node path to allow multiple matches. Matches:
      prop1
      sub::prop1
    <BLANKLINE>
    >>> ################ VALID USAGE #######################
    >>> atm_config_chg_impl(tree,'::prop1=two')
    True
    >>> atm_config_chg_impl(tree,'::prop1=two')
    False
    >>> atm_config_chg_impl(tree,'sub::prop1=one')
    True
    >>> atm_config_chg_impl(tree,'ANY::prop1=three')
    True
    >>> [item.text for item in get_xml_nodes(tree,'ANY::prop1')]
    ['three', 'three']
    >>> ################ TEST APPEND += #################
    >>> atm_config_chg_impl(tree,'a+=4')
    True
    >>> get_xml_nodes(tree,'a')[0].text
    '1,2,3, 4'
    >>> ################ ERROR, append to non-array and non-string
    >>> atm_config_chg_impl(tree,'c+=2')
    Traceback (most recent call last):
    SystemExit: ERROR: Error! Can only append with array and string types.
        - name: c
        - type: int
    >>> ################ Append to string ##################
    >>> atm_config_chg_impl(tree,'d+=two')
    True
    >>> get_xml_nodes(tree,'d')[0].text
    'onetwo'
    >>> ################ Append to array(string) ##################
    >>> atm_config_chg_impl(tree,'e+=two')
    True
    >>> get_xml_nodes(tree,'e')[0].text
    'one, two'
    >>> ################ TEST REMOVE -= #################
    >>> atm_config_chg_impl(tree,'a-=2')
    True
    >>> get_xml_nodes(tree,'a')[0].text
    '1,3,4'
    >>> atm_config_chg_impl(tree,'a-=1,3')
    True
    >>> get_xml_nodes(tree,'a')[0].text
    '4'
    >>> ################ ERROR, remove from non-array
    >>> atm_config_chg_impl(tree,'c-=1')
    Traceback (most recent call last):
    SystemExit: ERROR: Error! Can only remove with array types.
        - name: c
        - type: int
    >>> ################ ERROR, remove value not in list
    >>> atm_config_chg_impl(tree,'b-=9')
    Traceback (most recent call last):
    SystemExit: ERROR: Error! Value '9' not found in b. Current value: 1
    >>> ################ Test locked ##################
    >>> atm_config_chg_impl(tree, 'lprop2=yo')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot change lprop2, it is locked
    >>> atm_config_chg_impl(tree, 'lprop3=yo')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot change lprop3, it is locked
    >>> atm_config_chg_impl(tree, 'lprop4=yo')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot change lprop4, it is locked
    >>> ################ Test open nodes: add/rm leaves ##################
    >>> xml = '''
    ... <root>
    ...   <closed><x type="real">1.0</x></closed>
    ...   <open_empty open="true" leaf_type="array(real)"/>
    ...   <open_pos open="true" leaf_type="real" leaf_constraints="ge 0"/>
    ...   <open_full open="true" leaf_type="array(real)">
    ...     <f1 type="array(real)">1.0</f1>
    ...   </open_full>
    ...   <open_locked open="true" locked="true" leaf_type="real"/>
    ...   <open_dup1 open="true" leaf_type="real"/>
    ...   <open_dup2 open="true" leaf_type="real"/>
    ... </root>
    ... '''
    >>> tree = ET.fromstring(xml)
    >>> atm_config_chg_impl(tree,'open_empty::a:=1,2,3')
    True
    >>> [(c.tag,c.text,c.attrib) for c in get_xml_nodes(tree,'open_empty')[0]]
    [('a', '1,2,3', {'type': 'array(real)'})]
    >>> ###### The new leaf can be modified/appended to, as any other leaf
    >>> atm_config_chg_impl(tree,'open_empty::a=4')
    True
    >>> atm_config_chg_impl(tree,'open_empty::a+=5')
    True
    >>> get_xml_nodes(tree,'open_empty::a')[0].text
    '4, 5'
    >>> atm_config_chg_impl(tree,'open_full::f2:=7')
    True
    >>> [c.tag for c in get_xml_nodes(tree,'open_full')[0]]
    ['f1', 'f2']
    >>> ###### Leaves from the defaults can be removed, and added back
    >>> atm_config_chg_impl(tree,'open_full::f1~=')
    True
    >>> [c.tag for c in get_xml_nodes(tree,'open_full')[0]]
    ['f2']
    >>> atm_config_chg_impl(tree,'open_full::f1:=1')
    True
    >>> ###### Removing all leaves leaves the (open) node in place
    >>> atm_config_chg_impl(tree,'open_full::f1~=')
    True
    >>> atm_config_chg_impl(tree,'open_full::f2~=')
    True
    >>> len(get_xml_nodes(tree,'open_full')[0])
    0
    >>> ################ ERRORS: adding leaves #####################
    >>> atm_config_chg_impl(tree,'closed::y:=1')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot add leaf 'y': 'closed' is not an open node.
    Only nodes marked as open="true" in the defaults file allow adding/removing leaves.
    >>> atm_config_chg_impl(tree,'closed::x:=1')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot add leaf 'x': 'closed' is not an open node.
    Only nodes marked as open="true" in the defaults file allow adding/removing leaves.
    >>> atm_config_chg_impl(tree,'open_empy::a:=1')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot add leaf 'open_empy::a': 'open_empy' did not match any node
    >>> atm_config_chg_impl(tree,'a:=1')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot add/remove leaf 'a': the name must include the name of the parent open node (e.g., 'parent::a')
    >>> atm_config_chg_impl(tree,'open_empty::a:=1')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot add leaf 'a' to 'open_empty': a leaf with that name already exists.
      Use 'open_empty::a=<value>' to modify it.
    >>> atm_config_chg_impl(tree,'open_pos::b:=-1')
    Traceback (most recent call last):
    CIME.core.exceptions.CIMEError: ERROR: Value '-1.0' for entry 'b' violates constraint '-1.0 >= 0.0'
    >>> len(get_xml_nodes(tree,'open_pos')[0])
    0
    >>> atm_config_chg_impl(tree,'open_pos::b:=1')
    True
    >>> get_xml_nodes(tree,'open_pos::b')[0].attrib
    {'type': 'real', 'constraints': 'ge 0'}
    >>> atm_config_chg_impl(tree,'open_empty::b:=hi')
    Traceback (most recent call last):
    CIME.core.exceptions.CIMEError: ERROR: Could not refine 'hi' as type 'real':
    could not convert string to float: 'hi'
    >>> atm_config_chg_impl(tree,'open_empty::1b:=1')
    Traceback (most recent call last):
    SystemExit: ERROR: Invalid leaf name '1b'. Names must be alphanumeric (underscores allowed) and cannot start with a digit.
    >>> atm_config_chg_impl(tree,'open_locked::a:=1')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot change open_locked, it is locked
    >>> atm_config_chg_impl(tree,'ANY::a:=1')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot add leaf 'ANY::a': 'ANY' matches multiple nodes. Please, be more specific.
    >>> ################ ERRORS: removing leaves #####################
    >>> atm_config_chg_impl(tree,'closed::x~=')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot remove leaf 'x': 'closed' is not an open node.
    Only nodes marked as open="true" in the defaults file allow adding/removing leaves.
    >>> atm_config_chg_impl(tree,'open_empty::nope~=')
    Traceback (most recent call last):
    SystemExit: ERROR: open_empty::nope did not match any items
    >>> atm_config_chg_impl(tree,'open_empty~=')
    Traceback (most recent call last):
    SystemExit: ERROR: Cannot remove 'open_empty': it is not a leaf
    """
    node_name, new_value, op = parse_change(change)

    # Adding/removing leaves is dealt with separately: the leaf may not exist
    if op=="add":
        return add_leaf(xml_root, node_name, new_value)
    elif op=="rm":
        return rm_leaf(xml_root, node_name)

    append_this = op=="append"
    remove_this = op=="remove"
    matches = get_xml_nodes(xml_root, node_name)

    expect(len(matches) > 0, f"{node_name} did not match any items")

    any_change = False
    for node in matches:
        any_change |= apply_change(xml_root, node, new_value, append_this, remove_this)

    return any_change

###############################################################################
def create_parent_map(root):
###############################################################################
    pmap = {c: p for p in root.iter() for c in p}
    pmap[root] = None
    return pmap

###############################################################################
def get_parents(elem, parent_map):
###############################################################################
    """
    Return all parents of an elem in descending order (first item in list will
    be the furthest ancestor, last item will be direct parent)
    """
    results = []
    if not is_root(elem,parent_map):
        parent = parent_map[elem]
        results = get_parents(parent, parent_map)
        if not is_root(parent,parent_map):
            results.append(parent)

    return results

###############################################################################
def is_anchestor_of(parent,child,parent_map):
###############################################################################
    """
    >>> xml = '''
    ... <root>
    ...     <prop1>one</prop1>
    ...     <sub>
    ...         <prop1>two</prop1>
    ...         <prop2 type="integer" valid_values="1,2">2</prop2>
    ...     </sub>
    ... </root>
    ... '''
    >>> import xml.etree.ElementTree as ET
    >>> tree = ET.fromstring(xml)
    >>> parent_map = create_parent_map(tree)
    >>> sub = get_xml_nodes(tree,'sub')[0]
    >>> sub_prop1 = get_xml_nodes(tree,'sub::prop1')[0]
    >>> is_anchestor_of(sub,sub_prop1,parent_map)
    True
    >>> is_anchestor_of(sub_prop1,sub,parent_map)
    False
    """
    curr = child
    while curr is not None:
        if curr is parent:
            return True
        curr = parent_map[curr]
    return False

###############################################################################
def is_root (node,parent_map):
###############################################################################
    return parent_map[node] is None

###############################################################################
def print_var_impl(node,parent_map,full,dtype,value,valid_values,print_style="invalid",indent=""):
###############################################################################

    if not is_leaf(node):
        print (f"{indent}{node.tag}:")
        # This is not a leaf, so print all nested nodes.
        for child in node:
            # Since prints are nicely nested, use 'node-name' as print style
            print_var_impl(child,{},full,dtype,value,valid_values,'node-name',indent+"    ")
        return

    expect (print_style in ["node-name","full-scope","parent-scope"],
            f"Invalid print_style '{print_style}' for print_var_impl. Use 'full' or 'short'.")

    if print_style=="node-name":
        # Just the inner most name
        name = node.tag
    elif print_style=="parent-scope":
        parent = parent_map[node]
        name = node.tag if is_root(parent,parent_map) else parent.tag + "::" + node.tag
    else:
        parents = get_parents(node, parent_map)
        name = "::".join(e.tag for e in parents)
        name += "::" if parents else ""
        name += node.tag

    if full:
        expect ("type" in node.attrib.keys(),
                f"Error! Missing type information for {name}")
        print (f"{indent}{name}")
        print (f"{indent}    value: {node.text}")
        print (f"{indent}    type: {node.attrib['type']}")
        if "valid_values" not in node.attrib.keys():
            valid = []
        else:
            valid = node.attrib["valid_values"].split(",")
        print (f"{indent}    valid values: {valid}")
    elif dtype:
        expect ("type" in node.attrib.keys(),
                f"Error! Missing type information for {name}")
        print (f"{indent}{name}: {node.attrib['type']}")
    elif value:
        print (f"{indent}{node.text}")
    elif valid_values:
        if "valid_values" not in node.attrib.keys():
            valid = '<valid values not provided>'
        else:
            valid = node.attrib["valid_values"].split(",")
        print (f"{indent}{name}: {valid}")
    else:
        print (f"{indent}{name}: {node.text}")

###############################################################################
def print_var(xml_root,parent_map,var,full,dtype,value,valid_values,print_style="invalid",indent=""):
###############################################################################
    """
    >>> xml = '''
    ... <root>
    ...     <prop1>one</prop1>
    ...     <sub>
    ...         <prop1>two</prop1>
    ...         <prop2 type="integer" valid_values="1,2">2</prop2>
    ...     </sub>
    ... </root>
    ... '''
    >>> import xml.etree.ElementTree as ET
    >>> tree = ET.fromstring(xml)
    >>> parent_map = create_parent_map(tree)
    >>> ################ Missing type data #######################
    >>> print_var(tree,parent_map,'::prop1',False,True,False,False,"node-name")
    Traceback (most recent call last):
    SystemExit: ERROR: Error! Missing type information for prop1
    >>> print_var(tree,parent_map,'prop2',True,False,False,False,"node-name")
    prop2
        value: 2
        type: integer
        valid values: ['1', '2']
    >>> print_var(tree,parent_map,'prop2',False,True,False,False,"node-name")
    prop2: integer
    >>> print_var(tree,parent_map,'prop2',False,False,True,False,"node-name")
    2
    >>> print_var(tree,parent_map,'prop2',False,False,False,True,"node-name","    ")
        prop2: ['1', '2']
    """

    expect (print_style in ["node-name","full-scope","parent-scope"],
            f"Invalid print_style '{print_style}' for print_var. Use 'node-name', 'full-scope', or 'parent-scope'.")

    # Get matches
    matches = get_xml_nodes(xml_root,var)

    # If ANY is in the var name, we may hit the case where one of the matches
    # is a parent of another match. In this case, we want to get rid of 
    unique_matches = []
    for i in matches:
        add_this = True
        for j in matches:
            if i is not j and is_anchestor_of(j,i,parent_map):
                add_this = False
            
        if add_this:
            unique_matches.append(i)

    for node in unique_matches:
        print_var_impl(node,parent_map,full,dtype,value,valid_values,print_style,indent)

###############################################################################
def atm_query_impl(xml_root,variables,listall=False,full=False,value=False,
                   dtype=False, valid_values=False, grep=False,parent_map=None):
###############################################################################
    """
    >>> xml = '''
    ... <root>
    ...     <prop1>one</prop1>
    ...     <sub>
    ...         <prop1>two</prop1>
    ...         <prop2 type="integer" valid_values="1,2">2</prop2>
    ...     </sub>
    ... </root>
    ... '''
    >>> import xml.etree.ElementTree as ET
    >>> tree = ET.fromstring(xml)
    >>> vars = ['prop2','::prop1']
    >>> success = atm_query_impl(tree, vars)
        sub::prop2: 2
        prop1: one
    >>> success = atm_query_impl(tree, [], listall=True, valid_values=True)
        prop1: <valid values not provided>
        sub:
            prop1: <valid values not provided>
            prop2: ['1', '2']
    >>> success = atm_query_impl(tree,['prop1'], grep=True)
        prop1: one
        sub::prop1: two
    """

    if not parent_map:
        parent_map = create_parent_map(xml_root)
    if listall:
        print_var(xml_root,parent_map,'ANY',full,dtype,value,valid_values,"node-name","    ")

    elif grep:
        for regex in variables:
            expect("::" not in regex, "query --grep does not support including parent info")
            var_re = re.compile(f'{regex}')
            if var_re.search(xml_root.tag):
                if len(xml_root)>0:
                    parents = get_parents(xml_root,parent_map)
                    print (f"{'::'.join([p.tag for p in parents]) + '::' + xml_root.tag}:")
                print_var(xml_root,parent_map,'ANY',full,dtype,value,valid_values,"node-name","    ")
            else:
                for elem in xml_root:
                    if len(elem)>0:
                        atm_query_impl(elem,variables,listall,full,value,dtype,valid_values,grep,parent_map)
                    else:
                        if var_re.search(elem.tag):
                            nodes = get_xml_nodes(xml_root, "::"+elem.tag)
                            expect(len(nodes) == 1, "No matches?")
                            print_var_impl(nodes[0],parent_map,full,dtype,value,valid_values,"full-scope","    ")

    else:
        for var in variables:
            pmap = {} if var=='ANY' else parent_map
            print_var(xml_root,pmap,var,full,dtype,value,valid_values,"parent-scope","    ")

    return True
