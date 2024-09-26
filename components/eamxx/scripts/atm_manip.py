"""
Retrieve nodes from EAMxx XML config file.
"""

import sys, os, re

# Used for doctests
import xml.etree.ElementTree as ET # pylint: disable=unused-import

# Add path to cime_config folder
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "cime_config"))
from eamxx_buildnml_impl import check_value, is_array_type, get_child, find_node
from eamxx_buildnml_impl import gen_atm_proc_group
from utils import expect, run_cmd_no_fail

ATMCHANGE_SEP = "-ATMCHANGE_SEP-"
ATMCHANGE_BUFF_XML_NAME = "SCREAM_ATMCHANGE_BUFFER"

###############################################################################
def apply_atm_procs_list_changes_from_buffer(case, xml):
###############################################################################
    atmchg_buffer = case.get_value(ATMCHANGE_BUFF_XML_NAME)
    if atmchg_buffer:
        atmchgs = unbuffer_changes(case)

        for chg in atmchgs:
            if "atm_procs_list" in chg:
                atm_config_chg_impl(xml, chg)

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
    SystemExit: ERROR: Invalid xml node name format, 'sub::::prop1' contains ::::
    >>> ################ VALID USAGE #######################
    >>> get_xml_nodes(tree,'invalid::prop1')
    []
    >>> [item.text for item in get_xml_nodes(tree,'prop1')]
    ['one', 'two']
    >>> [item.text for item in get_xml_nodes(tree,'::prop1')]
    ['one']
    >>> [item.text for item in get_xml_nodes(tree,'prop2')]
    ['2']
    >>> item = get_xml_nodes(tree,'prop2')[0]
    >>> parent_map = create_parent_map(tree)
    >>> [p.tag for p in get_parents(item, parent_map)]
    ['root', 'sub']
    """
    expect('::::' not in name,
            "Badly formatted node name: found '::::'")

    tokens = name.split("::")
    expect (tokens[-1] != '', "Input query string ends with '::'. It should end with an actual node name")
    if 'ANY' in tokens:
        multiple_hits_ok = True

        # Check there's only ONE 'ANY' token
        expect(tokens.count('ANY') == 1, "Invalid xml node name format, multiple 'ANY' tokens found.")

        # Split tokens list into two parts: before and after 'ANY'
        before_any = tokens[:tokens.index('ANY')]
        after_any = tokens[tokens.index('ANY') + 1:]
        expect (after_any, "Input name should not end with ANY")

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
        name = '::'.join(after_any)
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
            name = "::".join(e.tag for e in parents) + "::" + node.tag
            error_str += "  " + name + "\n"

        expect(False, error_str)

    return result

###############################################################################
def modify_ap_list(xml_root, group, ap_list_str, append_this):
###############################################################################
    """
    Modify the atm_procs_list entry of this XML node (which is an atm proc group).
    This routine can only be used to add an atm proc group OR to remove some
    atm procs.
    >>> xml = '''
    ... <root>
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
    ... </root>
    ... '''
    >>> from eamxx_buildnml_impl import has_child
    >>> import xml.etree.ElementTree as ET
    >>> tree = ET.fromstring(xml)
    >>> node = ET.Element("my_group")
    >>> node.append(ET.Element("atm_procs_list"))
    >>> get_child(node,"atm_procs_list").text = ""
    >>> modify_ap_list(tree,node,"p1,p2",False)
    True
    >>> get_child(node,"atm_procs_list").text
    'p1,p2'
    >>> modify_ap_list(tree,node,"p1",True)
    True
    >>> get_child(node,"atm_procs_list").text
    'p1,p2,p1'
    >>> modify_ap_list(tree,node,"p1,p3",False)
    Traceback (most recent call last):
    SystemExit: ERROR: Unrecognized atm proc name 'p3'. To declare a new group, prepend and append '_' to the name.
    >>> modify_ap_list(tree,node,"p1,_my_group_",False)
    True
    >>> get_child(node,"atm_procs_list").text
    'p1,_my_group_'
    >>> defaults = get_child(tree,'atmosphere_processes_defaults')
    >>> has_child(defaults,'_my_group_')
    True
    """
    curr_apl = get_child(group,"atm_procs_list")
    if curr_apl.text==ap_list_str:
        return False

    ap_list = ap_list_str.split(",")
    expect (len(ap_list)==len(set(ap_list)),
            "Input list of atm procs contains repetitions")

    # If we're here b/c of a manual call of atmchange from command line, this will be None,
    # since we don't have this node in the genereated XML file. But in that case, we don't
    # have to actually add the new nodes, we can simply just modify the atm_procs_list entry
    # If, however, we're calling this from buildnml, then what we are passed in is the XML
    # tree from namelists_defaults_scream.xml, so this section *will* be present. And we
    # need to add the new atm procs group as children, so that buildnml knows how to build
    # them
    ap_defaults = find_node(xml_root,"atmosphere_processes_defaults")
    if ap_defaults is not None:

        # Figure out which aps in the list are new groups and which ones already
        # exist in the defaults
        add_aps = [n for n in ap_list if n not in curr_apl.text.split(',')]
        new_aps = [n for n in add_aps if find_node(ap_defaults,n) is None]

        for ap in new_aps:
            expect (ap[0]=="_" and ap[-1]=="_" and len(ap)>2,
                    f"Unrecognized atm proc name '{ap}'. To declare a new group, prepend and append '_' to the name.")
            group = gen_atm_proc_group("", ap_defaults)
            group.tag = ap

            ap_defaults.append(group)

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
def apply_change(xml_root, node, new_value, append_this):
###############################################################################
    any_change = False

    # User can change the list of atm procs in a group doing ./atmchange group_name=a,b,c
    # If we detect that this node is an atm proc group, don't modify the text, but do something els
    if node.tag=="atm_procs_list":
        parent_map = create_parent_map(xml_root)
        group = get_parents(node,parent_map)[-1]
        return modify_ap_list (xml_root,group,new_value,append_this)

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

    elif node.text != new_value:
        expect (not is_locked(xml_root, node), f"Cannot change {node.tag}, it is locked")
        check_value(node,new_value)
        node.text = new_value
        any_change = True

    return any_change

###############################################################################
def parse_change(change):
###############################################################################
    """
    >>> parse_change("a+=2")
    ('a', '2', True)
    >>> parse_change("a=hello")
    ('a', 'hello', False)
    """
    tokens = change.split('+=')
    if len(tokens)==2:
        append_this = True
    else:
        append_this = False
        tokens = change.split('=')

    expect (len(tokens)==2,
        f"Invalid change request '{change}'. Valid formats are:\n"
        f"  - A[::B[...]=value\n"
        f"  - A[::B[...]+=value  (implies append for this change)")
    node_name = tokens[0]
    new_value = tokens[1]

    return node_name,new_value,append_this

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
    >>> ################ INVALID TYPE #######################
    >>> atm_config_chg_impl(tree,'prop2=two')
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: Could not refine 'two' as type 'integer':
    could not convert string to float: 'two'
    >>> ################ INVALID VALUE #######################
    >>> atm_config_chg_impl(tree,'prop2=3')
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: Invalid value '3' for element 'prop2'. Value not in the valid list ('[1, 2]')
    >>> ################ AMBIGUOUS CHANGE #######################
    >>> atm_config_chg_impl(tree,'prop1=three')
    Traceback (most recent call last):
    SystemExit: ERROR: prop1 is ambiguous. Use ANY in the node path to allow multiple matches. Matches:
      root::prop1
      root::sub::prop1
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
    >>> [item.text for item in get_xml_nodes(tree,'prop1')]
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
    """
    node_name, new_value, append_this = parse_change(change)
    matches = get_xml_nodes(xml_root, node_name)

    expect(len(matches) > 0, f"{node_name} did not match any items")

    any_change = False
    for node in matches:
        any_change |= apply_change(xml_root, node, new_value, append_this)

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
    if elem in parent_map and parent_map[elem] is not None:
        parent = parent_map[elem]
        results = get_parents(parent, parent_map)
        if parent_map[parent] is not None:
            results.append(parent)

    return results

###############################################################################
def print_var_impl(node,parent_map,full,dtype,value,valid_values,print_style="invalid",indent=""):
###############################################################################

    expect (print_style in ["short","full"],
            f"Invalid print_style '{print_style}' for print_var_impl. Use 'full' or 'short'.")

    if print_style=="short":
        # Just the inner most name
        name = node.tag
    else:
        parents = get_parents(node, parent_map)
        name = "::".join(e.tag for e in parents) + "::" + node.tag

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
    >>> print_var(tree,parent_map,'::prop1',False,True,False,False,"short")
    Traceback (most recent call last):
    SystemExit: ERROR: Error! Missing type information for prop1
    >>> print_var(tree,parent_map,'prop2',True,False,False,False,"short")
    prop2
        value: 2
        type: integer
        valid values: ['1', '2']
    >>> print_var(tree,parent_map,'prop2',False,True,False,False,"short")
    prop2: integer
    >>> print_var(tree,parent_map,'prop2',False,False,True,False,"short")
    2
    >>> print_var(tree,parent_map,'prop2',False,False,False,True,"short","    ")
        prop2: ['1', '2']
    """

    expect (print_style in ["short","full"],
            f"Invalid print_style '{print_style}' for print_var. Use 'full' or 'short'.")

    # Get matches
    matches = get_xml_nodes(xml_root,var)

    for node in matches:
        print_var_impl(node,parent_map,full,dtype,value,valid_values,print_style,indent)

###############################################################################
def print_all_vars(xml_root,xml_node,parent_map,curr_namespace,full,dtype,value,valid_values,print_style,indent):
###############################################################################

    print (f"{indent}{xml_node.tag}")
    for c in xml_node:
        if len(c)>0:
            print_all_vars(xml_root,c,parent_map,curr_namespace+c.tag+"::",full,dtype,value,valid_values,print_style,indent+"    ")
        else:
            print_var(xml_root,parent_map,curr_namespace+c.tag,full,dtype,value,valid_values,print_style,indent+"    ")

###############################################################################
def atm_query_impl(xml_root,variables,listall=False,full=False,value=False,
                   dtype=False, valid_values=False, grep=False):
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
        root::sub::prop2: 2
        root::prop1: one
    >>> success = atm_query_impl(tree, [], listall=True, valid_values=True)
        root
            prop1: <valid values not provided>
            sub
                prop1: <valid values not provided>
                prop2: ['1', '2']
    >>> success = atm_query_impl(tree,['prop1'], grep=True)
        root::prop1: one
        sub::prop1: two
    """
    parent_map = create_parent_map(xml_root)
    if listall:
        print_all_vars(xml_root,xml_root,parent_map,"::",full,dtype,value,valid_values,"short","    ")

    elif grep:
        for regex in variables:
            expect("::" not in regex, "query --grep does not support including parent info")
            var_re = re.compile(f'{regex}')
            if var_re.search(xml_root.tag):
                print_all_vars(xml_root,xml_root,parent_map,"::",full,dtype,value,valid_values,"short","  ")
            else:
                for elem in xml_root:
                    if len(elem)>0:
                        atm_query_impl(elem,variables,listall,full,value,dtype,valid_values,grep)
                    else:
                        if var_re.search(elem.tag):
                            nodes = get_xml_nodes(xml_root, "::"+elem.tag)
                            expect(len(nodes) == 1, "No matches?")
                            print_var_impl(nodes[0],parent_map,full,dtype,value,valid_values,"full","    ")

    else:
        for var in variables:
            print_var(xml_root,parent_map,var,full,dtype,value,valid_values,"full","    ")

    return True
