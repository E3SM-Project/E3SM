"""
Retrieve nodes from EAMxx XML config file.
"""

import sys, os, re

# Used for doctests
import xml.etree.ElementTree as ET # pylint: disable=unused-import

# Add path to cime_config folder
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "cime_config"))
from eamxx_buildnml_impl import check_value, is_array_type
from utils import expect

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
    expect("::::" not in name, f"Invalid xml node name format, '{name}' contains ::::")

    if name.startswith("::"):
        prefix = "./"  # search immediate children only
        name = name[2:]
    else:
        prefix = ".//" # search entire tree

    try:
        xpath_str = prefix + name.replace("::", "/")
        result = xml_root.findall(xpath_str)
    except SyntaxError as e:
        expect(False, f"Invalid syntax '{name}' -> {e}")

    return result

###############################################################################
def apply_change(node, new_value, append_this):
###############################################################################
    any_change = False

    if append_this:
        expect ("type" in node.attrib.keys(),
                f"Error! Missing type information for {node.tag}")
        type_ = node.attrib["type"]
        expect (is_array_type(type_) or type_=="string",
                "Error! Can only append with array and string types.\n"
                f"    - name: {node.tag}\n"
                f"    - type: {type_}")
        if is_array_type(type_):
            node.text += ", " + new_value
        else:
            node.text += new_value

        any_change = True

    elif node.text != new_value:
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
def atm_config_chg_impl(xml_root, change, all_matches=False):
###############################################################################
    """
    >>> xml = '''
    ... <root>
    ...     <a type="array(int)">1,2,3</a>
    ...     <b type="array(int)">1</b>
    ...     <c type="int">1</c>
    ...     <d type="string">one</d>
    ...     <e type="array(string)">one</e>
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
    >>> atm_config_chg_impl(tree,'prop1->2')
    Traceback (most recent call last):
    SystemExit: ERROR: Invalid change request 'prop1->2'. Valid formats are:
      - A[::B[...]=value
      - A[::B[...]+=value  (implies append for this change)
    >>> ################ INVALID TYPE #######################
    >>> atm_config_chg_impl(tree,'prop2=two')
    Traceback (most recent call last):
    ValueError: Could not use 'two' as type 'integer'
    >>> ################ INVALID VALUE #######################
    >>> atm_config_chg_impl(tree,'prop2=3')
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: Invalid value '3' for element 'prop2'. Value not in the valid list ('[1, 2]')
    >>> ################ AMBIGUOUS CHANGE #######################
    >>> atm_config_chg_impl(tree,'prop1=three')
    Traceback (most recent call last):
    SystemExit: ERROR: prop1 is ambiguous, matches:
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
    >>> atm_config_chg_impl(tree,'prop1=three', all_matches=True)
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
    """
    node_name, new_value, append_this = parse_change(change)
    matches = get_xml_nodes(xml_root, node_name)

    expect(len(matches) > 0, f"{node_name} did not match any items")
    if len(matches) > 1 and not all_matches:
        parent_map = create_parent_map(xml_root)
        error_str = ""
        for node in matches:
            parents = get_parents(node, parent_map)
            name = "::".join(e.tag for e in parents) + "::" + node.tag
            error_str += "  " + name + "\n"

        expect(False, f"{node_name} is ambiguous, matches:\n{error_str}")

    any_change = False
    for node in matches:
        any_change |= apply_change(node, new_value, append_this)

    return any_change

###############################################################################
def create_parent_map(root):
###############################################################################
    return {c: p for p in root.iter() for c in p}

###############################################################################
def get_parents(elem, parent_map):
###############################################################################
    """
    Return all parents of an elem in descending order (first item in list will
    be the furthest ancestor, last item will be direct parent)
    """
    results = []
    if elem in parent_map:
        parent = parent_map[elem]
        results = get_parents(parent, parent_map) + [parent]

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

    # Get the shortest unique repr of the var name
    tokens = var.split("::")
    if tokens[0]=='':
        tokens.pop(0)

    while len(tokens)>1:
        new_name = "::".join(tokens[1:])
        matches = get_xml_nodes(xml_root, new_name)
        if len(matches) > 1:
            break
        else:
            tokens.pop(0)

    # Get node, along with all its parents (which might be used for 'full' print style)
    matches = get_xml_nodes(xml_root,var)
    expect(len(matches) == 1, f"Expected one match for {var}")
    node = matches[0]

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
