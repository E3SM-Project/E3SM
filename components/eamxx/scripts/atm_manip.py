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
class AmbiguousName (Exception):
    pass
###############################################################################

###############################################################################
def num_nodes_with_name (root,name,recurse=True):
###############################################################################
    """
    Count nodes with certain name in an XML tree

    >>> xml = '''
    ... <root>
    ...     <a/>
    ...     <sub>
    ...         <a/>
    ...     </sub>
    ... </root>
    ... '''
    >>> import xml.etree.ElementTree as ET
    >>> tree = ET.fromstring(xml)
    >>> num_nodes_with_name(tree,'a',recurse=False)
    1
    >>> num_nodes_with_name(tree,'a',recurse=True)
    2
    """

    count = 0

    for elem in root:
        if elem.tag==name:
            count += 1
        if recurse:
            count += num_nodes_with_name(elem,name)

    return count

###############################################################################
def find_node (root,name,recurse=True):
###############################################################################
    """
    >>> xml = '''
    ... <root>
    ...     <sub>
    ...         <a>2</a>
    ...     </sub>
    ... </root>
    ... '''
    >>> import xml.etree.ElementTree as ET
    >>> tree = ET.fromstring(xml)
    >>> a,parents = find_node(tree,'a')
    >>> print(f"{','.join(p.tag for p in parents)}")
    root,sub
    >>> print(a.text)
    2
    >>> print(len(parents))
    2
    >>> print(f"{','.join(p.tag for p in parents)}")
    root,sub
    >>> a,parents = find_node(tree,'a',recurse=False)
    >>> print(a)
    None
    """

    if root.tag==name:
        return root, []

    for elem in root:
        if elem.tag==name:
            return elem, [root]
        if len(elem)>0 and recurse:
            found, parents = find_node(elem,name,recurse=True)
            if found is not None:
                return found, [root] + parents

    return None, []

###############################################################################
def get_xml_node(xml_root,name,allow_not_found=False):
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
    >>> ################ INVALID SYNTAX #######################
    SystemExit: ERROR: Invalid change format. Expected A[::B[...]=value, got' prop1->2'
    >>> get_xml_node(tree,'sub::::prop1')
    Traceback (most recent call last):
    SystemExit: ERROR: Invalid xml node name format. Expected A[::B[...], got' sub::::prop1'
      Did you put two '::' in a row?
    >>> ################ INVALID NAMESPACE #######################
    >>> get_xml_node(tree,'invalid::prop1')
    Traceback (most recent call last):
    SystemExit: ERROR: Error! XML entry invalid not found in section root
    >>> ################ AMBIGUOUS ENTRY #######################
    >>> get_xml_node(tree,'prop1')
    Traceback (most recent call last):
    atm_manip.AmbiguousName: ERROR: Error! Multiple XML entries with name prop1 found in section root
    >>> ################ VALID USAGE #######################
    >>> n,p = get_xml_node(tree,'::prop1')
    >>> print(n.text)
    one
    >>> print(len(p))
    1
    >>> print(p[0].tag)
    root
    >>> n,p = get_xml_node(tree,'prop2')
    >>> print(n.text)
    2
    >>> m,p = get_xml_node(tree,'prop2')
    >>> print([k for k in n.attrib.keys()])
    ['type', 'valid_values']
    >>> print(len(p))
    2
    >>> print(f"{','.join(e.tag for e in p)}")
    root,sub
    """

    selectors = name.split("::")

    # Allow :: at the beginning (as in '::A::b'), but do not allow multiple :: operators
    expect('' not in selectors[1:],
        f"Invalid xml node name format. Expected A[::B[...], got' {name}'\n"
        "  Did you put two '::' in a row?")

    # Regardless of whether we have namespaces or not, the first selector must be unique through the whole XML tree
    s = selectors[0]
    if s == '':
        # User started with ::
        node = xml_root
        parents = []
    else:
        expect (allow_not_found or num_nodes_with_name(xml_root,s,recurse=True)>0,
            f"Error! XML entry {s} not found in section {xml_root.tag}")
        expect (allow_not_found or num_nodes_with_name(xml_root,s,recurse=True)==1,
            f"Error! Multiple XML entries with name {s} found in section {xml_root.tag}",
            AmbiguousName)

        node, parents = find_node(xml_root,s,recurse=True)

    if node is None:
        expect (allow_not_found,
                f"Error! XML entry {s} not found in section {xml_root.tag}")
        return None, []

    # If user specified selectors via namespace, recurse over them
    for s in selectors[1:]:
        expect (allow_not_found or num_nodes_with_name(node,s,recurse=False)>0,
            f"Error! XML entry {s} not found in section {node.tag}")
        expect (allow_not_found or num_nodes_with_name(node,s,recurse=False)==1,
            f"Error! Multiple XML entries with name {s} found in section {node.tag}")

        node, parents = find_node(node,s,recurse=False)

    return node, parents

###############################################################################
def apply_change (node, new_value, append_this):
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
def parse_change (change):
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
def atm_config_chg_impl(xml_parent, xml_node, change, all_matches=False):
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
    >>> atm_config_chg_impl(None,tree,'prop1->2')
    Traceback (most recent call last):
    SystemExit: ERROR: Invalid change request 'prop1->2'. Valid formats are:
      - A[::B[...]=value
      - A[::B[...]+=value  (implies append for this change)
    >>> ################ INVALID TYPE #######################
    >>> atm_config_chg_impl(None,tree,'prop2=two')
    Traceback (most recent call last):
    ValueError: Could not use 'two' as type 'integer'
    >>> ################ INVALID VALUE #######################
    >>> atm_config_chg_impl(None,tree,'prop2=3')
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: Invalid value '3' for element 'prop2'. Value not in the valid list ('[1, 2]')
    >>> ################ VALID USAGE #######################
    >>> atm_config_chg_impl(None,tree,'::prop1=two')
    (True, True)
    >>> atm_config_chg_impl(None,tree,'::prop1=two')
    (True, False)
    >>> atm_config_chg_impl(None,tree,'sub::prop1=one')
    (True, True)
    >>> ################ TEST APPEND += #################
    >>> atm_config_chg_impl(None,tree,'a+=4')
    (True, True)
    >>> get_xml_node(None,tree,'a')[0].text
    '1,2,3, 4'
    >>> ################ ERROR, append to non-array and non-string
    >>> atm_config_chg_impl(None,tree,'c+=2')
    Traceback (most recent call last):
    SystemExit: ERROR: Error! Can only append with array and string types.
        - name: c
        - type: int
    >>> ################ Append to string ##################
    >>> atm_config_chg_impl(None,tree,'d+=two')
    (True, True)
    >>> get_xml_node(None,tree,'d')[0].text
    'onetwo'
    >>> ################ Append to array(string) ##################
    >>> atm_config_chg_impl(None,tree,'e+=two')
    (True, True)
    >>> get_xml_node(None,tree,'e')[0].text
    'one, two'
    """

    any_change = False
    node_found = False
    if all_matches and len(xml_node)>0:
        # We have to go through the whole tree, since we need to apply
        # the changes to all nodes matching the node name
        for elem in xml_node:
            found_here, changed_here = atm_config_chg_impl(xml_node,elem,change, all_matches)
            any_change |= changed_here
            node_found |= found_here
    else:
        node_name,new_value,append_this = parse_change(change)
        node, __ = get_xml_node(xml_parent,node_name,allow_not_found=all_matches)

        node_found = node is not None
        if node is not None:
            any_change = apply_change (node,new_value,append_this)


    return node_found, any_change

###############################################################################
def print_var_impl(node,parents,full,dtype,value,valid_values,print_style="invalid",indent=""):
###############################################################################

    expect (print_style in ["short","full"],
            f"Invalid print_style '{print_style}' for print_var_impl. Use 'full' or 'short'.")

    if print_style=="short":
        # Just the inner most name
        name = node.tag
    else:
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
def print_var(xml_root,var,full,dtype,value,valid_values,print_style="invalid",indent=""):
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
    >>> ################ Missing type data #######################
    >>> print_var(tree,'::prop1',False,True,False,False,"short")
    Traceback (most recent call last):
    SystemExit: ERROR: Error! Missing type information for prop1
    >>> print_var(tree,'prop2',True,False,False,False,"short")
    prop2
        value: 2
        type: integer
        valid values: ['1', '2']
    >>> print_var(tree,'prop2',False,True,False,False,"short")
    prop2: integer
    >>> print_var(tree,'prop2',False,False,True,False,"short")
    2
    >>> print_var(tree,'prop2',False,False,False,True,"short","    ")
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
        try:
            get_xml_node(xml_root,new_name)
            tokens.pop(0)
        except AmbiguousName:
            # new_name was either "" or an ambiguous name, and get_xml_node failed
            break

    # Get node, along with all its parents (which might be used for 'full' print style)
    node, parents = get_xml_node(xml_root,var)

    print_var_impl(node,parents,full,dtype,value,valid_values,print_style,indent)

###############################################################################
def print_all_vars(xml_root,xml_node,curr_namespace,full,dtype,value,valid_values,print_style,indent):
###############################################################################

    print (f"{indent}{xml_node.tag}")
    for c in xml_node:
        if len(c)>0:
            print_all_vars(xml_root,c,curr_namespace+c.tag+"::",full,dtype,value,valid_values,print_style,indent+"    ")
        else:
            print_var(xml_root,curr_namespace+c.tag,full,dtype,value,valid_values,print_style,indent+"    ")

###############################################################################
def atm_query_impl(xml_root,variables,listall=False,full=False,value=False, \
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
    >>> success = atm_query_impl(tree, vars, False,False,False,False,False)
        root::sub::prop2: 2
        root::prop1: one
    >>> success = atm_query_impl(tree, [], True,False,False,False,True)
        root
            prop1: <valid values not provided>
            sub
                prop1: <valid values not provided>
                prop2: ['1', '2']
    >>> success = atm_query_impl(tree,['prop1'],False,False,False,False,False,True)
        root::prop1: one
        sub::prop1: two
    """

    if listall:
        print_all_vars(xml_root,xml_root,"::",full,dtype,value,valid_values,"short","    ")
    elif grep:
        for regex in variables:
            var_re = re.compile(f'{regex}')
            if var_re.search(xml_root.tag):
                print_all_vars(xml_root,xml_root,"::",full,dtype,value,valid_values,"short","  ")
            else:
                for elem in xml_root:
                    if len(elem)>0:
                        atm_query_impl(elem,variables,listall,full,value,dtype,valid_values,grep)
                    else:
                        #  print (f"checking {elem.tag}")
                        if var_re.search(elem.tag):
                            node, parents = find_node(xml_root,elem.tag,recurse=False)
                            print_var_impl(node,parents,full,dtype,value,valid_values,"full","    ")
    else:
        for var in variables:
            print_var(xml_root,var,full,dtype,value,valid_values,"full","    ")

    return True
