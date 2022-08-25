"""
Retrieve nodes from EAMxx XML config file.
"""

import sys, os
import xml.etree.ElementTree as ET

# Add path to cime_config folder
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), "cime_config"))
from eamxx_buildnml_impl import check_value

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

    for elem in root:
        if elem.tag==name:
            return elem, [root]
        if len(elem)>0 and recurse:
            found, parents = find_node(elem,name,recurse=True)
            if found is not None:
                return found, [root] + parents

    return None, []

###############################################################################
def expect(condition, error_msg, exc_type=SystemExit, error_prefix="ERROR:"):
###############################################################################
    """
    Similar to assert except doesn't generate an ugly stacktrace. Useful for
    checking user error, not programming error.
    """
    if not condition:
        msg = error_prefix + " " + error_msg
        raise exc_type(msg)

###############################################################################
def get_xml_node(xml_root,name):
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
        "Invalid xml node name format. Expected A[::B[...], got' {}'\n".format(name) +
        "  Did you put two '::' in a row?")

    # Regardless of whether we have namespaces or not, the first selector must be unique through the whole XML tree
    s = selectors[0]
    if s == '':
        # User started with ::
        node = xml_root
        parents = []
    else:
        expect (num_nodes_with_name(xml_root,s,recurse=True)>0,
            "Error! XML entry {} not found in section {}".format(s,xml_root.tag))
        expect (num_nodes_with_name(xml_root,s,recurse=True)==1,
            "Error! Multiple XML entries with name {} found in section {}"
            .format(s,xml_root.tag), AmbiguousName)

        node, parents = find_node(xml_root,s,recurse=True)

    # If user specified selectors via namespace, recurse over them
    for s in selectors[1:]:
        expect (num_nodes_with_name(node,s,recurse=False)>0,
            "Error! XML entry {} not found in section {}".format(s,node.tag))
        expect (num_nodes_with_name(node,s,recurse=False)==1,
            "Error! Multiple XML entries with name {} found in section {}"
            .format(s,node.tag))

        node, parents = find_node(node,s,recurse=False)

    return node, parents
