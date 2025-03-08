#!/usr/bin/env python3

"""
Used by buildnml. See buildnml for documetation.
"""

import os, sys, re, pwd, grp, stat, getpass
from collections import OrderedDict

import xml.etree.ElementTree as ET
import xml.dom.minidom as md

# Add path to scream libs
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))

# SCREAM imports
from eamxx_buildnml_impl import get_valid_selectors, get_child, refine_type, \
        resolve_all_inheritances, gen_atm_proc_group, check_all_values, find_node
from atm_manip import apply_atm_procs_list_changes_from_buffer, apply_non_atm_procs_list_changes_from_buffer

from utils import ensure_yaml # pylint: disable=no-name-in-module
ensure_yaml()
import yaml
from yaml_utils import Bools,Ints,Floats,Strings,array_representer,array_constructor

_CIMEROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..","..","..","cime")
sys.path.append(os.path.join(_CIMEROOT, "CIME", "Tools"))

# Cime imports
from standard_script_setup import * # pylint: disable=wildcard-import
from CIME.utils import expect, safe_copy, SharedArea

logger = logging.getLogger(__name__) # pylint: disable=undefined-variable

CIME_VAR_RE = re.compile(r'[$][{](\w+)[}]')

# These are special attributes used by buildnml and other scripts to
# perform some checks. In particular:
#  - type: allows to verify compatibility (e.g., can't assing 3.2 to an integer)
#  - valid_values: allows to specify a set of valid values
#  - locked: if set to true, the parameter cannot be modified (via atmchange)
#  - constraints: allows to specify constraints on values. Valid constraints
#    are lt, le, ne, gt, ge, and mod. Multiple constrained are separated by ';'.
#    Examples:
#      - constraints="ge 0; lt 4" means the value V must satisfy V>=0 && V<4.
#      - constraints="mod 2 eq 0" means the value V must be a multiple of 2.
METADATA_ATTRIBS = ("type", "valid_values", "locked", "constraints", "inherit", "doc", "append")

###############################################################################
def do_cime_vars(entry, case, refine=False, extra=None):
###############################################################################
    """
    For a parameter value, process any references to case values ("${CASEVAR}")

    >>> from eamxx_buildnml_impl import MockCase
    >>> case = MockCase({'foo':1, 'bar':2, 'baz_bar':'blah', 'baz2_bar2':4})
    >>> do_cime_vars('${foo}', case)
    '1'
    >>> do_cime_vars('hi ${foo} there', case)
    'hi 1 there'
    >>> do_cime_vars('hi ${foo}${bar} there', case)
    'hi 12 there'
    >>> do_cime_vars('hi ${baz_bar} there', case)
    'hi blah there'
    >>> do_cime_vars('hi ${baz2_bar2} there', case)
    'hi 4 there'
    >>> do_cime_vars('hi ${invalid} there', case)
    Traceback (most recent call last):
      ...
    CIME.utils.CIMEError: ERROR: Cannot resolve XML entry 'hi ${invalid} there', CIME has no value for 'invalid'
    >>> d = { 'foo' : '${foo}',
    ...      'subdict' : { 'bar' : 'foo', 'baz' : '${foo}' } }
    >>> do_cime_vars(d, case)
    {'foo': '1', 'subdict': {'bar': 'foo', 'baz': '1'}}
    """
    if type(entry) is dict:
        changed_keys = []
        for k,v in entry.items():
            # Apply replacement to the key.
            k1 = do_cime_vars(k,case,refine,extra)
            if (k1 != k): changed_keys.append((k,k1))
            # Apply replacement to the value.
            entry[k] = do_cime_vars(v,case,refine,extra)
        # Handle changed keys outside the original loop since keys can't change
        # in an items loop.
        for e in changed_keys:
            entry[e[1]] = entry[e[0]]
            del entry[e[0]]
    elif type(entry) is str:
        m = CIME_VAR_RE.search(entry)
        while m:
            cime_var = m.groups()[0]
            if type(extra) is dict and cime_var in extra:
                value = extra[cime_var]
            else:
                value = case.get_value(cime_var)
            expect(value is not None,
                   "Cannot resolve XML entry '{}', CIME has no value for '{}'".format(entry, cime_var))
            entry = entry.replace("${{{}}}".format(cime_var), str(value))
            m = CIME_VAR_RE.search(entry)

        if refine:
            entry = refine_type(entry)

    return entry

###############################################################################
def perform_consistency_checks(case, xml):
###############################################################################
    """
    There may be separate parts of the xml that must satisfy some consistency
    Here, we run any such check, so we can catch errors before submit time

    >>> from eamxx_buildnml_impl import MockCase
    >>> xml_str = '''
    ... <params>
    ...   <rrtmgp>
    ...     <rad_frequency type="integer">3</rad_frequency>
    ...   </rrtmgp>
    ... </params>
    ... '''
    >>> import xml.etree.ElementTree as ET
    >>> xml = ET.fromstring(xml_str)
    >>> case = MockCase({'ATM_NCPL':'24', 'REST_N':24, 'REST_OPTION':'nsteps'})
    >>> perform_consistency_checks(case,xml)
    >>> case = MockCase({'ATM_NCPL':'24', 'REST_N':2, 'REST_OPTION':'nsteps'})
    >>> perform_consistency_checks(case,xml)
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: rrtmgp::rad_frequency (3 steps) incompatible with restart frequency (2 steps).
     Please, ensure restart happens on a step when rad is ON
    >>> case = MockCase({'ATM_NCPL':'24', 'REST_N':10800, 'REST_OPTION':'nseconds'})
    >>> perform_consistency_checks(case,xml)
    >>> case = MockCase({'ATM_NCPL':'24', 'REST_N':7200, 'REST_OPTION':'nseconds'})
    >>> perform_consistency_checks(case,xml)
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: rrtmgp::rad_frequency incompatible with restart frequency.
     Please, ensure restart happens on a step when rad is ON
      rest_tstep: 7200
      rad_testep: 10800.0
    >>> case = MockCase({'ATM_NCPL':'24', 'REST_N':180, 'REST_OPTION':'nminutes'})
    >>> perform_consistency_checks(case,xml)
    >>> case = MockCase({'ATM_NCPL':'24', 'REST_N':120, 'REST_OPTION':'nminutes'})
    >>> perform_consistency_checks(case,xml)
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: rrtmgp::rad_frequency incompatible with restart frequency.
     Please, ensure restart happens on a step when rad is ON
      rest_tstep: 7200
      rad_testep: 10800.0
    >>> case = MockCase({'ATM_NCPL':'24', 'REST_N':6, 'REST_OPTION':'nhours'})
    >>> perform_consistency_checks(case,xml)
    >>> case = MockCase({'ATM_NCPL':'24', 'REST_N':8, 'REST_OPTION':'nhours'})
    >>> perform_consistency_checks(case,xml)
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: rrtmgp::rad_frequency incompatible with restart frequency.
     Please, ensure restart happens on a step when rad is ON
      rest_tstep: 28800
      rad_testep: 10800.0
    >>> case = MockCase({'ATM_NCPL':'12', 'REST_N':2, 'REST_OPTION':'ndays'})
    >>> perform_consistency_checks(case,xml)
    >>> case = MockCase({'ATM_NCPL':'10', 'REST_N':2, 'REST_OPTION':'ndays'})
    >>> perform_consistency_checks(case,xml)
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: rrtmgp::rad_frequency incompatible with restart frequency.
     Please, ensure restart happens on a step when rad is ON
     For daily (or less frequent) restart, rad_frequency must divide ATM_NCPL
    """

    # RRTMGP can be supercycled. Restarts cannot fall in the middle
    # of a rad superstep
    rrtmgp = find_node(xml,"rrtmgp")
    rest_opt = case.get_value("REST_OPTION")
    is_test = case.get_value("TEST")
    if rrtmgp is not None and rest_opt is not None and rest_opt not in ["never","none"]:
        rest_n = int(case.get_value("REST_N"))
        rad_freq = int(find_node(rrtmgp,"rad_frequency").text)
        atm_ncpl = int(case.get_value("ATM_NCPL"))
        atm_tstep = 86400 / atm_ncpl
        rad_tstep = atm_tstep * rad_freq


        # Some tests (ERS) make late (run-phase) changes, so we cannot validate restart
        # settings here.
        if rad_freq==1 or is_test:
            pass
        elif rest_opt in ["nsteps", "nstep"]:
            expect (rest_n % rad_freq == 0,
                    f"rrtmgp::rad_frequency ({rad_freq} steps) incompatible with "
                    f"restart frequency ({rest_n} steps).\n"
                    " Please, ensure restart happens on a step when rad is ON")
        elif rest_opt in ["nseconds", "nsecond", "nminutes", "nminute", "nhours", "nhour"]:
            if rest_opt in ["nseconds", "nsecond"]:
                factor = 1
            elif rest_opt in ["nminutes", "nminute"]:
                factor = 60
            else:
                factor = 3600

            rest_tstep = factor*rest_n
            expect (rest_tstep % rad_tstep == 0,
                    "rrtmgp::rad_frequency incompatible with restart frequency.\n"
                    " Please, ensure restart happens on a step when rad is ON\n"
                    f"  rest_tstep: {rest_tstep}\n"
                    f"  rad_testep: {rad_tstep}")

        else:
            # for "very infrequent" restarts, we request rad_freq to divide atm_ncpl
            expect (atm_ncpl % rad_freq ==0,
                    "rrtmgp::rad_frequency incompatible with restart frequency.\n"
                    " Please, ensure restart happens on a step when rad is ON\n"
                    " For daily (or less frequent) restart, rad_frequency must divide ATM_NCPL")

###############################################################################
def ordered_dump(data, item, Dumper=yaml.SafeDumper, **kwds):
###############################################################################
    """
    Copied from: https://stackoverflow.com/a/21912744
    Added ability to pass filename
    """

    class OrderedDumper(Dumper):
        pass
    def _dict_representer(dumper, data):
        return dumper.represent_mapping(
            yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
            data.items())
    OrderedDumper.add_representer(OrderedDict, _dict_representer)

    # These allow to dump arrays with a tag specifying the type
    OrderedDumper.add_representer(Bools,    array_representer)
    OrderedDumper.add_representer(Ints,     array_representer)
    OrderedDumper.add_representer(Floats,   array_representer)
    OrderedDumper.add_representer(Strings,  array_representer)

    if isinstance(item, str) and item.endswith(".yaml"):
        # Item is a filepath
        with open(item, "w") as fd:
            return yaml.dump(data, fd, OrderedDumper, **kwds)
    else:
        return yaml.dump(data, item, OrderedDumper, **kwds)

###############################################################################
def evaluate_selectors(element, case, ez_selectors):
###############################################################################
    """
    Evaluate and remove selectors from the unprocessed XML nml file in the repo.

    Elements with selectors are removed. If the selector evaulates to True, then
    the corresponding element text becomes the new text value of the original
    (no-selectors, AKA default) element.

    The metadata attributes are kept, to allow checks during calls to atmchange

    >>> from eamxx_buildnml_impl import MockCase
    >>> case = MockCase({'ATM_GRID':'ne4ne4', 'SCREAM_CMAKE_OPTIONS':'FOO=ON SCREAM_NUM_VERTICAL_LEV 128 BAR=OFF'})
    >>> xml_sel_good = '''
    ... <selectors_xml>
    ...   <selectors>
    ...     <selector name="grid" case_env="ATM_GRID"/>
    ...     <selector name="nlev" case_env="SCREAM_CMAKE_OPTIONS" regex=".*SCREAM_NUM_VERTICAL_LEV ([0-9]+).*"/>
    ...     <selector name="foo"  case_env="SCREAM_CMAKE_OPTIONS" regex=".*FOO=(ON|OFF).*"/>
    ...     <selector name="bar"  case_env="SCREAM_CMAKE_OPTIONS" regex=".*BAR=(ON|OFF).*"/>
    ...   </selectors>
    ... </selectors_xml>
    ... '''
    >>> xml_good = '''
    ... <namelist_defaults>
    ...   <var1>zero</var1>
    ...   <var1 grid="ne4ne4">one</var1>
    ...   <var2>zero</var2>
    ...   <var2 nlev="128">one</var2>
    ...   <var3>foo_off</var3>
    ...   <var3 foo="ON">foo_on</var3>
    ...   <var4>bar_off</var4>
    ...   <var4 bar="ON">bar_on</var4>
    ... </namelist_defaults>
    ... '''
    >>> import xml.etree.ElementTree as ET
    >>> good = ET.fromstring(xml_good)
    >>> selectors_good = get_valid_selectors(ET.fromstring(xml_sel_good))
    >>> ############## ALL GOOD #####################
    >>> evaluate_selectors(good,case,selectors_good)
    >>> get_child(good,'var1').text=="one"
    True
    >>> get_child(good,'var2').text=="one"
    True
    >>> get_child(good,'var3').text=="foo_on"
    True
    >>> get_child(good,'var4').text=="bar_off"
    True
    >>> ############## BAD SELECTOR DEFINITION #####################
    >>> xml_sel_bad1 = '''
    ... <selectors_xml>
    ...   <selectors>
    ...     <selector name="grid" case_env="BADENV"/>
    ...   </selectors>
    ... </selectors_xml>
    ... '''
    >>> selectors_bad1 = get_valid_selectors(ET.fromstring(xml_sel_bad1))
    >>> good = ET.fromstring(xml_good)
    >>> evaluate_selectors(good,case,selectors_bad1)
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: Bad easy selector 'grid' definition. Relies on unknown case value 'BADENV'
    >>> ############## BAD SELECTOR DEFINITION #####################
    >>> xml_sel_bad2 = '''
    ... <selectors_xml>
    ...   <selectors>
    ...     <selector name="grid" case_env="ATM_GRID" regex=".*"/>
    ...   </selectors>
    ... </selectors_xml>
    ... '''
    >>> selectors_bad2 = get_valid_selectors(ET.fromstring(xml_sel_bad2))
    >>> good = ET.fromstring(xml_good)
    >>> evaluate_selectors(good,case,selectors_bad2)
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: Selector 'grid' has invalid custom regex '.*' which does not capture exactly 1 group
    >>> ############## BAD SELECTOR NAME #####################
    >>> xml_bad1 = '''
    ... <namelist_defaults>
    ...   <var1>zero</var1>
    ...   <var1 my_grid="ne4ne4">one</var1>
    ... </namelist_defaults>
    ... '''
    >>> bad1 = ET.fromstring(xml_bad1)
    >>> evaluate_selectors(bad1,case,selectors_good)
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: Bad selector 'my_grid' for child 'var1'. 'my_grid' is not a valid case value or easy selector
    >>> ############## BAD DEFAULTS ORDERING #####################
    >>> xml_bad2 = '''
    ... <namelist_defaults>
    ...   <var1 grid="ne4ne4">one</var1>
    ...   <var1>zero</var1>
    ... </namelist_defaults>
    ... '''
    >>> bad2 = ET.fromstring(xml_bad2)
    >>> evaluate_selectors(bad2,case,selectors_good)
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: child 'var1' element without selectors occurred after other parameter elements for this parameter
    >>> ############## MULTIPLE MATCHES #####################
    >>> xml_bad3 = '''
    ... <namelist_defaults>
    ...   <var1>zero</var1>
    ...   <var1>one</var1>
    ... </namelist_defaults>
    ... '''
    >>> bad3 = ET.fromstring(xml_bad3)
    >>> evaluate_selectors(bad3,case,selectors_good)
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: child 'var1' element without selectors occurred after other parameter elements for this parameter
    """

    selected_child = {} # elem_name -> evaluated XML element
    children_to_remove = []
    child_base_value = {} # map elme name to values to be appended to if append=="base"
    child_type  = {} # map elme name to its type (since only first entry may have type specified)
    for child in element:
        # Note: in our system, an XML element is either a "node" (has children)
        # or a "leaf" (has a value).
        has_children = len(child) > 0
        if has_children:
            evaluate_selectors(child, case, ez_selectors)
        else:
            child_name = child.tag
            child.text = None if child.text is None else child.text.strip(' \n')
            child_val = child.text
            selectors = child.attrib

            if child_name not in child_type:
                child_type[child_name] = selectors["type"] if "type" in selectors.keys() else "unset"

            is_array = child_type[child_name].startswith("array")
            expect (is_array or "append" not in selectors.keys(),
                    "The 'append' metadata attribute is only supported for entries of array type\n"
                    f" param name: {child_name}\n"
                    f" param type: {child_type[child_name]}")

            append = selectors["append"] if "append" in selectors.keys() else "no"
            expect (append in ["no","base","last"],
                    "Unrecognized value for 'append' attribute\n" +
                    f"  param name  : {child_name}\n" +
                    f"  append value: {append}\n" +
                     "  valid values: base, last\n")
            if selectors:
                all_match = True
                had_case_selectors = False
                for k, v in selectors.items():
                    # Metadata attributes are used only when it's time to generate the input files
                    if k in METADATA_ATTRIBS:
                        if k=="type" and child_name in selected_child.keys():
                            if "type" in selected_child[child_name].attrib:
                                expect (v==selected_child[child_name].attrib["type"],
                                        f"The 'type' attribute of {child_name} is not consistent across different selectors")
                        continue

                    had_case_selectors = True
                    val_re = re.compile(v)

                    if k in ez_selectors:
                        ez_env, ez_regex = ez_selectors[k]
                        case_val = case.get_value(ez_env)
                        expect(case_val is not None,
                              "Bad easy selector '{}' definition. Relies on unknown case value '{}'".format(k, ez_env))

                        ez_regex_re = re.compile(ez_regex)
                        m = ez_regex_re.match(case_val)
                        if m:
                            groups = m.groups()
                            expect(len(groups) == 1,
                                    "Selector '{}' has invalid custom regex '{}' which does not capture exactly 1 group".format(k, ez_regex))
                            val = groups[0]
                        else:
                            # If the regex doesn't even match the case val, then we consider
                            # string below should ensure the selector will never match.
                            val = None

                    else:
                        val = case.get_value(k)
                        expect(val is not None,
                               "Bad selector '{0}' for child '{1}'. '{0}' is not a valid case value or easy selector".format(k, child_name))

                    if val is None or val_re.match(val) is None:
                        all_match = False
                        children_to_remove.append(child)
                        break

                if all_match:
                    if child_name in selected_child.keys():
                        orig_child = selected_child[child_name]
                        if append=="base":
                            orig_child.text = child_base_value[child_name] + "," + child.text
                        elif append=="last":
                            orig_child.text = orig_child.text + "," + child.text
                        else:
                            orig_child.text = child.text
                        children_to_remove.append(child)

                    else:
                        # If all selectors were the METADATA_ATTRIB ones, then this is the "base" value
                        if not had_case_selectors:
                            child_base_value[child_name] = child.text
                        selected_child[child_name] = child
                        # Make a copy of selectors.keys(), since selectors=child.attrib,
                        # and we might delete an entry, causing the error
                        #    RuntimeError: dictionary changed size during iteration

            else:
                expect(child_name not in selected_child,
                       "child '{}' element without selectors occurred after other parameter elements for this parameter".format(child_name))
                child_base_value[child_name] = child.text
                selected_child[child_name] = child
                child.text = do_cime_vars(child_val, case)

    for child_to_remove in children_to_remove:
        element.remove(child_to_remove)

###############################################################################
def expand_cime_vars(element, case):
###############################################################################
    """
    Expand all CIME variables inside an XML node text
    """

    for child in element:
        # Note: in our system, an XML element is either a "node" (has children)
        # or a "leaf" (has a value).
        has_children = len(child) > 0
        if has_children:
            expand_cime_vars(child, case)
        else:
            child.text = do_cime_vars(child.text, case)

###############################################################################
def write_pretty_xml(filepath, xml):
###############################################################################
    with open(filepath, "w") as fd:
        # dom has better pretty printing than ET in older python versions < 3.9
        dom = md.parseString(ET.tostring(xml, encoding="unicode"))
        pretty_xml = dom.toprettyxml(indent="  ")
        pretty_xml = os.linesep.join([s for s in pretty_xml.splitlines()
                                      if s.strip()])
        fd.write(pretty_xml)

###############################################################################
def _create_raw_xml_file_impl(case, xml, filepath=None):
###############################################################################
    """
    On input, xml contains the parsed content of namelist_defaults_eamxx.xml.
    On output, it contains the input parameters for this case.

    >>> from eamxx_buildnml_impl import MockCase
    >>> ############## BASIC TEST #####################
    >>> case = MockCase({'foo':1, 'bar':2, 'baz_bar':'blah', 'baz2_bar2':4})
    >>> xml = '''
    ... <namelist_defaults>
    ...     <selectors/>
    ...     <generated_files/>
    ...     <atmosphere_processes_defaults>
    ...         <atm_procs_list type="array(string)">P1,P2</atm_procs_list>
    ...         <atm_proc_base>
    ...             <prop1>zero</prop1>
    ...         </atm_proc_base>
    ...         <atm_proc_group inherit="atm_proc_base">
    ...             <atm_procs_list>NONE</atm_procs_list>
    ...             <prop2>one</prop2>
    ...         </atm_proc_group>
    ...         <P1 inherit="atm_proc_base">
    ...             <prop1>two</prop1>
    ...         </P1>
    ...         <P2 inherit="atm_proc_base">
    ...         </P2>
    ...     </atmosphere_processes_defaults>
    ... </namelist_defaults>
    ... '''
    >>> import xml.etree.ElementTree as ET
    >>> defaults = ET.fromstring(xml)
    >>> generated = _create_raw_xml_file_impl(case,defaults)
    >>> d = convert_to_dict(get_child(generated,'atmosphere_processes'))
    >>> import pprint
    >>> pp = pprint.PrettyPrinter(indent=4)
    >>> pp.pprint(d)
    OrderedDict([   ('atm_procs_list', 'P1,P2'),
                    ('prop2', 'one'),
                    ('prop1', 'zero'),
                    ('P1', OrderedDict([('prop1', 'two')])),
                    ('P2', OrderedDict([('prop1', 'zero')]))])
    >>> ############## INHERIT+CHILD SELECTOR #####################
    >>> case = MockCase({'ATM_GRID':'ne4ne4'})
    >>> xml = '''
    ... <namelist_defaults>
    ...     <selectors>
    ...       <selector name="grid" case_env="ATM_GRID"/>
    ...     </selectors>
    ...     <generated_files/>
    ...     <atmosphere_processes_defaults>
    ...         <atm_procs_list type="array(string)">P1,P2</atm_procs_list>
    ...         <atm_proc_base>
    ...             <prop1>zero</prop1>
    ...         </atm_proc_base>
    ...         <atm_proc_group inherit="atm_proc_base">
    ...             <atm_procs_list>NONE</atm_procs_list>
    ...             <prop2>one</prop2>
    ...         </atm_proc_group>
    ...         <P1 inherit="atm_proc_base">
    ...             <prop1>two</prop1>
    ...             <prop1 grid='ne4ne4'>two_selected</prop1>
    ...         </P1>
    ...         <P2 inherit="atm_proc_base">
    ...         </P2>
    ...     </atmosphere_processes_defaults>
    ... </namelist_defaults>
    ... '''
    >>> import xml.etree.ElementTree as ET
    >>> defaults = ET.fromstring(xml)
    >>> generated = _create_raw_xml_file_impl(case,defaults)
    >>> d = convert_to_dict(get_child(generated,'atmosphere_processes'))
    >>> import pprint
    >>> pp = pprint.PrettyPrinter(indent=4)
    >>> pp.pprint(d)
    OrderedDict([   ('atm_procs_list', 'P1,P2'),
                    ('prop2', 'one'),
                    ('prop1', 'zero'),
                    ('P1', OrderedDict([('prop1', 'two_selected')])),
                    ('P2', OrderedDict([('prop1', 'zero')]))])
    >>> ############## INHERIT+PARENT SELECTOR #####################
    >>> case = MockCase({'ATM_GRID':'ne4ne4'})
    >>> xml = '''
    ... <namelist_defaults>
    ...     <selectors>
    ...       <selector name='grid' case_env='ATM_GRID'/>
    ...     </selectors>
    ...     <generated_files/>
    ...     <atmosphere_processes_defaults>
    ...       <atm_procs_list type="array(string)">P1,P2</atm_procs_list>
    ...       <atm_proc_base>
    ...         <number_of_subcycles constraints='gt 0'>1</number_of_subcycles>
    ...         <enable_precondition_checks type='logical'>true</enable_precondition_checks>
    ...         <enable_postcondition_checks type='logical'>true</enable_postcondition_checks>
    ...       </atm_proc_base>
    ...       <physics_proc_base inherit='atm_proc_base'>
    ...         <Grid>Physics GLL</Grid>
    ...         <Grid grid='ne4ne4'>Physics PG2</Grid>
    ...       </physics_proc_base>
    ...       <atm_proc_group inherit="atm_proc_base">
    ...         <atm_procs_list>NONE</atm_procs_list>
    ...         <prop2>one</prop2>
    ...       </atm_proc_group>
    ...       <P1 inherit='physics_proc_base'>
    ...         <prop1>hi</prop1>
    ...       </P1>
    ...       <P2 inherit='atm_proc_base'>
    ...         <prop1>there</prop1>
    ...       </P2>
    ...     </atmosphere_processes_defaults>
    ... </namelist_defaults>
    ... '''
    >>> import xml.etree.ElementTree as ET
    >>> defaults = ET.fromstring(xml)
    >>> generated = _create_raw_xml_file_impl(case,defaults)
    >>> d = convert_to_dict(get_child(generated,'atmosphere_processes'))
    >>> import pprint
    >>> pp = pprint.PrettyPrinter(indent=4)
    >>> pp.pprint(d)
    OrderedDict([   ('atm_procs_list', 'P1,P2'),
                    ('prop2', 'one'),
                    ('number_of_subcycles', 1),
                    ('enable_precondition_checks', True),
                    ('enable_postcondition_checks', True),
                    (   'P1',
                        OrderedDict([   ('prop1', 'hi'),
                                        ('Grid', 'Physics PG2'),
                                        ('number_of_subcycles', 1),
                                        ('enable_precondition_checks', True),
                                        ('enable_postcondition_checks', True)])),
                    (   'P2',
                        OrderedDict([   ('prop1', 'there'),
                                        ('number_of_subcycles', 1),
                                        ('enable_precondition_checks', True),
                                        ('enable_postcondition_checks', True)]))])

    """

    # 0. Remove internal sections, that are not to appear in the
    #    processed xml file.
    #    Note: get_valid_selectors also removes the section
    get_child(xml,"generated_files",remove=True)
    selectors = get_valid_selectors(xml)

    # 1. Evaluate all selectors
    try:
        evaluate_selectors(xml, case, selectors)

        # 2. Apply all changes in the SCREAM_ATMCHANGE_BUFFER that may alter
        #    which atm processes are used
        apply_atm_procs_list_changes_from_buffer (case,xml)

        # 3. Resolve all inheritances
        resolve_all_inheritances(xml)

        # 4. Expand any CIME var that appears inside XML nodes text
        expand_cime_vars(xml,case)

        # 5. Grab the atmosphere_processes macro list, with all the defaults
        atm_procs_defaults = get_child(xml,"atmosphere_processes_defaults",remove=True)

        # 6. Get atm procs list
        atm_procs_list = get_child(atm_procs_defaults,"atm_procs_list",remove=True)

        # 7. Form the nested list of atm procs needed, append to atmosphere_driver section
        atm_procs = gen_atm_proc_group(atm_procs_list.text, atm_procs_defaults)
        atm_procs.tag = "atmosphere_processes"
        xml.append(atm_procs)

        # 8. Apply all changes in the SCREAM_ATMCHANGE_BUFFER that do not alter
        #    which atm processes are used
        apply_non_atm_procs_list_changes_from_buffer (case,xml)
    except BaseException as e:
        if filepath is not None:
            dbg_xml_path = filepath.replace(".xml", ".dbg.xml")
            write_pretty_xml(dbg_xml_path, xml)
            print(f"Error during XML creation, writing {dbg_xml_path}")

        raise e

    perform_consistency_checks (case, xml)

    return xml

###############################################################################
def create_raw_xml_file(case, caseroot):
###############################################################################
    """
    Create the $case/namelist_scream.xml file. This file is intended to be
    modified by users via the atmchange script if they want,
    to make tweaks to input files (yaml and/or nml).
    Note: users calls to atmchange do two things: 1) they add the change
          to the SCREAM_ATMCHANGE_BUFFER case variable, and 2) they
          call this function, which regenerates the scream xml file from
          the defaults, applyin all buffered changes.
    """
    raw_xml_file = os.path.join(caseroot, "namelist_scream.xml")
    if os.path.exists(raw_xml_file) and case.get_value("SCREAM_HACK_XML"):
        print("{} already exists and SCREAM_HACK_XML is on, will not overwrite. Remove to regenerate".format(raw_xml_file))

    else:
        print("Regenerating {}. Manual edits will be lost.".format(raw_xml_file))

        src = os.path.join(case.get_value("SRCROOT"), "components/eamxx/cime_config/namelist_defaults_eamxx.xml")

        # Some atmchanges will require structural changes to the XML file and must
        # be processed early by treating them as if they were made to the defaults file.
        with open(src, "r") as fd:
            defaults = ET.parse(fd).getroot()
            raw_xml = _create_raw_xml_file_impl(case, defaults, filepath=raw_xml_file)

        check_all_values(raw_xml)

        write_pretty_xml(raw_xml_file, raw_xml)

###############################################################################
def convert_to_dict(element):
###############################################################################
    """
    Convert an XML element to a dictonary where the tags are the keys
    >>> xml = '''
    ... <my_root>
    ...   <my__int>1</my__int>
    ...   <my__list>
    ...     <my_ints type="array(integer)">2,3</my_ints>
    ...     <my_strings type="array(string)">two,three</my_strings>
    ...   </my__list>
    ... </my_root>
    ... '''
    >>> import xml.etree.ElementTree as ET
    >>> import pprint
    >>> pp = pprint.PrettyPrinter(indent=4)
    >>> root = ET.fromstring(xml)
    >>> d = convert_to_dict(root)
    >>> pp.pprint(d)
    OrderedDict([   ('my int', 1),
                    (   'my list',
                        OrderedDict([   ('my_ints', [2, 3]),
                                        ('my_strings', ['two', 'three'])]))])
    """
    result = OrderedDict()
    for child in element:
        child_name = child.tag.replace("__", " ")

        has_children = len(child) > 0
        if has_children:
            result[child_name] = convert_to_dict(child)
        else:
            child_val = child.text
            force_type = None if "type" not in child.attrib.keys() else child.attrib["type"]
            result[child_name] = refine_type(child_val,force_type=force_type)

    return result

###############################################################################
def _dump_to_nml_impl(dict_contents):
###############################################################################
    """
    Given a dictionary representing namelists, dumps content to a string
    that can be directly written to file

    >>> ############# NON-NESTED DICT ##############
    >>> good1 = { 'bool': True, 'int': 1, 'real': 2.0, 'string': 's' }
    >>> print(_dump_to_nml_impl(good1))
    bool = .true.
    int = 1
    real = 2.0
    string = 's'
    <BLANKLINE>
    >>> ############# NESTED DICT ##############
    >>> good2 = {
    ...     'group1': {
    ...         'bool': True,
    ...         'int': 1
    ...     },
    ...     'group2': {
    ...         'real': 2.0,
    ...         'string': 's'
    ...     }
    ... }
    >>> print(_dump_to_nml_impl(good2))
    &group1
    bool = .true.
    int = 1
    /
    &group2
    real = 2.0
    string = 's'
    /
    <BLANKLINE>
    >>> ######### MIXED NESTED-NONNESTED DICT ##########
    >>> good2 = {
    ...     'group1': {
    ...         'bool': True,
    ...         'int': 1
    ...     },
    ...     'real': 2.0,
    ...     'string': 's'
    ... }
    >>> print(_dump_to_nml_impl(good2))
    Traceback (most recent call last):
    CIME.utils.CIMEError: ERROR: Error! _dump_to_nml_impl cannot mix nested and non-nested dicts.
    """

    result = ""
    nested = False
    not_nested = False
    for k, v in dict_contents.items():
        if isinstance(v, dict):
            result += "&{}\n".format(k)
            result += _dump_to_nml_impl(v)
            result += "/\n"
            nested = True

        elif isinstance(v, list):
            result += "{}(:) = {}\n".format(k, ", ".join([str(item) for item in v]))
            not_nested = True
        elif isinstance(v, bool):
            result += "{} = {}\n".format(k, ".true." if v else ".false.")
            not_nested = True
        elif isinstance(v, str):
            result += "{} = '{}'\n".format(k, v)
            not_nested = True
        else:
            result += "{} = {}\n".format(k, v)
            not_nested = True

    expect (not (nested and not_nested),
            "Error! _dump_to_nml_impl cannot mix nested and non-nested dicts.")

    return result

###############################################################################
def dump_to_nml(dict_contents, fd):
###############################################################################
    """
    Convert a dictionary to namelist file
    """
    fd.write(
"""!---------------------------------------------------------------
! Do NOT modify this file. It is generated by scream/buildnml using
! the data from $case/namelist_scream.xml. If you want to make some
! local changes, you can edit this XML file or use atmchange.
!---------------------------------------------------------------
""")

    fd.write(_dump_to_nml_impl(dict_contents))

###############################################################################
def create_input_files(caseroot, screamroot, rundir):
###############################################################################
    """
    Based on $case/namelist_scream.xml, create the actual input files.
    """
    raw_xml_file = os.path.join(caseroot, "namelist_scream.xml")
    with open(raw_xml_file, "r") as fd:
        tree = ET.parse(fd)
        raw_xml = tree.getroot()

    def_xml_file = os.path.join(screamroot, "cime_config/namelist_defaults_eamxx.xml")
    with open(def_xml_file, "r") as fd:
        tree = ET.parse(fd)
        generated_files = get_child(tree.getroot(),"generated_files")

    result = {}
    for file_ in generated_files:
        expect("name" in file_.attrib, "file element missing required 'name' attribute")
        expect("format" in file_.attrib, "file element missing required 'format' attribute")
        filename    = os.path.join(rundir, file_.attrib["name"])
        file_format = file_.attrib["format"]

        xml_sections = get_child(file_,"sections").text.split(",")

        dict_contents = {}
        for section_name in xml_sections:
            section = get_child(raw_xml,section_name)
            dict_contents[section_name.replace("__"," ")] = convert_to_dict(section)

        if file_format == "yaml":
            with open(filename, "w") as fd:
                fd.write(
"""################################################################
# Do NOT modify this file. It is generated by scream/buildnml using
# the data from $case/namelist_scream.xml. If you want to make some
# local changes, you can edit this XML file or use atmchange.
################################################################
""")
                ordered_dump(dict_contents, fd)

        elif file_format == "nml":
            with open(filename, "w") as fd:
                dump_to_nml(dict_contents, fd)

        else:
            expect(False, "Unknown input file format '{}'".format(file_format))

        result[file_.attrib["name"]] = dict_contents

    return result

###############################################################################
def get_file_parameters(caseroot):
###############################################################################
    """
    Find all XML elements with file="true" attribute. Returns a list
    of file paths.
    """
    raw_xml_file = os.path.join(caseroot, "namelist_scream.xml")
    with open(raw_xml_file, "r") as fd:
        tree    = ET.parse(fd)
        raw_xml = tree.getroot()

    result = []
    for item in raw_xml.findall('.//*[@type="file"]'):
        # Certain configurations may not need a file (e.g., a remap
        # file for SPA may not be needed if the model resolution
        # matches the data file resolution
        if item.text is None or item.text=="":
            continue
        result.append(item.text.strip())

    for item in raw_xml.findall('.//*[@type="array(file)"]'):
        result.extend(refine_type(item.text, force_type="array(file)"))

    # Remove duplicates. Not sure if an error would be warranted if dupes exist
    return list(OrderedDict.fromkeys(result))

###############################################################################
def create_input_data_list_file(case,caseroot):
###############################################################################
    """
    Create the scream.input_data_list file for this case. This will tell CIME
    what to download.
    """
    files_to_download = get_file_parameters(caseroot)

    # Add array parsing knowledge to yaml loader
    loader = yaml.SafeLoader
    loader.add_constructor("!bools",array_constructor)
    loader.add_constructor("!ints",array_constructor)
    loader.add_constructor("!floats",array_constructor)
    loader.add_constructor("!strings",array_constructor)

    # Grab all the output yaml files, open them, and check if horiz_remap_file or vertical_remap_file is used
    rundir   = case.get_value("RUNDIR")
    eamxx_xml_file = os.path.join(caseroot, "namelist_scream.xml")
    with open(eamxx_xml_file, "r") as fd:
        eamxx_xml = ET.parse(fd).getroot()

        scorpio = get_child(eamxx_xml,'Scorpio')
        out_files_xml = get_child(scorpio,"output_yaml_files",must_exist=False)
        #  out_files = out_files_xml.text.split(",") if (out_files_xml is not None and out_files_xml.text is not None) else []
        #  for fn in out_files:
        if (out_files_xml is not None and out_files_xml.text is not None):
            for fn in out_files_xml.text.split(","):
                # Get full name
                src_yaml = os.path.expanduser(os.path.join(fn.strip()))
                dst_yaml = os.path.expanduser(os.path.join(rundir,'data',os.path.basename(src_yaml)))

                # Load file, and look for the remap file entries
                content = yaml.load(open(dst_yaml,"r"),Loader=loader)
                if 'horiz_remap_file' in content.keys():
                    files_to_download += [content['horiz_remap_file']]
                if 'vertical_remap_file' in content.keys():
                    files_to_download += [content['vertical_remap_file']]

    input_data_list_file = "{}/Buildconf/scream.input_data_list".format(caseroot)
    if os.path.exists(input_data_list_file):
        os.remove(input_data_list_file)

    din_loc_root = case.get_value("DIN_LOC_ROOT")
    with open(input_data_list_file, "w") as fd:
        for idx, file_path in enumerate(list(set(files_to_download))):
            # Only add files whose full path starts with the CIME's input data location
            if file_path.startswith(din_loc_root):
                fd.write("scream_dl_input_{} = {}\n".format(idx, file_path))
                if os.path.exists(file_path):
                    if os.path.isdir(file_path):
                        raise IsADirectoryError(f"Input file '{file_path}' is a directory, not a regular file.")
                    if not os.path.isfile(file_path):
                        raise OSError(f"Input file '{file_path}' exists but is not a regular file.")
                    if not os.access(file_path,os.R_OK):
                        try:
                            file_stat = os.stat(file_path)

                            # Get owner and group names
                            owner = pwd.getpwuid(file_stat.st_uid).pw_name
                            group = grp.getgrgid(file_stat.st_gid).gr_name

                            # Get file permissions
                            permissions = stat.filemode(file_stat.st_mode)

                        except Exception as e:
                            raise RuntimeError(f"Error retrieving file info for '{file_path}': {e}") from e

                        curr_user = getpass.getuser()
                        user_info = pwd.getpwnam(curr_user)
                        group_ids = os.getgrouplist(curr_user, user_info.pw_gid)
                        curr_groups = [grp.getgrgid(gid).gr_name for gid in group_ids]

                        raise PermissionError ("Input file exists but it is not readable for current user\n"
                            f" - file name: {file_path}\n"
                            f" - file owner: {owner}\n"
                            f" - file group: {group}\n"
                            f" - permissions: {permissions}\n"
                            f" - current user: {curr_user}\n"
                            f" - current user groups: {curr_groups}\n")

###############################################################################
def do_cime_vars_on_yaml_output_files(case, caseroot):
###############################################################################

    rundir   = case.get_value("RUNDIR")
    eamxx_xml_file = os.path.join(caseroot, "namelist_scream.xml")

    with open(eamxx_xml_file, "r") as fd:
        eamxx_xml = ET.parse(fd).getroot()

    scorpio = get_child(eamxx_xml,'Scorpio')
    out_files_xml = get_child(scorpio,"output_yaml_files",must_exist=False)
    out_files = out_files_xml.text.split(",") if (out_files_xml is not None and out_files_xml.text is not None) else []

    # Add array parsing knowledge to yaml loader
    loader = yaml.SafeLoader
    loader.add_constructor("!bools",array_constructor)
    loader.add_constructor("!ints",array_constructor)
    loader.add_constructor("!floats",array_constructor)
    loader.add_constructor("!strings",array_constructor)

    # We will also change the 'output_yaml_files' entry in scream_input.yaml,
    # to point to the copied files in $rundir/data
    output_yaml_files = []
    scream_input_file = os.path.join(rundir,'data','scream_input.yaml')
    scream_input = yaml.load(open(scream_input_file,"r"),Loader=loader)

    # Determine the physics grid type for use in CIME-var substitution.
    pgt = 'GLL'
    atm_grid = case.get_value('ATM_GRID')
    if '.pg' in atm_grid:
        pgt = 'PG' + atm_grid[-1]

    for fn in out_files:
        # Get full name
        src_yaml = os.path.expanduser(os.path.join(fn.strip()))
        dst_yaml = os.path.expanduser(os.path.join(rundir,'data',os.path.basename(src_yaml)))

        # When expanding CIME vars, we load yaml, expand vars, and dump yaml.
        # In the process, comments are lost. So copy yaml to $rundir first,
        # to leave original file untouched.
        with SharedArea():
            safe_copy(src_yaml,dst_yaml)

        # Now load dst file, and process any CIME var present (if any)
        content = yaml.load(open(dst_yaml,"r"),Loader=loader)
        do_cime_vars(content,case,refine=True,
                     extra={'PHYSICS_GRID_TYPE': pgt})

        # ERS/ERP tests do not work well with INSTANT output. In particular, INSTANT
        # produces an output at t=0, which is not present in the restarted run, and
        # which also causes different timestamp in the file name.
        # Hence, change default output settings to perform a single AVERAGE step at the end of the run
        if case.get_value("TESTCASE") in ["ERP", "ERS"] and content['Averaging Type'].upper()=="INSTANT":
            hist_n = int(case.get_value("HIST_N",resolved=True))
            hist_opt = case.get_value("HIST_OPTION",resolved=True)
            content['output_control']['Frequency'] = hist_n
            content['output_control']['frequency_units'] = hist_opt
            content['output_control']['skip_t0_output'] = True
            print ("ERS/ERP test with INSTANT output detected. Adjusting output control specs:\n")
            print ("  - setting skip_t0_output=true\n")
            print ("  - setting freq and freq_units to HIST_N and HIST_OPTION respectively\n")

        # If frequency_units is not nsteps, verify that we don't request
        # a frequency faster than the model timestep
        if content['output_control']['frequency_units'] in ['nsecs','nmins','nhours']:
            freq  = content['output_control']['Frequency']
            units = content['output_control']['frequency_units']
            dt_out = 1 if units=="nsecs" else 60 if units=="nmins" else 3600
            dt_out = dt_out*int(freq)

            dt_atm = 86400 / case.get_value("ATM_NCPL")
            expect (dt_atm<=dt_out,
                   "Cannot have output frequency faster than atm timestep.\n"
                   f"   yaml file: {fn.strip()}\n"
                   f"   Frequency: {freq}\n"
                   f"   frequency_units: {units}\n"
                   f"   ATM_NCPL: {case.get_value('ATM_NCPL')}\n"
                   f" This yields dt_atm={dt_atm} > dt_output={dt_out}. Please, adjust 'Frequency' and/or 'frequency_units'\n")

        ordered_dump(content, open(dst_yaml, "w"))

        output_yaml_files.append(dst_yaml)

    # Now update the output yaml files entry, and dump the new content
    # of the scream input to YAML file
    scream_input["Scorpio"]["output_yaml_files"] = refine_type(",".join(output_yaml_files),"array(string)")
    with open(scream_input_file, "w") as fd:
        fd.write(
"""################################################################
# Do NOT modify this file. It is generated by scream/buildnml using
# the data from $case/namelist_scream.xml. If you want to make some
# local changes, you can edit this XML file or use atmchange.
################################################################
""")
        ordered_dump(scream_input, fd)
