import os, re, logging

from CIME.utils  import expect
logger=logging.getLogger(__name__)

# pragma pylint: disable=unsubscriptable-object

###############################################################################
def parse_namelists(namelist_lines, filename):
###############################################################################
    """
    Return data in form: {namelist -> {key -> value} }.
      value can be an int, string, list, or dict

    >>> teststr = '''&nml
    ...   val = 'foo'
    ...   aval = 'one','two', 'three'
    ...   maval = 'one', 'two',
    ...       'three', 'four'
    ...   dval = 'one->two', 'three -> four'
    ...   mdval = 'one   -> two',
    ...           'three -> four',
    ...           'five -> six'
    ...   nval = 1850
    ... /
    ...
    ... # Hello
    ...
    ...   &nml2
    ...   val2 = .false.
    ... /
    ... '''
    >>> parse_namelists(teststr.splitlines(), 'foo')
    {'nml': {'dval': {'three': 'four', 'one': 'two'}, 'val': "'foo'", 'maval': ["'one'", "'two'", "'three'", "'four'"], 'aval': ["'one'", "'two'", "'three'"], 'nval': '1850', 'mdval': {'five': 'six', 'three': 'four', 'one': 'two'}}, 'nml2': {'val2': '.false.'}}
    >>> parse_namelists('blah', 'foo')
    Traceback (most recent call last):
        ...
    SystemExit: ERROR: File 'foo' does not appear to be a namelist file, skipping

    >>> teststr = '''&nml
    ... val = 'one', 'two',
    ... val2 = 'three'
    ... /'''
    >>> parse_namelists(teststr.splitlines(), 'foo')
    Traceback (most recent call last):
        ...
    SystemExit: ERROR: In file 'foo', Incomplete multiline variable: 'val'

    >>> teststr = '''&nml
    ... val = 'one', 'two',
    ... /'''
    >>> parse_namelists(teststr.splitlines(), 'foo')
    Traceback (most recent call last):
        ...
    SystemExit: ERROR: In file 'foo', Incomplete multiline variable: 'val'

    >>> teststr = '''&nml
    ... val = 'one', 'two',
    ...       'three -> four'
    ... /'''
    >>> parse_namelists(teststr.splitlines(), 'foo')
    Traceback (most recent call last):
        ...
    SystemExit: ERROR: In file 'foo', multiline list variable 'val' had dict entries
    """

    comment_re = re.compile(r'^[#!]')
    namelist_re = re.compile(r'^&(\S+)$')
    name_re = re.compile(r"^([^\s=']+)\s*=\s*(.+)$")
   #name_re = re.compile(r"^([^\s=':]+)\s*=\s*(.+)$")
    rcline_re = re.compile(r"^([^&\s':]+)\s*:\s*(.+)$")
    dict_re = re.compile(r"^'(\S+)\s*->\s*(\S+)'")
    comma_re = re.compile(r'\s*,\s*')

    rv = {}
    current_namelist = None
    multiline_variable = None # (name, value)
    for line in namelist_lines:

        line = line.strip()

        logger.debug("Parsing line: '%s'" % line)

        if (line == "" or comment_re.match(line) is not None):
            logger.debug("  Line was whitespace or comment, skipping.")
            continue

        rcline = rcline_re.match(line)
        if (rcline is not None):
            # Defining a variable (AKA name)
            name, value = rcline.groups()

            value = value.replace('"',"'")
            logger.debug("  Parsing variable '%s' with data '%s'" % (name, value))

            if 'seq_maps.rc' not in rv:
                rv['seq_maps.rc'] = {}

            expect(name not in rv['seq_maps.rc'], "In file '%s', Duplicate name: '%s'" % (filename, name))
            rv['seq_maps.rc'][name] = value

        elif (current_namelist is None):
            # Must start a namelist
            expect(multiline_variable is None,
                   "In file '%s', Incomplete multiline variable: '%s'" % (filename, multiline_variable[0] if multiline_variable is not None else ""))

            # Unfortunately, other tools were using the old compare_namelists.pl script
            # to compare files that are not namelist files. We need a special error
            # to signify this event
            if (namelist_re.match(line) is None):
                expect(rv != {},
                       "File '%s' does not appear to be a namelist file, skipping" % filename)
                expect(False,
                       "In file '%s', Line '%s' did not begin a namelist as expected" % (filename, line))

            current_namelist = namelist_re.match(line).groups()[0]
            expect(current_namelist not in rv,
                   "In file '%s', Duplicate namelist '%s'" % (filename, current_namelist))

            rv[current_namelist] = {}

            logger.debug("  Starting namelist '%s'" % current_namelist)

        elif (line == "/"):
            # Ends a namelist
            logger.debug("  Ending namelist '%s'" % current_namelist)

            expect(multiline_variable is None,
                   "In file '%s', Incomplete multiline variable: '%s'" % (filename, multiline_variable[0] if multiline_variable is not None else ""))

            current_namelist = None

        elif (name_re.match(line)):
            # Defining a variable (AKA name)
            name, value = name_re.match(line).groups()

            value = value.replace('"',"'")
            logger.debug("  Parsing variable '%s' with data '%s'" % (name, value))

            expect(multiline_variable is None,
                   "In file '%s', Incomplete multiline variable: '%s'" % (filename, multiline_variable[0] if multiline_variable is not None else ""))
            expect(name not in rv[current_namelist], "In file '%s', Duplicate name: '%s'" % (filename, name))

            tokens = [item.strip() for item in comma_re.split(value) if item.strip() != ""]
            if ("->" in value):
                # dict
                rv[current_namelist][name] = {}
                for token in tokens:
                    m = dict_re.match(token)
                    expect(m is not None, "In file '%s', Dict entry '%s' does not match expected format" % (filename, token))
                    k, v = m.groups()
                    rv[current_namelist][name][k] = v
                    logger.debug("    Adding dict entry '%s' -> '%s'" % (k, v))

            elif ("," in value):
                # list
                rv[current_namelist][name] = tokens

                logger.debug("    Adding list entries: %s" % ", ".join(tokens))

            else:
                rv[current_namelist][name] = value

                logger.debug("    Setting to value '%s'" % value)

            if (line.endswith(",")):
                # Value will continue on in subsequent lines
                multiline_variable = (name, rv[current_namelist][name])

                logger.debug("    Var is multiline...")

        elif (multiline_variable is not None):
            # Continuation of list or dict variable
            current_value = multiline_variable[1]
            logger.debug("  Continuing multiline variable '%s' with data '%s'" % (multiline_variable[0], line))
            tokens = [item.strip() for item in comma_re.split(line) if item.strip() != ""]
            if (type(current_value) is list):
                expect("->" not in line, "In file '%s', multiline list variable '%s' had dict entries" % (filename, multiline_variable[0]))
                current_value.extend(tokens)
                logger.debug("    Adding list entries: %s" % ", ".join(tokens))
            elif (type(current_value) is dict):
                for token in tokens:
                    m = dict_re.match(token)
                    expect(m is not None, "In file '%s', Dict entry '%s' does not match expected format" % (filename, token))
                    k, v = m.groups()
                    current_value[k] = v
                    logger.debug("    Adding dict entry '%s' -> '%s'" % (k, v))
            else:
                expect(False, "In file '%s', Continuation should have been for list or dict, instead it was: '%s'" % (filename, type(current_value)))

            if (not line.endswith(",")):
                # Completed
                multiline_variable = None

                logger.debug("    Terminating multiline variable")

        else:
            expect(False, "In file '%s', Unrecognized line: '%s'" % (filename, line))

    return rv

###############################################################################
def normalize_string_value(name, value, case):
###############################################################################
    """
    Some of the string in namelists will contain data that's inherently prone
    to diffs, like file paths, etc. This function attempts to normalize that
    data so that it will not cause diffs.
    """
    # Any occurance of case must be normalized because test-ids might not match
    if (case is not None):
        case_re = re.compile(r'%s[.]([GC]+)[.]([^./\s]+)' % case)
        value = case_re.sub("%s.ACTION.TESTID" % case, value)

    if (name in ["runid", "model_version", "username"]):
        # Don't even attempt to diff these, we don't care
        return name.upper()
    elif (".log." in value):
        # Remove the part that's prone to diff
        components = value.split(".")
        return os.path.basename(".".join(components[0:-1]))
    elif (":" in value):
        items = value.split(":")
        items = [normalize_string_value(name, item, case) for item in items]
        return ":".join(items)
    elif ("/" in value):
        # File path, just return the basename unless its a seq_maps.rc mapping
        # mapname or maptype
        if "mapname" not in name and "maptype" not in name:
            return os.path.basename(value)
        else:
            return value
    else:
        return value

###############################################################################
def compare_values(name, gold_value, comp_value, case, do_print=False):
###############################################################################
    """
    Compare values for a specific variable in a namelist.
    """
    if (type(gold_value) != type(comp_value)):
        if do_print:
            print "  variable '%s' did not have expected type '%s', instead is type '%s'" % \
                (name, type(gold_value), type(comp_value))
        return False

    rv = True
    if (type(gold_value) is list):
        # Note, list values remain order sensitive
        for idx, gold_value_list_item in enumerate(gold_value):
            if (idx < len(comp_value)):
                rv &= compare_values("%s list item %d" % (name, idx),
                                     gold_value_list_item, comp_value[idx], case,
                                     do_print)
            else:
                rv = False
                if do_print:
                    print "  list variable '%s' missing value %s" % (name, gold_value_list_item)

        if (len(comp_value) > len(gold_value)):
            for comp_value_list_item in comp_value[len(gold_value):]:
                rv = False
                if do_print:
                    print "  list variable '%s' has extra value %s" % (name, comp_value_list_item)

    elif (type(gold_value) is dict):
        for key, gold_value_dict_item in gold_value.iteritems():
            if (key in comp_value):
                rv &= compare_values("%s dict item %s" % (name, key),
                                     gold_value_dict_item, comp_value[key], case,
                                     do_print)
            else:
                rv = False
                if do_print:
                    print "  dict variable '%s' missing key %s with value %s" \
                        % (name, key, gold_value_dict_item)

        for key in comp_value:
            if (key not in gold_value):
                rv = False
                if do_print:
                    print "  dict variable '%s' has extra key %s with value %s" \
                        % (name, key, comp_value[key])

    else:
        expect(type(gold_value) is str, "Unexpected type found: '%s'" % type(gold_value))
        norm_gold_value = normalize_string_value(name, gold_value, case)
        norm_comp_value = normalize_string_value(name, comp_value, case)

        if (norm_gold_value != norm_comp_value):
            rv = False
            if do_print:
                print "  BASE: %s = %s" % (name, norm_gold_value)
                print "  COMP: %s = %s" % (name, norm_comp_value)

    return rv

###############################################################################
def compare_namelists(gold_namelists, comp_namelists, case):
###############################################################################
    """
    Compare two namelists. Print diff information if any. Return true if
    equivalent.

    Expect args in form: {namelist -> {key -> value} }.
      value can be an int, string, list, or dict

    >>> teststr = '''&nml
    ...   val = 'foo'
    ...   aval = 'one','two', 'three'
    ...   maval = 'one', 'two', 'three', 'four'
    ...   dval = 'one -> two', 'three -> four'
    ...   mdval = 'one -> two', 'three -> four', 'five -> six'
    ...   nval = 1850
    ... /
    ... &nml2
    ...   val2 = .false.
    ... /
    ... '''
    >>> compare_namelists(parse_namelists(teststr.splitlines(), 'foo'), parse_namelists(teststr.splitlines(), 'bar'), None)
    True

    >>> teststr1 = '''&nml1
    ...   val11 = 'foo'
    ... /
    ... &nml2
    ...   val21 = 'foo'
    ...   val22 = 'foo', 'bar', 'baz'
    ...   val23 = 'baz'
    ...   val24 = '1 -> 2', '2 -> 3', '3 -> 4'
    ... /
    ... &nml3
    ...   val3 = .false.
    ... /'''
    >>> teststr2 = '''&nml01
    ...   val11 = 'foo'
    ... /
    ... &nml2
    ...   val21 = 'foo0'
    ...   val22 = 'foo', 'bar0', 'baz'
    ...   val230 = 'baz'
    ...   val24 = '1 -> 20', '2 -> 3', '30 -> 4'
    ... /
    ... &nml3
    ...   val3 = .false.
    ... /'''
    >>> compare_namelists(parse_namelists(teststr1.splitlines(), 'foo'), parse_namelists(teststr2.splitlines(), 'bar'), None)
    Differences in namelist 'nml2':
      BASE: val21 = 'foo'
      COMP: val21 = 'foo0'
      BASE: val22 list item 1 = 'bar'
      COMP: val22 list item 1 = 'bar0'
      missing variable: 'val23'
      BASE: val24 dict item 1 = 2
      COMP: val24 dict item 1 = 20
      dict variable 'val24' missing key 3 with value 4
      dict variable 'val24' has extra key 30 with value 4
      found extra variable: 'val230'
    Missing namelist: nml1
    Found extra namelist: nml01
    False

    >>> teststr1 = '''&rad_cnst_nl
    ... icecldoptics           = 'mitchell'
    ... logfile                = 'cpl.log.150514-001533'
    ... case_name              = 'ERB.f19_g16.B1850C5.skybridge_intel.C.150513-230221'
    ... runid                  = 'FOO'
    ... model_version          = 'cam5_3_36'
    ... username               = 'jgfouca'
    ... iceopticsfile          = '/projects/ccsm/inputdata/atm/cam/physprops/iceoptics_c080917.nc'
    ... liqcldoptics           = 'gammadist'
    ... liqopticsfile          = '/projects/ccsm/inputdata/atm/cam/physprops/F_nwvl200_mu20_lam50_res64_t298_c080428.nc'
    ... mode_defs              = 'mam3_mode1:accum:=', 'A:num_a1:N:num_c1:num_mr:+',
    ...   'A:so4_a1:N:so4_c1:sulfate:/projects/ccsm/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+', 'A:pom_a1:N:pom_c1:p-organic:/projects/ccsm/inputdata/atm/cam/physprops/ocpho_rrtmg_c101112.nc:+',
    ...   'A:soa_a1:N:soa_c1:s-organic:/projects/ccsm/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+', 'A:bc_a1:N:bc_c1:black-c:/projects/ccsm/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
    ...   'A:dst_a1:N:dst_c1:dust:/projects/ccsm/inputdata/atm/cam/physprops/dust4_rrtmg_c090521.nc:+', 'A:ncl_a1:N:ncl_c1:seasalt:/projects/ccsm/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc',
    ...   'mam3_mode2:aitken:=', 'A:num_a2:N:num_c2:num_mr:+',
    ...   'A:so4_a2:N:so4_c2:sulfate:/projects/ccsm/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+', 'A:soa_a2:N:soa_c2:s-organic:/projects/ccsm/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+',
    ...   'A:ncl_a2:N:ncl_c2:seasalt:/projects/ccsm/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc', 'mam3_mode3:coarse:=',
    ...   'A:num_a3:N:num_c3:num_mr:+', 'A:dst_a3:N:dst_c3:dust:/projects/ccsm/inputdata/atm/cam/physprops/dust4_rrtmg_c090521.nc:+',
    ...   'A:ncl_a3:N:ncl_c3:seasalt:/projects/ccsm/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+', 'A:so4_a3:N:so4_c3:sulfate:/projects/ccsm/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc'
    ... rad_climate            = 'A:Q:H2O', 'N:O2:O2', 'N:CO2:CO2',
    ...   'N:ozone:O3', 'N:N2O:N2O', 'N:CH4:CH4',
    ...   'N:CFC11:CFC11', 'N:CFC12:CFC12', 'M:mam3_mode1:/projects/ccsm/inputdata/atm/cam/physprops/mam3_mode1_rrtmg_c110318.nc',
    ...   'M:mam3_mode2:/projects/ccsm/inputdata/atm/cam/physprops/mam3_mode2_rrtmg_c110318.nc', 'M:mam3_mode3:/projects/ccsm/inputdata/atm/cam/physprops/mam3_mode3_rrtmg_c110318.nc'
    ... /'''
    >>> teststr2 = '''&rad_cnst_nl
    ... icecldoptics           = 'mitchell'
    ... logfile                = 'cpl.log.150514-2398745'
    ... case_name              = 'ERB.f19_g16.B1850C5.skybridge_intel.C.150513-1274213'
    ... runid                  = 'BAR'
    ... model_version          = 'cam5_3_36'
    ... username               = 'hudson'
    ... iceopticsfile          = '/something/else/inputdata/atm/cam/physprops/iceoptics_c080917.nc'
    ... liqcldoptics           = 'gammadist'
    ... liqopticsfile          = '/something/else/inputdata/atm/cam/physprops/F_nwvl200_mu20_lam50_res64_t298_c080428.nc'
    ... mode_defs              = 'mam3_mode1:accum:=', 'A:num_a1:N:num_c1:num_mr:+',
    ...   'A:so4_a1:N:so4_c1:sulfate:/something/else/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+', 'A:pom_a1:N:pom_c1:p-organic:/something/else/inputdata/atm/cam/physprops/ocpho_rrtmg_c101112.nc:+',
    ...   'A:soa_a1:N:soa_c1:s-organic:/something/else/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+', 'A:bc_a1:N:bc_c1:black-c:/something/else/inputdata/atm/cam/physprops/bcpho_rrtmg_c100508.nc:+',
    ...   'A:dst_a1:N:dst_c1:dust:/something/else/inputdata/atm/cam/physprops/dust4_rrtmg_c090521.nc:+', 'A:ncl_a1:N:ncl_c1:seasalt:/something/else/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc',
    ...   'mam3_mode2:aitken:=', 'A:num_a2:N:num_c2:num_mr:+',
    ...   'A:so4_a2:N:so4_c2:sulfate:/something/else/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc:+', 'A:soa_a2:N:soa_c2:s-organic:/something/else/inputdata/atm/cam/physprops/ocphi_rrtmg_c100508.nc:+',
    ...   'A:ncl_a2:N:ncl_c2:seasalt:/something/else/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc', 'mam3_mode3:coarse:=',
    ...   'A:num_a3:N:num_c3:num_mr:+', 'A:dst_a3:N:dst_c3:dust:/something/else/inputdata/atm/cam/physprops/dust4_rrtmg_c090521.nc:+',
    ...   'A:ncl_a3:N:ncl_c3:seasalt:/something/else/inputdata/atm/cam/physprops/ssam_rrtmg_c100508.nc:+', 'A:so4_a3:N:so4_c3:sulfate:/something/else/inputdata/atm/cam/physprops/sulfate_rrtmg_c080918.nc'
    ... rad_climate            = 'A:Q:H2O', 'N:O2:O2', 'N:CO2:CO2',
    ...   'N:ozone:O3', 'N:N2O:N2O', 'N:CH4:CH4',
    ...   'N:CFC11:CFC11', 'N:CFC12:CFC12', 'M:mam3_mode1:/something/else/inputdata/atm/cam/physprops/mam3_mode1_rrtmg_c110318.nc',
    ...   'M:mam3_mode2:/something/else/inputdata/atm/cam/physprops/mam3_mode2_rrtmg_c110318.nc', 'M:mam3_mode3:/something/else/inputdata/atm/cam/physprops/mam3_mode3_rrtmg_c110318.nc'
    ... /'''
    >>> compare_namelists(parse_namelists(teststr1.splitlines(), 'foo'), parse_namelists(teststr2.splitlines(), 'bar'), 'ERB.f19_g16.B1850C5.skybridge_intel')
    True
    """
    rv = True

    # We want to keep lists of differences and print results in a second pass,
    # in order to ensure that the order is not scrambled when we change Python
    # versions and/or parse_namelists implementation details.
    # (Would using an ordered dict in parse_namelists be easier?)
    different_namelists = []
    missing_namelists = []
    for namelist, gold_names in gold_namelists.iteritems():
        if (namelist not in comp_namelists):
            rv = False
            missing_namelists.append(namelist)
        else:
            comp_names = comp_namelists[namelist]
            namelists_equal = True
            for name, gold_value in gold_names.iteritems():
                if (name not in comp_names):
                    namelists_equal = False
                    break
                else:
                    comp_value = comp_names[name]
                    if not compare_values(name, gold_value, comp_value, case):
                        namelists_equal = False
                        break

            for name in comp_names:
                if (name not in gold_names):
                    namelists_equal = False
                    break

            if not namelists_equal:
                different_namelists.append(namelist)
                rv = False

    for namelist in sorted(different_namelists):
        print "Differences in namelist '%s':" % namelist
        gold_names = gold_namelists[namelist]
        comp_names = comp_namelists[namelist]
        for name in sorted(gold_names.keys()):
            if (name not in comp_names):
                print "  missing variable: '%s'" % name
            else:
                gold_value = gold_names[name]
                comp_value = comp_names[name]
                compare_values(name, gold_value, comp_value, case, do_print=True)
        for name in sorted(comp_names.keys()):
            if name not in gold_names:
                print "  found extra variable: '%s'" % name

    for namelist in sorted(missing_namelists):
        print "Missing namelist:", namelist

    extra_namelists = []
    for namelist in comp_namelists:
        if (namelist not in gold_namelists):
            rv = False
            extra_namelists.append(namelist)

    for namelist in sorted(extra_namelists):
        print "Found extra namelist:", namelist

    return rv

###############################################################################
def compare_namelist_files(gold_file, compare_file, case=None):
###############################################################################
    expect(os.path.exists(gold_file), "File not found: %s" % gold_file)
    expect(os.path.exists(compare_file), "File not found: %s" % compare_file)

    gold_namelists = parse_namelists(open(gold_file, "r").readlines(), gold_file)
    comp_namelists = parse_namelists(open(compare_file, "r").readlines(), compare_file)

    return compare_namelists(gold_namelists, comp_namelists, case)

###############################################################################
def is_namelist_file(file_path):
###############################################################################
    try:
        compare_namelist_files(file_path, file_path)
    except SystemExit as e:
        assert "does not appear to be a namelist file" in str(e), str(e)
        return False
    return True

