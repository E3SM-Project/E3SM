import os, re

from CIME.utils import expect

###############################################################################
def _normalize_string_value(value, case):
###############################################################################
    """
    Some of the strings are inherently prone to diffs, like file
    paths, etc. This function attempts to normalize that data so that
    it will not cause diffs.
    """
    # Any occurance of case must be normalized because test-ids might not match
    if (case is not None):
        case_re = re.compile(r'%s[.]([GC])[.]([^./\s]+)' % case)
        value = case_re.sub("%s.ACTION.TESTID" % case, value)

    if ("/" in value):
        # File path, just return the basename
        return os.path.basename(value)
    else:
        return value

###############################################################################
def _skip_comments_and_whitespace(lines, idx):
###############################################################################
    """
    Starting at idx, return next valid idx of lines that contains real data
    """
    if (idx == len(lines)):
        return idx

    comment_re = re.compile(r'^[#!]')

    lines_slice = lines[idx:]
    for line in lines_slice:
        line = line.strip()
        if (comment_re.match(line) is not None or line == ""):
            idx += 1
        else:
            return idx

    return idx

###############################################################################
def _compare_data(gold_lines, comp_lines, case):
###############################################################################
    """
    >>> teststr = '''
    ... data1
    ... data2 data3
    ... data4 data5 data6
    ...
    ... # Comment
    ... data7 data8 data9 data10
    ... '''
    >>> _compare_data(teststr.splitlines(), teststr.splitlines(), None)
    ''

    >>> teststr2 = '''
    ... data1
    ... data2 data30
    ... data4 data5 data6
    ... data7 data8 data9 data10
    ... data00
    ... '''
    >>> results = _compare_data(teststr.splitlines(), teststr2.splitlines(), None)
    >>> print results
    Inequivalent lines data2 data3 != data2 data30
      NORMALIZED: data2 data3 != data2 data30
    Found extra lines
    data00
    <BLANKLINE>
    """
    comments = ""
    gidx, cidx = 0, 0
    gnum, cnum = len(gold_lines), len(comp_lines)
    while (gidx < gnum or cidx < cnum):
        gidx = _skip_comments_and_whitespace(gold_lines, gidx)
        cidx = _skip_comments_and_whitespace(comp_lines, cidx)

        if (gidx == gnum):
            if (cidx == cnum):
                return comments
            else:
                comments += "Found extra lines\n"
                comments += "\n".join(comp_lines[cidx:]) + "\n"
                return comments
        elif (cidx == cnum):
            comments += "Missing lines\n"
            comments += "\n".join(gold_lines[gidx:1]) + "\n"
            return comments

        gold_value = gold_lines[gidx].strip()
        gold_value = gold_value.replace('"',"'")
        comp_value = comp_lines[cidx].strip()
        comp_value = comp_value.replace('"',"'")

        norm_gold_value = _normalize_string_value(gold_value, case)
        norm_comp_value = _normalize_string_value(comp_value, case)

        if (norm_gold_value != norm_comp_value):
            comments += "Inequivalent lines %s != %s\n" % (gold_value, comp_value)
            comments += "  NORMALIZED: %s != %s\n" % (norm_gold_value, norm_comp_value)

        gidx += 1
        cidx += 1

    return comments

###############################################################################
def compare_files(gold_file, compare_file, case=None):
###############################################################################
    """
    Returns true if files are the same, comments are returned too:
    (success, comments)
    """
    expect(os.path.exists(gold_file), "File not found: %s" % gold_file)
    expect(os.path.exists(compare_file), "File not found: %s" % compare_file)

    comments = _compare_data(open(gold_file, "r").readlines(),
                             open(compare_file, "r").readlines(),
                             case)
    return comments == "", comments
