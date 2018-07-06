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
        case_re = re.compile(r'{}[.]([GC])[.]([^./\s]+)'.format(case))
        value = case_re.sub("{}.ACTION.TESTID".format(case), value)

    if ("/" in value):
        # File path, just return the basename
        return os.path.basename(value)
    elif ("username" in value):
        return ''
    elif (".log." in value):
        # Remove the part that's prone to diff
        components = value.split(".")
        return os.path.basename(".".join(components[0:-1]))
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
def _compare_data(gold_lines, comp_lines, case, offset_method=False):
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
    ('', 0)

    >>> teststr2 = '''
    ... data1
    ... data2 data30
    ... data4 data5 data6
    ... data7 data8 data9 data10
    ... data00
    ... '''
    >>> results,_ = _compare_data(teststr.splitlines(), teststr2.splitlines(), None)
    >>> print(results)
    Inequivalent lines data2 data3 != data2 data30
      NORMALIZED: data2 data3 != data2 data30
    Found extra lines
    data00
    <BLANKLINE>
    >>> teststr3 = '''
    ... data1
    ... data4 data5 data6
    ... data7 data8 data9 data10
    ... data00
    ... '''
    >>> results,_ = _compare_data(teststr3.splitlines(), teststr2.splitlines(), None, offset_method=True)
    >>> print(results)
    Inequivalent lines data4 data5 data6 != data2 data30
      NORMALIZED: data4 data5 data6 != data2 data30
    <BLANKLINE>
    """
    comments = ""
    cnt = 0
    gidx, cidx = 0, 0
    gnum, cnum = len(gold_lines), len(comp_lines)
    while (gidx < gnum or cidx < cnum):
        gidx = _skip_comments_and_whitespace(gold_lines, gidx)
        cidx = _skip_comments_and_whitespace(comp_lines, cidx)

        if (gidx == gnum):
            if (cidx == cnum):
                return comments, cnt
            else:
                comments += "Found extra lines\n"
                comments += "\n".join(comp_lines[cidx:]) + "\n"
                return comments, cnt
        elif (cidx == cnum):
            comments += "Missing lines\n"
            comments += "\n".join(gold_lines[gidx:1]) + "\n"
            return comments, cnt

        gold_value = gold_lines[gidx].strip()
        gold_value = gold_value.replace('"',"'")
        comp_value = comp_lines[cidx].strip()
        comp_value = comp_value.replace('"',"'")

        norm_gold_value = _normalize_string_value(gold_value, case)
        norm_comp_value = _normalize_string_value(comp_value, case)

        if (norm_gold_value != norm_comp_value):
            comments += "Inequivalent lines {} != {}\n".format(gold_value, comp_value)
            comments += "  NORMALIZED: {} != {}\n".format(norm_gold_value, norm_comp_value)
            cnt += 1
        if offset_method and (norm_gold_value != norm_comp_value):
            if gnum > cnum:
                gidx += 1
            else:
                cidx += 1
        else:
            gidx += 1
            cidx += 1

    return comments, cnt

###############################################################################
def compare_files(gold_file, compare_file, case=None):
###############################################################################
    """
    Returns true if files are the same, comments are returned too:
    (success, comments)
    """
    expect(os.path.exists(gold_file), "File not found: {}".format(gold_file))
    expect(os.path.exists(compare_file), "File not found: {}".format(compare_file))

    comments, cnt = _compare_data(open(gold_file, "r").readlines(),
                             open(compare_file, "r").readlines(), case)

    if cnt > 0:
        comments2, cnt2 = _compare_data(open(gold_file, "r").readlines(),
                                        open(compare_file, "r").readlines(),
                                        case, offset_method=True)
        if cnt2 < cnt:
            comments = comments2

    return comments == "", comments
