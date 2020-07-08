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

###############################################################################
def compare_runconfigfiles(gold_file, compare_file, case=None):
###############################################################################
    """
    Returns true if files are the same, comments are returned too:
    (success, comments)
    """
    expect(os.path.exists(gold_file), "File not found: {}".format(gold_file))
    expect(os.path.exists(compare_file), "File not found: {}".format(compare_file))

    #create dictionary's of the runconfig files and compare them
    gold_dict = _parse_runconfig(gold_file)
    compare_dict = _parse_runconfig(compare_file)

    comments = findDiff(gold_dict, compare_dict, case=case)
    comments = comments.replace(" d1", " " + gold_file)
    comments = comments.replace(" d2", " " + compare_file)
    # this picks up the case that an entry in compare is not in gold
    if comments == "":
        comments = findDiff(compare_dict, gold_dict, case=case)
        comments = comments.replace(" d2", " " + gold_file)
        comments = comments.replace(" d1", " " + compare_file)

    return comments == "", comments

def _parse_runconfig(filename):
    runconfig = {}
    inrunseq = False
    insubsection = None
    subsection_re = re.compile(r'\s*(\S+)::')
    group_re = re.compile(r'\s*(\S+)\s*:\s*(\S+)')
    var_re = re.compile(r'\s*(\S+)\s*=\s*(\S+)')
    with open(filename, "r") as fd:
        for line in fd:
            # remove comments
            line = line.split('#')[0]
            subsection_match = subsection_re.match(line)
            group_match = group_re.match(line)
            var_match = var_re.match(line)
            if re.match(r'\s*runSeq\s*::', line):
                runconfig['runSeq'] = []
                inrunseq = True
            elif re.match(r'\s*::\s*', line):
                inrunseq = False
            elif inrunseq:
                runconfig['runSeq'].append(line)
            elif subsection_match:
                insubsection = subsection_match.group(1)
                runconfig[insubsection] = {}
            elif group_match:
                runconfig[group_match.group(1)] = group_match.group(2)
            elif insubsection and var_match:
                runconfig[insubsection][var_match.group(1)] = var_match.group(2)
    return runconfig

def findDiff(d1, d2, path="", case=None):
    comment = ""
    for k in d1.keys():
        if not k in d2:
            comment += path + ":\n"
            comment +=  k + " as key not in d2\n"
        else:
            if type(d1[k]) is dict:
                if path == "":
                    path = k
                else:
                    path = path + "->" + k
                comment += findDiff(d1[k],d2[k], path=path, case=case)
            else:
                if case in d1[k]:
                    pass
                elif "username" in k:
                    pass
                elif "logfile" in k:
                    pass
                elif d1[k] != d2[k]:
                    comment += path+":\n"
                    comment += " - {} : {}\n".format(k,d1[k])
                    comment += " + {} : {}\n".format(k,d2[k])
    return comment
