#!/usr/bin/env python

import re
import unittest
import collections

def flatten(l):
    "http://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists-in-python"
    if l:
        for el in l:
            if isinstance(el, collections.Iterable) and not isinstance(el, (str,bytes)):
# The following is incompatible with python 3.
#            if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
                for sub in flatten(el):
                    yield sub
            else:
                yield el

def parseDirectiveArguments(data):
    
    """Makes a list whose elements are delimited by commas in an input

    string. Commas inside a scope started with brackets, parens, single-
    or double-quotes, are skipped.  Scopes delimited by brackets or
    parentheses are assumed to be well formed.  I.e. we simply count
    opening and closing brackets and parens without regards to any other
    syntax or ordering rules. It would be nice to throw an exception or
    emit warnings when we detect suspicious syntax.
    """
    
    pos = 0; npos = len(data); str=''; maskCommas=False
    scopeCounts = {'[':0,'(':0,'"':0,"'":0} # Assume well formed scopes.
    scopeTerminators = {'[':']','(':')','"':'"',"'":"'"}
    scopeNames = {'(':'parens','[':'brackets','"':'double quotes',"'":'single quotes'}
    while (pos < npos):
        if data[pos] == ',' and not maskCommas:
            return [i for i in flatten([str,parseDirectiveArguments(data[pos+1:])])]
        else:
            for key in scopeCounts.keys():
                if data[pos] == key:
                    scopeCounts[key] = scopeCounts[key] + 1
                    maskCommas = True
                elif data[pos] == scopeTerminators[key]:
                    scopeCounts[key] = scopeCounts[key] - 1
                    if scopeCounts[key] < 0:
                        # Maybe try exceptions...
                        print('parseDirectiveArguments::error: mismatched '+scopeNames[key]+' parenCount < 0 "',str,'" from "',data,'"')
                        return None
                    else:
                        maskCommas = sum(map(abs,scopeCounts.values())) > 0
        str = str + data[pos]
        pos = pos + 1
    return [str]
        
    
class TestParseDirectiveArgs(unittest.TestCase):

    def test_args1(self):
        self.assertEqual(['a','b','c'],parseDirectiveArguments('a,b,c'))

    def test_args2(self):
        self.assertEqual(['a','b(1,2)','c((1,3,z(x,y(4))))'],parseDirectiveArguments('a,b(1,2),c((1,3,z(x,y(4))))'))

    def test_args3(self):
        self.assertEqual(['a','b','c[d,e,f(x,y)]'],parseDirectiveArguments('a,b,c[d,e,f(x,y)]'))

    def test_args4(self):
        self.assertEqual(['a','b','c[d,e,f(x,y]'],parseDirectiveArguments('a,b,c[d,e,f(x,y]'))

    def test_args5(self):
        self.assertEqual(['a','b="This, is, a, test."'],parseDirectiveArguments('a,b="This, is, a, test."'))
        self.assertEqual(["a","b='This, is, a, test.'"],parseDirectiveArguments("a,b='This, is, a, test.'"))

if __name__ == '__main__':
    print('starting')
    unittest.main()
    print('done')
