#!/usr/bin/env python

import unittest

import re
import string
import collections

def countBrackets(str,flagNegative=False):
    count = 0
    ok = True
    for iChar in str:
        if iChar == '[':
            count = count + 1
        elif iChar == ']':
            count = count -1
        if flagNegative:
            if ok:
                if count < 0:
                    ok = False
                    print('countBrackets::ERROR::mismatched brackets first detection',count)
                    print('str: ',str)
    if flagNegative:
        if count < 0:
            print('countBrackets::ERROR::mismatched brackets',count)
    return count

def rejoinBracketed_(textArray):
    if not textArray:
        return None
    if isinstance(textArray,str):
        return textArray
    if not isinstance(textArray,list):
        print("rejoinBracketed::ERROR::input is neither string or list thereof.")
        return None
    nBrackets = countBrackets(textArray[0],flagNegative=True)
    
    if len(textArray) == 1:
        return textArray
    elif nBrackets == 0:
        return [textArray[0],rejoinBracketed(textArray[1:])]
    elif nBrackets > 0:
        nextItem = textArray[0]
        lenTA = len(textArray); i = 1
        nBrackets0 = nBrackets # Save the original number of unclosed pairs.
        while i < lenTA:
            nextItem = nextItem+textArray[i]
            nBrackets = nBrackets + countBrackets(textArray[i],flagNegative=False)
            i = i + 1
            if nBrackets == 0:
                break # normal end
        if nBrackets != 0:
            print('rejoinBracketed::ERROR::bracket mismatch::nBrackets 0 !=',nBrackets)
        if i < lenTA:
            return [nextItem,rejoinBracketed(textArray[i:])]
        else:
            return [nextItem]
    return None

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

def rejoinBracketed(textArray):
    a = rejoinBracketed_(textArray)
    if isinstance(a,str):
        return a
    return [i for i in flatten(a)]

def splitKeepDelimiters(str, delimiter):
    return [e for e in re.split('([^,'+delimiter+']*)',str) if e != '']

def removeDelimiters(listStr, delimiter):
    return [e for e in listStr if e != delimiter]

def parseBrackets(str):
    """Splits a string on comma-delimited tokens while respecting text in nested brackets.
    E.g. parseBrackets("a,[b,[c,d]],f") => ['a','[b,[c,d]]','f']."""
    return removeDelimiters(rejoinBracketed(splitKeepDelimiters(str,',')),',')

class TestRejoinBracketed(unittest.TestCase):

    def testRejoinBracketed(self):
        self.assertEqual('a,b',rejoinBracketed_('a,b'))
        self.assertEqual('a,b',rejoinBracketed('a,b'))
        self.assertEqual(['a,b'],rejoinBracketed(['a,b']))
        self.assertEqual(['[a]','b'],rejoinBracketed(['[a]','b']))
        self.assertEqual(['[a,b]'],rejoinBracketed(['[a',',','b]']))
        self.assertEqual(['[a,b],c,d'],rejoinBracketed(['[a',',','b],c,d']))
        self.assertEqual(['a,[b,c],d'],rejoinBracketed(['a,[b,','c],d']))
        self.assertEqual(['a',',','[b,c]',',','d'],rejoinBracketed(splitKeepDelimiters('a,[b,c],d',',')))
        self.assertEqual(['a','[b,c]','d'],\
                             removeDelimiters(rejoinBracketed(splitKeepDelimiters('a,[b,c],d',',')),','))
        self.assertEqual(['[a,b]','[c,d]','e'],\
                             removeDelimiters(rejoinBracketed(splitKeepDelimiters('[a,b],[c,d],e',',')),','))
        self.assertEqual(['[[a,b],c]','d','e','[f,g]'],\
                             removeDelimiters(rejoinBracketed(splitKeepDelimiters('[[a,b],c],d,e,[f,g]',',')),','))

    def testParseBrackets(self):
        self.assertEqual(['[[a,b],c]','d','e','[f,g]'],parseBrackets('[[a,b],c],d,e,[f,g]'))
        self.assertEqual(['[[a,b],c,[x,y]]','d','e','[f,g]'],parseBrackets('[[a,b],c,[x,y]],d,e,[f,g]'))
        self.assertEqual(['a','[b,[c,d]]','f'],parseBrackets("a,[b,[c,d]],f"))

if __name__ == "__main__":
    unittest.main()   

