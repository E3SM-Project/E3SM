#!/usr/bin/env python

import unittest

debugLevel = 0
debug = 100

whitespace = ' \t\n'

def parseArgs (s) :
    """ Parse a string the form of a=1,b=[1,"c"], etc. Uses Python's 'eval' to bring.
    """
    args = {}
    pos = 0
    sstripped = s.strip()
    lhs0 = ''
    rhs0 = ''
    while pos < len(s) :
        c = s[pos]
        if c in whitespace :
            # ignore
            pass
        else:
            if c == '=' :
                if lhs0 == '' :
                    print 'syntax error: lhs empty'
                # use rhs0 as key in dictionary
                equalDone = False
                pos += 1
                firstCharRHS = ''
                firstCharQuote = False
                tripleQuote = False
                while pos < len(s) and not equalDone :
                    if debugLevel > debug : print 'rhs0-0: ',rhs0
                    c = s[pos]
                    if c in whitespace :
                        # ignore
                        pass
                    elif c != '[':
                        # not an array
                        if firstCharRHS == '' :
                            firstCharRHS=c
                            firstCharQuote = firstCharRHS in ['"',"'"]
                            tripleQuote = s[pos:pos+3] == firstCharRHS*3
                            if debugLevel > debug : print 'quote? - s: <',s[pos:pos+5],'> ...',' fc*3=<',firstCharRHS*3,'>?=<',s[pos:pos+3],'>'
                            if debugLevel > debug : print 'quote?      ',firstCharRHS,firstCharQuote,tripleQuote
                            if tripleQuote :
                                rhs0 = rhs0 + s[pos:pos+3]; pos = pos +3
                            else:
                                rhs0 = rhs0 + s[pos] ; pos = pos + 1
                        while pos < len(s) and not equalDone :
                            if debugLevel > debug : print 'rhs0-aa: ',rhs0
                            c = s[pos]
                            if not firstCharQuote :
                                if c == ',' :
                                    equalDone = True
                                elif c == '[' or c == ']' :
                                    print 'malformed statement'
                                else:
                                    rhs0 = rhs0 + c
                                pos += 1
                            else:
                                if debugLevel > debug : print 'c: ',c
                                if c == firstCharRHS :
                                    if tripleQuote :
                                        equalDone = s[pos:pos+3] == firstCharRHS*3
                                        if debugLevel > debug : print 'triple!'
                                        if equalDone :
                                            rhs0 = rhs0 + s[pos:pos+3]; pos = pos + 3
                                        else:
                                            rhs0 = rhs0 + c
                                    else:
                                        equalDone = True
                                        rhs0 = rhs0 + c; pos = pos + 1
                                        if debugLevel > debug : print 'firstChar done! -- ',c,' = ',firstCharRHS
                                    if equalDone :
                                        # skip white space
                                        p0 = pos + 1
                                        p0Done = not (p0 < len(s))
                                        while not p0Done :
                                            c = s[p0]; 
                                            p0Done = (p0 < len(s)) and not (c in whitespace)
                                            if not p0Done : p0 = p0 + 1
                                        pos = p0
                                else:
                                    rhs0 = rhs0 + c; pos = pos + 1
                            if debugLevel > debug : print 'rhs0-1: ',rhs0,equalDone,pos
                    elif c == '[' :
                        # an array
                        firstCharRHS=c
                        rhs0 = c; pos += 1
                        while pos < len(s) and not equalDone :
                            c = s[pos]
                            rhs0 = rhs0 + c
                            if c == ']':
                                equalDone = True
                                while pos < len(s) and s[pos] != ',' :
                                    pos += 1
                            pos += 1
                if debugLevel > debug : print 'l,r: ',lhs0,'---',rhs0,'---',len(rhs0)
                rhs0 = eval(rhs0)
                #                if not firstCharQuote : 
                #                    rhs0 = eval(rhs0)
                if debugLevel > debug : print 'rhs0-xx: ',rhs0
                args[lhs0]=rhs0
                rhs0=''; lhs0=''
            else:
                lhs0 = lhs0 + c
                pos += 1
    if debugLevel > debug : print 'args: ',args
    return args

def testParseArgs(s):
    print 'tpa: '+s,parseArgs(s)

# testParseArgs('arg1=1')
# testParseArgs('arg1=1,arg2=2')
# testParseArgs('arg1=[1,2]')

class TestParseArgs(unittest.TestCase):

    def test_ParseArgs_OneArgWithBrackets1(self):
        self.assertEqual(parseArgs('arg1="[1,2]"'), \
                         {'arg1':'[1,2]'})

    def test_ParseArgs_OneArgWithBrackets2(self):
        self.assertEqual(parseArgs("arg1='[1,2]'"), \
                         {'arg1':'[1,2]'})

    def test_ParseArgs_OneArgWithBrackets3(self):
        self.assertEqual(parseArgs("sq3='''[1,2]'''"), \
                         {'sq3':'[1,2]'})

    def test_ParseArgs_OneArgWithBrackets4(self):
        self.assertEqual(parseArgs('dq3="""[1,2]"""'), \
                         {'dq3':'[1,2]'})

    def test_ParseArgs_OneArgWithBrackets5(self):
        self.assertEqual(parseArgs('dq3="""[1,2]""",:test="timing"'), \
        {'dq3':'[1,2]',':test':'timing'})
        
    def test_ParseArgs_OneArgWithBrackets6(self):
        self.assertEqual(parseArgs('dq3="[1,2]",:test="timing"'), \
        {'dq3':'[1,2]',':test':'timing'})

    def test_ParseArgs_OneArgWithBrackets7(self):
        self.assertEqual(parseArgs('dq3="[1,2]",  :test="timing"'), \
        {'dq3':'[1,2]',':test':'timing'})

    def test_ParseArgs_oneArg(self):
        self.assertEqual(parseArgs('arg1=1'),{'arg1':1})
        
    def test_ParseArgs_twoArgs1(self):
        self.assertEqual(parseArgs('arg1=1,arg2=2'),{'arg1':1,'arg2':2})

    def test_ParseArgs_twoArgs2(self):
        self.assertEqual(parseArgs('arg1=1,:arg2=2'),{'arg1':1,':arg2':2})

    def test_ParseArgs_oneArgArray1(self):
        self.assertEqual(parseArgs('arg1=[1,2]'),{'arg1':[1,2]})

    def test_ParseArgs_TwoArgArray(self):
        self.assertEqual(parseArgs('arg1=[1,2],arg2=[3,4]'),{'arg1':[1,2],'arg2':[3,4]})

    def test_ParseArgs_ThreeArgs(self):
        self.assertEqual(parseArgs('arg1=[1,2],arg2=3,arg3=["a","b12"]'), \
                         {'arg1':[1,2],'arg2':3,'arg3':['a','b12']})


if __name__ == '__main__' :
    unittest.main()
