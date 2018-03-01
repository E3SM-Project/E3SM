#!/usr/bin/env python

import re
import unittest
from textwrap import dedent

debug = 10

class dataString() :
    "A stand-in until we get a stream going."
    def __init__(self,data):
        self.data = data
        self.position = 0
    def insert(self,pos,insertData):
        self.data = self.data[:pos] + insertData + self.data[pos:]
    def getLength(self):
        return len(self.data)
    def getPosition(self):
        return self.position
    def setPosition(self,pos):
        self.position = pos
    def getItem(self,pos):
        return self.data[pos]
    def getDataAtPosition(self,pos):
        return self.data[pos]
    def getData(self):
        return self.data
    def getSlice(self,i0,i1):
        return self.data[i0:i1]
    def getSliceForward(self,i0):
        return self.data[i0:]    
    def removeSlice(self,i0,i1) :
        self.data = self.data[:i0]+self.data[i1:]
        return self.data
    def getCurrentData(self):
        return self.data[self.position]
    def insertAtCurrent(self,includeData):
        self.insert(self.getPosition(),includeData)
    def append(self,appendData):
        # Maybe also need to check if appendData can be cat'd with str.
        if appendData :
            self.data = self.data + appendData
    def advanceAndGetNextData(self) :
        # better if we could throw exceptions here...
        self.position += 1
        return self.data[self.position]
    def validPosition(self,position):
        return ( 0 <= position ) and ( position <= self.getLength() )
    def findToEnd(self,start,s):
        return self.data[start:].find(s)
    def match(self,sre,start,end):
        return sre.match(self.data[start:end])
    def matchToEnd(self,sre,start):
        return sre.match(self.data[start:])
    def searchToEnd(self,sre,start):
        return sre.search(self.data[start:])
    def searchToPosition(self,sre,start,position):
        return sre.search(self.data[start:position])
    def finditerToEnd(self,sre,start):
        return sre.finditer(self.data[start:])
    def finditerToPosition(self,sre,start,position):
        return sre.finditer(self.data[start:position])
    

procDirectives = {
    }

class procDirective() :
    def __init__(self,name):
        self.name = name
        self.newPosition = -1
        self.tokens={}
        self.TokenREs={}
        procDirectives[self.name] = self
    def getLength(self) :
        return len(self.name)
    def match(self,name) :
        return self.name == name
    def evaluate(self,data,pos) :
        return 'data processed by '+ self.name
    def getNewPosition(self) :
        return self.newPosition

# Need to fix...
    def addTokenRE(self,args,key,defaultToken,prefix=r'''(?i)[ \t]*''',postfix='') :
        """
        Add a token/create an RE with a prefix that by default ignores preceding whitespace.
        Stores the RE in a dictionary for this directive.  Note this currently expects
        complex tokens like <EndToken> not something as overloaded as a close paren.
        """
        retToken = defaultToken
        # Delete the current token definition.
        if key in args.keys():
           retToken = args[key][:]
           del args[key]
        # Construct the RE and add the definition.
        tokenREString = prefix + retToken + postfix
        tokenRE = re.compile(tokenREString)
        # self.tokens[key] = 
        self.TokenREs[key] = tokenRE
        self.tokens[key]   = retToken[:]
        return tokenRE
    def searchTokenToEnd(self,key,data,start):
        m = data.searchToEnd(self.TokenREs[key],start)
        return m
    def searchTokenToPosition(self,key,data,start,end):
        m = data.searchToPosition(self.TokenREs[key],start,end)
        return m
    def finditerTokenToPosition(self,key,data,start,end):
        m = data.finditerToPosition(self.TokenREs[key],start,end)
        return m
    def makeTokenErrorMessage(self,msg,key) :
        retMessage = msg+' for '+str(key)+' -> '+self.tokens[key]
        print retMessage
        return retMessage

StringRE = re.compile(r'''(\".*\")|(\'.*\')''')
def matchString(s) :
    return StringRE.match(s)

StringRE_ErrorEOL = re.compile(r'''(\".*\n)|(\'.*\n)''')
def matchStringErrorEOL(s) :
    return StringRE_ErrorEOL.match(s)

StringParensRE = re.compile(r'''\(.*\)''')
def matchStringParensRE(s) :
    return StringParensRE.match(s)

StringParensAndCRRE = re.compile(r'''\((.*\n*)*\)''')
def matchStringParensAndCRRE(s) :
    return StringParensAndCRRE.match(s)

def removeQuotes(s) :
    r = s.replace('\'','')
    return r.replace('\"','')

def removeParens(s) :
    r = s.replace('(','')
    return r.replace(')','')

def getPCArgs(data):
    "PC = Parentheses and Commas"
    pos = 0
    if data[pos] != '(' :
        args = []
        len_arg = 0
    else:
        findex = data[pos:].find(')')
        if findex > -1 :
            m = matchStringParensAndCRRE(data[pos:pos+findex+1])
        else:
            print 'argumentError: closing parenthesis not found'
            m = None
        arg = m.group()
        len_arg = len(arg)
        args = removeParens(arg).split(',')
    return (args, len_arg)

def getPyArgs(data):
    "Get the string enclosed by the parens for interpretation via Python."
    pos = 0
    if data[pos] != '(' :
        args = []
        len_arg = 0
    else:
        findex = data[pos:].find(')')
        # print 'm ',data[pos:pos+findex+1]
        if findex > -1 :
            m = matchStringParensAndCRRE(data[pos:pos+findex+1])
        else:
            print 'argumentError: closing parenthesis not found'
            m = None
        arg = m.group()
        len_arg = len(arg)
        args = removeParens(arg)
    return (args, len_arg)

class includeDirective(procDirective) :
    def evaluate(self,data,pos) :
        includeData = ''
        self.startPosition = pos
        self.newPosition = pos
        # skip a space
        pos = pos + 1
        # Look for a quote-delimited string.
        if data.getItem(pos) != '"' and data.getItem(pos) != "'" :
            print 'includeDirective-error: expecting quote-delimeted string <'+data.data[pos-3:pos+3]+'>'
        else :
            # cheat on stream notion for data
            findex = data.findToEnd(pos,'\n')
            if findex > -1 :
                fileName = data.match(StringRE,pos,pos+findex).group()
            else:
                fileName = data.matchToEnd(StringRE,pos).group()
            if fileName :
                len_fileName = len(fileName)
                # add one for the skipped character
                self.newPosition = self.startPosition + len_fileName + 1
                fName = removeQuotes(fileName)
                with open(fName, 'r') as f:
                    includeData = f.read()
                # Insert the data ahead in the input stream.
                data.insert(self.newPosition,includeData)
                data.setPosition(self.newPosition)
            else:
                print 'includeDirective-error: no file found of name: <'+fileName+'>'
        return ''

proc_include = includeDirective('include')

class environment() :
    def __init__(self,input,output=''):
        self.input = dataString(input)
        self.output = dataString(output)

    def scan_and_proc(self) :
        # output = ''
        # outputPos = 0
        pos = 0
        while pos < self.input.getLength():
            if self.input.getItem(pos) == '@' :
                # processor directive
                # Note:  Need to put in lineNo counter
                # consume '@'
                pos = pos + 1
                directiveName = ''
                while ((not directiveName in procDirectives ) \
                       and ( pos < self.input.getLength() and self.input.getItem(pos) != '\n')) :
                    directiveName = directiveName + self.input.getItem(pos)
                    pos = pos + 1
                    if debug > 20 :
                        print 'scan_and_proc-directiveName: ',directiveName
                if directiveName in procDirectives :
                    # Note that side-effects are allowed in procDirectives!  Esp. w. INCLUDE.
                    p = procDirectives[directiveName]
                    self.output.append(p.evaluate(self.input,pos))
                    pos = p.getNewPosition()
                else:
                    if debug > 5 :
                        print 'scan_and_proc-debug-directiveNotFound: <'+directiveName+'>'
            else :
                self.output.append(self.input.getItem(pos))
                pos += 1
        return self.output.getData()

def pre(data = None, inFile = None, outFile = None) :
    if ( not data and not inFile ) :
        msg = 'pre: no input data found: --inFile used? For usage, please try -h.'
        print msg
        return msg
    if ( data and inFile ) :
        msg = 'pre: please only specify one of data or infile'
        print msg
        return msg
    if ( inFile ) :
        with open(inFile, 'r') as f:
            data = f.read()
    # Scan data
    # Process
    # Produce
    env = environment(data)
    result = env.scan_and_proc()

    if outFile :
        with open(outFile, 'w') as f:
            f.write(result)
        return 'Output written to '+outFile
    else :
        return result

class TestPreprocessor(unittest.TestCase) :
    
    def setUp(self):
        with open('test_include.inc','w') as f:
            f.write('<included text>\n')
        with open('test_include.txt','w') as f:
            f.write( dedent( \
                """
                This is a test of the include macro.
                The include macro will insert here: @include "test_include.inc".
                This is the last line of the file.
                """ ))
        with open('test_include_result.txt','w') as f:
            f.write( dedent( \
                """
                This is a test of the include macro.
                The include macro will insert here: <included text>
                .
                This is the last line of the file.
                """ ))
        pass

    def test_ALineOfText(self):
        self.assertEqual(pre(data="A line of 'text'.\n"),"A line of 'text'.\n")

    def test_AStringLiteral(self):
        self.assertEqual(pre(data='"AStringLiteral"'),'"AStringLiteral"')

    def test_INCLUDE(self):
        self.assertEqual( \
            pre(data='@include "test_include.inc"'), \
            "<included text>\n")

    def test_INCLUDEWithTextBefore(self):
        self.assertEqual( \
            pre(data='TextBefore@include "test_include.inc"'), \
            "TextBefore<included text>\n")

    def test_INCLUDEWithTextAfter(self):
        self.assertEqual( \
            pre(data='@include "test_include.inc"TextAfter'), \
            "<included text>\nTextAfter")

    def test_INCLUDEWithTextBeforeAndAfter(self):
        self.assertEqual( \
            pre(data='TextBefore@include "test_include.inc"TextAfter\n'), \
            "TextBefore<included text>\nTextAfter\n")

    def test_INCLUDEUsingFiles(self):
        pre(inFile='test_include.txt',outFile='test_include_output.txt')
        with open('test_include_output.txt','r') as f:
            actualResult = f.read()
        with open('test_include_result.txt','r') as f:
            expectedResult = f.read()
        self.assertEqual( actualResult, expectedResult )

    def test_getPCArgsWith3Args(self):
        self.assertEqual(getPCArgs('(a,b,c)'),(['a','b','c'],7))

    def test_getPCArgsWith5MixedArgs(self):
        self.assertEqual(getPCArgs('(ab,c,123,4,f)'),(['ab','c','123','4','f'],14))

    def test_getPCArgsWith5ArgsWithCRs(self):
        self.assertEqual(getPCArgs('(ab,c\n,123\n\n,4,f)'), \
                         (['ab','c\n','123\n\n','4','f'],17))

    def test_getPCArgsWith0Args(self):
        self.assertEqual(getPCArgs('()'),([''],2))

    def test_getPCArgsWith1Args(self):
        self.assertEqual(getPCArgs('(1)'),(['1'],3))
    
    # Much work to do if we need nested parens.    
    #def test_retainNestedPCArgs1(self):
    #    self.assertEqual(getPCArgs('(1,2,(3,4,(5,6)))'),())

    def test_addTokenRE1(self):
        m = proc_include.addTokenRE({':EndToken':'<EndToken>'},\
                                                 ':EndToken',\
                                                 '<DefaultEndToken>')\
                                                 .search('Testing\\n<EndToken>\\n')
        if m : 
            mSpan = m.span()
        else :
            mSpan = 'None'
        self.assertEqual(mSpan,(9,19))

    def test_addTokenRE2(self):
        m = proc_include.addTokenRE({':NotEndToken':'<NotEndTokenToken>'},\
                                                 ':EndToken',\
                                                 '<DefaultEndToken>')\
                                                 .search('Testing\\n<DefaultEndToken>\\n')
        if m : 
            mSpan = m.span()
        else :
            mSpan = None
        self.assertEqual(mSpan,(9,26))

    def test_addTokenRE3(self):
        m = proc_include.addTokenRE({':EndToken':'<EndToken>'},\
                                                 ':EndToken',\
                                                 '<DefaultEndToken>')\
                                                 .search('Testing\\n<DefaultToken>\\n')
        if m : 
            mSpan = m.span()
        else :
            mSpan = None
        self.assertEqual(mSpan,None)

    def test_addTokenRE4(self):
        RE1 = proc_include.addTokenRE({':EndToken':'<EndToken>', \
                                     ':BeginToken':'<BeginToken>'},\
                                                 ':EndToken',\
                                                 '<DefaultEndToken>')

        RE2 = proc_include.addTokenRE({':EndToken':'<EndToken>', \
                                       ':BeginToken':'<BeginToken>'}, \
                                                ':DefaultBeginToken',\
                                                 '<BeginToken>')
        
        testString = 'Testing\\n<BeginToken>\\nMiddle\\nTesting\\n<EndToken>\\n'
        m1 = RE1.search(testString); m2 = RE2.search(testString)
        if m1 :
            s1 = m1.span()
        else :
            s1 = None
            
        if m2 :
            s2 = m2.span()
        else :
            s2 = None

        self.assertEqual([s1,s2],[(40, 50),(9, 21)])
        #        self.assertEqual(str(s1)+''+str(s2),'(40, 50)(9, 21)')

    #def test_addTokenRE5(self):
    #    "Test for abilit to handle nested parens."
    #    RE1 = proc_include.addTokenRE({':EndToken':'\)', \
    #                                   ':BeginToken':'@assert\('},\
    #                                   ':EndToken','\)')
    #                                   #':BeginToken','@assert\(')
    #    testString = '@assert(a,(b),c)'
    #    m1 = RE1.search(testString)
    #    if m1 :
    #        s1 = m1.span()
    #    else :
    #        s1 = None
    #    self.assertEqual(s1,(15,16))
                
                                       
                                       



if __name__ == '__main__' :
    unittest.main()
