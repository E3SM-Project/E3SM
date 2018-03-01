#!/usr/bin/env python
# For python 2.6-2.7
from __future__ import print_function

from os.path import *
import re
# from parseBrackets import parseBrackets
from parseDirectiveArgs import parseDirectiveArguments

class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

assertVariants = 'Fail|Equal|True|False|LessThan|LessThanOrEqual|GreaterThan|GreaterThanOrEqual'
assertVariants += '|IsMemberOf|Contains|Any|All|NotAll|None|IsPermutationOf'
assertVariants += '|ExceptionRaised|SameShape|IsNaN|IsFinite'

def cppSetLineAndFile(line, file):
    return "#line " + str(line) + ' "' + file + '"\n'

def getSubroutineName(line):
    try:
        m = re.match('\s*subroutine\s+(\w*)\s*(\\([\w\s,]*\\))?\s*(!.*)*$', line, re.IGNORECASE)
        return m.groups()[0]
    except:
        raise MyError('Improper format in declaration of test procedure.')

def parseArgsFirstRest(directiveName,line):
    """If the @-directive has more than one argument, parse into first and rest strings.
    Added for assertAssociated.
    """

    argStr = ''; 
    if directiveName != '':
        m = re.match('\s*'+directiveName+'\s*\\((.*\w.*)\\)\s*$',line,re.IGNORECASE)
        if m:
            argStr = m.groups()[0]
        else:
            return None
    else:
        argStr = line

    args = parseDirectiveArguments(argStr)

    if args == []:
        returnArgs = None
    elif len(args) == 1:
        returnArgs = [args[0]]
    else:
        returnArgs = [args[0],','.join(args[1:])]
        
    return returnArgs


def parseArgsFirstSecondRest(directiveName,line):
    """If the @-directive must have at least two arguments, parse into first, second,
    and rest strings. Added for assertAssociated.
    """
    args1 = parseArgsFirstRest(directiveName,line)

    returnArgs = None

    if args1 != None:
        if len(args1) == 1:
            returnArgs = args1
        elif len(args1) == 2:
            args2 = parseArgsFirstRest('',args1[1])
            returnArgs = [args1[0]] + args2
        elif len(args1) == 3:
            print(-999,'parseArgsFirstSecondRest::error!')
            returnArgs = None

    return returnArgs


def getSelfObjectName(line):
    m = re.match('\s*subroutine\s+\w*\s*\\(\s*(\w+)\s*(,\s*\w+\s*)*\\)\s*$', line, re.IGNORECASE)
    if m:
        return m.groups()[0]
    else:
        return m

def getTypeName(line):
    m = re.match('\s*type(.*::\s*|\s+)(\w*)\s*$', line, re.IGNORECASE)
    return m.groups()[1]
 
class Action():
    def apply(self, line):
        m = self.match(line)
        if m: self.action(m, line)
        return m

#------------------
class AtTest(Action):
    def __init__(self, parser):
        self.parser = parser
        self.keyword = '@test'

    def match(self, line):
        m = re.match('\s*'+self.keyword+'(\s*(\\(.*\\))?\s*$)', line, re.IGNORECASE)
        return m

    def action(self, m, line):
        options = re.match('\s*'+self.keyword+'\s*\\((.*)\\)\s*$', line, re.IGNORECASE)
        method = {}

        if options:

            npesOption = re.search('npes\s*=\s*\\[([0-9,\s]+)\\]', options.groups()[0], re.IGNORECASE)
            if npesOption:
                npesString = npesOption.groups()[0]
                npes = map(int, npesString.split(','))
                method['npRequests'] = npes

            #ifdef is optional
            matchIfdef = re.match('.*ifdef\s*=\s*(\w+)', options.groups()[0], re.IGNORECASE)
            if matchIfdef: 
                ifdef = matchIfdef.groups()[0]
                method['ifdef'] = ifdef

            matchIfndef = re.match('.*ifndef\s*=\s*(\w+)', options.groups()[0], re.IGNORECASE)
            if matchIfndef: 
                ifndef = matchIfndef.groups()[0]
                method['ifndef'] = ifndef

            paramOption = re.search('testParameters\s*=\s*[{](.*)[}]', options.groups()[0], re.IGNORECASE)
            if paramOption:
                paramExpr = paramOption.groups()[0]
                method['testParameters'] = paramExpr

            casesOption = re.search('cases\s*=\s*(\\[[0-9,\s]+\\])', options.groups()[0], re.IGNORECASE)
            if casesOption:
                method['cases'] = casesOption.groups()[0]


        nextLine = self.parser.nextLine()
        method['name'] = getSubroutineName(nextLine)
        # save "self" name for use with @mpiAssert
        self.parser.currentSelfObjectName = getSelfObjectName(nextLine)

        # save "self" name for use with @mpiAssert
        dummyArgument = getSelfObjectName(nextLine)
        if dummyArgument:
            method['selfObjectName'] = dummyArgument

        self.parser.userTestMethods.append(method)
        self.parser.commentLine(line)
        self.parser.outputFile.write(nextLine)


#------------------
# deprecated - should now just use @test
class AtMpiTest(AtTest):
    def __init__(self, parser):
        self.parser = parser
        self.keyword = '@mpitest'

class AtTestCase(Action):
    def __init__(self, parser):
        self.parser = parser

    def match(self, line):
        m = re.match('\s*@testcase\s*(|\\(.*\\))\s*$', line, re.IGNORECASE)
        return m
    
    def action(self, m, line):
        options = re.match('\s*@testcase\s*\\((.*)\\)\s*$', line, re.IGNORECASE)
        if options:
            value = re.search('constructor\s*=\s*(\w*)', options.groups()[0], re.IGNORECASE)
            if value:
                self.parser.userTestCase['constructor'] = value.groups()[0]

            value = re.search('npes\s*=\s*\\[([0-9,\s]+)\\]', options.groups()[0], re.IGNORECASE)
            if value:
                npesString = value.groups()[0]

                npes = map(int, npesString.split(','))
                self.parser.userTestCase['npRequests'] = npes

            value = re.search('cases\s*=\s*(\\[[0-9,\s]+\\])', options.groups()[0], re.IGNORECASE)
            if value:
                cases = value.groups()[0]
                self.parser.userTestCase['cases'] = cases

            value = re.search('testParameters\s*=\s*[{](.*)[}]', options.groups()[0], re.IGNORECASE)
            if value:
                paramExpr = value.groups()[0]
                self.parser.userTestCase['testParameters'] = paramExpr

        nextLine = self.parser.nextLine()
        self.parser.userTestCase['type']=getTypeName(nextLine)
        self.parser.commentLine(line)
        self.parser.outputFile.write(nextLine)


class AtSuite(Action):
    def __init__(self, parser):
        self.parser = parser
    def match(self, line):
        nameRe = "'\w+'|" + """\w+"""
        m = re.match("\s*@suite\s*\\(\s*name\s*=\s*("+nameRe+")\s*\\)\s*$", line, re.IGNORECASE)
        return m

    def action(self, m, line):
        self.parser.suiteName=m.groups()[0][1:-1]
        self.parser.wrapModuleName = 'Wrap' + self.parser.suiteName


class AtBegin(Action):
    def __init__(self, parser):
        self.parser = parser

    def match(self, line):
        m = re.match('\s*module\s+(\w*)\s*$', line, re.IGNORECASE)
        return m

    def action(self, m, line):
        self.parser.userModuleName = m.groups()[0]
        self.parser.wrapModuleName = 'Wrap' + self.parser.userModuleName
        if not self.parser.suiteName:
            self.parser.suiteName = self.parser.userModuleName + "_suite"
        self.parser.outputFile.write(line)


class AtAssert(Action):
    def __init__(self, parser):
        self.parser = parser

    def match(self, line):
        m = re.match('\s*@assert('+assertVariants+')\s*\\((.*\w.*)\\)\s*$', line, re.IGNORECASE)
        return m

    def appendSourceLocation(self, fileHandle, fileName, lineNumber):
        fileHandle.write(" & location=SourceLocation( &\n")
        fileHandle.write(" & '" + str(basename(fileName)) + "', &\n")
        fileHandle.write(" & " + str(lineNumber) + ")")

    def action(self, m, line):
        p = self.parser
        
        p.outputFile.write(cppSetLineAndFile(p.currentLineNumber, p.fileName))
        p.outputFile.write("  call assert"+m.groups()[0]+"(" + m.groups()[1] + ", &\n")
        self.appendSourceLocation(p.outputFile, p.fileName, p.currentLineNumber)
        p.outputFile.write(" )\n")
        p.outputFile.write("  if (anyExceptions()) return\n")
        p.outputFile.write(cppSetLineAndFile(p.currentLineNumber+1, p.fileName))

class AtAssertAssociated(Action):
    def __init__(self,parser):
        self.parser = parser

    def match(self, line):
        m = re.match('\s*@assertassociated\s*\\((.*\w.*)\\)\s*$', line, re.IGNORECASE)

        if not m:
            m  = re.match( \
                '\s*@assertassociated\s*\\((\s*([^,]*\w.*),\s*([^,]*\w.*),(.*\w*.*))\\)\s*$', \
                line, re.IGNORECASE)

        # How to get both (a,b) and (a,b,c) to match?
        if not m:
            m  = re.match( \
                '\s*@assertassociated\s*\\((\s*([^,]*\w.*),\s*([^,]*\w.*))\\)\s*$', \
                line, re.IGNORECASE)
        return m

    def appendSourceLocation(self, fileHandle, fileName, lineNumber):
        fileHandle.write(" & location=SourceLocation( &\n")
        fileHandle.write(" & '" + str(basename(fileName)) + "', &\n")
        fileHandle.write(" & " + str(lineNumber) + ")")

    def action(self, m, line):
        p = self.parser

        # args = parseArgsFirstRest('@assertassociated',line)
        args = parseArgsFirstSecondRest('@assertassociated',line)

        # print(9000,line)
        # print(9001,args)
        
        p.outputFile.write(cppSetLineAndFile(p.currentLineNumber, p.fileName))
        if len(args) > 1:
            if re.match('.*message=.*',args[1],re.IGNORECASE):
                p.outputFile.write("  call assertTrue(associated(" + args[0] + "), " + args[1] + ", &\n")
            elif len(args) > 2:
                p.outputFile.write("  call assertTrue(associated(" + args[0] + "," + args[1] + "), " + args[2] + ", &\n")
            else:
                p.outputFile.write("  call assertTrue(associated(" + args[0] + "," + args[1] + "), &\n")
        else:
            p.outputFile.write("  call assertTrue(associated(" + args[0] + "), &\n")
        self.appendSourceLocation(p.outputFile, p.fileName, p.currentLineNumber)
        p.outputFile.write(" )\n")
        p.outputFile.write("  if (anyExceptions()) return\n")
        p.outputFile.write(cppSetLineAndFile(p.currentLineNumber+1, p.fileName))


class AtAssertNotAssociated(Action):
    def __init__(self,parser):
        self.parser = parser
        self.name='@assertnotassociated'

    def match(self, line):
        m = re.match('\s*@assert(not|un)associated\s*\\((.*\w.*)\\)\s*$', line, re.IGNORECASE)
        if m:
            self.name='@assert'+m.groups()[0]+'associated'
        else:
            self.name='@assertnotassociated'

        if not m:
            m  = re.match( \
                '\s*@assert(not|un)associated\s*\\((\s*([^,]*\w.*),\s*([^,]*\w.*),(.*\w*.*))\\)\s*$', \
                line, re.IGNORECASE)

        # How to get both (a,b) and (a,b,c) to match?
        if not m:
            m  = re.match( \
                '\s*@assert(not|un)associated\s*\\((\s*([^,]*\w.*),\s*([^,]*\w.*))\\)\s*$', \
                line, re.IGNORECASE)

        if m:
            self.name='@assert'+m.groups()[0]+'associated'
        else:
            self.name='@assertnotassociated'

                                
        return m

    def appendSourceLocation(self, fileHandle, fileName, lineNumber):
        fileHandle.write(" & location=SourceLocation( &\n")
        fileHandle.write(" & '" + str(basename(fileName)) + "', &\n")
        fileHandle.write(" & " + str(lineNumber) + ")")

    def action(self, m, line):
        p = self.parser

        #-- args = parseArgsFirstRest('@assertassociated',line)
        #ok args = parseArgsFirstSecondRest('@assertassociated',line)
        args = parseArgsFirstSecondRest(self.name,line)

        # print(9000,line)
        # print(9001,args)
        
        p.outputFile.write(cppSetLineAndFile(p.currentLineNumber, p.fileName))
        if len(args) > 1:
            if re.match('.*message=.*',args[1],re.IGNORECASE):
                p.outputFile.write("  call assertFalse(associated(" + args[0] + "), " + args[1] + ", &\n")
            elif len(args) > 2:
                p.outputFile.write("  call assertFalse(associated(" + args[0] + "," + args[1] + "), " + args[2] + ", &\n")
            else:
                p.outputFile.write("  call assertFalse(associated(" + args[0] + "," + args[1] + "), &\n")
        else:
            p.outputFile.write("  call assertFalse(associated(" + args[0] + "), &\n")
        self.appendSourceLocation(p.outputFile, p.fileName, p.currentLineNumber)
        p.outputFile.write(" )\n")
        p.outputFile.write("  if (anyExceptions()) return\n")
        p.outputFile.write(cppSetLineAndFile(p.currentLineNumber+1, p.fileName))


class AtAssertEqualUserDefined(Action):
    """Convenience directive replacing (a,b) with a call to assertTrue(a==b)
    and an error message, if none is provided when invoked.
    """
    def __init__(self,parser):
        self.parser = parser

    def match(self, line):
        m  = re.match( \
            '\s*@assertequaluserdefined\s*\\((\s*([^,]*\w.*),\s*([^,]*\w.*),(.*\w*.*))\\)\s*$', \
            line, re.IGNORECASE)

        # How to get both (a,b) and (a,b,c) to match?
        if not m:
            m  = re.match( \
                '\s*@assertequaluserdefined\s*\\((\s*([^,]*\w.*),\s*([^,]*\w.*))\\)\s*$', \
                line, re.IGNORECASE)
                    
        return m

    def appendSourceLocation(self, fileHandle, fileName, lineNumber):
        fileHandle.write(" & location=SourceLocation( &\n")
        fileHandle.write(" & '" + str(basename(fileName)) + "', &\n")
        fileHandle.write(" & " + str(lineNumber) + ")")

    def action(self, m, line):
        p = self.parser

        args = parseArgsFirstSecondRest('@assertequaluserdefined',line)
        
        p.outputFile.write(cppSetLineAndFile(p.currentLineNumber, p.fileName))
        if len(args) > 2:
            p.outputFile.write("  call assertTrue(" \
                               + args[0] + "==" + args[1] + ", " + args[2] + ", &\n")
        else:
            p.outputFile.write("  call assertTrue(" \
                               + args[0] + "==" + args[1] + ", &\n")
        if not re.match('.*message=.*',line,re.IGNORECASE):
            p.outputFile.write(" & message='<" + args[0] + "> not equal to <" + args[1] + ">', &\n")
        self.appendSourceLocation(p.outputFile, p.fileName, p.currentLineNumber)
        p.outputFile.write(" )\n")
        p.outputFile.write("  if (anyExceptions()) return\n")
        p.outputFile.write(cppSetLineAndFile(p.currentLineNumber+1, p.fileName))


class AtAssertEquivalent(Action):
    """Convenience directive replacing (a,b) with a call to assertTrue(a.eqv.b)
    and an error message, if none is provided when invoked.
    """
    def __init__(self,parser):
        self.parser = parser

    def match(self, line):
        m  = re.match( \
            '\s*@assertequivalent\s*\\((\s*([^,]*\w.*),\s*([^,]*\w.*),(.*\w*.*))\\)\s*$', \
            line, re.IGNORECASE)

        # How to get both (a,b) and (a,b,c) to match?
        if not m:
            m  = re.match( \
                '\s*@assertequivalent\s*\\((\s*([^,]*\w.*),\s*([^,]*\w.*))\\)\s*$', \
                line, re.IGNORECASE)
                    
        return m

    def appendSourceLocation(self, fileHandle, fileName, lineNumber):
        fileHandle.write(" & location=SourceLocation( &\n")
        fileHandle.write(" & '" + str(basename(fileName)) + "', &\n")
        fileHandle.write(" & " + str(lineNumber) + ")")

    def action(self, m, line):
        p = self.parser

        args = parseArgsFirstSecondRest('@assertequivalent',line)
        
        p.outputFile.write(cppSetLineAndFile(p.currentLineNumber, p.fileName))
        if len(args) > 2:
            p.outputFile.write("  call assertTrue(" \
                               + args[0] + ".eqv." + args[1] + ", " + args[2] + ", &\n")
        else:
            p.outputFile.write("  call assertTrue(" \
                               + args[0] + ".eqv." + args[1] + ", &\n")
        if not re.match('.*message=.*',line,re.IGNORECASE):
            p.outputFile.write(" & message='<" + args[0] + "> not equal to <" + args[1] + ">', &\n")
        self.appendSourceLocation(p.outputFile, p.fileName, p.currentLineNumber)
        p.outputFile.write(" )\n")
        p.outputFile.write("  if (anyExceptions()) return\n")
        p.outputFile.write(cppSetLineAndFile(p.currentLineNumber+1, p.fileName))
                            

class AtMpiAssert(Action):
    def __init__(self, parser):
        self.parser = parser

    def match(self, line):
        m = re.match('\s*@mpiassert('+assertVariants+')\s*\\((.*\w.*)\\)\s*$', line, re.IGNORECASE)
        return m

    def appendSourceLocation(self, fileHandle, fileName, lineNumber):
        fileHandle.write(" & location=SourceLocation( &\n")
        fileHandle.write(" & '" + str(basename(fileName)) + "', &\n")
        fileHandle.write(" & " + str(lineNumber) + ")")

    def action(self, m, line):
        p = self.parser
        
        p.outputFile.write(cppSetLineAndFile(p.currentLineNumber, p.fileName))
        p.outputFile.write("  call assert"+m.groups()[0]+"(" + m.groups()[1] + ", &\n")
        self.appendSourceLocation(p.outputFile, p.fileName, p.currentLineNumber)
        p.outputFile.write(" )\n")

        # 'this' object may not exist if test is commented out.
        if hasattr(p,'currentSelfObjectName'):
            p.outputFile.write("  if (anyExceptions("+p.currentSelfObjectName+"%context)) return\n")
        p.outputFile.write(cppSetLineAndFile(p.currentLineNumber+1, p.fileName))

class AtBefore(Action):
    def __init__(self, parser):
        self.parser = parser

    def match(self, line):
        m = re.match('\s*@before\s*$', line, re.IGNORECASE)
        return m 

    def action(self, m, line):
        nextLine = self.parser.nextLine()
        self.parser.userTestCase['setUp'] = getSubroutineName(nextLine)
        self.parser.commentLine(line)
        self.parser.outputFile.write(nextLine)

class AtAfter(Action):
    def __init__(self, parser):
        self.parser = parser

    def match(self, line):
        m = re.match('\s*@after\s*$', line, re.IGNORECASE)
        return m 

    def action(self, m, line):
        nextLine = self.parser.nextLine()
        self.parser.userTestCase['tearDown'] = getSubroutineName(nextLine)
        self.parser.commentLine(line)
        self.parser.outputFile.write(nextLine)

class AtTestParameter(Action):
    def __init__(self, parser):
        self.parser = parser

    def match(self, line):
        m = re.match('\s*@testParameter\s*(|.*)$', line, re.IGNORECASE)
        return m

    def action(self, m, line):
        options = re.match('\s*@testParameter\s*\\((.*)\\)\s*$', line, re.IGNORECASE)

        self.parser.commentLine(line)
        nextLine = self.parser.nextLine()
        if not 'testParameterType' in self.parser.userTestCase:
            self.parser.userTestCase['testParameterType'] = getTypeName(nextLine)
        self.parser.outputFile.write(nextLine)

        if options:
            value = re.search('constructor\s*=\s*(\w*)', options.groups()[0], re.IGNORECASE)
            if value:
                self.parser.userTestCase['testParameterConstructor'] = value.groups()[0]
            else:
                self.parser.userTestCase['testParameterConstructor'] = self.parser.userTestCase['testParameterType']


class Parser():
    def __init__(self, inputFileName, outputFileName):
        def getBaseName(fileName):
            from os.path import basename, splitext
            base = basename(fileName)
            return splitext(base)[0]

        self.fileName = inputFileName
        self.inputFile = open(inputFileName, 'r')
        self.outputFile = open(outputFileName, 'w')
        self.defaultSuiteName = getBaseName(inputFileName) + "_suite"
        self.suiteName = ''

        self.currentLineNumber = 0
        self.userModuleName = '' # if any

        self.userTestCase = {}
        self.userTestCase['setUpMethod'] = ''
        self.userTestCase['tearDownMethod'] = ''
        self.userTestCase['defaultTestParameterNpes'] = [] # is MPI if not empty
        self.userTestCase['defaultTestParametersExpr'] = ''
        self.userTestCase['defaultTestParameterCases'] = [] 

        self.userTestMethods = [] # each entry is a dictionary

        self.wrapModuleName = "Wrap" + getBaseName(inputFileName)
        self.currentLineNumber = 0

        self.actions=[]
        self.actions.append(AtTest(self))
        self.actions.append(AtMpiTest(self))
        self.actions.append(AtTestCase(self))
        self.actions.append(AtSuite(self))
        self.actions.append(AtBegin(self))

        self.actions.append(AtAssert(self))
        self.actions.append(AtAssertAssociated(self))
#        self.actions.append(AtAssertAssociatedWith(self))
        self.actions.append(AtAssertNotAssociated(self))
#        self.actions.append(AtAssertNotAssociatedWith(self))

        self.actions.append(AtAssertEqualUserDefined(self))
        self.actions.append(AtAssertEquivalent(self))
        
        self.actions.append(AtMpiAssert(self))
        self.actions.append(AtBefore(self))
        self.actions.append(AtAfter(self))
        self.actions.append(AtTestParameter(self))


    def commentLine(self, line):
        self.outputFile.write(re.sub('@','!@',line))

    def run(self):
        def parse(line):
            for action in self.actions:
                if (action.apply(line)): return
            self.outputFile.write(line)

        while True:
            line = self.nextLine()
            if  not line: break
            parse(line)

        if (not self.suiteName): self.suiteName = self.defaultSuiteName
        if ('testParameterType' in self.userTestCase and (not 'constructor' in self.userTestCase)):
            self.userTestCase['constructor'] = self.userTestCase['testParameterType']
        self.makeWrapperModule()

    def isComment(self, line):
        return re.match('\s*(!.*|)$', line)

    def nextLine(self):
        while True:
            self.currentLineNumber += 1
            line = self.inputFile.readline()
            if not line: break
            if (self.isComment(line)):
                self.outputFile.write(line)
                pass
            else:
                break
        return line


    def printHeader(self):
        self.outputFile.write('\n')
        self.outputFile.write('module ' + self.wrapModuleName + '\n')
        self.outputFile.write('   use pFUnit_mod\n')
        if (self.userModuleName): self.outputFile.write('   use ' + self.userModuleName + '\n')
        self.outputFile.write('   implicit none\n')
        self.outputFile.write('   private\n\n')



    def printTail(self):
        self.outputFile.write('\n')
        self.outputFile.write('end module ' + self.wrapModuleName + '\n\n')

    def printWrapUserTestCase(self):
        self.outputFile.write('   public :: WrapUserTestCase\n')
        self.outputFile.write('   public :: makeCustomTest\n')
        self.outputFile.write('   type, extends(' + self.userTestCase['type'] + ') :: WrapUserTestCase\n')
        self.outputFile.write('      procedure(userTestMethod), nopass, pointer :: testMethodPtr\n')
        self.outputFile.write('   contains\n')
        self.outputFile.write('      procedure :: runMethod\n')
        self.outputFile.write('   end type WrapUserTestCase\n\n')
        
        self.outputFile.write('   abstract interface\n')
        self.outputFile.write('     subroutine userTestMethod(this)\n')
        if self.userModuleName:
            self.outputFile.write('        use ' + self.userModuleName + '\n')
        if 'type' in self.userTestCase:
            self.outputFile.write('        class (' + self.userTestCase['type'] + '), intent(inout) :: this\n')
        self.outputFile.write('     end subroutine userTestMethod\n')
        self.outputFile.write('   end interface\n\n')

    def printRunMethod(self):
        self.outputFile.write('   subroutine runMethod(this)\n')
        self.outputFile.write('      class (WrapUserTestCase), intent(inout) :: this\n\n')
        self.outputFile.write('      call this%testMethodPtr(this)\n')
        self.outputFile.write('   end subroutine runMethod\n\n')

            
    def printParameterHeader(self, type):
        self.outputFile.write('   type (' + type + '), allocatable :: testParameters(:)\n')
        self.outputFile.write('   type (' + type + ') :: testParameter\n')
        self.outputFile.write('   integer :: iParam \n')
        self.outputFile.write('   integer, allocatable :: cases(:) \n')
        self.outputFile.write(' \n')


    def printMakeSuite(self):
        self.outputFile.write('function ' + self.suiteName + '() result(suite)\n')
        self.outputFile.write('   use pFUnit_mod\n')
        if (self.userModuleName): self.outputFile.write('   use ' + self.userModuleName + '\n')
        self.outputFile.write('   use '+ self.wrapModuleName + '\n')
        self.outputFile.write('   type (TestSuite) :: suite\n\n')

        if not self.userModuleName:
            for testMethod in self.userTestMethods:
                if ('ifdef' in testMethod):
                    self.outputFile.write('#ifdef ' + testMethod['ifdef'] + '\n')
                elif ('ifndef' in testMethod):
                    self.outputFile.write('#ifndef ' + testMethod['ifndef'] + '\n')
                self.outputFile.write('   external ' + testMethod['name'] + '\n')
                if ('ifdef' in testMethod or 'ifndef' in testMethod):
                    self.outputFile.write('#endif\n')
            self.outputFile.write('\n')
            if 'setUp' in self.userTestCase:
                self.outputFile.write('   external ' + self.userTestCase['setUp'] + '\n')
            if 'tearDown' in self.userTestCase:
                self.outputFile.write('   external ' + self.userTestCase['tearDown'] + '\n')
            self.outputFile.write('\n')

        self.outputFile.write('   integer, allocatable :: npes(:)\n\n')

        if 'testParameterType' in self.userTestCase:
            type = self.userTestCase['testParameterType']
            self.printParameterHeader(type)

        self.outputFile.write("   suite = newTestSuite('" + self.suiteName + "')\n\n")

        for testMethod in self.userTestMethods:
            if ('ifdef' in testMethod):
                self.outputFile.write('#ifdef ' + testMethod['ifdef'] + '\n')
            elif ('ifndef' in testMethod):
                self.outputFile.write('#ifndef ' + testMethod['ifndef'] + '\n')
            if 'type' in self.userTestCase:
                self.addUserTestMethod(testMethod)
            else:
                if 'npRequests' in testMethod:
                    self.addMpiTestMethod(testMethod)
                else: # vanilla
                    self.addSimpleTestMethod(testMethod)
            self.outputFile.write('\n')
            if ('ifdef' in testMethod or 'ifndef' in testMethod):
                self.outputFile.write('#endif\n')

        self.outputFile.write('\nend function ' + self.suiteName + '\n\n')

    def addSimpleTestMethod(self, testMethod):
        args = "'" + testMethod['name'] + "', " + testMethod['name']
        if 'setUp' in testMethod:
            args += ', ' + testMethod['setUp']
        elif 'setUp' in self.userTestCase:
            args += ', ' + self.userTestCase['setUp']

        if 'tearDown' in testMethod:
            args += ', ' + testMethod['tearDown']
        elif 'tearDown' in self.userTestCase:
            args += ', ' + self.userTestCase['tearDown']

        self.outputFile.write('   call suite%addTest(newTestMethod(' + args + '))\n')

    def addMpiTestMethod(self, testMethod):
        for npes in testMethod['npRequests']:
            args = "'" + testMethod['name'] + "', " + testMethod['name'] + ", " + str(npes)
            if 'setUp' in testMethod:
                args += ', ' + testMethod['setUp']
            elif 'setUp' in self.userTestCase:
                args += ', ' + self.userTestCase['setUp']

            if 'tearDown' in testMethod:
                args += ', ' + testMethod['tearDown']
            elif 'tearDown' in self.userTestCase:
                args += ', ' + self.userTestCase['tearDown']

            self.outputFile.write('   call suite%addTest(newMpiTestMethod(' + args + '))\n')
    
    def addUserTestMethod(self, testMethod):

        args = "'" + testMethod['name'] + "', " + testMethod['name']
        if 'npRequests' in testMethod:
            npRequests = testMethod['npRequests']
        else:
            if 'npRequests' in self.userTestCase:
                npRequests = self.userTestCase['npRequests']
            else:
                npRequests = [1]

        if 'cases' in testMethod:
            cases = testMethod['cases']
        elif 'cases' in self.userTestCase:
            cases = self.userTestCase['cases']

        testParameterArg = '' # unless

        if 'cases' in locals():
            testParameterArg = ', testParameter'
            self.outputFile.write('   cases = ' + testMethod['cases'] + '\n')
            self.outputFile.write('   testParameters = [(' + 
                                  self.userTestCase['testParameterConstructor'] + 
                                  '(cases(iCase)), iCase = 1, size(cases))]\n\n')

        if 'testParameterType' in self.userTestCase:
            if 'testParameters' in testMethod:
                testParameters = testMethod['testParameters']
            elif 'testParameters' in self.userTestCase:
                testParameters = self.userTestCase['testParameters']

        isMpiTestCase = 'npRequests' in self.userTestCase
        isMpiTestCase = isMpiTestCase or any('npRequests' in testMethod for testMethod in self.userTestMethods)

        if 'testParameters' in locals():
            testParameterArg = ', testParameter'
            self.outputFile.write('   testParameters = ' + testParameters + '\n\n')
        elif isMpiTestCase:
            testParameterArg = ', testParameter'
        

        for npes in npRequests:

            if 'testParameters' in locals() or 'cases' in locals():
                self.outputFile.write('   do iParam = 1, size(testParameters)\n')
                self.outputFile.write('      testParameter = testParameters(iParam)\n')

            if isMpiTestCase:
                self.outputFile.write('   call testParameter%setNumProcessesRequested(' + str(npes) + ')\n')

            self.outputFile.write('   call suite%addTest(makeCustomTest(' + 
                                  args + testParameterArg + '))\n')
            if 'cases' in locals() or 'testParameters' in locals():
                self.outputFile.write('   end do\n')

                

    def printMakeCustomTest(self, isMpiTestCase):
        args = 'methodName, testMethod'
        declareArgs =  '#ifdef INTEL_13\n'
        declareArgs +=  '      use pfunit_mod, only: testCase\n'
        declareArgs +=  '#endif\n'
        declareArgs +=  '      type (WrapUserTestCase) :: aTest\n'
        declareArgs +=  '#ifdef INTEL_13\n'
        declareArgs +=  '      target :: aTest\n'
        declareArgs +=  '      class (WrapUserTestCase), pointer :: p\n'
        declareArgs +=  '#endif\n'
        declareArgs += '      character(len=*), intent(in) :: methodName\n'
        declareArgs += '      procedure(userTestMethod) :: testMethod\n'
        
        if 'testParameterType' in self.userTestCase:
            args += ', testParameter'
            declareArgs += '      type (' + self.userTestCase['testParameterType'] + '), intent(in) :: testParameter\n'

        self.outputFile.write('   function makeCustomTest(' + args + ') result(aTest)\n')
        self.outputFile.write(declareArgs)

        if 'constructor' in self.userTestCase:
            if 'testParameterType' in self.userTestCase:
                constructor = self.userTestCase['constructor'] + '(testParameter)'
            else:
                constructor = self.userTestCase['constructor'] + '()'
            self.outputFile.write('      aTest%' + self.userTestCase['type'] + ' = ' + constructor + '\n\n')

        self.outputFile.write('      aTest%testMethodPtr => testMethod\n')
        
        self.outputFile.write('#ifdef INTEL_13\n')
        self.outputFile.write('      p => aTest\n')
        self.outputFile.write('      call p%setName(methodName)\n')
        self.outputFile.write('#else\n')
        self.outputFile.write('      call aTest%setName(methodName)\n')
        self.outputFile.write('#endif\n')

        if 'testParameterType' in self.userTestCase:
            self.outputFile.write('      call aTest%setTestParameter(testParameter)\n')
        

        self.outputFile.write('   end function makeCustomTest\n')

    def makeWrapperModule(self):
        #-----------------------------------------------------------
        # ! Start here
        self.printHeader()

        if 'type' in self.userTestCase:
            self.printWrapUserTestCase()
        
        self.outputFile.write('contains\n\n')

        if 'type' in self.userTestCase:
            self.printRunMethod()

        if 'type' in self.userTestCase:
            isMpiTestCase = 'npRequests' in self.userTestCase
            isMpiTestCase = isMpiTestCase or any('npRequests' in testMethod for testMethod in self.userTestMethods)
            if isMpiTestCase and not 'testParameterType' in self.userTestCase:
                self.userTestCase['testParameterType'] = 'MpiTestParameter'

            self.printMakeCustomTest(isMpiTestCase)

        self.printTail()
        self.printMakeSuite()

    def final(self):
        self.inputFile.close()
        self.outputFile.close()

if __name__ == "__main__":
    import sys
    print("Processing file", sys.argv[1])
    p = Parser(sys.argv[1], sys.argv[2])
    p.run()
    p.final()
    print(" ... Done.  Results in", sys.argv[2])


