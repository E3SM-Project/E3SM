#!/usr/bin/env python

import re
from parseArgs import parseArgs
from textwrap import dedent
from pre2 import *

debug = 0
debugLevel = 50



StringPyArg = re.compile(r'(.+=.+,?)+')
def matchStringPyArg(s) :
    return StringPyArg.findall(s)

class RepeatDirective(procDirective):
    def evaluate(self,data,pos):
        self.startPosition = pos; self.newPosition = pos

        # Args offset by parens
        argsAndLen = getPyArgs(data.getSliceForward(pos))
        args = argsAndLen[0]; lengthOfOriginalArgs = argsAndLen[1]

        args = parseArgs(args)
        
        if debug > debugLevel :
            print 'cs: ',args
        replaceString = args.keys()
        if debug > debugLevel :
            print 'cs-rs0: ',replaceString
        
        # Experiment with putting new text back onto the input...
        self.newPosition = pos + lengthOfOriginalArgs
        data.setPosition(self.newPosition)

        newText = ''; src = ''

        # Search for the end of the subroutine to determine how much to copy.
        # For the future:  How to deal with nested?  Maybe understand "contains"?
        # Also:  Consider the case of FUNCTION too.

        # allow caller to set end of block, but default to 'end subroutine'

        untilRE = self.addTokenRE(args,':until',r'end subroutine')
        
#mlr-        if not ':until' in replaceString :
#mlr-            untilEndToken = r'end subroutine'
#mlr-        else :
#mlr-            untilEndToken = args[':until'][:]
#mlr-            del args[':until']
#mlr-        RE=r'''(?i)[ \t]*'''+untilEndToken
#mlr-
#mlr-        if debug > debugLevel :
#mlr-            print 'args: ',args
#mlr-            print 'RE: ',RE
#mlr-
#mlr-        untilRE = re.compile(RE)
        
        m = data.searchToEnd(untilRE,self.newPosition)

        if not m :
            print 'Repeat syntax error:  untilEndToken: "'+untilEndToken+'" not found'
            return '<Repeat syntax error:  untilEndToken: "'+untilEndToken+'" not found>'

        findex = m.start()
        blockPos = self.newPosition
        blockPosEOL = blockPos + findex
        while data.getItem(blockPosEOL) != '\n' :
            blockPosEOL += 1
            if blockPosEOL > data.getLength() :
                print 'Repeat no EOL on end of subroutine!'
                break
        # copy a slice...

        self.newPosition = blockPosEOL
        if debug > debugLevel :
            print 'cs: rs  ',replaceString
            print 'cs: bp. ',blockPos,blockPosEOL

        # Copy the data
        origSrc = data.getSlice(blockPos,blockPosEOL)

        # When you reach a new case, make enough copies for the substitution.
        srcList = [origSrc[:]]
        for k in args.keys():
            p = '<'+k+'>'
            l = args[k]
            if not type(l) is list:
                l = [l]
            tmpList = []
            for s in srcList:
                for v in l:
                    # Make the substitution
                    tmpList.append(re.sub(p,str(v),s[:]))
            srcList = tmpList
        
        # Bring the copies together
        newText = newText + ''.join(srcList)

        # Experiment with putting data onto the input stream
        data.insert(self.newPosition,newText)
        return ''

        # Put onto the output
        # return newText

# Instantiate - also loads into the interpreter/preprocessor main loop.
proc_RepeatDirective = RepeatDirective('Repeat')

# Test...
class TestRepeatDirective(unittest.TestCase) :

    def test_copyBlock1(self) :
        self.assertEqual( \
            pre(data = dedent( \
                '''
                @Repeat(atype=1)
                ! A comment.
                subroutine foo<atype>(a)
                   Something(<atype>)
                end subroutine
                ''' )), dedent ( \
                '''
                
                ! A comment.
                subroutine foo1(a)
                   Something(1)
                end subroutine
                ''' ))

    def test_copyBlock2(self) :
        self.assertEqual( \
            pre(data = dedent( \
                '''
                @Repeat(atype=[1,2],:until="end subroutine")
                ! A comment.
                subroutine foo<atype>(a)
                   Something(<atype>)
                end subroutine
                ''' )), dedent ( \
                '''
                
                ! A comment.
                subroutine foo1(a)
                   Something(1)
                end subroutine
                ! A comment.
                subroutine foo2(a)
                   Something(2)
                end subroutine
                ''' ))

    def test_copyBlock2Vars(self) :
        self.assertEqual( \
            pre(data = dedent( \
                '''
                @Repeat(atype=1,btype=2)
                ! A comment.
                subroutine foo<atype><btype>(a)
                   Something(<atype>)
                   Something(<btype>)
                end subroutine
                ''' )), dedent ( \
                '''
                
                ! A comment.
                subroutine foo12(a)
                   Something(1)
                   Something(2)
                end subroutine
                ''' ))

    def test_copyBlock2VarsMulti(self) :
        self.assertEqual( \
            pre(data = dedent( \
                '''
                @Repeat(atype=[1,2],btype=[3,4,5],ctype=[6])
                ! A comment.
                subroutine foo<atype><btype><ctype>(a)
                   Something(<atype>)
                   Something(<btype>)
                   Something(<ctype>)
                end subroutine
                ''' )), dedent ( \
                '''
                
                ! A comment.
                subroutine foo136(a)
                   Something(1)
                   Something(3)
                   Something(6)
                end subroutine
                ! A comment.
                subroutine foo146(a)
                   Something(1)
                   Something(4)
                   Something(6)
                end subroutine
                ! A comment.
                subroutine foo156(a)
                   Something(1)
                   Something(5)
                   Something(6)
                end subroutine
                ! A comment.
                subroutine foo236(a)
                   Something(2)
                   Something(3)
                   Something(6)
                end subroutine
                ! A comment.
                subroutine foo246(a)
                   Something(2)
                   Something(4)
                   Something(6)
                end subroutine
                ! A comment.
                subroutine foo256(a)
                   Something(2)
                   Something(5)
                   Something(6)
                end subroutine
                ''' ))

    def test_copyBlock2VarsMultiWithStrings(self) :
        self.assertEqual( \
            pre(data = dedent( \
                '''
                @Repeat(atype=['Alt',"Pyg"],btype=["Orange",'Pear','Apple'],ctype=[6])
                ! A comment.
                subroutine foo<atype><btype><ctype>(a)
                   Something(<atype>)
                   Something(<btype>)
                   Something(<ctype>)
                end subroutine
                ''' )), dedent ( \
                '''
                
                ! A comment.
                subroutine fooAltOrange6(a)
                   Something(Alt)
                   Something(Orange)
                   Something(6)
                end subroutine
                ! A comment.
                subroutine fooAltPear6(a)
                   Something(Alt)
                   Something(Pear)
                   Something(6)
                end subroutine
                ! A comment.
                subroutine fooAltApple6(a)
                   Something(Alt)
                   Something(Apple)
                   Something(6)
                end subroutine
                ! A comment.
                subroutine fooPygOrange6(a)
                   Something(Pyg)
                   Something(Orange)
                   Something(6)
                end subroutine
                ! A comment.
                subroutine fooPygPear6(a)
                   Something(Pyg)
                   Something(Pear)
                   Something(6)
                end subroutine
                ! A comment.
                subroutine fooPygApple6(a)
                   Something(Pyg)
                   Something(Apple)
                   Something(6)
                end subroutine
                ''' ))

    def test_copyNaiveRecursion(self) :
        self.assertEqual( \
            pre(data = dedent( \
                '''
                @Repeat(atype=[1,2,3])
                @Repeat(z=['a','b','c'])
                subroutine NaiveRecurse<atype>(x)
                   Something(<z><atype>)
                   [<z>,<atype>]-[<atype>,<z>]
                end subroutine
                ''')), dedent( \
                '''
                
                
                subroutine NaiveRecurse1(x)
                   Something(a1)
                   [a,1]-[1,a]
                end subroutine
                subroutine NaiveRecurse1(x)
                   Something(b1)
                   [b,1]-[1,b]
                end subroutine
                subroutine NaiveRecurse1(x)
                   Something(c1)
                   [c,1]-[1,c]
                end subroutine
                
                subroutine NaiveRecurse2(x)
                   Something(a2)
                   [a,2]-[2,a]
                end subroutine
                subroutine NaiveRecurse2(x)
                   Something(b2)
                   [b,2]-[2,b]
                end subroutine
                subroutine NaiveRecurse2(x)
                   Something(c2)
                   [c,2]-[2,c]
                end subroutine
                
                subroutine NaiveRecurse3(x)
                   Something(a3)
                   [a,3]-[3,a]
                end subroutine
                subroutine NaiveRecurse3(x)
                   Something(b3)
                   [b,3]-[3,b]
                end subroutine
                subroutine NaiveRecurse3(x)
                   Something(c3)
                   [c,3]-[3,c]
                end subroutine
                '''))

    def test_copyNaiveRecursion1(self) :
        self.assertEqual( \
            pre(data = dedent( \
                '''
                @Repeat(type=['i','r'])
                @Repeat(kind=['32','64++<tol>'])
                @Repeat(tol=['0','32'])
                subroutine test_<type><kind>_tol<tol>(x,y,z)
                end subroutine
                ''' )), dedent( \
                '''                
                

                
                subroutine test_i32_tol0(x,y,z)
                end subroutine
                subroutine test_i32_tol32(x,y,z)
                end subroutine
                
                subroutine test_i64++0_tol0(x,y,z)
                end subroutine
                subroutine test_i64++32_tol32(x,y,z)
                end subroutine
                
                
                subroutine test_r32_tol0(x,y,z)
                end subroutine
                subroutine test_r32_tol32(x,y,z)
                end subroutine
                
                subroutine test_r64++0_tol0(x,y,z)
                end subroutine
                subroutine test_r64++32_tol32(x,y,z)
                end subroutine
                ''' ))

    def test_copyFunction1(self) :
        self.assertEqual( \
            pre(data = dedent( \
                '''
                @Repeat(atype=1,:until="end function")
                ! A comment.
                function foo<atype>(a)
                   Something(<atype>)
                end function
                ''' )), dedent ( \
                '''
                
                ! A comment.
                function foo1(a)
                   Something(1)
                end function
                ''' ))


if __name__ == '__main__' :
    unittest.main()
        
        
