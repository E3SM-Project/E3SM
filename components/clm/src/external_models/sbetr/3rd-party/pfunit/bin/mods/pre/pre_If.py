#!/usr/bin/env python

import re
from parseArgs import parseArgs
from textwrap import dedent
from pre2 import *
from interleavedp import interleavedp

debug = 00
debugLevel = 50

class interval :
    def __init__(self,start,end) :
        self.start = start
        self.end = end
        self.interval = [start,end]
    def getInterval(self):
        return self.interval
    def setInterval(self,start,end) :
        self.start = start
        self.end = end
        self.interval = [start,end]
    def getStart(self):
        return self.start
    def getEnd(self):
        return self.end

class IfDirective(procDirective) :
    def evaluate(self, data, pos) :
        self.startPosition = pos; self.newPosition = pos


        # Args offset by parens
        argsAndLen = getPyArgs(data.getSliceForward(pos))
        args = argsAndLen[0]; lengthOfOriginalArgs = argsAndLen[1]
        args = parseArgs(args)
        
        if debug > debugLevel :
            print 'cs: ',args
            print 'cs: pos: ',pos
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

        key = ':BeginIfToken'; defaultToken = '<BeginIf>'
        beginRE = self.addTokenRE(args,key,defaultToken,prefix='')

        key = ':ElseToken'; defaultToken = '<Else>'
        elseRE = self.addTokenRE(args,key,defaultToken,prefix='')

        key = ":EndIfToken"; defaultToken = '<EndIf>'
        endifRE = self.addTokenRE(args,key,defaultToken,prefix='')

        key = ":UntilEnd"; defaultToken = '<UntilEndIfs>'
        untilRE = self.addTokenRE(args,key,defaultToken,prefix='')

        mUntil = self.searchTokenToEnd(':UntilEnd',data,self.newPosition)
        # Or maybe findall?
        # We should signal an error if we detect more than one :UntilEnd token...
        # Get a pos and limit searches below.
        # But first, just find when the first stop is.
        if mUntil :
            # print ':UntilEnd found'
            mU_span = mUntil.span()
        else :
            if debug > debugLevel :
                print ':UntilEnd not found, searching to end'
            mU_span = [data.getLength()-1,data.getLength()-1]
        if debug > debugLevel :
            print 'mU_span: ',mU_span
            
        # until the first UntilEndifsToken
        mUntilPos = mU_span[0]
        if debug > debugLevel :
            # print 'mUntilPos: ',mUntilPos
            print 'data:    ',data.getSlice(self.newPosition,self.newPosition+mUntilPos+5)
            print 'data:    ',data.getSlice(self.newPosition+mUntilPos,self.newPosition+mUntilPos+5)

        mBeginIter = self.finditerTokenToPosition(':BeginIfToken',data,self.newPosition,self.newPosition+mUntilPos)
        if not mBeginIter :
            return self.makeTokenErrorMessage( \
                'Ifs-Syntax error... Maybe not found... ',':BeginIfToken')
        
        mEndIfIter = self.finditerTokenToPosition(':EndIfToken',data,self.newPosition,self.newPosition+mUntilPos)
        if not mEndIfIter :
            return self.makeTokenErrorMessage( \
                'Ifs-Syntax error... Maybe not found... ',':EndIfToken')

        mElseIter = self.finditerTokenToPosition(':ElseToken',data,self.newPosition,self.newPosition+mUntilPos)

        beginsStart = []; beginsEnd = []
        for i in mBeginIter :
            beginsStart.append(i.start())
            beginsEnd.append(i.end())
    
        elsesStart = []; elsesEnd = []
        if mElseIter :
            for i in mElseIter :
                elsesStart.append(i.start())
                elsesEnd.append(i.end())

        endsStart = []; endsEnd = []
        for i in mEndIfIter :
            endsStart.append(i.start())
            endsEnd.append(i.end())

        if debug > debugLevel :
            print 'be: ',beginsEnd
            print 'es: ',endsStart
            print 'ee: ',endsEnd

        pos = self.getNewPosition()
        # for i in range(len(beginsStart)) :
        #     print 'bs(',i,') : ',data.getSlice(pos+beginsStart[i],pos+beginsStart[i]+10)
        # for i in range(len(endsStart)) :
        #     print 'es(',i,') : ',data.getSlice(pos+endsStart[0],pos+endsStart[0]+10)

        ok = interleavedp(beginsEnd,endsStart,m=elsesStart)
        if not ok :
            return 'IF -- BeginIf-Else-EndIf interleaving error'

        ifIntervals = []; elseIntervals = []
        elseCount = 0
        #        for ii in range(len(beginsEnd)-1) :
        for ii in range(len(beginsEnd)) :
            beg = beginsEnd[ii]
            elsePos = None
            end = endsStart[ii]
            if elsesEnd :
                e = elsesEnd[elseCount]
                if (beg < e) and (e < end) :
                    elsePos = e #?
                    ifIntervals.append(interval(beg,elsesStart[elseCount]))
                    elseIntervals.append(interval(elsesEnd[elseCount],end))
                    elseCount = elseCount + 1
                else :
                    ifIntervals.append(interval(beg,end))
            else :
                ifIntervals.append(interval(beg,end))

        KeepIfClause = False
        if not ':test' in args.keys() :
            # print 'no :test in IF!'
            return 'no :test in IF!'
        else :
            #            print ':test = ',args[':test']
            p = eval(args[':test'])
            KeepIfClause = p
            #mlr ???            p = args[':test']
            #mlr-old KeepIfClause = p == 'True'

        if debug > debugLevel :
            print 'KeepIfClause, p: ',KeepIfClause,p
            print '0+++',data.getSliceForward(pos),'---'
        pos = self.newPosition
        if KeepIfClause :
            if debug > debugLevel :
                print 'removing else intervals...'
            offset = 0
            for i in elseIntervals :
                data.removeSlice(pos + i.getStart()+offset,pos + i.getEnd()+offset)
                offset = offset - ( i.getEnd() - i.getStart() )
        else:
            if debug > debugLevel :
                print 'removing if intervals...'
            offset = 0
            for i in ifIntervals :
                data.removeSlice(pos + i.getStart()+offset,pos + i.getEnd()+offset)
                offset = offset - ( i.getEnd() - i.getStart() )
        if debug > debugLevel :
            print '1+++',data.getSliceForward(pos),'---'

        # Remember data has been modified...
        pos = self.newPosition
        # Get the first one...
        mUntil = self.searchTokenToEnd(':UntilEnd',data,self.newPosition)
        if mUntil :
            # print ':UntilEnd found'
            mU_span = mUntil.span()
            # print 'mU_span0: ',mU_span
        else :
            # print ':UntilEnd found not found, searching to end'
            if debug > debugLevel :
                print 'mU_span1,endsEnd: ',endsEnd
            mU_span = [endsEnd[-1],endsEnd[-1]]
            # print 'mU_span1: ',mU_span
        # until the last :UntilEnd
        mUntilPos = mU_span[0]
        mUntilEnd = mU_span[1]

        clearTokens = True
        if ':keepTokens' in args.keys() :
            clearTokens = args[':keepTokens'] != 'True'
            
        if clearTokens :
            if debug > debugLevel :
                print 'clearing tokens in data from ',pos,' to ',pos + mUntilEnd
                print '+++',data.getSlice(pos,pos + mUntilEnd),'---'
            offsetEnd = 0
            for key in self.tokens.keys():
                # print 'clearing key = ',key,' = ',self.tokens[key]
                # print data.getSlice(pos,mUntilEnd+offset)
                clearIter = self.finditerTokenToPosition(key,data,pos,pos + mUntilEnd + offsetEnd )
                offset = 0
                for item in clearIter :
                    start = item.start()
                    end = item.end()
                    length = end - start
                    if debug > debugLevel :
                        print 'pos,start,length,offset',pos,start,length,offset
                        print '+++',data.getSliceForward(pos),'---'
                        print 'removing: ',data.getSlice(pos+start+offset,pos+end+offset)
                    data.removeSlice(pos+start+offset,pos+end+offset)
                    offset = offset - length
                    offsetEnd = offsetEnd - length
                    if debug > debugLevel :
                        print '+++',data.getSliceForward(pos),'---'



                                                         
        return ''
        


proc_IfDirective = IfDirective('If')

class TestIfDirective(unittest.TestCase) :

    def testTokenNotFound1(self):
        self.assertEqual(\
                         pre(data = dedent( \
                                            '''
                                            @If()
                                            Nobody but us chickens.
                                            No tokens here.
                                            ''')), dedent ( \
                                            '''
                                            no :test in IF!
                                            Nobody but us chickens.
                                            No tokens here.
                                            '''))

    def testNoTest(self):
        self.assertEqual(\
                         pre(data = dedent( \
                                            '''
                                            @If(:keepTokens)
                                            Nobody but us chickens.
                                            <BeginIf>
                                            No EndIf tokens here.
                                            ''')), dedent ( \
                                            '''
                                            IF -- BeginIf-Else-EndIf interleaving error
                                            Nobody but us chickens.
                                            <BeginIf>
                                            No EndIf tokens here.
                                            '''))

    def testIFTestFalse(self):
        self.assertEqual(\
                         pre(data = dedent( \
                                            '''
                                            @If(:test='False',:keepTokens='True')
                                            Nobody but us 
                                            <BeginIf>
                                            chickens. (If-Clause)
                                            <Else>
                                            roosters. (Else-Clause)
                                            <EndIf>
                                            ''')), dedent ( \
                                            '''

                                            Nobody but us 
                                            <BeginIf><Else>
                                            roosters. (Else-Clause)
                                            <EndIf>
                                            '''))

    def testIFTestTrue1(self):
        self.assertEqual(\
                         pre(data = dedent( \
                                            '''
                                            @If(:test='True',:keepTokens='True')
                                            Nobody but us 
                                            <BeginIf>
                                            chickens. (If-Clause)
                                            <Else>
                                            roosters. (Else-Clause)
                                            <EndIf>
                                            ''')), dedent ( \
                                            '''

                                            Nobody but us 
                                            <BeginIf>
                                            chickens. (If-Clause)
                                            <Else><EndIf>
                                            '''))


    def testIFTestTrue2(self):
        self.assertEqual(\
                         pre(data = dedent( \
                                            '''
                                            @If(:test='True',:keepTokens='True')
                                            Nobody but us 
                                            <BeginIf>
                                            chickens. (If-Clause)
                                            <Else>
                                            roosters. (Else-Clause)
                                            <EndIf>

                                            Foxes like
                                            <BeginIf>
                                            chickens. (If-Clause)
                                            <Else>
                                            roosters. (Else-Clause)
                                            <EndIf>

                                            ''')), dedent ( \
                                            '''

                                            Nobody but us 
                                            <BeginIf>
                                            chickens. (If-Clause)
                                            <Else><EndIf>

                                            Foxes like
                                            <BeginIf>
                                            chickens. (If-Clause)
                                            <Else><EndIf>

                                            '''))
        
    def testIFClearTokens(self):
        self.assertEqual(\
                         pre(data = dedent( \
                                            '''
                                            @If(:test='True')
                                            Nobody but us 
                                            <BeginIf>
                                            chickens. (If-Clause)
                                            <Else>
                                            roosters. (Else-Clause)
                                            <EndIf>

                                            Foxes like
                                            <BeginIf>
                                            chickens. (If-Clause)
                                            <Else>
                                            roosters. (Else-Clause)
                                            <EndIf>

                                            ''')), dedent ( \
                                            '''

                                            Nobody but us 

                                            chickens. (If-Clause)


                                            Foxes like

                                            chickens. (If-Clause)


                                            '''))

    def testIFClearTokensUntilEnd1(self):
        self.assertEqual(\
                         pre(data = dedent( \
                                            '''
                                            @If(:test='True')
                                            Nobody but us 
                                            <BeginIf>
                                            chickens. (If-Clause)
                                            <Else>
                                            roosters. (Else-Clause)
                                            <EndIf>

                                            Foxes like <BeginIf>chickens. (If-Clause)<Else>roosters. (Else-Clause)<EndIf><UntilEndIfs>
                                            <BeginIf>
                                            Protected from the clear...
                                            <EndIf>
                                            ''')), dedent ( \
                                            '''

                                            Nobody but us 

                                            chickens. (If-Clause)


                                            Foxes like chickens. (If-Clause)
                                            <BeginIf>
                                            Protected from the clear...
                                            <EndIf>
                                            '''))


if __name__ == '__main__' :
    unittest.main()
