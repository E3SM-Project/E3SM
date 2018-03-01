#!/usr/bin/env python

import unittest

def interleavedp(begins,ends,m=None) :
    
    if( len(begins) != len(ends) ) :
        print 'begin-end token number mismatch'
        # Should learn to throw...
        return False

    if m :
        if len(m) > len(begins) :
            print 'excess else tokens'
            return False

    ok = True
    for i in range(len(begins)-1) :
        ok = ok and begins[i] < ends[i] < begins[i+1] < ends[i+1]
        if not ok :
            print 'begin-end token order mismatch'
            return False

    if m :
        ok = True
        for i in range(len(m)-1) :
            ok = ok and m[i] < m[i+1]
            if not ok :
                print 'else tokens out of order'
        ok = True
        i = 0 ; j = 0 ; notDone = True
        if m[i] < begins[j] :
            print 'else token before if token'
            return False
        while j < len(begins) and notDone :
            k = 0
            while i < len(m) and m[i] < ends[j] :
                i = i + 1
                k = k + 1
                if k > 1 :
                    print 'too many else tokens'
                    return False
            notDone = i < len(m)
            j = j + 1
        if i < len(m):
            print 'else token after endifs token'
            return False
            
    return True


class TestInterleaved(unittest.TestCase):

    def test_InOrder(self):
        a = [1,3,5,7,9]
        b = [2,4,6,8,10]
        self.assertEqual(interleavedp(a,b),True)

    def test_NumberMismatch(self):
        a = [1,3,7,9]
        b = [2,4,6,8,10]
        self.assertEqual(interleavedp(a,b),False)
        
    def test_OrderMismatch1(self):
        a = [1,3,5,7,9]
        b = [2,4,6,10,12]
        self.assertEqual(interleavedp(a,b),False)

    def test_OrderMismatch2(self):
        a = [1,3,5,7,9]
        b = [4,6,10,12,14]
        self.assertEqual(interleavedp(a,b),False)

    def test_OrderMismatch3(self):
        a = [1,3,5,7,11]
        b = [2,4,6,8,10]
        self.assertEqual(interleavedp(a,b),False)

    def test_ElseMid1(self):
        a = [10]
        m = [15]
        b = [20]
        self.assertEqual(interleavedp(a,b,m=m),True)

    def test_ElseMid2(self):
        a = [10]
        m = [15]
        b = [15]
        self.assertEqual(interleavedp(a,b,m=m),False)

    def test_ElseMid3(self):
        a = [10]
        m = [15]
        b = [15]
        self.assertEqual(interleavedp(a,b,m=m),False)

    def test_ElseMid4(self):
        a = [10]
        m = [5]
        b = [20]
        self.assertEqual(interleavedp(a,b,m=m),False)

    def test_ElseMid5(self):
        a = [10,30]
        m = [15]
        b = [20,40]
        self.assertEqual(interleavedp(a,b,m=m),True)

    def test_ElseMid6(self):
        a = [10,30]
        m = [15,17]
        b = [20,40]
        self.assertEqual(interleavedp(a,b,m=m),False)

    def test_ElseMid7(self):
        a = [10,30]
        m = [25,30]
        b = [20,40]
        self.assertEqual(interleavedp(a,b,m=m),False)

    def test_ElseMid8(self):
        a = [10,30]
        m = [25,31,35]
        b = [20,40]
        self.assertEqual(interleavedp(a,b,m=m),False)

    def test_ElseMid9(self):
        a = [10,30,50]
        m = [25,55,35]
        b = [20,40,60]
        self.assertEqual(interleavedp(a,b,m=m),False)

    def test_ElseMid10(self):
        a = [10,30,50]
        m = [15,35,55]
        b = [20,40,60]
        self.assertEqual(interleavedp(a,b,m=m),True)

if __name__ == '__main__' :
    unittest.main()
