#!/bin/ksh

# *****************************COPYRIGHT****************************
# (c) British Crown Copyright 2009, the Met Office.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the
# following conditions are met:
# 
#     * Redistributions of source code must retain the above 
#       copyright  notice, this list of conditions and the following 
#       disclaimer.
#     * Redistributions in binary form must reproduce the above 
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials 
#       provided with the distribution.
#     * Neither the name of the Met Office nor the names of its 
#       contributors may be used to endorse or promote products
#       derived from this software without specific prior written 
#       permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR 
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
# 
# *****************************COPYRIGHT*******************************
# *****************************COPYRIGHT*******************************

# $Id: test_congvec.ksh,v 4.0 2009/02/13 08:21:07 hadmw Exp $ 

rm -f test_congvec congvec.out congvec.compare
make test_congvec || exit 1

#grep '^F77=' Makefile > congvec.out
# NEC sx6 rsh from front end
if [ $(hostname) = "tx01" ]
then
    dir=$(pwd)
    rsh sx601 "cd $dir ; ./test_congvec " >> congvec.out
else
    ./test_congvec >> congvec.out
fi

#tkdiff congvec.out congvec.expected

cat congvec.out | cut -c1-6 | sort | uniq -c | sort -k2 -n > congvec.compare

if diff -w congvec.compare congvec.expected
then
  echo congvec random number generator tests passed ok.
  exit 0
else
  echo there may be a problem with the test - files congvec.compare and congvec.expected do not match.
  set -x
  grep ^F77= Makefile
  head -20 congvec.out
  tkdiff congvec.compare congvec.expected
  exit 1
fi


