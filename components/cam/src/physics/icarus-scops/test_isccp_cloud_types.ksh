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

clear

make test_isccp_cloud_types || exit 1

# Test RNG gives correct results

./test_congvec.ksh || exit 1

echo Running ISCCP simulator tests - may take a few minutes....

rm -f stdout

./rcsid isccp_cloud_types.f77 > stdout
for top_height_direction in 1 2
do
  for top in 1 2 3
  do
    for overlap in 1 2 3
    do
      echo $top > stdin
      echo $overlap >> stdin
      echo $top_height_direction >> stdin
      (
        echo top=$top
        echo overlap=$overlap
        echo top_height_direction=$top_height_direction
        rm -f ftn09.*
         
        if [ $(hostname) = "tx01" ]
        then
            dir=$(pwd)
            rsh sx601 "cd $dir ; ./test_isccp_cloud_types < stdin"
        else
            ./test_isccp_cloud_types < stdin
        fi

      ) | sed 's/  \./ 0./g;s/ *$//' >> stdout
      for file in ftn09.*
      do
        export LC_ALL=C
        #sort < $file | uniq -c >> stdout
        sort < $file |sort >> stdout
      done
    done
  done
done

if diff stdout stdout.expected > stdout.diff
then
  echo tests passed ok.
  exit 0
else
  echo there may be a problem with the test - files stdout and stdout.expected do not match.
  echo tkdiff stdout stdout.expected
  less stdout.diff
  exit 1
fi

