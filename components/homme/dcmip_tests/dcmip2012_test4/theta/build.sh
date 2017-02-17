#!/bin/bash
cwd=`pwd`
cd ../../..
#  echo "make preqx-nlev30-interp"
#  make -j4 preqx-nlev30-interp
  echo "make theta-nlev30, interp"
  make -j4 theta-nlev30
cd $cwd
