#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make theta-l-nlev30, interp"
  make -j4 theta-l-nlev30
cd $cwd
