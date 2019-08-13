#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j4 theta-l-nlev30"
  make -j4 theta-l-nlev30
cd $cwd
