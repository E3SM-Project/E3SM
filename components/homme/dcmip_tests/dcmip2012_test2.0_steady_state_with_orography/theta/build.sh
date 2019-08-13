#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j4 theta-nlev30"
  make -j4 theta-nlev30
cd $cwd
