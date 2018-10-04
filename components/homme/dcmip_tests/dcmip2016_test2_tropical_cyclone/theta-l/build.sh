#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make theta-l-nlev30"
  make -j theta-l-nlev30
cd $cwd
