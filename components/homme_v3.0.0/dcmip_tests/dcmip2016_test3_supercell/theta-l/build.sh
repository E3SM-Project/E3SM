#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j theta-l-nlev40"
  make -j theta-l-nlev40
cd $cwd
