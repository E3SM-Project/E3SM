#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j theta-nlev40"
  make -j theta-nlev40
cd $cwd
