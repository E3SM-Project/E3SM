#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j theta-nlev30"
  make -j theta-nlev30
cd $cwd
