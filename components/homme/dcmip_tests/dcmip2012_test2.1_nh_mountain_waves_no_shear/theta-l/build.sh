#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j theta-l-nlev60"
  make -j theta-l-nlev60
cd $cwd
