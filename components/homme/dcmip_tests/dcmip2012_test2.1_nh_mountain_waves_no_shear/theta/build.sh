#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j theta-nlev60"
  make -j theta-nlev60
cd $cwd
