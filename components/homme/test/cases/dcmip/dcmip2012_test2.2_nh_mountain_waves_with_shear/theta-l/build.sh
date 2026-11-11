#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j theta-l-nlev60p"
  make -j theta-l-nlev60
cd $cwd
