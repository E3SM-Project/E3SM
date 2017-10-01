#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j theta-nlev60p"
  make -j theta-nlev60
cd $cwd
