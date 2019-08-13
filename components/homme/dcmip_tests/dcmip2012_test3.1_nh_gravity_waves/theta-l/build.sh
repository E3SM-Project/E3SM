#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j 4 theta-l-nlev20"
  make -j 4 theta-l-nlev20
cd $cwd
