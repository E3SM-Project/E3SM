#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j 4 theta-nlev20"
  make -j 4 theta-nlev20
cd $cwd
