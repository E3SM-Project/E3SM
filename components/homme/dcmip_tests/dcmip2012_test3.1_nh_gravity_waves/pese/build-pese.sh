#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j 4 pese-nlev20"
  make -j 4 pese-nlev20
cd $cwd