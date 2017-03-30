#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make pese-nlev30"
  make -j4 pese-nlev30
cd $cwd