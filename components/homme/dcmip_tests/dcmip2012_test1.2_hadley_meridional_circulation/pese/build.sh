#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make pese-nlev60"
  make -j4 pese-nlev60
cd $cwd