#!/bin/bash
cwd=`pwd`
cd ../..
  echo "make -j 4 preqx-nlev30"
  make -j 4 preqx-nlev30
cd $cwd