#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make -j 4 preqx-nlev20-interp"
  make -j 4 preqx-nlev20-interp
  make -j 4 theta-nlev20
cd $cwd