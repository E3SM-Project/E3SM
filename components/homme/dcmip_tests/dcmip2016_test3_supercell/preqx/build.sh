#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make preqx-nlev40-interp"
  make -j preqx-nlev40-interp
cd $cwd
