#!/bin/bash
cwd=`pwd`
cd ../../..
  echo "make preqx-nlev60-interp"
  make -j preqx-nlev60-interp theta-nlev60
cd $cwd