#!/bin/bash
cwd=`pwd`
cd ../..
  echo "make preqx-nlev60-interp"
  make -j4 preqx-nlev60-interp
cd $cwd