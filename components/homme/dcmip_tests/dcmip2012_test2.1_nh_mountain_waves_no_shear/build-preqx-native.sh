#!/bin/bash
cwd=`pwd`
cd ../..
  echo "make preqx-nlev60-native"
  make -j preqx-nlev60-native
cd $cwd