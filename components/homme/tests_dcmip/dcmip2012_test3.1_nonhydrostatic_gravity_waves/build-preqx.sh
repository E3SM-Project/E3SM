#!/bin/bash
cwd=`pwd`
cd ../..
  echo "make -j 4 preqx-nlev20"
  make -j 4 preqx-nlev20
cd $cwd