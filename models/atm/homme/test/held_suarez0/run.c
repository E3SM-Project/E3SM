#!/bin/sh -x
dir=`pwd`
cqsub -n 486 -t 06:00:00  -O ${dir}/run -C $dir $dir/../../build.Linux/preqx
