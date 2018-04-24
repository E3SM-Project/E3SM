#!/bin/bash

./diff_timeline ne4  /proj/imn/homme/preqx.cpu/jw-ne4-nlev72-qsize40-thread1-nodes1-tasks4/movies/asp_baroclinic2.nc \
                     /proj/imn/homme/preqx.openacc/jw-ne4-nlev72-qsize40-thread1-nodes1-tasks4/movies/asp_baroclinic2.nc diff

cat diff_*.dat
