#!/bin/bash

cd swtc2
rm *.png
./run_case.sh $1 new
./run_case.sh $1 orig
python3 make_plots.py
cd ..

cd swtc5
rm *.png
./run_case.sh $1 new
./run_case.sh $1 orig
python3 make_plots.py
cd ..

cd swtc6
rm *.png
./run_case.sh $1 new
./run_case.sh $1 orig
python3 make_plots.py
