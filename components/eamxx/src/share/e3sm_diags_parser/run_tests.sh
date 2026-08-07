#! /bin/bash

cmake --build build
./build/edp_tests -s
