====================
EISMINT-2 test cases
====================

To run a particular test case, first create a symbolic link (or copy of) your
executable, e.g.:

```sh
ln -s ../../../builds/mac-gnu/cism_driver/cism_driver
```

For a serial build, EISMINT-2 test A can be run using:

```sh
./cism_driver e2.a.config
```

For a parallel build, one would use:

```sh
mpirun -np 1 ./cism_driver e2.a.config
```

and likewise for other test cases.

Note that some tests use output "hotstart" files (marked with a "hot"
exentsion) from previous tests as input files, so the tests should be run
sequentially, starting with test A.

To view the results use ncview or another utility for viewing netCDF files.

More information on the EISMINT-2 test cases can be found here:
* Payne, A. J. and Coauthors, 2000: Results from the EISMINT model
  intercomparison: the effects of thermomechanical coupling. J Glaciol, 46,
  227–238. <http://homepages.vub.ac.be/~phuybrec/eismint.html>

(last edited on 2-6-14 by SFP)

