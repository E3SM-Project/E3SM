====================
EISMINT-1 test cases
====================

To run a particular test case, first create a symbolic link (or copy of) your executable, e.g.:

```sh
ln -s ../../../builds/mac-gnu/cism_driver/cism_driver
```

For a serial build, the EISMINT-1, moving margin, test case 3 can be run using:

```sh
./cism_driver e1-mm.3.config
```

and likewise for other test cases ("mm" = moving margin; "fm" = fixed margin).

For a parallel build, one would use:

```sh
mpirun -np 1 ./cism_driver e1-mm.3.config
```

To view the results use ncview or another utility for viewing netCDF files.

More information on the EISMINT-1 test cases can be found here:
* Huybrechts, P., A.J. Payne, 1996: The EISMINT benchmarks for testing
  ice-sheet models. Annals of Glaciology, 23, 1–12.
  <http://homepages.vub.ac.be/~phuybrec/eismint.html>

(last edited on 2-6-14 by SFP)
