===================
Isostasy test cases
===================

The test cases in this directory demonstrate the range of options for modeling
isostatic bedrock adjustment coupled to ice thickness evolution. The ice
evolution is based on a modified EISMINT-2 test case for an idealized,
axisymmetric ice sheet, which grows and then shrinks as a function of the mass
balance over the course of the model run.

Four configuration files allow for the four possible combinations of modeling
the bedrock response:

* llfa - local lithosphere, fluid asthenosphere
* llra - local lithosphere, relaxed asthenosphere
* elfa - elastic lithosphere, fluid asthenosphere
* elra - elastic lithosphere, relaxed asthenosphere

See the documentation for more discussion of these different models. 

To run a particular test case, first create a symbolic link (or copy of) your
executable, e.g.:

```sh
ln -s ../../../builds/mac-gnu/cism_driver/cism_driver
```

For a serial build, the EISMINT-1, moving margin, test case 3 can be run using:

```sh
./cism_driver isos.elra.config
```

and likewise for other test cases.

For a parallel build, one would use:

```sh
mpirun -np 1 ./cism_driver isos.elra.config 
```

To view the results use ncview or another utility for viewing netCDF files.

(last edited on 2-6-14 by SFP)
