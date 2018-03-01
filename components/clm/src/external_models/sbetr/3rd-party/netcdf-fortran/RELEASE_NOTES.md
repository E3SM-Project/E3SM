Release Notes {#nf_release_notes}
==============================

\brief Release notes file for the netcdf-fortran package.

This file contains a high-level description of this package's evolution.
Entries are in reverse chronological order (most recent first).

## 4.4.3 Released 2016-01-20

* Corrected a bug which would return a false-positive in `nf_test` when using netCDF-C `4.4.0`.

* Updated the `cfortran.doc` license document for the `cfortran.h` library.  The most recent version was pulled from http://cfortran.sourceforge.net.  The previous version did not reflect that the author had released cfortran under the LGPL.  See [Github Issue 27](https://github.com/Unidata/netcdf-fortran/issues/27) for more information.

## 4.4.2 Released 2015-02-02

* Added infrastructure to support the new `netcdf-c` option, `ENABLE_REMOTE_FORTRAN_BOOTSTRAP`.

* Incorporated changes submitted by Nico Schlomer which extends the cmake compatibility between `netcdf-c` and `netcdf-fortran`.

* Incorporated a patch submitted by Thomas Jahns which fixed `FC` being unconditionally overwritten by `F77` when `Fortran 90` was disabled.

## 4.4.1 Released 2014-09-09

* No significant changes from RC1.

### 4.4.1-RC1 Released 2014-08-05

* Added a new variable for cmake-based builds, `NC_EXTRA_DEPS`.  Use this to specify additional dependencies when linking against a static `netcdf-c` library, e.g.

```.fortran
netcdf-fortran/build$ cmake .. -DNC_EXTRA_DEPS="-lhdf5 -lhdf5_hl -lcurl"
```

* Fixed to build correctly with netCDF-3-only C library, for example C library configured with --disable-netcdf-4 (R. Weed).

## 4.4 Released 2014-07-08

* For 32-bit platforms fixed integer fill parameters, initialized potentially
  unitialized variables, and provided some missing defaults (R. Weed).

* Fixed CMake builds on 32-bit platforms.

* Added new `inq_path` and `rename_grps` functions analogous to
  corresponding C functions. Added associated tests (R. Weed).

* Added support for NF\_MPIIO, `NF_MPIPOSIX`, `NF_PNETCDF` flags and
  `NF_FILL_UINT`. (R. Weed)

* Fixed potential bug in attribute functions for integer values when
  Fortran `INTEGER*1` or `INTEGER*2` types are the same size as C
  long (R. Weed).

* Added test for compiler support of Fortran 2008 `ISO_FORTRAN_ENV`
  additions and TS29113 standard extension.

* Fixed `C_PTR_DIFF_T` issue reported by Orion Poplowski (R. Weed).

### 4.4-rc1 	Released 2013-10-06

* Added doxygen-generated documentation, using the `--enable-doxygen` and `-DENABLE_DOXYGEN` flags for autotools and cmake-based builds, respectively.

* Added missing error codes for DAP and some netCDF-4 errors

* Fixed some documentation for F77 API, added make rule for creating netcdf-f77 HTML files.

### 4.4-beta5 	Released 2013-08-27

* Added configuration files to github distribution.

### 4.4-beta4      

* Moved to GitHub from Subversion, the location of the new GitHub repository is at: http://github.com/Unidata/netCDF-Fortran

* Parallel-build portability fixes, particularly for
		OpenMPI and gcc/gfortran-4.8.x on the Mac.  Also added
		test from Reto Stöckli for NCF-250 bug, demonstrating
		it was fixed in previous commit.

* Add support for NF\_MPIIO, NF\_MPIPOSIX, NF\_PNETCDF, and
		NF\_FILL\_UINT in the data files.

* Add support for nf\_inq\_path.

* Add a pre-processor macro that can be used to bypass
		the home-brew C_PTRDIFF_T definition and use the
		standard one for compilers that support it.

* Fix a potential bug in nf\_attio to call the \_long
		version of some puts/gets instead of the \_int
		version. These were inside INT1\_IS\_C\_LONG and
		INT2\_IS\_C\_LONG ifdef blocks so they would have only
		showed up when those macros were true.

### 4.4-beta3	Released 2012-12-07

* Fixed bug that "make -j check" fails, but "make check" works fine.

* Fixed build problems resulting from syncing with separate C distribution.

* Synchronize with all changes made to version 4.2 since ts release.

### 4.4-beta2	Released 2012-06-29

* Made handling of --disable-f03 more transparent.

* Fixed adding flags for parallel I/O for MPI from David Warren.

* Removed all the old C code that's not needed for this separate distribution.

* Inadvertently broke the build until syncing with C distribution in later beta release.

### 4.4-beta1	Released 2012-03-02

* `Fortran 2003 Support`

    Version 4.4 is the first release to support fortran 2003 and to use the ISO C Bindings available in fortran 2003 to replace the older C code wrappers.

    Congratulations and thanks to Richard Weed at Mississippi State University, who is the author of new code.

    See the file `README_F03_MODS` for a more complete description of the changes. Many changes to the build structure have been made at the same time as the new 2003 code has been inserted.

    As part of the fortran 2003 refactor, the directory structure has been significantly modified.  All the previous F90 C wrapper code has been moved to the "libsrc" directory.

    All of the fortran code has been moved to the "fortran" directory. The directories names F77 and F90 have been removed. The most important consequence of this refactor is that pure Fortran77 compilers are no longer supported. It is assumed that the compiler supports at least Fortran 90 and also Fortran 77.  If it also supports the ISO C Bindings, then the new 2003 code is used instead of the older C wrappers.
