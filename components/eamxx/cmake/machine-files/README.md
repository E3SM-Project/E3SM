# Guide to making good scream/eamxx machine files.

### Basics

A machine file is a cmake macro preload file that makes life
better for developers because, instead of having to set loads of
values manually on the cmake command line, they can just load a
machine file via -C while will set up things in a way that are
known to work for that machine.

These machine files are used for full case (E3SM+CIME) builds and
can also be useful (but are optional) for eamxx standalone builds.
The eamxx standalone scripts infrastructure expects that a machine
file is available.

### Making a machine file

Machine files should have this at the top:
```
include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()
```
This will set SCREAM_MACHINE and look for an EKAT machine file
for this machine. An EKAT machine file is not necessary, but
you'll need to set up EKAT mpi info and kokkos stuff if EKAT
knows nothing about your machine.

You can enable kokkos execution spaces via:
`include (${EKAT_MACH_FILES_PATH}/kokkos/${exe_space}.cmake)`

You can enable kokkos architures via:
`include (${EKAT_MACH_FILES_PATH}/kokkos/${arch}.cmake)`

You can enable certain mpi implementations via:
`include (${EKAT_MACH_FILES_PATH}/kokkos/${mpi_impl}.cmake)`

If EKAT doesn't have presets for you architecture or MPI, you
can always set Kokkos and MPI stuff by hand.

NOTE: If CIME is not aware of this machine, `SCREAM_INPUT_ROOT`
must be set here.

### Making a machine+compiler file

You can have compiler-specific machine files if multiple compilers
are likely to be used on a machine. You'll want to name this file
`${mach}-{compiler}.cmake` and the first line should include the
generic machine file:

```
include(${CMAKE_CURRENT_LIST_DIR}/${mach}.cmake)
```

SCREAM_MACHINE will be set to `${mach}` not `${mach}-${compiler}`.

### Integration with E3SM/CIME

When eamxx is being build as a component within an E3SM case,
these machine files will still be used. It will prefer a
${mach}-${compiler}.cmake over a plain ${mach}.cmake file.

### Gotchas

Due to multiple layers of machine files and even build systems (in
the case of CIME), it is difficult to know what CACHE variables
have already been set and it's easy to find yourself in a situation
where CMake is ignoring your sets because the CACHE var has already
been set. To avoid this, I recommend adding FORCE to your cachevar
sets. This does remove the ability of the user to provide -DVAR=XXX
on the cmake command line, but the whole point of machine files is
to avoid having to set these things on the command line.
