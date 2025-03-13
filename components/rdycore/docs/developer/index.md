# RDycore Developer Guide - Overview

This document contains information for anyone interested in helping to
develop and/or maintain RDycore.

RDycore compriseѕ several related components:

* a C library that implements the actual river dynamical core model using
  [PETSc](https://petsc.org/release/), which provides performance portable
  numerical solvers, data structures like unstructured grids and vectors,
  and basic system-level utilities
* a Fortran library consisting of a thin wrapper around the C library, with
  exactly the same functionality
* standalone C and Fortran driver programs that demonstrate the capabilities of
  the RDycore library using RDycore's [YAML input format](../common/input.md)
* various other programs and tools, including C and Fortran drivers that
  can evaluate rates of convergence in the solution error for RDycore when
  applied to specific analytical problems

Most of these components are interesting only to RDycore developers. For those
who wish to incorporate RDycore within a larger framework like [E3SM](https://e3sm.org),
the C and Fortran RDycore libraries (`librdycore.a` and `librdycore_f90.a`, when
built) and their corresponding header/module files (`rdycore.h` and `rdycore.mod`,
respectively) are all that is needed.

### The RDycore library

## Contents

* [Building and Installing RDycore](../common/installation.md)
* [Development Process](development.md)
* [Code Structure and Organization](organization.md)
* [Developer Tools](tools.md)
* [Integrating RDycore with E3SM](e3sm.md)
* [Mesh Description](mesh.md)
