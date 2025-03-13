# RDycore User Guide - Overview

RDycore's primary purpose is to provide [E3SM](https://e3sm.org/) with the
capability to model coastal compound flooding. Accordingly, it has been
constructed as a performance-portable library that makes efficient use of
DOE's leadership-class computing facilities and is invoked by E3SM Fortran code.

Aside from the library, RDycore provides standalone drivers that can help you
test and evaluate the model's capabilities:

* standalone C and Fortran driver programs for running uncoupled flood
  simulations given appropriate initial and boundary conditions, source terms,
  etc.

* standalone C and Fortran verification programs that use the method of
  manufactured solutions (MMS) to evaluate the stability and accuracy of the
  underlying numerical methods by computing error norms of simulation results
  measured against analytical solutions. These MMS programs can also compute
  convergence rates, which are useful for identifying algorithmic and
  programming errors. Because these programs are more technical and used to
  identify defects, they are described in [the Developer Guide](../developer/index.md).

This guide also describes these standalone programs and their features. It also
explains how we integrate RDycore to E3SM to perform coupled simulations of
compound flooding.

## The Standalone C and Fortran Drivers

* [Standalone YAML input specification](../common/input.md)

