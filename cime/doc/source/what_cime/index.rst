.. _what-cime:

.. on documentation master file, created by
   sphinx-quickstart on Tue Jan 31 19:46:36 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#####################################
 What is CIME?  
#####################################

.. toctree::
   :maxdepth: 3
   :numbered:
      

CIME, pronounced "SEAM", contains the support scripts (configure,
build, run, test), data models, essential utility libraries, a “main”
and other tools that are needed to build a single-executable coupled
Earth System Model.  CIME is available in a stand-alone package that
can be compiled and tested without active prognostic components but is
typically included in the source of a climate model. CIME does not
contain: any active components, any intra-component coupling
capability (such as atmosphere physics-dynamics coupling).

*********
Overview
*********

CIME is comprised of:

1. A Case Control System to support configuration, compilation, execution, system testing and unit testing of a earth system model:

    i. Scripts to enable simple generation of model executables and associated input files for different scientific cases, component resolutions and combinations of full, data and stub components with a handful of commands.
    ii. Testing utilities to run defined system tests and report results for different configurations of the coupled system.

2. A default coupled model architecture:

    i. A programmer interface and libraries to implement a hub-and-spoke inter-component coupling architecture.
    ii. An implementation of a "hub" that needs 7 components (atm, ocn, lnd, sea-ice, land-ice, river, wave). a.k.a. “the driver”.
    iii. The ability to allow active and data components to be mixed in any combination as long as each component implements the coupling programmer interface.

3. Non-active Data and Stub components:

    i. “Data-only” versions of 6 of the 7 components that can replace active components at build-time.
    ii. “Stub” versions of all 7 components for building a complete system.

4. Source code for external libraries useful in scientific applications in general and climate models in particular.
    i.  Parallel I/O library.
    ii. The Model Coupling Toolkit.
    iii. Timing library.

5. Additional stand-alone tools:

    i. Parallel regridding weight generation program
    ii. Scripts to automate off-line load-balancing.
    iii. Scripts to conduct ensemble-based statistical consistency tests.
    iv. Netcdf file comparison program (for bit-for-bit).

*************************
Development
*************************

CIME is developed in an open-source, public repository hosted under the Earth
System Model Computational Infrastructure (ESMCI) organization on
Github at http://github.com/ESMCI/cime.



