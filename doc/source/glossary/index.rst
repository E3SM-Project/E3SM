.. _glossary:

########
Glossary
########

.. toctree::
   :maxdepth: 1
   :numbered:

*********
General
*********

.. glossary::

   active or prognostic component  
      Solves a complex set of equations to describe a sub-model’s behavior.

   case (CASE)
      An instance of a global climate model simulation. A case is defined by a component set, a model grid,
      a machine, a compiler, and any other additional customizations.

   component 
      A sub-model coupled with other components to constitute a global climate modeling system.
      Example components: atmosphere, ocean, land, etc.

   component set (COMPSET)
      A complete set of components to be linked together into a climate model to
      run a specific case.

   data component  
      Replacement for an active component. Sends and receives the same variables to and from other
      models (but ignores the variables received).

   grid (GRID)
      A set of numerical grids of a case. Each active component operates on its own numerical grid.

   resolution 
      Used to refer to a set of grids. Each grid within a set may have different resolution.

   stub component 
      Simply occupies the required place in the climate execution sequence and does send or receive
      any data.

*********
Coupling
*********

.. glossary::

   coupler 
      A component of the CIME infrastructure that is run from within the driver. It can be run on a
      subset of the total processors, and carries out mapping (interpolation), merging, diagnostics, and other
      calculations.

   driver 
      The hub that connects all components. CIME driver runs on all hardware processors, runs the top
      level instructions, and, executes the driver time loop.

   forcing
      An imposed perturbation of Earth's energy balance

   Model Coupling Toolkit or MCT 
      A library used by CIME for all data rearranging and mapping (interpolation)

   mask 
      Determines land/ocean boundaries in the model

   mapping 
      Interpolation of fields between components.

*********************
Files and Directories
*********************

.. glossary::

   archive directory (DOUT_S_ROOT)
      If short term archiving is activated (DOUT_S = TRUE), the restart files and run output files
      are copied to archive directory location (DOUT_S_ROOT).

   build directory (EXEROOT)
      Location where the case is built.

   case root (CASEROOT)
      The directory where the case is created. Includes namelist files, xml files, and scripts to setup,
      build, and run the case. Also, includes logs and timing output files. 

   CIME root (CIMEROOT)
      The directory where the CIME source code resides

   history files 
      NetCDF files that contain fields associated with the state of the model at a given time slice. 

   initial files
      Files required to start a file

   input data stream (DIN_LOC_ROOT)
      A time-series of input data files where all the fields in the stream are located in the
      same data file and all share the same spatial and temporal coordinates.

   namelist files 
      Each namelist file includes input parameters for a specific component.

   run directory (RUNDIR)
      Where the case is run.

   restart files 
      Written and read by each component in the RUNDIR to stop and subsequently restart in a bit-for-bit fashion.

   rpointer files 
      Text file written by the coupler in the RUNDIR with a list of necessary files required for model restart.

   XML files 
      Elements and attributes in these files configure a case. (building, running, batch, etc.) These files
      include env_archive.xml, env_batch.xml, env_build.xml, env_case.xml, env_mach_pes.xml, env_mach_specific.xml, env_run.xml
      in CASEROOT and can be queried and modifed using the xmlquery and xmlchange tools. 

***********
Development
***********

.. glossary::

   sandbox (SRCROOT)
      A checked out tag on a local or a remote machine. may be edited to create a new tag. or, may
      just be used for running cases.

   source modifications (CASEROOT/SourceMods)
      one or more source files that are modified by the user. Before building a case, CIME replaces
      the original source files with these files.

   tag 
      A snapshot of the source code. With each consecutive tag (one or more) answer-changing modifications
      to the source code of a component are introduced.

   user namelist files (CASEROOT/user_nl_*)
      User modifications for a given case can be specified in these files.

********
Testing
********

.. glossary::

   baseline 
      A set of test cases that is run using a tag which is complete, tested, and has no modifications
      in the source code. Used to assess the performance/accuracy of a case that is run using a sandbox.

   baseline failure
      A test that fails in its comparison with a baseline.

   blessing
      Part of the unit testing framework used by CIME scripts regression tests. 

   regression test
      A test that compares with baseline results to determine if any new errors have been introduced
      into the code base.

   unit testing 
      A fast, self-verifying test of a small piece of code.

*************
Miscellaneous
*************

.. glossary::

   ESP 
      External System Processing: handles data assimilation