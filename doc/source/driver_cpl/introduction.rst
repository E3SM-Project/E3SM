Introduction
============

The following provides an overview of the CIME driver/coupler. 
We will cover the top level driver implementation as well as the coupler component within the system. 
The driver runs on all hardware processors, runs the top level instructions, and, executes the driver time loop.
The coupler is a component of the CIME infrastructure that is run from within the driver.
It can be run on a subset of the total processors, and carries out mapping (interpolation), merging, diagnostics, and other calculations. 
The name cpl7 refers to the source code associated with both the driver and the coupler parts of the model. 
cpl7 code is located in the CIME source tree under driver_cpl/ and the main program of ``driver_cpl/driver/cesm_driver.F90``.

We also provide a general overview of the cpl7 design. 
Specific implementation issues are then discussed individually. 
Finally, there is a section summarizing all of the cpl7 namelist input. 
This document is written primarily to help users understand the inputs and controls within the cpl7 system, but to also provide some background about the associated implementation. 
`Coupler flow diagrams <http://www.cesm.ucar.edu/models/cesm2.0/cpl7/coupler_flow.pdf>`_ are provided in a separate document. 
Some additional documentation on how the coupler works can be found in Craig et al, `"A New Flexible Coupler for Earth System Modeling Developed for CCSM4 and CESM1" <http://hpc.sagepub.com/content/26/1/31>`_, International Journal of High Performance Computing Applications 2012 26: 31 DOI: 10.1177/1094342011428141.

