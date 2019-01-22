.. _ensemble-consistency-test:

==============================
Ensemble Consistency Test
==============================

CESM-ECT (CESM Ensemble Consistency Test) determines whether a new simulation set up (new machine, compiler, etc.) is statistically distinguishable from an accepted ensemble.
The verification tools in the CESM-ECT suite are:

* CAM-ECT - detects issues in CAM and CTSM (12 month runs)
* UF-CAM-ECT - detects issues in CAM and CTSM (9 time step runs)
* POP-ECT - detects issues in POP and CICE (12 month runs)

The ECT process involves comparing several runs (3) generated with the new scenario to an ensemble built on a trusted machine.
The python ECT tools are located in the pyCECT subdirectory or https://github.com/NCAR/PyCECT/releases.

Before the test can be run, a summary file is needed of the ensemble runs to which the comparison will be made.
Ensemble summary files (NetCDF) for existing tags for CAM-ECT, UF-CAM-ECT and POP-ECT that were created by CSEG are located (respectively) in the CESM input data
directories:

$CESMDATAROOT/inputdata/validation/ensembles
$CESMDATAROOT/inputdata/validation/uf_ensembles
$CESMDATAROOT/inputdata/validation/pop_ensembles

.. todo:: Add more content for ensemble consistency test
