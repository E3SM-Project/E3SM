.. _ensemble-consistency-test:

==============================
CESM-ECT (CESM Ensemble Consistency Test):
==============================

CESM-ECT is a suite of tests to determine whether a new
simulation set up (new machine, compiler, etc.) is statistically
distinguishable from an accepted ensemble.  The verification tools in
the CESM-ECT suite are:

CAM-ECT - detects issues in CAM and CLM (12 month runs)
UF-CAM-ECT - detects issues in CAM and CLM (9 time step runs)
POP-ECT - detects issues in POP and CICE (12 month runs)

The ECT process involves comparing runs generated with
the new scenario ( 3 for CAM-ECT and UF-CAM-ECT, and 1 for POP-ECT) 
to an ensemble built on a trusted machine (currently
cheyenne). The python ECT tools are located in the pyCECT
subdirectory or https://github.com/NCAR/PyCECT/releases.

-OR-
		
We now provide a web server for CAM-ECT and UF-CAM-ECT,  where 
you can upload the (3) generated runs for comparison to our ensemble.
Please see the webpage at http://www.cesm.ucar.edu/models/cesm2/verification/ 
for further instructions.

-----------------------------------
Creating or obtaining a summary file:
-----------------------------------

Before the test can be run, a summary file is needed of the ensemble
runs to which the comparison will be made. Ensemble summary files
(NetCDF) for existing tags for CAM-ECT, UF-CAM-ECT, and POP-ECT that
were created by CSEG are located (respectively) in the CESM input data
directories:

$CESMDATAROOT/inputdata/validation/ensembles
$CESMDATAROOT/inputdata/validation/uf_ensembles
$CESMDATAROOT/inputdata/validation/pop_ensembles

If none of our ensembles are suitable for your needs, then you may create
your own ensemble (and summary file) using the following instructions:
     
(1) To create a new ensemble, use the ensemble.py script in this directory. 
This script creates and compiles a case, then creates clones of the
original case, where the initial temperature perturbation is slightly modified
for each ensemble member.  At this time, cime includes functionality
to create ensembles for CAM-ECT, UF-CAM-ECT, and POP-ECT. 

(2) Use --ect <pop,cam> to specify whether ensemble is for CAM or POP.
(See 'python ensemble.py -h' for additional details).

(3) Use --ensemble <size> to specify the ensemble size.
Recommended ensemble sizes:
CAM-ECT: 151
UF-CAM-ECT: 350
POP-ECT 40

(4) Examples:

CAM-ECT:

python ensemble.py --case /glade/scratch/cesm_user/cesm_tag/ensemble/ensemble.cesm_tag.000 --mach cheyenne   --ensemble 151 --ect cam --project P99999999


UF-CAM-ECT:

python ensemble.py --case /glade/scratch/cesm_user/cesm_tag/uf_ensemble/ensemble.cesm_tag.uf.000 --mach cheyenne  --ensemble 350 --uf --ect cam --project P99999999

POP-ECT:

python ensemble.py --case /glade/scratch/cesm_user/cesm_tag/uf_ensemble/ensemble.cesm_tag.000 --mach cheyenne  --ensemble 40 --ect pop --project P99999999

Notes: 
       (a) ensemble.py accepts (most of) the argumenets of create_newcase

       (b) case name must end in ".000" and include the full path

       (c) ensemble size must be specified, and suggested defaults are listed
       	   above. Note that for CAM-ECT and UF-CAM-ECT, the ensemble size 
	   needs to be larger than the number of variables that ECT will evaluate.


(5) Once all ensemble simulations have run successfully, copy every cam history 
file (*.cam.h0.*) for CAM-ECT and UF-CAM-ECT) or monthly pop history file 
(*.pop.h.*) for POP-ECT from each ensemble run directory into a separate directory. 
Next create the ensemble summary using the pyCECT tool pyEnsSum.py (for CAM-ECT and
UF-CAM-ECT) or pyEnsSumPop.py (for POP-ECT).  For details see README_pyEnsSum.rst 
and README_pyEnsSumPop.rst with the pyCECT tools.

-------------------
Creating test runs:
-------------------

(1) Once an ensemble summary file has been created or chosen to
use from $CESMDATAROOT/inputdata/validation, the simulation
run(s) to be verified by ECT must be created via script ensemble.py.

NOTE: It is important that the **same** resolution and compset be used in the
individual runs as in the ensemble.  The NetCDF ensemble summary file global
attributes give this information.

(2) For example, for CAM-ECT:

python ensemble.py --case /glade/scratch/cesm_user/cesm_tag/camcase.cesm_tag.000 --ect cam --mach cheyenne --project P99999999
--compset   F2000climo --res f19_f19 
For example, for UF-CAM-ECT:

python ensemble.py --case /glade/scratch/cesm_user/cesm_tag/uf.camcase.cesm_tag.000 --ect cam --uf --mach cheyenne --project P99999999 --compset   F2000climo --res f19_f19 

For example, for POP-ECT:

python ensemble.py --case /glade/scratch/cesm_user/cesm_tag/popcase.cesm_tag.000 --ect pop --mach cheyenne  --project P99999999 --compset   G --res T62_g17 

(3) Next verify the new simulation(s) with the pyCECT tool pyCECT.py (see
README_pyCECT.rst with the pyCECT tools).
