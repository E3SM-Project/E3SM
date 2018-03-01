===============
README.pyEnsSum
===============

This package is used to create a summary file for a collection
(or ensemble) of runs. The user must specify the location of the input files.
The summary file calculates global means, RMSZ scores, PCA loadings, and max errors.
This summary file is required for pyCECT.py.

:AUTHORS: Haiying Xu, Allison Baker
:VERSION: 1.0.0
:COPYRIGHT: See the document entitled LICENSE.txt

Send questions and comments to Haiying Xu (haiyingx@ucar.edu).


This package includes:
----------------------
     	pyEnsSum.py
                            A script that generates an ensemble summary file
     		            from a collection of files.

        pyEnsLib.py
                            Library python script used by pyEnsSum.py

        pyEnsSum_test.sh
                            Bsub script to submit pyEnsSum.py to yellowstone

        ens_excluded_varlist.json
                            The variable list that will excluded from
                            reading and processing


Before you start to use the package, you need to load the following modules:
----------------------------------------------------------------------------
       - module load python
       - module load numpy
       - module load scipy
       - module load pynio
       - svn co https://proxy.subversion.ucar.edu/pubasap/pyTools/tags/v0.3 ASAPTool
       Note: need to install asaptools and setup PYTHONPATH by following
             the instruction at README.rst in ASAPTool, please make sure
             to set the correct python verion in the PYTHONPATH

To see all options (and defaults):
----------------------------------
       python pyEnsSum.py -h

Notes:
------
       For monthly average files, set tslice=0.

       For yearly average files, set tslice=1 (Because tslice ==0 is the initial conditions.)

       Esize can be less than or equal to the number of files in "--indir".

       Note that --res, --tag, --compset, and --mach only affect the metadata
       in the summary file.

       Recommended number of cores to use is one for each 3D variable (current
       default number of 3D variables is 42).

Examples for generating summary files:
--------------------------------------
	 (A) To generate (in parallel) a summary file for 151 simulations runs,

           we specify the size and data location:
	    --esize 151
	    --indir /glade/u/tdd/asap/verification/cesm1_3_beta11/sz151-yellowstone-intel/

           We also specify the name of file to create for the summary:
 	    --sumfile intel_summary.nc

	   Since these are yearly average files, we set
	    --tslice 1

	   We also specify the tag (cesm1_3_beta110 that will be written to the
	   metadata of intel_summary.nc):
	    --tag cesm1_3_beta11

           We can exclude some variables from the analysis by specifying them
	   in a json file:
            --jsonfile ens_excluded_varlist.json

           To generate only global_mean and related PCA loadings (i.e., exclude
	   RMSZ and max-error calculations.  This speeds up the calculation and
	   is useful for large ensemble sizes if RMSZ info is not needed.):
            --gmonly

           To enable parallel mode:
            --mpi_enable

	   This yields the following command:

           mpirun.lsf python  pyEnsSum.py --verbose --esize 151 --tslice 1 --indir /glade/u/tdd/asap/verification/cesm1_3_beta11/sz151-yellowstone-intel/ --tag cesm1_3_beta11 --sumfile intel_test.nc --jsonfile ens_excluded_varlist.json --gmonly --mpi_enable



	 (B) To generate (in serial) a summary file for 151 simulations runs,

           python  pyEnsSum.py --verbose --esize 151 --tslice 1 --indir /glade/u/tdd/asap/verification/cesm1_3_beta11/sz151-yellowstone-intel/ --tag cesm1_3_beta11 --sumfile intel_test.nc --jsonfile ens_excluded_varlist.json

