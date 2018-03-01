==================
README.pyEnsSumPop
==================

This package is used to create a summary file for a collection 
(or ensemble) of runs. The user must specify the location of the input files.
The summary file calculates averages, standard deviations, global means, and Z-scores.
This summary file is required for pyCECT.py.

:AUTHORS: Haiying Xu, Allison Baker
:VERSION: 1.0.0
:COPYRIGHT: See the document entitled LICENSE.txt

Send questions and comments to Haiying Xu (haiyingx@ucar.edu).


This package includes:  
----------------------
     	pyEnsSumPop.py             
                            A script that generates an ensemble summary file 
     		            from a collection of files.

        pyEnsLib.py     
                            Library python script used by pyEnsSum.py.

        pyEnsSumPop_test.sh        
                            Bsub script to submit pyEnsSum.py to yellowstone.

        pop_ensemble.json
                            The variable list that will be included for
                            reading and processing.


Before you start to use the package, you need to load the following modules: 
----------------------------------------------------------------------------
       - module load python 
       - module load numpy
       - module load scipy
       - module load pynio
       - Install asaptools
       git clone https://github.com/NCAR-CISL-ASAP/ASAPPyTools
       or
       pip install [--user] ASAPTools
      (Follow the instructions in README.rst in ASAPTools and please make sure
             to set the correct python verion in the PYTHONPATH)

       
To see all options (and defaults):
----------------------------------
       python pyEnsSumPop.py -h

Notes:
------
       For monthly average files, set tslice=0.

       For yearly average files, set tslice=1 (Because tslice ==0 is the initial conditions.)

       Note that --res, --tag, --compset, and --mach only affect the metadata 
       in the summary file.

       Recommended number of cores to use is one for each month. 

Examples for generating summary files:
--------------------------------------
	 (A) To generate (in parallel) a summary file for a set of simulation runs, 
       	 
           we specify the data location:
	    --indir /glade/scratch/haiyingx/pop_ensemble_data60/

           We also specify the name of file to create for the summary:
 	    --sumfile pop.ens.sum.nc

	   Since these are monthly average files, we set
	    --tslice 0

           We also specify the number of years, the number of months, 
           and the number of pertubations (e.g. ensemble size):
            --nyear 1
            --nmonth 12
            --npert 80

	   We also can specify the tag, resolution, machine and compset
	   information (that will be written to the
	   metadata of the summary file):
	    --tag cesm1_4_beta11
            --res T62_g16
            --mach yellowstone
            --compset G

           We can include some variables for the analysis by specifying them 
	   in a json file:
            --jsonfile pop_ensemble.json
       
           To generate only average, standard deviation and zscores (i.e., exclude 
	   global means calculations.  This speeds up the calculation and 
	   is useful for large ensemble sizes if global mean info is not needed.):
            --zscoreonly

           To enable parallel mode:
            --mpi_enable    

           On Yellowstone, please note to set span[ptile=1] in the job submission 
	   script if npert>=60.

	   This yields the following command:

           mpirun.lsf python  pyEnsSumPop.py --verbose --tslice 0 --indir /glade/scratch/haiyingx/pop_ensemble_data60/ --sumfile pop.ens.sum.nc --nyear 1 --nmonth 12 --npert 80 --jsonfile pop_ensemble.json --mpi_enable 




