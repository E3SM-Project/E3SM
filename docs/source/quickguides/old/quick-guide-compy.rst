
Quick guide for running in e3sm_unified environment for Compy and Others (v1)
=============================================================================

Activate the latest version of the e3sm_unified environment
-----------------------------------------------------------

Most of the E3SM analysis software is maintained with an Anaconda metapackage. To get all of the tools in the metapackage in your path, use one of the sets of commands below and the activation path from below for individual e3sm analysis machines. Shown below we also provide the paths where climatology observational data e3sm_diags uses, replace ``/climatology`` to ``/time-series`` for time-series data:


**Compy**
    ::

     source /compyfs/software/e3sm-unified/load_latest_e3sm_unified_py2.7.sh


obs at: ``/compyfs/e3sm_diags_data/obs_for_e3sm_diags/climatology/``

     

**For other analysis platforms**

**Cori**
    ::

     source /global/project/projectdirs/acme/software/anaconda_envs/load_latest_e3sm_unified_py2.7.sh
    
obs at: ``/global/project/projectdirs/acme/e3sm_diags/obs_for_e3sm_diags/climatology/``


**Anvil/blues**
    ::

     source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_py2.7.sh

obs at: ``/lcrc/soft/climate/e3sm_diags_data/obs_for_e3sm_diags/climatology/``


**Cooley**
    ::

     source /lus/theta-fs0/projects/ccsm/acme/tools/e3sm-unified/load_latest_e3sm_unified_py2.7.sh

obs at:``/lus/theta-fs0/projects/ClimateEnergy_3/e3sm_diags/obs_for_e3sm_diags/climatology/``


**acme1**
    ::

     source /usr/local/e3sm_unified/envs/load_latest_e3sm_unified_py2.7.sh

obs at:``/p/cscratch/acme/data/obs_for_e3sm_diags/climatology/``


**Rhea**
    ::

     source /ccs/proj/cli900/sw/rhea/e3sm-unified/load_latest_e3sm_unified_py2.7.sh
 
obs at:``/ccs/proj/cli115/e3sm_diags_data/obs_for_acme_diags/climatology/``


Change ``.sh`` to ``.csh`` for csh shells.


Running the entire annual latitude-longitude contour set
--------------------------------------------------------

Copy and paste the below code into ``myparams.py`` using your favorite text editor. Adjust any options as you like.

   **Tip:** Make a folder in the following directory ``/compyfs/www/`` based off your username.
   Then you can set ``results_dir`` to  ``/compyfs/www/<username>/lat_lon_demo`` in ``myparams.py`` below
   to view the results via a web browser here: https://compy-dtn.pnl.gov/<username>/lat_lon_demo


    .. code:: python

        reference_data_path = '/compyfs/e3sm_diags_data/obs_for_e3sm_diags/climatology/'
        test_data_path = '/compyfs/e3sm_diags_data/test_model_data_for_acme_diags/climatology/'

        test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

        sets = ["lat_lon"]
        seasons = ["ANN"]

        # 'mpl' and 'vcs' are for matplotlib or vcs plots respectively.
        backend = 'mpl'

        # Name of folder where all results will be stored.
        results_dir = 'lat_lon_demo'

        multiprocessing = True
        num_workers =  40

To enable multiprocessing rather than running in serial, the program will need to be ran in an
**interactive session** on compute nodes, or as a **batch job**.

Here are some hardware details for `compy`
   * 40 cores/node
   * 192 GB DRAM/node
   * 18400 total cores


Interactive session on compute nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, request an interactive session with a single node 
for one hour (running this example should take much less than this).

    ::

        salloc --nodes=1 --account=e3sm --time=00:30:00 


Once the session is available, launch E3SM Diagnostics.

    ::

        source /compyfs/software/e3sm-unified/load_latest_e3sm_unified_py2.7.sh
        e3sm_diags -p myparams.py


Batch job
^^^^^^^^^

Alternatively, you can also create a script and submit it to the batch system.
Copy and paste the code below into a file named ``diags.bash``.
Please remember to change what directory you're in to one accessible to you.

    .. code:: bash
    
        #!/bin/bash -l
        #SBATCH --job-name=diags
        #SBATCH --output=diags.o%j
        #SBATCH --account=e3sm
        #SBATCH --nodes=1
        #SBATCH --time=00:30:00

        # Please change the directory below.
        source /compyfs/software/e3sm-unified/load_latest_e3sm_unified_py2.7.sh
        e3sm_diags -p myparams.py

And then submit it

    ::

        sbatch diags.bash

View the status of your job with ``squeue -u <username>``.
Here's the meaning of some values under the State (``ST``) column:

* ``PD``: Pending
* ``R``: Running
* ``CA``: Cancelled
* ``CD``: Completed
* ``F``: Failed
* ``TO``: Timeout
* ``NF``: Node Failure


Back to running the latitude-longitude contour set
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
5. Once you ran the diagnostics in an interactive session or via a batch job, open the following webpage to view the results.


    ::

        lat_lon_demo/viewer/index.html

**Tip:** Once you're on the webpage for a specific plot, click on the
'Output Metadata' drop down menu to view the metadata for the displayed plot.
Running that command allows the displayed plot to be recreated.
Changing any of the options will modify the just that resulting figure.



Running all of the diagnostics sets
-----------------------------------

Copy and paste the following into ``all_sets.py`` using your
favorite text editor:

    .. code:: python

        reference_data_path = '/compyfs/e3sm_diags_data/obs_for_e3sm_diags/climatology/'
        test_data_path = '/compyfs/e3sm_diags_data/test_model_data_for_acme_diags/climatology/'

        test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

        # Not defining a sets parameter runs all of the default sets:
        # ['zonal_mean_xy', 'zonal_mean_2d', 'lat_lon', 'polar', 'cosp_histogram']
        sets = ['zonal_mean_xy', 'zonal_mean_2d', 'lat_lon', 'polar', 'cosp_histogram']

        # 'mpl' and 'vcs' are for matplotlib or vcs plots respectively.
        backend = 'mpl'

        # Name of folder where all results will be stored.
        results_dir = 'diag_demo'

        # Optional settings below:

        diff_title = 'Model - Obs'

        multiprocessing = True
        num_workers =  40


Advanced: Running custom diagnostics
------------------------------------
The following steps are for 'advanced' users, who want to run custom diagnostics.
So most users will not run the software like this.


By default, all of the E3SM diagnostics are ran for the sets that we defined above.
This takes some time, so instead we create our own diagnostics to be ran.


Copy and paste the code below in ``mydiags.cfg``.
Check :doc:`Available Parameters <../../available-parameters>`
for all available parameters.

For more examples of these types of files, look
`here <https://github.com/E3SM-Project/e3sm_diags/blob/master/e3sm_diags/driver/default_diags/lat_lon_model_vs_obs.cfg>`_
for the cfg file that was used to create all of the latitude-longitude sets.


    ::

        [#]
        sets = ["lat_lon"]
        case_id = "GPCP_v2.2"
        variables = ["PRECT"]
        ref_name = "GPCP_v2.2"
        reference_name = "GPCP (yrs1979-2014)"
        seasons = ["ANN", "DJF"]
        regions = ["global"]
        test_colormap = "WhiteBlueGreenYellowRed.rgb"
        reference_colormap = "WhiteBlueGreenYellowRed.rgb"
        diff_colormap = "BrBG"
        contour_levels = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16]
        diff_levels = [-5, -4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4, 5]

        [#]
        sets = ["lat_lon"]
        case_id = "SST_CL_HadISST"
        variables = ["SST"]
        ref_name = "HadISST_CL"
        reference_name = "HadISST/OI.v2 (Climatology) 1982-2001"
        seasons = ["ANN", "MAM"]
        contour_levels = [-1, 0, 1, 3, 6, 9, 12, 15, 18, 20, 22, 24, 26, 28, 29]
        diff_levels = [-5, -4, -3, -2, -1, -0.5, -0.2, 0.2, 0.5, 1, 2, 3, 4, 5]

Run E3SM diagnostics with the ``-d`` parameter.

    ::

        e3sm_diags -p myparams.py -d mydiags.cfg


