
Quick guide for NERSC Edison (v1)
=================================

Installation
------------

1. Log on to NERSC Edison

2. Make sure you're using ``bash``

::

    bash

3. Load the Anaconda module

::

    module load python/2.7-anaconda-4.4

4. Get the yml file to create an environment.

::

    wget https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/conda/e3sm_diags_env.yml


5. Allow Anaconda to download packages, even with a firewall.

::

    conda config --set ssl_verify false
    binstar config --set verify_ssl False


6. Use Anaconda to create a new environment with ``e3sm_diags`` installed.
Tip: You can change the name of the environment by adding ``-n new_env_name`` to the end of ``conda env create ...``.

::

    conda env create -f e3sm_diags_env.yml
    source activate e3sm_diags_env


Running the entire Latitude-longitude contour set
-------------------------------------------------

7. Copy and paste the below code into ``myparams.py`` using your
favorite text editor. Adjust any options as you like. 

.. code:: python

    reference_data_path = '/global/project/projectdirs/acme/acme_diags/obs_for_e3sm_diags/climatology/'
    test_data_path = '/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/climatology/'

    test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

    sets = ["lat_lon"]

    backend = 'mpl'  # 'vcs' is for vcs plots

    results_dir = 'lat_lon_demo'  # name of folder where all results will be stored


8. Run E3SM diags.

::

    e3sm_diags -p myparams.py

9. Open the following webpage to view the results:

::

    lat_lon_demo/viewer/index.html


Tip: Once you're on the webpage for a specific plot, click on the 'Output Metadata' 
drop down menu to view the metadata for the displayed plot.

* Running that command allows the displayed plot to be recreated. Changing any of the options will modify the resulting figure.

Running all of the diagnostics sets
-----------------------------------

Copy and paste the following into ``all_sets.py`` using your
favorite text editor:

.. code:: python

  reference_data_path = '/global/project/projectdirs/acme/acme_diags/obs_for_acme_diags/'
  test_data_path = '/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/'

  test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

  # Not defining a sets parameter runs all of the default sets:
  # ['zonal_mean_xy', 'zonal_mean_2d', 'lat_lon', 'polar', 'cosp_histogram']

  # optional settings below
  diff_title = 'Model - Obs'

  backend = 'mpl'  # 'mpl' is for the matplotlib plots.

  results_dir = 'diag_demo'  # name of folder where all results will be stored

  multiprocessing = True
  num_workers =  24

Compared to the previous short test above, note the following changes:

* Generate plots for all the available sets ('zonal_mean_xy', 'zonal_mean_2d', 
  'lat_lon', 'polar', 'cosp_histogram').
* Turn on multiprocessing with 24 workers.

Since the example above turns on multiprocessing, it should not be run interactively
on the Edison login nodes (NERSC would likely flag it and let you know about it).
Instead, it can be run either in an interactive session on compute nodes, or as a batch
job.


Interactive session on compute nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, request an interactive session with a single node (24 cores) for one hour
(running this example should take much less than this): ::

  salloc --nodes=1 --partition=regular --time=01:00:00

Once the session is available, launch E3SM Diags: ::

  source activate e3sm_diags_env
  e3sm_diags -p all_sets.py

Batch job
^^^^^^^^^

Alternatively, you can also create a script and submit it to the batch system.
Copy and paste the code below into a file named ``diags.bash``:

.. code:: bash
 
  #!/bin/bash -l
  #SBATCH --job-name=diags
  #SBATCH --output=diags.o%j
  #SBATCH --partition=regular
  #SBATCH --account=acme
  #SBATCH --nodes=1
  #SBATCH --time=01:00:00

  source activate e3sm_diags_env
  cd /global/cscratch1/sd/golaz/tmp
  e3sm_diags -p all_sets.py

And then submit it ::

  sbatch diags.bash

That's it!


Advanced: Running custom diagnostics
------------------------------------
The following steps are for 'advanced' users, who want to run custom diagnostics.
So most users will not run the software like this.

By default, all of the E3SM diagnostics are ran for the ``sets`` that
we defined above. This takes some time, so instead we create our own
diagnostics to be ran.

10. Copy and paste the code below in ``mydiags.cfg``.
Check :doc:`defining parameters <../../available-parameters>`
for all available parameters.

::

    [#]
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
    case_id = "SST_CL_HadISST"
    variables = ["SST"]
    ref_name = "HadISST_CL"
    reference_name = "HadISST/OI.v2 (Climatology) 1982-2001"
    seasons = ["ANN", "MAM"]
    contour_levels = [-1, 0, 1, 3, 6, 9, 12, 15, 18, 20, 22, 24, 26, 28, 29]
    diff_levels = [-5, -4, -3, -2, -1, -0.5, -0.2, 0.2, 0.5, 1, 2, 3, 4, 5]

11. Run E3SM diags.

::

    e3sm_diags -p myparams.py -d mydiags.cfg

12. Open the following webpage to view the results:

::

    lat_lon_demo/viewer/index.html
