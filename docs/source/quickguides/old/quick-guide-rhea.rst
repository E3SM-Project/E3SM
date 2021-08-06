Quick guide for OLCF Rhea (v1)
==============================

Running the software on Rhea shares similar steps as running on other machines. The path of datasets are different for different file systems.
To run ``e3sm_diags`` on ``Rhea`` (which shares the same file system as ``Titan``, and ``Rhea`` is more suitable for data analysis).

1. Log on to ``rhea``:

::

    ssh -Y rhea.ccs.ornl.gov

2. If you don't have Anaconda installed, follow `this
guide <https://docs.continuum.io/anaconda/install-linux>`__.

3. Make sure you are using ``bash``

::

    bash

Installing and creating an environment
--------------------------------------
The steps below detail how to create your own environment with ``e3sm_diags``.
However, it is possible to use the `E3SM Unified Environment <https://acme-climate.atlassian.net/wiki/spaces/EPWCD/pages/374407241/E3SM+Unified+Environment>`__ instead.
If you decide to use the unified environment, please do so and skip to step 5 (Note, as of April 4th, 2018, the lastest version of unified-env is 1.1.3).

4. Allow Anaconda to download packages, even with a firewall.

::

    conda config --set ssl_verify false
    binstar config --set verify_ssl False

5. Update Anaconda.

::

    conda update conda

6. Get the yml file to create an environment.

::

    wget https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/conda/e3sm_diags_env.yml

7. Remove any cached Anaconda packages. This will ensure that you always get the latest packages.

::

    conda clean --all

8. Use Anaconda to create a new environment with ``e3sm_diags`` installed.
**Tip:** You can change the name of the environment by adding ``-n new_env_name`` to the end of ``conda env create ...``.

::

    conda env create -f e3sm_diags_env.yml
    source activate e3sm_diags_env


Running all sets of diagnostics
-------------------------------------------------

9. Copy and paste the below code into ``myparams.py`` using your
favorite text editor. Adjust any options as you like.

.. code:: python


    reference_data_path = '/ccs/proj/cli115/e3sm_diags_data/obs_for_e3sm_diags/climatology/'
    test_data_path = '/ccs/proj/cli115/e3sm_diags_data/test_model_data_for_acme_diags/'

    test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

    #sets = ["lat_lon"] # To run only lat_lon countour diags 
    #without specifing sets, sets = ['zonal_mean_xy', 'zonal_mean_2d', 'lat_lon', 'polar', 'cosp_histogram'] 

    diff_title = 'Model - Obs'

    backend = 'mpl'  # 'vcs' is for vcs plotting backend.

    results_dir = 'lat_lon_demo'  # name of folder where all results will be stored


10. Run the diags.

::

    e3sm_diags -p myparams.py


11. Open the following webpage to view the results.

::

    firefox --no-remote lat_lon_demo/viewer/index.html &

-  The ``--no-remote`` option uses the Firefox installed on this machine,
   and not the one on your machine.

Above example runs the package in serial. Instead, it can be run much faster using multi-processing, either in an interactive session on compute nodes, or as a batch
job.

Adding below lines to ``myparams.py``:

.. code:: python

    multiprocessing = True
    num_workers = 16  # Number of processes to use

**Tip:** Once you're on the webpage for a specific plot, click on the 'Output Metadata' 
drop down menu to view the metadata for the displayed plot.

* Running that command allows the displayed plot to be recreated. Changing any of the options will modify the resulting figure.


Interactive session on compute nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, request an interactive session with a single node (16 cores) for an half hour
(running this example should take much less than this): ::


  qsub -I -A charging_project_name -q name_of_queue -V -l nodes=1 -l walltime=00:30:00

Once the session is available, launch E3SM Diags: ::

  source activate e3sm_diags_env
  e3sm_diags -p myparams.py

Batch job
^^^^^^^^^

Alternatively, you can also create a script and submit it to the batch system.
Copy and paste the code below into a file named ``diags.pbs`` and **change the following**:

* Change ``charging_project_name`` to a valid value
* Change ``$YOUR_WORKING_DIR`` to your working directory
* Get the path of your Anaconda binary

  * Run ``which conda``, and get a path like so:
    ``/ccs/home/zhang40/anaconda3/envs/e3sm_diags_env/bin/conda``
  * Copy everything from the beginning to 'anaconda2' (or 'anaconda3') put it in:
    ``export PATH="PASTE_HERE/bin:$PATH"``

    An example path is:

    ``export PATH="/ccs/home/zhang40/anaconda3/bin:$PATH"``

.. code:: bash

  #!/bin/bash -l
  # PLEASE CHANGE: charging_project_name
  #PBS -A charging_project_name
  #PBS -N e3sm_diags_test
  #PBS -j oe
  #PBS -l walltime=0:30:00,nodes=1
 
  # PLEASE CHANGE: the line below to your valid path
  export PATH="/ccs/home/zhang40/anaconda3/bin:$PATH"
  source activate e3sm_diags_env
  # PLEASE CHANGE: $YOUR_WORKING_DIR to a valid directory
  cd $YOUR_WORKING_DIR
  e3sm_diags -p myparams.py

And then submit it ::

  qsub diags.pbs

