
Quick guide for ANL ALCF Cooley using Singularity (v1)
======================================================

Obtaining the container image
-----------------------------

Singularity is the container runtime used at ALCF to run containers.
It's an alternative to Docker, but supports Docker containers.

1. Set the ``SINGULARITY_PULLFOLDER`` environmental variable to ``/lus/theta-fs0/projects/ccsm/acme/containers/``.
This is where the containers are stored.
If you're using bash, run the following.

    ::

        export SINGULARITY_PULLFOLDER="/lus/theta-fs0/projects/ccsm/acme/containers/"


2. Unlike Docker or Shifter (on NERSC), Singularity doesn't store the containers is a central repository.
So we just the directory below to store all downloaded containers.

View the ``e3sm_diags`` images available on ALCF.

    ::

        ls $SINGULARITY_PULLFOLDER | grep e3sm_diags

If the version you want to use is already available, then please continue to step 3.
Otherwise, you'll need to download the image you want, shown in step 2.


3. If the specific version you want or the ``latest`` image **is not shown**, download it.
You can view all of the images available on the 
`e3sm_diags Docker Hub <https://hub.docker.com/r/e3sm/e3sm_diags/tags/>`_.

Below, we are getting the image with the ``latest`` tag.
Images are stored in the ``SINGULARITY_PULLFOLDER`` directory.

    ::

        singularity pull docker://e3sm/e3sm_diags:latest 


4. ``wget`` the following script.

    ::

        wget https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/acme_diags/container/e3sm_diags_container.py



Running the entire annual latitude-longitude contour set
--------------------------------------------------------

4. Copy and paste the below code into ``myparams.py`` using your favorite text editor.
Adjust any options as you like.

    .. code:: python

        reference_data_path = '/lus/theta-fs0/projects/ClimateEnergy_3/e3sm_diags/obs_for_e3sm_diags/climatology/'
        test_data_path = '/lus/theta-fs0/projects/ClimateEnergy_3/e3sm_diags/test_model_data_for_e3sm_diags/climatology/'

        test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

        sets = ["lat_lon"]
        seasons = ["ANN"]

        # 'mpl' and 'vcs' are for matplotlib or vcs plots respectively.
        backend = 'mpl'

        # Name of folder where all results will be stored.
        results_dir = 'lat_lon_demo'

        # Optional parameters to run in parallel.
        multiprocessing = True
        num_workers =  12

Singularity can be ran on the login nodes, but you should run
either in an interactive session on the compute nodes, or as a batch job.

Interactive session on compute nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, request an interactive session with a single node (12 cores) for one hour
(running this example should take much less than this) on the ``default`` queue.
If obtaining a session takes too long, try to use the ``debug`` queue.
For more information,
`see here <https://www.alcf.anl.gov/user-guides/job-scheduling-policies-cooley>`_.

    ::

        qsub --interactive --nodecount=1 --queue=default --time=01:00:00


Once the session is available, reset the environmental variable.
In bash, use the below.

    ::

        export SINGULARITY_PULLFOLDER="/lus/theta-fs0/projects/ccsm/acme/containers/"

Then launch E3SM Diagnostics:

    ::
    
        python e3sm_diags_container.py --singularity -p myparams.py

**Tip:** You can select the version of the container you want to run with the ``--container_version`` argument.
If this argument isn't defined, it defaults to the ``latest`` container.

    ::

        python e3sm_diags_container.py --shifter --container_version v1.6.0 -p myparams.py

Use the command ``exit`` when the run is done.


Batch job
^^^^^^^^^

Alternatively, you can also create a script and submit it to the batch system.

First, create a temporary directory and make sure you download the script into it.
The reason we're doing this first is because the compute nodes
don't seem to not have Internet access.

    .. code:: bash
    
        mkdir ~/e3sm_diags_output
        wget -P ~/e3sm_diags_output https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/acme_diags/container/e3sm_diags_container.py


Also, make sure your ``myparams.py`` is also in the directory you've made.

Finally, copy and paste the code below into a file named ``diags.sh``.
You might need to make it executable with ``chmod u+x diags.sh``.

    .. code::bash

        #!/bin/bash -l
        #COBALT --jobname=diags
        #COBALT --output=diags.o%j
        #COBALT --queue=default
        #COBALT --nodecount=1
        #COBALT --time=01:00:00

        export SINGULARITY_PULLFOLDER="/lus/theta-fs0/projects/ccsm/acme/containers/"
        cd ~/e3sm_diags_output
        python e3sm_diags_container.py --singularity -p myparams.py

And then submit it

    ::

        qsub diags.bash

View the status of your job with ``qstat -u <username>``.


Back to running the latitude-longitude contour set
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
5. Once you ran the diagnostics in an interactive session or via a batch job,
open the following webpage to view the results.

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

        reference_data_path = '/lus/theta-fs0/projects/ClimateEnergy_3/e3sm_diags/obs_for_e3sm_diags/climatology/'
        test_data_path = '/lus/theta-fs0/projects/ClimateEnergy_3/e3sm_diags/test_model_data_for_e3sm_diags/climatology/'

        test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

        # Not defining a sets parameter runs all of the default sets:
        # ['zonal_mean_xy', 'zonal_mean_2d', 'lat_lon', 'polar', 'cosp_histogram']
        # Not defining a seasons parameter runs all of the seasons:
        # ['ANN', 'DJF', 'MAM', 'JJA', 'SON']

        # 'mpl' and 'vcs' are for matplotlib or vcs plots respectively.
        backend = 'mpl'

        # Name of folder where all results will be stored.
        results_dir = 'diag_demo'

        # Optional settings below:
        diff_title = 'Model - Obs'

        multiprocessing = True
        num_workers =  12


Compared to the previous short test above, note the following changes:

* Plots for all the available sets ('zonal_mean_xy', 'zonal_mean_2d',
  'lat_lon', 'polar', 'cosp_histogram') are generated.
* Plots for all of the seasons ('ANN', 'DJF', 'MAM', 'JJA', 'SON') are generated.


6. Again, run the diagnostics with this new parameter file (``all_sets.py``), either
   in an interactive session or via a batch job.


7. Open the following webpage to view the results.

    ::

        diags_demo/viewer/index.html



Advanced: Running custom diagnostics
------------------------------------
The following steps are for 'advanced' users, who want to run custom diagnostics.
So most users will not run the software like this.


By default, all of the E3SM diagnostics are ran for the sets that we defined above.
This takes some time, so instead we create our own diagnostics to be ran.


8. Copy and paste the code below in ``mydiags.cfg``.
Check :doc:`Available Parameters <../../available-parameters>`
for all available parameters.

For more examples of these types of files, look
`here <https://github.com/E3SM-Project/e3sm_diags/blob/master/acme_diags/driver/default_diags/lat_lon_model_vs_obs.cfg>`_
for the ``cfg`` file that was used to create all of the latitude-longitude sets.


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

9. Run E3SM diagnostics with the ``-d`` parameter.

    ::

        python e3sm_diags_container.py --singularity -p myparams.py -d mydiags.cfg
