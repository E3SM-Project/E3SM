
Quick guide for NERSC Edison using Shifter (v1)
===============================================

Obtaining the container image
-----------------------------

Shifter is the container runtime used at NERSC to run containers.
It's an alternative to Docker, but supports Docker containers.

1. View the ``e3sm_diags`` images available at NERSC.

    ::

        shifterimg images | grep e3sm_diags


If the version you want to use is already available, then please continue to step 3.

Otherwise, you'll need to download the image you want, shown in step 2.


2. If the specific version you want or the ``latest`` image **is not shown**, download it.
You can view all of the images available on the 
`e3sm_diags Docker Hub <https://hub.docker.com/r/e3sm/e3sm_diags/tags/>`_.
Below, we are getting the image with the latest tag.

    ::

        shifterimg -v pull docker:e3sm/e3sm_diags:latest 

* You'll see the same message with a timestamp print multiple times.
  This is normal and takes around 10 minutes or so.
  Something's just wrong with Shifter, we don't know why it does that.
* Once an image is downloaded from a public repo like this one, all users on NERSC can use it.
* You also cannot delete an image that you downloaded via Shifter for now.
  Please email NERSC support and they can do that for you.


3. ``wget`` the following script.

    ::

        wget https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/acme_diags/container/e3sm_diags_container.py



Running the entire annual latitude-longitude contour set
--------------------------------------------------------

4. Copy and paste the below code into ``myparams.py`` using your favorite text editor. Adjust any options as you like.

   **Tip:** Make a folder in the following directory ``/global/project/projectdirs/acme/www/`` based off your username.
   Then you can set ``results_dir`` to  ``/global/project/projectdirs/acme/www/<username>/lat_lon_demo`` in ``myparams.py`` below
   to view the results via a web browser here: http://portal.nersc.gov/project/acme/<username>/lat_lon_demo

    .. code:: python

        reference_data_path = '/global/project/projectdirs/acme/acme_diags/obs_for_e3sm_diags/climatology/'
        test_data_path = '/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/climatology/'

        test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

        sets = ["lat_lon"]
        seasons = ["ANN"]

        # 'mpl' and 'vcs' are for matplotlib or vcs plots respectively.
        backend = 'mpl'

        # Name of folder where all results will be stored.
        results_dir = 'lat_lon_demo'

Since Shifter cannot be ran on the login nodes, it must be ran either in an
**interactive session** on compute nodes, or as a **batch job**.


Interactive session on compute nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, request an interactive session with a single node (24 cores) for one hour (running this example should take much less than this).

If obtaining a session takes too long, try to use the ``debug`` partition.
Note that the maximum time allowed for this partition is ``00:30:00``.

    ::

        salloc --nodes=1 --partition=regular --time=01:00:00 -C haswell


Once the session is available, launch E3SM Diagnostics.

    ::

        python e3sm_diags_container.py --shifter -p myparams.py

**Tip:** You can select the version of the container you want to run with the ``--container_version argument``.
If this argument isn't defined, it defaults to the ``latest`` container.

    ::

        python e3sm_diags_container.py --shifter --container_version v1.5.0 -p myparams.py


Batch job
^^^^^^^^^

Alternatively, you can also create a script and submit it to the batch system.
Copy and paste the code below into a file named ``diags.bash``.
Please remember to change what directory you're in to one accessible to you.

    .. code:: bash
    
        #!/bin/bash -l
        #SBATCH --job-name=diags
        #SBATCH --output=diags.o%j
        #SBATCH --partition=regular
        #SBATCH --account=acme
        #SBATCH --nodes=1
        #SBATCH --time=01:00:00

        # Please change the directory below.
        cd /global/cscratch1/sd/golaz/tmp
        wget https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/acme_diags/container/e3sm_diags_container.py
        python e3sm_diags_container.py --shifter -p myparams.py

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

**Tip:** Once you're on the webpage for a specific plot, click on the 'Output Metadata' drop down menu to view the metadata for the displayed plot.
Running that command allows the displayed plot to be recreated.
Changing any of the options will modify the just that resulting figure.



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

        # 'mpl' and 'vcs' are for matplotlib or vcs plots respectively.
        backend = 'mpl'

        # Name of folder where all results will be stored.
        results_dir = 'diag_demo'

        # Optional settings below:

        diff_title = 'Model - Obs'

        multiprocessing = True
        num_workers =  24


Compared to the previous short test above, note the following changes:

* Plots for all the available sets ('zonal_mean_xy', 'zonal_mean_2d',
  'lat_lon', 'polar', 'cosp_histogram') are generated.
* Multiprocessing with 24 workers is enabled.


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

9. Run E3SM diagnostics with the ``-d`` parameter.

    ::

        python e3sm_diags_container.py --shifter -p myparams.py -d mydiags.cfg


