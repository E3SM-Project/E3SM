Acme1 quick guide for running e3sm_diags v2
=========================================================================

1. Installation
-----------------------------------------------------------

We will use the e3sm_unifed environment to install.
For the latest stable release or if you don't have access to e3sm analysis machines,
please instead refer to :ref:`Latest stable release <install_latest>`.

Most of the E3SM analysis software is maintained with an Anaconda metapackage
(E3SM unified environment).
If you have an account on Acme1,
then to get all of the tools in the metapackage in your path,
use the activation command below.
(Change ``.sh`` to ``.csh`` for csh shells.)

Below, we also provide the paths for observational data needed by ``e3sm_diags`` (<obs_path>),
and some sample model data for testing (<test_data_path>).
Both <obs_path> and <test_data_path> have two subdirectories:
``/climatology`` and ``/time-series`` for climatology and time-series data respectively.

Also listed below are paths where the HTML files (<html_path>) must be located to be displayed
at their corresponding web addresses (<web_address>).

<activation_command>: ``source /usr/local/e3sm_unified/envs/load_latest_e3sm_unified_acme1.sh``

<obs_path>: ``/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/``

<test_data_path>: ``/p/user_pub/e3sm/e3sm_diags_data/test_model_data_for_acme_diags/``

<html_path>: ``/var/www/acme/acme-diags/<username>/``

<web_address>: ``https://acme-viewer.llnl.gov/<username>/``
     


2. Config and run
--------------------------------------------------------

.. _Acme1_lat_lon:

Running the annual mean latitude-longitude contour set
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Copy and paste the below code into ``run_e3sm_diags.py`` using your favorite text editor.
Adjust any options as you like.

   **Tip:** Some of E3SM's analysis machines (**Acme1, Anvil, Compy, Cori**)
   have web servers setup to host html results.
   On Acme1,
   create the directory ``/var/www/acme/acme-diags/<username>/`` using your username.
   Set ``results_dir`` to ``/var/www/acme/acme-diags/<username>/doc_examples/lat_lon_demo``
   in ``run_e3sm_diags.py`` below. Then, you can view results via a web browser here:
   https://acme-viewer.llnl.gov/<username>/doc_examples/lat_lon_demo


    .. code:: python

        import os
        from e3sm_diags.parameter.core_parameter import CoreParameter
        from e3sm_diags.run import runner

        param = CoreParameter()

        param.reference_data_path = '/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/climatology/'
        param.test_data_path = '/p/user_pub/e3sm/e3sm_diags_data/test_model_data_for_acme_diags/climatology/'
        param.test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
        param.seasons = ["ANN"]   #all seasons ["ANN","DJF", "MAM", "JJA", "SON"] will run,if comment out"

        prefix = '/var/www/acme/acme-diags/<username>/doc_examples/'
        param.results_dir = os.path.join(prefix, 'lat_lon_demo')
        # Use the following if running in parallel:
        #param.multiprocessing = True
        #param.num_workers = 32
        
        # Use below to run all core sets of diags:
        #runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d']
        # Use below to run lat_lon map only:
        runner.sets_to_run = ['lat_lon']
        runner.run_diags([param])


Run in serial with:

    ::

        python run_e3sm_diags.py

The above run has the same results as running ``e3sm_diags -p lat_lon_params.py``
using the code below for ``lat_lon_params.py``:


    .. code:: python

        reference_data_path = '/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/climatology/'
        test_data_path = '/p/user_pub/e3sm/e3sm_diags_data/test_model_data_for_acme_diags/climatology/'

        test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

        sets = ["lat_lon"]
        seasons = ["ANN"]

        # 'mpl' for matplotlib plots
        backend = 'mpl'

        # Name of folder where all results will be stored.
        results_dir = '/var/www/acme/acme-diags/<username>/doc_examples/lat_lon_demo'

The new way of running (no ``-p``) is implemented in version 2.0.0,
preparing ``e3sm_diags`` to accomodate more diagnostics sets with set-specific parameters.



View results on the web
'''''''''''''''''''''''
Once the run is completed,
open  ``https://acme-viewer.llnl.gov/<username>/doc_examples/lat_lon_demo/viewer/index.html`` to view the results.
If you don't see the results, you may need to set proper permissions.
Run ``chmod -R 755 /var/www/acme/acme-diags/<username>/``.

**Tip:** Once you're on the webpage for a specific plot, click on the
'Output Metadata' drop down menu to view the metadata for the displayed plot.
Running that command allows the displayed plot to be recreated.
Changing any of the options will modify just that resulting figure.



Running all the core diagnostics sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Core diagnostics set includes:
**lat_lon**, **zonal_mean_xy**, **zonal_mean_2d**, **polar**, **cosp_histogram**,
**meridional_mean_2d**.
These diags share a common parameter space (core parameters).
To run all these sets without defining set-specific parameters
(e.g. **plev** for **zonal_mean_2d** and **meridional_mean_2d**.),
replace the ``runner.sets_to_run`` line in ``run_e3sm_diags.py`` with the one below:

 ::

   runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d']


Running area mean time series set
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In v2.0.0, the time series set was implemented to support regional averaged time series plotting
using monthly mean time series input.
This set is enabled if monthly mean time series is processed as documented
:doc:`here <../input-data-requirement>`.

A ``run_e3sm_diags.py`` example for running area mean time series alone:

    .. code:: python

        import os
        from e3sm_diags.parameter.core_parameter import CoreParameter
        from e3sm_diags.parameter.area_mean_time_series_parameter import AreaMeanTimeSeriesParameter
        from e3sm_diags.run import runner
        
        param = CoreParameter()
        
        param.reference_data_path = '/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/time-series/'
        param.test_data_path = '/p/user_pub/e3sm/e3sm_diags_data/test_model_data_for_acme_diags/time-series/E3SM_v1/'
        param.test_name = 'e3sm_v1'
        
        prefix = '/var/www/acme/acme-diags/<username>/doc_examples/'
        param.results_dir = os.path.join(prefix, 'area_mean_with_obs')
        # Use the following if running in parallel:
        #param.multiprocessing = True
        #param.num_workers =  40
        
        # We're passing in this new object as well, in
        # addition to the CoreParameter object.
        
        ts_param = AreaMeanTimeSeriesParameter()
        #ts_param.ref_names = ['none']   # Using this setting will plot only the model data, not the observation data
        ts_param.start_yr = '2002'
        ts_param.end_yr = '2008'
        
        runner.sets_to_run = ['area_mean_time_series']
        runner.run_diags([param, ts_param])


This set can also be ran with the core diagnostics sets,
so that all the plots are shown in one viewer.
The following is an example to run all sets:

    .. code:: python

        import os
        from e3sm_diags.parameter.core_parameter import CoreParameter
        from e3sm_diags.parameter.area_mean_time_series_parameter import AreaMeanTimeSeriesParameter
        from e3sm_diags.run import runner
        
        param = CoreParameter()
        
        param.reference_data_path = '/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/climatology/'
        param.test_data_path = '/p/user_pub/e3sm/e3sm_diags_data/test_model_data_for_acme_diags/climatology/'
        param.test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
        param.multiprocessing = True
        param.num_workers = 40
        prefix = '/var/www/acme/acme-diags/<username>/doc_examples'
        param.results_dir = os.path.join(prefix, 'all_sets')
        
        #
        ##Set specific parameters for new sets
        ts_param = AreaMeanTimeSeriesParameter()
        ts_param.reference_data_path = '/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/time-series/'
        ts_param.test_data_path = '/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/time-series/E3SM_v1/'
        ts_param.test_name = 'e3sm_v1'
        ts_param.start_yr = '2002'
        ts_param.end_yr = '2008'
        
        runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d', 'area_mean_time_series']
        runner.run_diags([param, ts_param])


Advanced: Running custom diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following steps are for 'advanced' users, who want to run custom diagnostics.
So, most users will not run the software like this.


By default, with ``e3sm_diags``,
a built in set of variables are defined for each diagonostics sets.
To do a short run, e.g. only running through a subset of variables,
a configuration file is needed to customize the run.


In the following example,
only precipitation and surface sea temperature are run to compare with
model and obs for lat_lon set.
Create ``mydiags.cfg`` file as below.

Check :doc:`Available Parameters <../available-parameters>` for all available parameters.

For a larger configuration file example, look
`here <https://github.com/E3SM-Project/e3sm_diags/blob/master/e3sm_diags/driver/default_diags/lat_lon_model_vs_obs.cfg>`_
for the cfg file that was used to create all of the latitude-longitude sets.


    ::

        [#]
        sets = ["lat_lon"]
        case_id = "GPCP_v2.3"
        variables = ["PRECT"]
        ref_name = "GPCP_v2.3"
        reference_name = "GPCP"
        seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
        regions = ["global"]
        test_colormap = "WhiteBlueGreenYellowRed.rgb"
        reference_colormap = "WhiteBlueGreenYellowRed.rgb"
        diff_colormap = "BrBG"
        contour_levels = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16]
        diff_levels = [-5, -4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4, 5]


Run E3SM diagnostics with the ``-d`` parameter.
Use the :ref:`above run script <Acme1_lat_lon>`. And run as following:

    ::

        python run_e3sm_diags.py -d mydiags.cfg


