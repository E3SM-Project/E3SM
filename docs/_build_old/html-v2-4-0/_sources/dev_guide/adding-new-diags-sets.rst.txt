How to Add Diagnostic Sets to E3SM Diagnostics
==============================================

Introduction
------------

When planning to expand the current diagnostic sets, efforts have
been made to structure the ``e3sm_diags`` code base in a more modularized
and expandable fashion. The developers will be able to add new
diagnostic sets easily and cleanly. This guide documents the essential
steps for guiding users to add custom diagnostics by providing
an example. If a user has a streamlined Python script for a
complete analysis, from reading in files, data manipulation,
computation, and visualization, it should be straightforward
to take the below steps to add the analysis to ``e3sm_diags``.

In this document, we will explain by example the process
of adding a diagnostics set into ``e3sm_diags``.
Some of the current diagnostics sets in ``e3sm_diags`` are
`the different kinds of plots seen here <../../../sample_output/modTS_vs_modTS_3years/viewer/index.html>`_.
If you have any questions or issues regarding this,
please make a Github issue on the e3sm_diags repo.

First off, to develop, you must have an 
:ref:`e3sm_diags development environment installed <dev-env>`.

Below are the components needed to add new diags:

1. **Parameters:** A set of variables that a user uses to define things pertaining
to a diagnostics run. These can be anything, from the path of the reference/test
data to parameters related to the plots created. As a reference, 
:doc:`look here for existing parameters and their description <../available-parameters>`.

2. **Driver:** The main code which takes a set of parameters and does the
diagnostics, including reading in data files, manipulating data sets, computing metrics,
and getting the data ready to be carried over to the plotting function.

3. **Adding Default Diagnostics:** Though users can choose what variables they
want to run diags on, we need to provide default variables.
Adding in derived variables is also covered in this guide.

4. **Plot:** A script which takes the output from the driver
and creates the plots.

5. **Viewer:** For a given plotset, create the htmls that host
the output of the plotting.


The Example Diagnostics Set
---------------------------

Say we have a diagnostics set that will simply take the difference between
some reference data and some test data (ex. regular lat-lon grid data in
climatology annual mean). Let's call it ``diff_diags``.
**From this point on,** `'diff_diags'` **is the name of this diagnostics set.**


Adding In The Parameters
------------------------

All of diagnostics sets in ``e3sm_diags`` share a set of core parameters,
which can be changed during runtime. The **default values** for these
parameters are defined in the
`acme_diags/parameter/core_parameter.py <https://github.com/E3SM-Project/e3sm_diags/blob/master/acme_diags/parameter/core_parameter.py>`_
folder.
:doc:`On the documentation website, there is more information explaining what each one does <../available-parameters>`.

First, open ``acme_diags/parameter/core_parameter.py`` and edit ``self.sets`` to include ``'diff_diags'``.

    .. code:: python

        self.sets = ['zonal_mean_xy', 'zonal_mean_2d', 'meridional_mean_2d',
            'lat_lon', 'polar', 'area_mean_time_series', 'cosp_histogram',
            'diff_diags',]

This ensures that when the user wants to run all of the plotsets,
the new one you made is included.

Now to make things more interesting, we want to add parameters
that are **just specific to this plotset**. Say we have parameters
called ``projection`` and ``central_lon``.

* ``projection`` changes the plot the that specific projection.
  The only values we support for now are ``'mercator'``, ``'platecarree'``, or ``'miller'``.
  The default value is ``'platecarree'``.
  For more information regarding these projections, `see here <https://scitools.org.uk/cartopy/docs/latest/crs/projections.html>`_.
* ``central_lon`` changes the central point of the plot. It's ``180`` be default.

In the ``acme_diags/parameter`` directory, create a file called ``diff_diags_parameter.py``.
Below is the code for it.

**When you create your own plotset, please change the name of the class.**
In this case, the name of the class in ``DiffDiagsParameter``.

    .. code:: python

        from .core_parameter import CoreParameter


        class DiffDiagsParameter(CoreParameter):
            def __init__(self):
                super().__init__()
                self.projection = 'platecarree'
                self.central_lon = 180

            def check_values(self):
                # Make sure that the core parameters are also valid.
                super().check_values()
                # Check that the user inputted values are valid.
                valid_projections = [
                    'mercator', 'platecarree',
                    'miller'
                ]
                if self.projection not in valid_projections:
                    msg = "Your projection ({}) isn't a valid value in: {}"
                    raise RuntimeError(msg.format(self.projection, valid_projections))
                if not (0 <= self.central_lon <= 360):
                    raise RuntimeError('central_lon must be between 0 and 360 inclusive.')

Note that we have a definition for the ``check_values()`` function.
This function makes sure that the user has inputted correct values for the parameters.
**This function is optional and isn't required.**

Also note that this Parameter class inherits from the ``CoreParameter``,
which contains the core parameters that all plotsets use. This means that
**this plotset will use all of the core parameters, in addition to the**
**set-specific ones define in the** ``DiffDiagsParameter``.


Letting The Parameter Class Be Used
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Open ``acme_diags/parameter/__init__.py`` and edit it like so.
Please read the comments.

    .. code:: python

        from .core_parameter import CoreParameter
        from .zonal_mean_2d_parameter import ZonalMean2dParameter
        # First, import the new Parameter class.
        from .diff_diags_parameter import DiffDiagsParameter

        SET_TO_PARAMETERS = {
            'zonal_mean_xy': CoreParameter,
            'zonal_mean_2d': ZonalMean2dParameter,
            'meridional_mean_2d': CoreParameter,
            'lat_lon': CoreParameter,
            'polar': CoreParameter,
            'cosp_histogram': CoreParameter,
            'area_mean_time_series': CoreParameter,
            # For the diff_diags plotset, we want to use the
            # below Parameter class.
            'diff_diags': DiffDiagsParameter,
        }

Adding In The Parser
^^^^^^^^^^^^^^^^^^^^

Notice that in the `acme_diags/parameter <https://github.com/E3SM-Project/e3sm_diags/tree/master/acme_diags/parameter>`_
folder, we have the following:

* ``core_parameter.py``, where ``CoreParameter`` is located.
* ``zonal_mean_2d_parameter.py``, where ``ZonalMean2dParameter`` is located.

  * This container parameters specific to the ``'zonal_mean_2d'`` plotset.
  * Again, since ``ZonalMean2dParameter`` inherits from ``CoreParameter``,
    it'll use all of the core parameters, in addition to it's set-specific ones.


**Every Parameter object needs a corresponding Parser object.**
The Parser is the command line parser which can take in the parameters
as command line arguments. The users have the option to run ``e3sm_diags``
via the command line. For example, this is done for the provenance, which
is the command shown in the bottom for each webpage with a plot.
One can reproduce or fine tune a single diagnostic by using the provenance command line.

In the ``acme_diags/parser`` directory, create a file called ``diff_diags_parser.py``.
Below is the code for it. Some points:

* Like how the ``DiffDiagsParameter`` inherits from ``CoreParameter``, our parser
  ``DiffDiagsParser`` inherits from ``CoreParser``.
* Remember to change the name of your class accordingly based on your plotset name.
* Please read the comments as well and implement the changes.

    .. code:: python

        from .core_parser import CoreParser
        # We need to import the corresponding Parameter object.
        from acme_diags.parameter.diff_diags_parameter import DiffDiagsParameter


        class DiffDiagsParser(CoreParser):
            def __init__(self, *args, **kwargs):
                if 'parameter_cls' in kwargs:
                    super().__init__(*args, **kwargs)
                else:
                    # We want this Parser to create objects of type DiffDiagsParameter.
                    super().__init__(parameter_cls=DiffDiagsParameter, *args, **kwargs)


            def load_default_args(self, files=[]):
                # This has '-p' and '--parameter' reserved.
                super().load_default_args(files)

                # The parameters unique to DiffDiagsParameter are added here.
                # For more information about adding arguments to the parser,
                # please search how to add arguments to Python's argument parser.
                # The way is exactly the same.
                self.add_argument(
                    '--projection',
                    type=str,
                    dest='projection',
                    help='The type of the projection '
                    + 'used when plotting.',
                    required=False)

                self.add_argument(
                    '--central_lon',
                    type=float,
                    dest='central_lon',
                    help='The central longitude of the plot.',
                    required=False)


Letting The Parser Class Be Used
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Open ``acme_diags/parser/__init__.py`` and edit it like so.
Please read the comments.

    .. code:: python

        from .core_parser import CoreParser
        from .zonal_mean_2d_parser import ZonalMean2dParser
        # First, import the new Parser class.
        from .diff_diags_parser import DiffDiagsParser


        SET_TO_PARSER = {
            'zonal_mean_xy': CoreParser,
            'zonal_mean_2d': ZonalMean2dParser,
            'meridional_mean_2d': CoreParser,
            'lat_lon': CoreParser,
            'polar': CoreParser,
            'cosp_histogram': CoreParser,
            'area_mean_time_series': CoreParser,
            # For the diff_diags plotset, we want to use the
            # below Parser class.
            'diff_diags': DiffDiagsParser,
        }


Adding In The Driver
--------------------

The driver is the main code which takes in a single Parameter object and does the diagnostics.
For each plotset, its corresponding driver is located in the
`acme_diags/driver <https://github.com/E3SM-Project/e3sm_diags/tree/master/acme_diags/driver>`_
folder. Please refer to these existing drivers, and if you need help creating
your driver, create a Github issue.

**This part of the code varies greatly based on the analysis.**
**There's no set way to do this.**

However, to get a variable based on the user's parameters (``reference_data_path``,
``test_data_path``, and more) and the way ``e3sm_diags`` input data is structured
(how the obs are named, the file naming conventions of the model files, etc.)
using the ``Dataset`` class is **highly recommended**.

* The ``diff_diags_driver.py`` below uses it.
* It's located in `acme_diags/driver/utils/dataset.py <https://github.com/E3SM-Project/e3sm_diags/blob/master/acme_diags/driver/utils/dataset.py>`_.
* With only two lines of code, here's how you get the variable PRECT from the test data with ANN climatology ran on it.

    .. code:: python

        test_data = utils.dataset.Dataset(parameter, test=True)
        prect_climo = test_data.get_climo_variable('PRECT', 'ANN')
* You can also get time-series data as well:

    .. code:: python

        test_data = utils.dataset.Dataset(parameter, test=True)
        prect_time_series = test_data.get_timeseries_variable('PRECT')

In ``acme_diags/driver``, create a file called ``diff_diags_driver.py``.
**Each Driver must have a** ``run_diags()`` **function which takes in a single Parameters object.**
**It also must return that Parameters object as well at the end of all of the for-loops.**

    .. code:: python

        from acme_diags.driver import utils
        from acme_diags.metrics import min_cdms, max_cdms, mean
        # The below will be defined in a future section.
        from acme_diags.plot.cartopy import diff_diags_plot

        def run_diag(parameter):
            variables = parameter.variables
            seasons = parameter.seasons

            test_data = utils.dataset.Dataset(parameter, test=True)
            ref_data = utils.dataset.Dataset(parameter, ref=True)

            for season in seasons:
                for var in variables:
                    test_var = test_data.get_climo_variable(var, season)
                    ref_var = ref_data.get_climo_variable(var, season)

                    # Only needed because our viewer (the final step)
                    # displays this data.
                    parameter.viewer_descr[var] = getattr(test_var, 'long_name', var)

                    # Regrid towards the lower resolution of the two
                    # variables for calculating the difference.
                    # The regrid_tool and regrid_method have default values.
                    test_var_reg, ref_var_reg = utils.general.regrid_to_lower_res(
                        test_var, ref_var, parameter.regrid_tool, parameter.regrid_method)

                    diff = test_var_reg - ref_var_reg

                    # We want to compute some metrics to plot as well.
                    metrics = {
                        'min': float(min_cdms(diff)),
                        'max': float(max_cdms(diff)),
                        'mean': float(mean(diff))
                    }

                    # This part will be defined in a forthcoming section.
                    diff_diags_plot.plot(diff, var, season, metrics, parameter)

            # Don't forget this.
            return parameter



Adding In The Config File For The Default Diagnostics
------------------------------------------------------------

When the user selects a certain number of plotsets to run, their
parameters are combined with default parameters for that plot set.
These are defined in
`acme_diags/driver/default_diags/ <https://github.com/E3SM-Project/e3sm_diags/tree/master/acme_diags/driver/default_diags>`_.
For the each plotset, we have two default files, one for ``model_vs_model`` runs and one for ``model_vs_obs`` runs.
The type of file used is determined by the ``run_type`` parameter in ``CoreParameter``, which is ``'model_vs_obs'`` by default.

Create a file ``diff_diags_model_vs_obs.cfg`` in the directory with the below contents.

    ::

        [#]
        sets = ["diff_diags"]
        case_id = "GPCP_v2.2"
        variables = ["NEW_PRECT"]
        ref_name = "GPCP_v2.2"
        seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
        diff_levels = [-5, -4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4, 5]

        [#]
        sets = ["diff_diags"]
        case_id = "SST_CL_HadISST"
        variables = ["SST"]
        ref_name = "HadISST_CL"
        seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
        diff_levels = [-5, -4, -3, -2, -1, -0.5, -0.2, 0.2, 0.5, 1, 2, 3, 4, 5]

        [#]
        sets = ["diff_diags"]
        case_id = "CERES-EBAF-TOA-v2.8"
        variables = ["SOLIN"]
        ref_name = "ceres_ebaf_toa_v2.8"
        seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
        diff_levels = [-5, -4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4, 5]

You must make sure the ``sets`` parameter in each section (a section starts with ``[#]``) is ``["diff_diags"]``.

In the above file, we have three sections. The result of this are three ``DiffDiagsParameter`` objects.
The user-inputted parameters will be added to each of these three objects.
Once combined with the user's input, each of the Parameter objects has the valid parameters to run the software.
Each is passed in to the ``run_diags()`` function in ``diff_diags_driver.py`` automatically.

Adding In Derived Variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^

This part only needs to be done if new variables are added for you new plotset.

A set of variables are defined in e3sm_diags in
`acme_diags/derivations/acme.py <https://github.com/E3SM-Project/e3sm_diags/blob/master/acme_diags/derivations/acme.py>`_.
Notice that in the above file, we had a variable ``NEW_PRECT``.
It's a new and derived variable, composed of two or more variables.
In this case, ``NEW_PRECT`` is composed of ``PRECT`` and ``PRECL``.
:doc:`Read more about derived variables here <../add-new-diagnostics>`.

Open `acme_diags/derivations/acme.py <https://github.com/E3SM-Project/e3sm_diags/blob/master/acme_diags/derivations/acme.py>`_
and add the function below. It handles what to do when we have ``NEW_PRECT`` as a variable.

    .. code:: python

        def new_prect(precc, precl):
            """
            Total precipitation flux = convective + large-scale.
            """
            var = precc + precl
            var = convert_units(var, "mm/day")
            var.long_name = "Total precipitation rate (convective + large-scale)"
            return var

We need to make sure that this function is actually called.
In the ``derived_variables`` dictionary in ``acme_diags/derivations/acme.py``, add the following entry.
It's basically the same as ``PRECT``, but with the ``new_prect()`` function being called.

    .. code:: python

        derived_variables = {
            'NEW_PRECT': OrderedDict([
                # This variable is 'PRECT' in newer versions of the obs data.
                # So we just get the variable and don't do anything.
                (('PRECT',), lambda prect: prect),
                # This variable is 'pr' in older versions of the obs data.
                (('pr',), lambda pr: qflxconvert_units(rename(pr))),
                # In the model data, it's composed of PRECC and PRECL.
                (('PRECC', 'PRECL'), lambda precc, precl: new_prect(precc, precl))
            ]),
            # Below is the old stuff. Don't insert the below.
            'PRECT': OrderedDict([
                (('pr',), lambda pr: qflxconvert_units(rename(pr))),
                (('PRECC', 'PRECL'), lambda precc, precl: prect(precc, precl))
            ]),


Installing This File
^^^^^^^^^^^^^^^^^^^^

Open ``setup.py`` in the root of this directory and add the following, somewhere before ``data_files`` is initialized.

    ::

        diff_diags_files = get_all_files_in_dir('acme_diags/driver/default_diags', 'diff_diags*')

Now modify ``data_files`` like below.

    .. code:: python

        (os.path.join(INSTALL_PATH, 'area_mean_time_series'),
        area_mean_time_series
        ),
        # We added in the below.
        (os.path.join(INSTALL_PATH, 'diff_diags'),
        diff_diags_files
        ),
        # The above was added in.
        (INSTALL_PATH,
        ['acme_diags/driver/acme_ne30_ocean_land_mask.nc',
        'misc/e3sm_logo.png'
        ])

**Every time you make a change to** ``diff_diags_model_vs_obs.cfg`` **and do a run, make sure you run**
``pip install .`` **as explained in the "Putting It All Together And Running On Cori At NERSC" section.**
**This is because when running the software, these files are obtained from where they're installed.**


Adding In The Plotting
----------------------

Data created in the driver needs to be carried over and needs to be plotted
using the plotting script. In our case, we only want to plot the ``diff`` data.
Other plotsets, like ``lat_lon``, plot more data.

In ``acme_diags/plot/cartopy/``, create a file called ``diff_diags_plot.py``.

The ``plot()`` function is what the driver, ``diff_diags_driver.py`` will call.
**Again, the code for this varies greatly based on the actual plot set.**
In this script, we're using the ``projection`` and ``central_lon`` parameters.

    .. code:: python

        import os
        import numpy as np
        import numpy.ma as ma
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        import cartopy.crs as ccrs
        from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
        from acme_diags.plot import get_colormap
        from acme_diags.driver.utils.general import get_output_dir


        def add_cyclic(var):
            lon = var.getLongitude()
            return var(longitude=(lon[0], lon[0] + 360.0, 'coe'))

        def get_ax_size(fig, ax):
            bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            width, height = bbox.width, bbox.height
            width *= fig.dpi
            height *= fig.dpi
            return width, height

        def plot(diff, var, season, metrics, parameter):
            # Create figure, projection
            fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)

            if parameter.projection == 'mercator':
                proj_cls = ccrs.Mercator
            elif parameter.projection == 'platecarree':
                proj_cls = ccrs.PlateCarree
            elif parameter.projection == 'miller':
                proj_cls = ccrs.Miller

            central_lon = parameter.central_lon
            proj = proj_cls(central_longitude=central_lon)

            diff = add_cyclic(diff)
            lon = diff.getLongitude()
            lat = diff.getLatitude()
            diff = ma.squeeze(diff.asma())

            # Contour levels
            clevels = parameter.diff_levels
            levels = None
            norm = None
            if len(clevels) > 0:
                levels = [-1.0e8] + clevels + [1.0e8]
                norm = colors.BoundaryNorm(boundaries=levels, ncolors=256)

            panel = (0.1691, 0.6810, 0.6465, 0.2258)
            # Contour plot
            ax = fig.add_axes(panel, projection=proj)
            # ax = fig.add_axes(panel[n], projection=proj)
            ax.set_global()
            cmap = get_colormap(parameter.diff_colormap, parameter)
            p1 = ax.contourf(lon, lat, diff,
                            transform=proj_cls(),
                            norm=norm,
                            levels=levels,
                            cmap=cmap,
                            extend='both',
                            )

            ax.set_aspect('auto')
            ax.coastlines(lw=0.3)
            if parameter.diff_title:
                ax.set_title(parameter.diff_title, fontdict={'fontsize': 11.5})
            ax.set_xticks([0, 60, 120, 180, 240, 300, 359.99], crs=proj_cls())
            ax.set_yticks([-90, -60, -30, 0, 30, 60, 90], crs=proj_cls())
            lon_formatter = LongitudeFormatter(
                zero_direction_label=True, number_format='.0f')
            lat_formatter = LatitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)
            ax.yaxis.set_major_formatter(lat_formatter)
            ax.tick_params(labelsize=8.0, direction='out', width=1)
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')

            # Color bar
            cbax = fig.add_axes(
                (panel[0] + 0.6635, panel[1] + 0.0215, 0.0326, 0.1792))
            cbar = fig.colorbar(p1, cax=cbax)
            w, h = get_ax_size(fig, cbax)

            if levels is None:
                cbar.ax.tick_params(labelsize=9.0, length=0)
            else:
                maxval = np.amax(np.absolute(levels[1:-1]))
                if maxval < 10.0:
                    fmt = "%5.2f"
                    pad = 25
                elif maxval < 100.0:
                    fmt = "%5.1f"
                    pad = 25
                else:
                    fmt = "%6.1f"
                    pad = 30
                cbar.set_ticks(levels[1:-1])
                labels = [fmt % l for l in levels[1:-1]]
                cbar.ax.set_yticklabels(labels, ha='right')
                cbar.ax.tick_params(labelsize=9.0, pad=pad, length=0)

            # Min, Mean, Max
            plotSideTitle = {'fontsize': 9.5}
            fig.text(panel[0] + 0.6635, panel[1] + 0.2107,
                    "Max\nMean\nMin", ha='left', fontdict=plotSideTitle)
            stats = metrics['min'], metrics['max'], metrics['mean']
            fig.text(panel[0] + 0.7635, panel[1] + 0.2107, "%.2f\n%.2f\n%.2f" %
                    stats[0:3], ha='right', fontdict=plotSideTitle)


            # Figure title
            if not parameter.main_title:
                fig.suptitle('{} {}'.format(var, season), x=0.5, y=0.96, fontsize=18)
            else:
                fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)

            # Save figure
            # Get the filename that the user has passed in and display that.
            # When running in a container, the paths are modified.
            file_name = '{}_{}.png'.format(var, season)
            path = os.path.join(get_output_dir('diff_diags', parameter,
                ignore_container=True), file_name)
            plt.savefig(path)
            print('Plot saved in: ' + path)
            plt.close()


Adding The Viewer
-----------------

Each plotset needs to have webpages generated for it that allow users to look at the resultant figures.
In ``e3sm_diags``, each of the plotset is **mapped to a function that takes in all of the Parameter objects**
**for that plotset, then creates the webpages and returns a** ``(display_name, url)`` **tuple of strings**.

First in `acme_diags/viewer/ <https://github.com/E3SM-Project/e3sm_diags/tree/master/acme_diags/viewer>`_
create a file ``diff_diags_viewer.py`` paste in the below code.

    .. code:: python

        import os
        from .utils import add_header, h1_to_h3
        from .default_viewer import create_metadata
        from cdp.cdp_viewer import OutputViewer


        def create_viewer(root_dir, parameters):
            """
            Given a set of parameters for a the diff_diags set,
            create a single webpage.

            Return the title and url for this page.
            """
            viewer = OutputViewer(path=root_dir)

            # The name that's displayed on the viewer.
            display_name = 'Diff Diagnostics'
            set_name = 'diff_diags'
            # The title of the colums on the webpage.
            cols = ['Description', 'Plot']
            viewer.add_page(display_name, short_name=set_name, columns=cols)
            viewer.add_group('Variable')

            for param in parameters:
                for var in param.variables:
                    for season in param.seasons:
                        viewer.add_row('{} {}'.format(var, season))
                        # Adding the description for this var to the current row.
                        # This was obtained and stored in the driver for this plotset.
                        viewer.add_col(param.viewer_descr[var])

                        file_name = '{}_{}.png'.format(var, season)
                        # We need to make sure we have relative paths, and not absolute ones.
                        # This is why we don't use get_output_dir() as in the plotting script
                        # to get the file name.
                        file_name = os.path.join('..', set_name, param.case_id, file_name)
                        viewer.add_col(file_name, is_file=True, title='Plot',
                            meta=create_metadata(param))

            url = viewer.generate_page()
            add_header(root_dir, os.path.join(root_dir, url), parameters)
            h1_to_h3(os.path.join(root_dir, url))

            return display_name, url


Now make sure that this ``create_viewer()`` function is actually called.
Open ``acme_diags/viewer/main.py`` and edit ``SET_TO_VIEWER``.

    .. code:: python

        # Import the newly created module.
        from . import diff_diags_viewer

        SET_TO_VIEWER = {
            'lat_lon': default_viewer.create_viewer,
            'polar': default_viewer.create_viewer,
            'zonal_mean_xy': default_viewer.create_viewer,
            'zonal_mean_2d': zonal_mean_2d_viewer.create_viewer,
            'meridional_mean_2d': default_viewer.create_viewer,
            'cosp_histogram': default_viewer.create_viewer,
            'area_mean_time_series': area_mean_time_series_viewer.create_viewer,
            # Add the below:
            'diff_diags': diff_diags_viewer.create_viewer,
        }


We use the CDP Viewer to create the webpages.
**This is not needed! Use whatever you want to create the webpages.**
Just make sure that your function that's mapped to the plotset in ``SET_TO_VIEWER``:

* Takes the ``root_dir`` (where the user wants the results outputted, so the ``results_dir`` parameter)
  and ``parameters`` (a list of Parameter objects for your plotset) as arguments.
* Returns a tuple ``(display_name, url)``, where ``display_name`` is the name
  displayed in the index and ``url`` is the URL of the webpage.


Putting It All Together And Running On Cori At NERSC
----------------------------------------------------

Go to the root of the repo where ``setup.py`` is located and run:

    ::

        pip install .


Some Examples To Run
^^^^^^^^^^^^^^^^^^^^
We are running ``e3sm_diags`` via the API. If you're not familar with running
the software that way, please see :doc:`this document <../examples/run-e3sm-diags-api>`.

Also, the paths to the reference and test data are for Cori at NERSC.
If on another machine, please change your paths accordingly.

Also, make sure to make sure your ``results_dir`` parameter is valid.


Running Without Any Set-Specific Parameters
"""""""""""""""""""""""""""""""""""""""""""

Call this script ``run_diff_diags_demo.py``.

    .. code:: python

        import os
        from acme_diags.parameter.core_parameter import CoreParameter
        from acme_diags.run import runner

        param = CoreParameter()

        param.reference_data_path = '/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/'
        param.test_data_path = '/global/cfs/cdirs/e3sm/acme_diags/test_model_data_for_acme_diags/climatology/'
        param.test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
        param.seasons = ["ANN"]
        prefix = '/global/cfs/cdirs/e3sm/www/shaheen2/runs_with_api'
        param.results_dir = os.path.join(prefix, 'diff_diags_demo')

        runner.sets_to_run = ['diff_diags']
        runner.run_diags([param])

Run it like so:

    ::

        python run_diff_diags_demo.py


Running With Set-Specific Parameters
""""""""""""""""""""""""""""""""""""

Call this script ``run_diff_diags_demo_specific.py``.

    .. code:: python

        import os
        from acme_diags.parameter.core_parameter import CoreParameter
        from acme_diags.parameter.diff_diags_parameter import DiffDiagsParameter
        from acme_diags.run import runner

        param = CoreParameter()

        param.reference_data_path = '/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/'
        param.test_data_path = '/global/cfs/cdirs/e3sm/acme_diags/test_model_data_for_acme_diags/climatology/'
        param.test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
        param.seasons = ["ANN"]
        prefix = '/global/cfs/cdirs/e3sm/www/shaheen2/runs_with_api'
        param.results_dir = os.path.join(prefix, 'diff_diags_demo_specific')

        # Set specific parameters.
        diff_diags_param = DiffDiagsParameter()
        diff_diags_param.projection = 'miller'
        diff_diags_param.central_lon = 30

        runner.sets_to_run = ['diff_diags']
        # We're passing in this new object as well, in
        # addtion to the CoreParameter object.
        runner.run_diags([param, diff_diags_param])

Run it like so:

    ::

        python run_diff_diags_demo_specific.py

The ticks on the plot won't match the newly changed ``central_lon`` but whatever.
This is just an example.


Running With Your Own Diags
"""""""""""""""""""""""""""

Create a file ``diags.cfg``.

    ::

        [#]
        sets = ["diff_diags"]
        case_id = "CERES-EBAF-TOA-v2.8"
        variables = ["ALBEDOC"]
        ref_name = "ceres_ebaf_toa_v2.8"
        seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
        diff_levels = [-0.25, -0.2, -0.15, -0.1, -0.07, -0.05, -0.02, 0.02, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25]


        [#]
        sets = ["diff_diags"]
        case_id = "CERES-EBAF-TOA-v2.8"
        variables = ["RESTOM"]
        ref_name = "ceres_ebaf_toa_v2.8"
        seasons = ["ANN"]
        diff_levels = [-50, -40, -30, -20, -10, -5, 5, 10, 20, 30, 40, 50]

        [#]
        sets = ["diff_diags"]
        case_id = "CERES-EBAF-TOA-v2.8"
        variables = ["FLUT"]
        ref_name = "ceres_ebaf_toa_v2.8"
        seasons = ["ANN"]
        diff_levels = [-50, -40, -30, -20, -10, -5, 5, 10, 20, 30, 40, 50]

Run it again with the same ``run_diff_diags_demo_specific.py`` previously defined.

    ::

        python run_diff_diags_demo_specific.py -d diags.cfg
