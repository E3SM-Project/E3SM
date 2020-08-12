.. _compass_ocean_isomip_plus:

ISOMIP+
=======

The second Ice Shelf-Ocean Model Intercomparison Project (ISOMIP+) is a set of
five experiments for ocean models with thermodynamically active ice shelf
cavities (i.e. with a parameterization of the sub-ice-shelf boundary layer
and of melt rates a the ice-ocean interface).  The experiments are described
in `Asay-Davis et al. (2016) <https://gmd.copernicus.org/articles/9/2471/2016/>`_.


ISOMIP+ is based on a earlier set of experiments,
`ISOMIP <https://ui.adsabs.harvard.edu/abs/2003AGUFM.C41A..05H/abstract>`_,
that was very idealized.  ISOMIP+ uses a more complex but still idealized ice
geometry.  It is a sister project of the third Marine Ice Sheet Model
Intercomparison Project
`(MISMIP+) <https://tc.copernicus.org/articles/14/2283/2020/tc-14-2283-2020.html>`_,
which dynamics of an ice sheet and ice shelf with significant variability in
two horizontal dimensions.  Topography for ISOMIP+ comes from MISMIP+
simulations with the BISICLES model and can be
`downloaded here <https://dataservices.gfz-potsdam.de/pik/showshort.php?id=escidoc:1487953>`_
(however, COMPASS will automatically download the topography, process it, and
regrid it to the MPAS mesh so users do not need to download it themselves).
ISOMIP+ uses far-field restoring over the full water column to force changes in
ocean conditions.  Restoring with a relatively warm far-field temperature
profile (the Ocean0 experiment) leads to a quasi-steady state within 1 to 2
year.  Two ISOMIP+ experiments prescribe dynamic topography, allowing models to
test their ability to handle moving boundaries and to see the effects that
moving topography has on ocean dynamics.

COM (2km) and TYP (5km) configurations
--------------------------------------

Each of the five experiments can be performed in a "common" (COM) or a "typical"
(TYP) setup.  The COM setup prescribes canonical horizontal (2 km) and vertical
(36 layers, ~20 m) resolutions as well as a linear equation of state, a
linear equation for the freezing point,and prescribed values for several
coefficients including: a root-mean-squared "tidal" velocity to include in the
computation of the friction velocity; top and bottom drag coefficients;
and horizontal and vertical viscosity and diffusivity.  The COM experiments also
prescribe a method for computing the value of the heat- and salt-transfer
coefficients (:math:`\Gamma_T` and :math:`\Gamma_S`) as part of the Ocean0
experiment described below.

The TYP experiments are intended to run at at resolutions typical for realistic
simulations performed with that ocean model and with typical parameter values.
The idea behind these experiments is to get a feel for how much these parameter
choices affect model spread compared with the COM choices.  For MPAS-Ocean, we
chose 5 km as our TYP resolution because, at the time, this was the target
resolution for our highest resolution simulations with ice-shelf cavities.
In practice, we have not run simulations with resolution this high.  But even
at 5 km resolution, the ISOMIP+ domain is not well resolved and the dynamics
are clearly affected by the poor resolution.  We are fully aware that we do not
capture small ice shelves like the one in ISOMIP+ in our typical MPAS-Ocean
simulations with ice-shelf cavities.

Ocean0
------

Ocean0 uses steady-state ice topography from the initial steady state of a
MISMIP+ Ice1 experiment (also described in
`Asay-Davis et al. (2016) <https://gmd.copernicus.org/articles/9/2471/2016/>`_).
produced with the BISICLES ice-sheet model. The ocean is initialized with a
relatively warm profile and restored the same profile in the far field. The
combination of warm initial conditions and restoring leads the system to reach a
quasi-equilibrium with strong melting over a few months to a few years. Because
Ocean0 is expected to reach a quasi-equilibrium within approximately 1 year,
this experiment is well suited to parameter studies. In particular, this
experiment is used to calibrate the values of the heat- and salt-transfer
coefficients, :math:`\Gamma_T` and :math:`\Gamma_S`, to achieve a target
melt rate of ~30 m/yr.


Ocean1
------

Ocean1 uses the same topography and restoring as Ocean0 but is initialized to a
colder, fresher profile that is expected to result in low melt rates during the
first several years of the simulation. Far-field restoring to the same warmer
profile as in Ocean0 leads to warmer and saltier water in the far field at
depth. The duration of the experiment is 20 years, which is typically sufficient
time to reach a quasi-steady state. Melt rates as well as the strengths of the
barotropic and overturning circulations toward the end of the simulation are
expected to be significantly larger than those within the first few years
because of the warming.

Ocean2
------

In Ocean2, the topography is from the end of the Ice1r simulation from BISICLES.
The ocean is initialized with the warm profile from Ocean0 and restored to the
cold profile from the initial condition of Ocean1. As with Ocean1, this
experiment runs for 20 years.

Ocean3 (not available)
----------------------

Ocean3 begins with the same topography as Ocean1, but in this experiment the
ice draft evolves over time according to a prescribed data set covering 100
years of ice retreat from Ice1r. Ocean3 is initialized and forced with warm
conditions, similar to Ocean0. Strong melting begins immediately as the
sub-ice-shelf circulation spins up, consistent with the conditions for Ice1r
used to generate the topography, and persist for the duration of the experiment.
Again, the topography is
`available here <https://dataservices.gfz-potsdam.de/pik/showshort.php?id=escidoc:1487953>`_.

MPAS-Ocean does not yet support the moving boundaries that this experiment
requires, so it has not yet been implemented in COMPASS.

Ocean4 (not available)
----------------------

Conceptually, Ocean4 is an extension of Ocean3. The ice- draft topography from
Ice1ra was produced by abruptly shutting off melting at year 100 and allowing
the ice to re-advance for 100 years. Thus, Ocean4 begins with the final
topography from Ocean3 (which is also the topography used in Ocean2). This time,
cold initial conditions and restorign are prescribed, leading to very low melt
rates, consistent with the lack of melting in the MISMIP+ run that produced the
ice topography.

Again, MPAS-Ocean does not yet support the moving boundaries, so this experiment
has not yet been implemented in COMPASS.

Setting up a run in COMPASS
---------------------------

In a local check-out of the ``MPAS-Dev/MPAS-Model/ocean/develop`` branch:

.. code-block:: bash

   cd testing_and_setup/compass/

To build test cases, you need to tell COMPASS where to find a few thing.
Open a file ``config.ocean`` and put the following, where we have used the
example path ``usr/projects/climate/username/mpas/model/ocean/develop`` as the
location where MPAS-Ocena has been checked out and compiled:

.. code-block:: ini

   # This file is the ocean core's configuration file. It is specific to the ocean
   # core, and a specific machine. Each machine will configure this file
   # differently, but it can be used to point on version of the testing
   # infrastructure at a different version of the model.


   # The namelists section defines paths to template namelists that will be used
   # to generate specific namelists. Typically these will point to the forward and
   # init namelists in the default_inputs directory after a successful build of
   # the ocean model.
   [namelists]
   forward = /usr/projects/climate/username/mpas/model/ocean/develop/namelist.ocean.forward
   init = /usr/projects/climate/username/mpas/model/ocean/develop/namelist.ocean.init


   # The streams section defines paths to template streams files that will be used
   # to generate specific streams files. Typically these will point to the forward and
   # init streams files in the default_inputs directory after a successful build of
   # the ocean model.
   [streams]
   forward = /usr/projects/climate/username/mpas/model/ocean/develop/streams.ocean.forward
   init = /usr/projects/climate/username/mpas/model/ocean/develop/streams.ocean.init


   # The executables section defines paths to required executables. These
   # executables are provided for use by specific test cases.
   # Full paths should be provided in order to access the executables from
   # anywhere on the machine.
   [executables]
   model = /usr/projects/climate/username/mpas/model/ocean/develop/ocean_model


   # The paths section describes paths that are used within the ocean core test
   # cases.
   [paths]

   # The mesh_database and the initial_condition_database are locations where
   # meshes / initial conditions might be found on a specific machine. They can be
   # the same directory, or different directory. Additionally, if they are empty
   # some test cases might download data into them, which will then be reused if
   # the test case is run again later.
   mpas_model = /usr/projects/climate/username/mpas/model/ocean/develop
   mesh_database = /usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/mesh_database
   initial_condition_database = /usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/initial_condition_database
   bathymetry_database = /usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/bathymetry_database

You can supply paths for the ``mesh_database``, ``initial_condition_database``
and ``bathymetry_database`` that are initially empty.  COMPASS will download
this data as needed.  If you are on an E3SM supported machine or LANL's
Institutional Computing, you should use the shared locations to avoid
unnecessary downloads.

List the available ISOMIP+ test cases:

.. code-block:: bash

   ./list_testcases.py -o ocean -c isomip_plus

There are 2 resolutions (2 km and 5 km) and 3 test cases at each resolution
(Ocean0, 1 and 2) plus a time-varying version of Ocean0 that is being used to
explore moving boundaries.  Pick the test case you are interested in running.
In this example, we choose COM (2km) Ocean0

Set up the test case as follows:

.. code-block:: bash

   ./setup_testcase.py -o ocean -c isomip_plus -r 2km -t Ocean0 -f config.ocean -m runtime_definitions/srun.xml --work_dir /lustre/scratch4/turquoise/username/isomip_plus_Ocean0

The directory ``/lustre/scratch4/turquoise/username/isomip_plus_Ocean0`` should
be some convenient path to set up and run the directory, typically in a scratch
space on a supercomputer or cluster.

Running the test case
---------------------

It is best to do a short test run (1 month) to make sure everything is working,
rather than jumping into a 2 or 20-year simulation.  How you do this will depend
on your system.  Here, we show starting an interactive job on LANL IC.  You will
need to figure out how set up your system with the appropriate compilers, MPI,
libraries, etc., which will be covered elsewhere in this documentation.

.. code-block:: bash

   cd /lustre/scratch4/turquoise/username/isomip_plus_Ocean0/ocean/isomip_plus/2km/Ocean0/
   salloc --nodes=1 --time=0:20:00 --account=e3sm

   <<load modules needed for MPI, python, etc.>>

   ./run_test.py

Running a full 2- or 20-year simulation
---------------------------------------

For a longer simulation, you will need a job script.

.. code-block:: bash

   cd /lustre/scratch4/turquoise/username/isomip_plus_Ocean0/ocean/isomip_plus/2km/Ocean0/forward
   vim job_script.bash

This is an example job script for Ocean0 on LANL IC:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --nodes=4
   #SBATCH --time=4:00:00
   #SBATCH --account=e3sm
   #SBATCH --job-name=Ocean0
   #SBATCH --output=Ocean0.o%j
   #SBATCH --error=Ocean0.e%j
   #SBATCH --qos=interactive

   # exit if there are any errors
   set -e

   <<load modules needed for MPI, python, etc.>>

   months_per_job=24
   end_date="0003-01-01_00:00:00"

   for month in `seq 0 $months_per_job`
   do
       ./check_progress.py -f namelist.ocean -e $end_date
       ./run.py
       ./setup_restart.py -f namelist.ocean
   done

For Ocean1 and Ocean2, the ``end_date`` should include the year ``0021``
(since the simulation starts at year ``0001``).

Submit the job:

.. code-block:: bash

   sbatch job_script.bash

Once it's running, monitor the progress with:

.. code-block:: bash

   tail log.ocean.0000.out

This writes a message for each time step (if all is going well).

The simulation runs one month at a time and then does some adjustment in a
python script to make sure sea level doesn't get out of control (there's a lot
of melting going on so we have to have a compensating level of "evaporation" at
the domain boundary).  It also will check to see if we've already reached year
2 and won't run again if so.

Some basic output is available in:

.. code-block:: none

   analysis_members/globalStats.0001-01-01_00.00.00.nc

To see the mean melt flux and how time is progressing there, do:

.. code-block:: bash

   ncdump -v xtime,landIceFreshwaterFluxAvg analysis_members/globalStats.0001-01-01_00.00.00.nc | tail -n 50

Keep in mind that the units are ``kg m^{-2} s^{-1}``, not m/yr, so not the most
intuitive output.  There is a ``viz`` package in the ``isomip_plus`` directory
that gets linked in each test case.  You can also look at output in paraview.

Visualization
-------------

There should be a link to the ``viz`` python package in the ``forward``
directory.  To run the ``viz`` package, you will need several dependencies as
part of a conda environment, including ``xarray``, ``dask``, ``numpy``,
``progressbar2``, ``netcdf4``, ``shapely``, ``matplotlib``, ``ffmpeg`` and
``scipy``.  These are all part of the ``compass`` environment defined by the
``compass`` package, so if you have already installed and activated the latest
``compass``, you're all set.  If you are on a supported machine, such as on
LANL IC, you can activate with a command like:

.. code-block:: bash

   source /usr/projects/climate/SHARED_CLIMATE/anaconda_envs/base/etc/profile.d/conda.sh
   conda activate compass_0.1.8

Once you have an appropriate conda environment, run:

.. code-block:: bash

   python -m viz

If you are running Ocean1 or Ocean2, you need to specify this.  Ocean0 is the
default experiment.  to see the other command-line flags, run:

.. code-block:: bash

   python -m viz --help

Running ``viz`` will take maybe 10 or 15 minutes (most of it on the overturning
streamfunction).  You should see something like:

.. code-block:: none

   barotropic streamfunction: 100% |##############################| Time:  0:00:15
   compute and caching transport on MPAS grid:
   [########################################] | 100% Completed |  7.2s
   interpolating tansport on z-level grid: 100% |#################| Time:  0:10:13
   caching transport on z-level grid:
   [########################################] | 100% Completed |  2.2s
   compute and caching vertical transport sum on z-level grid:
   [########################################] | 100% Completed |  2.4s
   bin overturning streamfunction: 100% |#########################| Time:  0:02:03
   plotting barotropic streamfunction: 100% |#####################| Time:  0:00:08
   plotting overturning streamfunction: 100% |####################| Time:  0:00:05
   plotting melt rate: 100% |#####################################| Time:  0:00:07
   plotting heat flux from ocean to ice-ocean interface: 100% |###| Time:  0:00:07
   plotting heat flux into ice at ice-ocean interface: 100% |#####| Time:  0:00:07
   plotting thermal driving: 100% |###############################| Time:  0:00:07
   plotting haline driving: 100% |################################| Time:  0:00:07
   plotting friction velocity: 100% |#############################| Time:  0:00:08
   plotting top temperature: 100% |###############################| Time:  0:00:09
   plotting bot temperature: 100% |###############################| Time:  0:00:08
   plotting temperature section: 100% |###########################| Time:  0:00:05
   plotting top salinity: 100% |##################################| Time:  0:00:08
   plotting bot salinity: 100% |##################################| Time:  0:00:08
   plotting salinity section: 100% |##############################| Time:  0:00:05
   plotting top potential density: 100% |#########################| Time:  0:00:10
   plotting bot potential density: 100% |#########################| Time:  0:00:08
   plotting potential density section: 100% |#####################| Time:  0:00:05
   running ffmpeg -y -r 30 -i ./plots/botPotRho/botPotRho_%04d.png -b:v 32000k -r 30 ./movies/botPotRho.mp4
   running ffmpeg -y -r 30 -i ./plots/botSalinity/botSalinity_%04d.png -b:v 32000k -r 30 ./movies/botSalinity.mp4
   running ffmpeg -y -r 30 -i ./plots/botTemp/botTemp_%04d.png -b:v 32000k -r 30 ./movies/botTemp.mp4
   running ffmpeg -y -r 30 -i ./plots/bsf/bsf_%04d.png -b:v 32000k -r 30 ./movies/bsf.mp4
   running ffmpeg -y -r 30 -i ./plots/frictionVelocity/frictionVelocity_%04d.png -b:v 32000k -r 30 ./movies/frictionVelocity.mp4
   running ffmpeg -y -r 30 -i ./plots/halineDriving/halineDriving_%04d.png -b:v 32000k -r 30 ./movies/halineDriving.mp4
   running ffmpeg -y -r 30 -i ./plots/iceHeatFlux/iceHeatFlux_%04d.png -b:v 32000k -r 30 ./movies/iceHeatFlux.mp4
   running ffmpeg -y -r 30 -i ./plots/meltRate/meltRate_%04d.png -b:v 32000k -r 30 ./movies/meltRate.mp4
   running ffmpeg -y -r 30 -i ./plots/oceanHeatFlux/oceanHeatFlux_%04d.png -b:v 32000k -r 30 ./movies/oceanHeatFlux.mp4
   running ffmpeg -y -r 30 -i ./plots/osf/osf_%04d.png -b:v 32000k -r 30 ./movies/osf.mp4
   running ffmpeg -y -r 30 -i ./plots/sectionPotRho/sectionPotRho_%04d.png -b:v 32000k -r 30 ./movies/sectionPotRho.mp4
   running ffmpeg -y -r 30 -i ./plots/sectionSalinity/sectionSalinity_%04d.png -b:v 32000k -r 30 ./movies/sectionSalinity.mp4
   running ffmpeg -y -r 30 -i ./plots/sectionTemp/sectionTemp_%04d.png -b:v 32000k -r 30 ./movies/sectionTemp.mp4
   running ffmpeg -y -r 30 -i ./plots/thermalDriving/thermalDriving_%04d.png -b:v 32000k -r 30 ./movies/thermalDriving.mp4
   running ffmpeg -y -r 30 -i ./plots/topPotRho/topPotRho_%04d.png -b:v 32000k -r 30 ./movies/topPotRho.mp4
   running ffmpeg -y -r 30 -i ./plots/topSalinity/topSalinity_%04d.png -b:v 32000k -r 30 ./movies/topSalinity.mp4
   running ffmpeg -y -r 30 -i ./plots/topTemp/topTemp_%04d.png -b:v 32000k -r 30 ./movies/topTemp.mp4

The more interesting results should be a series of movies in ``movies`` and 4
time series plots in ``plots`` (mean melt rate, total melt flux, mean thermal
driving and mean friction velocity) and the same plots in
``timeSeriesBelow300m``, but this time averaged only over the deepest part of
the ice shelf (where much of the action is).

More Information
----------------

.. toctree::
   :titlesonly:

   isomip_plus_at_lanl

