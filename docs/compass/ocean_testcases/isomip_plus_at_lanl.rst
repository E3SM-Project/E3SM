.. _compass_ocean_isomip_plus_at_lanl:

Instructions for setting up and running ISOMIP+ Ocean0 on LANL IC
=================================================================

In what follows, replace ``username`` with your user name.

1. SSH tricks
-------------

A couple of tricks for your laptop if you’re not already using them:
Save your SSH connections:

.. code-block:: bash

   vim ~/.ssh/config

Add the following:

.. code-block::

   Host *
       ControlMaster auto
       ControlPath ~/.ssh/connections/%r@%h:%p
       ServerAliveInterval 300

   Host wtrw
       Hostname wtrw.lanl.gov
       User username

.. code-block:: bash

   mkdir ~/.ssh/connections

Alias connections to LANL HPC machines:

.. code-block:: bash

   vim ~/.bashrc

Add:

.. code-block:: bash

   alias gr="ssh -t wtrw ssh gr-fe1"
   alias ba="ssh -t wtrw ssh ba-fe2"
   alias ko="ssh -t wtrw ssh ko-fe"

2. Making sure git is set up nicely
-----------------------------------

2.1 Storing your LANL IC SSH key on GitHub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It's useful to set up GitHub to know your public SSH key from LANL IC if you
haven't already done this.  It means you don’t have to type your password for
GitHub each time you git fetch, git push, etc.

I believe this is the right link for
`more details <https://help.github.com/en/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent>`_
If you haven't done this already and this gives you trouble, let me know and we
can work through it together.

2.2 git settings
^^^^^^^^^^^^^^^^

On IC or on your laptop, make sure you’ve got these settings defined:

.. code-block:: bash

   git config --global user.name "First Last"
   git config --global user.email user@domain.com
   git config --global core.editor vim
   git config --global push.default nothing
   git config --global color.ui true
   git config --global core.autocrlf input
   git config --global core.whitespace trailing-space
   git config --global alias.logg "log --graph --oneline --decorate"

I use ``git logg`` all the time so this last alias is particularly important.

2.3 git tab completion
^^^^^^^^^^^^^^^^^^^^^^

Download `git-completion.bash <https://raw.githubusercontent.com/git/git/master/contrib/completion/git-completion.bash>`_

.. code-block:: bash

   cd ~
   wget https://raw.githubusercontent.com/git/git/master/contrib/completion/git-completion.bash

Add this to your .bashrc

.. code-block::

   module load git
   source git-completion.bash

3. Forking and Cloning MPAS-Model
---------------------------------


* Go to: `https://github.com/MPAS-Dev/MPAS-Model <https://github.com/MPAS-Dev/MPAS-Model>`_
* Make your own fork by clicking “Fork” at the top right:
* Go to your new fork (e.g. `https://github.com/username/MPAS-Model <https://github.com/username/MPAS-Model>`_ )
* Whenever you ever need to know the link to clone your fork

  * Click on “Clone or download”
  * If it says “Clone with HTTPS”, click Use SSH (either works but SSH will use
    the SSH keys you’ve set up above and you never have to type my Git
    password.)
  * Copy the link with the clipboard icon

In a terminal window, log in to a LANL machine (I use Grizzly from here on
except where stated):

.. code-block:: bash

   ssh -t wtrw ssh gr-fe1

Make a directory for the code, e.g.:

.. code-block:: bash

   mkdir /usr/projects/climate/username
   cd /usr/projects/climate/username
   mkdir -p mpas/model
   cd mpas/model/

Clone the repo:

.. code-block:: bash

   git clone git@github.com:username/MPAS-Model.git repo
   cd repo

Rename your remote so it’s easier to not confuse it with other forks:

.. code-block:: bash

   git remote rename origin username/MPAS-Model

Add the main repo:

.. code-block:: bash

   git remote add MPAS-Dev/MPAS-Model git@github.com:MPAS-Dev/MPAS-Model.git

Add my fork (you can add other people’s forks in the same way):

.. code-block:: bash

   git remote add xylar/MPAS-Model git@github.com:xylar/MPAS-Model.git

Get the latest version of all the remotes (pruning anything that has been
deleted):

.. code-block:: bash

   git fetch --all -p

Let's store some settings you'll need to load every time you build MPAS.  The following
are only appropriate for Grizzly and we'll need a similar file with settings for
Badge and any other machines we might use in the future.

.. code-block:: bash

   vim ../setup_gr.bash

In this file, put:

.. code-block:: bash

   echo "Setting up grizzly intel environment for building and running MPAS"
   module purge
   module load git
   source /usr/projects/climate/SHARED_CLIMATE/anaconda_envs/base/etc/profile.d/conda.sh
   conda activate compass_py3.7
   module use /usr/projects/climate/SHARED_CLIMATE/modulefiles/all/
   module load intel/17.0.1 openmpi/1.10.5 netcdf/4.4.1 parallel-netcdf/1.5.0 pio/1.7.2
   export CORE=ocean

4. Checking out an MPAS branch and building the model
-----------------------------------------------------

**Note: this is a good place to come back to when you need to start over on
a new branch.**

Add a "worktree", a copy of the repo that we can point to a different branch.
We will work with my branch ``ocean/update_isomip_plus_viz``\ , where I have added some
new viz tools.  This is based off of the latest ``ocean/develop``. In general,
``ocean/develop`` is the place to start, since the ``master``  branch is updated only
rarely when we make releases:

.. code-block:: bash

   cd /usr/projects/climate/username/mpas/model/reop

Let's make sure we have the latest version of all the branches on all of the remotes

.. code-block:: bash

   git fetch --all -p

Okay, now we're ready to make a new folder to work from.

.. code-block:: bash

   git worktree add ../ocean/update_isomip_plus_viz -b ocean/update_isomip_plus_viz
   cd ../ocean/update_isomip_plus_viz

Take a look at which branch were on:

.. code-block::

   git logg

We don't start off on ``MPAS-Dev/MPAS-Model/ocean/update_isomip_plus_viz`` (even though
the name of the local branch might trick you into thinking you're there), so we need
to do a hard reset to put us there:

.. code-block:: bash

   git reset --hard xylar/MPAS-Model/ocean/update_isomip_plus_viz
   git logg

Now source the file with modules and settings for building MPAS on grizzly:

.. code-block:: bash

   source ../../setup_gr.bash

If all goes well, you should see ``comapss_py3.7`` as part of your command prompt and you should be read to build MPAS.

.. code-block:: bash

   make ifort

Take a coffee break, this will take some time.
...

5. Setting up a test case
-------------------------

Okay you're back and refreshed?  Let's set up a test case.

.. code-block:: bash

   cd testing_and_setup/compass/

COMPASS (COnfiguration of Model for Prediction Across Scales Setups -- yes, a litle tortured) is a set of python
scripts we use to set up and run our test cases.  To build test cases, you need to tell COMPASS where to find a few
thing on Grizzly.  Open a file ``config.ocean`` and put the following in it:

.. code-block::

   # This file is the ocean core's configuration file. It is specific to the ocean
   # core, and a specific machine. Each machine will configure this file
   # differently, but it can be used to point on version of the testing
   # infrastructure at a different version of the model.


   # The namelists section defines paths to template namelists that will be used
   # to generate specific namelists. Typically these will point to the forward and
   # init namelists in the default_inputs directory after a successful build of
   # the ocean model.
   [namelists]
   forward = /usr/projects/climate/username/mpas/model/ocean/update_isomip_plus_viz/namelist.ocean.forward
   init = /usr/projects/climate/username/mpas/model/ocean/update_isomip_plus_viz/namelist.ocean.init


   # The streams section defines paths to template streams files that will be used
   # to generate specific streams files. Typically these will point to the forward and
   # init streams files in the default_inputs directory after a successful build of
   # the ocean model.
   [streams]
   forward = /usr/projects/climate/username/mpas/model/ocean/update_isomip_plus_viz/streams.ocean.forward
   init = /usr/projects/climate/username/mpas/model/ocean/update_isomip_plus_viz/streams.ocean.init


   # The executables section defines paths to required executables. These
   # executables are provided for use by specific test cases.
   # Full paths should be provided in order to access the executables from
   # anywhere on the machine.
   [executables]
   model = /usr/projects/climate/username/mpas/model/ocean/update_isomip_plus_viz/ocean_model


   # The paths section describes paths that are used within the ocean core test
   # cases.
   [paths]

   # The mesh_database and the initial_condition_database are locations where
   # meshes / initial conditions might be found on a specific machine. They can be
   # the same directory, or different directory. Additionally, if they are empty
   # some test cases might download data into them, which will then be reused if
   # the test case is run again later.
   mpas_model = /usr/projects/climate/username/mpas/model/ocean/update_isomip_plus_viz
   mesh_database = /usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/mesh_database
   initial_condition_database = /usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/initial_condition_database
   bathymetry_database = /usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/bathymetry_database

In theory, you can point to default namelists, streams files and executables for other branches than
the one you're currently on but that's very rarely (if ever) going to be useful to you so you'll
just have to bear with all these redundant references to

.. code-block::

   /usr/projects/climate/username/mpas/model/ocean/update_isomip_plus_viz

If you want to set up a worktree for a different branch, the ``config.ocean`` looks the same except
that you would need to replace the above path with the one for your new worktree.

List the available test cases:

.. code-block:: bash

   ./list_testcases.py

At present, there are 107 of them!  Let's look at only the ISOMIP+ ones (component: ``ocean``\ , case: ``isomip_plus``\ ):

.. code-block:: bash

   ./list_testcases.py -o ocean -c isomip_plus

There are 2 resolutions (2 km and 5 km) and 3 test cases at each resolution (Ocean0, 1 and 2).  For now, we're
going to focus on Ocean0, which has boundary conditions and ocean properties consistent with a (very) warm
continental shelf.  This one spins up to a quasi-steady state in about 2 years (compared to several decades
for the other 2, which are purposefully designed as transient experiments) so it's a good starting point.
We'll use the 2 km version because the domain is only 80 km wide, so 5 km is really quite coarse.  Plus, this
is the "standard" resolution for ISOMIP+.

Set up the test case as follows:

.. code-block:: bash

   ./setup_testcase.py -o ocean -c isomip_plus -r 2km -t Ocean0 -f config.ocean -m runtime_definitions/srun.xml --work_dir /lustre/scratch4/turquoise/username/isomip_plus_Ocean0

6. Running the test case
------------------------

We'll do a short test run (1 month) to make sure everything is working, rathere than jumping into a 2-year simulation.

.. code-block:: bash

   cd /lustre/scratch4/turquoise/username/isomip_plus_Ocean0/ocean/isomip_plus/2km/Ocean0/
   salloc --nodes=1 --time=0:20:00 --account=e3sm

   module purge
   source /usr/projects/climate/SHARED_CLIMATE/anaconda_envs/base/etc/profile.d/conda.sh
   conda activate compass_py3.7
   module use /usr/projects/climate/SHARED_CLIMATE/modulefiles/all/
   module load intel/17.0.1 openmpi/1.10.5 netcdf/4.4.1 parallel-netcdf/1.5.0 pio/1.7.2

   ./run_test.py

If you don't have access to the ``e3sm`` account, ask Steve or Mark for help to get acces.  Somewhere on the
HPC website, there is a way to ask for access, but they may just be able to add you directly.

7. Running a full 2-year Ocean0 simulation
------------------------------------------

For this one, you should use a job script.

.. code-block::

   cd /lustre/scratch4/turquoise/username/isomip_plus_Ocean0/ocean/isomip_plus/2km/Ocean0/forward
   vim job_script.bash

Put this in the job script:

.. code-block::

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

   module purge
   source /usr/projects/climate/SHARED_CLIMATE/anaconda_envs/base/etc/profile.d/conda.sh
   conda activate compass_py3.7
   module use /usr/projects/climate/SHARED_CLIMATE/modulefiles/all/
   module load intel/17.0.1 openmpi/1.10.5 netcdf/4.4.1 parallel-netcdf/1.5.0 pio/1.7.2

   months_per_job=24
   end_date="0003-01-01_00:00:00"

   for month in `seq 0 $months_per_job`
   do
       ./check_progress.py -f namelist.ocean -e $end_date
       ./run.py
       ./setup_restart.py -f namelist.ocean
   done

Submit the job:

.. code-block::

   sbatch job_script.bash

Once it's running, monitor the progress with:

.. code-block::

   tail log.ocean.0000.out

This writes a message for each time step (if all is going well).

The simulation runs one month at a time and then does some adjustment in a python script to make sure sea level doesn't
get out of control (there's a lot of melting going on so we have to have a compensating level of "evaporation" at the
domain boundary).  It also will check to see if we've already reached year 2 and won't run again if so.

Some basic output is available in:

.. code-block::

   analysis_members/globalStats.0001-01-01_00.00.00.nc

To see the mean melt flux and how time is progressing there, do:

.. code-block::

   ncdump -v xtime,landIceFreshwaterFluxAvg analysis_members/globalStats.0001-01-01_00.00.00.nc | tail -n 50

Keep in mind that the units are ``kg m^{-2} s^{-1}``\ , not m/yr, so not the most intuitive output.  There are
some pretty outdated viz scripts in the ``viz`` directory linked there, but these might at least provide some
starting guidelines for how to do python viz.  You can also look at output in paraview.  I'll clean things
up and add instructions for viz in the near future as I have time.

8. Visualization
----------------

8.1 Running the default viz
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Viz should be light enough weight that you can run it on the login node but you could get an interactive job if you prefer.
It produces images, rather than anything interactive, so no need for x-windows or anything like that.

There should be a link to ``viz`` in the ``forward`` output directory.  This is a link to a python package (you can tell because
it contains a ``__init__.py`` (which is empty) and a ``__main__.py``\ , which is the main script for visualization.  To start
with, we'll run the default viz.  If you don't already have the compass conda environment loaded, do:

.. code-block:: bash

   source /usr/projects/climate/SHARED_CLIMATE/anaconda_envs/base/etc/profile.d/conda.sh
   conda activate compass_py3.7

Then, run:

.. code-block:: bash

   python -m viz

This will run the ``main()`` function in ``__main__.py``.  You could optionally set the input directory and the experiment
number but the defaults are the current directory and ``Ocean0``\ , respectively, so there's no need in this case.
This will take maybe 10 or 15 minutes (most of it on the overturning streamfunction).  You should see something like:

.. code-block::

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

The more interesting results should be a series of movies in ``movies`` and 4 time series plots in ``plots``
(mean melt rate, total melt flux, mean thermal driving and mean friction velocity) and the same plots in
``timeSeriesBelow300m``\ , but this time averaged only over the deepest part of the ice shelf (where much of the action is).

You'll likely need to scp or rsync them to your laptop to view them.  Let me know if it's not clear what these are.

8.2 Doing your own viz
^^^^^^^^^^^^^^^^^^^^^^

A starting point for doing your own viz is to make a local copy of ``__main__.py`` to edit:

.. code-block:: bash

   cp viz/__main__.py myviz.py
   vim myviz.py

You could, for example, take out the slow streamfunction stuff if you don't need that (it was added because I required it
as standard output in MISOMIP).

The script imports the following

.. code-block:: py

   from viz.streamfunction import compute_barotropic_streamfunction, \
       compute_overturning_streamfunction

These are functions for computing the stream functions and writing them to NetCDF files.

.. code-block:: py

   from viz.plot import MoviePlotter, TimeSeriesPlotter

These can be used to create "plotter" object that can then produce either time-series plots or a series of image for making movies.

.. code-block:: py

   from viz.misomip import compute_misomip_interp_coeffs, interp_misomip

These are used to write out MISOMIP standard output on a regular grid.

You can look at ``streamfunction.py``\ , ``plot.py`` and ``misomip.py`` to learn a bit more about what these do.  There's a bit
of commenting, particularly for the "public" functions that don't start with an underscore.

Maybe simplify it down to eliminate the streamfunction and MISOMIP stuff, and don't worry about the plots averaged over
the deeper part of the ice draft (none of this is probably all that relevant to you):

.. code-block:: py

   #!/usr/bin/env python

   import xarray
   import argparse

   from viz.plot import MoviePlotter, TimeSeriesPlotter

   def main():
       parser = argparse.ArgumentParser(
           description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
       parser.add_argument("-f", "--folder", dest="folder",
                           help="Folder for plots", default='.')
       parser.add_argument("-e", "--expt", dest="expt",
                           help="Experiment number (0, 1 or 2)", default=0)
       args = parser.parse_args()

       folder = args.folder
       expt = args.expt

       dsMesh = xarray.open_dataset('{}/init.nc'.format(folder))

       ds = xarray.open_mfdataset('{}/timeSeriesStatsMonthly*.nc'.format(folder),
                                  concat_dim='Time')

       tsPlotter = TimeSeriesPlotter(inFolder=folder,
                                     outFolder='{}/plots'.format(folder),
                                     expt=expt)
       tsPlotter.plot_melt_time_series()

       mPlotter = MoviePlotter(inFolder=folder,
                              outFolder='{}/plots'.format(folder),
                              expt=expt)

       mPlotter.plot_melt_rates()
       mPlotter.plot_ice_shelf_boundary_variables()
       mPlotter.plot_temperature()
       mPlotter.plot_salinity()
       mPlotter.plot_potential_density()

       mPlotter.images_to_movies(outFolder='{}/movies'.format(folder),
                                 framesPerSecond=30, extension='mp4')

   if __name__ == '__main__':
       main()

I've set things up to plot some of the more common fields by default.  The following plot either time series or
movies of some common fields related to the ice-ocean interface -- melt rates, thermal driving, friction velocity,
etc.

.. code-block:: py

       tsPlotter.plot_melt_time_series()
       ...
       mPlotter.plot_melt_rates()
       mPlotter.plot_ice_shelf_boundary_variables()

These functions plot 3D fields at the top of the ocean (either the ice draft or the sea surace), the sea floor
and in a transect through the middle of the domain:

.. code-block:: py

       mPlotter.plot_temperature()
       mPlotter.plot_salinity()
       mPlotter.plot_potential_density()

You could also add your own custom fields as long as they're available in the ``timeSeriesStatsMonthly*.nc`` files.

Here are a couple of examples:

.. code-block::

       # plot a time series of SST
       areaCell = tsPlotter.dsMesh.areaCell
       temperature = tsPlotter.ds.timeMonthly_avg_activeTracers_temperature
       sst = temperature.isel(nVertLevels=0)
       meanSST = (sst*areaCell).sum(dim='nCells')/areaCell.sum(dim='nCells')

       tsPlotter.plot_time_series(meanSST, 'mean sea-surface temperature',
                                  prefix='meanSST', units='deg C')
   ...
       # plot the x and y components of velocity at top, bottom and transect
       da = mPlotter.ds.timeMonthly_avg_velocityX
       mPlotter.plot_3d_field_top_bot_section(
           da, nameInTitle='x velocity', prefix='Vx', units='m/s',
           vmin=-0.2, vmax=0.2)
       da = mPlotter.ds.timeMonthly_avg_velocityY
       mPlotter.plot_3d_field_top_bot_section(
           da, nameInTitle='y velocity', prefix='Vy', units='m/s',
           vmin=-0.2, vmax=0.2)

Make sure any new plots with the movie plotter happen before movies are made (\ ``mPlotter.images_to_movies()``\ )
so they get included in the movies.

The data sets (\ ``ds``\ ) and data arrays (\ ``da``\ ) come from ``xarray``\ , which is a really handy package for working with
NetCDF-style files in memory in python.  It's a lot smarter about named dimensions than ``numpy`` and a lot more easy
to manipulate than default python ``NetCDF4`` data sets.  But there's a bit of a learning curve involving a lot of Googling
the documentation and StackOverflow.

Hopefully that's a start...

