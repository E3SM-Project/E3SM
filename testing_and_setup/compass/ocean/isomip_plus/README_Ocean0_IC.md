# Instructions for setting up and running ISOMIP+ Ocean0 on LANL IC

## 1. SSH tricks

A couple of tricks for your laptop if you’re not already using them:
Save your SSH connections:
```bash
vim ~/.ssh/config
```
Add the following:
```
Host *
    ControlMaster auto
    ControlPath ~/.ssh/connections/%r@%h:%p
    ServerAliveInterval 300

Host wtrw
    Hostname wtrw.lanl.gov
    User cbegeman
```
```bash
mkdir ~/.ssh/connections
```
Alias connections to LANL HPC machines:
```bash
vim ~/.bashrc
```
Add:
```bash
alias gr="ssh -t wtrw ssh gr-fe1"
alias ba="ssh -t wtrw ssh ba-fe2"
alias ko="ssh -t wtrw ssh ko-fe"
```
## 2. Making sure git is set up nicely

### 2.1 Storing your LANL IC SSH key on GitHub

It's useful to set up GitHub to know your public SSH key from LANL IC if you
haven't already done this.  It means you don’t have to type your password for
GitHub each time you git fetch, git push, etc.

I believe this is the right link for
[more details](https://help.github.com/en/articles/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)
If you haven't done this already and this gives you trouble, let me know and we
can work through it together.

### 2.2 git settings

On IC or on your laptop, make sure you’ve got these settings defined:
```bash
git config --global user.name "First Last"
git config --global user.email user@domain.com
git config --global core.editor vim
git config --global push.default nothing
git config --global color.ui true
git config --global core.autocrlf input
git config --global core.whitespace trailing-space
git config --global alias.logg "log --graph --oneline --decorate"
```
I use `git logg` all the time so this last alias is particularly important.

### 2.3 git tab completion

Download [git-completion.bash](https://raw.githubusercontent.com/git/git/master/contrib/completion/git-completion.bash)
```bash
cd ~
wget https://raw.githubusercontent.com/git/git/master/contrib/completion/git-completion.bash
```
Add this to your .bashrc
```
module load git
source git-completion.bash
```

## 3. Forking and Cloning MPAS-Model

* Go to: [https://github.com/MPAS-Dev/MPAS-Model](https://github.com/MPAS-Dev/MPAS-Model)
* Make your own fork by clicking “Fork” at the top right:
* Go to your new fork (e.g. [https://github.com/cbegeman/MPAS-Model](https://github.com/cbegeman/MPAS-Model))
* Whenever you ever need to know the link to clone your fork
  * Click on “Clone or download”
  * If it says “Clone with HTTPS”, click Use SSH (either works but SSH will use
    the SSH keys you’ve set up above and you never have to type my Git
    password.)
  * Copy the link with the clipboard icon

In a terminal window, log in to a LANL machine (I use Grizzly from here on
  except where stated): 
```bash
ssh -t wtrw ssh gr-fe1`
```
Make a directory for the code, e.g.:
```bash
mkdir /usr/projects/climate/cbegeman
cd /usr/projects/climate/cbegeman
mkdir -p mpas/model
cd mpas/model/
```
Clone the repo:
```bash
git clone git@github.com:cbegeman/MPAS-Model.git repo`
cd repo
```
Rename your remote so it’s easier to not confuse it with other forks:
```bash
git remote rename origin cbegeman/MPAS-Model
```
Add the main repo:
```bash
git remote add MPAS-Dev/MPAS-Model git@github.com:MPAS-Dev/MPAS-Model.git
```
Add my fork (you can add other people’s forks in the same way):
```bash
git remote add xylar/MPAS-Model git@github.com:xylar/MPAS-Model.git
```
Get the latest version of all the remotes (pruning anything that has been
deleted):
```bash
git fetch --all -p
```
Let's store some settings you'll need to load every time you build MPAS.  The following
are only appropriate for Grizzly and we'll need a similar file with settings for
Badge and any other machines we might use in the future.
```bash
vim ../setup_gr.bash
```
In this file, put:
```bash
echo "Setting up grizzly intel environment for building and running MPAS"
module purge
module load git
source /usr/projects/climate/SHARED_CLIMATE/anaconda_envs/base/etc/profile.d/conda.sh
conda activate compass_py3.7
module use /usr/projects/climate/SHARED_CLIMATE/modulefiles/all/
module load intel/17.0.1 openmpi/1.10.5 netcdf/4.4.1 parallel-netcdf/1.5.0 pio/1.7.2
export CORE=ocean
```

## 4. Cloning JIGSAW

There's one more bit of code that some test cases need (not Ocean0, but it's best to get
in in case you want to work with other test cases) so we'll clone that, too.
```bash
cd /usr/projects/climate/cbegeman/mpas
git clone git@github.com:dengwirda/jigsaw-geo-matlab.git
```

## 5. Checking out an MPAS branch and building the model

Add a “worktree”, a copy of the repo that we can point to a different branch.
We will work with my branch `ocean/develop`, where the latest MPAS-Ocean development 
is taking place.  In general, `ocean/develop` is the place to start, since the `master` 
branch is updated only rarely when we make releases:
```bash
cd /usr/projects/climate/cbegeman/mpas/model/reop
git worktree add ../ocean/develop -b ocean/develop
cd ../ocean/develop
```
Take a look at which branch were on:
```
git logg
```
We don't start off on `MPAS-Dev/MPAS-Model/ocean/develop` (even though
the name of the local branch might trick you into thinking you're there), so we need 
to do a hard reset to put us there:
```bash
git reset --hard MPAS-Dev/MPAS-Model/ocean/develop
git logg
```
Now source the file with modules and settings for building MPAS on grizzly:
```bash
source ../../setup_gr.bash
```
If all goes well, you should see `comapss_py3.7` as part of your command prompt and you should be read to build MPAS.
```bash
make ifort
```
There are currently 2 files involved in the build that aren't python 3 compatible.  They are downloaded by the build
process, so they can't easily be fixed before starting the build.  The way to handle them is to wait for the build to
fail (twice) and then run some intermediate commands to replace the broken files.  (I'm sorry for this, but it looks
like a proper solution is a little ways off).

After the first build error, which should mention something about `makedep.py` ins a directory with something about
CVMix, run the following:
```bash
wget https://raw.githubusercontent.com/CVMix/CVMix-src/master/src/shared/makedep.py -O src/core_ocean/.cvmix_all/src/shared/makedep.py
make ifort
```
Then, after a second, very similar error but which should mention something about BGC this time, run:
```bash
wget https://raw.githubusercontent.com/xylar/Ocean-BGC/python_3_support/makedep.py -O src/core_ocean/.BGC_all/makedep.py
make ifort
```
Take a coffee break, this will take some time.
...

## 6. Setting up a test case

Okay you're back and refreshed?  Let's set up a test case.
```bash
cd testing_and_setup/compass/
```
COMPASS (COnfiguration of Model for Prediction Across Scales Setups -- yes, a litle tortured) is a set of python 
scripts we use to set up and run our test cases.  To build test cases, you need to tell COMPASS where to find a few
thing on Grizzly.  Open a file `config.ocean` and put the following in it:
```
# This file is the ocean core's configuration file. It is specific to the ocean
# core, and a specific machine. Each machine will configure this file
# differently, but it can be used to point on version of the testing
# infrastructure at a different version of the model.


# The namelists section defines paths to template namelists that will be used
# to generate specific namelists. Typically these will point to the forward and
# init namelists in the default_inputs directory after a successful build of
# the ocean model.
[namelists]
forward = /usr/projects/climate/cbegeman/mpas/model/ocean/develop/namelist.ocean.forward
init = /usr/projects/climate/cbegeman/mpas/model/ocean/develop/namelist.ocean.init


# The streams section defines paths to template streams files that will be used
# to generate specific streams files. Typically these will point to the forward and
# init streams files in the default_inputs directory after a successful build of
# the ocean model.
[streams]
forward = /usr/projects/climate/cbegeman/mpas/model/ocean/develop/streams.ocean.forward
init = /usr/projects/climate/cbegeman/mpas/model/ocean/develop/streams.ocean.init


# The executables section defines paths to required executables. These
# executables are provided for use by specific test cases.
# Full paths should be provided in order to access the executables from
# anywhere on the machine.
[executables]
model = /usr/projects/climate/cbegeman/mpas/model/ocean/develop/ocean_model


# The paths section describes paths that are used within the ocean core test
# cases.
[paths]

# The mesh_database and the initial_condition_database are locations where
# meshes / initial conditions might be found on a specific machine. They can be
# the same directory, or different directory. Additionally, if they are empty
# some test cases might download data into them, which will then be reused if
# the test case is run again later.
mpas_model = /usr/projects/climate/cbegeman/mpas/model/ocean/develop
jigsaw-geo-matlab = /usr/projects/climate/cbegeman/mpas/jigsaw-geo-matlab
mesh_database = /usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/mesh_database
initial_condition_database = /usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/initial_condition_database
geometric_data = /usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/geometric_data_v0.1.1
bathymetry_database = /usr/projects/regionalclimate/COMMON_MPAS/ocean/grids/bathymetry_database
```
In theory, you can point to default namelists, streams files and executables for other branches than
the one you're currently on but that's very rarely (if ever) going to be useful to you so you'll
just have to bear with all these redundant references to 
```
/usr/projects/climate/cbegeman/mpas/model/ocean/develop
```
If you want to set up a worktree for a different branch, the `config.ocean` looks the same except
that you would need to replace the above path with the one for your new worktree.

List the available test cases:
```bash
./list_testcases.py
```
At present, there are 107 of them!  Let's look at only the ISOMIP+ ones (component: `ocean`, case: `isomip_plus`):
```bash
./list_testcases.py -o ocean -c isomip_plus
```
There are 2 resolutions (2 km and 5 km) and 3 test cases at each resolution (Ocean0, 1 and 2).  For now, we're
going to focus on Ocean0, which has boundary conditions and ocean properties consistent with a (very) warm 
continental shelf.  This one spins up to a quasi-steady state in about 2 years (compared to several decades
for the other 2, which are purposefully designed as transient experiments) so it's a good starting point.
We'll use the 2 km version because the domain is only 80 km wide, so 5 km is really quite coarse.  Plus, this
is the "standard" resolution for ISOMIP+.

Set up the test case as follows:
```bash
./setup_testcase.py -o ocean -c isomip_plus -r 2km -t Ocean0 -f config.ocean -m runtime_definitions/srun.xml --work_dir /lustre/scratch4/turquoise/cbegeman/isomip_plus_Ocean0
```

## 7. Running the test case

We'll do a short test run (1 month) to make sure everything is working, rathere than jumping into a 2-year simulation.
```bash
cd /lustre/scratch4/turquoise/cbegeman/isomip_plus_Ocean0/ocean/isomip_plus/2km/Ocean0/
salloc --nodes=1 --time=0:20:00 --account=e3sm

module purge
source /usr/projects/climate/SHARED_CLIMATE/anaconda_envs/base/etc/profile.d/conda.sh
conda activate compass_py3.7
module use /usr/projects/climate/SHARED_CLIMATE/modulefiles/all/
module load intel/17.0.1 openmpi/1.10.5 netcdf/4.4.1 parallel-netcdf/1.5.0 pio/1.7.2

./run_test.py
```
If you don't have access to the `e3sm` account, ask Steve or Mark for help to get acces.  Somewhere on the 
HPC website, there is a way to ask for access, but they may just be able to add you directly.


## 8. Running a full 2-year Ocean0 simulation

For this one, you should use a job script.
```
cd /lustre/scratch4/turquoise/cbegeman/isomip_plus_Ocean0/ocean/isomip_plus/2km/Ocean0/forward
vim job_script.bash
```
Put this in the job script:
```
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
```

Submit the job:
```
sbatch job_script.bash
```

Once it's running, monitor the progress with:
```
tail log.ocean.0000.out
```
This writes a message for each time step (if all is going well).

The simulation runs one month at a time and then does some adjustment in a python script to make sure sea level doesn't
get out of control (there's a lot of melting going on so we have to have a compensating level of "evaporation" at the
domain boundary).  It also will check to see if we've already reached year 2 and won't run again if so.

Some basic output is available in:
```
analysis_members/globalStats.0001-01-01_00.00.00.nc
```
To see the mean melt flux and how time is progressing there, do:
```
ncdump -v xtime,landIceFreshwaterFluxAvg analysis_members/globalStats.0001-01-01_00.00.00.nc | tail -n 50
```
Keep in mind that the units are `kg m^{-2} s^{-1}`, not m/yr, so not the most intuitive output.  There are
some pretty outdated viz scripts in the `viz` directory linked there, but these might at least provide some
starting guidelines for how to do python viz.  You can also look at output in paraview.  I'll clean things
up and add instructions for viz in the near future as I have time.


