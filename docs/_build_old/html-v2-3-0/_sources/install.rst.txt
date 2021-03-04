
Installation
============

The installation procedure depends on what version you'd like to install.

Activate **e3sm_unified** environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have an account on one of the E3SM supported machines (Cori, Compy, Acme1, Anvil, Cooley, Rhea), you
can access ``e3sm_diags`` by activating ``e3sm_unified``, which is a conda environment that pulls together Python
and other E3SM analysis tools such as ``e3sm_diags``, ``mpas-analysis``, ``NCO``, ``cdat`` and ``processflow``.

The paths to ``e3sm_unified`` activation scripts are machine dependent:

**Compy**
    ::

     source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified.sh


**Cori**
    ::

     source /global/cfs/cdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified.sh
    

**Anvil/blues**
    ::

     source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified.sh


**Cooley**
    ::

     source /lus/theta-fs0/projects/ccsm/acme/tools/e3sm-unified/load_latest_e3sm_unified.sh


**acme1**
    ::

     source /usr/local/e3sm_unified/envs/load_latest_e3sm_unified.sh


**Rhea**
    ::

     source /ccs/proj/cli900/sw/rhea/e3sm-unified/load_latest_e3sm_unified.sh
 

Change ``.sh`` to ``.csh`` for csh shells.
Note that ``e3sm_unified``'s development cycle is not in phase with ``e3sm_diags``,
therefore the version of ``e3sm_diags`` included may not be the latest.
To install latest stable releases, refer to following:

.. _conda_environment:

Installation in a conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the E3SM Unified environment doesn't serve your needs, you can alternatively
install the latest version in your own custom conda environment.

First, activate conda or install it if not available. Details vary on the machine.

**Compy**
    ::

     module load anaconda3/2019.03
     source /share/apps/anaconda3/2019.03/etc/profile.d/conda.sh


**Cori**
    ::

     module load python/3.7-anaconda-2019.10
     source /global/common/cori_cle7/software/python/3.7-anaconda-2019.10/etc/profile.d/conda.sh
    
    
**Others**
    ::

If the system doesn't come with conda pre-installed, install it following the
`official instructions <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_. 
Once installed, it is recommended to upgrade to the latest version:

   ::

       conda update conda

If you are working on a machine/network that intercepts SSL communications (such as acme1), you will get
an SSL error unless you disable the SSL verification:

   ::

       conda config --set ssl_verify false
       binstar config --set ssl_verify False


Once conda is properly working, you can either install the latest stable release or create a
development environment.

.. _install_latest:

Latest stable release
---------------------

Make sure conda is loaded or installed (see :ref:`above <conda_environment>`).

1. Get the yml file to create an environment. For Linux machines:

   ::

       wget https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/conda/e3sm_diags_env.yml

   For macOS:

   ::

       curl https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/conda/e3sm_diags_env_osx.yml

2. Change ``prefix`` in that file to be your conda prefix. Typically, this will be ``~/miniconda3/envs/e3sm_diags_env``.


3. Remove any cached conda packages. This will ensure that you always get the latest packages.

   ::

       conda clean --all

4. Use conda to create a new environment with ``e3sm_diags`` installed.

   Tip: You can change the name of the environment to anything you'd like using ``-n <my_env_name>``.

   ::

       conda env create -n e3sm_diags_env -f e3sm_diags_env.yml
       conda activate e3sm_diags_env


.. _dev-env:

Development environment
-----------------------

Unlike the latest stable release (i.e., the user environment), the development environment does not include E3SM Diags,
as the developer will ``pip install`` their changes to E3SM Diags (see step 6 below).

Make sure conda is loaded or installed (see :ref:`above <conda_environment>`).

1. Get the developmental yml file to create an environment. For Linux machines:

   ::

       wget https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/conda/e3sm_diags_env_dev.yml

   For macOS:

   ::

       curl https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/conda/e3sm_diags_env_dev_osx.yml

2. Change ``prefix`` in that file to be your conda prefix. Typically, this will be ``~/miniconda3/envs/e3sm_diags_env_dev``.


3. Remove any cached conda packages. This will ensure that you always get the latest packages.

   ::

       conda clean --all

4. Use conda to create a new environment. ``e3sm_diags`` **is not included in this environment.**

   ::

       conda env create -n e3sm_diags_env_dev -f e3sm_diags_env_dev.yml
       conda activate e3sm_diags_env_dev

5. Get the latest code from master

   ::

       git clone https://github.com/E3SM-Project/e3sm_diags.git


   or if you already have a clone of the repo, pull the latest code from master

   ::

       git pull origin master


   or checkout a new branch from master.

   ::

       git fetch origin master
       git checkout -b <branch-name> origin/master


6. Make any changes to E3SM Diags you want, then install with

   ::

       pip install .

7. Run a quick test which generates one of each plot type.

   ::

       cd tests/system
       python all_sets.py -d all_sets.cfg

8. Remember to view the generated html located here: ``all_sets/viewer/index.html``. These plots can be moved to the web
   for viewing by moving the generated directory ``all_sets`` to the ``html_path``. Each machine has a different
   ``html_path`` -- see :doc:`quick guide <quickguides/quick-guide-general>`. Files at the ``html_path`` can be viewed
   at ``web_address``. If you're not seeing the files there, you may need to change the permissions with ``chmod -R``
   (e.g., on Cori, ``chmod -R 755 /global/cfs/cdirs/e3sm/www/<username>/all_sets``).
