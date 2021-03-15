Installation
============

The installation procedure depends on what version you'd like to install.

Activate **e3sm_unified** environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have an account on one of the E3SM supported machines (NERSC, Compy, Acme1, LCRC, Cooley, Rhea), you
can access ``e3sm_diags`` by activating ``e3sm_unified``, which is a conda environment that pulls together Python
and other E3SM analysis tools such as ``e3sm_diags``, ``mpas-analysis``, ``NCO``, ``cdat`` and ``processflow``.

The paths to ``e3sm_unified`` activation scripts are machine dependent:

**Compy**
    ::

     source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified.sh


**NERSC**
    ::

     source /global/cfs/cdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified.sh


**LCRC**
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

Installation in a Conda Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the E3SM Unified environment doesn't serve your needs, you can alternatively
install the latest version in your own custom conda environment.

First, activate conda or install it if not available. Details vary on the machine.

Compy
~~~~~
    ::

     module load anaconda3/2019.03
     source /share/apps/anaconda3/2019.03/etc/profile.d/conda.sh


NERSC
~~~~~
    ::

     module load python/3.7-anaconda-2019.10
     source /global/common/cori_cle7/software/python/3.7-anaconda-2019.10/etc/profile.d/conda.sh

.. _conda_environment_others:

Others/Local
~~~~~~~~~~~~

If the system doesn't come with conda pre-installed, follow these instructions:

1. Download Conda

    Linux
        ::

            wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

    MacOS
        ::

            wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

2. Install Conda

    Linux
        ::

            bash ./Miniconda3-latest-Linux-x86_64.sh

    MacOS
        ::

            bash ./Miniconda3-latest-MacOSX-x86_64.sh

    - ``Do you wish the installer to initialize Miniconda3 by running conda init? [yes|no] yes``

3. If you are working on a machine/network that intercepts SSL communications (such as acme1), you will get
an SSL error unless you disable the SSL verification:

    ::

        conda config --set ssl_verify false
        binstar config --set ssl_verify False

4. Once conda is properly working, you can install the **(a) Latest Stable Release** or create a **(b) Development Environment**.

.. _install_latest:

(a) Latest Stable Release
-------------------------

1. Follow :ref:`"Others/Local" <conda_environment_others>` section for installing Conda.

2. Get the yml file to create an environment.

    ::

        wget https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/conda/e3sm_diags_env.yml


3. Change ``prefix`` in that file to be your conda prefix. Typically, this will be ``~/miniconda3/envs/e3sm_diags_env``.

4. Remove any cached conda packages. This will ensure that you always get the latest packages

    ::

        conda clean --all

5. Use conda to create a new environment with E3SM Diags (``e3sm_diags``) included.

    - Tip: Add the flag ``-n <name_of_env>`` to customize the name of the environment

    ::

        conda env create -f e3sm_diags_env.yml
        conda activate e3sm_diags_env

.. _dev-env:

(b) Development Environment
---------------------------

Unlike the latest stable release (i.e., the user environment), the development environment does not include E3SM Diags (``e3sm-diags``).
Instead, the developer will ``pip install .`` to build ``e3sm-diags`` with changes (see step 6 below).

Furthermore, the dev environment includes QA tools such as code formatters, linters, and ``pre-commit``.
**NOTE**: These QA tools are enforced using ``pre-commit`` checks in the CI/CD build, so you must use the dev environment for all contributions

1. Follow :ref:`"Others/Local" <conda_environment_others>` section for installing conda.

2. Clone your fork and keep in sync with upstream ``master``

    ::

        git clone https://github.com/<your-github-username>/e3sm_diags.git


   or if you already have a clone of the fork, rebase your fork with upstream to keep in sync

    ::

        # Add the remote, call it "upstream":
        git remote add upstream https://github.com/E3SM-Project/e3sm_diags.git

        # Fetch all the branches of that remote into remote-tracking branches
        git fetch upstream

        # Make sure that you're on your master branch:
        git checkout master

        # Rewrite your master branch so that any of your commits that
        # aren't already in upstream/master are replayed on top of the
        # other branch:
        git rebase upstream/master
        git push -f origin master


   Checkout a new branch from ``master``.

    ::

        git checkout -b <branch-name> master

3. Remove any cached conda packages. This will ensure that you always get the latest packages.

    ::

        conda clean --all

4. Enter the fork directory.

    ::

        cd e3sm_diags

5. Use conda to create a new dev environment (``e3sm_diags`` **is not included in this environment**).

    - Tip: Add the flag ``-n <name_of_env>`` to customize the name of the environment

    ::

        conda env create -f conda/e3sm_diags_env_dev.yml
        conda activate e3sm_diags_env_dev

6. Install ``pre-commit``.

    ::

        pre-commit install

7. Make the desired changes to E3SM Diags, then rebuild and install with:

    ::

        pip install .

8. Run a quick test which generates one of each plot type.

    ::

        cd tests/system
        python all_sets.py -d all_sets.cfg

9. Remember to view the generated html located here: ``all_sets/viewer/index.html``. These plots can be moved to the web
   for viewing by moving the generated directory ``all_sets`` to the ``html_path``. Each machine has a different
   ``html_path`` -- see :doc:`quick guide <quickguides/quick-guide-general>`. Files at the ``html_path`` can be viewed
   at ``web_address``. If you're not seeing the files there, you may need to change the permissions with ``chmod -R``
   (e.g., on NERSC, ``chmod -R 755 /global/cfs/cdirs/e3sm/www/<username>/all_sets``).

10. Commit changes and make sure ``pre-commit`` checks pass
    ::

        git commit -m "..."

    .. figure:: pre-commit-passing.png
       :alt: pre-commit Output

       ``pre-commit`` Output
