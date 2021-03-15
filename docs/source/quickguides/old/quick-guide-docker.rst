
Quick guide for using Docker (v1)
==================================

Docker allows software packages, called containers, to be ran.
E3SM Diagnostics has a Docker container that makes running the
software easier than traditional methods.
Any machine that has Docker installed can use the
instructions below to run E3SM Diagnostics.

Obtaining the container image
-----------------------------

1. Make sure that Docker is installed and running on your machine.
To check this, run the command below and you should see some output.

    ::

        docker info

2. View all of the ``e3sm_diags`` images.

    ::

        docker images | grep e3sm_diags

If the version you want to use is already available, then please continue to step 4.

Otherwise, you'll need to download the image you want, shown in step 3.

3. If the specific version you want or the ``latest`` image is **not** shown,
download it. You can view all of the images available on the
`e3sm_diags Docker Hub <https://hub.docker.com/r/e3sm/e3sm_diags/tags/>`_.

Below, we are getting the image with the latest tag:

    ::

        docker pull e3sm/e3sm_diags:latest

Running a quick test
--------------------

Since you probably don't have much sample data on your machine,
below are the steps to run a quick test. If you do have data,
please go to the "Running more diagnostics" section below.

4. Clone the ``e3sm_diags`` repo, and go to the test directory.

    ::

        git clone https://github.com/E3SM-Project/e3sm_diags.git
        cd e3sm_diags/tests/system/

5. ``wget`` or ``curl`` the script to run the container.

    ::

        wget https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/acme_diags/container/e3sm_diags_container.py

        # Or use this:
        curl -O https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/acme_diags/container/e3sm_diags_container.py

6. Run your diagnostics and examine the sample output.

    ::

        python e3sm_diags_container.py --docker -p all_sets.py -d all_sets.cfg

**Tip:** You can select the version of the container you want to run with the ``--container_version`` argument.
If this argument isn't defined, it defaults to the ``latest`` container.

    ::

        python e3sm_diags_container.py --docker --container_version v1.6.0 -p myparams.py


7. Open the html below to view the results.

    ::

        all_sets_results/viewer/index.html

Running more diagnostics
------------------------

To run other, more interesting diagnostics, you must download the data from one of our supported
machines (ALCF Cooley, NERSC Cori and others). 
For more information on the format of the input data, please see the
:doc:`input data requirements <../../input-data-requirement>`.

Below are the paths to the **observational data**:

    * Climatology data (6GB):
        * NERSC: ``/global/project/projectdirs/acme/acme_diags/obs_for_e3sm_diags/climatology/``
        * ALCF: ``/lus/theta-fs0/projects/ClimateEnergy_3/e3sm_diags/obs_for_e3sm_diags/climatology/``
    * Time-series data (145GB):
        * NERSC: ``/global/project/projectdirs/acme/acme_diags/obs_for_e3sm_diags/time-series/``
        * ALCF: ``/lus/theta-fs0/projects/ClimateEnergy_3/e3sm_diags/obs_for_e3sm_diags/time-series/``

We also have **sample model data** as well. You can use your own model data as well,
either climatology or time-series files created via ``nco``.
Again, if you want to use your own data, please see the
:doc:`input data requirements <../../input-data-requirement>`.

    * Climatology data (42GB):
        * We have data from three models in these directories, from about 11.5GB to 19GB.
        * NERSC: ``/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/climatology/``
        * ALCF: ``/lus/theta-fs0/projects/ClimateEnergy_3/e3sm_diags/test_model_data_for_e3sm_diags/climatology/``
    * Time-series data (107GB):
        * We have E3SM v1 data (94GB) and CESM1-CAM5 CMIP data (14GB).
        * NERSC: ``/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/time-series/``
        * ALCF: ``/lus/theta-fs0/projects/ClimateEnergy_3/e3sm_diags/test_model_data_for_e3sm_diags/time-series/``

Once the data is downloaded you can follow one of the
:doc:`many examples that we have <../../examples/index>`.

Some points to remember:

    * You must change the ``reference_data_path`` and ``test_data_path`` accordingly.

    * Every instance of ``e3sm_diags`` should be ``python e3sm_diags_container.py --docker``.
    
        * Ex: Use ``python e3sm_diags_container.py --docker -p myparams.py`` instead of ``e3sm_diags -p myparams.py``.
