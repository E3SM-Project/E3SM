.. _porting:

*********************************************
Porting and validating CIME on a new platform
*********************************************

One of the first steps for many users is getting CIME-based models running on their local machine.
This section describes that process.

===========================
Required libraries/packages
===========================

The machine needs to have:

- a functioning MPI environment (unless you plan to run on a single core with the CIME mpi-serial library).
- build tools gmake and cmake,
- a netcdf library version 4.3 or newer built with the same compiler you will use for CIME.

A pnetcdf library is optional.

If you are using MPI, make sure you can run a basic MPI parallel program on your machine before you attempt a CIME port. You can use this :ref:`MPI example <mpi-example>` to check.

=================
Steps for porting
=================

Porting CIME involves several steps in which you create, at a minimum, a **config_machines.xml** file in your **$HOME/.cime** directory.
In addition, if you have a batch system, you will also need to add a **config_batch.xml** file to your **$HOME/.cime** directory.

All files in **$HOME/.cime/** are appended to the xml objects that are read into memory from the **$CIME/config/[model]**, where **[model]** is either ``acme`` or ``cesm``.

Follow these steps:

#. Create a **$HOME/.cime** directory

#. Create a **config** file in that directory. It should contain the following two lines:
   ::

      [main]
      CIME_MODEL=cesm

   or

   ::

      [main]
      CIME_MODEL=acme

#. Create a **config_machines.xml** file in the same directory.

   This file contains all the information you need to set in order to configure a new machine to be CIME-compliant.

   Use one of the templates described here to create the file.

   - If you are on a MAC and have the libraries installed with either MacPorts or HomeBrew, copy
     **$CIME/config/cesm/machines/userdefined_laptop_template/config_machines.xml** to
     **$HOME/.cime/config_machines.xml**.

   - Otherwise, copy **$CIME/config/xml_schemas/config_machines_template.xml** to
     **$HOME/.cime/config_machines.xml**.

   Fill in the contents of **$HOME/.cime/config_machines.xml** that are specific to your machine. For more details see :ref:`customize the config_machines.xml file <customizing-machine-file>`.

   Check to ensure that your **config_machines.xml** file conforms to the CIME schema definition by doing the following:
   ::

      xmllint --noout --schema $CIME/config/xml_schemas/config_machines.xsd $HOME/.cime/config_machines.xml

#. If you have compiler settings that are specific to your machine, create a **$HOME/.cime/config_compilers.xml** file.

   The default compiler settings are set in **$CIME/config/[model]/machines/config_compilers.xml**, where **[model]** can be either ``acme`` or ``cesm``.

   There is no template for **config_compilers.xml**.

#.  If you have a batch system, create a **$HOME/.cime/config_batch.xml** file.

   Out-of-the-box batch settings are set in **$CIME/config/[model]/machines/config_batch.xml**, where **[model]** can be either ``acme`` or ``cesm``.

#. Once you have defined a basic configuration for your machine in your **$HOME/.cime** xml files, run **scripts_regression_test.py** interactively from the **$CIME/scripts/tests** directory.
   This performs a number of basic unit tests starting from the simplest and working toward more complicated ones.

After running those steps correctly, you are ready to try a case at your target compset and resolution.
   Once you have successfully created the required xml files in your .cime directory and are satisfied with the results you can merge them into the default files in the **config/$CIME_MODEL/machines** directory.
   If you would like to make this machine definition available generally you may then issue a pull request to add your changes to the git repository.

.. _customizing-machine-file:

===========================================
config_machines.xml - machine specific file
===========================================

The machine-specific file is defined in the model-specific :ref:`config_machines.xml <defining-machines>`.

To make a machine CIME-compatible, add the appropriate entries for the machine in **config_machines.xml**.

Each ``<machine>`` tag requires the following input:

- ``DESC``: a text description of the machine
- ``NODENAME_REGEX``: a regular expression used to identify the machine. It must work on compute nodes as well as login nodes. Use the ``machine`` option for **create_test** or **create_newcase** if this flag is not available.
- ``OS``: the machine's operating system
- ``PROXY``: optional http proxy for access to the internet
- ``COMPILERS``: compilers supported on the machine, in comma-separated list, default first
- ``MPILIBS``: mpilibs supported on the machine, in comma-separated list, default first
- ``PROJECT``: a project or account number used for batch jobs; can be overridden in environment or in **$HOME/.cime/config**
- ``SAVE_TIMING_DIR``: (ACME only) target directory for writing timing output
- ``CIME_OUTPUT_ROOT``: Base directory for case output; the **bld** and **run** directories are written below here
- ``DIN_LOC_ROOT``: location of the input data directory
- ``DIN_LOC_ROOT_CLMFORC``: optional input location for clm forcing data
- ``DOUT_S_ROOT``: root directory of short-term archive files
- ``DOUT_L_MSROOT``: root directory on mass store system for long-term archive files
- ``BASELINE_ROOT``: root directory for system test baseline files
- ``CCSM_CPRNC``: location of the cprnc tool, which compares model output in testing
- ``GMAKE``: gnu-compatible make tool; default is "gmake"
- ``GMAKE_J``: optional number of threads to pass to the gmake flag
- ``TESTS``: (ACME only) list of tests to run on the machine
- ``BATCH_SYSTEM``: batch system used on this machine (none is okay)
- ``SUPPORTED_BY``: contact information for support for this system
- ``MAX_TASKS_PER_NODE``: maximum number of threads/tasks per shared memory node on the machine
- ``PES_PER_NODE``: number of physical PES per shared node on the machine. In practice the MPI tasks per node will not exceed this value.
- ``PROJECT_REQUIRED``: Does this machine require a project to be specified to the batch system?
- ``mpirun``: The mpi exec to start a job on this machine.
  This is itself an element that has sub-elements that must be filled:

  * Must have a required ``<executable>`` element
  * May have optional attributes of ``compiler``, ``mpilib`` and/or ``threaded``
  * May have an optional ``<arguments>`` element which in turn contains one or more ``<arg>`` elements.
    These specify the arguments to the mpi executable and are dependent on your mpi library implementation.


- ``module_system``: How and what modules to load on this system. Module systems allow you to easily load multiple compiler environments on a machine. CIME provides support for two types of module tools: `module <http://www.tacc.utexas.edu/tacc-projects/mclay/lmod>`_ and `soft  <http://www.mcs.anl.gov/hs/software/systems/softenv/softenv-intro.html>`_. If neither of these is available on your machine, simply set ``<module_system type="none"\>``.

- ``environment_variables``: environment_variables to set on the system
   This contains sub-elements ``<env>`` with the ``name`` attribute specifying the environment variable name, and the element value specifying the corresponding environment variable value. If the element value is not set, the corresponding environment variable will be unset in your shell.

   For example, the following sets the environment variable ``OMP_STACKSIZE`` to 256M:
   ::

      <env name="OMP_STACKSIZE">256M</env>

   The following unsets this environment variable in the shell:
   ::

      <env name="OMP_STACKSIZE"></env>

   .. note:: These changes are **ONLY** activated for the CIME build and run environment, **BUT NOT** for your login shell. To activate them for your login shell, source either **$CASEROOT/.env_mach_specific.sh** or **$CASEROOT/.env_mach_specific.csh**, depending on your shell.

.. _customizing-compiler-file:

=================================================
config_compilers.xml - compiler paths and options
=================================================

The **config_compilers.xml** file defines compiler flags for building CIME (and also CESM and ACME prognostic CIME-driven components).

#. General compiler flags (e.g., for the gnu compiler) that are machine- and componen-independent are listed first.

#. Compiler flags specific to a particular operating system are listed next.

#. Compiler flags that are specific to particular machines are listed next.

#. Compiler flags that are specific to particular CIME-driven components are listed last.

The order of listing is a convention and not a requirement.

The possible elements and attributes that can exist in the file are documented in **$CIME/config/xml_schemas/config_compilers_v2.xsd**.

To clarify several conventions:

- The ``<append>`` element implies that any previous definition of that element's parent will be appended with the new element value.
  As an example, the following entry in **config_compilers.xml** would append the value of ``CPPDEFS`` with ``-D $OS`` where ``$OS`` is the environment value of ``OS``.

  ::

     <compiler>
        <CPPDEFS>
            <append> -D<env>OS</env> </append>
        </CPPDEFS>
     </compiler>

- The ``<base>`` element overwrites its parent element's value. For example, the following entry would overwrite the ``CONFIG_ARGS`` for machine ``melvin`` with a ``gnu`` compiler to be ``--host=Linux``.

  ::

     <compiler MACH="melvin" COMPILER="gnu">
        <CONFIG_ARGS>
           <base> --host=Linux </base>
        </CONFIG_ARGS>
     </compiler>


.. _customizing-batch-file:

===================================
config_batch.xml - batch directives
===================================

The **config_batch.xml** schema is defined in **$CIMEROOT/config/xml_schemas/config_batch.xsd**.

CIME supports these batch systems: pbs, cobalt, lsf and slurm.

As is the case for **config_compilers.xml**, the entries in **config_batch.xml** are hierarchical.

#. General configurations for each system are provided at the top of the file.

#. Specific modifications for a given machine are provided below.  In particular each machine should define its own queues.

#. Following is a machine-specific queue section.  This section details the parameters for each queue on the target machine.

#. The last section describes several things:

   - each job that will be submitted to the queue for a CIME workflow,

   - the template file that will be used to generate that job,

   - the prerequisites that must be met before the job is submitted, and

   - the dependencies that must be satisfied before the job is run.

By default the CIME workflow consists of two jobs (**case.run**, **case.st_archive**).

In addition, there is **case.test** job that is used by the CIME system test workflow.

====================================================
Validating your port
====================================================

The following port validation is recommended for any new machine.
Carrying out these steps does not guarantee the model is running properly in all cases nor that the model is scientifically valid on the new machine.

In addition to these tests, detailed validation should be carried out for any new production run.
That means verifying that model restarts are bit-for-bit identical with a baseline run, that the model is bit-for-bit reproducible when identical cases are run for several months, and that production cases are monitored carefully as they integrate forward to identify any potential problems as early as possible. Users are responsible for their own validation process, especially with respect to science validation.

These are the recommended steps for validating a port:

1. Verify functionality by performing these `functionality tests <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_:

::

   ERS_D.f19_g16.X
   ERS_D.T31_g37.A
   ERS_D.f19_g16.B1850CN
   ERI.ne30_g16.X
   ERI.T31_g37.A
   ERI.f19_g16.B1850CN
   ERS.ne30_ne30.F
   ERS.f19_g16.I
   ERS.T62_g16.C
   ERS.T62_g16.DTEST
   ERT.ne30_g16.B1850CN


2. Verify performance and scaling analysis.

   a. Create one or two `load-balanced <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ configurations to check into ``Machines/config_pes.xml`` for the new machine.

   b. Verify that performance and scaling are reasonable.

   c. Review timing summaries in ``$CASEROOT`` for load balance and throughput.

   d. Review coupler "daily" timing output for timing inconsistencies.
      As mentioned in `load balancing a case <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_, useful timing information is contained in a **cpl.log.$date** file that is produced for every run.
      The file contains the run time for each model day during the model run.
      This diagnostic is output as the model runs.
      Searc for ``tStamp`` in this file to see this information.
      The timing information is useful for tracking down temporal variability in model cost due to either inherent model variability cost (I/O, spin-up, seasonal, and so on) or hardware.
      The model daily cost generally is pretty constant unless I/O is written intermittently, such as at the end of the month.

3. Perform validation (both functional and scientific):

   a. Perform a new CIME validation test (**TODO: fill this in**)

   b. Follow the `CCSM4.0 CICE port-validation procedure <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.

   c. Follow the `CCSM4.0 POP2 port-validation procedure <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.

4. Perform two, one-year runs (using the expected load-balanced configuration) as separate job submissions and verify that atmosphere history files are BFB for the last month.
   Do this after some performance testing is complete; you can also combine this with the production test by running the first year as a single run and the second year as a multi-submission production run.
   This will test reproducibility, exact restart over the one-year timescale, and production capability all in one test.

5. Carry out a 20- to 30-year 1.9x2.5_gx1v6 resolution, B_1850_CN compset simulation and compare the results with the diagnostics plots for the 1.9x2.5_gx1v6 Pre-Industrial Control (see the `CCSM4.0 diagnostics <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_).
   Model output data for these runs will be available on the `Earth System Grid (ESG) <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ as well.




