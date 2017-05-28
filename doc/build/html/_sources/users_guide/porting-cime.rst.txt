.. _porting:

*********************************************
Porting and Validating CIME on a new Platform
*********************************************

One of the first steps many users will have to address is getting CIME based models running on their local machine. 
This section will describe that process. 

===========================
Required libraries/packages
===========================

- A functioning MPI environment (this is not required if you are only going to run on one core with the CIME mpi-serial library)
- Build tools gmake and cmake 
- A netcdf library version 4.3 or newer built with the same compiler you will use for CIME
- Optionally a pnetcdf library.
  
If you are using MPI, is usually very helpful to assure that you can run a basic mpi parallel program on your machine prior to attempting a CIME port. 
A simple example for assuring this functionality is in :ref:`MPI example <mpi-example>`.

=================
Steps for porting
=================

The crux of porting CIME is for you to create **at a minimum** the xml file, **config_machines.xml** and place these files in your **$HOME/.cime** directory.
In addition, if you have a batch system, you will also need to add a **config_batch.xml** file to your **$HOME/.cime** directory.

.. note::

   All files in **$HOME/.cime/** are appended to the xml objects read into memory from the **$CIME/config/[model]**, , where **[model]** can currently be either ``acme`` or ``cesm``.

A more detailed discussion of creating these files follows.

The following steps should be followed:

#. Create a **$HOME/.cime** directory

#. Create a **$HOME/.cime/config** file. At a minumum this should contain the following two lines:
   ::

      [main]
      CIME_MODEL=cesm

   or 

   ::

      [main]
      CIME_MODEL=acme

#. Create a **$HOME/.cime/config_machines.xml** file.

   This file contains all the information that a user needs to set in order to configure a new machine to be CIME complaint. 

   We provide several level of templates that you can use to create this file.

   - If you are on a MAC and have the libraries installed with either MacPorts or HomeBrew, copy
     **$CIME/config/cesm/machines/userdefined_laptop_template/config_machines.xml** to
     **$HOME/.cime/config_machines.xml**.

   - Otherwise, copy the template file
     **$CIME/config/xml_schemas/config_machines_template.xml** to
     **$HOME/.cime/config_machines.xml**.

   Fill in the contents of **$HOME/.cime/config_machines.xml** that are specific to your machine. 

   For more details see :ref:`customize the config_machines.xml file <customizing-machine-file>`. 

   Check that your **config_machines.xml** file conforms to the CIME schema definition by doing the following:
   ::

      xmllint --noout --schema $CIME/config/xml_schemas/config_machines.xsd $HOME/.cime/config_machines.xml

#. Optionally, create a **$HOME/.cime/config_compilers.xml** file, if you have compiler settings that are specific to your machine.

   The default compiler settings are set in **$CIME/config/[model]/machines/config_compilers.xml**, where **[model]** can currently be either ``acme`` or ``cesm``.

   There is no template for **config_compilers.xml**.

#. Optionally, create a **$HOME/.cime/config_batch.xml** file if you have a batch system.

   Out of the box batch settings are set in **$CIME/config/[model]/machines/config_batch.xml**, where **[model]** can currently be either ``acme`` or ``cesm``.

#. Once you have a basic configuration for your machine defined in your **$HOME/.cime** XML files, you should interactively run **scripts_regression_test.py** from the directory **$CIME/scripts/tests**. 
   This script will run a number of basic unit tests starting from the simplest tests and working toward more complicated ones.

#. Finally when all the previous steps have run correctly, you are ready to try a case at your target compset and resolution.
   Once you have successfully created the required xml files in your .cime directory and are satisfied with the results you can merge them into the default files in the **config/$CIME_MODEL/machines** directory.   
   If you would like to make this machine definition available generally you may then issue a pull request to add your changes to the git repository.  
   
.. _customizing-machine-file:

===========================================
config_machines.xml - machine specific file
===========================================

The machine specific files is defined in the model-specific :ref:`config_machines.xml <defining-machines>`.

The first step a user must take to make their machine CIME-compatible is to add the appropriate entries for their machine in **config_machines.xml**.

Each ``<machine>`` tag requires the following input: 

- ``DESC``: a text description of the machine, this field is current not used
- ``NODENAME_REGEX``: a regular expression used to identify this machine it must work on compute nodes as well as login nodes, use machine option to create_test or create_newcase if this flag is not available 
- ``OS``: the operating system of this machine. 
- ``PROXY``: optional http proxy for access to the internet
- ``COMPILERS``: compilers supported on this machine, comma seperated list, first is default 
- ``MPILIBS``: mpilibs supported on this machine, comma seperated list, first is default 
- ``PROJECT``: A project or account number used for batch jobs can be overridden in environment or $HOME/.cime/config 
- ``SAVE_TIMING_DIR``: (Acme only) directory to write timing output to 
- ``CIME_OUTPUT_ROOT``: Base directory for case output, the bld and run directories are written below here 
- ``DIN_LOC_ROOT``: location of the inputdata directory 
- ``DIN_LOC_ROOT_CLMFORC``: optional input location for clm forcing data  
- ``DOUT_S_ROOT``: root directory of short term archive files 
- ``DOUT_L_MSROOT``: root directory on mass store system of long term archive files
- ``BASELINE_ROOT``:  Root directory for system test baseline files 
- ``CCSM_CPRNC``: location of the cprnc tool, compares model output in testing
- ``GMAKE``: gnu compatible make tool, default is 'gmake' 
- ``GMAKE_J``: optional number of threads to pass to the gmake flag 
- ``TESTS``: (acme only) list of tests to run on this machine 
- ``BATCH_SYSTEM``: batch system used on this machine (none is okay) 
- ``SUPPORTED_BY``: contact information for support for this system 
- ``MAX_TASKS_PER_NODE``: maximum number of threads*tasks per shared memory node on this machine
- ``PES_PER_NODE``: number of physical PES per shared node on this machine, in practice the MPI tasks per node will not exceed this value 
- ``PROJECT_REQUIRED``: Does this machine require a project to be specified to the batch system? 
- ``mpirun``: The mpi exec to start a job on this machine. 
  This is itself an element that has sub elements that must be filled:

  * Must have a required ``<executable>`` element 
  * May have optional attributes of ``compiler``, ``mpilib`` and/or ``threaded``
  * May have an optional ``<arguments>`` element which in turn contain one or more ``<arg>`` elements. 
    These specify the arguments to the mpi executable and as a result are dependent on your mpi library implementation.


- ``module_system``: How and what modules to load on this system. Module systems allow you to easily load multiple compiler environments on a given machine. CIME provides support for two types of module tools: `module <http://www.tacc.utexas.edu/tacc-projects/mclay/lmod>`_ and `soft  <http://www.mcs.anl.gov/hs/software/systems/softenv/softenv-intro.html>`_.   If neither of these are available on your machine, the simply set ``<module_system type="none"\>``.
   
- ``environment_variables``: environment_variables to set on this system. 
   This contains sub elements, ``<env>`` with the ``name`` attribute specifying the environment variable name, and the element value specifying the corresponding environment variable value. If the element value is not set, then the corresponding environment variable will be unset in your shell. 

   As an example, the following sets the environment variable ``OMP_STACKSIZE`` to 256M.
   ::

      <env name="OMP_STACKSIZE">256M</env>

   and the following unsets this environment variable in the shell:
   ::

      <env name="OMP_STACKSIZE"></env>

   .. note:: These changes are **ONLY** activated for the CIME build and run environment, **BUT NOT** for your login shell. To activate them for your login shell, you would source either **$CASEROOT/.env_mach_specific.sh** or **$CASEROOT/.env_mach_specific.csh**, depending on your shell.

.. _customizing-compiler-file:

=================================================
config_compilers.xml - compiler paths and options
=================================================

The **config_compilers.xml** file defines compiler flags for building CIME (and also CESM and ACME prognostic CIME-driven components).  
There is a heirarchy that is intrinsic to this file.

.# General compiler flags (e.g. for the gnu compiler) that are machine and component independent are listed first.

.# Compiler flags specific to a particular operating systems are listed next.

.# Compiler flags that are specific to particular machines are listed next.

.# Compiler flags that are specific to particular CIME-driven components are listed last.

The order of listing is a convention and not a requirement.

The possible elements and attributes that can exist in this file are documented in **$CIME/config/xml_schemas/config_compilers_v2.xsd**.

Its useful to clarify several conventions:

- the ``<append>`` element implies that any previous defintion of that element's parent will be appended with the new element value.
  As an example, the following entry in **config_compilers.xml**, would append the value of ``CPPDEFS`` with ``-D $OS``
  where ``$OS`` is the environment value of ``OS``.
  
  ::
   
     <compiler>
        <CPPDEFS>
	    <append> -D<env>OS</env> </append>
	</CPPDEFS>
     </compiler>

- the ``<base>`` element overwrites its parent element's value with its element text.
  As an example, the following entry would overwrite the ``CONFIG_ARGS`` for machine ``melvin`` with a ``gnu`` compiler to be ``--host=Linux``.

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

CIME currently supports the following batch systems: pbs, cobalt, lsf and slurm.

As is the case for **config_compilers.xml**, the entries in **config_batch.xml** are heirarchical.

.# General configurations for each system are provided at the top of the file.

.# Specific modifications for a given machine are provided below.  In particular each machine should define its own queues.   

.# Following is a machine specific queue section.  This section details the parameters for each queue on the target machine.

.# The last section describes several things:

   - each job that will be submitted to the queue for a CIME workflow

   - the template file that will be used to generate that job

   - the prerequisets that must be met before the job is submitted 

   - the dependancies that must be satisfied before that job is run   

By default the CIME workflow consists of two jobs (**case.run**, **case.st_archive**). 

In addition, there is **case.test** job that is used by the CIME system test workflow. 

====================================================
Validating your port
====================================================

The following port validation is recommended for any new machine. 
Carrying out these steps does not guarantee the model is running properly in all cases nor that the model is scientifically valid on the new machine. 
In addition to these tests, detailed validation should be carried out for any new production run. 
That means verifying that model restarts are bit-for-bit identical with a baseline run, that the model is bit-for-bit reproducible when identical cases are run for several months, and that production cases are monitored very carefully as they integrate forward to identify any potential problems as early as possible. 
These are recommended steps for validating a port and are largely functional tests. 
Users are responsible for their own validation process, especially with respect to science validation.

1. Verify functionality by performing these `functionality tests <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.
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
      As has been mentioned in the section on `load balancing a case <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_, useful timing information is contained in cpl.log.$date file that is produced for every run. 
      The cpl.log file contains the run time for each model day during the model run. 
      This diagnostic is output as the model runs. 
      You can search for tStamp in this file to see this information. 
      This timing information is useful for tracking down temporal variability in model cost either due to inherent model variability cost (I/O, spin-up, seasonal, etc) or possibly due to variability due to hardware. 
      The model daily cost is generally pretty constant unless I/O is written intermittently such as at the end of the month.

3. Perform validation (both functional and scientific):

   a. Perform a new CIME validation test (**TODO: fill this in**)

   b. Follow the `CCSM4.0 CICE port-validation procedure <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.

   c. Follow the `CCSM4.0 POP2 port-validation procedure <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.

4. Perform two, one-year runs (using the expected load-balanced configuration) as separate job submissions and verify that atmosphere history files are bfb for the last month. 
   Do this after some performance testing is complete; you may also combine this with the production test by running the first year as a single run and the second year as a multi-submission production run. 
   This will test reproducibility, exact restart over the one-year timescale, and production capability all in one test.

5. Carry out a 20-30 year 1.9x2.5_gx1v6 resolution, B_1850_CN compset simulation and compare the results with the diagnostics plots for the 1.9x2.5_gx1v6 Pre-Industrial Control (see the `CCSM4.0 diagnostics <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_). 
   Model output data for these runs will be available on the `Earth System Grid (ESG) <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ as well.




