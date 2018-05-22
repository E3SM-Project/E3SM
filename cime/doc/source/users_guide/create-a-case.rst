.. _creating-a-case:

*********************************
Creating a Case
*********************************

Creating a CIME experiment or *case* requires, at a minimum, specifying a compset and a model grid and a case directory.

===================================
Calling **create_newcase**
===================================

The first step in creating a CIME-based experiment is to use **create_newcase**.

If you are not on an out-of-the box CIME-supported platform, you will need to :ref:`port <porting>` CIME to your system before proceeding.

Review the input options for **create_newcase** in the  **help** text.::

  > create_newcase --help

The only required arguments to **create_newcase** are shown here::

  > create_newcase --case [CASE] --compset [COMPSET] --res [GRID]

CIME supports out-of-the-box *component sets*, *model grids* and *hardware platforms* (machines).

The [CASE] argument must be a string and may not contain any of the following special characters
  > + * ? < > / { } [ ] ~ ` @ :

======================================
Results of calling **create_newcase**
======================================

Following is a simple example of using **create_newcase** with aliases for both compset and grid names.
The complete example appears in the :ref:`basic example <use-cases-basic-example>`.

Here, ``$CIMEROOT`` is the full pathname of the root directory of the CIME distribution::

  > cd $CIMEROOT/scripts
  > create_newcase --case ~/cime/example1 --compset A --res f09_g16_rx1

In the example, the command creates a ``$CASEROOT`` directory: **~/cime/example1**. If that directory already exists, a warning is printed and the command aborts. Additional details:

- ``$CASE`` can include letters, numbers,  ".", and "_". In the example, it is ``example1``.

- The compset is ``2000_DATM%NYF_SLND_DICE%SSMI_DOCN%DOM_DROF%NYF_SGLC_SWAV``.

- The model resolution is ``a%0.9x1.25_l%0.9x1.25_oi%gx1v6_r%r05_m%gx1v6_g%null_w%null``.

- **create_newcase** installs files in ``$CASEROOT`` to build and run the model and to optionally archive the case on the target platform.

Running **create_newcase** creates various scripts, files and directories ``$CASEROOT``, as shown here.

- ``user scripts``

   ====================  =====================================================================================================
   case.setup            Script used to set up the case (create the case.run script, the Macros file and user_nl_xxx files).
   case.build            Script to build component and utility libraries and model executable.
   case.submit           Script to submit the case to run using the machine's batch queueing system.
   case.cmpgen_namelist  Script to perform namelist baseline operations (compare, generate, or both).
   xmlchange             Script to modify values in the xml files.
   xmlquery              Script to query values in the xml files.
   preview_namelists     Script for users to see their component namelists in ``$CASEROOT/CaseDocs`` before running the model.
   preview_run           Script for users to see batch submit and mpirun command.
   check_input_data      Script for checking for various input data sets and moving them into place.
   check_case            Script to verify case is set up correctly
   pelayout              Script to query and modify the NTASKS, ROOTPE, and NTHRDS for each component model.  This a convenience script that can be used in place of xmlchange and xmlquery.
   ====================  =====================================================================================================

- ``XML files``

   =====================  ===============================================================================================================================
   env_archive.xml        Defines patters of files to be sent to the short-term archive.
   env_mach_specific.xml  Sets a number of machine-specific environment variables for building and/or running.

                          You can edit this file at any time.

   env_case.xml           Sets case specific variables (e.g. model components, model and case root directories).

                          Cannot be modified after a case has been created.

                          To make changes, your should re-run **create_newcase** with different options.
   env_build.xml          Sets model build settings.

                          This includes component resolutions and component compile-time configuration options.
                          You must run the case.build command after changing this file.

   env_mach_pes.xml       Sets component machine-specific processor layout (see :ref:`changing pe layout<changing-the-pe-layout>` ).

                          The settings in this are critical to a well-load-balanced simulation (see :ref:`load balancing <optimizing-processor-layout>`).
   env_run.xml            Sets runtime settings such as length of run, frequency of restarts, output of coupler diagnostics,
                          and short-term and long-term archiving.  This file can be edited at any time before a job starts.

   env_batch.xml          Sets batch system settings such as wallclock time and queue name.

   =====================  ===============================================================================================================================

- ``User Source Mods Directory``

   =====================  ===============================================================================================================================
   SourceMods             Top-level directory containing subdirectories for each compset component where
                          you can place modified source code for that component.  You may also place modified
			  buildnml and buildlib scripts here.
   =====================  ===============================================================================================================================

- ``Provenance``

   =====================  ===============================================================================================================================
   README.case            File detailing **create_newcase** usage. This is a good place to keep track of runtime problems and changes.
   CaseStatus             File containing a list of operations done in the current case.
   =====================  ===============================================================================================================================

- ``non-modifiable work directories``

   =====================  ===============================================================================================================================
   Buildconf/             Work directory containing scripts to generate component namelists and component and utility libraries
                          (PIO or MCT, for example). You should never have to edit the contents of this directory.
   LockedFiles/           Work directory that holds copies of files that should not be changed.

                          Certain xml files are *locked* after their variables have been used by should no longer be changed.

                          CIME does this by *locking* a file and not permitting you to modify that file unless, depending on the file,                              ``case.setup --clean`` or  ``case.build --clean`` is called.

   Tools/                 Work directory containing support utility scripts. You should never need to edit the contents of this directory.
   =====================  ===============================================================================================================================

The ``$CASEROOT`` xml files are organized so that variables can be locked at certain points after they have been resolved (used) in other parts of the scripts system.

CIME does the following:

- Locks variables in **env_case.xml** after **create_newcase**.

- Locks variables in **env_mach_pes.xml** after **case.setup**.

- Locks variables in **env_build.xml** after completion of **case.build**.

Variables in **env_run.xml**, **env_batch.xml** and **env_archive.xml** are never locked, and most can be changed at any time. There are some exceptions in the **env_batch.xml** file.

The **env_case.html** file can never be unlocked.

These other files can be "unlocked" as follows:

- To unlock **env_mach_pes.xml**, run ``case.setup --clean``.

- To unlock **env_build.xml**, run ``case.build --clean``.
