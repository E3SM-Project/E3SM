.. _building-a-case:

******************
Building a Case
******************

Once the case has been created and setup, its time to build the executable.
Several directories full of source code must be built all with the same compiler and flags.
**case.build** performs all build operations (setting dependecies, invoking Make,
creating the executable).

.. _building-the-model:

========================
Calling **case.build**
========================

After calling `case.setup <../Tools_user/case.setup.html>`_ , run `case.build <../Tools_user/case.build.html>`_  to build the model executable. Running this will:

1. Create the component namelists in ``$RUNDIR`` and ``$CASEROOT/CaseDocs``.
2. Create the necessary compiled libraries used by coupler and component models ``mct``, ``pio``, ``gptl`` and ``csm_share``.
   The libraries will be placed in a path below ``$SHAREDLIBROOT``.
3. Create the necessary compiled libraries for each component model. These are be placed in ``$EXEROOT/bld/lib``.
4. Create the model executable (``$MODEL.exe``), which is placed in ``$EXEROOT``.

You do not need to change the default build settings to create the executable, but it is useful to become familiar with them in order to make optimal use of the system. The CIME scripts provide you with a great deal of flexibility in customizing the build process.

The **env_build.xml** variables control various aspects of building the executable. Most of the variables should not be modified, but users can modify these:

- ``$BUILD_THREADED`` : if TRUE, the model will be built with OpenMP.

- ``$DEBUG`` : if TRUE, the model is compiled with debugging instead of optimization flags.

- ``$GMAKE_J`` : How many threads GNUMake should use while building.

The best way to see what xml variables are in your ``$CASEROOT`` directory is to use the `xmlquery <../Tools_user/xmlquery.html>`_  command. For usage information, run:
::

   > ./xmlquery --help

To build the model, change to your ``$CASEROOT`` directory and execute **case.build**.
::

   > cd $CASEROOT
   > ./case.build

Diagnostic comments appear as the build proceeds.

The `case.build <../Tools_user/case.build.html>`_  command generates the utility and component libraries and the model executable, and it generates build logs for each component.
Each log file is named form: **$component.bldlog.$datestamp**. They are located in ``$BLDDIR``. If they are compressed (as indicated by a .gz file extension), the build ran successfully.

Invoking `case.build <../Tools_user/case.build.html>`_  creates the following directory structure in ``$EXEROOT`` if the Intel compiler is used:
::

   atm/, cpl/, esp/, glc/, ice/, intel/, lib/, lnd/, ocn/, rof/, wav/

Except for **intel/** and **lib/**, each directory contains an **obj/** subdirectory for the target model component's compiled object files.

The *mct*, *pio*, *gptl* and *csm_share* libraries are placed in a directory tree that reflects their dependencies. See the **bldlog** for a given component to locate the library.

Special **include** modules are placed in **lib/include**. The model executable (**cesm.exe** or **e3sm.exe**, for example) is placed directly in ``$EXEROOT``.

Component namelists, component logs, output data sets, and restart files are placed in ``$RUNDIR``.
It is important to note that ``$RUNDIR`` and ``$EXEROOT`` are independent variables that are set in the **$CASEROOT/env_run.xml** file.

.. _rebuilding-the-model:

========================
Rebuilding the model
========================

Rebuild the model under the following circumstances:

If either **env_build.xml** or **Macros.make** has been modified, and/or if code is added to **SourceMods/src.**, it's safest to clean the build and rebuild from scratch as shown here:
::

   > cd $CASEROOT
   > ./case.build --clean-all

If you have ONLY modified the PE layout in **env_mach_pes.xml**, a clean may not be required.
::

   > cd $CASEROOT
   > ./case.build

If the threading has been changed (turned on or off) in any component since the previous build, the build script should fail with the following error and suggestion that the model be rebuilt from scratch:
::

   ERROR SMP STATUS HAS CHANGED
   SMP_BUILD = a0l0i0o0g0c0
   SMP_VALUE = a1l0i0o0g0c0
   A manual clean of your obj directories is strongly recommended.
   You should execute the following:
      ./case.build --clean
      ./case.build

    ---- OR ----

    You can override this error message at your own risk by executing:
      ./xmlchange SMP_BUILD=0
    Then rerun the build script interactively.

If there is any doubt, rebuild.

Run this to clean all of the model components (except for support libraries such as *mct* and *gptl*):
  ::

     > case.build --clean

Run this to clean everything associated with the build:
  ::

     > case.build --clean-all

You can also clean an individual component as shown here, where "compname" is the name of the component you want to clean (for example, atm, clm, pio and so on).
  ::

     > case.build --clean compname

Review the **help** text for more information.

.. _inputdata:

==========
Input data
==========

All active components and data components use input data sets. In order to run CIME and the CIME-compliant active components, a local disk needs the directory tree that is specified by the xml variable ``$DIN_LOC_ROOT`` to be populated with input data.

Input data is provided as part of the CIME release via data from a subversion input data server. It is downloaded from the server on an as-needed basis determined by the case. Data may already exist in the default local file system's input data area as specified by ``$DIN_LOC_ROOT``.

Input data can occupy significant space on a system, so users should share a common ``$DIN_LOC_ROOT`` directory on each system if possible.

The build process handles input data as follows:

- The **buildnml** scripts in the various component ``cime_config`` directories create listings of required component input data sets in the ``Buildconf/$component.input_data_list`` files.

- `check_input_data <../Tools_user/check_input_data.html>`_ , which is called by `case.build <../Tools_user/case.build.html>`_ , checks for the presence of the required input data files in the root directory ``$DIN_LOC_ROOT``.

- If all required data sets are found on the local disk, the build can proceed.

- If any of the required input data sets are not found locally, the
  files that are missing are listed. At this point, you must obtain
  the required data from the input data server with `check_input_data
  <../Tools_user/check_input_data.html>`_ as shown here: ::

     check_input_data --download

The **env_run.xml** variables ``$DIN_LOC_ROOT`` and ``$DIN_LOC_ROOT_CLMFORC`` determine where you should expect input data to reside on a local disk.
