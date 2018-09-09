.. _setting-up-a-case:

*********************************
Setting up a Case
*********************************

After creating a case, some aspects of the case are fixed (any variables in env_case.xml). Changing the pe-layout
(see :ref:`Changing Pes<defining-pes>`) or some aspects of the batch system you may be using must be modified before running
**case.setup**.

===================================
Calling **case.setup**
===================================

After creating a case or changing aspects of a case, such as the pe-layout, call the `case.setup <../Tools_user/case.setup.html>`_ command from ``$CASEROOT``.
This creates the following additional files and directories in ``$CASEROOT``:

   =============================   ===============================================================================================================================
   .case.run                       A (hidden) file with the commands that will be used to run the model (such as “mpirun”) and any batch directives needed. 
                                   The directive values are generated using the contents
                                   of **env_mach_pes.xml**. Running `case.setup --clean <../Tools_user/case.setup.html>`_  will remove this file. 
				   This file should not be edited directly and instead controlled through XML variables in **env_batch.xml**. It should also
				   *never* be run directly.

   Macros.make                     File containing machine-specific makefile directives for your target platform/compiler.
                                   This file is created if it does not already exist.

                                   The user can modify the file to change certain aspects of the build, such as compiler flags.
                                   Running `case.setup --clean <../Tools_user/case.setup.html>`_  will not remove the file once it has been created.
                                   However. if you remove or rename the Macros.make file, running `case.setup <../Tools_user/case.setup.html>`_ recreates it.

   user_nl_xxx[_NNNN]              Files where all user modifications to component namelists are made.

                                   **xxx** is any one of the set of components targeted for the case.
                                   For example, for a full active CESM compset, **xxx** is cam, clm or rtm, and so on.

                                   NNNN goes from 0001 to the number of instances of that component.
                                   (See :ref:`multiple instances<multi-instance>`)

                                   For a case with 1 instance of each component (default), NNNN will not appear
                                   in the user_nl file names.

                                   A user_nl file of a given name is created only once.

                                   Calling `case.setup --clean <../Tools_user/case.setup.html>`_ will *not remove* any user_nl files.

                                   Changing the number of instances in the **env_mach_pes.xml** file will cause only
                                   new user_nl files to be added to ``$CASEROOT``.

   CaseDocs/                       Directory that contains all the component namelists for the run.

                                   This is for reference only and files in this directory SHOULD NOT BE EDITED since they will
                                   be overwritten at build time and runtime.

   .env_mach_specific.*            Files summarizing the **module load** commands and environment variables that are set when
                                   the scripts in ``$CASEROOT`` are called. These files are not used by the case but can be
                                   useful for debugging **module load** and environment settings.

   software_environment.txt        This file records some aspects of the computing system on which the case is built,
                                   such as the shell environment.
   =============================   ===============================================================================================================================
