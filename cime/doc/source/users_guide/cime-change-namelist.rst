.. _namelist-gen:

Customizing your input variables
================================

CIME and CIME-compliant components primarily use Fortran namelists to control runtime options.  Some components use
other text-based files for runtime options.

All CIME-compliant components generate their input variable files using a **buildnml** script typically located in the
component's **cime_config** directory (or other location as set in **config_file.xml**).
**buildnml** may call other scripts to complete construction of the input file.

For example, the CIME data atmosphere model (DATM) generates namelists using the script **$CIMEROOT/components/data_comps/datm/cime_config/buildnml**.

You can customize a model's namelists in one of two ways:

1. by editing the **$CASEROOT/user_nl_xxx** files

  These files should be modified via keyword-value pairs that correspond to new namelist or input data settings.  They use the
  syntax of Fortran namelists.

2. by calling `xmlchange <../Tools_user/xmlchange.html>`_ to modify xml variables in your ``$CASEROOT``.

   Many of these variables are converted to Fortran namelist values for input by the models.  Variables that have
   to be coordinated between models in a coupled system (such as how many steps to run for) are usually in a CIME xml file.

You can generate the component namelists by running `preview_namelists <../Tools_user/preview_namelists.html>`_  from ``$CASEROOT``.

This results in the creation of component namelists (for example, atm_in, lnd_in, and so on) in ``$CASEROOT/CaseDocs/``.

.. warning:: The namelist files in ``CaseDocs`` are  there only for user reference and **SHOULD NOT BE EDITED** since they are overwritten every time `preview_namelists <../Tools_user/preview_namelists.html>`_ and `case.submit <../Tools_user/case.submit.html>`_ are called and the files read at runtime are not the ones in ``CaseDocs``.

.. _use-cases-modifying-driver-namelists:

Customizing driver input variables
-------------------------------------------

The driver input namelists/variables are contained in the files, **drv_in**, **drv_flds_in** and **seq_maps.rc**. Note that **seq_maps.rc** has a different file format than the other two input files.

All driver namelist variables are defined in the file **$CIMEROOT/src/drivers/mct/cime_config/namelist_definition_drv.xml**.

The variables that can be changed only by modifying xml variables appear with the *entry* attribute ``modify_via_xml="xml_variable_name"``.

All other driver namelist variables can be modified by by adding a keyword value pair at the end of ``user_nl_cpl``.

For example, to change the driver namelist value of ``eps_frac`` to ``1.0e-15``, add the following line to the end of the ``user_nl_cpl``:

::

   eps_frac = 1.0e-15

On the hand, to change the driver namelist value of the starting year/month/day, ``start_ymd`` to ``18500901``, use the command:

::

   ./xmlchange RUN_STARTDATE=1850-09-01

Note that

To see the result of change, call `preview_namelists <../Tools_user/preview_namelists.html>`_  and verify that the new value appears in **CaseDocs/drv_in**.

.. _basic_example:

Setting up a multi-year run
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This shows all of the steps necessary to do a multi-year simulation starting from a "cold start" for all components.  The
compset and resolution in this example are for a CESM fully-coupled case but the steps are similar for other models and cases.

1. Create a new case named EXAMPLE_CASE in your **$HOME** directory.

   ::

      > cd $CIME/scripts
      > ./create_newcase --case ~/EXAMPLE_CASE --compset B1850 --res f09_g17

2. Check the pe-layout by running **./pelayout**. Make sure it is suitable for your machine.
   If it is not use `xmlchange <../Tools_user/xmlchange.html>`_ or  `pelayout <../Tools_user/pelayout.html>`_ to modify your pe-layout.
   Then setup your case and build your executable.

   ::

      > cd ~/EXAMPLE_CASE
      > ./case.setup
      > ./case.build

   .. warning:: The case.build script can be compute intensive and may not be suitable to run on a login node. As an alternative you would submit this job to an interactive queue.
                For example, on the NCAR cheyenne platform, you would use **qcmd -- ./case.build** to do this.

3. In your case directory, set the job to run 12 model months, set the wallclock time, and submit the job.

   ::

      > ./xmlchange STOP_OPTION=nmonths
      > ./xmlchange STOP_N=12
      > ./xmlchange JOB_WALLCLOCK_TIME=06:00 --subgroup case.run
      > ./case.submit

4. Make sure the run succeeded.

   You should see the following line or similar at the end of the **cpl.log** file in your run directory or your short term archiving directory, set by ``$DOUT_S_ROOT``.

   ::

      (seq_mct_drv): ===============       SUCCESSFUL TERMINATION OF CPL7-cesm ===============

5. In the same case directory, Set the case to resubmit itself 10 times so it will run a total of 11 years (including the initial year), and resubmit the case. (Note that a resubmit will automatically change the run to be a continuation run).

   ::

      > ./xmlchange RESUBMIT=10
      > ./case.submit

   By default resubmitted runs are not submitted until the previous run is completed.  For 10 1-year runs as configured in this
   example, CIME will first submit a job for one year, then when that job completes it will submit a job for another year.  There will be
   only one job in the queue at a time.
   To change this behavior, and submit all jobs at once (with batch dependencies such that only one job is run at a time), use the command:

   ::

      > ./case.submit --resubmit-immediate

Setting up a branch or hybrid run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A branch or hybrid run uses initialization data from a previous run. Here is an example in which a valid load-balanced scenario is assumed.

1. The first step in setting up a branch or hybrid run is to create a new case. A CESM compset and resolution is assumed below.

   ::

      > cd $CIMEROOT/scripts
      > create_newcase --case ~/NEW_CASE --compset B1850 --res f09_g17
      > cd ~/NEW_CASE


2. For a branch run, use the following `xmlchange <../Tools_user/xmlchange.html>`_  commands to make **NEW_CASE** be a branch off of **EXAMPLE_CASE** at year 0001-02-01.

   ::

      > ./xmlchange RUN_TYPE=branch
      > ./xmlchange RUN_REFCASE=EXAMPLE_CASE
      > ./xmlchange RUN_REFDATE=0001-02-01

3. For a hybrid run, use the following `xmlchange <../Tools_user/xmlchange.html>`_  command to start **NEW_CASE** from **EXAMPLE_CASE** at year 0001-02-01.

   ::

      > ./xmlchange RUN_TYPE=hybrid
      > ./xmlchange RUN_REFCASE=EXAMPLE_CASE
      > ./xmlchange RUN_REFDATE=0001-02-01

   For a branch run, your **env_run.xml** file for **NEW_CASE** should be identical to the file for **EXAMPLE_CASE** except for the ``$RUN_TYPE`` setting.

   Also, modifications introduced into **user_nl_** files in **EXAMPLE_CASE** should be reintroduced in **NEW_CASE**.

4. Next, set up and build your case executable.
   ::

      > ./case.setup
      > ./case.build

5. Pre-stage the necessary restart/initial data in ``$RUNDIR``. Assume for this example that it was created in the **/rest/0001-02-01-00000** directory shown here:

   ::
      > cd $RUNDIR
      > cp /user/archive/EXAMPLE_CASE/rest/0001-02-01-00000/* .

   It is assumed that you already have a valid load-balanced scenario.
   Go back to the case directory, set the job to run 12 model months, and submit the job.
   ::

      > cd ~/NEW_CASE
      > ./xmlchange STOP_OPTION=nmonths
      > ./xmlchange STOP_N=12
      > ./xmlchange JOB_WALLCLOCK_TIME=06:00
      > ./case.submit

6.  Make sure the run succeeded (see above directions) and then change
    the run to a continuation run. Set it to resubmit itself 10 times
    so it will run a total of 11 years (including the initial year),
    then resubmit the case.
    ::

       > ./xmlchange CONTINUE_RUN=TRUE
       > ./xmlchange RESUMIT=10
       > ./case.submit

.. _changing-data-model-namelists:

Customizing data model input variable and stream files
------------------------------------------------------

Each data model can be runtime-configured with its own namelist.

Data Atmosphere (DATM)
~~~~~~~~~~~~~~~~~~~~~~

DATM is discussed in detail in :ref:`data atmosphere overview <data-atm>`.
DATM can be user-customized by changing either its  *namelist input files* or its *stream files*.
The namelist file for DATM is **datm_in** (or **datm_in_NNN** for multiple instances).

- To modify **datm_in** or **datm_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the **user_nl_datm** file or the **user_nl_datm_NNN** file in ``$CASEROOT``.

- To modify the contents of a DATM stream file, first run `preview_namelists <../Tools_user/preview_namelists.html>`_ to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp datm.streams.txt.[extension] user_datm.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_datm.streams.txt[extension``
  3. Edit the **user_datm.streams.txt.*** file.

**Example**

If the stream txt file is **datm.streams.txt.CORE2_NYF.GISS**, the modified copy should be **user_datm.streams.txt.CORE2_NYF.GISS**.
After calling `preview_namelists <../Tools_user/preview_namelists.html>`_ again, your edits should appear in **CaseDocs/datm.streams.txt.CORE2_NYF.GISS**.

Data Ocean (DOCN)
~~~~~~~~~~~~~~~~~~~~~~

DOCN is discussed in detail in :ref:`data ocean overview <data-ocean>`.
DOCN can be user-customized by changing either its namelist input or its stream files.
The namelist file for DOCN is **docn_in** (or **docn_in_NNN** for multiple instances).

- To modify **docn_in** or **docn_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the file in ``$CASEROOT``.

- To modify the contents of a DOCN stream file, first run `preview_namelists <../Tools_user/preview_namelists.html>`_ to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp docn.streams.txt.[extension] user_docn.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_docn.streams.txt[extension``
  3. Edit the **user_docn.streams.txt.*** file.

**Example**

As an example, if the stream text file is **docn.stream.txt.prescribed**, the modified copy should be **user_docn.streams.txt.prescribed**.
After changing this file and calling `preview_namelists <../Tools_user/preview_namelists.html>`_ again, your edits should appear in **CaseDocs/docn.streams.txt.prescribed**.

Data Sea-ice (DICE)
~~~~~~~~~~~~~~~~~~~~~~

DICE is discussed in detail in :ref:`data sea-ice overview <data-seaice>`.
DICE can be user-customized by changing either its namelist input or its stream files.
The namelist file for DICE is ``dice_in`` (or ``dice_in_NNN`` for multiple instances) and its values can be changed by editing the ``$CASEROOT`` file ``user_nl_dice`` (or ``user_nl_dice_NNN`` for multiple instances).

- To modify **dice_in** or **dice_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the file in ``$CASEROOT``.

- To modify the contents of a DICE stream file, first run `preview_namelists <../Tools_user/preview_namelists.html>`_ to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp dice.streams.txt.[extension] user_dice.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_dice.streams.txt[extension``
  3. Edit the **user_dice.streams.txt.*** file.

Data Land (DLND)
~~~~~~~~~~~~~~~~~~~~~~

DLND is discussed in detail in :ref:`data land overview <data-lnd>`.
DLND can be user-customized by changing either its namelist input or its stream files.
The namelist file for DLND is ``dlnd_in`` (or ``dlnd_in_NNN`` for multiple instances) and its values can be changed by editing the ``$CASEROOT`` file ``user_nl_dlnd`` (or ``user_nl_dlnd_NNN`` for multiple instances).

- To modify **dlnd_in** or **dlnd_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the file in ``$CASEROOT``.

- To modify the contents of a DLND stream file, first run `preview_namelists <../Tools_user/preview_namelists.html>`_ to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp dlnd.streams.txt.[extension] user_dlnd.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_dlnd.streams.txt[extension``
  3. Edit the **user_dlnd.streams.txt.*** file.

Data River (DROF)
~~~~~~~~~~~~~~~~~~~~~~

DROF is discussed in detail in :ref:`data river overview <data-river>`.
DROF can be user-customized by changing either its namelist input or its stream files.
The namelist file for DROF is ``drof_in`` (or ``drof_in_NNN`` for multiple instances) and its values can be changed by editing the ``$CASEROOT`` file ``user_nl_drof`` (or ``user_nl_drof_NNN`` for multiple instances).

- To modify **drof_in** or **drof_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the file in ``$CASEROOT``.

- To modify the contents of a DROF stream file, first run `preview_namelists <../Tools_user/preview_namelists.html>`_ to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp drof.streams.txt.[extension] user_drof.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_drof.streams.txt[extension``
  3. Edit the **user_drof.streams.txt.*** file.


Customizing CESM active component-specific namelist settings
------------------------------------------------------------

CAM
~~~

CIME calls **$SRCROOT/components/cam/cime_config/buildnml** to generate the CAM's namelist variables.

CAM-specific CIME xml variables are set in **$SRCROOT/components/cam/cime_config/config_component.xml** and are used by CAM's **buildnml** script to generate the namelist.

For complete documentation of namelist settings, see `CAM namelist variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.

To modify CAM namelist settings, add the appropriate keyword/value pair at the end of the **$CASEROOT/user_nl_cam** file. (See the documentation for each file at the top of that file.)

For example, to change the solar constant to 1363.27, modify **user_nl_cam** file to contain the following line at the end:
::

 solar_const=1363.27

To see the result, call `preview_namelists <../Tools_user/preview_namelists.html>`_ and verify that the new value appears in **CaseDocs/atm_in**.

CLM
~~~

CIME calls **$SRCROOT/components/clm/cime_config/buildnml** to generate the CLM namelist variables.

CLM-specific CIME xml variables are set in **$SRCROOT/components/clm/cime_config/config_component.xml** and are used by CLM's **buildnml** script to generate the namelist.

For complete documentation of namelist settings, see `CLM namelist variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.

To modify CLM namelist settings, add the appropriate keyword/value pair at the end of the **$CASEROOT/user_nl_clm** file.

To see the result, call `preview_namelists <../Tools_user/preview_namelists.html>`_ and verify that the changes appear correctly in **CaseDocs/lnd_in**.

MOSART
~~~~~~

CIME calls **$SRCROOT/components/mosart/cime_config/buildnml** to generate the MOSART namelist variables.

To modify MOSART namelist settings, add the appropriate keyword/value pair at the end of the **$CASEROOT/user_nl_rtm** file.

To see the result of your change, call `preview_namelists <../Tools_user/preview_namelists.html>`_ and verify that the changes appear correctly in **CaseDocs/rof_in**.

CICE
~~~~

CIME calls **$SRCROOT/components/cice/cime_config/buildnml** to generate the CICE namelist variables.

For complete documentation of namelist settings, see `CICE namelist variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.

To modify CICE namelist settings, add the appropriate keyword/value pair at the end of the **$CASEROOT/user_nl_cice** file.
(See the documentation for each file at the top of that file.)
To see the result of your change, call `preview_namelists <../Tools_user/preview_namelists.html>`_ and verify that the changes appear correctly in **CaseDocs/ice_in**.

In addition, `case.setup <../Tools_user/case.setup.html>`_  creates CICE's compile time `block decomposition variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ in **env_build.xml** as follows:

POP2
~~~~

CIME calls **$SRCROOT/components/pop2/cime_config/buildnml** to generate the POP2 namelist variables.

For complete documentation of namelist settings, see `POP2 namelist variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.

To modify POP2 namelist settings, add the appropriate keyword/value pair at the end of the **$CASEROOT/user_nl_pop2** file.
(See the documentation for each file at the top of that file.)
To see the result of your change, call `preview_namelists <../Tools_user/preview_namelists.html>`_ and verify that the changes appear correctly in **CaseDocs/ocn_in**.

In addition, `case.setup <../Tools_user/case.setup.html>`_ generates POP2's compile-time `block decomposition variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ in **env_build.xml** as shown here:

CISM
~~~~

See `CISM namelist variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ for a complete description of the CISM runtime namelist variables. This includes variables that appear both in **cism_in** and in **cism.config**.

To modify any of these settings, add the appropriate keyword/value pair at the end of the **user_nl_cism** file. (See the documentation for each file at the top of that file.)
Note that there is no distinction between variables that will appear in **cism_in** and those that will appear in **cism.config**: simply add a new variable setting in **user_nl_cism**, and it will be added to the appropriate place in **cism_in** or **cism.config**.
To see the result of your change, call `preview_namelists <../Tools_user/preview_namelists.html>`_ and verify that the changes appear correctly in **CaseDocs/cism_in** and **CaseDocs/cism.config**.

Some CISM runtime settings are sets via **env_run.xml**, as documented in `CISM runtime variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.
