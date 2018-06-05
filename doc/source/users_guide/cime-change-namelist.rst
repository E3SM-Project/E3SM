.. _namelist-gen:

Customizing your input variables
================================

CIME and current CIME-compliant components currently primarily uses fortran namelists.

All CIME-compliant components generate their input variable settings using a **buildnml** file located in the component's **cime_config** directory.

For example, the CIME data atmosphere model (DATM) generates namelists using the script **$CIMEROOT/components/data_comps/datm/cime_config/buildnml**.

You can customize your namelists  in one of two ways:

1. by editing the **$CASEROOT/user_nl_xxx** files

  These files should be modified via keyword-value pairs that correspond to new namelist or input data settings.

2. by calling `xmlchange <../Tools_user/xmlchange.html>`_ to modify xml variables in your ``$CASEROOT``.

You can preview the component namelists by running `preview_namelists <../Tools_user/preview_namelists.html>`_  from ``$CASEROOT``.

This results in the creation of component namelists (for example, atm_in, lnd_in, and so on) in ``$CASEROOT/CaseDocs/``.
The namelist files are there only for user reference and **SHOULD NOT BE EDITED** since they are overwritten every time `preview_namelists <../Tools_user/preview_namelists.html>`_ and `case.submit <../Tools_user/case.submit.html>`_ are called.

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
.. _changing-data-model-namelists:

Customizing data model input variable and stream files
------------------------------------------------------

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
