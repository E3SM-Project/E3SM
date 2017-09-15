.. _customizing-a-case:

**************************************************
Customizing a Case
**************************************************

All CIME-compliant components generate their namelist settings using the **cime_config/buildnml** file located in the component's directory tree.
For example, the CIME data atmosphere model (DATM) generates namelists using the script **$CIMEROOT/components/data_comps/datm/cime_config/buildnml**.

User-specific component namelist changes should be made only by:

- editing the **$CASEROOT/user_nl_xxx** files.

- using :ref:`xmlchange<modifying-an-xml-file>` to modify xml variables in **env_run.xml**, **env_build.xml** or **env_mach_pes.xml**.

You can preview the component namelists by running **preview_namelists** from ``$CASEROOT``.
This results in the creation of component namelists (for example, atm_in, lnd_in, and so on) in **$CASEROOT/CaseDocs/**. The namelist files are there only for user reference and SHOULD NOT BE EDITED since they are overwritten every time **preview_namelists**  and  **case.submit** are called.

The following sections summarize how to modify component-specific runtime settings:

.. _modifying-an-xml-file:

=================================================
Modifying an xml file
=================================================

Run the ``$CASEROOT`` script **xmlchange** to modify variables in xml files. The script performs variable error checking.

Here are two examples of how to invoke **xmlchange**:
::

   xmlchange <entry id>=<value>
   -- OR --
   xmlchange -id <entry id> -val <name> -file <filename>

The ``-id`` argument identifies the variable to be changed, and ``-val`` is the intended value of that variable. See the **help** text for more usage information:
::

   xmlchange --help


.. _changing-the-pe-layout:

=================================================
Customizing the PE layout
=================================================

Settings in the **env_mach_pes.xml** file determine:

- the number of MPI tasks and OpenMP threads for each component.
- the number of instances of each component.
- the layout of the components across the hardware processors.

Optimizing the throughput and efficiency of a CIME experiment often involves customizing the processor (PE) layout. (See :ref:`load balancing <optimizing-processor-layout>`.)
CIME provides significant flexibility with respect to the layout of components across different hardware processors. In general, the CIME components -- atm, lnd, ocn, and so on -- can run on overlapping or mutually unique processors. While each component is associated with a unique MPI communicator, the CIME driver runs on the union of all processors and controls the sequencing and hardware partitioning.

The component processor layout is determined by the following settings:

- the number of MPI tasks.
- the number of OpenMP threads per task.
- the root MPI task number from the global communicator.
- the maximum number of MPI tasks per node.

The entries in **env_mach_pes.xml** have the following meanings:

   ================== ================================================================================================
   MAX_TASKS_PER_MODE The total number of (MPI tasks) * (OpenMP threads) allowed on a node.
                         This is defined in **config_machines.xml** and therefore given a default setting, but
                         can be user modified.
   PES_PER_NODE       The maximum number of MPI tasks per node.
                         This is defined in **config_machines.xml** and therefore given a default setting, but
                         can be user modified.
   NTASKS             Total number of MPI tasks.
                          A negative value indicates nodes rather than tasks, where
                          PES_PER_NODE * -NTASKS equals the number of MPI tasks.
   NTHRDS             Number of OpenMP threads per MPI task.
   ROOTPE             The global MPI task of the component root task; if negative, indicates nodes rather than tasks.
   PSTRID             The stride of MPI tasks across the global set of pes (for now set to 1).
   NINST              The number of component instances, which are spread evenly across NTASKS.
   ================== ================================================================================================

----------------
**Example 1**
----------------

If a component has **NTASKS=16**, **NTHRDS=4** and **ROOTPE=32**, it will run on 64 hardware processors using 16 MPI tasks and 4 threads per task starting at global MPI task 32.

Each CIME component has corresponding entries for ``NTASKS``, ``NTHRDS``, ``ROOTPE`` and ``NINST`` in the **env_mach_pes.xml** file.

**Note:**

- ``NTASKS`` must be greater or equal to 1 even for inactive (stub) components.
- ``NTHRDS`` must be greater or equal to 1.
- If ``NTHRDS`` = 1, this generally means threading parallelization will be off for that component.
- ``NTHRDS`` should never be set to zero.
- The total number of hardware processors allocated to a component is ``NTASKS`` * ``NTHRDS``.
- The coupler processor inputs specify the pes used by coupler computation such as mapping, merging, diagnostics, and flux calculation. This is distinct from the driver, which automatically runs on the union of all processors to manage model concurrency and sequencing.
- The root processor is set relative to the MPI global communicator, not the hardware processors counts. An example of this is below.
- The layout of components on processors has no impact on the science.
- If all components have identical ``NTASKS``, ``NTHRDS``, and ``ROOTPE`` settings, all components will run sequentially on the same hardware processors.

The scientific sequencing is hardwired into the driver. Changing processor layouts does not change intrinsic coupling lags or coupling sequencing.

For a **fully active configuration**, the atmosphere component is hardwired in the driver to never run concurrently with the land or ice component. Performance improvements associated with processor layout concurrency therefore are constrained in this case such that there is never a performance reason not to overlap the atmosphere component with the land and ice components. Beyond that constraint, the land, ice, coupler and ocean models can run concurrently, and the ocean model can also run concurrently with the atmosphere model.

An important but often misunderstood point: The root processor for any given component is set relative to the MPI global communicator, not the hardware processor counts. For instance, in the following example, the atmosphere and ocean will run concurrently, each on 64 processors with the atmosphere running on MPI tasks 0-15 and the ocean running on MPI tasks 16-79.
::

   NTASKS(ATM)=6  NTHRRDS(ATM)=4  ROOTPE(ATM)=0
   NTASKS(OCN)=64 NTHRDS(OCN)=1   ROOTPE(OCN)=16

The first 16 tasks are each threaded 4 ways for the atmosphere. CIME ensures that the batch submission script (**$CASE.run**) automatically requests 128 hardware processors, and the first 16 MPI tasks will be laid out on the first 64 hardware processors with a stride of 4. The next 64 MPI tasks are laid out on the second set of 64 hardware processors.

If you had set ``ROOTPE_OCN`` to 64 in this example, a total of 176 processors would be requested, the atmosphere would be laid out on the first 64 hardware processors in 16x4 fashion, and the ocean model would be laid out on hardware processors 113-176. Hardware processors 65-112 would be allocated but completely idle.

----------------
**Example 2**
----------------

If a component has **NTASKS=-2**, **NTHRDS=4** and **ROOTPE=0**, **PES_PER_NODE=4**, **MAX_TASKS_PER_NODE=4**, it will run on (8 MPI tasks * 4 threads) = 32 hardware processors on 8 nodes.

If you intended 2 nodes INSTEAD of 8 nodes, then you would change **PES_PER_NODE=1** (using **xmlchange**).


**Note**: **env_mach_pes.xml** *cannot* be modified after **case.setup** has been invoked without first running the following:
::

   case.setup --clean

.. _changing-driver-namelists:

===================================================
Customizing driver namelists
===================================================

Driver namelist variables belong in two groups:

1. Those that are set directly from ``$CASEROOT`` xml variables.

2. Those that are set by the driver utility **$CIMEROOT/src/drivers/mct/cime_config/buildnml**.

All driver namelist variables are defined in the file **$CIMEROOT/src/drivers/mct/cime_config/namelist_definition_drv.xml**.
The variables that can be changed only by modifying xml variables appear with the *entry* attribute ``modify_via_xml="xml_variable_name"``.

All other variables that appear in the **namelist_definition_drv.xml** file can be modified by adding a keyword value pair at the end of ``user_nl_cpl``.
For example, to change the driver namelist value of ``eps_frac`` to ``1.0e-15``, add the following line to the end of the ``user_nl_cpl``:
::

   eps_frac = 1.0e-15

To see the result of change, call **preview_namelists** and verify that the new value appears in **CaseDocs/drv_in**.

.. _changing-data-model-namelists:

===================================================
Customizing data model namelists and stream files
===================================================
------------------------
Data Atmosphere (DATM)
------------------------

DATM is discussed in detail in :ref:`data atmosphere overview <data-atm>`.
DATM can be user-customized by changing either its  *namelist input files* or its *stream files*.
The namelist file for DATM is **datm_in** (or **datm_in_NNN** for multiple instances).

- To modify **datm_in** or **datm_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the **user_nl_datm** file or the **user_nl_datm_NNN** file in ``$CASEROOT``.

- To modify the contents of a DATM stream file, first run **preview_namelists** to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp datm.streams.txt.[extension] user_datm.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_datm.streams.txt[extension``
  3. Edit the **user_datm.streams.txt.*** file.

**Example**

If the stream txt file is **datm.streams.txt.CORE2_NYF.GISS**, the modified copy should be **user_datm.streams.txt.CORE2_NYF.GISS**.
After calling **preview_namelists** again, your edits should appear in **CaseDocs/datm.streams.txt.CORE2_NYF.GISS**.

------------------------
Data Ocean (DOCN)
------------------------

DOCN is discussed in detail in :ref:`data ocean overview <data-ocean>`.
DOCN can be user-customized by changing either its namelist input or its stream files.
The namelist file for DOCN is **docn_in** (or **docn_in_NNN** for multiple instances).

- To modify **docn_in** or **docn_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the file in ``$CASEROOT``.

- To modify the contents of a DOCN stream file, first run **preview_namelists** to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp docn.streams.txt.[extension] user_docn.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_docn.streams.txt[extension``
  3. Edit the **user_docn.streams.txt.*** file.

**Example**

As an example, if the stream text file is **docn.stream.txt.prescribed**, the modified copy should be **user_docn.streams.txt.prescribed**.
After changing this file and calling **preview_namelists** again, your edits should appear in **CaseDocs/docn.streams.txt.prescribed**.

------------------------
Data Sea-ice (DICE)
------------------------

DICE is discussed in detail in :ref:`data sea-ice overview <data-seaice>`.
DICE can be user-customized by changing either its namelist input or its stream files.
The namelist file for DICE is ``dice_in`` (or ``dice_in_NNN`` for multiple instances) and its values can be changed by editing the ``$CASEROOT`` file ``user_nl_dice`` (or ``user_nl_dice_NNN`` for multiple instances).

- To modify **dice_in** or **dice_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the file in ``$CASEROOT``.

- To modify the contents of a DICE stream file, first run **preview_namelists** to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp dice.streams.txt.[extension] user_dice.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_dice.streams.txt[extension``
  3. Edit the **user_dice.streams.txt.*** file.

------------------
Data Land (DLND)
------------------

DLND is discussed in detail in :ref:`data land overview <data-lnd>`.
DLND can be user-customized by changing either its namelist input or its stream files.
The namelist file for DLND is ``dlnd_in`` (or ``dlnd_in_NNN`` for multiple instances) and its values can be changed by editing the ``$CASEROOT`` file ``user_nl_dlnd`` (or ``user_nl_dlnd_NNN`` for multiple instances).

- To modify **dlnd_in** or **dlnd_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the file in ``$CASEROOT``.

- To modify the contents of a DLND stream file, first run **preview_namelists** to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp dlnd.streams.txt.[extension] user_dlnd.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_dlnd.streams.txt[extension``
  3. Edit the **user_dlnd.streams.txt.*** file.

------------------
Data River (DROF)
------------------

DROF is discussed in detail in :ref:`data river overview <data-river>`.
DROF can be user-customized by changing either its namelist input or its stream files.
The namelist file for DROF is ``drof_in`` (or ``drof_in_NNN`` for multiple instances) and its values can be changed by editing the ``$CASEROOT`` file ``user_nl_drof`` (or ``user_nl_drof_NNN`` for multiple instances).

- To modify **drof_in** or **drof_in_NNN**, add the appropriate keyword/value pair(s) for the namelist changes that you want at the end of the file in ``$CASEROOT``.

- To modify the contents of a DROF stream file, first run **preview_namelists** to list the *streams.txt* files in the **CaseDocs/** directory. Then, in the same directory:

  1. Make a *copy* of the file with the string *"user_"* prepended.
        ``> cp drof.streams.txt.[extension] user_drof.streams.txt[extension.``
  2. **Change the permissions of the file to be writeable.** (chmod 644)
        ``chmod 644 user_drof.streams.txt[extension``
  3. Edit the **user_drof.streams.txt.*** file.

=================================================================
Customizing CESM active component-specific namelist settings
=================================================================

---
CAM
---

CAM's `configure <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ and `build-namelist <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ utilities are called by ``Buildconf/cam.buildnml.csh``. The folllowing are used to set compset variables (for example, "-phys cam5" for CAM_CONFIG_OPTS) and in general should not be modified for supported compsets:
::

  `CAM_CONFIG_OPTS <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_
  `CAM_NAMELIST_OPTS <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_
  `CAM_NML_USECASE <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_

For complete documentation of namelist settings, see `CAM namelist variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.

To modify CAM namelist settings, add the appropriate keyword/value pair at the end of the **$CASEROOT/user_nl_cam** file. (See the documentation for each file at the top of that file.)

For example, to change the solar constant to 1363.27, modify **user_nl_cam** file to contain the following line at the end:
::

 solar_const=1363.27

To see the result, call **preview_namelists** and verify that the new value appears in **CaseDocs/atm_in**.

---
CLM
---

CIME calls **$SRCROOT/components/clm/cime_config/buildnml** to generate the CLM namelist variables.
CLM-specific CIME xml variables are set in **$SRCROOT/components/clm/cime_config/config_component.xml** and are used by CLM's **buildnml** script to generate the namelist.

For complete documentation of namelist settings, see `CLM namelist variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.

To modify CLM namelist settings, add the appropriate keyword/value pair at the end of the **$CASEROOT/user_nl_clm** file. To see the result, call **preview_namelists** and verify that the changes appear correctly in **CaseDocs/lnd_in**.

---
RTM
---

CIME calls **$SRCROOT/components/rtm/cime_config/buildnml** to generate the RTM namelist variables.

For complete documentation of namelist settings, see RTM namelist variables. //SHOULD THERE BE A LINK HERE?//

To modify RTM namelist settings, add the appropriate keyword/value pair at the end of the **$CASEROOT/user_nl_rtm** file. To see the result of your change, call **preview_namelists** and verify that the changes appear correctly in **CaseDocs/rof_in**.

---
CICE
---

The CICE `configure <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ and `build-namelist <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ utilities are called by **Buildconf/cice.buildnml.csh**. Note that `CICE_CONFIG_OPTS <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ and `CICE_NAMELIST_OPTS <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ are used to set compset-specific variables and in general should not be modified for supported compsets.

For complete documentation of namelist settings, see `CICE namelist variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.

To modify CICE namelist settings, add the appropriate keyword/value pair at the end of the **$CASEROOT/user_nl_cice** file. (See the documentation for each file at the top of that file.) To see the result of your change, call **preview_namelists** and verify that the changes appear correctly in **CaseDocs/ice_in**.

In addition, **case.setup** creates CICE's compile time `block decomposition variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ in **env_build.xml** as follows:
::

   ./case.setup
     ?
   Buildconf/cice.buildnml.csh and $NTASKS_ICE and $NTHRDS_ICE
     ?
   env_build.xml variables CICE_BLCKX, CICE_BLCKY, CICE_MXBLCKS, CICE_DECOMPTYPE
   CPP variables in cice.buildexe.csh

----
POP2
----
See `POP2 namelist variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ for complete description of the POP2 runtime namelist variables. Note that `OCN_COUPLING, OCN_ICE_FORCING andOCN_TRANSIENT <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ are normally used ONLY to set compset-specific variables and should not be edited. For complete documentation of namelist settings, see `CICE namelist variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_.

To modify POP2 namelist settings, add the appropriate keyword/value pair at the end of the **$CASEROOT/user_nl_pop2** file. (See the documentation for each file at the top of that file.) To see the result of your change, call **preview_namelists** and verify that the changes appear correctly in **CaseDocs/ocn_in**.

In addition, **cesm_setup** generates POP2's compile-time `block decomposition variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ in **env_build.xml** as shown here:
::

   ./cesm_setup
       ?
   Buildconf/pop2.buildnml.csh and $NTASKS_OCN and $NTHRDS_OCN
       ?
   env_build.xml variables POP2_BLCKX, POP2_BLCKY, POP2_MXBLCKS, POP2_DECOMPTYPE
   CPP variables in pop2.buildexe.csh

CISM
----
See `CISM namelist variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ for a complete description of the CISM runtime namelist variables. This includes variables that appear both in **cism_in** and in **cism.config**.

To modify any of these settings, add the appropriate keyword/value pair at the end of the **user_nl_cism** file. (See the documentation for each file at the top of that file.) To see the result of your change, call **preview_namelists** and verify that the changes appear correctly in **CaseDocs/cism_in** and **CaseDocs/cism.config**.

Some CISM runtime settings are sets via **env_run.xml**, as documented in `CISM runtime variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_. The model resolution, for example, is set via ``CISM_GRID``. The value of ``CISM_GRID`` determines the default value of a number of other namelist parameters.
