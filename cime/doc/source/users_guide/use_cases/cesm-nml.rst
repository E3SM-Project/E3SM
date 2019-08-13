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

----
CISM
----
See `CISM namelist variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_ for a complete description of the CISM runtime namelist variables. This includes variables that appear both in **cism_in** and in **cism.config**.

To modify any of these settings, add the appropriate keyword/value pair at the end of the **user_nl_cism** file. (See the documentation for each file at the top of that file.) To see the result of your change, call **preview_namelists** and verify that the changes appear correctly in **CaseDocs/cism_in** and **CaseDocs/cism.config**.

Some CISM runtime settings are sets via **env_run.xml**, as documented in `CISM runtime variables <http://www.cesm.ucar.edu/models/cesm2.0/external-link-here>`_. The model resolution, for example, is set via ``CISM_GRID``. The value of ``CISM_GRID`` determines the default value of a number of other namelist parameters.
