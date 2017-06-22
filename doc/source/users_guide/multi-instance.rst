.. _multi-instance:

Multi-instance component functionality
======================================

The CIME coupling infrastructure is capable of running multiple component instances (ensembles) under one model executable.  There are two modes of ensemble capability, single coupler in which all component instances are handled by a single coupler component or multi-coupler in which each instance includes a separate coupler component.  In the multi-coupler mode the entire model is duplicated for each instance while in the single coupler mode only active components need be duplicated.   In most cases the multi-coupler mode will give better performance and should be used.

The primary motivation for this development was to be able to run an ensemble Kalman-Filter for data assimilation and parameter estimation (UQ, for example).
However, it also provides the ability to run a set of experiments within a single model executable where each instance can have a different namelist, and to have all the output go to one directory.

An F compset is used in the following example. Using the multiple-instance code involves the following steps:

1. Create the case.
::

   > create_newcase --case Fmulti --compset F --res ne30_ne30_mg17
   > cd Fmulti

2. Assume this is the out-of-the-box pe-layout:
::

   Comp  NTASKS  NTHRDS  ROOTPE
   CPL :    144/     1;      0
   ATM :    144/     1;      0
   LND :    144/     1;      0
   ICE :    144/     1;      0
   OCN :    144/     1;      0
   ROF :    144/     1;      0
   GLC :    144/     1;      0
   WAV :    144/     1;      0
   ESP :      1/     1;      0

The atm, lnd and rof are active components in this compset. The ocn is a prescribed data component, cice is a mixed prescribed/active component (ice-coverage is prescribed), and glc, wav and esp are stub components.

Let's say we want to run two instances of CAM in this experiment.
We will also have to run two instances of CLM, CICE and RTM.
However, we can run either one or two instances of DOCN, and we can ignore the stub components since they do not do anything in this compset.

To run two instances of CAM, CLM, CICE, RTM and DOCN, invoke the following :ref: `xmlchange<modifying-an-xml-file>` commands in your **$CASEROOT** directory:
::

   > ./xmlchange NINST_ATM=2
   > ./xmlchange NINST_LND=2
   > ./xmlchange NINST_ICE=2
   > ./xmlchange NINST_ROF=2
   > ./xmlchange NINST_OCN=2

As a result, you will have two instances of CAM, CLM and CICE (prescribed), RTM, and DOCN, each running concurrently on 72 MPI tasks and all using the same coupler component.   In this single coupler mode the number of tasks for each component instance is NTASKS_COMPONENT/NINST_COMPONENT and the total number of tasks is the same as for the single instance case.

Now consider the multi coupler model.
To use this mode change the NINST values for the individual components back to 1 and the NINST_CPL to 2.
::
   > ./xmlchange NINST=1
   > ./xmlchange NINST_CPL=2

This configuration will run each component instance on the original 144 tasks but will generate two copies of the model (in the same executable) for a total of 288 tasks.

3. Set up the case
::

   > ./case.setup

A new **user_nl_xxx_NNNN** file (where NNNN is the number of the component instances) is generated when **case.setup** is called.
When calling **case.setup** with the **env_mach_pes.xml** file specifically, these files are created in **$CASEROOT**:
::

   user_nl_cam_0001,  user_nl_cam_0002
   user_nl_cice_0001, user_nl_cice_0002
   user_nl_clm_0001,  user_nl_clm_0002
   user_nl_rtm_0001,  user_nl_rtm_0002
   user_nl_docn_0001, user_nl_docn_0002
   user_nl_cpl

Also, **case.setup** creates the following ``*_in_*`` files and ``*txt*`` files in **$CASEROOT/CaseDocs**:
::

   atm_in_0001, atm_in_0002
   docn.streams.txt.prescribed_0001, docn.streams.txt.prescribed_0002
   docn_in_0001, docn_in_0002
   docn_ocn_in_0001, docn_ocn_in_0002
   drv_flds_in, drv_in
   ice_in_0001, ice_in_0002
   lnd_in_0001, lnd_in_0002
   rof_in_0001, rof_in_0002

The namelist for each component instance can be modified by changing the corresponding **user_nl_xxx_NNNN** file.
Modifying **user_nl_cam_0002** will result in your namelist changes being active ONLY for the second instance of CAM.
To change the DOCN stream txt file instance 0002, copy **docn.streams.txt.prescribed_0002** to your **$CASEROOT** directory with the name **user_docn.streams.txt.prescribed_0002** and modify it accordlingly.

Also keep these important points in mind:

#. Note that these changes can be made at create_newcase time with option --ninst # where # is a positive integer, use the additional logical option --ninst-couplers to invoke the multi-coupler mode.

#. **Multiple component instances can differ ONLY in namelist settings; they ALL use the same model executable.**

#. Calling **case.setup** with ``--clean`` *DOES NOT* remove the **user_nl_xxx_NN** (where xxx is the component name) files created by **case.setup**.

#. A special variable NINST_LAYOUT is provided for some experimental compsets, its value should be
   'concurrent' for all but a few special cases.

#. In **create_test** these options can be invoked with testname modifiers _N# for the single coupler mode and _C# for the multi-coupler mode.  These are mutually exclusive options, they cannot be combined.
