.. _multi-instance:

Multi-instance component functionality
======================================

The CIME coupling infrastructure is capable of running multiple component instances under one model executable. 
One caveat: If N multiple instances of any one active component are used, the same number of multiple instances of ALL active components are required. 
More details are discussed below.

The primary motivation for this development was to be able to run an ensemble Kalman-Filter for data assimilation and parameter estimation (UQ, for example). 
However, it also provides the ability to run a set of experiments within a single model executable where each instance can have a different namelist, and to have all the output go to one directory. 

An F compset is used in the following example. Using the multiple-instance code involves the following steps:

1. Create the case.
::

   > create_newcase --case Fmulti --compset F --res ne30_g16 
   > cd Fmulti

2. Assume this is the out-of-the-box pe-layout: 
::

   NTASKS(ATM)=128, NTHRDS(ATM)=1, ROOTPE(ATM)=0, NINST(ATM)=1
   NTASKS(LND)=128, NTHRDS(LND)=1, ROOTPE(LND)=0, NINST(LND)=1
   NTASKS(ICE)=128, NTHRDS(ICE)=1, ROOTPE(ICE)=0, NINST(ICE)=1
   NTASKS(OCN)=128, NTHRDS(OCN)=1, ROOTPE(OCN)=0, NINST(OCN)=1
   NTASKS(GLC)=128, NTHRDS(GLC)=1, ROOTPE(GLC)=0, NINST(GLC)=1
   NTASKS(WAV)=128, NTHRDS(WAV)=1, ROOTPE(WAV)=0, NINST(WAV)=1
   NTASKS(CPL)=128, NTHRDS(CPL)=1, ROOTPE(CPL)=0

The atm, lnd and rof are active components in this compset. The ocn is a prescribed data component, cice is a mixed prescribed/active component (ice-coverage is prescribed), and glc and wav are stub components.

Let's say we want to run two instances of CAM in this experiment. 
We will also have to run two instances of CLM, CICE and RTM. 
However, we can run either one or two instances of DOCN, and we can ignore glc and wav since they do not do anything in this compset as stub components.
 
To run two instances of CAM, CLM, CICE, RTM and DOCN, invoke the following commands in your **$CASEROOT** directory:
::

   > ./xmlchange NINST_ATM=2
   > ./xmlchange NINST_LND=2
   > ./xmlchange NINST_ICE=2
   > ./xmlchange NINST_ROF=2
   > ./xmlchange NINST_OCN=2

As a result, you will have two instances of CAM, CLM and CICE (prescribed), RTM, and DOCN, each running concurrently on 64 MPI tasks.

**TODO: put in reference to xmlchange".**

3. Set up the case
::

   > ./case.setup

A new **user_nl_xxx_NNNN** file (where NNNN is the number of the component instances) is generated when **case.setup** is called. 
When calling **case.setup** with the above ``env_mach_pes.xml`` file specifically, these files are created in **$CASEROOT**:
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

#. **Multiple component instances can differ ONLY in namelist settings; they ALL use the same model executable.**

#. Multiple-instance implementation supports only one coupler component.

#. Calling **case.setup** with ``--clean`` *DOES NOT* remove the **user_nl_xxx_NN** files created by **case.setup**.

#. Multiple instances generally should un concurrently, which is the default setting in **env_mach_pes.xml**. 
   The serial setting is only for EXPERT USERS in upcoming development code implementations.
