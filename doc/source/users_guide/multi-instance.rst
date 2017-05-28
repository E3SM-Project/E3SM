.. _multi-instance:

Multi-instance component functionality
======================================

The CIME coupling infrastructure has the capability to run multiple component instances under one model executable. 
The only caveat to this usage is that if N multiple instances of any one active component is used, then N multiple instances of ALL active components are required. 
More details are discussed below. 
The primary motivation for this development was to be able to run an ensemble Kalman-Filter for data assimilation and parameter estimation (e.g. UQ). 
However, it also provides you with the ability to run a set of experiments within a single model executable where each instance can have a different namelist, and have all the output go to one directory. 

In the following an F compset will be used as an illustration. Utilizing the multiple instance code involves the following steps:

1. create the case
::

   > create_newcase -case Fmulti -compset F -res ne30_g16 
   > cd Fmulti

2. Lets assume the following out of the box pe-layout 
::

   NTASKS(ATM)=128, NTHRDS(ATM)=1, ROOTPE(ATM)=0, NINST(ATM)=1
   NTASKS(LND)=128, NTHRDS(LND)=1, ROOTPE(LND)=0, NINST(LND)=1
   NTASKS(ICE)=128, NTHRDS(ICE)=1, ROOTPE(ICE)=0, NINST(ICE)=1
   NTASKS(OCN)=128, NTHRDS(OCN)=1, ROOTPE(OCN)=0, NINST(OCN)=1
   NTASKS(GLC)=128, NTHRDS(GLC)=1, ROOTPE(GLC)=0, NINST(GLC)=1
   NTASKS(WAV)=128, NTHRDS(WAV)=1, ROOTPE(WAV)=0, NINST(WAV)=1
   NTASKS(CPL)=128, NTHRDS(CPL)=1, ROOTPE(CPL)=0

In this F compset, the atm, lnd, rof are active components, the ocn is a prescribed data component, cice is a mixed prescribed/active component (ice-coverage is prescribed) and glc and wav are stub components. 
Lets say we want to run 2 instances of CAM in this experiment. 
The current implementation of multi-instances will also require you to run 2 instances of CLM, CICE and RTM. 
However, you have the flexibility to run either 1 or 2 instances of DOCN (we can ignore glc and wav since they do not do anything in this compset). 
To run 2 instances of CAM, CLM, CICE, RTM and DOCN, all you need to do is to invoke the following command in your ``$CASEROOT``:
::

   ./xmlchange NINST_ATM=2
   ./xmlchange NINST_LND=2
   ./xmlchange NINST_ICE=2
   ./xmlchange NINST_ROF=2
   ./xmlchange NINST_OCN=2

As a result of this, you will have 2 instances of CAM, CLM and CICE (prescribed), RTM, and DOCN,  each running concurrently on 64 MPI tasks  **TODO: put in reference to xmlchange".**

3. Setup the case
::

   > ./case.setup

New user_nl_xxx_NNNN file (where NNNN is the number of the component instances) will be generated when **case.setup** is called. 
In particular, calling **case.setup** with the above ``env_mach_pes.xml`` file will result in the following ``user_nl_*`` files in ``$CASEROOT``
::

   user_nl_cam_0001,  user_nl_cam_0002
   user_nl_cice_0001, user_nl_cice_0002
   user_nl_clm_0001,  user_nl_clm_0002
   user_nl_rtm_0001,  user_nl_rtm_0002
   user_nl_docn_0001, user_nl_docn_0002
   user_nl_cpl

and the following ``*_in_*`` files and ``*txt*`` files in $CASEROOT/CaseDocs:
::

   atm_in_0001, atm_in_0002
   docn.streams.txt.prescribed_0001, docn.streams.txt.prescribed_0002
   docn_in_0001, docn_in_0002
   docn_ocn_in_0001, docn_ocn_in_0002
   drv_flds_in, drv_in
   ice_in_0001, ice_in_0002
   lnd_in_0001, lnd_in_0002
   rof_in_0001, rof_in_0002

The namelist for each component instance can be modified by changing the corresponding user_nl_xxx_NNNN file for that component instance. 
Modifying the user_nl_cam_0002 will result in the namelist changes you put in to be active ONLY for instance 2 of CAM. 
To change the DOCN stream txt file instance 0002, you should place a copy of ``docn.streams.txt.prescribed_0002`` in ``$CASEROOT`` with the name ``user_docn.streams.txt.prescribed_0002`` and modify it accordlingly.

It is also important to stress the following points:

1. **Different component instances can ONLY differ by differences in namelist settings - they are ALL using the same model executable.**

2. Only 1 coupler component is supported currently in multiple instance implementation.

3. ``user_nl_xxx_NN`` files once they are created by **case.setup** *ARE NOT* removed by calling **case.setup -clean**. 

4. In general, you should run multiple instances concurrently (the default setting in ``env_mach_pes.xml``). 
   The serial setting is only for EXPERT USERS in upcoming development code implementations.
