!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_AccumulatorComms - Communication methods for the 
!          the Accumulator class.
!
! !DESCRIPTION:
!
! An {\em accumulator} is a data class used for computing running sums 
! and/or time averages of {\tt AttrVect} class data (see 
! {\tt m\_Accumulator} for details).  This module defines the 
! communications methods for the accumulator, employing both the 
! {\tt GlobalMap} and {\tt GlobalSegMap} decomposition descriptors.
!
! !INTERFACE:

 module m_AccumulatorComms
!
! !USES:
!
      use m_Accumulator, only : Accumulator
      use m_GlobalMap,   only : GlobalMap

      implicit none

      private	! except

! List of communications Methods for the Accumulator class

      public :: gather		! gather all local vectors to the root
      public :: scatter		! scatter from the root to all PEs
      public :: bcast		! bcast from root to all PEs

! Definition of interfaces for the communication methods for 
! the Accumulator:

    interface gather ; module procedure &
              GM_gather_, &
              GSM_gather_ 
    end interface
    interface scatter ; module procedure &
              GM_scatter_, &
              GSM_scatter_ 
    end interface
    interface bcast  ; module procedure bcast_  ; end interface

! !REVISION HISTORY:
! 	31Oct00 - Jay Larson <larson@mcs.anl.gov> - initial prototype--
!                 These routines were separated from the module 
!                 {\tt m\_Accumulator}
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - Specification of 
!                 APIs for the routines {\tt GSM_gather_()} and 
!                 {\tt GSM_scatter_()}.
!       10May01 - Jay Larson <larson@mcs.anl.gov> - Changes in the
!                 comms routine to match the MPI model for collective
!                 communications, and general clean-up of prologues.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_AccumulatorComms'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GM_gather_ - gather a vector using input GlobalMap.
!
! !DESCRIPTION:  {\tt GM\_gather()} takes a distributed (across the
! communicator associated with the handle {\tt comm}) input 
! {\tt Accumulator} argument {\tt iC} and gathers its data to the
! {\tt Accumulator} {\tt oC} on the {\tt root}.  The decomposition of 
! {\tt iC} is described by the input {\tt GlobalMap} argument {\tt Gmap}.
! The success (failure) of this operation is signified by the zero (nonzero) 
! value of the optional output argument {\tt stat}.
!
! !INTERFACE:

 subroutine GM_gather_(iC, oC, GMap, root, comm, stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90

      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_lsize => lsize
      use m_GlobalMap, only : GlobalMap_gsize => gsize

      use m_AttrVect,  only : AttrVect_lsize => lsize

      use m_Accumulator,   only : Accumulator_lsize => lsize
      use m_AttrVectComms, only : AttrVect_gather => gather

      implicit none

! !INPUT PARAMETERS: 
!
      type(Accumulator), intent(in)  :: iC
      type(GlobalMap) ,  intent(in)  :: GMap
      integer,           intent(in)  :: root
      integer,           intent(in)  :: comm

! !OUTPUT PARAMETERS: 
!
      type(Accumulator), intent(out) :: oC
      integer, optional,intent(out)  :: stat

! !REVISION HISTORY:
! 	13Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
! 	31Oct00 - Jay Larson <larson@mcs.anl.gov> - relocated to the
!                 module m_AccumulatorComms
! 	15Jan01 - Jay Larson <larson@mcs.anl.gov> - renamed GM_gather_
! 	10May01 - Jay Larson <larson@mcs.anl.gov> - revamped comms 
!                 model to match MPI comms model, and cleaned up prologue
!EOP ___________________________________________________________________

 character(len=*),parameter :: myname_=myname//'::GM_gather_'
 integer :: myID, ier

        ! Initialize status flag (if present)

  if(present(stat)) stat=0

  call MP_comm_rank(comm, myID, ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

       ! On the root, set oC%num_steps and oC%steps_done

  if(myID == root) then
     oC%num_steps = iC%num_steps
     oC%steps_done = iC%steps_done
  endif

       ! Gather distributed iC%av to oC%av on the root

  call AttrVect_gather(iC%av, oC%av, GMap, root, comm, ier)

  if(ier /= 0) then
    call MP_perr(myname_,'AttrVect_gather(iC%av, oC%av...',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

 end subroutine GM_gather_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GSM_gather_ - gather an Accumulator using input GlobalSegMap.
!
! !DESCRIPTION:  This routine takes the distrubuted (on the communcator
! associated with the handle {\tt comm}) input {\tt Accumulator} 
! argument {\tt iC} gathers it to the the {\tt Accumulator} argument 
! {\tt oC} (valid only on the {\tt root}).  The decompositon of {\tt iC}
! is contained in the input {\tt GlobalSegMap} argument {\tt GSMap}.  
! The success (failure) of this operation is signified by the zero 
! (nonzero) returned value of the {\tt INTEGER} flag {\tt stat}.
!
! !INTERFACE:

 subroutine GSM_gather_(iC, oC, GSMap, root, comm, stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90
      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_lsize => lsize
      use m_GlobalSegMap, only : GlobalSegMap_gsize => gsize

      use m_AttrVect,  only : AttrVect_lsize => lsize

      use m_Accumulator,   only : Accumulator_lsize => lsize

      use m_AttrVectComms, only : AttrVect_gather => gather

      implicit none

! !INPUT PARAMETERS: 
!
      type(Accumulator),  intent(in) :: iC
      type(GlobalSegMap), intent(in) :: GSMap
      integer,            intent(in) :: root
      integer,            intent(in) :: comm

! !OUTPUT PARAMETERS: 
!
      type(Accumulator), intent(out) :: oC
      integer, optional, intent(out) :: stat

! !REVISION HISTORY:
! 	15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
! 	10May01 - Jay Larson <larson@mcs.anl.gov> - Initial code and
!                 cleaned up prologue.
!EOP ___________________________________________________________________

 character(len=*),parameter :: myname_=myname//'::GSM_gather_'
 integer :: myID, ier

        ! Initialize status flag (if present)

  if(present(stat)) stat=0

  call MP_comm_rank(comm, myID, ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

       ! On the root, set oC%num_steps and oC%steps_done

  if(myID == root) then
     oC%num_steps = iC%num_steps
     oC%steps_done = iC%steps_done
  endif

       ! Gather distributed iC%av to oC%av on the root

  call AttrVect_gather(iC%av, oC%av, GSMap, root, comm, ier)

  if(ier /= 0) then
    call MP_perr(myname_,'AttrVect_gather(iC%av, oC%av...',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

 end subroutine GSM_gather_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GM_scatter_ - scatter an Accumulator using input GlobalMap.
!
! !DESCRIPTION:  This routine takes the input {\tt Accumulator} argument
! {\tt iC} (valid only on the {\tt root}), and scatters it to the 
! distributed {\tt Accumulator} argument {\tt oC} on the processes 
! associated with the communicator handle {\tt comm}.  The decompositon
! used to scatter the data is contained in the input {\tt GlobalMap} 
! argument {\tt GMap}.  The success (failure) of this operation is 
! signified by the zero (nonzero) returned value of the {\tt INTEGER} 
! flag {\tt stat}.
!
! !INTERFACE:

 subroutine GM_scatter_(iC, oC, GMap, root, comm, stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90

      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_lsize => lsize
      use m_GlobalMap, only : GlobalMap_gsize => gsize

      use m_Accumulator, only : Accumulator_lsize => lsize

      use m_AttrVectComms, only : AttrVect_scatter => scatter

      implicit none

! !INPUT PARAMETERS: 
!
      type(Accumulator), intent(in)  :: iC
      type(GlobalMap),   intent(in)  :: GMap
      integer,           intent(in)  :: root
      integer,           intent(in)  :: comm

! !OUTPUT PARAMETERS: 
!
      type(Accumulator), intent(out) :: oC
      integer, optional, intent(out) :: stat

! !REVISION HISTORY:
! 	14Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
! 	31Oct00 - Jay Larson <larson@mcs.anl.gov> - moved from the module
!                 m_Accumulator to m_AccumulatorComms
! 	15Jan01 - Jay Larson <larson@mcs.anl.gov> - renamed GM_scatter_.
!       10May01 - Jay Larson <larson@mcs.anl.gov> - revamped code to fit
!                 MPI-like comms model, and cleaned up prologue.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GM_scatter_'

  integer :: myID, ier
  integer :: StepBuff(2)

        ! Initialize status flag (if present)

  if(present(stat)) stat=0

  call MP_comm_rank(comm, myID, ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

        ! On the root, load up iC%num_steps and iC%steps_done

  StepBuff(1) = iC%num_steps
  StepBuff(2) = iC%steps_done

	! Broadcast the root value of StepBuff

  call MPI_BCAST(StepBuff, 2, MP_INTEGER, root, comm, ier)

  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(StepBuff...',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

        ! On all processes  unload iC%num_steps and 
        ! iC%steps_done from StepBuff

  oC%num_steps  = StepBuff(1)
  oC%steps_done = StepBuff(1)

	! Scatter the AttrVect component of aC

  call AttrVect_scatter(iC%av, oC%av, GMap, root, comm, ier)

  if(ier /= 0) then
    call MP_perr(myname_,'AttrVect_scatter(iC%av, oC%av...',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

 end subroutine GM_scatter_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: GSM_scatter_ - scatter Accumulator using input GlobalSegMap.
!
! !DESCRIPTION:  This routine takes the input {\tt Accumulator} argument
! {\tt iC} (valid only on the {\tt root}), and scatters it to the 
! distributed {\tt Accumulator} argument {\tt oC} on the processes 
! associated with the communicator handle {\tt comm}.  The decompositon
! used to scatter the data is contained in the input {\tt GlobalSegMap} 
! argument {\tt GSMap}.  The success (failure) of this operation is 
! signified by the zero (nonzero) returned value of the {\tt INTEGER} 
! flag {\tt stat}.
!
! !INTERFACE:

 subroutine GSM_scatter_(iC, oC, GSMap, root, comm, stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_lsize => lsize
      use m_GlobalSegMap, only : GlobalSegMap_gsize => gsize

      use m_Accumulator,  only : Accumulator_lsize => lsize

      use m_AttrVectComms, only : AttrVect_scatter => scatter

      implicit none

! !INPUT PARAMETERS: 
!
      type(Accumulator),  intent(in)  :: iC
      type(GlobalSegMap), intent(in)  :: GSMap
      integer,            intent(in)  :: root
      integer,            intent(in)  :: comm

! !OUTPUT PARAMETERS: 
!
      type(Accumulator),  intent(out) :: oC
      integer, optional,  intent(out) :: stat

! !REVISION HISTORY:
! 	15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
! 	10May01 - Jay Larson <larson@mcs.anl.gov> - Initial code/prologue
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::GSM_scatter_'

  integer :: myID, ier
  integer :: StepBuff(2)

        ! Initialize status flag (if present)

  if(present(stat)) stat=0

  call MP_comm_rank(comm, myID, ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

        ! On the root, load up iC%num_steps and iC%steps_done

  StepBuff(1) = iC%num_steps
  StepBuff(2) = iC%steps_done

	! Broadcast the root value of StepBuff

  call MPI_BCAST(StepBuff, 2, MP_INTEGER, root, comm, ier)

  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(StepBuff...',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

        ! On all processes  unload iC%num_steps and 
        ! iC%steps_done from StepBuff

  oC%num_steps  = StepBuff(1)
  oC%steps_done = StepBuff(1)

	! Scatter the AttrVect component of aC

  call AttrVect_scatter(iC%av, oC%av, GSMap, root, comm, ier)

  if(ier /= 0) then
    call MP_perr(myname_,'AttrVect_scatter(iC%av, oC%av...',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

 end subroutine GSM_scatter_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bcast_ - broadcast an Accumulator from the root to all PEs.
!
! !DESCRIPTION:  This routine takes the input {\tt Accumulator} argument
! {\tt aC} (on input valid only on the {\tt root}), and broadcasts it 
! to all the processes associated with the communicator handle 
! {\tt comm}.  The success (failure) of this operation is signified by 
! the zero (nonzero) returned value of the {\tt INTEGER} flag {\tt stat}.
!
! !INTERFACE:
!
 subroutine bcast_(aC, root, comm, stat)

!
! !USES:
!
      use m_die, only : die, perr
      use m_mpif90
      use m_AttrVectComms, only : AttrVect_bcast => bcast

      implicit none

! !INPUT PARAMETERS: 
!
      integer,intent(in) :: root
      integer,intent(in) :: comm

! !INPUT/OUTPUT PARAMETERS: 
!
      type(Accumulator), intent(inout) :: aC ! (IN) on root, (OUT) elsewhere

! !OUTPUT PARAMETERS: 
!
      integer, optional, intent(out)   :: stat

! !REVISION HISTORY:
! 	14Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
! 	31Oct00 - Jay Larson <larson@mcs.anl.gov> - moved from the module
!                 m_Accumulator to m_AccumulatorComms
!       09May01 - Jay Larson <larson@mcs.anl.gov> - cleaned up prologue
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::bcast_'

  integer :: myID
  integer :: ier

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

	! Broadcast the root value of num_steps

  call MPI_BCAST(aC%num_steps, 1, MP_INTEGER, root, comm, ier)

  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(aC%num_steps)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

	! Broadcast the root value of steps_done

  call MPI_BCAST(aC%steps_done, 1, MP_INTEGER, root, comm, ier)

  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(aC%steps_done)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

	! Broadcast the root value of aC%av

  call AttrVect_bcast(aC%av, root, comm, ier)

  if(ier /= 0) then
    call MP_perr(myname_,'AttrVect_bcast(aC%av)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

 end subroutine bcast_
 
 end module m_AccumulatorComms
