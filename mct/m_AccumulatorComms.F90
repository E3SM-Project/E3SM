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

    interface gather ; module procedure gather_ ; end interface
    interface scatter; module procedure scatter_; end interface
    interface bcast  ; module procedure bcast_  ; end interface

! !REVISION HISTORY:
! 	31Oct00 - Jay Larson <larson@mcs.anl.gov> - initial prototype--
!                 These routines were separated from the module 
!                 {\tt m\_Accumulator}
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_AccumulatorComms'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gather_ - gather a vector according to a given GlobalMap
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine gather_(iC,oC,Mp,root,comm,stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_lsize => lsize
      use m_GlobalMap, only : GlobalMap_gsize => gsize
      use m_AttrVect,  only : AttrVect_lsize => lsize
      use m_Accumulator,  only : Accumulator_lsize => lsize
      use m_AttrVectComms,  only : AttrVect_gather => gather
      use m_mpif90

      implicit none
      type(Accumulator),intent(in)  :: iC
      type(Accumulator),intent(out) :: oC
      type(GlobalMap) ,intent(in)   :: Mp
      integer, intent(in) :: root
      integer, intent(in) :: comm
      integer, optional,intent(out) :: stat

! !REVISION HISTORY:
! 	13Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
! 	31Oct00 - Jay Larson <larson@mcs.anl.gov> - relocated to the
!                 module m_AccumulatorComms
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gather_'

  integer :: nIA,nRA,niC,noC,ier
  integer :: myID
  integer :: sd_min, sd_max

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat = ier
    return
  endif

	! Verify the input: a _scatterd_ Accumulator.
        ! What we are checking here is that the GlobalMap
        ! size listed for iC is the same as the local size
        ! of iC

  niC = GlobalMap_lsize(Mp)
  noC = Accumulator_lsize(iC)
  if(niC /= noC) then
    write(stderr,'(2a,i4,a,i4)') myname_,	&
	': invalid input, lsize(Mp) =',niC,	&
	', lsize(iV) =',noC
    if(.not.present(stat)) call die(myname_)
    stat = -1
    return
  endif

        ! the gathered local size, for the output

  noC = GlobalMap_gsize(Mp)

        ! Check to be sure that all the accumulators have performed
        ! same number of time steps in the accumulation process:

  call MPI_REDUCE(iC%steps_done, sd_min, 1, MP_INTEGER, MP_MIN, &
                  root, comm, ier)

  call MPI_REDUCE(iC%steps_done, sd_max, 1, MP_INTEGER, MP_MAX, &
                  root, comm, ier)

  if(myID == root) then
     if(sd_min /= sd_max) then  ! i.e. timestepping mismatch
	write(stderr,'(2a,i4,a,i4)') myname_,	&
	': time accum. step mismatch, sd_min =',sd_min,	&
	', sd_max =',sd_max
	if(.not.present(stat)) call die(myname_)
	stat = -1
	return
     endif
  endif

        ! The Accumulator will be gathered onto the root processor,
        ! so oC will have zero storage space for non-root processors.

  if(myID /= root) noC = 0
  call initv_(oC,iC,noC,iC%num_steps,iC%steps_done) ! Create the output 
                                                    ! accumulator

        ! The only thing in the accumulator that needs to be
        ! gathered is its data, which are stored in its AttrVect
        ! component

  call AttrVect_gather(iC%av,oC%av,Mp,root,comm,ier)

  if(ier /= 0) then
    call MP_perr(myname_,'AttrVect_gather',ier)
    if(.not.present(stat)) call die(myname_)
    stat = ier
    return
  endif

 end subroutine gather_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: scatter_ - scatter an Accumulator using a GlobalMap.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine scatter_(iC,oC,Mp,root,comm,stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_GlobalMap, only : GlobalMap_lsize => lsize
      use m_GlobalMap, only : GlobalMap_gsize => gsize
      use m_Accumulator, only : Accumulator_lsize => lsize
      use m_mpif90
      use m_AttrVectComms, only : AttrVect_scatter => scatter

      implicit none
      type(Accumulator),intent(in)  :: iC
      type(Accumulator),intent(out) :: oC
      type(GlobalMap) ,intent(in)  :: Mp
      integer, intent(in) :: root
      integer, intent(in) :: comm
      integer, optional,intent(out) :: stat

! !REVISION HISTORY:
! 	14Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
! 	31Oct00 - Jay Larson <larson@mcs.anl.gov> - moved from the module
!                 m_Accumulator to m_AccumulatorComms
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::scatter_'

  integer :: niC,noC,ier
  integer :: num_steps, steps_done
  integer :: myID

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

	! Verify the input: a _gathered_ Accumulator.

  niC = GlobalMap_gsize(Mp)	! the _gathered_ local size
  if(myID /= root) niC=0

  noC = Accumulator_lsize(iC)
  if(niC /= noC) then
    write(stderr,'(2a,i4,a,i4)') myname_,	&
	': invalid input, rsize(Mp) =',niC,	&
	', lsize(iC) =',noC
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

  noC = GlobalMap_lsize(Mp)	! the _scatterd_ local size

 
	! Broadcast the root value of num_steps

  call MPI_BCAST(num_steps, 1, MP_INTEGER, root, comm, ier)

  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(num_steps)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

	! Broadcast the root value of steps_done

  call MPI_BCAST(steps_done, 1, MP_INTEGER, root, comm, ier)

  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(steps_done)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

	! Initialize the Accumulator to which data will be
        ! scattered:

  call initv_(oC,iC,noC,num_steps,steps_done)

	! Scatter the AttrVect component av using the GlobalMap Mp:

  call AttrVect_scatter(iC%av, oC%av, Mp, root, comm, ier)

 end subroutine scatter_



!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bcast_ - broadcast an Accumulator from the root to all PEs.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine bcast_(aC, root, comm, stat)
!
! !USES:
!
      use m_die, only : die, perr
      use m_mpif90
      use m_AttrVectComms, only : AttrVect_bcast => bcast

      implicit none

      type(Accumulator)  :: aC	! (IN) on the root, (OUT) elsewhere
      integer,intent(in) :: root
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	14Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
! 	31Oct00 - Jay Larson <larson@mcs.anl.gov> - moved from the module
!                 m_Accumulator to m_AccumulatorComms
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
