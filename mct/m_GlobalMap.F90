!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$ 
!BOP -------------------------------------------------------------------
!
! !MODULE: m_GlobalMap - a size description of a distributed array
!
! !DESCRIPTION:
!
! !INTERFACE:

 module m_GlobalMap
      implicit none
      private	! except

      public :: GlobalMap		! The class data structure
      public :: gsize
      public :: lsize
      public :: init
      public :: init_remote
      public :: clean
      public :: rank
      public :: bounds
      public :: comp_id

    type GlobalMap
      integer :: comp_id                        ! Component ID number
      integer :: gsize				! the Global size
      integer :: lsize				! my local size
      integer,dimension(:),pointer :: counts	! all local sizes
      integer,dimension(:),pointer :: displs	! PE ordered locations
    end type GlobalMap

    interface gsize; module procedure gsize_; end interface
    interface lsize; module procedure lsize_; end interface
    interface init ; module procedure	&
	initd_,	&	! initialize from all PEs
	initr_		! initialize from the root
    end interface
    interface init_remote; module procedure init_remote_; end interface
    interface clean; module procedure clean_; end interface
    interface rank ; module procedure rank_ ; end interface
    interface bounds; module procedure bounds_; end interface
    interface comp_id ; module procedure comp_id_ ; end interface

! !REVISION HISTORY:
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	09Nov00 - J.W. Larson <larson@mcs.anl.gov> - added init_remote
!                 interface.
!       26Jan01 - J.W. Larson <larson@mcs.anl.gov> - added storage for
!                 component ID number GlobalMap%comp_id, and associated
!                 method comp_id_()
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_GlobalMap'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initd_ - define the map from distributed data
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initd_(GMap, comp_id, ln, comm)
      use m_mpif90
      use m_die
      implicit none
      type(GlobalMap),intent(out) :: GMap
      integer,intent(in) :: comp_id ! Component
      integer,intent(in) :: ln	    ! the local size
      integer,intent(in) :: comm    ! f90 MPI communicator handle 

! !REVISION HISTORY:
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initd_'
  integer :: nPEs,myID,ier,l,i

  call MP_comm_size(comm,nPEs,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_size()',ier)

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)

  allocate(GMap%counts(0:nPEs-1),GMap%displs(0:nPEs-1),stat=ier)
  if(ier /= 0) call die(myname_,'allocate()',ier)

#ifdef MALL_ON
	call mall_ci(size(transfer(GMap%counts,(/1/))),myname_)
	call mall_ci(size(transfer(GMap%displs,(/1/))),myname_)
#endif

  call MPI_allgather(ln,1,MP_INTEGER,GMap%counts,1,MP_INTEGER,comm,ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_allgather()',ier)

  l=0
  do i=0,nPEs-1
    GMap%displs(i)=l
    l=l+GMap%counts(i)
  end do

  GMap%lsize=GMap%counts(myID)	! the local size
  GMap%gsize=l	! the global size
  GMap%comp_id = comp_id ! the component ID number

 end subroutine initd_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initr_ initialize the map from the root
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine initr_(GMap, comp_id, lns, root, comm)
      use m_mpif90
      use m_die
      use m_stdio
      implicit none
      type(GlobalMap),intent(out) :: GMap
      integer            :: comp_id             ! component ID number
      integer,dimension(:),intent(in) :: lns	! the distributed sizes
      integer,intent(in) :: root
      integer,intent(in) :: comm

! !REVISION HISTORY:
! 	29May98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initr_'
  integer :: nPEs,myID,ier,l,i

  call MP_comm_size(comm,nPEs,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_size()',ier)

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)

  allocate(GMap%counts(0:nPEs-1),GMap%displs(0:nPEs-1),stat=ier)
  if(ier /= 0) call die(myname_,'allocate()',ier)

#ifdef MALL_ON
	call mall_ci(size(transfer(GMap%counts,(/1/))),myname_)
	call mall_ci(size(transfer(GMap%displs,(/1/))),myname_)
#endif

  if(myID == root) then
    if(size(lns(:)) /= nPEs) then
      write(stderr,'(2a,2(a,i4))') myname_,	&
	': _root_ argument error',		&
	', size(lns) =',size(lns),		&
	', nPEs =',nPEs
      call die(myname_)
    endif

    GMap%counts(:)=lns(:)
  endif

  call MPI_bcast(GMap%counts,nPEs,MP_INTEGER,root,comm,ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_bcast()',ier)

  ! on each process, use GMap%counts(:) to compute GMap%displs(:)

  l=0
  do i=0,nPEs-1
    GMap%displs(i)=l
    l=l+GMap%counts(i)
  end do

  GMap%lsize=GMap%counts(myID)	! the local size
  GMap%gsize=l	! the global size

  ! finally, set and broadcast the component ID number GMap%comp_id

  if(myID == root) GMap%comp_id = comp_id

  call MPI_bcast(GMap%comp_id,1,MP_INTEGER,root,comm,ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_bcast()',ier)

 end subroutine initr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_remote_ initialize remote GlobalMap from the root
!
! !DESCRIPTION:
! This method takes data describing a GlobalMap for another communicator
!
! !INTERFACE:

 subroutine init_remote_(GMap, remote_lns, remote_npes, my_root, &
	                    my_comm, remote_comp_id)
      use m_mpif90
      use m_die
      use m_stdio
      implicit none
      type(GlobalMap),intent(out) :: GMap
      integer, dimension(:) :: remote_lns     ! distributed sizes on 
                                              ! the remote communicator
      integer            :: remote_npes       ! number of processes
                                              ! on remote communicator
      integer,intent(in) :: my_root           ! my root
      integer,intent(in) :: my_comm           ! my communicator
      integer,intent(in) :: remote_comp_id    ! remote component ID nubmer

! !REVISION HISTORY:
! 	08Nov00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
! 	26Jan01 - J.W. Larson <larson@mcs.anl.gov> - slight change--remote
!                 communicator is replaced by remote component ID number
!                 in argument remote_comp_id.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_remote_'
  integer :: nPEs,myID,ier,l,i


        ! Which processor am I on communicator my_comm?  Store
        ! the answer in myID:

  call MP_comm_rank(my_comm, myID, ier)
  if(ier /= 0) call MP_perr_die(myname_,'MP_comm_rank()',ier)

        ! allocate counts and displacements component arrays
        ! for the sake of compactness, store the value of remote_npes
        ! in the more tersely named variable nPEs.

  nPEs = remote_npes

  allocate(GMap%counts(0:nPEs-1),GMap%displs(0:nPEs-1),stat=ier)
  if(ier /= 0) call die(myname_,'allocate()',ier)

#ifdef MALL_ON
	call mall_ci(size(transfer(GMap%counts,(/1/))),myname_)
	call mall_ci(size(transfer(GMap%displs,(/1/))),myname_)
#endif

        ! On the Root processor, check the size of remote_lns(:)
        ! to see it is equal to nPEs, the number of remote processes,
        ! then store it as GMap%counts and broadcast it.

  if(myID == my_root) then
    if(size(remote_lns(:)) /= nPEs) then
      write(stderr,'(2a,2(a,i4))') myname_,	 &
	': _root_ argument error',		 &
	', size(remote_lns) =',size(remote_lns), &
	', nPEs =',nPEs
      call die(myname_)
    endif

    GMap%counts(:)=remote_lns(:)
  endif

  call MPI_bcast(GMap%counts, nPEs, MP_INTEGER, my_root, my_comm, ier)
  if(ier/=0) call MP_perr_die(myname_,'MPI_bcast()',ier)

        ! Now, on each processor of my_comm, compute from 
        ! GMap%counts(:) the entries of GMap%displs(:)

  l=0
  do i=0,nPEs-1
    GMap%displs(i)=l
    l=l+GMap%counts(i)
  end do

  GMap%lsize = GMap%counts(myID) ! the local size
  GMap%gsize = l      	         ! the global size
  GMap%comp_id = remote_comp_id  ! the remote component id

 end subroutine init_remote_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean the map
!
! !DESCRIPTION:
!
! !INTERFACE:

  subroutine clean_(GMap,stat)
      use m_die
      implicit none
      type(GlobalMap),intent(inout) :: GMap
      integer,optional,intent(out)  :: stat

! !REVISION HISTORY:
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	26Jan01 - J. Larson <larson@mcs.anl.gov> incorporated comp_id.
!       01Mar02 - E.T. Ong <eong@mcs.anl.gov> removed the die to prevent
!                 crashes and added stat argument.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  deallocate(GMap%counts,GMap%displs,stat=ier)

  if(present(stat)) then
     stat=ier
  else
     if(ier /= 0) call warn(myname_,'deallocate(GMap%...)',ier)
  endif
  
  if(ier == 0) then

#ifdef MALL_ON
	call mall_co(size(transfer(GMap%counts,(/1/))),myname_)
	call mall_co(size(transfer(GMap%displs,(/1/))),myname_)
#endif

  endif

  GMap%lsize = 0
  GMap%gsize = 0
  GMap%comp_id = 0

 end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - find the local size from the map
!
! !DESCRIPTION:
!
! !INTERFACE:

 function lsize_(GMap)
      implicit none
      type(GlobalMap),intent(in) :: GMap
      integer :: lsize_

! !REVISION HISTORY:
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'

  lsize_=GMap%lsize

 end function lsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gsize_ - find the global size from the map
!
! !DESCRIPTION:
!
! !INTERFACE:

 function gsize_(GMap)
      implicit none
      type(GlobalMap),intent(in) :: GMap
      integer :: gsize_

! !REVISION HISTORY:
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gsize_'

  gsize_=GMap%gsize

 end function gsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: rank_ - rank (which PE) of a given global index
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine rank_(GMap,i_g,rank)
      implicit none
      type(GlobalMap),intent(in) :: GMap
      integer, intent(in)  :: i_g	! a global index
      integer, intent(out) :: rank

! !REVISION HISTORY:
! 	05May98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::rank_'
  integer :: i,ilc,ile

  rank=-1	! if nowhere fits
  do i=0,size(GMap%displs)-1
    ilc=GMap%displs(i)
    ile=ilc+GMap%counts(i)

		! If i_g in (ilc,ile].  Note that i_g := [1:..]

    if(ilc < i_g .and. i_g <= ile) then
      rank=i
      return
    endif
  end do

 end subroutine rank_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bounds_ - upper/lower indices for a given process id.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine bounds_(GMap, pe_no, lbnd, ubnd)
      implicit none
      type(GlobalMap),intent(in) :: GMap
      integer,intent(in)  :: pe_no
      integer,intent(out) :: lbnd
      integer,intent(out) :: ubnd

! !REVISION HISTORY:
! 	30Jan01 - J. Larson <larson@mcs.anl.gov> - initial code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::bounds_'

  lbnd = GMap%displs(pe_no) + 1
  ubnd = lbnd + GMap%counts(pe_no) - 1

 end subroutine bounds_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: comp_id_ - the component id number in the GlobalMap.
!
! !DESCRIPTION:
!
! !INTERFACE:

    integer function comp_id_(GMap)
      implicit none
      type(GlobalMap),intent(in) :: GMap

! !REVISION HISTORY:
! 	25Jan02 - J. Larson <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::comp_id_'

  comp_id_ = GMap%comp_id

 end function comp_id_

 end module m_GlobalMap
