!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_AttrVectComms - Communications methods for the AttrVect
!
! !DESCRIPTION:
!
! In this module, we define communications methods specific to the 
! {\tt AttrVect} class (see m_AttrVect for more information about this
! class and its methods).
!
! !INTERFACE:
 module m_AttrVectComms
!
! !USES:
!
      use m_AttrVect ! AttrVect class and its methods


      implicit none

      private	! except

      public :: gather		! gather all local vectors to the root
      public :: scatter		! scatter from the root to all PEs
      public :: bcast		! bcast from root to all PEs

    interface gather ; module procedure gather_ ; end interface
    interface scatter; module procedure scatter_; end interface
    interface bcast  ; module procedure bcast_  ; end interface

! !REVISION HISTORY:
! 	27Oct00 - J.W. Larson <larson@mcs.anl.gov> - relocated routines
!                 from m_AttrVect to create this module.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_AttrVectComms'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: gather_ - gather a vector according to a given _map_
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine gather_(iV,oV,Mp,root,comm,stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90
      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_lsize => lsize
      use m_GlobalMap, only : GlobalMap_gsize => gsize
      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect,  only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr

      implicit none

      type(AttrVect),intent(in)  :: iV
      type(AttrVect),intent(out) :: oV
      type(GlobalMap) ,intent(in)  :: Mp
      integer, intent(in) :: root
      integer, intent(in) :: comm
      integer, optional,intent(out) :: stat

! !REVISION HISTORY:
! 	15Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	27Oct00 - J.W. Larson <larson@mcs.anl.gov> - relocated from
!                 m_AttrVect
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::gather_'
  integer :: nIA,nRA,niV,noV,ier
  integer :: myID

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

	! Verify the input: a _scatterd_ vector

  niV=GlobalMap_lsize(Mp)
  noV=AttrVect_lsize(iV)

  if(niV /= noV) then
    write(stderr,'(2a,i4,a,i4)') myname_,	&
	': invalid input, lsize(Mp) =',niV,	&
	', lsize(iV) =',noV
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

  noV=GlobalMap_gsize(Mp) ! the gathered local size, as for the output
  if(myID /= root) nov=0
  call AttrVect_init(oV,iV,noV)

  niV=GlobalMap_lsize(Mp) ! the scattered local size, as for the input

  nIA=AttrVect_nIAttr(oV)	! number of INTEGER attributes
  nRA=AttrVect_nRAttr(oV)	! number of REAL attributes

  call MPI_gatherv(iV%iAttr,niV*nIA,MP_INTEGER,			&
	oV%iAttr,Mp%counts*nIA,Mp%displs*nIA,MP_INTEGER,	&
	root,comm,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_gatherv(iAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  call MPI_gatherv(iV%rAttr,niV*nRA,MP_REAL,		&
	oV%rAttr,Mp%counts*nRA,Mp%displs*nRA,MP_REAL,	&
	root,comm,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_gatherv(rAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

 end subroutine gather_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: scatter_ - scatter a vecter according to a given _map_
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine scatter_(iV,oV,Mp,root,comm,stat)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90
      use m_GlobalMap, only : GlobalMap
      use m_GlobalMap, only : GlobalMap_lsize => lsize
      use m_GlobalMap, only : GlobalMap_gsize => gsize
      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect,  only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr

      implicit none

      type(AttrVect),intent(in)  :: iV
      type(AttrVect),intent(out) :: oV
      type(GlobalMap) ,intent(in)  :: Mp
      integer, intent(in) :: root
      integer, intent(in) :: comm
      integer, optional,intent(out) :: stat

! !REVISION HISTORY:
! 	21Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	27Oct00 - J.W. Larson <larson@mcs.anl.gov> - relocated from
!                 m_AttrVect
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::scatter_'
  integer :: nIA,nRA,niV,noV,ier
  integer :: myID

  if(present(stat)) stat=0

  call MP_comm_rank(comm,myID,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MP_comm_rank()',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

	! Verify the input: a _gathered_ vector

  niV = GlobalMap_gsize(Mp)  ! the _gathered_ local size
  if(myID /= root) niV=0

  noV = AttrVect_lsize(iV)
  if(niV /= noV) then
    write(stderr,'(2a,i4,a,i4)') myname_,	&
	': invalid input, rsize(Mp) =',niV,	&
	', lsize(iV) =',noV
    if(.not.present(stat)) call die(myname_)
    stat=-1
    return
  endif

  noV = GlobalMap_lsize(Mp) ! the _scatterd_ local size
  call AttrVect_initv(oV,iV,noV)

  nIA = AttrVect_nIAttr(iV)	! number of INTEGER attributes
  nRA = AttrVect_nRAttr(iV)	! number of REAL attributes

  call MPI_scatterv(iV%iAttr(1,1),Mp%counts*nIA,	&
	Mp%displs*nIA,MP_INTEGER,			&
	oV%iAttr(1,1),noV*nRA,MP_INTEGER,root,comm,ier )
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_scatterv(iAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  call MPI_scatterv(iV%rAttr(1,1),Mp%counts*nRA,	&
	Mp%displs*nRA,MP_REAL,				&
	oV%rAttr(1,1),noV*nRA,MP_REAL,root,comm,ier )
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_scatterv(rAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

 end subroutine scatter_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bcast_ - broadcast from the root to all PEs
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine bcast_(aV,root,comm,stat)
!
! !USES:
!
      use m_die, only : die, perr
      use m_mpif90
      use m_String, only : String,bcast,char
      use m_String, only : String_bcast => bcast
      use m_List, only : List_get => get
      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr

      implicit none

      type(AttrVect) :: aV	! (IN) on the root, (OUT) elsewhere
      integer,intent(in) :: root
      integer,intent(in) :: comm
      integer,optional,intent(out) :: stat

! !REVISION HISTORY:
! 	27Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	27Oct00 - J.W. Larson <larson@mcs.anl.gov> - relocated from
!                 m_AttrVect
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::bcast_'
  type(String) :: iLStr,rLStr
  integer :: nIA, nRA, lsize
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

	! Convert the two Lists to two Strings

  if(myID == root)	&
    call List_get(iLStr,aV%iList)

  call String_bcast(iLStr,root,comm,stat=ier)	! bcast.String()
  if(ier /= 0) then
    call perr(myname_,'bcast.String(iLstr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  if(myID == root)	&
    call List_get(rLStr,aV%rList)

  call String_bcast(rLStr,root,comm,stat=ier)	! bcast.String()
  if(ier /= 0) then
    call perr(myname_,'bcast.String(rLstr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  if(myID == root) lsize = AttrVect_lsize(aV)

	! Set lsize for all PEs

  call MPI_bcast(lsize,1,MP_INTEGER,root,comm,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(lsize)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  if(myID /= root) then
     call AttrVect_init(aV,iList=char(iLStr),rList=char(rLStr), &
	                lsize=lsize)
  endif

  nIA = AttrVect_nIAttr(aV)
  nRA = AttrVect_nRAttr(aV)

  call MPI_bcast(aV%iAttr,nIA*lsize,MP_INTEGER,root,comm,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(iAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

  call MPI_bcast(aV%rAttr,nRA*lsize,MP_REAL,   root,comm,ier)
  if(ier /= 0) then
    call MP_perr(myname_,'MPI_bcast(rAttr)',ier)
    if(.not.present(stat)) call die(myname_)
    stat=ier
    return
  endif

 end subroutine bcast_

 end module m_AttrVectComms
