!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Navigator - Array of pointers
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_Navigator
      implicit none
      private	! except

      public :: Navigator		! The class data structure
      public :: Navigator_init,init	! initialize an object
      public :: clean			! clean an object
      public :: lsize			! the true size
      public :: msize			! the maximum size
      public :: resize			! adjust the true size
      public :: get			! get an entry

      public :: ptr_displs		! referencing %displs(:)
      public :: ptr_counts		! referencing %counts(:)

    type Navigator
      integer :: lsize	! true size, if over-dimensioned
      integer,pointer,dimension(:) :: displs
      integer,pointer,dimension(:) :: counts
    end type Navigator

    interface Navigator_init; module procedure	&
	init_; end interface
    interface init  ; module procedure init_  ; end interface
    interface clean ; module procedure clean_ ; end interface
    interface lsize ; module procedure lsize_ ; end interface
    interface msize ; module procedure msize_ ; end interface
    interface resize; module procedure resize_; end interface
    interface get   ; module procedure get_   ; end interface
    interface ptr_displs; module procedure	&
	ptr_displs_; end interface
    interface ptr_counts; module procedure	&
	ptr_counts_; end interface

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
! 	21Oct00	- J.W. Larson <larson@mcs.anl.gov>
!                 minor modification (removal of 'private' 
!                 statement in declaration of Navigator.  This
!                 will be put back as the MCT matures.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Navigator'

contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize a navigator
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine init_(nav,lsize,stat)
      use m_mall,only : mall_ison,mall_mci
      use m_die ,only : die,perr
      implicit none
      type(Navigator),intent(out) :: nav	! the object
      integer,intent(in) :: lsize		! nominal size
      integer,optional,intent(out) :: stat	! status

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: ier

  allocate(nav%displs(lsize),nav%counts(lsize),stat=ier)
	if(ier/=0) then
	  call perr(myname_,'allocate()',ier)
	  if(.not.present(stat)) call die(myname_)
	  stat=ier
	  return
	endif
	if(mall_ison()) then
	  call mall_mci(nav%displs,myname)
	  call mall_mci(nav%counts,myname)
	endif

  nav%lsize=lsize
end subroutine init_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean a navigator
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(nav,stat)
      use m_mall,only : mall_ison,mall_mco
      use m_die ,only : die,perr
      implicit none
      type(Navigator),intent(inout) :: nav	! the object
      integer,optional,intent(out) :: stat	! status

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  deallocate(nav%displs,nav%counts,stat=ier)

  if(present(stat)) then
     stat=ier
  else
     if(ier /= 0) call warn(myname_,'deallocate(nav%...)',ier)
  endif

  if(ier == 0) then

     if(mall_ison()) then
	call mall_mco(nav%displs,myname)
	call mall_mco(nav%counts,myname)
     endif

  endif

  nav%lsize=0
end subroutine clean_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - return the true size
!
! !DESCRIPTION:
!
! !INTERFACE:

    function lsize_(nav)
      implicit none
      type(Navigator),intent(in) :: nav
      integer :: lsize_

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!       01Mar02 - E.T. Ong <eong@mcs.anl.gov> - removed die to prevent 
!                 crashes.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'

  lsize_=nav%lsize
end function lsize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: msize_ - return the maximum size
!
! !DESCRIPTION:
!
! !INTERFACE:

    function msize_(nav)
      implicit none
      type(Navigator),intent(in) :: nav
      integer :: msize_

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::msize_'

  msize_=size(nav%displs)
end function msize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: resize_ - adjust the true size
!
! !DESCRIPTION:
!
!       The user is responsibile for ensuring lsize is no greater than
!   the maximum size of the vector.
!
!       If lsize is not specified, the size of the vector is adjusted
!   to its original size.
!
!
! !INTERFACE:

    subroutine resize_(nav,lsize)
      implicit none
      type(Navigator),intent(inout) :: nav
      integer,optional,intent(in) :: lsize

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::resize_'
  integer :: m

  m=msize_(nav)
  nav%lsize=m
  if(present(lsize)) nav%lsize=lsize
end subroutine resize_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: get_ - get an entry according to user preference.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine get_(nav,inav,displ,count,lc,ln,le)
      implicit none
      type(Navigator),intent(in) :: nav
      integer,intent(in) :: inav
      integer,optional,intent(out) :: displ
      integer,optional,intent(out) :: count
      integer,optional,intent(out) :: lc
      integer,optional,intent(out) :: ln
      integer,optional,intent(out) :: le

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::get_'

	! No checking is done here to ensure the entry is a valid.
	! The user must ensure the index is valid.

  if(present(displ)) displ=nav%displs(inav)
  if(present(count)) count=nav%counts(inav)
  if(present(lc)) lc=nav%displs(inav)+1
  if(present(ln)) ln=nav%counts(inav)
  if(present(le)) le=nav%displs(inav)+nav%counts(inav)

end subroutine get_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_displs_ - returns pointer to displs(:) component.
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_displs_(nav,lbnd,ubnd)
      implicit none
      type(Navigator),intent(in) :: nav
      integer,optional,intent(in) :: lbnd
      integer,optional,intent(in) :: ubnd
      integer,pointer,dimension(:) :: ptr_displs_

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_displs_'
  integer :: lc,le

  if(present(lbnd).or.present(ubnd)) then
    lc=lbound(nav%displs,1)
    if(present(lbnd)) lc=lbnd
    le=ubound(nav%displs,1)
    if(present(ubnd)) le=ubnd
    ptr_displs_ => nav%displs(lc:le)
  else
    le=nav%lsize
    ptr_displs_ => nav%displs(1:le)
  endif

end function ptr_displs_
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ptr_counts_ - returns pointer to counts(:) component.
!
! !DESCRIPTION:
!
! !INTERFACE:

    function ptr_counts_(nav,lbnd,ubnd)
      implicit none
      type(Navigator),intent(in) :: nav
      integer,optional,intent(in) :: lbnd
      integer,optional,intent(in) :: ubnd
      integer,pointer,dimension(:) :: ptr_counts_

! !REVISION HISTORY:
! 	22May00	- Jing Guo <guo@dao.gsfc.nasa.gov>
!		- initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::ptr_counts_'
  integer :: lc,le

  if(present(lbnd).or.present(ubnd)) then
    lc=lbound(nav%counts,1)
    if(present(lbnd)) lc=lbnd
    le=ubound(nav%counts,1)
    if(present(ubnd)) le=ubnd
    ptr_counts_ => nav%counts(lc:le)
  else
    le=nav%lsize
    ptr_counts_ => nav%counts(1:le)
  endif

end function ptr_counts_

end module m_Navigator
