!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_AttrVect - a distributed Innovation vector
!
! !DESCRIPTION:
!
! An {\em attribute vector} is a scheme for storing bundles of integer 
! and real data vectors, indexed by lists of their respective attributes.
! The attribute vector is implemented in Fortran 90 using the 
! {\tt AttrVect} derived type.  This module contains the definition of
! {\tt AttrVect} class, and numerous methods that service it.
!
! !INTERFACE:

 module m_AttrVect
!
! !USES:
!
      use m_List, only : List   ! Support for rList and iList components.

      implicit none

      private	! except

      public :: AttrVect        ! The class data structure

      public :: init		! create a local vector
      public :: clean		! clean the local vector
      public :: zero		! zero the local vector
      public :: lsize		! size of the local vector
      public :: nIAttr		! number of integer attributes on local
      public :: nRAttr		! number of real attributes on local
      public :: indexIA		! index the integer attributes
      public :: indexRA		! index the real attributes
      public :: getIList        ! return list of integer attributes
      public :: getRList        ! return list of real attributes
      public :: Sort            ! sort entries, and return permutation
      public :: Permute         ! permute entries
      public :: SortPermute     ! sort and permute entries

    type AttrVect
      type(List) :: iList
      type(List) :: rList
      integer,dimension(:,:),pointer :: iAttr
      real   ,dimension(:,:),pointer :: rAttr
    end type AttrVect

    interface init   ; module procedure	&
	init_,	&
	initv_, &
	initl_
    end interface
    interface clean  ; module procedure clean_  ; end interface
    interface zero  ; module procedure zero_  ; end interface
    interface lsize  ; module procedure lsize_  ; end interface
    interface nIAttr ; module procedure nIAttr_ ; end interface
    interface nRAttr ; module procedure nRAttr_ ; end interface
    interface indexIA; module procedure indexIA_; end interface
    interface indexRA; module procedure indexRA_; end interface
    interface getIList; module procedure getIList_; end interface
    interface getRList; module procedure getRList_; end interface
    interface Sort    ; module procedure Sort_    ; end interface
    interface Permute ; module procedure Permute_ ; end interface
    interface SortPermute ; module procedure SortPermute_ ; end interface

! !REVISION HISTORY:
! 	10Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	10Oct00 - J.W. Larson <larson@mcs.anl.gov> - made getIList
!                 and getRList functions public and added appropriate
!                 interface definitions
!       20Oct00 - J.W. Larson <larson@mcs.anl.gov> - added Sort, 
!                 Permute, and SortPermute functions.
!       09May01 - J.W. Larson <larson@mcs.anl.gov> - added initl_().
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_AttrVect'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize with given iList, rList, and the size
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine init_(aV,iList,rList,lsize)
!
! !USES:
!
      use m_List, only : List
      use m_List, only : init,nitem
      use m_List, only : List_nullify => nullify
      use m_mall
      use m_die
      implicit none
      type(AttrVect),intent(out) :: aV
      character(len=*),optional,intent(in) :: iList
      character(len=*),optional,intent(in) :: rList
      integer,         optional,intent(in) :: lsize

! !REVISION HISTORY:
! 	09Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	09Oct01 - J.W. Larson <larson@mcs.anl.gov> - added feature to
!                 nullify all pointers before usage.  This was done to
!                 accomodate behavior of the f90 ASSOCIATED intrinsic 
!                 function on the AIX platform.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: nIA,nRA,n,ier

       ! Initially, nullify all pointers in the AttrVect aV:

  nullify(aV%iAttr)
  nullify(aV%rAttr)
  call List_nullify(aV%iList)
  call List_nullify(aV%rList)

  if(present(rList)) then
    call init(aV%rList,rList)	! init.List()
  endif

  if(present(iList)) then
    call init(aV%iList,iList)	! init.List()
  endif

  nIA=nitem(aV%iList)		! nitem.List()
  nRA=nitem(aV%rList)		! nitem.List()

  n=0
  if(present(lsize)) n=lsize

  allocate( aV%iAttr(nIA,n),aV%rAttr(nRA,n),	stat=ier)
  if(ier /= 0) call perr_die(myname_,'allocate()',ier)

#ifdef MALL_ON
	call mall_ci(size(transfer(aV%iAttr,(/1/)),myname_)
	call mall_ci(size(transfer(aV%rAttr,(/1/)),myname_)
#endif

 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initv_ - initialize on the vectors
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine initv_(aV,bV,lsize)
!
! !USES:
!
      use m_String, only : String,char
      use m_List,   only : get
      use m_die
      use m_stdio

      implicit none

      type(AttrVect),intent(out) :: aV
      type(AttrVect),intent(in)  :: bV
      integer,       intent(in)  :: lsize

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!       17May01 - R. Jacob <jacob@mcs.anl.gov> - add a check to see if
!                 input argument has been defined.  SGI will dump
!                 core if its not.
!       10Oct01 - J. Larson <larson@mcs.anl.gov> - Nullify all pointers
!                 in ouput AttrVect aV before initializing aV.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initv_'
  type(String) :: iLStr,rLStr

	! Step One:  Nullify all pointers in aV.  We will set
        ! only the pointers we really need for aV based on those
        ! currently ASSOCIATED in bV.

  call List_nullify(aV%iList)
  call List_nullify(aV%rList)
  nullify(aV%iAttr)
  nullify(aV%rAttr)

	! Convert the two Lists to two Strings

  if(.not.associated(bv%iList%bf) .and. & 
       .not.associated(bv%rList%bf)) then
     write(stderr,'(2a)')myname_, &
      'MCTERROR:  Trying to initialize a new AttrVect off an undefined AttrVect'
      call perr_die(myname_,'undefined input argument',0)
  endif

  if(associated(bv%iList%bf)) then
     call get(iLStr,bv%iList)
  endif

  if(associated(bv%rList%bf)) then
     call get(rLStr,bv%rList)
  endif

       ! Initialize the AttrVect aV depending on which parts of
       ! the input bV are valid:

  if(associated(bv%iList%bf) .and. associated(bv%rList%bf)) then
     call init_(aV,iList=char(iLStr),rList=char(rLStr),lsize=lsize)
  endif
  if(.not.associated(bv%iList%bf) .and. associated(bv%rList%bf)) then
     call init_(aV,rList=char(rLStr),lsize=lsize)
  endif
  if(associated(bv%iList%bf) .and. .not.associated(bv%rList%bf)) then
     call init_(aV,iList=char(iLStr),lsize=lsize)
  endif

 end subroutine initv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initl_ - initialize with Lists iList, rList, and the size
!
! !DESCRIPTION:  This routine initializes an {\tt AttrVect} directly 
! from input {\tt List} data type arguments {\tt iList} and {\tt rList} 
! (see the module {\tt m\_List} in mpeu for further details), and an
! input length {\tt lsize}.  The resulting {\tt AttrVect} is returned in
! the argument {\tt aV}.
!
! {\bf N.B.}:  The outcome of this routine, {\tt aV} represents 
! allocated memory.  When this {\tt AttrVect} is no longer needed, 
! it must be deallocated by invoking the routine {\tt AttrVect\_clean()}.  
! Failure to do so will spawn a memory leak.
!
! !INTERFACE:

 subroutine initl_(aV, iList, rList, lsize)

!
! !USES:
!
      use m_die
      use m_stdio
      use m_mall

      use m_String, only : String
      use m_String, only : String_clean => clean
      use m_String, only : String_toChar => toChar

      use m_List, only : List
      use m_List, only : List_init => init
      use m_List, only : List_nitem => nitem
      use m_List, only : List_copy => copy
      use m_List, only : List_get => get

      implicit none

! !INPUT PARAMETERS: 
!
      type(List), intent(in)  :: iList
      type(List), intent(in)  :: rList
      integer,    intent(in)  :: lsize

! !OUTPUT PARAMETERS:
!
      type(AttrVect), intent(out) :: aV

! !REVISION HISTORY:
! 	09May98 - J.W. Larson <larson@mcs.anl.gov> - initial version.
!       08Aug01 - E.T. Ong <eong@mcs.anl.gov> - change list assignment(=)
!                 to list copy to avoid compiler errors with pgf90.
!       10Oct01 - J. Larson <larson@mcs.anl.gov> - Nullify all pointers
!                 in ouput AttrVect aV before initializing aV.  Also, 
!                 greater caution taken regarding validity of input 
!                 arguments iList and rList.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::initl_'
  integer :: nIA, nRA, ier

	! Step One:  Nullify all pointers in aV.  We will set
        ! only the pointers we really need for aV based on those
        ! currently ASSOCIATED in bV.

  call List_nullify(aV%iList)
  call List_nullify(aV%rList)
  nullify(aV%iAttr)
  nullify(aV%rAttr)

       ! Assign iList to aV%iList and rList to aV%rList

  if(associated(iList%bf)) then
     call List_copy(aV%iList,iList)
  endif

  if(associated(rList%bf)) then
     call List_copy(aV%rList,rList)
  endif

       ! Count items in aV%iList and aV%rList: 

  if(associated(aV%iList%bf)) then
     nIA = List_nitem(aV%iList)
  else
     nIA = 0
  endif

  if(associated(aV%rList%bf)) then
     nRA = List_nitem(aV%rList)
  else
     nRA = 0
  endif

  if(lsize <= 0) then 
     write(stdout,*) myname_,":: Warning:  length argument lsize nonpositive."
  endif

  if(lsize > 0) then

     if(associated(aV%iList%bf) .and. associated(aV%rList%bf)) then
	allocate( aV%iAttr(nIA,lsize), aV%rAttr(nRA,lsize), stat=ier)
	if(ier /= 0) call perr_die(myname_,'allocate(aV%iAttr,aV%rAttr)',ier)
#ifdef MALL_ON
	call mall_ci(size(transfer(aV%iAttr,(/1/)),myname_)
	call mall_ci(size(transfer(aV%rAttr,(/1/)),myname_)
#endif
     endif

     if(associated(aV%iList%bf) .and. .not.associated(aV%rList%bf)) then
	allocate( aV%iAttr(nIA,lsize), stat=ier)
	if(ier /= 0) call perr_die(myname_,'allocate(aV%iAttr)',ier)
#ifdef MALL_ON
	call mall_ci(size(transfer(aV%iAttr,(/1/)),myname_)
#endif
     endif

     if(.not.associated(aV%iList%bf) .and. associated(aV%rList%bf)) then
	allocate( aV%rAttr(nRA,lsize), stat=ier)
	if(ier /= 0) call perr_die(myname_,'allocate(aV%rAttr)',ier)
#ifdef MALL_ON
	call mall_ci(size(transfer(aV%rAttr,(/1/)),myname_)
#endif
     endif

  endif ! if(lsize > 0)... block

 end subroutine initl_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean a vector
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine clean_(aV)
!
! !USES:
!
      use m_mall
      use m_stdio
      use m_die
      use m_List, only : List_clean => clean

      implicit none

      type(AttrVect),intent(inout) :: aV

! !REVISION HISTORY:
! 	09Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	10Oct01 - J. Larson <larson@mcs.anl.gov> - various fixes to 
!                 prevent deallocation of UNASSOCIATED pointers.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

	! Note that an undefined pointer may either crash the process
	! or return either .true. or .false. to the associated() test.
	! One should therefore avoid using the function on an 
	! undefined pointer.

        ! Clean up attribute lists:

  if(associated(aV%iList%bf)) call List_clean(aV%iList)
  if(associated(aV%rList%bf)) call List_clean(aV%rList)

        ! Clean up INTEGER attributes:

  if(associated(aV%iAttr)) then
     deallocate(aV%iAttr, stat=ier)
     if(ier /= 0) call perr_die(myname_,'deallocte(aV%iAttr)',ier)
#ifdef MALL_ON
     call mall_co(size(transfer(aV%iAttr,(/1/)),myname_)
#endif
  endif ! if(associated(aV%iAttr))...

        ! Clean up REAL attributes:

  if(associated(aV%rAttr)) then
     deallocate(aV%rAttr, stat=ier)
     if(ier /= 0) call perr_die(myname_,'deallocte(aV%rAttr)',ier)
#ifdef MALL_ON
     call mall_co(size(transfer(aV%rAttr,(/1/)),myname_)
#endif
  endif ! if(associated(aV%rAttr))...

 end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - the local size of the vector
!
! !DESCRIPTION:
!
! !INTERFACE:

 integer function lsize_(aV)

! !USES:

     use m_stdio, only : stderr
     use m_die
 
     implicit none

      type(AttrVect), intent(in) :: aV

! !REVISION HISTORY:
! 	09Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	10Oct01 - J. Larson <larson@mcs.anl.gov> - made code more robust
!                 to handle cases where the length of either aV%iAttr or
!                 aV%rAttr is zero, but the other is positive.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'
  integer :: iLength, rLength

	! One should try to avoid using this function on an undefined
	! or disassocated pointer.  However, it is understandable 
	! that an undefined or disassocated pointer has a size 0, if
	! the associated() test sucesses.

  lsize_=0

  if(associated(aV%iAttr)) then
     iLength = size(aV%iAttr,2)
  else
     iLength = 0
  endif

  if(associated(aV%rAttr)) then
     rLength = size(aV%rAttr,2)
  else
     rLength = 0
  endif

  if(iLength /= rLength) then

     if((rLength > 0) .and. (iLength > 0)) then
	call perr_die(myname_,'attribute array length mismatch', &
	     iLength-rLength)
     endif

     if((rLength > 0) .and. (iLength == 0)) then
	lsize_ = rLength
     endif

     if((iLength > 0) .and. (rLength == 0)) then
	lsize_ = iLength
     endif

  endif

  if(iLength == rLength) lsize_ = iLength

 end function lsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        Math and Computer Science Division, Argonne National Laboratory
!BOP -------------------------------------------------------------------
!
! !IROUTINE: zero_ - zero the local vector
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine zero_(aV)

! !USES:

     use m_die,only	: perr_die
     use m_stdio,only	: stderr
 
     implicit none

     type(AttrVect), intent(inout) :: aV

! !REVISION HISTORY:
! 	17May01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::zero_'
  integer :: i,j,size

  if(.not.associated(aV%iAttr) .and. .not.associated(aV%rAttr)) then
    write(stderr,'(2a)')myname_, &
      'MCTERROR:  Trying to zero an uninitialized AttrVect'
      call perr_die(myname_,'undefined input argument',0)
  endif

  size = lsize_(aV)
  do i=1,nIAttr_(aV)
   do j=1,size
    aV%iAttr(i,j)=0
   enddo
  enddo

  do i=1,nRAttr_(aV)
   do j=1,size
    aV%rAttr(i,j)=0.
   enddo
  enddo

 end subroutine zero_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nIAttr_ - number of INTEGER type attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

 function nIAttr_(aV)
!
! !USES:
!
      use m_List, only : nitem

      implicit none
      type(AttrVect),intent(in) :: aV
      integer :: nIAttr_

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	10Oct01 - J. Larson <larson@mcs.anl.gov> - made code more robust
!                 by checking status of pointers in aV%iList
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nIAttr_'

  if(associated(aV%iList%bf)) then
     nIAttr_ = nitem(aV%iList)
  else
     nIAttr_ = 0
  endif

 end function nIAttr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nRAttr_ - number of REAL type attributes
!
! !DESCRIPTION:
!
! !INTERFACE:

 function nRAttr_(aV)
!
! !USES:
!
      use m_List, only : nitem

      implicit none

      type(AttrVect),intent(in) :: aV
      integer :: nRAttr_

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	10Oct01 - J. Larson <larson@mcs.anl.gov> - made code more robust
!                 by checking status of pointers in aV%iList
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nRAttr_'

  if(associated(aV%rList%bf)) then
     nRAttr_ = nitem(aV%rList)
  else
     nRAttr_ = 0
  endif

 end function nRAttr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getIList_ - get an item from iList
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine getIList_(item,ith,aVect)
!
! !USES:
!
      use m_String, only : String
      use m_List,   only : get

      implicit none

      type(String),intent(out) :: item
      integer,     intent(in)  :: ith
      type(AttrVect),intent(in) :: aVect

! !REVISION HISTORY:
! 	24Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getIList_'

  call get(item,ith,aVect%iList)

 end subroutine getIList_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getRList_ - get an item from rList
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine getRList_(item,ith,aVect)
!
! !USES:
!
      use m_String, only : String
      use m_List,   only : get

      implicit none

      type(String),intent(out) :: item
      integer,     intent(in)  :: ith
      type(AttrVect),intent(in) :: aVect

! !REVISION HISTORY:
! 	24Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getRList_'

  call get(item,ith,aVect%rList)

 end subroutine getRList_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexIA_ - index the integer attribute List
!
! !DESCRIPTION:
!
! !INTERFACE:

 function indexIA_(aV,item,perrWith,dieWith)
!
! !USES:
!
      use m_List, only : index
      use m_die,  only : die
      use m_stdio,only : stderr

      implicit none

      type(AttrVect), intent(in) :: aV
      character(len=*),intent(in) :: item
      character(len=*),optional,intent(in) :: perrWith
      character(len=*),optional,intent(in) :: dieWith
      integer :: indexIA_

! !REVISION HISTORY:
! 	27Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexIA_'

  indexIA_=index(aV%iList,item)

	if(indexIA_==0) then
	  if(.not.present(dieWith)) then
	    if(present(perrWith)) write(stderr,'(4a)') perrWith, &
		'" indexIA_() error, not found "',trim(item),'"'
	  else
	    write(stderr,'(4a)') dieWith,	&
		'" indexIA_() error, not found "',trim(item),'"'
	    call die(dieWith)
	  endif
	endif

 end function indexIA_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexRA_ - index the integer attribute List
!
! !DESCRIPTION:
! 
!
! !INTERFACE:

 function indexRA_(aV,item,perrWith,dieWith)
!
! !USES:
!
      use m_List, only : index
      use m_die,  only : die
      use m_stdio,only : stderr

      implicit none

      type(AttrVect), intent(in) :: aV
      character(len=*),intent(in) :: item
      character(len=*),optional,intent(in) :: perrWith
      character(len=*),optional,intent(in) :: dieWith
      integer :: indexRA_

! !REVISION HISTORY:
! 	27Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexRA_'

  indexRA_=index(aV%rList,item)

	if(indexRA_==0) then
	  if(.not.present(dieWith)) then
	    if(present(perrWith)) write(stderr,'(4a)') perrWith, &
		'" indexRA_() error, not found "',trim(item),'"'
	  else
	    write(stderr,'(4a)') dieWith,	&
		'" indexRA_() error, not found "',trim(item),'"'
	    call die(dieWith)
	  endif
	endif

 end function indexRA_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Sort_ - return index permutation keyed by a list of
!            attributes
!
! !DESCRIPTION:
! The subroutine {\tt Sort\_()} uses a list of keys defined by the {\tt List} 
! {\tt sList}, searches for the appropriate integer or real attributes
! referenced by the items in {\tt sList} ( that is, it identifies the 
! appropriate entries in {aV\%iList} and {\tt aV\%rList}), and then 
! uses these keys to generate a permutation {\tt perm} that will put
! the entries of the attribute vector {\tt aV} in lexicographic order
! as defined by {\tt sList} (the ordering in {\tt sList} being from
! left to right.
!
! {\bf N.B.:}  This routine will fail if {\tt aV\%rList} and 
! {\tt aV\%rList} share one or more common entries. 
!
! {\bf N.B.:}  This routine will fail if {\tt aV\%rList} and 
!
! !INTERFACE:

 subroutine Sort_(aV, sList, perm, descend, perrWith, dieWith)
!
! !USES:
!
      use m_String,        only : String
      use m_String,        only : String_tochar => tochar
      use m_List ,         only : List_index => index
      use m_List ,         only : List_nitem => nitem
      use m_List ,         only : List_get   => get
      use m_die ,          only : die
      use m_stdio ,        only : stderr
      use m_SortingTools , only : IndexSet
      use m_SortingTools , only : IndexSort

      implicit none

      type(AttrVect), intent(in)                  :: aV
      type(List),     intent(in)                  :: sList
      integer, dimension(:), pointer              :: perm
      logical, dimension(:), optional, intent(in) :: descend
      character(len=*), optional, intent(in)      :: perrWith
      character(len=*), optional, intent(in)      :: dieWith


! !REVISION HISTORY:
! 	20Oct00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!       25Apr01 - R.L. Jacob <jacob@mcs.anl.gov> - add -1 to make a
!                 backwards loop go backwards
!       14Jun01 - J. Larson / E. Ong -- Fixed logic bug in REAL attribute
!                 sort (discovered by E. Ong), and cleaned up error / 
!                 shutdown logic.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::Sort_'

! local variables

        ! storage for key extracted from sList:

  type(String) :: key

        ! number of keys, loop index, error flag, and length:

  integer :: nkeys, n, ierr, length

        ! key indices for av%rAttr and av%iAttr, respectively:

  integer, dimension(:), allocatable :: rIndex, iIndex 

        ! count the sorting keys:

  nkeys = List_nitem(sList)

        ! allocate and initialize rIndex and iIndex to
        ! zero (the null return values from the functions
        ! indexRA_() and indexIA_() ).

  allocate(rIndex(nkeys), iIndex(nkeys), stat=ierr)
     
  rIndex = 0
  iIndex = 0

        ! Loop over the keys in the list, and identify the
        ! appropriate integer or real attribute, storing the
        ! attribute index in iIndex(:) or rIndex(:), respectively.

  do n = 1, nkeys

        ! grab the next key

     call List_get(key, n, sList)

        ! determine wheter this key refers to an
        ! integer or real attribute:

     rIndex(n) = indexRA_(aV, String_tochar(key))
     iIndex(n) = indexIA_(aV, String_tochar(key))

        ! If both rIndex(n) and iIndex(n) are greater than
        ! zero, then we have an integer attribute sharing 
        ! the same name as a real attribute, and there is 
        ! no clear path as to which one is the sort key.
        ! This is a fatal error that triggers shutdown.

     if ((rIndex(n) > 0) .and. (iIndex(n) > 0)) then
	if(.not.present(dieWith)) then
	   if(present(perrWith)) write(stderr,'(4a)') myname, &
		":: ambiguous key, ", perrWith, &
		" both iIndex(n) and rIndex(n) positive."
	    call die(myname_,":: both iIndex(n) and rIndex(n) > 0.")
	else
	   if(present(perrWith)) then
	       write(stderr,'(4a)') myname_,":: ", perrWith, &
		       " both iIndex(n) and rIndex(n) positive."
           endif
	   call die(myname_,dieWith)
	endif
     endif

        ! If both rIndex(n) and iIndex(n) are nonpositive,
        ! then the requested sort key is not present in either
        ! aV%rList or aV%iList, and we cannot perform the sort.
        ! This is a fatal error that triggers shutdown.

     if ((rIndex(n) <= 0) .and. (iIndex(n) <= 0)) then
	if(.not.present(dieWith)) then
	   if(present(perrWith)) write(stderr,'(4a)') myname,":: ", &
		perrWith, &
		" both iIndex(n) and rIndex(n) nonpositive"
	   call die(myname_,":: both iIndex(n) and rIndex(n) <= 0.")
	else
	   if(present(perrWith)) then
	      write(stderr,'(4a)') myname_,":: ", perrWith,	&
		   " both iIndex(n) and rIndex(n) nonpositive"
           endif
	   call die(myname_,dieWith)
	endif
     endif

        ! If only one of rIndex(n) or iIndex(n) is positive,
        ! set the other value to zero.

     if (iIndex(n) > 0) rIndex(n) = 0
     if (rIndex(n) > 0) iIndex(n) = 0

  enddo ! do n=1,nkeys

        ! Now we have the locations of the keys in the integer and
        ! real attribute storage areas aV%iAttr and aV%rAttr, respectively.
        ! our next step is to construct and initialize the permutation
        ! array perm.  First step--determine the length of aV using 
        ! lsize_():

  length = lsize_(aV)

  allocate(perm(length), stat=ierr)

        ! Initialize perm(i)=i, for i=1,length

  call IndexSet(perm)

        ! Now we can perform the stable successive keyed sorts to
        ! transform perm into the permutation that will place the
        ! entries of the attribute arrays in the lexicographic order
        ! defined by sList.  This is achieved by successive calls to
        ! IndexSort(), but in reverse order to the order of the keys
        ! as they appear in sList.

  do n=nkeys, 1, -1
     if(iIndex(n) > 0) then
	if(present(descend)) then
	   call IndexSort(length, perm, aV%iAttr(iIndex(n),:), &
		          descend(n))
	else
   	   call IndexSort(length, perm, aV%iAttr(iIndex(n),:), &
               		  descend=.false.)
	endif ! if(present(descend)...
     else
	if(rIndex(n) > 0) then
	   if(present(descend)) then
	      call IndexSort(length, perm, aV%rAttr(rIndex(n),:), &
		             descend(n))
	   else
	      call IndexSort(length, perm, aV%rAttr(rIndex(n),:), &
               		     descend=.false.)
	   endif ! if(present(descend)...
	endif ! if (rIndex(n) > 0)...
     endif ! if (iIndex(n) > 0)...
  enddo

        ! Now perm(1:length) is the transformation we seek--we are
        ! finished.

  deallocate(iIndex, rIndex, stat=ierr)  ! clean up allocated arrays.

 end subroutine Sort_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: Permute_ - return index permutation keyed by a list of
!            attributes
!
! !DESCRIPTION:
! The subroutine {\tt Permute\_()} uses a a permutation {\tt perm} (which can
! be generated by the routine {\tt Sort\_()} in this module) to rearrange
! the entries in the attribute integer and real storage areas of the
! input attribute vector {\tt aV}--{\tt aV\%iAttr} and {\tt aV\%rAttr}, 
! respectively.
!
! !INTERFACE:

 subroutine Permute_(aV, perm, perrWith, dieWith)
!
! !USES:
!
      use m_die ,          only : die
      use m_stdio ,        only : stderr
      use m_SortingTools , only : Permute

      implicit none

      type(AttrVect), intent(inout) :: aV
      integer, dimension(:), intent(in) :: perm
      character(len=*),optional,intent(in) :: perrWith
      character(len=*),optional,intent(in) :: dieWith

! !REVISION HISTORY:
! 	23Oct00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::Permute_'

! local variables

  integer, dimension(:,:), allocatable :: iAtmp
  real,    dimension(:,:), allocatable :: rAtmp
  integer :: i

        ! Check input arguments for compatibility--assure
        ! lsize_(aV) = size(perm); that is, make sure the
        ! index permutation is the same length as the vectors
        ! it will re-arrange.

  if (size(perm) /= lsize_(aV)) then
     if(.not.present(dieWith)) then
	if(present(perrWith)) write(stderr,'(4a,i8,a,i8)') myname, &
	  ":: size mismatch, ", perrWith, &
	  "size(perm)=",size(perm)," lsize_(aV)=",lsize_(aV)
     else
	write(stderr,'(4a,i8,a,i8)') myname, &
	 ":: size mismatch, ", dieWith,	&
	 "size(perm)=",size(perm)," lsize_(aV)=",lsize_(aV)
	call die(dieWith)
     endif
  endif

  if(size(perm) == lsize_(aV)) then

        ! Permute integer attributes:
     if(nIAttr_(aV) /= 0) then
	do i=1,nIAttr_(aV)
	   call Permute(aV%iAttr(i,:),perm,lsize_(aV))
	end do
     endif

        ! Permute real attributes:
     if(nRAttr_(aV) /= 0) then
	do i=1,nRAttr_(aV)
	   call Permute(aV%rAttr(i,:),perm,lsize_(aV))
	end do
     endif

  endif

 end subroutine Permute_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: SortPermute_ - Place AttrVect data in lexicographic order.
!
! !DESCRIPTION:
!
! The subroutine {\tt SortPermute\_()} uses the routine {\tt Sort\_()} 
! to create an index permutation {\tt perm} that will place the AttrVect
! entries in the lexicographic order defined by the keys in the List 
! variable {\tt key\_list}.  This permutation is then used by the routine
! {\tt Permute\_()} to place the AttreVect entries in lexicographic order.
!
! !INTERFACE:

 subroutine SortPermute_(aV, key_list, descend, perrWith, dieWith)
!
! !USES:
!
      use m_die ,          only : die
      use m_stdio ,        only : stderr

      implicit none

      type(AttrVect), intent(inout) :: aV
      type(List), intent(in)        :: key_list
      logical , dimension(:), optional, intent(in) :: descend
      character(len=*), optional, intent(in) :: perrWith
      character(len=*), optional, intent(in) :: dieWith


! !REVISION HISTORY:
! 	24Oct00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::Permute_'

! local variables

       ! Permutation array pointer perm(:)
  integer, dimension(:), pointer :: perm
       ! Error flag ierr
  integer :: ierr

       ! Step One: Generate the index permutation perm(:)

  if(present(descend)) then
     call Sort_(aV, key_list, perm, descend, perrWith, dieWith)
  else
     call Sort_(aV, key_list, perm, perrWith=perrWith, &
	        dieWith=dieWith)
  endif

       ! Step Two:  Apply the index permutation perm(:)

  call Permute_(aV, perm, perrWith, dieWith)

       ! Step Three:  deallocate temporary array used to 
       ! store the index permutation (this was allocated 
       ! in the routine Sort_()

  deallocate(perm, stat=ierr)

  end subroutine SortPermute_

! Sorting:
!
!	aV%iVect(:,:) =		&
!		aV%iVect((/(indx(i),i=1,lsize(aV))/),:)
!
!	aV%iVect((/(indx(i),i=1,lsize(aV))/),:) =	&
!		aV%iVect(:,:)
!
!	aV%iVect(:,ikx),aV%iVect(:,iks)

 end module m_AttrVect
!.

