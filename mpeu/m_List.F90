!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_List - a list manager
!
! !DESCRIPTION:
!
! !INTERFACE:

    module m_List
      implicit none
      private	! except

      public :: List		! The class data structure
      public :: init
      public :: clean
      public :: nullify
      public :: index
      public :: nitem
      public :: get
      public :: identical
      public :: assignment(=)
      public :: allocated
      public :: copy
      public :: exportToChar
      public :: CharBufferSize
      public :: concatenate
      public :: bcast
      public :: send
      public :: recv

    type List
      character(len=1),dimension(:),pointer :: bf
      integer,       dimension(:,:),pointer :: lc
    end type List

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	16May01 - J. Larson <larson@mcs.anl.gov> - Several changes / fixes:
!                 public interface for copy_(), corrected version of copy_(),
!                 corrected version of bcast_().
! 	15Oct01 - J. Larson <larson@mcs.anl.gov> - Added the LOGICAL 
!                 function identical_().
! 	14Dec01 - J. Larson <larson@mcs.anl.gov> - Added the LOGICAL 
!                 function allocated_().
! 	13Feb02 - J. Larson <larson@mcs.anl.gov> - Added the List query 
!                 functions exportToChar() and CharBufferLength().
!EOP ___________________________________________________________________

  interface init ; module procedure	&
	init_,		&
	initStr_,	&
	initstr1_
  end interface
  interface clean; module procedure clean_; end interface
  interface nullify; module procedure nullify_; end interface
  interface index; module procedure	&
	index_,		&
	indexStr_
  end interface
  interface nitem; module procedure nitem_; end interface
  interface get  ; module procedure	&
	get_,		&
	getall_,	&
	getrange_
  end interface
  interface identical; module procedure identical_; end interface
  interface assignment(=)
    module procedure copy_
  end interface
  interface allocated ; module procedure allocated_;  end interface
  interface copy ; module procedure copy_;  end interface
  interface exportToChar ; module procedure exportToChar_; end interface
  interface CharBufferSize 
     module procedure CharBufferSize_
  end interface
  interface concatenate ; module procedure concatenate_ ; end interface
  interface bcast; module procedure bcast_; end interface
  interface send; module procedure send_; end interface
  interface recv; module procedure recv_; end interface

  character(len=*),parameter :: myname='m_List'

contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialized a List from a character string
!
! !DESCRIPTION:
!
!	A list is a string in the form of ``\verb"cat:tiger:lion"'',
!   or ``\verb"lat:lon:lev"''.  Through the initialization call, the
!   items delimited by ``\verb":"'' are stored as an array of sub-
!   strings of a long string, accessible through an array of substring
!   indices.  The only constraints now on the valid list entries are,
!   (1) the value of an entry does not contain ``\verb":"'', and (2)
!   The leading and the trailing blanks are insignificant, although
!   any imbeded blanks are.
!
! !INTERFACE:

 subroutine init_(aList,Values)

      use m_die,only : die
      use m_mall,only : mall_mci,mall_ison
      implicit none
      type(List),intent(out)	  :: aList  ! an indexed string values
      character(len=*),intent(in) :: Values ! ":" delimited names

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::init_'
  character(len=1) :: c
  integer :: ib,ie,id,lb,le,ni,i,ier

	! Pass 1, getting the sizes
  le=0
  ni=0
  ib=1
  ie=0
  id=0
  do i=1,len(Values)
    c=Values(i:i)
    select case(c)
    case(' ')
      if(ib==i) ib=i+1	! moving ib up, starting from the next
    case(':')
      if(ib<=ie) then
	ni=ni+1
	id=1		! mark a ':'
      endif
      ib=i+1		! moving ib up, starting from the next
    case default
      ie=i
      if(id==1) then	! count an earlier marked ':'
	id=0
	le=le+1
      endif
      le=le+1
    end select
  end do
  if(ib<=ie) ni=ni+1

  allocate(aList%bf(le),aList%lc(0:1,ni),stat=ier)
  if(ier /= 0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(aList%bf,myname)
	  call mall_mci(aList%lc,myname)
	endif

	! Pass 2, copy the value and assign the pointers
  lb=1
  le=0
  ni=0
  ib=1
  ie=0
  id=0
  do i=1,len(Values)
    c=Values(i:i)

    select case(c)
    case(' ')
      if(ib==i) ib=i+1	! moving ib up, starting from the next
    case(':')
      if(ib<=ie) then
	ni=ni+1
	aList%lc(0:1,ni)=(/lb,le/)
	id=1		! mark a ':'
      endif

      ib=i+1		! moving ib up, starting from the next
      lb=le+2		! skip to the next non-':' and non-','
    case default
      ie=i
      if(id==1) then	! copy an earlier marked ':'
	id=0
	le=le+1
        aList%bf(le)=':'
      endif

      le=le+1
      aList%bf(le)=c
    end select
  end do
  if(ib<=ie) then
    ni=ni+1
    aList%lc(0:1,ni)=(/lb,le/)
  endif

 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initStr_ initialize with a String type
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine initStr_(aList,pstr)

      use m_String, only : String,toChar
      implicit none
      type(List),intent(out)	  :: aList  ! an indexed string values
      type(String),intent(in)	  :: pstr

! !REVISION HISTORY:
! 	23Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initStr_'

  call init_(aList,toChar(pstr))

 end subroutine initStr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initStr1_ initialize with an array of Strings
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine initStr1_(aList,strs)

      use m_String, only : String,toChar
      use m_String, only : len
      use m_String, only : ptr_chars
      use m_die,only : die
      implicit none
      type(List),intent(out)	  :: aList  ! an indexed string values
      type(String),dimension(:),intent(in)	  :: strs

! !REVISION HISTORY:
! 	23Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initStr1_'
  character(len=1),allocatable,dimension(:) :: ch1
  integer :: ier
  integer :: n,i,lc,le

  n=size(strs)
  le=0
  do i=1,n
    le=le+len(strs(i))
  end do
  le=le+n-1	! for n-1 ":"s

	allocate(ch1(le),stat=ier)
		if(ier/=0) call die(myname_,'allocate()',ier)

  le=0
  do i=1,n
    if(i>1) then
      le=le+1
      ch1(le)=':'
    endif

    lc=le+1
    le=le+len(strs(i))
    ch1(lc:le)=ptr_chars(strs(i))
  end do
    
  call init_(aList,toChar(ch1))

	deallocate(ch1,stat=ier)
		if(ier/=0) call die(myname_,'deallocate()',ier)

 end subroutine initStr1_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - clean a List variable
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine clean_(aList)

      use m_die,only : die
      use m_mall,only : mall_mco,mall_ison
      implicit none
      type(List),intent(inout) :: aList

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

	if(mall_ison()) then
	  call mall_mco(aList%bf,myname)
	  call mall_mco(aList%lc,myname)
	endif

  deallocate(aList%bf,aList%lc,stat=ier)
  if(ier /= 0) call die(myname_,'deallocate()',ier)

 end subroutine clean_

!BOP -------------------------------------------------------------------
!     Math + Computer Science Division / Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nullify_ - nullify a List variable
!
! !DESCRIPTION:  In Fortran 90, pointers may have three states of being:
! 1) {\tt ASSOCIATED}, that is the pointer is pointing at a target,2) 
! 2) {\tt UNASSOCIATED}, and 3) {\tt UNINITIALIZED}.  On some platforms, 
! the Fortran intrinsic function {\tt associated()} 
! will view uninitialized pointers as {\tt UNASSOCIATED} by default.
! This is not always the case.  It is good programming practice to 
! nullify pointers if they are not to be used.  This routine nullifies
! the pointers present in the {\tt List} datatype.
!
! !INTERFACE:

 subroutine nullify_(aList)

! !USES:
!
      use m_die,only : die

      implicit none

! !INPUT/OUTPUT PARAMETERS: 
!
      type(List),intent(inout) :: aList

! !REVISION HISTORY:
! 	18Jun01 - J.W. Larson - <larson@mcs.anl.gov> - initial version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nullify_'

  nullify(aList%bf)
  nullify(aList%lc)

 end subroutine nullify_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nitem_ - number of items in the list
!
! !DESCRIPTION:
!
! !INTERFACE:

 function nitem_(aList)

      implicit none
      type(List),intent(in) :: aList
      integer :: nitem_

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	10Oct01 - J.W. Larson <larson@mcs.anl.gov> - modified routine to
!                 check pointers aList%bf and aList%lc using  the f90 
!                 intrinsic ASSOCIATED before proceeding with the item
!                 count.  If these pointers are UNASSOCIATED, an item
!                 count of zero is returned.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nitem_'
  integer :: NumItems

       ! Initialize item count to zero

  NumItems = 0

       ! If the List pointers are ASSOCIATED, perform item count:

  if(ASSOCIATED(aList%bf) .and. ASSOCIATED(aList%lc)) then
     NumItems = size(aList%lc,2)
  endif

  nitem_ = NumItems

 end function nitem_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: index_ - lookup a list for a given item name
!
! !DESCRIPTION:
!
! !INTERFACE:

 function index_(aList,item)

      use m_String, only : toChar
      implicit none
      type(List),      intent(in) :: aList	! a List of names
      character(len=*),intent(in) :: item	! a given item name
      integer :: index_

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::index_'
  integer :: i,lb,le

  index_=0
  do i=1,size(aList%lc,2)		! == nitem_(aList)
    lb=aList%lc(0,i)
    le=aList%lc(1,i)
    if(item==toChar(aList%bf(lb:le))) then
      index_=i
      exit
    endif
  enddo

 end function index_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexStr_ - lookup a list for a given item name
!
! !DESCRIPTION:
!
! !INTERFACE:

 function indexStr_(aList,itemStr)

      use m_String,only : String,toChar
      implicit none
      type(List),      intent(in) :: aList	! a List of names
      type(String),    intent(in) :: itemStr
      integer :: indexStr_

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexStr_'
  integer :: i,lb,le

  indexStr_=0
  do i=1,size(aList%lc,2)		! == nitem_(aList)
    lb=aList%lc(0,i)
    le=aList%lc(1,i)
    if(toChar(itemStr)==toChar(aList%bf(lb:le))) then
      indexStr_=i
      exit
    endif
  enddo

 end function indexStr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: allocated_ - Check a list to see if it is allocated.
!
! !DESCRIPTION:
! This function checks the input {\tt List} argument {\tt inList} to 
! determine whether or not it has been allocated.  It does this by 
! invoking the Fortran90 intrinsic function {\tt associated()} on the 
! pointers {\tt inList\%bf} and {\tt inList\%lc}.  If both of these 
! pointers are associated, the return value is {\tt .TRUE.}.
!
! {\bf N.B.:}  In Fortran90, pointers have three different states of 
! existence:  {\tt ASSOCIATED}, {\tt UNASSOCIATED}, and {\tt UNDEFINED}.
!  If a pointer is {\tt UNDEFINED}, this function may return either 
! {\tt .TRUE.} or {\tt .FALSE.} values, depending on the Fortran90 
! compiler.  To avoid such problems, we advise that users invoke the 
! {\tt List} method {\tt nullify()} to nullify any {\tt List} pointers
! for {\tt List} variables that are not initialized.
!
! !INTERFACE:

 logical function allocated_(inList)

! !USES:

      use m_die,only : die

      implicit none

! !INPUT PARAMETERS:

      type(List), intent(in) :: inList

! !REVISION HISTORY:
! 	14Dec01 - J. Larson <larson@mcs.anl.gov> - inital version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::allocated_'

  allocated_ = associated(inList%bf) .and. associated(inList%lc)

 end function allocated_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: copy_ - Copy a List.  Pointers are copied as allocatables
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine copy_(yL,xL)	! yL=xL

      use m_die,only : die
      use m_String ,only : String
      use m_String ,only : String_clean
      use m_mall,only : mall_mci,mall_ison

      implicit none
      type(List),intent(out) :: yL
      type(List),intent(in)  :: xL

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
! 	16May01 - J. Larson <larson@mcs.anl.gov> - simpler, working 
!                 version that exploits the String datatype (see m_String)
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::copy_'
  type(String) DummStr

       ! Download input List info from xL to String DummStr

  call getall_(DummStr,xL)

       ! Initialize yL from DummStr

  call initStr_(yL,DummStr)

  call String_clean(DummStr)

 end subroutine copy_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportToChar_ - Export List to a CHARACTER.
!
! !DESCRIPTION:  This function returns the character buffer portion of
! the input {\tt List} argument {\tt inList}---that is, the contents of
! {\tt inList%bf}---as a {\tt CHARACTER} (suitable for printing).  An
! example of the use of this function is:
! \begin{verbatim}
!           write(*,*) exportToChar(inList) 
! \begin{verbatim}
!
! !INTERFACE:

 function exportToChar_(inList)

      use m_die,    only : die
      use m_stdio,  only : stderr
      use m_String, only : String
      use m_String, only : String_ToChar => toChar
      use m_String, only : String_clean

      implicit none

! ! INPUT PARAMETERS:

      type(List),         intent(in) :: inList

! ! OUTPUT PARAMETERS:

      character(len=size(inList%bf)) :: exportToChar_

! !REVISION HISTORY:
! 	13Feb02 - J. Larson <larson@mcs.anl.gov> - initial version.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::exportToChar_'
  type(String) DummStr

       ! Download input List info from inList to String DummStr
  if(allocated_(inList)) then
     call getall_(DummStr,inList)
     exportToChar_ = String_ToChar(DummStr)
     call String_clean(DummStr)
  else
     write(stderr,*) myname_,":: Argument inList not allocated."
     call die(myname_)
  endif

 end function exportToChar_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: CharBufferSize_ - Export List to a CHARACTER.
!
! !DESCRIPTION:  This function returns the length of the character 
! buffer portion of the input {\tt List} argument {\tt inList} (that 
! is, the number of characters stored in \tt inList%bf}) as an
! {\tt INTEGER}.  A usage example is presented below:
! \begin{verbatim}
!           integer :: BufferLength
!           BufferLength = CharBufferSize(inList) 
! \begin{verbatim}
!
! !INTERFACE:

 integer function CharBufferSize_(inList)

      use m_die,    only : die
      use m_stdio,  only : stderr

      implicit none

! ! INPUT PARAMETERS:

      type(List),         intent(in) :: inList

! !REVISION HISTORY:
! 	13Feb02 - J. Larson <larson@mcs.anl.gov> - initial version.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::CharBufferSize_'

  if(allocated_(inList)) then
     CharBufferSize_ = size(inList%bf)
  else
     write(stderr,*) myname_,":: Argument inList not allocated."
     call die(myname_)
  endif

 end function CharBufferSize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: get_ - return a numbered item from the List
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine get_(itemStr,ith,aList)
      use m_String, only : String,init,toChar
      implicit none
      type(String),intent(out) :: itemStr
      integer,     intent(in)  :: ith
      type(List),  intent(in)  :: aList

! !REVISION HISTORY:
! 	23Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::get_'
  integer :: lb,le

  if(ith>0 .and. ith <= size(aList%lc,2)) then
    lb=aList%lc(0,ith)
    le=aList%lc(1,ith)
    call init(itemStr,toChar(aList%bf(lb:le)))
  else
    call init(itemStr,'')
  endif

 end subroutine get_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getall_ - return all items from the List
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine getall_(itemStr,aList)
      use m_String, only : String,init,toChar
      implicit none
      type(String),intent(out) :: itemStr
      type(List),  intent(in)  :: aList

! !REVISION HISTORY:
! 	23Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getall_'
  integer :: lb,le,ni

  ni=size(aList%lc,2)
  lb=aList%lc(0,1)
  le=aList%lc(1,ni)
  call init(itemStr,toChar(aList%bf(lb:le)))

 end subroutine getall_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getrange_ - return a range of items from the List
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine getrange_(itemStr,i1,i2,aList)

      use m_String, only : String,init,toChar
      implicit none
      type(String),intent(out) :: itemStr
      integer,     intent(in)  :: i1
      integer,     intent(in)  :: i2
      type(List),  intent(in)  :: aList

! !REVISION HISTORY:
! 	23Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getrange_'
  integer :: lb,le,ni

  ni=size(aList%lc,2)
  lb=aList%lc(0,max(1,i1))
  le=aList%lc(1,min(ni,i2))
  call init(itemStr,toChar(aList%bf(lb:le)))

 end subroutine getrange_


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: identical_ - Compare two lists to see if they are the same.
!
! !DESCRIPTION:
!
! !INTERFACE:

 logical function identical_(yL,xL)

      use m_die,only : die
      use m_String ,only : String
      use m_String ,only : String_clean

      implicit none

      type(List),intent(in) :: yL
      type(List),intent(in) :: xL

! !REVISION HISTORY:
! 	14Oct01 - J. Larson <larson@mcs.anl.gov> - original version
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::identical_'

  logical :: myIdentical
  type(String) :: DummStr
  integer :: n, NumItems

       ! Compare the number of the items in the Lists xL and yL.
       ! If they differ, myIdentical is set to .FALSE. and we are
       ! finished.  If both Lists sport the same number of items,
       ! we must compare them one-by-one...

  if(nitem_(yL) == nitem_(xL)) then

     NumItems = nitem_(yL)

     COMPARE_LOOP:  do n=1,NumItems

	call get_(DummStr, n, yL)  ! retrieve nth tag as a String

	if( indexStr_(xL, Dummstr) /= n ) then ! a discrepency spotted.
	   call String_clean(Dummstr)
	   myIdentical = .FALSE.    
	   EXIT
	else
	   call String_clean(Dummstr)
	endif

	   myIdentical = .TRUE.   ! we survived the whole test process.

     end do COMPARE_LOOP

  else
     myIdentical = .FALSE.
  endif

  identical_ = myIdentical

 end function identical_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: set_indices_ - set the indices of given items
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine set_indices_(indices,aList,values)
      use m_String, only : String,clean
      implicit none
      integer,dimension(:),intent(out) :: indices
      type(List),intent(in)	  :: aList  ! an indexed string values
      character(len=*),intent(in) :: Values ! ":" delimited names

! !REVISION HISTORY:
! 	31May98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::set_indices_'
  type(List)   :: tList
  type(String) :: tStr
  integer :: n,i

  call init_(tList,values)
  n=min(nitem_(tList),size(indices))

  do i=1,n
    call get_(tStr,i,tList)
    indices(i)=indexStr_(aList,tStr)
  end do

  call clean_(tList)
  call clean(tStr)

 end subroutine set_indices_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: concatenate_ - Concatenates two Lists into a third List.
!
! !DESCRIPTION:  This routine takes two input {\tt List} arguments
! {\tt iList1} and {\tt iList2}, and concatenates them, producing an 
! output {\tt List} argument {\tt oList}.
!
! {\bf N.B.}:  The outcome of this routine is order dependent.  That is,
! the entries of {\tt iList2} will follow {\tt iList1}.
!
! {\bf N.B.}:  The outcome of this routine, {\tt oList} on non-root
! processes, represents allocated memory.  When this {\tt List} is
! no longer needed, it must be deallocated by invoking the routine 
! {\tt List\_clean()}.  Failure to do so will cause a memory leak.
!
! !INTERFACE:

    subroutine concatenate_(iList1, iList2, oList)
!
! !USES:
!
      use m_stdio
      use m_die, only : die

      use m_mpif90

      use m_String, only:  String
      use m_String, only:  String_toChar => toChar
      use m_String, only:  String_len
      use m_String, only:  String_clean => clean

      implicit none

! !INPUT PARAMETERS: 
!
      type(List),         intent(in)  :: iList1
      type(List),         intent(in)  :: iList2 

! !OUTPUT PARAMETERS: 
!
      type(List),         intent(out) :: oList

! !BUGS:  For now, the List concatenate algorithm relies on fixed-length
! CHARACTER variables as intermediate storage.  The lengths of these
! scratch variables is hard-wired to 10000, which should be large enough
! for most applications.  This undesirable feature should be corrected 
! ASAP.
!
! !REVISION HISTORY:
! 	08May01 - J.W. Larson - initial version.
! 	17May01 - J.W. Larson - Re-worked and tested successfully.
!EOP ___________________________________________________________________

 character(len=*),parameter :: myname_=myname//'::concatenate_'

 integer :: ilen1, ilen2, olen
 integer :: InNitems1, InNitems2, OutNitems
 integer :: ierr, n

 type(String) :: iStr1, iStr2
 character(10000) :: iChr1, iChr2, oChr

       ! First, handle the case of either iList1 and/or iList2 being
       ! null

  if((nitem_(iList1) == 0) .or. (nitem_(iList2) == 0)) then

     if((nitem_(iList1) == 0) .and. (nitem_(iList2) == 0)) then
	call init_(oList,'')
     endif

     if(nitem_(iList1) == 0) then
	call copy_(oList, iList2)
     endif

     if(nitem_(iList2) == 0) then
	call copy_(oList,iList1)
     endif

  else ! both lists are non-null

       ! Step one:  convert Lists to Strings

     call getall_(iStr1, iList1)
     call getall_(iStr2, iList2)

       ! Step two:  convert Strings to CHARACTER variables

     iChr1 = String_toChar(iStr1)
     iChr2 = String_toChar(iStr2)

     ilen1 = String_len(iStr1)
     ilen2 = String_len(iStr2)

       ! Step three:  concatenate CHARACTERs with the colon separator

     olen = ilen1 + ilen2 + 1

     oChr = trim(iChr1) // ':' // trim(iChr2)

       ! Step four:  initialize oList from a CHARACTER

     call init_(oList, trim(oChr))

       ! The concatenation is complete.  Now, clean up

     call String_clean(iStr1)
     call String_clean(iStr2)

  endif

 end subroutine concatenate_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: bcast_ - Broadcast a List variable...
!
! !DESCRIPTION:  This routine takes an input {\tt List} argument 
! {\tt iList} (on input, valid on the root only), and broadcasts it.
!
! {\bf N.B.}:  The outcome of this routine, {\tt ioList} on non-root
! processes, represents allocated memory.  When this {\tt List} is
! no longer needed, it must be deallocated by invoking the routine 
! {\tt List\_clean()}.  Failure to do so will cause a memory leak.
!
! !INTERFACE:

    subroutine bcast_(ioList, root, comm, status)
!
! !USES:
!
      use m_stdio
      use m_die, only : MP_perr_die

      use m_String, only:  String
      use m_String, only:  String_bcast => bcast
      use m_String, only:  String_clean => clean

      use m_mpif90

      implicit none

! !INPUT PARAMETERS: 
!
      integer,            intent(in)     :: root
      integer,            intent(in)     :: comm

! !INPUT/OUTPUT PARAMETERS: 
!
      type(List),         intent(inout)  :: ioList  


! !OUTPUT PARAMETERS: 
!
      integer, optional,  intent(out)    :: status

! !REVISION HISTORY:
! 	07May01 - J.W. Larson - initial version.
! 	14May01 - R.L. Jacob - fix error checking
! 	16May01 - J.W. Larson - new, simpler String-based algorigthm
!                 (see m_String for details), which works properly on
!                 the SGI platform.
!       13Jun01 - J.W. Larson <larson@mcs.anl.gov> - Initialize status
!                 (if present).
!EOP ___________________________________________________________________

 character(len=*),parameter :: myname_=myname//'::bcast_'
 integer :: myID, ierr
 type(String) :: DummStr

      ! Initialize status (if present)

  if(present(status)) status = 0

       ! Which process am I?

  call MPI_COMM_RANK(comm, myID, ierr)
  if(ierr /= 0) then
   if(present(status)) then
     status = ierr
     write(stderr,'(2a,i4)') myname_,":: MPI_COMM_RANK(), ierr=",ierr
     return
   else
     call MP_perr_die(myname_,"MPI_COMM_RANK()",ierr)
   endif
  endif

       ! on the root, convert ioList into the String variable DummStr

  if(myID == root) then
     call getall_(DummStr, ioList)
  endif

       ! Broadcast DummStr

  call String_bcast(DummStr, root, comm, ierr)
  if(ierr /= 0) then
   if(present(status)) then
     status = ierr
     write(stderr,'(2a,i4)') myname_,":: call String_bcast(), ierr=",ierr
     return
   else
     call MP_perr_die(myname_,"String_bcast() failed, stat=",ierr)
   endif
  endif

       ! Initialize ioList off the root using DummStr

  call initStr_(ioList, DummStr)

       ! And now, the List broadcast is complete.

  call String_clean(DummStr)

 end subroutine bcast_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: send_ - Point-to-point send of a List variable...
!
! !DESCRIPTION:  This routine takes an input {\tt List} argument 
! {\tt inList} and sends it to processor {\tt dest} on the communicator 
! associated with the fortran 90 {\tt INTEGER} handle {\tt comm}.  The 
! message is tagged by the input {\tt INTEGER} argument {\tt TagBase}.  
! The success (failure) of this operation is reported in the zero 
! (nonzero) optional output argument {\tt status}.
!
! {\bf N.B.}:  One must avoid assigning elsewhere the MPI tag values 
! {\tt TagBase} and {\tt TagBase+1}.  This is because {\tt send\_()} 
! performs the send of the {\tt List} as a pair of operations.  The 
! first send is the number of characters in {\tt inList\%bf}, and is 
! given MPI tag value {\tt TagBase}.  The second send is the 
! {\tt CHARACTER} data present in {\tt inList\%bf}, and is given MPI 
! tag value {\tt TagBase+1}.
!
! !INTERFACE:

    subroutine send_(inList, dest, TagBase, comm, status)
!
! !USES:
!
      use m_stdio
      use m_die, only : MP_perr_die

      use m_mpif90

      use m_String, only:  String
      use m_String, only:  String_toChar => toChar
      use m_String, only:  String_len
      use m_String, only:  String_clean => clean

      implicit none

! !INPUT PARAMETERS: 
!
      type(List),         intent(in)  :: inList  
      integer,            intent(in)  :: dest
      integer,            intent(in)  :: TagBase
      integer,            intent(in)  :: comm

! !OUTPUT PARAMETERS: 
!
      integer, optional,  intent(out) :: status

! !REVISION HISTORY:
! 	06Jun01 - J.W. Larson - initial version.
!       13Jun01 - J.W. Larson <larson@mcs.anl.gov> - Initialize status
!                 (if present).
!EOP ___________________________________________________________________

 character(len=*),parameter :: myname_=myname//'::send_'

 type(String) :: DummStr
 integer :: ierr, length

       ! Set status flag to zero (success) if present:

 if(present(status)) status = 0

       ! Step 1.  Extract CHARACTER buffer from inList and store it
       ! in String variable DummStr, determine its length.

 call getall_(DummStr, inList)
 length = String_len(DummStr)

       ! Step 2.  Send Length of String DummStr to process dest.

 call MPI_SEND(length, 1, MP_type(length), dest, TagBase, comm, ierr)
  if(ierr /= 0) then
     if(present(status)) then
        write(stderr,*) myname_,':: MPI_SEND(length...) failed.  ierror=',&
	     ierr
        status = ierr
        return
     else
        call MP_perr_die(myname_,':: MPI_SEND(length...) failed',ierr)
     endif
  endif

       ! Step 3.  Send CHARACTER portion of String DummStr 
       ! to process dest.

 call MPI_SEND(DummStr%c, length, MP_CHARACTER, dest, TagBase+1, &
               comm, ierr)
  if(ierr /= 0) then
     if(present(status)) then
        write(stderr,*) myname_,':: MPI_SEND(DummStr%c...) failed.  ierror=',&
	     ierr
        status = ierr
        return
     else
        call MP_perr_die(myname_,':: MPI_SEND(DummStr%c...) failed',ierr)
     endif
  endif

 end subroutine send_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: recv_ - Point-to-point receive of a List variable...
!
! !DESCRIPTION:  This routine receives the output {\tt List} argument 
! {\tt outList} from processor {\tt source} on the communicator associated
! with the fortran 90 {\tt INTEGER} handle {\tt comm}.  The message is
! tagged by the input {\tt INTEGER} argument {\tt TagBase}.  The success
! (failure) of this operation is reported in the zero (nonzero) optional
! output argument {\tt status}.
!
! {\bf N.B.}:  One must avoid assigning elsewhere the MPI tag values 
! {\tt TagBase} and {\tt TagBase+1}.  This is because {\tt recv\_()} 
! performs the receive of the {\tt List} as a pair of operations.  The 
! first receive is the number of characters in {\tt outList\%bf}, and 
! is given MPI tag value {\tt TagBase}.  The second receive is the 
! {\tt CHARACTER} data present in {\tt outList\%bf}, and is given MPI 
! tag value {\tt TagBase+1}.
!
! !INTERFACE:

    subroutine recv_(outList, source, TagBase, comm, status)
!
! !USES:
!

      use m_stdio
      use m_die, only : MP_perr_die

      use m_mpif90

      use m_String, only : String

      implicit none

! !INPUT PARAMETERS: 
!
      integer,            intent(in)  :: source
      integer,            intent(in)  :: TagBase
      integer,            intent(in)  :: comm

! !OUTPUT PARAMETERS: 
!
      type(List),         intent(out) :: outList  
      integer, optional,  intent(out) :: status

! !REVISION HISTORY:
! 	06Jun01 - J.W. Larson - initial version.
! 	11Jun01 - R. Jacob - small bug fix; status in MPI_RECV
!       13Jun01 - J.W. Larson <larson@mcs.anl.gov> - Initialize status
!                 (if present).
!EOP ___________________________________________________________________

 character(len=*),parameter :: myname_=myname//'::recv_'

 integer :: ierr, length
 integer :: MPstatus(MP_STATUS_SIZE)
 type(String) :: DummStr

       ! Initialize status to zero (success), if present.

  if(present(status)) status = 0

       ! Step 1.  Receive Length of String DummStr from process source.

 call MPI_RECV(length, 1, MP_type(length), source, TagBase, comm, &
               MPstatus, ierr)
  if(ierr /= 0) then
     if(present(status)) then
        write(stderr,*) myname_,':: MPI_RECV(length...) failed.  ierror=',&
	     ierr
        status = ierr
        return
     else
        call MP_perr_die(myname_,':: MPI_RECV(length...) failed',ierr)
     endif
  endif

 allocate(DummStr%c(length), stat=ierr)

       ! Step 2.  Send CHARACTER portion of String DummStr 
       ! to process dest.

 call MPI_RECV(DummStr%c, length, MP_CHARACTER, source, TagBase+1, &
               comm, MPstatus, ierr)
  if(ierr /= 0) then
     if(present(status)) then
        write(stderr,*) myname_,':: MPI_RECV(DummStr%c...) failed.  ierror=',&
	     ierr
        status = ierr
        return
     else
        call MP_perr_die(myname_,':: MPI_RECV(DummStr%c...) failed',ierr)
     endif
  endif

       ! Step 3.  Initialize outList.

 call initStr_(outList, DummStr)

 end subroutine recv_

 end module m_List
!.









