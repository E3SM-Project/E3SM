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
      public :: index
      public :: nitem
      public :: get
      public :: assignment(=)

    type List
      character(len=1),dimension(:),pointer :: bf
      integer,       dimension(:,:),pointer :: lc
    end type List

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  interface init ; module procedure	&
	init_,		&
	initStr_,	&
	initstr1_
  end interface
  interface clean; module procedure clean_; end interface
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
  interface assignment(=)
    module procedure copy_
  end interface

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
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nitem_'

  nitem_=size(aList%lc,2)

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
! !IROUTINE: copy_ - Copy a List.  Pointers are copied as allocatables
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine copy_(yL,xL)	! yL=xL
      use m_die,only : die
      use m_mall,only : mall_mci,mall_ison
      implicit none
      type(List),intent(out) :: yL
      type(List),intent(in)  :: xL

! !REVISION HISTORY:
! 	22Apr98 - Jing Guo <guo@thunder> - initial prototype/prolog/code
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::copy_'
  integer :: ln,ni,ier

  ln=len(xL%bf)
  ni=size(xL%lc,2)
  allocate(yL%bf(ln),yL%lc(0:1,ni),stat=ier)
  if(ier /= 0) call die(myname_,'allocate()',ier)

	if(mall_ison()) then
	  call mall_mci(yL%bf,myname)
	  call mall_mci(yL%lc,myname)
	endif

  yL%bf=xL%bf			! string copy
  yL%lc(0:1,:)=xL%lc(0:1,:)	! the locations are relative

	! Note that one may not be able to do this copy easily if
	! a pointer array is used at the place of %lc.  A pointer
	! to an array segment of %bf hides all location information
	! from programmers.  Pointer aliasing can only link the new
	! pointer to the old copy of %bf.  LBOUND() and UBOUND() of
	! a pointer will return only 1 and its size().

end subroutine copy_

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

end module m_List
!.
