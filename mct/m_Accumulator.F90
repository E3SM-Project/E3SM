!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Accumulator - a distributed accumulator for time-averaging.
!
! !DESCRIPTION:
!
! An {\em accumulator} is a data class used for computing running sums 
! and/or time averages of {\tt AttrVect} class data.  The fortran 90 
! implementation of this concept is the {\tt Accumulator} class, 
! which---along with its basic methods---is defined in this module.
!
! !INTERFACE:

 module m_Accumulator
!
! !USES:
!
      use m_List, only : List
      use m_AttrVect, only : AttrVect

      implicit none

      private	! except

! The class data structure

      public :: Accumulator     

! Defined constants

      public :: MCT_SUM
      public :: MCT_AVG


! List of Methods for the Accumulator class

      public :: init		! creation method
      public :: clean		! destruction method
      public :: lsize		! local length of the data arrays
      public :: nIAttr		! number of integer fields
      public :: nRAttr		! number of real fields
      public :: indexIA		! index the integer fields
      public :: indexRA		! index the real fields
      public :: getIList	! Return tag from INTEGER 
                                ! attribute list
      public :: getRList	! Return tag from REAL attribute
                                ! list

! Definition of the Accumulator class:

    type Accumulator
      integer :: num_steps      ! total number of accumulation steps
      integer :: steps_done     ! number of accumulation steps performed
      integer :: niAction(1:2)  ! number of integer actions
      integer :: nrAction(1:2)  ! number of real actions
      integer, pointer, dimension(:) :: iAction ! index of integer actions
      integer, pointer, dimension(:) :: rAction ! index of real actions
      type(AttrVect) :: av      ! accumulated sum field storage
      
    end type Accumulator

! Assignment of constants

    integer, parameter :: MCT_SUM = 1
    integer, parameter :: MCT_AVG = 2

! Definition of interfaces for the methods for the Accumulator:

    interface init   ; module procedure	&
	init_,	&
	initv_
    end interface
    interface clean  ; module procedure clean_  ; end interface
    interface lsize  ; module procedure lsize_  ; end interface
    interface nIAttr ; module procedure nIAttr_ ; end interface
    interface nRAttr ; module procedure nRAttr_ ; end interface
    interface indexIA; module procedure indexIA_; end interface
    interface indexRA; module procedure indexRA_; end interface
    interface getIList; module procedure getIList_; end interface
    interface getRList; module procedure getRList_; end interface

! !REVISION HISTORY:
! 	 7Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
! 	 7Feb01 - Jay Larson <larson@mcs.anl.gov> - Public interfaces
!                 to getIList() and getRList().
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Accumulator'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize with given iList, rList, length, 
! num\_steps, and steps\_done.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine init_(aC,iList,iAction,rList,rAction,lsize,num_steps,steps_done)
!
! !USES:
!
      use m_List, only : List_init=>init
      use m_List, only : List_nitem=>nitem
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_die

      implicit none

      type(Accumulator),intent(out)        :: aC
      character(len=*),optional,intent(in) :: iList
      integer,dimension(:),optional,intent(in) :: iaction
      character(len=*),optional,intent(in) :: rList
      integer,dimension(:),optional,intent(in) :: raction
      integer,         intent(in) :: lsize
      integer,         intent(in)          :: num_steps
      integer,         optional,intent(in) :: steps_done

! !REVISION HISTORY:
! 	11Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!       27JUL01 - E.T. Ong <eong@mcs.anl.gov> - added iAction, rAction,
!                 niAction, and nrAction to accumulator type. Also defined
!                 MCT_SUM and MCT_AVG for accumulator module.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::init_'
  integer :: n,i,ier
  integer :: steps_completed
  integer :: action(2)

  action(1) = MCT_SUM
  action(2) = MCT_AVG

        ! Argument dummy check
  
  ! Bug note: On denali, this portion will kill before MP_perr_die can execute
  if(.not.( (present(iList).and.present(iAction)) .or. &
            (present(rList).and.present(rAction)) ) ) then
     call MP_perr_die(myname_,"List and Action arguments must be paired",ier)
  endif

  if(lsize .le. 0) then
     call MP_perr_die(myname_,"lsize argument must be > 0",ier) 
  endif

        ! if the argument steps_done is not present, assume
        ! the accumulator is starting at step zero, that is,
        ! set st_complete to zero

  steps_completed = 0
  if(present(steps_done)) steps_completed = steps_done

        ! Set the stepping info:

  aC%num_steps = num_steps
  aC%steps_done = steps_completed


        ! Initialize the AttrVect component aC%av:


  call AttrVect_init(aC%av,iList,rList,lsize)



        ! Set indexing info

  aC%niAction = 0
  aC%nrAction = 0

  nullify(aC%iAction,aC%rAction)

  if(present(iAction)) then

      ! More argument checking

      if( size(iAction) /=  AttrVect_nIAttr(aC%av) ) then
	 call MP_perr_die(myname_,"size(iaction) /= size(iList)",ier)
      endif

      allocate(aC%iAction(1:AttrVect_nIAttr(aC%av)),stat=ier)
      if(ier /= 0) then
           call MP_perr_die(myname_,"iAction allocate",ier)
      endif
      
      do i=1,AttrVect_nIAttr(aC%av)

	 select case (iAction(i))
	 case (MCT_SUM)
	    aC%niAction(MCT_SUM) = aC%niAction(MCT_SUM) + 1
	 case(MCT_AVG)
	    aC%niAction(MCT_AVG) = aC%niAction(MCT_AVG) + 1
	 case default
	    call MP_perr_die(myname_,"illegal iAction assignment",ier)
	 end select

	 ! Safe? pointer copy
	 aC%iAction(i) = iAction(i)

      enddo

   endif

  if(present(rAction)) then

      ! More argument checking

      if( size(rAction) .ne.  AttrVect_nRAttr(aC%av) ) then
	 call MP_perr_die(myname_,"size(raction) /= size(rList)",ier)
      endif

      allocate(aC%rAction(1:AttrVect_nRAttr(aC%av)),stat=ier)
      if(ier /= 0) call MP_perr_die(myname_,"iAction allocate",ier)

      do i=1,AttrVect_nRAttr(aC%av)

	 select case (rAction(i))
	 case (MCT_SUM)
	    aC%nrAction(MCT_SUM) = aC%nrAction(MCT_SUM) + 1
	 case(MCT_AVG)
	    aC%nrAction(MCT_AVG) = aC%nrAction(MCT_AVG) + 1
	 case default
	    call MP_perr_die(myname_,"illegal rAction assignment",ier)
	 end select

	 ! Safe? pointer copy
	 aC%rAction(i) = rAction(i)

      enddo
	 
   endif

 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initv_-Initialize an accumulator using another Accumulator.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine initv_(aC,bC,lsize,num_steps,steps_done)
!
! !USES:
!
      use m_String, only : String
      use m_String, only : String_char => char
      use m_List,   only : List_get => get
      use m_die

      implicit none

      type(Accumulator),    intent(out) :: aC
      type(Accumulator),    intent(in)  :: bC
      integer,           intent(in)  :: lsize
      integer,           intent(in)  :: num_steps
      integer, optional, intent(in)  :: steps_done

! !REVISION HISTORY:
! 	11Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
! 	17May01 - R. Jacob <jacob@mcs.anl.gov> - change string_get to
!                 list_get
!       27JUL01 - E.T. Ong <eong@mcs.anl.gov> - added iaction,raction 
!                 compatibility
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::initv_'

  type(String) :: iLStr,rLStr
  integer :: steps_completed
  integer :: bC_iActions, bC_rActions 
  integer, dimension(:), allocatable :: iActionArray, rActionArray
  integer :: ier

        ! If the argument steps_done is present, set steps_completed
        ! to this value; otherwise, set it to zero

  steps_completed = 0
  if(present(steps_done)) steps_completed = steps_done 

	! Convert the two Lists to two Strings

  call List_get(iLStr,bC%av%iList)
  call List_get(rLStr,bC%av%rList)

  bC_iActions = size(bC%iAction)
  bC_rActions = size(bC%rAction)

        ! Convert the pointers to arrays

  allocate(iActionArray(bC_iActions),rActionArray(bC_rActions),stat=ier)
  if(ier /= 0) call MP_perr_die(myname_,"iActionArray/rActionArray allocate",ier)


        ! Call init with present arguments

  if( (bC_iActions > 0) .and. (bC_rActions > 0) ) then

     call init_(aC, iList=String_char(iLStr), iAction=bC%iAction, &
                rList=String_char(rLStr), rAction=bC%rAction, &
                lsize=lsize, num_steps=num_steps, steps_done=steps_completed)

  else 

     if( bC_iActions > 0 ) then
	call init_(aC, iList=String_char(iLStr), iAction=bC%iAction, &
                  lsize=lsize, num_steps=num_steps, steps_done=steps_completed)
     endif

     if( bC_rActions > 0 ) then
	call init_(aC, rList=String_char(rLStr), rAction=bC%rAction, &
                  lsize=lsize, num_steps=num_steps, steps_done=steps_completed)
     endif

  endif

  deallocate(iActionArray,rActionArray,stat=ier)
  if(ier /= 0) call MP_perr_die(myname_,"iActionArray/rActionArray deallocate",ier)

 end subroutine initv_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Destruction method for the Accumulator.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine clean_(aC)
!
! !USES:
!
      use m_mall
      use m_stdio
      use m_die
      use m_AttrVect, only : AttrVect_clean => clean

      implicit none

      type(Accumulator),intent(inout) :: aC

! !REVISION HISTORY:
! 	11Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!       27JUL01 - E.T. Ong <eong@mcs.anl.gov> - deallocate pointers iAction
!                 and rAction.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  call AttrVect_clean(aC%av)
  
  if( associated(aC%iAction) )  deallocate(aC%iAction,stat=ier)
  if(ier /= 0) then
     call MP_perr_die(myname_,"iAction deallocate",ier)
  endif

  if( associated(aC%rAction) ) deallocate(aC%rAction,stat=ier)
  if(ier /= 0) then
     call MP_perr_die(myname_,"rAction deallocate",ier)
  endif

 end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - local size of data storage in the Accumulator.
!
! !DESCRIPTION:
!
! !INTERFACE:

 function lsize_(aC)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_lsize => lsize

      implicit none

      type(Accumulator), intent(in) :: aC
      integer :: lsize_

! !REVISION HISTORY:
! 	12Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'


	! The function AttrVect_lsize is called to return
        ! its local size data

  lsize_=AttrVect_lsize(aC%aV)

 end function lsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nIAttr_ - number of INTEGER fields stored in the Accumulator.
!
! !DESCRIPTION:
!
! !INTERFACE:

    function nIAttr_(aC)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr

      implicit none

      type(Accumulator),intent(in) :: aC
      integer :: nIAttr_

! !REVISION HISTORY:
! 	12Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nIAttr_'

	! The function AttrVect_nIAttr is called to return the
        ! number of integer fields

  nIAttr_=AttrVect_nIAttr(aC%av)

 end function nIAttr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nRAttr_ - number of REAL fields stored in the Accumulator.
!
! !DESCRIPTION:
!
! !INTERFACE:

 function nRAttr_(aC)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
 
      implicit none

      type(Accumulator),intent(in) :: aC
      integer :: nRAttr_

! !REVISION HISTORY:
! 	12Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nRAttr_'

	! The function AttrVect_nRAttr is called to return the
        ! number of real fields

  nRAttr_=AttrVect_nRAttr(aC%aV)

 end function nRAttr_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getIList_ - get an item from the integer attribute list 
! in the accumulator data storage area (i.e. its AttrVect component).
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getIList_(item,ith,aC)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_getIList => getIList
      use m_String,   only : String

      implicit none
      type(String),intent(out)     :: item
      integer,     intent(in)      :: ith
      type(Accumulator),intent(in) :: aC

! !REVISION HISTORY:
! 	12Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getIList_'

  call AttrVect_getIList(item,ith,aC%av)

 end subroutine getIList_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: getRList_ - get an item from real attribute list in the
! accumulator data storage space (i.e. its AttrVect component).
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine getRList_(item,ith,aC)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_getRList => getRList
      use m_String,   only : String

      implicit none
      type(String),     intent(out) :: item
      integer,          intent(in)  :: ith
      type(Accumulator),intent(in)  :: aC

! !REVISION HISTORY:
! 	12Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::getRList_'

  call AttrVect_getRList(item,ith,aC%av)

 end subroutine getRList_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexIA_ - index the Accumulator integer attribute List.
!
! !DESCRIPTION:
!
! !INTERFACE:

 function indexIA_(aC,item,perrWith,dieWith)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_indexIA => indexIA
      use m_die,  only : die
      use m_stdio,only : stderr

      implicit none
 
     type(Accumulator), intent(in) :: aC
      character(len=*),intent(in) :: item
      character(len=*),optional,intent(in) :: perrWith
      character(len=*),optional,intent(in) :: dieWith
      integer :: indexIA_

! !REVISION HISTORY:
! 	14Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexIA_'

  indexIA_=AttrVect_indexIA(aC%aV,item)

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
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexRA_ - index the Accumulator real attribute list.
!
! !DESCRIPTION:
!
! !INTERFACE:

 function indexRA_(aC,item,perrWith,dieWith)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_die,  only : die
      use m_stdio,only : stderr

      implicit none

      type(Accumulator), intent(in) :: aC
      character(len=*),intent(in) :: item
      character(len=*),optional,intent(in) :: perrWith
      character(len=*),optional,intent(in) :: dieWith
      integer :: indexRA_

! !REVISION HISTORY:
! 	14Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexRA_'

  indexRA_=AttrVect_indexRA(aC%aV,item)

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

 end module m_Accumulator

