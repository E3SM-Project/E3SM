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
      public :: initialized     ! check if initialized
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
	initv_, &
        initp_ 
    end interface
    interface clean  ; module procedure clean_  ; end interface
    interface initialized; module procedure initialized_ ; end interface
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
!        09Aug01 - E.T. Ong <eong@mcs.anl.gov> - added initialized and
!                  initp_ routines. Added 'action' in Accumulator type.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Accumulator'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize an accumulator with given iList, rList, length, 
! num\_steps, and steps\_done.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine init_(aC,iList,iAction,rList,rAction,lsize,num_steps,steps_done)
!
! !USES:
!
      use m_List, only : List_init => init
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_zero => zero
      use m_die

      implicit none

      type(Accumulator),               intent(out) :: aC
      character(len=*),      optional, intent(in)  :: iList
      integer, dimension(:), optional, intent(in)  :: iAction
      character(len=*),      optional, intent(in)  :: rList
      integer, dimension(:), optional, intent(in)  :: rAction
      integer,                         intent(in)  :: lsize
      integer,                         intent(in)  :: num_steps
      integer,               optional, intent(in)  :: steps_done

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
  logical :: status

  action(1) = MCT_SUM
  action(2) = MCT_AVG

        ! Check that aC is not initialized

  status = initialized_(aC,die_flag=.true.,source_name=myname_)

        ! Initialize the Accumulator components.

  call initp_(aC,iAction,rAction,num_steps,steps_done,warning_flag=.true.)

        ! Initialize the AttrVect component aC:

  call AttrVect_init(aC%av,iList,rList,lsize)
  call AttrVect_zero(aC%av)

        ! Check that aC is initialized
  status = initialized_(aC,die_flag=.false.,source_name=myname_)

 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initp_ - initialize an accumulator with given iList, rList, 
!            length, num\_steps, and steps\_done, without initializing 
!            av component of accumulator
!
! !DESCRIPTION: This routine is meant to be used only by member modules, not
!               by the user. The warning flag must be present for the routine
!               to be called correctly. 
!
! !INTERFACE:

 subroutine initp_(aC,iAction,rAction,num_steps,steps_done,warning_flag)
!
! !USES:
!
      use m_die

      implicit none

      type(Accumulator),               intent(out) :: aC
      integer, dimension(:), optional, intent(in)  :: iAction
      integer, dimension(:), optional, intent(in)  :: rAction
      integer,                         intent(in)  :: num_steps
      integer,               optional, intent(in)  :: steps_done
      logical,                         intent(in)  :: warning_flag

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

        ! if the argument steps_done is not present, assume
        ! the accumulator is starting at step zero, that is,
        ! set steps_completed to zero

  steps_completed = 0
  if(present(steps_done)) steps_completed = steps_done

        ! Set the stepping info:

  aC%num_steps = num_steps
  aC%steps_done = steps_completed


        ! Set indexing info

  aC%niAction = 0
  aC%nrAction = 0

        ! Assign iAction and niAction components

  nullify(aC%iAction,aC%rAction)

  if(present(iAction)) then

      allocate(aC%iAction(1:size(iAction)),stat=ier)
      if(ier /= 0) call die(myname_,"iAction allocate",ier)
      
      do i=1,size(iAction)

	 select case (iAction(i))
	 case (MCT_SUM)
	    aC%niAction(MCT_SUM) = aC%niAction(MCT_SUM) + 1
	 case(MCT_AVG)
	    aC%niAction(MCT_AVG) = aC%niAction(MCT_AVG) + 1
	 case default
	    call die(myname_,"illegal iAction assignment")
	 end select

	 ! Safe? pointer copy
	 aC%iAction(i) = iAction(i)

      enddo

   endif

  if(present(rAction)) then

      allocate(aC%rAction(1:size(rAction)),stat=ier)
      if(ier /= 0) call die(myname_,"iAction allocate",ier)

      do i=1,size(rAction)

	 select case (rAction(i))
	 case (MCT_SUM)
	    aC%nrAction(MCT_SUM) = aC%nrAction(MCT_SUM) + 1
	 case(MCT_AVG)
	    aC%nrAction(MCT_AVG) = aC%nrAction(MCT_AVG) + 1
	 case default
	    call die(myname_,"illegal rAction assignment")
	 end select

	 ! Safe? pointer copy
	 aC%rAction(i) = rAction(i)

      enddo
	 
   endif

 end subroutine initp_



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
      integer, optional, intent(in)  :: lsize
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
  integer :: aC_lsize
  integer :: bC_iActions, bC_rActions 
  integer, dimension(:), allocatable :: iActionArray, rActionArray
  integer :: i,ier
  logical :: status

        ! Check aC and bC arguments

  status = initialized(aC,die_flag=.true.,source_name=myname_)
  status = initialized(bC,die_flag=.false.,source_name=myname_)

        ! If the argument steps_done is present, set steps_completed
        ! to this value; otherwise, set it to zero

  steps_completed = 0
  if(present(steps_done)) steps_completed = steps_done 

        ! If the argument lsize is present, 
        ! set aC_lsize to this value; otherwise, set it to the lsize of bC
 
  if(present(lsize)) then 
     aC_lsize = lsize     
  else
     aC_lsize = lsize_(bC)
  endif

	! Convert the two Lists to two Strings

  call List_get(iLStr,bC%av%iList)
  call List_get(rLStr,bC%av%rList)

  bC_iActions = size(bC%iAction)
  bC_rActions = size(bC%rAction)

        ! Convert the pointers to arrays

  allocate(iActionArray(bC_iActions),rActionArray(bC_rActions),stat=ier)
  if(ier /= 0) call die(myname_,"iActionArray/rActionArray allocate",ier)


        ! Call init with present arguments

  if( (bC_iActions > 0) .and. (bC_rActions > 0) ) then

     do i=1,bC_iActions
	iActionArray(i)=bC%iAction(i)
     enddo

     do i=1,bC_rActions
	rActionArray(i)=bC%rAction(i)
     enddo     

     call init_(aC, iList=String_char(iLStr), iAction=iActionArray, &
	        rList=String_char(rLStr), rAction=rActionArray, &
                lsize=aC_lsize, num_steps=num_steps, &
                steps_done=steps_completed)

  else 

     if( bC_iActions > 0 ) then

	do i=1,bC_iActions
	   iActionArray(i)=bC%iAction(i)
	enddo

	call init_(aC, iList=String_char(iLStr), iAction=iActionArray, &
                   lsize=aC_lsize, num_steps=num_steps, &
                   steps_done=steps_completed)
     endif

     if( bC_rActions > 0 ) then

	do i=1,bC_rActions
	   rActionArray(i)=bC%rAction(i)
	enddo

	call init_(aC, rList=String_char(rLStr), rAction=rActionArray, &
                   lsize=aC_lsize, num_steps=num_steps, &
                   steps_done=steps_completed)
     endif

  endif

  deallocate(iActionArray,rActionArray,stat=ier)
  if(ier /= 0) call die(myname_,"iActionArray/rActionArray deallocate",ier)

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

 subroutine clean_(aC,stat)
!
! !USES:
!
      use m_mall
      use m_stdio
      use m_die
      use m_AttrVect, only : AttrVect_clean => clean

      implicit none

      type(Accumulator), intent(inout) :: aC
      integer, optional, intent(out)   :: stat

! !REVISION HISTORY:
! 	11Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!       27JUL01 - E.T. Ong <eong@mcs.anl.gov> - deallocate pointers iAction
!                 and rAction.
!       01Mar02 - E.T. Ong <eong@mcs.anl.gov> removed the die to prevent
!                 crashes and added stat argument.       
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'
  integer :: ier

  if(present(stat)) then
     stat=0
     call AttrVect_clean(aC%av,stat)
  else
     call AttrVect_clean(aC%av)
  endif

  if( associated(aC%iAction) )  then

      deallocate(aC%iAction,stat=ier)

      if(ier /= 0) then
	 if(present(stat)) then
	    stat=ier
	 else
	    call warn(myname_,'deallocate(aC%iAction)',ier)
	 endif
      endif

   endif

  if( associated(aC%rAction) ) then

     deallocate(aC%rAction,stat=ier)

      if(ier /= 0) then
	 if(present(stat)) then
	    stat=ier
	 else
	    call warn(myname_,'deallocate(aC%rAction)',ier)
	 endif
      endif

  endif

 end subroutine clean_


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initialized_ - Check if an Accumulator is initialized.
!
! !DESCRIPTION: If an Accumulator is initialized correctly, the function
!               returns true. Argument die\_flag, if present, will err 
!               if die\_flag is true and aC is correctly initialized, and
!               if die\_flag is false and aC is incorrectly initialized.
!
! !INTERFACE:

 function initialized_(aC, die_flag, source_name)
!
! !USES:
!

   use m_die
   use m_AttrVect, only : AttrVect
   use m_AttrVect, only : Attr_nIAttr => nIAttr
   use m_AttrVect, only : Attr_nRAttr => nRAttr

   implicit none

   logical                                :: initialized_
   type(Accumulator),          intent(in) :: aC
   logical,          optional, intent(in) :: die_flag
   character(len=*), optional, intent(in) :: source_name

! !REVISION HISTORY:
!       07AUG01 - E.T. Ong <eong@mcs.anl.gov> - initital prototype
!
!EOP ___________________________________________________________________

   character(len=*),parameter :: myname_=myname//'::initialized_'
   integer :: i,ier
   logical :: init_kill,uninit_kill
   logical :: aC_associated

   if(present(die_flag)) then
      if(die_flag .eqv. .true.) then
	 init_kill = .true.
	 uninit_kill = .false.
      endif
      if(die_flag .eqv. .false.) then
	 init_kill = .false.
	 uninit_kill = .true.
      endif
   else
      init_kill = .false.
      uninit_kill = .false.
   endif

          ! Initial value
   initialized_ = .true.
   aC_associated = .true.

          ! If any of the pointers in the Accumulator are associated,
          ! then the Accumulator has been initialized

   if( .NOT. (associated(aC%iACtion) .and. associated(aC%rAction) .and. &
              associated(aC%av%iAttr) .and. associated(aC%av%rAttr)) ) then
      initialized_ = .false.
      aC_associated = .false.
      if(uninit_kill) then
	 if(present(source_name)) call perr(source_name,"Accumulator Initialization Error")
	 call die(myname_,"aC pointers are unassociated")
      endif

   endif
   
        ! More sanity checking

   if( aC_associated .eqv. .true. ) then

      if( (Attr_nIAttr(aC%av) == 0) .and. (Attr_nRAttr(aC%av) == 0) ) then
	 initialized_ = .false.
	 if(uninit_kill) then
	    if(present(source_name)) call perr(source_name,"Accumulator Initialization Error")
	    call die(myname_,"No attributes found in aC%av")
	 endif
      endif

      if(Attr_nIAttr(aC%av) > 0) then
	 if( size(aC%iAction) /=  Attr_nIAttr(aC%av) ) then
	    initialized_ = .false.
	    if(uninit_kill) then
	       if(present(source_name)) call perr(source_name,"Accumulator Initialization Error")
	       call die(myname_,"size(aC%iAction) /= nIAttr(aC%av)")
	    endif
	 endif

	 do i=1,Attr_nIAttr(aC%av)
	    if( (aC%iAction(i) /= MCT_SUM) .and. &
                (aC%iAction(i) /= MCT_AVG) ) then
	       initialized_ = .false.
	       if(uninit_kill) then
		  if(present(source_name)) call perr(source_name,"Accumulator Initialization Error")
		  call die(myname_,"Invalid aC%iAction")
	       endif
	    endif
	 enddo

      endif ! if(Attr_nIAttr(aC%av) > 0)

      if(Attr_nRAttr(aC%av) > 0) then

	 if( size(aC%rAction) /=  Attr_nRAttr(aC%av) ) then
	    initialized_ = .false.
	    if(uninit_kill) then
	      if(present(source_name)) call perr(source_name,"Accumulator Initialization Error")
	      call die(myname_,"size(aC%rAction) /= nRAttr(aC%av)")
	    endif
	 endif

	 do i=1,Attr_nRAttr(aC%av)
	    if( (aC%rAction(i) /= MCT_SUM) .and. &
                (aC%rAction(i) /= MCT_AVG) ) then
	       initialized_ = .false.
	       if(uninit_kill) then
		  if(present(source_name)) call perr(source_name,"Accumulator Initialization Error")
		  call die(myname_,"Invalid aC%rAction")
	       endif
	    endif
	 enddo

      endif ! if(Attr_nRAttr(aC%av) > 0)

   endif  ! if (aC_associated .eqv. .true.)

   if(init_kill) then
      if(initialized_ .eqv. .true.) then
	 if(present(source_name)) call perr(source_name,"Accumulator Initialization Error")
	 call die(myname_,"aC has been previously initialized")
      endif
   endif

 end function initialized_


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

