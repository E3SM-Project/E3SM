!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_Router -- Router class 
!
! !DESCRIPTION:
! The Router data type contains all the information needed
! to send an AttrVect from a component on M MPI-processes to a component
! on N MPI-processes.  
!
! !INTERFACE:

 module m_Router
!
! !USES:
!
      use m_Navigator
!     use MPH_module

      implicit none

      private   ! except

      public :: Router	        ! The class data structure
      public :: init            ! Create a Router
      public :: clean           ! Destroy a Router

!  if I am the sending model, then on return pe_list is the
!  processor ranks in the remote model to send to, num_segs
!  is how many segments of my GlobalSegMap to send to each, seg_starts
!  is the starting place for each segment and each destination processor
!  and seg_lengths is the total length of each segment, for each processor
!
!  if I am the receiving model, the on return pe_list is the
!  the processors to receive from, num_segs
!  
    type Router
      integer :: model1	       ! MPH component id of sending model
      integer :: model2	       ! MPH component id of receiving model
      character*4 :: type	       ! 'send' or 'recv'
      integer,dimension(:),pointer :: pe_list ! processor ranks of send/receive
      integer,dimension(:),pointer :: num_segs ! number of segments to send/receive
      integer,dimension(:,:),pointer :: seg_starts ! starting index
      integer,dimension(:,:),pointer :: seg_lengths ! total length

      real,dimension(:,:),pointer :: Rbuffer
      integer,dimension(:,:),pointer :: Ibuffer
      Type(Navigator) ::  buffer_nav
    end type Router


    interface init  ; module procedure init_  ; end interface
    interface clean ; module procedure clean_ ; end interface

! !REVISION HISTORY:
!      15Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_Router'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize a Router
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine init_(model1,model2,GSMap,Rout )
!
! !USES:
!
      use m_GlobalSegMap
      use m_die

      implicit none

      type(Router), intent(out)        :: Rout
      integer, intent(in)	       :: model1
      integer, intent(in)	       :: model2
      type(GlobalSegMap),intent(in)    :: GSMap  ! of the calling model

! !REVISION HISTORY:
!       15Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::init_'

 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Destroy a Router
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(Rout)
!
! !USES:
!

      implicit none

      type(Router), intent(inout) :: Rout

! !REVISION HISTORY:
!       15Jan01 - R. Jacob <jacob@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'


 end subroutine clean_

 end module m_Router

