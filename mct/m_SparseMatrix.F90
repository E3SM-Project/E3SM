!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_SparseMatrix -- Sparse Matrix class and methods.
!
! !DESCRIPTION:
! The SparseMatrix data type is a special case of the  AttrVect data 
! type (see m_AttrVect for details).  This data type has two storage 
! arrays, one for integer attributes (SparseMatrix\%iAttr) to hold row
! and column data, and one for real attributes (symMatx\%rAttr) which 
! holds the matrix element for that row and column.
!
! The set of attributes in each storage array are defined by a List; 
! SparseMatrix\%iList for integer attributes, and SparseMatrix\%rList 
! for real attributes.  The integer and real attribute tags are defined
! below:
!
! SparseMatrix\%iList components:
!    row : row index
!    col : column index
!
! SparseMatrix\%rList components:
!    weight : matrix element
!
! !INTERFACE:

 module m_SparseMatrix
!
! !USES:
!
      use m_AttrVect, only : SparseMatrix => AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVectComms, only : gather
      use m_AttrVectComms, only : scatter
      use m_AttrVectComms, only : bcast

      implicit none

      private   ! except

      public :: SparseMatrix    ! The class data structure
      public :: init            ! Create a SparseMatrix
      public :: clean           ! Destroy a SparseMatrix
      public :: gather          ! Gather a SparseMatrix
      public :: scatter         ! Scatter a SparseMatrix
      public :: bcast           ! Broadcast a SparseMatrix

    interface init  ; module procedure init_  ; end interface
    interface clean ; module procedure clean_ ; end interface

! !REVISION HISTORY:
!       19Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_SparseMatrix'

! SparseMatrix_iList components:
  character(len=*),parameter :: SparseMatrix_iList='row:col'
  integer,parameter :: SparseMatrix_irow=1
  integer,parameter :: SparseMatrix_icol=2

! SparseMatrix_rList components:
  character(len=*),parameter :: SparseMatrix_rList='weight'
  integer,parameter :: SparseMatrix_iweight=1

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: init_ - initialize a SparseMatrix
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine init_(sMat, lsize )
!
! !USES:
!
      use m_AttrVect, only : AttrVect_init => init
      use m_die

      implicit none

      type(SparseMatrix), intent(out)        :: sMat
      integer,         optional,intent(in)   :: lsize

! !REVISION HISTORY:
!       19Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::init_'

  integer :: n

        ! if lsize is present, use it to set n; if not, set n=0

  n = 0
  if(present(lsize)) n=lsize


        ! Initialize sMat using AttrVect_init

  call AttrVect_init(sMat, SparseMatrix_iList, SparseMatrix_rList, n)

 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Destroy a SparseMatrix.
!
! !DESCRIPTION:
!
! !INTERFACE:

    subroutine clean_(sMat)
!
! !USES:
!
      use m_AttrVect,only : AttrVect_clean => clean

      implicit none

      type(SparseMatrix), intent(inout) :: sMat

! !REVISION HISTORY:
!       19Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

  call AttrVect_clean(sMat)

 end subroutine clean_

 end module m_SparseMatrix

