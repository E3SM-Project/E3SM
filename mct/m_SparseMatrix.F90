!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_SparseMatrix -- Sparse Matrix class and methods.
!
! !DESCRIPTION:
! The SparseMatrix data type is a special case of the  AttrVect data 
! type (see the module {\tt m\_AttrVect} for more details).  This 
! data type has two storage arrays, one for integer attributes 
! (SparseMatrix\%iAttr) to hold row and column data, and one for 
! real attributes (SparseMatrix\%rAttr) which holds the matrix element 
! for that row and column.
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

      implicit none

      private   ! except

      public :: SparseMatrix     ! The class data structure
      public :: init             ! Create a SparseMatrix
      public :: clean            ! Destroy a SparseMatrix
      public :: local_row_range  ! Local (on-process) row range
      public :: global_row_range ! Local (on-process) row range
      public :: local_col_range  ! Local (on-process) column range
      public :: global_col_range ! Local (on-process) column range
      public :: row_sum          ! Return SparseMatrix row sums
      public :: row_sum_check    ! Check SparseMatrix row sums against
                                 ! input "valid" values

    interface init  ; module procedure init_  ; end interface
    interface clean ; module procedure clean_ ; end interface
    interface local_row_range ; module procedure local_row_range_ ; end interface
    interface global_row_range ; module procedure global_row_range_ ; end interface
    interface local_col_range ; module procedure local_col_range_ ; end interface
    interface global_col_range ; module procedure global_col_range_ ; end interface
    interface row_sum ; module procedure row_sum_ ; end interface
    interface row_sum_check ; module procedure row_sum_check_ ; end interface

! !REVISION HISTORY:
!       19Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!       15Jan01 - J.W. Larson <larson@mcs.anl.gov> - added numerous APIs
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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: local_row_range_ - range of rows with nonzero elements
!
! !DESCRIPTION: This routine examines the input distributed 
! {\tt SparseMatrix} variable {\tt sMat}, and returns local (on-process)
! range of rows having nonzero elements.  The first local row with 
! nonzero elements is returned in the {\tt INTEGER} argument 
! {\tt start\_row}, the last row in {\tt end\_row}.
!
! !INTERFACE:

 subroutine local_row_range_(sMat, start_row, end_row)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_init => init
      use m_die

      implicit none

      type(SparseMatrix), intent(in)  :: sMat
      integer,            intent(out) :: start_row
      integer,            intent(out) :: end_row

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::local_row_range_'

 end subroutine local_row_range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: global_row_range_ - range of rows with nonzero elements
!
! !DESCRIPTION:  This routine examines the input distributed 
! {\tt SparseMatrix} variable {\tt sMat}, and returns global range of 
! rows having nonzero ! elements.  The first local row with nonzero 
! elements is returned in the {\tt INTEGER} argument {\tt start\_row}, 
! the last row in {\tt end\_row}.  The global determination of these 
! values is facilitated by the input communicator handle {\tt comm}.
!
! !INTERFACE:

 subroutine global_row_range_(sMat, comm, start_row, end_row)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_init => init
      use m_die

      implicit none

      type(SparseMatrix), intent(in)  :: sMat
      integer,            intent(in)  :: comm
      integer,            intent(out) :: start_row
      integer,            intent(out) :: end_row

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::global_row_range_'

 end subroutine global_row_range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: local_col_range_ - range of cols with nonzero elements
!
! !DESCRIPTION: This routine examines the input distributed 
! {\tt SparseMatrix} variable {\tt sMat}, and returns local (on-process)
! range of columns having nonzero elements.  The first local column 
! with nonzero elements is returned in the {\tt INTEGER} argument 
! {\tt start\_col}, the last column in {\tt end\_col}.
!
! !INTERFACE:

 subroutine local_col_range_(sMat, start_col, end_col)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_init => init
      use m_die

      implicit none

      type(SparseMatrix), intent(in)  :: sMat
      integer,            intent(out) :: start_col
      integer,            intent(out) :: end_col

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::local_col_range_'

 end subroutine local_col_range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: global_col_range_ - range of cols with nonzero elements
!
! !DESCRIPTION:  This routine examines the input distributed 
! {\tt SparseMatrix} variable {\tt sMat}, and returns global range of 
! columns having nonzero ! elements.  The first global column with 
! nonzero elements is returned in the {\tt INTEGER} argument 
! {\tt start\_col}, the last column in {\tt end\_col}.  The global 
! determination of these values is facilitated by the input communicator 
! handle {\tt comm}.
!
! !INTERFACE:

 subroutine global_col_range_(sMat, comm, start_col, end_col)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_init => init
      use m_die

      implicit none

      type(SparseMatrix), intent(in)  :: sMat
      integer,            intent(in)  :: comm
      integer,            intent(out) :: start_col
      integer,            intent(out) :: end_col

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::global_col_range_'

 end subroutine global_col_range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: row_sum_ - Sum elements in each row of a sparse matrix
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine row_sum_(sMat, sums, start_row, num_rows)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_init => init
      use m_die

      implicit none

      type(SparseMatrix), intent(in)  :: sMat
      real, dimension(:), pointer     :: sums
      integer,            intent(out) :: start_row
      integer,            intent(out) :: num_rows

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::row_sum_'

 end subroutine row_sum_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: row_sum_check - Check row sums of a distributed SparseMatrix.
!
! !DESCRIPTION:  The routine {\tt row\_sum\_check()} sums the rows of 
! the input distributed (across the communicator identified by {\tt comm}) 
! {\tt SparseMatrix} variable {\tt sMat}.  It then compares these sums 
! with the {\tt num\_valid} input "valid" values stored in the array 
! {\tt valid\_sums}.  If all of the sums are within the absolute tolerence
! specified by the input argument {\tt abs\_tol} of any of the valid values,
! the output {\tt LOGICAL} flag {\tt valid} is set to {\tt .TRUE}.  
! Otherwise, this flag is returned with value {\tt .FALSE}.
!
! !INTERFACE:

 subroutine row_sum_check_(sMat, comm, num_valid, valid_sums, abs_tol, valid)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_init => init
      use m_die

      implicit none

      type(SparseMatrix), intent(in)  :: sMat
      integer,            intent(in)  :: comm
      integer,            intent(in)  :: num_valid
      real,               intent(in)  :: valid_sums(num_valid)
      real,               intent(in)  :: abs_tol
      logical,            intent(out) :: valid

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::row_sum_check_'

 end subroutine row_sum_check_

 end module m_SparseMatrix



