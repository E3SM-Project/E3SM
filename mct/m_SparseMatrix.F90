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
!    grow : global row index
!    gcol : global column index
!    lrow : local row index
!    lcol : local column index
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
!       25Feb01 - J.W. Larson <larson@mcs.anl.gov> - changed from row/column
!                 attributes to global and local row and column attributes
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_SparseMatrix'

! SparseMatrix_iList components:
  character(len=*),parameter :: SparseMatrix_iList='grow:gcol:lrow:lcol'
  integer,parameter :: SparseMatrix_igrow=1
  integer,parameter :: SparseMatrix_igcol=2
  integer,parameter :: SparseMatrix_ilrow=3
  integer,parameter :: SparseMatrix_ilcol=4

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
! {\tt SparseMatrix} variable {\tt sMat}, and returns the range of local 
! row values having nonzero elements.  The first local row with 
! nonzero elements is returned in the {\tt INTEGER} argument 
! {\tt start\_row}, the last row in {\tt end\_row}.
!
! !INTERFACE:

 subroutine local_row_range_(sMat, start_row, end_row)
!
! !USES:
!
      use m_die

      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_indexIA => indexIA

      implicit none

      type(SparseMatrix), intent(in)  :: sMat
      integer,            intent(out) :: start_row
      integer,            intent(out) :: end_row

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!       25Feb01 - Jay Larson <larson@mcs.anl.gov> - Initial prototype.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::local_row_range_'

  integer :: i, ilrow, lsize

  ilrow = AttrVect_indexIA(sMat, 'lrow')
  lsize = AttrVect_lsize(sMat)

       ! Initialize start_row and end_row:

  start_row = sMat%iAttr(ilrow,1)
  end_row = sMat%iAttr(ilrow,1)

  do i=1,lsize
     start_row = min(start_row, sMat%iAttr(ilrow,i))
     end_row = max(end_row, sMat%iAttr(ilrow,i))
  end do

 end subroutine local_row_range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: global_row_range_ - range of rows with nonzero elements
!
! !DESCRIPTION:  This routine examines the input distributed 
! {\tt SparseMatrix} variable {\tt sMat}, and returns the range of 
! global row values having nonzero elements.  The first local row with 
! nonzero elements is returned in the {\tt INTEGER} argument 
! {\tt start\_row}, the last row in {\tt end\_row}. 
!
! !INTERFACE:

 subroutine global_row_range_(sMat, comm, start_row, end_row)
!
! !USES:
!
      use m_die

      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_indexIA => indexIA

      implicit none

      type(SparseMatrix), intent(in)  :: sMat
      integer,            intent(in)  :: comm
      integer,            intent(out) :: start_row
      integer,            intent(out) :: end_row

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!       25Feb01 - Jay Larson <larson@mcs.anl.gov> - Initial prototype.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::global_row_range_'

  integer :: i, igrow, lsize

  igrow = AttrVect_indexIA(sMat, 'grow')
  lsize = AttrVect_lsize(sMat)

       ! Initialize start_row and end_row:

  start_row = sMat%iAttr(igrow,1)
  end_row = sMat%iAttr(igrow,1)

  do i=1,lsize
     start_row = min(start_row, sMat%iAttr(igrow,i))
     end_row = max(end_row, sMat%iAttr(igrow,i))
  end do

 end subroutine global_row_range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: local_col_range_ - range of cols with nonzero elements
!
! !DESCRIPTION: This routine examines the input distributed 
! {\tt SparseMatrix} variable {\tt sMat}, and returns the range of
! local column values having nonzero elements.  The first local column 
! with nonzero elements is returned in the {\tt INTEGER} argument 
! {\tt start\_col}, the last column in {\tt end\_col}.
!
! !INTERFACE:

 subroutine local_col_range_(sMat, start_col, end_col)
!
! !USES:
!
      use m_die

      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_indexIA => indexIA

      implicit none

      type(SparseMatrix), intent(in)  :: sMat
      integer,            intent(out) :: start_col
      integer,            intent(out) :: end_col

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!       25Feb01 - Jay Larson <larson@mcs.anl.gov> - Initial prototype.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::local_col_range_'

  integer :: i, ilcol, lsize

  ilcol = AttrVect_indexIA(sMat, 'lcol')
  lsize = AttrVect_lsize(sMat)

       ! Initialize start_col and end_col:

  start_col = sMat%iAttr(ilcol,1)
  end_col = sMat%iAttr(ilcol,1)

  do i=1,lsize
     start_col = min(start_col, sMat%iAttr(ilcol,i))
     end_col = max(end_col, sMat%iAttr(ilcol,i))
  end do

 end subroutine local_col_range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: global_col_range_ - range of cols with nonzero elements
!
! !DESCRIPTION:  This routine examines the input distributed 
! {\tt SparseMatrix} variable {\tt sMat}, and returns the range of 
! global column values having nonzero elements.  The first global 
! column with nonzero elements is returned in the {\tt INTEGER} argument 
! {\tt start\_col}, the last column in {\tt end\_col}.  
!
! !INTERFACE:

 subroutine global_col_range_(sMat, comm, start_col, end_col)
!
! !USES:
!
      use m_die

      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_indexIA => indexIA

      implicit none

      type(SparseMatrix), intent(in)  :: sMat
      integer,            intent(in)  :: comm
      integer,            intent(out) :: start_col
      integer,            intent(out) :: end_col

! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!       25Feb01 - Jay Larson <larson@mcs.anl.gov> - Initial prototype.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::global_col_range_'

  integer :: i, igcol, lsize

  igcol = AttrVect_indexIA(sMat, 'lcol')
  lsize = AttrVect_lsize(sMat)

       ! Initialize start_col and end_col:

  start_col = sMat%iAttr(igcol,1)
  end_col = sMat%iAttr(igcol,1)

  do i=1,lsize
     start_col = min(start_col, sMat%iAttr(igcol,i))
     end_col = max(end_col, sMat%iAttr(igcol,i))
  end do

 end subroutine global_col_range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: row_sum_ - Sum elements in each row of a sparse matrix
!
! !DESCRIPTION:
! Given an input {\tt SparseMatrix} argument {\tt sMat}, {\tt row\_sum\_()}
! returns the number of the rows {\tt num\_rows} in the sparse matrix and
! the sum of the elements in each row in the array {\tt sums}.  The input
! argument {\tt comm} is the Fortran 90 MPI communicator handle used to
! determine the number of rows and perform the sums.  The output arguments
! {\tt num\_rows} and {\tt sums} are valid on all processes.
!
! {\bf N.B.:  } This routine allocates an array {\tt sums}.  The user is
! responsible for deallocating this array when it is no longer needed.  
! Failure to do so will cause a memory leak.
!
! !INTERFACE:

 subroutine row_sum_(sMat, num_rows, sums, comm)

!
! !USES:
!
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_indexIA => indexIA
      use m_AttrVect, only : AttrVect_indexRA => indexRA

      implicit none

      type(SparseMatrix), intent(in)  :: sMat
      integer,            intent(out) :: num_rows
      real, dimension(:), pointer     :: sums
      integer,            intent(in)  :: comm


! !REVISION HISTORY:
!       15Jan01 - Jay Larson <larson@mcs.anl.gov> - API specification.
!       25Jan01 - Jay Larson <larson@mcs.anl.gov> - Prototype code.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::row_sum_'

  integer :: i, igrow, ierr, iwgt, lsize, myID
  integer :: start_row, end_row
  real, dimension(:), allocatable :: lsums

       ! Determine local rank

  call MP_COMM_RANK(comm, myID, ierr)

       ! Determine on each process the row of global row indices:

  call global_row_range_(sMat, comm, start_row, end_row)

       ! Determine across the communicator the _maximum_ value of
       ! end_row, which will be assigned to num_rows on each process:

  call MP_ALLREDUCE(end_row, num_rows, 1, MP_INTEGER, MP_MAX, &
                    comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,"MP_ALLREDUCE(end_row...",ierr)
  endif

       ! Allocate storage for the sums on each process.

  allocate(lsums(num_rows), sums(num_rows), stat=ierr)

  if(ierr /= 0) then
     call MP_perr_die(myname_,"allocate(lsums(...",ierr)
  endif

       ! Compute the local entries to lsum(1:num_rows) for each process:

  lsize = AttrVect_lsize(sMat)
  igrow = AttrVect_indexIA(sMat,'grow')
  iwgt = AttrVect_indexRA(sMat,'weight')

  lsums = 0.
  do i=1,lsize
     lsums(sMat%iAttr(igrow,i)) = lsums(sMat%iAttr(igrow,i)) + &
	                           sMat%rAttr(iwgt,i)
  end do

       ! Compute the global sum of the entries of lsums so that all
       ! processes own the global sums.

  call MP_ALLREDUCE(lsums, sums, num_rows, MP_REAL, MP_SUM, comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,"MP_ALLREDUCE(lsums...",ierr)
  endif

       ! Clean up...

  deallocate(lsums, stat=ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,"deallocate(lsums...",ierr)
  endif

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
!       25Feb01 - Jay Larson <larson@mcs.anl.gov> - Prototype code.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::row_sum_check_'

  integer :: i, j, num_invalid, num_rows
  real, dimension(:), pointer :: sums

       ! Compute row sums:

  call row_sum_(sMat, num_rows, sums, comm)

       ! Initialize for the scanning loop (assume the matrix row
       ! sums are valid):

  valid = .TRUE.
  i = 1

  SCAN_LOOP:  do 

       ! Count the number of elements in valid_sums(:) that
       ! are separated from sums(i) by more than abs_tol

     num_invalid = 0

     do j=1,num_valid
	if(abs(sums(i) - valid_sums(j)) > abs_tol) then
	   num_invalid = num_invalid + 1
	endif
     end do

       ! If num_invalid = num_valid, then we have failed to
       ! find a valid sum value within abs_tol of sums(i).  This
       ! one failure is enough to halt the process.

     if(num_invalid == num_valid) then
	valid = .FALSE.
	EXIT
     endif

       ! Prepare index i for the next element of sums(:)

     i = i + 1
     if( i > num_rows) EXIT

  end do SCAN_LOOP

 end subroutine row_sum_check_

 end module m_SparseMatrix



