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
! SparseMatrix\%data\%iList components:
!    grow : global row index
!    gcol : global column index
!    lrow : local row index
!    lcol : local column index
!
! SparseMatrix\%data\%rList components:
!    weight : matrix element
!
! !INTERFACE:

 module m_SparseMatrix
!
! !USES:
!
      use m_AttrVect, only : AttrVect

      private   ! except

! !PUBLIC TYPES:

      public :: SparseMatrix      ! The class data structure

      Type SparseMatrix
	 integer :: nrows
	 integer :: ncols
	 type(AttrVect) :: data
      End Type SparseMatrix

! !PUBLIC MEMBER FUNCTIONS:

      public :: init              ! Create a SparseMatrix
      public :: clean             ! Destroy a SparseMatrix
      public :: lsize             ! Local number of elements
      public :: indexIA           ! Index integer attribute
      public :: indexRA           ! Index real attribute
      public :: nRows             ! Total number of rows
      public :: nCols             ! Total number of columns

      public :: exportGlobalRowIndices    ! Return global row indices 
                                          ! for matrix elements
      public :: exportGlobalColumnIndices ! Return global column indices 
                                          ! for matrix elements
      public :: exportLocalRowIndices     ! Return local row indices 
                                          ! for matrix elements
      public :: exportLocalColumnIndices  ! Return local column indices 
                                          ! for matrix elements
      public :: exportMatrixElements      ! Return matrix elements

      public :: importGlobalRowIndices    ! Set global row indices 
                                          ! using 
      public :: importGlobalColumnIndices ! Return global column indices 
                                          ! for matrix elements
      public :: importLocalRowIndices     ! Return local row indices 
                                          ! for matrix elements
      public :: importLocalColumnIndices  ! Return local column indices 
                                          ! for matrix elements
      public :: importMatrixElements      ! Return matrix elements

      public :: GlobalNumElements ! Total number of nonzero elements
      public :: ComputeSparsity   ! Fraction of matrix that is nonzero
      public :: local_row_range   ! Local (on-process) row range
      public :: global_row_range  ! Local (on-process) row range
      public :: local_col_range   ! Local (on-process) column range
      public :: global_col_range  ! Local (on-process) column range
      public :: CheckBounds       ! Check row and column values
                                  ! for out-of-bounds values
      public :: row_sum           ! Return SparseMatrix row sums
      public :: row_sum_check     ! Check SparseMatrix row sums against
                                  ! input "valid" values
      public :: Sort              ! Sort matrix entries to generate an
                                  ! index permutation (to be used by
                                  ! Permute()
      public :: Permute           ! Permute matrix entries using index
                                  ! permutation gernerated by Sort()
      public :: SortPermute       ! Sort/Permute matrix entries

    interface init  ; module procedure init_  ; end interface
    interface clean ; module procedure clean_ ; end interface
    interface lsize ; module procedure lsize_ ; end interface
    interface indexIA ; module procedure indexIA_ ; end interface
    interface indexRA ; module procedure indexRA_ ; end interface
    interface nRows ; module procedure nRows_ ; end interface
    interface nCols ; module procedure nCols_ ; end interface

    interface exportGlobalRowIndices ; module procedure &
	 exportGlobalRowIndices_ 
    end interface

    interface exportGlobalColumnIndices ; module procedure &
	 exportGlobalColumnIndices_ 
    end interface

    interface exportLocalRowIndices ; module procedure &
	 exportLocalRowIndices_ 
    end interface

    interface exportLocalColumnIndices ; module procedure &
	 exportLocalColumnIndices_ 
    end interface

    interface exportMatrixElements ; module procedure &
	 exportMatrixElements_ 
    end interface

    interface importGlobalRowIndices ; module procedure &
	 importGlobalRowIndices_ 
    end interface

    interface importGlobalColumnIndices ; module procedure &
	 importGlobalColumnIndices_ 
    end interface

    interface importLocalRowIndices ; module procedure &
	 importLocalRowIndices_ 
    end interface

    interface importLocalColumnIndices ; module procedure &
	 importLocalColumnIndices_ 
    end interface

    interface importMatrixElements ; module procedure &
	 importMatrixElements_ 
    end interface

    interface GlobalNumElements ; module procedure &
	 GlobalNumElements_ 
    end interface

    interface ComputeSparsity ; module procedure &
	 ComputeSparsity_ 
    end interface

    interface local_row_range ; module procedure &
	 local_row_range_ 
    end interface

    interface global_row_range ; module procedure &
	 global_row_range_ 
    end interface

    interface local_col_range ; module procedure &
	 local_col_range_ 
    end interface

    interface global_col_range ; module procedure &
	 global_col_range_ 
    end interface

    interface CheckBounds; module procedure &
	 CheckBounds_ 
    end interface

    interface row_sum ; module procedure row_sum_ ; end interface

    interface row_sum_check ; module procedure &
	 row_sum_check_ 
    end interface

    interface Sort ; module procedure Sort_ ; end interface
    interface Permute ; module procedure Permute_ ; end interface
    interface SortPermute ; module procedure SortPermute_ ; end interface

! !REVISION HISTORY:
!       19Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!       15Jan01 - J.W. Larson <larson@mcs.anl.gov> - added numerous APIs
!       25Feb01 - J.W. Larson <larson@mcs.anl.gov> - changed from row/column
!                 attributes to global and local row and column attributes
!       23Apr01 - J.W. Larson <larson@mcs.anl.gov> - added number of rows
!                 and columns to the SparseMatrix type.  This means the
!                 SparseMatrix is no longer a straight AttrVect type.  This
!                 also made necessary the addition of lsize(), indexIA(),
!                 and indexRA().
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
! !DESCRIPTION:  This routine creates the storage space for the
! entries of a {\tt SparseMatrix}, and sets the number of rows and
! columns in it.  The input {\tt INTEGER} arguments {\tt nrows} and 
! {\tt ncols} specify the number of rows and columns respectively.
! The optional input argument {\tt lsize} specifies the number of 
! nonzero entries in the {\tt SparseMatrix}.  The initialized 
! {\tt SparseMatrix} is returned in the output argument {\tt sMat}.
!
! {\bf N.B.}:  This routine is allocating dynamical memory in the form
! of a {\tt SparseMatrix}.  The user must deallocate this space when
! the {\tt SparseMatrix} is no longer needed by invoking the routine
! {\tt clean\_()}.
!
! !INTERFACE:

 subroutine init_(sMat, nrows, ncols, lsize)
!
! !USES:
!
      use m_AttrVect, only : AttrVect_init => init
      use m_die

      implicit none

      type(SparseMatrix), intent(out)  :: sMat
      integer,            intent(in)   :: nrows
      integer,            intent(in)   :: ncols
      integer, optional,  intent(in)   :: lsize

! !REVISION HISTORY:
!       19Sep00 - Jay Larson <larson@mcs.anl.gov> - initial prototype
!       23Apr01 - Jay Larson <larson@mcs.anl.gov> - added arguments
!                 nrows and ncols--number of rows and columns in the
!                 SparseMatrix
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::init_'

  integer :: n

        ! if lsize is present, use it to set n; if not, set n=0

  n = 0
  if(present(lsize)) n=lsize

        ! Initialize number of rows and columns:

  sMat%nrows = nrows
  sMat%ncols = ncols

        ! Initialize sMat%data using AttrVect_init

  call AttrVect_init(sMat%data, SparseMatrix_iList, &
                     SparseMatrix_rList, n)

 end subroutine init_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Destroy a SparseMatrix.
!
! !DESCRIPTION:  This routine deallocates dynamical memory held by the
! input {\tt SparseMatrix} argument {\tt sMat}.  It also sets the number
! of rows and columns in the {\tt SparseMatrix} to zero.
!
! !INTERFACE:

    subroutine clean_(sMat,stat)
!
! !USES:
!
      use m_AttrVect,only : AttrVect_clean => clean

      implicit none

      type(SparseMatrix), intent(inout) :: sMat
      integer, optional,  intent(out)   :: stat

! !REVISION HISTORY:
!       19Sep00 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!       23Apr00 - J.W. Larson <larson@mcs.anl.gov> - added changes to
!                 accomodate clearing nrows and ncols.
!       01Mar02 - E.T. Ong <eong@mcs.anl.gov> Added stat argument.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::clean_'

       ! Deallocate memory held by sMat:

  if(present(stat)) then
     call AttrVect_clean(sMat%data,stat)
  else
     call AttrVect_clean(sMat%data)
  endif

       ! Set the number of rows and columns in sMat to zero:

  sMat%nrows = 0
  sMat%ncols = 0

 end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: lsize_ - Local number of elements in a SparseMatrix.
!
! !DESCRIPTION:  This {\tt INTEGER} function reports on-processor storage 
! of the number of nonzero elements in the input {\tt SparseMatrix} 
! argument {\tt sMat}.  
!
! !INTERFACE:

    integer function lsize_(sMat)
!
! !USES:
!
      use m_AttrVect,only : AttrVect_lsize => lsize

      implicit none

      type(SparseMatrix), intent(in) :: sMat

! !REVISION HISTORY:
!       23Apr00 - J.W. Larson <larson@mcs.anl.gov> - initial version.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::lsize_'

  lsize_ = AttrVect_lsize(sMat%data)

 end function lsize_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE:  GlobalNumElements_ - Number of nonzero elements.
!
! !DESCRIPTION:  This routine computes the number of nonzero elements 
! in a distributed {\tt SparseMatrix} variable {\tt sMat}.  The input 
! {\tt SparseMatrix} argument {\tt sMat} is examined on each process 
! to determine the number of nonzero elements it holds, and this value 
! is summed across the communicator associated with the input 
! {\tt INTEGER} handle {\tt comm}, with the total returned {\em on each
! process on the communicator}.
!
! !INTERFACE:

 integer function GlobalNumElements_(sMat, comm)

!
! !USES:
!
      use m_die
      use m_mpif90

      implicit none

! !INPUT PARAMETERS: 

      type(SparseMatrix), intent(in)  :: sMat
      integer, optional,  intent(in)  :: comm

! !REVISION HISTORY:
!       24Apr01 - Jay Larson <larson@mcs.anl.gov> - New routine.
!
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//':: GlobalNumElements_'

  integer :: MyNumElements, GNumElements, ierr

       ! Determine the number of locally held nonzero elements:

  MyNumElements = lsize_(sMat)

  call MPI_ALLREDUCE(MyNumElements, GNumElements, 1, MP_INTEGER, &
                     MP_SUM, comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,"MPI_ALLREDUCE(MyNumElements...",ierr)
  endif

 GlobalNumElements_ = GNumElements

 end function GlobalNumElements_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexIA_ - Index integer attribute of a SparseMatrix.
!
! !DESCRIPTION:  This {\tt INTEGER} function reports the row index 
! for a given {\tt INTEGER} attribute of the input {\tt SparseMatrix} 
! argument {\tt sMat}.  The attribute requested is represented by the 
! input {\tt CHARACTER} variable {\tt attribute}.  The list of integer 
! attributes one can request is defined in the description block of the 
! header of this module ({\tt m\_SparseMatrix}).
!
! Here is how {\tt indexIA\_} provides access to integer attribute data
! in a {\tt SparseMatrix} variable {\tt sMat}.  Suppose we wish to access
! global row information.  This attribute has associated with it the 
! string tag {\tt grow}.  The corresponding index returned ({\tt igrow}) 
! is determined by invoking {\tt indexIA\_}:
! \begin{verbatim}
! igrow = indexIA_(sMat, 'grow')
! \end{verbatim}
!
! Access to the global row index data in {\tt sMat} is thus obtained by 
! referencing {\tt sMat\%data\%iAttr(igrow,:)}.
!
!
! !INTERFACE:

    integer function indexIA_(sMat, attribute)
!
! !USES:
!
      use m_AttrVect,only : AttrVect_indexIA => indexIA

      implicit none

      type(SparseMatrix), intent(in) :: sMat
      character(len=*),   intent(in) :: attribute

! !REVISION HISTORY:
!       23Apr00 - J.W. Larson <larson@mcs.anl.gov> - initial version.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexIA_'

  indexIA_ = AttrVect_indexIA(sMat%data, attribute)

 end function indexIA_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: indexRA_ - Index real attribute of a SparseMatrix.
!
! !DESCRIPTION:  This {\tt INTEGER} function reports the row index 
! for a given {\tt REAL} attribute of the input {\tt SparseMatrix} 
! argument {\tt sMat}.  The attribute requested is represented by the 
! input {\tt CHARACTER} variable {\tt attribute}.  The list of real 
! attributes one can request is defined in the description block of the 
! header of this module ({\tt m\_SparseMatrix}).
!
! Here is how {\tt indexRA\_} provides access to integer attribute data
! in a {\tt SparseMatrix} variable {\tt sMat}.  Suppose we wish to access
! matrix element values.  This attribute has associated with it the 
! string tag {\tt weight}.  The corresponding index returned ({\tt iweight}) 
! is determined by invoking {\tt indexRA\_}:
! \begin{verbatim}
! iweight = indexRA_(sMat, 'weight')
! \end{verbatim}
!
! Access to the matrix element data in {\tt sMat} is thus obtained by 
! referencing {\tt sMat\%data\%rAttr(iweight,:)}.
!
! !INTERFACE:

    integer function indexRA_(sMat, attribute)
!
! !USES:
!
      use m_AttrVect,only : AttrVect_indexRA => indexRA

      implicit none

      type(SparseMatrix), intent(in) :: sMat
      character(len=*),   intent(in) :: attribute

! !REVISION HISTORY:
!       24Apr00 - J.W. Larson <larson@mcs.anl.gov> - initial version.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::indexRA_'

  indexRA_ = AttrVect_indexRA(sMat%data, attribute)

 end function indexRA_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nRows_ - Return the number of rows in a SparseMatrix.
!
! !DESCRIPTION:  This routine returns the {\em total} number of rows
! in the input {\tt SparseMatrix} argument {\tt sMat}.  This number of
! rows is a constant, and not dependent on the decomposition of the 
! {\tt SparseMatrix}.
!
! !INTERFACE:

    integer function nRows_(sMat)
!
! !USES:
!
      implicit none

      type(SparseMatrix), intent(in) :: sMat

! !REVISION HISTORY:
!       19Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nRows_'

  nRows_ = sMat%nrows

 end function nRows_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: nCols_ - Return the number of columns in a SparseMatrix.
!
! !DESCRIPTION:  This routine returns the {\em total} number of columns
! in the input {\tt SparseMatrix} argument {\tt sMat}.  This number of
! columns is a constant, and not dependent on the decomposition of the 
! {\tt SparseMatrix}.
!
! !INTERFACE:

    integer function nCols_(sMat)
!
! !USES:
!
      implicit none

      type(SparseMatrix), intent(in) :: sMat

! !REVISION HISTORY:
!       19Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::nCols_'

  nCols_ = sMat%ncols

 end function nCols_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportGlobalRowIndices_ - return Global Row Indices
!
! !DESCRIPTION:
! This routine extracts from the input {\tt SparseMatrix} argument 
! {\tt sMat} its global row indices, and returns them in the {\tt INTEGER} 
! output array {\tt GlobalRows}, and its length in the output {\tt INTEGER} 
! argument {\tt length}.
!
! {\bf N.B.:}  The flexibility of this routine regarding the pointer 
! association status of the output argument {\tt GlobalRows} means the
! user must invoke this routine with care.  If the user wishes this
! routine to fill a pre-allocated array, then obviously this array
! must be allocated prior to calling this routine.  If the user wishes
! that the routine {\em create} the output argument array {\tt GlobalRows},
! then the user must ensure this pointer is not allocated (i.e. the user
! must nullify this pointer) at the time this routine is invoked.
!
! {\bf N.B.:}  If the user has relied on this routine to allocate memory
! associated with the pointer {\tt GlobalRows}, then the user is responsible 
! for deallocating this array once it is no longer needed.  Failure to 
! do so will result in a memory leak.
!
! !INTERFACE:

 subroutine exportGlobalRowIndices_(sMat, GlobalRows, length)
!
! !USES:
!
      use m_die 
      use m_stdio

      use m_AttrVect,      only : AttrVect_exportIAttr => exportIAttr

      implicit none

! !INPUT PARAMETERS: 

      type(SparseMatrix),     intent(in)  :: sMat

! !OUTPUT PARAMETERS: 

      integer,  dimension(:), pointer     :: GlobalRows
      integer,                intent(out) :: length

! !REVISION HISTORY:
!        7May02 - J.W. Larson <larson@mcs.anl.gov> - initial version.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::exportGlobalRowIndices_'

       ! Export the data (inheritance from AttrVect)

  call AttrVect_exportIAttr(sMat%data, 'grow', GlobalRows, length)

 end subroutine exportGlobalRowIndices_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportGlobalColumnIndices_ - return Global Column Indices
!
! !DESCRIPTION:
! This routine extracts from the input {\tt SparseMatrix} argument 
! {\tt sMat} its global column indices, and returns them in the {\tt INTEGER} 
! output array {\tt GlobalColumns}, and its length in the output {\tt INTEGER} 
! argument {\tt length}.
!
! {\bf N.B.:}  The flexibility of this routine regarding the pointer 
! association status of the output argument {\tt GlobalColumns} means the
! user must invoke this routine with care.  If the user wishes this
! routine to fill a pre-allocated array, then obviously this array
! must be allocated prior to calling this routine.  If the user wishes
! that the routine {\em create} the output argument array {\tt GlobalColumns},
! then the user must ensure this pointer is not allocated (i.e. the user
! must nullify this pointer) at the time this routine is invoked.
!
! {\bf N.B.:}  If the user has relied on this routine to allocate memory
! associated with the pointer {\tt GlobalColumns}, then the user is responsible 
! for deallocating this array once it is no longer needed.  Failure to 
! do so will result in a memory leak.
!
! !INTERFACE:

 subroutine exportGlobalColumnIndices_(sMat, GlobalColumns, length)

!
! !USES:
!
      use m_die 
      use m_stdio

      use m_AttrVect,      only : AttrVect_exportIAttr => exportIAttr

      implicit none

! !INPUT PARAMETERS: 

      type(SparseMatrix),     intent(in)  :: sMat

! !OUTPUT PARAMETERS: 

      integer,  dimension(:), pointer     :: GlobalColumns
      integer,                intent(out) :: length

! !REVISION HISTORY:
!        7May02 - J.W. Larson <larson@mcs.anl.gov> - initial version.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::exportGlobalColumnIndices_'

       ! Export the data (inheritance from AttrVect)

  call AttrVect_exportIAttr(sMat%data, 'gcol', GlobalColumns, length)

 end subroutine exportGlobalColumnIndices_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportLocalRowIndices_ - return Local Row Indices
!
! !DESCRIPTION:
! This routine extracts from the input {\tt SparseMatrix} argument 
! {\tt sMat} its local row indices, and returns them in the {\tt INTEGER} 
! output array {\tt LocalRows}, and its length in the output {\tt INTEGER} 
! argument {\tt length}.
!
! {\bf N.B.:}  The flexibility of this routine regarding the pointer 
! association status of the output argument {\tt LocalRows} means the
! user must invoke this routine with care.  If the user wishes this
! routine to fill a pre-allocated array, then obviously this array
! must be allocated prior to calling this routine.  If the user wishes
! that the routine {\em create} the output argument array {\tt LocalRows},
! then the user must ensure this pointer is not allocated (i.e. the user
! must nullify this pointer) at the time this routine is invoked.
!
! {\bf N.B.:}  If the user has relied on this routine to allocate memory
! associated with the pointer {\tt LocalRows}, then the user is responsible 
! for deallocating this array once it is no longer needed.  Failure to 
! do so will result in a memory leak.
!
! !INTERFACE:

 subroutine exportLocalRowIndices_(sMat, LocalRows, length)
!
! !USES:
!
      use m_die 
      use m_stdio

      use m_AttrVect,      only : AttrVect_exportIAttr => exportIAttr

      implicit none

! !INPUT PARAMETERS: 

      type(SparseMatrix),     intent(in)  :: sMat

! !OUTPUT PARAMETERS: 

      integer,  dimension(:), pointer     :: LocalRows
      integer,                intent(out) :: length

! !REVISION HISTORY:
!        7May02 - J.W. Larson <larson@mcs.anl.gov> - initial version.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::exportLocalRowIndices_'

       ! Export the data (inheritance from AttrVect)

  call AttrVect_exportIAttr(sMat%data, 'lrow', LocalRows, length)

 end subroutine exportLocalRowIndices_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportLocalColumnIndices_ - return Local Column Indices
!
! !DESCRIPTION:
! This routine extracts from the input {\tt SparseMatrix} argument 
! {\tt sMat} its local column indices, and returns them in the {\tt INTEGER} 
! output array {\tt LocalColumns}, and its length in the output {\tt INTEGER} 
! argument {\tt length}.
!
! {\bf N.B.:}  The flexibility of this routine regarding the pointer 
! association status of the output argument {\tt LocalColumns} means the
! user must invoke this routine with care.  If the user wishes this
! routine to fill a pre-allocated array, then obviously this array
! must be allocated prior to calling this routine.  If the user wishes
! that the routine {\em create} the output argument array {\tt LocalColumns},
! then the user must ensure this pointer is not allocated (i.e. the user
! must nullify this pointer) at the time this routine is invoked.
!
! {\bf N.B.:}  If the user has relied on this routine to allocate memory
! associated with the pointer {\tt LocalColumns}, then the user is responsible 
! for deallocating this array once it is no longer needed.  Failure to 
! do so will result in a memory leak.
!
! !INTERFACE:

 subroutine exportLocalColumnIndices_(sMat, LocalColumns, length)

!
! !USES:
!
      use m_die 
      use m_stdio

      use m_AttrVect,      only : AttrVect_exportIAttr => exportIAttr

      implicit none

! !INPUT PARAMETERS: 

      type(SparseMatrix),     intent(in)  :: sMat

! !OUTPUT PARAMETERS: 

      integer,  dimension(:), pointer     :: LocalColumns
      integer,                intent(out) :: length

! !REVISION HISTORY:
!        7May02 - J.W. Larson <larson@mcs.anl.gov> - initial version.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::exportLocalColumnIndices_'

       ! Export the data (inheritance from AttrVect)

  call AttrVect_exportIAttr(sMat%data, 'lcol', LocalColumns, length)

 end subroutine exportLocalColumnIndices_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: exportMatrixElements_ - return MatrixElements as an array.
!
! !DESCRIPTION:
! This routine extracts the matrix elements from the input {\tt SparseMatrix} 
! argument {\tt sMat}, and returns them in the {\tt REAL} output array 
! {\tt MatrixElements}, and its length in the output {\tt INTEGER} 
! argument {\tt length}.
!
! {\bf N.B.:}  The flexibility of this routine regarding the pointer 
! association status of the output argument {\tt MatrixElements} means the
! user must invoke this routine with care.  If the user wishes this
! routine to fill a pre-allocated array, then obviously this array
! must be allocated prior to calling this routine.  If the user wishes
! that the routine {\em create} the output argument array {\tt MatrixElements},
! then the user must ensure this pointer is not allocated (i.e. the user
! must nullify this pointer) at the time this routine is invoked.
!
! {\bf N.B.:}  If the user has relied on this routine to allocate memory
! associated with the pointer {\tt MatrixElements}, then the user is responsible 
! for deallocating this array once it is no longer needed.  Failure to 
! do so will result in a memory leak.
!
! !INTERFACE:

 subroutine exportMatrixelements_(sMat, MatrixElements, length)

!
! !USES:
!
      use m_die 
      use m_stdio

      use m_AttrVect,      only : AttrVect_exportRAttr => exportRAttr

      implicit none

! !INPUT PARAMETERS: 

      type(SparseMatrix),     intent(in)  :: sMat

! !OUTPUT PARAMETERS: 

      real,  dimension(:),    pointer     :: MatrixElements
      integer,                intent(out) :: length

! !REVISION HISTORY:
!        7May02 - J.W. Larson <larson@mcs.anl.gov> - initial version.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::exportMatrixElements_'

       ! Export the data (inheritance from AttrVect)

  call AttrVect_exportRAttr(sMat%data, 'weight', MatrixElements, length)

 end subroutine exportMatrixElements_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: importGlobalRowIndices_ - Import Global Row Indices
!
! !DESCRIPTION:
! This routine imports global row index data into the {\tt SparseMatrix} 
! argument {\tt sMat}.  The user provides the index data in the input 
! {\tt INTEGER} vector {\tt inVect}.  The input {\tt INTEGER} argument 
! {\tt lsize} is used as a consistencey check to ensure the user is
! sufficient space in the {\tt SparseMatrix} to store the data.
!
! !INTERFACE:

 subroutine importGlobalRowIndices_(sMat, inVect, lsize)

!
! !USES:
!
      use m_die
      use m_stdio

      use m_AttrVect,      only : AttrVect_importIAttr => importIAttr

      implicit none

! !INPUT PARAMETERS: 

      integer,  dimension(:), pointer       :: inVect
      integer,                intent(in)    :: lsize

! !INPUT/OUTPUT PARAMETERS: 

      type(SparseMatrix),     intent(inout) :: sMat

! !REVISION HISTORY:
!        7May02 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::importGlobalRowIndices_'

       ! Argument Check:

  if(lsize > lsize_(sMat)) then
     write(stderr,*) myname_,':: ERROR, lsize > lsize_(sMat).', &
          'lsize = ',lsize,'lsize_(sMat) = ',lsize_(sMat)
     call die(myname_)
  endif

       ! Import the data (inheritance from AttrVect)

  call AttrVect_importIAttr(sMat%data, 'grow', inVect, lsize)

 end subroutine importGlobalRowIndices_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: importGlobalColumnIndices_ - Import Global Column Indices
!
! !DESCRIPTION:
! This routine imports global column index data into the {\tt SparseMatrix} 
! argument {\tt sMat}.  The user provides the index data in the input 
! {\tt INTEGER} vector {\tt inVect}.  The input {\tt INTEGER} argument 
! {\tt lsize} is used as a consistencey check to ensure the user is
! sufficient space in the {\tt SparseMatrix} to store the data.
!
! !INTERFACE:

 subroutine importGlobalColumnIndices_(sMat, inVect, lsize)

!
! !USES:
!
      use m_die
      use m_stdio

      use m_AttrVect,      only : AttrVect_importIAttr => importIAttr

      implicit none

! !INPUT PARAMETERS: 

      integer,  dimension(:), pointer       :: inVect
      integer,                intent(in)    :: lsize

! !INPUT/OUTPUT PARAMETERS: 

      type(SparseMatrix),     intent(inout) :: sMat

! !REVISION HISTORY:
!        7May02 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::importGlobalColumnIndices_'

       ! Argument Check:

  if(lsize > lsize_(sMat)) then
     write(stderr,*) myname_,':: ERROR, lsize > lsize_(sMat).', &
          'lsize = ',lsize,'lsize_(sMat) = ',lsize_(sMat)
     call die(myname_)
  endif

       ! Import the data (inheritance from AttrVect)

  call AttrVect_importIAttr(sMat%data, 'gcol', inVect, lsize)

 end subroutine importGlobalColumnIndices_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: importLocalRowIndices_ - Import Local Row Indices
!
! !DESCRIPTION:
! This routine imports local row index data into the {\tt SparseMatrix} 
! argument {\tt sMat}.  The user provides the index data in the input 
! {\tt INTEGER} vector {\tt inVect}.  The input {\tt INTEGER} argument 
! {\tt lsize} is used as a consistencey check to ensure the user is
! sufficient space in the {\tt SparseMatrix} to store the data.
!
! !INTERFACE:

 subroutine importLocalRowIndices_(sMat, inVect, lsize)

!
! !USES:
!
      use m_die
      use m_stdio

      use m_AttrVect,      only : AttrVect_importIAttr => importIAttr

      implicit none

! !INPUT PARAMETERS: 

      integer,  dimension(:), pointer       :: inVect
      integer,                intent(in)    :: lsize

! !INPUT/OUTPUT PARAMETERS: 

      type(SparseMatrix),     intent(inout) :: sMat

! !REVISION HISTORY:
!        7May02 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::importLocalRowIndices_'

       ! Argument Check:

  if(lsize > lsize_(sMat)) then
     write(stderr,*) myname_,':: ERROR, lsize > lsize_(sMat).', &
          'lsize = ',lsize,'lsize_(sMat) = ',lsize_(sMat)
     call die(myname_)
  endif

       ! Import the data (inheritance from AttrVect)

  call AttrVect_importIAttr(sMat%data, 'grow', inVect, lsize)

 end subroutine importLocalRowIndices_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: importLocalColumnIndices_ - Import Local Column Indices
!
! !DESCRIPTION:
! This routine imports local column index data into the {\tt SparseMatrix} 
! argument {\tt sMat}.  The user provides the index data in the input 
! {\tt INTEGER} vector {\tt inVect}.  The input {\tt INTEGER} argument 
! {\tt lsize} is used as a consistencey check to ensure the user is
! sufficient space in the {\tt SparseMatrix} to store the data.
!
! !INTERFACE:

 subroutine importLocalColumnIndices_(sMat, inVect, lsize)

!
! !USES:
!
      use m_die
      use m_stdio

      use m_AttrVect,      only : AttrVect_importIAttr => importIAttr

      implicit none

! !INPUT PARAMETERS: 

      integer,  dimension(:), pointer       :: inVect
      integer,                intent(in)    :: lsize

! !INPUT/OUTPUT PARAMETERS: 

      type(SparseMatrix),     intent(inout) :: sMat

! !REVISION HISTORY:
!        7May02 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::importLocalColumnIndices_'

       ! Argument Check:

  if(lsize > lsize_(sMat)) then
     write(stderr,*) myname_,':: ERROR, lsize > lsize_(sMat).', &
          'lsize = ',lsize,'lsize_(sMat) = ',lsize_(sMat)
     call die(myname_)
  endif

       ! Import the data (inheritance from AttrVect)

  call AttrVect_importIAttr(sMat%data, 'gcol', inVect, lsize)

 end subroutine importLocalColumnIndices_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: importMatrixElements_ - Import Local Column Indices
!
! !DESCRIPTION:
! This routine imports matrix elements index data into the 
! {\tt SparseMatrix} argument {\tt sMat}.  The user provides the index 
! data in the input {\tt REAL} vector {\tt inVect}.  The input 
! {\tt INTEGER} argument {\tt lsize} is used as a consistencey check 
! to ensure the user is sufficient space in the {\tt SparseMatrix} 
! to store the data.
!
! !INTERFACE:

 subroutine importMatrixElements_(sMat, inVect, lsize)

!
! !USES:
!
      use m_die
      use m_stdio

      use m_AttrVect,      only : AttrVect_importRAttr => importRAttr

      implicit none

! !INPUT PARAMETERS: 

      real,  dimension(:),    pointer       :: inVect
      integer,                intent(in)    :: lsize

! !INPUT/OUTPUT PARAMETERS: 

      type(SparseMatrix),     intent(inout) :: sMat

! !REVISION HISTORY:
!        7May02 - J.W. Larson <larson@mcs.anl.gov> - initial prototype.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::importMatrixElements_'

       ! Argument Check:

  if(lsize > lsize_(sMat)) then
     write(stderr,*) myname_,':: ERROR, lsize > lsize_(sMat).', &
          'lsize = ',lsize,'lsize_(sMat) = ',lsize_(sMat)
     call die(myname_)
  endif

       ! Import the data (inheritance from AttrVect)

  call AttrVect_importRAttr(sMat%data, 'weight', inVect, lsize)

 end subroutine importMatrixElements_

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
!       23Apr01 - Jay Larson <larson@mcs.anl.gov> - Modified to accomodate
!                 changes to the SparseMatrix type.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::local_row_range_'

  integer :: i, ilrow, lsize

  ilrow = AttrVect_indexIA(sMat%data, 'lrow')
  lsize = AttrVect_lsize(sMat%data)

       ! Initialize start_row and end_row:

  start_row = sMat%data%iAttr(ilrow,1)
  end_row = sMat%data%iAttr(ilrow,1)

  do i=1,lsize
     start_row = min(start_row, sMat%data%iAttr(ilrow,i))
     end_row = max(end_row, sMat%data%iAttr(ilrow,i))
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
!       23Apr01 - Jay Larson <larson@mcs.anl.gov> - Modified to accomodate
!                 changes to the SparseMatrix type.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::global_row_range_'

  integer :: i, igrow, lsize

  igrow = AttrVect_indexIA(sMat%data, 'grow')
  lsize = AttrVect_lsize(sMat%data)

       ! Initialize start_row and end_row:

  start_row = sMat%data%iAttr(igrow,1)
  end_row = sMat%data%iAttr(igrow,1)

  do i=1,lsize
     start_row = min(start_row, sMat%data%iAttr(igrow,i))
     end_row = max(end_row, sMat%data%iAttr(igrow,i))
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
!       23Apr01 - Jay Larson <larson@mcs.anl.gov> - Modified to accomodate
!                 changes to the SparseMatrix type.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::local_col_range_'

  integer :: i, ilcol, lsize

  ilcol = AttrVect_indexIA(sMat%data, 'lcol')
  lsize = AttrVect_lsize(sMat%data)

       ! Initialize start_col and end_col:

  start_col = sMat%data%iAttr(ilcol,1)
  end_col = sMat%data%iAttr(ilcol,1)

  do i=1,lsize
     start_col = min(start_col, sMat%data%iAttr(ilcol,i))
     end_col = max(end_col, sMat%data%iAttr(ilcol,i))
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
!       23Apr01 - Jay Larson <larson@mcs.anl.gov> - Modified to accomodate
!                 changes to the SparseMatrix type.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::global_col_range_'

  integer :: i, igcol, lsize

  igcol = AttrVect_indexIA(sMat%data, 'lcol')
  lsize = AttrVect_lsize(sMat%data)

       ! Initialize start_col and end_col:

  start_col = sMat%data%iAttr(igcol,1)
  end_col = sMat%data%iAttr(igcol,1)

  do i=1,lsize
     start_col = min(start_col, sMat%data%iAttr(igcol,i))
     end_col = max(end_col, sMat%data%iAttr(igcol,i))
  end do

 end subroutine global_col_range_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: ComputeSparsity_ - sparsity of a SparseMatrix
!
! !DESCRIPTION:  This routine computes the sparsity of a consolidated
! (all on one process) or distributed {\tt SparseMatrix}.  The input 
! {\tt SparseMatrix} argument {\tt sMat} is examined to determine the
! number of nonzero elements it holds, and this value is divided by the
! product of the number of rows and columns in {\tt sMat}.  If the 
! optional input argument {\tt comm} is given, then the distributed 
! elements are counted and the sparsity computed accordingly, and the 
! resulting value of {\tt sparsity} is returned {\em to all processes}.
!
! Given the inherent problems with multiplying and dividing large integers,
! the work in this routine is performed using floating point arithmetic on
! the logarithms of the number of rows, columns, and nonzero elements.
!
! !INTERFACE:

 subroutine ComputeSparsity_(sMat, sparsity, comm)

!
! !USES:
!
      use m_die
      use m_mpif90

      use m_AttrVect, only : AttrVect_lsize => lsize

      implicit none

! !INPUT PARAMETERS: 

      type(SparseMatrix), intent(in)  :: sMat
      integer, optional,  intent(in)  :: comm

! !OUTPUT PARAMETERS:

      real,               intent(out) :: sparsity

! !REVISION HISTORY:
!       23Apr01 - Jay Larson <larson@mcs.anl.gov> - New routine.
!
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::ComputeSparsity_'

  integer :: num_elements, num_rows, num_cols
  real    :: Lnum_elements, Lnum_rows, Lnum_cols, LMySparsity
  real    :: MySparsity
  integer :: ierr

       ! Extract number of nonzero elements and compute its logarithm

  num_elements = lsize_(sMat)
  Lnum_elements = log(float(num_elements))

       ! Extract number of rows and compute its logarithm

  num_rows = nRows_(sMat)
  Lnum_rows = log(float(num_rows))

       ! Extract number of columns and compute its logarithm

  num_cols = nCols_(sMat)
  Lnum_cols = log(float(num_cols))  

       ! Compute logarithm of the (local) sparsity

  LMySparsity = Lnum_elements - Lnum_rows - Lnum_cols

       ! Compute the (local) sparsity from its logarithm.

  MySparsity = exp(LMySparsity)

       ! If a communicator handle is present, sum up the
       ! distributed sparsity values to all processes.  If not,
       ! return the value of MySparsity computed above.

  if(present(comm)) then
     call MPI_ALLREDUCE(MySparsity, sparsity, 1, MP_INTEGER, &
                        MP_SUM, comm, ierr)
     if(ierr /= 0) then
	call MP_perr_die(myname_,"MPI_ALLREDUCE(MySparsity...",ierr)
     endif
  else
     sparsity = MySparsity
  endif

 end subroutine ComputeSparsity_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: CheckBounds_ - Check for out-of-bounds row/column values.
!
! !DESCRIPTION:  This routine examines the input distributed 
! {\tt SparseMatrix} variable {\tt sMat}, and examines the global row
! and column index for each element, comparing them with the known 
! maximum values for each (as returned by the routines {\tt nRows\_()}
! and {\tt nCols\_()}, respectively).  If global row or column entries 
! are non-positive, or greater than the defined maximum values, this
! routine stops execution with an error message.  If no out-of-bounds
! values are detected, the output {\tt INTEGER} status {\tt ierror} is 
! set to zero.
!
! !INTERFACE:

 subroutine CheckBounds_(sMat, ierror)
!
! !USES:
!
      use m_die

      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_indexIA => indexIA

      implicit none

      type(SparseMatrix), intent(in)  :: sMat
      integer,            intent(out) :: ierror

! !REVISION HISTORY:
!       24Apr01 - Jay Larson <larson@mcs.anl.gov> - Initial prototype.
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::CheckBounds_'

  integer :: MaxRow, MaxCol, NumElements
  integer :: igrow, igcol
  integer :: i

       ! Initially, set ierror to zero (success):

  ierror = 0

       ! Query sMat to find the number of rows and columns:

  MaxRow = nRows_(sMat)
  MaxCol = nCols_(sMat)

       ! Query sMat for the number of nonzero elements:

  NumElements = lsize_(sMat)

       ! Query sMat to index global row and column storage indices:

  igrow = indexIA_(sMat,'grow')
  igcol = indexIA_(sMat,'gcol')

       ! Scan the entries of sMat for row or column elements that
       ! are out-of-bounds.  Here, out-of-bounds means:  1) non-
       ! positive row or column indices; 2) row or column indices
       ! exceeding the stated number of rows or columns.

  do i=1,NumElements

       ! Row index out of bounds?

     if((sMat%data%iAttr(igrow,i) > MaxRow) .or. &
	  (sMat%data%iAttr(igrow,i) <= 0)) then
	ierror = 1
	call die(myname_,"Row index out of bounds",ierror)
     endif

       ! Column index out of bounds?

     if((sMat%data%iAttr(igcol,i) > MaxCol) .or. &
	  (sMat%data%iAttr(igcol,i) <= 0)) then
	ierror = 2
	call die(myname_,"Column index out of bounds",ierror)
     endif

  end do

 end subroutine CheckBounds_

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
!       23Apr01 - Jay Larson <larson@mcs.anl.gov> - Modified to accomodate
!                 changes to the SparseMatrix type.
!       18May01 - R. Jacob <jacob@mcs.anl.gov> - Use MP_TYPE function
!                 to set type in the mpi_allreduce
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::row_sum_'

  integer :: i, igrow, ierr, iwgt, lsize, myID
  integer :: start_row, end_row
  integer :: mp_Type_lsums
  real, dimension(:), allocatable :: lsums

       ! Determine local rank

  call MP_COMM_RANK(comm, myID, ierr)

       ! Determine on each process the row of global row indices:

  call global_row_range_(sMat, comm, start_row, end_row)

       ! Determine across the communicator the _maximum_ value of
       ! end_row, which will be assigned to num_rows on each process:

  call MPI_ALLREDUCE(end_row, num_rows, 1, MP_INTEGER, MP_MAX, &
                    comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,"MPI_ALLREDUCE(end_row...",ierr)
  endif

       ! Allocate storage for the sums on each process.

  allocate(lsums(num_rows), sums(num_rows), stat=ierr)

  if(ierr /= 0) then
     call die(myname_,"allocate(lsums(...",ierr)
  endif

       ! Compute the local entries to lsum(1:num_rows) for each process:

  lsize = AttrVect_lsize(sMat%data)
  igrow = AttrVect_indexIA(sMat%data,'grow')
  iwgt = AttrVect_indexRA(sMat%data,'weight')

  lsums = 0.
  do i=1,lsize
     lsums(sMat%data%iAttr(igrow,i)) = lsums(sMat%data%iAttr(igrow,i)) + &
	                           sMat%data%rAttr(iwgt,i)
  end do

       ! Compute the global sum of the entries of lsums so that all
       ! processes own the global sums.

  mp_Type_lsums=MP_Type(lsums)
  call MPI_ALLREDUCE(lsums, sums, num_rows, mp_Type_lsums, MP_SUM, comm, ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,"MPI_ALLREDUCE(lsums...",ierr)
  endif

       ! Clean up...

  deallocate(lsums, stat=ierr)
  if(ierr /= 0) then
     call die(myname_,"deallocate(lsums...",ierr)
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

!BOP -------------------------------------------------------------------
!
! !IROUTINE: Sort_ - return index permutation keyed by a list of
!            attributes
!
! !DESCRIPTION:
! The subroutine {\tt Sort\_()} uses a list of sorting keys defined by 
! the input {\tt List} argument {\tt key\_list}, searches for the appropriate 
! integer or real attributes referenced by the items in {\tt key\_list} 
! ( that is, it identifies the appropriate entries in {sMat\%data\%iList} 
! and {\tt sMat\%data\%rList}), and then uses these keys to generate an index 
! permutation {\tt perm} that will put the nonzero matrix entries of stored
! in {\tt sMat\%data} in lexicographic order as defined by {\tt key\_ist} 
! (the ordering in {\tt key\_list} being from left to right.  The optional 
! {\tt LOGICAL} array input argument {\tt descend} specifies whether or
! not to sort by each key in {\em descending} order or {\em ascending} 
! order.  Entries in {\tt descend} that have value {\tt .TRUE.} correspond 
! to a sort by the corresponding key in descending order.  If the argument
! {\tt descend} is not present, the sort is performed for all keys in
! ascending order.
!
! !INTERFACE:

 subroutine Sort_(sMat, key_list, perm, descend)

!
! !USES:
!
      use m_die ,          only : die
      use m_stdio ,        only : stderr

      use m_List ,         only : List

      use m_AttrVect, only: AttrVect_Sort => Sort

      implicit none
!
! !INPUT PARAMETERS: 

      type(SparseMatrix),              intent(in) :: sMat
      type(List),                      intent(in) :: key_list
      logical, dimension(:), optional, intent(in) :: descend
!
! !OUTPUT PARAMETERS: 

      integer, dimension(:), pointer              :: perm


! !REVISION HISTORY:
!       24Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::Sort_'

  if(present(descend)) then
     call AttrVect_Sort(sMat%data, key_list, perm, descend)
  else
     call AttrVect_Sort(sMat%data, key_list, perm)
  endif

 end Subroutine Sort_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: Permute_ - permute SparseMatrix entries
!
! !DESCRIPTION:
! The subroutine {\tt Permute\_()} uses an input index permutation 
! {\tt perm} to re-order the entries of the {\tt SparseMatrix} argument 
! {\tt sMat}.  The index permutation {\tt perm} is generated using the 
! routine {\tt Sort\_()} (in this module).
!
! !INTERFACE:

 subroutine Permute_(sMat, perm)

!
! !USES:
!
      use m_die ,          only : die
      use m_stdio ,        only : stderr

      use m_AttrVect, only: AttrVect_Permute => Permute

      implicit none
!
! !INPUT PARAMETERS: 


      integer, dimension(:), pointer               :: perm
!
! !INPUT/OUTPUT PARAMETERS: 

      type(SparseMatrix),            intent(inout) :: sMat


! !REVISION HISTORY:
!       24Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::Permute_'

  call AttrVect_Permute(sMat%data, perm)

 end Subroutine Permute_

!BOP -------------------------------------------------------------------
!
! !IROUTINE: SortPermute_ - sort/permute SparseMatrix entries.
!
! !DESCRIPTION:
! The subroutine {\tt SortPermute\_()} uses a list of sorting keys defined 
! by the input {\tt List} argument {\tt key\_list}, searches for the 
! appropriate integer or real attributes referenced by the items in 
! {\tt key\_ist} ( that is, it identifies the appropriate entries in 
! {sMat\%data\%iList} and {\tt sMat\%data\%rList}), and then uses these 
! keys to generate an index permutation that will put the nonzero matrix 
! entries of stored in {\tt sMat\%data} in lexicographic order as defined 
! by {\tt key\_list} (the ordering in {\tt key\_list} being from left to 
! right.  The optional {\tt LOGICAL} array input argument {\tt descend} 
! specifies whether or not to sort by each key in {\em descending} order 
! or {\em ascending} order.  Entries in {\tt descend} that have value 
! {\tt .TRUE.} correspond to a sort by the corresponding key in descending 
! order.  If the argument {\tt descend} is not present, the sort is 
! performed for all keys in ascending order.
!
! Once this index permutation is created, it is applied to re-order the 
! entries of the {\tt SparseMatrix} argument {\tt sMat} accordingly.
!
! !INTERFACE:

 subroutine SortPermute_(sMat, key_list, descend)

!
! !USES:
!
      use m_die ,          only : die
      use m_stdio ,        only : stderr

      use m_List ,         only : List

      implicit none
!
! !INPUT PARAMETERS: 

      type(List),                      intent(in)    :: key_list
      logical, dimension(:), optional, intent(in)    :: descend
!
! !INPUT/OUTPUT PARAMETERS: 

      type(SparseMatrix),              intent(inout) :: sMat

! !REVISION HISTORY:
!       24Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::SortPermute_'

  integer, dimension(:), pointer :: perm

       ! Create index permutation perm(:)

  if(present(descend)) then
     call Sort_(sMat, key_list, perm, descend)
  else
     call Sort_(sMat, key_list, perm)
  endif

       ! Apply index permutation perm(:) to re-order sMat:

  call Permute_(sMat, perm)

 end subroutine SortPermute_

 end module m_SparseMatrix



