!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_SparseMatrixComms -- sparse matrix communications methods.
!
! !DESCRIPTION:
! The {\tt SparseMatrix} datatype provides sparse matrix storage for 
! the parallel matrix-vector multiplication ${\bf y} = {\bf M} {\bf x}$.
! This module provides communications services for the {\tt SparseMatrix}
! type.  These services include scattering matrix elements based on row or 
! column decompositions, gathering of matrix elements to the root, and 
! broadcasting from the root.
!
! !INTERFACE:

 module m_SparseMatrixComms

      private   ! except

! !PUBLIC MEMBER FUNCTIONS:
!
      public :: ScatterByColumn
      public :: ScatterByRow
      public :: Gather
      public :: Bcast

    interface ScatterByColumn ; module procedure &
         ScatterByColumnGSMap_
    end interface

    interface ScatterByRow ; module procedure &
         ScatterByRowGSMap_
    end interface

    interface Gather ; module procedure Gather_ ; end interface
    interface Bcast ; module procedure Bcast_ ; end interface

! !REVISION HISTORY:
!       13Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial prototype
!                 and API specifications.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_SparseMatrixComms'

 contains

!-------------------------------------------------------------------------
!     Math + Computer Science Division / Argonne National Laboratory     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  ScatterByColumnGSMap_()--Column-based scatter for SparseMatrix.
! 
! !INTERFACE:

 subroutine ScatterByColumnGSMap_(columnGSMap, GsMat, LsMat, root, comm, stat)
!
! !USES:
!

   use m_die, only : MP_perr_die

   use m_mpif90

   use m_GlobalSegMap, only : GlobalSegMap

   use m_SparseMatrix, only : SparseMatrix
   use m_SparseMatrix, only : SparseMatrix_nRows => nRows
   use m_SparseMatrix, only : SparseMatrix_nCols => nCols

   use m_SparseMatrixDecomp, only : SparseMatrixDecompByColumn => ByColumn

   use m_AttrVectComms, only : AttrVect_Scatter => scatter

   implicit none

! !INPUT PARAMETERS: 
!
   type(GlobalSegMap), intent(in)    :: columnGSMap
   type(SparseMatrix), intent(inout) :: GsMat
   integer,            intent(in)    :: root
   integer,            intent(in)    :: comm

! !OUTPUT PARAMETERS:
!
   type(SparseMatrix), intent(out) :: LsMat
   integer,            intent(out) :: stat

! !DESCRIPTION: This routine scatters the input consolidated {\tt SparseMatrix} 
! argument {\tt GsMat} (valid only on the root) to a distributed 
! {\tt SparseMatrix} variable {\tt LsMat} across all the processes present on
! the communicator associated with the integer handle {\tt comm}.  The 
! decomposition defining the scatter is supplied by the input {\tt GlobalSegMap} 
! argument {\tt columnGSMap}.  The output integer flag {\tt stat} signifies a 
! successful operation if it is returned with value zero.
!
! {\bf N.B.:}  This routine returns an allocated {\tt SparseMatrix} 
! variable {\tt LsMat}.  The user must destroy this variable when it
! is no longer needed by invoking {\tt SparseMatrix\_Clean()}.
!
! !REVISION HISTORY: 
!
!       13Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial API spec.

!EOP
!-------------------------------------------------------------------------

  character(len=*),parameter :: myname_=myname//'ScatterByColumnGSMap_'
! GlobalSegMap used to create column decomposition of GsMat
  type(GlobalSegMap) :: MatGSMap
! Storage for the number of rows and columns in the SparseMatrix
  integer :: NumRowsColumns(2)
! Process ID
  integer :: myID
! Error flag
  integer :: ierr

       ! Which process are we?

  call MPI_COMM_RANK(comm, myID, ierr)

  if(ierr /= 0) then
     call MP_perr_die(myname_,"MPI_COMM_RANK() failed",ierr)
  endif

       ! Create from columnGSMap the corresponding GlobalSegMap
       ! that will decompose GsMat by column the same way.

  call SparseMatrixDecompByColumn(columnGSMap, GsMat, MatGSMap)

       ! Scatter the matrix element data GsMat%data accordingly

  call AttrVect_Scatter(GsMat%data, LsMat%data, MatGSMap, root, comm, stat)

       ! Now, distribute to all the processes the number of Rows and
       ! columns in GsMat (which are valid on the root only at this point)

  if(myID == root) then
     NumRowsColumns(1) = SparseMatrix_nRows(GsMat)
     NumRowsColumns(2) = SparseMatrix_nCols(GsMat)
  endif

  call MPI_Bcast(NumRowsColumns, 2, MP_INTEGER, root, comm, ierr)

  if(ierr /= 0) then
     call MP_perr_die(myname_,"MPI_Bcast(NumRowsColumns...",ierr)
  endif

       ! Unpack NumRowsColumns

  LsMat%nrows = NumRowsColumns(1)
  LsMat%ncols = NumRowsColumns(2)

 end subroutine ScatterByColumnGSMap_

!-------------------------------------------------------------------------
!     Math + Computer Science Division / Argonne National Laboratory     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  ScatterByRowGSMap_()--Row-based scatter for SparseMatrix.
! 
! !INTERFACE:

 subroutine ScatterByRowGSMap_(rowGSMap, GsMat, LsMat, root, comm, stat)
!
! !USES:
!
   use m_die, only : MP_perr_die

   use m_mpif90

   use m_GlobalSegMap, only : GlobalSegMap

   use m_SparseMatrix, only : SparseMatrix
   use m_SparseMatrix, only : SparseMatrix_nRows => nRows
   use m_SparseMatrix, only : SparseMatrix_nCols => nCols

   use m_SparseMatrixDecomp, only : SparseMatrixDecompByColumn => ByColumn

   use m_AttrVectComms, only : AttrVect_Scatter => scatter

   implicit none

! !INPUT PARAMETERS: 
!
   type(GlobalSegMap), intent(in)    :: rowGSMap
   type(SparseMatrix), intent(inout) :: GsMat
   integer,            intent(in)    :: root
   integer,            intent(in)    :: comm

! !OUTPUT PARAMETERS:
!
   type(SparseMatrix), intent(out) :: LsMat
   integer,            intent(out) :: stat

! !DESCRIPTION: This routine scatters the input consolidated {\tt SparseMatrix} 
! argument {\tt GsMat} (valid only on the root) to a distributed 
! {\tt SparseMatrix} variable {\tt LsMat} across all the processes present on
! the communicator associated with the integer handle {\tt comm}.  The 
! decomposition defining the scatter is supplied by the input {\tt GlobalSegMap} 
! argument {\tt rowGSMap}.  The output integer flag {\tt stat} signifies a 
! successful operation if it is returned with value zero.
!
! {\bf N.B.:}  This routine returns an allocated {\tt SparseMatrix} 
! variable {\tt LsMat}.  The user must destroy this variable when it
! is no longer needed by invoking {\tt SparseMatrix\_Clean()}.
!
! !REVISION HISTORY: 
!
!       13Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial API spec.

!EOP
!-------------------------------------------------------------------------

  character(len=*),parameter :: myname_=myname//'ScatterByRowGSMap_'
! GlobalSegMap used to create row decomposition of GsMat
  type(GlobalSegMap) :: MatGSMap
! Storage for the number of rows and columns in the SparseMatrix
  integer :: NumRowsColumns(2)
! Process ID
  integer :: myID
! Error flag
  integer :: ierr

       ! Which process are we?

  call MPI_COMM_RANK(comm, myID, ierr)

  if(ierr /= 0) then
     call MP_perr_die(myname_,"MPI_COMM_RANK() failed",ierr)
  endif

       ! Create from rowGSMap the corresponding GlobalSegMap
       ! that will decompose GsMat by row the same way.

  call SparseMatrixDecompByRow(rowGSMap, GsMat, MatGSMap)

       ! Scatter the matrix element data GsMat%data accordingly

  call AttrVect_Scatter(GsMat%data, LsMat%data, MatGSMap, root, comm, stat)

       ! Now, distribute to all the processes the number of rows and
       ! columns in GsMat (which are valid on the root only at this point)

  if(myID == root) then
     NumRowsColumns(1) = SparseMatrix_nRows(GsMat)
     NumRowsColumns(2) = SparseMatrix_nCols(GsMat)
  endif

  call MPI_Bcast(NumRowsColumns, 2, MP_INTEGER, root, comm, ierr)

  if(ierr /= 0) then
     call MP_perr_die(myname_,"MPI_Bcast(NumRowsColumns...",ierr)
  endif

       ! Unpack NumRowsColumns

  LsMat%nrows = NumRowsColumns(1)
  LsMat%ncols = NumRowsColumns(2)

 end subroutine ScatterByRowGSMap_

!-------------------------------------------------------------------------
!     Math + Computer Science Division / Argonne National Laboratory     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Gather_()--Gather a distributed SparseMatrix to one process.
! 
! !INTERFACE:

 subroutine Gather_(LsMat, GsMat, GSMap, root, comm, stat)
!
! !USES:
!
   use m_GlobalSegMap, only: GlobalSegMap

   use m_SparseMatrix, only: SparseMatrix
   use m_SparseMatrix, only: SparseMatrix_nRows => nRows
   use m_SparseMatrix, only: SparseMatrix_nCols => nCols

   use m_AttrVectComms, only : AttrVect_Gather => gather

   implicit none

! !INPUT PARAMETERS: 
!
   type(SparseMatrix), intent(in) :: LsMat
   type(GlobalSegMap), intent(in) :: GSMap
   integer,            intent(in) :: root
   integer,            intent(in) :: comm

! !OUTPUT PARAMETERS:
!
   type(SparseMatrix), intent(out) :: GsMat
   integer,            intent(out) :: stat

! !DESCRIPTION: This routine gathers the input distributed 
! {\tt SparseMatrix} argument {\tt LsMat} to a consolidated {\tt SparseMatrix} 
! variable {\tt GsMat} on the root.  The decomposition defining the gather is 
! supplied by the input {\tt GlobalSegMap} argument {\tt GSMap}.  The status 
! flag {\tt stat} has value zero if the operation has succeeded.
!
! {\bf N.B.:}  This routine returns an allocated {\tt SparseMatrix} 
! variable {\tt GsMat}.  The user must destroy this variable when it
! is no longer needed by invoking {\tt SparseMatrix\_Clean()}.
!
! !REVISION HISTORY: 
!
!       13Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial API spec.

!EOP
!-------------------------------------------------------------------------

  character(len=*),parameter :: myname_=myname//'Gather_'

  call AttrVect_Gather(LsMat%data, GsMat%data, GSMap, root, comm, stat)

       ! For now, the GsMat inherits the number of rows and columns from
       ! the corresponding values of LsMat on the root (this should be
       ! checked in future versions).

  GsMat%nrows = SparseMatrix_nRows(LsMat)
  GsMat%ncols = SparseMatrix_nCols(LsMat)

 end subroutine Gather_

!-------------------------------------------------------------------------
!     Math + Computer Science Division / Argonne National Laboratory     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Bcast_()--Broadcast a SparseMatrix.
! 
! !INTERFACE:

 subroutine Bcast_(sMat, root, comm, stat)

!
! !USES:
!

   use m_die, only : MP_perr_die

   use m_mpif90

   use m_GlobalSegMap, only: GlobalSegMap

   use m_AttrVectComms, only : AttrVect_Bcast => bcast

   use m_SparseMatrix, only: SparseMatrix
   use m_SparseMatrix, only: SparseMatrix_nRows => nRows
   use m_SparseMatrix, only: SparseMatrix_nCols => nCols

   implicit none

! !INPUT PARAMETERS: 
!
   integer,            intent(in) :: root
   integer,            intent(in) :: comm

! !INPUT/OUTPUT PARAMETERS: 
!
   type(SparseMatrix), intent(inout) :: sMat

! !OUTPUT PARAMETERS:
!
   integer,            intent(out) :: stat

! !DESCRIPTION: This routine broadcasts the {\tt SparseMatrix} argument 
! {\tt sMat} from the root to all processes on the communicator associated
! with the communicator handle {\tt comm}.  The status flag {\tt stat} 
! has value zero if the operation has succeeded.
!
! {\bf N.B.:}  This routine returns an allocated {\tt SparseMatrix} 
! variable {\tt sMat}.  The user must destroy this variable when it
! is no longer needed by invoking {\tt SparseMatrix\_Clean()}.
!
! !REVISION HISTORY: 
!
!       13Apr01 - J.W. Larson <larson@mcs.anl.gov> - initial API spec.

!EOP
!-------------------------------------------------------------------------

  character(len=*),parameter :: myname_=myname//'Bcast_'

! Storage for the number of rows and columns in the SparseMatrix
  integer :: NumRowsColumns(2)
! Process ID number
  integer :: myID
! Error flag
  integer :: ierr

  call AttrVect_Bcast(sMat%data, root, comm, stat)

  if(myID == root) then
     NumRowsColumns(1) = SparseMatrix_nRows(sMat)
     NumRowsColumns(2) = SparseMatrix_nCols(sMat)
  endif

  call MPI_Bcast(NumRowsColumns, 2, MP_INTEGER, root, comm, ierr)

  if(ierr /= 0) then
     call MP_perr_die(myname_,"MPI_Bcast(NumRowsColumns...",ierr)
  endif

       ! Unpack NumRowsColumns on broadcast destination processes

  if(myID /= root) then
     sMat%nrows = NumRowsColumns(1)
     sMat%ncols = NumRowsColumns(2)
  endif

 end subroutine Bcast_

 end module m_SparseMatrixComms
