!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_SparseMatrixPlus -- Class Parallel for Matrix-Vector Multiplication
!
! !DESCRIPTION:
! Matrix-vector multiplication is one of the MCT's core services, and is 
! used primarily for the interpolation of data fields from one physical 
! grid to another.  Let ${\bf x} \in \Re^{N_x}$ and 
! ${\bf y} \in \Re^{N_y}$ represent data fields on physical grids $A$ 
! and $B$, respectively.  Field data is interpolated from grid $A$ to grid 
! $B$ by
! $$ {\bf y} = {\bf M} {\bf x} , $$
! where {\bf M} is aa  ${N_y} \times {N_x}$ matrix.
! 
! Within MCT, the {\tt SparseMatrix} data type is MCT's object for 
! storing sparse matrices such as {\bf M} , and the {\tt AttrVect} data 
! type is MCT's field data storage object.  That is, {\bf x} and {\bf y} 
! are each stored in {\tt AttrVect} form, and {\bf M} is stored as a 
! {\tt SparseMatrix}.  
!
! For global address spaces (uniprocessor or shared-memory parallel), this 
! picture of matrix-vector multiplication is sufficient.  If one wishes 
! to perform {\em distributed-memory parallel} matrix-vector multiplication, 
! however, in addition to computation, one must consider {\em communication}.
!
! There are three basic message-passing parallel strategies for computing 
! ${\bf y} = {\bf M} {\bf x}$:
!
!\begin{enumerate}
! \item Decompose {\bf M} based on its {\em rows}, and corresponding to the 
! decomposition for the vector {\bf y}.  That is, if a given process owns 
! the $i^{\rm th}$ element of {\bf y}, then all the elements of row $i$ of 
! {\bf M} also reside on this process.  Then  ${\bf y} = {\bf M} {\bf x}$ is 
! implemented as follows:
! \begin{enumerate}
! \item Create an {\em intermediate vector} {\bf x'} that is the pre-image of 
! the elements of {\bf y} owned locally.
! \item Comunnicate with the appropriate processes on the local communicator to 
! gather from {\bf x} the elements of {\bf x'}.
! \item Compute ${\bf y} = {\bf M} {\bf x'}$.
! \item Destroy the data structure holding {\bf x'}.
! \end{enumerate}
! \item Decompose {\bf M} based on its {\em columns}, and corresponding to the 
! decomposition for the vector {\bf x}.  That is, if a given process owns 
! the $j^{\rm th}$ element of {\bf x}, then all the elements of column $j$ of 
! {\bf M} also reside on this process.  Then  ${\bf y} = {\bf M} {\bf x}$ is 
! implemented as follows:
! \begin{enumerate}
! \item Create an {\em intermediate vector} {\bf y'} that holds {\em partial sums}
! of elements of {\bf y} computed from {\bf x} and {\bf M}.
! \item Compute ${\bf y'} = {\bf M} {\bf x}$.
! \item Perform communications to route elements of {\bf y'} to their eventual 
! destinations in {\bf y}, where they will be summed, resulting in the distributed 
! vector {\bf y}.
! \item Destroy the data structure holding {\bf y'}.
! \end{enumerate}
! \item Decompose {\bf M} based on some arbitrary, user-supplied scheme.  This will 
! necessitate two intermediate vectors {\bf x'} and {\bf y'}. Then  
! ${\bf y} = {\bf M} {\bf x}$ is implemented as follows: 
! \begin{enumerate}
! \item Create {\em intermediate vectors} {\bf x'} and {\bf y'}.  The numbers of 
! elements in {\bf x'} and {\bf y'} are based {\bf M}, specifically its numbers of 
! {\em distinct} row and column index values, respectively.
! \item Comunnicate with the appropriate processes on the local communicator to 
! gather from {\bf x} the elements of {\bf x'}.
! \item Compute ${\bf y'} = {\bf M} {\bf x'}$.
! \item Perform communications to route elements of {\bf y'} to their eventual 
! destinations in {\bf y}, where they will be summed, resulting in the distributed 
! vector {\bf y}.
! \item Destroy the data structures holding {\bf x'} and {\bf y'}.
! \end{enumerate}
! \end{enumerate}
!
! These operations require information about many aspects of the multiplication 
! process.  These data are:
! \begin{itemize}
! \item The matrix-vector parallelization strategy, which is one of the following:
! \begin{enumerate}
! \item Distributed in {\bf x}, purely data local in {\bf y}, labeled by the 
! public data member {\tt Xonly}
! \item Purely data local {\bf x}, distributed in {\bf y}, labeled by the 
! public data member {\tt Yonly}
! \item Distributed in both {\bf x} and {\bf y}, labeled by the public data 
! member {\tt XandY}
! \end{enumerate}
! \item A communications scheduler to create {\bf x'} from {\bf x};
! \item A communications scheduler to deliver partial sums contained in {\bf y'} to 
! {\bf y}.
! \item Lengths of the intermediate vectors {\bf x'} and {\bf y'}.
! \end{itemize}
!
! In MCT, the above data are stored in a {\em master} class for {\tt SparseMatrix}-
! {\tt AttrVect} multiplication.  This master class is called a 
! {\tt SparseMatrixPlus}.
!
! This module contains the definition of the {\tt SparseMatrixPlus}, and a variety 
! of methods to support it.  These include initialization, destruction, query, and 
! data import/export.
!
! !INTERFACE:

 module m_SparseMatrixPlus

! !USES:

      use m_String, only : String
      use m_SparseMatrix, only : SparseMatrix
      use m_Rearranger, only : Rearranger

! !PUBLIC TYPES:

      public :: SparseMatrixPlus

      Type SparseMatrixPlus
        logical :: Initialized
        type(String) :: Strategy
        integer :: XPrimeLength
        type(Rearranger) :: XToXPrime
        integer :: YPrimeLength
        type(Rearranger) :: YPrimeToY
        type(SparseMatrix) :: Matrix
      End Type SparseMatrixPlus

! !PUBLIC MEMBER FUNCTIONS:

      public :: init
      public :: clean
      public :: initialized

      interface init ; module procedure initFromRoot_ ; end interface
      interface clean ; module procedure clean_ ; end interface
      interface initialized ; module procedure initialized_ ; end interface

! !PUBLIC DATA MEMBERS:

      public :: Xonly ! Matrix decomposed only by COLUMN (i.e., based 
                      ! on the decomposition of x)
      public :: Yonly ! Matrix decomposed only by ROW (i.e., based 
                      ! on the decomposition of y)
      public :: XandY ! Matrix has complex ROW/COLUMN decomposed

! !SEE ALSO:
! The MCT module m_SparseMatrix for more information about Sparse Matrices.
! The MCT module m_Rearranger for deatailed information about Communications 
! scheduling.
! The MCT module m_AttrVect for details regarding the Attribute Vector.
! The MCT module m_MatAttrVectMult for documentation of API's that use 
! the SparseMatrixPlus.
!
! !REVISION HISTORY:
! 29August 2002 - J. Larson <larson@mcs.anl.gov> - API specification.
!EOP -------------------------------------------------------------------

      character(len=*), parameter :: Xonly = 'Xonly'
      character(len=*), parameter :: Yonly = 'Yonly'
      character(len=*), parameter :: XandY = 'XandY'

      character(len=*), parameter :: myname = 'm_SparseMatrixPlus'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initFromRoot_ - Creation and Initializtion from the Root
!
! !DESCRIPTION:  
! This routine creates an {\tt SparseMatrixPlus} {\tt sMatPlus} using 
! the following elements:
! \begin{itemize}
! \item A {\tt SparseMatrix} (the input argument {\tt sMat}), whose 
! elements all reside only on the {\tt root} process of the MPI 
! communicator with an integer handle defined by the input {\tt INTEGER} 
! argument {\tt comm};
! \item A {\tt GlobalSegMap} (the input argument {\tt xGSMap}) describing
! the domain decomposition of the vector {\bf x} on the communicator 
! {\tt comm};
! \item A {\tt GlobalSegMap} (the input argument {\tt yGSMap}) describing
! the domain decomposition of the vector {\bf y} on the communicator 
! {\tt comm};
! \item The matrix-vector multiplication parallelization strategy.  This
! is set by the input {\tt CHARACTER} argument {\tt strategy}, which must 
! have value corresponding to one of the following public data members 
! defined in the declaration section of this module.  Acceptable values 
! for use in this routine are: {\tt Xonly} and {\tt Yonly}.
! \end{itemize}
!
! !INTERFACE:

 subroutine initFromRoot_(sMatPlus, sMat, xGSMap, yGSMap, strategy, &
                          root, comm, ComponentID)

! !USES:

      use m_die
      use m_stdio
      use m_mpif90

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_gsize => gsize
      use m_GlobalSegMap, only : GlobalSegMap_lsize => lsize
      use m_GlobalSegMap, only : GlobalSegMap_clean => clean

      use m_SparseMatrix, only : SparseMatrix
      use m_SparseMatrix, only : SparseMatrix_init => init
      use m_SparseMatrix, only : SparseMatrix_nRows => nRows
      use m_SparseMatrix, only : SparseMatrix_nCols => nCols

      use m_SparseMatrixComms, only : SparseMatrix_ScatterByRow => ScatterByRow
      use m_SparseMatrixComms, only : SparseMatrix_ScatterByColumn => &
                                                                ScatterByColumn

      use m_SparseMatrixToMaps, only : SparseMatrixToXGlobalSegMap
      use m_SparseMatrixToMaps, only : SparseMatrixToYGlobalSegMap

      use m_GlobalToLocal, only : GlobalToLocalMatrix

      use m_Rearranger, only : Rearranger
      use m_Rearranger, only : Rearranger_init => init

      implicit none

! !INPUT PARAMETERS:

      type(GlobalSegMap),     intent(in)    :: xGSMap
      type(GlobalSegMap),     intent(in)    :: yGSMap
      character(len=*),       intent(in)    :: strategy
      integer,                intent(in)    :: root
      integer,                intent(in)    :: comm
      integer,                intent(in)    :: ComponentID

! !INPUT/OUTPUT PARAMETERS:
      
      type(SparseMatrix),     intent(inout) :: sMat

! !OUTPUT PARAMETERS:

      type(SparseMatrixPlus), intent(out)   :: SMatPlus

! !REVISION HISTORY:
! 30Aug02 - Jay Larson <larson@mcs.anl.gov> - API Specification
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::initFromRoot_'

  type(GlobalSegMap) :: xPrimeGSMap, yPrimeGSMap

  integer :: myID, ierr

       ! Set initialization flag sMatPlus%Initialized to .FALSE.
       ! therefore, any return from a failure in creation will 
       ! be visible.

  SMatPlus%Initialized = .FALSE.

       ! Get local process ID number

  call MPI_COMM_RANK(myID, comm)

       ! Basic Input Argument Checks:

       ! On the root, where the matrix is stored, do its number of 
       ! rows and columns match the global lengths ofthe vectors y 
       ! and x, respectively? 

  if(myID == root) then

     if(GlobalSegMap_gsize(yGSMap) /= SparseMatrix_nRows(sMat)) then
	write(stderr,'(3a,i8,2a,i8)') myname, &
	     ':: FATAL--length of vector y different from row count of sMat.', &
	     'Length of y = ',GlobalSegMap_gsize(yGSMap),' Number of rows in ',&
	     'sMat = ',SparseMatrix_nRows(sMat)
	call die(myname_)
     endif

     if(GlobalSegMap_gsize(xGSMap) /= SparseMatrix_nCols(sMat)) then
	write(stderr,'(3a,i8,2a,i8)') myname, &
	     ':: FATAL--length of vector x different from column count of sMat.', &
	     'Length of y = ',GlobalSegMap_gsize(yGSMap),' Number of columns in ',&
	     'sMat = ',SparseMatrix_nCols(sMat)
	call die(myname_)
     endif

  endif ! if(myID == root) then...

       ! Check desired parallelization strategy name for validity.  
       ! If either of the strategies supported by this routine are
       ! provided, initialize the appropriate component of sMatPlus.

  select case(strategy)
  case(Xonly) ! decompose sMat by rows following decomposition of y
     call String_init(sMatPlus%Strategy, strategy)
  case(Yonly) ! decompose sMat by columns following decomposition of x
     call String_init(sMatPlus%Strategy, strategy)
  case(XandY) ! User has called the wrong routine.  Try initDistributed()
              ! instead.
     write(stderr,'(4a)') myname_, &
	  ':: ERROR--Strategy name = ',strategy,' not supported by this routine.'
  case default ! strategy name not recognized.
     write(stderr,'(5a)') myname_, &
	  ':: ERROR--Invalid parallelization strategy name = ',strategy,' not ', &
	  'recognized by this module.'
     call die(myname_)
  end select

       ! End Argument Sanity Checks.

       ! Based on the parallelization strategy, scatter sMat into 
       ! sMatPlus%Matrix accordingly.

  select case(strategy)
  case(Xonly)
       ! Scatter sMat by Row
     call SparseMatrix_ScatterByRow(yGSMap, sMat, sMatPlus%Matrix, root, &
	                            comm, ierr)
       ! Compute GlobalSegMap associated with intermediate vector x'
     call SparseMatrixToXGlobalSegMap(sMatPlus%Matrix, xPrimeGSMap, &
	                              root, comm, ComponentID)
       ! Determine length of x' from xPrimeGSMap:
     sMatPlus%XPrimeLength = GlobalSegMap_lsize(xPrimeGSMap, comm)
       ! Create Rearranger to assemble x' from x
     call Rearranger_init(xGSMap, xPrimeGSMap, comm, sMatPlus%XToXPrime)
       ! Create local column indices based on xPrimeGSMap
     call GlobalToLocalMatrix(sMatPlus%Matrix, xPrimeGSMap, 'column', comm)
       ! Destroy intermediate GlobalSegMap for x'
     call GlobalSegMap_clean(xPrimeGSMap)
  case(Yonly)
       ! Scatter sMat by Column
     call SparseMatrix_ScatterByColumn(xGSMap, sMat, sMatPlus%Matrix, root, &
	                               comm, ierr)
       ! Compute GlobalSegMap associated with intermediate vector y'
     call SparseMatrixToYGlobalSegMap(sMatPlus%Matrix, yPrimeGSMap, &
	                              root, comm, ComponentID)
       ! Determine length of y' from yPrimeGSMap:
     sMatPlus%YPrimeLength = GlobalSegMap_lsize(yPrimeGSMap, comm)
       ! Create Rearranger to assemble y from partial sums in y'
     call Rearranger_init(yPrimeGSMap, yGSMap, comm, sMatPlus%YPrimeToY)
       ! Create local row indices based on yPrimeGSMap
     call GlobalToLocalMatrix(sMatPlus%Matrix, yPrimeGSMap, 'row', comm)
       ! Destroy intermediate GlobalSegMap for y'
     call GlobalSegMap_clean(yPrimeGSMap)
  case default ! do nothing
  end select

       ! Since we've succeeded, set the SparseMatrixPlus initialization
       ! flag sMatPlus%Initialized to .TRUE.

  SMatPlus%Initialized = .TRUE.

 end subroutine initFromRoot_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initDistributed_ - Creation and Initializtion from the Root
!
! !DESCRIPTION:  
! This routine creates an {\tt SparseMatrixPlus} {\tt sMatPlus} using 
! the following elements:
! \begin{itemize}
! \item A {\tt SparseMatrix} (the input argument {\tt sMat}), whose 
! elements all reside only on the {\tt root} process of the MPI 
! communicator with an integer handle defined by the input {\tt INTEGER} 
! argument {\tt comm};
! \item A {\tt GlobalSegMap} (the input argument {\tt xGSMap}) describing
! the domain decomposition of the vector {\bf x} on the communicator 
! {\tt comm};
! \item A {\tt GlobalSegMap} (the input argument {\tt yGSMap}) describing
! the domain decomposition of the vector {\bf y} on the communicator 
! {\tt comm};
! \item The matrix-vector multiplication parallelization strategy.  This
! is set by the input {\tt CHARACTER} argument {\tt strategy}, which must 
! have value corresponding to one of the following public data members 
! defined in the declaration section of this module.  Acceptable values 
! for use in this routine are: {\tt Xonly} and {\tt Yonly}.
! \end{itemize}
!
! !INTERFACE:

 subroutine initDistributed_(sMatPlus, sMat, xGSMap, yGSMap, strategy, &
                          root, comm, ComponentID)

! !USES:

      use m_die
      use m_stdio
      use m_mpif90

      use m_GlobalSegMap, only : GlobalSegMap
      use m_GlobalSegMap, only : GlobalSegMap_gsize => gsize
      use m_GlobalSegMap, only : GlobalSegMap_lsize => lsize
      use m_GlobalSegMap, only : GlobalSegMap_clean => clean

      use m_SparseMatrix, only : SparseMatrix
      use m_SparseMatrix, only : SparseMatrix_init => init
      use m_SparseMatrix, only : SparseMatrix_nRows => nRows
      use m_SparseMatrix, only : SparseMatrix_nCols => nCols
      use m_SparseMatrix, only : SparseMatrix_Copy => Copy

      use m_SparseMatrixComms, only : SparseMatrix_ScatterByRow => ScatterByRow
      use m_SparseMatrixComms, only : SparseMatrix_ScatterByColumn => &
                                                                ScatterByColumn

      use m_SparseMatrixToMaps, only : SparseMatrixToXGlobalSegMap
      use m_SparseMatrixToMaps, only : SparseMatrixToYGlobalSegMap

      use m_GlobalToLocal, only : GlobalToLocalMatrix

      use m_Rearranger, only : Rearranger
      use m_Rearranger, only : Rearranger_init => init

      implicit none

! !INPUT PARAMETERS:

      type(GlobalSegMap),     intent(in)    :: xGSMap
      type(GlobalSegMap),     intent(in)    :: yGSMap
      character(len=*),       intent(in)    :: strategy
      integer,                intent(in)    :: root
      integer,                intent(in)    :: comm
      integer,                intent(in)    :: ComponentID

! !INPUT/OUTPUT PARAMETERS:
      
      type(SparseMatrix),     intent(inout) :: sMat

! !OUTPUT PARAMETERS:

      type(SparseMatrixPlus), intent(out)   :: SMatPlus

! !REVISION HISTORY:
! 30Aug02 - Jay Larson <larson@mcs.anl.gov> - API Specification
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::initDistributed_'

  type(GlobalSegMap) :: xPrimeGSMap, yPrimeGSMap

  integer :: myID, ierr

       ! Set initialization flag sMatPlus%Initialized to .FALSE.
       ! therefore, any return from a failure in creation will 
       ! be visible.

  SMatPlus%Initialized = .FALSE.

       ! Get local process ID number

  call MPI_COMM_RANK(myID, comm)

       ! Basic Input Argument Checks:

       ! A portion of sMat (even if there are no nonzero elements in 
       ! this local chunk) on each PE.  We must check to ensure the 
       ! number rows and columns match the global lengths ofthe 
       ! vectors y and x, respectively. 

  if(GlobalSegMap_gsize(yGSMap) /= SparseMatrix_nRows(sMat)) then
     write(stderr,'(3a,i8,2a,i8)') myname, &
	  ':: FATAL--length of vector y different from row count of sMat.', &
	  'Length of y = ',GlobalSegMap_gsize(yGSMap),' Number of rows in ',&
	  'sMat = ',SparseMatrix_nRows(sMat)
     call die(myname_)
  endif

  if(GlobalSegMap_gsize(xGSMap) /= SparseMatrix_nCols(sMat)) then
     write(stderr,'(3a,i8,2a,i8)') myname, &
	  ':: FATAL--length of vector x different from column count of sMat.', &
	  'Length of y = ',GlobalSegMap_gsize(yGSMap),' Number of columns in ',&
	  'sMat = ',SparseMatrix_nCols(sMat)
     call die(myname_)
  endif

       ! Check desired parallelization strategy name for validity.  
       ! If either of the strategies supported by this routine are
       ! provided, initialize the appropriate component of sMatPlus.

  select case(strategy)
  case(Xonly,Yonly) ! Wrong routine
     write(stderr,'(5a)') myname_, &
	  ':: FATAL--inappropriate strategy "',strategy,'" for this routine.', &
	  '  Try the MCT SparseMatrixPlus routine initFromRoot().'
     call die(myname_)
  case(XandY) ! This is the strategy this routine supports
     call String_init(sMatPlus%Strategy, strategy)
  case default ! strategy name not recognized.
     write(stderr,'(5a)') myname_, &
	  ':: ERROR--Invalid parallelization strategy name = ',strategy,' not ', &
	  'recognized by this module.'
     call die(myname_)
  end select

       ! End Argument Sanity Checks.

       ! Based on the XandY parallelization strategy, build SMatPlus
       ! First, copy Internals of sMat into sMatPlus%Matrix:
  call SparseMatrix_Copy(sMat, sMatPlus%Matrix)
       ! Compute GlobalSegMap associated with intermediate vector x'
  call SparseMatrixToXGlobalSegMap(sMatPlus%Matrix, xPrimeGSMap, &
                                   root, comm, ComponentID)
       ! Determine length of x' from xPrimeGSMap:
  sMatPlus%XPrimeLength = GlobalSegMap_lsize(xPrimeGSMap, comm)
       ! Create Rearranger to assemble x' from x
  call Rearranger_init(xGSMap, xPrimeGSMap, comm, sMatPlus%XToXPrime)
       ! Create local column indices based on xPrimeGSMap
  call GlobalToLocalMatrix(sMatPlus%Matrix, xPrimeGSMap, 'column', comm)
       ! Destroy intermediate GlobalSegMap for x'
  call GlobalSegMap_clean(xPrimeGSMap)
       ! Compute GlobalSegMap associated with intermediate vector y'
  call SparseMatrixToYGlobalSegMap(sMatPlus%Matrix, yPrimeGSMap, &
	                           root, comm, ComponentID)
       ! Determine length of y' from yPrimeGSMap:
  sMatPlus%YPrimeLength = GlobalSegMap_lsize(yPrimeGSMap, comm)
       ! Create Rearranger to assemble y from partial sums in y'
  call Rearranger_init(yPrimeGSMap, yGSMap, comm, sMatPlus%YPrimeToY)
       ! Create local row indices based on yPrimeGSMap
  call GlobalToLocalMatrix(sMatPlus%Matrix, yPrimeGSMap, 'row', comm)
       ! Destroy intermediate GlobalSegMap for y'
  call GlobalSegMap_clean(yPrimeGSMap)

       ! Since we've succeeded, set the SparseMatrixPlus initialization
       ! flag sMatPlus%Initialized to .TRUE.

  SMatPlus%Initialized = .TRUE.

 end subroutine initDistributed_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: clean_ - Destruction of a SparseMatrixPlus Object
!
! !DESCRIPTION:  
! This routine deallocates all allocated memory belonging to the 
! input/output {\tt SparseMatrixPlus} argument {\tt SMatP}, and sets 
! to zero its integer components describing intermediate vector length, 
! and sets its {\tt LOGICAL} flag signifying initialization to 
! {\tt .FALSE.}  The success (failure) of this operation is signified 
! by the zero (non-zero) value of the optional {\tt INTEGER} output 
! argument {\tt status}.  If the user does supply {\tt status} when 
! invoking this routine, failure of {\tt clean\_()} will lead to 
! termination of execution with an error message.
!
! !INTERFACE:

 subroutine clean_(SMatP, status)

! !USES:

      use m_die
      use m_stdio

      use m_String, only : String
      use m_String, only : String_init => init
      use m_String, only : String_ToChar => toChar
      use m_String, only : String_clean => clean

      use m_SparseMatrix, only : SparseMatrix
      use m_SparseMatrix, only : SparseMatrix_clean => clean

      use m_Rearranger, only : Rearranger
      use m_Rearranger, only : Rearranger_clean => clean

      implicit none

! !INPUT/OUTPUT PARAMETERS:

      type(SparseMatrixPlus), intent(inout)  :: SMatP

! !OUTPUT PARAMETERS:

      integer, optional,  intent(out)   :: status

! !REVISION HISTORY:
! 30Aug02 - Jay Larson <larson@mcs.anl.gov> - API Specification
!EOP ___________________________________________________________________
!
  character(len=*),parameter :: myname_=myname//'::clean_'

  integer :: myStatus
  type(String) :: SMatP_Strategy_copy

       ! If status was supplied, set it to zero (success)

  if(present(status)) status = 0

       ! The following string copy is superfluous. It is placed here
       ! to outwit a compiler bug in the SGI and SunOS compilers. 
       ! It occurs when a component of a derived type is used as an
       ! argument to String_ToChar. This bug crashes the compiler 
       ! with the error message:
       ! Error: Signal Segmentation fault in phase IR->WHIRL Conversion

  call String_init(SMatP_Strategy_copy,SMatP%Strategy)

       ! Use SMatP%Strategy to determine which Rearranger(s) need
       ! to be destroyed.  The CHARACTER parameters Xonly, Yonly, 
       ! and XandY are inherited from the declaration section of 
       ! this module.

  select case(String_ToChar(SMatP_Strategy_copy))
  case(Xonly) ! destroy X-rearranger only

     call Rearranger_clean(SMatP%XToXprime, myStatus)
     if(myStatus /= 0) then ! something went wrong
	write(stderr,'(3a,i8)') myname_, &
	     ':: ERROR - call to Rearranger_clean(SMatP%XToXprime) failed.', &
	     ' stat = ',myStatus
	if(present(status)) then
	   status = myStatus
	   return
	else
	   call die(myname_)
	endif
     endif

  case(Yonly) ! destroy Y-rearranger only

     call Rearranger_clean(SMatP%YprimeToY, myStatus)
     if(myStatus /= 0) then ! something went wrong
	write(stderr,'(3a,i8)') myname_, &
	     ':: ERROR - call to Rearranger_clean(SMatP%YPrimeToY) failed.', &
	     ' stat = ',myStatus
	if(present(status)) then
	   status = myStatus
	   return
	else
	   call die(myname_)
	endif
     endif

  case(XandY) ! destroy both X- and Y-rearrangers

     call Rearranger_clean(SMatP%XToXprime, myStatus)
     if(myStatus /= 0) then ! something went wrong
	write(stderr,'(3a,i8)') myname_, &
	     ':: ERROR - call to Rearranger_clean(SMatP%XToXprime) failed.', &
	     ' stat = ',myStatus
	if(present(status)) then
	   status = myStatus
	   return
	else
	   call die(myname_)
	endif
     endif

     call Rearranger_clean(SMatP%YprimeToY, myStatus)
     if(myStatus /= 0) then ! something went wrong
	write(stderr,'(3a,i8)') myname_, &
	     ':: ERROR - call to Rearranger_clean(SMatP%YPrimeToY) failed.', &
	     ' stat = ',myStatus
	if(present(status)) then
	   status = myStatus
	   return
	else
	   call die(myname_)
	endif
     endif

  case default ! do nothing--corresponds to purely data local case
  end select

       ! Zero out XPrimeLength and YPrimeLength

  SMatP%XPrimeLength = 0
  SMatP%YPrimeLength = 0

       ! Destroy the SparseMatrix component SMatP%Matrix

  call SparseMatrix_clean(SMatP%Matrix, myStatus)
  if(myStatus /= 0) then ! something went wrong
     write(stderr,'(2a,i8)') myname_, &
	  ':: ERROR - call to SparseMatrix_clean() failed with stat=',myStatus
     if(present(status)) then
	status = myStatus
	return
     else
	call die(myname_)
     endif
  endif

       ! Destroy the String SMatP%Strategy and its copy

  call String_clean(SMatP%Strategy)
  call String_clean(SMatP_Strategy_copy)

       ! Set initialization flag SMatP%Initialized to .FALSE.

  SMatP%Initialized = .FALSE.

 end subroutine clean_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!    Math and Computer Science Division, Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: initialized_() - Confirmation of Initialization
!
! !DESCRIPTION:
! This {\tt LOGICAL} query function tells the user if the input 
! {\tt SparseMatrixPlus} argument {\tt sMatPlus} has been initialized.
! The return value of {\tt initialized\_} is {\tt .TRUE.} if 
! {\tt sMatPlus} has been previously initialized, {\tt .FALSE.} if it
! has not.
!
! !INTERFACE:

 logical function initialized_(sMatPlus)
!
! !USES:
!
! No external modules are used by this function.
      

      implicit none

! !INPUT PARAMETERS: 
!
      type(SparseMatrixPlus), intent(in)  :: sMatPlus

! !REVISION HISTORY:
! 26Sep02 - Jay Larson <larson@mcs.anl.gov> - Implementation
!EOP ___________________________________________________________________
!
 character(len=*),parameter :: myname_=myname//'::initialized_'

  initialized_ = sMatPlus%initialized

 end function initialized_

 end module m_SparseMatrixPlus

