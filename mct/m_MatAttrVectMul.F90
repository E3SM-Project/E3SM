!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Math + Computer Science Division / Argonne National Laboratory   !
!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name$ 
!BOP -------------------------------------------------------------------
!
! !MODULE: m_MatAttrVectMul - Sparse Matrix AttrVect Multipication.
!
! !DESCRIPTION:
!
! This module contains routines supporting the sparse matrix-vector 
! multiplication
! $${\bf y} = {\bf M} {\bf x},$$
! where the vectors {\bf x} and {\bf y} are stored using the MCT 
! {\tt AttrVect} datatype, and {\bf M} is stored using either the MCT 
! {\tt SparseMatrix} or {\tt SparseMatrixPlus} type.  The {\tt SparseMatrix} 
! type is used to represent {\bf M} if the multiplication process is 
! purely data-local (e.g., in a global address space, or if the process
! has been rendered embarrasingly parallel by earlier or subsequent 
! vector data redistributions).  If the multiplication process is to 
! be explicitly distributed-memory parallel, then the {\tt SparseMatrixPlus}
! type is used to store the elements of {\bf M} and all information needed
! to coordinate data redistribution and reduction of partial sums.
!
! {\bf N.B.:} The matrix-vector multiplication routines in this module 
! process only the {\bf real} attributes of the {\tt AttrVect} arguments
! corresponding to {\bf x} and {\bf y}.  They ignore the integer attributes.
!
! !INTERFACE:

 module m_MatAttrVectMul

      private   ! except

      public :: sMatAvMult        ! The master Sparse Matrix -
                                  ! Attribute Vector multipy API

    interface sMatAvMult   ; module procedure &
        sMatAvMult_DataLocal_, &
        sMatAvMult_sMPlus_
    end interface

! !SEE ALSO:
! The MCT module m_AttrVect for more information about the AttrVect type.
! The MCT module m_SparseMatrix for more information about the SparseMatrix 
! type.
! The MCT module m_SparseMatrixPlus for more details about the master class 
! for parallel sparse matrix-vector multiplication, the SparseMatrixPlus.

! !REVISION HISTORY:
! 12Jan01 - J.W. Larson <larson@mcs.anl.gov> - initial module.
! 26Sep02 - J.W. Larson <larson@mcs.anl.gov> - added high-level, distributed
!           matrix-vector multiply routine using the SparseMatrixPlus class.
!
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname='m_MatAttrVectMul'

 contains

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Math + Computer Science Division / Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sMatAvMult_DataLocal -- Purely local matrix-vector multiply
!
! !DESCRIPTION:
!
! The sparse matrix-vector multiplication routine {\tt sMatAvMult\_DataLocal\_()} 
! operates on the assumption of total data locality, which is equivalent 
! to the following two conditions:
! \begin{enumerate}
! \item The input {\tt AttrVect} {\tt xAV} contains all the values referenced 
! by the local column indices stored in the input {\tt SparsMatrix} argument 
! {\tt sMat}; and
! \item The output {\tt AttrVect} {\tt yAV} contains all the values referenced 
! by the local row indices stored in the input {\tt SparsMatrix} argument 
! {\tt sMat}.
! \end{enumerate}
! The multiplication occurs for each of the common {\tt REAL} attributes 
! shared by {\tt xAV} and {\tt yAV}.  This routine is capable of 
! cross-indexing the attributes and performing the necessary multiplications.
!
! !INTERFACE:

 subroutine sMatAvMult_DataLocal_(xAV, sMat, yAV, Vector)
!
! !USES:
!
      use m_stdio, only : stderr
      use m_die,   only : MP_perr_die, die, warn

      use m_List, only : List_identical => identical
      use m_List, only : List_nitem => nitem

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_zero => zero
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : SharedAttrIndexList

      use m_SparseMatrix, only : SparseMatrix
      use m_SparseMatrix, only : SparseMatrix_lsize => lsize
      use m_SparseMatrix, only : SparseMatrix_indexIA => indexIA
      use m_SparseMatrix, only : SparseMatrix_indexRA => indexRA
      use m_SparseMatrix, only : SparseMatrix_vecinit => vecinit

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),     intent(in)    :: xAV
      logical,optional,   intent(in)    :: Vector

! !INPUT/OUTPUT PARAMETERS:

      type(SparseMatrix), intent(inout)    :: sMat
      type(AttrVect),     intent(inout) :: yAV

! !REVISION HISTORY:
! 15Jan01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
! 10Feb01 - J.W. Larson <larson@mcs.anl.gov> - Prototype code.
! 24Apr01 - J.W. Larson <larson@mcs.anl.gov> - Modified to accomodate
!           changes to the SparseMatrix datatype.
! 25Apr01 - J.W. Larson <larson@mcs.anl.gov> - Reversed loop order
!           for cache-friendliness
! 17May01 - R. Jacob <jacob@mcs.anl.gov> - Zero the output
!           attribute vector
! 10Oct01 - J. Larson <larson@mcs.anl.gov> - Added optional LOGICAL
!           input argument InterpInts to make application of the
!           multiply to INTEGER attributes optional
! 15Oct01 - J. Larson <larson@mcs.anl.gov> - Added feature to 
!           detect when attribute lists are identical, and cross-
!           indexing of attributes is not needed.
! 29Nov01 - E.T. Ong <eong@mcs.anl.gov> - Removed MP_PERR_DIE if
!           there are zero elements in sMat. This allows for
!           decompositions where a process may own zero points.
! 29Oct03 - R. Jacob <jacob@mcs.anl.gov> - add Vector argument to
!           optionally use the vector-friendly version provided by
!           Fujitsu
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sMatAvMult_DataLocal_'

! Matrix element count:
  integer :: num_elements

! Matrix row, column, and weight indices:
  integer :: icol, irow, iwgt

! Overlapping attribute index number
  integer :: num_indices

! Overlapping attribute index storage arrays:
  integer, dimension(:), pointer :: xAVindices, yAVindices

! Temporary variables for multiply do-loop
  integer :: row, col
  real :: wgt

! Error flag and loop indices
  integer :: ierr, i, m, n, l

! Character variable used as a data type flag:
  character*7 :: data_flag

! logical flag
  logical :: usevector 

  usevector = .false.
  if(present(Vector)) then
    if(Vector) usevector = .true.
  endif

       ! Retrieve the number of elements in sMat:

  num_elements = SparseMatrix_lsize(sMat)

       ! Indexing the sparse matrix sMat:

  irow = SparseMatrix_indexIA(sMat,'lrow')    ! local row index
  icol = SparseMatrix_indexIA(sMat,'lcol')    ! local column index
  iwgt = SparseMatrix_indexRA(sMat,'weight')  ! weight index

       ! zero the output AttributeVector

  call AttrVect_zero(yAV, zeroInts=.FALSE.)

       ! Multiplication sMat by REAL attributes in xAV:

  if(List_identical(xAV%rList, yAV%rList)) then ! no cross-indexing

     num_indices = List_nitem(xAV%rList)

     if(usevector) then

       if(.not.sMat%vecinit) then
            call SparseMatrix_vecinit(sMat)
       endif

       do m=1,num_indices
          do l=1,sMat%tbl_end
!CDIR NOLOOPCHG
             do i=sMat%row_s(l),sMat%row_e(l)
               col = sMat%tcol(i,l)
               wgt = sMat%twgt(i,l)
               if (col < 0) cycle
               yAV%rAttr(m,i) = yAV%rAttr(m,i) + wgt * xAV%rAttr(m,col)
             enddo
          enddo
       enddo
 
     else

       do n=1,num_elements

    	  row = sMat%data%iAttr(irow,n)
    	  col = sMat%data%iAttr(icol,n)
	  wgt = sMat%data%rAttr(iwgt,n)

         ! loop over attributes being regridded.

	  do m=1,num_indices

	     yAV%rAttr(m,row) = yAV%rAttr(m,row) + wgt * xAV%rAttr(m,col)

  	  end do ! m=1,num_indices

       end do ! n=1,num_elements

     endif

  else

     data_flag = 'REAL'
     call SharedAttrIndexList(xAV, yAV, data_flag, num_indices, &
	                      xAVindices, yAVindices)

       ! loop over matrix elements

     do n=1,num_elements

	row = sMat%data%iAttr(irow,n)
	col = sMat%data%iAttr(icol,n)
	wgt = sMat%data%rAttr(iwgt,n)

       ! loop over attributes being regridded.

	do m=1,num_indices

	   yAV%rAttr(yAVindices(m),row) = &
		yAV%rAttr(yAVindices(m),row) + &
		wgt * xAV%rAttr(xAVindices(m),col)

	end do ! m=1,num_indices

     end do ! n=1,num_elements

     deallocate(xAVindices, yAVindices, stat=ierr)
     if(ierr /= 0) call die(myname_,'first deallocate(xAVindices...',ierr)

  endif ! if(List_identical(xAV%rAttr, yAV%rAttr))...
        ! And we are finished!

 end subroutine sMatAvMult_DataLocal_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Math + Computer Science Division / Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sMatAvMult_SMPlus_ - Parallel Multiply Using SparseMatrixPlus
!
! !DESCRIPTION:
! This routine performs distributed parallel sparse matrix-vector 
! multiplication ${\bf y} = {\bf M} {\bf x}$, where {\bf y} and
! {\bf x} are represented by the {\tt AttrVect} arguments {\tt yAV} and
! {\tt xAV}, respectively.  The matrix {\bf M} is stored in the input 
! {\tt SparseMatrixPlus} argument {\tt sMatPlus}, which also contains 
! all the information needed to coordinate the communications required to 
! gather intermediate vectors used in the multiplication process, and to 
! reduce partial sums as needed.
!
! The reader should note that the vectors in this multiplication process
! are of the MCT {\tt AttrVect} type, which means that both {\tt xAV} 
! and {\tt yAV} may contain many different vectors of the same length, 
! all bundled together.  Each of these data vectors is known as an 
! {\em attribute}, and is indexible by matching the character string tag 
! for its name.  This routine is capable of cross indexing the attributes 
! in {\tt xAV} and {\tt yAV}, and will perform the matrix-vector multiply 
! on all of the attributes they share.
!
! !INTERFACE:

 subroutine sMatAvMult_SMPlus_(xAV, sMatPlus, yAV, Vector)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90

      use m_String, only : String
      use m_String, only : String_ToChar => ToChar

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_init => init
      use m_AttrVect, only : AttrVect_clean => clean

      use m_Rearranger, only : Rearranger
      use m_Rearranger, only : Rearrange

      use m_SparseMatrixPlus, only : SparseMatrixPlus
      use m_SparseMatrixPlus, only : Xonly
      use m_SparseMatrixPlus, only : Yonly
      use m_SparseMatrixPlus, only : XandY

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),         intent(in)    :: xAV
      logical, optional,      intent(in)    :: Vector

! !INPUT/OUTPUT PARAMETERS:

      type(AttrVect),         intent(inout) :: yAV
      type(SparseMatrixPlus), intent(inout) :: sMatPlus

! !SEE ALSO:
! The MCT module m_AttrVect for more information about the AttrVect type.
! The MCT module m_SparseMatrixPlus for more information about the 
! SparseMatrixPlus type.

! !REVISION HISTORY:
! 26Sep02 - J.W. Larson <larson@mcs.anl.gov> - API specification and
!           implementation.
! 29Oct03 - R. Jacob <jacob@mcs.anl.gov> - add vector argument to all
!           calls to Rearrange and DataLocal_.  Add optional input
!           argument to change value (assumed false)
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sMatAvMult_SMPlus_'
  type(AttrVect) :: xPrimeAV, yPrimeAV
  integer :: ierr
  logical ::  usevector
  character(len=5) :: strat


  usevector = .FALSE.
  if(present(Vector)) then
   if(Vector)usevector = .TRUE.
  endif
       ! Examine the parallelization strategy, and act accordingly

  strat = String_ToChar(sMatPlus%Strategy)
  select case( strat )
  case('Xonly')
       ! Create intermediate AttrVect for x'
     call AttrVect_init(xPrimeAV, xAV, sMatPlus%XPrimeLength)
       ! Rearrange data from x to get x'
     call Rearrange(xAV, xPrimeAV, sMatPlus%XToXPrime, &
                    sMatPlus%Tag ,vector=usevector)
       ! Perform perfectly data-local multiply y = Mx'
     call sMatAvMult_DataLocal_(xPrimeAV, sMatPlus%Matrix, yaV, Vector=usevector)
       ! Clean up space occupied by x'
     call AttrVect_clean(xPrimeAV, ierr)
  case('Yonly')
       ! Create intermediate AttrVect for y'
     call AttrVect_init(yPrimeAV, yAV, sMatPlus%YPrimeLength)
       ! Perform perfectly data-local multiply y' = Mx
     call sMatAvMult_DataLocal_(xAV, sMatPlus%Matrix, yPrimeAV, Vector=usevector)
       ! Rearrange/reduce partial sums in y' to get y
     call Rearrange(yPrimeAV, yAV, sMatPlus%YPrimeToY, sMatPlus%Tag, &
                    .TRUE., Vector=usevector)
       ! Clean up space occupied by y'
     call AttrVect_clean(yPrimeAV, ierr)
  case('XandY')
       ! Create intermediate AttrVect for x'
     call AttrVect_init(xPrimeAV, xAV, sMatPlus%XPrimeLength)
       ! Create intermediate AttrVect for y'
     call AttrVect_init(yPrimeAV, yAV, sMatPlus%YPrimeLength)
       ! Rearrange data from x to get x'
     call Rearrange(xAV, xPrimeAV, sMatPlus%XToXPrime, sMatPlus%Tag, &
                       Vector=usevector)
       ! Perform perfectly data-local multiply y' = Mx'
     call sMatAvMult_DataLocal_(xPrimeAV, sMatPlus%Matrix, yPrimeAV, &
                             Vector=usevector)
       ! Rearrange/reduce partial sums in y' to get y
     call Rearrange(yPrimeAV, yAV, sMatPlus%YPrimeToY, sMatPlus%Tag, &
                            .TRUE., Vector=usevector)
       ! Clean up space occupied by x'
     call AttrVect_clean(xPrimeAV, ierr)
       ! Clean up space occupied by y'
     call AttrVect_clean(yPrimeAV, ierr)
  case default
     write(stderr,'(4a)') myname_, &
	  ':: FATAL ERROR--parallelization strategy name ',&
	  String_ToChar(sMatPlus%Strategy),' not supported.'
     call die(myname_)
  end select

 end subroutine sMatAvMult_SMPlus_

 end module m_MatAttrVectMul




