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
! !INTERFACE:

 module m_MatAttrVectMul

      private   ! except

      public :: sMatAvMult        ! The master Sparse Matrix -
                                  ! Attribute Vector multipy API

    interface sMatAvMult   ; module procedure &
        sMatAvMult_xlyl_,  &
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
! !IROUTINE: sMatAvMult_xlyl_() -- Purely local matrix-vector multiply
!
! !DESCRIPTION:
!
! The sparse matrix-vector multiplication routine {\tt sMatAvMult\_xlyl\_()} 
! operates on the assumption of total data locality.  That is, the input 
! {\tt AttrVect} {\tt xaV} contains all the column data required by the 
! distributed {\tt SparseMatrix} {\tt sMat}, and the row data in {\tt sMat} 
! corresponds to the local index mapping of the output {\tt AttrVect} 
! {\tt yaV}.
!
! !INTERFACE:

 subroutine sMatAvMult_xlyl_(xaV, sMat, yaV, InterpInts)
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
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : AttrVect_indexIA => indexIA
      use m_AttrVect, only : SharedAttrIndexList

      use m_SparseMatrix, only : SparseMatrix
      use m_SparseMatrix, only : SparseMatrix_lsize => lsize
      use m_SparseMatrix, only : SparseMatrix_indexIA => indexIA
      use m_SparseMatrix, only : SparseMatrix_indexRA => indexRA

      implicit none

! !INPUT PARAMETERS:

      type(AttrVect),     intent(in)    :: xaV
      type(SparseMatrix), intent(in)    :: sMat
      logical, optional,  intent(in)    :: InterpInts

! !INPUT/OUTPUT PARAMETERS:

      type(AttrVect),     intent(inout) :: yaV

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
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sMatAvMult_xlyl_'

! Matrix element count:
  integer :: num_elements

! Matrix row, column, and weight indices:
  integer :: icol, irow, iwgt

! Overlapping attribute index number
  integer :: num_indices

! Overlapping attribute index storage arrays:
  integer, dimension(:), pointer :: xaVindices, yaVindices

! Temporary variables for multiply do-loop
  integer :: row, col
  real :: wgt

! Error flag and loop indices
  integer :: ierr, m, n

! Character variable used as a data type flag:
  character*7 :: data_flag

       ! Retrieve the number of elements in sMat:

  num_elements = SparseMatrix_lsize(sMat)

       ! Indexing the sparse matrix sMat:

  irow = SparseMatrix_indexIA(sMat,'lrow')    ! local row index
  icol = SparseMatrix_indexIA(sMat,'lcol')    ! local column index
  iwgt = SparseMatrix_indexRA(sMat,'weight')  ! weight index

       ! zero the output AttributeVector

  call AttrVect_zero(yaV)

       ! Regridding Operations:  First the REAL attributes:

  if(List_identical(xaV%rList, yaV%rList)) then ! no cross-indexing

     num_indices = List_nitem(xaV%rList)
  
       ! loop over matrix elements

     do n=1,num_elements

	row = sMat%data%iAttr(irow,n)
	col = sMat%data%iAttr(icol,n)
	wgt = sMat%data%rAttr(iwgt,n)

       ! loop over attributes being regridded.

	do m=1,num_indices

	   yaV%rAttr(m,row) = yaV%rAttr(m,row) + wgt * xaV%rAttr(m,col)

	end do ! m=1,num_indices

     end do ! n=1,num_elements

  else

     data_flag = 'REAL'
     call SharedAttrIndexList(xaV, yaV, data_flag, num_indices, &
	                      xaVindices, yaVindices)

       ! loop over matrix elements

     do n=1,num_elements

	row = sMat%data%iAttr(irow,n)
	col = sMat%data%iAttr(icol,n)
	wgt = sMat%data%rAttr(iwgt,n)

       ! loop over attributes being regridded.

	do m=1,num_indices

	   yaV%rAttr(yaVindices(m),row) = &
		yaV%rAttr(yaVindices(m),row) + &
		wgt * xaV%rAttr(xaVindices(m),col)

	end do ! m=1,num_indices

     end do ! n=1,num_elements

     deallocate(xaVindices, yaVindices, stat=ierr)
     if(ierr /= 0) call die(myname_,'first deallocate(xaVindices...',ierr)

  endif ! if(List_identical(xaV%rAttr, yaV%rAttr))...

       ! Regrid the INTEGER Attributes (if desired):

  if(present(InterpInts)) then  ! if this argument is not present, skip...

     if(InterpInts) then

	if(List_identical(xaV%iList, yaV%iList)) then ! no cross-indexing

	   num_indices = List_nitem(xaV%iList)
  
       ! loop over matrix elements

	   do n=1,num_elements

	      row = sMat%data%iAttr(irow,n)
	      col = sMat%data%iAttr(icol,n)
	      wgt = sMat%data%rAttr(iwgt,n)

	      ! loop over attributes being regridded.

	      do m=1,num_indices

		 yaV%iAttr(m,row) = yaV%iAttr(m,row) + wgt * xaV%iAttr(m,col)

	      end do ! m=1,num_indices

	   end do ! n=1,num_elements

	else ! must do attribute cross-indexing

	   data_flag = 'INTEGER'
	   call SharedAttrIndexList(xaV, yaV, data_flag, num_indices, &
		                    xaVindices, yaVindices)

       ! loop over matrix elements

	   do n=1,num_elements

	      row = sMat%data%iAttr(irow,n)
	      col = sMat%data%iAttr(icol,n)
	      wgt = sMat%data%rAttr(iwgt,n)

	      ! loop over attributes being regridded.

	      do m=1,num_indices

		 yaV%iAttr(yaVindices(m),row) = &
		      yaV%iAttr(yaVindices(m),row) + &
                      wgt * xaV%iAttr(xaVindices(m),col)

	      end do ! m=1,num_indices

	   end do ! n=1,num_elements

	   deallocate(xaVindices, yaVindices, stat=ierr)
	   if(ierr /= 0) call die(myname_,'second deallocate(xaVindices...',ierr)

	endif ! if(List_identical(xaV%iList, yaV%iList))...

     endif ! if(InterpInts)...

  endif ! if(present(InterpInts))...
     
 end subroutine sMatAvMult_xlyl_

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

 subroutine sMatAvMult_SMPlus_(xAV, sMatPlus, yAV)
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
      type(SparseMatrixPlus), intent(in)    :: sMatPlus

! !INPUT/OUTPUT PARAMETERS:

      type(AttrVect),         intent(inout) :: yAV

! !SEE ALSO:
! The MCT module m_AttrVect for more information about the AttrVect type.
! The MCT module m_SparseMatrixPlus for more information about the 
! SparseMatrixPlus type.
! !REVISION HISTORY:
! 26Sep02 - J.W. Larson <larson@mcs.anl.gov> - API specification and
!           implementation.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sMatAvMult_SMPlus_'
  type(AttrVect) :: xPrimeAV, yPrimeAV
  integer :: ierr

  select case(String_ToChar(sMatPlus%Strategy))
  case(Xonly)
       ! Create intermediate AttrVect for x'
     call AttrVect_init(xPrimeAV, xAV, sMatPlus%XPrimeLength)
       ! Rearrange data from x to get x'
     call Rearrange(xAV, xPrimeAV, sMatPlus%XToXPrime)
       ! Perform perfectly data-local multiply y = Mx'
     call sMatAvMult_xlyl_(xPrimeAV, sMatPlus%Matrix, yaV)
       ! Clean up space occupied by x'
     call AttrVect_clean(xPrimeAV, ierr)
  case(Yonly)
       ! Create intermediate AttrVect for y'
     call AttrVect_init(yPrimeAV, yAV, sMatPlus%YPrimeLength)
       ! Perform perfectly data-local multiply y' = Mx
     call sMatAvMult_xlyl_(xAV, sMatPlus%Matrix, yPrimeAV)
       ! Rearrange/reduce partial sums in y' to get y
     call Rearrange(yPrimeAV, yAV, sMatPlus%YPrimeToY, .TRUE.)
       ! Clean up space occupied by y'
     call AttrVect_clean(yPrimeAV, ierr)
  case(XandY)
       ! Create intermediate AttrVect for x'
     call AttrVect_init(xPrimeAV, xAV, sMatPlus%XPrimeLength)
       ! Create intermediate AttrVect for y'
     call AttrVect_init(yPrimeAV, yAV, sMatPlus%YPrimeLength)
       ! Rearrange data from x to get x'
     call Rearrange(xAV, xPrimeAV, sMatPlus%XToXPrime)
       ! Perform perfectly data-local multiply y' = Mx'
     call sMatAvMult_xlyl_(xPrimeAV, sMatPlus%Matrix, yPrimeAV)
       ! Rearrange/reduce partial sums in y' to get y
     call Rearrange(yPrimeAV, yAV, sMatPlus%YPrimeToY, .TRUE.)
       ! Clean up space occupied by x'
     call AttrVect_clean(xPrimeAV, ierr)
       ! Clean up space occupied by y'
     call AttrVect_clean(yPrimeAV, ierr)
  case default
     write(stderr,'(4a)') myname_, &
	  ':: ERROR--parallelization strategy name ',&
	  String_ToChar(sMatPlus%Strategy),' not supported.'
  end select

 end subroutine sMatAvMult_SMPlus_

 end module m_MatAttrVectMul




