!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Math + Computer Science Division / Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !MODULE: m_MatAttrVectMul - Sparse Matrix AttrVect Multipication.
!
! !DESCRIPTION:
!
! This module contains implementations of numerous parallel Sparse Matrix
! Attribute Vector multiplication routines.
!
! !INTERFACE:

 module m_MatAttrVectMul

      private   ! except

      public :: sMatAvMult        ! The master Sparse Matrix -
                                  ! Attribute Vector multipy API

    interface sMatAvMult   ; module procedure &
        sMatAvMult_xlyl_,  &
        sMatAvMult_gm_xdyl_,  &
        sMatAvMult_gsm_xdyl_,  &
        sMatAvMult_gm_xlyd_,  &
        sMatAvMult_gsm_xlyd_,  &
        sMatAvMult_gm_xdyd_,  &
        sMatAvMult_gsm_xdyd_ 
    end interface

! !REVISION HISTORY:
!       12Jan01 - J.W. Larson <larson@mcs.anl.gov> - initial module
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
      use m_die,   only : MP_perr_die

      use m_List, only : List_identical => identical
      use m_List, only : List_nitem => nitem

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_zero => zero
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : AttrVect_indexIA => indexIA

      use m_SparseMatrix, only : SparseMatrix
      use m_SparseMatrix, only : SparseMatrix_lsize => lsize
      use m_SparseMatrix, only : SparseMatrix_indexIA => indexIA
      use m_SparseMatrix, only : SparseMatrix_indexRA => indexRA

      use m_SharedAttrIndices, only : SharedAttrIndexList

      implicit none

      type(AttrVect),     intent(in)    :: xaV
      type(SparseMatrix), intent(in)    :: sMat
      logical, optional,  intent(in)    :: InterpInts

      type(AttrVect),     intent(inout) :: yaV

! !REVISION HISTORY:
!       15Jan01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!       10Feb01 - J.W. Larson <larson@mcs.anl.gov> - Prototype code.
!       24Apr01 - J.W. Larson <larson@mcs.anl.gov> - Modified to accomodate
!                 changes to the SparseMatrix datatype.
!       25Apr01 - J.W. Larson <larson@mcs.anl.gov> - Reversed loop order
!                 for cache-friendliness
!       17May01 - R. Jacob <jacob@mcs.anl.gov> - Zero the output
!                 attribute vector
!       10Oct01 - J. Larson <larson@mcs.anl.gov> - Added optional LOGICAL
!                 input argument InterpInts to make application of the
!                 multiply to INTEGER attributes optional
!       15Oct01 - J. Larson <larson@mcs.anl.gov> - Added feature to 
!                 detect when attribute lists are identical, and cross-
!                 indexing of attributes is not needed.
!       29Nov01 - E.T. Ong <eong@mcs.anl.gov> - Removed MP_PERR_DIE if
!                 there are zero elements in sMat. This allows for
!                 decompositions where a process may own zero points.
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

  if(num_elements == 0) then
     write(stderr,'(2a)') myname_, &
          ":: Warning: Zero elements in SparseMatrix sMat."
  endif

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
     if(ierr /= 0) then
	call MP_perr_die(myname_,'first deallocate(xaVindices...',ierr)
     endif

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
	   if(ierr /= 0) then
	      call MP_perr_die(myname_,'second deallocate(xaVindices...', &
		               ierr)
	   endif

	endif ! if(List_identical(xaV%iList, yaV%iList))...

     endif ! if(InterpInts)...

  endif ! if(present(InterpInts))...
     
 end subroutine sMatAvMult_xlyl_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Math + Computer Science Division / Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sMatAvMult_gm_xdyl_() -- Multiply, x GlobalMap distributed.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine sMatAvMult_gm_xdyl_(xaV, xGMap, sMat, yaV)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90
      use m_GlobalMap, only : GlobalMap
      use m_AttrVect, only : AttrVect
      use m_SparseMatrix, only : SparseMatrix

      implicit none

      type(AttrVect),     intent(in)    :: xaV
      type(GlobalMap),    intent(in)    :: xGMap
      type(SparseMatrix), intent(in)    :: sMat
      type(AttrVect),     intent(inout) :: yaV


! !REVISION HISTORY:
!       15Jan01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sMatAvMult_gm_xdyl_'

 end subroutine sMatAvMult_gm_xdyl_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Math + Computer Science Division / Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sMatAvMult_gsm_xdyl_() -- Multiply, x GlobalSegMap distributed.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine sMatAvMult_gsm_xdyl_(xaV, xGSMap, sMat, yaV)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90
      use m_GlobalSegMap, only : GlobalSegMap
      use m_AttrVect, only : AttrVect
      use m_SparseMatrix, only : SparseMatrix

      implicit none

      type(AttrVect)   ,  intent(in)    :: xaV
      type(GlobalSegMap), intent(in)    :: xGSMap
      type(SparseMatrix), intent(in)    :: sMat
      type(AttrVect),     intent(inout) :: yaV

! !REVISION HISTORY:
!       15Jan01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sMatAvMult_gsm_xdyl_'

 end subroutine sMatAvMult_gsm_xdyl_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Math + Computer Science Division / Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sMatAvMult_gm_xlyd_() -- Multiply, y GlobalMap distributed.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine sMatAvMult_gm_xlyd_(xaV, sMat, yaV, yGMap)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90
      use m_GlobalMap, only : GlobalMap
      use m_AttrVect, only : AttrVect
      use m_SparseMatrix, only : SparseMatrix

      implicit none

      type(AttrVect),     intent(in)    :: xaV
      type(SparseMatrix), intent(in)    :: sMat
      type(AttrVect),     intent(inout) :: yaV
      type(GlobalMap),    intent(in)    :: yGMap


! !REVISION HISTORY:
!       15Jan01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sMatAvMult_gm_xlyd_'

 end subroutine sMatAvMult_gm_xlyd_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Math + Computer Science Division / Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sMatAvMult_gsm_xlyd_() -- Multiply, x GlobalSegMap distributed.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine sMatAvMult_gsm_xlyd_(xaV, sMat, yaV, yGSMap)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90
      use m_GlobalSegMap, only : GlobalSegMap
      use m_AttrVect, only : AttrVect
      use m_SparseMatrix, only : SparseMatrix

      implicit none

      type(AttrVect)   ,  intent(in)    :: xaV
      type(SparseMatrix), intent(in)    :: sMat
      type(AttrVect),     intent(inout) :: yaV
      type(GlobalSegMap), intent(in)    :: yGSMap


! !REVISION HISTORY:
!       15Jan01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sMatAvMult_gsm_xlyd_'

 end subroutine sMatAvMult_gsm_xlyd_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Math + Computer Science Division / Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sMatAvMult_gm_xdyd_() -- Multiply, x, y GlobalMap distributed.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine sMatAvMult_gm_xdyd_(xaV, xGMap, sMat, yaV, yGMap)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90
      use m_GlobalMap, only : GlobalMap
      use m_AttrVect, only : AttrVect
      use m_SparseMatrix, only : SparseMatrix

      implicit none

      type(AttrVect),     intent(in)    :: xaV
      type(GlobalMap),    intent(in)    :: xGMap
      type(SparseMatrix), intent(in)    :: sMat
      type(AttrVect),     intent(inout) :: yaV
      type(GlobalMap),    intent(in)    :: yGMap


! !REVISION HISTORY:
!       15Jan01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sMatAvMult_gm_xdyd_'

 end subroutine sMatAvMult_gm_xdyd_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Math + Computer Science Division / Argonne National Laboratory   !
!BOP -------------------------------------------------------------------
!
! !IROUTINE: sMatAvMult_gsm_xdyd_() -- Mult., x,y GlobalSegMap distributed.
!
! !DESCRIPTION:
!
! !INTERFACE:

 subroutine sMatAvMult_gsm_xdyd_(xaV, xGSMap, sMat, yaV, yGSMap)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90
      use m_GlobalSegMap, only : GlobalSegMap
      use m_AttrVect, only : AttrVect
      use m_SparseMatrix, only : SparseMatrix

      implicit none

      type(AttrVect)   ,  intent(in)    :: xaV
      type(GlobalSegMap), intent(in)    :: xGSMap
      type(SparseMatrix), intent(in)    :: sMat
      type(AttrVect),     intent(inout) :: yaV
      type(GlobalSegMap), intent(in)    :: yGSMap


! !REVISION HISTORY:
!       15Jan01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sMatAvMult_gsm_xdyd_'

 end subroutine sMatAvMult_gsm_xdyd_

 end module m_MatAttrVectMul
