!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
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
!
! !USES:
!
      use m_AttrVect,     only : AttrVect
      use m_SparseMatrix, only : SparseMatrix

      use m_GlobalMap,    only : GlobalMap
      use m_GlobalSegMap, only : GlobalSegMap

      implicit none

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
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
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

 subroutine sMatAvMult_xlyl_(xaV, sMat, yaV)
!
! !USES:
!
      use m_stdio, only : stderr
      use m_die,   only : MP_perr_die

      use m_AttrVect, only : AttrVect
      use m_AttrVect, only : AttrVect_lsize => lsize
      use m_AttrVect, only : AttrVect_nIAttr => nIAttr
      use m_AttrVect, only : AttrVect_nRAttr => nRAttr
      use m_AttrVect, only : AttrVect_indexRA => indexRA
      use m_AttrVect, only : AttrVect_indexIA => indexIA

      use m_SparseMatrix, only : SparseMatrix
      use m_SharedAttrIndices, only : SharedAttrIndexList

      implicit none

      type(AttrVect),     intent(in  )  :: xaV
      type(AttrVect),     intent(inout) :: yaV
      type(SparseMatrix), intent(in)    :: sMat

! !REVISION HISTORY:
!       15Jan01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!       10Feb01 - J.W. Larson <larson@mcs.anl.gov> - Prototype code.
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
  integer :: xaVindex, yaVindex

! Error flag and loop indices
  integer :: ierr, m, n

! Character variable used as a data type flag:
  character*7 :: data_flag

       ! Retrieve the number of elements in sMat:

  num_elements = AttrVect_lsize(sMat)

  if(num_elements == 0) then
     write(stderr,'(2a)') myname_, &
          ":: Zero elements in SparseMatrix sMat."
     ierr = 1
     call MP_perr_die(myname_,'Zero steps in sMat.',ierr)
  endif

       ! Indexing the sparse matrix sMat:
  irow = AttrVect_indexIA(sMat,'row')  ! row index
  icol = AttrVect_indexIA(sMat,'col')  ! column index
  iwgt = AttrVect_indexRA(sMat,'col')  ! weight index

       ! Regridding Operations:  First the REAL attributes:

  data_flag = 'REAL'
  call SharedAttrIndexList(xaV, yaV, data_flag, num_indices, &
                              xaVindices, yaVindices)

       ! loop over attributes being regridded.

  do m=1,num_indices

     xaVindex = xaVindices(m)
     yaVindex = yaVindices(m)

       ! loop over matrix elements

     do n=1,num_elements

	yaV%rAttr(yaVindex,sMat%iAttr(irow,n)) = sMat%rAttr(iwgt, n) &
	     * xaV%rAttr(xaVindex,sMat%iAttr(icol,n))

     end do

  end do

  deallocate(xaVindices, yaVindices, stat=ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,'first deallocate(xaVindices...',ierr)
  endif

       ! Regrid the INTEGER Attributes:

  data_flag = 'INTEGER'
  call SharedAttrIndexList(xaV, yaV, data_flag, num_indices, &
                              xaVindices, yaVindices)

       ! loop over attributes being regridded.

  do m=1,num_indices

     xaVindex = xaVindices(m)
     yaVindex = yaVindices(m)

       ! loop over matrix elements

     do n=1,num_elements

	yaV%iAttr(yaVindex,sMat%iAttr(irow,n)) = sMat%rAttr(iwgt, n) &
	     * xaV%iAttr(xaVindex,sMat%iAttr(icol,n))

     end do

  end do

  deallocate(xaVindices, yaVindices, stat=ierr)
  if(ierr /= 0) then
     call MP_perr_die(myname_,'second deallocate(xaVindices...',ierr)
  endif
     

 end subroutine sMatAvMult_xlyl_

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
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
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
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
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
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
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
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
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
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
!       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
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
