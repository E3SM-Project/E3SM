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
! !INTERFACE:

 subroutine sMatAvMult_xlyl_(xaV, sMat, yaV)
!
! !USES:
!
      use m_stdio
      use m_die
      use m_mpif90
      use m_AttrVect, only : AttrVect
      use m_SparseMatrix, only : SparseMatrix

      implicit none

      type(AttrVect),     intent(in  )  :: xaV
      type(AttrVect),     intent(inout) :: yaV
      type(SparseMatrix), intent(in)    :: sMat

! !REVISION HISTORY:
!       15Jan01 - J.W. Larson <larson@mcs.anl.gov> - API specification.
!EOP ___________________________________________________________________

  character(len=*),parameter :: myname_=myname//'::sMatAvMult_xlyl_'

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
