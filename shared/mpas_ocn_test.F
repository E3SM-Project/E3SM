! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_test
!
!> \brief Driver for testing MPAS ocean core
!> \author Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date   October 2013
!> \details
!>  This module contains routines to test various components of
!>  the MPAS ocean core.
!
!-----------------------------------------------------------------------

module ocn_test

   use mpas_framework
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_timer
   use mpas_tensor_operations

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_test_suite

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: hmixOn
   type (timer_node), pointer :: del2Timer, del2TensorTimer, leithTimer, del4Timer, del4TensorTimer


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_test_suite
!
!> \brief   Call all internal start-up tests
!> \author  Mark Petersen, Doug Jacobsen, Todd Ringler
!> \date    October 2013
!> \details 
!>  Call all routines to test various MPAS-Ocean components.
!
!-----------------------------------------------------------------------

   subroutine ocn_test_suite(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: err1

      err=0

      call ocn_prep_test_tensor(domain,err1)
      err = ior(err1,err)

   end subroutine ocn_test_suite!}}}

!***********************************************************************
!
!  routine ocn_prep_test_tensor
!
!> \brief   set up scratch variables to test strain rate and tensor divergence operators
!> \author  Mark Petersen
!> \date    May 2013
!> \details 
!>  This routine sets up scratch variables to test strain rate and tensor divergence operators.
!
!-----------------------------------------------------------------------

   subroutine ocn_prep_test_tensor(domain,err)!{{{

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      type (mpas_pool_type), pointer :: scratchPool, meshPool
      logical, pointer :: config_test_tensors
      character (len=StrKIND), pointer :: config_tensor_test_function

      type (field2DInteger), pointer :: edgeSignOnCellField
      type (field2DReal), pointer :: edgeTangentVectorsField
      type (field2DReal), pointer :: normalVelocityTestField
      type (field2DReal), pointer :: tangentialVelocityTestField
      type (field3DReal), pointer :: strainRateR3CellField
      type (field3DReal), pointer :: strainRateR3CellSolutionField
      type (field3DReal), pointer :: strainRateR3EdgeField
      type (field3DReal), pointer :: strainRateLonLatRCellField
      type (field3DReal), pointer :: strainRateLonLatRCellSolutionField
      type (field3DReal), pointer :: strainRateLonLatREdgeField
      type (field3DReal), pointer :: divTensorR3CellField
      type (field3DReal), pointer :: divTensorR3CellSolutionField
      type (field3DReal), pointer :: divTensorLonLatRCellField
      type (field3DReal), pointer :: divTensorLonLatRCellSolutionField
      type (field3DReal), pointer :: outerProductEdgeField

      call mpas_pool_get_config(domain % configs, 'config_test_tensors', config_test_tensors)

      if (.not.config_test_tensors) return

      call mpas_pool_get_config(domain % configs, 'config_tensor_test_function', config_tensor_test_function)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)

      call mpas_pool_get_field(meshPool, 'edgeSignOnCell', edgeSignOnCellField)
      call mpas_pool_get_field(meshPool, 'edgeTangentVectors', edgeTangentVectorsField)

      call mpas_pool_get_field(scratchPool, 'normalVelocityTest', normalVelocityTestField)
      call mpas_pool_get_field(scratchPool, 'tangentialVelocityTest', tangentialVelocityTestField)
      call mpas_pool_get_field(scratchPool, 'strainRateR3Cell', strainRateR3CellField)
      call mpas_pool_get_field(scratchPool, 'strainRateR3CellSolution', strainRateR3CellSolutionField)
      call mpas_pool_get_field(scratchPool, 'strainRateR3Edge', strainRateR3EdgeField)
      call mpas_pool_get_field(scratchPool, 'strainRateLonLatRCell', strainRateLonLatRCellField)
      call mpas_pool_get_field(scratchPool, 'strainRateLonLatRCellSolution', strainRateLonLatRCellSolutionField)
      call mpas_pool_get_field(scratchPool, 'strainRateLonLatREdge', strainRateLonLatREdgeField)
      call mpas_pool_get_field(scratchPool, 'divTensorR3Cell', divTensorR3CellField)
      call mpas_pool_get_field(scratchPool, 'divTensorR3CellSolution', divTensorR3CellSolutionField)
      call mpas_pool_get_field(scratchPool, 'divTensorLonLatRCell', divTensorLonLatRCellField)
      call mpas_pool_get_field(scratchPool, 'divTensorLonLatRCellSolution', divTensorLonLatRCellSolutionField)
      call mpas_pool_get_field(scratchPool, 'outerProductEdge', outerProductEdgeField)

      call mpas_allocate_scratch_field(normalVelocityTestField, .false.)
      call mpas_allocate_scratch_field(tangentialVelocityTestField, .false.)
      call mpas_allocate_scratch_field(strainRateR3CellField, .false.)
      call mpas_allocate_scratch_field(strainRateR3CellSolutionField, .false.)
      call mpas_allocate_scratch_field(strainRateR3EdgeField, .false.)
      call mpas_allocate_scratch_field(strainRateLonLatRCellField, .false.)
      call mpas_allocate_scratch_field(strainRateLonLatRCellSolutionField, .false.)
      call mpas_allocate_scratch_field(strainRateLonLatREdgeField, .false.)
      call mpas_allocate_scratch_field(divTensorR3CellField, .false.)
      call mpas_allocate_scratch_field(divTensorR3CellSolutionField, .false.)
      call mpas_allocate_scratch_field(divTensorLonLatRCellField, .false.)
      call mpas_allocate_scratch_field(divTensorLonLatRCellSolutionField, .false.)
      call mpas_allocate_scratch_field(outerProductEdgeField, .false.)

      call mpas_test_tensor(domain, config_tensor_test_function, &
         edgeSignOnCellField, &
         edgeTangentVectorsField, &
         normalVelocityTestField, &
         tangentialVelocityTestField, &
         strainRateR3CellField, &
         strainRateR3CellSolutionField, &
         strainRateR3EdgeField, &
         strainRateLonLatRCellField, &
         strainRateLonLatRCellSolutionField, &
         strainRateLonLatREdgeField, &
         divTensorR3CellField, &
         divTensorR3CellSolutionField, &
         divTensorLonLatRCellField, &
         divTensorLonLatRCellSolutionField, &
         outerProductEdgeField )

      call mpas_deallocate_scratch_field(normalVelocityTestField, .false.)
      call mpas_deallocate_scratch_field(tangentialVelocityTestField, .false.)
      call mpas_deallocate_scratch_field(strainRateR3CellField, .false.)
      call mpas_deallocate_scratch_field(strainRateR3CellSolutionField, .false.)
      call mpas_deallocate_scratch_field(strainRateR3EdgeField, .false.)
      call mpas_deallocate_scratch_field(strainRateLonLatRCellField, .false.)
      call mpas_deallocate_scratch_field(strainRateLonLatRCellSolutionField, .false.)
      call mpas_deallocate_scratch_field(strainRateLonLatREdgeField, .false.)
      call mpas_deallocate_scratch_field(divTensorR3CellField, .false.)
      call mpas_deallocate_scratch_field(divTensorR3CellSolutionField, .false.)
      call mpas_deallocate_scratch_field(divTensorLonLatRCellField, .false.)
      call mpas_deallocate_scratch_field(divTensorLonLatRCellSolutionField, .false.)
      call mpas_deallocate_scratch_field(outerProductEdgeField, .false.)

      err = 0

   end subroutine ocn_prep_test_tensor!}}}

end module ocn_test

! vim: foldmethod=marker
