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

   use mpas_configure
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

      if (.not.config_test_tensors) return

      call mpas_allocate_scratch_field(domain % blocklist % scratch % normalVelocityTest, .false.)
      call mpas_allocate_scratch_field(domain % blocklist % scratch % tangentialVelocityTest, .false.)
      call mpas_allocate_scratch_field(domain % blocklist % scratch % strainRateR3Cell, .false.)
      call mpas_allocate_scratch_field(domain % blocklist % scratch % strainRateR3CellSolution, .false.)
      call mpas_allocate_scratch_field(domain % blocklist % scratch % strainRateR3Edge, .false.)
      call mpas_allocate_scratch_field(domain % blocklist % scratch % strainRateLonLatRCell, .false.)
      call mpas_allocate_scratch_field(domain % blocklist % scratch % strainRateLonLatRCellSolution, .false.)
      call mpas_allocate_scratch_field(domain % blocklist % scratch % strainRateLonLatREdge, .false.)
      call mpas_allocate_scratch_field(domain % blocklist % scratch % divTensorR3Cell, .false.)
      call mpas_allocate_scratch_field(domain % blocklist % scratch % divTensorR3CellSolution, .false.)
      call mpas_allocate_scratch_field(domain % blocklist % scratch % divTensorLonLatRCell, .false.)
      call mpas_allocate_scratch_field(domain % blocklist % scratch % divTensorLonLatRCellSolution, .false.)
      call mpas_allocate_scratch_field(domain % blocklist % scratch % outerProductEdge, .false.)

      call mpas_test_tensor(domain, config_tensor_test_function, &
         domain % blocklist % mesh % edgeSignOnCell, &
         domain % blocklist % mesh % edgeTangentVectors, &
         domain % blocklist % scratch % normalVelocityTest, &
         domain % blocklist % scratch % tangentialVelocityTest, &
         domain % blocklist % scratch % strainRateR3Cell, &
         domain % blocklist % scratch % strainRateR3CellSolution, &
         domain % blocklist % scratch % strainRateR3Edge, &
         domain % blocklist % scratch % strainRateLonLatRCell, &
         domain % blocklist % scratch % strainRateLonLatRCellSolution, &
         domain % blocklist % scratch % strainRateLonLatREdge, &
         domain % blocklist % scratch % divTensorR3Cell, &
         domain % blocklist % scratch % divTensorR3CellSolution, &
         domain % blocklist % scratch % divTensorLonLatRCell, &
         domain % blocklist % scratch % divTensorLonLatRCellSolution, &
         domain % blocklist % scratch % outerProductEdge )

      call mpas_deallocate_scratch_field(domain % blocklist % scratch % normalVelocityTest, .false.)
      call mpas_deallocate_scratch_field(domain % blocklist % scratch % tangentialVelocityTest, .false.)
      call mpas_deallocate_scratch_field(domain % blocklist % scratch % strainRateR3Cell, .false.)
      call mpas_deallocate_scratch_field(domain % blocklist % scratch % strainRateR3CellSolution, .false.)
      call mpas_deallocate_scratch_field(domain % blocklist % scratch % strainRateR3Edge, .false.)
      call mpas_deallocate_scratch_field(domain % blocklist % scratch % strainRateLonLatRCell, .false.)
      call mpas_deallocate_scratch_field(domain % blocklist % scratch % strainRateLonLatRCellSolution, .false.)
      call mpas_deallocate_scratch_field(domain % blocklist % scratch % strainRateLonLatREdge, .false.)
      call mpas_deallocate_scratch_field(domain % blocklist % scratch % divTensorR3Cell, .false.)
      call mpas_deallocate_scratch_field(domain % blocklist % scratch % divTensorR3CellSolution, .false.)
      call mpas_deallocate_scratch_field(domain % blocklist % scratch % divTensorLonLatRCell, .false.)
      call mpas_deallocate_scratch_field(domain % blocklist % scratch % divTensorLonLatRCellSolution, .false.)
      call mpas_deallocate_scratch_field(domain % blocklist % scratch % outerProductEdge, .false.)

      err = 0

   end subroutine ocn_prep_test_tensor!}}}

end module ocn_test

! vim: foldmethod=marker
