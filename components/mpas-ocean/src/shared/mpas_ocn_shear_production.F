! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_shear_production
!
!> \brief MPAS two equation turbulence model shear production
!> \author Luke Van Roekel 
!> \date   June 2021
!> \details
!>  This module contains the routines for computing shear production term 
!>  for the two equation turbulence closure
!>
!
!-----------------------------------------------------------------------

module ocn_shear_production

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_threading

   use ocn_constants
   use ocn_diagnostics_variables
   use ocn_config
   use ocn_mesh

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

   public :: ocn_shear_production_compute

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_shear_production_compute
!
!> \brief   Computes shear production term for the two equation model 
!> \author  Luke Van Roekel
!> \date    July 2021
!> \details
!>  This routine computes the shear production term for the  
!>  two equation model.
!
!-----------------------------------------------------------------------

   subroutine ocn_shear_production_compute(tkeTend, glsPsiTend, tke, glsPsi, &
                  c_psi1, surfaceFrictionVelocity, c_mu, err)
      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tkeTend,glsPsiTend          !< Input/Output:  tendencies

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         tke, glsPsi          !< Input: tke and psi variable at current timestep 

      real (kind=RKIND), dimension(:), intent(in) :: &
         surfaceFrictionVelocity 

      real (kind=RKIND), intent(in) :: c_psi1, c_mu

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

      integer :: iCell, iEdge, cell1, cell2
      integer :: i, k, nCells

      real (kind=RKIND) :: invAreaCell
      real (kind=RKIND) :: shearProduction, viscMid

      err = 0

      call mpas_timer_start("shear production compute")

      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(k, c_psi2, viscMid, shearProduction)
      do iCell = 1, nCellsOwned
        do k=1,maxLevelCell(iCell)
            viscMid = 0.5_RKIND*(vertViscTopOfCell(k,iCell) + vertViscTopOfCell(k+1,iCell))
            shearProduction = 2.0_RKIND*(viscMid*ugrad(1,k,iCell) - tke(k,iCell)/3.0_RKIND)* &
                              ugrad(1,k,iCell) + viscMid*(ugrad(2,k,iCell) + vgrad(1,k,iCell))* &
                              ugrad(2,k,iCell) + viscMid*(ugrad(3,k,iCell) + wgrad(1,k,iCell))* &
                              ugrad(3,k,iCell) + viscMid*(ugrad(2,k,iCell) + vgrad(1,k,iCell))* &
                              vgrad(1,k,iCell) + 2.0_RKIND*(viscMid*vgrad(2,k,iCell) -          &
                              tke(k,iCell)/3.0_RKIND)*vgrad(2,k,iCell) + viscMid*(vgrad(3,k,iCell) + &
                              wgrad(2,k,iCell))*vgrad(3,k,iCell) + viscMid*(ugrad(3,k,iCell) +  &
                              wgrad(1,k,iCell))*wgrad(1,k,iCell) + viscMid*(vgrad(3,k,iCell) +  &
                              wgrad(2,k,iCell))*wgrad(2,k,iCell) + 2.0_RKIND*(viscMid*wgrad(3,k,iCell) - &
                              tke(k,iCell)/3.0_RKIND)*wgrad(3,k,iCell)
            tkeTend(k,iCell) = tkeTend(k,iCell) + shearProduction
            glsPsiTend(k,iCell) = glsPsiTend(k,iCell) + glsPsi(k,iCell)/(tke(k,iCell) + 1.0E-15_RKIND) &
                                  *c_psi1*shearProduction
         end do 
         !TODO - need a bottom drag term here.
      end do
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("shear production compute")

   !--------------------------------------------------------------------

   end subroutine ocn_shear_production_compute

end module ocn_shear_production

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
