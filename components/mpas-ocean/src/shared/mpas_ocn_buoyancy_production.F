! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_buoyancy_production
!
!> \brief MPAS ocean buoyancy production calculation 
!> \author Luke Van Roekel
!> \date   July 2021
!> \details
!>  This module contains the computation for buoyancy production 
!>  for the two equation turbulence model
!>
!
!-----------------------------------------------------------------------

module ocn_buoyancy_production

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_threading
   use mpas_constants

   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables

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

   public :: ocn_buoyancy_production_compute

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_buoyancy_production_compute
!
!> \brief   Computes Buoyancy production term  
!> \author  Luke Van Roekel 
!> \date    July 2021
!> \details
!>  This routine computes the buoyancy production tendency for the 
!>  two equation turbulence model, applies turbulence for the tke and 
!>  dissipation in the same calculation.
!
!-----------------------------------------------------------------------

   subroutine ocn_buoyancy_production_compute(tkeTend, glsPsiTend, glsPsi, &
                                              tke, c_psi3, surfaceBuoyancyForcing, &
                                              err)!{{{
      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tkeTend, glsPsiTend          !< Input/Output: velocity tendency

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         glsPsi, tke

      real (kind=RKIND), dimension(:), intent(in) :: surfaceBuoyancyForcing

      real (kind=RKIND), intent(in) :: c_psi3
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
      integer, pointer :: nVertLevels

      real (kind=RKIND) :: invAreaCell
      real (kind=RKIND) :: buoyancyProduction, diffMid 

      err = 0

      call mpas_timer_start("buoyancy production compute")

      !
      ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
      !
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(k, diffMid, c_psi3, buoyancyProduction)
      do iCell = 1, nCellsOwned
         k = 1
         diffMid = 0.5_RKIND*(vertDiffTopOfCell(k+1,iCell) + vertDiffTopOfCell(k,iCell))
         buoyancyProduction = 0.5_RKIND*(gravity*thermExpCoeff(k,iCell)*tgrad(3,k,iCell)*diffMid - &
                              gravity*salineContractCoeff(k,iCell)*sgrad(3,k,iCell)*diffMid - &
                              surfaceBuoyancyForcing(iCell)) 


         tkeTend(k,iCell) = tkeTend(k,iCell) + buoyancyProduction
         glsPsiTend(k,iCell) = glsPsiTend(k,iCell) + glsPsi(k,iCell)/(tke(k,iCell) + 1.0E-15_RKIND)* &
                               c_psi3*buoyancyProduction

         do k=2,maxLevelCell(iCell)
            diffMid = 0.5_RKIND*(vertDiffTopOfCell(k+1,iCell) + vertDiffTopOfCell(k,iCell))
            buoyancyProduction = gravity*thermExpCoeff(k,iCell)*tgrad(3,k,iCell)*diffMid - &
                                 gravity*salineContractCoeff(k,iCell)*sgrad(3,k,iCell)*diffMid

            tkeTend(k,iCell) = tkeTend(k,iCell) + buoyancyProduction
            glsPsiTend(k,iCell) = glsPsiTend(k,iCell) + glsPsi(k,iCell)/(tke(k,iCell) + 1.0E-15_RKIND)* &
                                                  c_psi3*buoyancyProduction
         end do
      end do
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("buoyancy production compute")

   !--------------------------------------------------------------------

   end subroutine ocn_buoyancy_production_compute!}}}

end module ocn_buoyancy_production

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
