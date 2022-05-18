! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tke_glsPsi_transport
!
!> \brief MPAS ocean tke transport calculation 
!> \author Luke Van Roekel 
!> \date   July 2021
!> \details
!>  This module contains the main driver routine for computing
!>  the diffusive transport of TKE for the two equation model 
!>
!
!-----------------------------------------------------------------------

module ocn_tke_glsPsi_transport

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_threading

   use ocn_diagnostics_variables
   use ocn_mesh
   use ocn_constants
   use ocn_config

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

   public :: ocn_tke_transport_horizontal_compute, &
             ocn_tke_transport_vertical_compute,   &
             ocn_glsPsi_transport_horizontal_compute, &
             ocn_glsPsi_transport_vertical_compute

   !--------------------------------------------------------------------
   !
   ! Private member functions
   !
   !--------------------------------------------------------------------

   private :: turb_tridiagonal_solve

   real(kind=RKIND) :: eddyDiff2

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tke_transport_horizontal_compute
!
!> \brief   Computes horizontal Laplacian tendency term for tke 
!> \author  Luke Van Roekel
!> \date    July 2021
!> \details
!>  This routine computes the horizontal mixing tendency for tke 
!>  this is a parameterized form of the turbulent transport term
!
!-----------------------------------------------------------------------

   subroutine ocn_tke_transport_horizontal_compute(tke, tkeTend, sigma_k, err)
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        tke !< Input: tke at current timestep
      real (kind=RKIND), intent(in) :: sigma_k

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tkeTend          !< Input/Output: tke tendency

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
      integer :: j, i, k, iCellValid

      real (kind=RKIND) :: invAreaCell
      real (kind=RKIND) :: tke_turb_flux, flux, r_tmp

      err = 0

      call mpas_timer_start("tke horizontal turbulent transport")

      !
      ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
      !
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(invAreaCell, i, iEdge, cell1, cell2, r_tmp, k, &
      !$omp         tke_turb_flux, flux)
      do iCell = 1, nCellsOwned
        invAreaCell = 1.0_RKIND / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          r_tmp = dvEdge(iEdge) / dcEdge(iEdge)

          do k = 1, maxLevelEdgeTop(iEdge)
             ! \kappa_2 \nabla \phi on edge
             tke_turb_flux = tke(k, cell2) - tke(k, cell1)

             eddyDiff2 = 0.25_RKIND*(vertDiffTopOfCell(k,cell1) + vertDiffTopOfCell(k,cell2) + &
                                     vertDiffTopOfCell(k+1,cell1) + vertDiffTopOfCell(k+1,cell2))/sigma_k
             eddyDiff2 = eddyDiff2 + 0.25_RKIND*(dynamicViscosity(k,cell1) + dynamicViscosity(k,cell2) + &
                                                 dynamicViscosity(k+1,cell1) + dynamicViscosity(k+1,cell2))/rho_sw

             ! div(h \kappa_2 \nabla \phi) at cell center
             flux = eddyDiff2 * layerThickEdge(k, iEdge) * tke_turb_flux * r_tmp

             tkeTend(k, iCell) = tkeTend(k, iCell) - edgeSignOnCell(i, iCell) * flux * invAreaCell
          end do

        end do
      end do
      !$omp end do
      !$omp end parallel

      !Implements the horizontal boundary condition
      if(config_two_equation_model_choice == 'epsilon') then
         call mpas_log_write('BCs not yet implemented for k-eps scheme, choose a different model',MPAS_LOG_CRIT)
      else
         !Go through find boundary cell = 1 and add -tke(icell)*eddyDiff2*....

         do iCell=1,nCellsOwned
            do k=1,maxLevelCell(iCell)
               if(boundaryCellMask(k,iCell) == 1) then
                  do j=1,nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(j,iCell)
                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)
                     if(cell1 < 1) then
                        iCellValid = cell2
                     else
                        iCellValid = cell1
                     endif

                     eddyDiff2 = 0.25_RKIND*(vertDiffTopOfCell(k,iCellValid) + vertDiffTopOfCell(k+1,iCellValid))
                     eddyDiff2 = eddyDiff2+0.5_RKIND*(dynamicViscosity(k,iCellValid) + dynamicViscosity(k+1,iCellValid))/rho_sw
                     tke_turb_flux = -tke(k,iCellValid) !Assumes k=0 in boundary condition
                     ! div(h \kappa_2 \nabla \phi) at cell center
                     flux = eddyDiff2 * layerThickEdge(k, iEdge) * tke_turb_flux * dvEdge(iEdge)/dcEdge(iEdge)

                     tkeTend(k, iCellValid) = tkeTend(k, iCellValid) - edgeSignOnCell(i, iCellValid) * flux * invAreaCell
                  end do !nEdges loop
               end if !end cell mask if
            end do
         end do !nCells
      end if !config_two_equation_model_choice

      call mpas_timer_stop("tke horizontal turbulent transport")

   !--------------------------------------------------------------------

   end subroutine ocn_tke_transport_horizontal_compute

!***********************************************************************
!
!  routine ocn_tke_transport_vertical_compute
!
!> \brief   Compute vertical component of parameterized tke transport
!> \author  Luke Van Roekel
!> \date    July 2021
!> \details
!>  This routine initializes a variety of quantities related to
!>  Laplacian horizontal velocity mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tke_transport_vertical_compute(tke,layerThickness,sigma_k,dt)

      real(kind=RKIND), dimension(:,:), intent(inout) :: tke
      real(kind=RKIND), dimension(:,:), intent(in) :: layerThickness
      real(kind=RKIND), intent(in) :: dt
      integer :: iCell, k, N
      real (kind=RKIND), dimension(nVertLevels) :: A,B,C,tkeTmp,rhs
      real (kind=RKIND), intent(in) :: sigma_k

      call mpas_timer_start("tke vertical turbulent transport")

      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(N, A, B, C, rhs, tkeTmp, k)
      do iCell=1,nCellsOwned
         N = maxLevelCell(iCell)
         A(1) = 0.0_RKIND
         do k = 2, N
            A(k) = -2.0_RKIND*dt*(dynamicViscosity(k,iCell)/rho_sw+vertViscTopOfCell(k,iCell)/sigma_k) &
                  / (layerThickness(k-1,iCell) + layerThickness(k,iCell)) &
                  / layerThickness(k,iCell)
         end do

         do k = 1, N-1
            C(k) = -2.0_RKIND*dt*(dynamicViscosity(k+1,iCell)/rho_sw+vertViscTopOfCell(k+1,iCell)/sigma_k) &
                  / (layerThickness(k,iCell) + layerThickness(k+1,iCell)) &
                  / layerThickness(k,iCell)
         end do
         C(N) = 0.0_RKIND

         do k = 1, N
            B(k) = 1.0_RKIND - C(k) - A(k)
         end do

         rhs(:) = tke(:,iCell)
         tkeTmp(:) = 0.0_RKIND
         call turb_tridiagonal_solve(A(2:N), B, C(1:N-1), rhs(:), tkeTmp, &
                                N)

         tke(1:N,iCell) = tkeTmp(1:N)
         tke(N+1:nVertLevels,iCell) = 0.0_RKIND
      end do
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("tke vertical turbulent transport")

   end subroutine ocn_tke_transport_vertical_compute

!***********************************************************************
!
!  routine ocn_glsPsi_transport_horizontal_compute
!
!> \brief   Computes horizontal Laplacian tendency term for glsPsi 
!> \author  Luke Van Roekel
!> \date    July 2021
!> \details
!>  This routine computes the horizontal mixing tendency for glsPsi 
!>  this is a parameterized form of the turbulent transport term
!
!-----------------------------------------------------------------------

   subroutine ocn_glsPsi_transport_horizontal_compute(glsPsi, glsPsiTend, &
                    c_psi2, sigma_psi, err)
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        glsPsi !< Input: tke at current timestep
      real (kind=RKIND), intent(in) :: sigma_psi
      real (kind=RKIND) :: c_psi2

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         glsPsiTend          !< Input/Output: tke tendency

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

      integer :: iCell, iEdge, cell1, cell2, iCellValid
      integer :: j, i, k

      real (kind=RKIND) :: bcValue, invAreaCell
      real (kind=RKIND) :: glsPsi_turb_flux, flux, r_tmp

      err = 0

      call mpas_timer_start("glsPsi horizontal turbulent transport")

      !
      ! compute a boundary mask to enforce insulating boundary conditions in the horizontal
      !
      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(invAreaCell, i, iEdge, cell1, cell2, r_tmp, k, &
      !$omp         glsPsi_turb_flux, flux)
      do iCell = 1, nCellsOwned
        invAreaCell = 1.0_RKIND / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          cell1 = cellsOnEdge(1,iEdge)
          cell2 = cellsOnEdge(2,iEdge)

          r_tmp = dvEdge(iEdge) / dcEdge(iEdge)

          do k = 1, maxLevelEdgeTop(iEdge)
             ! \kappa_2 \nabla \phi on edge
             glsPsi_turb_flux = glsPsi(k, cell2) - glsPsi(k, cell1)

             eddyDiff2 = 0.25_RKIND*(vertDiffTopOfCell(k,cell1) + vertDiffTopOfCell(k,cell2) + &
                                     vertDiffTopOfCell(k+1,cell1) + vertDiffTopOfCell(k+1,cell2))/sigma_psi
             eddyDiff2 = eddyDiff2+0.25_RKIND*(dynamicViscosity(k,cell1)   + dynamicViscosity(k,cell2) + &
                                               dynamicViscosity(k+1,cell1) + dynamicViscosity(k+1,cell2))/rho_sw

             ! div(h \kappa_2 \nabla \phi) at cell center
             flux = eddyDiff2 * layerThickEdge(k, iEdge) * glsPsi_turb_flux * r_tmp

             glsPsiTend(k, iCell) = glsPsiTend(k, iCell) - edgeSignOnCell(i, iCell) * flux * invAreaCell
          end do

        end do
      end do
      !$omp end do
      !$omp end parallel

      !Implements the horizontal boundary condition
      if(config_two_equation_model_choice == 'epsilon') then
         call mpas_log_write('BCs not yet implemented for k-eps scheme, choose a different model',MPAS_LOG_CRIT)
      else if(config_two_equation_model_choice == 'omega') then
         !Go through find boundary cell = 1 and add -tke(icell)*eddyDiff2*....

         do iCell=1,nCellsOwned
            do k=1,maxLevelCell(iCell)
               if(boundaryCellMask(k,iCell) == 1) then
                  do j=1,nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(j,iCell)
                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)
                     if(cell1 < 1 .or. cell2 < 1) then
                        if(cell1 < 1) then
                           iCellValid = cell2
                        else
                           iCellValid = cell1
                        endif
                        eddyDiff2 = 0.25_RKIND*(vertDiffTopOfCell(k,iCellValid) + vertDiffTopOfCell(k+1,iCellValid))
                        eddyDiff2 = eddyDiff2+0.5_RKIND*(dynamicViscosity(k,iCellValid) + dynamicViscosity(k+1,iCellValid))/rho_sw
                        bcValue = 60.0_RKIND*dynamicViscosity(k,iCellValid)/rho_sw/(c_psi2*distanceToBoundary(k,iCellValid))
                        glsPsi_turb_flux = bcValue - glsPsi(k,iCellValid) !Assumes k=0 in boundary condition
                        ! div(h \kappa_2 \nabla \phi) at cell center
                        flux = eddyDiff2 * layerThickEdge(k, iEdge) * glsPsi_turb_flux * dvEdge(iEdge)/dcEdge(iEdge)

                        glsPsiTend(k, iCellValid) = glsPsiTend(k, iCellValid) - flux * invAreaCell
                     end if
                  end do !nEdges loop
               end if !end cell mask if
            end do
         end do !nCells
      else !kl and tau schemes
         !Go through find boundary cell = 1 and add -tke(icell)*eddyDiff2*....

         do iCell=1,nCellsOwned
            do k=1,maxLevelCell(iCell)
               if(boundaryCellMask(k,iCell) == 1) then
                  do j=1,nEdgesOnCell(iCell)
                     iEdge = edgesOnCell(j,iCell)
                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnEdge(2,iEdge)
                     if(cell1 < 1 .or. cell2 < 1) then
                        eddyDiff2 = 0.25_RKIND*(vertDiffTopOfCell(k,iCell) + vertDiffTopOfCell(k+1,iCell))
                        glsPsi_turb_flux = -glsPsi(k,iCell) !Assumes k=0 in boundary condition
                        ! div(h \kappa_2 \nabla \phi) at cell center
                        flux = eddyDiff2 * layerThickEdge(k, iEdge) * glsPsi_turb_flux * dvEdge(iEdge)/dcEdge(iEdge)

                        glsPsiTend(k, iCell) = glsPsiTend(k, iCell) - flux * invAreaCell
                     end if
                  end do !nEdges loop
               end if !end cell mask if
            end do
         end do !nCells
      end if !config_two_equation_model_choice

      call mpas_timer_stop("glsPsi horizontal turbulent transport")

   !--------------------------------------------------------------------

   end subroutine ocn_glsPsi_transport_horizontal_compute

!***********************************************************************
!
!  routine ocn_glsPsi_transport_vertical_compute
!
!> \brief   Compute vertical component of parameterized glsPsi transport
!> \author  Luke Van Roekel
!> \date    July 2021
!> \details
!>  This routine initializes a variety of quantities related to
!>  parameterized glsPsi transport in the horizontal in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_glsPsi_transport_vertical_compute(glsPsi,layerThickness,sigma_psi,dt)

      real (kind=RKIND), dimension(:,:), intent(inout) :: glsPsi
      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness
      real (kind=RKIND), intent(in) :: dt
      integer :: iCell, k, N
      real (kind=RKIND), dimension(nVertLevels) :: A,B,C,dissTmp,rhs
      real (kind=RKIND),intent(in) :: sigma_psi

      call mpas_timer_start("glsPsi vertical turbulent transport")

      !$omp parallel
      !$omp do schedule(runtime) &
      !$omp private(N, A, B, C, rhs, dissTmp, k)
      do iCell=1,nCellsOwned
         N = maxLevelCell(iCell)
         A(1) = 0.0_RKIND
         do k = 2, N
            A(k) = -2.0_RKIND*dt*(dynamicViscosity(k,iCell)/rho_sw+vertViscTopOfCell(k,iCell)/sigma_psi) &
                  / (layerThickness(k-1,iCell) + layerThickness(k,iCell)) &
                  / layerThickness(k,iCell)
         end do

         do k = 1, N-1
            C(k) = -2.0_RKIND*dt*(dynamicViscosity(k+1,iCell)/rho_sw+vertViscTopOfCell(k+1,iCell)/sigma_psi) &
                  / (layerThickness(k,iCell) + layerThickness(k+1,iCell)) &
                  / layerThickness(k,iCell)
         end do
         C(N) = 0.0_RKIND

         do k = 1, N
            B(k) = 1.0_RKIND - C(k) - A(k)
         end do

         rhs(:) = glsPsi(:,iCell)
         dissTmp(:) = 0.0_RKIND

         call turb_tridiagonal_solve(A(2:N), B, C(1:N-1), rhs(:), dissTmp, &
                                N)

         glsPsi(1:N,iCell) = dissTmp(1:N)
         glsPsi(N+1:nVertLevels,iCell) = 0.0_RKIND
      end do
      !$omp end do
      !$omp end parallel

      call mpas_timer_stop("glsPsi vertical turbulent transport")

   end subroutine ocn_glsPsi_transport_vertical_compute

!***********************************************************************

   subroutine turb_tridiagonal_solve(a,b,c,r,x,n) !{{{!
!
      !-----------------------------------------------------------------!
      !!
      ! input variables!
      !!
      !-----------------------------------------------------------------!
!
      integer,intent(in) :: n!
      real (KIND=RKIND), dimension(n), intent(in) :: a,b,c,r!
!
      !-----------------------------------------------------------------!
      !!
      ! output variables!
      !!
      !-----------------------------------------------------------------!
!
      real (KIND=RKIND), dimension(n), intent(out) :: x!
!
      !-----------------------------------------------------------------!
      !!
      ! local variables!
      !!
      !-----------------------------------------------------------------!

      real (KIND=RKIND), dimension(n) :: bTemp,rTemp!
      real (KIND=RKIND) :: m!
      integer i!
!
      ! Use work variables for b and r!
      bTemp(1) = b(1)!
      rTemp(1) = r(1)!
!
      ! First pass: set the coefficients!
      do i = 2,n!
         m = a(i-1)/bTemp(i-1)!
         bTemp(i) = b(i) - m*c(i-1)!
         rTemp(i) = r(i) - m*rTemp(i-1)!
      end do!
!
      x(n) = rTemp(n)/bTemp(n)!
       ! Second pass: back-substition!
      do i = n-1, 1, -1!
         x(i) = (rTemp(i) - c(i)*x(i+1))/bTemp(i)!
      end do!
!
   end subroutine turb_tridiagonal_solve

end module ocn_tke_glsPsi_transport

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
