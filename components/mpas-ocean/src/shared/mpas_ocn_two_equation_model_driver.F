! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_two_equation_model
!
!> \brief MPAS ocean two equation model driver
!> \author Luke Van Roekel
!> \date   June 2021
!> \details
!>    This module contains all necessary routines to compute 
!>    tendencies for the two equation turbulence closure 
!>    this is primarily meant to be for nonhydrostatic mode, but could be
!>    extended to the hydrostatic code as a GOTM replacement.
!!
!-----------------------------------------------------------------------

module ocn_two_equation_model

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_kind_types
   use mpas_dmpar
   use mpas_log
   use ocn_constants
   use ocn_mesh
   use ocn_tke_advection_mono
   use ocn_glsPsi_advection_mono
   use ocn_tke_glsPsi_transport
   use ocn_shear_production
   use ocn_buoyancy_production
   use ocn_diagnostics_variables
   use ocn_two_equation_structure_functions
   use ocn_config
   use ocn_tracer_advection_shared

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
   ! Public Member Functions
   !
   !--------------------------------------------------------------------

   public :: ocn_two_equation_init, &
             ocn_compute_two_equation_tendency, &
             ocn_compute_strain, &
             ocn_compute_two_equation_dissipation, &
             ocn_compute_two_equation_diff, &
             ocn_compute_two_equation_molecular_viscosity

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real(kind=RKIND) :: &
      c_mu, c_psi1, c_psi2, c_psi3, sigma_k, sigma_psi, c_p, c_m, c_n

   logical :: twoEquationMixingOn
!***********************************************************************

contains

!***********************************************************************
!
!   routine ocn_compute_two_equation_dissipation
!>  \brief    Compute the tendency terms for the k and psi terms 
!>  \author   Luke Van Roekel
!>  \date     June 2021
!>  \details
!>            
!>            
!>            
!>            
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_two_equation_dissipation(tke, glsPsi, tkeTend, glsPsiTend)

     real(kind=RKIND), dimension(:,:), intent(in) :: tke,        & ! tke at current time level
                                                     glsPsi        ! generalized variable from 
                                                                   ! Umlauf and Burchard 2003

     real(kind=RKIND), dimension(:,:), intent(inout) :: tkeTend,    & ! Tendency of TKE
                                                        glsPsiTend    ! tendency of Psi variable

     real(kind=RKIND) :: dissipation

     integer :: iCell, k

     do iCell = 1, nCellsOwned
        do k = 1, maxLevelCell(iCell)
           dissipation = c_mu**(3.0+c_p/c_n)*tke(k,iCell)**(1.5+c_m/c_n)* &
                         glsPsi(k,iCell)**(-1.0/c_n)
           tkeTend(k,iCell) = tkeTend(k,iCell) - dissipation
           glsPsiTend(k,iCell) = glsPsiTend(k,iCell) - c_psi2*glsPsi(k,iCell) / &
                                 (tke(k,iCell) + 1.0E-15_RKIND)*dissipation
        end do
     end do

   end subroutine ocn_compute_two_equation_dissipation

!***********************************************************************
!
!   routine ocn_compute_two_equation_tendency
!>  \brief    Compute the tendency terms for the k and psi terms 
!>  \author   Luke Van Roekel
!>  \date     June 2021
!>  \details
!>            
!>            
!>            
!>            
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_two_equation_tendency(tkeTend, glsPsiTend, tke, glsPsi, velocityZonal, &
                    velocityMeridional, normalVelocity, verticalVelocity, layerThickness,        &
                    activeTracers, dt, err)

     real(kind=RKIND),dimension(:,:), intent(in) :: velocityZonal, velocityMeridional, &
                                                    normalVelocity, layerThickness,    &
                                                    verticalVelocity, tke, glsPsi

     real(kind=RKIND),dimension(:,:), intent(inout) :: tkeTend, glsPsiTend

     real (kind=RKIND), dimension(:,:,:), intent(in) :: activeTracers
     real(kind=RKIND),dimension(:,:), allocatable :: normalThicknessFlux

     integer,intent(inout) :: err

     real(kind=RKIND), intent(in) :: dt

     integer :: iCell, k, nCells, nEdges, iEdge

     nCells = nCellsHalo(2)

     do iCell=1,nCells
        do k=1,nVertLevels
           tkeTend(k,iCell) = 0.0_RKIND
           glsPsiTend(k,iCell) = 0.0_RKIND
        end do
     end do

     nEdges = nEdgesAll
     allocate(normalThicknessFlux(nVertLevels,nEdges))
     do iEdge=1,nEdges
        do k=1,nVertLevels
           normalThicknessFlux(k,iEdge) = normalVelocity(k,iEdge)*layerThickEdge(k,iEdge)
        end do
     end do

     !Compute the strain for the remaining terms in the tendency
     call ocn_compute_strain(velocityZonal, velocityMeridional,  &
                verticalVelocity, activeTracers, layerThickness, err)
     call ocn_compute_two_equation_molecular_viscosity(activeTracers)

     call ocn_tke_advection_mono_tend(tkeTend, tke, layerThickness,    &
                                      normalThicknessFlux, verticalVelocity, dt,       &
                                      advCoefs, advCoefs3rd,            &
                                      nAdvCellsForEdge, advCellsForEdge,&
                                      advMaskHighOrder)

     call ocn_tke_transport_horizontal_compute(tke, tkeTend, sigma_k, err) 

     if (.not. config_LES_mode) then

        call ocn_glsPsi_advection_mono_tend(glsPsiTend, glsPsi, layerThickness,    &
                                            normalThicknessFlux, verticalVelocity, dt,       &
                                            advCoefs, advCoefs3rd,            &
                                            nAdvCellsForEdge, advCellsForEdge,&
                                            advMaskHighOrder)

         call ocn_glsPsi_transport_horizontal_compute(glsPsi, glsPsiTend, c_psi2, sigma_psi, err)

     endif

     call ocn_compute_two_equation_dissipation(tke, glsPsi, tkeTend, glsPsiTend)
     call ocn_buoyancy_production_compute(tkeTend, glsPsiTend, glsPsi, tke, c_psi3, &
                                          surfaceBuoyancyForcing, err)
     call ocn_shear_production_compute(tkeTend, glsPsiTend, tke, glsPsi, c_psi1, &
                                       surfaceFrictionVelocity, c_mu, err)
     deallocate(normalThicknessFlux)
  end subroutine ocn_compute_two_equation_tendency

!***********************************************************************
!
!   routine ocn_compute_two_equation_diff
!>  \brief    Computes the diffusivity for the two equation model
!>  \author   Luke Van Roekel
!>  \date     May 2021
!>  \details
!>            Computes the diffusivity/viscosity for the two equation model
!>            this is based on the tke and dissipation and diagnosed length 
!>            scale
!>
!>            The diffusivity and viscosity are placed in the vertDiffTopOfCell
!>            and vertViscTopOfCell arrays respectively
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_two_equation_diff(tke, glsPsi, layerThickness, dt, err)

     real(kind=RKIND),dimension(:,:), intent(inout) :: glsPsi, tke
     
     real(kind=RKIND),dimension(:,:), intent(in) :: layerThickness

     integer,intent(inout) :: err

     real(kind=RKIND), intent(in) :: dt

     real(kind=RKIND) :: length, dissipation, length_grid, length_stable, dcAv
     integer :: j, iCell, k, nCells, nEdges, iEdge, spot(2)

     !compute new diffusivity based on older values
     call ocn_canuto_structure_function_compute(Sm, Sh, tke, glsPsi, c_p, c_m, c_n, c_mu)
     nCells = nCellsHalo(1)

     if (config_LES_mode) then
        !$omp parallel
        !$omp do schedule(runtime) private(k,iEdge,dcAv,length_grid,length_stable,length)
        do iCell=1,nCells
           dcAv = 0.0_RKIND
           do j=1,nEdgesOnCell(iCell)
              iEdge = edgesOnCell(j, iCell)
              dcAv = dcAv + dcEdge(iEdge)
           end do
           dcAv = dcAv / (1.0E-15_RKIND + nEdgesOnCell(iCell))
           do k=2,maxLevelCell(iCell)
              !compute length scale
              !define the viscosity
              length_grid = sqrt(dcAv * 0.5_RKIND*(layerThickness(k-1,iCell) + layerThickness(k,iCell)))
              length_stable = 0.76_RKIND*sqrt(tke(k,iCell))/ &
                              (sqrt(max(BruntVaisalaFreqTop(k,iCell),0.0_RKIND)) + &
                                                             1.0E-15_RKIND)
              length = min(length_grid, length_stable)
              vertViscTopOfCell(k,iCell) = smagorinsky_coefficient*length*sqrt(tke(k,iCell))
              vertDiffTopOfCell(k,iCell) = vertViscTopOfCell(k,iCell) !Assume turbulent Pr of 1
           end do
        end do
        !$omp end do
        !$omp end parallel
     else
        !$omp parallel
        !$omp do schedule(runtime) private(k,length,dissipation)
        do iCell=1,nCells
           vertViscTopOfCell(1,iCell) = 0.0_RKIND
           vertViscTopOfCell(maxLevelCell(iCell)+1:nVertLevels,iCell) = 0.0_RKIND
           do k=2,maxLevelCell(iCell)
              dissipation = c_mu**(3.0+c_p/c_n)*tke(k,iCell)**(1.5+c_m/c_n)* &
                            glsPsi(k,iCell)**(-1.0/c_n)
              length = c_mu**3*tke(k,iCell)**1.5/(1.0E-15_RKIND + dissipation)
              vertViscTopOfCell(k,iCell) = length*Sm(k,iCell)*sqrt(tke(k,iCell))
           end do
        end do
        !$omp end do 
        !$omp end parallel
     end if

     call ocn_tke_transport_vertical_compute(tke,layerThickness,sigma_k,dt)
     if(.not. config_LES_mode) then
        call ocn_glsPsi_transport_vertical_compute(glsPsi,layerThickness,sigma_psi,dt)
     end if

     if (config_LES_mode) then
       if(minval(tke) < 0) then
         spot = minloc(tke)
         print *, 'tke = ', tke(:,spot(2))
         print *, 'vdiff = ',vertViscTopOfCell(:,spot(2))
         print *, 'dyn = ',dynamicViscosity(:,spot(2))
         print *, 'minval = ',minval(tke)
         print *, 'minloc = ',minloc(tke)
       endif
        !$omp parallel
        !$omp do schedule(runtime) private(k,dcAv,iEdge,length_grid,length_stable,length)
        do iCell=1,nCells
           dcAv = 0.0_RKIND
           do j=1,nEdgesOnCell(iCell)
              iEdge = edgesOnCell(j, iCell)
              dcAv = dcAv + dcEdge(iEdge)
           end do
           dcAv = dcAv / (1.0E-15_RKIND + nEdgesOnCell(iCell))
           do k=2,maxLevelCell(iCell)
              !compute length scale
              !define the viscosity
              length_grid = sqrt(dcAv * 0.5_RKIND*(layerThickness(k-1,iCell) + layerThickness(k,iCell)))
              length_stable = 0.76_RKIND*sqrt(max(tke(k,iCell),0.0_RKIND))/(sqrt(max(BruntVaisalaFreqTop(k,iCell),0.0_RKIND)) + &
                                                             1.0E-15_RKIND)
              length = min(length_grid, length_stable)
              vertViscTopOfCell(k,iCell) = smagorinsky_coefficient*length*sqrt(tke(k,iCell)) + &
                              dynamicViscosity(k,iCell)/rho_sw
              vertDiffTopOfCell(k,iCell) = vertViscTopOfCell(k,iCell) !Assume turbulent Pr of 1
           end do
        end do
        !$omp end do
        !$omp end parallel
     else
       ! Recompute now that everything is fully up to date
        call ocn_canuto_structure_function_compute(Sm, Sh, tke, glsPsi, c_p, c_m, c_n, c_mu)
        !$omp parallel
        !$omp do schedule(runtime) private(k,dissipation,length)
        do iCell=1,nCells
           vertDiffTopOfCell(1,iCell) = 0.0_RKIND
           vertViscTopOfCell(1,iCell) = 0.0_RKIND
           vertDiffTopOfCell(maxLevelCell(iCell)+1:nVertLevels,iCell) = 0.0_RKIND
           vertViscTopOfCell(maxLevelCell(iCell)+1:nVertLevels,iCell) = 0.0_RKIND
           do k=2,maxLevelCell(iCell)
              dissipation = c_mu**(3.0+c_p/c_n)*tke(k,iCell)**(1.5+c_m/c_n)* &
                            glsPsi(k,iCell)**(-1.0/c_n)
              length = c_mu**3*tke(k,iCell)**1.5/(1.0E-15_RKIND + dissipation)
              vertViscTopOfCell(k,iCell) = length*Sm(k,iCell)*sqrt(tke(k,iCell))
              vertDiffTopOfCell(k,iCell) = length*Sh(k,iCell)*sqrt(tke(k,iCell))
           end do
        end do
        !$omp end do 
        !$omp end parallel
     end if

   end subroutine ocn_compute_two_equation_diff

!***********************************************************************
!
!   routine ocn_compute_strain
!>  \brief    Computes the strain rate tensor magnitude for smagorinsky 
!>  \author   Luke Van Roekel
!>  \date     May 2021
!>  \details
!>            Computes the strain rate tensor at cell centers.  Uses
!>            the vector reconstruct coefficients to get x and y derivs
!>            For 2D need dudx, dudy, dvdx, dvdy
!>            For 3D adding dwdx and dwdy
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_strain(velocityZonal, velocityMeridional,  &
                verticalVelocity, activeTracers, layerThickness, err)

     real(kind=RKIND),dimension(:,:),intent(in) :: &
         velocityZonal, &  !Input: Zonal Velocity for gradients
         velocityMeridional, & !Input: Meridional Velocity for gradients
         layerThickness

     real(kind=RKIND),dimension(:,:),intent(in),optional :: &
         verticalVelocity !Input: optionally for nonhydrostaticMode include w

     real(kind=RKIND),dimension(:,:,:),intent(in) :: &
         activeTracers

     integer, intent(in) :: err

     integer :: iCell, iEdge, j, k, cell1, cell2

     ! initialize gradients

     do iCell=1,nCellsOwned
        ugrad(:,:,iCell) = 0.0_RKIND
        vgrad(:,:,iCell) = 0.0_RKIND
        tgrad(:,:,iCell) = 0.0_RKIND

        do j = 1,nEdgesOnCell(iCell) !for device may be better to move next lines inside k loop
           iEdge = edgesOnCell(j,iCell)
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
           do k = 1,maxLevelEdgeTop(iEdge)
              tgrad(1,k,iCell) = tgrad(1,k,iCell) + coeffs_reconstruct(1,j,iCell) * edgeSignOnCell(j,iCell) * &
                            (activeTracers(1,k,cell2) - activeTracers(1,k,cell1)) / dcEdge(iEdge) *   &
                            edgeMask(k,iEdge)
              tgrad(2,k,iCell) = tgrad(2,k,iCell) + coeffs_reconstruct(2,j,iCell) * edgeSignOnCell(j,iCell) * &
                            (activeTracers(1,k,cell2) - activeTracers(1,k,cell1)) / dcEdge(iEdge) *   &
                            edgeMask(k,iEdge)
              sgrad(1,k,iCell) = sgrad(1,k,iCell) + coeffs_reconstruct(1,j,iCell) * edgeSignOnCell(j,iCell) * &
                            (activeTracers(2,k,cell2) - activeTracers(2,k,cell1)) / dcEdge(iEdge) *   &
                            edgeMask(k,iEdge)
              sgrad(2,k,iCell) = sgrad(2,k,iCell) + coeffs_reconstruct(2,j,iCell) * edgeSignOnCell(j,iCell) * &
                            (activeTracers(2,k,cell2) - activeTracers(2,k,cell1)) / dcEdge(iEdge) *   &
                            edgeMask(k,iEdge)
              ugrad(1,k,iCell) = ugrad(1,k,iCell) + coeffs_reconstruct(1,j,iCell) * edgeSignOnCell(j,iCell) * &
                            (velocityZonal(k,cell2) - velocityZonal(k,cell1)) / dcEdge(iEdge) *   &
                            edgeMask(k,iEdge)
              ugrad(2,k,iCell) = ugrad(2,k,iCell) + coeffs_reconstruct(2,j,iCell) * edgeSignOnCell(j,iCell) * &
                            (velocityZonal(k,cell2) - velocityZonal(k,cell1)) / dcEdge(iEdge) *   &
                            edgeMask(k,iEdge)
              vgrad(1,k,iCell) = vgrad(1,k,iCell) + coeffs_reconstruct(1,j,iCell) * edgeSignOnCell(j,iCell) * &
                            (velocityMeridional(k,cell2) - velocityMeridional(k,cell1)) /         &
                            dcEdge(iEdge) * edgeMask(k,iEdge)
              vgrad(2,k,iCell) = vgrad(2,k,iCell) + coeffs_reconstruct(2,j,iCell) * edgeSignOnCell(j,iCell) * &
                            (velocityMeridional(k,cell2) - velocityMeridional(k,cell1)) /         &
                            dcEdge(iEdge) * edgeMask(k,iEdge)
          end do 
       end do !end n edges loop

       ! instead of doing this put in surface BCs for generation
       ! compute dudz dvdz
       ugrad(3,1,iCell) = (velocityZonal(1,iCell) - velocityZonal(2,iCell)) /   &
                            (0.5_RKIND*(layerThickness(1,iCell) + layerThickness(2,iCell)))
       vgrad(3,1,iCell) = (velocityMeridional(1,iCell) - velocityMeridional(2,iCell)) /   &
                            (0.5_RKIND*(layerThickness(1,iCell) + layerThickness(2,iCell)))
       tgrad(3,1,iCell) = (activeTracers(1,1,iCell) - activeTracers(1,2,iCell)) /   &
                            (0.5_RKIND*(layerThickness(1,iCell) + layerThickness(2,iCell)))
       sgrad(3,1,iCell) = (activeTracers(2,1,iCell) - activeTracers(2,2,iCell)) /   &
                            (0.5_RKIND*(layerThickness(1,iCell) + layerThickness(2,iCell)))

       do k=2,maxLevelCell(iCell)
          tgrad(3,k,iCell) = 0.5_RKIND*(tgrad(3,k-1,iCell) + (activeTracers(1,k-1,iCell) - activeTracers(1,k,iCell)) / &
                        (0.5_RKIND*(layerThickness(k-1,iCell) + layerThickness(k,iCell))))
          sgrad(3,k,iCell) = 0.5_RKIND*(sgrad(3,k-1,iCell) + (activeTracers(2,k-1,iCell) - activeTracers(2,k,iCell)) / &
                        (0.5_RKIND*(layerThickness(k-1,iCell) + layerThickness(k,iCell))))
          ugrad(3,k,iCell) = 0.5_RKIND*(ugrad(3,k-1,iCell) + (velocityZonal(k-1,iCell) - velocityZonal(k,iCell)) / &
                        (0.5_RKIND*(layerThickness(k-1,iCell) + layerThickness(k,iCell))))
          vgrad(3,k,iCell) = 0.5_RKIND*(vgrad(3,k-1,iCell) + (velocityMeridional(k-1,iCell) - velocityMeridional(k,iCell)) / &
                        (0.5_RKIND*(layerThickness(k-1,iCell) + layerThickness(k,iCell))))
       end do

    end do ! end cell loop

    ! add w part if requested 
    if(present(verticalVelocity)) then
       do iCell=1,nCellsOwned
          wgrad(:,:,iCell) = 0.0_RKIND

          do j = 1,nEdgesOnCell(iCell) !for device may be better to move next lines inside k loop
             iEdge = edgesOnCell(j,iCell)
             cell1 = cellsOnEdge(1,iEdge)
             cell2 = cellsOnEdge(2,iEdge)
         
             do k=1,maxLevelEdgeTop(iEdge)
                wgrad(1,k,iCell) = wgrad(1,k,iCell) + coeffs_reconstruct(1,j,iCell) * edgeSignOnCell(j,iCell) * &
                             (verticalVelocity(k,cell2) - verticalVelocity(k,cell1)) / dcEdge(iEdge) * &
                             edgeMask(k,iEdge)
                wgrad(2,k,iCell) = wgrad(2,k,iCell) + coeffs_reconstruct(2,j,iCell) * edgeSignOnCell(j,iCell) * &
                             (verticalVelocity(k,cell2) - verticalVelocity(k,cell1)) / dcEdge(iEdge) * &
                             edgeMask(k,iEdge)
             end do

          end do !end n edges loop

          do k=2,maxLevelCell(iCell)
             wgrad(3,k,iCell) = (verticalVelocity(k-1,iCell) - verticalVelocity(k,iCell)) / &
                         ((layerThickness(k,iCell)))
          end do

       end do !end icell loop
    end if !if present w 

 end subroutine ocn_compute_strain

!***********************************************************************
!
!   routine ocn_two_equation_init
!>  \brief    Configures options for the two equation model 
!>  \author   Luke Van Roekel
!>  \date     May 2021
!>  \details
!>            Initializes variables and allocates internal variables for 
!>            computing strain tensor and k and psi 
!
!----------------------------------------------------------------------

   subroutine ocn_two_equation_init( domain, err )

      type (domain_type), intent(inout) :: domain
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      integer, intent(out) :: err
      real(kind=RKIND), dimension(:,:), pointer :: layerThickness
      real(kind=RKIND), dimension(:,:), pointer :: distanceToBoundary
      real(kind=RKIND), dimension(:,:), pointer :: boundaryNormalAngle

      integer :: cell1, cell2, iEdge, k

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      twoEquationMixingOn = .false.
      if(config_use_two_equation_turbulence_model) then

         twoEquationMixingOn = .true.

         if (config_two_equation_model_choice == 'kl') then !MY 1982
            c_mu = 0.5544_RKIND
            sigma_k = 1.96_RKIND
            sigma_psi = 1.96_RKIND
            c_psi1 = 0.9_RKIND
            c_psi2 = 0.5_RKIND
            c_psi3 = 0.9_RKIND
            c_p = 0.0_RKIND
            c_m = 1.0_RKIND
            c_n = 1.0_RKIND
         elseif (config_two_equation_model_choice == 'epsilon') then ! Rotta 87
            c_mu = 0.5477_RKIND
            sigma_k = 1.0_RKIND
            sigma_psi = 1.3_RKIND
            c_psi1 = 1.44_RKIND
            c_psi2 = 1.92_RKIND
            c_psi3 = -0.629_RKIND
            c_p = 3.0_RKIND
            c_m = 3.0_RKIND / 2.0_RKIND
            c_n = -1.0_RKIND
        elseif (config_two_equation_model_choice == 'omega') then !Wilcox 88
            c_mu = 0.5477_RKIND
            sigma_k = 2.0_RKIND
            sigma_psi = 2.0_RKIND
            c_psi1 = 0.555_RKIND
            c_psi2 = 0.833_RKIND
            c_psi3 = -0.642_RKIND
            c_p = -1.0_RKIND
            c_m = 0.5_RKIND
            c_n = -1.0_RKIND
         elseif (config_two_equation_model_choice == 'ktau') then !Zilinkivech Wilcox 88
            c_mu = 0.5477_RKIND
            sigma_k = 1.46_RKIND
            sigma_psi = 10.8_RKIND
            c_psi1 = 0.173_RKIND
            c_psi2 = 0.225_RKIND
            c_psi3 = 0.0_RKIND
            c_p = -3.0_RKIND
            c_m = 0.5_RKIND
            c_n = 1.0_RKIND
         elseif (config_two_equation_model_choice == 'user-defined') then
            c_mu = config_user_c_mu
            sigma_k = config_user_sigma_k
            sigma_psi = config_user_sigma_psi
            c_psi1 = config_user_c_psi1
            c_psi2 = config_user_c_psi2
            c_psi3 = config_user_c_psi3
            c_p = config_user_c_p
            c_m = config_user_c_m
            c_n = config_user_c_n
         else !invalid choice 
            call mpas_log_write("Unknown option for config_two_equation_model_choice, valid values are: " // &
                                  " 'epsilon', 'kl', 'omega', 'ktau', 'user-defined'", MPAS_LOG_CRIT)
         endif

         if(config_use_cvmix) then
            call mpas_log_write("Cannot use two equation turbulence closure and cvmix at the same time: " // &
                    "set either config_use_cvmix or config_use_two_equation_model to .false.", MPAS_LOG_CRIT)
         end if

         !init advection of tke and glsPsi

         call ocn_tke_advection_mono_init(config_num_halos, config_two_equation_hor_adv_order,        &
                                                config_two_equation_vert_adv_order,              &
                                                config_two_equation_coef_3rd_order,              &
                                                config_two_equation_check_monotonicity, err)

         call ocn_glsPsi_advection_mono_init(config_num_halos, config_two_equation_hor_adv_order,        &
                                                config_two_equation_vert_adv_order,                 &
                                                config_two_equation_coef_3rd_order,                 &
                                                config_two_equation_check_monotonicity, err)

         call mpas_pool_get_array(meshPool, 'distanceToBoundary', distanceToBoundary)
         call mpas_pool_get_array(meshPool, 'boundaryNormalAngle', boundaryNormalAngle)
         if(.not. associated(distanceToBoundary)) then
            distanceToBoundary(:,:) = 1.0E50_RKIND
            ! make an assumption that the distance to boundary is just dcEdge/2.0 when edgeMask==0
            do iEdge=1,nEdgesAll
               do k=1,nVertLevels
                  if(edgeMask(k,iEdge) == 0) then
                     cell1 = cellsOnEdge(1,iEdge)
                     cell2 = cellsOnedge(2,iEdge)
                     distanceToBoundary(k,cell1) = dcEdge(iEdge)/2.0_RKIND
                     distanceToBoundary(k,cell2) = dcEdge(iEdge)/2.0_RKIND
                  endif
               end do
            end do
         end if

         if(.not. associated(boundaryNormalAngle)) then
            if(config_two_equation_model_choice == 'epsilon') then
               call mpas_log_write("Unable to use k-epsilon configuration without specifying the " // &
                        "boundaryNormalAngle array", MPAS_LOG_CRIT)
            end if
         end if
      end if

   end subroutine ocn_two_equation_init

   subroutine ocn_compute_two_equation_molecular_viscosity(activeTracers)

     !computes the molecular viscosity of water following Sharqawy et al (2010)
     real(kind=RKIND),dimension(:,:,:), intent(in) :: activeTracers
     real(kind=RKIND) :: A, B, C, dcEdgeAverage
     real(kind=RKIND) :: tAvg, sAvg, viscAvg, A1, B1, freshwaterMu

     real(kind=RKIND) :: boundaryTangentialVel
     integer :: iCell, j, k

     A = 2.414E-5_RKIND
     B = 247.8_RKIND
     C = 140.0_RKIND

     do iCell=1,nCellsAll
        !compute dcEdgeAverage
        dcEdgeAverage = 0.0_RKIND
        do j=1,nEdgesOnCell(iCell)
           dcEdgeAverage = dcEdgeAverage + dcEdge(edgesOnCell(j,iCell))
        end do
        dcEdgeAverage = dcEdgeAverage / nEdgesOnCell(iCell) + 1.0E-15_RKIND
        dynamicViscosity(1,iCell) = 0.0_RKIND
        do k = 1,maxLevelCell(iCell)-1
           tAvg = 0.5_RKIND*(activeTracers(1,k,iCell) + activeTracers(1,k+1,iCell))
           sAvg = 0.5_RKIND*(activeTracers(2,k,iCell) + activeTracers(2,k+1,iCell)) / 1000.0_RKIND
           viscAvg = 1.0E-15_RKIND + 0.5_RKIND*(dynamicViscosity(k,iCell) + dynamicViscosity(k+1,iCell))
!           freshwaterMu = A*10.0_RKIND**(B/(tAvg + 273.15_RKIND - C))
           freshwaterMu = 4.2844E-5 + (0.157_RKIND*(tAvg + 64.993_RKIND)**2.0 - 91.296_RKIND)**(-1.0)
           A1 = 1.541_RKIND + 0.01998*tAvg - 9.52E-5_RKIND*tAvg**2.0
           B1 = 7.974_RKIND - 0.07561*tAvg + 4.724E-4_RKIND*tAvg**2.0
           dynamicViscosity(k+1,iCell) = freshwaterMu*(1.0_RKIND + &
                                        A1*sAvg + &
                                        B1*sAvg**2)
          if(dynamicViscosity(k+1,iCell) < 0) then
            print *, activeTracers(1,:,iCell)
            print *, 'hi',k, tavg,savg, viscavg, freshwaterMu, A1, B1, dynamicViscosity(k+1,iCell)
            stop
          endif
           boundaryTangentialVel = velocityZonal(k,iCell)*sin(boundaryNormalAngle(k,iCell)) + &
                                  velocityMeridional(k,iCell)*cos(boundaryNormalAngle(k,iCell))
           !need to find the edge in question and then use that dcEdge -- maybe just do average for performance reasons?

           uStar(k,iCell) = sqrt(boundaryCellMask(k,iCell)*abs(boundaryTangentialVel/(0.5_RKIND*dcEdgeAverage)) &
                          * viscAvg/rho_sw)
           yPlus(k,iCell) = uStar(k,iCell)*distanceToBoundary(k,iCell) / (viscAvg/rho_sw)
        end do
        k=maxLevelCell(iCell)
        dynamicViscosity(maxLevelCell(iCell)+1:nVertLevels+1,iCell) = 0.0_RKIND

        viscAvg = 0.5_RKIND*(dynamicViscosity(k,iCell) + dynamicViscosity(k+1,iCell))
        boundaryTangentialVel = velocityZonal(k,iCell)*sin(boundaryNormalAngle(k,iCell)) + &
                                  velocityMeridional(k,iCell)*cos(boundaryNormalAngle(k,iCell))
           !need to find the edge in question and then use that dcEdge -- maybe just do average for performance reasons?

        uStar(k,iCell) = sqrt(boundaryCellMask(k,iCell)*abs(boundaryTangentialVel/(0.5_RKIND*dcEdgeAverage)) &
                          * viscAvg/rho_sw)
        yPlus(k,iCell) = uStar(k,iCell)*distanceToBoundary(k,iCell) / (viscAvg/rho_sw)

     end do

     !compute ustar

   end subroutine ocn_compute_two_equation_molecular_viscosity
end module ocn_two_equation_model

