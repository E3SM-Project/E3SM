! Copyright (c) 2021,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_fblts
!
!> \brief  MPAS barotropic ocean LTS Time integration scheme
!> \author Jeremy Lilly
!> \date   October 2023
!> \details
!>  This module contains the FB_LTS init routine and the FB_LTS
!>  barotropic ocean time integration scheme with splitting
!>  on the fast and slow tendency terms.
!
!-----------------------------------------------------------------------

module ocn_time_integration_fblts
       
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_timer

   use ocn_tendency
   use ocn_diagnostics
   use ocn_mesh
   use ocn_vmix
   use ocn_config
   use ocn_time_average_coupled

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------
   public :: ocn_time_integrator_fblts, &
             ocn_time_integration_fblts_init

   contains


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integrator_fblts
!
!> \brief  MPAS barotropic ocean FB_LTS time integration scheme
!> \author Jeremy Lilly
!> \date   October 2023
!> \details
!>  This routine integrates one timestep (dt) using an FB_LTS time
!>  integrator with a splitting of the fast and slow tendency terms
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integrator_fblts(domain, dt)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advance model state forward in time by the specified time step 
   ! using a local time stepping scheme with splitting of the fast and 
   ! slow tendencies
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
         dt  !< [in] time step (sec) to move forward

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain  !< [inout] model state to advance forward

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer :: &
         iCell, iEdge, iRegion, k, ic, ie, im, &  ! iterators
         M, &  ! M = dtCoarse / dtFine
         nRegions  ! number of interface regions (two)

      type (block_type), pointer :: &
         block  ! structure with subdomain data

      type (mpas_pool_type), pointer :: &
         tendPool, &  ! structure holding tendencies
         statePool, &  ! structure holding state variables
         meshPool, &  ! structure holding mesh variables
         verticalMeshPool, &  ! structure holding mesh variables
         forcingPool, &  ! structure holding forcing variables
         scratchPool, &  ! structure holding temporary variables
         tracersPool  ! structure holding tracers variables

      ! LTS Pools
      type (mpas_pool_type), pointer :: &
         LTSPool, &  ! structure holding LTS variables
         tendSlowPool, &  ! structure holding the slow tendency variables
         tendSum3rdPool, &  ! structure holding one of the correction terms for the interface
         prevTendSlowPool, nextTendSlowPool, &  ! structures containing intermediate data
         prevTendSum3rdPool, nextTendSum3rdPool  ! structures containing intermediate data

      ! Tend Array Pointers
      real (kind=RKIND), dimension(:,:), pointer ::  &
         normalVelocityTend, &  ! normal velocity fast tendency
         layerThicknessTend, &  ! layer thickness tendency
         normalVelocityTendSlow, &  ! normal velocity slow tendency
         normalVelocityTendSum3rd, &  ! one of the normal velocity correction terms for the interface
         layerThicknessTendSum3rd   ! one of the layer thickness correction terms for the interface

      ! State Array Pointers
      real (kind=RKIND), dimension(:,:), pointer :: &
         normalVelocityCur, &  ! normal velocity at time n
         normalVelocityNew, &  ! normal velocity at time n+1
         normalVelocityFirstStage, &  ! normal velocity at first stage of LTS
         normalVelocitySecondStage, &  ! normal velocity at second stage of LTS
         normalVelocityForTend, &  ! extra variable to store averages of data for tends calculations
         layerThicknessCur, &  ! layer thickness at time n
         layerThicknessNew, &  ! layer thickness at time n+1
         layerThicknessFirstStage, &  ! layer thickness at first stage of LTS
         layerThicknessSecondStage, &  ! layer thickness at second stage of LTS
         layerThicknessForTend  ! extra variable to store averages of data for tends calculations
      
      ! Local pointer for ssh
      real (kind=RKIND), dimension(:), pointer :: &
         ssh

      ! LTS objects
      real (kind=RKIND) :: &
         dtFine  ! fine dt, defined as dt / M
      integer, dimension(:,:), pointer :: &
         nCellsInLTSRegion, &  ! number of cells in a given LTS region
         nEdgesInLTSRegion  ! number of edges in a given LTS region
      integer, dimension(:,:,:), pointer :: &
         cellsInLTSRegion, &  ! list of cells in a given LTS region
         edgesInLTSRegion  ! list of edges in a given LTS region
      real (kind=RKIND) :: &
         weight1st, weight2nd  ! coefficients for each RK stage
      real (kind=RKIND) :: &
         weightTendSum3rd  ! coefficients for the interface correction

      integer err 
      err = 0
      
      !-----------------------------------------------------------------
      ! Begin routine
      !-----------------------------------------------------------------
      
      if (.not. config_disable_tr_all_tend) then
         call mpas_log_write("ERROR: tracers are not currently implemented for local time-stepping")
         call mpas_log_write("config_disable_tr_all_tend should be true in the namelist file")
         call abort
      end if
      
      
      call mpas_timer_start("FB_LTS time-step prep")
      
      ! LTS parameters 
      M = config_dt_scaling_LTS
      nRegions = 2
      dtFine = dt / M
      
      ! Weights for RK stages, weight for third stage is 1
      weight1st = 1.0_RKIND / 3.0_RKIND
      weight2nd = 1.0_RKIND / 2.0_RKIND
      
      block => domain % blocklist

      ! Retrieve model state, pools
      call mpas_pool_get_subpool(block%structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block%structs, 'verticalMesh', verticalMeshPool)
      call mpas_pool_get_subpool(block%structs, 'state', statePool)
      call mpas_pool_get_subpool(block%structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(block%structs, 'LTS', LTSPool)
      call mpas_pool_get_subpool(block%structs, 'tend', tendPool)
      call mpas_pool_get_subpool(block%structs, 'scratch', scratchPool)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      
      ! Retrieve state variables at necessary time levels
      call mpas_pool_get_array(statePool, 'normalVelocity', &
                                           normalVelocityCur, 1)
      call mpas_pool_get_array(statePool, 'normalVelocity', &
                                           normalVelocityNew, 2)
      call mpas_pool_get_array(statePool, 'normalVelocity', &
                                           normalVelocityFirstStage, 3)
      call mpas_pool_get_array(statePool, 'normalVelocity', &
                                           normalVelocitySecondStage, 4)
      call mpas_pool_get_array(statePool, 'normalVelocity', &
                                           normalVelocityForTend, 5)
      call mpas_pool_get_array(statePool, 'layerThickness', &
                                           layerThicknessCur, 1)
      call mpas_pool_get_array(statePool, 'layerThickness', &
                                           layerThicknessNew, 2)
      call mpas_pool_get_array(statePool, 'layerThickness', &
                                           layerThicknessFirstStage, 3)
      call mpas_pool_get_array(statePool, 'layerThickness', &
                                           layerThicknessSecondStage, 4)
      call mpas_pool_get_array(statePool, 'layerThickness', &
                                           layerThicknessForTend, 5)
      
      ! Retrieve tendency variables
      call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
      call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

      ! Retrieve LTS arrays
      call mpas_pool_get_array(LTSPool, 'cellsInLTSRegion', cellsInLTSRegion)
      call mpas_pool_get_array(LTSPool, 'nCellsInLTSRegion', nCellsInLTSRegion)
      call mpas_pool_get_array(LTSPool, 'edgesInLTSRegion', edgesInLTSRegion)
      call mpas_pool_get_array(LTSPool, 'nEdgesInLTSRegion', nEdgesInLTSRegion)
      
      ! Create and retrieve additional pools for LTS
      call mpas_pool_get_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)
      call mpas_pool_get_subpool(block % structs, 'tend_slow', tendSlowPool)

      call mpas_pool_copy_pool(tendPool, tendSum3rdPool, 1)
      call mpas_pool_copy_pool(tendPool, tendSlowPool, 1)


      call mpas_pool_get_array(tendSlowPool, 'normalVelocity', &
                                              normalVelocityTendSlow)
      call mpas_pool_get_array(tendSum3rdPool, 'normalVelocity', &
                                                normalVelocityTendSum3rd)
      call mpas_pool_get_array(tendSum3rdPool, 'layerThickness', &
                                                layerThicknessTendSum3rd)
      
      ! Init variables
      do iEdge = 1, nEdgesAll  ! do we need this?
         do k = 1, maxLevelEdgeTop(iEdge)
            normalVelocityNew(k, iEdge) = normalVelocityCur(k, iEdge)
            normalVelocityFirstStage(k, iEdge) = normalVelocityCur(k, iEdge)
            normalVelocitySecondStage(k, iEdge) = normalVelocityCur(k, iEdge)
         end do
      end do

      do iCell = 1, nCellsAll  ! do we need this?
         do k = 1, maxLevelCell(iCell)
            layerThicknessNew(k, iCell) = layerThicknessCur(k, iCell)
            layerThicknessFirstStage(k, iCell) = layerThicknessCur(k, iCell)
            layerThicknessSecondStage(k, iCell) = layerThicknessCur(k, iCell)
         end do
      end do

      normalVelocityTendSum3rd(:,:) = 0.0_RKIND
      layerThicknessTendSum3rd(:,:) = 0.0_RKIND
      

      call mpas_timer_stop("FB_LTS time-step prep")
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! BEGIN FB_LTS SCHEME
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      call mpas_timer_start("FB_LTS main loop")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN: Slow tendency calculation
      ! Calculate the slow tendencies on all LTS regions (only the 
      ! momentum equation contains slow terms)
      call mpas_timer_start("FB_LTS compute slow tendencies")
      
      call ocn_tend_vel(domain, tendSlowPool, statePool, forcingPool, 1, &
                        domain % dminfo, dt)

      call mpas_timer_stop("FB_LTS compute slow tendencies")
      ! END: Slow tendency calculation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN COARSE ADVANCEMENT
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN: Thickness stage 1
      ! Compute the first stage of the thickness on fine*, interface 1, 
      ! interface 2, and coarse
            
      ! Compute fast tendencies for thickness
      call mpas_timer_start("FB_LTS compute fast tendencies")
      call mpas_timer_start("FB_LTS compute fast thickness tendencies")

      call ocn_lts_thick_tend(LTSPool, tendPool, &
                              normalVelocityCur, layerThicknessCur, &
                              0, 1, 1, 1, 1)

      call mpas_timer_stop("FB_LTS compute fast thickness tendencies")
      call mpas_timer_stop("FB_LTS compute fast tendencies")

      ! Advance thickness solution with first stage
      call mpas_timer_start("FB_LTS advance solution")
      call mpas_timer_start("FB_LTS advance thickness solution")
      
      ! Interface layers
      do iRegion = 1,nRegions
         do ic = 1, nCellsInLTSRegion(iRegion,2)
            iCell = cellsInLTSRegion(iRegion,2,ic)
            layerThicknessFirstStage(:,iCell) = layerThicknessCur(:,iCell) &
                                                + weight1st * dt * layerThicknessTend(:,iCell)
         end do
      end do
      ! Coarse
      do ic = 1, nCellsInLTSRegion(2,1)
         iCell = cellsInLTSRegion(2,1,ic)
         layerThicknessFirstStage(:,iCell) = layerThicknessCur(:,iCell) &
                                             + weight1st * dt * layerThicknessTend(:,iCell)
      end do
      ! Fine layers close to interface layers
      do ic = 1, nCellsInLTSRegion(1,3)
         iCell = cellsInLTSRegion(1,3,ic)
         layerThicknessFirstStage(:,iCell) = layerThicknessCur(:,iCell) &
                                             + weight1st * dt * layerThicknessTend(:,iCell)
      end do

      call mpas_timer_stop("FB_LTS advance thickness solution")
      call mpas_timer_stop("FB_LTS advance solution")

      ! END: Thickness stage 1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN: Normal velocity stage 1
      ! Compute the first stage of the normal velocity on fine*, 
      ! interface 1, interface 2, and coarse
      
      ! Perform forward-backward average of thickness for stage one of 
      ! FB-RK(3,2):
      ! h = weight1*stage1 + (1-weight1)*current
      layerThicknessForTend(:,:) = config_fb_weight_1 * layerThicknessFirstStage(:,:) &
                                    + (1.0_RKIND - config_fb_weight_1) * layerThicknessCur(:,:)
      
      ! Halo update before tend calculation
      ! Might be able to remove this?
      call mpas_timer_start("FB_LTS prognostic halo update")

      ! Don't need to update normalVelocityCur, already updated

      call mpas_timer_start("FB_LTS thickness halo update")
      call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=5)
      call mpas_timer_stop("FB_LTS thickness halo update")

      call mpas_timer_stop("FB_LTS prognostic halo update")

      ! Compute fast tendencies for normal velocity
      call mpas_timer_start("FB_LTS compute fast tendencies")
      call mpas_timer_start("FB_LTS compute fast velocity tendencies")

      call ocn_lts_vel_tend(LTSPool, tendPool, &
                            normalVelocityCur, layerThicknessForTend, &
                            0, 1, 1, 1, 1)

      normalVelocityTend(:,:) = normalVelocityTend(:,:) &
                                + normalVelocityTendSlow(:,:)

      call mpas_timer_stop("FB_LTS compute fast velocity tendencies")
      call mpas_timer_stop("FB_LTS compute fast tendencies")

      ! Advance normal velocity solution with first stage
      call mpas_timer_start("FB_LTS advance solution")
      call mpas_timer_start("FB_LTS advance velocity solution")

      ! Interface regions
      do iRegion = 1,nRegions
         do ie = 1, nEdgesInLTSRegion(iRegion,2)
            iEdge = edgesInLTSRegion(iRegion,2,ie)
            normalVelocityFirstStage(:,iEdge) = normalVelocityCur(:,iEdge) &
                                                + weight1st * dt * normalVelocityTend(:,iEdge)
         end do
      end do
      ! Coarse
      do ie = 1, nEdgesInLTSRegion(2,1)
         iEdge = edgesInLTSRegion(2,1,ie)
         normalVelocityFirstStage(:,iEdge) = normalVelocityCur(:,iEdge) &
                                              + weight1st * dt * normalVelocityTend(:,iEdge)
      end do
      ! Fine layers close to interface layers
      do ie = 1, nEdgesInLTSRegion(1,3)
         iEdge = edgesInLTSRegion(1,3,ie)
         normalVelocityFirstStage(:,iEdge) = normalVelocityCur(:,iEdge) &
                                             + weight1st * dt * normalVelocityTend(:,iEdge)
      end do
      
      call mpas_timer_stop("FB_LTS advance velocity solution")
      call mpas_timer_stop("FB_LTS advance solution")

      ! END: Normal velocity stage 1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN: Thickness stage 2
      ! Compute the second stage of the thickness on fine*, interface 1, 
      ! interface 2, and coarse

      ! Halo update before tend calculation
      call mpas_timer_start("FB_LTS prognostic halo update")

      call mpas_timer_start("FB_LTS velocity halo update")
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=3)
      call mpas_timer_stop("FB_LTS velocity halo update")

      call mpas_timer_start("FB_LTS thickness halo update")
      call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=3)
      call mpas_timer_stop("FB_LTS thickness halo update")

      call mpas_timer_stop("FB_LTS prognostic halo update")

      ! Compute fast tendencies for thickness
      call mpas_timer_start("FB_LTS compute fast tendencies")
      call mpas_timer_start("FB_LTS compute fast thickness tendencies")

      call ocn_lts_thick_tend(LTSPool, tendPool, &
                              normalVelocityFirstStage, layerThicknessFirstStage, &
                              0, 1, 1, 1, 1)

      call mpas_timer_stop("FB_LTS compute fast thickness tendencies")
      call mpas_timer_stop("FB_LTS compute fast tendencies")

      ! Advance thickness solution with first stage
      call mpas_timer_start("FB_LTS advance solution")
      call mpas_timer_start("FB_LTS advance thickness solution")
      
      ! Interface layers
      do iRegion = 1,nRegions
         do ic = 1, nCellsInLTSRegion(iRegion,2)
            iCell = cellsInLTSRegion(iRegion,2,ic)
            layerThicknessSecondStage(:,iCell) = layerThicknessCur(:,iCell) &
                                                 + weight2nd * dt * layerThicknessTend(:,iCell)
         end do
      end do
      ! Coarse
      do ic = 1, nCellsInLTSRegion(2,1)
         iCell = cellsInLTSRegion(2,1,ic)
         layerThicknessSecondStage(:,iCell) = layerThicknessCur(:,iCell) &
                                              + weight2nd * dt * layerThicknessTend(:,iCell)
      end do
      ! Fine layers close to interface layers
      do ic = 1, nCellsInLTSRegion(1,3)
         iCell = cellsInLTSRegion(1,3,ic)
         layerThicknessSecondStage(:,iCell) = layerThicknessCur(:,iCell) &
                                              + weight2nd * dt * layerThicknessTend(:,iCell)
      end do

      call mpas_timer_stop("FB_LTS advance thickness solution")
      call mpas_timer_stop("FB_LTS advance solution")
      
      ! END: Thickness stage 2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN: Normal velocity stage 2
      ! Compute the second stage of the normal velocity on fine*, 
      ! interface 1, interface 2, and coarse
      
      ! Perform forward-backward average of thickness for stage two of 
      ! FB-RK(3,2):
      ! h = weight2*stage2 + (1-weight2)*current
      layerThicknessForTend(:,:) = config_fb_weight_2 * layerThicknessSecondStage(:,:) &
                                   + (1.0_RKIND - config_fb_weight_2) * layerThicknessCur(:,:)
      
      ! Halo update before tend calculation
      ! Might be able to remove this?
      call mpas_timer_start("FB_LTS prognostic halo update")

      ! Don't need to update normalVelocityFirstStage, previously updated

      call mpas_timer_start("FB_LTS thickness halo update")
      call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=5)
      call mpas_timer_stop("FB_LTS thickness halo update")

      call mpas_timer_stop("FB_LTS prognostic halo update")
 
      ! Compute fast tendencies for normal velocity
      call mpas_timer_start("FB_LTS compute fast tendencies")
      call mpas_timer_start("FB_LTS compute fast velocity tendencies")

      call ocn_lts_vel_tend(LTSPool, tendPool, &
                            normalVelocityFirstStage, layerThicknessForTend, &
                            0, 1, 1, 1, 1)

      normalVelocityTend(:,:) = normalVelocityTend(:,:) &
                                + normalVelocityTendSlow(:,:)

      call mpas_timer_stop("FB_LTS compute fast velocity tendencies")
      call mpas_timer_stop("FB_LTS compute fast tendencies")

      ! Advance normal velocity solution with first stage
      call mpas_timer_start("FB_LTS advance solution")
      call mpas_timer_start("FB_LTS advance velocity solution")
      
      ! Interface regions
      do iRegion = 1,nRegions
         do ie = 1, nEdgesInLTSRegion(iRegion,2)
            iEdge = edgesInLTSRegion(iRegion,2,ie)
            normalVelocitySecondStage(:,iEdge) = normalVelocityCur(:,iEdge) &
                                                 + weight2nd * dt * normalVelocityTend(:,iEdge)
         end do
      end do
      ! Coarse
      do ie = 1, nEdgesInLTSRegion(2,1)
         iEdge = edgesInLTSRegion(2,1,ie)
         normalVelocitySecondStage(:,iEdge) = normalVelocityCur(:,iEdge) &
                                              + weight2nd * dt * normalVelocityTend(:,iEdge)
      end do
      ! Fine layers close to interface layers
      do ie = 1, nEdgesInLTSRegion(1,3)
         iEdge = edgesInLTSRegion(1,3,ie)
         normalVelocitySecondStage(:,iEdge) = normalVelocityCur(:,iEdge) &
                                              + weight2nd * dt * normalVelocityTend(:,iEdge)
      end do
      
      call mpas_timer_stop("FB_LTS advance velocity solution")
      call mpas_timer_stop("FB_LTS advance solution")
      
      ! END: Normal velocity stage 2
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN: Thickness stage 3
      ! Compute the third stage of the thickness on fine*, interface 1, 
      ! interface 2, and coarse

      ! Halo update before tends calculation
      call mpas_timer_start("FB_LTS prognostic halo update")

      call mpas_timer_start("FB_LTS velocity halo update")
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=4)
      call mpas_timer_stop("FB_LTS velocity halo update")

      call mpas_timer_start("FB_LTS thickness halo update")
      call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=4)
      call mpas_timer_stop("FB_LTS thickness halo update")

      call mpas_timer_stop("FB_LTS prognostic halo update")

      ! Compute fast tendencies for thickness
      call mpas_timer_start("FB_LTS compute fast tendencies")
      call mpas_timer_start("FB_LTS compute fast thickness tendencies")

      call ocn_lts_thick_tend(LTSPool, tendPool, & 
                              normalVelocitySecondStage, layerThicknessSecondStage, &
                              0, 1, 1, 1, 1)

      call mpas_timer_stop("FB_LTS compute fast thickness tendencies")
      call mpas_timer_stop("FB_LTS compute fast tendencies")

      ! Advance thickness solution with first stage
      call mpas_timer_start("FB_LTS advance solution")
      call mpas_timer_start("FB_LTS advance thickness solution")
      
      ! Interface layers
      do iRegion = 1,nRegions
         do ic = 1, nCellsInLTSRegion(iRegion,2)
            iCell = cellsInLTSRegion(iRegion,2,ic)
            layerThicknessNew(:,iCell) = layerThicknessCur(:,iCell) &
                                         + dt * layerThicknessTend(:,iCell)
         end do
      end do
      ! Coarse
      do ic = 1, nCellsInLTSRegion(2,1)
         iCell = cellsInLTSRegion(2,1,ic)
         layerThicknessNew(:,iCell) = layerThicknessCur(:,iCell) &
                                      + dt * layerThicknessTend(:,iCell)
      end do
      ! Fine layers close to interface layers
      do ic = 1, nCellsInLTSRegion(1,3)
         iCell = cellsInLTSRegion(1,3,ic)
         layerThicknessNew(:,iCell) = layerThicknessCur(:,iCell) &
                                      + dt * layerThicknessTend(:,iCell)
      end do

      call mpas_timer_stop("FB_LTS advance thickness solution")
      call mpas_timer_stop("FB_LTS advance solution")

      ! END: Thickness stage 3
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN: Normal velocity stage 3
      ! Compute the third stage of the normal velocity on, interface 1, 
      ! and coarse

      ! Perform forward-backward average of thickness for stage three of 
      ! FB-RK(3,2):
      ! h = weight3*new + (1-2*weight3)*secondStage + weight3*current
      layerThicknessForTend(:,:) = config_fb_weight_3 * layerThicknessNew(:,:) &
                                   + (1.0_RKIND - 2.0_RKIND*config_fb_weight_3) * layerThicknessSecondStage(:,:) &
                                   + config_fb_weight_3 * layerThicknessCur(:,:)
      
      ! Halo update before tend calculation
      ! Might be able to remove this?
      call mpas_timer_start("FB_LTS prognostic halo update")

      ! Don't need to update normalVelocitySecondStage, previously updated

      call mpas_timer_start("FB_LTS thickness halo update")
      call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=5)
      call mpas_timer_stop("FB_LTS thickness halo update")

      call mpas_timer_stop("FB_LTS prognostic halo update")
      
      ! Compute fast tendencies for normal velocity
      call mpas_timer_start("FB_LTS compute fast tendencies")
      call mpas_timer_start("FB_LTS compute fast velocity tendencies")

      call ocn_lts_vel_tend(LTSPool, tendPool, &
                            normalVelocitySecondStage, layerThicknessForTend, &
                            0, 0, 1, 0, 1)

      normalVelocityTend(:,:) = normalVelocityTend(:,:) &
                                + normalVelocityTendSlow(:,:)

      call mpas_timer_stop("FB_LTS compute fast velocity tendencies")
      call mpas_timer_stop("FB_LTS compute fast tendencies")

      ! Advance normal velocity solution with first stage
      call mpas_timer_start("FB_LTS advance solution")
      call mpas_timer_start("FB_LTS advance velocity solution")
      
      ! Interface 1
      do ie = 1, nEdgesInLTSRegion(1,2)
         iEdge = edgesInLTSRegion(1,2,ie)
         normalVelocityNew(:,iEdge) = normalVelocityCur(:,iEdge) &
                                      + dt * normalVelocityTend(:,iEdge)
      end do
      ! Coarse
      do ie = 1, nEdgesInLTSRegion(2,1)
         iEdge = edgesInLTSRegion(2,1,ie)
         normalVelocityNew(:,iEdge) = normalVelocityCur(:,iEdge) &
                                      + dt * normalVelocityTend(:,iEdge)
      end do
      
      call mpas_timer_stop("FB_LTS advance velocity solution")
      call mpas_timer_stop("FB_LTS advance solution")

      ! END: Normal velocity stage 3
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END COARSE ADVANCEMENT
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN FINE ADVANCEMENT
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do im = 0, M-1


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! BEGIN: Thickness stage 1
         ! Compute the first stage of the thickness on fine

         ! Calculate predicted 'current' data at intermediate 
         ! time-level on interface 1:
         ! u,h = (im/M)*new + (1- (im/M))*current
         call mpas_timer_start('FB_LTS interface prediction')

         normalVelocityForTend(:,:) = normalVelocityCur(:,:)
         do ie = 1, nEdgesInLTSRegion(1,2)
            iEdge = edgesInLTSRegion(1,2,ie)
            normalVelocityForTend(:,iEdge) = (REAL(im)/M) * normalVelocityNew(:,iEdge) &
                                             + (1.0_RKIND - REAL(im)/M) * normalVelocityCur(:,iEdge)
         end do

         layerThicknessForTend(:,:) = layerThicknessCur(:,:)
         do ic = 1, nCellsInLTSRegion(1,2)
            iCell = cellsInLTSRegion(1,2,ic)
            layerThicknessForTend(:,iCell) = (REAL(im)/M) * layerThicknessNew(:,iCell) &
                                             + (1.0_RKIND - REAL(im)/M) * layerThicknessCur(:,iCell)
         end do

         call mpas_timer_stop('FB_LTS interface prediction')
         
         ! Halo update before tend calculation
         call mpas_timer_start("FB_LTS prognostic halo update")

         call mpas_timer_start("FB_LTS velocity halo update")
         call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=5)
         call mpas_timer_stop("FB_LTS velocity halo update")

         call mpas_timer_start("FB_LTS thickness halo update")
         call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=5)
         call mpas_timer_stop("FB_LTS thickness halo update")

         call mpas_timer_stop("FB_LTS prognostic halo update")
         
         ! Compute fast tendencies for thickness
         call mpas_timer_start("FB_LTS compute fast tendencies")
         call mpas_timer_start("FB_LTS compute fast thickness tendencies")

         call ocn_lts_thick_tend(LTSPool, tendPool, & 
                                 normalVelocityForTend, layerThicknessForTend, &
                                 1, 1, 0, 0, 0)

         call mpas_timer_stop("FB_LTS compute fast thickness tendencies")
         call mpas_timer_stop("FB_LTS compute fast tendencies")

         ! Advance thickness solution with first stage
         call mpas_timer_start("FB_LTS advance solution")
         call mpas_timer_start("FB_LTS advance thickness solution")
         
         ! Fine cells adjacent to interface 1
         do ic = 1, nCellsInLTSRegion(1,3)
            iCell = cellsInLTSRegion(1,3,ic)
            layerThicknessFirstStage(:,iCell) = layerThicknessCur(:,iCell) &
                                                + weight1st * dtFine * layerThicknessTend(:,iCell)
         end do
         ! Fine in the interior, away from interface 1
         do ic = 1, nCellsInLTSRegion(1,1)
            iCell = cellsInLTSRegion(1,1,ic)
            layerThicknessFirstStage(:,iCell) = layerThicknessCur(:,iCell) &
                                                + weight1st * dtFine * layerThicknessTend(:,iCell)
         end do

         call mpas_timer_stop("FB_LTS advance thickness solution")
         call mpas_timer_stop("FB_LTS advance solution")

         ! END: Thickness stage 1
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! BEGIN: Normal velocity stage 1
         ! Compute the first stage of the normal velocity on fine

         ! Perform forward-backward average of thickness for stage one of 
         ! FB-RK(3,2):
         ! h = weight1*stage1 + (1-weight1)*current
         layerThicknessForTend(:,:) = config_fb_weight_1 * layerThicknessFirstStage(:,:) &
                                      + (1.0_RKIND - config_fb_weight_1) * layerThicknessCur(:,:)

         ! Calculate predicted FB averaged thickness data at 
         ! intermediate time-level on interface 1:
         ! h = weight1*'firstStage' + (1-weight1)*'current'
         !   = weight1*( (im/M)*new + (1/M)*firstStage + (1-(im+1)/M)*current )
         !     + (1-weight1)*( (im/M)*new + (1-im/M)*current )
         call mpas_timer_start('FB_LTS interface prediction')
         
         do ic = 1, nCellsInLTSRegion(1,2)
            iCell = cellsInLTSRegion(1,2,ic)
            layerThicknessForTend(:,iCell) = config_fb_weight_1 &
                                             * ( (REAL(im)/M) * layerThicknessNew(:,iCell) &
                                               + (1.0_RKIND/M) * layerThicknessFirstStage(:,iCell) &
                                               + (1.0_RKIND - (REAL(im)+1.0_RKIND)/M) * layerThicknessCur(:,iCell) ) &
                                             + (1.0_RKIND - config_fb_weight_1) &
                                             * ( (REAL(im)/M) * layerThicknessNew(:,iCell) &
                                               + (1.0_RKIND - REAL(im)/M) * layerThicknessCur(:,iCell) )
         end do
         
         call mpas_timer_stop('FB_LTS interface prediction')

         ! Halo update before tend calculation
         call mpas_timer_start("FB_LTS prognostic halo update")

         ! Don't need to update normalVelocityForTend, previously updated

         call mpas_timer_start("FB_LTS thickness halo update")
         call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=5)
         call mpas_timer_stop("FB_LTS thickness halo update")

         call mpas_timer_stop("FB_LTS prognostic halo update")
         
         ! Compute fast tendencies for normal velocity
         call mpas_timer_start("FB_LTS compute fast tendencies")
         call mpas_timer_start("FB_LTS compute fast velocity tendencies")

         call ocn_lts_vel_tend(LTSPool, tendPool, &
                               normalVelocityForTend, layerThicknessForTend, &
                               1, 1, 0, 0, 0)

         normalVelocityTend(:,:) = normalVelocityTend(:,:) &
                                   + normalVelocityTendSlow(:,:)

         call mpas_timer_stop("FB_LTS compute fast velocity tendencies")
         call mpas_timer_stop("FB_LTS compute fast tendencies")
   
         ! Advance normal velocity solution with first stage
         call mpas_timer_start("FB_LTS advance solution")
         call mpas_timer_start("FB_LTS advance velocity solution")

         ! Fine cells adjacent to interface 1
         do ie = 1, nEdgesInLTSRegion(1,3)
            iEdge = edgesInLTSRegion(1,3,ie)
            normalVelocityFirstStage(:,iEdge) = normalVelocityCur(:,iEdge) &
                                                + weight1st * dtFine * normalVelocityTend(:,iEdge)
         end do
         ! Fine in the interior, away from interface 1
         do ie = 1, nEdgesInLTSRegion(1,1)
            iEdge = edgesInLTSRegion(1,1,ie)
            normalVelocityFirstStage(:,iEdge) = normalVelocityCur(:,iEdge) &
                                                + weight1st * dtFine * normalVelocityTend(:,iEdge)
         end do
         
         call mpas_timer_stop("FB_LTS advance velocity solution")
         call mpas_timer_stop("FB_LTS advance solution")

         ! END: Normal velocity stage 1
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! BEGIN: Thickness stage 2
         ! Compute the second stage of the thickness on fine

         ! Calculate predicted stage 1 data at intermediate 
         ! time-level on interface 1:
         ! u,h = (im/M)*new + (1/M)*firstStage + (1-(im+1)/M)*current
         call mpas_timer_start('FB_LTS interface prediction')

         normalVelocityForTend(:,:) = normalVelocityFirstStage(:,:)
         do ie = 1, nEdgesInLTSRegion(1,2)
            iEdge = edgesInLTSRegion(1,2,ie)
            normalVelocityForTend(:,iEdge) = (REAL(im)/M) * normalVelocityNew(:,iEdge) &
                                             + (1.0_RKIND/M) * normalVelocityFirstStage(:,iEdge) &
                                             + (1.0_RKIND - (REAL(im)+1.0_RKIND)/M) * normalVelocityCur(:,iEdge)
         end do

         layerThicknessForTend(:,:) = layerThicknessFirstStage(:,:)
         do ic = 1, nCellsInLTSRegion(1,2)
            iCell = cellsInLTSRegion(1,2,ic)
            layerThicknessForTend(:,iCell) = (REAL(im)/M) * layerThicknessNew(:,iCell) &
                                             + (1.0_RKIND/M) * layerThicknessFirstStage(:,iCell) &
                                             + (1.0_RKIND - (REAL(im)+1.0_RKIND)/M) * layerThicknessCur(:,iCell)
         end do

         call mpas_timer_stop('FB_LTS interface prediction')
         
         ! Halo update before tend calculation
         call mpas_timer_start("FB_LTS prognostic halo update")

         call mpas_timer_start("FB_LTS velocity halo update")
         call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=5)
         call mpas_timer_stop("FB_LTS velocity halo update")

         call mpas_timer_start("FB_LTS thickness halo update")
         call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=5)
         call mpas_timer_stop("FB_LTS thickness halo update")

         call mpas_timer_stop("FB_LTS prognostic halo update")
         
         ! Compute fast tendencies for thickness
         call mpas_timer_start("FB_LTS compute fast tendencies")
         call mpas_timer_start("FB_LTS compute fast thickness tendencies")

         call ocn_lts_thick_tend(LTSPool, tendPool, & 
                                 normalVelocityForTend, layerThicknessForTend, &
                                 1, 1, 0, 0, 0)

         call mpas_timer_stop("FB_LTS compute fast thickness tendencies")
         call mpas_timer_stop("FB_LTS compute fast tendencies")

         ! Advance thickness solution with second stage
         call mpas_timer_start("FB_LTS advance solution")
         call mpas_timer_start("FB_LTS advance thickness solution")
         
         ! Fine cells adjacent to interface 1
         do ic = 1, nCellsInLTSRegion(1,3)
            iCell = cellsInLTSRegion(1,3,ic)
            layerThicknessSecondStage(:,iCell) = layerThicknessCur(:,iCell) &
                                                 + weight2nd * dtFine * layerThicknessTend(:,iCell)
         end do
         ! Fine in the interior, away from interface 1
         do ic = 1, nCellsInLTSRegion(1,1)
            iCell = cellsInLTSRegion(1,1,ic)
            layerThicknessSecondStage(:,iCell) = layerThicknessCur(:,iCell) &
                                                 + weight2nd * dtFine * layerThicknessTend(:,iCell)
         end do

         call mpas_timer_stop("FB_LTS advance thickness solution")
         call mpas_timer_stop("FB_LTS advance solution")

         ! END: Thickness stage 2
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! BEGIN: Normal velocity stage 2
         ! Compute the second stage of the normal velocity on fine

         ! Perform forward-backward average of thickness for stage two of 
         ! FB-RK(3,2):
         ! h = weight2*stage2 + (1-weight2)*current
         layerThicknessForTend(:,:) = config_fb_weight_2 * layerThicknessSecondStage(:,:) &
                                      + (1.0_RKIND - config_fb_weight_2) * layerThicknessCur(:,:)

         ! Calculate predicted FB averaged thickness data at 
         ! intermediate time-level on interface 1:
         ! h = weight2*'secondStage' + (1-weight2)*'current'
         !   = weight2*( (im/M)*new + (1/M)*secondStage + (1-(im+1)/M)*current )
         !     + (1-weight2)*( (im/M)*new + (1-im/M)*current )
         call mpas_timer_start('FB_LTS interface prediction')
         
         do ic = 1, nCellsInLTSRegion(1,2)
            iCell = cellsInLTSRegion(1,2,ic)
            layerThicknessForTend(:,iCell) = config_fb_weight_2 &
                                             * ( (REAL(im)/M) * layerThicknessNew(:,iCell) &
                                               + (1.0_RKIND/M) * layerThicknessSecondStage(:,iCell) &
                                               + (1.0_RKIND - (REAL(im)+1.0_RKIND)/M) * layerThicknessCur(:,iCell) ) &
                                             + (1.0_RKIND - config_fb_weight_2) &
                                             * ( (REAL(im)/M) * layerThicknessNew(:,iCell) &
                                               + (1.0_RKIND - REAL(im)/M) * layerThicknessCur(:,iCell) )
         end do
         
         call mpas_timer_stop('FB_LTS interface prediction')

         ! Halo update before tend calculation
         call mpas_timer_start("FB_LTS prognostic halo update")

         ! Don't need to update normalVelocityForTend, previously updated

         call mpas_timer_start("FB_LTS thickness halo update")
         call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=5)
         call mpas_timer_stop("FB_LTS thickness halo update")

         call mpas_timer_stop("FB_LTS prognostic halo update")
         
         ! Compute fast tendencies for normal velocity
         call mpas_timer_start("FB_LTS compute fast tendencies")
         call mpas_timer_start("FB_LTS compute fast velocity tendencies")

         call ocn_lts_vel_tend(LTSPool, tendPool, &
                               normalVelocityForTend, layerThicknessForTend, &
                               1, 1, 0, 0, 0)

         normalVelocityTend(:,:) = normalVelocityTend(:,:) &
                                   + normalVelocityTendSlow(:,:)

         call mpas_timer_stop("FB_LTS compute fast velocity tendencies")
         call mpas_timer_stop("FB_LTS compute fast tendencies")

         ! Advance normal velocity solution with second stage
         call mpas_timer_start("FB_LTS advance solution")
         call mpas_timer_start("FB_LTS advance velocity solution")
         
         ! Fine cells adjacent to interface 1
         do ie = 1, nEdgesInLTSRegion(1,3)
            iEdge = edgesInLTSRegion(1,3,ie)
            normalVelocitySecondStage(:,iEdge) = normalVelocityCur(:,iEdge) &
                                                 + weight2nd * dtFine * normalVelocityTend(:,iEdge)
         end do
         ! Fine in the interior, away from interface 1
         do ie = 1, nEdgesInLTSRegion(1,1)
            iEdge = edgesInLTSRegion(1,1,ie)
            normalVelocitySecondStage(:,iEdge) = normalVelocityCur(:,iEdge) &
                                                 + weight2nd * dtFine * normalVelocityTend(:,iEdge)
         end do
         
         call mpas_timer_stop("FB_LTS advance velocity solution")
         call mpas_timer_stop("FB_LTS advance solution")

         ! END: Normal velocity stage 2
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! BEGIN: Thickness stage 3
         ! Compute the third stage of the thickness on fine and 
         ! increment the third stage interface correction terms

         ! Calculate predicted stage 2 data at intermediate 
         ! time-level on interface 1:
         ! u,h = (im/M)*new + (1/M)*secondStage + (1-(im+1)/M)*current
         call mpas_timer_start('FB_LTS interface prediction')

         normalVelocityForTend(:,:) = normalVelocitySecondStage(:,:)
         do ie = 1, nEdgesInLTSRegion(1,2)
            iEdge = edgesInLTSRegion(1,2,ie)
            normalVelocityForTend(:,iEdge) = (REAL(im)/M) * normalVelocityNew(:,iEdge) &
                                             + (1.0_RKIND/M) * normalVelocitySecondStage(:,iEdge) &
                                             + (1.0_RKIND - (REAL(im)+1.0_RKIND)/M) * normalVelocityCur(:,iEdge)
         end do

         layerThicknessForTend(:,:) = layerThicknessSecondStage(:,:)
         do ic = 1, nCellsInLTSRegion(1,2)
            iCell = cellsInLTSRegion(1,2,ic)
            layerThicknessForTend(:,iCell) = (REAL(im)/M) * layerThicknessNew(:,iCell) &
                                             + (1.0_RKIND/M) * layerThicknessSecondStage(:,iCell) &
                                             + (1.0_RKIND - (REAL(im)+1.0_RKIND)/M) * layerThicknessCur(:,iCell)
         end do

         call mpas_timer_stop('FB_LTS interface prediction')
         
         ! Halo update before tend calculation
         call mpas_timer_start("FB_LTS prognostic halo update")

         call mpas_timer_start("FB_LTS velocity halo update")
         call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=5)
         call mpas_timer_stop("FB_LTS velocity halo update")

         call mpas_timer_start("FB_LTS thickness halo update")
         call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=5)
         call mpas_timer_stop("FB_LTS thickness halo update")

         call mpas_timer_stop("FB_LTS prognostic halo update")
         
         ! Compute fast tendencies for thickness
         call mpas_timer_start("FB_LTS compute fast tendencies")
         call mpas_timer_start("FB_LTS compute fast thickness tendencies")

         call ocn_lts_thick_tend(LTSPool, tendPool, & 
                                 normalVelocityForTend, layerThicknessForTend, &
                                 1, 1, 1, 1, 0)

         call mpas_timer_stop("FB_LTS compute fast thickness tendencies")
         call mpas_timer_stop("FB_LTS compute fast tendencies")

         ! Advance thickness solution with third stage
         call mpas_timer_start("FB_LTS advance solution")
         call mpas_timer_start("FB_LTS advance thickness solution")

         ! Increment interface correction terms
         do iRegion = 1,nRegions
            do ic = 1, nCellsInLTSRegion(iRegion,2)
               iCell = cellsInLTSRegion(iRegion,2,ic)
               layerThicknessTendSum3rd(:,iCell) = layerThicknessTendSum3rd(:,iCell) &
                                                   + layerThicknessTend(:,iCell)
            end do
         end do
         ! Fine cells adjacent to interface 1
         do ic = 1, nCellsInLTSRegion(1,3)
            iCell = cellsInLTSRegion(1,3,ic)
            layerThicknessNew(:,iCell) = layerThicknessCur(:,iCell) &
                                          + dtFine * layerThicknessTend(:,iCell)
         end do
         ! Fine in the interior, away from interface 1
         do ic = 1, nCellsInLTSRegion(1,1)
            iCell = cellsInLTSRegion(1,1,ic)
            layerThicknessNew(:,iCell) = layerThicknessCur(:,iCell) &
                                          + dtFine * layerThicknessTend(:,iCell)
         end do

         call mpas_timer_stop("FB_LTS advance thickness solution")
         call mpas_timer_stop("FB_LTS advance solution")

         ! END: Thickness stage 3
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! BEGIN: Normal velocity stage 3
         ! Compute the third stage of the normal velocity on fine and 
         ! increment the third stage interface correction terms

         ! Perform forward-backward average of thickness for stage three of 
         ! FB-RK(3,2):
         ! h = weight3*new + (1-2*weight3)*stageTwo + weight3*current
         layerThicknessForTend(:,:) = config_fb_weight_3 * layerThicknessNew(:,:) &
                                      + (1.0_RKIND - 2.0_RKIND*config_fb_weight_3) * layerThicknessSecondStage(:,:) &
                                      + config_fb_weight_3 * layerThicknessCur(:,:)

         ! Calculate predicted FB averaged thickness data at 
         ! intermediate time-level on interface 1:
         ! h = weight3*new + (1-2*weight3)*'secondStage' + weight3*'current'
         !   = weight3*( ((im+1)/M)*new + (1-(im+1)/M)*current )
         !     + (1-2*weight3)*( (im/M)*new + (1/M)*secondStage + (1-(im+1)/M)*current )
         !     + weight3*( (im/M)*new + (1-im/M)*current )
         call mpas_timer_start('FB_LTS interface prediction')
         
         do ic = 1, nCellsInLTSRegion(1,2)
            iCell = cellsInLTSRegion(1,2,ic)
            layerThicknessForTend(:,iCell) = config_fb_weight_3 &
                                             * ( ((REAL(im)+1.0_RKIND)/M) * layerThicknessNew(:,iCell) &
                                               + (1.0_RKIND - (REAL(im)+1.0_RKIND)/M) * layerThicknessCur(:,iCell) ) &
                                             + (1.0_RKIND - 2.0_RKIND*config_fb_weight_3) &
                                             * ( (REAL(im)/M) * layerThicknessNew(:,iCell) &
                                               + (1.0_RKIND/M) * layerThicknessSecondStage(:,iCell) &
                                               + (1.0_RKIND - (REAL(im)+1.0_RKIND)/M) * layerThicknessCur(:,iCell) ) &
                                             + config_fb_weight_3 &
                                             * ( (REAL(im)/M) * layerThicknessNew(:,iCell) &
                                               + (1.0_RKIND - REAL(im)/M) * layerThicknessCur(:,iCell) )
         end do
         
         call mpas_timer_stop('FB_LTS interface prediction')

         ! Halo update before tend calculation
         call mpas_timer_start("FB_LTS prognostic halo update")
         
         ! Don't need to update normalVelocityForTend, previously updated

         call mpas_timer_start("FB_LTS thickness halo update")
         call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=5)
         call mpas_timer_stop("FB_LTS thickness halo update")

         call mpas_timer_stop("FB_LTS prognostic halo update")
         
         ! Compute fast tendencies for normal velocity
         call mpas_timer_start("FB_LTS compute fast tendencies")
         call mpas_timer_start("FB_LTS compute fast velocity tendencies")

         call ocn_lts_vel_tend(LTSPool, tendPool, &
                               normalVelocityForTend, layerThicknessForTend, &
                               1, 1, 1, 1, 0)

         normalVelocityTend(:,:) = normalVelocityTend(:,:) &
                                   + normalVelocityTendSlow(:,:)

         call mpas_timer_stop("FB_LTS compute fast velocity tendencies")
         call mpas_timer_stop("FB_LTS compute fast tendencies")

         ! Advance normal velocity solution with third stage
         ! Here, we write into normalVelocityCur, rather than
         ! normalVelocityNew, so that the New data is stored
         ! in normalVelocityCur to be used at the beginning of the
         ! next subcycled time-step
         call mpas_timer_start("FB_LTS advance solution")
         call mpas_timer_start("FB_LTS advance velocity solution")
         
         ! Increment interface correction terms
         do iRegion = 1,nRegions
            do ie = 1, nEdgesInLTSRegion(iRegion,2)
               iEdge = edgesInLTSRegion(iRegion,2,ie)
               normalVelocityTendSum3rd(:,iEdge) = normalVelocityTendSum3rd(:,iEdge) &
                                                   + normalVelocityTend(:,iEdge)
            end do
         end do
         ! Fine cells adjacent to interface 1
         do ie = 1, nEdgesInLTSRegion(1,3)
            iEdge = edgesInLTSRegion(1,3,ie)
            normalVelocityCur(:,iEdge) = normalVelocityCur(:,iEdge) &
                                         + dtFine * normalVelocityTend(:,iEdge)
         end do
         ! Fine in the interior, away from interface 1
         do ie = 1, nEdgesInLTSRegion(1,1)
            iEdge = edgesInLTSRegion(1,1,ie)
            normalVelocityCur(:,iEdge) = normalVelocityCur(:,iEdge) &
                                         + dtFine * normalVelocityTend(:,iEdge)
         end do
         
         call mpas_timer_stop("FB_LTS advance velocity solution")
         call mpas_timer_stop("FB_LTS advance solution")

         ! END: Normal velocity stage 3
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         
         ! Copy new solution into cur solution for use at the beginning
         ! of the next substep
         ! We only do this for layerThickness here because we write 
         ! into normalVelocityCur rather than normalVelocityNew above
         ! for stage 3 of the subcycled normal velocity
         call mpas_timer_start("FB_LTS copy soln")
         
         ! Fine adjacent to interface 1
         do ic = 1, nCellsInLTSRegion(1,3)
            iCell = cellsInLTSRegion(1,3,ic)
            layerThicknessCur(:,iCell) = layerThicknessNew(:,iCell)
         end do
         ! Fine in the interior, away from interface 1
         do ic = 1, nCellsInLTSRegion(1,1)
            iCell = cellsInLTSRegion(1,1,ic)
            layerThicknessCur(:,iCell) = layerThicknessNew(:,iCell)
         end do

         call mpas_timer_stop("FB_LTS copy soln")

      end do  ! end of fine time-step subcycling loop
      
      
      ! Copy data from normalVelocityCur into normalVelocityNew
      ! which currently holds the New data. 
      call mpas_timer_start("FB_LTS copy soln")

      ! Fine adjacent to interface 1
      do ie = 1, nEdgesInLTSRegion(1,3)
         iEdge = edgesInLTSRegion(1,3,ie)
         normalVelocityNew(:,iEdge) = normalVelocityCur(:,iEdge)
      end do
      ! Fine in the interior, away from interface 1
      do ie = 1, nEdgesInLTSRegion(1,1)
         iEdge = edgesInLTSRegion(1,1,ie)
         normalVelocityNew(:,iEdge) = normalVelocityCur(:,iEdge)
      end do
      
      call mpas_timer_stop("FB_LTS copy soln")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END FINE ADVANCEMENT
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN INTERFACE CORRECTION
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Calculate the corrected solution in the interface layers
      call mpas_timer_start("FB_LTS interface correction")
   
      ! Perform the correction on interface 1 and 2
      do iRegion = 1, nRegions
         ! Normal velocity correction
         do ie = 1, nEdgesInLTSRegion(iRegion,2)
            iEdge = edgesInLTSRegion(iRegion,2,ie)
            normalVelocityNew(:,iEdge) = normalVelocityCur(:,iEdge) &
                                         + dtFine * normalVelocityTendSum3rd(:,iEdge)
         end do
         ! Layer thickness correction
         do ic = 1, nCellsInLTSRegion(iRegion,2)
            iCell = cellsInLTSRegion(iRegion,2,ic)
            layerThicknessNew(:,iCell) = layerThicknessCur(:,iCell) &
                                         + dtFine * layerThicknessTendSum3rd(:,iCell)
         end do
      end do

      call mpas_timer_stop("FB_LTS interface correction")

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! END INTERFACE CORRECTION
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      
      ! Final prognostic halo update, is this needed? 
      call mpas_timer_start("FB_LTS prognostic halo update")

      call mpas_timer_start("FB_LTS velocity halo update")
      call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=2)
      call mpas_timer_stop("FB_LTS velocity halo update")
      
      call mpas_timer_start("FB_LTS thickness halo update")
      call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=2)
      call mpas_timer_stop("FB_LTS thickness halo update")

      call mpas_timer_stop("FB_LTS prognostic halo update")

      ! Implicit vmix, DO WE NEED THIS?
      call mpas_timer_start("FB_LTS implicit vmix")
      if (.not. config_disable_vel_vmix) then
         call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, &
                                   verticalMeshPool, scratchPool, tracersPool, 2)
         call ocn_vmix_implicit(dt, meshPool, statePool, forcingPool, scratchPool, err, 2)

         ! Halo update
         call mpas_timer_start("FB_LTS prognostic halo update")

         call mpas_timer_start("FB_LTS velocity halo update")
         call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=2)
         call mpas_timer_stop("FB_LTS velocity halo update")
      
         call mpas_timer_start("FB_LTS thickness halo update")
         call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=2)
         call mpas_timer_stop("FB_LTS thickness halo update")

         call mpas_timer_stop("FB_LTS prognostic halo update")
      end if   
      call mpas_timer_stop("FB_LTS implicit vmix")

      call mpas_timer_stop("FB_LTS main loop")
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! END FB_LTS SCHEME
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      

      call mpas_timer_start("FB_LTS cleanup")

      if (config_prescribe_velocity) then
         do iEdge = 1, nEdgesAll
            normalVelocityNew(:,iEdge) = normalVelocityCur(:,iEdge)
         end do
      end if

      if (config_prescribe_thickness) then
         do iCell = 1, nCellsAll
            layerThicknessNew(:,iCell) = layerThicknessCur(:,iCell)
         end do
      end if

      do iEdge = 1, nEdgesAll
         normalTransportVelocity(:,iEdge) = normalVelocityNew(:,iEdge)
      end do

      call mpas_reconstruct(meshPool, normalVelocityNew, &
                            velocityX, velocityY, velocityZ, &
                            velocityZonal, velocityMeridional, &
                            includeHalos = .true.)

      call mpas_reconstruct(meshPool, gradSSH, &
                            gradSSHX, gradSSHY, gradSSHZ, &
                            gradSSHZonal, gradSSHMeridional, &
                            includeHalos = .true.)

      do iCell = 1, nCellsAll
         surfaceVelocity(indexSurfaceVelocityZonal, iCell) = velocityZonal(1, iCell)
         surfaceVelocity(indexSurfaceVelocityMeridional, iCell) = velocityMeridional(1, iCell)

         SSHGradient(indexSSHGradientZonal, iCell) = gradSSHZonal(iCell)
         SSHGradient(indexSSHGradientMeridional, iCell) = gradSSHMeridional(iCell)
      end do

      call ocn_time_average_coupled_accumulate(statePool, forcingPool, 2)

      do iCell = 1, nCellsAll
         surfaceVelocity(indexSurfaceVelocityZonal, iCell) = velocityZonal(1, iCell)
         surfaceVelocity(indexSurfaceVelocityMeridional, iCell) = velocityMeridional(1, iCell)

         SSHGradient(indexSSHGradientZonal, iCell) = gradSSHZonal(iCell)
         SSHGradient(indexSSHGradientMeridional, iCell) = gradSSHMeridional(iCell)
      end do

      ! Update diagnostics
      call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, &
                                verticalMeshPool, scratchPool, tracersPool, 2)

      call mpas_timer_stop("FB_LTS cleanup")

   end subroutine ocn_time_integrator_fblts!}}}


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_fblts_init
!
!> \brief  Initialization for MPAS ocean FB_LTS time integration scheme
!> \author Giacomo Capodaglio
!> \date   November 2022
!> \details
!>  This routine computes the LTS arrays
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integration_fblts_init(domain)!{{{
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Sets the LTS arrays
   !
   ! Output: LTS arrays are written    
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      implicit none
      
      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: &
         domain  !< Input/output: model state
      
      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------
      
      type (block_type), pointer :: &
         block
      
      type (mpas_pool_type), pointer :: &
         LTSPool

      type (mpas_pool_type), pointer :: &
         tendSlowPool, &
         tendSum3rdPool, &
         prevTendSlowPool, nextTendSlowPool, &
         prevTendSum3rdPool, nextTendSum3rdPool

      type (mpas_pool_type), pointer :: &
         tendPool
      
      integer, dimension(:), allocatable :: &
         isLTSRegionEdgeAssigned
      
      integer, dimension(:), pointer :: &
         LTSRegion
      
      integer, dimension(:,:), pointer :: &
         nCellsInLTSRegion, &
         nEdgesInLTSRegion
      
      integer, dimension(:,:,:), pointer :: &
         cellsInLTSRegion, &
         edgesInLTSRegion
      
      integer, dimension(2) :: &
         minMaxLTSRegion
      
      integer :: &
         i, iCell, iEdge, iRegion, coarseRegions, &
         fineRegions, fineRegionsM1
      
      !-----------------------------------------------------------------
      ! Begin routine
      !-----------------------------------------------------------------
      
      minMaxLTSRegion(1) = 1
      minMaxLTSRegion(2) = 2

      block => domain % blocklist
      call mpas_pool_get_subpool(block%structs, 'tend', tendPool)

      call mpas_pool_create_pool(tendSum3rdPool)
      call mpas_pool_clone_pool(tendPool, tendSum3rdPool, 1)
      call mpas_pool_create_pool(tendSlowPool)
      call mpas_pool_clone_pool(tendPool, tendSlowPool, 1)

      call mpas_pool_add_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)
      call mpas_pool_add_subpool(block % structs, 'tend_slow', tendSlowPool)

      if (associated(block % prev)) then
         call mpas_pool_get_subpool(block % prev % structs, 'tend_sum_3rd', tendSum3rdPool)
         call mpas_pool_get_subpool(block % prev % structs, 'tend_slow', tendSlowPool)
      else
         nullify(prevTendSum3rdPool)
         nullify(prevTendSlowPool)
      end if

      if (associated(block % next)) then
         call mpas_pool_get_subpool(block % next % structs, 'tend_sum_3rd', nextTendSum3rdPool)
         call mpas_pool_get_subpool(block % next % structs, 'tend_slow', nextTendSlowPool)
      else
         nullify(nextTendSum3rdPool)
         nullify(nextTendSlowPool)
      end if

      call mpas_pool_get_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)
      call mpas_pool_get_subpool(block % structs, 'tend_slow', tendSlowPool)
   
      if (associated(prevTendSum3rdPool) .and. associated(nextTendSum3rdPool)) then
         call mpas_pool_link_pools(tendSum3rdPool, prevTendSum3rdPool, nextTendSum3rdPool)
      else if (associated(prevTendSum3rdPool)) then
         call mpas_pool_link_pools(tendSum3rdPool, prevTendSum3rdPool)
      else if (associated(nextTendSum3rdPool)) then
         call mpas_pool_link_pools(tendSum3rdPool,nextPool=nextTendSum3rdPool)
      else
         call mpas_pool_link_pools(tendSum3rdPool)
      end if

      if (associated(prevTendSlowPool) .and. associated(nextTendSlowPool)) then
         call mpas_pool_link_pools(tendSlowPool, prevTendSlowPool, nextTendSlowPool)
      else if (associated(prevTendSlowPool)) then
         call mpas_pool_link_pools(tendSlowPool, prevTendSlowPool)
      else if (associated(nextTendSlowPool)) then
         call mpas_pool_link_pools(tendSlowPool,nextPool=nextTendSlowPool)
      else
         call mpas_pool_link_pools(tendSlowPool)
      end if

      call mpas_pool_link_parinfo(block, tendSum3rdPool)
      call mpas_pool_link_parinfo(block, tendSlowPool)

      call mpas_pool_get_subpool(block % structs, 'LTS', LTSPool)

      call mpas_pool_get_array(LTSPool, 'LTSRegion', LTSRegion)
      call mpas_pool_get_array(LTSPool, 'cellsInLTSRegion', cellsInLTSRegion)
      call mpas_pool_get_array(LTSPool, 'nCellsInLTSRegion', nCellsInLTSRegion)
      call mpas_pool_get_array(LTSPool, 'edgesInLTSRegion', edgesInLTSRegion)
      call mpas_pool_get_array(LTSPool, 'nEdgesInLTSRegion', nEdgesInLTSRegion)

      ! LTS Regions code:
      ! 1 = fine
      ! 2 = coarse
      ! 3 = interface layer 1
      ! 4 = interface layer 2
      ! 5 = fine (to advance when doing 1st, 2nd, 3rd stage on interface)
      
      nCellsInLTSRegion(:,:) = 0
      nEdgesInLTSRegion(:,:) = 0

      ! This is a loop to build the lists of elements in the fine, coarse, 
      ! and interface regions. Only loops up to nCellsOwned because in the 
      ! time-stepping we only want to advance the cells owned by the MPI process
      do iCell = 1, nCellsOwned
         do iRegion = 1,2
            if (iRegion == minMaxLTSRegion(iRegion)) then
               if(LTSRegion(iCell) == minMaxLTSRegion(iRegion)) then
                   nCellsInLTSRegion(iRegion,1) = nCellsInLTSRegion(iRegion,1) + 1
                   cellsInLTSRegion(iRegion,1,nCellsInLTSRegion(iRegion,1)) = iCell
               end if
               if(LTSRegion(iCell) == (minMaxLTSRegion(iRegion) + 2) ) then
                  nCellsInLTSRegion(iRegion,2) = nCellsInLTSRegion(iRegion,2) + 1
                  cellsInLTSRegion(iRegion,2,nCellsInLTSRegion(iRegion,2)) = iCell
               end if
            end if
         end do
         if (LTSRegion(iCell) == 5) then
            nCellsInLTSRegion(1,3) = nCellsInLTSRegion(1,3) + 1
            cellsInLTSRegion(1,3,nCellsInLTSRegion(1,3)) = iCell
         end if
      end do
      
      ! Below we fill out the lists for the edges, according to the LTSRegion 
      ! that have been assigned to the cells.  We move from the fine to the 
      ! coarse (i.e. from the fine to the nearest LTS region in the direction 
      ! of the coarse). Note that edges shared between cells of different LTS 
      ! regions are owned by the cell in the LTS region closest to the fine 
      ! region, see Figure 3 in "Conservative explicit local time-stepping 
      ! schemes for the shallow water equations" by Hoang et al. (halo edges 
      ! however are owned by whatever processor they are initially assigned to).

      allocate(isLTSRegionEdgeAssigned(nEdgesOwned))
      isLTSRegionEdgeAssigned(:) = 0

      do iCell = 1, nCellsInLTSRegion(1,1)
         do i = 1, nEdgesOnCell(cellsInLTSRegion(1,1,iCell))
            iEdge = edgesOnCell(i,cellsInLTSRegion(1,1,iCell))
            if (iEdge .le. nEdgesOwned) then
               if (isLTSRegionEdgeAssigned(iEdge) == 0) then
                  nEdgesInLTSRegion(1,1) = nEdgesInLTSRegion(1,1) + 1
                  edgesInLTSRegion(1,1, nEdgesInLTSRegion(1,1)) = iEdge
                  isLTSRegionEdgeAssigned(iEdge) = 1
               end if
            end if
         end do
      end do

      fineRegions = 3
      fineRegionsM1 = 2
      do iRegion = 1, fineRegionsM1
         do iCell = 1, nCellsInLTSRegion(1, fineRegions - iRegion + 1)
            do i = 1, nEdgesOnCell(cellsInLTSRegion(1, fineRegions - iRegion + 1, iCell))
               iEdge = edgesOnCell(i,cellsInLTSRegion(1, fineRegions - iRegion + 1, iCell))
               if (iEdge .le. nEdgesOwned) then
                  if (isLTSRegionEdgeAssigned(iEdge) == 0) then
                     nEdgesInLTSRegion(1, fineRegions - iRegion + 1) &
                        = nEdgesInLTSRegion(1, fineRegions - iRegion + 1) + 1
                     edgesInLTSRegion(1, fineRegions - iRegion + 1, &
                                      nEdgesInLTSRegion(1, fineRegions - iRegion + 1)) &
                        = iEdge
                     isLTSRegionEdgeAssigned(iEdge) = 1
                  end if
              end if
            end do
         end do
      end do

      coarseRegions = 2
      do iRegion = 1, coarseRegions
         do iCell = 1, nCellsInLTSRegion(2, coarseRegions - iRegion + 1)
            do i = 1, nEdgesOnCell(cellsInLTSRegion(2,coarseRegions - iRegion + 1,iCell))
               iEdge = edgesOnCell(i,cellsInLTSRegion(2,coarseRegions - iRegion + 1,iCell))
               if (iEdge .le. nEdgesOwned) then
                  if (isLTSRegionEdgeAssigned(iEdge) == 0) then
                     nEdgesInLTSRegion(2, coarseRegions - iRegion + 1) &
                        = nEdgesInLTSRegion(2, coarseRegions - iRegion + 1) + 1
                     edgesInLTSRegion(2, coarseRegions - iRegion + 1, &
                                      nEdgesInLTSRegion(2, coarseRegions - iRegion + 1)) &
                        = iEdge
                     isLTSRegionEdgeAssigned(iEdge) = 1
                  end if
               end if
            end do
         end do
      end do

      deallocate(isLTSRegionEdgeAssigned)

   end subroutine ocn_time_integration_fblts_init!}}}


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_lts_thick_tend
!
!> \brief  Calculate thickness tendencies for FB_LTS
!> \author Jeremy Lilly
!> \date   October 2023 
!> \details
!>  This routine calculates the thickness tendency on different 
!>  LTS regions
!
!-----------------------------------------------------------------------

   subroutine ocn_lts_thick_tend(LTSPool, tendPool, &
                                 normalVelocity, layerThickness, &
                                 computeOnFineInterior, computeOnFineInterfaceAdjacent, &
                                 computeOnInterfaceOne, computeOnInterfaceTwo, &
                                 computeOnCoarse)!{{{
      
      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------
      
      integer, intent(in) :: &
         computeOnFineInterior, &
         computeOnFineInterfaceAdjacent, &
         computeOnInterfaceOne, &
         computeOnInterfaceTwo, &
         computeOnCoarse

      real (kind=RKIND), dimension(:,:), pointer, intent(in) :: &
         layerThickness, &
         normalVelocity

      type (mpas_pool_type), intent(in) :: &
         LTSPool  !< Input: LTS data

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------
      
      type (mpas_pool_type), intent(inout) :: &
         tendPool  !< Input/output: tendency variables
    
      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------
      
      integer, dimension(:,:), pointer :: &
         nCellsInLTSRegion
      
      integer, dimension(:,:,:), pointer :: &
         cellsInLTSRegion
      
      real (kind=RKIND), dimension(:,:), pointer :: &
         layerThicknessTend

      integer :: &
         iEdge, cell1, cell2, k, ie, &
         ic, i, iCell, kmin, kmax
      real (kind=RKIND) :: &
         invdcEdge, flux 
      
      !-----------------------------------------------------------------
      ! Begin routine
      !-----------------------------------------------------------------
      
      call mpas_pool_get_array(LTSPool, 'cellsInLTSRegion', cellsInLTSRegion)
      call mpas_pool_get_array(LTSPool, 'nCellsInLTSRegion', nCellsInLTSRegion)

      call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

      ! calculate layerThickEdgeFlux
      if (config_thickness_flux_type == 'centered') then

         do iEdge = 1, nEdgesAll
            kmin = minLevelEdgeBot(iEdge)
            kmax = maxLevelEdgeTop(iEdge)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k = 1,nVertLevels
               ! initialize layerThicknessEdgeMean to avoid divide by
               ! zero and NaN problems.
               layerThickEdgeFlux(k,iEdge) = -1.0e34_RKIND
            end do
            do k = kmin,kmax
               ! central differenced
               layerThickEdgeFlux(k,iEdge) = 0.5_RKIND * &
                                             ( layerThickness(k,cell1) + &
                                             layerThickness(k,cell2) )
            end do
         end do

      else if (config_thickness_flux_type == 'upwind') then

         do iEdge = 1, nEdgesAll
            kmin = minLevelEdgeBot(iEdge)
            kmax = maxLevelEdgeTop(iEdge)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            do k=1,nVertLevels
               ! initialize layerThicknessEdgeFlux to avoid divide by
               ! zero and NaN problems.
               layerThickEdgeFlux(k,iEdge) = -1.0e34_RKIND
            end do
            do k = kmin,kmax
               if (normalVelocity(k,iEdge) > 0.0_RKIND) then
                  layerThickEdgeFlux(k,iEdge) = layerThickness(k,cell1)
               elseif (normalVelocity(k,iEdge) < 0.0_RKIND) then
                  layerThickEdgeFlux(k,iEdge) = layerThickness(k,cell2)
               else
                  layerThickEdgeFlux(k,iEdge) = max(layerThickness(k,cell1), &
                                                    layerThickness(k,cell2))
               end if
            end do
         end do

      else

         call mpas_log_write('ERROR: config_thickness_flux_type selected not implemented for FB_LTS', &
                             messageType=MPAS_LOG_CRIT)
         call abort

      end if
   
      layerThicknessTend(:, :) = 0.0_RKIND

      ! Fine in the interior, away from interface 1
      if (computeOnFineInterior == 1) then
         do ic = 1, nCellsInLTSRegion(1,1)
            iCell = cellsInLTSRegion(1,1,ic)
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThickEdgeFlux(k, iEdge)
                  layerThicknessTend(k, iCell) = layerThicknessTend(k,iCell) &
                                                 + edgeSignOnCell(i, iCell) * flux * invAreaCell(iCell)
               end do
            end do
         end do
      end if

      ! Fine adjacent to interface 1
      if (computeOnFineInterfaceAdjacent == 1) then
         do ic = 1, nCellsInLTSRegion(1,3)
            iCell = cellsInLTSRegion(1,3,ic)
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThickEdgeFlux(k, iEdge)
                  layerThicknessTend(k, iCell) = layerThicknessTend(k,iCell) &
                                                 + edgeSignOnCell(i, iCell) * flux * invAreaCell(iCell)
               end do
            end do
         end do
      end if

      ! Interface 1
      if (computeOnInterfaceOne == 1) then
         do ic = 1, nCellsInLTSRegion(1,2)
            iCell = cellsInLTSRegion(1,2,ic)
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThickEdgeFlux(k, iEdge)
                  layerThicknessTend(k, iCell) = layerThicknessTend(k, iCell) & 
                                                 + edgeSignOnCell(i, iCell) * flux * invAreaCell(iCell)
               end do
            end do
         end do
      end if

      ! Interface 2
      if (computeOnInterfaceTwo == 1) then
         do ic = 1, nCellsInLTSRegion(2,2)
            iCell = cellsInLTSRegion(2,2,ic)
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThickEdgeFlux(k, iEdge)
                  layerThicknessTend(k, iCell) = layerThicknessTend(k, iCell) & 
                                                 + edgeSignOnCell(i, iCell) * flux * invAreaCell(iCell)
               end do
            end do
         end do
      end if

      ! Coarse
      if (computeOnCoarse == 1) then
         do ic = 1, nCellsInLTSRegion(2,1)
            iCell = cellsInLTSRegion(2,1,ic)
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThickEdgeFlux(k, iEdge)
                  layerThicknessTend(k, iCell) = layerThicknessTend(k,iCell) &
                                                 + edgeSignOnCell(i, iCell) * flux * invAreaCell(iCell)
               end do
            end do
         end do
      end if

   end subroutine ocn_lts_thick_tend!}}}


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_lts_vel_tend
!
!> \brief  Calculate velocity tendencies for FB_LTS
!> \author Jeremy Lilly
!> \date   October 2023 
!> \details
!>  This routine calculates the velocity tendency on different 
!>  LTS regions
!
!-----------------------------------------------------------------------

   subroutine ocn_lts_vel_tend(LTSPool, tendPool, &
                               normalVelocity, layerThickness, &
                               computeOnFineInterior, computeOnFineInterfaceAdjacent, &
                               computeOnInterfaceOne, computeOnInterfaceTwo, &
                               computeOnCoarse)!{{{
      
      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------
      
      integer, intent(in) :: &
         computeOnFineInterior, &
         computeOnFineInterfaceAdjacent, &
         computeOnInterfaceOne, &
         computeOnInterfaceTwo, &
         computeOnCoarse

      real (kind=RKIND), dimension(:,:), pointer, intent(in) :: &
         layerThickness, &
         normalVelocity

      type (mpas_pool_type), intent(in) :: &
         LTSPool  !< Input: LTS data

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------
      
      type (mpas_pool_type), intent(inout) :: &
         tendPool  !< Input/output: tendency variables
    
      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------
      
      integer, dimension(:,:), pointer :: &
         nEdgesInLTSRegion
      
      integer, dimension(:,:,:), pointer :: &
         edgesInLTSRegion
      
      real (kind=RKIND), dimension(nCellsAll) :: &
         ssh
      
      real (kind=RKIND), dimension(:,:), pointer :: &
         normalVelocityTend

      integer :: &
         iEdge, cell1, cell2, k, ie, &
         ic, i, iCell, kmin, kmax

      real (kind=RKIND) :: &
         invdcEdge, betaSelfAttrLoad, flux, ssh_sal_on, tidal_pot_for_on
      
      !-----------------------------------------------------------------
      ! Begin routine
      !-----------------------------------------------------------------
      
      betaSelfAttrLoad = config_self_attraction_and_loading_beta

      call mpas_pool_get_array(LTSPool, 'edgesInLTSRegion', edgesInLTSRegion)
      call mpas_pool_get_array(LTSPool, 'nEdgesInLTSRegion', nEdgesInLTSRegion)

      call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)

      if (config_use_self_attraction_loading) then
         ssh_sal_on = 1.0_RKIND
      else
         ssh_sal_on = 0.0_RKIND
      endif
      if (config_use_tidal_potential_forcing) then
         tidal_pot_for_on = 1.0_RKIND
      else
         tidal_pot_for_on = 0.0_RKIND
      end if
      
      ! ssh
      do iCell = 1, nCellsAll
         k = maxLevelCell(iCell)
         zTop(k:nVertLevels,iCell) = -bottomDepth(iCell) + layerThickness(k,iCell)
         do k = maxLevelCell(iCell)-1, minLevelCell(iCell), -1
            zTop(k,iCell) = zTop(k+1,iCell) + layerThickness(k,iCell)
         end do
         ! copy zTop(1,iCell) into sea-surface height array
         ssh(iCell) = zTop(minLevelCell(iCell),iCell)
      end do

      normalVelocityTend(:,:) = 0.0_RKIND

      ! Fine in the interior, away from interface 1
      if (computeOnFineInterior == 1) then
         do ie = 1, nEdgesInLTSRegion(1,1)
            iEdge = edgesInLTSRegion(1,1,ie)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)
            do k=kMin,kMax
               normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) &
                                             - edgeMask(k,iEdge) * invdcEdge &
                                             * ( gravity * ( (ssh(cell2) - ssh(cell1)) &
                                             - tidal_pot_for_on * (1.0_RKIND - ssh_sal_on) &
                                             * betaSelfAttrLoad * (ssh(cell2) - ssh(cell1)) ) )
            end do
         end do
      end if

      ! Fine adjacent to interface 1
      if (computeOnFineInterfaceAdjacent == 1) then
         do ie = 1, nEdgesInLTSRegion(1,3)
            iEdge = edgesInLTSRegion(1,3,ie)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)
            do k=kMin,kMax
               normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) &
                                             - edgeMask(k,iEdge) * invdcEdge &
                                             * ( gravity * ( (ssh(cell2) - ssh(cell1)) &
                                             - tidal_pot_for_on * (1.0_RKIND - ssh_sal_on) &
                                             * betaSelfAttrLoad * (ssh(cell2) - ssh(cell1)) ) )
            end do
         end do
      end if

      ! Interface 1
      if (computeOnInterfaceOne == 1) then
         do ie = 1, nEdgesInLTSRegion(1,2)
            iEdge = edgesInLTSRegion(1,2,ie)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)
            do k=kMin,kMax
               normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) &
                                             - edgeMask(k,iEdge) * invdcEdge &
                                             * ( gravity * ( (ssh(cell2) - ssh(cell1)) &
                                             - tidal_pot_for_on * (1.0_RKIND - ssh_sal_on) &
                                             * betaSelfAttrLoad * (ssh(cell2) - ssh(cell1)) ) )
            end do
         end do
      end if

      ! Interface 2
      if (computeOnInterfaceTwo == 1) then
         do ie = 1, nEdgesInLTSRegion(2,2)
            iEdge = edgesInLTSRegion(2,2,ie)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)
            do k=kMin,kMax
               normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) &
                                             - edgeMask(k,iEdge) * invdcEdge &
                                             * ( gravity * ( (ssh(cell2) - ssh(cell1)) &
                                             - tidal_pot_for_on * (1.0_RKIND - ssh_sal_on) &
                                             * betaSelfAttrLoad * (ssh(cell2) - ssh(cell1)) ) )
            end do
         end do
      end if

      ! Coarse
      if (computeOnCoarse == 1) then
         do ie = 1, nEdgesInLTSRegion(2,1)
            iEdge = edgesInLTSRegion(2,1,ie)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)
            do k=kMin,kMax
               normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) &
                                             - edgeMask(k,iEdge) * invdcEdge &
                                             * ( gravity * ( (ssh(cell2) - ssh(cell1)) &
                                             - tidal_pot_for_on * (1.0_RKIND - ssh_sal_on) & 
                                             * betaSelfAttrLoad * (ssh(cell2) - ssh(cell1)) ) )
            end do
         end do
      end if

   end subroutine ocn_lts_vel_tend!}}}


end module ocn_time_integration_fblts

