! Copyright (c) 2021,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_lts
!
!> \brief MPAS barotropic ocean LTS Time integration scheme
!> \author Giacomo Capodaglio
!> \date   November 2022
!> \details
!>  This module contains the LTS init routine and the LTS
!>  barotropic ocean time integration scheme with splitting
!>  on the fast and slow tendency terms.
!
!-----------------------------------------------------------------------

module ocn_time_integration_lts

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_constants
   use ocn_tendency
   use ocn_diagnostics
   use ocn_mesh
   use ocn_vmix
   use ocn_config
   use ocn_diagnostics_variables
   use ocn_equation_of_state
   use ocn_time_average_coupled
   use ocn_time_varying_forcing
   use ocn_wetting_drying

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
   public :: ocn_time_integrator_lts, &
             ocn_time_integration_lts_init

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integrator_lts
!
!> \brief MPAS barotropic ocean LTS time integration scheme
!> \author Giacomo Capodaglio
!> \date   November 2022
!> \details
!>  This routine integrates one timestep (dt) using an LTS time integrator
!>  with a splitting of the fast and slow tendency terms
!
!-----------------------------------------------------------------------

    subroutine ocn_time_integrator_lts(domain,dt)!{{{
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Advance model state forward in time by the specified time step using
    ! a local time stepping scheme with spltting of the fast and slow tendencies
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !-----------------------------------------------------------------
    ! Input variables
    !-----------------------------------------------------------------

    real (kind=RKIND), intent(in) :: &
       dt      !< [in] time step (sec) to move forward

    !-----------------------------------------------------------------
    ! Input/output variables
    !-----------------------------------------------------------------

    type (domain_type), intent(inout) :: &
       domain  !< [inout] model state to advance forward

    !-----------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------

    integer :: &
       iCell, iEdge, iRegion, k, ic, ie, im, & !iterators
       M, & ! M = dtCoarse / dtFine
       nRegions !number of interface regions (two)

    type (block_type), pointer :: &
       block ! structure with subdomain data

    type (mpas_pool_type), pointer :: &
       tendPool, & ! structure holding tendencies
       statePool, & ! structure holding state variables
       meshPool, & ! structure holding mesh variables
       verticalMeshPool, & ! structure holding mesh variables
       forcingPool, & ! structure holding forcing variables
       scratchPool, & ! structure holding temporary variables
       tracersPool  ! structure holding tracers variables

    ! LTS Pools
    type (mpas_pool_type), pointer :: &
       LTSPool, & ! structure holding LTS variables
       tendSlowPool, & ! structure holding the slow tendency variables
       tendSum1stPool, & ! structure holding one of the correction terms for the interface
       tendSum2ndPool, & ! structure holding one of the correction terms for the interface
       tendSum3rdPool, & ! structure holding one of the correction terms for the interface
       prevTendSlowPool, nextTendSlowPool, & ! structures containing intermediate data
       prevTendSum1stPool, nextTendSum1stPool, & ! structures containing intermediate data
       prevTendSum2ndPool, nextTendSum2ndPool, & ! structures containing intermediate data
       prevTendSum3rdPool, nextTendSum3rdPool    ! structures containing intermediate data

    ! Tend Array Pointers
    real (kind=RKIND), dimension(:,:), pointer ::  &
       normalVelocityTend, & ! normal velocity fast tendency
       layerThicknessTend, & ! layer thickness tendency
       normalVelocityTendSlow, & ! normal velocity slow tendency
       normalVelocityTendSum1st, & ! one of the normal velocity correction terms for the interface
       layerThicknessTendSum1st, & ! one of the layer thickness correction terms for the interface
       normalVelocityTendSum2nd, & ! one of the normal velocity correction terms for the interface
       layerThicknessTendSum2nd, & ! one of the layer thickness correction terms for the interface
       normalVelocityTendSum3rd, & ! one of the normal velocity correction terms for the interface
       layerThicknessTendSum3rd  ! one of the layer thickness correction terms for the interface

    ! State Array Pointers
    real (kind=RKIND), dimension(:,:), pointer :: &
       normalVelocityCur, & ! normal velocity at time n
       normalVelocityNew, & ! normal velocity at time n+1
       normalVelocityFirstStage, & ! normal velocity at first stage of LTS
       normalVelocitySecondStage, & ! normal velocity at second stage of LTS
       layerThicknessCur, & ! layer thickness at time n
       layerThicknessNew, & ! layer thickness at time n+1
       layerThicknessFirstStage, & ! layer thickness at first stage of LTS
       layerThicknessSecondStage ! layer thickness at second stage of LTS

    ! LTS Objects
    real (kind=RKIND) :: &
       dtFine, & ! fine dt, defined as dt / M
       alpha, alphaHat, beta, betaHat, gam, gamHat  ! interpolation coefficients for the predictor step of LTS
    integer, dimension(:,:), pointer :: &
       nCellsInLTSRegion, & ! number of cells in a given LTS region
       nEdgesInLTSRegion ! number of edges in a given LTS region
    integer, dimension(:,:,:), pointer :: &
       cellsInLTSRegion, & ! list of cells in a given LTS region
       edgesInLTSRegion ! list of edges in a given LTS region
    real (kind=RKIND) :: &
       weightOld, weightNew, weightTend ! coefficents for the time advancing
    real (kind=RKIND) :: &
       weightTendSum1st, weightTendSum2nd, weightTendSum3rd ! coefficients for the interface correction

    integer err 
    err = 0

    if (.not. config_disable_tr_all_tend) then

       call mpas_log_write("ERROR: tracers are not currently implemented for local time-stepping")
       call mpas_log_write("config_disable_tr_all_tend should be true in the namelist file")
       call abort

    end if


    call mpas_timer_start("lts time-step prep")

    M = config_dt_scaling_LTS
    nRegions = 2
    dtFine = dt / M

    weightOld = 0.75_RKIND
    weightNew = 0.25_RKIND
    weightTend = 0.25_RKIND

    weightTendSum1st = 1.0_RKIND / 6.0_RKIND
    weightTendSum2nd = 1.0_RKIND / 6.0_RKIND
    weightTendSum3rd = 2.0_RKIND / 3.0_RKIND

    block => domain % blocklist

    !--- Retrieve model state, pools
    !--- No longer support sub-blocks so this retrieves the only block

    call mpas_pool_get_subpool(block%structs, 'mesh', &
                                               meshPool)
    call mpas_pool_get_subpool(block%structs, 'verticalMesh', &
                                               verticalMeshPool)
    call mpas_pool_get_subpool(block%structs, 'state', &
                                               statePool)
    call mpas_pool_get_subpool(block%structs, 'forcing', &
                                               forcingPool)
    call mpas_pool_get_subpool(block%structs, 'LTS', &
                                               LTSPool)
    call mpas_pool_get_subpool(block%structs, 'tend', &
                                               tendPool)
    call mpas_pool_get_subpool(block%structs, 'scratch', &
                                               scratchPool)
    call mpas_pool_get_subpool(statePool, 'tracers', &
                                               tracersPool)

    !--- Retrieve state variables at necessary time levels

    call mpas_pool_get_array(statePool, 'normalVelocity', &
                                         normalVelocityCur, 1)
    call mpas_pool_get_array(statePool, 'normalVelocity', &
                                         normalVelocityNew, 2)
    call mpas_pool_get_array(statePool, 'normalVelocity', &
                                         normalVelocityFirstStage, 3)
    call mpas_pool_get_array(statePool, 'normalVelocity', &
                                         normalVelocitySecondStage, 4)

    call mpas_pool_get_array(statePool, 'layerThickness', &
                                         layerThicknessCur, 1)
    call mpas_pool_get_array(statePool, 'layerThickness', &
                                         layerThicknessNew, 2)
    call mpas_pool_get_array(statePool, 'layerThickness', &
                                         layerThicknessFirstStage, 3)
    call mpas_pool_get_array(statePool, 'layerThickness', &
                                         layerThicknessSecondStage, 4)

    !--- Retrieve tendency variables

    call mpas_pool_get_array(tendPool, 'normalVelocity', &
                                        normalVelocityTend)
    call mpas_pool_get_array(tendPool, 'layerThickness', &
                                        layerThicknessTend)

    !--- Retrieve LTS arrays
    call mpas_pool_get_array(LTSPool, 'cellsInLTSRegion', &
                                       cellsInLTSRegion)
    call mpas_pool_get_array(LTSPool, 'nCellsInLTSRegion', &
                                       nCellsInLTSRegion)
    call mpas_pool_get_array(LTSPool, 'edgesInLTSRegion', &
                                       edgesInLTSRegion)
    call mpas_pool_get_array(LTSPool, 'nEdgesInLTSRegion', &
                                       nEdgesInLTSRegion)

    !--- Update additional pools for LTS
    call mpas_pool_get_subpool(block % structs, 'tend_sum_1st', tendSum1stPool)
    call mpas_pool_get_subpool(block % structs, 'tend_sum_2nd', tendSum2ndPool)
    call mpas_pool_get_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)
    call mpas_pool_get_subpool(block % structs, 'tend_slow', tendSlowPool)


    call mpas_pool_copy_pool(tendPool, tendSum1stPool, 1)
    call mpas_pool_copy_pool(tendPool, tendSum2ndPool, 1)
    call mpas_pool_copy_pool(tendPool, tendSum3rdPool, 1)
    call mpas_pool_copy_pool(tendPool, tendSlowPool, 1)


    call mpas_pool_get_array(tendSlowPool, 'normalVelocity', normalVelocityTendSlow)

    call mpas_pool_get_array(tendSum1stPool, 'normalVelocity', normalVelocityTendSum1st)
    call mpas_pool_get_array(tendSum1stPool, 'layerThickness', layerThicknessTendSum1st)

    call mpas_pool_get_array(tendSum2ndPool, 'normalVelocity', normalVelocityTendSum2nd)
    call mpas_pool_get_array(tendSum2ndPool, 'layerThickness', layerThicknessTendSum2nd)

    call mpas_pool_get_array(tendSum3rdPool, 'normalVelocity', normalVelocityTendSum3rd)
    call mpas_pool_get_array(tendSum3rdPool, 'layerThickness', layerThicknessTendSum3rd)


    !--- Init variables
    do iEdge = 1, nEdgesAll
       do k = 1, maxLevelEdgeTop(iEdge)
          normalVelocityNew(k, iEdge) = normalVelocityCur(k, iEdge)
          normalVelocityFirstStage(k, iEdge) = normalVelocityCur(k, iEdge)
          normalVelocitySecondStage(k, iEdge) = normalVelocityCur(k, iEdge)
       end do
    end do

    do iCell = 1, nCellsAll
       do k = 1, maxLevelCell(iCell)
          layerThicknessNew(k, iCell) = layerThicknessCur(k, iCell)
          layerThicknessFirstStage(k, iCell) = layerThicknessCur(k, iCell)
          layerThicknessSecondStage(k, iCell) = layerThicknessCur(k, iCell)
       end do
    end do

    normalVelocityTendSum1st(:,:) = 0.0_RKIND
    layerThicknessTendSum1st(:,:) = 0.0_RKIND

    normalVelocityTendSum2nd(:,:) = 0.0_RKIND
    layerThicknessTendSum2nd(:,:) = 0.0_RKIND

    normalVelocityTendSum3rd(:,:) = 0.0_RKIND
    layerThicknessTendSum3rd(:,:) = 0.0_RKIND


    call mpas_timer_stop("lts time-step prep")

    call mpas_timer_start("lts main loop")

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BEGIN LTS SCHEME
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (config_use_wetting_drying) then
       wettingVelocityFactor(:, :) = 0.0_RKIND
    end if    

    ! -------------------------------- COMPUTE SLOW TEND FOR ALL LTS REGIONS  -------------------------------

    call mpas_timer_start("lts slow tend computation")

    ! SLOW TENDENCY COMPUTATION --- (this is only for the momentum equation)
    call ocn_tend_vel(domain,tendSlowPool, statePool, forcingPool, 1, domain % dminfo, dt)

    call mpas_timer_stop("lts slow tend computation")

    ! --------------------------------------  STEP  ------------------------------------------
    ! --- compute the first stage of LTS for interface layer 1, interface layer 2 and coarse

    ! --- compute tendencies for thickness and fast tendencies only for velocity
    call mpas_timer_start("lts compute tendencies")

    call ocn_lts_tends(statePool, LTSPool, tendPool, dt, 1, 0, 1, 1, 1)

    call mpas_timer_stop("lts compute tendencies")

    call mpas_timer_start("lts advance soln")

    ! --- advance solution for interface layer 1, interface layer 2 and coarse with first stage

    normalVelocityTend(:,:) = normalVelocityTend(:,:) + normalVelocityTendSlow(:,:)

    ! --- NORMAL VELOCITY

    ! --- interface layers
    do iRegion =1,nRegions
       do ie = 1, nEdgesInLTSRegion(iRegion,2)
          iEdge = edgesInLTSRegion(iRegion,2,ie)
          normalVelocityFirstStage(:,iEdge) =  ( normalVelocityCur(:,iEdge) + dt * normalVelocityTend(:,iEdge) ) &
                                               * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
       end do
    end do
    ! --- coarse
    do ie = 1, nEdgesInLTSRegion(2,1)
       iEdge = edgesInLTSRegion(2,1,ie)
       normalVelocityFirstStage(:,iEdge) =  ( normalVelocityCur(:,iEdge) + dt * normalVelocityTend(:,iEdge) ) &
                                            * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
    end do
    ! --- fine layers close to interface layers
    do ie = 1, nEdgesInLTSRegion(1,3)
       iEdge = edgesInLTSRegion(1,3,ie)
       normalVelocityFirstStage(:,iEdge) = ( normalVelocityCur(:,iEdge) + dt * normalVelocityTend(:,iEdge) ) &
                                           * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
    end do


    ! --- LAYER THICKNESS

    ! --- interface layers
    do iRegion =1,nRegions
       do ic = 1, nCellsInLTSRegion(iRegion,2)
          iCell = cellsInLTSRegion(iRegion,2,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessFirstStage(k,iCell) = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell)
          end do
       end do
    end do
    ! --- coarse
    do ic = 1, nCellsInLTSRegion(2,1)
       iCell = cellsInLTSRegion(2,1,ic)
       do k = 1, maxLevelCell(iCell)
          layerThicknessFirstStage(k,iCell) = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell)
       end do
    end do
    ! --- fine layers close to interface layers
    do ic = 1, nCellsInLTSRegion(1,3)
       iCell = cellsInLTSRegion(1,3,ic)
       do k = 1, maxLevelCell(iCell)
          layerThicknessFirstStage(k,iCell) = layerThicknessCur(k,iCell) + dt * layerThicknessTend(k,iCell)
       end do
    end do

    call mpas_timer_stop("lts advance soln")

    call mpas_timer_start("lts prognostic halo update")
    call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=3)
    call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=3)
    call mpas_timer_stop("lts prognostic halo update")

    ! --------------------------------------  STEP  ------------------------------------------
    ! --- compute the second stage of LTS for interface layer 1, interface layer 2 and coarse

    call mpas_timer_start("lts compute tendencies")

    ! TENDENCIES COMPUTATION ---
    call ocn_lts_tends(statePool, LTSPool, tendPool, dt, 3, 0, 0, 1, 1)

    call mpas_timer_stop("lts compute tendencies")


    call mpas_timer_start("lts advance soln")
    ! --- advance solution for interface layer 1, interface layer 2 and coarse with second stage

    normalVelocityTend(:,:) = normalVelocityTend(:,:) + normalVelocityTendSlow(:,:)

    ! --- NORMAL VELOCITY

    ! --- coarse
    do ie = 1, nEdgesInLTSRegion(2,1)
       iEdge = edgesInLTSRegion(2,1,ie)
       normalVelocitySecondStage(:,iEdge) = ( weightOld * normalVelocityCur(:,iEdge) &
                                             + weightNew * normalVelocityFirstStage(:,iEdge) &
                                             + weightTend * dt * normalVelocityTend(:,iEdge) ) &
                                            * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
    end do

    ! --- interface layers
    do iRegion = 1, nRegions
       do ie = 1, nEdgesInLTSRegion(iRegion,2)
          iEdge = edgesInLTSRegion(iRegion,2,ie)
          normalVelocitySecondStage(:,iEdge) = ( weightOld * normalVelocityCur(:,iEdge) &
                                                + weightNew * normalVelocityFirstStage(:,iEdge) &
                                                + weightTend * dt * normalVelocityTend(:,iEdge) ) &
                                               * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
       end do
    end do

    ! --- LAYER THICKNESS

    ! --- coarse
    do ic = 1, nCellsInLTSRegion(2,1)
       iCell = cellsInLTSRegion(2,1,ic)
       do k = 1, maxLevelCell(iCell)
          layerThicknessSecondStage(k,iCell) = weightOld * layerThicknessCur(k,iCell) &
                                               + weightNew * layerThicknessFirstStage(k,iCell) &
                                               + weightTend * dt * layerThicknessTend(k,iCell)
       end do
    end do

    ! --- interface layers
    do iRegion = 1, nRegions
       do ic = 1, nCellsInLTSRegion(iRegion,2)
          iCell = cellsInLTSRegion(iRegion,2,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessSecondStage(k,iCell) = weightOld * layerThicknessCur(k,iCell) &
                                                  + weightNew * layerThicknessFirstStage(k,iCell) &
                                                  + weightTend * dt * layerThicknessTend(k,iCell)
          end do
       end do
    end do

    call mpas_timer_stop("lts advance soln")


    call mpas_timer_start("lts prognostic halo update")

    call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=4)
    call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=4)

    call mpas_timer_stop("lts prognostic halo update")


    ! --------------------------------------  STEP  ------------------------------------------
    ! --- substepping for fine soln advancement begins here. we also compute
    ! the correction terms for the interface to be used after the substepping

    do im =1, M

       call mpas_timer_start("lts interface prediction")
       ! --- interpolate interface layer 1 values (this is the actual interface prediction) and store soln cur in soln new

       alpha = REAL((im - 1)) / M
       alphaHat = REAL((im - 1) * (im - 1)) / (M * M)
       do ie = 1, nEdgesInLTSRegion(1,2) !interface layer 1
          iEdge = edgesInLTSRegion(1,2,ie)
          normalVelocityNew(:,iEdge) = normalVelocityCur(:,iEdge)
          normalVelocityCur(:,iEdge) = (1.0_RKIND - alpha - alphaHat) * normalVelocityCur(:,iEdge) + (alpha - alphaHat) * normalVelocityFirstStage(:,iEdge) &
                                        + 2.0_RKIND * alphaHat *normalVelocitySecondStage(:,iEdge)
       end do
       do ic = 1, nCellsInLTSRegion(1,2) !interface layer 1
          iCell = cellsInLTSRegion(1,2,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
             layerThicknessCur(k,iCell) = (1.0_RKIND - alpha - alphaHat) * layerThicknessCur(k,iCell) + (alpha - alphaHat) * layerThicknessFirstStage(k,iCell) &
                                           + 2.0_RKIND * alphaHat * layerThicknessSecondStage(k,iCell)
          end do
       end do

       call mpas_timer_stop("lts interface prediction")

       call mpas_timer_start("lts prognostic halo update")

       call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=1)
       call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=1)

       call mpas_timer_stop("lts prognostic halo update")

       call mpas_timer_start("lts compute tendencies")

       ! --- compute tendencies for thickness and fast tendencies only for velocity
       call ocn_lts_tends(statePool, LTSPool, tendPool, dtFine, 1, 1, 1, 0, 1)

       call mpas_timer_stop("lts compute tendencies")

       call mpas_timer_start("lts advance soln")

       ! --- advance fine solution (first stage of LTS) AND sum up the tendecies above into tendSum2nd for interface layer 1 and interface
       ! layer 2 also restore soln cur by copying it back from soln new

       normalVelocityTend(:,:) = normalVelocityTend(:,:) + normalVelocityTendSlow(:,:)

       ! --- NORMAL VELOCITY

       ! --- interface layer 1 (correction)
       do ie = 1, nEdgesInLTSRegion(1,2)
          iEdge = edgesInLTSRegion(1,2,ie)
          normalVelocityCur(:,iEdge) = normalVelocityNew(:,iEdge)
          normalVelocityTendSum2nd(:,iEdge) = ( 1.0_RKIND / (3.0_RKIND * dt * weightTendSum2nd) * normalVelocityCur(:,iEdge) &
                                                + normalVelocityTendSum2nd(:,iEdge) + normalVelocityTend(:,iEdge) ) &
                                              * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
       end do
       ! --- interface layer 2 (correction)
       do ie = 1, nEdgesInLTSRegion(2,2)
          iEdge = edgesInLTSRegion(2,2,ie)
          normalVelocityTendSum2nd(:,iEdge) = ( 1.0_RKIND / (3.0_RKIND * dt * weightTendSum2nd) * normalVelocityCur(:,iEdge) &
                                                + normalVelocityTendSum2nd(:,iEdge) + normalVelocityTend(:,iEdge) ) &
                                              * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
       end do

       ! --- fine (soln advancement)
       do ie = 1, nEdgesInLTSRegion(1,1)
          iEdge = edgesInLTSRegion(1,1,ie)
          normalVelocityFirstStage(:,iEdge) = ( normalVelocityCur(:,iEdge) + dtFine * normalVelocityTend(:,iEdge) ) & !soln update for the fine
                                              * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
       end do
       do ie = 1, nEdgesInLTSRegion(1,3)
          iEdge = edgesInLTSRegion(1,3,ie)
          normalVelocityFirstStage(:,iEdge) = ( normalVelocityCur(:,iEdge) + dtFine * normalVelocityTend(:,iEdge) ) &!soln update for the fine
                                              * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
       end do

       ! --- LAYER THICKNESS

       ! --- interface layer 1 (correction)
       do ic = 1, nCellsInLTSRegion(1,2)
          iCell = cellsInLTSRegion(1,2,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessCur(k,iCell) = layerThicknessNew(k,iCell)
              layerThicknessTendSum2nd(k,iCell) = layerThicknessTendSum2nd(k,iCell) +  layerThicknessTend(k,iCell)
          end do
       end do

       ! --- interface layer 2 (correction)
       do ic = 1, nCellsInLTSRegion(2,2)
          iCell = cellsInLTSRegion(2,2,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessTendSum2nd(k,iCell) = layerThicknessTendSum2nd(k,iCell) +  layerThicknessTend(k,iCell)
          end do
       end do

       ! --- fine (soln advancement)
       do ic = 1, nCellsInLTSRegion(1,1)
          iCell = cellsInLTSRegion(1,1,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessFirstStage(k,iCell) = layerThicknessCur(k,iCell) + dtFine * layerThicknessTend(k,iCell) !soln update for the fine
          end do
       end do
       do ic = 1, nCellsInLTSRegion(1,3)
          iCell = cellsInLTSRegion(1,3,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessFirstStage(k,iCell) = layerThicknessCur(k,iCell) + dtFine * layerThicknessTend(k,iCell) !soln update for the fine
          end do
       end do

       call mpas_timer_stop("lts advance soln")


       call mpas_timer_start("lts interface prediction")
       ! --- interpolate interface layer 1 values (this is the actual interface prediction) and store soln cur in soln new

       beta = REAL(im) / M
       betaHat = REAL((im - 1) * (im + 1)) / (M * M)
       do ie = 1, nEdgesInLTSRegion(1,2) !interface layer 1
          iEdge = edgesInLTSRegion(1,2,ie)
          normalVelocityNew(:,iEdge) = normalVelocityFirstStage(:,iEdge)
          normalVelocityFirstStage(:,iEdge) = (1.0_RKIND - beta - betaHat) * normalVelocityCur(:,iEdge) + (beta - betaHat) * normalVelocityFirstStage(:,iEdge) &
                                               + 2.0_RKIND * betaHat *normalVelocitySecondStage(:,iEdge)
       end do
       do ic = 1, nCellsInLTSRegion(1,2) !interface layer 1
          iCell = cellsInLTSRegion(1,2,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessNew(k,iCell) = layerThicknessFirstStage(k,iCell)
             layerThicknessFirstStage(k,iCell) = (1.0_RKIND - beta - betaHat) * layerThicknessCur(k,iCell) + (beta - betaHat) * layerThicknessFirstStage(k,iCell) &
                                                  + 2.0_RKIND * betaHat * layerThicknessSecondStage(k,iCell)
          end do
       end do

       call mpas_timer_stop("lts interface prediction")


       call mpas_timer_start("lts prognostic halo update")
       call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=3)
       call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=3)
       call mpas_timer_stop("lts prognostic halo update")

       call mpas_timer_start("lts compute tendencies")

       ! TENDENCIES COMPUTATION ---
       call ocn_lts_tends(statePool, LTSPool, tendPool, dtFine, 3, 1, 1, 0, 1)

       call mpas_timer_stop("lts compute tendencies")


       call mpas_timer_start("lts advance soln")
       ! --- advance fine solution (second stage of LTS) AND sum up the tendecies above into tendSum1st for interface layer 1 and interface
       ! layer 2 also restore soln first stage by copying it back from soln new

       normalVelocityTend(:,:) = normalVelocityTend(:,:) + normalVelocityTendSlow(:,:)

       ! --- NORMAL VELOCITY

       ! --- interface layer 1 (correction)
       do ie = 1, nEdgesInLTSRegion(1,2)
          iEdge = edgesInLTSRegion(1,2,ie)
          normalVelocityFirstStage(:,iEdge) = normalVelocityNew(:,iEdge)
          normalVelocityTendSum1st(:,iEdge) = ( 1.0_RKIND / (3.0_RKIND * dt * weightTendSum1st) * normalVelocityCur(:,iEdge) & 
                                                + normalVelocityTendSum1st(:,iEdge) + normalVelocityTend(:,iEdge) ) &
                                              * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
       end do
       ! --- interface layer 2 (correction)
       do ie = 1, nEdgesInLTSRegion(2,2)
          iEdge = edgesInLTSRegion(2,2,ie)
          normalVelocityTendSum1st(:,iEdge) = ( 1.0_RKIND / (3.0_RKIND * dt * weightTendSum1st) * normalVelocityCur(:,iEdge) & 
                                                + normalVelocityTendSum1st(:,iEdge) + normalVelocityTend(:,iEdge) ) &
                                              * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
       end do

       ! --- fine (soln advancement)
       do ie = 1, nEdgesInLTSRegion(1,1)
          iEdge = edgesInLTSRegion(1,1,ie)
          normalVelocitySecondStage(:,iEdge) = ( weightOld * normalVelocityCur(:,iEdge) + weightNew * normalVelocityFirstStage(:,iEdge) &
                                                + weightTend *  dtFine * normalVelocityTend(:,iEdge) ) & !soln update for the fine
                                               * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
       end do
       do ie = 1, nEdgesInLTSRegion(1,3)
          iEdge = edgesInLTSRegion(1,3,ie)
          normalVelocitySecondStage(:,iEdge) = ( weightOld * normalVelocityCur(:,iEdge) + weightNew * normalVelocityFirstStage(:,iEdge) &
                                                + weightTend *  dtFine * normalVelocityTend(:,iEdge) ) & !soln update for the fine
                                               * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
       end do

       ! --- LAYER THICKNESS

       ! --- interface layer 1 (correction)
       do ic = 1, nCellsInLTSRegion(1,2)
          iCell = cellsInLTSRegion(1,2,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessFirstStage(k,iCell) = layerThicknessNew(k,iCell)
             layerThicknessTendSum1st(k,iCell) = layerThicknessTendSum1st(k,iCell) +  layerThicknessTend(k,iCell)
          end do
       end do

       ! --- interface layer 2 (correction)
       do ic = 1, nCellsInLTSRegion(2,2)
          iCell = cellsInLTSRegion(2,2,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessTendSum1st(k,iCell) = layerThicknessTendSum1st(k,iCell) +  layerThicknessTend(k,iCell)
          end do
       end do

       ! --- fine (soln advancement)
       do ic = 1, nCellsInLTSRegion(1,1)
          iCell = cellsInLTSRegion(1,1,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessSecondStage(k,iCell) = weightOld * layerThicknessCur(k,iCell) + weightNew * layerThicknessFirstStage(k,iCell) &
                                                  + weightTend * dtFine * layerThicknessTend(k,iCell) !soln update for the fine
          end do
       end do
       do ic = 1, nCellsInLTSRegion(1,3)
          iCell = cellsInLTSRegion(1,3,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessSecondStage(k,iCell) = weightOld * layerThicknessCur(k,iCell) + weightNew * layerThicknessFirstStage(k,iCell) &
                                                  + weightTend * dtFine * layerThicknessTend(k,iCell) !soln update for the fine
          end do
       end do

       call mpas_timer_stop("lts advance soln")


       call mpas_timer_start("lts interface prediction")
       ! --- interpolate interface layer 1 values (this is the actual interface prediction) and store soln cur in soln new

       gam = REAL(2 * (im - 1) + 1) / ( 2 * M )
       gamHat = REAL(2 * (im - 1) * (im - 1) + 2 * (im - 1) + 1) / (2 * M * M)
       do ie = 1, nEdgesInLTSRegion(1,2) !interface layer 1
          iEdge = edgesInLTSRegion(1,2,ie)
          normalVelocityNew(:,iEdge) = normalVelocitySecondStage(:,iEdge)
          normalVelocitySecondStage(:,iEdge) = (1.0_RKIND - gam - gamHat) * normalVelocityCur(:,iEdge) + (gam - gamHat) * normalVelocityFirstStage(:,iEdge) &
                                                + 2.0_RKIND * gamHat *normalVelocitySecondStage(:,iEdge)
       end do
       do ic = 1, nCellsInLTSRegion(1,2) !interface layer 1
          iCell = cellsInLTSRegion(1,2,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessNew(k,iCell) = layerThicknessSecondStage(k,iCell)
             layerThicknessSecondStage(k,iCell) = (1.0_RKIND - gam - gamHat) * layerThicknessCur(k,iCell) + (gam - gamHat) * layerThicknessFirstStage(k,iCell) &
                                                   + 2.0_RKIND * gamHat * layerThicknessSecondStage(k,iCell)
         end do
        end do

       call mpas_timer_stop("lts interface prediction")

       call mpas_timer_start("lts prognostic halo update")

       call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=4)
       call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=4)

       call mpas_timer_stop("lts prognostic halo update")

       call mpas_timer_start("lts compute tendencies")

       call ocn_lts_tends(statePool, LTSPool, tendPool, dtFine, 4, 1, 1, 0, 1)

       call mpas_timer_stop("lts compute tendencies")


       call mpas_timer_start("lts advance soln")
       ! --- advance fine solution (third stage of LTS) AND sum up the tendecies above into tendSum3rd for interface layer 1 and interface
       ! layer 2 also restore soln second stage by copying it back from soln new

       normalVelocityTend(:,:) = normalVelocityTend(:,:) + normalVelocityTendSlow(:,:)

       ! --- NORMAL VELOCITY

       ! --- interface layer 1 (correction)
       do ie = 1, nEdgesInLTSRegion(1,2)
          iEdge = edgesInLTSRegion(1,2,ie)
          normalVelocitySecondStage(:,iEdge) = normalVelocityNew(:,iEdge)
          normalVelocityTendSum3rd(:,iEdge) = ( 1.0_RKIND / (3.0_RKIND * dt * weightTendSum3rd) * normalVelocityCur(:,iEdge) &
                                                + normalVelocityTendSum3rd(:,iEdge) + normalVelocityTend(:,iEdge) ) &
                                              * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
       end do
       ! --- interface layer 2 (correction)
       do ie = 1, nEdgesInLTSRegion(2,2)
          iEdge = edgesInLTSRegion(2,2,ie)
          normalVelocityTendSum3rd(:,iEdge) = ( 1.0_RKIND / (3.0_RKIND * dt * weightTendSum3rd) * normalVelocityCur(:,iEdge) &
                                                + normalVelocityTendSum3rd(:,iEdge) + normalVelocityTend(:,iEdge) ) &
                                              * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
       end do

       ! --- fine (soln advancement)
       do ie = 1, nEdgesInLTSRegion(1,1)
          iEdge = edgesInLTSRegion(1,1,ie)
          normalVelocityCur(:,iEdge) =  ( (1.0_RKIND / 3.0_RKIND) * normalVelocityCur(:,iEdge) + (2.0_RKIND / 3.0_RKIND) * normalVelocitySecondStage(:,iEdge) &
                                         + (2.0_RKIND / 3.0_RKIND) *  dtFine * normalVelocityTend(:,iEdge) ) & !soln update for the fine
                                        * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
       end do
       do ie = 1, nEdgesInLTSRegion(1,3)
          iEdge = edgesInLTSRegion(1,3,ie)
          normalVelocityCur(:,iEdge) =  ( (1.0_RKIND / 3.0_RKIND) * normalVelocityCur(:,iEdge) + (2.0_RKIND / 3.0_RKIND) * normalVelocitySecondStage(:,iEdge) &
                                         + (2.0_RKIND / 3.0_RKIND) *  dtFine * normalVelocityTend(:,iEdge) ) & !soln update for the fine
                                        * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
       end do

       ! --- LAYER THICKNESS

       ! --- interface layer 1 (correction)
       do ic = 1, nCellsInLTSRegion(1,2)
          iCell = cellsInLTSRegion(1,2,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessSecondStage(k,iCell) = layerThicknessNew(k,iCell)
             layerThicknessTendSum3rd(k,iCell) = layerThicknessTendSum3rd(k,iCell) +  layerThicknessTend(k,iCell)
          end do
       end do

       ! --- interface layer 2 (correction)
       do ic = 1, nCellsInLTSRegion(2,2)
          iCell = cellsInLTSRegion(2,2,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessTendSum3rd(k,iCell) = layerThicknessTendSum3rd(k,iCell) +  layerThicknessTend(k,iCell)
          end do
       end do

       ! --- fine (soln advancement)
       do ic = 1, nCellsInLTSRegion(1,1)
          iCell = cellsInLTSRegion(1,1,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessCur(k,iCell) = (1.0_RKIND / 3.0_RKIND)  * layerThicknessCur(k,iCell) + (2.0_RKIND / 3.0_RKIND) * layerThicknessSecondStage(k,iCell) &
                                           + (2.0_RKIND / 3.0_RKIND) * dtFine * layerThicknessTend(k,iCell) !soln update for the fine
          end do
       end do
       do ic = 1, nCellsInLTSRegion(1,3)
          iCell = cellsInLTSRegion(1,3,ic)
          do k = 1, maxLevelCell(iCell)
             layerThicknessCur(k,iCell) = (1.0_RKIND / 3.0_RKIND)  * layerThicknessCur(k,iCell) + (2.0_RKIND / 3.0_RKIND) * layerThicknessSecondStage(k,iCell) &
                                           + (2.0_RKIND / 3.0_RKIND) * dtFine * layerThicknessTend(k,iCell) !soln update for the fine
          end do
       end do

       call mpas_timer_stop("lts advance soln")

    end do ! --- end substepping


    call mpas_timer_start("lts copy soln into new soln")
    ! --- copy old soln into new soln for the fine cells and edges

    ! --- normal velocity
    do ie = 1, nEdgesInLTSRegion(1,1)
       iEdge = edgesInLTSRegion(1,1,ie)
       normalVelocityNew(:,iEdge) = normalVelocityCur(:,iEdge)
    end do
    do ie = 1, nEdgesInLTSRegion(1,3)
       iEdge = edgesInLTSRegion(1,3,ie)
       normalVelocityNew(:,iEdge) = normalVelocityCur(:,iEdge)
    end do

    ! --- layer thickness
    do ic = 1, nCellsInLTSRegion(1,1)
       iCell = cellsInLTSRegion(1,1,ic)
       do k = 1, maxLevelCell(iCell)
          layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
       end do
    end do
    do ic = 1, nCellsInLTSRegion(1,3)
       iCell = cellsInLTSRegion(1,3,ic)
       do k = 1, maxLevelCell(iCell)
          layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell)
       end do
    end do

    call mpas_timer_stop("lts copy soln into new soln")

    ! --------------------------------------  STEP  ------------------------------------------
    ! --- compute the third stage of LTS for coarse
    !

    call mpas_timer_start("lts compute tendencies")

    call ocn_lts_tends(statePool, LTSPool, tendPool, dt, 4, 0, 0, 1, 0)

    call mpas_timer_stop("lts compute tendencies")

    call mpas_timer_start("lts advance soln")
    ! --- advance solution for coarse with third stage of LTS

    normalVelocityTend(:,:) = normalVelocityTend(:,:) + normalVelocityTendSlow(:,:)

    ! --- NORMAL VELOCITY

    ! --- coarse
    do ie = 1, nEdgesInLTSRegion(2,1)
       iEdge = edgesInLTSRegion(2,1,ie)
       normalVelocityNew(:,iEdge) = ( (1.0_RKIND / 3.0_RKIND)  * normalVelocityCur(:,iEdge) + (2.0_RKIND / 3.0_RKIND) * normalVelocitySecondStage(:,iEdge) &
                                     + (2.0_RKIND / 3.0_RKIND) * dt * normalVelocityTend(:,iEdge) ) &
                                    * (1.0_RKIND - wettingVelocityFactor(:, iEdge))
    end do

    ! --- LAYER THICKNESS
    ! --- coarse
    do ic = 1, nCellsInLTSRegion(2,1)
       iCell = cellsInLTSRegion(2,1,ic)
       do k = 1, maxLevelCell(iCell)
          layerThicknessNew(k,iCell) = (1.0_RKIND / 3.0_RKIND)  * layerThicknessCur(k,iCell) + (2.0_RKIND / 3.0_RKIND) * layerThicknessSecondStage(k,iCell) &
                                        + (2.0_RKIND / 3.0_RKIND) * dt * layerThicknessTend(k,iCell)
       end do
    end do

    call mpas_timer_stop("lts advance soln")

    ! --------------------------------------  STEP  ------------------------------------------
    ! --- interface correction
    !

    call mpas_timer_start("lts interface correction")
    ! --- correct solution for interface layer 1 and interface layer 2

    do iRegion = 1,nRegions

       ! --- NORMAL VELOCITY

       ! --- interface layers

       do ie = 1, nEdgesInLTSRegion(iRegion,2)
          iEdge = edgesInLTSRegion(iRegion,2,ie)
          ! The corrections are computed like this (taking normalVelocityTendSum2nd as an example):
          ! normalVelocityTendSum2nd(:,iEdge) = ( 1 / (3*dt*weightTendSum2nd) * normalVelocityCur(:,iEdge) &
          ! + normalVelocityTendSum2nd(:,iEdge) + normalVelocityTend(:,iEdge) )
          normalVelocityNew(:,iEdge) =   weightTendSum2nd * dtFine * normalVelocityTendSum2nd(:,iEdge) &
                                         + weightTendSum1st * dtFine * normalVelocityTendSum1st(:,iEdge) &
                                         + weightTendSum3rd * dtFine * normalVelocityTendSum3rd(:,iEdge) 
       end do

       ! --- LAYER THICKNESS

       ! --- interface layers

       do ic = 1, nCellsInLTSRegion(iRegion,2)
          iCell = cellsInLTSRegion(iRegion,2,ic)
          do k= 1, maxLevelCell(iCell)
             layerThicknessNew(k,iCell) = layerThicknessCur(k,iCell) + weightTendSum2nd * dtFine * layerThicknessTendSum2nd(k,iCell) &
                                          + weightTendSum1st * dtFine * layerThicknessTendSum1st(k,iCell) &
                                          + weightTendsum3rd * dtFine * layerThicknessTendSum3rd(k,iCell)
          end do
       end do
    end do

    call mpas_timer_stop("lts interface correction")

    call mpas_timer_start("lts prognostic halo update")
    call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=2)
    call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=2)
    call mpas_timer_stop("lts prognostic halo update")

    call mpas_timer_start("lts implicit vmix")
    if (.not. config_disable_vel_vmix) then
       call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, verticalMeshPool, scratchPool, tracersPool, 2)
       call ocn_vmix_implicit(dt, meshPool, statePool, forcingPool, scratchPool, err, 2)
       call mpas_dmpar_field_halo_exch(domain, 'normalVelocity', timeLevel=2)
       call mpas_dmpar_field_halo_exch(domain, 'layerThickness', timeLevel=2)
    end if   
    call mpas_timer_stop("lts implicit vmix")

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! END LTS SCHEME
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call mpas_timer_stop("lts main loop")

    !
    !  A little clean up at the end
    !

    call mpas_timer_start("lts cleanup phase")

    ! verify that cells are not dry at conclusion of time step
    if (config_use_wetting_drying) then
       call mpas_timer_start("lts check wet cells")

       ! ensure existing layerThickness is valid
       if (config_verify_not_dry) then
          call ocn_wetting_drying_verify(block, config_drying_min_cell_height, err)
       end if

       call mpas_timer_stop("lts check wet cells")
    end if

    ! direct application of tidal boundary condition
    if (config_use_tidal_forcing .and. trim(config_tidal_forcing_type) == 'direct') then
       call  ocn_lts_apply_direct_tidal_bc(statePool, forcingPool, verticalMeshPool, 2)
    end if

    if (config_prescribe_velocity) then
       do iEdge = 1, nEdgesAll
          normalVelocityNew(:, iEdge) = normalVelocityCur(:, iEdge)
       end do
    end if

    if (config_prescribe_thickness) then
       do iCell = 1, nCellsAll
          layerThicknessNew(:, iCell) = layerThicknessCur(:, iCell)
       end do
    end if

    do iEdge = 1, nEdgesAll
       normalTransportVelocity(:, iEdge) = normalVelocityNew(:, iEdge)
    end do

    call mpas_reconstruct(meshPool,  normalVelocityNew, &
                          velocityX, velocityY, velocityZ,   &
                          velocityZonal, velocityMeridional, &
                          includeHalos = .true.)

    call mpas_reconstruct(meshPool, gradSSH,          &
                          gradSSHX, gradSSHY, gradSSHZ,    &
                          gradSSHZonal, gradSSHMeridional, &
                          includeHalos = .true.)

    do iCell = 1, nCellsAll
       surfaceVelocity(indexSurfaceVelocityZonal, iCell) = velocityZonal(1, iCell)
       surfaceVelocity(indexSurfaceVelocityMeridional, iCell) = velocityMeridional(1, iCell)

       SSHGradient(indexSSHGradientZonal, iCell) = gradSSHZonal(iCell)
       SSHGradient(indexSSHGradientMeridional, iCell) = gradSSHMeridional(iCell)
    end do

    call ocn_time_average_coupled_accumulate( statePool, forcingPool, 2)

    do iCell = 1, nCellsAll
       surfaceVelocity(indexSurfaceVelocityZonal, iCell) = velocityZonal(1, iCell)
       surfaceVelocity(indexSurfaceVelocityMeridional, iCell) = velocityMeridional(1, iCell)

       SSHGradient(indexSSHGradientZonal, iCell) = gradSSHZonal(iCell)
       SSHGradient(indexSSHGradientMeridional, iCell) = gradSSHMeridional(iCell)
    end do

    ! DIAGNOSTICS UPDATE ---
    call ocn_diagnostic_solve(dt, statePool, forcingPool, meshPool, verticalMeshPool, scratchPool, tracersPool, 2)

    call mpas_timer_stop("lts cleanup phase")


  end subroutine ocn_time_integrator_lts!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_lts_init
!
!> \brief MPAS ocean LTS3 Time integration scheme
!> \author Giacomo Capodaglio
!> \date   November 2022
!> \details
!>  This routine computes the LTS arrays
!
!-----------------------------------------------------------------------

    subroutine ocn_time_integration_lts_init(domain)!{{{
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Sets the LTS arrays
    !
    ! Output: LTS arrays are written
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none

    type (domain_type), intent(inout) :: domain

    type (block_type), pointer :: block
    type (mpas_pool_type), pointer :: LTSPool
    type (mpas_pool_type), pointer :: tendPool
    type (mpas_pool_type), pointer :: &
       tendSlowPool, &
       tendSum1stPool, &
       tendSum2ndPool, &
       tendSum3rdPool, &
       prevTendSlowPool, nextTendSlowPool, &
       prevTendSum1stPool, nextTendSum1stPool, &
       prevTendSum2ndPool, nextTendSum2ndPool, &
       prevTendSum3rdPool, nextTendSum3rdPool
    integer, dimension(:), allocatable :: isLTSRegionEdgeAssigned
    integer :: i, iCell, iEdge, iRegion, coarseRegions, fineRegions, fineRegionsM1
    integer, dimension(:), pointer :: LTSRegion
    integer, dimension(:,:), pointer :: nCellsInLTSRegion, nEdgesInLTSRegion
    integer, dimension(:,:,:), pointer :: cellsInLTSRegion, edgesInLTSRegion
    integer, dimension(2) :: minMaxLTSRegion

    minMaxLTSRegion(1) = 1
    minMaxLTSRegion(2) = 2

    block => domain % blocklist

    ! Create additional pools 
    call mpas_pool_get_subpool(block%structs, 'tend', tendPool)

    call mpas_pool_create_pool(tendSum1stPool)
    call mpas_pool_clone_pool(tendPool, tendSum1stPool, 1)
    call mpas_pool_create_pool(tendSum2ndPool)
    call mpas_pool_clone_pool(tendPool, tendSum2ndPool, 1)
    call mpas_pool_create_pool(tendSum3rdPool)
    call mpas_pool_clone_pool(tendPool, tendSum3rdPool, 1)
    call mpas_pool_create_pool(tendSlowPool)
    call mpas_pool_clone_pool(tendPool, tendSlowPool, 1)

    call mpas_pool_add_subpool(block % structs, 'tend_sum_1st', tendSum1stPool)
    call mpas_pool_add_subpool(block % structs, 'tend_sum_2nd', tendSum2ndPool)
    call mpas_pool_add_subpool(block % structs, 'tend_sum_3rd', tendSum3rdPool)
    call mpas_pool_add_subpool(block % structs, 'tend_slow', tendSlowPool)

    if (associated(block % prev)) then
       call mpas_pool_get_subpool(block % prev % structs, 'tend_sum_1st', tendSum1stPool)
       call mpas_pool_get_subpool(block % prev % structs, 'tend_sum_2nd', tendSum2ndPool)
       call mpas_pool_get_subpool(block % prev % structs, 'tend_sum_3rd', tendSum3rdPool)
       call mpas_pool_get_subpool(block % prev % structs, 'tend_slow', tendSlowPool)
    else
       nullify(prevTendSum1stPool)
       nullify(prevTendSum2ndPool)
       nullify(prevTendSum3rdPool)
       nullify(prevTendSlowPool)
    end if

    if (associated(block % next)) then
         call mpas_pool_get_subpool(block % next % structs, 'tend_sum_1st', nextTendSum1stPool)
         call mpas_pool_get_subpool(block % next % structs, 'tend_sum_2nd', nextTendSum2ndPool)
         call mpas_pool_get_subpool(block % next % structs, 'tend_sum_3rd', nextTendSum3rdPool)
         call mpas_pool_get_subpool(block % next % structs, 'tend_slow', nextTendSlowPool)
    else
       nullify(nextTendSum1stPool)
       nullify(nextTendSum2ndPool)
       nullify(nextTendSum3rdPool)
       nullify(nextTendSlowPool)
    end if

    if (associated(prevTendSum1stPool) .and. associated(nextTendSum1stPool)) then
       call mpas_pool_link_pools(tendSum1stPool, prevTendSum1stPool, nextTendSum1stPool)
    else if (associated(prevTendSum1stPool)) then
       call mpas_pool_link_pools(tendSum1stPool, prevTendSum1stPool)
    else if (associated(nextTendSum1stPool)) then
       call mpas_pool_link_pools(tendSum1stPool,nextPool=nextTendSum1stPool)
    else
       call mpas_pool_link_pools(tendSum1stPool)
    end if

    if (associated(prevTendSum2ndPool) .and. associated(nextTendSum2ndPool)) then
       call mpas_pool_link_pools(tendSum2ndPool, prevTendSum2ndPool, nextTendSum2ndPool)
    else if (associated(prevTendSum2ndPool)) then
       call mpas_pool_link_pools(tendSum2ndPool, prevTendSum2ndPool)
    else if (associated(nextTendSum2ndPool)) then
       call mpas_pool_link_pools(tendSum2ndPool,nextPool=nextTendSum2ndPool)
    else
       call mpas_pool_link_pools(tendSum2ndPool)
    end if

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

    call mpas_pool_link_parinfo(block, tendSum1stPool)
    call mpas_pool_link_parinfo(block, tendSum2ndPool)
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
    ! 5 = fine (to advance when doing 1st stage on interface)

    nCellsInLTSRegion(:,:) = 0
    nEdgesInLTSRegion(:,:) = 0

    ! this is a loop to build the lists of elements in the fine, coarse, and interface regions
    do iCell = 1, nCellsOwned !this loops up to nCellsOwned because in the time stepping we only want to advance the cells owned by the MPI process
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


    ! below we fill out the lists for the edges, according to the
    ! LTSRegion that have been assigned to the cells.
    ! we move from the fine to the coarse (i.e. from the fine to the nearest
    ! LTS region in the direction of the coarse).
    ! Note that edges shared between cells of different LTS regions are owned
    ! by the cell in the LTS region closest to the fine region,
    ! see Figure 3 in "Conservative explicit local time-stepping schemes for
    ! the shallow water equations" by Hoang et al. (halo edges however
    ! are owned by whatever process they are initially assigned to)

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
                   nEdgesInLTSRegion(1, fineRegions - iRegion + 1) = nEdgesInLTSRegion(1, fineRegions - iRegion + 1) + 1
                   edgesInLTSRegion(1, fineRegions - iRegion + 1, nEdgesInLTSRegion(1, fineRegions - iRegion + 1)) = iEdge
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
                   nEdgesInLTSRegion(2, coarseRegions - iRegion + 1) = nEdgesInLTSRegion(2, coarseRegions - iRegion + 1) + 1
                   edgesInLTSRegion(2, coarseRegions - iRegion + 1, nEdgesInLTSRegion(2, coarseRegions - iRegion + 1)) = iEdge
                   isLTSRegionEdgeAssigned(iEdge) = 1
                end if
             end if
          end do
       end do
    end do

    deallocate(isLTSRegionEdgeAssigned)

    end subroutine ocn_time_integration_lts_init!}}}

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine ocn_lts_tends(statePool, LTSPool, tendPool, dt, timeLevelIn, computeOnFineBig, computeOnFineSmall, computeOnCoarse, computeOnInterface)!{{{

   integer, intent(in) :: timeLevelIn, computeOnFineBig, computeOnFineSmall, computeOnCoarse, computeOnInterface

   type (mpas_pool_type), intent(in) :: &
      statePool          !< Input: state variables

   type (mpas_pool_type), intent(in) :: &
      LTSPool          !< Input: LTS data

    real (kind=RKIND), intent(in) :: dt

   type (mpas_pool_type), intent(inout) :: &
      tendPool          !< Input: tendency variables

   integer, dimension(:,:), pointer :: nCellsInLTSRegion, nEdgesInLTSRegion
   integer, dimension(:,:,:), pointer :: cellsInLTSRegion, edgesInLTSRegion
   real (kind=RKIND), dimension(:), pointer :: ssh
   real (kind=RKIND), dimension(:,:), pointer :: normalVelocityTend, layerThicknessTend, normalVelocity, layerThickness

   integer :: iEdge, cell1, cell2, k, ie, iRegion, nRegions, ic, i, iCell, nCells, kmin, kmax
   real (kind=RKIND) :: invdcEdge, betaSelfAttrLoad, flux, ssh_sal_on, tidal_pot_for_on

   betaSelfAttrLoad = config_self_attraction_and_loading_beta
   nRegions = 2

   call mpas_pool_get_array(LTSPool, 'cellsInLTSRegion',cellsInLTSRegion)
   call mpas_pool_get_array(LTSPool, 'nCellsInLTSRegion',nCellsInLTSRegion)
   call mpas_pool_get_array(LTSPool, 'edgesInLTSRegion', edgesInLTSRegion)
   call mpas_pool_get_array(LTSPool, 'nEdgesInLTSRegion',nEdgesInLTSRegion)

   call mpas_pool_get_array(tendPool, 'normalVelocity', normalVelocityTend)
   call mpas_pool_get_array(tendPool, 'layerThickness', layerThicknessTend)

   call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevelIn)
   call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevelIn)
   call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevelIn)

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

   ! inline computation of the diagnostics

   ! layerThickEdgeFlux
   if(config_thickness_flux_type == 'centered') then

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
                                         (layerThickness(k,cell1) + &
                                          layerThickness(k,cell2))
         end do
      end do

   else if(config_thickness_flux_type == 'upwind') then

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
               layerThickEdgeFlux(k,iEdge) = max(layerThickness(k,cell1), layerThickness(k,cell2))
            end if
         end do
      end do

   else

      call mpas_log_write('ERROR: config_thickness_flux_type selected not implemented for LTS', messageType=MPAS_LOG_CRIT)
      call abort

   end if

   ! ssh
   do iCell = 1, nCellsAll
      k = maxLevelCell(iCell)
      zTop(k:nVertLevels,iCell) = -bottomDepth(iCell) + layerThickness(k,iCell)
      do k = maxLevelCell(iCell)-1, minLevelCell(iCell), -1
         zTop(k,iCell) = zTop(k+1,iCell) + layerThickness(k  ,iCell)
      end do
      ! copy zTop(1,iCell) into sea-surface height array
      ssh(iCell) = zTop(minLevelCell(iCell),iCell)
   end do

   if (config_use_wetting_drying) then
      call ocn_lts_wd(statePool, dt, timeLevelIn)
   end if

   ! inline computation of the tendencies
   normalVelocityTend(:,:) = 0.0_RKIND
   layerThicknessTend(:,:) = 0.0_RKIND

   ! interface
   if (computeOnInterface == 1) then
      do iRegion = 1, nRegions
         ! velocity tendency
         do ie = 1, nEdgesInLTSRegion(iRegion,2)
            iEdge = edgesInLTSRegion(iRegion,2,ie)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            invdcEdge = 1.0_RKIND / dcEdge(iEdge)
            kMin = minLevelEdgeBot(iEdge)
            kMax = maxLevelEdgeTop(iEdge)
            do k=kMin,kMax
               normalVelocityTend(k,iEdge) =  normalVelocityTend(k,iEdge) &
                                              - edgeMask(k,iEdge) * invdcEdge * ( gravity * ( (ssh(cell2) - ssh(cell1)) &
                                                - tidal_pot_for_on * (1.0_RKIND - ssh_sal_on) * betaSelfAttrLoad * (ssh(cell2) - ssh(cell1)) ) )
            end do
         end do
         ! thickness tendency
         do ic = 1, nCellsInLTSRegion(iRegion,2)
            iCell = cellsInLTSRegion(iRegion,2,ic)
            do i = 1, nEdgesOnCell(iCell)
               iEdge = edgesOnCell(i, iCell)
               do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                  flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThickEdgeFlux(k, iEdge) * (1.0_RKIND - wettingVelocityFactor(k, iEdge))
                  layerThicknessTend(k, iCell) = layerThicknessTend(k, iCell) + edgeSignOnCell(i, iCell) * flux * invAreaCell(iCell)
               end do
            end do
         end do
      end do
   end if

   if (computeOnFineBig == 1) then
      ! velocity tendency
      do ie = 1, nEdgesInLTSRegion(1,1)
         iEdge = edgesInLTSRegion(1,1,ie)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         invdcEdge = 1.0_RKIND / dcEdge(iEdge)
         kMin = minLevelEdgeBot(iEdge)
         kMax = maxLevelEdgeTop(iEdge)
         do k=kMin,kMax
            normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) &
                                          - edgeMask(k,iEdge) * invdcEdge * ( gravity * ( (ssh(cell2) - ssh(cell1))  &
                                            - tidal_pot_for_on * (1.0_RKIND - ssh_sal_on) * betaSelfAttrLoad * (ssh(cell2) - ssh(cell1)) ) )
         end do
      end do
      ! thickness tendency
      do ic = 1, nCellsInLTSRegion(1,1)
         iCell = cellsInLTSRegion(1,1,ic)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
               flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThickEdgeFlux(k, iEdge) * (1.0_RKIND - wettingVelocityFactor(k, iEdge))
               layerThicknessTend(k, iCell) = layerThicknessTend(k, iCell) + edgeSignOnCell(i, iCell) * flux * invAreaCell(iCell)
            end do
         end do
      end do
   end if

   if (computeOnFineSmall == 1) then
      ! velocity tendency
      do ie = 1, nEdgesInLTSRegion(1,3)
         iEdge = edgesInLTSRegion(1,3,ie)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         invdcEdge = 1.0_RKIND / dcEdge(iEdge)
         kMin = minLevelEdgeBot(iEdge)
         kMax = maxLevelEdgeTop(iEdge)
         do k=kMin,kMax
            normalVelocityTend(k,iEdge) = normalVelocityTend(k,iEdge) &
                                          - edgeMask(k,iEdge) * invdcEdge * ( gravity * ( (ssh(cell2) - ssh(cell1))  &
                                            - tidal_pot_for_on * (1.0_RKIND - ssh_sal_on) * betaSelfAttrLoad * (ssh(cell2) - ssh(cell1)) ) )
         end do
      end do
      ! thickness tendency
      do ic = 1, nCellsInLTSRegion(1,3)
         iCell = cellsInLTSRegion(1,3,ic)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
               flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThickEdgeFlux(k, iEdge) * (1.0_RKIND - wettingVelocityFactor(k, iEdge))
               layerThicknessTend(k, iCell) = layerThicknessTend(k, iCell) + edgeSignOnCell(i, iCell) * flux * invAreaCell(iCell)
            end do
         end do
      end do
   end if

   if (computeOnCoarse == 1) then
      ! velocity tendency
      do ie = 1, nEdgesInLTSRegion(2,1)
         iEdge = edgesInLTSRegion(2,1,ie)
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         invdcEdge = 1.0_RKIND / dcEdge(iEdge)
         kMin = minLevelEdgeBot(iEdge)
         kMax = maxLevelEdgeTop(iEdge)
         do k=kMin,kMax
            normalVelocityTend(k,iEdge) =  normalVelocityTend(k,iEdge) &
                                           - edgeMask(k,iEdge) * invdcEdge * ( gravity * ( (ssh(cell2) - ssh(cell1))  &
                                             - tidal_pot_for_on * (1.0_RKIND - ssh_sal_on) * betaSelfAttrLoad * (ssh(cell2) - ssh(cell1)) ) )
         end do
      end do
      ! thickness tendency
      do ic = 1, nCellsInLTSRegion(2,1)
         iCell = cellsInLTSRegion(2,1,ic)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            do k = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
               flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThickEdgeFlux(k, iEdge) * (1.0_RKIND - wettingVelocityFactor(k, iEdge))
               layerThicknessTend(k, iCell) = layerThicknessTend(k, iCell) + edgeSignOnCell(i, iCell) * flux * invAreaCell(iCell)
            end do
         end do
      end do
   end if

   end subroutine ocn_lts_tends!}}}

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine ocn_lts_wd(statePool, dt, timeLevelIn)

   type (mpas_pool_type), intent(in) :: &
      statePool          !< Input: state variables

   real (kind=RKIND), intent(in) :: dt

   integer, intent(in) :: timeLevelIn

   real (kind=RKIND), dimension(:, :), pointer :: layerThickness
   real (kind=RKIND), dimension(:, :), pointer :: normalVelocity
   real (kind=RKIND) :: divOutFlux, layerThicknessTmp

   integer :: iEdge, k, iCell, i

   call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevelIn)
   call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevelIn)

   do iCell = 1, nCellsAll
      do k = minLevelCell(iCell), maxLevelCell(iCell)
         divOutFlux = 0.0_RKIND
         layerThicknessTmp = layerThickness(k, iCell)
         do i = 1, nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i, iCell)
            if (k <= maxLevelEdgeTop(iEdge) .and. k >= minLevelEdgeBot(iEdge)) then
               ! only consider divergence flux leaving the cell
               if ( normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) < 0.0_RKIND ) then
                  divOutFlux = divOutFlux + &
                               normalVelocity(k, iEdge) * edgeSignOnCell(i, iCell) * &
                               layerThickEdgeFlux(k, iEdge) * dvEdge(iEdge) * &
                               invAreaCell(iCell)
                end if
             end if
          end do
          layerThicknessTmp = layerThicknessTmp + dt * divOutFlux

          call ocn_wetting_velocity_factor_on_cell_edges(wettingVelocityFactor, layerThicknessTmp, normalVelocity, iCell, k)

      end do
   end do

   call mpas_dmpar_exch_halo_field(wettingVelocityField)

   end subroutine ocn_lts_wd!}}}   

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine ocn_lts_apply_direct_tidal_bc(statePool, forcingPool, verticalMeshPool, timeLevelIn)!{{{

   type (mpas_pool_type), intent(inout) :: statePool

   type (mpas_pool_type), intent(in) :: forcingPool, verticalMeshPool

   integer, intent(in) :: timeLevelIn

   integer :: iCell, k

   real (kind=RKIND) :: totalDepth
   real (kind=RKIND), dimension(:,:), pointer :: layerThicknessNew, restingThickness
   real (kind=RKIND), dimension(:), pointer :: tidalInputMask, tidalBCValue

   call mpas_pool_get_array(statePool, 'layerThickness', layerThicknessNew, timeLevelIn)
   call mpas_pool_get_array(forcingPool, 'tidalInputMask', tidalInputMask)
   call mpas_pool_get_array(forcingPool, 'tidalBCValue', tidalBCValue)
   call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)

   do iCell=1, nCellsOwned
      if (tidalInputMask(iCell) == 1.0_RKIND) then
         ! compute total depth for relative thickness contribution
         totalDepth = 0.0_RKIND
         do k = minLevelCell(iCell), maxLevelCell(iCell)
            totalDepth = totalDepth + restingThickness(k,iCell)
         end do
         ! only modify layer thicknesses on tidal boundary
         do k = minLevelCell(iCell), maxLevelCell(iCell)
            layerThicknessNew(k, iCell) = tidalInputMask(iCell)*(tidalBCValue(iCell) + bottomDepth(iCell))*(restingThickness(k,iCell)/totalDepth)
         end do
      end if
   end do

   end subroutine ocn_lts_apply_direct_tidal_bc!}}}    

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   real function sphere_distance(lat1, lon1, lat2, lon2, radius)

   ! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on
   ! a sphere with given radius.

   implicit none

   real (kind=RKIND), intent(in) :: lat1, lon1, lat2, lon2, radius

   real (kind=RKIND) :: arg1

   arg1 = sqrt( sin(0.5*(lat2-lat1))**2 +  &
                cos(lat1)*cos(lat2)*sin(0.5*(lon2-lon1))**2 )
   sphere_distance = 2.*radius*asin(arg1)

   end function sphere_distance

end module ocn_time_integration_lts
