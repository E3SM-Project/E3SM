!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vmix_cvmix
!
!> \brief MPAS ocean vertical mixing interface to CVMix
!> \author Todd Ringler
!> \date   04 February 2013
!> \details
!>  This module contains the routines for calls into CVMix
!>
!
!-----------------------------------------------------------------------

module ocn_vmix_cvmix

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_constants
   use mpas_log

   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables
   use ocn_stokes_drift

   use cvmix_kinds_and_types
   use cvmix_put_get
   use cvmix_background
   use cvmix_ddiff
   use cvmix_convection
   use cvmix_shear
   use cvmix_tidal
   use cvmix_kpp

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

   public :: ocn_vmix_coefs_cvmix_build, &
             ocn_vmix_cvmix_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   type(cvmix_global_params_type) :: cvmix_global_params
   type(cvmix_bkgnd_params_type)  :: cvmix_background_params
   type(cvmix_shear_params_type)  :: cvmix_shear_params
   type(cvmix_tidal_params_type)  :: cvmix_tidal_params

   logical :: cvmixOn, cvmixConvectionOn, cvmixKPPOn
   real (kind=RKIND) :: backgroundVisc, backgroundDiff

   integer :: cvmixBackgroundChoice        ! user choice of cvmix background scheme

   integer, parameter ::                 & ! supported cvmix background schemes
      cvmixBackgroundTypeNone = 0,       & ! no background scheme
      cvmixBackgroundTypeConstant = 1,   & ! Constant background scheme
      cvmixBackgroundTypeBryanLewis = 2    ! Bryan Lewis background scheme

   integer :: cvmixShearMixingChoice       ! user choice of vmix shear mixing scheme

   integer, parameter ::                 & ! supported cvmix shear mixing schemes
      cvmixShearMixingTypePP = 1,        & ! PP mixing scheme
      cvmixShearMixingTypeKPP = 2          ! KPP mixing scheme

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_vmix_coefs_cmvix_build
!
!> \brief   Computes mixing coefficients using CVMix
!> \author  Todd Ringler
!> \date    04 February 2013
!> \details
!>  This routine computes the vertical mixing coefficients for momentum
!>  and tracers by calling CVMix routines.
!
!-----------------------------------------------------------------------

   subroutine ocn_vmix_coefs_cvmix_build(meshPool, statePool, forcingPool, err, timeLevelIn)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      integer, intent(in), optional :: timeLevelIn !< Input: time level for state pool

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         statePool         !< Input/Output: state information

      type (mpas_pool_type), intent(inout) :: &
         forcingPool   !< Input/Output: forcing information

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

      type(cvmix_data_type) :: cvmix_variables

      integer, dimension(:), pointer :: &
        maxLevelCell, minLevelCell, nEdgesOnCell, maxLevelEdgeTop, minLevelEdgeBot

      integer, dimension(:), pointer :: landIceMask
      integer :: landIceMaskValue

      real (kind=RKIND), dimension(:), pointer :: &
        latCell, lonCell, bottomDepth, fCell, &
        ssh, dcEdge, dvEdge, areaCell, iceFraction, &
        windSpeed10m, langmuirNumber

      real (kind=RKIND), dimension(:,:), pointer :: &
        layerThickness, &
        normalVelocity

      integer, dimension(:,:), pointer :: edgesOnCell, cellsOnEdge, cellsOnCell, cellMask

      integer :: k, i, iCell, iNeighbor, timeLevel, kIndexOBL, kav, iEdge, nCells
      integer :: edgeCount, nEdges, topIndex, nsmooth, kpp_stage
      integer, pointer :: nVertLevels, nVertLevelsP1
      integer, dimension(:), pointer :: nCellsArray
      integer, dimension(:), allocatable :: surfaceAverageIndex

      real (kind=RKIND) :: x, bulkRichardsonNumberStop, sfc_layer_depth
      real (kind=RKIND) :: normalVelocityAv, delU2, areaSum, blTemp,tangentialVelocityAv
      real (kind=RKIND) :: sigma
      real (kind=RKIND), dimension(:), allocatable :: Nsqr_iface, turbulentScalarVelocityScale, &
                                                      deltaVelocitySquared, normalVelocitySum, &
                                                      potentialDensitySum, RiTemp,              &
                                                      layerThicknessSum, layerThicknessEdgeSum, &
                                                      tangentialVelocitySum
      real (kind=RKIND), dimension(:), allocatable, target :: RiSmoothed, BVFSmoothed, OBLDepths, interfaceForcings
      logical :: bulkRichardsonFlag

      real (kind=RKIND) :: langmuirEnhancementFactor, alphaAngle

      real(kind=RKIND), dimension(:,:), pointer :: stokesDriftZonalWavenumber
      real(kind=RKIND), dimension(:,:), pointer :: stokesDriftMeridionalWavenumber
      real(kind=RKIND), dimension(:), pointer :: significantWaveHeight
      real(kind=RKIND), dimension(:), pointer :: windStressZonal, windStressMeridional



      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing mixing-related fields
      ! note that the user can choose multiple options and the
      !   mixing fields have to be added/merged together
      !
      !-----------------------------------------------------------------

      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      !
      ! only build up viscosity/diffusivity if CVMix is turned on
      !
      if ( .not. cvmixOn ) return

      !
      ! set parameters
      !
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'cellMask', cellMask)

      !
      ! set pointers for fields related to position on sphere
      !
      call mpas_pool_get_array(meshPool, 'latCell', latCell)
      call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
      call mpas_pool_get_array(meshPool, 'fCell', fCell)

      !
      ! set pointers for fields related to vertical mesh
      !
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'minLevelCell', minLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'minLevelEdgeBot', minLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

      !
      ! set pointers for fields related to ocean forcing state
      !
      call mpas_pool_get_array(forcingPool, 'iceFraction', iceFraction)
      call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)
      call mpas_pool_get_array(forcingPool, 'windSpeed10m', windSpeed10m)
      call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
      call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional)
      call mpas_pool_get_array(forcingPool, 'langmuirNumber', langmuirNumber)

      call mpas_pool_get_array(forcingPool, 'significantWaveHeight', significantWaveHeight)
      call mpas_pool_get_array(forcingPool, 'stokesDriftZonalWavenumber', stokesDriftZonalWavenumber)
      call mpas_pool_get_array(forcingPool, 'stokesDriftMeridionalWavenumber', stokesDriftMeridionalWavenumber)

      nCells = nCellsArray( size(nCellsArray) )

      !$omp parallel
      !$omp do schedule(runtime)
      do iCell = 1, nCells
         vertViscTopOfCell(:, iCell) = 0.0_RKIND
         vertDiffTopOfCell(:, iCell) = 0.0_RKIND
         vertNonLocalFlux(:, :, iCell) = 0.0_RKIND
      end do
      !$omp end do
      !$omp end parallel

      nCells = nCellsArray( 3 )

      !
      ! start by adding the mininum background values to the viscocity/diffusivity arrays
      !
      if (cvmixBackgroundChoice == cvmixBackgroundTypeConstant) then
         !$omp parallel
         !$omp do schedule(runtime) private(k)
         do iCell = 1, nCells
            do k = 1, nVertLevelsP1
               vertViscTopOfCell(k, iCell) = vertViscTopOfCell(k, iCell) + config_cvmix_background_viscosity
               vertDiffTopOfCell(k, iCell) = vertDiffTopOfCell(k, iCell) + config_cvmix_background_diffusion
            end do
         end do
         !$omp end do
         !$omp end parallel
      endif

      !
      ! allocate selected cvmix variables and loop over columns
      !
      cvmix_variables % max_nlev = nVertLevels
      allocate(cvmix_variables % Mdiff_iface(nVertLevels+1))
      allocate(cvmix_variables % Tdiff_iface(nVertLevels+1))
      allocate(cvmix_variables % Sdiff_iface(nVertLevels+1))
      allocate(cvmix_variables % zw_iface(nVertLevels+1))
      allocate(cvmix_variables % dzw(nVertLevels+1))
      allocate(cvmix_variables % zt_cntr(nVertLevels))
      allocate(cvmix_variables % dzt(nVertLevels))
      allocate(cvmix_variables % kpp_Tnonlocal_iface(nVertLevels+1))
      allocate(cvmix_variables % kpp_Snonlocal_iface(nVertLevels+1))

      ! Initialize some of the cvmix variables that are not set later.
      cvmix_variables % Mdiff_iface(1:nVertLevels+1) = 0.0_RKIND
      cvmix_variables % Tdiff_iface(1:nVertLevels+1) = 0.0_RKIND
      cvmix_variables % Sdiff_iface(1:nVertLevels+1) = 0.0_RKIND

      allocate(OBLDepths(nVertLevels))
      allocate(interfaceForcings(nVertLevels))

      allocate(Nsqr_iface(nVertLevels+1))
      allocate(turbulentScalarVelocityScale(nVertLevels))
      allocate(RiSmoothed(nVertLevels+1))
      allocate(BVFSmoothed(nVertLevels+1))
      allocate(RiTemp(nVertLevels+1))

      allocate(normalVelocitySum(nVertLevels))
      allocate(tangentialVelocitySum(nVertLevels))
      allocate(potentialDensitySum(nVertLevels))
      allocate(surfaceAverageIndex(nVertLevels))
      allocate(deltaVelocitySquared(nVertLevels))
      allocate(layerThicknessSum(nVertLevels))
      allocate(layerThicknessEdgeSum(nVertLevels))

      do k = 1, nVertLevels
         Nsqr_iface(k) = 0.0_RKIND
         turbulentScalarVelocityScale(k) = 0.0_RKIND
      end do
      Nsqr_iface(nVertLevelsP1) = 0.0_RKIND

      if (config_use_active_wave) then
        !call ocn_stokes_drift_reconstruct(forcingPool)
        call ocn_surface_stokes_drift(forcingPool)
      end if

      call mpas_timer_start('cvmix cell loop', .false.)
      do kpp_stage = 1,2
      ! TODO, mdt: determine where the openmp data race is coming from in this loop
!      !$omp parallel
!      !$omp do schedule(runtime) &
!      !$omp private(k, RiSmoothed, RiTemp, nsmooth, BVFSmoothed, alphaAngle, &
!      !$omp         langmuirEnhancementFactor, bulkRichardsonNumberStop, bulkRichardsonFlag, &
!      !$omp         topIndex, kIndexOBL, deltaVelocitySquared, sigma, OBLDepths, interfaceForcings, &
!      !$omp         surfaceAverageIndex, sfc_layer_depth, kav, i, iEdge, tangentialVelocitySum, &
!      !$omp         normalVelocitySum, layerThicknessEdgeSum, normalVelocityAv, delU2, potentialDensitySum, &
!      !$omp         layerThicknessSum, blTemp, zonalWavenumberCoeff, meridionalWavenumberCoeff, &
!      !$omp         tangentialVelocityAv) &
!      !$omp firstprivate(cvmix_variables, Nsqr_iface, turbulentScalarVelocityScale)
      do iCell = 1, nCells

         ! specify geometry/location
         cvmix_variables % SeaSurfaceHeight = ssh(iCell)
         cvmix_variables % Coriolis = fCell(iCell)
         cvmix_variables % lat = latCell(iCell) * 180.0_RKIND / pi
         cvmix_variables % lon = lonCell(iCell) * 180.0_RKIND / pi

         ! fill vertical position of column
         ! CVMix assume top of ocean is at z=0, so building all z-coordinate data based on layerThickness
         do k = 1, minLevelCell(iCell) - 1
            cvmix_variables % zw_iface(k) = 0.0_RKIND
            cvmix_variables % zt_cntr(k) = 0.0_RKIND
            cvmix_variables % dzw(k) = 0.0_RKIND
            cvmix_variables % dzt(k) = 0.0_RKIND
         enddo
         cvmix_variables % zw_iface(minLevelCell(iCell)) = 0.0_RKIND
         cvmix_variables % dzw(minLevelCell(iCell)) = layerThickness(minLevelCell(iCell),iCell)/2.0_RKIND
         cvmix_variables % zt_cntr(minLevelCell(iCell)) = -layerThickness(minLevelCell(iCell),iCell)/2.0_RKIND
         do k=minLevelCell(iCell)+1,maxLevelCell(iCell)
            cvmix_variables % zw_iface(k) = cvmix_variables % zw_iface(k-1) - layerThickness(k-1,iCell)
            cvmix_variables % zt_cntr(k) = cvmix_variables %  zw_iface(k) - layerThickness(k,iCell)/2.0_RKIND
            cvmix_variables % dzw(k) = cvmix_variables % zt_cntr(k-1) - cvmix_variables % zt_cntr(k)
            cvmix_variables % dzt(k) = layerThickness(k,iCell)
         enddo
         k = maxLevelCell(iCell)+1
         cvmix_variables % zw_iface(k) = cvmix_variables % zw_iface(k-1) - layerThickness(k-1,iCell)
         cvmix_variables % dzw(k) = cvmix_variables % zt_cntr(k-1) - cvmix_variables % zw_iface(k)
         do k = maxLevelCell(iCell) + 1, nVertLevels
            cvmix_variables % zw_iface(k+1) = cvmix_variables % zw_iface(maxLevelCell(iCell)+1)
            cvmix_variables % zt_cntr(k) = cvmix_variables % zw_iface(maxLevelCell(iCell)+1)
            cvmix_variables % dzw(k+1) = 0.0_RKIND
            cvmix_variables % dzt(k) = 0.0_RKIND
         enddo

         ! fill the intent(in) convective adjustment
         cvmix_variables % nlev = maxLevelCell(iCell)
         cvmix_variables % OceanDepth = bottomDepth(iCell)
         cvmix_variables % WaterDensity_cntr => density(1:maxLevelCell(iCell), iCell)
         cvmix_variables % AdiabWaterDensity_cntr => displacedDensity(1:maxLevelCell(iCell), iCell)

         if (kpp_stage == 2) then
         if (cvmixBackgroundChoice == cvmixBackgroundTypeBryanLewis) then
            ! See equation 9 in Bryan, K., and Lewis, L. J. (1979), A water mass
            ! model of the World Ocean, J. Geophys. Res., 84( C5), 2503– 2517
            ! https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/JC084iC05p02503
            ! The standard atan formulation is approximated following equation (13) of Girones et al (2012)
            ! https://doi.org/10.1109/MSP.2012.2219677
            do k = minLevelCell(iCell),maxLevelCell(iCell)
               x = (abs(cvmix_variables % zw_iface(k)) - config_cvmix_BryanLewis_transitionDepth) / &
                   config_cvmix_BryanLewis_transitionWidth
               vertDiffTopOfCell(k, iCell) = vertDiffTopOfCell(k, iCell) + config_cvmix_BryanLewis_bl1 + &
                        (config_cvmix_BryanLewis_bl2/2.0_RKIND)*(sign(1.0_RKIND, x)*(0.596227_RKIND*     &
                        abs(x)+x**2)/(1.0_RKIND+1.192454_RKIND*abs(x)+x**2))
               vertViscTopOfCell(k, iCell) = config_cvmix_prandtl_number*vertDiffTopOfCell(k, iCell)
            enddo
         endif

         ! fill Ri
         RiSmoothed(1:minLevelCell(iCell)-1) = RiTopOfCell(minLevelCell(iCell),iCell)
         RiSmoothed(minLevelCell(iCell):maxLevelCell(iCell)) = RiTopOfCell(minLevelCell(iCell):maxLevelCell(iCell),iCell)
         RiSmoothed(maxLevelCell(iCell)+1) = RiSmoothed(maxLevelCell(iCell))
         RiTemp(1:maxLevelCell(iCell)+1) = RiSmoothed(1:maxLevelCell(iCell)+1)

         ! Use a 1-2-1 filter to remove 2 dz noise in RiTopOfCell
         ! Note: loop limits here must not be minLevelCell or non-bfb changes are introduced
         do nsmooth=1,config_cvmix_num_ri_smooth_loops
            do k=2,maxLevelCell(iCell)
                 RiSmoothed(k) = (RiTemp(k-1) + 2.0_RKIND*RiTemp(k) + RiTemp(k+1))  / 4.0_RKIND
            enddo
            RiTemp(1:maxLevelCell(iCell)) = RiSmoothed(1:maxLevelCell(iCell))
         enddo

         cvmix_variables%ShearRichardson_iface => RiSmoothed(minLevelCell(iCell):maxLevelCell(iCell)+1)

         endif

         ! fill BVF
         BVFSmoothed(minLevelCell(iCell):maxLevelCell(iCell)) = max(0.0_RKIND, &
                        BruntVaisalaFreqTop(minLevelCell(iCell):maxLevelCell(iCell),iCell))
         BVFSmoothed(1:minLevelCell(iCell)-1) = max(0.0_RKIND,BVFSmoothed(minLevelCell(iCell)))
         BVFSmoothed(maxLevelCell(iCell)+1) = max(0.0_RKIND,BVFSmoothed(maxLevelCell(iCell)))
         cvmix_variables%SqrBuoyancyFreq_iface => BVFSmoothed(minLevelCell(iCell):maxLevelCell(iCell)+1)

         ! fill the intent(in) KPP
         cvmix_variables % SurfaceFriction = surfaceFrictionVelocity(iCell)
         cvmix_variables % SurfaceBuoyancyForcing = surfaceBuoyancyForcing(iCell)
         cvmix_variables % BulkRichardson_cntr => bulkRichardsonNumber(minLevelCell(iCell):maxLevelCell(iCell), iCell)

         if (kpp_stage == 2) then
         if (config_use_cvmix_shear) then

            do k = minLevelCell(iCell), maxLevelCell(iCell) + 1
               cvmix_variables % Mdiff_iface(k) = 0.0_RKIND
               cvmix_variables % Tdiff_iface(k) = 0.0_RKIND
            end do
            call cvmix_coeffs_shear( &
                 cvmix_variables, &
                 cvmix_shear_params)
            ! add shear mixing to vertical viscosity/diffusivity
            ! at present, shear mixing adds in background values when using PP, but background is
            ! accounted for seperately. so remove bac    kground from shear mixing values

            if(cvmixShearMixingChoice==cvmixShearMixingTypePP) then
               do k = minLevelCell(iCell), maxLevelCell(iCell) + 1
                  cvmix_variables % Mdiff_iface(k) = cvmix_variables % Mdiff_iface(k) - config_cvmix_background_viscosity
                  cvmix_variables % Tdiff_iface(k) = cvmix_variables % Tdiff_iface(k) - config_cvmix_background_diffusion
               end do
            endif

            do k = minLevelCell(iCell), maxLevelCell(iCell)
               vertViscTopOfCell(k, iCell) = vertViscTopOfCell(k, iCell) + cvmix_variables % Mdiff_iface(k)
               vertDiffTopOfCell(k, iCell) = vertDiffTopOfCell(k, iCell) + cvmix_variables % Tdiff_iface(k)
            end do

         endif ! if (config_use_cvmix_shear)
          endif ! stage 2 shear compute
         ! call kpp ocean mixed layer scheme
         if (cvmixKPPOn) then

          if (kpp_stage==1) then
            if (config_use_cvmix_fixed_boundary_layer) then
               cvmix_variables % BoundaryLayerDepth = config_cvmix_kpp_boundary_layer_depth
               cvmix_variables % kOBL_depth = cvmix_kpp_compute_kOBL_depth(  &
                          zw_iface = cvmix_variables%zw_iface(minLevelCell(iCell):maxLevelCell(iCell)+1), &
                          zt_cntr = cvmix_variables%zt_cntr(minLevelCell(iCell):maxLevelCell(iCell)),    &
                          OBL_depth = cvmix_variables % BoundaryLayerDepth )

            else

              ! set stratification
              do k=minLevelCell(iCell),maxLevelCell(iCell)
                  Nsqr_iface(k) = BVFSmoothed(k)
              enddo
              Nsqr_iface(1:minLevelCell(iCell)-1) = Nsqr_iface(minLevelCell(iCell))
              k=min(maxLevelCell(iCell)+1,nVertLevels)
              Nsqr_iface(k:maxLevelCell(iCell)+1) = Nsqr_iface(k-1)

              ! compute Langmuir number and Langmuir enhancement factor
              if (associated(landIceMask)) then
                  landIceMaskValue = landIceMask(iCell)
              else
                  landIceMaskValue = 0
              endif

              if ( landIceMaskValue .eq. 0 .and. iceFraction(iCell) .lt. 0.05_RKIND) then
                 if (config_cvmix_kpp_use_theory_wave) then
                    langmuirNumber(iCell) =  sqrt(surfaceFrictionVelocity(iCell) / ( &
                         cvmix_kpp_ustokes_SL_model(windSpeed10m(iCell), &
                                                    boundaryLayerDepth(iCell), &
                                                    cvmix_global_params)+1e-15_RKIND) )
                    langmuirEnhancementFactor =  &
                         cvmix_kpp_EFactor_model(windSpeed10m(iCell), &
                                                 surfaceFrictionVelocity(iCell), &
                                                 boundaryLayerDepth(iCell), &
                                                 cvmix_global_params)
                 else if (config_cvmix_kpp_use_active_wave) then
                    call ocn_stokes_drift_langmuir_number(windStressZonal(iCell), &
                                                       windStressMeridional(iCell), &
                                                       surfaceFrictionVelocity(iCell), &
                                                       significantWaveHeight(iCell), &
                                                       boundaryLayerDepth(iCell), &
                                                       stokesDriftZonalWavenumber(:,iCell), &
                                                       stokesDriftMeridionalWavenumber(:,iCell), &
                                                       alphaAngle, &
                                                       langmuirNumber(iCell))
                    call ocn_stokes_drift_kpp_enhancement_factor(alphaAngle, &
                                                              langmuirNumber(iCell), &
                                                              langmuirEnhancementFactor)
                 endif
              else
                 ! arbitrarily large Langmuir number
                 langmuirNumber(iCell) = 1.e10_RKIND
                 langmuirEnhancementFactor = 1.0_RKIND
              end if

              ! compute bulk Richardson number
              ! assume boundary layer depth is at bottom of every kIndexOBL cell
              bulkRichardsonNumberStop = config_cvmix_kpp_stop_OBL_search * config_cvmix_kpp_criticalBulkRichardsonNumber
              bulkRichardsonFlag = .false.
              topIndex = minLevelCell(iCell)
!             call mpas_timer_start('Bulk Richardson kIndexOBL loops')
              ! compute the index of the last cell in the KPP defined surface layer
              ! this index is used for the necessary surface layer averages of buoyancy and momentum
              do kIndexOBL = minLevelCell(iCell), maxLevelCell(iCell)

                 ! Reset deltaVelocitySquared and bulkRichardsonNumber at this layer for the later computation
                 deltaVelocitySquared(kIndexOBL) = 0.0_RKIND
                 bulkRichardsonNumber(kIndexOBL, iCell) = bulkRichardsonNumberStop - 1.0_RKIND

                 ! set OBL at bottom of kIndexOBL cell for computation of bulk Richardson number
                 cvmix_variables % BoundaryLayerDepth = abs(cvmix_variables % zw_iface(kIndexOBL+1))
                 sigma = -cvmix_variables % zt_cntr(kIndexOBL) / cvmix_variables % BoundaryLayerDepth

                 OBLDepths(kIndexOBL) = abs(cvmix_variables % zw_iface(kIndexOBL+1))
                 interfaceForcings(kIndexOBL) = cvmix_variables % SurfaceBuoyancyForcing

                 ! initialize the surfaceAverageIndex for cases when the if statement below is not true
                 surfaceAverageIndex(kIndexOBL) = minLevelCell(iCell)
                 ! move progressively downward to find the bottom most layer within the surface layer
                 sfc_layer_depth = cvmix_variables % BoundaryLayerDepth * config_cvmix_kpp_surface_layer_extent
                 do kav=topIndex,kIndexOBL
                    if(cvmix_variables%zw_iface(kav+1) < -sfc_layer_depth) then
                       surfaceAverageIndex(kIndexOBL) = kav
                       exit
                    end if
                 enddo
                 topIndex = max(1, surfaceAverageIndex(kIndexOBL) - 1)
              end do

              ! compute the turbulent scales in order to compute the bulk Richardson number
              call cvmix_kpp_compute_turbulent_scales( &
                   sigma_coord = config_cvmix_kpp_surface_layer_extent, &
                   OBL_depth = OBLDepths(minLevelCell(iCell):maxLevelCell(iCell)), &
                   surf_buoy_force = interfaceForcings(minLevelCell(iCell):maxLevelCell(iCell)), &
                   surf_fric_vel = cvmix_variables % SurfaceFriction, &
                   w_s = turbulentScalarVelocityScale(minLevelCell(iCell):maxLevelCell(iCell)))

              ! averaging over a surface layer assuming that BLdepth is cell bottom
              ! Build deltaVelocitySquared
              do i = 1, nEdgesOnCell(iCell)
                 iEdge = edgesOnCell(i, iCell)

                 deltaVelocitySquared(1:minLevelEdgeBot(iEdge)) = 0.0_RKIND

                 normalVelocitySum(minLevelEdgeBot(iEdge)) = normalVelocity(minLevelEdgeBot(iEdge), iEdge)* &
                                                          layerThickEdgeMean(minLevelEdgeBot(iEdge),iEdge)
                 tangentialVelocitySum(minLevelEdgeBot(iEdge)) = tangentialVelocity(minLevelEdgeBot(iEdge),iEdge)* &
                                                          layerThickEdgeMean(minLevelEdgeBot(iEdge),iEdge)
                 layerThicknessEdgeSum(minLevelEdgeBot(iEdge)) = layerThickEdgeMean(minLevelEdgeBot(iEdge),iEdge)

                 do kIndexOBL = minLevelEdgeBot(iEdge)+1, maxLevelEdgeTop(iEdge)
                    normalVelocitySum(kIndexOBL) = normalVelocitySum(kIndexOBL-1) + &
                                      layerThickEdgeMean(kIndexOBL, iEdge)*normalVelocity(kIndexOBL, iEdge)
                    tangentialVelocitySum(kIndexOBL) = tangentialVelocitySum(kIndexOBL-1) + &
                                      layerThickEdgeMean(kIndexOBL, iEdge)*tangentialVelocity(kIndexOBL, iEdge)
                    layerThicknessEdgeSum(kIndexOBL) = layerThicknessEdgeSum(kIndexOBL-1) + layerThickEdgeMean(kIndexOBL, iEdge)
                 end do

                 do kIndexOBL = minLevelEdgeBot(iEdge), maxLevelEdgeTop(iEdge)
                    normalVelocityAv = normalVelocitySum(surfaceAverageIndex(kIndexOBL )) / &
                        layerThicknessEdgeSum(surfaceAverageIndex(kIndexOBL))
                    tangentialVelocityAv = tangentialVelocitySum(surfaceAverageIndex(kIndexOBL)) / &
                        layerThicknessEdgeSum(surfaceAverageIndex(kIndexOBL))
                    delU2 = (normalVelocityAv - normalVelocity(kIndexOBL, iEdge))**2 + &
                            (tangentialVelocityAv - tangentialVelocity(kIndexOBL, iEdge))**2
                    deltaVelocitySquared(kIndexOBL) = deltaVelocitySquared(kIndexOBL) + edgeAreaFractionOfCell(i,iCell) * delU2
                 end do
              end do

              potentialDensitySum(minLevelCell(iCell)) = potentialDensity(minLevelCell(iCell), iCell)* &
                                                         layerThickness(minLevelCell(iCell), iCell)
              layerThicknessSum(minLevelCell(iCell)) = layerThickness(minLevelCell(iCell), iCell)
              do kIndexOBL = minLevelCell(iCell)+1, maxLevelCell(iCell)
                layerThicknessSum(kIndexOBL) = layerThicknessSum(kIndexOBL-1) + layerThickness(kIndexOBL, iCell)
                potentialDensitySum(kIndexOBL) = potentialDensitySum(kIndexOBL-1) + &
                    layerThickness(kIndexOBL, iCell)*potentialDensity(kIndexOBL, iCell)
              end do

              do kIndexOBL = minLevelCell(iCell), maxLevelCell(iCell)
                ! !compute shear contribution assuming BLdepth is cell bottom
                !Note that the factor of two is from averaging dot products to cell centers on a C-grid
                bulkRichardsonNumberShear(kIndexOBL,iCell) = max(deltaVelocitySquared(kIndexOBL), 1.0e-15_RKIND)

                bulkRichardsonNumberBuoy(kIndexOBL,iCell) = gravity * (potentialDensity(kIndexOBL, iCell) &
                                  - potentialDensitySum(surfaceAverageIndex(kIndexOBL)) &
                                  / layerThicknessSum(surfaceAverageIndex(kIndexOBL))) / rho_sw
              end do ! do kIndexOBL
!             call mpas_timer_stop('Bulk Richardson kIndexOBL loops')

              cvmix_variables % bulkRichardson_cntr(:) = cvmix_kpp_compute_bulk_Richardson( &
                   zt_cntr = cvmix_variables % zt_cntr(minLevelCell(iCell):maxLevelCell(iCell)), &
                   delta_buoy_cntr = bulkRichardsonNumberBuoy(minLevelCell(iCell):maxLevelCell(iCell),iCell), &
                   delta_Vsqr_cntr = bulkRichardsonNumberShear(minLevelCell(iCell):maxLevelCell(iCell),iCell), &
                   ws_cntr = turbulentScalarVelocityScale(:), &
                   Nsqr_iface = Nsqr_iface(minLevelCell(iCell):maxLevelCell(iCell)+1), &
                   EFactor = langmuirEnhancementFactor, &
                   LaSL = langmuirNumber(iCell), &
                   bfsfc = cvmix_variables%SurfaceBuoyancyForcing, &
                   ustar = cvmix_variables%SurfaceFriction )
              ! TODO: the surface buoyancy forcing here does not include penetrative solar radiation? <05-09-18, Qing Li> !

              ! each level of bulk Richardson is computed as if OBL resided at bottom of that level

              call cvmix_kpp_compute_OBL_depth(  &
                   Ri_bulk = bulkRichardsonNumber(minLevelCell(iCell):maxLevelCell(iCell),iCell), &
                   zw_iface = cvmix_variables % zw_iface(minLevelCell(iCell):maxLevelCell(iCell)+1), &
                   OBL_depth = cvmix_variables % BoundaryLayerDepth, &
                   kOBL_depth = cvmix_variables % kOBL_depth, &
                   zt_cntr = cvmix_variables % zt_cntr(minLevelCell(iCell):maxLevelCell(iCell)), &
                   surf_fric = cvmix_variables % SurfaceFriction, &
                   surf_buoy = cvmix_variables % SurfaceBuoyancyForcing, &
                   Coriolis = cvmix_variables % Coriolis)

             endif  ! if (config_use_cvmix_fixed_boundary_layer) then

             ! apply minimum limit to OBL
             if(cvmix_variables % BoundaryLayerDepth .lt. layerThickness(minLevelCell(iCell),iCell)/2.0_RKIND) then
                cvmix_variables % BoundaryLayerDepth = layerThickness(minLevelCell(iCell),iCell)/2.0_RKIND
                cvmix_variables % kOBL_depth = cvmix_kpp_compute_kOBL_depth(  &
                          zw_iface = cvmix_variables%zw_iface(minLevelCell(iCell):maxLevelCell(iCell)+1),&
                          zt_cntr = cvmix_variables%zt_cntr(minLevelCell(iCell):maxLevelCell(iCell)),    &
                          OBL_depth = cvmix_variables % BoundaryLayerDepth )
             endif

             ! apply minimum limit to OBL under sea-ice
             if(iceFraction(iCell).gt.0.15_RKIND) then
             if(cvmix_variables % BoundaryLayerDepth .lt. configure_cvmix_kpp_minimum_OBL_under_sea_ice) then
                cvmix_variables % BoundaryLayerDepth = configure_cvmix_kpp_minimum_OBL_under_sea_ice
                cvmix_variables % kOBL_depth = cvmix_kpp_compute_kOBL_depth(  &
                          zw_iface = cvmix_variables%zw_iface(minLevelCell(iCell):maxLevelCell(iCell)+1),&
                          zt_cntr = cvmix_variables%zt_cntr(minLevelCell(iCell):maxLevelCell(iCell)),    &
                          OBL_depth = cvmix_variables % BoundaryLayerDepth )
             endif
             endif

             ! apply maximum limit to OBL
             if(cvmix_variables % BoundaryLayerDepth .gt. abs(cvmix_variables%zt_cntr(maxLevelCell(iCell)))) then
                cvmix_variables % BoundaryLayerDepth = abs(cvmix_variables%zt_cntr(maxLevelCell(iCell)))
                cvmix_variables % kOBL_depth = cvmix_kpp_compute_kOBL_depth(  &
                          zw_iface = cvmix_variables%zw_iface(minLevelCell(iCell):maxLevelCell(iCell)+1), &
                          zt_cntr = cvmix_variables%zt_cntr(minLevelCell(iCell):maxLevelCell(iCell)),    &
                          OBL_depth = cvmix_variables % BoundaryLayerDepth )

             endif

             boundaryLayerDepth(iCell) = cvmix_variables % BoundaryLayerDepth
     endif !kpp stage 1 -- boundary layer compute

     if (kpp_stage == 2) then
            ! copy data into cvmix_variables
            do k = minLevelCell(iCell), maxLevelCell(iCell) + 1
              cvmix_variables % Mdiff_iface(k) = vertViscTopOfCell(k, iCell)
              cvmix_variables % Tdiff_iface(k) = vertDiffTopOfCell(k, iCell)
            end do

            !must reapply max and min limits to boundaryLayerDepth
            blTemp = max(boundaryLayerDepth(iCell), layerThickness(minLevelCell(iCell),iCell)/2.0_RKIND)
            boundaryLayerDepth(iCell) = min(blTemp, abs(cvmix_variables%zt_cntr(maxLevelCell(iCell))))

            cvmix_variables % kOBL_depth = cvmix_kpp_compute_kOBL_depth(  &
                          zw_iface = cvmix_variables%zw_iface(minLevelCell(iCell):maxLevelCell(iCell)+1), &
                          zt_cntr = cvmix_variables%zt_cntr(minLevelCell(iCell):maxLevelCell(iCell)),    &
                          OBL_depth = boundaryLayerDepth(iCell) )

            ! update Langmuir enhancement factor
            if (associated(landIceMask)) then
                  landIceMaskValue = landIceMask(iCell)
            else
                  landIceMaskValue = 0
            endif

            if (landIceMaskValue .eq. 0 .and. iceFraction(iCell) .lt. 0.05_RKIND) then
               if (config_cvmix_kpp_use_theory_wave) then
                  langmuirEnhancementFactor =  &
                       cvmix_kpp_EFactor_model(windSpeed10m(iCell), &
                                               surfaceFrictionVelocity(iCell), &
                                               boundaryLayerDepth(iCell), &
                                               cvmix_global_params)
               else if (config_cvmix_kpp_use_active_wave) then
                  call ocn_stokes_drift_langmuir_number(windStressZonal(iCell), &
                                                       windStressMeridional(iCell), &
                                                       surfaceFrictionVelocity(iCell), &
                                                       significantWaveHeight(iCell), &
                                                       boundaryLayerDepth(iCell), &
                                                       stokesDriftZonalWavenumber(:,iCell), &
                                                       stokesDriftMeridionalWavenumber(:,iCell), &
                                                       alphaAngle, &
                                                       langmuirNumber(iCell))
                  call ocn_stokes_drift_kpp_enhancement_factor(alphaAngle, &
                                                       langmuirNumber(iCell), &
                                                       langmuirEnhancementFactor)
               endif
            else
               langmuirEnhancementFactor = 1.0_RKIND
            endif
!           call mpas_timer_start('cvmix coeffs kpp', .false.)
            call cvmix_coeffs_kpp(                                              &
                          Mdiff_out = cvmix_variables % Mdiff_iface(minLevelCell(iCell):maxLevelCell(iCell)+1),       &
                          Tdiff_out = cvmix_variables % Tdiff_iface(minLevelCell(iCell):maxLevelCell(iCell)+1),       &
                          Sdiff_out = cvmix_variables % Sdiff_iface(minLevelCell(iCell):maxLevelCell(iCell)+1),       &
                          zw = cvmix_variables%zw_iface(minLevelCell(iCell):maxLevelCell(iCell)+1),            &
                          zt = cvmix_variables%zt_cntr(minLevelCell(iCell):maxLevelCell(iCell)),               &
                          old_Mdiff = cvmix_variables%Mdiff_iface(minLevelCell(iCell):maxLevelCell(iCell)+1),         &
                          old_Tdiff = cvmix_variables%Tdiff_iface(minLevelCell(iCell):maxLevelCell(iCell)+1),         &
                          old_Sdiff = cvmix_variables%Sdiff_iface(minLevelCell(iCell):maxLevelCell(iCell)+1),         &
                          OBL_depth = boundaryLayerDepth(iCell),                   &
                          kOBL_depth = cvmix_variables%kOBL_depth,                           &
                          Tnonlocal = cvmix_variables%kpp_Tnonlocal_iface(minLevelCell(iCell):maxLevelCell(iCell)+1), &
                          Snonlocal = cvmix_variables%kpp_Snonlocal_iface(minLevelCell(iCell):maxLevelCell(iCell)+1), &
                          surf_fric = cvmix_variables%SurfaceFriction,                      &
                          surf_buoy = cvmix_variables%SurfaceBuoyancyForcing,               &
                          nlev = maxLevelCell(iCell) - minLevelCell(iCell) + 1,             &
                          max_nlev = nVertLevels,         &
                          Langmuir_EFactor = langmuirEnhancementFactor )
!           call mpas_timer_stop('cvmix coeffs kpp')

            ! intent out of BoundaryLayerDepth is boundary layer depth measured in meters and vertical index
            indexBoundaryLayerDepth(iCell) = cvmix_variables % kOBL_depth

            do k = minLevelCell(iCell), maxLevelCell(iCell) + 1
                vertViscTopOfCell(k, iCell) = cvmix_variables % Mdiff_iface(k)
                vertDiffTopOfCell(k, iCell) = cvmix_variables % Tdiff_iface(k)
            end do

            ! store non-local flux terms
            ! these flux terms must be multiplied by the surfaceTracerFlux field
            ! the tracer tendency is then the vertical divergence of vertNonLocalFlux*surfaceTracerFlux
            ! both of these operations are done in ocn_tracer_nonlocalflux_tend routine
            vertNonLocalFlux(index_vertNonLocalFluxTemp,:,iCell) = cvmix_variables % kpp_Tnonlocal_iface(:)

    endif ! kpp stage 2
         endif !if (config_use_cvmix_kpp)

         if ( kpp_stage == 2) then
         ! call convective mixing scheme
         if (config_use_cvmix_convection) then
            do k = 1, maxLevelCell(iCell) + 1
               cvmix_variables % Mdiff_iface(k) = 0.0_RKIND
               cvmix_variables % Tdiff_iface(k) = 0.0_RKIND
            end do
            call cvmix_coeffs_conv( CVmix_vars = cvmix_variables )

            ! add convective mixing to vertical viscosity/diffusivity
            ! if using KPP, then do not apply convective mixing within the ocean boundary layer
            if(config_use_cvmix_kpp) then
               do k = ceiling(indexBoundaryLayerDepth(iCell)) + 1, maxLevelCell(iCell)
                  vertViscTopOfCell(k,iCell) = vertViscTopOfCell(k,iCell) + cvmix_variables % Mdiff_iface(k)
                  vertDiffTopOfCell(k,iCell) = vertDiffTopOfCell(k,iCell) + cvmix_variables % Tdiff_iface(k)
               enddo
            else
               vertViscTopOfCell(1:minLevelCell(iCell)-1, iCell) = 0.0_RKIND
               vertDiffTopOfCell(1:minLevelCell(iCell)-1, iCell) = 0.0_RKIND
               do k = minLevelCell(iCell), maxLevelCell(iCell) + 1
                  vertViscTopOfCell(k, iCell) = vertViscTopOfCell(k, iCell) + cvmix_variables % Mdiff_iface(k)
                  vertDiffTopOfCell(k, iCell) = vertDiffTopOfCell(k, iCell) + cvmix_variables % Tdiff_iface(k)
               end do
            endif
         endif  ! if (config_use_cvmix_convection)

         !
         ! put tidal mixing here
         !

         !
         ! put double diffusion mxing here
         !


         ! computation of viscosity/diffusivity complete
         ! impose no-flux boundary conditions at top and bottom by zero viscosity/diffusivity
         vertViscTopOfCell(1:minLevelCell(iCell), iCell) = 0.0_RKIND
         vertDiffTopOfCell(1:minLevelCell(iCell), iCell) = 0.0_RKIND
         do k = maxLevelCell(iCell)+1, nVertLevelsP1
            vertViscTopOfCell(k, iCell)=0.0_RKIND
            vertDiffTopOfCell(k, iCell)=0.0_RKIND
         end do

 endif ! kpp stage 2 convection calc
      end do  ! do iCell=1,mesh%nCells
!      !$omp end do
!      !$omp end parallel

      if (kpp_stage == 1 .and. config_cvmix_use_BLD_smoothing) then ! smooth boundary layer
         nCells = nCellsArray(2)

         !$omp parallel
         !$omp do schedule(runtime) private(nEdges, areaSum, edgeCount, iEdge, iNeighbor)
         do iCell=1,nCells
            nEdges = nEdgesOnCell(iCell)
            boundaryLayerDepthSmooth(iCell) = 0.0_RKIND
            areaSum = 0.0_RKIND
            edgeCount = 0

            do iEdge = 1,nEdges
               iNeighbor = cellsOnCell(iEdge,iCell)
               boundaryLayerDepthSmooth(iCell) = boundaryLayerDepthSmooth(iCell) +  &
                             2.0_RKIND * cellMask(minLevelCell(iNeighbor),iNeighbor) * areaCell(iNeighbor)    &
                             * boundaryLayerDepth(iNeighbor)
               areaSum = areaSum + 2.0_RKIND * areaCell(iNeighbor) * cellMask(minLevelCell(iNeighbor),iNeighbor)
               edgeCount = edgeCount + cellMask(minLevelCell(iNeighbor),iNeighbor)
            end do
            areaSum = areaSum + edgeCount * areaCell(iCell)
            boundaryLayerDepthSmooth(iCell) = boundaryLayerDepthSmooth(iCell) +   &
                         boundaryLayerDepth(iCell) * edgeCount * areaCell(iCell)
            boundaryLayerDepthSmooth(iCell) = boundaryLayerDepthSmooth(iCell) / areaSum
         end do
         !$omp end do

         !$omp do schedule(runtime)
          do iCell=1, nCells
             boundaryLayerDepth(iCell) = boundaryLayerDepthSmooth(iCell)
          enddo
         !$omp end do
         !$omp end parallel

       endif !stage 1 BLD smoothing
      end do ! kpp_stage
      call mpas_timer_stop('cvmix cell loop')

      ! dellocate cmvix variables
      deallocate(cvmix_variables % Mdiff_iface)
      deallocate(cvmix_variables % Tdiff_iface)
      deallocate(cvmix_variables % Sdiff_iface)
      deallocate(cvmix_variables % zw_iface)
      deallocate(cvmix_variables % dzw)
      deallocate(cvmix_variables % zt_cntr)
      deallocate(cvmix_variables % dzt)
      deallocate(cvmix_variables % kpp_Tnonlocal_iface)
      deallocate(cvmix_variables % kpp_Snonlocal_iface)

      deallocate(Nsqr_iface)
      deallocate(turbulentScalarVelocityScale)
      deallocate(RiSmoothed)
      deallocate(RiTemp)
      deallocate(BVFSmoothed)
      deallocate(normalVelocitySum)
      deallocate(tangentialVelocitySum)
      deallocate(potentialDensitySum)
      deallocate(surfaceAverageIndex)
      deallocate(deltaVelocitySquared)
      deallocate(layerThicknessEdgeSum)
      deallocate(layerThicknessSum)

      deallocate(OBLDepths)
      deallocate(interfaceForcings)

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_coefs_cvmix_build!}}}

!***********************************************************************
!
!  routine ocn_vmix_cvmix_init
!
!> \brief   Initializes ocean vertical mixing quantities by using
!> \ get and puts into CVMix
!> \author  Todd Ringler
!> \date    04 February 2013
!> \details
!>  This routine initializes a variety of quantities related to
!>  vertical mixing in the ocean. Parameters are set by calling into CVMix
!
!-----------------------------------------------------------------------


   subroutine ocn_vmix_cvmix_init(domain,err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      implicit none

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      integer, pointer :: nVertLevels

      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

      cvmixOn = config_use_cvmix
      backgroundVisc = config_cvmix_background_viscosity
      backgroundDiff = config_cvmix_background_diffusion
      cvmixConvectionOn = config_use_cvmix_convection
      cvmixKPPOn = config_use_cvmix_kpp

      !
      ! only initialize if CVMix is turned on
      !
      if (.not.config_use_cvmix) return

      !
      ! When CVMix is turned on, all other vertical mixing schemes should be off
      ! Test to make sure this is the case.
      !
      ! test here, err=1 if a problem

      !
      ! pull nVertLevels out of the mesh structure
      !
      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nVertLevels', nVertLevels)

      !
      ! put global parameters into CVMix
      !
      call cvmix_put(cvmix_global_params,  'max_nlev', nVertLevels)
      call cvmix_put(cvmix_global_params,  'prandtl',  config_cvmix_prandtl_number)
      call cvmix_put(cvmix_global_params,  'Gravity',  gravity)

      !
      ! set user choice for background scheme
      !
      select case (trim(config_cvmix_background_scheme))

      case ('none', 'None', 'NONE')

         cvmixBackgroundChoice = cvmixBackgroundTypeNone

      case ('constant', 'Constant', 'CONSTANT')

         cvmixBackgroundChoice = cvmixBackgroundTypeConstant

      case ('bryanlewis', 'bryanLewis', 'BryanLewis', 'BRYANLEWIS')

         cvmixBackgroundChoice = cvmixBackgroundTypeBryanLewis

      case default

         call mpas_log_write(&
            'Invalid choice for config_cvmix_background_scheme.  Choices are: none, &
             &constant, BryanLewis')

         err = 1

      end select

      !
      ! set user choice for background scheme
      !
      select case (trim(config_cvmix_shear_mixing_scheme))

      case ('pp', 'Pp', 'PP')

         cvmixShearMixingChoice = cvmixShearMixingTypePP

      case ('kpp', 'Kpp', 'KPP')

         cvmixShearMixingChoice = cvmixShearMixingTypeKPP

      case default

         call mpas_log_write(&
            'Invalid choice for config_cvmix_shear_mixing_scheme  Choices are: PP, KPP')

         err = 1

      end select

      !
      ! initialize background mixing
      !
      if ( cvmixBackgroundChoice == cvmixBackgroundTypeConstant ) then
           call cvmix_init_bkgnd( &
                  bkgnd_Tdiff = config_cvmix_background_diffusion, &
                  bkgnd_Mdiff = config_cvmix_background_viscosity, &
                  CVmix_bkgnd_params_user = cvmix_background_params)
      endif

      !
      ! initialize shear-based mixing
      !
      if (config_use_cvmix_shear) then
        if (cvmixBackgroundChoice == cvmixBackgroundTypeNone .and. cvmixShearMixingChoice == cvmixShearMixingTypePP) then
            call mpas_log_write("config_use_cvmix_shear cannot be used with with config_cvmix_shear_mixing_scheme = 'PP'" &
                 // "       without config_cvmix_background_scheme = 'constant' or 'BryanLewis'", MPAS_LOG_CRIT)
            err = 1
            return
        end if
        call cvmix_init_shear( &
               cvmix_shear_params, &
               mix_scheme = config_cvmix_shear_mixing_scheme, &
               PP_nu_zero = config_cvmix_shear_PP_nu_zero, &
               PP_alpha = config_cvmix_shear_PP_alpha, &
               PP_exp = config_cvmix_shear_PP_exp, &
               KPP_nu_zero = config_cvmix_shear_KPP_nu_zero, &
               KPP_Ri_zero = config_cvmix_shear_KPP_Ri_zero, &
               KPP_exp = config_cvmix_shear_KPP_exp)
      endif

      !
      ! initialize convective mixing
      !
      if (config_use_cvmix_convection) then

        ! config_cvmix_convective_basedOnBVF is not supported at this time
        if (.not.config_cvmix_convective_basedOnBVF) then
            call mpas_log_write("config_cvmix_convective_basedOnBVF = .false. is not supported. Change to true.", MPAS_LOG_CRIT)
            err = 1
            return
        endif

        call cvmix_init_conv( &
               convect_diff = config_cvmix_convective_diffusion,  &
               convect_visc = config_cvmix_convective_viscosity,  &
               lBruntVaisala = config_cvmix_convective_basedOnBVF, &
               BVsqr_convect = config_cvmix_convective_triggerBVF )
      endif

      !
      ! initialize tidal mixing
      !  (at present, tidal mixing can only use CVMix default parameter settings)
      !
      if (config_use_cvmix_tidal_mixing) then
        call cvmix_init_tidal(cvmix_tidal_params,'Simmons')
      endif

      !
      ! initialize double diffusion
      !  (at present, double diffusion can only use CVMix default parameter settings)
      !
      if (config_use_cvmix_double_diffusion) then
        call cvmix_init_ddiff( )
      endif

      !
      ! initialize KPP boundary layer scheme
      !
      if (config_use_cvmix_kpp) then
        if(config_cvmix_kpp_matching.eq."MatchBoth") then
           call mpas_log_write( &
              "Use of option MatchBoth is discouraged, use SimpleShapes instead", &
               MPAS_LOG_WARN)
        elseif(.not. config_cvmix_kpp_matching.eq."SimpleShapes") then
           call mpas_log_write( &
              "Unknown value for config_cvmix_kpp_matching., supported values are:" // &
              "         SimpleShapes or MatchBoth", &
              MPAS_LOG_CRIT)
           err = 1
           return
        endif

        if (trim(config_cvmix_kpp_langmuir_mixing_opt) .ne. "NONE" .and. &
            trim(config_cvmix_kpp_langmuir_mixing_opt) .ne. "LWF16" .and. &
            trim(config_cvmix_kpp_langmuir_mixing_opt) .ne. "RWHGK16") then
           call mpas_log_write( &
              "Unknown value for config_cvmix_kpp_langmuir_mixing_opt, supported values are:" // &
              "         NONE or LWF16 or RWHGK16", &
              MPAS_LOG_CRIT)
           err = 1
           return
        endif

        if (trim(config_cvmix_kpp_langmuir_entrainment_opt) .ne. "NONE" .and. &
            trim(config_cvmix_kpp_langmuir_entrainment_opt) .ne. "LWF16" .and. &
            trim(config_cvmix_kpp_langmuir_entrainment_opt) .ne. "LF17" .and. &
            trim(config_cvmix_kpp_langmuir_entrainment_opt) .ne. "RWHGK16") then
           call mpas_log_write( &
              "Unknown value for config_cvmix_kpp_langmuir_entrainment_opt, supported values are:" // &
              "         NONE or LWF16 or LF17 or RWHGK16", &
              MPAS_LOG_CRIT)
           err = 1
           return
        endif

        call cvmix_init_kpp ( &
               ri_crit = config_cvmix_kpp_criticalBulkRichardsonNumber, &
               interp_type = config_cvmix_kpp_interpolationOMLType, &
               interp_type2 = config_cvmix_kpp_interpolationOMLType, &
               lEkman = config_cvmix_kpp_EkmanOBL, &
               lMonOb = config_cvmix_kpp_MonObOBL, &
               MatchTechnique = config_cvmix_kpp_matching, &
               surf_layer_ext = config_cvmix_kpp_surface_layer_extent, &
               langmuir_mixing_str = config_cvmix_kpp_langmuir_mixing_opt, &
               langmuir_entrainment_str = config_cvmix_kpp_langmuir_entrainment_opt, &
               lenhanced_diff = config_cvmix_kpp_use_enhanced_diff)
      endif


   !--------------------------------------------------------------------

   end subroutine ocn_vmix_cvmix_init!}}}

!***********************************************************************

end module ocn_vmix_cvmix

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
