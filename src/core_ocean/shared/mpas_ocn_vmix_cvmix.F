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
   use mpas_io_units
   use mpas_constants
   use ocn_constants

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

   logical :: cvmixOn, cvmixBackgroundOn, cvmixConvectionOn, cvmixKPPOn
   real (kind=RKIND) :: backgroundVisc, backgroundDiff


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

   subroutine ocn_vmix_coefs_cvmix_build(meshPool, statePool, forcingPool, diagnosticsPool, err, timeLevelIn)!{{{

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
         diagnosticsPool   !< Input/Output: diagnostic information

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
        maxLevelCell, nEdgesOnCell

      real (kind=RKIND), dimension(:), pointer :: &
        latCell, lonCell, bottomDepth, surfaceBuoyancyForcing, surfaceFrictionVelocity, fCell, &
        boundaryLayerDepth, ssh, indexBoundaryLayerDepth, dcEdge, dvEdge, areaCell, iceFraction

      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfCell, vertDiffTopOfCell, layerThickness, &
        zMid, zTop, density, displacedDensity, potentialDensity, &
        bulkRichardsonNumber, RiTopOfCell, BruntVaisalaFreqTop, &
        bulkRichardsonNumberBuoy, bulkRichardsonNumberShear, unresolvedShear, normalVelocity

      real (kind=RKIND), dimension(:,:,:), pointer :: vertNonLocalFlux
      integer, pointer :: index_vertNonLocalFluxTemp
      integer, dimension(:,:), pointer :: edgesOnCell

      logical, pointer :: config_use_cvmix_shear, config_use_cvmix_convection, config_use_cvmix_kpp
      logical, pointer :: config_use_cvmix_fixed_boundary_layer
      real (kind=RKIND), pointer :: config_cvmix_kpp_stop_OBL_search, config_cvmix_kpp_criticalBulkRichardsonNumber
      real (kind=RKIND), pointer :: config_cvmix_kpp_boundary_layer_depth, config_cvmix_kpp_surface_layer_extent
      real (kind=RKIND), pointer :: configure_cvmix_kpp_minimum_OBL_under_sea_ice
      character (len=StrKIND), pointer :: config_cvmix_shear_mixing_scheme, config_cvmix_kpp_matching

      integer :: k, i, iCell, jCell, iNeighbor, iter, timeLevel, kIndexOBL, kav, iEdge, nCells
      integer :: topIndex
      integer, pointer :: nVertLevels, nVertLevelsP1
      integer, dimension(:), pointer :: nCellsArray
      integer, dimension(:), allocatable :: surfaceAverageIndex

      real (kind=RKIND) :: r, layerSum, bulkRichardsonNumberStop, sfc_layer_depth, invAreaCell
      real (kind=RKIND) :: normalVelocityAv, factor, delU2
      real (kind=RKIND) :: sigma, turbulentScalarVelocityScalePoint
      real (kind=RKIND), dimension(:), allocatable :: Nsqr_iface, turbulentScalarVelocityScale, &
                                                      deltaVelocitySquared, normalVelocitySum, &
                                                      potentialDensitySum
      real (kind=RKIND), dimension(:), allocatable, target :: RiSmoothed, BVFSmoothed, OBLDepths, interfaceForcings
      logical :: bulkRichardsonFlag

      real (kind=RKIND), pointer :: config_cvmix_background_viscosity, config_cvmix_background_diffusion

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
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_background_viscosity', config_cvmix_background_viscosity)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_background_diffusion', config_cvmix_background_diffusion)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_stop_OBL_search', config_cvmix_kpp_stop_OBL_search)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_criticalBulkRichardsonNumber', &
                                config_cvmix_kpp_criticalBulkRichardsonNumber)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_shear', config_use_cvmix_shear)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_convection', config_use_cvmix_convection)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_shear_mixing_scheme', config_cvmix_shear_mixing_scheme)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_matching', config_cvmix_kpp_matching)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_fixed_boundary_layer', config_use_cvmix_fixed_boundary_layer)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_boundary_layer_depth', config_cvmix_kpp_boundary_layer_depth)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_surface_layer_extent', config_cvmix_kpp_surface_layer_extent)
      call mpas_pool_get_config(ocnConfigs, 'configure_cvmix_kpp_minimum_OBL_under_sea_ice', &
                                configure_cvmix_kpp_minimum_OBL_under_sea_ice)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, timeLevel)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)

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
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)
      call mpas_pool_get_array(statePool, 'ssh', ssh, timeLevel)

      call mpas_pool_get_array(diagnosticsPool, 'zTop', zTop)
      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)

      !
      ! set pointers for fields related ocean state
      !
      call mpas_pool_get_array(diagnosticsPool, 'density', density)
      call mpas_pool_get_array(diagnosticsPool, 'displacedDensity', displacedDensity)
      call mpas_pool_get_array(diagnosticsPool, 'potentialDensity', potentialDensity)
      call mpas_pool_get_array(diagnosticsPool, 'bulkRichardsonNumber', bulkRichardsonNumber)
      call mpas_pool_get_array(diagnosticsPool, 'unresolvedShear', unresolvedShear)
      call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth', boundaryLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, 'RiTopOfCell', RiTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop',BruntVaisalaFreqTop)
      call mpas_pool_get_array(diagnosticsPool, 'bulkRichardsonNumberBuoy',bulkRichardsonNumberBuoy)
      call mpas_pool_get_array(diagnosticsPool, 'bulkRichardsonNumberShear',bulkRichardsonNumberShear)
      call mpas_pool_get_array(diagnosticsPool, 'indexBoundaryLayerDepth',indexBoundaryLayerDepth)

      !
      ! set pointers for fields related to ocean forcing state
      !
      call mpas_pool_get_array(forcingPool, 'iceFraction', iceFraction)

      !
      ! set pointers for fields related forcing at ocean surface
      !
      call mpas_pool_get_array(diagnosticsPool, 'surfaceFrictionVelocity', surfaceFrictionVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'surfaceBuoyancyForcing', surfaceBuoyancyForcing)

      !
      ! set pointers for viscosity/diffusivity and intialize to zero
      !
      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfCell', vertViscTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)


      !
      ! set pointers for nonlocal flux and intialize to zero
      !
      call mpas_pool_get_array(diagnosticsPool, 'vertNonLocalFlux', vertNonLocalFlux)
      call mpas_pool_get_dimension(diagnosticsPool, 'index_vertNonLocalFluxTemp', index_vertNonLocalFluxTemp)

      nCells = nCellsArray( size(nCellsArray) )

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         vertViscTopOfCell(:, iCell) = 0.0_RKIND
         vertDiffTopOfCell(:, iCell) = 0.0_RKIND
         vertNonLocalFlux(:, :, iCell) = 0.0_RKIND
      end do
      !$omp end do

      nCells = nCellsArray( 2 )

      !
      ! start by adding the mininum background values to the visocity/diffusivity arrays
      !
      if (cvmixBackgroundOn) then
         !$omp do schedule(runtime)
         do iCell = 1, nCells
            do k = 1, nVertLevelsP1
               vertViscTopOfCell(k, iCell) = vertViscTopOfCell(k, iCell) + config_cvmix_background_viscosity
               vertDiffTopOfCell(k, iCell) = vertDiffTopOfCell(k, iCell) + config_cvmix_background_diffusion
            end do
         end do
         !$omp end do
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

      allocate(normalVelocitySum(nVertLevels))
      allocate(potentialDensitySum(nVertLevels))
      allocate(surfaceAverageIndex(nVertLevels))
      allocate(deltaVelocitySquared(nVertLevels))

      do k = 1, nVertLevels
         Nsqr_iface(k) = 0.0_RKIND
         turbulentScalarVelocityScale(k) = 0.0_RKIND
      end do
      Nsqr_iface(nVertLevelsP1) = 0.0_RKIND

      call mpas_timer_start('cvmix cell loop', .false.)
      !$omp do schedule(runtime) private(k, bulkRichardsonNumberStop, kIndexOBL, bulkRichardsonFlag)
      do iCell = 1, nCells
         invAreaCell = 1.0_RKIND / areaCell(iCell)

         ! specify geometry/location
         cvmix_variables % SeaSurfaceHeight = ssh(iCell)
         cvmix_variables % Coriolis = fCell(iCell)
         cvmix_variables % lat = latCell(iCell) * 180.0_RKIND / 3.14_RKIND
         cvmix_variables % lon = lonCell(iCell) * 180.0_RKIND / 3.14_RKIND

         ! fill vertical position of column
         ! CVMix assume top of ocean is at z=0, so building all z-coordinate data based on layerThickness
         cvmix_variables % zw_iface(1) = 0.0_RKIND
         cvmix_variables % dzw(1) = layerThickness(1,iCell)/2.0_RKIND
         cvmix_variables % zt_cntr(1) = -layerThickness(1,iCell)/2.0_RKIND
         do k=2,maxLevelCell(iCell)
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

         ! fill Ri
         RiSmoothed(1:maxLevelCell(iCell)) = RiTopOfCell(1:maxLevelCell(iCell),iCell)
         RiSmoothed(maxLevelCell(iCell)+1) = RiSmoothed(maxLevelCell(iCell))
         do k=2,maxLevelCell(iCell)
            ! For input to cvmix, Richardson number should be positive.
             RiSmoothed(k)=max(0.0_RKIND, RiSmoothed(k))
         enddo
         cvmix_variables%ShearRichardson_iface => RiSmoothed(1:maxLevelCell(iCell)+1)

         ! fill BVF
         BVFSmoothed(1:maxLevelCell(iCell)) = max(0.0_RKIND,BruntVaisalaFreqTop(1:maxLevelCell(iCell),iCell))
         BVFSmoothed(maxLevelCell(iCell)+1) = max(0.0_RKIND,BVFSmoothed(maxLevelCell(iCell)))
         cvmix_variables%SqrBuoyancyFreq_iface => BVFSmoothed(1:maxLevelCell(iCell)+1)

         ! fill the intent(in) KPP
         cvmix_variables % SurfaceFriction = surfaceFrictionVelocity(iCell)
         cvmix_variables % SurfaceBuoyancyForcing = surfaceBuoyancyForcing(iCell)
         cvmix_variables % BulkRichardson_cntr => bulkRichardsonNumber(1:maxLevelCell(iCell), iCell)

         if (config_use_cvmix_shear) then

            do k = 1, maxLevelCell(iCell) + 1
               cvmix_variables % Mdiff_iface(k) = 0.0_RKIND
               cvmix_variables % Tdiff_iface(k) = 0.0_RKIND
            end do
            call cvmix_coeffs_shear( &
                 cvmix_variables, &
                 cvmix_shear_params)
            ! add shear mixing to vertical viscosity/diffusivity
            ! at present, shear mixing adds in background values when using PP, but background is
            ! accounted for seperately. so remove bac    kground from shear mixing values

            if(config_cvmix_shear_mixing_scheme=='PP') then
               do k = 1, maxLevelCell(iCell) + 1
                  cvmix_variables % Mdiff_iface(k) = cvmix_variables % Mdiff_iface(k) - config_cvmix_background_viscosity
                  cvmix_variables % Tdiff_iface(k) = cvmix_variables % Tdiff_iface(k) - config_cvmix_background_diffusion
               end do
            endif

            do k = 1, maxLevelCell(iCell)
               vertViscTopOfCell(k, iCell) = vertViscTopOfCell(k, iCell) + cvmix_variables % Mdiff_iface(k)
               vertDiffTopOfCell(k, iCell) = vertDiffTopOfCell(k, iCell) + cvmix_variables % Tdiff_iface(k)
            end do

         endif ! if (config_use_cvmix_shear)
         ! call kpp ocean mixed layer scheme
         if (cvmixKPPOn) then

            ! copy data into cvmix_variables
            do k = 1, maxLevelCell(iCell) + 1
              cvmix_variables % Mdiff_iface(k) = vertViscTopOfCell(k, iCell)
              cvmix_variables % Tdiff_iface(k) = vertDiffTopOfCell(k, iCell)
            end do

            if (config_use_cvmix_fixed_boundary_layer) then
               cvmix_variables % BoundaryLayerDepth = config_cvmix_kpp_boundary_layer_depth
               cvmix_variables % kOBL_depth = cvmix_kpp_compute_kOBL_depth(  &
                          zw_iface = cvmix_variables%zw_iface(1:maxLevelCell(iCell)+1), &
                          zt_cntr = cvmix_variables%zt_cntr(1:maxLevelCell(iCell)),    &
                          OBL_depth = cvmix_variables % BoundaryLayerDepth )

            else

              ! set stratification
              do k=1,maxLevelCell(iCell)
                  Nsqr_iface(k) = BVFSmoothed(k)
              enddo
              k=min(maxLevelCell(iCell)+1,nVertLevels)
              Nsqr_iface(k:maxLevelCell(iCell)+1) = Nsqr_iface(k-1)

              ! compute bulk Richardson number
              ! assume boundary layer depth is at bottom of every kIndexOBL cell
              bulkRichardsonNumberStop = config_cvmix_kpp_stop_OBL_search * config_cvmix_kpp_criticalBulkRichardsonNumber
              bulkRichardsonFlag = .false.
              topIndex = 1
!             call mpas_timer_start('Bulk Richardson kIndexOBL loops')
              ! compute the index of the last cell in the KPP defined surface layer
              ! this index is used for the necessary surface layer averages of buoyancy and momentum
              do kIndexOBL = 1, maxLevelCell(iCell)

                 ! Reset deltaVelocitySquared and bulkRichardsonNumber at this layer for the later computation
                 deltaVelocitySquared(kIndexOBL) = 0.0_RKIND
                 bulkRichardsonNumber(kIndexOBL, iCell) = bulkRichardsonNumberStop - 1.0_RKIND

                 ! set OBL at bottom of kIndexOBL cell for computation of bulk Richardson number
                 cvmix_variables % BoundaryLayerDepth = abs(cvmix_variables % zw_iface(kIndexOBL+1))
                 sigma = -cvmix_variables % zt_cntr(kIndexOBL) / cvmix_variables % BoundaryLayerDepth

                 OBLDepths(kIndexOBL) = abs(cvmix_variables % zw_iface(kIndexOBL+1))
                 interfaceForcings(kIndexOBL) = cvmix_variables % SurfaceBuoyancyForcing

                 ! initialize the surfaceAverageIndex for cases when the if statement below is not true
                 surfaceAverageIndex(kIndexOBL) = 1
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
                   OBL_depth = OBLDepths(1:maxLevelCell(iCell)), &
                   surf_buoy_force = interfaceForcings(1:maxLevelCell(iCell)), &
                   surf_fric_vel = cvmix_variables % SurfaceFriction, &
                   w_s = turbulentScalarVelocityScale(1:maxLevelCell(iCell)))

              ! averaging over a surface layer assuming that BLdepth is cell bottom
              ! Build deltaVelocitySquared
              do i = 1, nEdgesOnCell(iCell)
                 iEdge = edgesOnCell(i, iCell)
                 factor = 0.5_RKIND * dcEdge(iEdge) * dvEdge(iEdge) * invAreaCell

                 deltaVelocitySquared(1) = 0.0_RKIND

                 normalVelocitySum(1) = normalVelocity(1, iEdge)

                 do kIndexOBL = 2, maxLevelCell(iCell)
                    normalVelocitySum(kIndexOBL) = normalVelocitySum(kIndexOBL-1) + normalVelocity(kIndexOBL, iEdge)
                 end do

                 do kIndexOBL = 1, maxLevelCell(iCell)
                    normalVelocityAv = normalVelocitySum(surfaceAverageIndex(kIndexOBL )) / &
                        real( surfaceAverageIndex(kIndexOBL), kind=RKIND)
                    delU2 = ( normalVelocityAv - normalVelocity(kIndexOBL, iEdge) )**2
                    deltaVelocitySquared(kIndexOBL) = deltaVelocitySquared(kIndexOBL) + factor * delU2
                 end do
              end do

              potentialDensitySum(1) = potentialDensity(1, iCell)
              do kIndexOBL = 2, maxLevelCell(iCell)
                potentialDensitySum(kIndexOBL) = potentialDensitySum(kIndexOBL-1) + potentialDensity(kIndexOBL, iCell)
              end do

              do kIndexOBL = 1, maxLevelCell(iCell)
                ! !compute shear contribution assuming BLdepth is cell bottom
                bulkRichardsonNumberShear(kIndexOBL,iCell) = max(deltaVelocitySquared(kIndexOBL), 1.0e-15_RKIND)

                bulkRichardsonNumberBuoy(kIndexOBL,iCell) = gravity * (potentialDensity(kIndexOBL, iCell) &
                                  - potentialDensitySum(surfaceAverageIndex(kIndexOBL)) &
                                  / real(surfaceAverageIndex(kIndexOBL), kind=RKIND)) / rho_sw
              end do ! do kIndexOBL
!             call mpas_timer_stop('Bulk Richardson kIndexOBL loops')

              cvmix_variables % bulkRichardson_cntr(:) = cvmix_kpp_compute_bulk_Richardson( &
                   zt_cntr = cvmix_variables % zt_cntr(1:maxLevelCell(iCell)), &
                   delta_buoy_cntr = bulkRichardsonNumberBuoy(1:maxLevelCell(iCell),iCell), &
                   delta_Vsqr_cntr = bulkRichardsonNumberShear(1:maxLevelCell(iCell),iCell), &
                   ws_cntr = turbulentScalarVelocityScale(:), &
                   Nsqr_iface = Nsqr_iface(1:maxLevelCell(iCell)+1) )

              ! each level of bulk Richardson is computed as if OBL resided at bottom of that level

              call cvmix_kpp_compute_OBL_depth(  &
                   Ri_bulk = bulkRichardsonNumber(1:maxLevelCell(iCell),iCell), &
                   zw_iface = cvmix_variables % zw_iface(1:maxLevelCell(iCell)+1), &
                   OBL_depth = cvmix_variables % BoundaryLayerDepth, &
                   kOBL_depth = cvmix_variables % kOBL_depth, &
                   zt_cntr = cvmix_variables % zt_cntr(1:maxLevelCell(iCell)), &
                   surf_fric = cvmix_variables % SurfaceFriction, &
                   surf_buoy = cvmix_variables % SurfaceBuoyancyForcing, &
                   Coriolis = cvmix_variables % Coriolis)

             endif  ! if (config_use_cvmix_fixed_boundary_layer) then

             ! apply minimum limit to OBL
             if(cvmix_variables % BoundaryLayerDepth .lt. layerThickness(1,iCell)/2.0_RKIND) then
                cvmix_variables % BoundaryLayerDepth = layerThickness(1,iCell)/2.0_RKIND
                cvmix_variables % kOBL_depth = cvmix_kpp_compute_kOBL_depth(  &
                          zw_iface = cvmix_variables%zw_iface(1:maxLevelCell(iCell)+1),&
                          zt_cntr = cvmix_variables%zt_cntr(1:maxLevelCell(iCell)),    &
                          OBL_depth = cvmix_variables % BoundaryLayerDepth )
             endif

             ! apply minimum limit to OBL under sea-ice
             if(iceFraction(iCell).gt.0.15_RKIND) then
             if(cvmix_variables % BoundaryLayerDepth .lt. configure_cvmix_kpp_minimum_OBL_under_sea_ice) then
                cvmix_variables % BoundaryLayerDepth = configure_cvmix_kpp_minimum_OBL_under_sea_ice
                cvmix_variables % kOBL_depth = cvmix_kpp_compute_kOBL_depth(  &
                          zw_iface = cvmix_variables%zw_iface(1:maxLevelCell(iCell)+1),&
                          zt_cntr = cvmix_variables%zt_cntr(1:maxLevelCell(iCell)),    &
                          OBL_depth = cvmix_variables % BoundaryLayerDepth )
             endif
             endif

             ! apply maximum limit to OBL
             if(cvmix_variables % BoundaryLayerDepth .gt. abs(cvmix_variables%zt_cntr(maxLevelCell(iCell)))) then
                cvmix_variables % BoundaryLayerDepth = abs(cvmix_variables%zt_cntr(maxLevelCell(iCell)))
                cvmix_variables % kOBL_depth = cvmix_kpp_compute_kOBL_depth(  &
                          zw_iface = cvmix_variables%zw_iface(1:maxLevelCell(iCell)+1), &
                          zt_cntr = cvmix_variables%zt_cntr(1:maxLevelCell(iCell)),    &
                          OBL_depth = cvmix_variables % BoundaryLayerDepth )

             endif

!           call mpas_timer_start('cvmix coeffs kpp', .false.)
            call cvmix_coeffs_kpp(                                              &
                          Mdiff_out = cvmix_variables % Mdiff_iface(1:maxLevelCell(iCell)+1),       &
                          Tdiff_out = cvmix_variables % Tdiff_iface(1:maxLevelCell(iCell)+1),       &
                          Sdiff_out = cvmix_variables % Sdiff_iface(1:maxLevelCell(iCell)+1),       &
                          zw = cvmix_variables%zw_iface(1:maxLevelCell(iCell)+1),            &
                          zt = cvmix_variables%zt_cntr(1:maxLevelCell(iCell)),               &
                          old_Mdiff = cvmix_variables%Mdiff_iface(1:maxLevelCell(iCell)+1),         &
                          old_Tdiff = cvmix_variables%Tdiff_iface(1:maxLevelCell(iCell)+1),         &
                          old_Sdiff = cvmix_variables%Sdiff_iface(1:maxLevelCell(iCell)+1),         &
                          OBL_depth = cvmix_variables%BoundaryLayerDepth,                   &
                          kOBL_depth = cvmix_variables%kOBL_depth,                           &
                          Tnonlocal = cvmix_variables%kpp_Tnonlocal_iface(1:maxLevelCell(iCell)+1), &
                          Snonlocal = cvmix_variables%kpp_Snonlocal_iface(1:maxLevelCell(iCell)+1), &
                          surf_fric = cvmix_variables%SurfaceFriction,                      &
                          surf_buoy = cvmix_variables%SurfaceBuoyancyForcing,               &
                          nlev = maxLevelCell(iCell),                                  &
                          max_nlev = nVertLevels)
!           call mpas_timer_stop('cvmix coeffs kpp')

            ! intent out of BoundaryLayerDepth is boundary layer depth measured in meters and vertical index
            boundaryLayerDepth(iCell) = cvmix_variables % BoundaryLayerDepth
            indexBoundaryLayerDepth(iCell) = cvmix_variables % kOBL_depth



            if(config_cvmix_kpp_matching .eq. 'SimpleShapes') then
                do k = 1, int(indexBoundaryLayerDepth(iCell))
                   vertViscTopOfCell(k,iCell) = vertViscTopOfCell(k,iCell) + cvmix_variables % Mdiff_iface(k)
                   vertDiffTopOfCell(k,iCell) = vertDiffTopOfCell(k,iCell) + cvmix_variables % Tdiff_iface(k)
                end do
                do k = int(indexBoundaryLayerDepth(iCell))+1, maxLevelCell(iCell)+1
                   vertViscTopOfCell(k,iCell) = cvmix_variables % Mdiff_iface(k)
                   vertDiffTopOfCell(k,iCell) = cvmix_variables % Tdiff_iface(k)
                enddo
            else
                do k = 1, maxLevelCell(iCell) + 1
                   vertViscTopOfCell(k, iCell) = cvmix_variables % Mdiff_iface(k)
                   vertDiffTopOfCell(k, iCell) = cvmix_variables % Tdiff_iface(k)
                end do
            endif

            ! store non-local flux terms
            ! these flux terms must be multiplied by the surfaceTracerFlux field
            ! the tracer tendency is then the vertical divergence of vertNonLocalFlux*surfaceTracerFlux
            ! both of these operations are done in ocn_tracer_nonlocalflux_tend routine
            vertNonLocalFlux(index_vertNonLocalFluxTemp,:,iCell) = cvmix_variables % kpp_Tnonlocal_iface(:)

         endif !if (config_use_cvmix_kpp)

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
               do k = int(indexBoundaryLayerDepth(iCell)) + 1, maxLevelCell(iCell)
                  vertViscTopOfCell(k,iCell) = vertViscTopOfCell(k,iCell) + cvmix_variables % Mdiff_iface(k)
                  vertDiffTopOfCell(k,iCell) = vertDiffTopOfCell(k,iCell) + cvmix_variables % Tdiff_iface(k)
               enddo
            else
               do k = 1, maxLevelCell(iCell) + 1
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
         vertViscTopOfCell(1, iCell) = 0.0_RKIND
         vertDiffTopOfCell(1, iCell) = 0.0_RKIND
         do k = maxLevelCell(iCell)+1, nVertLevelsP1
            vertViscTopOfCell(k, iCell)=0.0_RKIND
            vertDiffTopOfCell(k, iCell)=0.0_RKIND
         end do

      end do  ! do iCell=1,mesh%nCells
      !$omp end do
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
      deallocate(BVFSmoothed)
      deallocate(normalVelocitySum)
      deallocate(potentialDensitySum)
      deallocate(surfaceAverageIndex)
      deallocate(deltaVelocitySquared)

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
      type (block_type), pointer :: block

      ! CVMix
      logical, pointer :: config_use_cvmix

      ! background
      logical, pointer :: config_use_cvmix_background
      real (kind=RKIND), pointer :: config_cvmix_background_viscosity, config_cvmix_background_diffusion
      real (kind=RKIND), pointer :: config_cvmix_prandtl_number

      ! Shear configs
      logical, pointer :: config_use_cvmix_shear
      character (len=StrKIND), pointer :: config_cvmix_shear_mixing_scheme
      real (kind=RKIND), pointer :: config_cvmix_shear_PP_nu_zero, config_cvmix_shear_PP_alpha, config_cvmix_shear_PP_exp, &
                                    config_cvmix_shear_KPP_nu_zero, config_cvmix_shear_KPP_Ri_zero, config_cvmix_shear_KPP_exp

      ! Convection configs
      logical, pointer ::  config_use_cvmix_convection
      real (kind=RKIND), pointer :: config_cvmix_convective_diffusion, config_cvmix_convective_viscosity, &
                                    config_cvmix_convective_triggerBVF
      logical, pointer :: config_cvmix_convective_basedOnBVF, config_cvmix_kpp_use_enhanced_diff

      ! Tidal mixing
      logical, pointer :: config_use_cvmix_tidal_mixing

      ! Double diffusion
      logical, pointer :: config_use_cvmix_double_diffusion

      ! KPP configs
      logical, pointer :: config_use_cvmix_kpp
      character (len=StrKIND), pointer :: config_cvmix_kpp_matching, config_cvmix_kpp_interpolationOMLType
      logical, pointer :: config_cvmix_kpp_EkmanOBL, config_cvmix_kpp_MonObOBL
      real (kind=RKIND), pointer :: config_cvmix_kpp_criticalBulkRichardsonNumber, &
                                    config_cvmix_kpp_surface_layer_extent, &
                                    config_cvmix_kpp_stop_OBL_search
      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix', config_use_cvmix)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_background', config_use_cvmix_background)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_background_viscosity', config_cvmix_background_viscosity)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_background_diffusion', config_cvmix_background_diffusion)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_shear', config_use_cvmix_shear)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_shear_mixing_scheme', config_cvmix_shear_mixing_scheme)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_shear_PP_nu_zero', config_cvmix_shear_PP_nu_zero)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_shear_PP_alpha', config_cvmix_shear_PP_alpha)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_shear_PP_exp', config_cvmix_shear_PP_exp)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_shear_KPP_nu_zero', config_cvmix_shear_KPP_nu_zero)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_shear_KPP_Ri_zero', config_cvmix_shear_KPP_Ri_zero)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_shear_KPP_exp', config_cvmix_shear_KPP_exp)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_convection', config_use_cvmix_convection)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_convective_basedOnBVF', config_cvmix_convective_basedOnBVF)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_convective_triggerBVF', config_cvmix_convective_triggerBVF)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_tidal_mixing', config_use_cvmix_tidal_mixing)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_double_diffusion', config_use_cvmix_double_diffusion)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_criticalBulkRichardsonNumber', &
                                config_cvmix_kpp_criticalBulkRichardsonNumber)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_stop_OBL_search', config_cvmix_kpp_stop_OBL_search)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_interpolationOMLType', config_cvmix_kpp_interpolationOMLType)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_interpolationOMLType', config_cvmix_kpp_interpolationOMLType)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_EkmanOBL', config_cvmix_kpp_EkmanOBL)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_MonObOBL', config_cvmix_kpp_MonObOBL)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_matching', config_cvmix_kpp_matching)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_surface_layer_extent', config_cvmix_kpp_surface_layer_extent)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_prandtl_number', config_cvmix_prandtl_number)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_convective_diffusion', config_cvmix_convective_diffusion)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_convective_viscosity', config_cvmix_convective_viscosity)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_use_enhanced_diff', config_cvmix_kpp_use_enhanced_diff)
      cvmixOn = config_use_cvmix
      cvmixBackgroundOn = config_use_cvmix_background
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

      !
      ! initialize background mixing
      !
      if (config_use_cvmix_background .or. config_use_cvmix_shear) then
        call cvmix_init_bkgnd( &
               bkgnd_Tdiff = config_cvmix_background_diffusion, &
               bkgnd_Mdiff = config_cvmix_background_viscosity, &
               CVmix_bkgnd_params_user = cvmix_background_params)
      endif

      !
      ! initialize shear-based mixing
      !
      if (config_use_cvmix_shear) then
        if (.not. config_use_cvmix_background .and. trim(config_cvmix_shear_mixing_scheme) == 'PP') then
            write(stderrUnit, *) "ERROR: config_use_cvmix_shear cannot be used with with config_cvmix_shear_mixing_scheme = 'PP'"
            write(stderrUnit, *) "       without config_use_cvmix_background = .true."
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
            write(stderrUnit, *) "ERROR: config_cvmix_convective_basedOnBVF = .false. is not supported. Change to true."
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
           write(stderrUnit,*) "WARNING: use of option MatchBoth is discouraged, use SimpleShapes instead"
        elseif(.not. config_cvmix_kpp_matching.eq."SimpleShapes") then
           write(stderrUnit,*) "ERROR: unknown value for config_cvmix_kpp_matching., supported values are:"
           write(stderrUnit,*) "         SimpleShapes or MatchBoth"
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
               lenhanced_diff = config_cvmix_kpp_use_enhanced_diff)
      endif


   !--------------------------------------------------------------------

   end subroutine ocn_vmix_cvmix_init!}}}

!***********************************************************************

end module ocn_vmix_cvmix

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
