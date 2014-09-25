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

   use mpas_grid_types
   use mpas_timer
   use mpas_io_units

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
   type(cvmix_data_type) :: cvmix_variables

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

   subroutine ocn_vmix_coefs_cvmix_build(meshPool, statePool, diagnosticsPool, err, timeLevelIn)!{{{

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

      integer, dimension(:), pointer :: &
        maxLevelCell

      real (kind=RKIND), dimension(:), pointer :: &
        latCell, lonCell, bottomDepth, surfaceBuoyancyForcing, surfaceFrictionVelocity, fCell, &
        boundaryLayerDepth, ssh, indexBoundaryLayerDepth
        
      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfCell, vertDiffTopOfCell, layerThickness, &
        zMid, zTop, density, displacedDensity, potentialDensity, &
        bulkRichardsonNumber, RiTopOfCell, BruntVaisalaFreqTop, &
        bulkRichardsonNumberBuoy, bulkRichardsonNumberShear

      real (kind=RKIND), dimension(:,:,:), pointer :: vertNonLocalFlux
      integer, pointer :: index_vertNonLocalFluxTemp

      integer, pointer :: config_cvmix_kpp_niterate
      logical, pointer :: config_use_cvmix_shear, config_use_cvmix_convection, config_use_cvmix_kpp
      character (len=StrKIND), pointer :: config_cvmix_shear_mixing_scheme, config_cvmix_kpp_matching

      integer :: k, iCell, jCell, iNeighbor, iter, timeLevel
      integer, pointer :: nVertLevels, nCells
      real (kind=RKIND) :: r, layerSum
      real (kind=RKIND), dimension(:), allocatable :: sigma, Nsqr_iface, turbulentScalarVelocityScale, tmp
      real (kind=RKIND), dimension(:), allocatable, target :: RiSmoothed, BVFSmoothed

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

!     write(stdoutUnit,*) 'TDR: ocn_vmix_coefs_cvmix_build enter'
      !
      ! only build up viscosity/diffusivity if CVMix is turned on
      !
      if ( .not. cvmixOn ) return

      !
      ! set parameters
      !
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_background_viscosity', config_cvmix_background_viscosity)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_background_diffusion', config_cvmix_background_diffusion)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_niterate', config_cvmix_kpp_niterate)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_shear', config_use_cvmix_shear)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_convection', config_use_cvmix_convection)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_shear_mixing_scheme', config_cvmix_shear_mixing_scheme)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_matching', config_cvmix_kpp_matching)

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

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
      call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth', boundaryLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, 'RiTopOfCell', RiTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'BruntVaisalaFreqTop',BruntVaisalaFreqTop)
      call mpas_pool_get_array(diagnosticsPool, 'bulkRichardsonNumberBuoy',bulkRichardsonNumberBuoy)
      call mpas_pool_get_array(diagnosticsPool, 'bulkRichardsonNumberShear',bulkRichardsonNumberShear)
      call mpas_pool_get_array(diagnosticsPool, 'indexBoundaryLayerDepth',indexBoundaryLayerDepth)

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
 
      vertViscTopOfCell = 0.0
      vertDiffTopOfCell = 0.0

      !
      ! set pointers for nonlocal flux and intialize to zero
      !
      call mpas_pool_get_array(diagnosticsPool, 'vertNonLocalFlux', vertNonLocalFlux)
      call mpas_pool_get_dimension(diagnosticsPool, 'index_vertNonLocalFluxTemp', index_vertNonLocalFluxTemp)

      vertNonLocalFlux = 0.0

      !
      ! start by adding the mininum background values to the visocity/diffusivity arrays
      !
      if (cvmixBackgroundOn) then
        vertViscTopOfCell(:,:) = vertViscTopOfCell(:,:) + config_cvmix_background_viscosity
        vertDiffTopOfCell(:,:) = vertDiffTopOfCell(:,:) + config_cvmix_background_diffusion
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
      allocate(cvmix_variables % BulkRichardson_cntr(nVertLevels))

      allocate(sigma(nVertLevels))
      allocate(Nsqr_iface(nVertLevels+1))
      allocate(turbulentScalarVelocityScale(nVertLevels))
      allocate(tmp(nVertLevels+1))
      allocate(RiSmoothed(nVertLevels+1))
      allocate(BVFSmoothed(nVertLevels+1))

      do iCell = 1, nCells

        ! specify geometry/location
        cvmix_variables % SeaSurfaceHeight = ssh(iCell)
        cvmix_variables % Coriolis = fCell(iCell)
        cvmix_variables % lat = latCell(iCell) * 180.0 / 3.14
        cvmix_variables % lon = lonCell(iCell) * 180.0 / 3.14

        ! fill vertical position of column
        ! CVMix assume top of ocean is at z=0, so building all z-coordinate data based on layerThickness
        cvmix_variables % zw_iface(1) = 0.0
        cvmix_variables % dzw(1) = layerThickness(1,iCell)/2.0
        cvmix_variables % zt_cntr(1) = -layerThickness(1,iCell)/2.0
        do k=2,maxLevelCell(iCell)
          cvmix_variables % zw_iface(k) = cvmix_variables % zw_iface(k-1) - layerThickness(k-1,iCell)
          cvmix_variables % zt_cntr(k) = cvmix_variables %  zw_iface(k) - layerThickness(k,iCell)/2.0
          cvmix_variables % dzw(k) = cvmix_variables % zt_cntr(k-1) - cvmix_variables % zt_cntr(k) 
          cvmix_variables % dzt(k) = layerThickness(k,iCell)
        enddo
        k = maxLevelCell(iCell)+1
        cvmix_variables % zw_iface(k) = cvmix_variables % zw_iface(k-1) - layerThickness(k-1,iCell)
        cvmix_variables % dzw(k) = cvmix_variables % zt_cntr(k-1) - cvmix_variables % zw_iface(k)
        do k = maxLevelCell(iCell) + 1, nVertLevels
          cvmix_variables % zw_iface(k+1) = cvmix_variables % zw_iface(maxLevelCell(iCell)+1)
          cvmix_variables % zt_cntr(k) = cvmix_variables % zw_iface(maxLevelCell(iCell)+1)
          cvmix_variables % dzw(k+1) = 0.0
          cvmix_variables % dzt(k) = 0.0
        enddo

        ! fill the intent(in) convective adjustment
        cvmix_variables % nlev = maxLevelCell(iCell)
        cvmix_variables % OceanDepth = bottomDepth(iCell)
        cvmix_variables % WaterDensity_cntr => density(:,iCell)
        cvmix_variables % AdiabWaterDensity_cntr => displacedDensity(:,iCell)

        ! eliminate 2dz mode from Ri
        RiSmoothed(1:nVertLevels) = RiTopOfCell(1:nVertLevels,iCell)
        RiSmoothed(nVertLevels+1) = RiSmoothed(nVertLevels)
        do k=2,maxLevelCell(iCell)
          tmp(k) = (RiSmoothed(k-1)+2*RiSmoothed(k)+RiSmoothed(k+1))/4.0
        enddo
        k=1
        tmp(k)=tmp(k+1)
        k=maxLevelCell(iCell)+1
        tmp(k:nVertLevels+1)=tmp(k-1)
        RiSmoothed(:) = tmp(:)
        cvmix_variables%ShearRichardson_iface => RiSmoothed

        ! eliminate 2dz mode from BVF
        BVFSmoothed(1:nVertLevels) = BruntVaisalaFreqTop(1:nVertLevels,iCell)
        BVFSmoothed(nVertLevels+1) = BVFSmoothed(nVertLevels)
        do k=2,maxLevelCell(iCell)
          tmp(k) = (BVFSmoothed(k-1)+2*BVFSmoothed(k)+BVFSmoothed(k+1))/4.0
        enddo
        k=1
        tmp(k)=tmp(k+1)
        k=maxLevelCell(iCell)+1
        tmp(k:nVertLevels+1)=tmp(k-1)
        BVFSmoothed(:) = tmp(:)
        cvmix_variables%SqrBuoyancyFreq_iface => BVFSmoothed

        ! fill the intent(in) KPP
        cvmix_variables % SurfaceFriction = surfaceFrictionVelocity(iCell)
        cvmix_variables % SurfaceBuoyancyForcing = surfaceBuoyancyForcing(iCell)

        ! call shear-based mixing scheme
        if (config_use_cvmix_shear) then
          cvmix_variables % Mdiff_iface(:)=0.0
          cvmix_variables % Tdiff_iface(:)=0.0
          call cvmix_coeffs_shear( &
                 cvmix_variables, &
                 cvmix_background_params, &
                 1, &
                 .false., &
                 cvmix_shear_params)

        ! add shear mixing to vertical viscosity/diffusivity
        ! at present, shear mixing adds in background values when using PP, but background is accounted for seperately. so remove background from shear mixing values
        if(config_cvmix_shear_mixing_scheme=='PP') then
          vertViscTopOfCell(:,iCell) = vertViscTopOfCell(:,iCell) + cvmix_variables % Mdiff_iface(:) - config_cvmix_background_viscosity
          vertDiffTopOfCell(:,iCell) = vertDiffTopOfCell(:,iCell) + cvmix_variables % Tdiff_iface(:) - config_cvmix_background_diffusion
        endif

        endif ! if (config_use_cvmix_shear)

        !
        ! put tidal mixing here
        !

        !
        ! put double diffusion mxing here
        !

        ! call kpp ocean mixed layer scheme
        if (cvmixKPPOn) then

        ! copy data into cvmix_variables, then iterate
        cvmix_variables % Mdiff_iface(:)=vertViscTopOfCell(:,iCell)
        cvmix_variables % Tdiff_iface(:)=vertDiffTopOfCell(:,iCell)
        cvmix_variables % BoundaryLayerDepth = boundaryLayerDepth(iCell)
        cvmix_variables % kOBL_depth = cvmix_kpp_compute_kOBL_depth(cvmix_variables%zw_iface, cvmix_variables%zt_cntr, cvmix_variables%BoundaryLayerDepth)

        do iter=1,config_cvmix_kpp_niterate

          if(cvmix_variables % BoundaryLayerDepth .lt. layerThickness(1,iCell)/2.0) then
             cvmix_variables % BoundaryLayerDepth = layerThickness(1,iCell)/2.0 + 0.01
          endif

          if(cvmix_variables % BoundaryLayerDepth .gt. abs(cvmix_variables%zw_iface(maxLevelCell(iCell)+1))) then
             cvmix_variables % BoundaryLayerDepth = abs(cvmix_variables%zt_cntr(maxLevelCell(iCell))) + 0.01
          endif

          ! compute ocean boundary layer depth 
          do k=1,maxLevelCell(iCell)
            sigma(k) = -cvmix_variables % zt_cntr(k) / cvmix_variables % BoundaryLayerDepth
            Nsqr_iface(k) = BVFSmoothed(k)
          enddo
          k=maxLevelCell(iCell)+1
          sigma(k:nVertLevels) = sigma(k-1)
          Nsqr_iface(k:nVertLevels+1) = Nsqr_iface(k-1)

          ! eliminate 2dz mode from N2
          tmp(:)=0.0
          do k=2,maxLevelCell(iCell)
            tmp(k) = (Nsqr_iface(k-1)+2*Nsqr_iface(k)+Nsqr_iface(k+1))/4.0
          enddo
          k=maxLevelCell(iCell)+1
          tmp(1)=tmp(2)
          tmp(k:nVertLevels+1)=tmp(k-1)
          Nsqr_iface(:)=tmp(:)

          ! compute the turbulent scales in order to compute the bulk Richardson number
          call cvmix_kpp_compute_turbulent_scales( &
            sigma(:), &
            cvmix_variables % BoundaryLayerDepth, &
            cvmix_variables % SurfaceBuoyancyForcing, &
            cvmix_variables % SurfaceFriction, &
            w_s = turbulentScalarVelocityScale(:) )

          cvmix_variables % BulkRichardson_cntr  = cvmix_kpp_compute_bulk_Richardson( &
            cvmix_variables % zt_cntr(1:nVertLevels), & 
            bulkRichardsonNumberBuoy(1:nVertLevels,iCell), &
            bulkRichardsonNumberShear(1:nVertLevels,iCell), &
            ws_cntr = turbulentScalarVelocityScale(:), &
            Nsqr_iface = Nsqr_iface(1:nVertLevels+1) )

          ! compute the boundary layer depth based on model state at current time step
          call cvmix_kpp_compute_OBL_depth( CVmix_vars = cvmix_variables)

          ! given current estimate of OBL and vertical profile of visc/diff, compute boundary layer mixing
          call cvmix_coeffs_kpp( CVmix_vars = cvmix_variables )

        end do ! iterate

        ! intent out of BoundaryLayerDepth is boundary layer depth measured in meters and vertical index
        boundaryLayerDepth(iCell) = cvmix_variables % BoundaryLayerDepth
        indexBoundaryLayerDepth(iCell) = cvmix_variables % kOBL_depth
        bulkRichardsonNumber(:,iCell) = cvmix_variables % BulkRichardson_cntr(:)

        ! if using KPP with "MatchBoth"   matching, then the output from KPP is the full viscosity/diffusivity
        ! if using KPP with "SimpleShape" matching, then the output from KPP needs to be added to current viscosity/diffusivity
        if(config_cvmix_kpp_matching.eq."MatchBoth") then
          vertViscTopOfCell(:,iCell) = cvmix_variables % Mdiff_iface(:)
          vertDiffTopOfCell(:,iCell) = cvmix_variables % Tdiff_iface(:)
        elseif(config_cvmix_kpp_matching.eq."SimpleShapes") then
          vertViscTopOfCell(:,iCell) = vertViscTopOfCell(:,iCell) + cvmix_variables % Mdiff_iface(:)
          vertDiffTopOfCell(:,iCell) = vertDiffTopOfCell(:,iCell) + cvmix_variables % Tdiff_iface(:)
        else
          stop
        endif

        ! store non-local flux terms
        ! these flux terms must be multiplied by the surfaceTracerFlux field
        ! the tracer tendency is then the vertical divergence of vertNonLocalFlux*surfaceTracerFlux
        ! both of these operations are done in ocn_tracer_nonlocalflux_tend routine
        vertNonLocalFlux(index_vertNonLocalFluxTemp,:,iCell) = cvmix_variables % kpp_Tnonlocal_iface(:)

        endif !if (config_use_cvmix_kpp)

        ! call convective mixing scheme
        if (config_use_cvmix_convection) then
          cvmix_variables % Mdiff_iface(:)=0.0
          cvmix_variables % Tdiff_iface(:)=0.0
          call cvmix_coeffs_conv( CVmix_vars = cvmix_variables )

          ! add convective mixing to vertical viscosity/diffusivity
          ! if using KPP, then do not apply convective mixing within the ocean boundary layer
          if(config_use_cvmix_kpp) then
            do k = int(indexBoundaryLayerDepth(iCell)) + 1, maxLevelCell(iCell)
              vertViscTopOfCell(k,iCell) = vertViscTopOfCell(k,iCell) + cvmix_variables % Mdiff_iface(k)
              vertDiffTopOfCell(k,iCell) = vertDiffTopOfCell(k,iCell) + cvmix_variables % Tdiff_iface(k)
            enddo
          else
            vertViscTopOfCell(:,iCell) = vertViscTopOfCell(:,iCell) + cvmix_variables % Mdiff_iface(:)
            vertDiffTopOfCell(:,iCell) = vertDiffTopOfCell(:,iCell) + cvmix_variables % Tdiff_iface(:)
          endif
        endif  ! if (config_use_cvmix_convection)

        ! computation of viscosity/diffusivity complete
        ! impose no-flux boundary conditions at top and bottom by zero viscosity/diffusivity
        vertViscTopOfCell(1,iCell) = 0.0
        vertDiffTopOfCell(1,iCell) = 0.0
        vertViscTopOfCell(maxLevelCell(iCell)+1:nVertLevels,iCell)=0.0
        vertDiffTopOfCell(maxLevelCell(iCell)+1:nVertLevels,iCell)=0.0
  
      end do  ! do iCell=1,mesh%nCells

      ! dellocate cmvix variables
      deallocate(cvmix_variables % Mdiff_iface)
      deallocate(cvmix_variables % Tdiff_iface)
      deallocate(cvmix_variables % zw_iface)
      deallocate(cvmix_variables % dzw)
      deallocate(cvmix_variables % zt_cntr)
      deallocate(cvmix_variables % dzt)
      deallocate(cvmix_variables % kpp_Tnonlocal_iface)
      deallocate(cvmix_variables % BulkRichardson_cntr)
  
      deallocate(sigma)
      deallocate(Nsqr_iface)
      deallocate(turbulentScalarVelocityScale)
      deallocate(tmp)
      deallocate(RiSmoothed)
      deallocate(BVFSmoothed)

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
      real (kind=RKIND), pointer :: config_cvmix_convective_diffusion, config_cvmix_convective_viscosity, config_cvmix_convective_triggerBVF
      logical, pointer :: config_cvmix_convective_basedOnBVF

      ! Tidal mixing
      logical, pointer :: config_use_cvmix_tidal_mixing

      ! Double diffusion
      logical, pointer :: config_use_cvmix_double_diffusion

      ! KPP configs
      logical, pointer :: config_use_cvmix_kpp
      character (len=StrKIND), pointer :: config_cvmix_kpp_matching, config_cvmix_kpp_interpolationOMLType
      logical, pointer :: config_cvmix_kpp_EkmanOBL, config_cvmix_kpp_MonObOBL
      real (kind=RKIND), pointer :: config_cvmix_kpp_criticalBulkRichardsonNumber, &
                                    config_cvmix_kpp_surface_layer_extent
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
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_criticalBulkRichardsonNumber', config_cvmix_kpp_criticalBulkRichardsonNumber)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_interpolationOMLType', config_cvmix_kpp_interpolationOMLType)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_interpolationOMLType', config_cvmix_kpp_interpolationOMLType)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_EkmanOBL', config_cvmix_kpp_EkmanOBL)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_MonObOBL', config_cvmix_kpp_MonObOBL)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_matching', config_cvmix_kpp_matching)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_kpp_surface_layer_extent', config_cvmix_kpp_surface_layer_extent)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_prandtl_number', config_cvmix_prandtl_number)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_convective_diffusion', config_cvmix_convective_diffusion)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_convective_viscosity', config_cvmix_convective_viscosity)

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
               bkgnd_diff = config_cvmix_background_diffusion, &
               bkgnd_visc = config_cvmix_background_viscosity, &
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
        call cvmix_init_kpp ( &
               ri_crit = config_cvmix_kpp_criticalBulkRichardsonNumber, &
               interp_type = config_cvmix_kpp_interpolationOMLType, &
               interp_type2 = config_cvmix_kpp_interpolationOMLType, &
               lEkman = config_cvmix_kpp_EkmanOBL, &
               lMonOb = config_cvmix_kpp_MonObOBL, &
               MatchTechnique = config_cvmix_kpp_matching, &
               surf_layer_ext = config_cvmix_kpp_surface_layer_extent) 
      endif


   !--------------------------------------------------------------------

   end subroutine ocn_vmix_cvmix_init!}}}

!***********************************************************************

end module ocn_vmix_cvmix

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
