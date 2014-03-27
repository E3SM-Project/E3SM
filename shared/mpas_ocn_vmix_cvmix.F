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
   type(cvmix_conv_params_type) :: cvmix_conv_params
   type(cvmix_kpp_params_type) :: cvmix_kpp_params
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
        maxLevelCell, indexBoundaryLayerDepth

      real (kind=RKIND), dimension(:), pointer :: &
        latCell, lonCell, bottomDepth, buoyancyForcingOBL, surfaceFrictionVelocity, fCell, ssh, &
        boundaryLayerDepth
        
      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfCell, vertDiffTopOfCell, layerThickness, &
        bulkRichardsonNumber, zMid, zTop, density, displacedDensity, potentialDensity

      integer :: iCell, timeLevel
      integer, pointer :: nVertLevels, nCellsSolve

      real (kind=RKIND), pointer :: config_cvmix_background_viscosity, config_cvmix_background_diffusion

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing mixing-related fields
      ! note that the user can choose multiple options and the 
      !   mixing fields have to be merged together
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

      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

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

      !
      ! set pointers for fields related ocean forcing at surface
      !
      call mpas_pool_get_array(diagnosticsPool, 'surfaceFrictionVelocity', surfaceFrictionVelocity)
      call mpas_pool_get_array(diagnosticsPool, 'buoyancyForcingOBL', buoyancyForcingOBL)
      call mpas_pool_get_array(diagnosticsPool, 'boundaryLayerDepth', boundaryLayerDepth)
      call mpas_pool_get_array(diagnosticsPool, 'indexBoundaryLayerDepth', indexBoundaryLayerDepth)

      !
      ! set pointers for viscosity/diffusivity and intialize to zero
      !
      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfCell', vertViscTopOfCell)
      call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)
 
      vertViscTopOfCell = 0.0
      vertDiffTopOfCell = 0.0

      !
      ! start by setting visocity/diffusivity to there mininum background values
      !
      if (cvmixBackgroundOn) then
        vertViscTopOfCell(:,:) = vertViscTopOfCell(:,:) + config_cvmix_background_viscosity
        vertDiffTopOfCell(:,:) = vertDiffTopOfCell(:,:) + config_cvmix_background_diffusion
!       write(stdoutUnit,*) 'TDR: config_use_cvmix_background',config_use_cvmix_background,maxval(vertViscTopOfCell(:,:))
      endif

      !
      ! allocate selected cvmix variables and loop over columns
      !
      allocate(cvmix_variables % visc_iface(nVertLevels+1))
      allocate(cvmix_variables % diff_iface(nVertLevels+1,2))
      allocate(cvmix_variables % zw_iface(nVertLevels+1))

      do iCell = 1, nCellsSolve

        ! zero the cvmix viscosity/diffusivity
        cvmix_variables % visc_iface(:)=0.0
        cvmix_variables % diff_iface(:,:)=0.0

        ! fill the intent(in) convective adjustment
        cvmix_variables % nlev = maxLevelCell(iCell)
        cvmix_variables % ocn_depth = bottomDepth(iCell)
        cvmix_variables % dens => displacedDensity(:,iCell)
        cvmix_variables % dens_lwr => density(:,iCell)

        ! call convective mixing scheme
        ! (NOTE: cvmix_coeffs_conv might be incorrect)
        if (cvmixConvectionOn) then
          call cvmix_coeffs_conv(cvmix_variables)

          ! add convective mixing to vertical viscosity/diffusivity
          vertViscTopOfCell(:,iCell) = vertViscTopOfCell(:,iCell) + cvmix_variables % visc_iface(:)
          vertDiffTopOfCell(:,iCell) = vertDiffTopOfCell(:,iCell) + cvmix_variables % diff_iface(:,1)

        endif

        ! call kpp ocean mixed layer scheme
        if (cvmixKPPOn) then

!         write(stdoutUnit,*) 'TDR: config_use_cvmix_kpp enter',config_use_cvmix_kpp

          ! set cvmix viscosity/diffusity to current total values to be used for matching
          cvmix_variables % visc_iface(1:nVertLevels) = vertViscTopOfCell(1:nVertLevels,iCell)
          cvmix_variables % diff_iface(1:nVertLevels,1) = vertDiffTopofCell(1:nVertLevels,iCell)
          cvmix_variables % diff_iface(1:nVertLevels,2) = vertDiffTopofCell(1:nVertLevels,iCell)

          ! set matching visc/diff to zero at upper (1) and lower (nVertlevels+1) boundaries
          cvmix_variables % visc_iface(1) = 0.0
          cvmix_variables % visc_iface(nVertLevels+1) = 0.0
          cvmix_variables % diff_iface(1,:) = 0.0
          cvmix_variables % diff_iface(nVertLevels+1,:) = 0.0

!         write(stdoutUnit,*) 'TDR: cvmix_variables % visc_iface',maxval(cvmix_variables % visc_iface)
!         write(stdoutUnit,*) 'TDR: cvmix_variables % diff_iface',maxval(cvmix_variables % diff_iface)

          ! set integer and real
          cvmix_variables % surf_hgt = ssh(iCell)
          cvmix_variables % Coriolis = fCell(iCell)
          cvmix_variables % lat = latCell(iCell) * 180.0 / 3.14
          cvmix_variables % lon = lonCell(iCell) * 180.0 / 3.14
          cvmix_variables % surf_fric = surfaceFrictionVelocity(iCell)
          cvmix_variables % surf_buoy = buoyancyForcingOBL(iCell)

!         write(stdoutUnit,*) 'TDR: zTop',maxval(diagnostics % zTop % array(1:nVertLevels,iCell))
!         write(stdoutUnit,*) 'TDR: bottomDepth',maxval(mesh % bottomDepth % array(:))

          ! fill zw_iface with interface coordinates
          cvmix_variables % zw_iface(1:nVertLevels) = zTop(1:nVertLevels,iCell)
          cvmix_variables % zw_iface(nVertLevels+1) = -bottomDepth(iCell)

!         write(stdoutUnit,*) 'TDR: zw_iface',cvmix_variables % zw_iface(:)

          ! point remainder of cvmix_variables to MPAS arrays
          cvmix_variables % zt =>  zMid(1:nVertLevels,iCell)
          cvmix_variables % dzt => layerThickness(1:nVertLevels,iCell)
          cvmix_variables % Rib => bulkRichardsonNumber(1:nVertLevels,iCell)

          ! compute the boundary layer depth
!         write(stdoutUnit,*) 'TDR: calling cvmix_kpp_compute_OBL_depth'
          call cvmix_kpp_compute_OBL_depth( cvmix_variables )
!         write(stdoutUnit,*) 'TDR: return cvmix_kpp_compute_OBL_depth'

!         write(stdoutUnit,*) 'TDR: OBL_depth, kOBL_depth',cvmix_variables % OBL_depth,cvmix_variables % kOBL_depth

          ! intent out of OBL_depth is boundary layer depth measured in meters and vertical index
          boundaryLayerDepth(iCell) = cvmix_variables % OBL_depth
          indexBoundaryLayerDepth(iCell) = cvmix_variables % kOBL_depth

!         write(stdoutUnit,*) 'TDR: OBL_depth',diagnostics % boundaryLayerDepth % array(iCell)

          ! given OBL and vertical profile of visc/diff, compute boundary layer mixing
!         write(stdoutUnit,*) 'TDR: calling cvmix_coeffs_kpp'
          call cvmix_coeffs_kpp( cvmix_variables )
!         write(stdoutUnit,*) 'TDR: return cvmix_coeffs_kpp'

          ! add convective mixing to vertical viscosity/diffusivity
          vertViscTopOfCell(:,iCell) = vertViscTopOfCell(:,iCell) + cvmix_variables % visc_iface(:)
          vertDiffTopOfCell(:,iCell) = vertDiffTopOfCell(:,iCell) + cvmix_variables % diff_iface(:,1)

        endif

      enddo
      deallocate(cvmix_variables % visc_iface)
      deallocate(cvmix_variables % diff_iface)
      deallocate(cvmix_variables % zw_iface)

!     write(stdoutUnit,*) 'TDR: exiting ocn_vmix_coefs_cvmix_build'

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

      logical, pointer :: config_use_cvmix, config_use_cvmix_background, config_use_cvmix_convection
      logical, pointer :: config_use_cvmix_kpp
      real (kind=RKIND), pointer :: config_cvmix_background_viscosity, config_cvmix_background_diffusion
      real (kind=RKIND), pointer :: config_cvmix_prandtl_number, config_cvmix_convective_diffusion, config_cvmix_convective_viscosity


      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix', config_use_cvmix)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_background', config_use_cvmix_background)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_background_viscosity', config_cvmix_background_viscosity)
      call mpas_pool_get_config(ocnConfigs, 'config_cvmix_background_diffusion', config_cvmix_background_diffusion)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_convection', config_use_cvmix_convection)
      call mpas_pool_get_config(ocnConfigs, 'config_use_cvmix_kpp', config_use_cvmix_kpp)
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
      ! if using background values, put these into CVMix
      !
      if (config_use_cvmix_background) then
        call cvmix_init_bkgnd( &
               bkgnd_diff = config_cvmix_background_diffusion, &
               bkgnd_visc = config_cvmix_background_viscosity, &
               CVmix_bkgnd_params_user = cvmix_background_params)
!       write(stdoutUnit,*) 'cvmix_init_bkgnd',config_use_cvmix_background
      endif

      !
      ! if using CVMix convection, put these into CVMix
      !
      if (config_use_cvmix_convection) then
        call cvmix_init_conv( &
               convect_diff = config_cvmix_convective_diffusion,  &
               convect_visc = config_cvmix_convective_viscosity,  &
               CVmix_conv_params_user = cvmix_conv_params)
!       write(stdoutUnit,*) 'cvmix_init_conv',config_use_cvmix_convection
      endif

      !
      ! if using CVMix convection, put these into CVMix
      !
      if (config_use_cvmix_kpp) then
        call cvmix_init_kpp ( CVmix_kpp_params_user = cvmix_kpp_params )
!       write(stdoutUnit,*) 'cvmix_init_kpp',config_use_cvmix_kpp
      endif

!     write(stdoutUnit,*) 'ocn_vmix_cvmix_init complete'

   !--------------------------------------------------------------------

   end subroutine ocn_vmix_cvmix_init!}}}

!***********************************************************************

end module ocn_vmix_cvmix

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
