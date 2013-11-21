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
   use mpas_configure
   use mpas_timer
   use mpas_io_units

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
   type(cvmix_shear_params_type) :: cvmix_shear_params
   type(cvmix_kpp_params_type) :: cvmix_kpp_params
   type(cvmix_data_type) :: cvmix_variables


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

   subroutine ocn_vmix_coefs_cvmix_build(mesh, state, diagnostics, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (state_type), intent(inout) :: &
         state         !< Input/Output: state information

      type (diagnostics_type), intent(inout) :: &
         diagnostics   !< Input/Output: diagnostic information

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
        latCell, lonCell, bottomDepth, buoyancyForcingOBL, surfaceFrictionVelocity, fCell 
        
      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfCell, vertDiffTopOfCell, layerThickness, &
        bulkRichardsonNumber, zMid, zTop, density, displacedDensity, potentialDensity

      real (kind=RKIND), dimension(:,:,:), pointer :: vertNonLocalTrans
      integer :: index_nonLocalTemp, index_nonLocalSalt, index_nonLocalZonalVel, index_nonLocalMeridVel

      integer :: k, iCell, nVertLevels

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

      !
      ! only build up viscosity/diffusivity if CVMix is turned on
      !
      if(.not.config_use_cvmix) return

      !
      ! set parameters
      !
      nVertLevels = mesh % nVertLevels

      !
      ! set pointers for fields related to position on sphere
      !
      latCell => mesh % latCell % array(:)
      lonCell => mesh % lonCell % array(:)
      fCell => mesh % fCell % array(:)

      !
      ! set pointers for fields related to vertical mesh
      !
      maxLevelCell => mesh % maxLevelCell % array(:)
      bottomDepth => mesh % bottomDepth % array(:)
      layerThickness => state % layerThickness % array(:,:)
      zTop => diagnostics % zTop % array(:,:)
      zMid => diagnostics % zMid % array(:,:)

      !
      ! set pointers for fields related ocean state
      !
      density => diagnostics % density % array(:,:)
      displacedDensity => diagnostics % displacedDensity % array(:,:)
      potentialDensity => diagnostics % potentialDensity % array(:,:)
      bulkRichardsonNumber => diagnostics % bulkRichardsonNumber % array(:,:)

      !
      ! set pointers for fields related ocean forcing at surface
      !
      surfaceFrictionVelocity => diagnostics % surfaceFrictionVelocity % array(:)
      buoyancyForcingOBL => diagnostics % buoyancyForcingOBL % array(:)

      !
      ! set pointers for viscosity/diffusivity and intialize to zero
      !
      vertViscTopOfCell => diagnostics % vertViscTopOfCell % array(:,:)
      vertDiffTopOfCell => diagnostics % vertDiffTopOfCell % array(:,:)
 
      vertViscTopOfCell = 0.0
      vertDiffTopOfCell = 0.0

      !
      ! set pointer for nonlocal transport and intialize to zero
      !
      vertNonLocalTrans = diagnostics % vertNonLocalTrans % array(:,:,:)
      index_nonLocalTemp = diagnostics % index_vertNonLocalTransTemp
      index_nonLocalSalt = diagnostics % index_vertNonLocalTransSalt
      index_nonLocalZonalVel = diagnostics % index_vertNonLocalTransZonalVel
      index_nonLocalMeridVel = diagnostics % index_vertNonLocalTransMeridVel

      vertNonLocalTrans = 0.0

      !
      ! start by setting visocity/diffusivity to there mininum background values
      !
      if (config_use_cvmix_background) then
        vertViscTopOfCell(:,:) = vertViscTopOfCell(:,:) + config_cvmix_background_viscosity
        vertDiffTopOfCell(:,:) = vertDiffTopOfCell(:,:) + config_cvmix_background_diffusion
      endif

      !
      ! allocate selected cvmix variables and loop over columns
      !
      allocate(cvmix_variables % visc_iface(nVertLevels+1))
      allocate(cvmix_variables % diff_iface(nVertLevels+1,2))
      allocate(cvmix_variables % zw_iface(nVertLevels+1))
      allocate(cvmix_variables % kpp_transport_iface(nVertLevels+1,4))

      do iCell=1,mesh%nCellsSolve

        ! specify geometry/location
        cvmix_variables % surf_hgt = state % ssh % array(iCell)
        cvmix_variables % Coriolis = mesh % fCell % array(iCell)
        cvmix_variables % lat = mesh % latCell % array(iCell) * 180.0 / 3.14
        cvmix_variables % lon = mesh % lonCell % array(iCell) * 180.0 / 3.14
        cvmix_variables % zw_iface(1:nVertLevels) =  diagnostics % zTop % array(1:nVertLevels,iCell)
        cvmix_variables % zw_iface(nVertLevels+1) = -mesh % bottomDepth % array(iCell)
        cvmix_variables % zt =>  diagnostics % zMid % array(1:nVertLevels,iCell)
        cvmix_variables % dzt => state % layerThickness % array(1:nVertLevels,iCell)

        ! fill the intent(in) for shear mixing
        cvmix_variables%Ri_iface => diagnostics % RiTopOfCell % array(:,iCell)

        ! fill the intent(in) convective adjustment
        ! TDR: check density
        cvmix_variables % nlev = mesh % maxLevelCell % array(iCell)
        cvmix_variables % ocn_depth = mesh % bottomDepth % array(iCell)
        cvmix_variables % dens => diagnostics % displacedDensity % array(:,iCell)
        cvmix_variables % dens_lwr => diagnostics % density % array(:,iCell)
        cvmix_variables % buoy_iface => diagnostics % BruntVaisalaFreqTop % array(:,iCell)

        ! fill the intent(in) KPP
        cvmix_variables % surf_fric = diagnostics % surfaceFrictionVelocity % array(iCell)
        cvmix_variables % surf_buoy = diagnostics % buoyancyForcingOBL % array(iCell)
        cvmix_variables % Rib => diagnostics % BulkRichardsonNumber % array(1:nVertLevels,iCell)

        ! call shear-based mixing scheme
        if (config_use_cvmix_shear) then
          cvmix_variables % visc_iface(:)=0.0
          cvmix_variables % diff_iface(:,:)=0.0
          call cvmix_coeffs_shear( &
                 CVmix_vars = cvmix_variables, &
                 CVmix_shear_params = cvmix_shear_params, &
                 CVmix_bkgnd_params = cvmix_background_params )

          ! add convective mixing to vertical viscosity/diffusivity
          vertViscTopOfCell(:,iCell) = vertViscTopOfCell(:,iCell) + cvmix_variables % visc_iface(:)
          vertDiffTopOfCell(:,iCell) = vertDiffTopOfCell(:,iCell) + cvmix_variables % diff_iface(:,1)

        endif

        ! call convective mixing scheme
        if (config_use_cvmix_convection) then
          cvmix_variables % visc_iface(:)=0.0
          cvmix_variables % diff_iface(:,:)=0.0
          call cvmix_coeffs_conv( & 
                 CVmix_vars = cvmix_variables,  &
                 CVmix_conv_params_user = cvmix_conv_params)

          ! add convective mixing to vertical viscosity/diffusivity
          vertViscTopOfCell(:,iCell) = vertViscTopOfCell(:,iCell) + cvmix_variables % visc_iface(:)
          vertDiffTopOfCell(:,iCell) = vertDiffTopOfCell(:,iCell) + cvmix_variables % diff_iface(:,1)
        endif

        ! call kpp ocean mixed layer scheme
        if (config_use_cvmix_kpp) then
          cvmix_variables % visc_iface(:)=vertViscTopOfCell(:,iCell)
          cvmix_variables % diff_iface(:,1)=vertDiffTopOfCell(:,iCell)
          cvmix_variables % diff_iface(:,2)=vertDiffTopOfCell(:,iCell)
          cvmix_variables % kpp_transport_iface(:,:)=0.0

          ! compute the boundary layer depth
          call cvmix_kpp_compute_OBL_depth( &
                 CVmix_vars = cvmix_variables, &
                 CVmix_kpp_params_user = cvmix_kpp_params )

          ! intent out of OBL_depth is boundary layer depth measured in meters and vertical index
          diagnostics % boundaryLayerDepth % array(iCell) = cvmix_variables % OBL_depth
          diagnostics % indexBoundaryLayerDepth % array(iCell) = cvmix_variables % kOBL_depth

          ! given OBL and vertical profile of visc/diff, compute boundary layer mixing
          call cvmix_coeffs_kpp( &
                 CVmix_vars = cvmix_variables, &
                 CVmix_kpp_params_user = cvmix_kpp_params )

          ! overwrite viscosity/diffusivity from surface through OBL depth
          do k=1,cvmix_variables % kOBL_depth
            vertViscTopOfCell(:,iCell) = cvmix_variables % visc_iface(:)
            vertDiffTopOfCell(:,iCell) = cvmix_variables % diff_iface(:,1)
          enddo

          vertNonLocalTrans(index_nonLocalTemp,:,iCell) = cvmix_variables % kpp_transport_iface(:,1)
          vertNonLocalTrans(index_nonLocalSalt,:,iCell) = cvmix_variables % kpp_transport_iface(:,2)
          vertNonLocalTrans(index_nonLocalZonalVel,:,iCell) = cvmix_variables % kpp_transport_iface(:,3)
          vertNonLocalTrans(index_nonLocalMeridVel,:,iCell) = cvmix_variables % kpp_transport_iface(:,4)

        endif

      enddo
      deallocate(cvmix_variables % visc_iface)
      deallocate(cvmix_variables % diff_iface)
      deallocate(cvmix_variables % zw_iface)
      deallocate(cvmix_variables % kpp_transport_iface)

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

      integer :: nVertLevels
      type (block_type), pointer :: block
      type (mesh_type), pointer :: mesh


      !
      ! assume no errors during initialization and set to 1 when error is encountered
      !
      err=0

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
      block => domain % blocklist
      do while (associated(block))
         mesh => block % mesh
         nVertLevels = mesh % nVertLevels
         block => block % next
      end do

      !
      ! put global parameters into CVMix
      !
      call cvmix_put(cvmix_global_params,  'max_nlev', nVertLevels)
      call cvmix_put(cvmix_global_params,  'prandtl',  config_cvmix_prandtl_number)

      !
      ! if using background mixing, put parameters into cvmix_background_params
      !
      if (config_use_cvmix_background) then
        call cvmix_init_bkgnd( &
               bkgnd_diff = config_cvmix_background_diffusion, &
               bkgnd_visc = config_cvmix_background_viscosity, &
               CVmix_bkgnd_params_user = cvmix_background_params)
      endif

      !
      ! if using shear-based mixing, put parameters into cvmix_shear_params
      !
      if (config_use_cvmix_shear) then
        call cvmix_init_shear( &
               CVmix_shear_params = cvmix_shear_params, &
               mix_scheme = config_cvmix_shear_mixing_scheme, &
               PP_nu_zero = config_cvmix_shear_PP_nu_zero, &
               PP_alpha = config_cvmix_shear_PP_alpha, &
               PP_exp = config_cvmix_shear_PP_exp, &
               KPP_nu_zero = config_cvmix_shear_KPP_nu_zero, &
               KPP_Ri_zero = config_cvmix_shear_KPP_Ri_zero, &
               KPP_exp = config_cvmix_shear_KPP_exp)
        write(6,*) 'cvmix_init_shear',config_use_cvmix_shear
      endif

      !
      ! if using CVMix convection, put parameters into cvmix_conv_params
      !
      if (config_use_cvmix_convection) then
        call cvmix_init_conv( &
               convect_diff = config_cvmix_convective_diffusion,  &
               convect_visc = config_cvmix_convective_viscosity,  &
               lBruntVaisala = config_cvmix_convective_basedOnBVF, &
               BVsqr_convect = config_cvmix_convective_triggerBVF, &
               CVmix_conv_params_user = cvmix_conv_params)
      endif

      !
      ! if using the KPP boundary layer scheme, put parameters into cvmix_kpp_params
      !
      if (config_use_cvmix_kpp) then
        call cvmix_init_kpp ( &
               ri_crit = config_cvmix_kpp_criticalBulkRichardsonNumber, &
               interp_type = config_cvmix_kpp_interpolationOMLType, &
               lEkman = config_cvmix_kpp_EkmanOBL, &
               lMonOb = config_cvmix_kpp_MonObOBL, &
               MatchTechnique = config_cvmix_kpp_matching, &
               CVmix_kpp_params_user = cvmix_kpp_params )
      endif


   !--------------------------------------------------------------------

   end subroutine ocn_vmix_cvmix_init!}}}

!***********************************************************************

end module ocn_vmix_cvmix

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker

