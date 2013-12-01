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
   type(cvmix_shear_params_type)  :: cvmix_shear_params
   type(cvmix_tidal_params_type)  :: cvmix_tidal_params
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
        latCell, lonCell, bottomDepth, surfaceBuoyancyForcing, surfaceFrictionVelocity, fCell 
        
      real (kind=RKIND), dimension(:,:), pointer :: &
        vertViscTopOfCell, vertDiffTopOfCell, layerThickness, &
        zMid, zTop, density, displacedDensity, potentialDensity, &
        bulkRichardsonNumber

      real (kind=RKIND), dimension(:,:,:), pointer :: vertNonLocalFlux
      integer :: index_nonLocalTemp, index_nonLocalSalt, index_nonLocalZonalVel, index_nonLocalMeridVel

      integer :: k, iCell, nVertLevels, jCell, iNeighbor
      real (kind=RKIND) :: r, layerSum, iter
      real (kind=RKIND), dimension(:), allocatable :: sigma, turbulentScalarVelocityScale

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
      ! set pointers for fields related forcing at ocean surface
      !
      surfaceFrictionVelocity => diagnostics % surfaceFrictionVelocity % array(:)
      surfaceBuoyancyForcing  => diagnostics % surfaceBuoyancyForcing % array(:)

      !
      ! set pointers for viscosity/diffusivity and intialize to zero
      !
      vertViscTopOfCell => diagnostics % vertViscTopOfCell % array(:,:)
      vertDiffTopOfCell => diagnostics % vertDiffTopOfCell % array(:,:)
 
      vertViscTopOfCell = 0.0
      vertDiffTopOfCell = 0.0

      !
      ! set pointers for nonlocal Flux and intialize to zero
      !
      vertNonLocalFlux => diagnostics % vertNonLocalFlux % array
      index_nonLocalTemp = diagnostics % index_vertNonLocalFluxTemp
      index_nonLocalSalt = diagnostics % index_vertNonLocalFluxSalt
      index_nonLocalZonalVel = diagnostics % index_vertNonLocalFluxZonalVel
      index_nonLocalMeridVel = diagnostics % index_vertNonLocalFluxMeridVel

      vertNonLocalFlux = 0.0

      !
      ! start by adding the mininum background values to the visocity/diffusivity arrays
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
      allocate(cvmix_variables % Rib(nVertLevels))

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
        cvmix_variables % nlev = mesh % maxLevelCell % array(iCell)
        cvmix_variables % ocn_depth = mesh % bottomDepth % array(iCell)
        cvmix_variables % dens => diagnostics % density % array(:,iCell)
        cvmix_variables % dens_lwr => diagnostics % displacedDensity % array(:,iCell)
        cvmix_variables % buoy_iface => diagnostics % BruntVaisalaFreqTop % array(:,iCell)

        ! fill the intent(in) KPP
        cvmix_variables % surf_fric = diagnostics % surfaceFrictionVelocity % array(iCell)
        cvmix_variables % surf_buoy = diagnostics % surfaceBuoyancyForcing % array(iCell)

        ! call shear-based mixing scheme
        if (config_use_cvmix_shear) then
          cvmix_variables % visc_iface(:)=0.0
          cvmix_variables % diff_iface(:,:)=0.0
          call cvmix_coeffs_shear( &
                 CVmix_vars = cvmix_variables, &
                 CVmix_shear_params = cvmix_shear_params, &
                 CVmix_bkgnd_params = cvmix_background_params )

          ! add shear mixing to vertical viscosity/diffusivity
          vertViscTopOfCell(:,iCell) = vertViscTopOfCell(:,iCell) + cvmix_variables % visc_iface(:)
          vertDiffTopOfCell(:,iCell) = vertDiffTopOfCell(:,iCell) + cvmix_variables % diff_iface(:,1)

          ! The following pertains to the 'PP' option of shear-generated mixing
          ! NOTE: at present, cvmix_coeffs_shear add in the background viscosity/diffusivity
          ! NOTE: the following two source code lines remove this double counting of the background
          ! NOTE: the following two source code lines should be deleted with cvmix_coeffs_shear is "fixed"
          if(config_cvmix_shear_mixing_scheme.eq.'PP') then
            vertViscTopOfCell(:,iCell) = vertViscTopOfCell(:,iCell) - config_cvmix_background_viscosity
            vertDiffTopOfCell(:,iCell) = vertDiffTopOfCell(:,iCell) - config_cvmix_background_diffusion
          endif

        endif ! if (config_use_cvmix_shear)

        ! TDR: put tidal here
        ! TDR: put ddiff here

        ! call kpp ocean mixed layer scheme
        if (config_use_cvmix_kpp) then
          cvmix_variables % visc_iface(:)=vertViscTopOfCell(:,iCell)
          cvmix_variables % diff_iface(:,1)=vertDiffTopOfCell(:,iCell)
          cvmix_variables % diff_iface(:,2)=vertDiffTopOfCell(:,iCell)
          cvmix_variables % kpp_transport_iface(:,:)=0.0

          ! copy boundary layer depth from previous time step into cvmix variables
          cvmix_variables % OBL_depth = diagnostics % boundaryLayerDepth % array(iCell)
          cvmix_variables % kOBL_depth = diagnostics % indexBoundaryLayerDepth % array(iCell)

          ! given (old time level) OBL and vertical profile of visc/diff, compute boundary layer mixing
          call cvmix_coeffs_kpp( CVmix_vars = cvmix_variables )

          ! if using KPP with "MatchBoth"   matching, then the output from KPP is the full viscosity/diffusivity
          ! if using KPP with "SimpleShape" matching, then the output from KPP needs to be added to current viscosity/diffusivity
          if(config_cvmix_kpp_matching.eq."MatchBoth") then
            vertViscTopOfCell(:,iCell) = cvmix_variables % visc_iface(:)
            vertDiffTopOfCell(:,iCell) = cvmix_variables % diff_iface(:,1)
          elseif(config_cvmix_kpp_matching.eq."SimpleShapes") then
            vertViscTopOfCell(:,iCell) = vertViscTopOfCell(:,iCell) + cvmix_variables % visc_iface(:)
            vertDiffTopOfCell(:,iCell) = vertDiffTopOfCell(:,iCell) + cvmix_variables % diff_iface(:,1)
          else
            stop
          endif

          ! store non-local flux terms
          ! these flux terms must be multiplied by the surfaceTracerFlux field
          ! the tracer tendency is then the vertical divergence of vertNonLocalFlux*surfaceTracerFlux
          ! both of these operations are done in ocn_tracer_nonlocalflux_tend routine
          ! NOTE: at present CVMix KPP stores four flux profiles (temperature, tracers, u velocity and v velocity)
          ! NOTE: at some point we might end up with only one flux profile (and thereby reduce storage of 3D arrays)
          vertNonLocalFlux(index_nonLocalTemp,:,iCell) = cvmix_variables % kpp_transport_iface(:,1)
          vertNonLocalFlux(index_nonLocalSalt,:,iCell) = cvmix_variables % kpp_transport_iface(:,2)
          vertNonLocalFlux(index_nonLocalZonalVel,:,iCell) = cvmix_variables % kpp_transport_iface(:,3)
          vertNonLocalFlux(index_nonLocalMeridVel,:,iCell) = cvmix_variables % kpp_transport_iface(:,4)

        endif ! if (config_use_cvmix_kpp)

        ! call convective mixing scheme
        if (config_use_cvmix_convection) then
          cvmix_variables % visc_iface(:)=0.0
          cvmix_variables % diff_iface(:,:)=0.0
          call cvmix_coeffs_conv( CVmix_vars = cvmix_variables )

          ! add convective mixing to vertical viscosity/diffusivity
          ! if using KPP, then do not apply convective mixing within the ocean boundary layer
          if(config_use_cvmix_kpp) then
            do k=int(diagnostics % indexBoundaryLayerDepth % array(iCell))+1,mesh % maxLevelCell % array(iCell)
              vertViscTopOfCell(k,iCell) = vertViscTopOfCell(k,iCell) + cvmix_variables % visc_iface(k)
              vertDiffTopOfCell(k,iCell) = vertDiffTopOfCell(k,iCell) + cvmix_variables % diff_iface(k,1)
            enddo
          else
            vertViscTopOfCell(:,iCell) = vertViscTopOfCell(:,iCell) + cvmix_variables % visc_iface(:)
            vertDiffTopOfCell(:,iCell) = vertDiffTopOfCell(:,iCell) + cvmix_variables % diff_iface(:,1)
          endif
        endif  ! if (config_use_cvmix_convection)

      enddo  ! do iCell=1,mesh%nCellsSolve

      ! 
      ! compute and smooth CVMix/KPP computed ocean boundary layer depth
      !
      if (config_use_cvmix_kpp) then

        allocate(sigma(nVertLevels))
        allocate(turbulentScalarVelocityScale(nVertLevels))

        do iCell=1,mesh%nCells
  
         ! the computation of the boundary layer depth is iterative
         do iter=1,config_cvmix_kpp_niteratorOBL_solver

           ! compute sigma at cell center
           sigma(1) = 0.0
           layerSum = 0.0
           do k=2,maxLevelCell(iCell)
             layerSum = layerSum + 0.5*(layerThickness(k-1,iCell)+layerThickness(k,iCell))
             sigma(k) = layerSum / diagnostics % boundaryLayerDepth % array(iCell)
           enddo
           do k=maxLevelCell(iCell),nVertLevels
             sigma(k) = sigma(maxLevelCell(iCell))
           enddo

           ! compute the turbulent scales in order to compute the bulk Richardson number
           call cvmix_kpp_compute_turbulent_scales( &
             sigma(:), &
             diagnostics % boundaryLayerDepth % array(iCell), &
             diagnostics % surfaceBuoyancyForcing % array(iCell), &
             diagnostics % surfaceFrictionVelocity % array(iCell),  &
             w_s = turbulentScalarVelocityScale(:) )

           cvmix_variables % Rib  = cvmix_kpp_compute_bulk_Richardson( &
             diagnostics % zMid % array(1:nVertLevels,iCell), & 
             diagnostics % bulkRichardsonNumberBuoy % array(1:nVertLevels,iCell), &
             diagnostics % bulkRichardsonNumberShear % array(1:nVertLevels,iCell), &
             ws_cntr = turbulentScalarVelocityScale(:), &
             Nsqr_iface = diagnostics % BruntVaisalaFreqTop % array(:,iCell) )

           ! compute the boundary layer depth based on model state at current time step
           ! boundaryLayerDepth is then smoothed and a halo update is applied elsewhere
           call cvmix_kpp_compute_OBL_depth( CVmix_vars = cvmix_variables)

           ! update boundary layer depth for sigma computation
           diagnostics % boundaryLayerDepth % array(iCell) = cvmix_variables % OBL_depth

         enddo

         ! intent out of OBL_depth is boundary layer depth measured in meters and vertical index
         diagnostics % boundaryLayerDepth % array(iCell) = cvmix_variables % OBL_depth
         diagnostics % indexBoundaryLayerDepth % array(iCell) = cvmix_variables % kOBL_depth
         diagnostics % bulkRichardsonNumber % array(:,iCell) = cvmix_variables % Rib(:)
  
        end do  ! do iCell=1,mesh%nCells
  
        do iCell=1,mesh%nCells
  
          ! loop over all cell neighbors and add up boundary layer depth
          r = 0.0
          do iNeighbor=1,mesh%nEdgesOnCell%array(iCell)
            jCell = mesh%CellsOnCell%array(iNeighbor,iCell)
            r = r + diagnostics % boundaryLayerDepth % array(jCell)
          end do

          ! replace boundary layer depth with average of itself and neigbhbors
          diagnostics % boundaryLayerDepth % array(iCell) = &
             0.5*(diagnostics % boundaryLayerDepth % array(iCell)+r/mesh%nEdgesOnCell%array(iCell))
  
        end do  ! do iCell=1,mesh%nCells

        deallocate(sigma)
        deallocate(turbulentScalarVelocityScale)

      end if ! if (config_use_cvmix_kpp)

      ! dellocate cmvix variables
      deallocate(cvmix_variables % visc_iface)
      deallocate(cvmix_variables % diff_iface)
      deallocate(cvmix_variables % zw_iface)
      deallocate(cvmix_variables % kpp_transport_iface)
      deallocate(cvmix_variables % Rib)

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
      ! initialize background mixing
      !
      if (config_use_cvmix_background) then
        call cvmix_init_bkgnd( &
               bkgnd_diff = config_cvmix_background_diffusion, &
               bkgnd_visc = config_cvmix_background_viscosity, &
               CVmix_bkgnd_params_user = cvmix_background_params)
      endif

      !
      ! initialize shear-based mixing
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

! vim: foldmethod=marker


