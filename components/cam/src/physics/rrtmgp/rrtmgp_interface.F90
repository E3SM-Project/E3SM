! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2017,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
!
! This module provides an interface to RRTMGP for a common use case --
!   users want to start from gas concentrations, pressures, and temperatures,
!   and compute clear-sky (aerosol plus gases) and all-sky fluxes.
! The routines here have the same names as those in mo_rrtmgp_[ls]w; normally users
!   will use either this module or the underling modules, but not both
module rrtmgp_interface
  use mo_rte_kind,   only: wp
  use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props, only: ty_optical_props_1scl, ty_optical_props_2str
  use mo_source_functions, only: ty_source_func_lw
  use mo_fluxes_byband, only: ty_fluxes_byband
  use mo_rte_lw, only: base_rte_lw => rte_lw
  use mo_rte_sw, only: base_rte_sw => rte_sw
  implicit none
  private
  public :: &
     rrtmgp_initialize, rrtmgp_finalize, &
     rte_lw, rte_sw, &
     get_nband, get_ngpt, get_band_lims_wavenumber, get_band_midpoints, &
     get_temp_min, get_temp_max, get_gpoint_bands

  public :: k_dist_lw, k_dist_sw
  type(ty_gas_optics_rrtmgp) :: k_dist_lw, k_dist_sw
  integer, parameter :: ERROR_MSG_LEN = 128

contains
   !----------------------------------------------------------------------------
   function rrtmgp_initialize(coefficients_file_sw, coefficients_file_lw, active_gases) result(error_msg)
      use rrtmgp_coefficients, only: load_and_init
      use mo_gas_concentrations, only: ty_gas_concs
      character(len=*), intent(in) :: coefficients_file_sw, coefficients_file_lw
      character(len=*), dimension(:), intent(in) :: active_gases
      character(len=ERROR_MSG_LEN) :: error_msg
      ! ty_gas_concs object that would normally hold volume mixing ratios for
      ! radiatively-important gases. Here, this is just used to provide the names
      ! of gases that are available in the model (needed by the kdist
      ! initialization routines that are called within the load_coefficients
      ! methods).
      type(ty_gas_concs) :: available_gases

      ! Need to initialize available_gases here! The only field of the
      ! available_gases type that is used in the kdist initialize is
      ! available_gases%gas_name, which gives the name of each gas that would be
      ! present in the ty_gas_concs object. So, we can just set this here, rather
      ! than trying to fully populate the ty_gas_concs object here, which would be
      ! impossible from this initialization routine because I do not think the
      ! rad_cnst objects are setup yet.
      ! TODO: This needs to be fixed to ONLY read in the data if masterproc, and then broadcast
      ! the other tasks!
      error_msg = set_available_gases(active_gases, available_gases)
      if (error_msg /= '') return

      ! Read coefficients from file and populate k-dist objects
      call load_and_init(k_dist_sw, coefficients_file_sw, available_gases)
      call load_and_init(k_dist_lw, coefficients_file_lw, available_gases)
   end function rrtmgp_initialize
   !----------------------------------------------------------------------------
   function set_available_gases(gases, gas_concentrations) result(error_msg)

      use mo_gas_concentrations, only: ty_gas_concs
      use mo_rrtmgp_util_string, only: lower_case

      type(ty_gas_concs), intent(inout) :: gas_concentrations
      character(len=*), intent(in) :: gases(:)
      character(len=32), dimension(size(gases)) :: gases_lowercase
      character(len=ERROR_MSG_LEN) :: error_msg
      integer :: igas

      error_msg = ''
      ! Initialize with lowercase gas names; we should work in lowercase
      ! whenever possible because we cannot trust string comparisons in RRTMGP
      ! to be case insensitive
      do igas = 1,size(gases)
         gases_lowercase(igas) = trim(lower_case(gases(igas)))
      end do
      error_msg = gas_concentrations%init(gases_lowercase)

   end function set_available_gases
   !----------------------------------------------------------------------------
   function rrtmgp_finalize() result(error_msg)
      character(len=ERROR_MSG_LEN) :: error_msg
      error_msg = ''
   end function rrtmgp_finalize
   !----------------------------------------------------------------------------
   function get_nband(band) result(nband)
      character(len=*), intent(in) :: band
      integer :: nband
      if (trim(band) == 'sw') then
         nband = k_dist_sw%get_nband()
      else if (trim(band) == 'lw') then
         nband = k_dist_lw%get_nband()
      else
         nband = -1
      end if
   end function get_nband
   !----------------------------------------------------------------------------
   function get_ngpt(band) result(ngpt)
      character(len=*), intent(in) :: band
      integer :: ngpt
      if (trim(band) == 'sw') then
         ngpt = k_dist_sw%get_ngpt()
      else if (trim(band) == 'lw') then
         ngpt = k_dist_lw%get_ngpt()
      else
         ngpt = -1
      end if
   end function get_ngpt
   !----------------------------------------------------------------------------
   function get_band_lims_wavenumber(band, nband) result(band_limits)
      character(len=*), intent(in) :: band
      integer, intent(in) :: nband
      real(wp), dimension(2,nband) :: band_limits
      if (trim(band) == 'sw') then
         band_limits = k_dist_sw%get_band_lims_wavenumber()
      else if (trim(band) == 'lw') then
         band_limits = k_dist_lw%get_band_lims_wavenumber()
      else
         band_limits = -1
      end if
   end function get_band_lims_wavenumber
   !----------------------------------------------------------------------------
   real(wp) function get_temp_min()
      get_temp_min = min(k_dist_sw%get_temp_min(), k_dist_lw%get_temp_min())
   end function get_temp_min
   !----------------------------------------------------------------------------
   real(wp) function get_temp_max()
      get_temp_max = max(k_dist_sw%get_temp_max(), k_dist_lw%get_temp_max())
   end function get_temp_max
   !----------------------------------------------------------------------------
   function get_gpoint_bands(band, ngpt)
      character(len=*), intent(in) :: band
      integer, intent(in) :: ngpt
      integer, dimension(ngpt) :: get_gpoint_bands
      if (trim(band) == 'sw') then
         get_gpoint_bands = k_dist_sw%get_gpoint_bands()
      else if (trim(band) == 'lw') then
         get_gpoint_bands = k_dist_lw%get_gpoint_bands()
      else
         get_gpoint_bands = -1
      end if
   end function
   !----------------------------------------------------------------------------
  ! --------------------------------------------------
  !
  ! Interfaces using clear (gas + aerosol) and all-sky categories, starting from
  !   pressures, temperatures, and gas amounts for the gas contribution
  !
  ! --------------------------------------------------
  function rte_lw( &
        gas_names, gas_vmr, p_lay, t_lay, p_lev,    &
        t_sfc, sfc_emis, cld_tau, aer_tau,          &
        allsky_flux_up, allsky_flux_dn, allsky_flux_net, &
        allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, &
        clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net, &
        clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net, &
        col_dry, t_lev, inc_flux, n_gauss_angles) result(error_msg)
    character(len=*), dimension(:), intent(in) :: gas_names
    real(wp), dimension(:,:,:),    intent(in   ) :: gas_vmr
    real(wp), dimension(:,:),    intent(in   ) :: p_lay, t_lay !< pressure [Pa], temperature [K] at layer centers (ncol,nlay)
    real(wp), dimension(:,:),    intent(in   ) :: p_lev        !< pressure at levels/interfaces [Pa] (ncol,nlay+1)
    real(wp), dimension(:),      intent(in   ) :: t_sfc     !< surface temperature           [K]  (ncol)
    real(wp), dimension(:,:),    intent(in   ) :: sfc_emis  !< emissivity at surface         []   (nband, ncol)
    real(wp), dimension(:,:,:),  intent(in   ) :: cld_tau   ! cloud absorption optical depth (ncol,nlay,ngpt)
    real(wp), dimension(:,:,:),  intent(in   ) :: aer_tau   ! cloud absorption optical depth (ncol,nlay,nband)
    real(wp), dimension(:,:), intent(inout), target :: &
       allsky_flux_up, allsky_flux_dn, allsky_flux_net, &
       clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net
    real(wp), dimension(:,:,:), intent(inout), target :: &
       allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, &
       clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net

    ! Optional inputs
    real(wp), dimension(:,:), &
              optional,       intent(in ) :: col_dry !< Molecular number density (ncol, nlay)
    real(wp), dimension(:,:), target, &
              optional,       intent(in ) :: t_lev     !< temperature at levels [K] (ncol, nlay+1)
    real(wp), dimension(:,:), target, &
              optional,       intent(in ) :: inc_flux   !< incident flux at domain top [W/m2] (ncol, ngpts)
    integer,  optional,       intent(in ) :: n_gauss_angles ! Number of angles used in Gaussian quadrature (no-scattering solution)
    character(len=128)                    :: error_msg

    ! Local variables
    type(ty_gas_concs) :: gas_concs    !< derived type encapsulating gas concentrations
    type(ty_optical_props_1scl) :: optical_props, cld_props, aer_props
    type(ty_source_func_lw) :: sources
    type(ty_fluxes_byband) :: allsky_fluxes, clrsky_fluxes
    integer :: ncol, nlay, ngpt, nband, nstr, igas
    logical :: top_at_1

    error_msg = ""

    ! Problem sizes
    ncol  = size(p_lay, 1)
    nlay  = size(p_lay, 2)
    ngpt  = k_dist_lw%get_ngpt()
    nband = k_dist_lw%get_nband()

    ! Flag for TOA->SFC or SFC->TOA
    top_at_1 = p_lay(1, 1) < p_lay(1, nlay)

    ! Associate output pointers
    allsky_fluxes%flux_up => allsky_flux_up
    allsky_fluxes%flux_dn => allsky_flux_dn
    allsky_fluxes%flux_net => allsky_flux_net
    clrsky_fluxes%flux_up => clrsky_flux_up
    clrsky_fluxes%flux_dn => clrsky_flux_dn
    clrsky_fluxes%flux_net => clrsky_flux_net
    allsky_fluxes%bnd_flux_up => allsky_bnd_flux_up
    allsky_fluxes%bnd_flux_dn => allsky_bnd_flux_dn
    allsky_fluxes%bnd_flux_net => allsky_bnd_flux_net
    clrsky_fluxes%bnd_flux_up => clrsky_bnd_flux_up
    clrsky_fluxes%bnd_flux_dn => clrsky_bnd_flux_dn
    clrsky_fluxes%bnd_flux_net => clrsky_bnd_flux_net

    ! Setup gas concentrations
    error_msg = gas_concs%init(gas_names)
    if (error_msg /= '') return
    do igas = 1,size(gas_names)
       error_msg = gas_concs%set_vmr(gas_names(igas), gas_vmr(igas,:,:))
    end do
    if (error_msg /= '') return

    ! Optical properties arrays
    error_msg = optical_props%alloc_1scl(ncol, nlay, k_dist_lw)
    if (error_msg /= '') return
    error_msg = cld_props%alloc_1scl(ncol, nlay, k_dist_lw)
    if (error_msg /= '') return
    cld_props%tau = cld_tau
    error_msg = aer_props%alloc_1scl(ncol, nlay, k_dist_lw%get_band_lims_wavenumber())
    if (error_msg /= '') return
    aer_props%tau = aer_tau

    ! Source function
    error_msg = sources%init(k_dist_lw)
    error_msg = sources%alloc(ncol, nlay)
    if (error_msg /= '') return

    ! Gas optical depth -- pressure need to be expressed as Pa
    error_msg = k_dist_lw%gas_optics(p_lay, p_lev, t_lay, t_sfc, gas_concs, &
                                  optical_props, sources,                &
                                  col_dry, t_lev)
    if (error_msg /= '') return

    ! Clear sky fluxes (gases + aerosols)
    error_msg = aer_props%increment(optical_props)
    if(error_msg /= '') return
    error_msg = base_rte_lw(optical_props, top_at_1, sources, &
                            sfc_emis, clrsky_fluxes,          &
                            inc_flux, n_gauss_angles)
    if(error_msg /= '') return

    ! All-sky fluxes (clear skies + clouds)
    error_msg = cld_props%increment(optical_props)
    if(error_msg /= '') return

    error_msg = base_rte_lw(optical_props, top_at_1, sources, &
                            sfc_emis, allsky_fluxes,          &
                            inc_flux, n_gauss_angles)

    call sources%finalize()
  end function rte_lw
  ! --------------------------------------------------
  function rte_sw( &
        gas_names, gas_vmr, p_lay, t_lay, p_lev, &
        mu0, sfc_alb_dir, sfc_alb_dif,  &
        cld_tau, cld_ssa, cld_asm,      &
        aer_tau, aer_ssa, aer_asm,      &
        allsky_flux_up, allsky_flux_dn, allsky_flux_net, &
        allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, allsky_bnd_flux_dn_dir, &
        clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net, &
        clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net, clrsky_bnd_flux_dn_dir, &
        col_dry, inc_flux, tsi_scaling  ) result(error_msg)
      character(len=*), dimension(:), intent(in   ) :: gas_names
      real(wp), dimension(:,:,:),  intent(in   ) :: gas_vmr
      real(wp), dimension(:,:),    intent(in   ) :: p_lay, t_lay !< pressure [Pa], temperature [K] at layer centers (ncol,nlay)
      real(wp), dimension(:,:),    intent(in   ) :: p_lev        !< pressure at levels/interfaces [Pa] (ncol,nlay+1)
      real(wp), dimension(:  ),    intent(in   ) :: mu0          !< cosine of solar zenith angle
      !  surface albedo for direct and diffuse radiation (band, col)
      real(wp), dimension(:,:),    intent(in   ) :: sfc_alb_dir, sfc_alb_dif
      real(wp), dimension(:,:,:),  intent(in   ) :: cld_tau, cld_ssa, cld_asm
      real(wp), dimension(:,:,:),  intent(in   ) :: aer_tau, aer_ssa, aer_asm
      real(wp), dimension(:,:),    intent(inout), target :: &
         allsky_flux_up, allsky_flux_dn, allsky_flux_net, &
         clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net
      real(wp), dimension(:,:,:), intent(inout), target :: &
         allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, allsky_bnd_flux_dn_dir, &
         clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net, clrsky_bnd_flux_dn_dir

      ! Optional inputs
      real(wp), dimension(:,:), &
                optional,       intent(in ) :: col_dry, &  !< Molecular number density (ncol, nlay)
                                               inc_flux    !< incident flux at domain top [W/m2] (ncol, ngpts)
      real(wp), optional,       intent(in ) :: tsi_scaling !< Optional scaling for total solar irradiance

      character(len=128)                    :: error_msg
      ! Local variables
      type(ty_optical_props_2str) :: optical_props, cloud_props, aer_props
      type(ty_gas_concs) :: gas_concs    !< derived type encapsulating gas concentrations
      type(ty_fluxes_byband) :: allsky_fluxes, clrsky_fluxes
      real(wp), dimension(:,:), allocatable :: toa_flux
      integer :: ncol, nlay, ngpt, nband, nstr, igas
      logical :: top_at_1

      error_msg = ""

      ! Problem sizes
      ncol  = size(p_lay, 1)
      nlay  = size(p_lay, 2)
      ngpt  = k_dist_sw%get_ngpt()
      nband = k_dist_sw%get_nband()
  
      ! Flag for TOA->SFC or SFC->TOA
      top_at_1 = p_lay(1, 1) < p_lay(1, nlay)

      ! Associate output pointers
      allsky_fluxes%flux_up => allsky_flux_up
      allsky_fluxes%flux_dn => allsky_flux_dn
      allsky_fluxes%flux_net => allsky_flux_net
      clrsky_fluxes%flux_up => clrsky_flux_up
      clrsky_fluxes%flux_dn => clrsky_flux_dn
      clrsky_fluxes%flux_net => clrsky_flux_net
      allsky_fluxes%bnd_flux_up => allsky_bnd_flux_up
      allsky_fluxes%bnd_flux_dn => allsky_bnd_flux_dn
      allsky_fluxes%bnd_flux_net => allsky_bnd_flux_net
      allsky_fluxes%bnd_flux_dn_dir => allsky_bnd_flux_dn_dir
      clrsky_fluxes%bnd_flux_up => clrsky_bnd_flux_up
      clrsky_fluxes%bnd_flux_dn => clrsky_bnd_flux_dn
      clrsky_fluxes%bnd_flux_net => clrsky_bnd_flux_net
      clrsky_fluxes%bnd_flux_dn_dir => clrsky_bnd_flux_dn_dir

      ! Setup gas concentrations
      error_msg = gas_concs%init(gas_names)
      if (error_msg /= '') return
      do igas = 1,size(gas_names)
         error_msg = gas_concs%set_vmr(gas_names(igas), gas_vmr(igas,:,:))
      end do
      if (error_msg /= '') return

      ! Optical properties arrays
      error_msg = optical_props%alloc_2str(ncol, nlay, k_dist_sw)
      error_msg = cloud_props%alloc_2str(ncol, nlay, k_dist_sw)
      error_msg = aer_props%alloc_2str(ncol, nlay, k_dist_sw%get_band_lims_wavenumber())
      if (error_msg /= '') return
  
      ! Populate optical properties
      cloud_props%tau = cld_tau
      cloud_props%ssa = cld_ssa
      cloud_props%g   = cld_asm
      aer_props%tau = aer_tau
      aer_props%ssa = aer_ssa
      aer_props%g   = aer_asm

      ! Delta scale
      error_msg = cloud_props%delta_scale()
      if (error_msg /= '') return
  
      ! TOA flux
      allocate(toa_flux(ncol, ngpt))
  
      ! Gas optical depth -- pressure need to be expressed as Pa
      error_msg = k_dist_sw%gas_optics(p_lay, p_lev, t_lay, gas_concs,  &
                                    optical_props, toa_flux,                          &
                                    col_dry)
      if (error_msg /= '') return
  
      ! If users have supplied an incident flux, use that
      if(present(inc_flux))    toa_flux(:,:) = inc_flux(:,:)
      if(present(tsi_scaling)) toa_flux(:,:) = toa_flux(:,:) * tsi_scaling
  
      ! Clear sky is gases + aerosols (if they're supplied)
      error_msg = aer_props%increment(optical_props)
      if(error_msg /= '') return
      error_msg = base_rte_sw(optical_props, top_at_1, &
                                 mu0, toa_flux,                   &
                                 sfc_alb_dir, sfc_alb_dif,        &
                                 clrsky_fluxes)
      if(error_msg /= '') return
  
      ! All-sky fluxes = clear skies + clouds
      error_msg = cloud_props%increment(optical_props)
      if(error_msg /= '') return
      error_msg = base_rte_sw(optical_props, top_at_1,  &
                                 mu0, toa_flux,                   &
                                 sfc_alb_dir, sfc_alb_dif,        &
                                 allsky_fluxes)
  
   end function rte_sw

   !----------------------------------------------------------------------------
   ! Function to calculate band midpoints from kdist band limits
   function get_band_midpoints(band, band_midpoints) result(error_msg)
      character(len=*), intent(in) :: band
      real(wp), intent(inout) :: band_midpoints(:)
      character(len=ERROR_MSG_LEN) :: error_msg
      real(wp) :: band_limits(2,size(band_midpoints))
      character(len=128), parameter :: subname = 'get_band_midpoints'
      integer :: i
      ! Get band limits
      error_msg = ''
      if (trim(band) == 'sw') then
         band_limits = k_dist_sw%get_band_lims_wavelength()
      else if (trim(band) == 'lw') then
         band_limits = k_dist_lw%get_band_lims_wavelength()
      else
         error_msg = trim(subname) // ': band ' // trim(band) // ' not known.'
         return
      end if
      ! Compute midpoints from band limits
      band_midpoints(:) = 0._wp
      do i = 1,size(band_midpoints)
         band_midpoints(i) = (band_limits(1,i) + band_limits(2,i)) / 2._wp
      end do
   end function get_band_midpoints
   !----------------------------------------------------------------------------

end module rrtmgp_interface
