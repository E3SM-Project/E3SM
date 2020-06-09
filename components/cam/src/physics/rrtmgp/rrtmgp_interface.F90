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

! Macro to handle string return codes from RRTMGP functions; if string is
! non-empty, we print the string to stdout, set a return value to 1, and return
! from the function. We do this as a macro because calling return from a fortran
! function would simply return from our error handler, not from the caller.
#define handle_error(error_msg) if (trim(error_msg) /= '') then; print *, trim(error_msg); error = 1; return; endif

module rrtmgp_interface
  use mo_rte_kind,   only: wp
#ifdef RRTMGPXX

  implicit none

  interface

    subroutine clear_gas_names() bind(C,name="clear_gas_names")
    end subroutine

    subroutine push_gas_name(gas_name) bind(C,name="push_gas_name")
      use iso_c_binding
      character(c_char), dimension(*), intent(in) :: gas_name
    end subroutine

    subroutine gas_names_to_string1d() bind(C,name="gas_names_to_string1d")
    end subroutine

    subroutine rrtmgp_initialize_cpp(coefficients_file_sw, coefficients_file_lw) bind(C,name="rrtmgp_initialize_cpp")
      use iso_c_binding
      character(c_char), dimension(*), intent(in) :: coefficients_file_sw
      character(c_char), dimension(*), intent(in) :: coefficients_file_lw
    end subroutine

    ! Called directly from radiation.F90
    function rrtmgp_finalize() result(ret) bind(C,name="rrtmgp_finalize")
      integer :: ret
    end function

    function get_nband_sw() result(ret) bind(C,name="get_nband_sw")
      integer :: ret
    end function

    function get_nband_lw() result(ret) bind(C,name="get_nband_lw")
      integer :: ret
    end function

    function get_ngpt_sw() result(ret) bind(C,name="get_ngpt_sw")
      integer :: ret
    end function

    function get_ngpt_lw() result(ret) bind(C,name="get_ngpt_lw")
      integer :: ret
    end function

    subroutine get_band_lims_wavenumber_sw(nbnd, band_limits_p) bind(C,name="get_band_lims_wavenumber_sw")
      integer, value :: nbnd
      real(8), dimension(*) :: band_limits_p
    end subroutine

    subroutine get_band_lims_wavenumber_lw(nbnd, band_limits_p) bind(C,name="get_band_lims_wavenumber_lw")
      integer, value :: nbnd
      real(8), dimension(*) :: band_limits_p
    end subroutine

    ! Called directly from radiation.F90
    function get_temp_min() result(ret) bind(C,name="get_temp_min")
      real(8) :: ret
    end function

    ! Called directly from radiation.F90
    function get_temp_max() result(ret) bind(C,name="get_temp_max")
      real(8) :: ret
    end function

    subroutine get_gpoint_bands_sw(nbnd, band_limits_p) bind(C,name="get_gpoint_bands_sw")
      integer, value :: nbnd
      integer, dimension(*) :: band_limits_p
    end subroutine

    subroutine get_gpoint_bands_lw(nbnd, band_limits_p) bind(C,name="get_gpoint_bands_lw")
      integer, value :: nbnd
      integer, dimension(*) :: band_limits_p
    end subroutine

    subroutine get_band_midpoints_sw(nbnd, band_midpoints_p) bind(C,name="get_band_midpoints_sw")
      integer, value :: nbnd
      real(8), dimension(*) :: band_midpoints_p
    end subroutine

    subroutine get_band_midpoints_lw(nbnd, band_midpoints_p) bind(C,name="get_band_midpoints_lw")
      integer, value :: nbnd
      real(8), dimension(*) :: band_midpoints_p
    end subroutine

    subroutine rrtmgp_run_lw_cpp(ngas, ncol, nlay, nbnd, ngpt, gas_vmr_p, p_lay_p, t_lay_p, p_lev_p, t_sfc_p,             &
                                 sfc_emis_p, cld_tau_p, aer_tau_p, allsky_flux_up_p, allsky_flux_dn_p, allsky_flux_net_p, &
                                 allsky_bnd_flux_up_p, allsky_bnd_flux_dn_p, allsky_bnd_flux_net_p, clrsky_flux_up_p,     &
                                 clrsky_flux_dn_p, clrsky_flux_net_p, clrsky_bnd_flux_up_p, clrsky_bnd_flux_dn_p,         &
                                 clrsky_bnd_flux_net_p, t_lev_p, n_gauss_angles) bind(C,name="rrtmgp_run_lw_cpp")
      integer, value :: ngas, ncol, nlay, nbnd, ngpt, n_gauss_angles
      real(8), dimension(*) :: gas_vmr_p, p_lay_p, t_lay_p, p_lev_p, t_sfc_p, sfc_emis_p, cld_tau_p, aer_tau_p,      &
                               allsky_flux_up_p, allsky_flux_dn_p, allsky_flux_net_p, allsky_bnd_flux_up_p,          &
                               allsky_bnd_flux_dn_p, allsky_bnd_flux_net_p, clrsky_flux_up_p, clrsky_flux_dn_p,      &
                               clrsky_flux_net_p, clrsky_bnd_flux_up_p, clrsky_bnd_flux_dn_p, clrsky_bnd_flux_net_p, &
                               t_lev_p
    end subroutine


    subroutine rrtmgp_run_sw_cpp(ngas, ncol, nlay, nbnd, ngpt, gas_vmr_p, p_lay_p, t_lay_p, p_lev_p,                          &
                                 mu0_p, sfc_alb_dir_p, sfc_alb_dif_p, cld_tau_p, cld_ssa_p, cld_asm_p,                        &
                                 aer_tau_p, aer_ssa_p, aer_asm_p, allsky_flux_up_p, allsky_flux_dn_p,                         &
                                 allsky_flux_net_p, allsky_bnd_flux_up_p, allsky_bnd_flux_dn_p, allsky_bnd_flux_net_p,        &
                                 allsky_bnd_flux_dn_dir_p, clrsky_flux_up_p, clrsky_flux_dn_p, clrsky_flux_net_p,             &
                                 clrsky_bnd_flux_up_p, clrsky_bnd_flux_dn_p, clrsky_bnd_flux_net_p, clrsky_bnd_flux_dn_dir_p, &
                                 tsi_scaling) bind(C,name="rrtmgp_run_sw_cpp")
      integer, value :: ngas, ncol, nlay, nbnd, ngpt
      real(8), value :: tsi_scaling
      real(8), dimension(*) :: gas_vmr_p, p_lay_p, t_lay_p, p_lev_p, mu0_p, sfc_alb_dir_p, sfc_alb_dif_p, cld_tau_p, cld_ssa_p, &
                               cld_asm_p, aer_tau_p, aer_ssa_p, aer_asm_p, allsky_flux_up_p, allsky_flux_dn_p,                  &
                               allsky_flux_net_p, allsky_bnd_flux_up_p, allsky_bnd_flux_dn_p, allsky_bnd_flux_net_p,            &
                               allsky_bnd_flux_dn_dir_p, clrsky_flux_up_p, clrsky_flux_dn_p, clrsky_flux_net_p,                 &
                               clrsky_bnd_flux_up_p, clrsky_bnd_flux_dn_p, clrsky_bnd_flux_net_p, clrsky_bnd_flux_dn_dir_p
    end subroutine

  end interface

contains

  function to_c_char(str) result(ret)
    use iso_c_binding
    character(len=*) :: str
    character(kind=c_char,len=len_trim(str)+1) :: ret
    ret = trim(str)//char(0)
  end function

  integer function rrtmgp_initialize(ngas, coefficients_file_sw, coefficients_file_lw, active_gases) result(error)
    integer, intent(in) :: ngas
    character(len=*), intent(in) :: coefficients_file_sw, coefficients_file_lw
    character(len=*), dimension(ngas), intent(in) :: active_gases
    integer :: igas
    call clear_gas_names()
    do igas = 1 , ngas
      call push_gas_name(to_c_char(active_gases(igas)))
    enddo
    call gas_names_to_string1d()
    call rrtmgp_initialize_cpp(to_c_char(coefficients_file_sw), to_c_char(coefficients_file_lw))
    error = 0
  end function

  integer function get_nband(band)
    character(len=*), intent(in) :: band
    if (trim(band) == 'sw') then
      get_nband = get_nband_sw()
    else if (trim(band) == 'lw') then
      get_nband = get_nband_lw()
    else
      get_nband = -1
    end if
  end function get_nband

  integer function get_ngpt(band)
    character(len=*), intent(in) :: band
    if (trim(band) == 'sw') then
      get_ngpt = get_ngpt_sw()
    else if (trim(band) == 'lw') then
      get_ngpt = get_ngpt_lw()
    else
      get_ngpt = -1
    end if
  end function get_ngpt

  integer function get_band_lims_wavenumber(nbnd, band, band_limits) result(error)
    integer, intent(in) :: nbnd
    character(len=*), intent(in) :: band
    real(wp), dimension(2,nbnd), intent(out) :: band_limits

    ! Initialize error status to 0
    error = 0

    if (trim(band) == 'sw') then
       call get_band_lims_wavenumber_sw(nbnd, band_limits)
    else if (trim(band) == 'lw') then
       call get_band_lims_wavenumber_lw(nbnd, band_limits)
    else
       band_limits = -1
       error = 1
       return
    end if
  end function get_band_lims_wavenumber

  integer function get_gpoint_bands(ngpt, band, gpoint_bands) result(error)
    integer, intent(in) :: ngpt
    character(len=*), intent(in) :: band
    integer, dimension(ngpt), intent(out) :: gpoint_bands

    ! Initialize error status to 0
    error = 0

    if (trim(band) == 'sw') then
      call get_gpoint_bands_sw(ngpt, gpoint_bands)
    else if (trim(band) == 'lw') then
      call get_gpoint_bands_lw(ngpt, gpoint_bands)
    else
      gpoint_bands = -1
      error = 1
      return
    end if
  end function

  integer function get_band_midpoints(nband, band, band_midpoints) result(error)
    integer, intent(in) :: nband
    character(len=*), intent(in) :: band
    real(wp), intent(inout) :: band_midpoints(nband)
    real(wp) :: band_limits(2,nband)
    character(len=128), parameter :: subname = 'get_band_midpoints'
    integer :: i

    ! Initialize error status to 0
    error = 0

    ! Get band limits
    if (trim(band) == 'sw') then
       call get_band_midpoints_sw(nband, band_midpoints)
    else if (trim(band) == 'lw') then
       call get_band_midpoints_lw(nband, band_midpoints)
    else
       print *, trim(subname) // ': band ' // trim(band) // ' not known.'
       error = 1
       return
    end if
  end function get_band_midpoints

  integer function rrtmgp_run_lw(ngas, ncol, nlay, nbnd, ngpt, &
                                 gas_names, gas_vmr, p_lay, t_lay, p_lev,    &
                                 t_sfc, sfc_emis, cld_tau, aer_tau,          &
                                 allsky_flux_up, allsky_flux_dn, allsky_flux_net, &
                                 allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, &
                                 clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net, &
                                 clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net, &
                                 col_dry, t_lev, inc_flux, n_gauss_angles) result(error)
    integer, intent(in) :: ngas, ncol, nlay, nbnd, ngpt
    character(len=*), dimension(ngas), intent(in) :: gas_names
    real(wp), dimension(ngas,ncol,nlay), intent(in   ) :: gas_vmr
    real(wp), dimension(ncol,nlay),      intent(in   ) :: p_lay, t_lay !< pressure [Pa], temperature [K] at layer centers (ncol,nlay)
    real(wp), dimension(ncol,nlay+1),    intent(in   ) :: p_lev        !< pressure at levels/interfaces [Pa] (ncol,nlay+1)
    real(wp), dimension(ncol),           intent(in   ) :: t_sfc     !< surface temperature           [K]  (ncol)
    real(wp), dimension(nbnd,ncol),      intent(in   ) :: sfc_emis  !< emissivity at surface         []   (nband, ncol)
    real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: cld_tau   ! cloud absorption optical depth (ncol,nlay,ngpt)
    real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: aer_tau   ! cloud absorption optical depth (ncol,nlay,nband)
    real(wp), dimension(ncol,nlay+1), intent(inout), target :: allsky_flux_up, allsky_flux_dn, allsky_flux_net, &
                                                               clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net
    real(wp), dimension(ncol,nlay+1,nbnd), intent(inout), target :: allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, &
                                                                    clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net
    ! Optional inputs
    real(wp), dimension(ncol,nlay)          , optional, intent(in ) :: col_dry !< Molecular number density (ncol, nlay)
    real(wp), dimension(ncol,nlay+1), target, optional, intent(in ) :: t_lev     !< temperature at levels [K] (ncol, nlay+1)
    real(wp), dimension(ncol,ngpt)  , target, optional, intent(in ) :: inc_flux   !< incident flux at domain top [W/m2] (ncol, ngpts)
    integer                                 , optional, intent(in ) :: n_gauss_angles ! Number of angles used in Gaussian quadrature (no-scattering solution)
    integer :: igas
    call clear_gas_names()
    do igas = 1 , ngas
      call push_gas_name(to_c_char(gas_names(igas)))
    enddo
    call gas_names_to_string1d()
    call rrtmgp_run_lw_cpp(ngas, ncol, nlay, nbnd, ngpt, gas_vmr, p_lay, t_lay, p_lev, t_sfc,           &
                           sfc_emis, cld_tau, aer_tau, allsky_flux_up, allsky_flux_dn, allsky_flux_net, &
                           allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, clrsky_flux_up, &
                           clrsky_flux_dn, clrsky_flux_net, clrsky_bnd_flux_up, clrsky_bnd_flux_dn,     &
                           clrsky_bnd_flux_net, t_lev, n_gauss_angles)

    error = 0;
  end function rrtmgp_run_lw

  integer function rrtmgp_run_sw(ngas, ncol, nlay, nbnd, ngpt, &
                                 gas_names, gas_vmr, p_lay, t_lay, p_lev, &
                                 mu0, sfc_alb_dir, sfc_alb_dif,  &
                                 cld_tau, cld_ssa, cld_asm,      &
                                 aer_tau, aer_ssa, aer_asm,      &
                                 allsky_flux_up, allsky_flux_dn, allsky_flux_net, &
                                 allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, allsky_bnd_flux_dn_dir, &
                                 clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net, &
                                 clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net, clrsky_bnd_flux_dn_dir, &
                                 col_dry, inc_flux, tsi_scaling  ) result(error)
    integer, intent(in) :: ngas, ncol, nlay, nbnd, ngpt
    character(len=*), dimension(ngas),     intent(in   ) :: gas_names
    real(wp), dimension(ngas,ncol,nlay),   intent(in   ) :: gas_vmr
    real(wp), dimension(ncol,nlay),        intent(in   ) :: p_lay, t_lay !< pressure [Pa], temperature [K] at layer centers (ncol,nlay)
    real(wp), dimension(ncol,nlay+1),      intent(in   ) :: p_lev        !< pressure at levels/interfaces [Pa] (ncol,nlay+1)
    real(wp), dimension(ncol),             intent(in   ) :: mu0          !< cosine of solar zenith angle
    !  surface albedo for direct and diffuse radiation (band, col)
    real(wp), dimension(nbnd,ncol),        intent(in   ) :: sfc_alb_dir, sfc_alb_dif
    real(wp), dimension(ncol,nlay,ngpt),   intent(in   ) :: cld_tau, cld_ssa, cld_asm
    real(wp), dimension(ncol,nlay,nbnd),   intent(in   ) :: aer_tau, aer_ssa, aer_asm
    real(wp), dimension(ncol,nlay+1),      intent(inout), target :: allsky_flux_up, allsky_flux_dn, allsky_flux_net, &
                                                                    clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net
    real(wp), dimension(ncol,nlay+1,nbnd), intent(inout), target :: allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, allsky_bnd_flux_dn_dir, &
                                                                    clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net, clrsky_bnd_flux_dn_dir
    ! Optional inputs
    real(wp), dimension(ncol,nlay), optional, intent(in ) :: col_dry     !< Molecular number density (ncol, nlay)
    real(wp), dimension(ncol,ngpt), optional, intent(in ) :: inc_flux    !< incident flux at domain top [W/m2] (ncol, ngpts)
    real(wp)                      , optional, intent(in ) :: tsi_scaling !< Optional scaling for total solar irradiance
    integer :: igas
    call clear_gas_names()
    do igas = 1 , ngas
      call push_gas_name(to_c_char(gas_names(igas)))
    enddo
    call gas_names_to_string1d()
    call rrtmgp_run_sw_cpp(ngas, ncol, nlay, nbnd, ngpt, gas_vmr, p_lay, t_lay, p_lev,                          &
                           mu0, sfc_alb_dir, sfc_alb_dif, cld_tau, cld_ssa, cld_asm,                            &
                           aer_tau, aer_ssa, aer_asm, allsky_flux_up, allsky_flux_dn,                           &
                           allsky_flux_net, allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net,        &
                           allsky_bnd_flux_dn_dir, clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net,             &
                           clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net, clrsky_bnd_flux_dn_dir, &
                           tsi_scaling)

    error = 0
  end function rrtmgp_run_sw 



#else  /*ifdef RRTMGPXX*/



  use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props, only: ty_optical_props_1scl, ty_optical_props_2str
  use mo_source_functions, only: ty_source_func_lw
  use mo_fluxes_byband, only: ty_fluxes_byband
  use mo_rte_lw, only: rte_lw
  use mo_rte_sw, only: rte_sw
  implicit none
  private
  public :: &
     rrtmgp_initialize, rrtmgp_finalize, &
     rrtmgp_run_lw, rrtmgp_run_sw, &
     get_nband, get_ngpt, get_band_lims_wavenumber, get_band_midpoints, &
     get_temp_min, get_temp_max, get_gpoint_bands

  public :: k_dist_lw, k_dist_sw
  type(ty_gas_optics_rrtmgp) :: k_dist_lw, k_dist_sw
  integer, parameter :: ERROR_MSG_LEN = 128

contains
   !----------------------------------------------------------------------------
   integer function rrtmgp_initialize(ngas, coefficients_file_sw, coefficients_file_lw, active_gases) result(error)
      use rrtmgp_coefficients, only: load_and_init
      use mo_gas_concentrations, only: ty_gas_concs
      integer, intent(in) :: ngas
      character(len=*), intent(in) :: coefficients_file_sw, coefficients_file_lw
      character(len=*), dimension(ngas), intent(in) :: active_gases
      ! ty_gas_concs object that would normally hold volume mixing ratios for
      ! radiatively-important gases. Here, this is just used to provide the names
      ! of gases that are available in the model (needed by the kdist
      ! initialization routines that are called within the load_coefficients
      ! methods).
      type(ty_gas_concs) :: available_gases

      ! Initialize error status to 0
      error = 0

      ! Need to initialize available_gases here! The only field of the
      ! available_gases type that is used in the kdist initialize is
      ! available_gases%gas_name, which gives the name of each gas that would be
      ! present in the ty_gas_concs object. So, we can just set this here, rather
      ! than trying to fully populate the ty_gas_concs object here, which would be
      ! impossible from this initialization routine because I do not think the
      ! rad_cnst objects are setup yet.
      ! TODO: This needs to be fixed to ONLY read in the data if masterproc, and then broadcast
      ! the other tasks!
      if (set_available_gases(ngas, active_gases, available_gases) /= 0) then
         error = 1
         return
      end if

      ! Read coefficients from file and populate k-dist objects
      call load_and_init(k_dist_sw, coefficients_file_sw, available_gases)
      call load_and_init(k_dist_lw, coefficients_file_lw, available_gases)

   end function rrtmgp_initialize
   !----------------------------------------------------------------------------
   integer function set_available_gases(ngas, gases, gas_concentrations) result(error)

      use mo_gas_concentrations, only: ty_gas_concs
      use mo_rrtmgp_util_string, only: lower_case

      integer, intent(in) :: ngas
      type(ty_gas_concs), intent(inout) :: gas_concentrations
      character(len=*), intent(in) :: gases(ngas)
      character(len=32), dimension(ngas) :: gases_lowercase
      integer :: igas

      ! Initialize error status to 0
      error = 0

      ! Initialize with lowercase gas names; we should work in lowercase
      ! whenever possible because we cannot trust string comparisons in RRTMGP
      ! to be case insensitive
      do igas = 1,ngas
         gases_lowercase(igas) = trim(lower_case(gases(igas)))
      end do
      handle_error(gas_concentrations%init(gases_lowercase))

   end function set_available_gases
   !----------------------------------------------------------------------------
   integer function rrtmgp_finalize() result(error)
      error = 0
   end function rrtmgp_finalize
   !----------------------------------------------------------------------------
   integer function get_nband(band)
      character(len=*), intent(in) :: band
      if (trim(band) == 'sw') then
         get_nband = k_dist_sw%get_nband()
      else if (trim(band) == 'lw') then
         get_nband = k_dist_lw%get_nband()
      else
         get_nband = -1
      end if
   end function get_nband
   !----------------------------------------------------------------------------
   integer function get_ngpt(band)
      character(len=*), intent(in) :: band
      if (trim(band) == 'sw') then
         get_ngpt = k_dist_sw%get_ngpt()
      else if (trim(band) == 'lw') then
         get_ngpt = k_dist_lw%get_ngpt()
      else
         get_ngpt = -1
      end if
   end function get_ngpt
   !----------------------------------------------------------------------------
   integer function get_band_lims_wavenumber(nbnd, band, band_limits) result(error)
      integer, intent(in) :: nbnd
      character(len=*), intent(in) :: band
      real(wp), dimension(2,nbnd), intent(out) :: band_limits

      ! Initialize error status to 0
      error = 0

      if (trim(band) == 'sw') then
         band_limits = k_dist_sw%get_band_lims_wavenumber()
      else if (trim(band) == 'lw') then
         band_limits = k_dist_lw%get_band_lims_wavenumber()
      else
         band_limits = -1
         error = 1
         return
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
   integer function get_gpoint_bands(ngpt, band, gpoint_bands) result(error)
      integer, intent(in) :: ngpt
      character(len=*), intent(in) :: band
      integer, dimension(ngpt), intent(out) :: gpoint_bands

      ! Initialize error status to 0
      error = 0

      if (trim(band) == 'sw') then
         gpoint_bands = k_dist_sw%get_gpoint_bands()
      else if (trim(band) == 'lw') then
         gpoint_bands = k_dist_lw%get_gpoint_bands()
      else
         gpoint_bands = -1
         error = 1
         return
      end if
   end function

   !----------------------------------------------------------------------------

   ! Function to calculate band midpoints from kdist band limits
   integer function get_band_midpoints(nband, band, band_midpoints) result(error)
      integer, intent(in) :: nband
      character(len=*), intent(in) :: band
      real(wp), intent(inout) :: band_midpoints(nband)
      real(wp) :: band_limits(2,nband)
      character(len=128), parameter :: subname = 'get_band_midpoints'
      integer :: i

      ! Initialize error status to 0
      error = 0

      ! Get band limits
      if (trim(band) == 'sw') then
         band_limits = k_dist_sw%get_band_lims_wavelength()
      else if (trim(band) == 'lw') then
         band_limits = k_dist_lw%get_band_lims_wavelength()
      else
         print *, trim(subname) // ': band ' // trim(band) // ' not known.'
         error = 1
         return
      end if

      ! Compute midpoints from band limits
      band_midpoints(:) = 0._wp
      do i = 1,nband
         band_midpoints(i) = (band_limits(1,i) + band_limits(2,i)) / 2._wp
      end do
   end function get_band_midpoints

   !----------------------------------------------------------------------------

   integer function rrtmgp_run_lw( &
        ngas, ncol, nlay, nbnd, ngpt, &
        gas_names, gas_vmr, p_lay, t_lay, p_lev,    &
        t_sfc, sfc_emis, cld_tau, aer_tau,          &
        allsky_flux_up, allsky_flux_dn, allsky_flux_net, &
        allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, &
        clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net, &
        clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net, &
        col_dry, t_lev, inc_flux, n_gauss_angles) result(error)
      integer, intent(in) :: ngas, ncol, nlay, nbnd, ngpt
      character(len=*), dimension(ngas), intent(in) :: gas_names
      real(wp), dimension(ngas,ncol,nlay), intent(in   ) :: gas_vmr
      real(wp), dimension(ncol,nlay),      intent(in   ) :: p_lay, t_lay !< pressure [Pa], temperature [K] at layer centers (ncol,nlay)
      real(wp), dimension(ncol,nlay+1),    intent(in   ) :: p_lev        !< pressure at levels/interfaces [Pa] (ncol,nlay+1)
      real(wp), dimension(ncol),           intent(in   ) :: t_sfc     !< surface temperature           [K]  (ncol)
      real(wp), dimension(nbnd,ncol),      intent(in   ) :: sfc_emis  !< emissivity at surface         []   (nband, ncol)
      real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: cld_tau   ! cloud absorption optical depth (ncol,nlay,ngpt)
      real(wp), dimension(ncol,nlay,nbnd),  intent(in   ) :: aer_tau   ! cloud absorption optical depth (ncol,nlay,nband)
      real(wp), dimension(ncol,nlay+1), intent(inout), target :: &
         allsky_flux_up, allsky_flux_dn, allsky_flux_net, &
         clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net
      real(wp), dimension(ncol,nlay+1,nbnd), intent(inout), target :: &
         allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, &
         clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net

      ! Optional inputs
      real(wp), dimension(ncol,nlay), &
                optional,       intent(in ) :: col_dry !< Molecular number density (ncol, nlay)
      real(wp), dimension(ncol,nlay+1), target, &
                optional,       intent(in ) :: t_lev     !< temperature at levels [K] (ncol, nlay+1)
      real(wp), dimension(ncol,ngpt), target, &
                optional,       intent(in ) :: inc_flux   !< incident flux at domain top [W/m2] (ncol, ngpts)
      integer,  optional,       intent(in ) :: n_gauss_angles ! Number of angles used in Gaussian quadrature (no-scattering solution)
  
      ! Local variables
      type(ty_gas_concs) :: gas_concs    !< derived type encapsulating gas concentrations
      type(ty_optical_props_1scl) :: optical_props, cld_props, aer_props
      type(ty_source_func_lw) :: sources
      type(ty_fluxes_byband) :: allsky_fluxes, clrsky_fluxes
      integer :: igas
      logical :: top_at_1
      character(128) :: error_msg

      ! Initialize error status to 0
      error = 0
  
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
      handle_error(gas_concs%init(gas_names))
      do igas = 1,ngas
         handle_error(gas_concs%set_vmr(gas_names(igas), gas_vmr(igas,:,:)))
      end do
  
      ! Optical properties arrays
      handle_error(optical_props%alloc_1scl(ncol, nlay, k_dist_lw))
      handle_error(cld_props%alloc_1scl(ncol, nlay, k_dist_lw))
      cld_props%tau = cld_tau
      handle_error(aer_props%alloc_1scl(ncol, nlay, k_dist_lw%get_band_lims_wavenumber()))
      aer_props%tau = aer_tau
  
      ! Source function
      handle_error(sources%init(k_dist_lw))
      handle_error(sources%alloc(ncol, nlay))
  
      ! Gas optical depth -- pressure need to be expressed as Pa
      error_msg = k_dist_lw%gas_optics(        &
         p_lay, p_lev, t_lay, t_sfc, gas_concs, &
         optical_props, sources,                &
         col_dry, t_lev                         &
      )
      handle_error(error_msg)
  
      ! Clear sky fluxes (gases + aerosols)
      handle_error(aer_props%increment(optical_props))
      error_msg = rte_lw(            &
         optical_props, top_at_1, sources, &
         sfc_emis, clrsky_fluxes,          &
         inc_flux, n_gauss_angles          &
      )
      handle_error(error_msg)
  
      ! All-sky fluxes (clear skies + clouds)
      handle_error(cld_props%increment(optical_props))
      error_msg = rte_lw(            &
         optical_props, top_at_1, sources, &
         sfc_emis, allsky_fluxes,          &
         inc_flux, n_gauss_angles          &
      )
      handle_error(error_msg)
  
      ! Clean up
      call sources%finalize()
   end function rrtmgp_run_lw
   ! --------------------------------------------------
   integer function rrtmgp_run_sw( &
          ngas, ncol, nlay, nbnd, ngpt, &
          gas_names, gas_vmr, p_lay, t_lay, p_lev, &
          mu0, sfc_alb_dir, sfc_alb_dif,  &
          cld_tau, cld_ssa, cld_asm,      &
          aer_tau, aer_ssa, aer_asm,      &
          allsky_flux_up, allsky_flux_dn, allsky_flux_net, &
          allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, allsky_bnd_flux_dn_dir, &
          clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net, &
          clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net, clrsky_bnd_flux_dn_dir, &
          col_dry, inc_flux, tsi_scaling  ) result(error)
      integer, intent(in) :: ngas, ncol, nlay, nbnd, ngpt
      character(len=*), dimension(ngas), intent(in   ) :: gas_names
      real(wp), dimension(ngas,ncol,nlay), intent(in   ) :: gas_vmr
      real(wp), dimension(ncol,nlay),      intent(in   ) :: p_lay, t_lay !< pressure [Pa], temperature [K] at layer centers (ncol,nlay)
      real(wp), dimension(ncol,nlay+1),    intent(in   ) :: p_lev        !< pressure at levels/interfaces [Pa] (ncol,nlay+1)
      real(wp), dimension(ncol),           intent(in   ) :: mu0          !< cosine of solar zenith angle
      !  surface albedo for direct and diffuse radiation (band, col)
      real(wp), dimension(nbnd,ncol),      intent(in   ) :: sfc_alb_dir, sfc_alb_dif
      real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: cld_tau, cld_ssa, cld_asm
      real(wp), dimension(ncol,nlay,nbnd), intent(in   ) :: aer_tau, aer_ssa, aer_asm
      real(wp), dimension(ncol,nlay+1),      intent(inout), target :: &
         allsky_flux_up, allsky_flux_dn, allsky_flux_net, &
         clrsky_flux_up, clrsky_flux_dn, clrsky_flux_net
      real(wp), dimension(ncol,nlay+1,nbnd), intent(inout), target :: &
         allsky_bnd_flux_up, allsky_bnd_flux_dn, allsky_bnd_flux_net, allsky_bnd_flux_dn_dir, &
         clrsky_bnd_flux_up, clrsky_bnd_flux_dn, clrsky_bnd_flux_net, clrsky_bnd_flux_dn_dir

      ! Optional inputs
      real(wp), dimension(ncol,nlay), &
                optional,       intent(in ) :: col_dry     !< Molecular number density (ncol, nlay)
      real(wp), dimension(ncol,ngpt), &
                optional,       intent(in ) :: inc_flux    !< incident flux at domain top [W/m2] (ncol, ngpts)
      real(wp), optional,       intent(in ) :: tsi_scaling !< Optional scaling for total solar irradiance

      ! Local variables
      type(ty_optical_props_2str) :: optical_props, cloud_props, aer_props
      type(ty_gas_concs) :: gas_concs    !< derived type encapsulating gas concentrations
      type(ty_fluxes_byband) :: allsky_fluxes, clrsky_fluxes
      real(wp), dimension(:,:), allocatable :: toa_flux
      integer :: igas
      logical :: top_at_1
      character(len=128) :: error_msg

      ! Initialize error status to 0
      error = 0

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
      handle_error(gas_concs%init(gas_names))
      do igas = 1,ngas
         handle_error(gas_concs%set_vmr(gas_names(igas), gas_vmr(igas,:,:)))
      end do

      ! Optical properties arrays
      handle_error(optical_props%alloc_2str(ncol, nlay, k_dist_sw))
      handle_error(cloud_props%alloc_2str(ncol, nlay, k_dist_sw))
      handle_error(aer_props%alloc_2str(ncol, nlay, k_dist_sw%get_band_lims_wavenumber()))
  
      ! Populate optical properties
      cloud_props%tau = cld_tau
      cloud_props%ssa = cld_ssa
      cloud_props%g   = cld_asm
      aer_props%tau = aer_tau
      aer_props%ssa = aer_ssa
      aer_props%g   = aer_asm

      ! Delta scale
      handle_error(cloud_props%delta_scale())
  
      ! TOA flux
      allocate(toa_flux(ncol, ngpt))
  
      ! Gas optical depth -- pressure need to be expressed as Pa
      error_msg = k_dist_sw%gas_optics( &
         p_lay, p_lev, t_lay, gas_concs, &
         optical_props, toa_flux,        &
         col_dry                         &
      )
      handle_error(error_msg)
  
      ! If users have supplied an incident flux, use that
      if(present(inc_flux))    toa_flux(:,:) = inc_flux(:,:)
      if(present(tsi_scaling)) toa_flux(:,:) = toa_flux(:,:) * tsi_scaling
  
      ! Clear sky is gases + aerosols (if they're supplied)
      handle_error(aer_props%increment(optical_props))
      error_msg = rte_sw(optical_props, top_at_1, &
                                 mu0, toa_flux,                   &
                                 sfc_alb_dir, sfc_alb_dif,        &
                                 clrsky_fluxes)
      handle_error(error_msg)
  
      ! All-sky fluxes = clear skies + clouds
      handle_error(cloud_props%increment(optical_props))
      error_msg = rte_sw(    &
         optical_props, top_at_1,  &
         mu0, toa_flux,            &
         sfc_alb_dir, sfc_alb_dif, &
         allsky_fluxes             &
      )
      handle_error(error_msg)
   end function rrtmgp_run_sw 

#endif  /*ifdef RRTMGPXX*/

end module rrtmgp_interface
