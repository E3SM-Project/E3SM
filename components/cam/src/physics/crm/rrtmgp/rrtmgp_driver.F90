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
!
module rrtmgp_driver
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
  public :: rte_lw, rte_sw

  public :: k_dist_lw, k_dist_sw
  type(ty_gas_optics_rrtmgp) :: k_dist_lw, k_dist_sw

contains
  ! --------------------------------------------------
  !
  ! Interfaces using clear (gas + aerosol) and all-sky categories, starting from
  !   pressures, temperatures, and gas amounts for the gas contribution
  !
  ! --------------------------------------------------
  function rte_lw( &
        gas_names, gas_vmr, p_lay, t_lay, p_lev,            &
        t_sfc, sfc_emis, cld_tau, aer_tau,                 &
        allsky_fluxes, clrsky_fluxes,                      &
        col_dry, t_lev, inc_flux, n_gauss_angles) result(error_msg)
    character(len=*), dimension(:), intent(in) :: gas_names
    real(wp), dimension(:,:,:),    intent(in   ) :: gas_vmr
    real(wp), dimension(:,:),    intent(in   ) :: p_lay, t_lay !< pressure [Pa], temperature [K] at layer centers (ncol,nlay)
    real(wp), dimension(:,:),    intent(in   ) :: p_lev        !< pressure at levels/interfaces [Pa] (ncol,nlay+1)
    real(wp), dimension(:),      intent(in   ) :: t_sfc     !< surface temperature           [K]  (ncol)
    real(wp), dimension(:,:),    intent(in   ) :: sfc_emis  !< emissivity at surface         []   (nband, ncol)
    real(wp), dimension(:,:,:),  intent(in   ) :: cld_tau   ! cloud absorption optical depth (ncol,nlay,ngpt)
    real(wp), dimension(:,:,:),  intent(in   ) :: aer_tau   ! cloud absorption optical depth (ncol,nlay,nband)
    type(ty_fluxes_byband),      intent(inout) :: allsky_fluxes, clrsky_fluxes

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
    integer :: ncol, nlay, ngpt, nband, nstr, igas
    logical :: top_at_1
    ! --------------------------------
    ! Problem sizes
    !
    error_msg = ""

    ncol  = size(p_lay, 1)
    nlay  = size(p_lay, 2)
    ngpt  = k_dist_lw%get_ngpt()
    nband = k_dist_lw%get_nband()

    top_at_1 = p_lay(1, 1) < p_lay(1, nlay)

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

    ! Clear sky is gases + aerosols (if they're supplied)
    error_msg = aer_props%increment(optical_props)
    if(error_msg /= '') return
    error_msg = base_rte_lw(optical_props, top_at_1, sources, &
                            sfc_emis, clrsky_fluxes,          &
                            inc_flux, n_gauss_angles)
    if(error_msg /= '') return

    ! All-sky fluxes = clear skies + clouds
    error_msg = cld_props%increment(optical_props)
    if(error_msg /= '') return

    error_msg = base_rte_lw(optical_props, top_at_1, sources, &
                            sfc_emis, allsky_fluxes,          &
                            inc_flux, n_gauss_angles)

    call sources%finalize()
  end function rte_lw
  ! --------------------------------------------------
  function rte_sw( &
        gas_concs, p_lay, t_lay, p_lev,     &
        mu0, sfc_alb_dir, sfc_alb_dif, cloud_props, &
        allsky_fluxes, clrsky_fluxes,               &
        aer_props, col_dry, inc_flux, tsi_scaling   ) result(error_msg)
    type(ty_gas_concs),          intent(in   ) :: gas_concs    !< derived type encapsulating gas concentrations
    real(wp), dimension(:,:),    intent(in   ) :: p_lay, t_lay !< pressure [Pa], temperature [K] at layer centers (ncol,nlay)
    real(wp), dimension(:,:),    intent(in   ) :: p_lev        !< pressure at levels/interfaces [Pa] (ncol,nlay+1)
    real(wp), dimension(:  ),    intent(in   ) :: mu0          !< cosine of solar zenith angle
    !  surface albedo for direct and diffuse radiation (band, col)
    real(wp), dimension(:,:),    intent(in   ) :: sfc_alb_dir, sfc_alb_dif
    type(ty_optical_props_2str), intent(in   ) :: cloud_props !< cloud optical properties (ncol,nlay,ngpt)
    type(ty_fluxes_byband),      intent(inout) :: allsky_fluxes, clrsky_fluxes

    ! Optional inputs
    type(ty_optical_props_2str), &
              optional,       intent(in ) :: aer_props   !< aerosol optical properties
    real(wp), dimension(:,:), &
              optional,       intent(in ) :: col_dry, &  !< Molecular number density (ncol, nlay)
                                             inc_flux    !< incident flux at domain top [W/m2] (ncol, ngpts)
    real(wp), optional,       intent(in ) :: tsi_scaling !< Optional scaling for total solar irradiance

    character(len=128)                    :: error_msg
    ! --------------------------------
    ! Local variables
    !
    type(ty_optical_props_2str) :: optical_props
    real(wp), dimension(:,:), allocatable :: toa_flux
    integer :: ncol, nlay, ngpt, nband, nstr
    logical :: top_at_1
    ! --------------------------------
    ! Problem sizes
    !
    error_msg = ""

    ncol  = size(p_lay, 1)
    nlay  = size(p_lay, 2)
    ngpt  = k_dist_sw%get_ngpt()
    nband = k_dist_sw%get_nband()

    top_at_1 = p_lay(1, 1) < p_lay(1, nlay)

    ! ------------------------------------------------------------------------------------
    !  Error checking
    !
    if(present(inc_flux) .and. present(tsi_scaling)) &
      error_msg = "rrtmpg_sw: only one of inc_flux, tsi_scaling may be supplied."

    if(present(aer_props)) then
      if(any([aer_props%get_ncol(), &
              aer_props%get_nlay()] /= [ncol, nlay])) &
        error_msg = "rrtmpg_sw: aerosol properties inconsistently sized"
      if(.not. any(aer_props%get_ngpt() /= [ngpt, nband])) &
        error_msg = "rrtmpg_sw: aerosol properties inconsistently sized"
    end if

    if(present(tsi_scaling)) then
      if(tsi_scaling <= 0._wp) &
        error_msg = "rrtmpg_sw: tsi_scaling is < 0"
    end if

    if(present(inc_flux)) then
      if(any([size(inc_flux, 1), &
              size(inc_flux, 2)] /= [ncol, ngpt])) &
        error_msg = "rrtmpg_sw: incident flux inconsistently sized"
    end if
    if(len_trim(error_msg) > 0) return

    ! Optical properties arrays
    error_msg = optical_props%alloc_2str(ncol, nlay, k_dist_sw)
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
    if(present(aer_props)) error_msg = aer_props%increment(optical_props)
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

end module rrtmgp_driver
