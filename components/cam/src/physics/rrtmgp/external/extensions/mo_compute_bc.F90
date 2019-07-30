module mo_compute_bc
  ! This code is part of RRTM for GCM Applications - Parallel (RRTMGP)
  !
  ! Contacts: Robert Pincus and Eli Mlawer
  ! email:  rrtmgp@aer.com
  !
  ! Copyright 2018,  Atmospheric and Environmental Research and
  ! Regents of the University of Colorado.  All right reserved.
  !
  ! Use and duplication is permitted under the terms of the
  !    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
  ! -------------------------------------------------------------------------------------------------
  ! This modules lets users determine upper boundary condition by
  !   computing the spectrally-resolved fluxes at the bottom of an isothermal layer
  !   extending from the lowest supplied pressure to the minimum pressure allowed by
  !   RRTMGP.
  ! This is only sensible if the user's domain extends nearly to the top of the atmosphere.
  ! Adding this thin extra layer makes heating rates in the top-most layer more reasonable
  !   especially in the longwave
  ! The boundary condition is on diffuse flux in the LW and direct flux in the SW
  ! -------------------------------------------------------------------------------------------------
  use mo_rte_kind,           only: wp, wl
  use mo_source_functions,   only: ty_source_func_lw
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props, ty_optical_props_arry, &
                                   ty_optical_props_1scl, ty_optical_props_2str
  use mo_gas_optics,         only: ty_gas_optics
  use mo_fluxes,             only: ty_fluxes
  use mo_rte_lw,             only: rte_lw
  use mo_rte_sw,             only: rte_sw
  implicit none
  private
  public :: compute_bc

  !
  ! Extend ty_fluxes to report spectrally-resolved downwelling flux at a single layer
  !
  type, extends(ty_fluxes) :: ty_fluxes_1lev
    real(wp), dimension(:,:), pointer :: gpt_flux_dn => NULL()    ! (ncol, nlev, nband)
  contains
    procedure :: reduce => reduce_1lev
    procedure :: are_desired => are_desired_1lev
  end type ty_fluxes_1lev
contains
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! The arguments to this routine follow those to the gas_optics routines
  !
  function compute_bc(k_dist,                      &
                      play, plev, tlay, gas_concs, &
                      flux_bc, mu0) result(error_msg)
    class(ty_gas_optics),     intent(in   ) :: k_dist
    real(wp), dimension(:,:), intent(in   ) :: play, &    ! layer pressures [Pa, mb]; (ncol,nlay)
                                               plev, &    ! level pressures [Pa, mb]; (ncol,nlay+1)
                                               tlay       ! layer temperatures [K]; (ncol,nlay)
    type(ty_gas_concs),       intent(in   ) :: gas_concs  ! Gas volume mixing ratios
    real(wp), dimension(:,:), target, &
                              intent(  out) :: flux_bc    ! Boundary condition to be applied (ncol,ngpt)
    real(wp), dimension(:), optional, &
                              intent(in   ) :: mu0        ! Must be provided for solar problems
    character(len=128)                      :: error_msg
    ! ----------------------------------------------------------
    !
    ! Local variables
    !
    logical :: top_at_1
    integer :: ncol, nlay, ngpt

    integer :: igas, ngases
    character(len=32), dimension(:), allocatable :: gas_names
    real(wp), dimension(size(play,1), size(play,2)) :: vmr

    integer :: top_lay
    real(wp), dimension(size(play,1), 1) :: play_1lay, tlay_1lay
    real(wp), dimension(size(play,1), 2) :: plev_1lay, tlev_1lay
    real(wp), dimension(k_dist%get_nband(),size(play,1)) &
                                          :: lower_bc ! emissivity or surface albedo
    type(ty_gas_concs)                    :: gas_concs_1lay  ! Gas volume mixing ratios
    class(ty_optical_props_arry), &
                              allocatable :: optical_props_1lay
    type(ty_fluxes_1lev)                  :: fluxes_1lev
    type(ty_source_func_lw)               :: lw_sources_1lay
    real(wp), dimension(size(play,1),k_dist%get_nband()) :: solar_src
    ! ----------------------------------------------------------
    !
    ! Problem extent
    !
    ncol  = size(play, dim=1)
    nlay  = size(play, dim=2)
    ngpt  = k_dist%get_ngpt()
    if(any([size(plev,1) /= ncol, size(plev,2) /= nlay+1])) then
      error_msg = "compute_bc: array plev has wrong dimensions"
      return
    end if
    if(any([size(tlay,1) /= ncol, size(tlay,2) /= nlay  ])) then
      error_msg = "compute_bc: array tlay has wrong dimensions"
      return
    end if
    if(present(mu0)) then
      if(size(mu0) /= ncol) then
        error_msg  = "compute bc: array mu0 has wrong dimensions"
        return
      end if
    end if

    !
    ! Vertical ordering?
    !
    top_at_1 = play(1, 1) < play(1, nlay)
    top_lay = merge(1, nlay, top_at_1)
    if(any(plev(:,top_lay) <= &
           k_dist%get_press_min() + 2._wp * spacing(k_dist%get_press_min()))) then
      error_msg = "compute_bc: pressures are too close to (or less than) min in gas optics "
      return
    end if
    !
    ! Make a single-layer isothermal atmosphere
    !
    tlay_1lay(1:ncol,1) = tlay(1:ncol, top_lay)
    tlev_1lay(1:ncol,1) = tlay(1:ncol, top_lay)
    tlev_1lay(1:ncol,2) = tlay(1:ncol, top_lay)
    plev_1lay(1:ncol,1) = k_dist%get_press_min()
    plev_1lay(1:ncol,2) = plev(1:ncol, top_lay+1)
    !
    ! Maybe there are better ways to interpolate pressure but the single layer
    !   should be thin enough that interpolation doesn't have much impact
    !
    play_1lay(1:ncol,1) = 0.5 * (plev_1lay(1:ncol,1) + plev_1lay(1:ncol,2))

    !
    ! Gas concentrations in the single layer are the same as in the top layer
    !
    ngases = gas_concs%get_num_gases()
    allocate(gas_names(ngases))
    gas_names = gas_concs%get_gas_names()
    do igas = 1, ngases
      error_msg = gas_concs%get_vmr(gas_names(igas), vmr)
      if(error_msg /= "") return
      error_msg = gas_concs_1lay%set_vmr(gas_names(igas), vmr(:, top_lay))
      if(error_msg /= "") return
    end do

    lower_bc(:,:) = 1._wp ! Value doesn't affect downward flux
    fluxes_1lev%gpt_flux_dn => flux_bc
    ! ---------------------------------------------------
    if(k_dist%source_is_internal()) then
      !
      ! Longwave specific variables
      !
      allocate(ty_optical_props_1scl::optical_props_1lay)
      select type (optical_props_1lay)
        type is (ty_optical_props_1scl)
          error_msg =  optical_props_1lay%alloc_1scl(ncol, 1, k_dist)
          if(error_msg /= "") return
      end select
      error_msg = lw_sources_1lay%alloc(ncol, 1, k_dist)
      if(error_msg /= "") return
      !
      ! Gas optics and sources
      !
      error_msg = k_dist%gas_optics(play_1lay, plev_1lay,               &
                                    tlay_1lay, tlay_1lay(1:ncol,1),     &
                                    gas_concs_1lay, optical_props_1lay, &
                                    lw_sources_1lay, tlev = tlev_1lay)
      !                                                                  !
      ! Compute fluxes
      !
      error_msg = rte_lw(optical_props_1lay, &
                         top_at_1,           &
                         lw_sources_1lay,    &
                         lower_bc, fluxes_1lev)
    else
      !
      ! Shortwave specific variables
      !
      if(.not. present(mu0)) then
        error_msg = "compute_bc: have to supply mu0 for solar calculations"
        return
      end if
      allocate(ty_optical_props_2str::optical_props_1lay)
      select type (optical_props_1lay)
        type is (ty_optical_props_2str)
          error_msg =  optical_props_1lay%alloc_2str(ncol, 1, k_dist)
          if(error_msg /= "") return
      end select
      !
      ! Gas optics and sources
      !
      error_msg = k_dist%gas_optics(play_1lay, plev_1lay,       &
                                    tlay_1lay,  gas_concs_1lay, &
                                    optical_props_1lay,         &
                                    solar_src)
      error_msg = rte_sw(optical_props_1lay, &
                         top_at_1, mu0,      &
                         solar_src,          &
                         lower_bc, lower_bc, fluxes_1lev)
    endif
  end function
  ! --------------------------------------------------------------------------------------
  function reduce_1lev(this, gpt_flux_up, gpt_flux_dn, spectral_disc, top_at_1, gpt_flux_dn_dir) result(error_msg)
    class(ty_fluxes_1lev),             intent(inout) :: this
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_up ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    real(kind=wp), dimension(:,:,:),   intent(in   ) :: gpt_flux_dn ! Fluxes by gpoint [W/m2](ncol, nlay+1, ngpt)
    class(ty_optical_props),           intent(in   ) :: spectral_disc  !< derived type with spectral information
    logical,                           intent(in   ) :: top_at_1
    real(kind=wp), dimension(:,:,:), optional, &
                                       intent(in   ) :: gpt_flux_dn_dir! Direct flux down
    character(len=128)                               :: error_msg
    ! ------
    integer :: ncol, nlev, ngpt, bottom_lev
    ! ------
    error_msg = ""
    ncol = size(gpt_flux_up, DIM=1)
    nlev = size(gpt_flux_up, DIM=2)
    ngpt = size(gpt_flux_up, DIM=3)
    if(nlev /= 2) then
      error_msg = "reduce: expecting only two layers when computing boundary condition fluxes"
      return
    end if
    bottom_lev = merge(2, 1, top_at_1)
   !
   ! Return the g-point flux at the bottomw of a two-layer domain
   !
    if(associated(this%gpt_flux_dn)) then
      if(any([size(this%gpt_flux_dn, 1) /= ncol,  &
              size(this%gpt_flux_dn, 2) /= ngpt])) then
        error_msg = "reduce: gpt_flux_dn array incorrectly sized"
      else
        if(present(gpt_flux_dn_dir)) then
          this%gpt_flux_dn(:,:) = gpt_flux_dn_dir(:,bottom_lev,:)
        else
          this%gpt_flux_dn(:,:) = gpt_flux_dn    (:,bottom_lev,:)
        end if
      end if
    end if
  end function reduce_1lev
  ! --------------------------------------------------------------------------------------
  ! Are any fluxes desired from this set of g-point fluxes? We can tell because memory will
  !   be allocated for output
  !
  function are_desired_1lev(this)
    class(ty_fluxes_1lev), intent(in   ) :: this
    logical                              :: are_desired_1lev

    are_desired_1lev = associated(this%gpt_flux_dn)
  end function are_desired_1lev
end module mo_compute_bc
