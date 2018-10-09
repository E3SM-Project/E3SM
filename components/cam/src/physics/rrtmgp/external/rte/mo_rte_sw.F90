! This code is part of Radiative Transfer for Energetics (RTE)
!
! Contacts: Robert Pincus and Eli Mlawer
! email:  rrtmgp@aer.com
!
! Copyright 2015-2018,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
! -------------------------------------------------------------------------------------------------
!
!  Contains a single routine to compute direct and diffuse fluxes of solar radiation given
!    atmospheric optical properties on a spectral grid
!    information about vertical ordering
!    boundary conditions
!      solar zenith angle, spectrally-resolved incident colimated flux, surface albedos for direct and diffuse radiation
!    optionally, a boundary condition for incident diffuse radiation
!
! It is the user's responsibility to ensure that boundary conditions (incident fluxes, surface albedos) are on the same
!   spectral grid as the optical properties.
!
! Final output is via user-extensible ty_fluxes which must reduce the detailed spectral fluxes to
!   whatever summary the user needs.
!
! The routine does error checking and choses which lower-level kernel to invoke based on
!   what kinds of optical properties are supplied
!
! -------------------------------------------------------------------------------------------------
module mo_rte_sw
  use mo_rte_kind,      only: wp, wl
  use mo_optical_props, only: ty_optical_props, &
                              ty_optical_props_arry, ty_optical_props_1scl, ty_optical_props_2str, ty_optical_props_nstr
  use mo_fluxes,        only: ty_fluxes
  use mo_rte_solver_kernels, &
                        only: apply_BC, sw_solver_noscat, sw_solver_2stream
  implicit none
  private

  public :: rte_sw

contains
  ! --------------------------------------------------
  function rte_sw(atmos, top_at_1,                 &
                  mu0, inc_flux,                   &
                  sfc_alb_dir, sfc_alb_dif,        &
                  fluxes, inc_flux_dif) result(error_msg)
    class(ty_optical_props_arry), intent(in   ) :: atmos           ! Optical properties provided as arrays
    logical,                      intent(in   ) :: top_at_1        ! Is the top of the domain at index 1?
                                                                   ! (if not, ordering is bottom-to-top)
    real(wp), dimension(:),       intent(in   ) :: mu0             ! cosine of solar zenith angle (ncol)
    real(wp), dimension(:,:),     intent(in   ) :: inc_flux,    &  ! incident flux at top of domain [W/m2] (ncol, ngpt)
                                                   sfc_alb_dir, &  ! surface albedo for direct and
                                                   sfc_alb_dif     ! diffuse radiation (nband, ncol)
    class(ty_fluxes),             intent(inout) :: fluxes          ! Class describing output calculations
    real(wp), dimension(:,:), optional, &
                                  intent(in   ) :: inc_flux_dif    ! incident diffuse flux at top of domain [W/m2] (ncol, ngpt)
    character(len=128)                          :: error_msg       ! If empty, calculation was successful
    ! --------------------------------
    !
    ! Local variables
    !
    integer :: ncol, nlay, ngpt, nband
    integer :: icol

    real(wp), dimension(:,:,:), allocatable :: gpt_flux_up, gpt_flux_dn, gpt_flux_dir
    real(wp), dimension(:,:),   allocatable :: sfc_alb_dir_gpt, sfc_alb_dif_gpt
    ! ------------------------------------------------------------------------------------
    ncol  = atmos%get_ncol()
    nlay  = atmos%get_nlay()
    ngpt  = atmos%get_ngpt()
    nband = atmos%get_nband()
    error_msg = ""

    allocate(gpt_flux_up (ncol, nlay+1, ngpt), gpt_flux_dn(ncol, nlay+1, ngpt), gpt_flux_dir(ncol, nlay+1, ngpt))
    allocate(sfc_alb_dir_gpt(ncol, ngpt), sfc_alb_dif_gpt(ncol, ngpt))

    ! ------------------------------------------------------------------------------------
    !
    ! Error checking -- consistency of sizes and validity of values
    !
    ! --------------------------------
    if(.not. fluxes%are_desired()) then
      error_msg = "rte_sw: no space allocated for fluxes"
      return
    end if

    !
    ! Sizes and values of input arrays
    !
    if(     size(mu0, 1)                                /=  ncol        ) &
      error_msg = "rte_sw: mu0 inconsistently sized"
    if(any(mu0 <= 0._wp)) &
      error_msg = "rte_sw: one or more mu0 <= 0"

    if(any([size(inc_flux, 1),    size(inc_flux, 2)]    /= [ncol, ngpt])) &
      error_msg = "rte_sw: inc_flux inconsistently sized"
    if(any(inc_flux <  0._wp)) &
      error_msg = "rte_sw: one or more inc_flux < 0"
    if(present(inc_flux_dif)) then
      if(any([size(inc_flux_dif, 1),    size(inc_flux_dif, 2)]    /= [ncol, ngpt])) &
        error_msg = "rte_sw: inc_flux_dif inconsistently sized"
      if(any(inc_flux_dif <  0._wp)) &
        error_msg = "rte_sw: one or more inc_flux_dif < 0"
    end if

    if(any([size(sfc_alb_dir, 1), size(sfc_alb_dir, 2)] /= [nband, ncol])) &
      error_msg = "rte_sw: sfc_alb_dir inconsistently sized"
    if(any(sfc_alb_dir < 0._wp .or. sfc_alb_dir > 1._wp)) &
      error_msg = "rte_sw: sfc_alb_dir out of bounds [0,1]"
    if(any([size(sfc_alb_dif, 1), size(sfc_alb_dif, 2)] /= [nband, ncol])) &
      error_msg = "rte_sw: sfc_alb_dif inconsistently sized"
    if(any(sfc_alb_dif < 0._wp .or. sfc_alb_dif > 1._wp)) &
      error_msg = "rte_sw: sfc_alb_dif out of bounds [0,1]"

    if(len_trim(error_msg) > 0) return

    !
    ! Ensure values of tau, ssa, and g are reasonable
    !
    error_msg =  atmos%validate()
    if(len_trim(error_msg) > 0) then
      if(len_trim(atmos%get_name()) > 0) &
        error_msg = trim(atmos%get_name()) // ': ' // trim(error_msg)
      return
    end if

    ! ------------------------------------------------------------------------------------
    ! Lower boundary condition -- expand surface albedos by band to gpoints
    !   and switch dimension ordering
    do icol = 1, ncol
      sfc_alb_dir_gpt(icol, 1:ngpt) = atmos%expand(sfc_alb_dir(:,icol))
    end do
    do icol = 1, ncol
      sfc_alb_dif_gpt(icol, 1:ngpt) = atmos%expand(sfc_alb_dif(:,icol))
    end do

    ! ------------------------------------------------------------------------------------
    !
    ! Compute the radiative transfer...
    !
    !
    ! Apply boundary conditions
    !   On input flux_dn is the diffuse component; the last action in each solver is to add
    !   direct and diffuse to represent the total, consistent with the LW
    !
    call apply_BC(ncol, nlay, ngpt, logical(top_at_1, wl),   inc_flux, mu0, gpt_flux_dir)
    if(present(inc_flux_dif)) then
      call apply_BC(ncol, nlay, ngpt, logical(top_at_1, wl), inc_flux_dif,  gpt_flux_dn )
    else
      call apply_BC(ncol, nlay, ngpt, logical(top_at_1, wl),                gpt_flux_dn )
    end if

    select type (atmos)
      class is (ty_optical_props_1scl)
        !
        ! Direct beam only
        !
        call sw_solver_noscat(ncol, nlay, ngpt, logical(top_at_1, wl), &
                              atmos%tau, mu0,                          &
                              gpt_flux_dir)
        !
        ! No diffuse flux
        !
        gpt_flux_up = 0._wp
        gpt_flux_dn = 0._wp
      class is (ty_optical_props_2str)
        !
        ! two-stream calculation with scattering
        !
        call sw_solver_2stream(ncol, nlay, ngpt, logical(top_at_1, wl), &
                               atmos%tau, atmos%ssa, atmos%g, mu0,      &
                               sfc_alb_dir_gpt, sfc_alb_dif_gpt,        &
                               gpt_flux_up, gpt_flux_dn, gpt_flux_dir)
      class is (ty_optical_props_nstr)
        !
        ! n-stream calculation
        !
        ! not yet implemented so fail
        !
        error_msg = 'sw_solver(...ty_optical_props_nstr...) not yet implemented'
    end select
    if (error_msg /= '') return
    !
    ! ...and reduce spectral fluxes to desired output quantities
    !
    error_msg = fluxes%reduce(gpt_flux_up, gpt_flux_dn, atmos, top_at_1, gpt_flux_dir)
  end function rte_sw
end module mo_rte_sw
