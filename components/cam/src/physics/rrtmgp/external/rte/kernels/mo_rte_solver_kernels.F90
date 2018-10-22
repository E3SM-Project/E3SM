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
! Numeric calculations for radiative transfer solvers.
!   Emission/absorption (no-scattering) calculations
!     solver for multi-angle Gaussian quadrature
!     solver for a single angle, calling
!       source function computation (linear-in-tau)
!       transport
!   Extinction-only calculation (direct solar beam)
!   Two-stream calculations
!     solvers for LW and SW with different boundary conditions and source functions
!       source function calculation for LW, SW
!       two-stream calculations for LW, SW (using different assumtions about phase function)
!       transport (adding)
!   Application of boundary conditions
!
! -------------------------------------------------------------------------------------------------
module mo_rte_solver_kernels
  use,  intrinsic :: iso_c_binding
  use mo_rte_kind, only: wp, wl
  implicit none
  private

  interface apply_BC
    module procedure apply_BC_gpt, apply_BC_factor, apply_BC_0
  end interface apply_BC

  public :: apply_BC, &
            lw_solver_noscat, lw_solver_noscat_GaussQuad, lw_solver_2stream, &
            sw_solver_noscat,                             sw_solver_2stream

  ! These routines don't really need to be visible but making them so is useful for testing.
  public :: lw_source_noscat, lw_combine_sources, &
            lw_source_2str, sw_source_2str, &
            lw_two_stream, sw_two_stream, &
            adding

  real(wp), parameter :: pi = acos(-1._wp)
contains
  ! -------------------------------------------------------------------------------------------------
  !
  ! Top-level longwave kernels
  !
  ! -------------------------------------------------------------------------------------------------
  !
  ! LW fluxes, no scattering, mu (cosine of integration angle) specified by column
  !   Does radiation calculation at user-supplied angles; converts radiances to flux
  !   using user-supplied weights
  !
  ! ---------------------------------------------------------------
  subroutine lw_solver_noscat(ncol, nlay, ngpt, top_at_1, D, weight,                             &
                              tau, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                              radn_up, radn_dn) bind(C, name="lw_solver_noscat")
    integer,                    intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                intent( in) :: top_at_1
    real(wp), dimension(ncol,       ngpt), intent( in) :: D            ! secant of propagation angle  []
    real(wp),                              intent( in) :: weight       ! quadrature weight
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: lay_source   ! Planck source at layer average temperature [W/m2]
    ! Planck source at layer edge for radiation in increasing/decreasing ilay direction
    ! lev_source_dec applies the mapping in layer i to the Planck function at layer i
    ! lev_source_inc applies the mapping in layer i to the Planck function at layer i+1
    real(wp), dimension(ncol,nlay,  ngpt), target, &
                                           intent( in) :: lev_source_inc, lev_source_dec
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_emis         ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_src          ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(out) :: radn_up, radn_dn ! Radiances [W/m2-str]
                                                                           ! Top level must contain incident flux boundary condition

    ! Local variables, no g-point dependency
    real(wp), dimension(ncol,nlay) :: tau_loc, &  ! path length (tau/mu)
                                        trans       ! transmissivity  = exp(-tau)
    real(wp), dimension(ncol,nlay) :: source_dn, source_up
    real(wp), dimension(ncol     ) :: source_sfc, sfc_albedo

    real(wp), dimension(:,:,:), pointer :: lev_source_up, lev_source_dn ! Mapping increasing/decreasing indicies to up/down

    real(wp), parameter :: pi = acos(-1._wp)
    integer             :: ilev, igpt, top_level
    ! ------------------------------------
    ! Which way is up?
    ! Level Planck sources for upward and downward radiation
    ! When top_at_1, lev_source_up => lev_source_dec
    !                lev_source_dn => lev_source_inc, and vice-versa
    if(top_at_1) then
      top_level = 1
      lev_source_up => lev_source_dec
      lev_source_dn => lev_source_inc
    else
      top_level = nlay+1
      lev_source_up => lev_source_inc
      lev_source_dn => lev_source_dec
    end if

    do igpt = 1, ngpt
      !
      ! Transport is for intensity
      !   convert flux at top of domain to intensity assuming azimuthal isotropy
      !
      radn_dn(:,top_level,igpt) = radn_dn(:,top_level,igpt)/(2._wp * pi * weight)

      !
      ! Optical path and transmission, used in source function and transport calculations
      !
      do ilev = 1, nlay
        tau_loc(:,ilev) = tau(:,ilev,igpt)*D(:,igpt)
        trans  (:,ilev) = exp(-tau_loc(:,ilev))
      end do
      !
      ! Source function for diffuse radiation
      !
      call lw_source_noscat(ncol, nlay, &
                            lay_source(:,:,igpt), lev_source_up(:,:,igpt), lev_source_dn(:,:,igpt), &
                            tau_loc, trans, source_dn, source_up)
      !
      ! Surface albedo, surface source function
      !
      sfc_albedo(:) = 1._wp - sfc_emis(:,igpt)
      source_sfc(:) = sfc_emis(:,igpt) * sfc_src(:,igpt)
      !
      ! Transport
      !
      call lw_transport_noscat(ncol, nlay, top_at_1,  &
                               tau_loc, trans, sfc_albedo, source_dn, source_up, source_sfc, &
                               radn_up(:,:,igpt), radn_dn(:,:,igpt))
      !
      ! Convert intensity to flux assuming azimuthal isotropy and quadrature weight
      !
      radn_dn(:,:,igpt) = 2._wp * pi * weight * radn_dn(:,:,igpt)
      radn_up(:,:,igpt) = 2._wp * pi * weight * radn_up(:,:,igpt)
    end do  ! g point loop

  end subroutine lw_solver_noscat
  ! -------------------------------------------------------------------------------------------------
  !
  ! LW transport, no scattering, multi-angle quadrature
  !   Users provide a set of weights and quadrature angles
  !   Routine sums over single-angle solutions for each sets of angles/weights
  !
  ! ---------------------------------------------------------------
  subroutine lw_solver_noscat_GaussQuad(ncol, nlay, ngpt, top_at_1, nmus, Ds, weights, &
                                   tau, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, flux_up, flux_dn) &
                                   bind(C, name="lw_solver_noscat_GaussQuad")
    integer,                    intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                intent( in) :: top_at_1
    integer,                    intent( in) :: nmus          ! number of quadrature angles
    real(wp), dimension(nmus),  intent( in) :: Ds, weights  ! quadrature secants, weights
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: lay_source   ! Planck source at layer average temperature [W/m2]
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: lev_source_inc
                                        ! Planck source at layer edge for radiation in increasing ilay direction [W/m2]
                                        ! Includes spectral weighting that accounts for state-dependent frequency to g-space mapping
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: lev_source_dec
                                               ! Planck source at layer edge for radiation in decreasing ilay direction [W/m2]
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_emis         ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_src          ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(out) :: flux_up, flux_dn ! Radiances [W/m2-str]
                                                                           ! Top level must contain incident flux boundary condition
    ! Local variables
    real(wp), dimension(ncol,nlay+1,ngpt) :: radn_dn, radn_up ! Fluxes per quad angle
    real(wp), dimension(ncol,       ngpt) :: Ds_ncol

    integer :: imu, top_level
    ! ------------------------------------
    !
    ! For the first angle output arrays store total flux
    !
    Ds_ncol(:,:) = Ds(1)
    call lw_solver_noscat(ncol, nlay, ngpt, &
                          top_at_1, Ds_ncol, weights(1), tau, &
                          lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                          flux_up, flux_dn)
    !
    ! For more than one angle use local arrays
    !
    top_level = MERGE(1, nlay+1, top_at_1)
    radn_dn(:,top_level,:) = flux_dn(:, top_level, :) ! Flux boundary condition
    do imu = 2, nmus
      Ds_ncol(:,:) = Ds(imu)
      call lw_solver_noscat(ncol, nlay, ngpt, &
                            top_at_1, Ds_ncol, weights(imu), tau, &
                            lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                            radn_up, radn_dn)
      flux_up(:,:,:) = flux_up(:,:,:) + radn_up(:,:,:)
      flux_dn(:,:,:) = flux_dn(:,:,:) + radn_dn(:,:,:)
    end do
  end subroutine lw_solver_noscat_GaussQuad
  ! -------------------------------------------------------------------------------------------------
  !
  ! Longwave two-stream calculation:
  !   combine RRTMGP-specific sources at levels
  !   compute layer reflectance, transmittance
  !   compute total source function at levels using linear-in-tau
  !   transport
  !
  ! -------------------------------------------------------------------------------------------------
   subroutine lw_solver_2stream (ncol, nlay, ngpt, top_at_1, &
                                 tau, ssa, g,                &
                                 lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                                 flux_up, flux_dn) bind(C, name="lw_solver_2stream")
    integer,                    intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                intent( in) :: top_at_1
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: tau, &  ! Optical thickness,
                                                          ssa, &  ! single-scattering albedo,
                                                          g       ! asymmetry parameter []
    real(wp), dimension(ncol,nlay,ngpt), intent( in) :: lay_source   ! Planck source at layer average temperature [W/m2]
    real(wp), dimension(ncol,nlay,ngpt), target, &
                                         intent( in) :: lev_source_inc, lev_source_dec
                                        ! Planck source at layer edge for radiation in increasing/decreasing ilay direction [W/m2]
                                        ! Includes spectral weighting that accounts for state-dependent frequency to g-space mapping
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_emis         ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_src          ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), &
                                           intent(out) :: flux_up, flux_dn  ! Fluxes [W/m2]
                                                                            ! Top level (= merge(1, nlay+1, top_at_1)
                                                                            ! must contain incident flux boundary condition
    ! ----------------------------------------------------------------------
    integer :: igpt
    real(wp), dimension(ncol,nlay  ) :: Rdif, Tdif, gamma1, gamma2
    real(wp), dimension(ncol       ) :: sfc_albedo
    real(wp), dimension(ncol,nlay+1) :: lev_source
    real(wp), dimension(ncol,nlay  ) :: source_dn, source_up
    real(wp), dimension(ncol       ) :: source_sfc
    ! ------------------------------------

    do igpt = 1, ngpt
      !
      ! RRTMGP provides source functions at each level using the spectral mapping
      !   of each adjacent layer. Combine these for two-stream calculations
      !
      call lw_combine_sources(ncol, nlay, top_at_1, &
                              lev_source_inc(:,:,igpt), lev_source_dec(:,:,igpt), &
                              lev_source)
      !
      ! Cell properties: reflection, transmission for diffuse radiation
      !   Coupling coefficients needed for source function
      !
      call lw_two_stream(ncol, nlay,                                 &
                         tau (:,:,igpt), ssa(:,:,igpt), g(:,:,igpt), &
                         gamma1, gamma2, Rdif, Tdif)
      !
      ! Source function for diffuse radiation
      !
      call lw_source_2str(ncol, nlay, top_at_1, &
                          sfc_emis(:,igpt), sfc_src(:,igpt), &
                          lay_source(:,:,igpt), lev_source, &
                          gamma1, gamma2, Rdif, Tdif, tau(:,:,igpt), &
                          source_dn, source_up, source_sfc)
      !
      ! Transport
      !
      sfc_albedo(1:ncol) = 1._wp - sfc_emis(:,igpt)
      call adding(ncol, nlay, top_at_1,              &
                  sfc_albedo,                        &
                  Rdif, Tdif,                        &
                  source_dn, source_up, source_sfc,  &
                  flux_up(:,:,igpt), flux_dn(:,:,igpt))
    end do

  end subroutine lw_solver_2stream
  ! -------------------------------------------------------------------------------------------------
  !
  !   Top-level shortwave kernels
  !
  ! -------------------------------------------------------------------------------------------------
  !
  !   Extinction-only i.e. solar direct beam
  !
  ! -------------------------------------------------------------------------------------------------
  pure subroutine sw_solver_noscat(ncol, nlay, ngpt, &
                              top_at_1, tau, mu0, flux_dir) bind(C, name="sw_solver_noscat")
    integer,                    intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                intent( in) :: top_at_1
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol            ), intent( in) :: mu0          ! cosine of solar zenith angle
    real(wp), dimension(ncol,nlay+1,ngpt), intent(out) :: flux_dir     ! Direct-beam flux, spectral [W/m2]
                                                                       ! Top level must contain incident flux boundary condition
    integer :: icol, ilev, igpt
    real(wp) :: mu0_inv(ncol)

    ! ------------------------------------
    mu0_inv(1:ncol) = 1._wp/mu0(1:ncol)
    ! Indexing into arrays for upward and downward propagation depends on the vertical
    !   orientation of the arrays (whether the domain top is at the first or last index)
    ! We write the loops out explicitly so compilers will have no trouble optimizing them.

    ! Downward propagation
    if(top_at_1) then
      ! For the flux at this level, what was the previous level, and which layer has the
      !   radiation just passed through?
      ! layer index = level index - 1
      ! previous level is up (-1)
      do igpt = 1, ngpt
        do ilev = 2, nlay+1
          flux_dir(:,ilev,igpt) = flux_dir(:,ilev-1,igpt) * exp(-tau(:,ilev,igpt)*mu0_inv(:))
        end do
      end do
    else
      ! layer index = level index
      ! previous level is up (+1)
      do igpt = 1, ngpt
        do ilev = nlay, 1, -1
          flux_dir(:,ilev,igpt) = flux_dir(:,ilev+1,igpt) * exp(-tau(:,ilev,igpt)*mu0_inv(:))
        end do
      end do
    end if
  end subroutine sw_solver_noscat
  ! -------------------------------------------------------------------------------------------------
  !
  ! Shortwave two-stream calculation:
  !   compute layer reflectance, transmittance
  !   compute solar source function for diffuse radiation
  !   transport
  !
  ! -------------------------------------------------------------------------------------------------
   subroutine sw_solver_2stream (ncol, nlay, ngpt, top_at_1, &
                                 tau, ssa, g, mu0,           &
                                 sfc_alb_dir, sfc_alb_dif,   &
                                 flux_up, flux_dn, flux_dir) bind(C, name="sw_solver_2stream")
    integer,                    intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                intent( in) :: top_at_1
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: tau, &  ! Optical thickness,
                                                          ssa, &  ! single-scattering albedo,
                                                          g       ! asymmetry parameter []
    real(wp), dimension(ncol            ), intent( in) :: mu0     ! cosine of solar zenith angle
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_alb_dir, sfc_alb_dif
                                                                  ! Spectral albedo of surface to direct and diffuse radiation
    real(wp), dimension(ncol,nlay+1,ngpt), &
                                           intent(out) :: flux_up, flux_dn, &  ! Fluxes [W/m2]
                                                          flux_dir             ! Downward direct
                                                                               ! Top level (= merge(1, nlay+1, top_at_1)
                                                                               ! must contain incident flux boundary condition
    integer :: igpt
    real(wp), dimension(ncol,nlay) :: Rdif, Tdif, Rdir, Tdir, Tnoscat
    real(wp), dimension(ncol,nlay) :: source_up, source_dn
    real(wp), dimension(ncol     ) :: source_srf
    ! ------------------------------------
    do igpt = 1, ngpt
      !
      ! Cell properties: transmittance and reflectance for direct and diffuse radiation
      !
      call sw_two_stream(ncol, nlay, mu0,                                &
                         tau (:,:,igpt), ssa (:,:,igpt), g   (:,:,igpt), &
                         Rdif, Tdif, Rdir, Tdir, Tnoscat)
      !
      ! Direct-beam and source for diffuse radiation
      !
      call sw_source_2str(ncol, nlay, top_at_1, Rdir, Tdir, Tnoscat, sfc_alb_dir(:,igpt),&
                          source_up, source_dn, source_srf, flux_dir(:,:,igpt))
      !
      ! Transport
      !
      call adding(ncol, nlay, top_at_1,            &
                     sfc_alb_dif(:,igpt), Rdif, Tdif, &
                     source_dn, source_up, source_srf, flux_up(:,:,igpt), flux_dn(:,:,igpt))
      !
      ! adding computes only diffuse flux; flux_dn is total
      !
      flux_dn(:,:,igpt) = flux_dn(:,:,igpt) + flux_dir(:,:,igpt)
    end do

  end subroutine sw_solver_2stream
  ! -------------------------------------------------------------------------------------------------
  !
  !   Lower-level longwave kernels
  !
  ! -------------------------------------------------------------------------------------------------
  !
  ! Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
  ! See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13
  !
  ! ---------------------------------------------------------------
  subroutine lw_source_noscat(ncol, nlay, lay_source, lev_source_up, lev_source_dn, tau, trans, &
                              source_dn, source_up) bind(C, name="lw_source_noscat")
    integer,                         intent(in) :: ncol, nlay
    real(wp), dimension(ncol, nlay), intent(in) :: lay_source, & ! Planck source at layer center
                                                   lev_source_up, & ! Planck source at levels (layer edges),
                                                   lev_source_dn, & !   increasing/decreasing layer index
                                                   tau,        & ! Optical path (tau/mu)
                                                   trans         ! Transmissivity (exp(-tau))
    real(wp), dimension(ncol, nlay), intent(out):: source_dn, source_up
                                                                   ! Source function at layer edges
                                                                   ! Down at the bottom of the layer, up at the top
    ! --------------------------------
    integer             :: icol, ilay
    real(wp)            :: fact
    real(wp), parameter :: tau_thresh = sqrt(epsilon(tau))
    ! ---------------------------------------------------------------
    do ilay = 1, nlay
      do icol = 1, ncol
      !
      ! Weighting factor. Use 2nd order series expansion when rounding error (~tau^2)
      !   is of order epsilon (smallest difference from 1. in working precision)
      !   Thanks to Peter Blossey
      !
      fact = merge((1._wp - trans(icol,ilay))/tau(icol,ilay) - trans(icol,ilay), &
                   tau(icol, ilay) * ( 0.5_wp - 1._wp/3._wp*tau(icol, ilay)   ), &
                   tau(icol, ilay) > tau_thresh)
      !
      ! Equation below is developed in Clough et al., 1992, doi:10.1029/92JD01419, Eq 13
      !
      source_dn(icol,ilay) = (1._wp - trans(icol,ilay)) * lev_source_dn(icol,ilay) + &
                              2._wp * fact * (lay_source(icol,ilay) - lev_source_dn(icol,ilay))
      source_up(icol,ilay) = (1._wp - trans(icol,ilay)) * lev_source_up(icol,ilay  ) + &
                              2._wp * fact * (lay_source(icol,ilay) - lev_source_up(icol,ilay))
      end do
    end do
  end subroutine lw_source_noscat
  ! -------------------------------------------------------------------------------------------------
  !
  ! Longwave no-scattering transport
  !
  ! -------------------------------------------------------------------------------------------------
  subroutine lw_transport_noscat(ncol, nlay, top_at_1, &
                                 tau, trans, sfc_albedo, source_dn, source_up, source_sfc, &
                                 radn_up, radn_dn) bind(C, name="lw_transport_noscat")
    integer,                          intent( in) :: ncol, nlay ! Number of columns, layers, g-points
    logical(wl),                      intent( in) :: top_at_1   !
    real(wp), dimension(ncol,nlay  ), intent( in) :: tau, &     ! Absorption optical thickness, pre-divided by mu []
                                                     trans      ! transmissivity = exp(-tau)
    real(wp), dimension(ncol       ), intent( in) :: sfc_albedo ! Surface albedo
    real(wp), dimension(ncol,nlay  ), intent( in) :: source_dn, &
                                                     source_up  ! Diffuse radiation emitted by the layer
    real(wp), dimension(ncol       ), intent( in) :: source_sfc ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1), intent(out) :: radn_up, radn_dn ! Radiances [W/m2-str]
                                                                      ! Top level must contain incident flux boundary condition
    ! Local variables
    integer :: ilev
    ! ---------------------------------------------------
    if(top_at_1) then
      !
      ! Top of domain is index 1
      !
      ! Downward propagation
      do ilev = 2, nlay+1
        radn_dn(:,ilev) = trans(:,ilev-1)*radn_dn(:,ilev-1) + source_dn(:,ilev-1)
      end do

      ! Surface reflection and emission
      radn_up(:,nlay+1) = radn_dn(:,nlay+1)*sfc_albedo(:) + source_sfc(:)

      ! Upward propagation
      do ilev = nlay, 1, -1
        radn_up(:,ilev) = trans(:,ilev  )*radn_up(:,ilev+1) + source_up(:,ilev)
      end do
    else
      !
      ! Top of domain is index nlay+1
      !
      ! Downward propagation
      do ilev = nlay, 1, -1
        radn_dn(:,ilev) = trans(:,ilev  )*radn_dn(:,ilev+1) + source_dn(:,ilev)
      end do

      ! Surface reflection and emission
      radn_up(:,     1) = radn_dn(:,     1)*sfc_albedo(:) + source_sfc(:)

      ! Upward propagation
      do ilev = 2, nlay+1
        radn_up(:,ilev) = trans(:,ilev-1) * radn_up(:,ilev-1) +  source_up(:,ilev-1)
      end do
    end if
  end subroutine lw_transport_noscat
  ! -------------------------------------------------------------------------------------------------
  !
  ! Longwave two-stream solutions to diffuse reflectance and transmittance for a layer
  !    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
  !
  ! Equations are developed in Meador and Weaver, 1980,
  !    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
  !
  ! -------------------------------------------------------------------------------------------------
  pure subroutine lw_two_stream(ncol, nlay, tau, w0, g, &
                                gamma1, gamma2, Rdif, Tdif) bind(C, name="lw_two_stream")
    integer,                        intent(in)  :: ncol, nlay
    real(wp), dimension(ncol,nlay), intent(in)  :: tau, w0, g
    real(wp), dimension(ncol,nlay), intent(out) :: gamma1, gamma2, Rdif, Tdif

    ! -----------------------
    integer  :: i, j

    ! Variables used in Meador and Weaver
    real(wp) :: k(ncol)

    ! Ancillary variables
    real(wp) :: RT_term(ncol)
    real(wp) :: exp_minusktau(ncol), exp_minus2ktau(ncol)

    real(wp), parameter :: LW_diff_sec = 1.66  ! 1./cos(diffusivity angle)
    ! ---------------------------------
    do j = 1, nlay
      do i = 1, ncol
        !
        ! Coefficients differ from SW implementation because the phase function is more isotropic
        !   Here we follow Fu et al. 1997, doi:10.1175/1520-0469(1997)054<2799:MSPITI>2.0.CO;2
        !   and use a diffusivity sec of 1.66
        !
        gamma1(i,j)= LW_diff_sec * (1._wp - 0.5_wp * w0(i,j) * (1._wp + g(i,j))) ! Fu et al. Eq 2.9
        gamma2(i,j)= LW_diff_sec *          0.5_wp * w0(i,j) * (1._wp - g(i,j))  ! Fu et al. Eq 2.10
      end do

      ! Written to encourage vectorization of exponential, square root
      ! Eq 18;  k = SQRT(gamma1**2 - gamma2**2), limited below to avoid div by 0.
      !   k = 0 for isotropic, conservative scattering; this lower limit on k
      !   gives relative error with respect to conservative solution
      !   of < 0.1% in Rdif down to tau = 10^-9
      k(1:ncol) = sqrt(max((gamma1(1:ncol,j) - gamma2(1:ncol,j)) * &
                           (gamma1(1:ncol,j) + gamma2(1:ncol,j)),  &
                           1.e-12_wp))
      exp_minusktau(1:ncol) = exp(-tau(1:ncol,j)*k(1:ncol))

      !
      ! Diffuse reflection and transmission
      !
      do i = 1, ncol
        exp_minus2ktau(i) = exp_minusktau(i) * exp_minusktau(i)

        ! Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
        RT_term(i) = 1._wp / (k     (i  ) * (1._wp + exp_minus2ktau(i))  + &
                              gamma1(i,j) * (1._wp - exp_minus2ktau(i)) )

        ! Equation 25
        Rdif(i,j) = RT_term(i) * gamma2(i,j) * (1._wp - exp_minus2ktau(i))

        ! Equation 26
        Tdif(i,j) = RT_term(i) * 2._wp * k(i) * exp_minusktau(i)
      end do

    end do
  end subroutine lw_two_stream
  ! -------------------------------------------------------------------------------------------------
  !
  ! Source function combination
  ! RRTMGP provides two source functions at each level
  !   using the spectral mapping from each of the adjascent layers.
  !   Need to combine these for use in two-stream calculation.
  !
  ! -------------------------------------------------------------------------------------------------
  subroutine lw_combine_sources(ncol, nlay, top_at_1, &
                                lev_src_inc, lev_src_dec, lev_source) bind(C, name="lw_combine_sources")
    integer,                           intent(in ) :: ncol, nlay
    logical(wl),                       intent(in ) :: top_at_1
    real(wp), dimension(ncol, nlay  ), intent(in ) :: lev_src_inc, lev_src_dec
    real(wp), dimension(ncol, nlay+1), intent(out) :: lev_source

    integer :: icol, ilay
    ! ---------------------------------------------------------------
    ilay = 1
    do icol = 1,ncol
      lev_source(icol, ilay) =        lev_src_dec(icol, ilay)
    end do
    do ilay = 2, nlay
      do icol = 1,ncol
        lev_source(icol, ilay) = sqrt(lev_src_dec(icol, ilay) * &
                                      lev_src_inc(icol, ilay-1))
      end do
    end do
    ilay = nlay+1
    do icol = 1,ncol
      lev_source(icol, ilay) =        lev_src_inc(icol, ilay-1)
    end do

  end subroutine lw_combine_sources
  ! ---------------------------------------------------------------
  !
  ! Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
  !   This version straight from ECRAD
  !   Source is provided as W/m2-str; factor of pi converts to flux units
  !
  ! ---------------------------------------------------------------
  subroutine lw_source_2str(ncol, nlay, top_at_1,   &
                            sfc_emis, sfc_src,      &
                            lay_source, lev_source, &
                            gamma1, gamma2, rdif, tdif, tau, source_dn, source_up, source_sfc) &
                            bind (C, name="lw_source_2str")
    integer,                         intent(in) :: ncol, nlay
    logical(wl),                     intent(in) :: top_at_1
    real(wp), dimension(ncol      ), intent(in) :: sfc_emis, sfc_src
    real(wp), dimension(ncol, nlay), intent(in) :: lay_source,    & ! Planck source at layer center
                                                   tau,           & ! Optical depth (tau)
                                                   gamma1, gamma2,& ! Coupling coefficients
                                                   rdif, tdif       ! Layer reflectance and transmittance
    real(wp), dimension(ncol, nlay+1), target, &
                                     intent(in)  :: lev_source       ! Planck source at layer edges
    real(wp), dimension(ncol, nlay), intent(out) :: source_dn, source_up
    real(wp), dimension(ncol      ), intent(out) :: source_sfc      ! Source function for upward radation at surface

    integer             :: icol, ilay
    real(wp)            :: Z, Zup_top, Zup_bottom, Zdn_top, Zdn_bottom
    real(wp), dimension(:), pointer :: lev_source_bot, lev_source_top
    ! ---------------------------------------------------------------
    do ilay = 1, nlay
      if(top_at_1) then
        lev_source_top => lev_source(:,ilay)
        lev_source_bot => lev_source(:,ilay+1)
      else
        lev_source_top => lev_source(:,ilay+1)
        lev_source_bot => lev_source(:,ilay)
      end if
      do icol = 1, ncol
        if (tau(icol,ilay) > 1.0e-8_wp) then
          !
          ! Toon et al. (JGR 1989) Eqs 26-27
          !
          Z = (lev_source_bot(icol)-lev_source_top(icol)) / (tau(icol,ilay)*(gamma1(icol,ilay)+gamma2(icol,ilay)))
          Zup_top        =  Z + lev_source_top(icol)
          Zup_bottom     =  Z + lev_source_bot(icol)
          Zdn_top        = -Z + lev_source_top(icol)
          Zdn_bottom     = -Z + lev_source_bot(icol)
          source_up(icol,ilay) = pi * (Zup_top    - rdif(icol,ilay) * Zdn_top    - tdif(icol,ilay) * Zup_bottom)
          source_dn(icol,ilay) = pi * (Zdn_bottom - rdif(icol,ilay) * Zup_bottom - tdif(icol,ilay) * Zdn_top)
        else
          source_up(icol,ilay) = 0._wp
          source_dn(icol,ilay) = 0._wp
        end if
      end do
    end do
    do icol = 1, ncol
      source_sfc(icol) = pi * sfc_emis(icol) * sfc_src(icol)
    end do
  end subroutine lw_source_2str
  ! -------------------------------------------------------------------------------------------------
  !
  !   Lower-level shortwave kernels
  !
  ! -------------------------------------------------------------------------------------------------
  !
  ! Two-stream solutions to direct and diffuse reflectance and transmittance for a layer
  !    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.
  !
  ! Equations are developed in Meador and Weaver, 1980,
  !    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
  !
  ! -------------------------------------------------------------------------------------------------
  pure subroutine sw_two_stream(ncol, nlay, mu0, tau, w0, g, &
                                Rdif, Tdif, Rdir, Tdir, Tnoscat) bind (C, name="sw_two_stream")
    integer,                        intent(in)  :: ncol, nlay
    real(wp), dimension(ncol),      intent(in)  :: mu0
    real(wp), dimension(ncol,nlay), intent(in)  :: tau, w0, g
    real(wp), dimension(ncol,nlay), intent(out) :: Rdif, Tdif, Rdir, Tdir, Tnoscat

    ! -----------------------
    integer  :: i, j

    ! Variables used in Meador and Weaver
    real(wp) :: gamma1(ncol), gamma2(ncol), gamma3(ncol), gamma4(ncol)
    real(wp) :: alpha1(ncol), alpha2(ncol), k(ncol)

    ! Ancillary variables
    real(wp) :: RT_term(ncol)
    real(wp) :: exp_minusktau(ncol), exp_minus2ktau(ncol)
    real(WP) :: k_mu, k_gamma3, k_gamma4
    real(wp) :: mu0_inv(ncol)
    ! ---------------------------------
    mu0_inv(1:ncol) = 1._wp/mu0(1:ncol)
    do j = 1, nlay
      do i = 1, ncol
        ! Zdunkowski Practical Improved Flux Method "PIFM"
        !  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
        !
        gamma1(i)= (8._wp - w0(i,j) * (5._wp + 3._wp * g(i,j))) * .25_wp
        gamma2(i)=  3._wp *(w0(i,j) * (1._wp -         g(i,j))) * .25_wp
        gamma3(i)= (2._wp - 3._wp * mu0(i) *           g(i,j) ) * .25_wp
        gamma4(i)=  1._wp - gamma3(i)

        alpha1(i) = gamma1(i) * gamma4(i) + gamma2(i) * gamma3(i)           ! Eq. 16
        alpha2(i) = gamma1(i) * gamma3(i) + gamma2(i) * gamma4(i)           ! Eq. 17
      end do

      ! Written to encourage vectorization of exponential, square root
      ! Eq 18;  k = SQRT(gamma1**2 - gamma2**2), limited below to avoid div by 0.
      !   k = 0 for isotropic, conservative scattering; this lower limit on k
      !   gives relative error with respect to conservative solution
      !   of < 0.1% in Rdif down to tau = 10^-9
      k(1:ncol) = sqrt(max((gamma1(1:ncol) - gamma2(1:ncol)) * &
                           (gamma1(1:ncol) + gamma2(1:ncol)),  &
                           1.e-12_wp))
      exp_minusktau(1:ncol) = exp(-tau(1:ncol,j)*k(1:ncol))

      !
      ! Diffuse reflection and transmission
      !
      do i = 1, ncol
        exp_minus2ktau(i) = exp_minusktau(i) * exp_minusktau(i)

        ! Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
        RT_term(i) = 1._wp / (k     (i) * (1._wp + exp_minus2ktau(i))  + &
                              gamma1(i) * (1._wp - exp_minus2ktau(i)) )

        ! Equation 25
        Rdif(i,j) = RT_term(i) * gamma2(i) * (1._wp - exp_minus2ktau(i))

        ! Equation 26
        Tdif(i,j) = RT_term(i) * 2._wp * k(i) * exp_minusktau(i)
      end do

      !
      ! Transmittance of direct, unscattered beam. Also used below
      !
      Tnoscat(1:ncol,j) = exp(-tau(1:ncol,j)*mu0_inv(1:ncol))

      !
      ! Direct reflect and transmission
      !
      do i = 1, ncol
        k_mu     = k(i) * mu0(i)
        k_gamma3 = k(i) * gamma3(i)
        k_gamma4 = k(i) * gamma4(i)

        !
        ! Equation 14, multiplying top and bottom by exp(-k*tau)
        !   and rearranging to avoid div by 0.
        !
        RT_term(i) =  w0(i,j) * RT_term(i)/merge(1._wp - k_mu*k_mu, &
                                                 epsilon(1._wp),    &
                                                 abs(1._wp - k_mu*k_mu) >= epsilon(1._wp))

        Rdir(i,j) = RT_term(i)  *                                        &
            ((1._wp - k_mu) * (alpha2(i) + k_gamma3)                     - &
             (1._wp + k_mu) * (alpha2(i) - k_gamma3) * exp_minus2ktau(i) - &
             2.0_wp * (k_gamma3 - alpha2(i) * k_mu)  * exp_minusktau (i) * Tnoscat(i,j))

        !
        ! Equation 15, multiplying top and bottom by exp(-k*tau),
        !   multiplying through by exp(-tau/mu0) to
        !   prefer underflow to overflow
        ! Omitting direct transmittance
        !
        Tdir(i,j) = -RT_term(i) *                                                                 &
                    ((1._wp + k_mu) * (alpha1(i) + k_gamma4)                     * Tnoscat(i,j) - &
                     (1._wp - k_mu) * (alpha1(i) - k_gamma4) * exp_minus2ktau(i) * Tnoscat(i,j) - &
                     2.0_wp * (k_gamma4 + alpha1(i) * k_mu)  * exp_minusktau (i))

      end do
    end do
  end subroutine sw_two_stream
  ! ---------------------------------------------------------------
  !
  ! Direct beam source for diffuse radiation in layers and at surface;
  !   report direct beam as a byproduct
  !
  ! -------------------------------------------------------------------------------------------------
  subroutine sw_source_2str(ncol, nlay, top_at_1, Rdir, Tdir, Tnoscat, sfc_albedo, &
                            source_up, source_dn, source_sfc, flux_dn_dir) bind(C, name="sw_source_2str")
    integer,                           intent(in   ) :: ncol, nlay
    logical(wl),                       intent(in   ) :: top_at_1
    real(wp), dimension(ncol, nlay  ), intent(in   ) :: Rdir, Tdir, Tnoscat ! Layer reflectance, transmittance for diffuse radiation
    real(wp), dimension(ncol        ), intent(in   ) :: sfc_albedo          ! surface albedo for direct radiation
    real(wp), dimension(ncol, nlay  ), intent(  out) :: source_dn, source_up
    real(wp), dimension(ncol        ), intent(  out) :: source_sfc          ! Source function for upward radation at surface
    real(wp), dimension(ncol, nlay+1), intent(inout) :: flux_dn_dir ! Direct beam flux
                                                                    ! intent(inout) because top layer includes incident flux

    integer :: ilev

    if(top_at_1) then
      do ilev = 1, nlay
        source_up(:,ilev) =        Rdir(:,ilev) * flux_dn_dir(:,ilev)
        source_dn(:,ilev) =        Tdir(:,ilev) * flux_dn_dir(:,ilev)
        flux_dn_dir(:,ilev+1) = Tnoscat(:,ilev) * flux_dn_dir(:,ilev)
      end do
      source_sfc(:) = flux_dn_dir(:,nlay+1)*sfc_albedo(:)
    else
      ! layer index = level index
      ! previous level is up (+1)
      do ilev = nlay, 1, -1
        source_up(:,ilev)   =    Rdir(:,ilev) * flux_dn_dir(:,ilev+1)
        source_dn(:,ilev)   =    Tdir(:,ilev) * flux_dn_dir(:,ilev+1)
        flux_dn_dir(:,ilev) = Tnoscat(:,ilev) * flux_dn_dir(:,ilev+1)
      end do
      source_sfc(:) = flux_dn_dir(:,     1)*sfc_albedo(:)
    end if
end subroutine sw_source_2str
! ---------------------------------------------------------------
!
! Transport of diffuse radiation through a vertically layered atmosphere.
!   Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
!   This routine is shared by longwave and shortwave
!
! -------------------------------------------------------------------------------------------------
subroutine adding(ncol, nlay, top_at_1, &
                  albedo_sfc,           &
                  rdif, tdif,           &
                  src_dn, src_up, src_sfc, &
                  flux_up, flux_dn) bind(C, name="adding")
  integer,                          intent(in   ) :: ncol, nlay
  logical(wl),                      intent(in   ) :: top_at_1
  real(wp), dimension(ncol       ), intent(in   ) :: albedo_sfc
  real(wp), dimension(ncol,nlay  ), intent(in   ) :: rdif, tdif
  real(wp), dimension(ncol,nlay  ), intent(in   ) :: src_dn, src_up
  real(wp), dimension(ncol       ), intent(in   ) :: src_sfc
  real(wp), dimension(ncol,nlay+1), intent(  out) :: flux_up
  ! intent(inout) because top layer includes incident flux
  real(wp), dimension(ncol,nlay+1), intent(inout) :: flux_dn
  ! ------------------
  integer :: ilev
  real(wp), dimension(ncol,nlay+1)  :: albedo, &  ! reflectivity to diffuse radiation below this level
                                                  ! alpha in SH08
                                       src        ! source of diffuse upwelling radiation from emission or
                                                  ! scattering of direct beam
                                                  ! G in SH08
  real(wp), dimension(ncol,nlay  )  :: denom      ! beta in SH08
  ! ------------------
  !
  ! Indexing into arrays for upward and downward propagation depends on the vertical
  !   orientation of the arrays (whether the domain top is at the first or last index)
  ! We write the loops out explicitly so compilers will have no trouble optimizing them.
  !
  if(top_at_1) then
    ilev = nlay + 1
    ! Albedo of lowest level is the surface albedo...
    albedo(:,ilev)  = albedo_sfc(:)
    ! ... and source of diffuse radiation is surface emission
    src(:,ilev) = src_sfc(:)

    !
    ! From bottom to top of atmosphere --
    !   compute albedo and source of upward radiation
    !
    do ilev = nlay, 1, -1
      denom(:, ilev) = 1._wp/(1._wp - rdif(:,ilev)*albedo(:,ilev+1))                 ! Eq 10
      albedo(:,ilev) = rdif(:,ilev) + &
                       tdif(:,ilev)*tdif(:,ilev) * albedo(:,ilev+1) * denom(:,ilev) ! Equation 9
      !
      ! Equation 11 -- source is emitted upward radiation at top of layer plus
      !   radiation emitted at bottom of layer,
      !   transmitted through the layer and reflected from layers below (tdiff*src*albedo)
      !
      src(:,ilev) =  src_up(:, ilev) + &
                     tdif(:,ilev) * denom(:,ilev) *       &
                       (src(:,ilev+1) + albedo(:,ilev+1)*src_dn(:,ilev))
    end do

    ! Eq 12, at the top of the domain upwelling diffuse is due to ...
    ilev = 1
    flux_up(:,ilev) = flux_dn(:,ilev) * albedo(:,ilev) + & ! ... reflection of incident diffuse and
                      src(:,ilev)                          ! emission from below

    !
    ! From the top of the atmosphere downward -- compute fluxes
    !
    do ilev = 2, nlay+1
      flux_dn(:,ilev) = (tdif(:,ilev-1)*flux_dn(:,ilev-1) + &  ! Equation 13
                         rdif(:,ilev-1)*src(:,ilev) +       &
                         src_dn(:,ilev-1)) * denom(:,ilev-1)
      flux_up(:,ilev) = flux_dn(:,ilev) * albedo(:,ilev) + & ! Equation 12
                        src(:,ilev)
    end do
  else
    ilev = 1
    ! Albedo of lowest level is the surface albedo...
    albedo(:,ilev)  = albedo_sfc(:)
    ! ... and source of diffuse radiation is surface emission
    src(:,ilev) = src_sfc(:)

    !
    ! From bottom to top of atmosphere --
    !   compute albedo and source of upward radiation
    !
    do ilev = 1, nlay
      denom(:, ilev  ) = 1._wp/(1._wp - rdif(:,ilev)*albedo(:,ilev))                ! Eq 10
      albedo(:,ilev+1) = rdif(:,ilev) + &
                         tdif(:,ilev)*tdif(:,ilev) * albedo(:,ilev) * denom(:,ilev) ! Equation 9
      !
      ! Equation 11 -- source is emitted upward radiation at top of layer plus
      !   radiation emitted at bottom of layer,
      !   transmitted through the layer and reflected from layers below (tdiff*src*albedo)
      !
      src(:,ilev+1) =  src_up(:, ilev) +  &
                       tdif(:,ilev) * denom(:,ilev) *       &
                       (src(:,ilev) + albedo(:,ilev)*src_dn(:,ilev))
    end do

    ! Eq 12, at the top of the domain upwelling diffuse is due to ...
    ilev = nlay+1
    flux_up(:,ilev) = flux_dn(:,ilev) * albedo(:,ilev) + & ! ... reflection of incident diffuse and
                      src(:,ilev)                          ! scattering by the direct beam below

    !
    ! From the top of the atmosphere downward -- compute fluxes
    !
    do ilev = nlay, 1, -1
      flux_dn(:,ilev) = (tdif(:,ilev)*flux_dn(:,ilev+1) + &  ! Equation 13
                         rdif(:,ilev)*src(:,ilev) + &
                         src_dn(:, ilev)) * denom(:,ilev)
      flux_up(:,ilev) = flux_dn(:,ilev) * albedo(:,ilev) + & ! Equation 12
                        src(:,ilev)

    end do
  end if
end subroutine adding
! -------------------------------------------------------------------------------------------------
!
! Upper boundary condition
!
! -------------------------------------------------------------------------------------------------
pure subroutine apply_BC_gpt(ncol, nlay, ngpt, top_at_1, inc_flux, flux_dn) bind (C, name="apply_BC_gpt")
  integer,                               intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
  logical(wl),                           intent( in) :: top_at_1
  real(wp), dimension(ncol,       ngpt), intent( in) :: inc_flux         ! Flux at top of domain
  real(wp), dimension(ncol,nlay+1,ngpt), intent(out) :: flux_dn          ! Flux to be used as input to solvers below

  !   Upper boundary condition
  if(top_at_1) then
    flux_dn(1:ncol,      1, 1:ngpt)  = inc_flux(1:ncol,1:ngpt)
  else
    flux_dn(1:ncol, nlay+1, 1:ngpt)  = inc_flux(1:ncol,1:ngpt)
  end if
end subroutine apply_BC_gpt
! ---------------------
pure subroutine apply_BC_factor(ncol, nlay, ngpt, top_at_1, inc_flux, factor, flux_dn) bind (C, name="apply_BC_factor")
  integer,                               intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
  logical(wl),                           intent( in) :: top_at_1
  real(wp), dimension(ncol,       ngpt), intent( in) :: inc_flux         ! Flux at top of domain
  real(wp), dimension(ncol            ), intent( in) :: factor           ! Factor to multiply incoming flux
  real(wp), dimension(ncol,nlay+1,ngpt), intent(out) :: flux_dn          ! Flux to be used as input to solvers below

  !   Upper boundary condition
  if(top_at_1) then
    flux_dn(1:ncol,      1, 1:ngpt)  = inc_flux(1:ncol,1:ngpt) * spread(factor, 2, ngpt)
  else
    flux_dn(1:ncol, nlay+1, 1:ngpt)  = inc_flux(1:ncol,1:ngpt) * spread(factor, 2, ngpt)
  end if
end subroutine apply_BC_factor
! ---------------------
pure subroutine apply_BC_0(ncol, nlay, ngpt, top_at_1, flux_dn) bind (C, name="apply_BC_0")
  integer,                               intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
  logical(wl),                           intent( in) :: top_at_1
  real(wp), dimension(ncol,nlay+1,ngpt), intent(out) :: flux_dn          ! Flux to be used as input to solvers below

  !   Upper boundary condition
  if(top_at_1) then
    flux_dn(1:ncol,      1, 1:ngpt)  = 0._wp
  else
    flux_dn(1:ncol, nlay+1, 1:ngpt)  = 0._wp
  end if
end subroutine apply_BC_0

end module mo_rte_solver_kernels
