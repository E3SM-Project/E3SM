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
    real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: radn_dn ! Radiances [W/m2-str]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: radn_up ! Radiances [W/m2-str]
                                                                           ! Top level must contain incident flux boundary condition

    ! Local variables, no g-point dependency
    real(wp), dimension(ncol,nlay,ngpt) :: tau_loc, &  ! path length (tau/mu)
                                           trans       ! transmissivity  = exp(-tau)
    real(wp), dimension(ncol,nlay,ngpt) :: source_dn, source_up
    real(wp), dimension(ncol,     ngpt) :: source_sfc, sfc_albedo

    real(wp), dimension(:,:,:), pointer :: lev_source_up, lev_source_dn ! Mapping increasing/decreasing indicies to up/down

    real(wp), parameter :: pi = acos(-1._wp)
    integer             :: icol, ilev, igpt, top_level
    ! ------------------------------------


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

    !$acc enter data copyin(d,tau,sfc_src,sfc_emis,lev_source_dec,lev_source_inc,lay_source,radn_dn) async
    !$acc enter data create(tau_loc,trans,source_dn,source_up,source_sfc,sfc_albedo,radn_up) async
    !$acc enter data attach(lev_source_up,lev_source_dn) async

    ! NOTE: This kernel produces small differences between GPU and CPU
    ! implementations on Ascent with PGI, we assume due to floating point
    ! differences in the exp() function. These differences are small in the
    ! RFMIP test case (10^-6).
    !$acc parallel loop collapse(3) async
    do igpt = 1, ngpt
      do ilev = 1, nlay
        do icol = 1, ncol
          !
          ! Optical path and transmission, used in source function and transport calculations
          !
          tau_loc(icol,ilev,igpt) = tau(icol,ilev,igpt)*D(icol,igpt)
          trans  (icol,ilev,igpt) = exp(-tau_loc(icol,ilev,igpt))
        end do
      end do
    end do

    !$acc parallel loop collapse(2) async
    do igpt = 1, ngpt
      do icol = 1, ncol
        !
        ! Transport is for intensity
        !   convert flux at top of domain to intensity assuming azimuthal isotropy
        !
        radn_dn(icol,top_level,igpt) = radn_dn(icol,top_level,igpt)/(2._wp * pi * weight)
        !
        ! Surface albedo, surface source function
        !
        sfc_albedo(icol,igpt) = 1._wp - sfc_emis(icol,igpt)
        source_sfc(icol,igpt) = sfc_emis(icol,igpt) * sfc_src(icol,igpt)
      end do
    end do

    !
    ! Source function for diffuse radiation
    !
    call lw_source_noscat(ncol, nlay, ngpt, &
                          lay_source, lev_source_up, lev_source_dn, &
                          tau_loc, trans, source_dn, source_up)

    !
    ! Transport
    !
    call lw_transport_noscat(ncol, nlay, ngpt, top_at_1,  &
                             tau_loc, trans, sfc_albedo, source_dn, source_up, source_sfc, &
                             radn_up, radn_dn)

    !
    ! Convert intensity to flux assuming azimuthal isotropy and quadrature weight
    !
    !$acc parallel loop collapse(3) async
    do igpt = 1, ngpt
      do ilev = 1, nlay+1
        do icol = 1, ncol
          radn_dn(icol,ilev,igpt) = 2._wp * pi * weight * radn_dn(icol,ilev,igpt)
          radn_up(icol,ilev,igpt) = 2._wp * pi * weight * radn_up(icol,ilev,igpt)
        end do
      end do
    end do

    !$acc exit data copyout(radn_dn,radn_up) async
    !$acc exit data delete(d,tau,sfc_src,sfc_emis,lev_source_dec,lev_source_inc,lay_source,tau_loc,trans,source_dn,source_up,source_sfc,sfc_albedo) async
    !$acc exit data detach(lev_source_up,lev_source_dn) async

  end subroutine lw_solver_noscat
  ! ---------------------------------------------------------------
  !
  ! LW transport, no scattering, multi-angle quadrature
  !   Users provide a set of weights and quadrature angles
  !   Routine sums over single-angle solutions for each sets of angles/weights
  !
  ! ---------------------------------------------------------------
  subroutine lw_solver_noscat_GaussQuad(ncol, nlay, ngpt, top_at_1, nmus, Ds, weights, &
                                   tau, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, flux_up, flux_dn) &
                                   bind (C, name="lw_solver_noscat_GaussQuad")
    integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                           intent(in   ) :: top_at_1
    integer,                               intent(in   ) :: nmus          ! number of quadrature angles
    real(wp), dimension(nmus),             intent(in   ) :: Ds, weights  ! quadrature secants, weights
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: lay_source   ! Planck source at layer average temperature [W/m2]
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: lev_source_inc
                                        ! Planck source at layer edge for radiation in increasing ilay direction [W/m2]
                                        ! Includes spectral weighting that accounts for state-dependent frequency to g-space mapping
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: lev_source_dec
                                               ! Planck source at layer edge for radiation in decreasing ilay direction [W/m2]
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_emis         ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_src          ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: flux_dn ! Radiances [W/m2-str]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: flux_up ! Radiances [W/m2-str]
                                                                           ! Top level must contain incident flux boundary condition
    ! Local variables
    real(wp), dimension(ncol,nlay+1,ngpt) :: radn_dn, radn_up ! Fluxes per quad angle
    real(wp), dimension(ncol,       ngpt) :: Ds_ncol
    real(wp), dimension(ncol,       ngpt) :: flux_top

    integer :: imu, top_level
    integer :: icol, ilev, igpt

    !$acc enter data copyin(Ds,weights,tau,lay_source,lev_source_inc,lev_source_dec,sfc_emis,sfc_src,flux_dn) async
    !$acc enter data create(flux_up,radn_dn,radn_up,Ds_ncol) async

    ! ------------------------------------
    ! ------------------------------------
    !$acc  parallel loop collapse(2) async
    do igpt = 1, ngpt
      do icol = 1, ncol
        Ds_ncol(icol, igpt) = Ds(1)
      end do
    end do

    call lw_solver_noscat(ncol, nlay, ngpt, &
                          top_at_1, Ds_ncol, weights(1), tau, &
                          lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                          flux_up, flux_dn)
    !
    ! For more than one angle use local arrays
    !
    top_level = MERGE(1, nlay+1, top_at_1)
    !$acc parallel loop collapse(2) async
    do igpt = 1,ngpt
      do icol = 1,ncol
        flux_top(icol,igpt) = flux_dn(icol,top_level,igpt)
      end do
    end do
    call apply_BC(ncol, nlay, ngpt, top_at_1, flux_top, radn_dn)

    do imu = 2, nmus
      !$acc  parallel loop collapse(2) async
      do igpt = 1, ngpt
        do icol = 1, ncol
          Ds_ncol(icol, igpt) = Ds(imu)
        end do
      end do
      call lw_solver_noscat(ncol, nlay, ngpt, &
                            top_at_1, Ds_ncol, weights(imu), tau, &
                            lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, &
                            radn_up, radn_dn)
      !$acc  parallel loop collapse(3) async
      do igpt = 1, ngpt
        do ilev = 1, nlay+1
          do icol = 1, ncol
            flux_up(icol,ilev,ngpt) = flux_up(icol,ilev,ngpt) + radn_up(icol,ilev,ngpt)
            flux_dn(icol,ilev,ngpt) = flux_dn(icol,ilev,ngpt) + radn_dn(icol,ilev,ngpt)
          end do
        end do
      end do

    end do ! imu loop

    !$acc exit data copyout(flux_up,flux_dn) async
    !$acc exit data delete(Ds,weights,tau,lay_source,lev_source_inc,lev_source_dec,sfc_emis,sfc_src,radn_dn,radn_up,Ds_ncol) async

    !$acc wait
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
    integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                           intent(in   ) :: top_at_1
    real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau, &  ! Optical thickness,
                                                            ssa, &  ! single-scattering albedo,
                                                            g       ! asymmetry parameter []
    real(wp), dimension(ncol,nlay,ngpt), intent(in   ) :: lay_source   ! Planck source at layer average temperature [W/m2]
    real(wp), dimension(ncol,nlay,ngpt), target, &
                                           intent(in   ) :: lev_source_inc, lev_source_dec
                                        ! Planck source at layer edge for radiation in increasing/decreasing ilay direction [W/m2]
                                        ! Includes spectral weighting that accounts for state-dependent frequency to g-space mapping
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_emis         ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_src          ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), &
                                           intent(  out) :: flux_up   ! Fluxes [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), &
                                           intent(inout) :: flux_dn  ! Fluxes [W/m2]
                                                                              ! Top level (= merge(1, nlay+1, top_at_1)
                                                                              ! must contain incident flux boundary condition
    ! ----------------------------------------------------------------------
    integer :: icol, igpt
    real(wp), dimension(ncol,nlay  ,ngpt) :: Rdif, Tdif, gamma1, gamma2
    real(wp), dimension(ncol       ,ngpt) :: sfc_albedo
    real(wp), dimension(ncol,nlay+1,ngpt) :: lev_source
    real(wp), dimension(ncol,nlay  ,ngpt) :: source_dn, source_up
    real(wp), dimension(ncol       ,ngpt) :: source_sfc
    ! ------------------------------------
    ! ------------------------------------
    !$acc enter data copyin(tau, ssa, g, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, flux_dn)
    !$acc enter data create(flux_up, Rdif, Tdif, gamma1, gamma2, sfc_albedo, lev_source, source_dn, source_up, source_sfc)
    !
    ! RRTMGP provides source functions at each level using the spectral mapping
    !   of each adjacent layer. Combine these for two-stream calculations
    !
    call lw_combine_sources(ncol, nlay, ngpt, top_at_1, &
                            lev_source_inc, lev_source_dec, &
                            lev_source)
    !
    ! Cell properties: reflection, transmission for diffuse radiation
    !   Coupling coefficients needed for source function
    !
    call lw_two_stream(ncol, nlay, ngpt, &
                       tau , ssa, g,     &
                       gamma1, gamma2, Rdif, Tdif)

    !
    ! Source function for diffuse radiation
    !
    call lw_source_2str(ncol, nlay, ngpt, top_at_1, &
                        sfc_emis, sfc_src, &
                        lay_source, lev_source, &
                        gamma1, gamma2, Rdif, Tdif, tau, &
                        source_dn, source_up, source_sfc)

    !$acc  parallel loop collapse(2)
    do igpt = 1, ngpt
      do icol = 1, ncol
        sfc_albedo(icol,igpt) = 1._wp - sfc_emis(icol,igpt)
      end do
    end do
    !
    ! Transport
    !
    call adding(ncol, nlay, ngpt, top_at_1,        &
                sfc_albedo,                        &
                Rdif, Tdif,                        &
                source_dn, source_up, source_sfc,  &
                flux_up, flux_dn)
    !$acc exit data delete(tau, ssa, g, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src)
    !$acc exit data delete(Rdif, Tdif, gamma1, gamma2, sfc_albedo, lev_source, source_dn, source_up, source_sfc)
    !$acc exit data copyout(flux_up, flux_dn)
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
    subroutine sw_solver_noscat(ncol, nlay, ngpt, &
                                top_at_1, tau, mu0, flux_dir) bind (C, name="sw_solver_noscat")
      integer,                    intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
      logical(wl),                intent(in   ) :: top_at_1
      real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau          ! Absorption optical thickness []
      real(wp), dimension(ncol            ), intent(in   ) :: mu0          ! cosine of solar zenith angle
      real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: flux_dir     ! Direct-beam flux, spectral [W/m2]
                                                                           ! Top level must contain incident flux boundary condition
      integer :: icol, ilev, igpt
      real(wp) :: mu0_inv(ncol)
      ! ------------------------------------
      ! ------------------------------------
      !$acc enter data copyin(tau, mu0) create(mu0_inv, flux_dir)
      !$acc parallel loop
      do icol = 1, ncol
        mu0_inv(icol) = 1._wp/mu0(icol)
      enddo
      ! Indexing into arrays for upward and downward propagation depends on the vertical
      !   orientation of the arrays (whether the domain top is at the first or last index)
      ! We write the loops out explicitly so compilers will have no trouble optimizing them.

      ! Downward propagation
      if(top_at_1) then
        ! For the flux at this level, what was the previous level, and which layer has the
        !   radiation just passed through?
        ! layer index = level index - 1
        ! previous level is up (-1)
        !$acc parallel loop collapse(2)
        do igpt = 1, ngpt
          do icol = 1, ncol
            do ilev = 2, nlay+1
              flux_dir(icol,ilev,igpt) = flux_dir(icol,ilev-1,igpt) * exp(-tau(icol,ilev,igpt)*mu0_inv(icol))
            end do
          end do
        end do
      else
        ! layer index = level index
        ! previous level is up (+1)
        !$acc parallel loop collapse(2)
        do igpt = 1, ngpt
          do icol = 1, ncol
            do ilev = nlay, 1, -1
              flux_dir(icol,ilev,igpt) = flux_dir(icol,ilev+1,igpt) * exp(-tau(icol,ilev,igpt)*mu0_inv(icol))
            end do
          end do
        end do
      end if
      !$acc exit data delete(tau, mu0, mu0_inv) copyout(flux_dir)
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
                                  flux_up, flux_dn, flux_dir) bind (C, name="sw_solver_2stream")
      integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
      logical(wl),                           intent(in   ) :: top_at_1
      real(wp), dimension(ncol,nlay,  ngpt), intent(in   ) :: tau, &  ! Optical thickness,
                                                              ssa, &  ! single-scattering albedo,
                                                              g       ! asymmetry parameter []
      real(wp), dimension(ncol            ), intent(in   ) :: mu0     ! cosine of solar zenith angle
      real(wp), dimension(ncol,       ngpt), intent(in   ) :: sfc_alb_dir, sfc_alb_dif
                                                                    ! Spectral albedo of surface to direct and diffuse radiation
      real(wp), dimension(ncol,nlay+1,ngpt), &
                                             intent(  out) :: flux_up ! Fluxes [W/m2]
      real(wp), dimension(ncol,nlay+1,ngpt), &                        ! Downward fluxes contain boundary conditions
                                             intent(inout) :: flux_dn, flux_dir
      ! -------------------------------------------
      integer :: icol, ilay, igpt
      real(wp), dimension(ncol,nlay,ngpt) :: Rdif, Tdif, Rdir, Tdir, Tnoscat
      real(wp), dimension(ncol,nlay,ngpt) :: source_up, source_dn
      real(wp), dimension(ncol     ,ngpt) :: source_srf
      ! ------------------------------------
      !
      ! Cell properties: transmittance and reflectance for direct and diffuse radiation
      !
      !$acc enter data copyin(tau, ssa, g, mu0, sfc_alb_dir, sfc_alb_dif, flux_dn, flux_dir)
      !$acc enter data create(Rdif, Tdif, Rdir, Tdir, Tnoscat, source_up, source_dn, source_srf, flux_up)
      call sw_two_stream(ncol, nlay, ngpt, mu0, &
                         tau , ssa , g   ,      &
                         Rdif, Tdif, Rdir, Tdir, Tnoscat)
      call sw_source_2str(ncol, nlay, ngpt, top_at_1,       &
                          Rdir, Tdir, Tnoscat, sfc_alb_dir, &
                          source_up, source_dn, source_srf, flux_dir)
      call adding(ncol, nlay, ngpt, top_at_1,   &
                  sfc_alb_dif, Rdif, Tdif,      &
                  source_dn, source_up, source_srf, flux_up, flux_dn)
      !
      ! adding computes only diffuse flux; flux_dn is total
      !
      !$acc  parallel loop collapse(3)
      do igpt = 1, ngpt
        do ilay = 1, nlay+1
          do icol = 1, ncol
            flux_dn(icol,ilay,igpt) = flux_dn(icol,ilay,igpt) + flux_dir(icol,ilay,igpt)
          end do
        end do
      end do
      !$acc exit data copyout(flux_up, flux_dn, flux_dir)
      !$acc exit data delete (tau, ssa, g, mu0, sfc_alb_dir, sfc_alb_dif, Rdif, Tdif, Rdir, Tdir, Tnoscat, source_up, source_dn, source_srf)

    end subroutine sw_solver_2stream

    ! -------------------------------------------------------------------------------------------------
    !
    !   Lower-level longwave kernels
    !
    ! ---------------------------------------------------------------
    !
    ! Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
    ! See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13
    !
    ! ---------------------------------------------------------------
    subroutine lw_source_noscat(ncol, nlay, ngpt, lay_source, lev_source_up, lev_source_dn, tau, trans, &
                                source_dn, source_up) bind(C, name="lw_source_noscat")
      integer,                               intent(in) :: ncol, nlay, ngpt
      real(wp), dimension(ncol, nlay, ngpt), intent(in) :: lay_source, & ! Planck source at layer center
                                                           lev_source_up, & ! Planck source at levels (layer edges),
                                                           lev_source_dn, & !   increasing/decreasing layer index
                                                           tau,        & ! Optical path (tau/mu)
                                                           trans         ! Transmissivity (exp(-tau))
      real(wp), dimension(ncol, nlay, ngpt), intent(out):: source_dn, source_up
                                                                     ! Source function at layer edges
                                                                     ! Down at the bottom of the layer, up at the top
      ! --------------------------------
      integer             :: icol, ilay, igpt
      real(wp)            :: fact
      real(wp), parameter :: tau_thresh = sqrt(epsilon(tau))
      ! ---------------------------------------------------------------
      ! ---------------------------------------------------------------
      !$acc  parallel loop collapse(3) async
      do igpt = 1, ngpt
        do ilay = 1, nlay
          do icol = 1, ncol
          !
          ! Weighting factor. Use 2nd order series expansion when rounding error (~tau^2)
          !   is of order epsilon (smallest difference from 1. in working precision)
          !   Thanks to Peter Blossey
          !
          fact = merge((1._wp - trans(icol,ilay,igpt))/tau(icol,ilay,igpt) - trans(icol,ilay,igpt), &
                       tau(icol,ilay,igpt) * ( 0.5_wp - 1._wp/3._wp*tau(icol,ilay,igpt) ), &
                       tau(icol,ilay,igpt) > tau_thresh)
          !
          ! Equation below is developed in Clough et al., 1992, doi:10.1029/92JD01419, Eq 13
          !
          source_dn(icol,ilay,igpt) = (1._wp - trans(icol,ilay,igpt)) * lev_source_dn(icol,ilay,igpt) + &
                                  2._wp * fact * (lay_source(icol,ilay, igpt) - lev_source_dn(icol,ilay,igpt))
          source_up(icol,ilay,igpt) = (1._wp - trans(icol,ilay,igpt)) * lev_source_up(icol,ilay,igpt) + &
                                  2._wp * fact * (lay_source(icol,ilay,igpt) - lev_source_up(icol,ilay,igpt))
          end do
        end do
      end do

    end subroutine lw_source_noscat
  ! ---------------------------------------------------------------
  !
  ! Longwave no-scattering transport
  !
  ! ---------------------------------------------------------------
  subroutine lw_transport_noscat(ncol, nlay, ngpt, top_at_1, &
                                 tau, trans, sfc_albedo, source_dn, source_up, source_sfc, &
                                 radn_up, radn_dn) bind(C, name="lw_transport_noscat")
    integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                           intent(in   ) :: top_at_1   !
    real(wp), dimension(ncol,nlay  ,ngpt), intent(in   ) :: tau, &     ! Absorption optical thickness, pre-divided by mu []
                                                            trans      ! transmissivity = exp(-tau)
    real(wp), dimension(ncol       ,ngpt), intent(in   ) :: sfc_albedo ! Surface albedo
    real(wp), dimension(ncol,nlay  ,ngpt), intent(in   ) :: source_dn, &
                                                            source_up  ! Diffuse radiation emitted by the layer
    real(wp), dimension(ncol       ,ngpt), intent(in   ) :: source_sfc ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: radn_dn ! Radiances [W/m2-str]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: radn_up ! Radiances [W/m2-str]
                                                                             ! Top level must contain incident flux boundary condition
    ! Local variables
    integer :: igpt, ilev, icol
    ! ---------------------------------------------------
    ! ---------------------------------------------------
    if(top_at_1) then
      !
      ! Top of domain is index 1
      !
      !$acc  parallel loop collapse(2) async
      do igpt = 1, ngpt
        do icol = 1, ncol
          ! Downward propagation
          do ilev = 2, nlay+1
            radn_dn(icol,ilev,igpt) = trans(icol,ilev-1,igpt)*radn_dn(icol,ilev-1,igpt) + source_dn(icol,ilev-1,igpt)
          end do

          ! Surface reflection and emission
          radn_up(icol,nlay+1,igpt) = radn_dn(icol,nlay+1,igpt)*sfc_albedo(icol,igpt) + source_sfc(icol,igpt)

          ! Upward propagation
          do ilev = nlay, 1, -1
            radn_up(icol,ilev,igpt) = trans(icol,ilev  ,igpt)*radn_up(icol,ilev+1,igpt) + source_up(icol,ilev,igpt)
          end do
        end do
      end do
    else
      !
      ! Top of domain is index nlay+1
      !
      !$acc  parallel loop collapse(2) async
      do igpt = 1, ngpt
        do icol = 1, ncol
          ! Downward propagation
          do ilev = nlay, 1, -1
            radn_dn(icol,ilev,igpt) = trans(icol,ilev  ,igpt)*radn_dn(icol,ilev+1,igpt) + source_dn(icol,ilev,igpt)
          end do

          ! Surface reflection and emission
          radn_up(icol,     1,igpt) = radn_dn(icol,     1,igpt)*sfc_albedo(icol,igpt) + source_sfc(icol,igpt)

          ! Upward propagation
          do ilev = 2, nlay+1
            radn_up(icol,ilev,igpt) = trans(icol,ilev-1,igpt) * radn_up(icol,ilev-1,igpt) +  source_up(icol,ilev-1,igpt)
          end do
        end do
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
  subroutine lw_two_stream(ncol, nlay, ngpt, tau, w0, g, &
                                gamma1, gamma2, Rdif, Tdif) bind(C, name="lw_two_stream")
    integer,                             intent(in)  :: ncol, nlay, ngpt
    real(wp), dimension(ncol,nlay,ngpt), intent(in)  :: tau, w0, g
    real(wp), dimension(ncol,nlay,ngpt), intent(out) :: gamma1, gamma2, Rdif, Tdif

    ! -----------------------
    integer  :: icol, ilay, igpt

    ! Variables used in Meador and Weaver
    real(wp) :: k

    ! Ancillary variables
    real(wp) :: RT_term
    real(wp) :: exp_minusktau, exp_minus2ktau

    real(wp), parameter :: LW_diff_sec = 1.66  ! 1./cos(diffusivity angle)
    ! ---------------------------------
    ! ---------------------------------
    !$acc enter data copyin(tau, w0, g)
    !$acc enter data create(gamma1, gamma2, Rdif, Tdif)

    !$acc  parallel loop collapse(3)
    do igpt = 1, ngpt
      do ilay = 1, nlay
        do icol = 1, ncol
          !
          ! Coefficients differ from SW implementation because the phase function is more isotropic
          !   Here we follow Fu et al. 1997, doi:10.1175/1520-0469(1997)054<2799:MSPITI>2.0.CO;2
          !   and use a diffusivity sec of 1.66
          !
          gamma1(icol,ilay,igpt)= LW_diff_sec * (1._wp - 0.5_wp * w0(icol,ilay,igpt) * (1._wp + g(icol,ilay,igpt))) ! Fu et al. Eq 2.9
          gamma2(icol,ilay,igpt)= LW_diff_sec *          0.5_wp * w0(icol,ilay,igpt) * (1._wp - g(icol,ilay,igpt))  ! Fu et al. Eq 2.10

          ! Written to encourage vectorization of exponential, square root
          ! Eq 18;  k = SQRT(gamma1**2 - gamma2**2), limited below to avoid div by 0.
          !   k = 0 for isotropic, conservative scattering; this lower limit on k
          !   gives relative error with respect to conservative solution
          !   of < 0.1% in Rdif down to tau = 10^-9
          k = sqrt(max((gamma1(icol,ilay,igpt) - gamma2(icol,ilay,igpt)) * &
                       (gamma1(icol,ilay,igpt) + gamma2(icol,ilay,igpt)),  &
                       1.e-12_wp))
          exp_minusktau = exp(-tau(icol,ilay,igpt)*k)

          !
          ! Diffuse reflection and transmission
          !
          exp_minus2ktau = exp_minusktau * exp_minusktau

          ! Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
          RT_term = 1._wp / (k * (1._wp + exp_minus2ktau)  + &
                    gamma1(icol,ilay,igpt) * (1._wp - exp_minus2ktau) )

          ! Equation 25
          Rdif(icol,ilay,igpt) = RT_term * gamma2(icol,ilay,igpt) * (1._wp - exp_minus2ktau)

          ! Equation 26
          Tdif(icol,ilay,igpt) = RT_term * 2._wp * k * exp_minusktau
        end do
      end do
    end do
    !$acc exit data delete (tau, w0, g)
    !$acc exit data copyout(gamma1, gamma2, Rdif, Tdif)
  end subroutine lw_two_stream
  ! -------------------------------------------------------------------------------------------------
  !
  ! Source function combination
  ! RRTMGP provides two source functions at each level
  !   using the spectral mapping from each of the adjascent layers.
  !   Need to combine these for use in two-stream calculation.
  !
  ! -------------------------------------------------------------------------------------------------
  subroutine lw_combine_sources(ncol, nlay, ngpt, top_at_1, &
                                lev_src_inc, lev_src_dec, lev_source) bind(C, name="lw_combine_sources")
    integer,                                 intent(in ) :: ncol, nlay, ngpt
    logical(wl),                             intent(in ) :: top_at_1
    real(wp), dimension(ncol, nlay  , ngpt), intent(in ) :: lev_src_inc, lev_src_dec
    real(wp), dimension(ncol, nlay+1, ngpt), intent(out) :: lev_source

    integer :: icol, ilay, igpt
    ! ---------------------------------------------------------------
    ! ---------------------------------
    !$acc enter data copyin(lev_src_inc, lev_src_dec)
    !$acc enter data create(lev_source)

    !$acc  parallel loop collapse(3)
    do igpt = 1, ngpt
      do ilay = 1, nlay+1
        do icol = 1,ncol
          if(ilay == 1) then
            lev_source(icol, ilay, igpt) =      lev_src_dec(icol, ilay,   igpt)
          else if (ilay == nlay+1) then
            lev_source(icol, ilay, igpt) =      lev_src_inc(icol, ilay-1, igpt)
          else
            lev_source(icol, ilay, igpt) = sqrt(lev_src_dec(icol, ilay, igpt) * &
                                                lev_src_inc(icol, ilay-1, igpt))
          end if
        end do
      end do
    end do
    !$acc exit data delete (lev_src_inc, lev_src_dec)
    !$acc exit data copyout(lev_source)
  end subroutine lw_combine_sources
  ! ---------------------------------------------------------------
  !
  ! Compute LW source function for upward and downward emission at levels using linear-in-tau assumption
  !   This version straight from ECRAD
  !   Source is provided as W/m2-str; factor of pi converts to flux units
  !
  ! ---------------------------------------------------------------
  subroutine lw_source_2str(ncol, nlay, ngpt, top_at_1,   &
                            sfc_emis, sfc_src,      &
                            lay_source, lev_source, &
                            gamma1, gamma2, rdif, tdif, tau, source_dn, source_up, source_sfc) &
                            bind (C, name="lw_source_2str")
    integer,                         intent(in) :: ncol, nlay, ngpt
    logical(wl),                     intent(in) :: top_at_1
    real(wp), dimension(ncol      , ngpt), intent(in) :: sfc_emis, sfc_src
    real(wp), dimension(ncol, nlay, ngpt), intent(in) :: lay_source,    & ! Planck source at layer center
                                                   tau,           & ! Optical depth (tau)
                                                   gamma1, gamma2,& ! Coupling coefficients
                                                   rdif, tdif       ! Layer reflectance and transmittance
    real(wp), dimension(ncol, nlay+1, ngpt), target, &
                                     intent(in)  :: lev_source       ! Planck source at layer edges
    real(wp), dimension(ncol, nlay, ngpt), intent(out) :: source_dn, source_up
    real(wp), dimension(ncol      , ngpt), intent(out) :: source_sfc      ! Source function for upward radation at surface

    integer             :: icol, ilay, igpt
    real(wp)            :: Z, Zup_top, Zup_bottom, Zdn_top, Zdn_bottom
    real(wp)            :: lev_source_bot, lev_source_top
    ! ---------------------------------------------------------------
    ! ---------------------------------
    !$acc enter data copyin(sfc_emis, sfc_src, lay_source, tau, gamma1, gamma2, rdif, tdif, lev_source)
    !$acc enter data create(source_dn, source_up, source_sfc)

    !$acc parallel loop collapse(3)
    do igpt = 1, ngpt
      do ilay = 1, nlay
        do icol = 1, ncol
          if (tau(icol,ilay,ngpt) > 1.0e-8_wp) then
            if(top_at_1) then
              lev_source_top = lev_source(icol,ilay  ,ngpt)
              lev_source_bot = lev_source(icol,ilay+1,ngpt)
            else
              lev_source_top = lev_source(icol,ilay+1,ngpt)
              lev_source_bot = lev_source(icol,ilay  ,ngpt)
            end if
            !
            ! Toon et al. (JGR 1989) Eqs 26-27
            !
            Z = (lev_source_bot-lev_source_top) / (tau(icol,ilay,igpt)*(gamma1(icol,ilay,igpt)+gamma2(icol,ilay,igpt)))
            Zup_top        =  Z + lev_source_top
            Zup_bottom     =  Z + lev_source_bot
            Zdn_top        = -Z + lev_source_top
            Zdn_bottom     = -Z + lev_source_bot
            source_up(icol,ilay,igpt) = pi * (Zup_top    - rdif(icol,ilay,igpt) * Zdn_top    - tdif(icol,ilay,igpt) * Zup_bottom)
            source_dn(icol,ilay,igpt) = pi * (Zdn_bottom - rdif(icol,ilay,igpt) * Zup_bottom - tdif(icol,ilay,igpt) * Zdn_top)
          else
            source_up(icol,ilay,igpt) = 0._wp
            source_dn(icol,ilay,igpt) = 0._wp
          end if
          if(ilay == 1) source_sfc(icol,igpt) = pi * sfc_emis(icol,igpt) * sfc_src(icol,igpt)
        end do
      end do
    end do
    !$acc exit data delete(sfc_emis, sfc_src, lay_source, tau, gamma1, gamma2, rdif, tdif, lev_source)
    !$acc exit data copyout(source_dn, source_up, source_sfc)

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
    subroutine sw_two_stream(ncol, nlay, ngpt, mu0, tau, w0, g, &
                                  Rdif, Tdif, Rdir, Tdir, Tnoscat) bind (C, name="sw_two_stream")
      integer,                             intent(in)  :: ncol, nlay, ngpt
      real(wp), dimension(ncol),           intent(in)  :: mu0
      real(wp), dimension(ncol,nlay,ngpt), intent(in)  :: tau, w0, g
      real(wp), dimension(ncol,nlay,ngpt), intent(out) :: Rdif, Tdif, Rdir, Tdir, Tnoscat

      ! -----------------------
      integer  :: icol,ilay,igpt

      ! Variables used in Meador and Weaver
      real(wp) :: gamma1, gamma2, gamma3, gamma4
      real(wp) :: alpha1, alpha2, k

      ! Ancillary variables
      real(wp) :: RT_term
      real(wp) :: exp_minusktau, exp_minus2ktau
      real(wp) :: k_mu, k_gamma3, k_gamma4
      real(wp) :: mu0_inv(ncol)
      ! ---------------------------------
      ! ---------------------------------
      !$acc enter data copyin (mu0, tau, w0, g)
      !$acc enter data create(Rdif, Tdif, Rdir, Tdir, Tnoscat, mu0_inv)

      !$acc parallel loop
      do icol = 1, ncol
        mu0_inv(icol) = 1._wp/mu0(icol)
      enddo

      ! NOTE: this kernel appears to cause small (10^-6) differences between GPU
      ! and CPU. This *might* be floating point differences in implementation of
      ! the exp function.
      !$acc  parallel loop collapse(3)
      do igpt = 1, ngpt
        do ilay = 1, nlay
          do icol = 1, ncol
            ! Zdunkowski Practical Improved Flux Method "PIFM"
            !  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
            !
            gamma1= (8._wp - w0(icol,ilay,igpt) * (5._wp + 3._wp * g(icol,ilay,igpt))) * .25_wp
            gamma2=  3._wp *(w0(icol,ilay,igpt) * (1._wp -         g(icol,ilay,igpt))) * .25_wp
            gamma3= (2._wp - 3._wp * mu0(icol)  *                  g(icol,ilay,igpt) ) * .25_wp
            gamma4=  1._wp - gamma3

            alpha1 = gamma1 * gamma4 + gamma2 * gamma3           ! Eq. 16
            alpha2 = gamma1 * gamma3 + gamma2 * gamma4           ! Eq. 17
            ! Written to encourage vectorization of exponential, square root
            ! Eq 18;  k = SQRT(gamma1**2 - gamma2**2), limited below to avoid div by 0.
            !   k = 0 for isotropic, conservative scattering; this lower limit on k
            !   gives relative error with respect to conservative solution
            !   of < 0.1% in Rdif down to tau = 10^-9
            k = sqrt(max((gamma1 - gamma2) * &
                         (gamma1 + gamma2),  &
                         1.e-12_wp))
            exp_minusktau = exp(-tau(icol,ilay,igpt)*k)
            !
            ! Diffuse reflection and transmission
            !
            exp_minus2ktau = exp_minusktau * exp_minusktau

            ! Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes
            RT_term = 1._wp / (k      * (1._wp + exp_minus2ktau)  + &
                               gamma1 * (1._wp - exp_minus2ktau) )

            ! Equation 25
            Rdif(icol,ilay,igpt) = RT_term * gamma2 * (1._wp - exp_minus2ktau)

            ! Equation 26
            Tdif(icol,ilay,igpt) = RT_term * 2._wp * k * exp_minusktau

            !
            ! Transmittance of direct, unscattered beam. Also used below
            !
            Tnoscat(icol,ilay,igpt) = exp(-tau(icol,ilay,igpt)*mu0_inv(icol))

            !
            ! Direct reflect and transmission
            !
            k_mu     = k * mu0(icol)
            k_gamma3 = k * gamma3
            k_gamma4 = k * gamma4

            !
            ! Equation 14, multiplying top and bottom by exp(-k*tau)
            !   and rearranging to avoid div by 0.
            !
            RT_term =  w0(icol,ilay,igpt) * RT_term/merge(1._wp - k_mu*k_mu, &
                                                         epsilon(1._wp),    &
                                                         abs(1._wp - k_mu*k_mu) >= epsilon(1._wp))

            Rdir(icol,ilay,igpt) = RT_term  *                                    &
               ((1._wp - k_mu) * (alpha2 + k_gamma3)                  - &
                (1._wp + k_mu) * (alpha2 - k_gamma3) * exp_minus2ktau - &
                2.0_wp * (k_gamma3 - alpha2 * k_mu)  * exp_minusktau  * Tnoscat(icol,ilay,igpt))

            !
            ! Equation 15, multiplying top and bottom by exp(-k*tau),
            !   multiplying through by exp(-tau/mu0) to
            !   prefer underflow to overflow
            ! Omitting direct transmittance
            !
            Tdir(icol,ilay,igpt) = &
                     -RT_term * ((1._wp + k_mu) * (alpha1 + k_gamma4) * Tnoscat(icol,ilay,igpt) - &
                                 (1._wp - k_mu) * (alpha1 - k_gamma4) * exp_minus2ktau * Tnoscat(icol,ilay,igpt) - &
                                  2.0_wp * (k_gamma4 + alpha1 * k_mu)  * exp_minusktau )

          end do
        end do
      end do
      !$acc exit data delete (mu0, tau, w0, g, mu0_inv)
      !$acc exit data copyout(Rdif, Tdif, Rdir, Tdir, Tnoscat)

    end subroutine sw_two_stream
  ! ---------------------------------------------------------------
  !
  ! Direct beam source for diffuse radiation in layers and at surface;
  !   report direct beam as a byproduct
  !
  subroutine sw_source_2str(ncol, nlay, ngpt, top_at_1, Rdir, Tdir, Tnoscat, sfc_albedo, &
                            source_up, source_dn, source_sfc, flux_dn_dir) bind(C, name="sw_source_2str")
    integer,                                 intent(in   ) :: ncol, nlay, ngpt
    logical(wl),                             intent(in   ) :: top_at_1
    real(wp), dimension(ncol, nlay  , ngpt), intent(in   ) :: Rdir, Tdir, Tnoscat ! Layer reflectance, transmittance for diffuse radiation
    real(wp), dimension(ncol        , ngpt), intent(in   ) :: sfc_albedo          ! surface albedo for direct radiation
    real(wp), dimension(ncol, nlay  , ngpt), intent(  out) :: source_dn, source_up
    real(wp), dimension(ncol        , ngpt), intent(  out) :: source_sfc          ! Source function for upward radation at surface
    real(wp), dimension(ncol, nlay+1, ngpt), intent(inout) :: flux_dn_dir ! Direct beam flux
                                                                    ! intent(inout) because top layer includes incident flux

    integer :: icol, ilev, igpt
    ! ---------------------------------
    ! ---------------------------------
    !$acc enter data copyin (Rdir, Tdir, Tnoscat, sfc_albedo, flux_dn_dir)
    !$acc enter data create(source_dn, source_up, source_sfc)

    if(top_at_1) then
      !$acc  parallel loop collapse(2)
      do igpt = 1, ngpt
        do icol = 1, ncol
          do ilev = 1, nlay
            source_up(icol,ilev,igpt)     =    Rdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt)
            source_dn(icol,ilev,igpt)     =    Tdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt)
            flux_dn_dir(icol,ilev+1,igpt) = Tnoscat(icol,ilev,igpt) * flux_dn_dir(icol,ilev,igpt)
            if(ilev == nlay) source_sfc(icol,igpt) = flux_dn_dir(icol,nlay+1,igpt)*sfc_albedo(icol,igpt)
          end do
        end do
      end do
    else
      ! layer index = level index
      ! previous level is up (+1)
      !$acc  parallel loop collapse(2)
      do igpt = 1, ngpt
        do icol = 1, ncol
          do ilev = nlay, 1, -1
            source_up(icol,ilev,igpt)   =    Rdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt)
            source_dn(icol,ilev,igpt)   =    Tdir(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt)
            flux_dn_dir(icol,ilev,igpt) = Tnoscat(icol,ilev,igpt) * flux_dn_dir(icol,ilev+1,igpt)
            if(ilev ==    1) source_sfc(icol,igpt) = flux_dn_dir(icol,    1,igpt)*sfc_albedo(icol,igpt)
          end do
        end do
      end do
    end if
    !$acc exit data copyout(source_dn, source_up, source_sfc, flux_dn_dir)
    !$acc exit data delete(Rdir, Tdir, Tnoscat, sfc_albedo)

  end subroutine sw_source_2str
! ---------------------------------------------------------------
!
! Transport of diffuse radiation through a vertically layered atmosphere.
!   Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
!   This routine is shared by longwave and shortwave
!
! -------------------------------------------------------------------------------------------------
  subroutine adding(ncol, nlay, ngpt, top_at_1, &
                    albedo_sfc,           &
                    rdif, tdif,           &
                    src_dn, src_up, src_sfc, &
                    flux_up, flux_dn) bind(C, name="adding")
    integer,                               intent(in   ) :: ncol, nlay, ngpt
    logical(wl),                           intent(in   ) :: top_at_1
    real(wp), dimension(ncol       ,ngpt), intent(in   ) :: albedo_sfc
    real(wp), dimension(ncol,nlay  ,ngpt), intent(in   ) :: rdif, tdif
    real(wp), dimension(ncol,nlay  ,ngpt), intent(in   ) :: src_dn, src_up
    real(wp), dimension(ncol       ,ngpt), intent(in   ) :: src_sfc
    real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: flux_up
    ! intent(inout) because top layer includes incident flux
    real(wp), dimension(ncol,nlay+1,ngpt), intent(inout) :: flux_dn
    ! ------------------
    integer :: icol, ilev, igpt
    real(wp), dimension(nlay+1) :: albedo, &  ! reflectivity to diffuse radiation below this level
                                              ! alpha in SH08
                                   src        ! source of diffuse upwelling radiation from emission or
                                              ! scattering of direct beam
                                              ! G in SH08
    real(wp), dimension(nlay  ) :: denom      ! beta in SH08
    ! ------------------
    ! ---------------------------------
    !
    ! Indexing into arrays for upward and downward propagation depends on the vertical
    !   orientation of the arrays (whether the domain top is at the first or last index)
    ! We write the loops out explicitly so compilers will have no trouble optimizing them.
    !
    !$acc enter data copyin(albedo_sfc, rdif, tdif, src_dn, src_up, src_sfc, flux_dn)
    !$acc enter data create(flux_up, albedo, src, denom)
    if(top_at_1) then
#ifdef __PGI
      !$acc parallel loop
#else
      !$acc parallel loop collapse(2) private(albedo, src, denom)
#endif
      do igpt = 1, ngpt
#ifdef __PGI
        !$acc loop private(albedo, src, denom)
#endif
        do icol = 1, ncol
          ilev = nlay + 1
          ! Albedo of lowest level is the surface albedo...
          albedo(ilev)  = albedo_sfc(icol,igpt)
          ! ... and source of diffuse radiation is surface emission
          src(ilev) = src_sfc(icol,igpt)

          !
          ! From bottom to top of atmosphere --
          !   compute albedo and source of upward radiation
          !
          do ilev = nlay, 1, -1
            denom(ilev) = 1._wp/(1._wp - rdif(icol,ilev,igpt)*albedo(ilev+1))                 ! Eq 10
            albedo(ilev) = rdif(icol,ilev,igpt) + &
                           tdif(icol,ilev,igpt)*tdif(icol,ilev,igpt) * albedo(ilev+1) * denom(ilev) ! Equation 9
            !
            ! Equation 11 -- source is emitted upward radiation at top of layer plus
            !   radiation emitted at bottom of layer,
            !   transmitted through the layer and reflected from layers below (tdiff*src*albedo)
            !
            src(ilev) =  src_up(icol, ilev, igpt) + &
                           tdif(icol,ilev,igpt) * denom(ilev) *       &
                             (src(ilev+1) + albedo(ilev+1)*src_dn(icol,ilev,igpt))
          end do

          ! Eq 12, at the top of the domain upwelling diffuse is due to ...
          ilev = 1
          flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(ilev) + & ! ... reflection of incident diffuse and
                                    src(ilev)                                  ! emission from below

          !
          ! From the top of the atmosphere downward -- compute fluxes
          !
          do ilev = 2, nlay+1
            flux_dn(icol,ilev,igpt) = (tdif(icol,ilev-1,igpt)*flux_dn(icol,ilev-1,igpt) + &  ! Equation 13
                               rdif(icol,ilev-1,igpt)*src(ilev) +       &
                               src_dn(icol,ilev-1,igpt)) * denom(ilev-1)
            flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(ilev) + & ! Equation 12
                              src(ilev)
          end do
        end do
      end do
    else
#ifdef __PGI
      !$acc parallel loop
#else
      !$acc parallel loop collapse(2) private(albedo, src, denom)
#endif
      do igpt = 1, ngpt
#ifdef __PGI
        !$acc loop private(albedo, src, denom)
#endif
        do icol = 1, ncol
          ilev = 1
          ! Albedo of lowest level is the surface albedo...
          albedo(ilev)  = albedo_sfc(icol,igpt)
          ! ... and source of diffuse radiation is surface emission
          src(ilev) = src_sfc(icol,igpt)

          !
          ! From bottom to top of atmosphere --
          !   compute albedo and source of upward radiation
          !
          do ilev = 1, nlay
            denom(ilev  ) = 1._wp/(1._wp - rdif(icol,ilev,igpt)*albedo(ilev))                ! Eq 10
            albedo(ilev+1) = rdif(icol,ilev,igpt) + &
                               tdif(icol,ilev,igpt)*tdif(icol,ilev,igpt) * albedo(ilev) * denom(ilev) ! Equation 9
            !
            ! Equation 11 -- source is emitted upward radiation at top of layer plus
            !   radiation emitted at bottom of layer,
            !   transmitted through the layer and reflected from layers below (tdiff*src*albedo)
            !
            src(ilev+1) =  src_up(icol, ilev, igpt) +  &
                             tdif(icol,ilev,igpt) * denom(ilev) *       &
                             (src(ilev) + albedo(ilev)*src_dn(icol,ilev,igpt))
          end do

          ! Eq 12, at the top of the domain upwelling diffuse is due to ...
          ilev = nlay+1
          flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(ilev) + & ! ... reflection of incident diffuse and
                            src(ilev)                          ! scattering by the direct beam below

          !
          ! From the top of the atmosphere downward -- compute fluxes
          !
          do ilev = nlay, 1, -1
            flux_dn(icol,ilev,igpt) = (tdif(icol,ilev,igpt)*flux_dn(icol,ilev+1,igpt) + &  ! Equation 13
                               rdif(icol,ilev,igpt)*src(ilev) + &
                               src_dn(icol, ilev, igpt)) * denom(ilev)
            flux_up(icol,ilev,igpt) = flux_dn(icol,ilev,igpt) * albedo(ilev) + & ! Equation 12
                              src(ilev)

          end do
        end do
      end do
    end if
    !$acc exit data delete(albedo_sfc, rdif, tdif, src_dn, src_up, src_sfc, albedo, src, denom)
    !$acc exit data copyout(flux_up, flux_dn)
  end subroutine adding
  ! ---------------------------------------------------------------
  !
  ! Upper boundary condition
  !
  ! ---------------------------------------------------------------
  subroutine apply_BC_gpt(ncol, nlay, ngpt, top_at_1, inc_flux, flux_dn) bind (C, name="apply_BC_gpt")
    integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                           intent(in   ) :: top_at_1
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: inc_flux         ! Flux at top of domain
    real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: flux_dn          ! Flux to be used as input to solvers below

    integer :: icol, igpt
    ! --------------
    ! --------------
    !   Upper boundary condition
    if(top_at_1) then
      !$acc  parallel loop collapse(2) async
      do igpt = 1, ngpt
        do icol = 1, ncol
          flux_dn(icol,      1, igpt)  = inc_flux(icol,igpt)
        end do
      end do
    else
      !$acc  parallel loop collapse(2) async
      do igpt = 1, ngpt
        do icol = 1, ncol
          flux_dn(icol, nlay+1, igpt)  = inc_flux(icol,igpt)
        end do
      end do
    end if
    !$acc wait
  end subroutine apply_BC_gpt
  ! ---------------------
  subroutine apply_BC_factor(ncol, nlay, ngpt, top_at_1, inc_flux, factor, flux_dn) bind (C, name="apply_BC_factor")
    integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                           intent(in   ) :: top_at_1
    real(wp), dimension(ncol,       ngpt), intent(in   ) :: inc_flux         ! Flux at top of domain
    real(wp), dimension(ncol            ), intent(in   ) :: factor           ! Factor to multiply incoming flux
    real(wp), dimension(ncol,nlay+1,ngpt), intent(out  ) :: flux_dn          ! Flux to be used as input to solvers below

    integer :: icol, igpt
    ! --------------
    ! --------------

    !   Upper boundary condition
    if(top_at_1) then
      !$acc  parallel loop collapse(2) async
      do igpt = 1, ngpt
        do icol = 1, ncol
          flux_dn(icol,      1, igpt)  = inc_flux(icol,igpt) * factor(icol)
        end do
      end do
    else
      !$acc  parallel loop collapse(2) async
      do igpt = 1, ngpt
        do icol = 1, ncol
          flux_dn(icol, nlay+1, igpt)  = inc_flux(icol,igpt) * factor(icol)
        end do
      end do
    end if
    !$acc wait
  end subroutine apply_BC_factor
  ! ---------------------
  subroutine apply_BC_0(ncol, nlay, ngpt, top_at_1, flux_dn) bind (C, name="apply_BC_0")
    integer,                               intent(in   ) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical(wl),                           intent(in   ) :: top_at_1
    real(wp), dimension(ncol,nlay+1,ngpt), intent(  out) :: flux_dn          ! Flux to be used as input to solvers below

    integer :: icol, igpt
    ! --------------
    ! --------------

    !   Upper boundary condition
    if(top_at_1) then
      !$acc  parallel loop collapse(2) async
      do igpt = 1, ngpt
        do icol = 1, ncol
          flux_dn(icol,      1, igpt)  = 0._wp
        end do
      end do
    else
      !$acc  parallel loop collapse(2) async
      do igpt = 1, ngpt
        do icol = 1, ncol
          flux_dn(icol, nlay+1, igpt)  = 0._wp
        end do
      end do
    end if
    !$acc wait
  end subroutine apply_BC_0

end module mo_rte_solver_kernels
