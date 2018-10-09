! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description: Numeric calculations for gas optics. Absorption and Rayleigh optical depths,
!   source functions.

module mo_gas_optics_kernels
  use mo_rte_kind,      only: wp
  use mo_util_string,   only : string_loc_in_array
  implicit none

  interface zero_array
    module procedure zero_array_3D, zero_array_4D
  end interface
contains
  ! --------------------------------------------------------------------------------------

  !
  ! Compute minor and major species opitcal depth from pre-computed interpolation coefficients
  !   (jeta,jtemp,jpress)
  !
  subroutine compute_tau_absorption(                &
                ncol,nlay,ngpt,ngas,nflav,          &  ! dimensions
                idx_h2o,                            &
                gpoint_flavor,                      &
                kmajor,                             &
                kminor_lower,                       &
                kminor_upper,                       &
                minor_limits_gpt_lower,             &
                minor_limits_gpt_upper,             &
                minor_scales_with_density_lower,    &
                minor_scales_with_density_upper,    &
                scale_by_complement_lower,          &
                scale_by_complement_upper,          &
                idx_minor_lower,                    &
                idx_minor_upper,                    &
                idx_minor_scaling_lower,            &
                idx_minor_scaling_upper,            &
                kminor_start_lower,                 &
                kminor_start_upper,                 &
                tropo,itropo_lower,itropo_upper,    &
                col_mix,fmajor,fminor,              &
                play,tlay,col_gas,                  &
                jeta,jtemp,jpress,                  &
                tau)
    ! ---------------------
    ! input dimensions
    integer,                                intent(in) :: ncol,nlay,ngpt,ngas,nflav
    integer,                                intent(in) :: idx_h2o
    ! ---------------------
    ! inputs from object
    integer,          dimension(2,ngpt),    intent(in) :: gpoint_flavor
    real(wp),         dimension(:,:,:,:),   intent(in) :: kmajor
    real(wp),         dimension(:,:,:),     intent(in) :: kminor_lower
    real(wp),         dimension(:,:,:),     intent(in) :: kminor_upper
    integer,          dimension(:,:),       intent(in) :: minor_limits_gpt_lower
    integer,          dimension(:,:),       intent(in) :: minor_limits_gpt_upper
    logical,          dimension(:),         intent(in) :: minor_scales_with_density_lower
    logical,          dimension(:),         intent(in) :: minor_scales_with_density_upper
    logical,          dimension(:),         intent(in) :: scale_by_complement_lower
    logical,          dimension(:),         intent(in) :: scale_by_complement_upper
    integer,          dimension(:),         intent(in) :: idx_minor_lower
    integer,          dimension(:),         intent(in) :: idx_minor_upper
    integer,          dimension(:),         intent(in) :: idx_minor_scaling_lower
    integer,          dimension(:),         intent(in) :: idx_minor_scaling_upper
    integer,          dimension(:),         intent(in) :: kminor_start_lower
    integer,          dimension(:),         intent(in) :: kminor_start_upper
    logical,          dimension(ncol,nlay), intent(in) :: tropo
    integer,          dimension(ncol,2),    intent(in) :: itropo_lower
    integer,          dimension(ncol,2),    intent(in) :: itropo_upper
    ! ---------------------
    ! inputs from profile or parent function
    real(wp), dimension(2,    nflav,ncol,nlay       ), intent(in) :: col_mix
    real(wp), dimension(2,2,2,nflav,ncol,nlay       ), intent(in) :: fmajor
    real(wp), dimension(2,2,  nflav,ncol,nlay       ), intent(in) :: fminor
    real(wp), dimension(            ncol,nlay       ), intent(in) :: play, tlay      ! pressure and temperature
    real(wp), dimension(            ncol,nlay,0:ngas), intent(in) :: col_gas
    integer,  dimension(2,    nflav,ncol,nlay       ), intent(in) :: jeta
    integer,  dimension(             ncol,nlay      ), intent(in) :: jtemp
    integer,  dimension(             ncol,nlay      ), intent(in) :: jpress
    ! ---------------------
    ! output - optical depth
    real(wp), dimension(ngpt,nlay,ncol), intent(inout) :: tau
    ! ---------------------

    ! ---------------------
    ! Major Species
    ! ---------------------
    call gas_optical_depths_major( &
          ncol,nlay,ngpt,nflav,           & ! dimensions
          gpoint_flavor,                  & ! inputs from object
          kmajor,                         &
          col_mix,fmajor,                 &
          jeta,tropo,jtemp,jpress,        & ! local input
          tau)
    ! ---------------------
    ! Minor Species - lower
    ! ---------------------
    call gas_optical_depths_minor(     &
           ncol,nlay,ngpt,ngas,nflav,  & ! dimensions
           idx_h2o,                    &
           gpoint_flavor(1,:),         &
           kminor_lower,               &
           minor_limits_gpt_lower,     &
           minor_scales_with_density_lower, &
           scale_by_complement_lower,  &
           idx_minor_lower,            &
           idx_minor_scaling_lower,    &
           kminor_start_lower,         &
           play, tlay,                 &
           col_gas,fminor,jeta,        &
           itropo_lower,jtemp,         &
           tau)
    ! ---------------------
    ! Minor Species - upper
    ! ---------------------
    call gas_optical_depths_minor(     &
           ncol,nlay,ngpt,ngas,nflav,  & ! dimensions
           idx_h2o,                    &
           gpoint_flavor(2,:),         &
           kminor_upper,               &
           minor_limits_gpt_upper,     &
           minor_scales_with_density_upper, &
           scale_by_complement_upper,  &
           idx_minor_upper,            &
           idx_minor_scaling_upper,    &
           kminor_start_upper,         &
           play, tlay,                 &
           col_gas,fminor,jeta,        &
           itropo_upper,jtemp,         &
           tau)

  end subroutine compute_tau_absorption
  ! --------------------------------------------------------------------------------------

  ! --------------------------------------------------------------------------------------
  !
  ! compute minor species optical depths
  !
  subroutine gas_optical_depths_major(ncol,nlay,ngpt,nflav,           & ! dimensions
                                      gpoint_flavor,                  & ! inputs from object
                                      kmajor,                         &
                                      col_mix,fmajor,                 &
                                      jeta,tropo,jtemp,jpress,        & ! local input
                                      tau)
    ! input dimensions
    integer, intent(in) :: ncol, nlay, ngpt, nflav ! dimensions

    ! inputs from object
    integer,  dimension(2,ngpt),  intent(in) :: gpoint_flavor
    real(wp), dimension(:,:,:,:), intent(in) :: kmajor

    ! inputs from profile or parent function
    real(wp), dimension(2,    nflav,ncol,nlay), intent(in) :: col_mix
    real(wp), dimension(2,2,2,nflav,ncol,nlay), intent(in) :: fmajor
    integer,  dimension(2,    nflav,ncol,nlay), intent(in) :: jeta
    logical,  dimension(ncol,nlay), intent(in) :: tropo
    integer,  dimension(ncol,nlay), intent(in) :: jtemp, jpress

    ! outputs
    real(wp), dimension(ngpt,nlay,ncol), intent(inout) :: tau
    ! -----------------
    ! local variables
    real(wp) :: tau_major ! major species optical depth
    ! local index
    integer :: icol, ilay, iflav, igpt, itropo

    ! -----------------

    do icol = 1, ncol
      do ilay = 1, nlay
        ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        itropo = merge(1,2,tropo(icol,ilay))
        ! optical depth calculation for major species
        do igpt = 1, ngpt
          ! binary species parameter (eta) and col_mix depend on band flavor
          iflav = gpoint_flavor(itropo, igpt)
          tau_major = &
            ! interpolation in temperature, pressure, and eta
            interpolate3D(col_mix(:,iflav,icol,ilay), &
                          fmajor(:,:,:,iflav,icol,ilay), kmajor, &
                          igpt, jeta(:,iflav,icol,ilay), jtemp(icol,ilay),jpress(icol,ilay)+itropo)
          tau(igpt,ilay,icol) = tau(igpt,ilay,icol) + tau_major
        end do ! igpt
      end do
    end do ! ilay
  end subroutine gas_optical_depths_major

  ! ----------------------------------------------------------
  !
  ! compute minor species optical depths
  !
  subroutine gas_optical_depths_minor(ncol, nlay,ngpt,ngas,nflav,   &
                                      idx_h2o,            &
                                      gpt_flv,            &
                                      kminor,             &
                                      minor_limits_gpt,   &
                                      minor_scales_with_density,    &
                                      scale_by_complement, &
                                      idx_minor, idx_minor_scaling, &
                                      kminor_start,        &
                                      play, tlay,          &
                                      col_gas,fminor,jeta, &
                                      itropo,jtemp,        &
                                      tau)
    integer,                                  intent(in ) :: ncol,nlay,ngpt,ngas,nflav
    integer,                                  intent(in ) :: idx_h2o
    integer,  dimension(:),                   intent(in ) :: gpt_flv
    real(wp), dimension(:,:,:),               intent(in ) :: kminor
    integer,  dimension(:,:),                 intent(in ) :: minor_limits_gpt
    logical,  dimension(:),                   intent(in ) :: minor_scales_with_density
    logical,  dimension(:),                   intent(in ) :: scale_by_complement
    integer,  dimension(:),                   intent(in ) :: kminor_start
    integer,  dimension(:),                   intent(in ) :: idx_minor, idx_minor_scaling
    real(wp), dimension(ncol,nlay),           intent(in ) :: play, tlay
    real(wp), dimension(ncol,nlay,0:ngas),    intent(in ) :: col_gas
    real(wp), dimension(2,2,nflav,ncol,nlay), intent(in ) :: fminor
    integer,  dimension(2,  nflav,ncol,nlay), intent(in ) :: jeta
    integer,  dimension(ncol, 2),             intent(in ) :: itropo
    integer,  dimension(ncol,nlay),           intent(in ) :: jtemp
    real(wp), dimension(ngpt,nlay,ncol),      intent(out) :: tau
    ! -----------------
    ! local variables
    real(wp), parameter :: PaTohPa = 0.01
    real(wp), dimension(ngas) :: vmr
    real(wp) :: scaling, kminor_loc, tau_minor ! minor species absorption coefficient, optical depth
    integer  :: icol, ilay, iflav, igpt, imnr
    integer  :: itl, itu, iml, imu
    integer  :: minor_start, minor_loc
    ! -----------------
    do imnr = 1, size(scale_by_complement,dim=1) ! loop over minor absorbers in each band
      do icol = 1, ncol
        ! Get layer range
        itl = itropo(icol,1)
        itu = itropo(icol,2)
        do ilay = itl,itu
          vmr(1:ngas) = col_gas(icol,ilay,1:ngas)/col_gas(icol,ilay,0)
          !
          ! Scaling of minor gas absortion coefficient begins with column amount of minor gas
          !
          scaling = col_gas(icol,ilay,idx_minor(imnr))
          !
          ! Density scaling (e.g. for h2o continuum, collision-induced absorption)
          !
          if (minor_scales_with_density(imnr)) then
            !
            ! NOTE: P needed in hPa to properly handle density scaling.
            !
            scaling = scaling * (PaTohPa*play(icol,ilay)/tlay(icol,ilay))
            if(idx_minor_scaling(imnr) > 0) then  ! there is a second gas that affects this gas's absorption
              ! scale by density of special gas
              if (scale_by_complement(imnr)) then ! scale by densities of all gases but the special one
                scaling = scaling * (1._wp - vmr(idx_minor_scaling(imnr)) / (1._wp+vmr(idx_h2o)) )
              else
                scaling = scaling *          vmr(idx_minor_scaling(imnr)) / (1._wp+vmr(idx_h2o))
              endif
            endif
          endif
          !
          ! Interpolation of absorption coefficient and calculation of optical depth
          !
          ! Which gpoint range does this minor gas affect?
          iml = minor_limits_gpt(1,imnr)
          imu = minor_limits_gpt(2,imnr)
          ! What is the starting point in the stored array of minor absorption coefficients?
          minor_start = kminor_start(imnr)
          do igpt = iml,imu
            tau_minor = 0._wp
            iflav = gpt_flv(igpt) ! eta interpolation depends on flavor
            minor_loc = minor_start + (igpt - iml) ! add offset to starting point
            kminor_loc = &
              interpolate2D(fminor(:,:,iflav,icol,ilay), &
                            kminor, &
                            minor_loc, jeta(:,iflav,icol,ilay), jtemp(icol,ilay))
              tau_minor = kminor_loc * scaling
            tau(igpt,ilay,icol) = tau(igpt,ilay,icol) + tau_minor
          enddo
        enddo
      enddo
    enddo
  end subroutine gas_optical_depths_minor
  ! ----------------------------------------------------------
  !
  ! compute Rayleigh scattering optical depths
  !
  subroutine compute_tau_rayleigh(ncol,nlay,ngpt,ngas,nflav, &
                                  gpoint_flavor,krayl,       &
                                  idx_h2o, col_dry,col_gas,  &
                                  fminor,jeta,tropo,jtemp,   &
                                  tau_rayleigh)
    integer,                                  intent(in ) :: ncol,nlay,ngpt,ngas,nflav
    integer,  dimension(:,:),                 intent(in ) :: gpoint_flavor
    real(wp), dimension(:,:,:,:),             intent(in ) :: krayl
    integer,                                  intent(in ) :: idx_h2o
    real(wp), dimension(ncol,nlay),           intent(in ) :: col_dry
    real(wp), dimension(ncol,nlay,0:ngas),    intent(in ) :: col_gas
    real(wp), dimension(2,2,nflav,ncol,nlay), intent(in ) :: fminor
    integer,  dimension(2,  nflav,ncol,nlay), intent(in ) :: jeta
    logical,  dimension(ncol,nlay),           intent(in ) :: tropo
    integer,  dimension(ncol,nlay),           intent(in ) :: jtemp
    ! outputs
    real(wp), dimension(ngpt,nlay,ncol),      intent(out) :: tau_rayleigh
    ! -----------------
    ! local variables
    real(wp) :: k ! rayleigh scattering coefficient
    integer  :: icol, ilay, iflav, igpt
    integer  :: itropo
    ! -----------------
    do ilay = 1, nlay
      do icol = 1, ncol
        itropo = merge(1,2,tropo(icol,ilay)) ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        do igpt = 1, ngpt
          iflav = gpoint_flavor(itropo, igpt)
          k = interpolate2D(fminor(:,:,iflav,icol,ilay), &
                            krayl(:,:,:,itropo),      &
                            igpt, jeta(:,iflav,icol,ilay), jtemp(icol,ilay))
          tau_rayleigh(igpt,ilay,icol) =  k * (col_gas(icol,ilay,idx_h2o)+col_dry(icol,ilay))
        end do ! igpt
      end do
    end do ! ilay
  end subroutine compute_tau_rayleigh

  ! ----------------------------------------------------------
  subroutine compute_Planck_source(ncol, nlay, ngpt, nbnd, nflav,  &
                    tlay, tlev, tsfc, sfc_lay,            &
                    fmajor, jeta, tropo, jtemp, jpress,   &
                    gpoint_bands, pfracin, temp_ref_min, totplnk_delta, totplnk, gpoint_flavor, &
                    sfc_src, lay_src, lev_src_inc, lev_src_dec)
    integer,                                    intent(in) :: ncol, nlay, ngpt, nbnd, nflav, sfc_lay
    real(wp), dimension(ncol,nlay  ),           intent(in) :: tlay
    real(wp), dimension(ncol,nlay+1),           intent(in) :: tlev
    real(wp), dimension(ncol       ),           intent(in) :: tsfc
    ! Interpolation variables
    real(wp), dimension(2,2,2,nflav,ncol,nlay), intent(in) :: fmajor
    integer,  dimension(2,    nflav,ncol,nlay), intent(in) :: jeta
    logical,  dimension(            ncol,nlay), intent(in) :: tropo
    integer,  dimension(            ncol,nlay), intent(in) :: jtemp, jpress
    ! Table-specific
    integer, dimension(ngpt),     intent(in) :: gpoint_bands  ! band number for each g-point
    real(wp),                     intent(in) :: temp_ref_min, totplnk_delta
    real(wp), dimension(:,:,:,:), intent(in) :: pfracin       ! change to provide size
    real(wp), dimension(:,:),     intent(in) :: totplnk       ! change to provide size
    integer,  dimension(:,:),     intent(in) :: gpoint_flavor ! change to provide size

    real(wp), dimension(ncol,     ngpt), intent(out) :: sfc_src
    real(wp), dimension(ncol,nlay,ngpt), intent(out) :: lay_src
    real(wp), dimension(ncol,nlay,ngpt), intent(out) :: lev_src_inc, lev_src_dec
    ! -----------------
    ! local
    integer  :: ilay, icol, igpt, itropo, iflav
    real(wp) :: pfrac          (ngpt,nlay,  ncol)
    real(wp) :: planck_function(nbnd,nlay+1,ncol)
    ! -----------------

    ! Calculation of fraction of band's Planck irradiance associated with each g-point
    do icol = 1, ncol
      do ilay = 1, nlay
        ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        itropo = merge(1,2,tropo(icol,ilay))
        do igpt = 1, ngpt
          iflav = gpoint_flavor(itropo, igpt) !eta interpolation depends on band's flavor
          pfrac(igpt,ilay,icol) = &
            ! interpolation in temperature, pressure, and eta
            interpolate3D((/1._wp,1._wp/), fmajor(:,:,:,iflav,icol,ilay), pfracin, &
                          igpt, jeta(:,iflav,icol,ilay), jtemp(icol,ilay),jpress(icol,ilay)+itropo)
        end do ! igpt
      end do ! col
    end do

    !
    ! Planck function by band for the surface
    ! Compute surface source irradiance for g-point, equals band irradiance x fraction for g-point
    !
    do icol = 1, ncol
      planck_function(1:nbnd,1,icol) = interpolate1D(tsfc(icol), temp_ref_min, totplnk_delta, totplnk)
      !
      ! Map to g-points
      !
      do igpt = 1, ngpt
        sfc_src(icol,igpt) = pfrac(igpt,sfc_lay,icol) * planck_function(gpoint_bands(igpt), 1, icol)
      end do
    end do ! icol

    do icol = 1, ncol
      do ilay = 1, nlay
        ! Compute layer source irradiance for g-point, equals band irradiance x fraction for g-point
        planck_function(1:nbnd,ilay,icol) = interpolate1D(tlay(icol,ilay), temp_ref_min, totplnk_delta, totplnk)
        !
        ! Map to g-points
        !
        do igpt = 1, ngpt
          lay_src(icol,ilay,igpt) = pfrac(igpt,ilay,icol) * planck_function(gpoint_bands(igpt),ilay,icol)
        end do
      end do ! ilay
    end do ! icol

    ! compute level source irradiances for each g-point, one each for upward and downward paths
    do icol = 1, ncol
      planck_function(1:nbnd,       1,icol) = interpolate1D(tlev(icol,     1), temp_ref_min, totplnk_delta, totplnk)
      do ilay = 1, nlay
        planck_function(1:nbnd,ilay+1,icol) = interpolate1D(tlev(icol,ilay+1), temp_ref_min, totplnk_delta, totplnk)
        !
        ! Map to g-points
        !
        do igpt = 1, ngpt
          lev_src_inc(icol,ilay,igpt) = pfrac(igpt,ilay,icol) * planck_function(gpoint_bands(igpt),ilay+1,icol)
          lev_src_dec(icol,ilay,igpt) = pfrac(igpt,ilay,icol) * planck_function(gpoint_bands(igpt),ilay,  icol)
        end do
      end do ! ilay
    end do ! icol

  end subroutine compute_Planck_source
  ! ----------------------------------------------------------
  !
  ! One dimensional interpolation -- return all values along second table dimension
  !
  pure function interpolate1D(val, offset, delta, table) result(res)
    ! input
    real(wp), intent(in) :: val,    & ! axis value at which to evaluate table
                            offset, & ! minimum of table axis
                            delta     ! step size of table axis
    real(wp), dimension(:,:), &
              intent(in) :: table ! dimensions (axis, values)
    ! output
    real(wp), dimension(size(table,dim=2)) :: res

    ! local
    real(wp) :: val0 ! fraction index adjusted by offset and delta
    integer :: index ! index term
    real(wp) :: frac ! fractional term
    ! -------------------------------------
    val0 = (val - offset) / delta
    frac = val0 - int(val0) ! get fractional part
    index = min(size(table,dim=1)-1, max(1, int(val0)+1)) ! limit the index range
    res(:) = table(index,:) + frac * (table(index+1,:) - table(index,:))
  end function interpolate1D
 ! ------------
 !   This function returns a single value from a subset (in gpoint) of the k table
 !
  pure function interpolate2D(fminor, k, igpt, jeta, jtemp) result(res)
    real(wp), dimension(2,2), intent(in) :: fminor ! interpolation fractions for minor species
                                       ! index(1) : reference eta level (temperature dependent)
                                       ! index(2) : reference temperature level
    real(wp), dimension(:,:,:), intent(in) :: k ! (g-point, eta, temp)
    integer,                    intent(in) :: igpt, jtemp ! interpolation index for temperature
    integer, dimension(2),      intent(in) :: jeta ! interpolation index for binary species parameter (eta)
    real(wp)                             :: res ! the result

    res =  &
      fminor(1,1) * k(igpt, jeta(1)  , jtemp  ) + &
      fminor(2,1) * k(igpt, jeta(1)+1, jtemp  ) + &
      fminor(1,2) * k(igpt, jeta(2)  , jtemp+1) + &
      fminor(2,2) * k(igpt, jeta(2)+1, jtemp+1)
  end function interpolate2D

  ! ----------------------------------------------------------
  ! interpolation in temperature, pressure, and eta
  pure function interpolate3D(scaling, fmajor, k, igpt, jeta, jtemp, jpress) result(res)
    real(wp), dimension(2),     intent(in) :: scaling
    real(wp), dimension(2,2,2), intent(in) :: fmajor ! interpolation fractions for major species
                                                     ! index(1) : reference eta level (temperature dependent)
                                                     ! index(2) : reference pressure level
                                                     ! index(3) : reference temperature level
    real(wp), dimension(:,:,:,:),intent(in) :: k ! (gpt, eta,temp,press)
    integer,                     intent(in) :: igpt
    integer, dimension(2),       intent(in) :: jeta ! interpolation index for binary species parameter (eta)
    integer,                     intent(in) :: jtemp ! interpolation index for temperature
    integer,                     intent(in) :: jpress ! interpolation index for pressure
    real(wp)                                :: res ! the result
    ! each code block is for a different reference temperature
    res =  &
      scaling(1) * &
      ( fmajor(1,1,1) * k(igpt, jeta(1)  , jpress-1, jtemp  ) + &
        fmajor(2,1,1) * k(igpt, jeta(1)+1, jpress-1, jtemp  ) + &
        fmajor(1,2,1) * k(igpt, jeta(1)  , jpress  , jtemp  ) + &
        fmajor(2,2,1) * k(igpt, jeta(1)+1, jpress  , jtemp  ) ) + &
      scaling(2) * &
      ( fmajor(1,1,2) * k(igpt, jeta(2)  , jpress-1, jtemp+1) + &
        fmajor(2,1,2) * k(igpt, jeta(2)+1, jpress-1, jtemp+1) + &
        fmajor(1,2,2) * k(igpt, jeta(2)  , jpress  , jtemp+1) + &
        fmajor(2,2,2) * k(igpt, jeta(2)+1, jpress  , jtemp+1) )
  end function interpolate3D

  ! ----------------------------------------------------------
  ! Compute interpolation coefficients
  ! for calculations of major optical depths, minor optical depths, Rayleigh,
  ! and Planck fractions
  subroutine interpolation(ncol,nlay,ngas,nflav,neta, &
    flavor,press_ref_log,temp_ref,press_ref_log_delta,temp_ref_min,temp_ref_delta,press_ref_trop_log,vmr_ref,nlay_ref, &
    play,tlay,col_gas, &
    jtemp,fmajor,fminor,col_mix,tropo,itropo_lower,itropo_upper,jeta,jpress)
    ! input dimensions
    integer, intent(in) :: ncol,nlay,ngas,nflav,neta

    ! inputs from object
    integer,  dimension(:,:),    intent(in) :: flavor
    real(wp), dimension(:),      intent(in) :: press_ref_log
    real(wp), dimension(:),      intent(in) :: temp_ref
    real(wp),                    intent(in) :: press_ref_log_delta, &
                                               temp_ref_min, temp_ref_delta, &
                                               press_ref_trop_log
    real(wp), dimension(:,0:,:), intent(in) :: vmr_ref
    integer,                     intent(in) :: nlay_ref

    ! inputs from profile or parent function
    real(wp), dimension(ncol,nlay),        intent(in) :: play, tlay
    real(wp), dimension(ncol,nlay,0:ngas), intent(in) :: col_gas

    ! outputs
    integer,  dimension(ncol,nlay), intent(out) :: jtemp, jpress
    logical,  dimension(ncol,nlay), intent(out) :: tropo
    integer,  dimension(ncol, 2),   intent(out) :: itropo_lower,itropo_upper
    integer,  dimension(2,    nflav,ncol,nlay), intent(out) :: jeta
    real(wp), dimension(2,    nflav,ncol,nlay), intent(out) :: col_mix
    real(wp), dimension(2,2,2,nflav,ncol,nlay), intent(out) :: fmajor
    real(wp), dimension(2,2,  nflav,ncol,nlay), intent(out) :: fminor
    ! -----------------
    ! local
    real(wp), dimension(ncol,nlay) :: ftemp, fpress ! interpolation fraction for temperature, pressure
    real(wp) :: locpress ! needed to find location in pressure grid
    real(wp) :: ratio_eta_half ! ratio of vmrs of major species that defines eta=0.5
                               ! for given flavor and reference temperature level
    real(wp) :: eta, feta   ! binary_species_parameter, interpolation variable for eta
    real(wp) :: loceta ! needed to find location in eta grid
    real(wp) :: ftemp_term
    ! -----------------
    ! local indexes
    integer :: icol, ilay, iflav, igases(2), itropo, itemp
    integer, dimension(ncol) :: itropo_last

    do ilay = 1, nlay
      do icol = 1, ncol
        ! index and factor for temperature interpolation
        jtemp(icol,ilay) = int((tlay(icol,ilay) - (temp_ref_min - temp_ref_delta)) / temp_ref_delta)
        jtemp(icol,ilay) = min(size(temp_ref) - 1, max(1, jtemp(icol,ilay))) ! limit the index range
        ftemp(icol,ilay) = (tlay(icol,ilay) - temp_ref(jtemp(icol,ilay))) / temp_ref_delta

        ! index and factor for pressure interpolation
        locpress = 1._wp + (log(play(icol,ilay)) - press_ref_log(1)) / press_ref_log_delta
        jpress(icol,ilay) = min(nlay_ref-2, max(1, int(locpress)))
        fpress(icol,ilay) = locpress - float(jpress(icol,ilay))

        ! determine if in lower or upper part of atmosphere
        tropo(icol,ilay) = log(play(icol,ilay)) > press_ref_trop_log

        ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        itropo = merge(1,2,tropo(icol,ilay))
        if (ilay .eq. 1) then
          itropo_last(icol) = itropo
        endif
        if (itropo_last(icol) .ne. itropo) then
          if (itropo .eq. 1) then ! layers go from TOA to surface
            itropo_lower(icol,1) = ilay
            itropo_lower(icol,2) = nlay
            itropo_upper(icol,1) = 1
            itropo_upper(icol,2) = ilay - 1
            itropo_last(icol) = itropo
          else ! layers go from surface to TOA
            itropo_lower(icol,1) = 1
            itropo_lower(icol,2) = ilay - 1
            itropo_upper(icol,1) = ilay
            itropo_upper(icol,2) = nlay
            itropo_last(icol) = itropo
          endif
        endif
        ! loop over implemented combinations of major species
        do iflav = 1, nflav
          igases(:) = flavor(:,iflav)
          do itemp = 1, 2
            ! compute interpolation fractions needed for lower, then upper reference temperature level
            ! compute binary species parameter (eta) for flavor and temperature and
            !  associated interpolation index and factors
            ratio_eta_half = vmr_ref(itropo,igases(1),(jtemp(icol,ilay)+itemp-1)) / &
                             vmr_ref(itropo,igases(2),(jtemp(icol,ilay)+itemp-1))
            col_mix(itemp,iflav,icol,ilay) = col_gas(icol,ilay,igases(1)) + ratio_eta_half * col_gas(icol,ilay,igases(2))
            eta = merge(col_gas(icol,ilay,igases(1)) / col_mix(itemp,iflav,icol,ilay), 0.5_wp, &
                        col_mix(itemp,iflav,icol,ilay) > 2._wp * tiny(col_mix))
            loceta = eta * float(neta-1)
            jeta(itemp,iflav,icol,ilay) = min(int(loceta)+1, neta-1)
            feta = mod(loceta, 1.0_wp)
            ! compute interpolation fractions needed for minor species
            ! ftemp_term = (1._wp-ftemp(icol,ilay)) for itemp = 1, ftemp(icol,ilay) for itemp=1
            ftemp_term = (real(2-itemp, wp) + real(2*itemp-3, wp) * ftemp(icol,ilay))
            fminor(1,itemp,iflav,icol,ilay) = (1._wp-feta) * ftemp_term
            fminor(2,itemp,iflav,icol,ilay) =        feta  * ftemp_term
            ! compute interpolation fractions needed for major species
            fmajor(1,1,itemp,iflav,icol,ilay) = (1._wp-fpress(icol,ilay)) * fminor(1,itemp,iflav,icol,ilay)
            fmajor(2,1,itemp,iflav,icol,ilay) = (1._wp-fpress(icol,ilay)) * fminor(2,itemp,iflav,icol,ilay)
            fmajor(1,2,itemp,iflav,icol,ilay) =        fpress(icol,ilay)  * fminor(1,itemp,iflav,icol,ilay)
            fmajor(2,2,itemp,iflav,icol,ilay) =        fpress(icol,ilay)  * fminor(2,itemp,iflav,icol,ilay)
          end do ! reference temperatures
        end do ! iflav
      end do ! icol,ilay
    end do

  end subroutine interpolation
  ! ----------------------------------------------------------
  !
  ! Combine absoprtion and Rayleigh optical depths for total tau, ssa, g
  !
  pure subroutine combine_and_reorder_2str(ncol, nlay, ngpt, tau_abs, tau_rayleigh, tau, ssa, g) &
      bind(C, name="combine_and_reorder_2str")
    integer,                             intent(in) :: ncol, nlay, ngpt
    real(wp), dimension(ngpt,nlay,ncol), intent(in   ) :: tau_abs, tau_rayleigh
    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau, ssa, g ! inout because components are allocated
    ! -----------------------
    integer  :: icol, ilay, igpt
    real(wp) :: t
    ! -----------------------
    do icol = 1, ncol
      do ilay = 1, nlay
        do igpt = 1, ngpt
           t = tau_abs(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol)
           tau(icol,ilay,igpt) = t
           g  (icol,ilay,igpt) = 0._wp
           if(t > 2._wp * tiny(t)) then
             ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / t
           else
             ssa(icol,ilay,igpt) = 0._wp
           end if
        end do
      end do
    end do
  end subroutine combine_and_reorder_2str
  ! ----------------------------------------------------------
  !
  ! Combine absoprtion and Rayleigh optical depths for total tau, ssa, p
  !   using Rayleigh scattering phase function
  !
  pure subroutine combine_and_reorder_nstr(ncol, nlay, ngpt, nmom, tau_abs, tau_rayleigh, tau, ssa, p) &
      bind(C, name="combine_and_reorder_nstr")
    integer, intent(in) :: ncol, nlay, ngpt, nmom
    real(wp), dimension(ngpt,nlay,ncol), intent(in ) :: tau_abs, tau_rayleigh
    real(wp), dimension(ncol,nlay,ngpt), intent(inout) :: tau, ssa
    real(wp), dimension(ncol,nlay,ngpt,nmom), &
                                         intent(inout) :: p
    ! -----------------------
    integer :: icol, ilay, igpt, imom
    real(wp) :: t
    ! -----------------------
    do icol = 1, ncol
      do ilay = 1, nlay
        do igpt = 1, ngpt
          t = tau_abs(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol)
          tau(icol,ilay,igpt) = t
          if(t > 2._wp * tiny(t)) then
            ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / t
          else
            ssa(icol,ilay,igpt) = 0._wp
          end if
          p(1,icol,ilay,igpt) = 0.0_wp
          p(2,icol,ilay,igpt) = 0.1_wp
          do imom = 1, nmom
            p(imom,icol,ilay,igpt) = 0.0_wp
          end do
        end do
      end do
    end do
  end subroutine combine_and_reorder_nstr
  ! ----------------------------------------------------------
  pure subroutine zero_array_3D(ni, nj, nk, array) bind(C, name="zero_array_3D")
    integer, intent(in) :: ni, nj, nk
    real(wp), dimension(ni, nj, nk), intent(out) :: array
    ! -----------------------
    integer :: i,j,k
    ! -----------------------
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          array(i,j,k) = 0.0_wp
        end do
      end do
    end do

  end subroutine zero_array_3D
  ! ----------------------------------------------------------
  pure subroutine zero_array_4D(ni, nj, nk, nl, array) bind(C, name="zero_array_4D")
    integer, intent(in) :: ni, nj, nk, nl
    real(wp), dimension(ni, nj, nk, nl), intent(out) :: array
    ! -----------------------
    integer :: i,j,k,l
    ! -----------------------
    do l = 1, nl
      do k = 1, nk
        do j = 1, nj
          do i = 1, ni
            array(i,j,k,l) = 0.0_wp
          end do
        end do
      end do
    end do

  end subroutine zero_array_4D
  ! ----------------------------------------------------------
end module mo_gas_optics_kernels
