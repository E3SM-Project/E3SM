module cam_rad_parameterizations
  use shr_orb_mod
  !
  ! CAM 3.0 parameterizations related to radiation
  !
  implicit none
  private
  
  real, parameter :: tmelt = 273.16  ! specify melting temperature
  !       Tabulated values of re(T) in the temperature interval
  !       180 K -- 274 K; hexagonal columns assumed
  !
  real, parameter :: iceSizeTableMinTemp = 180. 
  real, dimension(95), parameter :: &
    retab = (/ 5.92779, 6.26422, 6.61973, 6.99539, 7.39234,	&
        7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,	&
        10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319, &
        15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,	&
        20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,	&
        27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943, &
        31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601, &
        34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,	&
        38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,	&
        42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,	&
        50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,	&
        65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,	&
        93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424, &
        124.954, 130.630, 136.457, 142.446, 148.608, 154.956,	&
        161.503, 168.262, 175.248, 182.473, 189.952, 197.699,	&
        205.728, 214.055, 222.694, 231.661, 240.971, 250.639 /)


  public :: computeRe_Liquid, computeRe_Ice, albedo

contains
  !-----------------------------------------------------------------------
  elemental real function computeRe_Liquid(temperature, landfrac, icefrac, snowh) &
     result(rel)
    real,           intent(in) :: temperature, landfrac
    real, optional, intent(in) :: icefrac, snowh  ! Snow depth over land, water equivalent (m)

    real, parameter ::  rliqland  =  8.0, & ! liquid drop size if over land
                        rliqocean = 14.0, & ! liquid drop size if over ocean
                        rliqice   = 14.0    ! liquid drop size if over sea ice
 
    ! jrm Reworked effective radius algorithm
    ! Start with temperature-dependent value appropriate for continental air
    rel = rliqland + (rliqocean - rliqland) * min(1.0, max(0.0, (tmelt - temperature) * 0.05))
    
    if(present(snowh)) & ! Modify for snow depth over land
      rel = rel + (rliqocean - rel) * min(1.0, max(0.0, snowh*10.))
      
    ! Ramp between polluted value over land to clean value over ocean.
    rel = rel + (rliqocean-rel) * min(1.0, max(0.0, 1.0 - landfrac))
    
    if(present(icefrac)) & ! Ramp between the resultant value and a sea ice value in the presence of ice.
      rel = rel + (rliqice-rel) * min(1.0, max(0.0, icefrac))

  end function computeRe_Liquid
  !-----------------------------------------------------------------------
  elemental real function computeRe_Ice(temperature) result(rei) 
    real,           intent(in) :: temperature
  
    real    :: fraction
    integer :: index
    !
    !
    index = int(temperature - (iceSizeTableMinTemp - 1.))
    index = min(max(index, 1), 94)
    fraction = temperature - int(temperature)
    rei = retab(index) * (1. - fraction) + retab(index+1) * fraction
  end function computeRe_Ice
  !-----------------------------------------------------------------------
  subroutine albedo(ocean, coszrs, asdir, aldir, asdif, aldif)
    !-----------------------------------------------------------------------
    ! Computes surface albedos over ocean 
    ! and the surface (added by Marat Khairoutdinov)
  
    ! Two spectral surface albedos for direct (dir) and diffuse (dif)
    ! incident radiation are calculated. The spectral intervals are:
    !   s (shortwave)  = 0.2-0.7 micro-meters
    !   l (longwave)   = 0.7-5.0 micro-meters
    !
    ! Uses knowledge of surface type to specify albedo, as follows:
    !
    ! Ocean           Uses solar zenith angle to compute albedo for direct
    !                 radiation; diffuse radiation values constant; albedo
    !                 independent of spectral interval and other physical
    !                 factors such as ocean surface wind speed.
    !
    ! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
    ! Approximation for Solar Radiation in the NCAR Community Climate Model,
    ! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
    
    logical,            intent( in) :: ocean
    real, dimension(:), intent( in) :: coszrs   ! Cosine of solar zenith angle
    real, dimension(:), intent(out) :: asdir, & ! Srf alb for direct rad   0.2-0.7 micro-ms
                                       aldir, & ! Srf alb for direct rad   0.7-5.0 micro-ms
                                       asdif, & ! Srf alb for diffuse rad  0.2-0.7 micro-ms
                                       aldif    ! Srf alb for diffuse rad  0.7-5.0 micro-ms
    real, parameter :: adif = 0.06
    !-----------------------------------------------------------------------
    if (ocean) then
      !
      !
      ! Ice-free ocean albedos function of solar zenith angle only, and
      ! independent of spectral interval:
      !
      where(coszrs <= 0) 
        aldir(:) = 0; asdir(:) = 0; aldif(:) = 0; asdif(:) = 0
      elsewhere
        aldir(:) = ( .026 / (coszrs(:)**1.7 + .065)) + &
             (.15*(coszrs(:) - 0.10) * (coszrs(:) - 0.50) * (coszrs(:) - 1.00) )
        asdir(:) = aldir(:)
        aldif(:) = adif
        asdif(:) = adif
      end where
    else ! land
      where(coszrs <= 0) 
        aldir(:) = 0; asdir(:) = 0; aldif(:) = 0; asdif(:) = 0
      elsewhere
        ! Albedos for land type I (Briegleb)
        asdir(:) = 1.4 * 0.06 / ( 1. + 0.8 * coszrs(:))
        asdif(:) = 1.2 * 0.06
        aldir(:) = 1.4 * 0.24 / ( 1. + 0.8 * coszrs(:))
        aldif(:) = 1.2 * 0.24
      end where
    endif
  end subroutine albedo

end module cam_rad_parameterizations