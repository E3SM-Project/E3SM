!=======================================================================
! A New Algorithm for Calculation of Consine Solar Zenith Angle to average 
! over a radiation time step
! Author: Linjiong Zhou
! E-mail: linjiongzhou@hotmail.com
! Date  : 2015.02.22
! Ref.  : Zhou et al., GRL, 2015

! Ported to ACME by BSINGH (04/21/2015) 
! BSINGH Modification for porting:
!  1. Modified substantially to avoid reading atm_in namelist again here
!  2. Modified names of the subroutines to better represent their function
!  3. Added avg_slr_insol_init subroutine to initialize variables needed from radiation
!  4. Removed unnecessary variables and added 'only' to 'use' statements
!  5. Used shr_const_pi instead of computing pi in this subroutine
!  6. Rest of the code is ported as it is

! Usage:
! put this module file into cesm?_?_?/models/csm_share/shr
! replace "shr_orb_cosz = ..." in shr_orb_mod.F90 by:
! "call avg_slr_insol(lat, lon, declin, jday, shr_orb_cosz)"
!=======================================================================

module avg_slr_insol_mod
  
  use shr_kind_mod,   only: shr_kind_r8
  use shr_infnan_mod, only: shr_infnan_nan, shr_infnan_inf, assignment(=)
  
  implicit none
  private
  save
  
  real(shr_kind_r8) :: rdt  ! radiation time step (second)
  
  
  public :: avg_slr_insol
  public :: avg_slr_insol_init !initializes variables required from atm-radiation
  
contains
  !=======================================================================
  
  subroutine avg_slr_insol_init(dtime, iradsw)
    !-----------------------------------------------------------------------
    !BSINGH - This subroutine initializes radiation time step by using 
    !parameters from atm_in namelist (cam_inparm) file (read in 
    !runtime_opts.F90)
    !
    ! This subroutine is called from:
    !    shr_orb_mod_init [shr_orb_mod.F90]  
    !-----------------------------------------------------------------------
    
    !Args
    integer, intent(in) :: dtime           ! timestep size (s)
    integer, intent(in) :: iradsw          ! freq. of shortwave radiation calc in time steps (positive)
                                           ! or hours (negative).
    rdt = shr_infnan_inf
    !compute 'rdt' [a global module variable]
    if (iradsw .lt. 0) then
       rdt = - real(iradsw,shr_kind_r8) * real(3600.0,shr_kind_r8)
    else
       rdt =   real(iradsw,shr_kind_r8) * real(dtime,shr_kind_r8)
    end if

  end subroutine avg_slr_insol_init
  
  
  subroutine avg_slr_insol(lat, lon, dec, jday, cosz)

    use shr_const_mod, only: pi => shr_const_pi

    implicit none
    
    !-----------------------------------------------------------------------
    ! In/Out Arguements
    
    real(shr_kind_r8), intent(in) :: lat   ! latitude (radian)
    real(shr_kind_r8), intent(in) :: lon   ! longitude (radian)
    real(shr_kind_r8), intent(in) :: dec   ! solar declination (radian)
    real(shr_kind_r8), intent(in) :: jday  ! Julian calendar day (1.xx to 365.xx)
    
    real(shr_kind_r8), intent(out) :: cosz ! cosine solar zenith angle (1)
    
    !-----------------------------------------------------------------------
    ! Local Arguments

    real(shr_kind_r8) :: aa, bb
    real(shr_kind_r8) :: del, phi
    real(shr_kind_r8) :: cos_h, h
    real(shr_kind_r8) :: t1, t2, dt
    real(shr_kind_r8) :: tt1, tt2, tt3, tt4

    !-----------------------------------------------------------------------
    ! Compute Half-day Length   
    
    ! adjust latitude so that its tangent will be defined
    del = lat
    if (lat .eq.   pi / 2.0_shr_kind_r8) then
       del = lat - 1.0e-05_shr_kind_r8
    end if
    if (lat .eq. - pi / 2.0_shr_kind_r8) then
       del = lat + 1.0e-05_shr_kind_r8
    end if
    
    ! adjust declination so that its tangent will be defined
    phi = dec
    if (dec .eq.   pi / 2.0_shr_kind_r8) then
       phi = dec - 1.0e-05_shr_kind_r8
    end if
    if (dec .eq. - pi / 2.0_shr_kind_r8) then
       phi = dec + 1.0e-05_shr_kind_r8
    end if
    
    ! define the cosine of the half-day length
    ! adjust for cases of all daylight or all night
    cos_h = - tan(del) * tan(phi)
    if (cos_h .le. - 1.0_shr_kind_r8) then
       h = pi
    end if
    if (cos_h .ge.   1.0_shr_kind_r8) then
       h = 0.0_shr_kind_r8
    end if
    if (cos_h .gt. - 1.0_shr_kind_r8 .and. cos_h .lt. 1.0_shr_kind_r8) then
       h = acos(cos_h)
    end if
    
    !-----------------------------------------------------------------------
    ! Define Local Time t and t + dt
    
    ! adjust t to be between -pi and pi
    t1 = (jday - int(jday)) * 2.0_shr_kind_r8 * pi + lon - pi
    if (t1 .ge.   pi) t1 = t1 - 2.0_shr_kind_r8 * pi
    if (t1 .lt. - pi) t1 = t1 + 2.0_shr_kind_r8 * pi
    
    dt = rdt / 86400.0_shr_kind_r8 * 2.0_shr_kind_r8 * pi
    t2 = t1 + dt
    
    !-----------------------------------------------------------------------
    ! Comput Cosine Solar Zenith angle
    
    ! define terms needed in the cosine zenith angle equation
    aa = sin(lat) * sin(dec)
    bb = cos(lat) * cos(dec)
    
    ! define the hour angle
    ! force it to be between -h and h
    ! consider the situation when the night period is too short
    if (t2 .ge. pi .and. t1 .le. pi .and. pi - h .le. dt) then
       tt2 = h
       tt1 = min(max(t1,                      - h),                        h)
       tt4 = min(max(t2, 2.0_shr_kind_r8 * pi - h), 2.0_shr_kind_r8 * pi + h)
       tt3 = 2.0_shr_kind_r8 * pi - h
    else if (t2 .ge. - pi .and. t1 .le. - pi .and. pi - h .le. dt) then
       tt2 = - 2.0_shr_kind_r8 * pi + h
       tt1 = min(max(t1, - 2.0_shr_kind_r8 * pi - h), - 2.0_shr_kind_r8 * pi + h)
       tt4 = min(max(t2,                        - h),                          h)
       tt3 = - h
    else
       if (t2 .gt. pi) then
          tt2 = min(max(t2 - 2.0_shr_kind_r8 * pi, - h), h)
       else if (t2 .lt. - pi) then
          tt2 = min(max(t2 + 2.0_shr_kind_r8 * pi, - h), h)
       else
          tt2 = min(max(t2                       , - h), h)
       end if
       if (t1 .gt. pi) then
          tt1 = min(max(t1 - 2.0_shr_kind_r8 * pi, - h), h)
       else if (t1 .lt. - pi) then
          tt1 = min(max(t1 + 2.0_shr_kind_r8 * pi, - h), h)
       else
          tt1 = min(max(t1                       , - h), h)
       end if
       tt4 = 0.0_shr_kind_r8
       tt3 = 0.0_shr_kind_r8
    end if
    
    ! perform a time integration to obtain cosz if desired
    ! output is valid over the period from t to t + dt
    if (tt2 .gt. tt1 .or. tt4 .gt. tt3) then
       cosz = (aa * (tt2 - tt1) + bb * (sin(tt2) - sin(tt1))) / dt + &
            (aa * (tt4 - tt3) + bb * (sin(tt4) - sin(tt3))) / dt
    else
       cosz = 0.0_shr_kind_r8
    end if
    
  end subroutine avg_slr_insol
  
end module avg_slr_insol_mod
