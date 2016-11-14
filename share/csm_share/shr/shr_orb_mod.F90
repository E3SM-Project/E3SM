!===============================================================================
! SVN $Id: shr_orb_mod.F90 25434 2010-11-04 22:46:24Z tcraig $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/release_tags/cesm1_2_x_n02_share3_130715/shr/shr_orb_mod.F90 $
!===============================================================================

MODULE shr_orb_mod

   use shr_kind_mod, only: SHR_KIND_R8, SHR_KIND_IN
   use shr_sys_mod
   use shr_const_mod
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit

   IMPLICIT none

   !----------------------------------------------------------------------------
   ! PUBLIC: Interfaces and global data
   !----------------------------------------------------------------------------
   public :: shr_orb_cosz
   public :: shr_orb_params
   public :: shr_orb_decl
   public :: shr_orb_print

   real   (SHR_KIND_R8),public,parameter :: SHR_ORB_UNDEF_REAL = 1.e36_SHR_KIND_R8 ! undefined real
   integer(SHR_KIND_IN),public,parameter :: SHR_ORB_UNDEF_INT  = 2000000000        ! undefined int

   !----------------------------------------------------------------------------
   ! PRIVATE: by default everything else is private to this module
   !----------------------------------------------------------------------------
   private

   real   (SHR_KIND_R8),parameter :: pi                 = SHR_CONST_PI
   real   (SHR_KIND_R8),parameter :: SHR_ORB_ECCEN_MIN  =   0.0_SHR_KIND_R8 ! min value for eccen
   real   (SHR_KIND_R8),parameter :: SHR_ORB_ECCEN_MAX  =   0.1_SHR_KIND_R8 ! max value for eccen
   real   (SHR_KIND_R8),parameter :: SHR_ORB_OBLIQ_MIN  = -90.0_SHR_KIND_R8 ! min value for obliq
   real   (SHR_KIND_R8),parameter :: SHR_ORB_OBLIQ_MAX  = +90.0_SHR_KIND_R8 ! max value for obliq
   real   (SHR_KIND_R8),parameter :: SHR_ORB_MVELP_MIN  =   0.0_SHR_KIND_R8 ! min value for mvelp
   real   (SHR_KIND_R8),parameter :: SHR_ORB_MVELP_MAX  = 360.0_SHR_KIND_R8 ! max value for mvelp


!===============================================================================
CONTAINS
!===============================================================================

real(SHR_KIND_R8) pure FUNCTION shr_orb_cosz(jday,lat,lon,declin,dt_avg)

   !----------------------------------------------------------------------------
   !
   ! FUNCTION to return the cosine of the solar zenith angle.
   ! Assumes 365.0 days/year.
   !
   !--------------- Code History -----------------------------------------------
   !
   ! Original Author: Brian Kauffman
   ! Date:            Jan/98
   ! History:         adapted from statement FUNCTION in share/orb_cosz.h
   !
   !----------------------------------------------------------------------------

   real   (SHR_KIND_R8),intent(in) :: jday   ! Julian cal day (1.xx to 365.xx)
   real   (SHR_KIND_R8),intent(in) :: lat    ! Centered latitude (radians)
   real   (SHR_KIND_R8),intent(in) :: lon    ! Centered longitude (radians)
   real   (SHR_KIND_R8),intent(in) :: declin ! Solar declination (radians)
   real   (SHR_KIND_R8),intent(in), optional   :: dt_avg ! if present and set non-zero, then use in the
                                                         ! average cosz calculation
   logical :: use_dt_avg

   !----------------------------------------------------------------------------

   use_dt_avg = .false.
   if (present(dt_avg)) then
      if (dt_avg /= 0.0_shr_kind_r8) use_dt_avg = .true.
   end if


   ! If dt for the average cosz is specified, then call the shr_orb_avg_cosz
   if (use_dt_avg) then
      shr_orb_cosz =  shr_orb_avg_cosz(jday, lat, lon, declin, dt_avg)
   else
      shr_orb_cosz = sin(lat)*sin(declin) - &
      &              cos(lat)*cos(declin)*cos(jday*2.0_SHR_KIND_R8*pi + lon)
! RLJ old version above.  Restore version below more accurate version in separate PR 
!           cos(lat)*cos(declin) * &
!           cos((jday-floor(jday))*2.0_SHR_KIND_R8*pi + lon)
   end if

END FUNCTION shr_orb_cosz

!=======================================================================
! A New Algorithm for Calculation of Cosine Solar Zenith Angle
! Author: Linjiong Zhou
! E-mail: linjiongzhou@hotmail.com
! Date  : 2015.02.22
! Ref.  : Zhou et al., GRL, 2015
!=======================================================================

real (SHR_KIND_R8) pure function shr_orb_avg_cosz(jday, lat, lon, declin, dt_avg)

    use shr_const_mod, only : pi => shr_const_pi

    implicit none

!-----------------------------------------------------------------------
! In/Out Arguements

    real(SHR_KIND_R8), intent(in) :: jday   ! Julian calendar day (1.xx to 365.xx)
    real(SHR_KIND_R8), intent(in) :: lat    ! latitude (radian)
    real(SHR_KIND_R8), intent(in) :: lon    ! longitude (radian)
    real(SHR_KIND_R8), intent(in) :: declin ! solar declination (radian)
    real(SHR_KIND_R8), intent(in) :: dt_avg ! dt for averaged cosz calculation

!-----------------------------------------------------------------------
! Local Arguments

    real(SHR_KIND_R8),parameter :: piover2 = pi/2.0_SHR_KIND_R8
    real(SHR_KIND_R8),parameter :: twopi   = pi*2.0_SHR_KIND_R8

    real(SHR_KIND_R8) :: aa, bb
    real(SHR_KIND_R8) :: del, phi
    real(SHR_KIND_R8) :: cos_h, h
    real(SHR_KIND_R8) :: t1, t2, dt
    real(SHR_KIND_R8) :: tt1, tt2, tt3, tt4

!-----------------------------------------------------------------------
! Compute Half-day Length

    ! adjust latitude so that its tangent will be defined
    if (lat ==  piover2) then
        del = lat - 1.0e-05_SHR_KIND_R8
    else if (lat ==  -piover2) then
        del = lat + 1.0e-05_SHR_KIND_R8
    else
        del = lat
    end if

    ! adjust declination so that its tangent will be defined
    if (declin == piover2) then
        phi = declin - 1.0e-05_SHR_KIND_R8
    else if (declin == -piover2) then
        phi = declin + 1.0e-05_SHR_KIND_R8
    else
        phi = declin
    end if

    ! define the cosine of the half-day length
    ! adjust for cases of all daylight or all night
    cos_h = - tan(del) * tan(phi)
    if (cos_h <= -1.0_SHR_KIND_R8) then
        h = pi
    else if (cos_h >= 1.0_SHR_KIND_R8) then
        h = 0.0_SHR_KIND_R8
    else
        h = acos(cos_h)
    end if

!-----------------------------------------------------------------------
! Define Local Time t and t + dt

    ! adjust t to be between -pi and pi
    t1 = (jday - int(jday)) * twopi + lon - pi

    if (t1 >=  pi) then
       t1 = t1 - twopi
    else if (t1 < -pi) then
       t1 = t1 + twopi
    end if

    dt = dt_avg / 86400.0_SHR_KIND_R8 * twopi
    t2 = t1 + dt

!-----------------------------------------------------------------------
! Compute Cosine Solar Zenith angle

    ! define terms needed in the cosine zenith angle equation
    aa = sin(lat) * sin(declin)
    bb = cos(lat) * cos(declin)

    ! define the hour angle
    ! force it to be between -h and h
    ! consider the situation when the night period is too short
    if (t2 >= pi .and. t1 <= pi .and. pi - h <= dt) then
        tt2 = h
        tt1 = min(max(t1, -h)       ,         h)
        tt4 = min(max(t2, twopi - h), twopi + h)
        tt3 = twopi - h
    else if (t2 >= -pi .and. t1 <= -pi .and. pi - h <= dt) then
        tt2 = - twopi + h
        tt1 = min(max(t1, -twopi - h), -twopi + h)
        tt4 = min(max(t2, -h)        ,          h)
        tt3 = -h
    else
        if (t2 > pi) then
            tt2 = min(max(t2 - twopi, -h), h)
        else if (t2 < - pi) then
            tt2 = min(max(t2 + twopi, -h), h)
        else
            tt2 = min(max(t2 ,        -h), h)
        end if
        if (t1 > pi) then
            tt1 = min(max(t1 - twopi, -h), h)
        else if (t1 < - pi) then
            tt1 = min(max(t1 + twopi, -h), h)
        else
            tt1 = min(max(t1        , -h), h)
        end if
        tt4 = 0.0_SHR_KIND_R8
        tt3 = 0.0_SHR_KIND_R8
    end if

    ! perform a time integration to obtain cosz if desired
    ! output is valid over the period from t to t + dt
    if (tt2 > tt1 .or. tt4 > tt3) then
        shr_orb_avg_cosz = (aa * (tt2 - tt1) + bb * (sin(tt2) - sin(tt1))) / dt + &
               (aa * (tt4 - tt3) + bb * (sin(tt4) - sin(tt3))) / dt
    else
        shr_orb_avg_cosz = 0.0_SHR_KIND_R8
    end if

end function shr_orb_avg_cosz

!===============================================================================

SUBROUTINE shr_orb_params( iyear_AD , eccen  , obliq , mvelp     ,     &
           &               obliqr   , lambm0 , mvelpp, log_print )

!-------------------------------------------------------------------------------
!
! Calculate earths orbital parameters using Dave Threshers formula which
! came from Berger, Andre.  1978  "A Simple Algorithm to Compute Long-Term
! Variations of Daily Insolation".  Contribution 18, Institute of Astronomy
! and Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium
!
!------------------------------Code history-------------------------------------
!
! Original Author: Erik Kluzek
! Date:            Oct/97
!
!-------------------------------------------------------------------------------

   !----------------------------- Arguments ------------------------------------
   integer(SHR_KIND_IN),intent(in)    :: iyear_AD  ! Year to calculate orbit for
   real   (SHR_KIND_R8),intent(inout) :: eccen     ! orbital eccentricity
   real   (SHR_KIND_R8),intent(inout) :: obliq     ! obliquity in degrees
   real   (SHR_KIND_R8),intent(inout) :: mvelp     ! moving vernal equinox long
   real   (SHR_KIND_R8),intent(out)   :: obliqr    ! Earths obliquity in rad
   real   (SHR_KIND_R8),intent(out)   :: lambm0    ! Mean long of perihelion at
                                                   ! vernal equinox (radians)
   real   (SHR_KIND_R8),intent(out)   :: mvelpp    ! moving vernal equinox long
                                                   ! of perihelion plus pi (rad)
   logical             ,intent(in)    :: log_print ! Flags print of status/error

   !------------------------------ Parameters ----------------------------------
   integer(SHR_KIND_IN),parameter :: poblen =47 ! # of elements in series wrt obliquity
   integer(SHR_KIND_IN),parameter :: pecclen=19 ! # of elements in series wrt eccentricity
   integer(SHR_KIND_IN),parameter :: pmvelen=78 ! # of elements in series wrt vernal equinox
   real   (SHR_KIND_R8),parameter :: psecdeg = 1.0_SHR_KIND_R8/3600.0_SHR_KIND_R8 ! arc sec to deg conversion

   real   (SHR_KIND_R8) :: degrad = pi/180._SHR_KIND_R8   ! degree to radian conversion factor
   real   (SHR_KIND_R8) :: yb4_1950AD         ! number of years before 1950 AD

   character(len=*),parameter :: subname = '(shr_orb_params)'

   ! Cosine series data for computation of obliquity: amplitude (arc seconds),
   ! rate (arc seconds/year), phase (degrees).

   real   (SHR_KIND_R8), parameter :: obamp(poblen) =  & ! amplitudes for obliquity cos series
   &      (/   -2462.2214466_SHR_KIND_R8, -857.3232075_SHR_KIND_R8, -629.3231835_SHR_KIND_R8,   &
   &            -414.2804924_SHR_KIND_R8, -311.7632587_SHR_KIND_R8,  308.9408604_SHR_KIND_R8,   &
   &            -162.5533601_SHR_KIND_R8, -116.1077911_SHR_KIND_R8,  101.1189923_SHR_KIND_R8,   &
   &             -67.6856209_SHR_KIND_R8,   24.9079067_SHR_KIND_R8,   22.5811241_SHR_KIND_R8,   &
   &             -21.1648355_SHR_KIND_R8,  -15.6549876_SHR_KIND_R8,   15.3936813_SHR_KIND_R8,   &
   &              14.6660938_SHR_KIND_R8,  -11.7273029_SHR_KIND_R8,   10.2742696_SHR_KIND_R8,   &
   &               6.4914588_SHR_KIND_R8,    5.8539148_SHR_KIND_R8,   -5.4872205_SHR_KIND_R8,   &
   &              -5.4290191_SHR_KIND_R8,    5.1609570_SHR_KIND_R8,    5.0786314_SHR_KIND_R8,   &
   &              -4.0735782_SHR_KIND_R8,    3.7227167_SHR_KIND_R8,    3.3971932_SHR_KIND_R8,   &
   &              -2.8347004_SHR_KIND_R8,   -2.6550721_SHR_KIND_R8,   -2.5717867_SHR_KIND_R8,   &
   &              -2.4712188_SHR_KIND_R8,    2.4625410_SHR_KIND_R8,    2.2464112_SHR_KIND_R8,   &
   &              -2.0755511_SHR_KIND_R8,   -1.9713669_SHR_KIND_R8,   -1.8813061_SHR_KIND_R8,   &
   &              -1.8468785_SHR_KIND_R8,    1.8186742_SHR_KIND_R8,    1.7601888_SHR_KIND_R8,   &
   &              -1.5428851_SHR_KIND_R8,    1.4738838_SHR_KIND_R8,   -1.4593669_SHR_KIND_R8,   &
   &               1.4192259_SHR_KIND_R8,   -1.1818980_SHR_KIND_R8,    1.1756474_SHR_KIND_R8,   &
   &              -1.1316126_SHR_KIND_R8,    1.0896928_SHR_KIND_R8/)

   real   (SHR_KIND_R8), parameter :: obrate(poblen) = & ! rates for obliquity cosine series
   &        (/  31.609974_SHR_KIND_R8, 32.620504_SHR_KIND_R8, 24.172203_SHR_KIND_R8,   &
   &            31.983787_SHR_KIND_R8, 44.828336_SHR_KIND_R8, 30.973257_SHR_KIND_R8,   &
   &            43.668246_SHR_KIND_R8, 32.246691_SHR_KIND_R8, 30.599444_SHR_KIND_R8,   &
   &            42.681324_SHR_KIND_R8, 43.836462_SHR_KIND_R8, 47.439436_SHR_KIND_R8,   &
   &            63.219948_SHR_KIND_R8, 64.230478_SHR_KIND_R8,  1.010530_SHR_KIND_R8,   &
   &             7.437771_SHR_KIND_R8, 55.782177_SHR_KIND_R8,  0.373813_SHR_KIND_R8,   &
   &            13.218362_SHR_KIND_R8, 62.583231_SHR_KIND_R8, 63.593761_SHR_KIND_R8,   &
   &            76.438310_SHR_KIND_R8, 45.815258_SHR_KIND_R8,  8.448301_SHR_KIND_R8,   &
   &            56.792707_SHR_KIND_R8, 49.747842_SHR_KIND_R8, 12.058272_SHR_KIND_R8,   &
   &            75.278220_SHR_KIND_R8, 65.241008_SHR_KIND_R8, 64.604291_SHR_KIND_R8,   &
   &             1.647247_SHR_KIND_R8,  7.811584_SHR_KIND_R8, 12.207832_SHR_KIND_R8,   &
   &            63.856665_SHR_KIND_R8, 56.155990_SHR_KIND_R8, 77.448840_SHR_KIND_R8,   &
   &             6.801054_SHR_KIND_R8, 62.209418_SHR_KIND_R8, 20.656133_SHR_KIND_R8,   &
   &            48.344406_SHR_KIND_R8, 55.145460_SHR_KIND_R8, 69.000539_SHR_KIND_R8,   &
   &            11.071350_SHR_KIND_R8, 74.291298_SHR_KIND_R8, 11.047742_SHR_KIND_R8,   &
   &             0.636717_SHR_KIND_R8, 12.844549_SHR_KIND_R8/)

   real   (SHR_KIND_R8), parameter :: obphas(poblen) = & ! phases for obliquity cosine series
   &      (/    251.9025_SHR_KIND_R8, 280.8325_SHR_KIND_R8, 128.3057_SHR_KIND_R8,   &
   &            292.7252_SHR_KIND_R8,  15.3747_SHR_KIND_R8, 263.7951_SHR_KIND_R8,   &
   &            308.4258_SHR_KIND_R8, 240.0099_SHR_KIND_R8, 222.9725_SHR_KIND_R8,   &
   &            268.7809_SHR_KIND_R8, 316.7998_SHR_KIND_R8, 319.6024_SHR_KIND_R8,   &
   &            143.8050_SHR_KIND_R8, 172.7351_SHR_KIND_R8,  28.9300_SHR_KIND_R8,   &
   &            123.5968_SHR_KIND_R8,  20.2082_SHR_KIND_R8,  40.8226_SHR_KIND_R8,   &
   &            123.4722_SHR_KIND_R8, 155.6977_SHR_KIND_R8, 184.6277_SHR_KIND_R8,   &
   &            267.2772_SHR_KIND_R8,  55.0196_SHR_KIND_R8, 152.5268_SHR_KIND_R8,   &
   &             49.1382_SHR_KIND_R8, 204.6609_SHR_KIND_R8,  56.5233_SHR_KIND_R8,   &
   &            200.3284_SHR_KIND_R8, 201.6651_SHR_KIND_R8, 213.5577_SHR_KIND_R8,   &
   &             17.0374_SHR_KIND_R8, 164.4194_SHR_KIND_R8,  94.5422_SHR_KIND_R8,   &
   &            131.9124_SHR_KIND_R8,  61.0309_SHR_KIND_R8, 296.2073_SHR_KIND_R8,   &
   &            135.4894_SHR_KIND_R8, 114.8750_SHR_KIND_R8, 247.0691_SHR_KIND_R8,   &
   &            256.6114_SHR_KIND_R8,  32.1008_SHR_KIND_R8, 143.6804_SHR_KIND_R8,   &
   &             16.8784_SHR_KIND_R8, 160.6835_SHR_KIND_R8,  27.5932_SHR_KIND_R8,   &
   &            348.1074_SHR_KIND_R8,  82.6496_SHR_KIND_R8/)

   ! Cosine/sine series data for computation of eccentricity and fixed vernal
   ! equinox longitude of perihelion (fvelp): amplitude,
   ! rate (arc seconds/year), phase (degrees).

   real   (SHR_KIND_R8), parameter :: ecamp (pecclen) = & ! ampl for eccen/fvelp cos/sin series
   &      (/   0.01860798_SHR_KIND_R8,  0.01627522_SHR_KIND_R8, -0.01300660_SHR_KIND_R8,   &
   &           0.00988829_SHR_KIND_R8, -0.00336700_SHR_KIND_R8,  0.00333077_SHR_KIND_R8,   &
   &          -0.00235400_SHR_KIND_R8,  0.00140015_SHR_KIND_R8,  0.00100700_SHR_KIND_R8,   &
   &           0.00085700_SHR_KIND_R8,  0.00064990_SHR_KIND_R8,  0.00059900_SHR_KIND_R8,   &
   &           0.00037800_SHR_KIND_R8, -0.00033700_SHR_KIND_R8,  0.00027600_SHR_KIND_R8,   &
   &           0.00018200_SHR_KIND_R8, -0.00017400_SHR_KIND_R8, -0.00012400_SHR_KIND_R8,   &
   &           0.00001250_SHR_KIND_R8/)

   real   (SHR_KIND_R8), parameter :: ecrate(pecclen) = & ! rates for eccen/fvelp cos/sin series
   &      (/    4.2072050_SHR_KIND_R8,  7.3460910_SHR_KIND_R8, 17.8572630_SHR_KIND_R8,  &
   &           17.2205460_SHR_KIND_R8, 16.8467330_SHR_KIND_R8,  5.1990790_SHR_KIND_R8,  &
   &           18.2310760_SHR_KIND_R8, 26.2167580_SHR_KIND_R8,  6.3591690_SHR_KIND_R8,  &
   &           16.2100160_SHR_KIND_R8,  3.0651810_SHR_KIND_R8, 16.5838290_SHR_KIND_R8,  &
   &           18.4939800_SHR_KIND_R8,  6.1909530_SHR_KIND_R8, 18.8677930_SHR_KIND_R8,  &
   &           17.4255670_SHR_KIND_R8,  6.1860010_SHR_KIND_R8, 18.4174410_SHR_KIND_R8,  &
   &            0.6678630_SHR_KIND_R8/)

   real   (SHR_KIND_R8), parameter :: ecphas(pecclen) = & ! phases for eccen/fvelp cos/sin series
   &      (/    28.620089_SHR_KIND_R8, 193.788772_SHR_KIND_R8, 308.307024_SHR_KIND_R8,  &
   &           320.199637_SHR_KIND_R8, 279.376984_SHR_KIND_R8,  87.195000_SHR_KIND_R8,  &
   &           349.129677_SHR_KIND_R8, 128.443387_SHR_KIND_R8, 154.143880_SHR_KIND_R8,  &
   &           291.269597_SHR_KIND_R8, 114.860583_SHR_KIND_R8, 332.092251_SHR_KIND_R8,  &
   &           296.414411_SHR_KIND_R8, 145.769910_SHR_KIND_R8, 337.237063_SHR_KIND_R8,  &
   &           152.092288_SHR_KIND_R8, 126.839891_SHR_KIND_R8, 210.667199_SHR_KIND_R8,  &
   &            72.108838_SHR_KIND_R8/)

   ! Sine series data for computation of moving vernal equinox longitude of
   ! perihelion: amplitude (arc seconds), rate (arc sec/year), phase (degrees).

   real   (SHR_KIND_R8), parameter :: mvamp (pmvelen) = & ! amplitudes for mvelp sine series
   &      (/   7391.0225890_SHR_KIND_R8, 2555.1526947_SHR_KIND_R8, 2022.7629188_SHR_KIND_R8,  &
   &          -1973.6517951_SHR_KIND_R8, 1240.2321818_SHR_KIND_R8,  953.8679112_SHR_KIND_R8,  &
   &           -931.7537108_SHR_KIND_R8,  872.3795383_SHR_KIND_R8,  606.3544732_SHR_KIND_R8,  &
   &           -496.0274038_SHR_KIND_R8,  456.9608039_SHR_KIND_R8,  346.9462320_SHR_KIND_R8,  &
   &           -305.8412902_SHR_KIND_R8,  249.6173246_SHR_KIND_R8, -199.1027200_SHR_KIND_R8,  &
   &            191.0560889_SHR_KIND_R8, -175.2936572_SHR_KIND_R8,  165.9068833_SHR_KIND_R8,  &
   &            161.1285917_SHR_KIND_R8,  139.7878093_SHR_KIND_R8, -133.5228399_SHR_KIND_R8,  &
   &            117.0673811_SHR_KIND_R8,  104.6907281_SHR_KIND_R8,   95.3227476_SHR_KIND_R8,  &
   &             86.7824524_SHR_KIND_R8,   86.0857729_SHR_KIND_R8,   70.5893698_SHR_KIND_R8,  &
   &            -69.9719343_SHR_KIND_R8,  -62.5817473_SHR_KIND_R8,   61.5450059_SHR_KIND_R8,  &
   &            -57.9364011_SHR_KIND_R8,   57.1899832_SHR_KIND_R8,  -57.0236109_SHR_KIND_R8,  &
   &            -54.2119253_SHR_KIND_R8,   53.2834147_SHR_KIND_R8,   52.1223575_SHR_KIND_R8,  &
   &            -49.0059908_SHR_KIND_R8,  -48.3118757_SHR_KIND_R8,  -45.4191685_SHR_KIND_R8,  &
   &            -42.2357920_SHR_KIND_R8,  -34.7971099_SHR_KIND_R8,   34.4623613_SHR_KIND_R8,  &
   &            -33.8356643_SHR_KIND_R8,   33.6689362_SHR_KIND_R8,  -31.2521586_SHR_KIND_R8,  &
   &            -30.8798701_SHR_KIND_R8,   28.4640769_SHR_KIND_R8,  -27.1960802_SHR_KIND_R8,  &
   &             27.0860736_SHR_KIND_R8,  -26.3437456_SHR_KIND_R8,   24.7253740_SHR_KIND_R8,  &
   &             24.6732126_SHR_KIND_R8,   24.4272733_SHR_KIND_R8,   24.0127327_SHR_KIND_R8,  &
   &             21.7150294_SHR_KIND_R8,  -21.5375347_SHR_KIND_R8,   18.1148363_SHR_KIND_R8,  &
   &            -16.9603104_SHR_KIND_R8,  -16.1765215_SHR_KIND_R8,   15.5567653_SHR_KIND_R8,  &
   &             15.4846529_SHR_KIND_R8,   15.2150632_SHR_KIND_R8,   14.5047426_SHR_KIND_R8,  &
   &            -14.3873316_SHR_KIND_R8,   13.1351419_SHR_KIND_R8,   12.8776311_SHR_KIND_R8,  &
   &             11.9867234_SHR_KIND_R8,   11.9385578_SHR_KIND_R8,   11.7030822_SHR_KIND_R8,  &
   &             11.6018181_SHR_KIND_R8,  -11.2617293_SHR_KIND_R8,  -10.4664199_SHR_KIND_R8,  &
   &             10.4333970_SHR_KIND_R8,  -10.2377466_SHR_KIND_R8,   10.1934446_SHR_KIND_R8,  &
   &            -10.1280191_SHR_KIND_R8,   10.0289441_SHR_KIND_R8,  -10.0034259_SHR_KIND_R8/)

   real   (SHR_KIND_R8), parameter :: mvrate(pmvelen) = & ! rates for mvelp sine series
   &      (/    31.609974_SHR_KIND_R8, 32.620504_SHR_KIND_R8, 24.172203_SHR_KIND_R8,   &
   &             0.636717_SHR_KIND_R8, 31.983787_SHR_KIND_R8,  3.138886_SHR_KIND_R8,   &
   &            30.973257_SHR_KIND_R8, 44.828336_SHR_KIND_R8,  0.991874_SHR_KIND_R8,   &
   &             0.373813_SHR_KIND_R8, 43.668246_SHR_KIND_R8, 32.246691_SHR_KIND_R8,   &
   &            30.599444_SHR_KIND_R8,  2.147012_SHR_KIND_R8, 10.511172_SHR_KIND_R8,   &
   &            42.681324_SHR_KIND_R8, 13.650058_SHR_KIND_R8,  0.986922_SHR_KIND_R8,   &
   &             9.874455_SHR_KIND_R8, 13.013341_SHR_KIND_R8,  0.262904_SHR_KIND_R8,   &
   &             0.004952_SHR_KIND_R8,  1.142024_SHR_KIND_R8, 63.219948_SHR_KIND_R8,   &
   &             0.205021_SHR_KIND_R8,  2.151964_SHR_KIND_R8, 64.230478_SHR_KIND_R8,   &
   &            43.836462_SHR_KIND_R8, 47.439436_SHR_KIND_R8,  1.384343_SHR_KIND_R8,   &
   &             7.437771_SHR_KIND_R8, 18.829299_SHR_KIND_R8,  9.500642_SHR_KIND_R8,   &
   &             0.431696_SHR_KIND_R8,  1.160090_SHR_KIND_R8, 55.782177_SHR_KIND_R8,   &
   &            12.639528_SHR_KIND_R8,  1.155138_SHR_KIND_R8,  0.168216_SHR_KIND_R8,   &
   &             1.647247_SHR_KIND_R8, 10.884985_SHR_KIND_R8,  5.610937_SHR_KIND_R8,   &
   &            12.658184_SHR_KIND_R8,  1.010530_SHR_KIND_R8,  1.983748_SHR_KIND_R8,   &
   &            14.023871_SHR_KIND_R8,  0.560178_SHR_KIND_R8,  1.273434_SHR_KIND_R8,   &
   &            12.021467_SHR_KIND_R8, 62.583231_SHR_KIND_R8, 63.593761_SHR_KIND_R8,   &
   &            76.438310_SHR_KIND_R8,  4.280910_SHR_KIND_R8, 13.218362_SHR_KIND_R8,   &
   &            17.818769_SHR_KIND_R8,  8.359495_SHR_KIND_R8, 56.792707_SHR_KIND_R8,   &
   &            8.448301_SHR_KIND_R8,  1.978796_SHR_KIND_R8,  8.863925_SHR_KIND_R8,   &
   &             0.186365_SHR_KIND_R8,  8.996212_SHR_KIND_R8,  6.771027_SHR_KIND_R8,   &
   &            45.815258_SHR_KIND_R8, 12.002811_SHR_KIND_R8, 75.278220_SHR_KIND_R8,   &
   &            65.241008_SHR_KIND_R8, 18.870667_SHR_KIND_R8, 22.009553_SHR_KIND_R8,   &
   &            64.604291_SHR_KIND_R8, 11.498094_SHR_KIND_R8,  0.578834_SHR_KIND_R8,   &
   &             9.237738_SHR_KIND_R8, 49.747842_SHR_KIND_R8,  2.147012_SHR_KIND_R8,   &
   &             1.196895_SHR_KIND_R8,  2.133898_SHR_KIND_R8,  0.173168_SHR_KIND_R8/)

   real   (SHR_KIND_R8), parameter :: mvphas(pmvelen) = & ! phases for mvelp sine series
   &      (/    251.9025_SHR_KIND_R8, 280.8325_SHR_KIND_R8, 128.3057_SHR_KIND_R8,   &
   &            348.1074_SHR_KIND_R8, 292.7252_SHR_KIND_R8, 165.1686_SHR_KIND_R8,   &
   &            263.7951_SHR_KIND_R8,  15.3747_SHR_KIND_R8,  58.5749_SHR_KIND_R8,   &
   &             40.8226_SHR_KIND_R8, 308.4258_SHR_KIND_R8, 240.0099_SHR_KIND_R8,   &
   &            222.9725_SHR_KIND_R8, 106.5937_SHR_KIND_R8, 114.5182_SHR_KIND_R8,   &
   &            268.7809_SHR_KIND_R8, 279.6869_SHR_KIND_R8,  39.6448_SHR_KIND_R8,   &
   &            126.4108_SHR_KIND_R8, 291.5795_SHR_KIND_R8, 307.2848_SHR_KIND_R8,   &
   &             18.9300_SHR_KIND_R8, 273.7596_SHR_KIND_R8, 143.8050_SHR_KIND_R8,   &
   &            191.8927_SHR_KIND_R8, 125.5237_SHR_KIND_R8, 172.7351_SHR_KIND_R8,   &
   &            316.7998_SHR_KIND_R8, 319.6024_SHR_KIND_R8,  69.7526_SHR_KIND_R8,   &
   &            123.5968_SHR_KIND_R8, 217.6432_SHR_KIND_R8,  85.5882_SHR_KIND_R8,   &
   &            156.2147_SHR_KIND_R8,  66.9489_SHR_KIND_R8,  20.2082_SHR_KIND_R8,   &
   &            250.7568_SHR_KIND_R8,  48.0188_SHR_KIND_R8,   8.3739_SHR_KIND_R8,   &
   &             17.0374_SHR_KIND_R8, 155.3409_SHR_KIND_R8,  94.1709_SHR_KIND_R8,   &
   &            221.1120_SHR_KIND_R8,  28.9300_SHR_KIND_R8, 117.1498_SHR_KIND_R8,   &
   &            320.5095_SHR_KIND_R8, 262.3602_SHR_KIND_R8, 336.2148_SHR_KIND_R8,   &
   &            233.0046_SHR_KIND_R8, 155.6977_SHR_KIND_R8, 184.6277_SHR_KIND_R8,   &
   &            267.2772_SHR_KIND_R8,  78.9281_SHR_KIND_R8, 123.4722_SHR_KIND_R8,   &
   &            188.7132_SHR_KIND_R8, 180.1364_SHR_KIND_R8,  49.1382_SHR_KIND_R8,   &
   &            152.5268_SHR_KIND_R8,  98.2198_SHR_KIND_R8,  97.4808_SHR_KIND_R8,   &
   &            221.5376_SHR_KIND_R8, 168.2438_SHR_KIND_R8, 161.1199_SHR_KIND_R8,   &
   &             55.0196_SHR_KIND_R8, 262.6495_SHR_KIND_R8, 200.3284_SHR_KIND_R8,   &
   &            201.6651_SHR_KIND_R8, 294.6547_SHR_KIND_R8,  99.8233_SHR_KIND_R8,   &
   &            213.5577_SHR_KIND_R8, 154.1631_SHR_KIND_R8, 232.7153_SHR_KIND_R8,   &
   &            138.3034_SHR_KIND_R8, 204.6609_SHR_KIND_R8, 106.5938_SHR_KIND_R8,   &
   &            250.4676_SHR_KIND_R8, 332.3345_SHR_KIND_R8,  27.3039_SHR_KIND_R8/)

   !---------------------------Local variables----------------------------------
   integer(SHR_KIND_IN) :: i       ! Index for series summations
   real   (SHR_KIND_R8) :: obsum   ! Obliquity series summation
   real   (SHR_KIND_R8) :: cossum  ! Cos series summation for eccentricity/fvelp
   real   (SHR_KIND_R8) :: sinsum  ! Sin series summation for eccentricity/fvelp
   real   (SHR_KIND_R8) :: fvelp   ! Fixed vernal equinox long of perihelion
   real   (SHR_KIND_R8) :: mvsum   ! mvelp series summation
   real   (SHR_KIND_R8) :: beta    ! Intermediate argument for lambm0
   real   (SHR_KIND_R8) :: years   ! Years to time of interest ( pos <=> future)
   real   (SHR_KIND_R8) :: eccen2  ! eccentricity squared
   real   (SHR_KIND_R8) :: eccen3  ! eccentricity cubed

   !-------------------------- Formats -----------------------------------------
   character(len=*),parameter :: F00 = "('(shr_orb_params) ',4a)"
   character(len=*),parameter :: F01 = "('(shr_orb_params) ',a,i9)"
   character(len=*),parameter :: F02 = "('(shr_orb_params) ',a,f6.3)"
   character(len=*),parameter :: F03 = "('(shr_orb_params) ',a,es14.6)"

   !----------------------------------------------------------------------------
   ! radinp and algorithms below will need a degree to radian conversion factor

   if ( log_print .and. s_loglev > 0 ) then
     write(s_logunit,F00) 'Calculate characteristics of the orbit:'
   end if

   ! Check for flag to use input orbit parameters

   IF ( iyear_AD == SHR_ORB_UNDEF_INT ) THEN

      ! Check input obliq, eccen, and mvelp to ensure reasonable

      if( obliq == SHR_ORB_UNDEF_REAL )then
         write(s_logunit,F00) trim(subname)//' Have to specify orbital parameters:'
         write(s_logunit,F00) 'Either set: iyear_AD, OR [obliq, eccen, and mvelp]:'
         write(s_logunit,F00) 'iyear_AD is the year to simulate orbit for (ie. 1950): '
         write(s_logunit,F00) 'obliq, eccen, mvelp specify the orbit directly:'
         write(s_logunit,F00) 'The AMIP II settings (for a 1995 orbit) are: '
         write(s_logunit,F00) ' obliq =  23.4441'
         write(s_logunit,F00) ' eccen =   0.016715'
         write(s_logunit,F00) ' mvelp = 102.7'
         call shr_sys_abort(subname//' ERROR: unreasonable obliq')
      else if ( log_print ) then
         write(s_logunit,F00) 'Use input orbital parameters: '
      end if
      if( (obliq < SHR_ORB_OBLIQ_MIN).or.(obliq > SHR_ORB_OBLIQ_MAX) ) then
         write(s_logunit,F03) 'Input obliquity unreasonable: ', obliq
         call shr_sys_abort(subname//' ERROR: unreasonable obliq')
      end if
      if( (eccen < SHR_ORB_ECCEN_MIN).or.(eccen > SHR_ORB_ECCEN_MAX) ) then
         write(s_logunit,F03) 'Input eccentricity unreasonable: ', eccen
         call shr_sys_abort(subname//' ERROR: unreasonable eccen')
      end if
      if( (mvelp < SHR_ORB_MVELP_MIN).or.(mvelp > SHR_ORB_MVELP_MAX) ) then
         write(s_logunit,F03) 'Input mvelp unreasonable: ' , mvelp
         call shr_sys_abort(subname//' ERROR: unreasonable mvelp')
      end if
      eccen2 = eccen*eccen
      eccen3 = eccen2*eccen

   ELSE  ! Otherwise calculate based on years before present

      if ( log_print .and. s_loglev > 0) then
         write(s_logunit,F01) 'Calculate orbit for year: ' , iyear_AD
      end if
      yb4_1950AD = 1950.0_SHR_KIND_R8 - real(iyear_AD,SHR_KIND_R8)
      if ( abs(yb4_1950AD) .gt. 1000000.0_SHR_KIND_R8 )then
         write(s_logunit,F00) 'orbit only valid for years+-1000000'
         write(s_logunit,F00) 'Relative to 1950 AD'
         write(s_logunit,F03) '# of years before 1950: ',yb4_1950AD
         write(s_logunit,F01) 'Year to simulate was  : ',iyear_AD
         call shr_sys_abort(subname//' ERROR: unreasonable year')
      end if

      ! The following calculates the earths obliquity, orbital eccentricity
      ! (and various powers of it) and vernal equinox mean longitude of
      ! perihelion for years in the past (future = negative of years past),
      ! using constants (see parameter section) given in the program of:
      !
      ! Berger, Andre.  1978  A Simple Algorithm to Compute Long-Term Variations
      ! of Daily Insolation.  Contribution 18, Institute of Astronomy and
      ! Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium.
      !
      ! and formulas given in the paper (where less precise constants are also
      ! given):
      !
      ! Berger, Andre.  1978.  Long-Term Variations of Daily Insolation and
      ! Quaternary Climatic Changes.  J. of the Atmo. Sci. 35:2362-2367
      !
      ! The algorithm is valid only to 1,000,000 years past or hence.
      ! For a solution valid to 5-10 million years past see the above author.
      ! Algorithm below is better for years closer to present than is the
      ! 5-10 million year solution.
      !
      ! Years to time of interest must be negative of years before present
      ! (1950) in formulas that follow.

      years = - yb4_1950AD

      ! In the summations below, cosine or sine arguments, which end up in
      ! degrees, must be converted to radians via multiplication by degrad.
      !
      ! Summation of cosine series for obliquity (epsilon in Berger 1978) in
      ! degrees. Convert the amplitudes and rates, which are in arc secs, into
      ! degrees via multiplication by psecdeg (arc seconds to degrees conversion
      ! factor).  For obliq, first term is Berger 1978 epsilon star; second
      ! term is series summation in degrees.

      obsum = 0.0_SHR_KIND_R8
      do i = 1, poblen
         obsum = obsum + obamp(i)*psecdeg*cos((obrate(i)*psecdeg*years + &
         &       obphas(i))*degrad)
      end do
      obliq = 23.320556_SHR_KIND_R8 + obsum

      ! Summation of cosine and sine series for computation of eccentricity
      ! (eccen; e in Berger 1978) and fixed vernal equinox longitude of
      ! perihelion (fvelp; pi in Berger 1978), which is used for computation
      ! of moving vernal equinox longitude of perihelion.  Convert the rates,
      ! which are in arc seconds, into degrees via multiplication by psecdeg.

      cossum = 0.0_SHR_KIND_R8
      do i = 1, pecclen
        cossum = cossum+ecamp(i)*cos((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do

      sinsum = 0.0_SHR_KIND_R8
      do i = 1, pecclen
        sinsum = sinsum+ecamp(i)*sin((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do

      ! Use summations to calculate eccentricity

      eccen2 = cossum*cossum + sinsum*sinsum
      eccen  = sqrt(eccen2)
      eccen3 = eccen2*eccen

      ! A series of cases for fvelp, which is in radians.

      if (abs(cossum) .le. 1.0E-8_SHR_KIND_R8) then
        if (sinsum .eq. 0.0_SHR_KIND_R8) then
          fvelp = 0.0_SHR_KIND_R8
        else if (sinsum .lt. 0.0_SHR_KIND_R8) then
          fvelp = 1.5_SHR_KIND_R8*pi
        else if (sinsum .gt. 0.0_SHR_KIND_R8) then
          fvelp = .5_SHR_KIND_R8*pi
        endif
      else if (cossum .lt. 0.0_SHR_KIND_R8) then
        fvelp = atan(sinsum/cossum) + pi
      else if (cossum .gt. 0.0_SHR_KIND_R8) then
        if (sinsum .lt. 0.0_SHR_KIND_R8) then
          fvelp = atan(sinsum/cossum) + 2.0_SHR_KIND_R8*pi
        else
          fvelp = atan(sinsum/cossum)
        endif
      endif

      ! Summation of sin series for computation of moving vernal equinox long
      ! of perihelion (mvelp; omega bar in Berger 1978) in degrees.  For mvelp,
      ! first term is fvelp in degrees; second term is Berger 1978 psi bar
      ! times years and in degrees; third term is Berger 1978 zeta; fourth
      ! term is series summation in degrees.  Convert the amplitudes and rates,
      ! which are in arc seconds, into degrees via multiplication by psecdeg.
      ! Series summation plus second and third terms constitute Berger 1978
      ! psi, which is the general precession.

      mvsum = 0.0_SHR_KIND_R8
      do i = 1, pmvelen
        mvsum = mvsum + mvamp(i)*psecdeg*sin((mvrate(i)*psecdeg*years + &
        &       mvphas(i))*degrad)
      end do
      mvelp = fvelp/degrad + 50.439273_SHR_KIND_R8*psecdeg*years + 3.392506_SHR_KIND_R8 + mvsum

      ! Cases to make sure mvelp is between 0 and 360.

      do while (mvelp .lt. 0.0_SHR_KIND_R8)
        mvelp = mvelp + 360.0_SHR_KIND_R8
      end do
      do while (mvelp .ge. 360.0_SHR_KIND_R8)
        mvelp = mvelp - 360.0_SHR_KIND_R8
      end do

   END IF  ! end of test on whether to calculate or use input orbital params

   ! Orbit needs the obliquity in radians

   obliqr = obliq*degrad

   ! 180 degrees must be added to mvelp since observations are made from the
   ! earth and the sun is considered (wrongly for the algorithm) to go around
   ! the earth. For a more graphic explanation see Appendix B in:
   !
   ! A. Berger, M. Loutre and C. Tricot. 1993.  Insolation and Earth Orbital
   ! Periods.  J. of Geophysical Research 98:10,341-10,362.
   !
   ! Additionally, orbit will need this value in radians. So mvelp becomes
   ! mvelpp (mvelp plus pi)

   mvelpp = (mvelp + 180._SHR_KIND_R8)*degrad

   ! Set up an argument used several times in lambm0 calculation ahead.

   beta = sqrt(1._SHR_KIND_R8 - eccen2)

   ! The mean longitude at the vernal equinox (lambda m nought in Berger
   ! 1978; in radians) is calculated from the following formula given in
   ! Berger 1978.  At the vernal equinox the true longitude (lambda in Berger
   ! 1978) is 0.

   lambm0 = 2._SHR_KIND_R8*((.5_SHR_KIND_R8*eccen + .125_SHR_KIND_R8*eccen3)*(1._SHR_KIND_R8 + beta)*sin(mvelpp)  &
   &      - .250_SHR_KIND_R8*eccen2*(.5_SHR_KIND_R8    + beta)*sin(2._SHR_KIND_R8*mvelpp)            &
   &      + .125_SHR_KIND_R8*eccen3*(1._SHR_KIND_R8/3._SHR_KIND_R8 + beta)*sin(3._SHR_KIND_R8*mvelpp))

   if ( log_print ) then
     write(s_logunit,F03) '------ Computed Orbital Parameters ------'
     write(s_logunit,F03) 'Eccentricity      = ',eccen
     write(s_logunit,F03) 'Obliquity (deg)   = ',obliq
     write(s_logunit,F03) 'Obliquity (rad)   = ',obliqr
     write(s_logunit,F03) 'Long of perh(deg) = ',mvelp
     write(s_logunit,F03) 'Long of perh(rad) = ',mvelpp
     write(s_logunit,F03) 'Long at v.e.(rad) = ',lambm0
     write(s_logunit,F03) '-----------------------------------------'
   end if

END SUBROUTINE shr_orb_params

!===============================================================================

SUBROUTINE shr_orb_decl(calday ,eccen ,mvelpp ,lambm0 ,obliqr ,delta ,eccf)

!-------------------------------------------------------------------------------
!
! Compute earth/orbit parameters using formula suggested by
! Duane Thresher.
!
!---------------------------Code history----------------------------------------
!
! Original version:  Erik Kluzek
! Date:              Oct/1997
!
!-------------------------------------------------------------------------------

   !------------------------------Arguments--------------------------------
   real   (SHR_KIND_R8),intent(in)  :: calday ! Calendar day, including fraction
   real   (SHR_KIND_R8),intent(in)  :: eccen  ! Eccentricity
   real   (SHR_KIND_R8),intent(in)  :: obliqr ! Earths obliquity in radians
   real   (SHR_KIND_R8),intent(in)  :: lambm0 ! Mean long of perihelion at the
                                              ! vernal equinox (radians)
   real   (SHR_KIND_R8),intent(in)  :: mvelpp ! moving vernal equinox longitude
                                              ! of perihelion plus pi (radians)
   real   (SHR_KIND_R8),intent(out) :: delta  ! Solar declination angle in rad
   real   (SHR_KIND_R8),intent(out) :: eccf   ! Earth-sun distance factor (ie. (1/r)**2)

   !---------------------------Local variables-----------------------------
   real   (SHR_KIND_R8),parameter :: dayspy = 365.0_SHR_KIND_R8  ! days per year
   real   (SHR_KIND_R8),parameter :: ve     = 80.5_SHR_KIND_R8   ! Calday of vernal equinox
                                                     ! assumes Jan 1 = calday 1

   real   (SHR_KIND_R8) ::   lambm  ! Lambda m, mean long of perihelion (rad)
   real   (SHR_KIND_R8) ::   lmm    ! Intermediate argument involving lambm
   real   (SHR_KIND_R8) ::   lamb   ! Lambda, the earths long of perihelion
   real   (SHR_KIND_R8) ::   invrho ! Inverse normalized sun/earth distance
   real   (SHR_KIND_R8) ::   sinl   ! Sine of lmm

   ! Compute eccentricity factor and solar declination using
   ! day value where a round day (such as 213.0) refers to 0z at
   ! Greenwich longitude.
   !
   ! Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
   ! Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci.
   ! 35:2362-2367.
   !
   ! To get the earths true longitude (position in orbit; lambda in Berger
   ! 1978) which is necessary to find the eccentricity factor and declination,
   ! must first calculate the mean longitude (lambda m in Berger 1978) at
   ! the present day.  This is done by adding to lambm0 (the mean longitude
   ! at the vernal equinox, set as March 21 at noon, when lambda=0; in radians)
   ! an increment (delta lambda m in Berger 1978) that is the number of
   ! days past or before (a negative increment) the vernal equinox divided by
   ! the days in a model year times the 2*pi radians in a complete orbit.

   lambm = lambm0 + (calday - ve)*2._SHR_KIND_R8*pi/dayspy
   lmm   = lambm  - mvelpp

   ! The earths true longitude, in radians, is then found from
   ! the formula in Berger 1978:

   sinl  = sin(lmm)
   lamb  = lambm  + eccen*(2._SHR_KIND_R8*sinl + eccen*(1.25_SHR_KIND_R8*sin(2._SHR_KIND_R8*lmm)  &
   &     + eccen*((13.0_SHR_KIND_R8/12.0_SHR_KIND_R8)*sin(3._SHR_KIND_R8*lmm) - 0.25_SHR_KIND_R8*sinl)))

   ! Using the obliquity, eccentricity, moving vernal equinox longitude of
   ! perihelion (plus), and earths true longitude, the declination (delta)
   ! and the normalized earth/sun distance (rho in Berger 1978; actually inverse
   ! rho will be used), and thus the eccentricity factor (eccf), can be
   ! calculated from formulas given in Berger 1978.

   invrho = (1._SHR_KIND_R8 + eccen*cos(lamb - mvelpp)) / (1._SHR_KIND_R8 - eccen*eccen)

   ! Set solar declination and eccentricity factor

   delta  = asin(sin(obliqr)*sin(lamb))
   eccf   = invrho*invrho

   return

END SUBROUTINE shr_orb_decl

!===============================================================================

SUBROUTINE shr_orb_print( iyear_AD, eccen, obliq, mvelp )

!-------------------------------------------------------------------------------
!
! Print out the information on the Earths input orbital characteristics
!
!---------------------------Code history----------------------------------------
!
! Original version:  Erik Kluzek
! Date:              Oct/1997
!
!-------------------------------------------------------------------------------

   !---------------------------Arguments----------------------------------------
   integer(SHR_KIND_IN),intent(in) :: iyear_AD ! requested Year (AD)
   real   (SHR_KIND_R8),intent(in) :: eccen    ! eccentricity (unitless)
                                               ! (typically 0 to 0.1)
   real   (SHR_KIND_R8),intent(in) :: obliq    ! obliquity (-90 to +90 degrees)
                                               ! typically 22-26
   real   (SHR_KIND_R8),intent(in) :: mvelp    ! moving vernal equinox at perhel
                                               ! (0 to 360 degrees)
   !-------------------------- Formats -----------------------------------------
   character(len=*),parameter :: F00 = "('(shr_orb_print) ',4a)"
   character(len=*),parameter :: F01 = "('(shr_orb_print) ',a,i9.4)"
   character(len=*),parameter :: F02 = "('(shr_orb_print) ',a,f6.3)"
   character(len=*),parameter :: F03 = "('(shr_orb_print) ',a,es14.6)"
   !----------------------------------------------------------------------------

   if (s_loglev > 0) then
   if ( iyear_AD .ne. SHR_ORB_UNDEF_INT ) then
     if ( iyear_AD > 0 ) then
       write(s_logunit,F01) 'Orbital parameters calculated for year: AD ',iyear_AD
     else
       write(s_logunit,F01) 'Orbital parameters calculated for year: BC ',iyear_AD
     end if
   else if ( obliq /= SHR_ORB_UNDEF_REAL ) then
     write(s_logunit,F03) 'Orbital parameters: '
     write(s_logunit,F03) 'Obliquity (degree):              ', obliq
     write(s_logunit,F03) 'Eccentricity (unitless):         ', eccen
     write(s_logunit,F03) 'Long. of moving Perhelion (deg): ', mvelp
   else
     write(s_logunit,F03) 'Orbit parameters not set!'
   end if
   endif

END SUBROUTINE shr_orb_print
!===============================================================================

END MODULE shr_orb_mod
