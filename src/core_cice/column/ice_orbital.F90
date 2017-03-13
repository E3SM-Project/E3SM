!  SVN:$Id: ice_orbital.F90 1178 2017-03-08 19:24:07Z eclare $
!=======================================================================

! Orbital parameters computed from date
! author:  Bruce P. Briegleb, NCAR 
!
! 2006 ECH: Converted to free source form (F90)
! 2014 ECH: Moved routines from csm_share/shr_orb_mod.F90

      module ice_orbital

      use ice_kinds_mod
      use ice_constants_colpkg, only: c2, p5, pi, secday
      use ice_warnings, only: add_warning

      implicit none
      private
#ifdef CCSMCOUPLED
      public :: compute_coszen
#else
      public :: shr_orb_params, compute_coszen
#endif

!=======================================================================
 
      contains
 
!=======================================================================

! Uses orbital and lat/lon info to compute cosine solar zenith angle
! for the specified date.
!
! author:  Bruce P. Briegleb, NCAR 

      subroutine compute_coszen (tlat,          tlon,     &
                                 calendar_type, days_per_year, &
                                 nextsw_cday,   yday,  sec, &
                                 coszen,        dt)

      use ice_constants_colpkg, only: eccen, mvelpp, lambm0, obliqr, decln, eccf
#ifdef CCSMCOUPLED
      use shr_orb_mod, only: shr_orb_decl
#endif
 
      real (kind=dbl_kind), intent(in) :: &
         tlat, tlon          ! latitude and longitude (radians)

      character (len=char_len), intent(in) :: &
         calendar_type       ! differentiates Gregorian from other calendars

      integer (kind=int_kind), intent(in) :: &
         days_per_year, &    ! number of days in one year
         sec                 ! elapsed seconds into date

      real (kind=dbl_kind), intent(in) :: &
         nextsw_cday     , & ! julian day of next shortwave calculation
         yday                ! day of the year

      real (kind=dbl_kind), intent(inout) :: &
         coszen              ! cosine solar zenith angle 
                             ! negative for sun below horizon
 
      real (kind=dbl_kind), intent(in) :: &
         dt                  ! thermodynamic time step

      ! local variables

      real (kind=dbl_kind) :: ydayp1 ! day of year plus one time step
 
! Solar declination for next time step
 
#ifdef CCSMCOUPLED
      if (calendar_type == "GREGORIAN") then
         ydayp1 = min(nextsw_cday, real(days_per_year,kind=dbl_kind))
      else
         ydayp1 = nextsw_cday
      endif

      !--- update coszen when nextsw_cday valid
      if (ydayp1 > -0.5_dbl_kind) then
#else
      ydayp1 = yday + sec/secday
#endif
 
      call shr_orb_decl(ydayp1, eccen, mvelpp, lambm0, &
                        obliqr, decln, eccf)

      coszen = sin(tlat)*sin(decln) &
             + cos(tlat)*cos(decln) &
             *cos((sec/secday-p5)*c2*pi + tlon) !cos(hour angle)
 
#ifdef CCSMCOUPLED
      endif
#endif

      end subroutine compute_coszen
 
!===============================================================================

#ifndef CCSMCOUPLED
SUBROUTINE shr_orb_params( iyear_AD , eccen , obliq , mvelp    , &
           &               obliqr   , lambm0, mvelpp, log_print, &
                           l_stop, stop_label)

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
   integer(int_kind),intent(in)    :: iyear_AD  ! Year to calculate orbit for
   real   (dbl_kind),intent(inout) :: eccen     ! orbital eccentricity
   real   (dbl_kind),intent(inout) :: obliq     ! obliquity in degrees
   real   (dbl_kind),intent(inout) :: mvelp     ! moving vernal equinox long
   real   (dbl_kind),intent(out)   :: obliqr    ! Earths obliquity in rad
   real   (dbl_kind),intent(out)   :: lambm0    ! Mean long of perihelion at
                                                   ! vernal equinox (radians)
   real   (dbl_kind),intent(out)   :: mvelpp    ! moving vernal equinox long
                                                   ! of perihelion plus pi (rad)
   logical(log_kind),intent(in)    :: log_print ! Flags print of status/error

   logical(log_kind),intent(out)   :: l_stop    ! if true, abort model
   character (len=char_len), intent(out) :: stop_label

   !------------------------------ Parameters ----------------------------------
   real   (dbl_kind),parameter :: SHR_ORB_UNDEF_REAL = 1.e36_dbl_kind ! undefined real 
   integer(int_kind),parameter :: SHR_ORB_UNDEF_INT  = 2000000000        ! undefined int
   integer(int_kind),parameter :: poblen =47 ! # of elements in series wrt obliquity
   integer(int_kind),parameter :: pecclen=19 ! # of elements in series wrt eccentricity
   integer(int_kind),parameter :: pmvelen=78 ! # of elements in series wrt vernal equinox
   real   (dbl_kind),parameter :: psecdeg = 1.0_dbl_kind/3600.0_dbl_kind ! arc sec to deg conversion

   real   (dbl_kind) :: degrad = pi/180._dbl_kind   ! degree to radian conversion factor
   real   (dbl_kind) :: yb4_1950AD         ! number of years before 1950 AD

   real   (dbl_kind),parameter :: SHR_ORB_ECCEN_MIN  =   0.0_dbl_kind ! min value for eccen
   real   (dbl_kind),parameter :: SHR_ORB_ECCEN_MAX  =   0.1_dbl_kind ! max value for eccen
   real   (dbl_kind),parameter :: SHR_ORB_OBLIQ_MIN  = -90.0_dbl_kind ! min value for obliq
   real   (dbl_kind),parameter :: SHR_ORB_OBLIQ_MAX  = +90.0_dbl_kind ! max value for obliq
   real   (dbl_kind),parameter :: SHR_ORB_MVELP_MIN  =   0.0_dbl_kind ! min value for mvelp
   real   (dbl_kind),parameter :: SHR_ORB_MVELP_MAX  = 360.0_dbl_kind ! max value for mvelp

   character(len=*),parameter :: subname = '(shr_orb_params)'
 
   ! Cosine series data for computation of obliquity: amplitude (arc seconds),
   ! rate (arc seconds/year), phase (degrees).
 
   real   (dbl_kind), parameter :: obamp(poblen) =  & ! amplitudes for obliquity cos series
   &      (/   -2462.2214466_dbl_kind, -857.3232075_dbl_kind, -629.3231835_dbl_kind,   &
   &            -414.2804924_dbl_kind, -311.7632587_dbl_kind,  308.9408604_dbl_kind,   &
   &            -162.5533601_dbl_kind, -116.1077911_dbl_kind,  101.1189923_dbl_kind,   &
   &             -67.6856209_dbl_kind,   24.9079067_dbl_kind,   22.5811241_dbl_kind,   &
   &             -21.1648355_dbl_kind,  -15.6549876_dbl_kind,   15.3936813_dbl_kind,   &
   &              14.6660938_dbl_kind,  -11.7273029_dbl_kind,   10.2742696_dbl_kind,   &
   &               6.4914588_dbl_kind,    5.8539148_dbl_kind,   -5.4872205_dbl_kind,   &
   &              -5.4290191_dbl_kind,    5.1609570_dbl_kind,    5.0786314_dbl_kind,   &
   &              -4.0735782_dbl_kind,    3.7227167_dbl_kind,    3.3971932_dbl_kind,   &
   &              -2.8347004_dbl_kind,   -2.6550721_dbl_kind,   -2.5717867_dbl_kind,   &
   &              -2.4712188_dbl_kind,    2.4625410_dbl_kind,    2.2464112_dbl_kind,   &
   &              -2.0755511_dbl_kind,   -1.9713669_dbl_kind,   -1.8813061_dbl_kind,   &
   &              -1.8468785_dbl_kind,    1.8186742_dbl_kind,    1.7601888_dbl_kind,   &
   &              -1.5428851_dbl_kind,    1.4738838_dbl_kind,   -1.4593669_dbl_kind,   &
   &               1.4192259_dbl_kind,   -1.1818980_dbl_kind,    1.1756474_dbl_kind,   &
   &              -1.1316126_dbl_kind,    1.0896928_dbl_kind/)
 
   real   (dbl_kind), parameter :: obrate(poblen) = & ! rates for obliquity cosine series
   &        (/  31.609974_dbl_kind, 32.620504_dbl_kind, 24.172203_dbl_kind,   &
   &            31.983787_dbl_kind, 44.828336_dbl_kind, 30.973257_dbl_kind,   &
   &            43.668246_dbl_kind, 32.246691_dbl_kind, 30.599444_dbl_kind,   &
   &            42.681324_dbl_kind, 43.836462_dbl_kind, 47.439436_dbl_kind,   &
   &            63.219948_dbl_kind, 64.230478_dbl_kind,  1.010530_dbl_kind,   &
   &             7.437771_dbl_kind, 55.782177_dbl_kind,  0.373813_dbl_kind,   &
   &            13.218362_dbl_kind, 62.583231_dbl_kind, 63.593761_dbl_kind,   &
   &            76.438310_dbl_kind, 45.815258_dbl_kind,  8.448301_dbl_kind,   &
   &            56.792707_dbl_kind, 49.747842_dbl_kind, 12.058272_dbl_kind,   &
   &            75.278220_dbl_kind, 65.241008_dbl_kind, 64.604291_dbl_kind,   &
   &             1.647247_dbl_kind,  7.811584_dbl_kind, 12.207832_dbl_kind,   &
   &            63.856665_dbl_kind, 56.155990_dbl_kind, 77.448840_dbl_kind,   &
   &             6.801054_dbl_kind, 62.209418_dbl_kind, 20.656133_dbl_kind,   &
   &            48.344406_dbl_kind, 55.145460_dbl_kind, 69.000539_dbl_kind,   &
   &            11.071350_dbl_kind, 74.291298_dbl_kind, 11.047742_dbl_kind,   &
   &             0.636717_dbl_kind, 12.844549_dbl_kind/)
 
   real   (dbl_kind), parameter :: obphas(poblen) = & ! phases for obliquity cosine series
   &      (/    251.9025_dbl_kind, 280.8325_dbl_kind, 128.3057_dbl_kind,   &
   &            292.7252_dbl_kind,  15.3747_dbl_kind, 263.7951_dbl_kind,   &
   &            308.4258_dbl_kind, 240.0099_dbl_kind, 222.9725_dbl_kind,   &
   &            268.7809_dbl_kind, 316.7998_dbl_kind, 319.6024_dbl_kind,   &
   &            143.8050_dbl_kind, 172.7351_dbl_kind,  28.9300_dbl_kind,   &
   &            123.5968_dbl_kind,  20.2082_dbl_kind,  40.8226_dbl_kind,   &
   &            123.4722_dbl_kind, 155.6977_dbl_kind, 184.6277_dbl_kind,   &
   &            267.2772_dbl_kind,  55.0196_dbl_kind, 152.5268_dbl_kind,   &
   &             49.1382_dbl_kind, 204.6609_dbl_kind,  56.5233_dbl_kind,   &
   &            200.3284_dbl_kind, 201.6651_dbl_kind, 213.5577_dbl_kind,   &
   &             17.0374_dbl_kind, 164.4194_dbl_kind,  94.5422_dbl_kind,   &
   &            131.9124_dbl_kind,  61.0309_dbl_kind, 296.2073_dbl_kind,   &
   &            135.4894_dbl_kind, 114.8750_dbl_kind, 247.0691_dbl_kind,   &
   &            256.6114_dbl_kind,  32.1008_dbl_kind, 143.6804_dbl_kind,   &
   &             16.8784_dbl_kind, 160.6835_dbl_kind,  27.5932_dbl_kind,   &
   &            348.1074_dbl_kind,  82.6496_dbl_kind/)
 
   ! Cosine/sine series data for computation of eccentricity and fixed vernal 
   ! equinox longitude of perihelion (fvelp): amplitude, 
   ! rate (arc seconds/year), phase (degrees).
 
   real   (dbl_kind), parameter :: ecamp (pecclen) = & ! ampl for eccen/fvelp cos/sin series
   &      (/   0.01860798_dbl_kind,  0.01627522_dbl_kind, -0.01300660_dbl_kind,   &
   &           0.00988829_dbl_kind, -0.00336700_dbl_kind,  0.00333077_dbl_kind,   &
   &          -0.00235400_dbl_kind,  0.00140015_dbl_kind,  0.00100700_dbl_kind,   &
   &           0.00085700_dbl_kind,  0.00064990_dbl_kind,  0.00059900_dbl_kind,   &
   &           0.00037800_dbl_kind, -0.00033700_dbl_kind,  0.00027600_dbl_kind,   &
   &           0.00018200_dbl_kind, -0.00017400_dbl_kind, -0.00012400_dbl_kind,   &
   &           0.00001250_dbl_kind/)
 
   real   (dbl_kind), parameter :: ecrate(pecclen) = & ! rates for eccen/fvelp cos/sin series
   &      (/    4.2072050_dbl_kind,  7.3460910_dbl_kind, 17.8572630_dbl_kind,  &
   &           17.2205460_dbl_kind, 16.8467330_dbl_kind,  5.1990790_dbl_kind,  &
   &           18.2310760_dbl_kind, 26.2167580_dbl_kind,  6.3591690_dbl_kind,  &
   &           16.2100160_dbl_kind,  3.0651810_dbl_kind, 16.5838290_dbl_kind,  &
   &           18.4939800_dbl_kind,  6.1909530_dbl_kind, 18.8677930_dbl_kind,  &
   &           17.4255670_dbl_kind,  6.1860010_dbl_kind, 18.4174410_dbl_kind,  &
   &            0.6678630_dbl_kind/)
 
   real   (dbl_kind), parameter :: ecphas(pecclen) = & ! phases for eccen/fvelp cos/sin series
   &      (/    28.620089_dbl_kind, 193.788772_dbl_kind, 308.307024_dbl_kind,  &
   &           320.199637_dbl_kind, 279.376984_dbl_kind,  87.195000_dbl_kind,  &
   &           349.129677_dbl_kind, 128.443387_dbl_kind, 154.143880_dbl_kind,  &
   &           291.269597_dbl_kind, 114.860583_dbl_kind, 332.092251_dbl_kind,  &
   &           296.414411_dbl_kind, 145.769910_dbl_kind, 337.237063_dbl_kind,  &
   &           152.092288_dbl_kind, 126.839891_dbl_kind, 210.667199_dbl_kind,  &
   &            72.108838_dbl_kind/)
 
   ! Sine series data for computation of moving vernal equinox longitude of 
   ! perihelion: amplitude (arc seconds), rate (arc sec/year), phase (degrees).      
 
   real   (dbl_kind), parameter :: mvamp (pmvelen) = & ! amplitudes for mvelp sine series 
   &      (/   7391.0225890_dbl_kind, 2555.1526947_dbl_kind, 2022.7629188_dbl_kind,  &
   &          -1973.6517951_dbl_kind, 1240.2321818_dbl_kind,  953.8679112_dbl_kind,  &
   &           -931.7537108_dbl_kind,  872.3795383_dbl_kind,  606.3544732_dbl_kind,  &
   &           -496.0274038_dbl_kind,  456.9608039_dbl_kind,  346.9462320_dbl_kind,  &
   &           -305.8412902_dbl_kind,  249.6173246_dbl_kind, -199.1027200_dbl_kind,  &
   &            191.0560889_dbl_kind, -175.2936572_dbl_kind,  165.9068833_dbl_kind,  &
   &            161.1285917_dbl_kind,  139.7878093_dbl_kind, -133.5228399_dbl_kind,  &
   &            117.0673811_dbl_kind,  104.6907281_dbl_kind,   95.3227476_dbl_kind,  &
   &             86.7824524_dbl_kind,   86.0857729_dbl_kind,   70.5893698_dbl_kind,  &
   &            -69.9719343_dbl_kind,  -62.5817473_dbl_kind,   61.5450059_dbl_kind,  &
   &            -57.9364011_dbl_kind,   57.1899832_dbl_kind,  -57.0236109_dbl_kind,  &
   &            -54.2119253_dbl_kind,   53.2834147_dbl_kind,   52.1223575_dbl_kind,  &
   &            -49.0059908_dbl_kind,  -48.3118757_dbl_kind,  -45.4191685_dbl_kind,  &
   &            -42.2357920_dbl_kind,  -34.7971099_dbl_kind,   34.4623613_dbl_kind,  &
   &            -33.8356643_dbl_kind,   33.6689362_dbl_kind,  -31.2521586_dbl_kind,  &
   &            -30.8798701_dbl_kind,   28.4640769_dbl_kind,  -27.1960802_dbl_kind,  &
   &             27.0860736_dbl_kind,  -26.3437456_dbl_kind,   24.7253740_dbl_kind,  &
   &             24.6732126_dbl_kind,   24.4272733_dbl_kind,   24.0127327_dbl_kind,  &
   &             21.7150294_dbl_kind,  -21.5375347_dbl_kind,   18.1148363_dbl_kind,  &
   &            -16.9603104_dbl_kind,  -16.1765215_dbl_kind,   15.5567653_dbl_kind,  &
   &             15.4846529_dbl_kind,   15.2150632_dbl_kind,   14.5047426_dbl_kind,  &
   &            -14.3873316_dbl_kind,   13.1351419_dbl_kind,   12.8776311_dbl_kind,  &
   &             11.9867234_dbl_kind,   11.9385578_dbl_kind,   11.7030822_dbl_kind,  &
   &             11.6018181_dbl_kind,  -11.2617293_dbl_kind,  -10.4664199_dbl_kind,  &
   &             10.4333970_dbl_kind,  -10.2377466_dbl_kind,   10.1934446_dbl_kind,  &
   &            -10.1280191_dbl_kind,   10.0289441_dbl_kind,  -10.0034259_dbl_kind/)
 
   real   (dbl_kind), parameter :: mvrate(pmvelen) = & ! rates for mvelp sine series 
   &      (/    31.609974_dbl_kind, 32.620504_dbl_kind, 24.172203_dbl_kind,   &
   &             0.636717_dbl_kind, 31.983787_dbl_kind,  3.138886_dbl_kind,   &
   &            30.973257_dbl_kind, 44.828336_dbl_kind,  0.991874_dbl_kind,   &
   &             0.373813_dbl_kind, 43.668246_dbl_kind, 32.246691_dbl_kind,   &
   &            30.599444_dbl_kind,  2.147012_dbl_kind, 10.511172_dbl_kind,   &
   &            42.681324_dbl_kind, 13.650058_dbl_kind,  0.986922_dbl_kind,   &
   &             9.874455_dbl_kind, 13.013341_dbl_kind,  0.262904_dbl_kind,   &
   &             0.004952_dbl_kind,  1.142024_dbl_kind, 63.219948_dbl_kind,   &
   &             0.205021_dbl_kind,  2.151964_dbl_kind, 64.230478_dbl_kind,   &
   &            43.836462_dbl_kind, 47.439436_dbl_kind,  1.384343_dbl_kind,   &
   &             7.437771_dbl_kind, 18.829299_dbl_kind,  9.500642_dbl_kind,   &
   &             0.431696_dbl_kind,  1.160090_dbl_kind, 55.782177_dbl_kind,   &
   &            12.639528_dbl_kind,  1.155138_dbl_kind,  0.168216_dbl_kind,   &
   &             1.647247_dbl_kind, 10.884985_dbl_kind,  5.610937_dbl_kind,   &
   &            12.658184_dbl_kind,  1.010530_dbl_kind,  1.983748_dbl_kind,   &
   &            14.023871_dbl_kind,  0.560178_dbl_kind,  1.273434_dbl_kind,   &
   &            12.021467_dbl_kind, 62.583231_dbl_kind, 63.593761_dbl_kind,   &
   &            76.438310_dbl_kind,  4.280910_dbl_kind, 13.218362_dbl_kind,   &
   &            17.818769_dbl_kind,  8.359495_dbl_kind, 56.792707_dbl_kind,   &
   &            8.448301_dbl_kind,  1.978796_dbl_kind,  8.863925_dbl_kind,   &
   &             0.186365_dbl_kind,  8.996212_dbl_kind,  6.771027_dbl_kind,   &
   &            45.815258_dbl_kind, 12.002811_dbl_kind, 75.278220_dbl_kind,   &
   &            65.241008_dbl_kind, 18.870667_dbl_kind, 22.009553_dbl_kind,   &
   &            64.604291_dbl_kind, 11.498094_dbl_kind,  0.578834_dbl_kind,   &
   &             9.237738_dbl_kind, 49.747842_dbl_kind,  2.147012_dbl_kind,   &
   &             1.196895_dbl_kind,  2.133898_dbl_kind,  0.173168_dbl_kind/)

   real   (dbl_kind), parameter :: mvphas(pmvelen) = & ! phases for mvelp sine series
   &      (/    251.9025_dbl_kind, 280.8325_dbl_kind, 128.3057_dbl_kind,   &
   &            348.1074_dbl_kind, 292.7252_dbl_kind, 165.1686_dbl_kind,   &
   &            263.7951_dbl_kind,  15.3747_dbl_kind,  58.5749_dbl_kind,   &
   &             40.8226_dbl_kind, 308.4258_dbl_kind, 240.0099_dbl_kind,   &
   &            222.9725_dbl_kind, 106.5937_dbl_kind, 114.5182_dbl_kind,   &
   &            268.7809_dbl_kind, 279.6869_dbl_kind,  39.6448_dbl_kind,   &
   &            126.4108_dbl_kind, 291.5795_dbl_kind, 307.2848_dbl_kind,   &
   &             18.9300_dbl_kind, 273.7596_dbl_kind, 143.8050_dbl_kind,   &
   &            191.8927_dbl_kind, 125.5237_dbl_kind, 172.7351_dbl_kind,   &
   &            316.7998_dbl_kind, 319.6024_dbl_kind,  69.7526_dbl_kind,   &
   &            123.5968_dbl_kind, 217.6432_dbl_kind,  85.5882_dbl_kind,   &
   &            156.2147_dbl_kind,  66.9489_dbl_kind,  20.2082_dbl_kind,   &
   &            250.7568_dbl_kind,  48.0188_dbl_kind,   8.3739_dbl_kind,   &
   &             17.0374_dbl_kind, 155.3409_dbl_kind,  94.1709_dbl_kind,   &
   &            221.1120_dbl_kind,  28.9300_dbl_kind, 117.1498_dbl_kind,   &
   &            320.5095_dbl_kind, 262.3602_dbl_kind, 336.2148_dbl_kind,   &
   &            233.0046_dbl_kind, 155.6977_dbl_kind, 184.6277_dbl_kind,   &
   &            267.2772_dbl_kind,  78.9281_dbl_kind, 123.4722_dbl_kind,   &
   &            188.7132_dbl_kind, 180.1364_dbl_kind,  49.1382_dbl_kind,   &
   &            152.5268_dbl_kind,  98.2198_dbl_kind,  97.4808_dbl_kind,   &
   &            221.5376_dbl_kind, 168.2438_dbl_kind, 161.1199_dbl_kind,   &
   &             55.0196_dbl_kind, 262.6495_dbl_kind, 200.3284_dbl_kind,   &
   &            201.6651_dbl_kind, 294.6547_dbl_kind,  99.8233_dbl_kind,   &
   &            213.5577_dbl_kind, 154.1631_dbl_kind, 232.7153_dbl_kind,   &
   &            138.3034_dbl_kind, 204.6609_dbl_kind, 106.5938_dbl_kind,   &
   &            250.4676_dbl_kind, 332.3345_dbl_kind,  27.3039_dbl_kind/)
 
   !---------------------------Local variables----------------------------------
   integer(int_kind) :: i       ! Index for series summations
   real   (dbl_kind) :: obsum   ! Obliquity series summation
   real   (dbl_kind) :: cossum  ! Cos series summation for eccentricity/fvelp
   real   (dbl_kind) :: sinsum  ! Sin series summation for eccentricity/fvelp
   real   (dbl_kind) :: fvelp   ! Fixed vernal equinox long of perihelion
   real   (dbl_kind) :: mvsum   ! mvelp series summation
   real   (dbl_kind) :: beta    ! Intermediate argument for lambm0
   real   (dbl_kind) :: years   ! Years to time of interest ( pos <=> future)
   real   (dbl_kind) :: eccen2  ! eccentricity squared
   real   (dbl_kind) :: eccen3  ! eccentricity cubed
   integer (int_kind), parameter :: s_loglev    = 0         
   character(len=char_len_long) :: warning ! warning message

   !-------------------------- Formats -----------------------------------------
   character(*),parameter :: svnID  = "SVN " // &
   "$Id: ice_orbital.F90 1178 2017-03-08 19:24:07Z eclare $"
   character(*),parameter :: svnURL = "SVN <unknown URL>" 
!  character(*),parameter :: svnURL = "SVN " // &
!  "$URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_121022/shr/shr_orb_mod.F90 $"
   character(len=*),parameter :: F00 = "('(shr_orb_params) ',4a)"
   character(len=*),parameter :: F01 = "('(shr_orb_params) ',a,i9)"
   character(len=*),parameter :: F02 = "('(shr_orb_params) ',a,f6.3)"
   character(len=*),parameter :: F03 = "('(shr_orb_params) ',a,es14.6)"

   !----------------------------------------------------------------------------
   ! radinp and algorithms below will need a degree to radian conversion factor

   l_stop = .false.
   stop_label = ' '
 
   if ( log_print .and. s_loglev > 0 ) then
     write(warning,F00) 'Calculate characteristics of the orbit:'
     call add_warning(warning)
     write(warning,F00) svnID
     call add_warning(warning)
!    write(warning,F00) svnURL
!    call add_warning(warning)
   end if
 
   ! Check for flag to use input orbit parameters
 
   IF ( iyear_AD == SHR_ORB_UNDEF_INT ) THEN

      ! Check input obliq, eccen, and mvelp to ensure reasonable
 
      if( obliq == SHR_ORB_UNDEF_REAL )then
         write(warning,F00) trim(subname)//' Have to specify orbital parameters:'
         call add_warning(warning)
         write(warning,F00) 'Either set: iyear_AD, OR [obliq, eccen, and mvelp]:'
         call add_warning(warning)
         write(warning,F00) 'iyear_AD is the year to simulate orbit for (ie. 1950): '
         call add_warning(warning)
         write(warning,F00) 'obliq, eccen, mvelp specify the orbit directly:'
         call add_warning(warning)
         write(warning,F00) 'The AMIP II settings (for a 1995 orbit) are: '
         call add_warning(warning)
         write(warning,F00) ' obliq =  23.4441'
         call add_warning(warning)
         write(warning,F00) ' eccen =   0.016715'
         call add_warning(warning)
         write(warning,F00) ' mvelp = 102.7'
         call add_warning(warning)
         l_stop = .true.
         stop_label = 'unreasonable oblip'
      else if ( log_print ) then
         write(warning,F00) 'Use input orbital parameters: '
         call add_warning(warning)
      end if
      if( (obliq < SHR_ORB_OBLIQ_MIN).or.(obliq > SHR_ORB_OBLIQ_MAX) ) then
         write(warning,F03) 'Input obliquity unreasonable: ', obliq
         call add_warning(warning)
         l_stop = .true.
         stop_label = 'unreasonable obliq'
      end if
      if( (eccen < SHR_ORB_ECCEN_MIN).or.(eccen > SHR_ORB_ECCEN_MAX) ) then
         write(warning,F03) 'Input eccentricity unreasonable: ', eccen
         call add_warning(warning)
         l_stop = .true.
         stop_label = 'unreasonable eccen'
      end if
      if( (mvelp < SHR_ORB_MVELP_MIN).or.(mvelp > SHR_ORB_MVELP_MAX) ) then
         write(warning,F03) 'Input mvelp unreasonable: ' , mvelp
         call add_warning(warning)
         l_stop = .true.
         stop_label = 'unreasonable mvelp'
      end if
      eccen2 = eccen*eccen
      eccen3 = eccen2*eccen

   ELSE  ! Otherwise calculate based on years before present
 
      if ( log_print .and. s_loglev > 0) then
         write(warning,F01) 'Calculate orbit for year: ' , iyear_AD
         call add_warning(warning)
      end if
      yb4_1950AD = 1950.0_dbl_kind - real(iyear_AD,dbl_kind)
      if ( abs(yb4_1950AD) .gt. 1000000.0_dbl_kind )then
         write(warning,F00) 'orbit only valid for years+-1000000'
         call add_warning(warning)
         write(warning,F00) 'Relative to 1950 AD'
         call add_warning(warning)
         write(warning,F03) '# of years before 1950: ',yb4_1950AD
         call add_warning(warning)
         write(warning,F01) 'Year to simulate was  : ',iyear_AD
         call add_warning(warning)
         l_stop = .true.
         stop_label = 'unreasonable year'
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
  
      obsum = 0.0_dbl_kind
      do i = 1, poblen
         obsum = obsum + obamp(i)*psecdeg*cos((obrate(i)*psecdeg*years + &
         &       obphas(i))*degrad)
      end do
      obliq = 23.320556_dbl_kind + obsum
 
      ! Summation of cosine and sine series for computation of eccentricity 
      ! (eccen; e in Berger 1978) and fixed vernal equinox longitude of 
      ! perihelion (fvelp; pi in Berger 1978), which is used for computation 
      ! of moving vernal equinox longitude of perihelion.  Convert the rates, 
      ! which are in arc seconds, into degrees via multiplication by psecdeg.
 
      cossum = 0.0_dbl_kind
      do i = 1, pecclen
        cossum = cossum+ecamp(i)*cos((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do
 
      sinsum = 0.0_dbl_kind
      do i = 1, pecclen
        sinsum = sinsum+ecamp(i)*sin((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do
 
      ! Use summations to calculate eccentricity
 
      eccen2 = cossum*cossum + sinsum*sinsum
      eccen  = sqrt(eccen2)
      eccen3 = eccen2*eccen
 
      ! A series of cases for fvelp, which is in radians.
         
      if (abs(cossum) .le. 1.0E-8_dbl_kind) then
        if (sinsum .eq. 0.0_dbl_kind) then
          fvelp = 0.0_dbl_kind
        else if (sinsum .lt. 0.0_dbl_kind) then
          fvelp = 1.5_dbl_kind*pi
        else if (sinsum .gt. 0.0_dbl_kind) then
          fvelp = .5_dbl_kind*pi
        endif
      else if (cossum .lt. 0.0_dbl_kind) then
        fvelp = atan(sinsum/cossum) + pi
      else if (cossum .gt. 0.0_dbl_kind) then
        if (sinsum .lt. 0.0_dbl_kind) then
          fvelp = atan(sinsum/cossum) + 2.0_dbl_kind*pi
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
 
      mvsum = 0.0_dbl_kind
      do i = 1, pmvelen
        mvsum = mvsum + mvamp(i)*psecdeg*sin((mvrate(i)*psecdeg*years + &
        &       mvphas(i))*degrad)
      end do
      mvelp = fvelp/degrad + 50.439273_dbl_kind*psecdeg*years + 3.392506_dbl_kind + mvsum
 
      ! Cases to make sure mvelp is between 0 and 360.
 
      do while (mvelp .lt. 0.0_dbl_kind)
        mvelp = mvelp + 360.0_dbl_kind
      end do
      do while (mvelp .ge. 360.0_dbl_kind)
        mvelp = mvelp - 360.0_dbl_kind
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
 
   mvelpp = (mvelp + 180._dbl_kind)*degrad
 
   ! Set up an argument used several times in lambm0 calculation ahead.
 
   beta = sqrt(1._dbl_kind - eccen2)
 
   ! The mean longitude at the vernal equinox (lambda m nought in Berger
   ! 1978; in radians) is calculated from the following formula given in 
   ! Berger 1978.  At the vernal equinox the true longitude (lambda in Berger
   ! 1978) is 0.

   lambm0 = 2._dbl_kind*((.5_dbl_kind*eccen + .125_dbl_kind*eccen3)*(1._dbl_kind + beta)*sin(mvelpp)  &
   &      - .250_dbl_kind*eccen2*(.5_dbl_kind    + beta)*sin(2._dbl_kind*mvelpp)            &
   &      + .125_dbl_kind*eccen3*(1._dbl_kind/3._dbl_kind + beta)*sin(3._dbl_kind*mvelpp))
 
   if ( log_print ) then
     write(warning,F03) '------ Computed Orbital Parameters ------'
     call add_warning(warning)
     write(warning,F03) 'Eccentricity      = ',eccen
     call add_warning(warning)
     write(warning,F03) 'Obliquity (deg)   = ',obliq
     call add_warning(warning)
     write(warning,F03) 'Obliquity (rad)   = ',obliqr
     call add_warning(warning)
     write(warning,F03) 'Long of perh(deg) = ',mvelp
     call add_warning(warning)
     write(warning,F03) 'Long of perh(rad) = ',mvelpp
     call add_warning(warning)
     write(warning,F03) 'Long at v.e.(rad) = ',lambm0
     call add_warning(warning)
     write(warning,F03) '-----------------------------------------'
     call add_warning(warning)
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
   real   (dbl_kind),intent(in)  :: calday ! Calendar day, including fraction
   real   (dbl_kind),intent(in)  :: eccen  ! Eccentricity
   real   (dbl_kind),intent(in)  :: obliqr ! Earths obliquity in radians
   real   (dbl_kind),intent(in)  :: lambm0 ! Mean long of perihelion at the 
                                              ! vernal equinox (radians)
   real   (dbl_kind),intent(in)  :: mvelpp ! moving vernal equinox longitude
                                              ! of perihelion plus pi (radians)
   real   (dbl_kind),intent(out) :: delta  ! Solar declination angle in rad
   real   (dbl_kind),intent(out) :: eccf   ! Earth-sun distance factor (ie. (1/r)**2)
 
   !---------------------------Local variables-----------------------------
   real   (dbl_kind),parameter :: dayspy = 365.0_dbl_kind  ! days per year
   real   (dbl_kind),parameter :: ve     = 80.5_dbl_kind   ! Calday of vernal equinox
                                                     ! assumes Jan 1 = calday 1
 
   real   (dbl_kind) ::   lambm  ! Lambda m, mean long of perihelion (rad)
   real   (dbl_kind) ::   lmm    ! Intermediate argument involving lambm
   real   (dbl_kind) ::   lamb   ! Lambda, the earths long of perihelion
   real   (dbl_kind) ::   invrho ! Inverse normalized sun/earth distance
   real   (dbl_kind) ::   sinl   ! Sine of lmm
 
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
 
   lambm = lambm0 + (calday - ve)*2._dbl_kind*pi/dayspy
   lmm   = lambm  - mvelpp
 
   ! The earths true longitude, in radians, is then found from
   ! the formula in Berger 1978:
 
   sinl  = sin(lmm)
   lamb  = lambm  + eccen*(2._dbl_kind*sinl + eccen*(1.25_dbl_kind*sin(2._dbl_kind*lmm)  &
   &     + eccen*((13.0_dbl_kind/12.0_dbl_kind)*sin(3._dbl_kind*lmm) - 0.25_dbl_kind*sinl)))
 
   ! Using the obliquity, eccentricity, moving vernal equinox longitude of
   ! perihelion (plus), and earths true longitude, the declination (delta)
   ! and the normalized earth/sun distance (rho in Berger 1978; actually inverse
   ! rho will be used), and thus the eccentricity factor (eccf), can be 
   ! calculated from formulas given in Berger 1978.
 
   invrho = (1._dbl_kind + eccen*cos(lamb - mvelpp)) / (1._dbl_kind - eccen*eccen)
 
   ! Set solar declination and eccentricity factor
 
   delta  = asin(sin(obliqr)*sin(lamb))
   eccf   = invrho*invrho
 
   return
 
END SUBROUTINE shr_orb_decl
#endif

!=======================================================================
 
      end module ice_orbital
 
!=======================================================================
