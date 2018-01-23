!===============================================================================
! CVS: $Id: shr_orb_mod.f90,v 1.1 2004/09/08 19:46:04 samcvs Exp $
! CVS: $Source: /home/disk/eos1/samcvs/cvsroot/SAM/SRC/RAD_CAM/shr_orb_mod.f90,v $
! CVS: $Name:  $
!===============================================================================

MODULE shr_orb_mod

   use shr_kind_mod, only: SHR_KIND_R4, SHR_KIND_IN
!bloss   use shr_sys_mod, only: shr_sys_abort
!bloss   use shr_const_mod, only: SHR_CONST_PI

   IMPLICIT none

   !----------------------------------------------------------------------------
   ! PUBLIC: Interfaces and global data
   !----------------------------------------------------------------------------
   public shr_orb_cosz, shr_orb_params, shr_orb_decl, shr_orb_print
   real   (SHR_KIND_R4),public, parameter :: SHR_ORB_UNDEF_REAL = 1.e36      ! undefined real 
   integer(SHR_KIND_IN),public, parameter :: SHR_ORB_UNDEF_INT  = 2000000000 ! undefined int
   !----------------------------------------------------------------------------
   ! PRIVATE: by default everything else is private to this module
   !----------------------------------------------------------------------------
   private
!bloss   real   (SHR_KIND_R4),parameter :: pi = SHR_CONST_PI
   real   (SHR_KIND_R4), parameter :: pi = 3.14159265358979
   real   (SHR_KIND_R4),parameter :: SHR_ORB_ECCEN_MIN  =   0.0 ! min value for eccen
   real   (SHR_KIND_R4),parameter :: SHR_ORB_ECCEN_MAX  =   0.1 ! max value for eccen
   real   (SHR_KIND_R4),parameter :: SHR_ORB_OBLIQ_MIN  = -90.0 ! min value for obliq
   real   (SHR_KIND_R4),parameter :: SHR_ORB_OBLIQ_MAX  = +90.0 ! max value for obliq
   real   (SHR_KIND_R4),parameter :: SHR_ORB_MVELP_MIN  =   0.0 ! min value for mvelp
   real   (SHR_KIND_R4),parameter :: SHR_ORB_MVELP_MAX  = 360.0 ! max value for mvelp

CONTAINS

!===============================================================================

elemental real(SHR_KIND_R4) FUNCTION shr_orb_cosz(jday,lat,lon,declin)

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

   real   (SHR_KIND_R4),intent(in) :: jday   ! Julian cal day (1.xx to 365.xx)
   real   (SHR_KIND_R4),intent(in) :: lat    ! Centered latitude (radians)
   real   (SHR_KIND_R4),intent(in) :: lon    ! Centered longitude (radians)
   real   (SHR_KIND_R4),intent(in) :: declin ! Solar declination (radians)

   !----------------------------------------------------------------------------

   shr_orb_cosz = sin(lat)*sin(declin) - &
   &              cos(lat)*cos(declin)*cos(jday*2.0*pi + lon)

END FUNCTION shr_orb_cosz

!===============================================================================

SUBROUTINE shr_orb_params( iyear_AD , eccen  , obliq , mvelp     ,     &
           &               obliqr   , lambm0 , mvelpp, log_print )

!-------------------------------------------------------------------------------
!
! Calculate earths orbital parameters using Dave Threshers formula which 
! came from Berger, Andre.  1978  A Simple Algorithm to Compute Long-Term 
! Variations of Daily Insolation.  Contribution 18, Institute of Astronomy 
! and Geophysics, Universite Catholique de Louvain, Louvain-la-Neuve, Belgium
!
!------------------------------Code history-------------------------------------
!
! Original Author: Erik Kluzek
! Date:            Oct/97
!
!-------------------------------------------------------------------------------

   !----------------------------- Arguments ------------------------------------
   real   (SHR_KIND_R4),intent(inout) :: eccen     ! orbital eccentricity
   real   (SHR_KIND_R4),intent(inout) :: obliq     ! obliquity in degrees
   real   (SHR_KIND_R4),intent(inout) :: mvelp     ! moving vernal equinox long
   integer(SHR_KIND_IN),intent(in)    :: iyear_AD  ! Year to calculate orbit for
   logical             ,intent(in)    :: log_print ! Flags print of status/error
   real   (SHR_KIND_R4),intent(out)   :: obliqr    ! Earths obliquity in rad
   real   (SHR_KIND_R4),intent(out)   :: lambm0    ! Mean long of perihelion at
                                                   ! vernal equinox (radians)
   real   (SHR_KIND_R4),intent(out)   :: mvelpp    ! moving vernal equinox long
                                                   ! of perihelion plus pi (rad)

   !------------------------------ Parameters ----------------------------------
   integer(SHR_KIND_IN),parameter :: poblen =47 ! # of elements in series wrt obliquity
   integer(SHR_KIND_IN),parameter :: pecclen=19 ! # of elements in series wrt eccentricity
   integer(SHR_KIND_IN),parameter :: pmvelen=78 ! # of elements in series wrt vernal equinox
   real   (SHR_KIND_R4),parameter :: psecdeg = 1.0/3600.0 ! arc sec to deg conversion

   real   (SHR_KIND_R4) :: degrad = pi/180.   ! degree to radian conversion factor
   real   (SHR_KIND_R4) :: yb4_1950AD         ! number of years before 1950 AD
 
   ! Cosine series data for computation of obliquity: amplitude (arc seconds),
   ! rate (arc seconds/year), phase (degrees).
 
   real   (SHR_KIND_R4), parameter :: obamp(poblen) =  & ! amplitudes for obliquity cos series
   &      (/   -2462.2214466, -857.3232075, -629.3231835,   &
   &            -414.2804924, -311.7632587,  308.9408604,   &
   &            -162.5533601, -116.1077911,  101.1189923,   &
   &             -67.6856209,   24.9079067,   22.5811241,   &
   &             -21.1648355,  -15.6549876,   15.3936813,   &
   &              14.6660938,  -11.7273029,   10.2742696,   &
   &               6.4914588,    5.8539148,   -5.4872205,   &
   &              -5.4290191,    5.1609570,    5.0786314,   &
   &              -4.0735782,    3.7227167,    3.3971932,   &
   &              -2.8347004,   -2.6550721,   -2.5717867,   &
   &              -2.4712188,    2.4625410,    2.2464112,   &
   &              -2.0755511,   -1.9713669,   -1.8813061,   &
   &              -1.8468785,    1.8186742,    1.7601888,   &
   &              -1.5428851,    1.4738838,   -1.4593669,   &
   &               1.4192259,   -1.1818980,    1.1756474,   &
   &              -1.1316126,    1.0896928/)
 
   real   (SHR_KIND_R4), parameter :: obrate(poblen) = & ! rates for obliquity cosine series
   &        (/  31.609974, 32.620504, 24.172203,   &
   &            31.983787, 44.828336, 30.973257,   &
   &            43.668246, 32.246691, 30.599444,   &
   &            42.681324, 43.836462, 47.439436,   &
   &            63.219948, 64.230478,  1.010530,   &
   &             7.437771, 55.782177,  0.373813,   &
   &            13.218362, 62.583231, 63.593761,   &
   &            76.438310, 45.815258,  8.448301,   &
   &            56.792707, 49.747842, 12.058272,   &
   &            75.278220, 65.241008, 64.604291,   &
   &             1.647247,  7.811584, 12.207832,   &
   &            63.856665, 56.155990, 77.448840,   &
   &             6.801054, 62.209418, 20.656133,   &
   &            48.344406, 55.145460, 69.000539,   &
   &            11.071350, 74.291298, 11.047742,   &
   &             0.636717, 12.844549/)
 
   real   (SHR_KIND_R4), parameter :: obphas(poblen) = & ! phases for obliquity cosine series
   &      (/    251.9025, 280.8325, 128.3057,   &
   &            292.7252,  15.3747, 263.7951,   &
   &            308.4258, 240.0099, 222.9725,   &
   &            268.7809, 316.7998, 319.6024,   &
   &            143.8050, 172.7351,  28.9300,   &
   &            123.5968,  20.2082,  40.8226,   &
   &            123.4722, 155.6977, 184.6277,   &
   &            267.2772,  55.0196, 152.5268,   &
   &             49.1382, 204.6609,  56.5233,   &
   &            200.3284, 201.6651, 213.5577,   &
   &             17.0374, 164.4194,  94.5422,   &
   &            131.9124,  61.0309, 296.2073,   &
   &            135.4894, 114.8750, 247.0691,   &
   &            256.6114,  32.1008, 143.6804,   &
   &             16.8784, 160.6835,  27.5932,   &
   &            348.1074,  82.6496/)
 
   ! Cosine/sine series data for computation of eccentricity and fixed vernal 
   ! equinox longitude of perihelion (fvelp): amplitude, 
   ! rate (arc seconds/year), phase (degrees).
 
   real   (SHR_KIND_R4), parameter :: ecamp (pecclen) = & ! ampl for eccen/fvelp cos/sin series
   &      (/   0.01860798,  0.01627522, -0.01300660,   &
   &           0.00988829, -0.00336700,  0.00333077,   &
   &          -0.00235400,  0.00140015,  0.00100700,   &
   &           0.00085700,  0.00064990,  0.00059900,   &
   &           0.00037800, -0.00033700,  0.00027600,   &
   &           0.00018200, -0.00017400, -0.00012400,   &
   &           0.00001250/)
 
   real   (SHR_KIND_R4), parameter :: ecrate(pecclen) = & ! rates for eccen/fvelp cos/sin series
   &      (/    4.2072050,  7.3460910, 17.8572630,  &
   &           17.2205460, 16.8467330,  5.1990790,  &
   &           18.2310760, 26.2167580,  6.3591690,  &
   &           16.2100160,  3.0651810, 16.5838290,  &
   &           18.4939800,  6.1909530, 18.8677930,  &
   &           17.4255670,  6.1860010, 18.4174410,  &
   &            0.6678630/)
 
   real   (SHR_KIND_R4), parameter :: ecphas(pecclen) = & ! phases for eccen/fvelp cos/sin series
   &      (/    28.620089, 193.788772, 308.307024,  &
   &           320.199637, 279.376984,  87.195000,  &
   &           349.129677, 128.443387, 154.143880,  &
   &           291.269597, 114.860583, 332.092251,  &
   &           296.414411, 145.769910, 337.237063,  &
   &           152.092288, 126.839891, 210.667199,  &
   &            72.108838/)
 
   ! Sine series data for computation of moving vernal equinox longitude of 
   ! perihelion: amplitude (arc seconds), rate (arc sec/year), phase (degrees).      
 
   real   (SHR_KIND_R4), parameter :: mvamp (pmvelen) = & ! amplitudes for mvelp sine series 
   &      (/   7391.0225890, 2555.1526947, 2022.7629188,  &
   &          -1973.6517951, 1240.2321818,  953.8679112,  &
   &           -931.7537108,  872.3795383,  606.3544732,  &
   &           -496.0274038,  456.9608039,  346.9462320,  &
   &           -305.8412902,  249.6173246, -199.1027200,  &
   &            191.0560889, -175.2936572,  165.9068833,  &
   &            161.1285917,  139.7878093, -133.5228399,  &
   &            117.0673811,  104.6907281,   95.3227476,  &
   &             86.7824524,   86.0857729,   70.5893698,  &
   &            -69.9719343,  -62.5817473,   61.5450059,  &
   &            -57.9364011,   57.1899832,  -57.0236109,  &
   &            -54.2119253,   53.2834147,   52.1223575,  &
   &            -49.0059908,  -48.3118757,  -45.4191685,  &
   &            -42.2357920,  -34.7971099,   34.4623613,  &
   &            -33.8356643,   33.6689362,  -31.2521586,  &
   &            -30.8798701,   28.4640769,  -27.1960802,  &
   &             27.0860736,  -26.3437456,   24.7253740,  &
   &             24.6732126,   24.4272733,   24.0127327,  &
   &             21.7150294,  -21.5375347,   18.1148363,  &
   &            -16.9603104,  -16.1765215,   15.5567653,  &
   &             15.4846529,   15.2150632,   14.5047426,  &
   &            -14.3873316,   13.1351419,   12.8776311,  &
   &             11.9867234,   11.9385578,   11.7030822,  &
   &             11.6018181,  -11.2617293,  -10.4664199,  &
   &             10.4333970,  -10.2377466,   10.1934446,  &
   &            -10.1280191,   10.0289441,  -10.0034259/)
 
   real   (SHR_KIND_R4), parameter :: mvrate(pmvelen) = & ! rates for mvelp sine series 
   &      (/    31.609974, 32.620504, 24.172203,   &
   &             0.636717, 31.983787,  3.138886,   &
   &            30.973257, 44.828336,  0.991874,   &
   &             0.373813, 43.668246, 32.246691,   &
   &            30.599444,  2.147012, 10.511172,   &
   &            42.681324, 13.650058,  0.986922,   &
   &             9.874455, 13.013341,  0.262904,   &
   &             0.004952,  1.142024, 63.219948,   &
   &             0.205021,  2.151964, 64.230478,   &
   &            43.836462, 47.439436,  1.384343,   &
   &             7.437771, 18.829299,  9.500642,   &
   &             0.431696,  1.160090, 55.782177,   &
   &            12.639528,  1.155138,  0.168216,   &
   &             1.647247, 10.884985,  5.610937,   &
   &            12.658184,  1.010530,  1.983748,   &
   &            14.023871,  0.560178,  1.273434,   &
   &            12.021467, 62.583231, 63.593761,   &
   &            76.438310,  4.280910, 13.218362,   &
   &            17.818769,  8.359495, 56.792707,   &
   &            8.448301,  1.978796,  8.863925,   &
   &             0.186365,  8.996212,  6.771027,   &
   &            45.815258, 12.002811, 75.278220,   &
   &            65.241008, 18.870667, 22.009553,   &
   &            64.604291, 11.498094,  0.578834,   &
   &             9.237738, 49.747842,  2.147012,   &
   &             1.196895,  2.133898,  0.173168/)

   real   (SHR_KIND_R4), parameter :: mvphas(pmvelen) = & ! phases for mvelp sine series
   &      (/    251.9025, 280.8325, 128.3057,   &
   &            348.1074, 292.7252, 165.1686,   &
   &            263.7951,  15.3747,  58.5749,   &
   &             40.8226, 308.4258, 240.0099,   &
   &            222.9725, 106.5937, 114.5182,   &
   &            268.7809, 279.6869,  39.6448,   &
   &            126.4108, 291.5795, 307.2848,   &
   &             18.9300, 273.7596, 143.8050,   &
   &            191.8927, 125.5237, 172.7351,   &
   &            316.7998, 319.6024,  69.7526,   &
   &            123.5968, 217.6432,  85.5882,   &
   &            156.2147,  66.9489,  20.2082,   &
   &            250.7568,  48.0188,   8.3739,   &
   &             17.0374, 155.3409,  94.1709,   &
   &            221.1120,  28.9300, 117.1498,   &
   &            320.5095, 262.3602, 336.2148,   &
   &            233.0046, 155.6977, 184.6277,   &
   &            267.2772,  78.9281, 123.4722,   &
   &            188.7132, 180.1364,  49.1382,   &
   &            152.5268,  98.2198,  97.4808,   &
   &            221.5376, 168.2438, 161.1199,   &
   &             55.0196, 262.6495, 200.3284,   &
   &            201.6651, 294.6547,  99.8233,   &
   &            213.5577, 154.1631, 232.7153,   &
   &            138.3034, 204.6609, 106.5938,   &
   &            250.4676, 332.3345,  27.3039/)
 
   !---------------------------Local variables----------------------------------
   integer(SHR_KIND_IN) :: i       ! Index for series summations
   real   (SHR_KIND_R4) :: obsum   ! Obliquity series summation
   real   (SHR_KIND_R4) :: cossum  ! Cos series summation for eccentricity/fvelp
   real   (SHR_KIND_R4) :: sinsum  ! Sin series summation for eccentricity/fvelp
   real   (SHR_KIND_R4) :: fvelp   ! Fixed vernal equinox long of perihelion
   real   (SHR_KIND_R4) :: mvsum   ! mvelp series summation
   real   (SHR_KIND_R4) :: beta    ! Intermediate argument for lambm0
   real   (SHR_KIND_R4) :: years   ! Years to time of interest ( pos <=> future)
   real   (SHR_KIND_R4) :: eccen2  ! eccentricity squared
   real   (SHR_KIND_R4) :: eccen3  ! eccentricity cubed

   !-------------------------- Formats -----------------------------------------
   character(len=*),parameter :: F00 = "('(shr_orb_params) ',4a)"
   character(len=*),parameter :: F01 = "('(shr_orb_params) ',a,i9)"
   character(len=*),parameter :: F02 = "('(shr_orb_params) ',a,f6.3)"
   character(len=*),parameter :: F03 = "('(shr_orb_params) ',a,es14.6)"

   !----------------------------------------------------------------------------
   ! radinp and algorithms below will need a degree to radian conversion factor
 
   if ( log_print ) then
     write(6,F00) 'Calculate characteristics of the orbit:'
     write(6,F00) 'CVS revision: $Revision: 1.1 $'
     write(6,F00) 'CVS Tag     : $Name:  $'
   end if
 
   ! Check for flag to use input orbit parameters
 
   IF ( iyear_AD == SHR_ORB_UNDEF_INT ) THEN

      ! Check input obliq, eccen, and mvelp to ensure reasonable
 
      if( obliq == SHR_ORB_UNDEF_REAL )then
        if ( log_print ) then
          write(6,F00) 'Have to specify orbital parameters:'
          write(6,F00) 'Either set: iyear_AD, OR [obliq, eccen, and mvelp]:'
          write(6,F00) 'iyear_AD is the year to simulate orbit for (ie. 1950): '
          write(6,F00) 'obliq, eccen, mvelp specify the orbit directly:'
          write(6,F00) 'The AMIP II settings (for a 1995 orbit) are: '
          write(6,F00) ' obliq =  23.4441'
          write(6,F00) ' eccen =   0.016715'
          write(6,F00) ' mvelp = 102.7'
        end if
          call rad_error()
!tcram        call task_abort()
!bloss        call shr_sys_abort()
      else if ( log_print ) then
         write(6,F00) 'Use input orbital parameters: '
      end if
      if( (obliq < SHR_ORB_OBLIQ_MIN).or.(obliq > SHR_ORB_OBLIQ_MAX) ) then
         if ( log_print ) then
            write(6,F03) 'Input obliquity unreasonable: ', obliq
         end if
          call rad_error()
!tcram         call task_abort()
!bloss        call shr_sys_abort()
      end if
      if( (eccen < SHR_ORB_ECCEN_MIN).or.(eccen > SHR_ORB_ECCEN_MAX) ) then
         if ( log_print ) then
            write(6,F03) 'Input eccentricity unreasonable: ', eccen
         end if
          call rad_error()
!tcram         call task_abort()
!bloss        call shr_sys_abort()
      end if
      if( (mvelp < SHR_ORB_MVELP_MIN).or.(mvelp > SHR_ORB_MVELP_MAX) ) then
         if ( log_print ) then
            write(6,F03) 'Input mvelp unreasonable: ' , mvelp
         end if
          call rad_error()
!tcram         call task_abort()
!bloss        call shr_sys_abort()
      end if
      eccen2 = eccen*eccen
      eccen3 = eccen2*eccen

   ELSE  ! Otherwise calculate based on years before present
 
      if ( log_print ) then
         write(6,F01) 'Calculate orbit for year: ' , iyear_AD
      end if
      yb4_1950AD = 1950.0 - float(iyear_AD)
      if ( abs(yb4_1950AD) .gt. 1000000.0 )then
         if ( log_print ) then
            write(6,F00) 'orbit only valid for years+-1000000'
            write(6,F00) 'Relative to 1950 AD'
            write(6,F03) '# of years before 1950: ',yb4_1950AD
            write(6,F01) 'Year to simulate was  : ',iyear_AD
         end if
          call rad_error()
!tcram         call task_abort()
!bloss        call shr_sys_abort()
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
  
      obsum = 0.0
      do i = 1, poblen
         obsum = obsum + obamp(i)*psecdeg*cos((obrate(i)*psecdeg*years + &
         &       obphas(i))*degrad)
      end do
      obliq = 23.320556 + obsum
 
      ! Summation of cosine and sine series for computation of eccentricity 
      ! (eccen; e in Berger 1978) and fixed vernal equinox longitude of 
      ! perihelion (fvelp; pi in Berger 1978), which is used for computation 
      ! of moving vernal equinox longitude of perihelion.  Convert the rates, 
      ! which are in arc seconds, into degrees via multiplication by psecdeg.
 
      cossum = 0.0
      do i = 1, pecclen
        cossum = cossum+ecamp(i)*cos((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do
 
      sinsum = 0.0
      do i = 1, pecclen
        sinsum = sinsum+ecamp(i)*sin((ecrate(i)*psecdeg*years+ecphas(i))*degrad)
      end do
 
      ! Use summations to calculate eccentricity
 
      eccen2 = cossum*cossum + sinsum*sinsum
      eccen  = sqrt(eccen2)
      eccen3 = eccen2*eccen
 
      ! A series of cases for fvelp, which is in radians.
         
      if (abs(cossum) .le. 1.0E-8) then
        if (sinsum .eq. 0.0) then
          fvelp = 0.0
        else if (sinsum .lt. 0.0) then
          fvelp = 1.5*pi
        else if (sinsum .gt. 0.0) then
          fvelp = .5*pi
        endif
      else if (cossum .lt. 0.0) then
        fvelp = atan(sinsum/cossum) + pi
      else if (cossum .gt. 0.0) then
        if (sinsum .lt. 0.0) then
          fvelp = atan(sinsum/cossum) + 2.0*pi
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
 
      mvsum = 0.0
      do i = 1, pmvelen
        mvsum = mvsum + mvamp(i)*psecdeg*sin((mvrate(i)*psecdeg*years + &
        &       mvphas(i))*degrad)
      end do
      mvelp = fvelp/degrad + 50.439273*psecdeg*years + 3.392506 + mvsum
 
      ! Cases to make sure mvelp is between 0 and 360.
 
      do while (mvelp .lt. 0.0)
        mvelp = mvelp + 360.0
      end do
      do while (mvelp .ge. 360.0)
        mvelp = mvelp - 360.0
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
 
   mvelpp = (mvelp + 180.)*degrad
 
   ! Set up an argument used several times in lambm0 calculation ahead.
 
   beta = sqrt(1. - eccen2)
 
   ! The mean longitude at the vernal equinox (lambda m nought in Berger
   ! 1978; in radians) is calculated from the following formula given in 
   ! Berger 1978.  At the vernal equinox the true longitude (lambda in Berger
   ! 1978) is 0.

   lambm0 = 2.*((.5*eccen + .125*eccen3)*(1. + beta)*sin(mvelpp)  &
   &      - .250*eccen2*(.5    + beta)*sin(2.*mvelpp)            &
   &      + .125*eccen3*(1./3. + beta)*sin(3.*mvelpp))
 
   if ( log_print ) then
     write(6,F03) '------ Computed Orbital Parameters ------'
     write(6,F03) 'Eccentricity      = ',eccen
     write(6,F03) 'Obliquity (deg)   = ',obliq
     write(6,F03) 'Obliquity (rad)   = ',obliqr
     write(6,F03) 'Long of perh(deg) = ',mvelp
     write(6,F03) 'Long of perh(rad) = ',mvelpp
     write(6,F03) 'Long at v.e.(rad) = ',lambm0
     write(6,F03) '-----------------------------------------'
   end if
 
END SUBROUTINE shr_orb_params

!===============================================================================

pure SUBROUTINE shr_orb_decl(calday ,eccen ,mvelpp ,lambm0 ,obliqr ,delta ,eccf)

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
   real   (SHR_KIND_R4),intent(in)  :: calday ! Calendar day, including fraction
   real   (SHR_KIND_R4),intent(in)  :: eccen  ! Eccentricity
   real   (SHR_KIND_R4),intent(in)  :: obliqr ! Earths obliquity in radians
   real   (SHR_KIND_R4),intent(in)  :: lambm0 ! Mean long of perihelion at the 
                                              ! vernal equinox (radians)
   real   (SHR_KIND_R4),intent(in)  :: mvelpp ! moving vernal equinox longitude
                                              ! of perihelion plus pi (radians)
   real   (SHR_KIND_R4),intent(out) :: delta  ! Solar declination angle in rad
   real   (SHR_KIND_R4),intent(out) :: eccf   ! Earth-sun distance factor (ie. (1/r)**2)
 
   !---------------------------Local variables-----------------------------
   real   (SHR_KIND_R4),parameter :: dayspy = 365.0  ! days per year
   real   (SHR_KIND_R4),parameter :: ve     = 80.5   ! Calday of vernal equinox
                                                     ! assumes Jan 1 = calday 1
 
   real   (SHR_KIND_R4) ::   lambm  ! Lambda m, mean long of perihelion (rad)
   real   (SHR_KIND_R4) ::   lmm    ! Intermediate argument involving lambm
   real   (SHR_KIND_R4) ::   lamb   ! Lambda, the earths long of perihelion
   real   (SHR_KIND_R4) ::   invrho ! Inverse normalized sun/earth distance
   real   (SHR_KIND_R4) ::   sinl   ! Sine of lmm
 
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
 
   lambm = lambm0 + (calday - ve)*2.*pi/dayspy
   lmm   = lambm  - mvelpp
 
   ! The earths true longitude, in radians, is then found from
   ! the formula in Berger 1978:
 
   sinl  = sin(lmm)
   lamb  = lambm  + eccen*(2.*sinl + eccen*(1.25*sin(2.*lmm)  &
   &     + eccen*((13.0/12.0)*sin(3.*lmm) - 0.25*sinl)))
 
   ! Using the obliquity, eccentricity, moving vernal equinox longitude of
   ! perihelion (plus), and earths true longitude, the declination (delta)
   ! and the normalized earth/sun distance (rho in Berger 1978; actually inverse
   ! rho will be used), and thus the eccentricity factor (eccf), can be 
   ! calculated from formulas given in Berger 1978.
 
   invrho = (1. + eccen*cos(lamb - mvelpp)) / (1. - eccen*eccen)
 
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
   real   (SHR_KIND_R4),intent(in) :: eccen    ! eccentricity (unitless) 
                                               ! (typically 0 to 0.1)
   real   (SHR_KIND_R4),intent(in) :: obliq    ! obliquity (-90 to +90 degrees) 
                                               ! typically 22-26
   real   (SHR_KIND_R4),intent(in) :: mvelp    ! moving vernal equinox at perhel
                                               ! (0 to 360 degrees)
   !-------------------------- Formats -----------------------------------------
   character(len=*),parameter :: F00 = "('(shr_orb_print) ',4a)"
   character(len=*),parameter :: F01 = "('(shr_orb_print) ',a,i9.4)"
   character(len=*),parameter :: F02 = "('(shr_orb_print) ',a,f6.3)"
   character(len=*),parameter :: F03 = "('(shr_orb_print) ',a,es14.6)"
   !----------------------------------------------------------------------------
 
   if ( iyear_AD .ne. SHR_ORB_UNDEF_INT ) then
     if ( iyear_AD > 0 ) then
       write(6,F01) 'Orbital parameters calculated for year: AD ',iyear_AD
     else
       write(6,F01) 'Orbital parameters calculated for year: BC ',iyear_AD
     end if
   else if ( obliq /= SHR_ORB_UNDEF_REAL ) then
     write(6,F03) 'Orbital parameters: '
     write(6,F03) 'Obliquity (degree):              ', obliq
     write(6,F03) 'Eccentricity (unitless):         ', eccen
     write(6,F03) 'Long. of moving Perhelion (deg): ', mvelp
   else
     write(6,F03) 'Orbit parameters not set!'
   end if
 
END SUBROUTINE shr_orb_print
!===============================================================================

END MODULE shr_orb_mod
