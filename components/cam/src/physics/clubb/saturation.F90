!-------------------------------------------------------------------------
!$Id: saturation.F90 7991 2016-02-23 20:36:03Z raut@uwm.edu $
!===============================================================================
module saturation

! Description:
!   Contains functions that compute saturation with respect
!   to liquid or ice.
!-----------------------------------------------------------------------

#ifdef GFDL
    use model_flags, only: &  ! h1g, 2010-06-18
       I_sat_sphum
#endif

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  private  ! Change default so all items private

  public   :: sat_mixrat_liq, sat_mixrat_liq_lookup, sat_mixrat_ice, rcm_sat_adj, &
              sat_vapor_press_liq

  private  :: sat_vapor_press_liq_flatau, sat_vapor_press_liq_bolton
  private  :: sat_vapor_press_ice_flatau, sat_vapor_press_ice_bolton

  ! Lookup table of values for saturation 
  real( kind = core_rknd ), private, dimension(188:343) :: &
    svp_liq_lookup_table

  data svp_liq_lookup_table(188:343) / &
    0.049560547_core_rknd, 0.059753418_core_rknd, 0.070129395_core_rknd, 0.083618164_core_rknd, &
    0.09814453_core_rknd, 0.11444092_core_rknd, 0.13446045_core_rknd, 0.15686035_core_rknd,     &
    0.18218994_core_rknd, 0.21240234_core_rknd, 0.24725342_core_rknd, 0.28668213_core_rknd,     &
    0.33184814_core_rknd, 0.3826294_core_rknd, 0.4416504_core_rknd, 0.50775146_core_rknd,       &
    0.58343506_core_rknd, 0.6694946_core_rknd, 0.7668457_core_rknd, 0.87750244_core_rknd,       &
    1.0023804_core_rknd, 1.1434937_core_rknd, 1.3028564_core_rknd, 1.482544_core_rknd,          &
    1.6847534_core_rknd, 1.9118042_core_rknd, 2.1671143_core_rknd, 2.4535522_core_rknd,         &
    2.774231_core_rknd, 3.1330566_core_rknd, 3.5343628_core_rknd, 3.9819336_core_rknd,          &
    4.480713_core_rknd, 5.036072_core_rknd, 5.6540527_core_rknd, 6.340088_core_rknd,            &
    7.1015015_core_rknd, 7.9450684_core_rknd, 8.8793335_core_rknd, 9.91217_core_rknd,           &
    11.053528_core_rknd, 12.313049_core_rknd, 13.70166_core_rknd, 15.231018_core_rknd,          &
    16.91394_core_rknd, 18.764038_core_rknd, 20.795898_core_rknd, 23.025574_core_rknd,          &
    25.470093_core_rknd, 28.147766_core_rknd, 31.078003_core_rknd, 34.282043_core_rknd,         &
    37.782593_core_rknd, 41.60382_core_rknd, 45.771606_core_rknd, 50.31366_core_rknd,           &
    55.259644_core_rknd, 60.641174_core_rknd, 66.492004_core_rknd, 72.84802_core_rknd,          &
    79.74756_core_rknd, 87.23126_core_rknd, 95.34259_core_rknd, 104.12747_core_rknd,            &
    113.634796_core_rknd, 123.91641_core_rknd, 135.02725_core_rknd, 147.02563_core_rknd,        &
    159.97308_core_rknd, 173.93488_core_rknd, 188.97995_core_rknd, 205.18109_core_rknd,         &
    222.61517_core_rknd, 241.36334_core_rknd, 261.51108_core_rknd, 283.14853_core_rknd,         &
    306.37054_core_rknd, 331.27698_core_rknd, 357.97278_core_rknd, 386.56842_core_rknd,         &
    417.17978_core_rknd, 449.9286_core_rknd, 484.94254_core_rknd, 522.3556_core_rknd,           &
    562.30804_core_rknd, 604.947_core_rknd, 650.42645_core_rknd, 698.9074_core_rknd,            &
    750.55835_core_rknd, 805.55554_core_rknd, 864.0828_core_rknd, 926.3325_core_rknd,           &
    992.5052_core_rknd, 1062.8102_core_rknd, 1137.4657_core_rknd, 1216.6995_core_rknd,          &
    1300.7483_core_rknd, 1389.8594_core_rknd, 1484.2896_core_rknd, 1584.3064_core_rknd,         &
    1690.1881_core_rknd, 1802.224_core_rknd, 1920.7146_core_rknd, 2045.9724_core_rknd,          &
    2178.3218_core_rknd, 2318.099_core_rknd, 2465.654_core_rknd, 2621.3489_core_rknd,           &
    2785.5596_core_rknd, 2958.6758_core_rknd, 3141.101_core_rknd, 3333.2534_core_rknd,          &
    3535.5657_core_rknd, 3748.4863_core_rknd, 3972.4792_core_rknd, 4208.024_core_rknd,          &
    4455.616_core_rknd, 4715.7686_core_rknd, 4989.0127_core_rknd, 5275.8945_core_rknd,          &
    5576.9795_core_rknd, 5892.8535_core_rknd, 6224.116_core_rknd, 6571.3926_core_rknd,          &
    6935.3213_core_rknd, 7316.5674_core_rknd, 7715.8105_core_rknd, 8133.755_core_rknd,          &
    8571.125_core_rknd, 9028.667_core_rknd, 9507.15_core_rknd, 10007.367_core_rknd,             &
    10530.132_core_rknd, 11076.282_core_rknd, 11646.683_core_rknd, 12242.221_core_rknd,         &
    12863.808_core_rknd, 13512.384_core_rknd, 14188.913_core_rknd, 14894.385_core_rknd,         &
    15629.823_core_rknd, 16396.268_core_rknd, 17194.799_core_rknd, 18026.516_core_rknd,         &
    18892.55_core_rknd, 19794.07_core_rknd, 20732.262_core_rknd, 21708.352_core_rknd,           &
    22723.592_core_rknd, 23779.273_core_rknd, 24876.709_core_rknd, 26017.258_core_rknd,         &
    27202.3_core_rknd, 28433.256_core_rknd, 29711.578_core_rknd, 31038.766_core_rknd /

!$omp threadprivate( svp_liq_lookup_table )

  contains

!-------------------------------------------------------------------------
  elemental real( kind = core_rknd ) function sat_mixrat_liq( p_in_Pa, T_in_K )

! Description:
!   Used to compute the saturation mixing ratio of liquid water.

! References:
!   Formula from Emanuel 1994, 4.4.14
!-------------------------------------------------------------------------

    use constants_clubb, only: & 
      ep    ! Variable

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      p_in_Pa,  & ! Pressure    [Pa]
      T_in_K      ! Temperature [K]

   ! Local Variables
    real( kind = core_rknd ) :: esatv

    ! --- Begin Code ---

    ! Calculate the SVP for water vapor.
    esatv = sat_vapor_press_liq( T_in_K )

    ! If esatv exceeds the air pressure, then assume esatv~=0.5*pressure 
    !   and set rsat = ep = 0.622
    if ( p_in_Pa-esatv < 1.0_core_rknd ) then
      sat_mixrat_liq = ep
    else

#ifdef GFDL

    ! GFDL uses specific humidity
    ! Formula for Saturation Specific Humidity
     if ( I_sat_sphum )  then   ! h1g, 2010-06-18 begin mod
       sat_mixrat_liq = ep * ( esatv / ( p_in_Pa - (1.0_core_rknd-ep) * esatv ) )
     else
       sat_mixrat_liq = ep * ( esatv / ( p_in_Pa - esatv ) )
     endif                     ! h1g, 2010-06-18 end mod
#else
    ! Formula for Saturation Mixing Ratio:
    !
    ! rs = (epsilon) * [ esat / ( p - esat ) ];
    ! where epsilon = R_d / R_v
    sat_mixrat_liq = ep * ( esatv / ( p_in_Pa - esatv ) )
#endif

    end if

    return
  end function sat_mixrat_liq

!-------------------------------------------------------------------------
  elemental real( kind = core_rknd ) function sat_mixrat_liq_lookup( p_in_Pa, T_in_K )

! Description:
!   Used to compute the saturation mixing ratio of liquid water.
!   This function utilizes sat_vapor_press_liq_lookup; the SVP is found
!   using a lookup table rather than calculating it using various
!   approximations.

! References:
!   Formula from Emanuel 1994, 4.4.14
!-------------------------------------------------------------------------

    use constants_clubb, only: & 
      ep    ! Variable

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      p_in_Pa,  & ! Pressure    [Pa]
      T_in_K      ! Temperature [K]

   ! Local Variables
    real( kind = core_rknd ) :: esatv

    ! --- Begin Code ---

    ! Calculate the SVP for water vapor using a lookup table.
    esatv = sat_vapor_press_liq_lookup( T_in_K )

    ! If esatv exceeds the air pressure, then assume esatv~=0.5*pressure 
    !   and set rsat = ep = 0.622
    if ( p_in_Pa-esatv < 1.0_core_rknd ) then
      sat_mixrat_liq_lookup = ep
    else

#ifdef GFDL

    ! GFDL uses specific humidity
    ! Formula for Saturation Specific Humidity
     if( I_sat_sphum )  then   ! h1g, 2010-06-18 begin mod
           sat_mixrat_liq_lookup = ep * ( esatv / ( p_in_Pa - (1.0_core_rknd-ep) * esatv ) )
     else
           sat_mixrat_liq_lookup = ep * ( esatv / ( p_in_Pa - esatv ) )
     endif                     ! h1g, 2010-06-18 end mod
#else
    ! Formula for Saturation Mixing Ratio:
    !
    ! rs = (epsilon) * [ esat / ( p - esat ) ];
    ! where epsilon = R_d / R_v
    sat_mixrat_liq_lookup = ep * ( esatv / ( p_in_Pa - esatv ) )
#endif

    end if

    return
  end function sat_mixrat_liq_lookup

!-----------------------------------------------------------------
  elemental function sat_vapor_press_liq( T_in_K ) result ( esat )

! Description:
!   Computes SVP for water vapor. Calls one of the other functions
!   that calculate an approximation to SVP.

! References:
!   None

    use model_flags, only: &
      saturation_formula, & ! Variable
      saturation_bolton, &
      saturation_gfdl, &
      saturation_flatau

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: T_in_K     ! Temperature                          [K]

    ! Output Variables
    real( kind = core_rknd ) :: esat      ! Saturation Vapor Pressure over Water [Pa]

    ! Undefined approximation
    esat = -99999.999_core_rknd

    ! Saturation Vapor Pressure, esat, can be found to be approximated
    ! in many different ways.
    select case ( saturation_formula )
    case ( saturation_bolton )
      ! Using the Bolton 1980 approximations for SVP over vapor
      esat = sat_vapor_press_liq_bolton( T_in_K )

    case ( saturation_flatau )
      ! Using the Flatau, et al. polynomial approximation for SVP over vapor
      esat = sat_vapor_press_liq_flatau( T_in_K )

! ---> h1g
    case ( saturation_gfdl )
      ! Using GFDL polynomial approximation for SVP with respect to liquid
      esat = sat_vapor_press_liq_gfdl( T_in_K )
! <--- h1g

      ! Add new cases after this

    end select

    return

  end function sat_vapor_press_liq

!------------------------------------------------------------------------
  elemental function sat_vapor_press_liq_lookup( T_in_K ) result ( esat )

! Description:
!   Computes SVP for water vapor, using a lookup table.
!
!   The lookup table was constructed using the Flatau approximation.

! References:
!   ``Polynomial Fits to Saturation Vapor Pressure'' Falatau, Walko,
!     and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
!     pp. 1507--1513
!------------------------------------------------------------------------

    implicit none

    ! External
    intrinsic :: max, min, int, anint

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: T_in_K   ! Temperature   [K]

    ! Output Variables
    real( kind = core_rknd ) :: esat  ! Saturation vapor pressure over water [Pa]

    ! Local Variables
    integer :: T_in_K_int

    ! ---- Begin Code ----

    T_in_K_int = int( anint( T_in_K ) )

    ! Since this approximation is only good out to -85 degrees Celsius we
    ! truncate the result here
    T_in_K_int = min( max( T_in_K_int, 188 ), 343 )

    ! Use the lookup table to determine the saturation vapor pressure.
    esat = svp_liq_lookup_table( T_in_K_int )

    return
  end function sat_vapor_press_liq_lookup

!------------------------------------------------------------------------
  elemental function sat_vapor_press_liq_flatau( T_in_K ) result ( esat )

! Description:
!   Computes SVP for water vapor.

! References:
!   ``Polynomial Fits to Saturation Vapor Pressure'' Falatau, Walko,
!     and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
!     pp. 1507--1513
!------------------------------------------------------------------------

    use constants_clubb, only: T_freeze_K

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Constant parameters

    ! Relative error norm expansion (-50 to 50 deg_C) from
    ! Table 3 of pp. 1510 of Flatau et al. 1992 (Water Vapor)
    ! (The 100 coefficient converts from mb to Pa)
!   real, dimension(7), parameter :: a = & 
!   100.* (/ 6.11176750,      0.443986062,     0.143053301E-01, & 
!            0.265027242E-03, 0.302246994E-05, 0.203886313E-07, & 
!            0.638780966E-10 /)

    ! Relative error norm expansion (-85 to 70 deg_C) from
    ! Table 4 of pp. 1511 of Flatau et al.
    real( kind = core_rknd ), dimension(9), parameter :: a = & 
    100._core_rknd * &
             (/ 6.11583699_core_rknd,      0.444606896_core_rknd,     0.143177157E-01_core_rknd, &
             0.264224321E-03_core_rknd, 0.299291081E-05_core_rknd, 0.203154182E-07_core_rknd, & 
             0.702620698E-10_core_rknd, 0.379534310E-13_core_rknd,-0.321582393E-15_core_rknd /)

    real( kind = core_rknd ), parameter :: min_T_in_C = -85._core_rknd ! [deg_C]

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: T_in_K   ! Temperature   [K]

    ! Output Variables
    real( kind = core_rknd ) :: esat  ! Saturation vapor pressure over water [Pa]

    ! Local Variables
    real( kind = core_rknd ) :: T_in_C
!   integer :: i ! Loop index

    ! ---- Begin Code ----

    ! Determine deg K - 273.15
    T_in_C = T_in_K - T_freeze_K

    ! Since this approximation is only good out to -85 degrees Celsius we
    ! truncate the result here (Flatau, et al. 1992)
    T_in_C = max( T_in_C, min_T_in_C )

    ! Polynomial approx. (Flatau, et al. 1992)

    ! This is the generalized formula but is not computationally efficient.
    ! Based on Wexler's expressions(2.1)-(2.4) (See Flatau et al. p 1508)
    ! e_{sat} = a_1 + a_2 ( T - T_0 ) + ... + a_{n+1} ( T - T_0 )^n

!   esat = a(1)

!   do i = 2, size( a ) , 1
!     esat = esat + a(i) * ( T_in_C )**(i-1)
!   end do

    ! The 8th order polynomial fit.  When running deep 
    ! convective cases I noticed that absolute temperature often dips below
    ! -50 deg_C at higher altitudes, where the 6th order approximation is
    ! not accurate.  -dschanen 20 Nov 2008
    esat = a(1) + T_in_C*( a(2) + T_in_C*( a(3) + T_in_C*( a(4) + T_in_C &
    *( a(5) + T_in_C*( a(6) + T_in_C*( a(7) + T_in_C*( a(8) + T_in_C*( a(9) ) ) ) ) ) ) ) )

    return
  end function sat_vapor_press_liq_flatau


!------------------------------------------------------------------------
  elemental function sat_vapor_press_liq_bolton( T_in_K ) result ( esat )
! Description:
!   Computes SVP for water vapor.
! References:
!   Bolton 1980
!------------------------------------------------------------------------

    use constants_clubb, only: T_freeze_K

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: exp

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: T_in_K   ! Temperature   [K]

    ! Output Variables
    real( kind = core_rknd ) :: esat  ! Saturation vapor pressure over water [Pa]

    ! (Bolton 1980) approx.
    ! Generally this more computationally expensive than the Flatau polnomial expansion
    esat = 611.2_core_rknd * exp( (17.67_core_rknd*(T_in_K-T_freeze_K)) / &
      (T_in_K-29.65_core_rknd) ) ! Known magic number

    return
  end function sat_vapor_press_liq_bolton


! ---> h1g, 2010-06-16
!------------------------------------------------------------------------
  elemental function sat_vapor_press_liq_gfdl( T_in_K ) result ( esat )
! Description:
! copy from "GFDL polysvp.F90" 
!  Compute saturation vapor pressure with respect to liquid  by using 
! function from Goff and Gatch (1946)

!  Polysvp returned in units of pa.
!  T_in_K  is input in units of K.
!------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: T_in_K   ! Temperature   [K]

    ! Output Variables
    real( kind = core_rknd ) :: esat  ! Saturation vapor pressure over water [Pa]

    ! Goff Gatch equation, uncertain below -70 C
      
         esat = 10._core_rknd**(-7.90298_core_rknd*(373.16_core_rknd/T_in_K-1._core_rknd)+ &
             5.02808_core_rknd*log10(373.16_core_rknd/T_in_K)- &
             1.3816e-7_core_rknd*(10._core_rknd**(11.344_core_rknd &
               *(1._core_rknd-T_in_K/373.16_core_rknd))-1._core_rknd)+ &
             8.1328e-3_core_rknd*(10._core_rknd**(-3.49149_core_rknd &
               *(373.16_core_rknd/T_in_K-1._core_rknd))-1._core_rknd)+ &
             log10(1013.246_core_rknd))*100._core_rknd ! Known magic number

    return
  end function sat_vapor_press_liq_gfdl
! <--- h1g, 2010-06-16

!------------------------------------------------------------------------
  elemental real( kind = core_rknd ) function sat_mixrat_ice( p_in_Pa, T_in_K )

! Description:
!   Used to compute the saturation mixing ratio of ice.

! References:
!   Formula from Emanuel 1994, 4.4.15
!-------------------------------------------------------------------------

    use constants_clubb, only: & 
        ep ! Variable(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: trim

    ! Input Variables

    real( kind = core_rknd ), intent(in) :: &
      p_in_Pa, &          ! Pressure [Pa]
      T_in_K              ! Temperature [K]

    ! Local Variables

    real( kind = core_rknd ) :: esat_ice

    ! --- Begin Code ---

    ! Determine the SVP for the given temperature
    esat_ice = sat_vapor_press_ice( T_in_K )

    ! If esat_ice exceeds the air pressure, then assume esat_ice~=0.5*pressure 
    !   and set rsat = ep = 0.622
    if ( p_in_Pa-esat_ice < 1.0_core_rknd ) then
      sat_mixrat_ice = ep
    else

#ifdef GFDL
    ! GFDL uses specific humidity
    ! Formula for Saturation Specific Humidity
     if( I_sat_sphum )  then   ! h1g, 2010-06-18 begin mod
           sat_mixrat_ice = ep * ( esat_ice / ( p_in_Pa - (1.0_core_rknd-ep) * esat_ice ) )
     else
           sat_mixrat_ice  = ep * ( esat_ice / ( p_in_Pa - esat_ice ) )
     endif                     ! h1g, 2010-06-18 end mod
#else
    ! Formula for Saturation Mixing Ratio:
    !
    ! rs = (epsilon) * [ esat / ( p - esat ) ];
    ! where epsilon = R_d / R_v

    sat_mixrat_ice = ep * ( esat_ice / ( p_in_Pa - esat_ice ) )
#endif

    end if

    return
  end function sat_mixrat_ice

!------------------------------------------------------------------------
  elemental function sat_vapor_press_ice( T_in_K ) result ( esat_ice )
!
! Description:
!   Computes SVP for ice, using one of the various approximations.
!
! References:
!   None
!------------------------------------------------------------------------
 
    use model_flags, only: &
      saturation_formula, & ! Variable(s)
      saturation_bolton, &
      saturation_gfdl, &
      saturation_flatau

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variable
    real( kind = core_rknd ), intent(in) :: &
      T_in_K      ! Temperature     [K]

    ! Output Variable
    real( kind = core_rknd ) :: esat_ice    ! Saturation Vapor Pressure over Ice [Pa]

    ! Undefined approximation
    esat_ice = -99999.999_core_rknd

    select case ( saturation_formula )
    case ( saturation_bolton )
      ! Using the Bolton 1980 approximations for SVP over ice
      esat_ice = sat_vapor_press_ice_bolton( T_in_K )

    case ( saturation_flatau )
      ! Using the Flatau, et al. polynomial approximation for SVP over ice
      esat_ice = sat_vapor_press_ice_flatau( T_in_K )

! ---> h1g, 2010-06-16
    case ( saturation_gfdl )
      ! Using GFDL polynomial approximation for SVP with respect to ice
      esat_ice = sat_vapor_press_ice_gfdl( T_in_K )
! <--- h1g, 2010-06-16

      ! Add new cases after this

    end select

    return

  end function sat_vapor_press_ice

!------------------------------------------------------------------------
  elemental function sat_vapor_press_ice_flatau( T_in_K ) result ( esati )
!
! Description:
!   Computes SVP for ice.
!
! References:
!   ``Polynomial Fits to Saturation Vapor Pressure'' Falatau, Walko,
!     and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
!     pp. 1507--1513
!------------------------------------------------------------------------
    use constants_clubb, only: T_freeze_K

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: max

    ! Relative error norm expansion (-90 to 0 deg_C) from
    ! Table 4 of pp. 1511 of Flatau et al. 1992 (Ice)
    real( kind = core_rknd ), dimension(9), parameter :: a = & 
    100._core_rknd * (/ 6.09868993_core_rknd, 0.499320233_core_rknd, 0.184672631E-01_core_rknd, &
              0.402737184E-03_core_rknd, 0.565392987E-05_core_rknd, 0.521693933E-07_core_rknd, &
              0.307839583E-09_core_rknd, 0.105785160E-11_core_rknd, 0.161444444E-14_core_rknd /)

    real( kind = core_rknd ), parameter :: min_T_in_C = -90._core_rknd ! [deg_C]


    ! Input Variables
    real( kind = core_rknd ), intent(in) :: T_in_K   ! Temperature   [deg_K]

    ! Output Variables
    real( kind = core_rknd ) :: esati  ! Saturation vapor pressure over ice [Pa]

    ! Local Variables
    real( kind = core_rknd ) :: T_in_C ! Temperature [deg_C]
!   integer :: i

    ! ---- Begin Code ----

    ! Determine deg K - 273.15
    T_in_C = T_in_K - T_freeze_K

    ! Since this approximation is only good out to -90 degrees Celsius we
    ! truncate the result here (Flatau, et al. 1992)
    T_in_C = max( T_in_C, min_T_in_C )

    ! Polynomial approx. (Flatau, et al. 1992)
!   esati = a(1)

!   do i = 2, size( a ), 1
!     esati = esati + a(i) * ( T_in_C )**(i-1)
!   end do

    esati = a(1) + T_in_C*( a(2) + T_in_C*( a(3) + T_in_C*( a(4) + T_in_C &
    *( a(5) + T_in_C*( a(6) + T_in_C*( a(7) + T_in_C*( a(8) + T_in_C*( a(9) ) ) ) ) ) ) ) )

    return

  end function sat_vapor_press_ice_flatau

!------------------------------------------------------------------------
  elemental function sat_vapor_press_ice_bolton( T_in_K ) result ( esati )
!
! Description:
!   Computes SVP for ice.
!
! References:
!   Bolton 1980
!------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: exp, log

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: T_in_K   ! Temperature   [K]

    ! Output Variables
    real( kind = core_rknd ) :: esati  ! Saturation vapor pressure over ice [Pa]

    ! Exponential approx.
    esati = 100.0_core_rknd * exp( 23.33086_core_rknd - &
      (6111.72784_core_rknd/T_in_K) + (0.15215_core_rknd*log( T_in_K )) )

    return

  end function sat_vapor_press_ice_bolton


! ---> h1g, 2010-06-16
!------------------------------------------------------------------------
  elemental function sat_vapor_press_ice_gfdl( T_in_K ) result ( esati )
! Description:
! copy from "GFDL polysvp.F90" 
!  Compute saturation vapor pressure with respect to liquid  by using 
! function from Goff and Gatch (1946)
! 
!  Polysvp returned in units of pa.
!  T_in_K is input in units of K.
!------------------------------------------------------------------------
 
    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: T_in_K   ! Temperature   [K]

    ! Output Variables
    real( kind = core_rknd ) :: esati  ! Saturation vapor pressure over ice [Pa]

    ! Goff Gatch equation (good down to -100 C)

    esati = 10._core_rknd**(-9.09718_core_rknd* &
            (273.16_core_rknd/T_in_k-1._core_rknd)-3.56654_core_rknd* &
          log10(273.16_core_rknd/T_in_k)+0.876793_core_rknd* &
            (1._core_rknd-T_in_k/273.16_core_rknd)+ &
          log10(6.1071_core_rknd))*100._core_rknd ! Known magic number

    return

  end function sat_vapor_press_ice_gfdl
! <--- h1g, 2010-06-16

!-------------------------------------------------------------------------
  function rcm_sat_adj( thlm, rtm, p_in_Pa, exner ) result ( rcm )

    ! Description:
    !
    !   This function uses an iterative method to find the value of rcm
    !   from an initial profile that has saturation at some point.
    !
    ! References:
    !   None
    !-------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use constants_clubb, only: & 
        Cp,            & ! Variable(s)
        Lv,            &
        zero_threshold

    implicit none

    ! Local Constant(s)
    real( kind = core_rknd ), parameter :: &
      tolerance = 0.001_core_rknd ! Tolerance on theta calculation [K]

    integer, parameter :: &
      itermax = 1000000 ! Maximum interations

    ! External
    intrinsic :: max, abs

    ! Input Variable(s)
    real( kind = core_rknd ), intent(in) :: &
      thlm,    & ! Liquid Water Potential Temperature [K]
      rtm,     & ! Total Water Mixing Ratio       [kg/kg]
      p_in_Pa, & ! Pressure                          [Pa]
      exner      ! Exner function                     [-]

    ! Output Variable(s)
    real( kind = core_rknd ) :: rcm ! Cloud water mixing ratio      [kg/kg]

    ! Local Variable(s)
    real( kind = core_rknd ) :: &
      theta, answer, too_low, too_high ! [K]

    integer :: iteration

    ! ----- Begin Code -----

    ! Default initialization
    theta = thlm
    too_high = 0.0_core_rknd
    too_low  = 0.0_core_rknd

    do iteration = 1, itermax, 1

      answer = &
      theta - (Lv/(Cp*exner)) &
             *(MAX( rtm - sat_mixrat_liq(p_in_Pa,theta*exner), zero_threshold ))

      if ( ABS(answer - thlm) <= tolerance ) then
        exit
      else if ( answer - thlm > tolerance ) then
        too_high = theta
      else if ( thlm - answer > tolerance ) THEN
        too_low = theta
      end if

      ! For the first timestep, be sure to set a "too_high"
      ! that is "way too high."
      if ( iteration == 1 ) then
        too_high = theta + 20.0_core_rknd
      end if

      theta = (too_low + too_high)/2.0_core_rknd

    end do ! 1..itermax

    if ( iteration == itermax ) then
      ! Magic Eric Raut added to remove compiler warning (clearly this value is not used)
      rcm = 0.0_core_rknd
      
      stop "Error in rcm_sat_adj: could not determine rcm"
    else
      rcm = MAX( rtm - sat_mixrat_liq( p_in_Pa, theta*exner), zero_threshold )
      return
    end if

  end function rcm_sat_adj

end module saturation
