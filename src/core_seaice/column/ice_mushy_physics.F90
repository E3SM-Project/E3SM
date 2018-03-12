module ice_mushy_physics

  use ice_kinds_mod
  use ice_constants_colpkg, only: c0, c1, c2, c4, c8, c10, c1000, &
      p001, p01, p05, p1, p2, p5, pi, bignum, puny, ice_ref_salinity, &
      viscosity_dyn, rhow, rhoi, rhos, cp_ocn, cp_ice, Lfresh, gravit, &
      hs_min, ksno

  implicit none

  private
  public :: &
       conductivity_mush_array, &
       conductivity_snow_array, &
       enthalpy_snow, &
       enthalpy_brine, &
       enthalpy_mush, &
       enthalpy_mush_liquid_fraction, &
       enthalpy_of_melting, &
       temperature_snow, &
       temperature_mush, &
       temperature_mush_liquid_fraction, &
       liquidus_brine_salinity_mush, &
       liquidus_temperature_mush, &
       liquid_fraction, &
       density_brine
      
  !-----------------------------------------------------------------
  ! Constants for Liquidus relation from Assur (1958)
  !-----------------------------------------------------------------

  ! liquidus relation - higher temperature region
  real(kind=dbl_kind), parameter :: &
       az1_liq = -18.48_dbl_kind, &
       bz1_liq =    0.0_dbl_kind

  ! liquidus relation - lower temperature region
  real(kind=dbl_kind), parameter :: &
       az2_liq = -10.3085_dbl_kind, &
       bz2_liq =     62.4_dbl_kind

  ! liquidus break
  real(kind=dbl_kind), parameter :: &
       Tb_liq = -7.6362968855167352_dbl_kind, & ! temperature of liquidus break
       Sb_liq =  123.66702800276086_dbl_kind    ! salinity of liquidus break

  ! basic liquidus relation constants
  real(kind=dbl_kind), parameter :: &
       az1p_liq = az1_liq / c1000, &
       bz1p_liq = bz1_liq / c1000, &
       az2p_liq = az2_liq / c1000, &
       bz2p_liq = bz2_liq / c1000
  
  ! quadratic constants - higher temperature region
  real(kind=dbl_kind), parameter :: &
       AS1_liq = az1p_liq * (rhow * cp_ocn - rhoi * cp_ice)       , &
       AC1_liq = rhoi * cp_ice * az1_liq                          , & 
       BS1_liq = (c1 + bz1p_liq) * (rhow * cp_ocn - rhoi * cp_ice)  &
               + rhoi * Lfresh * az1p_liq                         , &
       BQ1_liq = -az1_liq                                         , &
       BC1_liq = rhoi * cp_ice * bz1_liq - rhoi * Lfresh * az1_liq, &
       CS1_liq = rhoi * Lfresh * (c1 + bz1p_liq)                  , &
       CQ1_liq = -bz1_liq                                         , &
       CC1_liq = -rhoi * Lfresh * bz1_liq
  
  ! quadratic constants - lower temperature region
  real(kind=dbl_kind), parameter :: &
       AS2_liq = az2p_liq * (rhow * cp_ocn - rhoi * cp_ice)       , &
       AC2_liq = rhoi * cp_ice * az2_liq                          , &
       BS2_liq = (c1 + bz2p_liq) * (rhow * cp_ocn - rhoi * cp_ice)  &
               + rhoi * Lfresh * az2p_liq                         , &
       BQ2_liq = -az2_liq                                         , &
       BC2_liq = rhoi * cp_ice * bz2_liq - rhoi * Lfresh * az2_liq, &
       CS2_liq = rhoi * Lfresh * (c1 + bz2p_liq)                  , &
       CQ2_liq = -bz2_liq                                         , &
       CC2_liq = -rhoi * Lfresh * bz2_liq
  
  ! break enthalpy constants
  real(kind=dbl_kind), parameter :: &
       D_liq = ((c1 + az1p_liq*Tb_liq + bz1p_liq) &
             / (       az1_liq*Tb_liq + bz1_liq)) &
             * ((cp_ocn*rhow - cp_ice*rhoi)*Tb_liq + Lfresh*rhoi), &
       E_liq = cp_ice*rhoi*Tb_liq - Lfresh*rhoi
  
  ! just fully melted enthapy constants
  real(kind=dbl_kind), parameter :: &
       F1_liq = (  -c1000 * cp_ocn * rhow) / az1_liq , &
       G1_liq =    -c1000                            , &
       H1_liq = (-bz1_liq * cp_ocn * rhow) / az1_liq , &
       F2_liq = (  -c1000 * cp_ocn * rhow) / az2_liq , &
       G2_liq =    -c1000                            , &
       H2_liq = (-bz2_liq * cp_ocn * rhow) / az2_liq
  
  ! warmer than fully melted constants
  real(kind=dbl_kind), parameter :: &
       I_liq = c1 / (cp_ocn * rhow)

  ! temperature to brine salinity
  real(kind=dbl_kind), parameter :: &
       J1_liq = bz1_liq / az1_liq         , &
       K1_liq = c1 / c1000                , &
       L1_liq = (c1 + bz1p_liq) / az1_liq , &
       J2_liq = bz2_liq  / az2_liq        , &
       K2_liq = c1 / c1000                , &
       L2_liq = (c1 + bz2p_liq) / az2_liq

  ! brine salinity to temperature
  real(kind=dbl_kind), parameter :: &
       M1_liq = az1_liq            , &
       N1_liq = -az1p_liq          , &
       O1_liq = -bz1_liq / az1_liq , &
       M2_liq = az2_liq            , &
       N2_liq = -az2p_liq          , &
       O2_liq = -bz2_liq / az2_liq

  !-----------------------------------------------------------------
  ! Other parameters
  !-----------------------------------------------------------------
  
  real(kind=dbl_kind), parameter :: &
       ki = 2.3_dbl_kind , & ! fresh ice conductivity (W m-1 K-1)
       kb = 0.5375_dbl_kind  ! brine conductivity (W m-1 K-1)

!=======================================================================

contains

!=======================================================================
! Physical Quantities
!=======================================================================

  subroutine conductivity_mush_array(nilyr, zqin, zSin, km)

    ! detemine the conductivity of the mush from enthalpy and salinity
    
    integer (kind=int_kind), intent(in) :: &
         nilyr   ! number of ice layers

    real(kind=dbl_kind), dimension(:), intent(in) :: &
         zqin, & ! ice layer enthalpy (J m-3) 
         zSin    ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         km      ! ice layer conductivity (W m-1 K-1)
    
    integer(kind=int_kind) :: &
         k       ! ice layer index

    real(kind=dbl_kind) :: Tmush
    
    do k = 1, nilyr
      
       Tmush = temperature_mush(zqin(k), zSin(k))
       
       km(k) = heat_conductivity(Tmush, zSin(k))
       
    enddo ! k

  end subroutine conductivity_mush_array

!=======================================================================

  function density_brine(Sbr) result(rho)
    
    ! density of brine from brine salinity

    real(kind=dbl_kind), intent(in) :: &
         Sbr ! brine salinity (ppt)

    real(kind=dbl_kind) :: &
         rho ! brine density (kg m-3)
    
    real(kind=dbl_kind), parameter :: &
         a = 1000.3_dbl_kind    , & ! zeroth empirical coefficient
         b = 0.78237_dbl_kind   , & ! linear empirical coefficient
         c = 2.8008e-4_dbl_kind     ! quadratic empirical coefficient
    
    rho = a + b * Sbr + c * Sbr**2
                
  end function density_brine

!=======================================================================
! Snow
!=======================================================================

  subroutine conductivity_snow_array(ks)

    ! heat conductivity of the snow

    real(kind=dbl_kind), dimension(:), intent(out) :: &
         ks ! snow layer conductivity (W m-1 K-1)

    ks = ksno

  end subroutine conductivity_snow_array

!=======================================================================
  
  function enthalpy_snow(zTsn) result(zqsn)
    
    ! enthalpy of snow from snow temperature

    real(kind=dbl_kind), intent(in) :: &
         zTsn ! snow layer temperature (C)

    real(kind=dbl_kind) :: &
         zqsn ! snow layer enthalpy (J m-3) 
    
    zqsn = -rhos * (-cp_ice * zTsn + Lfresh)
    
  end function enthalpy_snow

!=======================================================================
  
  function temperature_snow(zqsn) result(zTsn)
    
    ! temperature of snow from the snow enthalpy

    real(kind=dbl_kind), intent(in) :: &
         zqsn ! snow layer enthalpy (J m-3) 

    real(kind=dbl_kind) :: &
         zTsn ! snow layer temperature (C)
    
    real(kind=dbl_kind), parameter :: &
         A = c1 / (rhos * cp_ice) , &
         B = Lfresh / cp_ice
    
    zTsn = A * zqsn + B

  end function temperature_snow

!=======================================================================
! Mushy Layer Formulation - Assur (1958) liquidus
!=======================================================================

  function liquidus_brine_salinity_mush(zTin) result(Sbr)

    ! liquidus relation: equilibrium brine salinity as function of temperature
    ! based on empirical data from Assur (1958)

    real(kind=dbl_kind), intent(in) :: &
         zTin         ! ice layer temperature (C)

    real(kind=dbl_kind) :: &
         Sbr          ! ice brine salinity (ppt)

    real(kind=dbl_kind) :: &
         t_high   , & ! mask for high temperature liquidus region
         lsubzero     ! mask for sub-zero temperatures

    t_high   = merge(c1, c0, (zTin > Tb_liq))
    lsubzero = merge(c1, c0, (zTin <= c0))

    Sbr = ((zTin + J1_liq) / (K1_liq * zTin + L1_liq)) * t_high + &
          ((zTin + J2_liq) / (K2_liq * zTin + L2_liq)) * (c1 - t_high)

    Sbr = Sbr * lsubzero

  end function liquidus_brine_salinity_mush

!=======================================================================

  function liquidus_temperature_mush(Sbr) result(zTin)

    ! liquidus relation: equilibrium temperature as function of brine salinity
    ! based on empirical data from Assur (1958)

    real(kind=dbl_kind), intent(in) :: &
         Sbr    ! ice brine salinity (ppt)

    real(kind=dbl_kind) :: &
         zTin   ! ice layer temperature (C)

    real(kind=dbl_kind) :: &
         t_high ! mask for high temperature liquidus region

    t_high = merge(c1, c0, (Sbr <= Sb_liq))

    zTin = ((Sbr / (M1_liq + N1_liq * Sbr)) + O1_liq) * t_high + &
          ((Sbr / (M2_liq + N2_liq * Sbr)) + O2_liq) * (c1 - t_high)

  end function liquidus_temperature_mush

!=======================================================================

  function enthalpy_mush(zTin, zSin) result(zqin)

    ! enthalpy of mush from mush temperature and bulk salinity

    real(kind=dbl_kind), intent(in) :: &
         zTin, & ! ice layer temperature (C)
         zSin    ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         zqin    ! ice layer enthalpy (J m-3) 

    real(kind=dbl_kind) :: &
         phi     ! ice liquid fraction 

    phi = liquid_fraction(zTin, zSin)
    
    zqin = phi * (cp_ocn * rhow - cp_ice * rhoi) * zTin + &
           rhoi * cp_ice * zTin - (c1 - phi) * rhoi * Lfresh

  end function enthalpy_mush

!=======================================================================

  function enthalpy_mush_liquid_fraction(zTin, phi) result(zqin)

    ! enthalpy of mush from mush temperature and bulk salinity

    real(kind=dbl_kind), intent(in) :: &
         zTin, & ! ice layer temperature (C)
         phi     ! liquid fraction

    real(kind=dbl_kind) :: &
         zqin    ! ice layer enthalpy (J m-3) 

    zqin = phi * (cp_ocn * rhow - cp_ice * rhoi) * zTin + &
           rhoi * cp_ice * zTin - (c1 - phi) * rhoi * Lfresh

  end function enthalpy_mush_liquid_fraction

!=======================================================================

  function enthalpy_of_melting(zSin) result(qm)

    ! enthalpy of melting of mush
    ! energy needed to fully melt mush (T < 0)

    real(kind=dbl_kind), intent(in) :: &
         zSin ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         qm   ! melting ice enthalpy (J m-3)

    qm = cp_ocn * rhow * liquidus_temperature_mush(zSin)

  end function enthalpy_of_melting

!=======================================================================

  function enthalpy_brine(zTin) result(qbr)

    ! enthalpy of brine (fully liquid)

    real(kind=dbl_kind), intent(in) :: &
         zTin ! ice layer temperature (C)

    real(kind=dbl_kind) :: &
         qbr  ! brine enthalpy (J m-3)

    qbr = cp_ocn * rhow * zTin

  end function enthalpy_brine

!=======================================================================

  function temperature_mush(zqin, zSin) result(zTin)

    ! temperature of mush from mush enthalpy

    real(kind=dbl_kind), intent(in) :: &
         zqin   , & ! ice enthalpy (J m-3) 
         zSin       ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         zTin       ! ice layer temperature (C)

    real(kind=dbl_kind) :: &
         qb     , & ! liquidus break enthalpy
         q0     , & ! fully melted enthalpy
         A      , & ! quadratic equation A parameter
         B      , & ! quadratic equation B parameter
         C      , & ! quadratic equation C parameter
         S_low  , & ! mask for salinity less than the liquidus break salinity
         t_high , & ! mask for high temperature liquidus region
         t_low  , & ! mask for low temperature liquidus region
         q_melt     ! mask for all mush melted

    ! just melted enthalpy
    S_low = merge(c1, c0, (zSin < Sb_liq))
    q0 = ((F1_liq * zSin) / (G1_liq + zSin) + H1_liq) * S_low + &
         ((F2_liq * zSin) / (G2_liq + zSin) + H2_liq) * (c1 - S_low)
    q_melt = merge(c1, c0, (zqin > q0))

    ! break enthalpy
    qb = D_liq * zSin + E_liq
    t_high = merge(c1, c0, (zqin > qb))
    t_low = c1 - t_high

    ! quadratic values
    A = (AS1_liq * zSin                 + AC1_liq) * t_high + &
        (AS2_liq * zSin                 + AC2_liq) * t_low

    B = (BS1_liq * zSin + BQ1_liq * zqin + BC1_liq) * t_high + &
        (BS2_liq * zSin + BQ2_liq * zqin + BC2_liq) * t_low

    C = (CS1_liq * zSin + CQ1_liq * zqin + CC1_liq) * t_high + &
        (CS2_liq * zSin + CQ2_liq * zqin + CC2_liq) * t_low

    zTin = (-B + sqrt(max(B**2 - c4 * A * C,puny))) / (c2 * A)

    ! change T if all melted
    zTin = q_melt * zqin * I_liq + (c1 - q_melt) * zTin

  end function temperature_mush

!=======================================================================

  function temperature_mush_liquid_fraction(zqin, phi) result(zTin)

    ! temperature of mush from mush enthalpy

    real(kind=dbl_kind), intent(in) :: &
         zqin   , & ! ice enthalpy (J m-3) 
         phi        ! liquid fraction

    real(kind=dbl_kind) :: &
         zTin       ! ice layer temperature (C)

    zTin = (zqin + (c1 - phi) * rhoi * Lfresh) / &
          (phi * (cp_ocn * rhow - cp_ice * rhoi) + rhoi * cp_ice)

  end function temperature_mush_liquid_fraction

!=======================================================================

  function heat_conductivity(zTin, zSin) result(km)
    
    ! msuh heat conductivity from mush temperature and bulk salinity

    real(kind=dbl_kind), intent(in) :: &
         zTin              , & ! ice layer temperature (C)
         zSin                  ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         km                    ! ice layer conductivity (W m-1 K-1)
    
    real(kind=dbl_kind) :: &
         phi                   ! liquid fraction

    phi = liquid_fraction(zTin, zSin)

    km = phi * (kb - ki) + ki

  end function heat_conductivity

  !=======================================================================

  function liquid_fraction(zTin, zSin) result(phi)

    ! liquid fraction of mush from mush temperature and bulk salinity

    real(kind=dbl_kind), intent(in) :: &
         zTin, & ! ice layer temperature (C)
         zSin    ! ice layer bulk salinity (ppt)

    real(kind=dbl_kind) :: &
         phi , & ! liquid fraction
         Sbr     ! brine salinity (ppt)

    Sbr = max(liquidus_brine_salinity_mush(zTin),puny)
    phi = zSin / max(Sbr, zSin)

  end function liquid_fraction

!=======================================================================

end module ice_mushy_physics


