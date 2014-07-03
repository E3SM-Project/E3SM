MODULE co2calc

  !-----------------------------------------------------------------------------
  !   based upon OCMIP2 co2calc
  !
  !   CVS:$Id: co2calc.F90 941 2006-05-12 21:36:48Z klindsay $
  !   CVS:$Name$
  !-----------------------------------------------------------------------------

  USE constants
  USE blocks, ONLY : nx_block, ny_block, block, get_block
  USE domain, ONLY : blocks_clinic
  USE domain_size, ONLY : max_blocks_clinic
  USE kinds_mod
  USE state_mod, ONLY : ref_pressure
  USE io_types, ONLY : stdout
  USE time_management, ONLY : nsteps_run

#ifdef CCSMCOUPLED
   !*** ccsm
  USE shr_vmath_mod
#endif

  IMPLICIT NONE

  !-----------------------------------------------------------------------------
  !   public/private declarations
  !-----------------------------------------------------------------------------

  PRIVATE
  PUBLIC :: co2calc_row, comp_CO3terms, comp_co3_sat_vals

  !-----------------------------------------------------------------------------
  !   module parameters
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !   The current setting of xacc, a tolerance critera, will result in co2star
  !   being accurate to 3 significant figures (xx.y). Making xacc bigger will
  !   result in faster convergence also, but this is not recommended (xacc of
  !   10**-9 drops precision to 2 significant figures).
  !-----------------------------------------------------------------------------

  REAL(KIND=r8), PARAMETER :: xacc = 1e-10_r8
  INTEGER(KIND=int_kind), PARAMETER :: max_bracket_grow_it = 3
  INTEGER(KIND=int_kind), PARAMETER :: maxit = 100

  REAL(KIND=r8), PARAMETER :: salt_min = 0.1_r8
  REAL(KIND=r8), PARAMETER :: dic_min  = salt_min / 35.0_r8 * 1944.0_r8
  REAL(KIND=r8), PARAMETER :: alk_min  = salt_min / 35.0_r8 * 2225.0_r8

  !-----------------------------------------------------------------------------
  !   declarations for function coefficients & species concentrations
  !-----------------------------------------------------------------------------

  REAL(KIND=r8), DIMENSION(nx_block,max_blocks_clinic) :: &
       kw, kb, ks, kf, k1p, k2p, k3p, ksi, &
       bt, st, ft, dic, ta, pt, sit

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  SUBROUTINE co2calc_row(iblock, j, mask, locmip_k1_k2_bug_fix, lcomp_co3_coeffs, &
       temp, salt, dic_in, ta_in, pt_in, sit_in, phlo, phhi, ph, xco2_in, atmpres, &
       co2star, dco2star, pCO2surf, dpco2)

    !---------------------------------------------------------------------------
    !   SUBROUTINE co2calc_row
    !
    !   PURPOSE : Calculate delta co2*, etc. from total alkalinity, total CO2,
    !             temp, salinity (s), etc.
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: iblock, j
    LOGICAL(KIND=log_kind), DIMENSION(nx_block), INTENT(IN) :: mask
    LOGICAL(KIND=log_kind), INTENT(IN) :: locmip_k1_k2_bug_fix
    LOGICAL(KIND=log_kind), INTENT(IN) :: lcomp_co3_coeffs
    REAL(KIND=r8), DIMENSION(nx_block), INTENT(IN) :: &
         temp,     & ! temperature (degrees C)
         salt,     & ! salinity (PSU)
         dic_in,   & ! total inorganic carbon (nmol/cm^3)
         ta_in,    & ! total alkalinity (neq/cm^3)
         pt_in,    & ! inorganic phosphate (nmol/cm^3)
         sit_in,   & ! inorganic silicate (nmol/cm^3)
         xco2_in,  & ! atmospheric mole fraction CO2 in dry air (ppmv)
         atmpres     ! atmospheric pressure (atmosphere)

    !---------------------------------------------------------------------------
    !   input/output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=r8), DIMENSION(nx_block), INTENT(INOUT) :: &
         phlo,     & ! lower limit of pH range
         phhi        ! upper limit of pH range

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=r8), DIMENSION(nx_block), INTENT(OUT) :: &
         ph,       & ! computed ph values, for initial guess on next time step
         co2star,  & ! CO2*water (nmol/cm^3)
         dco2star, & ! delta CO2 (nmol/cm^3)
         pco2surf, & ! oceanic pCO2 (ppmv)
         dpco2       ! Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind) :: i
    INTEGER(KIND=int_kind) :: k

    REAL(KIND=r8) :: &
         mass_to_vol,  & ! (mol/kg) -> (mmol/m^3)
         vol_to_mass,  & ! (mmol/m^3) -> (mol/kg)
         co2starair,   & ! co2star saturation
         htotal2

    REAL(KIND=r8), DIMENSION(nx_block) :: &
         xco2,         & ! atmospheric CO2 (atm)
         htotal,       & ! free concentration of H ion
         k0,k1,k2,     & ! equilibrium constants for CO2 species
         ff              ! fugacity of CO2

    !---------------------------------------------------------------------------
    !   check for existence of ocean points
    !---------------------------------------------------------------------------

    IF (COUNT(mask) == 0) THEN
       ph          = c0
       co2star     = c0
       dco2star    = c0
       pCO2surf    = c0
       dpCO2       = c0
       RETURN
    END IF

    !---------------------------------------------------------------------------
    !   set unit conversion factors
    !---------------------------------------------------------------------------

    mass_to_vol = 1e6_r8 * rho_sw
    vol_to_mass = c1 / mass_to_vol

    k = 1

    !---------------------------------------------------------------------------
    !   compute thermodynamic CO3 coefficients
    !---------------------------------------------------------------------------

    IF (lcomp_co3_coeffs) THEN
       CALL comp_co3_coeffs(iblock, k, mask, temp, salt, k0, k1, k2, ff, &
                            k1_k2_pH_tot=locmip_k1_k2_bug_fix)
    END IF

    !---------------------------------------------------------------------------
    !   compute htotal
    !---------------------------------------------------------------------------

    CALL comp_htotal(iblock, j, k, mask, temp, dic_in, ta_in, pt_in, sit_in, &
                     k1, k2, phlo, phhi, htotal)

    !---------------------------------------------------------------------------
    !   convert xco2 from uatm to atm
    !---------------------------------------------------------------------------

    DO i = 1,nx_block
       IF (mask(i)) THEN
          xco2(i) = xco2_in(i) * 1e-6_r8
       END IF ! if mask
    END DO ! i loop

    !---------------------------------------------------------------------------
    !   Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
    !   ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
    !
    !   Compute co2starair
    !---------------------------------------------------------------------------

    DO i = 1,nx_block
       IF (mask(i)) THEN

          htotal2 = htotal(i) ** 2
          co2star(i) = dic(i,iblock) * htotal2 / &
               (htotal2 + k1(i) * htotal(i) + k1(i) * k2(i))
          co2starair = xco2(i) * ff(i) * atmpres(i)
          dco2star(i) = co2starair - co2star(i)
          ph(i) = -LOG10(htotal(i))

          !---------------------------------------------------------------------
          !   Add two output arguments for storing pCO2surf
          !   Should we be using K0 or ff for the solubility here?
          !---------------------------------------------------------------------

          pCO2surf(i) = co2star(i) / ff(i)
          dpCO2(i)    = pCO2surf(i) - xco2(i) * atmpres(i)

          !---------------------------------------------------------------------
          !   Convert units of output arguments
          !   Note: pCO2surf and dpCO2 are calculated in atm above.
          !---------------------------------------------------------------------

          co2star(i)  = co2star(i) * mass_to_vol
          dco2star(i) = dco2star(i) * mass_to_vol

          pCO2surf(i) = pCO2surf(i) * 1e6_r8
          dpCO2(i)    = dpCO2(i) * 1e6_r8

       ELSE ! if mask

          ph(i)       = c0
          co2star(i)  = c0
          dco2star(i) = c0
          pCO2surf(i) = c0
          dpCO2(i)    = c0

       END IF ! if mask
    END DO ! i loop

  END SUBROUTINE co2calc_row

  !*****************************************************************************

  SUBROUTINE comp_CO3terms(iblock, j, k, mask, lcomp_co3_coeffs, temp, salt, &
       dic_in, ta_in, pt_in, sit_in, phlo, phhi, ph, H2CO3, HCO3, CO3)

    !---------------------------------------------------------------------------
    !   SUBROUTINE comp_CO3terms
    !
    !   PURPOSE : Calculate H2CO3, HCO3, CO3 from
    !             total alkalinity, total CO2, temp, salinity (s), etc.
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: iblock
    INTEGER(KIND=int_kind), INTENT(IN) :: j, k
    LOGICAL(KIND=log_kind), DIMENSION(nx_block), INTENT(IN) :: mask
    LOGICAL(KIND=log_kind), INTENT(IN) :: lcomp_co3_coeffs
    REAL(KIND=r8), DIMENSION(nx_block), INTENT(IN) :: &
         temp,     & ! temperature (degrees C)
         salt,     & ! salinity (PSU)
         dic_in,   & ! total inorganic carbon (nmol/cm^3)
         ta_in,    & ! total alkalinity (neq/cm^3)
         pt_in,    & ! inorganic phosphate (nmol/cm^3)
         sit_in      ! inorganic silicate (nmol/cm^3)

    !---------------------------------------------------------------------------
    !   input/output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=r8), DIMENSION(nx_block), INTENT(INOUT) :: &
         phlo,     & ! lower limit of pH range
         phhi        ! upper limit of pH range

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=r8), DIMENSION(nx_block), INTENT(OUT) :: &
         pH,         & ! computed ph values, for initial guess on next time step
         H2CO3,      & ! Carbonic Acid Concentration
         HCO3,       & ! Bicarbonate Ion Concentration
         CO3           ! Carbonate Ion Concentration

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind) :: i

    REAL(KIND=r8) :: &
         mass_to_vol,  & ! (mol/kg) -> (mmol/m^3)
         vol_to_mass,  & ! (mmol/m^3) -> (mol/kg)
         htotal2, denom

    REAL(KIND=r8), DIMENSION(nx_block) :: &
         htotal,       & ! free concentration of H ion
         k0,k1,k2,     & ! equilibrium constants for CO2 species
         ff              ! fugacity of CO2

    !---------------------------------------------------------------------------
    !   check for existence of ocean points
    !---------------------------------------------------------------------------

    IF (COUNT(mask) == 0) THEN
       ph         = c0
       H2CO3      = c0
       HCO3       = c0
       CO3        = c0
       RETURN
    END IF

    !---------------------------------------------------------------------------
    !   set unit conversion factors
    !---------------------------------------------------------------------------

    mass_to_vol = 1e6_r8 * rho_sw
    vol_to_mass = c1 / mass_to_vol

    !------------------------------------------------------------------------
    !   compute thermodynamic CO3 coefficients
    !------------------------------------------------------------------------

    IF (lcomp_co3_coeffs) THEN
       CALL comp_co3_coeffs(iblock, k, mask, temp, salt, k0, k1, k2, ff, k1_k2_pH_tot=.true.)
    END IF

    !------------------------------------------------------------------------
    !   compute htotal
    !------------------------------------------------------------------------

    CALL comp_htotal(iblock, j, k, mask, temp, dic_in, &
                     ta_in, pt_in, sit_in, k1, k2, &
                     phlo, phhi, htotal)

    !------------------------------------------------------------------------
    !   Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
    !   ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49-51)
    !------------------------------------------------------------------------

    DO i = 1,nx_block
       IF (mask(i)) THEN

          htotal2  = htotal(i) ** 2
          denom    = c1 / (htotal2 + k1(i) * htotal(i) + k1(i) * k2(i))
          H2CO3(i) = dic(i,iblock) * htotal2 * denom
          HCO3(i)  = dic(i,iblock) * k1(i) * htotal(i) * denom
          CO3(i)   = dic(i,iblock) * k1(i) * k2(i) * denom
          ph(i)    = -LOG10(htotal(i))

          !------------------------------------------------------------------
          !   Convert units of output arguments
          !------------------------------------------------------------------

          H2CO3(i) = H2CO3(i) * mass_to_vol
          HCO3(i)  = HCO3(i) * mass_to_vol
          CO3(i)   = CO3(i) * mass_to_vol

       ELSE ! if mask

          ph(i)    = c0
          H2CO3(i) = c0
          HCO3(i)  = c0
          CO3(i)   = c0

       END IF ! if mask
    END DO ! i loop

  END SUBROUTINE comp_CO3terms

  !*****************************************************************************

  SUBROUTINE comp_co3_coeffs(iblock, k, mask, temp, salt, k0, k1, k2, ff, k1_k2_pH_tot)

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: iblock
    INTEGER(KIND=int_kind), INTENT(IN) :: k
    LOGICAL(KIND=log_kind), DIMENSION(nx_block), INTENT(IN) :: mask
    REAL(KIND=r8), DIMENSION(nx_block), INTENT(IN) :: &
         temp,     & ! temperature (degrees C)
         salt        ! salinity (PSU)
    LOGICAL(KIND=log_kind), INTENT(IN) :: k1_k2_pH_tot

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=r8), DIMENSION(nx_block), INTENT(OUT) :: &
         k0,k1,k2,     & ! equilibrium constants for CO2 species
         ff              ! fugacity of CO2

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind) :: i

    REAL(KIND=r8) :: &
         press_bar       ! pressure at level k [bars]

    REAL(KIND=r8), DIMENSION(nx_block) :: &
         salt_lim,     & ! bounded salt
         tk,           & ! temperature (K)
         is,           & ! ionic strength
         scl,          & ! chlorinity
         tk100, tk1002, invtk, dlogtk, is2, sqrtis, &
         s2, sqrts, s15, invRtk, arg, &
         deltaV,Kappa,lnKfac,Kfac, & ! pressure correction terms
         log_1_m_1p005em3_s, &
         log_1_p_tot_sulfate_div_ks

    !---------------------------------------------------------------------------

    press_bar = ref_pressure(k)

    !---------------------------------------------------------------------------
    !   Calculate all constants needed to convert between various
    !   measured carbon species. References for each equation are
    !   noted in the code.  Once calculated, the constants are stored
    !   and passed in the common block "const". The original version
    !   of this code was based on the code by Dickson in Version 2 of
    !   "Handbook of Methods for the Analysis of the Various Parameters
    !   of the Carbon Dioxide System in Seawater", DOE, 1994 (SOP No. 3,
    !   p25-26).
    !   Derive simple terms used more than once
    !---------------------------------------------------------------------------

    salt_lim = max(salt,salt_min)
    tk       = T0_Kelvin + temp
    tk100    = tk * 1e-2_r8
    tk1002   = tk100 * tk100
    invtk    = c1 / tk
#ifdef CCSMCOUPLED
    CALL shr_vmath_log(tk, dlogtk, nx_block)
#else
    dlogtk   = LOG(tk)
#endif
    invRtk   = (c1 / 83.1451_r8) * invtk

    is       = 19.924_r8 * salt_lim / (c1000 - 1.005_r8 * salt_lim)
    is2      = is * is
#ifdef CCSMCOUPLED
    CALL shr_vmath_sqrt(is, sqrtis, nx_block)
    CALL shr_vmath_sqrt(salt_lim, sqrts, nx_block)
#else
    sqrtis   = SQRT(is)
    sqrts    = SQRT(salt_lim)
#endif
    s2       = salt_lim * salt_lim
    scl      = salt_lim / 1.80655_r8

    arg = c1 - 0.001005_r8 * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_log(arg, log_1_m_1p005em3_s, nx_block)
#else
    log_1_m_1p005em3_s = LOG(arg)
#endif

    !---------------------------------------------------------------------------
    !   f = k0(1-pH2O)*correction term for non-ideality
    !   Weiss & Price (1980, Mar. Chem., 8, 347-359;
    !                 Eq 13 with table 6 values)
    !---------------------------------------------------------------------------

    arg = -162.8301_r8 + 218.2968_r8 / tk100 + &
          90.9241_r8 * (dlogtk + LOG(1e-2_r8)) - 1.47696_r8 * tk1002 + &
          salt_lim * (.025695_r8 - .025225_r8 * tk100 + 0.0049867_r8 * tk1002)
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, ff, nx_block)
#else
    ff = EXP(arg)
#endif

    !---------------------------------------------------------------------------
    !   K0 from Weiss 1974
    !---------------------------------------------------------------------------

    arg = 93.4517_r8 / tk100 - 60.2409_r8 + 23.3585_r8 * (dlogtk + LOG(1e-2_r8)) + &
          salt_lim * (.023517_r8 - 0.023656_r8 * tk100 + 0.0047036_r8 * tk1002)
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k0, nx_block)
#else
    k0 = EXP(arg)
#endif

    !---------------------------------------------------------------------------
    !   k1 = [H][HCO3]/[H2CO3]
    !   k2 = [H][CO3]/[HCO3]
    !   if k1_k2_pH_tot == .true., then use
    !      Lueker, Dickson, Keeling (2000) using Mehrbach et al. data on total scale
    !   otherwise, use
    !      Millero p.664 (1995) using Mehrbach et al. data on seawater scale
    !      this is only present to be consistent w/ OCMIP2 code
    !      it should not be used for new runs
    !      the only reason to use it is to be compatible with prior
    !      long spun up runs that had used it
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    IF (k1_k2_pH_tot) THEN
       ! total pH scale
       arg = 3633.86_r8 * invtk - 61.2172_r8 + &
             9.67770_r8 * dlogtk - 0.011555_r8 * salt_lim + &
             0.0001152_r8 * s2
    ELSE
       ! seawater pH scale, see comment above
       arg = 3670.7_r8 * invtk - 62.008_r8 + &
             9.7944_r8 * dlogtk - 0.0118_r8 * salt_lim + &
             0.000116_r8 * s2
    END IF
    arg = -LOG(c10) * arg
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k1, nx_block)
#else
    k1 = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -25.5_r8 + 0.1271_r8 * temp
       Kappa  = (-3.08_r8 + 0.0877_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, nx_block)
#else
       Kfac = EXP(lnKfac)
#endif
       k1 = k1 * Kfac
    END IF

    IF (k1_k2_pH_tot) THEN
       ! total pH scale
       arg = 471.78_r8 * invtk + 25.9290_r8 - &
             3.16967_r8 * dlogtk - 0.01781_r8 * salt_lim + 0.0001122_r8 * s2
    ELSE
       ! seawater pH scale, see comment above
       arg = 1394.7_r8 * invtk + 4.777_r8 - &
             0.0184_r8 * salt_lim + 0.000118_r8 * s2
    END IF
    arg = -LOG(c10) * arg
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k2, nx_block)
#else
    k2 = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -15.82_r8 - 0.0219_r8 * temp
       Kappa  = (1.13_r8 - 0.1475_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, nx_block)
#else
       Kfac = EXP(lnKfac)
#endif
       k2 = k2 * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   kb = [H][BO2]/[HBO2]
    !   Millero p.669 (1995) using data from Dickson (1990)
    !   CO2SYS states that this in on total pH scale
    !   pressure correction from Millero 1979, p. 1657
    !      omitting salinity contribution
    !---------------------------------------------------------------------------

    arg = (-8966.90_r8 - 2890.53_r8 * sqrts - &
           77.942_r8 * salt_lim + 1.728_r8 * salt_lim * sqrts - &
           0.0996_r8 * s2) * invtk + &
          (148.0248_r8 + 137.1942_r8 * sqrts + 1.62142_r8 * salt_lim) + &
          (-24.4344_r8 - 25.085_r8 * sqrts - 0.2474_r8 * salt_lim) * dlogtk + &
          0.053105_r8 * sqrts * tk
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, kb(:,iblock), nx_block)
#else
    kb(:,iblock) = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -29.48_r8 + (0.1622_r8 - 0.002608_r8 * temp) * temp
       Kappa  = -2.84_r8 * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, nx_block)
#else
       Kfac = EXP(lnKfac)
#endif
       kb(:,iblock) = kb(:,iblock) * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   k1p = [H][H2PO4]/[H3PO4]
    !   DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -4576.752_r8 * invtk + 115.525_r8 - &
          18.453_r8 * dlogtk + &
          (-106.736_r8 * invtk + 0.69171_r8) * sqrts + &
          (-0.65643_r8 * invtk - 0.01844_r8) * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k1p(:,iblock), nx_block)
#else
    k1p(:,iblock) = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -14.51_r8 + (0.1211_r8 - 0.000321_r8 * temp) * temp
       Kappa  = (-2.67_r8 + 0.0427_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, nx_block)
#else
       Kfac = EXP(lnKfac)
#endif
       k1p(:,iblock) = k1p(:,iblock) * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   k2p = [H][HPO4]/[H2PO4]
    !   DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -8814.715_r8 * invtk + 172.0883_r8 - &
          27.927_r8 * dlogtk + &
          (-160.340_r8 * invtk + 1.3566_r8) * sqrts + &
          (0.37335_r8 * invtk - 0.05778_r8) * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k2p(:,iblock), nx_block)
#else
    k2p(:,iblock) = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -23.12_r8 + (0.1758_r8 - 0.002647_r8 * temp) * temp
       Kappa  = (-5.15_r8 + 0.09_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, nx_block)
#else
       Kfac = EXP(lnKfac)
#endif
       k2p(:,iblock) = k2p(:,iblock) * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   k3p = [H][PO4]/[HPO4]
    !   DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -3070.75_r8 * invtk - 18.141_r8 + &
          (17.27039_r8 * invtk + 2.81197_r8) * sqrts + &
          (-44.99486_r8 * invtk - 0.09984_r8) * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, k3p(:,iblock), nx_block)
#else
    k3p(:,iblock) = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -26.57_r8 + (0.202_r8 - 0.003042_r8 * temp) * temp
       Kappa  = (-4.08_r8 + 0.0714_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, nx_block)
#else
       Kfac = EXP(lnKfac)
#endif
       k3p(:,iblock) = k3p(:,iblock) * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   ksi = [H][SiO(OH)3]/[Si(OH)4]
    !   Millero p.671 (1995) using data from Yao and Millero (1995)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !      apply boric acid values
    !---------------------------------------------------------------------------

    arg = -8904.2_r8 * invtk + 117.385_r8 - &
          19.334_r8 * dlogtk + &
          (-458.79_r8 * invtk + 3.5913_r8) * sqrtis + &
          (188.74_r8 * invtk - 1.5998_r8) * is + &
          (-12.1652_r8 * invtk + 0.07871_r8) * is2 + &
          log_1_m_1p005em3_s
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, ksi(:,iblock), nx_block)
#else
    ksi(:,iblock) = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -29.48_r8 + (0.1622_r8 - 0.002608_r8 * temp) * temp
       Kappa  = -2.84_r8 * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, nx_block)
#else
       Kfac = EXP(lnKfac)
#endif
       ksi(:,iblock) = ksi(:,iblock) * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   kw = [H][OH]
    !   Millero p.670 (1995) using composite data
    !   following DOE Handbook, 0.015 substracted from constant to
    !   approximately convert from SWS pH scale to total pH scale
    !   pressure correction from Millero 1983
    !      note that deltaV coeffs in Millero 1995 are those actually
    !      freshwater deltaV coeffs from Millero 1983
    !---------------------------------------------------------------------------

    arg = -13847.26_r8 * invtk + 148.9652_r8 - 23.6521_r8 * dlogtk + &
          (118.67_r8 * invtk - 5.977_r8 + 1.0495_r8 * dlogtk) * sqrts - &
          0.01615_r8 * salt_lim
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, kw(:,iblock), nx_block)
#else
    kw = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -20.02_r8 + (0.1119_r8 - 0.001409_r8 * temp) * temp
       Kappa  = (-5.13_r8 + 0.0794_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, nx_block)
#else
       Kfac = EXP(lnKfac)
#endif
       kw(:,iblock) = kw(:,iblock) * Kfac
    END IF

    !---------------------------------------------------------------------------
    !   ks = [H][SO4]/[HSO4], free pH scale
    !   Dickson (1990, J. chem. Thermodynamics 22, 113)
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------------

    arg = -4276.1_r8 * invtk + 141.328_r8 - 23.093_r8 * dlogtk + &
          (-13856.0_r8 * invtk + 324.57_r8 - 47.986_r8 * dlogtk) * sqrtis + &
          (35474.0_r8 * invtk - 771.54_r8 + 114.723_r8 * dlogtk) * is - &
          2698.0_r8 * invtk * is * sqrtis + &
          1776.0_r8 * invtk * is2 + &
          log_1_m_1p005em3_s
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, ks(:,iblock), nx_block)
#else
    ks(:,iblock) = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -18.03_r8 + (0.0466_r8 + 0.000316_r8 * temp) * temp
       Kappa  = (-4.53_r8 + 0.09_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, nx_block)
#else
       Kfac = EXP(lnKfac)
#endif
       ks(:,iblock) = ks(:,iblock) * Kfac
    END IF

    !---------------------------------------------------------------------
    !   kf = [H][F]/[HF]
    !   Dickson and Riley (1979) -- change pH scale to total
    !   pressure correction from Millero 1995, p. 675
    !      w/ typo corrections from CO2SYS
    !---------------------------------------------------------------------

    arg = c1 + (0.1400_r8 / 96.062_r8) * (scl) / ks(:,iblock)
#ifdef CCSMCOUPLED
       CALL shr_vmath_log(arg, log_1_p_tot_sulfate_div_ks, nx_block)
#else
    log_1_p_tot_sulfate_div_ks = LOG(arg)
#endif
    arg = 1590.2_r8 * invtk - 12.641_r8 + 1.525_r8 * sqrtis + &
          log_1_m_1p005em3_s + log_1_p_tot_sulfate_div_ks
#ifdef CCSMCOUPLED
    CALL shr_vmath_exp(arg, kf(:,iblock), nx_block)
#else
    kf(:,iblock) = EXP(arg)
#endif

    IF (k > 1) THEN
       deltaV = -9.78_r8 - (0.009_r8 + 0.000942_r8 * temp) * temp
       Kappa  = (-3.91_r8 + 0.054_r8 * temp) * p001
       lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(lnKfac, Kfac, nx_block)
#else
       Kfac = EXP(lnKfac)
#endif
       kf(:,iblock) = kf(:,iblock) * Kfac
    END IF

    !---------------------------------------------------------------------
    !   Calculate concentrations for borate, sulfate, and fluoride
    !   bt : Uppstrom (1974)
    !   st : Morris & Riley (1966)
    !   ft : Riley (1965)
    !---------------------------------------------------------------------

    bt(:,iblock) = 0.000232_r8 / 10.811_r8 * scl
    st(:,iblock) = 0.14_r8 / 96.062_r8 * scl
    ft(:,iblock) = 0.000067_r8 / 18.9984_r8 * scl

  END SUBROUTINE comp_co3_coeffs

  !*****************************************************************************

  SUBROUTINE comp_htotal(iblock, j, k, mask, temp, dic_in, ta_in, pt_in, sit_in, &
                         k1, k2, phlo, phhi, htotal)

    !---------------------------------------------------------------------------
    !   SUBROUTINE comp_htotal
    !
    !   PURPOSE : Calculate htotal from total alkalinity, total CO2,
    !             temp, salinity (s), etc.
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: iblock, j
    INTEGER(KIND=int_kind), INTENT(IN) :: k
    LOGICAL(KIND=log_kind), DIMENSION(nx_block), INTENT(IN) :: mask
    REAL(KIND=r8), DIMENSION(nx_block), INTENT(IN) :: &
         temp,     & ! temperature (degrees C)
         dic_in,   & ! total inorganic carbon (nmol/cm^3)
         ta_in,    & ! total alkalinity (neq/cm^3)
         pt_in,    & ! inorganic phosphate (nmol/cm^3)
         sit_in,   & ! inorganic silicate (nmol/cm^3)
         k1,k2       ! equilibrium constants for CO2 species

    !---------------------------------------------------------------------------
    !   input/output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=r8), DIMENSION(nx_block), INTENT(INOUT) :: &
         phlo,     & ! lower limit of pH range
         phhi        ! upper limit of pH range

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=r8), DIMENSION(nx_block), INTENT(OUT) :: &
         htotal      ! free concentration of H ion

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind) :: i

    REAL(KIND=r8) :: &
         mass_to_vol,  & ! (mol/kg) -> (mmol/m^3)
         vol_to_mass     ! (mmol/m^3) -> (mol/kg)

    REAL(KIND=r8), DIMENSION(nx_block) :: &
         x1, x2          ! bounds on htotal for solver

    !---------------------------------------------------------------------------
    !   check for existence of ocean points
    !---------------------------------------------------------------------------

    IF (COUNT(mask) == 0) THEN
       htotal = c0
       RETURN
    END IF

    !---------------------------------------------------------------------------
    !   set unit conversion factors
    !---------------------------------------------------------------------------

    mass_to_vol = 1e6_r8 * rho_sw
    vol_to_mass = c1 / mass_to_vol

    !---------------------------------------------------------------------------
    !   convert tracer units to per mass
    !---------------------------------------------------------------------------

    DO i = 1,nx_block
       IF (mask(i)) THEN
          dic(i,iblock)  = max(dic_in(i),dic_min) * vol_to_mass
          ta(i,iblock)   = max(ta_in(i),alk_min)  * vol_to_mass
          pt(i,iblock)   = max(pt_in(i),c0)       * vol_to_mass
          sit(i,iblock)  = max(sit_in(i),c0)      * vol_to_mass

          x1(i) = c10 ** (-phhi(i))
          x2(i) = c10 ** (-phlo(i))
       END IF ! if mask
    END DO ! i loop

    !---------------------------------------------------------------------------
    !   If DIC and TA are known then either a root finding or iterative
    !   method must be used to calculate htotal. In this case we use
    !   the Newton-Raphson "safe" method taken from "Numerical Recipes"
    !   (function "rtsafe.f" with error trapping removed).
    !
    !   As currently set, this procedure iterates about 12 times. The
    !   x1 and x2 values set below will accomodate ANY oceanographic
    !   values. If an initial guess of the pH is known, then the
    !   number of iterations can be reduced to about 5 by narrowing
    !   the gap between x1 and x2. It is recommended that the first
    !   few time steps be run with x1 and x2 set as below. After that,
    !   set x1 and x2 to the previous value of the pH +/- ~0.5.
    !---------------------------------------------------------------------------

    CALL drtsafe_row(iblock, j, k, mask, k1, k2, x1, x2, xacc, htotal)

  END SUBROUTINE comp_htotal

  !*****************************************************************************

  SUBROUTINE drtsafe_row(iblock, j, k, mask_in, k1, k2, x1, x2, xacc, soln)

    !---------------------------------------------------------------------------
    !   Vectorized version of drtsafe, which was a modified version of
    !      Numerical Recipes algorithm.
    !   Keith Lindsay, Oct 1999
    !
    !   Algorithm comment :
    !      Iteration from Newtons method is used unless it leaves
    !      bracketing interval or the dx is > 0.5 the previous dx.
    !      In that case, bisection method is used.
    !---------------------------------------------------------------------------

#ifdef CCSMCOUPLED
    USE shr_sys_mod, ONLY : shr_sys_abort
#endif

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: iblock, j, k
    LOGICAL(KIND=log_kind), DIMENSION(nx_block), INTENT(IN) :: mask_in
    REAL(KIND=r8), DIMENSION(nx_block), INTENT(IN) :: k1, k2
    REAL(KIND=r8), INTENT(IN) :: xacc

    !---------------------------------------------------------------------------
    !   input/output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=r8), DIMENSION(nx_block), INTENT(INOUT) :: x1, x2

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=r8), DIMENSION(nx_block), INTENT(OUT) :: soln

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    LOGICAL(KIND=log_kind) :: leave_bracket, dx_decrease
    LOGICAL(KIND=log_kind), DIMENSION(nx_block) :: mask
    INTEGER(KIND=int_kind) ::  i, it
    REAL(KIND=r8) :: temp
    REAL(KIND=r8), DIMENSION(nx_block) :: xlo, xhi, flo, fhi, f, df, dxold, dx
    TYPE (block) :: this_block

    !---------------------------------------------------------------------------
    !   bracket root at each location and set up first iteration
    !---------------------------------------------------------------------------

    mask = mask_in

    it = 0

    DO
       CALL talk_row(iblock, mask, k1, k2, x1, flo, df)
       CALL talk_row(iblock, mask, k1, k2, x2, fhi, df)

       WHERE ( mask )
          mask = (flo > c0 .AND. fhi > c0) .OR. &
                 (flo < c0 .AND. fhi < c0)
       END WHERE

       IF (.NOT. ANY(mask)) EXIT

       it = it + 1

       DO i = 1,nx_block
          IF (mask(i)) THEN
             this_block = get_block(blocks_clinic(iblock), iblock)
             WRITE(stdout,*) '(co2calc.F90:drtsafe_row) ', &
                'i_glob = ', this_block%i_glob(i), &
                ', j_glob = ', this_block%j_glob(j), ', k = ', k, &
                ', nsteps_run = ', nsteps_run, ', it = ', it
             WRITE(stdout,*) '(co2calc.F90:drtsafe_row) ', &
                '   x1,f = ', x1(i), flo(i)
             WRITE(stdout,*) '(co2calc.F90:drtsafe_row) ', &
                '   x2,f = ', x2(i), fhi(i)
          END IF
       END DO

       IF (it > max_bracket_grow_it) THEN
          CALL shr_sys_abort('bounding bracket for pH solution not found')
       END IF

       WHERE ( mask )
          dx = sqrt(x2 / x1)
          x2 = x2 * dx
          x1 = x1 / dx
       END WHERE
    END DO

    mask = mask_in

    DO i = 1,nx_block
       IF (mask(i)) THEN
          IF (flo(i) .LT. c0) THEN
             xlo(i) = x1(i)
             xhi(i) = x2(i)
          ELSE
             xlo(i) = x2(i)
             xhi(i) = x1(i)
             temp = flo(i)
             flo(i) = fhi(i)
             fhi(i) = temp
          END IF
          soln(i) = p5 * (xlo(i) + xhi(i))
          dxold(i) = ABS(xlo(i) - xhi(i))
          dx(i) = dxold(i)
       END IF
    END DO

    CALL talk_row(iblock, mask, k1, k2, soln, f, df)

    !---------------------------------------------------------------------------
    !   perform iterations, zeroing mask when a location has converged
    !---------------------------------------------------------------------------

    DO it = 1,maxit
       DO i = 1,nx_block
          IF (mask(i)) THEN
             leave_bracket = ((soln(i) - xhi(i)) * df(i) - f(i)) * &
                  ((soln(i) - xlo(i)) * df(i) - f(i)) .GE. 0
             dx_decrease = ABS(c2 * f(i)) .LE. ABS(dxold(i) * df(i))
             IF (leave_bracket .OR. .NOT. dx_decrease) THEN
                dxold(i) = dx(i)
                dx(i) = p5 * (xhi(i) - xlo(i))
                soln(i) = xlo(i) + dx(i)
                IF (xlo(i) .EQ. soln(i)) mask(i) = .FALSE.
             ELSE
                dxold(i) = dx(i)
                dx(i) = -f(i) / df(i)
                temp = soln(i)
                soln(i) = soln(i) + dx(i)
                IF (temp .EQ. soln(i)) mask(i) = .FALSE.
             END IF
             IF (ABS(dx(i)) .LT. xacc) mask(i) = .FALSE.
          END IF
       END DO

       IF (.NOT. ANY(mask)) RETURN

       CALL talk_row(iblock, mask, k1, k2, soln, f, df)

       DO i = 1,nx_block
          IF (mask(i)) THEN
             IF (f(i) .LT. c0) THEN
                xlo(i) = soln(i)
                flo(i) = f(i)
             ELSE
                xhi(i) = soln(i)
                fhi(i) = f(i)
             END IF
          END IF
       END DO

    END DO ! iteration loop

#ifdef CCSMCOUPLED
    CALL shr_sys_abort('lack of convergence in drtsafe_row')
#endif

  END SUBROUTINE drtsafe_row

  !*****************************************************************************

  SUBROUTINE talk_row(iblock, mask, k1, k2, x, fn, df)

    !---------------------------------------------------------------------------
    !   This routine computes TA as a function of DIC, htotal and constants.
    !   It also calculates the derivative of this function with respect to
    !   htotal. It is used in the iterative solution for htotal. In the call
    !   "x" is the input value for htotal, "fn" is the calculated value for
    !   TA and "df" is the value for dTA/dhtotal.
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: iblock
    LOGICAL(KIND=log_kind), DIMENSION(nx_block), INTENT(IN) :: mask
    REAL(KIND=r8), DIMENSION(nx_block), INTENT(IN) :: k1, k2, x

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=r8), DIMENSION(nx_block), INTENT(OUT) :: fn, df

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind) :: i

    REAL(KIND=r8) :: &
         x1, x1_r, x2, x2_r, x3, k12, k12p, k123p, &
         a, a_r, a2_r, da, b, b_r, b2_r, db, c, c_r, &
         kb_p_x1_r, ksi_p_x1_r, c1_p_c_ks_x1_r_r, c1_p_kf_x1_r_r

    !---------------------------------------------------------------------------

    DO i = 1,nx_block
       IF (mask(i)) THEN
          x1 = x(i)
          x1_r = c1 / x1
          x2 = x1 * x1
          x2_r = x1_r * x1_r
          x3 = x2 * x1
          k12 = k1(i) * k2(i)
          k12p = k1p(i,iblock) * k2p(i,iblock)
          k123p = k12p * k3p(i,iblock)
          a = x3 + k1p(i,iblock) * x2 + k12p * x1 + k123p
          a_r = c1 / a
          a2_r = a_r * a_r
          da = c3 * x2 + c2 * k1p(i,iblock) * x1 + k12p
          b = x2 + k1(i) * x1 + k12
          b_r = c1 / b
          b2_r = b_r * b_r
          db = c2 * x1 + k1(i)
          c = c1 + st(i,iblock) / ks(i,iblock)
          c_r = c1 / c
          kb_p_x1_r = c1 / (kb(i,iblock) + x1)
          ksi_p_x1_r = c1 / (ksi(i,iblock) + x1)
          c1_p_c_ks_x1_r_r = c1 / (c1 + c * ks(i,iblock) * x1_r)
          c1_p_kf_x1_r_r = c1 / (c1 + kf(i,iblock) * x1_r)

          !---------------------------------------------------------------------
          !   fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
          !---------------------------------------------------------------------

          fn(i) = k1(i) * dic(i,iblock) * x1 * b_r &
               + c2 * dic(i,iblock) * k12 * b_r &
               + bt(i,iblock) * kb(i,iblock) * kb_p_x1_r &
               + kw(i,iblock) * x1_r &
               + pt(i,iblock) * k12p * x1 * a_r &
               + c2 * pt(i,iblock) * k123p * a_r &
               + sit(i,iblock) * ksi(i,iblock) * ksi_p_x1_r &
               - x1 * c_r &
               - st(i,iblock) * c1_p_c_ks_x1_r_r &
               - ft(i,iblock) * c1_p_kf_x1_r_r &
               - pt(i,iblock) * x3 * a_r &
               - ta(i,iblock)

          !---------------------------------------------------------------------
          !   df = d(fn)/dx
          !---------------------------------------------------------------------

          df(i) = k1(i) * dic(i,iblock) * (b - x1 * db) * b2_r &
               - c2 * dic(i,iblock) * k12 * db * b2_r &
               - bt(i,iblock) * kb(i,iblock) * kb_p_x1_r * kb_p_x1_r &
               - kw(i,iblock) * x2_r &
               + (pt(i,iblock) * k12p * (a - x1 * da)) * a2_r &
               - c2 * pt(i,iblock) * k123p * da * a2_r &
               - sit(i,iblock) * ksi(i,iblock) * ksi_p_x1_r * ksi_p_x1_r &
               - c1 * c_r &
               - st(i,iblock) * c1_p_c_ks_x1_r_r * c1_p_c_ks_x1_r_r * (c * ks(i,iblock) * x2_r) &
               - ft(i,iblock) * c1_p_kf_x1_r_r * c1_p_kf_x1_r_r * kf(i,iblock) * x2_r &
               - pt(i,iblock) * x2 * (c3 * a - x1 * da) * a2_r

       END IF ! if mask
    END DO ! i loop

  END SUBROUTINE talk_row

  !*****************************************************************************

  SUBROUTINE comp_co3_sat_vals(k, mask, temp, salt, co3_sat_calc, co3_sat_arag)

    !---------------------------------------------------------------------------
    !   SUBROUTINE comp_co3_sat_vals
    !
    !   PURPOSE : Calculate co3 concentration at calcite and aragonite saturation
    !             from temp, salinity (s), press
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   input arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: k
    LOGICAL(KIND=log_kind), DIMENSION(nx_block, ny_block), INTENT(IN) :: mask
    REAL(KIND=r8), DIMENSION(nx_block,ny_block), INTENT(IN) :: &
         temp,     & ! temperature (degrees C)
         salt        ! salinity (PSU)

    !---------------------------------------------------------------------------
    !   output arguments
    !---------------------------------------------------------------------------

    REAL(KIND=r8), DIMENSION(nx_block, ny_block), INTENT(OUT) :: &
         co3_sat_calc,&! co3 concentration at calcite saturation
         co3_sat_arag  ! co3 concentration at aragonite saturation

    !---------------------------------------------------------------------------
    !   local variable declarations
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind) :: i, j

    REAL(KIND=r8) :: &
         mass_to_vol,  & ! (mol/kg) -> (mmol/m^3)
         press_bar       ! pressure at level k [bars]

    REAL(KIND=r8), DIMENSION(nx_block) :: &
         salt_lim,     & ! bounded salt
         tk,           & ! temperature (K)
         log10tk,invtk,sqrts,s15,invRtk,arg,&
         K_calc,       & ! thermodynamic constant for calcite
         K_arag,       & ! thermodynamic constant for aragonite
         deltaV,Kappa, & ! pressure correction terms
         lnKfac,Kfac,  & ! pressure correction terms
         inv_Ca          ! inverse of Calcium concentration (mol/kg)

    !---------------------------------------------------------------------------
    !   check for existence of ocean points
    !---------------------------------------------------------------------------

    IF (COUNT(mask) == 0) THEN
       co3_sat_calc = c0
       co3_sat_arag = c0
       RETURN
    END IF

    !---------------------------------------------------------------------------
    !   set unit conversion factors
    !---------------------------------------------------------------------------

    mass_to_vol = 1e6_r8 * rho_sw

    !---------------------------------------------------------------------------

    press_bar = ref_pressure(k)

    DO j = 1,ny_block

       !------------------------------------------------------------------------
       !   check for existence of ocean points on this row
       !------------------------------------------------------------------------

       IF (COUNT(mask(:,j)) == 0) THEN
          co3_sat_calc(:,j) = c0
          co3_sat_arag(:,j) = c0
          CYCLE
       END IF

       salt_lim = max(salt(:,j),salt_min)
       tk       = T0_Kelvin + temp(:,j)
#ifdef CCSMCOUPLED
       CALL shr_vmath_log(tk, log10tk, nx_block)
#else
       log10tk  = LOG(tk)
#endif
       log10tk  = log10tk/LOG(c10)
       invtk    = c1 / tk
       invRtk   = (c1 / 83.1451_r8) * invtk

#ifdef CCSMCOUPLED
       CALL shr_vmath_sqrt(salt_lim, sqrts, nx_block)
#else
       sqrts    = SQRT(salt_lim)
#endif
       s15      = sqrts * salt_lim

       !------------------------------------------------------------------------
       !   Compute K_calc, K_arag, apply pressure factor
       !   Mucci, Amer. J. of Science 283:781-799, 1983 & Millero 1979
       !------------------------------------------------------------------------

       arg = -171.9065_r8 - 0.077993_r8 * tk + 2839.319_r8 * invtk + 71.595_r8 * log10tk + &
             (-0.77712_r8 + 0.0028426_r8 * tk + 178.34_r8 * invtk) * sqrts - &
             0.07711_r8 * salt_lim + 0.0041249_r8 * s15
       arg = LOG(c10) * arg
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(arg, K_calc, nx_block)
#else
       K_calc = EXP(arg)
#endif

       IF (k > 1) THEN
          deltaV = -48.76_r8 + 0.5304_r8 * temp(:,j)
          Kappa  = (-11.76_r8 + 0.3692_r8 * temp(:,j)) * p001
          lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
          CALL shr_vmath_exp(lnKfac, Kfac, nx_block)
#else
          Kfac = EXP(lnKfac)
#endif
          K_calc = K_calc * Kfac
       END IF

       arg = -171.945_r8 - 0.077993_r8 * tk + 2903.293_r8 * invtk + 71.595_r8 * log10tk + &
            (-0.068393_r8 + 0.0017276_r8 * tk + 88.135_r8 * invtk) * sqrts - &
            0.10018_r8 * salt_lim + 0.0059415_r8 * s15
       arg = LOG(c10) * arg
#ifdef CCSMCOUPLED
       CALL shr_vmath_exp(arg, K_arag, nx_block)
#else
       K_arag = EXP(arg)
#endif

       IF (k > 1) THEN
          deltaV = deltaV + 2.8_r8
          lnKfac = (-deltaV + p5 * Kappa * press_bar) * press_bar * invRtk
#ifdef CCSMCOUPLED
          CALL shr_vmath_exp(lnKfac, Kfac, nx_block)
#else
          Kfac = EXP(lnKfac)
#endif
          K_arag = K_arag * Kfac
       END IF

       WHERE ( mask(:,j) )

          !------------------------------------------------------------------
          !   Compute CO3 concentration at calcite & aragonite saturation
          !------------------------------------------------------------------

          inv_Ca = (35.0_r8 / 0.01028_r8) / salt_lim 
          co3_sat_calc(:,j) = K_calc * inv_Ca
          co3_sat_arag(:,j) = K_arag * inv_Ca

          !------------------------------------------------------------------
          !   Convert units of output arguments
          !------------------------------------------------------------------

          co3_sat_calc(:,j) = co3_sat_calc(:,j) * mass_to_vol
          co3_sat_arag(:,j) = co3_sat_arag(:,j) * mass_to_vol

       ELSEWHERE

          co3_sat_calc(:,j) = c0
          co3_sat_arag(:,j) = c0

       END WHERE

    END DO ! j loop

  END SUBROUTINE comp_co3_sat_vals

  !*****************************************************************************

END MODULE co2calc
