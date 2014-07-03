! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  Calculates particle production rates due to nucleation <rhompe>:
!!  binary homogeneous nucleation of sulfuric acid and water only
!!  Numerical method follows Zhao & Turco, JAS, V.26, No.5, 1995.
!!
!!  This was moved from sulfnuc to make the code more manageable.
!!
!!  @author Mike Mills, Chuck Bardeen
!!  @version Sep-2011
subroutine sulfnucrate(carma,cstate, iz, igroup, nucbin, nucrate, rc) 
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  use sulfate_utils
  
  implicit none
  
  type(carma_type), intent(in)         :: carma       !! the carma object
  type(carmastate_type), intent(inout) :: cstate      !! the carma state object
  integer, intent(in)                  :: iz          !! level index
  integer, intent(in)                  :: igroup      !! group index
  integer, intent(out)                 :: nucbin      !! bin in which nucleation occurs
  real(kind=f), intent(out)            :: nucrate     !! nucleation rate #/x/y/z/s
  integer, intent(inout)               :: rc          !! return code, negative indicates failure

  !  Local declarations     
  integer           :: i, ibin, ie    
  real(kind=f)      :: dens(46)
  real(kind=f)      :: pa(46)
  real(kind=f)      :: pb(46)
  real(kind=f)      :: c1(46)
  real(kind=f)      :: c2(46)
  real(kind=f)      :: fct(46)
  real(kind=f)      :: wtmolr         ! molecular weight ration of H2SO4/H2O
  real(kind=f)      :: h2o_cgs        ! H2O densities in g/cm3
  real(kind=f)      :: h2so4_cgs      ! H2SO4 densities in g/cm3
  real(kind=f)      :: h2o            ! H2O concentrations in molec/cm3 
  real(kind=f)      :: h2so4          ! H2SO4 concentrations in molec/cm3
  real(kind=f)      :: h2oln          ! H2O ambient vapor pressures [dynes/cm2]
  real(kind=f)      :: h2so4ln        ! H2SO4 ambient vapor pressures [dynes/cm2]
  real(kind=f)      :: rh             ! relative humidity of water wrt liquid water
  real(kind=f)      :: SA             ! total surface area of pre-existing wet particles
  real(kind=f)      :: SAbin          ! bin surface area of pre-existing wet particles
  real(kind=f)      :: cw
  real(kind=f)      :: dw
  real(kind=f)      :: wvp            ! water eq.vp over solution
  real(kind=f)      :: wvpln
  real(kind=f)      :: t0_kulm
  real(kind=f)      :: seqln
  real(kind=f)      :: t_crit_kulm
  real(kind=f)      :: factor_kulm
  real(kind=f)      :: dw1, dw2
  real(kind=f)      :: dens1
  real(kind=f)      :: dens11
  real(kind=f)      :: dens12
  real(kind=f)      :: xfrac
  real(kind=f)      :: wstar
  real(kind=f)      :: dstar
  real(kind=f)      :: rhln
  real(kind=f)      :: raln
  real(kind=f)      :: wfstar
  real(kind=f)      :: sigma
  real(kind=f)      :: ystar
  real(kind=f)      :: rstar
  real(kind=f)      :: r2
  real(kind=f)      :: gstar
  real(kind=f)      :: rb
  real(kind=f)      :: beta1
  real(kind=f)      :: beta2
  real(kind=f)      :: rpr
  real(kind=f)      :: rpre
  real(kind=f)      :: fracmol
  real(kind=f)      :: zphi
  real(kind=f)      :: zeld
  real(kind=f)      :: cfac
  real(kind=f)      :: ahom
  real(kind=f)      :: exhom
  real(kind=f)      :: rmstar
  real(kind=f)      :: frac_h2so4
  real(kind=f)      :: rhomlim  
  real(kind=f)      :: dnpot(46), dnwf(46)
  real(kind=f)      :: rho_H2SO4_wet

 5 format(/,'microfast::WARNING - nucleation rate exceeds 5.e1: ie=', i2,', iz=',i4,',lat=', &
              f7.2,',lon=',f7.2, ', rhompe=', e10.3)	      
  
  !  Parameterized fit developed by Mike Mills in 1994 to the partial molal
  !  Gibbs energies (F2|o-F2) vs. weight percent H2SO4 table in Giauque et al.,
  !  J. Am. Chem. Soc, 82, 62-70, 1960.  The parameterization gives excellent
  !  agreement.  Ayers (GRL, 7, 433-436, 1980) refers to F2|o-F2 as mu - mu_0 
  !  (chemical potential).  This parameterization may be replaced by a lookup
  !  table, as was done ultimately in the Garcia-Solomon sulfate code.
  do i = 1, 46
    dnpot(i) = 4.184_f * (23624.8_f - 1.14208e8_f / ((dnwtp(i) - 105.318_f)**2 + 4798.69_f))
    dnwf(i) = dnwtp(i) / 100._f
  end do

  ! Molecular weight ratio of H2SO4 / H2O:
  wtmolr = gwtmol(igash2so4) / gwtmol(igash2o)

  ! Compute H2O and H2SO4 densities in g/cm3      
  h2o_cgs   = gc(iz, igash2o)   / (zmet(iz) * xmet(iz) * ymet(iz))
  h2so4_cgs = gc(iz, igash2so4) / (zmet(iz) * xmet(iz) * ymet(iz))

  ! Compute H2O and H2SO4 concentrations in molec/cm3      
  h2o   = h2o_cgs   * AVG / gwtmol(igash2o)
  h2so4 = h2so4_cgs * AVG / gwtmol(igash2so4)

  ! Compute relative humidity of water wrt liquid water       
  rh = (supsatl(iz, igash2o) + 1._f) * 100._f

  ! Compute ln of H2O and H2SO4 ambient vapor pressures [dynes/cm2]
  h2oln   = log(h2o_cgs   * (RGAS / gwtmol(igash2o))   * t(iz))
  h2so4ln = log(h2so4_cgs * (RGAS / gwtmol(igash2so4)) * t(iz))

  ! loop through wt pcts and calculate vp/composition for each
  do i = 1, 46
    dens(i) = dnc0(i) + dnc1(i) * t(iz)

    ! Calc. water eq.vp over solution using (Lin & Tabazadeh eqn 5, JGR, 2001)  
    cw = 22.7490_f + 0.0424817_f * dnwtp(i) - 0.0567432_f * dnwtp(i)**0.5_f - 0.000621533_f * dnwtp(i)**2
    dw = -5850.24_f + 21.9744_f * dnwtp(i) - 44.5210_f * dnwtp(i)**0.5_f - 0.384362_f * dnwtp(i)**2
    
    ! pH20 | eq[mb]
    wvp   = exp(cw + dw / t(iz))
    
    ! Ln(pH2O | eq [dynes/cm2])
    wvpln = log(wvp * 1013250._f / 1013.25_f)
        
    ! Save the water eq.vp over solution at each wt pct into this array:
    !
    ! Ln(pH2O/pH2O|eq) with both terms in dynes/cm2
    pb(i) = h2oln - wvpln

    ! Calc. sulfuric acid eq.vp over solution using (Ayers et. al., GRL, V.7, No.6, June 1980)
    !
    ! T0 set in the low end of the Ayers measurement range (338-445K)
    t0_kulm = 340._f
    seqln   = -10156._f / t0_kulm + 16.259_f

    ! Now calc. Kulmala correction (J. CHEM. PHYS. V.93, No.1, 1 July 1990)
    !
    ! Critical temperature = 1.5 * Boiling point
    t_crit_kulm = 905._f
    factor_kulm = -1._f / t(iz) + 1._f / t0_kulm + 0.38_f / (t_crit_kulm - t0_kulm) * &
      (1.0_f + log(t0_kulm / t(iz)) - t0_kulm / t(iz))

    ! For pure sulfuric acid
    seqln = seqln + 10156._f * factor_kulm

    ! Now adjust vp based on weight % composition using parameterization of Giauque 1960
    !
    ! Adjust for WTPCT composition
    seqln = seqln - dnpot(i) / (8.3143_f * t(iz))

    ! Convert atmospheres => dynes/cm2
    seqln = seqln + log(1013250._f)
  
    ! Save the sulfuric acid eq.vp over solution at each wt pct into this array:
    !
    ! Ln(pH2SO4/pH2SO4|eq) with both terms in dynes/cm2
    pa(i) = h2so4ln - seqln

    ! Create 2-component solutions of varying composition c1 and c2
    c1(i) = pa(i) - pb(i) * wtmolr
    c2(i) = pa(i) * dnwf(i) + pb(i) * (1._f - dnwf(i)) * wtmolr
  end do  ! end of loop through wtpcts

  ! Now loop through until we find the c1+c2 combination with minimum Gibbs free energy
  dw2     = dnwtp(46) - dnwtp(45)
  dens1   = (dens(46) - dens(45)) / dw2
  fct(46) = c1(46) + c2(46) * 100._f * dens1 / dens(46)
  dens12 = dens1
    
  do i = 45, 2, -1
    dw1    = dw2
    dens11 = dens12
    dw2    = dnwtp(i) - dnwtp(i-1)
    dens12 = (dens(i) - dens(i-1)) / dw2
    dens1  = (dens11 * dw2 + dens12 * dw1) / (dw1 + dw2)

    fct(i) = c1(i) + c2(i) * 100._f * dens1 / dens(i)

    ! Find saddle where fct(i)<0<fct(i+1)
    if (fct(i) * fct(i+1) <= 0._f) exit
  end do
  
  if (i == 1) then
    dens1  = (dens(2) - dens(1)) / (dnwtp(2) - dnwtp(1))
    fct(1) = c1(1) + c2(1) * 100._f * dens1 / dens(1)
  end if
      
  ! Possibility 1: loop finds no saddle, so no nucleation occurs:
  if (fct(i) * fct(i+1) > 0._f) then
    nucbin  = 0 
    nucrate = 0.0_f
    
    return

  ! Possibility 2: loop crossed the saddle; interpolate to find exact value:
  else if (fct(i) * fct(i+1) < 0._f) then
    xfrac = fct(i+1) / (fct(i+1) - fct(i))
    wstar = dnwtp(i+1) * (1.0_f - xfrac) + dnwtp(i) * xfrac ! critical wtpct
    dstar = dens(i+1)  * (1.0_f - xfrac) + dens(i)  * xfrac
    rhln  = pb(i+1) * (1.0_f - xfrac) + pb(i) * xfrac
    raln  = pa(i+1) * (1.0_f - xfrac) + pa(i) * xfrac

  ! Possibility 3: loop found the saddle point exactly
  else
    dstar = dens(i)
  
    ! critical wtpct
    wstar = dnwtp(i)
    rhln  = pb(i)
    raln  = pa(i)
  end if

  ! Critical weight fraction
  wfstar = wstar / 100._f

  if ((wfstar < 0._f) .or. (wfstar > 1._f)) then
    write(LUNOPRT,*)'sulfnuc: wstar out of bounds!'
    rc = RC_ERROR
    return
  end if
      
  ! Critical surface tension  [erg/cm2]
  sigma = sulfate_surf_tens(carma, wstar, t(iz), rc)  
      
  ! Critical Y (eqn 13 in Zhao & Turco 1993) [erg/cm3]
  ystar = dstar * RGAS * t(iz) * (wfstar / gwtmol(igash2so4) &
      * raln + (1._f - wfstar) / gwtmol(igash2o) * rhln)
  if (ystar < 1.e-20_f) then
    nucbin  = 0 
    nucrate = 0.0_f

    return
  end if

  ! Critical cluster radius [cm]        
  rstar = 2._f * sigma / ystar 
  rstar = max(rstar, 0.0_f)
  r2    = rstar * rstar

  ! Critical Gibbs free energy [erg]
  gstar = (4._f * PI / 3._f) * r2 * sigma 
 
  !   kT/(2*Pi*M) = [erg/mol/K]*[K]/[g/mol] = [erg/g] = [cm2/s2]
  !   RB[erg/mol] = RGAS[erg/mol/K] * T[K] / (2Pi)
  rb = RGAS * t(iz) / 2._f / PI
      
  ! Beta[cm/s] = sqrt(RB[erg/mol] / WTMOL[g/mol])
  beta1 = sqrt(rb / gwtmol(igash2so4)) ! H2SO4
  beta2 = sqrt(rb / gwtmol(igash2o))   ! H2O

  ! RPR[molecules/s] = 4Pi * R2[cm2] * H2O[molecules/cm3] * Beta[cm/s]
  rpr = 4._f * PI * r2 * h2o * beta1

  ! RPRE[/cm3/s] = RPR[/s] * H2SO4[/cm3]; first part of Zhao & Turco eqn 16
  rpre = rpr * h2so4

  ! Zeldovitch non-equilibrium correction factor [unitless]
  ! Jaecker-Voirol & Mirabel, 1988 (not considered in Zhao & Turco) 
  fracmol = 1._f /(1._f + wtmolr * (1._f - wfstar) / wfstar)
  zphi    = atan(fracmol)
  zeld    = 0.25_f / (sin(zphi))**2

  ! Empirical correction factor:
  cfac = 0.0_f

  ! Gstar exponential term in Zhao & Turco eqn 16 [unitless]
  ahom = (-gstar / BK / t(iz)) + cfac
  if (ahom .lt. -500._f) then
    exhom=0.0_f
  else
    exhom = exp(min(ahom, 28.0_f))
  endif

  !   Calculate mass of critical nucleus
  rho_H2SO4_wet = sulfate_density(carma, wtpct(iz),t(iz), rc)
  rmstar = (4._f * PI / 3._f) * rho_H2SO4_wet * r2 * rstar
     
  ! Calculate dry mass of critical nucleus
  rmstar = rmstar * wfstar

  !   Calc bin # of crit nucleus
  if (rmstar.lt.rmassup(1,igroup)) then
    nucbin = 1
  else
    nucbin = 2 + int(log(rmstar / rmassup(1,igroup)) / log(rmrat(igroup)))
  endif
 
  ! If none of the bins are large enough for the critical radius, then
  ! no nucleation will occur.
  if (nucbin > NBIN) then
    nucbin  = 0 
    nucrate = 0.0_f
  else
    ! Calculate the nucleation rate [#/cm3/s], Zhao & Turco eqn 16.
    nucrate = rpre * zeld * exhom
        
    ! Scale to #/x/y/z/s
    nucrate = nucrate * zmet(iz) * xmet(iz) * ymet(iz)
  endif
    
  return
end subroutine sulfnucrate
     
