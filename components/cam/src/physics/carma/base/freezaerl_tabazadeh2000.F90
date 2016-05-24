! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluates particle loss rates due to nucleation <rnuclg>:
!! aerosol freezing only.
!!
!! The parameterization described by Tabazadeh et al. [GRL, 27, 1111, 2000.] is
!! used.
!! 
!! The loss rates for all particle elements in a particle group are equal.
!!
!! @author Eric Jensen, Chuck Bardeen
!! @version Mar-1995, Nov-2009
subroutine freezaerl_tabazadeh2000(carma, cstate, iz, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(in)                  :: iz      !! z index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  !  Local declarations
  ! Define parameters needed for freezing nucleation calculations.
!  real(kind=f), parameter ::  adelf = 1.29e-12_f
!  real(kind=f), parameter ::  bdelf = 0.05_f
  real(kind=f), parameter ::  prenuc = 2.075e33_f * RHO_W / RHO_I
!  real(kind=f), parameter ::  rmiv   = 0.6_f
  
  integer                              :: igas    !! gas index
  integer                              :: igroup  !! group index
  integer                              :: ibin    !! bin index
  integer                              :: iepart  !! element for condensing group index
  integer                              :: inuc    !! nucleating element index
  integer                              :: ienucto !! index of target nucleation element
  integer                              :: ignucto !! index of target nucleation group
  integer                              :: isol
  real(kind=f)                         :: A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10
  real(kind=f)                         :: c0, C1, C2, C3, C4, c5
  real(kind=f)                         :: d0, d1, d2, d3, d4, d5
  real(kind=f)                         :: e0, e1, e2, e3, e4, e5
  real(kind=f)                         :: sifreeze
  real(kind=f)                         :: rhoibar
  real(kind=f)                         :: rlhbar
  real(kind=f)                         :: act
  real(kind=f)                         :: CONTL
  real(kind=f)                         :: CONTH
  real(kind=f)                         :: H2SO4m
  real(kind=f)                         :: WT
  real(kind=f)                         :: vrat
  real(kind=f)                         :: wtfrac
  real(kind=f)                         :: den
  real(kind=f)                         :: diffact
  real(kind=f)                         :: S260, S220, S180
  real(kind=f)                         :: sigma
  real(kind=f)                         :: sigsula
  real(kind=f)                         :: sigicea
  real(kind=f)                         :: sigsulice
  real(kind=f)                         :: ag
  real(kind=f)                         :: delfg
  real(kind=f)                         :: expon
  real(kind=f)                         :: ssl
  real(kind=f)                         :: fkelv


  ! Loop over particle groups.
  do igroup = 1,NGROUP

    igas   = inucgas(igroup)
    iepart = ienconc(igroup)
    isol   = isolelem(iepart)

    if( igas .ne. 0 )then

      ! Calculate nucleation loss rates.  Do not allow nucleation into
      ! an evaporating bin.
!     if( nnuc2elem(iepart) .gt. 1 )then
        do inuc = 1,nnuc2elem(iepart)
  
          ienucto = inuc2elem(inuc,iepart)
          if( ienucto .ne. 0 )then
            ignucto = igelem( ienucto )
    
            ! Only compute nucleation rate for aerosol freezing.
            !
            ! NOTE: If heterogeneous nucleation of glassy aerosols is being used
            ! as a nucleation mechanism, then both the heterogeneous nucleation and
            ! the homogeneous freezing need to be considered.
              if ((iand(inucproc(iepart,ienucto), I_AF_TABAZADEH_2000) /= 0)) then
      
              !  Loop over particle bins.  Loop from largest to smallest for 
              !  evaluation of index of smallest bin nucleated during time step <inucstep>.
              do ibin =NBIN,1,-1
    
                ! Bypass calculation if few particles are present
                if( pconmax(iz,igroup) .gt. FEW_PC )then
      
                  ! Calculate approximate critical saturation needed for homogeneous freezing
                  ! of sulfate aerosols (see Jensen and Toon, GRL, 1994).
                  sifreeze = 0.3_f
                  
                  ! NOTE: The wieght percent can become negative from this parameterization,
                  ! which is not physicsal. With small supersaturations, the water activity 
                  ! becomes postive (>1.013) the weight percent becomes negative. Don't allow
                  ! the the supsatl to be greater than 0.
                  ssl = max(-1.0_f, min(0._f, supsatl(iz,igas)))
                  
    
                  ! Homogeneous freezing of sulfate aerosols should only occur of SL < Scrit
                  ! and SI > <sifreeze>.
                  if( supsati(iz,igas) .gt. sifreeze)then
    
                    ! Calculate mean ice density and latent heat of freezing over temperature
                    ! interval [T0,T]
        
                    rhoibar = ( 0.916_f * (t(iz)-T0) - &
                      1.75e-4_f/2._f * ((t(iz)-T0)**2) - &
                      5.e-7_f * ((t(iz)-T0)**3)/3._f ) / (t(iz)-T0)
        
                    rlhbar = ( 79.7_f * (t(iz)-T0) + &
                      0.485_f/2._f * (t(iz)-T0)**2 - &
                      2.5e-3_f/3._f * (t(iz)-T0)**3 ) &
                      / (t(iz)-T0) * 4.186e7*18._f
        
                    ! Equilibrium H2SO4 weight percent for fixed water activity
                    act = min(1.0_f, ssl + 1._f)

                    ! Kelvin effect on water activity
                    fkelv = exp(akelvin(iz,igas) / r(ibin,igroup))                                ! ?
                    act   = act / fkelv

                    IF(act .LT. 0.05_f) THEN
                      CONTL = 12.37208932_f * (act**(-0.16125516114_f)) - &
                              30.490657554_f * act  - 2.1133114241_f
                      CONTH = 13.455394705_f * (act**(-0.1921312255_f))  - &
                              34.285174604_f * act - 1.7620073078_f
                    END IF
                    IF(act .GE. 0.05_f .and. act .LE. 0.85_f) THEN
                      CONTL = 11.820654354_f * (act**(-0.20786404244_f)) - &
                              4.807306373_f * act  - 5.1727540348_f
                      CONTH = 12.891938068_f * (act**(-0.23233847708_f)) - &
                              6.4261237757_f * act - 4.9005471319_f
                    END IF
                    IF(act .GT. 0.85_f) THEN
                      CONTL = -180.06541028_f * (act**(-0.38601102592_f)) - &
                              93.317846778_f * act + 273.88132245_f
                       CONTH = -176.95814097_f * (act**(-0.36257048154_f)) - &
                              90.469744201_f * act + 267.45509988_f
                    END IF
                    H2SO4m = CONTL + ((CONTH - CONTL) * (t(iz) -190._f)/70._f)
                    WT = (98.0_f * H2SO4m)/(1000._f + 98._f * H2SO4m)
                    WT = 100._f * WT
        
                    ! Volume ratio of wet/dry aerosols.
                    vrat = rhosol(isol)/RHO_W * ((100._f-wt)/wt) + 1._f
        
                    ! Calculation sulfate solution density from Myhre et al. (1998).
                    wtfrac = WT/100._f
                    C1      = t(iz) - 273.15_f
                    C2      = C1**2
                    C3      = C1**3
                    C4      = C1**4
                    A0 = 999.8426_f + 334.5402e-4_f*C1 - 569.1304e-5_f*C2
                    A1 = 547.2659_f - 530.0445e-2_f*C1 + 118.7671e-4_f*C2 &
                         + 599.0008e-6_f*C3
                    A2 = 526.295e+1_f + 372.0445e-1_f*C1 + 120.1909e-3_f*C2 &
                         - 414.8594e-5_f*C3 + 119.7973e-7_f*C4
                    A3 = -621.3958e+2_f - 287.7670_f*C1 - 406.4638e-3_f*C2 &
                         + 111.9488e-4_f*C3 + 360.7768e-7_f*C4
                    A4 = 409.0293e+3_f + 127.0854e+1_f*C1 + 326.9710e-3_f*C2 &
                         - 137.7435e-4*C3 - 263.3585e-7*C4
                    A5 = -159.6989e+4_f - 306.2836e+1_f*C1 + 136.6499e-3_f*C2 &
                         + 637.3031e-5_f*C3
                    A6 = 385.7411e+4_f + 408.3717e+1_f*C1 - 192.7785e-3_f*C2
                    A7 = -580.8064e+4_f - 284.4401e+1_f*C1
                    A8 = 530.1976e+4_f + 809.1053_f*C1
                    A9 = -268.2616e+4_f
                    A10 = 576.4288e+3_f
                    den = A0 + wtfrac*A1 + wtfrac**2 * A2 + &
                          wtfrac**3 * A3 + wtfrac**4 * A4
                    den = den + wtfrac**5 * A5 + wtfrac**6 * A6 + &
                          wtfrac**7 * A7
                    den = den + wtfrac**8 * A8 + wtfrac**9 * A9 + &
                          wtfrac**10 * A10
        
                    ! Activation energy is based on Koop's lab data.
                    IF(t(iz) .GT. 220._f) then
                      A0      = 104525.93058_f
                      A1      = -1103.7644651_f
                      A2      = 1.070332702_f
                      A3      = 0.017386254322_f
                      A4      = -1.5506854268e-06_f
                      A5      = -3.2661912497e-07_f
                      A6      = 6.467954459e-10_f
                    ELSE
                      A0      = -17459.516183_f
                      A1      = 458.45827551_f
                      A2      = -4.8492831317_f
                      A3      = 0.026003658878_f
                      A4      = -7.1991577798e-05_f
                      A5      = 8.9049094618e-08_f
                      A6      = -2.4932257419e-11_f
                    END IF
                    
                    diffact = ( A0 + A1*t(iz) + A2*t(iz)**2 + &
                        A3*t(iz)**3 + A4*t(iz)**4 + &
                        A5*t(iz)**5 + A6*t(iz)**6 ) * 1.0e-13_f
        
                    ! Surface energy
        
                    ! Weight percent function for T = 260 K
                    c0      = 77.40682664_f
                    c1      = -0.006963123274_f
                    c2      = -0.009682499074_f
                    c3      = 0.00088797988_f
                    c4      = -2.384669516e-05_f
                    c5      = 2.095358048e-07_f
                    S260 = c0 + c1*wt + c2*wt**2 + c3*wt**3 + &
                           c4*wt**4 + c5*wt**5
        
                    ! Weight percent function for T = 220 K
                    d0      = 82.01197792_f
                    d1      = 0.5312072092_f
                    d2      = -0.1050692123_f
                    d3      = 0.005415260617_f
                    d4      = -0.0001145573827_f
                    d5      = 8.969257061e-07_f
                    S220 = d0 + d1*wt + d2*wt**2 + d3*wt**3 + &
                           d4*wt**4 + d5*wt**5
        
                    ! Weight percent function for T = 180K
                    e0      = 85.75507114_f
                    e1      = 0.09541966318_f
                    e2      = -0.1103647657_f
                    e3      = 0.007485866933_f
                    e4      = -0.0001912224154_f
                    e5      = 1.736789787e-06_f
                    S180 = e0 + e1*wt + e2*wt**2 + e3*wt**3 + &
                           e4*wt**4 + e5*wt**5
        
                    if( t(iz) .GE. 220._f ) then
                      sigma = S260 + ((260._f-t(iz))*(S220-S260))/40._f
                    else
                      sigma = S220 + ((220._f-t(iz))*(S180-S220))/40._f
                    endif
        
                    sigsula = sigma
                    sigicea = 105._f
                    sigsulice = abs( sigsula - sigicea )
        
                    ! Critical ice germ radius formed in the sulfate solution
                    ag = 2._f*gwtmol(igas)*sigsulice / &
                          ( rlhbar * rhoibar * log(T0/t(iz)) + &
                          rhoibar * rgas * 0.5_f * (T0+t(iz)) * &
                         log(ssl+1._f) )
                          
                    if( ag .lt. 0._f ) ag = 1.e10_f
        
                    ! Gibbs free energy of ice germ formation in the ice/sulfate solution
                    delfg = 4._f/3._f*PI * sigsulice * (ag**2)
        
                    ! Ice nucleation rate in a 0.2 micron aerosol (/sec)
                    expon = ( -diffact - delfg ) / BK / t(iz)
                    expon = max( -100._f*ONE, expon )
                    rnuclg(ibin,igroup,ignucto) = prenuc * &
                           sqrt(sigsulice*t(iz)) * &
                           vrat*vol(ibin,igroup) * exp( expon )
                           
                    ! This parameterizations has problems that sometimes yield negative nucleation
                    ! rates. It would be best to fix the parameterization, but at least keep negative
                    ! values from being return.
                    if (rnuclg(ibin,igroup,ignucto) < 0._f) then
                      rnuclg(ibin,igroup,ignucto) = 0._f
                    end if
                  
                  
    !             xh = 0.1 * r(ibin,igroup) / ag
    !             phih = sqrt( 1. - 2.*rmiv*xh + xh**2 )
    !             rath = (xh-rmiv) / phih
    !             fv3h = xh**3 * ( 2.*ONE - 3.*rath + rath**3 )
    !             fv4h = 3. * rmiv * xh**2 * (rath-1.)
    !             if( abs(rath) .gt. 1.e0-1.e-8 ) fv3h = 0.
    !             if( abs(rath) .gt. 1.e0-1.e-10 ) fv4h = 0.
    ! 
    !             fh = 0.5 * ( ONE + ((ONE-rmiv*xh)/phih)**3 +
    !     $                              fv3h + fv4h )
    !
    !             expon = ( -delfwat2ice - delfg ) / BK / t3(ixyz)
    !             expon = max( -POWMAX, expon )
                  endif
                endif   ! pconmax(ixyz,igroup) .gt. FEW_PC 
              enddo      ! ibin = 1,NBIN
            endif       ! inucproc(iepart,ienucto) .eq. I_DROPACT
          endif
        enddo        ! inuc = 1,nnuc2elem(iepart)
!     endif       ! (nnuc2elem(iepart) .gt. 1)
    endif        ! (igas = inucgas(igroup) .ne. 0)
  enddo         ! igroup = 1,NGROUP

  ! Return to caller with particle loss rates due to nucleation evaluated.
  return
end
