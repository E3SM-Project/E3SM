! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluates particle loss rates due to nucleation <rnuclg>:
!! aerosol freezing only.
!!
!! The parameterization described by Koop et al., Nature 406, 611-614, 2000
!! is used.
!! 
!! The loss rates for all particle elements in a particle group are equal.
!!
!! To avoid nucleation into an evaporating bin, this subroutine must
!! be called after growp, which evaluates evaporation loss rates <evaplg>.
!!
!! @author Eric Jensen, Chuck Bardeen
!! @version Dec-2003, Apr-2010
subroutine freezaerl_koop2000(carma, cstate, iz, rc)

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
  real(kind=f), parameter ::  prenuc = 2.075e33_f * RHO_W / RHO_I
  real(kind=f), parameter ::  kt0    = 1.6e0_f
  real(kind=f), parameter ::  dkt0dp = -8.8e0_f
  real(kind=f), parameter ::  kti    = 0.22e0_f
  real(kind=f), parameter ::  dktidp = -0.17e0_f
  
  logical                              :: evapfrom_nucto
  integer                              :: igas    ! gas index
  integer                              :: igroup  ! group index
  integer                              :: ibin    ! bin index
  integer                              :: iepart  ! element for condensing group index
  integer                              :: inuc    ! nucleating element index
  integer                              :: isol    ! solute index of freezing particle
  integer                              :: ienucto ! index of target nucleation element
  integer                              :: ignucto ! index of target nucleation group
  integer                              :: inucto  ! index of target nucleation bin
  real(kind=f)                         :: sifreeze
  real(kind=f)                         :: aw
  real(kind=f)                         :: CONTL
  real(kind=f)                         :: CONTH
  real(kind=f)                         :: H2SO4m
  real(kind=f)                         :: WT
  real(kind=f)                         :: volrat
  real(kind=f)                         :: ssi
  real(kind=f)                         :: ssl
  real(kind=f)                         :: rjj
  real(kind=f)                         :: rlogj
  real(kind=f)                         :: daw
  real(kind=f)                         :: riv
  real(kind=f)                         :: vw0
  real(kind=f)                         :: awi
  real(kind=f)                         :: rsi
  real(kind=f)                         :: dmy
  real(kind=f)                         :: rlnt
  real(kind=f)                         :: td
  real(kind=f)                         :: pp
  real(kind=f)                         :: pp2
  real(kind=f)                         :: pp3
  real(kind=f)                         :: vi
  real(kind=f)                         :: fkelv
  real(kind=f)                         :: fkelvi


  rc = RC_OK
  
  !  Aerosol freezing limited to T < 240K
  if (t(iz) <= 240._f) then

    !  Loop over particle groups.
    do igroup = 1,NGROUP

      igas   = inucgas(igroup)
      iepart = ienconc(igroup)
      isol   = isolelem(iepart)

      if (igas .ne. 0) then
      
        !  Bypass calculation if few particles are present
        if (pconmax(iz,igroup) .gt. FEW_PC) then

          ! Calculate nucleation loss rates.  Do not allow nucleation into
          ! an evaporating bin.
          do inuc = 1, nnuc2elem(iepart)

            ienucto = inuc2elem(inuc,iepart)
            if (ienucto /= 0) then
              ignucto = igelem( ienucto )

              ! Only compute nucleation rate for aerosol freezing
              !
              ! NOTE: If heterogeneous nucleation of glassy aerosols is being used
              ! as a nucleation mechanism, then both the heterogeneous nucleation and
              ! the homogeneous freezing need to be considered.
              if ((iand(inucproc(iepart,ienucto), I_AF_KOOP_2000) /= 0)) then
                  
                !  Loop over particle bins.
                do ibin = 1, NBIN
  
                  ssi = supsati(iz,igas)
                  ssl = supsatl(iz,igas)
  
                  ! Calculate approximate critical saturation needed for homogeneous freezing
                  ! of sulfate aerosols (see Jensen and Toon, GRL, 1994).
                  sifreeze = 0.3_f
             
                  ! Homogeneous freezing of sulfate aerosols should only occur if SL < Scrit
                  ! and SI > <sifreeze>.
                  if (ssi > sifreeze) then
  
                    !  Koop et al. nucleation rate parameterization
                    td    = t(iz)
                    rlnt  = log(td)
                    dmy   = 210368._f + 131.438_f * td - (3.32373e6_f / td) - 41729.1_f * rlnt    ! eqn 2, potential difference [J mol-1]
                    rsi   = RGAS / 1.e7_f       ! gas constant [J mol-1 K-1]
                    awi   = exp(dmy / (rsi * td))                                                 ! Notes (p: ambient vs. at pressure) ?
  
                    vw0   = -230.76_f - 0.1478_f * td + (4099.2_f / td) + 48.8341_f * rlnt        ! eqn 4
                    vi    = 19.43_f - 2.2e-3_f * td + 1.08e-5_f * td * td                         ! eqn 5
                    
                    pp    = 1.e-10_f * p(iz)  ! pressure [GPa]
                    pp2   = pp * pp * 0.5_f
                    pp3   = pp2 * pp / 3._f
                    riv   = vw0 * (pp - kt0 * pp2 - dkt0dp * pp3) - vi * (pp - kti * pp2 - dktidp * pp3) ! eqn 3
      
                    riv   = riv * 1.e3_f      ! [GPa cm3 mol-1] to [Pa m3 mol-1] 
  
                    ! NOTE: The wieght percent can become negative from this parameterization,
                    ! which is not physicsal. With small supersaturations, the water activity 
                    ! becomes postive (>1.013) the weight percent becomes negative. Don't allow
                    ! the the supsatl to be greater than 0.
                    ssl = max(-1.0_f, min(0._f, ssl))
  
                    ! Water activity
                    aw    = 1._f + ssl                                                            ! ?

                    ! Kelvin effect on water activity
                    fkelv = exp(akelvin(iz,igas) / r(ibin,igroup))                                ! ?
                    aw    = aw / fkelv

                    ! Nucleation rate
                    !
                    ! NOTE: This formulation is only valid for daw in the range of
                    ! 0.26 < daw < 0.34, so limit daw to that range.
                    daw   = aw * exp(riv / (rsi*td)) - awi                                        ! eqn 6
                    daw = min(0.34_f, max(daw, 0.26_f))                                           ! eqn 7
                    
                    rlogj = ((29180._f * daw - 26924._f) * daw + 8502._f) * daw - 906.7_f         ! eqn 7
                    rlogj = min(rlogj, POWMAX*0.3_f)
                    rjj   = 10._f**(rlogj)                                                        ! [cm-3 s-1]
    

                    ! Calculate volume ratio of wet/dry aerosols
                    if (aw < 0.05_f) then
                      CONTL =  12.37208932_f  * (aw**(-0.16125516114_f)) - 30.490657554_f * aw - 2.1133114241_f
                      CONTH =  13.455394705_f * (aw**(-0.1921312255_f))  - 34.285174604_f * aw - 1.7620073078_f
                    elseif (aw <= 0.85) then
                      CONTL =  11.820654354_f * (aw**(-0.20786404244_f)) - 4.807306373_f  * aw - 5.1727540348_f
                      CONTH =  12.891938068_f * (aw**(-0.23233847708_f)) - 6.4261237757_f * aw - 4.9005471319_f
                    else
                      CONTL = -180.06541028_f * (aw**(-0.38601102592_f)) - 93.317846778_f * aw + 273.88132245_f
                      CONTH = -176.95814097_f * (aw**(-0.36257048154_f)) - 90.469744201_f * aw + 267.45509988_f
                    endif
                    
                    H2SO4m = CONTL + ((CONTH - CONTL) * (t(iz) - 190._f) / 70._f)
                    WT = (98.0_f * H2SO4m) / (1000._f + 98._f * H2SO4m)
                    WT = max(0._f, min(1._f, WT))
                    WT = 100._f * WT
  
                    ! Volume ratio of wet/dry aerosols.
                    if (WT <= 0._f) then
                      volrat = 1.e10_f
                    else
                      volrat = rhosol(isol) / RHO_W * ((100._f - WT) / WT) + 1._f
                    endif

                    rnuclg(ibin,igroup,ignucto) = rnuclg(ibin,igroup,ignucto) + rjj * volrat * vol(ibin,igroup)                ! [s-1]
                  endif   ! ssi > sifreeze .and. target droplets not evaporating
                enddo    ! ibin = 1,NBIN
              endif     ! inucproc(iepart,ienucto) .eq. I_DROPACT
            endif
          enddo      ! inuc = 1,nnuc2elem(iepart)
        endif       ! pconmax .gt. FEW_PC
      endif        ! (igas = inucgas(igroup) .ne. 0)
    enddo         ! igroup = 1,NGROUP
  endif

  ! Return to caller with particle loss rates due to nucleation evaluated.
  return
end subroutine
