! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluate particle loss rates due to particle heating.
!!
!! The net energy absorbed by each particle is calculatated as <qrad>, and
!! this heating rate is included in the caclulation of <dmdt> in growevapl. The
!! particle temperature perturbation realtive to atmospheric temperature <dtpart>
!! and the radiative heating of the atmosphere by particles <partheat>
!! are also calculated.
!!
!! This algorithm is based upon the model described in the appendix of
!!   Toon et al., J. Geophys. Res., 94, 11359-11380, 1989.
!!
!! This routine assumes that the following variable/tables have already been
!! set up:
!!
!!   <radint> intensity of incoming radiance (solar+ir) [erg/cm2/sr/s/cm]
!!   <wave>   wavelengths used for integration [cm]
!!   <dwave>  width of wavelength bands for integration [cm]
!!   <do_wave_emit> whether planck emission should be doen for the band
!!   <qext>   extinction [cm2]
!!   <ssa>    single scattering albedo
!!
!! @author Chuck Bardeen
!! @version Jan-2010
subroutine pheat(carma, cstate, iz, igroup, iepart, ibin, igas, dmdt, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  
  use planck, only         : planckIntensity, planckBandIntensity, planckBandIntensityWidger1976, planckBandIntensityConley2011

  implicit none

    
  type(carma_type), intent(inout)      :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(in)                  :: iz      !! vertical index
  integer, intent(in)                  :: igroup  !! group index
  integer, intent(in)                  :: iepart  !! group's concentration element index
  integer, intent(in)                  :: ibin    !! bin index
  integer, intent(in)                  :: igas    !! gas index
  real(kind=f), intent(out)            :: dmdt    !! particle growth rate (g/s)
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer, parameter                   :: MAX_ITER      = 10      ! Maximum number of iterations
  real(kind=f), parameter              :: DDTP_LIMIT    = 0.01_f   ! Convergence criteria for iteration.
  
  integer                              :: iter                    ! iteration
  integer                              :: iwvl                    ! wavelength band index
  integer                              :: ieother(NELEM)
  integer                              :: nother
  integer                              :: ieoth_rel
  integer                              :: ieoth_abs
  integer                              :: jother
  integer                              :: isol
  real(kind=f)                         :: otherm(NELEM)
  real(kind=f)                         :: argsol
  real(kind=f)                         :: othermtot
  real(kind=f)                         :: condm
  real(kind=f)                         :: akas
  real(kind=f)                         :: expon
  real(kind=f)                         :: g0
  real(kind=f)                         :: g1
  real(kind=f)                         :: g2
  real(kind=f)                         :: ss
  real(kind=f)                         :: pvap
  real(kind=f)                         :: qrad                    ! particle net radiation (erg/s)
!  real(kind=f)                         :: qrad0                   ! particle net radiation (Tp=Ta) (erg/s)
  real(kind=f)                         :: rlh                     ! latent heat (erg/g)
  real(kind=f)                         :: tp                      ! particle temperature (K)
  real(kind=f)                         :: dtp                     ! change in particle temperature (K)
  real(kind=f)                         :: dtpl                    ! last change in particle temperature (K)
  real(kind=f)                         :: ddtp                    ! change in particle temperature in last iteration (K)
  real(kind=f)                         :: plkint                  ! planck intensity
  
  ! <akas> is combined kelvin (curvature) and solute factors.
  !
  ! Ignore solute factor for ice particles.
  if( is_grp_ice(igroup) )then
    expon = akelvini(iz,igas) / rup_wet(iz,ibin,igroup)
  else
  
    argsol = 0._f
  
    ! Consider growth of average particle at radius <rup(ibin,igroup)>.
    ! 
    ! Treat solute effect first: <asol> is solute factor.
    !
    ! Only need to treat solute effect if <nelemg(igroup)> > 1
    if( nelemg(igroup) .gt. 1 )then
  
      ! <condm> is mass concentration of condensed gas <igas> in particle.
      ! <nother> is number of other elements in group having mass.
      ! <otherm> are mass concentrations of other elements in particle group.
      ! <othermtot> is total mass concentrations of other elements in particle.
      nother = 0
      othermtot = 0._f
  
      ! <ieoth_rel> is relative element number of other element in group.
      do ieoth_rel  = 2,nelemg(igroup)       
  
        ! <ieoth_abs> is absolute element number of other element.
        ieoth_abs = iepart + ieoth_rel - 1    
  
        if( itype(ieoth_abs) .eq. I_COREMASS )then
          nother = nother + 1
          ieother(nother) = ieoth_abs
          otherm(nother) = pc(iz,ibin,ieoth_abs)
          othermtot = othermtot + otherm(nother)
        endif
  
      enddo
  
      condm = rmass(ibin,igroup) * pc(iz,ibin,iepart) - othermtot
  
      if( condm .le. 0._f )then
  
        ! Zero mass for the condensate -- <asol> is a small value << 1
        argsol = 1e6_f     
  
      else
  
        ! Sum over masses of other elements in group for argument of solute factor.
        do jother = 1,nother
          isol = isolelem(ieother(jother))
          
          ! Some elements aren't soluble, so skip them.
          if(isol .gt. 0 ) argsol = argsol + sol_ions(isol)*otherm(jother)/solwtmol(isol)
        enddo 
       
        argsol = argsol*gwtmol(igas)/condm
      endif 
    endif    ! nelemg(igroup) > 1

    expon = akelvin(iz,igas)  / rup_wet(iz,ibin,igroup) - argsol 
  endif
  
  expon = max(-POWMAX, expon)
  akas  = exp( expon )

  ! Trick for removing haze droplets from droplet bins:
  ! allows haze droplets to exist under supersaturated conditions;
  ! when below supersaturation, haze droplets will evaporate.
!          if( (.not. is_grp_ice(igroup)) .and. (akas .lt. 1._f) .and. &
!              (supsatl(iz,igas) .lt. 0._f) ) akas = 1._f

  ! <dmdt> is growth rate in mass space [g/s].
  g0 =  gro(iz,ibin+1,igroup)
  g1 = gro1(iz,ibin+1,igroup)
  g2 = gro2(iz,igroup)

  if( is_grp_ice(igroup) )then
    ss   = supsati(iz,igas)
    pvap = pvapi(iz,igas)
  else
    ss = supsatl(iz,igas)
    pvap = pvapl(iz,igas)
  endif


  ! If particle heating is being considered, then determine qrad and tpart to
  ! determine dmdt.
  !
  ! NOTE: If no optical properties, then can't do the particle heating calculation.
  if ((.not. do_pheat) .or. (.not. do_mie(igroup))) then

    ! Ignore the qrad term.
    dmdt = pvap * ( ss + 1._f - akas ) * g0 / ( 1._f + g0 * g1 * pvap )
                     
  else
  
    ! Latent heat of condensing gas 
    if( is_grp_ice(igroup) )then
      rlh = rlhe(iz,igas) + rlhm(iz,igas)
    else
      rlh = rlhe(iz,igas)
    endif
  
    ! The particle temperature must be solved for by iterating, with an
    ! initial guess that the particle temperature is the ambient temperature.
    !
    ! NOTE: We could also try a guest that is based upon an equilibrium
    ! between upwelling IR and collisonal heating, which was identified by
    ! Jensen [1989] as the dominant terms.
    !
    !     radp = 0.d0
    !      
    !     do iwvl = 1, Nwave
    !       radp = radp + (4.0d0*PI * absk(iwvl,ibin+1,igroup) *
    !    $    radint3(ixyz,iwvl) * dwave(iwvl))
    !     end do
    !      
    !     dtp2 = radp /
    !    $  (4.d0*PI*rlow(ibin+1,igroup)*thcondnc(iz)*ft(iz,ibin+1,igroup))
    tp   = t(iz)
    dtp  = 0._f
    dtpl = 0._f
        
    do iter = 1, MAX_ITER
  
      ! Calculate the net radiative flux on the particle, which requires
      ! integrating the incoming and outgoing flux over the spectral
      ! interval.
      qrad = 0._f

      do iwvl = 1, NWAVE

        ! There may be overlap between bands, so only do the emission
        ! for each range of wavelengths once.
        if (do_wave_emit(iwvl)) then
        
          ! Get an integral across the entire band. There are several
          ! techniques for doing this that vary in accuracy and
          ! performance. Comments below are based on the CAM RRTMG
          ! band structure.
          
          ! Just use the band center.
          !
          ! NOTE: This generates about a 20% error, but is the fastest
!         plkint = planckIntensity(wave(iwvl), tp)
          
          ! Brute Force integral
          !
          ! The slowest technique, and not as accurate as either Widger
          ! and Woodall or Conley, even at 100 iterations.
!         plkint = planckBandIntensity(wave(iwvl), dwave(iwvl), tp, 60)
          
          ! Integral using Widger and Woodall, 1976.
          ! 
          ! NOTE: One of the fastest technique at 2 iterations, but yields errors
          ! of about 2%. Can handle wide rage of band sizes.
!          plkint = planckBandIntensityWidger1976(wave(iwvl), dwave(iwvl), tp, 2)

          ! Using method developed by Andrew Conley.
          !
          ! This is similar in performance to Widger and Woodall, but is more
          ! accurate with errors of about 0.3%. It had trouble with SW bands that
          ! are very large, but the latest version has improved performance and
          ! it does work with the RRTMG band structure.
          plkint = planckBandIntensityConley2011(wave(iwvl), dwave(iwvl), tp, 1)
          
        else
          plkint = 0._f
        end if

        qrad = qrad + 4.0_f * PI * (1._f - ssa(iwvl,ibin+1,igroup)) * qext(iwvl,ibin+1,igroup) * PI * (rlow_wet(iz,ibin+1,igroup) ** 2) * arat(ibin+1,igroup) * &
             (radint(iz,iwvl) - plkint) * dwave(iwvl)
      end do
      
      ! Save of the Qrad association with the ambient air temperature.
!      if (iter == 0) then
!        qrad0 = qrad
!      end if

      ! Calculate the change in mass using eq. A3 from Toon et al. [1989].
      dmdt = pvap * ( ss + 1._f - akas * (1._f + qrad * g1 * g2 )) * &
             g0 / ( 1._f + g0 * g1 * pvap )
  
      ! Calculate a new particle temperature based upon the loss of mass and
      ! energy being absorbed.
      if ((dmdt * dtime) .le. (- rmass(ibin+1, igroup))) then
        dtp = ((rlh * (- rmass(ibin+1, igroup) / dtime)) + qrad) / &
               (4._f * PI * rlow_wet(iz,ibin+1,igroup) * thcondnc(iz,ibin+1,igroup) * ft(iz,ibin+1,igroup))
      else
        dtp = ((rlh * dmdt) + qrad) / &
               (4._f * PI * rlow_wet(iz,ibin+1,igroup) * thcondnc(iz,ibin+1,igroup) * ft(iz,ibin+1,igroup))
      end if

      tp = t(iz) + dtp
          
      ddtp = dtp - dtpl
      dtpl = dtp

      if (abs(ddtp) .le. DDTP_LIMIT) then
        exit
      end if
          
      if ((iter .gt. 1) .and. (ddtp .gt. dtpl)) then
        exit
      end if
    end do

    dtpart(iz,ibin,igroup) = dtp
  
    ! Calculate the contribution of this bin to the heating of the atmosphere. CARMA does
    ! not actually apply this heating to change the temperature.
    !
    ! From Pruppacher & Klett [2000], eq. 13-19, the heat transfer to
    ! one particle is:
    !
    !   dq/dt = 4*pi*r*thcondnc*Ft(r)*(T - Tp(r))
    !
    ! so the total heating rate of the air by the particle is:
    !
    !   dT/dt = -Sum((4*pi*r*thcondnc*Ft(r)*(T-Tp(r))*pc(r))) / (Cp,air*arho)
    !
    ! or
    !
    !   dT/dt = Sum((4*pi*r*thcondnc*Ft(r)*dtp*pc(r))) / (Cp,air*arho)
    !
    ! where dtp = Tp(r) - T
    ! 
    ! NOTE: Using these terms will cause the model parent model to go out of
    ! energy balance, since qrad difference is not being communicated to the
    ! other layers.
    if (do_pheatatm) then

      ! NOTE: If the particle is going to evaporate entirely during the timestep,
      ! then assume that there is no particle heating.
      if ((dmdt * dtime) .gt. (- rmass(ibin+1, igroup))) then
    
        ! If the particles are radiatively active, then the parent model's radiation
        ! code is calculated based upon Ta, not Tp. Adjust for this error in Qrad.
!        phprod = phprod + (qrad - qrad0) * pc(iz,ibin+1,iepart) / CP / rhoa(iz)

        ! Now add in the heating from thermal conduction.
        phprod = phprod + 4._f * PI * rlow_wet(iz,ibin+1,igroup) * thcondnc(iz,ibin+1,igroup) * &
                 ft(iz,ibin+1,igroup) * dtp * pc(iz,ibin+1,iepart) / (CP * rhoa(iz))
      end if
    end if
  end if

  !  Return to caller with particle loss rates for growth and evaporation
  !  evaluated.
  return
end
