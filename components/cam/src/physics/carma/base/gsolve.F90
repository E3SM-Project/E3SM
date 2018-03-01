! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine calculates new gas concentrations.
!!
!! @author Andy Ackerman, Bill McKie, Chuck Bardeen
!! @version Dec-1995, Sep-1997, Nov-2009
subroutine gsolve(carma, cstate, iz, previous_ice, previous_liquid, scale_threshold, rc)

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
  real(kind=f), intent(in)             :: previous_ice(NGAS)      !! total ice at the start of substep
  real(kind=f), intent(in)             :: previous_liquid(NGAS)   !! total liquid at the start of substep
  real(kind=f)                         :: scale_threshold !! Scaling factor for convergence thresholds
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local Variables
  integer                              :: igas    !! gas index
  real(kind=f)                         :: gc_cgs
  real(kind=f)                         :: rvap
  real(kind=f)                         :: total_ice(NGAS)      ! total ice
  real(kind=f)                         :: total_liquid(NGAS)   ! total liquid
  real(kind=f)                         :: threshold            ! convergence threshold
  
  
  1 format(/,'gsolve::ERROR - negative gas concentration for ',a,' : iz=',i4,',lat=', &
              f7.2,',lon=',f7.2,',gc=',e10.3,',gasprod=',e10.3,',supsati=',e10.3, &
              ',supsatl=',e10.3,',t=',f6.2)
  2 format('gsolve::ERROR - conditions at beginning of the step : gc=',e10.3,',supsati=',e17.10, &
              ',supsatl=',e17.10,',t=',f6.2,',d_gc=',e10.3,',d_t=',f6.2)
  3 format(/,'microfast::WARNING - gas concentration change exceeds threshold: ',a,' : iz=',i4,',lat=', &
              f7.2,',lon=',f7.2, ', (gc-gcl)/gcl=', e10.3)
  

  ! Determine the total amount of condensate for each gas.
  call totalcondensate(carma, cstate, iz, total_ice, total_liquid, rc)
  
  do igas = 1,NGAS
  
    ! We do not seem to be conserving mass and energy, so rather than relying upon gasprod
    ! and rlheat, recalculate the total change in condensate to determine the change
    ! in gas and energy.
    !
    ! This is because in the old scheme, the particles were solved for implicitly, but the
    ! gas and latent heat were solved for explicitly using the same rates.
    gasprod(igas) = ((previous_ice(igas) - total_ice(igas)) + (previous_liquid(igas) - total_liquid(igas))) / dtime
    rlprod        = rlprod - ((previous_ice(igas) - total_ice(igas)) * (rlhe(iz,igas) + rlhm(iz,igas)) + &
                       (previous_liquid(igas) - total_liquid(igas)) * (rlhe(iz,igas))) / (CP * rhoa(iz) * dtime)     

    ! Don't let the gas concentration go negative.
    gc(iz,igas) = gc(iz,igas) + dtime * gasprod(igas)

    if (gc(iz,igas) < 0.0_f) then
      if (do_substep) then
        if (nretries == maxretries) then 
          if (do_print) write(LUNOPRT,1) trim(gasname(igas)), iz, lat, lon, gc(iz,igas), gasprod(igas), &
            supsati(iz,igas), supsatl(iz,igas), t(iz)
          if (do_print) write(LUNOPRT,2) gcl(iz,igas), supsatiold(iz,igas), supsatlold(iz,igas), told(iz), d_gc(iz,igas), d_t(iz)
        end if
      else
        if (do_print) write(LUNOPRT,1) trim(gasname(igas)), iz, lat, lon, gc(iz,igas), gasprod(igas), &
          supsati(iz, igas), supsatl(iz,igas), t(iz)
      end if

      rc = RC_WARNING_RETRY
    end if
    
    ! If gas changes by too much, then retry the calculation.
    threshold = dgc_threshold(igas) / scale_threshold
    
    if (threshold /= 0._f) then
      if ((dtime * gasprod(igas) / gc(iz,igas)) > threshold) then
        if (do_substep) then
          if (nretries == maxretries) then 
            if (do_print) write(LUNOPRT,3) trim(gasname(igas)), iz, lat, lon, dtime * gasprod(igas) / gc(iz,igas)
            if (do_print) write(LUNOPRT,2) gcl(iz,igas), supsatiold(iz,igas), supsatlold(iz,igas), told(iz), d_gc(iz,igas), d_t(iz)
          end if
        else
          if (do_print) write(LUNOPRT,3) trim(gasname(igas)), iz, lat, lon, dtime * gasprod(igas) / gc(iz,igas)
        end if
  
        rc = RC_WARNING_RETRY
      end if
    end if
  end do

  ! Return to caller with new gas concentrations.
  return
end
