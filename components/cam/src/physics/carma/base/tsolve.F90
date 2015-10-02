! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates new potential temperature concentration 
!! (and updates temperature) due to microphysical and radiative forcings.
!! The equation solved (the first law of thermodynamics in flux form) is
!!
!!   d(rhostar theta)     rhostar theta       d(qv)    1  dF
!!   ---------------  = - ------------- * ( L ----- + --- -- )
!!         dt                 Cp T              dt    rho dz
!!
!! where
!!   rhostar = scaled air density
!!   theta   = potential temperature
!!   t       = time
!!   Cp      = specific heat (at constant pressure) of air
!!   T       = air temperature
!!   qv      = water vapor mixing ratio
!!   L       = latent heat
!!   F       = net radiative flux
!!   z       = unscaled altitude
!!
!! @author Andy Ackerman
!! @version Oct-1997
subroutine tsolve(carma, cstate, iz, scale_threshold, rc)

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
  real(kind=f)                         :: scale_threshold !! Scaling factor for convergence thresholds
  integer, intent(inout)               :: rc      !! return code, negative indicates failure
  
  1 format(/,'tsolve::ERROR - negative temperature for : iz=',i4,',lat=',&
              f7.2,',lon=',f7.2,',T=',e10.3,',dT=',e10.3,',t_old=',e10.3,',d_gc=',e10.3,',dT_adv=',e10.3)
  2 format(/,'tsolve::ERROR - temperature change to large for : iz=',i4,',lat=',&
              f7.2,',lon=',f7.2,',T=',e10.3,',dT_rlh=',e10.3,',dT_pth=',e10.3,',t_old=',e10.3,',d_gc=',e10.3,',dT_adv=',e10.3)
  3 format(/,'tsolve::ERROR - temperature change to large for : iz=',i4,',lat=',&
              f7.2,',lon=',f7.2,',T=',e10.3,',dT_rlh=',e10.3,',dT_pth=',e10.3,',t_old=',e10.3)
  4 format(/,'tsolve::ERROR - negative temperature for : iz=',i4,',lat=',&
              f7.2,',lon=',f7.2,',T=',e10.3,',dT=',e10.3,',t_old=',e10.3)
      
  real(kind=f)      :: dt           ! delta temperature
  real(kind=f)      :: threshold    ! convergence threshold
  
      
  ! Solve for the new <t> due to latent heat exchange and radiative heating.
  ! Latent and radiative heating rates are in units of [deg_K/s].
  !
  ! NOTE: In the embedded model rhoa and p are handled by the parent model and
  ! won't change during one time step.
  !
  ! NOTE: Radiative heating by the particles is handled by the parent model, so
  ! that term does not need to be added here.
  dt         = dtime * rlprod
  rlheat(iz) = rlheat(iz) + rlprod * dtime   
  
  ! With particle heating, you must also include the impact of heat
  ! conduction from the particle
  !
  ! NOTE: We are ignoring the energy to heat the particle, since we
  ! are not tracking the particle temperature. Thus ...
  if (do_pheatatm) then
    dt           = dt + dtime * phprod
    partheat(iz) = partheat(iz) + phprod * dtime
  end if
  
  t(iz) = t(iz) + dt
 
  
  ! Don't let the temperature go negative.
  if (t(iz) < 0._f) then
    if (do_substep) then
      if (nretries == maxretries) then 
        if (do_print) write(LUNOPRT,1) iz, lat, lon, t(iz), dt, told(iz), d_gc(iz, 1), d_t(iz)
      end if
    else
      if (do_print) write(LUNOPRT,4) iz, lat, lon, t(iz), dt, told(iz)
    end if
    
    rc = RC_WARNING_RETRY
  end if

  ! Don't let the temperature change by more than the threshold in any given substep,
  ! to prevent overshooting that doesn't result in negative gas concentrations, but
  ! does result in excessive temperature swings.
  threshold = dt_threshold / scale_threshold
  
  if (threshold /= 0._f) then
    if (abs(abs(dt)) > threshold) then
      if (do_substep) then
        if (nretries == maxretries) then 
          if (do_print) write(LUNOPRT,2) iz, lat, lon, t(iz), rlprod*dtime, dtime*partheat(iz), told(iz), d_gc(iz, 1), d_t(iz)
        end if
      else
        if (do_print) write(LUNOPRT,3) iz, lat, lon, t(iz), rlprod*dtime, dtime*partheat(iz), told(iz)
      end if
  
      rc = RC_WARNING_RETRY
    end if
  end if

  ! Return to caller with new temperature.
  return
end
