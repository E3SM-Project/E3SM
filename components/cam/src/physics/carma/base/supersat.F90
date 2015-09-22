! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine evaluates supersaturations <supsatl> and <supsati> for all gases.
!!
!! @author Andy Ackerman, Chuck Bardeen
!! @version Dec-1995, Aug-2010
subroutine supersat(carma, cstate, iz, igas, rc)

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
  integer, intent(in)                  :: igas    !! gas index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  real(kind=f)  :: rvap
  real(kind=f)  :: gc_cgs
  real(kind=f)  :: alpha

  ! Calculate vapor pressures.
  call vaporp(carma, cstate, iz, igas, rc)

  ! Define gas constant for this gas
  rvap = RGAS / gwtmol(igas)

  gc_cgs = gc(iz,igas) / (zmet(iz)*xmet(iz)*ymet(iz))

  supsatl(iz,igas) = (gc_cgs * rvap * t(iz) - pvapl(iz,igas)) / pvapl(iz,igas)
  supsati(iz,igas) = (gc_cgs * rvap * t(iz) - pvapi(iz,igas)) / pvapi(iz,igas)

  ! For subgrid scale clouds, the supersaturation needs to be increased be scaled
  ! based upon cloud fraction. This approach is similar to Wilson and Ballard (1999),
  ! except that only the water vapor (no liquid water) is used to determine the available
  ! water.
  !
  ! NOTE: This assumes that the cloud is an ice cloud.
  if (do_incloud) then
    alpha = rhcrit(iz) * (1._f - cldfrc(iz)) + cldfrc(iz)
    
    supsatl(iz,igas) = (gc_cgs * rvap * t(iz) - alpha * pvapl(iz,igas)) / pvapl(iz,igas)
    supsati(iz,igas) = (gc_cgs * rvap * t(iz) - alpha * pvapi(iz,igas)) / pvapi(iz,igas)
    
    ! Limit supersaturation to liquid saturation.
    supsatl(iz,igas) = min(supsatl(iz,igas), 0._f)
    supsati(iz,igas) = min(supsati(iz,igas), (pvapl(iz,igas) - alpha * pvapi(iz,igas)) / pvapi(iz,igas))        
  end if

  return
end


!! This routine evaluates supersaturations <supsatl> and <supsati> for all gases, but
!! thus version of the routine does not scale the supersaturation based on the cloud
!! fraction. It also assumes that vaporp has already been called.
!!
!! @author Andy Ackerman, Chuck Bardeen
!! @version Dec-1995, Aug-2010
subroutine supersat_nocldf(carma, cstate, iz, igas, ssi, ssl, rc)

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
  integer, intent(in)                  :: igas    !! gas index
  real(kind=f), intent(out)            :: ssl
  real(kind=f), intent(out)            :: ssi
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  real(kind=f)  :: rvap
  real(kind=f)  :: gc_cgs
  real(kind=f)  :: alpha

  ! Calculate vapor pressures.
  call vaporp(carma, cstate, iz, igas, rc)

  ! Define gas constant for this gas
  rvap = RGAS / gwtmol(igas)

  gc_cgs = gc(iz,igas) / (zmet(iz)*xmet(iz)*ymet(iz))

  ssl = (gc_cgs * rvap * t(iz) - pvapl(iz,igas)) / pvapl(iz,igas)
  ssi = (gc_cgs * rvap * t(iz) - pvapi(iz,igas)) / pvapi(iz,igas)

  return
end
