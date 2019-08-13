! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates vertrical transport rates.
!! Currently treats diffusion only.
!! Not necessarily generalized for irregular grid.
!!
!! @author Eric Jensen
!! @version Dec-1996
subroutine vertdif(carma, cstate, igroup, ibin, itbnd, ibbnd, vertdifu, vertdifd, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)         :: carma           !! the carma object
  type(carmastate_type), intent(inout) :: cstate          !! the carma state object
  integer, intent(in)                  :: igroup          !! particle group index
  integer, intent(in)                  :: ibin            !! particle bin index
  integer, intent(in)                  :: itbnd           !! top boundary condition
  integer, intent(in)                  :: ibbnd           !! bottom boundary condition
  real(kind=f), intent(out)            :: vertdifu(NZP1)  !! upward vertical diffusion rate into level k from level k-1 [cm/s]
  real(kind=f), intent(out)            :: vertdifd(NZP1)  !! downward vertical diffusion rate into level k from level k-1 [cm/s]
  integer, intent(inout)               :: rc              !! return code, negative indicates failure
  
  ! Local Variables
  integer      :: k
  integer      :: nzm1
  integer      :: itwo
  real(kind=f) :: dz_avg
  real(kind=f) :: rhofact
  real(kind=f) :: xex
  real(kind=f) :: ttheta


  ! Set some constants
  nzm1 = max( 1, NZ-1 )
  itwo = min( 2, NZ   )

  !  Loop over vertical levels.
  do k = 2, NZ

    dz_avg = dz(k)                            ! layer thickness

    !  Check the vertical coordinate

    if( igridv .eq. I_CART ) then
      rhofact = log(  rhoa(k)/rhoa(k-1) &
                    * zmet(k-1)/zmet(k) )
      xex = rhoa(k-1)/rhoa(k) * &
            zmet(k)/zmet(k-1)
      vertdifu(k) = ( rhofact * dkz(k, ibin, igroup) / dz_avg ) / ( 1._f - xex )

      vertdifd(k) = vertdifu(k) * xex


    !  ...else you're in sigma or hybrid coordinates...
    elseif(( igridv .eq. I_SIG ) .or. ( igridv .eq. I_HYBRID )) then
      vertdifu(k) = dkz(k, ibin, igroup) / dz_avg
      vertdifd(k) = dkz(k, ibin, igroup) / dz_avg

    !  ...else write an error (maybe redundant)...
    else
      if (do_print) write(LUNOPRT,*) 'vertdif::ERROR - Invalid vertical grid type (', igridv, ').'
      rc = -1
      return
    endif
  enddo

  ! Fluxes at boundaries specified by user
  if( ibbnd .eq. I_FLUX_SPEC ) then
    vertdifu(1) = 0._f
    vertdifd(1) = 0._f
  endif

  if( itbnd .eq. I_FLUX_SPEC ) then
    vertdifu(NZ+1) = 0._f
    vertdifd(NZ+1) = 0._f
  endif

  ! Diffusion across boundaries using fixed boundary concentration:
  if( ibbnd .eq. I_FIXED_CONC ) then
    dz_avg = dz(1)                            ! layer thickness
    rhofact = log( rhoa(itwo)/rhoa(1) )
    ttheta = rhofact
    if( ttheta .ge. 0._f ) then
      ttheta = min(ttheta,POWMAX)
    else
      ttheta = max(ttheta,-POWMAX)
    endif

    xex = exp(-ttheta)
    if( abs(ONE - xex) .lt. ALMOST_ZERO ) xex = ALMOST_ONE

    vertdifu(1) = ( rhofact * dkz(1, ibin, igroup) / dz_avg ) / ( 1._f - xex )
    vertdifd(1) = vertdifu(1) * xex
  endif

  if( itbnd .eq. I_FIXED_CONC ) then
    dz_avg = dz(NZ)                            ! layer thickness
    rhofact = log( rhoa(NZ)/rhoa(nzm1) )
    ttheta = rhofact
    if( ttheta .ge. 0._f ) then
      ttheta = min(ttheta,POWMAX)
    else
      ttheta = max(ttheta,-POWMAX)
    endif

    xex = exp(-ttheta)
    if( abs(ONE - xex) .lt. ALMOST_ZERO ) xex = ALMOST_ONE

    vertdifu(NZ+1) = ( rhofact * dkz(NZ+1, ibin, igroup) / dz_avg ) / ( 1._f - xex )
    vertdifd(NZ+1) = vertdifu(NZ+1) * xex
  endif

  ! Return to caller with vertical diffusion rates.
  return
end
