! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine solves the vertical transport equation.
!!  <cvert> is temporary storage for concentrations (particles,
!!  gas, potential temperature) being transported.
!!  New values of <cvert> are calculated.
!!
!! @author Eric Jensen
!! @version Dec-1996
subroutine versol(carma, cstate, cvert, itbnd, ibbnd, ftop, fbot, cvert_tbnd, cvert_bbnd, &
  vertadvu, vertadvd, vertdifu, vertdifd, rc)

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
  real(kind=f), intent(inout)          :: cvert(NZ)       !! quantity being transported
  integer, intent(in)                  :: itbnd           !! top boundary condition
  integer, intent(in)                  :: ibbnd           !! bottom boundary condition
  real(kind=f), intent(in)             :: ftop            !! flux at top boundary
  real(kind=f), intent(in)             :: fbot            !! flux at bottom boundary
  real(kind=f), intent(in)             :: cvert_tbnd      !! quantity at top boundary
  real(kind=f), intent(in)             :: cvert_bbnd      !! quantity at bottom boundary
  real(kind=f), intent(in)             :: vertadvu(NZP1)  !! upward vertical transport rate into level k from level k-1 [cm/s]
  real(kind=f), intent(in)             :: vertadvd(NZP1)  !! downward vertical transport rate into level k from level k-1 [cm/s]
  real(kind=f), intent(in)             :: vertdifu(NZP1)  !! upward vertical diffusion rate into level k from level k-1 [cm/s]
  real(kind=f), intent(in)             :: vertdifd(NZP1)  !! downward vertical diffusion rate into level k from level k-1 [cm/s]
  integer, intent(inout)               :: rc              !! return code, negative indicates failure

!  Declare local variables
!
  integer                        :: k
  real(kind=f)                   :: al(NZ)
  real(kind=f)                   :: bl(NZ)
  real(kind=f)                   :: dl(NZ)
  real(kind=f)                   :: el(NZ)
  real(kind=f)                   :: fl(NZ)
  real(kind=f)                   :: ul(NZ)
  real(kind=f)                   :: ctempl(NZ)
  real(kind=f)                   :: ctempu(NZ)
  real(kind=f)                   :: divcor(NZ)
  real(kind=f)                   :: uc
  real(kind=f)                   :: cour
  real(kind=f)                   :: denom
     
  ! Divergence adjustments are not being generated.
  divcor(:) = 0._f

  ! Determine whether transport should be solved explicitly (uc=0)
  ! or implicitly (uc=1).
  uc = 0._f
  do k = 1,NZ
    cour = dz(k)/dtime - &
    ( vertdifu(k+1) + vertdifd(k) + vertadvu(k+1) + vertadvd(k) )

    if( cour .lt. 0._f .and. uc .ne. 1._f )then
      uc = 1.0_f

      ! NOTE: This can happen a lot and clutters up the log. Should we print it out or not?     
!      write(LUNOPRT,'(a,i3,7(1x,1pe8.1))') &
!         'in versol: k dz/dt vdifd vdifu vadvd vadvu cour uc = ', &
!          k, dz(k)/dtime, vertdifd(k), vertdifu(k+1), &
!          vertadvd(k), vertadvu(k+1), cour, uc
    endif
  enddo

  ! Store concentrations in local variables (shifted up and down
  ! a vertical level).
  do k = 2,NZ
    ctempl(k) = cvert(k-1)
    ctempu(k-1) = cvert(k)
  enddo

  if( ibbnd .eq. I_FIXED_CONC ) then
    ctempl(1) = cvert_bbnd
  else
    ctempl(1) = 0._f
  endif

  if( itbnd .eq. I_FIXED_CONC ) then
    ctempu(NZ) = cvert_tbnd
  else
    ctempu(NZ) = 0._f
  endif
  
  ! Calculate coefficients of the transport equation:
  !   al(k)c(k+1) + bl(k)c(k) + ul(k)c(k-1) = dl(k)

  do k = 1,NZ
    al(k) = uc * ( vertdifd(k+1) + vertadvd(k+1) )
    bl(k) = -( uc*(vertdifd(k)+vertdifu(k+1)+ &
                   vertadvd(k)+vertadvu(k+1)) &
               + dz(k)/dtime )
    ul(k) = uc * ( vertdifu(k) + vertadvu(k) )
    dl(k) = cvert(k) * &
      ( (1._f - uc)*(vertdifd(k)+vertdifu(k+1)+ &
                 vertadvd(k)+vertadvu(k+1)) &
      - dz(k)/dtime ) - &
      (1._f - uc) * ( (vertdifu(k)+vertadvu(k))*ctempl(k) + &
      (vertdifd(k+1)+vertadvd(k+1))*ctempu(k) ) - &
      divcor(k) * dz(k)
  enddo

  ! Boundary fluxes: <ftop> is the downward flux across the
  ! upper boundary; <fbot> is the upward flux across the
  ! lower boundary.
  if(( igridv .eq. I_SIG ) .or. ( igridv .eq. I_HYBRID )) then
    if( itbnd .eq. I_FLUX_SPEC ) dl(1) = dl(1) - ftop
    if( ibbnd .eq. I_FLUX_SPEC ) dl(NZ) = dl(NZ) - fbot
  else
    if( itbnd .eq. I_FLUX_SPEC ) dl(NZ) = dl(NZ) - ftop
    if( ibbnd .eq. I_FLUX_SPEC ) dl(1) = dl(1) - fbot
  endif

  ! Calculate recursion relations.
  el(1) = dl(1)/bl(1)
  fl(1) = al(1)/bl(1)
  do k = 2,NZ
    denom = bl(k) - ul(k) * fl(k-1)
    el(k) = ( dl(k) - ul(k)*el(k-1) ) / denom
    fl(k) = al(k) / denom
  enddo

  ! Calculate new concentrations.

  cvert(NZ) = el(NZ)
  do k = NZ-1,1,-1
    cvert(k) = el(k) - fl(k)*cvert(k+1)
  enddo

  ! Return to caller with new concentrations.
  return
end
