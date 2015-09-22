! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine solves for sedimentation using an explicit substepping approach. It
!!  is faster and handles large cfl and irregular grids better than the normal PPM
!!  solver (versol), but it is more diffusive.
!!
!! @author Andy Ackerman, Chuck Bardeen 
!! version Aug 2010
subroutine versub(carma, cstate, pcmax, cvert, itbnd, ibbnd, ftop, fbot, cvert_tbnd, cvert_bbnd, &
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
  real(kind=f), intent(in)             :: pcmax(NZ)       !! maximum particle concentration (#/x/y/z)
  real(kind=f), intent(inout)          :: cvert(NZ)       !! quantity being transported (#/x/y/z)
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
  integer                        :: iz
  integer                        :: istep
  integer                        :: nstep_sed
  real(kind=f)                   :: fvert(NZ)
  real(kind=f)                   :: up(NZP1)
  real(kind=f)                   :: dn(NZP1)
  real(kind=f)                   :: cfl_max
  real(kind=f)                   :: fvert_1
  real(kind=f)                   :: fvert_nz

  ! Determine the total upward and downward velocities.
  up(:) = vertadvu(:) + vertdifu(:)
  dn(:) = vertadvd(:) + vertdifd(:)

  ! Compute the maximum CFL for each bin that has a significant concentration
  ! of particles.
  cfl_max = 0._f

  do iz = 1, NZ
    if (pcmax(iz) > SMALL_PC) then
      cfl_max = max(cfl_max, max(abs(up(iz)), abs(up(iz+1)), abs(dn(iz)), abs(dn(iz+1))) * dtime / dz(iz))
    end if
  end do

  ! Use the maximum CFL determined above to figure out how much substepping is
  ! needed to sediment explicitly without violating the CFL anywhere in the column.
  if (cfl_max >= 0._f) then
    nstep_sed = int(1._f + cfl_max)
  else
    nstep_sed = 0
  endif
  
  ! If velocities are in both directions, then more steps are needed to make sure
  ! that no more than half of the concentration can be transported in either direction.
  if (maxval(up(:) * dn(:)) > 0._f) then
    nstep_sed = nstep_sed * 2
  end if

  ! Determine the top and bottom boundary fluxes, keeping in mind that
  ! the velocities and grid coordinates are reversed in sigma or hybrid
  ! coordinates
  if ((igridv .eq. I_SIG) .or. (igridv .eq. I_HYBRID)) then
    if (itbnd .eq. I_FLUX_SPEC) then
      fvert_nz = -fbot
    else
      fvert_nz = cvert_bbnd*dn(NZ+1)
    end if

    if (ibbnd .eq. I_FLUX_SPEC) then
      fvert_1 = -ftop
    else
      fvert_1 = cvert_tbnd*up(1)
    end if
  
  else
    if (itbnd .eq. I_FLUX_SPEC) then
      fvert_nz = ftop
    else
      fvert_nz = cvert_tbnd*dn(NZ+1)
    end if

    if (ibbnd .eq. I_FLUX_SPEC) then
      fvert_1 = fbot
    else
      fvert_1 = cvert_bbnd*up(1)
    end if
  endif

  ! Sediment the particles using multiple iterations to satisfy the CFL.
  do istep = 1, nstep_sed

    ! Determine the net particle flux at each gridbox. The first and last levels
    ! need special treatment to handle to bottom and top boundary conditions.
    fvert(1) = (-cvert(1)*dn(1) + fvert_1 + cvert(2)*dn(2) - cvert(1)*up(2))
    
    do iz = 2, NZ-1
      fvert(iz) = (-cvert(iz)*dn(iz) + cvert(iz-1)*up(iz) + cvert(iz+1)*dn(iz+1) - cvert(iz)*up(iz+1))
    end do

    fvert(NZ) = (-cvert(NZ)*dn(NZ) + cvert(NZ-1)*up(NZ) + fvert_nz - cvert(NZ)*up(NZ+1))
    
    ! Now update the actual concentrations.
    cvert(:) = cvert(:) + fvert(:) * dtime / nstep_sed / dz(:)
  enddo

  return
end subroutine versub
