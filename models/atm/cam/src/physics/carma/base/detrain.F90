! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine moves condensate from the detrained bins (pcd) to the
!!  particle bins.
!!
!! @author Chuck Bardeen
!! @version May 2010
subroutine detrain(carma, cstate, rc)

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
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer                              :: iz      ! z index
  integer                              :: ibin    ! bin index
  integer                              :: ielem   ! element index

  rc = RC_OK
  
  ! Add the detrained condensate to the particle bins.
  !
  ! NOTE: For now, do this all prior to the fast microphysics, but eventually it may
  ! be better to move it into microfast and substep the detrained condensate.
  pc(:,:,:)  = pc(:,:,:) + pcd(:,:,:)
  pcd(:,:,:) = 0._f

  ! Prevent particle concentrations from dropping below SMALL_PC
  do iz = 1, NZ
    do ibin = 1, NBIN
      do ielem = 1, NELEM
        call smallconc(carma, cstate, iz, ibin, ielem, rc)
      end do
    end do
  end do
  
  !  Return to caller with new particle number concentrations.
  return
end
