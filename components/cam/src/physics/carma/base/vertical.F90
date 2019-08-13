! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine drives the vertical transport calculations.
!!
!!  NOTE: Since this is only for sedimentation and brownian diffusion of a column within
!! a parent model, the advection of air density, gases and potential temperature have
!! been removed. Also, the divergence corrections (divcor) for 1D transport are not
!! applied, since these columns exist within a parent model that is responsible for the
!! advection.
!!
!! @author Eric Jensen
!! version Mar-1995
subroutine vertical(carma, cstate, rc)

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

  ! Declare local variables
  integer        :: ielem
  integer        :: ibin
  integer        :: ig
  real(kind=f)   :: vertadvu(NZP1)
  real(kind=f)   :: vertadvd(NZP1)
  real(kind=f)   :: vertdifu(NZP1)
  real(kind=f)   :: vertdifd(NZP1)
  real(kind=f)   :: vtrans(NZP1)
  real(kind=f)   :: old_pc(NZ)
  
  rc = RC_OK
  
  do ielem = 1,NELEM          ! Loop over particle elements
    ig = igelem(ielem)        ! particle group
    
    ! Should this group participate in sedimentation?
    if (grp_do_vtran(ig)) then
    
      ! Are there enough particles in the column to bother?
      if (maxval(pconmax(:,ig)) .gt. FEW_PC) then

        do ibin = 1,NBIN          ! Loop over particle mass bins          
          vtrans(:) = -vf(:,ibin,ig)
    
          ! If dry deposition is enabled for this group, then set
          ! the deposition velocity at the surface.
          if (grp_do_drydep(ig)) then
            if (igridv .eq. I_CART) then
              vtrans(1) = -vd(ibin, ig)
            else
              vtrans(NZP1) = -vd(ibin, ig)
            end if
          end if
          
          !  Calculate particle transport rates due to vertical advection
          !  and vertical diffusion, and solve for concentrations at end of time step.
          call vertadv(carma, cstate, vtrans, pc(:,ibin,ielem), itbnd_pc, ibbnd_pc, &
            pc_topbnd(ibin,ielem), pc_botbnd(ibin,ielem), vertadvu, vertadvd, rc)
          if (rc < RC_OK) return
            
          call vertdif(carma, cstate, ig, ibin, itbnd_pc, ibbnd_pc, vertdifu, vertdifd, rc)
          if (rc < RC_OK) return
   
          old_pc(:) = pc(:,ibin,ielem)
          
          ! There are 2 different solvers, versol with uses a PPM scheme and versub
          ! which using an explicit substepping approach.
          if (do_explised) then
            call versub(carma, cstate, pconmax(:,ig)*xmet(:)*ymet(:)*zmet(:), pc(:,ibin,ielem), itbnd_pc, ibbnd_pc, &
              ftoppart(ibin,ielem), fbotpart(ibin,ielem), &
              pc_topbnd(ibin,ielem), pc_botbnd(ibin,ielem), &
              vertadvu, vertadvd, vertdifu, vertdifd, rc)
            if (rc < RC_OK) return
          else
            call versol(carma, cstate, pc(:,ibin,ielem), itbnd_pc, ibbnd_pc, &
              ftoppart(ibin,ielem), fbotpart(ibin,ielem), &
              pc_topbnd(ibin,ielem), pc_botbnd(ibin,ielem), &
              vertadvu, vertadvd, vertdifu, vertdifd, rc)
            if (rc < RC_OK) return
          end if
          
          ! A clunky way to get the mass flux to the surface and to conserve mass
          ! is to determine the total before and after. Anything lost went to the
          ! surface.
          !
          ! NOTE: This only works if you assume nothing is lost out the top. It would be
          ! better to figure out how to get this directly from versol.
          pc_surf(ibin,ielem) = pc_surf(ibin, ielem) + sum(old_pc(:) * dz(:) / xmet(:) / ymet(:)) - &
            sum(pc(:,ibin,ielem) * dz(:) / xmet(:) / ymet(:))
          sedimentationflux(ibin,ielem) = ( sum(old_pc(:) * dz(:) / xmet(:) / ymet(:)) - &
            sum(pc(:,ibin,ielem) * dz(:) / xmet(:) / ymet(:)) ) / dtime
        enddo  ! ibin
      endif
    endif
  enddo  ! ielem

  ! Return to caller with new particle concentrations.
  return
end
