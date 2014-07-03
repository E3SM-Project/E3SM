! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  Calculates particle production rates due to nucleation <rhompe>:
!!  binary homogeneous nucleation of sulfuric acid and water only
!!  Numerical method follows Zhao & Turco, JAS, V.26, No.5, 1995.
!!
!!  @author Mike Mills
!!  @version Jul-2001
subroutine sulfnuc(carma,cstate, iz, rc) 
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  use sulfate_utils
  
  implicit none
  
  type(carma_type), intent(in)         :: carma       !! the carma object
  type(carmastate_type), intent(inout) :: cstate      !! the carma state object
  integer, intent(in)                  :: iz          !! level index
  integer, intent(inout)               :: rc          !! return code, negative indicates failure

  !  Local declarations     
  integer           :: igroup     ! group index
  integer           :: ielem      ! concentration element index
  integer           :: nucbin     ! bin in which nucleation takes place    
  real(kind=f)      :: nucrate    ! nucleation rate (#/x/y/z/s)
    

  ! Cycle through each group, only proceed if BHN
  do igroup = 1 , NGROUP
    
    ielem = ienconc(igroup)      
    
    if (inucproc(ielem,ielem) .eq. I_HOMNUC) then
    
      ! This is where all of the pre calculation needs to go, so that it isn't
      ! done when the model is not configured for homogeneous nucleation of
      ! sulfates.
      call sulfnucrate(carma,cstate, iz, igroup, nucbin, nucrate, rc)
      if (rc /= RC_OK) return
                    
      ! Do further calculations only if nucleation occurred
      if (nucrate .gt. 0._f) then

        rhompe(nucbin, ielem) = rhompe(nucbin, ielem) + nucrate
        
        ! Since homogeneous nucleation doesn't go through upgxfer or downgxfer, then
        ! then the effects of latent heat need to be accounted for here.
!        rlprod = rlprod + rhompe(nucbin, ielem) * rmass(nucbin,igroup) * rlh_nuc(ielem,ielem) / (CP * rhoa(iz))
      end if
    end if
  end do

  return
end
