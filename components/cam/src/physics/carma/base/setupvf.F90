! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine calculates fall velocities for particles. Since there are
!!  several different approaches, this routine dispatches the call to the
!!  proper subordinate routine based upon the setup routine defined in the
!!  particle group.
!!
!!
!! @author Andy Ackerman
!! @version Mar-2010
subroutine setupvf(carma, cstate, rc)

  !     types
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
  integer                 :: igroup, i, j, k, k1, k2, ibin, iz, nzm1

  !  Define formats
  2 format(/,'Fall velocities and Reynolds number (prior to interpolation)')
  3 format(/,'Particle group ',i3,' using algorithm ',i3,/, &
         ' bin   lev  p [dyne/cm2]  T [K]       r [cm]    wet r [cm]       bpm', &
      '   vf [cm/s]        re'/)
  4 format(i3,4x,i3,7(1pe11.3,4x)) 
  
  !  Loop over all groups.
  do igroup = 1, NGROUP
      
    ! There are different implementations of the fall velocity calculation. Some of
    ! these routines may be more appropriate for certain types of partciles.
    select case(ifallrtn(igroup))

      case (I_FALLRTN_STD)
        call setupvf_std(carma, cstate, igroup, rc)
      
      case(I_FALLRTN_STD_SHAPE)
        call setupvf_std_shape(carma, cstate, igroup, rc)
      
      case(I_FALLRTN_HEYMSFIELD2010)
        call setupvf_heymsfield2010(carma, cstate, igroup, rc)
      
      case default
        if (do_print) write(LUNOPRT,*) "setupvf:: ERROR - Unknown fall velocity routine  (", ifallrtn(igroup), &
          ") for group (", igroup, ")."
        rc = -1
        return
    end select
  enddo
    
  ! Constant value if <ifall> = 0
  if (ifall .eq. 0) then
    vf(:,:,:) = vf_const
  end if

  ! Print out fall velocities and reynolds' numbers.
#ifdef DEBUG
  if (do_print_init) then
   
    write(LUNOPRT,2)

    do j = 1, NGROUP
        
      write(LUNOPRT,3) j, ifallrtn(j)
      
      do i = 1,NBIN
    
        do k = NZ, 1, -1
          write(LUNOPRT,4) i,k,p(k),t(k),r(i,j),r_wet(k,i,j),bpm(k,i,j),vf(k,i,j),re(k,i,j)
        end do
      enddo
    enddo
  
    write(LUNOPRT,*) ""
  end if
#endif

  ! Interpolate <vf> from layer mid-pts to layer boundaries.
  ! <vf(k)> is the fall velocity at the lower edge of the layer
  nzm1 = max(1, NZ-1)

  ! Set upper boundary before averaging
  vf(NZP1,:,:) = vf(NZ,:,:)

  if (NZ .gt. 1) then
    vf(NZ,:,:) = sqrt(vf(nzm1,:,:) * vf(NZ,:,:))

    if (NZ .gt. 2) then
      do iz = NZ-1, 2, -1
        vf(iz,:,:) = sqrt(vf(iz-1,:,:) * vf(iz,:,:))
      enddo
    endif
  endif
  
  ! Scale cartesian fallspeeds to the appropriate vertical coordinate system.
  ! Non--cartesion coordinates are assumed to be positive downward, but
  ! vertical velocities in this model are always assumed to be positive upward. 
  if( igridv .ne. I_CART )then
    do igroup=1,NGROUP
      do ibin=1,NBIN
        vf(:,ibin,igroup) = -vf(:,ibin,igroup) / zmetl(:)
      enddo
    enddo
  endif

  ! Return to caller with fall velocities evaluated.
  return
end
