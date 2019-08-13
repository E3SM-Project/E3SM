! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluates particle loss rates due to nucleation <rnuclg>:
!! droplet freezing only.
!! 
!! The loss rates for all particle elements in a particle group are equal.
!!
!! @author Eric Jensen, Chuck Bardeen
!! @version Jan-2000, Nov-2009
subroutine freezdropl(carma, cstate, iz, rc)

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
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  !  Local declarations
  integer                              :: igroup  !! group index
  integer                              :: ibin    !! bin index
  integer                              :: iepart  !! element for condensing group index
  integer                              :: inuc    !! nucleating element index
  integer                              :: ienucto !! index of target nucleation element
  integer                              :: ignucto !! index of target nucleation group


  ! Loop over particle groups.
  do igroup = 1,NGROUP
  
    iepart = ienconc( igroup )            ! particle number density element
  
    ! Calculate nucleation loss rates.
    do inuc = 1,nnuc2elem(iepart)
    
      ienucto = inuc2elem(inuc,iepart)
      
      if( ienucto .ne. 0 )then
        ignucto = igelem( ienucto )
  
        ! Only compute nucleation rate for droplet freezing
        if( inucproc(iepart,ienucto) .eq. I_DROPFREEZE ) then
    
          ! Loop over particle bins.  
          do ibin = 1,NBIN
      
            ! Bypass calculation if few particles are present 
            if( pc(iz,ibin,iepart) .gt. FEW_PC )then
      
              ! Temporary simple kludge: Set <rnuclg> to 1.e2 if T < -40C
              if( t(iz) .lt. T0-40._f ) then
                rnuclg(ibin,igroup,ignucto) = 1.e2_f
              endif
      
            endif     ! pc(source particles) .gt. FEW_PC
          enddo      ! ibin = 1,NBIN
        endif       ! inucproc(iepart,ienucto) .eq. I_DROPFREEZE
      endif
    enddo        ! inuc = 1,nnuc2elem(iepart)
  enddo         ! igroup = 1,NGROUP

  ! Return to caller with particle loss rates due to nucleation evaluated.
  return
end
