! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluates particle loss rates due to nucleation <rnuclg>:
!! Ice crystal melting only.
!! 
!! The loss rates for all particle elements in a particle group are equal.
!!
!! @author Eric Jensen, Chuck Bardeen
!! @version Jan-2000, Nov-2009
subroutine melticel(carma, cstate, iz, rc)

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
  integer                              :: ienucto !! index of target nucleation element
  integer                              :: ignucto !! index of target nucleation group
  integer                              :: inuc    !! nucleating element index


  ! Loop over particle groups.
  do igroup = 1,NGROUP

    iepart = ienconc( igroup )            ! particle number density element

    ! Calculate nucleation loss rates.
    do inuc = 1,nnuc2elem(iepart)

      ienucto = inuc2elem(inuc,iepart)
      
      if( ienucto .ne. 0 )then
        ignucto = igelem( ienucto )

        ! Only compute nucleation rate for ice crystal melting
        if( inucproc(iepart,ienucto) .eq. I_ICEMELT ) then
  
          ! Loop over particle bins.  Loop from largest to smallest for 
          ! evaluation of index of smallest bin nucleated during time step <inucstep>.
          do ibin = NBIN,1,-1
  
            ! Bypass calculation if few particles are present 
            if( pconmax(iz,igroup) .gt. FEW_PC )then
  
              ! Temporary simple kludge: Set <rnuclg> to 1.e2 if T > 0C
              if( t(iz) .gt. T0 ) then
                rnuclg(ibin,igroup,ignucto) = 1.e2_f
              endif
            endif   ! pconmax(ixyz,igroup) .gt. FEW_PC
          enddo      ! ibin = 1,NBIN
        endif       ! inucproc(iepart,ienucto) .eq. I_DROPFREEZE
      endif
    enddo       ! inuc = 1,nnuc2elem(iepart)
  enddo         ! igroup = 1,NGROUP

  ! Return to caller with particle loss rates due to nucleation evaluated.
  return
end
