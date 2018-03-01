! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluates particle loss rates due to nucleation <rnuclg>:
!! droplet activation only.
!!  
!! The loss rates for all particle elements in a particle group are equal.
!!
!! To avoid nucleation into an evaporating bin, this subroutine must
!! be called after growp, which evaluates evaporation loss rates <evaplg>.
!!
!! @author Andy Ackerman
!! @version Dec-1995
subroutine actdropl(carma, cstate, iz, rc)

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
  integer                              :: igas    !! gas index
  integer                              :: igroup  !! group index
  integer                              :: ibin    !! bin index
  integer                              :: iepart  !! element for condensing group index
  integer                              :: inuc    !! nucleating element index
  integer                              :: ienucto !! index of target nucleation element
  integer                              :: ignucto !! index of target nucleation group
  integer                              :: inucto  !! index of target nucleation bin
  logical                              :: evapfrom_nucto !! .true. when target droplets are evaporating
  

  ! This calculation is only necessary for temperatures greater
  ! than -40C.
  if( t(iz) .ge. (T0 - 40._f) ) then

    ! Loop over particle groups.
    do igroup = 1,NGROUP
  
      ! Bypass calculation if few particles are present
      if( pconmax(iz,igroup) .gt. FEW_PC )then
  
        igas = inucgas(igroup)                ! condensing gas
        iepart = ienconc( igroup )            ! particle number density element
  
        if( igas .ne. 0 )then
  
          ! Calculate nucleation loss rates.  Do not allow nucleation into
          ! an evaporating bin.
          do inuc = 1,nnuc2elem(iepart)
          
            ienucto = inuc2elem(inuc,iepart)
            if( ienucto .ne. 0 )then
              ignucto = igelem( ienucto )
            else
              ignucto = 0
            endif
      
            ! Only compute nucleation rate for droplet activation
            if( inucproc(iepart,ienucto) .eq. I_DROPACT ) then
      
              ! Loop over particle bins.  Loop from largest to smallest for 
              ! evaluation of index of smallest bin nucleated during time step <inucstep>.
              do ibin = NBIN, 1, -1
      
                if( ignucto .ne. 0 )then
                  inucto = inuc2bin(ibin,igroup,ignucto)
                else
                  inucto = 0
                endif
      
                ! Set <evapfrom_nucto> to .true. when target droplets are evaporating
                if( inucto .ne. 0 )then
                 evapfrom_nucto = evaplg(inucto,ignucto) .gt. 0._f
                else
                 evapfrom_nucto = .false.
                endif
                  
                if( (supsatl(iz,igas) .gt. scrit(iz,ibin,igroup)) .and. &
                    (.not. evapfrom_nucto) .and. &
                    (pc(iz,ibin,iepart) .gt. SMALL_PC) )then
      
                  rnuclg(ibin,igroup,ignucto) = 1.e3_f
                endif
              enddo   ! ibin = 1,NBIN
            endif    ! inucproc(iepart,ienucto) .eq. I_DROPACT
          enddo     ! inuc = 1,nnuc2elem(iepart)
        endif      ! (igas = inucgas(igroup)) .ne. 0 
      endif       ! pconmax(iz,igroup) .gt. FEW_PC
    enddo        ! igroup = 1,NGROUP
  endif         ! t(iz) .ge. T0-40.

  ! Return to caller with particle loss rates due to nucleation evaluated.
  return
end
