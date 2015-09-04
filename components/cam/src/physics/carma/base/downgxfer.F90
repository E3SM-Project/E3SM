! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates particle source terms <rnucpe> due to particle
!! element transfer processes for which the source element number is larger
!! than the target element number.
!!
!! @author Andy Ackerman
!! @version Dec-1995
subroutine downgxfer(carma, cstate, iz, rc)

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

  ! Local declarations
  integer                              :: igroup  ! group index
  integer                              :: iepart
  integer                              :: ibin    !! bin index
  integer                              :: ielem   !! element index
  integer                              :: i
  integer                              :: jefrom
  integer                              :: iefrom
  integer                              :: igfrom
  integer                              :: ipow_from
  integer                              :: ipow_to
  integer                              :: ipow
  integer                              :: jfrom
  integer                              :: ifrom
  integer                              :: ic
  integer                              :: iecore
  real(kind=f)                         :: elemass
  real(kind=f)                         :: totmass
  real(kind=f)                         :: rmasscore
  real(kind=f)                         :: fracmass
  real(kind=f)                         :: rnucprod
  

  ! Calculate nucleation source terms for which the source element
  ! number is greater than the target element number

  ! Set nucleation production rates to zero to avoid double-application
  ! of rates calculated in upgxfer.f
  rnucpe(:,:) = 0._f

  ! Loop over particle elements and bins
  do ielem = 1, NELEM
    do ibin = 1, NBIN

      ! Define group & particle # concentration indices for current element
      igroup = igelem(ielem)      ! target particle group
      iepart = ienconc(igroup)    ! target particle number concentration element
      
      ! First calculate production terms due to nucleation <rnucpe>.

      ! Loop over elements that nucleate to element <ielem>.
      do jefrom = 1,nnucelem(ielem)

        iefrom = inucelem(jefrom,ielem)    ! source particle element

        ! Only calculate production rates here if <ielem> is less than
        ! <iefrom>.  Otherwise, production is calculated in upgxfer.f
        if( ielem .lt. iefrom ) then
          igfrom = igelem(iefrom)            ! source particle group

          ! <ipow> is the power to which the source particle mass must be taken
          ! to match the type of the target element.  This ugliness could be
          ! handled much more slickly in setupnuc()
          if( itype(iefrom) .eq. I_INVOLATILE .or. itype(iefrom) .eq. I_VOLATILE )then
            ipow_from = 0
          elseif ( itype(iefrom) .eq. I_COREMASS .or. itype(iefrom) .eq. I_VOLCORE )then
            ipow_from = 1
          else
            ipow_from = 2
          endif

          if( itype(ielem) .eq. I_INVOLATILE .or. itype(ielem) .eq. I_VOLATILE )then
            ipow_to = 0
          elseif ( itype(ielem) .eq. I_COREMASS .or. itype(ielem) .eq. I_VOLCORE )then
            ipow_to = 1
          else
            ipow_to = 2
          endif

          ipow = ipow_to - ipow_from

          ! Loop over bins that nucleate to bin <ibin>.
          do jfrom = 1,nnucbin(igfrom,ibin,igroup)

            ifrom = inucbin(jfrom,igfrom,ibin,igroup)    ! bin of source

            ! Bypass calculation if few source particles are present
            if( pconmax(iz,igfrom) .gt. FEW_PC )then
              if( rnuclg(ifrom,igfrom,igroup) .gt. 0._f )then
           
                ! First calculate mass associated with the source element <elemass>
                ! (this is <rmass> for all source elements except particle number
                ! concentration in a multicomponent particle group).
                if( ncore(igfrom) .eq. 0 .or. itype(iefrom) .gt. I_VOLATILE )then
                  elemass = rmass(ifrom,igfrom)
                else
                  totmass  = pc(iz,ifrom,iefrom) * rmass(ifrom,igfrom)
                  rmasscore = pc(iz,ifrom,icorelem(1,igfrom))
                  
                  do ic = 2,ncore(igfrom)
                    iecore = icorelem(ic,igfrom)
                    rmasscore = rmasscore + pc(iz,ifrom,iecore)
                  enddo
                  
                  fracmass = 1._f - rmasscore/totmass
                  elemass  = fracmass * rmass(ifrom,igfrom)
                endif

                rnucprod = rnuclg(ifrom,igfrom,igroup) * &
                  pc(iz,ifrom,iefrom) * elemass**ipow

                rnucpe(ibin,ielem) = rnucpe(ibin,ielem) + rnucprod

                ! Calculate latent heat associated with nucleation to <ibin,ielem>
                ! from <ifrom,iefrom>
!                rlprod = rlprod + rnucprod * rlh_nuc(iefrom,ielem) / &
!                  (CP * rhoa(iz)) * elemass

              endif  ! (rnuclg > 0.)
            endif   ! (pconmax > FEW_PC)
          enddo    ! (jfrom = 1,nnucbin)
        endif     ! (ielem < iefrom)
      enddo      ! (jefrom = 1,nnucelem)
    enddo       ! (ibin = 1, NBIN)
  enddo       ! (ielem = 1, NELEM)

  ! Return to caller with down-grid production terms evaluated.
  return
end
