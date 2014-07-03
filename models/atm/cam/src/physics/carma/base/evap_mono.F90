! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates particle source terms <evappe> due to total
!! evaporation from bin <ibin> group <ig> into a monodisperse
!! distribution.
!!
!! Distinct evaporation of cores has not been treated.
!!
!! @author Andy Ackerman
!! @version Aug-2001
subroutine evap_mono(carma,cstate,iz,ibin,ig,iavg,ieto,igto,rc)

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
  integer, intent(in)                  :: ibin    !! bin index
  integer, intent(in)                  :: ig      !! group index
  integer, intent(in)                  :: iavg
  integer, intent(in)                  :: ieto
  integer, intent(in)                  :: igto
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer                              :: ic  
  integer                              :: iecore
  integer                              :: ie2cn
  integer                              :: jbin
  logical                              :: conserve_mass
  real(kind=f)                         :: factor
  real(kind=f)                         :: fracmass


  ! Define option to conserve mass or number when a choice must be made
  ! during monodisperse total evaporation beyond CN grid -- should be done in setupaer()
  conserve_mass = .true.

  ! Set automatic flag for total evaporation used in gasexchange()
  totevap(ibin,ig) = .true.

  ! Possibly put all of core mass into largest, smallest, or
  ! smallest nucelated CN bin 
  if( too_big .or. too_small .or. nuc_small )then

    if( too_big )then
      jbin = NBIN
    elseif( too_small )then
      jbin = 1
    else
      jbin = 1
    endif

    if( conserve_mass )then
      factor = coreavg/rmass(jbin,igto)
    else
      factor = ONE
    endif

    ! First the CN number concentration element
    evappe(jbin,ieto) = evappe(jbin,ieto) + factor*evdrop

    ! Now the CN cores
    do ic = 2, ncore(ig)
      iecore = icorelem(ic,ig)
      ie2cn  = ievp2elem(iecore)
      evappe(jbin,ie2cn) = evappe(jbin,ie2cn) + &
        factor*evcore(ic)*rmass(jbin,igto)
    enddo
  else

    ! Partition core mass between two CN bins, conserving total core mass
    ! and number.  The number will be subdivided into bins <iavg> and <iavg>-1.
    if( iavg .le. 1 .or. iavg .gt. NBIN )then
      if (do_print) write(LUNOPRT, *) "evap_mono: bad iavg = , ", iavg
      rc = RC_ERROR
      return
    endif

    fracmass = ( rmass(iavg,igto) - coreavg ) / diffmass(iavg,igto,iavg-1,igto)
!    fracmass = max( 0._f, min( ONE, fracmass ) )
  
    ! First the CN number concentration element
    evappe(iavg-1,ieto) = evappe(iavg-1,ieto) + evdrop*fracmass
    evappe(iavg,ieto) = evappe(iavg,ieto) + evdrop*( ONE - fracmass )
  
    ! Now the cores
    do ic = 2, ncore(ig)
      iecore = icorelem(ic,ig)
      ie2cn  = ievp2elem(iecore)
      evappe(iavg-1,ie2cn) = evappe(iavg-1,ie2cn) + &
          rmass(iavg-1,igto)*evcore(ic)*fracmass
      evappe(iavg,ie2cn) = evappe(iavg,ie2cn) + &
          rmass(iavg,igto)*evcore(ic)*( ONE - fracmass )
    enddo
  endif

  return
end
