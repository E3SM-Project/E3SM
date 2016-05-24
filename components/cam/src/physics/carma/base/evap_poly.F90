! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates particle source terms <evappe> due to
!! total evaporation into a polydisperse CN distribution by assuming
!! that the pdf of core mass is log-normal skewed by mass raised to
!! the -3/2 power (which guarantees average core mass from pdf is the
!! same as average core mass).
!!
!! Distinct evaporation of cores has not been treated.
!!
!! @author Andy Ackerman
!! @version Aug-2001
subroutine evap_poly(carma,cstate,iz,ibin,ig,iavg,ieto,igto,rc)

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
  integer                              :: ito  
  integer                              :: kount_s
  integer                              :: kount_l
  integer                              :: iecore
  integer                              :: ie2cn
  real(kind=f)                         :: prob(NBIN)
  real(kind=f)                         :: rn_norms
  real(kind=f)                         :: rn_norml
  real(kind=f)                         :: rm_norms
  real(kind=f)                         :: rm_norml
  real(kind=f)                         :: expon
  real(kind=f)                         :: rmassto
  real(kind=f)                         :: dmto
  real(kind=f)                         :: weightl
  real(kind=f)                         :: weights


  ! Treat total evaporation from a polydisperse core mass distribution:
  ! assume a log-normal CN size distribution and conserve number and mass as
  ! described by Turco (NASA Technical Paper 1362).
  !  
  ! Set automatic flag for total evaporation used in gasexchange()
  totevap(ibin,ig) = .true.
  
  ! Calculate number <rn_norms,rn_norml> and mass <rm_norms,rm_norml>
  ! normalization factors for cores smaller and larger than <rmass(m,igto)>.
  rn_norms = 0._f
  rn_norml = 0._f
  rm_norms = 0._f
  rm_norml = 0._f
  kount_s = 0
  kount_l = 0

  do ito = 1, NBIN

    rmassto = rmass(ito,igto)           
    dmto = dm(ito,igto)

    ! <prob> is probability that core mass is in CN bin <ito>.
    if( coreavg .gt. 0._f .and. coresig .gt. 0._f )then
       expon = -log( rmassto/coreavg )**2 / ( 2.*coresig )
       expon = max(-POWMAX, expon)
    else
      expon = 0._f
    endif
        
    prob(ito) = rmassto**(-1.5_f) * exp( expon )

    if( ito .lt. iavg )then
      rn_norms = rn_norms + prob(ito)*dmto
      rm_norms = rm_norms + prob(ito)*dmto*rmassto
      kount_s = kount_s + 1
    else
      rn_norml = rn_norml + prob(ito)*dmto
      rm_norml = rm_norml + prob(ito)*dmto*rmassto
      kount_l = kount_l + 1
    endif
  enddo
  
  ! Calculate mass weighting factors <weights,weightl> for small and
  ! large cores.
  if( kount_s .eq. 0 )then
    weightl = ONE
  elseif( kount_l .eq. 0 )then
    weightl = 0._f
  else
    rm_norms = rm_norms/rn_norms
    rm_norml = rm_norml/rn_norml
    weightl = (coreavg - rm_norms) / (rm_norml - rm_norms)
    if( weightl .gt. ALMOST_ONE )then
      weightl = ONE
    elseif( weightl .lt. ALMOST_ZERO )then
      weightl = 0._f
    endif
  endif
  
  weights = ONE - weightl

  ! Renormalize probability distribution function and evaluate the CN
  ! evaporation source term <evappe>.
  do ito = 1, NBIN

!    if( ito .le. iavg )then
    if( ito .lt. iavg )then      ! Kevin M
      prob(ito) = prob(ito)*weights/rn_norms
    else
      prob(ito) = prob(ito)*weightl/rn_norml
    endif
    
    ! First the CN number concentration element
    evappe(ito,ieto) = evappe(ito,ieto) + evdrop*prob(ito)*dm(ito,igto)
  
    ! Now the CN core elements
    do ic = 2, ncore(ig)
      iecore = icorelem(ic,ig)
      ie2cn  = ievp2elem(iecore)
      evappe(ito,ie2cn) = evappe(ito,ie2cn) + &
        rmass(ito,igto)*evcore(ic)*prob(ito)*dm(ito,igto)
    enddo
  enddo
 
  return
end
