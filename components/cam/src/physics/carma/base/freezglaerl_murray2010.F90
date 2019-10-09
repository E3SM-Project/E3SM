! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine evaluates particle loss rates due to nucleation <rnuclg>:
!! heterogeneous nucleation of glassy aerosols only,.
!!
!! The parameterization of glass aerosols is described in Murray et al.
!! [Nature Geosciences, 2010], and is based upon measurements of the nucleation of
!! citric acid aerosols at cold temperatures.
!!
!! NOTE: This implementation assumes that the aerosol being nucleated is the total
!! aerosol population and not just the fraction of aerosols that are glassy. To
!! account for homogenous freezing of the aerosol population, the routine freezaerl
!! also needs to be called and the overall nucleation rate is the sum of
!! the rates for homogeneous freezing and for heterogenous nucleation.
!!
!! The parameter fglass is the fraction of the total aerosol population that will be
!! in a glassy state for T <= 212K.
!! 
!! The loss rates for all particle elements in a particle group are equal.
!!
!! @author Chuck Bardeen, Eric Jensen
!! @version Apr-2010
subroutine freezglaerl_murray2010(carma, cstate, iz, rc)

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
	! Define parameters needed for freezing nucleation calculations.
  real(kind=f), parameter              :: kice1   = 7.7211e-5_f  ! Fit constant from Murray et al.
  real(kind=f), parameter              :: kice2   = 9.2688e-3_f  ! Fit constant from Murray et al.
  real(kind=f), parameter              :: ssmin   = 0.21_f       ! Minimum supersaturation for nucleation
  real(kind=f), parameter              :: ssmax   = 0.7_f        ! Maximum supersaturation for nucleation
  real(kind=f), parameter              :: tglass  = 212._f       ! Maximum temperature for glassy state
  real(kind=f), parameter              :: fglass  = 0.5_f        ! Fraction of aerosols that can become glassy
	
  integer                              :: igas    ! gas index
  integer                              :: igroup  ! group index
  integer                              :: ibin    ! bin index
  integer                              :: iepart  ! element for condensing group index
  integer                              :: inuc    ! nucleating element index
  integer                              :: ienucto ! index of target nucleation element
  integer                              :: ignucto ! index of target nucleation group
  integer                              :: inucto  ! index of target nucleation bin
  real(kind=f)                         :: dfice   ! difference in fraction of aerosol nucleated
  real(kind=f)                         :: ssi, ssiold
  
  ! Assume success.
  rc = RC_OK

  ! Loop over particle groups.
  do igroup = 1,NGROUP

    igas = inucgas(igroup)                ! condensing gas
    iepart = ienconc( igroup )            ! particle number density element

    if (igas /= 0) then

      ! Calculate nucleation loss rates.  Do not allow nucleation into
      ! an evaporating bin.
      do inuc = 1, nnuc2elem(iepart)

        ienucto = inuc2elem(inuc,iepart)
        if (ienucto /= 0) then
          ignucto = igelem(ienucto)
  
          ! Only compute nucleation rate for glassy aerosol freezing.
              if ((iand(inucproc(iepart,ienucto), I_AF_MURRAY_2010) /= 0)) then
          
            ! Is it cold enough for aerosols to be in a glassy state.
            if (t(iz) <= tglass) then
          
              !  Loop over particle bins.  Loop from largest to smallest for 
              !  evaluation of index of smallest bin nucleated during time step <inucstep>.
              do ibin = NBIN, 1, -1
      
                ! Bypass calculation if few particles are present or if it isn't cold enough
                ! for the aerosols to be present in a glassy state.
                if (pconmax(iz,igroup) > FEW_PC) then
                
                  ! Murray et al. [2010] doesn't really give a nucleation rate. Instead it gives
                  ! a fraction of glassy aerosol particles that have been nucleated as a function
                  ! of ice supersaturation.
                  !
                  ! Since CARMA really wants to work with rates, use the difference in relative
                  ! humidity and the length of the timestep to come up with an approximation to
                  ! a nucleation rate.
                  
                  ! The supersaturation must be greater than .21 for heterogeneous nucleation to
                  ! commence. The fraction of glassy aerosol nucleated is:
                  !
                  !   fice = 7.7211e-5 * RHi(%) - 9.2688e-3   for 121 % < RHi < 170 %
                  !
                  ! To get a pseudo production rate, use
                  !
                  !   rnuclg = (fice(RHi) - fice(RHi_old)) / dtime
                  !
                  ssi    = supsati(iz,igas)
                  ssiold = supsatiold(iz,igas)
                  
                  if ((ssi >= ssmin) .and. (ssi > ssiold)) then
                    dfice = kice1 * (1._f + min(ssmax, ssi)) * 100._f - kice2
                    
                    if (ssiold >= ssmin) then
                      dfice = dfice - (kice1 * (1._f + min(ssmax, ssiold)) * 100._f - kice2)
                    endif
                    
                    ! Add the rate of heterogenous freezing to the rate of homogeneous
                    rnuclg(ibin,igroup,ignucto) = rnuclg(ibin,igroup,ignucto) + fglass * dfice / dtime
                  endif    
                endif
              enddo
            endif
          endif
        endif
      enddo 
    endif
  enddo

  ! Return to caller with particle loss rates due to nucleation evaluated.
  return
end
