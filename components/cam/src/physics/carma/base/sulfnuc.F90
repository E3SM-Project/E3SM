! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  Calculates particle production rates due to nucleation <rhompe>:
!!  binary homogeneous nucleation of sulfuric acid and water only
!!  Numerical method follows Zhao & Turco, JAS, V.26, No.5, 1995.
!!
!!  @author Mike Mills, Chuck Bardeen
!!  @version Jun-2013
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
  integer           :: ibin       ! bin index
  integer           :: igas       ! gas index
  integer           :: iepart     ! concentration element index
  integer           :: nucbin     ! bin in which nucleation takes place
  integer           :: ignucto    ! index of target nucleation group
  integer           :: ienucto    ! index of target nucleation element
  integer           :: inuc
  real(kind=f)      :: nucrate    ! nucleation rate (#/x/y/z/s)
  real(kind=f)      :: h2o        ! H2O concentrations in molec/cm3 
  real(kind=f)      :: h2so4      ! H2SO4 concentrations in molec/cm3
  real(kind=f)      :: beta1
  real(kind=f)      :: beta2
  real(kind=f)      :: ftry
  real(kind=f)      :: rstar      ! critical radius (cm)

  ! Cycle through each group, only proceed if BHN
  rstar = -1._f
  
  do igroup = 1 , NGROUP
    
    igas = inucgas(igroup)                ! condensing gas
    
    if (igas .ne. 0) then

      iepart = ienconc(igroup)              ! particle number density element

      if (inucproc(iepart,iepart) .eq. I_HOMNUC) then
    
        ! This is where all of the pre calculation needs to go, so that it isn't
        ! done when the model is not configured for homogeneous nucleation of
        ! sulfates.
        call sulfnucrate(carma, cstate, iz, igroup, h2so4, h2o, beta1, beta2, ftry, rstar, nucbin, nucrate, rc)
        if (rc /= RC_OK) return
        
        ! Do further calculations only if nucleation occurred
        if (nucrate .gt. 0._f) then

          rhompe(nucbin, iepart) = rhompe(nucbin, iepart) + nucrate
        
          ! Since homogeneous nucleation doesn't go through upgxfer or downgxfer, then
          ! then the effects of latent heat need to be accounted for here.
  !        rlprod = rlprod + rhompe(nucbin, ielem) * rmass(nucbin,igroup) * rlh_nuc(ielem,ielem) / (CP * rhoa(iz))
        end if
      end if
    end if
  end do

  ! Cycle through each group, only proceed if heterogeneous nucleation
  !
  ! NOTE: Only do heterogeneous nucleation if an rstar was determined by homogeneous
  ! nucleation.
  if (rstar > 0._f) then
    do igroup = 1 , NGROUP
    
      igas = inucgas(igroup)                ! condensing gas
    
      if (igas .ne. 0) then

        iepart = ienconc(igroup)              ! particle number density element

        ! Calculate heterogeneous nucleation loss rates.  Do not allow nucleation into
        ! an evaporating bin.
        !
        ! NOTE: Heterogeneous nucleation assumes that homogeneous nucleation was called
        ! first to determine the critical cluster size.
        !
        ! <ienucto> is index of target nucleation element;
        ! <ignucto> is index of target nucleation group.
        do inuc = 1, nnuc2elem(iepart)

          ienucto = inuc2elem(inuc,iepart)
        
          if (ienucto .ne. 0) then
            ignucto = igelem(ienucto)
          else
            ignucto = 0
          endif
        
          if (inucproc(iepart,ienucto) .eq. I_HETNUCSULF) then
    
            do ibin = NBIN, 1, -1

              ! Bypass calculation if few particles are present
              if (pconmax(iz,igroup) .gt. FEW_PC) then

                ! This is where all of the pre calculation needs to go, so that it isn't
                ! done when the model is not configured for homogeneous nucleation of
                ! sulfates.
                call sulfhetnucrate(carma, cstate, iz, igroup, ibin, h2so4, h2o, beta1, beta2, ftry, rstar, nucrate, rc)
                if (rc /= RC_OK) return
                    
                rnuclg(ibin, igroup, ignucto) = rnuclg(ibin, igroup, ignucto) + nucrate
              end if
            end do
          end if
        end do
      end if

    end do
  end if
  
  return
end
