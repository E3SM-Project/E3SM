! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates the dry deposition velocity, vd  [cm s^-1]
!! Method: Zhang et al., 2001
!! vd = vf(pver) + 1./ (rs + ra)
!! rs is the surface resistance, which is calculated in here
!! ra is the aerodynamic resistance, which is from parent dynamic model, like CAM
!! use carma_do_drydep flag optionally to decide if the CARMA or the parent model does the dry deposition
!! @author Tianyi Fan
!! @version Nov-2010
subroutine setupvdry(carma, cstate, lndfv, ocnfv, icefv, lndram, ocnram, iceram, lndfrac, ocnfrac, icefrac, rc)
  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod

  implicit none

  type(carma_type), intent(in)                :: carma    !! the carma object
  type(carmastate_type), intent(inout)        :: cstate   !! the carma state object
  real(kind=f), intent(in)                    :: lndfv    !! the surface friction velocity over land  [cm/s]
  real(kind=f), intent(in)                    :: ocnfv    !! the surface friction velocity over ocean  [cm/s]
  real(kind=f), intent(in)                    :: icefv    !! the surface friction velocity over ice  [cm/s]
  real(kind=f), intent(in)                    :: lndram   !! the aerodynamic resistance over land [s/cm]
  real(kind=f), intent(in)                    :: ocnram   !! the aerodynamic resistance over ocean [s/cm]
  real(kind=f), intent(in)                    :: iceram   !! the aerodynamic resistance over ice [s/cm]
  real(kind=f), intent(in)                    :: lndfrac  !! land fraction
  real(kind=f), intent(in)                    :: ocnfrac  !! ocn fraction
  real(kind=f), intent(in)                    :: icefrac  !! ice fraction
  integer, intent(inout)                      :: rc       !! return code, negative indicates failure

  ! Local declarations
  integer         :: ielem, igroup, ibin, icnst, k
  real(kind=f)    :: vd_lnd, vd_ocn, vd_ice  ! the deposition velocity of land,ocean and sea ice
  real(kind=f)    :: rs                      ! surface resistance  [s/m]
  real(kind=f)    :: vfall(NBIN, NGROUP)     ! fall velocity [m/s]
  integer         :: cnsttype                ! if constituent is prognostic
  integer         :: maxbin                  ! last prognostic bin
  integer         :: ibot, ibotp1            ! index of bottom layer


  if (do_drydep) then

    if (igridv .eq. I_CART) then
      ibot = 1
      ibotp1 = 1
      vfall(:,:) = vf(ibotp1, :, :)  ![cm/s]
    else
      ibot = NZ
      ibotp1 = NZP1
      vfall(:,:) = -vf(ibotp1, :, :) * zmetl(ibotp1)  ! [z_unit/s] -> [cm/s]
    end if

    do ielem = 1, NELEM
      igroup = igelem(ielem)

      if (grp_do_drydep(igroup)) then
        do ibin = 1, NBIN
          vd_lnd = 0._f
          vd_ocn = 0._f
          vd_ice = 0._f

          ! land
          if (lndfrac > 0._f) then
            call calcrs(carma, cstate, lndfv, t(ibot), r_wet(ibot, ibin, igroup), &
                 bpm(ibot, ibin, igroup), vfall(ibin,igroup), rs, 1, rc)
            vd_lnd = vfall(ibin, igroup) + 1._f / (lndram + rs)
          end if

          ! ocean
          if (ocnfrac > 0._f) then
            call calcrs(carma, cstate, ocnfv, t(ibot), r_wet(ibot, ibin, igroup), &
                 bpm(ibot, ibin, igroup), vfall(ibin,igroup), rs, 2, rc)
            vd_ocn = vfall(ibin, igroup) + 1._f / (ocnram + rs)
          end if

          ! sea ice
          if (icefrac > 0._f) then
            call calcrs(carma, cstate, icefv, t(ibot), r_wet(ibot, ibin, igroup), &
                 bpm(ibot, ibin, igroup), vfall(ibin,igroup), rs, 3, rc)
            vd_ice = vfall(ibin, igroup) + 1._f / (iceram + rs)
          end if

          vd(ibin, igroup) = (lndfrac * vd_lnd + ocnfrac * vd_ocn + icefrac * vd_ice)   ![cm/s]
        end do   ! ibin
      else
        vd(:, igroup) = vfall(:, igroup)   ! [cm/s]
      end if  ! if grp_do_drydep
    end do  ! ielem

    ! change scale for non-catesian vertical coordinate
    ! Scale cartesian fallspeeds to the appropriate vertical coordinate system.
    ! Non--cartesion coordinates are assumed to be positive downward, but
    ! vertical velocities in this model are always assumed to be positive upward.
    if( igridv /= I_CART )then
      vd(:,:) = -vd(:,:) / zmetl(NZP1)
    end if
  end if

  return
end
