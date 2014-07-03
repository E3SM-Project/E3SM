! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine calculates new average particle densities.
!!
!!  The particle mass density can change at each time-step due to
!!  changes in the core mass fraction.
!!
!!  For particles that are hydrophilic and whose particle size changes based
!!  upon the relative humidity, and wet radius and density are also calculated.
!!  For particles that do not swell, the wet and dry radius and densities are
!!  the same.
!!
!! @author   Chuck Bardeen Eric Jensen
!! @ version May-2009; Oct-1995
!!
!! @see wetr
subroutine rhopart(carma, cstate, rc)

  ! types
  use carma_precision_mod
  use carma_enums_mod
  use carma_constants_mod
  use carma_types_mod
  use carmastate_mod
  use carma_mod
  use sulfate_utils
  use wetr

  implicit none

  type(carma_type), intent(in)         :: carma   !! the carma object
  type(carmastate_type), intent(inout) :: cstate  !! the carma state object
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local declarations
  integer         :: iz             !! z index
  integer         :: igroup         !! group index
  integer         :: ibin           !! bin index
  integer         :: iepart         !! element in group containing the particle concentration
  integer         :: jcore
  integer         :: iecore
  real(kind=f)    :: vcore(NBIN)
  real(kind=f)    :: mcore(NBIN)
  real(kind=f)    :: r_ratio
  real(kind=f)    :: gc_cgs

  1 format(/,'rhopart::WARNING - core mass > total mass, truncating : iz=',i4,',igroup=',&
              i4,',ibin=',i4,',total mass=',e10.3,',core mass=',e10.3,',using rhop=',f9.4)

  ! Calculate average particle mass density for each group
  do igroup = 1,NGROUP

    ! Define particle # concentration element index for current group
    iepart = ienconc(igroup)     ! element of particle number concentration 
    
    do iz = 1, NZ
    
      ! If there are no cores, than the density of the particle is just the density
      ! of the element.
      if (ncore(igroup) < 1) then
        rhop(iz,:,igroup) = rhoelem(:,iepart)
      
      ! Otherwise, the density changes depending on the amount of core and volatile
      ! components.
      else

        ! Calculate volume of cores and the mass of shell material
        ! <vcore> is the volume of core material and <rmshell> is the
        ! mass of shell material.
        vcore(:) = 0._f
        mcore(:) = 0._f
  
        do jcore = 1,ncore(igroup)
          iecore = icorelem(jcore,igroup)    ! core element
  
          mcore(:) = mcore(:) + pc(iz,:,iecore)
          vcore(:) = vcore(:) + pc(iz,:,iecore) / rhoelem(:,iecore)
        enddo
      
        ! Calculate average density
        do ibin = 1,NBIN
                  
          ! If there is no core, the the density is that of the volatile element.
          if (mcore(ibin) == 0._f) then
            rhop(iz,ibin,igroup) = rhoelem(ibin,iepart)
          else
          
            ! Since core mass and particle number (i.e. total mass) are advected separately,
            ! numerical diffusion during advection can cause problems where the core mass
            ! becomes greater than the total mass. To prevent adevction errors from making the
            ! group inconsistent, we will truncate core mass if it is larger than the total
            ! mass.
            if (mcore(ibin) > (rmass(ibin,igroup) * pc(iz,ibin,iepart))) then
            
              ! Calculate the density.
              rhop(iz,ibin,igroup) = mcore(ibin) / vcore(ibin)
              
              ! NOTE: This error happens a lot, so this error message is commented out
              ! by default.
!              if (do_print) write(LUNOPRT,1) iz, igroup, ibin, pc(iz,ibin,iepart)*rmass(ibin,igroup), &
!                mcore(ibin), rhop(iz,ibin,igroup)
!              rc = RC_WARNING

              ! Repair total mass.
              pc(iz,ibin,iepart) = mcore(ibin) / rmass(ibin,igroup)
            else 
              rhop(iz,ibin,igroup) = (rmass(ibin,igroup) * pc(iz,ibin,iepart)) / &
              ((pc(iz,ibin,iepart)*rmass(ibin,igroup) - mcore(ibin))/rhoelem(ibin,iepart) + vcore(ibin))
            end if
          end if
        enddo
      endif
    
      ! If these particles are hygroscopic and grow in response to the relative
      ! humidity, then caclulate a wet radius and wet density. Otherwise the wet
      ! and dry radius are the same.
    
      ! Determine the weight percent of sulfate, and store it for later use.
      if (irhswell(igroup) == I_WTPCT_H2SO4) then
        gc_cgs     = gc(iz, igash2so4) / (xmet(iz) * ymet(iz) * zmet(iz))
        wtpct(iz)  = wtpct_tabaz(carma, t(iz), gc_cgs, pvapl(iz, igash2o), rc)
        if (rc < 0) return
      end if
          
      ! Loop over particle size bins.
      do ibin = 1,NBIN
      
        ! If humidty affects the particle, then determine the equilbirium
        ! radius and density based upon the relative humidity.
        if (irhswell(igroup) == I_WTPCT_H2SO4) then
        
          ! rlow
          call getwetr(carma, igroup, relhum(iz), rlow(ibin,igroup), rlow_wet(iz,ibin,igroup), &
            rhop(iz,ibin,igroup), rhop_wet(iz,ibin,igroup), rc, wgtpct=wtpct(iz), temp=t(iz))
          if (rc < 0) return

          ! rup
          call getwetr(carma, igroup, relhum(iz), rup(ibin,igroup), rup_wet(iz,ibin,igroup), &
            rhop(iz,ibin,igroup), rhop_wet(iz,ibin,igroup), rc, wgtpct=wtpct(iz), temp=t(iz))
          if (rc < 0) return

          ! r
          call getwetr(carma, igroup, relhum(iz), r(ibin,igroup), r_wet(iz,ibin,igroup), &
            rhop(iz,ibin,igroup), rhop_wet(iz,ibin,igroup), rc, wgtpct=wtpct(iz), temp=t(iz))
          if (rc < 0) return

        else
          ! rlow
          call getwetr(carma, igroup, relhum(iz), rlow(ibin,igroup), rlow_wet(iz,ibin,igroup), &
            rhop(iz,ibin,igroup), rhop_wet(iz,ibin,igroup), rc)
          if (rc < 0) return

          ! rup
          call getwetr(carma, igroup, relhum(iz), rup(ibin,igroup), rup_wet(iz,ibin,igroup), &
            rhop(iz,ibin,igroup), rhop_wet(iz,ibin,igroup), rc)
          if (rc < 0) return

          ! r
          call getwetr(carma, igroup, relhum(iz), r(ibin,igroup), r_wet(iz,ibin,igroup), &
            rhop(iz,ibin,igroup), rhop_wet(iz,ibin,igroup), rc)
          if (rc < 0) return
        end if
      end do
    end do
  enddo

  !  Return to caller with new particle number densities.
  return
end
