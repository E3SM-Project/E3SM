! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!! This routine calculates particle source terms due to evaporation <evappe>.
!!
!! @author Andy Ackerman
!! @version Aug-2001
subroutine evapp(carma, cstate, iz, rc)

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
  integer                              :: ibin    !! bin index
  integer                              :: ielem   !! element index
  integer                              :: ig      !! source group index
  integer                              :: ip      !! source number concentration element
  integer                              :: ic  
  integer                              :: ic1     !! element of first core mass in group
  integer                              :: iecore
  integer                              :: ieto  
  integer                              :: igto  
  integer                              :: iavg
  logical                              :: evap_total
  real(kind=f)                         :: sig_mono
  real(kind=f)                         :: coretot
  real(kind=f)                         :: coremom
  real(kind=f)                         :: smf
  integer                              :: nbin

  
  ! Define criterion for monodisperse core mass distributions
  sig_mono = sqrt( ALMOST_ZERO )

  ! Loop over source groups (from which evaporation is being treated)
  do ig = 1, NGROUP

    ip = ienconc(ig)

    ! No evaporation unless particles are volatile
    if( itype(ip) .eq. I_VOLATILE )then

      ! Make sure that these always get intializaed, since they can
      ! cause problems in other parts of the code if they aren't.     
      totevap(:,ig) = .false.
      cmf(:,ig)     = 0._f

      if (pconmax(iz, ig) > FEW_PC) then
 
        ic1 = icorelem(1,ig)

        ! Loop over source bins and calculate temporary evaporation source
        ! for droplets in next smaller bin assuming no total evaporation <evdrop>
        do ibin = 1, NBIN
          evdrop = pc(iz,ibin,ip)*evaplg(ibin,ig)

          ! Check for evaporation of a sufficient number of droplets
!          if( evdrop .gt. 0._f .and. pc(iz,ibin,ip) .gt. SMALL_PC )then
          if( evdrop .gt. 0._f )then
 
            ! No cores: transfer droplets within group
            if( ic1 .eq. 0 )then
              call evap_ingrp(carma,cstate,iz,ibin,ig,ip,rc)
            else
          
              ! First core is not involatile (therefore none are)
              ! -- this is a hack until enforced/checked in setupbins() --
              ! transfer droplets within group
              !
              if( itype(ic1) .ne. I_COREMASS )then
                call evap_ingrp(carma,cstate,iz,ibin,ig,ip,rc)
              else

                ! Have cores: calculate <evcore> the amount of the source term 
                ! by number <evdrop> associated with total evaporation of secondary cores
                coretot = pc(iz,ibin,ic1)
                do ic = 2, ncore(ig)
                  iecore = icorelem(ic,ig)
                  if( itype(iecore) .eq. I_COREMASS )then
                    coretot = coretot + pc(iz,ibin,iecore)
                  endif
                enddo
                do ic = 2, ncore(ig)
                  iecore = icorelem(ic,ig)
                  if( itype(iecore) .eq. I_COREMASS )then
                    evcore(ic) = evdrop*pc(iz,ibin,iecore)/coretot  
                  endif
                enddo 
  
                ! Calculate average particle core mass and fraction
                coreavg = coretot / pc(iz,ibin,ip) 
                coreavg = min( rmass(ibin,ig), coreavg )
                cmf(ibin,ig) = coreavg / rmass(ibin,ig)
  !                 cmf(ibin,ig) = max( 0., min( ONE, cmf(ibin,ig) ) )
  
                ! Get target number concentration element and group for total evaporation
                ! and evaluate logical flags regarding position on CN bin and index of
                ! target CN bin
                ieto = ievp2elem(ic1)
                
                ! To treat internal mixtures, it is possible for the condensate to
                ! totally evaporate and have core mass, but for there not to be another
                ! group to which the core mass should go. So allow no evp2elem, but
                ! always use the in group evaporation.
                if (ieto == 0) then
                  nuc_small = .false.
                else
                  igto = igelem(ieto)
    
                  too_small = coreavg .lt. rmass(1,igto)
                  nbin = NBIN
                  too_big   = coreavg .gt. rmass(nbin,igto)
    
                  if( .not. (too_small .or. too_big) )then
                    iavg = log( coreavg / rmassmin(igto) ) / &
                           log( rmrat(igto) ) + 2
                    iavg = min( iavg, NBIN )
                  endif
    
                  ! Only consider size of evaporating cores relative to nuc_small
                  ! when treating core second moment for this particle group
                  if( if_sec_mom(ig) )then
                    nuc_small = coreavg .lt. rmass(1,igto)
                  else
                    nuc_small = .false.
                  endif
                end if
  
                ! Want total evaporation when 
                !  cores smaller than smallest nucleated 
                !  OR evaporating droplets are in bin 1
                !  OR droplets will be created with core mass fraction > 1
                evap_total = nuc_small .or. ibin .eq. 1 .or. &
                    rmrat(ig)*cmf(ibin,ig) .gt. ONE
  
                ! No core second moment: evaporate to monodisperse CN cores or within group.!
                if( .not. if_sec_mom(ig) )then
  
                  if( evap_total .and. (ieto /= 0) )then
                    call evap_mono(carma,cstate,iz,ibin,ig,iavg,ieto,igto,rc)
                  else
                    call evap_ingrp(carma,cstate,iz,ibin,ig,ip,rc)
                  endif
  
                ! Have core second moments: evaporate to mono- or polydisperse CN cores
                ! or within group.  First calculate average core second moment <coremom>, 
                ! second moment fraction <smf>, and square of the logarithm of the geometric
                ! standard deviation of the assumed core mass distribution <coresig>.
                else
  
                  coremom = pc(iz,ibin,imomelem(ig)) /  pc(iz,ibin,ip)
                  smf = coremom / rmass(ibin,ig)**2
                  coresig = log( smf / cmf(ibin,ig)**2 )
  
                  ! Want total evaporation for above reasons 
                  !  OR droplets will be created with core moment fraction > 1
                  evap_total = evap_total .or.  rmrat(ig)**2*smf .gt. ONE 
  
                  if( evap_total  .and. (ieto /= 0) )then
                
                    ! Want monodisperse total evaporation when
                    !  cores smaller than smallest nucleated 
                    !  OR evaporating core distribution is narrow
                    ! Otherwise want polydisperse total evaporation 
                    if( nuc_small .or. coresig .le. sig_mono )then
                      call evap_mono(carma,cstate,iz,ibin,ig,iavg,ieto,igto,rc)
                    else
                      call evap_poly(carma,cstate,iz,ibin,ig,iavg,ieto,igto,rc)
                    endif
  
                  ! Droplet evaporation within group
                  else
                    call evap_ingrp(carma,cstate,iz,ibin,ig,ip,rc)
                  endif
                endif      ! if_sec_mom(ig)
              endif        ! itype(ic1)
            endif          ! ic1=0 
          endif            ! evaplg > 0
        enddo              ! ibin=1,NBIN
      endif                ! enough particles
    endif                  ! volatile particles
  enddo                    ! ig=1,NGROUP

  ! Return to caller with evaporation production terms evaluated.
  return
end
