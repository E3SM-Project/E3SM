! Include shortname defintions, so that the F77 code does not have to be modified to
! reference the CARMA structure.
#include "carma_globaer.h"

!!  This routine calculates the number of sub-timesteps <ntsubsteps>
!!  for the current model spatial point.
!!
!! @author Eric Jensen
!! @version Apr-2000
subroutine nsubsteps(carma, cstate, iz, dtime_save, ntsubsteps, rc)

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
  integer, intent(in)                  :: iz          !! z index
  real(kind=f), intent(in)             :: dtime_save  !! original (not substepped) dtime
  integer, intent(inout)               :: ntsubsteps  !! suggested number of substeps
  integer, intent(inout)               :: rc          !! return code, negative indicates failure

  ! Local declarations
  integer                        :: ig      ! group index
  integer                        :: igas    ! gas index
  integer                        :: ibin    ! bin index
  integer                        :: iepart
  integer                        :: inuc
  integer                        :: ienucto
  integer                        :: ibin_small(NGROUP)
  real(kind=f)                   :: g0
  real(kind=f)                   :: g1
  real(kind=f)                   :: dmdt
  real(kind=f)                   :: dt_adv
  real(kind=f)                   :: ss
  real(kind=f)                   :: ssold
  real(kind=f)                   :: pvap
  real(kind=f)                   :: vf_max


  ! If substepping is disabled, then use one substep
  if (.not. do_substep) then
    ntsubsteps = 1
  else 
    ! Set default values
    ntsubsteps = minsubsteps
    
    ! Find the bin number of the smallest particle bin that
    ! contains a significant number of particles.
    ! Also check for significant activation of water droplets.
  
    if( ntsubsteps .lt. maxsubsteps )then
  
      do ig = 1, NGROUP
  
        if( pconmax(iz,ig) .gt. FEW_PC) then
  
          ibin_small(ig) = NBIN

          ! element of particle number concentration  
          iepart = ienconc(ig)
          
          if( itype(iepart) .eq. I_INVOLATILE ) then
  
            ! condensing gas
            igas = inucgas(ig)
  
            if (igas /= 0) then
          
              ss = max( supsatl(iz,igas), supsatlold(iz,igas) )
  
              do inuc = 1,nnuc2elem(iepart)
                ienucto = inuc2elem(inuc,iepart)
    
                if( inucproc(iepart,ienucto) .eq. I_DROPACT ) then
                  do ibin = 1, NBIN
                    if( pc(iz,ibin,iepart) / xmet(iz) / ymet(iz) / zmet(iz) .gt. conmax * pconmax(iz,ig) .and. &
                        ss .gt. scrit(iz,ibin,ig) )then
                      ntsubsteps = maxsubsteps
                    endif
                  enddo
                endif
              enddo
            endif
  
          elseif( itype(iepart) .eq. I_VOLATILE ) then
          
            do ibin = NBIN-1, 1, -1
              if( pc(iz,ibin,iepart) / xmet(iz) / ymet(iz) / zmet(iz) .gt. conmax * pconmax(iz,ig) )then
                ibin_small(ig) = ibin
              endif
            enddo
  
          endif
        endif
      enddo
    endif
  
    ! Calculate the growth rate of a particle with the mode radius for
    ! each volatile group.  The maximum time-step to use is then the
    ! mass growth rate divided by the mass bin width / 2.
    if( ntsubsteps .lt. maxsubsteps )then
  
      dt_adv = dtime_save
      do ig = 1, NGROUP
      
        ! element of particle number concentration
        iepart = ienconc(ig)
        
        ! condensing gas
        igas = igrowgas(iepart)

        if (igas /= 0) then
  
          if( pconmax(iz,ig) .gt. FEW_PC ) then
    
            if( itype(iepart) .eq. I_VOLATILE ) then
    
              if( is_grp_ice(ig) )then
                ss = supsati(iz,igas)
                pvap = pvapi(iz,igas)
              else
                ss = supsatl(iz,igas)
                pvap = pvapl(iz,igas)
              endif
    
              g0 = gro(iz,ibin_small(ig),ig)
              g1 = gro1(iz,ibin_small(ig),ig)
              dmdt = abs( pvap * ss * g0 / ( 1._f + g0*g1*pvap ) )
              
              if (dmdt /= 0._f) then
                dt_adv = min( dt_adv, dm(ibin_small(ig),ig)/dmdt )
              end if
            endif
          endif
        endif
      enddo
  
      ntsubsteps = nint(min(real(maxsubsteps, kind=f), real(dtime_save, kind=f) / dt_adv))
      ntsubsteps = max( minsubsteps, ntsubsteps )
    endif

    ! If the ice supersaturation is large enough for homogeneous freezing
    ! of sulfate aerosols, then use maximum number of substeps
    if( ntsubsteps .lt. (maxsubsteps) )then
      do ig = 1, NGROUP
      
        ! element of particle number concentration  
        iepart = ienconc(ig)

        ! condensing gas
        igas = inucgas(ig)

        if (igas /= 0) then
          
          do inuc = 1,nnuc2elem(iepart)
            ienucto = inuc2elem(inuc,iepart)
    
            if (iand(inucproc(iepart,ienucto), I_AERFREEZE) .ne. 0) then
              if( (supsati(iz,igas) .gt. 0.4_f) .and. (t(iz) .lt. 233.16_f) ) then
                ntsubsteps = maxsubsteps
              endif
            endif
          enddo
        endif
      enddo
    endif
  endif

  ! Return to caller with number of sub-timesteps evaluated.
  return
end
