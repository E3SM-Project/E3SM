#include "carma_globaer.h"

!! This routine calculates coagulation production terms <coagpe>.
!! See [Jacobson, et al., Atmos. Env., 28, 1327, 1994] for details
!! on the coagulation algorithm.
!!
!! @author Eric Jensen
!! @version Oct-1995
subroutine coagp(carma, cstate, ibin, ielem, rc)

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
  integer, intent(in)                  :: ibin    !! bin index
  integer, intent(in)                  :: ielem   !! element index
  integer, intent(inout)               :: rc      !! return code, negative indicates failure

  ! Local Variables
  integer                        :: iz
  integer                        :: igrp
  integer                        :: i_pkern
  integer                        :: iquad
  integer                        :: ig
  integer                        :: jg
  integer                        :: i
  integer                        :: j
  integer                        :: iefrom
  integer                        :: iefrom_cm
  integer                        :: je
  integer                        :: je_cm
  integer                        :: ic
  integer                        :: iecore
  real(kind=f)                   :: totmass
  real(kind=f)                   :: rmasscore
  real(kind=f)                   :: fracmass
  real(kind=f)                   :: elemass
  real(kind=f)                   :: rmi
  real(kind=f)                   :: rmj


  ! Definition of i,j,k,n used in comments: colision between i and j bins
  ! yields particle between bins k and k+1.  Production in bin n is calculated.
  
  
  ! Determine group that particles are produced in
  igrp = igelem(ielem)
  
  ! Particle number production
  !
  ! Coagulation between particle in group <ig> bin <i> with particle in
  ! group <jg> bin <j> results in particle with mass between bins k and k+1.
  ! First, loop over group-bin quads <ig,i,jg,j> resulting in production in
  ! bin <ibin> = k+1.  The set of quads <igup,jgup,iup,jup> is
  ! defined in setupcoag.
  
  do iquad = 1, npairu(igrp,ibin)

    ig = igup(igrp,ibin,iquad)           ! source group
    jg = jgup(igrp,ibin,iquad)           ! source group
    i  = iup(igrp,ibin,iquad)            ! source bin
    j  = jup(igrp,ibin,iquad)            ! source bin
  
    iefrom = icoagelem(ielem,ig)         ! source element for <i> particle
    
    if( if_sec_mom(igrp) )then
      iefrom_cm = icoagelem_cm(ielem,ig)        ! core mass moment source element
    endif

    ! If <iefrom> = 0 then there is no contribution to production
    if( iefrom .ne. 0 ) then

      je = ienconc(jg)                   ! source element for <j> particle
  
      if( if_sec_mom(igrp) )then
        je_cm = icoagelem_cm(ielem,jg)           ! core mass moment source element
      endif
        
      ! If ielem is core mass type and <ig> is a CN type and <ig> is different
      ! from <igrp>, then we must multiply production by mass
      ! per particle (<rmi>) of element <icoagelem>.  (this is <rmass> for all source
      ! elements except particle number concentration in a multicomponent CN group).
      do iz = 1, NZ

        ! Bypass calculation if few source particles present
        if( pconmax(iz,ig) .gt. FEW_PC .and. &
            pconmax(iz,jg) .gt. FEW_PC )then

          rmi = 1._f
          i_pkern = 1

          if( itype(ielem) .eq. I_COREMASS .or. &
              itype(ielem) .eq. I_VOLCORE )then          ! core mass

            i_pkern = 3        ! Use different kernel for core mass prod.

            if( ( itype(ienconc(ig)) .eq. I_INVOLATILE .or. &
                  itype(ienconc(ig)) .eq. I_VOLATILE ) &
                .and. ig .ne. igrp ) then 

                ! CN source and ig different from igrp

              if( ncore(ig) .eq. 0 )then        ! No cores in source group

                if(icomp(ienconc(ig)) .eq. icomp(ielem)) then
                  rmi = rmass(i,ig)
                else
                  rmi = 0._f
                endif

              elseif( itype(iefrom) .eq. I_INVOLATILE .or. &
                      itype(iefrom) .eq. I_VOLATILE ) then

                !  Source element is number concentration elem of mixed CN group
                totmass  = pc(iz,i,iefrom) * rmass(i,ig)
                rmasscore = pc(iz,i,icorelem(1,ig))
                
                do ic = 2,ncore(ig)
                  iecore = icorelem(ic,ig)
                  rmasscore = rmasscore + pc(iz,i,iecore)
                enddo
                
                fracmass = 1._f - rmasscore/totmass
                elemass  = fracmass * rmass(i,ig)
                rmi = elemass
              endif
            endif  ! ig is a CCN and not igrp

          elseif( itype(ielem) .eq. I_CORE2MOM )then      ! core mass^2

            i_pkern = 5    ! Use different kernel for core mass^2 production
            rmj = 1._f

            if( itype(ienconc(ig)) .eq. I_INVOLATILE ) then
              rmi = rmass(i,ig)
              rmj = rmass(j,jg)
            endif

          endif  ! itype(ielem) is a coremass or core2mom

          ! For each spatial grid point, sum up coagulation production
          ! contributions from each quad.
          if( itype(ielem) .ne. I_CORE2MOM )then
            coagpe(iz,ibin,ielem) = coagpe(iz,ibin,ielem) + &
              pc(iz,i,iefrom)*pcl(iz,j,je)*rmi * &
              ckernel(iz,i,j,ig,jg) * &
              pkernel(i,j,ig,jg,igrp,i_pkern) 
          else
            coagpe(iz,ibin,ielem) = coagpe(iz,ibin,ielem) + &
              ( pc(iz,i,iefrom)*pcl(iz,j,je)*rmi**2 + &
              pc(iz,i,iefrom_cm)*rmi* &
              pcl(iz,j,je_cm)*rmj ) * &
              ckernel(iz,i,j,ig,jg) * &
              pkernel(i,j,ig,jg,igrp,i_pkern) 
          endif
        endif    ! end of ( pconmax .gt. FEW_PC )
      enddo    ! iz = 1, NZ
    endif  ! iefrom .ne. 0
  enddo  ! iquad 

  ! Next, loop over group-bin quads for production in bin <ibin> = k from
  ! bin <i> due to collision between bins <i> and <j>.
  ! Production will only occur if either k != <i> or igrp != <ig>
  do iquad = 1, npairl(igrp,ibin)

    ig = iglow(igrp,ibin,iquad)
    jg = jglow(igrp,ibin,iquad)
    i  = ilow(igrp,ibin,iquad)
    j  = jlow(igrp,ibin,iquad)
  
    iefrom = icoagelem(ielem,ig)          ! source element for <i> particle
  
    if( if_sec_mom(igrp) )then
      iefrom_cm = icoagelem_cm(ielem,ig)        ! core mass moment source element
    endif

    if( iefrom .ne. 0 ) then
  
      je = ienconc(jg)                    ! source element for <j> particle
  
      if( if_sec_mom(igrp) )then
        je_cm = icoagelem_cm(ielem,jg)           ! core mass moment source element
      endif
        
      do iz = 1, NZ

        ! Bypass calculation if few particles present
        if( pconmax(iz,ig) .gt. FEW_PC .and. &
            pconmax(iz,jg) .gt. FEW_PC )then

          rmi = 1._f
          i_pkern = 2

          if( itype(ielem) .eq. I_COREMASS .or. &
              itype(ielem) .eq. I_VOLCORE )then          ! core mass

            i_pkern = 4     ! Use different kernel for core mass production

            if( ( itype(ienconc(ig)) .eq. I_INVOLATILE .or. &
                  itype(ienconc(ig)) .eq. I_VOLATILE ) &
                .and. ig .ne. igrp ) then 

              ! CN source and ig different from igrp

              if( ncore(ig) .eq. 0 )then          ! No cores in source group
                rmi = rmass(i,ig)

              elseif( itype(iefrom) .eq. I_INVOLATILE .or. &
                      itype(iefrom) .eq. I_VOLATILE ) then

                ! Source element is number concentration elem of mixed CN group

                totmass  = pc(iz,i,iefrom) * rmass(i,ig)
                rmasscore = pc(iz,i,icorelem(1,ig))
                do ic = 2,ncore(ig)
                  iecore = icorelem(ic,ig)
                  rmasscore = rmasscore + pc(iz,i,iecore)
                enddo
                fracmass = 1._f - rmasscore/totmass
                elemass  = fracmass * rmass(i,ig)
                rmi = elemass

              endif  ! pure CN group or CN group w/ cores

            endif  ! src group is CN and different from the target group

          elseif( itype(ielem) .eq. I_CORE2MOM )then      ! core mass^2

            i_pkern = 6   ! Use different kernel for core mass^2 production
            rmj = 1._f
            if( itype(ienconc(ig)) .eq. I_INVOLATILE ) then
              rmi = rmass(i,ig)
              rmj = rmass(j,jg)
            endif
          endif  ! itype(ielem)

          if( itype(ielem) .ne. I_CORE2MOM )then
           coagpe(iz,ibin,ielem) = coagpe(iz,ibin,ielem) + &
                pc(iz,i,iefrom)*pcl(iz,j,je)*rmi * &
                ckernel(iz,i,j,ig,jg) * &
                pkernel(i,j,ig,jg,igrp,i_pkern)
          else
           coagpe(iz,ibin,ielem) = coagpe(iz,ibin,ielem) + &
                ( pc(iz,i,iefrom)*pcl(iz,j,je)*rmi**2 + &
                  pc(iz,i,iefrom_cm)*rmi* &
                  pcl(iz,j,je_cm)*rmj ) * &
                ckernel(iz,i,j,ig,jg) * &
                pkernel(i,j,ig,jg,igrp,i_pkern)
          endif
        endif    ! end of ( pconmax .gt. FEW_PC )
      enddo    ! iz = 1, NZ
    endif      ! end of (iefrom .ne. 0)
  enddo  ! iquad

  ! Return to caller with coagulation production terms evaluated.
 
  return
end
