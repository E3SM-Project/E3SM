!=======================================================================
!
!BOP
!
! !MODULE: ice_aerosol - Aerosol tracer within sea ice
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$$
!
! authors Marika Holland, NCAR
!         David Bailey, NCAR
!
! !INTERFACE:
!
      module ice_aerosol
!
! !USES:
!
      use ice_kinds_mod
      use ice_constants
      use ice_fileunits
      use ice_restart, only: lenstr, restart_dir, restart_file, &
                             pointer_file, runtype
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
!
!EOP
!
      implicit none

      logical (kind=log_kind) :: & 
         restart_aero      ! if .true., read aerosol tracer restart file

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !ROUTINE: init_aerosol
!
! !DESCRIPTION:
!
!  Initialize ice aerosol tracer (call prior to reading restart data)
! 
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine init_aerosol
!
! !USES:
!
      use ice_state, only: filename_aero 
!
!EOP
!

      if (trim(filename_aero) /= 'none') restart_aero = .true.

      if (restart_aero) then
         if (trim(runtype) == 'continue') then
            call read_restart_aero
         else
            call read_restart_aero(filename_aero)
         endif
      endif

      end subroutine init_aerosol

!=======================================================================

!BOP
!
! !ROUTINE: update_aerosol
!
! !DESCRIPTION:
!
!  Increase aerosol in ice or snow surface due to deposition
!
! !REVISION HISTORY: same as module
!
! !INTERFACE:
!
      subroutine update_aerosol (nx_block, ny_block,  &
                                dt,       icells,     &
                                indxi,    indxj,      &
                                meltt,    melts,      &
                                meltb,    congel,     &
                                snoice,               &
                                fsnow,                &
                                trcrn,                &
                                aice_old,             &
                                vice_old, vsno_old,   &
                                vicen, vsnon, aicen,  &
                                faero, fsoot)
!
! !USES:
!
      use ice_domain_size, only: max_ntrcr, nilyr, nslyr, n_aero, n_aeromx
      use ice_state, only: nt_aero 
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells                ! number of cells with ice present

      integer (kind=int_kind), dimension (nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj     ! compressed indices for cells with ice

      real (kind=dbl_kind), intent(in) :: &
         dt                    ! time step

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(in) :: &
         meltt,    &
         melts,    &
         meltb,    &
         congel,   &
         snoice,   &
         fsnow,    &
         vicen,    &
         vsnon,    &
         aicen,    &
         aice_old, &
         vice_old, &
         vsno_old 

      real (kind=dbl_kind), dimension(nx_block,ny_block,n_aeromx), &
         intent(in) :: &
         faero

      real (kind=dbl_kind), dimension(nx_block,ny_block,n_aeromx), &
         intent(inout) :: &
         fsoot

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_ntrcr), &
         intent(inout) :: &
         trcrn

!
!  local variables
!
      integer (kind=int_kind) :: i, j, ij, k
      integer (kind=int_kind) :: n  ! print_points
!
      real (kind=dbl_kind), dimension(icells) :: &
         dzssl, &
         dzint, &
         dzssli, &
         dzinti

      real (kind=dbl_kind), dimension(icells) :: &
        dhs_evap, dhi_evap,  &
        dhs_melts, dhs_snoice, dhi_meltt, dhi_snoice, &
        dhi_congel, dhi_meltb 
      real (kind=dbl_kind), dimension(icells,n_aeromx) :: &
        aerotot, aerotot0   ! for diagnostics 

      real (kind=dbl_kind) :: &
         dzssl_new, &
         dzint_new, &
         dzssli_new, &
         dzinti_new, &
         dznew

      real (kind=dbl_kind), dimension(nx_block,ny_block,n_aeromx,2) :: &
         aerosno, aeroice, &
         aerosno0, aeroice0  ! for diagnostic prints

      real (kind=dbl_kind), dimension(nx_block,ny_block,n_aeromx) :: &
         fsoot_old

      real (kind=dbl_kind) :: &
         hs_old, hi_old, hslyr_old, hilyr_old, dhs, dhi, hs, hi, &
         hslyr, hilyr, sloss1, sloss2 
      real (kind=dbl_kind), dimension(n_aeromx) :: &
         kscav, kscavsi

!MH These need to be the same as in the DE code. Put in a common place?
      real (kind=dbl_kind) :: &
        hi_ssl, hs_ssl

      data hs_ssl  / .040_dbl_kind /
      data hi_ssl  / .050_dbl_kind /
      data kscav   / .03_dbl_kind, .20_dbl_kind,&
           .02_dbl_kind,.02_dbl_kind,.01_dbl_kind,.01_dbl_kind /
      data kscavsi / .03_dbl_kind, .20_dbl_kind,&
           .02_dbl_kind,.02_dbl_kind,.01_dbl_kind,.01_dbl_kind /

      aerosno(:,:,:,:) = c0
      aeroice(:,:,:,:) = c0
      aerosno0(:,:,:,:) = c0
      aeroice0(:,:,:,:) = c0
      fsoot_old(:,:,:) = fsoot(:,:,:)

      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

         hs_old=vsno_old(i,j)/aice_old(i,j)
         hi_old=vice_old(i,j)/aice_old(i,j)
         hslyr_old=hs_old/real(nslyr,kind=dbl_kind)
         hilyr_old=hi_old/real(nilyr,kind=dbl_kind)

         dzssl(ij)=min(hslyr_old/c2,hs_ssl)
         dzint(ij)=hs_old-dzssl(ij)
         dzssli(ij)=min(hilyr_old/c2,hi_ssl)
         dzinti(ij)=hi_old-dzssli(ij)

         if (aicen(i,j) > c0) then
            hs = vsnon(i,j)/aicen(i,j)
            hi = vicen(i,j)/aicen(i,j)
            dhs_melts(ij)=-melts(i,j)/aicen(i,j)
            dhi_snoice(ij)=snoice(i,j)/aicen(i,j)
            dhs_snoice(ij)=dhi_snoice(ij)*rhoi/rhos
            dhi_meltt(ij)=-meltt(i,j)/aicen(i,j)
            dhi_meltb(ij)=-meltb(i,j)/aicen(i,j)
            dhi_congel(ij)=congel(i,j)/aicen(i,j)
         else
            hs = vsnon(i,j)/aice_old(i,j)
            hi = vicen(i,j)/aice_old(i,j)
            dhs_melts(ij)=-melts(i,j)/aice_old(i,j)
            dhi_snoice(ij)=snoice(i,j)/aice_old(i,j)
            dhs_snoice(ij)=dhi_snoice(ij)*rhoi/rhos
            dhi_meltt(ij)=-meltt(i,j)/aice_old(i,j)
            dhi_meltb(ij)=-meltb(i,j)/aice_old(i,j)
            dhi_congel(ij)=congel(i,j)/aice_old(i,j)
         endif

         dhs_evap(ij)=hs-(hs_old+dhs_melts(ij)-dhs_snoice(ij)+&
             fsnow(i,j)/rhos*dt)
         dhi_evap(ij)=hi-(hi_old+dhi_meltt(ij)+dhi_meltb(ij)+ &
             dhi_congel(ij)+dhi_snoice(ij))
      enddo

      do ij = 1, icells
       i = indxi(ij)
       j = indxj(ij)
       do k=1,n_aero
         aerosno(i,j,k,:)=&
          trcrn(i,j,nt_aero+(k-1)*4  :nt_aero+(k-1)*4+1)*vsno_old(i,j)    ! aerosol in snow
         aeroice(i,j,k,:)=&
          trcrn(i,j,nt_aero+(k-1)*4+2:nt_aero+(k-1)*4+3)*vice_old(i,j)    ! aerosol in ice
         aerosno0(i,j,k,:)=aerosno(i,j,k,:)
         aeroice0(i,j,k,:)=aeroice(i,j,k,:)
         aerotot0(ij,k)=aerosno(i,j,k,2)+aerosno(i,j,k,1) &
           +aeroice(i,j,k,2)+aeroice(i,j,k,1)
       enddo
      enddo

! apply evaporation
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         dzint(ij)=dzint(ij) + min(dzssl(ij)+dhs_evap(ij),c0)
         dzssl(ij)=max(dzssl(ij)+dhs_evap(ij),c0)
         dzinti(ij)=dzinti(ij) + min(dzssli(ij)+dhi_evap(ij),c0)
         dzssli(ij)=max(dzssli(ij)+dhi_evap(ij),c0)
      enddo

!     basal ice growth
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         dzinti(ij)=dzinti(ij)+dhi_congel(ij)
      enddo

!     surface snow melt
      do ij = 1, icells
        i = indxi(ij)
        j = indxj(ij)
        if (-dhs_melts(ij) > puny) then
         do k=1,n_aero
          sloss1=c0
          sloss2=c0
          if (dzssl(ij) > puny)  &
           sloss1=kscav(k)*aerosno(i,j,k,1)  &
                 *min(-dhs_melts(ij),dzssl(ij))/dzssl(ij)
          aerosno(i,j,k,1)=aerosno(i,j,k,1)-sloss1
          if (dzint(ij) > puny)  &
           sloss2=kscav(k)*aerosno(i,j,k,2) &
                 *max(-dhs_melts(ij)-dzssl(ij),c0)/dzint(ij)
          aerosno(i,j,k,2)=aerosno(i,j,k,2)-sloss2
          fsoot(i,j,k)=fsoot(i,j,k)+(sloss1+sloss2)/dt
         enddo  ! n_aero

!      update snow thickness
         dzint(ij)=dzint(ij)+min(dzssl(ij)+dhs_melts(ij),c0)
         dzssl(ij)=max(dzssl(ij)+dhs_melts(ij),c0)

         if ( dzssl(ij) <= puny ) then ! ssl melts away
          aerosno(i,j,:,2)=aerosno(i,j,:,1)+aerosno(i,j,:,2)
          aerosno(i,j,:,1)=c0
          dzssl(ij)=max(dzssl(ij),c0)
         endif
         if (dzint(ij) <= puny ) then   ! all snow melts away
          aeroice(i,j,:,1)=&
             aeroice(i,j,:,1)+aerosno(i,j,:,1)+aerosno(i,j,:,2)
          aerosno(i,j,:,:)=c0
          dzint(ij)=max(dzint(ij),c0)
         endif
        endif
       enddo

!     surface ice melt
      do ij = 1, icells
        i = indxi(ij)
        j = indxj(ij)
        if (-dhi_meltt(ij) > puny) then
         do k=1,n_aero
          sloss1=c0
          sloss2=c0
          if (dzssli(ij) > puny)  &
           sloss1=kscav(k)*aeroice(i,j,k,1)  &
                 *min(-dhi_meltt(ij),dzssli(ij))/dzssli(ij)
          aeroice(i,j,k,1)=aeroice(i,j,k,1)-sloss1
          if (dzinti(ij) > puny)  &
           sloss2=kscav(k)*aeroice(i,j,k,2)  &
                 *max(-dhi_meltt(ij)-dzssli(ij),c0)/dzinti(ij)
          aeroice(i,j,k,2)=aeroice(i,j,k,2)-sloss2
          fsoot(i,j,k)=fsoot(i,j,k)+(sloss1+sloss2)/dt
         enddo

         dzinti(ij)=dzinti(ij)+min(dzssli(ij)+dhi_meltt(ij),c0)
         dzssli(ij)=max(dzssli(ij)+dhi_meltt(ij),c0)
         if (dzssli(ij) <= puny) then   ! ssl ice melts away
          do k=1,n_aero
           aeroice(i,j,k,2)=aeroice(i,j,k,1)+aeroice(i,j,k,2)
           aeroice(i,j,k,1)=c0
          enddo
          dzssli(ij)=max(dzssli(ij),c0)
         endif
         if (dzinti(ij) <= puny) then   ! all ice melts away
          do k=1,n_aero
           fsoot(i,j,k)=fsoot(i,j,k)  &
                       +(aeroice(i,j,k,1)+aeroice(i,j,k,2))/dt
           aeroice(i,j,k,:)=c0
          enddo
          dzinti(ij)=max(dzinti(ij),c0)
         endif
        endif
       enddo

!     basal ice melt.  Assume all soot lost in basal melt
       do ij = 1, icells
        i = indxi(ij)
        j = indxj(ij)
        if (-dhi_meltb(ij) > puny) then
         do k=1,n_aero
          sloss1=c0
          sloss2=c0
          if (dzssli(ij) > puny)  &
           sloss1=max(-dhi_meltb(ij)-dzinti(ij),c0)  &
                 *aeroice(i,j,k,1)/dzssli(ij)
          aeroice(i,j,k,1)=aeroice(i,j,k,1)-sloss1
          if (dzinti(ij) > puny)  &
           sloss2=min(-dhi_meltb(ij),dzinti(ij))  &
                 *aeroice(i,j,k,2)/dzinti(ij)
          aeroice(i,j,k,2)=aeroice(i,j,k,2)-sloss2
          fsoot(i,j,k)=fsoot(i,j,k)+(sloss1+sloss2)/dt
         enddo

         dzssli(ij) = dzssli(ij)+min(dzinti(ij)+dhi_meltb(ij), c0)
         dzinti(ij) = max(dzinti(ij)+dhi_meltb(ij), c0)           
        endif
       enddo

!     snowfall
       do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if (fsnow(i,j) > c0) &
           dzssl(ij)=dzssl(ij)+fsnow(i,j)/rhos*dt
      enddo 

!     snoice formation
      do ij = 1, icells
        i = indxi(ij)
        j = indxj(ij)
        if (dhs_snoice(ij) > puny) then
         do k=1,n_aero
          sloss1=c0
          sloss2=c0
          if (dzint(ij) > puny)  &
           sloss2 = min(dhs_snoice(ij),dzint(ij))  &
                  *aerosno(i,j,k,2)/dzint(ij)
          aerosno(i,j,k,2) = aerosno(i,j,k,2) - sloss2
          if (dzssl(ij) > puny)  &
           sloss1 = max(dhs_snoice(ij)-dzint(ij),c0)  &
                  *aerosno(i,j,k,1)/dzssl(ij)
          aerosno(i,j,k,1) = aerosno(i,j,k,1) - sloss1
          aeroice(i,j,k,1) = aeroice(i,j,k,1) &
                           + (c1-kscavsi(k))*(sloss2+sloss1)
          fsoot(i,j,k)=fsoot(i,j,k)+kscavsi(k)*(sloss2+sloss1)/dt
         enddo
         dzssl(ij)=dzssl(ij)-max(dhs_snoice(ij)-dzint(ij),c0)
         dzint(ij)=max(dzint(ij)-dhs_snoice(ij),c0)
         dzssli(ij)=dzssli(ij)+dhi_snoice(ij)
        endif
      enddo

!     aerosol deposition
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if (aicen(i,j) > c0) then
            hs = vsnon(i,j) / aicen(i,j)
         else
            hs = c0
         endif
         if (hs > hsmin) then    ! should this really be hsmin or 0? 
                                 ! should use same hsmin value as in radiation
          do k=1,n_aero
           aerosno(i,j,k,1)=aerosno(i,j,k,1) &
                           + faero(i,j,k)*dt*aicen(i,j)
          enddo
         else
          do k=1,n_aero
           aeroice(i,j,k,1)=aeroice(i,j,k,1) &
                           + faero(i,j,k)*dt*aicen(i,j)
          enddo
         endif
      enddo

!     redistribute aerosol within vertical layers
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         if (aicen(i,j) > c0) then
            hs = vsnon(i,j) / aicen(i,j)     ! new snow thickness
            hi = vicen(i,j) / aicen(i,j)     ! new ice thickness
         else
            hs = c0
            hi = c0
         endif
         if (dzssl(ij) <= puny) then   ! nothing in SSL
          do k=1,n_aero
           aerosno(i,j,k,2)=aerosno(i,j,k,2)+aerosno(i,j,k,1)
           aerosno(i,j,k,1)=c0
          enddo
         endif
         if (dzint(ij) <= puny) then   ! nothing in Snow Int
          do k=1,n_aero
           aeroice(i,j,k,1)=aeroice(i,j,k,1)+aerosno(i,j,k,2)
           aerosno(i,j,k,2)=c0
          enddo
         endif
         if (dzssli(ij) <= puny) then   ! nothing in Ice SSL
          do k=1,n_aero
           aeroice(i,j,k,2)=aeroice(i,j,k,2)+aeroice(i,j,k,1)
           aeroice(i,j,k,1)=c0
          enddo
         endif

         if (dzinti(ij) <= puny) then   ! nothing in Ice INT
          do k=1,n_aero
           fsoot(i,j,k)=fsoot(i,j,k)+&
               (aeroice(i,j,k,1)+aeroice(i,j,k,2))/dt
           aeroice(i,j,k,:)=c0
          enddo
         endif

         hslyr=hs/real(nslyr,kind=dbl_kind)
         hilyr=hi/real(nilyr,kind=dbl_kind)
         dzssl_new=min(hslyr/c2,hs_ssl)         ! ssl for snow
         dzint_new=hs-dzssl_new
         dzssli_new=min(hilyr/c2,hi_ssl)        ! ssl for ice
         dzinti_new=hi-dzssli_new
 
         if (hs > hsmin) then
          do k=1,n_aero
           dznew=min(dzssl_new-dzssl(ij),c0)
           sloss1=c0
           if (dzssl(ij) > puny) &
            sloss1=dznew*aerosno(i,j,k,1)/dzssl(ij) ! not neccesarily a loss term
           dznew=max(dzssl_new-dzssl(ij),c0)
           if (dzint(ij) > puny) &
            sloss1=sloss1+aerosno(i,j,k,2)*dznew/dzint(ij) ! not really a loss term
           aerosno(i,j,k,1) =aerosno(i,j,k,1)+sloss1 
           aerosno(i,j,k,2) =aerosno(i,j,k,2)-sloss1
          enddo
         else
          aeroice(i,j,:,1)=aeroice(i,j,:,1)  &
                          +aerosno(i,j,:,1)+aerosno(i,j,:,2)
          aerosno(i,j,:,:) = c0
         endif

         if (vicen(i,j) > puny) then ! may want a limit on hi instead?
          do k=1,n_aero
           sloss2=c0
           dznew=min(dzssli_new-dzssli(ij),c0)
           if (dzssli(ij) > puny) & 
            sloss2=dznew*aeroice(i,j,k,1)/dzssli(ij)
           dznew=max(dzssli_new-dzssli(ij),c0)
           if (dzinti(ij) > puny) & 
            sloss2=sloss2+aeroice(i,j,k,2)*dznew/dzinti(ij) ! not really a loss term
           aeroice(i,j,k,1) =aeroice(i,j,k,1)+sloss2 
           aeroice(i,j,k,2) =aeroice(i,j,k,2)-sloss2
          enddo
         else
          fsoot(i,j,:)=fsoot(i,j,:)+(aeroice(i,j,:,1)+aeroice(i,j,:,2))/dt
          aeroice(i,j,:,:) = c0
         endif

         do k=1,n_aero
          aerotot(ij,k)=aerosno(i,j,k,2)+aerosno(i,j,k,1) &
           +aeroice(i,j,k,2)+aeroice(i,j,k,1)
          if ( ( (aerotot(ij,k)-aerotot0(ij,k)) &
             - ( faero(i,j,k)*aicen(i,j)     &
             -  (fsoot(i,j,k)-fsoot_old(i,j,k)) )*dt ) > 0.00001) then
            write(nu_diag,*) 'aerosol tracer:      ',k
            write(nu_diag,*) 'aerotot-aerotot0     ',aerotot(ij,k)-aerotot0(ij,k) 
            write(nu_diag,*) 'faero-fsoot          ',faero(i,j,k)*aicen(i,j)*dt &
                                                   -(fsoot(i,j,k)-fsoot_old(i,j,k))*dt
          endif
         enddo
      enddo

!     reload tracers
      do ij = 1, icells
       i = indxi(ij)
       j = indxj(ij)
       if (vicen(i,j) > puny) &
        aeroice(i,j,:,:)=aeroice(i,j,:,:)/vicen(i,j)
       if (vsnon(i,j) > puny) &
        aerosno(i,j,:,:)=aerosno(i,j,:,:)/vsnon(i,j)
       do k=1,n_aero
        do n=1,2
         trcrn(i,j,nt_aero+(k-1)*4+n-1)=aerosno(i,j,k,n)
         trcrn(i,j,nt_aero+(k-1)*4+n+1)=aeroice(i,j,k,n)
        enddo
!       do n=1,4
!        if (trcrn(i,j,nt_aero+(k-1)*4+n-1) < puny) then
!          fsoot(i,j,k)=fsoot(i,j,k)+  &
!               trcrn(i,j,nt_aero+(k-1)*4+n-1)/dt
!          trcrn(i,j,nt_aero+(k-1)*4+n-1)=c0
!        endif
!       enddo
       enddo
      enddo

      do ij = 1, icells
       i = indxi(ij)
       j = indxj(ij)
       if (trcrn(i,j,nt_aero) < -puny .or. trcrn(i,j,nt_aero+1) < -puny    &
       .or. trcrn(i,j,nt_aero+2) < -puny .or. trcrn(i,j,nt_aero+3) < -puny) then
           write(nu_diag,*) 'MH aerosol negative in aerosol code'
           write(nu_diag,*) 'MH INT neg in aerosol my_task = ',&
                               my_task &
                               ,' printing point = ',n &
                               ,' i and j = ',i,j
           write(nu_diag,*) 'MH Int Neg aero snowssl= ',aerosno0(i,j,1,1)
           write(nu_diag,*) 'MH Int Neg aero new snowssl= ',aerosno(i,j,1,1)
           write(nu_diag,*) 'MH Int Neg aero snowint= ',aerosno0(i,j,1,2)
           write(nu_diag,*) 'MH Int Neg aero new snowint= ',aerosno(i,j,1,2)
           write(nu_diag,*) 'MH Int Neg aero ice_ssl= ',aeroice0(i,j,1,1)
           write(nu_diag,*) 'MH Int Neg aero new ice_ssl= ',aeroice(i,j,1,1)
           write(nu_diag,*) 'MH Int Neg aero ice_int= ',aeroice0(i,j,1,2)
           write(nu_diag,*) 'MH Int Neg aero new ice_int= ',aeroice(i,j,1,2)
           write(nu_diag,*) 'MH Int Neg aero aicen= ',aicen(i,j)
           write(nu_diag,*) 'MH Int Neg aero vicen= ',vicen(i,j)
           write(nu_diag,*) 'MH Int Neg aero vsnon= ',vsnon(i,j)
           write(nu_diag,*) 'MH Int Neg aero viceold= ',vice_old(i,j)
           write(nu_diag,*) 'MH Int Neg aero vsnoold= ',vsno_old(i,j)
           write(nu_diag,*) 'MH Int Neg aero melts= ',melts(i,j)
           write(nu_diag,*) 'MH Int Neg aero meltt= ',meltt(i,j)
           write(nu_diag,*) 'MH Int Neg aero meltb= ',meltb(i,j)
           write(nu_diag,*) 'MH Int Neg aero congel= ',congel(i,j)
           write(nu_diag,*) 'MH Int Neg aero snoice= ',snoice(i,j)
           write(nu_diag,*) 'MH Int Neg aero evap sno?= ',dhs_evap(ij)
           write(nu_diag,*) 'MH Int Neg aero evap ice?= ',dhi_evap(ij)
           write(nu_diag,*) 'MH Int Neg aero fsnow= ',fsnow(i,j)
           write(nu_diag,*) 'MH Int Neg aero faero= ',faero(i,j,1)
           write(nu_diag,*) 'MH Int Neg aero fsoot= ',fsoot(i,j,1)
           trcrn(i,j,nt_aero)=max(trcrn(i,j,nt_aero),c0)
           trcrn(i,j,nt_aero+1)=max(trcrn(i,j,nt_aero+1),c0)
           trcrn(i,j,nt_aero+2)=max(trcrn(i,j,nt_aero+2),c0)
           trcrn(i,j,nt_aero+3)=max(trcrn(i,j,nt_aero+3),c0)
          endif
      enddo

      end subroutine update_aerosol



!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================
!
!BOP
!
! !IROUTINE: write_restart_aero - dumps all fields required for restart
!
! !INTERFACE:
!
      subroutine write_restart_aero(filename_spec)
!
! !DESCRIPTION:
!
! Dumps all values needed for restarting
!
! !REVISION HISTORY:
!
! authors Elizabeth Hunke, LANL (original version)
!         David Bailey, NCAR
!         Marika Holland, NCAR
!
! !USES:
!
      use ice_domain_size
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_state
      use ice_read_write
      use ice_restart, only: lenstr, restart_dir, restart_file, pointer_file
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: filename

      logical (kind=log_kind) :: diag

      ! construct path/file
      if (present(filename_spec)) then
         filename = trim(filename_spec)
      else
         iyear = nyr + year_init - 1
         imonth = month
         iday = mday
         
         write(filename,'(a,a,a,i4.4,a,i2.2,a,i2.2,a,i5.5)') &
              restart_dir(1:lenstr(restart_dir)), &
              restart_file(1:lenstr(restart_file)),'.aero.', &
              iyear,'-',month,'-',mday,'-',sec
      end if
         
      ! begin writing restart data
      call ice_open(nu_dump_aero,filename,0)

      if (my_task == master_task) then
        write(nu_dump_aero) istep1,time,time_forc
        write(nu_diag,*) 'Writing ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------

      do k = 1, n_aero
      do n = 1, ncat
       call ice_write(nu_dump_aero,0,trcrn(:,:,nt_aero  +(k-1)*4,n,:),'ruf8',diag)
       call ice_write(nu_dump_aero,0,trcrn(:,:,nt_aero+1+(k-1)*4,n,:),'ruf8',diag)
       call ice_write(nu_dump_aero,0,trcrn(:,:,nt_aero+2+(k-1)*4,n,:),'ruf8',diag)
       call ice_write(nu_dump_aero,0,trcrn(:,:,nt_aero+3+(k-1)*4,n,:),'ruf8',diag)
      enddo
      enddo

      if (my_task == master_task) close(nu_dump_aero)

      end subroutine write_restart_aero

!=======================================================================
!BOP
!
! !IROUTINE: read_restart_aero - reads all fields required for restart
!
! !INTERFACE:
!
      subroutine read_restart_aero(filename_spec)
!
! !DESCRIPTION:
!
! Reads all values needed for an ice aerosol restart
!
! !REVISION HISTORY:
!
! authors Elizabeth Hunke, LANL (original version)
!         David Bailey, NCAR
!         Marika Holland, NCAR
!
! !USES:
!
      use ice_domain_size
      use ice_calendar, only: sec, month, mday, nyr, istep1, &
                              time, time_forc, idate, year_init
      use ice_state
      use ice_read_write
      use ice_restart, only: lenstr, restart_dir, restart_file, pointer_file
!
! !INPUT/OUTPUT PARAMETERS:
!
      character(len=char_len_long), intent(in), optional :: filename_spec

!EOP
!
      integer (kind=int_kind) :: &
          i, j, k, n, it, iblk, & ! counting indices
          iyear, imonth, iday     ! year, month, day

      character(len=char_len_long) :: &
         filename, filename0, string1, string2

      logical (kind=log_kind) :: &
         diag

      if (my_task == master_task) then
         ! reconstruct path/file
         if (present(filename_spec)) then
            filename = filename_spec
         else
            open(nu_rst_pointer,file=pointer_file)
            read(nu_rst_pointer,'(a)') filename0
            filename = trim(filename0)
            close(nu_rst_pointer)

            n = index(filename0,trim(restart_file))
            if (n == 0) call abort_ice('soot restart: filename discrepancy')
            string1 = trim(filename0(1:n-1))
            string2 = trim(filename0(n+lenstr(restart_file):lenstr(filename0)))
            write(filename,'(a,a,a,a)') &
               string1(1:lenstr(string1)), &
               restart_file(1:lenstr(restart_file)),'.aero', &
               string2(1:lenstr(string2))
         endif
      endif ! master_task

      call ice_open(nu_restart_aero,filename,0)

      if (my_task == master_task) then
        read(nu_restart_aero) istep1,time,time_forc
        write(nu_diag,*) 'Reading ',filename(1:lenstr(filename))
      endif

      diag = .true.

      !-----------------------------------------------------------------

      do k = 1, n_aero
      do n = 1, ncat
       call ice_read(nu_restart_aero,0,trcrn(:,:,nt_aero  +(k-1)*4,n,:),'ruf8',&
            diag,field_type=field_type_scalar,field_loc=field_loc_center)
       call ice_read(nu_restart_aero,0,trcrn(:,:,nt_aero+1+(k-1)*4,n,:),'ruf8',&
            diag,field_type=field_type_scalar,field_loc=field_loc_center)
       call ice_read(nu_restart_aero,0,trcrn(:,:,nt_aero+2+(k-1)*4,n,:),'ruf8',&
            diag,field_type=field_type_scalar,field_loc=field_loc_center)
       call ice_read(nu_restart_aero,0,trcrn(:,:,nt_aero+3+(k-1)*4,n,:),'ruf8',&
            diag,field_type=field_type_scalar,field_loc=field_loc_center)
      enddo
      enddo

      if (my_task == master_task) close(nu_restart_aero)

      end subroutine read_restart_aero

!=======================================================================

      end module ice_aerosol

!=======================================================================
