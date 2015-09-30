!  SVN:$Id: ice_aerosol.F90 1012 2015-06-26 12:34:09Z eclare $
!=======================================================================

! Aerosol tracer within sea ice
!
! authors Marika Holland, NCAR
!         David Bailey, NCAR

      module ice_aerosol

      use ice_kinds_mod
      use ice_constants_colpkg, only: c0, c1, c2, puny, rhoi, rhos, hs_min, nspint
      use ice_constants, only: field_loc_center, field_type_scalar
      use ice_domain_size, only: max_aero
      use ice_fileunits, only: nu_diag

      implicit none

      private
      public :: faero_default, update_aerosol

      ! aerosol optical properties   -> band  |
      !                                       v aerosol
      ! for combined dust category, use category 4 properties
      real (kind=dbl_kind), dimension(nspint,max_aero), public :: & 
         kaer_tab, & ! aerosol mass extinction cross section (m2/kg)
         waer_tab, & ! aerosol single scatter albedo (fraction)
         gaer_tab    ! aerosol asymmetry parameter (cos(theta))

!=======================================================================

      contains

!=======================================================================

! constant values for atmospheric aerosols
!
! authors: Elizabeth Hunke, LANL

      subroutine faero_default

      use ice_flux, only: faero_atm

      faero_atm(:,:,1,:) = 1.e-15_dbl_kind ! W/m^2 s
      faero_atm(:,:,2,:) = 1.e-13_dbl_kind
      faero_atm(:,:,3,:) = 1.e-11_dbl_kind

      ! this data is used in dEdd radiation
         kaer_tab = reshape((/ &      ! aerosol mass extinction cross section (m2/kg)
          11580.61872,   5535.41835,   2793.79690, &
          25798.96479,  11536.03871,   4688.24207, &
            196.49772,    204.14078,    214.42287, &
           2665.85867,   2256.71027,    820.36024, &
            840.78295,   1028.24656,   1163.03298, &
            387.51211,    414.68808,    450.29814/), &
            (/nspint,max_aero/))
         waer_tab = reshape((/ &      ! aerosol single scatter albedo (fraction)
              0.29003,      0.17349,      0.06613, &
              0.51731,      0.41609,      0.21324, &
              0.84467,      0.94216,      0.95666, &
              0.97764,      0.99402,      0.98552, &
              0.94146,      0.98527,      0.99093, &
              0.90034,      0.96543,      0.97678/), &
              (/nspint,max_aero/))
         gaer_tab = reshape((/ &      ! aerosol asymmetry parameter (cos(theta))
              0.35445,      0.19838,      0.08857, &
              0.52581,      0.32384,      0.14970, &
              0.83162,      0.78306,      0.74375, &
              0.68861,      0.70836,      0.54171, &
              0.70239,      0.66115,      0.71983, &
              0.78734,      0.73580,      0.64411/), &
              (/nspint,max_aero/))

      end subroutine faero_default

!=======================================================================

! read atmospheric aerosols
!
! authors: Elizabeth Hunke, LANL

      subroutine faero_data

      use ice_calendar, only: month, mday, istep, sec
      use ice_domain_size, only: max_blocks
      use ice_blocks, only: nx_block, ny_block
      use ice_flux, only: faero_atm
      use ice_forcing, only: interp_coeff_monthly, read_clim_data_nc, interpolate_data

#ifdef ncdf 
      ! local parameters

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), &
         save :: &
         aero1_data    , & ! field values at 2 temporal data points
         aero2_data    , & ! field values at 2 temporal data points
         aero3_data        ! field values at 2 temporal data points

      character (char_len_long) :: & 
         aero_file,   &   ! netcdf filename
         fieldname        ! field name in netcdf file

      integer (kind=int_kind) :: & 
         ixm,ixp     , & ! record numbers for neighboring months
         maxrec      , & ! maximum record number
         recslot     , & ! spline slot for current record
         midmonth        ! middle day of month

      logical (kind=log_kind) :: readm

    !-------------------------------------------------------------------
    ! monthly data 
    !
    ! Assume that monthly data values are located in the middle of the 
    ! month.
    !-------------------------------------------------------------------

      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(month)))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      ixm  = mod(month+maxrec-2,maxrec) + 1
      ixp  = mod(month,         maxrec) + 1
      if (mday >= midmonth) ixm = 99  ! other two points will be used
      if (mday <  midmonth) ixp = 99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)

      ! Read 2 monthly values 
      readm = .false.
      if (istep==1 .or. (mday==midmonth .and. sec==0)) readm = .true.

!      aero_file = trim(atm_data_dir)//'faero.nc'   
      aero_file = '/usr/projects/climate/eclare/DATA/gx1v3/faero.nc'   

      fieldname='faero_atm001'
      call read_clim_data_nc (readm, 0,  ixm, month, ixp, &
                              aero_file, fieldname, aero1_data, &
                              field_loc_center, field_type_scalar)

      fieldname='faero_atm002'
      call read_clim_data_nc (readm, 0,  ixm, month, ixp, &
                              aero_file, fieldname, aero2_data, &
                              field_loc_center, field_type_scalar)

      fieldname='faero_atm003'
      call read_clim_data_nc (readm, 0,  ixm, month, ixp, &
                              aero_file, fieldname, aero3_data, &
                              field_loc_center, field_type_scalar)

      call interpolate_data (aero1_data, faero_atm(:,:,1,:)) ! W/m^2 s
      call interpolate_data (aero2_data, faero_atm(:,:,2,:))
      call interpolate_data (aero3_data, faero_atm(:,:,3,:))

      where (faero_atm(:,:,:,:) > 1.e20) faero_atm(:,:,:,:) = c0

#endif

      end subroutine faero_data

!=======================================================================

!  Increase aerosol in ice or snow surface due to deposition
!  and vertical cycling

      subroutine update_aerosol(dt,                  &
                                meltt,    melts,      &
                                meltb,    congel,     &
                                snoice,               &
                                fsnow,                &
                                aerosno,  aeroice,    &
                                aice_old,             &
                                vice_old, vsno_old,   &
                                vicen, vsnon, aicen,  &
                                faero_atm, faero_ocn)

      use ice_communicate, only: my_task
      use ice_domain_size, only: nilyr, nslyr, n_aero, max_aero
      use ice_shortwave, only: hi_ssl, hs_ssl
      use ice_colpkg_tracers, only: nt_aero 

      real (kind=dbl_kind), intent(in) :: &
         dt                    ! time step

      real (kind=dbl_kind), &
         intent(in) :: &
         meltt,    & ! thermodynamic melt/growth rates
         melts,    &
         meltb,    &
         congel,   &
         snoice,   &
         fsnow,    &
         vicen,    & ! ice volume (m)
         vsnon,    & ! snow volume (m)
         aicen,    & ! ice area fraction
         aice_old, & ! values prior to thermodynamic changes
         vice_old, &
         vsno_old 

      real (kind=dbl_kind), dimension(:), &
         intent(in) :: &
         faero_atm   ! aerosol deposition rate (W/m^2 s)

      real (kind=dbl_kind), dimension(:), &
         intent(inout) :: &
         faero_ocn   ! aerosol flux to ocean (W/m^2 s)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         aerosno,  aeroice    ! kg/m^2

      !  local variables
      integer (kind=int_kind) :: k, n

      real (kind=dbl_kind) :: &
         dzssl,  dzssl_new,      & ! snow ssl thickness
         dzint,  dzint_new,      & ! snow interior thickness
         dzssli, dzssli_new,     & ! ice ssl thickness
         dzinti, dzinti_new,     & ! ice interior thickness
         dznew,                  & ! tracks thickness changes
         hs, hi,                 & ! snow/ice thickness (m)
         dhs_evap, dhi_evap,     & ! snow/ice thickness change due to evap
         dhs_melts, dhi_meltt,   & ! ... due to surface melt
         dhs_snoice, dhi_snoice, & ! ... due to snow-ice formation
         dhi_congel, dhi_meltb,  & ! ... due to bottom growth, melt
         hslyr, hilyr,           & ! snow, ice layer thickness (m)
         hslyr_old, hilyr_old,   & ! old snow, ice layer thickness (m)
         hs_old, hi_old,         & ! old snow, ice thickness (m)
         sloss1, sloss2,         & ! aerosol mass loss (kg/m^2)
         ar                        ! 1/aicen(i,j)

      real (kind=dbl_kind), dimension(max_aero) :: &
         kscav, kscavsi       ! scavenging by melt water

      real (kind=dbl_kind), dimension(n_aero) :: &
         aerotot, aerotot0, & ! for conservation check
         focn_old             ! for conservation check

      real (kind=dbl_kind), dimension(n_aero,2) :: &
         aerosno0, aeroice0   ! for diagnostic prints

      ! echmod:  this assumes max_aero=6
      data kscav   / .03_dbl_kind, .20_dbl_kind, .02_dbl_kind, &
                     .02_dbl_kind, .01_dbl_kind, .01_dbl_kind / 
      data kscavsi / .03_dbl_kind, .20_dbl_kind, .02_dbl_kind, &
                     .02_dbl_kind, .01_dbl_kind, .01_dbl_kind / 

    !-------------------------------------------------------------------
    ! initialize
    !-------------------------------------------------------------------
      focn_old(:) = faero_ocn(:)
      aerosno0(:,:) = c0
      aeroice0(:,:) = c0
      
      hs_old    = vsno_old/aice_old
      hi_old    = vice_old/aice_old
      hslyr_old = hs_old/real(nslyr,kind=dbl_kind)
      hilyr_old = hi_old/real(nilyr,kind=dbl_kind)
      
      dzssl  = min(hslyr_old/c2, hs_ssl)
      dzssli = min(hilyr_old/c2, hi_ssl)
      dzint  = hs_old - dzssl
      dzinti = hi_old - dzssli
      
      if (aicen > c0) then
         ar = c1/aicen
         hs = vsnon*ar
         hi = vicen*ar
         dhs_melts  = -melts*ar
         dhi_snoice = snoice*ar
         dhs_snoice = dhi_snoice*rhoi/rhos
         dhi_meltt  = -meltt*ar
         dhi_meltb  = -meltb*ar
         dhi_congel = congel*ar
      else ! ice disappeared during time step
         hs = vsnon/aice_old
         hi = vicen/aice_old
         dhs_melts  = -melts/aice_old
         dhi_snoice = snoice/aice_old
         dhs_snoice = dhi_snoice*rhoi/rhos
         dhi_meltt  = -meltt/aice_old
         dhi_meltb  = -meltb/aice_old
         dhi_congel = congel/aice_old
      endif
      
      dhs_evap = hs - (hs_old + dhs_melts - dhs_snoice &
                              + fsnow/rhos*dt)
      dhi_evap = hi - (hi_old + dhi_meltt + dhi_meltb &
                              + dhi_congel + dhi_snoice)

      do k = 1, n_aero
         aerosno0(k,:) = aerosno(k,:)
         aeroice0(k,:) = aeroice(k,:)
         aerotot0(k) = aerosno(k,2) + aerosno(k,1) &
                     + aeroice(k,2) + aeroice(k,1)
      enddo
      
    !-------------------------------------------------------------------
    ! evaporation
    !-------------------------------------------------------------------
      dzint  = dzint  + min(dzssl  + dhs_evap, c0)
      dzinti = dzinti + min(dzssli + dhi_evap, c0)
      dzssl  = max(dzssl  + dhs_evap, c0)
      dzssli = max(dzssli + dhi_evap, c0)

    !-------------------------------------------------------------------
    ! basal ice growth
    !-------------------------------------------------------------------
      dzinti = dzinti + dhi_congel

    !-------------------------------------------------------------------
    ! surface snow melt
    !-------------------------------------------------------------------
      if (-dhs_melts > puny) then
         do k = 1, n_aero
            sloss1 = c0
            sloss2 = c0
            if (dzssl > puny)  &
                 sloss1 = kscav(k)*aerosno(k,1)  &
                                  *min(-dhs_melts,dzssl)/dzssl
            aerosno(k,1) = aerosno(k,1) - sloss1
            if (dzint > puny)  &
                 sloss2 = kscav(k)*aerosno(k,2) &
                                  *max(-dhs_melts-dzssl,c0)/dzint
            aerosno(k,2) = aerosno(k,2) - sloss2
            faero_ocn(k) = faero_ocn(k) + (sloss1+sloss2)/dt
         enddo  ! n_aero

         ! update snow thickness
         dzint=dzint+min(dzssl+dhs_melts, c0)
         dzssl=max(dzssl+dhs_melts, c0)

         if ( dzssl <= puny ) then ! ssl melts away
            aerosno(:,2) = aerosno(:,1) + aerosno(:,2)
            aerosno(:,1) = c0
            dzssl = max(dzssl, c0)
         endif
         if (dzint <= puny ) then  ! all snow melts away
            aeroice(:,1) = aeroice(:,1) &
                         + aerosno(:,1) + aerosno(:,2)
            aerosno(:,:) = c0
            dzint = max(dzint, c0)
         endif
      endif

    !-------------------------------------------------------------------
    ! surface ice melt
    !-------------------------------------------------------------------
      if (-dhi_meltt > puny) then
         do k = 1, n_aero
            sloss1 = c0
            sloss2 = c0
            if (dzssli > puny)  &
                 sloss1 = kscav(k)*aeroice(k,1)  &
                                  *min(-dhi_meltt,dzssli)/dzssli
            aeroice(k,1) = aeroice(k,1) - sloss1
            if (dzinti > puny)  &
                 sloss2 = kscav(k)*aeroice(k,2)  &
                                  *max(-dhi_meltt-dzssli,c0)/dzinti
            aeroice(k,2) = aeroice(k,2) - sloss2
            faero_ocn(k) = faero_ocn(k) + (sloss1+sloss2)/dt
         enddo
         
         dzinti = dzinti + min(dzssli+dhi_meltt, c0)
         dzssli = max(dzssli+dhi_meltt, c0)
         if (dzssli <= puny) then   ! ssl ice melts away
            do k = 1, n_aero
               aeroice(k,2) = aeroice(k,1) + aeroice(k,2)
               aeroice(k,1) = c0
            enddo
            dzssli = max(dzssli, c0)
         endif
         if (dzinti <= puny) then   ! all ice melts away
            do k = 1, n_aero
               faero_ocn(k) = faero_ocn(k)  &
                            + (aeroice(k,1)+aeroice(k,2))/dt
               aeroice(k,:)=c0
            enddo
            dzinti = max(dzinti, c0)
         endif
      endif

    !-------------------------------------------------------------------
    ! basal ice melt.  Assume all aero lost in basal melt
    !-------------------------------------------------------------------
      if (-dhi_meltb > puny) then
         do k=1,n_aero
            sloss1=c0
            sloss2=c0
            if (dzssli > puny)  &
                 sloss1 = max(-dhi_meltb-dzinti, c0)  &
                        *aeroice(k,1)/dzssli
            aeroice(k,1) = aeroice(k,1) - sloss1
            if (dzinti > puny)  &
                 sloss2 = min(-dhi_meltb, dzinti)  &
                        *aeroice(k,2)/dzinti
            aeroice(k,2) = aeroice(k,2) - sloss2
            faero_ocn(k) = faero_ocn(k) + (sloss1+sloss2)/dt
         enddo

         dzssli = dzssli + min(dzinti+dhi_meltb, c0)
         dzinti = max(dzinti+dhi_meltb, c0)           
      endif

    !-------------------------------------------------------------------
    ! snowfall
    !-------------------------------------------------------------------
      if (fsnow > c0) dzssl = dzssl + fsnow/rhos*dt

    !-------------------------------------------------------------------
    ! snow-ice formation
    !-------------------------------------------------------------------
      if (dhs_snoice > puny) then
         do k = 1, n_aero
            sloss1 = c0
            sloss2 = c0
            if (dzint > puny)  &
                 sloss2 = min(dhs_snoice, dzint)  &
                        *aerosno(k,2)/dzint
            aerosno(k,2) = aerosno(k,2) - sloss2
            if (dzssl > puny)  &
                 sloss1 = max(dhs_snoice-dzint, c0)  &
                        *aerosno(k,1)/dzssl
            aerosno(k,1) = aerosno(k,1) - sloss1
            aeroice(k,1) = aeroice(k,1) &
                         + (c1-kscavsi(k))*(sloss2+sloss1)
            faero_ocn(k) = faero_ocn(k) &
                         + kscavsi(k)*(sloss2+sloss1)/dt
         enddo
         dzssl  = dzssl - max(dhs_snoice-dzint, c0)
         dzint  = max(dzint-dhs_snoice, c0)
         dzssli = dzssli + dhi_snoice
      endif

    !-------------------------------------------------------------------
    ! aerosol deposition
    !-------------------------------------------------------------------
      if (aicen > c0) then
         hs = vsnon * ar
      else
         hs = c0
      endif
      if (hs > hs_min) then    ! should this really be hs_min or 0? 
         ! should use same hs_min value as in radiation
         do k=1,n_aero
            aerosno(k,1) = aerosno(k,1) &
                         + faero_atm(k)*dt*aicen
         enddo
      else
         do k=1,n_aero
            aeroice(k,1) = aeroice(k,1) &
                         + faero_atm(k)*dt*aicen
         enddo
      endif

    !-------------------------------------------------------------------
    ! redistribute aerosol within vertical layers
    !-------------------------------------------------------------------
      if (aicen > c0) then
         hs = vsnon  * ar     ! new snow thickness
         hi = vicen  * ar     ! new ice thickness
      else
         hs = c0
         hi = c0
      endif
      if (dzssl <= puny) then   ! nothing in SSL
         do k=1,n_aero
            aerosno(k,2) = aerosno(k,2) + aerosno(k,1)
            aerosno(k,1) = c0
         enddo
      endif
      if (dzint <= puny) then   ! nothing in Snow Int
         do k = 1, n_aero
            aeroice(k,1) = aeroice(k,1) + aerosno(k,2)
            aerosno(k,2) = c0
         enddo
      endif
      if (dzssli <= puny) then  ! nothing in Ice SSL
         do k = 1, n_aero
            aeroice(k,2) = aeroice(k,2) + aeroice(k,1)
            aeroice(k,1) = c0
         enddo
      endif
      
      if (dzinti <= puny) then  ! nothing in Ice INT
         do k = 1, n_aero
            faero_ocn(k) = faero_ocn(k) &
                         + (aeroice(k,1)+aeroice(k,2))/dt
            aeroice(k,:)=c0
         enddo
      endif
      
      hslyr      = hs/real(nslyr,kind=dbl_kind)
      hilyr      = hi/real(nilyr,kind=dbl_kind)
      dzssl_new  = min(hslyr/c2, hs_ssl)
      dzssli_new = min(hilyr/c2, hi_ssl)
      dzint_new  = hs - dzssl_new
      dzinti_new = hi - dzssli_new

      if (hs > hs_min) then
         do k = 1, n_aero
            dznew = min(dzssl_new-dzssl, c0)
            sloss1 = c0
            if (dzssl > puny) &
                 sloss1 = dznew*aerosno(k,1)/dzssl ! not neccesarily a loss
            dznew = max(dzssl_new-dzssl, c0)
            if (dzint > puny) &
                 sloss1 = sloss1 + aerosno(k,2)*dznew/dzint
            aerosno(k,1) = aerosno(k,1) + sloss1 
            aerosno(k,2) = aerosno(k,2) - sloss1
         enddo
      else
         aeroice(:,1) = aeroice(:,1)  &
                      + aerosno(:,1) + aerosno(:,2)
         aerosno(:,:) = c0
      endif
      
      if (vicen > puny) then ! may want a limit on hi instead?
         do k = 1, n_aero
            sloss2 = c0
            dznew = min(dzssli_new-dzssli, c0)
            if (dzssli > puny) & 
                 sloss2 = dznew*aeroice(k,1)/dzssli
            dznew = max(dzssli_new-dzssli, c0)
            if (dzinti > puny) & 
                 sloss2 = sloss2 + aeroice(k,2)*dznew/dzinti
            aeroice(k,1) = aeroice(k,1) + sloss2 
            aeroice(k,2) = aeroice(k,2) - sloss2
         enddo
      else
         faero_ocn(:) = faero_ocn(:) + (aeroice(:,1)+aeroice(:,2))/dt
         aeroice(:,:) = c0
      endif
      
    !-------------------------------------------------------------------
    ! check conservation
    !-------------------------------------------------------------------
      do k = 1, n_aero
         aerotot(k) = aerosno(k,2) + aerosno(k,1) &
                    + aeroice(k,2) + aeroice(k,1)
         if ((aerotot(k)-aerotot0(k)) &
              - (   faero_atm(k)*aicen &
              - (faero_ocn(k)-focn_old(k)) )*dt  > puny) then
            
            write(nu_diag,*) 'aerosol tracer:  ',k
            write(nu_diag,*) 'aerotot-aerotot0 ',aerotot(k)-aerotot0(k) 
            write(nu_diag,*) 'faero_atm-faero_ocn      ', &
                 (faero_atm(k)*aicen-(faero_ocn(k)-focn_old(k)))*dt
         endif
      enddo

    !-------------------------------------------------------------------
    ! check for negative values
    !-------------------------------------------------------------------

!echmod:  note that this does not test or fix all aero tracers         
      if (aeroice(1,1) < -puny .or. &
          aeroice(1,2) < -puny .or. &
          aerosno(1,1) < -puny .or. &
          aerosno(1,2) < -puny) then

         write(nu_diag,*) 'MH aerosol negative in aerosol code'
         write(nu_diag,*) 'MH INT neg in aerosol my_task = ',&
              my_task, &
              ' printing point = ',n
         write(nu_diag,*) 'MH Int Neg aero snowssl= '    ,aerosno0(1,1)
         write(nu_diag,*) 'MH Int Neg aero new snowssl= ',aerosno (1,1)
         write(nu_diag,*) 'MH Int Neg aero snowint= '    ,aerosno0(1,2)
         write(nu_diag,*) 'MH Int Neg aero new snowint= ',aerosno (1,2)
         write(nu_diag,*) 'MH Int Neg aero ice_ssl= '    ,aeroice0(1,1)
         write(nu_diag,*) 'MH Int Neg aero new ice_ssl= ',aeroice (1,1)
         write(nu_diag,*) 'MH Int Neg aero ice_int= '    ,aeroice0(1,2)
         write(nu_diag,*) 'MH Int Neg aero new ice_int= ',aeroice (1,2)
         write(nu_diag,*) 'MH Int Neg aero aicen= '      ,aicen   
         write(nu_diag,*) 'MH Int Neg aero vicen= '      ,vicen   
         write(nu_diag,*) 'MH Int Neg aero vsnon= '      ,vsnon   
         write(nu_diag,*) 'MH Int Neg aero viceold= '    ,vice_old
         write(nu_diag,*) 'MH Int Neg aero vsnoold= '    ,vsno_old
         write(nu_diag,*) 'MH Int Neg aero melts= '      ,melts   
         write(nu_diag,*) 'MH Int Neg aero meltt= '      ,meltt   
         write(nu_diag,*) 'MH Int Neg aero meltb= '      ,meltb   
         write(nu_diag,*) 'MH Int Neg aero congel= '     ,congel  
         write(nu_diag,*) 'MH Int Neg aero snoice= '     ,snoice  
         write(nu_diag,*) 'MH Int Neg aero evap sno?= '  ,dhs_evap
         write(nu_diag,*) 'MH Int Neg aero evap ice?= '  ,dhi_evap
         write(nu_diag,*) 'MH Int Neg aero fsnow= '      ,fsnow
         write(nu_diag,*) 'MH Int Neg aero faero_atm= '  ,faero_atm(1)
         write(nu_diag,*) 'MH Int Neg aero faero_ocn= '  ,faero_ocn(1)

!echmod:  note that this does not test or fix all aero tracers         
         aeroice(1,1) = max(aeroice(1,1), c0)
         aeroice(1,2) = max(aeroice(1,2), c0)
         aerosno(1,1) = max(aerosno(1,1), c0)
         aerosno(1,2) = max(aerosno(1,2), c0)

      endif

      end subroutine update_aerosol

!=======================================================================

      end module ice_aerosol

!=======================================================================
