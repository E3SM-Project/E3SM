!  SVN:$Id: ice_algae.F90 1012 2015-06-26 12:34:09Z eclare $
!=======================================================================
!
! Compute biogeochemistry in the skeletal layer
!
! authors: Nicole Jeffery, LANL
!          Scott Elliot,   LANL
!          Elizabeth C. Hunke, LANL
!
      module ice_algae

      use ice_kinds_mod
      use ice_domain_size, only: nblyr, nilyr, max_blocks
      use ice_blocks, only: nx_block, ny_block
      use ice_fileunits, only: nu_diag, nu_restart_bgc, nu_rst_pointer, &
          nu_dump_bgc, flush_fileunit
      use ice_read_write, only: ice_open, ice_read, ice_write
      use ice_communicate, only: my_task, master_task
      use ice_exit, only: abort_ice
      use ice_zbgc_shared ! everything
      use ice_state, only: vicen, vice, trcr
      use ice_colpkg_tracers, only: ntrcr, nt_bgc_Am_sk, &
          nt_bgc_C_sk, nt_bgc_chl_sk, nt_bgc_DMS_sk, nt_bgc_DMSPd_sk, &
          nt_bgc_DMSPp_sk, nt_bgc_N_sk, nt_bgc_Nit_sk, nt_bgc_Sil_sk

      implicit none

      private
      public :: get_forcing_bgc, bgc_diags, write_restart_bgc, &
                algal_dyn, read_restart_bgc, &
                skl_biogeochemistry

      real (kind=dbl_kind), parameter, private :: &
         R_Si2N   = 1.5_dbl_kind ! algal Si to N (mole/mole)

!=======================================================================

      contains

!=======================================================================
!
! Read and interpolate annual climatologies of silicate and nitrate.
! Restore model quantities to data if desired.
!
! author: Elizabeth C. Hunke, LANL

      subroutine get_forcing_bgc

      use ice_calendar, only: dt, istep, mday, month, sec
      use ice_constants, only: field_loc_center, field_type_scalar
      use ice_domain, only: nblocks
      use ice_flux, only: sss
      use ice_forcing, only: trestore, trest, &
          read_clim_data, interpolate_data, interp_coeff_monthly

      integer (kind=int_kind) :: &
         i, j, iblk  , & ! horizontal indices
         ixm,ixp     , & ! record numbers for neighboring months
         maxrec      , & ! maximum record number
         recslot     , & ! spline slot for current record
         midmonth        ! middle day of month

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
         nitdat      , & ! data value toward which nitrate is restored
         sildat          ! data value toward which silicate is restored

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks) :: &
         nit_data, & ! field values at 2 temporal data points
         sil_data

      logical (kind=log_kind) :: readm

      if (trim(nit_data_type) == 'clim'.or. &
          trim(sil_data_type) == 'clim') then

         nit_file = trim(bgc_data_dir)//'nitrate_WOA2005_surface_monthly' ! gx1 only
         sil_file = trim(bgc_data_dir)//'silicate_WOA2005_surface_monthly' ! gx1 only

         if (my_task == master_task .and. istep == 1) then
         if (trim(sil_data_type)=='clim' .AND. tr_bgc_Sil_sk) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'silicate data interpolated to timestep:'
            write (nu_diag,*) trim(sil_file)
         endif
         if (trim(nit_data_type)=='clim' .AND. tr_bgc_Nit_sk) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'nitrate data interpolated to timestep:'
            write (nu_diag,*) trim(nit_file)
            if (restore_bgc) write (nu_diag,*) &
              'bgc restoring timescale (days) =', trestore
         endif
         endif                     ! my_task, istep

    !-------------------------------------------------------------------
    ! monthly data
    !
    ! Assume that monthly data values are located in the middle of the
    ! month.
    !-------------------------------------------------------------------

         midmonth = 15          ! data is given on 15th of every month
!!!      midmonth = fix(p5 * real(daymo(month)))  ! exact middle

         ! Compute record numbers for surrounding months
         maxrec = 12
         ixm  = mod(month+maxrec-2,maxrec) + 1
         ixp  = mod(month,         maxrec) + 1
         if (mday >= midmonth) ixm = -99 ! other two points will be used
         if (mday <  midmonth) ixp = -99

         ! Determine whether interpolation will use values 1:2 or 2:3
         ! recslot = 2 means we use values 1:2, with the current value (2)
         !  in the second slot
         ! recslot = 1 means we use values 2:3, with the current value (2)
         !  in the first slot
         recslot = 1            ! latter half of month
         if (mday < midmonth) recslot = 2 ! first half of month

         ! Find interpolation coefficients
         call interp_coeff_monthly (recslot)

         readm = .false.
         if (istep==1 .or. (mday==midmonth .and. sec==0)) readm = .true.

      endif   ! sil/nit_data_type

    !-------------------------------------------------------------------
    ! Read two monthly silicate values and interpolate.
    ! Restore toward interpolated value.
    !-------------------------------------------------------------------

      if (trim(sil_data_type)=='clim'  .AND. tr_bgc_Sil_sk) then
         call read_clim_data (readm, 0, ixm, month, ixp, &
                              sil_file, sil_data, &
                              field_loc_center, field_type_scalar)
         call interpolate_data (sil_data, sildat)

         if (restore_bgc) then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sil(i,j,iblk) = sil(i,j,iblk)  &
                             + (sildat(i,j,iblk)-sil(i,j,iblk))*dt/trest
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         endif

      endif

    !-------------------------------------------------------------------
    ! Read two monthly nitrate values and interpolate.
    ! Restore toward interpolated value.
    !-------------------------------------------------------------------

      if (trim(nit_data_type)=='clim' .AND. tr_bgc_Nit_sk) then
         call read_clim_data (readm, 0, ixm, month, ixp, &
                              nit_file, nit_data, &
                              field_loc_center, field_type_scalar)
         call interpolate_data (nit_data, nitdat)

         if (restore_bgc) then
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               nit(i,j,iblk) = nit(i,j,iblk)  &
                             + (nitdat(i,j,iblk)-nit(i,j,iblk))*dt/trest
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
         endif

      elseif (trim(nit_data_type) == 'sss' .AND. tr_bgc_Nit_sk) then
          
         !$OMP PARALLEL DO PRIVATE(iblk,i,j)
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               nit(i,j,iblk) =  sss(i,j,iblk)        
            enddo
            enddo
         enddo
         !$OMP END PARALLEL DO
                 
      endif

      end subroutine get_forcing_bgc

!=======================================================================
!
! skeletal layer biochemistry
! 
      subroutine skl_biogeochemistry (nx_block, ny_block,  &
                                      icells,   dt,        &
                                      indxi,    indxj,     &
                                      nbtrcr,              &
                                      flux_bio, ocean_bio, &
                                      hmix,     aicen,     &
                                      meltb,    congel,    &
                                      fswthru,  first_ice, &
                                      trcrn,    grow_Cn)

      use ice_constants, only: p5, p05, p1, c1, c0, puny

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of cells with aicen > puny 
         nbtrcr                ! number of bgc tracers

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj    ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), intent(in) :: &
         dt         ! time step 

      logical (kind=log_kind), dimension (nx_block,ny_block), &
         intent(in) :: &
         first_ice  ! initialized values should be used

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         hmix   , & ! mixed layer depth (m)
         aicen  , & ! ice area 
         meltb  , & ! bottom ice melt (m)
         congel , & ! bottom ice growth (m)
         fswthru    ! shortwave passing through ice to ocean

      real (kind=dbl_kind), dimension(nx_block,ny_block,ntrcr), &
         intent(inout) :: &
         trcrn      ! tracer array
   
      ! history variables

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(out) :: &
         grow_Cn    ! specific growth (1/s)  

      real (kind=dbl_kind), dimension (nx_block,ny_block,nbtrcr), &
         intent(inout) :: &
         flux_bio,& ! ocean tracer flux (mmol/m^2/s) positive into ocean
         ocean_bio  ! ocean tracer concentration (mmol/m^3)

      !  local variables

      integer (kind=int_kind) :: i, j, ij, nn

      real (kind=dbl_kind), dimension(icells,nbtrcr):: &
         react  , & ! biological sources and sinks (mmol/m^3)
         cinit  , & ! initial brine concentration (mmol/m^2)
         congel_alg ! congelation flux contribution to ice algae (mmol/m^2 s) 

      real (kind=dbl_kind), dimension (nbtrcr):: &
         flux_bio_temp, & ! tracer flux to ocean (mmol/m^2 s)
         PVflag       , & ! 1 for tracers that flow with the brine, 0 otherwise
         cling            ! 1 for tracers that cling, 0 otherwise

      real (kind=dbl_kind), parameter :: &
         PVc = 1.e-6_dbl_kind   , & ! type 'constant' piston velocity (m/s) 
         PV_scale_growth = p5   , & ! scale factor in Jin PV during ice growth
         PV_scale_melt = p05    , & ! scale factor in Jin PV during ice melt
         MJ1 = 9.667e-9_dbl_kind, & ! coefficients in Jin 2008 (m/s)
         MJ2 = 38.8_dbl_kind    , & ! 4.49e-4_dbl_kind*secday (unitless)
         MJ3 = 1.04e7_dbl_kind  , & ! 1.39e-3_dbl_kind*secday^2  (s/m)
         PV_frac_max = 0.9_dbl_kind ! max piston velocity coefficient

      real (kind=dbl_kind), dimension(icells) :: &
         PVt       , & ! type 'Jin2006' piston velocity (m/s) 
         ice_growth, & ! Jin2006 defn: congel rate or bottom melt rate (m/s)
         f_meltn       ! vertical melt fraction of skeletal layer in dt

      real (kind=dbl_kind):: &
         rphi_sk   , & ! 1 / skeletal layer porosity
         cinit_tmp     ! temporary variable for concentration (mmol/m^2)

      !-------------------------------------------------------------------
      ! Initialize 
      !-------------------------------------------------------------------
     
      do nn = 1, nbtrcr          
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            cinit     (ij,nn) = c0
            congel_alg(ij,nn) = c0
            react     (ij,nn) = c0
         enddo 
         PVflag(nn) = c1
         if (bgc_tracer_type(nn) < p5) PVflag(nn) = c0
         cling(nn) = c0
      enddo 
      cling(nlt_bgc_N) = c1

      rphi_sk = c1/phi_sk

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
 
         PVt       (ij) = c0
         f_meltn   (ij) = c0
         ice_growth(ij) = (congel(i,j)-meltb(i,j))/dt

         if (first_ice(i,j)) then     
               trcrn(i,j,nt_bgc_N_sk)     = ocean_bio(i,j,nlt_bgc_N)    *sk_l*rphi_sk
            if (tr_bgc_Nit_sk) &
               trcrn(i,j,nt_bgc_Nit_sk)   = ocean_bio(i,j,nlt_bgc_NO)   *sk_l*rphi_sk
            if (tr_bgc_Am_sk)  &
               trcrn(i,j,nt_bgc_Am_sk)    = ocean_bio(i,j,nlt_bgc_NH)   *sk_l*rphi_sk
            if (tr_bgc_Sil_sk) &
               trcrn(i,j,nt_bgc_Sil_sk)   = ocean_bio(i,j,nlt_bgc_Sil)  *sk_l*rphi_sk
            if (tr_bgc_C_sk)   &
               trcrn(i,j,nt_bgc_C_sk)     = ocean_bio(i,j,nlt_bgc_C)    *sk_l*rphi_sk
            if (tr_bgc_chl_sk) &
               trcrn(i,j,nt_bgc_chl_sk)   = ocean_bio(i,j,nlt_bgc_chl)  *sk_l*rphi_sk
            if (tr_bgc_DMSPp_sk) &
               trcrn(i,j,nt_bgc_DMSPp_sk) = ocean_bio(i,j,nlt_bgc_DMSPp)*sk_l*rphi_sk
            if (tr_bgc_DMSPd_sk) &
               trcrn(i,j,nt_bgc_DMSPd_sk) = ocean_bio(i,j,nlt_bgc_DMSPd)*sk_l*rphi_sk
            if (tr_bgc_DMS_sk) &
               trcrn(i,j,nt_bgc_DMS_sk)   = ocean_bio(i,j,nlt_bgc_DMS)  *sk_l*rphi_sk
         endif ! first_ice

                              cinit(ij,nlt_bgc_N)     = trcrn(i,j,nt_bgc_N_sk)    
         if (tr_bgc_Nit_sk)   cinit(ij,nlt_bgc_NO)    = trcrn(i,j,nt_bgc_Nit_sk)  
         if (tr_bgc_Am_sk)    cinit(ij,nlt_bgc_NH)    = trcrn(i,j,nt_bgc_Am_sk)   
         if (tr_bgc_Sil_sk)   cinit(ij,nlt_bgc_Sil)   = trcrn(i,j,nt_bgc_Sil_sk)  
         if (tr_bgc_C_sk)     cinit(ij,nlt_bgc_C)     = trcrn(i,j,nt_bgc_C_sk)    
         if (tr_bgc_chl_sk)   cinit(ij,nlt_bgc_chl)   = trcrn(i,j,nt_bgc_chl_sk)  
         if (tr_bgc_DMSPp_sk) cinit(ij,nlt_bgc_DMSPp) = trcrn(i,j,nt_bgc_DMSPp_sk)
         if (tr_bgc_DMSPd_sk) cinit(ij,nlt_bgc_DMSPd) = trcrn(i,j,nt_bgc_DMSPd_sk)
         if (tr_bgc_DMS_sk)   cinit(ij,nlt_bgc_DMS)   = trcrn(i,j,nt_bgc_DMS_sk)  
            
      enddo     ! ij

      do nn = 1, nbtrcr          
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            if (cinit(ij,nn) < c0) then
               write(nu_diag,*)'initial sk_bgc < 0, ij,nn,nbtrcr,cinit(ij,nn)', &
                    ij,nn,nbtrcr,cinit(ij,nn)
               call abort_ice ('ice_bgc.F90: BGC error1')
            endif
         enddo 
      enddo 

      !-------------------------------------------------------------------
      ! 'Jin2006':
      ! 1. congel/melt dependent piston velocity (PV) for growth and melt
      ! 2. If congel > melt use 'congel'; if melt > congel use 'melt'
      ! 3. For algal N, PV for ice growth only provides a seeding concentration 
      ! 4. Melt affects nutrients and algae in the same manner through PV(melt)
      !-------------------------------------------------------------------

      if (trim(bgc_flux_type) == 'Jin2006') then  

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)
 
            if (ice_growth(ij) > c0) then  ! ice_growth(ij) = congel(i,j)/dt
               PVt(ij) = -min(abs(PV_scale_growth*(MJ1 + MJ2*ice_growth(ij) &
                                                       - MJ3*ice_growth(ij)**2)), &
                              PV_frac_max*sk_l/dt)  
            else                           ! ice_growth(ij) = -meltb(i,j)/dt
               PVt(ij) =  min(abs(PV_scale_melt  *(      MJ2*ice_growth(ij) &
                                                       - MJ3*ice_growth(ij)**2)), &
                              PV_frac_max*sk_l/dt)
            endif
   
            if (ice_growth(ij) < c0) then ! flux from ice to ocean
               ! Algae melt like nutrients
               f_meltn(ij) = PVt(ij)*cinit(ij,nlt_bgc_N)/sk_l  ! for algae only
            elseif (ice_growth(ij) > c0 .AND. &
                   cinit(ij,nlt_bgc_N) < ocean_bio(i,j,nlt_bgc_N)*sk_l/phi_sk) then
               ! Growth only contributes to seeding from ocean 
               congel_alg(ij,nlt_bgc_N) = (ocean_bio(i,j,nlt_bgc_N)*sk_l/phi_sk &
                                        - cinit(ij,nlt_bgc_N))/dt
            endif ! PVt > c0       
         enddo    ! ij

      !----------------------------------------------------------------------
      ! 'constant':
      ! 1. Constant PV for congel > melt
      ! 2. For algae, PV for ice growth only provides a seeding concentration 
      ! 3. Melt loss (f_meltn) affects algae only and is proportional to melt
      !-----------------------------------------------------------------------

      else   ! bgc_flux_type = 'constant'

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
         do ij = 1, icells
            i = indxi(ij)
            j = indxj(ij)

            if (ice_growth(ij) > c0) PVt(ij) = -PVc
            if (ice_growth(ij) >= c0 .AND. &
                   cinit(ij,nlt_bgc_N)/sk_l < ocean_bio(i,j,nlt_bgc_N)/phi_sk) then
               congel_alg(ij,nlt_bgc_N) = (ocean_bio(i,j,nlt_bgc_N)*sk_l/phi_sk &
                                        - cinit(ij,nlt_bgc_N))/dt
            elseif (ice_growth(ij) < c0) then
               f_meltn(ij) = min(c1, meltb(i,j)/sk_l)*cinit(ij,nlt_bgc_N)/dt
           endif
         enddo  ! ij

      endif  ! bgc_flux_type

      !-----------------------------------------------------------------------
      ! begin building biogeochemistry terms
      !-----------------------------------------------------------------------

      call algal_dyn (nx_block,        ny_block,        &
                      icells,          dt,              &
                      indxi,           indxj,           &
                      fswthru,         react,           & 
                      cinit,           nbtrcr,          &
                      grow_Cn,                          &
                      tr_bgc_N_sk,     tr_bgc_Nit_sk,   &
                      tr_bgc_Am_sk,    tr_bgc_Sil_sk,   &
                      tr_bgc_C_sk,     tr_bgc_chl_sk,   &
                      tr_bgc_DMSPp_sk, tr_bgc_DMSPd_sk, & 
                      tr_bgc_DMS_sk)

      !-----------------------------------------------------------------------
      ! compute new tracer concencentrations
      !-----------------------------------------------------------------------
  
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
        
         do nn = 1, nbtrcr

      !-----------------------------------------------------------------------
      ! if PVt(ij) > 0, ie melt, then ocean_bio term drops out (MJ2006)
      ! Combine boundary fluxes
      !-----------------------------------------------------------------------
           
         PVflag(nn) = SIGN(PVflag(nn),PVt(ij))
         cinit_tmp = max(c0, cinit(ij,nn) + react(ij,nn)*sk_l)
         flux_bio_temp(nn) = (PVflag(nn)*PVt(ij)*cinit_tmp/sk_l &
                           -  PVflag(nn)*min(c0,PVt(ij))*ocean_bio(i,j,nn)) &
                           + f_meltn(ij)*cling(nn) - congel_alg(ij,nn)

         if (cinit_tmp < flux_bio_temp(nn)*dt) then
             flux_bio_temp(nn) = cinit_tmp/dt*(c1-puny)
         endif

         cinit(ij,nn) = cinit_tmp - flux_bio_temp(nn)*dt
         flux_bio(i,j,nn) = flux_bio(i,j,nn) + flux_bio_temp(nn)*phi_sk  

         ! Uncomment to update ocean concentration
         ! Currently not coupled with ocean biogeochemistry
         !ocean_bio(i,j,nn) = ocean_bio(i,j,nn) &
         !                  + flux_bio(i,j,nn)/hmix(i,j)*aicen(i,j)

         if (cinit(ij,nn) < c0) then
              write(nu_diag,*) 'sk_bgc < 0 after algal fluxes, ij,nn,cinit,flux_bio',&
                               ij,nn,cinit(ij,nn),flux_bio(i,j,nn)
              write(nu_diag,*) 'cinit_tmp,flux_bio_temp,f_meltn,congel_alg,PVt,PVflag: '
              write(nu_diag,*) cinit_tmp,flux_bio_temp(nn),f_meltn(ij), &
                               congel_alg(ij,nn),PVt(ij),PVflag(nn)
              write(nu_diag,*) 'congel, meltb: ',congel(i,j),meltb(i,j)
              call abort_ice ('ice_bgc.F90: BGC error3')
         endif
         
         enddo  ! nbtrcr

      !-----------------------------------------------------------------------
      ! reload tracer array
      !-----------------------------------------------------------------------

                              trcrn(i,j,nt_bgc_N_sk)     = cinit(ij,nlt_bgc_N)
         if (tr_bgc_Nit_sk)   trcrn(i,j,nt_bgc_Nit_sk)   = cinit(ij,nlt_bgc_NO)
         if (tr_bgc_Am_sk)    trcrn(i,j,nt_bgc_Am_sk)    = cinit(ij,nlt_bgc_NH)
         if (tr_bgc_Sil_sk)   trcrn(i,j,nt_bgc_Sil_sk)   = cinit(ij,nlt_bgc_Sil)
         if (tr_bgc_C_sk)     trcrn(i,j,nt_bgc_C_sk)     = trcrn(i,j,nt_bgc_N_sk)*R_C2N
         if (tr_bgc_chl_sk)   trcrn(i,j,nt_bgc_chl_sk)   = trcrn(i,j,nt_bgc_N_sk)*R_chl2N
         if (tr_bgc_DMSPp_sk) trcrn(i,j,nt_bgc_DMSPp_sk) = cinit(ij,nlt_bgc_DMSPp)
         if (tr_bgc_DMSPd_sk) trcrn(i,j,nt_bgc_DMSPd_sk) = cinit(ij,nlt_bgc_DMSPd)
         if (tr_bgc_DMS_sk)   trcrn(i,j,nt_bgc_DMS_sk)   = cinit(ij,nlt_bgc_DMS)
        
      enddo ! icells

      end subroutine skl_biogeochemistry

!=======================================================================
!
! Do biogeochemistry from subroutine algal_dynamics
! authors: Scott Elliott, LANL
!          Nicole Jeffery, LANL

      subroutine algal_dyn (nx_block,     ny_block,     &
                            icells,       dt,           &
                            indxi,        indxj,        &
                            fswthru,      reactb,       & 
                            ltrcrn,       nbtrcr,       &
                            growN,                      &
                            tr_bio_N,     tr_bio_NO,    &
                            tr_bio_NH,    tr_bio_Sil,   &
                            tr_bio_C,     tr_bio_chl,   &
                            tr_bio_DMSPp, tr_bio_DMSPd, &
                            tr_bio_DMS)

      use ice_constants, only: p1, p5, c0, c1, secday

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         icells            , & ! number of cells with aicen > puny
         nbtrcr                ! number of layer tracers

      integer (kind=int_kind), dimension(nx_block*ny_block), &
         intent(in) :: &
         indxi, indxj ! compressed indices for cells with aicen > puny

      real (kind=dbl_kind), intent(in) :: &
         dt         ! time step 

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(in) :: &
         fswthru    ! average shortwave passing through current ice layer (W/m^2)

      real (kind=dbl_kind), dimension (nx_block,ny_block), intent(inout) :: &
         growN      !  algal specific growth rate   (1/s)

      real (kind=dbl_kind), dimension(icells,nbtrcr), intent(inout) :: &
         reactb     ! biological reaction terms (mmol/m^3)

      real (kind=dbl_kind), dimension(icells,nbtrcr), intent(in) :: &
         ltrcrn     ! concentrations in layer * sk_l

      ! tracer flags for vertical or skeletal layer bgc
      logical (kind=log_kind), intent(in):: & 
         tr_bio_N    , & ! algal nitrogen
         tr_bio_NO   , & ! algal nitrate  
         tr_bio_NH   , & ! algal ammonium
         tr_bio_Sil  , & ! algal silicate
         tr_bio_C    , & ! algal carbon
         tr_bio_chl  , & ! algal chlorophyll
         tr_bio_DMSPp, & ! DMSP particulate
         tr_bio_DMSPd, & ! DMSP dissolved
         tr_bio_DMS      ! DMS

      !  local variables

      real (kind=dbl_kind), parameter :: & 
         T_bot      = -1.8_dbl_kind , & ! interface close to freezing (C)
         chlabs     = 9.e-4_dbl_kind, & ! chlorophyll absorption (1/(mg/m^3))
                                        ! assuming skeletal layer thickness = 3 cm
         mu_max     = 1.5_dbl_kind,   & ! maximum growth rate 'Jin2006' (1/day) 
         T_max      = -1.8_dbl_kind,  & ! maximum growth at Tmax (C)
         op_dep_min = 0.1          ,  & ! optical depth above which light attentuates
         grow_Tdep  = 0.0633_dbl_kind,& ! and its T dependence (1/C)
         fr_graze   = p1           ,  & ! A93 val for S, set to zero in Jin 06
         fr_graze_s = 0.5_dbl_kind ,  & ! fraction of grazing spilled or slopped
         fr_graze_a = 0.5_dbl_kind ,  & ! fraction of grazing assimilated
         fr_graze_e = 0.5_dbl_kind ,  & ! fraction of assimilation excreted  
         alpha2max  = 0.8_dbl_kind,   & ! light limitation (1/(W/m^2))
         !beta2max   = 0.018_dbl_kind, & ! corresponding light inhibition (1/W/m^2)
         K_Nit      = c1           ,  & ! nitrate half saturation (mmol/m^3) 
         K_Am       = c1           ,  & ! ammonium half saturation (mmol/m^3) 
         K_Sil      = 4.0_dbl_kind ,  & ! silicon half saturation (mmol/m^3)
         mort_pre   = 0.0208_dbl_kind,& ! prefix to mortality (1/day) 
         mort_Tdep  = 0.03_dbl_kind , & ! T dependence of mortality (1/C) 
         fr_mort2min= c1            , & ! fractionation to remin
         !t_nitrif   = 67.0_dbl_kind , & ! nitrification time scale (days)
         max_loss   = 0.9               ! restrict uptake to 90% of remaining value 

      real (kind=dbl_kind), parameter :: &
         fr_excrt_2S= c1           , &  ! excretion is efficient initially
         y_sk_DMS   = c1           , &  ! and conversion given high yield
         t_sk_conv  = 10.0_dbl_kind, &  ! at a Stefels rate (days)
         t_sk_ox    = 10.0_dbl_kind     ! DMS in turn oxidizes slowly (days)

     ! real (kind=dbl_kind), parameter :: &
     !    pr_l       = 10.0_dbl_kind, & ! product layer thickness (m) 
     !    chl_pr_v   = 0.1_dbl_kind , & ! fixed nondiatom chlorophyll in ml (mg/m^3)
     !    R_chl2N_nd = 3.0_dbl_kind , & ! shade adaptation below (mg/millimole)
     !    R_C2N_nd   = 7.0_dbl_kind , & ! open water ratio (mole/mole)
     !    t_pr_dsrp  = 10.0_dbl_kind    ! disruption time scale (days)
         
     ! real (kind=dbl_kind), parameter :: &
     !    R_S2N_nd   = 0.03_dbl_kind, & ! open water ratio nondiatoms (mole/mole)
     !    y_pr_DMS   = c1           , & ! but we begin again with unit yield
     !    t_pr_conv  = 10.0_dbl_kind, & ! and a similar conversion (days)
     !    t_pr_ox    = 10.0_dbl_kind    ! plus round final time scale (days)

      integer (kind=int_kind) :: i, j, ij

      real (kind=dbl_kind) :: &
         Nin        , &     ! algal nitrogen concentration on volume (mmol/m^3) 
         Cin        , &     ! algal carbon concentration on volume (mmol/m^3)
         chlin      , &     ! algal chlorophyll concentration on volume (mg/m^3)
         NOin       , &     ! nitrate concentration on volume (mmol/m^3) 
         NHin       , &     ! ammonia/um concentration on volume (mmol/m^3) 
         Silin      , &     ! silicon concentration on volume (mmol/m^3) 
         DMSPpin    , &     ! DMSPp concentration on volume (mmol/m^3)
         DMSPdin    , &     ! DMSPd concentration on volume (mmol/m^3)
         DMSin              ! DMS concentration on volume (mmol/m^3)

      real (kind=dbl_kind) :: &
         op_dep     , &  ! bottom layer attenuation exponent (optical depth)
         Iavg_loc        ! bottom layer attenuated Fswthru (W/m^2)

      real (kind=dbl_kind) :: &
         L_lim    , &  ! overall light limitation 
         Nit_lim  , &  ! overall nitrate limitation
         Am_lim   , &  ! overall ammonium limitation
         N_lim    , &  ! overall nitrogen species limitation
         Sil_lim  , &  ! overall silicon limitation
         fr_Nit   , &  ! fraction of local ecological growth as nitrate
         fr_Am    , &  ! fraction of local ecological growth as ammonia
         growmax_N, &  ! maximum growth rate in N currency (mmol/m^3/s)
         grow_N   , &  ! true growth rate in N currency (mmol/m^3/s)
         potU_Nit , &  ! potential nitrate uptake (mmol/m^3/s)
         potU_Am  , &  ! potential ammonium uptake (mmol/m^3/s)
         U_Nit    , &  ! actual nitrate uptake (mmol/m^3/s)
         U_Am     , &  ! actual ammonium uptake (mmol/m^3/s)
         U_Sil         ! actual silicon uptake (mmol/m^3/s)

      real (kind=dbl_kind) :: &
         resp     , &  ! respiration (mmol/m^3/s)
         graze    , &  ! grazing (mmol/m^3/s)
         mort     , &  ! sum of mortality and excretion (mmol/m^3/s)
         nitrif        ! nitrification (mmol/m^3/s)

      ! source terms underscore s, removal underscore r

      real (kind=dbl_kind) :: &
         N_s_p     , &  ! algal nitrogen photosynthesis (mmol/m^3)
         N_r_g     , &  ! algal nitrogen losses to grazing (mmol/m^3)
         N_r_r     , &  ! algal nitrogen losses to respiration (mmol/m^3)
         N_r_mo    , &  ! algal nitrogen losses to mortality (mmol/m^3)
         N_s       , &  ! net algal nitrogen sources (mmol/m^3)
         N_r       , &  ! net algal nitrogen removal (mmol/m^3)
         C_s       , &  ! net algal carbon sources (mmol/m^3)
         C_r       , &  ! net algal carbon removal (mmol/m^3)
         NO_s_n    , &  ! skl nitrate from nitrification (mmol/m^3)
         NO_r_p    , &  ! skl nitrate uptake by algae (mmol/m^3)
         NO_s      , &  ! net skl nitrate sources (mmol/m^3)
         NO_r      , &  ! net skl nitrate removal (mmol/m^3)
         NH_s_e    , &  ! skl ammonium source from excretion (mmol/m^3)
         NH_s_r    , &  ! skl ammonium source from respiration (mmol/m^3)
         NH_s_mo   , &  ! skl ammonium source from mort/remin (mmol/m^3) 
         NH_r_p    , &  ! skl ammonium uptake by algae (mmol/m^3)
         NH_r_n    , &  ! skl ammonium removal to nitrification (mmol/m^3)
         NH_s      , &  ! net skl ammonium sources (mmol/m^3)
         NH_r      , &  ! net skl ammonium removal (mmol/m^3)
         Sil_r_p   , &  ! skl silicon uptake by algae (mmol/m^3)
         Sil_s     , &  ! net skl silicon sources (mmol/m^3)
         Sil_r          ! net skl silicon removal (mmol/m^3)

      real (kind=dbl_kind) :: &
         DMSPd_s_s , &  ! skl dissolved DMSP from grazing spillage (mmol/m^3)
         DMSPd_s_e , &  ! skl dissolved DMSP from zooplankton excretion (mmol/m^3)
         DMSPd_s_mo, &  ! skl dissolved DMSP from MBJ algal mortexc (mmol/m^3)
         DMSPd_r_c , &  ! skl dissolved DMSP conversion (mmol/m^3)
         DMSPd_s   , &  ! net skl dissolved DMSP sources (mmol/m^3)
         DMSPd_r   , &  ! net skl dissolved DMSP removal (mmol/m^3)
         DMS_s_c   , &  ! skl DMS source via conversion (mmol/m^3)
         DMS_r_o   , &  ! skl DMS losses due to oxidation (mmol/m^3)
         DMS_s     , &  ! net skl DMS sources (mmol/m^3)
         DMS_r          ! net skl DMS removal (mmol/m^3)

     ! real (kind=dbl_kind) :: &
     !    DMSP_pr_s_nd , &  ! product layer dissolved DMSP from local bio (mmol/m^3)
     !    DMSP_pr_s_me , &  ! product layer dissolved DMSP from melting (mmol/m^3)
     !    DMSP_pr_r_c  , &  ! product layer dissolved DMSP conversion (mmol/m^3)
     !    DMSP_pr_s    , &  ! net product dissolved DMSP sources (mmol/m^3)
     !    DMSP_pr_r    , &  ! net product dissolved DMSP removal (mmol/m^3)
     !    DMS_pr_s_c   , &  ! product layer DMS source via conversion (mmol/m^3)
     !    DMS_pr_r_o   , &  ! product layer DMS losses due to oxidation (mmol/m^3)
     !    DMS_pr_s     , &  ! net product DMS sources (mmol/m^3)
     !    DMS_pr_r          ! net product DMS removal (mmol/m^3)

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)

      !-----------------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------------

         Nin     = c0
         Cin     = c0
         chlin   = c0
         NOin    = c0
         NHin    = c0
         Silin   = c0
         DMSPpin = c0
         DMSPdin = c0
         DMSin   = c0

                               NOin     = ltrcrn(ij,nlt_bgc_NO)   /sk_l
         if (tr_bio_N)         Nin      = ltrcrn(ij,nlt_bgc_N)    /sk_l
         if (tr_bio_C)         Cin      = ltrcrn(ij,nlt_bgc_C)    /sk_l
         if (tr_bio_NH)        NHin     = ltrcrn(ij,nlt_bgc_NH)   /sk_l
         if (tr_bio_Sil)       Silin    = ltrcrn(ij,nlt_bgc_Sil)  /sk_l
         if (tr_bio_DMSPp)     DMSPpin  = ltrcrn(ij,nlt_bgc_DMSPp)/sk_l
         if (tr_bio_DMSPd)     DMSPdin  = ltrcrn(ij,nlt_bgc_DMSPd)/sk_l
         if (tr_bio_DMS)       DMSin    = ltrcrn(ij,nlt_bgc_DMS)  /sk_l
         chlin = R_chl2N * Nin  

      !-----------------------------------------------------------------------
      ! Light limitation
      !-----------------------------------------------------------------------

         op_dep   = chlabs * chlin   

         ! Okhotsk maxima causes a reevaluation.
         ! The new concept is, late algae at the bottom of the bottom strongly 
         ! attenuated.  Since k about 0.03 1/m(mg/m3), efold at about 30 mg/m3
         ! values of order hundreds will be shut down...
         ! Newest algae experience super aborption above because they sit low.
         ! More than perhaps two efolds and light falls below half value.

         if (op_dep > op_dep_min) then
            Iavg_loc = fswthru(i,j) * (c1 - exp(-op_dep)) / op_dep
         else
            Iavg_loc = fswthru(i,j)
         endif

         ! With light inhibition
         ! L_lim = (c1 - exp(-alpha2max*Iavg_loc)) * exp(-beta2max*Iavg_loc)      

         ! Without light inhibition
         L_lim = (c1 - exp(-alpha2max*Iavg_loc)) 

      !-----------------------------------------------------------------------
      ! Nutrient limitation
      !-----------------------------------------------------------------------

         Nit_lim = NOin/(NOin + K_Nit)
         Am_lim  = c0
         N_lim = Nit_lim
         if (tr_bio_NH) then
            Am_lim = NHin/(NHin + K_Am)
            N_lim  = min(c1, Nit_lim + Am_lim)  
         endif
         Sil_lim = c1
         if (tr_bio_Sil) Sil_lim = Silin/(Silin + K_Sil)
 
         ! Growth and uptake computed within the bottom layer 
         ! Note here per A93 discussions and MBJ model, salinity is a universal 
         ! restriction.  Comparison with available column nutrients inserted 
         ! but in tests had no effect.
         ! Primary production reverts to SE form, see MBJ below and be careful

         growmax_N = mu_max / secday * exp(grow_Tdep * (T_bot - T_max))* Nin 
         grow_N    = min(L_lim, N_lim, Sil_lim) * growmax_N
         potU_Nit  = Nit_lim                    * growmax_N
         potU_Am   = Am_lim                     * growmax_N 
         U_Am      = min(grow_N, potU_Am)
         U_Nit     = grow_N - U_Am
         U_Sil     = R_Si2N * grow_N

         if (tr_bio_Sil) U_Sil = min(U_Sil, max_loss * Silin/dt)
         U_Nit = min(U_Nit, max_loss * NOin/dt)  
         U_Am  = min(U_Am,  max_loss * NHin/dt)    
 
         grow_N = min(U_Sil/R_Si2N,U_Nit + U_Am)
         fr_Am = c0
         if (tr_bio_NH) then
            fr_Am = p5
            if (grow_N > c0) fr_Am = min(U_Am/grow_N, c1)
         endif
         fr_Nit = c1 - fr_Am
         U_Nit  = fr_Nit * grow_N
         U_Am   = fr_Am  * grow_N
         U_Sil  = R_Si2N * grow_N
    
         resp   = fr_resp  * grow_N
         graze  = fr_graze * grow_N
         mort   = mort_pre * exp(mort_Tdep*(T_bot-T_max)) * Nin / secday 
         nitrif = c0 ! (NHin / t_nitrif) / secday
 
         ! history variables

         growN(i,j) = grow_N 
         if (Nin > c0) growN(i,j) = grow_N/Nin  ! specific growth rate (per s)

      !-----------------------------------------------------------------------
      ! Define reaction terms
      !-----------------------------------------------------------------------

         ! Since the framework remains incomplete at this point,
         ! it is assumed as a starting expedient that 
         ! DMSP loss to melting results in 10% conversion to DMS
         ! which is then given a ten day removal constant.
         ! Grazing losses are channeled into rough spillage and assimilation
         ! then following ammonia there is some recycling.

         !--------------------------------------------------------------------
         ! Algal reaction term
         ! N_react = (grow_N*(c1 -fr_resp - fr_graze) - mort)*dt
         !--------------------------------------------------------------------

         N_s_p  = grow_N * dt  
         N_r_g  = graze  * dt 
         N_r_r  = resp   * dt
         N_r_mo = mort   * dt
         N_s    = N_s_p
         N_r    = N_r_g + N_r_r + N_r_mo 

         !--------------------------------------------------------------------
         ! Carbon chemistry 
         ! C_react = R_C2N * N_react 
         !--------------------------------------------------------------------

         C_s = R_C2N * N_s
         C_r = R_C2N * N_r

         !--------------------------------------------------------------------
         ! Nitrate reaction term
         ! NO_react = (nitrif - fr_Nit*grow_N)*dt
         !--------------------------------------------------------------------

         NO_s_n = nitrif * dt               
         NO_r_p = U_Nit  * dt                
         NO_s   = NO_s_n
         NO_r   = NO_r_p

         !--------------------------------------------------------------------
         ! Ammonium reaction term
         ! NH_react = (-nitrif - (c1-fr_Nit - fr_resp 
         !     - fr_graze*fr_graze_e*fr_graze_a)*grow_N + mort*fr_mort2min)*dt  
         !--------------------------------------------------------------------

         NH_s_r  = N_r_r 
         NH_s_e  = fr_graze_e * fr_graze_a * N_r_g
         NH_s_mo = fr_mort2min * N_r_mo
         NH_r_p  = U_Am   * dt
         NH_r_n  = nitrif * dt
         NH_s    = NH_s_r + NH_s_e + NH_s_mo
         NH_r    = NH_r_p + NH_r_n 

         !--------------------------------------------------------------------
         ! Silica reaction term
         ! Sil_react = - R_Si2N * grow_N * dt
         !--------------------------------------------------------------------
     
         Sil_r_p = U_Sil * dt
         Sil_s   = c0
         Sil_r   = Sil_r_p 

         !--------------------------------------------------------------------
         ! Sulfur cycle begins here
         !--------------------------------------------------------------------
         ! Grazing losses are channeled into rough spillage and assimilation
         ! then onward and the MBJ mortality channel is included
         ! It is assumed as a starting expedient that 
         ! DMSP loss to melting gives partial conversion to DMS in product layer
         ! which then undergoes Stefels removal.

         !--------------------------------------------------------------------
         ! DMSPd  reaction term
         ! DMSPd_react = R_S2N*((fr_graze_s+fr_excrt_2S*fr_graze_e*fr_graze_a)
         !                      *fr_graze*grow_N + fr_mort2min*mort)*dt
         !             - [\DMSPd]/t_sk_conv*dt
         !--------------------------------------------------------------------

         DMSPd_s_s  = fr_graze_s * R_S2N * N_r_g
         DMSPd_s_e  = fr_excrt_2S * fr_graze_e * fr_graze_a * R_S2N * N_r_g
         DMSPd_s_mo = fr_mort2min * R_S2N * N_r_mo
         DMSPd_r_c  = DMSPdin * dt / (t_sk_conv * secday) 
         DMSPd_s    = DMSPd_s_s + DMSPd_s_e + DMSPd_s_mo 
         DMSPd_r    = DMSPd_r_c

         !--------------------------------------------------------------------
         ! DMS reaction term
         ! DMS_react = ([\DMSPd]*y_sk_DMS/t_sk_conv - c1/t_sk_ox *[\DMS])*dt
         !--------------------------------------------------------------------
 
         DMS_s_c = y_sk_DMS * DMSPd_r_c 
         DMS_r_o = DMSin * dt / (t_sk_ox * secday)
         DMS_s   = DMS_s_c
         DMS_r   = DMS_r_o

         ! for mixed layer sulfur chemistry, fluxes kept separate for ice 
         ! area weighting
         ! units are different here, but not sure if they need to be changed
         ! no fluxes into the product layer here

         ! DMSP_pr_s_nd = chl_pr_v*pr_l * R_S2N_nd/R_chl2N_nd * dt &
         !              / (t_pr_dsrp * secday)
         ! DMSP_pr_s_me = fr_melt_2S * DMSPp_sk_r_me
         ! DMSP_pr_r_c  = dmsp(i,j,iblk)*pr_l * dt / (t_pr_conv * secday) 
         ! DMSP_pr_f    = F_DMSP * dt
         ! DMSP_pr_s    = DMSP_pr_s_nd                ! + DMSP_pr_s_me + DMSP_pr_f
         ! DMSP_pr_r    = DMSP_pr_r_c

         ! DMS_pr_s_c = y_pr_DMS * DMSP_pr_r_c
         ! DMS_pr_r_o = dms(i,j,iblk)*pr_l * dt / (t_pr_ox * secday)
         ! DMS_pr_f   = F_DMS * dt 
         ! DMS_pr_s   = DMS_pr_s_c                   ! + DMS_pr_f
         ! DMS_pr_r   = DMS_pr_r_o

      !-----------------------------------------------------------------------
      ! Load reaction array
      !-----------------------------------------------------------------------

                           reactb(ij,nlt_bgc_NO)    = NO_s    - NO_r
         if (tr_bio_N)     reactb(ij,nlt_bgc_N)     = N_s     - N_r 
         if (tr_bio_C)     reactb(ij,nlt_bgc_C)     = C_s     - C_r
         if (tr_bio_NH)    reactb(ij,nlt_bgc_NH)    = NH_s    - NH_r
         if (tr_bio_Sil)   reactb(ij,nlt_bgc_Sil)   = Sil_s   - Sil_r
         if (tr_bio_DMSPd) reactb(ij,nlt_bgc_DMSPd) = DMSPd_s - DMSPd_r
         if (tr_bio_DMS)   reactb(ij,nlt_bgc_DMS)   = DMS_s   - DMS_r

      enddo

      end subroutine algal_dyn

!=======================================================================
!
! Writes diagnostic info (max, min, global sums, etc) to standard out
!
! authors: Elizabeth C. Hunke, LANL
!          Nicole Jeffery, LANL

      subroutine bgc_diags (dt)

      use ice_broadcast, only: broadcast_scalar
      use ice_constants, only: c0, mps_to_cmpdy, c100
      use ice_diagnostics, only: npnt, print_points, pmloc, piloc, pjloc, pbloc
      use ice_timers, only: timer_bgc, ice_timer_start, ice_timer_stop

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      ! local variables

      integer (kind=int_kind) :: &
         i, j, n, iblk

      ! fields at diagnostic points
      real (kind=dbl_kind), dimension(npnt) :: &
         pN_sk, pNit_sk, pAm_sk, pSil_sk, &
         pDMSPp_sk, pDMSPd_sk, pDMS_sk, pN_ac, &
         pNit_ac, pAm_ac, pSil_ac, pDMSP_ac, pDMS_ac, &
         pflux_NO, pflux_N, pflux_Sil, pflux_NH

      call ice_timer_start(timer_bgc) ! biogeochemistry

      if (print_points) then

      !-----------------------------------------------------------------
      ! biogeochemical 
      ! state of the ice and associated fluxes for 2 defined points
      ! NOTE these are computed for the last timestep only (not avg)
      !-----------------------------------------------------------------

         do n = 1, npnt
            if (my_task == pmloc(n)) then
               i = piloc(n)
               j = pjloc(n)
               iblk = pbloc(n)
               pAm_ac(n)   = c0
               pSil_ac(n)  = c0
               pDMSP_ac(n) = c0
               pDMS_ac(n)  = c0

               pN_ac(n)    = ocean_bio(i,j,nlt_bgc_N,iblk)    ! algalN(i,j,iblk)
               pNit_ac(n)  = ocean_bio(i,j,nlt_bgc_NO,iblk)   ! nit(i,j,iblk)
               if (tr_bgc_Am_sk) &
               pAm_ac(n)   = ocean_bio(i,j,nlt_bgc_NH,iblk)   ! amm(i,j,iblk)
               if (tr_bgc_Sil_sk) &
               pSil_ac(n)  = ocean_bio(i,j,nlt_bgc_Sil,iblk)  ! sil(i,j,iblk)
               if (tr_bgc_DMS_sk) then
               pDMSP_ac(n) = ocean_bio(i,j,nlt_bgc_DMSPp,iblk)! dmsp(i,j,iblk)
               pDMS_ac(n)  = ocean_bio(i,j,nlt_bgc_DMS,iblk)  ! dms(i,j,iblk)
               endif

               ! fluxes in mmol/m^2/d
               ! concentrations are bulk in mmol/m^3

               if (skl_bgc) then
                  pNit_sk(n)   = c0
                  pAm_sk(n)    = c0
                  pSil_sk(n)   = c0
                  pDMSPp_sk(n) = c0
                  pDMSPd_sk(n) = c0
                  pDMS_sk(n)   = c0
               
                  pN_sk(n)       = trcr    (i,j, nt_bgc_N_sk,   iblk)*phi_sk/sk_l
                  pflux_N(n)     = flux_bio(i,j,nlt_bgc_N,      iblk)*mps_to_cmpdy/c100 
                  if (tr_bgc_Nit_sk) then
                     pNit_sk(n)  = trcr    (i,j, nt_bgc_Nit_sk, iblk)*phi_sk/sk_l   
                     pflux_NO(n) = flux_bio(i,j,nlt_bgc_NO,     iblk)*mps_to_cmpdy/c100 
                  endif
                  if (tr_bgc_Am_sk) then
                     pAm_sk(n)   = trcr    (i,j, nt_bgc_Am_sk,  iblk)*phi_sk/sk_l
                     pflux_NH(n) = flux_bio(i,j,nlt_bgc_NH,     iblk)*mps_to_cmpdy/c100 
                  endif
                  if (tr_bgc_Sil_sk) then
                     pSil_sk(n)  = trcr    (i,j, nt_bgc_Sil_sk, iblk)*phi_sk/sk_l       
                     pflux_Sil(n)= flux_bio(i,j,nlt_bgc_Sil,    iblk)*mps_to_cmpdy/c100 
                  endif
                  if (tr_bgc_DMSPp_sk) &
                     pDMSPp_sk(n) = trcr   (i,j,nt_bgc_DMSPp_sk,iblk)*phi_sk/sk_l  
                  if (tr_bgc_DMSPd_sk) &
                     pDMSPd_sk(n) = trcr   (i,j,nt_bgc_DMSPd_sk,iblk)*phi_sk/sk_l 
                  if (tr_bgc_DMS_sk)   &
                     pDMS_sk  (n) = trcr   (i,j,nt_bgc_DMS_sk,  iblk)*phi_sk/sk_l  
               endif

            endif                 ! my_task = pmloc
 
            call broadcast_scalar(pN_ac    (n), pmloc(n))     
            call broadcast_scalar(pNit_ac  (n), pmloc(n))             
            call broadcast_scalar(pAm_ac   (n), pmloc(n))             
            call broadcast_scalar(pSil_ac  (n), pmloc(n))             
            call broadcast_scalar(pDMSP_ac (n), pmloc(n))             
            call broadcast_scalar(pDMS_ac  (n), pmloc(n))
            call broadcast_scalar(pflux_N  (n), pmloc(n))     
            call broadcast_scalar(pflux_NO (n), pmloc(n))             
            call broadcast_scalar(pflux_NH (n), pmloc(n))             
            call broadcast_scalar(pflux_Sil(n), pmloc(n))

            if (skl_bgc) then   ! skl_bgc
               call broadcast_scalar(pN_sk    (n), pmloc(n))            
               call broadcast_scalar(pNit_sk  (n), pmloc(n))             
               call broadcast_scalar(pAm_sk   (n), pmloc(n))             
               call broadcast_scalar(pSil_sk  (n), pmloc(n))             
               call broadcast_scalar(pDMSPp_sk(n), pmloc(n))             
               call broadcast_scalar(pDMSPd_sk(n), pmloc(n))             
               call broadcast_scalar(pDMS_sk  (n), pmloc(n))        
           endif

         enddo                  ! npnt
      endif                     ! print_points

      !-----------------------------------------------------------------
      ! start spewing
      !-----------------------------------------------------------------

      if (my_task == master_task) then

      call flush_fileunit(nu_diag)

      !-----------------------------------------------------------------
      ! diagnostics for Arctic and Antarctic points
      !-----------------------------------------------------------------

      if (print_points) then
      if (skl_bgc) then  

      write(nu_diag,*) '----------BGC----------'

      write(nu_diag,*) '------bulk skl bgc-----'
      write(nu_diag,900) 'nitrogen    (mmol/m^3) = ',pN_sk(1),pN_sk(2)
      write(nu_diag,900) 'nitrate     (mmol/m^3) = ',pNit_sk(1),pNit_sk(2)
      if (tr_bgc_Am_sk) &
      write(nu_diag,900) 'ammonia/um  (mmol/m^3) = ',pAm_sk(1),pAm_sk(2)
      if (tr_bgc_Sil_sk) &
      write(nu_diag,900) 'silicon     (mmol/m^3) = ',pSil_sk(1),pSil_sk(2)
      if (tr_bgc_DMS_sk) then
      write(nu_diag,900) 'DMSPp       (mmol/m^3) = ',pDMSPp_sk(1),pDMSPp_sk(2)
      write(nu_diag,900) 'DMSPd       (mmol/m^3) = ',pDMSPd_sk(1),pDMSPd_sk(2)
      write(nu_diag,900) 'DMS         (mmol/m^3) = ',pDMS_sk(1),pDMS_sk(2)
      endif

      write(nu_diag,*) '---ice-ocean fluxes----'
      write(nu_diag,900) 'algalN flx(mmol/m^2/d) = ',pflux_N(1),pflux_N(2)
      write(nu_diag,900) 'nit. flux (mmol/m^2/d) = ',pflux_NO(1),pflux_NO(2)
      if (tr_bgc_Am_sk) &
      write(nu_diag,900) 'amm. flux (mmol/m^2/d) = ',pflux_NH(1),pflux_NH(2)
      if (tr_bgc_Sil_sk) &
      write(nu_diag,900) 'sil. flux (mmol/m^2/d) = ',pflux_Sil(1),pflux_Sil(2)

      write(nu_diag,*) '---ocean mixed layer---'
      write(nu_diag,900) 'algal N     (mmol/m^3) = ',pN_ac(1),pN_ac(2)
      write(nu_diag,900) 'nitrate     (mmol/m^3) = ',pNit_ac(1),pNit_ac(2)
      if (tr_bgc_Am_sk) &
      write(nu_diag,900) 'ammonia/um  (mmol/m^3) = ',pAm_ac(1),pAm_ac(2)
      if (tr_bgc_Sil_sk) &
      write(nu_diag,900) 'silicon     (mmol/m^3) = ',pSil_ac(1),pSil_ac(2)
      if (tr_bgc_DMS_sk) then
      write(nu_diag,900) 'DMSP        (mmol/m^3) = ',pDMSP_ac(1),pDMSP_ac(2)
      write(nu_diag,900) 'DMS         (mmol/m^3) = ',pDMS_ac(1),pDMS_ac(2)
      endif

      endif                   ! skl_bgc
      endif                   ! print_points

      endif                   ! my_task = master_task 

      call ice_timer_stop(timer_bgc) ! biogeochemistry

  900 format (a25,2x,f24.17,2x,f24.17)

      end subroutine bgc_diags

!=======================================================================
!---! these subroutines write/read Fortran unformatted data files ..
!=======================================================================
!
! Dumps all values needed for a bgc restart
!
! author Elizabeth C. Hunke, LANL

      subroutine write_restart_bgc()

      use ice_domain_size, only: ncat
      use ice_state, only: trcrn
      use ice_restart,only: write_restart_field

      ! local variables

      logical (kind=log_kind) :: diag

      diag = .true.

      !-----------------------------------------------------------------
      ! Skeletal layer BGC
      !-----------------------------------------------------------------

         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_N_sk,:,:), &
                                  'ruf8','bgc_N_sk',ncat,diag)
         if (tr_bgc_C_sk) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_C_sk,:,:), &
                                  'ruf8','bgc_C_sk',ncat,diag)
         if (tr_bgc_chl_sk) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_chl_sk,:,:), &
                                  'ruf8','bgc_chl_sk',ncat,diag)
         if (tr_bgc_Nit_sk) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Nit_sk,:,:), &
                                  'ruf8','bgc_Nit_sk',ncat,diag)
         if (tr_bgc_Am_sk) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Am_sk,:,:), &
                                  'ruf8','bgc_Am_sk',ncat,diag)
         if (tr_bgc_Sil_sk) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_Sil_sk,:,:), &
                                  'ruf8','bgc_Sil_sk',ncat,diag)
         if (tr_bgc_DMSPp_sk) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPp_sk,:,:), &
                                  'ruf8','bgc_DMSPp_sk',ncat,diag)
         if (tr_bgc_DMSPd_sk) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMSPd_sk,:,:), &
                                  'ruf8','bgc_DMSPd_sk',ncat,diag)
         if (tr_bgc_DMS_sk) &
         call write_restart_field(nu_dump_bgc,0,trcrn(:,:,nt_bgc_DMS_sk,:,:), &
                                  'ruf8','bgc_DMS_sk',ncat,diag)
      
      !-----------------------------------------------------------------
      ! Ocean BGC
      !-----------------------------------------------------------------

      if (tr_bgc_N_sk) &
      call write_restart_field(nu_dump_bgc,0,algalN,'ruf8','algalN',1,diag)
      if (tr_bgc_Nit_sk) &
      call write_restart_field(nu_dump_bgc,0,nit,   'ruf8','nit',   1,diag)
      if (tr_bgc_Am_sk) &
      call write_restart_field(nu_dump_bgc,0,amm,   'ruf8','amm',   1,diag)
      if (tr_bgc_Sil_sk) &
      call write_restart_field(nu_dump_bgc,0,sil,   'ruf8','sil',   1,diag)
      if (tr_bgc_DMSPp_sk) &
      call write_restart_field(nu_dump_bgc,0,dmsp,  'ruf8','dmsp',  1,diag)
      if (tr_bgc_DMS_sk) &
      call write_restart_field(nu_dump_bgc,0,dms,   'ruf8','dms',   1,diag)

      end subroutine write_restart_bgc

!=======================================================================
!
! Reads all values needed for a bgc restart
!
! author Elizabeth C. Hunke, LANL

      subroutine read_restart_bgc()

      use ice_constants, only: field_loc_center, field_type_scalar
      use ice_domain_size, only: ncat
      use ice_state, only: trcrn
      use ice_restart,only: read_restart_field

      ! local variables

      logical (kind=log_kind) :: diag

      diag = .true.

      !-----------------------------------------------------------------
      ! Skeletal Layer BGC
      !-----------------------------------------------------------------

      if (my_task == master_task) write(nu_diag,*) 'skl bgc restart'

      call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_N_sk,:,:), &
           'ruf8','bgc_N_sk',ncat,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_C_sk) &
      call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_C_sk,:,:), &
           'ruf8','bgc_C_sk',ncat,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_chl_sk) &
      call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_chl_sk,:,:), &
           'ruf8','bgc_chl_sk',ncat,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_Nit_sk) &
      call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Nit_sk,:,:), &
           'ruf8','bgc_Nit_sk',ncat,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_Am_sk) &
      call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Am_sk,:,:), &
           'ruf8','bgc_Am_sk',ncat,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_Sil_sk) &
      call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_Sil_sk,:,:), &
           'ruf8','bgc_Sil_sk',ncat,diag,field_loc_center,field_type_scalar)
      if(tr_bgc_DMSPp_sk) &
      call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMSPp_sk,:,:), &
           'ruf8','bgc_DMSPp_sk',ncat,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_DMSPd_sk) &
      call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMSPd_sk,:,:), &
           'ruf8','bgc_DMSPd_sk',ncat,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_DMS_sk) &
      call read_restart_field(nu_restart_bgc,0,trcrn(:,:,nt_bgc_DMS_sk,:,:), &
           'ruf8','bgc_DMS_sk',ncat,diag,field_loc_center,field_type_scalar)

      !-----------------------------------------------------------------
      ! Ocean BGC
      !-----------------------------------------------------------------

      if (my_task == master_task) write(nu_diag,*) 'mixed layer ocean bgc restart'

      if (tr_bgc_N_sk) &
      call read_restart_field(nu_restart_bgc,0,algalN,'ruf8','algalN',&
                              1,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_Nit_sk) &
      call read_restart_field(nu_restart_bgc,0,nit   ,'ruf8','nit'   ,&
                              1,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_Am_sk) &
      call read_restart_field(nu_restart_bgc,0,amm   ,'ruf8','amm'   ,&
                              1,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_Sil_sk) &
      call read_restart_field(nu_restart_bgc,0,sil   ,'ruf8','sil'   ,&
                              1,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_DMSPp_sk) &
      call read_restart_field(nu_restart_bgc,0,dmsp  ,'ruf8','dmsp'  ,&
                              1,diag,field_loc_center,field_type_scalar)
      if (tr_bgc_DMS_sk) &
      call read_restart_field(nu_restart_bgc,0,dms   ,'ruf8','dms'   ,&
                              1,diag,field_loc_center,field_type_scalar)

      end subroutine read_restart_bgc

!=======================================================================

      end module ice_algae

!=======================================================================
