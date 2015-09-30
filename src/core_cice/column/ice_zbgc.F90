!  SVN:$Id: ice_zbgc.F90 1012 2015-06-26 12:34:09Z eclare $
!=======================================================================
!
! Biogeochemistry driver
!
! authors: Nicole Jeffery, LANL
!          Scott Elliot,   LANL
!          Elizabeth C. Hunke, LANL
!
      module ice_zbgc

      use ice_kinds_mod
      use ice_zbgc_shared ! everything

      implicit none 

      private
      public :: add_new_ice_bgc, init_zbgc, init_bgc, &
                init_history_bgc, biogeochemistry

!=======================================================================

      contains

!=======================================================================

! Namelist variables, set to default values; may be altered at run time
! 
! author Elizabeth C. Hunke, LANL
!        Nicole Jeffery, LANL

      subroutine init_zbgc

      use ice_broadcast, only: broadcast_scalar
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c1, p5, c0, rhos, rhoi
      use ice_domain_size, only: max_ntrcr, max_nbtrcr, nblyr
      use ice_exit, only: abort_ice
      use ice_fileunits, only: nu_nml, nml_filename, get_fileunit, &
                               release_fileunit, nu_diag
      use ice_restart_shared, only: runtype
      use ice_state, only: trcr_base, trcr_depend
      use ice_colpkg_tracers, only: tr_brine, nt_fbri, ntrcr, nbtrcr, &
          nt_bgc_N_sk, nt_bgc_Nit_sk, nt_bgc_chl_sk, nt_bgc_Am_sk, &
          nt_bgc_Sil_sk, nt_bgc_DMSPp_sk, nt_bgc_DMSPd_sk, &
          nt_bgc_DMS_sk, nt_bgc_C_sk

      integer (kind=int_kind) :: &
        nml_error ! namelist i/o error flag

      !-----------------------------------------------------------------
      ! namelist variables
      !-----------------------------------------------------------------

      namelist /zbgc_nml/  &
         tr_brine, bgc_data_dir, sil_data_type, nit_data_type, &
         restore_bgc, skl_bgc, &
         tr_bgc_C_sk, tr_bgc_chl_sk, tr_bgc_Am_sk, tr_bgc_Sil_sk, &
         tr_bgc_DMSPp_sk, tr_bgc_DMSPd_sk, tr_bgc_DMS_sk, &
         restart_bgc, restart_hbrine, phi_snow, bgc_flux_type

      !-----------------------------------------------------------------
      ! default values
      !-----------------------------------------------------------------

      tr_brine        = .false.  ! brine height differs from ice height
      restore_bgc     = .false.  ! restore bgc if true
      skl_bgc         = .false.  ! solve skeletal biochemistry in diffuse bio
      bgc_data_dir    = 'unknown_bgc_data_dir'
      sil_data_type   = 'default'
      nit_data_type   = 'default'
      tr_bgc_C_sk     = .false.  ! biogeochemistry, 
      tr_bgc_chl_sk   = .false.  ! biogeochemistry,
      tr_bgc_Am_sk    = .false.  ! biogeochemistry, 
      tr_bgc_Sil_sk   = .false.  ! biogeochemistry,
      tr_bgc_DMSPp_sk = .false.  ! biogeochemistry, trace gases (skeletal)
      tr_bgc_DMSPd_sk = .false.  ! biogeochemistry, trace gases (skeletal)
      tr_bgc_DMS_sk   = .false.  ! biogeochemistry, trace gases (skeletal) 
      restart_bgc     = .false.  ! biogeochemistry restart
      restart_hbrine  = .false.  ! hbrine restart
      phi_snow        = p5       ! snow porosity
      bgc_flux_type   = 'Jin2006'! type of ocean-ice poston velocity ('constant')

      !-----------------------------------------------------------------
      ! read from input file
      !-----------------------------------------------------------------

      call get_fileunit(nu_nml)

      if (my_task == master_task) then
         open (nu_nml, file=trim(nml_filename), status='old',iostat=nml_error)
         if (nml_error /= 0) then
            nml_error = -1
         else
            nml_error =  1
         endif 

         print*,'Reading zbgc_nml'
         do while (nml_error > 0)
            read(nu_nml, nml=zbgc_nml,iostat=nml_error)
         end do
         if (nml_error == 0) close(nu_nml)
      endif
      call broadcast_scalar(nml_error, master_task)
      if (nml_error /= 0) then
         call abort_ice('ice: error reading zbgc namelist')
      endif
      call release_fileunit(nu_nml)

      !-----------------------------------------------------------------
      ! brine
      !-----------------------------------------------------------------

      if (trim(runtype) == 'continue') restart_hbrine = .true.

      call broadcast_scalar(tr_brine,           master_task)
      call broadcast_scalar(restart_hbrine,     master_task)
      call broadcast_scalar(phi_snow,           master_task)

      nt_fbri = c0
      if (tr_brine) then
          nt_fbri = ntrcr + 1   ! ice volume fraction with salt
          ntrcr = ntrcr + 1
          trcr_depend(nt_fbri)   = 1  ! volume-weighted
          trcr_base  (nt_fbri,1) = 0  ! volume-weighted
          trcr_base  (nt_fbri,2) = 1  ! volume-weighted
          trcr_base  (nt_fbri,3) = 0  ! volume-weighted
      endif

      if (phi_snow .le. c0) phi_snow = c1-rhos/rhoi

      if (my_task == master_task) then
         write(nu_diag,1010) ' tr_brine                  = ', tr_brine
         if (tr_brine) then
         write(nu_diag,1010) ' restart_hbrine            = ', restart_hbrine
         write(nu_diag,1005) ' phi_snow                  = ', phi_snow
         endif
         write(nu_diag,1010) ' skl_bgc                   = ', skl_bgc
      endif

      !-----------------------------------------------------------------
      ! skeletal layer biogeochemistry
      !-----------------------------------------------------------------

      if (TRBGCS == 0 .and. skl_bgc) then
         write(nu_diag,*) &
            'WARNING: skl_bgc=T but 0 bgc tracers compiled'
         write(nu_diag,*) &
            'WARNING: setting skl_bgc = F'
         skl_bgc = .false.
      endif

      if (trim(runtype) == 'continue') restart_bgc = .true.

      call broadcast_scalar(skl_bgc,            master_task)
      call broadcast_scalar(restart_bgc,        master_task)

      if (skl_bgc) then
            tr_bgc_N_sk      = .true.   ! minimum NP biogeochemistry
            tr_bgc_Nit_sk    = .true.
      else
            tr_bgc_N_sk      = .false.
            tr_bgc_C_sk      = .false.
            tr_bgc_chl_sk    = .false.
            tr_bgc_Nit_sk    = .false.
            tr_bgc_Am_sk     = .false.
            tr_bgc_Sil_sk    = .false.
            tr_bgc_DMSPp_sk  = .false.
            tr_bgc_DMSPd_sk  = .false.
            tr_bgc_DMS_sk    = .false.
      endif

      call broadcast_scalar(bgc_flux_type,      master_task)
      call broadcast_scalar(restore_bgc,        master_task)
      call broadcast_scalar(bgc_data_dir,       master_task)
      call broadcast_scalar(sil_data_type,      master_task)
      call broadcast_scalar(nit_data_type,      master_task)
      call broadcast_scalar(tr_bgc_N_sk,        master_task)
      call broadcast_scalar(tr_bgc_C_sk,        master_task)
      call broadcast_scalar(tr_bgc_chl_sk,      master_task)
      call broadcast_scalar(tr_bgc_Nit_sk,      master_task)
      call broadcast_scalar(tr_bgc_Am_sk,       master_task)
      call broadcast_scalar(tr_bgc_Sil_sk,      master_task)
      call broadcast_scalar(tr_bgc_DMSPp_sk,    master_task)
      call broadcast_scalar(tr_bgc_DMSPd_sk,    master_task)
      call broadcast_scalar(tr_bgc_DMS_sk,      master_task)

      if (skl_bgc) then

      if (my_task == master_task) then

         write(nu_diag,1030) ' bgc_flux_type             = ', bgc_flux_type
         write(nu_diag,1010) ' restart_bgc               = ', restart_bgc
         write(nu_diag,1010) ' restore_bgc               = ', restore_bgc
         write(nu_diag,*)    ' bgc_data_dir              = ', &
                               trim(bgc_data_dir)
         write(nu_diag,*)    ' sil_data_type             = ', &
                               trim(sil_data_type)
         write(nu_diag,*)    ' nit_data_type             = ', &
                               trim(nit_data_type)
         write(nu_diag,1010) ' tr_bgc_N_sk               = ', tr_bgc_N_sk
         write(nu_diag,1010) ' tr_bgc_C_sk               = ', tr_bgc_C_sk
         write(nu_diag,1010) ' tr_bgc_chl_sk             = ', tr_bgc_chl_sk
         write(nu_diag,1010) ' tr_bgc_Nit_sk             = ', tr_bgc_Nit_sk
         write(nu_diag,1010) ' tr_bgc_Am_sk              = ', tr_bgc_Am_sk
         write(nu_diag,1010) ' tr_bgc_Sil_sk             = ', tr_bgc_Sil_sk
         write(nu_diag,1010) ' tr_bgc_DMSPp_sk           = ', tr_bgc_DMSPp_sk
         write(nu_diag,1010) ' tr_bgc_DMSPd_sk           = ', tr_bgc_DMSPd_sk
         write(nu_diag,1010) ' tr_bgc_DMS_sk             = ', tr_bgc_DMS_sk
        
      endif   ! master_task

      !-----------------------------------------------------------------
      ! assign tracer indices and dependencies
      !-----------------------------------------------------------------

         nbtrcr = 0
         nlt_bgc_NO = 0
         nlt_bgc_N = 0
         nlt_bgc_C = 0
         nlt_bgc_chl = 0
         nlt_bgc_NH = 0
         nlt_bgc_Sil = 0
         nlt_bgc_DMSPp = 0
         nlt_bgc_DMSPd = 0
         nlt_bgc_DMS = 0

         ntrcr = ntrcr + 1              ! algalN, required tracer
         nt_bgc_N_sk = ntrcr
         nbtrcr = nbtrcr + 1
         nlt_bgc_N = nbtrcr
      
         ntrcr = ntrcr + 1              ! nitrate, required tracer 
         nt_bgc_Nit_sk = ntrcr
         nbtrcr = nbtrcr + 1
         nlt_bgc_NO = nbtrcr

         if (tr_bgc_C_sk) then
             ntrcr = ntrcr + 1
             nt_bgc_C_sk = ntrcr
             nbtrcr = nbtrcr + 1
             nlt_bgc_C = nbtrcr
         endif    
         if (tr_bgc_chl_sk)then
             ntrcr = ntrcr + 1
             nt_bgc_chl_sk = ntrcr
             nbtrcr = nbtrcr + 1
             nlt_bgc_chl = nbtrcr
         endif 
         if (tr_bgc_Am_sk)then
             ntrcr = ntrcr + 1
             nt_bgc_Am_sk = ntrcr
             nbtrcr = nbtrcr + 1
             nlt_bgc_NH = nbtrcr
         endif    
         if (tr_bgc_Sil_sk)then
             ntrcr = ntrcr + 1
             nt_bgc_Sil_sk = ntrcr
             nbtrcr = nbtrcr + 1
             nlt_bgc_Sil = nbtrcr
         endif    
         if (tr_bgc_DMSPp_sk)then
             ntrcr = ntrcr + 1
             nt_bgc_DMSPp_sk = ntrcr
             nbtrcr = nbtrcr + 1
             nlt_bgc_DMSPp = nbtrcr
         endif    
         if (tr_bgc_DMSPd_sk)then
             ntrcr = ntrcr + 1
             nt_bgc_DMSPd_sk = ntrcr
             nbtrcr = nbtrcr + 1
             nlt_bgc_DMSPd = nbtrcr
         endif    
         if (tr_bgc_DMS_sk)then
             ntrcr = ntrcr + 1
             nt_bgc_DMS_sk = ntrcr
             nbtrcr = nbtrcr + 1
             nlt_bgc_DMS = nbtrcr
         endif  
      endif  ! skl_bgc

      if (nbtrcr > max_nbtrcr) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'nbtrcr > max_nbtrcr'
         write (nu_diag,*) 'nbtrcr, max_nbtrcr:',nbtrcr, max_nbtrcr
         call abort_ice ('ice: ice_zbgc error')
      endif

      if (ntrcr > max_ntrcr-1) then
         write(nu_diag,*) 'max_ntrcr-1 < number of namelist tracers'
         write(nu_diag,*) 'max_ntrcr-1 = ',max_ntrcr-1,' ntrcr = ',ntrcr
         call abort_ice('max_ntrcr-1 < number of namelist tracers')
      endif                               

      if (skl_bgc .and. TRBGCS < 2) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'comp_ice must have number of bgc tracers >= 2'
         write (nu_diag,*) 'number of bgc tracers compiled:',TRBGCS
         call abort_ice ('ice: ice_zbgc error')
      endif

      if (my_task == master_task) then
         if (skl_bgc) then
            write(nu_diag,1020)'nt_bgc_N_sk = ', nt_bgc_N_sk
            write(nu_diag,1020)'nt_bgc_Nit_sk = ', nt_bgc_Nit_sk
         endif
         if (tr_brine .or. skl_bgc) then
            write(nu_diag,1020)'nblyr = ', nblyr
            write(nu_diag,1020) 'ntrcr (w/ bgc) = ', ntrcr
         endif
      endif

      ! BGC layer model (on bottom "skeletal" layer)
      if (tr_bgc_N_sk)     trcr_depend(nt_bgc_N_sk)     = 0 ! algae  (skeletal)
      if (tr_bgc_C_sk)     trcr_depend(nt_bgc_C_sk)     = 0 ! 
      if (tr_bgc_chl_sk)   trcr_depend(nt_bgc_chl_sk)   = 0 ! 
      if (tr_bgc_Nit_sk)   trcr_depend(nt_bgc_Nit_sk)   = 0 ! nutrients
      if (tr_bgc_Am_sk)    trcr_depend(nt_bgc_Am_sk)    = 0 ! 
      if (tr_bgc_Sil_sk)   trcr_depend(nt_bgc_Sil_sk)   = 0 ! 
      if (tr_bgc_DMSPp_sk) trcr_depend(nt_bgc_DMSPp_sk) = 0 ! trace gases
      if (tr_bgc_DMSPd_sk) trcr_depend(nt_bgc_DMSPd_sk) = 0 !
      if (tr_bgc_DMS_sk)   trcr_depend(nt_bgc_DMS_sk)   = 0 !

      if (tr_bgc_N_sk)     bgc_tracer_type(nlt_bgc_N)     = c0 ! algae
      if (tr_bgc_C_sk)     bgc_tracer_type(nlt_bgc_C)     = c0 ! 
      if (tr_bgc_chl_sk)   bgc_tracer_type(nlt_bgc_chl)   = c0 ! 
      if (tr_bgc_Nit_sk)   bgc_tracer_type(nlt_bgc_NO)    = c1 ! nutrients
      if (tr_bgc_Am_sk)    bgc_tracer_type(nlt_bgc_NH)    = c1 ! 
      if (tr_bgc_Sil_sk)   bgc_tracer_type(nlt_bgc_Sil)   = c1 ! 
      if (tr_bgc_DMSPp_sk) bgc_tracer_type(nlt_bgc_DMSPp) = c0 ! trace gases
      if (tr_bgc_DMSPd_sk) bgc_tracer_type(nlt_bgc_DMSPd) = c1 !
      if (tr_bgc_DMS_sk)   bgc_tracer_type(nlt_bgc_DMS)   = c1 !
   
 1000    format (a30,2x,f9.2)  ! a30 to align formatted, unformatted statements
 1005    format (a30,2x,f9.6)  ! float
 1010    format (a30,2x,l6)    ! logical
 1020    format (a30,2x,i6)    ! integer
 1030    format (a30,   a8)    ! character

      end subroutine init_zbgc

!=======================================================================

!  Initialize vertical profile of biogeochemistry

      subroutine init_bgc

      use ice_algae, only: read_restart_bgc
      use ice_blocks, only: nx_block, ny_block
      use ice_communicate, only: my_task, master_task
      use ice_constants, only: c1, c0, c10, c5, p15, &
              field_type_scalar, field_loc_center
      use ice_domain, only: nblocks
      use ice_domain_size, only: max_blocks
      use ice_fileunits, only: nu_diag, nu_forcing
      use ice_flux, only: sss
      use ice_calendar, only: month
      use ice_read_write, only: ice_read, ice_open
      use ice_state, only: trcrn, aicen
      use ice_colpkg_tracers, only: &
          nt_bgc_N_sk, nt_bgc_Nit_sk, nt_bgc_chl_sk, nt_bgc_Am_sk, &
          nt_bgc_Sil_sk, nt_bgc_DMSPp_sk, nt_bgc_DMSPd_sk, &
          nt_bgc_DMS_sk, nt_bgc_C_sk

      integer (kind=int_kind) :: &
         i, j, iblk       , & ! horizontal indices
         nbits

      logical (kind=log_kind) :: &
         dbug             ! prints debugging output if true

      real (kind=dbl_kind), dimension (nx_block,ny_block,max_blocks) :: &
         work1

      dbug = .true.

      if (.not. skl_bgc) return
  
      if (restart_bgc) then       

         call read_restart_bgc

      else ! not restarting
         
      !-----------------------------------------------------------------
      ! default ocean values
      !-----------------------------------------------------------------

         sil   (:,:,:) = c10       ! ocean silicate       (mmol/m^3)
         nit   (:,:,:) = c5        ! ocean nitrate        (mmol/m^3)
         amm   (:,:,:) = c1        ! ocean ammonia        (mmol/m^3)
         dmsp  (:,:,:) = R_S2N*p15 ! sulfur cycle product (mmol/m^3)
         dms   (:,:,:) = c0        ! sulfur cycle product (mmol/m^3)
         algalN(:,:,:) = p15       ! algal concentration  (mmol/m^3)

      !-----------------------------------------------------------------
      !  skeletal layer model
      !-----------------------------------------------------------------

         if (tr_bgc_N_sk)     trcrn(:,:,nt_bgc_N_sk,    :,:) = &
                                                         p15/phi_sk*sk_l
         if (tr_bgc_C_sk)     trcrn(:,:,nt_bgc_C_sk,    :,:) = &
                                                   R_C2N*p15/phi_sk*sk_l
         if (tr_bgc_chl_sk)   trcrn(:,:,nt_bgc_chl_sk,  :,:) = &
                                                 R_chl2N*p15/phi_sk*sk_l
         if (tr_bgc_Nit_sk)   trcrn(:,:,nt_bgc_Nit_sk,  :,:) = &
                                                          c5/phi_sk*sk_l
         if (tr_bgc_Am_sk)    trcrn(:,:,nt_bgc_Am_sk,   :,:) = &
                                                          c1/phi_sk*sk_l
         if (tr_bgc_Sil_sk)   trcrn(:,:,nt_bgc_Sil_sk,  :,:) = &
                                                         c10/phi_sk*sk_l
         if (tr_bgc_DMSPp_sk) trcrn(:,:,nt_bgc_DMSPp_sk,:,:) = &
                                                   R_S2N*p15/phi_sk*sk_l
         if (tr_bgc_DMSPd_sk) trcrn(:,:,nt_bgc_DMSPd_sk,:,:) = c0
         if (tr_bgc_DMS_sk)   trcrn(:,:,nt_bgc_DMS_sk,  :,:) = c0
 
      !-----------------------------------------------------------------
      ! silicate
      !-----------------------------------------------------------------

         nbits = 64                ! double precision data

         if (tr_bgc_Sil_sk) then
         if (trim(sil_data_type) == 'clim') then  
            ! gx1 only
            sil_file = trim(bgc_data_dir)//'silicate_WOA2005_surface_monthly'

            if (my_task == master_task) then
               write (nu_diag,*) ' '
               write (nu_diag,*) 'silicate initialized from:'
               write (nu_diag,*) trim(sil_file)
            endif

            if (my_task == master_task) &
               call ice_open (nu_forcing, sil_file, nbits)

            call ice_read (nu_forcing, month, work1, 'rda8', dbug, &
                           field_loc_center, field_type_scalar)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  sil(i,j,iblk) = work1(i,j,iblk)
               enddo
               enddo
            enddo

            if (my_task == master_task) close(nu_forcing)

         else ! default
        
            ! use WOA2005_surface (winter or spring) for a specific location
            ! Bering (60, 180), Okhotsk (55, 150E),  Chukchi (70, 170W) 
            ! Labrador Sea (56, 50W), central(0,86)
            !          March:                      (25, 50, 30, 2.5, 20)
            ! mmol/m^3 Apr, May, Jun spring range: (20, 40, 10, 2.5, 20)
            !          Jan, Feb, Mar winter range: (20, 60, 25, 2.5, 20)
          
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  sil(i,j,iblk) = 30.0_dbl_kind  !chukchi, march
               enddo
               enddo
            enddo

         endif ! sil_data_type
         endif ! tr_bgc_Sil_sk

      !-----------------------------------------------------------------
      ! nitrate
      !-----------------------------------------------------------------

         if (tr_bgc_Nit_sk) then
         if (trim(nit_data_type) == 'clim') then
            ! gx1 only
            nit_file = trim(bgc_data_dir)//'nitrate_WOA2005_surface_monthly'

            if (my_task == master_task) then
               write (nu_diag,*) ' '
               write (nu_diag,*) 'nitrate initialized from:'
               write (nu_diag,*) trim(nit_file)
            endif

            if (my_task == master_task) &
               call ice_open (nu_forcing, nit_file, nbits)

            call ice_read (nu_forcing, month, work1, 'rda8', dbug, &
                           field_loc_center, field_type_scalar)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  nit(i,j,iblk) = work1(i,j,iblk)                 
               enddo
               enddo
            enddo
            
            if (my_task == master_task) close(nu_forcing)

         elseif (trim(nit_data_type) == 'sss') then
           
            if (my_task == master_task) then
               write (nu_diag,*) ' '
               write (nu_diag,*) 'nitrate initialized from salinity'
            endif
             
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  nit(i,j,iblk) = sss(i,j,iblk)        
               enddo
               enddo
            enddo

         else ! default
            
            if (my_task == master_task) then
               write (nu_diag,*) ' '
               write (nu_diag,*) 'nitrate initialized from March, Chukchi Sea'
            endif
             
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  nit(i,j,iblk) = c10
               enddo
               enddo
            enddo

         endif ! nit_data_type
         endif ! tr_bgc_Nit_sk
              
      endif    ! restart_bgc

      end subroutine init_bgc

!=======================================================================

      subroutine biogeochemistry (dt, iblk)

      use ice_algae, only: skl_biogeochemistry
      use ice_blocks, only: nx_block, ny_block, block, get_block
      use ice_brine, only: preflushing_changes, compute_microS_mushy, &
                           update_hbrine
      use ice_constants, only: c0, c1, puny
      use ice_domain, only: blocks_ice
      use ice_domain_size, only: ncat, nblyr
      use ice_flux, only: hin_old, meltbn, melttn, congeln, snoicen, &
                          sss, sst, meltsn, hmix
      use ice_arrays_column, only:  fswthrun
      use ice_state, only: aicen_init, vicen_init, aicen, vicen, vsnon, &
          trcrn
      use ice_colpkg_tracers, only: nt_fbri, tr_brine, ntrcr, nbtrcr, &
          nt_bgc_N_sk, nt_bgc_Nit_sk, nt_bgc_chl_sk, nt_bgc_Am_sk, &
          nt_bgc_Sil_sk, nt_bgc_DMSPp_sk, nt_bgc_DMSPd_sk, &
          nt_bgc_DMS_sk, nt_bgc_C_sk
      use ice_timers, only: timer_bgc, ice_timer_start, ice_timer_stop

      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step

      integer (kind=int_kind), intent(in) :: &
         iblk    ! block index

      ! local variables

      integer (kind=int_kind) :: &
         i, j, ij    , & ! horizontal indices
         ilo,ihi,jlo,jhi, & ! beginning and end of physical domain
         n               ! thickness category index

      integer (kind=int_kind) :: &
         icells          ! number of cells with aicen > puny

      integer (kind=int_kind), dimension(nx_block*ny_block) :: &
         indxi, indxj    ! indirect indices for cells with aicen > puny

      real (kind=dbl_kind), dimension (nx_block*ny_block) :: &
         hin         , & ! new ice thickness
         hsn         , & ! snow thickness  (m)
         hbr_old     , & ! old brine thickness before growh/melt
         kavg        , & ! average ice permeability (m^2)
         zphi_o      , & ! surface ice porosity 
         hbrin           ! brine height

      real (kind=dbl_kind), dimension (nx_block*ny_block,nblyr+2) :: &
      ! Defined on Bio Grid points
         bSin        , & ! salinity on the bio grid  (ppt)
         brine_sal   , & ! brine salinity (ppt)
         brine_rho   , & ! brine_density (kg/m^3)
      ! Defined on Bio Grid interfaces
         iphin       , & ! porosity 
         ibrine_sal  , & ! brine salinity  (ppt)
         ibrine_rho      ! brine_density (kg/m^3)

      real (kind=dbl_kind), dimension (nx_block,ny_block) :: &
         grow_Cn         ! C growth

      real (kind=dbl_kind), dimension (nx_block,ny_block,nbtrcr) :: &
         flux_bion       ! tracer flux to ocean

      type (block) :: &
         this_block      ! block information for current block

      if (tr_brine .or. skl_bgc) then

         call ice_timer_start(timer_bgc) ! biogeochemistry

      !-----------------------------------------------------------------
      ! initialize
      !-----------------------------------------------------------------

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         ! Define ocean tracer concentration
         do j = 1, ny_block
         do i = 1, nx_block
            if (tr_bgc_Nit_sk)   ocean_bio(i,j,nlt_bgc_NO   ,iblk) = nit   (i,j,iblk)
            if (tr_bgc_chl_sk)   ocean_bio(i,j,nlt_bgc_chl  ,iblk) = algalN(i,j,iblk)*R_chl2N
            if (tr_bgc_Am_sk)    ocean_bio(i,j,nlt_bgc_NH   ,iblk) = amm   (i,j,iblk)
            if (tr_bgc_C_sk)     ocean_bio(i,j,nlt_bgc_C    ,iblk) = algalN(i,j,iblk)*R_C2N
            if (tr_bgc_Sil_sk)   ocean_bio(i,j,nlt_bgc_Sil  ,iblk) = sil   (i,j,iblk)
            if (tr_bgc_DMSPp_sk) ocean_bio(i,j,nlt_bgc_DMSPp,iblk) = dmsp  (i,j,iblk)
            if (tr_bgc_DMSPd_sk) ocean_bio(i,j,nlt_bgc_DMSPd,iblk) = dmsp  (i,j,iblk)
            if (tr_bgc_DMS_sk)   ocean_bio(i,j,nlt_bgc_DMS  ,iblk) = dms   (i,j,iblk)
            if (tr_bgc_N_sk)     ocean_bio(i,j,nlt_bgc_N    ,iblk) = algalN(i,j,iblk)
         enddo
         enddo

         do n = 1, ncat
            
            hin_old(:,:,n,iblk) = c0
            flux_bion(:,:,:) = c0
            do j = jlo, jhi
            do i = ilo, ihi
               if (aicen_init(i,j,n,iblk) > puny) then 
                  hin_old(i,j,n,iblk) = vicen_init(i,j,n,iblk) &
                                      / aicen_init(i,j,n,iblk)
               else
                  first_ice(i,j,n,iblk) = .true.
                  if (tr_brine) trcrn(i,j,nt_fbri,n,iblk) = c1
               endif
            enddo
            enddo

            icells = 0
            do j = jlo, jhi
            do i = ilo, ihi
               if (aicen(i,j,n,iblk) > puny) then
                  icells = icells + 1
                  indxi(icells) = i
                  indxj(icells) = j
               endif
            enddo               ! i
            enddo               ! j

            if (icells > 0) then

      !-----------------------------------------------------------------
      ! brine dynamics
      !-----------------------------------------------------------------

            if (tr_brine) then 

               call preflushing_changes (nx_block,   ny_block,            &
                                icells,              n,                   &
                                indxi,               indxj,               &
                                aicen  (:,:,n,iblk),                      &
                                vicen  (:,:,n,iblk), vsnon  (:,:,n,iblk), &
                                meltbn (:,:,n,iblk), melttn (:,:,n,iblk), &
                                congeln(:,:,n,iblk), snoicen(:,:,n,iblk), &
                                hin_old(:,:,n,iblk),                      & 
                                trcrn  (:,:,nt_fbri,n,iblk),              &
                                dhbr_top(:,:,n,iblk),dhbr_bot(:,:,n,iblk),&
                                hbr_old,             hin,                 &
                                hsn,                 first_ice(:,:,n,iblk))

               ! Requires the average ice permeability = kavg(:)
               ! and the surface ice porosity = zphi_o(:)
               ! computed in "compute_microS" or from "thermosaline_vertical"

               call compute_microS_mushy (nx_block,  ny_block,            &
                                icells,              n,                   &
                                indxi,               indxj,               &
                                trcrn(:,:,:,n,iblk), hin_old(:,:,n,iblk), &
                                hbr_old,                                  &
                                sss  (:,:,iblk),     sst(:,:,iblk),       & 
                                bTiz (:,:,:,n,iblk), bphi(:,:,:,n,iblk),  &
                                kavg,                zphi_o,              &
                                bSin,                brine_sal,           &
                                brine_rho,           iphin,               &
                                ibrine_rho,          ibrine_sal)
 
!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
               do ij = 1, icells
                  i = indxi(ij)
                  j = indxj(ij)

               call update_hbrine (meltbn  (i,j,n,iblk), melttn(i,j,n,iblk),  &
                                   meltsn  (i,j,n,iblk), dt,                  &
                                   hin     (ij),         hsn   (ij),          &
                                   hin_old (i,j,n,iblk), hbrin (ij),          &
                                   hbr_old (ij),                              &
                                   trcrn   (i,j,nt_fbri,n,iblk),              &
                                   dhbr_top(i,j,n,iblk), dhbr_bot(i,j,n,iblk),&
                                   kavg    (ij),         zphi_o(ij),          &
                                   darcy_V (i,j,n,iblk))
                    
               hbri(i,j,iblk) = hbri(i,j,iblk) + hbrin(ij)*aicen_init(i,j,n,iblk)  
               enddo                     ! ij

            endif ! tr_brine

      !-----------------------------------------------------------------
      ! biogeochemistry
      !-----------------------------------------------------------------

            if (skl_bgc) then
               call skl_biogeochemistry (nx_block, ny_block,            &
                                         icells,   dt,                  &
                                         indxi,    indxj,               &  
                                         nbtrcr,                        &
                                         flux_bion(:,:,1:nbtrcr),       &
                                         ocean_bio(:,:,1:nbtrcr, iblk), &
                                         hmix     (:,:,          iblk), &
                                         aicen    (:,:,        n,iblk), & 
                                         meltbn   (:,:,        n,iblk), &
                                         congeln  (:,:,        n,iblk), &
                                         fswthrun (:,:,        n,iblk), &
                                         first_ice(:,:,        n,iblk), &
                                         trcrn    (:,:,1:ntrcr,n,iblk), &
                                         grow_Cn)

               call merge_bgc_fluxes_skl(nx_block, ny_block,                 &
                                         icells,                             &
                                         indxi,    indxj,                    &
                                         nbtrcr,                             &
                                         aicen_init(:,:,            n,iblk), &
                                         trcrn     (:,:,nt_bgc_N_sk,n,iblk), &
                                         flux_bion (:,:,1:nbtrcr),           &
                                         flux_bio  (:,:,1:nbtrcr,     iblk), &
                                         PP_net    (:,:,              iblk), &
                                         grow_net  (:,:,              iblk), &
                                         grow_Cn)
            endif  

            do ij = 1, icells
               i = indxi(ij)
               j = indxj(ij)
               first_ice(i,j,n,iblk) = .false.
            enddo
   
            endif               ! icells      
         enddo                  ! ncat

         call ice_timer_stop(timer_bgc) ! biogeochemistry

      endif  ! tr_brine .or. skl_bgc

      end subroutine biogeochemistry

!=======================================================================

! Aggregate flux information from all ice thickness categories
! for skeletal layer biogeochemistry
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL

      subroutine merge_bgc_fluxes_skl (nx_block, ny_block,   &
                               icells,               &
                               indxi,    indxj,      &
                               nbtrcr, &
                               aicen,  algal_N,        &
                               flux_bion, flux_bio,  &
                               PP_net, grow_net,     &
                               grow_Cn)

      use ice_constants, only: c1

      integer (kind=int_kind), intent(in) :: &
          nx_block, ny_block, & ! block dimensions
          icells            , & ! number of cells with aicen > puny
          nbtrcr                ! number of bgc tracers

      integer (kind=int_kind), dimension(nx_block*ny_block), &
          intent(in) :: &
          indxi, indxj          ! compressed indices for cells with aicen > puny

      ! single category fluxes
      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in):: &          
          aicen                 ! category ice area fraction

      real (kind=dbl_kind), dimension (nx_block,ny_block), &
          intent(in) :: &
          algal_N               ! (mmol N/m^2)
     
      real (kind=dbl_kind), dimension(nx_block,ny_block,nbtrcr), &
          intent(in):: &
          flux_bion             ! all bio fluxes to ocean, on categories

      real (kind=dbl_kind), dimension(nx_block,ny_block,nbtrcr), &
          intent(inout):: &
          flux_bio              ! all bio fluxes to ocean, aggregated

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(in):: & 
          grow_Cn               ! specific growth (/s) 

      ! history output
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
          intent(inout):: & 
          PP_net , &            ! Bulk net PP (mg C/m^2/s)
          grow_net              ! net specific growth (/s)
      
      ! local variables

      integer (kind=int_kind) :: &
          ij, i, j, &   ! horizontal indices
          k             ! tracer indice
    
      !-----------------------------------------------------------------
      ! Merge fluxes
      !-----------------------------------------------------------------

!DIR$ CONCURRENT !Cray
!cdir nodep      !NEC
!ocl novrec      !Fujitsu
      do ij = 1, icells
         i = indxi(ij)
         j = indxj(ij)
         do k = 1,nbtrcr
            flux_bio (i,j,k) = flux_bio(i,j,k) + flux_bion(i,j,k)*aicen(i,j)
         enddo
        
         PP_net   (i,j) = PP_net  (i,j) &
                        + algal_N(i,j)*phi_sk*grow_Cn(i,j)*(c1-fr_resp) &
                        * R_C2N*R_gC2molC * aicen(i,j)
         grow_net (i,j) = grow_net(i,j) + grow_Cn(i,j) * phi_sk*aicen(i,j) 
      enddo                     ! ij
      
      end subroutine merge_bgc_fluxes_skl

!=======================================================================

! Initialize bgc fields written to history files
!
! authors: Nicole Jeffery, LANL
!          Elizabeth C. Hunke, LANL

      subroutine init_history_bgc

      use ice_constants, only: c0

      PP_net     (:,:,:)   = c0
      grow_net   (:,:,:)   = c0
      hbri       (:,:,:)   = c0
      flux_bio   (:,:,:,:) = c0
      flux_bio_ai(:,:,:,:) = c0

      end subroutine init_history_bgc

!=======================================================================

! Adjust biogeochemical tracers when new frazil ice forms

      subroutine add_new_ice_bgc (dt,                               &
                                  aicen_init, vicen_init, vi0_init, &
                                  aicen,      vicen,      vi0new,   &
                                  trcrn,      nbtrcr,               &
                                  sss,        ocean_bio,  flux_bio, &
                                  hsurp,      l_stop,     stop_label)

      use ice_constants, only: c0, c1, puny
      use ice_domain_size, only: ncat
      use ice_fileunits, only: nu_diag
      use ice_itd, only: column_sum, &
                         column_conservation_check
      use ice_colpkg_tracers, only: tr_brine, nt_fbri
      use ice_timers, only: timer_bgc, ice_timer_start, ice_timer_stop

      real (kind=dbl_kind), intent(in) :: &
         dt              ! time step (s)

      real (kind=dbl_kind), dimension (:), intent(in) :: &
         aicen_init  , & ! initial concentration of ice
         vicen_init  , & ! intiial volume per unit area of ice  (m)
         aicen       , & ! concentration of ice
         vicen           ! volume per unit area of ice          (m)

      real (kind=dbl_kind), dimension (:,:), intent(inout) :: &
         trcrn           ! ice tracers

      real (kind=dbl_kind), intent(in) :: &
         sss         , & ! sea surface salinity (ppt)
         vi0_init    , & ! volume of new ice added to cat 1 (initial)
         vi0new      , & ! volume of new ice added to cat 1
         hsurp           ! thickness of new ice added to each cat

      integer (kind=int_kind), intent(in) :: &
         nbtrcr         ! number of biology tracers

      real (kind=dbl_kind), dimension (:), intent(inout) :: &
         flux_bio   ! tracer flux to ocean from biology (mmol/m^2/s) 
        
       real (kind=dbl_kind), dimension (:), intent(in) :: &
         ocean_bio   ! ocean concentration of biological tracer

      logical (kind=log_kind), intent(out) :: &
         l_stop          ! if true, abort on return

      character (char_len), intent(out) :: stop_label

      ! local variables

      integer (kind=int_kind) :: &
         n           , & ! ice category index
         k               ! ice layer index

      real (kind=dbl_kind) :: &
         vbri1       , & ! starting volume of existing brine
         vbri_init   , & ! brine volume summed over categories
         vbri_final  , & ! brine volume summed over categories
         vsurp       , & ! volume of new ice added to each cat
         vtmp            ! total volume of new and old ice
        
      real (kind=dbl_kind), dimension (ncat) :: &
         vbrin           ! trcrn(i,j,nt_fbri,n)*vicen(i,j,n) 

      character (len=char_len) :: &
         fieldid         ! field identifier

      call ice_timer_start(timer_bgc) ! biogeochemistry

      !-----------------------------------------------------------------     
      ! brine
      !-----------------------------------------------------------------     
      do n = 1, ncat
         vbrin(n) = vicen_init(n)
         if (tr_brine) vbrin(n) =  trcrn(nt_fbri,n)*vicen_init(n)
      enddo

      call column_sum (ncat, vbrin(:), vbri_init)

      vbri_init = vbri_init + vi0_init

      !-----------------------------------------------------------------     
      ! ocean flux
      !-----------------------------------------------------------------     
      do k = 1, nbtrcr  ! only correct for dissolved tracers
         flux_bio(k) = flux_bio(k) - vi0_init/dt*ocean_bio(k) & 
                     * (bgc_tracer_type(k)*initbio_frac &
                           + (c1-bgc_tracer_type(k)))
      enddo

      !-----------------------------------------------------------------
      ! Distribute bgc in new ice volume among all ice categories by 
      ! increasing ice thickness, leaving ice area unchanged.
      !-----------------------------------------------------------------

      if (hsurp > c0) then   ! add ice to all categories

         ! Diffuse_bio handles concentration changes from ice growth/melt
         ! ice area does not change
         ! add salt to the bottom 

         do n = 1,ncat
            vbrin(n) = vbrin(n) + hsurp * aicen_init(n) 
            if (tr_brine) then
               trcrn(nt_fbri,n) = c1
               if (vicen(n) > c0) trcrn(nt_fbri,n) = vbrin(n)/vicen(n)
            endif
         enddo              ! n

      endif ! hsurp > 0

      !-----------------------------------------------------------------
      ! Combine bgc in new ice grown in open water with category 1 ice.
      !-----------------------------------------------------------------

      if (vi0new > c0) then  ! add ice to category 1

         ! Diffuse_bio handles concentration changes from ice growth/melt
         ! ice area changes
         ! add salt throughout

         vbri1    = vbrin(1) 
         vbrin(1) = vbrin(1) + vi0new
         if (tr_brine) then
            trcrn(nt_fbri,1) = c1
            if (vicen(1) > c0) trcrn(nt_fbri,1) = vbrin(1)/vicen(1)
         endif

      endif ! vi0new > 0
        
      if (tr_brine) then
         call column_sum (ncat, vbrin(:), vbri_final)
         fieldid = 'vbrin, add_new_ice'
         call column_conservation_check (fieldid,               &
                                         vbri_init, vbri_final, &
                                         puny,                  &
                                         l_stop,    nu_diag)
         if (l_stop) then
            stop_label = 'add_new_ice_bgc: Column conservation error'
            return
         endif
      endif

      call ice_timer_stop(timer_bgc) ! biogeochemistry

      end subroutine add_new_ice_bgc

!=======================================================================

      end module ice_zbgc

!=======================================================================
