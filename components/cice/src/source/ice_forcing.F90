!=======================================================================
!
!BOP
!
! !MODULE: ice_forcing - reads and interpolates input forcing data
!
! !DESCRIPTION:
!
! Reads and interpolates forcing data for atmosphere and ocean quantities.
!
! !REVISION HISTORY:
!  SVN:$Id: ice_forcing.F90 140 2008-07-25 20:15:53Z eclare $
!
! authors: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! 2004 WHL: Block structure added
! 2005 WHL: ECMWF option added
! 2006 ECH: LY option added
! 2006 WHL: Module name changed from ice_flux_in
! 2006 ECH: Fixed bugs, rearranged routines, edited comments, etc.
!           Added NCAR ocean forcing file
!           Converted to free source form (F90)
! 2007: netcdf version of read_data added by Alison McLaren, Met Office
!
! !INTERFACE:
!
      module ice_forcing
!
! !USES:
!
      use ice_kinds_mod
      use ice_blocks, only: nx_block, ny_block
      use ice_domain_size
      use ice_communicate, only: my_task, master_task
      use ice_constants
      use ice_calendar, only: istep, istep1, time, time_forc, year_init, &
                              sec, mday, month, nyr, yday, daycal, dayyr, &
                              daymo, days_per_year
      use ice_fileunits
      use ice_atmo, only: calc_strair
      use ice_exit
      use ice_timers
!
!EOP
!
      implicit none
      save

      integer (kind=int_kind) :: &
         ycycle          , & ! number of years in forcing cycle
         fyear_init      , & ! first year of data in forcing cycle
         fyear           , & ! current year in forcing cycle
         fyear_final         ! last year in cycle

      character (char_len_long) :: &        ! input data file names
         height_file, &
          uwind_file, &
          vwind_file, &
           wind_file, &
          strax_file, &
          stray_file, &
           potT_file, &
           tair_file, &
          humid_file, &
           rhoa_file, &
            fsw_file, &
            flw_file, &
           rain_file, &
            sst_file, &
            sss_file, &
           pslv_file, &
         sublim_file, &
           snow_file  

      character (char_len_long), dimension(ncat) :: &  ! input data file names
        topmelt_file, &
        botmelt_file

      real (kind=dbl_kind) :: &
           c1intp, c2intp , & ! interpolation coefficients
           ftime              ! forcing time (for restart)

      integer (kind=int_kind) :: &
           oldrecnum = 0  , & ! old record number (save between steps)
           oldrecslot = 1     ! old record slot (save between steps)

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
          cldf                ! cloud fraction

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks) :: &
            fsw_data, & ! field values at 2 temporal data points
           cldf_data, &
          fsnow_data, &
           Tair_data, &
           uatm_data, &
           vatm_data, &
           wind_data, &
          strax_data, &
          stray_data, &
             Qa_data, &
           rhoa_data, &
           potT_data, &
           zlvl_data, &
            flw_data, &
            sst_data, &
            sss_data, & 
           uocn_data, &
           vocn_data, &
         sublim_data, &
          frain_data

      real (kind=dbl_kind), & 
           dimension(nx_block,ny_block,2,max_blocks,ncat) :: &
        topmelt_data, &
        botmelt_data

      character(char_len) :: & 
         atm_data_format, & ! 'bin'=binary or 'nc'=netcdf
         ocn_data_format, & ! 'bin'=binary or 'nc'=netcdf
         atm_data_type, & ! 'default', 'monthly', 'ncar', 'ecmwf', 
                          ! 'LYq' or 'hadgem'
         sss_data_type, & ! 'default', 'clim', or 'ncar'
         sst_data_type, & ! 'default', 'clim', 'ncar', 
                          !     'hadgem_sst' or 'hadgem_sst_uvocn'
         precip_units     ! 'mm_per_month', 'mm_per_sec', 'mks'
 
      character(char_len_long) :: & 
         atm_data_dir , & ! top directory for atmospheric data
         ocn_data_dir , & ! top directory for ocean data
         oceanmixed_file  ! file name for ocean forcing data

      integer (kind=int_kind), parameter :: & 
         nfld = 8    ! number of fields to search for in forcing file

      ! as in the dummy atm (latm)
      real (kind=dbl_kind), parameter :: &
         frcvdr = 0.28_dbl_kind, & ! frac of incoming sw in vis direct band
         frcvdf = 0.24_dbl_kind, & ! frac of incoming sw in vis diffuse band
         frcidr = 0.31_dbl_kind, & ! frac of incoming sw in near IR direct band
         frcidf = 0.17_dbl_kind    ! frac of incoming sw in near IR diffuse band

      real (kind=dbl_kind), &
       dimension (nx_block,ny_block,max_blocks,nfld,12) :: & 
         ocn_frc_m   ! ocn data for 12 months

      logical (kind=log_kind) :: &
         restore_sst                 ! restore sst if true

      integer (kind=int_kind) :: &
         trestore                    ! restoring time scale (days)

      real (kind=dbl_kind) :: & 
         trest                       ! restoring time scale (sec)

      logical (kind=log_kind) :: &
         dbug             ! prints debugging output if true

!=======================================================================

      contains

!=======================================================================
!
!BOP
!
! !IROUTINE: init_forcing_atmo - initialize atmospheric forcing
!
! !INTERFACE:
!
      subroutine init_forcing_atmo
!
! !DESCRIPTION:
!
! Determine the current and final year of the forcing cycle based on
! namelist input; initialize the forcing data filenames.
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      fyear       = fyear_init + mod(nyr-1,ycycle) ! current year
      fyear_final = fyear_init + ycycle - 1 ! last year in forcing cycle

      if (trim(atm_data_type) /= 'default' .and. &
                          my_task == master_task) then
         write (nu_diag,*) ' Initial forcing data year = ',fyear_init
         write (nu_diag,*) ' Final   forcing data year = ',fyear_final
      endif

    !-------------------------------------------------------------------
    ! Get filenames for input forcing data     
    !-------------------------------------------------------------------

      ! default forcing values from init_flux_atm
      if (trim(atm_data_type) == 'ncar') then
         call NCAR_files(fyear)
      elseif (trim(atm_data_type) == 'ecmwf') then
         call ecmwf_files(fyear)    
      elseif (trim(atm_data_type) == 'LYq') then
         call LY_files(fyear)
      elseif (trim(atm_data_type) == 'hadgem') then
         call hadgem_files(fyear)
      elseif (trim(atm_data_type) == 'monthly') then
         call monthly_files(fyear)
      endif

      end subroutine init_forcing_atmo

!=======================================================================
!BOP
!
! !IROUTINE: init_forcing_ocn - initialize sss and sst
!
! !INTERFACE:
!
      subroutine init_forcing_ocn(dt)
!
! !DESCRIPTION:
!
! Set sea surface salinity and freezing temperature to annual mean value 
!  using a 12-month climatology.
! Read sst data for current month, and adjust sst based on freezing
! temperature.  No interpolation in time.

! Note: SST is subsequently prognosed if CICE is run with a mixed layer
! ocean (oceanmixed\_ice = T), and can be restored to data 
! (restore\_sst = T). SSS is not prognosed by CICE. 
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      use ice_domain, only: nblocks
      use ice_flux, only: sss, sst, Tf, Tfrzpt
      use ice_work, only:  work1
      use ice_read_write
      use ice_therm_vertical, only: ustar_scale
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
           dt                   ! time step
!
!EOP
!
      integer (kind=int_kind) :: &
           i, j, iblk       , & ! horizontal indices
           k                , & ! month index
           fid              , & ! file id for netCDF file 
          nbits

      logical (kind=log_kind) :: diag

      character (char_len) :: & 
            fieldname    ! field name in netcdf file

      nbits = 64                ! double precision data

    !-------------------------------------------------------------------
    ! Sea surface salinity (SSS)
    ! initialize to annual climatology created from monthly data
    !-------------------------------------------------------------------

      if (trim(sss_data_type) == 'clim') then

!         sss_file = trim(ocn_data_dir)//'sss_Lev.mm'
            sss_file = trim(ocn_data_dir)//'sss.mm.100x116.da' ! gx3 only
!!!         sss_file = trim(ocn_data_dir)//'sss_12.r'

         if (my_task == master_task) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'SSS climatology computed from:'
            write (nu_diag,*) trim(sss_file)
         endif

         if (my_task == master_task) &
              call ice_open (nu_forcing, sss_file, nbits)

         sss(:,:,:) = c0

         do k = 1,12            ! loop over 12 months
            call ice_read (nu_forcing, k, work1, 'rda8', dbug, &
                           field_loc_center, field_type_scalar)
            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block
                  sss(i,j,iblk) = sss(i,j,iblk) + work1(i,j,iblk)
               enddo
               enddo
            enddo
         enddo                  ! k

         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sss(i,j,iblk) = sss(i,j,iblk) / c12   ! annual average
               sss(i,j,iblk) = max(sss(i,j,iblk),c0)
               if (trim(Tfrzpt) == 'constant') then
                  Tf (i,j,iblk) = -1.8_dbl_kind ! deg C
               else ! default:  Tfrzpt = 'linear_S'
                  Tf (i,j,iblk) = -depressT * sss(i,j,iblk) ! deg C
               endif
            enddo
            enddo
         enddo

         ! close file
         if (my_task == master_task) close(nu_forcing)

      endif                     ! sss_data_type

    !-------------------------------------------------------------------
    ! Sea surface temperature (SST)
    ! initialize to data for current month
    !-------------------------------------------------------------------

      if (restore_sst) then
         if (trestore == 0) then
            trest = dt          ! use data instantaneously
         else
            trest = real(trestore,kind=dbl_kind) * secday ! seconds
         endif
      endif

      if (trim(sst_data_type) == 'clim') then

         if (nx_global == 320) then ! gx1
            sst_file = trim(ocn_data_dir)//'sst_clim_hurrell.dat'
         else                   ! gx3
!            sst_file = trim(ocn_data_dir)//'sst_Lev.mm'
            sst_file = trim(ocn_data_dir)//'sst.mm.100x116.da'
!!!            sst_file = trim(ocn_data_dir)//'sst_12.r'
         endif

         if (my_task == master_task) then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'Initial SST file:', trim(sst_file)
         endif

         if (my_task == master_task) &
              call ice_open (nu_forcing, sst_file, nbits)

         call ice_read (nu_forcing, month, sst, 'rda8', dbug, &
                        field_loc_center, field_type_scalar)

         if (my_task == master_task) close(nu_forcing)

         ! Make sure sst is not less than freezing temperature Tf
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sst(i,j,iblk) = max(sst(i,j,iblk),Tf(i,j,iblk))
            enddo
            enddo
         enddo

      endif                     ! init_sst_data


      if (trim(sst_data_type) == 'hadgem_sst' .or.  &
          trim(sst_data_type) == 'hadgem_sst_uvocn') then

       	 diag = .true.   ! write diagnostic information 

         sst_file = trim (ocn_data_dir)//'MONTHLY/sst.1997.nc'

       	 if (my_task == master_task) then

             write (nu_diag,*) ' '
             write (nu_diag,*) 'Initial SST file:', trim(sst_file)

             call ice_open_nc(sst_file,fid)

         endif
 
         fieldname='sst'
         call ice_read_nc(fid,month,fieldname,sst,diag)

         if (my_task == master_task) call ice_close_nc(fid)  

         ! Make sure sst is not less than freezing temperature Tf
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sst(i,j,iblk) = max(sst(i,j,iblk),Tf(i,j,iblk))
            enddo
            enddo
         enddo

      endif                        ! sst_data_type


      if (trim(sst_data_type) == 'ncar' .or.  &
          trim(sss_data_type) == 'ncar') then
         call ocn_data_ncar_init
      endif

      ! set ustar_scale for case of zero currents
      ! default value of c1 (for nonzero currents) set in init_thermo_vertical

      if (trim(sst_data_type) /= 'ncar' .or.  &
          trim(sss_data_type) /= 'ncar' .or.  & 
          trim(sst_data_type) /= 'hadgem_sst_uvocn') then
         ustar_scale = c10            ! for zero currents
      endif

      end subroutine init_forcing_ocn

!=======================================================================
!BOP
!
! !IROUTINE: get_forcing_atmo - Get atmospheric forcing data and interpolate
!
! !INTERFACE:
!
      subroutine get_forcing_atmo
!
! !DESCRIPTION:
!
! Get atmospheric forcing data and interpolate as necessary
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      use ice_boundary, only: ice_HaloUpdate
      use ice_domain
      use ice_blocks
      use ice_flux
      use ice_state
      use ice_grid, only: ANGLET, hm
!
!EOP
!
      integer (kind=int_kind) :: &
         iblk, &              ! block index
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      type (block) :: &
         this_block           ! block information for current block
      
      fyear = fyear_init + mod(nyr-1,ycycle)  ! current year
      if (trim(atm_data_type) /= 'default' .and. istep <= 1 &
                   .and. my_task == master_task) then
         write (nu_diag,*) ' Current forcing data year = ',fyear
      endif

      ftime = time         ! forcing time
      time_forc = ftime    ! for restarting

    !-------------------------------------------------------------------
    ! Read and interpolate atmospheric data
    !-------------------------------------------------------------------

      if (trim(atm_data_type) == 'ncar') then
         call ncar_data
      elseif (trim(atm_data_type) == 'ecmwf') then
         call ecmwf_data
      elseif (trim(atm_data_type) == 'LYq') then
         call LY_data
      elseif (trim(atm_data_type) == 'hadgem') then
         call hadgem_data
      elseif (trim(atm_data_type) == 'monthly') then
         call monthly_data
      else    ! default values set in init_flux
         return
      endif

    !-------------------------------------------------------------------
    ! Convert forcing data to fields needed by ice model
    !-------------------------------------------------------------------

      do iblk = 1, nblocks

         this_block = get_block(blocks_ice(iblk),iblk)         
         ilo = this_block%ilo
         ihi = this_block%ihi
         jlo = this_block%jlo
         jhi = this_block%jhi

         call prepare_forcing (nx_block, ny_block, &
                               ilo, ihi, jlo, jhi, &
                               hm    (:,:,iblk),   &
                               Tair  (:,:,iblk),   &
                               fsw   (:,:,iblk),   &   
                               cldf  (:,:,iblk),   &
                               flw   (:,:,iblk),   &
                               frain (:,:,iblk),   &
                               fsnow (:,:,iblk),   &
                               Qa    (:,:,iblk),   &
                               rhoa  (:,:,iblk),   &
                               uatm  (:,:,iblk),   &
                               vatm  (:,:,iblk),   &
                               strax (:,:,iblk),   &
                               stray (:,:,iblk),   &
                               zlvl  (:,:,iblk),   &
                               wind  (:,:,iblk),   &
                               swvdr (:,:,iblk),   &
                               swvdf (:,:,iblk),   &
                               swidr (:,:,iblk),   &
                               swidf (:,:,iblk),   &
                               potT  (:,:,iblk),   &
                               ANGLET(:,:,iblk),   &
                               trcr  (:,:,nt_Tsfc,iblk), &
                               sst   (:,:,iblk),   &
                               aice  (:,:,iblk) )

      enddo                     ! iblk

      call ice_timer_start(timer_bound)
      call ice_HaloUpdate (swvdr,             halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate (swvdf,             halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate (swidr,             halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_HaloUpdate (swidf,             halo_info, &
                           field_loc_center,  field_type_scalar)
      call ice_timer_stop(timer_bound)

      end subroutine get_forcing_atmo

!=======================================================================
!BOP
!
! !IROUTINE: get_forcing_ocn - interpolate sss, sst; restore sst
!
! !INTERFACE:
!
      subroutine get_forcing_ocn (dt)
!
! !DESCRIPTION:
!
! Read and interpolate annual climatologies of SSS and SST.
! Restore model SST to data if desired.
! Interpolate ocean fields to U grid if necessary.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain, only: nblocks
      use ice_ocean
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      if (trim(sst_data_type) == 'clim' .or.  &
          trim(sss_data_type) == 'clim') then
         call ocn_data_clim(dt)
      elseif (trim(sst_data_type) == 'ncar' .or.  &
              trim(sss_data_type) == 'ncar') then
         call ocn_data_ncar(dt)      
      elseif (trim(sst_data_type) == 'hadgem_sst' .or.  &
              trim(sst_data_type) == 'hadgem_sst_uvocn') then
         call ocn_data_hadgem(dt)  
      endif

      end subroutine get_forcing_ocn

!=======================================================================
!
!BOP
!
! !IROUTINE: read_data - Read data needed for interpolation
!
! !INTERFACE:
!
      subroutine read_data (flag, recd, yr, ixm, ixx, ixp, &
                            maxrec, data_file, field_data, &
                            field_loc, field_type)
!
! !DESCRIPTION:
!
! If data is at the beginning of a one-year record, get data from
!  the previous year.
! If data is at the end of a one-year record, get data from the
!  following year.
! If no earlier data exists (beginning of fyear_init), then
!  (1) For monthly data, get data from the end of fyear_final.
!  (2) For more frequent data, let the ixm value equal the
!      first value of the year.
! If no later data exists (end of fyear_final), then
!  (1) For monthly data, get data from the beginning of fyear_init.
!  (2) For more frequent data, let the ipx value
!      equal the last value of the year.
! In other words, we assume persistence when daily or 6-hourly
!   data is missing, and we assume periodicity when monthly data
!   is missing.
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      use ice_read_write
      use ice_diagnostics, only: check_step
!
! !INPUT/OUTPUT PARAMETERS:
!
      logical (kind=log_kind), intent(in) :: flag

      integer (kind=int_kind), intent(in) :: &
         recd                , & ! baseline record number
         yr                  , & ! year of forcing data
         ixm, ixx, ixp       , & ! record numbers of 3 data values
                                 ! relative to recd
         maxrec                  ! maximum record value

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), &
         intent(out) :: &
         field_data              ! 2 values needed for interpolation

      integer (kind=int_kind), intent(in) :: &
           field_loc, &      ! location of field on staggered grid
           field_type        ! type of field (scalar, vector, angle)
!
!EOP
!
      character (char_len_long) :: &
         data_file               ! data file to be read

      integer (kind=int_kind) :: &
         nbits            , & ! = 32 for single precision, 64 for double
         nrec             , & ! record number to read
         n2, n4           , & ! like ixm and ixp, but
                              ! adjusted at beginning and end of data
         arg                  ! value of time argument in field_data

      call ice_timer_start(timer_readwrite)  ! reading/writing

      nbits = 64              ! double precision data

      if (istep1 > check_step) dbug = .true.  !! debugging

      if (my_task==master_task .and. (dbug)) then
         write(nu_diag,*) '  ', trim(data_file)
      endif

      if (flag) then

      !-----------------------------------------------------------------
      ! Initialize record counters
      ! (n2, n4 will change only at the very beginning or end of
      !  a forcing cycle.)
      !-----------------------------------------------------------------
         n2 = ixm
         n4 = ixp
         arg = 0

      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------

         if (ixm /= 99) then
         ! currently in first half of data interval
            if (ixx <= 1) then
               if (yr > fyear_init) then ! get data from previous year
                  call file_year (data_file, yr-1)
               else             ! yr = fyear_init, no prior data exists
                  if (maxrec > 12) then ! extrapolate from first record
                     if (ixx == 1) n2 = ixx
                  else          ! go to end of fyear_final
                     call file_year (data_file, fyear_final)
                  endif
               endif            ! yr > fyear_init
            endif               ! ixx <= 1

            call ice_open (nu_forcing, data_file, nbits)

            arg = 1
            nrec = recd + n2
            call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
                           'rda8', dbug, field_loc, field_type)

            if (ixx==1 .and. my_task == master_task) close(nu_forcing)
         endif                  ! ixm ne 99

         ! always read ixx data from data file for current year
         call file_year (data_file, yr)
         call ice_open (nu_forcing, data_file, nbits)

         arg = arg + 1
         nrec = recd + ixx
         call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
                        'rda8', dbug, field_loc, field_type)

         if (ixp /= 99) then
         ! currently in latter half of data interval
            if (ixx==maxrec) then
               if (yr < fyear_final) then ! get data from following year
                  if (my_task == master_task) close(nu_forcing)
                  call file_year (data_file, yr+1)
                  call ice_open (nu_forcing, data_file, nbits)
               else             ! yr = fyear_final, no more data exists
                  if (maxrec > 12) then ! extrapolate from ixx
                     n4 = ixx
                  else          ! go to beginning of fyear_init
                     if (my_task == master_task) close(nu_forcing)
                     call file_year (data_file, fyear_init)

                     call ice_open (nu_forcing, data_file, nbits)

                  endif
               endif            ! yr < fyear_final
            endif               ! ixx = maxrec

            arg = arg + 1
            nrec = recd + n4
            call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
                           'rda8', dbug, field_loc, field_type)
         endif                  ! ixp /= 99

         if (my_task == master_task) close(nu_forcing)

      endif                     ! flag

      call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine read_data

!=======================================================================
!
!BOP
!
! !IROUTINE: read_data_nc - Read netcdf data needed for interpolation
!
! !INTERFACE:
!
      subroutine read_data_nc (flag, recd, yr, ixm, ixx, ixp, &
                            maxrec, data_file, fieldname, field_data, &
                            field_loc, field_type)
!
! !DESCRIPTION:
!
! If data is at the beginning of a one-year record, get data from
!  the previous year.
! If data is at the end of a one-year record, get data from the
!  following year.
! If no earlier data exists (beginning of fyear_init), then
!  (1) For monthly data, get data from the end of fyear_final.
!  (2) For more frequent data, let the ixm value equal the
!      first value of the year.
! If no later data exists (end of fyear_final), then
!  (1) For monthly data, get data from the beginning of fyear_init.
!  (2) For more frequent data, let the ipx value
!      equal the last value of the year.
! In other words, we assume persistence when daily or 6-hourly
!   data is missing, and we assume periodicity when monthly data
!   is missing.
!
! !REVISION HISTORY:
!
! Adapted by Alison McLaren, Met Office from read_data
!
! !USES:
!
      use ice_read_write
      use ice_diagnostics, only: check_step
!
! !INPUT/OUTPUT PARAMETERS:
!
      logical (kind=log_kind), intent(in) :: flag

      integer (kind=int_kind), intent(in) :: &
         recd                , & ! baseline record number
         yr                  , & ! year of forcing data
         ixm, ixx, ixp       , & ! record numbers of 3 data values
                                 ! relative to recd
         maxrec                  ! maximum record value

      character (char_len_long) :: &
         data_file               ! data file to be read

      character (char_len), intent(in) :: &
         fieldname               ! field name in netCDF file

      integer (kind=int_kind), intent(in) :: &
           field_loc, &      ! location of field on staggered grid
           field_type        ! type of field (scalar, vector, angle)

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), &
         intent(out) :: &
         field_data              ! 2 values needed for interpolation
!
!EOP
!
#ifdef ncdf 
      integer (kind=int_kind) :: &
         nrec             , & ! record number to read
         n2, n4           , & ! like ixm and ixp, but
                              ! adjusted at beginning and end of data
         arg              , & ! value of time argument in field_data
         fid                  ! file id for netCDF routines


      call ice_timer_start(timer_readwrite)  ! reading/writing

      if (istep1 > check_step) dbug = .true.  !! debugging

      if (my_task==master_task .and. (dbug)) then
         write(nu_diag,*) '  ', trim(data_file)
      endif

      if (flag) then

      !-----------------------------------------------------------------
      ! Initialize record counters
      ! (n2, n4 will change only at the very beginning or end of
      !  a forcing cycle.)
      !-----------------------------------------------------------------
         n2 = ixm
         n4 = ixp
         arg = 0

      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------

         if (ixm /= 99) then
         ! currently in first half of data interval
            if (ixx <= 1) then
               if (yr > fyear_init) then ! get data from previous year
                  call file_year (data_file, yr-1)
               else             ! yr = fyear_init, no prior data exists
                  if (maxrec > 12) then ! extrapolate from first record
                     if (ixx == 1) n2 = ixx
                  else          ! go to end of fyear_final
                     call file_year (data_file, fyear_final)
                  endif
               endif            ! yr > fyear_init
            endif               ! ixx <= 1

            call ice_open_nc (data_file, fid)

            arg = 1
            nrec = recd + n2

            call ice_read_nc & 
                 (fid, nrec, fieldname, field_data(:,:,arg,:), dbug, &
                  field_loc, field_type)

            if (ixx==1) call ice_close_nc(fid)
         endif                  ! ixm ne 99

         ! always read ixx data from data file for current year
         call file_year (data_file, yr)
         call ice_open_nc (data_file, fid)

         arg = arg + 1
         nrec = recd + ixx

         call ice_read_nc & 
              (fid, nrec, fieldname, field_data(:,:,arg,:), dbug, &
               field_loc, field_type)

         if (ixp /= 99) then
         ! currently in latter half of data interval
            if (ixx==maxrec) then
               if (yr < fyear_final) then ! get data from following year
                  call ice_close_nc(fid)
                  call file_year (data_file, yr+1)
                  call ice_open_nc (data_file, fid)
               else             ! yr = fyear_final, no more data exists
                  if (maxrec > 12) then ! extrapolate from ixx
                     n4 = ixx
                  else          ! go to beginning of fyear_init
                     call ice_close_nc(fid)
                     call file_year (data_file, fyear_init)
                     call ice_open_nc (data_file, fid)

                  endif
               endif            ! yr < fyear_final
            endif               ! ixx = maxrec

            arg = arg + 1
            nrec = recd + n4

            call ice_read_nc & 
                 (fid, nrec, fieldname, field_data(:,:,arg,:), dbug, &
                  field_loc, field_type)
         endif                  ! ixp /= 99

         call ice_close_nc(fid)

      endif                     ! flag

      call ice_timer_stop(timer_readwrite)  ! reading/writing

#else
      field_data = c0 ! to satisfy intent(out) attribute
#endif
      end subroutine read_data_nc

!=======================================================================
!
!BOP
!
! !IROUTINE: read_clim_data - read annual climatological data
!
! !INTERFACE:
!
      subroutine read_clim_data (readflag, recd, ixm, ixx, ixp, &
                                 data_file, field_data, &
                                 field_loc, field_type)
!
! !DESCRIPTION:
!
! Read data needed for interpolation, as in read_data.
! Assume a one-year cycle of climatological data, so that there is
!  no need to get data from other years or to extrapolate data beyond
!  the forcing time period.
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      use ice_read_write
      use ice_diagnostics, only: check_step
!
! !INPUT/OUTPUT PARAMETERS:
!
      logical (kind=log_kind),intent(in) :: readflag

      integer (kind=int_kind), intent(in) :: &
        recd            , & ! baseline record number
        ixm,ixx,ixp         ! record numbers of 3 data values
                            ! relative to recd

      character (char_len_long), intent(in) ::  data_file

      integer (kind=int_kind), intent(in) :: &
           field_loc, &      ! location of field on staggered grid
           field_type        ! type of field (scalar, vector, angle)

      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), &
        intent(out) :: &
        field_data         ! 2 values needed for interpolation
!
!EOP
!
      integer (kind=int_kind) :: &
        nbits          , & ! = 32 for single precision, 64 for double
        nrec           , & ! record number to read
        arg                ! value of time argument in field_data

      call ice_timer_start(timer_readwrite)  ! reading/writing

      nbits = 64                ! double precision data

      if (istep1 > check_step) dbug = .true.  !! debugging

      if (my_task==master_task .and. (dbug)) &
        write(nu_diag,*) '  ', trim(data_file)

      if (readflag) then

      !-----------------------------------------------------------------
      ! read data
      !-----------------------------------------------------------------

         call ice_open (nu_forcing, data_file, nbits)

         arg = 0
         if (ixm /= 99) then
            arg = 1
            nrec = recd + ixm
            call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
                           'rda8', dbug, field_loc, field_type)
         endif

         arg = arg + 1
         nrec = recd + ixx
         call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
                        'rda8', dbug, field_loc, field_type)

         if (ixp /= 99) then
            arg = arg + 1
            nrec = recd + ixp
            call ice_read (nu_forcing, nrec, field_data(:,:,arg,:), &
                           'rda8', dbug, field_loc, field_type)
         endif

         if (my_task == master_task) close (nu_forcing)
      endif                     ! readflag

      call ice_timer_stop(timer_readwrite)  ! reading/writing

      end subroutine read_clim_data

!=======================================================================
!
!BOP
!
! !IROUTINE: interp_coeff_monthly - Compute monthly data interpolation coefficients
!
! !INTERFACE:
!
      subroutine interp_coeff_monthly (recslot)
!
! !DESCRIPTION:
!
! Compute coefficients for interpolating monthly data to current time step.
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
          recslot         ! slot (1 or 2) for current record
!
!EOP
!
      real (kind=dbl_kind) :: &
          tt           , & ! seconds elapsed in current year
          t1, t2           ! seconds elapsed at month midpoint

      real (kind=dbl_kind) :: &
          daymid(0:13)     ! month mid-points

      daymid(1:13) = 14._dbl_kind   ! time frame ends 0 sec into day 15
      daymid(0)    = 14._dbl_kind - daymo(12)  ! Dec 15, 0 sec

      ! make time cyclic
      tt = mod(ftime/secday,dayyr)

      ! Find neighboring times

      if (recslot==2) then      ! first half of month
        t2 = daycal(month) + daymid(month)   ! midpoint, current month
        if (month == 1) then
          t1 = daymid(0)                 ! Dec 15 (0 sec)
        else
          t1 = daycal(month-1) + daymid(month-1) ! midpoint, previous month
        endif
      else                      ! second half of month
        t1 = daycal(month) + daymid(month)    ! midpoint, current month
        t2 = daycal(month+1) + daymid(month+1)! day 15 of next month (0 sec)
      endif

      ! Compute coefficients
      c1intp = (t2 - tt) / (t2 - t1)
      c2intp =  c1 - c1intp

      end subroutine interp_coeff_monthly

!=======================================================================
!
!BOP
!
! !IROUTINE: interp_coeff
!
! !INTERFACE:
!
      subroutine interp_coeff (recnum, recslot, secint, dataloc)
!
! !DESCRIPTION:
!
! Compute coefficients for interpolating data to current time step.
! Works for any data interval that divides evenly into a
!  year (daily, 6-hourly, etc.)
! Use interp_coef_monthly for monthly data.
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      integer (kind=int_kind), intent(in) :: &
          recnum      , & ! record number for current data value
          recslot     , & ! spline slot for current record
          dataloc         ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval

      real (kind=dbl_kind), intent(in) :: &
          secint                    ! seconds in data interval
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      real (kind=dbl_kind) :: &
          secyr            ! seconds in a year

      real (kind=dbl_kind) :: &
          tt           , & ! seconds elapsed in current year
          t1, t2       , & ! seconds elapsed at data points
          rcnum            ! recnum => dbl_kind

      secyr = dayyr * secday         ! seconds in a year
      tt = mod(ftime,secyr)

      ! Find neighboring times
      rcnum = real(recnum,kind=dbl_kind)
      if (recslot==2) then           ! current record goes in slot 2
         if (dataloc==1) then        ! data located at middle of interval
            t2 = (rcnum-p5)*secint
         else                        !  data located at end of interval
            t2 = rcnum*secint
         endif
         t1 = t2 - secint            !  - 1 interval
      else                           ! recslot = 1
         if (dataloc==1) then        ! data located at middle of interval
            t1 = (rcnum-p5)*secint
         else                        
            t1 = rcnum*secint        ! data located at end of interval
         endif
         t2 = t1 + secint            !  + 1 interval
      endif

      ! Compute coefficients
      c1intp =  abs((t2 - tt) / (t2 - t1))
      c2intp =  c1 - c1intp

      end subroutine interp_coeff

!=======================================================================
!BOP
!
! !IROUTINE: interpolate_data
!
! !INTERFACE:
!
      subroutine interpolate_data (field_data, field)
!
! !DESCRIPTION:
!
! Linear interpolation
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
      use ice_domain, only: nblocks
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), dimension(nx_block,ny_block,2,max_blocks), &
        intent(in) :: &
        field_data    ! 2 values used for interpolation

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks), &
        intent(out) :: &
        field         ! interpolated field
!
!EOP
!
      integer (kind=int_kind) :: i,j, iblk

      do iblk = 1, nblocks
         do j = 1, ny_block
         do i = 1, nx_block
            field(i,j,iblk) = c1intp * field_data(i,j,1,iblk) &
                            + c2intp * field_data(i,j,2,iblk)
         enddo
         enddo
      enddo

      end subroutine interpolate_data

!=======================================================================
!
!BOP
!
! !IROUTINE: file_year - construct name of atmospheric data file
!
! !INTERFACE:
!
      subroutine file_year (data_file, yr)
!
! !DESCRIPTION:
!
! Construct the correct name of the atmospheric data file
! to be read, given the year and assuming the naming convention
! that filenames end with 'yyyy.dat' or 'yyyy.r' or 'yyyy.nc'.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      character (char_len_long), intent(inout) ::  data_file
!
!EOP
!
      integer (kind=int_kind), intent(in) :: yr

      character (char_len_long) :: tmpname

      integer (kind=int_kind) :: i

      if (trim(atm_data_type) == 'ecmwf') then ! NPS/ECMWF naming convention
         i = index(data_file,'.r') - 5
         tmpname = data_file
         write(data_file,'(a,i4.4,a)') tmpname(1:i), yr, '.r'
      elseif (trim(atm_data_type) == 'hadgem') then ! netcdf
         i = index(data_file,'.nc') - 5
         tmpname = data_file
         write(data_file,'(a,i4.4,a)') tmpname(1:i), yr, '.nc'
      else                                     ! LANL/NCAR naming convention
         i = index(data_file,'.dat') - 5
         tmpname = data_file
         write(data_file,'(a,i4.4,a)') tmpname(1:i), yr, '.dat'
      endif

      end subroutine file_year

!=======================================================================
!
!BOP
!
! !IROUTINE: prepare_forcing - finish manipulating forcing
!
! !INTERFACE:
!
      subroutine prepare_forcing (nx_block, ny_block, &
                                  ilo, ihi, jlo, jhi, &
                                  hm,                 &
                                  Tair,     fsw,      &    
                                  cldf,     flw,      &
                                  frain,    fsnow,    &
                                  Qa,       rhoa,     &
                                  uatm,     vatm,     &
                                  strax,    stray,    &
                                  zlvl,     wind,     &
                                  swvdr,    swvdf,    &
                                  swidr,    swidf,    &
                                  potT,     ANGLET,   &
                                  Tsfc,     sst,      &
                                  aice)
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         Tair    , & ! air temperature  (K)
         ANGLET  , & ! ANGLE converted to T-cells
         Tsfc    , & ! ice skin temperature
         sst     , & ! sea surface temperature
         aice    , & ! ice area fraction
         hm          ! land mask
     
      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: &
         fsw     , & ! incoming shortwave radiation (W/m^2)
         cldf    , & ! cloud fraction
         frain   , & ! rainfall rate (kg/m^2 s)
         fsnow   , & ! snowfall rate (kg/m^2 s)
         Qa      , & ! specific humidity (kg/kg)
         rhoa    , & ! air density (kg/m^3)
         uatm    , & ! wind velocity components (m/s)
         vatm    , &
         strax   , & ! wind stress components (N/m^2)
         stray   , &
         zlvl    , & ! atm level height (m)
         wind    , & ! wind speed (m/s)
         flw     , & ! incoming longwave radiation (W/m^2)
         swvdr   , & ! sw down, visible, direct  (W/m^2)
         swvdf   , & ! sw down, visible, diffuse (W/m^2)
         swidr   , & ! sw down, near IR, direct  (W/m^2)
         swidf   , & ! sw down, near IR, diffuse (W/m^2)
         potT        ! air potential temperature  (K)
!
!EOP
!
      integer (kind=int_kind) :: &
         i, j

      real (kind=dbl_kind) :: workx, worky, &
         fcc, sstk, rtea, ptem, qlwm, precip_factor

      do j = jlo, jhi
      do i = ilo, ihi

      !-----------------------------------------------------------------
      ! make sure interpolated values are physically realistic
      !-----------------------------------------------------------------
         cldf (i,j) = max(min(cldf(i,j),c1),c0)
         fsw  (i,j) = max(fsw(i,j),c0)
         fsnow(i,j) = max(fsnow(i,j),c0)
         rhoa (i,j) = max(rhoa(i,j),c0)
         Qa   (i,j) = max(Qa(i,j),c0)

      enddo                     ! i
      enddo                     ! j

      !-----------------------------------------------------------------
      ! calculations specific to datasets
      !-----------------------------------------------------------------

      if (trim(atm_data_type) == 'ncar') then
         do j = jlo, jhi
         do i = ilo, ihi

      !-----------------------------------------------------------------
      ! correct known biases in NCAR data (as in CCSM latm)
      !-----------------------------------------------------------------

            Qa (i,j) = Qa (i,j) * 0.94_dbl_kind
            fsw(i,j) = fsw(i,j) * 0.92_dbl_kind

      !-----------------------------------------------------------------
      ! compute downward longwave as in Parkinson and Washington (1979)
      ! (for now)
      !-----------------------------------------------------------------

            flw(i,j) = stefan_boltzmann*Tair(i,j)**4 &
                     * (c1 - 0.261_dbl_kind &
                      *exp(-7.77e-4_dbl_kind*(Tffresh - Tair(i,j))**2)) &
                     * (c1 + 0.275_dbl_kind*cldf(i,j))

      ! precip is in mm/month; converted to mks below

         enddo
         enddo

      elseif (trim(atm_data_type) == 'ecmwf') then
         do j = jlo, jhi
         do i = ilo, ihi

      !-----------------------------------------------------------------
      ! The following assumes that the input Qa is really dew point temp
      ! (deg K) and need to be converted to specific humidity (kg/kg).
      ! Cf. ice_atmo module.
      !-----------------------------------------------------------------

         Qa (i,j) = (qqqocn/rhoa(i,j)) * exp(-TTTocn/Qa(i,j))
         Qa (i,j) = max(Qa(i,j),c0)

      !-----------------------------------------------------------------
      ! compute downward longwave as in Parkinson and Washington (1979)
      ! (for now)
      !-----------------------------------------------------------------

         flw(i,j) = stefan_boltzmann*Tair(i,j)**4 &
                  * (c1 - 0.261_dbl_kind &
                   *exp(-7.77e-4_dbl_kind*(Tffresh - Tair(i,j))**2)) &
                  * (c1 + 0.275_dbl_kind*cldf(i,j))

      ! precip is in mm/month; converted to mks below

         enddo
         enddo

      elseif (trim(atm_data_type) == 'LYq') then

      !-----------------------------------------------------------------
      ! longwave, Rosati and Miyakoda, JPO 18, p. 1607 (1988) - sort of 
      !-----------------------------------------------------------------

         do j = jlo, jhi
         do i = ilo, ihi
            fcc = c1 - 0.8_dbl_kind * cldf(i,j)
            sstk = (Tsfc(i,j) * aice(i,j) &
                 + sst(i,j) * (c1 - aice(i,j))) + Tffresh
            rtea = sqrt(c1000*Qa(i,j) /  &
                  (0.622_dbl_kind+0.378_dbl_kind*Qa(i,j)))
            ptem = Tair(i,j)    ! get this from stability?
            qlwm = ptem * ptem * ptem  &
                 * ( ptem*(0.39_dbl_kind-0.05_dbl_kind*rtea)*fcc  &
                                      + c4*(sstk-ptem) )
            flw(i,j) = emissivity*stefan_boltzmann * ( sstk**4 - qlwm )
            flw(i,j) = flw(i,j) * hm(i,j) ! land mask
         enddo
         enddo

      endif                     ! atm_data_type

      !-----------------------------------------------------------------
      ! Compute other fields needed by model
      !-----------------------------------------------------------------

      ! convert precipitation units to kg/m^2 s
      if (trim(precip_units) == 'mm_per_month') then
         precip_factor = c12/(secday*days_per_year) 
      elseif (trim(precip_units) == 'mm_per_day') then
         precip_factor = c1/secday
      elseif (trim(precip_units) == 'mm_per_sec' .or. &
              trim(precip_units) == 'mks') then 
         precip_factor = c1    ! mm/sec = kg/m^2 s
      endif

      do j = jlo, jhi
      do i = ilo, ihi

         zlvl(i,j) = c10
         potT(i,j) = Tair(i,j)

        ! divide shortwave into spectral bands
         swvdr(i,j) = fsw(i,j)*frcvdr        ! visible direct
         swvdf(i,j) = fsw(i,j)*frcvdf        ! visible diffuse
         swidr(i,j) = fsw(i,j)*frcidr        ! near IR direct
         swidf(i,j) = fsw(i,j)*frcidf        ! near IR diffuse
                 
        ! convert precipitation units to kg/m^2 s
         fsnow(i,j) = fsnow(i,j) * precip_factor

      enddo                     ! i
      enddo                     ! j

      ! determine whether precip is rain or snow
      ! HadGEM forcing provides separate snowfall and rainfall rather 
      ! than total precipitation
      if (trim(atm_data_type) /= 'hadgem') then

        do j = jlo, jhi
        do i = ilo, ihi
           frain(i,j) = c0                     
           if (Tair(i,j) >= Tffresh) then
               frain(i,j) = fsnow(i,j)
               fsnow(i,j) = c0
           endif
        enddo                     ! i
        enddo                     ! j

      endif

      if (calc_strair) then

        do j = jlo, jhi
        do i = ilo, ihi

            wind(i,j) = sqrt(uatm(i,j)**2 + vatm(i,j)**2)

      !-----------------------------------------------------------------
      ! Rotate zonal/meridional vectors to local coordinates.
      ! Velocity comes in on T grid, but is oriented geographically ---
      ! need to rotate to pop-grid FIRST using ANGLET
      ! then interpolate to the U-cell centers  (otherwise we
      ! interpolate across the pole).
      ! Use ANGLET which is on the T grid !
      ! Atmo variables are needed in T cell centers in subroutine 
      ! atmo_boundary_layer, and are interpolated to the U grid later as 
      ! necessary.
      !-----------------------------------------------------------------
           workx      = uatm(i,j) ! wind velocity, m/s
           worky      = vatm(i,j)
           uatm (i,j) = workx*cos(ANGLET(i,j)) & ! convert to POP grid
                      + worky*sin(ANGLET(i,j))   ! note uatm, vatm, wind
           vatm (i,j) = worky*cos(ANGLET(i,j)) & !  are on the T-grid here
                      - workx*sin(ANGLET(i,j))

        enddo                     ! i
        enddo                     ! j

      else  ! strax, stray, wind are read from files

        do j = jlo, jhi
        do i = ilo, ihi

           workx      = strax(i,j) ! wind stress
           worky      = stray(i,j)
           strax(i,j) = workx*cos(ANGLET(i,j)) & ! convert to POP grid
                      + worky*sin(ANGLET(i,j))   ! note strax, stray, wind
           stray(i,j) = worky*cos(ANGLET(i,j)) & !  are on the T-grid here
                      - workx*sin(ANGLET(i,j))

        enddo                     ! i
        enddo                     ! j

      endif                   ! calc_strair

      end subroutine prepare_forcing

!=======================================================================
! NCAR atmospheric forcing
!=======================================================================
!
!BOP
!
! !IROUTINE: ncar_files - construct filenames for NCAR bulk atmospheric data
!
! !INTERFACE:
!
      subroutine ncar_files (yr)
!
! !DESCRIPTION:
!
! Construct filenames based on the LANL naming conventions for NCAR data.
! Edit for other directory structures or filenames.
! Note: The year number in these filenames does not matter, because
!       subroutine file\_year will insert the correct year.
!
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           yr                   ! current forcing year
!
!EOP
!
      fsw_file = &
           trim(atm_data_dir)//'ISCCPM/MONTHLY/RADFLX/swdn.1996.dat'
      call file_year(fsw_file,yr)

      flw_file = &
           trim(atm_data_dir)//'ISCCPM/MONTHLY/RADFLX/cldf.1996.dat'
      call file_year(flw_file,yr)

      rain_file = &
           trim(atm_data_dir)//'MXA/MONTHLY/PRECIP/prec.1996.dat'
      call file_year(rain_file,yr)

      uwind_file = &
           trim(atm_data_dir)//'NCEP/4XDAILY/STATES/u_10.1996.dat'
      call file_year(uwind_file,yr)

      vwind_file = &
           trim(atm_data_dir)//'NCEP/4XDAILY/STATES/v_10.1996.dat'
      call file_year(vwind_file,yr)

      tair_file = &
           trim(atm_data_dir)//'NCEP/4XDAILY/STATES/t_10.1996.dat'
      call file_year(tair_file,yr)

      humid_file = &
           trim(atm_data_dir)//'NCEP/4XDAILY/STATES/q_10.1996.dat'
      call file_year(humid_file,yr)

      rhoa_file = &
           trim(atm_data_dir)//'NCEP/4XDAILY/STATES/dn10.1996.dat'
      call file_year(rhoa_file,yr)

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Forcing data year =', fyear
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(fsw_file)
         write (nu_diag,*) trim(flw_file)
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(uwind_file)
         write (nu_diag,*) trim(vwind_file)
         write (nu_diag,*) trim(tair_file)
         write (nu_diag,*) trim(humid_file)
         write (nu_diag,*) trim(rhoa_file)
      endif                     ! master_task

      end subroutine ncar_files

!=======================================================================
!
!BOP
!
! !IROUTINE: ncar_data - read NCAR bulk atmospheric data
!
! !INTERFACE:
!
      subroutine ncar_data
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      use ice_flux
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      integer (kind=int_kind) :: &
          i, j        , &
          ixm,ixx,ixp , & ! record numbers for neighboring months
          recnum      , & ! record number
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          dataloc     , & ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
          midmonth        ! middle day of month

      real (kind=dbl_kind) :: &
          sec6hr              ! number of seconds in 6 hours

      logical (kind=log_kind) :: readm, read6

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

      if (trim(atm_data_format) == 'bin') then
         call read_data (readm, 0, fyear, ixm, month, ixp, &
                         maxrec, fsw_file, fsw_data, &
                         field_loc_center, field_type_scalar)
         call read_data (readm, 0, fyear, ixm, month, ixp, &
                         maxrec, flw_file, cldf_data, &
                         field_loc_center, field_type_scalar)
         call read_data (readm, 0, fyear, ixm, month, ixp, &
                         maxrec, rain_file, fsnow_data, &
                         field_loc_center, field_type_scalar)
      else
         call abort_ice ('nonbinary atm_data_format unavailable')
!        The routine exists, for example:  
!         call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
!                            maxrec, fsw_file, 'fsw', fsw_data, &
!                            field_loc_center, field_type_scalar)
!         call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
!                            maxrec, flw_file, 'cldf',cldf_data, &
!                            field_loc_center, field_type_scalar)
!         call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
!                            maxrec, rain_file,'prec',fsnow_data, &
!                            field_loc_center, field_type_scalar)
      endif

      ! Interpolate to current time step
      call interpolate_data (fsw_data,   fsw)
      call interpolate_data (cldf_data,  cldf)
      call interpolate_data (fsnow_data, fsnow)

    !-------------------------------------------------------------------
    ! 6-hourly data
    !
    ! Assume that the 6-hourly value is located at the end of the
    !  6-hour period.  This is the convention for NCEP reanalysis data.
    !  E.g. record 1 gives conditions at 6 am GMT on 1 January.
    !-------------------------------------------------------------------

      dataloc = 2               ! data located at end of interval
      sec6hr = secday/c4        ! seconds in 6 hours
      maxrec = 1460             ! 365*4

      ! current record number
      recnum = 4*int(yday) - 3 + int(real(sec,kind=dbl_kind)/sec6hr)

      ! Compute record numbers for surrounding data

      ixm = mod(recnum+maxrec-2,maxrec) + 1
      ixx = mod(recnum-1,       maxrec) + 1
!      ixp = mod(recnum,         maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record always goes in slot 2.

      recslot = 2
      ixp = 99
      call interp_coeff (recnum, recslot, sec6hr, dataloc)

      ! Read
      read6 = .false.
      if (istep==1 .or. oldrecnum /= recnum) read6 = .true.

      if (trim(atm_data_format) == 'bin') then
         call read_data (read6, 0, fyear, ixm, ixx, ixp, &
                         maxrec, tair_file, Tair_data, &
                         field_loc_center, field_type_scalar)
         call read_data (read6, 0, fyear, ixm, ixx, ixp, &
                         maxrec, uwind_file, uatm_data, &
                         field_loc_center, field_type_scalar)
         call read_data (read6, 0, fyear, ixm, ixx, ixp, &
                         maxrec, vwind_file, vatm_data, &
                         field_loc_center, field_type_scalar)
         call read_data (read6, 0, fyear, ixm, ixx, ixp, &
                         maxrec, rhoa_file, rhoa_data, &
                         field_loc_center, field_type_scalar)
         call read_data (read6, 0, fyear, ixm, ixx, ixp, &
                         maxrec, humid_file, Qa_data, &
                         field_loc_center, field_type_scalar)
      else
         call abort_ice ('nonbinary atm_data_format unavailable')
      endif

      ! Interpolate
      call interpolate_data (Tair_data, Tair)
      call interpolate_data (uatm_data, uatm)
      call interpolate_data (vatm_data, vatm)
      call interpolate_data (rhoa_data, rhoa)
      call interpolate_data (Qa_data,   Qa)

      ! Save record number for next time step
      oldrecnum = recnum

      end subroutine ncar_data

!=======================================================================
! ECMWF atmospheric forcing
!=======================================================================
!BOP
!
! !IROUTINE: ecmwf_files - construct filenames for ECMWF atmospheric data
!
! !INTERFACE:
!
      subroutine ecmwf_files (yr)
!
! !DESCRIPTION:
!
! Construct filenames based on naming conventions used by Wieslaw Maslowski
!  for reading (mostly) ECMWF atmospheric data. 
! Edit for other directory structures or filenames.
! Note: The year number in these filenames does not matter, because
!       subroutine file\_year will insert the correct year.
!
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL (based on ncar_files)
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           yr                   ! current forcing year
!
!EOP
!
      fsw_file = &
           trim(atm_data_dir)//'sol_2002.r'
      call file_year(fsw_file,yr)

      flw_file = &
           trim(atm_data_dir)//'flo_2002.r'
      call file_year(flw_file,yr)

      rain_file = &
           trim(atm_data_dir)//'prec_lanl_12.r'
      ! Comment out the file_year call if rain file is from climatology
!!!      call file_year(rain_file,yr)
 
      uwind_file = &
           trim(atm_data_dir)//'ucmp_2002.r'
      call file_year(uwind_file,yr)

      vwind_file = &
           trim(atm_data_dir)//'vcmp_2002.r'
      call file_year(vwind_file,yr)

      tair_file = &
           trim(atm_data_dir)//'tair_2002.r'
      call file_year(tair_file,yr)

      humid_file = &
           trim(atm_data_dir)//'qa_2002.r'
      call file_year(humid_file,yr)

      rhoa_file = &
           trim(atm_data_dir)//'rhoa_ncar85-88_12.r'
      ! Comment out the file_year call if rhoa file is from climatology
!!!      call file_year(rhoa_file,yr)

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Forcing data year = ', fyear
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(fsw_file)
         write (nu_diag,*) trim(flw_file)
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(uwind_file)
         write (nu_diag,*) trim(vwind_file)
         write (nu_diag,*) trim(tair_file)
         write (nu_diag,*) trim(humid_file)
         write (nu_diag,*) trim(rhoa_file)
      endif                     ! master_task

      end subroutine ecmwf_files

!=======================================================================
!BOP
!
! !IROUTINE: ecmwf_data - read ECMWF atmospheric data
!
! !INTERFACE:
!
      subroutine ECMWF_data
!
! !DESCRIPTION:
!
! Read ECMWF atmospheric data.
!
! !REVISION HISTORY:
!
! authors: William H. Lipscomb, LANL
!          Wieslaw Maslowski, NPS
!          Based on ncar_data
!
! !USES:
!
      use ice_flux
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      integer (kind=int_kind) :: &
          ixm,ixx,ixp , & ! record numbers for neighboring months
          recnum      , & ! record number
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          dataloc     , & ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
          midmonth        ! middle day of month

      logical (kind=log_kind) :: readm, readd

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

      if (trim(atm_data_format) == 'bin') then
         call read_clim_data (readm, 0,  ixm, month, ixp, &
              rhoa_file, rhoa_data, field_loc_center, field_type_scalar)
         call read_clim_data (readm, 0,  ixm, month, ixp, &
              rain_file, fsnow_data, field_loc_center, field_type_scalar)
      else
         call abort_ice ('nonbinary atm_data_format unavailable')
      endif

      ! Interpolate to current time step
      call interpolate_data (fsnow_data, fsnow)
      call interpolate_data (rhoa_data, rhoa)  

    !-------------------------------------------------------------------
    ! Daily data
    !
    ! Assume that the daily value is located in the middle of the
    !  24-hour period.
    !-------------------------------------------------------------------

      dataloc = 1          ! data located in middle of interval
      maxrec = 365         ! days in a year (no leap years)

      ! current record number
      recnum = int(yday)   ! current record number

      ! Compute record numbers for surrounding data

      ixm = mod(recnum+maxrec-2,maxrec) + 1
      ixx = mod(recnum-1,       maxrec) + 1
      ixp = mod(recnum,         maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record goes in slot 2

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      if (real(sec,kind=dbl_kind) < p5*secday-puny) then  ! first half of day
         recslot = 2
         ixp = 99
      else                             ! second half of day
         recslot = 1
         ixm = 99
      endif

      call interp_coeff (recnum, recslot, secday, dataloc)

      ! Read new data at midpoint of day

      readd = .false.
      if (istep==1 .or. (recslot==1 .and. oldrecslot==2)) &
           readd = .true.

      if (trim(atm_data_format) == 'bin') then
         call read_data (readd, 0, fyear, ixm, ixx, ixp, maxrec, &
                         tair_file, Tair_data, &
                         field_loc_center, field_type_scalar)
         call read_data (readd, 0, fyear, ixm, ixx, ixp, maxrec, &
                         uwind_file, uatm_data, &
                         field_loc_center, field_type_scalar)
         call read_data (readd, 0, fyear, ixm, ixx, ixp, maxrec, &
                         vwind_file, vatm_data, &
                         field_loc_center, field_type_scalar)
         call read_data (readd, 0, fyear, ixm, ixx, ixp, maxrec, &
                         fsw_file, fsw_data, &
                         field_loc_center, field_type_scalar)
         call read_data (readd, 0, fyear, ixm, ixx, ixp, maxrec, &
                         flw_file, flw_data, &
                         field_loc_center, field_type_scalar)
         call read_data (readd, 0, fyear, ixm, ixx, ixp, maxrec, &
                         humid_file, Qa_data, &
                         field_loc_center, field_type_scalar)
      else
         call abort_ice ('nonbinary atm_data_format unavailable')
      endif

      ! Interpolate
      call interpolate_data (Tair_data, Tair)
      call interpolate_data (uatm_data, uatm)
      call interpolate_data (vatm_data, vatm)
      call interpolate_data ( fsw_data, fsw)
      call interpolate_data ( flw_data, flw)
      call interpolate_data (  Qa_data, Qa)

      ! Save recslot for next time step
      oldrecslot = recslot

      end subroutine ecmwf_data

!=======================================================================
! Large and Yeager forcing (AOMIP style)
!=======================================================================
!
!BOP
!
! !IROUTINE: LY_files - construct filenames for Large and Yeager data
!     note:  includes AOMIP (OMIP) cldf climatology 
!
! !INTERFACE:
!
      subroutine LY_files (yr)
!
! !DESCRIPTION:
!
! Construct filenames based on the LANL naming conventions for NCAR data.
! Edit for other directory structures or filenames.
! Note: The year number in these filenames does not matter, because
!       subroutine file_year will insert the correct year.
!
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           yr                   ! current forcing year
!
!EOP
!
      flw_file = &
           trim(atm_data_dir)//'MONTHLY/cldf.omip.dat'

      rain_file = &
           trim(atm_data_dir)//'MONTHLY/prec.nmyr.dat'

      uwind_file = &
           trim(atm_data_dir)//'4XDAILY/u_10.1996.dat'
      call file_year(uwind_file,yr)

      vwind_file = &
           trim(atm_data_dir)//'4XDAILY/v_10.1996.dat'
      call file_year(vwind_file,yr)

      tair_file = &
           trim(atm_data_dir)//'4XDAILY/t_10.1996.dat'
      call file_year(tair_file,yr)

      humid_file = &
           trim(atm_data_dir)//'4XDAILY/q_10.1996.dat'
      call file_year(humid_file,yr)

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Forcing data year = ', fyear         
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(flw_file)
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(uwind_file)
         write (nu_diag,*) trim(vwind_file)
         write (nu_diag,*) trim(tair_file)
         write (nu_diag,*) trim(humid_file)
      endif                     ! master_task

      end subroutine LY_files

!=======================================================================
!
!BOP
!
! !IROUTINE: LY_data - read Large and Yeager atmospheric data
!        note:  also uses AOMIP protocol, in part
!
! !INTERFACE:
!
      subroutine LY_data
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      use ice_global_reductions
      use ice_domain, only: nblocks, distrb_info, blocks_ice
      use ice_flux 
      use ice_grid, only: hm, tlon, tlat, tmask, umask
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: & 
          i, j        , &
          imx,ixx,ipx , & ! record numbers for neighboring months
          recnum      , & ! record number
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          midmonth    , & ! middle day of month
          dataloc     , & ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
          iblk        , & ! block index
          ilo,ihi,jlo,jhi ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
          sec6hr          , & ! number of seconds in 6 hours
          vmin, vmax

      logical (kind=log_kind) :: readm, read6

      type (block) :: &
         this_block           ! block information for current block

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
      imx  = mod(month+maxrec-2,maxrec) + 1
      ipx  = mod(month,         maxrec) + 1
      if (mday >= midmonth) imx = 99  ! other two points will be used
      if (mday <  midmonth) ipx = 99

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

      call read_clim_data (readm, 0, imx, month, ipx,  &
             flw_file, cldf_data, field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, imx, month, ipx,  &
             rain_file, fsnow_data, field_loc_center, field_type_scalar)

      call interpolate_data (cldf_data, cldf)
      call interpolate_data (fsnow_data, fsnow)  ! units mm/s = kg/m^2/s

    !-------------------------------------------------------------------
    ! 6-hourly data
    ! 
    ! Assume that the 6-hourly value is located at the end of the
    !  6-hour period.  This is the convention for NCEP reanalysis data.
    !  E.g. record 1 gives conditions at 6 am GMT on 1 January.
    !-------------------------------------------------------------------

      dataloc = 2               ! data located at end of interval
      sec6hr = secday/c4        ! seconds in 6 hours
      maxrec = 1460             ! 365*4

      ! current record number
      recnum = 4*int(yday) - 3 + int(real(sec,kind=dbl_kind)/sec6hr)

      ! Compute record numbers for surrounding data (2 on each side)

      imx = mod(recnum+maxrec-2,maxrec) + 1
      ixx = mod(recnum-1,       maxrec) + 1
!     ipx = mod(recnum,         maxrec) + 1

      ! Compute interpolation coefficients
      ! If data is located at the end of the time interval, then the
      !  data value for the current record goes in slot 2

      recslot = 2
      ipx = 99
      call interp_coeff (recnum, recslot, sec6hr, dataloc)

      ! Read
      read6 = .false.
      if (istep==1 .or. oldrecnum .ne. recnum) read6 = .true.

      if (trim(atm_data_format) == 'bin') then
         call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec, &
                         tair_file, Tair_data, &
                         field_loc_center, field_type_scalar)
         call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec, &
                         uwind_file, uatm_data, &
                         field_loc_center, field_type_scalar)
         call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec, &
                         vwind_file, vatm_data, &
                         field_loc_center, field_type_scalar)
         call read_data (read6, 0, fyear, imx, ixx, ipx, maxrec, &
                         humid_file, Qa_data, &
                         field_loc_center, field_type_scalar)
      else
         call abort_ice ('nonbinary atm_data_format unavailable')
      endif

      ! Interpolate
      call interpolate_data (Tair_data, Tair)
      call interpolate_data (uatm_data, uatm)
      call interpolate_data (vatm_data, vatm)
      call interpolate_data (Qa_data, Qa)

      do iblk = 1, nblocks
        call Qa_fixLY(nx_block,  ny_block, &
                                 Tair (:,:,iblk), &
                                 Qa   (:,:,iblk))

        do j = 1, ny_block
          do i = 1, nx_block
            Qa  (i,j,iblk) = Qa  (i,j,iblk) * hm(i,j,iblk)
            Tair(i,j,iblk) = Tair(i,j,iblk) * hm(i,j,iblk)
            uatm(i,j,iblk) = uatm(i,j,iblk) * hm(i,j,iblk)
            vatm(i,j,iblk) = vatm(i,j,iblk) * hm(i,j,iblk)
          enddo
        enddo

      ! AOMIP
        this_block = get_block(blocks_ice(iblk),iblk)         
        ilo = this_block%ilo
        ihi = this_block%ihi
        jlo = this_block%jlo
        jhi = this_block%jhi

        call compute_shortwave(nx_block, ny_block, &
                               ilo, ihi, jlo, jhi, &
                               TLON (:,:,iblk), &
                               TLAT (:,:,iblk), &
                               hm   (:,:,iblk), &
                               Qa   (:,:,iblk), &
                               cldf (:,:,iblk), &
                               fsw  (:,:,iblk))

      enddo  ! iblk

      ! Save record number
      oldrecnum = recnum

         if (dbug) then
           if (my_task == master_task) write (nu_diag,*) 'LY_bulk_data'
           vmin = global_minval(fsw,distrb_info,tmask)
                               
           vmax = global_maxval(fsw,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'fsw',vmin,vmax 
           vmin = global_minval(cldf,distrb_info,tmask)
           vmax = global_maxval(cldf,distrb_info,tmask)
           if (my_task.eq.master_task) & 
               write (nu_diag,*) 'cldf',vmin,vmax
           vmin =global_minval(fsnow,distrb_info,tmask)
           vmax =global_maxval(fsnow,distrb_info,tmask)
           if (my_task.eq.master_task) & 
               write (nu_diag,*) 'fsnow',vmin,vmax
           vmin = global_minval(Tair,distrb_info,tmask)
           vmax = global_maxval(Tair,distrb_info,tmask)
           if (my_task.eq.master_task) & 
               write (nu_diag,*) 'Tair',vmin,vmax
           vmin = global_minval(uatm,distrb_info,umask)
           vmax = global_maxval(uatm,distrb_info,umask)
           if (my_task.eq.master_task) & 
               write (nu_diag,*) 'uatm',vmin,vmax
           vmin = global_minval(vatm,distrb_info,umask)
           vmax = global_maxval(vatm,distrb_info,umask)
           if (my_task.eq.master_task) & 
               write (nu_diag,*) 'vatm',vmin,vmax
           vmin = global_minval(Qa,distrb_info,tmask)
           vmax = global_maxval(Qa,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'Qa',vmin,vmax

        endif                   ! dbug

      end subroutine LY_data

!=======================================================================

      subroutine compute_shortwave(nx_block,  ny_block, &
                                   ilo, ihi, jlo, jhi, &
                                   TLON, TLAT, hm, Qa, cldf, fsw)

!---!-------------------------------------------------------------------
!---! AOMIP shortwave forcing
!---!-------------------------------------------------------------------

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block, & ! block dimensions
         ilo,ihi,jlo,jhi       ! beginning and end of physical domain

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         TLON, TLAT     , & ! longitude, latitude
         Qa             , & ! specific humidity
         cldf           , & ! cloud fraction
         hm                 ! land mask

      real (kind=dbl_kind), dimension(nx_block,ny_block),  &
         intent(inout) :: &
         fsw                ! shortwave

      real (kind=dbl_kind) :: &
         hour_angle, &
         solar_time, &
         declin    , &
         cosZ      , &
         e, d      , &
         sw0       , &
         deg2rad   

      integer (kind=int_kind) :: &
         i, j

      do j=jlo,jhi
       do i=ilo,ihi
        deg2rad = pi/c180
        solar_time = mod(real(sec,kind=dbl_kind),secday)/c3600 &
                   + c12*sin(p5*TLON(i,j))
        hour_angle = (c12 - solar_time)*pi/c12
        declin = 23.44_dbl_kind*cos((172._dbl_kind-yday) &
                 * c2*pi/c365)*deg2rad     ! use dayyr instead of c365???
        cosZ = sin(TLAT(i,j))*sin(declin) &
             + cos(TLAT(i,j))*cos(declin)*cos(hour_angle)
        cosZ = max(cosZ,c0)
        e = 1.e5*Qa(i,j)/(0.622_dbl_kind + 0.378_dbl_kind*Qa(i,j))
        d = (cosZ+2.7_dbl_kind)*e*1.e-5_dbl_kind+1.085_dbl_kind*cosZ+p1
        sw0 = 1353._dbl_kind*cosZ**2/d
        sw0 = max(sw0,c0)

        ! total downward shortwave for cice
        Fsw(i,j) = sw0*(c1-p6*cldf(i,j)**3) 
        Fsw(i,j) = Fsw(i,j)*hm(i,j)
       enddo
      enddo

      end subroutine compute_shortwave

!=======================================================================

      subroutine Qa_fixLY(nx_block, ny_block, Tair, Qa)

      use ice_work, only: worka

      integer (kind=int_kind), intent(in) :: &
         nx_block, ny_block ! block dimensions

      real (kind=dbl_kind), dimension(nx_block,ny_block), intent(in) :: &
         Tair               ! air temperature

      real (kind=dbl_kind), dimension(nx_block,ny_block), &
         intent(inout) :: &
         Qa                 ! specific humidity

      worka = Tair - Tffresh
      worka = c2 + (0.7859_dbl_kind + 0.03477_dbl_kind*worka) &
                     /(c1 + 0.00412_dbl_kind*worka) & ! 2+ converts ea mb -> Pa
                + 0.00422_dbl_kind*worka              ! for ice
      ! vapor pressure
      worka = (c10**worka)      ! saturated 
      worka = max(worka,puny)   ! puny over land to prevent division by zero
      ! specific humidity
      worka = 0.622_dbl_kind*worka/(1.e5_dbl_kind-0.378_dbl_kind*worka)

      Qa = min(Qa, worka)

      end subroutine Qa_fixLY

!=======================================================================
! HadGEM or HadGAM atmospheric forcing
!=======================================================================
!
!BOP
!
! !IROUTINE: hadgem_files - construct filenames for HadGEM or HadGAM files
!
! !INTERFACE:
!
      subroutine hadgem_files (yr)
!
! !DESCRIPTION:
!
! Construct filenames based on selected model options
!
! Note: The year number in these filenames does not matter, because
!       subroutine file\_year will insert the correct year.
!
! !REVISION HISTORY:
!
! author: Alison McLaren, Met Office
!
! !USES:
!
      use ice_therm_vertical, only: calc_Tsfc 
      use ice_ocean, only: oceanmixed_ice
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           yr                   ! current forcing year
!
!EOP
!
      integer (kind=int_kind) :: &
           n           ! thickness category index

      ! -----------------------------------------------------------
      ! Rainfall and snowfall
      ! -----------------------------------------------------------

      snow_file = &
           trim(atm_data_dir)//'MONTHLY/snowfall.1996.nc'
           call file_year(snow_file,yr)

      rain_file = &
           trim(atm_data_dir)//'MONTHLY/rainfall.1996.nc'
           call file_year(rain_file,yr)

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(snow_file)
      endif

      if (calc_strair) then

         ! --------------------------------------------------------
         ! Wind velocity
         ! --------------------------------------------------------

         uwind_file = &
           trim(atm_data_dir)//'MONTHLY/u_10.1996.nc'
           call file_year(uwind_file,yr)

         vwind_file = &
           trim(atm_data_dir)//'MONTHLY/v_10.1996.nc'
           call file_year(vwind_file,yr)

         if (my_task == master_task) then
            write (nu_diag,*) trim(uwind_file)
            write (nu_diag,*) trim(vwind_file)
         endif

      else

         ! --------------------------------------------------------
         ! Wind stress
         ! --------------------------------------------------------

         strax_file = &
              trim(atm_data_dir)//'MONTHLY/taux.1996.nc'
         call file_year(strax_file,yr)

         stray_file = &
              trim(atm_data_dir)//'MONTHLY/tauy.1996.nc'
         call file_year(stray_file,yr)

         if (my_task == master_task) then
            write (nu_diag,*) trim(strax_file)
            write (nu_diag,*) trim(stray_file)
         endif

         if (calc_Tsfc .or. oceanmixed_ice) then

            ! --------------------------------------------------
            ! Wind speed
            ! --------------------------------------------------

            wind_file = &
               trim(atm_data_dir)//'MONTHLY/wind_10.1996.nc'
            call file_year(wind_file,yr)

            if (my_task == master_task) then
               write (nu_diag,*) trim(wind_file)
            endif

         endif   ! calc_Tsfc or oceanmixed_ice

      endif  ! calc_strair

      ! --------------------------------------------------------------
      ! Atmosphere properties.  Even if these fields are not 
      ! being used to force the ice (i.e. calc_Tsfc=.false.), they
      ! are still needed to generate forcing for mixed layer model or
      ! to calculate wind stress
      ! --------------------------------------------------------------

       if (calc_Tsfc .or. oceanmixed_ice .or. calc_strair) then  

         fsw_file = &
           trim(atm_data_dir)//'MONTHLY/SW_incoming.1996.nc'
           call file_year(fsw_file,yr)

         flw_file = &
           trim(atm_data_dir)//'MONTHLY/LW_incoming.1996.nc'
           call file_year(flw_file,yr)

         tair_file = &
           trim(atm_data_dir)//'MONTHLY/t_10.1996.nc'
           call file_year(tair_file,yr)

         humid_file = &
           trim(atm_data_dir)//'MONTHLY/q_10.1996.nc'
           call file_year(humid_file,yr)

         rhoa_file = &
           trim(atm_data_dir)//'MONTHLY/rho_10.1996.nc'
           call file_year(rhoa_file,yr)

         if (my_task == master_task) then
            write (nu_diag,*) trim(fsw_file)
            write (nu_diag,*) trim(flw_file)
            write (nu_diag,*) trim(tair_file)
            write (nu_diag,*) trim(humid_file)
            write (nu_diag,*) trim(rhoa_file)
         endif                     ! master_task

      endif ! calc_Tsfc or oceanmixed_ice  or calc_strair

      if (.not. calc_Tsfc) then

         ! ------------------------------------------------------
         ! Sublimation, topmelt and botmelt
         ! ------------------------------------------------------

         do n = 1, ncat

            ! 'topmelt' = fsurf - fcondtop.
            write(topmelt_file(n), '(a,i1,a)')  &
              trim(atm_data_dir)//'MONTHLY/topmeltn',n,'.1996.nc'
              call file_year(topmelt_file(n),yr)

            ! 'botmelt' = fcondtop. 
            write(botmelt_file(n), '(a,i1,a)')  &
              trim(atm_data_dir)//'MONTHLY/botmeltn',n,'.1996.nc'
              call file_year(botmelt_file(n),yr)

         enddo

         ! 'sublim' = - flat / Lsub. 
         sublim_file = &
           trim(atm_data_dir)//'MONTHLY/sublim.1996.nc'
           call file_year(sublim_file,yr)

         if (my_task == master_task) then
            do n = 1, ncat
               write (nu_diag,*) trim(topmelt_file(n))
               write (nu_diag,*) trim(botmelt_file(n))
            enddo
            write (nu_diag,*) trim(sublim_file)

         endif

      endif  ! .not. calc_Tsfc

      end subroutine hadgem_files

!=======================================================================
!
!BOP
!
! !IROUTINE: hadgem_data - read HadGEM or HadGAM atmospheric data
!
! !INTERFACE:
!
      subroutine hadgem_data
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: Alison McLaren, Met Office
!
! !USES:
!
      use ice_domain, only: nblocks
      use ice_flux
      use ice_state, only: aice,aicen
      use ice_ocean, only: oceanmixed_ice
      use ice_therm_vertical, only: calc_Tsfc
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      integer (kind=int_kind) :: &
          i, j        , & ! horizontal indices
          n           , & ! thickness category index
          iblk        , & ! block index
          ixm,ixx,ixp , & ! record numbers for neighboring months
          recnum      , & ! record number
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          dataloc     , & ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
          midmonth        ! middle day of month

      logical (kind=log_kind) :: readm

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
            topmelt, & ! temporary fields
            botmelt, &
            sublim

      character (char_len) :: & 
            fieldname    ! field name in netcdf file

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

      ! -----------------------------------------------------------
      ! Rainfall and snowfall
      ! -----------------------------------------------------------

      fieldname='rainfall'
      call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, rain_file, fieldname, frain_data, &
                      field_loc_center, field_type_scalar)
      fieldname='snowfall'
      call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, snow_file, fieldname, fsnow_data, &
                      field_loc_center, field_type_scalar)

      ! Interpolate to current time step
      call interpolate_data (fsnow_data, fsnow)
      call interpolate_data (frain_data, frain)

      if (calc_strair) then

         ! --------------------------------------------------------
         ! Wind velocity
         ! --------------------------------------------------------

         fieldname='u_10'
         call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, uwind_file, fieldname, uatm_data, &
                      field_loc_center, field_type_vector)
         fieldname='v_10'
         call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, vwind_file, fieldname, vatm_data, &
                      field_loc_center, field_type_vector)

         ! Interpolate to current time step
         call interpolate_data (uatm_data, uatm)
         call interpolate_data (vatm_data, vatm)

      else

         ! --------------------------------------------------------
         ! Wind stress
         ! --------------------------------------------------------

         fieldname='taux'
         call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, strax_file, fieldname, strax_data, &
                      field_loc_center, field_type_vector)
         fieldname='tauy'
         call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, stray_file, fieldname, stray_data, &
                      field_loc_center, field_type_vector)

         ! Interpolate to current time step
         call interpolate_data (strax_data, strax)
         call interpolate_data (stray_data, stray)

         if (calc_Tsfc .or. oceanmixed_ice) then

            ! --------------------------------------------------
            ! Wind speed
            ! --------------------------------------------------

            fieldname='wind_10'
            call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, wind_file, fieldname, wind_data, &
                      field_loc_center, field_type_scalar)

            ! Interpolate to current time step
            call interpolate_data (wind_data, wind)

         endif   ! calc_Tsfc or oceanmixed_ice

      endif      ! calc_strair

      ! -----------------------------------------------------------
      ! SW incoming, LW incoming, air temperature, density and 
      ! humidity at 10m.  
      !
      ! Even if these fields are not being used to force the ice 
      ! (i.e. calc_Tsfc=.false.), they are still needed to generate 
      ! forcing for mixed layer model or to calculate wind stress
      ! -----------------------------------------------------------

      if (calc_Tsfc .or. oceanmixed_ice .or. calc_strair) then  

         fieldname='SW_incoming'
         call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, fsw_file, fieldname, fsw_data, &
                      field_loc_center, field_type_scalar)
         fieldname='LW_incoming'
         call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, flw_file, fieldname, flw_data, &
                      field_loc_center, field_type_scalar)
         fieldname='t_10'
         call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, tair_file, fieldname, Tair_data, &
                      field_loc_center, field_type_scalar)
         fieldname='rho_10'
         call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, rhoa_file, fieldname, rhoa_data, &
                      field_loc_center, field_type_scalar)
         fieldname='q_10'
         call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, humid_file, fieldname, Qa_data, &
                      field_loc_center, field_type_scalar)

         ! Interpolate onto current timestep

         call interpolate_data (fsw_data,   fsw)
         call interpolate_data (flw_data,  flw)
         call interpolate_data (Tair_data, Tair)
         call interpolate_data (rhoa_data, rhoa)
         call interpolate_data (Qa_data,   Qa)

      endif       ! calc_Tsfc or oceanmixed_ice or calc_strair

      if (.not. calc_Tsfc) then

         ! ------------------------------------------------------
         ! Sublimation, topmelt and botmelt
         ! ------------------------------------------------------

         fieldname='sublim'
         call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, sublim_file, fieldname, sublim_data, &
                      field_loc_center, field_type_scalar)

         ! Interpolate to current time step
         call interpolate_data (sublim_data, sublim)

         do n = 1, ncat
            write(fieldname, '(a,i1)') 'topmeltn',n
            call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
              maxrec, topmelt_file(n), fieldname, topmelt_data(:,:,:,:,n), &
              field_loc_center, field_type_scalar)

            write(fieldname, '(a,i1)') 'botmeltn',n
            call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
              maxrec, botmelt_file(n), fieldname, botmelt_data(:,:,:,:,n), &
              field_loc_center, field_type_scalar)

            call interpolate_data (topmelt_data(:,:,:,:,n), topmelt)
            call interpolate_data (botmelt_data(:,:,:,:,n), botmelt)

            !--------------------------------------------------------
            ! Convert from UM variables to CICE variables
            !  topmelt = fsurf - fcondtop
            !  botmelt = fcondtop  (as zero layer)
            !
            ! Convert UM sublimation data into CICE LH flux
            ! (sublim = - flatn / Lsub) and have same value for all 
            ! categories
            !--------------------------------------------------------

            do iblk = 1, nblocks
               do j = 1, ny_block
               do i = 1, nx_block

                  fcondtopn_f(i,j,n,iblk) = botmelt(i,j,iblk)
                  fsurfn_f(i,j,n,iblk)    = topmelt(i,j,iblk) & 
                                            + botmelt(i,j,iblk)
                  flatn_f(i,j,n,iblk)    = - sublim(i,j,iblk)*Lsub
                  
               enddo
               enddo

            enddo

         enddo  ! ncat

      endif   ! .not. calc_Tsfc 

      end subroutine hadgem_data

!=======================================================================
! monthly forcing 
!=======================================================================
!
!BOP
!
! !IROUTINE: monthly_files - construct filenames for monthly data
!
! !INTERFACE:
!
      subroutine monthly_files (yr)
!
! !DESCRIPTION:
!
! Construct filenames based on the LANL naming conventions for NCAR data.
! Edit for other directory structures or filenames.
! Note: The year number in these filenames does not matter, because
!       subroutine file_year will insert the correct year.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer (kind=int_kind), intent(in) :: &
           yr                   ! current forcing year
!
!EOP
!
      flw_file = &
           trim(atm_data_dir)//'MONTHLY/cldf.omip.dat'

      rain_file = &
           trim(atm_data_dir)//'MONTHLY/prec.nmyr.dat'

      tair_file = &
           trim(atm_data_dir)//'MONTHLY/t_10.1996.dat'
      call file_year(tair_file,yr)

      humid_file = &
           trim(atm_data_dir)//'MONTHLY/q_10.1996.dat'
      call file_year(humid_file,yr)

      ! stress/speed is used instead of wind components
      strax_file = &
           trim(atm_data_dir)//'MONTHLY/strx.1996.dat'
      call file_year(strax_file,yr)

      stray_file = &
           trim(atm_data_dir)//'MONTHLY/stry.1996.dat'
      call file_year(stray_file,yr)

      wind_file = &
           trim(atm_data_dir)//'MONTHLY/wind.1996.dat'
      call file_year(wind_file,yr)

      if (my_task == master_task) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'Forcing data year = ', fyear         
         write (nu_diag,*) 'Atmospheric data files:'
         write (nu_diag,*) trim(flw_file)
         write (nu_diag,*) trim(rain_file)
         write (nu_diag,*) trim(tair_file)
         write (nu_diag,*) trim(humid_file)
         write (nu_diag,*) trim(uwind_file)
         write (nu_diag,*) trim(vwind_file)
      endif                     ! master_task

      end subroutine monthly_files

!=======================================================================
!
!BOP
!
! !IROUTINE: monthly_data - read monthly atmospheric data
!
! !INTERFACE:
!
      subroutine monthly_data
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      use ice_global_reductions
      use ice_domain, only: nblocks, distrb_info, blocks_ice
      use ice_flux 
      use ice_grid, only: hm, tlon, tlat, tmask, umask
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: & 
          i, j        , &
          imx,ipx     , & ! record numbers for neighboring months
          recnum      , & ! record number
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          midmonth    , & ! middle day of month
          iblk        , & ! block index
          ilo,ihi,jlo,jhi ! beginning and end of physical domain

      real (kind=dbl_kind) :: &
          vmin, vmax

      logical (kind=log_kind) :: readm

      type (block) :: &
         this_block           ! block information for current block
      
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
      imx  = mod(month+maxrec-2,maxrec) + 1
      ipx  = mod(month,         maxrec) + 1
      if (mday >= midmonth) imx = 99  ! other two points will be used
      if (mday <  midmonth) ipx = 99

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

      call read_clim_data (readm, 0, imx, month, ipx,  &
             flw_file, cldf_data, &
             field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, imx, month, ipx,  &
             rain_file, fsnow_data, &
             field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, imx, month, ipx,  &
             tair_file, Tair_data, &
             field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, imx, month, ipx,  &
             humid_file, Qa_data, &
             field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, imx, month, ipx,  &
             wind_file, wind_data, &
             field_loc_center, field_type_scalar)
      call read_clim_data (readm, 0, imx, month, ipx,  &
             strax_file, strax_data, &
             field_loc_center, field_type_vector)
      call read_clim_data (readm, 0, imx, month, ipx,  &
             stray_file, stray_data, &
             field_loc_center, field_type_vector)

      call interpolate_data (cldf_data, cldf)
      call interpolate_data (fsnow_data, fsnow)  ! units mm/s = kg/m^2/s
      call interpolate_data (Tair_data, Tair)
      call interpolate_data (Qa_data, Qa)
      call interpolate_data (wind_data, wind)
      call interpolate_data (strax_data, strax)
      call interpolate_data (stray_data, stray)

      do iblk = 1, nblocks
        call Qa_fixLY(nx_block,  ny_block, &
                                 Tair (:,:,iblk), &
                                 Qa   (:,:,iblk))

        do j = 1, ny_block
          do i = 1, nx_block
            Qa   (i,j,iblk) = Qa   (i,j,iblk) * hm(i,j,iblk)
            Tair (i,j,iblk) = Tair (i,j,iblk) * hm(i,j,iblk)
            wind (i,j,iblk) = wind (i,j,iblk) * hm(i,j,iblk)
            strax(i,j,iblk) = strax(i,j,iblk) * hm(i,j,iblk)
            stray(i,j,iblk) = stray(i,j,iblk) * hm(i,j,iblk)
          enddo
        enddo

      ! AOMIP
      this_block = get_block(blocks_ice(iblk),iblk)         
      ilo = this_block%ilo
      ihi = this_block%ihi
      jlo = this_block%jlo
      jhi = this_block%jhi

      call compute_shortwave(nx_block, ny_block, &
                             ilo, ihi, jlo, jhi, &
                             TLON (:,:,iblk), &
                             TLAT (:,:,iblk), &
                             hm   (:,:,iblk), &
                             Qa   (:,:,iblk), &
                             cldf (:,:,iblk), &
                             fsw  (:,:,iblk))

      enddo  ! iblk

      ! Save record number
      oldrecnum = recnum

         if (dbug) then
           if (my_task == master_task) write (nu_diag,*) 'LY_bulk_data'
           vmin = global_minval(fsw,distrb_info,tmask)
           vmax = global_maxval(fsw,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'fsw',vmin,vmax 
           vmin = global_minval(cldf,distrb_info,tmask)
           vmax = global_maxval(cldf,distrb_info,tmask)
           if (my_task.eq.master_task) & 
               write (nu_diag,*) 'cldf',vmin,vmax
           vmin =global_minval(fsnow,distrb_info,tmask)
           vmax =global_maxval(fsnow,distrb_info,tmask)
           if (my_task.eq.master_task) & 
               write (nu_diag,*) 'fsnow',vmin,vmax
           vmin = global_minval(Tair,distrb_info,tmask)
           vmax = global_maxval(Tair,distrb_info,tmask)
           if (my_task.eq.master_task) & 
               write (nu_diag,*) 'Tair',vmin,vmax
           vmin = global_minval(wind,distrb_info,umask)
           vmax = global_maxval(wind,distrb_info,umask)
           if (my_task.eq.master_task) & 
               write (nu_diag,*) 'wind',vmin,vmax
           vmin = global_minval(strax,distrb_info,umask)
           vmax = global_maxval(strax,distrb_info,umask)
           if (my_task.eq.master_task) & 
               write (nu_diag,*) 'strax',vmin,vmax
           vmin = global_minval(stray,distrb_info,umask)
           vmax = global_maxval(stray,distrb_info,umask)
           if (my_task.eq.master_task) & 
               write (nu_diag,*) 'stray',vmin,vmax
           vmin = global_minval(Qa,distrb_info,tmask)
           vmax = global_maxval(Qa,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'Qa',vmin,vmax

        endif                   ! dbug

      end subroutine monthly_data

!=======================================================================
! Climatological ocean forcing
!=======================================================================
!BOP
!
! !IROUTINE: ocn_data_clim - interpolate sss, sst; restore sst
!
! !INTERFACE:
!
      subroutine ocn_data_clim (dt)
!
! !DESCRIPTION:
!
! Interpolate monthly sss, sst data to timestep.
! Restore prognostic sst to data.
! Interpolate fields from U grid to T grid if necessary.
!
! !REVISION HISTORY:
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL
!
! !USES:
!
      use ice_domain, only: nblocks
      use ice_ocean
      use ice_flux, only: Tf, sss, sst, uocn, vocn, ss_tltx, ss_tlty, Tfrzpt
      use ice_grid, only: t2ugrid_vector
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: &
          i, j, iblk  , & ! horizontal indices
          ixm,ixp     , & ! record numbers for neighboring months
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          midmonth        ! middle day of month

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
          sstdat              ! data value toward which SST is restored

      logical (kind=log_kind) :: readm

      if (my_task == master_task .and. istep == 1) then
         if (trim(sss_data_type)=='clim') then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'SSS data interpolated to timestep:'
            write (nu_diag,*) trim(sss_file)
         endif
         if (trim(sst_data_type)=='clim') then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'SST data interpolated to timestep:'
            write (nu_diag,*) trim(sst_file)
            if (restore_sst) write (nu_diag,*) &
              'SST restoring timescale (days) =', trestore
         endif
      endif                     ! my_task, istep

    !-------------------------------------------------------------------
    ! monthly data
    !
    ! Assume that monthly data values are located in the middle of the
    ! month.
    !-------------------------------------------------------------------

      if (trim(sss_data_type)=='clim' .or.  &
          trim(sst_data_type)=='clim') then

         midmonth = 15          ! data is given on 15th of every month
!!!      midmonth = fix(p5 * real(daymo(month)))  ! exact middle

         ! Compute record numbers for surrounding months
         maxrec = 12
         ixm  = mod(month+maxrec-2,maxrec) + 1
         ixp  = mod(month,         maxrec) + 1
         if (mday >= midmonth) ixm = 99 ! other two points will be used
         if (mday <  midmonth) ixp = 99

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

      endif   ! sss/sst_data_type

    !-------------------------------------------------------------------
    ! Read two monthly SSS values and interpolate.
    ! Note: SSS is restored instantaneously to data.
    !-------------------------------------------------------------------

      if (trim(sss_data_type)=='clim') then
         call read_clim_data (readm, 0, ixm, month, ixp, &
                              sss_file, sss_data, &
                              field_loc_center, field_type_scalar)
         call interpolate_data (sss_data, sss)

         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sss(i,j,iblk) = max(sss(i,j,iblk), c0)
               if (trim(Tfrzpt) == 'constant') then
                  Tf (i,j,iblk) = -1.8_dbl_kind ! deg C
               else ! default:  Tfrzpt = 'linear_S'
                  Tf (i,j,iblk) = -depressT * sss(i,j,iblk) ! deg C
               endif
            enddo
            enddo
         enddo
      endif

    !-------------------------------------------------------------------
    ! Read two monthly SST values and interpolate.
    ! Restore toward interpolated value.
    !-------------------------------------------------------------------

      if (trim(sst_data_type)=='clim') then
         call read_clim_data (readm, 0, ixm, month, ixp, &
                              sst_file, sst_data, &
                              field_loc_center, field_type_scalar)
         call interpolate_data (sst_data, sstdat)

         if (restore_sst) then
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sst(i,j,iblk) = sst(i,j,iblk)  &
                         + (sstdat(i,j,iblk)-sst(i,j,iblk))*dt/trest
            enddo
            enddo
         enddo
         endif
      endif

      end subroutine ocn_data_clim

!=======================================================================
! NCAR CCSM M-configuration (AIO) ocean forcing
!=======================================================================
!
!BOP
!
! !IROUTINE: ocn_data_ncar_init - reads data set
!
! !INTERFACE:
!
      subroutine ocn_data_ncar_init
!
! !DESCRIPTION:
!
! Reads NCAR pop ocean forcing data set 'pop_frc_gx1v3_010815.nc'
! 
! List of ocean forcing fields: Note that order is important!
! (order is determined by field list in vname).
! 
! For ocean mixed layer-----------------------------units 
! 
! 1  sst------temperature---------------------------(C)   \\
! 2  sss------salinity------------------------------(ppt) \\
! 3  hbl------depth---------------------------------(m)   \\ 
! 4  u--------surface u current---------------------(m/s) \\
! 5  v--------surface v current---------------------(m/s) \\
! 6  dhdx-----surface tilt x direction--------------(m/m) \\
! 7  dhdy-----surface tilt y direction--------------(m/m) \\
! 8  qdp------ocean sub-mixed layer heat flux-------(W/m2)\\ 
!
! Fields 4, 5, 6, 7 are on the U-grid; 1, 2, 3, and 8 are
! on the T-grid.
!
! !REVISION HISTORY:
!
! authors: Bruce Briegleb, NCAR
!          Elizabeth Hunke, LANL
!
! !USES:
!
      use ice_domain, only: nblocks, distrb_info
      use ice_gather_scatter
      use ice_exit
      use ice_work, only: work1
      use ice_read_write
#ifdef ncdf
      use netcdf
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!
      integer (kind=int_kind) :: & 
        n   , & ! field index
        m   , & ! month index
        nrec, & ! record number for direct access
        nbits

      character(char_len) :: &
        vname(nfld) ! variable names to search for in file
      data vname /  &
           'T',      'S',      'hblt',  'U',     'V', &
           'dhdx',   'dhdy',   'qdp' /

      integer (kind=int_kind) :: &
        fid        , & ! file id 
        dimid          ! dimension id 

      integer (kind=int_kind) :: &
        status  , & ! status flag
        nlat    , & ! number of longitudes of data
        nlon        ! number of latitudes  of data

      if (my_task == master_task) then

         write (nu_diag,*) 'WARNING: evp_prep calculates surface tilt'
         write (nu_diag,*) 'WARNING: stress from geostrophic currents,'
         write (nu_diag,*) 'WARNING: not data from ocean forcing file.'
         write (nu_diag,*) 'WARNING: Alter ice_dyn_evp.F if desired.'

         if (restore_sst) write (nu_diag,*)  &
             'SST restoring timescale = ',trestore,' days' 

         sst_file = trim(ocn_data_dir)//oceanmixed_file ! not just sst

        !---------------------------------------------------------------
        ! Read in ocean forcing data from an existing file
        !---------------------------------------------------------------
        write (nu_diag,*) 'ocean mixed layer forcing data file = ', &
                           sst_file

      endif ! master_task

      if (trim(ocn_data_format) == 'nc') then
#ifdef ncdf
        if (my_task == master_task) then
          call ice_open_nc(sst_file, fid)

!          status = nf90_inq_dimid(fid,'nlon',dimid)
          status = nf90_inq_dimid(fid,'ni',dimid)
          status = nf90_inquire_dimension(fid,dimid,len=nlon)
  
!          status = nf90_inq_dimid(fid,'nlat',dimid)
          status = nf90_inq_dimid(fid,'nj',dimid)
          status = nf90_inquire_dimension(fid,dimid,len=nlat)

          if( nlon .ne. nx_global ) then
            call abort_ice ('ice: ocn frc file nlon ne nx_global')
          endif
          if( nlat .ne. ny_global ) then
            call abort_ice ('ice: ocn frc file nlat ne ny_global')
          endif

        endif ! master_task

        ! Read in ocean forcing data for all 12 months
        do n=1,nfld
          do m=1,12
                
            ! Note: netCDF does single to double conversion if necessary
            if (n >= 4 .and. n <= 7) then
               call ice_read_nc(fid, m, vname(n), work1, dbug, &
                                field_loc_NEcorner, field_type_vector)
            else
               call ice_read_nc(fid, m, vname(n), work1, dbug, &
                                field_loc_center, field_type_scalar)
            endif
            ocn_frc_m(:,:,:,n,m) = work1(:,:,:)

          enddo               ! month loop
        enddo               ! field loop

        if (my_task == master_task) status = nf90_close(fid)
#endif

      else  ! binary format

        nbits = 64
        call ice_open (nu_forcing, sst_file, nbits)

        nrec = 0
        do n=1,nfld
           do m=1,12
              nrec = nrec + 1
              if (n >= 4 .and. n <= 7) then
                call ice_read (nu_forcing, nrec, work1, 'rda8', dbug, &
                               field_loc_NEcorner, field_type_vector)
              else
                call ice_read (nu_forcing, nrec, work1, 'rda8', dbug, &
                               field_loc_center, field_type_scalar)
              endif
              ocn_frc_m(:,:,:,n,m) = work1(:,:,:)
           enddo               ! month loop
        enddo               ! field loop
        close (nu_forcing)

      endif

!echmod - currents cause Fram outflow to be too large
              ocn_frc_m(:,:,:,4,:) = c0
              ocn_frc_m(:,:,:,5,:) = c0
!echmod

      end subroutine ocn_data_ncar_init

!=======================================================================
!
!BOP
!
! !IROUTINE: ocn_data_ncar - interpolates data to timestep
!
! !INTERFACE:
!
      subroutine ocn_data_ncar(dt)
!
! !DESCRIPTION:
!
! Interpolate monthly ocean data to timestep.
! Restore sst if desired. sst is updated with surface fluxes in ice_ocean.F.
! 
! !REVISION HISTORY:
!
! authors: same as module
!
! !USES:
!
      use ice_global_reductions
      use ice_domain, only: nblocks, distrb_info
      use ice_gather_scatter
      use ice_exit
      use ice_work, only: work_g1, work1
      use ice_flux, only: sss, sst, Tf, uocn, vocn, ss_tltx, ss_tlty, &
            qdp, hmix, Tfrzpt
!      use ice_ocean
      use ice_restart, only: restart
      use ice_grid, only: hm, tmask, umask
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
!
!EOP
!
      integer (kind=int_kind) :: & 
          i, j, n, iblk   , &
          imx,ipx         , & ! record numbers for neighboring months
          maxrec          , & ! maximum record number
          recslot         , & ! spline slot for current record
          midmonth            ! middle day of month

      real (kind=dbl_kind) :: &
          vmin, vmax

    !-------------------------------------------------------------------
    ! monthly data 
    !
    ! Assume that monthly data values are located in the middle of the 
    ! month.
    !-------------------------------------------------------------------
      
      midmonth = 15  ! data is given on 15th of every month
!      midmonth = fix(p5 * real(daymo(month),kind=dbl_kind))  ! exact middle

      ! Compute record numbers for surrounding months
      maxrec = 12
      imx  = mod(month+maxrec-2,maxrec) + 1
      ipx  = mod(month,         maxrec) + 1
      if (mday >= midmonth) imx = 99  ! other two points will be used
      if (mday <  midmonth) ipx = 99

      ! Determine whether interpolation will use values 1:2 or 2:3
      ! recslot = 2 means we use values 1:2, with the current value (2)
      !  in the second slot
      ! recslot = 1 means we use values 2:3, with the current value (2)
      !  in the first slot
      recslot = 1                             ! latter half of month
      if (mday < midmonth) recslot = 2        ! first half of month

      ! Find interpolation coefficients
      call interp_coeff_monthly (recslot)

      do n = nfld, 1, -1
        do iblk = 1, nblocks
        ! use sst_data arrays as temporary work space until n=1
        if (imx /= 99) then  ! first half of month
          sst_data(:,:,1,iblk) = ocn_frc_m(:,:,iblk,n,imx)
          sst_data(:,:,2,iblk) = ocn_frc_m(:,:,iblk,n,month)
        else                 ! second half of month
          sst_data(:,:,1,iblk) = ocn_frc_m(:,:,iblk,n,month)
          sst_data(:,:,2,iblk) = ocn_frc_m(:,:,iblk,n,ipx)
        endif
        enddo
        call interpolate_data (sst_data,work1)
        ! masking by hm is necessary due to NaNs in the data file
        do j = 1, ny_block 
          do i = 1, nx_block 
            if (n == 2) sss    (i,j,:) = c0
            if (n == 3) hmix   (i,j,:) = c0
            if (n == 4) uocn   (i,j,:) = c0
            if (n == 5) vocn   (i,j,:) = c0
            if (n == 6) ss_tltx(i,j,:) = c0
            if (n == 7) ss_tlty(i,j,:) = c0
            if (n == 8) qdp    (i,j,:) = c0
            do iblk = 1, nblocks
              if (hm(i,j,iblk) == c1) then
                if (n == 2) sss    (i,j,iblk) = work1(i,j,iblk)
                if (n == 3) hmix   (i,j,iblk) = work1(i,j,iblk)
                if (n == 4) uocn   (i,j,iblk) = work1(i,j,iblk)
                if (n == 5) vocn   (i,j,iblk) = work1(i,j,iblk)
                if (n == 6) ss_tltx(i,j,iblk) = work1(i,j,iblk)
                if (n == 7) ss_tlty(i,j,iblk) = work1(i,j,iblk)
                if (n == 8) qdp    (i,j,iblk) = work1(i,j,iblk)
              endif
            enddo
          enddo
        enddo
      enddo

      do j = 1, ny_block 
        do i = 1, nx_block 
          sss (i,j,:) = max (sss(i,j,:), c0) 
            if (trim(Tfrzpt) == 'constant') then
               Tf (i,j,:) = -1.8_dbl_kind ! deg C
            else ! default:  Tfrzpt = 'linear_S'
               Tf (i,j,:) = -depressT * sss(i,j,:) ! deg C
            endif
          hmix(i,j,:) = max(hmix(i,j,:), c0) 
        enddo 
      enddo 

      if (restore_sst) then
        do j = 1, ny_block 
         do i = 1, nx_block 
           sst(i,j,:) = sst(i,j,:) + (work1(i,j,:)-sst(i,j,:))*dt/trest 
         enddo 
        enddo 
!     else sst is only updated in ice_ocean.F
      endif

      ! initialize sst properly on first step
      if (istep1 <= 1 .and. .not. (restart)) then
        call interpolate_data (sst_data,sst)
        do iblk = 1, nblocks
         do j = 1, ny_block 
          do i = 1, nx_block 
            if (hm(i,j,iblk) == c1) then
              sst(i,j,iblk) =  max (sst(i,j,iblk), Tf(i,j,iblk)) 
            else
              sst(i,j,iblk) = c0
            endif
          enddo 
         enddo 
        enddo 
      endif

      if (dbug) then
         if (my_task == master_task)  &
               write (nu_diag,*) 'ocn_data_ncar'
           vmin = global_minval(Tf,distrb_info,tmask)
           vmax = global_maxval(Tf,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'Tf',vmin,vmax
           vmin = global_minval(sst,distrb_info,tmask)
           vmax = global_maxval(sst,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'sst',vmin,vmax
           vmin = global_minval(sss,distrb_info,tmask)
           vmax = global_maxval(sss,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'sss',vmin,vmax
           vmin = global_minval(hmix,distrb_info,tmask)
           vmax = global_maxval(hmix,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'hmix',vmin,vmax
           vmin = global_minval(uocn,distrb_info,umask)
           vmax = global_maxval(uocn,distrb_info,umask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'uocn',vmin,vmax
           vmin = global_minval(vocn,distrb_info,umask)
           vmax = global_maxval(vocn,distrb_info,umask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'vocn',vmin,vmax
           vmin = global_minval(ss_tltx,distrb_info,umask)
           vmax = global_maxval(ss_tltx,distrb_info,umask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'ss_tltx',vmin,vmax
           vmin = global_minval(ss_tlty,distrb_info,umask)
           vmax = global_maxval(ss_tlty,distrb_info,umask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'ss_tlty',vmin,vmax
           vmin = global_minval(qdp,distrb_info,tmask)
           vmax = global_maxval(qdp,distrb_info,tmask)
           if (my_task.eq.master_task)  &
               write (nu_diag,*) 'qdp',vmin,vmax
      endif

      end subroutine ocn_data_ncar

!=======================================================================
!
!BOP
!
! !IROUTINE: ocn_data_hadgem - read HadGEM ocean data
!
! !INTERFACE:
!
      subroutine ocn_data_hadgem(dt)
!
! !DESCRIPTION:
!  Reads in HadGEM ocean forcing data as required from netCDF files
!  Current options (selected by sst_data_type)
!  hadgem_sst: 		read in sst only 
!  hadgem_sst_uvocn:	read in sst plus uocn and vocn	
!        
!
! !REVISION HISTORY:
!
! authors: Ann Keen, Met Office
!
! !USES:
!
      use ice_domain, only: nblocks
      use ice_flux, only: sst, uocn, vocn
      use ice_grid, only: t2ugrid_vector, ANGLET
!
! !INPUT/OUTPUT PARAMETERS:
!
!
!EOP
!
      real (kind=dbl_kind), intent(in) :: &
         dt      ! time step
 
     integer (kind=int_kind) :: &
          i, j        , & ! horizontal indices
          n           , & ! thickness category index
          iblk        , & ! block index
          ixm,ixx,ixp , & ! record numbers for neighboring months
          recnum      , & ! record number
          maxrec      , & ! maximum record number
          recslot     , & ! spline slot for current record
          dataloc     , & ! = 1 for data located in middle of time interval
                          ! = 2 for date located at end of time interval
          midmonth        ! middle day of month

      real (kind=dbl_kind), dimension(nx_block,ny_block,max_blocks) :: &
          sstdat              ! data value toward which SST is restored

      real (kind=dbl_kind) :: workx, worky

      logical (kind=log_kind) :: readm

      character (char_len) :: & 
            fieldname    	! field name in netcdf file

      character (char_len_long) :: & 
            filename    	! name of netCDF file

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


      if (my_task == master_task .and. istep == 1) then
         write (nu_diag,*) ' '
         write (nu_diag,*) 'SST data interpolated to timestep:'
         write (nu_diag,*) trim(ocn_data_dir)//'MONTHLY/sst.1997.nc'
         if (restore_sst) write (nu_diag,*) &
              'SST restoring timescale (days) =', trestore
         if (trim(sst_data_type)=='hadgem_sst_uvocn') then
            write (nu_diag,*) ' '
            write (nu_diag,*) 'uocn and vocn interpolated to timestep:'
            write (nu_diag,*) trim(ocn_data_dir)//'MONTHLY/uocn.1997.nc'
            write (nu_diag,*) trim(ocn_data_dir)//'MONTHLY/vocn.1997.nc'
         endif
      endif                     ! my_task, istep


      ! -----------------------------------------------------------
      ! SST
      ! -----------------------------------------------------------
      sst_file = trim(ocn_data_dir)//'MONTHLY/sst.1997.nc'	
      fieldname='sst'
      call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, sst_file, fieldname, sst_data, &
                      field_loc_center, field_type_scalar)
      
      ! Interpolate to current time step
      call interpolate_data (sst_data, sstdat)

      ! Restore SSTs if required
        if (restore_sst) then
         do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block
               sst(i,j,iblk) = sst(i,j,iblk)  &
                         + (sstdat(i,j,iblk)-sst(i,j,iblk))*dt/trest
            enddo
            enddo
         enddo
         endif      


      ! -----------------------------------------------------------
      ! Ocean currents
      ! --------------
      ! Values read in are on T grid and oriented geographically, hence 
      ! vectors need to be rotated to model grid and then interpolated
      ! to U grid.   
      ! Also need to be converted from cm s-1 (UM) to m s-1 (CICE)
      ! -----------------------------------------------------------

      if (trim(sst_data_type)=='hadgem_sst_uvocn') then

      	filename = trim(ocn_data_dir)//'MONTHLY/uocn.1997.nc'	
      	fieldname='uocn'
      	call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, filename, fieldname, uocn_data, &
                      field_loc_center, field_type_vector)
      
      	! Interpolate to current time step
      	call interpolate_data (uocn_data, uocn)

      	filename = trim(ocn_data_dir)//'MONTHLY/vocn.1997.nc'	
      	fieldname='vocn'
      	call read_data_nc (readm, 0, fyear, ixm, month, ixp, &
                      maxrec, filename, fieldname, vocn_data, &
                      field_loc_center, field_type_vector)
      
      	! Interpolate to current time step
      	call interpolate_data (vocn_data, vocn)

     !----------------------------------------------------------------- 
     ! Rotate zonal/meridional vectors to local coordinates, 
     ! and change  units
     !----------------------------------------------------------------- 

        do iblk = 1, nblocks
            do j = 1, ny_block
            do i = 1, nx_block

               workx      = uocn(i,j,iblk) 
               worky      = vocn(i,j,iblk)
               uocn(i,j,iblk) = workx*cos(ANGLET(i,j,iblk)) & 
                                  + worky*sin(ANGLET(i,j,iblk))   
               vocn(i,j,iblk) = worky*cos(ANGLET(i,j,iblk)) & 
                                  - workx*sin(ANGLET(i,j,iblk))

		uocn(i,j,iblk) = uocn(i,j,iblk) * cm_to_m
		vocn(i,j,iblk) = vocn(i,j,iblk) * cm_to_m

            enddo		! i
            enddo		! j
        enddo		! nblocks

     !----------------------------------------------------------------- 
     ! Interpolate to U grid 
     !----------------------------------------------------------------- 

	call t2ugrid_vector(uocn)
	call t2ugrid_vector(vocn)


     endif    !   sst_data_type = hadgem_sst_uvocn


     end subroutine ocn_data_hadgem


!=======================================================================

      end module ice_forcing

!=======================================================================
