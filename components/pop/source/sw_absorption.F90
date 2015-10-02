!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module sw_absorption

!BOP
! !MODULE: sw_absorption
!
! !DESCRIPTION:
!  This module computes shortwave absorption
!
!   sw_absorption_type 
!     = 'top-layer'   :  all sw absorbed in top ocean layer
!     = 'jerlov'   :  sw absorption by jerlov water type
!     = 'chlorophyll' :  sw absorption by chlorophyll amount
!
!   For chlorophyll based absorption, transmission from surface 
!   (z=0) to level z = A_1 exp(-B_1 z) + A_2 exp(-B_2 z) 
!
!   A,B coefficients from Table 1a, J. Carter Ohlmann, 2002, Ocean Radiant 
!   Heating in Climate Models, Journal of Climate, 2003, in press. 
!   Chlorophyll concentration in mg/m^3. Values for chl < .01 mg/m^3 are 
!   linear extrpolations from two smallest chlorophyll amounts. Values 
!   over 3 mg/m^3 up to 10 mg/m^3 were provided by Carter specifically 
!   for this application. Chlorophyll amounts limited from .001 to 10 mg/m^3
!
!   Note that A,B values are relative to the net shortwave flux just below 
!   the ocean surface [E_d(0-) in Carter notation]. This means ocean
!   albedo effects have been accounted for (i.e. E_d(0-) = E_d(0+)(1-alpha_o)
!   where E_d(0+) is the down incident shortwave flux at the ocean surface
!   and alpha_o is the ocean albedo).
!
!   Note also that the second order effects of solar zenith angle and
!   cloud cover have been ignored in this implementation.
!
!   Chlorophyll monthly distribution is provided over all ocean, and used to
!   compute regional distribution of shortwave absorption in the upper ocean.
!
!   Chlorophyll based transmissions are calculated once a month and
!   held fixed over the month in this implementation.

! !REVISION HISTORY:
!  SVN:$Id: sw_absorption.F90 26603 2011-01-28 23:09:02Z njn01 $

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod

   use kinds_mod
   use domain_size
   use domain
   use constants
   use io
   use io_types
   use grid
   use forcing_tools
   use time_management
   use prognostic
   use forcing_shf
   use tavg, only: define_tavg_field, accumulate_tavg_field
   use exit_mod
   use registry, only: registry_match

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_sw_absorption, &
             add_sw_absorb,      &
             sw_absorb_frac,     &
             set_chl,            &
             sw_trans_chl

! !PUBLIC DATA MEMBERS:

   character (char_len), public ::       &
      sw_absorption_type          ! type of sw_absorption

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic),public :: &
      TRANS,                     &! transmission (1 to 0) from surface
                                  ! to depth ZTRANS
      TRANSKM1,                  &! transmission from surface to level k-1
      ZTRANS                      ! depth of transmission      cm

!EOP
!BOC


   character (char_len) ::       &
      chl_option,                &! chlorophyll option ['file', 'model']
      chl_filename,              &! chlorophyll data file name
      chl_file_fmt,              &! chlorophyll data file format
      chl_data_name               ! chlorophyll short name (eg, 'CHL')
 
   integer (int_kind) ::         &
      chl_bndy_loc,              &! location and field types for ghost
      chl_bndy_type               !    cell update routines

   integer (int_kind) ::         &
      jerlov_water_type           ! jerlov water type from 1 to 5

   real (r8), dimension(0:km) :: &
      sw_absorb                   ! sw transmission for jerlov water type

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      CHL                         ! chlorophyll amount      mg/m^3

   real (r8) ::   &
      chlmin,     &! minimum chlorophyll amount    mg/m^3
      chlmax,     &! maximum chlorophyll amount    mg/m^3
      dlogchl      ! delta log chl           

   integer (int_kind), dimension(nx_block,ny_block,max_blocks_clinic) :: &
     CHLINDX       ! chlorophyll index from table

   real (r8), allocatable, dimension(:,:,:,:) :: &
      CHL_DATA    ! global chlorophyll climatology by month

   integer (int_kind), parameter :: &
      nsub = 400,  &! total number chlorophyll concentrations in look up table
      ksol = 2*km, &! (0:ksol) transmission levels (layer interface + mid-point)
      nchl = 31     ! number of chlorophyll absorption coefficients

   real (r8), dimension(0:ksol,0:nsub) :: & 
      Tr       ! table look-up transmissions (level x chlorophyll)

   real (r8), dimension(0:ksol) ::        &
      ztr      ! transmission array levels (layer interface + mid-point)

   real (r8) ::  &
      A1,B1,A2,B2 ! A,B coefficient interpolation values

   real (r8), parameter, dimension(nchl) :: &!chl for table look-up
     chlcnc = (/                            &
     .001_r8, .005_r8, .01_r8,  .02_r8,     &
     .03_r8,  .05_r8,  .10_r8,  .15_r8,     &
     .20_r8,  .25_r8,  .30_r8,  .35_r8,     &
     .40_r8,  .45_r8,  .50_r8,  .60_r8,     &
     .70_r8,  .80_r8,  .90_r8, 1.00_r8,     &
     1.50_r8, 2.00_r8, 2.50_r8, 3.00_r8,    &
     4.00_r8, 5.00_r8, 6.00_r8, 7.00_r8,    &
     8.00_r8, 9.00_r8,10.00_r8  /)

   real (r8), parameter, dimension(nchl) :: &
      A_1 = (/                              &
      0.4421_r8, 0.4451_r8, 0.4488_r8,      &
                            0.4563_r8,      &
      0.4622_r8, 0.4715_r8, 0.4877_r8,      &
                            0.4993_r8,      &
      0.5084_r8, 0.5159_r8, 0.5223_r8,      &
                            0.5278_r8,      &
      0.5326_r8, 0.5369_r8, 0.5408_r8,      &
                            0.5474_r8,      &
      0.5529_r8, 0.5576_r8, 0.5615_r8,      &
                            0.5649_r8,      &
      0.5757_r8, 0.5802_r8, 0.5808_r8,      &
                            0.5788_r8,      &
      0.56965_r8,0.55638_r8,0.54091_r8,     &
                            0.52442_r8,     &
      0.50766_r8,0.49110_r8,0.47505_r8  /)

   real (r8), parameter, dimension(nchl) :: &
      A_2 = (/                              &
      0.2981_r8, 0.2963_r8, 0.2940_r8,      &
                            0.2894_r8,      &
      0.2858_r8, 0.2800_r8, 0.2703_r8,      &
                            0.2628_r8,      &
      0.2571_r8, 0.2523_r8, 0.2481_r8,      &
                            0.2444_r8,      &
      0.2411_r8, 0.2382_r8, 0.2356_r8,      &
                            0.2309_r8,      &
      0.2269_r8, 0.2235_r8, 0.2206_r8,      &
                            0.2181_r8,      &
      0.2106_r8, 0.2089_r8, 0.2113_r8,      &
                            0.2167_r8,      &
      0.23357_r8,0.25504_r8,0.27829_r8,     &
                            0.30274_r8,     &
      0.32698_r8,0.35056_r8,0.37303_r8 /)

   real (r8), parameter, dimension(nchl) :: &
      B_1 = (/                              &
      0.0287_r8, 0.0301_r8, 0.0319_r8,      &
                            0.0355_r8,      &
      0.0384_r8, 0.0434_r8, 0.0532_r8,      &
                            0.0612_r8,      &
      0.0681_r8, 0.0743_r8, 0.0800_r8,      &
                            0.0853_r8,      &
      0.0902_r8, 0.0949_r8, 0.0993_r8,      &
                            0.1077_r8,      &
      0.1154_r8, 0.1227_r8, 0.1294_r8,      &
                            0.1359_r8,      &
      0.1640_r8, 0.1876_r8, 0.2082_r8,      &
                            0.2264_r8,      &
      0.25808_r8,0.28498_r8,0.30844_r8,     &
                           0.32932_r8,      &
      0.34817_r8,0.36540_r8,0.38132_r8 /)

   real (r8), parameter, dimension(nchl) ::  &
      B_2 = (/                               &
      0.3192_r8, 0.3243_r8, 0.3306_r8,       &
      0.3433_r8,                             &
      0.3537_r8, 0.3705_r8, 0.4031_r8,       &
                            0.4262_r8,       &
      0.4456_r8, 0.4621_r8, 0.4763_r8,       &
                            0.4889_r8,       &
      0.4999_r8, 0.5100_r8, 0.5191_r8,       &
                            0.5347_r8,       &
      0.5477_r8, 0.5588_r8, 0.5682_r8,       &
                            0.5764_r8,       &
      0.6042_r8, 0.6206_r8, 0.6324_r8,       &
                            0.6425_r8,       &
      0.66172_r8,0.68144_r8,0.70086_r8,      &
                            0.72144_r8,      &
      0.74178_r8,0.76190_r8,0.78155_r8 /) 

   integer (int_kind) :: &
      tavg_QSW_HTP,     & ! tavg id for QSW_HTP (solar short-wave heat flux in top layer)
      tavg_QSW_3D         ! tavg id for 3D QSW at top of cell

!-----------------------------------------------------------------------
!     named field indices
!-----------------------------------------------------------------------

      integer (kind=int_kind) :: chl_nf_ind = 0

!EOC
!***********************************************************************

   contains

!***********************************************************************

   subroutine init_sw_absorption
   
!-----------------------------------------------------------------------
!
!     initialize shortwave absorption.
!
!-----------------------------------------------------------------------

   use named_field_mod, only: named_field_get_index

   implicit none

   integer (int_kind) :: &
      k,n,               &! level index
      nu,                &! unit for input dataset
      nml_error,         &! namelist error flag
      bid,               &! local block address for this block
      iblock

   real (r8), allocatable, dimension(:,:,:,:) :: &
      TEMP_DATA   ! temporary data array
 
   type (datafile) ::        &
      forcing_file            !data file structure for input forcing file

   type (io_field_desc) ::   &
      io_chl                  ! io field descriptor for input sss field
             
   type (io_dim) ::          &
      i_dim, j_dim,          &! dimension descriptors for horiz dimensions
      month_dim               ! dimension descriptor  for monthly data

   namelist /sw_absorption_nml/ &
        sw_absorption_type,     &! 'top-layer', 'jerlov' or 'chlorophyll'
        jerlov_water_type,      &! jerlov water type from 1 to 5
        chl_option,             &! chlorophyll option ['file', 'model']
        chl_filename,           &! include local filepath in name
        chl_file_fmt             ! include local filepath in name

!-----------------------------------------------------------------------

   if (.not. registry_match('init_passive_tracers')) then
      call exit_POP(sigAbort, 'init_passive_tracers not called ' /&
         &/ 'before init_sw_absorption. This is necessary to ' /&
         &/ 'allow for registry of model_chlorophyll.')
   end if

!-----------------------------------------------------------------------
!
!     read shortwave absorption namelist input after setting default values.
!
!-----------------------------------------------------------------------

 
   sw_absorption_type   = 'jerlov'
   jerlov_water_type    =    3
   chl_option           = 'file'
   chl_filename         = 'unknown-chl'
   chl_file_fmt         = 'bin'

   if (my_task == master_task) then
   open (nml_in, file=nml_filename, status='old',iostat=nml_error)
   if (nml_error /= 0) then
      nml_error = -1
   else
      nml_error =  1
   endif
   do while (nml_error > 0)
      read(nml_in, nml=sw_absorption_nml,iostat=nml_error)
   end do
   if (nml_error == 0) close(nml_in)
   endif
 

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
     call exit_POP(sigAbort,'ERROR reading sw_absorption_nml')
   endif

   if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,ndelim_fmt)
       write(stdout,blank_fmt)
       write(stdout,*) ' Short-wave absorption:'
       write(stdout,blank_fmt)
       write(stdout,*) ' sw_absorption namelist  '
       write(stdout,blank_fmt)
       write(stdout, sw_absorption_nml)
       write(stdout,blank_fmt)
   endif

   call broadcast_scalar(sw_absorption_type,    master_task)
   call broadcast_scalar(jerlov_water_type,     master_task)
   call broadcast_scalar(chl_option,            master_task)
   call broadcast_scalar(chl_filename,          master_task)

   if (sw_absorption_type .ne. 'top-layer'.and.  &
          sw_absorption_type .ne. 'jerlov'.and.  &
          sw_absorption_type .ne. 'chlorophyll' ) then
     call exit_POP(sigAbort,'ERROR sw_absorption_type unknown')
   endif  

   if (sw_absorption_type .eq. 'chlorophyll') then
      if (chl_option .ne. 'file' .and. chl_option .ne. 'model') then
         call exit_POP(sigAbort,'ERROR chl_option unknown')
      endif
      if (chl_option .eq. 'model') then
         call named_field_get_index('model_chlorophyll', chl_nf_ind, &
            exit_on_err=.false.)
         if (chl_nf_ind == 0) then
            call exit_POP(sigAbort,'chl_option==model, but ' /&
               &/ 'model_chlorophyll is not registered')
         endif
      endif
   endif

!-----------------------------------------------------------------------
!     define shortwave solar absorption model.
!-----------------------------------------------------------------------

   select case (sw_absorption_type)

   case ('top-layer')

     sw_absorb(0)    = c1
     do k = 1, km
       sw_absorb(k) = c0
     enddo

   case ('jerlov')

     sw_absorb(0)  = c1
     sw_absorb(km) = c0
     do k = 1, km-1
       call sw_absorb_frac(zw(k),sw_absorb(k))
     enddo

   case ('chlorophyll')

     call set_chl_trn   ! compute transmission table

     if (chl_option == 'file') then
!-----------------------------------------------------------------------
!
!   Read in chlorophyll data
!
!-----------------------------------------------------------------------
 
     chl_bndy_loc  = field_loc_center
     chl_bndy_type = field_type_scalar
     chl_data_name = 'CHL'

     allocate( CHL_DATA (nx_block,ny_block,max_blocks_clinic,12))
     allocate( TEMP_DATA(nx_block,ny_block,12,max_blocks_clinic))

     CHL_DATA = c0

     forcing_file = construct_file(chl_file_fmt,                 &
                                   full_name=trim(chl_filename), &
                                   record_length = rec_type_dbl,  &
                                   recl_words=nx_global*ny_global)

      call data_set(forcing_file,'open_read')

      i_dim     = construct_io_dim('i',nx_global)
      j_dim     = construct_io_dim('j',ny_global)
      month_dim = construct_io_dim('month',12)

      io_chl = construct_io_field( &
                    trim(chl_data_name),                        &
                    dim1=i_dim, dim2=j_dim, dim3=month_dim,     &
                    field_loc  = chl_bndy_loc,                  &
                    field_type = chl_bndy_type,                 &
                    d3d_array=TEMP_DATA(:,:,:,:))
      call data_set(forcing_file,'define',io_chl)
      call data_set(forcing_file,'read'  ,io_chl)
      call destroy_io_field(io_chl)


      !*** re-order data 

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
         do n=1,12
            CHL_DATA(:,:,iblock,n) = TEMP_DATA(:,:,n,iblock)
         end do
      end do
      !$OMP END PARALLEL DO

      deallocate(TEMP_DATA)
      call data_set(forcing_file,'close')
      call destroy_file(forcing_file)


     if (my_task.eq.master_task) then
       write(stdout,blank_fmt)
       write(stdout,*) ' CHL monthly file read: ',chl_filename
     endif

!-----------------------------------------------------------------------
!
!   Set initial chlorophyll amount
!
!-----------------------------------------------------------------------

     CHL(:,:,:) = CHL_DATA(:,:,:,imonth) ! constant across month
     CHL = max(CHL,chlmin)               ! set lowest allowed limit
     CHL = min(CHL,chlmax)               ! set highest allowed limit
     CHLINDX = log10(CHL/chlmin)/dlogchl ! compute chlorphyll index
     CHLINDX = max(CHLINDX,0)            ! minimum limit for chl index
     CHLINDX = min(CHLINDX,nsub)         ! maximum limit for chl index

     end if ! chl_option == 'file'

   end select

   call define_tavg_field(tavg_QSW_HTP,'QSW_HTP',2,                    &
                          long_name='Solar Short-Wave Heat Flux in top layer', &
                          units='watt/m^2', grid_loc='2110',           &
                          coordinates='TLONG TLAT time')

   call define_tavg_field(tavg_QSW_3D,'QSW_3D',3,                      &
                          long_name='Solar Short-Wave Heat Flux', &
                          units='watt/m^2', grid_loc='3112',           &
                          coordinates='TLONG TLAT z_w_top time')

!-----------------------------------------------------------------------

   end subroutine init_sw_absorption

!***********************************************************************

   subroutine set_chl

   use named_field_mod, only: named_field_get

   logical (kind=log_kind) :: first_call = .true.

   integer (int_kind) ::     &
      iblock

!-----------------------------------------------------------------------
!
!  check that set_sflux_passive_tracers has been called
!
!-----------------------------------------------------------------------

   if (first_call .and. chl_option == 'model') then
      if (.not. registry_match('set_sflux_passive_tracers')) then
         call exit_POP(sigAbort, 'set_sflux_passive_tracers not ' /&
            &/ 'called before set_chl. This is necessary for ' /&
            &/ 'exact restart when Chl is prognostic.')
      end if
   endif
 
!-----------------------------------------------------------------------
!
!   Update chlorophyll amount; update every new month
!
!-----------------------------------------------------------------------

   if( sw_absorption_type .eq. 'chlorophyll' ) then
      if (imonth .ne. imonth_last .or. &    ! if new month, update chl 
          chl_option == 'model' ) then      ! update every timestep for 'model'

         !$OMP PARALLEL DO PRIVATE(iblock)
         do iblock=1,nblocks_clinic
            if (chl_option == 'file') then
               CHL(:,:,iblock) = CHL_DATA(:,:,iblock,imonth)            ! constant across month
            else
               call named_field_get(chl_nf_ind, iblock, CHL(:,:,iblock))
            endif
   
            CHL(:,:,iblock) = max(CHL(:,:,iblock),chlmin)               ! set lowest allowed limit
            CHL(:,:,iblock) = min(CHL(:,:,iblock),chlmax)               ! set highest allowed limit
            CHLINDX(:,:,iblock) = log10(CHL(:,:,iblock)/chlmin)/dlogchl ! compute chlorphyll index
            CHLINDX(:,:,iblock) = max(CHLINDX(:,:,iblock),0)            ! minimum limit for chl index
            CHLINDX(:,:,iblock) = min(CHLINDX(:,:,iblock),nsub)         ! maximum limit for chl index
         end do
         !$OMP END PARALLEL DO

      endif
   endif

   first_call = .false.

   end subroutine set_chl

!***********************************************************************

   subroutine set_chl_trn

!-----------------------------------------------------------------------
!
!   Set chlorophyll transmissions for table look-up.
!
!   Look-up table coefficients computed over the range .001
!     to 10 mg/m^3, by a constant step in log(chl). We
!     choose to have 100 steps per decade, so increment
!     in chl is around 2.3% .
!
!   Note that table coefficients require depth in m, so cm to m
!   unit conversion is done where appropriate.
!
!-----------------------------------------------------------------------

   implicit none

   real (r8) ::               &
      chlamnt,                &! chlorophyll amount (density) in mg/m^3
      w1,                     &! interpolation weight
      w2                       ! interpolation weight

   real (r8) ::               &
      logchl,                 &! log of chlorophyll
      zprnt                    ! z level for print

   integer (int_kind) ::      &
      n,                      &! index for absorber amount
      m,                      &! index for absorber amount
      mc,                     &! specific absorber index
      k,                      &! vertical index
      kint,                   &! vertical index for layer interface
      kmid                     ! vertical index for layer mid-point

   real (r8) ::               &
      Trprnt,                 &! chlorophyll transmissions for diagnostic print
      Trjrlv,                 &! Jerlov transmissions for diagnostic print
      arg                      ! argument for exponential

   real (r8), parameter ::    &
      maxarg = 35.0_r8         ! maximum allowed exponential argument

   logical (log_kind)      :: &
      prnt = .false.           ! if true, generates diagnostic prints

!-----------------------------------------------------------------------
!
!   Diagnostic prints if desired
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then

    if (prnt) then
    write(stdout,5)
 5     format(                                                          &
       ' Transmissions for selected chlorphyll amounts'/                &
       ' Interface   Mid-point  Transmission   Jerlov Tr (type 3 = IB)') 
    do n=1,nchl
     if (n .eq.  1 .or. n .eq.  3 .or.          &
         n .eq.  6 .or. n .eq.  7 .or.          &
         n .eq. 15 .or. n .eq. 20 .or.          &
         n .eq. 24 .or. n .eq. 31) then
     write(stdout,8) chlcnc(n) 
 8   format(/' chlorophyll amount = ',f6.3/)
     do k=1,km

       kint = k-1
       if (kint.eq.0 ) then
      Trprnt = c1
       else
      arg = min(B_1(n)*zw(kint)*mpercm,maxarg)
      Trprnt = A_1(n)*exp(-arg)
      arg = min(B_2(n)*zw(kint)*mpercm,maxarg)
      Trprnt = Trprnt + A_2(n)*exp(-arg)
       endif
       arg = min(zw(kint)*mpercm/1.0_r8,maxarg)
       Trjrlv = .67_r8*exp(-arg)
       arg = min(zw(kint)*mpercm/17.0_r8,maxarg)
       Trjrlv = Trjrlv + .33_r8*exp(-arg)
       if (kint .eq. 0 ) then
      zprnt  = c0
      Trjrlv = c1
      if (Trprnt .gt. .00005 )  &
             write(stdout,10) kint,zprnt,Trprnt,Trjrlv
 10      format(1x,i3,1x,f6.1,15x,f6.4,8x,f6.4)
       else
      if (Trprnt .gt. .00005 )  &
             write(stdout,11) kint,zw(kint)*mpercm,Trprnt,Trjrlv
 11      format(1x,i3,1x,f6.1,16x,f5.4,8x,f5.4)
       endif

       if (k .le. km-1 ) then
      kmid = kint + 1

      arg = min(B_1(n)*zt(kmid)*mpercm,maxarg)
      Trprnt = A_1(n)*exp(-arg)
      arg = min(B_2(n)*zt(kmid)*mpercm,maxarg)
      Trprnt = Trprnt + A_2(n)*exp(-arg)

      arg = min(zt(kmid)*mpercm/1.0_r8,maxarg)
      Trjrlv = .67_r8*exp(-arg)
      arg = min(zt(kmid)*mpercm/17.0_r8,maxarg)
      Trjrlv = Trjrlv + .33_r8*exp(-arg)

      if (Trprnt .gt. .00005 ) &
             write(stdout,15) kmid,zt(kmid)*mpercm,Trprnt,Trjrlv
 15      format(12x,i3,1x,f6.1,5x,f5.4,8x,f5.4)
       endif
     enddo
     endif
    enddo
    endif

   endif

!-----------------------------------------------------------------------
!
!   Compute transmission array
!
!-----------------------------------------------------------------------

   ztr(0) = c0
   do k=1,km
     ztr(2*k-1) = zt(k)
     ztr(2*k)   = zw(k)
   enddo

   chlmin  = chlcnc(1)
   chlmax  = chlcnc(nchl) 
   dlogchl = (log10(chlmax)-log10(chlmin))/real(nsub)

   if (my_task == master_task) then
     if (prnt) then 
       write(stdout,20) chlmin,chlmax,nsub,dlogchl
 20    format(/1x,' minimum chlorophyll  = ',1pe11.4,' mg/m^3'/  &
            1x,' maximum chlorophyll  = ',1pe11.4,' mg/m^3'/     &
            1x,' number of intervals  = ',i5/                    &
            1x,' log chlorophyll step = ',1pe11.4)
     endif
   endif

   logchl = log10(chlmin) - dlogchl
   do n=0,nsub
     logchl = logchl + dlogchl
     chlamnt = 10**(logchl)
     do m=1,nchl-1
       if (chlcnc(m) .le. chlamnt .and. &
           chlamnt .le. chlcnc(m+1) ) then
      mc = m
      goto 30
       endif
     enddo
     if (my_task == master_task) then
       write(stdout,25) chlamnt
 25    format(/' Could not find range for chlamnt = ',1pe11.4)
       call exit_POP(sigAbort,' set_chl_trans range error for chlamnt ')
     endif
 30     continue

     w2 = (chlamnt-chlcnc(mc))/(chlcnc(mc+1)-chlcnc(mc))
     w1 = c1 - w2
     A1 = A_1(mc)*w1 + A_1(mc+1)*w2
     A2 = A_2(mc)*w1 + A_2(mc+1)*w2
     B1 = B_1(mc)*w1 + B_1(mc+1)*w2
     B2 = B_2(mc)*w1 + B_2(mc+1)*w2
     if (my_task == master_task) then
       if (prnt) then
      write(stdout,35) n,chlamnt,A1,A2,B1,B2
 35         format(/' n, chl = ',i5,2x,1pe11.4/  &
         ' A1,A2,B1,B2 = ',4(f6.4,2x)/           &
         ' Level  Depth(m)  Transmission ')
       endif
     endif

     Tr(0,n) = c1
     do k=1,ksol
       arg = min(B1*ztr(k)*mpercm,maxarg)
       Tr(k,n) = A1*exp(-arg)
       arg = min(B2*ztr(k)*mpercm,maxarg)
       Tr(k,n) = Tr(k,n) + A2*exp(-arg)
     enddo

     if (my_task == master_task) then
       if (prnt) then
      do k=0,ksol
        if (Tr(k,n) .gt. .00005 ) &
            write(stdout,40) k,ztr(k)*mpercm,Tr(k,n)
 40        format(1x,i3,4x,f6.1,5x,f6.4)
      enddo
       endif
     endif

   enddo

   if (my_task.eq.master_task) then
     write(stdout,blank_fmt)
     write(stdout,*) ' Chlorophyll transmission table computed'
     call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
   endif


   end subroutine set_chl_trn


!***********************************************************************
!BOP
! !IROUTINE: sw_absorb_frac
! !INTERFACE:

 subroutine sw_absorb_frac( depth, sw_absorb_fraction )

! !DESCRIPTION:
!  Computes fraction of solar short-wave flux penetrating to
!  specified depth due to exponential decay in Jerlov water type.
!  Reference : two band solar absorption model of Simpson and
!     Paulson (1977)
!  Note: below 200m the solar penetration gets set to zero,
!     otherwise the limit for the exponent ($+/- 5678$) needs to be 
!     taken care of.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8) :: &
      depth     ! vertical depth (cm, >0.) for desired sw fraction

! !OUTPUT PARAMETERS:

   real (r8) :: &
     sw_absorb_fraction     ! short wave (radiation) fractional decay

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      num_water_types = 5  ! max number of different water types

   real (r8), parameter :: &
      depth_cutoff = -200.0_r8

   real (r8) :: &
      depth_neg_meters

!-----------------------------------------------------------------------
!
!   define Jerlov water properties with rfac, depth1, depth2
!     Jerlov water type :  I       IA      IB      II      III
!     jerlov_water_type :  1       2       3       4       5
!
!-----------------------------------------------------------------------

   real (r8), dimension(num_water_types) ::                       &
      rfac   = (/ 0.58_r8, 0.62_r8, 0.67_r8, 0.77_r8, 0.78_r8 /), &
      depth1 = (/ 0.35_r8, 0.60_r8, 1.00_r8, 1.50_r8, 1.40_r8 /), &
      depth2 = (/ 23.0_r8, 20.0_r8, 17.0_r8, 14.0_r8, 7.90_r8 /)

!-----------------------------------------------------------------------
!
!  compute absorption fraction
!
!-----------------------------------------------------------------------

   depth_neg_meters =  -depth*mpercm  ! convert from cm to m and
                                      ! change sign

   if (depth_neg_meters < depth_cutoff) then
      sw_absorb_fraction = c0
   else
      sw_absorb_fraction =     rfac(jerlov_water_type)*            &
                 exp(depth_neg_meters/depth1(jerlov_water_type)) +  &
                     (c1 - rfac(jerlov_water_type))*                &
                 exp(depth_neg_meters/depth2(jerlov_water_type))
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine sw_absorb_frac

!***********************************************************************
!BOP
! !IROUTINE: add_sw_absorb
! !INTERFACE:

 subroutine add_sw_absorb(T_SOURCE, SHF_QSW, k, this_block)

! !DESCRIPTION:
!  If surface short wave heat flux is available, this routine caculates
!  the flux which passes through the top layer and enters lower vertical
!  depths in the ocean.  This flux is added as a source term in the
!  baroclinic equations.
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,nt), intent(inout) :: &
      T_SOURCE     ! source terms for all tracers (to avoid copies)
                   ! sw absorption added only to other potential
                   ! temperature tracers

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k            ! vertical level index

   type (block), intent(in) :: &
      this_block   ! block info for this block
 
   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SHF_QSW      
 

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bid                ! local block index
 
   real (r8), dimension(nx_block,ny_block) :: &
      WORK    ! temporary work space
 

!-----------------------------------------------------------------------
!
!  calculate short wave absorption if available
!  absorption profile pre-calculated in init_shf routine
!
!-----------------------------------------------------------------------

   if (lsw_absorb .or. registry_match('lcoupled')) then  ! short wave flux is available

      bid = this_block%local_id

      WORK = max(SHF_QSW,c0) !*** insure no neg QSW - store in work

      select case (sw_absorption_type)
 
      case ('top-layer', 'jerlov')

        if (partial_bottom_cells) then
           where (k < KMT(:,:,bid))
              T_SOURCE(:,:,1) = T_SOURCE(:,:,1) +                      &
                                WORK*(sw_absorb(k-1) - sw_absorb(k))*  &
                                dzr(k)
           elsewhere  !  do not allow energy absorption by the ground
              T_SOURCE(:,:,1) = T_SOURCE(:,:,1) + WORK*                &
                                (sw_absorb(k-1))/DZT(:,:,k,bid)
           endwhere
        else
           where (k < KMT(:,:,bid))
              T_SOURCE(:,:,1) = T_SOURCE(:,:,1) +                      &
                                WORK*(sw_absorb(k-1) - sw_absorb(k))*  &
                                dzr(k)
           elsewhere  !  do not allow energy absorption by the ground
              T_SOURCE(:,:,1) = T_SOURCE(:,:,1) +                      &
                                WORK *(sw_absorb(k-1))*dzr(k)
           endwhere
        endif
        
        if (k == 1) then
             call accumulate_tavg_field(WORK*(sw_absorb(0)-sw_absorb(1))/hflux_factor, &
                                        tavg_QSW_HTP,bid,k)
        endif
           call accumulate_tavg_field(WORK*sw_absorb(k-1)/hflux_factor, &
                                      tavg_QSW_3D,bid,k)
 
      case ('chlorophyll')
 
        if( k.eq.1 ) TRANSKM1(:,:,bid) = c1         ! surface z=0
        call sw_trans_chl(2*k,this_block)  ! use explicit index and
                                           ! chl to get transmission
 
        if (partial_bottom_cells) then
           where (k < KMT(:,:,bid))
              T_SOURCE(:,:,1) = T_SOURCE(:,:,1) + &
                                WORK*(TRANSKM1(:,:,bid)-TRANS(:,:,bid))*dzr(k)

           elsewhere  !  do not allow energy absorption by the ground
              T_SOURCE(:,:,1) = T_SOURCE(:,:,1) + &
                                WORK*TRANSKM1(:,:,bid)/DZT(:,:,k,bid)
           endwhere
        else
           where (k < KMT(:,:,bid))
              T_SOURCE(:,:,1) = T_SOURCE(:,:,1) + WORK*(TRANSKM1(:,:,bid)- &
                                TRANS(:,:,bid))*dzr(k)
           elsewhere  !  do not allow energy absorption by the ground
              T_SOURCE(:,:,1) = T_SOURCE(:,:,1) + WORK*TRANSKM1(:,:,bid)*dzr(k)
           endwhere
        endif

        if (k == 1) then
             call accumulate_tavg_field(WORK*(TRANSKM1(:,:,bid)-TRANS(:,:,bid))/hflux_factor, &
                                        tavg_QSW_HTP,bid,k)
        endif
           call accumulate_tavg_field(WORK*TRANSKM1(:,:,bid)/hflux_factor, &
                                      tavg_QSW_3D,bid,k)
  
        TRANSKM1(:,:,bid) = TRANS(:,:,bid)
 
      end select

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine add_sw_absorb
!***********************************************************************


   subroutine sw_trans_chl(kin,this_block)

!-----------------------------------------------------------------------
!
!   Given level kin (>0) or depth (positive) and chlorophyll amount 
!     chl, find transmission from table look-up values.
!
!-----------------------------------------------------------------------

   implicit none
  
! !INPUT PARAMETERS:

   type (block), intent(in) :: &
      this_block   ! block info for this block

   integer (int_kind), intent(in) :: &
      kin                           ! input index (if > 0), or implied ZTRANS (kin = 0)  
!-----------------------------------------------------------------------
!
!   local variables
!
!-----------------------------------------------------------------------
   integer (int_kind) ::             &
      i,j,k,                         &! horizontal, vertical indices
      bid

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) ::  &
      w1,                            &! weight for level interpolation
      w2                              ! weight for level interpolation

   integer (int_kind), dimension(nx_block,ny_block) :: &
      kindx                           ! interpolation level index

   bid = this_block%local_id


!-----------------------------------------------------------------------
!
!   Explicit level known
!
!-----------------------------------------------------------------------

   if (kin .gt. 0 ) then

     do j=1,ny_block
     do i=1,nx_block
       TRANS(i,j,bid)  =  Tr(kin ,CHLINDX(i,j,bid))
     enddo
     enddo

   else

!-----------------------------------------------------------------------
!
!   Find level indices and interpolate (chlorophyll index
!   has previously been computed when chl was updated)
!
!-----------------------------------------------------------------------

     kindx = ksol - 1

!-----------------------------------------------------------------------
!   F90 compiler requires conditional to be a scalar;
!   use explicit indices
!-----------------------------------------------------------------------

     do k=1,ksol
       do j=1,ny_block
       do i=1,nx_block

       if (ztr(k-1).le.ZTRANS(i,j,bid) .and. ZTRANS(i,j,bid).lt.ztr(k) ) then
         w2(i,j,bid)  = (ZTRANS(i,j,bid)-ztr(k-1))/(ztr(k)-ztr(k-1))
         w1(i,j,bid)  = c1 - w2(i,j,bid)
         kindx(i,j) = k - 1
       endif 

      enddo 
      enddo

     enddo

!-----------------------------------------------------------------------
!   F90 compiler requires a vector subscript to be a rank 1 expression;
!   use explicit indices
!-----------------------------------------------------------------------

     do j=1,ny_block
     do i=1,nx_block
     TRANS(i,j,bid) = w1(i,j,bid)*Tr(kindx(i,j)  ,CHLINDX(i,j,bid)) +  &
                      w2(i,j,bid)*Tr(kindx(i,j)+1,CHLINDX(i,j,bid))
     enddo 
     enddo

   endif

   end subroutine sw_trans_chl

!***********************************************************************



   end module sw_absorption

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
