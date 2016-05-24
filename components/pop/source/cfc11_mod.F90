!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module cfc11_mod

!BOP
! !MODULE: cfc11_mod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$Id: cfc11_mod.F90 56176 2013-12-20 18:35:46Z mlevy@ucar.edu $

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_IOUnitsMod

   use blocks, only: nx_block, ny_block, block
   use domain_size, only: max_blocks_clinic, km, nx_global, ny_global
   use domain, only: nblocks_clinic, distrb_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use prognostic, only: tracer_field
   use kinds_mod
   use constants, only: c0, c1, char_blank, delim_fmt
   use io, only: data_set
   use io_types, only: stdout, datafile, io_field_desc, io_dim,       &
       nml_in, nml_filename, construct_file, construct_io_dim,        &
       construct_io_field, rec_type_dbl, destroy_file,                &
       destroy_io_field, get_unit, release_unit
   use io_tools, only: document
   use tavg, only: define_tavg_field, accumulate_tavg_field
   use timers, only: get_timer
   use passive_tracer_tools, only: forcing_monthly_every_ts,          &
       init_forcing_monthly_every_ts, ind_name_pair, tracer_read,     &
       rest_read_tracer_block, file_read_tracer_block, read_field

   implicit none

  !-----------------------------------------------------------------------------
  !   public/private declarations
  !-----------------------------------------------------------------------------

   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: &
       cfc11_tracer_cnt, &
       cfc11_init, &
       cfc11_set_sflux,  &
       cfc11_tavg_forcing, &
       SCHMIDT_CFC

!EOP
!BOC
!maltrud
     character(char_len) :: mmessage

  !-----------------------------------------------------------------------------
  !   module variables required by forcing_passive_tracer
  !-----------------------------------------------------------------------------

  integer(int_kind), parameter :: &
       cfc11_tracer_cnt = 1

!-----------------------------------------------------------------------
!     index bounds of passive tracer module variables in TRACER
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
     num_cfc11_years

  !-----------------------------------------------------------------------------
  !   relative tracer indices
  !-----------------------------------------------------------------------------

  integer(int_kind), parameter :: &
       cfc11_ind = 1     ! cfc11 index

  !-----------------------------------------------------------------------------
  !   derived type & parameter for tracer index lookup
  !-----------------------------------------------------------------------------

  type(ind_name_pair), dimension(cfc11_tracer_cnt) :: &
       ind_name_table = (/ ind_name_pair(cfc11_ind, 'CFC11') /)

  !-----------------------------------------------------------------------------
  !   derived type for tracer initialization
  !-----------------------------------------------------------------------------

    type(tracer_read) :: &
         gas_flux_fice,       & ! ice fraction for gas fluxes
         gas_flux_ws,         & ! wind speed for gas fluxes
         gas_flux_ap            ! atmospheric pressure for gas fluxes

  type(forcing_monthly_every_ts), save :: &
       fice_file,              & ! ice fraction, if read from file
       xkw_file,               & ! a * wind-speed ** 2, if read from file
       ap_file                   ! atmoshperic pressure, if read from file

  !-----------------------------------------------------------------------------
  !   define tavg id for 2d fields related to surface fluxes
  !-----------------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_FICE,         &! tavg id for ice fraction
      tavg_XKW,          &! tavg id for xkw
      tavg_ATM_PRESS,    &! tavg id for atmospheric pressure
      tavg_pCFC11,       &! tavg id for cfc11 partial pressure
      tavg_PV,           &! tavg id for piston velocity
      tavg_SCHMIDT_CFC11,&! tavg id for cfc11 schmidt number
      tavg_CFC11SAT       ! tavg id for cfc11 saturation

  !-----------------------------------------------------------------------------
  !   define array for holding flux-related quantities that need to be time-averaged
  !   this is necessary since the forcing routines are called before tavg flags
  !-----------------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable ::   &
      CFC11_SFLUX_TAVG

   integer (int_kind) :: cfc11_sflux_timer

   character(char_len) :: &
         cfc11_formulation,           & ! how to calculate flux (ocmip or bulk)
         cfc11_interp_weight_file,    & ! filename for hemispheric weights
         cfc11_interp_weight_file_fmt,& ! file format for weights file
         cfc11_time_series_file         ! filename for ascii time series of atm cfc11

  logical(log_kind), dimension(:,:,:), allocatable, save :: &
       LAND_MASK

   real (r8), dimension(:,:,:), allocatable :: &
      CFC11_INTERP_WEIGHT,  &! latitudinally dependent weight for atm CFC11
      pCFC11                 ! surface patial pressure of CFC11

   real (r8), dimension(:,:), allocatable :: &
      cfc11_time_series ! time series of atm CFC11

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: cfc11_init
! !INTERFACE:

 subroutine cfc11_init(init_ts_file_fmt, read_restart_filename, &
                       tracer_d_module, TRACER_MODULE, errorCode)

! !DESCRIPTION:
!  Initialize cfc11 tracer module. This involves setting metadata, reading
!  the module namelist and setting initial conditions.

! !REVISION HISTORY:
!  same as module

! !USES:

   use broadcast, only: broadcast_scalar, broadcast_array
   use constants, only: c0, field_loc_center, blank_fmt, &
       field_type_scalar
   use prognostic, only: nx_global, ny_global, curtime, oldtime
   use grid, only: KMT, zt, zw, n_topo_smooth, fill_points
   use forcing_tools, only: find_forcing_times
   use time_management, only: freq_opt_nyear, freq_opt_nmonth

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(cfc11_tracer_cnt), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real(r8), dimension(nx_block,ny_block,km,cfc11_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

    character(*), parameter :: subname = 'cfc11_mod:cfc11_init'

    character(char_len) :: &
         init_cfc11_option,           & ! option for initialization of bgc
         init_cfc11_init_file,        & ! filename for option 'file'
         init_cfc11_init_file_fmt       ! file format for option 'file'

    integer(int_kind) :: &
         n,                   & ! index for looping over tracers
         k,                   & ! index for looping over depth levels
         l,                   & ! index for looping over time levels
         ind,                 & ! tracer index for tracer name from namelist
         iblock,              & ! index for looping over blocks
         nml_error              ! namelist i/o error flag

    type(tracer_read), dimension(cfc11_tracer_cnt) :: &
         tracer_init_int,      &! namelist variable for initializing tracers
         tracer_init_ext        ! namelist variable for initializing tracers

    namelist /cfc11_nml/ &
         init_cfc11_option, init_cfc11_init_file, tracer_init_ext, &
         init_cfc11_init_file_fmt, cfc11_interp_weight_file,       &
         cfc11_interp_weight_file_fmt, cfc11_time_series_file,     &
         gas_flux_fice, gas_flux_ws, gas_flux_ap, &
         cfc11_formulation

   character (char_len) ::  &
      cfc11_restart_filename  ! modified file name for restart file

!-----------------------------------------------------------------------
!  initialize forcing_monthly_every_ts variables
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call init_forcing_monthly_every_ts(fice_file)
   call init_forcing_monthly_every_ts(xkw_file)
   call init_forcing_monthly_every_ts(ap_file)

!-----------------------------------------------------------------------
!  initialize tracer_d values
!-----------------------------------------------------------------------

   tracer_d_module(cfc11_ind)%short_name = 'CFC11'
   tracer_d_module(cfc11_ind)%long_name  = 'CFC11'
   tracer_d_module(cfc11_ind)%units      = 'fmol/cm^3'
   tracer_d_module(cfc11_ind)%tend_units = 'fmol/cm^3/s'
   tracer_d_module(cfc11_ind)%flux_units = 'fmol/cm^2/s'

    !---------------------------------------------------------------------------
    !   default namelist settings
    !---------------------------------------------------------------------------

    init_cfc11_option = 'unknown'
    init_cfc11_init_file = 'unknown'
    init_cfc11_init_file_fmt = 'bin'
    cfc11_interp_weight_file = 'unknown'
    cfc11_interp_weight_file_fmt = 'bin'
    cfc11_time_series_file = 'unknown'
    cfc11_formulation = 'unknown'

    do n = 1,cfc11_tracer_cnt
       tracer_init_ext(n)%mod_varname  = 'unknown'
       tracer_init_ext(n)%filename     = 'unknown'
       tracer_init_ext(n)%file_varname = 'unknown'
       tracer_init_ext(n)%scale_factor = c1
       tracer_init_ext(n)%default_val  = c0
       tracer_init_ext(n)%file_fmt     = 'bin'
    end do

    gas_flux_fice%filename     = 'unknown'
    gas_flux_fice%file_varname = 'FICE'
    gas_flux_fice%scale_factor = c1
    gas_flux_fice%default_val  = c0
    gas_flux_fice%file_fmt     = 'bin'

    gas_flux_ws%filename     = 'unknown'
    gas_flux_ws%file_varname = 'XKW'
    gas_flux_ws%scale_factor = c1
    gas_flux_ws%default_val  = c0
    gas_flux_ws%file_fmt     = 'bin'

    gas_flux_ap%filename     = 'unknown'
    gas_flux_ap%file_varname = 'P'
    gas_flux_ap%scale_factor = c1
    gas_flux_ap%default_val  = c0
    gas_flux_ap%file_fmt     = 'bin'

    if (my_task == master_task) then
       open (nml_in, file=nml_filename, status='old',iostat=nml_error)
       if (nml_error /= 0) then
          nml_error = -1
       else
          nml_error =  1
       endif
       do while (nml_error > 0)
          read(nml_in, nml=cfc11_nml,iostat=nml_error)
       end do
       if (nml_error == 0) close(nml_in)
    endif

    call broadcast_scalar(nml_error, master_task)
    if (nml_error /= 0) then
       call document(subname, 'cfc11_nml not found')
       mmessage = 'ERROR : stopping in '/&
                                         &/ subname
       call exit_POP(sigAbort,mmessage)
    end if

    !---------------------------------------------------------------------------
    !   broadcast all namelist variables
    !---------------------------------------------------------------------------

    call broadcast_scalar(init_cfc11_option, master_task)
    call broadcast_scalar(init_cfc11_init_file, master_task)
    call broadcast_scalar(init_cfc11_init_file_fmt, master_task)
    call broadcast_scalar(cfc11_interp_weight_file, master_task)
    call broadcast_scalar(cfc11_interp_weight_file_fmt, master_task)
    call broadcast_scalar(cfc11_time_series_file, master_task)
    call broadcast_scalar(cfc11_formulation, master_task)

    call broadcast_scalar(gas_flux_fice%filename, master_task)
    call broadcast_scalar(gas_flux_fice%file_varname, master_task)
    call broadcast_scalar(gas_flux_fice%scale_factor, master_task)
    call broadcast_scalar(gas_flux_fice%default_val, master_task)
    call broadcast_scalar(gas_flux_fice%file_fmt, master_task)

    fice_file%input = gas_flux_fice

    call broadcast_scalar(gas_flux_ws%filename, master_task)
    call broadcast_scalar(gas_flux_ws%file_varname, master_task)
    call broadcast_scalar(gas_flux_ws%scale_factor, master_task)
    call broadcast_scalar(gas_flux_ws%default_val, master_task)
    call broadcast_scalar(gas_flux_ws%file_fmt, master_task)

    xkw_file%input = gas_flux_ws

    call broadcast_scalar(gas_flux_ap%filename, master_task)
    call broadcast_scalar(gas_flux_ap%file_varname, master_task)
    call broadcast_scalar(gas_flux_ap%scale_factor, master_task)
    call broadcast_scalar(gas_flux_ap%default_val, master_task)
    call broadcast_scalar(gas_flux_ap%file_fmt, master_task)

    ap_file%input = gas_flux_ap

    do n = 1,cfc11_tracer_cnt
       call broadcast_scalar(tracer_init_ext(n)%mod_varname, master_task)
       call broadcast_scalar(tracer_init_ext(n)%filename, master_task)
       call broadcast_scalar(tracer_init_ext(n)%file_varname, master_task)
       call broadcast_scalar(tracer_init_ext(n)%scale_factor, master_task)
       call broadcast_scalar(tracer_init_ext(n)%default_val, master_task)
       call broadcast_scalar(tracer_init_ext(n)%file_fmt, master_task)
    end do

    !---------------------------------------------------------------------------
    !   initialize tracers
    !---------------------------------------------------------------------------

    select case (init_cfc11_option)

    case ('ccsm_startup', 'zero', 'ccsm_startup_spunup')
      TRACER_MODULE = c0
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d CFC11 set to all zeros'
          write(stdout,delim_fmt)
          call POP_IOUnitsFlush(POP_stdout) ; call POP_IOUnitsFlush(stdout)
      endif

    case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid' )

       cfc11_restart_filename = char_blank

       if (init_cfc11_init_file == 'same_as_TS') then
          if (read_restart_filename == 'undefined') then
             call exit_POP(sigAbort, &
                'no restart file to read cfc11 from')
          end if
          cfc11_restart_filename = read_restart_filename
          init_cfc11_init_file_fmt = init_ts_file_fmt

       else  ! do not read from TS restart file

          cfc11_restart_filename = trim(init_cfc11_init_file)

       end if

       call rest_read_tracer_block(init_cfc11_init_file_fmt, &
                                   cfc11_restart_filename,   &
                                   tracer_d_module,          &
                                   TRACER_MODULE)

    case ('file')
       call document(subname, 'cfc11 being read from separate file')

       call file_read_tracer_block(init_cfc11_init_file_fmt, &
                                   init_cfc11_init_file,     &
                                   tracer_d_module,          &
                                   ind_name_table,           &
                                   tracer_init_ext,          &
                                   TRACER_MODULE)

       if (n_topo_smooth > 0) then
          do k=1,km
             call fill_points(k,TRACER_MODULE(:,:,k,1,curtime,:), &
                              errorCode)

             if (errorCode /= POP_Success) then
                call POP_ErrorSet(errorCode, &
                   'cfc11_init: error in fill_points')
                return
             endif
          enddo
       endif

    case default
       call document(subname, 'init_cfc11_option = ', init_cfc11_option)
!      call exit_POP('ERROR: stopping in ' // subname)

    end select

    !---------------------------------------------------------------------------
    !   apply land mask to tracers
    !---------------------------------------------------------------------------

    do iblock=1,nblocks_clinic
    do n = 1,cfc11_tracer_cnt
       do k = 1,km
          where (k > KMT(:,:,iblock))
             TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
             TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
          end where
       end do
    end do
    enddo

    !---------------------------------------------------------------------------
    !   allocate and initialize CFC11 partial pressure array
    !---------------------------------------------------------------------------

    allocate( pCFC11(nx_block,ny_block,max_blocks_clinic) )
    pCFC11 = c0

    !---------------------------------------------------------------------------
    !   allocate and initialize LAND_MASK (true for ocean points)
    !---------------------------------------------------------------------------

    allocate( LAND_MASK(nx_block,ny_block,max_blocks_clinic) )
    LAND_MASK = (KMT.gt.0)

    call get_timer(cfc11_sflux_timer, 'CFC11_SFLUX',1, &
                                          distrb_clinic%nprocs)

!-------------------------------------------------------------------------------
!  call other initialization subroutines
!-------------------------------------------------------------------------------

    call init_cfc11_tavg
    call cfc11_init_sflux

!-----------------------------------------------------------------------
!EOC

 end subroutine cfc11_init

!***********************************************************************
!BOP
! !IROUTINE: init_cfc11_tavg
! !INTERFACE:

 subroutine init_cfc11_tavg

! !DESCRIPTION:
!  Define tavg fields not automatically handled by the base model.

! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      var_cnt             ! how many tavg variables are defined

!-----------------------------------------------------------------------

   var_cnt = 0

   call define_tavg_field(tavg_FICE,'CFC11_FICE',2,               &
                          long_name='CFC11 Ice Fraction',         &
                          units='fraction', grid_loc='2111')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_XKW,'CFC11_XKW',2,                 &
                          long_name='CFC11 XKW',                  &
                          units='m/s', grid_loc='2111')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_pCFC11,'pCFC11',2,                 &
                          long_name='CFC11 partial pressure',     &
                          units='atmospheres', grid_loc='2111')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_PV,'CFC11_PV',2,                   &
                          long_name='CFC11 piston velocity',      &
                          units='m/s', grid_loc='2111')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_ATM_PRESS,'CFC11_ATM_PRESS',2,     &
                          long_name='CFC11 Atmospheric Pressure', &
                          units='atmospheres', grid_loc='2111')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_SCHMIDT_CFC11,'CFC11_SCHMIDT',2,   &
                          long_name='CFC11 Schmidt Number',       &
                          units='none', grid_loc='2111')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC11SAT,'CFC11_SAT',2,            &
                          long_name='CFC11 Saturation',           &
                          units='mmol/m^3', grid_loc='2111')
   var_cnt = var_cnt+1

!-----------------------------------------------------------------------

   allocate(CFC11_SFLUX_TAVG(nx_block,ny_block,var_cnt,max_blocks_clinic))

!-----------------------------------------------------------------------
!EOC

 end subroutine init_cfc11_tavg

!***********************************************************************
!BOP
! !IROUTINE: cfc11_init_sflux
! !INTERFACE:

 subroutine cfc11_init_sflux

   use broadcast, only: broadcast_scalar, broadcast_array
   use constants, only: field_loc_center, blank_fmt, field_type_scalar
   use grid, only: KMT, zt, zw
   use forcing_tools, only: find_forcing_times

   type (datafile) ::    &
      in_file             ! data file type for init ts file

   type (io_field_desc) :: &
      io_tracer           ! io field descriptors for input Tracer

   type (io_dim) :: &
      i_dim, j_dim, k_dim, month_dim ! dimension descriptors

   real (r8), dimension(:,:,:,:), allocatable :: &
      TEMP_DATA           ! temp array for reading Tracer data

    integer(int_kind) :: &
         io_error,            & ! io error status
         nu,                  & ! io unit number
         n,                   & ! index for looping over tracers
         k,                   & ! index for looping over levels
         nc,                  & ! index for looping over columns in time series file
         ny,                  & ! index for looping over years in time series file
         iblock                 ! index for looping over blocks

    !---------------------------------------------------------------------------
    !   read gas flux forcing (if required)
    !   otherwise, use values passed in
    !---------------------------------------------------------------------------

   select case (cfc11_formulation)

     case ('ocmip')

       allocate(TEMP_DATA(nx_block,ny_block,12,max_blocks_clinic))

    !---------------------------------------------------------------------------
    !   first, read ice file
    !---------------------------------------------------------------------------

       allocate(fice_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))
       in_file = construct_file(fice_file%input%file_fmt,             &
                         full_name=trim(fice_file%input%filename),    &
                         record_length = rec_type_dbl,                &
                         recl_words=nx_global*ny_global)
       call data_set(in_file,'open_read')

       i_dim = construct_io_dim('i',nx_global)
       j_dim = construct_io_dim('j',ny_global)
       month_dim = construct_io_dim('month',12)

       io_tracer = &
           construct_io_field(trim(fice_file%input%file_varname), &
           dim1=i_dim, dim2=j_dim, dim3=month_dim,                &
           field_loc = field_loc_center,                          &
           field_type = field_type_scalar,                        &
           d3d_array=TEMP_DATA)

       call data_set(in_file,'define',io_tracer)
       call data_set(in_file,'read'  ,io_tracer)
       call destroy_io_field(io_tracer)

!untested       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1,nblocks_clinic
       do n=1,12
          fice_file%DATA(:,:,iblock,1,n) = &
             TEMP_DATA(:,:,n,iblock)*fice_file%input%scale_factor
          where (.not. LAND_MASK(:,:,iblock))
             fice_file%DATA(:,:,iblock,1,n) = c0
          end where
       end do
       end do
!untested       !$OMP END PARALLEL DO

       call data_set(in_file,'close')
       call destroy_file(in_file)

       if (my_task == master_task) then
          write(stdout,blank_fmt)
          write(stdout,'(a12,a)') ' file read: ', &
             trim(fice_file%input%filename)
       endif

       call find_forcing_times(fice_file%data_time, &
            fice_file%data_inc, fice_file%interp_type, &
            fice_file%data_next, fice_file%data_time_min_loc, &
            fice_file%data_update, fice_file%data_type)

    !---------------------------------------------------------------------------
    !   next, read piston velocity file
    !---------------------------------------------------------------------------

       allocate(xkw_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))
       in_file = construct_file(xkw_file%input%file_fmt,                 &
                         full_name=trim(xkw_file%input%filename), &
                         record_length = rec_type_dbl,                &
                         recl_words=nx_global*ny_global)
       call data_set(in_file,'open_read')

       i_dim = construct_io_dim('i',nx_global)
       j_dim = construct_io_dim('j',ny_global)
       month_dim = construct_io_dim('month',12)

       io_tracer = &
           construct_io_field(trim(xkw_file%input%file_varname), &
           dim1=i_dim, dim2=j_dim, dim3=month_dim,               &
           field_loc = field_loc_center,                         &
           field_type = field_type_scalar,                       &
           d3d_array=TEMP_DATA)

       call data_set(in_file,'define',io_tracer)
       call data_set(in_file,'read'  ,io_tracer)
       call destroy_io_field(io_tracer)

!untested       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1,nblocks_clinic
       do n=1,12
          xkw_file%DATA(:,:,iblock,1,n) = &
             TEMP_DATA(:,:,n,iblock)*xkw_file%input%scale_factor
          where (.not. LAND_MASK(:,:,iblock))
             xkw_file%DATA(:,:,iblock,1,n) = c0
          end where
       end do
       end do
!untested       !$OMP END PARALLEL DO

       call data_set(in_file,'close')
       call destroy_file(in_file)

       if (my_task == master_task) then
          write(stdout,blank_fmt)
          write(stdout,'(a12,a)') ' file read: ', &
             trim(xkw_file%input%filename)
       endif

       call find_forcing_times(xkw_file%data_time, &
            xkw_file%data_inc, xkw_file%interp_type, &
            xkw_file%data_next, xkw_file%data_time_min_loc, &
            xkw_file%data_update, xkw_file%data_type)

    !---------------------------------------------------------------------------
    !   last, read atmospheric pressure file
    !---------------------------------------------------------------------------

       allocate(ap_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))
       in_file = construct_file(ap_file%input%file_fmt,                &
                         full_name=trim(ap_file%input%filename),       &
                         record_length = rec_type_dbl,                 &
                         recl_words=nx_global*ny_global)
       call data_set(in_file,'open_read')

       i_dim = construct_io_dim('i',nx_global)
       j_dim = construct_io_dim('j',ny_global)
       month_dim = construct_io_dim('month',12)

       io_tracer = &
           construct_io_field(trim(ap_file%input%file_varname), &
           dim1=i_dim, dim2=j_dim, dim3=month_dim,              &
           field_loc = field_loc_center,                        &
           field_type = field_type_scalar,                      &
           d3d_array=TEMP_DATA)

       call data_set(in_file,'define',io_tracer)
       call data_set(in_file,'read'  ,io_tracer)
       call destroy_io_field(io_tracer)

!untested       !$OMP PARALLEL DO PRIVATE(iblock, n)
       do iblock=1,nblocks_clinic
       do n=1,12
          ap_file%DATA(:,:,iblock,1,n) = &
             TEMP_DATA(:,:,n,iblock)*ap_file%input%scale_factor
          where (.not. LAND_MASK(:,:,iblock))
             ap_file%DATA(:,:,iblock,1,n) = c0
          end where
       end do
       end do
!untested       !$OMP END PARALLEL DO

       call data_set(in_file,'close')
       call destroy_file(in_file)

       if (my_task == master_task) then
          write(stdout,blank_fmt)
          write(stdout,'(a12,a)') ' file read: ', &
             trim(ap_file%input%filename)
       endif

       call find_forcing_times(ap_file%data_time, &
            ap_file%data_inc, ap_file%interp_type, &
            ap_file%data_next, ap_file%data_time_min_loc, &
            ap_file%data_update, ap_file%data_type)

       deallocate(TEMP_DATA)

     case ('bulk')

       if (my_task == master_task) then
          write(stdout,*)  &
             ' Using fields from bulk forcing for calculating CFC11 flux'
       endif

     case default
      call exit_POP(sigAbort, &
                    'cfc11_init_sflux: Unknown value for cfc11_formulation')

   end select

    !---------------------------------------------------------------------------
    !   now read latitudinally dependent weights file
    !---------------------------------------------------------------------------

       allocate(CFC11_INTERP_WEIGHT(nx_block,ny_block,max_blocks_clinic))
       in_file = construct_file(cfc11_interp_weight_file_fmt,                 &
                         full_name=trim(cfc11_interp_weight_file), &
                         record_length = rec_type_dbl,                &
                         recl_words=nx_global*ny_global)
       call data_set(in_file,'open_read')

       i_dim = construct_io_dim('i',nx_global)
       j_dim = construct_io_dim('j',ny_global)

       io_tracer = &
           construct_io_field(trim(ap_file%input%file_varname), &
           dim1=i_dim, dim2=j_dim,                              &
           field_loc = field_loc_center,                        &
           field_type = field_type_scalar,                      &
           d2d_array=CFC11_INTERP_WEIGHT)

       call data_set(in_file,'define',io_tracer)
       call data_set(in_file,'read'  ,io_tracer)
       call destroy_io_field(io_tracer)
       call data_set(in_file,'close')
       call destroy_file(in_file)

       if (my_task == master_task) then
          write(stdout,blank_fmt)
          write(stdout,'(a12,a)') ' file read: ', &
             trim(cfc11_interp_weight_file)
       endif

    !---------------------------------------------------------------------------
    !   now read in hemispherically averaged CFC11 time series file
    !   first, read until end of file to determine how big to allocate
    !      the time series array
    !---------------------------------------------------------------------------

      call get_unit(nu)
      if (my_task == master_task) then
        write(stdout,*) 'Trying to read CFC11 time series from file:',  &
                           trim(cfc11_time_series_file)
        open(nu, file=cfc11_time_series_file, status='old', iostat= io_error)
        if (io_error /= 0) then
           io_error = -1
        else
           io_error =  1
        endif
        allocate( cfc11_time_series(1,3) )   ! do this temporarily
        num_cfc11_years = 0
        do while (io_error >= 0)
           read(nu,*,iostat = io_error) ( cfc11_time_series(1,nc), nc = 1, 3 )
           num_cfc11_years = num_cfc11_years + 1
        end do
        num_cfc11_years = num_cfc11_years - 1
        close (nu)
        write(stdout,*)num_cfc11_years,' years to be read from file: ',  &
                       trim(cfc11_time_series_file)
        deallocate( cfc11_time_series )
      endif

    !---------------------------------------------------------------------------
    !   check that num_cfc11_years is > 0
    !---------------------------------------------------------------------------

    call broadcast_scalar(num_cfc11_years, master_task)
    if (num_cfc11_years == 0) then
       call exit_POP(sigAbort,'cfc11_init_sflux: num_cfc11_years = 0')
    end if

    !---------------------------------------------------------------------------
    !   allocate time series array on all processors, then read in on
    !      master_task, then broadcast to all
    !   time series file has 3 columns:  year, Northern Hemisphere, Southern Hem
    !---------------------------------------------------------------------------

      allocate( cfc11_time_series(num_cfc11_years,3) )

      if (my_task == master_task) then
        open(nu, file=cfc11_time_series_file, status='old', iostat= io_error)
        do ny = 1, num_cfc11_years
           read(nu,*,iostat = io_error) ( cfc11_time_series(ny,nc), nc = 1, 3 )
        end do
        close (nu)
        write(stdout,*)' done reading file: ', trim(cfc11_time_series_file)
      endif
      call release_unit(nu)

      call broadcast_array (cfc11_time_series, master_task)

!-----------------------------------------------------------------------
!EOC

 end subroutine cfc11_init_sflux

!***********************************************************************
!BOP
! !IROUTINE: cfc11_set_sflux
! !INTERFACE:

 subroutine cfc11_set_sflux(WIND_VEL,IFRAC,PRESS,SST,SSS,SURF_VALS,  &
                              STF_MODULE)

! !DESCRIPTION:
!  Compute CFC11 surface flux and store related tavg fields for
!  subsequent accumulating.

! !REVISION HISTORY:
!  same as module

! !USES:

   use constants, only: field_loc_center, field_type_scalar, &
       mpercm, eps, p5, cmperm
   use time_management, only: thour00, check_time_flag, iyear,  &
       seconds_this_year, seconds_in_year, tday
   use broadcast, only: broadcast_scalar
   use forcing_tools, only: update_forcing_data, interpolate_forcing

!-----------------------------------------------------------------------
!   NOTE that variables are a mish-mash of model/non-model units
!   WIND_VEL should be m/s
!   XKW should be cm/s
!   salinity should be psu
!   pressure should be dyne/cm**2 (NOT atmospheres or Pascals)
!-----------------------------------------------------------------------

! !INPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block,2,max_blocks_clinic), intent(in) :: &
         WIND_VEL ! surface wind velocity (m/s)

   real(r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
         IFRAC, & ! sea ice fraction (non-dimensional)
         PRESS, & ! sea level atmospheric pressure (dyne/cm**2)
         SST,   & ! sea surface temperature (C)
         SSS      ! sea surface salinity (psu)

   real(r8), dimension(nx_block,ny_block,cfc11_tracer_cnt,max_blocks_clinic), &
         intent(in) :: SURF_VALS ! module tracers

! !OUTPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block,cfc11_tracer_cnt,max_blocks_clinic), &
         intent(inout) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

    character(*), parameter :: subname = 'cfc11_mod:cfc11_set_sflux'

    integer(int_kind) :: &
         j, iblock,    & ! loop indices
         cfc11_first_year, cfc11_last_year  !  first, last year in cfc11 time series

    real(r8) :: &  !  used for linear interpolation of time series values
          pcfc11_north, pcfc11_south, del_year

    real(r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
         FICE_USED,    & ! used ice fraction (non-dimensional)
         XKW,          & ! part of piston velocity (cm/s)
         AP_USED         ! used atm pressure (converted from dyne/cm**2 to atm)

    real(r8), dimension(nx_block,ny_block) :: &
         CFC11_loc,    & ! local copy of surface cfc11, forced to be non-negative
         XKW_ICE,      & ! common portion of piston vel., (1-fice)*xkw (cm/s)
         SCHMIDT_USED, & ! used Schmidt number
         PV,           & ! piston velocity (cm/s)
         CFC11SAT_1atm,& ! calculated CFC11 saturation (mmol/m^3)
         CFC11SAT_USED,& ! used CFC11 saturation (mmol/m^3)
         FLUX            ! tracer flux (nmol/cm^2/s)

   character (char_len) :: &
      tracer_data_label          ! label for what is being updated

   character (char_len), dimension(:), allocatable :: &
      tracer_data_names          ! short names for input data fields

   integer (int_kind), dimension(:), allocatable :: &
      tracer_bndy_loc,          &! location and field type for ghost
      tracer_bndy_type           !    cell updates

   real (r8), dimension(:,:,:,:), allocatable :: &
      TEMP_DATA           ! temp array for reading Tracer data

    real(r8), dimension(nx_block,ny_block) :: &
         WORK1, WORK2 ! temporaries for averages

    !---------------------------------------------------------------------------
    !   local parameters
    !---------------------------------------------------------------------------

    real(r8), parameter :: &
         xkw_coeff = 8.6e-9_r8       ! xkw_coeff = 0.31 cm/hr s^2/m^2 in (s/cm)

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    STF_MODULE = c0

    !---------------------------------------------------------------------------
    !   Interpolate gas flux forcing data if necessary
    !---------------------------------------------------------------------------

    select case (cfc11_formulation)

    case ('ocmip')

       if (thour00 >= fice_file%data_update) then
          allocate(tracer_data_names(1), &
                   tracer_bndy_loc  (1), &
                   tracer_bndy_type (1))
          tracer_data_names = fice_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Ice Fraction'
          call update_forcing_data(          fice_file%data_time,      &
               fice_file%data_time_min_loc,  fice_file%interp_type,    &
               fice_file%data_next,          fice_file%data_update,    &
               fice_file%data_type,          fice_file%data_inc,       &
               fice_file%DATA(:,:,:,:,1:12), fice_file%data_renorm,    &
               tracer_data_label,            tracer_data_names,        &
               tracer_bndy_loc,              tracer_bndy_type,         &
               fice_file%filename,           fice_file%input%file_fmt)
          deallocate(tracer_data_names, &
                     tracer_bndy_loc  , &
                     tracer_bndy_type )
       end if
       allocate(TEMP_DATA(nx_block,ny_block,max_blocks_clinic,1))
       call interpolate_forcing(TEMP_DATA, &
            fice_file%DATA(:,:,:,:,1:12), &
            fice_file%data_time,         fice_file%interp_type, &
            fice_file%data_time_min_loc, fice_file%interp_freq, &
            fice_file%interp_inc,        fice_file%interp_next, &
            fice_file%interp_last,       0)
       FICE_USED = TEMP_DATA(:,:,:,1)
       deallocate(TEMP_DATA)

       if (thour00 >= xkw_file%data_update) then
          allocate(tracer_data_names(1), &
                   tracer_bndy_loc  (1), &
                   tracer_bndy_type (1))
          tracer_data_names = xkw_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Piston Velocity'
          call update_forcing_data(         xkw_file%data_time,      &
               xkw_file%data_time_min_loc,  xkw_file%interp_type,    &
               xkw_file%data_next,          xkw_file%data_update,    &
               xkw_file%data_type,          xkw_file%data_inc,       &
               xkw_file%DATA(:,:,:,:,1:12), xkw_file%data_renorm,    &
               tracer_data_label,           tracer_data_names,        &
               tracer_bndy_loc,             tracer_bndy_type,         &
               xkw_file%filename,           xkw_file%input%file_fmt)
          deallocate(tracer_data_names, &
                     tracer_bndy_loc  , &
                     tracer_bndy_type )
       end if
       allocate(TEMP_DATA(nx_block,ny_block,max_blocks_clinic,1))
       call interpolate_forcing(TEMP_DATA,     &
            xkw_file%DATA(:,:,:,:,1:12), &
            xkw_file%data_time,         xkw_file%interp_type, &
            xkw_file%data_time_min_loc, xkw_file%interp_freq, &
            xkw_file%interp_inc,        xkw_file%interp_next, &
            xkw_file%interp_last,       0)
       XKW = TEMP_DATA(:,:,:,1)
       deallocate(TEMP_DATA)

       if (thour00 >= ap_file%data_update) then
          allocate(tracer_data_names(1), &
                   tracer_bndy_loc  (1), &
                   tracer_bndy_type (1))
          tracer_data_names = ap_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Atmospheric Pressure'
          call update_forcing_data(        ap_file%data_time,      &
               ap_file%data_time_min_loc,  ap_file%interp_type,    &
               ap_file%data_next,          ap_file%data_update,    &
               ap_file%data_type,          ap_file%data_inc,       &
               ap_file%DATA(:,:,:,:,1:12), ap_file%data_renorm,    &
               tracer_data_label,          tracer_data_names,        &
               tracer_bndy_loc,            tracer_bndy_type,         &
               ap_file%filename,           ap_file%input%file_fmt)
          deallocate(tracer_data_names, &
                     tracer_bndy_loc  , &
                     tracer_bndy_type )
       end if
       allocate(TEMP_DATA(nx_block,ny_block,max_blocks_clinic,1))
       call interpolate_forcing(TEMP_DATA, &
            ap_file%DATA(:,:,:,:,1:12), &
            ap_file%data_time,         ap_file%interp_type, &
            ap_file%data_time_min_loc, ap_file%interp_freq, &
            ap_file%interp_inc,        ap_file%interp_next, &
            ap_file%interp_last,       0)
       AP_USED = TEMP_DATA(:,:,:,1)
       deallocate(TEMP_DATA)

   case ('bulk')

!untested      !$OMP PARALLEL DO PRIVATE(iblock)
      do iblock = 1, nblocks_clinic

         FICE_USED(:,:,iblock) = IFRAC(:,:,iblock)
         XKW(:,:,iblock) = cmperm*cmperm * xkw_coeff*  &  ! wind is m/s, xkw is cm/s
            (WIND_VEL(:,:,1,iblock)**2 + WIND_VEL(:,:,2,iblock)**2)
         AP_USED(:,:,iblock) = PRESS(:,:,iblock)
      enddo
!untested      !$OMP END PARALLEL DO

   end select

       !------------------------------------------------------------------------
       !   assume PRESS is in cgs units (dyne/cm**2) since that is what is
       !     required for pressure forcing in barotropic
       !   want units to be atmospheres
       !   convertion from dyne/cm**2 to Pascals is P(mks) = P(cgs)/10.
       !   convertion from Pascals to atm is P(atm) = P(Pa)/101.325e+3_r8
       !
       !   Set bad AP values to 1. This is necessary for runs restarting off
       !   a run in which the flux coupler didnot restart on AP correctly.
       !------------------------------------------------------------------------

!untested   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1, nblocks_clinic

      AP_USED(:,:,iblock) = AP_USED(:,:,iblock) / 101.325e+4_r8

      AP_USED(:,:,iblock) = merge( c1, AP_USED(:,:,iblock),   &
                                   (AP_USED(:,:,iblock) > 1.5_r8 .or.   &
                                    AP_USED(:,:,iblock) < 0.5_r8) )
   enddo
!untested   !$OMP END PARALLEL DO

    !---------------------------------------------------------------------------
    !   make sure current model year is within range of cfc11 time series
    !---------------------------------------------------------------------------

    cfc11_first_year = int(cfc11_time_series(1,1) + eps)
    cfc11_last_year  = int(cfc11_time_series(num_cfc11_years,1) + eps)
    if (iyear < cfc11_first_year .or. iyear > cfc11_last_year ) then
       call exit_POP(sigAbort,  &
                'model year out of range of CFC11 years')
    end if

    !---------------------------------------------------------------------------
    !   use linear interpolation of time series to current time
    !---------------------------------------------------------------------------

    del_year = seconds_this_year/seconds_in_year
    j = iyear - cfc11_first_year + 1

    if (iyear == cfc11_first_year .and. del_year <= p5) then
       pcfc11_north = cfc11_time_series(1,2)
       pcfc11_south = cfc11_time_series(1,3)
    elseif (iyear == cfc11_last_year .and. del_year >= p5) then
       pcfc11_north = cfc11_time_series(num_cfc11_years,2)
       pcfc11_south = cfc11_time_series(num_cfc11_years,3)
    else
       if (del_year >= p5) then
          pcfc11_north =  &
             (1.5_r8 - del_year)*cfc11_time_series(j    ,2) +  &
             (del_year - p5    )*cfc11_time_series(j + 1,2)
          pcfc11_south =  &
             (1.5_r8 - del_year)*cfc11_time_series(j    ,3) +  &
             (del_year - p5    )*cfc11_time_series(j + 1,3)
       else
          pcfc11_north =  &
             (p5 + del_year)*cfc11_time_series(j    ,2) +  &
             (p5 - del_year)*cfc11_time_series(j - 1,2)
          pcfc11_south =  &
             (p5 + del_year)*cfc11_time_series(j    ,3) +  &
             (p5 - del_year)*cfc11_time_series(j - 1,3)
       endif
    endif

!maltrud debug
!   if (my_task == master_task) write(stdout,*)  &
!      tday, pcfc11_north, pcfc11_south, 'pCFC11_NS'

    pCFC11 = CFC11_INTERP_WEIGHT*pcfc11_north +   &
             (c1 - CFC11_INTERP_WEIGHT)*pcfc11_south

    !---------------------------------------------------------------------------
    !   Compute XKW_ICE. XKW is zero over land, so XKW_ICE is too.
    !---------------------------------------------------------------------------

!untested    !$OMP PARALLEL DO PRIVATE(iblock,XKW_ICE,SCHMIDT_USED,CFC11SAT_1atm, &
!untested    !$OMP                     CFC11SAT_USED,PV,FLUX)

    do iblock = 1, nblocks_clinic

       XKW_ICE = XKW(:,:,iblock)
       where (FICE_USED(:,:,iblock) > 0.2_r8   &
        .and. FICE_USED(:,:,iblock) < 0.9999_r8)
          XKW_ICE = (c1 - FICE_USED(:,:,iblock)) * XKW_ICE
       end where
       where (FICE_USED(:,:,iblock) >= 0.9999_r8)
          XKW_ICE = c0
       end where

    !---------------------------------------------------------------------------
    !   compute CFC11 flux
    !---------------------------------------------------------------------------

       SCHMIDT_USED = SCHMIDT_CFC(SST(:,:,iblock), LAND_MASK(:,:,iblock), 11)

       CFC11SAT_1atm = pCFC11(:,:,iblock)*   &
          SOLUBILITY_CFC(SST(:,:,iblock),  &
                         SSS(:,:,iblock),   &
                         LAND_MASK(:,:,iblock), 11)

       where (LAND_MASK(:,:,iblock))
          PV = XKW_ICE * SQRT(660.0_r8 / SCHMIDT_USED)
          CFC11SAT_USED = AP_USED(:,:,iblock) * CFC11SAT_1atm
!         CFC11_loc = max(c0, SURF_VALS(:,:,cfc11_ind,iblock))
!         FLUX = PV * (CFC11SAT_USED - CFC11_loc)
          FLUX = PV * (CFC11SAT_USED - SURF_VALS(:,:,cfc11_ind,iblock))
          STF_MODULE(:,:,cfc11_ind,iblock) = FLUX
       end where

       CFC11_SFLUX_TAVG(:,:,1,iblock) = FICE_USED(:,:,iblock)
       CFC11_SFLUX_TAVG(:,:,2,iblock) = XKW_ICE(:,:)
       CFC11_SFLUX_TAVG(:,:,3,iblock) = AP_USED(:,:,iblock)
       CFC11_SFLUX_TAVG(:,:,4,iblock) = pCFC11(:,:,iblock)
       CFC11_SFLUX_TAVG(:,:,5,iblock) = PV * mpercm
       CFC11_SFLUX_TAVG(:,:,6,iblock) = SCHMIDT_USED
       CFC11_SFLUX_TAVG(:,:,7,iblock) = CFC11SAT_USED
!      CFC11_SFLUX_TAVG(:,:,7,iblock) = CFC11_INTERP_WEIGHT(:,:,iblock)

    enddo

!untested    !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine cfc11_set_sflux

!***********************************************************************
!BOP
! !IROUTINE: cfc11_tavg_forcing
! !INTERFACE:

 subroutine cfc11_tavg_forcing

! !DESCRIPTION:
!  Make accumulation calls for forcing related tavg fields. This is
!  necessary because the forcing routines are called before tavg flags
!  are set.

! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock              ! block loop index

!-----------------------------------------------------------------------

!untested   !$OMP PARALLEL DO PRIVATE(iblock)

   do iblock = 1,nblocks_clinic
      call accumulate_tavg_field(CFC11_SFLUX_TAVG(:,:,1,iblock),tavg_FICE,iblock,1)
      call accumulate_tavg_field(CFC11_SFLUX_TAVG(:,:,2,iblock),tavg_XKW,iblock,1)
      call accumulate_tavg_field(CFC11_SFLUX_TAVG(:,:,3,iblock),tavg_ATM_PRESS,iblock,1)
      call accumulate_tavg_field(CFC11_SFLUX_TAVG(:,:,4,iblock),tavg_pCFC11,iblock,1)
      call accumulate_tavg_field(CFC11_SFLUX_TAVG(:,:,5,iblock),tavg_PV,iblock,1)
      call accumulate_tavg_field(CFC11_SFLUX_TAVG(:,:,6,iblock),tavg_SCHMIDT_CFC11,iblock,1)
      call accumulate_tavg_field(CFC11_SFLUX_TAVG(:,:,7,iblock),tavg_CFC11SAT,iblock,1)
   end do

!untested   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine cfc11_tavg_forcing

!***********************************************************************
!BOP
! !IROUTINE: SCHMIDT_CFC
! !INTERFACE:

 function SCHMIDT_CFC(PT_2D,LAND_MASK,kn)

! !DESCRIPTION:
!  CFC 11 and 12 schmidt number
!  as a fonction of temperature.
!
!  ref: Zheng et al (1998), JGR, vol 103,No C1
!
!  Original author: J-C Dutay - LSCE
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block) :: &
      PT_2D              ! temperature (degree Celsius)

   logical (log_kind), dimension(nx_block,ny_block) :: LAND_MASK

   integer(int_kind) ::  &
      kn                 ! 11 = CFC-11, 12 = CFC-12

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block) :: &
      SCHMIDT_CFC        ! returned value, non-dimensional

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), dimension(11:12) :: a1, a2, a3, a4

!-----------------------------------------------------------------------
!   coefficients with t in degree Celsius
!-----------------------------------------------------------------------

   a1(11) = 3501.8_r8
   a2(11) = -210.31_r8
   a3(11) =    6.1851_r8
   a4(11) =   -0.07513_r8

   a1(12) = 3845.4_r8
   a2(12) = -228.95_r8
   a3(12) =    6.1908_r8
   a4(12) =   -0.067430_r8

   where (LAND_MASK)
      SCHMIDT_CFC = a1(kn) + a2(kn) * PT_2D + a3(kn) *PT_2D*PT_2D  &
                     + a4(kn) *PT_2D*PT_2D*PT_2D
   elsewhere
      SCHMIDT_CFC = c0
   endwhere

!-----------------------------------------------------------------------
!EOC

 end function SCHMIDT_CFC

!***********************************************************************
!BOP
! !IROUTINE: SOLUBILITY_CFC
! !INTERFACE:

 function SOLUBILITY_CFC(PT,PS,LAND_MASK,kn)

! !DESCRIPTION:
!  CFC 11 and 12 Solubilities in seawater
!  ref: Warner & Weiss (1985) , Deep Sea Research, vol32
!
!  Original author: J-C Dutay - LSCE
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use constants, only: T0_Kelvin, c1000

! !INPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block) :: &
      PT,              & ! temperature (degree Celsius)
      PS                 ! salinity    (o/oo) (multiply by 1000 to get psu)

   logical (log_kind), dimension(nx_block,ny_block) :: LAND_MASK

   integer(int_kind) ::  &
      kn                 ! 11 = CFC-11, 12 = CFC-12

! !OUTPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block) :: &
      SOLUBILITY_CFC     ! returned value, in mol/m3/pptv
                         ! 1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real(r8), dimension(nx_block,ny_block) :: WORK

   real(r8), dimension(11:12) ::  &
      a1, a2, a3, a4, b1, b2, b3

!-----------------------------------------------------------------------
!  for CFC 11
!-----------------------------------------------------------------------

   a1 ( 11) = -229.9261_r8
   a2 ( 11) =  319.6552_r8
   a3 ( 11) =  119.4471_r8
   a4 ( 11) =   -1.39165_r8
   b1 ( 11) =   -0.142382_r8
   b2 ( 11) =    0.091459_r8
   b3 ( 11) =   -0.0157274_r8

!-----------------------------------------------------------------------
!  for CFC/12
!-----------------------------------------------------------------------

   a1 ( 12) = -218.0971_r8
   a2 ( 12) =  298.9702_r8
   a3 ( 12) =  113.8049_r8
   a4 ( 12) =   -1.39165_r8
   b1 ( 12) =   -0.143566_r8
   b2 ( 12) =    0.091015_r8
   b3 ( 12) =   -0.0153924_r8

   WORK = merge( ((PT + T0_Kelvin)* 0.01_r8), c1, LAND_MASK)

!-----------------------------------------------------------------------
!  coefficient for solubility in  mol/l/atm
!-----------------------------------------------------------------------

   where (LAND_MASK)
      SOLUBILITY_CFC  &
         = exp ( a1 ( kn)                &
               + a2 ( kn)/ WORK          &
               + a3 ( kn)* log ( WORK )  &
               + a4 ( kn)* WORK * WORK   &
               + PS*( ( b3( kn)*WORK + b2( kn) )*WORK + b1(kn) ) )
   elsewhere
      SOLUBILITY_CFC = c0
   endwhere

!-----------------------------------------------------------------------
!  conversion from mol/(l * atm) to mol/(m^3 * atm) to mol/(m3 * pptv)
!-----------------------------------------------------------------------

   SOLUBILITY_CFC = c1000 * SOLUBILITY_CFC
   SOLUBILITY_CFC = 1.0e-12_r8 * SOLUBILITY_CFC

!-----------------------------------------------------------------------
!EOC

 end function SOLUBILITY_CFC

!***********************************************************************

end module cfc11_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
