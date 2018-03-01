!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module cfc_mod

!BOP
! !MODULE: cfc_mod
!
!  Module for Chlorofluorocarbons (CFCs)
!
!  The units of concentration for these tracers are
!     fmol/cm^3 == nmol/m^3 == pmol/l ~= pmol/kg.
!  These units are chosen because ship measurements are typically
!  given in units of pmol/kg.
!
!  The units of surface fluxes for these tracers are
!     fmol/cm^3 * cm/s == fmol/cm^2/s.
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!  SVN:$Id: $

! !USES:

   use POP_KindsMod
   use POP_ErrorMod

   use kinds_mod
   use blocks, only: nx_block, ny_block, block
   use domain_size, only: max_blocks_clinic, km
   use domain, only: nblocks_clinic, distrb_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use constants, only: c0, c1
   use io_types, only: stdout
   use io_tools, only: document
   use tavg, only: define_tavg_field, accumulate_tavg_field
   use passive_tracer_tools, only: forcing_monthly_every_ts, &
       ind_name_pair, tracer_read, read_field
   use broadcast
   use netcdf

   implicit none
   save

!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: &
       cfc_tracer_cnt, &
       cfc_init, &
       cfc_set_sflux,  &
       cfc_tavg_forcing

!EOP
!BOC

!-----------------------------------------------------------------------
!  module variables required by passive_tracers
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
       cfc_tracer_cnt = 2

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      cfc11_ind =  1,  & ! CFC11
      cfc12_ind =  2     ! CFC12

!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(cfc_tracer_cnt) :: &
      ind_name_table = (/ &
      ind_name_pair(cfc11_ind, 'CFC11'), &
      ind_name_pair(cfc12_ind, 'CFC12') /)

!-----------------------------------------------------------------------
!  mask that eases avoidance of computation over land
!-----------------------------------------------------------------------

   logical (log_kind), dimension(:,:,:), allocatable :: &
      LAND_MASK

!-----------------------------------------------------------------------
!  forcing related variables
!-----------------------------------------------------------------------

   character(char_len) :: &
      cfc_formulation,     & ! how to calculate flux (ocmip or model)
      pcfc_file              ! filename for ascii time series of atm cfc11

   integer (int_kind) ::  &
      model_year,          & ! arbitrary model year
      data_year,           & ! year in data that corresponds to model_year
      pcfc_data_len          ! length of atmospheric pcfc record

   real (r8), parameter :: &
      max_pcfc_extension = 2.0_r8
      ! maximum number of years that pcfc record will be extrapolated

   real (r8), dimension(:), allocatable :: &
      pcfc_date,           & ! date for atmospheric pcfc record (years)
      pcfc11_nh,           & ! pcfc11 data for northern hemisphere (pmol/mol)
      pcfc11_sh,           & ! pcfc11 data for southern hemisphere (pmol/mol)
      pcfc12_nh,           & ! pcfc12 data for northern hemisphere (pmol/mol)
      pcfc12_sh              ! pcfc12 data for southern hemisphere (pmol/mol)

   real (r8), dimension(:,:,:,:), allocatable :: &
      INTERP_WORK            ! temp array for interpolate_forcing output

   type(forcing_monthly_every_ts) :: &
      fice_file,           & ! ice fraction, if read from file
      xkw_file,            & ! a * wind-speed ** 2, if read from file
      ap_file                ! atmoshperic pressure, if read from file

!-----------------------------------------------------------------------
!  define tavg id for 2d fields related to surface fluxes
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable ::   &
      CFC_SFLUX_TAVG

   integer (int_kind) :: &
      tavg_CFC_IFRAC,      & ! tavg id for ice fraction
      tavg_CFC_XKW,        & ! tavg id for xkw
      tavg_CFC_ATM_PRESS,  & ! tavg id for atmospheric pressure
      tavg_pCFC11,         & ! tavg id for cfc11 partial pressure
      tavg_pCFC12,         & ! tavg id for cfc12 partial pressure
      tavg_CFC11_SCHMIDT,  & ! tavg id for cfc11 Schmidt number
      tavg_CFC12_SCHMIDT,  & ! tavg id for cfc12 Schmidt number
      tavg_CFC11_PV,       & ! tavg id for cfc11 piston velocity
      tavg_CFC11_surf_sat, & ! tavg id for cfc11 surface saturation
      tavg_CFC12_PV,       & ! tavg id for cfc12 piston velocity
      tavg_CFC12_surf_sat    ! tavg id for cfc12 surface saturation

!-----------------------------------------------------------------------
!  data_ind is the index into data for current timestep, i.e
!  data_ind is largest integer less than pcfc_data_len s.t.
!  pcfc_date(i) <= iyear + (iday_of_year-1+frac_day)/days_in_year
!                  - model_year + data_year
!  Note that data_ind is always strictly less than pcfc_data_len.
!  To enable OpenMP parallelism, duplicating data_ind for each block
!-----------------------------------------------------------------------

   integer (int_kind), dimension(:), allocatable :: &
      data_ind

!-----------------------------------------------------------------------
!  timers
!-----------------------------------------------------------------------

   integer (int_kind) :: cfc_sflux_timer

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: cfc_init
! !INTERFACE:

 subroutine cfc_init(init_ts_file_fmt, read_restart_filename, &
                     tracer_d_module, TRACER_MODULE, errorCode)

! !DESCRIPTION:
!  Initialize cfc tracer module. This involves setting metadata, reading
!  the modules namelist and setting initial conditions.

! !REVISION HISTORY:
!  same as module

! !USES:

   use constants, only: char_blank, delim_fmt
   use prognostic, only: curtime, oldtime
   use grid, only: KMT, n_topo_smooth, fill_points
   use grid, only: REGION_MASK
   use io_types, only: nml_in, nml_filename
   use prognostic, only: tracer_field
   use timers, only: get_timer
   use passive_tracer_tools, only: init_forcing_monthly_every_ts, &
       rest_read_tracer_block, file_read_tracer_block

! !INPUT PARAMETERS:

   character (*), intent(in) ::  &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(cfc_tracer_cnt), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

   real (r8), dimension(nx_block,ny_block,km,cfc_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: sub_name = 'cfc_mod:cfc_init'

   character(char_len) :: &
      init_cfc_option,        & ! option for initialization of bgc
      init_cfc_init_file,     & ! filename for option 'file'
      init_cfc_init_file_fmt    ! file format for option 'file'

   integer (int_kind) :: &
      n,                      & ! index for looping over tracers
      k,                      & ! index for looping over depth levels
      iblock,                 & ! index for looping over blocks
      nml_error                 ! namelist i/o error flag

   type(tracer_read), dimension(cfc_tracer_cnt) :: &
      tracer_init_ext           ! namelist variable for initializing tracers

   type(tracer_read) :: &
      gas_flux_fice,          & ! ice fraction for gas fluxes
      gas_flux_ws,            & ! wind speed for gas fluxes
      gas_flux_ap               ! atmospheric pressure for gas fluxes

   namelist /cfc_nml/ &
      init_cfc_option, init_cfc_init_file, init_cfc_init_file_fmt, &
      tracer_init_ext, pcfc_file, model_year, data_year, &
      cfc_formulation, gas_flux_fice, gas_flux_ws, gas_flux_ap

   character (char_len) ::  &
      cfc_restart_filename      ! modified file name for restart file

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

   do n = 1, cfc_tracer_cnt
      tracer_d_module(n)%short_name = ind_name_table(n)%name
      tracer_d_module(n)%long_name  = ind_name_table(n)%name
      tracer_d_module(n)%units      = 'fmol/cm^3'
      tracer_d_module(n)%tend_units = 'fmol/cm^3/s'
      tracer_d_module(n)%flux_units = 'fmol/cm^2/s'
   end do

!-----------------------------------------------------------------------
!  default namelist settings
!-----------------------------------------------------------------------

   init_cfc_option        = 'unknown'
   init_cfc_init_file     = 'unknown'
   init_cfc_init_file_fmt = 'bin'

   do n = 1, cfc_tracer_cnt
      tracer_init_ext(n)%mod_varname  = 'unknown'
      tracer_init_ext(n)%filename     = 'unknown'
      tracer_init_ext(n)%file_varname = 'unknown'
      tracer_init_ext(n)%scale_factor = c1
      tracer_init_ext(n)%default_val  = c0
      tracer_init_ext(n)%file_fmt     = 'bin'
   end do

   pcfc_file       = 'unknown'
   model_year      = 1
   data_year       = 1931
   cfc_formulation = 'model'

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
         read(nml_in, nml=cfc_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(sub_name, 'cfc_nml not found')
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ sub_name)
   endif

!-----------------------------------------------------------------------
!  broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_cfc_option, master_task)
   call broadcast_scalar(init_cfc_init_file, master_task)
   call broadcast_scalar(init_cfc_init_file_fmt, master_task)

   do n = 1, cfc_tracer_cnt
      call broadcast_scalar(tracer_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%filename, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(tracer_init_ext(n)%default_val, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_fmt, master_task)
   end do

   call broadcast_scalar(pcfc_file, master_task)
   call broadcast_scalar(model_year, master_task)
   call broadcast_scalar(data_year, master_task)
   call broadcast_scalar(cfc_formulation, master_task)

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

!-----------------------------------------------------------------------
!   initialize tracers
!-----------------------------------------------------------------------

   select case (init_cfc_option)

   case ('ccsm_startup', 'zero', 'ccsm_startup_spunup')
      TRACER_MODULE = c0
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d CFCs set to all zeros'
          write(stdout,delim_fmt)
      endif

   case ('restart', 'ccsm_continue', 'ccsm_branch', 'ccsm_hybrid' )

      cfc_restart_filename = char_blank

      if (init_cfc_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(sub_name, 'no restart file to read CFCs from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ sub_name)
         endif
         cfc_restart_filename = read_restart_filename
         init_cfc_init_file_fmt = init_ts_file_fmt

      else  ! do not read from TS restart file

         cfc_restart_filename = trim(init_cfc_init_file)

      endif

      call rest_read_tracer_block(init_cfc_init_file_fmt, &
                                  cfc_restart_filename,   &
                                  tracer_d_module,        &
                                  TRACER_MODULE)

   case ('file')

      call document(sub_name, 'CFCs being read from separate file')

      call file_read_tracer_block(init_cfc_init_file_fmt, &
                                  init_cfc_init_file,     &
                                  tracer_d_module,        &
                                  ind_name_table,         &
                                  tracer_init_ext,        &
                                  TRACER_MODULE)

      if (n_topo_smooth > 0) then
         do n = 1, cfc_tracer_cnt
            do k = 1, km
               call fill_points(k,TRACER_MODULE(:,:,k,n,curtime,:), &
                                errorCode)

               if (errorCode /= POP_Success) then
                  call POP_ErrorSet(errorCode, &
                     'cfc_init: error in fill_points')
                  return
               endif
            end do
         end do
      endif

   case default
      call document(sub_name, 'init_cfc_option', init_cfc_option)
      call exit_POP(sigAbort, 'unknown init_cfc_option')

   end select

!-----------------------------------------------------------------------
!  apply land mask to tracers
!-----------------------------------------------------------------------

   do iblock = 1, nblocks_clinic
   do n = 1, cfc_tracer_cnt
      do k = 1, km
         where (k > KMT(:,:,iblock))
            TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
            TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
         end where
      end do
   end do
   end do

!-----------------------------------------------------------------------
!  allocate and initialize LAND_MASK (true for ocean points)
!-----------------------------------------------------------------------

   allocate( LAND_MASK(nx_block,ny_block,max_blocks_clinic) )
   LAND_MASK = (KMT.gt.0)

   call get_timer(cfc_sflux_timer, 'CFC_SFLUX', 1, distrb_clinic%nprocs)

!-----------------------------------------------------------------------
!  call other initialization subroutines
!-----------------------------------------------------------------------

   call cfc_init_tavg
   call cfc_init_sflux

!-----------------------------------------------------------------------
!EOC

 end subroutine cfc_init

!***********************************************************************
!BOP
! !IROUTINE: cfc_init_tavg
! !INTERFACE:

 subroutine cfc_init_tavg

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

   call define_tavg_field(tavg_CFC_IFRAC,'CFC_IFRAC',2,           &
                          long_name='Ice Fraction for CFC fluxes',&
                          units='fraction', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC_XKW,'CFC_XKW',2,               &
                          long_name='XKW for CFC fluxes',         &
                          units='cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC_ATM_PRESS,'CFC_ATM_PRESS',2,   &
                          long_name='Atmospheric Pressure for CFC fluxes',&
                          units='atmospheres', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_pCFC11,'pCFC11',2,                 &
                          long_name='CFC11 atmospheric partial pressure',&
                          units='pmol/mol', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_pCFC12,'pCFC12',2,                 &
                          long_name='CFC12 atmospheric partial pressure',&
                          units='pmol/mol', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC11_SCHMIDT,'CFC11_SCHMIDT',2,   &
                          long_name='CFC11 Schmidt Number',       &
                          units='none', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC12_SCHMIDT,'CFC12_SCHMIDT',2,   &
                          long_name='CFC12 Schmidt Number',       &
                          units='none', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC11_PV,'CFC11_PV',2,             &
                          long_name='CFC11 piston velocity',      &
                          units='cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC11_surf_sat,'CFC11_surf_sat',2, &
                          long_name='CFC11 Saturation',           &
                          units='fmol/cm^3', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC12_PV,'CFC12_PV',2,             &
                          long_name='CFC12 piston velocity',      &
                          units='cm/s', grid_loc='2110')
   var_cnt = var_cnt+1

   call define_tavg_field(tavg_CFC12_surf_sat,'CFC12_surf_sat',2, &
                          long_name='CFC12 Saturation',           &
                          units='fmol/cm^3', grid_loc='2110')
   var_cnt = var_cnt+1

!-----------------------------------------------------------------------

   allocate(CFC_SFLUX_TAVG(nx_block,ny_block,var_cnt,max_blocks_clinic))
   CFC_SFLUX_TAVG = c0

!-----------------------------------------------------------------------
!EOC

 end subroutine cfc_init_tavg

!***********************************************************************
!BOP
! !IROUTINE: cfc_init_sflux
! !INTERFACE:

 subroutine cfc_init_sflux

! !USES:

   use forcing_tools, only: find_forcing_times

! !DESCRIPTION:
!  Initialize surface flux computations for cfc tracer module.
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: sub_name = 'cfc_mod:cfc_init_sflux'

   integer (int_kind) :: &
      n,                 & ! index for looping over tracers
      iblock               ! index for looping over blocks

   real (r8), dimension (nx_block,ny_block) :: WORK

   real (r8), dimension (nx_block,ny_block,12,max_blocks_clinic), target :: &
      WORK_READ            ! temporary space to read in fields

!-----------------------------------------------------------------------

   call read_pcfc_data

!-----------------------------------------------------------------------
!  read gas flux forcing (if required)
!  otherwise, use values passed in
!-----------------------------------------------------------------------

   select case (cfc_formulation)

   case ('ocmip')

!-----------------------------------------------------------------------
!  allocate space for interpolate_forcing
!-----------------------------------------------------------------------

      allocate(INTERP_WORK(nx_block,ny_block,max_blocks_clinic,1))

!-----------------------------------------------------------------------
!  first, read ice file
!-----------------------------------------------------------------------

      allocate(fice_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(fice_file%input%file_fmt, &
                      fice_file%input%filename, &
                      fice_file%input%file_varname, &
                      WORK_READ)
      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         fice_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            fice_file%DATA(:,:,iblock,1,n) = c0
         fice_file%DATA(:,:,iblock,1,n) = &
            fice_file%DATA(:,:,iblock,1,n) * fice_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(fice_file%data_time, &
                              fice_file%data_inc, fice_file%interp_type, &
                              fice_file%data_next, fice_file%data_time_min_loc, &
                              fice_file%data_update, fice_file%data_type)

!-----------------------------------------------------------------------
!  next, read piston velocity file
!-----------------------------------------------------------------------

      allocate(xkw_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(xkw_file%input%file_fmt, &
                      xkw_file%input%filename, &
                      xkw_file%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         xkw_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            xkw_file%DATA(:,:,iblock,1,n) = c0
         xkw_file%DATA(:,:,iblock,1,n) = &
            xkw_file%DATA(:,:,iblock,1,n) * xkw_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(xkw_file%data_time, &
                              xkw_file%data_inc, xkw_file%interp_type, &
                              xkw_file%data_next, xkw_file%data_time_min_loc, &
                              xkw_file%data_update, xkw_file%data_type)

!-----------------------------------------------------------------------
!  last, read atmospheric pressure file
!-----------------------------------------------------------------------

      allocate(ap_file%DATA(nx_block,ny_block,max_blocks_clinic,1,12))

      call read_field(ap_file%input%file_fmt, &
                      ap_file%input%filename, &
                      ap_file%input%file_varname, &
                      WORK_READ)

      !$OMP PARALLEL DO PRIVATE(iblock, n)
      do iblock=1,nblocks_clinic
      do n=1,12
         ap_file%DATA(:,:,iblock,1,n) = WORK_READ(:,:,n,iblock)
         where (.not. LAND_MASK(:,:,iblock)) &
            ap_file%DATA(:,:,iblock,1,n) = c0
         ap_file%DATA(:,:,iblock,1,n) = &
            ap_file%DATA(:,:,iblock,1,n) * ap_file%input%scale_factor
      end do
      end do
      !$OMP END PARALLEL DO

      call find_forcing_times(ap_file%data_time, &
                              ap_file%data_inc, ap_file%interp_type, &
                              ap_file%data_next, ap_file%data_time_min_loc, &
                              ap_file%data_update, ap_file%data_type)

   case ('model')

      if (my_task == master_task) then
         write(stdout,*)  &
            ' Using fields from model forcing for calculating CFC flux'
      endif

   case default
      call document(sub_name, 'cfc_formulation', cfc_formulation)

      call exit_POP(sigAbort, &
                    'cfc_init_sflux: Unknown value for cfc_formulation')

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine cfc_init_sflux

!***********************************************************************
!BOP
! !IROUTINE: read_pcfc_data
! !INTERFACE:

 subroutine read_pcfc_data

! !DESCRIPTION:
!  subroutine to read in atmospheric pcfc data
!
!  Have the master_task do the following :
!     1) get length of data
!     2) allocate memory for data
!     3) read in data, checking for consistent lengths
!  Then, outside master_task conditional
!     1) broadcast length of data
!     2) have non-mastertasks allocate memory for data
!     3) broadcast data

! !USES:

!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   character(*), parameter :: sub_name = 'cfc_mod:read_pcfc_data'

   character (len=char_len) :: &
      varname           ! name of variable being processed

   integer (int_kind) :: &
      stat,           & ! status of netCDF call
      ncid,           & ! netCDF file id
      varid,          & ! netCDF variable id
      ndims             ! number of dimensions for varid

   integer (int_kind), dimension(1) :: &
      data_dimid        ! netCDF dimension id that all data should have

!-----------------------------------------------------------------------
!  perform netCDF I/O on master_task
!  jump out of master_task conditional if an error is encountered
!-----------------------------------------------------------------------

   if (my_task == master_task) then

      stat = nf90_open(pcfc_file, 0, ncid)
      if (stat /= 0) then
         write(stdout,*) 'error from nf_open: ', nf90_strerror(stat)
         go to 99
      endif

!-----------------------------------------------------------------------
!  get length of data by examining pcfc11_nh
!  keep track of dimid for later comparison when reading in data
!-----------------------------------------------------------------------

      varname = 'pcfc11_nh'

      stat = nf90_inq_varid(ncid, varname, varid)

      if (stat /= 0) then
         write(stdout,*) 'error from nf_inq_varid for pcfc11_nh: ', nf90_strerror(stat)
         go to 99
      endif

      stat = nf90_inquire_variable(ncid, varid, ndims=ndims)
      if (stat /= 0) then
         write(stdout,*) 'nf_inq_varndims for pcfc11_nh: ', nf90_strerror(stat)
         go to 99
      endif
      if (ndims /= 1) then
         write(stdout,*) 'ndims /= 1 for pcfc11_nh'
         go to 99
      endif

      stat = nf90_inquire_variable(ncid, varid, dimids=data_dimid)
      if (stat /= 0) then
         write(stdout,*) 'nf_inq_vardimid for pcfc11_nh: ', nf90_strerror(stat)
         go to 99
      endif

      stat = nf90_inquire_dimension(ncid, data_dimid(1), len=pcfc_data_len)
      if (stat /= 0) then
         write(stdout,*) 'nf_inq_dimlen for pcfc11_nh: ', nf90_strerror(stat)
         go to 99
      endif

      call document(sub_name, 'pcfc_data_len', pcfc_data_len)

      allocate(pcfc_date(pcfc_data_len))
      allocate(pcfc11_nh(pcfc_data_len))
      allocate(pcfc11_sh(pcfc_data_len))
      allocate(pcfc12_nh(pcfc_data_len))
      allocate(pcfc12_sh(pcfc_data_len))

      stat = nf90_inquire_dimension(ncid, data_dimid(1), name=varname)
      if (stat /= 0) then
         write(stdout,*) 'nf_inq_dimname for dim of pcfc11_nh: ', nf90_strerror(stat)
         go to 99
      endif

      call read_1dvar_cdf(ncid, data_dimid, varname,     pcfc_date, stat)
      if (stat /= 0) go to 99
      call read_1dvar_cdf(ncid, data_dimid, 'pcfc11_nh', pcfc11_nh, stat)
      if (stat /= 0) go to 99
      call read_1dvar_cdf(ncid, data_dimid, 'pcfc11_sh', pcfc11_sh, stat)
      if (stat /= 0) go to 99
      call read_1dvar_cdf(ncid, data_dimid, 'pcfc12_nh', pcfc12_nh, stat)
      if (stat /= 0) go to 99
      call read_1dvar_cdf(ncid, data_dimid, 'pcfc12_sh', pcfc12_sh, stat)
      if (stat /= 0) go to 99

      stat = nf90_close(ncid)
      if (stat /= 0) then
         write(stdout,*) 'nf_close: ', nf90_strerror(stat)
         go to 99
      endif

      call document(sub_name, 'pcfc_data_len', pcfc_data_len)
      call document(sub_name, 'pcfc_date(end)', pcfc_date(pcfc_data_len))
      call document(sub_name, 'pcfc11_nh(end)', pcfc11_nh(pcfc_data_len))
      call document(sub_name, 'pcfc11_sh(end)', pcfc11_sh(pcfc_data_len))
      call document(sub_name, 'pcfc12_nh(end)', pcfc12_nh(pcfc_data_len))
      call document(sub_name, 'pcfc12_sh(end)', pcfc12_sh(pcfc_data_len))

   endif ! my_task == master_task

99 call broadcast_scalar(stat, master_task)
   if (stat /= 0) call exit_POP(sigAbort, 'stopping in ' /&
                                                          &/ sub_name)

   call broadcast_scalar(pcfc_data_len, master_task)

   if (my_task /= master_task) then
      allocate(pcfc_date(pcfc_data_len))
      allocate(pcfc11_nh(pcfc_data_len))
      allocate(pcfc11_sh(pcfc_data_len))
      allocate(pcfc12_nh(pcfc_data_len))
      allocate(pcfc12_sh(pcfc_data_len))
   endif

   call broadcast_array(pcfc_date, master_task)
   call broadcast_array(pcfc11_nh, master_task)
   call broadcast_array(pcfc11_sh, master_task)
   call broadcast_array(pcfc12_nh, master_task)
   call broadcast_array(pcfc12_sh, master_task)

!-----------------------------------------------------------------------
!EOC

 end subroutine read_pcfc_data

!***********************************************************************
!BOP
! !IROUTINE: read_1dvar_cdf
! !INTERFACE:

 subroutine read_1dvar_cdf(ncid, data_dimid, varname, data, stat)

! !DESCRIPTION:
!  Subroutine to read in a single 1D variable from a netCDF file
!  that is supposed to be on a particular dimension
!
! !REVISION HISTORY:
!  same as module

! !USES:

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      ncid              ! netCDF file id

   integer (int_kind), dimension(1), intent(in) :: &
      data_dimid        ! netCDF dimension id that all data should have

   character (len=*), intent(in) :: &
      varname           ! name of variable being read

! !OUTPUT PARAMETERS:

   real (r8), dimension(:), intent(out) :: &
      data              ! where data is going

   integer (int_kind), intent(out) :: &
      stat              ! status of netCDF call
!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      varid,          & ! netCDF variable id
      ndims             ! number of dimensions for varid

   integer (int_kind), dimension(1) :: &
      dimid             ! netCDF dimension id

!-----------------------------------------------------------------------

   stat = nf90_inq_varid(ncid, varname, varid)
   if (stat /= 0) then
      write(stdout,*) 'nf_inq_varid for ', trim(varname), ' : ', nf90_strerror(stat)
      return
   endif

   stat = nf90_inquire_variable(ncid, varid, ndims=ndims)
   if (stat /= 0) then
      write(stdout,*) 'nf_inq_varndims for ', trim(varname), ' : ', nf90_strerror(stat)
      return
   endif
   if (ndims /= 1) then
      write(stdout,*) 'ndims /= 1 for ', trim(varname)
      return
   endif

   stat = nf90_inquire_variable(ncid, varid, dimids=dimid)
   if (stat /= 0) then
      write(stdout,*) 'nf_inq_vardimid for ', trim(varname), ' : ', nf90_strerror(stat)
      return
   endif
   if (dimid(1) /= data_dimid(1)) then
      write(stdout,*) 'dimid mismatch for ', trim(varname)
      return
   endif

   stat = nf90_get_var(ncid, varid, data)
   if (stat /= 0) then
      write(stdout,*) 'nf_get_var_double for ', trim(varname), ' : ', nf90_strerror(stat)
      return
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine read_1dvar_cdf

!***********************************************************************
!BOP
! !IROUTINE: cfc_set_sflux
! !INTERFACE:

 subroutine cfc_set_sflux(U10_SQR,IFRAC,PRESS,SST,SSS, &
                          SURF_VALS_OLD,SURF_VALS_CUR,STF_MODULE)

! !DESCRIPTION:
!  Compute CFC11 surface flux and store related tavg fields for
!  subsequent accumulating.

! !REVISION HISTORY:
!  same as module

! !USES:

   use constants, only: field_loc_center, field_type_scalar, p5
   use time_management, only: thour00
   use forcing_tools, only: update_forcing_data, interpolate_forcing
   use timers, only: timer_start, timer_stop

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      U10_SQR,   & ! 10m wind speed squared (cm/s)**2
      IFRAC,     & ! sea ice fraction (non-dimensional)
      PRESS,     & ! sea level atmospheric pressure (dyne/cm**2)
      SST,       & ! sea surface temperature (C)
      SSS          ! sea surface salinity (psu)

   real (r8), dimension(nx_block,ny_block,cfc_tracer_cnt,max_blocks_clinic), &
         intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR ! module tracers

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,cfc_tracer_cnt,max_blocks_clinic), &
         intent(inout) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      iblock             ! block index

   integer (int_kind) :: i, j

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      IFRAC_USED,      & ! used ice fraction (non-dimensional)
      XKW_USED,        & ! part of piston velocity (cm/s)
      AP_USED            ! used atm pressure (converted from dyne/cm**2 to atm)

   real (r8), dimension(nx_block,ny_block) :: &
      SURF_VALS,       & ! filtered surface tracer values
      pCFC11,          & ! atmospheric CFC11 mole fraction (pmol/mol)
      pCFC12,          & ! atmospheric CFC11 mole fraction (pmol/mol)
      CFC11_SCHMIDT,   & ! CFC11 Schmidt number
      CFC12_SCHMIDT,   & ! CFC12 Schmidt number
      CFC11_SOL_0,     & ! solubility of CFC11 at 1 atm (mol/l/atm)
      CFC12_SOL_0,     & ! solubility of CFC12 at 1 atm (mol/l/atm)
      XKW_ICE,         & ! common portion of piston vel., (1-fice)*xkw (cm/s)
      PV,              & ! piston velocity (cm/s)
      CFC_surf_sat       ! CFC surface saturation (either CFC11 or CFC12) (fmol/cm^3)

   character (char_len) :: &
      tracer_data_label          ! label for what is being updated

   character (char_len), dimension(1) :: &
      tracer_data_names          ! short names for input data fields

   integer (int_kind), dimension(1) :: &
      tracer_bndy_loc,          &! location and field type for ghost
      tracer_bndy_type           !    cell updates

   logical (log_kind), save :: &
      first = .true.

!-----------------------------------------------------------------------
!  local parameters
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      xkw_coeff = 8.6e-9_r8      ! xkw_coeff = 0.31 cm/hr s^2/m^2 in (s/cm)

!-----------------------------------------------------------------------

   call timer_start(cfc_sflux_timer)

   if (first) then
      allocate( data_ind(max_blocks_clinic) )
      data_ind = -1
      first = .false.
   endif

   do iblock = 1, nblocks_clinic
      IFRAC_USED(:,:,iblock) = c0
      XKW_USED(:,:,iblock) = c0
      AP_USED(:,:,iblock) = c0
   end do

!-----------------------------------------------------------------------
!  Interpolate gas flux forcing data if necessary
!-----------------------------------------------------------------------

   if (cfc_formulation == 'ocmip') then
       if (thour00 >= fice_file%data_update) then
          tracer_data_names = fice_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Ice Fraction'
          call update_forcing_data(          fice_file%data_time,   &
               fice_file%data_time_min_loc,  fice_file%interp_type, &
               fice_file%data_next,          fice_file%data_update, &
               fice_file%data_type,          fice_file%data_inc,    &
               fice_file%DATA(:,:,:,:,1:12), fice_file%data_renorm, &
               tracer_data_label,            tracer_data_names,     &
               tracer_bndy_loc,              tracer_bndy_type,      &
               fice_file%filename,           fice_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            fice_file%DATA(:,:,:,:,1:12), &
            fice_file%data_time,         fice_file%interp_type, &
            fice_file%data_time_min_loc, fice_file%interp_freq, &
            fice_file%interp_inc,        fice_file%interp_next, &
            fice_file%interp_last,       0)
       IFRAC_USED = INTERP_WORK(:,:,:,1)

       if (thour00 >= xkw_file%data_update) then
          tracer_data_names = xkw_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Piston Velocity'
          call update_forcing_data(         xkw_file%data_time,   &
               xkw_file%data_time_min_loc,  xkw_file%interp_type, &
               xkw_file%data_next,          xkw_file%data_update, &
               xkw_file%data_type,          xkw_file%data_inc,    &
               xkw_file%DATA(:,:,:,:,1:12), xkw_file%data_renorm, &
               tracer_data_label,           tracer_data_names,    &
               tracer_bndy_loc,             tracer_bndy_type,     &
               xkw_file%filename,           xkw_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            xkw_file%DATA(:,:,:,:,1:12), &
            xkw_file%data_time,         xkw_file%interp_type, &
            xkw_file%data_time_min_loc, xkw_file%interp_freq, &
            xkw_file%interp_inc,        xkw_file%interp_next, &
            xkw_file%interp_last,       0)
       XKW_USED = INTERP_WORK(:,:,:,1)

       if (thour00 >= ap_file%data_update) then
          tracer_data_names = ap_file%input%file_varname
          tracer_bndy_loc   = field_loc_center
          tracer_bndy_type  = field_type_scalar
          tracer_data_label = 'Atmospheric Pressure'
          call update_forcing_data(        ap_file%data_time,   &
               ap_file%data_time_min_loc,  ap_file%interp_type, &
               ap_file%data_next,          ap_file%data_update, &
               ap_file%data_type,          ap_file%data_inc,    &
               ap_file%DATA(:,:,:,:,1:12), ap_file%data_renorm, &
               tracer_data_label,          tracer_data_names,   &
               tracer_bndy_loc,            tracer_bndy_type,    &
               ap_file%filename,           ap_file%input%file_fmt)
       endif
       call interpolate_forcing(INTERP_WORK, &
            ap_file%DATA(:,:,:,:,1:12), &
            ap_file%data_time,         ap_file%interp_type, &
            ap_file%data_time_min_loc, ap_file%interp_freq, &
            ap_file%interp_inc,        ap_file%interp_next, &
            ap_file%interp_last,       0)
       AP_USED = INTERP_WORK(:,:,:,1)
   endif

   !$OMP PARALLEL DO PRIVATE(iblock,SURF_VALS,pCFC11,pCFC12,CFC11_SCHMIDT, &
   !$OMP                     CFC12_SCHMIDT,CFC11_SOL_0,CFC12_SOL_0,XKW_ICE,&
   !$OMP                     PV,CFC_surf_sat)
   do iblock = 1, nblocks_clinic

      if (cfc_formulation == 'ocmip') then
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) < 0.2000_r8) &
            IFRAC_USED(:,:,iblock) = 0.2000_r8
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) > 0.9999_r8) &
            IFRAC_USED(:,:,iblock) = 0.9999_r8
      endif

      if (cfc_formulation == 'model') then
         where (LAND_MASK(:,:,iblock))
            IFRAC_USED(:,:,iblock) = IFRAC(:,:,iblock)

            XKW_USED(:,:,iblock) = xkw_coeff * U10_SQR(:,:,iblock)

            AP_USED(:,:,iblock) = PRESS(:,:,iblock)
         endwhere
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) < c0) &
            IFRAC_USED(:,:,iblock) = c0
         where (LAND_MASK(:,:,iblock) .and. IFRAC_USED(:,:,iblock) > c1) &
            IFRAC_USED(:,:,iblock) = c1
      endif

!-----------------------------------------------------------------------
!  assume PRESS is in cgs units (dyne/cm**2) since that is what is
!    required for pressure forcing in barotropic
!  want units to be atmospheres
!  convertion from dyne/cm**2 to Pascals is P(mks) = P(cgs)/10.
!  convertion from Pascals to atm is P(atm) = P(Pa)/101.325e+3_r8
!-----------------------------------------------------------------------

      AP_USED(:,:,iblock) = AP_USED(:,:,iblock) * (c1 / 1013.25e+3_r8)

      call comp_pcfc(iblock, LAND_MASK(:,:,iblock), data_ind(iblock), &
                     pCFC11, pCFC12)

      call comp_cfc_schmidt(LAND_MASK(:,:,iblock), SST(:,:,iblock), &
                            CFC11_SCHMIDT, CFC12_SCHMIDT)

      call comp_cfc_sol_0(LAND_MASK(:,:,iblock), SST(:,:,iblock), SSS(:,:,iblock), &
                          CFC11_SOL_0, CFC12_SOL_0)

      where (LAND_MASK(:,:,iblock))
         CFC_SFLUX_TAVG(:,:,1,iblock) = IFRAC_USED(:,:,iblock)
         CFC_SFLUX_TAVG(:,:,2,iblock) = XKW_USED(:,:,iblock)
         CFC_SFLUX_TAVG(:,:,3,iblock) = AP_USED(:,:,iblock)
         CFC_SFLUX_TAVG(:,:,4,iblock) = pCFC11
         CFC_SFLUX_TAVG(:,:,5,iblock) = pCFC12
         CFC_SFLUX_TAVG(:,:,6,iblock) = CFC11_SCHMIDT
         CFC_SFLUX_TAVG(:,:,7,iblock) = CFC12_SCHMIDT

         XKW_ICE = (c1 - IFRAC_USED(:,:,iblock)) * XKW_USED(:,:,iblock)

         PV = XKW_ICE * sqrt(660.0_r8 / CFC11_SCHMIDT)
         CFC_SFLUX_TAVG(:,:,8,iblock) = PV
         CFC_surf_sat = AP_USED(:,:,iblock) * CFC11_SOL_0 * pCFC11
         CFC_SFLUX_TAVG(:,:,9,iblock) = CFC_surf_sat
         SURF_VALS = p5*(SURF_VALS_OLD(:,:,cfc11_ind,iblock) + &
                         SURF_VALS_CUR(:,:,cfc11_ind,iblock))
         STF_MODULE(:,:,cfc11_ind,iblock) = &
            PV * (CFC_surf_sat - SURF_VALS)

         PV = XKW_ICE * sqrt(660.0_r8 / CFC12_SCHMIDT)
         CFC_SFLUX_TAVG(:,:,10,iblock) = PV
         CFC_surf_sat = AP_USED(:,:,iblock) * CFC12_SOL_0 * pCFC12
         CFC_SFLUX_TAVG(:,:,11,iblock) = CFC_surf_sat
         SURF_VALS = p5*(SURF_VALS_OLD(:,:,cfc12_ind,iblock) + &
                         SURF_VALS_CUR(:,:,cfc12_ind,iblock))
         STF_MODULE(:,:,cfc12_ind,iblock) = &
            PV * (CFC_surf_sat - SURF_VALS)
      elsewhere
         STF_MODULE(:,:,cfc11_ind,iblock) = c0
         STF_MODULE(:,:,cfc12_ind,iblock) = c0
      endwhere

   end do
   !$OMP END PARALLEL DO

   call timer_stop(cfc_sflux_timer)

!-----------------------------------------------------------------------
!EOC

 end subroutine cfc_set_sflux

!***********************************************************************
!BOP
! !IROUTINE: comp_pcfc
! !INTERFACE:

 subroutine comp_pcfc(iblock, LAND_MASK, data_ind, pCFC11, pCFC12)

! !DESCRIPTION:
!  Compute atmospheric mole fractions of CFCs
!  Linearly interpolate hemispheric values to current time step
!  Spatial pattern is determined by :
!     Northern Hemisphere value is used North of 10N
!     Southern Hemisphere value is used North of 10S
!     Linear Interpolation (in latitude) is used between 10N & 10S

! !REVISION HISTORY:
!  same as module

! !USES:

   use grid, only : TLATD
   use constants, only : c10
   use time_management, only : iyear, iday_of_year, frac_day, days_in_year

! !INPUT PARAMETERS:

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK          ! land mask for this block

   integer (int_kind) :: &
      iblock          ! block index

! !INPUT/OUTPUT PARAMETERS:

   integer (int_kind) :: &
      data_ind  ! data_ind is the index into data for current timestep, 
                ! i.e data_ind is largest integer less than pcfc_data_len s.t.
                !  pcfc_date(i) <= iyear + (iday_of_year-1+frac_day)/days_in_year
                !                  - model_year + data_year
                !  note that data_ind is always strictly less than pcfc_data_len
                !  and is initialized to -1 before the first call

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      pCFC11,       & ! atmospheric CFC11 mole fraction (pmol/mol)
      pCFC12          ! atmospheric CFC11 mole fraction (pmol/mol)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      i, j              ! loop indices

   real (r8) :: &
      mapped_date,    & ! date of current model timestep mapped to data timeline
      weight,         & ! weighting for temporal interpolation
      pcfc11_nh_curr, & ! pcfc11_nh for current time step (pmol/mol)
      pcfc11_sh_curr, & ! pcfc11_sh for current time step (pmol/mol)
      pcfc12_nh_curr, & ! pcfc12_nh for current time step (pmol/mol)
      pcfc12_sh_curr    ! pcfc12_sh for current time step (pmol/mol)

!-----------------------------------------------------------------------
!  Generate mapped_date and check to see if it is too large.
!  The check for mapped_date being too small only needs to be done
!  on the first time step.
!-----------------------------------------------------------------------


   mapped_date = iyear + (iday_of_year-1+frac_day)/days_in_year &
                 - model_year + data_year

   if (mapped_date >= pcfc_date(pcfc_data_len) + max_pcfc_extension) &
      call exit_POP(sigAbort, 'model date maps too far beyond pcfc_date(end)')

!-----------------------------------------------------------------------
!  Assume atmospheric concentrations are zero before record.
!-----------------------------------------------------------------------

   if (mapped_date < pcfc_date(1)) then
      pCFC11 = c0
      pCFC12 = c0
      data_ind = 1
      return
   endif

!-----------------------------------------------------------------------
!  On first time step, perform linear search to find data_ind.
!-----------------------------------------------------------------------

   if (data_ind == -1) then
      do data_ind = pcfc_data_len-1,1,-1
         if (mapped_date >= pcfc_date(data_ind)) exit
      end do
   endif

!-----------------------------------------------------------------------
!  See if data_ind need to be updated,
!  but do not set it to pcfc_data_len.
!-----------------------------------------------------------------------

   if (data_ind < pcfc_data_len-1) then
      if (mapped_date >= pcfc_date(data_ind+1)) data_ind = data_ind + 1
   endif

!-----------------------------------------------------------------------
!  Generate hemisphere values for current time step.
!-----------------------------------------------------------------------

   weight = (mapped_date - pcfc_date(data_ind)) &
            / (pcfc_date(data_ind+1) - pcfc_date(data_ind))

   pcfc11_nh_curr = &
      weight * pcfc11_nh(data_ind+1) + (c1-weight) * pcfc11_nh(data_ind)

   pcfc11_sh_curr = &
      weight * pcfc11_sh(data_ind+1) + (c1-weight) * pcfc11_sh(data_ind)

   pcfc12_nh_curr = &
      weight * pcfc12_nh(data_ind+1) + (c1-weight) * pcfc12_nh(data_ind)

   pcfc12_sh_curr = &
      weight * pcfc12_sh(data_ind+1) + (c1-weight) * pcfc12_sh(data_ind)

!-----------------------------------------------------------------------
!     Merge hemisphere values.
!-----------------------------------------------------------------------

   do j = 1, ny_block
      do i = 1, nx_block
         if (LAND_MASK(i,j)) then
            if (TLATD(i,j,iblock) < -c10) then
               pCFC11(i,j) = pcfc11_sh_curr
               pCFC12(i,j) = pcfc12_sh_curr
            else if (TLATD(i,j,iblock) > c10) then
               pCFC11(i,j) = pcfc11_nh_curr
               pCFC12(i,j) = pcfc12_nh_curr
            else
               pCFC11(i,j) = pcfc11_sh_curr + (TLATD(i,j,iblock)+c10) &
                  * 0.05_r8 * (pcfc11_nh_curr - pcfc11_sh_curr)
               pCFC12(i,j) = pcfc12_sh_curr + (TLATD(i,j,iblock)+c10) &
                  * 0.05_r8 * (pcfc12_nh_curr - pcfc12_sh_curr)
            endif
         endif
      end do
   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_pcfc

!***********************************************************************
!BOP
! !IROUTINE: comp_cfc_schmidt
! !INTERFACE:

 subroutine comp_cfc_schmidt(LAND_MASK, SST_in, CFC11_SCHMIDT, CFC12_SCHMIDT)

! !DESCRIPTION:
!  Compute Schmidt numbers of CFCs.
!  Ref : Zheng et al. (1998), JGR, Vol. 103, No. C1, pp 1375-1379
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   logical (log_kind), dimension(nx_block,ny_block), intent(in) :: &
      LAND_MASK          ! land mask for this block

   real (r8), dimension(nx_block,ny_block), intent(in) :: &
      SST_in             ! sea surface temperature (C)

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      CFC11_SCHMIDT,   & ! Schmidt number of CFC11 (non-dimensional)
      CFC12_SCHMIDT      ! Schmidt number of CFC12 (non-dimensional)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      a1_11 = 3501.8_r8,        a1_12 = 3845.4_r8, &
      a2_11 = -210.31_r8,       a2_12 = -228.95_r8, &
      a3_11 =    6.1851_r8,     a3_12 =    6.1908_r8, &
      a4_11 =   -0.07513_r8,    a4_12 =   -0.067430_r8

   real (r8), dimension(nx_block,ny_block) :: &
      SST                ! sea surface temperature (C)

!-----------------------------------------------------------------------
! Zheng's fit only uses data up to 30
! when temp exceeds 35, use 35
! CFC11 fit goes negative for temp > 42.15
!-----------------------------------------------------------------------

   SST = merge(SST_in, 35.0_r8, SST_in < 35.0_r8)

   where (LAND_MASK)
      CFC11_SCHMIDT = a1_11 + SST * (a2_11 + SST * (a3_11 + a4_11 * SST))
      CFC12_SCHMIDT = a1_12 + SST * (a2_12 + SST * (a3_12 + a4_12 * SST))
   elsewhere
      CFC11_SCHMIDT = c0
      CFC12_SCHMIDT = c0
   endwhere

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_cfc_schmidt

!***********************************************************************
!BOP
! !IROUTINE: comp_cfc_sol_0
! !INTERFACE:

 subroutine comp_cfc_sol_0(LAND_MASK, SST, SSS, CFC11_SOL_0, CFC12_SOL_0)

! !DESCRIPTION:
!  Compute solubilities of CFCs at 1 atm.
!  Ref : Warner & Weiss (1985), Deep Sea Reasearch,
!        Vol 32, No. 12, pp. 1485-1497 (Table 5)
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use constants, only: T0_Kelvin

! !INPUT PARAMETERS:

   logical (log_kind), dimension(nx_block,ny_block) :: &
      LAND_MASK          ! land mask for this block

   real (r8), dimension(nx_block,ny_block) :: &
      SST,             & ! sea surface temperature (C)
      SSS                ! sea surface salinity (psu)

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block), intent(out) :: &
      CFC11_SOL_0,     & ! solubility of CFC11 at 1 atm (mol/l/atm)
      CFC12_SOL_0        ! solubility of CFC12 at 1 atm (mol/l/atm)

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------

   real (r8), parameter :: &
      a1_11 = -229.9261_r8,     a1_12 = -218.0971_r8, &
      a2_11 =  319.6552_r8,     a2_12 =  298.9702_r8, &
      a3_11 =  119.4471_r8,     a3_12 =  113.8049_r8, &
      a4_11 =   -1.39165_r8,    a4_12 =   -1.39165_r8, &
      b1_11 =   -0.142382_r8,   b1_12 =   -0.143566_r8, &
      b2_11 =    0.091459_r8,   b2_12 =    0.091015_r8, &
      b3_11 =   -0.0157274_r8,  b3_12 =   -0.0153924_r8

   real (r8), dimension(nx_block,ny_block) :: &
      SSTKp01  ! .01 * sea surface temperature (in Kelvin)

!-----------------------------------------------------------------------

   SSTKp01 = merge( ((SST + T0_Kelvin)* 0.01_r8), c1, LAND_MASK)

   where (LAND_MASK)
      CFC11_SOL_0 = EXP(a1_11 + a2_11 / SSTKp01 &
                        + a3_11 * LOG(SSTKp01) + a4_11 * SSTKp01 ** 2 &
                        + SSS * (b1_11 + SSTKp01 * (b2_11 + b3_11 * SSTKp01)))
      CFC12_SOL_0 = EXP(a1_12 + a2_12 / SSTKp01 &
                        + a3_12 * LOG(SSTKp01) + a4_12 * SSTKp01 ** 2 &
                        + SSS * (b1_12 + SSTKp01 * (b2_12 + b3_12 * SSTKp01)))
   elsewhere
      CFC11_SOL_0 = c0
      CFC12_SOL_0 = c0
   endwhere

!-----------------------------------------------------------------------
!EOC

 end subroutine comp_cfc_sol_0

!***********************************************************************
!BOP
! !IROUTINE: cfc_tavg_forcing
! !INTERFACE:

 subroutine cfc_tavg_forcing

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

   !$OMP PARALLEL DO PRIVATE(iblock)

   do iblock = 1, nblocks_clinic
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,1,iblock),tavg_CFC_IFRAC,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,2,iblock),tavg_CFC_XKW,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,3,iblock),tavg_CFC_ATM_PRESS,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,4,iblock),tavg_pCFC11,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,5,iblock),tavg_pCFC12,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,6,iblock),tavg_CFC11_SCHMIDT,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,7,iblock),tavg_CFC12_SCHMIDT,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,8,iblock),tavg_CFC11_PV,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,9,iblock),tavg_CFC11_surf_sat,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,10,iblock),tavg_CFC12_PV,iblock,1)
         call accumulate_tavg_field(CFC_SFLUX_TAVG(:,:,11,iblock),tavg_CFC12_surf_sat,iblock,1)
   end do

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!EOC

 end subroutine cfc_tavg_forcing

!***********************************************************************

end module cfc_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
