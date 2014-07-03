!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module moby_mod

!BOP
! !MODULE: moby_mod
!
! !DESCRIPTION:
!

! !REVISION HISTORY:

!  SVN:$Id:  $

!-----------------------------------------------------------------------
!  variables/subroutines/function used from other modules
!  The following are used extensively in this moby module, so are used at
!  the module level. The use statements for variables that are only needed
!  locally are located at the module subprogram level.
!-----------------------------------------------------------------------

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_GridHorzMod
   use POP_FieldMod
   use POP_HaloMod

   use kinds_mod
   use constants
   use communicate
   use broadcast
   use global_reductions
   use blocks
   use domain_size
   use domain
   use exit_mod
   use prognostic
   use grid
   use io
   use io_types
   use io_tools
   use moby_parms
   use tavg
   use timers
   use passive_tracer_tools
   use named_field_mod
   use forcing_tools
   use time_management
   use registry
   use co2calc
#ifdef CCSMCOUPLED
   use POP_MCT_vars_mod
   use shr_strdata_mod
   use shr_sys_mod
#endif


   implicit none
   save
   private

!-----------------------------------------------------------------------
!  public/private declarations
!-----------------------------------------------------------------------

   public :: &
      moby_init,                  &
      moby_reset,                 &
      moby_set_sflux,             &
      moby_set_interior,          &
      moby_set_interior_3D,       &
      moby_tavg_forcing,          &
      moby_tracer_cnt,            &
      moby_tracer_ref_val,        &
      moby_write_restart,         &
      moby_qsw_distrb_const,      &
      lmoby,                      &
      POP_mobyFinal,              &
      POP_mobySendTime
      
!-----------------------------------------------------------------------
!  module variables 
!-----------------------------------------------------------------------

   integer ::                     &
      myThid,                     &
      moby_NumPTRACERS 

   integer (int_kind), parameter :: &
      moby_tracer_cnt = 47
   integer (int_kind), parameter :: &
      moby_max_np = 100,            &
      moby_max_nz = 100           

!-----------------------------------------------------------------------
!  Massive array for storing 3D passive-tracer tendencies computed in
!  darwin_forcing
!-----------------------------------------------------------------------
   real (r8), dimension(:,:,:,:), allocatable :: & 
      PT3D_MOBY
!-----------------------------------------------------------------------
!  Not-so-massive array for storing CO2 flux computed in darwin_forcing
!-----------------------------------------------------------------------
   real (r8), dimension(:,:,:), allocatable :: & 
      FLXCO2_MOBY

!-----------------------------------------------------------------------
!  flags controlling which portion of code are executed
!  usefull for debugging
!-----------------------------------------------------------------------

  logical (log_kind) :: &
     moby_lsource_sink, &
     moby_lflux_gas_o2, &
     moby_lflux_gas_co2

  logical (log_kind), dimension(:,:,:), allocatable :: &
     LAND

!-----------------------------------------------------------------------
!  relative tracer indices
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!  derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(moby_tracer_cnt) :: &
      moby_ind_name_table

 
!-----------------------------------------------------------------------
!  options for forcing of gas fluxes
!-----------------------------------------------------------------------

   integer (int_kind), parameter :: &
      gas_flux_forcing_iopt_drv   = 1,   &
      gas_flux_forcing_iopt_file  = 2,   &
      atm_co2_iopt_const          = 1,   &
      atm_co2_iopt_drv_prog       = 2,   &
      atm_co2_iopt_drv_diag       = 3

   integer (int_kind) :: &
      gas_flux_forcing_iopt,             &
      atm_co2_iopt

   real (r8)       :: &
      atm_co2_const  ! value of atmospheric co2 (ppm, dry-air, 1 atm)

   character(char_len) :: &
      gas_flux_forcing_file    ! file containing gas flux forcing fields

!-----------------------------------------------------------------------

   real (r8) :: &
      rest_time_inv_surf,  & ! inverse restoring timescale at surface
      rest_time_inv_deep,  & ! inverse restoring timescale at depth
      rest_z0,             & ! shallow end of transition regime
      rest_z1                ! deep end of transition regime

   type(tracer_read) :: &
      gas_flux_fice,       & ! ice fraction for gas fluxes
      gas_flux_ws,         & ! wind speed for gas fluxes
      gas_flux_ap,         & ! atmospheric pressure for gas fluxes
      fesedflux_input        ! namelist input for iron_fl
  
!-----------------------------------------------------------------------
!  module variables related to ph computations
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:), allocatable, target :: &
      PH_PREV           ! computed ph from previous time step


!-----------------------------------------------------------------------
!  restoring climatologies for nutrients
!-----------------------------------------------------------------------

   logical (log_kind) :: &
      moby_lrest_po4,  & ! restoring on po4 
      moby_lrest_no3,  & ! restoring on no3 
      moby_lrest_sio3    ! restoring on sio3 

   real (r8), dimension(km) :: &
      moby_nutr_rest_time_inv ! inverse restoring time scale for nutrients (1/secs)

   real (r8), dimension(:,:,:,:), allocatable, target :: &
      PO4_CLIM, NO3_CLIM, SiO3_CLIM
   real (r8), dimension(:,:,:,:), allocatable, target :: &
      FESEDFLUX

   character(char_len) :: &
      moby_nutr_rest_file               ! file containing nutrient fields

   real (r8), dimension(:,:,:), allocatable, target :: &
      NUTR_RESTORE_RTAU            ! inverse restoring timescale for variable
                                   ! interior restoring

   integer (int_kind), dimension(:,:,:), allocatable :: &
      NUTR_RESTORE_MAX_LEVEL       ! maximum level for applying variable
                                   ! interior restoring

   real (r8), dimension(:,:,:,:), allocatable :: &
      INTERP_WORK                  ! temp array for interpolate_forcing output

   type(forcing_monthly_every_ts) :: &
      moby_iron_flux,                 & ! iron component of surface dust flux
      moby_fice_file,                 & ! ice fraction, if read from file
      moby_xkw_file,                  & ! a * wind-speed ** 2, if read from file
      moby_ap_file                      ! atmoshperic pressure, if read from file

   character(char_len) :: &
      ndep_data_type               ! type of ndep forcing

   integer (int_kind) :: &
      moby_ndep_shr_stream_year_first, & ! first year in stream to use
      moby_ndep_shr_stream_year_last,  & ! last year in stream to use
      moby_ndep_shr_stream_year_align    ! align moby_ndep_shr_stream_year_first with this model year

   integer (int_kind), parameter :: &
      moby_ndep_shr_stream_var_cnt = 2, & ! number of variables in ndep shr_stream
      moby_ndep_shr_stream_no_ind  = 1, & ! index for NO forcing
      moby_ndep_shr_stream_nh_ind  = 2    ! index for NH forcing

   character(char_len) :: &
      moby_ndep_shr_stream_file          ! file containing domain and input data

   real (r8) :: &
      moby_ndep_shr_stream_scale_factor  ! unit conversion factor

#ifdef CCSMCOUPLED
   type(shr_strdata_type) :: ndep_sdat ! input data stream for ndep
#endif

!-----------------------------------------------------------------------
!  derived type for implicit handling of sinking particulate matter
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!  define tavg id for 2d fields related to surface fluxes
!-----------------------------------------------------------------------

   integer (int_kind) ::      &
      tavg_MOBY_IFRAC,        &! tavg id for ice fraction
      tavg_MOBY_IFRAC_2,      &! tavg id for duplicate of ice fraction
      tavg_MOBY_XKW,          &! tavg id for xkw
      tavg_MOBY_XKW_2,        &! tavg id for duplicate of xkw
      tavg_MOBY_ATM_PRESS,    &! tavg id for atmospheric pressure
      tavg_MOBY_DIC_GAS_FLUX, &! tavg id for dic flux
      tavg_IRON_FLUX           ! tavg id for dust flux

!-----------------------------------------------------------------------
!  define tavg id for nonstandard 3d fields
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_PAR_avg        ! tavg id for available radiation avg over mixed layer


!-----------------------------------------------------------------------
!  define tavg id for nonstandard 2d fields
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      tavg_O2_ZMIN,      &! tavg id for vertical minimum of O2
      tavg_O2_ZMIN_DEPTH  ! tavg id for depth of vertical minimum of O2

!-----------------------------------------------------------------------
!  define tavg id for nonstandard 3d fields
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     tavg_P_iron_FLUX_IN   ! tavg id for p_iron flux into cell

!-----------------------------------------------------------------------
!  define tavg id for MORE nonstandard 3d fields
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!  define tavg id for MORE nonstandard 3d fields
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!  define array for holding flux-related quantities that need to be time-averaged
!  this is necessary since the forcing routines are called before tavg flags
!-----------------------------------------------------------------------

   real (r8), dimension(:,:,:,:), allocatable :: &
      MOBY_SFLUX_TAVG

!-----------------------------------------------------------------------
!  average surface tracer value related variables
!  used as reference value for virtual flux computations
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      alk_ind,           &! relative index for ALK
      dic_ind             ! relative index for DIC

   logical (log_kind), dimension(moby_tracer_cnt) :: &
      vflux_flag                ! which tracers get virtual fluxes applied

   integer (int_kind) :: &
      moby_comp_surf_avg_flag   ! time flag id for computing average
                                ! surface tracer values

   real (r8), dimension(moby_tracer_cnt)  :: &
      surf_avg                                ! average surface tracer values

   logical (log_kind) :: &
      moby_qsw_distrb_const

!-----------------------------------------------------------------------
!  iron patch fertilization
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!  timers
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      moby_shr_strdata_advance_timer,        &
      moby_interior_timer,                   &
      moby_interior_3D_timer,                &
      moby_sflux_timer

!-----------------------------------------------------------------------
!  named field indices
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      sflux_co2_nf_ind   = 0       ! air-sea co2 gas flux

!-----------------------------------------------------------------------

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic) :: &
      PAR_out           ! photosynthetically available radiation (W/m^2)

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  parameters for CESM/MOBY scaling
!-----------------------------------------------------------------------
   real (r8), parameter :: m3percm3 = 1.0e-6_r8
   real (r8), parameter :: cm3perm3 = 1.0e6_r8

!-----------------------------------------------------------------------
!  darwin-side passive-tracer indices
!-----------------------------------------------------------------------
   integer (int_kind), dimension(moby_max_nz) ::   &
      moby_iZooP,         &
      moby_iZooN,         &                
      moby_iZooFe,        &
      moby_iZooSi

   real(r8), dimension(moby_max_np) ::   &
      moby_R_NP,         &
      moby_R_FeP,        &                
      moby_R_SiP,        &
      moby_R_PC

   integer (int_kind) ::   &
      moby_nzmax,         &
      moby_npmax,         &
      moby_iphy,          &
      moby_iPO4,          &
      moby_iDOP,          &
      moby_iPOP,          &
      moby_iNO3,          &
      moby_iNH4,          &
      moby_iNO2,          &
      moby_iDON,          &
      moby_iPON,          &
      moby_iFeT,          &
      moby_iDOFe,         &
      moby_iPOFe,         &
      moby_iPOSi,         &
      moby_iSi

!-----------------------------------------------------------------------
!  other darwin-side variables
!-----------------------------------------------------------------------

   real(r8)  ::   &
      moby_parconv         !conversion from W/m2 to uEin/m2/s
  
!-----------------------------------------------------------------------
!  time output strings
!-----------------------------------------------------------------------

   character (char_len) ::  &
      time_string,          &
      date_string

!-----------------------------------------------------------------------
!  misc
!-----------------------------------------------------------------------
 
   logical (log_kind) ::  &
      lmoby,              &!is moby used?
      lmoby_continue       !communicate restarts with MOBY models

   integer (int_kind) ::  &
      ishare_val_put = 0, &
      ishare_val_get = 1

!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: moby_init
! !INTERFACE:

 subroutine moby_init(init_ts_file_fmt, read_restart_filename, &
                      tracer_d_module, TRACER_MODULE, tadvect_ctype, moby_ind_begin, errorCode)

! !DESCRIPTION:
!  Initialize moby tracer module. This involves setting metadata, reading
!  the module namelist, setting initial conditions, setting up forcing,
!  and defining additional tavg variables.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      init_ts_file_fmt,    & ! format (bin or nc) for input file
      read_restart_filename  ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

  !type (tracer_field), dimension(moby_tracer_cnt), intent(inout) :: &
   type (tracer_field), dimension(:), intent(inout) :: &
      tracer_d_module   ! descriptors for each tracer

  !real (r8), dimension(nx_block,ny_block,km,moby_tracer_cnt,3,max_blocks_clinic), &
   real (r8), dimension(:,:,:,:,:,:), &
      intent(inout) :: TRACER_MODULE

   integer (int_kind), intent(in) :: moby_ind_begin

! !OUTPUT PARAMETERS:

  !character (char_len), dimension(moby_tracer_cnt), intent(out) :: &
   character (char_len), dimension(:), intent(out) :: &
      tadvect_ctype     ! advection method for moby tracers

   integer (POP_i4), intent(out) :: &
      errorCode

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


 end subroutine moby_init

!***********************************************************************
!BOP
! !IROUTINE: moby_extract_surf_avg
! !INTERFACE:

 subroutine moby_extract_surf_avg(init_moby_init_file_fmt, &
                                  moby_restart_filename)

! !DESCRIPTION:
!  Extract average surface values from restart file.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      init_moby_init_file_fmt, & ! file format (bin or nc)
      moby_restart_filename      ! file name for restart file

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


 end subroutine moby_extract_surf_avg

!***********************************************************************
!BOP
! !IROUTINE: moby_init_tavg
! !INTERFACE:

 subroutine moby_init_tavg

! !DESCRIPTION:
!  call define_tavg_field for nonstandard tavg fields
!
! !REVISION HISTORY:
!  same as module
!

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


 end subroutine moby_init_tavg

!***********************************************************************
!BOP
! !IROUTINE: moby_set_interior_3D
! !INTERFACE:

 subroutine moby_set_interior_3D(TEMP_OLD, TEMP_CUR, SALT_OLD, SALT_CUR, &
    TRACER_MODULE_OLD, TRACER_MODULE_CUR)


! !DESCRIPTION:
!  Compute time derivatives for moby state variables. Because of the 3D nature
!  of the moby darwin models' design, this routine is called for all k-levels
!  at once, and tendencies are stored for later access level-by-level in baroclinic.F90.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,max_blocks_clinic), intent(in) :: &
      TEMP_OLD,          &! old potential temperature (C)
      TEMP_CUR,          &! current potential temperature (C)
      SALT_OLD,          &! old salinity (msu)
      SALT_CUR            ! current salinity (msu)

  !real (r8), dimension(nx_block,ny_block,km,moby_tracer_cnt,max_blocks_clinic), intent(in) :: &
   real (r8), dimension(:,:,:,:,:), intent(in) :: &
      TRACER_MODULE_OLD, &! old moby tracer values
      TRACER_MODULE_CUR   ! current moby tracer values

! !OUTPUT PARAMETERS:


!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


 end subroutine moby_set_interior_3D

!***********************************************************************
!BOP
! !IROUTINE: moby_set_interior
! !INTERFACE:

 subroutine moby_set_interior(k, DTRACER_MODULE, this_block)

! !DESCRIPTION:
!  Access the previously computed moby time derivatives 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      k                   ! vertical level index

   type (block), intent(in) :: &
      this_block          ! block info for the current block

! !OUTPUT PARAMETERS:

  !real (r8), dimension(nx_block,ny_block,moby_tracer_cnt), intent(out) :: &
   real (r8), dimension(:,:,:), intent(out) :: &
      DTRACER_MODULE      ! previously computed source/sink terms for this k level

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


 end subroutine moby_set_interior


!***********************************************************************
!BOP
! !IROUTINE: moby_reset
! !INTERFACE:

 subroutine moby_reset(TRACER_MODULE, bid)

! !DESCRIPTION:
!  Reset iron based upon darwin_fe_chem computations. 
!
!  This subroutine performs the function of the second call to darwin_fe_chem 
!  in darwin_forcing.  The routine is necessary because in CESM, darwin_forcing
!  is modified so that it only computes and returns tracer tendencies -- it
!  does not update the tracers.  So in CESM, a second call to darwin_fe_chem
!  in darwin_forcing would act upon the the *old* value of the iFeT tracer, 
!  not the new one, which is the intended action. 
! 
!  Instead, the second call to darwin_fe_chem is done in this routine, moby_reset,
!  after the iFeT tracer has been updated by the base POP2 model.
!
!
! !REVISION HISTORY:
!  same as module

! !USES:

   use time_management, only: mix_pass, c2dtt
   use prognostic, only:  newtime

! !INPUT PARAMETERS:

   integer(int_kind), intent(in) :: bid

! !INPUT/OUTPUT PARAMETERS:

  !real(r8), dimension(nx_block,ny_block,km,moby_tracer_cnt), intent(inout) :: &
   real(r8), dimension(:,:,:,:), intent(inout) :: &
      TRACER_MODULE      ! moby tracers

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


 end subroutine moby_reset

!***********************************************************************
!BOP
! !IROUTINE: moby_init_sflux
! !INTERFACE:

 subroutine moby_init_sflux

! !DESCRIPTION:
!  Initialize surface flux computations for moby tracer module.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


 end subroutine moby_init_sflux

!***********************************************************************
!BOP
! !IROUTINE: moby_init_interior_restore
! !INTERFACE:

 subroutine moby_init_interior_restore

! !DESCRIPTION:
!  Initialize interior restoring computations for moby tracer module.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


 end subroutine moby_init_interior_restore

!***********************************************************************
!BOP
! !IROUTINE: moby_set_sflux
! !INTERFACE:

 subroutine moby_set_sflux(SHF_QSW_RAW, SHF_QSW,                 &
                           U10_SQR,IFRAC,PRESS,                  &
                           SURF_VALS_OLD,SURF_VALS_CUR,STF_MODULE)
   use shr_pio_mod, only : shr_pio_getiotype, shr_pio_getiosys

! !DESCRIPTION:
!  Compute surface fluxes for moby tracer module.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      SHF_QSW_RAW,  &! penetrative solar heat flux, from coupler (degC*cm/s)
      SHF_QSW,      &! SHF_QSW used by physics, may have diurnal cylce imposed (degC*cm/s)
      U10_SQR,      &! 10m wind speed squared (cm/s)**2
      IFRAC,        &! sea ice fraction (non-dimensional)
      PRESS          ! sea level atmospheric pressure (dyne/cm**2)

  !real (r8), dimension(nx_block,ny_block,moby_tracer_cnt,max_blocks_clinic), &
   real (r8), dimension(:,:,:,:), &
      intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR ! module tracers

! !INPUT/OUTPUT PARAMETERS:

  !real (r8), dimension(nx_block,ny_block,moby_tracer_cnt,max_blocks_clinic), &
   real (r8), dimension(:,:,:,:), &
      intent(inout) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


 end subroutine moby_set_sflux

!*****************************************************************************
!BOP
! !IROUTINE: moby_tavg_forcing
! !INTERFACE:

 subroutine moby_tavg_forcing(STF_MODULE)

! !DESCRIPTION:
!  Accumulate non-standard forcing related tavg variables.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

 !real (r8), dimension(nx_block,ny_block,moby_tracer_cnt,max_blocks_clinic), &
  real (r8), dimension(:,:,:,:), &
     intent(in) :: STF_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


 end subroutine moby_tavg_forcing

!*****************************************************************************
!BOP
! !IROUTINE: moby_comp_surf_avg
! !INTERFACE:

 subroutine moby_comp_surf_avg(SURF_VALS_OLD,SURF_VALS_CUR)

! !DESCRIPTION:
!  compute average surface tracer values
!
!  avg = sum(SURF_VAL*TAREA) / sum(TAREA)
!  with the sum taken over ocean points only
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

 !real (r8), dimension(nx_block,ny_block,moby_tracer_cnt,max_blocks_clinic), &
  real (r8), dimension(:,:,:,:), &
      intent(in) :: SURF_VALS_OLD, SURF_VALS_CUR

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


 end subroutine moby_comp_surf_avg

!*****************************************************************************
!BOP
! !IROUTINE: moby_write_restart
! !INTERFACE:

 subroutine moby_write_restart(restart_file, action)

! !DESCRIPTION:
!  write auxiliary fields & scalars to restart files
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character(*), intent(in) :: action

! !INPUT/OUTPUT PARAMETERS:

   type (datafile), intent (inout)  :: restart_file

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


 end subroutine moby_write_restart

!*****************************************************************************
!BOP
! !IROUTINE: moby_tracer_ref_val
! !INTERFACE:

 function moby_tracer_ref_val(ind)

! !DESCRIPTION:
!  return reference value for tracers using virtual fluxes
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: ind

! !OUTPUT PARAMETERS:

   real (r8) :: moby_tracer_ref_val

!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


! dummy value -- needed for moby_mod "stub" to compile
      moby_tracer_ref_val = c0

!-----------------------------------------------------------------------
!EOC

 end function moby_tracer_ref_val



!***********************************************************************
!BOP
! !IROUTINE: POP_mobyInit1
! !INTERFACE:

 subroutine POP_mobyInit1 (log_filename)


! !DESCRIPTION:
!  This routine initializes the pop2 MITgcm_moby interface
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:

   character (char_len), intent(in) ::  &
      log_filename  ! output log filename for moby

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


 end subroutine POP_mobyInit1

!***********************************************************************
!BOP
! !IROUTINE: POP_mobyInit2
! !INTERFACE:

 subroutine POP_mobyInit2


! !DESCRIPTION:
!  This routine initializes the pop2 MITgcm_moby interface
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


 end subroutine POP_mobyInit2


!***********************************************************************
!BOP
! !IROUTINE: POP_mobySurfaceForcingSet
! !INTERFACE:

 subroutine POP_mobySurfaceForcingSet (PAR_out,U10_SQR,IFRAC,PRESS,IRON)


! !DESCRIPTION:
!  This routine sends POP2 surface forcing information to the MITgcm_moby interface
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      PAR_out,      &! photosynthetically available radiation (convert to W/m^2)
      U10_SQR,      &! 10m wind speed squared (cm/s)**2
      IFRAC,        &! sea ice fraction (non-dimensional)
      PRESS,        &! sea level atmospheric pressure (dyne/cm**2)
      IRON           ! bio-available iron flux  (mol/m**2/s)

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


 end subroutine POP_mobySurfaceForcingSet


!***********************************************************************
!BOP
! !IROUTINE: POP_mobyMeanArea
! !INTERFACE:

 subroutine POP_mobyMeanArea(TEMP_OLD, TEMP_CUR, SALT_OLD, SALT_CUR, &
    TRACER_MODULE_OLD, TRACER_MODULE_CUR)

! !DESCRIPTION:
!  Compute surface mean of tracers and send to darwin model.
!  Called by moby_set_interior_3D.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,max_blocks_clinic), intent(in) :: &
      TEMP_OLD,          &! old potential temperature (C)
      TEMP_CUR,          &! current potential temperature (C)
      SALT_OLD,          &! old salinity (msu)
      SALT_CUR            ! current salinity (msu)

  !real (r8), dimension(nx_block,ny_block,km,moby_tracer_cnt,max_blocks_clinic), intent(in) :: &
   real (r8), dimension(:,:,:,:,:), intent(in) :: &
      TRACER_MODULE_OLD, &! old moby tracer values
      TRACER_MODULE_CUR   ! current moby tracer values

! !OUTPUT PARAMETERS:


!EOP
!BOC
!-----------------------------------------------------------------------
!  local variables
!-----------------------------------------------------------------------


 end subroutine POP_mobyMeanArea

!***********************************************************************
!BOP
! !IROUTINE: POP_mobyCons
! !INTERFACE:

 subroutine POP_mobyCons 


! !DESCRIPTION:
!  This routine calls the moby conservation subroutines; it replaces darwin_cons.F and
!   its quota equivalent
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:

! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


 end subroutine POP_mobyCons

!***********************************************************************

!BOP
! !IROUTINE: POP_mobySendTime
! !INTERFACE:

 subroutine POP_mobySendTime


! !DESCRIPTION:
!  This routine sends pop2 time information to the MITgcm_moby interface
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


 end subroutine POP_mobySendTime

!***********************************************************************
!BOP
! !IROUTINE: POP_mobyFinal
! !INTERFACE:

 subroutine POP_mobyFinal


! !DESCRIPTION:
!  This routine sends the finalize signal to the MITgcm_moby interface
!
! !REVISION HISTORY:
!  same as module
!
! !OUTPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


 end subroutine POP_mobyFinal

!***********************************************************************
!BOP
! !IROUTINE: POP_mobyConsistencyChecks
! !INTERFACE:

 subroutine POP_mobyConsistencyChecks(myThid)


! !DESCRIPTION:
!  This routine checks for incompatible/unsupported CESM/MOBY options
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:
   integer (int_kind) ::  &
      myThid

! !OUTPUT PARAMETERS:


!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


 end subroutine POP_mobyConsistencyChecks

!***********************************************************************
!BOP
! !IROUTINE: moby_tracer_mean_area
! !INTERFACE:

!function moby_tracer_mean_area (MY_TRACER, k)

! end function moby_tracer_mean_area



  subroutine moby_global_tracer_volume (tracer_mean, modified_volume_t)

! !DESCRIPTION:
!
!  This subroutine computes the global sum of volume * tracer
!  (returned in tracer_mean) and the global ocean volume at T points
!  at curtime (returned in modified_volume_t).
!  This routine is nearly identical to diag_for_tracer_budgets, which
!  cannot be used in this module without sorting out a circular dependency.

! !REVISION HISTORY:
!  based on subroutine diag_for_tracer_budgets in budget_diagnostics.F90


! !INPUT PARAMETERS:


! !OUTPUT PARAMETERS:

   real (r8), dimension(nt), intent(out) :: tracer_mean
   real (r8),                intent(out) :: modified_volume_t

!EOP
!BOC


   end subroutine moby_global_tracer_volume
!*****************************************************************************

!BOP
! !IROUTINE: moby_lookup_longname
! !INTERFACE:

!function moby_lookup_longname(name)


!end function moby_lookup_longname


 end module moby_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
