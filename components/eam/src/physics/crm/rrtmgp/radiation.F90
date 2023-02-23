#define _IDX321(i, j, k, nx, ny, nz) (nx * (ny * (k - 1) + (j - 1)) + i)
module radiation
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to RRTMGP
!
! Revision history:
! Nov 2017  B. Hillman  Initial version, adapted from RRTMG interface.
!---------------------------------------------------------------------------------

   ! E3SM-specific modules that are used throughout this module. An effort was made
   ! to keep imports as local as possible, so we only load a few of these at the
   ! module (rather than the subroutine) level.
   use iso_c_binding
   use shr_kind_mod,     only: r8=>shr_kind_r8, cl=>shr_kind_cl
   use ppgrid,           only: pcols, pver, pverp, begchunk, endchunk
   use cam_abortutils,   only: endrun
   use scamMod,          only: scm_crm_mode, single_column
   use rad_constituents, only: N_DIAG
   use radconstants,     only: &
      nswbands, nlwbands, &
      get_band_index_sw, get_band_index_lw, test_get_band_index, &
      get_sw_spectral_midpoints, get_lw_spectral_midpoints, &
      get_sw_spectral_boundaries, &
      rrtmg_to_rrtmgp_swbands
   use cam_history_support, only: add_hist_coord
   use physconst, only: cpair, cappa, stebol

   ! RRTMGP gas optics object to store coefficient information. This is imported
   ! here so that we can make the k_dist objects module data and only load them
   ! once.
   use rrtmgp_interface, only: &
      rrtmgp_initialize, rrtmgp_finalize, &
      rrtmgp_run_sw, rrtmgp_run_lw, &
      get_min_temperature, get_max_temperature, &
      get_gpoint_bands_sw, get_gpoint_bands_lw, &
      nswgpts, nlwgpts

   ! Use my assertion routines to perform sanity checks
   use assertions, only: assert, assert_valid, assert_range

   use radiation_state, only: ktop, kbot, nlev_rad
   use radiation_utils, only: handle_error, clip_values, &
      fluxes_t, initialize_fluxes, free_fluxes, reset_fluxes, &
      expand_day_fluxes, get_gas_vmr

   ! For MMF
   use crmdims, only: crm_nx_rad, crm_ny_rad, crm_nz

   implicit none
   private
   save

   ! Public routines provided by this module.
   ! TODO: radiation_defaultopts, radiation_setops, and radiation_printops exist
   ! only because the radiation namelist has traditionally been read at the driver
   ! level by runtime_opts. I am not sure this is the best solution, and in fact
   ! CESM seems to have gone away from this practice altogether, opting instead for
   ! individual modules to be responsible for reading their own namelists. This
   ! might be a better practice going forward, and would simplify the logic here
   ! somewhat. To do this, we would need to use the radiation_readnl routine
   ! (copied below), and add a line in runtime_opts that calls this routine, and
   ! then remove the code in runtime_opts that reads the radiation namelist
   ! variables from the CAM namelist. The calls to radiation_defaultopts,
   ! radiation_setops, and radiation_printops would also then need to be removed
   ! from runtime_ops.
   public :: &
      radiation_register,    &! registers radiation physics buffer fields
      radiation_nextsw_cday, &! calendar day of next radiation calculation
      radiation_do,          &! query which radiation calcs are done this timestep
      radiation_init,        &! calls radini
      radiation_final,       &! deallocate
      radiation_readnl,      &! read radiation namelist
      radiation_tend          ! moved from radctl.F90

   ! Counter variables for use with the CFMIP Observation Simulator Package (COSP).
   ! TODO: This seems like somewhat of an awkward way of determining when to run
   ! COSP, and it might make more sense to implement a more elegant routine to call
   ! that determines whether or not to run COSP at a given timestep (similar to the
   ! radiation_do routine in this module).
   ! TODO: move these to COSP interface instead
   integer, public, allocatable :: cosp_cnt(:)       ! counter for cosp
   integer, public              :: cosp_cnt_init = 0 ! initial value for cosp counter

   ! Declare namelist variables as module data. This also sets default values for
   ! namelist variables.
   integer :: iradsw = -1  ! freq. of shortwave radiation calc in time steps (positive)
                           ! or hours (negative).
   integer :: iradlw = -1  ! frequency of longwave rad. calc. in time steps (positive)
                           ! or hours (negative).
   integer :: irad_always = 0  ! Specifies length of time in timesteps (positive)
                               ! or hours (negative) SW/LW radiation will be
                               ! run continuously from the start of an
                               ! initial or restart run

   ! The spectralflux flag determines whether or not spectral (per band) fluxes are
   ! calculated. If true, upward and downward fluxes are calculated per band,
   ! otherwise just broadband "shortwave" and "longwave" fluxes are calculated.
   ! TODO: while it seems that setting spectralflux = .true. add spectral fluxes to
   ! the physics buffer, I do not see where these fields are added to the history
   ! buffer. It might be that the output of these fields are just handled by the
   ! physics buffer output routines (they are declared with "global" scope, so
   ! written on restart files at least), but it would be good to include them on
   ! the history buffer too so that we can annotate them with a description of what
   ! they are, rather than expecting the user to know what they are from the pbuf
   ! fields.
   logical :: spectralflux  = .false.  ! calculate fluxes (up and down) per band.

   ! Flag to indicate whether or not to use the radiation timestep for solar zenith
   ! angle calculations. If true, use the radiation timestep for all solar zenith 
   ! angle (cosz) calculations.
   ! TODO: How does this differ if value is .false.?
   logical :: use_rad_dt_cosz  = .false. 

   ! Flag to indicate whether to do aerosol optical calculations. This
   ! zeroes out the aerosol optical properties if False
   logical :: do_aerosol_rad = .true.

   ! Value for prescribing an invariant solar constant (i.e. total solar
   ! irradiance at TOA). Used for idealized experiments such as RCE.
   ! Disabled when value is less than 0.
   real(r8) :: fixed_total_solar_irradiance = -1

   ! The RRTMGP warnings are printed when the state variables need to be limted,
   ! such as when the temperature drops too low. This is not normally an issue,
   ! but in aquaplanet and RCE configurations these situations occur much more 
   ! frequently, so this flag was added to be able to disable those messages.
   logical :: rrtmgp_enable_temperature_warnings = .true.

   ! Model data that is not controlled by namelist fields specifically follows
   ! below.

   ! TODO: it seems CAM allows for multiple passes through the radiation (and other
   ! physics?) for different diagnostic purposes, but this does not seem to be well
   ! documented anywhere. We should add some documentation here on this if we are
   ! going to keep this functionality.
   character(len=4) :: diag(0:N_DIAG) =(/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ', &
                                         '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

   ! Output diagnostic brightness temperatures at the top of the
   ! atmosphere for 7 TOVS/HIRS channels (2,4,6,8,10,11,12) and 4 TOVS/MSU 
   ! channels (1,2,3,4).
   ! TODO: where are these options set?
   logical :: dohirs = .false. 
   integer :: ihirsfq = 1      ! frequency (timesteps) of brightness temperature calcs

   ! time step to use for the shr_orb_cosz calculation, if use_rad_dt_cosz set to true
   ! TODO: where is this set, and what is shr_orb_cosz? Alternative solar zenith
   ! angle calculation? What is the other behavior?
   real(r8) :: dt_avg = 0.0_r8  

   ! k-distribution coefficients files to read from. These are set via namelist
   ! variables.
   character(len=cl) :: rrtmgp_coefficients_file_sw, rrtmgp_coefficients_file_lw

   ! Band midpoints; these need to be module variables because of how cam_history works;
   ! add_hist_coord sets up pointers to these, so they need to persist.
   real(r8), target :: sw_band_midpoints(nswbands), lw_band_midpoints(nlwbands)

   ! Set name for this module (for purpose of writing output and log files)
   character(len=*), parameter :: module_name = 'radiation'

   ! Gases that we want to use in the radiative calculations. These need to be set
   ! here, because the RRTMGP kdist initialization needs to know the names of the
   ! gases before these are available via the rad_cnst interface. So, if we want to
   ! change which gases are used at some point, this needs to be changed here, or
   ! we need to come up with a more clever way of getting these from the model,
   ! since this is kind of a hack to hard-code them in here.
   character(len=3), dimension(8) :: active_gases = (/ &
      'H2O', 'CO2', 'O3 ', 'N2O', &
      'CO ', 'CH4', 'O2 ', 'N2 ' &
   /)

   ! Null-terminated C-compatible version of active gases for use with C++
   ! routines.
   character(len=len(active_gases)+1), dimension(size(active_gases)), target :: active_gases_c

   ! Stuff to generate random numbers for perturbation growth tests. This needs to
   ! be public module data because restart_physics needs to read it to write it to
   ! restart files (I think). Making this public module data may not be the best
   ! way of doing this, so maybe this should be reconsidered in the future.
   !
   ! kiss_seed_num is number of seed values to use?
   integer, public, parameter :: kiss_seed_num = 4

   ! TODO: what does this mean? This seems to be part of the additions for the
   ! perturbation growth functionality.
   integer, public, allocatable :: rad_randn_seedrst(:,:,:)

   ! Total number of physics chunks until this processor (TODO: what exactly does
   ! this mean?)
   integer, public, allocatable :: tot_chnk_till_this_prc(:,:,:)

   !============================================================================

contains

   !============================================================================

   subroutine radiation_readnl(nlfile, dtime_in)
   !-------------------------------------------------------------------------------
   ! Purpose: Read radiation_nl namelist group.
   !-------------------------------------------------------------------------------

      use namelist_utils,  only: find_group_name
      use units,           only: getunit, freeunit
      use spmd_utils,      only: mpicom, mstrid=>masterprocid, &
                                 mpi_integer, mpi_logical, &
                                 mpi_character, masterproc, &
                                 mpi_real8
      use time_manager,    only: get_step_size
      use cam_logfile,     only: iulog

      ! File containing namelist input
      character(len=*), intent(in) :: nlfile
      integer, intent(in), optional :: dtime_in

      ! Local variables
      integer :: unitn, ierr
      integer :: dtime  ! timestep size
      character(len=*), parameter :: subroutine_name = 'radiation_readnl'

      ! Variables defined in namelist
      namelist /radiation_nl/ rrtmgp_coefficients_file_lw,     &
                              rrtmgp_coefficients_file_sw,     &
                              iradsw, iradlw, irad_always,     &
                              use_rad_dt_cosz, spectralflux,   &
                              do_aerosol_rad,                  &
                              fixed_total_solar_irradiance,    &
                              rrtmgp_enable_temperature_warnings

      ! Read the namelist, only if called from master process
      ! TODO: better documentation and cleaner logic here?
      if (masterproc) then
         unitn = getunit()
         open(unitn, file=trim(nlfile), status='old')
         call find_group_name(unitn, 'radiation_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, radiation_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subroutine_name // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)
         call freeunit(unitn)
      end if

#ifdef SPMD
      ! Broadcast namelist variables
      call mpibcast(rrtmgp_coefficients_file_lw, cl, mpi_character, mstrid, mpicom, ierr)
      call mpibcast(rrtmgp_coefficients_file_sw, cl, mpi_character, mstrid, mpicom, ierr)
      call mpibcast(iradsw, 1, mpi_integer, mstrid, mpicom, ierr)
      call mpibcast(iradlw, 1, mpi_integer, mstrid, mpicom, ierr)
      call mpibcast(irad_always, 1, mpi_integer, mstrid, mpicom, ierr)
      call mpibcast(use_rad_dt_cosz, 1, mpi_logical, mstrid, mpicom, ierr)
      call mpibcast(spectralflux, 1, mpi_logical, mstrid, mpicom, ierr)
      call mpibcast(do_aerosol_rad, 1, mpi_logical, mstrid, mpicom, ierr)
      call mpibcast(fixed_total_solar_irradiance, 1, mpi_real8, mstrid, mpicom, ierr)
      call mpibcast(rrtmgp_enable_temperature_warnings, 1, mpi_logical, mstrid, mpicom, ierr)
#endif

      ! Convert iradsw, iradlw and irad_always from hours to timesteps if necessary
      if (present(dtime_in)) then
         dtime = dtime_in
      else
         dtime = get_step_size()
      end if
      if (iradsw      < 0) iradsw      = nint((-iradsw     *3600._r8)/dtime)
      if (iradlw      < 0) iradlw      = nint((-iradlw     *3600._r8)/dtime)
      if (irad_always < 0) irad_always = nint((-irad_always*3600._r8)/dtime)

      ! Print runtime options to log.
      if (masterproc) then
         write(iulog,*) 'RRTMGP radiation scheme parameters:'
         write(iulog,10) trim(rrtmgp_coefficients_file_lw), trim(rrtmgp_coefficients_file_sw), &
                         iradsw, iradlw, irad_always, &
                         use_rad_dt_cosz, spectralflux, &
                         do_aerosol_rad, fixed_total_solar_irradiance, &
                         rrtmgp_enable_temperature_warnings
      end if
   10 format('  LW coefficents file: ',                                a/, &
             '  SW coefficents file: ',                                a/, &
             '  Frequency (timesteps) of Shortwave Radiation calc:  ',i5/, &
             '  Frequency (timesteps) of Longwave Radiation calc:   ',i5/, &
             '  SW/LW calc done every timestep for first N steps. N=',i5/, &
             '  Use average zenith angle:                           ',l5/, &
             '  Output spectrally resolved fluxes:                  ',l5/, &
             '  Do aerosol radiative calculations:                  ',l5/, &
             '  Fixed solar consant (disabled with -1):             ',f10.4/, &
             '  Enable temperature warnings:                        ',l5/ )

   end subroutine radiation_readnl


   !================================================================================================

   subroutine radiation_register()

      !----------------------------------------------------------------------------
      ! Register radiation fields in the physics buffer
      !----------------------------------------------------------------------------

      use physics_buffer, only: pbuf_add_field, dtype_r8

      integer :: idx  ! dummy index for adding fields to physics buffer

      ! Heating rate profiles; QRS is the shortwave radiative heating rate, and QRL
      ! is the longwave radiative heating rate
      ! TODO: Do QRS and QRL need to be set to "global" scope on the physics
      ! buffer? Doing so forces them to be written to restart files, but is this
      ! needed? It does not look like their values get read anywhere in this
      ! module, just overwritten, so they may not need to be written to restarts.
      call pbuf_add_field('QRS', 'global', dtype_r8, (/pcols,pver/), idx)
      call pbuf_add_field('QRL', 'global', dtype_r8, (/pcols,pver/), idx)

      ! If the namelist has been configured for preserving the spectral fluxes, then create
      ! physics buffer variables to store the results. These are fluxes per
      ! spectral band, as follows:
      !     SU: shortwave up
      !     SD: shortwave down
      !     LU: longwave up
      !     LD: longwave down
      ! TODO: Do these need to be "global"? Why not "physpkg"?
      if (spectralflux) then
         call pbuf_add_field('SU', 'global', dtype_r8, (/pcols,pverp,nswbands/), idx)
         call pbuf_add_field('SD', 'global', dtype_r8, (/pcols,pverp,nswbands/), idx)
         call pbuf_add_field('LU', 'global', dtype_r8, (/pcols,pverp,nlwbands/), idx)
         call pbuf_add_field('LD', 'global', dtype_r8, (/pcols,pverp,nlwbands/), idx)
      end if

   end subroutine radiation_register

   !===============================================================================

   function radiation_do(op, timestep)

      !---------------------------------------------------------------------------- 
      ! Purpose: Returns true if the specified operation is done this timestep.
      !
      ! History: Copied from RRTMG implementation.
      !----------------------------------------------------------------------------

      use time_manager, only: get_nstep

      ! Arguments
      character(len=*), intent(in)  :: op             ! name of operation
      integer, intent(in), optional :: timestep       ! model timestep
      logical                       :: radiation_do   ! return value

      ! Local variables
      integer :: nstep             ! current timestep number

      ! If passed a timestep explicitly, use that. Otherwise, get timestep from
      ! time_manager routine.
      if (present(timestep)) then
         nstep = timestep
      else
         nstep = get_nstep()
      end if

      ! Figure out whether or not to do radiation this timestep. In general, we
      ! allow different frequencies of call the radiation code for the shortwave
      ! and the longwave. In practice, these are probably usually the same.
      select case (op)
         case ('sw') ! do a shortwave heating calc this timestep?
            if (iradsw==0) then
               radiation_do = .false.
            else
               radiation_do = nstep == 0 .or. iradsw == 1                     &
                             .or. (mod(nstep-1,iradsw) == 0 .and. nstep /= 1) &
                             .or. nstep <= irad_always
            end if
         case ('lw') ! do a longwave heating calc this timestep?
            if (iradlw==0) then
               radiation_do = .false.
            else
               radiation_do = nstep == 0 .or. iradlw == 1                     &
                             .or. (mod(nstep-1,iradlw) == 0 .and. nstep /= 1) &
                             .or. nstep <= irad_always
            end if
         case default
            call endrun('radiation_do: unknown operation:'//op)
      end select
   end function radiation_do

   !================================================================================================

   real(r8) function radiation_nextsw_cday()
        
      !----------------------------------------------------------------------- 
      ! Purpose: Returns calendar day of next sw radiation calculation
      !
      ! History: Copied from RRTMG implementation.
      !-----------------------------------------------------------------------

      use time_manager, only: get_curr_calday, get_nstep, get_step_size

      ! Local variables
      integer :: nstep      ! timestep counter
      logical :: dosw       ! true => do shosrtwave calc   
      integer :: offset     ! offset for calendar day calculation
      integer :: dTime      ! integer timestep size
      real(r8):: calday     ! calendar day of 
      !-----------------------------------------------------------------------

      radiation_nextsw_cday = -1._r8
      dosw   = .false.
      nstep  = get_nstep()
      dtime  = get_step_size()
      offset = 0
      if (iradsw/=0) then
         do while (.not. dosw)
            nstep = nstep + 1
            offset = offset + dtime
            if (radiation_do('sw', nstep)) then
               radiation_nextsw_cday = get_curr_calday(offset=offset) 
               dosw = .true.
            end if
         end do
      end if
  
   end function radiation_nextsw_cday

   !================================================================================================

   subroutine radiation_init(state)
   !-------------------------------------------------------------------------------
   ! Purpose: Initialize the radiation parameterization and add fields to the 
   ! history buffer
   ! 
   ! History: Copied from RRTMG implemenation.
   !-------------------------------------------------------------------------------
      use physics_buffer,     only: pbuf_get_index
      use cam_history,        only: addfld, horiz_only, add_default
      use constituents,       only: cnst_get_ind
      use phys_control,       only: phys_getopts
      use rad_constituents,   only: N_DIAG, rad_cnst_get_call_list, rad_cnst_get_info
      use cospsimulator_intr, only: docosp, cospsimulator_intr_init
      use hirsbt,             only: hirsbt_init
      use hirsbtpar,          only: hirsname, msuname
      use modal_aer_opt,      only: modal_aer_opt_init
      use time_manager,       only: get_nstep, get_step_size, is_first_restart_step
      use radiation_data,     only: init_rad_data
      use physics_types,      only: physics_state

      ! For optics
      use cloud_rad_props, only: cloud_rad_props_init
      use ebert_curry, only: ec_rad_props_init
      use slingo, only: slingo_rad_props_init

      ! Physics state is going to be needed for perturbation growth tests, but we
      ! have yet to implement this in RRTMGP. It is a vector because at the point
      ! where this subroutine is called, physics_state has not subset for a
      ! specific chunk (i.e., state is a vector of all chunks, indexed by lchnk)
      type(physics_state), intent(in) :: state(:)

      integer :: icall, nmodes
      logical :: active_calls(0:N_DIAG)
      integer :: nstep                       ! current timestep number
      logical :: history_amwg                ! output the variables used by the AMWG diag package
      logical :: history_vdiag               ! output the variables used by the AMWG variability diag package
      logical :: history_budget              ! output tendencies and state variables for CAM4
                                             ! temperature, water vapor, cloud ice and cloud
                                             ! liquid budgets.
      integer :: history_budget_histfile_num ! output history file number for budget fields
      integer :: err
      integer :: dtime  ! time step
      integer :: cldfsnow_idx = 0 

      character(len=128) :: error_message

      character(len=32) :: subname = 'radiation_init'

      character(len=10), dimension(3) :: dims_crm_rad = (/'crm_nx_rad','crm_ny_rad','crm_nz    '/)

      !-----------------------------------------------------------------------

      ! Initialize cloud optics
      call cloud_rad_props_init()
      call ec_rad_props_init()
      call slingo_rad_props_init()

      ! Initialize output fields for offline driver.
      ! TODO: do we need to keep this functionality? Where is the offline driver?
      ! Do we need to write a new offline driver for RRTMGP?
      call init_rad_data()

      ! Do initialization for perturbation growth tests
      call perturbation_growth_init()

      ! Setup the RRTMGP interface
      call rrtmgp_initialize(size(active_gases), active_gases, rrtmgp_coefficients_file_sw, rrtmgp_coefficients_file_lw)

      ! Set number of levels used in radiation calculations
#ifdef NO_EXTRA_RAD_LEVEL
      nlev_rad = pver
#else
      nlev_rad = pver + 1
#endif

      ! Indices on radiation grid that correspond to top and bottom of the model
      ! grid...this is for cases where we want to add an extra layer above the
      ! model top to deal with radiative heating above the model top.
      ktop = nlev_rad - pver + 1
      kbot = nlev_rad

      ! Set the radiation timestep for cosz calculations if requested using the 
      ! adjusted iradsw value from radiation
      if (use_rad_dt_cosz)  then
         dtime  = get_step_size()
         dt_avg = iradsw*dtime
      end if

      ! Get options for history file outputs
      call phys_getopts(history_amwg_out = history_amwg,    &
                        history_vdiag_out = history_vdiag,   &
                        history_budget_out = history_budget,  &
                        history_budget_histfile_num_out = history_budget_histfile_num)
      
      ! Determine whether modal aerosols are affecting the climate, and if so
      ! then initialize the modal aerosol optics module
      call rad_cnst_get_info(0, nmodes=nmodes)
      if (nmodes > 0) call modal_aer_opt_init()

      ! TODO: what is this?
      call hirsbt_init()

      ! "irad_always" is number of time steps to execute radiation continuously from start of
      ! initial OR restart run
      nstep = get_nstep()
      if (irad_always > 0) then
         irad_always = irad_always + nstep
      end if

      ! Initialize the satellite simulator package (COSP). 
      ! TODO: Should this be moved to a higher level? 
      ! This should probably not be specific to a given radiation
      ! package. Also move most of this to cospsimulator package to handle itself,
      ! rather than relying on radiation driver handling this logic. Too much
      ! duplicate code.
      if (docosp) then
         ! Initialization for the simulator package.
         call cospsimulator_intr_init
      
         ! Allocate and initialize the counter that is used to determine when to
         ! call the simulator package.
         allocate(cosp_cnt(begchunk:endchunk))
         if (is_first_restart_step()) then
            cosp_cnt(begchunk:endchunk) = cosp_cnt_init
         else
            cosp_cnt(begchunk:endchunk) = 0     
         end if
      end if

      !
      ! Add fields to history buffer
      !
      
      ! Register new dimensions
      call get_sw_spectral_midpoints(sw_band_midpoints, 'cm-1')
      call get_lw_spectral_midpoints(lw_band_midpoints, 'cm-1')
      call add_hist_coord('swband', nswbands, 'Shortwave wavenumber', 'cm-1', sw_band_midpoints)
      call add_hist_coord('lwband', nlwbands, 'Longwave wavenumber', 'cm-1', lw_band_midpoints)

      ! Shortwave radiation
      call addfld('TOT_CLD_VISTAU', (/ 'lev' /), 'A',   '1', &
                  'Total gridbox cloud visible optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('TOT_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', &
                  'Total in-cloud visible optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('LIQ_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', &
                  'Liquid in-cloud visible optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('ICE_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', &
                  'Ice in-cloud visible optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)

      call add_default('TOT_CLD_VISTAU',  1, ' ')
      call add_default('TOT_ICLD_VISTAU', 1, ' ')

      ! Band-by-band cloud optical properties
      call addfld('CLOUD_TAU_SW', (/'lev   ','swband'/), 'I', '1', &
                  'Cloud shortwave extinction optical depth', sampling_seq='rad_lwsw', flag_xyfill=.true.) 
      call addfld('CLOUD_SSA_SW', (/'lev   ','swband'/), 'I', '1', &
                  'Cloud shortwave single scattering albedo', sampling_seq='rad_lwsw', flag_xyfill=.true.) 
      call addfld('CLOUD_G_SW', (/'lev   ','swband'/), 'I', '1', &
                  'Cloud shortwave assymmetry parameter', sampling_seq='rad_lwsw', flag_xyfill=.true.) 
      call addfld('CLOUD_TAU_LW', (/'lev   ','lwband'/), 'I', '1', &
                  'Cloud longwave absorption optical depth', sampling_seq='rad_lwsw', flag_xyfill=.true.) 

      ! Band-by-band shortwave albedos
      call addfld('SW_ALBEDO_DIR', (/'swband'/), 'I', '1', &
                  'Shortwave direct-beam albedo', &
                  flag_xyfill=.true., sampling_seq='rad_lwsw') 
      call addfld('SW_ALBEDO_DIF', (/'swband'/), 'I', '1', &
                  'Shortwave diffuse-beam albedo', &
                  flag_xyfill=.true., sampling_seq='rad_lwsw') 

      ! get list of active radiation calls
      call rad_cnst_get_call_list(active_calls)

      ! TODO: what is active_calls and diag?
      do icall = 0,N_DIAG
         if (active_calls(icall)) then
            call addfld('SOLIN'//diag(icall), horiz_only, 'A',   'W/m2', &
                        'Solar insolation', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('SOLL'//diag(icall),  horiz_only, 'A',   'W/m2', &
                        'Solar downward near infrared direct to surface', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('SOLS'//diag(icall),  horiz_only, 'A',   'W/m2', &
                        'Solar downward visible direct to surface', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('SOLLD'//diag(icall), horiz_only, 'A',   'W/m2', &
                        'Solar downward near infrared diffuse to surface', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('SOLSD'//diag(icall), horiz_only, 'A',   'W/m2', &
                        'Solar downward visible diffuse to surface', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('QRS'//diag(icall),   (/ 'lev' /), 'A', 'K/s', &
                        'Solar heating rate', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('QRSC'//diag(icall),  (/ 'lev' /), 'A', 'K/s', &
                        'Clearsky solar heating rate', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSNS'//diag(icall),   horiz_only,  'A', 'W/m2', &
                        'Net solar flux at surface', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSNT'//diag(icall),   horiz_only,  'A', 'W/m2', &
                        'Net solar flux at top of model', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSNTOA'//diag(icall), horiz_only,  'A', 'W/m2', &
                        'Net solar flux at top of atmosphere', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSUTOA'//diag(icall), horiz_only,  'A',  'W/m2', &
                        'Upwelling solar flux at top of atmosphere', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSUT'//diag(icall), horiz_only,  'A',  'W/m2', &
                        'Upwelling solar flux at top of model', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSUTC'//diag(icall), horiz_only,  'A',  'W/m2', &
                        'Clearsky upwelling solar flux at top of model', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSNTOAC'//diag(icall), horiz_only, 'A',  'W/m2', &
                        'Clearsky net solar flux at top of atmosphere', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSUTOAC'//diag(icall), horiz_only, 'A',  'W/m2', &
                        'Clearsky upwelling solar flux at top of atmosphere', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSN200'//diag(icall), horiz_only,  'A',  'W/m2', &
                        'Net shortwave flux at 200 mb', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSN200C'//diag(icall), horiz_only, 'A',  'W/m2', &
                        'Clearsky net shortwave flux at 200 mb', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSNTC'//diag(icall), horiz_only,   'A',  'W/m2', &
                        'Clearsky net solar flux at top of model', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSNSC'//diag(icall), horiz_only,   'A',  'W/m2', &
                        'Clearsky net solar flux at surface', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSDSC'//diag(icall), horiz_only,   'A',  'W/m2', &
                        'Clearsky downwelling solar flux at surface', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSDS'//diag(icall), horiz_only,    'A',  'W/m2', &
                        'Downwelling solar flux at surface', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FUS'//diag(icall),  (/ 'ilev' /),  'A',  'W/m2', &
                        'Shortwave upward flux', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FDS'//diag(icall),  (/ 'ilev' /),  'A',  'W/m2', &
                        'Shortwave downward flux', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FDS_DIR'//diag(icall),  (/ 'ilev' /),  'A',  'W/m2', &
                        'Shortwave direct-beam downward flux', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FNS'//diag(icall),  (/ 'ilev' /),  'A',  'W/m2', &
                        'Shortwave net flux', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FUSC'//diag(icall),  (/ 'ilev' /), 'A',  'W/m2', &
                        'Shortwave clear-sky upward flux', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FDSC'//diag(icall),  (/ 'ilev' /), 'A',  'W/m2', &
                        'Shortwave clear-sky downward flux', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FDSC_DIR'//diag(icall),  (/ 'ilev' /), 'A',  'W/m2', &
                        'Shortwave clear-sky direct-beam downward flux', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FNSC'//diag(icall),  (/ 'ilev' /), 'A',  'W/m2', &
                        'Shortwave clear-sky net flux', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSNIRTOA'//diag(icall), horiz_only, 'A', 'W/m2', &
                        'Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSNRTOAC'//diag(icall), horiz_only, 'A', 'W/m2', &
                        'Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FSNRTOAS'//diag(icall), horiz_only, 'A', 'W/m2', &
                        'Net near-infrared flux (>= 0.7 microns) at top of atmosphere', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('SWCF'//diag(icall),     horiz_only, 'A', 'W/m2', &
                        'Shortwave cloud forcing', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)

            if (history_amwg) then
               call add_default('SOLIN'//diag(icall),   1, ' ')
               call add_default('QRS'//diag(icall),     1, ' ')
               call add_default('FSNS'//diag(icall),    1, ' ')
               call add_default('FSNT'//diag(icall),    1, ' ')
               call add_default('FSNTOA'//diag(icall),  1, ' ')
               call add_default('FSUTOA'//diag(icall),  1, ' ')
               call add_default('FSNTOAC'//diag(icall), 1, ' ')
               call add_default('FSUTOAC'//diag(icall), 1, ' ')
               call add_default('FSNTC'//diag(icall),   1, ' ')
               call add_default('FSNSC'//diag(icall),   1, ' ')
               call add_default('FSDSC'//diag(icall),   1, ' ')
               call add_default('FSDS'//diag(icall),    1, ' ')
               call add_default('SWCF'//diag(icall),    1, ' ')
            end if  ! history_amwg
         end if  ! active_calls(icall)
      end do  ! icall = 0,N_DIAG

      ! Cosine of solar zenith angle (primarily for debugging)
      call addfld('COSZRS', horiz_only, 'A', 'None', &
                  'Cosine of solar zenith angle', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)

      if (single_column .and. scm_crm_mode) then
         call add_default ('FUS     ', 1, ' ')
         call add_default ('FUSC    ', 1, ' ')
         call add_default ('FDS     ', 1, ' ')
         call add_default ('FDSC    ', 1, ' ')
      end if

      ! Longwave radiation
      do icall = 0,N_DIAG
         if (active_calls(icall)) then
            call addfld('QRL'//diag(icall),     (/'lev'/),  'A', 'K/s',  &
                        'Longwave heating rate',                         &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('QRLC'//diag(icall),    (/'lev'/),  'A', 'K/s',  &
                        'Clearsky longwave heating rate',                &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FLDS'//diag(icall),    horiz_only, 'A', 'W/m2', &
                        'Downwelling longwave flux at surface',          &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FLDSC'//diag(icall),   horiz_only, 'A', 'W/m2', &
                        'Clearsky Downwelling longwave flux at surface', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FLNS'//diag(icall),    horiz_only, 'A', 'W/m2', &
                        'Net longwave flux at surface',                  &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FLNT'//diag(icall),    horiz_only, 'A', 'W/m2', &
                        'Net longwave flux at top of model',             &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FLUT'//diag(icall),    horiz_only, 'A', 'W/m2', &
                        'Upwelling longwave flux at top of model',       &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FLUTC'//diag(icall),   horiz_only, 'A', 'W/m2', &
                        'Clearsky upwelling longwave flux at top of model', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FLNTC'//diag(icall),   horiz_only, 'A', 'W/m2', &
                        'Clearsky net longwave flux at top of model',    &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('LWCF'//diag(icall),    horiz_only, 'A', 'W/m2', &
                        'Longwave cloud forcing',                        &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FLN200'//diag(icall),  horiz_only, 'A', 'W/m2', &
                        'Net longwave flux at 200 mb',                   &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FLN200C'//diag(icall), horiz_only, 'A', 'W/m2', &
                        'Clearsky net longwave flux at 200 mb',          &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FLNSC'//diag(icall),   horiz_only, 'A', 'W/m2', &
                        'Clearsky net longwave flux at surface',         &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FUL'//diag(icall),     (/'ilev'/), 'A', 'W/m2', &
                        'Longwave upward flux', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FDL'//diag(icall),     (/'ilev'/), 'A', 'W/m2', &
                        'Longwave downward flux', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FNL'//diag(icall),     (/'ilev'/), 'A', 'W/m2', &
                        'Longwave net flux', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FULC'//diag(icall),    (/'ilev'/), 'A', 'W/m2', &
                        'Longwave clear-sky upward flux', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FDLC'//diag(icall),    (/'ilev'/), 'A', 'W/m2', &
                        'Longwave clear-sky downward flux', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)
            call addfld('FNLC'//diag(icall),    (/'ilev'/), 'A', 'W/m2', &
                        'Longwave clear-sky net flux', &
                        sampling_seq='rad_lwsw', flag_xyfill=.true.)

            call add_default('QRL'//diag(icall),   1, ' ')
            if (history_amwg) then
               call add_default('QRL'//diag(icall),   1, ' ')
               call add_default('FLNS'//diag(icall),  1, ' ')
               call add_default('FLDS'//diag(icall),  1, ' ')
               call add_default('FLNT'//diag(icall),  1, ' ')
               call add_default('FLUT'//diag(icall),  1, ' ')
               call add_default('FLUTC'//diag(icall), 1, ' ')
               call add_default('FLNTC'//diag(icall), 1, ' ')
               call add_default('FLNSC'//diag(icall), 1, ' ')
               call add_default('LWCF'//diag(icall),  1, ' ')
            endif  ! history_amwg
         end if
      end do

      ! Add cloud-scale radiative quantities
      call addfld('CRM_QRAD', dims_crm_rad, 'A', 'K/s', &
                  'Radiative heating tendency', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('CRM_QRS ', dims_crm_rad, 'I', 'K/s', &
                  'CRM Shortwave radiative heating rate', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('CRM_QRSC', dims_crm_rad, 'I', 'K/s', &
                  'CRM clear-sky shortwave radiative heating rate', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('CRM_QRL ', dims_crm_rad, 'I', 'K/s', &
                  'CRM Longwave radiative heating rate', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('CRM_QRLC', dims_crm_rad, 'I', 'K/s', &
                  'CRM clear-sky longwave radiative heating rate', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('CRM_CLD_RAD', dims_crm_rad, 'I', 'fraction', &
                  'CRM cloud fraction', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)

      call addfld('EMIS', (/ 'lev' /), 'A', '1', 'Cloud longwave emissivity')

      ! Add default fields for single column mode
      if (single_column.and.scm_crm_mode) then
         call add_default ('FUL     ', 1, ' ')
         call add_default ('FULC    ', 1, ' ')
         call add_default ('FDL     ', 1, ' ')
         call add_default ('FDLC    ', 1, ' ')
      endif

      ! HIRS/MSU diagnostic brightness temperatures
      if (dohirs) then
         call addfld(hirsname(1),horiz_only,'A','K','HIRS CH2 infra-red brightness temperature', &
                     sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld(hirsname(2),horiz_only,'A','K','HIRS CH4 infra-red brightness temperature', &
                     sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld(hirsname(3),horiz_only,'A','K','HIRS CH6 infra-red brightness temperature', &
                     sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld(hirsname(4),horiz_only,'A','K','HIRS CH8 infra-red brightness temperature', &
                     sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld(hirsname(5),horiz_only,'A','K','HIRS CH10 infra-red brightness temperature', &
                     sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld(hirsname(6),horiz_only,'A','K','HIRS CH11 infra-red brightness temperature', &
                     sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld(hirsname(7),horiz_only,'A','K','HIRS CH12 infra-red brightness temperature', &
                     sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld(msuname(1),horiz_only,'A','K','MSU CH1 microwave brightness temperature', &
                     sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld(msuname(2),horiz_only,'A','K','MSU CH2 microwave brightness temperature', &
                     sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld(msuname(3),horiz_only,'A','K','MSU CH3 microwave brightness temperature', &
                     sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld(msuname(4),horiz_only,'A','K','MSU CH4 microwave brightness temperature', &
                     sampling_seq='rad_lwsw', flag_xyfill=.true.)

         ! Add HIRS/MSU fields to default history files
         call add_default (hirsname(1), 1, ' ')
         call add_default (hirsname(2), 1, ' ')
         call add_default (hirsname(3), 1, ' ')
         call add_default (hirsname(4), 1, ' ')
         call add_default (hirsname(5), 1, ' ')
         call add_default (hirsname(6), 1, ' ')
         call add_default (hirsname(7), 1, ' ')
         call add_default (msuname(1), 1, ' ')
         call add_default (msuname(2), 1, ' ')
         call add_default (msuname(3), 1, ' ')
         call add_default (msuname(4), 1, ' ')
      end if

      ! Heating rate needed for d(theta)/dt computation
      ! TODO: why include this here? Can't this be calculated from QRS and QRL,
      ! which are already output?
      call addfld('HR',(/ 'lev' /), 'A','K/s', &
                  'Heating rate needed for d(theta)/dt computation', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)

      if (history_budget .and. history_budget_histfile_num > 1) then
         call add_default ('QRL     ', history_budget_histfile_num, ' ')
         call add_default ('QRS     ', history_budget_histfile_num, ' ')
      end if

      if (history_vdiag) then
         call add_default('FLUT', 2, ' ')
         call add_default('FLUT', 3, ' ')
      end if

      cldfsnow_idx = 0
      cldfsnow_idx = pbuf_get_index('CLDFSNOW',errcode=err)
      if (cldfsnow_idx > 0) then
         call addfld('CLDFSNOW',(/ 'lev' /),'I','1','CLDFSNOW',flag_xyfill=.true.)
         call addfld('SNOW_ICLD_VISTAU', (/ 'lev' /), 'A', '1', &
                     'Snow in-cloud extinction visible sw optical depth', &
                     sampling_seq='rad_lwsw', flag_xyfill=.true.)
      endif

   end subroutine radiation_init

   subroutine radiation_final()
      call rrtmgp_finalize()
   end subroutine radiation_final

   subroutine perturbation_growth_init()

      use cam_logfile, only: iulog

      character(len=32) :: subname = 'perturbation_growth_init'

      ! Wrap modes in an ifdef for now since this is not implemented here, and I
      ! have some work to do to figure out what the heck this is meant to do.
#ifdef DO_PERGRO_MODS
      ! Modification needed by pergro_mods for generating random numbers
      if (pergro_mods) then
         max_chnks_in_blk = maxval(npchunks(:))  !maximum of the number for chunks in each procs
         allocate(clm_rand_seed(pcols,kiss_seed_num,max_chnks_in_blk), stat=astat)
         if (astat /= 0) then
            write(iulog,*) 'radiation.F90(rrtmg)-radiation_init: failed to allocate clm_rand_seed; error = ',astat
            call endrun
         end if

         allocate(tot_chnk_till_this_prc(0:npes-1), stat=astat)
         if (astat /= 0) then
            write(errstr,*) 'radiation.F90(rrtmg)-radiation_init: failed to allocate tot_chnk_till_this_prc variable; error = ',astat
            call endrun(errstr)
         end if
          
         !BSINGH - Build lat lon relationship to chunk and column
         !Compute maximum number of chunks each processor have
         if (masterproc) then
             tot_chnk_till_this_prc(0:npes-1) = huge(1)
             do ipes = 0, npes - 1
                tot_chnk_till_this_prc(ipes) = 0
                do ipes_tmp = 0, ipes-1
                   tot_chnk_till_this_prc(ipes) = tot_chnk_till_this_prc(ipes) + npchunks(ipes_tmp)
                enddo
             enddo
          endif
#ifdef SPMD
          !BSINGH - Ideally we should use mpi_scatter but we are using this variable
          !in "if(masterproc)" below in phys_run1, so broadcast is iused here
          call mpibcast(tot_chnk_till_this_prc,npes, mpi_integer, 0, mpicom)
#endif
          call get_block_bounds_d(firstblock,lastblock)
          
          allocate(clm_id(pcols,max_chnks_in_blk), stat=astat)
          if( astat /= 0 ) then
             write(errstr,*) 'radiation.F90(rrtmg)-radiation_init: failed to allocate clm_id; error = ',astat
             call endrun(errstr)
          end if
          
          allocate(clm_id_mstr(pcols,max_chnks_in_blk,npes), stat=astat)
          if( astat /= 0 ) then
             write(errstr,*) 'radiation.F90(rrtmg)-radiation_init: failed to allocate clm_id_mstr; error = ',astat
             call endrun(errstr)
          end if
          !compute all clm ids on masterproc and then scatter it ....
          if(masterproc) then
             do igcol = 1, ngcols_p
                imap = latlon_to_dyn_gcol_map(igcol)
                chunkid  = knuhcs(imap)%chunkid
                icol = knuhcs(imap)%col
                iown  = chunks(chunkid)%owner
                ilchnk = (chunks(chunkid)%lcid - lastblock) - tot_chnk_till_this_prc(iown)
                clm_id_mstr(icol,ilchnk,iown+1) = igcol
             enddo
          endif
          
#ifdef SPMD
          !Scatter
          tot_cols = pcols*max_chnks_in_blk
          call MPI_Scatter( clm_id_mstr, tot_cols,  mpi_integer, &
               clm_id,    tot_cols,  mpi_integer, 0,             &
               MPI_COMM_WORLD,ierr)
#else
          !BSINGH - Haven't tested it.....               
          call endrun('radiation.F90(rrtmg)-radiation_init: non-mpi compiles are not tested yet for pergro test...')
#endif       
       endif
          
       if (is_first_restart_step()) then
          cosp_cnt(begchunk:endchunk)=cosp_cnt_init
          if (pergro_mods) then
             !--------------------------------------
             !Read seeds from restart file
             !--------------------------------------
             !For restart runs, rad_randn_seedrst array  will already be allocated in the restart_physics.F90
             
             do ilchnk = 1, max_chnks_in_blk
                lchnk = begchunk + (ilchnk -1)
                ncol = phys_state(lchnk)%ncol
                do iseed = 1, kiss_seed_num
                   do icol = 1, ncol                
                      clm_rand_seed(icol,iseed,ilchnk) = rad_randn_seedrst(icol,iseed,lchnk)
                   enddo
                enddo
             enddo
          endif
       else
          cosp_cnt(begchunk:endchunk)=0           
          if (pergro_mods) then
             !---------------------------------------
             !create seeds based off of column ids
             !---------------------------------------
             !allocate array rad_randn_seedrst for initial run for  maintaining exact restarts
             !For restart runs, it will already be allocated in the restart_physics.F90
             allocate(rad_randn_seedrst(pcols,kiss_seed_num,begchunk:endchunk), stat=astat)
             if( astat /= 0 ) then
                write(iulog,*) 'radiation.F90(rrtmg)-radiation_init: failed to allocate rad_randn_seedrst; error = ',astat
                call endrun
             end if
             do ilchnk = 1, max_chnks_in_blk
                lchnk = begchunk + (ilchnk -1)
                ncol = phys_state(lchnk)%ncol
                do iseed = 1, kiss_seed_num
                   do icol = 1, ncol
                      id = clm_id(icol,ilchnk)
                      clm_rand_seed(icol,iseed,ilchnk) = id + (iseed -1)
                   enddo
                enddo
             enddo
          endif
       end if

#endif /* DO_PERGRO_MODS */

   end subroutine perturbation_growth_init


   !===============================================================================
        
   !----------------------------------------------------------------------------
   ! Main driver for radiation computation. Calls the shortwave and longwave
   ! drivers, and calculates the radiative heating from the resulting fluxes.
   ! Primary output from this routine is the heating tendency due to radiative
   ! transfer, as a ptend object.
   !
   ! Notes for MMF:
   ! We are currently using a separate driver when compiling the code using the
   ! MMF configuration, and right now this code is quite a bit more complicated
   ! than it needs to be. The plan with this is to factor pbuf and state objects
   ! out of the lower-level routines so that we can call most of these routines
   ! with arbitrary numbers of columns, and then pack the MMF data into arrays
   ! dimensioned (ncol*crm_nx_rad*crm_ny_rad, nlev, ...). The majority of this
   ! routine could then be run as-is, with just a "packing" routine at the top,
   ! and then a domain averaging routine at the end (to get the domain-averaged
   ! fluxes the non-CRM code expects, like surface fluxes to exchange with the
   ! surface models). This will expose a lot more parallelism and make the code
   ! a lot more GPU-friendly, but we still need to refactor the aerosol optics
   ! code to accomplish this, which is a big task.
   subroutine radiation_tend(state_in,ptend,    pbuf,          cam_out, cam_in,  &
                             landfrac,icefrac,  snowh,                           &
                             fsns,    fsnt,     flns,          flnt,             &
                             fsds,    net_flux, is_cmip6_volc, dt,               &
                             clear_rh                                            )

      ! Performance module needed for timing functions
      use perf_mod, only: t_startf, t_stopf

      ! CAM derived types; needed for surface exchange fields, physics state, and
      ! tendency fields
      use camsrfexch, only: cam_out_t, cam_in_t
      use physics_types, only: physics_state, physics_ptend, physics_state_copy

      ! Utilities for interacting with the physics buffer
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
                                pbuf_get_index

      ! For calculating radiative heating tendencies
      use radheat, only: radheat_tend

      ! For getting radiative constituent gases
      use rad_constituents, only: N_DIAG, rad_cnst_get_call_list
      use constituents,     only: cnst_get_ind

      ! Index to visible channel for diagnostic outputs
      use radconstants, only: idx_sw_diag

      ! CAM history module provides subroutine to send output data to the history
      ! buffer to be aggregated and written to disk
      use cam_history, only: outfld

      use pkg_cldoptics, only: cldefr  ! for sam1mom microphysics

      use physconst,    only: gravit, cpair
      use phys_control, only: phys_getopts

      use radiation_state, only: set_rad_state
      use radiation_utils, only: calculate_heating_rate
      use cam_optics, only: set_aerosol_optics_lw, set_aerosol_optics_sw
      use mmf_optics, only: get_cloud_optics_sw, sample_cloud_optics_sw, &
                            get_cloud_optics_lw, sample_cloud_optics_lw

#ifdef MODAL_AERO
      use modal_aero_data, only: ntot_amode
#endif

      ! For running CFMIP Observation Simulator Package (COSP)
      use cospsimulator_intr, only: docosp, cospsimulator_intr_run, cosp_nradsteps

      ! ---------------------------------------------------------------------------
      ! Arguments
      ! ---------------------------------------------------------------------------

      ! Physics state variables
      type(physics_state), intent(in), target :: state_in

      ! Heating tendencies calculated in this subroutine
      type(physics_ptend), intent(out) :: ptend

      ! Fields from other parameterizations that persist across timesteps
      type(physics_buffer_desc), pointer :: pbuf(:)

      ! Surface fluxes?
      type(cam_out_t), intent(inout) :: cam_out
      type(cam_in_t), intent(in) :: cam_in

      ! Net flux calculated in this routine; used to check energy conservation in
      ! the physics package driver?
      real(r8), intent(inout) :: net_flux(pcols)

      ! This should be module data or something specific to aerosol where it is
      ! used?
      logical,  intent(in)    :: is_cmip6_volc    ! true if cmip6 style volcanic file is read otherwise false 

      real(r8),  intent(in)   :: dt               ! time step(s) - needed for aerosol optics call

      real(r8), optional,  intent(in)    :: clear_rh(pcols,pver) ! optional clear air relative humidity
                                                                 ! that gets passed to modal_aero_wateruptake_dr

      ! These are not used anymore and exist only because the radiation call is
      ! inflexible
      real(r8), intent(in)    :: landfrac(pcols)  ! land fraction
      real(r8), intent(in)    :: icefrac(pcols)   ! land fraction
      real(r8), intent(in)    :: snowh(pcols)     ! Snow depth (liquid water equivalent)

      ! Surface and top fluxes
      real(r8), intent(inout) :: fsns(pcols)      ! Surface solar absorbed flux
      real(r8), intent(inout) :: fsnt(pcols)      ! Net column abs solar flux at model top
      real(r8), intent(inout) :: flns(pcols)      ! Srf longwave cooling (up-down) flux
      real(r8), intent(inout) :: flnt(pcols)      ! Net outgoing lw flux at model top
      real(r8), intent(inout) :: fsds(pcols)      ! Surface solar down flux

      ! ---------------------------------------------------------------------------
      ! Local variables
      ! ---------------------------------------------------------------------------

      ! Pointers to heating rates on physics buffer
      real(r8), pointer :: qrs(:,:)  ! shortwave radiative heating rate 
      real(r8), pointer :: qrl(:,:)  ! longwave  radiative heating rate 

      ! Clear-sky heating rates are not on the physics buffer, and we have no
      ! reason to put them there, so declare these as regular arrays here
      real(r8) :: qrsc(pcols,pver)
      real(r8) :: qrlc(pcols,pver)

      ! Temporary variable for heating rate output
      real(r8) :: hr(pcols,pver)

      ! Number of columns; ncol is number of GCM columns in current chunk,
      ! ncol_tot is product of GCM and CRM columns (for packing GCM and CRM
      ! columns into a single dimension)
      integer :: ncol, ncol_tot

      ! Loop indices
      integer :: icall, ix, iy, iz, ilev, j, iband, igpt

      ! Everyone needs a name
      character(*), parameter :: subname = 'radiation_tend'

      ! Radiative fluxes
      type(fluxes_t) :: fluxes_allsky, fluxes_clrsky, fluxes_allsky_packed, fluxes_clrsky_packed, fluxes_allsky_day, fluxes_clrsky_day

      ! Copy of state
      type(physics_state) :: state

      ! For loops over diagnostic calls (TODO: what does this mean?)
      logical :: active_calls(0:N_DIAG)

      ! State fields that are passed into RRTMGP. Some of these may need to
      ! modified from what exists in the physics_state object, i.e. to clip
      ! temperatures to make sure they are within the valid range.
      real(r8), dimension(pcols*crm_nx_rad*crm_ny_rad,nlev_rad) :: tmid_packed, pmid_packed
      real(r8), dimension(pcols*crm_nx_rad*crm_ny_rad,nlev_rad+1) :: tint_packed, pint_packed
      real(r8), dimension(pcols,nlev_rad) :: tmid_col, pmid_col
      real(r8), dimension(pcols,nlev_rad+1) :: tint_col, pint_col

      ! Daytime subsets
      real(r8), dimension(pcols*crm_nx_rad*crm_ny_rad,nlev_rad) :: tmid_day, pmid_day
      real(r8), dimension(pcols*crm_nx_rad*crm_ny_rad,nlev_rad+1) :: tint_day, pint_day
      real(r8), dimension(nswbands,pcols*crm_nx_rad*crm_ny_rad) :: albedo_dir_day, albedo_dif_day

      ! Surface emissivity needed for longwave
      real(r8) :: surface_emissivity(nlwbands,pcols*crm_nx_rad*crm_ny_rad)

      ! Temporary heating rates on radiation vertical grid
      real(r8), dimension(pcols*crm_nx_rad*crm_ny_rad,pver) :: qrl_packed, qrlc_packed

      ! Packed optics for COSP
      real(r8), dimension(pcols*crm_nx_rad*crm_ny_rad,pver) :: dtau_packed, dems_packed

      ! Cloud properties for optics
      real(r8), dimension(pcols,pver) :: rel, rei, cld, iclwp, iciwp
      ! PACKED cloud physical properties for optics
      real(r8), dimension(pcols*crm_nx_rad*crm_ny_rad,pver) :: &
         rel_p, rei_p, cld_p, iclwp_p, iciwp_p, &
         rel_d, rei_d, cld_d, iclwp_d, iciwp_d

      ! Pointers to fields on physics buffer
      real(r8), pointer, dimension(:,:,:,:) :: crm_t, crm_qv, crm_qc, crm_qi, crm_cld, crm_qrad
      real(r8), pointer, dimension(:,:,:,:) :: crm_rel, crm_rei

      ! Indices to water tracers
      integer, parameter :: ixwatvap = 1

      ! Arrays to hold gas volume mixing ratios
      real(r8), dimension(size(active_gases),pcols,nlev_rad) :: vmr_col
      real(r8), dimension(size(active_gases),pcols*crm_nx_rad*crm_ny_rad,nlev_rad) :: vmr_packed, vmr_day

      ! CRM heating rates
      real(r8), dimension(pcols,crm_nx_rad,crm_ny_rad,crm_nz) :: crm_qrs, crm_qrsc, crm_qrl, crm_qrlc

      ! Albedo for shortwave calculations
      real(r8), dimension(nswbands,pcols) :: albedo_dir_col, albedo_dif_col
      real(r8), dimension(nswbands,pcols*crm_nx_rad*crm_ny_rad) :: albedo_dir_packed, albedo_dif_packed

      ! Cosine solar zenith angle
      real(r8) :: coszrs(pcols), coszrs_packed(pcols*crm_nx_rad*crm_ny_rad)
      real(r8) :: coszrs_day(pcols*crm_nx_rad*crm_ny_rad)

      ! Gathered indicies of day and night columns 
      ! chunk_column_index = day_indices(daylight_column_index)
      integer :: nday, nnight     ! Number of daylight columns
      integer :: day_indices(pcols), night_indices(pcols)   ! Indicies of daylight coumns
      integer, dimension(pcols*crm_nx_rad*crm_ny_rad) :: day_indices_packed, night_indices_packed

      ! Scaling factor for total sky irradiance; used to account for orbital
      ! eccentricity, and could be used to scale total sky irradiance for different
      ! climates as well (i.e., paleoclimate simulations)
      real(r8) :: tsi_scaling

      ! Temporary packed heating rates
      real(r8), dimension(pcols*crm_nx_rad*crm_ny_rad,pver) :: qrs_packed, qrsc_packed

      ! Working variables for optics
      real(r8), dimension(pcols,pver,nlwbands) :: aer_tau_bnd_lw
      real(r8), dimension(pcols*crm_nx_rad*crm_ny_rad,nlev_rad,nlwgpts) :: cld_tau_gpt_lw
      real(r8), dimension(pcols*crm_nx_rad*crm_ny_rad,nlev_rad,nlwbands) :: aer_tau_bnd_lw_packed, cld_tau_bnd_lw
      real(r8), dimension(pcols,pver,nswbands) :: &
         aer_tau_bnd_sw, aer_ssa_bnd_sw, aer_asm_bnd_sw
      real(r8), dimension(pcols*crm_nx_rad*crm_ny_rad,nlev_rad,nswbands) :: &
         cld_tau_bnd_day, cld_ssa_bnd_day, cld_asm_bnd_day
      real(r8), dimension(pcols*crm_nx_rad*crm_ny_rad,nlev_rad,nswgpts) :: &
         cld_tau_gpt_day, cld_ssa_gpt_day, cld_asm_gpt_day
      real(r8), dimension(pcols*crm_nx_rad*crm_ny_rad,nlev_rad,nswbands) :: &
         aer_tau_bnd_sw_packed, aer_ssa_bnd_sw_packed, aer_asm_bnd_sw_packed, &
         aer_tau_bnd_day, aer_ssa_bnd_day, aer_asm_bnd_day
      ! NOTE: these are diagnostic only (and only output in non-MMF version, but needed for optics calls)
      real(r8), dimension(pcols,pver,nswbands) :: liq_tau_bnd_sw, ice_tau_bnd_sw, snw_tau_bnd_sw
      real(r8), dimension(pcols,pver,nlwbands) :: liq_tau_bnd_lw, ice_tau_bnd_lw, snw_tau_bnd_lw

      integer, dimension(nswgpts) :: gpoint_bands_sw
      integer, dimension(nlwgpts) :: gpoint_bands_lw

      integer :: cosp_lwband  ! index to LW band used by COSP passive simulators
      integer :: cosp_swband  ! index to SW band used by COSP passive simulators

      ! Loop variables
      integer :: icol, ilay, iday

      !----------------------------------------------------------------------

      ! Copy state so we can use CAM routines with arrays replaced with data
      ! from CRM fields
      call physics_state_copy(state_in, state)

      ! Number of physics columns in this "chunk"
      ncol = state%ncol

      ! Total number of columns to operate on in one call to the radiative
      ! transfer solvers. When using superparameterization, this will be ncol
      ! times the total number of CRM columns so that we can pack all of the
      ! data together into one call.
      ncol_tot = ncol * crm_nx_rad * crm_ny_rad

      ! Set pointers to heating rates stored on physics buffer. These will be
      ! modified in this routine.
      call pbuf_get_field(pbuf, pbuf_get_index('QRS'), qrs)
      call pbuf_get_field(pbuf, pbuf_get_index('QRL'), qrl)

      ! Calculate derived quantities and pack data. This includes cloud
      ! optical properties. This is done separately from the blocks of
      ! code where shortwave and longwave fluxes are calculated so that we
      ! only have to pack/copy data once, rather than once for each of
      ! shortwave and longwave.
      if (radiation_do('sw') .or. radiation_do('lw')) then

         ! CRM fields
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_T_RAD'  ), crm_t   )
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_QV_RAD' ), crm_qv  )
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_QC_RAD' ), crm_qc  )
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_QI_RAD' ), crm_qi  )
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_CLD_RAD'), crm_cld )
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_REL')    , crm_rel )
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_REI')    , crm_rei )

         ! Output CRM cloud fraction
         call outfld('CRM_CLD_RAD', crm_cld(1:ncol,:,:,:), state%ncol, state%lchnk)

         ! Set surface emissivity to 1 here. There is a note in the RRTMG
         ! implementation that this is treated in the land model, but the old
         ! RRTMG implementation also sets this to 1. This probably does not make
         ! a lot of difference either way, but if a more intelligent value
         ! exists or is assumed in the model we should use it here as well.
         surface_emissivity = 1.0_r8

         ! Loop over "diagnostic calls"; these are additional configurations of
         ! gases used to calculate radiative effects for diagnostic purposes only,
         ! without affecting the climate. The only version that affects the
         ! simulation is the first (index 0).
         call rad_cnst_get_call_list(active_calls)
         do icall = N_DIAG,0,-1
            if (active_calls(icall)) then
               ! Calculate effective radius for optics used with single-moment microphysics and for COSP
               call cldefr(state%lchnk, ncol, state%t, rel, rei, state%ps, state%pmid, landfrac, icefrac, snowh)
               ! Populate CRM effective radii (used for COSP)
               do iz = 1,crm_nz
                  do iy = 1,crm_ny_rad
                     do ix = 1,crm_nx_rad
                        do icol = 1,ncol
                           ilay = pver - iz + 1
                           crm_rel(icol,ix,iy,iz) = rel(icol,ilay)
                           crm_rei(icol,ix,iy,iz) = rei(icol,ilay)
                        end do
                     end do
                  end do
               end do

               ! Get albedo. This uses CAM routines internally and just provides a
               ! wrapper to improve readability of the code here.
               call set_albedo(cam_in, albedo_dir_col(1:nswbands,1:ncol), albedo_dif_col(1:nswbands,1:ncol))

               ! Send albedos to history buffer (useful for debugging)
               call outfld('SW_ALBEDO_DIR', transpose(albedo_dir_col(1:nswbands,1:ncol)), ncol, state%lchnk)
               call outfld('SW_ALBEDO_DIF', transpose(albedo_dif_col(1:nswbands,1:ncol)), ncol, state%lchnk)

               ! Get cosine solar zenith angle for current time step, and send to
               ! history buffer
               call set_cosine_solar_zenith_angle(state, dt_avg, coszrs(1:ncol))
               call outfld('COSZRS', coszrs(1:ncol), ncol, state%lchnk)

               ! Gather night/day column indices for subsetting SW inputs; we only want to
               ! do the shortwave radiative transfer during the daytime to save
               ! computational cost (and because RRTMGP will fail for cosine solar zenith
               ! angles less than or equal to zero)
               day_indices_packed = 0
               call set_daynight_indices(coszrs(1:ncol), day_indices(1:ncol), night_indices(1:ncol))
               nnight = count(night_indices(1:ncol) > 0)

               ! Do aerosol optics; this was moved outside the CRM loop to 
               ! optimize performance. The impact was estimated to be negligible.
               aer_tau_bnd_sw = 0
               aer_ssa_bnd_sw = 0
               aer_asm_bnd_sw = 0
               aer_tau_bnd_lw = 0
               if (do_aerosol_rad) then
                  if (radiation_do('sw')) then
                     call t_startf('rad_aerosol_optics_sw')
                     call set_aerosol_optics_sw( &
                          icall, dt, state, pbuf, night_indices(1:nnight), is_cmip6_volc, &
                          aer_tau_bnd_sw, aer_ssa_bnd_sw, aer_asm_bnd_sw,  &
                          clear_rh=clear_rh)
                     ! Now reorder bands to be consistent with RRTMGP
                     ! TODO: fix the input files themselves!
                     do icol = 1,size(aer_tau_bnd_sw,1)
                        do ilay = 1,size(aer_tau_bnd_sw,2)
                           aer_tau_bnd_sw(icol,ilay,:) = reordered(aer_tau_bnd_sw(icol,ilay,:), rrtmg_to_rrtmgp_swbands)
                           aer_ssa_bnd_sw(icol,ilay,:) = reordered(aer_ssa_bnd_sw(icol,ilay,:), rrtmg_to_rrtmgp_swbands)
                           aer_asm_bnd_sw(icol,ilay,:) = reordered(aer_asm_bnd_sw(icol,ilay,:), rrtmg_to_rrtmgp_swbands)
                        end do
                     end do
                     call t_stopf('rad_aerosol_optics_sw')
                  end if ! radiation_do('sw')
                  if (radiation_do('lw')) then
                     call t_startf('rad_aerosol_optics_lw')
                     call set_aerosol_optics_lw(icall, dt, state, pbuf, is_cmip6_volc, aer_tau_bnd_lw)
                     call t_stopf('rad_aerosol_optics_lw')
                  end if ! radiation_do('lw')
               end if ! do_aerosol_rad

               aer_tau_bnd_lw_packed = 0._r8
               aer_tau_bnd_sw_packed = 0._r8
               aer_ssa_bnd_sw_packed = 0._r8
               aer_asm_bnd_sw_packed = 0._r8

               ! make sure water path variables are zeroed out
               iclwp = 0_r8
               iciwp = 0_r8
               cld   = 0_r8

               ! Populate fields with CRM data and pack into arrays dimensioned ncol_tot = ncol*crm_nx_rad*crm_ny_rad
               do iy = 1,crm_ny_rad
                  do ix = 1,crm_nx_rad
                     ! Fill GCM columns with CRM data
                     call t_startf('rad_set_state')
                     do iz = 1,crm_nz
                        ilev = pver - iz + 1
                        do icol = 1,ncol
                           state%q(icol,ilev,ixwatvap) = crm_qv(icol,ix,iy,iz)
                           state%t(icol,ilev)          = crm_t (icol,ix,iy,iz)
                           cld(icol,ilev) = crm_cld(icol,ix,iy,iz)
                           if (cld(icol,ilev) > 0) then
                              iclwp(icol,ilev) = crm_qc(icol,ix,iy,iz) * state%pdel(icol,ilev) / gravit / max(0.01_r8, cld(icol,ilev))
                              iciwp(icol,ilev) = crm_qi(icol,ix,iy,iz) * state%pdel(icol,ilev) / gravit / max(0.01_r8, cld(icol,ilev))
                           else
                              iclwp(icol,ilev) = 0
                              iciwp(icol,ilev) = 0
                           end if
                        end do  ! icol = 1,ncol
                     end do  ! iz = 1,crm_nz
                     ! Setup state arrays, which contain an extra level above model top to handle heating above the model
                     call set_rad_state(                                            &
                        state                      , cam_in                       , &
                        tmid_col(1:ncol,1:nlev_rad), tint_col(1:ncol,1:nlev_rad+1), &
                        pmid_col(1:ncol,1:nlev_rad), pint_col(1:ncol,1:nlev_rad+1)  &
                     )
                     call t_stopf('rad_set_state')

                     ! Set gas concentrations
                     call t_startf('rad_gas_concentrations')
                     call get_gas_vmr(icall, state, pbuf, active_gases, vmr_col(:,:,ktop:kbot))
                     vmr_col(:,:,1) = vmr_col(:,:,2)  ! Extra layer above model top
                     call t_stopf('rad_gas_concentrations')

                     ! Pack data
                     call t_startf('rad_pack_columns')
                     do icol = 1,ncol
                        j = _IDX321(icol, ix, iy, ncol, crm_nx_rad, crm_ny_rad)
                        cld_p(j,:) = cld(icol,:)
                        iclwp_p(j,:) = iclwp(icol,:)
                        iciwp_p(j,:) = iciwp(icol,:)
                        rel_p(j,:) = rel(icol,:)
                        rei_p(j,:) = rei(icol,:)
                        coszrs_packed(j) = coszrs(icol)
                        albedo_dir_packed(:,j) = albedo_dir_col(:,icol)
                        albedo_dif_packed(:,j) = albedo_dif_col(:,icol)
                        pmid_packed(j,:) = pmid_col(icol,:)
                        tmid_packed(j,:) = tmid_col(icol,:)
                        pint_packed(j,:) = pint_col(icol,:)
                        tint_packed(j,:) = tint_col(icol,:)
                        ! Note: leave top level empty for optics
                        aer_tau_bnd_lw_packed(j,ktop:kbot,:) = aer_tau_bnd_lw(icol,:,:)
                        aer_tau_bnd_sw_packed(j,ktop:kbot,:) = aer_tau_bnd_sw(icol,:,:)
                        aer_ssa_bnd_sw_packed(j,ktop:kbot,:) = aer_ssa_bnd_sw(icol,:,:)
                        aer_asm_bnd_sw_packed(j,ktop:kbot,:) = aer_asm_bnd_sw(icol,:,:)
                        vmr_packed(:,j,:) = vmr_col(:,icol,:)
                        day_indices_packed(j) = day_indices(icol)
                     end do
                     call t_stopf('rad_pack_columns')

                  end do  ! ix = 1,crm_nx_rad
               end do  ! iy = 1,crm_ny_rad

               ! Check (and possibly clip) values before passing to RRTMGP driver
               call handle_error(clip_values(tint_packed(1:ncol_tot,1:nlev_rad+1), get_min_temperature(), get_max_temperature(), trim(subname) // ' tint'), &
                  fatal=.false., warn=rrtmgp_enable_temperature_warnings)
               call handle_error(clip_values(aer_tau_bnd_sw_packed,  0._r8, huge(aer_tau_bnd_sw_packed), trim(subname) // ' aer_tau_bnd_sw', tolerance=1e-10_r8))
               call handle_error(clip_values(aer_ssa_bnd_sw_packed,  0._r8,                       1._r8, trim(subname) // ' aer_ssa_bnd_sw', tolerance=1e-10_r8))
               call handle_error(clip_values(aer_asm_bnd_sw_packed, -1._r8,                       1._r8, trim(subname) // ' aer_asm_bnd_sw', tolerance=1e-10_r8))
               call handle_error(clip_values(aer_tau_bnd_lw_packed,  0._r8, huge(aer_tau_bnd_lw_packed), trim(subname) // ' aer_tau_bnd_lw', tolerance=1e-10_r8))

               ! Do shortwave stuff...
               if (radiation_do('sw')) then

                  if (fixed_total_solar_irradiance<0) then
                     ! Get orbital eccentricity factor to scale total sky irradiance
                     tsi_scaling = get_eccentricity_factor()
                  else
                     ! For fixed TSI we divide by the default solar constant of 1360.9
                     ! At some point we will want to replace this with a method that 
                     ! retrieves the solar constant
                     tsi_scaling = fixed_total_solar_irradiance / 1360.9_r8
                  end if

                  ! Allocate shortwave fluxes (allsky and clearsky)
                  ! NOTE: fluxes defined at interfaces, so initialize to have vertical
                  ! dimension nlev_rad+1, while we initialized the RRTMGP input variables to
                  ! have vertical dimension nlev_rad (defined at midpoints).
                  call t_startf('rad_initialize_fluxes_sw')
                  call initialize_fluxes(ncol    , nlev_rad+1, nswbands, fluxes_allsky    , do_direct=.true.)
                  call initialize_fluxes(ncol    , nlev_rad+1, nswbands, fluxes_clrsky    , do_direct=.true.)
                  call initialize_fluxes(ncol_tot, nlev_rad+1, nswbands, fluxes_allsky_packed, do_direct=.true.)
                  call initialize_fluxes(ncol_tot, nlev_rad+1, nswbands, fluxes_clrsky_packed, do_direct=.true.)
                  call t_stopf('rad_initialize_fluxes_sw')

                  ! Calculate shortwave fluxes
                  ! If no daytime columns in this chunk, then we return zeros
                  call set_daynight_indices(coszrs_packed(1:ncol_tot), day_indices_packed(1:ncol_tot), night_indices_packed(1:ncol_tot))
                  nday = count(day_indices_packed(1:ncol_tot) > 0)
                  if (nday > 0) then
                     ! Compress to daytime-only arrays
                     call t_startf('rad_compress_day_columns')
                     do iday = 1,nday
                        icol = day_indices_packed(iday)
                        tmid_day(iday,:) = tmid_packed(icol,:)
                        pmid_day(iday,:) = pmid_packed(icol,:)
                        pint_day(iday,:) = pint_packed(icol,:)
                        albedo_dir_day(:,iday) = albedo_dir_packed(:,icol)
                        albedo_dif_day(:,iday) = albedo_dif_packed(:,icol)
                        coszrs_day(iday) = coszrs_packed(icol)
                        vmr_day(:,iday,:) = vmr_packed(:,icol,:)
                        cld_d(iday,:) = cld_p(icol,:)
                        iclwp_d(iday,:) = iclwp_p(icol,:)
                        iciwp_d(iday,:) = iciwp_p(icol,:)
                        rel_d(iday,:) = rel_p(icol,:)
                        rei_d(iday,:) = rei_p(icol,:)
                        aer_tau_bnd_day(iday,:,:) = aer_tau_bnd_sw_packed(icol,:,:)
                        aer_ssa_bnd_day(iday,:,:) = aer_ssa_bnd_sw_packed(icol,:,:)
                        aer_asm_bnd_day(iday,:,:) = aer_asm_bnd_sw_packed(icol,:,:)
                     end do
                     call t_stopf('rad_compress_day_columns')

                     ! Do cloud optics
                     call t_startf('rad_cloud_optics_sw')
                     cld_tau_bnd_day = 0._r8
                     cld_ssa_bnd_day = 0._r8
                     cld_asm_bnd_day = 0._r8
                     cld_tau_gpt_day = 0._r8
                     cld_ssa_gpt_day = 0._r8
                     cld_asm_gpt_day = 0._r8
                     call get_cloud_optics_sw( &
                        nday, pver, nswbands, &
                        cld_d, iclwp_d, iciwp_d, &
                        rel_d, rei_d, &
                        cld_tau_bnd_day(:,ktop:kbot,:), cld_ssa_bnd_day(:,ktop:kbot,:), cld_asm_bnd_day(:,ktop:kbot,:) &
                     )
                     ! Now reorder bands to be consistent with RRTMGP
                     ! We need to fix band ordering because the old input files assume RRTMG
                     ! band ordering, but this has changed in RRTMGP.
                     ! TODO: fix the input files themselves!
                     do icol = 1,nday
                        do ilay = 1,nlev_rad
                           cld_tau_bnd_day(icol,ilay,:) = reordered(cld_tau_bnd_day(icol,ilay,:), rrtmg_to_rrtmgp_swbands)
                           cld_ssa_bnd_day(icol,ilay,:) = reordered(cld_ssa_bnd_day(icol,ilay,:), rrtmg_to_rrtmgp_swbands)
                           cld_asm_bnd_day(icol,ilay,:) = reordered(cld_asm_bnd_day(icol,ilay,:), rrtmg_to_rrtmgp_swbands)
                        end do
                     end do
                     ! MCICA sampling to get cloud optical properties by gpoint/cloud state
                     call get_gpoint_bands_sw(gpoint_bands_sw)
                     call sample_cloud_optics_sw( &
                        nday, pver, nswgpts, gpoint_bands_sw, &
                        pmid_day(:,ktop:kbot), cld_d, &
                        cld_tau_bnd_day(:,ktop:kbot,:), cld_ssa_bnd_day(:,ktop:kbot,:), cld_asm_bnd_day(:,ktop:kbot,:), &
                        cld_tau_gpt_day(:,ktop:kbot,:), cld_ssa_gpt_day(:,ktop:kbot,:), cld_asm_gpt_day(:,ktop:kbot,:) &
                     )
                     call handle_error(clip_values(cld_tau_gpt_day,  0._r8, huge(cld_tau_gpt_day), trim(subname) // ' cld_tau_gpt_day', tolerance=1e-10_r8))
                     call handle_error(clip_values(cld_ssa_gpt_day,  0._r8,                1._r8, trim(subname) // ' cld_ssa_gpt_day', tolerance=1e-10_r8))
                     call handle_error(clip_values(cld_asm_gpt_day, -1._r8,                1._r8, trim(subname) // ' cld_asm_gpt_day', tolerance=1e-10_r8))
                     call t_stopf('rad_cloud_optics_sw')

                     ! Allocate shortwave fluxes (allsky and clearsky)
                     ! NOTE: fluxes defined at interfaces, so initialize to have vertical
                     ! dimension nlev_rad+1, while we initialized the RRTMGP input variables to
                     ! have vertical dimension nlev_rad (defined at midpoints).
                     call initialize_fluxes(nday, nlev_rad+1, nswbands, fluxes_allsky_day, do_direct=.true.)
                     call initialize_fluxes(nday, nlev_rad+1, nswbands, fluxes_clrsky_day, do_direct=.true.)

                     ! Do shortwave radiative transfer calculations
                     call t_startf('rad_rrtmgp_run_sw')
                     call rrtmgp_run_sw( &
                        size(active_gases), nday, nlev_rad, &
                        vmr_day(:,1:nday,1:nlev_rad), &
                        pmid_day(1:nday,1:nlev_rad), &
                        tmid_day(1:nday,1:nlev_rad), &
                        pint_day(1:nday,1:nlev_rad+1), &
                        coszrs_day(1:nday), &
                        albedo_dir_day(1:nswbands,1:nday), &
                        albedo_dif_day(1:nswbands,1:nday), &
                        cld_tau_gpt_day(1:nday,1:nlev_rad,1:nswgpts), cld_ssa_gpt_day(1:nday,1:nlev_rad,1:nswgpts), cld_asm_gpt_day(1:nday,1:nlev_rad,1:nswgpts), &
                        aer_tau_bnd_day(1:nday,1:nlev_rad,1:nswbands), aer_ssa_bnd_day(1:nday,1:nlev_rad,1:nswbands), aer_asm_bnd_day(1:nday,1:nlev_rad,1:nswbands), &
                        fluxes_allsky_day%flux_up    , fluxes_allsky_day%flux_dn    , fluxes_allsky_day%flux_net    , fluxes_allsky_day%flux_dn_dir    , &
                        fluxes_allsky_day%bnd_flux_up, fluxes_allsky_day%bnd_flux_dn, fluxes_allsky_day%bnd_flux_net, fluxes_allsky_day%bnd_flux_dn_dir, &
                        fluxes_clrsky_day%flux_up    , fluxes_clrsky_day%flux_dn    , fluxes_clrsky_day%flux_net    , fluxes_clrsky_day%flux_dn_dir    , &
                        fluxes_clrsky_day%bnd_flux_up, fluxes_clrsky_day%bnd_flux_dn, fluxes_clrsky_day%bnd_flux_net, fluxes_clrsky_day%bnd_flux_dn_dir, &
                        tsi_scaling &
                     )
                     call t_stopf('rad_rrtmgp_run_sw')

                     ! Expand fluxes from daytime-only arrays to full chunk arrays
                     call t_startf('rad_expand_fluxes_sw')
                     call expand_day_fluxes(fluxes_allsky_day, fluxes_allsky_packed, day_indices_packed(1:nday))
                     call expand_day_fluxes(fluxes_clrsky_day, fluxes_clrsky_packed, day_indices_packed(1:nday))
                     call t_stopf('rad_expand_fluxes_sw')

                     ! Clean up after ourselves
                     call free_fluxes(fluxes_allsky_day)
                     call free_fluxes(fluxes_clrsky_day)

                  else

                     call reset_fluxes(fluxes_allsky_packed)
                     call reset_fluxes(fluxes_clrsky_packed)

                  end if  ! nday > 0

                  ! Calculate heating rates
                  call t_startf('rad_heating_rate_sw')
                  call calculate_heating_rate(fluxes_allsky_packed%flux_up(1:ncol_tot,ktop:kbot+1), &
                     fluxes_allsky_packed%flux_dn(1:ncol_tot,ktop:kbot+1), &
                     pint_packed(1:ncol_tot,ktop:kbot+1), &
                     qrs_packed(1:ncol_tot,1:pver) &
                  )
                  call calculate_heating_rate(fluxes_clrsky_packed%flux_up(1:ncol_tot,ktop:kbot+1), &
                     fluxes_clrsky_packed%flux_dn(1:ncol_tot,ktop:kbot+1), &
                     pint_packed(1:ncol_tot,ktop:kbot+1), &
                     qrsc_packed(1:ncol_tot,1:pver) &
                  )
                  call t_stopf('rad_heating_rate_sw')

                  ! Calculate CRM domain averages
                  call t_startf('rad_average_fluxes_sw')
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, qrs_packed (1:ncol_tot,:), qrs (1:ncol,:))
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, qrsc_packed(1:ncol_tot,:), qrsc(1:ncol,:))
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_allsky_packed%flux_up, fluxes_allsky%flux_up)
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_allsky_packed%flux_dn, fluxes_allsky%flux_dn)
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_allsky_packed%flux_net   , fluxes_allsky%flux_net   )
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_allsky_packed%flux_dn_dir, fluxes_allsky%flux_dn_dir)
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_clrsky_packed%flux_up    , fluxes_clrsky%flux_up    )
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_clrsky_packed%flux_dn    , fluxes_clrsky%flux_dn    )
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_clrsky_packed%flux_net   , fluxes_clrsky%flux_net   )
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_clrsky_packed%flux_dn_dir, fluxes_clrsky%flux_dn_dir)
                  do iband = 1,nswbands
                     call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_allsky_packed%bnd_flux_up    (:,:,iband), fluxes_allsky%bnd_flux_up    (:,:,iband))
                     call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_allsky_packed%bnd_flux_dn    (:,:,iband), fluxes_allsky%bnd_flux_dn    (:,:,iband))
                     call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_allsky_packed%bnd_flux_net   (:,:,iband), fluxes_allsky%bnd_flux_net   (:,:,iband))
                     call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_allsky_packed%bnd_flux_dn_dir(:,:,iband), fluxes_allsky%bnd_flux_dn_dir(:,:,iband))
                     call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_clrsky_packed%bnd_flux_up    (:,:,iband), fluxes_clrsky%bnd_flux_up    (:,:,iband))
                     call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_clrsky_packed%bnd_flux_dn    (:,:,iband), fluxes_clrsky%bnd_flux_dn    (:,:,iband))
                     call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_clrsky_packed%bnd_flux_net   (:,:,iband), fluxes_clrsky%bnd_flux_net   (:,:,iband))
                     call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_clrsky_packed%bnd_flux_dn_dir(:,:,iband), fluxes_clrsky%bnd_flux_dn_dir(:,:,iband))
                  end do
                  call t_stopf('rad_average_fluxes_sw')

                  ! Map to CRM columns
                  do iy = 1,crm_ny_rad
                     do ix = 1,crm_nx_rad
                        do icol = 1,ncol
                           j = _IDX321(icol, ix, iy, ncol, crm_nx_rad, crm_ny_rad)
                           do iz = 1,crm_nz
                              ilev = pver - iz + 1
                              crm_qrs (icol,ix,iy,iz) = qrs_packed(j,ilev)
                              crm_qrsc(icol,ix,iy,iz) = qrsc_packed(j,ilev)
                           end do
                        end do
                     end do
                  end do

                  ! Send fluxes to history buffer
                  call output_fluxes_sw(icall, state, fluxes_allsky, fluxes_clrsky, qrs,  qrsc)
                  call outfld('CRM_QRS' , crm_qrs (1:ncol,:,:,:)/cpair, ncol, state%lchnk)
                  call outfld('CRM_QRSC', crm_qrsc(1:ncol,:,:,:)/cpair, ncol, state%lchnk)

                  ! Set net fluxes used by other components (land?) 
                  call set_net_fluxes_sw(fluxes_allsky, fsds, fsns, fsnt)

                  ! Set surface fluxes that are used by the land model
                  call export_surface_fluxes(fluxes_allsky, cam_out, 'shortwave')
                  
                  ! Free memory allocated for shortwave fluxes
                  call free_fluxes(fluxes_allsky)
                  call free_fluxes(fluxes_clrsky)
                  call free_fluxes(fluxes_allsky_packed)
                  call free_fluxes(fluxes_clrsky_packed)

               else
                  ! Back heating out of pdel scaled heating from previous step
                  qrs(1:ncol,1:pver) = qrs(1:ncol,1:pver) / state%pdel(1:ncol,1:pver)
               end if  ! dosw

               ! Do longwave stuff...
               if (radiation_do('lw')) then
                  ! Compute cloud optics
                  call t_startf('rad_cloud_optics_lw')
                  cld_tau_gpt_lw = 0._r8
                  cld_tau_bnd_lw = 0._r8
                  call get_cloud_optics_lw(ncol_tot, pver, nlwbands, cld_p, iclwp_p, iciwp_p, rei_p, cld_tau_bnd_lw(:,ktop:kbot,:))
                  ! Do mcica sampling of cloud optics
                  call get_gpoint_bands_lw(gpoint_bands_lw)
                  call sample_cloud_optics_lw( &
                     ncol_tot, pver, nlwgpts, gpoint_bands_lw, &
                     pmid_packed(:,ktop:kbot), cld_p, &
                     cld_tau_bnd_lw(:,ktop:kbot,:), cld_tau_gpt_lw(:,ktop:kbot,:) &
                  )
                  ! Save CRM cloud optics for cosp
                  call handle_error(clip_values(cld_tau_gpt_lw,  0._r8, huge(cld_tau_gpt_lw), trim(subname) // ' cld_tau_gpt_lw', tolerance=1e-10_r8))
                  call t_stopf('rad_cloud_optics_lw')

                  ! Allocate longwave outputs
                  ! NOTE: fluxes defined at interfaces, so initialize to have vertical dimension nlev_rad+1
                  call t_startf('rad_initialize_fluxes_lw')
                  call initialize_fluxes(ncol, nlev_rad+1, nlwbands, fluxes_allsky)
                  call initialize_fluxes(ncol, nlev_rad+1, nlwbands, fluxes_clrsky)
                  call initialize_fluxes(ncol_tot, nlev_rad+1, nlwbands, fluxes_allsky_packed)
                  call initialize_fluxes(ncol_tot, nlev_rad+1, nlwbands, fluxes_clrsky_packed)
                  call t_stopf('rad_initialize_fluxes_lw')

                  ! Calculate longwave fluxes
                  call t_startf('rad_fluxes_lw')
                  call radiation_driver_lw(                                                &
                     size(active_gases), ncol_tot, nlev_rad, vmr_packed(:,1:ncol_tot,1:nlev_rad),                                     &
                     surface_emissivity(1:nlwbands,1:ncol_tot),                            &
                     pmid_packed(1:ncol_tot,1:nlev_rad  ), tmid_packed(1:ncol_tot,1:nlev_rad  ), &
                     pint_packed(1:ncol_tot,1:nlev_rad+1), tint_packed(1:ncol_tot,1:nlev_rad+1), &
                     cld_tau_gpt_lw                      , aer_tau_bnd_lw_packed,                &
                     fluxes_allsky_packed                , fluxes_clrsky_packed                  &
                  )
                  call t_stopf('rad_fluxes_lw')

                  ! Calculate heating rates
                  call t_startf('rad_heating_lw')
                  call calculate_heating_rate(fluxes_allsky_packed%flux_up(:,ktop:kbot+1), &
                                              fluxes_allsky_packed%flux_dn(:,ktop:kbot+1), &
                                              pint_packed(1:ncol_tot,ktop:kbot+1)        , &
                                              qrl_packed(1:ncol_tot,1:pver)                )
                  call calculate_heating_rate(fluxes_clrsky_packed%flux_up(:,ktop:kbot+1), &
                                              fluxes_clrsky_packed%flux_dn(:,ktop:kbot+1), &
                                              pint_packed(1:ncol_tot,ktop:kbot+1)        , &
                                              qrlc_packed(1:ncol_tot,1:pver)               )
                  call t_stopf('rad_heating_lw')

                  ! Calculate domain averages
                  call t_startf('rad_average_fluxes_lw')
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, qrl_packed (1:ncol_tot,:)    , qrl (1:ncol,:))
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, qrlc_packed(1:ncol_tot,:)    , qrlc(1:ncol,:))
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_allsky_packed%flux_up , fluxes_allsky%flux_up )
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_allsky_packed%flux_dn , fluxes_allsky%flux_dn )
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_allsky_packed%flux_net, fluxes_allsky%flux_net)
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_clrsky_packed%flux_up , fluxes_clrsky%flux_up )
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_clrsky_packed%flux_dn , fluxes_clrsky%flux_dn )
                  call average_packed_array(ncol, crm_nx_rad, crm_ny_rad, fluxes_clrsky_packed%flux_net, fluxes_clrsky%flux_net)
                  call t_stopf('rad_average_fluxes_lw')
                            
                  ! Map to CRM columns
                  do iy = 1,crm_ny_rad
                     do ix = 1,crm_nx_rad
                        do icol = 1,ncol
                           j = _IDX321(icol, ix, iy, ncol, crm_nx_rad, crm_ny_rad)
                           do iz = 1,crm_nz
                              ilev = pver - iz + 1
                              crm_qrl(icol,ix,iy,iz) = qrl_packed(j,ilev)
                              crm_qrlc(icol,ix,iy,iz) = qrlc_packed(j,ilev)
                           end do
                        end do
                     end do
                  end do

                  ! Send fluxes and heating rates to history buffer
                  call output_fluxes_lw(icall, state, fluxes_allsky, fluxes_clrsky, qrl, qrlc)
                  call outfld('CRM_QRL', crm_qrl(1:ncol,:,:,:)/cpair, ncol, state%lchnk)
                  call outfld('CRM_QRLC', crm_qrlc(1:ncol,:,:,:)/cpair, ncol, state%lchnk)

                  ! Set net fluxes used in other components
                  call set_net_fluxes_lw(fluxes_allsky, flns, flnt)

                  ! Export surface fluxes that are used by the land model
                  call export_surface_fluxes(fluxes_allsky, cam_out, 'longwave')

                  ! Free memory allocated for fluxes
                  call free_fluxes(fluxes_allsky)
                  call free_fluxes(fluxes_clrsky)
                  call free_fluxes(fluxes_allsky_packed)
                  call free_fluxes(fluxes_clrsky_packed)

               else
                  ! Back heating out of pdel scaled heating rate from previous step
                  qrl(1:ncol,1:pver) = qrl(1:ncol,1:pver) / state%pdel(1:ncol,1:pver)
               end if
            end if  ! active calls
         end do  ! loop over diagnostic calls

         ! Update net CRM heating tendency
         ! TODO: should we be updating this every timestep, even if using
         ! previous heating tendency since pdel would have changed?
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_QRAD'), crm_qrad)
         crm_qrad = 0
         do iz = 1,crm_nz
            do iy = 1,crm_ny_rad
               do ix = 1,crm_nx_rad
                  do icol = 1,ncol
                     ilev = pver - iz + 1
                     crm_qrad(icol,ix,iy,iz) = (crm_qrs(icol,ix,iy,iz) + crm_qrl(icol,ix,iy,iz)) / cpair
                     crm_qrad(icol,ix,iy,iz) = crm_qrad(icol,ix,iy,iz) * state%pdel(icol,ilev)
                  end do
               end do
            end do
         end do
         call outfld('CRM_QRAD', crm_qrad(1:ncol,:,:,:), ncol, state%lchnk)

      end if  ! sw or lw

      ! If we ran radiation this timestep, check if we should run COSP
      if (radiation_do('sw') .and. radiation_do('lw')) then
         if (docosp) then
            ! Advance counter and run COSP if new count value is equal to cosp_nradsteps
            cosp_cnt(state%lchnk) = cosp_cnt(state%lchnk) + 1
            if (cosp_cnt(state%lchnk) == cosp_nradsteps) then
               ! Extract LW and SW bands we need
               cosp_lwband = get_band_index_lw(10.5_r8, 'micron')
               cosp_swband = get_band_index_sw(0.67_r8, 'micron')
               do ilay = 1,pver
                  do j = 1,ncol_tot
                     dems_packed(j,ilay) = 1._r8 - exp(-cld_tau_bnd_lw(j,ilay+1,cosp_lwband))
                  end do
               end do
               do ilay = 1,pver
                  do j = 1,nday
                     icol = day_indices_packed(j)
                     dtau_packed(j,ilay) = cld_tau_bnd_day(icol,ilay+1,cosp_swband)
                  end do
               end do
               ! Call cosp
               call t_startf('cospsimulator_intr_run')
               call cospsimulator_intr_run( &
                  ncol_tot, state, pbuf, cam_in, coszrs, &
                  dtau_packed, dems_packed &
               )
               call t_stopf('cospsimulator_intr_run')
               ! Reset counter
               cosp_cnt(state%lchnk) = 0
            end if
         end if
      end if

      ! Compute net radiative heating tendency
      call t_startf('radheat_tend')
      call radheat_tend(state, pbuf, ptend, &
                        qrl, qrs, &
                        fsns, fsnt, flns, flnt, &
                        cam_in%asdir, net_flux)
      call t_stopf('radheat_tend')

      ! Compute heating rate for dtheta/dt
      call t_startf ('rad_heating_rate')
      do ilay = 1, pver
         do icol = 1, ncol
            hr(icol,ilay) = (qrs(icol,ilay) + qrl(icol,ilay)) / cpair * (1.e5_r8 / state%pmid(icol,ilay))**cappa
         end do
      end do
      call outfld('HR', hr(1:ncol,:), ncol, state%lchnk)
      call t_stopf ('rad_heating_rate')

      ! convert radiative heating rates to Q*dp to carry across timesteps
      ! for energy conservation
      qrs(1:ncol,1:pver) = qrs(1:ncol,1:pver) * state%pdel(1:ncol,1:pver)
      qrl(1:ncol,1:pver) = qrl(1:ncol,1:pver) * state%pdel(1:ncol,1:pver)

   end subroutine radiation_tend

   subroutine radiation_driver_sw(ncol, &
                                  gas_vmr, &
                                  pmid, pint, tmid, albedo_dir, albedo_dif, coszrs, &
                                  cld_tau_gpt, cld_ssa_gpt, cld_asm_gpt, &
                                  aer_tau_bnd, aer_ssa_bnd, aer_asm_bnd, &
                                  fluxes_allsky, fluxes_clrsky, tsi_scaling)
     
      use perf_mod, only: t_startf, t_stopf

      ! Inputs
      integer, intent(in) :: ncol
      type(fluxes_t), intent(inout) :: fluxes_allsky, fluxes_clrsky
      real(r8), intent(in), dimension(:,:,:) :: gas_vmr
      real(r8), intent(in), dimension(:,:) :: pmid, pint, tmid
      real(r8), intent(in), dimension(:,:) :: albedo_dir, albedo_dif
      real(r8), intent(in), dimension(:) :: coszrs
      real(r8), intent(in), dimension(:,:,:) :: cld_tau_gpt, cld_ssa_gpt, cld_asm_gpt
      real(r8), intent(in), dimension(:,:,:) :: aer_tau_bnd, aer_ssa_bnd, aer_asm_bnd
      ! Scaling factor for total sky irradiance; used to account for orbital
      ! eccentricity, and could be used to scale total sky irradiance for different
      ! climates as well (i.e., paleoclimate simulations)
      real(r8), intent(in) :: tsi_scaling

      ! Compressed daytime-only arrays
      real(r8), dimension(ncol) :: coszrs_day
      real(r8), dimension(nswbands,ncol) :: albedo_dir_day, albedo_dif_day
      real(r8), dimension(ncol,nlev_rad) :: pmid_day, tmid_day
      real(r8), dimension(ncol,nlev_rad+1) :: pint_day
      real(r8), dimension(size(gas_vmr, 1),ncol,nlev_rad) :: gas_vmr_day
      real(r8), dimension(ncol,pver,nswgpts) :: cld_tau_gpt_day, cld_ssa_gpt_day, cld_asm_gpt_day
      real(r8), dimension(ncol,pver,nswbands) :: aer_tau_bnd_day, aer_ssa_bnd_day, aer_asm_bnd_day
      type(fluxes_t) :: fluxes_allsky_day, fluxes_clrsky_day

      real(r8), dimension(size(cld_tau_gpt,1),nlev_rad,size(cld_tau_gpt,3)) :: cld_tau_gpt_rad
      real(r8), dimension(size(cld_ssa_gpt,1),nlev_rad,size(cld_ssa_gpt,3)) :: cld_ssa_gpt_rad
      real(r8), dimension(size(cld_asm_gpt,1),nlev_rad,size(cld_asm_gpt,3)) :: cld_asm_gpt_rad
      real(r8), dimension(size(aer_tau_bnd,1),nlev_rad,size(aer_tau_bnd,3)) :: aer_tau_bnd_rad
      real(r8), dimension(size(aer_ssa_bnd,1),nlev_rad,size(aer_ssa_bnd,3)) :: aer_ssa_bnd_rad
      real(r8), dimension(size(aer_asm_bnd,1),nlev_rad,size(aer_asm_bnd,3)) :: aer_asm_bnd_rad

      ! Gathered indicies of day and night columns 
      ! chunk_column_index = day_indices(daylight_column_index)
      integer :: iday, icol
      integer :: nday  ! Number of daylight columns
      integer :: day_indices(ncol), night_indices(ncol)   ! Indicies of daylight coumns


      ! Everybody needs a name
      character(*), parameter :: subroutine_name = 'radiation_driver_sw'


      ! Gather night/day column indices for subsetting SW inputs; we only want to
      ! do the shortwave radiative transfer during the daytime to save
      ! computational cost (and because RRTMGP will fail for cosine solar zenith
      ! angles less than or equal to zero)
      call set_daynight_indices(coszrs(1:ncol), day_indices(1:ncol), night_indices(1:ncol))
      nday = count(day_indices(1:ncol) > 0)

      ! If no daytime columns in this chunk, then we return zeros
      if (nday == 0) then
         call reset_fluxes(fluxes_allsky)
         call reset_fluxes(fluxes_clrsky)
         return
      end if

      ! Compress to daytime-only arrays
      ! TODO: do this BEFORE computing optics, because packing cloud optics is expensive
      call t_startf('rad_compress_day_columns')
      do iday = 1,nday
         icol = day_indices(iday)
         tmid_day(iday,:) = tmid(icol,:)
         pmid_day(iday,:) = pmid(icol,:)
         pint_day(iday,:) = pint(icol,:)
         albedo_dir_day(:,iday) = albedo_dir(:,icol)
         albedo_dif_day(:,iday) = albedo_dif(:,icol)
         coszrs_day(iday) = coszrs(icol)
         gas_vmr_day(:,iday,:) = gas_vmr(:,icol,:)
         cld_tau_gpt_rad(iday,:,:) = cld_tau_gpt(icol,:,:)
         cld_ssa_gpt_rad(iday,:,:) = cld_ssa_gpt(icol,:,:)
         cld_asm_gpt_rad(iday,:,:) = cld_asm_gpt(icol,:,:)
         aer_tau_bnd_rad(iday,:,:) = aer_tau_bnd(icol,:,:)
         aer_ssa_bnd_rad(iday,:,:) = aer_ssa_bnd(icol,:,:)
         aer_asm_bnd_rad(iday,:,:) = aer_asm_bnd(icol,:,:)
      end do
      call t_stopf('rad_compress_day_columns')

      ! Allocate shortwave fluxes (allsky and clearsky)
      ! NOTE: fluxes defined at interfaces, so initialize to have vertical
      ! dimension nlev_rad+1, while we initialized the RRTMGP input variables to
      ! have vertical dimension nlev_rad (defined at midpoints).
      call initialize_fluxes(nday, nlev_rad+1, nswbands, fluxes_allsky_day, do_direct=.true.)
      call initialize_fluxes(nday, nlev_rad+1, nswbands, fluxes_clrsky_day, do_direct=.true.)

      ! Do shortwave radiative transfer calculations
      call t_startf('rad_rrtmgp_run_sw')
      call rrtmgp_run_sw( &
         size(active_gases), nday, nlev_rad, &
         gas_vmr_day(:,1:nday,1:nlev_rad), &
         pmid_day(1:nday,1:nlev_rad), &
         tmid_day(1:nday,1:nlev_rad), &
         pint_day(1:nday,1:nlev_rad+1), &
         coszrs_day(1:nday), &
         albedo_dir_day(1:nswbands,1:nday), &
         albedo_dif_day(1:nswbands,1:nday), &
         cld_tau_gpt_rad(1:nday,1:nlev_rad,1:nswgpts), cld_ssa_gpt_rad(1:nday,1:nlev_rad,1:nswgpts), cld_asm_gpt_rad(1:nday,1:nlev_rad,1:nswgpts), &
         aer_tau_bnd_rad(1:nday,1:nlev_rad,1:nswbands), aer_ssa_bnd_rad(1:nday,1:nlev_rad,1:nswbands), aer_asm_bnd_rad(1:nday,1:nlev_rad,1:nswbands), &
         fluxes_allsky_day%flux_up    , fluxes_allsky_day%flux_dn    , fluxes_allsky_day%flux_net    , fluxes_allsky_day%flux_dn_dir    , &
         fluxes_allsky_day%bnd_flux_up, fluxes_allsky_day%bnd_flux_dn, fluxes_allsky_day%bnd_flux_net, fluxes_allsky_day%bnd_flux_dn_dir, &
         fluxes_clrsky_day%flux_up    , fluxes_clrsky_day%flux_dn    , fluxes_clrsky_day%flux_net    , fluxes_clrsky_day%flux_dn_dir    , &
         fluxes_clrsky_day%bnd_flux_up, fluxes_clrsky_day%bnd_flux_dn, fluxes_clrsky_day%bnd_flux_net, fluxes_clrsky_day%bnd_flux_dn_dir, &
         tsi_scaling &
      )
      call t_stopf('rad_rrtmgp_run_sw')

      ! Expand fluxes from daytime-only arrays to full chunk arrays
      call t_startf('rad_expand_fluxes_sw')
      call expand_day_fluxes(fluxes_allsky_day, fluxes_allsky, day_indices(1:nday))
      call expand_day_fluxes(fluxes_clrsky_day, fluxes_clrsky, day_indices(1:nday))
      call t_stopf('rad_expand_fluxes_sw')

      ! Clean up after ourselves
      call free_fluxes(fluxes_allsky_day)
      call free_fluxes(fluxes_clrsky_day)

   end subroutine radiation_driver_sw

   !----------------------------------------------------------------------------

   subroutine radiation_driver_lw(ngas, ncol, nlev, gas_vmr, surface_emissivity, &
                                  pmid, tmid, pint, tint, &
                                  cld_tau_gpt, aer_tau_bnd, &
                                  fluxes_allsky, fluxes_clrsky)

      use perf_mod, only: t_startf, t_stopf

      integer, intent(in) :: ngas, ncol, nlev
      real(r8), intent(in) :: gas_vmr(:,:,:)
      real(r8), intent(in) :: surface_emissivity(:,:)
      real(r8), intent(in) :: pmid(:,:), tmid(:,:), pint(:,:), tint(:,:)
      real(r8), intent(in) :: cld_tau_gpt(:,:,:), aer_tau_bnd(:,:,:)
      type(fluxes_t), intent(inout) :: fluxes_allsky, fluxes_clrsky

      ! Compute fluxes
      call t_startf('rrtmgp_run_lw')
      call rrtmgp_run_lw( &
         ngas, ncol, nlev, &
         gas_vmr(:,1:ncol,:), &
         pmid(1:ncol,1:nlev), tmid(1:ncol,1:nlev), pint(1:ncol,1:nlev+1), tint(1:ncol,1:nlev+1), &
         surface_emissivity(1:nlwbands,1:ncol), &
         cld_tau_gpt(1:ncol,:,:)  , aer_tau_bnd(1:ncol,:,:)  , &
         fluxes_allsky%flux_up    , fluxes_allsky%flux_dn    , fluxes_allsky%flux_net    , &
         fluxes_allsky%bnd_flux_up, fluxes_allsky%bnd_flux_dn, fluxes_allsky%bnd_flux_net, &
         fluxes_clrsky%flux_up    , fluxes_clrsky%flux_dn    , fluxes_clrsky%flux_net    , &
         fluxes_clrsky%bnd_flux_up, fluxes_clrsky%bnd_flux_dn, fluxes_clrsky%bnd_flux_net  &
         )
      call t_stopf('rrtmgp_run_lw')
   end subroutine radiation_driver_lw

   !----------------------------------------------------------------------------

   ! Utility routine to compute domain averages of CRM data that has been
   ! "packed" into a single dimension to hold both global GCM columns and CRM
   ! "subcolumns" within each GCM column. Input will be 2D arrays dimensioned
   ! (ncol_tot,nlev) where ncol_tot is the total number of CRM columns ncol *
   ! crm_nx_rad * crm_ny_rad and nlev is number of vertical levels. Output will
   ! be 2D arrays dimensioned (ncol,nlev), where the averaging has been done
   ! over the CRM columns.
   subroutine average_packed_array(ncol, nx, ny, array_packed, array_avg)
      integer, intent(in) :: ncol, nx, ny
      real(r8), intent(in) :: array_packed(:,:)
      real(r8), intent(out) :: array_avg(:,:)
      integer :: ncol_tot
      real(r8) :: area_factor
      integer :: icol, ix, iy, iz, j

      if (nx * ny > 1) then
         area_factor = 1._r8 / (nx * ny)
      else
         area_factor = 1
      end if
      array_avg = 0
      ncol_tot = ncol * nx * ny
      call assert(size(array_packed,1) == ncol_tot, 'size(array_packed,1) /= ncol_tot')
      call assert(size(array_packed,2) == size(array_avg,2), 'size(array_packed,2) /= size(array_avg,2)')
      do iz = 1,size(array_packed,2)
         do iy = 1,ny
            do ix = 1,nx
               do icol = 1,ncol
                  j = _IDX321(icol, ix, iy, ncol, nx, ny)
                  array_avg(icol,iz) = array_avg(icol,iz) + array_packed(j,iz) * area_factor
               end do
            end do
         end do 
      end do
   end subroutine average_packed_array

   !----------------------------------------------------------------------------

   subroutine export_surface_fluxes(fluxes, cam_out, band)
      use camsrfexch, only: cam_out_t

      type(fluxes_t), intent(in) :: fluxes
      type(cam_out_t), intent(inout) :: cam_out
      character(len=*), intent(in) :: band
      integer :: icol
      real(r8), dimension(size(fluxes%bnd_flux_dn,1), &
                          size(fluxes%bnd_flux_dn,2), &
                          size(fluxes%bnd_flux_dn,3)) :: flux_dn_diffuse

      ! TODO: check this code! This seems to differ from what is in RRTMG. 
      !
      ! Export surface fluxes
      !
      ! To break the fluxes into the UV/vis and near-IR bands use the same scheme 
      ! as for the albedos which is hardcoded for 14 spectral bands.
      !
      ! sols(pcols)      Direct solar rad on surface (< 0.7)
      ! soll(pcols)      Direct solar rad on surface (>= 0.7)
      !
      ! Near-IR bands (1-10), 820-16000 cm-1, 0.625-12.195 microns
      !
      ! Put half of band 10 in each of the UV/visible and near-IR values,
      ! since this band straddles 0.7 microns:
      !
      ! UV/visible bands 11-14, 16000-50000 cm-1, 0.200-0.625 micron
      !
      ! NOTE: bands are shifted relative to RRTMG! Band 10 used to be band 9, band
      ! 1 used to be band 14.
      !
      ! TODO: this hard-coded band mapping is BAD! We need to do this more
      ! intelligently.
      if (trim(band) == 'shortwave') then

         ! Reset fluxes
         cam_out%soll = 0
         cam_out%sols = 0
         cam_out%solld = 0
         cam_out%solsd = 0

         ! Calculate diffuse flux from total and direct
         flux_dn_diffuse = fluxes%bnd_flux_dn - fluxes%bnd_flux_dn_dir

         ! Calculate broadband surface solar fluxes (UV/visible vs near IR) for
         ! each column.
         do icol = 1,size(fluxes%bnd_flux_dn, 1)

            ! Direct fluxes
            cam_out%soll(icol) &
               = sum(fluxes%bnd_flux_dn_dir(icol,kbot+1,1:9)) &
               + 0.5_r8 * fluxes%bnd_flux_dn_dir(icol,kbot+1,10)
            cam_out%sols(icol) &
               = 0.5_r8 * fluxes%bnd_flux_dn_dir(icol,kbot+1,10) &
               + sum(fluxes%bnd_flux_dn_dir(icol,kbot+1,11:14))

            ! Diffuse fluxes
            cam_out%solld(icol) &
               = sum(flux_dn_diffuse(icol,kbot+1,1:9)) &
               + 0.5_r8 * flux_dn_diffuse(icol,kbot+1,10)
            cam_out%solsd(icol) &
               = 0.5_r8 * flux_dn_diffuse(icol,kbot+1,10) &
               + sum(flux_dn_diffuse(icol,kbot+1,11:14))

            ! Net shortwave flux at surface
            cam_out%netsw(icol) = fluxes%flux_net(icol,kbot+1)
         end do
      else if (trim(band) == 'longwave') then
         do icol = 1,size(fluxes%flux_dn, 1)
            cam_out%flwds(icol) = fluxes%flux_dn(icol,kbot+1)
         end do
      else
         call endrun('flux_type ' // band // ' not known.')
      end if

   end subroutine export_surface_fluxes

   !----------------------------------------------------------------------------

   subroutine set_net_fluxes_sw(fluxes, fsds, fsns, fsnt)

      type(fluxes_t), intent(in) :: fluxes
      real(r8), intent(inout) :: fsds(:)
      real(r8), intent(inout) :: fsns(:)
      real(r8), intent(inout) :: fsnt(:)

      integer :: ncol

      ! Copy data
      ncol = size(fluxes%flux_up, 1)
      fsds(1:ncol) = fluxes%flux_dn(1:ncol,kbot+1)
      fsns(1:ncol) = fluxes%flux_dn(1:ncol,kbot+1) - fluxes%flux_up(1:ncol,kbot+1)
      fsnt(1:ncol) = fluxes%flux_dn(1:ncol,ktop) - fluxes%flux_up(1:ncol,ktop)

   end subroutine set_net_fluxes_sw

   !----------------------------------------------------------------------------

   subroutine set_net_fluxes_lw(fluxes, flns, flnt)

      type(fluxes_t), intent(in) :: fluxes
      real(r8), intent(inout) :: flns(:)
      real(r8), intent(inout) :: flnt(:)
      integer :: ncol

      ! Copy data
      ncol = size(fluxes%flux_up, 1)
      flns(1:ncol) = fluxes%flux_up(1:ncol,kbot+1) - fluxes%flux_dn(1:ncol,kbot+1)
      flnt(1:ncol) = fluxes%flux_up(1:ncol,ktop) - fluxes%flux_dn(1:ncol,ktop)

   end subroutine set_net_fluxes_lw

   !----------------------------------------------------------------------------

   subroutine set_daynight_indices(coszrs, day_indices, night_indices)
      ! Input: cosine of solar zenith angle
      real(r8), intent(in) :: coszrs(:)

      ! Output: array of indices to daytime columns
      integer, intent(inout) :: day_indices(:), night_indices(:)

      ! Loop indices; icol is index to physics columns in current chunk and iday is
      ! index to daytime indices
      integer :: icol, iday, inight

      ! Subroutine name for error messages
      character(len=128) :: subname = 'set_daynight_indices'

      ! Initialize array of daytime indices to be all zero. If any zeros exist when
      ! we are done, something went wrong.
      day_indices(:) = 0
      night_indices(:) = 0

      ! Loop over columns and identify daytime columns as those where the cosine
      ! solar zenith angle exceeds zero. Note that we wrap the setting of
      ! day_indices in an if-then to make sure we are not accesing day_indices out
      ! of bounds, and stopping with an informative error message if we do for some
      ! reason.
      iday = 0
      inight = 0
      do icol = 1,size(coszrs)
         if (coszrs(icol) > 0._r8) then
            iday = iday + 1
            day_indices(iday) = icol
         else
            inight = inight + 1
            night_indices(inight) = icol
         end if
      end do

      ! Check indices
      call assert(count(day_indices > 0) == count(coszrs > 0), &
                  trim(subname) // ': count(day_indices > 0) != count(coszrs > 0)')
      call assert(count(night_indices > 0) == count(coszrs <= 0), &
                  trim(subname) // ': count(night_indices > 0) != count(coszrs <= 0)')
      call assert(count(day_indices > 0) + count(night_indices > 0) == size(coszrs), &
                  trim(subname) // ': day_indices + night_indices != size(coszrs)')
      call assert_range(day_indices, 0, size(coszrs), trim(subname) // 'day_indices')

   end subroutine set_daynight_indices


   ! Subroutine to calculate the solar insolation, accounting for orbital
   ! eccentricity and solar variability
   function get_eccentricity_factor() result(eccentricity_factor)
      use cam_control_mod, only: eccen, mvelpp, lambm0, obliqr
      use shr_orb_mod, only: shr_orb_decl
      use time_manager, only: get_curr_calday

      ! Variables needed to use shr_orb_decl to get orbital eccentricity factor 
      ! (earth-sun distance factor)
      real(r8) :: calday ! Calendar day, including fraction
      real(r8) :: solar_declination    ! Solar declination angle in rad
      real(r8) :: eccentricity_factor  ! Earth-sun distance factor (ie. (1/r)**2)

      ! Get orbital eccentricity factor for the current calendar day
      calday = get_curr_calday()
      call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, &
                        solar_declination, eccentricity_factor)

   end function get_eccentricity_factor


   ! Get and set cosine of the solar zenith angle for all columns in a physics chuck
   ! based on input physics_state object and timestep. This routine serves mainly
   ! as a wrapper for the "zenith" subroutine that handles the task of grabbing the
   ! appropriate calendar day for this timestep and the latitude and longitude 
   ! values for the columns in the current chunk, and containing addition CAM
   ! module use to improve readability of the main radiation_tend routine.
   subroutine set_cosine_solar_zenith_angle(state, dt, coszrs)
      use physics_types,   only: physics_state
      use phys_grid,       only: get_rlat_all_p, get_rlon_all_p
      use time_manager,    only: get_curr_calday
      use orbit,           only: zenith

      ! Inputs
      type(physics_state), intent(in) :: state
      real(r8), intent(in) :: dt

      ! Output cosine solar zenith angle. Note that this is declared as an
      ! assumed-shape array, but below we will require coszrs to be at least as
      ! large as state%ncol (because latitudes and longitudes are only defined for
      ! the state%ncol columns that exist in the current physics "chunk")
      real(r8), intent(inout) :: coszrs(:)

      ! Local variables
      real(r8) :: calday  ! current calendar day
      real(r8) :: clat(size(coszrs))  ! current latitudes(radians)
      real(r8) :: clon(size(coszrs))  ! current longitudes(radians)

      ! Make sure desired coszrs has the correct shape. The "zenith" subroutine
      ! expects input arrays to have shape state%ncol, although this should
      ! probably be relaxed to take assumed-shape arrays.
      call assert(size(coszrs) >= state%ncol, 'size(coszrs) < ncol')

      ! Get solar zenith angle from CAM utility routine. The "zenith" routine needs
      ! the current calendar day, and the latitude and longitudes for all columns
      ! in the current physics "chunk", so we use CAM routines to retrieve those
      ! values here.
      calday = get_curr_calday()
      call get_rlat_all_p(state%lchnk, state%ncol, clat(1:state%ncol))
      call get_rlon_all_p(state%lchnk, state%ncol, clon(1:state%ncol))

      ! Call zenith to calculate cosine solar zenith angle. Note we only pass the
      ! first state%ncol columns in case coszrs was allocated to be larger than
      ! ncol in the calling routine (i.e, pcols).
      call zenith(calday, clat(1:state%ncol), clon(1:state%ncol), &
                  coszrs(1:state%ncol), state%ncol, dt)
   end subroutine set_cosine_solar_zenith_angle


   ! Set surface albedos from cam surface exchange object for direct and diffuse
   ! beam radiation. This code was copied from the RRTMG implementation, but moved
   ! to a subroutine with some better variable names.
   subroutine set_albedo(cam_in, albedo_dir, albedo_dif)
      use camsrfexch, only: cam_in_t
      use radiation_utils, only: clip_values

      type(cam_in_t), intent(in) :: cam_in
      real(r8), intent(inout) :: albedo_dir(:,:)   ! surface albedo, direct radiation
      real(r8), intent(inout) :: albedo_dif(:,:)  ! surface albedo, diffuse radiation

      ! Local namespace
      real(r8), dimension(nswbands) :: lower_bounds, upper_bounds
      integer :: ncol, iband
      character(len=10) :: subname = 'set_albedo'

      ! Check dimension sizes of output arrays.
      ! albedo_dir and albedo_dif should have sizes nswbands,ncol, but ncol
      ! can change so we just check that it is less than or equal to pcols (the
      ! maximum size ncol is ever allowed to be).
      call assert(size(albedo_dir, 1) == nswbands, &
                  'set_albedo: size(albedo_dir, 1) /= nswbands')
      call assert(size(albedo_dir, 2) <= pcols, &
                  'set_albedo: size(albedo_dir, 2) > pcols')
      call assert(all(shape(albedo_dir) == shape(albedo_dif)), &
                  'set_albedo: albedo_dir and albedo_dif have inconsistent shapes')
      
      ncol = size(albedo_dir, 2)

      ! Initialize albedo
      albedo_dir(:,:) = 0._r8
      albedo_dif(:,:) = 0._r8

      ! Albedos are input as broadband (visible, and near-IR), and we need to map
      ! these to appropriate bands. Bands are categorized broadly as "visible" or
      ! "infrared" based on wavenumber, so we get the wavenumber limits here
      call get_sw_spectral_boundaries(lower_bounds, upper_bounds, 'cm^-1')

      ! We need to reorder the spectral bounds since we store them in RRTMG
      ! order in radconstants!
      lower_bounds = reordered(lower_bounds, rrtmg_to_rrtmgp_swbands)
      upper_bounds = reordered(upper_bounds, rrtmg_to_rrtmgp_swbands)

      ! Loop over bands, and determine for each band whether it is broadly in the
      ! visible or infrared part of the spectrum (visible or "not visible")
      do iband = 1,nswbands
         if (is_visible(lower_bounds(iband)) .and. &
             is_visible(upper_bounds(iband))) then

            ! Entire band is in the visible
            albedo_dir(iband,1:ncol) = cam_in%asdir(1:ncol)
            albedo_dif(iband,1:ncol) = cam_in%asdif(1:ncol)

         else if (.not.is_visible(lower_bounds(iband)) .and. &
                  .not.is_visible(upper_bounds(iband))) then

            ! Entire band is in the longwave (near-infrared)
            albedo_dir(iband,1:ncol) = cam_in%aldir(1:ncol)
            albedo_dif(iband,1:ncol) = cam_in%aldif(1:ncol)

         else

            ! Band straddles the visible to near-infrared transition, so we take
            ! the albedo to be the average of the visible and near-infrared
            ! broadband albedos
            albedo_dir(iband,1:ncol) = 0.5 * (cam_in%aldir(1:ncol) + cam_in%asdir(1:ncol))
            albedo_dif(iband,1:ncol) = 0.5 * (cam_in%aldif(1:ncol) + cam_in%asdif(1:ncol))

         end if
      end do

      ! Check values and clip if necessary (albedos should not be larger than 1)
      ! NOTE: this does actually issue warnings for albedos larger than 1, but this
      ! was never checked for RRTMG, so albedos will probably be slight different
      ! than the implementation in RRTMG!
      call handle_error(clip_values( &
         albedo_dir, 0._r8, 1._r8, trim(subname) // ': albedo_dir', tolerance=0.01_r8) &
      )
      call handle_error(clip_values( &
         albedo_dif, 0._r8, 1._r8, trim(subname) // ': albedo_dif', tolerance=0.01_r8) &
      )


   end subroutine set_albedo

   !-------------------------------------------------------------------------------

   ! Function to check if a wavenumber is in the visible or IR
   logical function is_visible(wavenumber)

      ! Input wavenumber; this needs to be input in inverse cm (cm^-1)
      real(r8), intent(in) :: wavenumber

      ! Threshold between visible and infrared is 0.7 micron, or 14286 cm^-1
      real(r8), parameter :: visible_wavenumber_threshold = 14286._r8  ! cm^-1

      ! Wavenumber is in the visible if it is above the visible threshold
      ! wavenumber, and in the infrared if it is below the threshold
      if (wavenumber > visible_wavenumber_threshold) then
         is_visible = .true.
      else
         is_visible = .false.
      end if

   end function is_visible

   !-------------------------------------------------------------------------------

   ! Send shortwave fluxes and heating rates to history buffer
   subroutine output_fluxes_sw(icall, state, flux_all, flux_clr, qrs, qrsc)
      use physconst, only: cpair
      use physics_types, only: physics_state
      use cam_history, only: outfld
      
      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(fluxes_t), intent(in) :: flux_all
      type(fluxes_t), intent(in) :: flux_clr
      real(r8), intent(in) :: qrs(:,:), qrsc(:,:)

      ! SW cloud radiative effect
      real(r8) :: cloud_radiative_effect(state%ncol)

      ! Working variables
      integer :: ncol
      integer :: ktop_rad = 1

      ncol = state%ncol

      ! All-sky flux_all%fluxes at model interfaces
      call outfld('FDS'//diag(icall), flux_all%flux_dn(1:ncol,ktop:kbot+1), ncol, state%lchnk)
      call outfld('FUS'//diag(icall), flux_all%flux_up(1:ncol,ktop:kbot+1), ncol, state%lchnk)
      call outfld('FNS'//diag(icall), flux_all%flux_net(1:ncol,ktop:kbot+1), ncol, state%lchnk)
      call outfld('FDS_DIR'//diag(icall), flux_all%flux_dn_dir(1:ncol,ktop:kbot+1), ncol, state%lchnk)

      ! Clear-sky fluxes at model interfaces
      call outfld('FDSC'//diag(icall), flux_clr%flux_dn(1:ncol,ktop:kbot+1), ncol, state%lchnk)
      call outfld('FUSC'//diag(icall), flux_clr%flux_up(1:ncol,ktop:kbot+1), ncol, state%lchnk)
      call outfld('FNSC'//diag(icall), flux_clr%flux_net(1:ncol,ktop:kbot+1), ncol, state%lchnk)
      call outfld('FDSC_DIR'//diag(icall), flux_clr%flux_dn_dir(1:ncol,ktop:kbot+1), ncol, state%lchnk)

      ! All-sky fluxes
      call outfld('FSNT'//diag(icall), flux_all%flux_net(1:ncol,ktop), ncol, state%lchnk)
      call outfld('FSNS'//diag(icall), flux_all%flux_net(1:ncol,kbot+1), ncol, state%lchnk)
      call outfld('FSDS'//diag(icall), flux_all%flux_dn(1:ncol,kbot+1), ncol, state%lchnk)
      call outfld('FSUT'//diag(icall), flux_all%flux_up(1:ncol,ktop), ncol, state%lchnk)

      ! TOA fluxes (above model top, use index to rad top)
      call outfld('FSUTOA'//diag(icall), flux_all%flux_up(1:ncol,ktop_rad), ncol, state%lchnk)
      call outfld('FSNTOA'//diag(icall), flux_all%flux_net(1:ncol,ktop_rad), ncol, state%lchnk)

      ! Clear-sky fluxes
      call outfld('FSNTC'//diag(icall), flux_clr%flux_net(1:ncol,ktop), ncol, state%lchnk)
      call outfld('FSNSC'//diag(icall), flux_clr%flux_net(1:ncol,kbot+1), ncol, state%lchnk)
      call outfld('FSDSC'//diag(icall), flux_clr%flux_dn(1:ncol,kbot+1), ncol, state%lchnk)
      call outfld('FSUTOAC'//diag(icall), flux_clr%flux_up(1:ncol,ktop_rad), ncol, state%lchnk)
      call outfld('FSNTOAC'//diag(icall), flux_clr%flux_net(1:ncol,ktop_rad), ncol, state%lchnk)
      call outfld('FSUTC'//diag(icall), flux_clr%flux_up(1:ncol,ktop), ncol, state%lchnk)

      ! Calculate and output the shortwave cloud radiative effect (SWCF in history)
      cloud_radiative_effect(1:ncol) = flux_all%flux_net(1:ncol,ktop_rad) - flux_clr%flux_net(1:ncol,ktop_rad)
      call outfld('SWCF'//diag(icall), cloud_radiative_effect, ncol, state%lchnk)

      ! Send solar insolation to history buffer
      call outfld('SOLIN'//diag(icall), flux_clr%flux_dn(1:ncol,1), ncol, state%lchnk)
                        
      ! Send heating rates to history buffer
      call outfld('QRS'//diag(icall), qrs(1:ncol,1:pver)/cpair, ncol, state%lchnk)
      call outfld('QRSC'//diag(icall), qrsc(1:ncol,1:pver)/cpair, ncol, state%lchnk)

   end subroutine output_fluxes_sw


   ! Send longwave fluxes and heating rates to history buffer
   subroutine output_fluxes_lw(icall, state, flux_all, flux_clr, qrl, qrlc)
      use physconst, only: cpair
      use physics_types, only: physics_state
      use cam_history, only: outfld
      
      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(fluxes_t), intent(in) :: flux_all
      type(fluxes_t), intent(in) :: flux_clr

      ! Heating rates
      real(r8), intent(in) :: qrl(:,:), qrlc(:,:)

      ! Cloud radiative effect for output files
      real(r8) :: cloud_radiative_effect(state%ncol)

      ! Working arrays
      real(r8), dimension(pcols,pver+1) :: flux_up, flux_dn, flux_net
      integer :: ncol

      ncol = state%ncol

      ! Do all-sky fluxes, map first
      flux_up(1:ncol,1:pver+1) = flux_all%flux_up(1:ncol,ktop:kbot+1)
      flux_dn(1:ncol,1:pver+1) = flux_all%flux_dn(1:ncol,ktop:kbot+1)
      flux_net(1:ncol,1:pver+1) = flux_up(1:ncol,1:pver+1) - flux_dn(1:ncol,1:pver+1)

      ! All-sky fluxes at model interfaces
      call outfld('FDL'//diag(icall), flux_dn(1:ncol,1:pver+1), ncol, state%lchnk)
      call outfld('FUL'//diag(icall), flux_up(1:ncol,1:pver+1), ncol, state%lchnk)
      call outfld('FNL'//diag(icall), flux_net(1:ncol,1:pver+1), ncol, state%lchnk)

      ! Clear-sky fluxes at model interfaces
      flux_up(1:ncol,1:pver+1) = flux_clr%flux_up(1:ncol,ktop:kbot+1)
      flux_dn(1:ncol,1:pver+1) = flux_clr%flux_dn(1:ncol,ktop:kbot+1)
      flux_net(1:ncol,1:pver+1) = flux_up(1:ncol,1:pver+1) - flux_dn(1:ncol,1:pver+1)
      call outfld('FDLC'//diag(icall), flux_dn(1:ncol,1:pver+1), ncol, state%lchnk)
      call outfld('FULC'//diag(icall), flux_up(1:ncol,1:pver+1), ncol, state%lchnk)
      call outfld('FNLC'//diag(icall), flux_net(1:ncol,1:pver+1), ncol, state%lchnk)

      ! All-sky fluxes
      ! NOTE: sign change on net fluxes because internally net fluxes are assumed
      ! to down minus up, but for longwave outputs we want upward positive
      call outfld('FLNT'//diag(icall), -flux_all%flux_net(1:ncol,ktop), ncol, state%lchnk)
      call outfld('FLNS'//diag(icall), -flux_all%flux_net(1:ncol,kbot+1), ncol, state%lchnk)
      call outfld('FLUT'//diag(icall), flux_all%flux_up(1:ncol,ktop), ncol, state%lchnk)
      call outfld('FLDS'//diag(icall), flux_all%flux_dn(1:ncol,kbot+1), ncol, state%lchnk)

      ! Clear-sky fluxes
      ! NOTE: sign change on net fluxes because internally net fluxes are assumed
      ! to down minus up, but for longwave outputs we want upward positive
      call outfld('FLNTC'//diag(icall), -flux_clr%flux_net(1:ncol,ktop), ncol, state%lchnk)
      call outfld('FLNSC'//diag(icall), -flux_clr%flux_net(1:ncol,kbot+1), ncol, state%lchnk)
      call outfld('FLUTC'//diag(icall), flux_clr%flux_up(1:ncol,ktop), ncol, state%lchnk)
      call outfld('FLDSC'//diag(icall), flux_clr%flux_dn(1:ncol,kbot+1), ncol, state%lchnk)

      ! Calculate and output the cloud radiative effect (LWCF in history)
      cloud_radiative_effect(1:ncol) = flux_all%flux_net(1:ncol,ktop) - flux_clr%flux_net(1:ncol,ktop)
      call outfld('LWCF'//diag(icall), cloud_radiative_effect, ncol, state%lchnk)

      ! Send heating rates to history buffer
      call outfld('QRL'//diag(icall), qrl(1:ncol,1:pver)/cpair, ncol, state%lchnk)
      call outfld('QRLC'//diag(icall), qrlc(1:ncol,1:pver)/cpair, ncol, state%lchnk)

   end subroutine output_fluxes_lw

   !----------------------------------------------------------------------------

   ! Utility function to reorder an array given a new indexing
   function reordered(array_in, new_indexing) result(array_out)

      ! Inputs
      real(r8), intent(in) :: array_in(:)
      integer, intent(in) :: new_indexing(:)

      ! Output, reordered array
      real(r8), dimension(size(array_in)) :: array_out

      ! Loop index
      integer :: ii

      ! Check inputs
      call assert(size(array_in) == size(new_indexing), 'reorder_array: sizes inconsistent')

      ! Reorder array based on input index mapping, which maps old indices to new
      do ii = 1,size(new_indexing)
         array_out(ii) = array_in(new_indexing(ii))
      end do

   end function reordered

   !----------------------------------------------------------------------------

   subroutine output_cloud_optics_sw(state, tau, ssa, asm)
      use ppgrid, only: pver
      use physics_types, only: physics_state
      use cam_history, only: outfld
      use radconstants, only: idx_sw_diag

      type(physics_state), intent(in) :: state
      real(r8), intent(in), dimension(:,:,:) :: tau, ssa, asm
      integer :: sw_band_index
      character(len=*), parameter :: subname = 'output_cloud_optics_sw'

      ! Check values
      call assert_valid(tau(1:state%ncol,1:pver,1:nswbands), &
                        trim(subname) // ': optics%optical_depth')
      call assert_valid(ssa(1:state%ncol,1:pver,1:nswbands), &
                        trim(subname) // ': optics%single_scattering_albedo')
      call assert_valid(asm(1:state%ncol,1:pver,1:nswbands), &
                        trim(subname) // ': optics%assymmetry_parameter')

      ! Send outputs to history buffer
      call outfld('CLOUD_TAU_SW', &
                  tau(1:state%ncol,1:pver,1:nswbands), &
                  state%ncol, state%lchnk)
      call outfld('CLOUD_SSA_SW', &
                  ssa(1:state%ncol,1:pver,1:nswbands), &
                  state%ncol, state%lchnk)
      call outfld('CLOUD_G_SW', &
                  asm(1:state%ncol,1:pver,1:nswbands), &
                  state%ncol, state%lchnk)
      sw_band_index = get_band_index_sw(550._r8, 'nm')
      call outfld('TOT_ICLD_VISTAU', &
                  tau(1:state%ncol,1:pver,sw_band_index), &
                  state%ncol, state%lchnk)
   end subroutine output_cloud_optics_sw

   !----------------------------------------------------------------------------

   subroutine output_cloud_optics_lw(state, tau)

      use ppgrid, only: pver
      use physics_types, only: physics_state
      use cam_history, only: outfld

      type(physics_state), intent(in) :: state
      real(r8), intent(in), dimension(:,:,:) :: tau

      ! Check values
      call assert_valid(tau(1:state%ncol,1:pver,1:nlwbands), 'cld_tau_lw')

      ! Output
      call outfld('CLOUD_TAU_LW', &
                  tau(1:state%ncol,1:pver,1:nlwbands), &
                  state%ncol, state%lchnk)

   end subroutine output_cloud_optics_lw

   !----------------------------------------------------------------------------

   ! Should we do snow optics? Check for existence of "cldfsnow" variable
   logical function do_snow_optics()
      use physics_buffer, only: physics_buffer_desc, pbuf_get_index
      use phys_control, only: phys_getopts
      use cam_abortutils, only: endrun
      real(r8), pointer :: pbuf(:)
      integer :: err, idx
      logical :: use_MMF
      character(len=16) :: MMF_microphysics_scheme

      idx = pbuf_get_index('CLDFSNOW', errcode=err)
      if (idx > 0) then
         do_snow_optics = .true.
      else
         do_snow_optics = .false.
      end if

      ! Reset to false if using MMF with 1-mom scheme
      call phys_getopts(use_MMF_out           = use_MMF          )
      call phys_getopts(MMF_microphysics_scheme_out = MMF_microphysics_scheme)
      if (use_MMF .and. (trim(MMF_microphysics_scheme) == 'sam1mom')) then
         do_snow_optics = .false.
      end if

      return
   end function do_snow_optics

   !----------------------------------------------------------------------------

end module radiation
