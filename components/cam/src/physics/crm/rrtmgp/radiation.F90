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
   use shr_kind_mod,     only: r8=>shr_kind_r8, cl=>shr_kind_cl
   use ppgrid,           only: pcols, pver, pverp, begchunk, endchunk
   use cam_abortutils,   only: endrun
   use scamMod,          only: scm_crm_mode, single_column, swrad_off
   use rad_constituents, only: N_DIAG

   ! RRTMGP gas optics object to store coefficient information. This is imported
   ! here so that we can make the k_dist objects module data and only load them
   ! once.
   use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp

   ! Use my assertion routines to perform sanity checks
   use assertions, only: assert, assert_valid, assert_range

   use radiation_state, only: ktop, kbot, nlev_rad
   use radiation_utils, only: compress_day_columns, expand_day_columns, &
                              handle_error

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

   ! k-distribution coefficients. These will be populated by reading from the
   ! RRTMGP coefficients files, specified by coefficients_file_sw and
   ! coefficients_file_lw in the radiation namelist. They exist as module data
   ! because we only want to load those files once.
   type(ty_gas_optics_rrtmgp) :: k_dist_sw, k_dist_lw

   ! k-distribution coefficients files to read from. These are set via namelist
   ! variables.
   character(len=cl) :: coefficients_file_sw, coefficients_file_lw

   ! Set name for this module (for purpose of writing output and log files)
   character(len=*), parameter :: module_name = 'radiation'

   ! Number of shortwave and longwave bands in use by the RRTMGP radiation code.
   ! This information will be stored in the k_dist_sw and k_dist_lw objects and
   ! may
   ! be retrieved using the k_dist_sw%get_nband() and k_dist_lw%get_nband()
   ! methods, but I think we need to save these as private module data so that
   ! we
   ! can automatically allocate arrays later in subroutine headers, i.e.:
   !
   !     real(r8) :: cld_tau(pcols,pver,nswbands)
   !
   ! and so forth. Previously some of this existed in radconstants.F90, but I do
   ! not think we need to use that.
   ! EDIT: maybe these JUST below in radconstants.F90?
   integer :: nswbands, nlwbands

   ! Also, save number of g-points as private module data
   integer :: nswgpts, nlwgpts

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
      character(len=cl) :: rrtmgp_coefficients_file_lw, rrtmgp_coefficients_file_sw

      ! Variables defined in namelist
      namelist /radiation_nl/ rrtmgp_coefficients_file_lw,     &
                              rrtmgp_coefficients_file_sw,     &
                              iradsw, iradlw, irad_always,     &
                              use_rad_dt_cosz, spectralflux,   &
                              do_aerosol_rad, fixed_total_solar_irradiance

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
#endif

      ! Set module data
      coefficients_file_lw = rrtmgp_coefficients_file_lw
      coefficients_file_sw = rrtmgp_coefficients_file_sw

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
         write(iulog,10) trim(coefficients_file_lw), trim(coefficients_file_sw), &
                         iradsw, iradlw, irad_always, &
                         use_rad_dt_cosz, spectralflux, &
                         do_aerosol_rad, fixed_total_solar_irradiance
      end if
   10 format('  LW coefficents file: ',                                a/, &
             '  SW coefficents file: ',                                a/, &
             '  Frequency (timesteps) of Shortwave Radiation calc:  ',i5/, &
             '  Frequency (timesteps) of Longwave Radiation calc:   ',i5/, &
             '  SW/LW calc done every timestep for first N steps. N=',i5/, &
             '  Use average zenith angle:                           ',l5/, &
             '  Output spectrally resolved fluxes:                  ',l5/, &
             '  Do aerosol radiative calculations:                  ',l5/, &
             '  Fixed solar consant (disabled with -1):             ',f10.4/ )

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
            radiation_do = nstep == 0 .or. iradsw == 1                     &
                          .or. (mod(nstep-1,iradsw) == 0 .and. nstep /= 1) &
                          .or. nstep <= irad_always
         case ('lw') ! do a longwave heating calc this timestep?
            radiation_do = nstep == 0 .or. iradlw == 1                     &
                          .or. (mod(nstep-1,iradlw) == 0 .and. nstep /= 1) &
                          .or. nstep <= irad_always
         case ('aeres') ! write absorptivity/emissivity to restart file this timestep?
            ! for RRTMG there is no abs/ems restart file
            radiation_do = .false.
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
      do while (.not. dosw)
         nstep = nstep + 1
         offset = offset + dtime
         if (radiation_do('sw', nstep)) then
            radiation_nextsw_cday = get_curr_calday(offset=offset) 
            dosw = .true.
         end if
      end do
      if(radiation_nextsw_cday == -1._r8) then
         call endrun('error in radiation_nextsw_cday')
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
      use cam_history_support, only: add_hist_coord
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

      ! RRTMGP modules
      use rrtmgp_coefficients, only: rrtmgp_load_coefficients=>load_and_init
      use mo_gas_concentrations, only: ty_gas_concs

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

      logical :: use_MMF  ! MMF flag

      character(len=128) :: error_message

      ! ty_gas_concs object that would normally hold volume mixing ratios for
      ! radiatively-important gases. Here, this is just used to provide the names
      ! of gases that are available in the model (needed by the kdist
      ! initialization routines that are called within the load_coefficients
      ! methods).
      type(ty_gas_concs) :: available_gases

      real(r8), target :: sw_band_midpoints(nswbands), lw_band_midpoints(nlwbands)
      character(len=32) :: subname = 'radiation_init'

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

      ! Read gas optics coefficients from file
      ! Need to initialize available_gases here! The only field of the
      ! available_gases type that is used int he kdist initialize is
      ! available_gases%gas_name, which gives the name of each gas that would be
      ! present in the ty_gas_concs object. So, we can just set this here, rather
      ! than trying to fully populate the ty_gas_concs object here, which would be
      ! impossible from this initialization routine because I do not thing the
      ! rad_cnst objects are setup yet.
      call set_available_gases(active_gases, available_gases)
      call rrtmgp_load_coefficients(k_dist_sw, coefficients_file_sw, available_gases)
      call rrtmgp_load_coefficients(k_dist_lw, coefficients_file_lw, available_gases)

      ! Get number of bands used in shortwave and longwave and set module data
      ! appropriately so that these sizes can be used to allocate array sizes.
      nswbands = k_dist_sw%get_nband()
      nlwbands = k_dist_lw%get_nband()

      ! Set number of g-points for used for correlated-k. These are determined
      ! by the absorption coefficient data.
      nswgpts = k_dist_sw%get_ngpt()
      nlwgpts = k_dist_lw%get_ngpt()

      ! Set number of levels used in radiation calculations
#ifdef NO_EXTRA_RAD_LEVEL
      nlev_rad = pver
#else
      nlev_rad = pver + 1
#endif

      ! Indices on radiation grid that correspond to top and bottom of the model
      ! grid...this is for cases where we want to add an extra layer above the
      ! model top to deal with radiative heating above the model top...why???
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
      !if (docosp) call cospsimulator_intr_init()
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
      sw_band_midpoints(:) = get_band_midpoints(nswbands, k_dist_sw)
      lw_band_midpoints(:) = get_band_midpoints(nlwbands, k_dist_lw)
      call assert(all(sw_band_midpoints > 0), subname // ': negative sw_band_midpoints')
      call add_hist_coord('swband', nswbands, 'Shortwave band', 'wavelength', sw_band_midpoints)
      call add_hist_coord('lwband', nlwbands, 'Longwave band', 'wavelength', lw_band_midpoints)

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
                  'Cloud shortwave extinction optical depth', sampling_seq='rad_lwsw') 
      call addfld('CLOUD_SSA_SW', (/'lev   ','swband'/), 'I', '1', &
                  'Cloud shortwave single scattering albedo', sampling_seq='rad_lwsw') 
      call addfld('CLOUD_G_SW', (/'lev   ','swband'/), 'I', '1', &
                  'Cloud shortwave assymmetry parameter', sampling_seq='rad_lwsw') 
      call addfld('CLOUD_TAU_LW', (/'lev   ','lwband'/), 'I', '1', &
                  'Cloud longwave absorption optical depth', sampling_seq='rad_lwsw') 

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
                        sampling_seq='rad_lwsw')
            call addfld('SOLL'//diag(icall),  horiz_only, 'A',   'W/m2', &
                        'Solar downward near infrared direct to surface', &
                        sampling_seq='rad_lwsw')
            call addfld('SOLS'//diag(icall),  horiz_only, 'A',   'W/m2', &
                        'Solar downward visible direct to surface', &
                        sampling_seq='rad_lwsw')
            call addfld('SOLLD'//diag(icall), horiz_only, 'A',   'W/m2', &
                        'Solar downward near infrared diffuse to surface', &
                        sampling_seq='rad_lwsw')
            call addfld('SOLSD'//diag(icall), horiz_only, 'A',   'W/m2', &
                        'Solar downward visible diffuse to surface', &
                        sampling_seq='rad_lwsw')
            call addfld('QRS'//diag(icall),   (/ 'lev' /), 'A', 'K/s', &
                        'Solar heating rate', &
                        sampling_seq='rad_lwsw')
            call addfld('QRSC'//diag(icall),  (/ 'lev' /), 'A', 'K/s', &
                        'Clearsky solar heating rate', &
                        sampling_seq='rad_lwsw')
            call addfld('FSNS'//diag(icall),   horiz_only,  'A', 'W/m2', &
                        'Net solar flux at surface', &
                        sampling_seq='rad_lwsw')
            call addfld('FSNT'//diag(icall),   horiz_only,  'A', 'W/m2', &
                        'Net solar flux at top of model', &
                        sampling_seq='rad_lwsw')
            call addfld('FSNTOA'//diag(icall), horiz_only,  'A', 'W/m2', &
                        'Net solar flux at top of atmosphere', &
                        sampling_seq='rad_lwsw')
            call addfld('FSUTOA'//diag(icall), horiz_only,  'A',  'W/m2', &
                        'Upwelling solar flux at top of atmosphere', &
                        sampling_seq='rad_lwsw')
            call addfld('FSUT'//diag(icall), horiz_only,  'A',  'W/m2', &
                        'Upwelling solar flux at top of model', &
                        sampling_seq='rad_lwsw')
            call addfld('FSUTC'//diag(icall), horiz_only,  'A',  'W/m2', &
                        'Clearsky upwelling solar flux at top of model', &
                        sampling_seq='rad_lwsw')
            call addfld('FSNTOAC'//diag(icall), horiz_only, 'A',  'W/m2', &
                        'Clearsky net solar flux at top of atmosphere', &
                        sampling_seq='rad_lwsw')
            call addfld('FSUTOAC'//diag(icall), horiz_only, 'A',  'W/m2', &
                        'Clearsky upwelling solar flux at top of atmosphere', &
                        sampling_seq='rad_lwsw')
            call addfld('FSN200'//diag(icall), horiz_only,  'A',  'W/m2', &
                        'Net shortwave flux at 200 mb', &
                        sampling_seq='rad_lwsw')
            call addfld('FSN200C'//diag(icall), horiz_only, 'A',  'W/m2', &
                        'Clearsky net shortwave flux at 200 mb', &
                        sampling_seq='rad_lwsw')
            call addfld('FSNTC'//diag(icall), horiz_only,   'A',  'W/m2', &
                        'Clearsky net solar flux at top of model', &
                        sampling_seq='rad_lwsw')
            call addfld('FSNSC'//diag(icall), horiz_only,   'A',  'W/m2', &
                        'Clearsky net solar flux at surface', &
                        sampling_seq='rad_lwsw')
            call addfld('FSDSC'//diag(icall), horiz_only,   'A',  'W/m2', &
                        'Clearsky downwelling solar flux at surface', &
                        sampling_seq='rad_lwsw')
            call addfld('FSDS'//diag(icall), horiz_only,    'A',  'W/m2', &
                        'Downwelling solar flux at surface', &
                        sampling_seq='rad_lwsw')
            call addfld('FUS'//diag(icall),  (/ 'ilev' /),  'A',  'W/m2', &
                        'Shortwave upward flux')
            call addfld('FDS'//diag(icall),  (/ 'ilev' /),  'A',  'W/m2', &
                        'Shortwave downward flux')
            call addfld('FDS_DIR'//diag(icall),  (/ 'ilev' /),  'A',  'W/m2', &
                        'Shortwave direct-beam downward flux')
            call addfld('FNS'//diag(icall),  (/ 'ilev' /),  'A',  'W/m2', &
                        'Shortwave net flux')
            call addfld('FUSC'//diag(icall),  (/ 'ilev' /), 'A',  'W/m2', &
                        'Shortwave clear-sky upward flux')
            call addfld('FDSC'//diag(icall),  (/ 'ilev' /), 'A',  'W/m2', &
                        'Shortwave clear-sky downward flux')
            call addfld('FDSC_DIR'//diag(icall),  (/ 'ilev' /), 'A',  'W/m2', &
                        'Shortwave clear-sky direct-beam downward flux')
            call addfld('FNSC'//diag(icall),  (/ 'ilev' /), 'A',  'W/m2', &
                        'Shortwave clear-sky net flux')
            call addfld('FSNIRTOA'//diag(icall), horiz_only, 'A', 'W/m2', &
                        'Net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
                        sampling_seq='rad_lwsw')
            call addfld('FSNRTOAC'//diag(icall), horiz_only, 'A', 'W/m2', &
                        'Clearsky net near-infrared flux (Nimbus-7 WFOV) at top of atmosphere', &
                        sampling_seq='rad_lwsw')
            call addfld('FSNRTOAS'//diag(icall), horiz_only, 'A', 'W/m2', &
                        'Net near-infrared flux (>= 0.7 microns) at top of atmosphere', &
                        sampling_seq='rad_lwsw')
            call addfld('SWCF'//diag(icall),     horiz_only, 'A', 'W/m2', &
                        'Shortwave cloud forcing', &
                        sampling_seq='rad_lwsw')

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
                  sampling_seq='rad_lwsw')

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
                        sampling_seq='rad_lwsw')
            call addfld('QRLC'//diag(icall),    (/'lev'/),  'A', 'K/s',  &
                        'Clearsky longwave heating rate',                &
                        sampling_seq='rad_lwsw')
            call addfld('FLDS'//diag(icall),    horiz_only, 'A', 'W/m2', &
                        'Downwelling longwave flux at surface',          &
                        sampling_seq='rad_lwsw')
            call addfld('FLDSC'//diag(icall),   horiz_only, 'A', 'W/m2', &
                        'Clearsky Downwelling longwave flux at surface', &
                        sampling_seq='rad_lwsw')
            call addfld('FLNS'//diag(icall),    horiz_only, 'A', 'W/m2', &
                        'Net longwave flux at surface',                  &
                        sampling_seq='rad_lwsw')
            call addfld('FLNT'//diag(icall),    horiz_only, 'A', 'W/m2', &
                        'Net longwave flux at top of model',             &
                        sampling_seq='rad_lwsw')
            call addfld('FLUT'//diag(icall),    horiz_only, 'A', 'W/m2', &
                        'Upwelling longwave flux at top of model',       &
                        sampling_seq='rad_lwsw')
            call addfld('FLUTC'//diag(icall),   horiz_only, 'A', 'W/m2', &
                        'Clearsky upwelling longwave flux at top of model', &
                        sampling_seq='rad_lwsw')
            call addfld('FLNTC'//diag(icall),   horiz_only, 'A', 'W/m2', &
                        'Clearsky net longwave flux at top of model',    &
                        sampling_seq='rad_lwsw')
            call addfld('LWCF'//diag(icall),    horiz_only, 'A', 'W/m2', &
                        'Longwave cloud forcing',                        &
                        sampling_seq='rad_lwsw')
            call addfld('FLN200'//diag(icall),  horiz_only, 'A', 'W/m2', &
                        'Net longwave flux at 200 mb',                   &
                        sampling_seq='rad_lwsw')
            call addfld('FLN200C'//diag(icall), horiz_only, 'A', 'W/m2', &
                        'Clearsky net longwave flux at 200 mb',          &
                        sampling_seq='rad_lwsw')
            call addfld('FLNSC'//diag(icall),   horiz_only, 'A', 'W/m2', &
                        'Clearsky net longwave flux at surface',         &
                        sampling_seq='rad_lwsw')
            call addfld('FUL'//diag(icall),     (/'ilev'/), 'A', 'W/m2', &
                        'Longwave upward flux')
            call addfld('FDL'//diag(icall),     (/'ilev'/), 'A', 'W/m2', &
                        'Longwave downward flux')
            call addfld('FNL'//diag(icall),     (/'ilev'/), 'A', 'W/m2', &
                        'Longwave net flux')
            call addfld('FULC'//diag(icall),    (/'ilev'/), 'A', 'W/m2', &
                        'Longwave clear-sky upward flux')
            call addfld('FDLC'//diag(icall),    (/'ilev'/), 'A', 'W/m2', &
                        'Longwave clear-sky downward flux')
            call addfld('FNLC'//diag(icall),    (/'ilev'/), 'A', 'W/m2', &
                        'Longwave clear-sky net flux')

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
      call phys_getopts(use_MMF_out=use_MMF)
      if (use_MMF) then
         call addfld ('CRM_QRAD', (/'crm_nx_rad','crm_ny_rad','crm_nz    '/), 'A', 'K/s', &
                      'Radiative heating tendency')
         call addfld ('CRM_QRS ', (/'crm_nx_rad','crm_ny_rad','crm_nz    '/), 'I', 'K/s', &
                      'CRM Shortwave radiative heating rate')
         call addfld ('CRM_QRSC', (/'crm_nx_rad','crm_ny_rad','crm_nz    '/), 'I', 'K/s', &
                      'CRM clear-sky shortwave radiative heating rate')
         call addfld ('CRM_QRL ', (/'crm_nx_rad','crm_ny_rad','crm_nz    '/), 'I', 'K/s', &
                      'CRM Longwave radiative heating rate' )
         call addfld ('CRM_QRLC', (/'crm_nx_rad','crm_ny_rad','crm_nz    '/), 'I', 'K/s', &
                      'CRM clear-sky longwave radiative heating rate' )
         call addfld ('CRM_CLD_RAD', (/'crm_nx_rad','crm_ny_rad','crm_nz    '/), 'I', 'fraction', &
                      'CRM cloud fraction')
      end if

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
         call addfld (hirsname(1),horiz_only,'A','K','HIRS CH2 infra-red brightness temperature')
         call addfld (hirsname(2),horiz_only,'A','K','HIRS CH4 infra-red brightness temperature')
         call addfld (hirsname(3),horiz_only,'A','K','HIRS CH6 infra-red brightness temperature')
         call addfld (hirsname(4),horiz_only,'A','K','HIRS CH8 infra-red brightness temperature')
         call addfld (hirsname(5),horiz_only,'A','K','HIRS CH10 infra-red brightness temperature')
         call addfld (hirsname(6),horiz_only,'A','K','HIRS CH11 infra-red brightness temperature')
         call addfld (hirsname(7),horiz_only,'A','K','HIRS CH12 infra-red brightness temperature')
         call addfld (msuname(1),horiz_only,'A','K','MSU CH1 microwave brightness temperature')
         call addfld (msuname(2),horiz_only,'A','K','MSU CH2 microwave brightness temperature')
         call addfld (msuname(3),horiz_only,'A','K','MSU CH3 microwave brightness temperature')
         call addfld (msuname(4),horiz_only,'A','K','MSU CH4 microwave brightness temperature')

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
      call addfld ('HR',(/ 'lev' /), 'A','K/s','Heating rate needed for d(theta)/dt computation')

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


   subroutine perturbation_growth_init()

      use cam_logfile, only: iulog

      character(len=32) :: subname = 'perturbation_growth_init'

      ! Wrap modes in an ifdef for now since this is not implemented here, and I
      ! have some work to do to figure out what the heck this is meant to do.
#ifdef DO_PERGRO_MODS
      !Modification needed by pergro_mods for generating random numbers
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


   subroutine set_available_gases(gases, gas_concentrations)

      use mo_gas_concentrations, only: ty_gas_concs
      use mo_rrtmgp_util_string, only: lower_case

      type(ty_gas_concs), intent(inout) :: gas_concentrations
      character(len=*), intent(in) :: gases(:)
      character(len=32), dimension(size(gases)) :: gases_lowercase
      integer :: igas

      ! Initialize with lowercase gas names; we should work in lowercase
      ! whenever possible because we cannot trust string comparisons in RRTMGP
      ! to be case insensitive
      do igas = 1,size(gases)
         gases_lowercase(igas) = trim(lower_case(gases(igas)))
      end do
      call handle_error(gas_concentrations%init(gases_lowercase))

   end subroutine set_available_gases


   ! Function to calculate band midpoints from kdist band limits
   function get_band_midpoints(nbands, kdist) result(band_midpoints)
      integer, intent(in) :: nbands
      type(ty_gas_optics_rrtmgp), intent(in) :: kdist
      real(r8) :: band_midpoints(nbands)
      real(r8) :: band_limits(2,nbands)
      integer :: i
      character(len=32) :: subname = 'get_band_midpoints'

      call assert(kdist%get_nband() == nbands, trim(subname) // ': kdist%get_nband() /= nbands')

      band_limits = kdist%get_band_lims_wavelength()

      call assert(size(band_limits, 1) == size(kdist%get_band_lims_wavelength(), 1), &
                  subname // ': band_limits and kdist inconsistently sized')
      call assert(size(band_limits, 2) == size(kdist%get_band_lims_wavelength(), 2), &
                  subname // ': band_limits and kdist inconsistently sized')
      call assert(size(band_limits, 2) == size(band_midpoints), &
                  subname // ': band_limits and band_midpoints inconsistently sized')
      call assert(all(band_limits > 0), subname // ': negative band limit wavelengths!')

      band_midpoints(:) = 0._r8
      do i = 1,nbands
         band_midpoints(i) = (band_limits(1,i) + band_limits(2,i)) / 2._r8
      end do
      call assert(all(band_midpoints > 0), subname // ': negative band_midpoints!')

   end function get_band_midpoints


   !===============================================================================
        
   !----------------------------------------------------------------------------
   ! Main driver for radiation computation. Calls the shortwave and longwave
   ! drivers, and calculates the radiative heating from the resulting fluxes.
   ! Primary output from this routine is the heating tendency due to radiative
   ! transfer, as a ptend object.
   subroutine radiation_tend(state_in,ptend,    pbuf,          cam_out, cam_in,  &
                             landfrac,landm,    icefrac,       snowh,            &
                             fsns,    fsnt,     flns,          flnt,             &
                             fsds,    net_flux, is_cmip6_volc                    )

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

      ! RRTMGP radiation drivers and derived types
      use mo_gas_concentrations, only: ty_gas_concs
      use mo_optical_props, only: ty_optical_props_1scl, & 
                                  ty_optical_props_2str
      use mo_fluxes_byband, only: ty_fluxes_byband
      use mo_rrtmgp_util_string, only: lower_case

      ! CAM history module provides subroutine to send output data to the history
      ! buffer to be aggregated and written to disk
      use cam_history, only: outfld

      use pkg_cldoptics, only: cldefr  ! for sam1mom microphysics

      use physconst,    only: gravit, cpair
      use phys_control, only: phys_getopts

      use radiation_state, only: set_rad_state
      use radiation_utils, only: calculate_heating_rate
      use cam_optics, only: free_optics_sw       , free_optics_lw       , &
                            set_cloud_optics_lw  , set_cloud_optics_sw  , &
                            set_aerosol_optics_lw, set_aerosol_optics_sw

#ifdef MODAL_AERO
      use modal_aero_data, only: ntot_amode
#endif


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

      ! These are not used anymore and exist only because the radiation call is
      ! inflexible
      real(r8), intent(in)    :: landfrac(pcols)  ! land fraction
      real(r8), intent(in)    :: landm(pcols)     ! land fraction ramp
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

      ! Need to set effective radii for single-moment microphysics optics
      real(r8), pointer :: rel(:,:), rei(:,:)

      ! Options for MMF/SP
      logical :: use_MMF
      character(len=16) :: MMF_microphysics_scheme

      ! Flag to carry (QRS,QRL)*dp across time steps. 
      ! TODO: what does this mean?
      logical :: conserve_energy = .true.

      ! Number of columns; ncol is number of GCM columns in current chunk,
      ! ncol_tot is product of GCM and CRM columns (for packing GCM and CRM
      ! columns into a single dimension)
      integer :: ncol, ncol_tot

      ! Loop indices
      integer :: icall, ic, ix, iy, iz, ilev, j, igas, iband, igpt

      ! Everyone needs a name
      character(*), parameter :: subroutine_name = 'radiation_tend'

      ! Radiative fluxes
      type(ty_fluxes_byband) :: fluxes_allsky    , fluxes_clrsky, &
                                fluxes_allsky_all, fluxes_clrsky_all

      ! Copy of state
      type(physics_state) :: state

      ! For loops over diagnostic calls (TODO: what does this mean?)
      logical :: active_calls(0:N_DIAG)

      ! State fields that are passed into RRTMGP. Some of these may need to
      ! modified from what exists in the physics_state object, i.e. to clip
      ! temperatures to make sure they are within the valid range.
      real(r8), dimension(pcols * crm_nx_rad * crm_ny_rad, nlev_rad) :: tmid, pmid
      real(r8), dimension(pcols * crm_nx_rad * crm_ny_rad, nlev_rad+1) :: tint, pint
      real(r8), dimension(pcols, nlev_rad) :: tmid_col, pmid_col
      real(r8), dimension(pcols, nlev_rad+1) :: tint_col, pint_col

      ! Surface emissivity needed for longwave
      real(r8) :: surface_emissivity(nlwbands,pcols * crm_nx_rad * crm_ny_rad)

      ! Temporary heating rates on radiation vertical grid
      real(r8), dimension(pcols * crm_nx_rad * crm_ny_rad,pver) :: qrl_all, qrlc_all

      ! Pointers to CRM fields on physics buffer
      real(r8), pointer, dimension(:,:) :: iclwp, iciwp, cld
      real(r8), pointer, dimension(:,:,:,:) :: crm_t, crm_qv, crm_qc, crm_qi, crm_cld, crm_qrad
      real(r8), pointer, dimension(:,:,:,:,:) :: crm_qaerwat, crm_dgnumwet
      real(r8), dimension(pcols,pver) :: iclwp_save, iciwp_save, cld_save
      integer :: ixwatvap, ixcldliq, ixcldice

      real(r8), pointer, dimension(:,:,:) :: qaerwat, dgnumwet
#ifdef MODAL_AERO
      real(r8), dimension(pcols,pver,ntot_amode) :: qaerwat_save, dgnumwet_save
#endif
      real(r8), pointer :: dei(:,:)
      real(r8), dimension(pcols,pver) :: dei_save, rei_save, rel_save

      ! Arrays to hold gas volume mixing ratios
      real(r8), dimension(size(active_gases),pcols,nlev_rad) :: vmr_col
      real(r8), dimension(size(active_gases),pcols * crm_nx_rad * crm_ny_rad,nlev_rad) :: vmr_all

      ! Clear-sky CRM heating (for debugging)
      real(r8), dimension(pcols,crm_nx_rad,crm_ny_rad,crm_nz) :: &
               crm_qrs, crm_qrsc, crm_qrl, crm_qrlc

      ! Albedo for shortwave calculations
      real(r8), dimension(nswbands,pcols) :: albedo_direct_col, albedo_diffuse_col, &
                                             albedo_direct_day, albedo_diffuse_day
      real(r8), dimension(nswbands,pcols * crm_nx_rad * crm_ny_rad) :: &
                                             albedo_direct_all, albedo_diffuse_all

      ! Cosine solar zenith angle for all columns in chunk
      real(r8) :: coszrs(pcols), coszrs_all(pcols * crm_nx_rad * crm_ny_rad)

      ! Incoming solar radiation, scaled for solar zenith angle 
      ! and earth-sun distance
      real(r8) :: solar_irradiance_by_gpt(pcols,nswgpts)

      ! Gathered indicies of day and night columns 
      ! chunk_column_index = day_indices(daylight_column_index)
      integer :: nday, nnight     ! Number of daylight columns
      integer :: day_indices(pcols), night_indices(pcols)   ! Indicies of daylight coumns


      ! Scaling factor for total sky irradiance; used to account for orbital
      ! eccentricity, and could be used to scale total sky irradiance for different
      ! climates as well (i.e., paleoclimate simulations)
      real(r8) :: tsi_scaling

      ! Cloud and aerosol optics
      type(ty_optical_props_2str) :: &
         aer_optics_sw_col, aer_optics_sw_all, &
         cld_optics_sw_col, cld_optics_sw_all
      type(ty_optical_props_1scl) :: &
         aer_optics_lw_col, aer_optics_lw_all, &
         cld_optics_lw_col, cld_optics_lw_all

      real(r8), dimension(pcols * crm_nx_rad * crm_ny_rad,pver) :: qrs_all, qrsc_all

      ! Area factor for calculating CRM domain averages
      real(r8), parameter :: area_factor = 1._r8 / (crm_nx_rad * crm_ny_rad)

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

      ! Set flags for MMF/SP
      call phys_getopts(use_MMF_out           = use_MMF          )
      call phys_getopts(MMF_microphysics_scheme_out = MMF_microphysics_scheme)

      ! Set pointers to heating rates stored on physics buffer. These will be
      ! modified in this routine.
      call pbuf_get_field(pbuf, pbuf_get_index('QRS'), qrs)
      call pbuf_get_field(pbuf, pbuf_get_index('QRL'), qrl)

      ! Set surface emissivity to 1 here. There is a note in the RRTMG
      ! implementation that this is treated in the land model, but the old
      ! RRTMG implementation also sets this to 1. This probably does not make
      ! a lot of difference either way, but if a more intelligent value
      ! exists or is assumed in the model we should use it here as well.
      surface_emissivity = 1.0_r8

      ! pbuf fields we need to overwrite with CRM fields to work with optics
      call pbuf_get_field(pbuf, pbuf_get_index('ICLWP'), iclwp)
      call pbuf_get_field(pbuf, pbuf_get_index('ICIWP'), iciwp)
      call pbuf_get_field(pbuf, pbuf_get_index('CLD'  ), cld  )
      call pbuf_get_field(pbuf, pbuf_get_index('DEI'  ), dei  )
      call pbuf_get_field(pbuf, pbuf_get_index('REL'), rel)
      call pbuf_get_field(pbuf, pbuf_get_index('REI'), rei)
#ifdef MODAL_AERO
      call pbuf_get_field(pbuf, pbuf_get_index('QAERWAT'), qaerwat)
      call pbuf_get_field(pbuf, pbuf_get_index('DGNUMWET'), dgnumwet)
#endif

      ! CRM fields
      call phys_getopts(use_MMF_out=use_MMF)
      if (use_MMF) then
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_T_RAD'  ), crm_t  )
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_QV_RAD' ), crm_qv )
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_QC_RAD' ), crm_qc )
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_QI_RAD' ), crm_qi )
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_CLD_RAD'), crm_cld)
#ifdef MODAL_AERO
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_QAERWAT' ), crm_qaerwat )
         call pbuf_get_field(pbuf, pbuf_get_index('CRM_DGNUMWET'), crm_dgnumwet)
#endif

         ! Output CRM cloud fraction
         call outfld('CRM_CLD_RAD', crm_cld, state%ncol, state%lchnk)
      end if

      ! Save pbuf things to restore when we are done working with them. This is
      ! needed because the CAM optics routines bury pbuf dependencies down deep
      ! in the call stack, so we need to overwrite with CRM state here in order
      ! to use CRM information to calculate optics. This can go away if we
      ! refactor the optics routines to take array arguments instead of the high
      ! level pbuf and state data structures.
      iclwp_save = iclwp
      iciwp_save = iciwp
      cld_save   = cld
      dei_save   = dei
      rei_save   = rei
      rel_save   = rel
#ifdef MODAL_AERO
      qaerwat_save = qaerwat
      dgnumwet_save = dgnumwet
#endif
 
      ! Indices into rad constituents arrays
      ixwatvap = 1
      call cnst_get_ind('CLDLIQ', ixcldliq)
      call cnst_get_ind('CLDICE', ixcldice)

      ! Loop over "diagnostic calls"; these are additional configurations of
      ! gases used to calculate radiative effects for diagnostic purposes only,
      ! without affecting the climate. The only version that affects the
      ! simulation is the first (index 0).
      call rad_cnst_get_call_list(active_calls)
      do icall = N_DIAG,0,-1
         if (active_calls(icall)) then

            ! Calculate derived quantities and pack data. This includes cloud
            ! optical properties. This is done separately from the blocks of
            ! code where shortwave and longwave fluxes are calculated so that we
            ! only have to pack/copy data once, rather than once for each of
            ! shortwave and longwave.
            if (radiation_do('sw') .or. radiation_do('lw')) then

               ! Calculate effective radius for optics used with single-moment microphysics
               if (use_MMF .and. (trim(MMF_microphysics_scheme) .eq. 'sam1mom')) then 
                  call cldefr(state%lchnk, ncol, landfrac, state%t, rel, rei, state%ps, state%pmid, landm, icefrac, snowh)
               end if

               ! Get albedo. This uses CAM routines internally and just provides a
               ! wrapper to improve readability of the code here.
               call set_albedo(cam_in, albedo_direct_col(1:nswbands,1:ncol), albedo_diffuse_col(1:nswbands,1:ncol))

               ! Send albedos to history buffer (useful for debugging)
               call outfld('SW_ALBEDO_DIR', transpose(albedo_direct_col(1:nswbands,1:ncol)), ncol, state%lchnk)
               call outfld('SW_ALBEDO_DIF', transpose(albedo_diffuse_col(1:nswbands,1:ncol)), ncol, state%lchnk)

               ! Get cosine solar zenith angle for current time step, and send to
               ! history buffer
               call set_cosine_solar_zenith_angle(state, dt_avg, coszrs(1:ncol))
               call outfld('COSZRS', coszrs(1:ncol), ncol, state%lchnk)

               ! If the swrad_off flag is set, meaning we should not do SW radiation, then 
               ! we just set coszrs to zero everywhere. TODO: why not just set dosw false 
               ! and skip the loop?
               if (swrad_off) coszrs(:) = 0._r8

               ! Gather night/day column indices for subsetting SW inputs; we only want to
               ! do the shortwave radiative transfer during the daytime to save
               ! computational cost (and because RRTMGP will fail for cosine solar zenith
               ! angles less than or equal to zero)
               call set_daynight_indices(coszrs(1:ncol), day_indices(1:ncol), night_indices(1:ncol))
               nday = count(day_indices(1:ncol) > 0)
               nnight = count(night_indices(1:ncol) > 0)

               ! Initialize cloud optics objects
               call handle_error(cld_optics_sw_col%alloc_2str( &
                  ncol, nlev_rad, k_dist_sw, name='cld_optics_sw' &
               ))
               call handle_error(cld_optics_sw_all%alloc_2str( &
                  ncol_tot, nlev_rad, k_dist_sw, name='cld_optics_sw' &
               ))
               call handle_error(cld_optics_lw_col%alloc_1scl( &
                  ncol, nlev_rad, k_dist_lw, name='cld_optics_lw' &
               ))
               call handle_error(cld_optics_lw_all%alloc_1scl( &
                  ncol_tot, nlev_rad, k_dist_lw, name='cld_optics_lw' &
               ))

               ! Initialize aerosol optics; passing only the wavenumber bounds for each
               ! "band" rather than passing the full spectral discretization object, and
               ! omitting the "g-point" mapping forces the optics to be indexed and
               ! stored by band rather than by g-point. This is most consistent with our
               ! treatment of aerosol optics in the model, and prevents us from having to
               ! map bands to g-points ourselves since that will all be handled by the
               ! private routines internal to the optics class.
               call handle_error(aer_optics_sw_all%alloc_2str(              &
                  ncol_tot, nlev_rad, k_dist_sw%get_band_lims_wavenumber(), &
                  name='aer_optics_sw_all'                                  &
               ))
               call handle_error(aer_optics_lw_all%alloc_1scl(          &
                  ncol_tot, nlev_rad, k_dist_lw%get_band_lims_wavenumber(), &
                  name='aer_optics_lw_all'                              &
               ))
               call handle_error(aer_optics_sw_col%alloc_2str(              &
                  ncol    , nlev_rad, k_dist_sw%get_band_lims_wavenumber(), &
                  name='aer_optics_sw_col'                                  &
               ))
               call handle_error(aer_optics_lw_col%alloc_1scl(          &
                  ncol    , nlev_rad, k_dist_lw%get_band_lims_wavenumber(), &
                  name='aer_optics_lw_col'                              &
               ))

               ! Loop over CRM columns; call routines designed to work with
               ! pbuf/state over ncol columns for each CRM column index, and pack
               ! into arrays dimensioned ncol_tot = ncol * ncrms
               j = 1
               do iy = 1,crm_ny_rad
                  do ix = 1,crm_nx_rad

                     ! Overwrite state and pbuf with CRM
                     if (use_MMF) then
                        call t_startf('rad_overwrite_state')
                        do iz = 1,crm_nz
                           ilev = pver - iz + 1
                           do ic = 1,ncol
                              ! NOTE: I do not think these are used by optics
                              state%q(ic,ilev,ixcldliq) = crm_qc(ic,ix,iy,iz)
                              state%q(ic,ilev,ixcldice) = crm_qi(ic,ix,iy,iz)
                              state%q(ic,ilev,ixwatvap) = crm_qv(ic,ix,iy,iz)
                              state%t(ic,ilev)          = crm_t (ic,ix,iy,iz)

                              ! In-cloud liquid and ice water paths (used by cloud optics)
                              cld(ic,ilev) = crm_cld(ic,ix,iy,iz)
                              if (cld(ic,ilev) > 0) then
                                 iclwp(ic,ilev) = crm_qc(ic,ix,iy,iz)          &
                                                * state%pdel(ic,ilev) / gravit &
                                                / max(0.01_r8, cld(ic,ilev))
                                 iciwp(ic,ilev) = crm_qi(ic,ix,iy,iz)          &
                                                * state%pdel(ic,ilev) / gravit &
                                                / max(0.01_r8, cld(ic,ilev))
                              else
                                 iclwp(ic,ilev) = 0
                                 iciwp(ic,ilev) = 0
                              end if
#ifdef MODAL_AERO
                              ! Use CRM scale aerosol water to calculate aerosol optical depth. Here we assume no
                              ! aerosol water uptake for cloudy sky on CRM grids. This is not really physically
                              ! correct, but if we assume 100% of relative humidity for aerosol water uptake, this
                              ! will bias 'AODVIS' to be large, since 'AODVIS' is used to compare with observed
                              ! clear sky AOD. In the future, AODVIS should be calculated from clear sky CRM AOD
                              ! only. But before this is done, we will assume no water uptake on CRM grids for
                              ! cloudy conditions (The radiative effects of this assumption will be small, since
                              ! aerosol effects are small relative to cloud effects for cloudy sky anyway.
                              ! -Minghuai Wang (minghuai.wang@pnl.gov)
                              qaerwat (ic,ilev,1:ntot_amode) =  crm_qaerwat(ic,ix,iy,iz,1:ntot_amode)
                              dgnumwet(ic,ilev,1:ntot_amode) = crm_dgnumwet(ic,ix,iy,iz,1:ntot_amode)
#endif
                              ! DEI is used for 2-moment optics
                              dei(1:ncol,1:pver) = 2._r8 * rei(1:ncol,1:pver)
                           end do  ! ic = 1,ncol
                        end do  ! iz = 1,crm_nz
                        call t_stopf('rad_overwrite_state')
                     end if  ! use_MMF

                     ! Setup state arrays, which may contain an extra level
                     ! above model top to handle heating above the model
                     call t_startf('rad_set_state')
                     call set_rad_state(                                            &
                        state                      , cam_in                       , &
                        tmid_col(1:ncol,1:nlev_rad), tint_col(1:ncol,1:nlev_rad+1), &
                        pmid_col(1:ncol,1:nlev_rad), pint_col(1:ncol,1:nlev_rad+1)  &
                     )
                     call t_stopf('rad_set_state')

                     ! Set gas concentrations
                     call t_startf('rad_gas_concentrations')
                     do igas = 1,size(active_gases)
                        ! Get volume mixing ratio for this gas
                        call get_gas_vmr(icall, state, pbuf, trim(active_gases(igas)), vmr_col(igas,1:ncol,ktop:kbot))

                        ! Copy top model level to level above model top
                        vmr_col(igas,1:ncol,1) = vmr_col(igas,1:ncol,ktop)
                     end do
                     call t_stopf('rad_gas_concentrations')

                     ! Do shortwave cloud optics calculations
                     ! TODO: refactor the set_cloud_optics codes to allow passing arrays
                     ! rather than state/pbuf so that we can use this for superparameterized
                     ! simulations...or alternatively add logic within the set_cloud_optics
                     ! routines to handle this.
                     if (radiation_do('sw')) then
                        call t_startf('rad_cloud_optics_sw')
                        call set_cloud_optics_sw(state, pbuf, k_dist_sw, cld_optics_sw_col)
                        call t_stopf('rad_cloud_optics_sw')

                        ! Get shortwave aerosol optics
                        call t_startf('rad_aerosol_optics_sw')
                        call set_aerosol_optics_sw(icall, state, pbuf, &
                                                   night_indices(1:nnight), &
                                                   is_cmip6_volc, &
                                                   aer_optics_sw_col)
                        call t_stopf('rad_aerosol_optics_sw')
                     end if

                     ! Do optics; TODO: refactor to take array arguments?
                     if (radiation_do('lw')) then
                        call t_startf('rad_cloud_optics_lw')
                        call set_cloud_optics_lw(state, pbuf, k_dist_lw, cld_optics_lw_col)
                        call t_stopf('rad_cloud_optics_lw')

                        call t_startf('rad_aerosol_optics_lw')
                        call set_aerosol_optics_lw(icall, state, pbuf, is_cmip6_volc, aer_optics_lw_col)
                        call t_stopf('rad_aerosol_optics_lw')
                     end if

                     ! Pack data
                     call t_startf('rad_pack_columns')
                     do ic = 1,ncol
                        coszrs_all(j) = coszrs(ic)
                        albedo_direct_all(:,j) = albedo_direct_col(:,ic)
                        albedo_diffuse_all(:,j) = albedo_diffuse_col(:,ic)
                        pmid(j,:) = pmid_col(ic,:)
                        tmid(j,:) = tmid_col(ic,:)
                        pint(j,:) = pint_col(ic,:)
                        tint(j,:) = tint_col(ic,:)
                        cld_optics_lw_all%tau(j,:,:) = cld_optics_lw_col%tau(ic,:,:)
                        aer_optics_lw_all%tau(j,:,:) = aer_optics_lw_col%tau(ic,:,:)
                        cld_optics_sw_all%tau(j,:,:) = cld_optics_sw_col%tau(ic,:,:)
                        cld_optics_sw_all%ssa(j,:,:) = cld_optics_sw_col%ssa(ic,:,:)
                        cld_optics_sw_all%g  (j,:,:) = cld_optics_sw_col%g  (ic,:,:)
                        aer_optics_sw_all%tau(j,:,:) = aer_optics_sw_col%tau(ic,:,:)
                        aer_optics_sw_all%ssa(j,:,:) = aer_optics_sw_col%ssa(ic,:,:)
                        aer_optics_sw_all%g  (j,:,:) = aer_optics_sw_col%g  (ic,:,:)
                        vmr_all(:,j,:) = vmr_col(:,ic,:)
                        j = j + 1
                     end do  ! ic = 1,ncol
                     call t_stopf('rad_pack_columns')

                  end do  ! ix = 1,crm_nx_rad
               end do  ! iy = 1,crm_ny_rad

            end if  !(radiation_do('sw') .or. radiation_do('lw')) then

            ! Do shortwave stuff...
            if (radiation_do('sw')) then

               ! Get orbital eccentricity factor to scale total sky irradiance
               tsi_scaling = get_eccentricity_factor()

               ! Allocate shortwave fluxes (allsky and clearsky)
               ! NOTE: fluxes defined at interfaces, so initialize to have vertical
               ! dimension nlev_rad+1, while we initialized the RRTMGP input variables to
               ! have vertical dimension nlev_rad (defined at midpoints).
               call initialize_rrtmgp_fluxes(ncol    , nlev_rad+1, nswbands, fluxes_allsky    , do_direct=.true.)
               call initialize_rrtmgp_fluxes(ncol    , nlev_rad+1, nswbands, fluxes_clrsky    , do_direct=.true.)
               call initialize_rrtmgp_fluxes(ncol_tot, nlev_rad+1, nswbands, fluxes_allsky_all, do_direct=.true.)
               call initialize_rrtmgp_fluxes(ncol_tot, nlev_rad+1, nswbands, fluxes_clrsky_all, do_direct=.true.)

               ! Do shortwave radiative transfer calculations
               call t_startf('rad_calculate_fluxes_sw')
               call calculate_fluxes_sw( &
                  active_gases(:), vmr_all(:,1:ncol_tot,1:nlev_rad), &
                  pmid(1:ncol_tot,1:nlev_rad), tmid(1:ncol_tot,1:nlev_rad), &
                  pint(1:ncol_tot,1:nlev_rad+1), &
                  coszrs_all(1:ncol_tot), &
                  albedo_direct_all(1:nswbands,1:ncol_tot), &
                  albedo_diffuse_all(1:nswbands,1:ncol_tot), &
                  cld_optics_sw_all, aer_optics_sw_all, &
                  fluxes_allsky_all, fluxes_clrsky_all, tsi_scaling &
               )
               call t_stopf('rad_calculate_fluxes_sw')

               ! Calculate heating rates
               call t_startf('rad_heating_rate_sw')
               call calculate_heating_rate(fluxes_allsky_all%flux_up(1:ncol_tot,ktop:kbot+1), &
                                           fluxes_allsky_all%flux_dn(1:ncol_tot,ktop:kbot+1), &
                                           pint(1:ncol_tot,ktop:kbot+1), &
                                           qrs_all(1:ncol_tot,1:pver))
               call calculate_heating_rate(fluxes_clrsky_all%flux_up(1:ncol_tot,ktop:kbot+1), &
                                           fluxes_clrsky_all%flux_dn(1:ncol_tot,ktop:kbot+1), &
                                           pint(1:ncol_tot,ktop:kbot+1), &
                                           qrsc_all(1:ncol_tot,1:pver))
               call t_stopf('rad_heating_rate_sw')

               ! Calculate CRM domain averages
               call t_startf('rad_average_fluxes_sw')
               call average_packed_array(qrs_all (1:ncol_tot,:)                     , qrs (1:ncol,:)                     )
               call average_packed_array(qrsc_all(1:ncol_tot,:)                     , qrsc(1:ncol,:)                     )
               call average_packed_array(fluxes_allsky_all%flux_up    (1:ncol_tot,:), fluxes_allsky%flux_up    (1:ncol,:))
               call average_packed_array(fluxes_allsky_all%flux_dn    (1:ncol_tot,:), fluxes_allsky%flux_dn    (1:ncol,:))
               call average_packed_array(fluxes_allsky_all%flux_net   (1:ncol_tot,:), fluxes_allsky%flux_net   (1:ncol,:))
               call average_packed_array(fluxes_allsky_all%flux_dn_dir(1:ncol_tot,:), fluxes_allsky%flux_dn_dir(1:ncol,:))
               do iband = 1,nswbands
                  call average_packed_array(fluxes_allsky_all%bnd_flux_up    (1:ncol_tot,:,iband), fluxes_allsky%bnd_flux_up    (1:ncol,:,iband))
                  call average_packed_array(fluxes_allsky_all%bnd_flux_dn    (1:ncol_tot,:,iband), fluxes_allsky%bnd_flux_dn    (1:ncol,:,iband))
                  call average_packed_array(fluxes_allsky_all%bnd_flux_net   (1:ncol_tot,:,iband), fluxes_allsky%bnd_flux_net   (1:ncol,:,iband))
                  call average_packed_array(fluxes_allsky_all%bnd_flux_dn_dir(1:ncol_tot,:,iband), fluxes_allsky%bnd_flux_dn_dir(1:ncol,:,iband))
               end do
               call average_packed_array(fluxes_clrsky_all%flux_up    (1:ncol_tot,:), fluxes_clrsky%flux_up    (1:ncol,:))
               call average_packed_array(fluxes_clrsky_all%flux_dn    (1:ncol_tot,:), fluxes_clrsky%flux_dn    (1:ncol,:))
               call average_packed_array(fluxes_clrsky_all%flux_net   (1:ncol_tot,:), fluxes_clrsky%flux_net   (1:ncol,:))
               call average_packed_array(fluxes_clrsky_all%flux_dn_dir(1:ncol_tot,:), fluxes_clrsky%flux_dn_dir(1:ncol,:))
               do iband = 1,nswbands
                  call average_packed_array(fluxes_clrsky_all%bnd_flux_up    (1:ncol_tot,:,iband), fluxes_clrsky%bnd_flux_up    (1:ncol,:,iband))
                  call average_packed_array(fluxes_clrsky_all%bnd_flux_dn    (1:ncol_tot,:,iband), fluxes_clrsky%bnd_flux_dn    (1:ncol,:,iband))
                  call average_packed_array(fluxes_clrsky_all%bnd_flux_net   (1:ncol_tot,:,iband), fluxes_clrsky%bnd_flux_net   (1:ncol,:,iband))
                  call average_packed_array(fluxes_clrsky_all%bnd_flux_dn_dir(1:ncol_tot,:,iband), fluxes_clrsky%bnd_flux_dn_dir(1:ncol,:,iband))
               end do
               call t_stopf('rad_average_fluxes_sw')

               ! Send fluxes to history buffer
               call output_fluxes_sw(icall, state, fluxes_allsky, fluxes_clrsky, qrs,  qrsc)

               ! Map to CRM columns
               if (use_MMF) then
                  j = 1
                  do iy = 1,crm_ny_rad
                     do ix = 1,crm_nx_rad
                        do ic = 1,ncol
                           do iz = 1,crm_nz
                              ilev = pver - iz + 1
                              crm_qrs (ic,ix,iy,iz) = qrs_all(j,ilev)
                              crm_qrsc(ic,ix,iy,iz) = qrsc_all(j,ilev)
                           end do
                           j = j + 1
                        end do
                     end do
                  end do
                  call outfld('CRM_QRS' , crm_qrs (1:ncol,:,:,:)/cpair, ncol, state%lchnk)
                  call outfld('CRM_QRSC', crm_qrsc(1:ncol,:,:,:)/cpair, ncol, state%lchnk)
               end if

               ! Set net fluxes used by other components (land?) 
               call set_net_fluxes_sw(fluxes_allsky, fsds, fsns, fsnt)

               ! Set surface fluxes that are used by the land model
               call export_surface_fluxes(fluxes_allsky, cam_out, 'shortwave')
               
               ! Free memory allocated for shortwave fluxes
               call free_fluxes(fluxes_allsky)
               call free_fluxes(fluxes_clrsky)
               call free_fluxes(fluxes_allsky_all)
               call free_fluxes(fluxes_clrsky_all)

               ! Free optical properties
               call free_optics_sw(cld_optics_sw_all)
               call free_optics_sw(aer_optics_sw_all)
               call free_optics_sw(cld_optics_sw_col)
               call free_optics_sw(aer_optics_sw_col)

            else

               ! Conserve energy
               if (conserve_energy) then
                  qrs(1:ncol,1:pver) = qrs(1:ncol,1:pver) / state%pdel(1:ncol,1:pver)
               end if

            end if  ! dosw

            ! Do longwave stuff...
            if (radiation_do('lw')) then

               ! Allocate longwave outputs; why is this not part of the
               ! ty_fluxes_byband object?
               ! NOTE: fluxes defined at interfaces, so initialize to have vertical
               ! dimension nlev_rad+1
               call initialize_rrtmgp_fluxes(ncol, nlev_rad+1, nlwbands, fluxes_allsky)
               call initialize_rrtmgp_fluxes(ncol, nlev_rad+1, nlwbands, fluxes_clrsky)
               call initialize_rrtmgp_fluxes(ncol_tot, nlev_rad+1, nlwbands, fluxes_allsky_all)
               call initialize_rrtmgp_fluxes(ncol_tot, nlev_rad+1, nlwbands, fluxes_clrsky_all)

               ! Calculate longwave fluxes
               call t_startf('rad_fluxes_lw')
               call calculate_fluxes_lw(                                          &
                  active_gases, vmr_all(:,1:ncol_tot,1:nlev_rad),               &
                  surface_emissivity(1:nlwbands,1:ncol_tot),                    &
                  pmid(1:ncol_tot,1:nlev_rad  ), tmid(1:ncol_tot,1:nlev_rad  ), &
                  pint(1:ncol_tot,1:nlev_rad+1), tint(1:ncol_tot,1:nlev_rad+1), &
                  cld_optics_lw_all            , aer_optics_lw_all,             &
                  fluxes_allsky_all            , fluxes_clrsky_all              &
               )
               call t_stopf('rad_fluxes_lw')

               ! Calculate heating rates
               call t_startf('rad_heating_lw')
               call calculate_heating_rate(fluxes_allsky_all%flux_up(:,ktop:kbot+1), &
                                           fluxes_allsky_all%flux_dn(:,ktop:kbot+1), &
                                           pint(1:ncol_tot,ktop:kbot+1)            , &
                                           qrl_all(1:ncol_tot,1:pver)                )
               call calculate_heating_rate(fluxes_clrsky_all%flux_up(:,ktop:kbot+1), &
                                           fluxes_clrsky_all%flux_dn(:,ktop:kbot+1), &
                                           pint(1:ncol_tot,ktop:kbot+1)            , &
                                           qrlc_all(1:ncol_tot,1:pver)               )
               call t_stopf('rad_heating_lw')

               ! Calculate domain averages
               call t_startf('rad_average_fluxes_lw')
               call average_packed_array(qrl_all (1:ncol_tot,:)                  , qrl (1:ncol,:))
               call average_packed_array(qrlc_all(1:ncol_tot,:)                  , qrlc(1:ncol,:))
               call average_packed_array(fluxes_allsky_all%flux_up (1:ncol_tot,:), fluxes_allsky%flux_up (1:ncol,:))
               call average_packed_array(fluxes_allsky_all%flux_dn (1:ncol_tot,:), fluxes_allsky%flux_dn (1:ncol,:))
               call average_packed_array(fluxes_allsky_all%flux_net(1:ncol_tot,:), fluxes_allsky%flux_net(1:ncol,:))
               call average_packed_array(fluxes_clrsky_all%flux_up (1:ncol_tot,:), fluxes_clrsky%flux_up (1:ncol,:))
               call average_packed_array(fluxes_clrsky_all%flux_dn (1:ncol_tot,:), fluxes_clrsky%flux_dn (1:ncol,:))
               call average_packed_array(fluxes_clrsky_all%flux_net(1:ncol_tot,:), fluxes_clrsky%flux_net(1:ncol,:))
               call t_stopf('rad_average_fluxes_lw')
                         
               ! Send fluxes to history buffer
               call output_fluxes_lw(icall, state, fluxes_allsky, fluxes_clrsky, qrl, qrlc)

               ! Map to CRM columns
               if (use_MMF) then
                  j = 1
                  do iy = 1,crm_ny_rad
                     do ix = 1,crm_nx_rad
                        do ic = 1,ncol
                           do iz = 1,crm_nz
                              ilev = pver - iz + 1
                              crm_qrl(ic,ix,iy,iz) = qrl_all(j,ilev)
                              crm_qrlc(ic,ix,iy,iz) = qrlc_all(j,ilev)
                           end do
                           j = j + 1
                        end do
                     end do
                  end do
                  call outfld('CRM_QRL', crm_qrl(1:ncol,:,:,:)/cpair, ncol, state%lchnk)
                  call outfld('CRM_QRLC', crm_qrlc(1:ncol,:,:,:)/cpair, ncol, state%lchnk)
               end if

               ! Set net fluxes used in other components
               call set_net_fluxes_lw(fluxes_allsky, flns, flnt)

               ! Export surface fluxes that are used by the land model
               call export_surface_fluxes(fluxes_allsky, cam_out, 'longwave')

               ! Free memory allocated for fluxes
               call free_fluxes(fluxes_allsky)
               call free_fluxes(fluxes_clrsky)
               call free_fluxes(fluxes_allsky_all)
               call free_fluxes(fluxes_clrsky_all)
               call free_optics_lw(cld_optics_lw_col)
               call free_optics_lw(cld_optics_lw_all)
               call free_optics_lw(aer_optics_lw_col)
               call free_optics_lw(aer_optics_lw_all)

            else

               ! Conserve energy (what does this mean exactly?)
               if (conserve_energy) then
                  qrl(1:ncol,1:pver) = qrl(1:ncol,1:pver) / state%pdel(1:ncol,1:pver)
               end if

            end if  ! radiation_do('lw')

         end if  ! active calls
      end do  ! loop over diagnostic calls

      ! Restore pbuf fields
      iclwp = iclwp_save
      iciwp = iciwp_save
      cld   = cld_save
      dei   = dei_save
      rei   = rei_save
      rel   = rel_save
#ifdef MODAL_AERO
      qaerwat = qaerwat_save
      dgnumwet = dgnumwet_save
#endif

      ! Compute net radiative heating tendency
      call t_startf('radheat_tend')
      call radheat_tend(state, pbuf, ptend, &
                        qrl, qrs, &
                        fsns, fsnt, flns, flnt, &
                        cam_in%asdir, net_flux)
      call t_stopf('radheat_tend')

      ! convert radiative heating rates to Q*dp to carry across timesteps
      ! for energy conservation
      if (conserve_energy) then
         qrs(1:ncol,1:pver) = qrs(1:ncol,1:pver) * state%pdel(1:ncol,1:pver)
         qrl(1:ncol,1:pver) = qrl(1:ncol,1:pver) * state%pdel(1:ncol,1:pver)
      end if

      ! Update net CRM heating tendency, IF doing radiation this timestep
      if (use_MMF) then
         call t_startf('rad_update_crm_heating')
         if (radiation_do('sw') .or. radiation_do('lw')) then
            call pbuf_get_field(pbuf, pbuf_get_index('CRM_QRAD'), crm_qrad)
            crm_qrad = 0
            do iz = 1,crm_nz
               do iy = 1,crm_ny_rad
                  do ix = 1,crm_nx_rad
                     do ic = 1,ncol
                        crm_qrad(ic,ix,iy,iz) = (crm_qrs(ic,ix,iy,iz) + crm_qrl(ic,ix,iy,iz)) / cpair
                     end do
                  end do
               end do
            end do
            call outfld('CRM_QRAD', crm_qrad(1:ncol,:,:,:), ncol, state%lchnk)
         end if
         call t_stopf('rad_update_crm_heating')
      end if  ! use_MMF

   end subroutine radiation_tend

   !----------------------------------------------------------------------------

   subroutine calculate_fluxes_sw(gas_names, gas_vmr, &
                                  pmid, tmid, pint, &
                                  coszrs, alb_dir, alb_dif, &
                                  cld_optics, aer_optics, &
                                  fluxes_allsky, fluxes_clrsky, &
                                  tsi_scaling)

      use perf_mod, only: t_startf, t_stopf
      use mo_fluxes_byband, only: ty_fluxes_byband
      use mo_optical_props, only: ty_optical_props_2str
      use mo_gas_concentrations, only: ty_gas_concs
      use mo_rrtmgp_util_string, only: lower_case
      use mo_rrtmgp_clr_all_sky, only: rte_sw
      use cam_optics, only: free_optics_sw

      character(len=*), intent(in) :: gas_names(:)
      real(r8), intent(in), dimension(:,:,:) :: gas_vmr
      real(r8), intent(in), dimension(:,:) :: pmid, tmid, pint
      real(r8), intent(in), dimension(:) :: coszrs
      real(r8), intent(in), dimension(:,:) :: alb_dir, alb_dif
      type(ty_optical_props_2str), intent(in) :: cld_optics, aer_optics
      type(ty_fluxes_byband), intent(inout) :: fluxes_allsky, fluxes_clrsky
      real(r8), intent(in) :: tsi_scaling

      ! For day-only arrays/objects
      type(ty_optical_props_2str) :: cld_optics_day, aer_optics_day
      type(ty_fluxes_byband) :: fluxes_allsky_day, fluxes_clrsky_day
      real(r8), dimension(size(pmid,1), size(pmid,2)) :: pmid_day, tmid_day
      real(r8), dimension(size(pint,1), size(pint,2)) :: pint_day
      real(r8), dimension(size(coszrs)) :: coszrs_day
      real(r8), dimension(size(alb_dir,1),size(alb_dir,2)) :: alb_dir_day, alb_dif_day
      real(r8), dimension(size(gas_vmr,1),size(gas_vmr,2),size(gas_vmr,3)) :: gas_vmr_day

      type(ty_gas_concs) :: gas_concs
      integer :: ncol, nday, nlev, igas, iday, icol
      integer, dimension(size(coszrs)) :: day_indices, night_indices

      ! Character array to hold lowercase gas names
      character(len=32), allocatable :: gas_names_lower(:)


      ncol = size(pmid,1)
      nlev = size(pmid,2)

      ! Get day-night indices
      call set_daynight_indices(coszrs, day_indices, night_indices)
      nday = count(day_indices > 0)

      ! If we don't have any daytime indices, zero fluxes and return
      if (nday <= 0) then
         call reset_fluxes(fluxes_allsky)
         call reset_fluxes(fluxes_clrsky)
         return
      end if

      ! Allocate daytime-only optics
      call handle_error(cld_optics_day%alloc_2str(nday, nlev, k_dist_sw, name='sw day-time cloud optics'))
      call handle_error(aer_optics_day%alloc_2str(nday, nlev, k_dist_sw%get_band_lims_wavenumber(), name='sw day-time aerosol optics'))

      ! Allocate fluxes
      ! TODO: use pointers to stack arrays to save from allocating/deallocating
      ! each step
      call initialize_rrtmgp_fluxes(nday, nlev+1, nswbands, fluxes_allsky_day, do_direct=.true.)
      call initialize_rrtmgp_fluxes(nday, nlev+1, nswbands, fluxes_clrsky_day, do_direct=.true.)

      ! Check incoming optical properties
      call handle_error(cld_optics%validate())
      call handle_error(aer_optics%validate())

      ! Compress to day-time only
      do iday = 1,nday
         icol = day_indices(iday)

         ! Compress state
         pmid_day(iday,:) = pmid(icol,:)
         tmid_day(iday,:) = tmid(icol,:)
         pint_day(iday,:) = pint(icol,:)

         ! Compress coszrs and albedo
         coszrs_day(iday) = coszrs(icol)
         alb_dir_day(:,iday) = alb_dir(:,icol)
         alb_dif_day(:,iday) = alb_dif(:,icol)

         ! Compress gases
         gas_vmr_day(:,iday,:) = gas_vmr(:,icol,:)

         ! Compress optics
         cld_optics_day%tau(iday,:,:) = cld_optics%tau(icol,:,:)
         cld_optics_day%ssa(iday,:,:) = cld_optics%ssa(icol,:,:)
         cld_optics_day%g  (iday,:,:) = cld_optics%g  (icol,:,:)
         aer_optics_day%tau(iday,:,:) = aer_optics%tau(icol,:,:)
         aer_optics_day%ssa(iday,:,:) = aer_optics%ssa(icol,:,:)
         aer_optics_day%g  (iday,:,:) = aer_optics%g  (icol,:,:)
      end do

      ! Check incoming optical properties
      call handle_error(cld_optics_day%validate())
      call handle_error(aer_optics_day%validate())

      ! Initialize gas concentrations with lower case names
      allocate(gas_names_lower(size(gas_names)))
      do igas = 1,size(gas_names)
         gas_names_lower(igas) = trim(lower_case(gas_names(igas)))
      end do
      call handle_error(gas_concs%init(gas_names_lower))

      ! Populate gas concentrations
      do igas = 1,size(gas_names)
         call handle_error(gas_concs%set_vmr( &
            gas_names_lower(igas), gas_vmr_day(igas,1:nday,:) &
         ))
      end do
      deallocate(gas_names_lower)

      call handle_error(cld_optics_day%validate())
      call handle_error(aer_optics_day%validate())
        
      ! Compute fluxes
      call t_startf('rad_rte_sw')
      call handle_error(rte_sw(k_dist_sw, gas_concs, &
                               pmid_day(1:nday,1:nlev), &
                               tmid_day(1:nday,1:nlev), &
                               pint_day(1:nday,1:nlev+1), &
                               coszrs_day(1:nday), &
                               alb_dir_day(1:nswbands,1:nday), &
                               alb_dif_day(1:nswbands,1:nday), &
                               cld_optics_day, &
                               fluxes_allsky_day, fluxes_clrsky_day, &
                               aer_props=aer_optics_day, &
                               tsi_scaling=tsi_scaling))
      call t_stopf('rad_rte_sw')

      ! Expand daytime-only fluxes
      call reset_fluxes(fluxes_allsky)
      call reset_fluxes(fluxes_clrsky)
      do iday = 1,nday
         icol = day_indices(iday)
         fluxes_allsky%flux_up(icol,:) = fluxes_allsky_day%flux_up(iday,:)
         fluxes_allsky%flux_dn(icol,:) = fluxes_allsky_day%flux_dn(iday,:)
         fluxes_allsky%flux_net(icol,:) = fluxes_allsky_day%flux_net(iday,:)
         fluxes_allsky%bnd_flux_up(icol,:,:) = fluxes_allsky_day%bnd_flux_up(iday,:,:)
         fluxes_allsky%bnd_flux_dn(icol,:,:) = fluxes_allsky_day%bnd_flux_dn(iday,:,:)
         fluxes_allsky%bnd_flux_net(icol,:,:) = fluxes_allsky_day%bnd_flux_net(iday,:,:)
         fluxes_allsky%bnd_flux_dn_dir(icol,:,:) = fluxes_allsky_day%bnd_flux_dn_dir(iday,:,:)
         fluxes_clrsky%flux_up(icol,:) = fluxes_clrsky_day%flux_up(iday,:)
         fluxes_clrsky%flux_dn(icol,:) = fluxes_clrsky_day%flux_dn(iday,:)
         fluxes_clrsky%flux_net(icol,:) = fluxes_clrsky_day%flux_net(iday,:)
         fluxes_clrsky%bnd_flux_up(icol,:,:) = fluxes_clrsky_day%bnd_flux_up(iday,:,:)
         fluxes_clrsky%bnd_flux_dn(icol,:,:) = fluxes_clrsky_day%bnd_flux_dn(iday,:,:)
         fluxes_clrsky%bnd_flux_net(icol,:,:) = fluxes_clrsky_day%bnd_flux_net(iday,:,:)
         fluxes_clrsky%bnd_flux_dn_dir(icol,:,:) = fluxes_clrsky_day%bnd_flux_dn_dir(iday,:,:)
      end do

      ! Free memory
      call free_optics_sw(cld_optics_day)
      call free_optics_sw(aer_optics_day)
      call free_fluxes(fluxes_allsky_day)
      call free_fluxes(fluxes_clrsky_day)

   end subroutine calculate_fluxes_sw

   !----------------------------------------------------------------------------

   subroutine calculate_fluxes_lw(gas_names, gas_vmr, emis_sfc, &
                                pmid, tmid, pint, tint, &
                                cld_optics, aer_optics, &
                                fluxes_allsky, fluxes_clrsky)

      use perf_mod, only: t_startf, t_stopf
      use mo_rrtmgp_clr_all_sky, only: rte_lw
      use mo_fluxes_byband, only: ty_fluxes_byband
      use mo_optical_props, only: ty_optical_props_1scl
      use mo_gas_concentrations, only: ty_gas_concs
      use mo_rrtmgp_util_string, only: lower_case

      character(len=*), intent(in) :: gas_names(:)
      real(r8), intent(in) :: gas_vmr(:,:,:)
      real(r8), intent(in) :: emis_sfc(:,:)
      real(r8), intent(in) :: pmid(:,:), tmid(:,:), pint(:,:), tint(:,:)
      type(ty_optical_props_1scl), intent(in) :: cld_optics, aer_optics
      type(ty_fluxes_byband), intent(inout) :: fluxes_allsky, fluxes_clrsky

      type(ty_gas_concs) :: gas_concs
      integer :: ncol, nlev, igas

      ! Character array to hold lowercase gas names
      character(len=32), allocatable :: gas_names_lower(:)

      ncol = size(pmid,1)
      nlev = size(pmid,2)

      ! Initialize gas concentrations with lower case names
      allocate(gas_names_lower(size(gas_names)))
      do igas = 1,size(gas_names)
         gas_names_lower(igas) = trim(lower_case(gas_names(igas)))
      end do
      call handle_error(gas_concs%init(gas_names_lower))

      ! Populate gas concentrations
      do igas = 1,size(gas_names)
         call handle_error(gas_concs%set_vmr( &
            gas_names_lower(igas), gas_vmr(igas,:,:) &
         ))
      end do
      deallocate(gas_names_lower)


      ! Do longwave radiative transfer calculations
      call t_startf('rad_rte_lw')
      call handle_error(rte_lw(k_dist_lw, gas_concs, &
                               pmid(1:ncol,1:nlev), tmid(1:ncol,1:nlev), &
                               pint(1:ncol,1:nlev+1), tint(1:ncol,nlev+1), &
                               emis_sfc(1:nlwbands,1:ncol), &
                               cld_optics, &
                               fluxes_allsky, fluxes_clrsky, &
                               aer_props=aer_optics, &
                               t_lev=tint(1:ncol,1:nlev+1), &
                               n_gauss_angles=1))  ! Set to 3 for consistency with RRTMG
      call t_stopf('rad_rte_lw')

   end subroutine calculate_fluxes_lw

   !----------------------------------------------------------------------------

   ! Utility routine to compute domain averages of CRM data that has been
   ! "packed" into a single dimension to hold both global GCM columns and CRM
   ! "subcolumns" within each GCM column. Input will be 2D arrays dimensioned
   ! (ncol_tot,nlev) where ncol_tot is the total number of CRM columns ncol *
   ! crm_nx_rad * crm_ny_rad and nlev is number of vertical levels. Output will
   ! be 2D arrays dimensioned (ncol,nlev), where the averaging has been done
   ! over the CRM columns.
   subroutine average_packed_array(array_packed, array_avg)
      real(r8), intent(in) :: array_packed(:,:)
      real(r8), intent(out) :: array_avg(:,:)
      integer :: ncol, ncol_tot
      real(r8) :: area_factor
      integer :: ic, ix, iy, iz, j

      if (crm_nx_rad * crm_ny_rad > 1) then
         area_factor = 1._r8 / (crm_nx_rad * crm_ny_rad)
      else
         area_factor = 1
      end if
      array_avg = 0
      ncol = size(array_avg, 1)
      ncol_tot = ncol * crm_nx_rad * crm_ny_rad
      call assert(size(array_packed,1) == ncol_tot, 'size(array_packed,1) /= ncol_tot')
      call assert(size(array_packed,2) == size(array_avg,2), 'size(array_packed,2) /= size(array_avg,2)')
      do iz = 1,size(array_packed,2)
         j = 1
         do iy = 1,crm_ny_rad
            do ix = 1,crm_nx_rad
               do ic = 1,ncol
                  array_avg(ic,iz) = array_avg(ic,iz) + array_packed(j,iz) * area_factor
                  j = j + 1
               end do
            end do
         end do 
      end do
   end subroutine average_packed_array

   !----------------------------------------------------------------------------

   subroutine export_surface_fluxes(fluxes, cam_out, band)
      use mo_fluxes_byband, only: ty_fluxes_byband
      use camsrfexch, only: cam_out_t

      type(ty_fluxes_byband), intent(in) :: fluxes
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

      use mo_fluxes_byband, only: ty_fluxes_byband
      type(ty_fluxes_byband), intent(in) :: fluxes
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

      use mo_fluxes_byband, only: ty_fluxes_byband
      type(ty_fluxes_byband), intent(in) :: fluxes
      real(r8), intent(inout) :: flns(:)
      real(r8), intent(inout) :: flnt(:)
      integer :: ncol

      ! Copy data
      ncol = size(fluxes%flux_up, 1)
      flns(1:ncol) = fluxes%flux_up(1:ncol,kbot+1) - fluxes%flux_dn(1:ncol,kbot+1)
      flnt(1:ncol) = fluxes%flux_up(1:ncol,ktop) - fluxes%flux_dn(1:ncol,ktop)

   end subroutine set_net_fluxes_lw

   !----------------------------------------------------------------------------

   subroutine reset_fluxes(fluxes)

      use mo_rte_kind, only: wp
      use mo_fluxes_byband, only: ty_fluxes_byband
      type(ty_fluxes_byband), intent(inout) :: fluxes

      ! Reset broadband fluxes
      fluxes%flux_up(:,:) = 0._wp
      fluxes%flux_dn(:,:) = 0._wp
      fluxes%flux_net(:,:) = 0._wp
      if (associated(fluxes%flux_dn_dir)) fluxes%flux_dn_dir(:,:) = 0._wp

      ! Reset band-by-band fluxes
      fluxes%bnd_flux_up(:,:,:) = 0._wp
      fluxes%bnd_flux_dn(:,:,:) = 0._wp
      fluxes%bnd_flux_net(:,:,:) = 0._wp
      if (associated(fluxes%bnd_flux_dn_dir)) fluxes%bnd_flux_dn_dir(:,:,:) = 0._wp

   end subroutine reset_fluxes


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
   subroutine set_albedo(cam_in, albedo_direct, albedo_diffuse)
      use camsrfexch, only: cam_in_t
      use radiation_utils, only: clip_values

      type(cam_in_t), intent(in) :: cam_in
      real(r8), intent(inout) :: albedo_direct(:,:)   ! surface albedo, direct radiation
      real(r8), intent(inout) :: albedo_diffuse(:,:)  ! surface albedo, diffuse radiation

      ! Local namespace
      real(r8) :: wavenumber_limits(2,nswbands)
      integer :: ncol, iband

      ! Check dimension sizes of output arrays.
      ! albedo_direct and albedo_diffuse should have sizes nswbands,ncol, but ncol
      ! can change so we just check that it is less than or equal to pcols (the
      ! maximum size ncol is ever allowed to be).
      call assert(size(albedo_direct, 1) == nswbands, &
                  'set_albedo: size(albedo_direct, 1) /= nswbands')
      call assert(size(albedo_direct, 2) <= pcols, &
                  'set_albedo: size(albedo_direct, 2) > pcols')
      call assert(all(shape(albedo_direct) == shape(albedo_diffuse)), &
                  'set_albedo: albedo_direct and albedo_diffuse have inconsistent shapes')
      
      ncol = size(albedo_direct, 2)

      ! Initialize albedo
      albedo_direct(:,:) = 0._r8
      albedo_diffuse(:,:) = 0._r8

      ! Albedos are input as broadband (visible, and near-IR), and we need to map
      ! these to appropriate bands. Bands are categorized broadly as "visible" or
      ! "infrared" based on wavenumber, so we get the wavenumber limits here
      wavenumber_limits = k_dist_sw%get_band_lims_wavenumber()

      ! Loop over bands, and determine for each band whether it is broadly in the
      ! visible or infrared part of the spectrum (visible or "not visible")
      do iband = 1,nswbands
         if (is_visible(wavenumber_limits(1,iband)) .and. &
             is_visible(wavenumber_limits(2,iband))) then

            ! Entire band is in the visible
            albedo_direct(iband,1:ncol) = cam_in%asdir(1:ncol)
            albedo_diffuse(iband,1:ncol) = cam_in%asdif(1:ncol)

         else if (.not.is_visible(wavenumber_limits(1,iband)) .and. &
                  .not.is_visible(wavenumber_limits(2,iband))) then

            ! Entire band is in the longwave (near-infrared)
            albedo_direct(iband,1:ncol) = cam_in%aldir(1:ncol)
            albedo_diffuse(iband,1:ncol) = cam_in%aldif(1:ncol)

         else

            ! Band straddles the visible to near-infrared transition, so we take
            ! the albedo to be the average of the visible and near-infrared
            ! broadband albedos
            albedo_direct(iband,1:ncol) = 0.5 * (cam_in%aldir(1:ncol) + cam_in%asdir(1:ncol))
            albedo_diffuse(iband,1:ncol) = 0.5 * (cam_in%aldif(1:ncol) + cam_in%asdif(1:ncol))

         end if
      end do

      ! Check values and clip if necessary (albedos should not be larger than 1)
      ! NOTE: this does actually issue warnings for albedos larger than 1, but this
      ! was never checked for RRTMG, so albedos will probably be slight different
      ! than the implementation in RRTMG!
      call clip_values(albedo_direct, 0._r8, 1._r8, varname='albedo_direct')
      call clip_values(albedo_diffuse, 0._r8, 1._r8, varname='albedo_diffuse')

   end subroutine set_albedo

   !-------------------------------------------------------------------------------

   ! Function to check if a wavenumber is in the visible or IR
   logical function is_visible(wavenumber)

      use mo_rte_kind, only: wp
      
      ! Input wavenumber; this needs to be input in inverse cm (cm^-1)
      real(wp), intent(in) :: wavenumber

      ! Threshold between visible and infrared is 0.7 micron, or 14286 cm^-1
      real(wp), parameter :: visible_wavenumber_threshold = 14286._wp  ! cm^-1

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
      use mo_fluxes_byband, only: ty_fluxes_byband
      
      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(ty_fluxes_byband), intent(in) :: flux_all
      type(ty_fluxes_byband), intent(in) :: flux_clr
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
      use mo_fluxes_byband, only: ty_fluxes_byband
      
      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(ty_fluxes_byband), intent(in) :: flux_all
      type(ty_fluxes_byband), intent(in) :: flux_clr

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


   ! For some reason the RRTMGP flux objects do not include initialization
   ! routines, but rather expect the user to associate the individual fluxes (which
   ! are pointers) with appropriate targets. Instead, this routine treats those
   ! pointers as allocatable members and allocates space for them. TODO: is this
   ! appropriate use?
   subroutine initialize_rrtmgp_fluxes(ncol, nlevels, nbands, fluxes, do_direct)

      use mo_fluxes_byband, only: ty_fluxes_byband
      integer, intent(in) :: ncol, nlevels, nbands
      type(ty_fluxes_byband), intent(inout) :: fluxes
      logical, intent(in), optional :: do_direct

      logical :: do_direct_local

      if (present(do_direct)) then
         do_direct_local = .true.
      else
         do_direct_local = .false.
      end if

      ! Allocate flux arrays
      ! NOTE: fluxes defined at interfaces, so need to either pass nlevels as
      ! number of model levels plus one, or allocate as nlevels+1 if nlevels
      ! represents number of model levels rather than number of interface levels.

      ! Broadband fluxes
      allocate(fluxes%flux_up(ncol, nlevels))
      allocate(fluxes%flux_dn(ncol, nlevels))
      allocate(fluxes%flux_net(ncol, nlevels))
      if (do_direct_local) allocate(fluxes%flux_dn_dir(ncol, nlevels))

      ! Fluxes by band
      allocate(fluxes%bnd_flux_up(ncol, nlevels, nbands))
      allocate(fluxes%bnd_flux_dn(ncol, nlevels, nbands))
      allocate(fluxes%bnd_flux_net(ncol, nlevels, nbands))
      if (do_direct_local) allocate(fluxes%bnd_flux_dn_dir(ncol, nlevels, nbands))

      ! Initialize
      call reset_fluxes(fluxes)

   end subroutine initialize_rrtmgp_fluxes

   subroutine free_fluxes(fluxes)
      use mo_fluxes_byband, only: ty_fluxes_byband
      type(ty_fluxes_byband), intent(inout) :: fluxes
      if (associated(fluxes%flux_up)) deallocate(fluxes%flux_up)
      if (associated(fluxes%flux_dn)) deallocate(fluxes%flux_dn)
      if (associated(fluxes%flux_net)) deallocate(fluxes%flux_net)
      if (associated(fluxes%flux_dn_dir)) deallocate(fluxes%flux_dn_dir)
      if (associated(fluxes%bnd_flux_up)) deallocate(fluxes%bnd_flux_up)
      if (associated(fluxes%bnd_flux_dn)) deallocate(fluxes%bnd_flux_dn)
      if (associated(fluxes%bnd_flux_net)) deallocate(fluxes%bnd_flux_net)
      if (associated(fluxes%bnd_flux_dn_dir)) deallocate(fluxes%bnd_flux_dn_dir)
   end subroutine free_fluxes


   subroutine get_gas_vmr(icall, state, pbuf, gas_name, vmr)

      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc
      use rad_constituents, only: rad_cnst_get_gas
      use mo_rrtmgp_util_string, only: lower_case, string_loc_in_array

      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      character(len=*), intent(in) :: gas_name
      real(r8), intent(out) :: vmr(:,:)

      ! Mass mixing ratio
      real(r8), pointer :: mass_mix_ratio(:,:)

      ! Gases and molecular weights. Note that we do NOT have CFCs yet (I think
      ! this is coming soon in RRTMGP). RRTMGP also allows for absorption due to
      ! CO and N2, which RRTMG did not have.
      character(len=3), dimension(8) :: gas_species = (/ &
         'H2O', 'CO2', 'O3 ', 'N2O', &
         'CO ', 'CH4', 'O2 ', 'N2 ' &
      /)
      real(r8), dimension(8) :: mol_weight_gas = (/ &
         18.01528, 44.0095, 47.9982, 44.0128, &
         28.0101, 16.04246, 31.998, 28.0134 &
      /)  ! g/mol

      ! Molar weight of air
      real(r8), parameter :: mol_weight_air = 28.97  ! g/mol
                                       
      ! Defaults for gases that are not available (TODO: is this still accurate?)
      real(r8), parameter :: co_vmr = 1.0e-7_r8
      real(r8), parameter :: n2_vmr = 0.7906_r8

      ! Loop indices
      integer :: igas

      ! Number of columns
      integer :: ncol

      ! Name of routine
      character(len=32) :: subname = 'get_gas_vmr'

      ! Number of columns in chunk
      ncol = state%ncol

      ! Get index into gas names we define above that we know the molecular 
      ! weights for; if this gas is not in list of gases we know about, skip
      igas = string_loc_in_array(gas_name, gas_species)
      if (igas <= 0) then
         call endrun('Gas name ' // trim(gas_name) // ' not recognized.')
      end if
         
      ! initialize
      vmr(:,:) = 0._r8

      select case(trim(gas_species(igas)))

         case('CO')

            ! CO not available, use default
            vmr(1:ncol,1:pver) = co_vmr

         case('N2')

            ! N2 not available, use default
            vmr(1:ncol,1:pver) = n2_vmr

         case('H2O')

            ! Water vapor is represented as specific humidity in CAM, so we
            ! need to handle water a little differently; first, read water vapor
            ! specific humidity into mass mixing ratio array
            call rad_cnst_get_gas(icall, trim(gas_species(igas)), state, pbuf, &
                                  mass_mix_ratio)
            
            ! Convert to volume mixing ratio by multiplying by the ratio of
            ! molecular weight of dry air to molecular weight of gas. Note that
            ! first specific humidity (held in the mass_mix_ratio array read
            ! from rad_constituents) is converted to an actual mass mixing
            ! ratio.
            vmr(1:ncol,1:pver) = mass_mix_ratio(1:ncol,1:pver) / ( &
               1._r8 - mass_mix_ratio(1:ncol,1:pver) &
            )  * mol_weight_air / mol_weight_gas(igas)

         case DEFAULT

            ! Get mass mixing ratio from the rad_constituents interface
            call rad_cnst_get_gas(icall, trim(gas_species(igas)), state, pbuf,  mass_mix_ratio)

            ! Convert to volume mixing ratio by multiplying by the ratio of
            ! molecular weight of dry air to molecular weight of gas
            vmr(1:ncol,1:pver) = mass_mix_ratio(1:ncol,1:pver) &
                               * mol_weight_air / mol_weight_gas(igas)

      end select

   end subroutine get_gas_vmr
   

end module radiation
