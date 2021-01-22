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
   use cam_history_support, only: add_hist_coord
   use scamMod,          only: scm_crm_mode, single_column, swrad_off
   use rad_constituents, only: N_DIAG
   use radconstants,     only: nswbands, nlwbands, &
      get_band_index_sw, get_band_index_lw, test_get_band_index, &
      get_sw_spectral_midpoints, get_lw_spectral_midpoints, &
      rrtmg_to_rrtmgp_swbands
   use physconst, only: cpair, cappa

   ! RRTMGP gas optics object to store coefficient information. This is imported
   ! here so that we can make the k_dist objects module data and only load them
   ! once.
   use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
   use mo_rte_kind, only: wp

   ! Use my assertion routines to perform sanity checks
   use assertions, only: assert, assert_valid, assert_range

   use radiation_state, only: ktop, kbot, nlev_rad
   use radiation_utils, only: compress_day_columns, expand_day_columns, &
                              handle_error

   use pio,            only: file_desc_t, pio_nowrite
   use cam_pio_utils,    only: cam_pio_openfile,cam_pio_closefile
   use cam_grid_support, only: cam_grid_check, cam_grid_id,cam_grid_get_dim_names
   use ncdio_atm,       only: infld
   use time_manager,   only: get_curr_date
   use cam_logfile,     only: iulog

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

   ! Logical to indicate if aerosol optical properties are read in directly from
   ! input file i.e., SPA:
   logical :: do_SPA_optics = .false.

   character(len=16), parameter :: unset_str = 'UNSET'

   character(len=128) :: SPA_lookup_dir     = unset_str ! location of SPA input files

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

   ! k-distribution coefficients. These will be populated by reading from the
   ! RRTMGP coefficients files, specified by coefficients_file_sw and
   ! coefficients_file_lw in the radiation namelist. They exist as module data
   ! because we only want to load those files once.
   type(ty_gas_optics_rrtmgp) :: k_dist_sw, k_dist_lw

   ! k-distribution coefficients files to read from. These are set via namelist
   ! variables.
   character(len=cl) :: rrtmgp_coefficients_file_sw, rrtmgp_coefficients_file_lw

   ! Number of g-points in k-distribution (set based on absorption coefficient inputdata)
   integer :: nswgpts, nlwgpts

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

   ! Indices to pbuf fields
   integer :: cldfsnow_idx = 0

   !needed for SPA
   integer :: aer_tau_bnd_lw_mon_1_idx,aer_tau_bnd_lw_mon_2_idx, aer_tau_bnd_sw_mon_1_idx, aer_tau_bnd_sw_mon_2_idx,&
              aer_ssa_bnd_sw_mon_1_idx,aer_ssa_bnd_sw_mon_2_idx, aer_asm_bnd_sw_mon_1_idx, aer_asm_bnd_sw_mon_2_idx 
   integer :: current_month

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
                              rrtmgp_enable_temperature_warnings,&
                              do_SPA_optics, SPA_lookup_dir

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
      call mpibcast(do_SPA_optics, 1, mpi_logical, mstrid, mpicom, ierr)
      call mpibcast(SPA_lookup_dir, cl, mpi_character, mstrid, mpicom, ierr)
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

    !for SPA
      call pbuf_add_field('ATBLM_1', 'global', dtype_r8, (/pcols,pverp,nlwbands/), aer_tau_bnd_lw_mon_1_idx)
      call pbuf_add_field('ATBLM_2', 'global', dtype_r8,(/pcols,pverp,nlwbands/), aer_tau_bnd_lw_mon_2_idx)
      call pbuf_add_field('ATBSM_1', 'global', dtype_r8, (/pcols,pverp,nlwbands/), aer_tau_bnd_sw_mon_1_idx)
      call pbuf_add_field('ATBSM_2', 'global', dtype_r8,(/pcols,pverp,nlwbands/), aer_tau_bnd_sw_mon_2_idx)
      call pbuf_add_field('ASBSM_1', 'global', dtype_r8, (/pcols,pverp,nlwbands/), aer_ssa_bnd_sw_mon_1_idx)
      call pbuf_add_field('ASBSM_2', 'global', dtype_r8,(/pcols,pverp,nlwbands/), aer_ssa_bnd_sw_mon_2_idx)
      call pbuf_add_field('AABSM_1', 'global', dtype_r8, (/pcols,pverp,nlwbands/), aer_asm_bnd_sw_mon_1_idx)
      call pbuf_add_field('AABSM_2', 'global', dtype_r8,(/pcols,pverp,nlwbands/), aer_asm_bnd_sw_mon_2_idx)

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

   subroutine get_aerosol_optical_property_from_file(month_int,aerosol_optical_property_str,band_identifier,nbands,aerosol_optical_property)
      use pio,            only: file_desc_t, pio_nowrite
      use cam_pio_utils,    only: cam_pio_openfile,cam_pio_closefile
      use cam_grid_support, only: cam_grid_check,cam_grid_id,cam_grid_get_dim_names
      use ncdio_atm,       only: infld
      use time_manager,   only: get_curr_date
      use cam_logfile,     only: iulog

      integer, intent(in) :: month_int, nbands !nbands is either nswbands or nlwbands
      character*(*), intent(in) :: aerosol_optical_property_str, band_identifier
      real(r8), intent(inout) :: aerosol_optical_property(pcols,pver,nbands,begchunk:endchunk)

     !internal variables
      character(len=100) :: base_file_name, optics_lookup_dir
      character(len=500) :: filename
      character(len=20) :: dim1name, dim2name
      character(len=20) :: mon_str
      type(file_desc_t) :: nccn_ncid
      integer :: year, month, day, tod, next_month, grid_id
      logical :: found = .false.

      write(mon_str,*) month_int

      !assign base_file_name the name of the CCN file being used 
      base_file_name = "spa_optics_file_"

      mon_str = adjustl(mon_str)

     !retrieve the name of the relevant file by combining base file name !with
     !!month and full file path:
      filename = trim(SPA_lookup_dir)//'/'//trim(base_file_name)//trim(mon_str)//'.nc'

      grid_id = cam_grid_id('physgrid')

      call cam_grid_get_dim_names(grid_id, dim1name, dim2name)

      call cam_pio_openfile(nccn_ncid,filename,PIO_NOWRITE)

      call infld(aerosol_optical_property_str,nccn_ncid,dim1name,'lev',band_identifier,1,pcols,1,pver,1,nbands,begchunk,endchunk,&
           aerosol_optical_property, found, gridname='physgrid')

      call cam_pio_closefile(nccn_ncid)

    end subroutine get_aerosol_optical_property_from_file

   subroutine radiation_init(state,pbuf)
   !-------------------------------------------------------------------------------
   ! Purpose: Initialize the radiation parameterization and add fields to the
   ! history buffer
   !
   ! History: Copied from RRTMG implemenation.
   !-------------------------------------------------------------------------------
      use physics_buffer,     only: pbuf_get_index,pbuf_set_field, physics_buffer_desc
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
      use physics_types, only: physics_state
      use cam_logfile,     only: iulog
      use spmd_utils,       only: iam

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
   
      type(physics_buffer_desc), pointer :: pbuf(:,:)

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
      logical :: use_MMF  ! SPCAM flag

      character(len=128) :: error_message

      ! ty_gas_concs object that would normally hold volume mixing ratios for
      ! radiatively-important gases. Here, this is just used to provide the names
      ! of gases that are available in the model (needed by the kdist
      ! initialization routines that are called within the load_coefficients
      ! methods).
      type(ty_gas_concs) :: available_gases

      character(len=32) :: subname = 'radiation_init'

      !needed for SPA
      integer :: year, month, day, tod, next_month
      real(r8), pointer :: aerosol_optical_property(:,:,:,:)
      real(r8), pointer :: aerosol_optical_property_lw(:,:,:,:)
      real(r8), pointer :: aerosol_optical_property_sw(:,:,:,:)

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
      ! the other tasks!
      ! TODO: This needs to be fixed to ONLY read in the data if masterproc, and then broadcast
      call set_available_gases(active_gases, available_gases)
      call rrtmgp_load_coefficients(k_dist_sw, rrtmgp_coefficients_file_sw, available_gases)
      call rrtmgp_load_coefficients(k_dist_lw, rrtmgp_coefficients_file_lw, available_gases)

      ! Make sure number of bands in absorption coefficient files matches what we expect
      call assert(nswbands == k_dist_sw%get_nband(), 'nswbands does not match absorption coefficient data')
      call assert(nlwbands == k_dist_lw%get_nband(), 'nlwbands does not match absorption coefficient data')

      ! Number of gpoints depend on inputdata, so initialize here
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
      call get_sw_spectral_midpoints(sw_band_midpoints, 'cm-1')
      call get_lw_spectral_midpoints(lw_band_midpoints, 'cm-1')
      call add_hist_coord('swband', nswbands, 'Shortwave wavenumber', 'cm-1', sw_band_midpoints)
      call add_hist_coord('lwband', nlwbands, 'Longwave wavenumber', 'cm-1', lw_band_midpoints)

      ! Shortwave radiation
      call addfld('TOT_CLD_VISTAU', (/ 'lev' /), 'A',   '1', &
                  'Total gridbox cloud visible (550 nm) optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('TOT_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', &
                  'Total in-cloud visible (550 nm) optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('LIQ_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', &
                  'Liquid in-cloud visible (550 nm) optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('ICE_ICLD_VISTAU', (/ 'lev' /), 'A',  '1', &
                  'Ice in-cloud visible (550 nm) optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)

      call add_default('TOT_CLD_VISTAU',  1, ' ')
      call add_default('TOT_ICLD_VISTAU', 1, ' ')

      ! Band-by-band cloud optical properties
      call addfld('CLOUD_TAU_SW', (/'lev   ','swband'/), 'I', '1', &
                  'Cloud shortwave extinction optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('CLOUD_SSA_SW', (/'lev   ','swband'/), 'I', '1', &
                  'Cloud shortwave single scattering albedo', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('CLOUD_G_SW', (/'lev   ','swband'/), 'I', '1', &
                  'Cloud shortwave assymmetry parameter', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('CLOUD_TAU_LW', (/'lev   ','lwband'/), 'I', '1', &
                  'Cloud longwave absorption optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)

      ! Band-by-band shortwave albedos
      call addfld('SW_ALBEDO_DIR', (/'swband'/), 'I', '1', &
                  'Shortwave direct-beam albedo', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
      call addfld('SW_ALBEDO_DIF', (/'swband'/), 'I', '1', &
                  'Shortwave diffuse-beam albedo', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)

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

      call addfld('EMIS', (/ 'lev' /), 'A', '1', 'Cloud longwave (10.5 micron) emissivity')

      ! Add default fields for single column mode
      if (single_column.and.scm_crm_mode) then
         call add_default ('FUL     ', 1, ' ')
         call add_default ('FULC    ', 1, ' ')
         call add_default ('FDL     ', 1, ' ')
         call add_default ('FDLC    ', 1, ' ')
      endif

      ! HIRS/MSU diagnostic brightness temperatures
      if (dohirs) then
         call addfld (hirsname(1),horiz_only,'A','K','HIRS CH2 infra-red brightness temperature', &
                      sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld (hirsname(2),horiz_only,'A','K','HIRS CH4 infra-red brightness temperature', &
                      sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld (hirsname(3),horiz_only,'A','K','HIRS CH6 infra-red brightness temperature', &
                      sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld (hirsname(4),horiz_only,'A','K','HIRS CH8 infra-red brightness temperature', &
                      sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld (hirsname(5),horiz_only,'A','K','HIRS CH10 infra-red brightness temperature', &
                      sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld (hirsname(6),horiz_only,'A','K','HIRS CH11 infra-red brightness temperature', &
                      sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld (hirsname(7),horiz_only,'A','K','HIRS CH12 infra-red brightness temperature', &
                      sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld (msuname(1),horiz_only,'A','K','MSU CH1 microwave brightness temperature', &
                      sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld (msuname(2),horiz_only,'A','K','MSU CH2 microwave brightness temperature', &
                      sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld (msuname(3),horiz_only,'A','K','MSU CH3 microwave brightness temperature', &
                      sampling_seq='rad_lwsw', flag_xyfill=.true.)
         call addfld (msuname(4),horiz_only,'A','K','MSU CH4 microwave brightness temperature', &
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
      call addfld ('HR',(/ 'lev' /), 'A','K/s','Heating rate needed for d(theta)/dt computation', &
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

      if (do_SPA_optics) then !initialize SPA
         !find current_month
         call get_curr_date(year,month,day,tod)
         current_month = month
         if (month==12) then
            next_month = 1
         else
            next_month = month + 1
         end if

         allocate(aerosol_optical_property_lw(pcols,pver,nlwbands,begchunk:endchunk))
         
         call get_aerosol_optical_property_from_file(current_month,'AER_TAU_LW','lwband',nlwbands,aerosol_optical_property_lw)
         call pbuf_set_field(pbuf, aer_tau_bnd_lw_mon_1_idx,aerosol_optical_property_lw)
         call get_aerosol_optical_property_from_file(next_month,'AER_TAU_LW','lwband',nlwbands,aerosol_optical_property_lw)
         call pbuf_set_field(pbuf,aer_tau_bnd_lw_mon_2_idx,aerosol_optical_property_lw)

         deallocate(aerosol_optical_property_lw) 

         allocate(aerosol_optical_property_sw(pcols,pver,nswbands,begchunk:endchunk))
         
         call get_aerosol_optical_property_from_file(current_month,'AER_TAU_SW','swband',nswbands,aerosol_optical_property_sw)
         call pbuf_set_field(pbuf, aer_tau_bnd_sw_mon_1_idx,aerosol_optical_property_sw)
         call get_aerosol_optical_property_from_file(next_month,'AER_TAU_SW','swband',nswbands,aerosol_optical_property_sw)
         call pbuf_set_field(pbuf,aer_tau_bnd_sw_mon_2_idx,aerosol_optical_property_sw)

         call get_aerosol_optical_property_from_file(current_month,'AER_SSA_SW','swband',nswbands,aerosol_optical_property_sw)
         call pbuf_set_field(pbuf, aer_ssa_bnd_sw_mon_1_idx,aerosol_optical_property_sw)
         call get_aerosol_optical_property_from_file(next_month,'AER_SSA_SW','swband',nswbands,aerosol_optical_property_sw)
         call pbuf_set_field(pbuf,aer_ssa_bnd_sw_mon_2_idx,aerosol_optical_property_sw)

         call get_aerosol_optical_property_from_file(current_month,'AER_G_SW','swband',nswbands,aerosol_optical_property_sw)
         call pbuf_set_field(pbuf,aer_asm_bnd_sw_mon_1_idx,aerosol_optical_property_sw)
         call get_aerosol_optical_property_from_file(next_month,'AER_G_SW','swband',nswbands,aerosol_optical_property_sw)
         call pbuf_set_field(pbuf,aer_asm_bnd_sw_mon_2_idx,aerosol_optical_property_sw)

         deallocate(aerosol_optical_property_sw)

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


   !===============================================================================

   !----------------------------------------------------------------------------
   ! Main driver for radiation computation. Calls the shortwave and longwave
   ! drivers, and calculates the radiative heating from the resulting fluxes.
   ! Primary output from this routine is the heating tendency due to radiative
   ! transfer, as a ptend object.
   subroutine radiation_tend(state,   ptend,    pbuf,          cam_out, cam_in,  &
                             landfrac,landm,    icefrac,       snowh,            &
                             fsns,    fsnt,     flns,          flnt,             &
                             fsds,    net_flux, is_cmip6_volc, dt                )

      ! Performance module needed for timing functions
      use perf_mod, only: t_startf, t_stopf

      ! CAM derived types; needed for surface exchange fields, physics state, and
      ! tendency fields
      use camsrfexch, only: cam_out_t, cam_in_t
      use physics_types, only: physics_state, physics_ptend

      ! Utilities for interacting with the physics buffer
      use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
                                pbuf_get_index

      ! For calculating radiative heating tendencies
      use radheat, only: radheat_tend

      ! For getting radiative constituent gases
      use rad_constituents, only: N_DIAG, rad_cnst_get_call_list

      ! RRTMGP radiation drivers and derived types
      use mo_fluxes_byband, only: ty_fluxes_byband

      ! CAM history module provides subroutine to send output data to the history
      ! buffer to be aggregated and written to disk
      use cam_history, only: outfld

      use radiation_state, only: set_rad_state
      use radiation_utils, only: clip_values
      use cam_control_mod, only: aqua_planet

      use cam_optics, only: get_cloud_optics_sw, sample_cloud_optics_sw, &
                            get_cloud_optics_lw, sample_cloud_optics_lw, &
                            set_aerosol_optics_sw
      use aer_rad_props, only: aer_rad_props_lw

      ! For running CFMIP Observation Simulator Package (COSP)
      use cospsimulator_intr, only: docosp, cospsimulator_intr_run, cosp_nradsteps

      ! ---------------------------------------------------------------------------
      ! Arguments
      ! ---------------------------------------------------------------------------

      ! Physics state variables
      type(physics_state), intent(in), target :: state

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

      ! Pointers to fields on the physics buffer
      real(r8), pointer, dimension(:,:) :: &
         cld, cldfsnow, &
         iclwp, iciwp, icswp, dei, des, lambdac, mu, &
         rei, rel

      ! Clear-sky heating rates are not on the physics buffer, and we have no
      ! reason to put them there, so declare these are regular arrays here
      real(r8) :: qrsc(pcols,pver)
      real(r8) :: qrlc(pcols,pver)

      ! Temporary variable for heating rate output
      real(r8) :: hr(pcols,pver)

      real(r8), dimension(pcols,nlev_rad) :: tmid, pmid
      real(r8), dimension(pcols,nlev_rad+1) :: pint, tint
      real(r8), dimension(nswbands,pcols) :: albedo_dir, albedo_dif

      ! Cosine solar zenith angle for all columns in chunk
      real(r8) :: coszrs(pcols)

      ! Cloud, snow, and aerosol optical properties
      real(r8), dimension(pcols,pver,nswgpts ) :: cld_tau_gpt_sw, cld_ssa_gpt_sw, cld_asm_gpt_sw
      real(r8), dimension(pcols,pver,nswbands) :: cld_tau_bnd_sw, cld_ssa_bnd_sw, cld_asm_bnd_sw
      real(r8), dimension(pcols,pver,nswbands) :: aer_tau_bnd_sw, aer_ssa_bnd_sw, aer_asm_bnd_sw
      real(r8), dimension(pcols,pver,nlwbands) :: cld_tau_bnd_lw, aer_tau_bnd_lw
      real(r8), dimension(pcols,pver,nlwgpts ) :: cld_tau_gpt_lw
      ! NOTE: these are diagnostic only
      real(r8), dimension(pcols,pver,nswbands) :: liq_tau_bnd_sw, ice_tau_bnd_sw, snw_tau_bnd_sw
      real(r8), dimension(pcols,pver,nlwbands) :: liq_tau_bnd_lw, ice_tau_bnd_lw, snw_tau_bnd_lw

      ! Needed for COSP: cloud lw emissivity and gridbox mean snow sw and lw optical depth
      real(r8), dimension(pcols,pver) :: cld_emis_lw, gb_snow_tau_sw, gb_snow_tau_lw

      ! COSP simulator bands
      integer :: cosp_swband, cosp_lwband

      ! Gas volume mixing ratios
      real(r8), dimension(size(active_gases),pcols,pver) :: gas_vmr

      ! Needed for shortwave aerosol; TODO: remove this dependency
      integer :: nday, nnight     ! Number of daylight columns
      integer :: day_indices(pcols), night_indices(pcols)   ! Indicies of daylight coumns

      ! Flag to carry (QRS,QRL)*dp across time steps.
      ! TODO: what does this mean?
      logical :: conserve_energy = .true.

      ! Number of columns
      integer :: ncol

      ! For loops over diagnostic calls
      integer :: icall
      logical :: active_calls(0:N_DIAG)

      integer :: icol, ilay

      ! Everyone needs a name
      character(*), parameter :: subname = 'radiation_tend'

      ! Radiative fluxes
      type(ty_fluxes_byband) :: fluxes_allsky, fluxes_clrsky

      ! Zero-array for cloud properties if not diagnosed by microphysics
      real(r8), target, dimension(pcols,pver) :: zeros
      real(r8), dimension(pcols,pver) :: c_cldf  ! Combined cloud/snow fraction
 
      !needed for SPA
      integer :: year, month, day, tod, next_month
      real(r8) :: fraction_of_month
      real(r8), pointer :: aerosol_optical_property_lw(:,:,:,:)
      real(r8), pointer :: aerosol_optical_property_sw(:,:,:,:)
      real(r8), pointer :: aer_tau_bnd_lw_mon_1(:,:,:)
      real(r8), pointer :: aer_tau_bnd_lw_mon_2(:,:,:)
      real(r8), pointer :: aer_tau_bnd_sw_mon_1(:,:,:)
      real(r8), pointer :: aer_tau_bnd_sw_mon_2(:,:,:)
      real(r8), pointer :: aer_ssa_bnd_sw_mon_1(:,:,:)
      real(r8), pointer :: aer_ssa_bnd_sw_mon_2(:,:,:)
      real(r8), pointer :: aer_asm_bnd_sw_mon_1(:,:,:)
      real(r8), pointer :: aer_asm_bnd_sw_mon_2(:,:,:)

      real(r8), dimension(12):: days_per_month
      nullify(aerosol_optical_property_lw)
      nullify(aerosol_optical_property_sw)
      !fill days_per_month, SPA doesn't recognize leap year
      days_per_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)

      !----------------------------------------------------------------------

      ! Number of physics columns in this "chunk"
      ncol = state%ncol

      ! Set pointers to heating rates stored on physics buffer. These will be
      ! modified in this routine.
      call pbuf_get_field(pbuf, pbuf_get_index('QRS'), qrs)
      call pbuf_get_field(pbuf, pbuf_get_index('QRL'), qrl)

      ! Get fields from pbuf for optics
      call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cld)
      call pbuf_get_field(pbuf, pbuf_get_index('ICLWP'), iclwp)
      call pbuf_get_field(pbuf, pbuf_get_index('ICIWP'), iciwp)
      call pbuf_get_field(pbuf, pbuf_get_index('DEI'), dei)
      call pbuf_get_field(pbuf, pbuf_get_index('REL'), rel)
      call pbuf_get_field(pbuf, pbuf_get_index('REI'), rei)
      call pbuf_get_field(pbuf, pbuf_get_index('LAMBDAC'), lambdac)
      call pbuf_get_field(pbuf, pbuf_get_index('MU'), mu)
      ! May or may not have snow properties depending on microphysics
      if (do_snow_optics()) then
         call pbuf_get_field(pbuf, pbuf_get_index('CLDFSNOW'), cldfsnow)
         call pbuf_get_field(pbuf, pbuf_get_index('ICSWP'), icswp)
         call pbuf_get_field(pbuf, pbuf_get_index('DES'), des)
      else
         zeros = 0
         cldfsnow => zeros
         icswp => zeros
         des => zeros
      end if

      ! Compute combined cloud and snow fraction
      do icol = 1,size(cld,1)
         do ilay = 1,size(cld,2)
            c_cldf(icol,ilay) = max(cld(icol,ilay), cldfsnow(icol,ilay))
         end do
      end do

      ! Initialize clearsky-heating rates to make sure we do not get garbage
      ! for columns beyond ncol or nday
      qrsc(:,:) = 0
      qrlc(:,:) = 0

      if (radiation_do('sw') .or. radiation_do('lw')) then
         ! Make copies of state variables because we may need to modify in-place, and because we need
         ! to add a level above model top to account for heating
         call set_rad_state( &
            state, cam_in, &
            tmid(1:ncol,1:nlev_rad), tint(1:ncol,1:nlev_rad+1), &
            pmid(1:ncol,1:nlev_rad), pint(1:ncol,1:nlev_rad+1) &
         )

         ! Check temperatures to make sure they are within the bounds of the
         ! absorption coefficient look-up tables. If out of bounds, clip
         ! values to min/max specified
         call t_startf('rrtmgp_check_temperatures')
         call handle_error(clip_values( &
            tmid(1:ncol,1:nlev_rad), k_dist_lw%get_temp_min(), k_dist_lw%get_temp_max(), &
            trim(subname) // ' tmid' &
         ), fatal=.false., warn=rrtmgp_enable_temperature_warnings)
         call handle_error(clip_values( &
            tint(1:ncol,1:nlev_rad+1), k_dist_lw%get_temp_min(), k_dist_lw%get_temp_max(), &
            trim(subname) // ' tint' &
         ), fatal=.false., warn=rrtmgp_enable_temperature_warnings)
         call t_stopf('rrtmgp_check_temperatures')
      end if

      ! Do shortwave stuff...
      if (radiation_do('sw')) then

         ! Allocate shortwave fluxes
         call initialize_rrtmgp_fluxes(ncol, nlev_rad+1, nswbands, fluxes_allsky, do_direct=.true.)
         call initialize_rrtmgp_fluxes(ncol, nlev_rad+1, nswbands, fluxes_clrsky, do_direct=.true.)

         ! Get albedo. This uses CAM routines internally and just provides a
         ! wrapper to improve readability of the code here.
         call set_albedo(cam_in, albedo_dir(1:nswbands,1:ncol), albedo_dif(1:nswbands,1:ncol))

         ! Send albedos to history buffer (useful for debugging)
         call outfld('SW_ALBEDO_DIR', transpose(albedo_dir(1:nswbands,1:ncol)), ncol, state%lchnk)
         call outfld('SW_ALBEDO_DIF', transpose(albedo_dif(1:nswbands,1:ncol)), ncol, state%lchnk)

         ! Get cosine solar zenith angle for current time step.
         call set_cosine_solar_zenith_angle(state, dt_avg, coszrs(1:ncol))
         call outfld('COSZRS', coszrs(1:ncol), ncol, state%lchnk)
         ! If the swrad_off flag is set, meaning we should not do SW radiation, then
         ! we just set coszrs to zero everywhere. TODO: why not just set dosw false
         ! and skip the loop?
         if (swrad_off) coszrs(:) = 0._r8

         ! Do shortwave cloud optics calculations
         call t_startf('rad_cld_optics_sw')
         cld_tau_gpt_sw = 0._r8
         cld_ssa_gpt_sw = 0._r8
         cld_asm_gpt_sw = 0._r8
         call get_cloud_optics_sw( &
            ncol, pver, nswbands, do_snow_optics(), cld, cldfsnow, iclwp, iciwp, icswp, &
            lambdac, mu, dei, des, rel, rei, &
            cld_tau_bnd_sw, cld_ssa_bnd_sw, cld_asm_bnd_sw, &
            liq_tau_bnd_sw, ice_tau_bnd_sw, snw_tau_bnd_sw &
         )
         ! Output the band-by-band cloud optics BEFORE we reorder bands, because
         ! we hard-coded the indices for diagnostic bands in radconstants.F90 to
         ! correspond to the optical property look-up tables.
         call output_cloud_optics_sw( &
            state, coszrs, c_cldf, cld_tau_bnd_sw, cld_ssa_bnd_sw, cld_asm_bnd_sw, &
            liq_tau_bnd_sw, ice_tau_bnd_sw &
         )
         ! Now reorder bands to be consistent with RRTMGP
         ! We need to fix band ordering because the old input files assume RRTMG
         ! band ordering, but this has changed in RRTMGP.
         ! TODO: fix the input files themselves!
         do icol = 1,size(cld_tau_bnd_sw,1)
            do ilay = 1,size(cld_tau_bnd_sw,2)
               cld_tau_bnd_sw(icol,ilay,:) = reordered(cld_tau_bnd_sw(icol,ilay,:), rrtmg_to_rrtmgp_swbands)
               cld_ssa_bnd_sw(icol,ilay,:) = reordered(cld_ssa_bnd_sw(icol,ilay,:), rrtmg_to_rrtmgp_swbands)
               cld_asm_bnd_sw(icol,ilay,:) = reordered(cld_asm_bnd_sw(icol,ilay,:), rrtmg_to_rrtmgp_swbands)
            end do
         end do
         ! And now do the MCICA sampling to get cloud optical properties by
         ! gpoint/cloud state
         call sample_cloud_optics_sw( &
            ncol, pver, nswgpts, k_dist_sw%get_gpoint_bands(), &
            state%pmid, cld, cldfsnow, &
            cld_tau_bnd_sw, cld_ssa_bnd_sw, cld_asm_bnd_sw, &
            cld_tau_gpt_sw, cld_ssa_gpt_sw, cld_asm_gpt_sw &
         )
         call t_stopf('rad_cld_optics_sw')

         ! Aerosol needs night indices
         ! TODO: remove this dependency, it's just used to mask aerosol outputs
         call set_daynight_indices(coszrs(1:ncol), day_indices(1:ncol), night_indices(1:ncol))
         nday = count(day_indices(1:ncol) > 0)
         nnight = count(night_indices(1:ncol) > 0)

         ! Loop over diagnostic calls
         call rad_cnst_get_call_list(active_calls)
         do icall = N_DIAG,0,-1
            if (active_calls(icall)) then
               ! Get gas concentrations
               call t_startf('rad_gas_concentrations_sw')
               call get_gas_vmr(icall, state, pbuf, active_gases, gas_vmr)
               call t_stopf('rad_gas_concentrations_sw')
               ! Get aerosol optics
               if (do_aerosol_rad) then
                  call t_startf('rad_aer_optics_sw')
                  aer_tau_bnd_sw = 0._r8
                  aer_ssa_bnd_sw = 0._r8
                  aer_asm_bnd_sw = 0._r8
                  if (do_SPA_optics) then
                     !get current time step's date
                     call get_curr_date(year,month,day,tod)
                     !populate with corresponding pbuf
                     !variables
                     call pbuf_get_field(pbuf,aer_tau_bnd_sw_mon_1_idx, aer_tau_bnd_sw_mon_1)
                     call pbuf_get_field(pbuf,aer_tau_bnd_sw_mon_2_idx, aer_tau_bnd_sw_mon_2)
                     call pbuf_get_field(pbuf,aer_ssa_bnd_sw_mon_1_idx, aer_ssa_bnd_sw_mon_1)
                     call pbuf_get_field(pbuf,aer_ssa_bnd_sw_mon_2_idx, aer_ssa_bnd_sw_mon_2)
                     call pbuf_get_field(pbuf,aer_asm_bnd_sw_mon_1_idx, aer_asm_bnd_sw_mon_1)
                     call pbuf_get_field(pbuf,aer_asm_bnd_sw_mon_2_idx, aer_asm_bnd_sw_mon_2)

                     if (current_month .ne. month) then

                        aer_tau_bnd_sw_mon_1 = aer_tau_bnd_sw_mon_2
                        aer_ssa_bnd_sw_mon_1 = aer_ssa_bnd_sw_mon_2
                        aer_asm_bnd_sw_mon_1 = aer_asm_bnd_sw_mon_2

                        if (month==12) then
                           next_month = 1
                        else
                           next_month = month + 1
                        end if

                        allocate(aerosol_optical_property_sw(pcols,pver,nswbands,begchunk:endchunk))
                        call get_aerosol_optical_property_from_file(next_month,'AER_TAU_SW','swband',nswbands,aerosol_optical_property_sw)
                        aer_tau_bnd_sw_mon_2 = aerosol_optical_property_sw(:,:,:,state%lchnk)
                        call get_aerosol_optical_property_from_file(next_month,'AER_SSA_SW','swband',nswbands,aerosol_optical_property_sw)
                        aer_ssa_bnd_sw_mon_2 = aerosol_optical_property_sw(:,:,:,state%lchnk)
                        call get_aerosol_optical_property_from_file(next_month,'AER_G_SW','swband',nswbands,aerosol_optical_property_sw)
                        aer_asm_bnd_sw_mon_2 = aerosol_optical_property_sw(:,:,:,state%lchnk)

                        deallocate(aerosol_optical_property_sw)

                        current_month = month
                     end if
                     !interpolate between two months to calculate prescribed
                     !aerosol optical property based on current date
                     fraction_of_month = (day*3600.0*24.0 + tod)/(3600*24*days_per_month(current_month)) !tod is in seconds
                     aer_tau_bnd_sw = aer_tau_bnd_sw_mon_1*(1-fraction_of_month) + aer_tau_bnd_sw_mon_2*(fraction_of_month)
                     aer_ssa_bnd_sw = aer_ssa_bnd_sw_mon_1*(1-fraction_of_month) + aer_ssa_bnd_sw_mon_2*(fraction_of_month)
                     aer_asm_bnd_sw = aer_asm_bnd_sw_mon_1*(1-fraction_of_month) + aer_asm_bnd_sw_mon_2*(fraction_of_month)

                  else

                     call set_aerosol_optics_sw( &
                          icall, dt, state, pbuf, &
                          night_indices(1:nnight), &
                          is_cmip6_volc, &
                          aer_tau_bnd_sw, aer_ssa_bnd_sw, aer_asm_bnd_sw &
                  )
                  end if !SPA optics
                  ! Now reorder bands to be consistent with RRTMGP
                  ! TODO: fix the input files themselves!
                  do icol = 1,size(aer_tau_bnd_sw,1)
                     do ilay = 1,size(aer_tau_bnd_sw,2)
                        aer_tau_bnd_sw(icol,ilay,:) = reordered(aer_tau_bnd_sw(icol,ilay,:), rrtmg_to_rrtmgp_swbands)
                        aer_ssa_bnd_sw(icol,ilay,:) = reordered(aer_ssa_bnd_sw(icol,ilay,:), rrtmg_to_rrtmgp_swbands)
                        aer_asm_bnd_sw(icol,ilay,:) = reordered(aer_asm_bnd_sw(icol,ilay,:), rrtmg_to_rrtmgp_swbands)
                     end do
                  end do
                  call t_stopf('rad_aer_optics_sw')
               else
                  aer_tau_bnd_sw = 0
                  aer_ssa_bnd_sw = 0
                  aer_asm_bnd_sw = 0
               end if

               ! Check (and possibly clip) values before passing to RRTMGP driver
               call handle_error(clip_values(cld_tau_gpt_sw,  0._r8, huge(cld_tau_gpt_sw), trim(subname) // ' cld_tau_gpt_sw', tolerance=1e-10_r8))
               call handle_error(clip_values(cld_ssa_gpt_sw,  0._r8,                1._r8, trim(subname) // ' cld_ssa_gpt_sw', tolerance=1e-10_r8))
               call handle_error(clip_values(cld_asm_gpt_sw, -1._r8,                1._r8, trim(subname) // ' cld_asm_gpt_sw', tolerance=1e-10_r8))
               call handle_error(clip_values(aer_tau_bnd_sw,  0._r8, huge(aer_tau_bnd_sw), trim(subname) // ' aer_tau_bnd_sw', tolerance=1e-10_r8))
               call handle_error(clip_values(aer_ssa_bnd_sw,  0._r8,                1._r8, trim(subname) // ' aer_ssa_bnd_sw', tolerance=1e-10_r8))
               call handle_error(clip_values(aer_asm_bnd_sw, -1._r8,                1._r8, trim(subname) // ' aer_asm_bnd_sw', tolerance=1e-10_r8))

               ! Call the shortwave radiation driver
               call radiation_driver_sw( &
                  ncol, active_gases, gas_vmr, &
                  pmid, pint, tmid, albedo_dir, albedo_dif, coszrs, &
                  cld_tau_gpt_sw, cld_ssa_gpt_sw, cld_asm_gpt_sw, &
                  aer_tau_bnd_sw, aer_ssa_bnd_sw, aer_asm_bnd_sw, &
                  fluxes_allsky, fluxes_clrsky, qrs, qrsc &
               )
               ! Send fluxes to history buffer
               call output_fluxes_sw(icall, state, fluxes_allsky, fluxes_clrsky, qrs,  qrsc)
            end if
         end do

         ! Set net fluxes used by other components (land?)
         call set_net_fluxes_sw(fluxes_allsky, fsds, fsns, fsnt)

         ! Set surface fluxes that are used by the land model
         call export_surface_fluxes(fluxes_allsky, cam_out, 'shortwave')

         ! Free memory allocated for shortwave fluxes
         call free_fluxes(fluxes_allsky)
         call free_fluxes(fluxes_clrsky)

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

         call t_startf('rad_cld_optics_lw')
         cld_tau_gpt_lw = 0._r8
         call get_cloud_optics_lw( &
            ncol, pver, nlwbands, do_snow_optics(), cld, cldfsnow, iclwp, iciwp, icswp, &
            lambdac, mu, dei, des, rei, &
            cld_tau_bnd_lw, liq_tau_bnd_lw, ice_tau_bnd_lw, snw_tau_bnd_lw &
         )
         call sample_cloud_optics_lw( &
            ncol, pver, nlwgpts, k_dist_lw%get_gpoint_bands(), &
            state%pmid, cld, cldfsnow, &
            cld_tau_bnd_lw, cld_tau_gpt_lw &
         )
         call output_cloud_optics_lw(state, cld_tau_bnd_lw)
         call t_stopf('rad_cld_optics_lw')

         ! Loop over diagnostic calls
         call rad_cnst_get_call_list(active_calls)
         do icall = N_DIAG,0,-1
            if (active_calls(icall)) then
               ! Get gas concentrations
               call t_startf('rad_gas_concentrations_sw')
               call get_gas_vmr(icall, state, pbuf, active_gases, gas_vmr)
               call t_stopf('rad_gas_concentrations_sw')
               ! Get aerosol optics
               if (do_aerosol_rad) then
                  call t_startf('rad_aer_optics_lw')
                  aer_tau_bnd_lw = 0._r8
                 if (do_SPA_optics) then
                     !get current time step's date
                     call get_curr_date(year,month,day,tod)
                     !populate with corresponding pbuf
                     !variables
                     call pbuf_get_field(pbuf,aer_tau_bnd_lw_mon_1_idx, aer_tau_bnd_lw_mon_1)
                     call pbuf_get_field(pbuf,aer_tau_bnd_lw_mon_2_idx, aer_tau_bnd_lw_mon_2)
                     if (current_month .ne. month) then
                        aer_tau_bnd_lw_mon_1 = aer_tau_bnd_lw_mon_2
                        if (month==12) then
                           next_month = 1
                        else
                           next_month = month + 1
                        end if
                        allocate(aerosol_optical_property_lw(pcols,pver,nlwbands,begchunk:endchunk))
                        call get_aerosol_optical_property_from_file(next_month,'AER_TAU_LW','lwband',nlwbands,aerosol_optical_property_lw)
                        aer_tau_bnd_lw_mon_2 = aerosol_optical_property_lw(:,:,:,state%lchnk)
                        deallocate(aerosol_optical_property_lw)
                        current_month = month
                     end if
                     !interpolate between two months to calculate prescribed
                     !aerosol optical property based on current date
                     fraction_of_month = (day*3600.0*24.0 + tod)/(3600*24*days_per_month(current_month)) !tod is in seconds
                     aer_tau_bnd_lw = aer_tau_bnd_lw_mon_1*(1-fraction_of_month) + aer_tau_bnd_lw_mon_2*(fraction_of_month)
                  else

                     call aer_rad_props_lw(is_cmip6_volc, icall, dt, state, pbuf, aer_tau_bnd_lw)
                  end if
                  call t_stopf('rad_aer_optics_lw')
               else
                  aer_tau_bnd_lw = 0
               end if

               ! Check (and possibly clip) values before passing to RRTMGP driver
               call handle_error(clip_values(cld_tau_gpt_lw,  0._r8, huge(cld_tau_gpt_lw), trim(subname) // ': cld_tau_gpt_lw', tolerance=1e-10_r8))
               call handle_error(clip_values(aer_tau_bnd_lw,  0._r8, huge(aer_tau_bnd_lw), trim(subname) // ': aer_tau_bnd_lw', tolerance=1e-10_r8))

               ! Call the longwave radiation driver to calculate fluxes and heating rates
               call radiation_driver_lw( &
                  ncol, active_gases, gas_vmr, &
                  pmid, pint, tmid, tint, &
                  cld_tau_gpt_lw, aer_tau_bnd_lw, &
                  fluxes_allsky, fluxes_clrsky, qrl, qrlc &
               )
               ! Send fluxes to history buffer
               call output_fluxes_lw(icall, state, fluxes_allsky, fluxes_clrsky, qrl, qrlc)
            end if
         end do

         ! Set net fluxes used in other components
         call set_net_fluxes_lw(fluxes_allsky, flns, flnt)

         ! Export surface fluxes that are used by the land model
         call export_surface_fluxes(fluxes_allsky, cam_out, 'longwave')

         ! Free memory allocated for shortwave fluxes
         call free_fluxes(fluxes_allsky)
         call free_fluxes(fluxes_clrsky)

      else

         ! Conserve energy (what does this mean exactly?)
         if (conserve_energy) then
            qrl(1:ncol,1:pver) = qrl(1:ncol,1:pver) / state%pdel(1:ncol,1:pver)
         end if

      end if  ! dolw

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
      call t_stopf ('rad_heating_rate')
      call outfld('HR', hr(1:ncol,:), ncol, state%lchnk)

      ! convert radiative heating rates to Q*dp for energy conservation
      if (conserve_energy) then
         qrs(1:ncol,1:pver) = qrs(1:ncol,1:pver) * state%pdel(1:ncol,1:pver)
         qrl(1:ncol,1:pver) = qrl(1:ncol,1:pver) * state%pdel(1:ncol,1:pver)
      end if

      ! If we ran radiation this timestep, check if we should run COSP
      if (radiation_do('sw') .or. radiation_do('lw')) then
         if (docosp) then
            ! Advance counter and run COSP if new count value is equal to cosp_nradsteps
            cosp_cnt(state%lchnk) = cosp_cnt(state%lchnk) + 1
            if (cosp_cnt(state%lchnk) == cosp_nradsteps) then

               ! Find bands used for cloud simulator calculations
               cosp_lwband = get_band_index_lw(10.5_r8, 'micron')
               cosp_swband = get_band_index_sw(0.67_r8, 'micron')

               ! Compute quantities needed for COSP
               cld_emis_lw = 1._r8 - exp(-cld_tau_bnd_lw(:,:,cosp_lwband))
               gb_snow_tau_sw = cldfsnow * snw_tau_bnd_sw(:,:,cosp_swband)
               gb_snow_tau_lw = cldfsnow * snw_tau_bnd_lw(:,:,cosp_lwband)

               ! Run COSP; note that "snow" is passed as LW absorption optical depth, even though
               ! the argument name in "snow_emis_in". Optical depth is apparently converted to
               ! emissivity somewhere in the COSP routines. This should probably be changed in the
               ! future for clarity, but is consistent with the implementation in RRTMG.
               call t_startf('cospsimulator_intr_run')
               call cospsimulator_intr_run( &
                  state, pbuf, cam_in, cld_emis_lw, coszrs, &
                  cld_swtau_in=cld_tau_bnd_sw(:,:,cosp_swband), &
                  snow_tau_in=gb_snow_tau_sw, snow_emis_in=gb_snow_tau_lw &
               )
               call t_stopf('cospsimulator_intr_run')

               ! Reset counter
               cosp_cnt(state%lchnk) = 0
            end if
         end if
      end if

   end subroutine radiation_tend

   !----------------------------------------------------------------------------

   subroutine radiation_driver_sw(ncol, &
                                  gas_names, gas_vmr, &
                                  pmid, pint, tmid, albedo_dir, albedo_dif, coszrs, &
                                  cld_tau_gpt, cld_ssa_gpt, cld_asm_gpt, &
                                  aer_tau_bnd, aer_ssa_bnd, aer_asm_bnd, &
                                  fluxes_allsky, fluxes_clrsky, qrs, qrsc)

      use perf_mod, only: t_startf, t_stopf
      use mo_rrtmgp_clr_all_sky, only: rte_sw
      use mo_fluxes_byband, only: ty_fluxes_byband
      use mo_optical_props, only: ty_optical_props_2str
      use mo_gas_concentrations, only: ty_gas_concs
      use radiation_utils, only: calculate_heating_rate
      use cam_optics, only: get_cloud_optics_sw, sample_cloud_optics_sw, &
                            compress_optics_sw, set_aerosol_optics_sw

      ! Inputs
      integer, intent(in) :: ncol
      type(ty_fluxes_byband), intent(inout) :: fluxes_allsky, fluxes_clrsky
      real(r8), intent(inout) :: qrs(:,:), qrsc(:,:)
      character(len=*), intent(in), dimension(:) :: gas_names
      real(r8), intent(in), dimension(:,:,:) :: gas_vmr
      real(r8), intent(in), dimension(:,:) :: pmid, pint, tmid
      real(r8), intent(in), dimension(:,:) :: albedo_dir, albedo_dif
      real(r8), intent(in), dimension(:) :: coszrs
      real(r8), intent(in), dimension(:,:,:) :: cld_tau_gpt, cld_ssa_gpt, cld_asm_gpt
      real(r8), intent(in), dimension(:,:,:) :: aer_tau_bnd, aer_ssa_bnd, aer_asm_bnd

      ! Temporary fluxes compressed to daytime only arrays
      type(ty_fluxes_byband) :: fluxes_allsky_day, fluxes_clrsky_day

      ! Temporary heating rates on radiation vertical grid (and daytime only)
      real(r8), dimension(ncol,nlev_rad) :: qrs_rad, qrsc_rad

      ! Albedo for shortwave calculations
      real(r8), dimension(nswbands,ncol) :: albedo_dir_day, albedo_dif_day

      ! Cloud and aerosol optics
      type(ty_optical_props_2str) :: aer_optics_sw, cld_optics_sw

      ! Gas concentrations
      type(ty_gas_concs) :: gas_concentrations

      ! Incoming solar radiation, scaled for solar zenith angle
      ! and earth-sun distance
      real(r8) :: solar_irradiance_by_gpt(ncol,nswgpts)

      ! Gathered indicies of day and night columns
      ! chunk_column_index = day_indices(daylight_column_index)
      integer :: nday, nnight     ! Number of daylight columns
      integer :: day_indices(ncol), night_indices(ncol)   ! Indicies of daylight coumns

      ! Cosine solar zenith angle for daytime columns
      real(r8) :: coszrs_day(ncol)  ! cosine solar zenith angle

      ! Scaling factor for total sky irradiance; used to account for orbital
      ! eccentricity, and could be used to scale total sky irradiance for different
      ! climates as well (i.e., paleoclimate simulations)
      real(r8) :: tsi_scaling

      ! Loop indices
      integer :: iband

      ! State fields that are passed into RRTMGP. Some of these may need to
      ! modified from what exist in the physics_state object, i.e. to clip
      ! temperatures to make sure they are within the valid range.
      real(r8), dimension(ncol,nlev_rad) :: tmid_day, pmid_day
      real(r8), dimension(ncol,nlev_rad+1) :: pint_day, tint_day

      real(r8), dimension(size(active_gases),ncol,pver) :: gas_vmr_day
      integer :: igas, iday, icol

      ! Everybody needs a name
      character(*), parameter :: subroutine_name = 'radiation_driver_sw'


      if (fixed_total_solar_irradiance<0) then
         ! Get orbital eccentricity factor to scale total sky irradiance
         tsi_scaling = get_eccentricity_factor()
      else
         ! For fixed TSI we divide by the default solar constant of 1360.9
         ! At some point we will want to replace this with a method that
         ! retrieves the solar constant
         tsi_scaling = fixed_total_solar_irradiance / 1360.9_r8
      end if

      ! Gather night/day column indices for subsetting SW inputs; we only want to
      ! do the shortwave radiative transfer during the daytime to save
      ! computational cost (and because RRTMGP will fail for cosine solar zenith
      ! angles less than or equal to zero)
      call set_daynight_indices(coszrs(1:ncol), day_indices(1:ncol), night_indices(1:ncol))
      nday = count(day_indices(1:ncol) > 0)
      nnight = count(night_indices(1:ncol) > 0)

      ! If no daytime columns in this chunk, then we return zeros
      if (nday == 0) then
         call reset_fluxes(fluxes_allsky)
         call reset_fluxes(fluxes_clrsky)
         qrs(1:ncol,1:pver) = 0
         qrsc(1:ncol,1:pver) = 0
         return
      end if

      ! Compress state to daytime-only
      do iday = 1,nday
         icol = day_indices(iday)
         tmid_day(iday,:) = tmid(icol,:)
         pmid_day(iday,:) = pmid(icol,:)
         pint_day(iday,:) = pint(icol,:)
      end do

      ! Compress to daytime-only arrays
      do iband = 1,nswbands
         call compress_day_columns(albedo_dir(iband,1:ncol), albedo_dir_day(iband,1:nday), day_indices(1:nday))
         call compress_day_columns(albedo_dif(iband,1:ncol), albedo_dif_day(iband,1:nday), day_indices(1:nday))
      end do
      call compress_day_columns(coszrs(1:ncol), coszrs_day(1:nday), day_indices(1:nday))

      ! Allocate shortwave fluxes (allsky and clearsky)
      ! TODO: why do I need to provide my own routines to do this? Why is
      ! this not part of the ty_fluxes_byband object?
      !
      ! NOTE: fluxes defined at interfaces, so initialize to have vertical
      ! dimension nlev_rad+1, while we initialized the RRTMGP input variables to
      ! have vertical dimension nlev_rad (defined at midpoints).
      call initialize_rrtmgp_fluxes(nday, nlev_rad+1, nswbands, fluxes_allsky_day, do_direct=.true.)
      call initialize_rrtmgp_fluxes(nday, nlev_rad+1, nswbands, fluxes_clrsky_day, do_direct=.true.)

      ! Populate RRTMGP optics
      call handle_error(cld_optics_sw%alloc_2str(nday, nlev_rad, k_dist_sw, name='shortwave cloud optics'))
      cld_optics_sw%tau = 0
      cld_optics_sw%ssa = 1
      cld_optics_sw%g   = 0
      call compress_optics_sw( &
         day_indices, cld_tau_gpt, cld_ssa_gpt, cld_asm_gpt, &
         cld_optics_sw%tau(1:nday,2:nlev_rad,:), &
         cld_optics_sw%ssa(1:nday,2:nlev_rad,:), &
         cld_optics_sw%g  (1:nday,2:nlev_rad,:)  &
      )
      ! Apply delta scaling to account for forward-scattering
      call handle_error(cld_optics_sw%delta_scale())

      ! Initialize aerosol optics; passing only the wavenumber bounds for each
      ! "band" rather than passing the full spectral discretization object, and
      ! omitting the "g-point" mapping forces the optics to be indexed and
      ! stored by band rather than by g-point. This is most consistent with our
      ! treatment of aerosol optics in the model, and prevents us from having to
      ! map bands to g-points ourselves since that will all be handled by the
      ! private routines internal to the optics class.
      call handle_error(aer_optics_sw%alloc_2str( &
         nday, nlev_rad, k_dist_sw%get_band_lims_wavenumber(), &
         name='shortwave aerosol optics' &
      ))

      ! Populate aerosol optics (and compress to daytime only)
      if (do_aerosol_rad) then
         aer_optics_sw%tau = 0
         aer_optics_sw%ssa = 1
         aer_optics_sw%g   = 0
         call compress_optics_sw( &
            day_indices, &
            aer_tau_bnd(1:ncol,1:pver,:), &
            aer_ssa_bnd(1:ncol,1:pver,:), &
            aer_asm_bnd(1:ncol,1:pver,:), &
            aer_optics_sw%tau(1:nday,2:nlev_rad,:), &
            aer_optics_sw%ssa(1:nday,2:nlev_rad,:), &
            aer_optics_sw%g  (1:nday,2:nlev_rad,:)  &
         )
         ! Apply delta scaling to account for forward-scattering
         call handle_error(aer_optics_sw%delta_scale())
      else
         aer_optics_sw%tau(:,:,:) = 0
         aer_optics_sw%ssa(:,:,:) = 0
         aer_optics_sw%g  (:,:,:) = 0
      end if

      ! Compress gases to daytime-only
      call t_startf('rad_set_gases_sw')
      do igas = 1,size(active_gases)
         call compress_day_columns(gas_vmr(igas,1:ncol,1:pver), &
                                   gas_vmr_day(igas,1:nday,1:pver), &
                                   day_indices(1:nday))
      end do
      call set_gas_concentrations(nday, gas_names, gas_vmr_day, gas_concentrations)
      call t_stopf('rad_set_gases_sw')

      ! Do shortwave radiative transfer calculations
      call t_startf('rad_calculations_sw')
      call handle_error(rte_sw( &
         k_dist_sw, gas_concentrations, &
         pmid_day(1:nday,1:nlev_rad), &
         tmid_day(1:nday,1:nlev_rad), &
         pint_day(1:nday,1:nlev_rad+1), &
         coszrs_day(1:nday), &
         albedo_dir_day(1:nswbands,1:nday), &
         albedo_dif_day(1:nswbands,1:nday), &
         cld_optics_sw, &
         fluxes_allsky_day, fluxes_clrsky_day, &
         aer_props=aer_optics_sw, &
         tsi_scaling=tsi_scaling &
      ))
      call t_stopf('rad_calculations_sw')

      ! Calculate heating rates on the DAYTIME columns
      call t_startf('rad_heating_rate_sw')
      call calculate_heating_rate(      &
         fluxes_allsky_day%flux_up,     &
         fluxes_allsky_day%flux_dn,     &
         pint_day(1:nday,1:nlev_rad+1), &
         qrs_rad(1:nday,1:nlev_rad)     &
      )
      call calculate_heating_rate(      &
         fluxes_clrsky_day%flux_up,     &
         fluxes_clrsky_day%flux_dn,     &
         pint_day(1:nday,1:nlev_rad+1), &
         qrsc_rad(1:nday,1:nlev_rad)    &
      )
      call t_stopf('rad_heating_rate_sw')

      ! Expand fluxes from daytime-only arrays to full chunk arrays
      call t_startf('rad_expand_fluxes_sw')
      call expand_day_fluxes(fluxes_allsky_day, fluxes_allsky, day_indices(1:nday))
      call expand_day_fluxes(fluxes_clrsky_day, fluxes_clrsky, day_indices(1:nday))
      call t_stopf('rad_expand_fluxes_sw')

      ! Expand heating rates to all columns and map back to CAM levels
      call t_startf('rad_expand_heating_rate_sw')
      call expand_day_columns(qrs_rad(1:nday,ktop:kbot), qrs(1:ncol,1:pver), day_indices(1:nday))
      call expand_day_columns(qrsc_rad(1:nday,ktop:kbot), qrsc(1:ncol,1:pver), day_indices(1:nday))
      call t_stopf('rad_expand_heating_rate_sw')

      ! Free fluxes and optical properties
      call free_optics_sw(cld_optics_sw)
      call free_optics_sw(aer_optics_sw)
      call free_fluxes(fluxes_allsky_day)
      call free_fluxes(fluxes_clrsky_day)

   end subroutine radiation_driver_sw

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

   subroutine output_cloud_optics_sw(state, coszrs, cldf, tau, ssa, asm, tau_liq, tau_ice)
      use ppgrid, only: pver
      use physics_types, only: physics_state
      use cam_history, only: outfld
      use radconstants, only: idx_sw_diag
      use cam_history_support, only: fillvalue

      type(physics_state), intent(in) :: state
      real(r8), intent(in), dimension(:) :: coszrs
      real(r8), intent(in), dimension(:,:) :: cldf
      real(r8), intent(in), dimension(:,:,:) :: tau, ssa, asm, tau_liq, tau_ice
      real(r8), dimension(size(tau,1), size(tau,2)) :: &
         tot_cld_vistau, tot_icld_vistau, liq_icld_vistau, ice_icld_vistau
      character(len=*), parameter :: subname = 'output_cloud_optics_sw'
      integer :: sw_band_index, ncol, nlay, icol, ilay

      ncol = state%ncol
      nlay = pver

      ! Check values
      call assert_valid(tau(1:ncol,1:nlay,1:nswbands), &
                        trim(subname) // ': optics%optical_depth')
      call assert_valid(ssa(1:ncol,1:nlay,1:nswbands), &
                        trim(subname) // ': optics%single_scattering_albedo')
      call assert_valid(asm(1:ncol,1:nlay,1:nswbands), &
                        trim(subname) // ': optics%assymmetry_parameter')

      ! Output band-by-band optical properties
      call outfld('CLOUD_TAU_SW', &
                  tau(1:ncol,1:nlay,1:nswbands), &
                  ncol, state%lchnk)
      call outfld('CLOUD_SSA_SW', &
                  ssa(1:ncol,1:nlay,1:nswbands), &
                  ncol, state%lchnk)
      call outfld('CLOUD_G_SW', &
                  asm(1:ncol,1:nlay,1:nswbands), &
                  ncol, state%lchnk)

      ! Compute and output diagnostics for single (550 nm) band
      sw_band_index = get_band_index_sw(550._r8, 'nm')
      do icol = 1,size(tau, 1)
         do ilay = 1,size(tau, 2)
            ! NOTE: this logical is inside the loop to make this play nice with GPU parallelism
            if (coszrs(icol) > 0) then
               tot_icld_vistau(icol,ilay) = tau    (icol,ilay,sw_band_index)
               liq_icld_vistau(icol,ilay) = tau_liq(icol,ilay,sw_band_index)
               ice_icld_vistau(icol,ilay) = tau_ice(icol,ilay,sw_band_index)
               tot_cld_vistau(icol,ilay)  = cldf(icol,ilay) * tot_icld_vistau(icol,ilay)
            else
               tot_icld_vistau(icol,ilay) = fillvalue
               liq_icld_vistau(icol,ilay) = fillvalue
               ice_icld_vistau(icol,ilay) = fillvalue
               tot_cld_vistau(icol,ilay)  = fillvalue
            end if
         end do
      end do
      call outfld('TOT_CLD_VISTAU' ,  tot_cld_vistau(1:ncol,1:nlay), ncol, state%lchnk)
      call outfld('TOT_ICLD_VISTAU', tot_icld_vistau(1:ncol,1:nlay), ncol, state%lchnk)
      call outfld('LIQ_ICLD_VISTAU', liq_icld_vistau(1:ncol,1:nlay), ncol, state%lchnk)
      call outfld('ICE_ICLD_VISTAU', ice_icld_vistau(1:ncol,1:nlay), ncol, state%lchnk)
   end subroutine output_cloud_optics_sw

   !----------------------------------------------------------------------------

   subroutine output_cloud_optics_lw(state, tau)

      use ppgrid, only: pver
      use physics_types, only: physics_state
      use cam_history, only: outfld
      use radconstants, only: rrtmg_lw_cloudsim_band

      type(physics_state), intent(in) :: state
      real(r8), intent(in), dimension(:,:,:) :: tau
      real(r8), dimension(size(tau, 1), size(tau, 2)) :: emis
      integer :: lw_band_index, ncol, nlay

      ncol = state%ncol
      nlay = pver

      ! Check values
      call assert_valid(tau(1:state%ncol,1:pver,1:nlwbands), 'cld_tau_lw')

      ! Output
      call outfld('CLOUD_TAU_LW', tau(1:ncol,1:pver,1:nlwbands), ncol, state%lchnk)
      lw_band_index = get_band_index_lw(10.5_r8, 'micron')
      emis(1:ncol,:) = 1._r8 - exp(-tau(1:ncol,:,lw_band_index))
      call outfld('EMIS', emis(1:ncol,1:pver), ncol, state%lchnk)

   end subroutine output_cloud_optics_lw

   !----------------------------------------------------------------------------

   subroutine radiation_driver_lw(ncol, &
                                  gas_names, gas_vmr, &
                                  pmid, pint, tmid, tint, &
                                  cld_tau_gpt, aer_tau_bnd, &
                                  fluxes_allsky, fluxes_clrsky, qrl, qrlc)

      use perf_mod, only: t_startf, t_stopf
      use mo_rrtmgp_clr_all_sky, only: rte_lw
      use mo_fluxes_byband, only: ty_fluxes_byband
      use mo_optical_props, only: ty_optical_props_1scl
      use mo_gas_concentrations, only: ty_gas_concs
      use radiation_utils, only: calculate_heating_rate

      ! Inputs
      integer, intent(in) :: ncol
      type(ty_fluxes_byband), intent(inout) :: fluxes_allsky, fluxes_clrsky
      real(r8), intent(inout) :: qrl(:,:), qrlc(:,:)
      character(len=*), intent(in), dimension(:) :: gas_names
      real(r8), intent(in), dimension(:,:,:) :: gas_vmr
      real(r8), intent(in), dimension(:,:) :: pmid, pint, tmid, tint
      real(r8), intent(in), dimension(:,:,:) :: cld_tau_gpt, aer_tau_bnd

      ! Everybody needs a name
      character(*), parameter :: subroutine_name = 'radiation_driver_lw'

      ! Surface emissivity needed for longwave
      real(r8) :: surface_emissivity(nlwbands,ncol)

      ! Temporary heating rates on radiation vertical grid
      real(r8), dimension(ncol,nlev_rad) :: qrl_rad, qrlc_rad

      ! RRTMGP types
      type(ty_gas_concs) :: gas_concentrations
      type(ty_optical_props_1scl) :: aer_optics_lw
      type(ty_optical_props_1scl) :: cld_optics_lw

      ! Set surface emissivity to 1 here. There is a note in the RRTMG
      ! implementation that this is treated in the land model, but the old
      ! RRTMG implementation also sets this to 1. This probably does not make
      ! a lot of difference either way, but if a more intelligent value
      ! exists or is assumed in the model we should use it here as well.
      ! TODO: set this more intelligently?
      surface_emissivity(1:nlwbands,1:ncol) = 1.0_r8

      ! Populate RRTMGP optics
      call t_startf('longwave cloud optics')
      call handle_error(cld_optics_lw%alloc_1scl(ncol, nlev_rad, k_dist_lw, name='longwave cloud optics'))
      cld_optics_lw%tau = 0.0
      cld_optics_lw%tau(1:ncol,2:nlev_rad,:) = cld_tau_gpt(1:ncol,1:pver,:)
      call handle_error(cld_optics_lw%delta_scale())
      call t_stopf('longwave cloud optics')

      ! Initialize aerosol optics; passing only the wavenumber bounds for each
      ! "band" rather than passing the full spectral discretization object, and
      ! omitting the "g-point" mapping forces the optics to be indexed and
      ! stored by band rather than by g-point. This is most consistent with our
      ! treatment of aerosol optics in the model, and prevents us from having to
      ! map bands to g-points ourselves since that will all be handled by the
      ! private routines internal to the optics class.
      call handle_error(aer_optics_lw%alloc_1scl(ncol, nlev_rad, k_dist_lw%get_band_lims_wavenumber()))
      call aer_optics_lw%set_name('longwave aerosol optics')

      ! Set gas concentrations (I believe the active gases may change
      ! for different values of icall, which is why we do this within
      ! the loop).
      call t_startf('rad_gas_concentrations_lw')
      call set_gas_concentrations(ncol, gas_names, gas_vmr, gas_concentrations)
      call t_stopf('rad_gas_concentrations_lw')

      if (do_aerosol_rad) then
         ! Get longwave aerosol optics
         call t_startf('rad_aer_optics_lw')
         aer_optics_lw%tau(:,:,:) = 0
         aer_optics_lw%tau(1:ncol,ktop:kbot,1:nlwbands) = aer_tau_bnd(1:ncol,1:pver,1:nlwbands)
         ! Apply delta scaling to account for forward-scattering
         call handle_error(aer_optics_lw%delta_scale())
         call t_stopf('rad_aer_optics_lw')
      else
         aer_optics_lw%tau(:,:,:) = 0
      end if

      ! Do longwave radiative transfer calculations
      call t_startf('rad_calculations_lw')
      call handle_error(rte_lw( &
         k_dist_lw, gas_concentrations, &
         pmid(1:ncol,1:nlev_rad), tmid(1:ncol,1:nlev_rad), &
         pint(1:ncol,1:nlev_rad+1), tint(1:ncol,nlev_rad+1), &
         surface_emissivity(1:nlwbands,1:ncol), &
         cld_optics_lw, &
         fluxes_allsky, fluxes_clrsky, &
         aer_props=aer_optics_lw, &
         t_lev=tint(1:ncol,1:nlev_rad+1), &
         n_gauss_angles=1 & ! Set to 3 for consistency with RRTMG
      ))
      call t_stopf('rad_calculations_lw')

      ! Calculate heating rates
      call calculate_heating_rate(  &
         fluxes_allsky%flux_up,     &
         fluxes_allsky%flux_dn,     &
         pint(1:ncol,1:nlev_rad+1), &
         qrl_rad(1:ncol,1:nlev_rad) &
      )
      call calculate_heating_rate(   &
         fluxes_clrsky%flux_up,      &
         fluxes_clrsky%flux_dn,      &
         pint(1:ncol,1:nlev_rad+1),  &
         qrlc_rad(1:ncol,1:nlev_rad) &
      )

      ! Map heating rates to CAM columns and levels
      qrl(1:ncol,1:pver) = qrl_rad(1:ncol,ktop:kbot)
      qrlc(1:ncol,1:pver) = qrlc_rad(1:ncol,ktop:kbot)

      ! Free fluxes and optical properties
      call free_optics_lw(cld_optics_lw)
      call free_optics_lw(aer_optics_lw)

   end subroutine radiation_driver_lw

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

   !----------------------------------------------------------------------------

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

   !----------------------------------------------------------------------------

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

   !----------------------------------------------------------------------------

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
      real(r8) :: wavenumber_limits(2,nswbands)
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
      wavenumber_limits(:,:) = k_dist_sw%get_band_lims_wavenumber()

      ! Loop over bands, and determine for each band whether it is broadly in the
      ! visible or infrared part of the spectrum (visible or "not visible")
      do iband = 1,nswbands
         if (is_visible(wavenumber_limits(1,iband)) .and. &
             is_visible(wavenumber_limits(2,iband))) then

            ! Entire band is in the visible
            albedo_dir(iband,1:ncol) = cam_in%asdir(1:ncol)
            albedo_dif(iband,1:ncol) = cam_in%asdif(1:ncol)

         else if (.not.is_visible(wavenumber_limits(1,iband)) .and. &
                  .not.is_visible(wavenumber_limits(2,iband))) then

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
      ! was never checked for RRTMG, so albedos will probably be slightly different
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

      ! TOA fluxes (above model top, use index to rad top)
      call outfld('FSUTOA'//diag(icall), flux_all%flux_up(1:ncol,ktop_rad), ncol, state%lchnk)
      call outfld('FSNTOA'//diag(icall), flux_all%flux_net(1:ncol,ktop_rad), ncol, state%lchnk)

      ! Clear-sky fluxes
      call outfld('FSNTC'//diag(icall), flux_clr%flux_net(1:ncol,ktop), ncol, state%lchnk)
      call outfld('FSNSC'//diag(icall), flux_clr%flux_net(1:ncol,kbot+1), ncol, state%lchnk)
      call outfld('FSDSC'//diag(icall), flux_clr%flux_dn(1:ncol,kbot+1), ncol, state%lchnk)
      call outfld('FSUTOAC'//diag(icall), flux_clr%flux_up(1:ncol,ktop_rad), ncol, state%lchnk)
      call outfld('FSNTOAC'//diag(icall), flux_clr%flux_net(1:ncol,ktop_rad), ncol, state%lchnk)

      ! Calculate and output the shortwave cloud radiative effect (SWCF in history)
      cloud_radiative_effect(1:ncol) = flux_all%flux_net(1:ncol,ktop_rad) - flux_clr%flux_net(1:ncol,ktop_rad)
      call outfld('SWCF'//diag(icall), cloud_radiative_effect, ncol, state%lchnk)

      ! Send solar insolation to history buffer
      call outfld('SOLIN'//diag(icall), flux_clr%flux_dn(1:ncol,1), ncol, state%lchnk)

      ! Send heating rates to history buffer
      call outfld('QRS'//diag(icall), qrs(1:ncol,1:pver)/cpair, ncol, state%lchnk)
      call outfld('QRSC'//diag(icall), qrsc(1:ncol,1:pver)/cpair, ncol, state%lchnk)

   end subroutine output_fluxes_sw

   !----------------------------------------------------------------------------

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

   !----------------------------------------------------------------------------

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

   !----------------------------------------------------------------------------

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

   !----------------------------------------------------------------------------

   subroutine free_optics_sw(optics)
      use mo_optical_props, only: ty_optical_props_2str
      type(ty_optical_props_2str), intent(inout) :: optics
      if (allocated(optics%tau)) deallocate(optics%tau)
      if (allocated(optics%ssa)) deallocate(optics%ssa)
      if (allocated(optics%g)) deallocate(optics%g)
      call optics%finalize()
   end subroutine free_optics_sw

   !----------------------------------------------------------------------------

   subroutine free_optics_lw(optics)
      use mo_optical_props, only: ty_optical_props_1scl
      type(ty_optical_props_1scl), intent(inout) :: optics
      if (allocated(optics%tau)) deallocate(optics%tau)
      call optics%finalize()
   end subroutine free_optics_lw

   !----------------------------------------------------------------------------

   subroutine get_gas_vmr(icall, state, pbuf, gas_names, gas_vmr)

      use physics_types, only: physics_state
      use physics_buffer, only: physics_buffer_desc
      use rad_constituents, only: rad_cnst_get_gas

      integer, intent(in) :: icall
      type(physics_state), intent(in) :: state
      type(physics_buffer_desc), pointer :: pbuf(:)
      character(len=*), intent(in), dimension(:) :: gas_names
      real(r8), intent(out), dimension(:,:,:) :: gas_vmr

      ! Mass mixing ratio
      real(r8), pointer :: mmr(:,:)

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
      real(r8), parameter :: co_vol_mix_ratio = 1.0e-7_r8
      real(r8), parameter :: n2_vol_mix_ratio = 0.7906_r8

      ! Loop indices
      integer :: igas

      ! Number of columns
      integer :: ncol

      ! Name of subroutine for error messages
      character(len=32) :: subname = 'get_gas_vmr'

      ! Number of columns in chunk
      ncol = state%ncol

      ! initialize
      gas_vmr(:,:,:) = 0._r8

      ! For each gas species needed for RRTMGP, read the mass mixing ratio from the
      ! CAM rad_constituents interface, convert to volume mixing ratios, and
      ! subset for daytime-only indices if needed.
      do igas = 1,size(gas_names)

         select case(trim(gas_names(igas)))

            case('CO')

               ! CO not available, use default
               gas_vmr(igas,1:ncol,1:pver) = co_vol_mix_ratio

            case('N2')

               ! N2 not available, use default
               gas_vmr(igas,1:ncol,1:pver) = n2_vol_mix_ratio

            case('H2O')

               ! Water vapor is represented as specific humidity in CAM, so we
               ! need to handle water a little differently
               call rad_cnst_get_gas(icall, trim(gas_species(igas)), state, pbuf, mmr)

               ! Convert to volume mixing ratio by multiplying by the ratio of
               ! molecular weight of dry air to molecular weight of gas. Note that
               ! first specific humidity (held in the mass_mix_ratio array read
               ! from rad_constituents) is converted to an actual mass mixing
               ! ratio.
               gas_vmr(igas,1:ncol,1:pver) = mmr(1:ncol,1:pver) / ( &
                  1._r8 - mmr(1:ncol,1:pver) &
               )  * mol_weight_air / mol_weight_gas(igas)

            case DEFAULT

               ! Get mass mixing ratio from the rad_constituents interface
               call rad_cnst_get_gas(icall, trim(gas_species(igas)), state, pbuf, mmr)

               ! Convert to volume mixing ratio by multiplying by the ratio of
               ! molecular weight of dry air to molecular weight of gas
               gas_vmr(igas,1:ncol,1:pver) = mmr(1:ncol,1:pver) &
                                            * mol_weight_air / mol_weight_gas(igas)

         end select

      end do  ! igas

   end subroutine get_gas_vmr

   !----------------------------------------------------------------------------

   subroutine set_gas_concentrations(ncol, gas_names, gas_vmr, gas_concentrations)
      use mo_gas_concentrations, only: ty_gas_concs
      use mo_rrtmgp_util_string, only: lower_case

      integer, intent(in) :: ncol
      character(len=*), intent(in), dimension(:) :: gas_names
      real(r8), intent(in), dimension(:,:,:) :: gas_vmr
      type(ty_gas_concs), intent(out) :: gas_concentrations

      ! Local variables
      real(r8), dimension(pcols,nlev_rad) :: vol_mix_ratio_out

      ! Loop indices
      integer :: igas

      ! Character array to hold lowercase gas names
      character(len=32), allocatable :: gas_names_lower(:)

      ! Name of subroutine for error messages
      character(len=32) :: subname = 'set_gas_concentrations'

      ! Initialize gas concentrations with lower case names
      allocate(gas_names_lower(size(gas_names)))
      do igas = 1,size(gas_names)
         gas_names_lower(igas) = trim(lower_case(gas_names(igas)))
      end do
      call handle_error(gas_concentrations%init(gas_names_lower))

      ! For each gas, add level above model top and set values in RRTMGP object
      do igas = 1,size(gas_names)
         vol_mix_ratio_out = 0
         ! Map to radiation grid
         vol_mix_ratio_out(1:ncol,ktop:kbot) = gas_vmr(igas,1:ncol,1:pver)
         ! Copy top-most model level to top-most rad level (which could be above
         ! the top of the model)
         vol_mix_ratio_out(1:ncol,1) = gas_vmr(igas,1:ncol,1)
         ! Set volumn mixing ratio in gas concentration object for just columns
         ! in this chunk
         call handle_error(gas_concentrations%set_vmr( &
            trim(lower_case(gas_names(igas))), vol_mix_ratio_out(1:ncol,1:nlev_rad)) &
         )
      end do

   end subroutine set_gas_concentrations

   !----------------------------------------------------------------------------

   logical function string_in_list(string, list)
      character(len=*), intent(in) :: string
      character(len=*), intent(in) :: list(:)
      integer :: list_index

      ! Set to false initially
      string_in_list = .false.

      ! Loop over items in list, and if we find a match then set flag to true and
      ! exit loop
      do list_index = 1,size(list)
         if (trim(list(list_index)) == trim(string)) then
            string_in_list = .true.
            exit
         end if
      end do
   end function string_in_list

   !----------------------------------------------------------------------------

   subroutine expand_day_fluxes(daytime_fluxes, expanded_fluxes, day_indices)
      use mo_rte_kind, only: wp
      use mo_fluxes_byband, only: ty_fluxes_byband
      type(ty_fluxes_byband), intent(in) :: daytime_fluxes
      type(ty_fluxes_byband), intent(inout) :: expanded_fluxes
      integer, intent(in) :: day_indices(:)

      integer :: nday, iday, icol

      ! Reset fluxes in expanded_fluxes object to zero
      call reset_fluxes(expanded_fluxes)

      ! Number of daytime columns is number of indices greater than zero
      nday = count(day_indices > 0)

      ! Loop over daytime indices and map daytime fluxes into expanded arrays
      do iday = 1,nday

         ! Map daytime index to proper column index
         icol = day_indices(iday)

         ! Expand broadband fluxes
         expanded_fluxes%flux_up(icol,:) = daytime_fluxes%flux_up(iday,:)
         expanded_fluxes%flux_dn(icol,:) = daytime_fluxes%flux_dn(iday,:)
         expanded_fluxes%flux_net(icol,:) = daytime_fluxes%flux_net(iday,:)
         if (associated(daytime_fluxes%flux_dn_dir)) then
            expanded_fluxes%flux_dn_dir(icol,:) = daytime_fluxes%flux_dn_dir(iday,:)
         end if

         ! Expand band-by-band fluxes
         expanded_fluxes%bnd_flux_up(icol,:,:) = daytime_fluxes%bnd_flux_up(iday,:,:)
         expanded_fluxes%bnd_flux_dn(icol,:,:) = daytime_fluxes%bnd_flux_dn(iday,:,:)
         expanded_fluxes%bnd_flux_net(icol,:,:) = daytime_fluxes%bnd_flux_net(iday,:,:)
         if (associated(daytime_fluxes%bnd_flux_dn_dir)) then
            expanded_fluxes%bnd_flux_dn_dir(icol,:,:) = daytime_fluxes%bnd_flux_dn_dir(iday,:,:)
         end if

      end do

   end subroutine expand_day_fluxes

   ! Should we do snow optics? Check for existence of "cldfsnow" variable
   logical function do_snow_optics()
      use physics_buffer, only: pbuf_get_index
      use cam_abortutils, only: endrun
      real(r8), pointer :: pbuf(:)
      integer :: err, idx

      idx = pbuf_get_index('CLDFSNOW', errcode=err)
      if (idx > 0) then
         do_snow_optics = .true.
      else
         do_snow_optics = .false.
      end if

      return
   end function do_snow_optics

end module radiation
