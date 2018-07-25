module radiation
!---------------------------------------------------------------------------------
! Purpose:
!
! CAM interface to RRTMGP
!
! Revision history:
! Nov 2017  B. Hillman  Initial version, adapted from RRTMG interface.
!
! # NOTES
! 
! ## How to clean up this module
! This is a bit of a mess right now, but the code runs and computes both clear
! and cloudy sky fluxes using RRTMGP. There are some differences to sort out
! with RRTMG, but this can be done while cleaning up the interface to enable
! adding support for the super-parameterized configuration. In doing so, we
! would like to make this interface flexible enough to allow using for
! calculated fluxes at the CRM level...that is, calling some form of this driver
! from within the CRM itself to calculate the radiative heating on both the CRM
! spatial and time domain so that we can update the CRM heating at each CRM
! timestep. Instead, the current behavior (implemented in RRTMG) uses
! the time-averaged CRM state to calculate some kind of averaged heating rate,
! and then applies that same heating rate at each CRM timestep the next time the
! CRM is called.
!
! In order to allow using this with the CRM, either at this level or at the CRM
! level, we can try to contain the RRTMGP interface to a single subroutine that
! takes as arguments either just simple arrays, or derived types we know about
! (i.e., cam_fluxes, cam_optics), and doing whatever conversions are necessary
! within that routine. This could all be contained within the cam_fluxes class,
! so maybe this would look something like:
!
!    cam_fluxes_sw%initialize()
!    cam_cloud_optics_sw%initialize()
!    cam_aerosol_optics_sw%initialize()
!    cam_fluxes_sw%calculate(k_dist_sw, cam_cloud_optics_sw, cam_aerosol_optics_sw)
! 
! It would then be straightforward to implement routines to populate a
! cam_optics object with CRM-scale inputs, and use cam_fluxes%calculate to
! calculate fluxes on the CRM domain.
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
use mo_gas_optics, only: ty_gas_optics
use mo_rte_kind, only: wp

! Use my assertion routines to perform sanity checks
use assertions, only: assert, assert_valid, assert_range

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

! Private module data
! TODO: remove these?
integer :: cldfsnow_idx = 0 

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
type(ty_gas_optics) :: k_dist_sw, k_dist_lw

! k-distribution coefficients files to read from. These are set via namelist
! variables.
character(len=cl) :: coefficients_file_sw, coefficients_file_lw

! Number of shortwave and longwave bands in use by the RRTMGP radiation code.
! This information will be stored in the k_dist_sw and k_dist_lw objects and may
! be retrieved using the k_dist_sw%get_nband() and k_dist_lw%get_nband()
! methods, but I think we need to save these as private module data so that we
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

! Number of vertical levels in radiation calculations. This will generally
! include an extra layer above the model top (num_rad_layers = pver + 1).
integer :: num_rad_layers

! Indices to map levels on the radiation grid to levels on the cam grid
integer, allocatable :: cam_layers(:), cam_interfaces(:)

! Set name for this module (for purpose of writing output and log files)
character(len=*), parameter :: module_name = 'radiation'

! Gases that we want to use in the radiative calculations. These need to be set
! here, because the RRTMGP kdist initialization needs to know the names of the
! gases before these are available via the rad_cnst interface. So, if we want to
! change which gases are used at some point, this needs to be changed here, or
! we need to come up with a more clever way of getting these from the model,
! since this is kind of a hack to hard-code them in here.
!character(len=3), dimension(8) :: active_gases = (/ &
!   'H2O', 'CO2', 'O3 ', 'N2O', &
!   'CO ', 'CH4', 'O2 ', 'N2 ' &
!/)
character(len=3), dimension(6) :: active_gases = (/ &
   'H2O', 'CO2', 'O3 ', 'N2O', &
   'CH4', 'O2 ' &
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

! Mapping from old RRTMG sw bands to new band ordering in RRTMGP
integer, dimension(14) :: map_rrtmg_to_rrtmgp_swbands = (/ &
   14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 &
/)

! Interface blocks for overloaded procedures
interface clip_values
   module procedure clip_values_1d, clip_values_2d
end interface clip_values

interface compress_day_columns
   module procedure compress_day_columns_1d, compress_day_columns_2d
end interface compress_day_columns
!===============================================================================

contains

!===============================================================================

subroutine radiation_readnl(nlfile, dtime_in)
!-------------------------------------------------------------------------------
! Purpose: Read radiation_nl namelist group.
!-------------------------------------------------------------------------------

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, &
                              mpi_integer, mpi_logical, &
                              mpi_character, masterproc
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
   namelist /radiation_nl/ rrtmgp_coefficients_file_lw, &
                           rrtmgp_coefficients_file_sw, &
                           iradsw, iradlw, irad_always, &
                           use_rad_dt_cosz, spectralflux

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
   ! TODO: do we need to broadcast? Only if masterproc?
   call mpibcast(rrtmgp_coefficients_file_lw, cl, mpi_character, mstrid, mpicom, ierr)
   call mpibcast(rrtmgp_coefficients_file_sw, cl, mpi_character, mstrid, mpicom, ierr)
   call mpibcast(iradsw, 1, mpi_integer, mstrid, mpicom, ierr)
   call mpibcast(iradlw, 1, mpi_integer, mstrid, mpicom, ierr)
   call mpibcast(irad_always, 1, mpi_integer, mstrid, mpicom, ierr)
   call mpibcast(use_rad_dt_cosz, 1, mpi_logical, mstrid, mpicom, ierr)
   call mpibcast(spectralflux, 1, mpi_logical, mstrid, mpicom, ierr)
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
                      use_rad_dt_cosz, spectralflux
   end if
10 format('  LW coefficents file: ',                                a/, &
          '  SW coefficents file: ',                                a/, &
          '  Frequency (timesteps) of Shortwave Radiation calc:  ',i5/, &
          '  Frequency (timesteps) of Longwave Radiation calc:   ',i5/, &
          '  SW/LW calc done every timestep for first N steps. N=',i5/, &
          '  Use average zenith angle:                           ',l5/, &
          '  Output spectrally resolved fluxes:                  ',l5/)

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

   ! Chemistry interface needs shortwave down at surface
   ! TODO: why add the others? Are these needed?
   call pbuf_add_field('FSDS', 'global', dtype_r8, (/pcols/), idx)
   call pbuf_add_field('FSNS', 'global', dtype_r8, (/pcols/), idx)
   call pbuf_add_field('FSNT', 'global', dtype_r8, (/pcols/), idx)
   call pbuf_add_field('FLNS', 'global', dtype_r8, (/pcols/), idx)
   call pbuf_add_field('FLNT', 'global', dtype_r8, (/pcols/), idx)

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
   use physconst,          only: gravit, epsilo, stebol, &
                                 pstd, mwdry, mwco2, mwo3
   use phys_control,       only: phys_getopts
   use rad_constituents,   only: N_DIAG, rad_cnst_get_call_list, rad_cnst_get_info
   use cospsimulator_intr, only: docosp, cospsimulator_intr_init
   use hirsbt,             only: hirsbt_init
   use hirsbtpar,          only: hirsname, msuname
   use modal_aer_opt,      only: modal_aer_opt_init
   use time_manager,       only: get_nstep, get_step_size, is_first_restart_step
   use radiation_data,     only: init_rad_data
   use physics_types, only: physics_state

   ! RRTMGP modules
   use mo_load_coefficients, only: rrtmgp_load_coefficients=>load_and_init
   use mo_gas_concentrations, only: ty_gas_concs

   ! For optics
   use cloud_rad_props, only: cloud_rad_props_init

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

   logical :: use_SPCAM  ! SPCAM flag

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

   ! Likewise for g-points
   nswgpts = k_dist_sw%get_ngpt()
   nlwgpts = k_dist_lw%get_ngpt()

   ! Set number of levels used in radiation calculations
   num_rad_layers = pver + 1

   ! Setup a mapping from the RRTMGP grid back to the CAM grid
   allocate(cam_layers(num_rad_layers), cam_interfaces(num_rad_layers+1))
   cam_layers(:) = get_cam_layers(num_rad_layers, pver)
   cam_interfaces(:) = get_cam_layers(num_rad_layers+1, pver+1)

   ! Set the radiation timestep for cosz calculations if requested using the 
   ! adjusted iradsw value from radiation
   if (use_rad_dt_cosz)  then
      dtime  = get_step_size()
      dt_avg = iradsw*dtime
   end if

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
   call add_hist_coord('rad_lay', num_rad_layers, 'Radiation layers', '1')
   call add_hist_coord('rad_lev', num_rad_layers+1, 'Radiation layers', '1')

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
   call addfld('CLOUD_TAU_SW', (/'lev','swband'/), 'I', '1', &
               'Cloud extinction optical depth', sampling_seq='rad_lwsw') 
   call addfld('CLOUD_TAU_LW', (/'lev','lwband'/), 'I', '1', &
               'Cloud absorption optical depth', sampling_seq='rad_lwsw') 

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
         call addfld('FUS'//diag(icall),  (/ 'ilev' /),  'I',  'W/m2', &
                     'Shortwave upward flux')
         call addfld('FDS'//diag(icall),  (/ 'ilev' /),  'I',  'W/m2', &
                     'Shortwave downward flux')
         call addfld('FDS_DIR'//diag(icall),  (/ 'ilev' /),  'I',  'W/m2', &
                     'Shortwave direct-beam downward flux')
         call addfld('FNS'//diag(icall),  (/ 'ilev' /),  'I',  'W/m2', &
                     'Shortwave net flux')
         call addfld('FUSC'//diag(icall),  (/ 'ilev' /), 'I',  'W/m2', &
                     'Shortwave clear-sky upward flux')
         call addfld('FDSC'//diag(icall),  (/ 'ilev' /), 'I',  'W/m2', &
                     'Shortwave clear-sky downward flux')
         call addfld('FDSC_DIR'//diag(icall),  (/ 'ilev' /), 'I',  'W/m2', &
                     'Shortwave clear-sky direct-beam downward flux')
         call addfld('FNSC'//diag(icall),  (/ 'ilev' /), 'I',  'W/m2', &
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
         call addfld('FUL'//diag(icall),     (/'ilev'/), 'I', 'W/m2', &
                     'Longwave upward flux')
         call addfld('FDL'//diag(icall),     (/'ilev'/), 'I', 'W/m2', &
                     'Longwave downward flux')
         call addfld('FNL'//diag(icall),     (/'ilev'/), 'I', 'W/m2', &
                     'Longwave net flux')
         call addfld('FULC'//diag(icall),    (/'ilev'/), 'I', 'W/m2', &
                     'Longwave clear-sky upward flux')
         call addfld('FDLC'//diag(icall),    (/'ilev'/), 'I', 'W/m2', &
                     'Longwave clear-sky downward flux')
         call addfld('FNLC'//diag(icall),    (/'ilev'/), 'I', 'W/m2', &
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

   cldfsnow_idx = pbuf_get_index('CLDFSNOW',errcode=err)
   if (cldfsnow_idx > 0) then
      call addfld('CLDFSNOW',(/ 'lev' /),'I','1','CLDFSNOW',flag_xyfill=.true.)
      call addfld('SNOW_ICLD_VISTAU', (/ 'lev' /), 'A', '1', &
                  'Snow in-cloud extinction visible sw optical depth', &
                  sampling_seq='rad_lwsw', flag_xyfill=.true.)
   endif

   ! Inputs to radiation for debugging
   call addfld('RAD_TMID', (/'rad_lay'/), 'I', 'K', &
               'Temperature input to radiation', &
               flag_xyfill=.true., sampling_seq='rad_lwsw')
   call addfld('RAD_PMID', (/'rad_lay'/), 'I', 'Pa', &
               'Pressure input to radiation', &
               flag_xyfill=.true., sampling_seq='rad_lwsw')
   call addfld('RAD_TINT', (/'rad_lev'/), 'I', 'K', &
               'Temperature input to radiation', &
               flag_xyfill=.true., sampling_seq='rad_lwsw')
   call addfld('RAD_PINT', (/'rad_lev'/), 'I', 'Pa', &
               'Pressure input to radiation', &
               flag_xyfill=.true., sampling_seq='rad_lwsw')
   call addfld('RAD_COLDRY', (/'rad_lay'/), 'I', '1', &
               'COLDRY input to radiation', &
               flag_xyfill=.true., sampling_seq='rad_lwsw')
   call addfld('RAD_H2OVMR', (/'rad_lay'/), 'I', '1', &
               'H2O volumn mixing ratio input to radiation', &
               flag_xyfill=.true., sampling_seq='rad_lwsw')
   
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

#else
   write(iulog,*) subname // ': PERGRO not implemented for RRTMGP, doing nothing.'
#endif /* DO_PERGRO_MODS */

end subroutine perturbation_growth_init


subroutine set_available_gases(gases, gas_concentrations)

   use mo_gas_concentrations, only: ty_gas_concs
   use mo_util_string, only: lower_case
   use mo_rte_kind, only: wp

   type(ty_gas_concs), intent(inout) :: gas_concentrations
   character(len=*), intent(in) :: gases(:)
   integer :: i

   ! Use set_vmr method to set gas names. This just sets the vmr to zero, since
   ! this routine is just used to build a list of available gas names.
   do i = 1,size(gases)
      call handle_error(gas_concentrations%set_vmr(trim(lower_case(gases(i))), 0._wp))
   end do

end subroutine set_available_gases


! Function to calculate band midpoints from kdist band limits
function get_band_midpoints(nbands, kdist) result(band_midpoints)
   integer, intent(in) :: nbands
   type(ty_gas_optics), intent(in) :: kdist
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
  
subroutine radiation_tend(state, ptend, pbuf, cam_out, cam_in, net_flux)
!-------------------------------------------------------------------------------
! 
! Purpose: 
! Driver for radiation computation.
! 
! Method: 
! Radiation uses cgs units, so conversions must be done from
! model fields to radiation fields.
!
! Revision history:
! May 2004    D.B. Coleman     Merge of code from radctl.F90 and parts of tphysbc.F90.
! 2004-08-09  B. Eaton         Add pointer variables for constituents.
! 2004-08-24  B. Eaton         Access O3 and GHG constituents from chem_get_cnst.
! 2004-08-30  B. Eaton         Replace chem_get_cnst by rad_constituent_get.
! 2007-11-05  M. Iacono        Install rrtmg_lw and sw as radiation model.
! 2007-12-27  M. Iacono        Modify to use CAM cloud optical properties with rrtmg.
! 2009-06     Minghuai Wang,   add treatments for MMF CAM
!              These modifications are based on the spcam3.5, which was developled
!              by Marat Khairoutdinov (mkhairoutdin@ms.cc.sunysb.edu). The spcam3.5 
!              only have one radiation package (camrt). 
!              Short wave and long wave radiation codes are called for every column
!              of the CRM domain. CRM fields are named as *_crm, and domain-averaged fields are
!              named as *_m. The domain-averaged fields and CRM fields are outputed at the end
!              of the CRM domain loop (last_column=.true.).
!              Several variables in state are updated with those from CRM output
!              (liquid water, qc_rad; ice water, qi_rad; water vapor, qv_rad; 
!              and temperature, t_rad). Several variables in pbuf are also updated 
!              with those in CRM domain (cld, cicewp, cliqwp, csnowp, cldfsnow).  
!              At the end of the radiation calculation, state and those in pbuf are
!              restored to the old values. 
!              Finally, a new cloud simulator are called, which takes account of cloud fileds 
!              from the CRM domain. 
! 2009-07-13, Minghuai Wang: MMF CAM
!             When Morrison's two momenent microphysics is used in SAM, droplet and ice crystal effective radius
!             used in this radiation code are calcualted at each CRM column by calling m2005_effradius
! 2009-10-21, Minghuai Wang: MMF CAM
!             CRM-scale aerosol water is used to calculate aerosol optical depth
! 2017-11     Ben Hillman: Adapted for use with RRTMGP
!----------------------------------------------------------------------------------------------

   ! Performance module needed for timing functions
   use perf_mod, only: t_startf, t_stopf

   ! CAM derived types; needed for surface exchange fields, physics state, and
   ! tendency fields
   use camsrfexch, only: cam_out_t, cam_in_t
   use physics_types, only: physics_state, physics_ptend

   ! Utilities for interacting with the physics buffer
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
                             pbuf_get_index, pbuf_old_tim_idx

   ! For calculating radiative heating tendencies
   use radheat, only: radheat_tend

   ! For getting radiative constituent gases
   use rad_constituents, only: N_DIAG, rad_cnst_get_call_list

   ! Index to visible channel for diagnostic outputs
   use radconstants, only: idx_sw_diag

   ! RRTMGP radiation drivers and derived types
   use mo_rrtmgp_clr_all_sky, only: rte_sw, rte_lw
   use mo_gas_concentrations, only: ty_gas_concs
   use mo_optical_props, only: ty_optical_props, &
                               ty_optical_props_1scl, &
                               ty_optical_props_2str
   use mo_fluxes_byband, only: ty_fluxes_byband

   ! CAM history module provides subroutine to send output data to the history
   ! buffer to be aggregated and written to disk
   use cam_history, only: outfld

   ! CAM optical properties; includes cam_optics_type class for holding optical
   ! properties, and subroutines to get CAM aerosol and cloud optical properties
   ! via CAM parameterizations
   use cam_optics, only: cam_optics_type
   use physconst, only: cpair, stebol

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

   ! ---------------------------------------------------------------------------
   ! Local variables
   ! ---------------------------------------------------------------------------
   integer :: nstep  ! current timestep number

   ! Pointers to heating rates on physics buffer
   real(r8), pointer :: qrs(:,:) => null()  ! shortwave radiative heating rate 
   real(r8), pointer :: qrl(:,:) => null()  ! longwave  radiative heating rate 
   real(r8), pointer :: fsns(:,:) => null()
   real(r8), pointer :: fsnt(:,:) => null()
   real(r8), pointer :: flns(:,:) => null()
   real(r8), pointer :: flnt(:,:) => null()

   ! Clear-sky heating rates are not on the physics buffer, and we have no
   ! reason to put them there, so declare these are regular arrays here
   real(r8) :: qrsc(pcols,pver) = 0._r8
   real(r8) :: qrlc(pcols,pver) = 0._r8

   ! Pointer to hold gridbox-mean cloud fraction (nday,nlev)
   real(r8), pointer :: cloud_fraction(:,:)

   ! Flag to carry (QRS,QRL)*dp across time steps. 
   ! TODO: what does this mean?
   logical :: conserve_energy = .true.

   ! Loop variables; i is used to loop over columns, and k is used to loop over
   ! levels
   integer :: i, k, k_cam
   integer :: icall

   ! nlay is number of layers used in radiation calculation, nlev is number of
   ! interface levels between those layers
   integer :: ncol, nlay, nlev

   ! Gathered indicies of day and night columns 
   ! chunk_column_index = day_indices(daylight_column_index)
   integer :: nday, nnight     ! Number of daylight columns
   integer :: day_indices(pcols), night_indices(pcols)   ! Indicies of daylight coumns

   integer :: model_top

   character(*), parameter :: subroutine_name = 'radiation_tend'

   ! For loops over diagnostic calls (TODO: what does this mean?)
   logical :: active_calls(0:N_DIAG)

   ! Cosine solar zenith angle for all columns in chunk
   real(r8) :: coszrs(pcols)

   ! Incoming solar radiation, scaled for solar zenith angle 
   ! and earth-sun distance
   real(r8) :: solar_irradiance_by_gpt(pcols,nswgpts)

   ! Blackbody temperature
   real(r8) :: t_sfc(pcols)

   ! State fields that are passed into RRTMGP. Some of these may need to
   ! modified from what exist in the physics_state object, i.e. to clip
   ! temperatures to make sure they are within the valid range.
   real(r8) :: coszrs_day(pcols)  ! cosine solar zenith angle
   real(r8), dimension(pcols,num_rad_layers) :: tmid, tmid_day, pmid, pmid_day, tint_cam
   real(r8), dimension(pcols,num_rad_layers+1) :: pint, pint_day, tint
   real(r8) :: surface_emissivity(nlwbands,pcols)

   ! Scaling factor for total sky irradiance; used to account for orbital
   ! eccentricity, and could be used to scale total sky irradiance for different
   ! climates as well (i.e., paleoclimate simulations)
   real(r8) :: tsi_scaling

   ! Temporary heating rates on radiation vertical grid
   real(r8), dimension(pcols,num_rad_layers) :: qrs_rad, qrsc_rad, &
                                                qrl_rad, qrlc_rad

   ! Albedo for shortwave calculations
   real(r8) :: albedo_direct(nswbands,pcols), albedo_direct_day(nswbands,pcols)
   real(r8) :: albedo_diffuse(nswbands,pcols), albedo_diffuse_day(nswbands,pcols)

   ! RRTMGP types
   type(ty_gas_concs) :: gas_concentrations_sw, gas_concentrations_lw
   type(ty_optical_props_1scl) :: aerosol_optics_lw
   type(ty_optical_props_1scl) :: cloud_optics_lw
   type(ty_optical_props_2str) :: aerosol_optics_sw
   type(ty_optical_props_2str) :: cloud_optics_sw
   type(ty_fluxes_byband) :: flux_lw_allsky, flux_lw_clearsky
   type(ty_fluxes_byband) :: flux_sw_allsky, flux_sw_clearsky
   type(ty_fluxes_byband) :: flux_sw_allsky_day, flux_sw_clearsky_day

   !----------------------------------------------------------------------
   
   ! Number of physics columns in this "chunk"; used in multiple places
   ! throughout this subroutine, so set once for convenience
   ncol = state%ncol

   ! Set pointers to heating rates stored on physics buffer. These will be
   ! modified in this routine.
   call pbuf_get_field(pbuf, pbuf_get_index('QRS'), qrs)
   call pbuf_get_field(pbuf, pbuf_get_index('QRL'), qrl)
  
   ! Set pointer to cloud fraction; this is used by McICA routines
   ! TODO: why the extra arguments to pbuf_get_field here? Are these necessary?
   call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction, &
                       start=(/1,1,pbuf_old_tim_idx()/), &
                       kount=(/pcols,pver,1/))

   ! Figure out which vertical levels to use; only use those between the min and
   ! max of press_ref read from the coefficients file. I think the implementation 
   ! that Brian Eaton at NCAR used was to look at the reference pressure pref. We
   ! could also look at the min and max over all columns in this chunk. That
   ! might be the most transparent thing to do. Note that this code here finds
   ! the indices to the top and bottom interfaces, not model midpoints. Note
   ! also that in the CESM implementation (and also in the RRTMG implementation)
   ! an extra layer is added above the RRTMGP topmost layer. I am not 100% sure
   ! why this is, but I think maybe it is to get "top of atmosphere" vs "top of
   ! model" heating and fluxes? TODO.
   !
   ! Size of grid for radiation; note kbot > ktop because CAM grid goes from TOA
   ! to surface. RRTMGP grid will go from surface to TOA internally, but this is
   ! handled automatically so we do not worry about this here. Note that nlay is
   ! the number of *layers* used, so the number of *interfaces* should be nlay+1
   ! TODO: set more intelligently
   nlay = num_rad_layers !pver + 1

   ! Get cosine solar zenith angle for current time step. 
   call set_cosine_solar_zenith_angle(state, dt_avg, coszrs(:ncol))

   ! Send values for this chunk to history buffer
   call outfld('COSZRS', coszrs(:ncol), ncol, state%lchnk)

   ! Get orbital eccentricity factor to scale total sky irradiance
   tsi_scaling = get_eccentricity_factor()

   ! If the swrad_off flag is set, meaning we should not do SW radiation, then 
   ! we just set coszrs to zero everywhere. TODO: why not just set dosw false 
   ! and skip the loop?
   if (swrad_off) coszrs(:) = 0._r8

   ! Gather night/day column indices for subsetting SW inputs; we only want to
   ! do the shortwave radiative transfer during the daytime to save
   ! computational cost (and because RRTMGP will fail for cosine solar zenith
   ! angles less than or equal to zero)
   call set_daynight_indices(coszrs(:ncol), day_indices(:ncol), night_indices(:ncol))
   nday = count(day_indices(:ncol) > 0)
   nnight = count(night_indices(:ncol) > 0)

   ! Do shortwave stuff...
   if (radiation_do('sw') .and. nday > 0) then

      ! Populate RRTMGP input variables. Use the day_indices index array to
      ! map CAM variables on all columns to the daytime-only arrays, and take
      ! only the ktop:kbot vertical levels (mapping CAM vertical grid to
      ! RRTMGP vertical grid).

      ! Map CAM to rad levels
      call map_cam_to_rad_levels(state%t(:ncol,:pver), &
                                 tmid(:ncol,:nlay), &
                                 cam_layers(:nlay)) 
      call map_cam_to_rad_levels(state%pmid(:ncol,:pver), &
                                 pmid(:ncol,:nlay), &
                                 cam_layers(:nlay)) 
      call map_cam_to_rad_levels(state%pint(:ncol,:pver+1), &
                                 pint(:ncol,:nlay+1), &
                                 cam_interfaces(:nlay+1)) 

      ! Fill layer above model top; this is done 
      ! consistent with the RRTMG implementation
      if (nlay > pver) then
         tmid(:ncol,1) = state%t(:ncol,1)
         pmid(:ncol,1) = 0.5_r8 * state%pint(:ncol,1)
         pint(:ncol,1) = 1.01_r8
      end if

      ! Get albedo. This uses CAM routines internally and just provides a
      ! wrapper to improve readability of the code here.
      call set_albedo(cam_in, albedo_direct(:nswbands,:ncol), albedo_diffuse(:nswbands,:ncol))

      ! Send albedos to history buffer (useful for debugging)
      call outfld('SW_ALBEDO_DIR', transpose(albedo_direct(:nswbands,:ncol)), ncol, state%lchnk)
      call outfld('SW_ALBEDO_DIF', transpose(albedo_diffuse(:nswbands,:ncol)), ncol, state%lchnk)

      ! Compress to daytime-only arrays
      do i = 1,nswbands
         call compress_day_columns(albedo_direct(i,:ncol), albedo_direct_day(i,:nday), day_indices(:nday))
         call compress_day_columns(albedo_diffuse(i,:ncol), albedo_diffuse_day(i,:nday), day_indices(:nday))
      end do
      call compress_day_columns(tmid(:ncol,:nlay), tmid_day(:nday,:nlay), day_indices(:nday))
      call compress_day_columns(pmid(:ncol,:nlay), pmid_day(:nday,:nlay), day_indices(:nday))
      call compress_day_columns(pint(:ncol,:nlay+1), pint_day(:nday,:nlay+1), day_indices(:nday))
      call compress_day_columns(coszrs(:ncol), coszrs_day(:nday), day_indices(:nday))

      ! Allocate shortwave fluxes (allsky and clearsky)
      ! TODO: why do I need to provide my own routines to do this? Why is 
      ! this not part of the ty_fluxes_byband object?
      !
      ! NOTE: fluxes defined at interfaces, so initialize to have vertical
      ! dimension nlay+1, while we initialized the RRTMGP input variables to
      ! have vertical dimension nlay (defined at midpoints).
      call initialize_rrtmgp_fluxes(nday, nlay+1, nswbands, flux_sw_allsky_day, do_direct=.true.)
      call initialize_rrtmgp_fluxes(nday, nlay+1, nswbands, flux_sw_clearsky_day, do_direct=.true.)
      call initialize_rrtmgp_fluxes(ncol, nlay+1, nswbands, flux_sw_allsky, do_direct=.true.)
      call initialize_rrtmgp_fluxes(ncol, nlay+1, nswbands, flux_sw_clearsky, do_direct=.true.)

      ! Do shortwave cloud optics calculations
      call t_startf('shortwave cloud optics')
      call handle_error(cloud_optics_sw%alloc_2str(nday, nlay, k_dist_sw))
      call cloud_optics_sw%set_name('shortwave cloud optics')

      ! TODO: refactor the set_cloud_optics codes to allow passing arrays
      ! rather than state/pbuf so that we can use this for superparameterized
      ! simulations...or alternatively add logic within the set_cloud_optics
      ! routines to handle this.
      call set_cloud_optics_sw(state, pbuf, &
                               day_indices(:nday), cam_layers(:nlay), &
                               k_dist_sw, cloud_optics_sw)
      call t_stopf('shortwave cloud optics')

      ! Initialize aerosol optics; passing only the wavenumber bounds for each
      ! "band" rather than passing the full spectral discretization object, and
      ! omitting the "g-point" mapping forces the optics to be indexed and
      ! stored by band rather than by g-point. This is most consistent with our
      ! treatment of aerosol optics in the model, and prevents us from having to
      ! map bands to g-points ourselves since that will all be handled by the
      ! private routines internal to the optics class.
      call handle_error(aerosol_optics_sw%alloc_2str(nday, nlay, k_dist_sw%get_band_lims_wavenumber()))
      call aerosol_optics_sw%set_name('shortwave aerosol optics')

      ! Loop over diagnostic calls 
      ! TODO: more documentation on what this means
      !
      ! NOTE: the climate (icall==0) calculation must occur last, so we loop
      ! backwards.
      call rad_cnst_get_call_list(active_calls)
      do icall = N_DIAG,0,-1
         if (active_calls(icall)) then

            ! Get shortwave aerosol optics
            ! NOTE: this gets vertical dimension nlay since cloud optics are
            ! defined at model midpoints (representing optics through the whole
            ! layer)
            call t_startf('shortwave aerosol optics')
            call set_aerosol_optics_sw(icall, state, pbuf, &
                                       day_indices(:nday), &
                                       night_indices(:nnight), &
                                       cam_layers(:nlay), &
                                       aerosol_optics_sw)
            call t_stopf('shortwave aerosol optics')

            ! Set gas concentrations (I believe the gases may change for
            ! different values of icall, which is why we do this within the
            ! loop)
            call set_gas_concentrations(icall, state, pbuf, nlay, &
                                        gas_concentrations_sw, &
                                        cam_layers(:nlay), &
                                        day_indices=day_indices(:nday))

            ! Do shortwave radiative transfer calculations
            call t_startf('shortwave radiation calculations')
            call handle_error(rte_sw( &
               k_dist_sw, gas_concentrations_sw, &
               pmid_day(:nday,:nlay), &
               tmid_day(:nday,:nlay), &
               pint_day(:nday,:nlay+1), &
               coszrs_day(:nday), &
               albedo_direct_day(:nswbands,:nday), &
               albedo_diffuse_day(:nswbands,:nday), &
               cloud_optics_sw, &
               flux_sw_allsky_day, flux_sw_clearsky_day, &
               aer_props=aerosol_optics_sw, &
               tsi_scaling=tsi_scaling &
            ))
            call t_stopf('shortwave radiation calculations')

            ! Expand fluxes from daytime-only arrays to full chunk arrays
            call expand_daytime_fluxes(flux_sw_allsky_day, flux_sw_allsky, day_indices(:nday))
            call expand_daytime_fluxes(flux_sw_clearsky_day, flux_sw_clearsky, day_indices(:nday))

            ! Calculate heating rates
            call calculate_heating_rate(flux_sw_allsky, pint(:ncol,:nlay+1), qrs_rad(:ncol,:nlay))
            call calculate_heating_rate(flux_sw_clearsky, pint(:ncol,:nlay+1), qrsc_rad(:ncol,:nlay))

            ! Map heating rates to CAM columns and levels
            call map_rad_to_cam_levels(qrs_rad(:ncol,:nlay), qrs(:ncol,:pver), cam_layers(:nlay))
            call map_rad_to_cam_levels(qrsc_rad(:ncol,:nlay), qrsc(:ncol,:pver), cam_layers(:nlay))

            ! Send solar insolation to history buffer
            call outfld('SOLIN'//diag(icall), flux_sw_clearsky%flux_dn(:ncol,1), ncol, state%lchnk)
                        
            ! Send heating rates to history buffer
            call outfld('QRS'//diag(icall), qrs(:ncol,:pver)/cpair, ncol, state%lchnk)
            call outfld('QRSC'//diag(icall), qrsc(:ncol,:pver)/cpair, ncol, state%lchnk)

            ! Send fluxes to history buffer
            call output_fluxes_sw(icall, state, &
                                  flux_sw_allsky, flux_sw_clearsky, &
                                  cam_interfaces(:nlay+1))

            ! Copy outputs to pbuf fields that are used by other model
            ! components (i.e., land)
            call copy_fluxes_to_pbuf(pbuf, flux_sw_allsky, 'shortwave')

         end if
      end do

   else
      ! Conserve energy
      if (conserve_energy) then
         qrs(:ncol,:pver) = qrs(:ncol,:pver) / state%pdel(:ncol,:pver)
      end if
   end if  ! dosw

   ! Do longwave stuff...
   if (radiation_do('lw')) then

      ! Initialize aerosol optics; passing only the wavenumber bounds for each
      ! "band" rather than passing the full spectral discretization object, and
      ! omitting the "g-point" mapping forces the optics to be indexed and
      ! stored by band rather than by g-point. This is most consistent with our
      ! treatment of aerosol optics in the model, and prevents us from having to
      ! map bands to g-points ourselves since that will all be handled by the
      ! private routines internal to the optics class.
      call handle_error(aerosol_optics_lw%alloc_1scl(ncol, nlay, k_dist_lw%get_band_lims_wavenumber()))
      call aerosol_optics_lw%set_name('longwave aerosol optics')

      ! Initialize cloud optics object; cloud_optics_lw will be indexed by
      ! g-point, rather than by band, and subcolumn routines will associate each
      ! g-point with a stochastically-sampled cloud state
      call handle_error(cloud_optics_lw%alloc_1scl(ncol, nlay, k_dist_lw))
      call cloud_optics_lw%set_name('longwave cloud optics')

      ! Allocate longwave outputs; why is this not part of the
      ! ty_fluxes_byband object?
      ! NOTE: fluxes defined at interfaces, so initialize to have vertical
      ! dimension nlay+1
      call initialize_rrtmgp_fluxes( &
         ncol, nlay+1, nlwbands, &
         flux_lw_allsky &
      )
      call initialize_rrtmgp_fluxes( &
         ncol, nlay+1, nlwbands, &
         flux_lw_clearsky &
      )
        
      ! Calculate interface temperature explicitly
      call set_interface_temperature(state, cam_in, tint_cam)

      ! Populate state variables
      call map_cam_to_rad_levels(state%t(:ncol,:pver), &
                                 tmid(:ncol,:nlay), &
                                 cam_layers(:nlay)) 
      call map_cam_to_rad_levels(state%pmid(:ncol,:pver), &
                                 pmid(:ncol,:nlay), &
                                 cam_layers(:nlay)) 
      call map_cam_to_rad_levels(state%pint(:ncol,:pver+1), &
                                 pint(:ncol,:nlay+1), &
                                 cam_interfaces(:nlay+1)) 
      call map_cam_to_rad_levels(tint_cam(:ncol,:pver+1), &
                                 tint(:ncol,:nlay+1), &
                                 cam_interfaces(:nlay+1)) 

      ! Special case for top-most rad level, which is an extra level above
      ! model top; fill temperatures in extra level with temperatures from
      ! model top, set interface pressure at top of extra level to be a small
      ! number, and set midpoint pressure to be half of the pressure at the
      ! top most model interface.
      if (nlay > pver) then
         tmid(:ncol,1) = state%t(:ncol,1)
         pmid(:ncol,1) = 0.5 * state%pint(:ncol,1)
         pint(:ncol,1) = 1.01_r8
      end if

      ! Blackbody temperature of surface
      t_sfc(:ncol) = sqrt(sqrt(cam_in%lwup(:ncol)/stebol))

      ! Do longwave cloud optics calculations
      call t_startf('longwave cloud optics')
      call set_cloud_optics_lw(state, pbuf, cam_layers(:nlay), k_dist_lw, &
                               cloud_optics_lw)
      call t_stopf('longwave cloud optics')

      ! Set surface emissivity to 1 here. There is a note in the RRTMG
      ! implementation that this is treated in the land model, but the old
      ! RRTMG implementation also sets this to 1. This probably does not make
      ! a lot of difference either way, but if a more intelligent value
      ! exists or is assumed in the model we should use it here as well.
      ! TODO: set this more intelligently?
      surface_emissivity(:nlwbands,:ncol) = 1.0_r8

      ! get list of active radiation calls
      call rad_cnst_get_call_list(active_calls)

      ! Loop over diagnostic calls (what does this mean?)
      do icall = N_DIAG,0,-1
         if (active_calls(icall)) then

            ! Set gas concentrations (I believe the active gases may change
            ! for different values of icall, which is why we do this within
            ! the loop).
            call t_startf('longwave gas concentrations')
            call set_gas_concentrations(icall, state, pbuf, nlay, &
                                        gas_concentrations_lw, &
                                        cam_layers(:nlay))
            call t_stopf('longwave gas concentrations')

            ! Get longwave aerosol optics
            call t_startf('longwave aerosol optics')
            call set_aerosol_optics_lw(icall, state, pbuf, cam_layers(:nlay), aerosol_optics_lw)
            call t_stopf('longwave aerosol optics')

            ! Do longwave radiative transfer calculations
            call t_startf('longwave radiation calculations')
            call handle_error(rte_lw( &
               k_dist_lw, gas_concentrations_lw, &
               pmid(1:ncol,1:nlay), tmid(1:ncol,1:nlay), &
               pint(1:ncol,1:nlay+1), t_sfc(1:ncol), &
               surface_emissivity(1:nlwbands,1:ncol), &
               cloud_optics_lw, &
               flux_lw_allsky, flux_lw_clearsky, &
               aer_props=aerosol_optics_lw, &
               t_lev=tint(1:ncol,1:nlay+1) &
            ))
            call t_stopf('longwave radiation calculations')

            ! Calculate heating rates
            call calculate_heating_rate(flux_lw_allsky, pint(:ncol,:nlay+1), &
                                        qrl_rad(:ncol,:nlay))
            call calculate_heating_rate(flux_lw_clearsky, pint(:ncol,:nlay+1), &
                                        qrlc_rad(:ncol,:nlay))

            ! Map heating rates to CAM columns and levels
            call map_rad_to_cam_levels(qrl_rad(:ncol,:nlay), qrl(:ncol,:pver), cam_layers(:nlay))
            call map_rad_to_cam_levels(qrlc_rad(:ncol,:nlay), qrlc(:ncol,:pver), cam_layers(:nlay))
                        
            ! Send heating rates to history buffer
            call outfld('QRL'//diag(icall), qrl(:ncol,:pver)/cpair, ncol, state%lchnk)
            call outfld('QRLC'//diag(icall), qrlc(:ncol,:pver)/cpair, ncol, state%lchnk)

            ! Copy fluxes to pbuf
            call copy_fluxes_to_pbuf(pbuf, flux_lw_allsky, 'longwave')

            ! Send fluxes to history buffer
            call output_fluxes_lw(icall, state, flux_lw_allsky, flux_lw_clearsky, cam_interfaces(:nlay+1))

         end if  ! active calls
      end do  ! loop over diagnostic calls
            
   else
      ! Conserve energy (what does this mean exactly?)
      if (conserve_energy) then
         qrl(:ncol,:pver) = qrl(:ncol,:pver) / state%pdel(:ncol,:pver)
      end if
   end if  ! dolw

   ! Compute net radiative heating tendency
   call t_startf('radheat_tend')
   call pbuf_get_field(pbuf, pbuf_get_index('FSNS'), fsns)
   call pbuf_get_field(pbuf, pbuf_get_index('FSNT'), fsnt)
   call pbuf_get_field(pbuf, pbuf_get_index('FLNS'), flns)
   call pbuf_get_field(pbuf, pbuf_get_index('FLNT'), flnt)
   call radheat_tend(state, pbuf, ptend, &
                     qrl, qrs, &
                     fsns, fsnt, flns, flnt, &
                     cam_in%asdir, net_flux)
   call t_stopf('radheat_tend')

   ! convert radiative heating rates to Q*dp for energy conservation
   ! TODO: this gets converted back here? what is going on here?
   if (conserve_energy) then
      do k = 1,pver
         do i = 1,ncol
            qrs(i,k) = qrs(i,k)*state%pdel(i,k)
            qrl(i,k) = qrl(i,k)*state%pdel(i,k)
         end do
      end do
   end if

   ! copy surface fluxes to cam_out structure
   call export_surface_fluxes(flux_sw_allsky, cam_out, 'shortwave')
   call export_surface_fluxes(flux_lw_allsky, cam_out, 'longwave')

end subroutine radiation_tend


subroutine output_rad_state(state, tmid, tint, pmid, pint, coldry, h2ovmr)

   use physics_types, only: physics_state
   use cam_history, only: outfld
   type(physics_state), intent(in) :: state
   real(r8), intent(in) :: tmid(:,:), tint(:,:), pmid(:,:), pint(:,:)
   real(r8), intent(in), optional :: coldry(:,:), h2ovmr(:,:)
   integer :: ncol, lchnk

   ncol = state%ncol
   lchnk = state%lchnk
   call outfld('RAD_TMID', tmid(:ncol,:num_rad_layers), ncol, lchnk)
   call outfld('RAD_PMID', pmid(:ncol,:num_rad_layers), ncol, lchnk)
   call outfld('RAD_TINT', tint(:ncol,:num_rad_layers+1), ncol, lchnk)
   call outfld('RAD_PINT', pint(:ncol,:num_rad_layers+1), ncol, lchnk)
   if (present(coldry)) then
      call outfld('RAD_COLDRY', coldry(:ncol,:num_rad_layers), ncol, lchnk)
   end if
   if (present(h2ovmr)) then
      call outfld('RAD_H2OVMR', h2ovmr(:ncol,:num_rad_layers), ncol, lchnk)
   end if

end subroutine output_rad_state


! Old calculation from RRTMG...seems to be inconsistent with internal RRTMGP
! function to calculation dry column amount
subroutine set_coldry(p_lev, h2ovmr, coldry)

   use shr_const_mod, only: shr_const_avogad
   use physconst, only: grav=>gravit

   real(r8), intent(in) :: p_lev(:,:)
   real(r8), intent(in) :: h2ovmr(:,:)
   real(r8), intent(inout) :: coldry(:,:)

   real(r8), parameter :: amd = 28.9660_r8      ! Effective molecular weight of dry air (g/mol)
   real(r8), parameter :: amw = 18.0160_r8      ! Molecular weight of water vapor (g/mol)
   real(r8) :: amm

   ! CIME Avogadro's constant in units of molecules / kmole ... convert to the
   ! normal molecules / mole here
   real(r8), parameter :: avogad = shr_const_avogad * 1.e-3_r8

   real(r8) :: dp
   integer :: i, l

   call assert(all(shape(h2ovmr) == shape(coldry)), 'set_coldry: inconsistent shapes')
   call assert(all(shape(p_lev) == (/size(h2ovmr, 1), size(h2ovmr, 2)+1/)), &
               'set_coldry: p_lev inconsistent shape')
   call assert_range(p_lev, 1._r8, 1.e10_r8, 'set_coldry: p_lev out of bounds.')

   coldry(:,:) = 0._r8
   do i = 1,size(h2ovmr, 1)
      do l = 1,size(h2ovmr, 2)

         ! Mass of air (mass of dry air plus mass of water vapor)?
         ! (g/mol)
         amm = (1._r8 - h2ovmr(i,l)) * amd + h2ovmr(i,l) * amw

         ! Pressure thickness (Pa)
         dp = abs(p_lev(i,l+1) - p_lev(i,l))

         ! Hydrostatic equation
         coldry(i,l) = 1.e-2_r8 * dp * 1.e3_r8 * avogad / (1.e2_r8 * grav * amm * (1._r8 + h2ovmr(i,l)))
         if (coldry(i,l) < 0 .or. coldry(i,l) >= huge(coldry)) then
            print *, 'coldry invalid at ', i, l, coldry(i,l)
            print *, '    dp = ', dp
            print *, '    h2ovmr = ', h2ovmr(i,l)
            print *, '    10._r8 * dp * avogad = ', 10._r8 * dp * avogad
            print *, '    amm = ', amm
            call endrun('coldry')
         end if
      enddo

      ! This is weird...top level uses interface pressure from the bottom
      ! interface of the top level, and the h2o volumn mixing ratio of the
      ! second level? Why? This seems to handle the extra level above model
      ! top...dp here is taken to be the pressure difference between the last
      ! real model interface and ZERO, which is why p_lev appears by itself.
      ! h2ovmr should be the same for the top TWO levels, based on how we set
      ! h2ovmr in set_gas_concentrations, so we should be able to take
      ! h2ovmr(i,1) or h2ovmr(i,2) with no difference here. Note that we take
      ! p_lev(i,2) instead of p_lev(i,1) here because p_lev(i,1) is the "fake"
      ! top interface level that we added above model top.
      coldry(i,1) = 1.e-2_r8 * (p_lev(i,2)) * 1.e3_r8 * avogad / (1.e2_r8 * grav * amm * (1._r8 + h2ovmr(i,2)))
   end do

   call assert(all(coldry >= 0), 'set_coldry: coldry has negative values')
   call assert_valid(coldry, 'set_coldry: coldry has inf or NaNs')
end subroutine set_coldry


subroutine set_interface_temperature(state, cam_in, tint)

   use physics_types, only: physics_state
   use camsrfexch, only: cam_in_t
   use physconst, only: stebol

   type(physics_state), intent(in) :: state
   type(cam_in_t), intent(in) :: cam_in
   real(r8), intent(inout) :: tint(:,:)
   real(r8) :: dy
   integer :: i, k

   ! Calculate interface temperatures (following method used in radtpl for the
   ! longwave), using surface upward flux and stebol constant in mks units.
   ! NOTE: this code copied from RRTMG implementation, and DIFFERS from what is
   ! done in RRTMGP if interface temperatures are omitted! Letting RRTMGP handle
   ! this leads to large differences in near-surface fluxes with RRTMG. TODO:
   ! why is this? Are we doing something wrong here? Maybe it's just the surface
   ! temperature (tint(:,pverp)) that is to blame?
   do i = 1,state%ncol
      tint(i,1) = state%t(i,1)
      tint(i,pverp) = sqrt(sqrt(cam_in%lwup(i)/stebol))
      do k = 2,pver
         dy = (state%lnpint(i,k) - state%lnpmid(i,k)) / (state%lnpmid(i,k-1) - state%lnpmid(i,k))
         tint(i,k) = state%t(i,k) - dy * (state%t(i,k) - state%t(i,k-1))
      end do
   end do

end subroutine set_interface_temperature


subroutine export_surface_fluxes(fluxes, cam_out, band)
   use mo_fluxes_byband, only: ty_fluxes_byband
   use camsrfexch, only: cam_out_t

   type(ty_fluxes_byband), intent(in) :: fluxes
   type(cam_out_t), intent(inout) :: cam_out
   character(len=*), intent(in) :: band
   integer :: i, kbot

   ! TODO: check this code! This seems to differ from what is in RRTMG. Also, I
   ! think RRTMGP lets us break output into direct and diffuse beam now, so much
   ! of this might be able to be simplified.
   !
   ! Export surface fluxes
   !
   ! The RRTMGP output isn't broken down into direct and diffuse.  For first cut
   ! put entire flux into direct and leave the diffuse set to zero.  To break the
   ! fluxes into the UV/vis and near-IR bands use the same scheme as for the albedos
   ! which is hardcoded for 14 spectral bands.
   !
   ! sols(pcols)      Direct solar rad on surface (< 0.7)
   ! soll(pcols)      Direct solar rad on surface (>= 0.7)
   !
   ! Near-IR bands (1-9 and 14), 820-16000 cm-1, 0.625-12.195 microns
   !
   ! Put half of band 9 in each of the UV/visible and near-IR values,
   ! since this band straddles 0.7 microns:
   !
   ! UV/visible bands 10-13, 16000-50000 cm-1, 0.200-0.625 micron
   kbot = size(fluxes%flux_dn, 2)
   if (trim(band) == 'shortwave') then
      cam_out%soll = 0
      cam_out%sols = 0
      cam_out%solld = 0
      cam_out%solsd = 0
      do i = 1,size(fluxes%bnd_flux_dn, 1)
         cam_out%soll(i) &
            = sum(fluxes%bnd_flux_dn(i,kbot,1:8)) &
            + 0.5_r8 * fluxes%bnd_flux_dn(i,kbot,9) &
            + fluxes%bnd_flux_dn(i,kbot,14)

         cam_out%sols(i) &
            = 0.5_r8 * fluxes%bnd_flux_dn(i,kbot,9) &
            + sum(fluxes%bnd_flux_dn(i,kbot,10:13))

         cam_out%netsw(i) = fluxes%flux_net(i,kbot)
      end do
   else if (trim(band) == 'longwave') then
      do i = 1,size(fluxes%flux_dn, 1)
         cam_out%flwds(i) = fluxes%flux_dn(i,kbot)
      end do
   else
      call endrun('flux_type ' // band // ' not known.')
   end if

end subroutine export_surface_fluxes


subroutine copy_fluxes_to_pbuf(pbuf, fluxes, band)
   use physics_buffer, only: physics_buffer_desc, pbuf_get_field, &
                             pbuf_get_index
   use mo_fluxes_byband, only: ty_fluxes_byband

   type(physics_buffer_desc), pointer :: pbuf(:)
   type(ty_fluxes_byband), intent(in) :: fluxes
   character(len=*), intent(in) :: band

   ! Pointers to pbuf fields
   real(r8), pointer :: fsds(:)
   real(r8), pointer :: fsns(:)
   real(r8), pointer :: fsnt(:)
   real(r8), pointer :: flns(:)
   real(r8), pointer :: flnt(:)

   integer :: ncol, nlevels, ktop, kbot
   character(len=32) :: subname = 'copy_fluxes_to_pbuf'

   ncol = size(fluxes%flux_up, 1)
   nlevels = size(fluxes%flux_up, 2)
   ktop = 1
   kbot = nlevels
   if (trim(band) == 'shortwave') then
      ! Associate pointers
      call pbuf_get_field(pbuf, pbuf_get_index('FSDS'), fsds)
      call pbuf_get_field(pbuf, pbuf_get_index('FSNS'), fsns)
      call pbuf_get_field(pbuf, pbuf_get_index('FSNT'), fsnt)

      ! Copy data
      fsds(:ncol) = fluxes%flux_dn(:ncol,kbot)
      fsns(:ncol) = fluxes%flux_dn(:ncol,kbot) - fluxes%flux_up(:ncol,kbot) !fluxes%flux_net(:ncol,kbot)
      fsnt(:ncol) = fluxes%flux_dn(:ncol,ktop) - fluxes%flux_up(:ncol,kbot) !fluxes%flux_net(:ncol,ktop)

      ! Check values
      call assert_valid(fsds(:ncol), 'fsds')
      call assert_valid(fsns(:ncol), 'fsns')
      call assert_valid(fsnt(:ncol), 'fsnt')
   else if (trim(band) == 'longwave') then
      ! Associate pointers
      call pbuf_get_field(pbuf, pbuf_get_index('FLNS'), flns)
      call pbuf_get_field(pbuf, pbuf_get_index('FLNT'), flnt)

      ! Copy data
      flns(:ncol) = fluxes%flux_up(:ncol,kbot) - fluxes%flux_dn(:ncol,kbot) !fluxes%flux_net(:ncol,kbot)
      flnt(:ncol) = fluxes%flux_up(:ncol,ktop) - fluxes%flux_dn(:ncol,ktop) !fluxes%flux_net(:ncol,ktop)
   else
      call endrun(trim(subname) // ': flux_type ' // trim(band) // ' not known.')
   end if

end subroutine copy_fluxes_to_pbuf

!-------------------------------------------------------------------------------
! Map fields from CAM vertical grid to radiation vertical grid, which may
! contain more or less vertical levels than CAM (more when we add an extra layer
! above the top of the model, less when not all pressures are above the minimum
! threshold required by RRTMGP lookup tables).
subroutine map_cam_to_rad_levels(xcam, xrad, indices)
   real(r8), intent(in) :: xcam(:,:)
   real(r8), intent(inout) :: xrad(:,:)
   integer, intent(in) :: indices(:)
   integer :: krad, kcam
   character(len=32) :: subname = 'map_cam_to_rad_levels'

   ! Check sizes
   call assert(size(xcam, 1) == size(xrad, 1), &
               trim(subname) // ': dim 1 inconsistent.')
   call assert(size(xrad, 2) == size(indices), &
               trim(subname) // ': dim 2 inconsistent.')

   do krad = 1,size(indices)
      kcam = max(1, indices(krad))
      xrad(:,krad) = xcam(:,kcam)
   end do
end subroutine map_cam_to_rad_levels
!-------------------------------------------------------------------------------
! Map fields from radiation vertical grid back to CAM vertical grid. Uses the
! same indices as above subroutine, which specify the CAM levels corresponding
! to the radiation levels. I.e., indices(1) = 0 means that the first radiation
! level corresponds to a level above the top of the CAM grid.
subroutine map_rad_to_cam_levels(xrad, xcam, indices)
   real(r8), intent(in) :: xrad(:,:)
   real(r8), intent(inout) :: xcam(:,:)
   integer, intent(in) :: indices(:)
   integer :: krad, kcam
   character(len=32) :: subname = 'map_rad_to_cam_levels'

   ! Check sizes
   call assert(size(xcam, 1) == size(xrad, 1), &
               trim(subname) // ': dim 1 inconsistent.')
   call assert(size(xrad, 2) == size(indices), &
               trim(subname) // ': dim 2 inconsistent.')

   do krad = 1,size(indices)
      kcam = max(1, indices(krad))
      xcam(:,kcam) = xrad(:,krad)
   end do
end subroutine map_rad_to_cam_levels
!-------------------------------------------------------------------------------

subroutine calculate_heating_rate(fluxes, pint, heating_rate)

   use physconst, only: gravit
   use mo_fluxes_byband, only: ty_fluxes_byband

   ! Inputs
   type(ty_fluxes_byband), intent(in) :: fluxes
   real(r8), intent(in) :: pint(:,:)

   ! Output heating rate; same size as pint with one fewer vertical level
   real(r8), intent(out) :: heating_rate(:,:)

   ! Loop indices
   integer :: i, k

   ! Everyone needs a name
   character(len=32) :: subname = 'calculate_heating_rate'

   ! Get dimension sizes and make sure arrays conform
   call assert(size(pint,1) == size(fluxes%flux_up,1), subname // ': sizes do not conform.')
   call assert(size(pint,2) == size(fluxes%flux_up,2), subname // ': sizes do not conform.')
   call assert(size(heating_rate,1) == size(fluxes%flux_up,1), subname // ': sizes do not conform.')
   call assert(size(heating_rate,2) == size(fluxes%flux_up,2)-1, subname // ': sizes do not conform.')

   ! Loop over levels and calculate heating rates; note that the fluxes *should*
   ! be defined at interfaces, so the loop ktop,kbot and grabbing the current
   ! and next value of k should be safe. ktop should be the top interface, and
   ! kbot + 1 should be the bottom interface.
   !
   ! NOTE: to get heating rate in K/day, normally we would use:
   !
   !     H = dF / dp * g * (sec/day) * (1e-5) / (cpair)
   !
   ! Here we just use
   !
   !     H = dF / dp * g
   !
   ! Why? Something to do with convenience with applying the fluxes to the
   ! heating tendency?
   do k = 1,size(pint,2)-1
      do i = 1,size(pint,1)
         heating_rate(i,k) = ( &
            fluxes%flux_up(i,k+1) - fluxes%flux_up(i,k) - &
            fluxes%flux_dn(i,k+1) + fluxes%flux_dn(i,k) &
         ) * gravit / (pint(i,k+1) - pint(i,k))
      end do
   end do

end subroutine calculate_heating_rate
    

subroutine set_aerosol_optics_lw(icall, state, pbuf, cam_layers, optics_out)
   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc, pbuf_get_index, &
                             pbuf_get_field
   use mo_optical_props, only: ty_optical_props_1scl
   use aer_rad_props, only: aer_rad_props_lw
   use radconstants, only: nlwbands

   integer, intent(in) :: icall
   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer, intent(in) :: cam_layers(:)
   type(ty_optical_props_1scl), intent(inout) :: optics_out

   ! Subroutine name for error messages
   character(len=*), parameter :: subroutine_name = 'set_aerosol_optics_lw'

   ! Output from CAM routines; expected to have dimension pcols
   real(r8) :: absorption_tau(pcols,pver,nlwbands)

   ! Loop variable
   integer :: k

   ! Get aerosol absorption optical depth from CAM routine
   absorption_tau = 0.0
   call aer_rad_props_lw(icall, state, pbuf, absorption_tau)

   ! Populate the RRTMGP optical properties object with CAM optical depth
   optics_out%tau(:,:,:) = 0.0
   do k = 1,size(cam_layers)
      if (cam_layers(k) > 0) then
         optics_out%tau(:,k,:) = absorption_tau(:,cam_layers(k),:)
      end if
   end do

end subroutine set_aerosol_optics_lw


subroutine set_aerosol_optics_sw(icall, state, pbuf, &
                                 day_indices, night_indices, cam_layers, &
                                 optics_out)
   use ppgrid, only: pcols, pver
   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc
   use aer_rad_props, only: aer_rad_props_sw
   use radconstants, only: nswbands
   use mo_optical_props, only: ty_optical_props_2str
   integer, intent(in) :: icall
   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer, intent(in) :: day_indices(:), night_indices(:)
   integer, intent(in) :: cam_layers(:)
   type(ty_optical_props_2str), intent(inout) :: optics_out

   ! NOTE: aer_rad_props expects 0:pver indexing on these! It appears this is to
   ! account for the extra layer added above model top, but it is not entirely
   ! clear. This is not done for the longwave, and it is not really documented
   ! anywhere that I can find. Regardless, optical properties for the zero index
   ! are set to zero in aer_rad_props_sw as far as I can tell.
   !
   ! NOTE: dimension ordering is different than for cloud optics!
   real(r8), dimension(pcols,0:pver,nswbands) :: tau, tau_w, tau_w_g, tau_w_f

   integer :: ncol, nbnd, nlay
   integer :: iday, icol, ilay
   integer :: k_rad, k_cam

   ! Everyone needs a name
   character(len=*), parameter :: subroutine_name = 'set_aerosol_optics_sw'

   ncol = state%ncol
   nbnd = nswbands

   ! Get aerosol absorption optical depth from CAM routine
   tau = 0.0
   tau_w = 0.0
   tau_w_g = 0.0
   tau_w_f = 0.0
   call aer_rad_props_sw(icall, state, pbuf, &
                         count(night_indices > 0), night_indices, &
                         tau, tau_w, tau_w_g, tau_w_f)

   ! Reset outputs (also handles case where radiation grid contains an extra
   ! layer above CAM grid)
   optics_out%tau = 0
   optics_out%ssa = 1
   optics_out%g = 0

   ! Assign daytime columns
   do iday = 1,count(day_indices > 0)

      ! Get index into full chunk-wide array for this daytime index
      icol = day_indices(iday)

      ! Loop over levels to map CAM grid to radiation grid (radiation grid might
      ! contain extra levels, or vice versa)
      do k_rad = 1,size(cam_layers)

         ! Get location of rad level on CAM grid. If this level is ABOVE the CAM
         ! grid, then skip this index (leave optical properties at initial
         ! values, 0 for tau and g, 1 for ssa).
         k_cam = cam_layers(k_rad)
         if (k_cam < 1) cycle

         ! Copy cloud optical depth over directly
         optics_out%tau(iday,k_rad,1:nbnd) = tau(icol,k_cam,1:nbnd)

         ! Extract single scattering albedo from the product-defined fields
         where (tau(icol,k_cam,1:nbnd) > 0)
            optics_out%ssa(iday,k_rad,1:nbnd) &
               = tau_w(icol,k_cam,1:nbnd) / tau(icol,k_cam,1:nbnd)
         elsewhere
            optics_out%ssa(iday,k_rad,1:nbnd) = 1
         endwhere

         ! Extract assymmetry parameter from the product-defined fields
         where (tau_w(icol,k_cam,1:nbnd) > 0)
            optics_out%g(iday,k_rad,1:nbnd) &
               = tau_w_g(icol,k_cam,1:nbnd) / tau_w(icol,k_cam,1:nbnd)
         elsewhere
            optics_out%g(iday,k_rad,1:nbnd) = 0
         endwhere

      end do

   end do

   ! We need to fix band ordering because the old input files assume RRTMG band
   ! ordering, but this has changed in RRTMGP.
   ! TODO: fix the input files themselves!
   do icol = 1,size(optics_out%tau,1)
      do ilay = 1,size(optics_out%tau,2)
         optics_out%tau(icol,ilay,:) = reordered(optics_out%tau(icol,ilay,:), map_rrtmg_to_rrtmgp_swbands)
         optics_out%ssa(icol,ilay,:) = reordered(optics_out%ssa(icol,ilay,:), map_rrtmg_to_rrtmgp_swbands)
         optics_out%g(icol,ilay,:) = reordered(optics_out%g(icol,ilay,:), map_rrtmg_to_rrtmgp_swbands)
      end do
   end do

   ! Check values
   call handle_error(optics_out%validate())
   
end subroutine set_aerosol_optics_sw


subroutine compress_day_columns_1d(xcol, xday, day_indices)

   real(r8), intent(in) :: xcol(:)
   real(r8), intent(inout) :: xday(:)
   integer, intent(in) :: day_indices(:)
   integer :: icol, iday
   character(len=32) :: subname = 'compress_day_columns'

   ! Check dimensions
   call assert(size(xday, 1) == size(day_indices, 1), trim(subname) // ': inconsistent sizes')
   call assert(size(xcol, 1) >= size(day_indices, 1), trim(subname) // ': inconsistent sizes')
   
   ! Do mapping
   do iday = 1,count(day_indices>0)
      icol = day_indices(iday)
      xday(iday) = xcol(icol)
   end do

end subroutine
!-------------------------------------------------------------------------------
subroutine compress_day_columns_2d(xcol, xday, day_indices)

   real(r8), intent(in) :: xcol(:,:)
   real(r8), intent(inout) :: xday(:,:)
   integer, intent(in) :: day_indices(:)
   integer :: icol, iday
   character(len=32) :: subname = 'compress_day_columns'

   ! Check dimensions
   call assert(size(xday, 1) == size(day_indices, 1), trim(subname) // ': inconsistent sizes')
   call assert(size(xcol, 1) >= size(day_indices, 1), trim(subname) // ': inconsistent sizes')
   call assert(size(xday, 2) == size(xcol, 2), trim(subname) // ': inconsistent sizes')

   ! Do mapping
   do iday = 1,count(day_indices>0)
      icol = day_indices(iday)
      xday(iday,:) = xcol(icol,:)
   end do

end subroutine
!-------------------------------------------------------------------------------


subroutine expand_daytime_fluxes(daytime_fluxes, expanded_fluxes, day_indices)
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

end subroutine expand_daytime_fluxes


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


subroutine output_cloud_optics_lw(state, optics)
   use physics_types, only: physics_state
   use cam_optics, only: cam_optics_type
   use cam_history, only: outfld

   type(physics_state), intent(in) :: state
   type(cam_optics_type), intent(in) :: optics
   
   ! Output 
   call assert_valid(optics%optical_depth(:state%ncol,:pver,:nlwbands), 'cloud_tau_lw')
   call outfld('CLOUD_TAU_LW', optics%optical_depth(:state%ncol,:pver,:nlwbands), state%ncol, state%lchnk)

end subroutine output_cloud_optics_lw


subroutine output_cloud_optics_sw(state, optics)
   use physics_types, only: physics_state
   use cam_optics, only: cam_optics_type
   use cam_history, only: outfld
   use radconstants, only: idx_sw_diag

   type(physics_state), intent(in) :: state
   type(cam_optics_type), intent(in) :: optics

   ! Output 
   call assert_valid(optics%optical_depth(:state%ncol,:pver,:nswbands), 'cloud_tau_sw')
   call outfld('CLOUD_TAU_SW', &
               optics%optical_depth(:state%ncol,:pver,:nswbands), &
               state%ncol, state%lchnk)
   call outfld('TOT_ICLD_VISTAU', &
               optics%optical_depth(:state%ncol,:pver,idx_sw_diag), &
               state%ncol, state%lchnk)
end subroutine output_cloud_optics_sw


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
         if (iday <= size(day_indices)) then
            day_indices(iday) = icol
         else
            call endrun(trim(subname) // ': iday > size(day_indices)')
         end if
      else
         inight = inight + 1
         if (inight <= size(night_indices)) then
            night_indices(inight) = icol
         else
            call endrun(trim(subname) // ': inight > size(night_indices)')
         end if
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
   call get_rlat_all_p(state%lchnk, state%ncol, clat(:state%ncol))
   call get_rlon_all_p(state%lchnk, state%ncol, clon(:state%ncol))

   ! Call zenith to calculate cosine solar zenith angle. Note we only pass the
   ! first state%ncol columns in case coszrs was allocated to be larger than
   ! ncol in the calling routine (i.e, pcols).
   call zenith(calday, clat(:state%ncol), clon(:state%ncol), &
               coszrs(:state%ncol), state%ncol, dt)
end subroutine set_cosine_solar_zenith_angle


! Set surface albedos from cam surface exchange object for direct and diffuse
! beam radiation. This code was copied from the RRTMG implementation, but moved
! to a subroutine with some better variable names.
subroutine set_albedo(cam_in, albedo_direct, albedo_diffuse)
   use camsrfexch, only: cam_in_t
   type(cam_in_t), intent(in) :: cam_in
   real(r8), intent(inout) :: albedo_direct(:,:)   ! surface albedo, direct radiation
   real(r8), intent(inout) :: albedo_diffuse(:,:)  ! surface albedo, diffuse radiation

   ! Local namespace
   real(r8) :: wavenumber_limits(2,nswbands)
   integer :: ncol, ibnd

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
   do ibnd = 1,nswbands
      if (is_visible(wavenumber_limits(1,ibnd)) .and. &
          is_visible(wavenumber_limits(2,ibnd))) then

         ! Entire band is in the visible
         albedo_direct(ibnd,:ncol) = cam_in%asdir(:ncol)
         albedo_diffuse(ibnd,:ncol) = cam_in%asdif(:ncol)

      else if (.not.is_visible(wavenumber_limits(1,ibnd)) .and. &
               .not.is_visible(wavenumber_limits(2,ibnd))) then

         ! Entire band is in the longwave
         albedo_direct(ibnd,:ncol) = cam_in%aldir(:ncol)
         albedo_diffuse(ibnd,:ncol) = cam_in%aldif(:ncol)

      else

         ! Band straddles the visible to near-infrared transition, so we take
         ! the albedo to be the average of the visible and near-infrared
         ! broadband albedos
         albedo_direct(ibnd,:ncol) = 0.5 * (cam_in%aldir(:ncol) + cam_in%asdir(:ncol))
         albedo_diffuse(ibnd,:ncol) = 0.5 * (cam_in%aldif(:ncol) + cam_in%asdif(:ncol))

      end if
   end do

#ifdef DO_FRAGILE_BAND_MAPPING
   ! Surface albedo (band mapping is hardcoded for RRTMG(P) code)
   ! This mapping assumes nswbands=14.
   ! TODO: we should handle this more intelligently
   if (nswbands /= 14) then
      call endrun(module_name //': ERROR: albedo band mapping assumes nswbands=14')
   end if

   ! Code lifted from RRTMG implementation
   do icol = 1,size(albedo_direct, 2)

      do ibnd = 1,9
         albedo_direct(ibnd,icol) = cam_in%aldir(icol)
         albedo_diffuse(ibnd,icol) = cam_in%aldif(icol)
      enddo

      ! Set band 25 (or, band 10 counting from 1) to use linear average of UV/visible
      ! and near-IR values, since this band straddles 0.7 microns:
      albedo_direct(10,icol) = 0.5*(cam_in%aldir(icol) + cam_in%asdir(icol))
      albedo_diffuse(10,icol) = 0.5*(cam_in%aldif(icol) + cam_in%asdif(icol))

      ! UV/visible bands 25-28 (11-14), 16000-50000 cm-1, 0.200-0.625 micron
      do ibnd = 11,nswbands
         albedo_direct(ibnd,icol) = cam_in%asdir(icol)
         albedo_diffuse(ibnd,icol) = cam_in%asdif(icol)
      enddo
   end do
#endif

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

subroutine clip_values_1d(x, min_x, max_x, varname, warn)
   real(r8), intent(inout) :: x(:)
   real(r8), intent(in) :: min_x
   real(r8), intent(in) :: max_x
   character(len=*), intent(in), optional :: varname
   logical, intent(in), optional :: warn

   logical :: warn_local

   warn_local = .false.
   if (present(warn)) then
      warn_local = warn
   end if

   ! Look for values less than threshold
   if (any(x < min_x)) then
      ! Raise warning?
      if (warn_local) then
         if (present(varname)) then
            print *, module_name // ' warning: ', &
                     count(x < min_x), ' values are below threshold for variable ', &
                     trim(varname), '; min = ', minval(x)
         else
            print *, module_name // ' warning: ', &
                     count(x < min_x), ' values are below threshold; min = ', minval(x)
         end if
      end if

      ! Clip values
      where (x < min_x)
         x = min_x
      endwhere
   end if

   ! Look for values greater than threshold 
   if (any(x > max_x)) then 
      ! Raise warning?
      if (warn_local) then
         if (present(varname)) then
            print *, module_name // ' warning: ', &
                     count(x > max_x), ' values are above threshold for variable ', &
                     trim(varname), '; max = ', maxval(x)
         else
            print *, module_name // ' warning: ', &
                     count(x > max_x), ' values are above threshold; max = ', maxval(x)
         end if
      end if

      ! Clip values
      where (x > max_x)
         x = max_x
      end where
   end if
end subroutine clip_values_1d

subroutine clip_values_2d(x, min_x, max_x, varname, warn)
   real(r8), intent(inout) :: x(:,:)
   real(r8), intent(in) :: min_x
   real(r8), intent(in) :: max_x
   character(len=*), intent(in), optional :: varname
   logical, intent(in), optional :: warn

   integer :: i, j, nbad
   logical :: warn_local

   warn_local = .false.
   if (present(warn)) then
      warn_local = warn
   end if

   ! look for values less than threshold
   if (any(x < min_x)) then
      ! Raise warning?
      if (warn_local) then
         if (present(varname)) then
            print *, module_name // ' warning: ', &
                     count(x < min_x), ' values are below threshold for variable ', &
                     trim(varname), '; min = ', minval(x)
         else
            print *, module_name // ' warning: ', &
                     count(x < min_x), ' values are below threshold; min = ', minval(x)
         end if
      end if

      ! Clip values
      where (x < min_x)
         x = min_x
      endwhere
   end if

   ! Look for values greater than threshold
   if (any(x > max_x)) then 
      ! Raise warning?
      if (warn_local) then
         if (present(varname)) then
            print *, module_name // ' warning: ', &
                     count(x > max_x), ' values are above threshold for variable ', &
                     trim(varname), '; max = ', maxval(x)
         else
            print *, module_name // ' warning: ', &
                     count(x > max_x), ' values are above threshold; max = ', maxval(x)
         end if
      end if

      ! Clip values
      where (x > max_x)
         x = max_x
      end where
   end if

end subroutine clip_values_2d


! Send shortwave fluxes and heating rates to history buffer
subroutine output_fluxes_sw(icall, state, flux_all, flux_clr, lev_indices)
   use physconst, only: cpair
   use physics_types, only: physics_state
   use cam_history, only: outfld
   use mo_fluxes_byband, only: ty_fluxes_byband
   
   integer, intent(in) :: icall
   type(physics_state), intent(in) :: state
   type(ty_fluxes_byband), intent(in) :: flux_all
   type(ty_fluxes_byband), intent(in) :: flux_clr
   integer, intent(in) :: lev_indices(:)

   ! Working variables
   integer :: ncol, nlev
   integer :: ktop_cam, kbot_cam
   integer :: ktop_rad, kbot_rad
   integer :: k

   ncol = state%ncol
   nlev = size(flux_all%flux_up,2)
   ktop_rad = 1
   kbot_rad = nlev

   ! Find index to top of model
   ktop_cam = 0
   do k = 1,size(lev_indices)
      if (lev_indices(k) == 1) then
         ktop_cam = k
         exit
      end if
   end do
   call assert(ktop_cam > 0, ': no layer top found.')
   kbot_cam = kbot_rad

   ! All-sky flux_all%fluxes at model interfaces
   call outfld('FDS'//diag(icall), flux_all%flux_dn(:ncol,ktop_cam:), ncol, state%lchnk)
   call outfld('FUS'//diag(icall), flux_all%flux_up(:ncol,ktop_cam:), ncol, state%lchnk)
   call outfld('FNS'//diag(icall), flux_all%flux_net(:ncol,ktop_cam:), ncol, state%lchnk)
   call outfld('FDS_DIR'//diag(icall), flux_all%flux_dn_dir(:ncol,ktop_cam:), ncol, state%lchnk)

   ! Clear-sky flux_all%fluxes at model interfaces
   call outfld('FDSC'//diag(icall), flux_clr%flux_dn(:ncol,ktop_cam:), ncol, state%lchnk)
   call outfld('FUSC'//diag(icall), flux_clr%flux_up(:ncol,ktop_cam:), ncol, state%lchnk)
   call outfld('FNSC'//diag(icall), flux_clr%flux_net(:ncol,ktop_cam:), ncol, state%lchnk)
   call outfld('FDSC_DIR'//diag(icall), flux_clr%flux_dn_dir(:ncol,ktop_cam:), ncol, state%lchnk)

   ! All-sky flux_all%fluxes
   call outfld('FSNT'//diag(icall), flux_all%flux_net(:ncol,ktop_cam), ncol, state%lchnk)
   call outfld('FSNS'//diag(icall), flux_all%flux_net(:ncol,kbot_cam), ncol, state%lchnk)
   call outfld('FSDS'//diag(icall), flux_all%flux_dn(:ncol,kbot_cam), ncol, state%lchnk)

   ! TOA flux_all%fluxes (above model top, use index to rad top)
   call outfld('FSUTOA'//diag(icall), flux_all%flux_up(:ncol,ktop_rad), ncol, state%lchnk)
   call outfld('FSNTOA'//diag(icall), flux_all%flux_net(:ncol,ktop_rad), ncol, state%lchnk)

   ! Clear-sky flux_all%fluxes
   call outfld('FSNTC'//diag(icall), flux_clr%flux_net(:ncol,ktop_cam), ncol, state%lchnk)
   call outfld('FSNSC'//diag(icall), flux_clr%flux_net(:ncol,kbot_cam), ncol, state%lchnk)
   call outfld('FSDSC'//diag(icall), flux_clr%flux_dn(:ncol,kbot_cam), ncol, state%lchnk)
   call outfld('FSUTOAC'//diag(icall), flux_clr%flux_up(:ncol,ktop_rad), ncol, state%lchnk)
   call outfld('FSNTOAC'//diag(icall), flux_clr%flux_net(:ncol,ktop_rad), ncol, state%lchnk)

end subroutine output_fluxes_sw


! Send longwave fluxes and heating rates to history buffer
subroutine output_fluxes_lw(icall, state, flux_all, flux_clr, lev_indices)
   use physconst, only: cpair
   use physics_types, only: physics_state
   use cam_history, only: outfld
   use mo_fluxes_byband, only: ty_fluxes_byband
   
   integer, intent(in) :: icall
   type(physics_state), intent(in) :: state
   type(ty_fluxes_byband), intent(in) :: flux_all
   type(ty_fluxes_byband), intent(in) :: flux_clr
   integer, intent(in) :: lev_indices(:)

   ! Working arrays
   real(r8), dimension(pcols,pver+1) :: flux_up, flux_dn, flux_net

   integer :: ncol, nlev
   integer :: ktop_cam, kbot_cam
   integer :: ktop_rad, kbot_rad
   integer :: k

   ncol = state%ncol
   nlev = size(flux_all%flux_up, 2)

   ! Do error checking here

   ktop_rad = 1
   kbot_rad = nlev

   ! Find index to top of model
   do k = 1,size(lev_indices)
      if (lev_indices(k) == 1) then
         ktop_cam = k
         exit
      end if
   end do
   kbot_cam = nlev

   ! Do all-sky fluxes, map first
   call map_rad_to_cam_levels(flux_all%flux_up(:ncol,:nlev), flux_up(:ncol,:pver+1), lev_indices(:nlev))
   call map_rad_to_cam_levels(flux_all%flux_dn(:ncol,:nlev), flux_dn(:ncol,:pver+1), lev_indices(:nlev))
   flux_net(:ncol,:pver+1) = flux_up(:ncol,:pver+1) - flux_dn(:ncol,:pver+1)

   ! All-sky fluxes at model interfaces
   call outfld('FDL'//diag(icall), flux_dn(:ncol,:pver+1), ncol, state%lchnk)
   call outfld('FUL'//diag(icall), flux_up(:ncol,:pver+1), ncol, state%lchnk)
   call outfld('FNL'//diag(icall), flux_net(:ncol,:pver+1), ncol, state%lchnk)

   ! Clear-sky fluxes at model interfaces
   call map_rad_to_cam_levels(flux_clr%flux_up(:ncol,:nlev), flux_up(:ncol,:pver+1), lev_indices(:nlev))
   call map_rad_to_cam_levels(flux_clr%flux_dn(:ncol,:nlev), flux_dn(:ncol,:pver+1), lev_indices(:nlev))
   flux_net(:ncol,:pver+1) = flux_up(:ncol,:pver+1) - flux_dn(:ncol,:pver+1)
   call outfld('FDLC'//diag(icall), flux_dn(:ncol,:pver+1), ncol, state%lchnk)
   call outfld('FULC'//diag(icall), flux_up(:ncol,:pver+1), ncol, state%lchnk)
   call outfld('FNLC'//diag(icall), flux_net(:ncol,:pver+1), ncol, state%lchnk)

   ! All-sky fluxes
   ! NOTE: sign change on net fluxes because internally net fluxes are assumed
   ! to down minus up, but for longwave outputs we want upward positive
   call outfld('FLNT'//diag(icall), -flux_all%flux_net(:ncol,ktop_rad), ncol, state%lchnk)
   call outfld('FLNS'//diag(icall), -flux_all%flux_net(:ncol,kbot_rad), ncol, state%lchnk)
   call outfld('FLUT'//diag(icall), flux_all%flux_up(:ncol,ktop_rad), ncol, state%lchnk)
   call outfld('FLDS'//diag(icall), flux_all%flux_dn(:ncol,kbot_rad), ncol, state%lchnk)

   ! Clear-sky fluxes
   ! NOTE: sign change on net fluxes because internally net fluxes are assumed
   ! to down minus up, but for longwave outputs we want upward positive
   call outfld('FLNTC'//diag(icall), -flux_clr%flux_net(:ncol,ktop_rad), ncol, state%lchnk)
   call outfld('FLNSC'//diag(icall), -flux_clr%flux_net(:ncol,kbot_rad), ncol, state%lchnk)
   call outfld('FLUTC'//diag(icall), flux_clr%flux_up(:ncol,ktop_rad), ncol, state%lchnk)
   call outfld('FLDSC'//diag(icall), flux_clr%flux_dn(:ncol,kbot_rad), ncol, state%lchnk)

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


subroutine set_gas_concentrations(icall, state, pbuf, nlayers, &
                                  gas_concentrations, &
                                  lev_indices, &
                                  day_indices)

   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc
   use rad_constituents, only: rad_cnst_get_gas
   use mo_gas_concentrations, only: ty_gas_concs
   use mo_util_string, only: lower_case

   integer, intent(in) :: icall
   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer, intent(in) :: nlayers
   type(ty_gas_concs), intent(inout) :: gas_concentrations
   integer, intent(in) :: lev_indices(:)
   integer, intent(in), optional :: day_indices(:)

   ! Local variables
   real(r8), dimension(pcols,nlayers) :: vol_mix_ratio, &
                                         vol_mix_ratio_day, &
                                         vol_mix_ratio_out
   real(r8), pointer :: mass_mix_ratio(:,:)
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
   integer :: igas, iday, icol

   ! Number of columns and daytime columns
   integer :: ncol, nday

   character(len=32) :: subname = 'set_gas_concentrations'

   ! Indices to model top on the RADIATION VERTICAL GRID
   integer :: k, k_cam
   
   ! Check sizes
   call assert(size(lev_indices) == nlayers, trim(subname) // ': inconsistent sizes')

   ! Number of columns in chunk
   ncol = state%ncol

   ! For each gas species needed for RRTMGP, read the mass mixing ratio from the
   ! CAM rad_constituents interface, convert to volume mixing ratios, and
   ! subset for daytime-only indices if needed.
   do igas = 1,size(gas_species)

      ! If this gas is not in list of active gases, then skip
      if (.not.string_in_list(gas_species(igas), active_gases)) cycle

      ! initialize
      vol_mix_ratio = 0._r8

      select case(trim(gas_species(igas)))

         case('CO')

            ! CO not available, use default
            vol_mix_ratio(:,:) = co_vol_mix_ratio

         case('N2')

            ! N2 not available, use default
            vol_mix_ratio(:,:) = n2_vol_mix_ratio

         case('H2O')

            ! Water vapor is represented as specific humidity in CAM, so we
            ! need to handle water a little differently
            call rad_cnst_get_gas(icall, trim(gas_species(igas)), state, pbuf, &
                                  mass_mix_ratio)
            
            ! Convert to volume mixing ratio by multiplying by the ratio of
            ! molecular weight of dry air to molecular weight of gas. Note that
            ! first specific humidity (held in the mass_mix_ratio array read
            ! from rad_constituents) is converted to an actual mass mixing
            ! ratio.
            vol_mix_ratio(:ncol,:pver) = mass_mix_ratio(:ncol,:pver) / ( &
               1._r8 - mass_mix_ratio(:ncol,:pver) &
            )  * mol_weight_air / mol_weight_gas(igas)

         case DEFAULT

            ! Get mass mixing ratio from the rad_constituents interface
            call rad_cnst_get_gas(icall, trim(gas_species(igas)), state, pbuf, &
                                  mass_mix_ratio)

            ! Convert to volume mixing ratio by multiplying by the ratio of
            ! molecular weight of dry air to molecular weight of gas
            vol_mix_ratio(:ncol,:pver) = mass_mix_ratio(:ncol,:pver) &
                                       * mol_weight_air / mol_weight_gas(igas)

      end select

      ! If this gas is not in list of active gases, then zero out (useful for
      ! debugging)
      if (.not.string_in_list(gas_species(igas), active_gases)) then
         vol_mix_ratio(:,:) = 0._r8
      end if

      ! Make sure we do not have any negative volume mixing ratios
      call assert(all(vol_mix_ratio(:ncol,:pver) >= 0), &
                  trim(subname) // ': invalid gas concentration for ' // &
                  trim(gas_species(igas)))

      ! Map to radiation grid
      call map_cam_to_rad_levels(vol_mix_ratio(:ncol,:pver), &
                                 vol_mix_ratio_out(:ncol,:nlayers), &
                                 lev_indices(:nlayers))

      ! Populate the RRTMGP gas concentration object with values for this gas.
      ! NOTE: RRTMGP makes some assumptions about gas names internally, so we
      ! need to pass the gas names as LOWER case here!
      if (present(day_indices)) then

         ! Number of daytime columns in this chunk is number of indices in
         ! day_indices array that are greater than zero
         nday = count(day_indices > 0)

         ! Populate compressed array with just daytime values
         call compress_day_columns(vol_mix_ratio_out(:ncol,:nlayers), &
                                   vol_mix_ratio_day(:nday,:nlayers), &
                                   day_indices(:nday))
         
         ! Set volumn mixing ratio in gas concentration object for just daytime
         ! columns in this chunk
         call handle_error( &
            gas_concentrations%set_vmr(trim(lower_case(gas_species(igas))),vol_mix_ratio_day(:nday,:nlayers)) & 
         )
      else
         ! Set volumn mixing ratio in gas concentration object for just columns
         ! in this chunk
         call handle_error( &
            gas_concentrations%set_vmr(trim(lower_case(gas_species(igas))),vol_mix_ratio_out(:ncol,:nlayers)) & 
         )
      end if

   end do

end subroutine set_gas_concentrations


function get_cam_layers(nlay_rad, nlay_cam) result(cam_layers)

   integer, intent(in) :: nlay_rad, nlay_cam
   integer :: cam_layers(nlay_rad)
   integer :: k

   ! CAM layers that correspond to RRTMGP layers; note that an extra layer above
   ! the model top may be used for the RRTMGP grid.
   cam_layers(:) = 0
   do k = 1,nlay_rad
      cam_layers(k) = nlay_cam-nlay_rad+k
   end do

end function get_cam_layers


logical function string_in_list(string, list)
   character(len=*), intent(in) :: string
   character(len=*), intent(in) :: list(:)
   integer :: i
   
   ! Set to false initially
   string_in_list = .false.

   ! Loop over items in list, and if we find a match then set flag to true and
   ! exit loop
   do i = 1,size(list)
      if (trim(list(i)) == trim(string)) then
         string_in_list = .true.
         exit 
      end if
   end do
end function string_in_list


subroutine handle_error(error_message, stop_on_error)
   character(len=*), intent(in) :: error_message
   logical, intent(in), optional :: stop_on_error
   logical :: stop_on_error_local = .true.

   ! Allow passing of an optional flag to not stop the run if an error is
   ! encountered. This allows this subroutine to be used when inquiring if a
   ! variable exists without failing.
   if (present(stop_on_error)) then
      stop_on_error_local = stop_on_error
   else
      stop_on_error_local = .true.
   end if

   ! If we encounter an error, fail if we require success. Otherwise do
   ! nothing and return silently.
   if (len(trim(error_message)) > 0) then
      if (stop_on_error_local) then
         call endrun(module_name // ': ' // error_message)
      end if
   end if
end subroutine handle_error


subroutine set_cloud_optics_sw(state, pbuf, col_indices, lev_indices, kdist, optics_out)
   
   use ppgrid, only: pcols, pver, pverp
   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc, &
                             pbuf_get_field, &
                             pbuf_get_index
   use mo_optical_props, only: ty_optical_props_2str
   use mo_gas_optics, only: ty_gas_optics
   use mcica_subcol_gen, only: mcica_subcol_mask
   use cam_optics, only: cam_optics_type, get_cloud_optics_sw

   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer, intent(in) :: col_indices(:), lev_indices(:)
   type(ty_gas_optics), intent(in) :: kdist
   type(ty_optical_props_2str), intent(inout) :: optics_out

   ! Type to hold optics on CAM grid
   type(cam_optics_type) :: optics_cam

   ! Pointer to cloud fraction on physics buffer
   real(r8), pointer :: cloud_fraction(:,:), snow_fraction(:,:)

   ! Combined cloud and snow fraction
   real(r8) :: combined_cloud_fraction(pcols,pver)

   ! For MCICA sampling routine
   integer, parameter :: changeseed = 1

   ! Dimension sizes
   integer :: ngpt, ncol, nday, nlay

   ! McICA subcolumn cloud flag
   logical, dimension(nswgpts,pcols,pver) :: iscloudy

   ! Loop variables
   integer :: icol, ilev, igpt, ibnd, icol_cam, ilev_cam

   character(len=32) :: subname = 'set_cloud_optics_sw'

   ncol = state%ncol
   ngpt = kdist%get_ngpt()

   ! Retrieve the mean in-cloud optical properties via CAM cloud radiative
   ! properties interface (cloud_rad_props). This retrieves cloud optical
   ! properties by *band* -- these will be mapped to g-points when doing
   ! the subcolumn sampling to account for cloud overlap.
   call optics_cam%initialize(nswbands, ncol, pver)
   call get_cloud_optics_sw(state, pbuf, optics_cam)

   ! We need to fix band ordering because the old input files assume RRTMG band
   ! ordering, but this has changed in RRTMGP.
   ! TODO: fix the input files themselves!
   do icol = 1,size(optics_cam%optical_depth,1)
      do ilev = 1,size(optics_cam%optical_depth,2)
         optics_cam%optical_depth(icol,ilev,:) = reordered(optics_cam%optical_depth(icol,ilev,:), map_rrtmg_to_rrtmgp_swbands)
         optics_cam%single_scattering_albedo(icol,ilev,:) = reordered(optics_cam%single_scattering_albedo(icol,ilev,:), map_rrtmg_to_rrtmgp_swbands)
         optics_cam%assymmetry_parameter(icol,ilev,:) = reordered(optics_cam%assymmetry_parameter(icol,ilev,:), map_rrtmg_to_rrtmgp_swbands)
      end do
   end do

   ! Send in-cloud optical depth for visible band to history buffer
   ! TODO: send single scattering albedo and assymmetry parameter for
   ! diagnostic purposes as well?
   call output_cloud_optics_sw(state, optics_cam)

   ! Initialize (or reset) output cloud optics object
   optics_out%tau = 0.0
   optics_out%ssa = 1.0
   optics_out%g = 0.0

   ! Get cloud and snow fractions, and combine
   call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction)
   call pbuf_get_field(pbuf, pbuf_get_index('CLDFSNOW'), snow_fraction)
   combined_cloud_fraction(1:ncol,1:pver) = max(cloud_fraction(1:ncol,1:pver), &
                                                snow_fraction(1:ncol,1:pver))

   ! Do MCICA sampling of optics here. This will map bands to gpoints,
   ! while doing stochastic sampling of cloud state
   call mcica_subcol_mask(ngpt, ncol, pver, changeseed, &
                          state%pmid(1:ncol,1:pver), &
                          combined_cloud_fraction(1:ncol,1:pver), &
                          iscloudy(1:ngpt,1:ncol,1:pver))
   
   ! -- generate subcolumns for homogeneous clouds -----
   ! where there is a cloud, set the subcolumn cloud properties;
   optics_out%tau(:,:,:) = 0
   optics_out%ssa(:,:,:) = 1
   optics_out%g(:,:,:) = 0
   do ilev = 1,size(lev_indices)

      ! Map vertical levels; if this level is above model top, skip it (leave set
      ! to initial values of tau = 0, ssa = 1, g = 0).
      ilev_cam = lev_indices(ilev)
      if (ilev_cam < 1) cycle

      ! Loop over columns and map CAM columns to those on radiation grid
      ! (daytime-only columns)
      do icol = 1,size(col_indices)
         icol_cam = col_indices(icol)
         call assert(icol_cam > 0, trim(subname) // ': col_indices < 0')

         ! Loop over g-points and map bands to g-points; each subcolumn
         ! corresponds to a single g-point. This is how this code implements the
         ! McICA assumptions: simultaneously sampling over cloud state and
         ! g-point.
         do igpt = 1,ngpt
            if (iscloudy(igpt,icol_cam,ilev_cam)) then
               ibnd = kdist%convert_gpt2band(igpt)
               optics_out%tau(icol,ilev,igpt) = optics_cam%optical_depth(icol_cam,ilev_cam,ibnd)
               optics_out%ssa(icol,ilev,igpt) = optics_cam%single_scattering_albedo(icol_cam,ilev_cam,ibnd)
               optics_out%g(icol,ilev,igpt) = optics_cam%assymmetry_parameter(icol_cam,ilev_cam,ibnd)
            else
               optics_out%tau(icol,ilev,igpt) = 0._r8
               optics_out%ssa(icol,ilev,igpt) = 1._r8
               optics_out%g(icol,ilev,igpt) = 0._r8
            end if
         end do
      end do
   end do

   ! Apply delta scaling to account for forward-scattering
   ! TODO: delta_scale takes the forward scattering fraction as an optional
   ! parameter. In the current cloud optics_sw scheme, forward scattering is taken
   ! just as g^2, which delta_scale assumes if forward scattering fraction is
   ! omitted in the function call. In the future, we should explicitly pass
   ! this. This just requires modifying the get_cloud_optics_sw procedures to also
   ! pass the foward scattering fraction that the CAM cloud optics_sw assumes.
   !call handle_error(optics_out%delta_scale())

   ! Check cloud optics_sw
   call handle_error(optics_out%validate())

end subroutine set_cloud_optics_sw


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


subroutine set_cloud_optics_lw(state, pbuf, lev_indices, kdist, optics_out)
   
   use ppgrid, only: pcols, pver
   use physics_types, only: physics_state
   use physics_buffer, only: physics_buffer_desc, &
                             pbuf_get_field, pbuf_get_index
   use mo_optical_props, only: ty_optical_props_1scl
   use mo_gas_optics, only: ty_gas_optics
   use mcica_subcol_gen, only: mcica_subcol_mask
   use cam_optics, only: cam_optics_type, get_cloud_optics_lw

   type(physics_state), intent(in) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   integer, intent(in) :: lev_indices(:)
   type(ty_gas_optics), intent(in) :: kdist
   type(ty_optical_props_1scl), intent(inout) :: optics_out

   type(cam_optics_type) :: optics_cam
   real(r8), pointer :: cloud_fraction(:,:)
   real(r8), pointer :: snow_fraction(:,:)
   real(r8) :: combined_cloud_fraction(pcols,pver)

   ! For MCICA sampling routine
   integer, parameter :: changeseed = 1

   ! Dimension sizes
   integer :: ngpt, ncol

   ! Temporary arrays to hold mcica-sampled cloud optics (ngpt,ncol,pver)
   real(r8), dimension(pcols,pver,nlwgpts) :: optical_depth_gpt
   logical, dimension(nlwgpts,pcols,pver) :: iscloudy

   ! Loop variables
   integer :: icol, ilev, igpt, ibnd, ilev_cam

   ! Initialize (or reset) output cloud optics object
   optics_out%tau = 0.0

   ! Set dimension size working variables
   ngpt = kdist%get_ngpt()
   ncol = state%ncol

   ! Get cloud optics using CAM routines. This should combine cloud with snow
   ! optics, if "snow clouds" are being considered
   call optics_cam%initialize(nlwbands, ncol, pver)
   call get_cloud_optics_lw(state, pbuf, optics_cam)

   ! Check values
   call assert_range(optics_cam%optical_depth, 0._r8, 1e20_r8, &
                     'set_cloud_optics_lw: optics_cam%optical_depth')

   ! Send cloud optics to history buffer
   call output_cloud_optics_lw(state, optics_cam)

   ! Get cloud and snow fractions, and combine
   call pbuf_get_field(pbuf, pbuf_get_index('CLD'), cloud_fraction)
   call pbuf_get_field(pbuf, pbuf_get_index('CLDFSNOW'), snow_fraction)

   ! Combine cloud and snow fractions for MCICA sampling
   combined_cloud_fraction(1:ncol,1:pver) = max(cloud_fraction(1:ncol,1:pver), &
                                                snow_fraction(1:ncol,1:pver))

   ! Do MCICA sampling of optics here. This will map bands to gpoints,
   ! while doing stochastic sampling of cloud state
   !
   ! First, just get the stochastic subcolumn cloudy mask...
   call mcica_subcol_mask(ngpt, ncol, pver, changeseed, &
                          state%pmid(1:ncol,1:pver), &
                          combined_cloud_fraction(1:ncol,1:pver), &
                          iscloudy(1:ngpt,1:ncol,1:pver))

   ! ... and now map optics to g-points, selecting a single subcolumn for each
   ! g-point. This implementation generates homogeneous clouds, but it would be
   ! straightforward to extend this to handle horizontally heterogeneous clouds
   ! as well.
   ! NOTE: incoming optics should be in-cloud quantites and not grid-averaged 
   ! quantities!
   optics_out%tau = 0
   do ilev = 1,size(lev_indices)

      ! Get level index on CAM grid (i.e., the index that this rad level
      ! corresponds to in CAM fields). If this index is above the model top
      ! (index less than 0) then skip setting the optical properties for this
      ! level (leave set to zero)
      ilev_cam = lev_indices(ilev)
      if (ilev_cam < 1) cycle

      do icol = 1,ncol
         do igpt = 1,ngpt
            if (iscloudy(igpt,icol,ilev_cam) .and. (combined_cloud_fraction(icol,ilev_cam) > 0._r8) ) then
               ibnd = kdist%convert_gpt2band(igpt)
               optics_out%tau(icol,ilev,igpt) = optics_cam%optical_depth(icol,ilev_cam,ibnd)
            else
               optics_out%tau(icol,ilev,igpt) = 0._r8
            end if
         end do
      end do
   end do

   ! Apply delta scaling to account for forward-scattering
   ! TODO: delta_scale takes the forward scattering fraction as an optional
   ! parameter. In the current cloud optics_lw scheme, forward scattering is taken
   ! just as g^2, which delta_scale assumes if forward scattering fraction is
   ! omitted in the function call. In the future, we should explicitly pass
   ! this. This just requires modifying the get_cloud_optics_lw procedures to also
   ! pass the foward scattering fraction that the CAM cloud optics_lw assumes.
   !call handle_error(optics_out%delta_scale())

   ! Check values
   call assert_range(optics_out%tau, 0._r8, 1e20_r8, &
                     'set_cloud_optics_lw: optics_out%tau')

   ! Check cloud optics
   call handle_error(optics_out%validate())

   call optics_cam%finalize()

end subroutine set_cloud_optics_lw


! Get index to highest vertical level above a minimum pressure. This function is
! provided to determine the highest CAM level to use in the radiation
! calculation for a given chunk, such that all levels are above some minimum
! pressure and thus within range of the RRTMGP lookup tables.
integer function min_pressure_index(pressure, min_pressure)
   real(r8), intent(in) :: pressure(:,:)
   real(r8), intent(in) :: min_pressure
   integer :: k

   do k = 1,size(pressure,2)
      if (all(pressure(:,k) > min_pressure)) then
         min_pressure_index = k
         exit
      end if
   end do
end function

! Get index to lowest vertical level with a pressure below a max pressure value.
! This function is provided to determine the lowest CAM level to use in the
! radiation calculation for a given chunk, such that all levels are below some
! maximum pressure and thus within range of the RRTMGP lookup tables.
integer function max_pressure_index(pressure, max_pressure)
   real(r8), intent(in) :: pressure(:,:)
   real(r8), intent(in) :: max_pressure
   integer :: k

   do k = size(pressure,2),1,-1
      if (all(pressure(:,k) < max_pressure)) then
         max_pressure_index = k
         exit
      end if
   end do
end function


end module radiation
