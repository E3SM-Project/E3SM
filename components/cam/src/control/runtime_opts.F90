module runtime_opts

!----------------------------------------------------------------------- 
! 
! Purpose: This module is responsible for reading CAM namelist cam_inparm 
!          and broadcasting namelist values if needed.  
! 
! Author:
!   Original routines:  CMS
!   Module:             T. Henderson, September 2003
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
use shr_kind_mod,    only: r8 => shr_kind_r8, SHR_KIND_CL
use spmd_utils,      only: masterproc
use namelist_utils,  only: find_group_name
use pmgrid,          only: plon
use cam_instance,    only: inst_suffix
use cam_history
use cam_control_mod
use cam_logfile,     only: iulog
use pspect
use units
use constituents,    only: pcnst, readtrace
use tracers,         only: tracers_flag
use time_manager,    only: dtime
use filenames,       only: ncdata, bnd_topo, &
                           absems_data, &
                           caseid, &
                           brnch_retain_casename
use dycore,          only: dycore_is
use cam_abortutils,      only: endrun
use rayleigh_friction, only: rayk0, raykrange, raytau0

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
implicit none
private
save


!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
public read_namelist        ! Set and/or get all runtime options

!-----------------------------------------------------------------------
! Private data ---------------------------------------------------------
!-----------------------------------------------------------------------

character(len=SHR_KIND_CL), private :: nlfilename = 'atm_in' ! Namelist filename

!-----------------------------------------------------------------------
!
! SOMEWHAT ALPHABETICAL listing of variables in the cam_inparm namelist:
!
! variable                description
! --------             -----------------
!
! bnd_topo             Path and filename of topography dataset
! 
! absems_data          Dataset with absorption and emissivity factors.
!
! dtime = nnnn,        Model time step in seconds. Default is dycore dependent.
! 
! nlvdry = nn,         Number of layers over which to do dry
!                      adjustment. Defaults to 3.
! 
! cam_branch_file      Filepath of restart file to branch from (nsrest=3)
!                      Full pathname required.
character(len=256) :: cam_branch_file = ' '
!

!------------------------------------------------------------------
! The following 3 are specific to Rayleigh friction
! integer rayk0         vertical level at which rayleigh friction term is centered
! real(r8) raykrange    range of rayleigh friction profile; if 0, range is set automatically
! real(r8) raytau0      approximate value of decay time at model top (days);
!                       if 0., no rayleigh friction is applied
!------------------------------------------------------------------
!
!
! pertlim = n.n        Max size of perturbation to apply to initial
!                      temperature field.
!
! new_random           logical: if .true., use RNG in dynamics/se/random_xgc.F90
!                      instead of the fortran intrinsic.
!
! seed_custom          integer: if > 0, use new seeding mechanism that uses a
!                      custom seed rather than a custom limit. Default 0
!
! seed_clock           logical: if .true., XOR the system_clock with the seed,
!                      wheter it includes a custom seed or not. Default .false.
!
! phys_alltoall        Dynamics/physics transpose option. See phys_grid module.
!
integer :: phys_alltoall
! 
! phys_loadbalance     Load balance option for performance tuning of 
!                      physics chunks.  See phys_grid module.  
integer :: phys_loadbalance
! 
! phys_twin_algorithm  Load balance option for performance tuning of 
!                      physics chunks.  See phys_grid module.  
integer :: phys_twin_algorithm
! 
! phys_chnk_per_thd    Performance tuning option for physics chunks.  See 
!                      phys_grid module.  
integer :: phys_chnk_per_thd
! 
! tracers_flag = .F.    If true, implement tracer test code. Number of tracers determined
!                      in tracers_suite.F90 must agree with PCNST
!
! readtrace = .T.      If true, tracer initial conditions obtained from 
!                      initial file. 
!
! print_step_cost      true => print per timestep cost info
!
!
!   logical indirect     
!                    ! true => include indirect radiative effects of
!                    ! sulfate aerosols.  Default is false.
!
!
! met_data_file        name of file that contains the offline meteorology data
! met_data_path        name of directory that contains the offline meteorology data
!
! met_filenames_list   name of file that contains names of the offline 
!                      meteorology data files
!
! met_remove_file      true => the offline meteorology file will be removed
!
! met_cell_wall_winds  true => the offline meteorology winds are defined on the model
!                      grid cell walls
! Physics buffer
logical :: pbuf_global_allocate       ! allocate all buffers as global (default: .true.)


! Conservation checks

logical            :: print_energy_errors ! switch for diagnostic output from check_energy module

! SCM Options
logical  :: single_column
real(r8) :: scmlat,scmlon
real(r8) :: scm_relaxation_low
real(r8) :: scm_relaxation_high
integer, parameter :: max_chars = 128
character(len=max_chars) iopfile
character(len=200) :: scm_clubb_iop_name
logical  :: scm_iop_srf_prop
logical  :: scm_relaxation
logical  :: scm_diurnal_avg
logical  :: scm_crm_mode
logical  :: scm_observed_aero
logical  :: swrad_off
logical  :: lwrad_off
logical  :: precip_off

contains

!=======================================================================

  subroutine read_namelist(single_column_in, scmlon_in, scmlat_in, nlfilename_in )

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Read data from namelist cam_inparm to define the run. Process some of the
   ! namelist variables to determine history and restart/branch file path 
   ! names.  Check input namelist variables for validity and print them
   ! to standard output. 
   ! 
   ! Method: 
   ! Important Note for running on SUN systems: "implicit automatic (a-z)"
   ! will not work because namelist data must be static.
   !
   ! Author: 
   ! Original version:  CCM1
   ! Standardized:      L. Bath, June 1992
   !                    T. Acker, March 1996
   !     
   !-----------------------------------------------------------------------

   ! Note that the following interfaces are prototypes proposed by Henderson 
   ! and Eaton.  They minimize coupling with other modules.  Design of these 
   ! interfaces should be refined via review by other CAM developers.  
   ! Interface *_defaultopts() gets default values from the responsible 
   ! module (Expert) prior to namelist read.  
   ! Interface *_setopts() sends values to the responsible module (Expert) 
   ! after namelist read.  Erroneous values are handled by Experts.  
   ! TBH  9/8/03 
   !
   use phys_grid,        only: phys_grid_defaultopts, phys_grid_setopts

   use chem_surfvals,    only: chem_surfvals_readnl
   use check_energy,     only: check_energy_defaultopts, check_energy_setopts
   use cam_restart,      only: restart_defaultopts, restart_setopts, restart_printopts
   use carma_flags_mod,  only: carma_readnl
   use co2_cycle,        only: co2_cycle_readnl
   use shr_string_mod,   only: shr_string_toUpper
   use scamMod,          only: scam_setopts,scam_default_opts

   ! Some modules read their own namelist input.
   use spmd_utils,          only: spmd_utils_readnl
   use physconst,           only: physconst_readnl
   use phys_control,        only: phys_ctl_readnl
   use wv_saturation,       only: wv_sat_readnl
   use ref_pres,            only: ref_pres_readnl
   use cam3_aero_data,      only: cam3_aero_data_readnl
   use cam3_ozone_data,     only: cam3_ozone_data_readnl
   use macrop_driver,       only: macrop_driver_readnl
   use microp_driver,       only: microp_driver_readnl
   use microp_aero,         only: microp_aero_readnl
   use subcol,              only: subcol_readnl
   use cloud_fraction,      only: cldfrc_readnl
   use cldfrc2m,            only: cldfrc2m_readnl
   use unicon_cam,          only: unicon_cam_readnl
   use cldwat,              only: cldwat_readnl
   use zm_conv,             only: zmconv_readnl
   use hk_conv,             only: hkconv_readnl
   use uwshcu,              only: uwshcu_readnl
   use pkg_cld_sediment,    only: cld_sediment_readnl
   use gw_drag,             only: gw_drag_readnl
   use qbo,                 only: qbo_readnl
   use iondrag,             only: iondrag_readnl
   use phys_debug_util,     only: phys_debug_readnl
   use rad_constituents,    only: rad_cnst_readnl
   use radiation_data,      only: rad_data_readnl
   use modal_aer_opt,       only: modal_aer_opt_readnl
   use clubb_intr,          only: clubb_readnl
   use chemistry,           only: chem_readnl
   use prescribed_volcaero, only: prescribed_volcaero_readnl
   use aerodep_flx,         only: aerodep_flx_readnl
   use solar_data,          only: solar_data_readnl
   use tropopause,          only: tropopause_readnl
   use aoa_tracers,         only: aoa_tracers_readnl
   use prescribed_ozone,    only: prescribed_ozone_readnl
   use prescribed_aero,     only: prescribed_aero_readnl
   use prescribed_ghg,      only: prescribed_ghg_readnl
   use aircraft_emit,       only: aircraft_emit_readnl
   use cospsimulator_intr,  only: cospsimulator_intr_readnl
   use sat_hist,            only: sat_hist_readnl
   ! Needed by sat_hist_readnl
   use cam_history,         only: hfilename_spec, mfilt, fincl, nhtfrq, avgflag_pertape
   use vertical_diffusion,  only: vd_readnl
   use cam_history_support, only: fieldname_len, fieldname_lenp2
   use cam_diagnostics,     only: diag_readnl
   use nudging,             only: nudging_readnl
   use radheat,             only: radheat_readnl
#if ( defined OFFLINE_DYN )
   use metdata,             only: metdata_readnl
#endif
   use radiation,           only: radiation_readnl

!---------------------------Arguments-----------------------------------

   logical , intent(in), optional :: single_column_in 
   real(r8), intent(in), optional :: scmlon_in
   real(r8), intent(in), optional :: scmlat_in
   character(len=*)    , optional :: nlfilename_in
!-----------------------------------------------------------------------

   include 'netcdf.inc'

!---------------------------Local variables-----------------------------
   character(len=*), parameter ::  subname = "read_namelist"

   integer ntspdy         ! number of timesteps per day
   integer t              ! history tape index
   integer lastchar       ! index to last char of a char variable
   integer ierr           ! error code
   integer unitn          ! namelist unit number

   integer f, i
   integer, parameter :: max_chars = 128


   ! Define the cam_inparm namelist
   ! ***NOTE*** If a namelist option is not described in the CAM Users Guide,
   !            it is not supported.
   namelist /cam_inparm/ ncdata, bnd_topo, &
                     cam_branch_file  , &
                     absems_data, &
                     dtime, &
                     nlvdry,  &
                     pertlim ,&
                     new_random ,&
                     seed_custom ,&
                     seed_clock ,&
                     readtrace, rayk0, raykrange, raytau0, &
                     tracers_flag, &
                     indirect, &
                     print_step_cost,  &
                     phys_alltoall, phys_loadbalance, phys_twin_algorithm, &
                     phys_chnk_per_thd

   ! physics buffer
   namelist /cam_inparm/ pbuf_global_allocate

   ! conservation checks
   namelist /cam_inparm/ print_energy_errors

   ! scam
   namelist /cam_inparm/ iopfile,scm_iop_srf_prop,scm_relaxation, &
                         scm_relaxation_low, scm_relaxation_high, &
                         scm_diurnal_avg,scm_crm_mode,scm_clubb_iop_name, &
                         scm_observed_aero,swrad_off,lwrad_off, precip_off

!-----------------------------------------------------------------------

   if (present(nlfilename_in)) then
      nlfilename = nlfilename_in
   end if

   ! Determine preset values (this is currently being phased out)
   call preset ()

   ! Preset sulfate aerosol related variables
   indirect  = .false.

   ! restart write interval
   call restart_defaultopts( &
      cam_branch_file_out          =cam_branch_file            )

   ! Get default values of runtime options for physics chunking.
   call phys_grid_defaultopts(                      &
      phys_loadbalance_out    =phys_loadbalance,    &
      phys_twin_algorithm_out =phys_twin_algorithm, &
      phys_alltoall_out       =phys_alltoall,       &
      phys_chnk_per_thd_out   =phys_chnk_per_thd    )

   ! conservation
   call check_energy_defaultopts( &
      print_energy_errors_out = print_energy_errors )

   ! Set default options for single column model
   if (present(single_column_in)) then
      call scam_default_opts(scmlat_out=scmlat,scmlon_out=scmlon, &
        single_column_out=single_column, &
        scm_iop_srf_prop_out=scm_iop_srf_prop,&
        scm_relaxation_out=scm_relaxation, &
        scm_relaxation_low_out=scm_relaxation_low, &
        scm_relaxation_high_out=scm_relaxation_high, &
        scm_diurnal_avg_out=scm_diurnal_avg, &
        scm_crm_mode_out=scm_crm_mode, &
        scm_observed_aero_out=scm_observed_aero, &
        swrad_off_out=swrad_off, &
        lwrad_off_out=lwrad_off, &
        precip_off_out=precip_off, &
        scm_clubb_iop_name_out=scm_clubb_iop_name)
   end if

   ! Read in the cam_inparm namelist from input filename
   if (masterproc) then
      write(iulog,*) 'Read in cam_inparm namelist from: ', trim(nlfilename)
      unitn = getunit()
      open( unitn, file=trim(nlfilename), status='old' )

      ! Look for cam_inparm group name in the input file.  If found, leave the
      ! file positioned at that namelist group.
      call find_group_name(unitn, 'cam_inparm', status=ierr)
      if (ierr == 0) then  ! found cam_inparm
         read(unitn, cam_inparm, iostat=ierr)  ! read the cam_inparm namelist group
         if (ierr /= 0) then
            call endrun( subname//':: namelist read returns an'// &
                          ' error condition for cam_inparm' )
         end if
      else
         call endrun(subname // ':: can''t find cam_inparm in file ' // trim(nlfilename))
      end if
      close( unitn )
      call freeunit( unitn )

      ! Check CASE namelist variable
      if (caseid==' ') then
         call endrun ('READ_NAMELIST: Namelist variable CASEID must be set')
      end if

      lastchar = len(caseid)
      if (caseid(lastchar:lastchar) /= ' ') then
         write(iulog,*)'READ_NAMELIST: CASEID must not exceed ', len(caseid)-1, ' characters'
         call endrun
      end if
   end if

#if ( defined SPMD )
   ! Scatter namelist data to all processes
   call distnl ( )
#endif

   ! restart write interval
   call restart_setopts( nsrest,            &
      cam_branch_file_in          =cam_branch_file            )


   ! Set runtime options for physics chunking.
   call phys_grid_setopts(                          &
       phys_loadbalance_in    =phys_loadbalance,    &
       phys_twin_algorithm_in =phys_twin_algorithm, &
       phys_alltoall_in       =phys_alltoall,       &
       phys_chnk_per_thd_in   =phys_chnk_per_thd    )

   ! conservation
   call check_energy_setopts( &
      print_energy_errors_in = print_energy_errors )

   ! Set runtime options for single column mode
   if (present(single_column_in) .and. present(scmlon_in) .and. present(scmlat_in)) then 
      if (single_column_in) then
         single_column = single_column_in
         scmlon = scmlon_in
         scmlat = scmlat_in
         call scam_setopts( scmlat_in=scmlat,scmlon_in=scmlon, &
                            iopfile_in=iopfile,single_column_in=single_column,&
                            scm_iop_srf_prop_in=scm_iop_srf_prop,&
                            scm_relaxation_in=scm_relaxation, &
                            scm_relaxation_low_in=scm_relaxation_low, &
                            scm_relaxation_high_in=scm_relaxation_high, &
                            scm_diurnal_avg_in=scm_diurnal_avg, &
                            scm_crm_mode_in=scm_crm_mode, &
                            scm_observed_aero_in=scm_observed_aero, &
                            swrad_off_in=swrad_off, &
                            lwrad_off_in=lwrad_off, &
                            precip_off_in=precip_off, &
                            scm_clubb_iop_name_in=scm_clubb_iop_name)
      end if
   endif

   ! Call subroutines for modules to read their own namelist.
   ! In some cases namelist default values may depend on settings from
   ! other modules, so there may be an order dependence in the following
   ! calls.
   ! ***N.B.*** In particular, physconst_readnl should be called before
   !            the other readnl methods in case that method is used to set
   !            physical constants, some of which are set at runtime
   !            by the physconst_readnl method.
   ! Modules that read their own namelist are responsible for making sure
   ! all processes receive the values.

   call spmd_utils_readnl(nlfilename)
   call history_readnl(nlfilename, dtime)
   call physconst_readnl(nlfilename)
   call chem_surfvals_readnl(nlfilename)
   call phys_ctl_readnl(nlfilename)
   call wv_sat_readnl(nlfilename)
   call ref_pres_readnl(nlfilename)
   call cam3_aero_data_readnl(nlfilename)
   call cam3_ozone_data_readnl(nlfilename)
   call macrop_driver_readnl(nlfilename)
   call microp_driver_readnl(nlfilename)
   call microp_aero_readnl(nlfilename)
   call clubb_readnl(nlfilename)
   call subcol_readnl(nlfilename)
   call cldfrc_readnl(nlfilename)
   call cldfrc2m_readnl(nlfilename)
   call unicon_cam_readnl(nlfilename)
   call zmconv_readnl(nlfilename)
   call cldwat_readnl(nlfilename)
   call hkconv_readnl(nlfilename)
   call uwshcu_readnl(nlfilename)
   call cld_sediment_readnl(nlfilename)
   call gw_drag_readnl(nlfilename)
   call qbo_readnl(nlfilename)
   call iondrag_readnl(nlfilename)
   call phys_debug_readnl(nlfilename)
   call rad_cnst_readnl(nlfilename)
   call rad_data_readnl(nlfilename)
   call modal_aer_opt_readnl(nlfilename)
   call chem_readnl(nlfilename)
   call prescribed_volcaero_readnl(nlfilename)
   call solar_data_readnl(nlfilename)
   call carma_readnl(nlfilename)
   call tropopause_readnl(nlfilename)
   call aoa_tracers_readnl(nlfilename)
   call aerodep_flx_readnl(nlfilename)
   call prescribed_ozone_readnl(nlfilename)
   call prescribed_aero_readnl(nlfilename)
   call prescribed_ghg_readnl(nlfilename)
   call co2_cycle_readnl(nlfilename)
   call aircraft_emit_readnl(nlfilename)
   call cospsimulator_intr_readnl(nlfilename)
   call sat_hist_readnl(nlfilename, hfilename_spec, mfilt, fincl, nhtfrq, avgflag_pertape)
   call diag_readnl(nlfilename)
   call nudging_readnl(nlfilename)
   call radheat_readnl(nlfilename)
   call vd_readnl(nlfilename)
#if ( defined OFFLINE_DYN )
   call metdata_readnl(nlfilename)
#endif

   ! Read radiation namelist
   call radiation_readnl(nlfilename, dtime_in=dtime)

   ! Print cam_inparm input variables to standard output
   if (masterproc) then
      write(iulog,*)' ------------------------------------------'
      write(iulog,*)'     *** INPUT VARIABLES (CAM_INPARM) ***'
      write(iulog,*)' ------------------------------------------'
      if (nsrest/=0) then
         write(iulog,*) '  Continuation of an earlier run'
      else
         write(iulog,*) '         Initial run'
      end if
      write(iulog,*) ' ********** CASE = ',trim(caseid),' **********'
      write(iulog,'(1x,a)') ctitle
      if (len_trim(ncdata) > 0) then
         write(iulog,*) 'Initial dataset is: ',trim(ncdata)
      end if
      write(iulog,*)'Topography dataset is: ', trim(bnd_topo)
      write(iulog,*)'Time-invariant (absorption/emissivity) factor dataset is: ', trim(absems_data)

      ! Type of run
      write(iulog,*)'Run type flag (NSREST) 0=initial, 1=restart, 3=branch ',nsrest

      call restart_printopts()

      ! Write physics variables from namelist cam_inparm to std. output
      write(iulog,9108) nlvdry
9108 format('Lowest level for dry adiabatic adjust (NLVDRY)',i10)

      if ( (adiabatic .and. ideal_phys) .or. (adiabatic .and. aqua_planet) .or. &
           (ideal_phys .and. aqua_planet) ) then
         call endrun ('READ_NAMELIST: Only one of ADIABATIC, IDEAL_PHYS, or AQUA_PLANET can be .true.')
      end if

#ifdef COUP_SOM
      if (adiabatic .or. ideal_phys .or. aqua_planet )then
         call endrun ('READ_NAMELIST: adiabatic, ideal_phys or aqua_planet can not be used with SOM')
      end if
#else
      if (adiabatic)   write(iulog,*) 'Model will run ADIABATICALLY (i.e. no physics)'
      if (ideal_phys)  write(iulog,*) 'Run ONLY the "idealized" dynamical core of the ', &
                                  'model  (dynamics + Held&Suarez-specified physics)'
      if (aqua_planet) then
         write(iulog,*) 'Running model in "AQUA_PLANET" mode'
      else
         write(iulog,*) 'NOT Running model in "AQUA_PLANET" mode'
      end if
#endif


   end if

   ! set public data in cam_control_mod
   moist_physics = (.not. adiabatic) .and. (.not. ideal_phys)

#ifdef PERGRO
   if (masterproc) then
      write(iulog,*)'pergro for cloud water is true'
   end if
#endif

   ntspdy = nint(86400._r8/dtime) ! no. timesteps per day


end subroutine read_namelist


!=======================================================================

#ifdef SPMD
subroutine distnl
!-----------------------------------------------------------------------
!     
! Purpose:     
! Distribute namelist data to all processors.
!
! The cpp SPMD definition provides for the funnelling of all program i/o
! through the master processor. Processor 0 either reads restart/history
! data from the disk and distributes it to all processors, or collects
! data from all processors and writes it to disk.
!     
!---------------------------Code history-------------------------------
!
! Original version:  CCM2
! Standardized:      J. Rosinski, Oct 1995
!                    J. Truesdale, Feb. 1996
!
!-----------------------------------------------------------------------
   use mpishorthand
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
! 
   call mpibcast (dtime,       1,mpiint,0,mpicom)
   call mpibcast (nsrest  ,1,mpiint,0,mpicom)
   call mpibcast (nlvdry  ,1,mpiint,0,mpicom)

   call mpibcast (rayk0    ,1,mpiint,0,mpicom)
   call mpibcast (raykrange,1,mpir8,0,mpicom)
   call mpibcast (raytau0  ,1,mpir8,0,mpicom)

   call mpibcast (tracers_flag,1,mpilog,0,mpicom)
   call mpibcast (readtrace   ,1,mpilog,0,mpicom)
   call mpibcast (adiabatic   ,1,mpilog,0,mpicom)
   call mpibcast (ideal_phys  ,1,mpilog,0,mpicom)
   call mpibcast (aqua_planet ,1,mpilog,0,mpicom)

   call mpibcast (print_step_cost,1,mpilog,0,mpicom)
   call mpibcast (pertlim     ,1, mpir8 , 0, mpicom )
   call mpibcast (new_random  ,1, mpilog, 0, mpicom )
   call mpibcast (seed_custom ,1, mpiint, 0, mpicom )
   call mpibcast (seed_clock  ,1, mpilog, 0, mpicom )

   call mpibcast (caseid  ,len(caseid) ,mpichar,0,mpicom)
   call mpibcast (ctitle  ,len(ctitle),mpichar,0,mpicom)
   call mpibcast (ncdata  ,len(ncdata) ,mpichar,0,mpicom)
   call mpibcast (bnd_topo  ,len(bnd_topo) ,mpichar,0,mpicom)
   call mpibcast (absems_data,len(absems_data),mpichar,0,mpicom)
   call mpibcast (cam_branch_file  ,len(cam_branch_file) ,mpichar,0,mpicom)

   call mpibcast (indirect     , 1 ,mpilog, 0,mpicom)

   ! Physics chunk tuning
   call mpibcast (phys_loadbalance   ,1,mpiint,0,mpicom)
   call mpibcast (phys_twin_algorithm,1,mpiint,0,mpicom)
   call mpibcast (phys_alltoall      ,1,mpiint,0,mpicom)
   call mpibcast (phys_chnk_per_thd  ,1,mpiint,0,mpicom)

   ! Physics buffer
   call mpibcast (pbuf_global_allocate, 1, mpilog, 0, mpicom)

   ! Conservation
   call mpibcast (print_energy_errors, 1, mpilog, 0, mpicom)

end subroutine distnl
#endif



subroutine preset
!----------------------------------------------------------------------- 
! 
! Purpose: Preset namelist CAM_INPARM input variables and initialize some other variables
! 
! Method: Hardwire the values
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
   use rgrid
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!-----------------------------------------------------------------------
!
!
! Flags
!
   print_step_cost = .false.   ! print per timestep cost info
!
! rgrid: set default to full grid
!
   nlon(:) = plon
!!
!! Unit numbers: set to invalid
!!
!   ncid_ini = -1
!   ncid_sst = -1
!   ncid_trc = -1
!
   return
end subroutine preset

end module runtime_opts
