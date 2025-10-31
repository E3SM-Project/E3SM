module zm_conv_types
   !----------------------------------------------------------------------------
   ! Purpose: utility methods for ZM deep convection scheme
   !----------------------------------------------------------------------------
#ifdef SCREAM_CONFIG_IS_CMAKE
   use zm_eamxx_bridge_params, only: r8, masterproc
#else
   use shr_kind_mod,     only: r8=>shr_kind_r8
   use spmd_utils,       only: masterproc
   use mpishorthand
#endif
   use cam_logfile,      only: iulog
   use shr_sys_mod,      only: shr_sys_flush

   public :: zm_const_set_to_global    ! set constant values using global values from physconst/shr_const_mod
   public :: zm_const_set_for_testing  ! set constant values consistent with shr_const_mod for testing
   public :: zm_const_print            ! print constant values for log file
   public :: zm_param_set_for_testing  ! set parameter values for testing
   public :: zm_param_mpi_broadcast    ! broadcast parameter values to all MPI ranks
   public :: zm_param_print            ! print parameter values for log file
   
   ! ZM derived types
   public :: zm_const_t ! derived type to hold ZM constants
   public :: zm_param_t ! derived type to hold ZM tunable parameters

   ! invalid values for parameter initialization
   real(r8), parameter :: unset_r8  = huge(1.0_r8)
   integer , parameter :: unset_int = huge(1)

!===================================================================================================

type :: zm_const_t
   !----------------------------------------------------------------------------
   ! Purpose: derived type to hold ZM constant values
   !----------------------------------------------------------------------------
   real(r8) :: pi      ! 3.14159...
   real(r8) :: grav    ! gravitational acceleration      [m/s2]
   real(r8) :: rgrav   ! 1/grav                          [s2/m]
   real(r8) :: rdair   ! gas constant for dry air        [J/K/kg]
   real(r8) :: rh2o    ! gas constant for water vapor    [J/K/kg]
   real(r8) :: zvir    ! virtual temperature parameter   []
   real(r8) :: cpair   ! specific heat of dry air        [J/K/kg]
   real(r8) :: cpwv    ! specific heat of water vapor    [J/K/kg]
   real(r8) :: cpliq   ! specific heat of liquid water   [J/K/kg]
   real(r8) :: tfreez  ! freezing point of water         [K]
   real(r8) :: latvap  ! latent heat of vaporization     [J/kg]
   real(r8) :: latice  ! latent heat of fusion           [J/kg]
   real(r8) :: epsilo  ! ratio of h2o to dry air molecular weights
end type zm_const_t

!===================================================================================================

type :: zm_param_t
   !----------------------------------------------------------------------------
   ! Purpose: derived type to hold ZM tunable parameters
   !----------------------------------------------------------------------------
   real(r8) :: tau             = unset_r8    ! convective adjustment time scale
   real(r8) :: alfa            = unset_r8    ! max downdraft mass flux fraction
   real(r8) :: ke              = unset_r8    ! evaporation efficiency
   real(r8) :: dmpdz           = unset_r8    ! fractional mass entrainment rate [1/m]
   logical  :: tpert_fix       = .false.     ! flag to disable using applying tpert to PBL-rooted convection
   real(r8) :: tpert_fac       = 0           ! tunable temperature perturbation factor
   real(r8) :: tiedke_add      = unset_r8    ! tunable temperature perturbation
   real(r8) :: c0_lnd          = unset_r8    ! autoconversion coefficient over land
   real(r8) :: c0_ocn          = unset_r8    ! autoconversion coefficient over ocean
   integer  :: num_cin         = 0           ! num of neg buoyancy regions allowed before the conv top and CAPE calc are completed
   integer  :: limcnv          = 0           ! upper pressure interface level to limit deep convection
   integer  :: mx_bot_lyr_adj  = unset_int   ! bot layer index adjustment for launch level search
   logical  :: trig_dcape      = .false.     ! true if to using DCAPE trigger - based on CAPE generation by the dycor
   logical  :: trig_ull        = .false.     ! true if to using the "unrestricted launch level" (ULL) mode
   logical  :: clos_dyn_adj    = .false.     ! flag for mass flux adjustment to CAPE closure
   logical  :: no_deep_pbl     = .false.     ! flag to eliminate deep convection within PBL
   ! ZM micro parameters
   logical  :: zm_microp       = .false.     ! switch for convective microphysics
   real(r8) :: auto_fac        = unset_r8    ! ZM microphysics enhancement factor for droplet-rain autoconversion
   real(r8) :: accr_fac        = unset_r8    ! ZM microphysics enhancement factor for droplet-rain accretion
   real(r8) :: micro_dcs       = unset_r8    ! ZM microphysics size threshold for cloud ice to snow autoconversion [m]
   ! MCSP parameters
   logical  :: mcsp_enabled    = .false.     ! flag for mesoscale coherent structure parameterization (MSCP)
   real(r8) :: mcsp_t_coeff    = 0           ! MCSP coefficient for temperature tendencies
   real(r8) :: mcsp_q_coeff    = 0           ! MCSP coefficient for specific humidity tendencies
   real(r8) :: mcsp_u_coeff    = 0           ! MCSP coefficient for zonal momentum tendencies
   real(r8) :: mcsp_v_coeff    = 0           ! MCSP coefficient for meridional momentum tendencies
end type zm_param_t

!===================================================================================================
contains
!===================================================================================================

subroutine zm_const_set_to_global(zm_const)
   !----------------------------------------------------------------------------
   ! Purpose: set constant values using global values from physconst/shr_const_mod
   !----------------------------------------------------------------------------
#ifdef SCREAM_CONFIG_IS_CMAKE
   use zm_eamxx_bridge_physconst, only: pi,gravit,rair,cpair,cpwv,cpliq,rh2o,tmelt,latvap,latice,epsilo
#else
   use physconst, only: pi,gravit,rair,cpair,cpwv,cpliq,rh2o,tmelt,latvap,latice,epsilo
#endif
   !----------------------------------------------------------------------------
   type(zm_const_t), intent(inout) :: zm_const
   !----------------------------------------------------------------------------
   zm_const%pi       = pi
   zm_const%grav     = gravit
   zm_const%rgrav    = 1.0_r8/gravit
   zm_const%rdair    = rair
   zm_const%cpair    = cpair
   zm_const%cpwv     = cpwv
   zm_const%cpliq    = cpliq
   zm_const%rh2o     = rh2o
   zm_const%tfreez   = tmelt
   zm_const%latvap   = latvap
   zm_const%latice   = latice
   zm_const%epsilo   = epsilo
   zm_const%zvir     = 1.608_r8 ! use this instead of physconst value to avoid non-BFB diffs
end subroutine zm_const_set_to_global

!===================================================================================================

subroutine zm_const_set_for_testing(zm_const)
   !----------------------------------------------------------------------------
   ! Purpose: set constant values consistent with shr_const_mod for testing
   ! Note - normal model operation should use zm_const_set_to_global()
   !----------------------------------------------------------------------------
   type(zm_const_t), intent(inout) :: zm_const
   !----------------------------------------------------------------------------
   ! Local variables
   real(r8) :: boltz    ! Boltzmann's constant            [J/K/molecule]
   real(r8) :: avogad   ! Avogadro's number               [molecules/kmole]
   real(r8) :: rgas     ! Universal gas constant          [J/K/kmole]
   real(r8) :: mwdair   ! molecular weight dry air        [kg/kmole]
   real(r8) :: mwwv     ! molecular weight water vapor    [kg/kmole]
   !----------------------------------------------------------------------------
   boltz    = 1.38065e-23_r8
   avogad   = 6.02214e26_r8
   rgas     = avogad*boltz
   mwdair   = 28.966_r8
   mwwv     = 18.016_r8
   zm_const%pi       = 3.14159265358979323846_R8
   zm_const%grav     = 9.80616_r8
   zm_const%rgrav    = 1.0_r8/9.80616_r8
   zm_const%rdair    = rgas/mwdair
   zm_const%rh2o     = rgas/mwwv
   zm_const%zvir     = zm_const%rh2o/zm_const%rdair - 1.0_r8
   zm_const%cpair    = 1.00464e3_r8
   zm_const%cpwv     = 1.810e3_r8
   zm_const%cpliq    = 4.188e3_r8
   zm_const%tfreez   = 273.15_r8
   zm_const%latvap   = 2.501e6_r8
   zm_const%latice   = 3.337e5_r8
   zm_const%epsilo   = mwwv/mwdair
end subroutine zm_const_set_for_testing

!===================================================================================================

subroutine zm_const_print(zm_const)
   !----------------------------------------------------------------------------
   ! Purpose: print constant values for log file
   !----------------------------------------------------------------------------
   type(zm_const_t), intent(in) :: zm_const
   !----------------------------------------------------------------------------
   character(len=32) :: indent = '  '
   !----------------------------------------------------------------------------
   if (masterproc) then
      write(iulog,*) ''
      write(iulog,*) 'ZM deep convection constant values:'
      write(iulog,*) indent,'pi              : ',zm_const%pi
      write(iulog,*) indent,'grav            : ',zm_const%grav
      write(iulog,*) indent,'rgrav           : ',zm_const%rgrav
      write(iulog,*) indent,'rdair           : ',zm_const%rdair
      write(iulog,*) indent,'rh2o            : ',zm_const%rh2o
      write(iulog,*) indent,'zvir            : ',zm_const%zvir
      write(iulog,*) indent,'cpair           : ',zm_const%cpair
      write(iulog,*) indent,'cpwv            : ',zm_const%cpwv
      write(iulog,*) indent,'cpliq           : ',zm_const%cpliq
      write(iulog,*) indent,'tfreez          : ',zm_const%tfreez
      write(iulog,*) indent,'latvap          : ',zm_const%latvap
      write(iulog,*) indent,'latice          : ',zm_const%latice
      write(iulog,*) indent,'epsilo          : ',zm_const%epsilo
      write(iulog,*) ''
      call shr_sys_flush(iulog)
   end if ! masterproc
end subroutine zm_const_print

!===================================================================================================

subroutine zm_param_set_for_testing(zm_param)
   !----------------------------------------------------------------------------
   ! Purpose: set parameter values for testing
   !----------------------------------------------------------------------------
   type(zm_param_t), intent(inout) :: zm_param
   !----------------------------------------------------------------------------
   zm_param%tau             = 3600.0_r8
   zm_param%alfa            = 0.14D0
   zm_param%ke              = 2.5E-6
   zm_param%dmpdz           = -0.7e-3
   zm_param%tpert_fix       = .true.
   zm_param%tpert_fac       = 2.0D0
   zm_param%tiedke_add      = 0.8D0
   zm_param%c0_lnd          = 0.0020D0
   zm_param%c0_ocn          = 0.0020D0
   zm_param%num_cin         = 1
   zm_param%limcnv          = 23 ! note - default for E3SMv3 => ne30pg2 w/ L80 
   zm_param%mx_bot_lyr_adj  = 1
   zm_param%trig_dcape      = .true.
   zm_param%trig_ull        = .true.
   zm_param%clos_dyn_adj    = .true.
   zm_param%no_deep_pbl     = .false.
   ! ZM micro parameters
   zm_param%zm_microp       = .true.
   zm_param%auto_fac        = 7.0D0
   zm_param%accr_fac        = 1.5D0
   zm_param%micro_dcs       = 150.E-6
   ! MCSP parameters
   zm_param%mcsp_enabled    = .true.
   zm_param%mcsp_t_coeff    = 0.3
   zm_param%mcsp_q_coeff    = 0
   zm_param%mcsp_u_coeff    = 0
   zm_param%mcsp_v_coeff    = 0
end subroutine zm_param_set_for_testing

!===================================================================================================

subroutine zm_param_mpi_broadcast(zm_param)
   !----------------------------------------------------------------------------
   ! Purpose: broadcast parameter values to all MPI ranks
   !----------------------------------------------------------------------------
   type(zm_param_t), intent(inout) :: zm_param
   !----------------------------------------------------------------------------
#ifdef SCREAM_CONFIG_IS_CMAKE
   ! The EAMXx bridge to ZM will rely on zm_param_set_for_testing()
   ! so MPI broadcasting of ZM parameter is unnecessary
#else
   call mpibcast(zm_param%tau,               1, mpir8,  0, mpicom)
   call mpibcast(zm_param%alfa,              1, mpir8,  0, mpicom)
   call mpibcast(zm_param%ke,                1, mpir8,  0, mpicom)
   call mpibcast(zm_param%dmpdz,             1, mpir8,  0, mpicom)
   call mpibcast(zm_param%tpert_fix,         1, mpilog, 0, mpicom)
   call mpibcast(zm_param%tpert_fac,         1, mpir8,  0, mpicom)
   call mpibcast(zm_param%tiedke_add,        1, mpir8,  0, mpicom)
   call mpibcast(zm_param%c0_lnd,            1, mpir8,  0, mpicom)
   call mpibcast(zm_param%c0_ocn,            1, mpir8,  0, mpicom)
   call mpibcast(zm_param%num_cin,           1, mpiint, 0, mpicom)
   call mpibcast(zm_param%limcnv,            1, mpiint, 0, mpicom)
   call mpibcast(zm_param%mx_bot_lyr_adj,    1, mpiint, 0, mpicom)
   call mpibcast(zm_param%trig_dcape,        1, mpilog, 0, mpicom)
   call mpibcast(zm_param%trig_ull,          1, mpilog, 0, mpicom)
   call mpibcast(zm_param%clos_dyn_adj,      1, mpilog, 0, mpicom)
   call mpibcast(zm_param%no_deep_pbl,       1, mpilog, 0, mpicom)
   call mpibcast(zm_param%zm_microp,         1, mpilog, 0, mpicom) ! ZM micro parameters
   call mpibcast(zm_param%auto_fac,          1, mpir8,  0, mpicom)
   call mpibcast(zm_param%accr_fac,          1, mpir8,  0, mpicom)
   call mpibcast(zm_param%micro_dcs,         1, mpir8,  0, mpicom)
   call mpibcast(zm_param%mcsp_enabled,      1, mpilog, 0, mpicom) ! MCSP parameters
   call mpibcast(zm_param%mcsp_t_coeff,      1, mpir8,  0, mpicom)
   call mpibcast(zm_param%mcsp_q_coeff,      1, mpir8,  0, mpicom)
   call mpibcast(zm_param%mcsp_u_coeff,      1, mpir8,  0, mpicom)
   call mpibcast(zm_param%mcsp_v_coeff,      1, mpir8,  0, mpicom)
#endif
end subroutine zm_param_mpi_broadcast

!===================================================================================================

subroutine zm_param_print(zm_param)
   !----------------------------------------------------------------------------
   ! Purpose: print parameter values for log file
   !----------------------------------------------------------------------------
   type(zm_param_t), intent(in) :: zm_param
   !----------------------------------------------------------------------------
   character(len=32) :: indent = '  '
   !----------------------------------------------------------------------------
   if (masterproc) then
      write(iulog,*) ''
      write(iulog,*) 'ZM deep convection parameter values:'
      write(iulog,*) indent,'tau             : ',zm_param%tau
      write(iulog,*) indent,'alfa            : ',zm_param%alfa
      write(iulog,*) indent,'ke              : ',zm_param%ke
      write(iulog,*) indent,'dmpdz           : ',zm_param%dmpdz
      write(iulog,*) indent,'tpert_fix       : ',zm_param%tpert_fix
      write(iulog,*) indent,'tpert_fac       : ',zm_param%tpert_fac
      write(iulog,*) indent,'tiedke_add      : ',zm_param%tiedke_add
      write(iulog,*) indent,'c0_lnd          : ',zm_param%c0_lnd
      write(iulog,*) indent,'c0_ocn          : ',zm_param%c0_ocn
      write(iulog,*) indent,'num_cin         : ',zm_param%num_cin
      write(iulog,*) indent,'limcnv          : ',zm_param%limcnv
      write(iulog,*) indent,'mx_bot_lyr_adj  : ',zm_param%mx_bot_lyr_adj
      write(iulog,*) indent,'trig_dcape      : ',zm_param%trig_dcape
      write(iulog,*) indent,'trig_ull        : ',zm_param%trig_ull
      write(iulog,*) indent,'clos_dyn_adj    : ',zm_param%clos_dyn_adj
      write(iulog,*) indent,'no_deep_pbl     : ',zm_param%no_deep_pbl
      ! ZM micro parameters
      write(iulog,*) indent,'zm_microp       : ',zm_param%zm_microp
      write(iulog,*) indent,'auto_fac        : ',zm_param%auto_fac
      write(iulog,*) indent,'accr_fac        : ',zm_param%accr_fac
      write(iulog,*) indent,'micro_dcs       : ',zm_param%micro_dcs
      ! MCSP parameters
      write(iulog,*) indent,'mcsp_enabled    : ',zm_param%mcsp_enabled
      write(iulog,*) indent,'mcsp_t_coeff    : ',zm_param%mcsp_t_coeff
      write(iulog,*) indent,'mcsp_q_coeff    : ',zm_param%mcsp_q_coeff
      write(iulog,*) indent,'mcsp_u_coeff    : ',zm_param%mcsp_u_coeff
      write(iulog,*) indent,'mcsp_v_coeff    : ',zm_param%mcsp_v_coeff
      write(iulog,*) ''
      call shr_sys_flush(iulog)
   end if ! masterproc
end subroutine zm_param_print

!===================================================================================================

end module zm_conv_types
