module zm_conv_types
   !----------------------------------------------------------------------------
   ! Purpose: utility methods for ZM deep convection scheme
   !----------------------------------------------------------------------------
   use shr_kind_mod,     only: r8=>shr_kind_r8

   public :: zm_const_set_to_global    ! set zm_const using global values from physconst/shr_const_mod
   public :: zm_const_set_for_testing  ! set zm_const consistent with shr_const_mod for testing
   public :: zm_param_mpi_broadcast    ! broadcast parameter values to all MPI ranks
   
   ! ZM derived types
   public :: zm_const_t ! derived type to hold ZM constants
   public :: zm_param_t ! derived type to hold ZM tunable parameters

   ! invalid values for parameter initialization
   real(r8), parameter :: unset_r8  = huge(1.0_r8)
   integer , parameter :: unset_int = huge(1)

!===================================================================================================

type :: zm_const_t
   real(r8) :: grav    ! gravitational acceleration      [m/s2]
   real(r8) :: boltz   ! Boltzmann's constant            [J/K/molecule]
   real(r8) :: avogad  ! Avogadro's number               [molecules/kmole]
   real(r8) :: rgas    ! Universal gas constant          [J/K/kmole]
   real(r8) :: mwdair  ! molecular weight dry air        [kg/kmole]
   real(r8) :: mwwv    ! molecular weight water vapor    [kg/kmole]
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
   logical  :: trig_dcape      = .false.     ! true if to using DCAPE trigger - based on CAPE generation by the dycor
   logical  :: trig_ull        = .false.     ! true if to using the "unrestricted launch level" (ULL) mode
   integer  :: num_cin         = 0           ! num of neg buoyancy regions allowed before the conv top and CAPE calc are completed
   integer  :: limcnv          = 0           ! upper pressure interface level to limit deep convection
   logical  :: tpert_fix       = .false.     ! flag to disable using applying tpert to PBL-rooted convection
   real(r8) :: tpert_fac       = 0           ! tunable temperature perturbation factor
   real(r8) :: dmpdz           = unset_r8    ! fractional mass entrainment rate [1/m]
   real(r8) :: tiedke_add      = unset_r8    ! tunable temperature perturbation
   integer  :: mx_bot_lyr_adj  = unset_int   ! bot layer index adjustment for launch level search
end type zm_param_t

!===================================================================================================
contains
!===================================================================================================

subroutine zm_const_set_to_global(zm_const)
   !----------------------------------------------------------------------------
   ! Purpose: set zm_const using global values from physconst/shr_const_mod
   !----------------------------------------------------------------------------
   use physconst, only: gravit,rair,cpair,cpwv,cpliq,rh2o,tmelt,latvap,latice,epsilo
   !----------------------------------------------------------------------------
   type(zm_const_t), intent(inout) :: zm_const
   !----------------------------------------------------------------------------
   zm_const%grav     = gravit
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
   ! Purpose: set zm_const consistent with shr_const_mod for testing
   ! Note - normal model operation should use zm_const_set_to_global()
   !----------------------------------------------------------------------------
   type(zm_const_t), intent(inout) :: zm_const
   !----------------------------------------------------------------------------
   zm_const%grav     = 9.80616_r8
   zm_const%boltz    = 1.38065e-23_r8
   zm_const%avogad   = 6.02214e26_r8
   zm_const%rgas     = zm_const%avogad*zm_const%boltz
   zm_const%mwdair   = 28.966_r8
   zm_const%mwwv     = 18.016_r8
   zm_const%rdair    = zm_const%rgas/zm_const%mwdair
   zm_const%rh2o     = zm_const%rgas/zm_const%mwwv
   zm_const%zvir     = zm_const%rh2o/zm_const%rdair - 1.0_r8
   zm_const%cpair    = 1.00464e3_r8
   zm_const%cpwv     = 1.810e3_r8
   zm_const%cpliq    = 4.188e3_r8
   zm_const%tfreez   = 273.15_r8
   zm_const%latvap   = 2.501e6_r8
   zm_const%latice   = 3.337e5_r8
   zm_const%epsilo   = zm_const%mwwv/zm_const%mwdair
end subroutine zm_const_set_for_testing

!===================================================================================================

subroutine zm_param_mpi_broadcast(zm_param)
   !----------------------------------------------------------------------------
   ! Purpose: broadcast parameter values to all MPI ranks
   !----------------------------------------------------------------------------
   use mpishorthand
   !----------------------------------------------------------------------------
   type(zm_param_t), intent(inout) :: zm_param
   !----------------------------------------------------------------------------
   call mpibcast(zm_param%trig_dcape,        1, mpilog, 0, mpicom)
   call mpibcast(zm_param%trig_ull,          1, mpilog, 0, mpicom)
   call mpibcast(zm_param%tiedke_add,        1, mpir8,  0, mpicom)
   call mpibcast(zm_param%dmpdz,             1, mpir8,  0, mpicom)
   call mpibcast(zm_param%num_cin,           1, mpiint, 0, mpicom)
   call mpibcast(zm_param%mx_bot_lyr_adj,    1, mpiint, 0, mpicom)
   call mpibcast(zm_param%tpert_fix,         1, mpilog, 0, mpicom)
   call mpibcast(zm_param%tpert_fac,         1, mpir8,  0, mpicom)
end subroutine zm_param_mpi_broadcast

!===================================================================================================

end module zm_conv_types
