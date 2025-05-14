module zm_aero_type
   !----------------------------------------------------------------------------
   ! Purpose: aerosol derived type for ZM microphysics
   ! Original Author: Xialiang Song and Guang Zhang, June 2010
   !----------------------------------------------------------------------------
#ifdef SCREAM_CONFIG_IS_CMAKE
   use zm_eamxx_bridge_params, only: r8
#else
   use shr_kind_mod,     only: r8=>shr_kind_r8
#endif

   public :: zm_aero_t     ! structure to hold aerosol state information for ZM microphysics

!===================================================================================================

! generic 2D pointer type for zm_aero_t
type, public :: ptr2d
   real(r8), pointer :: val(:,:)
end type ptr2d

! structure to hold aerosol state information for ZM microphysics
type :: zm_aero_t

   ! Aerosol treatment
   character(len=5) :: scheme  ! either 'bulk' or 'modal'

   ! Bulk aerosols
   integer :: nbulk    =  0 ! number of bulk aerosols affecting climate
   integer :: idxsul   = -1 ! index in aerosol list for sulfate
   integer :: idxdst1  = -1 ! index in aerosol list for dust1
   integer :: idxdst2  = -1 ! index in aerosol list for dust2
   integer :: idxdst3  = -1 ! index in aerosol list for dust3
   integer :: idxdst4  = -1 ! index in aerosol list for dust4
   integer :: idxbcphi = -1 ! index in aerosol list for Soot (BCPHI)

   real(r8),    allocatable :: num_to_mass_aer(:)  ! conversion of mmr to number conc for bulk aerosols
   type(ptr2d), allocatable :: mmr_bulk(:)         ! array of pointers to bulk aerosol mmr
   real(r8),    allocatable :: mmrg_bulk(:,:,:)    ! gathered bulk aerosol mmr

   ! Modal aerosols
   integer                  :: nmodes = 0      ! number of modes
   integer,     allocatable :: nspec(:)        ! number of species in each mode
   type(ptr2d), allocatable :: num_a(:)        ! number mixing ratio of modes (interstitial phase)
   type(ptr2d), allocatable :: mmr_a(:,:)      ! species mmr in each mode (interstitial phase)
   real(r8),    allocatable :: numg_a(:,:,:)   ! gathered number mixing ratio of modes (interstitial phase)
   real(r8),    allocatable :: mmrg_a(:,:,:,:) ! gathered species mmr in each mode (interstitial phase)
   real(r8),    allocatable :: voltonumblo(:)  ! volume to number conversion (lower bound) for each mode
   real(r8),    allocatable :: voltonumbhi(:)  ! volume to number conversion (upper bound) for each mode
   real(r8),    allocatable :: specdens(:,:)   ! density of modal species
   real(r8),    allocatable :: spechygro(:,:)  ! hygroscopicity of modal species

   integer :: mode_accum_idx  = -1  ! index of accumulation mode
   integer :: mode_aitken_idx = -1  ! index of aitken mode
   integer :: mode_coarse_idx = -1  ! index of coarse mode
   integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
   integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode

   type(ptr2d), allocatable :: dgnum(:)        ! mode dry radius
   real(r8),    allocatable :: dgnumg(:,:,:)   ! gathered mode dry radius

   real(r8) :: sigmag_aitken

end type zm_aero_t

!===================================================================================================

end module zm_aero_type
