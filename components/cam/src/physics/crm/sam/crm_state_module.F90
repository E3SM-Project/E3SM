module crm_state_module

   use params,       only: crm_rknd
   use crmdims,      only: crm_nx, crm_ny, crm_nz

   implicit none
   private

   public crm_state_type

   !------------------------------------------------------------------------------------------------
   type crm_state_type
      ! Purpose: Define a class that will encapulate the CRM-level data that needs to be passed
      ! between the GCM and the CRM, and to encapulate operations on that data (i.e.,
      ! initialization, writing fields to netCDF files, displaying info, et.c).

      ! TODO: should these be allocatable, and initialized once per node (at, maybe,
      ! crm_physics_init)?
      !---------------------------------------------------------------------------------------------
      ! CRM-scale fields

      ! NOTE: these were intent(inout) before, so these need to persist across calls; pointers so
      ! they can be used without making a bunch of temporary arrays. Dimensions should be
      ! (ncol,crm_nx,crm_ny,crm_nz)
      real(crm_rknd), pointer :: u_wind(:,:,:,:)       ! CRM u-wind component
      real(crm_rknd), pointer :: v_wind(:,:,:,:)       ! CRM v-wind component
      real(crm_rknd), pointer :: w_wind(:,:,:,:)       ! CRM w-wind component
      real(crm_rknd), pointer :: temperature(:,:,:,:)  ! CRM temperuture

      ! Microphysics
      real(crm_rknd), pointer :: qt(:,:,:,:) 

      ! NOTE: These are terrible variable names...replace with more descriptive names.
      ! for m2005...
      real(crm_rknd), pointer :: nc(:,:,:,:)
      real(crm_rknd), pointer :: qr(:,:,:,:)
      real(crm_rknd), pointer :: nr(:,:,:,:)
      real(crm_rknd), pointer :: qi(:,:,:,:)
      real(crm_rknd), pointer :: ni(:,:,:,:)
      real(crm_rknd), pointer :: qs(:,:,:,:)
      real(crm_rknd), pointer :: ns(:,:,:,:)
      real(crm_rknd), pointer :: qg(:,:,:,:)
      real(crm_rknd), pointer :: ng(:,:,:,:)
      real(crm_rknd), pointer :: qc(:,:,:,:)

      ! for sam1mom...
      real(crm_rknd), pointer :: qp(:,:,:,:)
      real(crm_rknd), pointer :: qn(:,:,:,:)

   contains
      ! Type-bound procedures. Initialization should nullify fields
      procedure, public :: initialize=>crm_state_initialize
      procedure, public :: finalize=>crm_state_finalize
      !procedure, public :: dump=>crm_state_dump

   end type crm_state_type
   !------------------------------------------------------------------------------------------------
contains

   !------------------------------------------------------------------------------------------------
   ! Type-bound procedures for crm_state_type
   subroutine crm_state_initialize(this)
      class(crm_state_type), intent(inout) :: this

      ! Nullify pointers
      this%u_wind => null()
      this%v_wind => null()
      this%w_wind => null()
      this%temperature => null()

      this%qt => null()
      this%qc => null()
      this%qi => null()
      this%qr => null()
      this%qs => null()
      this%qg => null()
      this%nc => null()
      this%ni => null()
      this%nr => null()
      this%ns => null()
      this%ng => null()

      this%qp => null()
      this%qn => null()

   end subroutine crm_state_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_state_finalize(this)
      class(crm_state_type), intent(inout) :: this

      ! Nullify pointers
      this%u_wind => null()
      this%v_wind => null()
      this%w_wind => null()
      this%temperature => null()

#ifdef m2005
      this%qt => null()
      this%qc => null()
      this%qi => null()
      this%qr => null()
      this%qs => null()
      this%qg => null()
      this%nc => null()
      this%ni => null()
      this%nr => null()
      this%ns => null()
      this%ng => null()
#else
      this%qt => null()
      this%qp => null()
      this%qn => null()
#endif

   end subroutine crm_state_finalize
end module crm_state_module
