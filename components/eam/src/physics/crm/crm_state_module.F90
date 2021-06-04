module crm_state_module
   use openacc_utils
   use params_kind,  only: crm_rknd

   implicit none
   private

   public crm_state_type
   public crm_state_initialize
   public crm_state_finalize

   !------------------------------------------------------------------------------------------------
   type crm_state_type
      ! Purpose: Define a type that will encapsulate the CRM-level data that needs to be passed
      ! between the GCM and the CRM, and to encapulate operations on that data (i.e.,
      ! initialization, writing fields to netCDF files, displaying info, et.c).

      !---------------------------------------------------------------------------------------------
      ! CRM-scale fields

      ! NOTE: these were intent(inout) before, so these need to persist across calls; pointers so
      ! they can be used without making a bunch of temporary arrays. Dimensions should be
      ! (ncol,crm_nx,crm_ny,crm_nz)
      real(crm_rknd), allocatable :: u_wind(:,:,:,:)       ! CRM u-wind component
      real(crm_rknd), allocatable :: v_wind(:,:,:,:)       ! CRM v-wind component
      real(crm_rknd), allocatable :: w_wind(:,:,:,:)       ! CRM w-wind component
      real(crm_rknd), allocatable :: temperature(:,:,:,:)  ! CRM temperuture
      real(crm_rknd), allocatable :: qt(:,:,:,:)           ! CRM total water

      ! 2-moment microphsics variables
      real(crm_rknd), allocatable :: qc(:,:,:,:)   ! mass mixing ratio of cloud water
      real(crm_rknd), allocatable :: nc(:,:,:,:)   ! number concentration of cloud water
      real(crm_rknd), allocatable :: qr(:,:,:,:)   ! mass mixing ratio of rain
      real(crm_rknd), allocatable :: nr(:,:,:,:)   ! number concentration of rain
      real(crm_rknd), allocatable :: qi(:,:,:,:)   ! mass mixing ratio of cloud ice
      real(crm_rknd), allocatable :: ni(:,:,:,:)   ! number concentration of cloud ice
      real(crm_rknd), allocatable :: qs(:,:,:,:)   ! mass mixing ratio of snow
      real(crm_rknd), allocatable :: ns(:,:,:,:)   ! number concentration of snow
      real(crm_rknd), allocatable :: qg(:,:,:,:)   ! mass mixing ratio of graupel
      real(crm_rknd), allocatable :: ng(:,:,:,:)   ! number concentration of graupel
      
      ! 1-moment microphsics variables
      real(crm_rknd), allocatable :: qp(:,:,:,:)   ! mass mixing ratio of precipitating condensate
      real(crm_rknd), allocatable :: qn(:,:,:,:)   ! mass mixing ratio of cloud condensate

   end type crm_state_type
   !------------------------------------------------------------------------------------------------
contains

   !------------------------------------------------------------------------------------------------
   ! Type-bound procedures for crm_state_type
   subroutine crm_state_initialize(state,ncrms,crm_nx,crm_ny,crm_nz,MMF_microphysics_scheme)
      type(crm_state_type), intent(inout) :: state
      integer,              intent(in   ) :: ncrms, crm_nx, crm_ny, crm_nz
      character(len=*),    intent(in   ) :: MMF_microphysics_scheme    ! CRM microphysics scheme

      ! Allocate memory
      if (.not. allocated(state%u_wind))      allocate(state%u_wind(ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(state%v_wind))      allocate(state%v_wind(ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(state%w_wind))      allocate(state%w_wind(ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(state%temperature)) allocate(state%temperature(ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(state%qt))          allocate(state%qt(ncrms,crm_nx,crm_ny,crm_nz))

      call prefetch(state%u_wind)
      call prefetch(state%v_wind)
      call prefetch(state%w_wind)
      call prefetch(state%temperature)
      call prefetch(state%qt)

      if (trim(MMF_microphysics_scheme) .eq. 'm2005') then
         if (.not. allocated(state%qc))          allocate(state%qc(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%qi))          allocate(state%qi(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%qr))          allocate(state%qr(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%qs))          allocate(state%qs(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%qg))          allocate(state%qg(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%nc))          allocate(state%nc(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%ni))          allocate(state%ni(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%nr))          allocate(state%nr(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%ns))          allocate(state%ns(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%ng))          allocate(state%ng(ncrms,crm_nx,crm_ny,crm_nz))
         call prefetch(state%qc)
         call prefetch(state%qi)
         call prefetch(state%qr)
         call prefetch(state%qs)
         call prefetch(state%qg)
         call prefetch(state%nc)
         call prefetch(state%ni)
         call prefetch(state%nr)
         call prefetch(state%ns)
         call prefetch(state%ng)
      end if
      if (trim(MMF_microphysics_scheme) .eq. 'sam1mom') then
         if (.not. allocated(state%qp))          allocate(state%qp(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%qn))          allocate(state%qn(ncrms,crm_nx,crm_ny,crm_nz))
         call prefetch(state%qp)
         call prefetch(state%qn)
      end if

   end subroutine crm_state_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_state_finalize(state, MMF_microphysics_scheme)
      type(crm_state_type), intent(inout) :: state
      character(len=*), intent(in) :: MMF_microphysics_scheme    ! CRM microphysics scheme

      ! Nullify pointers
      if (allocated(state%u_wind))      deallocate(state%u_wind)
      if (allocated(state%v_wind))      deallocate(state%v_wind)
      if (allocated(state%w_wind))      deallocate(state%w_wind)
      if (allocated(state%temperature)) deallocate(state%temperature)
      if (allocated(state%qt))          deallocate(state%qt)

      if (trim(MMF_microphysics_scheme) .eq. 'm2005') then
         if (allocated(state%qc)) deallocate(state%qc)
         if (allocated(state%qi)) deallocate(state%qi)
         if (allocated(state%qr)) deallocate(state%qr)
         if (allocated(state%qs)) deallocate(state%qs)
         if (allocated(state%qg)) deallocate(state%qg)
         if (allocated(state%nc)) deallocate(state%nc)
         if (allocated(state%ni)) deallocate(state%ni)
         if (allocated(state%nr)) deallocate(state%nr)
         if (allocated(state%ns)) deallocate(state%ns)
         if (allocated(state%ng)) deallocate(state%ng)
      end if
      if (trim(MMF_microphysics_scheme) .eq. 'sam1mom') then
         if (allocated(state%qp)) deallocate(state%qp)
         if (allocated(state%qn)) deallocate(state%qn)
      end if

   end subroutine crm_state_finalize
end module crm_state_module
