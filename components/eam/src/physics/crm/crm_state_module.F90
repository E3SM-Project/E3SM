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
      real(crm_rknd), allocatable :: temperature(:,:,:,:)  ! CRM temperature
      real(crm_rknd), allocatable :: qv(:,:,:,:)           ! CRM water vapor

      ! 1-moment microphsics variables
      real(crm_rknd), allocatable :: qp(:,:,:,:)   ! mass mixing ratio of precipitating condensate
      real(crm_rknd), allocatable :: qn(:,:,:,:)   ! mass mixing ratio of cloud condensate

      ! 2-moment microphysics variables (p3)
      real(crm_rknd), allocatable :: qc(:,:,:,:)   ! mass mixing ratio of cloud water
      real(crm_rknd), allocatable :: nc(:,:,:,:)   ! number concentration of cloud water
      real(crm_rknd), allocatable :: qr(:,:,:,:)   ! mass mixing ratio of rain
      real(crm_rknd), allocatable :: nr(:,:,:,:)   ! number concentration of rain
      real(crm_rknd), allocatable :: qi(:,:,:,:)   ! mass mixing ratio of cloud ice
      real(crm_rknd), allocatable :: ni(:,:,:,:)   ! number concentration of cloud ice

      ! p3 microphysics variables not included above
      real(crm_rknd), allocatable :: qm(:,:,:,:) ! averaged riming density
      real(crm_rknd), allocatable :: bm(:,:,:,:) ! averaged riming volume

      ! "previous" state variables needed for P3
      real(crm_rknd), allocatable :: t_prev(:,:,:,:)  ! previous CRM time step temperature
      real(crm_rknd), allocatable :: q_prev(:,:,:,:) ! previous CRM time step water vapor

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
      if (.not. allocated(state%qv))          allocate(state%qv(ncrms,crm_nx,crm_ny,crm_nz))

      call prefetch(state%u_wind)
      call prefetch(state%v_wind)
      call prefetch(state%w_wind)
      call prefetch(state%temperature)
      call prefetch(state%qv)

      if (trim(MMF_microphysics_scheme) .eq. 'sam1mom') then
         if (.not. allocated(state%qp))          allocate(state%qp(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%qn))          allocate(state%qn(ncrms,crm_nx,crm_ny,crm_nz))
         call prefetch(state%qp)
         call prefetch(state%qn)
      end if

      if (trim(MMF_microphysics_scheme).eq.'p3') then
         if (.not. allocated(state%qc))          allocate(state%qc(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%qi))          allocate(state%qi(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%qr))          allocate(state%qr(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%nc))          allocate(state%nc(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%ni))          allocate(state%ni(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%nr))          allocate(state%nr(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%qm))          allocate(state%qm(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%bm))          allocate(state%bm(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%t_prev))      allocate(state%t_prev(ncrms,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(state%q_prev))      allocate(state%q_prev(ncrms,crm_nx,crm_ny,crm_nz))
         call prefetch(state%qc)
         call prefetch(state%qi)
         call prefetch(state%qr)
         call prefetch(state%nc)
         call prefetch(state%ni)
         call prefetch(state%nr)
         call prefetch(state%qm)
         call prefetch(state%bm)
         call prefetch(state%t_prev)
         call prefetch(state%q_prev)
      end if

   end subroutine crm_state_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_state_finalize(state, MMF_microphysics_scheme)
      type(crm_state_type), intent(inout) :: state
      character(len=*), intent(in) :: MMF_microphysics_scheme    ! CRM microphysics scheme

      if (allocated(state%u_wind))      deallocate(state%u_wind)
      if (allocated(state%v_wind))      deallocate(state%v_wind)
      if (allocated(state%w_wind))      deallocate(state%w_wind)
      if (allocated(state%temperature)) deallocate(state%temperature)
      if (allocated(state%qv))          deallocate(state%qv)

      if (allocated(state%qp)) deallocate(state%qp)
      if (allocated(state%qn)) deallocate(state%qn)

      if (allocated(state%qc)) deallocate(state%qc)
      if (allocated(state%qi)) deallocate(state%qi)
      if (allocated(state%qr)) deallocate(state%qr)
      if (allocated(state%nc)) deallocate(state%nc)
      if (allocated(state%ni)) deallocate(state%ni)
      if (allocated(state%nr)) deallocate(state%nr)

      if (allocated(state%qm)) deallocate(state%qm)
      if (allocated(state%bm)) deallocate(state%bm)
      if (allocated(state%t_prev)) deallocate(state%t_prev)
      if (allocated(state%q_prev)) deallocate(state%q_prev)
      
   end subroutine crm_state_finalize
end module crm_state_module
