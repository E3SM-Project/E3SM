module crm_rad_module
   ! Module to encapsulate data and methods specific to radiation. This exists
   ! as a separate module with separate derived types because radiation is
   ! handled in a special way when interacting between the GCM and CRM.
   ! The radiation calculations may be on a separate grid than the CRM, and the
   ! quantities may be aggregated and updated between calls in a different way than
   ! the other diagnostic quantities. This module should also contain methods to
   ! update the radiation in the future, should we choose to put the radiation
   ! calculations on the CRM.
   use params_kind, only: crm_rknd, r8
   use openacc_utils

   implicit none

   public crm_rad_type

   type crm_rad_type
      ! Radiative heating
      real(crm_rknd), allocatable :: qrad(:,:,:,:)

      ! Quantities used by the radiation code. Note that these are strange in that they are 
      ! time-averages, but spatially-resolved.
      real(crm_rknd), allocatable :: temperature(:,:,:,:) ! rad temperature
      real(crm_rknd), allocatable :: qv (:,:,:,:) ! rad vapor
      real(crm_rknd), allocatable :: qc (:,:,:,:) ! rad cloud water
      real(crm_rknd), allocatable :: qi (:,:,:,:) ! rad cloud ice
      real(crm_rknd), allocatable :: cld(:,:,:,:) ! rad cloud fraction

      ! Only relevant when using 2-moment microphysics (ex. P3)
      real(crm_rknd), allocatable :: nc(:,:,:,:) ! rad cloud droplet number (#/kg)
      real(crm_rknd), allocatable :: ni(:,:,:,:) ! rad cloud ice crystal number (#/kg)
      real(crm_rknd), allocatable :: qs(:,:,:,:) ! rad cloud snow (kg/kg)
      real(crm_rknd), allocatable :: ns(:,:,:,:) ! rad cloud snow crystal number (#/kg)
   end type crm_rad_type

contains

   !------------------------------------------------------------------------------------------------
   ! Type-bound procedures for crm_rad_type
   subroutine crm_rad_initialize(rad, ncrms, crm_nx_rad, crm_ny_rad, crm_nz, MMF_microphysics_scheme)
      type(crm_rad_type),  intent(inout) :: rad
      integer,             intent(in   ) :: ncrms, crm_nx_rad, crm_ny_rad, crm_nz
      character(len=*),    intent(in   ) :: MMF_microphysics_scheme    ! CRM microphysics scheme

      if (.not. allocated(rad%qrad))        allocate(rad%qrad       (ncrms, crm_nx_rad, crm_ny_rad, crm_nz))
      if (.not. allocated(rad%temperature)) allocate(rad%temperature(ncrms, crm_nx_rad, crm_ny_rad, crm_nz))
      if (.not. allocated(rad%qv))          allocate(rad%qv         (ncrms, crm_nx_rad, crm_ny_rad, crm_nz))
      if (.not. allocated(rad%qc))          allocate(rad%qc         (ncrms, crm_nx_rad, crm_ny_rad, crm_nz))
      if (.not. allocated(rad%qi))          allocate(rad%qi         (ncrms, crm_nx_rad, crm_ny_rad, crm_nz))
      if (.not. allocated(rad%cld))         allocate(rad%cld        (ncrms, crm_nx_rad, crm_ny_rad, crm_nz))
      if (.not. allocated(rad%nc))          allocate(rad%nc         (ncrms, crm_nx_rad, crm_ny_rad, crm_nz))
      if (.not. allocated(rad%ni))          allocate(rad%ni         (ncrms, crm_nx_rad, crm_ny_rad, crm_nz))
      if (.not. allocated(rad%qs))          allocate(rad%qs         (ncrms, crm_nx_rad, crm_ny_rad, crm_nz))
      if (.not. allocated(rad%ns))          allocate(rad%ns         (ncrms, crm_nx_rad, crm_ny_rad, crm_nz))

      call prefetch(rad%qrad)
      call prefetch(rad%temperature)
      call prefetch(rad%qv)
      call prefetch(rad%qc)
      call prefetch(rad%qi)
      call prefetch(rad%cld)
      call prefetch(rad%nc)
      call prefetch(rad%ni)
      call prefetch(rad%qs)
      call prefetch(rad%ns)

      rad%qrad        = 0
      rad%temperature = 0
      rad%qv          = 0
      rad%qc          = 0
      rad%qi          = 0
      rad%cld         = 0
      rad%nc          = 0
      rad%ni          = 0
      rad%qs          = 0
      rad%ns          = 0

   end subroutine crm_rad_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_rad_finalize(rad, MMF_microphysics_scheme)
      type(crm_rad_type), intent(inout) :: rad
      character(len=*), intent(in) :: MMF_microphysics_scheme    ! CRM microphysics scheme

      if (allocated(rad%qrad))        deallocate(rad%qrad)
      if (allocated(rad%temperature)) deallocate(rad%temperature)
      if (allocated(rad%qv))          deallocate(rad%qv)
      if (allocated(rad%qc))          deallocate(rad%qc)
      if (allocated(rad%qi))          deallocate(rad%qi)
      if (allocated(rad%cld))         deallocate(rad%cld)
      if (allocated(rad%nc))       deallocate(rad%nc)
      if (allocated(rad%ni))       deallocate(rad%ni)
      if (allocated(rad%qs))       deallocate(rad%qs)
      if (allocated(rad%ns))       deallocate(rad%ns)

   end subroutine crm_rad_finalize
   !------------------------------------------------------------------------------------------------

end module crm_rad_module
