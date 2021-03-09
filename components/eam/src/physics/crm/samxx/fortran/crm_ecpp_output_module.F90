module crm_ecpp_output_module

   use params, only: crm_rknd

   implicit none

   private

   public crm_ecpp_output_type

   !------------------------------------------------------------------------------------------------
   type crm_ecpp_output_type
      ! Purpose: Derived type to encapsulate the CRM output relevant for ECPP
      real(crm_rknd), allocatable :: abnd         (:,:,:,:,:)  ! cloud fraction for each sub-sub class for full time period
      real(crm_rknd), allocatable :: abnd_tf      (:,:,:,:,:)  ! cloud fraction for end-portion of time period
      real(crm_rknd), allocatable :: massflxbnd   (:,:,:,:,:)  ! sub-class vertical mass flux (kg/m2/s) at layer bottom boundary.
      real(crm_rknd), allocatable :: acen         (:,:,:,:,:)  ! cloud fraction for each sub-sub class for full time period
      real(crm_rknd), allocatable :: acen_tf      (:,:,:,:,:)  ! cloud fraction for end-portion of time period
      real(crm_rknd), allocatable :: rhcen        (:,:,:,:,:)  ! relative humidity (0-1)
      real(crm_rknd), allocatable :: qcloudcen    (:,:,:,:,:)  ! cloud water (kg/kg)
      real(crm_rknd), allocatable :: qicecen      (:,:,:,:,:)  ! cloud ice (kg/kg)
      real(crm_rknd), allocatable :: qlsinkcen    (:,:,:,:,:)  ! cloud water loss rate from precipitation (/s??)
      real(crm_rknd), allocatable :: precrcen     (:,:,:,:,:)  ! liquid (rain) precipitation rate (kg/m2/s)
      real(crm_rknd), allocatable :: precsolidcen (:,:,:,:,:)  ! solid (rain) precipitation rate (kg/m2/s)
      real(crm_rknd), allocatable :: qlsink_afcen (:,:,:,:,:)  ! cld water loss rate from precip calc from cloud water after precipitating (/s)
      real(crm_rknd), allocatable :: qlsink_bfcen (:,:,:,:,:)  ! cld water loss rate from precip calc from cloud water before precipitating (/s)
      real(crm_rknd), allocatable :: qlsink_avgcen(:,:,:,:,:)  ! cld water loss rate from precip calc from praincen and qlcoudcen averaged over ntavg1_ss time step (/s??)
      real(crm_rknd), allocatable :: praincen     (:,:,:,:,:)  ! cld water loss rate from precip (kg/kg/s)
      real(crm_rknd), allocatable :: wupthresh_bnd      (:,:)  ! vert velocity threshold for updraft (m/s)
      real(crm_rknd), allocatable :: wdownthresh_bnd    (:,:)  ! vert velocity threshold for downdraft (m/s)
      real(crm_rknd), allocatable :: wwqui_cen          (:,:)  ! vert velocity variance in quiescent class (m2/s2) at layer center
      real(crm_rknd), allocatable :: wwqui_bnd          (:,:)  ! vert velocity variance in quiescent class (m2/s2) at layer boundary
      real(crm_rknd), allocatable :: wwqui_cloudy_cen   (:,:)  ! vert velocity variance in quiescent and cloudy class (m2/s2) at layer center
      real(crm_rknd), allocatable :: wwqui_cloudy_bnd   (:,:)  ! vert velocity variance in quiescent and cloudy class (m2/s2) at layer boundary
#if defined( ECPP )
   contains
      procedure, public :: initialize=>crm_ecpp_initialize
      procedure, public :: finalize=>crm_ecpp_finalize
#endif /* ECPP */
   end type crm_ecpp_output_type
   !------------------------------------------------------------------------------------------------

contains

#if defined( ECPP )

   !------------------------------------------------------------------------------------------------
   subroutine crm_ecpp_initialize(this, ncol, nlev)
      use ecppvars, only: NCLASS_CL, ncls_ecpp_in, NCLASS_PR
      class(crm_ecpp_output_type), intent(inout) :: this
      integer, intent(in) :: ncol, nlev
      if (.not.allocated(this%acen            )) allocate(this%acen            (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%acen_tf         )) allocate(this%acen_tf         (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%rhcen           )) allocate(this%rhcen           (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%qcloudcen       )) allocate(this%qcloudcen       (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%qicecen         )) allocate(this%qicecen         (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%qlsinkcen       )) allocate(this%qlsinkcen       (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%precrcen        )) allocate(this%precrcen        (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%precsolidcen    )) allocate(this%precsolidcen    (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%qlsink_afcen    )) allocate(this%qlsink_afcen    (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%qlsink_bfcen    )) allocate(this%qlsink_bfcen    (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%qlsink_avgcen   )) allocate(this%qlsink_avgcen   (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%praincen        )) allocate(this%praincen        (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%wwqui_cen       )) allocate(this%wwqui_cen       (ncol,nlev) )
      if (.not.allocated(this%wwqui_cloudy_cen)) allocate(this%wwqui_cloudy_cen(ncol,nlev) )
      if (.not.allocated(this%abnd            )) allocate(this%abnd            (ncol,nlev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%abnd_tf         )) allocate(this%abnd_tf         (ncol,nlev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%massflxbnd      )) allocate(this%massflxbnd      (ncol,nlev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(this%wupthresh_bnd   )) allocate(this%wupthresh_bnd   (ncol,nlev+1) )
      if (.not.allocated(this%wdownthresh_bnd )) allocate(this%wdownthresh_bnd (ncol,nlev+1) )
      if (.not.allocated(this%wwqui_bnd       )) allocate(this%wwqui_bnd       (ncol,nlev+1) )
      if (.not.allocated(this%wwqui_cloudy_bnd)) allocate(this%wwqui_cloudy_bnd(ncol,nlev+1) )
   end subroutine crm_ecpp_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_ecpp_finalize(this)
      class(crm_ecpp_output_type), intent(inout) :: this
      deallocate(this%acen            )
      deallocate(this%acen_tf         )
      deallocate(this%rhcen           )
      deallocate(this%qcloudcen       )
      deallocate(this%qicecen         )
      deallocate(this%qlsinkcen       )
      deallocate(this%precrcen        )
      deallocate(this%precsolidcen    )
      deallocate(this%qlsink_afcen    )
      deallocate(this%qlsink_bfcen    )
      deallocate(this%qlsink_avgcen   )
      deallocate(this%praincen        )
      deallocate(this%abnd            )
      deallocate(this%abnd_tf         )
      deallocate(this%massflxbnd      )
      deallocate(this%wupthresh_bnd   )
      deallocate(this%wdownthresh_bnd )
      deallocate(this%wwqui_cen       )
      deallocate(this%wwqui_cloudy_cen)
      deallocate(this%wwqui_bnd       )
      deallocate(this%wwqui_cloudy_bnd)
   end subroutine crm_ecpp_finalize
   !------------------------------------------------------------------------------------------------

#endif /* ECPP */

end module crm_ecpp_output_module