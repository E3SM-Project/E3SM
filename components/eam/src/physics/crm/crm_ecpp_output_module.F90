module crm_ecpp_output_module

   use params_kind, only: crm_rknd

   implicit none

   private

   public crm_ecpp_output_type
   public crm_ecpp_output_initialize
   public crm_ecpp_output_finalize
   public crm_ecpp_output_copy

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
   end type crm_ecpp_output_type
   !------------------------------------------------------------------------------------------------

contains

   !------------------------------------------------------------------------------------------------
   subroutine crm_ecpp_output_initialize(output, ncol, nlev)
      use ecppvars, only: NCLASS_CL, ncls_ecpp_in, NCLASS_PR
      type(crm_ecpp_output_type), intent(inout) :: output
      integer, intent(in) :: ncol, nlev
      if (.not.allocated(output%acen            )) allocate(output%acen            (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%acen_tf         )) allocate(output%acen_tf         (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%rhcen           )) allocate(output%rhcen           (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%qcloudcen       )) allocate(output%qcloudcen       (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%qicecen         )) allocate(output%qicecen         (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%qlsinkcen       )) allocate(output%qlsinkcen       (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%precrcen        )) allocate(output%precrcen        (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%precsolidcen    )) allocate(output%precsolidcen    (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%qlsink_afcen    )) allocate(output%qlsink_afcen    (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%qlsink_bfcen    )) allocate(output%qlsink_bfcen    (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%qlsink_avgcen   )) allocate(output%qlsink_avgcen   (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%praincen        )) allocate(output%praincen        (ncol,nlev,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%wwqui_cen       )) allocate(output%wwqui_cen       (ncol,nlev) )
      if (.not.allocated(output%wwqui_cloudy_cen)) allocate(output%wwqui_cloudy_cen(ncol,nlev) )
      if (.not.allocated(output%abnd            )) allocate(output%abnd            (ncol,nlev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%abnd_tf         )) allocate(output%abnd_tf         (ncol,nlev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%massflxbnd      )) allocate(output%massflxbnd      (ncol,nlev+1,NCLASS_CL,ncls_ecpp_in,NCLASS_PR) )
      if (.not.allocated(output%wupthresh_bnd   )) allocate(output%wupthresh_bnd   (ncol,nlev+1) )
      if (.not.allocated(output%wdownthresh_bnd )) allocate(output%wdownthresh_bnd (ncol,nlev+1) )
      if (.not.allocated(output%wwqui_bnd       )) allocate(output%wwqui_bnd       (ncol,nlev+1) )
      if (.not.allocated(output%wwqui_cloudy_bnd)) allocate(output%wwqui_cloudy_bnd(ncol,nlev+1) )
   end subroutine crm_ecpp_output_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_ecpp_output_finalize(output)
      type(crm_ecpp_output_type), intent(inout) :: output
      if (allocated(output%acen            )) deallocate(output%acen            )
      if (allocated(output%acen_tf         )) deallocate(output%acen_tf         )
      if (allocated(output%rhcen           )) deallocate(output%rhcen           )
      if (allocated(output%qcloudcen       )) deallocate(output%qcloudcen       )
      if (allocated(output%qicecen         )) deallocate(output%qicecen         )
      if (allocated(output%qlsinkcen       )) deallocate(output%qlsinkcen       )
      if (allocated(output%precrcen        )) deallocate(output%precrcen        )
      if (allocated(output%precsolidcen    )) deallocate(output%precsolidcen    )
      if (allocated(output%qlsink_afcen    )) deallocate(output%qlsink_afcen    )
      if (allocated(output%qlsink_bfcen    )) deallocate(output%qlsink_bfcen    )
      if (allocated(output%qlsink_avgcen   )) deallocate(output%qlsink_avgcen   )
      if (allocated(output%praincen        )) deallocate(output%praincen        )
      if (allocated(output%abnd            )) deallocate(output%abnd            )
      if (allocated(output%abnd_tf         )) deallocate(output%abnd_tf         )
      if (allocated(output%massflxbnd      )) deallocate(output%massflxbnd      )
      if (allocated(output%wupthresh_bnd   )) deallocate(output%wupthresh_bnd   )
      if (allocated(output%wdownthresh_bnd )) deallocate(output%wdownthresh_bnd )
      if (allocated(output%wwqui_cen       )) deallocate(output%wwqui_cen       )
      if (allocated(output%wwqui_cloudy_cen)) deallocate(output%wwqui_cloudy_cen)
      if (allocated(output%wwqui_bnd       )) deallocate(output%wwqui_bnd       )
      if (allocated(output%wwqui_cloudy_bnd)) deallocate(output%wwqui_cloudy_bnd)
   end subroutine crm_ecpp_output_finalize
   !------------------------------------------------------------------------------------------------
   subroutine crm_ecpp_output_copy(output, output_copy, col_beg, col_end)
      type(crm_ecpp_output_type), intent(in   ) :: output
      type(crm_ecpp_output_type), intent(inout) :: output_copy
      integer,                    intent(in   ) :: col_beg
      integer,                    intent(in   ) :: col_end
      integer :: ncol_copy
      ncol_copy = col_end - col_beg + 1
      output_copy%abnd            (1:ncol_copy,:,:,:,:) = output%abnd            (col_beg:col_end,:,:,:,:)
      output_copy%abnd_tf         (1:ncol_copy,:,:,:,:) = output%abnd_tf         (col_beg:col_end,:,:,:,:)
      output_copy%massflxbnd      (1:ncol_copy,:,:,:,:) = output%massflxbnd      (col_beg:col_end,:,:,:,:)
      output_copy%acen            (1:ncol_copy,:,:,:,:) = output%acen            (col_beg:col_end,:,:,:,:)
      output_copy%acen_tf         (1:ncol_copy,:,:,:,:) = output%acen_tf         (col_beg:col_end,:,:,:,:)
      output_copy%rhcen           (1:ncol_copy,:,:,:,:) = output%rhcen           (col_beg:col_end,:,:,:,:)
      output_copy%qcloudcen       (1:ncol_copy,:,:,:,:) = output%qcloudcen       (col_beg:col_end,:,:,:,:)
      output_copy%qicecen         (1:ncol_copy,:,:,:,:) = output%qicecen         (col_beg:col_end,:,:,:,:)
      output_copy%qlsinkcen       (1:ncol_copy,:,:,:,:) = output%qlsinkcen       (col_beg:col_end,:,:,:,:)
      output_copy%precrcen        (1:ncol_copy,:,:,:,:) = output%precrcen        (col_beg:col_end,:,:,:,:)
      output_copy%precsolidcen    (1:ncol_copy,:,:,:,:) = output%precsolidcen    (col_beg:col_end,:,:,:,:)
      output_copy%qlsink_afcen    (1:ncol_copy,:,:,:,:) = output%qlsink_afcen    (col_beg:col_end,:,:,:,:)
      output_copy%qlsink_bfcen    (1:ncol_copy,:,:,:,:) = output%qlsink_bfcen    (col_beg:col_end,:,:,:,:)
      output_copy%qlsink_avgcen   (1:ncol_copy,:,:,:,:) = output%qlsink_avgcen   (col_beg:col_end,:,:,:,:)
      output_copy%praincen        (1:ncol_copy,:,:,:,:) = output%praincen        (col_beg:col_end,:,:,:,:)
      output_copy%wupthresh_bnd   (1:ncol_copy,:)       = output%wupthresh_bnd   (col_beg:col_end,:)
      output_copy%wdownthresh_bnd (1:ncol_copy,:)       = output%wdownthresh_bnd (col_beg:col_end,:)
      output_copy%wwqui_cen       (1:ncol_copy,:)       = output%wwqui_cen       (col_beg:col_end,:)
      output_copy%wwqui_bnd       (1:ncol_copy,:)       = output%wwqui_bnd       (col_beg:col_end,:)
      output_copy%wwqui_cloudy_cen(1:ncol_copy,:)       = output%wwqui_cloudy_cen(col_beg:col_end,:)
      output_copy%wwqui_cloudy_bnd(1:ncol_copy,:)       = output%wwqui_cloudy_bnd(col_beg:col_end,:)
   end subroutine crm_ecpp_output_copy
   !------------------------------------------------------------------------------------------------

end module crm_ecpp_output_module
