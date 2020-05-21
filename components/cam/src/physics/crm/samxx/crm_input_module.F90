module crm_input_module
   use params, only: crm_rknd, crm_iknd
   implicit none
   private
   public crm_input_type
   type crm_input_type

      real(crm_rknd), pointer, contiguous :: zmid(:,:)           ! Global grid height (m)
      real(crm_rknd), pointer, contiguous :: zint(:,:)           ! Global grid interface height (m)
      real(crm_rknd), pointer, contiguous :: tl(:,:)             ! Global grid temperature (K)
      real(crm_rknd), pointer, contiguous :: ql(:,:)             ! Global grid water vapor (g/g)
      real(crm_rknd), pointer, contiguous :: qccl(:,:)           ! Global grid cloud liquid water (g/g)
      real(crm_rknd), pointer, contiguous :: qiil(:,:)           ! Global grid cloud ice (g/g)
      real(crm_rknd), pointer, contiguous :: ps(:)               ! Global grid surface pressure (Pa)
      real(crm_rknd), pointer, contiguous :: pmid(:,:)           ! Global grid pressure (Pa)
      real(crm_rknd), pointer, contiguous :: pint(:,:)           ! Global grid pressure (Pa)
      real(crm_rknd), pointer, contiguous :: pdel(:,:)           ! Layer's pressure thickness (Pa)
      real(crm_rknd), pointer, contiguous :: phis(:)             ! Global grid surface geopotential (m2/s2)
      real(crm_rknd), pointer, contiguous :: ul(:,:)             ! Global grid u (m/s)
      real(crm_rknd), pointer, contiguous :: vl(:,:)             ! Global grid v (m/s)
      real(crm_rknd), pointer, contiguous :: ocnfrac(:)          ! area fraction of the ocean
      real(crm_rknd), pointer, contiguous :: tau00  (:)          ! large-scale surface stress (N/m2)
      real(crm_rknd), pointer, contiguous :: wndls  (:)          ! large-scale surface wind (m/s)
      real(crm_rknd), pointer, contiguous :: bflxls (:)          ! large-scale surface buoyancy flux (K m/s)
      real(crm_rknd), pointer, contiguous :: fluxu00(:)          ! surface momenent fluxes [N/m2]
      real(crm_rknd), pointer, contiguous :: fluxv00(:)          ! surface momenent fluxes [N/m2]
      real(crm_rknd), pointer, contiguous :: fluxt00(:)          ! surface sensible heat fluxes [K Kg/ (m2 s)]
      real(crm_rknd), pointer, contiguous :: fluxq00(:)          ! surface latent heat fluxes [ kg/(m2 s)]

   contains
      procedure, public :: initialize=>crm_input_initialize
      procedure, public :: finalize=>crm_input_finalize
   end type crm_input_type
   !------------------------------------------------------------------------------------------------

contains
   !------------------------------------------------------------------------------------------------
   ! Type-bound procedures for crm_input_type
   subroutine crm_input_initialize(this, nlev, ncrms)
      class(crm_input_type), intent(inout) :: this
      integer(crm_iknd), intent(in) :: nlev, ncrms
      
      allocate(this%zmid   (ncrms,nlev  ))
      allocate(this%zint   (ncrms,nlev+1))
      allocate(this%tl     (ncrms,nlev  ))
      allocate(this%ql     (ncrms,nlev  ))
      allocate(this%qccl   (ncrms,nlev  ))
      allocate(this%qiil   (ncrms,nlev  ))
      allocate(this%ps     (ncrms       ))
      allocate(this%pmid   (ncrms,nlev  ))
      allocate(this%pint   (ncrms,nlev+1))
      allocate(this%pdel   (ncrms,nlev  ))
      allocate(this%phis   (ncrms       ))
      allocate(this%ul     (ncrms,nlev  ))
      allocate(this%vl     (ncrms,nlev  ))
      allocate(this%ocnfrac(ncrms       ))
      allocate(this%tau00  (ncrms       ))
      allocate(this%wndls  (ncrms       ))
      allocate(this%bflxls (ncrms       ))
      allocate(this%fluxu00(ncrms       ))
      allocate(this%fluxv00(ncrms       ))
      allocate(this%fluxt00(ncrms       ))
      allocate(this%fluxq00(ncrms       ))

      ! Initialize
      this%zmid = 0
      this%zint = 0
      this%tl = 0
      this%ql = 0
      this%qccl = 0
      this%qiil = 0
      this%ps = 0
      this%pmid = 0
      this%pint = 0
      this%pdel = 0
      this%phis = 0
      this%ul = 0
      this%vl = 0
      this%ocnfrac = 0
      this%tau00   = 0
      this%wndls   = 0
      this%bflxls  = 0
      this%fluxu00 = 0
      this%fluxv00 = 0
      this%fluxt00 = 0
      this%fluxq00 = 0
   end subroutine crm_input_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_input_finalize(this)
      class(crm_input_type), intent(inout) :: this

      deallocate(this%zmid)
      deallocate(this%zint)
      deallocate(this%tl)
      deallocate(this%ql)
      deallocate(this%qccl)
      deallocate(this%qiil)
      deallocate(this%ps)
      deallocate(this%pmid)
      deallocate(this%pint)
      deallocate(this%pdel)
      deallocate(this%phis)
      deallocate(this%ul)
      deallocate(this%vl)
      deallocate(this%ocnfrac)
      deallocate(this%tau00)
      deallocate(this%wndls)
      deallocate(this%bflxls)
      deallocate(this%fluxu00)
      deallocate(this%fluxv00)
      deallocate(this%fluxt00)
      deallocate(this%fluxq00)

   end subroutine crm_input_finalize 

end module crm_input_module
