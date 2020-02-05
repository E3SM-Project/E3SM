module crm_input_module
   use params, only: crm_rknd
#ifdef MODAL_AERO
   use modal_aero_data, only: ntot_amode
#endif
   use openacc_utils
   implicit none
   private
   public crm_input_type
   type crm_input_type

      real(crm_rknd), allocatable :: zmid(:,:)           ! Global grid height (m)
      real(crm_rknd), allocatable :: zint(:,:)           ! Global grid interface height (m)
      real(crm_rknd), allocatable :: tl(:,:)             ! Global grid temperature (K)
      real(crm_rknd), allocatable :: ql(:,:)             ! Global grid water vapor (g/g)
      real(crm_rknd), allocatable :: qccl(:,:)           ! Global grid cloud liquid water (g/g)
      real(crm_rknd), allocatable :: qiil(:,:)           ! Global grid cloud ice (g/g)
      real(crm_rknd), allocatable :: ps(:)               ! Global grid surface pressure (Pa)
      real(crm_rknd), allocatable :: pmid(:,:)           ! Global grid pressure (Pa)
      real(crm_rknd), allocatable :: pint(:,:)           ! Global grid pressure (Pa)
      real(crm_rknd), allocatable :: pdel(:,:)           ! Layer's pressure thickness (Pa)
      real(crm_rknd), allocatable :: phis(:)             ! Global grid surface geopotential (m2/s2)
      real(crm_rknd), allocatable :: ul(:,:)             ! Global grid u (m/s)
      real(crm_rknd), allocatable :: vl(:,:)             ! Global grid v (m/s)
      real(crm_rknd), allocatable :: ocnfrac(:)          ! area fraction of the ocean
      real(crm_rknd), allocatable :: tau00  (:)          ! large-scale surface stress (N/m2)
      real(crm_rknd), allocatable :: wndls  (:)          ! large-scale surface wind (m/s)
      real(crm_rknd), allocatable :: bflxls (:)          ! large-scale surface buoyancy flux (K m/s)
      real(crm_rknd), allocatable :: fluxu00(:)          ! surface momenent fluxes [N/m2]
      real(crm_rknd), allocatable :: fluxv00(:)          ! surface momenent fluxes [N/m2]
      real(crm_rknd), allocatable :: fluxt00(:)          ! surface sensible heat fluxes [K Kg/ (m2 s)]
      real(crm_rknd), allocatable :: fluxq00(:)          ! surface latent heat fluxes [ kg/(m2 s)]

#if defined( m2005 ) && defined( MODAL_AERO )
      real(crm_rknd), allocatable :: naermod (:,:,:)     ! Aerosol number concentration [/m3]
      real(crm_rknd), allocatable :: vaerosol(:,:,:)     ! aerosol volume concentration [m3/m3]
      real(crm_rknd), allocatable :: hygro   (:,:,:)     ! hygroscopicity of aerosol mode 
#endif

#if defined( MMF_ESMT )
      real(crm_rknd), allocatable :: ul_esmt(:,:)        ! input u for ESMT
      real(crm_rknd), allocatable :: vl_esmt(:,:)        ! input v for ESMT
#endif

   contains
      procedure, public :: initialize=>crm_input_initialize
      procedure, public :: finalize=>crm_input_finalize
   end type crm_input_type
   !------------------------------------------------------------------------------------------------

contains
   !------------------------------------------------------------------------------------------------
   ! Type-bound procedures for crm_input_type
   subroutine crm_input_initialize(this, ncrms, nlev)
      class(crm_input_type), intent(inout) :: this
      integer, intent(in) :: ncrms, nlev
      
      if (.not. allocated(this%zmid))     allocate(this%zmid(ncrms,nlev))
      if (.not. allocated(this%zint))     allocate(this%zint(ncrms,nlev+1))
      if (.not. allocated(this%tl))       allocate(this%tl(ncrms,nlev))
      if (.not. allocated(this%ql))       allocate(this%ql(ncrms,nlev))
      if (.not. allocated(this%qccl))     allocate(this%qccl(ncrms,nlev))
      if (.not. allocated(this%qiil))     allocate(this%qiil(ncrms,nlev))
      if (.not. allocated(this%ps))       allocate(this%ps(ncrms))
      if (.not. allocated(this%pmid))     allocate(this%pmid(ncrms,nlev))
      if (.not. allocated(this%pint))     allocate(this%pint(ncrms,nlev+1))
      if (.not. allocated(this%pdel))     allocate(this%pdel(ncrms,nlev))
      if (.not. allocated(this%phis))     allocate(this%phis(ncrms))
      if (.not. allocated(this%ul))       allocate(this%ul(ncrms,nlev))
      if (.not. allocated(this%vl))       allocate(this%vl(ncrms,nlev))
      if (.not. allocated(this%ocnfrac))  allocate(this%ocnfrac(ncrms))
      if (.not. allocated(this%tau00))    allocate(this%tau00(ncrms))
      if (.not. allocated(this%wndls))    allocate(this%wndls(ncrms))
      if (.not. allocated(this%bflxls))   allocate(this%bflxls(ncrms))
      if (.not. allocated(this%fluxu00))  allocate(this%fluxu00(ncrms))
      if (.not. allocated(this%fluxv00))  allocate(this%fluxv00(ncrms))
      if (.not. allocated(this%fluxt00))  allocate(this%fluxt00(ncrms))
      if (.not. allocated(this%fluxq00))  allocate(this%fluxq00(ncrms))

      call prefetch(this%zmid)
      call prefetch(this%zint)
      call prefetch(this%tl)
      call prefetch(this%ql)
      call prefetch(this%qccl)
      call prefetch(this%qiil)
      call prefetch(this%ps)
      call prefetch(this%pmid)
      call prefetch(this%pint)
      call prefetch(this%pdel)
      call prefetch(this%phis)
      call prefetch(this%ul)
      call prefetch(this%vl)
      call prefetch(this%ocnfrac)
      call prefetch(this%tau00)
      call prefetch(this%wndls)
      call prefetch(this%bflxls)
      call prefetch(this%fluxu00)
      call prefetch(this%fluxv00)
      call prefetch(this%fluxt00)
      call prefetch(this%fluxq00)

#if defined( m2005 ) && defined( MODAL_AERO )
      if (.not. allocated(this%naermod))  allocate(this%naermod(ncrms,nlev,ntot_amode))
      if (.not. allocated(this%vaerosol)) allocate(this%vaerosol(ncrms,nlev,ntot_amode))
      if (.not. allocated(this%hygro))    allocate(this%hygro(ncrms,nlev,ntot_amode))
      call prefetch(this%naermod)
      call prefetch(this%vaerosol)
      call prefetch(this%hygro)
#endif

#if defined(MMF_ESMT)
      if (.not. allocated(this%ul_esmt))  allocate(this%ul_esmt(ncrms,nlev))
      if (.not. allocated(this%vl_esmt))  allocate(this%vl_esmt(ncrms,nlev))
#endif

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
#if defined( m2005 ) && defined( MODAL_AERO )
      this%naermod  = 0
      this%vaerosol = 0
      this%hygro    = 0
#endif
#if defined( MMF_ESMT )
      this%ul_esmt = 0
      this%vl_esmt = 0
#endif

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

#if defined( m2005 ) && defined( MODAL_AERO )
      deallocate(this%naermod)
      deallocate(this%vaerosol)
      deallocate(this%hygro)
#endif

#if defined(MMF_ESMT)
      deallocate(this%ul_esmt)
      deallocate(this%vl_esmt)
#endif

   end subroutine crm_input_finalize 

end module crm_input_module
