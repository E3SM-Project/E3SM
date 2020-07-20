
module crm_input_module
   use params, only: crm_rknd
#ifdef MODAL_AERO
   use modal_aero_data, only: ntot_amode
#endif
#if defined(_OPENACC)
   use openacc_utils
#endif

   implicit none
   public 
   real(crm_rknd), allocatable :: crm_input_zmid(:,:)           ! Global grid height (m)
   real(crm_rknd), allocatable :: crm_input_zint(:,:)           ! Global grid interface height (m)
   real(crm_rknd), allocatable :: crm_input_tl(:,:)             ! Global grid temperature (K)
   real(crm_rknd), allocatable :: crm_input_ql(:,:)             ! Global grid water vapor (g/g)
   real(crm_rknd), allocatable :: crm_input_qccl(:,:)           ! Global grid cloud liquid water (g/g)
   real(crm_rknd), allocatable :: crm_input_qiil(:,:)           ! Global grid cloud ice (g/g)
   real(crm_rknd), allocatable :: crm_input_ps(:)               ! Global grid surface pressure (Pa)
   real(crm_rknd), allocatable :: crm_input_pmid(:,:)           ! Global grid pressure (Pa)
   real(crm_rknd), allocatable :: crm_input_pint(:,:)           ! Global grid pressure (Pa)
   real(crm_rknd), allocatable :: crm_input_pdel(:,:)           ! Layer's pressure thickness (Pa)
   real(crm_rknd), allocatable :: crm_input_phis(:)             ! Global grid surface geopotential (m2/s2)
   real(crm_rknd), allocatable :: crm_input_ul(:,:)             ! Global grid u (m/s)
   real(crm_rknd), allocatable :: crm_input_vl(:,:)             ! Global grid v (m/s)
   real(crm_rknd), allocatable :: crm_input_ocnfrac(:)          ! area fraction of the ocean
   real(crm_rknd), allocatable :: crm_input_tau00(:)            ! large-scale surface stress (N/m2)
   real(crm_rknd), allocatable :: crm_input_wndls(:)            ! large-scale surface wind (m/s)
   real(crm_rknd), allocatable :: crm_input_bflxls(:)           ! large-scale surface buoyancy flux (K m/s)
   real(crm_rknd), allocatable :: crm_input_fluxu00(:)          ! surface momenent fluxes [N/m2]
   real(crm_rknd), allocatable :: crm_input_fluxv00(:)          ! surface momenent fluxes [N/m2]
   real(crm_rknd), allocatable :: crm_input_fluxt00(:)          ! surface sensible heat fluxes [K Kg/ (m2 s)]
   real(crm_rknd), allocatable :: crm_input_fluxq00(:)          ! surface latent heat fluxes [ kg/(m2 s)]

#if defined( m2005 ) && defined( MODAL_AERO )
   real(crm_rknd), allocatable :: crm_input_naermod(:,:,:)      ! Aerosol number concentration [/m3]
   real(crm_rknd), allocatable :: crm_input_vaerosol(:,:,:)     ! aerosol volume concentration [m3/m3]
   real(crm_rknd), allocatable :: crm_input_hygro(:,:,:)        ! hygroscopicity of aerosol mode 
#endif

#if defined( SP_ESMT )
   real(crm_rknd), allocatable :: crm_input_ul_esmt(:,:)        ! input u for ESMT
   real(crm_rknd), allocatable :: crm_input_vl_esmt(:,:)        ! input v for ESMT
#endif

   !------------------------------------------------------------------------------------------------
   public :: crm_input_initialize
   public :: crm_input_finalize

#if defined(_OPENMP)
   public :: update_device_input
#endif
contains
   !------------------------------------------------------------------------------------------------
   ! Type-bound procedures for crm_input_type
   subroutine crm_input_initialize(ncrms, nlev)
      integer, intent(in) :: ncrms, nlev
      integer :: istat, icrm, i, j
      integer, save :: icall
      icall = icall + 1

      if (.not. allocated(crm_input_zmid))     allocate(crm_input_zmid(ncrms,nlev))
      if (.not. allocated(crm_input_zint))     allocate(crm_input_zint(ncrms,nlev+1))
      if (.not. allocated(crm_input_tl))       allocate(crm_input_tl(ncrms,nlev))
      if (.not. allocated(crm_input_ql))       allocate(crm_input_ql(ncrms,nlev))
      if (.not. allocated(crm_input_qccl))     allocate(crm_input_qccl(ncrms,nlev))
      if (.not. allocated(crm_input_qiil))     allocate(crm_input_qiil(ncrms,nlev))
      if (.not. allocated(crm_input_ps))       allocate(crm_input_ps(ncrms))
      if (.not. allocated(crm_input_pmid))     allocate(crm_input_pmid(ncrms,nlev))
      if (.not. allocated(crm_input_pint))     allocate(crm_input_pint(ncrms,nlev+1))
      if (.not. allocated(crm_input_pdel))     allocate(crm_input_pdel(ncrms,nlev))
      if (.not. allocated(crm_input_phis))     allocate(crm_input_phis(ncrms))
      if (.not. allocated(crm_input_ul))       allocate(crm_input_ul(ncrms,nlev))
      if (.not. allocated(crm_input_vl))       allocate(crm_input_vl(ncrms,nlev))
      if (.not. allocated(crm_input_ocnfrac))  allocate(crm_input_ocnfrac(ncrms))
      if (.not. allocated(crm_input_tau00))    allocate(crm_input_tau00(ncrms))
      if (.not. allocated(crm_input_wndls))    allocate(crm_input_wndls(ncrms))
      if (.not. allocated(crm_input_bflxls))   allocate(crm_input_bflxls(ncrms))
      if (.not. allocated(crm_input_fluxu00))  allocate(crm_input_fluxu00(ncrms))
      if (.not. allocated(crm_input_fluxv00))  allocate(crm_input_fluxv00(ncrms))
      if (.not. allocated(crm_input_fluxt00))  allocate(crm_input_fluxt00(ncrms))
      if (.not. allocated(crm_input_fluxq00))  allocate(crm_input_fluxq00(ncrms))

#if defined( m2005 ) && defined( MODAL_AERO )
      if (.not. allocated(crm_input_naermod))  allocate(crm_input_naermod(ncrms,nlev,ntot_amode))
      if (.not. allocated(crm_input_vaerosol)) allocate(crm_input_vaerosol(ncrms,nlev,ntot_amode))
      if (.not. allocated(crm_input_hygro))    allocate(crm_input_hygro(ncrms,nlev,ntot_amode))
#endif

#if defined(SP_ESMT)
      if (.not. allocated(crm_input_ul_esmt))  allocate(crm_input_ul_esmt(ncrms,nlev))
      if (.not. allocated(crm_input_vl_esmt))  allocate(crm_input_vl_esmt(ncrms,nlev))
#endif

#if defined(_OPENACC)
      call prefetch(crm_input_zmid)
      call prefetch(crm_input_zint)
      call prefetch(crm_input_tl)
      call prefetch(crm_input_ql)
      call prefetch(crm_input_qccl)
      call prefetch(crm_input_qiil)
      call prefetch(crm_input_ps)
      call prefetch(crm_input_pmid)
      call prefetch(crm_input_pint)
      call prefetch(crm_input_pdel)
      call prefetch(crm_input_phis)
      call prefetch(crm_input_ul)
      call prefetch(crm_input_vl)
      call prefetch(crm_input_ocnfrac)
      call prefetch(crm_input_tau00)
      call prefetch(crm_input_wndls)
      call prefetch(crm_input_bflxls)
      call prefetch(crm_input_fluxu00)
      call prefetch(crm_input_fluxv00)
      call prefetch(crm_input_fluxt00)
      call prefetch(crm_input_fluxq00)
#if defined( m2005 ) && defined( MODAL_AERO )
      call prefetch(crm_input_naermod)
      call prefetch(crm_input_vaerosol)
      call prefetch(crm_input_hygro)
#endif
#if defined(SP_ESMT)
      call prefetch(crm_input_ul_esmt)
      call prefetch(crm_input_vl_esmt)
#endif

#elif defined(_OPENMP)
      !$omp target enter data map(alloc: crm_input_zmid) 
      !$omp target enter data map(alloc: crm_input_zint) 
      !$omp target enter data map(alloc: crm_input_tl) 
      !$omp target enter data map(alloc: crm_input_ql) 
      !$omp target enter data map(alloc: crm_input_qccl) 
      !$omp target enter data map(alloc: crm_input_qiil) 
      !$omp target enter data map(alloc: crm_input_ps) 
      !$omp target enter data map(alloc: crm_input_pmid) 
      !$omp target enter data map(alloc: crm_input_pint) 
      !$omp target enter data map(alloc: crm_input_pdel) 
      !$omp target enter data map(alloc: crm_input_phis) 
      !$omp target enter data map(alloc: crm_input_ul) 
      !$omp target enter data map(alloc: crm_input_vl) 
      !$omp target enter data map(alloc: crm_input_ocnfrac) 
      !$omp target enter data map(alloc: crm_input_tau00) 
      !$omp target enter data map(alloc: crm_input_wndls) 
      !$omp target enter data map(alloc: crm_input_bflxls) 
      !$omp target enter data map(alloc: crm_input_fluxu00) 
      !$omp target enter data map(alloc: crm_input_fluxv00)
      !$omp target enter data map(alloc: crm_input_fluxt00) 
      !$omp target enter data map(alloc: crm_input_fluxq00) 
#if defined( m2005 ) && defined( MODAL_AERO )
      !$omp target enter data map(alloc: crm_input_naermod)
      !$omp target enter data map(alloc: crm_input_vaerosol)
      !$omp target enter data map(alloc: crm_input_hygro)
#endif
#if defined(SP_ESMT)
      !$omp target enter data map(alloc: crm_input_ul_esmt)
      !$omp target enter data map(alloc: crm_input_vl_esmt)
#endif
#endif

      ! Initialize
      crm_input_zmid = 0
      crm_input_zint = 0
      crm_input_tl = 0
      crm_input_ql = 0
      crm_input_qccl = 0
      crm_input_qiil = 0
      crm_input_ps = 0
      crm_input_pmid = 0
      crm_input_pint = 0
      crm_input_pdel = 0
      crm_input_phis = 0
      crm_input_ul = 0
      crm_input_vl = 0
      crm_input_ocnfrac = 0
      crm_input_tau00   = 0
      crm_input_wndls   = 0
      crm_input_bflxls  = 0
      crm_input_fluxu00 = 0
      crm_input_fluxv00 = 0
      crm_input_fluxt00 = 0
      crm_input_fluxq00 = 0
#if defined( m2005 ) && defined( MODAL_AERO )
      crm_input_naermod  = 0
      crm_input_vaerosol = 0
      crm_input_hygro    = 0
#endif
#if defined( SP_ESMT )
      crm_input_ul_esmt = 0
      crm_input_vl_esmt = 0
#endif
   end subroutine crm_input_initialize

   !\-------------------------
   ! update device data
   !/
#if defined(_OPENMP)
   subroutine update_device_input()
      !$omp target update to (crm_input_zmid)
      !$omp target update to (crm_input_zint)
      !$omp target update to (crm_input_tl)
      !$omp target update to (crm_input_ql)
      !$omp target update to (crm_input_qccl)
      !$omp target update to (crm_input_qiil)
      !$omp target update to (crm_input_ps)
      !$omp target update to (crm_input_pmid)
      !$omp target update to (crm_input_pint)
      !$omp target update to (crm_input_pdel)
      !$omp target update to (crm_input_phis)
      !$omp target update to (crm_input_ul)
      !$omp target update to (crm_input_vl)
      !$omp target update to (crm_input_ocnfrac)
      !$omp target update to (crm_input_tau00)
      !$omp target update to (crm_input_wndls)
      !$omp target update to (crm_input_bflxls)
      !$omp target update to (crm_input_fluxu00)
      !$omp target update to (crm_input_fluxv00)
      !$omp target update to (crm_input_fluxt00)
      !$omp target update to (crm_input_fluxq00)
#if defined( m2005 ) && defined( MODAL_AERO )
      !$omp target update to (crm_input_naermod)
      !$omp target update to (crm_input_vaerosol)
      !$omp target update to (crm_input_hygro)
#endif
#if defined(SP_ESMT)
      !$omp target update to (crm_input_ul_esmt)
      !$omp target update to (crm_input_vl_esmt)
#endif
   end subroutine update_device_input
#endif

   subroutine crm_input_finalize()
#if defined(_OPENMP)
      !$omp target exit data map(delete: crm_input_zmid)
      !$omp target exit data map(delete: crm_input_zint)
      !$omp target exit data map(delete: crm_input_tl)
      !$omp target exit data map(delete: crm_input_ql)
      !$omp target exit data map(delete: crm_input_qccl)
      !$omp target exit data map(delete: crm_input_qiil)
      !$omp target exit data map(delete: crm_input_ps)
      !$omp target exit data map(delete: crm_input_pmid)
      !$omp target exit data map(delete: crm_input_pint)
      !$omp target exit data map(delete: crm_input_pdel)
      !$omp target exit data map(delete: crm_input_phis)
      !$omp target exit data map(delete: crm_input_ul)
      !$omp target exit data map(delete: crm_input_vl)
      !$omp target exit data map(delete: crm_input_ocnfrac)
      !$omp target exit data map(delete: crm_input_tau00)
      !$omp target exit data map(delete: crm_input_wndls)
      !$omp target exit data map(delete: crm_input_bflxls)
      !$omp target exit data map(delete: crm_input_fluxu00)
      !$omp target exit data map(delete: crm_input_fluxv00)
      !$omp target exit data map(delete: crm_input_fluxt00)
      !$omp target exit data map(delete: crm_input_fluxq00)
#if defined( m2005 ) && defined( MODAL_AERO )
      !$omp target exit data map(delete: crm_input_naermod)
      !$omp target exit data map(delete: crm_input_vaerosol)
      !$omp target exit data map(delete: crm_input_hygro)
#endif
#if defined(SP_ESMT)
      !$omp target exit data map(delete: crm_input_ul_esmt)
      !$omp target exit data map(delete: crm_input_vl_esmt)
#endif
#endif
      deallocate(crm_input_zmid)
      deallocate(crm_input_zint)
      deallocate(crm_input_tl)
      deallocate(crm_input_ql)
      deallocate(crm_input_qccl)
      deallocate(crm_input_qiil)
      deallocate(crm_input_ps)
      deallocate(crm_input_pmid)
      deallocate(crm_input_pint)
      deallocate(crm_input_pdel)
      deallocate(crm_input_phis)
      deallocate(crm_input_ul)
      deallocate(crm_input_vl)

      deallocate(crm_input_ocnfrac)
      deallocate(crm_input_tau00)
      deallocate(crm_input_wndls)
      deallocate(crm_input_bflxls)
      deallocate(crm_input_fluxu00)
      deallocate(crm_input_fluxv00)
      deallocate(crm_input_fluxt00)
      deallocate(crm_input_fluxq00)
#if defined( m2005 ) && defined( MODAL_AERO )
      deallocate(crm_input_naermod)
      deallocate(crm_input_vaerosol)
      deallocate(crm_input_hygro)
#endif
#if defined(SP_ESMT)
      deallocate(crm_input_ul_esmt)
      deallocate(crm_input_vl_esmt)
#endif
   end subroutine crm_input_finalize 
end module crm_input_module
