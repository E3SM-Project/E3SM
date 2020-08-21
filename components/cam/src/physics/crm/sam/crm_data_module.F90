module crm_data_mod
   use params, only: crm_rknd
#ifdef MODAL_AERO
   use modal_aero_data, only: ntot_amode
#endif
   use crmdims, only: crm_nx_rad, crm_ny_rad, &
                      crm_nx, crm_ny, crm_nz
#if defined(_OPENACC)
   use openacc_utils
#endif
   implicit none

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
   real(crm_rknd), allocatable :: crm_input_tau00  (:)          ! large-scale surface stress (N/m2)
   real(crm_rknd), allocatable :: crm_input_wndls  (:)          ! large-scale surface wind (m/s)
   real(crm_rknd), allocatable :: crm_input_bflxls (:)          ! large-scale surface buoyancy flux (K m/s)
   real(crm_rknd), allocatable :: crm_input_fluxu00(:)          ! surface momenent fluxes [N/m2]
   real(crm_rknd), allocatable :: crm_input_fluxv00(:)          ! surface momenent fluxes [N/m2]
   real(crm_rknd), allocatable :: crm_input_fluxt00(:)          ! surface sensible heat fluxes [K Kg/ (m2 s)]
   real(crm_rknd), allocatable :: crm_input_fluxq00(:)          ! surface latent heat fluxes [ kg/(m2 s)]

#if defined( m2005 ) && defined( MODAL_AERO )
   real(crm_rknd), allocatable :: crm_input_naermod (:,:,:)     ! Aerosol number concentration [/m3]
   real(crm_rknd), allocatable :: crm_input_vaerosol(:,:,:)     ! aerosol volume concentration [m3/m3]
   real(crm_rknd), allocatable :: crm_input_hygro   (:,:,:)     ! hygroscopicity of aerosol mode 
#endif

#if defined( MMF_ESMT )
   real(crm_rknd), allocatable :: crm_input_ul_esmt(:,:)        ! input u for ESMT
   real(crm_rknd), allocatable :: crm_input_vl_esmt(:,:)        ! input v for ESMT
#endif

   real(crm_rknd), allocatable :: crm_rad_temperature  (:,:,:,:)
   real(crm_rknd), allocatable :: crm_rad_qv           (:,:,:,:)
   real(crm_rknd), allocatable :: crm_rad_qc           (:,:,:,:)
   real(crm_rknd), allocatable :: crm_rad_qi           (:,:,:,:)
   real(crm_rknd), allocatable :: crm_rad_cld          (:,:,:,:)
   real(crm_rknd), allocatable :: crm_rad_qrad         (:,:,:,:)

   real(crm_rknd), allocatable :: crm_state_u_wind     (:,:,:,:)
   real(crm_rknd), allocatable :: crm_state_v_wind     (:,:,:,:)
   real(crm_rknd), allocatable :: crm_state_w_wind     (:,:,:,:)
   real(crm_rknd), allocatable :: crm_state_temperature(:,:,:,:)
   real(crm_rknd), allocatable :: crm_state_qt         (:,:,:,:)
   real(crm_rknd), allocatable :: crm_state_qp         (:,:,:,:)
   real(crm_rknd), allocatable :: crm_state_qn         (:,:,:,:)

   ! These are copies of the SAM cloud and precip liquid and ice, previously
   ! passed in and out of crm() via qc_crm, qi_crm, etc.
   real(crm_rknd), allocatable :: crm_output_qcl(:,:,:,:)
   real(crm_rknd), allocatable :: crm_output_qci(:,:,:,:)
   real(crm_rknd), allocatable :: crm_output_qpl(:,:,:,:)
   real(crm_rknd), allocatable :: crm_output_qpi(:,:,:,:)

   real(crm_rknd), allocatable :: crm_output_tk (:,:,:,:)
   real(crm_rknd), allocatable :: crm_output_tkh(:,:,:,:)
   real(crm_rknd), allocatable :: crm_output_prec_crm(:,:,:) ! CRM precipiation rate (surface)

   ! 2-moment process rates
   real(crm_rknd), allocatable :: crm_output_wvar(:,:,:,:) ! vertical velocity variance (m/s)
   real(crm_rknd), allocatable :: crm_output_aut (:,:,:,:) ! cloud water autoconversion (1/s)
   real(crm_rknd), allocatable :: crm_output_acc (:,:,:,:) ! cloud water accretion (1/s)
   real(crm_rknd), allocatable :: crm_output_evpc(:,:,:,:) ! cloud water evaporation (1/s)
   real(crm_rknd), allocatable :: crm_output_evpr(:,:,:,:) ! rain evaporation (1/s)
   real(crm_rknd), allocatable :: crm_output_mlt (:,:,:,:) ! ice, snow, graupel melting (1/s)
   real(crm_rknd), allocatable :: crm_output_sub (:,:,:,:) ! ice, snow, graupel sublimation (1/s)
   real(crm_rknd), allocatable :: crm_output_dep (:,:,:,:) ! ice, snow, graupel deposition (1/s)
   real(crm_rknd), allocatable :: crm_output_con (:,:,:,:) ! cloud water condensation(1/s)

   ! Cloud area fractions
   real(crm_rknd), allocatable :: crm_output_cltot(:)  ! shaded cloud fraction
   real(crm_rknd), allocatable :: crm_output_clhgh(:)  ! shaded cloud fraction
   real(crm_rknd), allocatable :: crm_output_clmed(:)  ! shaded cloud fraction
   real(crm_rknd), allocatable :: crm_output_cllow(:)  ! shaded cloud fraction

   real(crm_rknd), allocatable :: crm_output_cldtop(:,:)  ! cloud top ... pressure???
   real(crm_rknd), allocatable :: crm_output_precc(:)   ! convective precipitation rate
   real(crm_rknd), allocatable :: crm_output_precl(:)   ! stratiform precipitation rate
   real(crm_rknd), allocatable :: crm_output_precsc(:)   ! convective snow precipitation rate
   real(crm_rknd), allocatable :: crm_output_precsl(:)   ! stratiform snow precipitation rate

   ! TODO: These diagnostics are currently all on the GCM vertical grid. I
   ! think this is
   ! misleading though, and overly complicates crm_module. I think the better
   ! thing to do would
   ! be to define everything within crm_module on the CRM grid, and then do
   ! the
   ! mapping/interpolation at the GCM (crm_physics) level. For now though, I
   ! am just copying
   ! these over directly to minimize chances of me making a mistake. Also,
   ! some of these probably
   ! do not need to be calculated here, and might make more sense to
   ! calculate at the
   ! crm_physics_tend level. For example, I think tendencies should be
   ! calculated in
   ! crm_physics_tend, from, for example, something like crm_crm_output_uwind
   ! - crm_input%uwind.
   real(crm_rknd), allocatable :: crm_output_qc_mean(:,:)  ! mean cloud water
   real(crm_rknd), allocatable :: crm_output_qi_mean(:,:)  ! mean cloud ice
   real(crm_rknd), allocatable :: crm_output_qs_mean(:,:)  ! mean snow
   real(crm_rknd), allocatable :: crm_output_qg_mean(:,:)  ! mean graupel
   real(crm_rknd), allocatable :: crm_output_qr_mean(:,:)  ! mean rain
#ifdef m2005
   real(crm_rknd), allocatable :: crm_output_nc_mean(:,:)  ! mean cloud water (#/kg)
   real(crm_rknd), allocatable :: crm_output_ni_mean(:,:)  ! mean cloud ice (#/kg)
   real(crm_rknd), allocatable :: crm_output_ns_mean(:,:)  ! mean snow (#/kg)
   real(crm_rknd), allocatable :: crm_output_ng_mean(:,:)  ! mean graupel (#/kg)
   real(crm_rknd), allocatable :: crm_output_nr_mean(:,:)  ! mean rain (#/kg)

   ! Time and domain averaged process rates
   real(crm_rknd), allocatable :: crm_output_aut_a (:,:)  ! cloud water autoconversion (1/s)
   real(crm_rknd), allocatable :: crm_output_acc_a (:,:)  ! cloud water accretion (1/s)
   real(crm_rknd), allocatable :: crm_output_evpc_a(:,:)  ! cloud water evaporation (1/s)
   real(crm_rknd), allocatable :: crm_output_evpr_a(:,:)  ! rain evaporation (1/s)
   real(crm_rknd), allocatable :: crm_output_mlt_a (:,:)  ! ice, snow, graupel melting (1/s)
   real(crm_rknd), allocatable :: crm_output_sub_a (:,:)  ! ice, snow, graupel sublimation (1/s)
   real(crm_rknd), allocatable :: crm_output_dep_a (:,:)  ! ice, snow, graupel deposition (1/s)
   real(crm_rknd), allocatable :: crm_output_con_a (:,:)  ! cloud water condensation(1/s)
#endif /* m2005 */

#if defined( MMF_MOMENTUM_FEEDBACK )
   real(crm_rknd), allocatable :: crm_output_ultend(:,:)            ! tendency of ul
   real(crm_rknd), allocatable :: crm_output_vltend(:,:)            ! tendency of vl
#endif

#if defined( MMF_ESMT )
   real(crm_rknd), allocatable :: crm_output_u_tend_esmt(:,:)       ! CRM scalar u-momentum tendency
   real(crm_rknd), allocatable :: crm_output_v_tend_esmt(:,:)       ! CRM scalar v-momentum tendency
#endif

   real(crm_rknd), allocatable :: crm_output_sltend  (:,:)          ! CRM output tendency of static energy
   real(crm_rknd), allocatable :: crm_output_qltend  (:,:)          ! CRM output tendency of water vapor
   real(crm_rknd), allocatable :: crm_output_qcltend (:,:)          ! CRM output tendency of cloud liquid water
   real(crm_rknd), allocatable :: crm_output_qiltend (:,:)          ! CRM output tendency of cloud ice

   ! These are all time and spatial averages, on the GCM grid
   real(crm_rknd), allocatable :: crm_output_cld   (:,:)      ! cloud fraction
   real(crm_rknd), allocatable :: crm_output_gicewp(:,:)      ! ice water path
   real(crm_rknd), allocatable :: crm_output_gliqwp(:,:)      ! ice water path
   real(crm_rknd), allocatable :: crm_output_mctot (:,:)      ! cloud mass flux
   real(crm_rknd), allocatable :: crm_output_mcup  (:,:)      ! updraft cloud mass flux
   real(crm_rknd), allocatable :: crm_output_mcdn  (:,:)      ! downdraft cloud mass flux
   real(crm_rknd), allocatable :: crm_output_mcuup (:,:)      ! unsat updraft cloud mass flux
   real(crm_rknd), allocatable :: crm_output_mcudn (:,:)      ! unsat downdraft cloud mass flux

   ! For convective transport
   real(crm_rknd), allocatable :: crm_output_mu_crm(:,:)      ! mass flux up
   real(crm_rknd), allocatable :: crm_output_md_crm(:,:)      ! mass flux down
   real(crm_rknd), allocatable :: crm_output_du_crm(:,:)      ! mass detrainment from updraft
   real(crm_rknd), allocatable :: crm_output_eu_crm(:,:)      ! mass entrainment from updraft
   real(crm_rknd), allocatable :: crm_output_ed_crm(:,:)      ! mass detrainment from downdraft
   real(crm_rknd), allocatable :: crm_output_jt_crm(:)        ! index of cloud (convection) top
   real(crm_rknd), allocatable :: crm_output_mx_crm(:)        ! index of cloud (convection) bottom

   ! Other stuff...
   real(crm_rknd), allocatable :: crm_output_flux_qt      (:,:)  ! nonprecip water flux        [kg/m2/s]
   real(crm_rknd), allocatable :: crm_output_fluxsgs_qt   (:,:)  ! sgs non-precip water flux   [kg/m2/s]
   real(crm_rknd), allocatable :: crm_output_tkez         (:,:)  ! tke profile                 [kg/m/s2]
   real(crm_rknd), allocatable :: crm_output_tkesgsz      (:,:)  ! sgs tke profile             [kg/m/s2]
   real(crm_rknd), allocatable :: crm_output_tkz          (:,:)  ! tk profile [m2/s]
   real(crm_rknd), allocatable :: crm_output_flux_u       (:,:)  ! x-momentum flux             [m2/s2]
   real(crm_rknd), allocatable :: crm_output_flux_v       (:,:)  ! y-momentum flux             [m2/s2]
   real(crm_rknd), allocatable :: crm_output_flux_qp      (:,:)  ! precipitating water flux    [kg/m2/s or mm/s]
   real(crm_rknd), allocatable :: crm_output_precflux     (:,:)  ! precipitation flux          [m/s]
   real(crm_rknd), allocatable :: crm_output_qt_ls        (:,:)  ! tend of nonprec water due to large-scale   [kg/kg/s]
   real(crm_rknd), allocatable :: crm_output_qt_trans     (:,:)  ! tend of nonprec water due to transport     [kg/kg/s]
   real(crm_rknd), allocatable :: crm_output_qp_trans     (:,:)  ! tend of prec water due to transport     [kg/kg/s]
   real(crm_rknd), allocatable :: crm_output_qp_fall      (:,:)  ! tend of prec water due to fall-out      [kg/kg/s]
   real(crm_rknd), allocatable :: crm_output_qp_src       (:,:)  ! tend of prec water due to conversion    [kg/kg/s]
   real(crm_rknd), allocatable :: crm_output_qp_evp       (:,:)  ! tend of prec water due to evp           [kg/kg/s]
   real(crm_rknd), allocatable :: crm_output_t_ls         (:,:)  ! tend of lwse  due to large-scale           [kg/kg/s] ???
   real(crm_rknd), allocatable :: crm_output_prectend     (:)    ! column integrated tend in precip water+ice [kg/m2/s]
   real(crm_rknd), allocatable :: crm_output_precstend    (:)    ! column integrated tend in precip ice       [kg/m2/s]
   real(crm_rknd), allocatable :: crm_output_taux     (:)    ! zonal CRM surface stress perturbation      [N/m2]
   real(crm_rknd), allocatable :: crm_output_tauy     (:)    ! merid CRM surface stress perturbation      [N/m2]
   real(crm_rknd), allocatable :: crm_output_z0m          (:)    ! surface stress                             [N/m2]
   real(crm_rknd), allocatable :: crm_output_timing_factor(:)    ! crm cpu efficiency
#ifdef MAML
   ! MAML variables
   real(crm_rknd), allocatable :: crm_output_crm_pcp(ncrms,crm_nx,crm_ny) ! CRM precip rate for MAML (m/s)
   real(crm_rknd), allocatable :: crm_output_crm_snw(ncrms,crm_nx,crm_ny) ! CRM snow rate for MAML (m/s)
#endif
   !------------------------------------------------------------------------------------------------

contains
   !------------------------------------------------------------------------------------------------
   ! Type-bound procedures for crm_input_type
   subroutine allocate_input(ncrms, nlev)
      integer, intent(in) :: ncrms, nlev
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

#if defined(MMF_ESMT)
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
#if defined(MMF_ESMT)
      call prefetch(crm_input_ul_esmt)
      call prefetch(crm_input_vl_esmt)
#endif
#elif defined(_OPENMP)


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
#if defined( MMF_ESMT )
      crm_input_ul_esmt = 0
      crm_input_vl_esmt = 0
#endif
   end subroutine allocate_input

   subroutine allocate_states(ncrms)
     integer, intent(in) :: ncrms

     if (.not. allocated(crm_state_u_wind))       allocate( crm_state_u_wind(ncrms,crm_nx,crm_ny,crm_nz) )
     if (.not. allocated(crm_state_v_wind))       allocate( crm_state_v_wind(ncrms,crm_nx,crm_ny,crm_nz) )
     if (.not. allocated(crm_state_w_wind))       allocate( crm_state_w_wind(ncrms,crm_nx,crm_ny,crm_nz) )
     if (.not. allocated(crm_state_temperature))  allocate( crm_state_temperature(ncrms,crm_nx,crm_ny,crm_nz) )
     if (.not. allocated(crm_state_qt))           allocate( crm_state_qt(ncrms,crm_nx,crm_ny,crm_nz) )
     if (.not. allocated(crm_state_qp))           allocate( crm_state_qp(ncrms,crm_nx,crm_ny,crm_nz) )
     if (.not. allocated(crm_state_qn))           allocate( crm_state_qn(ncrms,crm_nx,crm_ny,crm_nz) )

#if defined(_OPENACC)
     call prefetch(crm_state_u_wind)
     call prefetch(crm_state_v_wind)
     call prefetch(crm_state_w_wind)
     call prefetch(crm_state_temperature)
     call prefetch(crm_state_qt)
     call prefetch(crm_state_qp)
     call prefetch(crm_state_qn)
#elif defined(_OPENMP)

#endif
   end subroutine allocate_states

   subroutine allocate_output(ncrms, nlev)
      integer, intent(in) :: ncrms, nlev

      ! Allocate instantaneous outputs
      if (.not. allocated(crm_output_qcl)) allocate(crm_output_qcl(ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(crm_output_qci)) allocate(crm_output_qci(ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(crm_output_qpl)) allocate(crm_output_qpl(ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(crm_output_qpi)) allocate(crm_output_qpi(ncrms,crm_nx,crm_ny,crm_nz))

      if (.not. allocated(crm_output_tk )) allocate(crm_output_tk (ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(crm_output_tkh)) allocate(crm_output_tkh(ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(crm_output_prec_crm)) allocate(crm_output_prec_crm(ncrms,crm_nx,crm_ny))

      if (.not. allocated(crm_output_wvar)) allocate(crm_output_wvar(ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(crm_output_aut))  allocate(crm_output_aut (ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(crm_output_acc))  allocate(crm_output_acc (ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(crm_output_evpc)) allocate(crm_output_evpc(ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(crm_output_evpr)) allocate(crm_output_evpr(ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(crm_output_mlt))  allocate(crm_output_mlt (ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(crm_output_sub))  allocate(crm_output_sub (ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(crm_output_dep))  allocate(crm_output_dep (ncrms,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(crm_output_con))  allocate(crm_output_con (ncrms,crm_nx,crm_ny,crm_nz))

      ! Allocate domain and time-averaged fields
      if (.not. allocated(crm_output_cltot)) allocate(crm_output_cltot(ncrms))
      if (.not. allocated(crm_output_cllow)) allocate(crm_output_cllow(ncrms))
      if (.not. allocated(crm_output_clmed)) allocate(crm_output_clmed(ncrms))
      if (.not. allocated(crm_output_clhgh)) allocate(crm_output_clhgh(ncrms))

      if (.not. allocated(crm_output_precc)) allocate(crm_output_precc(ncrms))
      if (.not. allocated(crm_output_precl)) allocate(crm_output_precl(ncrms))
      if (.not. allocated(crm_output_precsc)) allocate(crm_output_precsc(ncrms))
      if (.not. allocated(crm_output_precsl)) allocate(crm_output_precsl(ncrms))

      ! NOTE: this output had a bug in the previous implementation
      if (.not. allocated(crm_output_cldtop)) allocate(crm_output_cldtop(ncrms,nlev))

      if (.not. allocated(crm_output_qc_mean)) allocate(crm_output_qc_mean(ncrms,nlev))
      if (.not. allocated(crm_output_qi_mean)) allocate(crm_output_qi_mean(ncrms,nlev))
      if (.not. allocated(crm_output_qs_mean)) allocate(crm_output_qs_mean(ncrms,nlev))
      if (.not. allocated(crm_output_qg_mean)) allocate(crm_output_qg_mean(ncrms,nlev))
      if (.not. allocated(crm_output_qr_mean)) allocate(crm_output_qr_mean(ncrms,nlev))

#ifdef m2005
      if (.not. allocated(crm_output_nc_mean)) allocate(crm_output_nc_mean(ncrms,nlev))
      if (.not. allocated(crm_output_ni_mean)) allocate(crm_output_ni_mean(ncrms,nlev))
      if (.not. allocated(crm_output_ns_mean)) allocate(crm_output_ns_mean(ncrms,nlev))
      if (.not. allocated(crm_output_ng_mean)) allocate(crm_output_ng_mean(ncrms,nlev))
      if (.not. allocated(crm_output_nr_mean)) allocate(crm_output_nr_mean(ncrms,nlev))

      if (.not. allocated(crm_output_aut_a )) allocate(crm_output_aut_a (ncrms,nlev))
      if (.not. allocated(crm_output_acc_a )) allocate(crm_output_acc_a (ncrms,nlev))
      if (.not. allocated(crm_output_evpc_a)) allocate(crm_output_evpc_a(ncrms,nlev))
      if (.not. allocated(crm_output_evpr_a)) allocate(crm_output_evpr_a(ncrms,nlev))
      if (.not. allocated(crm_output_mlt_a )) allocate(crm_output_mlt_a (ncrms,nlev))
      if (.not. allocated(crm_output_sub_a )) allocate(crm_output_sub_a (ncrms,nlev))
      if (.not. allocated(crm_output_dep_a )) allocate(crm_output_dep_a (ncrms,nlev))
      if (.not. allocated(crm_output_con_a )) allocate(crm_output_con_a (ncrms,nlev))
#endif /* m2005 */

#if defined( MMF_MOMENTUM_FEEDBACK )
      if (.not. allocated(crm_output_ultend )) allocate(crm_output_ultend (ncrms,nlev))
      if (.not. allocated(crm_output_vltend )) allocate(crm_output_vltend (ncrms,nlev))
#endif

#if defined( MMF_ESMT )
      if (.not. allocated(crm_output_u_tend_esmt )) allocate(crm_output_u_tend_esmt (ncrms,nlev))
      if (.not. allocated(crm_output_v_tend_esmt )) allocate(crm_output_v_tend_esmt (ncrms,nlev))
#endif

#ifdef MAML
      if (.not. allocated(crm_output_crm_pcp)) allocate(crm_output_crm_pcp(ncrms,crm_nx,crm_ny))
      if (.not. allocated(crm_output_crm_snw)) allocate(crm_output_crm_snw(ncrms,crm_nx,crm_ny))
#endif
      if (.not. allocated(crm_output_sltend )) allocate(crm_output_sltend (ncrms,nlev))
      if (.not. allocated(crm_output_qltend )) allocate(crm_output_qltend (ncrms,nlev))
      if (.not. allocated(crm_output_qcltend)) allocate(crm_output_qcltend(ncrms,nlev))
      if (.not. allocated(crm_output_qiltend)) allocate(crm_output_qiltend(ncrms,nlev))
      if (.not. allocated(crm_output_cld   )) allocate(crm_output_cld (ncrms,nlev))  ! cloud fraction
      if (.not. allocated(crm_output_gicewp)) allocate(crm_output_gicewp(ncrms,nlev))  ! ice water path
      if (.not. allocated(crm_output_gliqwp)) allocate(crm_output_gliqwp(ncrms,nlev))  ! ice water path
      if (.not. allocated(crm_output_mctot )) allocate(crm_output_mctot (ncrms,nlev))  ! cloud mass flux
      if (.not. allocated(crm_output_mcup  )) allocate(crm_output_mcup (ncrms,nlev))  ! updraft cloud mass flux
      if (.not. allocated(crm_output_mcdn  )) allocate(crm_output_mcdn (ncrms,nlev))  ! downdraft cloud mass flux
      if (.not. allocated(crm_output_mcuup )) allocate(crm_output_mcuup (ncrms,nlev))  ! unsat updraft cloud mass flux
      if (.not. allocated(crm_output_mcudn )) allocate(crm_output_mcudn (ncrms,nlev))  ! unsat downdraft cloud mass flux

      if (.not. allocated(crm_output_mu_crm)) allocate(crm_output_mu_crm(ncrms,nlev))  ! mass flux up
      if (.not. allocated(crm_output_md_crm)) allocate(crm_output_md_crm(ncrms,nlev))  ! mass flux down
      if (.not. allocated(crm_output_du_crm)) allocate(crm_output_du_crm(ncrms,nlev))  ! mass detrainment from updraft
      if (.not. allocated(crm_output_eu_crm)) allocate(crm_output_eu_crm(ncrms,nlev))  ! mass entrainment from updraft
      if (.not. allocated(crm_output_ed_crm)) allocate(crm_output_ed_crm(ncrms,nlev))  ! mass detrainment from downdraft
      if (.not. allocated(crm_output_jt_crm)) allocate(crm_output_jt_crm(ncrms))       ! index of cloud (convection) top
      if (.not. allocated(crm_output_mx_crm)) allocate(crm_output_mx_crm(ncrms))       ! index of cloud (convection) bottom

      if (.not. allocated(crm_output_flux_qt      )) allocate(crm_output_flux_qt      (ncrms,nlev))
      if (.not. allocated(crm_output_fluxsgs_qt   )) allocate(crm_output_fluxsgs_qt   (ncrms,nlev))
      if (.not. allocated(crm_output_tkez         )) allocate(crm_output_tkez (ncrms,nlev))
      if (.not. allocated(crm_output_tkesgsz      )) allocate(crm_output_tkesgsz      (ncrms,nlev))
      if (.not. allocated(crm_output_tkz          )) allocate(crm_output_tkz (ncrms,nlev))
      if (.not. allocated(crm_output_flux_u       )) allocate(crm_output_flux_u       (ncrms,nlev))
      if (.not. allocated(crm_output_flux_v       )) allocate(crm_output_flux_v       (ncrms,nlev))
      if (.not. allocated(crm_output_flux_qp      )) allocate(crm_output_flux_qp      (ncrms,nlev))
      if (.not. allocated(crm_output_precflux     )) allocate(crm_output_precflux     (ncrms,nlev))
      if (.not. allocated(crm_output_qt_ls        )) allocate(crm_output_qt_ls        (ncrms,nlev))
      if (.not. allocated(crm_output_qt_trans     )) allocate(crm_output_qt_trans     (ncrms,nlev))
      if (.not. allocated(crm_output_qp_trans     )) allocate(crm_output_qp_trans     (ncrms,nlev))
      if (.not. allocated(crm_output_qp_fall      )) allocate(crm_output_qp_fall      (ncrms,nlev))
      if (.not. allocated(crm_output_qp_src       )) allocate(crm_output_qp_src       (ncrms,nlev))
      if (.not. allocated(crm_output_qp_evp       )) allocate(crm_output_qp_evp       (ncrms,nlev))
      if (.not. allocated(crm_output_t_ls         )) allocate(crm_output_t_ls (ncrms,nlev))
      if (.not. allocated(crm_output_prectend     )) allocate(crm_output_prectend     (ncrms))
      if (.not. allocated(crm_output_precstend    )) allocate(crm_output_precstend    (ncrms))
      if (.not. allocated(crm_output_taux         )) allocate(crm_output_taux (ncrms))
      if (.not. allocated(crm_output_tauy         )) allocate(crm_output_tauy (ncrms))
      if (.not. allocated(crm_output_z0m          )) allocate(crm_output_z0m (ncrms))
      if (.not. allocated(crm_output_timing_factor)) allocate(crm_output_timing_factor(ncrms))

#if defined(_OPENACC)
      call prefetch(crm_output_qcl)
      call prefetch(crm_output_qci)
      call prefetch(crm_output_qpl)
      call prefetch(crm_output_qpi)
      call prefetch(crm_output_tk )
      call prefetch(crm_output_tkh)
      call prefetch(crm_output_prec_crm)
      call prefetch(crm_output_wvar)
      call prefetch(crm_output_aut )
      call prefetch(crm_output_acc )
      call prefetch(crm_output_evpc)
      call prefetch(crm_output_evpr)
      call prefetch(crm_output_mlt )
      call prefetch(crm_output_sub )
      call prefetch(crm_output_dep )
      call prefetch(crm_output_con )
      call prefetch(crm_output_cltot)
      call prefetch(crm_output_cllow)
      call prefetch(crm_output_clmed)
      call prefetch(crm_output_clhgh)
      call prefetch(crm_output_precc)
      call prefetch(crm_output_precl)
      call prefetch(crm_output_precsc)
      call prefetch(crm_output_precsl)
      call prefetch(crm_output_cldtop)
      call prefetch(crm_output_qc_mean)
      call prefetch(crm_output_qi_mean)
      call prefetch(crm_output_qs_mean)
      call prefetch(crm_output_qg_mean)
      call prefetch(crm_output_qr_mean)
#ifdef MAML
      call prefetch(crm_output_crm_pcp)
      call prefetch(crm_output_crm_snw)
#endif
      call prefetch(crm_output_sltend  )
      call prefetch(crm_output_qltend  )
      call prefetch(crm_output_qcltend )
      call prefetch(crm_output_qiltend )
      call prefetch(crm_output_cld    )
      call prefetch(crm_output_gicewp )
      call prefetch(crm_output_gliqwp )
      call prefetch(crm_output_mctot  )
      call prefetch(crm_output_mcup   )
      call prefetch(crm_output_mcdn   )
      call prefetch(crm_output_mcuup  )
      call prefetch(crm_output_mcudn  )
      call prefetch(crm_output_mu_crm )
      call prefetch(crm_output_md_crm )
      call prefetch(crm_output_du_crm )
      call prefetch(crm_output_eu_crm )
      call prefetch(crm_output_ed_crm )
      call prefetch(crm_output_jt_crm )
      call prefetch(crm_output_mx_crm )
      call prefetch(crm_output_flux_qt       )
      call prefetch(crm_output_fluxsgs_qt    )
      call prefetch(crm_output_tkez          )
      call prefetch(crm_output_tkesgsz       )
      call prefetch(crm_output_tkz           )
      call prefetch(crm_output_flux_u        )
      call prefetch(crm_output_flux_v        )
      call prefetch(crm_output_flux_qp       )
      call prefetch(crm_output_precflux      )
      call prefetch(crm_output_qt_ls         )
      call prefetch(crm_output_qt_trans      )
      call prefetch(crm_output_qp_trans      )
      call prefetch(crm_output_qp_fall       )
      call prefetch(crm_output_qp_src        )
      call prefetch(crm_output_qp_evp        )
      call prefetch(crm_output_t_ls          )
      call prefetch(crm_output_prectend      )
      call prefetch(crm_output_precstend     )
      call prefetch(crm_output_taux          )
      call prefetch(crm_output_tauy          )
      call prefetch(crm_output_z0m           )
#endif
      ! Initialize 
      crm_output_qcl = 0
      crm_output_qci = 0
      crm_output_qpl = 0
      crm_output_qpi = 0

      crm_output_tk = 0
      crm_output_tkh = 0
      crm_output_prec_crm = 0

      ! 2-moment process rates
      crm_output_wvar = 0
      crm_output_aut  = 0
      crm_output_acc  = 0
      crm_output_evpc = 0
      crm_output_evpr = 0
      crm_output_mlt  = 0
      crm_output_sub  = 0
      crm_output_dep  = 0
      crm_output_con  = 0

      crm_output_cltot = 0
      crm_output_cllow = 0
      crm_output_clmed = 0
      crm_output_clhgh = 0

      crm_output_cldtop = 0
      crm_output_precc = 0
      crm_output_precl = 0
      crm_output_precsc = 0
      crm_output_precsl = 0

      crm_output_qc_mean = 0
      crm_output_qi_mean = 0
      crm_output_qs_mean = 0
      crm_output_qg_mean = 0
      crm_output_qr_mean = 0
#ifdef m2005
      crm_output_nc_mean = 0
      crm_output_ni_mean = 0
      crm_output_ns_mean = 0
      crm_output_ng_mean = 0
      crm_output_nr_mean = 0

      crm_output_aut_a = 0
      crm_output_acc_a = 0
      crm_output_evpc_a = 0
      crm_output_evpr_a = 0
      crm_output_mlt_a = 0
      crm_output_sub_a = 0
      crm_output_dep_a = 0
      crm_output_con_a = 0
#endif

#if defined( MMF_MOMENTUM_FEEDBACK )
      crm_output_ultend = 0
      crm_output_vltend = 0
#endif

#if defined( MMF_ESMT )
      crm_output_u_tend_esmt = 0
      crm_output_v_tend_esmt = 0
#endif

#ifdef MAML
      crm_output_crm_pcp = 0
      crm_output_crm_snw = 0
#endif

      crm_output_sltend  = 0
      crm_output_qltend  = 0
      crm_output_qcltend = 0
      crm_output_qiltend = 0
      crm_output_cld    = 0
      crm_output_gicewp = 0
      crm_output_gliqwp = 0
      crm_output_mctot  = 0
      crm_output_mcup   = 0
      crm_output_mcdn   = 0
      crm_output_mcuup  = 0
      crm_output_mcudn  = 0

      ! Convective transport
      crm_output_mu_crm = 0
      crm_output_md_crm = 0
      crm_output_eu_crm = 0
      crm_output_du_crm = 0
      crm_output_ed_crm = 0
      crm_output_jt_crm = 0
      crm_output_mx_crm = 0

      ! Other stuff...
      crm_output_flux_qt       = 0
      crm_output_fluxsgs_qt    = 0
      crm_output_tkez          = 0
      crm_output_tkesgsz       = 0
      crm_output_tkz           = 0
      crm_output_flux_u        = 0
      crm_output_flux_v        = 0
      crm_output_flux_qp       = 0
      crm_output_precflux      = 0
      crm_output_qt_ls         = 0
      crm_output_qt_trans      = 0
      crm_output_qp_trans      = 0
      crm_output_qp_fall       = 0
      crm_output_qp_src        = 0
      crm_output_qp_evp        = 0
      crm_output_t_ls          = 0
      crm_output_prectend      = 0
      crm_output_precstend     = 0
      crm_output_taux      = 0
      crm_output_tauy      = 0
      crm_output_z0m           = 0
      crm_output_timing_factor = 0
   end subroutine allocate_output

   subroutine allocate_rad(ncrms)
     integer, intent(in) :: ncrms

     if (.not. allocated(crm_rad_temperature)) allocate( crm_rad_temperature(ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
     if (.not. allocated(crm_rad_qv))          allocate( crm_rad_qv(ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
     if (.not. allocated(crm_rad_qc))          allocate( crm_rad_qc(ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
     if (.not. allocated(crm_rad_qi))          allocate( crm_rad_qi(ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
     if (.not. allocated(crm_rad_cld))         allocate( crm_rad_cld(ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
     if (.not. allocated(crm_rad_qrad))        allocate( crm_rad_qrad(ncrms,crm_nx_rad,crm_ny_rad,crm_nz) )
#if defined(_OPENACC)
     call prefetch(crm_rad_temperature)
     call prefetch(crm_rad_qv)
     call prefetch(crm_rad_qc)
     call prefetch(crm_rad_qi)
     call prefetch(crm_rad_cld)
     call prefetch(crm_rad_qrad)
#elif defined(_OPENMP)

#endif   
   end subroutine allocate_rad

   subroutine get_input_data(input, ncrms, nlev)
      use crm_input_module
      integer, intent(in) :: ncrms, nlev
      type(crm_input_type), intent(in) :: input

      crm_input_zmid     = input%zmid(1:ncrms,1:nlev)
      crm_input_zint     = input%zint(1:ncrms,1:nlev+1)
      crm_input_tl       = input%tl(1:ncrms,1:nlev)
      crm_input_ql       = input%ql(1:ncrms,1:nlev)
      crm_input_qccl     = input%qccl(1:ncrms,1:nlev)
      crm_input_qiil     = input%qiil(1:ncrms,1:nlev)
      crm_input_ps       = input%ps(1:ncrms)
      crm_input_pmid     = input%pmid(1:ncrms,1:nlev)
      crm_input_pint     = input%pint(1:ncrms,1:nlev+1)
      crm_input_pdel     = input%pdel(1:ncrms,1:nlev)
      crm_input_phis     = input%phis(1:ncrms)
      crm_input_ul       = input%ul(1:ncrms,1:nlev)
      crm_input_vl       = input%vl(1:ncrms,1:nlev)
      crm_input_ocnfrac  = input%ocnfrac(1:ncrms)
      crm_input_tau00    = input%tau00(1:ncrms)
      crm_input_wndls    = input%wndls(1:ncrms)
      crm_input_bflxls   = input%bflxls(1:ncrms)
      crm_input_fluxu00  = input%fluxu00(1:ncrms)
      crm_input_fluxv00  = input%fluxv00(1:ncrms)
      crm_input_fluxt00  = input%fluxt00(1:ncrms)
      crm_input_fluxq00  = input%fluxq00(1:ncrms)

#if defined( m2005 ) && defined( MODAL_AERO )
      crm_input_naermod  = input%naermod(1:ncrms,1:nlev,1:ntot_amode)
      crm_input_vaerosol = input%vaerosol(1:ncrms,1:nlev,1:ntot_amode)
      crm_input_hygro    = input%hygro(1:ncrms,1:nlev,1:ntot_amode)
#endif

#if defined(MMF_ESMT)
      crm_input_ul_esmt = input%ul_esmt(1:ncrms,1:nlev)
      crm_input_vl_esmt = input%vl_esmt(1:ncrms,1:nlev)
#endif
   end subroutine get_input_data

   subroutine get_rad_data(rad, ncrms)
    use crm_rad_module,only: crm_rad_type
     implicit none
     integer, intent(in) :: ncrms
     type(crm_rad_type), intent(out) :: rad
     crm_rad_temperature = rad%temperature(1:ncrms,1:crm_nx_rad,1:crm_ny_rad,1:crm_nz)
     crm_rad_qv          = rad%qv(1:ncrms,1:crm_nx_rad,1:crm_ny_rad,1:crm_nz)
     crm_rad_qc          = rad%qc(1:ncrms,1:crm_nx_rad,1:crm_ny_rad,1:crm_nz)
     crm_rad_qi          = rad%qi(1:ncrms,1:crm_nx_rad,1:crm_ny_rad,1:crm_nz)
     crm_rad_cld         = rad%cld(1:ncrms,1:crm_nx_rad,1:crm_ny_rad,1:crm_nz)
     crm_rad_qrad        = rad%qrad(1:ncrms,1:crm_nx_rad,1:crm_ny_rad,1:crm_nz)
   end subroutine get_rad_data

   subroutine get_state_data(state, ncrms)
     use crm_state_module,only: crm_state_type
     implicit none
     integer, intent(in) :: ncrms
     type(crm_state_type), intent(in) :: state
     crm_state_u_wind      = state%u_wind(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz)
     crm_state_v_wind      = state%v_wind(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz)
     crm_state_w_wind      = state%w_wind(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz)
     crm_state_temperature = state%temperature(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz)
     crm_state_qt          = state%qt(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz)
     crm_state_qp          = state%qp(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz)
     crm_state_qn          = state%qn(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz)
   end subroutine get_state_data

   subroutine load_rad_data(rad, ncrms)
     use crm_rad_module,only: crm_rad_type
     implicit none
     integer, intent(in) :: ncrms
     type(crm_rad_type) :: rad
     rad%temperature(1:ncrms,1:crm_nx_rad,1:crm_ny_rad,1:crm_nz) = crm_rad_temperature
     rad%qv(1:ncrms,1:crm_nx_rad,1:crm_ny_rad,1:crm_nz) = crm_rad_qv
     rad%qc(1:ncrms,1:crm_nx_rad,1:crm_ny_rad,1:crm_nz) = crm_rad_qc
     rad%qi(1:ncrms,1:crm_nx_rad,1:crm_ny_rad,1:crm_nz) = crm_rad_qi
     rad%cld(1:ncrms,1:crm_nx_rad,1:crm_ny_rad,1:crm_nz) = crm_rad_cld
     rad%qrad(1:ncrms,1:crm_nx_rad,1:crm_ny_rad,1:crm_nz) = crm_rad_qrad
   end subroutine load_rad_data

   subroutine load_state_data(state, ncrms)
     use crm_state_module,only: crm_state_type
     implicit none
     integer, intent(in) :: ncrms
     type(crm_state_type) :: state
     state%u_wind(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_state_u_wind
     state%v_wind(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_state_v_wind
     state%w_wind(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_state_w_wind
     state%temperature(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_state_temperature
     state%qt(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_state_qt
     state%qp(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_state_qp
     state%qn(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_state_qn
   end subroutine load_state_data
  
   subroutine load_output_data(output, ncrms, nlev)
    use crm_output_module, only: crm_output_type
    implicit none
    integer, intent(in) :: ncrms, nlev
    type(crm_output_type) :: output
    integer i, j, k, icrm

    output%qcl(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_output_qcl
    output%qci(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_output_qci
    output%qpl(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_output_qpl
    output%qpi(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_output_qpi

    output%tk(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz)  = crm_output_tk
    output%tkh(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_output_tkh
    output%prec_crm(1:ncrms,1:crm_nx,1:crm_ny)   = crm_output_prec_crm

    output%wvar(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_output_wvar
    output%aut(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_output_aut
    output%acc(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_output_acc
    output%evpc(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_output_evpc
    output%evpr(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_output_evpr
    output%mlt(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_output_mlt
    output%sub(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_output_sub
    output%dep(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_output_dep
    output%con(1:ncrms,1:crm_nx,1:crm_ny,1:crm_nz) = crm_output_con

    ! Allocate domain and time-averaged fields
    output%cltot(1:ncrms) = crm_output_cltot
    output%cllow(1:ncrms) = crm_output_cllow
    output%clmed(1:ncrms) = crm_output_clmed 
    output%clhgh(1:ncrms) = crm_output_clhgh

    output%precc(1:ncrms) = crm_output_precc
    output%precl(1:ncrms) = crm_output_precl
    output%precsc(1:ncrms) = crm_output_precsc
    output%precsl(1:ncrms) = crm_output_precsl 

    ! NOTE: this output had a bug in the previous implementation
    output%cldtop(1:ncrms,1:nlev) = crm_output_cldtop

    output%qc_mean(1:ncrms,1:nlev) = crm_output_qc_mean
    output%qi_mean(1:ncrms,1:nlev) = crm_output_qi_mean
    output%qs_mean(1:ncrms,1:nlev) = crm_output_qs_mean
    output%qg_mean(1:ncrms,1:nlev) = crm_output_qg_mean
    output%qr_mean(1:ncrms,1:nlev) = crm_output_qr_mean

#ifdef m2005
    output%nc_mean(1:ncrms,1:nlev) = crm_output_nc_mean
    output%ni_mean(1:ncrms,1:nlev) = crm_output_ni_mean
    output%ns_mean(1:ncrms,1:nlev) = crm_output_ns_mean
    output%ng_mean(1:ncrms,1:nlev) = crm_output_ng_mean
    output%nr_mean(1:ncrms,1:nlev) = crm_output_nr_mean

    output%aut_a(1:ncrms,1:nlev) = crm_output_aut_a
    output%acc_a(1:ncrms,1:nlev) = crm_output_acc_a
    output%evpc_a(1:ncrms,1:nlev) = crm_output_evpc_a
    output%evpr_a(1:ncrms,1:nlev) = crm_output_evpr_a
    output%mlt_a(1:ncrms,1:nlev) = crm_output_mlt_a
    output%sub_a(1:ncrms,1:nlev) = crm_output_sub_a
    output%dep_a(1:ncrms,1:nlev) = crm_output_dep_a
    output%con_a(1:ncrms,1:nlev) = crm_output_con_a
#endif /* m2005 */

#if defined( MMF_MOMENTUM_FEEDBACK )
    output%ultend(1:ncrms,1:nlev) = crm_output_ultend
    output%vltend(1:ncrms,1:nlev) = crm_output_vltend
#endif

#if defined( MMF_ESMT )
    output%u_tend_esmt(1:ncrms,1:nlev) = crm_output_u_tend_esmt
    output%v_tend_esmt(1:ncrms,1:nlev) = crm_output_v_tend_esmt
#endif

#ifdef MAML
    output%crm_pcp(1:ncrms,1:crm_nx,1:crm_ny) = crm_output_crm_pcp
    output%crm_snw(1:ncrms,1:crm_nx,1:crm_ny) = crm_output_crm_snw
#endif

    output%sltend(1:ncrms,1:nlev) = crm_output_sltend
    output%qltend(1:ncrms,1:nlev) = crm_output_qltend
    output%qcltend(1:ncrms,1:nlev) = crm_output_qcltend
    output%qiltend(1:ncrms,1:nlev) = crm_output_qiltend
    output%cld(1:ncrms,1:nlev) = crm_output_cld
    output%gicewp(1:ncrms,1:nlev) = crm_output_gicewp
    output%gliqwp(1:ncrms,1:nlev) = crm_output_gliqwp
    output%mctot(1:ncrms,1:nlev) = crm_output_mctot
    output%mcup(1:ncrms,1:nlev) = crm_output_mcup
    output%mcdn(1:ncrms,1:nlev) = crm_output_mcdn
    output%mcuup(1:ncrms,1:nlev) = crm_output_mcuup
    output%mcudn(1:ncrms,1:nlev) = crm_output_mcudn

    output%mu_crm(1:ncrms,1:nlev) = crm_output_mu_crm
    output%md_crm(1:ncrms,1:nlev) = crm_output_md_crm
    output%du_crm(1:ncrms,1:nlev) = crm_output_du_crm
    output%eu_crm(1:ncrms,1:nlev) = crm_output_eu_crm
    output%ed_crm(1:ncrms,1:nlev) = crm_output_ed_crm
    output%jt_crm(1:ncrms) = crm_output_jt_crm
    output%mx_crm(1:ncrms) = crm_output_mx_crm

    output%flux_qt(1:ncrms,1:nlev) = crm_output_flux_qt
    output%fluxsgs_qt(1:ncrms,1:nlev) = crm_output_fluxsgs_qt
    output%tkez(1:ncrms,1:nlev) = crm_output_tkez
    output%tkesgsz(1:ncrms,1:nlev) = crm_output_tkesgsz
    output%tkz(1:ncrms,1:nlev) = crm_output_tkz
    output%flux_u(1:ncrms,1:nlev) = crm_output_flux_u
    output%flux_v(1:ncrms,1:nlev) = crm_output_flux_v
    output%flux_qp(1:ncrms,1:nlev) = crm_output_flux_qp
    output%precflux(1:ncrms,1:nlev) = crm_output_precflux
    output%qt_ls(1:ncrms,1:nlev) = crm_output_qt_ls
    output%qt_trans(1:ncrms,1:nlev) = crm_output_qt_trans
    output%qp_trans(1:ncrms,1:nlev) = crm_output_qp_trans
    output%qp_fall(1:ncrms,1:nlev) = crm_output_qp_fall
    output%qp_src(1:ncrms,1:nlev) = crm_output_qp_src
    output%qp_evp(1:ncrms,1:nlev) = crm_output_qp_evp
    output%t_ls(1:ncrms,1:nlev) = crm_output_t_ls
    output%prectend(1:ncrms) = crm_output_prectend
    output%precstend(1:ncrms) = crm_output_precstend
    output%taux(1:ncrms) = crm_output_taux
    output%tauy(1:ncrms) = crm_output_tauy
    output%z0m(1:ncrms) = crm_output_z0m
    output%timing_factor(1:ncrms) = crm_output_timing_factor
   end subroutine load_output_data 

   !------------------------------------------------------------------------------
   subroutine deallocate_input()
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

#if defined(MMF_ESMT)
     deallocate(crm_input_ul_esmt)
     deallocate(crm_input_vl_esmt)
#endif
   end subroutine deallocate_input

   subroutine deallocate_states()
     deallocate(crm_state_u_wind)
     deallocate(crm_state_v_wind)
     deallocate(crm_state_w_wind)
     deallocate(crm_state_temperature)
     deallocate(crm_state_qt)
     deallocate(crm_state_qp)
     deallocate(crm_state_qn)
   end subroutine deallocate_states

   subroutine deallocate_rad()
     deallocate(crm_rad_temperature)
     deallocate(crm_rad_qv)
     deallocate(crm_rad_qc)
     deallocate(crm_rad_qi)
     deallocate(crm_rad_cld)
     deallocate(crm_rad_qrad)
   end subroutine deallocate_rad

   subroutine deallocate_output()
     if (allocated(crm_output_qcl)) deallocate(crm_output_qcl)
     if (allocated(crm_output_qci)) deallocate(crm_output_qci)
     if (allocated(crm_output_qpl)) deallocate(crm_output_qpl)
     if (allocated(crm_output_qpi)) deallocate(crm_output_qpi)
     if (allocated(crm_output_tk )) deallocate(crm_output_tk )
     if (allocated(crm_output_tkh)) deallocate(crm_output_tkh)
     if (allocated(crm_output_prec_crm)) deallocate(crm_output_prec_crm)

     if (allocated(crm_output_wvar)) deallocate(crm_output_wvar)
     if (allocated(crm_output_aut)) deallocate(crm_output_aut)
     if (allocated(crm_output_acc)) deallocate(crm_output_acc)
     if (allocated(crm_output_evpc)) deallocate(crm_output_evpc)
     if (allocated(crm_output_evpr)) deallocate(crm_output_evpr)
     if (allocated(crm_output_mlt)) deallocate(crm_output_mlt)
     if (allocated(crm_output_sub)) deallocate(crm_output_sub)
     if (allocated(crm_output_dep)) deallocate(crm_output_dep)
     if (allocated(crm_output_con)) deallocate(crm_output_con)

     if (allocated(crm_output_cltot)) deallocate(crm_output_cltot)
     if (allocated(crm_output_cllow)) deallocate(crm_output_cllow)
     if (allocated(crm_output_clmed)) deallocate(crm_output_clmed)
     if (allocated(crm_output_clhgh)) deallocate(crm_output_clhgh)
     if (allocated(crm_output_cldtop)) deallocate(crm_output_cldtop)
     if (allocated(crm_output_precc)) deallocate(crm_output_precc)
     if (allocated(crm_output_precl)) deallocate(crm_output_precl)
     if (allocated(crm_output_precsc)) deallocate(crm_output_precsc)
     if (allocated(crm_output_precsl)) deallocate(crm_output_precsl)

     if (allocated(crm_output_qc_mean)) deallocate(crm_output_qc_mean)
     if (allocated(crm_output_qi_mean)) deallocate(crm_output_qi_mean)
     if (allocated(crm_output_qs_mean)) deallocate(crm_output_qs_mean)
     if (allocated(crm_output_qg_mean)) deallocate(crm_output_qg_mean)
     if (allocated(crm_output_qr_mean)) deallocate(crm_output_qr_mean)
#ifdef m2005
     if (allocated(crm_output_nc_mean)) deallocate(crm_output_nc_mean)
     if (allocated(crm_output_ni_mean)) deallocate(crm_output_ni_mean)
     if (allocated(crm_output_ns_mean)) deallocate(crm_output_ns_mean)
     if (allocated(crm_output_ng_mean)) deallocate(crm_output_ng_mean)
     if (allocated(crm_output_nr_mean)) deallocate(crm_output_nr_mean)

     ! Time and domain-averaged process rates
     if (allocated(crm_output_aut_a)) deallocate(crm_output_aut_a)
     if (allocated(crm_output_acc_a)) deallocate(crm_output_acc_a)
     if (allocated(crm_output_evpc_a)) deallocate(crm_output_evpc_a)
     if (allocated(crm_output_evpr_a)) deallocate(crm_output_evpr_a)
     if (allocated(crm_output_mlt_a)) deallocate(crm_output_mlt_a)
     if (allocated(crm_output_sub_a)) deallocate(crm_output_sub_a)
     if (allocated(crm_output_dep_a)) deallocate(crm_output_dep_a)
     if (allocated(crm_output_con_a)) deallocate(crm_output_con_a)
#endif

#if defined( MMF_MOMENTUM_FEEDBACK )
     if (allocated(crm_output_ultend)) deallocate(crm_output_ultend)
     if (allocated(crm_output_vltend)) deallocate(crm_output_vltend)
#endif

#if defined( MMF_ESMT )
     if (allocated(crm_output_u_tend_esmt)) deallocate(crm_output_u_tend_esmt)
     if (allocated(crm_output_v_tend_esmt)) deallocate(crm_output_v_tend_esmt)
#endif

#ifdef MAML
     if (allocated(crm_output_crm_pcp)) deallocate(crm_output_crm_pcp)
     if (allocated(crm_output_crm_snw)) deallocate(crm_output_crm_snw)
#endif
     if (allocated(crm_output_sltend)) deallocate(crm_output_sltend)
     if (allocated(crm_output_qltend)) deallocate(crm_output_qltend)
     if (allocated(crm_output_qcltend)) deallocate(crm_output_qcltend)
     if (allocated(crm_output_qiltend)) deallocate(crm_output_qiltend)

     if (allocated(crm_output_cld)) deallocate(crm_output_cld)
     if (allocated(crm_output_gicewp)) deallocate(crm_output_gicewp)
     if (allocated(crm_output_gliqwp)) deallocate(crm_output_gliqwp)
     if (allocated(crm_output_mctot)) deallocate(crm_output_mctot)
     if (allocated(crm_output_mcup)) deallocate(crm_output_mcup)
     if (allocated(crm_output_mcdn)) deallocate(crm_output_mcdn)
     if (allocated(crm_output_mcuup)) deallocate(crm_output_mcuup)
     if (allocated(crm_output_mcudn)) deallocate(crm_output_mcudn)

     if (allocated(crm_output_mu_crm)) deallocate(crm_output_mu_crm)
     if (allocated(crm_output_md_crm)) deallocate(crm_output_md_crm)
     if (allocated(crm_output_du_crm)) deallocate(crm_output_du_crm)
     if (allocated(crm_output_eu_crm)) deallocate(crm_output_eu_crm)
     if (allocated(crm_output_ed_crm)) deallocate(crm_output_ed_crm)
     if (allocated(crm_output_jt_crm)) deallocate(crm_output_jt_crm)
     if (allocated(crm_output_mx_crm)) deallocate(crm_output_mx_crm)

     if (allocated(crm_output_flux_qt)) deallocate(crm_output_flux_qt)
     if (allocated(crm_output_fluxsgs_qt)) deallocate(crm_output_fluxsgs_qt)
     if (allocated(crm_output_tkez)) deallocate(crm_output_tkez)
     if (allocated(crm_output_tkesgsz)) deallocate(crm_output_tkesgsz)
     if (allocated(crm_output_tkz)) deallocate(crm_output_tkz)
     if (allocated(crm_output_flux_u)) deallocate(crm_output_flux_u)
     if (allocated(crm_output_flux_v)) deallocate(crm_output_flux_v)
     if (allocated(crm_output_flux_qp)) deallocate(crm_output_flux_qp)
     if (allocated(crm_output_precflux)) deallocate(crm_output_precflux)
     if (allocated(crm_output_qt_ls)) deallocate(crm_output_qt_ls)
     if (allocated(crm_output_qt_trans)) deallocate(crm_output_qt_trans)
     if (allocated(crm_output_qp_trans)) deallocate(crm_output_qp_trans)
     if (allocated(crm_output_qp_fall)) deallocate(crm_output_qp_fall)
     if (allocated(crm_output_qp_src)) deallocate(crm_output_qp_src)
     if (allocated(crm_output_qp_evp)) deallocate(crm_output_qp_evp)
     if (allocated(crm_output_t_ls)) deallocate(crm_output_t_ls)
     if (allocated(crm_output_prectend)) deallocate(crm_output_prectend)
     if (allocated(crm_output_precstend)) deallocate(crm_output_precstend)
     if (allocated(crm_output_taux)) deallocate(crm_output_taux)
     if (allocated(crm_output_tauy)) deallocate(crm_output_tauy)
     if (allocated(crm_output_z0m)) deallocate(crm_output_z0m)
     if (allocated(crm_output_timing_factor)) deallocate(crm_output_timing_factor)
   end subroutine deallocate_output

end module crm_data_mod
