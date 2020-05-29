
module crm_output_module
   use params,       only: crm_rknd
   use crmdims,      only: crm_nx, crm_ny, crm_nz
#if defined(_OPENACC)
   use openacc_utils
#endif
   implicit none
   ! Derived type to encapsulate CRM output fields (things that are output
   ! only, intent(out) in the original implementation)

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

   ! TODO: These diagnostics are currently all on the GCM vertical grid. I think this is
   ! misleading though, and overly complicates crm_module. I think the better thing to do would
   ! be to define everything within crm_module on the CRM grid, and then do the
   ! mapping/interpolation at the GCM (crm_physics) level. For now though, I am just copying
   ! these over directly to minimize chances of me making a mistake. Also, some of these probably
   ! do not need to be calculated here, and might make more sense to calculate at the
   ! crm_physics_tend level. For example, I think tendencies should be calculated in
   ! crm_physics_tend, from, for example, something like crm_output%uwind - crm_input%uwind.
   real(crm_rknd), allocatable :: crm_output_qc_mean(:,:)  ! mean cloud water
   real(crm_rknd), allocatable :: crm_output_qi_mean(:,:)  ! mean cloud ice
   real(crm_rknd), allocatable :: crm_output_qs_mean(:,:)  ! mean snow
   real(crm_rknd), allocatable :: crm_output_qg_mean(:,:)  ! mean graupel
   real(crm_rknd), allocatable :: crm_output_qr_mean(:,:)  ! mean rain
#ifdef m2005
   real(crm_rknd), allocatable :: crm_output_nc_mean(:,:)  ! mean cloud water  (#/kg)
   real(crm_rknd), allocatable :: crm_output_ni_mean(:,:)  ! mean cloud ice    (#/kg)
   real(crm_rknd), allocatable :: crm_output_ns_mean(:,:)  ! mean snow         (#/kg)
   real(crm_rknd), allocatable :: crm_output_ng_mean(:,:)  ! mean graupel      (#/kg)
   real(crm_rknd), allocatable :: crm_output_nr_mean(:,:)  ! mean rain         (#/kg)

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

#if defined( SPMOMTRANS )
    real(crm_rknd), allocatable :: crm_output_ultend(:,:)            ! tendency of ul
    real(crm_rknd), allocatable :: crm_output_vltend(:,:)            ! tendency of vl
#endif

#if defined( SP_ESMT )
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
    real(crm_rknd), allocatable :: crm_output_tkz          (:,:)  ! tk profile                  [m2/s]
    real(crm_rknd), allocatable :: crm_output_flux_u       (:,:)  ! x-momentum flux             [m2/s2]
    real(crm_rknd), allocatable :: crm_output_flux_v       (:,:)  ! y-momentum flux             [m2/s2]
    real(crm_rknd), allocatable :: crm_output_flux_qp      (:,:)  ! precipitating water flux    [kg/m2/s or mm/s]
    real(crm_rknd), allocatable :: crm_output_precflux     (:,:)  ! precipitation flux          [m/s]
    real(crm_rknd), allocatable :: crm_output_qt_ls        (:,:)  ! tend of nonprec water due to large-scale   [kg/kg/s]
    real(crm_rknd), allocatable :: crm_output_qt_trans     (:,:)  ! tend of nonprec water due to transport     [kg/kg/s]
    real(crm_rknd), allocatable :: crm_output_qp_trans     (:,:)  ! tend of    prec water due to transport     [kg/kg/s]
    real(crm_rknd), allocatable :: crm_output_qp_fall      (:,:)  ! tend of    prec water due to fall-out      [kg/kg/s]
    real(crm_rknd), allocatable :: crm_output_qp_src       (:,:)  ! tend of    prec water due to conversion    [kg/kg/s]
    real(crm_rknd), allocatable :: crm_output_qp_evp       (:,:)  ! tend of    prec water due to evp           [kg/kg/s]
    real(crm_rknd), allocatable :: crm_output_t_ls         (:,:)  ! tend of lwse  due to large-scale           [kg/kg/s] ???
    real(crm_rknd), allocatable :: crm_output_prectend     (:)    ! column integrated tend in precip water+ice [kg/m2/s]
    real(crm_rknd), allocatable :: crm_output_precstend    (:)    ! column integrated tend in precip ice       [kg/m2/s]
    real(crm_rknd), allocatable :: crm_output_taux     (:)    ! zonal CRM surface stress perturbation      [N/m2]
    real(crm_rknd), allocatable :: crm_output_tauy     (:)    ! merid CRM surface stress perturbation      [N/m2]
    real(crm_rknd), allocatable :: crm_output_z0m          (:)    ! surface stress                             [N/m2]
    real(crm_rknd), allocatable :: crm_output_timing_factor(:)    ! crm cpu efficiency

    public :: crm_output_initialize
    public :: crm_output_finalize
#if defined(_OPENMP)
    public :: update_host_output
#endif

contains
   !--------------------------------------------------------------------
   subroutine crm_output_initialize(ncol, nlev) 
      integer, intent(in), optional :: ncol, nlev
      ! Allocate arrays if dimensions are passed as input
      if (present(ncol)) then
         ! Allocate instantaneous outputs
         if (.not. allocated(crm_output_qcl)) allocate(crm_output_qcl(ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(crm_output_qci)) allocate(crm_output_qci(ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(crm_output_qpl)) allocate(crm_output_qpl(ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(crm_output_qpi)) allocate(crm_output_qpi(ncol,crm_nx,crm_ny,crm_nz))

         if (.not. allocated(crm_output_tk )) allocate(crm_output_tk (ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(crm_output_tkh)) allocate(crm_output_tkh(ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(crm_output_prec_crm)) allocate(crm_output_prec_crm(ncol,crm_nx,crm_ny))

         if (.not. allocated(crm_output_wvar)) allocate(crm_output_wvar(ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(crm_output_aut))  allocate(crm_output_aut (ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(crm_output_acc))  allocate(crm_output_acc (ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(crm_output_evpc)) allocate(crm_output_evpc(ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(crm_output_evpr)) allocate(crm_output_evpr(ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(crm_output_mlt))  allocate(crm_output_mlt (ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(crm_output_sub))  allocate(crm_output_sub (ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(crm_output_dep))  allocate(crm_output_dep (ncol,crm_nx,crm_ny,crm_nz))
         if (.not. allocated(crm_output_con))  allocate(crm_output_con (ncol,crm_nx,crm_ny,crm_nz))


         ! Allocate domain and time-averaged fields
         if (.not. allocated(crm_output_cltot)) allocate(crm_output_cltot(ncol))
         if (.not. allocated(crm_output_cllow)) allocate(crm_output_cllow(ncol))
         if (.not. allocated(crm_output_clmed)) allocate(crm_output_clmed(ncol))
         if (.not. allocated(crm_output_clhgh)) allocate(crm_output_clhgh(ncol))

         if (.not. allocated(crm_output_precc))  allocate(crm_output_precc(ncol))
         if (.not. allocated(crm_output_precl))  allocate(crm_output_precl(ncol))
         if (.not. allocated(crm_output_precsc)) allocate(crm_output_precsc(ncol))
         if (.not. allocated(crm_output_precsl)) allocate(crm_output_precsl(ncol))

         ! NOTE: this output had a bug in the previous implementation
         if (.not. allocated(crm_output_cldtop)) allocate(crm_output_cldtop(ncol,nlev))

         if (.not. allocated(crm_output_qc_mean)) allocate(crm_output_qc_mean(ncol,nlev))
         if (.not. allocated(crm_output_qi_mean)) allocate(crm_output_qi_mean(ncol,nlev))
         if (.not. allocated(crm_output_qs_mean)) allocate(crm_output_qs_mean(ncol,nlev))
         if (.not. allocated(crm_output_qg_mean)) allocate(crm_output_qg_mean(ncol,nlev))
         if (.not. allocated(crm_output_qr_mean)) allocate(crm_output_qr_mean(ncol,nlev))
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
#elif defined(_OPENMP)
         !$omp target enter data map(alloc:crm_output_qcl)
         !$omp target enter data map(alloc:crm_output_qci)
         !$omp target enter data map(alloc:crm_output_qpl)
         !$omp target enter data map(alloc:crm_output_qpi)
         !$omp target enter data map(alloc:crm_output_tk)
         !$omp target enter data map(alloc:crm_output_tkh)
         !$omp target enter data map(alloc:crm_output_prec_crm)
         !$omp target enter data map(alloc:crm_output_wvar)
         !$omp target enter data map(alloc:crm_output_aut)
         !$omp target enter data map(alloc:crm_output_acc)
         !$omp target enter data map(alloc:crm_output_evpc)
         !$omp target enter data map(alloc:crm_output_evpr)
         !$omp target enter data map(alloc:crm_output_mlt)
         !$omp target enter data map(alloc:crm_output_sub)
         !$omp target enter data map(alloc:crm_output_dep)
         !$omp target enter data map(alloc:crm_output_con)
         !$omp target enter data map(alloc:crm_output_cltot)
         !$omp target enter data map(alloc:crm_output_cllow)
         !$omp target enter data map(alloc:crm_output_clmed)
         !$omp target enter data map(alloc:crm_output_clhgh)
         !$omp target enter data map(alloc:crm_output_precc)
         !$omp target enter data map(alloc:crm_output_precl)
         !$omp target enter data map(alloc:crm_output_precsc)
         !$omp target enter data map(alloc:crm_output_precsl)
         !$omp target enter data map(alloc:crm_output_cldtop)
         !$omp target enter data map(alloc:crm_output_qc_mean)
         !$omp target enter data map(alloc:crm_output_qi_mean)
         !$omp target enter data map(alloc:crm_output_qs_mean)
         !$omp target enter data map(alloc:crm_output_qg_mean)
         !$omp target enter data map(alloc:crm_output_qr_mean)
#endif
#ifdef m2005
         if (.not. allocated(crm_output_nc_mean)) allocate(crm_output_nc_mean(ncol,nlev))
         if (.not. allocated(crm_output_ni_mean)) allocate(crm_output_ni_mean(ncol,nlev))
         if (.not. allocated(crm_output_ns_mean)) allocate(crm_output_ns_mean(ncol,nlev))
         if (.not. allocated(crm_output_ng_mean)) allocate(crm_output_ng_mean(ncol,nlev))
         if (.not. allocated(crm_output_nr_mean)) allocate(crm_output_nr_mean(ncol,nlev))

         if (.not. allocated(crm_output_aut_a )) allocate(crm_output_aut_a (ncol,nlev))
         if (.not. allocated(crm_output_acc_a )) allocate(crm_output_acc_a (ncol,nlev))
         if (.not. allocated(crm_output_evpc_a)) allocate(crm_output_evpc_a(ncol,nlev))
         if (.not. allocated(crm_output_evpr_a)) allocate(crm_output_evpr_a(ncol,nlev))
         if (.not. allocated(crm_output_mlt_a )) allocate(crm_output_mlt_a (ncol,nlev))
         if (.not. allocated(crm_output_sub_a )) allocate(crm_output_sub_a (ncol,nlev))
         if (.not. allocated(crm_output_dep_a )) allocate(crm_output_dep_a (ncol,nlev))
         if (.not. allocated(crm_output_con_a )) allocate(crm_output_con_a (ncol,nlev))
#endif /* m2005 */
#if defined( SPMOMTRANS )
         if (.not. allocated(crm_output_ultend )) allocate(crm_output_ultend (ncol,nlev))
         if (.not. allocated(crm_output_vltend )) allocate(crm_output_vltend (ncol,nlev))
#endif
#if defined( SP_ESMT )
         if (.not. allocated(crm_output_u_tend_esmt )) allocate(crm_output_u_tend_esmt (ncol,nlev))
         if (.not. allocated(crm_output_v_tend_esmt )) allocate(crm_output_v_tend_esmt (ncol,nlev))
#endif   
         if (.not. allocated(crm_output_sltend ))  allocate(crm_output_sltend (ncol,nlev))
         if (.not. allocated(crm_output_qltend ))  allocate(crm_output_qltend (ncol,nlev))
         if (.not. allocated(crm_output_qcltend))  allocate(crm_output_qcltend(ncol,nlev))
         if (.not. allocated(crm_output_qiltend))  allocate(crm_output_qiltend(ncol,nlev))

         if (.not. allocated(crm_output_cld   )) allocate(crm_output_cld   (ncol,nlev))  ! cloud fraction
         if (.not. allocated(crm_output_gicewp)) allocate(crm_output_gicewp(ncol,nlev))  ! ice water path
         if (.not. allocated(crm_output_gliqwp)) allocate(crm_output_gliqwp(ncol,nlev))  ! ice water path
         if (.not. allocated(crm_output_mctot )) allocate(crm_output_mctot (ncol,nlev))  ! cloud mass flux
         if (.not. allocated(crm_output_mcup  )) allocate(crm_output_mcup  (ncol,nlev))  ! updraft cloud mass flux
         if (.not. allocated(crm_output_mcdn  )) allocate(crm_output_mcdn  (ncol,nlev))  ! downdraft cloud mass flux
         if (.not. allocated(crm_output_mcuup )) allocate(crm_output_mcuup (ncol,nlev))  ! unsat updraft cloud mass flux
         if (.not. allocated(crm_output_mcudn )) allocate(crm_output_mcudn (ncol,nlev))  ! unsat downdraft cloud mass flux

         if (.not. allocated(crm_output_mu_crm)) allocate(crm_output_mu_crm(ncol,nlev))  ! mass flux up
         if (.not. allocated(crm_output_md_crm)) allocate(crm_output_md_crm(ncol,nlev))  ! mass flux down
         if (.not. allocated(crm_output_du_crm)) allocate(crm_output_du_crm(ncol,nlev))  ! mass detrainment from updraft
         if (.not. allocated(crm_output_eu_crm)) allocate(crm_output_eu_crm(ncol,nlev))  ! mass entrainment from updraft
         if (.not. allocated(crm_output_ed_crm)) allocate(crm_output_ed_crm(ncol,nlev))  ! mass detrainment from downdraft
         if (.not. allocated(crm_output_jt_crm)) allocate(crm_output_jt_crm(ncol))       ! index of cloud (convection) top
         if (.not. allocated(crm_output_mx_crm)) allocate(crm_output_mx_crm(ncol))       ! index of cloud (convection) bottom

         if (.not. allocated(crm_output_flux_qt      )) allocate(crm_output_flux_qt      (ncol,nlev))
         if (.not. allocated(crm_output_fluxsgs_qt   )) allocate(crm_output_fluxsgs_qt   (ncol,nlev))
         if (.not. allocated(crm_output_tkez         )) allocate(crm_output_tkez         (ncol,nlev))
         if (.not. allocated(crm_output_tkesgsz      )) allocate(crm_output_tkesgsz      (ncol,nlev))
         if (.not. allocated(crm_output_tkz          )) allocate(crm_output_tkz          (ncol,nlev))
         if (.not. allocated(crm_output_flux_u       )) allocate(crm_output_flux_u       (ncol,nlev))
         if (.not. allocated(crm_output_flux_v       )) allocate(crm_output_flux_v       (ncol,nlev))
         if (.not. allocated(crm_output_flux_qp      )) allocate(crm_output_flux_qp      (ncol,nlev))
         if (.not. allocated(crm_output_precflux     )) allocate(crm_output_precflux     (ncol,nlev))
         if (.not. allocated(crm_output_qt_ls        )) allocate(crm_output_qt_ls        (ncol,nlev))
         if (.not. allocated(crm_output_qt_trans     )) allocate(crm_output_qt_trans     (ncol,nlev))
         if (.not. allocated(crm_output_qp_trans     )) allocate(crm_output_qp_trans     (ncol,nlev))
         if (.not. allocated(crm_output_qp_fall      )) allocate(crm_output_qp_fall      (ncol,nlev))
         if (.not. allocated(crm_output_qp_src       )) allocate(crm_output_qp_src       (ncol,nlev))
         if (.not. allocated(crm_output_qp_evp       )) allocate(crm_output_qp_evp       (ncol,nlev))
         if (.not. allocated(crm_output_t_ls         )) allocate(crm_output_t_ls         (ncol,nlev))
         if (.not. allocated(crm_output_prectend     )) allocate(crm_output_prectend     (ncol))
         if (.not. allocated(crm_output_precstend    )) allocate(crm_output_precstend    (ncol))
         if (.not. allocated(crm_output_taux         )) allocate(crm_output_taux         (ncol))
         if (.not. allocated(crm_output_tauy         )) allocate(crm_output_tauy         (ncol))
         if (.not. allocated(crm_output_z0m          )) allocate(crm_output_z0m          (ncol))
         if (.not. allocated(crm_output_timing_factor)) allocate(crm_output_timing_factor(ncol))
#if defined(_OPENACC)
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
         call prefetch(crm_output_timing_factor )
#elif defined(_OPENMP)
         !$omp target enter data map(alloc: crm_output_sltend)
         !$omp target enter data map(alloc: crm_output_qltend)
         !$omp target enter data map(alloc: crm_output_qcltend)
         !$omp target enter data map(alloc: crm_output_qiltend)
         !$omp target enter data map(alloc: crm_output_cld)
         !$omp target enter data map(alloc: crm_output_gicewp)
         !$omp target enter data map(alloc: crm_output_gliqwp)
         !$omp target enter data map(alloc: crm_output_mctot)
         !$omp target enter data map(alloc: crm_output_mcup)
         !$omp target enter data map(alloc: crm_output_mcdn)
         !$omp target enter data map(alloc: crm_output_mcuup)
         !$omp target enter data map(alloc: crm_output_mcudn)
         !$omp target enter data map(alloc: crm_output_mu_crm)
         !$omp target enter data map(alloc: crm_output_md_crm)
         !$omp target enter data map(alloc: crm_output_du_crm)
         !$omp target enter data map(alloc: crm_output_eu_crm)
         !$omp target enter data map(alloc: crm_output_ed_crm)
         !$omp target enter data map(alloc: crm_output_jt_crm)
         !$omp target enter data map(alloc: crm_output_mx_crm)
         !$omp target enter data map(alloc: crm_output_flux_qt)
         !$omp target enter data map(alloc: crm_output_fluxsgs_qt)
         !$omp target enter data map(alloc: crm_output_tkez)
         !$omp target enter data map(alloc: crm_output_tkesgsz)
         !$omp target enter data map(alloc: crm_output_tkz)
         !$omp target enter data map(alloc: crm_output_flux_u)
         !$omp target enter data map(alloc: crm_output_flux_v)
         !$omp target enter data map(alloc: crm_output_flux_qp)
         !$omp target enter data map(alloc: crm_output_precflux)
         !$omp target enter data map(alloc: crm_output_qt_ls)
         !$omp target enter data map(alloc: crm_output_qt_trans)
         !$omp target enter data map(alloc: crm_output_qp_trans)
         !$omp target enter data map(alloc: crm_output_qp_fall)
         !$omp target enter data map(alloc: crm_output_qp_src)
         !$omp target enter data map(alloc: crm_output_qp_evp)
         !$omp target enter data map(alloc: crm_output_t_ls)
         !$omp target enter data map(alloc: crm_output_prectend)
         !$omp target enter data map(alloc: crm_output_precstend)
         !$omp target enter data map(alloc: crm_output_taux)
         !$omp target enter data map(alloc: crm_output_tauy)
         !$omp target enter data map(alloc: crm_output_z0m)
         !$omp target enter data map(alloc: crm_output_timing_factor)
#endif
      end if ! present(ncol)

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

#if defined( SPMOMTRANS )
      crm_output_ultend = 0
      crm_output_vltend = 0
#endif

#if defined( SP_ESMT )
      crm_output_u_tend_esmt = 0
      crm_output_v_tend_esmt = 0
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
      crm_output_taux          = 0
      crm_output_tauy          = 0
      crm_output_z0m           = 0
      crm_output_timing_factor = 0
   end subroutine crm_output_initialize
   !\
   ! update host data from device
   !/
#if defined(_OPENMP)
   subroutine update_host_output() 
      ! Allocate arrays if dimensions are passed as input
      !$omp target update from (crm_output_qcl)
      !$omp target update from (crm_output_qci)
      !$omp target update from (crm_output_qpl)
      !$omp target update from (crm_output_qpi)
      !$omp target update from (crm_output_tk)
      !$omp target update from (crm_output_tkh)
      !$omp target update from (crm_output_prec_crm)
      !$omp target update from (crm_output_wvar)
      !$omp target update from (crm_output_aut)
      !$omp target update from (crm_output_acc)
      !$omp target update from (crm_output_evpc)
      !$omp target update from (crm_output_evpr)
      !$omp target update from (crm_output_mlt)
      !$omp target update from (crm_output_sub)
      !$omp target update from (crm_output_dep)
      !$omp target update from (crm_output_con)
      !$omp target update from (crm_output_cltot)
      !$omp target update from (crm_output_cllow)
      !$omp target update from (crm_output_clmed)
      !$omp target update from (crm_output_clhgh)
      !$omp target update from (crm_output_precc)
      !$omp target update from (crm_output_precl)
      !$omp target update from (crm_output_precsc)
      !$omp target update from (crm_output_precsl)
      !$omp target update from (crm_output_cldtop)
      !$omp target update from (crm_output_qc_mean)
      !$omp target update from (crm_output_qi_mean)
      !$omp target update from (crm_output_qs_mean)
      !$omp target update from (crm_output_qg_mean)
      !$omp target update from (crm_output_qr_mean)
      !$omp target update from(crm_output_sltend)
      !$omp target update from(crm_output_qltend)
      !$omp target update from(crm_output_qcltend)
      !$omp target update from(crm_output_qiltend)
      !$omp target update from(crm_output_cld)
      !$omp target update from(crm_output_gicewp)
      !$omp target update from(crm_output_gliqwp)
      !$omp target update from(crm_output_mctot)
      !$omp target update from(crm_output_mcup)
      !$omp target update from(crm_output_mcdn)
      !$omp target update from(crm_output_mcuup)
      !$omp target update from(crm_output_mcudn)
      !$omp target update from(crm_output_mu_crm)
      !$omp target update from(crm_output_md_crm)
      !$omp target update from(crm_output_du_crm)
      !$omp target update from(crm_output_eu_crm)
      !$omp target update from(crm_output_ed_crm)
      !$omp target update from(crm_output_jt_crm)
      !$omp target update from(crm_output_mx_crm)
      !$omp target update from(crm_output_flux_qt)
      !$omp target update from(crm_output_fluxsgs_qt)
      !$omp target update from(crm_output_tkez)
      !$omp target update from(crm_output_tkesgsz)
      !$omp target update from(crm_output_tkz)
      !$omp target update from(crm_output_flux_u)
      !$omp target update from(crm_output_flux_v)
      !$omp target update from(crm_output_flux_qp)
      !$omp target update from(crm_output_precflux)
      !$omp target update from(crm_output_qt_ls)
      !$omp target update from(crm_output_qt_trans)
      !$omp target update from(crm_output_qp_trans) 
      !$omp target update from(crm_output_qp_fall)
      !$omp target update from(crm_output_qp_src)
      !$omp target update from(crm_output_qp_evp)
      !$omp target update from(crm_output_t_ls)
      !$omp target update from(crm_output_prectend)
      !$omp target update from(crm_output_precstend)
      !$omp target update from(crm_output_taux)
      !$omp target update from(crm_output_tauy)
      !$omp target update from(crm_output_z0m)
      !$omp target update from(crm_output_timing_factor)
   end subroutine update_host_output
#endif
   !\
   ! finalize the device data, clean up the map
   !/
   subroutine crm_output_finalize()
#if defined(_OPENMP)
      !$omp target exit data map(delete: crm_output_qcl)
      !$omp target exit data map(delete: crm_output_qci)
      !$omp target exit data map(delete: crm_output_qpl)
      !$omp target exit data map(delete: crm_output_qpi)
      !$omp target exit data map(delete: crm_output_tk )
      !$omp target exit data map(delete: crm_output_tkh)
      !$omp target exit data map(delete: crm_output_prec_crm)
      !$omp target exit data map(delete: crm_output_wvar)
      !$omp target exit data map(delete: crm_output_aut )
      !$omp target exit data map(delete: crm_output_acc )
      !$omp target exit data map(delete: crm_output_evpc)
      !$omp target exit data map(delete: crm_output_evpr)
      !$omp target exit data map(delete: crm_output_mlt )
      !$omp target exit data map(delete: crm_output_sub )
      !$omp target exit data map(delete: crm_output_dep )
      !$omp target exit data map(delete: crm_output_con )
      !$omp target exit data map(delete: crm_output_cltot)
      !$omp target exit data map(delete: crm_output_cllow)
      !$omp target exit data map(delete: crm_output_clmed)
      !$omp target exit data map(delete: crm_output_clhgh)
      !$omp target exit data map(delete: crm_output_precc)
      !$omp target exit data map(delete: crm_output_precl)
      !$omp target exit data map(delete: crm_output_precsc)
      !$omp target exit data map(delete: crm_output_precsl)
      !$omp target exit data map(delete: crm_output_cldtop)
      !$omp target exit data map(delete: crm_output_qc_mean)
      !$omp target exit data map(delete: crm_output_qi_mean)
      !$omp target exit data map(delete: crm_output_qs_mean)
      !$omp target exit data map(delete: crm_output_qg_mean)
      !$omp target exit data map(delete: crm_output_qr_mean)
      !$omp target exit data map(delete: crm_output_sltend  )
      !$omp target exit data map(delete: crm_output_qltend  )
      !$omp target exit data map(delete: crm_output_qcltend )
      !$omp target exit data map(delete: crm_output_qiltend )
      !$omp target exit data map(delete: crm_output_cld    )
      !$omp target exit data map(delete: crm_output_gicewp )
      !$omp target exit data map(delete: crm_output_gliqwp )
      !$omp target exit data map(delete: crm_output_mctot  )
      !$omp target exit data map(delete: crm_output_mcup   )
      !$omp target exit data map(delete: crm_output_mcdn   )
      !$omp target exit data map(delete: crm_output_mcuup  )
      !$omp target exit data map(delete: crm_output_mcudn  )
      !$omp target exit data map(delete: crm_output_mu_crm )
      !$omp target exit data map(delete: crm_output_md_crm )
      !$omp target exit data map(delete: crm_output_du_crm )
      !$omp target exit data map(delete: crm_output_eu_crm )
      !$omp target exit data map(delete: crm_output_ed_crm )
      !$omp target exit data map(delete: crm_output_jt_crm )
      !$omp target exit data map(delete: crm_output_mx_crm )
      !$omp target exit data map(delete: crm_output_flux_qt       )
      !$omp target exit data map(delete: crm_output_fluxsgs_qt    )
      !$omp target exit data map(delete: crm_output_tkez          )
      !$omp target exit data map(delete: crm_output_tkesgsz       )
      !$omp target exit data map(delete: crm_output_tkz           )
      !$omp target exit data map(delete: crm_output_flux_u        )
      !$omp target exit data map(delete: crm_output_flux_v        )
      !$omp target exit data map(delete: crm_output_flux_qp       )
      !$omp target exit data map(delete: crm_output_precflux      )
      !$omp target exit data map(delete: crm_output_qt_ls         )
      !$omp target exit data map(delete: crm_output_qt_trans      )
      !$omp target exit data map(delete: crm_output_qp_trans      )
      !$omp target exit data map(delete: crm_output_qp_fall       )
      !$omp target exit data map(delete: crm_output_qp_src        )
      !$omp target exit data map(delete: crm_output_qp_evp        )
      !$omp target exit data map(delete: crm_output_t_ls          )
      !$omp target exit data map(delete: crm_output_prectend      )
      !$omp target exit data map(delete: crm_output_precstend     )
      !$omp target exit data map(delete: crm_output_taux          )
      !$omp target exit data map(delete: crm_output_tauy          )
      !$omp target exit data map(delete: crm_output_z0m           )
      !$omp target exit data map(delete: crm_output_timing_factor )
#endif
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

#if defined( SPMOMTRANS )
      if (allocated(crm_output_ultend)) deallocate(crm_output_ultend)
      if (allocated(crm_output_vltend)) deallocate(crm_output_vltend)
#endif

#if defined( SP_ESMT )
      if (allocated(crm_output_u_tend_esmt)) deallocate(crm_output_u_tend_esmt)
      if (allocated(crm_output_v_tend_esmt)) deallocate(crm_output_v_tend_esmt)
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
   end subroutine crm_output_finalize
end module crm_output_module
