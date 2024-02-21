module crm_output_module
   use params_kind,  only: crm_rknd
   use openacc_utils
   implicit none
   public crm_output_type
   type crm_output_type
      ! Derived type to encapsulate CRM output fields (things that are output
      ! only, intent(out) in the original implementation)

      ! These are copies of the SAM cloud and precip liquid and ice, previously
      ! passed in and out of crm() via qc_crm, qi_crm, etc.
      real(crm_rknd), allocatable :: qcl(:,:,:,:)
      real(crm_rknd), allocatable :: qci(:,:,:,:)
      real(crm_rknd), allocatable :: qpl(:,:,:,:)
      real(crm_rknd), allocatable :: qpi(:,:,:,:)

      real(crm_rknd), allocatable :: tk (:,:,:,:)
      real(crm_rknd), allocatable :: tkh(:,:,:,:)
      real(crm_rknd), allocatable :: prec_crm(:,:,:) ! CRM precipiation rate (surface)

      ! 2-moment process rates
      real(crm_rknd), allocatable :: wvar(:,:,:,:) ! vertical velocity variance (m/s)
      real(crm_rknd), allocatable :: aut (:,:,:,:) ! cloud water autoconversion (1/s)
      real(crm_rknd), allocatable :: acc (:,:,:,:) ! cloud water accretion (1/s)
      real(crm_rknd), allocatable :: evpc(:,:,:,:) ! cloud water evaporation (1/s)
      real(crm_rknd), allocatable :: evpr(:,:,:,:) ! rain evaporation (1/s)
      real(crm_rknd), allocatable :: mlt (:,:,:,:) ! ice, snow, graupel melting (1/s)
      real(crm_rknd), allocatable :: sub (:,:,:,:) ! ice, snow, graupel sublimation (1/s)
      real(crm_rknd), allocatable :: dep (:,:,:,:) ! ice, snow, graupel deposition (1/s)
      real(crm_rknd), allocatable :: con (:,:,:,:) ! cloud water condensation(1/s)

      ! Cloud area fractions
      real(crm_rknd), allocatable :: cltot(:)  ! shaded cloud fraction
      real(crm_rknd), allocatable :: clhgh(:)  ! shaded cloud fraction
      real(crm_rknd), allocatable :: clmed(:)  ! shaded cloud fraction
      real(crm_rknd), allocatable :: cllow(:)  ! shaded cloud fraction

      real(crm_rknd), allocatable :: cldtop(:,:)  ! cloud top ... pressure???
      real(crm_rknd), allocatable :: precc(:)   ! convective precipitation rate
      real(crm_rknd), allocatable :: precl(:)   ! stratiform precipitation rate
      real(crm_rknd), allocatable :: precsc(:)   ! convective snow precipitation rate
      real(crm_rknd), allocatable :: precsl(:)   ! stratiform snow precipitation rate

      ! TODO: These diagnostics are currently all on the GCM vertical grid. I think this is
      ! misleading though, and overly complicates crm_module. I think the better thing to do would
      ! be to define everything within crm_module on the CRM grid, and then do the
      ! mapping/interpolation at the GCM (crm_physics) level. For now though, I am just copying
      ! these over directly to minimize chances of me making a mistake. Also, some of these probably
      ! do not need to be calculated here, and might make more sense to calculate at the
      ! crm_physics_tend level. For example, I think tendencies should be calculated in
      ! crm_physics_tend, from, for example, something like crm_output%uwind - crm_input%uwind.
      real(crm_rknd), allocatable :: qv_mean(:,:)  ! mean cloud water
      real(crm_rknd), allocatable :: qc_mean(:,:)  ! mean cloud water
      real(crm_rknd), allocatable :: qi_mean(:,:)  ! mean cloud ice
      real(crm_rknd), allocatable :: qr_mean(:,:)  ! mean rain
      real(crm_rknd), allocatable :: qs_mean(:,:)  ! mean snow
      real(crm_rknd), allocatable :: qg_mean(:,:)  ! mean graupel
      real(crm_rknd), allocatable :: qm_mean(:,:)  ! mean ice rime mass
      real(crm_rknd), allocatable :: bm_mean(:,:)  ! mean ice rime volume
      real(crm_rknd), allocatable :: rho_d_mean(:,:)  ! mean dry density
      real(crm_rknd), allocatable :: rho_v_mean(:,:)  ! mean vapor density

      real(crm_rknd), allocatable :: nc_mean(:,:)  ! mean cloud water  (#/kg)
      real(crm_rknd), allocatable :: ni_mean(:,:)  ! mean cloud ice    (#/kg)
      real(crm_rknd), allocatable :: nr_mean(:,:)  ! mean rain         (#/kg)

      real(crm_rknd), allocatable :: ultend  (:,:)          ! CRM output tendency of zonal wind
      real(crm_rknd), allocatable :: vltend  (:,:)          ! CRM output tendency of meridional wind
      real(crm_rknd), allocatable :: sltend  (:,:)          ! CRM output tendency of static energy
      real(crm_rknd), allocatable :: qltend  (:,:)          ! CRM output tendency of water vapor
      real(crm_rknd), allocatable :: qcltend (:,:)          ! CRM output tendency of cloud liquid water
      real(crm_rknd), allocatable :: qiltend (:,:)          ! CRM output tendency of cloud ice

      real(crm_rknd), allocatable :: t_vt_tend (:,:)       ! CRM output tendency for LSE variance transport
      real(crm_rknd), allocatable :: q_vt_tend (:,:)       ! CRM output tendency for QT  variance transport
      real(crm_rknd), allocatable :: u_vt_tend (:,:)       ! CRM output tendency for U variance transport
      real(crm_rknd), allocatable :: t_vt_ls   (:,:)       ! large-scale LSE variance transport tendency from GCM
      real(crm_rknd), allocatable :: q_vt_ls   (:,:)       ! large-scale QT  variance transport tendency from GCM
      real(crm_rknd), allocatable :: u_vt_ls   (:,:)       ! large-scale U variance transport tendency from GCM

      ! These are all time and spatial averages, on the GCM grid
      real(crm_rknd), allocatable :: cld   (:,:)      ! cloud fraction
      real(crm_rknd), allocatable :: gicewp(:,:)      ! ice water path
      real(crm_rknd), allocatable :: gliqwp(:,:)      ! ice water path

      real(crm_rknd), allocatable :: liq_ice_exchange(:,:) ! P3 liq-ice phase change tendency
      real(crm_rknd), allocatable :: vap_liq_exchange(:,:) ! P3 vap-liq phase change tendency
      real(crm_rknd), allocatable :: vap_ice_exchange(:,:) ! P3 vap-ice phase change tendency

      real(crm_rknd), allocatable :: mctot (:,:)      ! cloud mass flux
      real(crm_rknd), allocatable :: mcup  (:,:)      ! updraft cloud mass flux
      real(crm_rknd), allocatable :: mcdn  (:,:)      ! downdraft cloud mass flux
      real(crm_rknd), allocatable :: mcuup (:,:)      ! unsat updraft cloud mass flux
      real(crm_rknd), allocatable :: mcudn (:,:)      ! unsat downdraft cloud mass flux

      ! For convective transport
      real(crm_rknd), allocatable :: mu_crm(:,:)      ! mass flux up
      real(crm_rknd), allocatable :: md_crm(:,:)      ! mass flux down
      real(crm_rknd), allocatable :: du_crm(:,:)      ! mass detrainment from updraft
      real(crm_rknd), allocatable :: eu_crm(:,:)      ! mass entrainment from updraft
      real(crm_rknd), allocatable :: ed_crm(:,:)      ! mass detrainment from downdraft
      real(crm_rknd), allocatable :: jt_crm(:)        ! index of cloud (convection) top
      real(crm_rknd), allocatable :: mx_crm(:)        ! index of cloud (convection) bottom

      ! Other stuff...
      real(crm_rknd), allocatable :: flux_qt      (:,:)  ! nonprecip water flux        [kg/m2/s]
      real(crm_rknd), allocatable :: fluxsgs_qt   (:,:)  ! sgs non-precip water flux   [kg/m2/s]
      real(crm_rknd), allocatable :: tkez         (:,:)  ! tke profile                 [kg/m/s2]
      real(crm_rknd), allocatable :: tkew         (:,:)  ! vertical velocity variance  [kg/m/s2]
      real(crm_rknd), allocatable :: tkesgsz      (:,:)  ! sgs tke profile             [kg/m/s2]
      real(crm_rknd), allocatable :: tkz          (:,:)  ! tk profile                  [m2/s]
      real(crm_rknd), allocatable :: flux_u       (:,:)  ! x-momentum flux             [m2/s2]
      real(crm_rknd), allocatable :: flux_v       (:,:)  ! y-momentum flux             [m2/s2]
      real(crm_rknd), allocatable :: flux_qp      (:,:)  ! precipitating water flux    [kg/m2/s or mm/s]
      real(crm_rknd), allocatable :: precflux     (:,:)  ! precipitation flux          [m/s]
      real(crm_rknd), allocatable :: qt_ls        (:,:)  ! tend of nonprec water due to large-scale   [kg/kg/s]
      real(crm_rknd), allocatable :: qt_trans     (:,:)  ! tend of nonprec water due to transport     [kg/kg/s]
      real(crm_rknd), allocatable :: qp_trans     (:,:)  ! tend of    prec water due to transport     [kg/kg/s]
      real(crm_rknd), allocatable :: qp_fall      (:,:)  ! tend of    prec water due to fall-out      [kg/kg/s]
      real(crm_rknd), allocatable :: qp_src       (:,:)  ! tend of    prec water due to conversion    [kg/kg/s]
      real(crm_rknd), allocatable :: qp_evp       (:,:)  ! tend of    prec water due to evp           [kg/kg/s]
      real(crm_rknd), allocatable :: t_ls         (:,:)  ! tend of lwse  due to large-scale           [kg/kg/s] ???
      real(crm_rknd), allocatable :: prectend     (:)    ! column integrated tend in precip water+ice [kg/m2/s]
      real(crm_rknd), allocatable :: precstend    (:)    ! column integrated tend in precip ice       [kg/m2/s]
      real(crm_rknd), allocatable :: taux         (:)    ! zonal CRM surface stress perturbation      [N/m2]
      real(crm_rknd), allocatable :: tauy         (:)    ! merid CRM surface stress perturbation      [N/m2]
      real(crm_rknd), allocatable :: z0m          (:)    ! surface stress                             [N/m2]
      real(crm_rknd), allocatable :: subcycle_factor(:)    ! crm cpu efficiency

      real(crm_rknd), allocatable :: dt_sgs       (:,:)  ! CRM temperature tendency from SGS   [K/s]
      real(crm_rknd), allocatable :: dqv_sgs      (:,:)  ! CRM water vapor tendency from SGS   [kg/kg/s]
      real(crm_rknd), allocatable :: dqc_sgs      (:,:)  ! CRM cloud water tendency from SGS   [kg/kg/s]
      real(crm_rknd), allocatable :: dqi_sgs      (:,:)  ! CRM cloud ice tendency from SGS     [kg/kg/s]
      real(crm_rknd), allocatable :: dqr_sgs      (:,:)  ! CRM liquid rain tendency from SGS   [kg/kg/s]

      real(crm_rknd), allocatable :: dt_micro     (:,:)  ! CRM temperature tendency from micro [K/s]
      real(crm_rknd), allocatable :: dqv_micro    (:,:)  ! CRM water vapor tendency from micro [kg/kg/s]
      real(crm_rknd), allocatable :: dqc_micro    (:,:)  ! CRM cloud water tendency from micro [kg/kg/s]
      real(crm_rknd), allocatable :: dqi_micro    (:,:)  ! CRM cloud ice tendency from micro   [kg/kg/s]
      real(crm_rknd), allocatable :: dqr_micro    (:,:)  ! CRM liquid rain tendency from micro  [kg/kg/s]

      real(crm_rknd), allocatable :: dt_dycor  (:,:)
      real(crm_rknd), allocatable :: dqv_dycor (:,:)
      real(crm_rknd), allocatable :: dqc_dycor (:,:)
      real(crm_rknd), allocatable :: dqi_dycor (:,:)
      real(crm_rknd), allocatable :: dqr_dycor (:,:)

      real(crm_rknd), allocatable :: dt_sponge (:,:)
      real(crm_rknd), allocatable :: dqv_sponge(:,:)
      real(crm_rknd), allocatable :: dqc_sponge(:,:)
      real(crm_rknd), allocatable :: dqi_sponge(:,:)
      real(crm_rknd), allocatable :: dqr_sponge(:,:)

      real(crm_rknd), allocatable :: rho_d_ls     (:,:)  ! large-scale forcing of dry density   [kg/m3/s]
      real(crm_rknd), allocatable :: rho_v_ls     (:,:)  ! large-scale forcing of vapor density [kg/m3/s]
      real(crm_rknd), allocatable :: rho_l_ls     (:,:)  ! large-scale forcing of vapor density [kg/m3/s]
      real(crm_rknd), allocatable :: rho_i_ls     (:,:)  ! large-scale forcing of vapor density [kg/m3/s]

   end type crm_output_type

contains

   !------------------------------------------------------------------------------------------------
   subroutine crm_output_initialize(output, ncol, nlev, crm_nx, crm_ny, crm_nz, MMF_microphysics_scheme)
      type(crm_output_type), intent(inout) :: output
      integer,               intent(in   ) :: ncol, nlev, crm_nx, crm_ny, crm_nz
      character(len=*),      intent(in   ) :: MMF_microphysics_scheme    ! CRM microphysics scheme

      ! Allocate instantaneous outputs
      if (.not. allocated(output%qcl)) allocate(output%qcl(ncol,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(output%qci)) allocate(output%qci(ncol,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(output%qpl)) allocate(output%qpl(ncol,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(output%qpi)) allocate(output%qpi(ncol,crm_nx,crm_ny,crm_nz))

      if (.not. allocated(output%tk )) allocate(output%tk (ncol,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(output%tkh)) allocate(output%tkh(ncol,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(output%prec_crm)) allocate(output%prec_crm(ncol,crm_nx,crm_ny))

      if (.not. allocated(output%wvar)) allocate(output%wvar(ncol,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(output%aut))  allocate(output%aut (ncol,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(output%acc))  allocate(output%acc (ncol,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(output%evpc)) allocate(output%evpc(ncol,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(output%evpr)) allocate(output%evpr(ncol,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(output%mlt))  allocate(output%mlt (ncol,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(output%sub))  allocate(output%sub (ncol,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(output%dep))  allocate(output%dep (ncol,crm_nx,crm_ny,crm_nz))
      if (.not. allocated(output%con))  allocate(output%con (ncol,crm_nx,crm_ny,crm_nz))


      ! Allocate domain and time-averaged fields
      if (.not. allocated(output%cltot)) allocate(output%cltot(ncol))
      if (.not. allocated(output%cllow)) allocate(output%cllow(ncol))
      if (.not. allocated(output%clmed)) allocate(output%clmed(ncol))
      if (.not. allocated(output%clhgh)) allocate(output%clhgh(ncol))

      if (.not. allocated(output%precc))  allocate(output%precc(ncol))
      if (.not. allocated(output%precl))  allocate(output%precl(ncol))
      if (.not. allocated(output%precsc)) allocate(output%precsc(ncol))
      if (.not. allocated(output%precsl)) allocate(output%precsl(ncol))

      ! NOTE: this output had a bug in the previous implementation
      if (.not. allocated(output%cldtop)) allocate(output%cldtop(ncol,nlev))

      if (.not. allocated(output%qv_mean)) allocate(output%qv_mean(ncol,nlev))
      if (.not. allocated(output%qc_mean)) allocate(output%qc_mean(ncol,nlev))
      if (.not. allocated(output%qi_mean)) allocate(output%qi_mean(ncol,nlev))
      if (.not. allocated(output%qr_mean)) allocate(output%qr_mean(ncol,nlev))
      if (.not. allocated(output%qs_mean)) allocate(output%qs_mean(ncol,nlev))
      if (.not. allocated(output%qg_mean)) allocate(output%qg_mean(ncol,nlev))
      if (.not. allocated(output%qm_mean)) allocate(output%qm_mean(ncol,nlev))
      if (.not. allocated(output%bm_mean)) allocate(output%bm_mean(ncol,nlev))
      if (.not. allocated(output%rho_d_mean)) allocate(output%rho_d_mean(ncol,nlev))
      if (.not. allocated(output%rho_v_mean)) allocate(output%rho_v_mean(ncol,nlev))

      call prefetch(output%qcl)
      call prefetch(output%qci)
      call prefetch(output%qpl)
      call prefetch(output%qpi)
      call prefetch(output%tk )
      call prefetch(output%tkh)
      call prefetch(output%prec_crm)
      call prefetch(output%wvar)
      call prefetch(output%aut )
      call prefetch(output%acc )
      call prefetch(output%evpc)
      call prefetch(output%evpr)
      call prefetch(output%mlt )
      call prefetch(output%sub )
      call prefetch(output%dep )
      call prefetch(output%con )
      call prefetch(output%cltot)
      call prefetch(output%cllow)
      call prefetch(output%clmed)
      call prefetch(output%clhgh)
      call prefetch(output%precc)
      call prefetch(output%precl)
      call prefetch(output%precsc)
      call prefetch(output%precsl)
      call prefetch(output%cldtop)
      call prefetch(output%qv_mean)
      call prefetch(output%qc_mean)
      call prefetch(output%qi_mean)
      call prefetch(output%qr_mean)
      call prefetch(output%qs_mean)
      call prefetch(output%qg_mean)
      call prefetch(output%qm_mean)
      call prefetch(output%bm_mean)
      call prefetch(output%rho_d_mean)
      call prefetch(output%rho_v_mean)

      if (.not. allocated(output%nc_mean)) allocate(output%nc_mean(ncol,nlev))
      if (.not. allocated(output%ni_mean)) allocate(output%ni_mean(ncol,nlev))
      if (.not. allocated(output%nr_mean)) allocate(output%nr_mean(ncol,nlev))

      if (.not. allocated(output%ultend ))  allocate(output%ultend (ncol,nlev))
      if (.not. allocated(output%vltend ))  allocate(output%vltend (ncol,nlev))      
      if (.not. allocated(output%sltend ))  allocate(output%sltend (ncol,nlev))
      if (.not. allocated(output%qltend ))  allocate(output%qltend (ncol,nlev))
      if (.not. allocated(output%qcltend))  allocate(output%qcltend(ncol,nlev))
      if (.not. allocated(output%qiltend))  allocate(output%qiltend(ncol,nlev))

      if (.not. allocated(output%t_vt_tend))  allocate(output%t_vt_tend(ncol,nlev))
      if (.not. allocated(output%q_vt_tend))  allocate(output%q_vt_tend(ncol,nlev))
      if (.not. allocated(output%u_vt_tend))  allocate(output%u_vt_tend(ncol,nlev))
      if (.not. allocated(output%t_vt_ls  ))  allocate(output%t_vt_ls  (ncol,nlev))
      if (.not. allocated(output%q_vt_ls  ))  allocate(output%q_vt_ls  (ncol,nlev))
      if (.not. allocated(output%u_vt_ls  ))  allocate(output%u_vt_ls  (ncol,nlev))

      if (.not. allocated(output%cld   )) allocate(output%cld   (ncol,nlev))  ! cloud fraction
      if (.not. allocated(output%gicewp)) allocate(output%gicewp(ncol,nlev))  ! ice water path
      if (.not. allocated(output%gliqwp)) allocate(output%gliqwp(ncol,nlev))  ! ice water path

      if (.not. allocated(output%liq_ice_exchange)) allocate(output%liq_ice_exchange(ncol,nlev)) ! P3 liq-ice phase change tendency
      if (.not. allocated(output%vap_liq_exchange)) allocate(output%vap_liq_exchange(ncol,nlev)) ! P3 vap-liq phase change tendency
      if (.not. allocated(output%vap_ice_exchange)) allocate(output%vap_ice_exchange(ncol,nlev)) ! P3 vap-ice phase change tendency
      
      if (.not. allocated(output%mctot )) allocate(output%mctot (ncol,nlev))  ! cloud mass flux
      if (.not. allocated(output%mcup  )) allocate(output%mcup  (ncol,nlev))  ! updraft cloud mass flux
      if (.not. allocated(output%mcdn  )) allocate(output%mcdn  (ncol,nlev))  ! downdraft cloud mass flux
      if (.not. allocated(output%mcuup )) allocate(output%mcuup (ncol,nlev))  ! unsat updraft cloud mass flux
      if (.not. allocated(output%mcudn )) allocate(output%mcudn (ncol,nlev))  ! unsat downdraft cloud mass flux

      if (.not. allocated(output%mu_crm)) allocate(output%mu_crm(ncol,nlev))  ! mass flux up
      if (.not. allocated(output%md_crm)) allocate(output%md_crm(ncol,nlev))  ! mass flux down
      if (.not. allocated(output%du_crm)) allocate(output%du_crm(ncol,nlev))  ! mass detrainment from updraft
      if (.not. allocated(output%eu_crm)) allocate(output%eu_crm(ncol,nlev))  ! mass entrainment from updraft
      if (.not. allocated(output%ed_crm)) allocate(output%ed_crm(ncol,nlev))  ! mass detrainment from downdraft
      if (.not. allocated(output%jt_crm)) allocate(output%jt_crm(ncol))       ! index of cloud (convection) top
      if (.not. allocated(output%mx_crm)) allocate(output%mx_crm(ncol))       ! index of cloud (convection) bottom

      if (.not. allocated(output%flux_qt      )) allocate(output%flux_qt      (ncol,nlev))
      if (.not. allocated(output%fluxsgs_qt   )) allocate(output%fluxsgs_qt   (ncol,nlev))
      if (.not. allocated(output%tkez         )) allocate(output%tkez         (ncol,nlev))
      if (.not. allocated(output%tkew         )) allocate(output%tkew         (ncol,nlev))
      if (.not. allocated(output%tkesgsz      )) allocate(output%tkesgsz      (ncol,nlev))
      if (.not. allocated(output%tkz          )) allocate(output%tkz          (ncol,nlev))
      if (.not. allocated(output%flux_u       )) allocate(output%flux_u       (ncol,nlev))
      if (.not. allocated(output%flux_v       )) allocate(output%flux_v       (ncol,nlev))
      if (.not. allocated(output%flux_qp      )) allocate(output%flux_qp      (ncol,nlev))
      if (.not. allocated(output%precflux     )) allocate(output%precflux     (ncol,nlev))
      if (.not. allocated(output%qt_ls        )) allocate(output%qt_ls        (ncol,nlev))
      if (.not. allocated(output%qt_trans     )) allocate(output%qt_trans     (ncol,nlev))
      if (.not. allocated(output%qp_trans     )) allocate(output%qp_trans     (ncol,nlev))
      if (.not. allocated(output%qp_fall      )) allocate(output%qp_fall      (ncol,nlev))
      if (.not. allocated(output%qp_src       )) allocate(output%qp_src       (ncol,nlev))
      if (.not. allocated(output%qp_evp       )) allocate(output%qp_evp       (ncol,nlev))
      if (.not. allocated(output%t_ls         )) allocate(output%t_ls         (ncol,nlev))
      if (.not. allocated(output%prectend     )) allocate(output%prectend     (ncol))
      if (.not. allocated(output%precstend    )) allocate(output%precstend    (ncol))
      if (.not. allocated(output%taux         )) allocate(output%taux         (ncol))
      if (.not. allocated(output%tauy         )) allocate(output%tauy         (ncol))
      if (.not. allocated(output%z0m          )) allocate(output%z0m          (ncol))
      if (.not. allocated(output%subcycle_factor)) allocate(output%subcycle_factor(ncol))

      if (.not. allocated(output%dt_sgs       )) allocate(output%dt_sgs       (ncol,nlev))
      if (.not. allocated(output%dqv_sgs      )) allocate(output%dqv_sgs      (ncol,nlev))
      if (.not. allocated(output%dqc_sgs      )) allocate(output%dqc_sgs      (ncol,nlev))
      if (.not. allocated(output%dqi_sgs      )) allocate(output%dqi_sgs      (ncol,nlev))
      if (.not. allocated(output%dqr_sgs      )) allocate(output%dqr_sgs      (ncol,nlev))

      if (.not. allocated(output%dt_micro     )) allocate(output%dt_micro     (ncol,nlev))
      if (.not. allocated(output%dqv_micro    )) allocate(output%dqv_micro    (ncol,nlev))
      if (.not. allocated(output%dqc_micro    )) allocate(output%dqc_micro    (ncol,nlev))
      if (.not. allocated(output%dqi_micro    )) allocate(output%dqi_micro    (ncol,nlev))
      if (.not. allocated(output%dqr_micro    )) allocate(output%dqr_micro    (ncol,nlev))

      if (.not. allocated(output%dt_dycor     )) allocate(output%dt_dycor     (ncol,nlev))
      if (.not. allocated(output%dqv_dycor    )) allocate(output%dqv_dycor    (ncol,nlev))
      if (.not. allocated(output%dqc_dycor    )) allocate(output%dqc_dycor    (ncol,nlev))
      if (.not. allocated(output%dqi_dycor    )) allocate(output%dqi_dycor    (ncol,nlev))
      if (.not. allocated(output%dqr_dycor    )) allocate(output%dqr_dycor    (ncol,nlev))

      if (.not. allocated(output%dt_sponge    )) allocate(output%dt_sponge    (ncol,nlev))
      if (.not. allocated(output%dqv_sponge   )) allocate(output%dqv_sponge   (ncol,nlev))
      if (.not. allocated(output%dqc_sponge   )) allocate(output%dqc_sponge   (ncol,nlev))
      if (.not. allocated(output%dqi_sponge   )) allocate(output%dqi_sponge   (ncol,nlev))
      if (.not. allocated(output%dqr_sponge   )) allocate(output%dqr_sponge   (ncol,nlev))

      if (.not. allocated(output%rho_d_ls     )) allocate(output%rho_d_ls     (ncol,nlev))
      if (.not. allocated(output%rho_v_ls     )) allocate(output%rho_v_ls     (ncol,nlev))
      if (.not. allocated(output%rho_l_ls     )) allocate(output%rho_l_ls     (ncol,nlev))
      if (.not. allocated(output%rho_i_ls     )) allocate(output%rho_i_ls     (ncol,nlev))

      call prefetch(output%sltend  )
      call prefetch(output%qltend  )
      call prefetch(output%qcltend )
      call prefetch(output%qiltend )

      call prefetch(output%t_vt_tend )
      call prefetch(output%q_vt_tend )
      call prefetch(output%u_vt_tend )
      call prefetch(output%t_vt_ls   )
      call prefetch(output%q_vt_ls   )
      call prefetch(output%u_vt_ls   )

      call prefetch(output%cld    )
      call prefetch(output%gicewp )
      call prefetch(output%gliqwp )
      
      call prefetch(output%liq_ice_exchange)
      call prefetch(output%vap_liq_exchange)
      call prefetch(output%vap_ice_exchange)

      call prefetch(output%mctot  )
      call prefetch(output%mcup   )
      call prefetch(output%mcdn   )
      call prefetch(output%mcuup  )
      call prefetch(output%mcudn  )

      call prefetch(output%mu_crm )
      call prefetch(output%md_crm )
      call prefetch(output%du_crm )
      call prefetch(output%eu_crm )
      call prefetch(output%ed_crm )
      call prefetch(output%jt_crm )
      call prefetch(output%mx_crm )
      call prefetch(output%flux_qt       )
      call prefetch(output%fluxsgs_qt    )
      call prefetch(output%tkez          )
      call prefetch(output%tkew          )
      call prefetch(output%tkesgsz       )
      call prefetch(output%tkz           )
      call prefetch(output%flux_u        )
      call prefetch(output%flux_v        )
      call prefetch(output%flux_qp       )
      call prefetch(output%precflux      )
      call prefetch(output%qt_ls         )
      call prefetch(output%qt_trans      )
      call prefetch(output%qp_trans      )
      call prefetch(output%qp_fall       )
      call prefetch(output%qp_src        )
      call prefetch(output%qp_evp        )
      call prefetch(output%t_ls          )
      call prefetch(output%prectend      )
      call prefetch(output%precstend     )
      call prefetch(output%taux          )
      call prefetch(output%tauy          )
      call prefetch(output%z0m           )
      call prefetch(output%subcycle_factor )

      call prefetch(output%dt_sgs)
      call prefetch(output%dqv_sgs)
      call prefetch(output%dqc_sgs)
      call prefetch(output%dqi_sgs)
      call prefetch(output%dt_micro)
      call prefetch(output%dqv_micro)
      call prefetch(output%dqc_micro)
      call prefetch(output%dqi_micro)

      call prefetch(output%dt_dycor  )
      call prefetch(output%dqv_dycor )
      call prefetch(output%dqc_dycor )
      call prefetch(output%dqi_dycor )
      call prefetch(output%dt_sponge )
      call prefetch(output%dqv_sponge)
      call prefetch(output%dqc_sponge)
      call prefetch(output%dqi_sponge)

      call prefetch(output%rho_d_ls)
      call prefetch(output%rho_v_ls)
      call prefetch(output%rho_l_ls)
      call prefetch(output%rho_i_ls)

      ! Initialize 
      output%qcl = 0
      output%qci = 0
      output%qpl = 0
      output%qpi = 0

      output%tk = 0
      output%tkh = 0
      output%prec_crm = 0

      ! 2-moment process rates
      output%wvar = 0
      output%aut  = 0 
      output%acc  = 0
      output%evpc = 0
      output%evpr = 0
      output%mlt  = 0
      output%sub  = 0
      output%dep  = 0
      output%con  = 0

      output%cltot = 0
      output%cllow = 0
      output%clmed = 0
      output%clhgh = 0

      output%cldtop = 0
      output%precc = 0
      output%precl = 0
      output%precsc = 0
      output%precsl = 0

      output%qv_mean = 0
      output%qc_mean = 0
      output%qi_mean = 0
      output%qr_mean = 0
      output%qs_mean = 0
      output%qg_mean = 0
      output%qm_mean = 0
      output%bm_mean = 0
      output%rho_d_mean = 0
      output%rho_v_mean = 0

      output%nc_mean = 0
      output%ni_mean = 0
      output%nr_mean = 0

      output%ultend  = 0
      output%vltend  = 0
      output%sltend  = 0
      output%qltend  = 0
      output%qcltend = 0
      output%qiltend = 0

      output%t_vt_tend = 0
      output%q_vt_tend = 0
      output%u_vt_tend = 0
      output%t_vt_ls   = 0
      output%q_vt_ls   = 0
      output%u_vt_ls   = 0

      output%cld    = 0
      output%gicewp = 0
      output%gliqwp = 0

      output%liq_ice_exchange = 0
      output%vap_liq_exchange = 0
      output%vap_ice_exchange = 0

      output%mctot  = 0
      output%mcup   = 0
      output%mcdn   = 0
      output%mcuup  = 0
      output%mcudn  = 0

      ! Convective transport
      output%mu_crm = 0
      output%md_crm = 0
      output%eu_crm = 0
      output%du_crm = 0
      output%ed_crm = 0
      output%jt_crm = 0
      output%mx_crm = 0

      ! Other stuff...
      output%flux_qt       = 0
      output%fluxsgs_qt    = 0
      output%tkez          = 0
      output%tkew          = 0
      output%tkesgsz       = 0
      output%tkz           = 0
      output%flux_u        = 0
      output%flux_v        = 0
      output%flux_qp       = 0
      output%precflux      = 0
      output%qt_ls         = 0
      output%qt_trans      = 0
      output%qp_trans      = 0
      output%qp_fall       = 0
      output%qp_src        = 0
      output%qp_evp        = 0
      output%t_ls          = 0
      output%prectend      = 0
      output%precstend     = 0
      output%taux          = 0
      output%tauy          = 0
      output%z0m           = 0
      output%subcycle_factor = 0

      output%dt_sgs    = 0
      output%dqv_sgs   = 0
      output%dqc_sgs   = 0
      output%dqi_sgs   = 0
      output%dt_micro  = 0
      output%dqv_micro = 0
      output%dqc_micro = 0
      output%dqi_micro = 0

      output%dt_dycor   = 0
      output%dqv_dycor  = 0
      output%dqc_dycor  = 0
      output%dqi_dycor  = 0
      output%dt_sponge  = 0
      output%dqv_sponge = 0
      output%dqc_sponge = 0
      output%dqi_sponge = 0

      output%rho_d_ls = 0
      output%rho_v_ls = 0
      output%rho_l_ls = 0
      output%rho_i_ls = 0

   end subroutine crm_output_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_output_finalize(output, MMF_microphysics_scheme)
      type(crm_output_type), intent(inout) :: output
      character(len=*), intent(in) :: MMF_microphysics_scheme    ! CRM microphysics scheme

      if (allocated(output%qcl)) deallocate(output%qcl)
      if (allocated(output%qci)) deallocate(output%qci)
      if (allocated(output%qpl)) deallocate(output%qpl)
      if (allocated(output%qpi)) deallocate(output%qpi)
      if (allocated(output%tk )) deallocate(output%tk )
      if (allocated(output%tkh)) deallocate(output%tkh)
      if (allocated(output%prec_crm)) deallocate(output%prec_crm)

      if (allocated(output%wvar)) deallocate(output%wvar)
      if (allocated(output%aut)) deallocate(output%aut)
      if (allocated(output%acc)) deallocate(output%acc)
      if (allocated(output%evpc)) deallocate(output%evpc)
      if (allocated(output%evpr)) deallocate(output%evpr)
      if (allocated(output%mlt)) deallocate(output%mlt)
      if (allocated(output%sub)) deallocate(output%sub)
      if (allocated(output%dep)) deallocate(output%dep)
      if (allocated(output%con)) deallocate(output%con)

      if (allocated(output%cltot)) deallocate(output%cltot)
      if (allocated(output%cllow)) deallocate(output%cllow)
      if (allocated(output%clmed)) deallocate(output%clmed)
      if (allocated(output%clhgh)) deallocate(output%clhgh)
      if (allocated(output%cldtop)) deallocate(output%cldtop)
      if (allocated(output%precc)) deallocate(output%precc)
      if (allocated(output%precl)) deallocate(output%precl)
      if (allocated(output%precsc)) deallocate(output%precsc)
      if (allocated(output%precsl)) deallocate(output%precsl)

      if (allocated(output%qv_mean)) deallocate(output%qv_mean)
      if (allocated(output%qc_mean)) deallocate(output%qc_mean)
      if (allocated(output%qi_mean)) deallocate(output%qi_mean)
      if (allocated(output%qr_mean)) deallocate(output%qr_mean)
      if (allocated(output%qs_mean)) deallocate(output%qs_mean)
      if (allocated(output%qg_mean)) deallocate(output%qg_mean)
      if (allocated(output%qm_mean)) deallocate(output%qm_mean)
      if (allocated(output%bm_mean)) deallocate(output%bm_mean)
      if (allocated(output%rho_d_mean)) deallocate(output%rho_d_mean)
      if (allocated(output%rho_v_mean)) deallocate(output%rho_v_mean)
      
      if (allocated(output%nc_mean)) deallocate(output%nc_mean)
      if (allocated(output%ni_mean)) deallocate(output%ni_mean)
      if (allocated(output%nr_mean)) deallocate(output%nr_mean)

      if (allocated(output%ultend))  deallocate(output%ultend)
      if (allocated(output%vltend))  deallocate(output%vltend)
      if (allocated(output%sltend))  deallocate(output%sltend)
      if (allocated(output%qltend))  deallocate(output%qltend)
      if (allocated(output%qcltend)) deallocate(output%qcltend)
      if (allocated(output%qiltend)) deallocate(output%qiltend)

      if (allocated(output%t_vt_tend)) deallocate(output%t_vt_tend)
      if (allocated(output%q_vt_tend)) deallocate(output%q_vt_tend)
      if (allocated(output%u_vt_tend)) deallocate(output%u_vt_tend)
      if (allocated(output%t_vt_ls))   deallocate(output%t_vt_ls)
      if (allocated(output%q_vt_ls))   deallocate(output%q_vt_ls)
      if (allocated(output%u_vt_ls))   deallocate(output%u_vt_ls)

      if (allocated(output%cld)) deallocate(output%cld)
      if (allocated(output%gicewp)) deallocate(output%gicewp)
      if (allocated(output%gliqwp)) deallocate(output%gliqwp)

      if (allocated(output%liq_ice_exchange)) deallocate(output%liq_ice_exchange)
      if (allocated(output%vap_liq_exchange)) deallocate(output%vap_liq_exchange)
      if (allocated(output%vap_ice_exchange)) deallocate(output%vap_ice_exchange)

      if (allocated(output%mctot)) deallocate(output%mctot)
      if (allocated(output%mcup)) deallocate(output%mcup)
      if (allocated(output%mcdn)) deallocate(output%mcdn)
      if (allocated(output%mcuup)) deallocate(output%mcuup)
      if (allocated(output%mcudn)) deallocate(output%mcudn)

      if (allocated(output%mu_crm)) deallocate(output%mu_crm)
      if (allocated(output%md_crm)) deallocate(output%md_crm)
      if (allocated(output%du_crm)) deallocate(output%du_crm)
      if (allocated(output%eu_crm)) deallocate(output%eu_crm)
      if (allocated(output%ed_crm)) deallocate(output%ed_crm)
      if (allocated(output%jt_crm)) deallocate(output%jt_crm)
      if (allocated(output%mx_crm)) deallocate(output%mx_crm)

      if (allocated(output%flux_qt)) deallocate(output%flux_qt)
      if (allocated(output%fluxsgs_qt)) deallocate(output%fluxsgs_qt)
      if (allocated(output%tkez)) deallocate(output%tkez)
      if (allocated(output%tkew)) deallocate(output%tkew)
      if (allocated(output%tkesgsz)) deallocate(output%tkesgsz)
      if (allocated(output%tkz)) deallocate(output%tkz)
      if (allocated(output%flux_u)) deallocate(output%flux_u)
      if (allocated(output%flux_v)) deallocate(output%flux_v)
      if (allocated(output%flux_qp)) deallocate(output%flux_qp)
      if (allocated(output%precflux)) deallocate(output%precflux)
      if (allocated(output%qt_ls)) deallocate(output%qt_ls)
      if (allocated(output%qt_trans)) deallocate(output%qt_trans)
      if (allocated(output%qp_trans)) deallocate(output%qp_trans)
      if (allocated(output%qp_fall)) deallocate(output%qp_fall)
      if (allocated(output%qp_src)) deallocate(output%qp_src)
      if (allocated(output%qp_evp)) deallocate(output%qp_evp)
      if (allocated(output%t_ls)) deallocate(output%t_ls)
      if (allocated(output%prectend)) deallocate(output%prectend)
      if (allocated(output%precstend)) deallocate(output%precstend)
      if (allocated(output%taux)) deallocate(output%taux)
      if (allocated(output%tauy)) deallocate(output%tauy)
      if (allocated(output%z0m)) deallocate(output%z0m)
      if (allocated(output%subcycle_factor)) deallocate(output%subcycle_factor)

      if (allocated(output%dt_sgs   )) deallocate(output%dt_sgs)
      if (allocated(output%dqv_sgs  )) deallocate(output%dqv_sgs)
      if (allocated(output%dqc_sgs  )) deallocate(output%dqc_sgs)
      if (allocated(output%dqi_sgs  )) deallocate(output%dqi_sgs)
      if (allocated(output%dt_micro )) deallocate(output%dt_micro)
      if (allocated(output%dqv_micro)) deallocate(output%dqv_micro)
      if (allocated(output%dqc_micro)) deallocate(output%dqc_micro)
      if (allocated(output%dqi_micro)) deallocate(output%dqi_micro)

      if (allocated(output%dt_dycor  )) deallocate(output%dt_dycor  )
      if (allocated(output%dqv_dycor )) deallocate(output%dqv_dycor )
      if (allocated(output%dqc_dycor )) deallocate(output%dqc_dycor )
      if (allocated(output%dqi_dycor )) deallocate(output%dqi_dycor )
      if (allocated(output%dt_sponge )) deallocate(output%dt_sponge )
      if (allocated(output%dqv_sponge)) deallocate(output%dqv_sponge)
      if (allocated(output%dqc_sponge)) deallocate(output%dqc_sponge)
      if (allocated(output%dqi_sponge)) deallocate(output%dqi_sponge)

      if (allocated(output%rho_d_ls)) deallocate(output%rho_d_ls)
      if (allocated(output%rho_v_ls)) deallocate(output%rho_v_ls)
      if (allocated(output%rho_l_ls)) deallocate(output%rho_l_ls)
      if (allocated(output%rho_i_ls)) deallocate(output%rho_i_ls)

   end subroutine crm_output_finalize
   !------------------------------------------------------------------------------------------------

end module crm_output_module
