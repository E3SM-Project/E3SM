module crm_output_module
   use params,       only: crm_rknd
   use crmdims,      only: crm_nx, crm_ny, crm_nz
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
      real(crm_rknd), allocatable :: qc_mean(:,:)  ! mean cloud water
      real(crm_rknd), allocatable :: qi_mean(:,:)  ! mean cloud ice
      real(crm_rknd), allocatable :: qs_mean(:,:)  ! mean snow
      real(crm_rknd), allocatable :: qg_mean(:,:)  ! mean graupel
      real(crm_rknd), allocatable :: qr_mean(:,:)  ! mean rain
#ifdef m2005
      real(crm_rknd), allocatable :: nc_mean(:,:)  ! mean cloud water  (#/kg)
      real(crm_rknd), allocatable :: ni_mean(:,:)  ! mean cloud ice    (#/kg)
      real(crm_rknd), allocatable :: ns_mean(:,:)  ! mean snow         (#/kg)
      real(crm_rknd), allocatable :: ng_mean(:,:)  ! mean graupel      (#/kg)
      real(crm_rknd), allocatable :: nr_mean(:,:)  ! mean rain         (#/kg)

      ! Time and domain averaged process rates
      real(crm_rknd), allocatable :: aut_a (:,:)  ! cloud water autoconversion (1/s)
      real(crm_rknd), allocatable :: acc_a (:,:)  ! cloud water accretion (1/s)
      real(crm_rknd), allocatable :: evpc_a(:,:)  ! cloud water evaporation (1/s)
      real(crm_rknd), allocatable :: evpr_a(:,:)  ! rain evaporation (1/s)
      real(crm_rknd), allocatable :: mlt_a (:,:)  ! ice, snow, graupel melting (1/s)
      real(crm_rknd), allocatable :: sub_a (:,:)  ! ice, snow, graupel sublimation (1/s)
      real(crm_rknd), allocatable :: dep_a (:,:)  ! ice, snow, graupel deposition (1/s)
      real(crm_rknd), allocatable :: con_a (:,:)  ! cloud water condensation(1/s)
#endif /* m2005 */

#if defined( MMF_MOMENTUM_FEEDBACK )
      real(crm_rknd), allocatable :: ultend(:,:)            ! tendency of ul
      real(crm_rknd), allocatable :: vltend(:,:)            ! tendency of vl
#endif

#if defined( MMF_ESMT )
      real(crm_rknd), allocatable :: u_tend_esmt(:,:)       ! CRM scalar u-momentum tendency
      real(crm_rknd), allocatable :: v_tend_esmt(:,:)       ! CRM scalar v-momentum tendency
#endif

      real(crm_rknd), allocatable :: sltend  (:,:)          ! CRM output tendency of static energy
      real(crm_rknd), allocatable :: qltend  (:,:)          ! CRM output tendency of water vapor
      real(crm_rknd), allocatable :: qcltend (:,:)          ! CRM output tendency of cloud liquid water
      real(crm_rknd), allocatable :: qiltend (:,:)          ! CRM output tendency of cloud ice

      ! These are all time and spatial averages, on the GCM grid
      real(crm_rknd), allocatable :: cld   (:,:)      ! cloud fraction
      real(crm_rknd), allocatable :: gicewp(:,:)      ! ice water path
      real(crm_rknd), allocatable :: gliqwp(:,:)      ! ice water path
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
      real(crm_rknd), allocatable :: taux     (:)    ! zonal CRM surface stress perturbation      [N/m2]
      real(crm_rknd), allocatable :: tauy     (:)    ! merid CRM surface stress perturbation      [N/m2]
      real(crm_rknd), allocatable :: z0m          (:)    ! surface stress                             [N/m2]
      real(crm_rknd), allocatable :: timing_factor(:)    ! crm cpu efficiency

#ifdef MAML
      ! MAML variables
      real(crm_rknd), allocatable :: crm_pcp(ncrms,crm_nx,crm_ny) ! CRM precip rate for MAML (m/s)
      real(crm_rknd), allocatable :: crm_snw(ncrms,crm_nx,crm_ny) ! CRM snow rate for MAML (m/s)
#endif

   end type crm_output_type

contains

   !------------------------------------------------------------------------------------------------
   subroutine crm_output_initialize(output, ncol, nlev)
      type(crm_output_type), intent(inout) :: output
      integer, intent(in), optional :: ncol, nlev

      ! Allocate arrays if dimensions are passed as input
      if (present(ncol)) then

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

         if (.not. allocated(output%qc_mean)) allocate(output%qc_mean(ncol,nlev))
         if (.not. allocated(output%qi_mean)) allocate(output%qi_mean(ncol,nlev))
         if (.not. allocated(output%qs_mean)) allocate(output%qs_mean(ncol,nlev))
         if (.not. allocated(output%qg_mean)) allocate(output%qg_mean(ncol,nlev))
         if (.not. allocated(output%qr_mean)) allocate(output%qr_mean(ncol,nlev))

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
         call prefetch(output%qc_mean)
         call prefetch(output%qi_mean)
         call prefetch(output%qs_mean)
         call prefetch(output%qg_mean)
         call prefetch(output%qr_mean)

#ifdef m2005
         if (.not. allocated(output%nc_mean)) allocate(output%nc_mean(ncol,nlev))
         if (.not. allocated(output%ni_mean)) allocate(output%ni_mean(ncol,nlev))
         if (.not. allocated(output%ns_mean)) allocate(output%ns_mean(ncol,nlev))
         if (.not. allocated(output%ng_mean)) allocate(output%ng_mean(ncol,nlev))
         if (.not. allocated(output%nr_mean)) allocate(output%nr_mean(ncol,nlev))

         if (.not. allocated(output%aut_a )) allocate(output%aut_a (ncol,nlev))
         if (.not. allocated(output%acc_a )) allocate(output%acc_a (ncol,nlev))
         if (.not. allocated(output%evpc_a)) allocate(output%evpc_a(ncol,nlev))
         if (.not. allocated(output%evpr_a)) allocate(output%evpr_a(ncol,nlev))
         if (.not. allocated(output%mlt_a )) allocate(output%mlt_a (ncol,nlev))
         if (.not. allocated(output%sub_a )) allocate(output%sub_a (ncol,nlev))
         if (.not. allocated(output%dep_a )) allocate(output%dep_a (ncol,nlev))
         if (.not. allocated(output%con_a )) allocate(output%con_a (ncol,nlev))
#endif /* m2005 */

#if defined( MMF_MOMENTUM_FEEDBACK )
         if (.not. allocated(output%ultend )) allocate(output%ultend (ncol,nlev))
         if (.not. allocated(output%vltend )) allocate(output%vltend (ncol,nlev))
#endif

#if defined( MMF_ESMT )
         if (.not. allocated(output%u_tend_esmt )) allocate(output%u_tend_esmt (ncol,nlev))
         if (.not. allocated(output%v_tend_esmt )) allocate(output%v_tend_esmt (ncol,nlev))
#endif

#ifdef MAML
         if (.not. allocated(output%crm_pcp)) allocate(output%crm_pcp(ncol,crm_nx,crm_ny))
         if (.not. allocated(output%crm_snw)) allocate(output%crm_snw(ncol,crm_nx,crm_ny))
         call prefetch(output%crm_pcp)
         call prefetch(output%crm_snw)
#endif
         
         if (.not. allocated(output%sltend ))  allocate(output%sltend (ncol,nlev))
         if (.not. allocated(output%qltend ))  allocate(output%qltend (ncol,nlev))
         if (.not. allocated(output%qcltend))  allocate(output%qcltend(ncol,nlev))
         if (.not. allocated(output%qiltend))  allocate(output%qiltend(ncol,nlev))

         if (.not. allocated(output%cld   )) allocate(output%cld   (ncol,nlev))  ! cloud fraction
         if (.not. allocated(output%gicewp)) allocate(output%gicewp(ncol,nlev))  ! ice water path
         if (.not. allocated(output%gliqwp)) allocate(output%gliqwp(ncol,nlev))  ! ice water path
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
         if (.not. allocated(output%timing_factor)) allocate(output%timing_factor(ncol))

         call prefetch(output%sltend  )
         call prefetch(output%qltend  )
         call prefetch(output%qcltend )
         call prefetch(output%qiltend )
         call prefetch(output%cld    )
         call prefetch(output%gicewp )
         call prefetch(output%gliqwp )
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
         call prefetch(output%timing_factor )

      end if ! present(ncol)

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

      output%qc_mean = 0
      output%qi_mean = 0
      output%qs_mean = 0
      output%qg_mean = 0
      output%qr_mean = 0
#ifdef m2005
      output%nc_mean = 0
      output%ni_mean = 0
      output%ns_mean = 0
      output%ng_mean = 0
      output%nr_mean = 0

      output%aut_a = 0
      output%acc_a = 0
      output%evpc_a = 0
      output%evpr_a = 0
      output%mlt_a = 0
      output%sub_a = 0
      output%dep_a = 0
      output%con_a = 0
#endif

#if defined( MMF_MOMENTUM_FEEDBACK )
      output%ultend = 0
      output%vltend = 0
#endif

#if defined( MMF_ESMT )
      output%u_tend_esmt = 0
      output%v_tend_esmt = 0
#endif

#ifdef MAML
      output%crm_pcp = 0
      output%crm_snw = 0
#endif

      output%sltend  = 0
      output%qltend  = 0
      output%qcltend = 0
      output%qiltend = 0

      output%cld    = 0
      output%gicewp = 0
      output%gliqwp = 0
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
      output%taux      = 0
      output%tauy      = 0
      output%z0m           = 0
      output%timing_factor = 0

   end subroutine crm_output_initialize
   !------------------------------------------------------------------------------------------------
   subroutine crm_output_finalize(output)
      type(crm_output_type), intent(inout) :: output
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

      if (allocated(output%qc_mean)) deallocate(output%qc_mean)
      if (allocated(output%qi_mean)) deallocate(output%qi_mean)
      if (allocated(output%qs_mean)) deallocate(output%qs_mean)
      if (allocated(output%qg_mean)) deallocate(output%qg_mean)
      if (allocated(output%qr_mean)) deallocate(output%qr_mean)
#ifdef m2005
      if (allocated(output%nc_mean)) deallocate(output%nc_mean)
      if (allocated(output%ni_mean)) deallocate(output%ni_mean)
      if (allocated(output%ns_mean)) deallocate(output%ns_mean)
      if (allocated(output%ng_mean)) deallocate(output%ng_mean)
      if (allocated(output%nr_mean)) deallocate(output%nr_mean)

      ! Time and domain-averaged process rates
      if (allocated(output%aut_a)) deallocate(output%aut_a)
      if (allocated(output%acc_a)) deallocate(output%acc_a)
      if (allocated(output%evpc_a)) deallocate(output%evpc_a)
      if (allocated(output%evpr_a)) deallocate(output%evpr_a)
      if (allocated(output%mlt_a)) deallocate(output%mlt_a)
      if (allocated(output%sub_a)) deallocate(output%sub_a)
      if (allocated(output%dep_a)) deallocate(output%dep_a)
      if (allocated(output%con_a)) deallocate(output%con_a)
#endif

#if defined( MMF_MOMENTUM_FEEDBACK )
      if (allocated(output%ultend)) deallocate(output%ultend)
      if (allocated(output%vltend)) deallocate(output%vltend)
#endif

#if defined( MMF_ESMT )
      if (allocated(output%u_tend_esmt)) deallocate(output%u_tend_esmt)
      if (allocated(output%v_tend_esmt)) deallocate(output%v_tend_esmt)
#endif

#ifdef MAML
      if (allocated(output%crm_pcp)) deallocate(output%crm_pcp)
      if (allocated(output%crm_snw)) deallocate(output%crm_snw)
#endif

      if (allocated(output%sltend)) deallocate(output%sltend)
      if (allocated(output%qltend)) deallocate(output%qltend)
      if (allocated(output%qcltend)) deallocate(output%qcltend)
      if (allocated(output%qiltend)) deallocate(output%qiltend)

      if (allocated(output%cld)) deallocate(output%cld)
      if (allocated(output%gicewp)) deallocate(output%gicewp)
      if (allocated(output%gliqwp)) deallocate(output%gliqwp)
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
      if (allocated(output%timing_factor)) deallocate(output%timing_factor)

   end subroutine crm_output_finalize
   !------------------------------------------------------------------------------------------------

end module crm_output_module
