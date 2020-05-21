module crm_output_module
   use params,       only: crm_rknd, crm_iknd
   use crmdims,      only: crm_nx, crm_ny, crm_nz
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
   end type crm_output_type

contains

   !------------------------------------------------------------------------------------------------
   subroutine crm_output_initialize(output, ncol, nlev)
      type(crm_output_type), intent(inout) :: output
      integer(crm_iknd), intent(in), optional :: ncol, nlev

      ! Allocate arrays if dimensions are passed as input
      if (present(ncol)) then

         ! Allocate instantaneous outputs
         allocate(output%qcl     (ncol,crm_nx,crm_ny,crm_nz))
         allocate(output%qci     (ncol,crm_nx,crm_ny,crm_nz))
         allocate(output%qpl     (ncol,crm_nx,crm_ny,crm_nz))
         allocate(output%qpi     (ncol,crm_nx,crm_ny,crm_nz))

         allocate(output%tk      (ncol,crm_nx,crm_ny,crm_nz))
         allocate(output%tkh     (ncol,crm_nx,crm_ny,crm_nz))
         allocate(output%prec_crm(ncol,crm_nx,crm_ny       ))

         allocate(output%wvar    (ncol,crm_nx,crm_ny,crm_nz))
         allocate(output%aut     (ncol,crm_nx,crm_ny,crm_nz))
         allocate(output%acc     (ncol,crm_nx,crm_ny,crm_nz))
         allocate(output%evpc    (ncol,crm_nx,crm_ny,crm_nz))
         allocate(output%evpr    (ncol,crm_nx,crm_ny,crm_nz))
         allocate(output%mlt     (ncol,crm_nx,crm_ny,crm_nz))
         allocate(output%sub     (ncol,crm_nx,crm_ny,crm_nz))
         allocate(output%dep     (ncol,crm_nx,crm_ny,crm_nz))
         allocate(output%con     (ncol,crm_nx,crm_ny,crm_nz))


         ! Allocate domain and time-averaged fields
         allocate(output%cltot (ncol))
         allocate(output%cllow (ncol))
         allocate(output%clmed (ncol))
         allocate(output%clhgh (ncol))

         allocate(output%precc (ncol))
         allocate(output%precl (ncol))
         allocate(output%precsc(ncol))
         allocate(output%precsl(ncol))

         ! NOTE: this output had a bug in the previous implementation
         allocate(output%cldtop       (ncol,nlev))

         allocate(output%qc_mean      (ncol,nlev))
         allocate(output%qi_mean      (ncol,nlev))
         allocate(output%qs_mean      (ncol,nlev))
         allocate(output%qg_mean      (ncol,nlev))
         allocate(output%qr_mean      (ncol,nlev))

         allocate(output%sltend       (ncol,nlev))
         allocate(output%qltend       (ncol,nlev))
         allocate(output%qcltend      (ncol,nlev))
         allocate(output%qiltend      (ncol,nlev))

         allocate(output%cld          (ncol,nlev))  ! cloud fraction
         allocate(output%gicewp       (ncol,nlev))  ! ice water path
         allocate(output%gliqwp       (ncol,nlev))  ! ice water path
         allocate(output%mctot        (ncol,nlev))  ! cloud mass flux
         allocate(output%mcup         (ncol,nlev))  ! updraft cloud mass flux
         allocate(output%mcdn         (ncol,nlev))  ! downdraft cloud mass flux
         allocate(output%mcuup        (ncol,nlev))  ! unsat updraft cloud mass flux
         allocate(output%mcudn        (ncol,nlev))  ! unsat downdraft cloud mass flux

         allocate(output%mu_crm       (ncol,nlev))  ! mass flux up
         allocate(output%md_crm       (ncol,nlev))  ! mass flux down
         allocate(output%du_crm       (ncol,nlev))  ! mass detrainment from updraft
         allocate(output%eu_crm       (ncol,nlev))  ! mass entrainment from updraft
         allocate(output%ed_crm       (ncol,nlev))  ! mass detrainment from downdraft
         allocate(output%jt_crm       (ncol     ))       ! index of cloud (convection) top
         allocate(output%mx_crm       (ncol     ))       ! index of cloud (convection) bottom

         allocate(output%flux_qt      (ncol,nlev))
         allocate(output%fluxsgs_qt   (ncol,nlev))
         allocate(output%tkez         (ncol,nlev))
         allocate(output%tkesgsz      (ncol,nlev))
         allocate(output%tkz          (ncol,nlev))
         allocate(output%flux_u       (ncol,nlev))
         allocate(output%flux_v       (ncol,nlev))
         allocate(output%flux_qp      (ncol,nlev))
         allocate(output%precflux     (ncol,nlev))
         allocate(output%qt_ls        (ncol,nlev))
         allocate(output%qt_trans     (ncol,nlev))
         allocate(output%qp_trans     (ncol,nlev))
         allocate(output%qp_fall      (ncol,nlev))
         allocate(output%qp_src       (ncol,nlev))
         allocate(output%qp_evp       (ncol,nlev))
         allocate(output%t_ls         (ncol,nlev))
         allocate(output%prectend     (ncol     ))
         allocate(output%precstend    (ncol     ))
         allocate(output%taux         (ncol     ))
         allocate(output%tauy         (ncol     ))
         allocate(output%z0m          (ncol     ))
         allocate(output%timing_factor(ncol     ))

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
      deallocate(output%qcl)
      deallocate(output%qci)
      deallocate(output%qpl)
      deallocate(output%qpi)
      deallocate(output%tk )
      deallocate(output%tkh)
      deallocate(output%prec_crm)

      deallocate(output%wvar)
      deallocate(output%aut)
      deallocate(output%acc)
      deallocate(output%evpc)
      deallocate(output%evpr)
      deallocate(output%mlt)
      deallocate(output%sub)
      deallocate(output%dep)
      deallocate(output%con)

      deallocate(output%cltot)
      deallocate(output%cllow)
      deallocate(output%clmed)
      deallocate(output%clhgh)
      deallocate(output%cldtop)
      deallocate(output%precc)
      deallocate(output%precl)
      deallocate(output%precsc)
      deallocate(output%precsl)

      deallocate(output%qc_mean)
      deallocate(output%qi_mean)
      deallocate(output%qs_mean)
      deallocate(output%qg_mean)
      deallocate(output%qr_mean)
      deallocate(output%sltend)
      deallocate(output%qltend)
      deallocate(output%qcltend)
      deallocate(output%qiltend)

      deallocate(output%cld)
      deallocate(output%gicewp)
      deallocate(output%gliqwp)
      deallocate(output%mctot)
      deallocate(output%mcup)
      deallocate(output%mcdn)
      deallocate(output%mcuup)
      deallocate(output%mcudn)

      deallocate(output%mu_crm)
      deallocate(output%md_crm)
      deallocate(output%du_crm)
      deallocate(output%eu_crm)
      deallocate(output%ed_crm)
      deallocate(output%jt_crm)
      deallocate(output%mx_crm)

      deallocate(output%flux_qt)
      deallocate(output%fluxsgs_qt)
      deallocate(output%tkez)
      deallocate(output%tkesgsz)
      deallocate(output%tkz)
      deallocate(output%flux_u)
      deallocate(output%flux_v)
      deallocate(output%flux_qp)
      deallocate(output%precflux)
      deallocate(output%qt_ls)
      deallocate(output%qt_trans)
      deallocate(output%qp_trans)
      deallocate(output%qp_fall)
      deallocate(output%qp_src)
      deallocate(output%qp_evp)
      deallocate(output%t_ls)
      deallocate(output%prectend)
      deallocate(output%precstend)
      deallocate(output%taux)
      deallocate(output%tauy)
      deallocate(output%z0m)
      deallocate(output%timing_factor)

   end subroutine crm_output_finalize
   !------------------------------------------------------------------------------------------------

end module crm_output_module
