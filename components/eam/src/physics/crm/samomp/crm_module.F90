!---------------------------------------------------------------------------------------------------
! MMF main driver (original version by Marat Khairoutdinov)
!---------------------------------------------------------------------------------------------------
module crm_module
  use perf_mod
  use task_init_mod, only: task_init
  use params_kind, only: crm_rknd, r8
  use abcoefs_mod, only: abcoefs
  use kurant_mod, only: kurant
  use setperturb_mod, only: setperturb
  use boundaries_mod, only: boundaries
  use forcing_mod, only: forcing
  use advect_mom_mod, only: advect_mom
  use adams_mod, only: adams
  use advect_all_scalars_mod, only: advect_all_scalars
  use sat_mod
  use crm_surface_mod
#ifdef sam1mom
  use precip_init_mod
#endif
  use zero_mod
  use buoyancy_mod
  use pressure_mod
  use uvw_mod
  use diagnose_mod
  use damping_mod
  use ice_fall_mod
  use coriolis_mod
  use setparm_mod,            only: setparm
  use crm_state_module,       only: crm_state_type
  use crm_rad_module,         only: crm_rad_type
  use crm_input_module,       only: crm_input_type
  use crm_output_module,      only: crm_output_type
  use crm_ecpp_output_module, only: crm_ecpp_output_type
#ifdef ECPP  
  use module_ecpp_crm_driver, only: ecpp_crm_stat, ecpp_crm_init, ecpp_crm_cleanup
#endif

contains

subroutine crm( ncrms, dt_gl, plev,       &
                crm_input, crm_state, crm_rad, &
                crm_ecpp_output, crm_output, crm_clear_rh, &
                latitude0, longitude0, gcolp, igstep, &
                use_VT, VT_wn_max, &
                use_crm_accel_in, crm_accel_factor_in, crm_accel_uv_in)
  !-----------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------
  use shr_const_mod         , only: SHR_CONST_TKFRZ
  use vars
  use params
  use microphysics
  use sgs
  use crmtracers
  use scalar_momentum_mod
#ifdef MODAL_AERO
  use modal_aero_data       , only: ntot_amode
#endif
  use crmdims               , only: crm_nx_rad, crm_ny_rad
#ifdef ECPP
  use ecppvars              , only: qlsink, precr, precsolid, &
                                    area_bnd_final, area_bnd_sum, area_cen_final, area_cen_sum, &
                                    mass_bnd_final, mass_bnd_sum, rh_cen_sum, qcloud_cen_sum, qice_cen_sum, &
                                    qlsink_cen_sum, precr_cen_sum, precsolid_cen_sum, xkhvsum, wup_thresh, wdown_thresh, &
                                    wwqui_cen_sum, wwqui_bnd_sum, wwqui_cloudy_cen_sum, wwqui_cloudy_bnd_sum, &
                                    qlsink_bf_cen_sum, qlsink_avg_cen_sum, prain_cen_sum, qlsink_bf, prain
  use ecppvars              , only: NCLASS_CL, ncls_ecpp_in, NCLASS_PR
#endif /* ECPP */
  use accelerate_crm_mod    , only: use_crm_accel, crm_accel_factor, crm_accel_nstop, accelerate_crm, crm_accel_uv
#ifndef MMF_STANDALONE
  use cam_abortutils        , only: endrun
#endif

  implicit none

  !-----------------------------------------------------------------------------------------------
  ! Interface variable declarations
  !-----------------------------------------------------------------------------------------------
  integer , intent(in   ) :: ncrms                            ! Number of "vector" GCM columns to push down into CRM for SIMD vectorization / more threading
  integer , intent(in   ) :: plev                             ! number of levels in parent model
  real(r8), intent(in   ) :: dt_gl                            ! global model's time step
  type(crm_input_type), target,intent(in   ) :: crm_input
  type(crm_state_type), target,intent(inout) :: crm_state
  type(crm_rad_type),   target,intent(inout) :: crm_rad
  type(crm_ecpp_output_type), intent(inout) :: crm_ecpp_output
  type(crm_output_type), target, intent(inout) :: crm_output
  real(r8), dimension(ncrms,nzm), intent(  out) :: crm_clear_rh
  real(crm_rknd), intent(in) :: latitude0(:)
  real(crm_rknd), intent(in) :: longitude0(:)
  integer       , intent(in) :: igstep
  integer       , intent(in) :: gcolp(:)
  logical       , intent(in) :: use_VT
  integer       , intent(in) :: VT_wn_max
  logical       , intent(in) :: use_crm_accel_in
  real(crm_rknd), intent(in) :: crm_accel_factor_in
  logical       , intent(in) :: crm_accel_uv_in

  !-----------------------------------------------------------------------------------------------
  ! Local variable declarations
  !-----------------------------------------------------------------------------------------------

  real(r8),       parameter :: umax = 0.5D0*crm_dx/crm_dt       ! maxumum ampitude of the l.s. wind
  real(r8),       parameter :: wmin = 2.                      ! minimum up/downdraft velocity for stat
  real(crm_rknd), parameter :: cwp_threshold = 0.001D0          ! threshold for cloud condensate for shaded fraction calculation
  integer,        parameter :: perturb_seed_scale = 1000      ! scaling value for setperturb() seed value (seed = gcol * perturb_seed_scale)
  real(r8)        :: crm_run_time                             ! length of CRM integration
  real(r8)        :: icrm_run_time                            ! = 1 / crm_run_time
  real(r8)        :: factor_xy, factor_xyt, idt_gl
  real(crm_rknd)  :: tmp1, tmp2, tmp
  real(crm_rknd)  :: u2z,v2z,w2z
  integer         :: i,j,k,l,ptop,nn,icyc,icrm
  integer         :: kx
  real(crm_rknd)  :: qsat, omg, rh_tmp
  real(crm_rknd), allocatable  :: colprec(:), colprecs(:)
  real(crm_rknd), allocatable  :: ustar(:), bflx(:), wnd(:)
  real(r8)      , allocatable  :: qtot (:,:)    ! Total water for water conservation check

  ! These should all be inputs
  integer         :: iseed             ! seed for random perturbation
  ! variables for radiation grouping method
  real(crm_rknd) :: crm_nx_rad_fac
  real(crm_rknd) :: crm_ny_rad_fac
  integer        :: i_rad
  integer        :: j_rad
  logical :: crm_accel_ceaseflag   ! indicates if accelerate_crm needs to be aborted for remainder of crm call

  ! Arrays
  real(crm_rknd), allocatable :: t00(:,:)
  real(crm_rknd), allocatable :: tln  (:,:)
  real(crm_rknd), allocatable :: qln  (:,:)
  real(crm_rknd), allocatable :: qccln(:,:)
  real(crm_rknd), allocatable :: qiiln(:,:)
  real(crm_rknd), allocatable :: uln  (:,:)
  real(crm_rknd), allocatable :: vln  (:,:)
#if defined(MMF_ESMT)
  real(crm_rknd), allocatable  :: uln_esmt(:,:)
  real(crm_rknd), allocatable  :: vln_esmt(:,:)     ! tempoerary variables for expliciit scalar momentum transport
#endif
  real(crm_rknd), allocatable  :: cwp     (:,:,:)
  real(crm_rknd), allocatable  :: cwph    (:,:,:)
  real(crm_rknd), allocatable  :: cwpm    (:,:,:)
  real(crm_rknd), allocatable  :: cwpl    (:,:,:)
  logical       , allocatable  :: flag_top(:,:,:)
  real(crm_rknd), allocatable  :: cltemp  (:,:,:)
  real(crm_rknd), allocatable  :: cmtemp  (:,:,:)
  real(crm_rknd), allocatable  :: chtemp  (:,:,:)
  real(crm_rknd), allocatable  :: cttemp  (:,:,:)

  real(r8), allocatable :: dd_crm (:,:)     ! mass entraiment from downdraft
  real(r8), allocatable :: mui_crm(:,:)     ! mass flux up at the interface
  real(r8), allocatable :: mdi_crm(:,:)     ! mass flux down at the interface

  real(r8), allocatable :: crm_clear_rh_cnt(:,:) ! counter for clear air relative humidity

  real(crm_rknd), pointer :: crm_input_zmid(:,:)
  real(crm_rknd), pointer :: crm_input_zint(:,:)
  real(crm_rknd), pointer :: crm_input_tl(:,:)
  real(crm_rknd), pointer :: crm_input_ql(:,:)
  real(crm_rknd), pointer :: crm_input_qccl(:,:)
  real(crm_rknd), pointer :: crm_input_qiil(:,:)
  real(crm_rknd), pointer :: crm_input_ps(:)
  real(crm_rknd), pointer :: crm_input_pmid(:,:)
  real(crm_rknd), pointer :: crm_input_pint(:,:)
  real(crm_rknd), pointer :: crm_input_pdel(:,:)
  real(crm_rknd), pointer :: crm_input_phis(:)
  real(crm_rknd), pointer :: crm_input_ul(:,:)
  real(crm_rknd), pointer :: crm_input_vl(:,:)
  real(crm_rknd), pointer :: crm_input_ocnfrac(:)
  real(crm_rknd), pointer :: crm_input_tau00(:)
  real(crm_rknd), pointer :: crm_input_wndls(:)
  real(crm_rknd), pointer :: crm_input_bflxls(:)
  real(crm_rknd), pointer :: crm_input_fluxu00(:)
  real(crm_rknd), pointer :: crm_input_fluxv00(:)
  real(crm_rknd), pointer :: crm_input_fluxt00(:)
  real(crm_rknd), pointer :: crm_input_fluxq00(:)

#if defined( m2005 ) && defined( MODAL_AERO )
  real(crm_rknd), pointer :: crm_input_naermod(:,:,:) 
  real(crm_rknd), pointer :: crm_input_vaerosol(:,:,:)
  real(crm_rknd), pointer :: crm_input_hygro(:,:,:)
#endif

#if defined( SP_ESMT )
  real(crm_rknd), pointer :: crm_input_ul_esmt(:,:) 
  real(crm_rknd), pointer :: crm_input_vl_esmt(:,:) 
#endif

  real(crm_rknd), pointer :: crm_output_qcl(:,:,:,:)
  real(crm_rknd), pointer :: crm_output_qci(:,:,:,:)
  real(crm_rknd), pointer :: crm_output_qpl(:,:,:,:)
  real(crm_rknd), pointer :: crm_output_qpi(:,:,:,:)

  real(crm_rknd), pointer :: crm_output_tk (:,:,:,:)
  real(crm_rknd), pointer :: crm_output_tkh(:,:,:,:)
  real(crm_rknd), pointer :: crm_output_prec_crm(:,:,:)

  ! 2-moment process rates
  real(crm_rknd), pointer :: crm_output_wvar(:,:,:,:)
  real(crm_rknd), pointer :: crm_output_aut (:,:,:,:)
  real(crm_rknd), pointer :: crm_output_acc (:,:,:,:)
  real(crm_rknd), pointer :: crm_output_evpc(:,:,:,:)
  real(crm_rknd), pointer :: crm_output_evpr(:,:,:,:)
  real(crm_rknd), pointer :: crm_output_mlt (:,:,:,:)
  real(crm_rknd), pointer :: crm_output_sub (:,:,:,:)
  real(crm_rknd), pointer :: crm_output_dep (:,:,:,:)
  real(crm_rknd), pointer :: crm_output_con (:,:,:,:)

  ! Cloud area fractions
  real(crm_rknd), pointer :: crm_output_cltot(:)
  real(crm_rknd), pointer :: crm_output_clhgh(:)
  real(crm_rknd), pointer :: crm_output_clmed(:)
  real(crm_rknd), pointer :: crm_output_cllow(:)

  real(crm_rknd), pointer :: crm_output_cldtop(:,:)
  real(crm_rknd), pointer :: crm_output_precc(:)
  real(crm_rknd), pointer :: crm_output_precl(:)
  real(crm_rknd), pointer :: crm_output_precsc(:)
  real(crm_rknd), pointer :: crm_output_precsl(:)

  real(crm_rknd), pointer :: crm_output_qc_mean(:,:)
  real(crm_rknd), pointer :: crm_output_qi_mean(:,:)
  real(crm_rknd), pointer :: crm_output_qs_mean(:,:)
  real(crm_rknd), pointer :: crm_output_qg_mean(:,:)
  real(crm_rknd), pointer :: crm_output_qr_mean(:,:)
#ifdef m2005
  real(crm_rknd), pointer :: crm_output_nc_mean(:,:)
  real(crm_rknd), pointer :: crm_output_ni_mean(:,:)
  real(crm_rknd), pointer :: crm_output_ns_mean(:,:)
  real(crm_rknd), pointer :: crm_output_ng_mean(:,:)
  real(crm_rknd), pointer :: crm_output_nr_mean(:,:)

  ! Time and domain averaged process rates
  real(crm_rknd), pointer :: crm_output_aut_a (:,:)
  real(crm_rknd), pointer :: crm_output_acc_a (:,:)
  real(crm_rknd), pointer :: crm_output_evpc_a(:,:)
  real(crm_rknd), pointer :: crm_output_evpr_a(:,:)
  real(crm_rknd), pointer :: crm_output_mlt_a (:,:)
  real(crm_rknd), pointer :: crm_output_sub_a (:,:)
  real(crm_rknd), pointer :: crm_output_dep_a (:,:)
  real(crm_rknd), pointer :: crm_output_con_a (:,:)
#endif /* m2005 */

#if defined( SPMOMTRANS )
  real(crm_rknd), pointer :: crm_output_ultend(:,:)
  real(crm_rknd), pointer :: crm_output_vltend(:,:)
#endif

#if defined( SP_ESMT )
  real(crm_rknd), pointer :: crm_output_u_tend_esmt(:,:)
  real(crm_rknd), pointer :: crm_output_v_tend_esmt(:,:)
#endif
  real(crm_rknd), pointer :: crm_output_sltend  (:,:)
  real(crm_rknd), pointer :: crm_output_qltend  (:,:)
  real(crm_rknd), pointer :: crm_output_qcltend (:,:)
  real(crm_rknd), pointer :: crm_output_qiltend (:,:)

  ! These are all time and spatial averages, on the GCM grid
  real(crm_rknd), pointer :: crm_output_cld   (:,:)
  real(crm_rknd), pointer :: crm_output_gicewp(:,:)
  real(crm_rknd), pointer :: crm_output_gliqwp(:,:)
  real(crm_rknd), pointer :: crm_output_mctot (:,:)
  real(crm_rknd), pointer :: crm_output_mcup  (:,:)
  real(crm_rknd), pointer :: crm_output_mcdn  (:,:)
  real(crm_rknd), pointer :: crm_output_mcuup (:,:)
  real(crm_rknd), pointer :: crm_output_mcudn (:,:)

  ! For convective transport
  real(crm_rknd), pointer :: crm_output_mu_crm(:,:)
  real(crm_rknd), pointer :: crm_output_md_crm(:,:)
  real(crm_rknd), pointer :: crm_output_du_crm(:,:)
  real(crm_rknd), pointer :: crm_output_eu_crm(:,:)
  real(crm_rknd), pointer :: crm_output_ed_crm(:,:)
  real(crm_rknd), pointer :: crm_output_jt_crm(:)
  real(crm_rknd), pointer :: crm_output_mx_crm(:)

  ! Other stuff...
  real(crm_rknd), pointer :: crm_output_flux_qt(:,:)
  real(crm_rknd), pointer :: crm_output_fluxsgs_qt(:,:)
  real(crm_rknd), pointer :: crm_output_tkez(:,:)
  real(crm_rknd), pointer :: crm_output_tkew(:,:)
  real(crm_rknd), pointer :: crm_output_tkesgsz(:,:)
  real(crm_rknd), pointer :: crm_output_tkz(:,:)
  real(crm_rknd), pointer :: crm_output_flux_u(:,:)
  real(crm_rknd), pointer :: crm_output_flux_v(:,:)
  real(crm_rknd), pointer :: crm_output_flux_qp(:,:)
  real(crm_rknd), pointer :: crm_output_precflux(:,:)
  real(crm_rknd), pointer :: crm_output_qt_ls(:,:)
  real(crm_rknd), pointer :: crm_output_qt_trans(:,:)
  real(crm_rknd), pointer :: crm_output_qp_trans(:,:)
  real(crm_rknd), pointer :: crm_output_qp_fall(:,:)
  real(crm_rknd), pointer :: crm_output_qp_src(:,:)
  real(crm_rknd), pointer :: crm_output_qp_evp(:,:)
  real(crm_rknd), pointer :: crm_output_t_ls(:,:)
  real(crm_rknd), pointer :: crm_output_prectend (:)
  real(crm_rknd), pointer :: crm_output_precstend(:)
  real(crm_rknd), pointer :: crm_output_taux(:)
  real(crm_rknd), pointer :: crm_output_tauy(:)
  real(crm_rknd), pointer :: crm_output_z0m(:)
  real(crm_rknd), pointer :: crm_output_subcycle_factor(:) 

#ifdef MAML
  ! MAML variables
  real(crm_rknd), pointer :: crm_output_crm_pcp(:,:,:)
  real(crm_rknd), pointer :: crm_output_crm_snw(:,:,:)
#endif

  real(crm_rknd), pointer :: crm_rad_temperature  (:,:,:,:)
  real(crm_rknd), pointer :: crm_rad_qv           (:,:,:,:)
  real(crm_rknd), pointer :: crm_rad_qc           (:,:,:,:)
  real(crm_rknd), pointer :: crm_rad_qi           (:,:,:,:)
  real(crm_rknd), pointer :: crm_rad_cld          (:,:,:,:)
  real(crm_rknd), pointer :: crm_rad_qrad         (:,:,:,:)
  real(crm_rknd), pointer :: crm_rad_nc           (:,:,:,:)
  real(crm_rknd), pointer :: crm_rad_ni           (:,:,:,:)
  real(crm_rknd), pointer :: crm_rad_qs           (:,:,:,:)
  real(crm_rknd), pointer :: crm_rad_ns           (:,:,:,:)

  real(crm_rknd), pointer :: crm_state_u_wind     (:,:,:,:)
  real(crm_rknd), pointer :: crm_state_v_wind     (:,:,:,:)
  real(crm_rknd), pointer :: crm_state_w_wind     (:,:,:,:)
  real(crm_rknd), pointer :: crm_state_temperature(:,:,:,:)
  real(crm_rknd), pointer :: crm_state_qt         (:,:,:,:)
  real(crm_rknd), pointer :: crm_state_qp         (:,:,:,:)
  real(crm_rknd), pointer :: crm_state_qn         (:,:,:,:)
  !\--------------------------------------------------------------
  ! Execute code starts here...
  !/

  ! Variance transport not yet implemented in OpenMP, so abort if someone tries
  ! to run with it
  if (use_VT) then
#ifdef MMF_STANDALONE
     write(0,*) 'Variance transport not supported for OpenMP build.'
     stop
#else
     call endrun('Variance transport not supported for OpenMP build.')
#endif
  end if

  use_crm_accel    = use_crm_accel_in   
  crm_accel_factor = crm_accel_factor_in
  crm_accel_uv = crm_accel_uv_in

  call allocate_crm()
  call allocate_params(ncrms)
  call allocate_vars(ncrms)
  call allocate_grid(ncrms)
  call allocate_tracers(ncrms)
  call allocate_sgs(ncrms)
  call allocate_micro(ncrms)
#ifdef sam1mom
  call allocate_micro_params(ncrms)
#endif
#if defined(MMF_ESMT)
  call allocate_scalar_momentum(ncrms)
#endif

  crm_accel_ceaseflag = .false.

  !-----------------------------------------------

  dostatis  = .false.    ! no statistics are collected.
  idt_gl    = 1._r8/dt_gl
  ptop      = plev-nzm+1
  factor_xy = 1._r8/dble(nx*ny)
  crm_rad_temperature = 0.
  crm_rad_qv  = 0.
  crm_rad_qc  = 0.
  crm_rad_qi  = 0.
  crm_rad_cld = 0.
#ifdef m2005
  crm_rad_nc = 0.0
  crm_rad_ni = 0.0
  crm_rad_qs = 0.0
  crm_rad_ns = 0.0
#endif /* m2005 */
  do icrm = 1 , ncrms
    bflx(icrm) = crm_input_bflxls(icrm)
    wnd (icrm) = crm_input_wndls (icrm)
  enddo

!-----------------------------------------

  call task_init ()
  call setparm()

  do icrm = 1 , ncrms
    fcor(icrm)= 4*pi/86400.D0*sin(latitude0(icrm)*pi/180.D0)
    fcorz(icrm) = sqrt(4.D0*(2*pi/(3600.D0*24.D0))**2-fcor(icrm)**2)
    fcory(icrm,:) = fcor(icrm)
    fcorzy(icrm,:) = fcorz(icrm)
    do j=1,ny
      do i=1,nx
        latitude(icrm,i,j) = latitude0(icrm)
        longitude(icrm,i,j) = longitude0(icrm)
      end do
    end do

    if(crm_input_ocnfrac(icrm).gt.0.5D0) then
       OCEAN(icrm) = .true.
    else
       LAND(icrm) = .true.
    end if

    ! Create CRM vertical grid and initialize some vertical reference arrays:
    do k = 1, nzm
      z(icrm,k) = crm_input_zmid(icrm,plev-k+1) - crm_input_zint(icrm,plev+1)
      zi(icrm,k) = crm_input_zint(icrm,plev-k+2)- crm_input_zint(icrm,plev+1)
      pres(icrm,k) = crm_input_pmid(icrm,plev-k+1)/100.D0
      presi(icrm,k) = crm_input_pint(icrm,plev-k+2)/100.D0
      prespot(icrm,k)=(1000.D0/pres(icrm,k))**(rgas/cp)
      bet(icrm,k) = ggr/crm_input_tl(icrm,plev-k+1)
      gamaz(icrm,k)=ggr/cp*z(icrm,k)
    end do ! k
   ! zi(icrm,nz) =  crm_input_zint(plev-nz+2)
    zi(icrm,nz) = crm_input_zint(icrm,plev-nz+2)-crm_input_zint(icrm,plev+1) !+++mhwang, 2012-02-04
    presi(icrm,nz) = crm_input_pint(icrm, plev-nz+2)/100.D0

    dz(icrm) = 0.5D0*(z(icrm,1)+z(icrm,2))
    do k=2,nzm
      adzw(icrm,k) = (z(icrm,k)-z(icrm,k-1))/dz(icrm)
    end do
    adzw(icrm,1)  = 1.
    adzw(icrm,nz) = adzw(icrm,nzm)
    !+++mhwang fix the adz bug. (adz needs to be consistent with zi)
    !2012-02-04 Minghuai Wang (minghuai.wang@pnnl.gov)
    do k=1, nzm
      adz(icrm,k)=(zi(icrm,k+1)-zi(icrm,k))/dz(icrm)
    end do

    do k = 1,nzm
      rho(icrm,k) = crm_input_pdel(icrm,plev-k+1)/ggr/(adz(icrm,k)*dz(icrm))
    end do
    do k=2,nzm
    ! rhow(icrm,k) = 0.5*(rho(icrm,k)+rho(icrm,k-1))
    !+++mhwang fix the rhow bug (rhow needes to be consistent with crm_input_pmid)
    !2012-02-04 Minghuai Wang (minghuai.wang@pnnl.gov)
      rhow(icrm,k) = (crm_input_pmid(icrm,plev-k+2)-crm_input_pmid(icrm,plev-k+1))/ggr/(adzw(icrm,k)*dz(icrm))
    end do
    rhow(icrm,1) = 2.D0*rhow(icrm,2) - rhow(icrm,3)
    rhow(icrm,nz)= 2.D0*rhow(icrm,nzm) - rhow(icrm,nzm-1)
  enddo

  call t_startf('crm_gpu_region')
  !$omp taskwait
  call update_device_crm()
  call update_device_grid()
  call update_device_params()
  call update_device_vars()

  ! Initialize clear air relative humidity for aerosol water uptake
  !$omp target teams distribute parallel do collapse(2)
  do k = 1, nzm
    do icrm = 1, ncrms
      crm_clear_rh(icrm,k) = 0
      crm_clear_rh_cnt(icrm,k) = 0
    end do
  end do

  !  Initialize CRM fields:
  !$omp target teams distribute parallel do collapse(4)
  do k = 1 , nzm
    do j = 1 , ny
      do i = 1 , nx
        do icrm = 1 , ncrms
          u   (icrm,i,j,k) = crm_state_u_wind     (icrm,i,j,k)
#ifdef MAML
          !open the crm v component
          v   (icrm,i,j,k) = crm_state_v_wind     (icrm,i,j,k)
#else       
          v   (icrm,i,j,k) = crm_state_v_wind     (icrm,i,j,k)*YES3D
#endif
          w   (icrm,i,j,k) = crm_state_w_wind     (icrm,i,j,k)
          tabs(icrm,i,j,k) = crm_state_temperature(icrm,i,j,k)
#if defined(MMF_ESMT)
          u_esmt(icrm,i,j,k) = crm_input_ul_esmt(icrm,plev-k+1)
          v_esmt(icrm,i,j,k) = crm_input_vl_esmt(icrm,plev-k+1)
#endif /* MMF_ESMT */
        enddo
      enddo
    enddo
  enddo

  ! limit the velocity at the very first step:
  if(u(1,1,1,1).eq.u(1,2,1,1).and.u(1,3,1,2).eq.u(1,4,1,2)) then
    !$omp target teams distribute parallel do collapse(4)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm=1,ncrms
            u(icrm,i,j,k) = min( umax, max(-umax,u(icrm,i,j,k)) )
#ifdef MAML
            !open the crm v component
            v(icrm,i,j,k) = min( umax, max(-umax,v(icrm,i,j,k)) ) 
#else     
            v(icrm,i,j,k) = min( umax, max(-umax,v(icrm,i,j,k)) )*YES3D
#endif
          enddo
        enddo
      enddo
    enddo
  endif
 
  ! Populate microphysics array from crm_state
  !$omp target teams distribute parallel do collapse(4)
  do k = 1 , nzm
    do j = 1 , ny
      do i = 1 , nx
        do icrm = 1 , ncrms
#ifdef m2005
          micro_field(icrm,i,j,k,1 )  = crm_state%qt(icrm,i,j,k)
          micro_field(icrm,i,j,k,2 )  = crm_state%nc(icrm,i,j,k)
          micro_field(icrm,i,j,k,3 )  = crm_state%qr(icrm,i,j,k)
          micro_field(icrm,i,j,k,4 )  = crm_state%nr(icrm,i,j,k)
          micro_field(icrm,i,j,k,5 )  = crm_state%qi(icrm,i,j,k)
          micro_field(icrm,i,j,k,6 )  = crm_state%ni(icrm,i,j,k)
          micro_field(icrm,i,j,k,7 )  = crm_state%qs(icrm,i,j,k)
          micro_field(icrm,i,j,k,8 )  = crm_state%ns(icrm,i,j,k)
          micro_field(icrm,i,j,k,9 )  = crm_state%qg(icrm,i,j,k)
          micro_field(icrm,i,j,k,10)  = crm_state%ng(icrm,i,j,k)
          cloudliq   (icrm,i,j,k)     = crm_state%qc(icrm,i,j,k)
#else
          micro_field(icrm,i,j,k,1) = crm_state_qt(icrm,i,j,k)
          micro_field(icrm,i,j,k,2) = crm_state_qp(icrm,i,j,k)
          qn         (icrm,i,j,k)   = crm_state_qn(icrm,i,j,k)
#endif
        enddo
      enddo
    enddo
  enddo

#ifdef m2005
  do icrm = 1 , ncrms
    do k=1, nzm
#ifdef MODAL_AERO
      ! set aerosol data
      l=plev-k+1
      naer(icrm,k, 1:ntot_amode) = crm_input_naermod (icrm,l, 1:ntot_amode)
      vaer(icrm,k, 1:ntot_amode) = crm_input_vaerosol(icrm,l, 1:ntot_amode)
      hgaer(icrm,k, 1:ntot_amode) = crm_input_hygro   (icrm,l, 1:ntot_amode)
#endif /* MODAL_AERO */
      do j=1, ny
        do i=1, nx
          if(cloudliq(icrm,i,j,k).gt.0) then
            if(dopredictNc) then
              if( micro_field(icrm,i,j,k,incl).eq.0) micro_field(icrm,i,j,k,incl) = 1.0D6*Nc0/rho(icrm,k)
            endif
          endif
        enddo
      enddo
    enddo
  enddo
#endif /* m2005 */

  call micro_init(ncrms)

  ! initialize sgs fields
  call sgs_init(ncrms)

  !$omp target teams distribute parallel do
  do icrm = 1 , ncrms
    colprec (icrm)=0
    colprecs(icrm)=0
  enddo
  !$omp target teams distribute parallel do collapse(2)
  do k = 1 , nzm
    do icrm = 1 , ncrms
      u0   (icrm,k)=0.
      v0   (icrm,k)=0.
      t0   (icrm,k)=0.
      t00  (icrm,k)=0.
      tabs0(icrm,k)=0.
      q0   (icrm,k)=0.
      qv0  (icrm,k)=0.
      qn0  (icrm,k)=0.0
      qp0  (icrm,k)=0.0
      tke0 (icrm,k)=0.0
    enddo
  enddo

  !$omp target teams distribute parallel do collapse(4)
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          t(icrm,i,j,k) = tabs(icrm,i,j,k)+gamaz(icrm,k)-fac_cond*qcl(icrm,i,j,k)-fac_sub*qci(icrm,i,j,k) &
                                                        -fac_cond*qpl(icrm,i,j,k)-fac_sub*qpi(icrm,i,j,k)

          tmp = (qpl(icrm,i,j,k)+qpi(icrm,i,j,k))*crm_input_pdel(icrm,plev-k+1)
          !$omp atomic update
          colprec (icrm)=colprec (icrm)+tmp

          tmp = qpi(icrm,i,j,k)*crm_input_pdel(icrm,plev-k+1)

          !$omp atomic update
          colprecs(icrm)=colprecs(icrm)+tmp
          !$omp atomic update
          u0   (icrm,k)=u0   (icrm,k)+u(icrm,i,j,k)
          !$omp atomic update
          v0   (icrm,k)=v0   (icrm,k)+v(icrm,i,j,k)
          !$omp atomic update
          t0   (icrm,k)=t0   (icrm,k)+t(icrm,i,j,k)

          tmp = t(icrm,i,j,k)+fac_cond*qpl(icrm,i,j,k)+fac_sub*qpi(icrm,i,j,k)
          !$omp atomic update
          t00  (icrm,k)=t00  (icrm,k)+tmp
          !$omp atomic update
          tabs0(icrm,k)=tabs0(icrm,k)+tabs(icrm,i,j,k)

          tmp = qv(icrm,i,j,k)+qcl(icrm,i,j,k)+qci(icrm,i,j,k)
          !$omp atomic update
          q0   (icrm,k)=q0   (icrm,k)+tmp
          !$omp atomic update
          qv0  (icrm,k)=qv0  (icrm,k)+qv(icrm,i,j,k)

          tmp = qcl(icrm,i,j,k) + qci(icrm,i,j,k)
          !$omp atomic update
          qn0  (icrm,k)=qn0  (icrm,k)+tmp

          tmp = qpl(icrm,i,j,k) + qpi(icrm,i,j,k)
          !$omp atomic update
          qp0  (icrm,k)=qp0  (icrm,k)+tmp
          !$omp atomic update
          tke0 (icrm,k)=tke0 (icrm,k)+sgs_field(icrm,i,j,k,1)
        enddo
      enddo
    enddo
  enddo

  !$omp target teams distribute parallel do collapse(2)
  do k=1,nzm
    do icrm = 1 , ncrms
      u0   (icrm,k) = u0   (icrm,k) * factor_xy
      v0   (icrm,k) = v0   (icrm,k) * factor_xy
      t0   (icrm,k) = t0   (icrm,k) * factor_xy
      t00  (icrm,k) = t00  (icrm,k) * factor_xy
      tabs0(icrm,k) = tabs0(icrm,k) * factor_xy
      q0   (icrm,k) = q0   (icrm,k) * factor_xy
      qv0  (icrm,k) = qv0  (icrm,k) * factor_xy
      qn0  (icrm,k) = qn0  (icrm,k) * factor_xy
      qp0  (icrm,k) = qp0  (icrm,k) * factor_xy
      tke0 (icrm,k) = tke0 (icrm,k) * factor_xy
      l = plev-k+1
      uln  (icrm,l) = min( umax, max(-umax,crm_input_ul(icrm,l)) )
#ifdef MAML
      !open the crm v component
      vln  (icrm,l) = min( umax, max(-umax,crm_input_vl(icrm,l)) )
#else
      vln  (icrm,l) = min( umax, max(-umax,crm_input_vl(icrm,l)) )*YES3D
#endif
      ttend(icrm,k) = (crm_input_tl(icrm,l)+gamaz(icrm,k)- fac_cond*(crm_input_qccl(icrm,l)+crm_input_qiil(icrm,l))-fac_fus*crm_input_qiil(icrm,l)-t00(icrm,k))*idt_gl
      qtend(icrm,k) = (crm_input_ql(icrm,l)+crm_input_qccl(icrm,l)+crm_input_qiil(icrm,l)-q0(icrm,k))*idt_gl
      utend(icrm,k) = (uln(icrm,l)-u0(icrm,k))*idt_gl
      vtend(icrm,k) = (vln(icrm,l)-v0(icrm,k))*idt_gl
      ug0  (icrm,k) = uln(icrm,l)
      vg0  (icrm,k) = vln(icrm,l)
      tg0  (icrm,k) = crm_input_tl(icrm,l)+gamaz(icrm,k)-fac_cond*crm_input_qccl(icrm,l)-fac_sub*crm_input_qiil(icrm,l)
      qg0  (icrm,k) = crm_input_ql(icrm,l)+crm_input_qccl(icrm,l)+crm_input_qiil(icrm,l)
    end do ! k
  end do ! icrm

  !$omp target teams distribute parallel do 
  do icrm = 1 , ncrms
    uhl(icrm) = u0(icrm,1)
    vhl(icrm) = v0(icrm,1)
    ! estimate roughness length assuming logarithmic profile of velocity near the surface:
    ustar(icrm) = sqrt(crm_input_tau00(icrm)/rho(icrm,1))
    z0(icrm) = z0_est(z(icrm,1),bflx(icrm),wnd(icrm),ustar(icrm))
    z0(icrm) = max(real(0.00001D0,crm_rknd),min(real(1.,crm_rknd),z0(icrm)))
    crm_output_subcycle_factor(icrm) = 0.
    crm_output_prectend (icrm)=colprec (icrm)
    crm_output_precstend(icrm)=colprecs(icrm)
  enddo

!---------------------------------------------------
#ifdef m2005
  crm_output_nc_mean = 0.
  crm_output_ni_mean = 0.
  crm_output_ns_mean = 0.
  crm_output_ng_mean = 0.
  crm_output_nr_mean = 0.
  crm_output_aut_a  = 0.
  crm_output_acc_a  = 0.
  crm_output_evpc_a = 0.
  crm_output_evpr_a = 0.
  crm_output_mlt_a  = 0.
  crm_output_sub_a  = 0.
  crm_output_dep_a  = 0.
  crm_output_con_a  = 0.
  aut1a  = 0.
  acc1a  = 0.
  evpc1a = 0.
  evpr1a = 0.
  mlt1a  = 0.
  sub1a  = 0.
  dep1a  = 0.
  con1a  = 0.
#endif /* m2005 */
  !$omp target teams distribute parallel do collapse(2)
  do k = 1 , plev+1
    do icrm = 1 , ncrms
      if (k <= plev) then
        crm_output_cld       (icrm,k) = 0.
        crm_output_cldtop    (icrm,k) = 0.
        crm_output_gicewp    (icrm,k) = 0
        crm_output_gliqwp    (icrm,k) = 0
        crm_output_mctot     (icrm,k) = 0.
        crm_output_mcup      (icrm,k) = 0.
        crm_output_mcdn      (icrm,k) = 0.
        crm_output_mcuup     (icrm,k) = 0.
        crm_output_mcudn     (icrm,k) = 0.
        crm_output_qc_mean   (icrm,k) = 0.
        crm_output_qi_mean   (icrm,k) = 0.
        crm_output_qs_mean   (icrm,k) = 0.
        crm_output_qg_mean   (icrm,k) = 0.
        crm_output_qr_mean   (icrm,k) = 0.
        crm_output_mu_crm    (icrm,k) = 0.
        crm_output_md_crm    (icrm,k) = 0.
        crm_output_eu_crm    (icrm,k) = 0.
        crm_output_du_crm    (icrm,k) = 0.
        crm_output_ed_crm    (icrm,k) = 0.
        crm_output_flux_qt   (icrm,k) = 0.
        crm_output_flux_u    (icrm,k) = 0.
        crm_output_flux_v    (icrm,k) = 0.
        crm_output_fluxsgs_qt(icrm,k) = 0.
        crm_output_tkez      (icrm,k) = 0.
        crm_output_tkew      (icrm,k) = 0.
        crm_output_tkesgsz   (icrm,k) = 0.
        crm_output_tkz       (icrm,k) = 0.
        crm_output_flux_qp   (icrm,k) = 0.
        crm_output_precflux  (icrm,k) = 0.
        crm_output_qt_trans  (icrm,k) = 0.
        crm_output_qp_trans  (icrm,k) = 0.
        crm_output_qp_fall   (icrm,k) = 0.
        crm_output_qp_evp    (icrm,k) = 0.
        crm_output_qp_src    (icrm,k) = 0.
        crm_output_qt_ls     (icrm,k) = 0.
        crm_output_t_ls      (icrm,k) = 0.
        dd_crm               (icrm,k) = 0.
      endif
      mui_crm(icrm,k) = 0.
      mdi_crm(icrm,k) = 0.
    enddo
  enddo
  !$omp target teams distribute parallel do 
  do icrm = 1 , ncrms
    crm_output_jt_crm(icrm) = 0.
    crm_output_mx_crm(icrm) = 0.
  enddo

!--------------------------------------------------
#ifdef sam1mom
  if(doprecip) call precip_init(ncrms)
#endif
  !$omp taskwait

  do icrm = 1, ncrms
    if ( igstep <= 1 ) then
        iseed = gcolp(icrm) * perturb_seed_scale
        call setperturb(ncrms,icrm,iseed)
    end if

    !--------------------------
    ! sanity check for method to calculate radiation
    ! over averaged groups of columns instead of each individually
    if ( mod(nx,crm_nx_rad)==0 .or. mod(nx,crm_nx_rad)==0  ) then
      crm_nx_rad_fac = real(crm_nx_rad,crm_rknd)/real(nx,crm_rknd)
      crm_ny_rad_fac = real(crm_ny_rad,crm_rknd)/real(ny,crm_rknd)
    else
      write(0,*) "crm_nx_rad and crm_ny_rad need to be divisible by nx and ny"
#ifdef MMF_STANDALONE
      stop
#else
      call endrun('crm main')
#endif
    end if
  enddo

#ifdef MAML
  if(crm_nx_rad.NE.crm_nx .or. crm_ny_rad.NE.crm_ny) then 
     write(0,*) "crm_nx_rad and crm_ny_rad have to be equal to crm_nx and crm_ny in the MAML configuration"
     call endrun('crm main')
  end if
#endif

#ifdef ECPP
  call ecpp_crm_init(ncrms,dt_gl)

  qlsink    = 0.0
  qlsink_bf = 0.0
  prain     = 0.0
  precr     = 0.0
  precsolid = 0.0
#endif /* ECPP */

  nstop = dt_gl/dt
  dt = dt_gl/nstop

  crm_run_time  = dt_gl
  icrm_run_time = 1._r8/crm_run_time

  if (use_crm_accel) then
    call crm_accel_nstop(nstop)  ! reduce nstop by factor of (1 + crm_accel_factor)
  end if

  !========================================================================================
  !----------------------------------------------------------------------------------------
  !   Main time loop
  !----------------------------------------------------------------------------------------
  !========================================================================================
  nstep = 0
  do while (nstep < nstop)
    nstep = nstep + 1

    !------------------------------------------------------------------
    !  Check if the dynamical time step should be decreased
    !  to handle the cases when the flow being locally linearly unstable
    !------------------------------------------------------------------
    call kurant(ncrms)
    !$omp taskwait

    do icyc=1,ncycle
      icycle = icyc
      dtn = dt/ncycle
      dt3(na) = dtn
      dtfactor = dtn/dt

      !$omp target teams distribute parallel do 
      do icrm = 1 , ncrms
        crm_output_subcycle_factor(icrm) = crm_output_subcycle_factor(icrm)+1
      enddo

      !---------------------------------------------
      !  	the Adams-Bashforth scheme in time
      call abcoefs(ncrms)

      !---------------------------------------------
      !  	initialize stuff:
      call zero(ncrms)

      !-----------------------------------------------------------
      !       Buoyancy term:
      call buoyancy(ncrms)

      !------------------------------------------------------------
      !       Large-scale and surface forcing:
      call forcing(ncrms)

      ! Apply radiative tendency
      !omp target teams distribute parallel do collapse(4)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              i_rad = (i-1) / (nx/crm_nx_rad) + 1
              j_rad = (j-1) / (ny/crm_ny_rad) + 1
              t(icrm,i,j,k) = t(icrm,i,j,k) + crm_rad_qrad(icrm,i_rad,j_rad,k)*dtn
            enddo
          enddo
        enddo
      enddo

      !----------------------------------------------------------
      !   	suppress turbulence near the upper boundary (spange):
      if (dodamping) call damping(ncrms)

      !---------------------------------------------------------
      !   Ice fall-out
      if(docloud) then
        call ice_fall(ncrms)
      endif

      !----------------------------------------------------------
      !     Update scalar boundaries after large-scale processes:
      call boundaries(ncrms,3)

      !---------------------------------------------------------
      !     Update boundaries for velocities:
      call boundaries(ncrms,0)

      !-----------------------------------------------
      !     surface fluxes:
      if (dosurface) call crm_surface(ncrms,bflx)

      !-----------------------------------------------------------
      !  SGS physics:
      if (dosgs) call sgs_proc(ncrms)

      !----------------------------------------------------------
      !     Fill boundaries for SGS diagnostic fields:
      call boundaries(ncrms,4)

      !-----------------------------------------------
      !       advection of momentum:
      call advect_mom(ncrms)

      !----------------------------------------------------------
      !	SGS effects on momentum:
      if(dosgs) call sgs_mom(ncrms)

      !-----------------------------------------------------------
      !       Coriolis force:
      if (docoriolis) call coriolis(ncrms)

      !---------------------------------------------------------
      !       compute rhs of the Poisson equation and solve it for pressure.
      call pressure(ncrms)

      !---------------------------------------------------------
      !       find velocity field at n+1/2 timestep needed for advection of scalars:
      !  Note that at the end of the call, the velocities are in nondimensional form.
      call adams(ncrms)

      !----------------------------------------------------------
      !     Update boundaries for all prognostic scalar fields for advection:
      call boundaries(ncrms,2)

      !---------------------------------------------------------
      !      advection of scalars :
      call advect_all_scalars(ncrms)

      !-----------------------------------------------------------
      !    Convert velocity back from nondimensional form:
      call uvw(ncrms)

      !----------------------------------------------------------
      !     Update boundaries for scalars to prepare for SGS effects:
      call boundaries(ncrms,3)

      !---------------------------------------------------------
      !      SGS effects on scalars :
      if (dosgs) call sgs_scalars(ncrms)

      !-----------------------------------------------------------
      !       Calculate PGF for scalar momentum tendency
#if defined( MMF_ESMT ) 
      call scalar_momentum_tend(ncrms)
#endif

      !-----------------------------------------------------------
      !       Cloud condensation/evaporation and precipitation processes:
      if(docloud) call micro_proc(ncrms)

      !-----------------------------------------------------------
      !       Apply mean-state acceleration
      if (use_crm_accel .and. .not. crm_accel_ceaseflag) then
        ! Use Jones-Bretherton-Pritchard methodology to accelerate
        ! CRM horizontal mean evolution artificially.
        call accelerate_crm(ncrms, nstep, nstop, crm_accel_ceaseflag)
      endif

      !-----------------------------------------------------------
      !    Compute diagnostics fields:
      call diagnose(ncrms)

      !----------------------------------------------------------
      ! Rotate the dynamic tendency arrays for Adams-bashforth scheme:
      nn=na
      na=nc
      nc=nb
      nb=nn
    enddo ! icycle

#ifdef ECPP
    ! Here ecpp_crm_stat is called every CRM time step (dt), not every subcycle time step (dtn).
    ! This is what the original MMF model did (crm_rad_temperature, crm_rad_qv, ...). Do we want to call ecpp_crm_stat
    ! every subcycle time step??? +++mhwang
    call ecpp_crm_stat(ncrms)
#endif
    !$omp target teams distribute parallel do collapse(3)
    do j = 1 , ny
      do i = 1 , nx
        do icrm = 1 , ncrms
          cwp(icrm,i,j) = 0.
          cwph(icrm,i,j) = 0.
          cwpm(icrm,i,j) = 0.
          cwpl(icrm,i,j) = 0.

          flag_top(icrm,i,j) = .true.

          cltemp(icrm,i,j) = 0.0; cmtemp(icrm,i,j) = 0.0
          chtemp(icrm,i,j) = 0.0; cttemp(icrm,i,j) = 0.0
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(3)
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          do k=1,nzm
            l = plev-k+1
            tmp1 = rho(icrm,nz-k)*adz(icrm,nz-k)*dz(icrm)*(qcl(icrm,i,j,nz-k)+qci(icrm,i,j,nz-k))
            !$omp atomic update
            cwp(icrm,i,j) = cwp(icrm,i,j)+tmp1
            cttemp(icrm,i,j) = max(cf3d(icrm,i,j,nz-k), cttemp(icrm,i,j))
            if(cwp(icrm,i,j).gt.cwp_threshold.and.flag_top(icrm,i,j)) then
                !$omp atomic update
                crm_output_cldtop(icrm,l) = crm_output_cldtop(icrm,l) + 1
                flag_top(icrm,i,j) = .false.
            endif
            if(pres(icrm,nz-k).ge.700.D0) then
                !$omp atomic update
                cwpl(icrm,i,j) = cwpl(icrm,i,j)+tmp1
                cltemp(icrm,i,j) = max(cf3d(icrm,i,j,nz-k), cltemp(icrm,i,j))
            else if(pres(icrm,nz-k).lt.400.D0) then
                !$omp atomic update
                cwph(icrm,i,j) = cwph(icrm,i,j)+tmp1
                chtemp(icrm,i,j) = max(cf3d(icrm,i,j,nz-k), chtemp(icrm,i,j))
            else
                !$omp atomic update
                cwpm(icrm,i,j) = cwpm(icrm,i,j)+tmp1
                cmtemp(icrm,i,j) = max(cf3d(icrm,i,j,nz-k), cmtemp(icrm,i,j))
            endif
            tmp1 = rho(icrm,k)*adz(icrm,k)*dz(icrm)
            if(tmp1*(qcl(icrm,i,j,k)+qci(icrm,i,j,k)).gt.cwp_threshold) then
                 !$omp atomic update
                 crm_output_cld(icrm,l) = crm_output_cld(icrm,l) + cf3d(icrm,i,j,k)
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).gt.2*wmin) then
                   tmp = rho(icrm,k)*0.5D0*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * cf3d(icrm,i,j,k)
                   !$omp atomic update
                   crm_output_mcup (icrm,l) = crm_output_mcup (icrm,l) + tmp
                   tmp = rho(icrm,k)*0.5D0*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * (1.0D0 - cf3d(icrm,i,j,k))
                   !$omp atomic update
                   crm_output_mcuup(icrm,l) = crm_output_mcuup(icrm,l) + tmp
                 endif
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).lt.-2*wmin) then
                   tmp = rho(icrm,k)*0.5D0*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * cf3d(icrm,i,j,k)
                   !$omp atomic update
                   crm_output_mcdn (icrm,l) = crm_output_mcdn (icrm,l) + tmp
                   tmp = rho(icrm,k)*0.5D0*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * (1.D0 - cf3d(icrm,i,j,k))
                   !$omp atomic update
                   crm_output_mcudn(icrm,l) = crm_output_mcudn(icrm,l) + tmp
                 endif
            else
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).gt.2*wmin) then
                   tmp = rho(icrm,k)*0.5D0*(w(icrm,i,j,k+1)+w(icrm,i,j,k))
                   !$omp atomic update
                   crm_output_mcuup(icrm,l) = crm_output_mcuup(icrm,l) + tmp
                 endif
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).lt.-2*wmin) then
                    tmp = rho(icrm,k)*0.5D0*(w(icrm,i,j,k+1)+w(icrm,i,j,k))
                   !$omp atomic update
                   crm_output_mcudn(icrm,l) = crm_output_mcudn(icrm,l) + tmp
                 endif
            endif

            !$omp atomic update
            crm_output_gliqwp(icrm,l) = crm_output_gliqwp(icrm,l) + qcl(icrm,i,j,k)
            !$omp atomic update
            crm_output_gicewp(icrm,l) = crm_output_gicewp(icrm,l) + qci(icrm,i,j,k)
          enddo
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(4)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            ! Reduced radiation method allows for fewer radiation calculations
            ! by collecting statistics and doing radiation over column groups
            i_rad = (i-1) / (nx/crm_nx_rad) + 1
            j_rad = (j-1) / (ny/crm_ny_rad) + 1

            !$omp atomic update
            crm_rad_temperature(icrm,i_rad,j_rad,k) = crm_rad_temperature(icrm,i_rad,j_rad,k) + tabs(icrm,i,j,k)
            tmp = max(real(0.,crm_rknd),qv(icrm,i,j,k))
            !$omp atomic update
            crm_rad_qv         (icrm,i_rad,j_rad,k) = crm_rad_qv         (icrm,i_rad,j_rad,k) + tmp
            !$omp atomic update
            crm_rad_qc         (icrm,i_rad,j_rad,k) = crm_rad_qc         (icrm,i_rad,j_rad,k) + qcl(icrm,i,j,k)
            !$omp atomic update
            crm_rad_qi         (icrm,i_rad,j_rad,k) = crm_rad_qi         (icrm,i_rad,j_rad,k) + qci(icrm,i,j,k)
            if (qcl(icrm,i,j,k) + qci(icrm,i,j,k) > 0) then
               !$omp atomic update
               crm_rad_cld     (icrm,i_rad,j_rad,k) = crm_rad_cld        (icrm,i_rad,j_rad,k) + cf3d(icrm,i,j,k)
            else
               rh_tmp = qv(icrm,i,j,k)/qsatw_crm(tabs(icrm,i,j,k),pres(icrm,k))
               !$omp atomic update
               crm_clear_rh(icrm,k) = crm_clear_rh(icrm,k) + rh_tmp
               !$omp atomic update
               crm_clear_rh_cnt(icrm,k) = crm_clear_rh_cnt(icrm,k) + 1
            endif
#ifdef m2005
            !$omp atomic update
            crm_rad_nc         (icrm,i_rad,j_rad,k) = crm_rad_nc         (icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,incl)
            !$omp atomic update
            crm_rad_ni         (icrm,i_rad,j_rad,k) = crm_rad_ni         (icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,inci)
            !$omp atomic update
            crm_rad_qs         (icrm,i_rad,j_rad,k) = crm_rad_qs         (icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,iqs )
            !$omp atomic update
            crm_rad_ns         (icrm,i_rad,j_rad,k) = crm_rad_ns         (icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,ins )
#endif
          enddo
        enddo
      enddo
    enddo

    ! Diagnose mass fluxes to drive CAM's convective transport of tracers.
    ! definition of mass fluxes is taken from Xu et al., 2002, QJRMS.
    !$omp target teams distribute parallel do collapse(3)
    do j=1, ny
      do i=1, nx
        do icrm = 1 , ncrms
          do k=1, nzm+1
            l=plev+1-k+1
            if(w(icrm,i,j,k).gt.0.) then
              kx=max(1, k-1)
              qsat = qsatw_crm(tabs(icrm,i,j,kx),pres(icrm,kx))
              if(qcl(icrm,i,j,kx)+qci(icrm,i,j,kx).gt.min(real(1.D-5,crm_rknd),0.01D0*qsat)) then
                tmp = rhow(icrm,k)*w(icrm,i,j,k)
                !$omp atomic update
                mui_crm(icrm,l) = mui_crm(icrm,l)+tmp
              endif
            else if (w(icrm,i,j,k).lt.0.) then
              kx=min(k+1, nzm)
              qsat = qsatw_crm(tabs(icrm,i,j,kx),pres(icrm,kx))
              if(qcl(icrm,i,j,kx)+qci(icrm,i,j,kx).gt.min(real(1.D-5,crm_rknd),0.01D0*qsat)) then
                tmp = rhow(icrm,k)*w(icrm,i,j,k)
                !$omp atomic update
                mdi_crm(icrm,l) = mdi_crm(icrm,l)+tmp
              else if(qpl(icrm,i,j,kx)+qpi(icrm,i,j,kx).gt.1.0D-4) then
                tmp = rhow(icrm,k)*w(icrm,i,j,k)
                !$omp atomic update
                mdi_crm(icrm,l) = mdi_crm(icrm,l)+tmp
              endif
            endif
          enddo
        enddo
      enddo
    enddo

    !$omp target teams distribute parallel do collapse(3)
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          if(cwp(icrm,i,j).gt.cwp_threshold) then
            !$omp atomic update
            crm_output_cltot(icrm) = crm_output_cltot(icrm) + cttemp(icrm,i,j)
          endif
          if(cwph(icrm,i,j).gt.cwp_threshold) then
            !$omp atomic update
            crm_output_clhgh(icrm) = crm_output_clhgh(icrm) + chtemp(icrm,i,j)
          endif
          if(cwpm(icrm,i,j).gt.cwp_threshold) then
            !$omp atomic update
            crm_output_clmed(icrm) = crm_output_clmed(icrm) + cmtemp(icrm,i,j)
          endif
          if(cwpl(icrm,i,j).gt.cwp_threshold) then
            !$omp atomic update
            crm_output_cllow(icrm) = crm_output_cllow(icrm) + cltemp(icrm,i,j)
          endif
        enddo
      enddo
    enddo

  enddo ! nstep

  ! for time-averaging crm output statistics
  factor_xyt = factor_xy / real(nstop,crm_rknd) 

  !========================================================================================
  !----------------------------------------------------------------------------------------
  ! End main time loop
  !----------------------------------------------------------------------------------------
  !========================================================================================

  tmp1 = crm_nx_rad_fac * crm_ny_rad_fac / real(nstop,crm_rknd)

  !$omp target teams distribute parallel do collapse(4)
  do k=1,nzm
    do j=1,crm_ny_rad
      do i=1,crm_nx_rad
        do icrm = 1 , ncrms
          crm_rad_temperature(icrm,i,j,k) = crm_rad_temperature(icrm,i,j,k) * tmp1
          crm_rad_qv         (icrm,i,j,k) = crm_rad_qv         (icrm,i,j,k) * tmp1
          crm_rad_qc         (icrm,i,j,k) = crm_rad_qc         (icrm,i,j,k) * tmp1
          crm_rad_qi         (icrm,i,j,k) = crm_rad_qi         (icrm,i,j,k) * tmp1
          crm_rad_cld        (icrm,i,j,k) = crm_rad_cld        (icrm,i,j,k) * tmp1
#ifdef m2005
          crm_rad_nc         (icrm,i,j,k) = crm_rad_nc         (icrm,i,j,k) * tmp1
          crm_rad_ni         (icrm,i,j,k) = crm_rad_ni         (icrm,i,j,k) * tmp1
          crm_rad_qs         (icrm,i,j,k) = crm_rad_qs         (icrm,i,j,k) * tmp1
          crm_rad_ns         (icrm,i,j,k) = crm_rad_ns         (icrm,i,j,k) * tmp1
#endif /* m2005 */
        enddo
      enddo
    enddo
  enddo

  !$omp target teams distribute parallel do collapse(2)
  do k=1,nzm
    do icrm = 1 , ncrms
      if (crm_clear_rh_cnt(icrm,k)>0) then
        crm_clear_rh(icrm,k) = crm_clear_rh(icrm,k) / crm_clear_rh_cnt(icrm,k)
      endif
    enddo
  enddo

  ! no CRM tendencies above its top
  !$omp target teams distribute parallel do collapse(2)
  do k = 1 , ptop-1
    do icrm = 1 , ncrms
      tln  (icrm,k) = crm_input_tl  (icrm,k)
      qln  (icrm,k) = crm_input_ql  (icrm,k)
      qccln(icrm,k) = crm_input_qccl(icrm,k)
      qiiln(icrm,k) = crm_input_qiil(icrm,k)
      uln  (icrm,k) = crm_input_ul  (icrm,k)
      vln  (icrm,k) = crm_input_vl  (icrm,k)
    enddo
  enddo

  !  Compute tendencies due to CRM:
  !$omp target teams distribute parallel do collapse(2)
  do k = ptop,plev
    do icrm = 1 , ncrms
      tln  (icrm,k) = 0.
      qln  (icrm,k) = 0.
      qccln(icrm,k) = 0.
      qiiln(icrm,k) = 0.
      uln  (icrm,k) = 0.
      vln  (icrm,k) = 0.
    enddo
  enddo
  !$omp target teams distribute parallel do
  do icrm = 1 , ncrms
    colprec (icrm)=0
    colprecs(icrm)=0
  enddo

#if defined( MMF_ESMT )
  ! Initialize updated domain mean momentum scalars to input values above CRM 
  !$omp target teams distribute parallel do collapse(2)
  do k = 1,ptop-1
    do icrm = 1 , ncrms
      uln_esmt(icrm,k)  = crm_input_ul_esmt(icrm,k)
      vln_esmt(icrm,k)  = crm_input_vl_esmt(icrm,k)
    end do
  end do
  ! Initialize updated domain mean momentum scalars to zero for CRM levels
  !$omp target teams distribute parallel do collapse(2)
  do k = ptop,plev
    do icrm = 1 , ncrms
      uln_esmt(icrm,k) = 0.
      vln_esmt(icrm,k) = 0.
    end do
  end do
#endif /* MMF_ESMT */

  !$omp target teams distribute parallel do collapse(4)
  do k = 1,nzm
    do i=1,nx
      do j=1,ny
        do icrm=1,ncrms
          l = plev-k+1

          tmp = (qpl(icrm,i,j,k)+qpi(icrm,i,j,k))*crm_input_pdel(icrm,plev-k+1)
          !$omp atomic update
          colprec (icrm)= colprec (icrm)+tmp

          tmp = qpi(icrm,i,j,k)*crm_input_pdel(icrm,plev-k+1)
          !$omp atomic update
          colprecs(icrm)= colprecs(icrm)+tmp
          !$omp atomic update
          tln(icrm,l)  = tln(icrm,l)  +tabs(icrm,i,j,k)
          !$omp atomic update
          qln(icrm,l)  = qln(icrm,l)  +qv(icrm,i,j,k)
          !$omp atomic update
          qccln(icrm,l)= qccln(icrm,l)+qcl(icrm,i,j,k)
          !$omp atomic update
          qiiln(icrm,l)= qiiln(icrm,l)+qci(icrm,i,j,k)
          !$omp atomic update
          uln(icrm,l)  = uln(icrm,l)  +u(icrm,i,j,k)
          !$omp atomic update
          vln(icrm,l)  = vln(icrm,l)  +v(icrm,i,j,k)

#if defined(MMF_ESMT)
          !$omp atomic update
          uln_esmt(icrm,l) = uln_esmt(icrm,l)+u_esmt(icrm,i,j,k)
          !$omp atomic update
          vln_esmt(icrm,l) = vln_esmt(icrm,l)+v_esmt(icrm,i,j,k)
#endif
        enddo ! j
      enddo ! i
    enddo ! k
  enddo ! icrm

#if defined(MMF_ESMT)
  !$omp taskwait
  do icrm=1,ncrms
    uln_esmt(icrm,ptop:plev) = uln_esmt(icrm,ptop:plev) * factor_xy
    vln_esmt(icrm,ptop:plev) = vln_esmt(icrm,ptop:plev) * factor_xy

    crm_output_u_tend_esmt(icrm,:) = (uln_esmt(icrm,:) - crm_input_ul_esmt(icrm,:))*icrm_run_time
    crm_output_v_tend_esmt(icrm,:) = (vln_esmt(icrm,:) - crm_input_vl_esmt(icrm,:))*icrm_run_time

    ! don't use tendencies from two top CRM levels
    crm_output_u_tend_esmt(icrm,ptop:ptop+1) = 0.
    crm_output_v_tend_esmt(icrm,ptop:ptop+1) = 0.
  enddo
#endif

  !$omp target teams distribute parallel do collapse(2)
  do k = ptop , plev
    do icrm = 1 , ncrms
      tln  (icrm,k) = tln  (icrm,k) * factor_xy
      qln  (icrm,k) = qln  (icrm,k) * factor_xy
      qccln(icrm,k) = qccln(icrm,k) * factor_xy
      qiiln(icrm,k) = qiiln(icrm,k) * factor_xy
      uln  (icrm,k) = uln  (icrm,k) * factor_xy
      vln  (icrm,k) = vln  (icrm,k) * factor_xy
    enddo
  enddo

  !$omp target teams distribute parallel do collapse(2)
  do k = 1 , plev
    do icrm = 1 , ncrms
      crm_output_sltend (icrm,k) = cp * (tln  (icrm,k) - crm_input_tl  (icrm,k)) * icrm_run_time
      crm_output_qltend (icrm,k) =      (qln  (icrm,k) - crm_input_ql  (icrm,k)) * icrm_run_time
      crm_output_qcltend(icrm,k) =      (qccln(icrm,k) - crm_input_qccl(icrm,k)) * icrm_run_time
      crm_output_qiltend(icrm,k) =      (qiiln(icrm,k) - crm_input_qiil(icrm,k)) * icrm_run_time
#if defined(MMF_MOMENTUM_FEEDBACK)
      crm_output_ultend(icrm,k) =       (uln  (icrm,k) - crm_input_ul  (icrm,k)) * icrm_run_time
      crm_output_vltend(icrm,k) =       (vln  (icrm,k) - crm_input_vl  (icrm,k)) * icrm_run_time
#endif /* MMF_MOMENTUM_FEEDBACK */
    enddo
  enddo
  !$omp target teams distribute parallel do
  do icrm = 1 , ncrms
    crm_output_prectend (icrm) = (colprec (icrm)-crm_output_prectend (icrm))/ggr*factor_xy * icrm_run_time
    crm_output_precstend(icrm) = (colprecs(icrm)-crm_output_precstend(icrm))/ggr*factor_xy * icrm_run_time
  enddo

  ! don't use CRM tendencies from two crm top levels
  ! radiation tendencies are added back after the CRM call (see crm_physics_tend)
  !$omp target teams distribute parallel do collapse(2)
  do k = ptop,ptop+1
    do icrm = 1 , ncrms
      crm_output_sltend (icrm,k) = 0.
      crm_output_qltend (icrm,k) = 0.
      crm_output_qcltend(icrm,k) = 0.
      crm_output_qiltend(icrm,k) = 0.
#if defined(MMF_MOMENTUM_FEEDBACK)
      crm_output_ultend (icrm,k) = 0.
      crm_output_vltend (icrm,k) = 0.
#endif /* MMF_MOMENTUM_FEEDBACK */
    enddo
  enddo

  !-------------------------------------------------------------
  !
  ! Save the last step to the permanent core:
  !$omp target teams distribute parallel do collapse(4)
  do k = 1 , nzm
    do j = 1 , ny
      do i = 1 , nx
        do icrm = 1 , ncrms
          crm_state_u_wind     (icrm,i,j,k) = u   (icrm,i,j,k)
          crm_state_v_wind     (icrm,i,j,k) = v   (icrm,i,j,k)
          crm_state_w_wind     (icrm,i,j,k) = w   (icrm,i,j,k)
          crm_state_temperature(icrm,i,j,k) = tabs(icrm,i,j,k)
#ifdef m2005
          crm_state%qt(icrm,i,j,k) = micro_field(icrm,i,j,k,1 )
          crm_state%nc(icrm,i,j,k) = micro_field(icrm,i,j,k,2 )
          crm_state%qr(icrm,i,j,k) = micro_field(icrm,i,j,k,3 )
          crm_state%nr(icrm,i,j,k) = micro_field(icrm,i,j,k,4 )
          crm_state%qi(icrm,i,j,k) = micro_field(icrm,i,j,k,5 )
          crm_state%ni(icrm,i,j,k) = micro_field(icrm,i,j,k,6 )
          crm_state%qs(icrm,i,j,k) = micro_field(icrm,i,j,k,7 )
          crm_state%ns(icrm,i,j,k) = micro_field(icrm,i,j,k,8 )
          crm_state%qg(icrm,i,j,k) = micro_field(icrm,i,j,k,9 )
          crm_state%ng(icrm,i,j,k) = micro_field(icrm,i,j,k,10)
          crm_state%qc(icrm,i,j,k) = cloudliq   (icrm,i,j,k)
#else
          crm_state_qt(icrm,i,j,k) = micro_field(icrm,i,j,k,1)
          crm_state_qp(icrm,i,j,k) = micro_field(icrm,i,j,k,2)
          crm_state_qn(icrm,i,j,k) = qn         (icrm,i,j,k)
#endif
          crm_output_tk (icrm,i,j,k) = sgs_field_diag(icrm,i,j,k,1)
          crm_output_tkh(icrm,i,j,k) = sgs_field_diag(icrm,i,j,k,2)
          crm_output_qcl (icrm,i,j,k) = qcl  (icrm,i,j,k)
          crm_output_qci (icrm,i,j,k) = qci  (icrm,i,j,k)
          crm_output_qpl (icrm,i,j,k) = qpl  (icrm,i,j,k)
          crm_output_qpi (icrm,i,j,k) = qpi  (icrm,i,j,k)
#ifdef m2005
          crm_output_wvar(icrm,i,j,k) = wvar (icrm,i,j,k)
          crm_output_aut (icrm,i,j,k) = aut1 (icrm,i,j,k)
          crm_output_acc (icrm,i,j,k) = acc1 (icrm,i,j,k)
          crm_output_evpc(icrm,i,j,k) = evpc1(icrm,i,j,k)
          crm_output_evpr(icrm,i,j,k) = evpr1(icrm,i,j,k)
          crm_output_mlt (icrm,i,j,k) = mlt1 (icrm,i,j,k)
          crm_output_sub (icrm,i,j,k) = sub1 (icrm,i,j,k)
          crm_output_dep (icrm,i,j,k) = dep1 (icrm,i,j,k)
          crm_output_con (icrm,i,j,k) = con1 (icrm,i,j,k)
#endif /* m2005 */
        enddo
      enddo
    enddo
  enddo
  !$omp target teams distribute parallel do
  do icrm = 1 , ncrms
    crm_output_z0m (icrm) = z0   (icrm)
    crm_output_taux(icrm) = taux0(icrm) / dble(nstop)
    crm_output_tauy(icrm) = tauy0(icrm) / dble(nstop)
  enddo

  !---------------------------------------------------------------
  !  Diagnostics:

  ! hm add 9/7/11, change from GCM-time step avg to end-of-timestep
  !$omp target teams distribute parallel do collapse(4)
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        do icrm=1,ncrms
          l = plev-k+1
          !$omp atomic update
          crm_output_qc_mean(icrm,l) = crm_output_qc_mean(icrm,l) + qcl(icrm,i,j,k)
          !$omp atomic update
          crm_output_qi_mean(icrm,l) = crm_output_qi_mean(icrm,l) + qci(icrm,i,j,k)
          !$omp atomic update
          crm_output_qr_mean(icrm,l) = crm_output_qr_mean(icrm,l) + qpl(icrm,i,j,k)
#ifdef sam1mom
          omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tgrmin)*a_gr))

          tmp = qpi(icrm,i,j,k)*omg
          !$omp atomic update
          crm_output_qg_mean(icrm,l) = crm_output_qg_mean(icrm,l) + tmp

          tmp = qpi(icrm,i,j,k)*(1.-omg)
          !$omp atomic update
          crm_output_qs_mean(icrm,l) = crm_output_qs_mean(icrm,l) + tmp
#else
          !$omp atomic update
          crm_output_qg_mean(icrm,l) = crm_output_qg_mean(icrm,l) + micro_field(icrm,i,j,k,iqg)
          !$omp atomic update
          crm_output_qs_mean(icrm,l) = crm_output_qs_mean(icrm,l) + micro_field(icrm,i,j,k,iqs)

          !$omp atomic update
          crm_output_nc_mean(icrm,l) = crm_output_nc_mean(icrm,l) + micro_field(icrm,i,j,k,incl)
          !$omp atomic update
          crm_output_ni_mean(icrm,l) = crm_output_ni_mean(icrm,l) + micro_field(icrm,i,j,k,inci)
          !$omp atomic update
          crm_output_nr_mean(icrm,l) = crm_output_nr_mean(icrm,l) + micro_field(icrm,i,j,k,inr )
          !$omp atomic update
          crm_output_ng_mean(icrm,l) = crm_output_ng_mean(icrm,l) + micro_field(icrm,i,j,k,ing )
          !$omp atomic update
          crm_output_ns_mean(icrm,l) = crm_output_ns_mean(icrm,l) + micro_field(icrm,i,j,k,ins )
#endif /* sam1mom */
        enddo
      enddo
    enddo
  enddo

  !$omp target teams distribute parallel do collapse(2)
  do k = 1 , plev
    do icrm = 1 , ncrms
      crm_output_cld   (icrm,k) = min( 1._r8, crm_output_cld   (icrm,k) * factor_xyt )
      crm_output_cldtop(icrm,k) = min( 1._r8, crm_output_cldtop(icrm,k) * factor_xyt )
      crm_output_gicewp(icrm,k) = crm_output_gicewp(icrm,k)*crm_input_pdel(icrm,k)*1000.D0/ggr * factor_xyt
      crm_output_gliqwp(icrm,k) = crm_output_gliqwp(icrm,k)*crm_input_pdel(icrm,k)*1000.D0/ggr * factor_xyt
      crm_output_mcup  (icrm,k) = crm_output_mcup (icrm,k) * factor_xyt
      crm_output_mcdn  (icrm,k) = crm_output_mcdn (icrm,k) * factor_xyt
      crm_output_mcuup (icrm,k) = crm_output_mcuup(icrm,k) * factor_xyt
      crm_output_mcudn (icrm,k) = crm_output_mcudn(icrm,k) * factor_xyt
      crm_output_mctot (icrm,k) = crm_output_mcup(icrm,k) + crm_output_mcdn(icrm,k) + crm_output_mcuup(icrm,k) + crm_output_mcudn(icrm,k)

      crm_output_qc_mean(icrm,k) = crm_output_qc_mean(icrm,k) * factor_xy
      crm_output_qi_mean(icrm,k) = crm_output_qi_mean(icrm,k) * factor_xy
      crm_output_qs_mean(icrm,k) = crm_output_qs_mean(icrm,k) * factor_xy
      crm_output_qg_mean(icrm,k) = crm_output_qg_mean(icrm,k) * factor_xy
      crm_output_qr_mean(icrm,k) = crm_output_qr_mean(icrm,k) * factor_xy
    enddo
  enddo

#ifdef m2005
  do icrm=1,ncrms
    crm_output_nc_mean(icrm,:) = crm_output_nc_mean(icrm,:) * factor_xy
    crm_output_ni_mean(icrm,:) = crm_output_ni_mean(icrm,:) * factor_xy
    crm_output_ns_mean(icrm,:) = crm_output_ns_mean(icrm,:) * factor_xy
    crm_output_ng_mean(icrm,:) = crm_output_ng_mean(icrm,:) * factor_xy
    crm_output_nr_mean(icrm,:) = crm_output_nr_mean(icrm,:) * factor_xy

    ! hm 8/31/11 new output, gcm-grid- and time-step avg
    ! add loop over i,j do get horizontal avg, and flip vertical array
    do k=1,nzm
      l = plev-k+1
      do j=1,ny
        do i=1,nx
          crm_output_aut_a (icrm,l) = crm_output_aut_a (icrm,l) + aut1a(icrm,i,j,k)
          crm_output_acc_a (icrm,l) = crm_output_acc_a (icrm,l) + acc1a(icrm,i,j,k)
          crm_output_evpc_a(icrm,l) = crm_output_evpc_a(icrm,l) + evpc1a(icrm,i,j,k)
          crm_output_evpr_a(icrm,l) = crm_output_evpr_a(icrm,l) + evpr1a(icrm,i,j,k)
          crm_output_mlt_a (icrm,l) = crm_output_mlt_a (icrm,l) + mlt1a(icrm,i,j,k)
          crm_output_sub_a (icrm,l) = crm_output_sub_a (icrm,l) + sub1a(icrm,i,j,k)
          crm_output_dep_a (icrm,l) = crm_output_dep_a (icrm,l) + dep1a(icrm,i,j,k)
          crm_output_con_a (icrm,l) = crm_output_con_a (icrm,l) + con1a(icrm,i,j,k)
        enddo
      enddo
    enddo

    ! note, rates are divded by dt to get mean rate over step
    crm_output_aut_a (icrm,:) = crm_output_aut_a (icrm,:) * factor_xyt / dt
    crm_output_acc_a (icrm,:) = crm_output_acc_a (icrm,:) * factor_xyt / dt
    crm_output_evpc_a(icrm,:) = crm_output_evpc_a(icrm,:) * factor_xyt / dt
    crm_output_evpr_a(icrm,:) = crm_output_evpr_a(icrm,:) * factor_xyt / dt
    crm_output_mlt_a (icrm,:) = crm_output_mlt_a (icrm,:) * factor_xyt / dt
    crm_output_sub_a (icrm,:) = crm_output_sub_a (icrm,:) * factor_xyt / dt
    crm_output_dep_a (icrm,:) = crm_output_dep_a (icrm,:) * factor_xyt / dt
    crm_output_con_a (icrm,:) = crm_output_con_a (icrm,:) * factor_xyt / dt
  enddo
#endif /* m2005 */

  !$omp target teams distribute parallel do
  do icrm = 1 , ncrms
    crm_output_precc (icrm) = 0.
    crm_output_precl (icrm) = 0.
    crm_output_precsc(icrm) = 0.
    crm_output_precsl(icrm) = 0.
  enddo

  !$omp target teams distribute parallel do collapse(3)
  do j=1,ny
    do i=1,nx
      do icrm = 1 , ncrms
#ifdef sam1mom
        precsfc(icrm,i,j) = precsfc(icrm,i,j)*dz(icrm)/dt/dble(nstop)
        precssfc(icrm,i,j) = precssfc(icrm,i,j)*dz(icrm)/dt/dble(nstop)
#endif /* sam1mom */
#ifdef m2005
        precsfc(icrm,i,j) = precsfc(icrm,i,j)*dz(icrm)/dt/dble(nstop)     !mm/s/dz --> mm/s
        precssfc(icrm,i,j) = precssfc(icrm,i,j)*dz(icrm)/dt/dble(nstop)   !mm/s/dz --> mm/s
#endif /* m2005 */

#ifdef MAML
        ! output CRM precip and snow to pass down individual CLM instances
        crm_output_crm_pcp(icrm,i,j) = precsfc(icrm,i,j)/1000.D0      ! mm/s --> m/s
        crm_output_crm_snw(icrm,i,j) = precssfc(icrm,i,j)/1000.D0     ! mm/s --> m/s
#endif

        if(precsfc(icrm,i,j).gt.10.D0/86400.D0) then
           !$omp atomic update
           crm_output_precc (icrm) = crm_output_precc (icrm) + precsfc(icrm,i,j)
           !$omp atomic update
           crm_output_precsc(icrm) = crm_output_precsc(icrm) + precssfc(icrm,i,j)
        else
           !$omp atomic update
           crm_output_precl (icrm) = crm_output_precl (icrm) + precsfc(icrm,i,j)
           !$omp atomic update
           crm_output_precsl(icrm) = crm_output_precsl(icrm) + precssfc(icrm,i,j)
        endif
      enddo
    enddo
  enddo

  !$omp target teams distribute parallel do collapse(3)
  do j = 1 , ny
    do i = 1 , nx
      do icrm = 1 , ncrms
        crm_output_prec_crm(icrm,i,j) = precsfc(icrm,i,j)/1000.D0           !mm/s --> m/s
      enddo
    enddo
  enddo

  !$omp target teams distribute parallel do
  do icrm = 1 , ncrms
    crm_output_precc (icrm) = crm_output_precc (icrm)*factor_xy/1000.D0
    crm_output_precl (icrm) = crm_output_precl (icrm)*factor_xy/1000.D0
    crm_output_precsc(icrm) = crm_output_precsc(icrm)*factor_xy/1000.D0
    crm_output_precsl(icrm) = crm_output_precsl(icrm)*factor_xy/1000.D0

    crm_output_cltot(icrm) = crm_output_cltot(icrm) * factor_xyt
    crm_output_clhgh(icrm) = crm_output_clhgh(icrm) * factor_xyt
    crm_output_clmed(icrm) = crm_output_clmed(icrm) * factor_xyt
    crm_output_cllow(icrm) = crm_output_cllow(icrm) * factor_xyt

    crm_output_jt_crm(icrm) = plev * 1.0
    crm_output_mx_crm(icrm) = 1.0
  enddo

  !$omp target teams distribute parallel do collapse(2)
  do k=1, plev
    do icrm = 1 , ncrms
      crm_output_mu_crm(icrm,k)=0.5D0*(mui_crm(icrm,k)+mui_crm(icrm,k+1))
      crm_output_md_crm(icrm,k)=0.5D0*(mdi_crm(icrm,k)+mdi_crm(icrm,k+1))
      crm_output_mu_crm(icrm,k)=crm_output_mu_crm(icrm,k)*ggr/100.D0          !kg/m2/s --> mb/s
      crm_output_md_crm(icrm,k)=crm_output_md_crm(icrm,k)*ggr/100.D0          !kg/m2/s --> mb/s
      crm_output_eu_crm(icrm,k) = 0.
      if(mui_crm(icrm,k)-mui_crm(icrm,k+1).gt.0) then
        crm_output_eu_crm(icrm,k)=(mui_crm(icrm,k)-mui_crm(icrm,k+1))*ggr/crm_input_pdel(icrm,k)    !/s
      else
        crm_output_du_crm(icrm,k)=-1.0*(mui_crm(icrm,k)-mui_crm(icrm,k+1))*ggr/crm_input_pdel(icrm,k)   !/s
      endif
      if(mdi_crm(icrm,k+1)-mdi_crm(icrm,k).lt.0) then
        crm_output_ed_crm(icrm,k)=(mdi_crm(icrm,k)-mdi_crm(icrm,k+1))*ggr/crm_input_pdel(icrm,k) ! /s
      else
        dd_crm(icrm,k)=-1.*(mdi_crm(icrm,k)-mdi_crm(icrm,k+1))*ggr/crm_input_pdel(icrm,k)   !/s
      endif
      if(abs(crm_output_mu_crm(icrm,k)).gt.1.0D-15.or.abs(crm_output_md_crm(icrm,k)).gt.1.0D-15) then
        tmp = k
        !$omp atomic update
        crm_output_jt_crm(icrm) = min( tmp , crm_output_jt_crm(icrm) )

        tmp = k
        !$omp atomic update
        crm_output_mx_crm(icrm) = max( tmp , crm_output_mx_crm(icrm) )
      endif
    enddo
  enddo

  !-------------------------------------------------------------
  !       Fluxes and other stat:
  !-------------------------------------------------------------
  !$omp target teams distribute parallel do collapse(2)
  do k=1,nzm
    do icrm = 1 , ncrms
      u2z = 0.
      v2z = 0.
      w2z = 0.
      do j=1,ny
        do i=1,nx
          u2z = u2z+(u(icrm,i,j,k)-u0(icrm,k))**2
          v2z = v2z+(v(icrm,i,j,k)-v0(icrm,k))**2
          w2z = w2z+0.5D0*(w(icrm,i,j,k+1)**2+w(icrm,i,j,k)**2)
        enddo
      enddo
      !+++mhwang
      ! mkwsb, mkle, mkadv, mkdiff (also crm_output_flux_u, crm_output_flux_v,icrm) seem not calculted correclty in the spcam3.5 codes.
      ! Only values at the last time step are calculated, but is averaged over the entire GCM
      ! time step.
      !---mhwang

      tmp1 = dz(icrm)/rhow(icrm,k)
      tmp2 = tmp1/dtn                        ! dtn is calculated inside of the icyc loop.
                                             ! It seems wrong to use it here ???? +++mhwang
      mkwsb(icrm,k,:) = mkwsb(icrm,k,:) * tmp1*rhow(icrm,k) * factor_xy/nstop     !kg/m3/s --> kg/m2/s
      mkwle(icrm,k,:) = mkwle(icrm,k,:) * tmp2*rhow(icrm,k) * factor_xy/nstop     !kg/m3   --> kg/m2/s
      mkadv(icrm,k,:) = mkadv(icrm,k,:) * factor_xy*icrm_run_time     ! kg/kg  --> kg/kg/s
      mkdiff(icrm,k,:) = mkdiff(icrm,k,:) * factor_xy*icrm_run_time   ! kg/kg  --> kg/kg/s

      ! qpsrc, qpevp, qpfall in M2005 are calculated in micro_flux.
      qpsrc(icrm,k) = qpsrc(icrm,k) * factor_xy*icrm_run_time
      qpevp(icrm,k) = qpevp(icrm,k) * factor_xy*icrm_run_time
      qpfall(icrm,k) = qpfall(icrm,k) * factor_xy*icrm_run_time   ! kg/kg in M2005 ---> kg/kg/s
      precflux(icrm,k) = precflux(icrm,k) * factor_xy*dz(icrm)/dt/nstop  !kg/m2/dz in M2005 -->kg/m2/s or mm/s (idt_gl=1/dt/nstop)

      l = plev-k+1
      crm_output_flux_u    (icrm,l) = (uwle(icrm,k) + uwsb(icrm,k))*tmp1*factor_xy/nstop
      crm_output_flux_v    (icrm,l) = (vwle(icrm,k) + vwsb(icrm,k))*tmp1*factor_xy/nstop
#ifdef sam1mom
      crm_output_flux_qt   (icrm,l) = mkwle(icrm,k,1) + mkwsb(icrm,k,1)
      crm_output_fluxsgs_qt(icrm,l) = mkwsb(icrm,k,1)
      crm_output_flux_qp   (icrm,l) = mkwle(icrm,k,2) + mkwsb(icrm,k,2)
      crm_output_qt_trans  (icrm,l) = mkadv(icrm,k,1) + mkdiff(icrm,k,1)
      crm_output_qp_trans  (icrm,l) = mkadv(icrm,k,2) + mkdiff(icrm,k,2)
#endif /* sam1mom */
#ifdef m2005
      crm_output_flux_qt   (icrm,l) = mkwle(icrm,k,1   ) + mkwsb(icrm,k,1   ) +  &
                         mkwle(icrm,k,iqci) + mkwsb(icrm,k,iqci)
      crm_output_fluxsgs_qt(icrm,l) = mkwsb(icrm,k,1   ) + mkwsb(icrm,k,iqci)
      crm_output_flux_qp   (icrm,l) = mkwle(icrm,k,iqr) + mkwsb(icrm,k,iqr) +  &
                         mkwle(icrm,k,iqs) + mkwsb(icrm,k,iqs) + mkwle(icrm,k,iqg) + mkwsb(icrm,k,iqg)
      crm_output_qt_trans  (icrm,l) = mkadv(icrm,k,1) + mkadv(icrm,k,iqci) + &
                         mkdiff(icrm,k,1) + mkdiff(icrm,k,iqci)
      crm_output_qp_trans  (icrm,l) = mkadv(icrm,k,iqr) + mkadv(icrm,k,iqs) + mkadv(icrm,k,iqg) + &
                         mkdiff(icrm,k,iqr) + mkdiff(icrm,k,iqs) + mkdiff(icrm,k,iqg)
#endif /* m2005 */
      tmp = 0
      do j = 1 , ny
        do i = 1 , nx
          tmp = tmp + sgs_field(icrm,i,j,k,1)
        enddo
      enddo
      crm_output_tkesgsz   (icrm,l)= rho(icrm,k)*tmp*factor_xy
      crm_output_tkez      (icrm,l)= rho(icrm,k)*0.5*(u2z+v2z*YES3D+w2z)*factor_xy + crm_output_tkesgsz(icrm,l)
      crm_output_tkew      (icrm,l)= rho(icrm,k)*0.5*w2z*factor_xy
      tmp = 0
      do j = 1 , ny
        do i = 1 , nx
          tmp = tmp + sgs_field_diag(icrm,i,j,k,1)
        enddo
      enddo
      crm_output_tkz       (icrm,l) = tmp * factor_xy
      crm_output_precflux  (icrm,l) = precflux(icrm,k)/1000.D0       !mm/s  -->m/s

      crm_output_qp_fall   (icrm,l) = qpfall(icrm,k)
      crm_output_qp_evp    (icrm,l) = qpevp(icrm,k)
      crm_output_qp_src    (icrm,l) = qpsrc(icrm,k)

      crm_output_qt_ls     (icrm,l) = qtend(icrm,k)
      crm_output_t_ls      (icrm,l) = ttend(icrm,k)
    enddo
  enddo

  call update_host_crm()
  !$omp taskwait

  call t_stopf('crm_gpu_region')

#ifdef ECPP
  do icrm = 1 , ncrms
    crm_ecpp_output%abnd         (icrm,:,:,:,:)=0.0
    crm_ecpp_output%abnd_tf      (icrm,:,:,:,:)=0.0
    crm_ecpp_output%massflxbnd   (icrm,:,:,:,:)=0.0
    crm_ecpp_output%acen         (icrm,:,:,:,:)=0.0
    crm_ecpp_output%acen_tf      (icrm,:,:,:,:)=0.0
    crm_ecpp_output%rhcen        (icrm,:,:,:,:)=0.0
    crm_ecpp_output%qcloudcen    (icrm,:,:,:,:)=0.0
    crm_ecpp_output%qicecen      (icrm,:,:,:,:)=0.0
    crm_ecpp_output%qlsinkcen    (icrm,:,:,:,:)=0.0
    crm_ecpp_output%precrcen     (icrm,:,:,:,:)=0.0
    crm_ecpp_output%precsolidcen (icrm,:,:,:,:)=0.0
    crm_ecpp_output%qlsink_bfcen (icrm,:,:,:,:)=0.0
    crm_ecpp_output%qlsink_avgcen(icrm,:,:,:,:)=0.0
    crm_ecpp_output%praincen     (icrm,:,:,:,:)=0.0

    crm_ecpp_output%wupthresh_bnd   (icrm,:)=0.0
    crm_ecpp_output%wdownthresh_bnd (icrm,:)=0.0
    crm_ecpp_output%wwqui_cen       (icrm,:)=0.0
    crm_ecpp_output%wwqui_bnd       (icrm,:)=0.0
    crm_ecpp_output%wwqui_cloudy_cen(icrm,:)=0.0
    crm_ecpp_output%wwqui_cloudy_bnd(icrm,:)=0.0

    ! default is clear, non-precipitating, and quiescent class
    crm_ecpp_output%abnd   (icrm,:,1,1,1)=1.0
    crm_ecpp_output%abnd_tf(icrm,:,1,1,1)=1.0
    crm_ecpp_output%acen   (icrm,:,1,1,1)=1.0
    crm_ecpp_output%acen_tf(icrm,:,1,1,1)=1.0

    do k=1, nzm
      l=plev-k+1
      crm_ecpp_output%acen            (icrm,l,:,:,:) = area_cen_sum        (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%acen_tf         (icrm,l,:,:,:) = area_cen_final      (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%rhcen           (icrm,l,:,:,:) = rh_cen_sum          (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%qcloudcen       (icrm,l,:,:,:) = qcloud_cen_sum      (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%qicecen         (icrm,l,:,:,:) = qice_cen_sum        (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%qlsinkcen       (icrm,l,:,:,:) = qlsink_cen_sum      (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%precrcen        (icrm,l,:,:,:) = precr_cen_sum       (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%precsolidcen    (icrm,l,:,:,:) = precsolid_cen_sum   (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%wwqui_cen       (icrm,l)       = wwqui_cen_sum       (k,icrm)
      crm_ecpp_output%wwqui_cloudy_cen(icrm,l)       = wwqui_cloudy_cen_sum(k,icrm)
      crm_ecpp_output%qlsink_bfcen    (icrm,l,:,:,:) = qlsink_bf_cen_sum   (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%qlsink_avgcen   (icrm,l,:,:,:) = qlsink_avg_cen_sum  (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%praincen        (icrm,l,:,:,:) = prain_cen_sum       (k,:,1:ncls_ecpp_in,:,icrm)
    enddo
    do k=1, nzm+1
      l=plev+1-k+1
      crm_ecpp_output%abnd            (icrm,l,:,:,:) = area_bnd_sum        (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%abnd_tf         (icrm,l,:,:,:) = area_bnd_final      (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%massflxbnd      (icrm,l,:,:,:) = mass_bnd_sum        (k,:,1:ncls_ecpp_in,:,icrm)
      crm_ecpp_output%wupthresh_bnd   (icrm,l)       = wup_thresh          (k,icrm)
      crm_ecpp_output%wdownthresh_bnd (icrm,l)       = wdown_thresh        (k,icrm)
      crm_ecpp_output%wwqui_bnd       (icrm,l)       = wwqui_bnd_sum       (k,icrm)
      crm_ecpp_output%wwqui_cloudy_bnd(icrm,l)       = wwqui_cloudy_bnd_sum(k,icrm)
    enddo
  enddo
#endif /* ECPP */

  crm_output_subcycle_factor(:) = crm_output_subcycle_factor(:) / nstop

#ifdef ECPP
  ! Deallocate ECPP variables
  call ecpp_crm_cleanup ()
#endif
  call deallocate_params()
  call deallocate_grid()
  call deallocate_tracers()
  call deallocate_sgs()
  call deallocate_vars()
  call deallocate_micro()
#ifdef sam1mom
  call deallocate_micro_params()
#endif
#if defined( MMF_ESMT )
  call deallocate_scalar_momentum()
#endif
  call deallocate_crm()

  CONTAINS
  subroutine allocate_crm()
    implicit none
    allocate( t00(ncrms,nz) )
    allocate( tln(ncrms,plev) )
    allocate( qln(ncrms,plev) )
    allocate( qccln(ncrms,plev) )
    allocate( qiiln(ncrms,plev) )
    allocate( uln(ncrms,plev) )
    allocate( vln(ncrms,plev) )
#if defined(MMF_ESMT)
    allocate( uln_esmt(ncrms,plev) )
    allocate( vln_esmt(ncrms,plev) )
#endif
    allocate( cwp(ncrms,nx,ny) )
    allocate( cwph(ncrms,nx,ny) )
    allocate( cwpm(ncrms,nx,ny) )
    allocate( cwpl(ncrms,nx,ny) )
    allocate( flag_top(ncrms,nx,ny) )
    allocate( cltemp(ncrms,nx,ny) )
    allocate( cmtemp(ncrms,nx,ny) )
    allocate( chtemp(ncrms,nx,ny) )
    allocate( cttemp(ncrms,nx,ny) )
    allocate( dd_crm (ncrms,plev)   )
    allocate( mui_crm(ncrms,plev+1) )
    allocate( mdi_crm(ncrms,plev+1) )
    allocate( ustar(ncrms) )
    allocate( bflx(ncrms) )
    allocate( wnd(ncrms) )
    allocate( qtot (ncrms,20) )
    allocate( colprec (ncrms) )
    allocate( colprecs(ncrms) )
    allocate( crm_clear_rh_cnt(ncrms,nzm) )

    crm_input_zmid    => crm_input%zmid(:,:)
    crm_input_zint    => crm_input%zint(:,:)
    crm_input_tl      => crm_input%tl(:,:)
    crm_input_ql      => crm_input%ql(:,:)
    crm_input_qccl    => crm_input%qccl(:,:)
    crm_input_qiil    => crm_input%qiil(:,:)
    crm_input_ps      => crm_input%ps(:)
    crm_input_pmid    => crm_input%pmid(:,:)
    crm_input_pint    => crm_input%pint(:,:)
    crm_input_pdel    => crm_input%pdel(:,:)
    crm_input_phis    => crm_input%phis(:)
    crm_input_ul      => crm_input%ul(:,:)
    crm_input_vl      => crm_input%vl(:,:)
    crm_input_ocnfrac => crm_input%ocnfrac(:)
    crm_input_tau00   => crm_input%tau00(:)
    crm_input_wndls   => crm_input%wndls(:)
    crm_input_bflxls  => crm_input%bflxls(:)
#if defined( m2005 ) && defined( MODAL_AERO )
    crm_input_naermod  => crm_input%naermod(:,:,:)
    crm_input_vaerosol => crm_input%vaerosol(:,:,:)
    crm_input_hygro    => crm_input%hygro(:,:,:)
#endif
#if defined( MMF_ESMT )
    crm_input_ul_esmt  => crm_input%ul_esmt(:,:)
    crm_input_vl_esmt  => crm_input%vl_esmt(:,:)
#endif
    crm_input_fluxu00 => crm_input%fluxu00(:)
    crm_input_fluxv00 => crm_input%fluxv00(:)
    crm_input_fluxt00 => crm_input%fluxt00(:)
    crm_input_fluxq00 => crm_input%fluxq00(:)

    crm_output_qcl    => crm_output%qcl(:,:,:,:)
    crm_output_qci    => crm_output%qci(:,:,:,:)
    crm_output_qpl    => crm_output%qpl(:,:,:,:)
    crm_output_qpi    => crm_output%qpi(:,:,:,:)

    crm_output_tk       => crm_output%tk (:,:,:,:)
    crm_output_tkh      => crm_output%tkh(:,:,:,:)
    crm_output_prec_crm => crm_output%prec_crm(:,:,:)
    crm_output_wvar     => crm_output%wvar(:,:,:,:)
    crm_output_aut      => crm_output%aut (:,:,:,:)
    crm_output_acc      => crm_output%acc (:,:,:,:)
    crm_output_evpc     => crm_output%evpc(:,:,:,:)
    crm_output_evpr     => crm_output%evpr(:,:,:,:)
    crm_output_mlt      => crm_output%mlt (:,:,:,:)
    crm_output_sub      => crm_output%sub (:,:,:,:)
    crm_output_dep      => crm_output%dep (:,:,:,:)
    crm_output_con      => crm_output%con (:,:,:,:)

    ! Cloud area fractions
    crm_output_cltot    => crm_output%cltot(:)
    crm_output_clhgh    => crm_output%clhgh(:)
    crm_output_clmed    => crm_output%clmed(:)
    crm_output_cllow    => crm_output%cllow(:)

    crm_output_cldtop   => crm_output%cldtop(:,:)
    crm_output_precc    => crm_output%precc(:)
    crm_output_precl    => crm_output%precl(:)
    crm_output_precsc   => crm_output%precsc(:)
    crm_output_precsl   => crm_output%precsl(:)

    crm_output_qc_mean  => crm_output%qc_mean(:,:)
    crm_output_qi_mean  => crm_output%qi_mean(:,:)
    crm_output_qs_mean  => crm_output%qs_mean(:,:)
    crm_output_qg_mean  => crm_output%qg_mean(:,:)
    crm_output_qr_mean  => crm_output%qr_mean(:,:)

#ifdef m2005
    crm_output_nc_mean  => crm_output%nc_mean(:,:)
    crm_output_ni_mean  => crm_output%ni_mean(:,:)
    crm_output_ns_mean  => crm_output%ns_mean(:,:)
    crm_output_ng_mean  => crm_output%ng_mean(:,:)
    crm_output_nr_mean  => crm_output%nr_mean(:,:)

    ! Time and domain averaged process rates
    crm_output_aut_a    => crm_output%aut_a (:,:)
    crm_output_acc_a    => crm_output%acc_a (:,:)
    crm_output_evpc_a   => crm_output%evpc_a(:,:)
    crm_output_evpr_a   => crm_output%evpr_a(:,:)
    crm_output_mlt_a    => crm_output%mlt_a (:,:)

    crm_output_sub_a    => crm_output%sub_a (:,:)
    crm_output_dep_a    => crm_output%dep_a (:,:)
    crm_output_con_a    => crm_output%con_a (:,:)
#endif /* m2005 */

#if defined( MMF_MOMENTUM_FEEDBACK )
    crm_output_ultend  => crm_output%ultend(:,:)
    crm_output_vltend  => crm_output%vltend(:,:)
#endif

#if defined( MMF_ESMT )
    crm_output_t_tend   => crm_output%u_tend_esmt(:,:)
    crm_output_v_tend   => crm_output%v_tend_esmt(:,:)
#endif
    crm_output_sltend  => crm_output%sltend  (:,:)
    crm_output_qltend  => crm_output%qltend  (:,:)
    crm_output_qcltend => crm_output%qcltend (:,:)
    crm_output_qiltend => crm_output%qiltend (:,:)

    ! These are all time and spatial averages, on the GCM grid
    crm_output_cld    => crm_output%cld   (:,:)
    crm_output_gicewp => crm_output%gicewp(:,:)
    crm_output_gliqwp => crm_output%gliqwp(:,:)
    crm_output_mctot  => crm_output%mctot (:,:)
    crm_output_mcup   => crm_output%mcup  (:,:)
    crm_output_mcdn   => crm_output%mcdn  (:,:)
    crm_output_mcuup  => crm_output%mcuup (:,:)
    crm_output_mcudn  => crm_output%mcudn (:,:)

    ! For convective transport
    crm_output_mu_crm  => crm_output%mu_crm(:,:)
    crm_output_md_crm  => crm_output%md_crm(:,:)
    crm_output_du_crm  => crm_output%du_crm(:,:)
    crm_output_eu_crm  => crm_output%eu_crm(:,:)
    crm_output_ed_crm  => crm_output%ed_crm(:,:)
    crm_output_jt_crm  => crm_output%jt_crm(:)
    crm_output_mx_crm  => crm_output%mx_crm(:)

    ! Other stuff...
    crm_output_flux_qt    => crm_output%flux_qt(:,:)
    crm_output_fluxsgs_qt => crm_output%fluxsgs_qt(:,:)
    crm_output_tkez       => crm_output%tkez(:,:)
    crm_output_tkew       => crm_output%tkew(:,:)
    crm_output_tkesgsz    => crm_output%tkesgsz(:,:)
    crm_output_tkz        => crm_output%tkz(:,:)
    crm_output_flux_u     => crm_output%flux_u(:,:)
    crm_output_flux_v     => crm_output%flux_v(:,:)
    crm_output_flux_qp    => crm_output%flux_qp(:,:)
    crm_output_precflux   => crm_output%precflux(:,:)
    crm_output_qt_ls      => crm_output%qt_ls(:,:)
    crm_output_qt_trans   => crm_output%qt_trans(:,:)
    crm_output_qp_trans   => crm_output%qp_trans(:,:)
    crm_output_qp_fall    => crm_output%qp_fall(:,:)
    crm_output_qp_src     => crm_output%qp_src(:,:)
    crm_output_qp_evp     => crm_output%qp_evp(:,:)
    crm_output_t_ls       => crm_output%t_ls(:,:)
    crm_output_prectend   => crm_output%prectend(:)
    crm_output_precstend  => crm_output%precstend(:)
    crm_output_taux       => crm_output%taux(:)
    crm_output_tauy       => crm_output%tauy(:)
    crm_output_z0m        => crm_output%z0m(:)
#ifdef MAML
    ! MAML variables
    crm_output_crm_pcp   => crm_output%crm_pcp(:,:,:)
    crm_output_crm_snw   => crm_output%crm_snw(:,:,:)
#endif
    crm_output_subcycle_factor => crm_output%subcycle_factor(:)

    crm_rad_temperature => crm_rad%temperature(1:ncrms,:,:,:)
    crm_rad_qv          => crm_rad%qv         (1:ncrms,:,:,:)
    crm_rad_qc          => crm_rad%qc         (1:ncrms,:,:,:)
    crm_rad_qi          => crm_rad%qi         (1:ncrms,:,:,:)
    crm_rad_cld         => crm_rad%cld        (1:ncrms,:,:,:)
    crm_rad_qrad        => crm_rad%qrad       (1:ncrms,:,:,:)

#ifdef m2005
    crm_rad_nc => crm_rad%nc(1:ncrms,:,:,:)
    crm_rad_ni => crm_rad%ni(1:ncrms,:,:,:)
    crm_rad_qs => crm_rad%qs(1:ncrms,:,:,:)
    crm_rad_ns => crm_rad%ns(1:ncrms,:,:,:)
#endif /* m2005 */

    crm_state_u_wind      => crm_state%u_wind     (1:ncrms,:,:,:)
    crm_state_v_wind      => crm_state%v_wind     (1:ncrms,:,:,:)
    crm_state_w_wind      => crm_state%w_wind     (1:ncrms,:,:,:)
    crm_state_temperature => crm_state%temperature(1:ncrms,:,:,:)
    crm_state_qt          => crm_state%qt         (1:ncrms,:,:,:)
    crm_state_qp          => crm_state%qp         (1:ncrms,:,:,:)
    crm_state_qn          => crm_state%qn         (1:ncrms,:,:,:)

    !$omp target enter data map(alloc: t00      )
    !$omp target enter data map(alloc: tln      )
    !$omp target enter data map(alloc: qln      )
    !$omp target enter data map(alloc: qccln    )
    !$omp target enter data map(alloc: qiiln    )
    !$omp target enter data map(alloc: uln      )
    !$omp target enter data map(alloc: vln      )
    !$omp target enter data map(alloc: cwp      )
    !$omp target enter data map(alloc: cwph     )
    !$omp target enter data map(alloc: cwpm     )
    !$omp target enter data map(alloc: cwpl     )
    !$omp target enter data map(alloc: flag_top )
    !$omp target enter data map(alloc: cltemp   )
    !$omp target enter data map(alloc: cmtemp   )
    !$omp target enter data map(alloc: chtemp   )
    !$omp target enter data map(alloc: cttemp   )
    !$omp target enter data map(alloc: dd_crm   )
    !$omp target enter data map(alloc: mui_crm  )
    !$omp target enter data map(alloc: mdi_crm  )
    !$omp target enter data map(alloc: ustar    )
    !$omp target enter data map(alloc: bflx     )
    !$omp target enter data map(alloc: wnd      )
    !$omp target enter data map(alloc: qtot     )
    !$omp target enter data map(alloc: colprec  )
    !$omp target enter data map(alloc: colprecs )
    !$omp target enter data map(alloc: crm_clear_rh_cnt )
    !$omp target enter data map(alloc: crm_clear_rh )

    !$omp target enter data map(alloc: crm_rad_temperature )
    !$omp target enter data map(alloc: crm_rad_qv )
    !$omp target enter data map(alloc: crm_rad_qc )
    !$omp target enter data map(alloc: crm_rad_qi )
    !$omp target enter data map(alloc: crm_rad_cld )
    !$omp target enter data map(alloc: crm_rad_qrad )

    !$omp target enter data map(alloc: crm_state_u_wind )
    !$omp target enter data map(alloc: crm_state_v_wind )
    !$omp target enter data map(alloc: crm_state_w_wind )
    !$omp target enter data map(alloc: crm_state_temperature )
    !$omp target enter data map(alloc: crm_state_qt )
    !$omp target enter data map(alloc: crm_state_qp )
    !$omp target enter data map(alloc: crm_state_qn )

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
    !$omp target enter data map(alloc: crm_output_tkew)
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
    !$omp target enter data map(alloc: crm_output_subcycle_factor)

  end subroutine allocate_crm

  subroutine deallocate_crm()
    implicit none
    !$omp target exit data map(delete: t00      )
    !$omp target exit data map(delete: tln      )
    !$omp target exit data map(delete: qln      )
    !$omp target exit data map(delete: qccln    )
    !$omp target exit data map(delete: qiiln    )
    !$omp target exit data map(delete: uln      )
    !$omp target exit data map(delete: vln      )
    !$omp target exit data map(delete: cwp      )
    !$omp target exit data map(delete: cwph     )
    !$omp target exit data map(delete: cwpm     )
    !$omp target exit data map(delete: cwpl     )
    !$omp target exit data map(delete: flag_top )
    !$omp target exit data map(delete: cltemp   )
    !$omp target exit data map(delete: cmtemp   )
    !$omp target exit data map(delete: chtemp   )
    !$omp target exit data map(delete: cttemp   )
    !$omp target exit data map(delete: dd_crm   )
    !$omp target exit data map(delete: mui_crm  )
    !$omp target exit data map(delete: mdi_crm  )
    !$omp target exit data map(delete: ustar    )
    !$omp target exit data map(delete: bflx     )
    !$omp target exit data map(delete: wnd      )
    !$omp target exit data map(delete: qtot     )
    !$omp target exit data map(delete: colprec  )
    !$omp target exit data map(delete: colprecs )
    !$omp target exit data map(delete: crm_clear_rh_cnt )
    !$omp target exit data map(delete: crm_clear_rh )

    !$omp target exit data map(delete: crm_rad_temperature )
    !$omp target exit data map(delete: crm_rad_qv )
    !$omp target exit data map(delete: crm_rad_qc )
    !$omp target exit data map(delete: crm_rad_qi )
    !$omp target exit data map(delete: crm_rad_cld )
    !$omp target exit data map(delete: crm_rad_qrad )

    !$omp target exit data map(delete: crm_state_u_wind )
    !$omp target exit data map(delete: crm_state_v_wind )
    !$omp target exit data map(delete: crm_state_w_wind )
    !$omp target exit data map(delete: crm_state_temperature )
    !$omp target exit data map(delete: crm_state_qt )
    !$omp target exit data map(delete: crm_state_qp )
    !$omp target exit data map(delete: crm_state_qn )

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
    !$omp target exit data map(delete: crm_output_tkew          )
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
    !$omp target exit data map(delete: crm_output_subcycle_factor )

    deallocate( t00)
    deallocate( tln)
    deallocate( qln)
    deallocate( qccln)
    deallocate( qiiln)
    deallocate( uln)
    deallocate( vln)
#if defined(MMF_ESMT)
    deallocate( uln_esmt)
    deallocate( vln_esmt)
#endif
    deallocate( cwp)
    deallocate( cwph)
    deallocate( cwpm)
    deallocate( cwpl)
    deallocate( flag_top)
    deallocate( cltemp)
    deallocate( cmtemp)
    deallocate( chtemp)
    deallocate( cttemp)
    deallocate( dd_crm  )
    deallocate( mui_crm )
    deallocate( mdi_crm )
    deallocate( ustar )
    deallocate( bflx )
    deallocate( wnd )
    deallocate( qtot )
    deallocate( colprec  )
    deallocate( colprecs )
    deallocate( crm_clear_rh_cnt )
  end subroutine deallocate_crm

  subroutine update_device_crm()
    implicit none
    !$omp target update to( t00      )
    !$omp target update to( tln      )
    !$omp target update to( qln      )
    !$omp target update to( qccln    )
    !$omp target update to( qiiln    )
    !$omp target update to( uln      )
    !$omp target update to( vln      )
    !$omp target update to( cwp      )
    !$omp target update to( cwph     )
    !$omp target update to( cwpm     )
    !$omp target update to( cwpl     )
    !$omp target update to( flag_top )
    !$omp target update to( cltemp   )
    !$omp target update to( cmtemp   )
    !$omp target update to( chtemp   )
    !$omp target update to( cttemp   )
    !$omp target update to( dd_crm   )
    !$omp target update to( mui_crm  )
    !$omp target update to( mdi_crm  )
    !$omp target update to( ustar    )
    !$omp target update to( bflx     )
    !$omp target update to( wnd      )
    !$omp target update to( qtot     )
    !$omp target update to( colprec  )
    !$omp target update to( colprecs )
    !$omp target update to( crm_clear_rh )
    !$omp target update to( crm_clear_rh_cnt )

    !$omp target update to( crm_rad_temperature )
    !$omp target update to( crm_rad_qv )
    !$omp target update to( crm_rad_qc )
    !$omp target update to( crm_rad_qi )
    !$omp target update to( crm_rad_cld )
    !$omp target update to( crm_rad_qrad )

    !$omp target update to( crm_state_u_wind )
    !$omp target update to( crm_state_v_wind )
    !$omp target update to( crm_state_w_wind )
    !$omp target update to( crm_state_temperature )
    !$omp target update to( crm_state_qt )
    !$omp target update to( crm_state_qp )
    !$omp target update to( crm_state_qn )

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

  end subroutine update_device_crm
  !\
  ! update host from device
  !/
  subroutine update_host_crm()
    implicit none
    !$omp target update from( t00      )
    !$omp target update from( tln      )
    !$omp target update from( qln      )
    !$omp target update from( qccln    )
    !$omp target update from( qiiln    )
    !$omp target update from( uln      )
    !$omp target update from( vln      )
    !$omp target update from( cwp      )
    !$omp target update from( cwph     )
    !$omp target update from( cwpm     )
    !$omp target update from( cwpl     )
    !$omp target update from( flag_top )
    !$omp target update from( cltemp   )
    !$omp target update from( cmtemp   )
    !$omp target update from( chtemp   )
    !$omp target update from( cttemp   )
    !$omp target update from( dd_crm   )
    !$omp target update from( mui_crm  )
    !$omp target update from( mdi_crm  )
    !$omp target update from( ustar    )
    !$omp target update from( bflx     )
    !$omp target update from( wnd      )
    !$omp target update from( qtot     )
    !$omp target update from( colprec  )
    !$omp target update from( colprecs )
    !$omp target update from( crm_clear_rh )
    !$omp target update from( crm_clear_rh_cnt )

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
    !$omp target update from (crm_output_sltend)
    !$omp target update from (crm_output_qltend)
    !$omp target update from (crm_output_qcltend)
    !$omp target update from (crm_output_qiltend)
    !$omp target update from (crm_output_cld)
    !$omp target update from (crm_output_gicewp)
    !$omp target update from (crm_output_gliqwp)
    !$omp target update from (crm_output_mctot)
    !$omp target update from (crm_output_mcup)
    !$omp target update from (crm_output_mcdn)

    !$omp target update from (crm_output_mcuup)
    !$omp target update from (crm_output_mcudn)
    !$omp target update from (crm_output_mu_crm)
    !$omp target update from (crm_output_md_crm)
    !$omp target update from (crm_output_du_crm)
    !$omp target update from (crm_output_eu_crm)
    !$omp target update from (crm_output_ed_crm)
    !$omp target update from (crm_output_jt_crm)
    !$omp target update from (crm_output_mx_crm)
    !$omp target update from (crm_output_flux_qt)
    !$omp target update from (crm_output_fluxsgs_qt)
    !$omp target update from (crm_output_tkez)
    !$omp target update from (crm_output_tkew)
    !$omp target update from (crm_output_tkesgsz)
    !$omp target update from (crm_output_tkz)
    !$omp target update from (crm_output_tkew)
    !$omp target update from (crm_output_flux_u)
    !$omp target update from (crm_output_flux_v)
    !$omp target update from (crm_output_flux_qp)
    !$omp target update from (crm_output_precflux)
    !$omp target update from (crm_output_qt_ls)
    !$omp target update from (crm_output_qt_trans)
    !$omp target update from (crm_output_qp_trans) 
    !$omp target update from (crm_output_qp_fall)
    !$omp target update from (crm_output_qp_src)
    !$omp target update from (crm_output_qp_evp)
    !$omp target update from (crm_output_t_ls)
    !$omp target update from (crm_output_prectend)
    !$omp target update from (crm_output_precstend)
    !$omp target update from (crm_output_taux)
    !$omp target update from (crm_output_tauy)
    !$omp target update from (crm_output_z0m)
    !$omp target update from (crm_output_subcycle_factor)

    !$omp target update from( crm_rad_temperature )
    !$omp target update from( crm_rad_qv )
    !$omp target update from( crm_rad_qc )
    !$omp target update from( crm_rad_qi )
    !$omp target update from( crm_rad_cld )
    !$omp target update from( crm_rad_qrad )

    !$omp target update from( crm_state_u_wind )
    !$omp target update from( crm_state_v_wind )
    !$omp target update from( crm_state_w_wind )
    !$omp target update from( crm_state_temperature )
    !$omp target update from( crm_state_qt )
    !$omp target update from( crm_state_qp )
    !$omp target update from( crm_state_qn )
  end subroutine update_host_crm
end subroutine crm
end module crm_module
