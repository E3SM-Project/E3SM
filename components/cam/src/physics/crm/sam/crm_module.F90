
module crm_module
  use openacc_utils, only: prefetch
  use perf_mod
  use task_init_mod, only: task_init
  use abcoefs_mod, only: abcoefs
  use kurant_mod, only: kurant
  use setperturb_mod, only: setperturb
  use boundaries_mod, only: boundaries
  use forcing_mod, only: forcing
  use advect_mom_mod, only: advect_mom
  use adams_mod, only: adams
  use advect_all_scalars_mod, only: advect_all_scalars
  use sat_mod
  use crmsurface_mod
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

  use crm_state_module,       only: crm_state_type
  use crm_rad_module,         only: crm_rad_type
  use crm_input_module,       only: crm_input_type
  use crm_output_module,      only: crm_output_type
  use crm_ecpp_output_module, only: crm_ecpp_output_type
  use phys_grid             , only: get_rlon_p, get_rlat_p, get_gcol_p  
#ifdef ECPP  
  use module_ecpp_crm_driver, only: ecpp_crm_stat, ecpp_crm_init, ecpp_crm_cleanup
#endif
!---------------------------------------------------------------
!  Super-parameterization's main driver
!  Marat Khairoutdinov, 2001-2009
!---------------------------------------------------------------
use setparm_mod, only : setparm

contains

subroutine crm(lchnk, icol, ncrms, dt_gl, plev, &
                crm_input, crm_state, crm_rad,  &
#ifdef CLUBB_CRM
                clubb_buffer,           &
                crm_cld, clubb_tk,      &
                clubb_tkh, relvar,      &
                accre_enhan, qclvar,    &
#endif
#ifdef MAML
                crm_pcp,     crm_snw,              &
#endif
                crm_ecpp_output, crm_output )
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    use shr_kind_mod          , only: r8 => shr_kind_r8
    use ppgrid                , only: pcols
    use vars
    use params
    use microphysics
    use sgs
    use crmtracers
    use scalar_momentum_mod
#ifdef MODAL_AERO
    use modal_aero_data       , only: ntot_amode
#endif
    use crmdims               , only: nclubbvars, crm_nx_rad, crm_ny_rad
#ifdef CLUBB_CRM
    use clubb_sgs             , only: advance_clubb_sgs, clubb_sgs_setup, clubb_sgs_cleanup, apply_clubb_sgs_tndcy, apply_clubb_sgs_tndcy_scalars, &
                                      apply_clubb_sgs_tndcy_mom, t2thetal, total_energy
    use clubb_precision       , only: time_precision, core_rknd
    use clubbvars             , only: up2, vp2, wprtp, wpthlp, wp2, wp3, rtp2, thlp2, rtpthlp, upwp, vpwp, cloud_frac, t_tndcy, qc_tndcy, qv_tndcy, &
                                      u_tndcy, v_tndcy, lrestart_clubb, rho_ds_zt, rho_ds_zm, thv_ds_zt, thv_ds_zm, invrs_rho_ds_zt, invrs_rho_ds_zm, &
                                      tracer_tndcy, sclrp2, sclrprtp, sclrpthlp, wpsclrp, relvarg, accre_enhang, qclvarg, edsclr_dim, sclr_dim, rho_ds_zt, &
                                      rho_ds_zm, rtm_spurious_source, thlm_spurious_source
    use fill_holes            , only: vertical_integral
    use numerical_check       , only: calculate_spurious_source
    use grid_class            , only: gr
#endif
#ifdef ECPP
    use ecppvars              , only: qlsink, precr, precsolid, &
                                      area_bnd_final, area_bnd_sum, area_cen_final, area_cen_sum, &
                                      mass_bnd_final, mass_bnd_sum, rh_cen_sum, qcloud_cen_sum, qice_cen_sum, &
                                      qlsink_cen_sum, precr_cen_sum, precsolid_cen_sum, xkhvsum, wup_thresh, wdown_thresh, &
                                      wwqui_cen_sum, wwqui_bnd_sum, wwqui_cloudy_cen_sum, wwqui_cloudy_bnd_sum, &
                                      qlsink_bf_cen_sum, qlsink_avg_cen_sum, prain_cen_sum, qlsink_bf, prain
    use ecppvars              , only: NCLASS_CL, ncls_ecpp_in, NCLASS_PR
#endif /* ECPP */
    use accelerate_crm_mod    , only: use_crm_accel, crm_accel_factor, crm_accel_nstop, accelerate_crm
    use cam_abortutils        , only: endrun
    use time_manager          , only: get_nstep

    implicit none

    !-----------------------------------------------------------------------------------------------
    ! Interface variable declarations
    !-----------------------------------------------------------------------------------------------

    integer , intent(in   ) :: lchnk                            ! chunk identifier (only for lat/lon and random seed)
    integer , intent(in   ) :: ncrms                            ! Number of "vector" GCM columns to push down into CRM for SIMD vectorization / more threading
    integer , intent(in   ) :: plev                             ! number of levels in parent model
    real(r8), intent(in   ) :: dt_gl                            ! global model's time step
    integer , intent(in   ) :: icol                (ncrms)      ! column identifier (only for lat/lon and random seed)
    type(crm_input_type),      intent(in   ) :: crm_input
    type(crm_state_type),      intent(inout) :: crm_state
    type(crm_rad_type), target,intent(inout) :: crm_rad
#ifdef CLUBB_CRM
    real(r8), intent(inout), target :: clubb_buffer(ncrms,crm_nx, crm_ny, crm_nz+1,1:nclubbvars)
    real(r8), intent(  out) :: crm_cld             (ncrms,crm_nx, crm_ny, crm_nz+1)
    real(r8), intent(  out) :: clubb_tk            (ncrms,crm_nx, crm_ny, crm_nz)
    real(r8), intent(  out) :: clubb_tkh           (ncrms,crm_nx, crm_ny, crm_nz)
    real(r8), intent(  out) :: relvar              (ncrms,crm_nx, crm_ny, crm_nz)
    real(r8), intent(  out) :: accre_enhan         (ncrms,crm_nx, crm_ny, crm_nz)
    real(r8), intent(  out) :: qclvar              (ncrms,crm_nx, crm_ny, crm_nz)
#endif /* CLUBB_CRM */
#ifdef MAML
    real(r8), intent(inout) :: crm_pcp(ncrms, crm_nx,crm_ny)  ! CRM precip rate (m/s)
    real(r8), intent(inout) :: crm_snw(ncrms,crm_nx,crm_ny) ! CRM snow rate (m/s)
#endif
    type(crm_ecpp_output_type),intent(inout) :: crm_ecpp_output
    type(crm_output_type), target,     intent(inout) :: crm_output

    !-----------------------------------------------------------------------------------------------
    ! Local variable declarations
    !-----------------------------------------------------------------------------------------------

    real(r8),       parameter :: umax = 0.5*crm_dx/crm_dt       ! maxumum ampitude of the l.s. wind
    real(r8),       parameter :: wmin = 2.                      ! minimum up/downdraft velocity for stat
    real(crm_rknd), parameter :: cwp_threshold = 0.001          ! threshold for cloud condensate for shaded fraction calculation
    integer,        parameter :: perturb_seed_scale = 1000      ! scaling value for setperturb() seed value (seed = gcol * perturb_seed_scale)
    real(r8)        :: crm_run_time                             ! length of CRM integration
    real(r8)        :: icrm_run_time                            ! = 1 / crm_run_time
    real(r8)        :: factor_xy, factor_xyt, idt_gl
    real(crm_rknd)  :: tmp1, tmp2, tmp
    real(crm_rknd)  :: u2z,v2z,w2z
    integer         :: i,j,k,l,ptop,nn,icyc,icrm
    integer         :: kx
    real(crm_rknd)  :: qsat, omg
    real(crm_rknd), allocatable  :: colprec(:), colprecs(:)
    real(crm_rknd), allocatable  :: ustar(:), bflx(:), wnd(:)
    real(r8)      , allocatable  :: qtot (:,:)    ! Total water for water conservation check

    !!! These should all be inputs
    integer         :: igstep            ! GCM time steps
    integer         :: iseed             ! seed for random perturbation
    !!! variables for radiation grouping method
    real(crm_rknd) :: crm_nx_rad_fac
    real(crm_rknd) :: crm_ny_rad_fac
    integer        :: i_rad
    integer        :: j_rad
    logical :: crm_accel_ceaseflag   ! indicates if accelerate_crm needs to be aborted for remainder of crm call

    !!! Arrays
    real(crm_rknd), allocatable :: t00(:,:)
    real(crm_rknd), allocatable :: tln  (:,:)
    real(crm_rknd), allocatable :: qln  (:,:)
    real(crm_rknd), allocatable :: qccln(:,:)
    real(crm_rknd), allocatable :: qiiln(:,:)
    real(crm_rknd), allocatable :: uln  (:,:)
    real(crm_rknd), allocatable :: vln  (:,:)
#if defined(SP_ESMT)
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

    real(crm_rknd), pointer :: crm_rad_temperature  (:,:,:,:)
    real(crm_rknd), pointer :: crm_rad_qv           (:,:,:,:)
    real(crm_rknd), pointer :: crm_rad_qc           (:,:,:,:)
    real(crm_rknd), pointer :: crm_rad_qi           (:,:,:,:)
    real(crm_rknd), pointer :: crm_rad_cld          (:,:,:,:)
    real(crm_rknd), pointer :: crm_rad_qrad         (:,:,:,:)
    real(crm_rknd), pointer :: crm_state_u_wind     (:,:,:,:)
    real(crm_rknd), pointer :: crm_state_v_wind     (:,:,:,:)
    real(crm_rknd), pointer :: crm_state_w_wind     (:,:,:,:)
    real(crm_rknd), pointer :: crm_state_temperature(:,:,:,:)
    real(crm_rknd), pointer :: crm_state_qt         (:,:,:,:)
    real(crm_rknd), pointer :: crm_state_qp         (:,:,:,:)
    real(crm_rknd), pointer :: crm_state_qn         (:,:,:,:)

  !-----------------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------------

  allocate( t00(ncrms,nz) )
  allocate( tln(ncrms,plev) )
  allocate( qln(ncrms,plev) )
  allocate( qccln(ncrms,plev) )
  allocate( qiiln(ncrms,plev) )
  allocate( uln(ncrms,plev) )
  allocate( vln(ncrms,plev) )
#if defined(SP_ESMT)
  allocate( uln_esmt(plev,ncrms) )
  allocate( vln_esmt(plev,ncrms) )
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

  call prefetch( t00      )
  call prefetch( tln      )
  call prefetch( qln      )
  call prefetch( qccln    )
  call prefetch( qiiln    )
  call prefetch( uln      )
  call prefetch( vln      )
  call prefetch( cwp      ) 
  call prefetch( cwph     ) 
  call prefetch( cwpm     ) 
  call prefetch( cwpl     ) 
  call prefetch( flag_top ) 
  call prefetch( cltemp   ) 
  call prefetch( cmtemp   ) 
  call prefetch( chtemp   ) 
  call prefetch( cttemp   ) 
  call prefetch( dd_crm   ) 
  call prefetch( mui_crm  ) 
  call prefetch( mdi_crm  ) 
  call prefetch( ustar    ) 
  call prefetch( bflx     ) 
  call prefetch( wnd      ) 
  call prefetch( qtot     ) 
  call prefetch( colprec  ) 
  call prefetch( colprecs ) 

  call allocate_params(ncrms)
  call allocate_vars(ncrms)
  call allocate_grid(ncrms)
  call allocate_tracers(ncrms)
  call allocate_sgs(ncrms)
  call allocate_micro(ncrms)
#ifdef sam1mom
  call allocate_micro_params(ncrms)
#endif
#if defined(SP_ESMT)
  call allocate_scalar_momentum(ncrms)
#endif

  crm_rad_temperature => crm_rad%temperature(1:ncrms,:,:,:)
  crm_rad_qv          => crm_rad%qv         (1:ncrms,:,:,:)
  crm_rad_qc          => crm_rad%qc         (1:ncrms,:,:,:)
  crm_rad_qi          => crm_rad%qi         (1:ncrms,:,:,:)
  crm_rad_cld         => crm_rad%cld        (1:ncrms,:,:,:)
  crm_rad_qrad        => crm_rad%qrad       (1:ncrms,:,:,:)

  crm_state_u_wind      => crm_state%u_wind     (1:ncrms,:,:,:)
  crm_state_v_wind      => crm_state%v_wind     (1:ncrms,:,:,:)
  crm_state_w_wind      => crm_state%w_wind     (1:ncrms,:,:,:)
  crm_state_temperature => crm_state%temperature(1:ncrms,:,:,:)
  crm_state_qt          => crm_state%qt         (1:ncrms,:,:,:)
  crm_state_qp          => crm_state%qp         (1:ncrms,:,:,:)
  crm_state_qn          => crm_state%qn         (1:ncrms,:,:,:)
  
  crm_accel_ceaseflag = .false.

  !Loop over "vector columns"
  do icrm = 1 , ncrms
    latitude0 (icrm) = get_rlat_p(lchnk, icol(icrm)) * 57.296_r8
    longitude0(icrm) = get_rlon_p(lchnk, icol(icrm)) * 57.296_r8
  enddo

  igstep = get_nstep()

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
  crm_rad%nc = 0.0
  crm_rad%ni = 0.0
  crm_rad%qs = 0.0
  crm_rad%ns = 0.0
#endif /* m2005 */
  do icrm = 1 , ncrms
    bflx(icrm) = crm_input%bflxls(icrm)
    wnd (icrm) = crm_input%wndls (icrm)
  enddo

!-----------------------------------------

  call task_init ()
  call setparm()

  do icrm = 1 , ncrms
    fcor(icrm)= 4*pi/86400.*sin(latitude0(icrm)*pi/180.)
    fcorz(icrm) = sqrt(4.*(2*pi/(3600.*24.))**2-fcor(icrm)**2)
    fcory(icrm,:) = fcor(icrm)
    fcorzy(icrm,:) = fcorz(icrm)
    do j=1,ny
      do i=1,nx
        latitude(icrm,i,j) = latitude0(icrm)
        longitude(icrm,i,j) = longitude0(icrm)
      end do
    end do

    if(crm_input%ocnfrac(icrm).gt.0.5) then
       OCEAN(icrm) = .true.
    else
       LAND(icrm) = .true.
    end if

    ! Create CRM vertical grid and initialize some vertical reference arrays:
    do k = 1, nzm
      z(icrm,k) = crm_input%zmid(icrm,plev-k+1) - crm_input%zint(icrm,plev+1)
      zi(icrm,k) = crm_input%zint(icrm,plev-k+2)- crm_input%zint(icrm,plev+1)
      pres(icrm,k) = crm_input%pmid(icrm,plev-k+1)/100.
      presi(icrm,k) = crm_input%pint(icrm,plev-k+2)/100.
      prespot(icrm,k)=(1000./pres(icrm,k))**(rgas/cp)
      bet(icrm,k) = ggr/crm_input%tl(icrm,plev-k+1)
      gamaz(icrm,k)=ggr/cp*z(icrm,k)
    end do ! k
   ! zi(icrm,nz) =  crm_input%zint(plev-nz+2)
    zi(icrm,nz) = crm_input%zint(icrm,plev-nz+2)-crm_input%zint(icrm,plev+1) !+++mhwang, 2012-02-04
    presi(icrm,nz) = crm_input%pint(icrm, plev-nz+2)/100.

    dz(icrm) = 0.5*(z(icrm,1)+z(icrm,2))
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
      rho(icrm,k) = crm_input%pdel(icrm,plev-k+1)/ggr/(adz(icrm,k)*dz(icrm))
    end do
    do k=2,nzm
    ! rhow(icrm,k) = 0.5*(rho(icrm,k)+rho(icrm,k-1))
    !+++mhwang fix the rhow bug (rhow needes to be consistent with crm_input%pmid)
    !2012-02-04 Minghuai Wang (minghuai.wang@pnnl.gov)
      rhow(icrm,k) = (crm_input%pmid(icrm,plev-k+2)-crm_input%pmid(icrm,plev-k+1))/ggr/(adzw(icrm,k)*dz(icrm))
    end do
    rhow(icrm,1) = 2.*rhow(icrm,2) - rhow(icrm,3)
    rhow(icrm,nz)= 2.*rhow(icrm,nzm) - rhow(icrm,nzm-1)
  enddo

  call t_startf('crm_gpu_region')

  !  Initialize CRM fields:
  !$acc parallel loop collapse(4) async(asyncid)
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
        enddo
      enddo
    enddo
  enddo

  ! limit the velocity at the very first step:
  if(u(1,1,1,1).eq.u(1,2,1,1).and.u(1,3,1,2).eq.u(1,4,1,2)) then
    !$acc parallel loop collapse(4) async(asyncid)
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

#if defined(SP_ESMT)
  do k=1,nzm
    u_esmt(:,:,:,k) = crm_input%ul_esmt(icrm,plev-k+1)
    v_esmt(:,:,:,k) = crm_input%vl_esmt(icrm,plev-k+1)
  end do
#endif

  ! Populate microphysics array from crm_state
  !$acc parallel loop collapse(4) async(asyncid)
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
      naer(icrm,k, 1:ntot_amode) = crm_input%naermod (icrm,l, 1:ntot_amode)
      vaer(icrm,k, 1:ntot_amode) = crm_input%vaerosol(icrm,l, 1:ntot_amode)
      hgaer(icrm,k, 1:ntot_amode) = crm_input%hygro   (icrm,l, 1:ntot_amode)
#endif /* MODAL_AERO */
      do j=1, ny
        do i=1, nx
          if(cloudliq(icrm,i,j,k).gt.0) then
            if(dopredictNc) then
              if( micro_field(icrm,i,j,k,incl).eq.0) micro_field(icrm,i,j,k,incl) = 1.0e6*Nc0/rho(icrm,k)
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

  !$acc parallel loop async(asyncid)
  do icrm = 1 , ncrms
    colprec (icrm)=0
    colprecs(icrm)=0
  enddo
  !$acc parallel loop collapse(2) async(asyncid)
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

  !$acc parallel loop collapse(4) async(asyncid)
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          t(icrm,i,j,k) = tabs(icrm,i,j,k)+gamaz(icrm,k)-fac_cond*qcl(icrm,i,j,k)-fac_sub*qci(icrm,i,j,k) &
                                                        -fac_cond*qpl(icrm,i,j,k)-fac_sub*qpi(icrm,i,j,k)

          tmp = (qpl(icrm,i,j,k)+qpi(icrm,i,j,k))*crm_input%pdel(icrm,plev-k+1)
          !$acc atomic update
          colprec (icrm)=colprec (icrm)+tmp

          tmp = qpi(icrm,i,j,k)*crm_input%pdel(icrm,plev-k+1)
          !$acc atomic update
          colprecs(icrm)=colprecs(icrm)+tmp
          !$acc atomic update
          u0   (icrm,k)=u0   (icrm,k)+u(icrm,i,j,k)
          !$acc atomic update
          v0   (icrm,k)=v0   (icrm,k)+v(icrm,i,j,k)
          !$acc atomic update
          t0   (icrm,k)=t0   (icrm,k)+t(icrm,i,j,k)

          tmp = t(icrm,i,j,k)+fac_cond*qpl(icrm,i,j,k)+fac_sub*qpi(icrm,i,j,k)
          !$acc atomic update
          t00  (icrm,k)=t00  (icrm,k)+tmp
          !$acc atomic update
          tabs0(icrm,k)=tabs0(icrm,k)+tabs(icrm,i,j,k)

          tmp = qv(icrm,i,j,k)+qcl(icrm,i,j,k)+qci(icrm,i,j,k)
          !$acc atomic update
          q0   (icrm,k)=q0   (icrm,k)+tmp
          !$acc atomic update
          qv0  (icrm,k)=qv0  (icrm,k)+qv(icrm,i,j,k)

          tmp = qcl(icrm,i,j,k) + qci(icrm,i,j,k)
          !$acc atomic update
          qn0  (icrm,k)=qn0  (icrm,k)+tmp

          tmp = qpl(icrm,i,j,k) + qpi(icrm,i,j,k)
          !$acc atomic update
          qp0  (icrm,k)=qp0  (icrm,k)+tmp
          !$acc atomic update
          tke0 (icrm,k)=tke0 (icrm,k)+sgs_field(icrm,i,j,k,1)
        enddo
      enddo
    enddo
  enddo

  !$acc parallel loop collapse(2) async(asyncid)
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
      uln  (icrm,l) = min( umax, max(-umax,crm_input%ul(icrm,l)) )
#ifdef MAML
      !open the crm v component
      vln  (icrm,l) = min( umax, max(-umax,crm_input%vl(icrm,l)) )
#else
      vln  (icrm,l) = min( umax, max(-umax,crm_input%vl(icrm,l)) )*YES3D
#endif
      ttend(icrm,k) = (crm_input%tl(icrm,l)+gamaz(icrm,k)- fac_cond*(crm_input%qccl(icrm,l)+crm_input%qiil(icrm,l))-fac_fus*crm_input%qiil(icrm,l)-t00(icrm,k))*idt_gl
      qtend(icrm,k) = (crm_input%ql(icrm,l)+crm_input%qccl(icrm,l)+crm_input%qiil(icrm,l)-q0(icrm,k))*idt_gl
      utend(icrm,k) = (uln(icrm,l)-u0(icrm,k))*idt_gl
      vtend(icrm,k) = (vln(icrm,l)-v0(icrm,k))*idt_gl
      ug0  (icrm,k) = uln(icrm,l)
      vg0  (icrm,k) = vln(icrm,l)
      tg0  (icrm,k) = crm_input%tl(icrm,l)+gamaz(icrm,k)-fac_cond*crm_input%qccl(icrm,l)-fac_sub*crm_input%qiil(icrm,l)
      qg0  (icrm,k) = crm_input%ql(icrm,l)+crm_input%qccl(icrm,l)+crm_input%qiil(icrm,l)
    end do ! k
  end do ! icrm

  !$acc parallel loop async(asyncid)
  do icrm = 1 , ncrms
    uhl(icrm) = u0(icrm,1)
    vhl(icrm) = v0(icrm,1)
    ! estimate roughness length assuming logarithmic profile of velocity near the surface:
    ustar(icrm) = sqrt(crm_input%tau00(icrm)/rho(icrm,1))
    z0(icrm) = z0_est(z(icrm,1),bflx(icrm),wnd(icrm),ustar(icrm))
    z0(icrm) = max(real(0.00001,crm_rknd),min(real(1.,crm_rknd),z0(icrm)))
    crm_output%timing_factor(icrm) = 0.
    crm_output%prectend (icrm)=colprec (icrm)
    crm_output%precstend(icrm)=colprecs(icrm)
  enddo

!---------------------------------------------------
#ifdef m2005
  crm_output%nc_mean = 0.
  crm_output%ni_mean = 0.
  crm_output%ns_mean = 0.
  crm_output%ng_mean = 0.
  crm_output%nr_mean = 0.
  crm_output%aut_a  = 0.
  crm_output%acc_a  = 0.
  crm_output%evpc_a = 0.
  crm_output%evpr_a = 0.
  crm_output%mlt_a  = 0.
  crm_output%sub_a  = 0.
  crm_output%dep_a  = 0.
  crm_output%con_a  = 0.
  aut1a  = 0.
  acc1a  = 0.
  evpc1a = 0.
  evpr1a = 0.
  mlt1a  = 0.
  sub1a  = 0.
  dep1a  = 0.
  con1a  = 0.
#endif /* m2005 */
  !$acc parallel loop collapse(2) async(asyncid)
  do k = 1 , plev+1
    do icrm = 1 , ncrms
      if (k <= plev) then
        crm_output%cld       (icrm,k) = 0.
        crm_output%cldtop    (icrm,k) = 0.
        crm_output%gicewp    (icrm,k) = 0
        crm_output%gliqwp    (icrm,k) = 0
        crm_output%mctot     (icrm,k) = 0.
        crm_output%mcup      (icrm,k) = 0.
        crm_output%mcdn      (icrm,k) = 0.
        crm_output%mcuup     (icrm,k) = 0.
        crm_output%mcudn     (icrm,k) = 0.
        crm_output%qc_mean   (icrm,k) = 0.
        crm_output%qi_mean   (icrm,k) = 0.
        crm_output%qs_mean   (icrm,k) = 0.
        crm_output%qg_mean   (icrm,k) = 0.
        crm_output%qr_mean   (icrm,k) = 0.
        crm_output%mu_crm    (icrm,k) = 0.
        crm_output%md_crm    (icrm,k) = 0.
        crm_output%eu_crm    (icrm,k) = 0.
        crm_output%du_crm    (icrm,k) = 0.
        crm_output%ed_crm    (icrm,k) = 0.
        crm_output%flux_qt   (icrm,k) = 0.
        crm_output%flux_u    (icrm,k) = 0.
        crm_output%flux_v    (icrm,k) = 0.
        crm_output%fluxsgs_qt(icrm,k) = 0.
        crm_output%tkez      (icrm,k) = 0.
        crm_output%tkesgsz   (icrm,k) = 0.
        crm_output%tkz       (icrm,k) = 0.
        crm_output%flux_qp   (icrm,k) = 0.
        crm_output%precflux  (icrm,k) = 0.
        crm_output%qt_trans  (icrm,k) = 0.
        crm_output%qp_trans  (icrm,k) = 0.
        crm_output%qp_fall   (icrm,k) = 0.
        crm_output%qp_evp    (icrm,k) = 0.
        crm_output%qp_src    (icrm,k) = 0.
        crm_output%qt_ls     (icrm,k) = 0.
        crm_output%t_ls      (icrm,k) = 0.
        dd_crm               (icrm,k) = 0.
      endif
      mui_crm(icrm,k) = 0.
      mdi_crm(icrm,k) = 0.
    enddo
  enddo
  !$acc parallel loop async(asyncid)
  do icrm = 1 , ncrms
    crm_output%jt_crm(icrm) = 0.
    crm_output%mx_crm(icrm) = 0.
  enddo

!--------------------------------------------------
#ifdef sam1mom
  if(doprecip) call precip_init(ncrms)
#endif
  !$acc wait(asyncid)

  do icrm = 1 , ncrms
    if ( igstep <= 1 ) then
        iseed = get_gcol_p(lchnk,icol(icrm)) * perturb_seed_scale
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
      call endrun('crm main')
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

    !$acc parallel loop async(asyncid)
    do icrm = 1 , ncrms
      crm_output%timing_factor(icrm) = crm_output%timing_factor(icrm)+1
    enddo

    !------------------------------------------------------------------
    !  Check if the dynamical time step should be decreased
    !  to handle the cases when the flow being locally linearly unstable
    !------------------------------------------------------------------
    call kurant(ncrms)
    !$acc wait(asyncid)

    do icyc=1,ncycle
      icycle = icyc
      dtn = dt/ncycle
      dt3(na) = dtn
      dtfactor = dtn/dt

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

      !!! Apply radiative tendency
      !$acc parallel loop collapse(4) async(asyncid)
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
#ifdef CLUBB_CRM
      if ( docloud .or. doclubb ) then
        call ice_fall(ncrms)
      endif
#else
      if(docloud) then
        call ice_fall(ncrms)
      endif
#endif

      !----------------------------------------------------------
      !     Update scalar boundaries after large-scale processes:
      call boundaries(ncrms,3)

      !---------------------------------------------------------
      !     Update boundaries for velocities:
      call boundaries(ncrms,0)

      !-----------------------------------------------
      !     surface fluxes:
      if (dosurface) call crmsurface(ncrms,bflx)

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
#if defined( SP_ESMT ) 
      call scalar_momentum_tend(ncrms)
#endif

      !-----------------------------------------------------------
      !       Cloud condensation/evaporation and precipitation processes:
#ifdef CLUBB_CRM
      if(docloud.or.dosmoke.or.doclubb) call micro_proc(ncrms)
#else
      if(docloud.or.dosmoke) call micro_proc(ncrms)
#endif /*CLUBB_CRM*/

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
    !$acc parallel loop collapse(3) async(asyncid)
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

    !$acc parallel loop gang vector collapse(3) async(asyncid)
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          do k=1,nzm
            l = plev-k+1
            tmp1 = rho(icrm,nz-k)*adz(icrm,nz-k)*dz(icrm)*(qcl(icrm,i,j,nz-k)+qci(icrm,i,j,nz-k))
            cwp(icrm,i,j) = cwp(icrm,i,j)+tmp1
            cttemp(icrm,i,j) = max(cf3d(icrm,i,j,nz-k), cttemp(icrm,i,j))
            if(cwp(icrm,i,j).gt.cwp_threshold.and.flag_top(icrm,i,j)) then
                !$acc atomic update
                crm_output%cldtop(icrm,l) = crm_output%cldtop(icrm,l) + 1
                flag_top(icrm,i,j) = .false.
            endif
            if(pres(icrm,nz-k).ge.700.) then
                cwpl(icrm,i,j) = cwpl(icrm,i,j)+tmp1
                cltemp(icrm,i,j) = max(cf3d(icrm,i,j,nz-k), cltemp(icrm,i,j))
            else if(pres(icrm,nz-k).lt.400.) then
                cwph(icrm,i,j) = cwph(icrm,i,j)+tmp1
                chtemp(icrm,i,j) = max(cf3d(icrm,i,j,nz-k), chtemp(icrm,i,j))
            else
                cwpm(icrm,i,j) = cwpm(icrm,i,j)+tmp1
                cmtemp(icrm,i,j) = max(cf3d(icrm,i,j,nz-k), cmtemp(icrm,i,j))
            endif
            tmp1 = rho(icrm,k)*adz(icrm,k)*dz(icrm)
            if(tmp1*(qcl(icrm,i,j,k)+qci(icrm,i,j,k)).gt.cwp_threshold) then
                 !$acc atomic update
                 crm_output%cld(icrm,l) = crm_output%cld(icrm,l) + cf3d(icrm,i,j,k)
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).gt.2*wmin) then
                   tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * cf3d(icrm,i,j,k)
                   !$acc atomic update
                   crm_output%mcup (icrm,l) = crm_output%mcup (icrm,l) + tmp
                   tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * (1.0 - cf3d(icrm,i,j,k))
                   !$acc atomic update
                   crm_output%mcuup(icrm,l) = crm_output%mcuup(icrm,l) + tmp
                 endif
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).lt.-2*wmin) then
                   tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * cf3d(icrm,i,j,k)
                   !$acc atomic update
                   crm_output%mcdn (icrm,l) = crm_output%mcdn (icrm,l) + tmp
                   tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k)) * (1. - cf3d(icrm,i,j,k))
                   !$acc atomic update
                   crm_output%mcudn(icrm,l) = crm_output%mcudn(icrm,l) + tmp
                 endif
            else
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).gt.2*wmin) then
                   tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k))
                   !$acc atomic update
                   crm_output%mcuup(icrm,l) = crm_output%mcuup(icrm,l) + tmp
                 endif
                 if(w(icrm,i,j,k+1)+w(icrm,i,j,k).lt.-2*wmin) then
                    tmp = rho(icrm,k)*0.5*(w(icrm,i,j,k+1)+w(icrm,i,j,k))
                   !$acc atomic update
                   crm_output%mcudn(icrm,l) = crm_output%mcudn(icrm,l) + tmp
                 endif
            endif

            !$acc atomic update
            crm_output%gliqwp(icrm,l) = crm_output%gliqwp(icrm,l) + qcl(icrm,i,j,k)
            !$acc atomic update
            crm_output%gicewp(icrm,l) = crm_output%gicewp(icrm,l) + qci(icrm,i,j,k)
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop gang vector collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            !!! Reduced radiation method allows for fewer radiation calculations
            !!! by collecting statistics and doing radiation over column groups
            i_rad = (i-1) / (nx/crm_nx_rad) + 1
            j_rad = (j-1) / (ny/crm_ny_rad) + 1

            !$acc atomic update
            crm_rad_temperature(icrm,i_rad,j_rad,k) = crm_rad_temperature(icrm,i_rad,j_rad,k) + tabs(icrm,i,j,k)
            tmp = max(real(0.,crm_rknd),qv(icrm,i,j,k))
            !$acc atomic update
            crm_rad_qv         (icrm,i_rad,j_rad,k) = crm_rad_qv         (icrm,i_rad,j_rad,k) + tmp
            !$acc atomic update
            crm_rad_qc         (icrm,i_rad,j_rad,k) = crm_rad_qc         (icrm,i_rad,j_rad,k) + qcl(icrm,i,j,k)
            !$acc atomic update
            crm_rad_qi         (icrm,i_rad,j_rad,k) = crm_rad_qi         (icrm,i_rad,j_rad,k) + qci(icrm,i,j,k)
            if (qcl(icrm,i,j,k) + qci(icrm,i,j,k) > 0) then
               !$acc atomic update
               crm_rad_cld     (icrm,i_rad,j_rad,k) = crm_rad_cld        (icrm,i_rad,j_rad,k) + cf3d(icrm,i,j,k)
            endif
#ifdef m2005
            !$acc atomic update
            crm_rad%nc         (icrm,i_rad,j_rad,k) = crm_rad%nc         (icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,incl)
            !$acc atomic update
            crm_rad%ni         (icrm,i_rad,j_rad,k) = crm_rad%ni         (icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,inci)
            !$acc atomic update
            crm_rad%qs         (icrm,i_rad,j_rad,k) = crm_rad%qs         (icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,iqs )
            !$acc atomic update
            crm_rad%ns         (icrm,i_rad,j_rad,k) = crm_rad%ns         (icrm,i_rad,j_rad,k) + micro_field(icrm,i,j,k,ins )
#endif
          enddo
        enddo
      enddo
    enddo

    ! Diagnose mass fluxes to drive CAM's convective transport of tracers.
    ! definition of mass fluxes is taken from Xu et al., 2002, QJRMS.
    !$acc parallel loop collapse(3) async(asyncid)
    do j=1, ny
      do i=1, nx
        do icrm = 1 , ncrms
          do k=1, nzm+1
            l=plev+1-k+1
            if(w(icrm,i,j,k).gt.0.) then
              kx=max(1, k-1)
              qsat = qsatw_crm(tabs(icrm,i,j,kx),pres(icrm,kx))
              if(qcl(icrm,i,j,kx)+qci(icrm,i,j,kx).gt.min(real(1.e-5,crm_rknd),0.01*qsat)) then
                tmp = rhow(icrm,k)*w(icrm,i,j,k)
                !$acc atomic update
                mui_crm(icrm,l) = mui_crm(icrm,l)+tmp
              endif
            else if (w(icrm,i,j,k).lt.0.) then
              kx=min(k+1, nzm)
              qsat = qsatw_crm(tabs(icrm,i,j,kx),pres(icrm,kx))
              if(qcl(icrm,i,j,kx)+qci(icrm,i,j,kx).gt.min(real(1.e-5,crm_rknd),0.01*qsat)) then
                tmp = rhow(icrm,k)*w(icrm,i,j,k)
                !$acc atomic update
                mdi_crm(icrm,l) = mdi_crm(icrm,l)+tmp
              else if(qpl(icrm,i,j,kx)+qpi(icrm,i,j,kx).gt.1.0e-4) then
                tmp = rhow(icrm,k)*w(icrm,i,j,k)
                !$acc atomic update
                mdi_crm(icrm,l) = mdi_crm(icrm,l)+tmp
              endif
            endif
          enddo
        enddo
      enddo
    enddo

    !$acc parallel loop collapse(3) async(asyncid)
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          if(cwp(icrm,i,j).gt.cwp_threshold) then
            !$acc atomic update
            crm_output%cltot(icrm) = crm_output%cltot(icrm) + cttemp(icrm,i,j)
          endif
          if(cwph(icrm,i,j).gt.cwp_threshold) then
            !$acc atomic update
            crm_output%clhgh(icrm) = crm_output%clhgh(icrm) + chtemp(icrm,i,j)
          endif
          if(cwpm(icrm,i,j).gt.cwp_threshold) then
            !$acc atomic update
            crm_output%clmed(icrm) = crm_output%clmed(icrm) + cmtemp(icrm,i,j)
          endif
          if(cwpl(icrm,i,j).gt.cwp_threshold) then
            !$acc atomic update
            crm_output%cllow(icrm) = crm_output%cllow(icrm) + cltemp(icrm,i,j)
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

  !$acc parallel loop collapse(4) async(asyncid)
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
          crm_rad%nc         (icrm,i,j,k) = crm_rad%nc         (icrm,i,j,k) * tmp1
          crm_rad%ni         (icrm,i,j,k) = crm_rad%ni         (icrm,i,j,k) * tmp1
          crm_rad%qs         (icrm,i,j,k) = crm_rad%qs         (icrm,i,j,k) * tmp1
          crm_rad%ns         (icrm,i,j,k) = crm_rad%ns         (icrm,i,j,k) * tmp1
#endif /* m2005 */
        enddo
      enddo
    enddo
  enddo

  ! no CRM tendencies above its top
  !$acc parallel loop collapse(2) async(asyncid)
  do k = 1 , ptop-1
    do icrm = 1 , ncrms
      tln  (icrm,k) = crm_input%tl  (icrm,k)
      qln  (icrm,k) = crm_input%ql  (icrm,k)
      qccln(icrm,k) = crm_input%qccl(icrm,k)
      qiiln(icrm,k) = crm_input%qiil(icrm,k)
      uln  (icrm,k) = crm_input%ul  (icrm,k)
      vln  (icrm,k) = crm_input%vl  (icrm,k)
    enddo
  enddo

  !  Compute tendencies due to CRM:
  !$acc parallel loop collapse(2) async(asyncid)
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
  !$acc parallel loop async(asyncid)
  do icrm = 1 , ncrms
    colprec (icrm)=0
    colprecs(icrm)=0
  enddo

#if defined( SP_ESMT )
  uln_esmt(1:ptop-1,:)  = crm_input%ul_esmt(:,1:ptop-1)
  vln_esmt(1:ptop-1,:)  = crm_input%vl_esmt(:,1:ptop-1)
  uln_esmt(ptop:plev,:) = 0.
  vln_esmt(ptop:plev,:) = 0.
#endif /* SP_ESMT */

  !$acc parallel loop collapse(4) async(asyncid)
  do k = 1,nzm
    do i=1,nx
      do j=1,ny
        do icrm=1,ncrms
          l = plev-k+1

          tmp = (qpl(icrm,i,j,k)+qpi(icrm,i,j,k))*crm_input%pdel(icrm,plev-k+1)
          !$acc atomic update
          colprec (icrm)= colprec (icrm)+tmp

          tmp = qpi(icrm,i,j,k)*crm_input%pdel(icrm,plev-k+1)
          !$acc atomic update
          colprecs(icrm)= colprecs(icrm)+tmp
          !$acc atomic update
          tln(icrm,l)  = tln(icrm,l)  +tabs(icrm,i,j,k)
          !$acc atomic update
          qln(icrm,l)  = qln(icrm,l)  +qv(icrm,i,j,k)
          !$acc atomic update
          qccln(icrm,l)= qccln(icrm,l)+qcl(icrm,i,j,k)
          !$acc atomic update
          qiiln(icrm,l)= qiiln(icrm,l)+qci(icrm,i,j,k)
          !$acc atomic update
          uln(icrm,l)  = uln(icrm,l)  +u(icrm,i,j,k)
          !$acc atomic update
          vln(icrm,l)  = vln(icrm,l)  +v(icrm,i,j,k)

#if defined(SP_ESMT)
          !$acc atomic update
          uln_esmt(l,icrm) = uln_esmt(l,icrm)+u_esmt(icrm,i,j,k)
          !$acc atomic update
          vln_esmt(l,icrm) = vln_esmt(l,icrm)+v_esmt(icrm,i,j,k)
#endif
        enddo ! j
      enddo ! i
    enddo ! k
  enddo ! icrm

#if defined(SP_ESMT)
  !$acc wait(asyncid)
  do icrm=1,ncrms
    uln_esmt(ptop:plev,icrm) = uln_esmt(ptop:plev,icrm) * factor_xy
    vln_esmt(ptop:plev,icrm) = vln_esmt(ptop:plev,icrm) * factor_xy

    crm_output%u_tend_esmt(icrm,:) = (uln_esmt(:,icrm) - crm_input%ul_esmt(icrm,:))*icrm_run_time
    crm_output%v_tend_esmt(icrm,:) = (vln_esmt(:,icrm) - crm_input%vl_esmt(icrm,:))*icrm_run_time

    ! don't use tendencies from two top levels,
    crm_output%u_tend_esmt(icrm,ptop:ptop+1) = 0.
    crm_output%v_tend_esmt(icrm,ptop:ptop+1) = 0.
  enddo
#endif

#if defined(SPMOMTRANS)
  !$acc wait(asyncid)
  do icrm=1,ncrms
    !!! resolved convective momentum transport (CMT) tendencies
    crm_output%ultend(icrm,:) = (uln(icrm,:) - crm_input%ul(icrm,:))*icrm_run_time
    crm_output%vltend(icrm,:) = (vln(icrm,:) - crm_input%vl(icrm,:))*icrm_run_time

    !!! don't use tendencies from two top levels
    crm_output%ultend(icrm,ptop:ptop+1) = 0.
    crm_output%vltend(icrm,ptop:ptop+1) = 0.
  enddo
#endif /* SPMOMTRANS */


  !$acc parallel loop collapse(2) async(asyncid)
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

  !$acc parallel loop collapse(2) async(asyncid)
  do k = 1 , plev
    do icrm = 1 , ncrms
      crm_output%sltend (icrm,k) = cp * (tln  (icrm,k) - crm_input%tl  (icrm,k)) * icrm_run_time
      crm_output%qltend (icrm,k) =      (qln  (icrm,k) - crm_input%ql  (icrm,k)) * icrm_run_time
      crm_output%qcltend(icrm,k) =      (qccln(icrm,k) - crm_input%qccl(icrm,k)) * icrm_run_time
      crm_output%qiltend(icrm,k) =      (qiiln(icrm,k) - crm_input%qiil(icrm,k)) * icrm_run_time
    enddo
  enddo
  !$acc parallel loop async(asyncid)
  do icrm = 1 , ncrms
    crm_output%prectend (icrm) = (colprec (icrm)-crm_output%prectend (icrm))/ggr*factor_xy * icrm_run_time
    crm_output%precstend(icrm) = (colprecs(icrm)-crm_output%precstend(icrm))/ggr*factor_xy * icrm_run_time
  enddo

  !!! don't use CRM tendencies from two crm top levels
  !!! radiation tendencies are added back after the CRM call (see crm_physics_tend)
  !$acc parallel loop collapse(2) async(asyncid)
  do k = ptop,ptop+1
    do icrm = 1 , ncrms
      crm_output%sltend (icrm,k) = 0.
      crm_output%qltend (icrm,k) = 0.
      crm_output%qcltend(icrm,k) = 0.
      crm_output%qiltend(icrm,k) = 0.
    enddo
  enddo

  !-------------------------------------------------------------
  !
  ! Save the last step to the permanent core:
  !$acc parallel loop collapse(4) async(asyncid)
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
          crm_output%tk (icrm,i,j,k) = sgs_field_diag(icrm,i,j,k,1)
          crm_output%tkh(icrm,i,j,k) = sgs_field_diag(icrm,i,j,k,2)
          crm_output%qcl (icrm,i,j,k) = qcl  (icrm,i,j,k)
          crm_output%qci (icrm,i,j,k) = qci  (icrm,i,j,k)
          crm_output%qpl (icrm,i,j,k) = qpl  (icrm,i,j,k)
          crm_output%qpi (icrm,i,j,k) = qpi  (icrm,i,j,k)
#ifdef m2005
          crm_output%wvar(icrm,i,j,k) = wvar (icrm,i,j,k)
          crm_output%aut (icrm,i,j,k) = aut1 (icrm,i,j,k)
          crm_output%acc (icrm,i,j,k) = acc1 (icrm,i,j,k)
          crm_output%evpc(icrm,i,j,k) = evpc1(icrm,i,j,k)
          crm_output%evpr(icrm,i,j,k) = evpr1(icrm,i,j,k)
          crm_output%mlt (icrm,i,j,k) = mlt1 (icrm,i,j,k)
          crm_output%sub (icrm,i,j,k) = sub1 (icrm,i,j,k)
          crm_output%dep (icrm,i,j,k) = dep1 (icrm,i,j,k)
          crm_output%con (icrm,i,j,k) = con1 (icrm,i,j,k)
#endif /* m2005 */
        enddo
      enddo
    enddo
  enddo
  !$acc parallel loop async(asyncid)
  do icrm = 1 , ncrms
    crm_output%z0m (icrm) = z0   (icrm)
    crm_output%taux(icrm) = taux0(icrm) / dble(nstop)
    crm_output%tauy(icrm) = tauy0(icrm) / dble(nstop)
  enddo

  !---------------------------------------------------------------
  !  Diagnostics:

  ! hm add 9/7/11, change from GCM-time step avg to end-of-timestep
  !$acc parallel loop collapse(4) async(asyncid)
  do k=1,nzm
    do j=1,ny
      do i=1,nx
        do icrm=1,ncrms
          l = plev-k+1
          !$acc atomic update
          crm_output%qc_mean(icrm,l) = crm_output%qc_mean(icrm,l) + qcl(icrm,i,j,k)
          !$acc atomic update
          crm_output%qi_mean(icrm,l) = crm_output%qi_mean(icrm,l) + qci(icrm,i,j,k)
          !$acc atomic update
          crm_output%qr_mean(icrm,l) = crm_output%qr_mean(icrm,l) + qpl(icrm,i,j,k)
#ifdef sam1mom
          omg = max(real(0.,crm_rknd),min(real(1.,crm_rknd),(tabs(icrm,i,j,k)-tgrmin)*a_gr))

          tmp = qpi(icrm,i,j,k)*omg
          !$acc atomic update
          crm_output%qg_mean(icrm,l) = crm_output%qg_mean(icrm,l) + tmp

          tmp = qpi(icrm,i,j,k)*(1.-omg)
          !$acc atomic update
          crm_output%qs_mean(icrm,l) = crm_output%qs_mean(icrm,l) + tmp
#else
          !$acc atomic update
          crm_output%qg_mean(icrm,l) = crm_output%qg_mean(icrm,l) + micro_field(icrm,i,j,k,iqg)
          !$acc atomic update
          crm_output%qs_mean(icrm,l) = crm_output%qs_mean(icrm,l) + micro_field(icrm,i,j,k,iqs)

          !$acc atomic update
          crm_output%nc_mean(icrm,l) = crm_output%nc_mean(icrm,l) + micro_field(icrm,i,j,k,incl)
          !$acc atomic update
          crm_output%ni_mean(icrm,l) = crm_output%ni_mean(icrm,l) + micro_field(icrm,i,j,k,inci)
          !$acc atomic update
          crm_output%nr_mean(icrm,l) = crm_output%nr_mean(icrm,l) + micro_field(icrm,i,j,k,inr )
          !$acc atomic update
          crm_output%ng_mean(icrm,l) = crm_output%ng_mean(icrm,l) + micro_field(icrm,i,j,k,ing )
          !$acc atomic update
          crm_output%ns_mean(icrm,l) = crm_output%ns_mean(icrm,l) + micro_field(icrm,i,j,k,ins )
#endif /* sam1mom */
        enddo
      enddo
    enddo
  enddo

  !$acc parallel loop collapse(2) async(asyncid)
  do k = 1 , plev
    do icrm = 1 , ncrms
      crm_output%cld   (icrm,k) = min( 1._r8, crm_output%cld   (icrm,k) * factor_xyt )
      crm_output%cldtop(icrm,k) = min( 1._r8, crm_output%cldtop(icrm,k) * factor_xyt )
      crm_output%gicewp(icrm,k) = crm_output%gicewp(icrm,k)*crm_input%pdel(icrm,k)*1000./ggr * factor_xyt
      crm_output%gliqwp(icrm,k) = crm_output%gliqwp(icrm,k)*crm_input%pdel(icrm,k)*1000./ggr * factor_xyt
      crm_output%mcup  (icrm,k) = crm_output%mcup (icrm,k) * factor_xyt
      crm_output%mcdn  (icrm,k) = crm_output%mcdn (icrm,k) * factor_xyt
      crm_output%mcuup (icrm,k) = crm_output%mcuup(icrm,k) * factor_xyt
      crm_output%mcudn (icrm,k) = crm_output%mcudn(icrm,k) * factor_xyt
      crm_output%mctot (icrm,k) = crm_output%mcup(icrm,k) + crm_output%mcdn(icrm,k) + crm_output%mcuup(icrm,k) + crm_output%mcudn(icrm,k)

      crm_output%qc_mean(icrm,k) = crm_output%qc_mean(icrm,k) * factor_xy
      crm_output%qi_mean(icrm,k) = crm_output%qi_mean(icrm,k) * factor_xy
      crm_output%qs_mean(icrm,k) = crm_output%qs_mean(icrm,k) * factor_xy
      crm_output%qg_mean(icrm,k) = crm_output%qg_mean(icrm,k) * factor_xy
      crm_output%qr_mean(icrm,k) = crm_output%qr_mean(icrm,k) * factor_xy
    enddo
  enddo

#ifdef m2005
  do icrm=1,ncrms
    crm_output%nc_mean(icrm,:) = crm_output%nc_mean(icrm,:) * factor_xy
    crm_output%ni_mean(icrm,:) = crm_output%ni_mean(icrm,:) * factor_xy
    crm_output%ns_mean(icrm,:) = crm_output%ns_mean(icrm,:) * factor_xy
    crm_output%ng_mean(icrm,:) = crm_output%ng_mean(icrm,:) * factor_xy
    crm_output%nr_mean(icrm,:) = crm_output%nr_mean(icrm,:) * factor_xy

    ! hm 8/31/11 new output, gcm-grid- and time-step avg
    ! add loop over i,j do get horizontal avg, and flip vertical array
    do k=1,nzm
      l = plev-k+1
      do j=1,ny
        do i=1,nx
          crm_output%aut_a (icrm,l) = crm_output%aut_a (icrm,l) + aut1a(icrm,i,j,k)
          crm_output%acc_a (icrm,l) = crm_output%acc_a (icrm,l) + acc1a(icrm,i,j,k)
          crm_output%evpc_a(icrm,l) = crm_output%evpc_a(icrm,l) + evpc1a(icrm,i,j,k)
          crm_output%evpr_a(icrm,l) = crm_output%evpr_a(icrm,l) + evpr1a(icrm,i,j,k)
          crm_output%mlt_a (icrm,l) = crm_output%mlt_a (icrm,l) + mlt1a(icrm,i,j,k)
          crm_output%sub_a (icrm,l) = crm_output%sub_a (icrm,l) + sub1a(icrm,i,j,k)
          crm_output%dep_a (icrm,l) = crm_output%dep_a (icrm,l) + dep1a(icrm,i,j,k)
          crm_output%con_a (icrm,l) = crm_output%con_a (icrm,l) + con1a(icrm,i,j,k)
        enddo
      enddo
    enddo

    ! note, rates are divded by dt to get mean rate over step
    crm_output%aut_a (icrm,:) = crm_output%aut_a (icrm,:) * factor_xyt / dt
    crm_output%acc_a (icrm,:) = crm_output%acc_a (icrm,:) * factor_xyt / dt
    crm_output%evpc_a(icrm,:) = crm_output%evpc_a(icrm,:) * factor_xyt / dt
    crm_output%evpr_a(icrm,:) = crm_output%evpr_a(icrm,:) * factor_xyt / dt
    crm_output%mlt_a (icrm,:) = crm_output%mlt_a (icrm,:) * factor_xyt / dt
    crm_output%sub_a (icrm,:) = crm_output%sub_a (icrm,:) * factor_xyt / dt
    crm_output%dep_a (icrm,:) = crm_output%dep_a (icrm,:) * factor_xyt / dt
    crm_output%con_a (icrm,:) = crm_output%con_a (icrm,:) * factor_xyt / dt
  enddo
#endif /* m2005 */

  !$acc parallel loop async(asyncid)
  do icrm = 1 , ncrms
    crm_output%precc (icrm) = 0.
    crm_output%precl (icrm) = 0.
    crm_output%precsc(icrm) = 0.
    crm_output%precsl(icrm) = 0.
  enddo

  !$acc parallel loop collapse(3) async(asyncid)
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
        !  output CRM level precip  so that individual CRM precip values are passed down
        !  to CLM.
        crm_pcp(icrm,i,j) = precsfc(icrm,i,j)/1000.      ! mm/s --> m/s
        crm_snw(icrm,i,j) = precssfc(icrm,i,j)/1000.     ! mm/s --> m/s
#endif
        if(precsfc(icrm,i,j).gt.10./86400.) then
           !$acc atomic update
           crm_output%precc (icrm) = crm_output%precc (icrm) + precsfc(icrm,i,j)
           !$acc atomic update
           crm_output%precsc(icrm) = crm_output%precsc(icrm) + precssfc(icrm,i,j)
        else
           !$acc atomic update
           crm_output%precl (icrm) = crm_output%precl (icrm) + precsfc(icrm,i,j)
           !$acc atomic update
           crm_output%precsl(icrm) = crm_output%precsl(icrm) + precssfc(icrm,i,j)
        endif
      enddo
    enddo
  enddo

  !$acc parallel loop collapse(3) async(asyncid)
  do j = 1 , ny
    do i = 1 , nx
      do icrm = 1 , ncrms
        crm_output%prec_crm(icrm,i,j) = precsfc(icrm,i,j)/1000.           !mm/s --> m/s
      enddo
    enddo
  enddo

  !$acc parallel loop async(asyncid)
  do icrm = 1 , ncrms
    crm_output%precc (icrm) = crm_output%precc (icrm)*factor_xy/1000.
    crm_output%precl (icrm) = crm_output%precl (icrm)*factor_xy/1000.
    crm_output%precsc(icrm) = crm_output%precsc(icrm)*factor_xy/1000.
    crm_output%precsl(icrm) = crm_output%precsl(icrm)*factor_xy/1000.

    crm_output%cltot(icrm) = crm_output%cltot(icrm) * factor_xyt
    crm_output%clhgh(icrm) = crm_output%clhgh(icrm) * factor_xyt
    crm_output%clmed(icrm) = crm_output%clmed(icrm) * factor_xyt
    crm_output%cllow(icrm) = crm_output%cllow(icrm) * factor_xyt

    crm_output%jt_crm(icrm) = plev * 1.0
    crm_output%mx_crm(icrm) = 1.0
  enddo

  !$acc parallel loop collapse(2) async(asyncid)
  do k=1, plev
    do icrm = 1 , ncrms
      crm_output%mu_crm(icrm,k)=0.5*(mui_crm(icrm,k)+mui_crm(icrm,k+1))
      crm_output%md_crm(icrm,k)=0.5*(mdi_crm(icrm,k)+mdi_crm(icrm,k+1))
      crm_output%mu_crm(icrm,k)=crm_output%mu_crm(icrm,k)*ggr/100.          !kg/m2/s --> mb/s
      crm_output%md_crm(icrm,k)=crm_output%md_crm(icrm,k)*ggr/100.          !kg/m2/s --> mb/s
      crm_output%eu_crm(icrm,k) = 0.
      if(mui_crm(icrm,k)-mui_crm(icrm,k+1).gt.0) then
        crm_output%eu_crm(icrm,k)=(mui_crm(icrm,k)-mui_crm(icrm,k+1))*ggr/crm_input%pdel(icrm,k)    !/s
      else
        crm_output%du_crm(icrm,k)=-1.0*(mui_crm(icrm,k)-mui_crm(icrm,k+1))*ggr/crm_input%pdel(icrm,k)   !/s
      endif
      if(mdi_crm(icrm,k+1)-mdi_crm(icrm,k).lt.0) then
        crm_output%ed_crm(icrm,k)=(mdi_crm(icrm,k)-mdi_crm(icrm,k+1))*ggr/crm_input%pdel(icrm,k) ! /s
      else
        dd_crm(icrm,k)=-1.*(mdi_crm(icrm,k)-mdi_crm(icrm,k+1))*ggr/crm_input%pdel(icrm,k)   !/s
      endif
      if(abs(crm_output%mu_crm(icrm,k)).gt.1.0e-15.or.abs(crm_output%md_crm(icrm,k)).gt.1.0e-15) then
        tmp = k
        !$acc atomic update
        crm_output%jt_crm(icrm) = min( tmp , crm_output%jt_crm(icrm) )

        tmp = k
        !$acc atomic update
        crm_output%mx_crm(icrm) = max( tmp , crm_output%mx_crm(icrm) )
      endif
    enddo
  enddo

  !-------------------------------------------------------------
  !       Fluxes and other stat:
  !-------------------------------------------------------------
  !$acc parallel loop collapse(2) async(asyncid)
  do k=1,nzm
    do icrm = 1 , ncrms
      u2z = 0.
      v2z = 0.
      w2z = 0.
      do j=1,ny
        do i=1,nx
          u2z = u2z+(u(icrm,i,j,k)-u0(icrm,k))**2
          v2z = v2z+(v(icrm,i,j,k)-v0(icrm,k))**2
          w2z = w2z+0.5*(w(icrm,i,j,k+1)**2+w(icrm,i,j,k)**2)
        enddo
      enddo
      !+++mhwang
      ! mkwsb, mkle, mkadv, mkdiff (also crm_output%flux_u, crm_output%flux_v,icrm) seem not calculted correclty in the spcam3.5 codes.
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
      crm_output%flux_u    (icrm,l) = (uwle(icrm,k) + uwsb(icrm,k))*tmp1*factor_xy/nstop
      crm_output%flux_v    (icrm,l) = (vwle(icrm,k) + vwsb(icrm,k))*tmp1*factor_xy/nstop
#ifdef sam1mom
      crm_output%flux_qt   (icrm,l) = mkwle(icrm,k,1) + mkwsb(icrm,k,1)
      crm_output%fluxsgs_qt(icrm,l) = mkwsb(icrm,k,1)
      crm_output%flux_qp   (icrm,l) = mkwle(icrm,k,2) + mkwsb(icrm,k,2)
      crm_output%qt_trans  (icrm,l) = mkadv(icrm,k,1) + mkdiff(icrm,k,1)
      crm_output%qp_trans  (icrm,l) = mkadv(icrm,k,2) + mkdiff(icrm,k,2)
#endif /* sam1mom */
#ifdef m2005
      crm_output%flux_qt   (icrm,l) = mkwle(icrm,k,1   ) + mkwsb(icrm,k,1   ) +  &
                         mkwle(icrm,k,iqci) + mkwsb(icrm,k,iqci)
      crm_output%fluxsgs_qt(icrm,l) = mkwsb(icrm,k,1   ) + mkwsb(icrm,k,iqci)
      crm_output%flux_qp   (icrm,l) = mkwle(icrm,k,iqr) + mkwsb(icrm,k,iqr) +  &
                         mkwle(icrm,k,iqs) + mkwsb(icrm,k,iqs) + mkwle(icrm,k,iqg) + mkwsb(icrm,k,iqg)
      crm_output%qt_trans  (icrm,l) = mkadv(icrm,k,1) + mkadv(icrm,k,iqci) + &
                         mkdiff(icrm,k,1) + mkdiff(icrm,k,iqci)
      crm_output%qp_trans  (icrm,l) = mkadv(icrm,k,iqr) + mkadv(icrm,k,iqs) + mkadv(icrm,k,iqg) + &
                         mkdiff(icrm,k,iqr) + mkdiff(icrm,k,iqs) + mkdiff(icrm,k,iqg)
#endif /* m2005 */
      tmp = 0
      do j = 1 , ny
        do i = 1 , nx
          tmp = tmp + sgs_field(icrm,i,j,k,1)
        enddo
      enddo
      crm_output%tkesgsz   (icrm,l)= rho(icrm,k)*tmp*factor_xy
      crm_output%tkez      (icrm,l)= rho(icrm,k)*0.5*(u2z+v2z*YES3D+w2z)*factor_xy + crm_output%tkesgsz(icrm,l)
      tmp = 0
      do j = 1 , ny
        do i = 1 , nx
          tmp = tmp + sgs_field_diag(icrm,i,j,k,1)
        enddo
      enddo
      crm_output%tkz       (icrm,l) = tmp * factor_xy
      crm_output%precflux  (icrm,l) = precflux(icrm,k)/1000.       !mm/s  -->m/s

      crm_output%qp_fall   (icrm,l) = qpfall(icrm,k)
      crm_output%qp_evp    (icrm,l) = qpevp(icrm,k)
      crm_output%qp_src    (icrm,l) = qpsrc(icrm,k)

      crm_output%qt_ls     (icrm,l) = qtend(icrm,k)
      crm_output%t_ls      (icrm,l) = ttend(icrm,k)
    enddo
  enddo

  !$acc wait(asyncid)

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

  crm_output%timing_factor(:) = crm_output%timing_factor(:) / nstop

#ifdef ECPP
  !!! Deallocate ECPP variables
  call ecpp_crm_cleanup ()
#endif

  deallocate( t00)
  deallocate( tln)
  deallocate( qln)
  deallocate( qccln)
  deallocate( qiiln)
  deallocate( uln)
  deallocate( vln)
#if defined(SP_ESMT)
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
#ifdef CLUBB_CRM
  deallocate( rtm_integral_before )
  deallocate( rtm_integral_after )
  deallocate( thlm_integral_before)
  deallocate( thlm_integral_after)
  deallocate( thlm_before)
  deallocate( thlm_after)
  deallocate( rtm_column)
#endif /* CLUBB_CRM */
  deallocate( dd_crm  )
  deallocate( mui_crm )
  deallocate( mdi_crm )
  deallocate( ustar )
  deallocate( bflx )
  deallocate( wnd )
  deallocate( qtot )
  deallocate( colprec  )
  deallocate( colprecs )

  call deallocate_params()
  call deallocate_grid()
  call deallocate_tracers()
  call deallocate_sgs()
  call deallocate_vars()
  call deallocate_micro()
#ifdef sam1mom
  call deallocate_micro_params()
#endif
#if defined( SP_ESMT )
  call deallocate_scalar_momentum()
#endif

end subroutine crm

end module crm_module
