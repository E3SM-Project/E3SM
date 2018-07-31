module shoc_intr

  !------------------------------------------------------------------- !
  ! Module to interface E3SM with Simplified Higher                    !
  ! Order Closure (SHOC), developed                                    !
  !    by Peter Bogenschutz (Bogenschutz and Krueger 2013).            !
  !                                                                    !
  ! SHOC replaces the exisiting turbulence, shallow convection, and    !
  !   macrophysics in E3SM                                             !  
  !                                                                    !  
  !                                                                    !
  !---------------------------Code history---------------------------- !
  ! Authors:  P. Bogenschutz                                           ! 
  !                                                                    ! 
  !------------------------------------------------------------------- !

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pver, pverp
  use physconst,     only: rair, cpair, gravit, latvap, latice, zvir, &
                           rh2o, karman, tms_orocnst, tms_z0fac  
  use constituents,  only: pcnst, cnst_add
  
  implicit none			   

  ! define physics buffer indicies here
  integer :: tke_idx, &     ! turbulent kinetic energy
  integer :: wthv_idx,&     ! buoyancy flux
  
  ! integer, public :: &
    ixtke = 0
    
  real, parameter :: tke_tol = (2.e-2_r8)**2
  
  real(r8), parameter :: &
      theta0   = 300._r8, &             ! Reference temperature                     [K]
      ts_nudge = 86400._r8, &           ! Time scale for u/v nudging (not used)     [s]
      p0_shoc = 100000._r8, &
      shoc_tk1 = 268.15_r8, &
      shoc_tk2 = 238.15_r8, &
      shoc_liq_deep = 8.e-6, &
      shoc_liq_sh = 10.e-6, &
      shoc_ice_deep = 25.e-6, &
      shoc_ice_sh = 50.e-6  
      
  logical      :: lq(pcnst)
  logical      :: lq2(pcnst)
  
  contains
  
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  
  subroutine shoc_register_e3sm()

#ifdef SHOC_SGS
    ! Add SHOC fields to pbuf
    use physics_buffer,  only: pbuf_add_field, dtype_r8, dyn_time_lvls
    use ppgrid,          only: pver, pverp, pcols  
  
    ! TKE is prognostic in SHOC and should be advected by dynamics
    call cnst_add('TKE',0._r8,0._r8,0._r8,longname='turbulent kinetic energy',cam_outfield=.true.)
  
    ! Fields that are not prognostic should be added to PBUF
    call pbuf_add_field('WTHV', 'global', dtype_r8, (/pcols,pverp,dyn_time_lvls/), wthv_idx) 
  
#endif
  
  end shoc_register_e3sm
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  
  subroutine shoc_init_cnst(name, q, gcid)
  
  !------------------------------------------------------------------- !
  ! Initialize the state for SHOC's prognostic variable                !
  !------------------------------------------------------------------- !
  
    character(len=*), intent(in)  :: name     ! constituent name
    real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
    integer,          intent(in)  :: gcid(:)  ! global column id
    
#ifdef SHOC_SGS
    if (trim(name) == trim('TKE')) q = tke_tol
#endif  

  end subroutine shoc_init_cnst
  
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  
  subroutine shoc_readnl(nlfile)
  
  !------------------------------------------------------------------- !
  ! Read in any namelist parameters here                               !
  !   (currently none)                                                 !
  !------------------------------------------------------------------- !  
  
  end subroutine shoc_readnl   
  
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  
  subroutine shoc_init_e3sm(pbuf2d)

  !------------------------------------------------------------------- !
  ! Initialize SHOC for E3SM                                           !
  !------------------------------------------------------------------- !  

    use physics_types,          only: physics_state, physics_ptend
    use ppgrid,                 only: pver, pverp, pcols
    use time_manager,              only: is_first_step
    
    implicit none
    !  Input Variables
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    
    lq(1:pcnst) = .true.
    edsclr_dim = pcnst
    
    ! Define physics buffers indexes
    cld_idx     = pbuf_get_index('CLD')         ! Cloud fraction
    concld_idx  = pbuf_get_index('CONCLD')      ! Convective cloud cover
    ast_idx     = pbuf_get_index('AST')         ! Stratiform cloud fraction
    alst_idx    = pbuf_get_index('ALST')        ! Liquid stratiform cloud fraction
    aist_idx    = pbuf_get_index('AIST')        ! Ice stratiform cloud fraction
    qlst_idx    = pbuf_get_index('QLST')        ! Physical in-stratus LWC 
    qist_idx    = pbuf_get_index('QIST')        ! Physical in-stratus IWC
    dp_frac_idx = pbuf_get_index('DP_FRAC')     ! Deep convection cloud fraction
    icwmrdp_idx = pbuf_get_index('ICWMRDP')     ! In-cloud deep convective mixing ratio
    sh_frac_idx = pbuf_get_index('SH_FRAC')     ! Shallow convection cloud fraction
    relvar_idx  = pbuf_get_index('RELVAR')      ! Relative cloud water variance
    accre_enhan_idx = pbuf_get_index('ACCRE_ENHAN') ! accretion enhancement for MG
    prer_evap_idx   = pbuf_get_index('PRER_EVAP')
    qrl_idx         = pbuf_get_index('QRL')
    cmfmc_sh_idx    = pbuf_get_index('CMFMC_SH')
    
    if (is_first_step()) then
      call pbuf_set_field(pbuf2, wthv_idx, 0.0_r8) 
    endif
    
    if (prog_modal_aero) then
       ! Turn off modal aerosols and decrement edsclr_dim accordingly
       call rad_cnst_get_info(0, nmodes=nmodes)
 
       do m = 1, nmodes
          call rad_cnst_get_mode_num_idx(m, lptr)
          lq(lptr)=.false.
          edsclr_dim = edsclr_dim-1
 
          call rad_cnst_get_info(0, m, nspec=nspec)
          do l = 1, nspec
             call rad_cnst_get_mam_mmr_idx(m, l, lptr)
             lq(lptr)=.false.
             edsclr_dim = edsclr_dim-1
          end do
       end do
 
       !  In addition, if running with MAM, droplet number is transported
       !  in dropmixnuc, therefore we do NOT want SHOC to apply transport
       !  tendencies to avoid double counted.  Else, we apply tendencies.
       call cnst_get_ind('NUMLIQ',ixnumliq)
       lq(ixnumliq) = .false.
       edsclr_dim = edsclr_dim-1
    endif    
    
    ! --------------- !
    ! End             !
    ! Initialization  !
    ! --------------- !    
  
  end subroutine shoc_init_e3sm   
  
  ! =============================================================================== !
  !                                                                                 !
  ! =============================================================================== !
  
  subroutine shoc_tend_e3sm( &
                             state, ptend_all, pbuf, hdtime, &
			     cmfmc, cam_in, sgh30, &
			     macmic_it, cld_macmic_num_steps, &
			     dlf, det_s, det_ice, alst_o)
  
  !------------------------------------------------------------------- !
  ! Provide tendencies of shallow convection , turbulence, and         !
  !   macrophysics from SHOC to E3SM                                   !
  !------------------------------------------------------------------- !  
  
    use physics_types,  only: physics_state, physics_ptend, &
                              physics_state_copy, physics_ptend_init, &
                              physics_ptend_sum 
			      
    use physics_update_mod, only: physics_update

    use physics_buffer, only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                              pbuf_set_field, physics_buffer_desc	
			      
    use ppgrid,         only: pver, pverp, pcols
    use constituents,   only: cnst_get_ind
    use camsrfexch,     only: cam_in_t
    use ref_pres,       only: top_lev => trop_cloud_top_lev  
    use time_manager,   only: is_first_step   
    use cam_abortutils, only: endrun
    use wv_saturation,  only: qsat
    use micro_mg_cam,   only: micro_mg_version    	
    
    implicit none
    
   ! --------------- !
   ! Input Auguments !
   ! --------------- !

   type(physics_state), intent(in)    :: state                    ! Physics state variables                 [vary]
   type(cam_in_t),      intent(in)    :: cam_in
   real(r8),            intent(in)    :: hdtime                   ! Host model timestep                     [s]
   real(r8),            intent(in)    :: dlf(pcols,pver)          ! Detraining cld H20 from deep convection [kg/ks/s]
   real(r8),            intent(in)    :: cmfmc(pcols,pverp)       ! convective mass flux--m sub c           [kg/m2/s]
   real(r8),            intent(in)    :: sgh30(pcols)             ! std deviation of orography              [m]
   integer,             intent(in)    :: cld_macmic_num_steps     ! number of mac-mic iterations
   integer,             intent(in)    :: macmic_it                ! number of mac-mic iterations    		      		  
   ! ---------------------- !
   ! Input-Output Auguments !
   ! ---------------------- !
    
   type(physics_buffer_desc), pointer :: pbuf(:)

   ! ---------------------- !
   ! Output Auguments !
   ! ---------------------- !

   type(physics_ptend), intent(out)   :: ptend_all                 ! package tendencies

   ! These two variables are needed for energy check    
   real(r8),            intent(out)   :: det_s(pcols)              ! Integral of detrained static energy from ice
   real(r8),            intent(out)   :: det_ice(pcols)            ! Integral of detrained ice for energy check

   real(r8), intent(out) :: alst_o(pcols,pver)  ! H. Wang: for old liquid status fraction 
   
   ! --------------- !
   ! Local Variables !
   ! --------------- !
   
#ifdef SHOC_SGS

   type(physics_state) :: state1                ! Local copy of state variable
   type(physics_ptend) :: ptend_loc             ! Local tendency from processes, added up to return as ptend_all

   integer :: i, j, k, t, ixind, nadv
   integer :: ixcldice, ixcldliq, ixnumliq, ixnumice, ixq
   integer :: itim_old
   integer :: ncol, lchnk                       ! # of columns, and chunk identifier
   integer :: err_code                          ! Diagnostic, for if some calculation goes amiss.
   integer :: begin_height, end_height
   integer :: icnt, clubbtop
   
   real(r8) :: dtime                            ! SHOC time step                              [s]   
   real(r8) :: edsclr_in(pverp,edsclr_dim)      ! Scalars to be diffused through SHOC         [units vary]   
   real(r8) :: tke_in(pcols,pverp)
   real(r8) :: thlm_in(pcols,pverp)
   real(r8) :: qv_in(pcols,pverp)
   real(r8) :: rcm_in(pcols,pverp)
   real(r8) :: pres_in(pcols,pverp)
   real(r8) :: um_in(pcols,pverp)
   real(r8) :: vm_in(pcols,pverp)
   real(r8) :: cloudfrac_shoc(pcols,pverp)
   real(r8) :: rcm_shoc(pcols,pverp)
   
   ! Variables below are needed to compute energy integrals for conservation
   real(r8) :: ke_a(pcols), ke_b(pcols), te_a(pcols), te_b(pcols)
   real(r8) :: wv_a(pcols), wv_b(pcols), wl_b(pcols), wl_a(pcols)
   real(r8) :: se_dis, se_a(pcols), se_b(pcols), clubb_s(pver)
   
   ! --------------- !
   ! Pointers        !
   ! --------------- !
   
   real(r8), pointer, dimension(:,:) :: wthv ! buoyancy flux
   real(r8), pointer, dimension(:,:) :: um       ! mean east-west wind                          [m/s]
   real(r8), pointer, dimension(:,:) :: vm       ! mean north-south wind                        [m/s]
   real(r8), pointer, dimension(:,:) :: cld      ! cloud fraction                               [fraction]
   real(r8), pointer, dimension(:,:) :: concld   ! convective cloud fraction                    [fraction]
   real(r8), pointer, dimension(:,:) :: ast      ! stratiform cloud fraction                    [fraction]
   real(r8), pointer, dimension(:,:) :: alst     ! liquid stratiform cloud fraction             [fraction]
   real(r8), pointer, dimension(:,:) :: aist     ! ice stratiform cloud fraction                [fraction]         
   real(r8), pointer, dimension(:,:) :: qlst     ! Physical in-stratus LWC                      [kg/kg]
   real(r8), pointer, dimension(:,:) :: qist     ! Physical in-stratus IWC                      [kg/kg]
   real(r8), pointer, dimension(:,:) :: deepcu   ! deep convection cloud fraction               [fraction]
   real(r8), pointer, dimension(:,:) :: shalcu   ! shallow convection cloud fraction            [fraction]    
   real(r8), pointer, dimension(:,:) :: khzt     ! eddy diffusivity on thermo levels            [m^2/s]
   real(r8), pointer, dimension(:,:) :: khzm     ! eddy diffusivity on momentum levels          [m^2/s]
   real(r8), pointer, dimension(:) :: pblh     ! planetary boundary layer height                [m]

#endif
   det_s(:)   = 0.0_r8
   det_ice(:) = 0.0_r8
#ifdef SHOC_SGS
   
   !------------------------------------------------------------------!
   !------------------------------------------------------------------!
   !------------------------------------------------------------------!
   !       MAIN COMPUTATION BEGINS HERE                               !
   !------------------------------------------------------------------!
   !------------------------------------------------------------------!
   !------------------------------------------------------------------!   
   
 !  Get indicees for cloud and ice mass and cloud and ice number

   call cnst_get_ind('Q',ixq)
   call cnst_get_ind('CLDLIQ',ixcldliq)
   call cnst_get_ind('CLDICE',ixcldice)
   call cnst_get_ind('NUMLIQ',ixnumliq)
   call cnst_get_ind('NUMICE',ixnumice)
   
   call physics_ptend_init(ptend_loc,state%psetcols, 'clubb_ice1', ls=.true., lu=.true., lv=.true., lq=lq)
   
   call physics_state_copy(state,state1)
   
   !  Determine number of columns and which chunk computation is to be performed on
   ncol = state%ncol
   lchnk = state%lchnk    
   
   !  Determine time step of physics buffer 
   itim_old = pbuf_old_tim_idx()     
   
   !  Establish associations between pointers and physics buffer fields   
   call pbuf_get_field(pbuf, wthv_idx,     wthv,     start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))  
   call pbuf_get_field(pbuf, cld_idx,     cld,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, concld_idx,  concld,  start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, ast_idx,     ast,     start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, alst_idx,    alst,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, aist_idx,    aist,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, qlst_idx,    qlst,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   call pbuf_get_field(pbuf, qist_idx,    qist,    start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
   
   call pbuf_get_field(pbuf, prer_evap_idx, prer_evap)
   call pbuf_get_field(pbuf, accre_enhan_idx, accre_enhan)
   call pbuf_get_field(pbuf, cmeliq_idx,  cmeliq)
   call pbuf_get_field(pbuf, relvar_idx,  relvar)
   call pbuf_get_field(pbuf, dp_frac_idx, deepcu)
   call pbuf_get_field(pbuf, sh_frac_idx, shalcu)
   call pbuf_get_field(pbuf, kvm_idx,     khzt)
   call pbuf_get_field(pbuf, kvh_idx,     khzm)
   call pbuf_get_field(pbuf, pblh_idx,    pblh)
   call pbuf_get_field(pbuf, icwmrdp_idx, dp_icwmr)
   call pbuf_get_field(pbuf, cmfmc_sh_idx, cmfmc_sh)  
   
   !  Determine SHOC time step and make it sub-step friendly
   !  For now we want SHOC time step to be 5 min.  However, there are certain
   !  instances when a 5 min time step will not be possible (based on 
   !  host model time step or on macro-micro sub-stepping   
   
   dtime = 300._r8  !+DPAB, probably want to make this a namelist variable
   
   !  Now check to see if dtime is greater than the host model 
   !    (or sub stepped) time step.  If it is, then simply 
   !    set it equal to the host (or sub step) time step.  
   !    This section is mostly to deal with small host model
   !    time steps (or small sub-steps)
   
   if (dtime .gt. hdtime) then
     dtime = hdtime
   endif
   
   !  Now check to see if SHOC time step divides evenly into
   !    the host model time step.  If not, force it to divide evenly.
   !    We also want it to be 5 minutes or less.  This section is
   !    mainly for host model time steps that are not evenly divisible
   !    by 5 minutes  
   
   if (mod(hdtime,dtime) .ne. 0) then
     dtime = hdtime/2._r8
     do while (dtime .gt. 300._r8) 
       dtime = dtime/2._r8
     end do
   endif  
   
   !  If resulting host model time step and SHOC time step do not divide evenly
   !    into each other, have model throw a fit.  

   if (mod(hdtime,dtime) .ne. 0) then
     call endrun('shoc_tend_e3sm:  SHOC time step and HOST time step NOT compatible')
   endif      
   
   minqn = 0._r8
   newfice(:,:) = 0._r8
   where(state1%q(:ncol,:pver,3) .gt. minqn) &
       newfice(:ncol,:pver) = state1%q(:ncol,:pver,3)/(state1%q(:ncol,:pver,2)+state1%q(:ncol,:pver,3))  
       
   do k=1,pver
     do i=1,ncol
       exner(i,k) = 1._r8/((state1%pmid(i,k)/p0_shoc)**(rair/cpair))
     enddo
   enddo       
       
   !  At each SHOC call, initialize mean momentum  and thermo SHOC state 
   !  from the E3SM state
   
   do k=1,pver   ! loop over levels
     do i=1,ncol ! loop over columns
     
       rvm(i,k) = state1%q(i,k,ixq)
       rcm(i,k) = state1%q(i,k,ixcldliq)
       um(i,k) = state1%u(i,k)
       vm(i,k) = state1%v(i,k)
       thlm(i,k) = state1%t(i,k)*exner(i,k)-(latvap/cpair)*state1%q(i,k,ixcldliq)
       
       if (macmic_it .eq. 1) then
         tke(i,k) = state1%q(i,k,ixtke)
       endif
     
     enddo
   enddo    
   
   rvm(1:ncol,pverp) = rtm(1:ncol,pver)
   rcm(1:ncol,pverp) = rcm(1:ncol,pver)
   um(1:ncol,pverp) = um(1:ncol,pver)
   vm(1:ncol,pverp) = vm(1:ncol,pver)
   thlm(1:ncol,pverp) = thlm(1:ncol,pver) 
   tke(1:col,pverp) = tke(1:ncol,pver)     
   
   ! Compute integrals of static energy, kinetic energy, water vapor, and liquid water
   ! for the computation of total energy before CLUBB is called.  This is for an 
   ! effort to conserve energy since liquid water potential temperature (which CLUBB 
   ! conserves) and static energy (which CAM conserves) are not exactly equal.   
   se_b = 0._r8
   ke_b = 0._r8
   wv_b = 0._r8
   wl_b = 0._r8
   do k=1,pver
     do i=1,ncol
       se_b(i) = se_b(i) + state1%s(i,k)*state1%pdel(i,k)/gravit
       ke_b(i) = ke_b(i) + 0.5_r8*(um(i,k)**2+vm(i,k)**2)*state1%pdel(i,k)/gravit
       wv_b(i) = wv_b(i) + state1%q(i,k,ixq)*state1%pdel(i,k)/gravit
       wl_b(i) = wl_b(i) + state1%q(i,k,ixcldliq)*state1%pdel(i,k)/gravit
     enddo
   enddo
   
   ! ------------------------------------------------- !
   ! Begin module to compute turbulent mountain stress !
   ! ------------------------------------------------- !
    if ( do_tms) then
       call t_startf('compute_tms')
       call compute_tms( pcols,        pver,      ncol,                   &
                     state1%u,     state1%v,  state1%t,  state1%pmid, &
                     state1%exner, state1%zm, sgh30,     ksrftms,     &
                     tautmsx,      tautmsy,   cam_in%landfrac ) 
       call t_stopf('compute_tms')
    endif
    
   ! ------------------------------------------------- !
   ! End module to compute turbulent mountain stress   !
   ! ------------------------------------------------- !  
   
   ! ------------------------------------------------- !
   ! Prepare inputs for SHOC call                      !
   ! ------------------------------------------------- ! 
   
   !  Define the CLUBB thermodynamic grid (in units of m)
   wm_zt(:,1) = 0._r8
   do k=1,pver
     do i=1,ncol
       zt_g(i,k+1) = state1%zm(i,pver-k+1)-state1%zi(i,pver+1)
       wm_zt(i,k+1) = -1._r8*state1%omega(i,pver-k+1)/(rho(i,k+1)*gravit)
     enddo
   enddo
     
   do k=1,pverp
     do i=1,ncol
       zi_g(k) = state1%zi(i,pverp-k+1)-state1%zi(i,pver+1)
     enddo
   enddo       
   
   zt_g(:,1) = -1._r8*zt_g(:,2)
   
   ! Figure out way to deal with surface fluxes
   ! and TMS
   
      !  Surface fluxes provided by host model
!      wpthlp_sfc = cam_in%shf(i)/(cpair*rho_ds_zm(1))       ! Sensible heat flux
!      wprtp_sfc  = cam_in%cflx(i,1)/(rho_ds_zm(1))      ! Latent heat flux
!      upwp_sfc   = cam_in%wsx(i)/rho_ds_zm(1)               ! Surface meridional momentum flux
!      vpwp_sfc   = cam_in%wsy(i)/rho_ds_zm(1)               ! Surface zonal momentum flux  
      
      ! ------------------------------------------------- !
      ! Apply TMS                                         !
      ! ------------------------------------------------- !    
!       if ( do_tms) then
!          upwp_sfc = upwp_sfc-((ksrftms(i)*state1%u(i,pver))/rho_ds_zm(1))
!          vpwp_sfc = vpwp_sfc-((ksrftms(i)*state1%v(i,pver))/rho_ds_zm(1)) 
!	endif           

    ! Need to flip arrays around for SHOC 
    do k=1,pverp
      do i=1,ncol
        um_in(i,k)      = um(i,pverp-k+1)
        vm_in(i,k)      = vm(i,pverp-k+1)   
        rvm_in(i,k)     = rvm(i,pverp-k+1)
        rcm_in(i,k)     = rcm(i,pverp-k+1)
        thlm_in(i,k)    = thlm(i,pverp-k+1)
        tke_in(i,k)     = tke(i,pverp-k+1)
	wthv_in(i,k)    = wthv(i,pverp-k+1)
        pres_in(i,k)    = state1%pmid(i,pver-k+1)
      enddo  
    enddo   
    
   !  Do the same for tracers 
   icnt=0
   do ixind=1,pcnst
     if (lq(ixind))  then 
       icnt=icnt+1
       do k=1,pver
         edsclr_in(k+1,icnt) = state1%q(i,pver-k+1,ixind)
       enddo
       edsclr_in(1,icnt) = edsclr_in(2,icnt)
     end if
   enddo    
    
   ! ------------------------------------------------- !
   ! Actually call SHOC                                !
   ! ------------------------------------------------- !   
   
   do t=1,nadv
   
     call shoc_main &
          (pcols, pverp, dtime, &
	  zt_g, zi_g, pres_in, &
	  tke_in, thlm_in, rvm_in, wm_zt, &
	  um_in, vm_in rcm_in, edsclr_in, &
	  edsclr_dim, wthv_in, &
	  cloudfrac_shoc, rcm_shoc) 
   
   enddo  ! end time loop
   
   ! Arrays need to be "flipped" to CAM grid
   
   do k=1,pverp
     do i=1,ncol 
       um(i,k) = um_in(pverp-k+1)
       vm(i,k) = vm_in(pverp-k+1)
       thlm(i,k) = thlm_in(pverp-k+1)
       rvm(i,k) = rvm_in(pverp-k+1)
       rcm(i,k) = rcm_shoc(pverp-k+1)
       cloud_frac(i,k) = min(cloudfrac_shoc(pverp-k+1),1._r8)
       wthv(i,k) = wthv_in(pverp-k+1)
       tke(i,k) = tke_in(pverp-k+1)
       
       do ixind=1,edsclr_dim
         edsclr_out(k,ixind) = edsclr_in(pverp-k+1,ixind)
       enddo       
       
     enddo
   enddo
   
   ! Compute integrals for static energy, kinetic energy, water vapor, and liquid water
   ! after SHOC is called.  This is for energy conservation purposes. 
   se_a = 0._r8
   ke_a = 0._r8
   wv_a = 0._r8
   wl_a = 0._r8
   do k=1,pver
     do i=1,ncol
       clubb_s(i,k) = cpair*((thlm(i,k)+(latvap/cpair)*rcm(i,k))/exner(i,k))+ &
                      gravit*state1%zm(i,k)+state1%phis(i)
       se_a(i) = se_a(i) + clubb_s(i,k)*state1%pdel(i,k)/gravit
       ke_a(i) = ke_a(i) + 0.5_r8*(um(i,k)**2+vm(i,k)**2)*state1%pdel(i,k)/gravit
       wv_a(i) = wv_a(i) + (rvm(i,k)*state1%pdel(i,k)/gravit
       wl_a(i) = wl_a(i) + (rcm(i,k))*state1%pdel(i,k)/gravit
     enddo    
   enddo     
   
   ! Based on these integrals, compute the total energy before and after CLUBB call
   do i=1,ncol
     te_a(i) = se_a(i) + ke_a(i) + (latvap+latice)*wv_a(i)+latice*wl_a(i)
     te_b(i) = se_b(i) + ke_b(i) + (latvap+latice)*wv_b(i)+latice*wl_b(i)
     te_b(i) = te_b(i)+(cam_in%shf(i)+(cam_in%cflx(i,1))*(latvap+latice))*hdtime
   enddo  
   
   ! Limit the energy fixer to find highest layer where CLUBB is active
   ! Find first level where wp2 is higher than lowest threshold
   do i=1,ncol
      shoctop(i) = 1
     do while (tke(i,clubbtop) .eq. tke_tol .and. shoctop .lt. pver-1)
       shoctop(i) = shoctop(i) + 1
     enddo   
   
     ! Compute the disbalance of total energy, over depth where CLUBB is active
     se_dis(i) = (te_a(i) - te_b(i))/(state1%pint(i,pverp)-state1%pint(i,shoctop(i)))  
   enddo    
   
   do k=1,pver
     do i=1,ncol
       shoc_s(i,k) = shoc_s(i,k) - se_dis(i)*gravit
     enddo
   enddo
   
   !  Now compute the tendencies of CLUBB to CAM, note that pverp is the ghost point
   !  for all variables and therefore is never called in this loop
   do k=1,pver
     do i=1,ncol
       
       ptend_loc%u(i,k) = (um(i,k)-state1%u(i,k))/hdtime
       ptend_loc%v(i,k)   = (vm(i,k)-state1%v(i,k))/hdtime           
       ptend_loc%q(i,k,ixq) = (rvm(i,k)-state1%q(i,k,ixq))/hdtime ! water vapor
       ptend_loc%q(i,k,ixcldliq) = (rcm(i,k)-state1%q(i,k,ixcldliq))/hdtime   ! Tendency of liquid water
       ptend_loc%s(i,k) = (shoc_s(i,k)-state1%s(i,k))/hdtime
       
       if (macmic_it .eq. cld_macmic_num_steps) then
         ptend_loc%q(i,k,ixtke)=(tke(i,k)-state1%q(i,k,ixtke))/hdtime ! TKE
       endif
       
     enddo
   enddo   
   
   !  Apply tendencies to ice mixing ratio, liquid and ice number, and aerosol constituents.
   !  Loading up this array doesn't mean the tendencies are applied.  
   ! edsclr_out is compressed with just the constituents being used, ptend and state are not compressed

   icnt=0
   do ixind=1,pcnst
     if (lq(ixind)) then
       icnt=icnt+1
       if ((ixind /= ixq)       .and. (ixind /= ixcldliq) .and.&
          (ixind /= ixtke) )
          ptend_loc%q(i,k,ixind) = (edsclr_out(k,icnt)-state1%q(i,k,ixind))/hdtime ! transported constituents 
       end if
     end if
   enddo   
   
   cmeliq(:,:) = ptend_loc%q(:,:,ixcldliq)
   
   ! Update physics tendencies
   call physics_ptend_init(ptend_all, state%psetcols, 'shoc')
   call physics_ptend_sum(ptend_loc,ptend_all,ncol)
   call physics_update(state1,ptend_loc,hdtime)
   
   ! ------------------------------------------------------------ !
   ! ------------------------------------------------------------ ! 
   ! ------------------------------------------------------------ !
   ! The rest of the code deals with diagnosing variables         !
   ! for microphysics/radiation computation and macrophysics      !
   ! ------------------------------------------------------------ !
   ! ------------------------------------------------------------ !
   ! ------------------------------------------------------------ !   
   
   ! --------------------------------------------------------------------------------- !  
   !  COMPUTE THE ICE CLOUD DETRAINMENT                                                !
   !  Detrainment of convective condensate into the environment or stratiform cloud    !
   ! --------------------------------------------------------------------------------- !
   
   call physics_ptend_init(ptend_loc,state%psetcols, 'clubb_det', ls=.true., lq=lqice)   
   do k=1,pver
      do i=1,ncol
         if( state1%t(i,k) > shoc_tk1 ) then
            dum1 = 0.0_r8
         elseif ( state1%t(i,k) < shoc_tk2 ) then
            dum1 = 1.0_r8
         else
            !Note: Denominator is changed from 30.0_r8 to (shoc_tk1 - shoc_tk2),
            !(clubb_tk1 - clubb_tk2) is also 30.0 but it introduced a non-bfb change
            dum1 = ( shoc_tk1 - state1%t(i,k) ) /(shoc_tk1 - shoc_tk2)
         endif
        
         ptend_loc%q(i,k,ixcldliq) = dlf(i,k) * ( 1._r8 - dum1 )
         ptend_loc%q(i,k,ixcldice) = dlf(i,k) * dum1
         ptend_loc%q(i,k,ixnumliq) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) * ( 1._r8 - dum1 ) ) &
                                     / (4._r8*3.14_r8* shoc_liq_deep**3*997._r8) + & ! Deep    Convection
                                     3._r8 * (                         dlf2(i,k)    * ( 1._r8 - dum1 ) ) &
                                     / (4._r8*3.14_r8*shoc_liq_sh**3*997._r8)     ! Shallow Convection 
         ptend_loc%q(i,k,ixnumice) = 3._r8 * ( max(0._r8, ( dlf(i,k) - dlf2(i,k) )) *  dum1 ) &
                                     / (4._r8*3.14_r8*shoc_ice_deep**3*500._r8) + & ! Deep    Convection
                                     3._r8 * (                         dlf2(i,k)    *  dum1 ) &
                                     / (4._r8*3.14_r8*shoc_ice_sh**3*500._r8)     ! Shallow Convection
         ptend_loc%s(i,k)          = dlf(i,k) * dum1 * latice
 
         ! Only rliq is saved from deep convection, which is the reserved liquid.  We need to keep
         !   track of the integrals of ice and static energy that is effected from conversion to ice
         !   so that the energy checker doesn't complain.
         det_s(i)                  = det_s(i) + ptend_loc%s(i,k)*state1%pdel(i,k)/gravit
         det_ice(i)                = det_ice(i) - ptend_loc%q(i,k,ixcldice)*state1%pdel(i,k)/gravit
 
      enddo
    enddo

    det_ice(:ncol) = det_ice(:ncol)/1000._r8  ! divide by density of water
   
    call physics_ptend_sum(ptend_loc,ptend_all,ncol)
    call physics_update(state1,ptend_loc,hdtime)
   
    ! For purposes of this implementaiton, just set relvar and accre_enhan to 1
    relvar(:,:) = 1._r8   
    accre_enhan(:,:) = 1._r8  
   
    ! --------------------------------------------------------------------------------- ! 
    !  Diagnose some quantities that are computed in macrop_tend here.                  !
    !  These are inputs required for the microphysics calculation.                      !
    !                                                                                   !
    !  FIRST PART COMPUTES THE STRATIFORM CLOUD FRACTION FROM SHOC CLOUD FRACTION       !
    ! --------------------------------------------------------------------------------- ! 
   
    ! HW: set alst to alst_o before getting updated
    if(liqcf_fix) then
      if(.not.is_first_step()) alst_o(:ncol,:pver) = alst(:ncol,:pver)
    endif

    !  initialize variables 
    alst(:,:) = 0.0_r8
    qlst(:,:) = 0.0_r8 
 
    do k=1,pver
      do i=1,ncol
        alst(i,k) = cloud_frac(i,k)   
        qlst(i,k) = rcm(i,k)/max(0.01_r8,alst(i,k))  ! Incloud stratus condensate mixing ratio
      enddo
    enddo

    ! HW
    if(liqcf_fix) then
      if(is_first_step()) alst_o(:ncol,:pver) = alst(:ncol,:pver)
    endif  
   
    ! --------------------------------------------------------------------------------- !  
    !  THIS PART COMPUTES CONVECTIVE AND DEEP CONVECTIVE CLOUD FRACTION                 !
    ! --------------------------------------------------------------------------------- ! 
 
    deepcu(:,pver) = 0.0_r8
    shalcu(:,pver) = 0.0_r8
 
    do k=1,pver-1
      do i=1,ncol
        !  diagnose the deep convective cloud fraction, as done in macrophysics based on the 
        !  deep convective mass flux, read in from pbuf.  Since shallow convection is never 
        !  called, the shallow convective mass flux will ALWAYS be zero, ensuring that this cloud
        !  fraction is purely from deep convection scheme.  
        deepcu(i,k) = max(0.0_r8,min(dp1*log(1.0_r8+500.0_r8*(cmfmc(i,k+1)-cmfmc_sh(i,k+1))),0.6_r8))
        shalcu(i,k) = 0._r8
       
        if (deepcu(i,k) <= frac_limit .or. dp_icwmr(i,k) < ic_limit) then
          deepcu(i,k) = 0._r8
        endif
             
        !  using the deep convective cloud fraction, and SHOC cloud fraction (variable 
        !  "cloud_frac"), compute the convective cloud fraction.  This follows the formulation
        !  found in macrophysics code.  Assumes that convective cloud is all nonstratiform cloud 
        !  from CLUBB plus the deep convective cloud fraction
        concld(i,k) = min(cloud_frac(i,k)-alst(i,k)+deepcu(i,k),0.80_r8)
      enddo
    enddo   
   
    ! --------------------------------------------------------------------------------- !  
    !  COMPUTE THE ICE CLOUD FRACTION PORTION                                           !
    !  use the aist_vector function to compute the ice cloud fraction                   !
    ! --------------------------------------------------------------------------------- !
   
    do k=1,pver
      call aist_vector(state1%q(:,k,ixq),state1%t(:,k),state1%pmid(:,k),state1%q(:,k,ixcldice), &
           state1%q(:,k,ixnumice),cam_in%landfrac(:),cam_in%snowhland(:),aist(:,k),ncol)
    enddo
   
    ! --------------------------------------------------------------------------------- !  
    !  THIS PART COMPUTES THE LIQUID STRATUS FRACTION                                   !
    !                                                                                   !
    !  For now leave the computation of ice stratus fraction from macrop_driver intact  !
    !  because SHOC does nothing with ice.  Here I simply overwrite the liquid stratus ! 
    !  fraction that was coded in macrop_driver                                         !
    ! --------------------------------------------------------------------------------- !  
 
    !  Recompute net stratus fraction using maximum over-lapping assumption, as done
    !  in macrophysics code, using alst computed above and aist read in from physics buffer            

    cldthresh=1.e-18_r8

    do k=1,pver
      do i=1,ncol

        ast(i,k) = max(alst(i,k),aist(i,k))

        qist(i,k) = state1%q(i,k,ixcldice)/max(0.01_r8,aist(i,k)) 
      enddo
    enddo
   
    !  Probably need to add deepcu cloud fraction to the cloud fraction array, else would just 
    !  be outputting the shallow convective cloud fraction 

    do k=1,pver
      do i=1,ncol
        cloud_frac(i,k) = min(ast(i,k)+deepcu(i,k),1.0_r8)
      enddo
    enddo
   
    ! --------------------------------------------------------------------------------- !  
    !  DIAGNOSE THE PBL DEPTH                                                           !
    !  this is needed for aerosol code                                                  !
    ! --------------------------------------------------------------------------------- ! 

    do i=1,ncol
      do k=1,pver
        th(i,k) = state1%t(i,k)*state1%exner(i,k)
        thv(i,k) = th(i,k)*(1.0_r8+zvir*state1%q(i,k,ixq))
      enddo
    enddo
 
    ! diagnose surface friction and obukhov length (inputs to diagnose PBL depth)
    do i=1,ncol
      rrho = (1._r8/gravit)*(state1%pdel(i,pver)/dz_g(pver))
      call calc_ustar( state1%t(i,pver), state1%pmid(i,pver), cam_in%wsx(i), cam_in%wsy(i), &
                       rrho, ustar2(i) )
      call calc_obklen( th(i,pver), thv(i,pver), cam_in%cflx(i,1), cam_in%shf(i), rrho, ustar2(i), &
                        kinheat(i), kinwat(i), kbfs(i), obklen(i) )  
    enddo
   
    dummy2(:) = 0._r8
    dummy3(:) = 0._r8
   
    where (kbfs .eq. -0.0_r8) kbfs = 0.0_r8

    !  Compute PBL depth according to Holtslag-Boville Scheme
    call pblintd(ncol, thv, state1%zm, state1%u, state1%v, &
                ustar2, obklen, kbfs, pblh, dummy2, &
                state1%zi, cloud_frac(:,1:pver), 1._r8-cam_in%landfrac, dummy3)  
		
    ! Assign the first pver levels of cloud_frac back to cld
    cld(:,1:pver) = cloud_frac(:,1:pver)		 
    return
#endif             
  end subroutine shoc_tend_e3sm      

end module shoc_intr
