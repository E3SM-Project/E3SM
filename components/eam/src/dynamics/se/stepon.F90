!BOP
!
! !MODULE: stepon -- FV Dynamics specific time-stepping
!
! !INTERFACE:
module stepon

! !USES:
! from cam
   use shr_kind_mod,   only: r8 => shr_kind_r8
   use shr_sys_mod,    only: shr_sys_flush
   use pmgrid,         only: plev, plevp, plat
   use spmd_utils,     only: iam, masterproc, mpicom
   use constituents,   only: pcnst, cnst_name, cnst_longname
   use cam_abortutils, only: endrun
   use ppgrid,         only: begchunk, endchunk
   use physconst,      only: zvir, cappa
   use physics_types,  only: physics_state, physics_tend
   use dyn_comp,       only: dyn_import_t, dyn_export_t
   use perf_mod,       only: t_startf, t_stopf, t_barrierf
   use time_manager,   only: get_step_size
! from SE
   use derivative_mod, only: derivinit, derivative_t
   use viscosity_mod, only : compute_zeta_C0, compute_div_C0
   use quadrature_mod, only: gauss, gausslobatto, quadrature_t
   use edge_mod,       only: edge_g, edgeVpack_nlyr, edgeVunpack_nlyr
   use parallel_mod,   only : par
   use scamMod,        only: use_iop, doiopupdate, single_column, &
                             setiopupdate, readiopdata
   use element_mod,    only: element_t
   use shr_const_mod,       only: SHR_CONST_PI

   implicit none
   private
   save

!
! !PUBLIC MEMBER FUNCTIONS: 
!
  public stepon_init   ! Initialization
  public stepon_run1    ! run method phase 1
  public stepon_run2    ! run method phase 2
  public stepon_run3    ! run method phase 3
  public stepon_final  ! Finalization

!----------------------------------------------------------------------
!
! !DESCRIPTION: Module for dynamics time-stepping.
!
! !REVISION HISTORY:
!
! 2006.05.31  JPE    Created
! 2019.06.28  MT     Updated to output vorticity/divergence, new DSS interface
!
!EOP
!----------------------------------------------------------------------
!BOC
!
! !PRIVATE DATA MEMBERS:
!

  logical :: iop_update_phase1

  type (derivative_t)   :: deriv           ! derivative struct
  type (quadrature_t)   :: gv,gp           ! quadratures on velocity and pressure grids
!-----------------------------------------------------------------------


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_init --- Time stepping initialization
!
! !INTERFACE:
subroutine stepon_init(dyn_in, dyn_out )
! !USES:
  use dimensions_mod,         only: nlev, nelemd, npsq
  use dyn_grid,               only: fv_nphys
  use cam_history,            only: addfld, add_default, horiz_only
  use cam_history,            only: register_vector_field
  use gravity_waves_sources,  only: gws_init
  use phys_control,           only: use_gw_front
  use cam_history_support,    only: max_fieldname_len

! !OUTPUT PARAMETERS
!
  type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
  type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container

  integer :: m
  character(len=max_fieldname_len) :: grid_name
! !DESCRIPTION:
!
! Allocate data, initialize values, setup grid locations and other
! work needed to prepare the dynamics to run. Return weights and 
! vertical coords to atmosphere calling program.
!
!EOP
!-----------------------------------------------------------------------
!BOC

  ! This is not done in dyn_init due to a circular dependency issue.
  if(par%dynproc) then
     if (use_gw_front)  call gws_init(dyn_in%elem)
  end if

  ! fields that are written by the dycore
  ! these calls cant be done in dyn_init() because physics grid
  ! is not initialized at that point if making a restart runs
  !
  ! Forcing from physics
  call addfld ('FU',  (/ 'lev' /), 'A', 'm/s2', 'Zonal wind forcing term',     gridname='GLL')
  call addfld ('FV',  (/ 'lev' /), 'A', 'm/s2', 'Meridional wind forcing term',gridname='GLL')
  call register_vector_field('FU', 'FV')
  call addfld ('VOR', (/ 'lev' /), 'A', '1/s',  'Relative Vorticity (2D)',     gridname='GLL')
  call addfld ('DIV', (/ 'lev' /), 'A', '1/s',  'Divergence (2D)',             gridname='GLL')
  call addfld ('DIV_Qflux', (/ 'lev' /), 'A', '1/s kg/m2',  'Divergence of Qdp (2D)', gridname='GLL') !(zhang73)

  call addfld ('ETADOT', (/ 'ilev' /), 'A', '1/s', 'Vertical (eta) velocity', gridname='physgrid')

  if (fv_nphys > 0) then
    grid_name = 'GLL'
  else
    grid_name = 'physgrid'
  end if 

  call addfld ('PS&IC',horiz_only,'I','Pa', 'Surface pressure',gridname=grid_name)
  call addfld ('U&IC', (/'lev'/), 'I','m/s','Zonal wind',      gridname=grid_name)
  call addfld ('V&IC', (/'lev'/), 'I','m/s','Meridional wind', gridname=grid_name)
  call addfld ('T&IC', (/'lev'/), 'I','K',  'Temperature',     gridname=grid_name)
  do m = 1,pcnst
    call addfld (trim(cnst_name(m))//'&IC',(/'lev'/),'I','kg/kg',cnst_longname(m),gridname=grid_name)
  end do
  
  call add_default ('U&IC',0, 'I')
  call add_default ('V&IC',0, 'I')
  call add_default ('PS&IC      ',0, 'I')
  call add_default ('T&IC       ',0, 'I')
  do m = 1,pcnst
    call add_default(trim(cnst_name(m))//'&IC',0, 'I')
  end do

  call addfld('DYN_T'    ,(/ 'lev' /), 'A', 'K',    'Temperature (dyn grid)', gridname='GLL')
  call addfld('DYN_Q'    ,(/ 'lev' /), 'A', 'kg/kg','Water Vapor (dyn grid',  gridname='GLL' )
  call addfld('DYN_U'    ,(/ 'lev' /), 'A', 'm/s',  'Zonal Velocity',         gridname='GLL')
  call addfld('DYN_V'    ,(/ 'lev' /), 'A', 'm/s',  'Meridional Velocity',    gridname='GLL')
  call addfld('DYN_OMEGA',(/ 'lev' /), 'A', 'Pa/s', 'Vertical Velocity',      gridname='GLL' )
  call addfld('DYN_PS'   ,horiz_only,  'A', 'Pa',   'Surface pressure',       gridname='GLL')

end subroutine stepon_init

!-----------------------------------------------------------------------
!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_run1 -- Phase 1 of dynamics run method.
!
! !INTERFACE:
subroutine stepon_run1( dtime_out, phys_state, phys_tend,               &
                        pbuf2d, dyn_in, dyn_out )
  
  use dp_coupling, only: d_p_coupling
  use time_mod,    only: tstep      ! dynamics timestep
  use time_manager, only: is_last_step
  use control_mod, only: ftype
  use physics_buffer, only : physics_buffer_desc
  use hycoef,      only: hyam, hybm
  use se_single_column_mod, only: scm_setfield, scm_setinitial
  implicit none
!
! !OUTPUT PARAMETERS:
!

   real(r8), intent(out) :: dtime_out   ! Time-step
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend), intent(inout) :: phys_tend(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
   type (physics_buffer_desc), pointer :: pbuf2d(:,:)
   type (element_t), pointer :: elem(:)

!-----------------------------------------------------------------------

   elem => dyn_out%elem

   ! NOTE: dtime_out computed here must match formula below
   dtime_out = get_step_size()

   if(par%dynproc) then
      if(tstep <= 0)  call endrun( 'bad tstep')
      if(dtime_out <= 0)  call endrun( 'bad dtime')
   end if
   !----------------------------------------------------------
   ! Move data into phys_state structure.
   !----------------------------------------------------------
   
  ! Determine whether it is time for an IOP update;
  ! doiopupdate set to true if model time step > next available IOP 
  if (use_iop .and. .not. is_last_step()) then
    call setiopupdate
  end if
  
  if (single_column) then
    iop_update_phase1 = .true. 
    if (doiopupdate) call readiopdata( iop_update_phase1,hyam,hybm )
    call scm_setfield(elem,iop_update_phase1)       
  endif 
  
   call t_barrierf('sync_d_p_coupling', mpicom)
   call t_startf('d_p_coupling')
   call d_p_coupling (phys_state, phys_tend,  pbuf2d, dyn_out )
   call t_stopf('d_p_coupling') 
   
end subroutine stepon_run1

subroutine stepon_run2(phys_state, phys_tend, dyn_in, dyn_out )
   use bndry_mod,       only: bndry_exchangeV
   use dimensions_mod,  only: nlev, nlevp, nelemd, np, npsq
   use dyn_grid,        only: fv_nphys
   use dp_coupling,     only: p_d_coupling
   use parallel_mod,    only: par
   use dyn_comp,        only: TimeLevel, hvcoord
   use time_mod,        only: tstep, TimeLevel_Qdp   !  dynamics typestep
   use control_mod,     only: ftype, qsplit
   use hycoef,          only: hyai, hybi
   use cam_history,     only: outfld, hist_fld_active
   use prim_driver_base,only: applyCAMforcing_tracers
   use prim_advance_mod,only: applyCAMforcing_dynamics
   use element_ops,     only: get_temperature

   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
   real(r8) :: temperature(np,np,nlev)   ! Temperature from dynamics
   integer :: kptr, ie, ic, m, i, j, k, tl_f, tl_fQdp, velcomp
   real(r8) :: rec2dt
   real(r8) :: dp(np,np,nlev),fq,fq0,qn0, ftmp(npsq,nlev,2)
   real(r8) :: tmp_dyn(np,np,nlev,nelemd)
   real(r8) :: tmp_dyn_qf(np,np,nlev,nelemd) !(zhang73)
   real(r8) :: fmtmp(np,np,nlev)
   real(r8) :: p_m(np,np,nlev)    ! temporary midpoint pressure for DYN_OMEGA output
   real(r8) :: p_i(np,np,nlevp)   ! temporary interface pressure for DYN_OMEGA output
   real(r8) :: omega(np,np,nlev)  ! temporary omega for DYN_OMEGA output
   real(r8) :: dtime
   integer :: nlev_tot
   nlev_tot=(3+pcnst)*nlev


   dtime = get_step_size()


   ! copy from phys structures -> dynamics structures
   call t_barrierf('sync_p_d_coupling', mpicom)
   call t_startf('p_d_coupling')
   call p_d_coupling(phys_state, phys_tend,  dyn_in)
   call t_stopf('p_d_coupling')

   if(.not.par%dynproc) return

   call t_startf('stepon_bndry_exch')
   ! do boundary exchange
   if (.not. single_column) then 
      do ie=1,nelemd

         if (fv_nphys>0) then
            ! We need to apply mass matrix weighting when FV physics grid is used
            do k = 1,nlev
               dyn_in%elem(ie)%derived%FT(:,:,k)   = dyn_in%elem(ie)%derived%FT(:,:,k)   * dyn_in%elem(ie)%spheremp(:,:)
               do m = 1,2
                  dyn_in%elem(ie)%derived%FM(:,:,m,k) = dyn_in%elem(ie)%derived%FM(:,:,m,k) * dyn_in%elem(ie)%spheremp(:,:)
               end do
               do ic = 1,pcnst
                  dyn_in%elem(ie)%derived%FQ(:,:,k,ic) = dyn_in%elem(ie)%derived%FQ(:,:,k,ic) * dyn_in%elem(ie)%spheremp(:,:)
               end do
            end do ! k = 1, nlev
         end if ! fv_nphys>0

         kptr=0
         ! fmtmp can be removed if theta and preqx model had the same size FM array
         fmtmp=dyn_in%elem(ie)%derived%FM(:,:,1,:)
         call edgeVpack_nlyr(edge_g,dyn_in%elem(ie)%desc,fmtmp,nlev,kptr,nlev_tot)
         kptr=kptr+nlev
         
         fmtmp=dyn_in%elem(ie)%derived%FM(:,:,2,:)
         call edgeVpack_nlyr(edge_g,dyn_in%elem(ie)%desc,fmtmp,nlev,kptr,nlev_tot)
         kptr=kptr+nlev
         
         call edgeVpack_nlyr(edge_g,dyn_in%elem(ie)%desc,dyn_in%elem(ie)%derived%FT(:,:,:),nlev,kptr,nlev_tot)
         kptr=kptr+nlev
         call edgeVpack_nlyr(edge_g,dyn_in%elem(ie)%desc,dyn_in%elem(ie)%derived%FQ(:,:,:,:),nlev*pcnst,kptr,nlev_tot)

      end do ! ie
   endif ! .not. single_column

   call bndry_exchangeV(par, edge_g)

   ! NOTE: rec2dt MUST be 1/dtime_out as computed above

   rec2dt = 1._r8/dtime

   do ie=1,nelemd
      if (.not. single_column) then

         kptr=0

         fmtmp=dyn_in%elem(ie)%derived%FM(:,:,1,:)
         call edgeVunpack_nlyr(edge_g,dyn_in%elem(ie)%desc,fmtmp,nlev,kptr,nlev_tot)
         dyn_in%elem(ie)%derived%FM(:,:,1,:)=fmtmp
         kptr=kptr+nlev
         
         fmtmp=dyn_in%elem(ie)%derived%FM(:,:,2,:)
         call edgeVunpack_nlyr(edge_g,dyn_in%elem(ie)%desc,fmtmp,nlev,kptr,nlev_tot)
         dyn_in%elem(ie)%derived%FM(:,:,2,:)=fmtmp
         kptr=kptr+nlev
         
         call edgeVunpack_nlyr(edge_g,dyn_in%elem(ie)%desc,dyn_in%elem(ie)%derived%FT(:,:,:),nlev,kptr,nlev_tot)
         kptr=kptr+nlev
         
         call edgeVunpack_nlyr(edge_g,dyn_in%elem(ie)%desc,dyn_in%elem(ie)%derived%FQ(:,:,:,:),nlev*pcnst,kptr,nlev_tot)

         if (fv_nphys>0) then
            ! We need to apply inverse mass matrix weighting when FV physics grid is used
            do k = 1,nlev
               dyn_in%elem(ie)%derived%FT(:,:,k)   = dyn_in%elem(ie)%derived%FT(:,:,k)   * dyn_in%elem(ie)%rspheremp(:,:)
               do m = 1,2
                  dyn_in%elem(ie)%derived%FM(:,:,m,k) = dyn_in%elem(ie)%derived%FM(:,:,m,k) * dyn_in%elem(ie)%rspheremp(:,:)
               end do
               do ic = 1,pcnst
                  dyn_in%elem(ie)%derived%FQ(:,:,k,ic) = dyn_in%elem(ie)%derived%FQ(:,:,k,ic) * dyn_in%elem(ie)%rspheremp(:,:)
               end do
            end do ! k = 1, nlev
         end if ! fv_nphys>0
         
      endif ! .not. single_column

      tl_f = TimeLevel%n0   ! timelevel which was adjusted by physics

      call TimeLevel_Qdp(TimeLevel, qsplit, tl_fQdp)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ftype=2,4:  apply forcing to Q,ps.  Return dynamics tendencies
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( (ftype==2) .or. (ftype==4) ) then
         ! apply forcing to states tl_f 
         ! requires forward-in-time timestepping, checked in namelist_mod.F90

         call applyCAMforcing_tracers(dyn_in%elem(ie),hvcoord,tl_f,tl_fQdp,dtime,.true.)

      endif ! if ftype == 2 or == 4

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ftype=1:  apply all forcings as an adjustment
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ftype==1) then
         ! apply forcing to state tl_f
         ! requires forward-in-time timestepping, checked in namelist_mod.F90i

         !ftype1 also requires a call ty applycamforcing_dynamics, below
         call applyCAMforcing_tracers(dyn_in%elem(ie),hvcoord,tl_f,tl_fQdp,dtime,.true.)
      endif !ftype=1 

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ftype=0 and ftype<0 (debugging options):  just return tendencies to dynamics
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ftype<=0) then

         do ic=1,pcnst
            ! Q  =  data used for forcing, at timelevel nm1   t
            ! FQ =  adjusted Q returned by forcing,  at time  t+dt
            ! tendency = (FQ*dp - Q*dp) / dt 
            ! Convert this to a tendency on Qdp:  CAM physics did not change dp3d
!$omp parallel do private(k, j, i)
            do k=1,nlev
               do j=1,np
                  do i=1,np
                     dyn_in%elem(ie)%derived%FQ(i,j,k,ic)=&
                          (  dyn_in%elem(ie)%derived%FQ(i,j,k,ic) - &
                          dyn_in%elem(ie)%state%Q(i,j,k,ic))*rec2dt* &  
                          dyn_in%elem(ie)%state%dp3d(i,j,k,tl_f) 
                  end do
               end do
            end do
         end do
      endif ! ftype<0
   end do   ! ie loop

   if (ftype==1) then
      call applyCAMforcing_dynamics(dyn_in%elem,hvcoord,tl_f,dtime,1,nelemd)
   endif

   call t_stopf('stepon_bndry_exch')



   ! Most output is done by physics.  We pass to the physics state variables
   ! at timelevel "tl_f".  
   ! we will output dycore variables here to ensure they are always at the same
   ! time as what the physics is writing.  
   ! in single_column mode, dycore is not fully initialized so dont use
   ! outfld() on dycore decompositions
   if (.not. single_column) then 
   if (hist_fld_active('VOR')) then
      call compute_zeta_C0(tmp_dyn,dyn_in%elem,par,tl_f)
      do ie=1,nelemd
         call outfld('VOR',tmp_dyn(1,1,1,ie),npsq,ie)
      enddo
   endif
   ! output DIV_Qflux: 1) +return: tmp_dyn_qf, 2) +input: tl_fQdp for <compute_div_C0> (zhang73)
   if (hist_fld_active('DIV')) then
      !call compute_div_C0(tmp_dyn,dyn_in%elem,par,tl_f)
      call compute_div_C0(tmp_dyn,tmp_dyn_qf,dyn_in%elem,par,tl_f,tl_fQdp)
      do ie=1,nelemd
         call outfld('DIV',tmp_dyn(1,1,1,ie),npsq,ie)
      enddo
   endif
   if (hist_fld_active('DIV_Qflux')) then
      call compute_div_C0(tmp_dyn,tmp_dyn_qf,dyn_in%elem,par,tl_f,tl_fQdp)
      do ie=1,nelemd
         call outfld('DIV_Qflux',tmp_dyn_qf(1,1,1,ie),npsq,ie)
      enddo
   endif
   if (hist_fld_active('FU') .or. hist_fld_active('FV') ) then
      do ie=1,nelemd
         do j=1,np
            do i=1,np
               ftmp(i+(j-1)*np,:,1) = dyn_in%elem(ie)%derived%FM(i,j,1,:)
               ftmp(i+(j-1)*np,:,2) = dyn_in%elem(ie)%derived%FM(i,j,2,:)
            end do
         end do
         call outfld('FU',ftmp(:,:,1),npsq,ie)
         call outfld('FV',ftmp(:,:,2),npsq,ie)
      end do
   endif
   endif
   
   do ie = 1,nelemd
      call get_temperature(dyn_in%elem(ie),temperature,hvcoord,tl_f)
      call outfld('DYN_T'     ,temperature                            ,npsq,ie)
      call outfld('DYN_Q'     ,dyn_in%elem(ie)%state%Q(:,:,:,1)       ,npsq,ie)
      call outfld('DYN_U'     ,dyn_in%elem(ie)%state%V(:,:,1,:,tl_f)  ,npsq,ie)
      call outfld('DYN_V'     ,dyn_in%elem(ie)%state%V(:,:,2,:,tl_f)  ,npsq,ie)

      ! at this point in the code we are on floating lagrangian levels
      ! need to compute ps for output
      p_i(:,:,nlevp) = hvcoord%hyai(1)*hvcoord%ps0 + &
           sum(dyn_in%elem(ie)%state%dp3d(:,:,:,tl_f),3)
      call outfld('DYN_PS'    ,p_i(:,:,nlevp),npsq,ie)
#ifdef MODEL_THETA_L
      ! omega_p is just omega when using the theta dycore
      call outfld('DYN_OMEGA',dyn_in%elem(ie)%derived%omega_p(:,:,:)  ,npsq,ie)
#else
      ! Multiply omega_p by pressure to get omega for output
      p_i(:,:,1) = hvcoord%hyai(1)*hvcoord%ps0
      do k = 1,nlev
         p_i(:,:,k+1) = p_i(:,:,k) + dyn_in%elem(ie)%state%dp3d(:,:,k,tl_f)
         p_m(:,:,k)   = ( p_i(:,:,k+1) + p_i(:,:,k) )/2
         omega(:,:,k) = dyn_in%elem(ie)%derived%omega_p(:,:,k) * p_m(:,:,k)
      enddo
      call outfld('DYN_OMEGA',omega(:,:,:),npsq,ie)
#endif
   end do
   
   
   
   end subroutine stepon_run2
   

subroutine stepon_run3(dtime, cam_out, phys_state, dyn_in, dyn_out)
   use camsrfexch,  only: cam_out_t     
   use dyn_comp,    only: dyn_run
   use time_mod,    only: tstep
   use hycoef,      only: hyam, hybm
   use dimensions_mod, only: nlev, nelemd, np, npsq
   use se_single_column_mod, only: scm_setfield, scm_setinitial
   use dyn_comp, only: TimeLevel
   use cam_history,     only: outfld   
   use cam_logfile, only: iulog
   real(r8), intent(in) :: dtime   ! Time-step
   real(r8) :: ftmp_temp(np,np,nlev,nelemd), ftmp_q(np,np,nlev,pcnst,nelemd)
   real(r8) :: out_temp(npsq,nlev), out_q(npsq,nlev), out_u(npsq,nlev), &
               out_v(npsq,nlev), out_psv(npsq)  
   real(r8), parameter :: rad2deg = 180.0 / SHR_CONST_PI
   real(r8), parameter :: fac = 1000._r8	     
   type(cam_out_t),     intent(inout) :: cam_out(:) ! Output from CAM to surface
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
   type (element_t), pointer :: elem(:)
   integer :: rc, i, j, k, p, ie, tl_f
#if defined (E3SM_SCM_REPLAY)
   real(r8) :: forcing_temp(npsq,nlev), forcing_q(npsq,nlev,pcnst)
#endif   
   
   elem => dyn_out%elem
   
#if (defined E3SM_SCM_REPLAY)   

   tl_f = TimeLevel%n0   ! timelevel which was adjusted by physics
   
   ! Save ftmp stuff to get state before dynamics is called
   do ie=1,nelemd
     ftmp_temp(:,:,:,ie) = dyn_in%elem(ie)%state%T(:,:,:,tl_f)
     ftmp_q(:,:,:,:,ie) = dyn_in%elem(ie)%state%Q(:,:,:,:)
   enddo

#endif   
   
   if (single_column) then
     
     ! Update IOP properties e.g. omega, divT, divQ
     
     iop_update_phase1 = .false. 
     if (doiopupdate) then
       call scm_setinitial(elem)
       call readiopdata(iop_update_phase1,hyam,hybm)
       call scm_setfield(elem,iop_update_phase1)
     endif   

   endif   

   call t_barrierf('sync_dyn_run', mpicom)
   call t_startf ('dyn_run')
   call dyn_run(dyn_out,rc)	
   call t_stopf  ('dyn_run')
   
   ! Update to get tendency 
#if (defined E3SM_SCM_REPLAY) 

   tl_f = TimeLevel%n0
   
   do ie=1,nelemd
     do k=1,nlev
       do j=1,np
         do i=1,np
	 
           ! Note that this calculation will not provide b4b results with 
	   !  an E3SM because the dynamics tendency is not computed in the exact
	   !  same way as an E3SM run, introducing error with roundoff 
	   forcing_temp(i+(j-1)*np,k) = (dyn_in%elem(ie)%state%T(i,j,k,tl_f) - &
	        ftmp_temp(i,j,k,ie))/dtime - dyn_in%elem(ie)%derived%FT(i,j,k)	
           out_temp(i+(j-1)*np,k) = dyn_in%elem(ie)%state%T(i,j,k,tl_f)
	   out_u(i+(j-1)*np,k) = dyn_in%elem(ie)%state%v(i,j,1,k,tl_f)
	   out_v(i+(j-1)*np,k) = dyn_in%elem(ie)%state%v(i,j,2,k,tl_f)
	   out_q(i+(j-1)*np,k) = dyn_in%elem(ie)%state%Q(i,j,k,1)
	   out_psv(i+(j-1)*np) = dyn_in%elem(ie)%state%ps_v(i,j,tl_f)

	   do p=1,pcnst	 	
	     forcing_q(i+(j-1)*np,k,p) = (dyn_in%elem(ie)%state%Q(i,j,k,p) - &
	        ftmp_q(i,j,k,p,ie))/dtime
	   enddo
	   
	 enddo
       enddo
     enddo
     
     call outfld('Ps',out_psv,npsq,ie)
     call outfld('t',out_temp,npsq,ie)
     call outfld('q',out_q,npsq,ie)
     call outfld('u',out_u,npsq,ie)
     call outfld('v',out_v,npsq,ie)
     call outfld('divT3d',forcing_temp,npsq,ie)
     do p=1,pcnst
       call outfld(trim(cnst_name(p))//'_dten',forcing_q(:,:,p),npsq,ie)   
     enddo
     
   enddo

#endif   

end subroutine stepon_run3


!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_final --- Dynamics finalization
!
! !INTERFACE:
subroutine stepon_final(dyn_in, dyn_out)
  use dyn_grid,         only: fv_physgrid_final, fv_nphys
  use cam_logfile, only: iulog
  use prim_driver_base,only: prim_finalize
! !PARAMETERS:
  ! WARNING: intent(out) here means that pointers in dyn_in and dyn_out
  ! are nullified. Unless this memory is released in some other routine,
  ! this is a memory leak.
  ! These currently seem to point to dyn_grid global data.
  type (dyn_import_t), intent(out) :: dyn_in  ! Dynamics import container
  type (dyn_export_t), intent(out) :: dyn_out ! Dynamics export container
!
! !DESCRIPTION:
!
! Deallocate data needed for dynamics. Finalize any dynamics specific
! files or subroutines.
!
!EOP
!-----------------------------------------------------------------------
!BOC

   ! Deallocate variables needed for the FV physics grid
   if (fv_nphys > 0) then
     if (par%masterproc)  write(iulog,*) "stepon: phygrid finalization..."
     call fv_physgrid_final()
   end if ! fv_nphys > 0
   if (par%masterproc)  write(iulog,*) "stepon: HOMME finalization..."
   call prim_finalize
   if (par%masterproc)  write(iulog,*) "stepon: End of finalization"
!EOC
end subroutine stepon_final

!-----------------------------------------------------------------------


end module stepon
