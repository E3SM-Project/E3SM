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
   use quadrature_mod, only: gauss, gausslobatto, quadrature_t
   use edgetype_mod,       only: EdgeBuffer_t
   use edge_mod,       only: initEdgeBuffer, FreeEdgeBuffer, edgeVpack, edgeVunpack
   use parallel_mod,   only : par
   use scamMod,        only: use_iop, doiopupdate, single_column, &
                             setiopupdate, readiopdata
   use element_mod,    only: element_t

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
!
!EOP
!----------------------------------------------------------------------
!BOC
!
! !PRIVATE DATA MEMBERS:
!

  logical :: iop_update_surface

  type (derivative_t)   :: deriv           ! derivative struct
  type (quadrature_t)   :: gv,gp           ! quadratures on velocity and pressure grids
  type (EdgeBuffer_t) :: edgebuf              ! edge buffer
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
  use dimensions_mod, only: nlev, nelemd, npsq
  use cam_history,    only: addfld, add_default, horiz_only
  use cam_history,    only: register_vector_field
  use control_mod,    only: smooth_phis_numcycle
  use gravity_waves_sources, only: gws_init
  use phys_control,   only: use_gw_front

! !OUTPUT PARAMETERS
!
  type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
  type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container

  integer :: m
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
     call initEdgeBuffer(par, edgebuf, dyn_in%elem, (3+pcnst)*nlev)
     if (use_gw_front)  call gws_init(dyn_in%elem)
  end if

  ! fields that are written by the dycore
  ! these calls cant be done in dyn_init() because physics grid
  ! is not initialized at that point if making a restart runs
  !
  ! Forcing from physics
  ! FU, FV, other dycores, doc, says "m/s" but I think that is m/s^2
  call addfld ('FU',  (/ 'lev' /), 'A', 'm/s2', 'Zonal wind forcing term',     gridname='GLL')
  call addfld ('FV',  (/ 'lev' /), 'A', 'm/s2', 'Meridional wind forcing term',gridname='GLL')
  call register_vector_field('FU', 'FV')
  call addfld ('VOR', (/ 'lev' /), 'A', '1/s',  'Vorticity',                   gridname='GLL')
  call addfld ('DIV', (/ 'lev' /), 'A', '1/s',  'Divergence',                  gridname='GLL')

  if (smooth_phis_numcycle>0) then
     call addfld ('PHIS_SM',  horiz_only, 'I', 'm2/s2', 'Surface geopotential (smoothed)',                gridname='GLL')
     call addfld ('SGH_SM',   horiz_only, 'I', 'm',     'Standard deviation of orography (smoothed)',     gridname='GLL')
     call addfld ('SGH30_SM', horiz_only, 'I', 'm',     'Standard deviation of 30s orography (smoothed)', gridname='GLL')
  endif

  call addfld ('CONVU   ', (/ 'ilev' /),'A', 'm/s2    ','Zonal component IE->KE conversion term',      gridname='physgrid')
  call addfld ('CONVV   ', (/ 'ilev' /),'A', 'm/s2    ','Meridional component IE->KE conversion term', gridname='physgrid')
  call register_vector_field('CONVU', 'CONVV')
  call addfld ('DIFFU   ', (/ 'ilev' /),'A', 'm/s2    ','U horizontal diffusion',                      gridname='physgrid')
  call addfld ('DIFFV   ', (/ 'ilev' /),'A', 'm/s2    ','V horizontal diffusion',                      gridname='physgrid')
  call register_vector_field('DIFFU', 'DIFFV')
  
  call addfld ('ETADOT', (/ 'ilev' /), 'A', '1/s', 'Vertical (eta) velocity', gridname='physgrid')
  call addfld ('U&IC',   (/ 'lev' /),  'I', 'm/s', 'Zonal wind',              gridname='physgrid' )
  call addfld ('V&IC',   (/ 'lev' /),  'I', 'm/s', 'Meridional wind',         gridname='physgrid' )
  ! Don't need to register U&IC V&IC since we don't interpolate IC files
  call add_default ('U&IC',0, 'I')
  call add_default ('V&IC',0, 'I')

  call addfld ('PS&IC', horiz_only,  'I', 'Pa', 'Surface pressure',gridname='physgrid')
  call addfld ('T&IC',  (/ 'lev' /), 'I', 'K',  'Temperature',     gridname='physgrid')

  call add_default ('PS&IC      ',0, 'I')
  call add_default ('T&IC       ',0, 'I')
  do m = 1,pcnst
     call addfld (trim(cnst_name(m))//'&IC', (/ 'lev' /), 'I', 'kg/kg', cnst_longname(m), gridname='physgrid')
  end do
  do m = 1,pcnst
     call add_default(trim(cnst_name(m))//'&IC',0, 'I')
  end do


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
  use time_mod,    only: tstep, phys_tscale       ! dynamics timestep
  use time_manager, only: is_last_step
  use control_mod, only: ftype
  use physics_buffer, only : physics_buffer_desc
  use hycoef,      only: hyam, hybm
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

!-----------------------------------------------------------------------

   ! NOTE: dtime_out computed here must match formula below
   dtime_out = get_step_size()
   if (phys_tscale/=0) then
      dtime_out=phys_tscale  ! set by user in namelist
   endif

   if(par%dynproc) then
      if(tstep <= 0)  call endrun( 'bad tstep')
      if(dtime_out <= 0)  call endrun( 'bad dtime')
   end if
   !----------------------------------------------------------
   ! Move data into phys_state structure.
   !----------------------------------------------------------
   
   call t_barrierf('sync_d_p_coupling', mpicom)
   call t_startf('d_p_coupling')
   call d_p_coupling (phys_state, phys_tend,  pbuf2d, dyn_out )
   call t_stopf('d_p_coupling')
   
  ! Determine whether it is time for an IOP update;
  ! doiopupdate set to true if model time step > next available IOP 
  if (use_iop .and. .not. is_last_step()) then
    call setiopupdate
  end if
  
  if (single_column) then
    iop_update_surface = .true. 
    if (doiopupdate) call readiopdata( iop_update_surface,hyam,hybm )
  endif   
   
end subroutine stepon_run1

subroutine stepon_run2(phys_state, phys_tend, dyn_in, dyn_out )
   use bndry_mod,      only: bndry_exchangeV
   use dimensions_mod, only: nlev, nelemd, np, npsq
   use dp_coupling,    only: p_d_coupling
   use parallel_mod,   only: par
   use dyn_comp,       only: TimeLevel
   
   use time_mod,        only: tstep, phys_tscale, TimeLevel_Qdp   !  dynamics typestep
   use control_mod,     only: ftype, qsplit, smooth_phis_numcycle
   use hycoef,          only: hyai, hybi, ps0
   use cam_history,     only: outfld, hist_fld_active
   use nctopo_util_mod, only: phisdyn,sghdyn,sgh30dyn


   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend), intent(inout) :: phys_tend(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
   integer :: kptr, ie, ic, i, j, k, tl_f, tl_fQdp
   real(r8) :: rec2dt, dyn_ps0
   real(r8) :: dp(np,np,nlev),dp_tmp,fq,fq0,qn0, ftmp(npsq,nlev,2)
   real(r8) :: dtime

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
       kptr=0
       call edgeVpack(edgebuf,dyn_in%elem(ie)%derived%FM(:,:,:,:),2*nlev,kptr,ie)
       kptr=kptr+2*nlev

       call edgeVpack(edgebuf,dyn_in%elem(ie)%derived%FT(:,:,:),nlev,kptr,ie)
       kptr=kptr+nlev
       call edgeVpack(edgebuf,dyn_in%elem(ie)%derived%FQ(:,:,:,:),nlev*pcnst,kptr,ie)
     end do
   endif

   call bndry_exchangeV(par, edgebuf)

   ! NOTE: rec2dt MUST be 1/dtime_out as computed above

   rec2dt = 1._r8/dtime
   if (phys_tscale/=0) then
      rec2dt = 1._r8/phys_tscale
   endif


   do ie=1,nelemd
     if (.not. single_column) then
       kptr=0

       call edgeVunpack(edgebuf,dyn_in%elem(ie)%derived%FM(:,:,:,:),2*nlev,kptr,ie)
       kptr=kptr+2*nlev

       call edgeVunpack(edgebuf,dyn_in%elem(ie)%derived%FT(:,:,:),nlev,kptr,ie)
       kptr=kptr+nlev

       call edgeVunpack(edgebuf,dyn_in%elem(ie)%derived%FQ(:,:,:,:),nlev*pcnst,kptr,ie)
     endif

      tl_f = TimeLevel%n0   ! timelevel which was adjusted by physics

      call TimeLevel_Qdp(TimeLevel, qsplit, tl_fQdp)

      dyn_ps0=ps0

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ftype=2,3,4:  apply forcing to Q,ps.  Return dynamics tendencies
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( (ftype==2) .or. (ftype==3) .or. (ftype==4) ) then
         ! apply forcing to states tl_f 
         ! requires forward-in-time timestepping, checked in namelist_mod.F90
!$omp parallel do private(k)
         do k=1,nlev
            dp(:,:,k) = ( hyai(k+1) - hyai(k) )*dyn_ps0 + &
                 ( hybi(k+1) - hybi(k) )*dyn_in%elem(ie)%state%ps_v(:,:,tl_f)
         enddo

         if (ftype == 3) then ! ftype == 3, scale tendencies with current dp
           do k=1,nlev
             do j=1,np
               do i=1,np
                  dyn_in%elem(ie)%derived%FT(i,j,k) = dyn_in%elem(ie)%derived%FT(i,j,k)*dp(i,j,k)
                  dyn_in%elem(ie)%derived%FM(i,j,1:2,k) = dyn_in%elem(ie)%derived%FM(i,j,1:2,k)*dp(i,j,k)
               end do
             end do
           end do
         endif !ftype 3

         do k=1,nlev
            do j=1,np
               do i=1,np
                  do ic=1,pcnst
                     ! apply forcing to Qdp
                     ! dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_fQdp) = &
                     !        dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_fQdp) + fq 
                     dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_fQdp) = &
                          dp(i,j,k)*dyn_in%elem(ie)%derived%FQ(i,j,k,ic)

! BEWARE critical region if using OpenMP over k (AAM)
                     if (ic==1) then
                        fq = dp(i,j,k)*(  dyn_in%elem(ie)%derived%FQ(i,j,k,ic) - &
                             dyn_in%elem(ie)%state%Q(i,j,k,ic))
                        ! force ps_v to conserve mass:  
                        dyn_in%elem(ie)%state%ps_v(i,j,tl_f)= &
                             dyn_in%elem(ie)%state%ps_v(i,j,tl_f) + fq
                     endif
                  enddo
               end do
            end do
         end do

!$omp parallel do private(k, j, i, ic, dp_tmp)
         do k=1,nlev
          do ic=1,pcnst
            do j=1,np
               do i=1,np
                  ! make Q consistent now that we have updated ps_v above
                  ! recompute dp, since ps_v was changed above
                  dp_tmp = ( hyai(k+1) - hyai(k) )*dyn_ps0 + &
                       ( hybi(k+1) - hybi(k) )*dyn_in%elem(ie)%state%ps_v(i,j,tl_f)
                  dyn_in%elem(ie)%state%Q(i,j,k,ic)= &
                       dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_fQdp)/dp_tmp
               end do
            end do
          end do
         end do

      endif ! if ftype == 2 or == 3 or == 4

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ftype=1:  apply all forcings as an adjustment
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ftype==1) then
         ! apply forcing to state tl_f
         ! requires forward-in-time timestepping, checked in namelist_mod.F90
!$omp parallel do private(k)
         do k=1,nlev
            dp(:,:,k) = ( hyai(k+1) - hyai(k) )*dyn_ps0 + &
                 ( hybi(k+1) - hybi(k) )*dyn_in%elem(ie)%state%ps_v(:,:,tl_f)
         enddo
         do k=1,nlev
            do j=1,np
               do i=1,np

                  do ic=1,pcnst
                     ! back out tendency: Qdp*dtime 
                     fq = dp(i,j,k)*(  dyn_in%elem(ie)%derived%FQ(i,j,k,ic) - &
                          dyn_in%elem(ie)%state%Q(i,j,k,ic))
                     
                     ! apply forcing to Qdp
                     dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_fQdp) = &
                          dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_fQdp) + fq 

                     ! BEWARE critical region if using OpenMP over k (AAM)
                     if (ic==1) then
                        ! force ps_v to conserve mass:  
                        dyn_in%elem(ie)%state%ps_v(i,j,tl_f)= &
                             dyn_in%elem(ie)%state%ps_v(i,j,tl_f) + fq
                     endif
                  enddo

                  ! force V, T, both timelevels
                  dyn_in%elem(ie)%state%v(i,j,:,k,tl_f)= &
                       dyn_in%elem(ie)%state%v(i,j,:,k,tl_f) +  &
                       dtime*dyn_in%elem(ie)%derived%FM(i,j,:,k)
                  
                  dyn_in%elem(ie)%state%T(i,j,k,tl_f)= &
                       dyn_in%elem(ie)%state%T(i,j,k,tl_f) + &
                       dtime*dyn_in%elem(ie)%derived%FT(i,j,k)
                  
               end do
            end do
         end do

!$omp parallel do private(k, j, i, ic, dp_tmp)
         do k=1,nlev
            do ic=1,pcnst
               do j=1,np
                  do i=1,np
                     ! make Q consistent now that we have updated ps_v above
                     dp_tmp = ( hyai(k+1) - hyai(k) )*dyn_ps0 + &
                          ( hybi(k+1) - hybi(k) )*dyn_in%elem(ie)%state%ps_v(i,j,tl_f)
                     dyn_in%elem(ie)%state%Q(i,j,k,ic)= &
                          dyn_in%elem(ie)%state%Qdp(i,j,k,ic,tl_fQdp)/dp_tmp
                  end do
               end do
            end do
         end do
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ftype=0 and ftype<0 (debugging options):  just return tendencies to dynamics
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ftype<=0) then

         do ic=1,pcnst
            ! Q  =  data used for forcing, at timelevel nm1   t
            ! FQ =  adjusted Q returned by forcing,  at time  t+dt
            ! tendency = (FQ*dp - Q*dp) / dt 
            ! Convert this to a tendency on Qdp:  CAM physics does not change ps
            ! so use ps_v at t.  (or, if physics worked with Qdp
            ! and returned FQdp, tendency = (FQdp-Qdp)/2dt and since physics
            ! did not change dp, dp would be the same in both terms)
!$omp parallel do private(k, j, i, dp_tmp)
            do k=1,nlev
               do j=1,np
                  do i=1,np
                     dp_tmp = ( hyai(k+1) - hyai(k) )*dyn_ps0 + &
                          ( hybi(k+1) - hybi(k) )*dyn_in%elem(ie)%state%ps_v(i,j,tl_f)
                     
                     dyn_in%elem(ie)%derived%FQ(i,j,k,ic)=&
                          (  dyn_in%elem(ie)%derived%FQ(i,j,k,ic) - &
                          dyn_in%elem(ie)%state%Q(i,j,k,ic))*rec2dt*dp_tmp
                  end do
               end do
            end do
         end do
      endif
   end do
   call t_stopf('stepon_bndry_exch')



   ! Most output is done by physics.  We pass to the physics state variables
   ! at timelevel "tl_f".  
   ! we will output dycore variables here to ensure they are always at the same
   ! time as what the physics is writing.  
#if 0
   if (hist_fld_active('VOR')) then
      call compute_zeta_C0(tmp_dyn,elem,hybrid,1,nelemd,tl_f,k)
      do ie=1,nelemd
         do j=1,np
            do i=1,np
               ftmp(i+(j-1)*np,1:pver,1) = tmp_dyn(i,j,1:pver)
            end do
         end do
         if (.not. single_column) call outfld('VOR',ftmp(:,:,1),npsq,ie)
      enddo
   endif
   if (hist_fld_active('DIV')) then
      call compute_div_C0(tmp_dyn,elem,hybrid,1,nelemd,tl_f,k)
      do ie=1,nelemd
         do j=1,np
            do i=1,np
               ftmp(i+(j-1)*np,1:pver,1) = tmp_dyn(i,j,1:pver)
            end do
         end do
         if (.not. single_column) call outfld('DIV',ftmp(:,:,1),npsq,ie)
      enddo
   endif
#endif
   if (smooth_phis_numcycle>0) then
      if (hist_fld_active('PHIS_SM')) then
         do ie=1,nelemd
            do j=1,np
               do i=1,np
                  ftmp(i+(j-1)*np,1,1) = phisdyn(i,j,ie)
               end do
            end do
            if (.not. single_column) call outfld('PHIS_SM',ftmp(:,1,1),npsq,ie)
         enddo
      endif
      if (hist_fld_active('SGH_SM')) then
         do ie=1,nelemd
            do j=1,np
               do i=1,np
                  ftmp(i+(j-1)*np,1,1) = sghdyn(i,j,ie)
               end do
            end do
            if (.not. single_column) call outfld('SGH_SM',ftmp(:,1,1),npsq,ie)
         enddo
      endif
      if (hist_fld_active('SGH30_SM')) then
         do ie=1,nelemd
            do j=1,np
               do i=1,np
                  ftmp(i+(j-1)*np,1,1) = sgh30dyn(i,j,ie)
               end do
            end do
            if (.not. single_column) call outfld('SGH30_SM',ftmp(:,1,1),npsq,ie)
         enddo
      endif
   end if
   
   if (hist_fld_active('FU') .or. hist_fld_active('FV') .and. .not. single_column) then
      do ie=1,nelemd
         do k=1,nlev
            do j=1,np
               do i=1,np
                  ftmp(i+(j-1)*np,k,1) = dyn_in%elem(ie)%derived%FM(i,j,1,k)
                  ftmp(i+(j-1)*np,k,2) = dyn_in%elem(ie)%derived%FM(i,j,2,k)
               end do
            end do
         end do
         
         call outfld('FU',ftmp(:,:,1),npsq,ie)
         call outfld('FV',ftmp(:,:,2),npsq,ie)
      end do
   endif
   
   
   
   
   end subroutine stepon_run2
   

subroutine stepon_run3(dtime, cam_out, phys_state, dyn_in, dyn_out)
   use camsrfexch,  only: cam_out_t     
   use dyn_comp,    only: dyn_run
   use time_mod,    only: tstep
   use hycoef,      only: hyam, hybm
   use se_single_column_mod, only: scm_setfield, scm_setinitial
   real(r8), intent(in) :: dtime   ! Time-step
   type(cam_out_t),     intent(inout) :: cam_out(:) ! Output from CAM to surface
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
   type (element_t), pointer :: elem(:)
   integer :: rc
   
   elem => dyn_out%elem
   
   if (single_column) then
     
     ! Update IOP properties e.g. omega, divT, divQ
     
     iop_update_surface = .false. 
     if (doiopupdate) then
       call scm_setinitial(elem)
       call readiopdata(iop_update_surface,hyam,hybm)
       call scm_setfield(elem)
     endif   

   endif   

   call t_barrierf('sync_dyn_run', mpicom)
   call t_startf ('dyn_run')
   call dyn_run(dyn_out,rc)	
   call t_stopf  ('dyn_run')

end subroutine stepon_run3


!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  stepon_final --- Dynamics finalization
!
! !INTERFACE:
subroutine stepon_final(dyn_in, dyn_out)

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


!EOC
end subroutine stepon_final

!-----------------------------------------------------------------------


end module stepon
