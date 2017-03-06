module stepon
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Module for time-stepping of the CAM Eulerian Spectral dynamics.
! 
!-----------------------------------------------------------------------
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use shr_sys_mod,      only: shr_sys_flush
  use pmgrid,           only: plev, plat, plevp, plon, beglat, endlat
  use spmd_utils,       only: masterproc
  use scanslt,          only: advection_state
  use prognostics,      only: ps, u3, v3, t3, q3, qminus, div, &
                              dpsl, dpsm, omga, phis, n3, n3m2, n3m1
  use camsrfexch,       only: cam_out_t     
  use ppgrid,           only: begchunk, endchunk
  use physics_types,    only: physics_state, physics_tend
  use time_manager,     only: is_first_step, get_step_size
  use iop,              only: setiopupdate, readiopdata
  use scamMod,          only: use_iop,doiopupdate,use_pert_frc,wfld,wfldh,single_column
  use perf_mod

  implicit none
  private
  save

  public stepon_init     ! Initialization
  public stepon_run1     ! Run method phase 1
  public stepon_run2     ! Run method phase 2
  public stepon_run3     ! Run method phase 3
  public stepon_final    ! Finalization
!
! Private module data
!
  type(physics_state), pointer :: phys_state(:)   ! Physics state data
  type(physics_tend ), pointer :: phys_tend(:)    ! Physics tendency data

  real(r8) :: detam(plev)               ! intervals between vert full levs.
  real(r8) :: cwava(plat)               ! weight applied to global integrals
  real(r8), allocatable :: t2(:,:,:)    ! temp tendency
  real(r8), allocatable :: fu(:,:,:)    ! u wind tendency
  real(r8), allocatable :: fv(:,:,:)    ! v wind tendency
  real(r8), allocatable :: flx_net(:,:) ! net flux from physics
  real(r8), allocatable :: fq(:,:,:,:)  ! Q tendencies,for eul_nsplit>1
  real(r8), allocatable :: t2_save(:,:,:)    ! temp tendency
  real(r8), allocatable :: fu_save(:,:,:)    ! u wind tendency
  real(r8), allocatable :: fv_save(:,:,:)    ! v wind tendency
  real(r8) :: coslat(plon)              ! cosine of latitude
  real(r8) :: rcoslat(plon)             ! Inverse of coseine of latitude
  real(r8) :: rpmid(plon,plev)          ! inverse of midpoint pressure
  real(r8) :: pdel(plon,plev)           ! Pressure depth of layer
  real(r8) :: pint(plon,plevp)          ! Pressure at interfaces
  real(r8) :: pmid(plon,plev)           ! Pressure at midpoint
  type(advection_state) :: adv_state    ! Advection state data

  real(r8) :: etamid(plev)              ! vertical coords at midpoints or pmid if single_column

!======================================================================= 
contains
!======================================================================= 

subroutine stepon_init(dyn_in, dyn_out)
!----------------------------------------------------------------------- 
! 
! Purpose:  Initialization, primarily of dynamics.
!
!----------------------------------------------------------------------- 
   use dyn_comp,       only: dyn_import_t, dyn_export_t
   use scanslt,        only: scanslt_initial
   use commap,         only: clat
   use constituents,   only: pcnst
   use physconst,      only: gravit
   use rgrid,          only: nlon
   use eul_control_mod,only: eul_nsplit
#if ( defined BFB_CAM_SCAM_IOP )
   use iop,            only:init_iop_fields
#endif
!-----------------------------------------------------------------------
! Arguments
!
  type(dyn_import_t) :: dyn_in                       ! included for compatibility
  type(dyn_export_t) :: dyn_out                      ! included for compatibility
!-----------------------------------------------------------------------
!  Local variables
!
   integer :: k, lat, i
   !-----------------------------------------------------------------------

   call t_startf ('stepon_startup')

   call scanslt_initial(adv_state, etamid, gravit, detam, cwava)
   !
   ! Initial guess for trajectory midpoints in spherical coords.
   ! nstep = 0:  use arrival points as initial guess for trajectory midpoints.
   ! nstep > 0:  use calculated trajectory midpoints from previous time 
   ! step as first guess.
   ! NOTE:  reduce number of iters necessary for convergence after nstep = 1.
   !
   if (is_first_step()) then
      do lat=beglat,endlat
	 if (.not. single_column) then
            do i=1,nlon(lat)
               coslat(i) = cos(clat(lat))
               rcoslat(i) = 1._r8/coslat(i)
            end do
         endif
         !     
         ! Set current time pressure arrays for model levels etc.
         !
         call plevs0(nlon(lat), plon, plev, ps(1,lat,n3), pint, pmid, pdel)
         !
         do k=1,plev
            do i=1,nlon(lat)
               rpmid(i,k) = 1._r8/pmid(i,k)
            end do
         end do

	 if (.not. single_column) then
         !
         ! Calculate vertical motion field
         !
            call omcalc (rcoslat, div(1,1,lat,n3), u3(1,1,lat,n3), v3(1,1,lat,n3), dpsl(1,lat), &
                      dpsm(1,lat), pmid, pdel, rpmid   ,pint(1,plevp), &
                      omga(1,1,lat), nlon(lat))
         else
         
            omga(1,:,lat)=wfld(:)
         endif
      end do
   end if

   allocate(t2(plon,plev,beglat:endlat))
   allocate(fu(plon,plev,beglat:endlat))
   allocate(fv(plon,plev,beglat:endlat))
   allocate( flx_net(plon,beglat:endlat))
   if (eul_nsplit>1) then
      allocate(fq(plon,plev,pcnst,beglat:endlat))
      allocate(t2_save(plon,plev,beglat:endlat))
      allocate(fu_save(plon,plev,beglat:endlat))
      allocate(fv_save(plon,plev,beglat:endlat))
   endif
   !
   ! Beginning of basic time step loop
   !
   call t_stopf ('stepon_startup')


#if ( defined BFB_CAM_SCAM_IOP )
   if (is_first_step()) then
      call init_iop_fields()
   endif
#endif
end subroutine stepon_init

!
!======================================================================= 
!

subroutine stepon_run1( ztodt, phys_state, phys_tend , pbuf2d, dyn_in, dyn_out)
!----------------------------------------------------------------------- 
! 
! Purpose:  Phase 1 run method of dynamics. Set the time-step
!           to use for physics. And couple from dynamics to physics.
!
!----------------------------------------------------------------------- 
  use dyn_comp,       only: dyn_import_t, dyn_export_t
  use time_manager,   only: get_nstep
  use prognostics,    only: pdeld
  
  use dp_coupling,    only: d_p_coupling
  use eul_control_mod,only: eul_nsplit
  use physics_buffer, only : physics_buffer_desc
  real(r8), intent(out) :: ztodt            ! twice time step unless nstep=0
  type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
  type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  type(dyn_import_t) :: dyn_in                       ! included for compatibility
  type(dyn_export_t) :: dyn_out                      ! included for compatibility

  real(r8) :: dtime                     ! timestep size
  !----------------------------------------------------------------------- 

  dtime = get_step_size()

  ztodt = 2.0_r8*dtime

  ! If initial time step adjust dt
  if (is_first_step()) ztodt = dtime

  ! subcycling case, physics dt is always dtime
  if (eul_nsplit>1) ztodt = dtime	

  ! Dump state variables to IC file
  call t_startf ('diag_dynvar_ic')
  call diag_dynvar_ic (phis, ps(:,beglat:endlat,n3m1), t3(:,:,beglat:endlat,n3m1), u3(:,:,beglat:endlat,n3m1), &
                       v3(:,:,beglat:endlat,n3m1), q3(:,:,:,beglat:endlat,n3m1) )
  call t_stopf ('diag_dynvar_ic')
  !
  !----------------------------------------------------------
  ! Couple from dynamics to physics
  !----------------------------------------------------------
  !
  call t_startf ('d_p_coupling')
  call d_p_coupling (ps(:,beglat:endlat,n3m2), t3(:,:,beglat:endlat,n3m2), u3(:,:,beglat:endlat,n3m2), &
                     v3(:,:,beglat:endlat,n3m2), q3(:,:,:,beglat:endlat,n3m2), &
                     omga, phis, phys_state, phys_tend,  pbuf2d, pdeld(:,:,:,n3m2))
  call t_stopf  ('d_p_coupling')
end subroutine stepon_run1

!
!======================================================================= 
!

subroutine stepon_run2( phys_state, phys_tend, dyn_in, dyn_out )
!----------------------------------------------------------------------- 
! 
! Purpose:  Phase 2 run method of dynamics. Couple from physics
!           to dynamics.
!
!----------------------------------------------------------------------- 
  use dyn_comp,       only: dyn_import_t, dyn_export_t
  use dp_coupling,    only: p_d_coupling
  type(physics_state), intent(in):: phys_state(begchunk:endchunk)
  type(physics_tend), intent(in):: phys_tend(begchunk:endchunk)
  type(dyn_import_t) :: dyn_in                       ! included for compatibility
  type(dyn_export_t) :: dyn_out                      ! included for compatibility

  call t_startf ('p_d_coupling')
  call p_d_coupling (phys_state, phys_tend, t2, fu, fv, flx_net, &
                     qminus )
  call t_stopf  ('p_d_coupling')
end subroutine stepon_run2

!
!======================================================================= 
!

subroutine stepon_run3( ztodt, cam_out, phys_state, dyn_in, dyn_out )
!----------------------------------------------------------------------- 
! 
! Purpose:  Final phase of dynamics run method. Run the actual dynamics.
!
!----------------------------------------------------------------------- 
  use dyn_comp,       only: dyn_import_t, dyn_export_t
  use eul_control_mod,only: eul_nsplit
  real(r8), intent(in) :: ztodt            ! twice time step unless nstep=0
  type(cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)
  type(physics_state), intent(in):: phys_state(begchunk:endchunk)
  type(dyn_import_t) :: dyn_in                       ! included for compatibility
  type(dyn_export_t) :: dyn_out                      ! included for compatibility
  real(r8) :: dt_dyn0,dt_dyn 
  integer :: stage
  if (single_column) then
     
     ! Determine whether it is time for an IOP update;
     ! doiopupdate set to true if model time step > next available IOP
     if (use_iop) then
        call setiopupdate
     end if
     
     ! Update IOP properties e.g. omega, divT, divQ
     
     if (doiopupdate) call readiopdata()
     
  endif

  !----------------------------------------------------------
  ! DYNPKG Call the Dynamics Package
  !----------------------------------------------------------
  call t_startf ('dynpkg')

  if (eul_nsplit==1) then	
     call dynpkg(adv_state, t2      ,fu      ,fv      ,etamid  ,  &
       cwava   ,detam   ,flx_net ,ztodt)
  else
     dt_dyn0 = ztodt/eul_nsplit
     dt_dyn = dt_dyn0
     if (is_first_step()) dt_dyn = 2*dt_dyn0

     ! convert q adjustment to a tendency
     fq = (qminus(:,:,:,:) - q3(:,:,:,:,n3m2))/ztodt
     ! save a copy of t2,fu,fv
     t2_save=t2		
     fu_save=fu
     fv_save=fv

     call apply_fq(qminus,q3(:,:,:,:,n3m2),fq,dt_dyn0)
     call dynpkg(adv_state, t2      ,fu      ,fv      ,etamid  ,  &
       cwava   ,detam   ,flx_net ,dt_dyn0)

     do stage=2,eul_nsplit	
        t2=t2_save		
        fu=fu_save
        fv=fv_save
        call apply_fq(qminus,q3(:,:,:,:,n3m2),fq,dt_dyn)
        call dynpkg(adv_state, t2      ,fu      ,fv      ,etamid  ,  &
             cwava   ,detam   ,flx_net ,dt_dyn)
     enddo
  endif

  call t_stopf  ('dynpkg')
end subroutine stepon_run3



subroutine apply_fq(qminus,q3,fq,dt)
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plat, plev, plevp, beglat, endlat
   use rgrid,        only: nlon
   use constituents,   only: pcnst

   real(r8), intent(in) :: q3(plon,plev,beglat:endlat,pcnst)
   real(r8), intent(in) :: fq(plon,plev,beglat:endlat,pcnst)
   real(r8), intent(out) :: qminus(plon,plev,beglat:endlat,pcnst)
   real(r8), intent(in) :: dt

   !local 	
   real(r8) :: q_tmp,fq_tmp
   integer :: q,c,k,i

   do q=1,pcnst 
   do c=beglat,endlat
   do k=1,plev
   do i=1,nlon(c)
      fq_tmp = dt*fq(i,k,c,q)
      q_tmp  = q3(i,k,c,q)
      ! if forcing is > 0, do nothing (it makes q less negative)
      if (fq_tmp<0 .and. q_tmp+fq_tmp<0 ) then
         ! reduce magnitude of forcing so it wont drive q negative 
         ! but we only reduce the magnitude of the forcing, dont increase
         ! its magnitude or change the sign
         
         ! if q<=0, then this will set fq=0  (q already negative)
         ! if q>0, then we know from above that fq < -q < 0, so we 
         ! can reduce the magnitive of fq by setting fq = -q:
         fq_tmp = min(-q_tmp,0._r8)
      endif
      qminus(i,k,c,q) = q_tmp + fq_tmp
   enddo
   enddo
   enddo
   enddo
   	
end subroutine


!
!======================================================================= 
!

subroutine stepon_final(dyn_in, dyn_out)
!----------------------------------------------------------------------- 
! 
! Purpose:  Stepon finalization.
!
!----------------------------------------------------------------------- 
   use dyn_comp,       only: dyn_import_t, dyn_export_t
   use scanslt, only: scanslt_final
   type(dyn_import_t) :: dyn_in                       ! included for compatibility
   type(dyn_export_t) :: dyn_out                      ! included for compatibility

   call scanslt_final( adv_state )
   deallocate(t2)
   deallocate(fu)
   deallocate(fv)
   deallocate(flx_net)

end subroutine stepon_final
!
!======================================================================= 
!

end module stepon
