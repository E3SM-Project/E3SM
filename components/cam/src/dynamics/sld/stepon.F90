module stepon
!-----------------------------------------------------------------------
!
! Purpose:
! Module for time-stepping of the Spectral Semi-Lagrangian Dynamical core.
!
!-----------------------------------------------------------------------
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use shr_sys_mod,      only: shr_sys_flush
  use pmgrid,           only: plev, plat, plevp, plon, beglat, endlat
  use spmd_utils,       only: masterproc
  use scanslt,          only: advection_state
  use prognostics,      only: u3, v3, t3, q3, n3, div, dpsl, dpsm, ps, omga, &
                              phis
  use camsrfexch,       only: cam_out_t     
  use ppgrid,           only: begchunk, endchunk
  use physics_types,    only: physics_state, physics_tend
  use time_manager,     only: is_first_step
  use perf_mod
  use dyn_comp,         only: dyn_import_t, dyn_export_t
  implicit none

  private   ! By default make all data and methods private to this module
!
! Public methods
!
  public stepon_init     ! Initialization
  public stepon_run1     ! Run method phase 1
  public stepon_run2     ! Run method phase 2
  public stepon_run3     ! Run method phase 3
  public stepon_final    ! Finalization
!
! Private module data
!
  save
  type(physics_state), pointer :: phys_state(:)   ! Physics state data
  type(physics_tend ), pointer :: phys_tend(:)    ! Physics tendency data

  real(r8) :: detam(plev)               ! intervals between vert full levs.
  real(r8) :: cwava(plat)               ! weight applied to global integrals
  real(r8), allocatable :: t2(:,:,:)    ! temp tendency
  real(r8), allocatable :: fu(:,:,:)    ! u wind tendency
  real(r8), allocatable :: fv(:,:,:)    ! v wind tendency
  real(r8), allocatable :: flx_net(:,:) ! net flux from physics
  real(r8) :: coslat(plon)              ! cosine of latitude 
  real(r8) :: rcoslat(plon)             ! Inverse of coseine of latitude
  real(r8) :: rpmid(plon,plev)          ! inverse of midpoint pressure
  real(r8) :: pdel(plon,plev)           ! Pressure depth of layer
  real(r8) :: pint(plon,plevp)          ! Pressure at interfaces
  real(r8) :: pmid(plon,plev)           ! Pressure at midpoint
  type(advection_state) :: adv_state    ! Advection state data

!=======================================================================
contains
!=======================================================================

subroutine stepon_init( gw, etamid , dyn_in, dyn_out)
!-----------------------------------------------------------------------
!
! Purpose:
! Stepon initialization, primarily dynamics initialization.
!
!-----------------------------------------------------------------------
   use rgrid,          only: nlon
   use commap,         only: clat
   use time_manager,   only: is_first_step
   use scanslt,        only: slt_initial
   use hycoef,         only: hyam, hybm
!-----------------------------------------------------------------------
! Output arguments
!
  real(r8), intent(out) :: gw(plat)                  ! Gaussian weights
  real(r8), intent(out) :: etamid(plev)              ! vertical coords at midpoints
  type(dyn_import_t) :: dyn_in                       ! included for compatibility
  type(dyn_export_t) :: dyn_out                      ! included for compatibility
	
! Local variables
!
   integer i,k,lat                    ! longitude,level,latitude indices
!-----------------------------------------------------------------------
   !
   ! Define eta coordinates on midpoints
   !
   call t_startf ('stepon_startup')

   do k=1,plev
      etamid(k) = hyam(k) + hybm(k)
   end do

   if (is_first_step()) then
      do lat=beglat,endlat
         do i=1,nlon(lat)
            coslat(i) = cos(clat(lat))
            rcoslat(i) = 1._r8/coslat(i)
         end do
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
         !
         ! Calculate vertical motion field
         !
         call omcalc (rcoslat, div(1,1,lat,n3), u3(1,1,lat,n3), v3(1,1,lat,n3), dpsl(1,lat), &
                      dpsm(1,lat), pmid, pdel, rpmid, pint(1,plevp), &
                      omga(1,1,lat), nlon(lat))
      end do
   end if

   call slt_initial( detam, gw, cwava, etamid, adv_state )

   allocate(t2(plon,plev,beglat:endlat))
   allocate(fu(plon,plev,beglat:endlat))
   allocate(fv(plon,plev,beglat:endlat))
   allocate(flx_net(plon,beglat:endlat))
   call t_stopf ('stepon_startup')

end subroutine stepon_init

!
!=======================================================================
!

subroutine stepon_run1( ztodt, phys_state, phys_tend, pbuf2d, dyn_in, dyn_out )
!-----------------------------------------------------------------------
!
! Purpose:	First phase of run method for SLD. Return
!		time-step to atmosphere and couple to physics.
!
!-----------------------------------------------------------------------
   use time_manager,   only: get_step_size
   use dp_coupling,    only: d_p_coupling, p_d_coupling
   
   use physics_buffer, only : physics_buffer_desc
!
! Output arguments:
!
   real(r8), intent(out) :: ztodt      ! time step
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   type(dyn_import_t) :: dyn_in                       ! included for compatibility
   type(dyn_export_t) :: dyn_out                      ! included for compatibility

   ztodt = get_step_size()

   !
   ! Dump state variables to IC file
   !
   call t_startf ('diag_dynvar_ic')
   call diag_dynvar_ic (phis, ps(:,:,n3), t3(:,:,:,n3), u3(:,:,:,n3), &
                        v3(:,:,:,n3), q3(:,:,:,:,n3) )
   call t_stopf ('diag_dynvar_ic')
   !
   !----------------------------------------------------------
   ! PHYS_RUN Call the Physics package
   !----------------------------------------------------------
   !
   call t_startf ('d_p_coupling')
   call d_p_coupling (ps(:,:,n3), t3(:,:,:,n3), u3(:,:,:,n3), &
                      v3(:,:,:,n3), q3(:,:,:,:,n3), &
                      omga, phis, phys_state, phys_tend,  pbuf2d)
   call t_stopf  ('d_p_coupling')
end subroutine stepon_run1

!
!=======================================================================
!

subroutine stepon_run2( phys_state, phys_tend, dyn_in, dyn_out )
!-----------------------------------------------------------------------
!
! Purpose:	Second phase of run method for SLD. Couple from physics
!		to dynamics.
!
!-----------------------------------------------------------------------
  use dp_coupling,    only: p_d_coupling
!
! Input arguments:
!
  type(physics_state), intent(in):: phys_state(begchunk:endchunk)
  type(physics_tend), intent(in):: phys_tend(begchunk:endchunk)
  type(dyn_import_t) :: dyn_in                       ! included for compatibility
  type(dyn_export_t) :: dyn_out                      ! included for compatibility

   call t_startf ('p_d_coupling')
   call p_d_coupling (phys_state, phys_tend, t2, fu, fv, flx_net,q3(:,:,:,:,n3))
   call t_stopf  ('p_d_coupling')
end subroutine stepon_run2

!
!=======================================================================
!

subroutine stepon_run3( ztodt, etamid, cam_out, phys_state, dyn_in, dyn_out )
!-----------------------------------------------------------------------
!
! Purpose:	Final phase of run method for SLD. Run the actual
!		dynamics.
!
!-----------------------------------------------------------------------
!
! Input arguments:
!
   real(r8), intent(in) :: ztodt            ! time step
   real(r8), intent(in) :: etamid(plev)     ! vertical coords at midpoints
   type(physics_state), intent(in):: phys_state(begchunk:endchunk)
  type(dyn_import_t) :: dyn_in                       ! included for compatibility
  type(dyn_export_t) :: dyn_out                      ! included for compatibility
!
! Input/output argument:
!
   type(cam_out_t), intent(inout) :: cam_out(begchunk:endchunk)
   !----------------------------------------------------------
   ! DYNPKG Call the Dynamics Package
   !----------------------------------------------------------

   call t_startf ('dynpkg')
   call dynpkg(t2      ,fu      ,fv      ,etamid  , &
                cwava   ,detam   ,flx_net , ztodt  ,adv_state)
   call t_stopf  ('dynpkg')
end subroutine stepon_run3

!
!=======================================================================
!

subroutine stepon_final(dyn_in, dyn_out)
!-----------------------------------------------------------------------
!
! Purpose:  Stepon finalization.
!
!-----------------------------------------------------------------------
   use scanslt, only: slt_final
   type(dyn_import_t) :: dyn_in                       ! included for compatibility
   type(dyn_export_t) :: dyn_out                      ! included for compatibility

   call slt_final( adv_state )
   deallocate(t2)
   deallocate(fu)
   deallocate(fv)
end subroutine stepon_final

!
!=======================================================================
!

end module stepon
