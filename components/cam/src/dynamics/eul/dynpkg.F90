
subroutine dynpkg (adv_state, t2      ,fu      ,fv      ,etamid  ,          &
                   cwava   ,detam   ,flx_net ,ztodt   )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Driving routines for dynamics and transport.
! 
! Method: 
! 
! Author: 
! Original version:  CCM3
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plat, plev, plevp, beglat, endlat
   use pspect
   use comspe
   use scanslt,      only: scanslt_run, plond, platd, advection_state
   use scan2,        only: scan2run
   use scamMod,      only: single_column,scm_crm_mode,switch,wfldh
#if ( defined BFB_CAM_SCAM_IOP )
   use iop, only: t2sav
   use rgrid, only: nlon
#endif
   use perf_mod
!-----------------------------------------------------------------------
   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   type(advection_state), intent(inout) :: adv_state              ! Advection state data
   real(r8), intent(inout) :: t2(plon,plev,beglat:endlat)         ! temp tendency
   real(r8), intent(inout) :: fu(plon,plev,beglat:endlat)         ! u wind tendency
   real(r8), intent(inout) :: fv(plon,plev,beglat:endlat)         ! v wind tendency

   real(r8), intent(in) :: etamid(plev)                ! vertical coords at midpoints 
   real(r8), intent(inout) :: cwava(plat)                 ! weight applied to global integrals
   real(r8), intent(inout) :: detam(plev)                 ! intervals between vert full levs.
   real(r8), intent(in) :: flx_net(plon,beglat:endlat) ! net flux from physics
   real(r8), intent(in) :: ztodt                       ! twice time step unless nstep=0
!
!---------------------------Local workspace-----------------------------
!
   real(r8) etadot(plon,plevp,beglat:endlat)     ! Vertical motion (slt)
!
! Fourier coefficient arrays which have a latitude index on them for
! multitasking. These arrays are defined in LINEMSAC and used in QUAD
! to compute spectral coefficients. They contain a latitude index so
! that the sums over latitude can be performed in a specified order.
!
   real(r8) grlps1(2*maxm,plat/2)      ! ------------------------------
   real(r8) grlps2(2*maxm,plat/2)      ! |
   real(r8) grt1(2*maxm,plev,plat/2)   ! |
   real(r8) grt2(2*maxm,plev,plat/2)   ! |
   real(r8) grz1(2*maxm,plev,plat/2)   ! |
   real(r8) grz2(2*maxm,plev,plat/2)   ! |
   real(r8) grd1(2*maxm,plev,plat/2)   ! |
   real(r8) grd2(2*maxm,plev,plat/2)   ! |
   real(r8) grfu1(2*maxm,plev,plat/2)  ! |- see quad for definitions
   real(r8) grfu2(2*maxm,plev,plat/2)  ! | 
   real(r8) grfv1(2*maxm,plev,plat/2)  ! |
   real(r8) grfv2(2*maxm,plev,plat/2)  ! |
   real(r8) grut1(2*maxm,plev,plat/2)  ! |
   real(r8) grut2(2*maxm,plev,plat/2)  ! |
   real(r8) grvt1(2*maxm,plev,plat/2)  ! |
   real(r8) grvt2(2*maxm,plev,plat/2)  ! |
   real(r8) grrh1(2*maxm,plev,plat/2)  ! |
   real(r8) grrh2(2*maxm,plev,plat/2)  ! ------------------------------
   real(r8) :: vcour(plev,plat)        ! maximum Courant number in slice
   real(r8) :: vmax2d(plev,plat)       ! max. wind at each level, latitude
   real(r8) :: vmax2dt(plev,plat)      ! max. truncated wind at each lvl,lat
   integer c

   call settau(ztodt/2)
   if(single_column.and.scm_crm_mode) return
!----------------------------------------------------------
! SCANDYN Dynamics scan
!----------------------------------------------------------
!
#if ( defined BFB_CAM_SCAM_IOP )
do c=beglat,endlat
   t2sav(:nlon(c),:,c)= t2(:nlon(c),:,c)
enddo
#endif

if ( single_column ) then
   etadot(1,:,1)=wfldh(:)
else
   call t_startf('scandyn')
   call scandyn(ztodt   ,etadot  ,etamid  ,grlps1  ,grt1    ,  &
                grz1    ,grd1    ,grfu1   ,grfv1   ,grut1   ,  &
                grvt1   ,grrh1   ,grlps2  ,grt2    ,grz2    ,  &
                grd2    ,grfu2   ,grfv2   ,grut2   ,grvt2   ,  &
                grrh2   ,vcour   ,vmax2d,  vmax2dt ,detam   ,  &
                cwava   ,flx_net ,t2      ,fu      ,fv      )
   call t_stopf('scandyn')
endif
!
!----------------------------------------------------------
! SLT scan from south to north
!----------------------------------------------------------
!
   call t_startf('sltrun')
   call scanslt_run(adv_state, ztodt   ,etadot  , detam, etamid, cwava  )
   call t_stopf('sltrun')

   if ( single_column ) then
   call scan2run (ztodt,   cwava,   etamid ,t2      ,fu      ,fv    )
   else
!
!----------------------------------------------------------
! Accumulate spectral coefficients
!----------------------------------------------------------
!
   call t_startf('dynpkg_alloc')
   allocate( vz  (2*lpspt,plev) )
   allocate( d   (2*lpspt,plev) )
   allocate( t   (2*lpspt,plev) )
   allocate( alps(2*lpspt) )
   call t_stopf('dynpkg_alloc')

   call t_startf('dyndrv')
   call dyndrv(grlps1  ,grt1    ,grz1    ,grd1    ,grfu1   ,  &
               grfv1   ,grut1   ,grvt1   ,grrh1   ,grlps2  ,  &
               grt2    ,grz2    ,grd2    ,grfu2   ,grfv2   ,  &
               grut2   ,grvt2   ,grrh2   ,vmax2d  ,vmax2dt ,  &
               vcour, ztodt   )
   call t_stopf('dyndrv')
!
!----------------------------------------------------------
! Second gaussian scan (spectral -> grid)
!----------------------------------------------------------
!
   call t_startf('scan2')
   call scan2run (ztodt,   cwava,   etamid)
   call t_stopf('scan2')

   call t_startf('dynpkg_dealloc')
   deallocate( vz )
   deallocate( d )
   deallocate( t )
   deallocate( alps )
   call t_stopf('dynpkg_dealloc')
endif

   return
end subroutine dynpkg

