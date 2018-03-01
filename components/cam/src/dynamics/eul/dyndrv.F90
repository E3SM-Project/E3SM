subroutine dyndrv(grlps1,  grt1,    grz1,    grd1,    grfu1,    &
                  grfv1,   grut1,   grvt1,   grrh1,   grlps2,   &
                  grt2,    grz2,    grd2,    grfu2,   grfv2,    &
                  grut2,   grvt2,   grrh2,   vmax2d,  vmax2dt,  &
                  vcour, ztodt   )
!-----------------------------------------------------------------------
!
! Driving routine for Gaussian quadrature, semi-implicit equation
! solution and linear part of horizontal diffusion.
! The need for this interface routine is to have a multitasking
! driver for the spectral space routines it invokes.
!
!---------------------------Code history--------------------------------
!
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, B. Boville, J. Hack, August 1992
! Reviewed:          D. Williamson, March 1996
! Modified:          P. Worley, September 2002
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use commap
!   use time_manager, only: get_step_size, is_first_step
   use spmd_utils, only: iam
   use perf_mod

   implicit none

!
! Input arguments
!
   real(r8), intent(in) :: grlps1(2*maxm,(plat+1)/2)       ! ----------------------------
   real(r8), intent(inout) :: grt1(2*maxm,plev,(plat+1)/2) ! |
   real(r8), intent(inout) :: grz1(2*maxm,plev,(plat+1)/2) ! |
   real(r8), intent(inout) :: grd1(2*maxm,plev,(plat+1)/2) ! |
   real(r8), intent(in) :: grfu1(2*maxm,plev,(plat+1)/2)   ! |
   real(r8), intent(in) :: grfv1(2*maxm,plev,(plat+1)/2)   ! |
   real(r8), intent(in) :: grut1(2*maxm,plev,(plat+1)/2)   ! |
   real(r8), intent(in) :: grvt1(2*maxm,plev,(plat+1)/2)   ! |
   real(r8), intent(in) :: grrh1(2*maxm,plev,(plat+1)/2)   ! |- see linems and quad for
   real(r8), intent(in) :: grlps2(2*maxm,(plat+1)/2)       ! |  definitions: these variables are
   real(r8), intent(inout) :: grt2(2*maxm,plev,(plat+1)/2) ! |  declared here for data scoping
   real(r8), intent(inout) :: grz2(2*maxm,plev,(plat+1)/2) ! |
   real(r8), intent(inout) :: grd2(2*maxm,plev,(plat+1)/2) ! |
   real(r8), intent(in) :: grfu2(2*maxm,plev,(plat+1)/2)   ! |
   real(r8), intent(in) :: grfv2(2*maxm,plev,(plat+1)/2)   ! |
   real(r8), intent(in) :: grut2(2*maxm,plev,(plat+1)/2)   ! |
   real(r8), intent(in) :: grvt2(2*maxm,plev,(plat+1)/2)   ! |
   real(r8), intent(in) :: grrh2(2*maxm,plev,(plat+1)/2)   ! ----------------------------
   real(r8), intent(inout) :: vmax2d(plev,plat)        ! max. wind at each level, latitude
   real(r8), intent(inout) :: vmax2dt(plev,plat)       ! max. truncated wind at each lvl,lat
   real(r8), intent(inout) :: vcour(plev,plat)         ! maximum Courant number in slice
   real(r8), intent(in) :: ztodt               
!
!---------------------------Local workspace-----------------------------
!
   real(r8) ztdtsq(pnmax)                ! 2dt*(n(n+1)/a^2)
   real(r8) zdt                          ! dt unless nstep = 0
   real(r8) ztdt                         ! 2*zdt (2dt)
   integer irow                      ! latitude pair index
   integer lm                        ! local longitudinal wavenumber index
   integer n                         ! total wavenumber index
   integer k                         ! level index

   call t_startf('dyn')

!$OMP PARALLEL DO PRIVATE (IROW)
   do irow=1,plat/2
      call dyn(irow,   grlps1(:,irow),   grt1(:,:,irow),    &
               grz1(:,:,irow),   grd1(:,:,irow),   &
               grfu1(:,:,irow),  grfv1(:,:,irow),  &
               grut1(:,:,irow),  grvt1(:,:,irow),  &
               grrh1(:,:,irow),  &
               grlps2(:,irow),   grt2(:,:,irow),   &
               grz2(:,:,irow),   grd2(:,:,irow),   &
               grfu2(:,:,irow),  &
               grfv2(:,:,irow),  grut2(:,:,irow),  &
               grvt2(:,:,irow),  grrh2(:,:,irow),ztodt  )
   end do

   call t_stopf('dyn')
!
!-----------------------------------------------------------------------
!
! Build vector with del^2 response function
!

   ztdt = ztodt
   zdt = ztdt/2       	
!   zdt = get_step_size()
!   if (is_first_step()) zdt = .5_r8*zdt
!   ztdt = 2._r8*zdt

   
   do n=1,pnmax
      ztdtsq(n) = ztdt*sq(n)
   end do

   call t_startf ('quad-tstep')

#ifdef OUTER_OMP
!$OMP PARALLEL DO PRIVATE(LM)
#endif
   do lm=1,numm(iam)
!
! Perform Gaussian quadrature
!
      call quad(lm,     zdt,     ztdtsq,  grlps1,  grlps2,  &
                grt1,   grz1,    grd1,    grfu1,   grfv1,   &
                grvt1,  grrh1,   grt2,    grz2,    grd2,   &
                grfu2,  grfv2,   grvt2,   grrh2   )
!
! Complete time advance, solve vertically coupled semi-implicit system
!
      call tstep(lm,zdt,ztdtsq)
   end do
   call t_stopf  ('quad-tstep')
!
! Find out if courant limit has been exceeded.  If so, the limiter will be
! applied in HORDIF
!
   call t_startf('courlim')
   call courlim(vmax2d,  vmax2dt, vcour   )
   call t_stopf('courlim')
!
! Linear part of horizontal diffusion
!
   call t_startf('hordif')

!$OMP PARALLEL DO PRIVATE(K)
   do k=1,plev
      call hordif(k,ztdt)
   end do

   call t_stopf('hordif')

   return
end subroutine dyndrv
