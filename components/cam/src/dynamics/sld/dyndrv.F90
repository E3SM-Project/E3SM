subroutine dyndrv(grlps1  ,grt1    ,grq1    ,grz1    ,grd1    , &
                  grfu1   ,grfv1   ,grlps2  ,grt2    ,grq2    , &
                  grz2    ,grd2    ,grfu2   ,grfv2   ,vmax2d  , &
                  vmax2dt ,vcour   )
!-----------------------------------------------------------------------
!
! Purpose:
! Driving routine for Gaussian quadrature, semi-implicit equation
! solution and linear part of horizontal diffusion.
! The need for this interface routine is to have a multitasking
! driver for the spectral space routines it invokes.
!
! Author:  J. Rosinski
! Modified: P. Worley, October 2002
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe
  use commap
  use time_manager, only: get_step_size
  use spmd_utils, only : iam
  use perf_mod
  implicit none

!------------------------------Arguments--------------------------------
!
  real(r8), intent(inout)   :: grlps1(2*maxm,plat/2)      ! ----------------------------
  real(r8), intent(inout)   :: grt1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grq1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grz1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grd1  (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grfu1 (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grfv1 (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grlps2(2*maxm,plat/2)      ! | definitions: these
  real(r8), intent(inout)   :: grt2  (2*maxm,plev,plat/2) ! | variables are declared here
  real(r8), intent(inout)   :: grq2  (2*maxm,plev,plat/2) ! | for data scoping
  real(r8), intent(inout)   :: grz2  (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grd2  (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grfu2 (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: grfv2 (2*maxm,plev,plat/2) ! |
  real(r8), intent(inout)   :: vmax2d (plev,plat)         ! max. wind at each level, latitude
  real(r8), intent(inout)   :: vmax2dt(plev,plat)         ! max. truncated wind at each lvl,lat
  real(r8), intent(inout)   :: vcour  (plev,plat)         ! maximum Courant number in slice
!
!---------------------------Local workspace-----------------------------
!
  real(r8) ztdtsq(pnmax)    ! dt*(n(n+1)/a^2)
  real(r8) zdt              ! dt
  real(r8) ztdt             ! dt
  integer irow              ! latitude pair index
  integer lm                ! local longitudinal wavenumber index
  integer n                 ! spectral index
  integer k                 ! level index
!
!-----------------------------------------------------------------------
!
  call t_startf ('dyn')

!$OMP PARALLEL DO PRIVATE (IROW)

  do irow=1,plat/2
     call dyn(irow    ,grlps1(1,irow),grt1(1,1,irow),grq1(1,1,irow),grz1(1,1,irow), &
              grd1(1,1,irow),grfu1(1,1,irow),grfv1(1,1,irow),grlps2(1,irow)  , &
                 grt2(1,1,irow), &
              grq2(1,1,irow),grz2 (1,1,irow),grd2 (1,1,irow),grfu2 (1,1,irow), &
                 grfv2(1,1,irow)  )
  end do

  call t_stopf  ('dyn')
!
!-----------------------------------------------------------------------
!
! Build vector with del^2 response function
!
  zdt  = get_step_size()*0.5_r8
  ztdt = 2._r8*zdt
  do n=1,pnmax
     ztdtsq(n) = ztdt*sq(n)
  end do

  call t_startf ('quad')

!$OMP PARALLEL DO PRIVATE(LM)

  do lm=1,numm(iam)
!
! Perform Gaussian quadrature
!
     call quad(lm      ,grlps1  ,grlps2  ,grt1    ,grq1    , &
               grz1    ,grd1    ,grfu1   ,grfv1   ,grt2    , &
               grq2    ,grz2    ,grd2    ,grfu2   ,grfv2   )
   end do

   call t_stopf  ('quad')

   call t_startf ('tstep')

!$OMP PARALLEL DO PRIVATE(LM)

   do lm=1,numm(iam)
!
! Complete time advance, solve vertically coupled semi-implicit system
!
     call tstep(lm      ,ztdtsq  )
!
! Solve vertically de-coupled semi-implicit divergence/vorticity system
! in normal mode space
!
     call vertnm(lm)
!
! Transform divergence and vorticity back from normal mode space and
! complete time advance.
!
     call tstep1(lm      ,zdt     )
  end do

  call t_stopf  ('tstep')
!
! Find out if courant limit has been exceeded.  If so, the limiter will
! be applied in HORDIF
!
  call t_startf('courlim')
  call courlim(vmax2d  ,vmax2dt ,vcour   )
  call t_stopf('courlim')
!
! Linear part of horizontal diffusion
!
  call t_startf('hordif')

!$OMP PARALLEL DO PRIVATE(K)

  do k=1,plev
     call hordif(k       ,ztdt    )
  end do

  call t_stopf('hordif')

  return
end subroutine dyndrv
