!-----------------------------------------------------------------------
!BOP
! !ROUTINE:  epvd --- Calculate absolute potential vorticity
!
! !INTERFACE:
      subroutine epvd( grid, u, v, pt, delp, grav, ae, omega, epv )
! !USES:
      use shr_kind_mod, only : r8 => shr_kind_r8
      use mapz_module, only  : ppme
      use dynamics_vars, only : T_FVDYCORE_GRID
#if defined( SPMD )
      use parutilitiesmodule, only: sumop,  parcollective
      use mod_comm, only :  gid, mp_send3d, mp_recv3d
#endif
      implicit none

! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid   ! grid (for XY decomp)

      real (r8) :: u(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) 
      real (r8) :: v(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) 
      real (r8) :: pt(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km) 
      real (r8) :: delp(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)
      real(r8), intent(in)        :: GRAV    ! Constants, passed as arguments to 
      real(r8), intent(in)        :: AE      ! ensure portability between 
      real(r8), intent(in)        :: OMEGA   ! CAM and GEOS5

! !OUTPUT PARAMETERS:
      real(r8) epv(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)

! !DESCRIPTION:
!     Compute absolute vorticity on the D grid
!        epv = -g * (vort+f0)*dpt/dp
!
! !REVISION HISTORY:
!   WS  99.11.02   Documentation; indentation; jfirstxy:jlastxy
!   WS  00.07.08   Use precision module; Kevin's ghost indices
!   WS  05.02.16   Rewritten for FVdycore_GridCompMod, XY decomposition 
!   WS  05.05.25   Add constants to avoid dependencies on GEOS_Mod
!
! !BUGS:
!   Not yet tested...
!
!EOP
!---------------------------------------------------------------------
!BOC
      real(r8), parameter ::  D0_0                    =  0.0_r8
      real(r8), parameter ::  D1_0                    =  1.0_r8
      real(r8), parameter ::  D2_0                    =  2.0_r8

      real(r8) :: te(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km+1)
      real(r8) :: te2(grid%ifirstxy:grid%ilastxy,grid%km+1)
      real(r8) :: t2(grid%ifirstxy:grid%ilastxy,grid%km)
      real(r8) :: delp2(grid%ifirstxy:grid%ilastxy,grid%km)
      real(r8) :: fx(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy+1)
      real(r8) :: fy(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy)

! Geometric arrays
      real(r8) :: rdx(grid%jfirstxy:grid%jlastxy)  ! 1 / ae*cos(\theta)* dtheta
      real(r8) :: cy(grid%jfirstxy:grid%jlastxy)   ! 1 / ae*cos(\theta)* dlam

      integer  :: i, j, k,  js2g0, jn2g0
      integer  :: iam, myidxy_y, nprxy_x, nprxy_y, dest, src    ! SPMD related
      integer  :: im, jm, km                                    ! problem dimensions
      integer  :: ifirstxy, ilastxy, jfirstxy, jlastxy          ! This PE's intervals
      real(r8) :: c1, c2, rdy

      real(r8), allocatable :: veast(:,:)     ! East halo
      real(r8), allocatable :: unorth(:,:)    ! North halo
      real(r8), allocatable :: fx_sp(:,:), fx_np(:,:)
      real(r8), allocatable :: f0(:)          ! Coriolis force
      real(r8), allocatable :: vort(:,:)      ! Relative vorticity

      im       = grid%im
      jm       = grid%jm
      km       = grid%km

      ifirstxy = grid%ifirstxy
      ilastxy  = grid%ilastxy
      jfirstxy = grid%jfirstxy
      jlastxy  = grid%jlastxy

      iam      = grid%iam
      myidxy_y = grid%myidxy_y
      nprxy_x  = grid%nprxy_x
      nprxy_y  = grid%nprxy_y


      js2g0 = max(2,jfirstxy)
      jn2g0 = min(jm-1,jlastxy)

      allocate(veast(jfirstxy:jlastxy,km))     ! East halo
      allocate(unorth(ifirstxy:ilastxy,km))    ! North halo
      allocate(fx_sp(im,km), fx_np(im,km) )
      allocate(f0(jfirstxy:jlastxy))           ! Coriolis force
      allocate(vort(ifirstxy:ilastxy,jfirstxy:jlastxy))  ! Relative vorticity


! Geometric factors

      do j=jfirstxy,jlastxy
         f0(j)  = D2_0*omega*grid%sinp(j)
      enddo
      rdy   = D1_0/(ae*grid%dp)
      do j=js2g0,jn2g0
         rdx(j) = D1_0/(grid%dl*ae*grid%cosp(j))
         cy(j) =  rdy / grid%cosp(j)
      enddo

      unorth = D0_0
! Periodic boundary  (for the case of no decomposition in X)
      do k=1,km
         do j=jfirstxy,jlastxy
            veast(j,k) = v(ifirstxy,j,k)
         enddo
      enddo

#if defined( SPMD )
      if (nprxy_y > 1) then
! Nontrivial y decomposition
        call mp_send3d( grid%commxy, iam-nprxy_x, iam+nprxy_x, im, jm, km,    &
                        ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km,         &
                        ifirstxy, ilastxy, jfirstxy, jfirstxy, 1, km, u )
      endif
      if (nprxy_x > 1) then
! Nontrivial x decomposition
        dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
        src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
        call mp_send3d( grid%commxy, dest, src, im, jm, km,                   &
                        ifirstxy, ilastxy, jfirstxy, jlastxy, 1,km,          &
                        ifirstxy, ifirstxy, jfirstxy, jlastxy, 1, km, v )
      endif
#endif

! Compute PT at layer edges.

!$omp  parallel do                 &
!$omp  default(shared)             &
!$omp  private(i,j,k,t2,delp2,te2)

      do 1000 j=jfirstxy,jlastxy

        do k=1,km
          do i=ifirstxy,ilastxy
	    t2(i,k) =   pt(i,j,k)
            delp2(i,k) = delp(i,j,k)
          enddo
        enddo

        call ppme(t2,te2,delp2,ilastxy-ifirstxy+1,km)

        do k=1,km+1
          do i=ifirstxy,ilastxy
	    te(i,j,k) = te2(i,k)
          enddo
        enddo

1000  continue


!
! Prepare sum of U-winds for vorticities at pole
!
      fx_sp = D0_0
      fx_np = D0_0
!$omp  parallel do                  &
!$omp  default(shared)              &
!$omp  private(i,k)
      do k=1,km
        if ( jfirstxy == 1 ) then  ! SP
          do i=ifirstxy,ilastxy
            fx_sp(i,k) = u(i,2,k)*grid%cose(2)
          enddo
        endif
        if ( jlastxy == jm ) then  ! NP
          do i=ifirstxy,ilastxy
            fx_np(i,k) = u(i,jm,k)*grid%cose(jm)
          enddo
        endif
      enddo


#if defined( SPMD )
      if ( nprxy_y > 1 ) then
! Non-trivial Y decomposition
        call mp_recv3d( grid%commxy, iam+nprxy_x, im, jm, km,                  &
                        ifirstxy, ilastxy, jlastxy+1, jlastxy+1, 1, km,       &
                        ifirstxy, ilastxy, jlastxy+1, jlastxy+1, 1, km, unorth )
      endif
      if ( nprxy_x > 1 ) then
! Non-trivial X decomposition
        call mp_recv3d( grid%commxy, src, im, jm, km,                          &
                        ilastxy+1, ilastxy+1, jfirstxy, jlastxy, 1, km,       &
                        ilastxy+1, ilastxy+1, jfirstxy, jlastxy, 1, km, veast )
      endif
#endif

#if defined( SPMD )
!
! Collect on all PETs the weighted U-winds at both poles
!
      if (nprxy_x > 1) then
        call parcollective(grid%commxy_x, sumop, im, km, fx_sp)
        call parcollective(grid%commxy_x, sumop, im, km, fx_np)
      endif
#endif

!$omp  parallel do                  &
!$omp  default(shared)              &
!$omp  private(i,j,k,fx,fy,vort,c1,c2)

      do 2000 k=1,km
! Compute relative vorticity
        do j=js2g0,jlastxy
          do i=ifirstxy,ilastxy
            fx(i,j) = u(i,j,k)*grid%cose(j)
          enddo
        enddo
        if ( jlastxy < jm ) then
          do i=ifirstxy,ilastxy
            fx(i,jlastxy+1) = unorth(i,k)*grid%cose(jlastxy+1)
          enddo
        endif

        do j=js2g0,jn2g0
          do i=ifirstxy,ilastxy-1
            fy(i,j) =  v(i+1,j,k) - v(i,j,k)
          enddo
        enddo
        do j=js2g0,jn2g0
          fy(ilastxy,j) = veast(j,k) - v(ilastxy,j,k)
        enddo

        do j=js2g0,jn2g0
          do i=ifirstxy,ilastxy
            vort(i,j) = (fx(i,j)-fx(i,j+1))*cy(j) + fy(i,j)*rdx(j)
          enddo
        enddo

! Vort at poles computed by circulation theorem

        if ( jfirstxy == 1 ) then
          c1 = -SUM(fx_sp(1:im,k))*rdy*grid%rcap
          do i=ifirstxy,ilastxy
            vort(i,  1) = c1
          enddo
        endif 
        if ( jlastxy == jm )  then
          c2 = SUM(fx_np(1:im,k))*rdy*grid%rcap
          do i=ifirstxy,ilastxy
            vort(i,jm) = c2
          enddo
        endif
       
        do j=jfirstxy,jlastxy
          do i=ifirstxy,ilastxy
! Entropy is the thermodynamic variable in the following formulation.
            epv(i,j,k) = grav*(vort(i,j)+f0(j))*(te(i,j,k)-te(i,j,k+1))  &
                       / (pt(i,j,k)*delp(i,j,k))
          enddo
        enddo
!!!        write(iulog,*) "k", k, ifirstxy, jfirstxy, "minmax epv", minval(epv(:,:,k)), &
!!!              maxval(epv(:,:,k)), minloc(epv(:,:,k)), maxloc(epv(:,:,k))
2000  continue

      deallocate(veast)
      deallocate(unorth)
      deallocate(fx_sp,fx_np)
      deallocate(f0)
      deallocate(vort)

      return
      end
