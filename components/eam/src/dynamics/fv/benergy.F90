!-----------------------------------------------------------------------
!BOP
! !ROUTINE: benergy --- Calculate the total energy (based on GFDL)
!
! !INTERFACE:

      subroutine benergy(grid,  u,    v,    t3,   delp,              &
                         qqq,   pe,   peln, phis,                    &
                         r_vir, cp,   rg,   tte,  te0 )

! !USES:

      use shr_kind_mod,  only: r8 => shr_kind_r8
      use dynamics_vars, only: T_FVDYCORE_GRID
      use cam_logfile,   only: iulog

#if defined( SPMD )
      use mod_comm,      only: mp_send3d, mp_recv3d
#endif
      implicit none

! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid         ! YZ decomposition

! U-winds
      real(r8), intent(in) ::  u(grid%ifirstxy:grid%ilastxy,         &
                                 grid%jfirstxy:grid%jlastxy,         &
                                 grid%km)
! V-winds
      real(r8), intent(in) ::  v(grid%ifirstxy:grid%ilastxy,         &
                                 grid%jfirstxy:grid%jlastxy,         &
                                 grid%km)

! Temperature (K)
      real(r8), intent(in) :: t3(grid%ifirstxy:grid%ilastxy,         &
                                 grid%jfirstxy:grid%jlastxy,         &
                                 grid%km)

! Delta pressure
      real(r8), intent(in) :: delp(grid%ifirstxy:grid%ilastxy,       &
                                   grid%jfirstxy:grid%jlastxy,       &
                                   grid%km)

! Specific humidity
      real(r8), intent(in) :: qqq(grid%ifirstxy:grid%ilastxy,        &
                                   grid%jfirstxy:grid%jlastxy,       &
                                   grid%km)

! Edge pressure
      real(r8), intent(in) :: pe(grid%ifirstxy:grid%ilastxy,         &
                                 grid%km+1,                          &
                                 grid%jfirstxy:grid%jlastxy)

! Edge pressure
      real(r8), intent(in) :: peln(grid%ifirstxy:grid%ilastxy,       &
                                   grid%km+1,                        &
                                   grid%jfirstxy:grid%jlastxy)

! Surface heights
      real(r8), intent(in) :: phis(grid%ifirstxy:grid%ilastxy,       &
                                   grid%jfirstxy:grid%jlastxy)

      real(r8), intent(in) :: r_vir  ! Virtual effect constant ( rwv/rg-1 )
      real(r8), intent(in) :: cp     ! C_p ( = rg / cappa )
      real(r8), intent(in) :: rg     ! Gas constant for dry air

! !OUTPUT PARAMETERS:

! column integrated Total Energy
      real(r8), intent(out)   :: tte(grid%jm) 
! globally integrated total energy
      real(r8), intent(out)   :: te0

! !DESCRIPTION:
!    Determines the column and globally integrated total energy
!
! !REVISION HISTORY:
!
! SJL 99.04.13 : Delivered as release 0.9.8
! WS  99.05.18 : Added im, jm, km, te, dz as arguments
! WS  99.05.25 : Replaced IMR by IM, JMR by JM-1; removed fvcore.h
! WS  99.10.11 : Ghosted U, now fully limited to jfirst:jlast
! WS  99.11.23 : Pruned te, additional cleaning
! WS  00.05.14 : Renamed ghost indices as per Kevin's definitions
! WS  00.07.13 : Changed PILGRIM API
! WS  00.08.28 : Cosmetic changes
! AAM 00.08.28 : Added kfirst,klast
! WS  00.12.01 : Replaced MPI_ON with SPMD; hs now distributed
! AAM 01.06.15 : Changes for zero diff
! WS  01.12.10 : Ghosted PT
! WS  01.12.31 : Ghosted U,V
! WS  02.07.04 : Fixed 2D decomposition bug dest/src for mp_send3d
! WS  03.10.22 : pmgrid removed (now spmd_dyn)
! WS  03.12.03 : added grid as input argument
! WS  04.10.07 : Removed dependency on spmd_dyn; info now in GRID
! WS  06.05.02 : Rewritten for XY decomposition based on GFDL-code 
! WS  06.06.21 : Extensive debugging of revised version
!
!EOP
!---------------------------------------------------------------------
!BOC

! Local
      real (r8), parameter :: D0_0        = 0.0_r8
      real (r8), parameter :: D0_25       = 0.25_r8
      real (r8), parameter :: D0_5        = 0.5_r8
      real (r8), parameter :: D1_0        = 1.0_r8

      integer   :: im, jm, km, ifirstxy, ilastxy, jfirstxy, jlastxy
      integer   :: iam, myidxy_x, myidxy_y, nprxy_x, nprxy_y, dest, src    ! SPMD related
      integer   :: i, j, k, js1g0, js2g0, jn1g0, jn1g1, jn2g0, ktot, jtot, itot

      real (r8) :: u2(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy+1)
      real (r8) :: v2(grid%ifirstxy:grid%ilastxy+1,grid%jfirstxy:grid%jlastxy)

      real (r8) :: tm(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy)
      real (r8) :: bte(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy)
      real (r8) :: te_sp(grid%ifirstxy:grid%ilastxy,grid%km)
      real (r8) :: te_np(grid%ifirstxy:grid%ilastxy,grid%km)
      real (r8) :: gztop(grid%ifirstxy:grid%ilastxy)
      real (r8) :: xsum(grid%jfirstxy:grid%jlastxy)
      real (r8) :: sp_sum(grid%km), np_sum(grid%km)
      real (r8) :: tm_sp(grid%km), tm_np(grid%km)
      real (r8) :: tmp

      real (r8) :: te(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy, &
                      grid%km)
      real (r8) :: dz(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy, &
                      grid%km)
      real(r8)  :: veast(grid%jfirstxy:grid%jlastxy,grid%km)     ! East halo
      real(r8)  :: unorth(grid%ifirstxy:grid%ilastxy,grid%km)    ! North halo

      im     = grid%im
      jm     = grid%jm
      km     = grid%km

      ifirstxy = grid%ifirstxy
      ilastxy  = grid%ilastxy
      jfirstxy = grid%jfirstxy
      jlastxy  = grid%jlastxy

      iam      = grid%iam
      myidxy_x = grid%myidxy_x
      myidxy_y = grid%myidxy_y
      nprxy_x  = grid%nprxy_x
      nprxy_y  = grid%nprxy_y
      
      js1g0  = max(1,jfirstxy)
      js2g0  = max(2,jfirstxy)
      jn2g0  = min(jm-1,jlastxy)
      jn1g0  = min(jm,jlastxy)
      jn1g1  = min(jm,jlastxy+1)

      itot   = ilastxy - ifirstxy + 1
      jtot   = jlastxy - jfirstxy + 1

#if defined(SPMD)
      call mp_send3d( grid%commxy, iam-nprxy_x, iam+nprxy_x, im, jm, km,      &
                      ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km,           &
                      ifirstxy, ilastxy, jfirstxy, jfirstxy, 1, km, u )
      call mp_recv3d( grid%commxy, iam+nprxy_x, im, jm, km,                   &
                      ifirstxy, ilastxy, jlastxy+1, jlastxy+1, 1, km,        &
                      ifirstxy, ilastxy, jlastxy+1, jlastxy+1, 1, km, unorth )

      if (itot .ne. im) then
         dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
         src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
         call mp_send3d( grid%commxy, dest, src, im, jm, km,                  &
                         ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km,        &
                         ifirstxy, ifirstxy, jfirstxy, jlastxy, 1, km, v )
         call mp_recv3d( grid%commxy, src, im, jm, km,                        &
                         ilastxy+1, ilastxy+1, jfirstxy, jlastxy, 1, km,     &
                         ilastxy+1, ilastxy+1, jfirstxy, jlastxy, 1, km, veast )
      else
!$omp parallel do private(j, k)
         do k = 1,km
            do j=jfirstxy,jlastxy
               veast(j,k) = v(1,j,k)
            enddo
         enddo
      endif
#else
      !$omp parallel do private(j, k)
      do k = 1,km
         do j=1,jm
            veast(j,k) = v(1,j,k)
         enddo
      enddo
#endif


!-----------------------------------------------------------------------------------------------


!$omp parallel do private(i, j, k, u2, v2, tm)
  do k=1,km
!
! Check the poles for consistent values

      do j=js2g0,jlastxy
         do i=ifirstxy,ilastxy
            u2(i,j) = grid%cose(j) * u(i,j,k)**2
         enddo
      enddo

      if ( jlastxy /= jm ) then    ! Pull information out of northern halo
         do i=ifirstxy,ilastxy
            u2(i,jlastxy+1) = grid%cose(jlastxy+1) * unorth(i,k)**2
         enddo
      endif

      do j=js2g0,jn2g0
         do i=ifirstxy,ilastxy
            v2(i,j) = v(i,j,k)**2
         enddo
         v2(ilastxy+1,j) = veast(j,k)**2  ! eastern halo
      enddo

      do j=js2g0,jn2g0
         do i=ifirstxy,ilastxy
            te(i,j,k) = D0_25*((u2(i,j)+u2(i,j+1))*grid%acosu(j) + &
                                v2(i,j) + v2(i+1,j))
         enddo
      enddo

      do j=jfirstxy,jlastxy
         do i=ifirstxy,ilastxy
            tm(i,j) = t3(i,j,k)*(D1_0+r_vir*qqq(i,j,k))
         enddo
      enddo

      do j=js2g0,jn2g0
         do i=ifirstxy, ilastxy
            te(i,j,k) = delp(i,j,k) * ( te(i,j,k) + cp*tm(i,j) )
         enddo
      enddo

      if ( jfirstxy == 1 ) then
        do i=ifirstxy,ilastxy
           te_sp(i,k) = D0_5*u2(i,2)/grid%cose(2)
        enddo
        tm_sp(k) = tm(ifirstxy,1)    ! All tm(:,1) should be the same
      endif

      if ( jlastxy == jm ) then
        do i=ifirstxy,ilastxy
           te_np(i,k)= D0_5*u2(i,jm)/grid%cose(jm)
        enddo
        tm_np(k) = tm(ifirstxy,jm)   ! All tm(:,jm) should be the same
      endif

      do j=jfirstxy,jlastxy
         do i=ifirstxy,ilastxy
            dz(i,j,k) = rg*tm(i,j)
         enddo
      enddo
  enddo


  if ( jfirstxy == 1 ) then
     call par_xsum( grid, te_sp, km, sp_sum )
!$omp parallel do private(i, k, tmp)
     do k=1,km
        tmp = delp(ifirstxy,1,k) * (D0_5*sp_sum(k)/real(im,r8) +  &
                                    cp*tm_sp(k))
        do i=ifirstxy,ilastxy
           te(i,1,k)  = tmp
        enddo
     enddo
  endif 
  if ( jlastxy == jm ) then
     call par_xsum( grid, te_np, km, np_sum )
!$omp parallel do private(i, k, tmp)
     do k=1,km
        tmp = delp(ifirstxy,jm,k) * (D0_5*np_sum(k)/real(im,r8) +& 
                                     cp*tm_np(k))
        do i=ifirstxy,ilastxy
           te(i,jm,k) = tmp
        enddo
     enddo
  endif

  bte = D0_0
!$omp parallel do private(i,j,k,gztop)
  do j=jfirstxy,jlastxy
! Perform vertical integration
     do i=ifirstxy,ilastxy
        gztop(i) = phis(i,j)
        do k=1,km
           gztop(i) = gztop(i) + dz(i,j,k)*(peln(i,k+1,j)-peln(i,k,j))
        enddo
     enddo

     if (j == 1) then
!       gztop(:) should all have identical values   WS 2006.06.22: this checks out
! SP
        tte(1) = pe(ifirstxy,km+1,1)*phis(ifirstxy,1) - pe(ifirstxy,1,1)*gztop(ifirstxy)
        do k=1,km
           tte(1) = tte(1) + te(ifirstxy,1,k)
        enddo
        tte(1)  = grid%acap * tte(1)
     elseif (j == jm) then
!       gztop(:) should all have identical values   WS 2006.06.22: this checks out
! NP
        tte(jm) = pe(ifirstxy,km+1,jm)*phis(ifirstxy,jm) - pe(ifirstxy,1,jm)*gztop(ifirstxy)
        do k=1,km
           tte(jm) = tte(jm) + te(ifirstxy,jm,k)
        enddo
        tte(jm) = grid%acap * tte(jm)
     else
! Interior

        do i=ifirstxy,ilastxy
           bte(i,j) = pe(i,km+1,j)*phis(i,j) - pe(i,1,j)*gztop(i)
        enddo

        do k=1,km
           do i=ifirstxy,ilastxy
              bte(i,j) = bte(i,j) + te(i,j,k)
           enddo
        enddo
     endif
  enddo

  call par_xsum(grid, bte, jtot, xsum)
  
!$omp parallel do private(j)
  do j=js2g0,jn2g0
     tte(j) = xsum(j)*grid%cosp(j)
  enddo

  call par_vecsum(jm, jfirstxy, jlastxy, tte, te0, grid%commxy_y, grid%nprxy_y)

  write(iulog,*) "myidxy_x/y:", myidxy_x, myidxy_y, "The total energy is", te0

!EOC
  end subroutine benergy
!-----------------------------------------------------------------------
