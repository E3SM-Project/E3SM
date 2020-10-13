!-----------------------------------------------------------------------
!BOP
! !ROUTINE: d2a3ikj -- Generalized 2nd order D-to-A grid transform (3D)
!                      Output array is i,k,j
!
! !INTERFACE:

      subroutine d2a3dikj( grid, u, v, ua, va )

! !USES:

      use shr_kind_mod, only: r8 => shr_kind_r8
      use dynamics_vars, only : T_FVDYCORE_GRID

#if defined( SPMD )
      use parutilitiesmodule, only : parcollective, sumop
      use mod_comm, only: mp_send3d, mp_recv3d
#endif

      use shr_reprosum_mod, only : shr_reprosum_calc, &
           shr_reprosum_tolExceeded, &
           shr_reprosum_reldiffmax, &
           shr_reprosum_recompute
      use cam_logfile,   only : iulog
      use perf_mod

      implicit none
! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid
      real(r8), intent(in) :: u(grid%ifirstxy:grid%ilastxy,                  &
                                grid%jfirstxy:grid%jlastxy,grid%km) ! U-Wind
      real(r8), intent(in) :: v(grid%ifirstxy:grid%ilastxy,                  &
                                grid%jfirstxy:grid%jlastxy,grid%km) ! V-Wind

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: ua(grid%ifirstxy:grid%ilastxy,grid%km,      &
                                    grid%jfirstxy:grid%jlastxy)     ! U-Wind
      real(r8), intent(inout) :: va(grid%ifirstxy:grid%ilastxy,grid%km,      &
                                    grid%jfirstxy:grid%jlastxy)     ! V-Wind

! !DESCRIPTION:
!
!     This routine performs a second order 
!     interpolation of three-dimensional wind
!     fields on a D grid to an A grid.  !
!
! !REVISION HISTORY:
!     WS  00.12.22 : Creation from d2a3d
!     AAM 01.06.13 : Generalized to 2D decomposition
!     WS  02.04.25 : New mod_comm interfaces
!     WS  03.08.13 : Use unorth for ghosting U (aligned with d2a3dijk)
!     WS  05.07.06 : Simplified interface with grid
!     WS  06.09.08 : isolated magic numbers as F90 parameters
!     PW  08.07.03 : introduced reprosum logic
!     SS  12.10.29 : reprosum is now in csm_share
!
!EOP
!-----------------------------------------------------------------------
!BOC
      real(r8), parameter ::  D0_0                    =  0.0_r8
      real(r8), parameter ::  D0_5                    =  0.5_r8

      integer  :: imh, i, j, k, m, itot, jtot, ltot, ik
      real(r8) :: veast(grid%jfirstxy:grid%jlastxy,grid%km)
      real(r8) :: unorth(grid%ifirstxy:grid%ilastxy,grid%km)

      real(r8) :: uva(grid%ifirstxy:grid%ilastxy,grid%km,2)
      real(r8) :: uvn(grid%km,2), uvs(grid%km,2)
      real(r8) :: rel_diff(2,grid%km,2)
      real(r8),allocatable :: uva_tmp(:)

      integer  :: ifirstxy, ilastxy, jfirstxy, jlastxy, im, jm, km
      integer  :: myidxy_y, myidxy_x, nprxy_x, iam

      logical  :: write_warning

      real(r8), pointer :: coslon(:),sinlon(:)  ! Sine and cosine in longitude

#if defined( SPMD )
      integer dest, src, incount, outcount
#endif

      myidxy_y = grid%myidxy_y
      myidxy_x = grid%myidxy_x
      nprxy_x  = grid%nprxy_x
      iam      = grid%iam

      im       = grid%im
      jm       = grid%jm
      km       = grid%km

      ifirstxy = grid%ifirstxy
      ilastxy  = grid%ilastxy
      jfirstxy = grid%jfirstxy
      jlastxy  = grid%jlastxy

      coslon   => grid%coslon
      sinlon   => grid%sinlon

      itot = ilastxy-ifirstxy+1
      jtot = jlastxy-jfirstxy+1
      imh = im/2

#if defined( SPMD )
! Set ua on A-grid
      call mp_send3d( grid%commxy, iam-nprxy_x, iam+nprxy_x, im, jm, km,      &
                      ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km,           &
                      ifirstxy, ilastxy, jfirstxy, jfirstxy, 1, km, u )
      call mp_recv3d( grid%commxy, iam+nprxy_x, im, jm, km,                   &
                      ifirstxy, ilastxy, jlastxy+1, jlastxy+1, 1, km,        &
                      ifirstxy, ilastxy, jlastxy+1, jlastxy+1, 1, km, unorth )

      if ( jlastxy .lt. jm ) then
!$omp  parallel do private(i, k)

         do k=1,km
            do i=ifirstxy,ilastxy
               ua(i,k,jlastxy) = D0_5 * ( u(i,jlastxy,k) + unorth(i,k) )
            enddo
         enddo
      endif
#endif

!$omp  parallel do private(i,j,k)
      do k=1,km
        do j=jfirstxy, jlastxy-1
          do i=ifirstxy,ilastxy
            ua(i,k,j) = D0_5*(u(i,j,k) + u(i,j+1,k))
          enddo
        enddo
      enddo

! Set va on A-grid

!$omp  parallel do private(j,k)

      do k = 1,km
         do j=jfirstxy,jlastxy
            veast(j,k) = v(ifirstxy,j,k)
         enddo
      enddo

#if defined( SPMD )
      if (itot .ne. im) then
         dest = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
         src  = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
         call mp_send3d( grid%commxy, dest, src, im, jm, km,               &
                         ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km,     &
                         ifirstxy, ifirstxy, jfirstxy, jlastxy, 1, km, v )
         call mp_recv3d( grid%commxy, src, im, jm, km,                     &
                         ilastxy+1, ilastxy+1, jfirstxy, jlastxy, 1, km,  &
                         ilastxy+1, ilastxy+1, jfirstxy, jlastxy, 1, km, veast )
      endif
#endif

!$omp  parallel do private(i,j,k)

      do k=1,km
         do j=jfirstxy, jlastxy
            do i=ifirstxy,ilastxy-1
               va(i,k,j) = D0_5*(v(i,j,k) + v(i+1,j,k))
            enddo
            va(ilastxy,k,j) = D0_5*(v(ilastxy,j,k) + veast(j,k))
         enddo
      enddo

      if (jfirstxy .eq. 1) then
! Was (something like) ...
!            do i=1,imh
!               us(k) = us(k) + (uvaglob(i+imh,k,1)-uvaglob(i,k,1))*sinlon(i)  &
!                     + (uvaglob(i,k,2)-uvaglob(i+imh,k,2))*coslon(i)
!               vs(k) = vs(k) + (uvaglob(i+imh,k,1)-uvaglob(i,k,1))*coslon(i)  & 
!                     + (uvaglob(i+imh,k,2)-uvaglob(i,k,2))*sinlon(i)
!            enddo

!$omp  parallel do private(i,k)
         do k = 1,km
            do i=ifirstxy,min(imh,ilastxy)
               uva(i,k,1) = -ua(i,k,2)*sinlon(i) + va(i,k,2)*coslon(i)
               uva(i,k,2) = -ua(i,k,2)*coslon(i) - va(i,k,2)*sinlon(i)
            enddo
            do i=max(imh+1,ifirstxy),ilastxy
               uva(i,k,1) =  ua(i,k,2)*sinlon(i-imh) - va(i,k,2)*coslon(i-imh)
               uva(i,k,2) =  ua(i,k,2)*coslon(i-imh) + va(i,k,2)*sinlon(i-imh)
            enddo
         enddo

         call t_startf("d2a3dikj_reprosum")
         call shr_reprosum_calc(uva, uvs, itot, itot, 2*km, gbl_count=im, &
                        commid=grid%commxy_x, rel_diff=rel_diff)
         call t_stopf("d2a3dikj_reprosum")

         ! check that "fast" reproducible sum is accurate enough. If not, calculate
         ! using old method
         write_warning = .false.
         if (myidxy_x == 0) write_warning = .true.
         if ( shr_reprosum_tolExceeded('d2a3dikj/South Pole', 2*km, write_warning, &
              iulog, rel_diff) ) then
            if ( shr_reprosum_recompute ) then
               call t_startf("d2a3dikj_sumfix")
               allocate( uva_tmp(im) )
               do m = 1,2
                  do k = 1,km
                     if (rel_diff(1,k,m) > shr_reprosum_reldiffmax) then
                        uva_tmp(:) = D0_0
                        do i = ifirstxy,ilastxy
                           uva_tmp(i) = uva(i,k,m)
                        enddo
#if defined(SPMD)
                        call parcollective(grid%commxy_x,sumop,im,uva_tmp)
#endif
                        uvs(k,m) = D0_0
                        do i = 1,im
                           uvs(k,m) = uvs(k,m) + uva_tmp(i)
                        enddo
                     endif
                  enddo
               enddo
               deallocate( uva_tmp )
               call t_stopf("d2a3dikj_sumfix")
            endif
         endif

!$omp  parallel do private(i,k)
         do k = 1,km
            uvs(k,1) = uvs(k,1)/im
            uvs(k,2) = uvs(k,2)/im
            do i=ifirstxy,min(imh,ilastxy)
               ua(i,k,1) = -uvs(k,1)*sinlon(i) - uvs(k,2)*coslon(i)
               va(i,k,1) =  uvs(k,1)*coslon(i) - uvs(k,2)*sinlon(i)
            enddo
            do i=max(imh+1,ifirstxy),ilastxy
               ua(i,k,1) =  uvs(k,1)*sinlon(i-imh) + uvs(k,2)*coslon(i-imh)
               va(i,k,1) = -uvs(k,1)*coslon(i-imh) + uvs(k,2)*sinlon(i-imh)
            enddo
         enddo

      endif

      if (jlastxy .eq. jm) then
! Was (something like) ...
!            do i=1,imh
!               un(k) = un(k) + (uaglob(i+imh,k,jm-1)-uaglob(i,k,jm-1))*sinlon(i) &
!                     + (vaglob(i+imh,k,jm-1)-vaglob(i,k,jm-1))*coslon(i)
!               vn(k) = vn(k) + (uaglob(i,k,jm-1)-uaglob(i+imh,k,jm-1))*coslon(i) &
!                     + (vaglob(i+imh,k,jm-1)-vaglob(i,k,jm-1))*sinlon(i)
!            enddo

!$omp  parallel do private(i,k)
         do k = 1,km
            do i=ifirstxy,min(imh,ilastxy)
               uva(i,k,1) = -ua(i,k,jm-1)*sinlon(i) - va(i,k,jm-1)*coslon(i)
               uva(i,k,2) =  ua(i,k,jm-1)*coslon(i) - va(i,k,jm-1)*sinlon(i)
            enddo
            do i=max(imh+1,ifirstxy),ilastxy
               uva(i,k,1) =  ua(i,k,jm-1)*sinlon(i-imh) + va(i,k,jm-1)*coslon(i-imh)
               uva(i,k,2) = -ua(i,k,jm-1)*coslon(i-imh) + va(i,k,jm-1)*sinlon(i-imh)
            enddo
         enddo

         call t_startf("d2a3dikj_reprosum")
         call shr_reprosum_calc(uva, uvn, itot, itot, 2*km, gbl_count=im, &
                        commid=grid%commxy_x, rel_diff=rel_diff)
         call t_stopf("d2a3dikj_reprosum")

         ! check that "fast" reproducible sum is accurate enough. If not, calculate
         ! using old method
         write_warning = .false.
         if (myidxy_x == 0) write_warning = .true.
         if ( shr_reprosum_tolExceeded('d2a3dikj/Nouth Pole', 2*km, write_warning, &
              iulog, rel_diff) ) then
            if ( shr_reprosum_recompute ) then
               call t_startf("d2a3dikj_sumfix")
               allocate( uva_tmp(im) )
               do m = 1,2
                  do k = 1,km
                     if (rel_diff(1,k,m) > shr_reprosum_reldiffmax) then
                        uva_tmp(:) = D0_0
                        do i = ifirstxy,ilastxy
                           uva_tmp(i) = uva(i,k,m)
                        enddo
#if defined(SPMD)
                        call parcollective(grid%commxy_x,sumop,im,uva_tmp)
#endif
                        uvn(k,m) = D0_0
                        do i = 1,im
                           uvn(k,m) = uvn(k,m) + uva_tmp(i)
                        enddo
                     endif
                  enddo
               enddo
               deallocate( uva_tmp )
               call t_stopf("d2a3dikj_sumfix")
            endif
         endif

!$omp  parallel do private(i,k)
         do k = 1,km
            uvn(k,1) = uvn(k,1)/im
            uvn(k,2) = uvn(k,2)/im
            do i=ifirstxy,min(imh,ilastxy)
               ua(i,k,jm) = -uvn(k,1)*sinlon(i) + uvn(k,2)*coslon(i)
               va(i,k,jm) = -uvn(k,1)*coslon(i) - uvn(k,2)*sinlon(i)
            enddo
            do i=max(imh+1,ifirstxy),ilastxy
               ua(i,k,jm) =  uvn(k,1)*sinlon(i-imh) - uvn(k,2)*coslon(i-imh)
               va(i,k,jm) =  uvn(k,1)*coslon(i-imh) + uvn(k,2)*sinlon(i-imh)
            enddo
         enddo

      endif

      return
!EOC
      end
!-----------------------------------------------------------------------
