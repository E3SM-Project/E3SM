!-----------------------------------------------------------------------
!BOP
! !ROUTINE: d2a3ijk -- Generalized 2nd order D-to-A grid transform (3D)
!                      Output array is i,j,k
!
! !INTERFACE:

      subroutine d2a3dijk(grid, u, v, ua, va )

! !USES:

      use shr_kind_mod, only: r8 => shr_kind_r8
      use dynamics_vars, only : T_FVDYCORE_GRID

#if defined( SPMD )
      use parutilitiesmodule, only : parcollective, sumop
      use mod_comm, only: mp_send3d, mp_recv3d
#endif

      implicit none
! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid
      real(r8), intent(in) :: u(grid%ifirstxy:grid%ilastxy,            &
                                grid%jfirstxy:grid%jlastxy,grid%km) ! U-Wind ghosted N1
      real(r8), intent(in) :: v(grid%ifirstxy:grid%ilastxy,            &
                                grid%jfirstxy:grid%jlastxy,grid%km) ! V-Wind

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: ua(grid%ifirstxy:grid%ilastxy,        &
                                    grid%jfirstxy:grid%jlastxy,grid%km) ! U-Wind
      real(r8), intent(inout) :: va(grid%ifirstxy:grid%ilastxy,        &
                                    grid%jfirstxy:grid%jlastxy,grid%km) ! V-Wind


! !DESCRIPTION:
!
!     This routine performs a second order 
!     interpolation of three-dimensional wind
!     fields on a D grid to an A grid.  !
!
! !REVISION HISTORY:
!     WS  00.12.22 : Creation from d2a3d
!     AAM 01.06.13 : Generalized to 2D decomposition
!     WS  02.04.25 : Newest mod_comm interfaces
!     WS  05.07.06 : Simplified interface with grid
!     WS  06.09.08 : Isolated magic numbers as F90 parameters
!
!EOP
!-----------------------------------------------------------------------
!BOC

      real(r8), parameter ::  D0_0                    =  0.0_r8
      real(r8), parameter ::  D0_5                    =  0.5_r8

      integer  :: imh, i, j, k, itot, jtot, ltot, lbegin, lend, ik
      real(r8) :: un(grid%km), vn(grid%km), us(grid%km), vs(grid%km)
      real(r8) :: veast(grid%jfirstxy:grid%jlastxy,grid%km)
      real(r8) :: unorth(grid%ifirstxy:grid%ilastxy,grid%km)
      real(r8) :: uvaglob(grid%im,grid%km,4)
      real(r8) :: uvaloc(grid%ifirstxy:grid%ilastxy,grid%km,4)
      real(r8) :: uaglob(grid%im),vaglob(grid%im)

      integer  :: im, jm, km, ifirstxy, ilastxy, jfirstxy, jlastxy
      integer  :: myidxy_y, myidxy_x, nprxy_x, iam

      real(r8), pointer :: coslon(:), sinlon(:)

#if defined( SPMD )
      integer  :: dest, src, incount, outcount
#endif

      im   = grid%im
      jm   = grid%jm
      km   = grid%km

      ifirstxy = grid%ifirstxy
      ilastxy  = grid%ilastxy
      jfirstxy = grid%jfirstxy
      jlastxy  = grid%jlastxy

      myidxy_x = grid%myidxy_x
      myidxy_y = grid%myidxy_y
      nprxy_x  = grid%nprxy_x
      iam      = grid%iam

      coslon  => grid%coslon
      sinlon  => grid%sinlon

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
               ua(i,jlastxy,k) = D0_5 * ( u(i,jlastxy,k) + unorth(i,k) )
            enddo
         enddo
      endif
#endif

!$omp  parallel do private(i,j,k)

      do k=1,km
        do j=jfirstxy, jlastxy-1
          do i=ifirstxy,ilastxy
            ua(i,j,k) = D0_5*(u(i,j,k) + u(i,j+1,k))
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
         call mp_send3d( grid%commxy, dest, src, im, jm, km,             &
                         ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km,   &
                         ifirstxy, ifirstxy, jfirstxy, jlastxy, 1, km, v )
         call mp_recv3d( grid%commxy, src, im, jm, km,                                 & 
                         ilastxy+1, ilastxy+1, jfirstxy, jlastxy, 1, km,&
                         ilastxy+1, ilastxy+1, jfirstxy, jlastxy, 1, km, veast )
      endif
#endif

!$omp  parallel do private(i,j,k)

      do k=1,km
         do j=jfirstxy, jlastxy
            do i=ifirstxy,ilastxy-1
               va(i,j,k) = D0_5*(v(i,j,k) + v(i+1,j,k))
            enddo
            va(ilastxy,j,k) = D0_5*(v(ilastxy,j,k) + veast(j,k))
         enddo
      enddo

!$omp  parallel do private(i,ik,k)

      do ik=1,4
         do k=1,km
            do i=1,im
               uvaglob(i,k,ik) = D0_0
            enddo
         enddo
      enddo

      lbegin = 0
      lend = 0
      if (jfirstxy .eq. 1) then
!$omp  parallel do private(i,k)
         do k = 1,km
            do i=ifirstxy,ilastxy
               uvaloc(i,k,1) = ua(i,2,k)
               uvaloc(i,k,2) = va(i,2,k)
               uvaglob(i,k,1) = ua(i,2,k)
               uvaglob(i,k,2) = va(i,2,k)
            enddo
         enddo
         lbegin = 1
         lend = 2
      endif

      if (jlastxy .eq. jm) then
!$omp  parallel do private(i,k)
         do k = 1,km
            do i=ifirstxy,ilastxy
               uvaloc(i,k,3) = ua(i,jm-1,k)
               uvaloc(i,k,4) = va(i,jm-1,k)
               uvaglob(i,k,3) = ua(i,jm-1,k)
               uvaglob(i,k,4) = va(i,jm-1,k)
            enddo
         enddo
         lbegin = 3
         lend = 4
      endif
      if (jtot .eq. jm) lbegin=1

#if defined( SPMD )
      if (itot .ne. im) then
         ltot = lend-lbegin+1
         if (jfirstxy .eq. 1 .or. jlastxy .eq. jm) then
            call parcollective(grid%commxy_x, sumop, im, km, ltot, &
                               uvaglob(:,:,lbegin:lend))
         endif
      endif
#endif

      if ( jfirstxy .eq. 1 ) then
! Projection at SP
!$omp  parallel do private(i,k,uaglob,vaglob)
         do k=1,km
            us(k) = D0_0
            vs(k) = D0_0
            do i=1,imh
               us(k) = us(k) + (uvaglob(i+imh,k,1)-uvaglob(i,k,1))*sinlon(i)  &
                     + (uvaglob(i,k,2)-uvaglob(i+imh,k,2))*coslon(i)
               vs(k) = vs(k) + (uvaglob(i+imh,k,1)-uvaglob(i,k,1))*coslon(i)  & 
                     + (uvaglob(i+imh,k,2)-uvaglob(i,k,2))*sinlon(i)
            enddo

            us(k) = us(k)/im
            vs(k) = vs(k)/im
            do i=1,imh
               uaglob(i)   = -us(k)*sinlon(i) - vs(k)*coslon(i)
               vaglob(i)   =  us(k)*coslon(i) - vs(k)*sinlon(i)
               uaglob(i+imh) = -uaglob(i)
               vaglob(i+imh) = -vaglob(i)
            enddo
            do i=ifirstxy,ilastxy
               ua(i,1,k) = uaglob(i)
               va(i,1,k) = vaglob(i)
            enddo
         enddo
      endif

      if ( jlastxy .eq. jm ) then
! Projection at NP
!$omp  parallel do private(i,k,uaglob,vaglob)
         do k=1,km
            un(k) = D0_0
            vn(k) = D0_0
            do i=1,imh
               un(k) = un(k) + (uvaglob(i+imh,k,3)-uvaglob(i,k,3))*sinlon(i) &
                     + (uvaglob(i+imh,k,4)-uvaglob(i,k,4))*coslon(i)
               vn(k) = vn(k) + (uvaglob(i,k,3)-uvaglob(i+imh,k,3))*coslon(i) &
                     + (uvaglob(i+imh,k,4)-uvaglob(i,k,4))*sinlon(i)
            enddo

            un(k) = un(k)/im
            vn(k) = vn(k)/im
            do i=1,imh
               uaglob(i) = -un(k)*sinlon(i) + vn(k)*coslon(i)
               vaglob(i) = -un(k)*coslon(i) - vn(k)*sinlon(i)
               uaglob(i+imh) = -uaglob(i)
               vaglob(i+imh) = -vaglob(i)
            enddo
            do i=ifirstxy,ilastxy
               ua(i,jm,k) = uaglob(i)
               va(i,jm,k) = vaglob(i)
            enddo
         enddo
      endif

      return
!EOC
      end
!-----------------------------------------------------------------------
