!-----------------------------------------------------------------------
!BOP
! !ROUTINE: uv3s_update --  update u3s, v3s (XY decomposition)
!
! !INTERFACE:

   subroutine uv3s_update(grid, dua, u3s, dva, v3s, dt5)

! !USES:

      use shr_kind_mod, only: r8 => shr_kind_r8

#if defined( SPMD )
      use parutilitiesmodule, only : pargatherreal
      use mod_comm, only : mp_send3d, mp_recv3d
#endif
      use cam_history,   only: outfld

      use dynamics_vars, only: T_FVDYCORE_GRID

      implicit none
! !INPUT PARAMETERS:
      type (T_FVDYCORE_GRID), intent(in) :: grid
! dudt on A-grid 
      real(r8),intent(in)  :: dua(grid%ifirstxy:grid%ilastxy,grid%km,grid%jfirstxy:grid%jlastxy)
! dvdt on A-grid 
      real(r8),intent(in)  :: dva(grid%ifirstxy:grid%ilastxy,grid%km,grid%jfirstxy:grid%jlastxy)
      real(r8),intent(in)  :: dt5     ! weighting factor

! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: u3s(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy, &
                                     grid%km)          ! U-Wind on D Grid
      real(r8), intent(inout) :: v3s(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy, &
                                     grid%km)          ! V-Wind on D Grid

! !DESCRIPTION:
!
!     This routine performs the update for the N-S staggered u-wind
!       and the E-W staggered v-wind
!
! !REVISION HISTORY:
!    WS   00.12.22 : Creation from d2a3d
!    SJL  01.01.20 : modifications
!    AAM  01.06.08 : Name change; folding in of v3s update and outfld calls
!    WS   02.04.25 : New mod_comm interfaces
!    WS   02.07.04 : Fixed 2D decomposition bug dest/src for mp_send3d
!    WS   03.07.22 : Removed strip3zatyt4 from use list (no longer used)
!    WS   05.07.14 : Simplified interface with grid argument
!    WS   05.09.23 : Modified for XY decomposition
!
!EOP
!-----------------------------------------------------------------------
!BOC

   integer  :: i, j, k
   integer  :: im, jm, km, ifirstxy, ilastxy, jfirstxy, jlastxy, idim

#if defined( SPMD )
   real(r8) :: duasouth(grid%ifirstxy:grid%ilastxy,grid%km)
   real(r8) :: dvawest(grid%km,grid%jfirstxy:grid%jlastxy)
   integer  :: dest, src
   integer  :: iam, nprxy_x, myidxy_y
#endif
   real(r8) :: tmp
   real(r8) :: u3s_tmp (grid%ifirstxy:grid%ilastxy,grid%km)
   real(r8) :: v3s_tmp (grid%ifirstxy:grid%ilastxy,grid%km)
   real(r8) :: fu3s    (grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)
   real(r8) :: fv3s    (grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)
   real(r8) :: fu3s_tmp(grid%ifirstxy:grid%ilastxy,grid%km)
   real(r8) :: fv3s_tmp(grid%ifirstxy:grid%ilastxy,grid%km)

   fu3s(:,:,:) = 0._r8
   fv3s(:,:,:) = 0._r8

   im     =  grid%im
   jm     =  grid%jm
   km     =  grid%km

   ifirstxy =  grid%ifirstxy
   ilastxy  =  grid%ilastxy
   jfirstxy =  grid%jfirstxy
   jlastxy  =  grid%jlastxy

#if defined( SPMD )
      iam      = grid%iam
      nprxy_x  = grid%nprxy_x
      myidxy_y = grid%myidxy_y
!
! Transfer dua(:,jlast) to the node directly to the north; dva(ifirst, to east)
!
      call mp_send3d( grid%commxy, iam+nprxy_x, iam-nprxy_x, im, km, jm,     &
                      ifirstxy, ilastxy, 1, km, jfirstxy, jlastxy,          &
                      ifirstxy, ilastxy, 1, km, jlastxy, jlastxy, dua )
      call mp_recv3d( grid%commxy, iam-nprxy_x, im, km, jm,                  &
                      ifirstxy, ilastxy, 1, km, jfirstxy-1, jfirstxy-1,     &
                      ifirstxy, ilastxy, 1, km, jfirstxy-1, jfirstxy-1, duasouth )

      dest = myidxy_y*nprxy_x + MOD(iam+1,nprxy_x)
      src  = myidxy_y*nprxy_x + MOD(iam+nprxy_x-1,nprxy_x)
      call mp_send3d( grid%commxy, dest, src, im, km, jm,                    &
                      ifirstxy, ilastxy, 1, km, jfirstxy, jlastxy,          &
                      ilastxy, ilastxy, 1, km, jfirstxy, jlastxy, dva )
      call mp_recv3d( grid%commxy, src, im, km, jm,                          &
                      ifirstxy-1, ifirstxy-1, 1, km, jfirstxy, jlastxy, &
                      ifirstxy-1, ifirstxy-1, 1, km, jfirstxy, jlastxy, dvawest )
#endif

!$omp parallel do private (i, j, k)

      do k = 1, km

!
! Adjust D-grid winds by interpolating A-grid tendencies.
!

        do j = jfirstxy+1, jlastxy
          do i = ifirstxy, ilastxy
             tmp         =  u3s(i,j,k)
             u3s (i,j,k) =  u3s(i,j,k) + dt5*(dua(i,k,j)+dua(i,k,j-1))
             fu3s(i,j,k) = (u3s(i,j,k) - tmp)/(2._r8*dt5)
          enddo
        enddo

        do j = max(jfirstxy,2), min(jlastxy,jm-1)
           do i=ifirstxy+1,ilastxy
              tmp         =  v3s(i,j,k)
              v3s (i,j,k) =  v3s(i,j,k) + dt5*(dva(i,k,j)+dva(i-1,k,j))
              fv3s(i,j,k) = (v3s(i,j,k) - tmp)/(2._r8*dt5)
           enddo
        enddo

#if defined( SPMD )
        if ( jfirstxy .gt. 1 ) then
          do i = ifirstxy, ilastxy
             tmp                =  u3s(i,jfirstxy,k)
             u3s (i,jfirstxy,k) =  u3s(i,jfirstxy,k) +                         &
                         dt5*( dua(i,k,jfirstxy) + duasouth(i,k) )
             fu3s(i,jfirstxy,k) = (u3s(i,jfirstxy,k) - tmp)/(2._r8*dt5)
          enddo
        endif
        do j = max(jfirstxy,2), min(jlastxy,jm-1)
           tmp                =  v3s(ifirstxy,j,k)
           v3s (ifirstxy,j,k) =  v3s(ifirstxy,j,k) + dt5*(dva(ifirstxy,k,j)+dvawest(k,j))
           fv3s(ifirstxy,j,k) = (v3s(ifirstxy,j,k) - tmp)/(2._r8*dt5)
        enddo
#else
        do j = max(jfirstxy,2), min(jlastxy,jm-1)
           tmp         =  v3s(1,j,k)
           v3s (1,j,k) =  v3s(1,j,k) + dt5*(dva(1,k,j)+dva(im,k,j))
           fv3s(1,j,k) = (v3s(1,j,k) - tmp)/(2._r8*dt5)
        enddo
#endif

      enddo

      idim = ilastxy - ifirstxy + 1

!$omp parallel do private (i, j, k, u3s_tmp, v3s_tmp, fu3s_tmp, fv3s_tmp)

      do j = jfirstxy, jlastxy
         do k = 1, km
            do i = ifirstxy, ilastxy
               u3s_tmp (i,k) = u3s (i,j,k)
               v3s_tmp (i,k) = v3s (i,j,k)
               fu3s_tmp(i,k) = fu3s(i,j,k)
               fv3s_tmp(i,k) = fv3s(i,j,k)
            enddo
         enddo

         call outfld ('FU      ', dua(:,:,j), idim, j )
         call outfld ('FV      ', dva(:,:,j), idim, j )
         call outfld ('US      ', u3s_tmp   , idim, j )
         call outfld ('VS      ', v3s_tmp   , idim, j )
         call outfld ('FU_S    ', fu3s_tmp  , idim, j )
         call outfld ('FV_S    ', fv3s_tmp  , idim, j )

      enddo

      return
!EOC
      end subroutine uv3s_update
!-----------------------------------------------------------------------
