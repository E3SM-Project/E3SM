
#define _IDX(l1,u1,i1,l2,u2,i2,l3,u3,i3,l4,u4,i4) \
           ( ((i4)-(l4))*((u3)-(l3)+1)*((u2)-(l2)+1)*((u1)-(l1)+1) + \
             ((i3)-(l3))              *((u2)-(l2)+1)*((u1)-(l1)+1) + \
             ((i2)-(l2))                            *((u1)-(l1)+1) + \
             ((i1)-(l1)) + 1 )

module bound_exchange_mod
  use params, only: asyncid
  implicit none

contains

  subroutine bound_exchange(ncrms,f,dimx1,dimx2,dimy1,dimy2,dimz,i_1, i_2, j_1, j_2, id)
    ! periodic boundary exchange
    use grid
    use params, only: crm_rknd
    use openacc_utils
    implicit none
    integer dimx1, dimx2, dimy1, dimy2, dimz, ncrms
    integer i_1, i_2, j_1, j_2
    real(crm_rknd) f(ncrms,dimx1:dimx2, dimy1:dimy2, dimz)
    integer id   ! id of the sent field (dummy variable)
    real(crm_rknd), allocatable :: buffer(:)  ! buffer for sending data
    integer i, j, k, n, icrm
    integer i1, i2, j1, j2

    allocate(buffer((nx+ny)*3*nz*ncrms))
    call prefetch( buffer )

    i1 = i_1 - 1
    i2 = i_2 - 1
    j1 = j_1 - 1
    j2 = j_2 - 1

    !----------------------------------------------------------------------
    !  Send buffers to neighbors
    !----------------------------------------------------------------------

    if(RUN3D) then
      ! "North" -> "South":
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,dimz
        do j=ny-j1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              n = _IDX( 1,ncrms,icrm , 1,nx,i , ny-j1,ny,j , 1,dimz,k )
              buffer(n) = f(icrm,i,j,k)
            end do
          end do
        end do
      end do
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,dimz
        do j=-j1,0
          do i=1,nx
            do icrm = 1 , ncrms
              n = _IDX( 1,ncrms,icrm , 1,nx,i , -j1,0,j , 1,dimz,k )
              f(icrm,i,j,k) = buffer(n)
            end do
          end do
        end do
      end do

      ! "North-East" -> "South-West":
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,dimz
        do j=ny-j1,ny
          do i=nx-i1,nx
            do icrm = 1 , ncrms
              n = _IDX( 1,ncrms,icrm , nx-i1,nx,i , ny-j1,ny,j , 1,dimz,k )
              buffer(n) = f(icrm,i,j,k)
            end do
          end do
        end do
      end do
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,dimz
        do j=-j1,0
          do i=-i1,0
            do icrm = 1 , ncrms
              n = _IDX( 1,ncrms,icrm , -i1,0,i , -j1,0,j , 1,dimz,k )
              f(icrm,i,j,k) = buffer(n)
            end do
          end do
        end do
      end do

      ! "South-East" -> "North-West":
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,dimz
        do j=1,1+j2
          do i=nx-i1,nx
            do icrm = 1 , ncrms
              n = _IDX( 1,ncrms,icrm , nx-i1,nx,i , 1,1+j2,j , 1,dimz,k )
              buffer(n) = f(icrm,i,j,k)
            end do
          end do
        end do
      end do
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,dimz
        do j=nyp1,nyp1+j2
          do i=-i1,0
            do icrm = 1 , ncrms
              n = _IDX( 1,ncrms,icrm , -i1,0,i , nyp1,nyp1+j2,j , 1,dimz,k )
              f(icrm,i,j,k) = buffer(n)
            end do
          end do
        end do
      end do

      ! "South" -> "North":
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,dimz
        do j=1,1+j2
          do i=1,nx
            do icrm = 1 , ncrms
              n = _IDX( 1,ncrms,icrm , 1,nx,i , 1,1+j2,j , 1,dimz,k )
              buffer(n) = f(icrm,i,j,k)
            end do
          end do
        end do
      end do
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,dimz
        do j=nyp1,nyp1+j2
          do i=1,nx
            do icrm = 1 , ncrms
              n = _IDX( 1,ncrms,icrm , 1,nx,i , nyp1,nyp1+j2,j , 1,dimz,k )
              f(icrm,i,j,k) = buffer(n)
            end do
          end do
        end do
      end do

      ! "South-West" -> "North-East":
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,dimz
        do j=1,1+j2
          do i=1,1+i2
            do icrm = 1 , ncrms
              n = _IDX( 1,ncrms,icrm , 1,1+i2,i , 1,1+j2,j , 1,dimz,k )
              buffer(n) = f(icrm,i,j,k)
            end do
          end do
        end do
      end do
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,dimz
        do j=nyp1,nyp1+j2
          do i=nxp1,nxp1+i2
            do icrm = 1 , ncrms
              n = _IDX( 1,ncrms,icrm , nxp1,nxp1+i2,i , nyp1,nyp1+j2,j , 1,dimz,k )
              f(icrm,i,j,k) = buffer(n)
            end do
          end do
        end do
      end do


      ! To "North-West" -> "South-East":
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,dimz
        do j=ny-j1,ny
          do i=1,1+i2
            do icrm = 1 , ncrms
              n = _IDX( 1,ncrms,icrm , 1,1+i2,i , ny-j1,ny,j , 1,dimz,k )
              buffer(n) = f(icrm,i,j,k)
            end do
          end do
        end do
      end do
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,dimz
        do j=-j1,0
          do i=nxp1,nxp1+i2
            do icrm = 1 , ncrms
              n = _IDX( 1,ncrms,icrm , nxp1,nxp1+i2,i , -j1,0,j , 1,dimz,k )
              f(icrm,i,j,k) = buffer(n)
            end do
          end do
        end do
      end do


    endif

    !  "East" -> "West":
    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,dimz
      do j=1,ny
        do i=nx-i1,nx
          do icrm = 1 , ncrms
            n = _IDX( 1,ncrms,icrm , nx-i1,nx,i , 1,ny,j , 1,dimz,k )
            buffer(n) = f(icrm,i,j,k)
          end do
        end do
      end do
    end do
    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,dimz
      do j=1,ny
        do i=-i1,0
          do icrm = 1 , ncrms
            n = _IDX( 1,ncrms,icrm , -i1,0,i , 1,ny,j , 1,dimz,k )
            f(icrm,i,j,k) = buffer(n)
          end do
        end do
      end do
    end do

    ! "West" -> "East":
    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,dimz
      do j=1,ny
        do i=1,1+i2
          do icrm = 1 , ncrms
            n = _IDX( 1,ncrms,icrm , 1,1+i2,i , 1,ny,j , 1,dimz,k )
            buffer(n) = f(icrm,i,j,k)
          end do
        end do
      end do
    end do
    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,dimz
      do j=1,ny
        do i=nxp1,nxp1+i2
          do icrm = 1 , ncrms
            n = _IDX( 1,ncrms,icrm , nxp1,nxp1+i2,i , 1,ny,j , 1,dimz,k )
            f(icrm,i,j,k) = buffer(n)
          end do
        end do
      end do
    end do

    deallocate(buffer)

  end subroutine bound_exchange


end module bound_exchange_mod
