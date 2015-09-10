module integxz_mod
contains     
      subroutine integxz(a,b,n1,n2)
   implicit none
    integer, intent(in) :: n1, n2

    real, intent(inout) :: a (n1, n2)
    real, intent(out) :: b (n1, n2) 

    integer :: i, k, ip, im


!cc z smoothing:
      do k=2,n2-1
      do i=1,n1
      b(i,k)=0.25*(a(i,k+1)+2.*a(i,k)+a(i,k-1))
      enddo
      enddo

      do i=1,n1
      b(i,1 )=0.5*(a(i, 2)+a(i,   1))
      b(i,n2)=0.5*(a(i,n2)+a(i,n2-1))
      enddo

      do k=1,n2
      do i=1,n1
      a(i,k)=b(i,k)
      enddo
      enddo

!cc x smoothing:
      do i=1,n1
      im=i-1
      if(im.eq.0) im=n1-1
      ip=i+1
      if(ip.eq.n1+1) ip=2
      do k=1,n2
      b(i,k)=0.25*(a(im,k)+2.*a(i,k)+a(ip,k))
      enddo
      enddo

      do k=1,n2
      do i=1,n1
      a(i,k)=b(i,k)
      enddo
      enddo

      return
    end subroutine integxz

    subroutine integxz_noise(a,b,n1,n2,nz_noise)
    implicit none
    integer, intent(in) :: n1, n2, nz_noise
    real a (n1, n2), b (n1, n2) 
    integer :: k, i, ip, im

!cc z smoothing:
      do k=2,nz_noise
      do i=1,n1
      b(i,k)=0.25*(a(i,k+1)+2.*a(i,k)+a(i,k-1))
      enddo
      enddo

      do k=2,nz_noise
      do i=1,n1
      a(i,k)=b(i,k)
      enddo
      enddo

!c x smoothing:
      do i=1,n1
      im=i-1
      if(im.eq.0) im=n1-1
      ip=i+1
      if(ip.eq.n1+1) ip=2
      do k=2,nz_noise
      b(i,k)=0.25*(a(im,k)+2.*a(i,k)+a(ip,k))
      enddo
      enddo

      do k=2,nz_noise
      do i=1,n1
      a(i,k)=b(i,k)
      enddo
      enddo

      return
    end subroutine integxz_noise

  end module integxz_mod
