module zero_mean_mod
contains  
      subroutine zero_mean(thn,nx,nz,nz_noise,hg)
    implicit none
    integer, intent(in) :: nx, nz, nz_noise
    real, intent(inout) :: thn (nx, nz) 
    real, intent(in) ::  hg (nz) 
    integer :: i, k
    real :: sum, sum1, anorm, del
      anorm=real(nx-1)
      sum=0.
      sum1=0.
      do k=2,nz_noise
      del=.5*(hg(k+1)-hg(k-1))
      sum1=sum1+del
      do i=1,nx-1
      sum=sum+del*thn(i,k)/anorm
      enddo
      enddo
      sum=sum/sum1
      do i=1,nx
      do k=2,nz_noise
      thn(i,k)=thn(i,k)-sum
      enddo
      enddo

      return
    end subroutine zero_mean
  end module zero_mean_mod
