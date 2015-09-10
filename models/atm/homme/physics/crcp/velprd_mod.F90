module velprd_mod
  use grid_init_mod, only : dt, dxi, dzi 
contains  
! ux,uxp,uxa,uz,uzp,uza
  subroutine velprd_1(ux,uxp,uxa,uz,uzp,uza,rho, nx, nz)

    implicit none
    integer, intent(in) :: nx,nz
    real :: ux(nx,nz), uz(nx,nz) 
    real ::  uxp(nx,nz), uzp(nx,nz) 
    real ::  uxa(nx+1,nz), uza(nx,nz+1), rho(nx,nz) 
    integer :: i, k
    real :: dtdx, dtdz

      dtdx = dt*dxi
      dtdz = dt*dzi

       do k=1,nz
       do i=2,nx
       uxa(i,k) =(0.75*(ux(i-1,k)*rho(i-1,k)+ux(i,k)*rho(i,k))  &
            - 0.25*(uxp(i-1,k)*rho(i-1,k)+uxp(i,k)*rho(i,k)) )*dtdx 
       enddo
!cc cyclic in horizontal
       uxa(1,k) = uxa(nx,k)
       uxa(nx+1,k) = uxa(2,k)
       enddo
          
       do k=2,nz
       do i=1,nx
       uza(i,k) =(0.75*(uz(i,k-1)*rho(i,k-1)+uz(i,k)*rho(i,k)) &
               - 0.25*(uzp(i,k-1)*rho(i,k-1)+uzp(i,k)*rho(i,k)) )*dtdz 
       enddo
       enddo
!cc zero flux in vertical
       do i=1,nx
       uza(i,1) = - uza(i,2)
       uza(i,nz+1) = - uza(i,nz)
       enddo

       return
     end subroutine velprd_1
   end module velprd_mod
