module press_grad_mod
  use task_util_mod
  use bound_duvdt_mod
  implicit none

contains

  subroutine press_grad(ncrms)
    !       pressure term of the momentum equations
    use vars
    use params, only: dowallx, dowally
    implicit none
    integer, intent(in) :: ncrms
    real *8 rdx,rdy,rdz,dpdx,dpdy,dpdz
    integer i,j,k,kb,jb,ib, icrm

    rdx=1./dx
    rdy=1./dy
#if defined(_OPENACC)
    !$acc parallel loop collapse(4) async(asyncid)
#elif defined(_OPENMP)
    !$omp target teams distribute parallel do collapse(4)
#endif
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            kb=max(1,k-1)
            rdz = 1./(dz(icrm)*adzw(icrm,k))
            jb=j-YES3D
            ib=i-1
            dpdx = (p(icrm,i,j,k)-p(icrm,ib,j,k))*rdx
            dpdy = (p(icrm,i,j,k)-p(icrm,i,jb,k))*rdy
            dpdz = (p(icrm,i,j,k)-p(icrm,i,j,kb))*rdz
#if defined(_OPENACC)
            !$acc atomic update
#elif defined(_OPENMP)
            !$omp atomic update
#endif
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-dpdx
#if defined(_OPENACC)
            !$acc atomic update
#elif defined(_OPENMP)
            !$omp atomic update
#endif      
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-dpdy
#if defined(_OPENACC)
            !$acc atomic update
#elif defined(_OPENMP)
            !$omp atomic update
#endif      
            dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na)-dpdz
          end do ! i
        end do ! j
      end do ! k
    enddo
#if defined(_OPENACC)
    !$acc parallel loop collapse(4) async(asyncid)
#elif defined(_OPENMP)
    !$omp target teams distribute parallel do collapse(4)
#endif
    do k=1,nzm
      do j=1-YES3D,ny !bloss: 0,n* fixes computation of dp/d* in stats.
        do i=0,nx
          do icrm = 1 , ncrms
#if defined(_OPENACC)
            !$acc atomic update
#elif defined(_OPENMP)
            !$omp atomic update
#endif      
            p(icrm,i,j,k)=p(icrm,i,j,k)*rho(icrm,k)  ! convert p'/rho to p'
          end do
        end do
      end do
    enddo

    if(dowallx.and.mod(rank,nsubdomains_x).eq.0) then
#if defined(_OPENACC)
      !$acc parallel loop collapse(3) async(asyncid)
#elif defined(_OPENMP)
      !$omp target teams distribute parallel do collapse(3)
#endif
      do k=1,nzm
        do j=1,ny
            do icrm = 1 , ncrms
            dudt(icrm,1,j,k,na) = 0.
          end do
        end do
      enddo

    end if

    if(dowally.and.RUN3D.and.rank.lt.nsubdomains_x) then
#if defined(_OPENACC)
      !$acc parallel loop collapse(3) async(asyncid)
#elif defined(_OPENMP)
      !$omp target teams distribute parallel do collapse(3)
#endif
      do k=1,nzm
        do i=1,nx
          do icrm = 1 , ncrms
            dvdt(icrm,i,1,k,na) = 0.
          end do
        end do
      enddo

    end if

    call bound_duvdt(ncrms)

  end subroutine press_grad

end module press_grad_mod
