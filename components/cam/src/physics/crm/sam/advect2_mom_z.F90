module advect2_mom_z_mod
  use params, only: asyncid
  implicit none

contains

  subroutine advect2_mom_z(ncrms)
    !       momentum tendency due to the 2nd-order-central vertical advection
    use vars
    use params, only: crm_rknd
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd), allocatable :: fuz(:,:,:,:)
    real(crm_rknd), allocatable :: fvz(:,:,:,:)
    real(crm_rknd), allocatable :: fwz(:,:,:,:)
    integer i, j, k, kc, kb,icrm
    real(crm_rknd) dz25, www, rhoi, dfuz, dfvz, dfwz

    allocate( fuz(ncrms,nx,ny,nz ) )
    allocate( fvz(ncrms,nx,ny,nz ) )
    allocate( fwz(ncrms,nx,ny,nzm) )
#if defined(_OPENACC)
    call prefetch( fuz )
    call prefetch( fvz )
    call prefetch( fwz )
#elif defined(_OPENMP)
    !$omp target enter data map(alloc: fuz)
    !$omp target enter data map(alloc: fvz)
    !$omp target enter data map(alloc: fwz)
#endif
#if defined(_OPENACC)
    !$acc parallel loop collapse(2) async(asyncid)
#elif defined(_OPENMP)
    !$omp target teams distribute parallel do collapse(2)
#endif
    do k = 1 , nz
      do icrm = 1 , ncrms
        uwle(icrm,k) = 0.
        vwle(icrm,k) = 0.
      enddo
    enddo
#if defined(_OPENACC)
    !$acc parallel loop collapse(3) async(asyncid)
#elif defined(_OPENMP)
    !$omp target teams distribute parallel do collapse(2)
#endif
    do j=1,ny
      do i=1,nx
        do icrm = 1 , ncrms
          dz25=1./(4.*dz(icrm))
          fuz(icrm,i,j,1  ) = 0.
          fuz(icrm,i,j,nz ) = 0.
          fvz(icrm,i,j,1  ) = 0.
          fvz(icrm,i,j,nz ) = 0.
          fwz(icrm,i,j,1  ) = 0.
          fwz(icrm,i,j,nzm) = 0.
        end do
      end do
    enddo

    if(RUN3D) then
#if defined(_OPENACC)
      !$acc parallel loop collapse(4) async(asyncid)
#elif defined(_OPENMP)
      !$omp target teams distribute parallel do collapse(4)
#endif
      do k=2,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              dz25=1./(4.*dz(icrm))
              kb = k-1
              rhoi = dz25 * rhow(icrm,k)
              fuz(icrm,i,j,k) = rhoi*(w(icrm,i,j,k)+w(icrm,i-1,j  ,k))*(u(icrm,i,j,k)+u(icrm,i,j,kb))
              fvz(icrm,i,j,k) = rhoi*(w(icrm,i,j,k)+w(icrm,i  ,j-1,k))*(v(icrm,i,j,k)+v(icrm,i,j,kb))
#if defined(_OPENACC)
              !$acc atomic update
#elif defined(_OPENMP)
              !$omp atomic update
#endif
              uwle(icrm,k) = uwle(icrm,k)+fuz(icrm,i,j,k)
#if defined(_OPENACC)
              !$acc atomic update
#elif defined(_OPENMP)
              !$omp atomic update
#endif
              vwle(icrm,k) = vwle(icrm,k)+fvz(icrm,i,j,k)
            end do
          end do
        end do
      end do

    else
#if defined(_OPENACC)
      !$acc parallel loop collapse(4) async(asyncid)
#elif defined(_OPENMP)
      !$omp target teams distribute parallel do collapse(4)
#endif
      do k=2,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              dz25=1./(4.*dz(icrm))
              kb = k-1
              rhoi = dz25 * rhow(icrm,k)
              www = rhoi*(w(icrm,i,j,k)+w(icrm,i-1,j,k))
              fuz(icrm,i,j,k) = www*(u(icrm,i,j,k)+u(icrm,i,j,kb))
              fvz(icrm,i,j,k) = www*(v(icrm,i,j,k)+v(icrm,i,j,kb))
#if defined(_OPENACC)
              !$acc atomic update
#elif defined(_OPENMP)
              !$omp atomic update
#endif
              uwle(icrm,k) = uwle(icrm,k)+fuz(icrm,i,j,k)
#if defined(_OPENACC)
              !$acc atomic update
#elif defined(_OPENMP)
              !$omp atomic update
#endif
              vwle(icrm,k) = vwle(icrm,k)+fvz(icrm,i,j,k)
            end do
          end do
        end do
      end do

    endif
#if defined(_OPENACC)
    !$acc parallel loop collapse(4) async(asyncid)
#elif defined(_OPENMP)
    !$omp target teams distribute parallel do collapse(4)
#endif
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            dz25=1./(4.*dz(icrm))
            kc = k+1
            rhoi = 1./(rho(icrm,k)*adz(icrm,k))
            dfuz = (fuz(icrm,i,j,kc)-fuz(icrm,i,j,k))*rhoi
            dfvz = (fvz(icrm,i,j,kc)-fvz(icrm,i,j,k))*rhoi
#if defined(_OPENACC)
            !$acc atomic update
#elif defined(_OPENMP)
            !$omp atomic update
#endif
            dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-dfuz
#if defined(_OPENACC)
            !$acc atomic update
#elif defined(_OPENMP)
            !$omp atomic update
#endif
            dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-dfvz
          !  dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)-(fuz(icrm,i,j,kc)-fuz(icrm,i,j,k))*rhoi
          !  dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-(fvz(icrm,i,j,kc)-fvz(icrm,i,j,k))*rhoi
            fwz(icrm,i,j,k)=dz25*(w(icrm,i,j,kc)*rhow(icrm,kc)+w(icrm,i,j,k)*rhow(icrm,k))*(w(icrm,i,j,kc)+w(icrm,i,j,k))
          end do
        end do
      end do
    end do
#if defined(_OPENACC)
    !$acc parallel loop collapse(4) async(asyncid)
#elif defined(_OPENMP)
    !$omp target teams distribute parallel do collapse(4)
#endif
    do k=2,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            kb=k-1
            rhoi = 1./(rhow(icrm,k)*adzw(icrm,k))
            dfwz = (fwz(icrm,i,j,k)-fwz(icrm,i,j,kb))*rhoi
#if defined(_OPENACC)
            !$acc atomic update
#elif defined(_OPENMP)
            !$omp atomic update
#endif
            dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na)-dfwz
          !  dwdt(icrm,i,j,k,na)=dwdt(icrm,i,j,k,na)-(fwz(icrm,i,j,k)-fwz(icrm,i,j,kb))*rhoi
          end do
        end do
      end do ! k
    end do
#if defined(_OPENMP)
    !$omp target exit data map(delete: fuz)
    !$omp target exit data map(delete: fvz)
    !$omp target exit data map(delete: fwz)
#endif
    deallocate( fuz )
    deallocate( fvz )
    deallocate( fwz )

  end subroutine advect2_mom_z

end module advect2_mom_z_mod
