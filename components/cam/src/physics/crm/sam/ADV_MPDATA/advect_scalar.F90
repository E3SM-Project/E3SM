module advect_scalar_mod
  use advect_scalar2D_mod
  use advect_scalar3D_mod
  use openacc_utils
  implicit none

contains

  subroutine advect_scalar (ncrms,f,fadv,flux)

    !     positively definite monotonic advection with non-oscillatory option

    use grid
    use vars, only: u, v, w, rho, rhow
    use params, only: crm_rknd

    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) f(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    real(crm_rknd) flux(ncrms,nz), fadv(ncrms,nz)
    real(crm_rknd), allocatable :: f0(:,:,:,:)
    real(crm_rknd) tmp
    integer i,j,k,icrm

    allocate( f0(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) )
    call prefetch(f0)

    !$acc parallel loop collapse(4) async(asyncid)
    do k = 1 , nzm
      do j = dimy1_s,dimy2_s
        do i = dimx1_s,dimx2_s
          do icrm = 1 , ncrms
            f0(icrm,i,j,k) = f(icrm,i,j,k)
          enddo
        enddo
      enddo
    enddo

    if(RUN3D) then
      call advect_scalar3D(ncrms, f, u, v, w, rho, rhow, flux)
    else
      call advect_scalar2D(ncrms, f, u, w, rho, rhow, flux)
    endif

    !$acc parallel loop collapse(2) async(asyncid)
    do k=1,nzm
      do icrm = 1 , ncrms
        fadv(icrm,k)=0.
      enddo
    enddo
    !$acc parallel loop collapse(4) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            tmp = f(icrm,i,j,k)-f0(icrm,i,j,k)
            !$acc atomic update
            fadv(icrm,k)=fadv(icrm,k)+tmp
          end do
        end do
      end do
    enddo

    deallocate( f0 )

  end subroutine advect_scalar

end module advect_scalar_mod
