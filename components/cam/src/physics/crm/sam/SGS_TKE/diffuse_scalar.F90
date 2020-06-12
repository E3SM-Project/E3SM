module diffuse_scalar_mod
  use diffuse_scalar2D_mod
  use diffuse_scalar3D_mod
  implicit none

contains

  subroutine diffuse_scalar (ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,tkh,f,fluxb,fluxt,fdiff,flux)
    use grid
    use vars, only: rho, rhow
    use params
    use openacc_utils
    implicit none
    integer, intent(in) :: ncrms
    ! input:
    integer :: dimx1_d,dimx2_d,dimy1_d,dimy2_d
    real(crm_rknd) grdf_x(ncrms,nzm)! grid factor for eddy diffusion in x
    real(crm_rknd) grdf_y(ncrms,nzm)! grid factor for eddy diffusion in y
    real(crm_rknd) grdf_z(ncrms,nzm)! grid factor for eddy diffusion in z
    real(crm_rknd) tkh(ncrms,dimx1_d:dimx2_d, dimy1_d:dimy2_d, nzm) ! SGS eddy conductivity
    real(crm_rknd) f(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) ! scalar
    real(crm_rknd) fluxb(ncrms,nx,ny)   ! bottom flux
    real(crm_rknd) fluxt(ncrms,nx,ny)   ! top flux
    real(crm_rknd) fdiff(ncrms,nz)
    real(crm_rknd) flux(ncrms,nz)
    ! Local
    real(crm_rknd), allocatable :: df(:,:,:,:)  ! scalar
    real(crm_rknd) :: tmp
    integer i,j,k,icrm

    allocate( df(ncrms,dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) )
    call prefetch(df)

    !$acc parallel loop collapse(4) async(asyncid)
    do k = 1 , nzm
      do j = dimy1_s , dimy2_s
        do i = dimx1_s , dimx2_s
          do icrm = 1 , ncrms
            df(icrm,i,j,k) = f(icrm,i,j,k)
          enddo
        enddo
      enddo
    enddo

    if(RUN3D) then
      call diffuse_scalar3D (ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,grdf_y,grdf_z,f,fluxb,fluxt,tkh,rho,rhow,flux)
    else
      call diffuse_scalar2D (ncrms,dimx1_d,dimx2_d,dimy1_d,dimy2_d,grdf_x,       grdf_z,f,fluxb,fluxt,tkh,rho,rhow,flux)
    endif

    !$acc parallel loop collapse(2) async(asyncid)
    do k=1,nzm
      do icrm = 1 , ncrms
        fdiff(icrm,k)=0.
      enddo
    enddo
    !$acc parallel loop collapse(2) async(asyncid)
    do k=1,nzm
      do j=1,ny
        do i=1,nx
          do icrm = 1 , ncrms
            tmp = f(icrm,i,j,k)-df(icrm,i,j,k)
            !$acc atomic update
            fdiff(icrm,k)=fdiff(icrm,k)+tmp
          end do
        end do
      end do
    enddo

    deallocate( df )

  end subroutine diffuse_scalar

end module diffuse_scalar_mod
