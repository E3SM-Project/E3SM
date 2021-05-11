module advect2_mom_xy_mod
  use params, only: asyncid
  implicit none

contains

  subroutine advect2_mom_xy(ncrms)

    !        momentum tendency due to 2nd-order-central horizontal advection

    use vars
    use params, only: crm_rknd

    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) dx25, dy25, irho, fu1, fu2, fv1, fv2, fw1, fw2
    integer i, j, k, kc, kcu, ic, jb, ib, jc,icrm

    dx25 = 0.25D0 / dx
    dy25 = 0.25D0 / dy

    if(RUN3D) then

      !$acc parallel loop collapse(4) async(asyncid)
      do k = 1,nzm
        do j = 1, ny
          do i = 1, nx
            do icrm = 1 , ncrms
              kc= k+1
              kcu =min(kc, nzm)
              irho = 1./(rhow(icrm,kc)*adzw(icrm,kc))
              jb = j-1
              ic = i+1
              fu1 = dx25*(u(icrm,ic-1,j,k)+u(icrm,i-1,j,k))*(u(icrm,i-1,j,k)+u(icrm,ic-1,j,k))
              fu2 = dx25*(u(icrm,ic  ,j,k)+u(icrm,i  ,j,k))*(u(icrm,i  ,j,k)+u(icrm,ic  ,j,k))
              dudt(icrm,i,j,k,na)  = dudt(icrm,i,j,k,na)  - (fu2-fu1)
              fv1 = dx25*(u(icrm,ic-1,j,k)+u(icrm,ic-1,jb,k))*(v(icrm,i-1,j,k)+v(icrm,ic-1,j,k))
              fv2 = dx25*(u(icrm,ic  ,j,k)+u(icrm,ic  ,jb,k))*(v(icrm,i  ,j,k)+v(icrm,ic  ,j,k))
              dvdt(icrm,i,j,k,na)  = dvdt(icrm,i,j,k,na)  - (fv2-fv1)
              fw1 = dx25*(u(icrm,ic-1,j,k)*rho(icrm,k)*adz(icrm,k)+u(icrm,ic-1,j,kcu)*rho(icrm,kcu)*adz(icrm,kcu))*(w(icrm,i-1,j,kc)+w(icrm,ic-1,j,kc))
              fw2 = dx25*(u(icrm,ic  ,j,k)*rho(icrm,k)*adz(icrm,k)+u(icrm,ic  ,j,kcu)*rho(icrm,kcu)*adz(icrm,kcu))*(w(icrm,i  ,j,kc)+w(icrm,ic  ,j,kc))
              dwdt(icrm,i,j,kc,na) = dwdt(icrm,i,j,kc,na)-irho*(fw2-fw1)

              jc = j+1
              ib = i-1
              fu1 = dy25*(v(icrm,i,jc-1,k)+v(icrm,ib,jc-1,k))*(u(icrm,i,j-1,k)+u(icrm,i,jc-1,k))
              fu2 = dy25*(v(icrm,i,jc  ,k)+v(icrm,ib,jc  ,k))*(u(icrm,i,j  ,k)+u(icrm,i,jc  ,k))
              dudt(icrm,i,j,k,na) = dudt(icrm,i,j,k,na) - (fu2-fu1)
              fv1 = dy25*(v(icrm,i,jc-1,k)+v(icrm,i,j-1,k))*(v(icrm,i,j-1,k)+v(icrm,i,jc-1,k))
              fv2 = dy25*(v(icrm,i,jc  ,k)+v(icrm,i,j  ,k))*(v(icrm,i,j  ,k)+v(icrm,i,jc  ,k))
              dvdt(icrm,i,j,k,na) = dvdt(icrm,i,j,k,na) - (fv2-fv1)
              fw1 = dy25*(v(icrm,i,jc-1,k)*rho(icrm,k)*adz(icrm,k)+v(icrm,i,jc-1,kcu)*rho(icrm,kcu)*adz(icrm,kcu))*(w(icrm,i,j-1,kc)+w(icrm,i,jc-1,kc))
              fw2 = dy25*(v(icrm,i,jc  ,k)*rho(icrm,k)*adz(icrm,k)+v(icrm,i,jc  ,kcu)*rho(icrm,kcu)*adz(icrm,kcu))*(w(icrm,i,j  ,kc)+w(icrm,i,jc  ,kc))
              dwdt(icrm,i,j,kc,na)= dwdt(icrm,i,j,kc,na)-irho*(fw2-fw1)
            end do
          end do
        end do ! k
      end do

    else

      j=1
      !$acc parallel loop collapse(3) async(asyncid)
      do k = 1,nzm
        do i = 1, nx
          do icrm = 1 , ncrms
            kc= k+1
            kcu =min(kc, nzm)
            irho = 1./(rhow(icrm,kc)*adzw(icrm,kc))
            ic = i+1
            fu1 = dx25*(u(icrm,ic-1,j,k)+u(icrm,i-1,j,k))*(u(icrm,i-1,j,k)+u(icrm,ic-1,j,k))
            fu2 = dx25*(u(icrm,ic  ,j,k)+u(icrm,i  ,j,k))*(u(icrm,i  ,j,k)+u(icrm,ic  ,j,k))
            dudt(icrm,i,j,k,na)  = dudt(icrm,i,j,k,na)  - (fu2-fu1)
            fv1 = dx25*(u(icrm,ic-1,j,k)+u(icrm,i-1,j,k))*(v(icrm,i-1,j,k)+v(icrm,ic-1,j,k))
            fv2 = dx25*(u(icrm,ic  ,j,k)+u(icrm,i  ,j,k))*(v(icrm,i  ,j,k)+v(icrm,ic  ,j,k))
            dvdt(icrm,i,j,k,na)  = dvdt(icrm,i,j,k,na)  - (fv2-fv1)
            fw1 = dx25*(u(icrm,ic-1,j,k)*rho(icrm,k)*adz(icrm,k)+u(icrm,ic-1,j,kcu)*rho(icrm,kcu)*adz(icrm,kcu))*(w(icrm,i-1,j,kc)+w(icrm,ic-1,j,kc))
            fw2 = dx25*(u(icrm,ic  ,j,k)*rho(icrm,k)*adz(icrm,k)+u(icrm,ic  ,j,kcu)*rho(icrm,kcu)*adz(icrm,kcu))*(w(icrm,i  ,j,kc)+w(icrm,ic  ,j,kc))
            dwdt(icrm,i,j,kc,na) = dwdt(icrm,i,j,kc,na)-irho*(fw2-fw1)
          end do
        end do ! k
      end do

    endif

  end subroutine advect2_mom_xy

end module advect2_mom_xy_mod
