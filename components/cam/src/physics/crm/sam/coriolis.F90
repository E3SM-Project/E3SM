module coriolis_mod
  use params, only: asyncid
  implicit none

contains

  subroutine coriolis(ncrms)
    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) u_av, v_av, w_av
    integer i,j,k,ib,ic,jb,jc,kc,icrm

    if(RUN3D) then
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              kc=k+1
              jb=j-1
              jc=j+1
              ib=i-1
              ic=i+1
              v_av=0.25*(v(icrm,i,j,k)+v(icrm,i,jc,k)+v(icrm,ib,j,k)+v(icrm,ib,jc,k))
              w_av=0.25*(w(icrm,i,j,kc)+w(icrm,ib,j,kc)+w(icrm,i,j,k)+w(icrm,ib,j,k))
              dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)+fcory(icrm,j)*(v_av-vg0(icrm,k))-fcorzy(icrm,j)*w_av
              u_av=0.25*(u(icrm,i,j,k)+u(icrm,ic,j,k)+u(icrm,i,jb,k)+u(icrm,ic,jb,k))
              dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-0.5*(fcory(icrm,j)+fcory(icrm,jb))*(u_av-ug0(icrm,k))
            end do ! i
          end do ! j
        end do ! k
      end do
    else
      !$acc parallel loop collapse(4) async(asyncid)
      do k=1,nzm
        do j=1,ny
          do i=1,nx
            do icrm = 1 , ncrms
              kc=k+1
              ib=i-1
              ic=i+1
              w_av=0.25*(w(icrm,i,j,kc)+w(icrm,ib,j,kc)+w(icrm,i,j,k)+w(icrm,ib,j,k))
              dudt(icrm,i,j,k,na)=dudt(icrm,i,j,k,na)+fcory(icrm,j)*(v(icrm,i,j,k)-vg0(icrm,k))-fcorzy(icrm,j)*w_av
              dvdt(icrm,i,j,k,na)=dvdt(icrm,i,j,k,na)-fcory(icrm,j)*(u(icrm,i,j,k)-ug0(icrm,k))
            end do ! i
          end do ! i
        end do ! k
      end do
    endif

  end subroutine coriolis

end module coriolis_mod
