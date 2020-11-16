module shear_prod2D_mod
  use params, only: asyncid
  implicit none

contains

  subroutine shear_prod2D(ncrms,def2)
    use vars
    use params, only: crm_rknd
    implicit none
    integer, intent(in) :: ncrms
    real(crm_rknd) def2(ncrms,nx,ny,nzm)
    real(crm_rknd) rdx0,rdx,rdx_up,rdx_dn
    real(crm_rknd) rdz,rdzw_up,rdzw_dn
    integer i,j,k,ib,ic,kb,kc,icrm

    rdx0=1./dx
    j=1

    !$acc parallel loop collapse(3) async(asyncid)
    do k=1,nzm
      do i=1,nx
        do icrm = 1 , ncrms
          if ( k >= 2 .and. k <= nzm-1) then
            kb=k-1
            kc=k+1
            rdz = 1./(dz(icrm)*adz(icrm,k))
            rdzw_up = 1./(dz(icrm)*adzw(icrm,kc))
            rdzw_dn = 1./(dz(icrm)*adzw(icrm,k))
            rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
            rdx_up=rdx0 * sqrt(dx*rdzw_up)
            rdx_dn=rdx0 * sqrt(dx*rdzw_dn)
            ib=i-1
            ic=i+1
            def2(icrm,i,j,k)=2.* ( &
            ( (u(icrm,ic,j,k)-u(icrm,i,j,k))*rdx)**2+ &
            ( (w(icrm,i,j,kc)-w(icrm,i,j,k))*rdz)**2 ) &
            + 0.5 * ( &
            ( (v(icrm,ic,j ,k)-v(icrm,i ,j ,k))*rdx )**2 +  &
            ( (v(icrm,i ,j ,k)-v(icrm,ib,j ,k))*rdx )**2 +   &
            ( (u(icrm,ic,j,kc)-u0(icrm,kc)-u(icrm,ic,j, k)+u0(icrm,k))*rdzw_up+ &
            (w(icrm,ic,j,kc)-w(icrm,i ,j,kc))*rdx_up )**2 + &
            ( (u(icrm,i ,j,kc)-u0(icrm,kc)-u(icrm,i ,j, k)+u0(icrm,k))*rdzw_up+ &
            (w(icrm,i ,j,kc)-w(icrm,ib,j,kc))*rdx_up )**2 + &
            ( (u(icrm,ic,j,k )-u0(icrm,k)-u(icrm,ic,j,kb)+u0(icrm,kb))*rdzw_dn+ &
            (w(icrm,ic,j,k )-w(icrm,i ,j,k ))*rdx_dn )**2 + &
            ( (u(icrm,i ,j,k )-u0(icrm,k)-u(icrm,i ,j,kb)+u0(icrm,kb))*rdzw_dn+ &
            (w(icrm,i ,j,k )-w(icrm,ib,j,k ))*rdx_dn )**2 + &
            ( (v(icrm,i,j ,kc)-v0(icrm,kc)-v(icrm,i,j , k)+v0(icrm,k))*rdzw_up )**2 + &
            ( (v(icrm,i,j ,k )-v0(icrm,k)-v(icrm,i,j ,kb)+v0(icrm,kb))*rdzw_dn )**2 )
          elseif (k == 1) then
            kc=k+1
            rdz = 1./(dz(icrm)*adz(icrm,k))
            rdzw_up = 1./(dz(icrm)*adzw(icrm,kc))
            rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
            rdx_up=rdx0 * sqrt(dx*rdzw_up)
            ib=i-1
            ic=i+1
            def2(icrm,i,j,k)=2.* ( &
            ( (u(icrm,ic,j,k)-u(icrm,i,j,k))*rdx)**2+ &
            ( (w(icrm,i,j,kc)-w(icrm,i,j,k))*rdz)**2 ) &
            + 0.5 * ( &
            ( (v(icrm,ic,j ,k)-v(icrm,i ,j ,k))*rdx )**2 + &
            ( (v(icrm,i ,j ,k)-v(icrm,ib,j ,k))*rdx )**2 ) &
            +( (v(icrm,i,j ,kc)-v0(icrm,kc)-v(icrm,i,j,k)+v0(icrm,k))*rdzw_up )**2 &
            + 0.5 * ( &
            ( (u(icrm,ic,j,kc)-u0(icrm,kc)-u(icrm,ic,j, k)+u0(icrm,k))*rdzw_up+ &
            (w(icrm,ic,j,kc)-w(icrm,i ,j,kc))*rdx_up )**2 + &
            ( (u(icrm,i ,j,kc)-u0(icrm,kc)-u(icrm,i ,j, k)+u0(icrm,k))*rdzw_up+ &
            (w(icrm,i ,j,kc)-w(icrm,ib,j,kc))*rdx_up )**2 )
          elseif (k == nzm) then
            kc=k+1
            kb=k-1
            rdz = 1./(dz(icrm)*adz(icrm,k))
            rdzw_dn = 1./(dz(icrm)*adzw(icrm,k))
            rdx=rdx0 * sqrt(dx*rdz) ! take into account grid anisotropy
            rdx_dn=rdx0 * sqrt(dx*rdzw_dn)
            ib=i-1
            ic=i+1
            def2(icrm,i,j,k)=2.* ( &
            ( (u(icrm,ic,j,k)-u(icrm,i,j,k))*rdx)**2+ &
            ( (w(icrm,i,j,kc)-w(icrm,i,j,k))*rdz)**2 ) &
            + 0.5 * ( &
            ( (v(icrm,ic,j ,k)-v(icrm,i ,j ,k))*rdx )**2 +  &
            ( (v(icrm,i ,j ,k)-v(icrm,ib,j ,k))*rdx )**2 )   &
            + ( (v(icrm,i,j ,k )-v0(icrm,k)-v(icrm,i,j ,kb)+v0(icrm,kb))*rdzw_dn )**2 &
            + 0.5 * ( &
            ( (u(icrm,ic,j,k )-u0(icrm,k)-u(icrm,ic,j,kb)+u0(icrm,kb))*rdzw_dn+ &
            (w(icrm,ic,j,k )-w(icrm,i ,j,k ))*rdx_dn )**2 + &
            ( (u(icrm,i ,j,k )-u0(icrm,k)-u(icrm,i ,j,kb)+u0(icrm,kb))*rdzw_dn+ &
            (w(icrm,i ,j,k )-w(icrm,ib,j,k ))*rdx_dn )**2 )
          endif
        end do
      end do ! k
    end do
  end subroutine shear_prod2D

end module shear_prod2D_mod
