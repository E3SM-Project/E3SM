
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module derivative_mod
  use derivative_mod_base, only: derivative_t, subcell_integration, subcell_dss_fluxes, subcell_div_fluxes, subcell_Laplace_fluxes, allocate_subcell_integration_matrix,   &
                                 derivinit, &
                                 gradient_sphere_wk_testcov, gradient_sphere_wk_testcontra, ugradv_sphere, vorticity_sphere, vorticity_sphere_diag, curl_sphere,     &
                                 curl_sphere_wk_testcov, vlaplace_sphere_wk, element_boundary_integral, edge_flux_u_cg, limiter_optim_iter_full, &
                                 laplace_sphere_wk, divergence_sphere_wk, gradient_sphere, divergence_sphere
  use kinds, only : real_kind, longdouble_kind
  use dimensions_mod, only : np, nelemd, nlev
  use quadrature_mod, only : quadrature_t, gauss, gausslobatto,legendre, jacobi
  use parallel_mod, only : abortmp
  ! needed for spherical differential operators:
  use physical_constants, only : rrearth 
  use element_mod, only : element_t
  use control_mod, only : hypervis_scaling, hypervis_power
  implicit none
  private

  public ::  derivative_t, subcell_integration, subcell_dss_fluxes, subcell_div_fluxes, subcell_Laplace_fluxes, allocate_subcell_integration_matrix,   &
             derivinit, &
             gradient_sphere_wk_testcov, gradient_sphere_wk_testcontra, ugradv_sphere, vorticity_sphere, vorticity_sphere_diag, curl_sphere,     &
             curl_sphere_wk_testcov, vlaplace_sphere_wk, element_boundary_integral, edge_flux_u_cg, limiter_optim_iter_full, &
             laplace_sphere_wk, divergence_sphere_wk, gradient_sphere, divergence_sphere
  public :: laplace_sphere_wk_openacc
  public :: divergence_sphere_wk_openacc
  public :: gradient_sphere_openacc
  public :: divergence_sphere_openacc

contains

  subroutine laplace_sphere_wk_openacc(s,grads,deriv,elem,var_coef,laplace,len,nets,nete,ntl,tl)
    use element_mod, only: element_t
    use control_mod, only: hypervis_scaling, hypervis_power
    implicit none
    !input:  s = scalar
    !ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
    !note: for this form of the operator, grad(s) does not need to be made C0
    real(kind=real_kind) , intent(in   ) :: s(np,np,len,ntl,nelemd)
    real(kind=real_kind) , intent(inout) :: grads(np,np,2,len,nelemd)
    type (derivative_t)  , intent(in   ) :: deriv
    type (element_t)     , intent(in   ) :: elem(:)
    logical              , intent(in   ) :: var_coef
    real(kind=real_kind) , intent(  out) :: laplace(np,np,len,ntl,nelemd)
    integer              , intent(in   ) :: len,nets,nete,ntl,tl
    integer :: i,j,k,ie
    ! Local
    real(kind=real_kind) :: oldgrads(2)
    call gradient_sphere_openacc(s,deriv,elem(:),grads,len,nets,nete,ntl,tl)
    !$acc parallel loop gang vector collapse(4) present(grads,elem(:)) private(oldgrads)
    do ie = nets , nete
      do k = 1 , len
        do j = 1 , np
          do i = 1 , np
            if (var_coef) then
              if (hypervis_power/=0 ) then
                ! scalar viscosity with variable coefficient
                grads(i,j,1,k,ie) = grads(i,j,1,k,ie)*elem(ie)%variable_hyperviscosity(i,j)
                grads(i,j,2,k,ie) = grads(i,j,2,k,ie)*elem(ie)%variable_hyperviscosity(i,j)
              else if (hypervis_scaling /=0 ) then
                oldgrads = grads(i,j,:,k,ie)
                grads(i,j,1,k,ie) = sum(oldgrads(:)*elem(ie)%tensorVisc(i,j,1,:))
                grads(i,j,2,k,ie) = sum(oldgrads(:)*elem(ie)%tensorVisc(i,j,2,:))
              endif
            endif
          enddo
        enddo
      enddo
    enddo
    ! note: divergnece_sphere and divergence_sphere_wk are identical *after* bndry_exchange
    ! if input is C_0.  Here input is not C_0, so we should use divergence_sphere_wk().  
    call divergence_sphere_wk_openacc(grads,deriv,elem(:),laplace,len,nets,nete,ntl,tl)
  end subroutine laplace_sphere_wk_openacc

  subroutine divergence_sphere_wk_openacc(v,deriv,elem,div,len,nets,nete,ntl,tl)
    use element_mod, only: element_t
    use physical_constants, only: rrearth
    implicit none
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v, integrated by parts
!   Computes  -< grad(psi) dot v > 
!   (the integrated by parts version of < psi div(v) > )
!   note: after DSS, divergence_sphere () and divergence_sphere_wk() 
!   are identical to roundoff, as theory predicts.
    real(kind=real_kind), intent(in) :: v(np,np,2,len,ntl,nelemd)  ! in lat-lon coordinates
    type (derivative_t) , intent(in) :: deriv
    type (element_t)    , intent(in) :: elem(:)
    real(kind=real_kind), intent(out):: div(np,np,len,ntl,nelemd)
    integer             , intent(in) :: len
    integer             , intent(in) :: nets , nete , ntl , tl
    ! Local
    integer, parameter :: kchunk = 8
    integer :: i,j,l,k,ie,kc,kk
    real(kind=real_kind) :: vtemp(np,np,2,kchunk), tmp, deriv_tmp(np,np)
    ! latlon- > contra
    !$acc parallel loop gang collapse(2) present(v,elem(:),div,deriv) private(vtemp,deriv_tmp)
    do ie = nets , nete
      do kc = 1 , len/kchunk+1
        !$acc cache(vtemp,deriv_tmp)
        !$acc loop vector collapse(3)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= len) then
                vtemp(i,j,1,kk)=elem(ie)%spheremp(i,j)*(elem(ie)%Dinv(i,j,1,1)*v(i,j,1,k,tl,ie) + elem(ie)%Dinv(i,j,1,2)*v(i,j,2,k,tl,ie))
                vtemp(i,j,2,kk)=elem(ie)%spheremp(i,j)*(elem(ie)%Dinv(i,j,2,1)*v(i,j,1,k,tl,ie) + elem(ie)%Dinv(i,j,2,2)*v(i,j,2,k,tl,ie))
              endif
              if (kk == 1) deriv_tmp(i,j) = deriv%Dvv(i,j)
            enddo
          enddo
        enddo
        !$acc loop vector collapse(3) private(tmp)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= len) then
                tmp = 0.
                do l = 1 , np
                  tmp = tmp - ( vtemp(l,j,1,kk)*deriv_tmp(i,l) + vtemp(i,l,2,kk)*deriv_tmp(j,l) )
                enddo
                div(i,j,k,tl,ie) = tmp * rrearth
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine divergence_sphere_wk_openacc

  subroutine gradient_sphere_openacc(s,deriv,elem,ds,len,nets,nete,ntl,tl)
    use element_mod, only: element_t
    use physical_constants, only: rrearth
    implicit none
    !   input s:  scalar
    !   output  ds: spherical gradient of s, lat-lon coordinates
    real(kind=real_kind), intent(in) :: s(np,np,len,ntl,nelemd)
    type(derivative_t)  , intent(in) :: deriv
    type(element_t)     , intent(in) :: elem(:)
    real(kind=real_kind), intent(out):: ds(np,np,2,len,ntl,nelemd)
    integer             , intent(in) :: len
    integer             , intent(in) :: nets,nete,ntl,tl
    integer, parameter :: kchunk = 8
    integer :: i, j, l, k, ie, kc, kk
    real(kind=real_kind) :: dsdx00, dsdy00
    real(kind=real_kind) :: stmp(np,np,kchunk), deriv_tmp(np,np)
    !$acc parallel loop gang collapse(2) present(ds,elem(:),s,deriv%Dvv) private(stmp,deriv_tmp)
    do ie = nets , nete
      do kc = 1 , len/kchunk+1
        !$acc cache(stmp,deriv_tmp)
        !$acc loop vector collapse(3)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k > len) k = len
              stmp(i,j,kk) = s(i,j,k,tl,ie)
              if (kk == 1) deriv_tmp(i,j) = deriv%Dvv(i,j)
            enddo
          enddo
        enddo
        !$acc loop vector collapse(3)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= len) then
                dsdx00=0.0d0
                dsdy00=0.0d0
                do l = 1 , np
                  dsdx00 = dsdx00 + deriv_tmp(l,i)*stmp(l,j,kk)
                  dsdy00 = dsdy00 + deriv_tmp(l,j)*stmp(i,l,kk)
                enddo
                ds(i,j,1,k,tl,ie) = ( elem(ie)%Dinv(i,j,1,1)*dsdx00 + elem(ie)%Dinv(i,j,2,1)*dsdy00 ) * rrearth
                ds(i,j,2,k,tl,ie) = ( elem(ie)%Dinv(i,j,1,2)*dsdx00 + elem(ie)%Dinv(i,j,2,2)*dsdy00 ) * rrearth
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine gradient_sphere_openacc

  subroutine divergence_sphere_openacc(v,deriv,elem,div,len,nets,nete,ntl,tl)
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v
    use element_mod   , only: element_t
    use physical_constants, only: rrearth
    implicit none
    real(kind=real_kind), intent(in   ) :: v(np,np,2,len,ntl,nelemd)  ! in lat-lon coordinates
    type(derivative_t)  , intent(in   ) :: deriv
    type(element_t)     , intent(in   ) :: elem(:)
    real(kind=real_kind), intent(  out) :: div(np,np,len,ntl,nelemd)
    integer             , intent(in   ) :: len , nets , nete , ntl , tl
    ! Local
    integer, parameter :: kchunk = 8
    integer :: i, j, l, k, ie, kc, kk
    real(kind=real_kind) ::  dudx00, dvdy00, gv(np,np,kchunk,2), deriv_tmp(np,np)
    ! convert to contra variant form and multiply by g
    !$acc parallel loop gang collapse(2) private(gv,deriv_tmp) present(v,deriv,div,elem(:))
    do ie = nets , nete
      do kc = 1 , len/kchunk+1
        !$acc cache(gv,deriv_tmp)
        !$acc loop vector collapse(3) private(k)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= len) then
                gv(i,j,kk,1)=elem(ie)%metdet(i,j)*(elem(ie)%Dinv(i,j,1,1)*v(i,j,1,k,tl,ie) + elem(ie)%Dinv(i,j,1,2)*v(i,j,2,k,tl,ie))
                gv(i,j,kk,2)=elem(ie)%metdet(i,j)*(elem(ie)%Dinv(i,j,2,1)*v(i,j,1,k,tl,ie) + elem(ie)%Dinv(i,j,2,2)*v(i,j,2,k,tl,ie))
              endif
              if (kk == 1) deriv_tmp(i,j) = deriv%dvv(i,j)
            enddo
          enddo
        enddo
        ! compute d/dx and d/dy         
        !$acc loop vector collapse(3) private(dudx00,dvdy00,k)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= len) then
                dudx00=0.0d0
                dvdy00=0.0d0
                do l = 1 , np
                  dudx00 = dudx00 + deriv_tmp(l,i)*gv(l,j,kk,1)
                  dvdy00 = dvdy00 + deriv_tmp(l,j)*gv(i,l,kk,2)
                enddo
                div(i,j,k,tl,ie)=(dudx00+dvdy00)*(elem(ie)%rmetdet(i,j)*rrearth)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine divergence_sphere_openacc

end module derivative_mod

