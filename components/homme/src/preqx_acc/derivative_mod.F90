
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module derivative_mod
  use derivative_mod_base, only: derivative_t, subcell_integration, subcell_dss_fluxes, subcell_div_fluxes, subcell_Laplace_fluxes, allocate_subcell_integration_matrix,   &
                                 derivinit, &
                                 gradient_sphere_wk_testcov, gradient_sphere_wk_testcontra, ugradv_sphere, vorticity_sphere, vorticity_sphere_diag, curl_sphere,     &
                                 curl_sphere_wk_testcov, vlaplace_sphere_wk, element_boundary_integral, edge_flux_u_cg, limiter_optim_iter_full, &
                                 laplace_sphere_wk, divergence_sphere_wk, gradient_sphere, divergence_sphere, laplace_z
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
             laplace_sphere_wk, divergence_sphere_wk, gradient_sphere, divergence_sphere, laplace_z
  public :: laplace_sphere_wk_openacc
  public :: divergence_sphere_wk_openacc
  public :: gradient_sphere_openacc
  public :: divergence_sphere_openacc
  public :: vorticity_sphere_openacc
  public :: vlaplace_sphere_wk_openacc

contains


  subroutine vlaplace_sphere_wk_openacc(v,vor,div,deriv,elem,var_coef,len,nets,nete,ntl_in,tl_in,ntl_out,tl_out,laplace,nu_ratio,klim_in)
    !   input:  v = vector in lat-lon coordinates
    !   ouput:  weak laplacian of v, in lat-lon coordinates
    !   logic:
    !      tensorHV:     requires cartesian
    !      nu_div/=nu:   requires contra formulatino
    !   One combination NOT supported:  tensorHV and nu_div/=nu then abort
    implicit none
    real(kind=real_kind), intent(in   ) :: v(np,np,2,len,ntl_in,nets:nete)
    real(kind=real_kind), intent(  out) :: vor(np,np,len,nets:nete)
    real(kind=real_kind), intent(  out) :: div(np,np,len,nets:nete)
    logical             , intent(in   ) :: var_coef
    type (derivative_t) , intent(in   ) :: deriv
    type (element_t)    , intent(in   ) :: elem(:)
    real(kind=real_kind), intent(in   ) :: nu_ratio
    integer             , intent(in   ) :: len,nets,nete,ntl_in,tl_in,ntl_out,tl_out
    real(kind=real_kind), intent(  out) :: laplace(np,np,2,len,ntl_out,nets:nete)
    integer, optional   , intent(in   ) :: klim_in
    if (hypervis_scaling/=0 .and. var_coef) then
      call abortmp('hypervis_scaling/=0 .and. var_coef not supported in OpenACC!')
    else
      ! all other cases, use contra formulation:
      call vlaplace_sphere_wk_contra(v,vor,div,deriv,elem,var_coef,len,nets,nete,ntl_in,tl_in,ntl_out,tl_out,laplace,nu_ratio,klim_in)
    endif
  end subroutine vlaplace_sphere_wk_openacc

  subroutine vlaplace_sphere_wk_contra(v,vor,div,deriv,elem,var_coef,len,nets,nete,ntl_in,tl_in,ntl_out,tl_out,laplace,nu_ratio,klim_in)
    !   input:  v = vector in lat-lon coordinates
    !   ouput:  weak laplacian of v, in lat-lon coordinates
    implicit none
    real(kind=real_kind), intent(in   ) :: v(np,np,2,len,ntl_in,nets:nete)
    real(kind=real_kind), intent(  out) :: vor(np,np,len,nets:nete)
    real(kind=real_kind), intent(  out) :: div(np,np,len,nets:nete)
    logical             , intent(in   ) :: var_coef
    type (derivative_t) , intent(in   ) :: deriv
    type (element_t)    , intent(in   ) :: elem(:)
    real(kind=real_kind), intent(in   ) :: nu_ratio
    integer             , intent(in   ) :: len,nets,nete,ntl_in,tl_in,ntl_out,tl_out
    real(kind=real_kind), intent(  out) :: laplace(np,np,2,len,ntl_out,nets:nete)
    integer, optional   , intent(in   ) :: klim_in
    ! Local
    integer :: i,j,l,m,n,k,ie,klim
    klim = len
    if (present(klim_in)) klim = klim_in
    call divergence_sphere_openacc(v,deriv,elem,div,len,nets,nete,ntl_in,tl_in,1,1,klim)
    call vorticity_sphere_openacc (v,deriv,elem,vor,len,nets,nete,ntl_in,tl_in,1,1,klim)
    !$acc parallel loop gang vector collapse(4) present(div,vor,elem)
    do ie = nets , nete
      do k = 1 , klim
        do j = 1 , np
          do i = 1 , np
            if (var_coef .and. hypervis_power/=0 ) then
              ! scalar viscosity with variable coefficient
              div(i,j,k,ie) = div(i,j,k,ie)*elem(ie)%variable_hyperviscosity(i,j)
              vor(i,j,k,ie) = vor(i,j,k,ie)*elem(ie)%variable_hyperviscosity(i,j)
            endif
            div(i,j,k,ie) = nu_ratio*div(i,j,k,ie)
          enddo
        enddo
      enddo
    enddo
    call gradient_minus_curl_sphere_wk_testcov_openacc(div,vor,deriv,elem,len,nets,nete,1,1,ntl_out,tl_out,laplace,klim)
    !$acc parallel loop gang vector collapse(4) present(laplace,elem,v)
    do ie = nets , nete
      do k = 1 , klim
        do n=1,np
          do m=1,np
            ! add in correction so we dont damp rigid rotation
            laplace(m,n,1,k,tl_out,ie)=laplace(m,n,1,k,tl_out,ie) + 2*elem(ie)%spheremp(m,n)*v(m,n,1,k,tl_in,ie)*(rrearth**2)
            laplace(m,n,2,k,tl_out,ie)=laplace(m,n,2,k,tl_out,ie) + 2*elem(ie)%spheremp(m,n)*v(m,n,2,k,tl_in,ie)*(rrearth**2)
          enddo
        enddo
      enddo
    enddo
  end subroutine vlaplace_sphere_wk_contra

  !TODO: make this efficient with shared memory!
  subroutine gradient_minus_curl_sphere_wk_testcov_openacc(s1,s2,deriv,elem,len,nets,nete,ntl_in,tl_in,ntl_out,tl_out,ds,klim_in)
    !   integrated-by-parts gradient, w.r.t. COVARIANT test functions
    !   input s:  scalar
    !   output  ds: weak gradient, lat/lon coordinates
    implicit none
    type (derivative_t) , intent(in   ) :: deriv
    type (element_t)    , intent(in   ) :: elem(:)
    real(kind=real_kind), intent(in   ) :: s1(np,np  ,len,ntl_in ,nets:nete)
    real(kind=real_kind), intent(in   ) :: s2(np,np  ,len,ntl_in ,nets:nete)
    real(kind=real_kind), intent(  out) :: ds(np,np,2,len,ntl_out,nets:nete)
    integer             , intent(in   ) :: len,nets,nete,ntl_in,tl_in,ntl_out,tl_out
    integer, optional   , intent(in   ) :: klim_in
    integer :: i,j,l,k,ie,klim
    real(kind=real_kind) :: dscontra1, dscontra2
    klim = len
    if (present(klim_in)) klim = klim_in
    !$acc parallel loop gang vector collapse(4) private(dscontra1, dscontra2) present(elem,deriv,s1,s2,ds)
    do ie=nets,nete
      do k=1,klim
        do j=1,np
          do i=1,np
            dscontra1 = 0
            dscontra2 = 0
            do l=1,np
              dscontra1=dscontra1-( (elem(ie)%mp(l,j)*elem(ie)%metinv(i,j,1,1)*elem(ie)%metdet(i,j)*s1(l,j,k,tl_in,ie)*deriv%Dvv(i,l)) + &
                                    (elem(ie)%mp(i,l)*elem(ie)%metinv(i,j,2,1)*elem(ie)%metdet(i,j)*s1(i,l,k,tl_in,ie)*deriv%Dvv(j,l)) - &
                                    (elem(ie)%mp(i,l)*s2(i,l,k,tl_in,ie)*deriv%Dvv(j,l)) ) *rrearth
              dscontra2=dscontra2-( (elem(ie)%mp(l,j)*elem(ie)%metinv(i,j,1,2)*elem(ie)%metdet(i,j)*s1(l,j,k,tl_in,ie)*deriv%Dvv(i,l)) + &
                                    (elem(ie)%mp(i,l)*elem(ie)%metinv(i,j,2,2)*elem(ie)%metdet(i,j)*s1(i,l,k,tl_in,ie)*deriv%Dvv(j,l)) + &
                                    (elem(ie)%mp(l,j)*s2(l,j,k,tl_in,ie)*deriv%Dvv(i,l)) ) *rrearth
            enddo
            ds(i,j,1,k,tl_out,ie)=(elem(ie)%D(i,j,1,1)*dscontra1 + elem(ie)%D(i,j,1,2)*dscontra2)
            ds(i,j,2,k,tl_out,ie)=(elem(ie)%D(i,j,2,1)*dscontra1 + elem(ie)%D(i,j,2,2)*dscontra2)
          enddo
        enddo
      enddo
    enddo
  end subroutine gradient_minus_curl_sphere_wk_testcov_openacc


  subroutine vorticity_sphere_openacc(v,deriv,elem,vort,len,nets,nete,ntl_in,tl_in,ntl_out,tl_out,klim_in)
    implicit none
!   input:  v = velocity in lat-lon coordinates
!   ouput:  spherical vorticity of v
    type (derivative_t) , intent(in   ) :: deriv
    type (element_t)    , intent(in   ) :: elem(:)
    real(kind=real_kind), intent(in   ) :: v(np,np,2,len,ntl_in,nets:nete)
    real(kind=real_kind), intent(  out) :: vort(np,np,len,ntl_out,nets:nete)
    integer             , intent(in   ) :: len,nets,nete,ntl_in,tl_in,ntl_out,tl_out
    integer, optional   , intent(in   ) :: klim_in
    integer, parameter :: kchunk = 8
    integer :: i, j, l, k, ie, kc, kk, klim
    real(kind=real_kind) :: dvdx00,dudy00
    real(kind=real_kind) :: vco(np,np,2,kchunk)
    klim = len
    if (present(klim_in)) klim = klim_in
    ! convert to covariant form
    !$acc parallel loop gang collapse(2) private(vco) present(elem,v,vort,deriv)
    do ie = nets , nete
      do kc = 1 , len/kchunk+1
        !$acc cache(vco)
        !$acc loop vector collapse(3) private(k)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= klim) then
                vco(i,j,1,kk)=(elem(ie)%D(i,j,1,1)*v(i,j,1,k,tl_in,ie) + elem(ie)%D(i,j,2,1)*v(i,j,2,k,tl_in,ie))
                vco(i,j,2,kk)=(elem(ie)%D(i,j,1,2)*v(i,j,1,k,tl_in,ie) + elem(ie)%D(i,j,2,2)*v(i,j,2,k,tl_in,ie))
              endif
            enddo
          enddo
        enddo
        !$acc loop vector collapse(3) private(k,dudy00,dvdx00)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= klim) then
                dudy00=0.0d0
                dvdx00=0.0d0
                do l = 1 , np
                  dvdx00 = dvdx00 + deriv%Dvv(l,i)*vco(l,j,2,kk)
                  dudy00 = dudy00 + deriv%Dvv(l,j)*vco(i,l,1,kk)
                enddo
                vort(i,j,k,tl_out,ie)=(dvdx00-dudy00)*(elem(ie)%rmetdet(i,j)*rrearth)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine vorticity_sphere_openacc

  subroutine laplace_sphere_wk_openacc(s,grads,deriv,elem,var_coef,laplace,len,nets,nete,ntl_in,tl_in,ntl_out,tl_out,klim_in)
    use element_mod, only: element_t
    use control_mod, only: hypervis_scaling, hypervis_power
    implicit none
    !input:  s = scalar
    !ouput:  -< grad(PHI), grad(s) >   = weak divergence of grad(s)
    !note: for this form of the operator, grad(s) does not need to be made C0
    real(kind=real_kind) , intent(in   ) :: s(np,np,len,ntl_in,nelemd)
    real(kind=real_kind) , intent(  out) :: grads(np,np,2,len,nelemd)
    type (derivative_t)  , intent(in   ) :: deriv
    type (element_t)     , intent(in   ) :: elem(:)
    logical              , intent(in   ) :: var_coef
    real(kind=real_kind) , intent(  out) :: laplace(np,np,len,ntl_out,nelemd)
    integer              , intent(in   ) :: len,nets,nete,ntl_in,tl_in,ntl_out,tl_out
    integer, optional    , intent(in   ) :: klim_in
    integer :: i,j,k,ie,klim
    ! Local
    real(kind=real_kind) :: oldgrads(2)
    klim = len
    if (present(klim_in)) klim = klim_in
    call gradient_sphere_openacc(s,deriv,elem(:),grads,len,nets,nete,ntl_in,tl_in,1,1,klim)
    !$acc parallel loop gang vector collapse(4) private(oldgrads)
    do ie = nets , nete
      do k = 1 , klim
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
    call divergence_sphere_wk_openacc(grads,deriv,elem(:),laplace,len,nets,nete,1,1,ntl_out,tl_out,klim)
  end subroutine laplace_sphere_wk_openacc

  subroutine divergence_sphere_wk_openacc(v,deriv,elem,div,len,nets,nete,ntl_in,tl_in,ntl_out,tl_out,klim_in)
    use element_mod, only: element_t
    use physical_constants, only: rrearth
    implicit none
    !   input:  v = velocity in lat-lon coordinates
    !   ouput:  div(v)  spherical divergence of v, integrated by parts
    !   Computes  -< grad(psi) dot v >
    !   (the integrated by parts version of < psi div(v) > )
    !   note: after DSS, divergence_sphere () and divergence_sphere_wk()
    !   are identical to roundoff, as theory predicts.
    real(kind=real_kind), intent(in) :: v(np,np,2,len,ntl_in,nelemd)  ! in lat-lon coordinates
    type (derivative_t) , intent(in) :: deriv
    type (element_t)    , intent(in) :: elem(:)
    real(kind=real_kind), intent(out):: div(np,np,len,ntl_out,nelemd)
    integer             , intent(in) :: len
    integer             , intent(in) :: nets , nete , ntl_in , tl_in , ntl_out , tl_out
    integer, optional    , intent(in   ) :: klim_in
    ! Local
    integer, parameter :: kchunk = 8
    integer :: i,j,l,k,ie,kc,kk,klim
    real(kind=real_kind) :: vtemp(np,np,2,kchunk), tmp, deriv_tmp(np,np)
    klim = len
    if (present(klim_in)) klim = klim_in
    ! latlon- > contra
    !$acc parallel loop gang collapse(2) private(vtemp,deriv_tmp)
    do ie = nets , nete
      do kc = 1 , len/kchunk+1
        !$acc cache(vtemp,deriv_tmp)
        !$acc loop vector collapse(3)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= klim) then
                vtemp(i,j,1,kk)=elem(ie)%spheremp(i,j)*(elem(ie)%Dinv(i,j,1,1)*v(i,j,1,k,tl_in,ie) + elem(ie)%Dinv(i,j,1,2)*v(i,j,2,k,tl_in,ie))
                vtemp(i,j,2,kk)=elem(ie)%spheremp(i,j)*(elem(ie)%Dinv(i,j,2,1)*v(i,j,1,k,tl_in,ie) + elem(ie)%Dinv(i,j,2,2)*v(i,j,2,k,tl_in,ie))
                if (kk == 1) deriv_tmp(i,j) = deriv%Dvv(i,j)
              endif
            enddo
          enddo
        enddo
        !$acc loop vector collapse(3) private(tmp)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= klim) then
                tmp = 0.
                do l = 1 , np
                  tmp = tmp - ( vtemp(l,j,1,kk)*deriv_tmp(i,l) + vtemp(i,l,2,kk)*deriv_tmp(j,l) )
                enddo
                div(i,j,k,tl_out,ie) = tmp * rrearth
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine divergence_sphere_wk_openacc

  subroutine gradient_sphere_openacc(s,deriv,elem,ds,len,nets,nete,ntl_in,tl_in,ntl_out,tl_out,klim_in)
    use element_mod, only: element_t
    use physical_constants, only: rrearth
    implicit none
    !   input s:  scalar
    !   output  ds: spherical gradient of s, lat-lon coordinates
    real(kind=real_kind), intent(in) :: s(np,np,len,ntl_in,nelemd)
    type(derivative_t)  , intent(in) :: deriv
    type(element_t)     , intent(in) :: elem(:)
    real(kind=real_kind), intent(out):: ds(np,np,2,len,ntl_out,nelemd)
    integer             , intent(in) :: len
    integer             , intent(in) :: nets,nete,ntl_in,ntl_out,tl_in,tl_out
    integer, optional    , intent(in   ) :: klim_in
    integer, parameter :: kchunk = 8
    integer :: i, j, l, k, ie, kc, kk,klim
    real(kind=real_kind) :: dsdx00, dsdy00
    real(kind=real_kind) :: stmp(np,np,kchunk), deriv_tmp(np,np)
    klim = len
    if (present(klim_in)) klim = klim_in
    !$acc parallel loop gang collapse(2) private(stmp,deriv_tmp)
    do ie = nets , nete
      do kc = 1 , len/kchunk+1
        !$acc cache(stmp,deriv_tmp)
        !$acc loop vector collapse(3)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k > klim) k = len
              stmp(i,j,kk) = s(i,j,k,tl_in,ie)
              if (kk == 1) deriv_tmp(i,j) = deriv%Dvv(i,j)
            enddo
          enddo
        enddo
        !$acc loop vector collapse(3)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= klim) then
                dsdx00=0.0d0
                dsdy00=0.0d0
                do l = 1 , np
                  dsdx00 = dsdx00 + deriv_tmp(l,i)*stmp(l,j,kk)
                  dsdy00 = dsdy00 + deriv_tmp(l,j)*stmp(i,l,kk)
                enddo
                ds(i,j,1,k,tl_out,ie) = ( elem(ie)%Dinv(i,j,1,1)*dsdx00 + elem(ie)%Dinv(i,j,2,1)*dsdy00 ) * rrearth
                ds(i,j,2,k,tl_out,ie) = ( elem(ie)%Dinv(i,j,1,2)*dsdx00 + elem(ie)%Dinv(i,j,2,2)*dsdy00 ) * rrearth
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine gradient_sphere_openacc

  subroutine divergence_sphere_openacc(v,deriv,elem,div,len,nets,nete,ntl_in,tl_in,ntl_out,tl_out,klim_in)
!   input:  v = velocity in lat-lon coordinates
!   ouput:  div(v)  spherical divergence of v
    use element_mod   , only: element_t
    use physical_constants, only: rrearth
    implicit none
    real(kind=real_kind), intent(in   ) :: v(np,np,2,len,ntl_in,nelemd)  ! in lat-lon coordinates
    type(derivative_t)  , intent(in   ) :: deriv
    type(element_t)     , intent(in   ) :: elem(:)
    real(kind=real_kind), intent(  out) :: div(np,np,len,ntl_out,nelemd)
    integer             , intent(in   ) :: len , nets , nete , ntl_in, ntl_out , tl_in , tl_out
    integer, optional   , intent(in   ) :: klim_in
    ! Local
    integer, parameter :: kchunk = 8
    integer :: i, j, l, k, ie, kc, kk, klim
    real(kind=real_kind) ::  dudx00, dvdy00, gv(np,np,kchunk,2), deriv_tmp(np,np)
    klim = len
    if (present(klim_in)) klim = klim_in
    ! convert to contra variant form and multiply by g
    !$acc parallel loop gang collapse(2) private(gv,deriv_tmp)
    do ie = nets , nete
      do kc = 1 , len/kchunk+1
        !$acc cache(gv,deriv_tmp)
        !$acc loop vector collapse(3) private(k)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= klim) then
                gv(i,j,kk,1)=elem(ie)%metdet(i,j)*(elem(ie)%Dinv(i,j,1,1)*v(i,j,1,k,tl_in,ie) + elem(ie)%Dinv(i,j,1,2)*v(i,j,2,k,tl_in,ie))
                gv(i,j,kk,2)=elem(ie)%metdet(i,j)*(elem(ie)%Dinv(i,j,2,1)*v(i,j,1,k,tl_in,ie) + elem(ie)%Dinv(i,j,2,2)*v(i,j,2,k,tl_in,ie))
                if (kk == 1) deriv_tmp(i,j) = deriv%dvv(i,j)
              endif
            enddo
          enddo
        enddo
        ! compute d/dx and d/dy
        !$acc loop vector collapse(3) private(dudx00,dvdy00,k)
        do kk = 1 , kchunk
          do j = 1 , np
            do i = 1 , np
              k = (kc-1)*kchunk+kk
              if (k <= klim) then
                dudx00=0.0d0
                dvdy00=0.0d0
                do l = 1 , np
                  dudx00 = dudx00 + deriv_tmp(l,i)*gv(l,j,kk,1)
                  dvdy00 = dvdy00 + deriv_tmp(l,j)*gv(i,l,kk,2)
                enddo
                div(i,j,k,tl_out,ie)=(dudx00+dvdy00)*(elem(ie)%rmetdet(i,j)*rrearth)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  end subroutine divergence_sphere_openacc

end module derivative_mod
