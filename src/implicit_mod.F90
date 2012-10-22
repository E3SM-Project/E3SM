#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


module implicit_mod

    use ,intrinsic :: iso_c_binding 

contains

  subroutine advance_imp_nonstag(elem, edge1, edge2, edge3, red, deriv, cg,   &
       hybrid, blkjac, lambdasq, dt, pmean, tl, nets, nete, xstate)

    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar
    use physical_constants, only : g
    use element_mod, only : element_t
    use parallel_mod, only : parallel_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk
    use time_mod, only : timelevel_t, smooth
    use control_mod, only :  topology, test_case
    use shallow_water_mod, only : tc1_velocity
    use bndry_mod, only : bndry_exchangev
    use cg_mod, only : cg_t
    use precon_mod, only : pcg_presolver_nonstag
    use solver_mod, only : blkjac_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use derived_type_mod ,only : derived_type, precon_type, initialize, init_precon
    use ,intrinsic :: iso_c_binding 
    use global_norms_mod
    use perf_mod, only : t_startf, t_stopf

    implicit none

    type (element_t) ,intent(inout)   :: elem(:)
    type (EdgeBuffer_t)  , intent(in) :: edge1
    type (EdgeBuffer_t)  , intent(in) :: edge2
    type (EdgeBuffer_t)  , intent(in) :: edge3
    type (reductionbuffer_ordered_1d_t)  ,intent(in) :: red
    type (derivative_t)  , intent(in) :: deriv
    type (cg_t)                       :: cg
    type (hybrid_t)      , intent(in) :: hybrid
    type (blkjac_t)      , intent(inout) :: blkjac(:)
    real (kind=real_kind), intent(inout) :: lambdasq(:)
    real (kind=real_kind), intent(in) :: dt
    real (kind=real_kind), intent(in) :: pmean
    type (timeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete

    ! =================
    ! Local
    ! =================

    ! Thread private working set ...

    type (parallel_t)          :: par
    real (c_double) ,allocatable ,dimension(:) :: xstate

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx, lx
    integer    :: comm
    integer    :: state_size

! state_object is a derived data type passed from doloca thru cpp to calc_f
! in the form of a pointer
    type(derived_type) ,target         :: state_object
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    type(precon_type) ,target          :: pre_object
    type(precon_type) ,pointer         :: pptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_pre

  interface 
   subroutine noxsolve(vectorSize,vector,v_container,p_container) &
     bind(C,name='noxsolve')
    use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
      integer(c_int)                :: vectorSize
      real(c_double)  ,dimension(*) :: vector
      type(c_ptr)                   :: v_container 
      type(c_ptr)                   :: p_container  !precon ptr
    end subroutine noxsolve
  end interface

    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    comm = par%comm
    call t_startf('implicit header')

! n+1 guess is result at n - need this to match up with time_mod bootstrap style
       do ie=nets,nete
             elem(ie)%state%v(:,:,1,:,np1)= elem(ie)%state%v(:,:,1,:,n0)
             elem(ie)%state%v(:,:,2,:,np1)= elem(ie)%state%v(:,:,2,:,n0)
             elem(ie)%state%p(:,:,:,np1)= elem(ie)%state%p(:,:,:,n0)
       end do !ie

       lx = 0
       do n=1,nvar
        do ie=nets,nete
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) xstate(lx) = elem(ie)%state%v(i,j,1,k,np1)
             if (n==2) xstate(lx) = elem(ie)%state%v(i,j,2,k,np1)
             if (n==3) xstate(lx) = elem(ie)%state%p(i,j,k,np1)
            end do  !np
          end do  !np
        end do  !nlev
       end do !ie
       end do !nvar

    call initialize(state_object, lenx, elem, pmean, edge3, &
        hybrid, deriv, dt, tl, nets, nete)

    call init_precon(pre_object, lenx, elem, blkjac, edge1, edge2, edge3, &
        red, deriv, cg, lambdasq, dt, pmean, tl, nets, nete)

    fptr => state_object
    c_ptr_to_object =  c_loc(fptr)
    pptr => pre_object
    c_ptr_to_pre =  c_loc(pptr)

! ForTrilinos interface to use nox and loca, and returns xstate(n+1)
   call noxsolve(size(xstate), xstate, c_ptr_to_object, c_ptr_to_pre)

       lx = 0
       do n=1,nvar
        do ie=nets,nete
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) elem(ie)%state%v(i,j,1,k,np1)=xstate(lx)
             if (n==2) elem(ie)%state%v(i,j,2,k,np1)=xstate(lx)
             if (n==3) elem(ie)%state%p(i,j,k,np1)=xstate(lx)
            end do  !np
          end do  !np
        end do  !nlev
       end do !ie
       end do !nvar

       ! ==========================================
       ! If SW test case 1, velocity is prescribed.
       ! reset v back to initial values
       ! ==========================================

! TODO update with vortex and swirl possibly using set_prescribed_velocity
       do ie=nets,nete
       if (topology == "cube" .and. test_case=="swtc1") then
          do k=1,nlev
             elem(ie)%state%v(:,:,:,k,np1)=tc1_velocity(elem(ie)%spherep,elem(ie)%Dinv)
             elem(ie)%state%v(:,:,:,k,n0 )=elem(ie)%state%v(:,:,:,k,np1)
             elem(ie)%state%v(:,:,:,k,nm1)=elem(ie)%state%v(:,:,:,k,np1)
          end do
       end if
       end do !ie

    call t_stopf('implicit header')

    !$OMP BARRIER

  end subroutine advance_imp_nonstag

#ifdef SPHEREW
  subroutine residual(xstate, fx, nelemd, c_ptr_to_object) bind(C,name='calc_f')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf

    implicit none

    real (c_double) ,intent(in)        :: xstate(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,ps
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,E_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

    call t_startf('FE implicit')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xstate(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xstate(lx)
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xstate(lx)
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne
      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)
       end do
      end do
     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

          v1       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! contra
          v2       = fptr%base(ie)%state%v(i,j,2,k,n0)   ! contra 
          v1p      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! contra
          v2p      = fptr%base(ie)%state%v(i,j,2,k,np1)   ! contra 
          v1m      = fptr%base(ie)%state%v(i,j,1,k,nm1)   ! contra
          v2m      = fptr%base(ie)%state%v(i,j,2,k,nm1)   ! contra 

          ulatlon(i,j,1)=fptr%base(ie)%D(1,1,i,j)*v1 + fptr%base(ie)%D(1,2,i,j)*v2 ! contra->latlon
          ulatlon(i,j,2)=fptr%base(ie)%D(2,1,i,j)*v1 + fptr%base(ie)%D(2,2,i,j)*v2 ! contra->latlon
          up(i,j,1)=fptr%base(ie)%D(1,1,i,j)*v1p + fptr%base(ie)%D(1,2,i,j)*v2p    ! contra->latlon
          up(i,j,2)=fptr%base(ie)%D(2,1,i,j)*v1p + fptr%base(ie)%D(2,2,i,j)*v2p    ! contra->latlon
          um(i,j,1)=fptr%base(ie)%D(1,1,i,j)*v1m + fptr%base(ie)%D(1,2,i,j)*v2m    ! contra->latlon
          um(i,j,2)=fptr%base(ie)%D(2,1,i,j)*v1m + fptr%base(ie)%D(2,2,i,j)*v2m    ! contra->latlon

          E_n(i,j) = 0.5D0*(ulatlon(i,j,1)**2 + ulatlon(i,j,2)**2)  + &
                         fptr%base(ie)%state%p(i,j,k,n0) + ps(i,j)
          E(i,j)   = 0.5D0*(up(i,j,1)**2 + up(i,j,2)**2)  + &
                         fptr%base(ie)%state%p(i,j,k,np1) + ps(i,j)

          pv_n(i,j,1) = ulatlon(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          pv_n(i,j,2) = ulatlon(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          pv(i,j,1) = up(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,np1) )
          pv(i,j,2) = up(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,np1) )

       end do
      end do
! residual time level n
     grade_n = gradient_sphere(E_n,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     div_n = divergence_sphere(pv_n,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 

     grade = gradient_sphere(E,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     zeta = vorticity_sphere(up,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     div = divergence_sphere(pv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
          ! ==============================================
          ! Compute residual terms
          ! ==============================================
          do j=1,np
             do i=1,np

            vtens_n(i,j,1,k,ie)=spheremp(i,j)* &
          (ulatlon(i,j,2)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,1))
            vtens_n(i,j,2,k,ie)=spheremp(i,j)* &
          (-ulatlon(i,j,1)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,2))

            vtens(i,j,1,k,ie)=spheremp(i,j)* &
          (up(i,j,2)*(fcor(i,j) + zeta(i,j)) - grade(i,j,1))
            vtens(i,j,2,k,ie)=spheremp(i,j)* &
          (-up(i,j,1)*(fcor(i,j) + zeta(i,j)) - grade(i,j,2))

            ptens_n(i,j,k,ie) =  -spheremp(i,j)*div_n(i,j)
            ptens(i,j,k,ie) =  -spheremp(i,j)*div(i,j)

! calculate nonlinear residual using Crank-Nicolson

       ptens(i,j,k,ie) = spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,np1)- &
          fptr%base(ie)%state%p(i,j,k,n0))*dti - 0.5*ptens(i,j,k,ie) &
          - 0.5*ptens_n(i,j,k,ie)

      vtens(i,j,1,k,ie) = spheremp(i,j)* &
        (up(i,j,1)-ulatlon(i,j,1))*dti - 0.5*vtens(i,j,1,k,ie) &
          - 0.5*vtens_n(i,j,1,k,ie)

       vtens(i,j,2,k,ie) = spheremp(i,j)* &
         (up(i,j,2)-ulatlon(i,j,2))*dti - 0.5*vtens(i,j,2,k,ie) &
          - 0.5*vtens_n(i,j,2,k,ie)

            end do
          end do
       end do

      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
   end do

   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge3)
   !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, fptr%base(ie)%desc)

       kptr=nlev
       call edgeVunpack(fptr%edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

     end do !ie
       ! ===========================================================
       ! Compute velocity and pressure tendencies for all levels
       ! ===========================================================
      lx = 0
    do n=1,nvar
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
              lx = lx+1
! TODO have if statement include swirl and vortex, possibly using set_prescribed_velocity
           if (topology == "cube" .and. test_case=="swtc1") then
              if (n==1) fx(lx) = 0.0
              if (n==2) fx(lx) = 0.0
           else 
                vtens1=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens2=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

              if (n==1) fx(lx) = fptr%base(ie)%Dinv(1,1,i,j)*vtens1 &
                     + fptr%base(ie)%Dinv(1,2,i,j)*vtens2
              if (n==2) fx(lx) = fptr%base(ie)%Dinv(2,1,i,j)*vtens1 &
                     + fptr%base(ie)%Dinv(2,2,i,j)*vtens2
           end if
              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
            end do
          end do
       end do
     end do !ie
     end do

    call t_stopf('FE implicit')

  end subroutine residual

  subroutine precon_gmres(vv, z, nelemd, xstate, c_ptr_to_object, c_ptr_to_pre) &
!   bind(C,name='precon_gmres')
   bind(C,name='precon')

    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use physical_constants, only : rearth, g
    use element_mod, only : element_t
!    use parallel_mod, only : parallel_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
!    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use shallow_water_mod, only : tc1_velocity
    use bndry_mod, only : bndry_exchangev
    use cg_mod, only : cg_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use derived_type_mod ,only : derived_type, precon_type
!    use global_norms_mod
    use perf_mod, only : t_startf, t_stopf

    implicit none

    ! =================
    ! Local
    ! =================

    ! Thread private working set ...

!    type (parallel_t)          :: par

    integer    :: i,j,k,n,ie
    integer    :: nm1,n0,np1
    integer     :: lx, nn
    integer     :: ns, ne
!    integer    :: comm

    real (c_double) ,intent(out)       :: vv(nelemd) ! combo of dv and dp 
    real (c_double) ,intent(in)        :: z(nelemd)  ! rhs from lin solver
    real (c_double) ,intent(in)        :: xstate(nelemd)
    integer(c_int)  ,intent(in) ,value :: nelemd
    type(derived_type) ,pointer         :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    type(precon_type) ,pointer         :: pptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_pre

    real (c_double)       :: zt(nelemd) ! tmp vv to prevent overwrite
    real (c_double)       :: xt(nelemd) ! current xstate to prevent overwrite

  interface 
    subroutine init_prec(vectorSize,vector,rhs,comm,v_container,p_container) &
        bind(C,name='init_prec')
    use ,intrinsic :: iso_c_binding
      integer(c_int)                :: vectorSize,comm
      real(c_double)  ,dimension(*) :: vector
      real(c_double)  ,dimension(*) :: rhs
      type(c_ptr)                   :: v_container
      type(c_ptr)                   :: p_container
    end subroutine init_prec

    subroutine finish_prec() bind(C,name='finish_prec')
    use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
    end subroutine finish_prec
  
   subroutine precon_solve(vectorSize,vector,rhs,v_container,p_container) &
     bind(C,name='precon_solve')
    use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
      integer(c_int)                :: vectorSize
      real(c_double)  ,dimension(*) :: vector
      real(c_double)  ,dimension(*) :: rhs
      type(c_ptr)                   :: v_container 
      type(c_ptr)                   :: p_container 
    end subroutine precon_solve
  end interface

    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr
    call c_f_pointer(c_ptr_to_pre, pptr) ! convert C ptr to F ptr

    call t_startf('precon_gmres')

    zt = z
!    xt = xstate

!    call init_prec(size(xstate), xt, zt, 1, c_ptr_to_object, c_ptr_to_pre)

! ForTrilinos interface to use applyJacobianInverse with F_l(x)=residual_lin(x)
!    call precon_solve(size(xstate), xt, zt, c_ptr_to_object, c_ptr_to_pre)

    vv = zt

!    call finish_prec()

    call t_stopf('precon_gmres')

    !$OMP BARRIER

  end subroutine precon_gmres

  subroutine residual_lin(xs, fx, nelemd, c_ptr_to_object) bind(C,name='calc_f_lin')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf

    implicit none

    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,ps
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,E_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

    call t_startf('residual lin')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne
      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)
       end do
      end do
     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

          v1       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! contra
          v2       = fptr%base(ie)%state%v(i,j,2,k,n0)   ! contra 
          v1p      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! contra
          v2p      = fptr%base(ie)%state%v(i,j,2,k,np1)   ! contra 
          v1m      = fptr%base(ie)%state%v(i,j,1,k,nm1)   ! contra
          v2m      = fptr%base(ie)%state%v(i,j,2,k,nm1)   ! contra 

          ulatlon(i,j,1)=fptr%base(ie)%D(1,1,i,j)*v1 + fptr%base(ie)%D(1,2,i,j)*v2 ! contra->latlon
          ulatlon(i,j,2)=fptr%base(ie)%D(2,1,i,j)*v1 + fptr%base(ie)%D(2,2,i,j)*v2 ! contra->latlon
          up(i,j,1)=fptr%base(ie)%D(1,1,i,j)*v1p + fptr%base(ie)%D(1,2,i,j)*v2p    ! contra->latlon
          up(i,j,2)=fptr%base(ie)%D(2,1,i,j)*v1p + fptr%base(ie)%D(2,2,i,j)*v2p    ! contra->latlon
          um(i,j,1)=fptr%base(ie)%D(1,1,i,j)*v1m + fptr%base(ie)%D(1,2,i,j)*v2m    ! contra->latlon
          um(i,j,2)=fptr%base(ie)%D(2,1,i,j)*v1m + fptr%base(ie)%D(2,2,i,j)*v2m    ! contra->latlon

          E_n(i,j) = 0.5D0*(ulatlon(i,j,1)**2 + ulatlon(i,j,2)**2)  + &
                         fptr%base(ie)%state%p(i,j,k,n0) + ps(i,j)
          E(i,j)   = 0.5D0*(up(i,j,1)**2 + up(i,j,2)**2)  + &
                         fptr%base(ie)%state%p(i,j,k,np1) + ps(i,j)

          pv_n(i,j,1) = ulatlon(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          pv_n(i,j,2) = ulatlon(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          pv(i,j,1) = up(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,np1) )
          pv(i,j,2) = up(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,np1) )

       end do
      end do
! residual time level n
     grade_n = gradient_sphere(E_n,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     div_n = divergence_sphere(pv_n,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 

     grade = gradient_sphere(E,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     zeta = vorticity_sphere(up,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     div = divergence_sphere(pv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
          ! ==============================================
          ! Compute residual terms
          ! ==============================================
          do j=1,np
             do i=1,np

            vtens_n(i,j,1,k,ie)=spheremp(i,j)* &
          (ulatlon(i,j,2)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,1))
            vtens_n(i,j,2,k,ie)=spheremp(i,j)* &
          (-ulatlon(i,j,1)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,2))

            vtens(i,j,1,k,ie)=spheremp(i,j)* &
          (up(i,j,2)*(fcor(i,j) + zeta(i,j)) - grade(i,j,1))
            vtens(i,j,2,k,ie)=spheremp(i,j)* &
          (-up(i,j,1)*(fcor(i,j) + zeta(i,j)) - grade(i,j,2))

            ptens_n(i,j,k,ie) =  -spheremp(i,j)*div_n(i,j)
            ptens(i,j,k,ie) =  -spheremp(i,j)*div(i,j)

! calculate nonlinear residual 

       ptens(i,j,k,ie) = spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,np1)- &
          fptr%base(ie)%state%p(i,j,k,n0))*dti - 0.5*ptens(i,j,k,ie) &
          - 0.5*ptens_n(i,j,k,ie)

      vtens(i,j,1,k,ie) = spheremp(i,j)* &
        (up(i,j,1)-ulatlon(i,j,1))*dti - 0.5*vtens(i,j,1,k,ie) &
          - 0.5*vtens_n(i,j,1,k,ie)

       vtens(i,j,2,k,ie) = spheremp(i,j)* &
         (up(i,j,2)-ulatlon(i,j,2))*dti - 0.5*vtens(i,j,2,k,ie) &
          - 0.5*vtens_n(i,j,2,k,ie)

            end do
          end do
       end do

      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
   end do

   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge3)
   !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, fptr%base(ie)%desc)

       kptr=nlev
       call edgeVunpack(fptr%edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

     end do !ie
       ! ===========================================================
       ! Compute velocity and pressure tendencies for all levels
       ! ===========================================================
      lx = 0
    do n=1,nvar
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
              lx = lx+1
           if (topology == "cube" .and. test_case=="swtc1") then
              if (n==1) fx(lx) = 0.0
              if (n==2) fx(lx) = 0.0
           else 
                vtens1=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens2=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

              if (n==1) fx(lx) = fptr%base(ie)%Dinv(1,1,i,j)*vtens1 &
                     + fptr%base(ie)%Dinv(1,2,i,j)*vtens2
              if (n==2) fx(lx) = fptr%base(ie)%Dinv(2,1,i,j)*vtens1 &
                     + fptr%base(ie)%Dinv(2,2,i,j)*vtens2
           end if
              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
! for testing, ID preconditioner
!              if (n==1) fx(lx) = xs(lx)
!              if (n==2) fx(lx) = xs(lx)
!              if (n==3) fx(lx) = xs(lx)
            end do
          end do
       end do
     end do !ie
     end do

    call t_stopf('residual lin')

  end subroutine residual_lin

#endif  !SPHEREW

! precon_si() is a preconditioner for the fully implicit solver based on 
! the semi-implicit solver advance_si_nonstag

  subroutine precon_si(vv, z, nelemd, xstate, c_ptr_to_object, c_ptr_to_pre) &
!   bind(C,name='precon')
   bind(C,name='precon_si')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use physical_constants, only : rearth, g
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use derivative_mod, only : derivative_t,  gradient_wk, divergence, &
         vorticity_sphere, gradient_sphere, divergence_sphere
    use time_mod, only : timelevel_t
    use control_mod, only : topology, test_case
    use shallow_water_mod, only : tc1_velocity
    use cg_mod, only : cg_t
    use solver_mod, only : blkjac_t
    use precon_mod, only : pcg_presolver_nonstag
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod ,only : derived_type, precon_type
    use perf_mod, only : t_startf, t_stopf

    implicit none

    real (c_double) ,intent(out)       :: vv(nelemd) ! combo of dv and dp 
    real (c_double) ,intent(in)        :: z(nelemd)  ! rhs from lin solver
    real (c_double) ,intent(in)        :: xstate(nelemd)
    integer(c_int)  ,intent(in) ,value :: nelemd
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    type(precon_type) ,pointer         :: pptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_pre

    ! =================
    ! Local Variables
    ! =================

    real (kind=real_kind), dimension(np,np)     :: mp
    real (kind=real_kind), dimension(np,np)     :: metdetp
    real (kind=real_kind), dimension(np,np)     :: fcor
    real (kind=real_kind), dimension(np,np)     :: rmv
    real (kind=real_kind), dimension(np,np)     :: mv
    real (kind=real_kind), dimension(np,np)     :: rspheremp
    real (kind=real_kind), dimension(np,np)     :: spheremp
    real (kind=real_kind), dimension(np,np)     :: metdet
    real (kind=real_kind), dimension(2,2,np,np) :: met
    real (kind=real_kind), dimension(2,2,np,np) :: metinv

    ! Thread private working set ...

    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv, zv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: Ru
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: grad_dp
    real (kind=real_kind), dimension(np,np,nlev,nelem)    :: vgradp 
    real (kind=real_kind), dimension(np,np,nlev,nelem)    :: Rs   
    real (kind=real_kind), dimension(np,np,nlev,nelem)    :: dp, zp 
    real (kind=real_kind), dimension(np,np,2) :: gradp,grappm1 ! weak p grad, time level n,n-1
    real (kind=real_kind), dimension(np,np,2) :: grade   ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2) :: gv,gvm1 ! metdet*v(n,n-1), v contravar
    real (kind=real_kind), dimension(np,np,2) :: ulatlon     !(cov) latlon velocity
    real (kind=real_kind), dimension(np,np)   :: E       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)   :: zeta    ! relative vorticity
    real (kind=real_kind), dimension(np,np)   :: div,divm1 ! tstep n,n-1 velocity div

    real (kind=real_kind) ::  v1,v2
    real (kind=real_kind) ::  gradp1,gradp2
    real (kind=real_kind) ::  grad_dp1,grad_dp2
    real (kind=real_kind) ::  Ru1,Ru2

    real (kind=real_kind) ::  lenscale
    real (kind=real_kind) ::  pmean
    real (kind=real_kind) ::  time_adv

    real*8  :: et,st
    integer i,j,k,ie, lenx
    integer kptr
    integer point
    integer iptr, ieptr
    integer nm1,n0,np1
    integer ns, ne
    integer nn, lx

    real (kind=real_kind) :: p0sum,v0sum

    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr
    call c_f_pointer(c_ptr_to_pre,pptr) ! convert C ptr to F ptr

    n0         = pptr%tl%n0
    np1        = pptr%tl%np1
    nm1        = pptr%tl%nm1
    lenx       = pptr%n 
    ns         = pptr%nets 
    ne         = pptr%nete 
    pmean      = pptr%pmean

    lenscale=rearth
    call t_startf('precon_si')

    lx = 0
     do nn=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (nn==1) zv(i,j,1,k,ie) = z(lx)
             if (nn==2) zv(i,j,2,k,ie) = z(lx)
             if (nn==3) zp(i,j,k,ie)   = z(lx)
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    !$OMP BARRIER

    do ie=ns,ne

        rspheremp(:,:)         = pptr%base(ie)%rspheremp(:,:)

       do k=1,nlev

          ! =========================================================
          ! Scale velocity tendency and v.grad(p) term by inverse mass
          ! matrix, scale velocity by metric g factor.
          !     11 np*np Flops
          ! =========================================================

          do j=1,np
             do i=1,np

             gvm1(i,j,1) = pptr%base(ie)%D(1,1,i,j)*zv(i,j,1,k,ie) + &
  pptr%base(ie)%D(1,2,i,j)*zv(i,j,2,k,ie)
             gvm1(i,j,2) = pptr%base(ie)%D(2,1,i,j)*zv(i,j,1,k,ie) + &
  pptr%base(ie)%D(2,2,i,j)*zv(i,j,2,k,ie)

             end do
          end do

          ! ==========================================================
          ! Compute divergence 
          ! ==========================================================

           divm1  = divergence_sphere(gvm1,pptr%deriv,pptr%base(ie))

          do j=1,np
             do i=1,np
                Rs(i,j,k,ie) =        &
                 +  pptr%base(ie)%spheremp(i,j)*(zp(i,j,k,ie)  &
                -  pptr%dt*pmean*divm1(i,j) )
             end do
          end do

       end do
       !DBG print *,'advance_si: point #14'
       kptr=0
       call edgeVpack(pptr%edge1, Rs(1,1,1,ie),nlev,kptr,pptr%base(ie)%desc)

    end do

    call bndry_exchangeV(pptr%cg%hybrid,pptr%edge1)
    !$OMP BARRIER

    do ie=ns,ne
       kptr=0
       call edgeVunpack(pptr%edge1, Rs(1,1,1,ie), nlev, kptr, pptr%base(ie)%desc)
       do k=1,nlev
          do j=1,np
             do i=1,np
                Rs(i,j,k,ie) = pptr%base(ie)%rspheremp(i,j)*Rs(i,j,k,ie)
             enddo
          enddo
       enddo
    enddo

    ! ======================================================
    ! Invoke the solver!
    ! ======================================================

    !DBG print *,'advance_si: before call to pcg_presolver'
    point = 3
#ifdef _HTRACE
    !JMD    call EVENT_POINT(point)
#endif
    !$OMP BARRIER

    dp(:,:,:,ns:ne) = pcg_presolver_nonstag(pptr, &
         Rs(:,:,:,ns:ne) )     ! rhs of Helmholtz problem

    point = 4 
#ifdef _HTRACE
    !JMD   call EVENT_POINT(point)
#endif

    do ie=ns,ne

       metinv(:,:,:,:)  = pptr%base(ie)%metinv(:,:,:,:)
       spheremp(:,:)  = pptr%base(ie)%spheremp(:,:)
       rspheremp(:,:)  = pptr%base(ie)%rspheremp(:,:)

       do k=1,nlev

          ! compute grad dp needed to back substitute for du

          grad_dp(:,:,:,k,ie)=gradient_sphere(dp(:,:,k,ie),pptr%deriv,pptr%base(ie)%Dinv)

          ! ==================================================
          ! Rotate grad_dp to form contravariant object:
          !        6 np*np Flops
          ! ==================================================

          do j=1,np
             do i=1,np
                grad_dp(i,j,1,k,ie) = spheremp(i,j)*grad_dp(i,j,1,k,ie) 
                grad_dp(i,j,2,k,ie) = spheremp(i,j)*grad_dp(i,j,2,k,ie) 
             end do
          end do
       enddo

       kptr=0
       call edgeVpack(pptr%edge2, grad_dp(1,1,1,1,ie),2*nlev,kptr,pptr%base(ie)%desc)
    end do

    !$OMP BARRIER
    call bndry_exchangeV(pptr%cg%hybrid,pptr%edge2)
    !$OMP BARRIER
    do ie=ns,ne

       kptr=0      
       call edgeVunpack(pptr%edge2, grad_dp(1,1,1,1,ie), 2*nlev, kptr, pptr%base(ie)%desc)

       do k=1,nlev

          ! ==============================
          ! Update velocity
          !    16 np*np Flops
          ! ==============================

          do j=1,np
             do i=1,np
               grad_dp1 = pptr%base(ie)%rspheremp(i,j)*grad_dp(i,j,1,k,ie)
               grad_dp2 = pptr%base(ie)%rspheremp(i,j)*grad_dp(i,j,2,k,ie)
           grad_dp(i,j,1,k,ie)   = pptr%base(ie)%Dinv(1,1,i,j)*grad_dp1 + &
                      pptr%base(ie)%Dinv(1,2,i,j)*grad_dp2
           grad_dp(i,j,2,k,ie)   = pptr%base(ie)%Dinv(2,1,i,j)*grad_dp1 + &
                      pptr%base(ie)%Dinv(2,2,i,j)*grad_dp2

               dv(i,j,1,k,ie) = zv(i,j,1,k,ie) - pptr%dt*grad_dp(i,j,1,k,ie)
               dv(i,j,2,k,ie) = zv(i,j,2,k,ie) - pptr%dt*grad_dp(i,j,2,k,ie)  

             end do
          end do
       end do

    end do ! ie loop

    lx = 0
     do nn=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
! for testing, ID preconditioner 
!              vv(lx) = z(lx)
           if (topology == "cube" .and. test_case=="swtc1") then
             if (nn.ne.3) vv(lx) = zv(i,j,nn,k,ie)
             if (nn==3)   vv(lx) = dp(i,j,k,ie)
           else 
             if (nn.ne.3) vv(lx) = dv(i,j,nn,k,ie)
             if (nn==3)   vv(lx) = dp(i,j,k,ie)
           end if
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar
    !$OMP BARRIER
    call t_stopf('precon_si')

  end subroutine precon_si

end module implicit_mod

