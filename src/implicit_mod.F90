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

    type (element_t)  :: pc_elem(size(elem))
    type (element_t)  :: jac_elem(size(elem))

!    type (element_t)  :: pc_elem(nets-nete+1)
!    type (element_t)  :: jac_elem(nets-nete+1)

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
!    type(c_funptr)                     :: res_calc
    !type(precon_type) ,target          :: pre_object
    !type(precon_type) ,pointer         :: pptr=>NULL()

    type(derived_type) ,target          :: pre_object
    type(derived_type) ,pointer         :: pptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_pre

    type(derived_type) ,target          :: ajac_object
    type(derived_type) ,pointer         :: jptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_jac


!    type(c_funptr)                     :: precon

  interface 
   subroutine noxsolve(vectorSize,vector,v_container,p_container,j_container) &
     bind(C,name='noxsolve')
    use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
      integer(c_int)                :: vectorSize
      real(c_double)  ,dimension(*) :: vector
      type(c_ptr)                   :: v_container 
      type(c_ptr)                   :: p_container  !precon ptr
      type(c_ptr)                   :: j_container  !precon ptr
    end subroutine noxsolve
  end interface

    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    comm = par%comm
! n+1 guess is result at n - need this to match up with time_mod bootstrap style
     call t_startf('implicit header')
     call t_startf('implicit_pre_noxsolve')

!Moved conditionals out of do loops
     lx = 0
     do ie=nets,nete
       do k=1,nlev
         do j=1,np
           do i=1,np
             lx = lx+1
             xstate(lx) = elem(ie)%state%v(i,j,1,k,n0)
           end do  !np
         end do  !np
       end do  !nlev
     end do !ie

     do ie=nets,nete
       do k=1,nlev
         do j=1,np
           do i=1,np
             lx = lx+1
             xstate(lx) = elem(ie)%state%v(i,j,2,k,n0)
           end do  !np
         end do  !np
       end do  !nlev
     end do !ie

     do ie=nets,nete
       do k=1,nlev
         do j=1,np
           do i=1,np
             lx = lx+1
             xstate(lx) = elem(ie)%state%p(i,j,k,n0)
           end do  !np
         end do  !np
       end do  !nlev
     end do !ie

     call t_stopf('implicit_pre_noxsolve')

     call t_startf('implicit_init')

     pc_elem=elem
     jac_elem=elem

    call initialize(state_object, lenx, elem, pmean, edge1,edge2,edge3, &
        hybrid, deriv, dt, tl, nets, nete)

    call initialize(pre_object, lenx, pc_elem, pmean, edge1,edge2,edge3, &
        hybrid, deriv, dt, tl, nets, nete)

    call initialize(ajac_object, lenx, jac_elem, pmean, edge1,edge2,edge3, &
        hybrid, deriv, dt, tl, nets, nete)

    
   ! call init_precon(pre_object, lenx, elem, blkjac, edge1, edge2, edge3, &
   !     red, deriv, cg, lambdasq, dt, pmean, tl, nets, nete)

    fptr => state_object
    c_ptr_to_object =  c_loc(fptr)

    pptr => pre_object
    c_ptr_to_pre =  c_loc(pptr)

    jptr => ajac_object
    c_ptr_to_jac =  c_loc(jptr)
    call t_stopf('implicit_init')

! ForTrilinos interface to use nox and loca, and returns xstate(n+1)

    call t_startf('noxsolve')
    call noxsolve(size(xstate), xstate, c_ptr_to_object, c_ptr_to_pre,c_ptr_to_jac)
    call t_stopf('noxsolve')
    call t_startf('implicit_post_noxsolve')

!Moved conditionals out of do loops
      lx = 0
      do ie=nets,nete
        do k=1,nlev
          do j=1,np
            do i=1,np
              lx = lx+1
              elem(ie)%state%v(i,j,1,k,np1)=xstate(lx)
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie

      do ie=nets,nete
        do k=1,nlev
          do j=1,np
            do i=1,np
              lx = lx+1
              elem(ie)%state%v(i,j,2,k,np1)=xstate(lx)
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie

      do ie=nets,nete
        do k=1,nlev
          do j=1,np
            do i=1,np
              lx = lx+1
              elem(ie)%state%p(i,j,k,np1)=xstate(lx)
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie

     call t_stopf('implicit_post_noxsolve')

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
    use control_mod, only :  topology, test_case, tstep_type
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
    call t_startf('FE_implicit_init')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

!Moved conditionals out of do loops
      lx = 0
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
              lx = lx+1
              fptr%base(ie)%state%v(i,j,1,k,np1) = xstate(lx)
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie

      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
              lx = lx+1
              fptr%base(ie)%state%v(i,j,2,k,np1) = xstate(lx)
            end do  !np
          end do  !np
        end do  !nlev
       end do !ie

       do ie=ns,ne
         do k=1,nlev
           do j=1,np
             do i=1,np
               lx = lx+1
               fptr%base(ie)%state%p(i,j,k,np1)  = xstate(lx)
             end do  !np
           end do  !np
         end do  !nlev
       end do !ie

       call t_stopf('FE_implicit_init')

       call t_startf('FE_implicit_KE_resid_calc')

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

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   !
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !

! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   !
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   ! 

! set un-1
          um(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,nm1)   ! 
          um(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,nm1)   !  


          E(i,j)   = 0.5D0*(up(i,j,1)**2 + up(i,j,2)**2)  + &
                         fptr%base(ie)%state%p(i,j,k,np1) + ps(i,j)
          if (tstep_type==12) then  !compute extra terms needed for CN
             E_n(i,j) = 0.5D0*(ulatlon(i,j,1)**2 + ulatlon(i,j,2)**2)  + &
                         fptr%base(ie)%state%p(i,j,k,n0) + ps(i,j)
             pv_n(i,j,1) = ulatlon(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
             pv_n(i,j,2) = ulatlon(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          endif

          pv(i,j,1) = up(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,np1) )
          pv(i,j,2) = up(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,np1) )

       end do
      end do

!call flush(6)
!stop

! residual time level n
     if (tstep_type==12) then  !compute extra terms needed for CN
        grade_n = gradient_sphere(E_n,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
        zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
        div_n = divergence_sphere(pv_n,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     endif 
     grade = gradient_sphere(E,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     zeta = vorticity_sphere(up,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     div = divergence_sphere(pv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
          ! ==============================================
          ! Compute residual terms
          ! ==============================================
          do j=1,np
             do i=1,np

           if (tstep_type==12) then  !compute extra terms needed for CN
               vtens_n(i,j,1,k,ie)=spheremp(i,j)* &
                  (ulatlon(i,j,2)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,1))
               vtens_n(i,j,2,k,ie)=spheremp(i,j)* &
                  (-ulatlon(i,j,1)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,2))
               ptens_n(i,j,k,ie) =  -spheremp(i,j)*div_n(i,j)
           endif

            vtens(i,j,1,k,ie)=spheremp(i,j)* &
               (up(i,j,2)*(fcor(i,j) + zeta(i,j)) - grade(i,j,1))
            vtens(i,j,2,k,ie)=spheremp(i,j)* &
               (-up(i,j,1)*(fcor(i,j) + zeta(i,j)) - grade(i,j,2))

            ptens(i,j,k,ie) =  -spheremp(i,j)*div(i,j)

! calculate nonlinear residual 
    if (tstep_type==12) then !Crank Nicolson 2nd order

       ptens(i,j,k,ie) = spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,np1)- &
          fptr%base(ie)%state%p(i,j,k,n0))*dti - 0.5*ptens(i,j,k,ie) &
          - 0.5*ptens_n(i,j,k,ie)

      vtens(i,j,1,k,ie) = spheremp(i,j)* &
        (up(i,j,1)-ulatlon(i,j,1))*dti - 0.5*vtens(i,j,1,k,ie) &
          - 0.5*vtens_n(i,j,1,k,ie)

       vtens(i,j,2,k,ie) = spheremp(i,j)* &
         (up(i,j,2)-ulatlon(i,j,2))*dti - 0.5*vtens(i,j,2,k,ie) &
          - 0.5*vtens_n(i,j,2,k,ie)

    else if (tstep_type==13) then !BDF2 2nd order

      if (nstep==0) then ! BE bootstrap
       ptens(i,j,k,ie)   = spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,np1)- &
          fptr%base(ie)%state%p(i,j,k,n0))*dti - ptens(i,j,k,ie)

       vtens(i,j,1,k,ie) = spheremp(i,j)* &
         (up(i,j,1)-ulatlon(i,j,1))*dti - vtens(i,j,1,k,ie)

       vtens(i,j,2,k,ie) = spheremp(i,j)* &
         (up(i,j,2)-ulatlon(i,j,2))*dti - vtens(i,j,2,k,ie)

!      if (nstep==0) then ! CN bootstrap
!       ptens(i,j,k,ie) = spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,np1)- &
!          fptr%base(ie)%state%p(i,j,k,n0))*dti - 0.5*ptens(i,j,k,ie) &
!          - 0.5*ptens_n(i,j,k,ie)
!
!       vtens(i,j,1,k,ie) = spheremp(i,j)* &
!         (up(i,j,1)-ulatlon(i,j,1))*dti - 0.5*vtens(i,j,1,k,ie) &
!          - 0.5*vtens_n(i,j,1,k,ie)
!
!       vtens(i,j,2,k,ie) = spheremp(i,j)* &
!         (up(i,j,2)-ulatlon(i,j,2))*dti - 0.5*vtens(i,j,2,k,ie) &
!          - 0.5*vtens_n(i,j,2,k,ie)

      else 

       ptens(i,j,k,ie) = &
       (1+gam)*spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,np1) - &
          fptr%base(ie)%state%p(i,j,k,n0))*dti - &
          gam*spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,n0) - &
          fptr%base(ie)%state%p(i,j,k,nm1))*dti - &
            ptens(i,j,k,ie) 

       vtens(i,j,1,k,ie) = &
       (1+gam)*spheremp(i,j)*(up(i,j,1)-ulatlon(i,j,1))*dti - &
          gam*spheremp(i,j)*(ulatlon(i,j,1)-um(i,j,1))*dti - &
          vtens(i,j,1,k,ie)

       vtens(i,j,2,k,ie) = &
       (1+gam)*spheremp(i,j)*(up(i,j,2)-ulatlon(i,j,2))*dti - &
          gam*spheremp(i,j)*(ulatlon(i,j,2)-um(i,j,2))*dti - &
          vtens(i,j,2,k,ie)
       end if 

    else ! Backward Euler 1st order method
       ptens(i,j,k,ie)   = spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,np1)- &
          fptr%base(ie)%state%p(i,j,k,n0))*dti - ptens(i,j,k,ie)

       vtens(i,j,1,k,ie) = spheremp(i,j)* &
         (up(i,j,1)-ulatlon(i,j,1))*dti - vtens(i,j,1,k,ie)

       vtens(i,j,2,k,ie) = spheremp(i,j)* &
         (up(i,j,2)-ulatlon(i,j,2))*dti - vtens(i,j,2,k,ie)

    end if
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

   !pw++
   call t_stopf('FE_implicit_KE_resid_calc')
   !pw--

   !$OMP BARRIER
   !pw++
   call t_startf('FE_implicit_bndry_ex')
   !pw--
   !$OMP BARRIER
   call bndry_exchangeV(fptr%hybrid,fptr%edge3)
   !$OMP BARRIER
   !pw++
   call t_stopf('FE_implicit_bndry_ex')
   !pw--
   !$OMP BARRIER

   !pw++
   call t_startf('FE_implicit_bndry_unpack')
   !pw--
   

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, fptr%base(ie)%desc)

       kptr=nlev
       call edgeVunpack(fptr%edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

     end do !ie

         call t_stopf('FE_implicit_bndry_unpack')

         ! ===========================================================
         ! Compute velocity and pressure tendencies for all levels
         ! ===========================================================

         call t_startf('FE_implicit_vel_pres')
 

!Moved conditionals out of do loops
       lx = 0
 if (topology == "cube" .and. test_case=="swtc1") then
       do ie=ns,ne
         do k=1,nlev
           do j=1,np
             do i=1,np
               lx = lx+1
               fx(lx) = 0.0
             end do
           end do
         end do
       end do !ie

       do ie=ns,ne
         do k=1,nlev
           do j=1,np
             do i=1,np
               lx = lx+1
               fx(lx) = 0.0
             end do
           end do
         end do
       end do !ie

       do ie=ns,ne
         do k=1,nlev
           do j=1,np
             do i=1,np
               lx = lx+1
               fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
             end do
           end do
         end do
       end do !ie

 else

       do ie=ns,ne
         do k=1,nlev
           do j=1,np
             do i=1,np
               lx = lx+1
               fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
             end do
           end do
         end do
       end do !ie

       do ie=ns,ne
         do k=1,nlev
           do j=1,np
             do i=1,np
               lx = lx+1
               fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
             end do
           end do
         end do
       end do !ie

       do ie=ns,ne
         do k=1,nlev
           do j=1,np
             do i=1,np
               lx = lx+1
               fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
             end do
           end do
         end do
       end do !ie
     end if

     call t_stopf('FE_implicit_vel_pres')

     call t_stopf('FE implicit')

  end subroutine residual

  subroutine precon_gmres(vv, z, nelemd, xstate, c_ptr_to_object, c_ptr_to_pre) &
   bind(C,name='precon')

    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use physical_constants, only : g
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

    integer    :: i,j,k,n,ie
    integer    :: nm1,n0,np1
    integer     :: lx, nn
    integer     :: ns, ne

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
  
!   subroutine precon_solve(vectorSize,vector,rhs,v_container,p_container) &
!     bind(C,name='precon_solve')
!    use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
!      integer(c_int)                :: vectorSize
!      real(c_double)  ,dimension(*) :: vector
!      real(c_double)  ,dimension(*) :: rhs
!      type(c_ptr)                   :: v_container 
!      type(c_ptr)                   :: p_container 
!    end subroutine precon_solve

! Aaron! 
!   subroutine precon_solve_wgmres(vectorSize,vector,rhs,v_container,p_container) &
!     bind(C,name='precon_wgmres')
!    use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
!      integer(c_int)                :: vectorSize
!      real(c_double)  ,dimension(*) :: vector
!      real(c_double)  ,dimension(*) :: rhs
!      type(c_ptr)                   :: v_container 
!      type(c_ptr)                   :: p_container 
!    end subroutine precon_solve

  end interface

    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr
    call c_f_pointer(c_ptr_to_pre, pptr) ! convert C ptr to F ptr

    call t_startf('precon_gmres')

    zt = z
    xt = xstate

! ForTrilinos interface to use applyJacobianInverse with F_l(x)=residual_lin(x)
!    call precon_solve(size(xstate), xt, zt, c_ptr_to_object, c_ptr_to_pre)
! Aaron! this will point to a GMRES solve by trilinos
!    call precon_precon_wgmres(size(xstate), xt, zt, c_ptr_to_object, c_ptr_to_pre)

    vv = zt

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
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case, tstep_type
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none

    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

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
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)
       end do
      end do

     do k=1,nlev

!       usum1       = 0.0d0
!       usum2       = 0.0d0
!       psum        = 0.0d0
!       tot_nv      = 1.0d0/float(nv*nv)

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  

! set p=gh=phi-ghs=phi-pmean 
          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)-pmean
! set dp
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)

!          usum1 = usum1 + ulatlon(i,j,1)
!          usum2 = usum2 + ulatlon(i,j,2)
!          psum = psum + fptr%base(ie)%state%p(i,j,k,n0)

       end do
      end do

!          usum1 = usum1*tot_np
!          usum2 = usum2*tot_np
!          psum =  psum*tot_np
!      do j=1,np
!       do i=1,np
!          E_n(i,j) = 0.5D0*(usum1*ulatlon(i,j,1) + usum2*ulatlon(i,j,2))  + &
!                         fptr%base(ie)%state%p(i,j,k,n0) + ps(i,j)
!          E(i,j)   = 0.5D0*(usum1*up(i,j,1) + usum2*up(i,j,2))  + &
!                         fptr%base(ie)%state%p(i,j,k,np1) + ps(i,j)
! Aaron! figure out what E should be. u_old*u_new? or some form?
!          E_n(i,j) = 0.5D0*(ulatlon(i,j,1)**2 + ulatlon(i,j,2)**2)  + &
!                         fptr%base(ie)%state%p(i,j,k,n0) + ps(i,j)
!          E(i,j)   = 0.5D0*(up(i,j,1)**2 + up(i,j,2)**2)  + &
!                         fptr%base(ie)%state%p(i,j,k,np1) + ps(i,j)
! orig
!          pv_n(i,j,1) = ulatlon(i,j,1)*( pmean + fptr%base(ie)%state%p(i,j,k,n0) )
!          pv_n(i,j,2) = ulatlon(i,j,2)*( pmean + fptr%base(ie)%state%p(i,j,k,n0) )
!          pv(i,j,1) = up(i,j,1)*( pmean + fptr%base(ie)%state%p(i,j,k,np1) )
!          pv(i,j,2) = up(i,j,2)*( pmean + fptr%base(ie)%state%p(i,j,k,np1) )

! this is element (3,3) of precon matrix, "G" in the notes
          !pv_n(i,j,1) = usum1*( fptr%base(ie)%state%p(i,j,k,n0) )
          !pv_n(i,j,2) = usum2*( fptr%base(ie)%state%p(i,j,k,n0) )
          !pv(i,j,1) = usum1*( fptr%base(ie)%state%p(i,j,k,np1) )
          !pv(i,j,2) = usum2*( fptr%base(ie)%state%p(i,j,k,np1) )
! set dh
!          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)
!         dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)
!       end do
!      end do




! create operators for full residual calc
     !grade_n = gradient_sphere(E_n,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector
     grade_n1 = gradient_sphere(up(:,:,1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     grade_n2 = gradient_sphere(up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
!     div_n = divergence_sphere(pv_n,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 

     grade = gradient_sphere(E,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
!     zeta = vorticity_sphere(up,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
!     div = divergence_sphere(pv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 

          ! ==============================================
          ! Compute G-D(diagF)^-1 D^T(dh) 
          ! ==============================================

! compute terms in F_diag
!     gradef_n = gradient_sphere_diag(Ef_n,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector
!     gradef = gradient_sphere_diag(Ef,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector
!     zetaf_n = vorticity_sphere_diag(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
!     zetaf = vorticity_sphere_diag(up,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
! adjust these to be like v_s_diag but with only one component 
!     zetaf1 = vorticity_sphere_vx(up,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
!     zetaf2 = halfvorticity_sphere_uy(up,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 

! compute diagF that just has dt term
       do j=1,np
        do i=1,np
! zetaf here should be the one half
       diagf(i,j,1,k,ie) = spheremp(i,j)*dti 
       !+ rspheremp(i,j)*(ulatlon(i,j,2)*zetaf(i,j))
!       diagf(i,j,1,k,ie) = rspheremp(i,j)*dti + rspheremp(i,j)*(ulatlon(i,j,2)*zetaf1(i,j))
! zetaf here should be the other half
       diagf(i,j,2,k,ie) = spheremp(i,j)*dti 
       !+ rspheremp(i,j)*(-ulatlon(i,j,1)*zetaf(i,j))
!       diagf(i,j,2,k,ie) = rspheremp(i,j)*dti + rspheremp(i,j)*(-ulatlon(i,j,1)*zetaf2(i,j))
! compute (diagF)^-1
         diagf(i,j,1,k,ie)   =  1.0d0/(diagf(i,j,1,k,ie))
         diagf(i,j,2,k,ie)   =  1.0d0/(diagf(i,j,2,k,ie))

         !schur(i,j,1)   =  fptr%base(ie)%state%p(i,j,k,n0)*diagf(i,j,1,k,ie)*grade_n3(i,j,1)
         !schur(i,j,2)   =  fptr%base(ie)%state%p(i,j,k,n0)*diagf(i,j,2,k,ie)*grade_n3(i,j,2)

         schur(i,j,1)   =  diagf(i,j,1,k,ie)*grade_n3(i,j,1)
         schur(i,j,2)   =  diagf(i,j,2,k,ie)*grade_n3(i,j,2)

       end do 
      end do

! compute D (dh*diagF^-1*Dtranspose)
     div_h =  dh(i,j)*divergence_sphere(schur,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 


          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

            vtens(i,j,1,k,ie)= spheremp(i,j)*&
            (ulatlon(i,j,1)*grade_n1(i,j,1)+&
             ulatlon(i,j,2)*grade_n1(i,j,2)+&
             up(i,j,2)*(fcor(i,j) - zeta_n(i,j)) + grade(i,j,1))

           vtens(i,j,2,k,ie)= spheremp(i,j)*&
            (ulatlon(i,j,1)*grade_n2(i,j,1)+&
             ulatlon(i,j,2)*grade_n2(i,j,2)+&
             up(i,j,2)*(fcor(i,j) + zeta_n(i,j)) + grade(i,j,1))


           ! vtens(i,j,2,k,ie)=spheremp(i,j)* &
          !(-up(i,j,1)*(fcor(i,j) + zeta(i,j)) - grade(i,j,2))
          !  vtens_n(i,j,1,k,ie)=spheremp(i,j)* &
          !(ulatlon(i,j,2)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,1))
          !  vtens_n(i,j,2,k,ie)=spheremp(i,j)* &
          !(-ulatlon(i,j,1)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,2))

!            vtens(i,j,1,k,ie)=spheremp(i,j)* &
!          (usum2*(fcor(i,j) + zeta(i,j)) - grade(i,j,1))
!            vtens(i,j,2,k,ie)=spheremp(i,j)* &
!          (-usum1*(fcor(i,j) + zeta(i,j)) - grade(i,j,2))
!            vtens_n(i,j,1,k,ie)=spheremp(i,j)* &
!          (ulatlon(i,j,2)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,1))
!            vtens_n(i,j,2,k,ie)=spheremp(i,j)* &
!          (-ulatlon(i,j,1)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,2))

         ptens(i,j,k,ie) = spheremp(i,j)*&
            ( ulatlon(i,j,1)*grade_n3(i,j,1)+&
             ulatlon(i,j,2)*grade_n3(i,j,2)) +spheremp(i,j)*dh_n(i,j)*div_h(i,j)
!  - ( approx schur complement + G )
      ! ptens(i,j,k,ie) =  -spheremp(i,j)*div_h(i,j) + ptens(i,j,k,ie)



! calculate nonlinear residual 
    if (tstep_type==12) then !Crank Nicolson 2nd order

      vtens(i,j,1,k,ie) = spheremp(i,j)* &
        (up(i,j,1)-ulatlon(i,j,1))*dti - 0.5*vtens(i,j,1,k,ie) &
          - 0.5*vtens_n(i,j,1,k,ie)

       vtens(i,j,2,k,ie) = spheremp(i,j)* &
         (up(i,j,2)-ulatlon(i,j,2))*dti - 0.5*vtens(i,j,2,k,ie) &
          - 0.5*vtens_n(i,j,2,k,ie)

       ptens(i,j,k,ie) = spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,np1)- &
          fptr%base(ie)%state%p(i,j,k,n0))*dti - 0.5*ptens(i,j,k,ie) &
          - 0.5*ptens_n(i,j,k,ie)

    else if (tstep_type==13) then !BDF2 2nd order

      if (nstep==0) then ! CN bootstrap

       vtens(i,j,1,k,ie) = spheremp(i,j)* &
         (up(i,j,1)-ulatlon(i,j,1))*dti - 0.5*vtens(i,j,1,k,ie) &
          - 0.5*vtens_n(i,j,1,k,ie)

       vtens(i,j,2,k,ie) = spheremp(i,j)* &
         (up(i,j,2)-ulatlon(i,j,2))*dti - 0.5*vtens(i,j,2,k,ie) &
          - 0.5*vtens_n(i,j,2,k,ie)

       ptens(i,j,k,ie) = spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,np1)- &
         fptr%base(ie)%state%p(i,j,k,n0))*dti - 0.5*ptens(i,j,k,ie) &
          - 0.5*ptens_n(i,j,k,ie)

      else 

       vtens(i,j,1,k,ie) = &
       (1+gam)*spheremp(i,j)*(up(i,j,1)-ulatlon(i,j,1))*dti - &
          gam*spheremp(i,j)*(ulatlon(i,j,1)-um(i,j,1))*dti - &
          vtens(i,j,1,k,ie)

       vtens(i,j,2,k,ie) = &
       (1+gam)*spheremp(i,j)*(up(i,j,2)-ulatlon(i,j,2))*dti - &
          gam*spheremp(i,j)*(ulatlon(i,j,2)-um(i,j,2))*dti - &
          vtens(i,j,2,k,ie)

       ptens(i,j,k,ie) =  &
       (1+gam)*spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,np1) - &
          fptr%base(ie)%state%p(i,j,k,n0))*dti - &
          gam*spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,n0) - &
          fptr%base(ie)%state%p(i,j,k,nm1))*dti - &
            ptens(i,j,k,ie) 

       end if 

    else ! Backward Euler 1st order method

!       vtens(i,j,1,k,ie) = spheremp(i,j)* &
!         (up(i,j,1)-ulatlon(i,j,1))*dti - vtens(i,j,1,k,ie)
       vtens(i,j,1,k,ie) = spheremp(i,j)*up(i,j,1)*dti + vtens(i,j,1,k,ie)

!       vtens(i,j,2,k,ie) = spheremp(i,j)* &
!         (up(i,j,2)-ulatlon(i,j,2))*dti - vtens(i,j,2,k,ie)
       vtens(i,j,2,k,ie) = spheremp(i,j)*up(i,j,2)*dti + vtens(i,j,2,k,ie)

!       ptens(i,j,k,ie)   = spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,np1)- &
!          fptr%base(ie)%state%p(i,j,k,n0))*dti - ptens(i,j,k,ie)
       ptens(i,j,k,ie)   = spheremp(i,j)*dh(i,j)*dti + ptens(i,j,k,ie)

    end if
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

              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

!              if (n==1) fx(lx) = xs(lx)
!              if (n==2) fx(lx) = xs(lx)
           end if
              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
!              if (n==3) fx(lx) = xs(lx)
            end do
          end do
       end do
     end do !ie
     end do

    call t_stopf('residual lin')

  end subroutine residual_lin

! precon_si() is a preconditioner for the fully implicit solver based on 
! the semi-implicit solver advance_si_nonstag

  subroutine precon_si(vv, z, nelemd, xstate, c_ptr_to_object, c_ptr_to_pre) &
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


             gvm1(i,j,1) = pptr%base(ie)%D(1,1,i,j)*zv(i,j,1,k,ie) + pptr%base(ie)%D(1,2,i,j)*zv(i,j,2,k,ie)
             gvm1(i,j,2) = pptr%base(ie)%D(2,1,i,j)*zv(i,j,1,k,ie) + pptr%base(ie)%D(2,2,i,j)*zv(i,j,2,k,ie)

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
! identity precondioner test
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


  subroutine sw_picard_fdiag_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_fdiag')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0

!write(6,*)'dti=',dti,'\n'
!xs=1.0d0
!write(6,*)'xs=[',xs,']\n'

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   ! 


          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)



! set dh_n=p_notes=g(h+hmean)=phi
          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)


          pv1(i,j,1) = up(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          pv1(i,j,2) = 0.d00

          pv2(i,j,1) = 0.d00
          pv2(i,j,2) = up(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )

          dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
          dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )

       end do !np
      end do !np

! create operators for full residual calc
     grade_n1 = gradient_sphere(up(:,:,1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_1
     grade_n2 = gradient_sphere(up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_2
     grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi

     !newton_grade_n1 = gradient_sphere(ulatlon(:,:,1)*up(:,:,1) ,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_1
     !newton_grade_n2 = gradient_sphere(ulatlon(:,:,2)*up(:,:,2),fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_2
     !newton_grade_n3 = gradient_sphere(dh_n,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector        (grad)phi

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
    ! newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
    ! newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w

!for diagonal approximation to f
     newton_vortdelv2 = vorticity_sphere_diag(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere_diag(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w



     divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     !divdpv = divergence_sphere(dpv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

! Jacobian


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*dti - &                 !delv1/dt
             ulatlon(i,j,2)*newton_vortdelv1(i,j)+& !u2*d(delv1)/dy & dx
!!            newton_grade_n1(i,j,1)+ &      !delv1*d(v1)/dx
!            -up(i,j,2)*(zeta_n(i,j)+fcor(i,j))+& !delv2*(w+f)
!            -ulatlon(i,j,2)*newton_vortdelv2(i,j)+& !v2*d(delv2)/dx &dy
!!            newton_grade_n2(i,j,1)+ &      !-delv2*d(v2)/dx
            grade_n3(i,j,1)) ! + d(delp)/dx



           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
!             up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
!             ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
!!             newton_grade_n2(i,j,2)+&  !delu1*d(u2)/dy
             up(i,j,2)*dti +&                !delu2/dt
             ulatlon(i,j,1)*newton_vortdelv2(i,j)+& !u1*d(delu2)/dy & dx
!!            newton_grade_n1(i,j,2)+ &      !delu1*d(u1)/dy
             grade_n3(i,j,2))!d(delp)/dy


         ptens(i,j,k,ie) = spheremp(i,j)*(&
             dh(i,j)*dti+&                         !delp/dt
             divpv1(i,j)+divpv2(i,j))

            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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

              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

           end if
              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
            end do
          end do
       end do
     end do !ie
     end do


  end subroutine sw_picard_fdiag_op





  subroutine sw_picard_simple_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_simple')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv,divgrad
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0

!write(6,*)'dti=',dti,'\n'
!xs=1.0d0
!write(6,*)'xs=[',xs,']\n'

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   !
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   ! 


          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)



! set dh_n=p_notes=g(h+hmean)=phi
          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)


          pv1(i,j,1) = up(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          pv1(i,j,2) = 0.d00

          pv2(i,j,1) = 0.d00
          pv2(i,j,2) = up(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )

          dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
          dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )

       end do !np
      end do !np

! create operators for full residual calc
     grade_n1 = gradient_sphere(up(:,:,1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_1
     grade_n2 = gradient_sphere(up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_2
     grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi

     !newton_grade_n1 = gradient_sphere(ulatlon(:,:,1)*up(:,:,1) ,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_1
     !newton_grade_n2 = gradient_sphere(ulatlon(:,:,2)*up(:,:,2),fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_2
     !newton_grade_n3 = gradient_sphere(dh_n,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector        (grad)phi

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
    ! newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
    ! newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w

!for diagonal approximation to f
     newton_vortdelv2 = vorticity_sphere_diag(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere_diag(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w



     divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divdpv = divergence_sphere(dpv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 

     divgrad= divergence_sphere((1.0d0/dti)*grade_n3,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

! Jacobian


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*dti - &                 !delv1/dt
             ulatlon(i,j,2)*newton_vortdelv1(i,j)& !u2*d(delv1)/dy & dx
!!            newton_grade_n1(i,j,1)+ &      !delv1*d(v1)/dx
!            -up(i,j,2)*(zeta_n(i,j)+fcor(i,j))+& !delv2*(w+f)
!            -ulatlon(i,j,2)*newton_vortdelv2(i,j)+& !v2*d(delv2)/dx &dy
!!            newton_grade_n2(i,j,1)+ &      !-delv2*d(v2)/dx
!            grade_n3(i,j,1)) ! + d(delp)/dx
)



           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
!             up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
!             ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
!!             newton_grade_n2(i,j,2)+&  !delu1*d(u2)/dy
             up(i,j,2)*dti +&                !delu2/dt
             ulatlon(i,j,1)*newton_vortdelv2(i,j)& !u1*d(delu2)/dy & dx
!!            newton_grade_n1(i,j,2)+ &      !delu1*d(u1)/dy
!             grade_n3(i,j,2))!d(delp)/dy
)


         ptens(i,j,k,ie) = spheremp(i,j)*(&
             dh(i,j)*dti+divpv1(i,j)+divpv2(i,j)+&                         !delp/dt
 divdpv(i,j)-divgrad(i,j))

            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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


              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)


           end if
              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
            end do
          end do
       end do
     end do !ie
     end do


  end subroutine sw_picard_simple_op

  subroutine sw_picard_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0

!write(6,*)'dti=',dti,'\n'
!xs=1.0d0
!write(6,*)'xs=[',xs,']\n'

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  


          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)



! set dh_n=p_notes=g(h+hmean)=phi
          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)

          pv1(i,j,1) = up(i,j,1)
          pv1(i,j,2) = 0.d00
          pv2(i,j,1) = 0.d00
          pv2(i,j,2) = up(i,j,2)

          !pv1(i,j,1) = up(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          !pv1(i,j,2) = 0.d00
          !pv2(i,j,1) = 0.d00
          !pv2(i,j,2) = up(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )

          dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
          dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )

       end do !np
      end do !np

! create operators for full residual calc
     grade_n1 = gradient_sphere(up(:,:,1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_1
     grade_n2 = gradient_sphere(up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_2
     grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi

     !newton_grade_n1 = gradient_sphere(ulatlon(:,:,1)*up(:,:,1) ,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_1
     !newton_grade_n2 = gradient_sphere(ulatlon(:,:,2)*up(:,:,2),fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_2
     !newton_grade_n3 = gradient_sphere(dh_n,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector        (grad)phi

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w

     divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divdpv = divergence_sphere(dpv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*dti - &                 !delv1/dt
            ulatlon(i,j,2)*newton_vortdelv1(i,j)-& !u2*d(delv1)/dy & dx
!            newton_grade_n1(i,j,1)+ &      !delv1*d(v1)/dx
            up(i,j,2)*(zeta_n(i,j)+fcor(i,j))-& !delv2*(w+f)
            ulatlon(i,j,2)*newton_vortdelv2(i,j)+& !v2*d(delv2)/dx &dy
!            newton_grade_n2(i,j,1)+ &      !-delv2*d(v2)/dx
            grade_n3(i,j,1)) ! + d(delp)/dx



           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
             ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
!             newton_grade_n2(i,j,2)+&  !delu1*d(u2)/dy
             up(i,j,2)*dti +&                !delu2/dt
             ulatlon(i,j,1)*newton_vortdelv2(i,j)+& !u1*d(delu2)/dy & dx
!            newton_grade_n1(i,j,2)+ &      !delu1*d(u1)/dy
             grade_n3(i,j,2))!d(delp)/dy


         ptens(i,j,k,ie) = spheremp(i,j)*(&
             dh(i,j)*dti+&                         !delp/dt
             !divpv1(i,j)+divpv2(i,j)+&
             ( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )*&
             (divpv1(i,j)+divpv2(i,j))+&
             divdpv(i,j))

            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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
              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
            end if
            if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
          end do
        end do
      end do
    end do !ie
  end do

  end subroutine sw_picard_op

  subroutine sw_picard_alphaschur_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_alphaschur')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only :  np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    real (c_double)       :: fxtemp(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(pv,pv)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dt,dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv
    real (kind=real_kind) ::  alpha




    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx


alpha=1.0d0

!    call t_startf('sw jacobian op')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    dt         = fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

    vtens=0.0d0
    ptens=0.0d0

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

  do ie=ns,ne

      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np

     do k=1,nlev

     grade_n3 = gradient_sphere(fptr%base(ie)%state%p(:,:,k,np1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

     do j=1,np
        do i=1,np
!vtens is in latlon coordinates
          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
            grade_n3(i,j,1)) ! + d(delp)/dx
           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             grade_n3(i,j,2))!d(delp)/dy
        end do !np
      end do !np

     end do !nlev

      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge2, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie

   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge2)
   !$OMP BARRIER

     do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge2, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

     end do !ie


!!apply rspheremp
!
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
                vtens(i,j,1,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
            end do
          end do
       end do
     end do !ie
! end of Bt cal


       ! ===========================================================
       ! Compute velocity and pressure tendencies for all levels
       ! ===========================================================
!F_diag_inv(B'p) (B'p is in vtens and matches along element interfaces so F_diag_inv can be applied locally )

!  ApplyFDiagInv (the long way)
      lx = 0

    do n=1,nvar
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
              lx = lx+1
           if (topology == "cube" .and. test_case=="swtc1") then
!should we take this out of the preconditioner or leave this in?
              if (n==1) fxtemp(lx) = 0.0
              if (n==2) fxtemp(lx) = 0.0

           else
              !if(n==1) fxtemp(lx)=dt*fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              !if(n==2) fxtemp(lx)=dt*fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
              if(n==1) fxtemp(lx)=dt*vtens(i,j,1,k,ie)
              if(n==2) fxtemp(lx)=dt*vtens(i,j,2,k,ie)

           end if
              if (n==3) fxtemp(lx) =0.0
            end do
          end do
       end do
     end do !ie
     end do

  lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = fxtemp(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = fxtemp(lx)
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

! we now need to apply B to vtens, as well as G to dh (a.k.a. fptr%base(ie)%state%p(:,:,k,np1), and sum these two quantities
! now compute Gdp -BF_diag_inv(B'p) 
!
     do ie=ns,ne
         do j=1,np
           do i=1,np
             spheremp(i,j)=fptr%base(ie)%spheremp(i,j)
             rspheremp(i,j)=fptr%base(ie)%rspheremp(i,j)
           end do !np
         end do !np

       do k=1,nlev
!
!         ! ==============================================
!         ! Compute kinetic energy term at each time level
!         ! ==============================================
!
         do j=1,np
           do i=1,np
! set u
             ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
             ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
             up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
             up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  

             dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
             dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)

             pv1(i,j,1) = up(i,j,1)
             pv1(i,j,2) = 0.d00
             pv2(i,j,1) = 0.d00
             pv2(i,j,2) = up(i,j,2)

             dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )  
             dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )   

        !  pv1(i,j,1) = vtens(i,j,1,k,ie)
        !  pv1(i,j,2) = 0.d00
!
!          pv2(i,j,1) = 0.d00
!          pv2(i,j,2) = vtens(i,j,2,k,ie)

             dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
             dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )
           end do !np
         end do !np
!! create operators for full residual calc
!
     divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divdpv = divergence_sphere(dpv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

         ptens(i,j,k,ie) = spheremp(i,j)*(&
             !alpha*G
             alpha*( &
             dh(i,j)*dti+&                         !delp/dt
             divdpv(i,j)) + &
             !(1-alpha)*B(inv(diagF))B'
             (1.0d0-alpha )*(  ( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )*&
             (divpv1(i,j)+divpv2(i,j))))

            end do !np
          end do !np
       end do !nlev

!
!      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge1, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
  end do !ie



   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge1)
   !$OMP BARRIER


     do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge1, ptens(1,1,1,ie), nlev, kptr, fptr%base(ie)%desc)
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
              if (n==1) fx(lx) = 0.0d0
              if (n==2) fx(lx) = 0.0d0
              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
!              !if (n==3) fx(lx) = ptens(i,j,k,ie) 
!              if (n==3) fx(lx) = 0.0d0
            end do
          end do
       end do
     end do !ie
     end do
!

  end subroutine sw_picard_alphaschur_op


  subroutine sw_picard_schur_op_clip(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_schur_clip')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case,tstep_type
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !real (c_double)       :: fxtemp(nelemd)
    real (c_double)       :: fxtemp(2*nelemd)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dt,dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

    call t_startf('Precon Schur')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    dt         = fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

    if (tstep_type==12) then !crank nicholson
            dti=2*dti !crank nicholson has a factor of 2 in the time step coefficient
    endif

vtens=0.0d0
ptens=0.0d0

    lx = 0
!     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
!             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
!             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
!             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
!     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
       end do !np
      end do !np



     do k=1,nlev

     grade_n3 = gradient_sphere(fptr%base(ie)%state%p(:,:,k,np1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
            grade_n3(i,j,1)) ! + d(delp)/dx

           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             grade_n3(i,j,2))!d(delp)/dy


            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge2, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge2)
   !$OMP BARRIER

     do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge2, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

     end do !ie


!!apply rspheremp
!
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
                vtens(i,j,1,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
            end do
          end do
       end do
     end do !ie
! end of Bt cal


       ! ===========================================================
       ! Compute velocity and pressure tendencies for all levels
       ! ===========================================================
!F_diag_inv(B'p) (B'p is in vtens and matches along element interfaces so F_diag_inv can be applied locally )

!  ApplyFDiagInv 
      lx = 0

    do n=1,2
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
              lx = lx+1
           if (topology == "cube" .and. test_case=="swtc1") then
              if (n==1) fxtemp(lx) = 0.0
              if (n==2) fxtemp(lx) = 0.0

           else
              if(n==1) fxtemp(lx)=dt*vtens(i,j,1,k,ie)
              if(n==2) fxtemp(lx)=dt*vtens(i,j,2,k,ie)
           end if
            end do
          end do
       end do
     end do !ie
     end do



  lx = 0
     do n=1,2
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = fxtemp(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = fxtemp(lx)
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar





! we now need to apply B to vtens, as well as G to dh (a.k.a. fptr%base(ie)%state%p(:,:,k,np1), and sum these two quantities
! now compute Gdp -BF_diag_inv(B'p) 
!
   do ie=ns,ne
     do j=1,np
       do i=1,np
         spheremp(i,j)=fptr%base(ie)%spheremp(i,j)
         rspheremp(i,j)=fptr%base(ie)%rspheremp(i,j)
       end do !np
     end do !np

     do k=1,nlev
!
!         ! ==============================================
!         ! Compute kinetic energy term at each time level
!         ! ==============================================
!
      do j=1,np
       do i=1,np


! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  


          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)


          pv1(i,j,1) = up(i,j,1)
          pv1(i,j,2) = 0.d00
          pv2(i,j,1) = 0.d00
          pv2(i,j,2) = up(i,j,2)



          dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )  
          dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )   


          dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
          dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )
       end do !np
      end do !np
!! create operators for full residual calc
!
     divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divdpv = divergence_sphere(dpv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
!
!
!          ! ==============================================
!          ! Compute residual terms
!          ! ==============================================
!
!
!
      do j=1,np
        do i=1,np
!
         ptens(i,j,k,ie) = spheremp(i,j)*(&
             !G
             dh(i,j)*dti+&                         !delp/dt
             divdpv(i,j) - &
             !-BdiagFinvB'
             (  ( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )*&
             (divpv1(i,j)+divpv2(i,j))))

            end do !np
          end do !np
       end do !nlev

!
!      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge1, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
  end do !ie



   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge1)
   !$OMP BARRIER


     do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge1, ptens(1,1,1,ie), nlev, kptr, fptr%base(ie)%desc)
     end do !ie

       ! ===========================================================
       ! Compute velocity and pressure tendencies for all levels
       ! ===========================================================

      lx = 0

!    do n=1,nvar
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
              lx = lx+1
              
        
!              if (n==1) fx(lx) = 0.0d0
!              if (n==2) fx(lx) = 0.0d0
              fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
            end do
          end do
       end do
     end do !ie
!     end do
!


    call t_stopf('Precon Schur')


  end subroutine sw_picard_schur_op_clip



  subroutine sw_picard_schur_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_schur')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    real (c_double)       :: fxtemp(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dt,dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    dt         = fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0


    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

     grade_n3 = gradient_sphere(fptr%base(ie)%state%p(:,:,k,np1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!vtens is in latlon coordinates

          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
            grade_n3(i,j,1)) ! + d(delp)/dx

           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             grade_n3(i,j,2))!d(delp)/dy


            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge2, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge2)
   !$OMP BARRIER

     do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge2, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

     end do !ie


!!apply rspheremp
!
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
                vtens(i,j,1,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
            end do
          end do
       end do
     end do !ie
! end of Bt cal


       ! ===========================================================
       ! Compute velocity and pressure tendencies for all levels
       ! ===========================================================
!F_diag_inv(B'p) (B'p is in vtens and matches along element interfaces so F_diag_inv can be applied locally )

!  ApplyFDiagInv (the long way)
      lx = 0

    do n=1,nvar
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
              lx = lx+1
           if (topology == "cube" .and. test_case=="swtc1") then
!should we take this out of the preconditioner or leave this in?
              if (n==1) fxtemp(lx) = 0.0
              if (n==2) fxtemp(lx) = 0.0

           else
              !if(n==1) fxtemp(lx)=dt*fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              !if(n==2) fxtemp(lx)=dt*fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
              if(n==1) fxtemp(lx)=dt*vtens(i,j,1,k,ie)
              if(n==2) fxtemp(lx)=dt*vtens(i,j,2,k,ie)

           end if
              if (n==3) fxtemp(lx) =0.0
            end do
          end do
       end do
     end do !ie
     end do



  lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = fxtemp(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = fxtemp(lx)
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar





! we now need to apply B to vtens, as well as G to dh (a.k.a. fptr%base(ie)%state%p(:,:,k,np1), and sum these two quantities
! now compute Gdp -BF_diag_inv(B'p) 
!
   do ie=ns,ne
     do j=1,np
       do i=1,np
         spheremp(i,j)=fptr%base(ie)%spheremp(i,j)
         rspheremp(i,j)=fptr%base(ie)%rspheremp(i,j)
       end do !np
     end do !np

     do k=1,nlev
!
!         ! ==============================================
!         ! Compute kinetic energy term at each time level
!         ! ==============================================
!
      do j=1,np
       do i=1,np


! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  


          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)


          pv1(i,j,1) = up(i,j,1)
          pv1(i,j,2) = 0.d00
          pv2(i,j,1) = 0.d00
          pv2(i,j,2) = up(i,j,2)



          dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )  
          dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )   

        !  pv1(i,j,1) = vtens(i,j,1,k,ie)
        !  pv1(i,j,2) = 0.d00
!
!          pv2(i,j,1) = 0.d00
!          pv2(i,j,2) = vtens(i,j,2,k,ie)

          dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
          dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )
       end do !np
      end do !np
!! create operators for full residual calc
!
     divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divdpv = divergence_sphere(dpv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
!
!
!          ! ==============================================
!          ! Compute residual terms
!          ! ==============================================
!
!
!
      do j=1,np
        do i=1,np
!
         ptens(i,j,k,ie) = spheremp(i,j)*(&
             !G
             dh(i,j)*dti+&                         !delp/dt
             divdpv(i,j) - &
             !-BdiagFinvB'
             (  ( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )*&
             (divpv1(i,j)+divpv2(i,j))))

            end do !np
          end do !np
       end do !nlev

!
!      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge1, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
  end do !ie



   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge1)
   !$OMP BARRIER


     do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge1, ptens(1,1,1,ie), nlev, kptr, fptr%base(ie)%desc)
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
              
        
              if (n==1) fx(lx) = 0.0d0
              if (n==2) fx(lx) = 0.0d0
              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
!              !if (n==3) fx(lx) = ptens(i,j,k,ie) 
!              if (n==3) fx(lx) = 0.0d0
            end do
          end do
       end do
     end do !ie
     end do
!




  end subroutine sw_picard_schur_op












  subroutine sw_picard_FDFinvBt2_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_FDFinvBt2')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    real (c_double)       :: fxtemp(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dt,dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    dt         = fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0


    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

     grade_n3 = gradient_sphere(fptr%base(ie)%state%p(:,:,k,np1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!vtens is in latlon coordinates

          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
            grade_n3(i,j,1)) ! + d(delp)/dx

           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             grade_n3(i,j,2))!d(delp)/dy


            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge2, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge2)
   !$OMP BARRIER

     do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge2, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

     end do !ie


!!apply rspheremp
!
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
                vtens(i,j,1,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
            end do
          end do
       end do
     end do !ie
! end of Bt cal


       ! ===========================================================
       ! Compute velocity and pressure tendencies for all levels
       ! ===========================================================
!F_diag_inv(B'p) (B'p is in vtens and matches along element interfaces so F_diag_inv can be applied locally )

!  ApplyFDiagInv (the long way)
      lx = 0

    do n=1,nvar
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
              lx = lx+1
           if (topology == "cube" .and. test_case=="swtc1") then
!should we take this out of the preconditioner or leave this in?
              if (n==1) fxtemp(lx) = 0.0
              if (n==2) fxtemp(lx) = 0.0

           else
              if(n==1) fxtemp(lx)=dt*vtens(i,j,1,k,ie)
              if(n==2) fxtemp(lx)=dt*vtens(i,j,2,k,ie)

           end if
              if (n==3) fxtemp(lx) =0.0
            end do
          end do
       end do
     end do !ie
     end do



  lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = fxtemp(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = fxtemp(lx)
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar





! we now need to apply B to vtens, as well as G to dh (a.k.a. fptr%base(ie)%state%p(:,:,k,np1), and sum these two quantities
! now compute Gdp -BF_diag_inv(B'p) 
!
   do ie=ns,ne
     do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
         spheremp(i,j)=fptr%base(ie)%spheremp(i,j)
         rspheremp(i,j)=fptr%base(ie)%rspheremp(i,j)
       end do !np
     end do !np







    do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  


          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)


       end do !np
      end do !np

! create operators for full residual calc

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!picard linearization of 11 block 
!note this 11 block is representing a 2x2 system for up1 and up2 

!          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
!             up(i,j,1)*dti )
!
!           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
!             up(i,j,2)*dti )
!


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*dti - &                 !delv1/dt
            ulatlon(i,j,2)*newton_vortdelv1(i,j)-& !u2*d(delv1)/dy & dx
            up(i,j,2)*(zeta_n(i,j)+fcor(i,j))-& !delv2*(w+f)
            ulatlon(i,j,2)*newton_vortdelv2(i,j)) !v2*d(delv2)/dx &dy


           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
             ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
             up(i,j,2)*dti +&                !delu2/dt
             ulatlon(i,j,1)*newton_vortdelv2(i,j)) !u1*d(delu2)/dy & dx



            end do !np
          end do !np
       end do !nlev



      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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
              
              if (n==1) fx(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie) 
              if (n==2) fx(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie) 
              !if (n==1) fx(lx) = vtens(i,j,1,k,ie) 
              !if (n==2) fx(lx) = vtens(i,j,2,k,ie) 
              if (n==3) fx(lx) = 0.0d0
!              !if (n==3) fx(lx) = ptens(i,j,k,ie) 
!              if (n==3) fx(lx) = 0.0d0
            end do
          end do
       end do
     end do !ie
     end do
!




  end subroutine sw_picard_FDFinvBt2_op


 subroutine sw_picard_DFinvBt_op_clip(xs, nelemd,fx,nret, c_ptr_to_object) bind(C,name='sw_picard_DFinvBt_clip')

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
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    integer(c_int) ,intent(in) ,value  :: nret
    real (c_double) ,intent(out)       :: fx(nret)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dt,dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

    call t_startf('Precon DFinvBt_clip')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    dt         = fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0


    lx = 0
!     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
      !       if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
      !       if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
      !       if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
!     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
       end do !np
      end do !np



     do k=1,nlev

     grade_n3 = gradient_sphere(fptr%base(ie)%state%p(:,:,k,np1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
            grade_n3(i,j,1)) ! + d(delp)/dx

           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             grade_n3(i,j,2))!d(delp)/dy


            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge2, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge2)
   !$OMP BARRIER

     do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge2, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

     end do !ie


!!apply rspheremp
!
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
                vtens(i,j,1,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

            end do
          end do
       end do
     end do !ie
! end of Bt cal


       ! ===========================================================
       ! Compute velocity and pressure tendencies for all levels
       ! ===========================================================


      lx = 0

    do n=1,2
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
              lx = lx+1
           if (topology == "cube" .and. test_case=="swtc1") then
              if (n==1) fx(lx) = 0.0
              if (n==2) fx(lx) = 0.0

           else

              if(n==1) fx(lx)=dt*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=dt*vtens(i,j,2,k,ie)
           end if

            end do
          end do
       end do
     end do !ie
     end do

    call t_stopf('Precon DFinvBt_clip')

  end subroutine sw_picard_DFinvBt_op_clip




  subroutine sw_picard_DFinvBt_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_DFinvBt')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dt,dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    dt         = fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0


    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

     grade_n3 = gradient_sphere(fptr%base(ie)%state%p(:,:,k,np1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!vtens is in latlon coordinates

          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
            grade_n3(i,j,1)) ! + d(delp)/dx

           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             grade_n3(i,j,2))!d(delp)/dy


            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge2, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge2)
   !$OMP BARRIER

     do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge2, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

     end do !ie


!!apply rspheremp
!
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
                vtens(i,j,1,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

              !  vtens(i,j,1,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              !  vtens(i,j,2,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
            end do
          end do
       end do
     end do !ie
! end of Bt cal


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
!should we take this out of the preconditioner or leave this in?
              if (n==1) fx(lx) = 0.0
              if (n==2) fx(lx) = 0.0

           else

              if(n==1) fx(lx)=dt*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=dt*vtens(i,j,2,k,ie)

           end if
              if (n==3) fx(lx) =0.0
            end do
          end do
       end do
     end do !ie
     end do


  end subroutine sw_picard_DFinvBt_op






  subroutine sw_picard_FDFinvBt_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_FDFinvBt')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dt,dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    dt         = fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0

!write(6,*)'dti=',dti,'\n'
!xs=1.0d0
!write(6,*)'xs=[',xs,']\n'

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

!      do j=1,np
!       do i=1,np
!
!          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)
!
!       end do !np
!      end do !np

! create operators for full residual calc
     !grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector        (grad)del phi
     grade_n3 = gradient_sphere(fptr%base(ie)%state%p(:,:,k,np1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi




          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!vtens is in latlon coordinates

          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
            grade_n3(i,j,1)) ! + d(delp)/dx

           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             grade_n3(i,j,2))!d(delp)/dy


            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge2, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge2)
   !$OMP BARRIER

     do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge2, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

     end do !ie




!!apply rspheremp

    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
                vtens(i,j,1,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
            end do
          end do
       end do
     end do !ie
!!
!end of B' application


!apply inverse of diag of F, which we approximate as dt*inv(Mass)

    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
                vtens(i,j,1,k,ie)=dt*fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie)=dt*fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
            end do
          end do
       end do
     end do !ie



  do ie=ns,ne
   do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np


          delv1(i,j,1)=vtens(i,j,1,k,ie)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=vtens(i,j,2,k,ie)


       end do !np
      end do !np

! create operators for full residual calc

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!picard linearization of 11 block 

!
          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
            vtens(i,j,1,k,ie)*dti - &                 !delv1/dt
            ulatlon(i,j,2)*newton_vortdelv1(i,j)-& !u2*d(delv1)/dy & dx
            vtens(i,j,2,k,ie)*(zeta_n(i,j)+fcor(i,j))-& !delv2*(w+f)
            ulatlon(i,j,2)*newton_vortdelv2(i,j)) !v2*d(delv2)/dx &dy



           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             vtens(i,j,1,k,ie)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
             ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
             vtens(i,j,2,k,ie)*dti +&                !delu2/dt
             ulatlon(i,j,1)*newton_vortdelv2(i,j)) !u1*d(delu2)/dy & dx


!test
!          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
!             vtens(i,j,1,k,ie)*dti )                  !delv1/dt
!           
!          vtens(i,j,2,k,ie)= spheremp(i,j)*( &
!             vtens(i,j,2,k,ie)*dti )                !delu2/dt
!


            end do !np
          end do !np
       end do !nlev





   kptr=0
       call edgeVpack(fptr%edge2, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
!       kptr=0
!       call edgerotate(pptr%edge2,2*nlev,kptr,pptr%base(ie)%desc)
    end do !ie

    !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge2)
    !$OMP BARRIER
    do ie=ns,ne

       kptr=0
       call edgeVunpack(fptr%edge2, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

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
!should we take this out of the preconditioner or leave this in?
              if (n==1) fx(lx) = 0.0
              if (n==2) fx(lx) = 0.0

           else
!check with Mark on the role of the rspheremp in this part of the code, he alluded to it during our conversation
! make sure this is consistent with the equations we are solving

              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

           end if
              if (n==3) fx(lx) =0.0
            end do
          end do
       end do
     end do !ie
     end do


  end subroutine sw_picard_FDFinvBt_op








  subroutine sw_picard_FBt_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_FBt')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only :  np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dt,dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    dt         = fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0


    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

!      do j=1,np
!       do i=1,np
!
!          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)
!
!       end do !np
!      end do !np

! create operators for full residual calc
     !grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector        (grad)del phi
     grade_n3 = gradient_sphere(fptr%base(ie)%state%p(:,:,k,np1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi




          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!vtens is in latlon coordinates

          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
            grade_n3(i,j,1)) ! + d(delp)/dx

           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             grade_n3(i,j,2))!d(delp)/dy


            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge2, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge2)
   !$OMP BARRIER

     do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge2, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

     end do !ie









!!apply rspheremp
!
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
                vtens(i,j,1,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
            end do
          end do
       end do
     end do !ie
!


  do ie=ns,ne
   do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np


          delv1(i,j,1)=vtens(i,j,1,k,ie)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=vtens(i,j,2,k,ie)


       end do !np
      end do !np

! create operators for full residual calc

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!picard linearization of 11 block 


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
            vtens(i,j,1,k,ie)*dti - &                 !delv1/dt
            ulatlon(i,j,2)*newton_vortdelv1(i,j)-& !u2*d(delv1)/dy & dx
            vtens(i,j,2,k,ie)*(zeta_n(i,j)+fcor(i,j))-& !delv2*(w+f)
            ulatlon(i,j,2)*newton_vortdelv2(i,j)) !v2*d(delv2)/dx &dy



           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             vtens(i,j,1,k,ie)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
             ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
             vtens(i,j,2,k,ie)*dti +&                !delu2/dt
             ulatlon(i,j,1)*newton_vortdelv2(i,j)) !u1*d(delu2)/dy & dx


            end do !np
          end do !np
       end do !nlev





   kptr=0
       call edgeVpack(fptr%edge2, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
!       kptr=0
!       call edgerotate(pptr%edge2,2*nlev,kptr,pptr%base(ie)%desc)
    end do !ie

    !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge2)
    !$OMP BARRIER
    do ie=ns,ne

       kptr=0
       call edgeVunpack(fptr%edge2, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

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
!should we take this out of the preconditioner or leave this in?
              if (n==1) fx(lx) = 0.0
              if (n==2) fx(lx) = 0.0

           else
              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

           end if
              if (n==3) fx(lx) =0.0
            end do
          end do
       end do
     end do !ie
     end do


  end subroutine sw_picard_FBt_op







  subroutine sw_picard_block_FDiagInv(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_block_FDiagInv')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only :  np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dt, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dt         = fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0

!write(6,*)'dti=',dti,'\n'
!xs=1.0d0
!write(6,*)'xs=[',xs,']\n'
!write(6,*)'nvar=',nvar,';'
!write(6,*)'ne=',ne,';'
!write(6,*)'nlev=',nlev,';'
!write(6,*)'nv=',nv,';'

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  


          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)


       end do !np
      end do !np

! create operators for full residual calc

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!Diag F inv


!          vtens(i,j,1,k,ie)= ( &
!             up(i,j,1)*fptr%dt)                  !delv1*dt
!
!           vtens(i,j,2,k,ie)= ( &
!             up(i,j,2)*fptr%dt)                !delu2*dt



          vtens(i,j,1,k,ie)= (rspheremp(i,j)*( &
             up(i,j,1)*dt )          )       !delv1/dt
             !up(i,j,1) )          )       !delv1/dt
!
           vtens(i,j,2,k,ie)=(rspheremp(i,j)*( &
             up(i,j,2)*dt)           )     !delu2/dt
             !up(i,j,2))           )     !delu2/dt





            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
!       kptr=0
!       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
!       kptr=nlev
!       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


   !$OMP BARRIER
!    call bndry_exchangeV(fptr%hybrid,fptr%edge3)
   !$OMP BARRIER

!    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
!       kptr=0
!       call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, fptr%base(ie)%desc)
!       kptr=nlev
!       call edgeVunpack(fptr%edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

!     end do !ie
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

              if(n==1) fx(lx)=vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=vtens(i,j,2,k,ie)


           end if
              if (n==3) fx(lx) =0.0 
            end do
          end do
       end do
     end do !ie
     end do



  end subroutine sw_picard_block_FDiagInv







  subroutine sw_picard_block_FDiag(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_block_FDiag')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0

!write(6,*)'dti=',dti,'\n'
!xs=1.0d0
!write(6,*)'xs=[',xs,']\n'
!write(6,*)'nvar=',nvar,';'
!write(6,*)'ne=',ne,';'
!write(6,*)'nlev=',nlev,';'
!write(6,*)'nv=',nv,';'

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  


          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)


       end do !np
      end do !np

! create operators for full residual calc

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!picard linearization of 11 block 

      !    vtens(i,j,1,k,ie)= ( &
      !       up(i,j,1)*dti )                 !delv1/dt
!
!           vtens(i,j,2,k,ie)= ( &
!             up(i,j,2)*dti)                !delu2/dt



          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*dti )                 !delv1/dt
           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             up(i,j,2)*dti)                !delu2/dt


            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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

              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
              !if(n==1) fx(lx)=vtens(i,j,1,k,ie)
              !if(n==2) fx(lx)=vtens(i,j,2,k,ie)

           end if
              if (n==3) fx(lx) =0.0 
            end do
          end do
       end do
     end do !ie
     end do



  end subroutine sw_picard_block_FDiag








  subroutine sw_picard_block_11alt(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_block_11alt')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0

!write(6,*)'dti=',dti,'\n'
!xs=1.0d0
!write(6,*)'xs=[',xs,']\n'
!write(6,*)'nvar=',nvar,';'
!write(6,*)'ne=',ne,';'
!write(6,*)'nlev=',nlev,';'
!write(6,*)'np=',np,';'

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
!write(6,*)'fcordiag(',ie,',',i,',',j,')=',fcor(i,j),'\n'
!call flush(6)

        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  


          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)


       end do !np
      end do !np

! create operators for full residual calc

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!picard linearization of 11 block 



          vtens(i,j,1,k,ie)= 0.0d0

!          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
!             up(i,j,1)*dti + &                 !delv1/dt
!            !-ulatlon(i,j,2)*newton_vortdelv1(i,j)+& !u2*d(delv1)/dy & dx
!            -up(i,j,2)*(zeta_n(i,j)+fcor(i,j))) !delv2*(w+f)
!            !-up(i,j,2)*(fcor(i,j))) !delv2*(w+f)
!            !-up(i,j,2)*(zeta_n(i,j)+fcor(i,j))+& !delv2*(w+f)
!            !-ulatlon(i,j,2)*newton_vortdelv2(i,j)) !v2*d(delv2)/dx &dy

!           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
!             up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
!            ! up(i,j,1)*(fcor(i,j))+& !delu1*(w+f)
!             !ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
!             up(i,j,2)*dti )                !delu2/dt
!             !up(i,j,2)*dti +&                !delu2/dt
!             !ulatlon(i,j,1)*newton_vortdelv2(i,j)) !u1*d(delu2)/dy & dx

!full picard
!          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
!             up(i,j,1)*dti + &                 !delv1/dt
!            -ulatlon(i,j,2)*newton_vortdelv1(i,j)+& !u2*d(delv1)/dy & dx
!            -up(i,j,2)*(zeta_n(i,j)+fcor(i,j))+& !delv2*(w+f)
!            -ulatlon(i,j,2)*newton_vortdelv2(i,j)) !v2*d(delv2)/dx &dy
!
!           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
!             up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
!             ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
!             up(i,j,2)*dti +&                !delu2/dt
!             ulatlon(i,j,1)*newton_vortdelv2(i,j)) !u1*d(delu2)/dy & dx
!end full picard

          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*dti - &                 !delv1/dt
         !   -ulatlon(i,j,2)*newton_vortdelv1(i,j)+& !u2*d(delv1)/dy & dx
            !-up(i,j,2)*(zeta_n(i,j)+fcor(i,j))) !delv2*(w+f)
            up(i,j,2)*(fcor(i,j))) !delv2*(w+f)
         !   -ulatlon(i,j,2)*newton_vortdelv2(i,j)) !v2*d(delv2)/dx &dy

           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             !up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
             up(i,j,1)*(fcor(i,j))+& !delu1*(w+f)
         !    ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
             up(i,j,2)*dti )                !delu2/dt
         !    ulatlon(i,j,1)*newton_vortdelv2(i,j)) !u1*d(delu2)/dy & dx







!          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
!             up(i,j,1)*dti + &                 !delv1/dt
!            !-ulatlon(i,j,2)*newton_vortdelv1(i,j)+& !u2*d(delv1)/dy & dx
!            -up(i,j,2)*(zeta_n(i,j)+fcor(i,j))) !delv2*(w+f)
!            !-up(i,j,2)*(fcor(i,j))) !delv2*(w+f)
!            !-up(i,j,2)*(zeta_n(i,j)+fcor(i,j))+& !delv2*(w+f)
!            !-ulatlon(i,j,2)*newton_vortdelv2(i,j)) !v2*d(delv2)/dx &dy
!
!           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
!             up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
!            ! up(i,j,1)*(fcor(i,j))+& !delu1*(w+f)
!             !ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
!             up(i,j,2)*dti )                !delu2/dt
!             !up(i,j,2)*dti +&                !delu2/dt
!             !ulatlon(i,j,1)*newton_vortdelv2(i,j)) !u1*d(delu2)/dy & dx


            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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

              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

           end if
              if (n==3) fx(lx) =0.0 
            end do
          end do
       end do
     end do !ie
     end do



  end subroutine sw_picard_block_11alt






 subroutine sw_picard_block_11_clip(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_block_11_clip')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case,tstep_type
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

    call t_startf('precon 11')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5
    
    if (tstep_type==12) then !crank nicholson

            dti=2*dti !crank nicholson has a factor of 2 in the time step coefficient
    endif


vtens=0.0d0
ptens=0.0d0


    lx = 0
     !do n=1,nvar
     do n=1,2 !only need the u and v blocks
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
            ! if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  


          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)


       end do !np
      end do !np

! create operators for full residual calc

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!picard linearization of 11 block 
!note this 11 block is representing a 2x2 system for up1 and up2 


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
            up(i,j,1)*dti - &                 !delv1/dt
            ulatlon(i,j,2)*newton_vortdelv1(i,j)-& !u2*d(delv1)/dy & dx
            up(i,j,2)*(zeta_n(i,j)+fcor(i,j))-& !delv2*(w+f)
            ulatlon(i,j,2)*newton_vortdelv2(i,j)) !v2*d(delv2)/dx &dy

           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
             ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
             up(i,j,2)*dti +&                !delu2/dt
             ulatlon(i,j,1)*newton_vortdelv2(i,j)) !u1*d(delu2)/dy & dx



            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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

    do n=1,2 !only u and v blocks
    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
              lx = lx+1
           if (topology == "cube" .and. test_case=="swtc1") then
              if (n==1) fx(lx) = 0.0
              if (n==2) fx(lx) = 0.0
              
           else 

              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

           end if
            end do
          end do
       end do
     end do !ie
     end do


    call t_stopf('precon 11')

  end subroutine sw_picard_block_11_clip

















  subroutine sw_picard_block_11(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_block_11')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0

!write(6,*)'dti=',dti,'\n'
!xs=1.0d0
!write(6,*)'xs=[',xs,']\n'
!write(6,*)'nvar=',nvar,';'
!write(6,*)'ne=',ne,';'
!write(6,*)'nlev=',nlev,';'
!write(6,*)'nv=',nv,';'

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  


          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)


       end do !np
      end do !np

! create operators for full residual calc

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!picard linearization of 11 block 
!note this 11 block is representing a 2x2 system for up1 and up2 


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
            up(i,j,1)*dti - &                 !delv1/dt
            ulatlon(i,j,2)*newton_vortdelv1(i,j)-& !u2*d(delv1)/dy & dx
            up(i,j,2)*(zeta_n(i,j)+fcor(i,j))-& !delv2*(w+f)
            ulatlon(i,j,2)*newton_vortdelv2(i,j)) !v2*d(delv2)/dx &dy

           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
             ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
             up(i,j,2)*dti +&                !delu2/dt
             ulatlon(i,j,1)*newton_vortdelv2(i,j)) !u1*d(delu2)/dy & dx



            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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

              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

           end if
              if (n==3) fx(lx) =0.0 
            end do
          end do
       end do
     end do !ie
     end do



  end subroutine sw_picard_block_11


  subroutine sw_picard_block_12(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_block_12')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5d0

vtens=0.0d0
ptens=0.0d0

!write(6,*)'dti=',dti,'\n'
!xs=1.0d0
!write(6,*)'xs=[',xs,']\n'

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  


          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)



! set dh_n=p_notes=g(h+hmean)=phi
          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)


          pv1(i,j,1) = up(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          pv1(i,j,2) = 0.d00

          pv2(i,j,1) = 0.d00
          pv2(i,j,2) = up(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )

          dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
          dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )

       end do !np
      end do !np

! create operators for full residual calc
     grade_n1 = gradient_sphere(up(:,:,1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_1
     grade_n2 = gradient_sphere(up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_2
     grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi

     !newton_grade_n1 = gradient_sphere(ulatlon(:,:,1)*up(:,:,1) ,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_1
     !newton_grade_n2 = gradient_sphere(ulatlon(:,:,2)*up(:,:,2),fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_2
     !newton_grade_n3 = gradient_sphere(dh_n,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector        (grad)phi

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w

     divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     !divdpv = divergence_sphere(dpv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
            grade_n3(i,j,1)) ! + d(delp)/dx


           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             grade_n3(i,j,2))!d(delp)/dy



            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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
              if (n==1) fx(lx) = 0.0d0
              if (n==2) fx(lx) = 0.0d0
              
           else 

              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

           end if
              if (n==3) fx(lx) = 0.0d0
            end do
          end do
       end do
     end do !ie
     end do


  end subroutine sw_picard_block_12




  subroutine sw_picard_block_21_clip(xs, nelemd,fx, nret, c_ptr_to_object) bind(C,name='sw_picard_block_21_clip')

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
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    integer(c_int) ,intent(in) ,value  :: nret
    real (c_double) ,intent(out)       :: fx(nret)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx


    call t_startf('picard 21_clip')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0


    lx = 0
     do n=1,2
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  



! set dh_n=p_notes=g(h+hmean)=phi
!          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
!          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)


          pv1(i,j,1) = up(i,j,1)
          pv1(i,j,2) = 0.d00
          pv2(i,j,1) = 0.d00
          pv2(i,j,2) = up(i,j,2)


       end do !np
      end do !np

! create operators for full residual calc

     divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

         ptens(i,j,k,ie) = spheremp(i,j)*(&
             ( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )*&
             (divpv1(i,j)+divpv2(i,j)))


            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
!       kptr=nlev
!       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge3)
   !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, fptr%base(ie)%desc)
!       kptr=nlev
!       call edgeVunpack(fptr%edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)

     end do !ie
      ! ===========================================================
       ! Compute velocity and pressure tendencies for all levels
       ! ===========================================================

      lx = 0

    do ie=ns,ne
      do k=1,nlev
         do j=1,np
            do i=1,np
              lx = lx+1
!           if (topology == "cube" .and. test_case=="swtc1") then
!              if (n==1) fx(lx) = 0.0
!              if (n==2) fx(lx) = 0.0
!              
!           else 
!
!             if (n==1) fx(lx) = 0.0
!             if (n==2) fx(lx) = 0.0
!           end if
              !if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
               fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
            end do
          end do
       end do
     end do !ie

    call t_stopf('picard 21_clip')

  end subroutine sw_picard_block_21_clip




  subroutine sw_picard_block_21(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_block_21')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0

!write(6,*)'dti=',dti,'\n'
!xs=1.0d0
!write(6,*)'xs=[',xs,']\n'

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  



! set dh_n=p_notes=g(h+hmean)=phi
          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)


          !pv1(i,j,1) = up(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          !pv1(i,j,2) = 0.d00
          !pv2(i,j,1) = 0.d00
          !pv2(i,j,2) = up(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
!

          pv1(i,j,1) = up(i,j,1)
          pv1(i,j,2) = 0.d00
          pv2(i,j,1) = 0.d00
          pv2(i,j,2) = up(i,j,2)


       end do !np
      end do !np

! create operators for full residual calc

     divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     !divdpv = divergence_sphere(dpv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

         ptens(i,j,k,ie) = spheremp(i,j)*(&
             ( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )*&
             (divpv1(i,j)+divpv2(i,j)))




!         ptens(i,j,k,ie) = spheremp(i,j)*(&
!             divpv1(i,j)+divpv2(i,j))

            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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

             if (n==1) fx(lx) = 0.0
             if (n==2) fx(lx) = 0.0
           end if
              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
            end do
          end do
       end do
     end do !ie
     end do


  end subroutine sw_picard_block_21


  subroutine sw_picard_block_22(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_block_22')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
!          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
!          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  



          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)

          dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
          dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )

       end do !np
      end do !np



     divdpv = divergence_sphere(dpv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

         ptens(i,j,k,ie) = spheremp(i,j)*(&
             dh(i,j)*dti+&                         !delp/dt
             divdpv(i,j))


!         ptens(i,j,k,ie) = spheremp(i,j)*(&
!             dh(i,j)*dti)                         !delp/dt

            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge1, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
  end do !ie


   !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge1)
   !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge1, ptens(1,1,1,ie), nlev, kptr, fptr%base(ie)%desc)

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
              if (n==1) fx(lx) = 0.0
              if (n==2) fx(lx) = 0.0
              !if (n==3) fx(lx) = 0.0 !testing the saddlepoint matrix
              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
            end do
          end do
       end do
     end do !ie
     end do


  end subroutine sw_picard_block_22


















  subroutine sw_jacobian_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_jacobian')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0

!write(6,*)'dti=',dti,'\n'
!xs=1.0d0
!write(6,*)'xs=[',xs,']\n'

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   !
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   !
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  

          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)



!test
!         up(i,j,1)=v1p
!          up(i,j,2)=v2p
!end test


! where is ps???? shouldn't ps be included in dp? ps=g*h_s in comment above, ANSWER: NO, also no pmean, otherwise we would be adding these contributions at each nonlinear update 
!what is p 
!what is pmean
!ps= surface height g*h_s in notes


! set dh_n=p_notes=g(h+hmean)=phi
          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
!test unit vector
!          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)

! set dp
          !dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)+ps(i,j)+pmean
          !dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)+pmean
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)


          !pv1(i,j,1) = up(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          pv1(i,j,1) = up(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          pv1(i,j,2) = 0.d00

          pv2(i,j,1) = 0.d00
          !pv2(i,j,2) = up(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          pv2(i,j,2) = up(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )

          dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
          dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )

          !dpv(i,j,1) = ulatlon(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,np1))
          !dpv(i,j,2) = ulatlon(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,np1))


       end do !np
      end do !np

! create operators for full residual calc
     grade_n1 = gradient_sphere(up(:,:,1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_1
     grade_n2 = gradient_sphere(up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_2
     grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi
     !grade_n3 = gradient_sphere(dh+ps,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector        (grad)del phi
     !grade_n3 = gradient_sphere(dh+ps+pmean,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector        (grad)del phi


     !newton_grade_n1 = gradient_sphere(ulatlon(:,:,1),fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_1
     !newton_grade_n2 = gradient_sphere(ulatlon(:,:,2),fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_2

! .5*(uco(:,:,1)*upco(:,:,1)+upco(:,:,1)*uco(:,:,1))=uco(:,:,1)*upco(:,:,1)  .5*(u_1*delu^1 + delu_1*u^1)
! .5*(uco(:,:,2)*upco(:,:,2)+upco(:,:,2)*uco(:,:,2))=uco(:,:,2)*upco(:,:,2)  .5*(u_2*delu^2 + delu_2*u^2)
!     newton_grade_n1 = gradient_sphere(uco(:,:,1)*upco(:,:,1) ,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_1
!     newton_grade_n2 = gradient_sphere(uco(:,:,2)*upco(:,:,2),fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_2

! .5*(uco(:,:,1)*upco(:,:,1)+upco(:,:,1)*uco(:,:,1))=uco(:,:,1)*upco(:,:,1)  .5*(u_1*delu^1 + delu_1*u^1)
! .5*(uco(:,:,2)*upco(:,:,2)+upco(:,:,2)*uco(:,:,2))=uco(:,:,2)*upco(:,:,2)  .5*(u_2*delu^2 + delu_2*u^2)
     newton_grade_n1 = gradient_sphere(ulatlon(:,:,1)*up(:,:,1) ,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)u_1
     newton_grade_n2 = gradient_sphere(ulatlon(:,:,2)*up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)u_2




     newton_grade_n3 = gradient_sphere(dh_n,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)phi

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w

     divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divdpv = divergence_sphere(dpv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

! Jacobian


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*dti - &                 !delv1/dt
            ulatlon(i,j,2)*newton_vortdelv1(i,j)+& !u2*d(delv1)/dy & dx
            newton_grade_n1(i,j,1)- &      !delv1*d(v1)/dx
            up(i,j,2)*(zeta_n(i,j)+fcor(i,j))-& !delv2*(w+f)
            ulatlon(i,j,2)*newton_vortdelv2(i,j)+& !v2*d(delv2)/dx &dy
            newton_grade_n2(i,j,1)+ &      !-delv2*d(v2)/dx
            grade_n3(i,j,1)) ! + d(delp)/dx



           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
             ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
             newton_grade_n2(i,j,2)+&  !delu1*d(u2)/dy
             up(i,j,2)*dti +&                !delu2/dt
             ulatlon(i,j,1)*newton_vortdelv2(i,j)+& !u1*d(delu2)/dy & dx
            newton_grade_n1(i,j,2)+ &      !delu1*d(u1)/dy
             grade_n3(i,j,2))!d(delp)/dy


         ptens(i,j,k,ie) = spheremp(i,j)*(&
             dh(i,j)*dti+&                         !delp/dt
             divpv1(i,j)+divpv2(i,j)+divdpv(i,j))





            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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
              
!test unit base vector
           else 

              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

           end if
              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
            end do
          end do
       end do
     end do !ie
     end do


!write(6,*)'M=[',spheremp*rspheremp,']'
!    call t_stopf('sw jacobian op')

!  write(6,*)'fx=[',fx,']\n'
!fx=xs
!fx=base_vector

  end subroutine sw_jacobian_op



















  subroutine update_prec_state(xs, nelemd, c_ptr_to_object) bind(C,name='update_prec_state')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none

    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne
   
    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: nm1,n0,np1
    integer    :: lenx,lenp, lx

!write(6,*)'update prec state\n'
!    call t_startf('residual lin')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             !hack to put current point
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,n0) = xs(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,n0) = xs(lx)
             if (n==3) fptr%base(ie)%state%p(i,j,k,n0)  = xs(lx)
           !  write(6,*)'xs=',xs(lx),''
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar
 end subroutine update_prec_state



  subroutine get_jac_vector(xs, nelemd, c_ptr_to_object) bind(C,name='get_jac_vector')

    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only :  np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none

    real (c_double)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne


    integer    :: i,j,k,n,ie
    integer    :: nm1,n0,np1
    integer    :: lenx,lenp, lx

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0         = fptr%tl%n0
    ns         = fptr%nets
    ne         = fptr%nete

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             !hack to put current point
             if (n==1) xs(lx)= fptr%base(ie)%state%v(i,j,1,k,n0)
             if (n==2) xs(lx)= fptr%base(ie)%state%v(i,j,2,k,n0)
             if (n==3) xs(lx)= fptr%base(ie)%state%p(i,j,k,n0)  
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar
 end subroutine get_jac_vector








  subroutine sw_mass_test(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_mass_test')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only :  np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0

!write(6,*)'dti=',dti,'\n'
!xs=1.0d0
!write(6,*)'xs=[',xs,']\n'

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  


          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)



! set dh_n=p_notes=g(h+hmean)=phi
          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)


          pv1(i,j,1) = up(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          pv1(i,j,2) = 0.d00

          pv2(i,j,1) = 0.d00
          pv2(i,j,2) = up(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )

          dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
          dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )

       end do !np
      end do !np

! create operators for full residual calc
     grade_n1 = gradient_sphere(up(:,:,1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_1
     grade_n2 = gradient_sphere(up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_2
     grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi

     !newton_grade_n1 = gradient_sphere(ulatlon(:,:,1)*up(:,:,1) ,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_1
     !newton_grade_n2 = gradient_sphere(ulatlon(:,:,2)*up(:,:,2),fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_2
     !newton_grade_n3 = gradient_sphere(dh_n,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector        (grad)phi

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w

     divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     !divdpv = divergence_sphere(dpv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

! Mass


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
             up(i,j,1))


           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             up(i,j,2))

         ptens(i,j,k,ie) = spheremp(i,j)*(&
             dh(i,j))

            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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

              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
           end if
            end do
          end do
       end do
     end do !ie
     end do


  end subroutine sw_mass_test





  subroutine sw_picard_block_ID11(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_block_ID11')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0

!write(6,*)'dti=',dti,'\n'
!xs=1.0d0
!write(6,*)'xs=[',xs,']\n'
!write(6,*)'nvar=',nvar,';'
!write(6,*)'ne=',ne,';'
!write(6,*)'nlev=',nlev,';'
!write(6,*)'nv=',nv,';'

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  


          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)


       end do !np
      end do !np

! create operators for full residual calc

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!picard linearization of 11 block 


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
            up(i,j,1)*dti - &                 !delv1/dt
            ulatlon(i,j,2)*newton_vortdelv1(i,j)-& !u2*d(delv1)/dy & dx
            up(i,j,2)*(zeta_n(i,j)+fcor(i,j))-& !delv2*(w+f)
            ulatlon(i,j,2)*newton_vortdelv2(i,j)) !v2*d(delv2)/dx &dy



           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
             ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
             up(i,j,2)*dti +&                !delu2/dt
             ulatlon(i,j,1)*newton_vortdelv2(i,j)) !u1*d(delu2)/dy & dx


            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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
              if (n==1) fx(lx) = 0.0d0
              if (n==2) fx(lx) = 0.0d0
              
           else 
              if (n==1) fx(lx) = xs(lx)
              if (n==2) fx(lx) = xs(lx)

           end if
              if (n==3) fx(lx) =0.0d0 
              !if (n==3) fx(lx) =xs(lx)
            end do
          end do
       end do
     end do !ie
     end do



  end subroutine sw_picard_block_ID11













  subroutine sw_picard_block_ID22(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_block_ID22')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only :  np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0

!write(6,*)'dti=',dti,'\n'
!xs=1.0d0
!write(6,*)'xs=[',xs,']\n'
!write(6,*)'nvar=',nvar,';'
!write(6,*)'ne=',ne,';'
!write(6,*)'nlev=',nlev,';'
!write(6,*)'np=',np,';'

    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  


          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)


       end do !np
      end do !np

! create operators for full residual calc

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

!picard linearization of 11 block 


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*dti - &                 !delv1/dt
            ulatlon(i,j,2)*newton_vortdelv1(i,j)-& !u2*d(delv1)/dy & dx
            up(i,j,2)*(zeta_n(i,j)+fcor(i,j))-& !delv2*(w+f)
            ulatlon(i,j,2)*newton_vortdelv2(i,j)) !v2*d(delv2)/dx &dy



           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
             ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
             up(i,j,2)*dti +&                !delu2/dt
             ulatlon(i,j,1)*newton_vortdelv2(i,j)) !u1*d(delu2)/dy & dx


            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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
              if (n==1) fx(lx) = 0.0d0
              if (n==2) fx(lx) = 0.0d0
              
           else 
              if (n==1) fx(lx) = 0.0d0
              if (n==2) fx(lx) = 0.0d0

           end if
              if (n==3) fx(lx) =xs(lx)
            end do
          end do
       end do
     end do !ie
     end do



  end subroutine sw_picard_block_ID22


  subroutine test_id(xs, nelemd, fx, c_ptr_to_object) bind(C,name='test_id')
  use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only :  np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(in)        :: xs(nelemd)
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'test id'
!call flush(6)

!    call t_startf('test id')
!call flush(6)
!write(6,*)'startf'
!call flush(6)

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
!call flush(6)
!write(6,*)'cf pointer'
!call flush(6)

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

!write(6,*)'setting'
!call flush(6)
    fx=xs

!    call t_stopf('test id')

!write(6,*)'returning'
!call flush(6)
  end subroutine test_id


  subroutine simple_prec_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='simple_prec_op')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none

    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

    call t_startf('simple prec op')

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
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)
       end do
      end do

     do k=1,nlev

!       usum1       = 0.0d0
!       usum2       = 0.0d0
!       psum        = 0.0d0
!       tot_nv      = 1.0d0/float(nv*nv)

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   !
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   ! 
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   !
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   ! 

! set p=gh=phi-ghs=phi-pmean 
          !dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)-pmean
          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean 
! set dp
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)

!          usum1 = usum1 + ulatlon(i,j,1)
!          usum2 = usum2 + ulatlon(i,j,2)
!          psum = psum + fptr%base(ie)%state%p(i,j,k,n0)

       end do
      end do

!          usum1 = usum1*tot_nv
!          usum2 = usum2*tot_nv
!          psum =  psum*tot_nv
!      do j=1,np
!       do i=1,np
!          E_n(i,j) = 0.5D0*(usum1*ulatlon(i,j,1) + usum2*ulatlon(i,j,2))  + &
!                         fptr%base(ie)%state%p(i,j,k,n0) + ps(i,j)
!          E(i,j)   = 0.5D0*(usum1*up(i,j,1) + usum2*up(i,j,2))  + &
!                         fptr%base(ie)%state%p(i,j,k,np1) + ps(i,j)
! Aaron! figure out what E should be. u_old*u_new? or some form?
!          E_n(i,j) = 0.5D0*(ulatlon(i,j,1)**2 + ulatlon(i,j,2)**2)  + &
!                         fptr%base(ie)%state%p(i,j,k,n0) + ps(i,j)
!          E(i,j)   = 0.5D0*(up(i,j,1)**2 + up(i,j,2)**2)  + &
!                         fptr%base(ie)%state%p(i,j,k,np1) + ps(i,j)
! orig
!          pv_n(i,j,1) = ulatlon(i,j,1)*( pmean + fptr%base(ie)%state%p(i,j,k,n0) )
!          pv_n(i,j,2) = ulatlon(i,j,2)*( pmean + fptr%base(ie)%state%p(i,j,k,n0) )
!          pv(i,j,1) = up(i,j,1)*( pmean + fptr%base(ie)%state%p(i,j,k,np1) )
!          pv(i,j,2) = up(i,j,2)*( pmean + fptr%base(ie)%state%p(i,j,k,np1) )

! this is element (3,3) of precon matrix, "G" in the notes
          !pv_n(i,j,1) = usum1*( fptr%base(ie)%state%p(i,j,k,n0) )
          !pv_n(i,j,2) = usum2*( fptr%base(ie)%state%p(i,j,k,n0) )
          !pv(i,j,1) = usum1*( fptr%base(ie)%state%p(i,j,k,np1) )
          !pv(i,j,2) = usum2*( fptr%base(ie)%state%p(i,j,k,np1) )
! set dh
!          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)
!         dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)
!       end do
!      end do




! create operators for full residual calc
     !grade_n = gradient_sphere(E_n,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector
     grade_n1 = gradient_sphere(up(:,:,1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     grade_n2 = gradient_sphere(up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
!     div_n = divergence_sphere(pv_n,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 

!     grade = gradient_sphere(E,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector
!     zeta = vorticity_sphere(up,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
!     div = divergence_sphere(pv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 

          ! ==============================================
          ! Compute G-D(diagF)^-1 D^T(dh) 
          ! ==============================================

! compute terms in F_diag
!     gradef_n = gradient_sphere_diag(Ef_n,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector
!     gradef = gradient_sphere_diag(Ef,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector
!     zetaf_n = vorticity_sphere_diag(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
!     zetaf = vorticity_sphere_diag(up,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
! adjust these to be like v_s_diag but with only one component 
!     zetaf1 = vorticity_sphere_vx(up,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
!     zetaf2 = halfvorticity_sphere_uy(up,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 

! compute diagF that just has dt term
       do j=1,np
        do i=1,np
! zetaf here should be the one half
       diagf(i,j,1,k,ie) = spheremp(i,j)*dti 
       !+ rspheremp(i,j)*(ulatlon(i,j,2)*zetaf(i,j))
!       diagf(i,j,1,k,ie) = rspheremp(i,j)*dti + rspheremp(i,j)*(ulatlon(i,j,2)*zetaf1(i,j))
! zetaf here should be the other half
       diagf(i,j,2,k,ie) = spheremp(i,j)*dti 
       !+ rspheremp(i,j)*(-ulatlon(i,j,1)*zetaf(i,j))
!       diagf(i,j,2,k,ie) = rspheremp(i,j)*dti + rspheremp(i,j)*(-ulatlon(i,j,1)*zetaf2(i,j))
! compute (diagF)^-1
         diagf(i,j,1,k,ie)   =  1.0d0/(diagf(i,j,1,k,ie))
         diagf(i,j,2,k,ie)   =  1.0d0/(diagf(i,j,2,k,ie))

         !schur(i,j,1)   =  fptr%base(ie)%state%p(i,j,k,n0)*diagf(i,j,1,k,ie)*grade_n3(i,j,1)
         !schur(i,j,2)   =  fptr%base(ie)%state%p(i,j,k,n0)*diagf(i,j,2,k,ie)*grade_n3(i,j,2)

         schur(i,j,1)   =  diagf(i,j,1,k,ie)*grade_n3(i,j,1)
         schur(i,j,2)   =  diagf(i,j,2,k,ie)*grade_n3(i,j,2)

       end do 
      end do

! compute D (dh*diagF^-1*Dtranspose)
     div_h =  dh(i,j)*divergence_sphere(schur,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 


          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

            vtens(i,j,1,k,ie)= spheremp(i,j)*&
            (ulatlon(i,j,1)*grade_n1(i,j,1)+&
             ulatlon(i,j,2)*grade_n1(i,j,2)+&
             up(i,j,2)*(fcor(i,j) - zeta_n(i,j)) + grade_n3(i,j,1))

           vtens(i,j,2,k,ie)= spheremp(i,j)*&
            (ulatlon(i,j,1)*grade_n2(i,j,1)+&
             ulatlon(i,j,2)*grade_n2(i,j,2)+&
             up(i,j,1)*(fcor(i,j) + zeta_n(i,j)) + grade_n3(i,j,2))


           ! vtens(i,j,2,k,ie)=spheremp(i,j)* &
          !(-up(i,j,1)*(fcor(i,j) + zeta(i,j)) - grade(i,j,2))
          !  vtens_n(i,j,1,k,ie)=spheremp(i,j)* &
          !(ulatlon(i,j,2)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,1))
          !  vtens_n(i,j,2,k,ie)=spheremp(i,j)* &
          !(-ulatlon(i,j,1)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,2))

!            vtens(i,j,1,k,ie)=spheremp(i,j)* &
!          (usum2*(fcor(i,j) + zeta(i,j)) - grade(i,j,1))
!            vtens(i,j,2,k,ie)=spheremp(i,j)* &
!          (-usum1*(fcor(i,j) + zeta(i,j)) - grade(i,j,2))
!            vtens_n(i,j,1,k,ie)=spheremp(i,j)* &
!          (ulatlon(i,j,2)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,1))
!            vtens_n(i,j,2,k,ie)=spheremp(i,j)* &
!          (-ulatlon(i,j,1)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,2))

         ptens(i,j,k,ie) = spheremp(i,j)*&
            ( ulatlon(i,j,1)*grade_n3(i,j,1)+&
             ulatlon(i,j,2)*grade_n3(i,j,2) +dh_n(i,j)*div_h(i,j))
!  - ( approx schur complement + G )
      ! ptens(i,j,k,ie) =  -spheremp(i,j)*div_h(i,j) + ptens(i,j,k,ie)




!       vtens(i,j,1,k,ie) = spheremp(i,j)* &
!         (up(i,j,1)-ulatlon(i,j,1))*dti - vtens(i,j,1,k,ie)
       vtens(i,j,1,k,ie) = spheremp(i,j)*up(i,j,1)*dti + vtens(i,j,1,k,ie)

!       vtens(i,j,2,k,ie) = spheremp(i,j)* &
!         (up(i,j,2)-ulatlon(i,j,2))*dti - vtens(i,j,2,k,ie)
       vtens(i,j,2,k,ie) = spheremp(i,j)*up(i,j,2)*dti + vtens(i,j,2,k,ie)

!       ptens(i,j,k,ie)   = spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,np1)- &
!          fptr%base(ie)%state%p(i,j,k,n0))*dti - ptens(i,j,k,ie)
       ptens(i,j,k,ie)   = spheremp(i,j)*dh(i,j)*dti + ptens(i,j,k,ie)

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
               if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
               if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

!              if (n==1) fx(lx) = xs(lx)
!              if (n==2) fx(lx) = xs(lx)
           end if
              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
!              if (n==3) fx(lx) = xs(lx)
            end do
          end do
       end do
     end do !ie
     end do

    call t_stopf('simple prec op')

  end subroutine simple_prec_op




  subroutine test_picard_lin(xs, nelemd,fx, c_ptr_to_object) bind(C,name='test_picard_lin')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only :  np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none

    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(precon_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

  !  call t_startf('picard prec op')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

    write(6,*)'setting state=xs'
    call flush(6)
    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    write(6,*)'state=xs'
    call flush(6)
    do ie=ns,ne


      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
       end do
      end do

    write(6,*)'defining up and ulatlon '
    call flush(6)
     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   !
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   ! 
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   !
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   ! 

! set dh_n=p_notes=g(h+hmean)
          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean 

! set dp
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)


       end do
      end do

    write(6,*)'up and ulatlon defined'
    call flush(6)



! create operators for full residual calc
     grade_n1 = gradient_sphere(up(:,:,1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     grade_n2 = gradient_sphere(up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
    write(6,*)'grade terms defined'
    call flush(6)

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np
!
            vtens(i,j,1,k,ie)= spheremp(i,j)*&
            (ulatlon(i,j,1)*grade_n1(i,j,1)+&
             ulatlon(i,j,2)*grade_n1(i,j,2)+&
             up(i,j,2)*(fcor(i,j) - zeta_n(i,j)) + grade_n3(i,j,1))

           vtens(i,j,2,k,ie)= spheremp(i,j)*&
            (ulatlon(i,j,1)*grade_n2(i,j,1)+&
             ulatlon(i,j,2)*grade_n2(i,j,2)+&
             up(i,j,1)*(fcor(i,j) + zeta_n(i,j)) + grade_n3(i,j,2))


         ptens(i,j,k,ie) = spheremp(i,j)*&
            ( ulatlon(i,j,1)*grade_n3(i,j,1)+&
             ulatlon(i,j,2)*grade_n3(i,j,2)+ &
             dh_n(i,j)*(grade_n1(i,j,1) +grade_n2(i,j,2)))
!
!
!!!
!!!code is blowing up here for some reason
!       vtens(i,j,1,k,ie) = spheremp(i,j)*up(i,j,1)*dti + vtens(i,j,1,k,ie)
!!!


!
!       vtens(i,j,2,k,ie) = spheremp(i,j)*up(i,j,2)*dti + vtens(i,j,2,k,ie)
!
!       ptens(i,j,k,ie)   = spheremp(i,j)*dh(i,j)*dti + ptens(i,j,k,ie)
!
            end do
          end do
       end do


    write(6,*)'vtens and ptens defined'
    call flush(6)
!      ! ===================================================
!      ! Pack cube edges of tendencies, rotate velocities
!      ! ===================================================
!       kptr=0
!       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
!       kptr=nlev
!       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
   end do
!!
!!   !$OMP BARRIER
!!    call bndry_exchangeV(fptr%hybrid,fptr%edge3)
!!   !$OMP BARRIER
!!
    do ie=ns,ne
!
!       ! ===========================================================
!       ! Unpack the edges for vgradp and vtens
!       ! ===========================================================
!       kptr=0
!       call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, fptr%base(ie)%desc)
!       kptr=nlev
!       call edgeVunpack(fptr%edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, fptr%base(ie)%desc)
!
     end do !ie
       ! ===========================================================
       ! Compute velocity and pressure tendencies for all levels
       ! ===========================================================

!      lx = 0
!
!    do n=1,nvar
!    do ie=ns,ne
!      do k=1,nlev
!         do j=1,np
!            do i=1,np
!              lx = lx+1
!           if (topology == "cube" .and. test_case=="swtc1") then
!              if (n==1) fx(lx) = 0.0
!              if (n==2) fx(lx) = 0.0
!           else 
!                vtens1=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
!                vtens2=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
!
!              if (n==1) fx(lx) = fptr%base(ie)%Dinv(1,1,i,j)*vtens1 &
!                     + fptr%base(ie)%Dinv(1,2,i,j)*vtens2
!              if (n==2) fx(lx) = fptr%base(ie)%Dinv(2,1,i,j)*vtens1 &
!                     + fptr%base(ie)%Dinv(2,2,i,j)*vtens2
!           end if
!              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
!            end do
!          end do
!       end do
!     end do !ie
!     end do





!    call t_stopf('picard prec op')

  end subroutine test_picard_lin



  subroutine picard_prec_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='picard_lin')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none

    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

    call t_startf('picard prec op')

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
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne


      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
       end do
      end do

     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  

! set dh_n=p_notes=g(h+hmean)
          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean 

! set dp
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)


       end do
      end do




! create operators for full residual calc
     grade_n1 = gradient_sphere(up(:,:,1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     grade_n2 = gradient_sphere(up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

            vtens(i,j,1,k,ie)= spheremp(i,j)*&
            (ulatlon(i,j,1)*grade_n1(i,j,1)+&
             ulatlon(i,j,2)*grade_n1(i,j,2)+&
             up(i,j,2)*(fcor(i,j) - zeta_n(i,j)) + grade_n3(i,j,1))

           vtens(i,j,2,k,ie)= spheremp(i,j)*&
            (ulatlon(i,j,1)*grade_n2(i,j,1)+&
             ulatlon(i,j,2)*grade_n2(i,j,2)+&
             up(i,j,1)*(fcor(i,j) + zeta_n(i,j)) + grade_n3(i,j,2))


         ptens(i,j,k,ie) = spheremp(i,j)*&
            ( ulatlon(i,j,1)*grade_n3(i,j,1)+&
             ulatlon(i,j,2)*grade_n3(i,j,2)+ &
             dh_n(i,j)*(grade_n1(i,j,1) +grade_n2(i,j,2)))


! calculate nonlinear residual 



!       vtens(i,j,1,k,ie) = spheremp(i,j)* &
!         (up(i,j,1)-ulatlon(i,j,1))*dti - vtens(i,j,1,k,ie)
       vtens(i,j,1,k,ie) = spheremp(i,j)*up(i,j,1)*dti + vtens(i,j,1,k,ie)

!       vtens(i,j,2,k,ie) = spheremp(i,j)* &
!         (up(i,j,2)-ulatlon(i,j,2))*dti - vtens(i,j,2,k,ie)
       vtens(i,j,2,k,ie) = spheremp(i,j)*up(i,j,2)*dti + vtens(i,j,2,k,ie)

!       ptens(i,j,k,ie)   = spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,np1)- &
!          fptr%base(ie)%state%p(i,j,k,n0))*dti - ptens(i,j,k,ie)
       ptens(i,j,k,ie)   = spheremp(i,j)*dh(i,j)*dti + ptens(i,j,k,ie)

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

              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

!              if (n==1) fx(lx) = xs(lx)
!              if (n==2) fx(lx) = xs(lx)
           end if
              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
!              if (n==3) fx(lx) = xs(lx)
            end do
          end do
       end do
     end do !ie
     end do

    call t_stopf('picard prec op')


  end subroutine picard_prec_op



  subroutine picard_diag_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='picard_diag')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none

    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx


write(6,*)'picard diag\n'

    call t_startf('picard diag op')

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
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes

             !write(6,*)'xs=',xs(lx),''
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne


      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
       end do
      end do

     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   !  

! set dh_n=p_notes=g(h+hmean)
          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean 

! set dp
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)

           !  write(6,*)'ulatlon1=',ulatlon(i,j,1),''
           !  write(6,*)'ulatlon2=',ulatlon(i,j,2),''
           !  write(6,*)'pn=',dh_n(i,j),''

       end do
      end do

!stop


! create operators for full residual calc
     grade_n1 = gradient_sphere(up(:,:,1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     grade_n2 = gradient_sphere(up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector
     grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

          !  vtens(i,j,1,k,ie)= spheremp(i,j)*&
          !  (ulatlon(i,j,1)*grade_n1(i,j,1)+&
          !   ulatlon(i,j,2)*grade_n1(i,j,2)+&
          !   up(i,j,2)*(fcor(i,j) - zeta_n(i,j)) + grade_n3(i,j,1))
!
!           vtens(i,j,2,k,ie)= spheremp(i,j)*&
!            (ulatlon(i,j,1)*grade_n2(i,j,1)+&
!             ulatlon(i,j,2)*grade_n2(i,j,2)+&
!             up(i,j,1)*(fcor(i,j) + zeta_n(i,j)) + grade_n3(i,j,2))

        ! ptens(i,j,k,ie) = spheremp(i,j)*&
        !    ( ulatlon(i,j,1)*grade_n3(i,j,1)+&
        !     ulatlon(i,j,2)*grade_n3(i,j,2))+ &
        !     dh_n(i,j)*(grade_n1(i,j,1) +grade_n2(i,j,2))


              vtens(i,j,1,k,ie)= spheremp(i,j)*&
            (ulatlon(i,j,1)*grade_n1(i,j,1)+&
             ulatlon(i,j,2)*grade_n1(i,j,2))
             !up(i,j,2)*(fcor(i,j) - zeta_n(i,j)) + grade_n3(i,j,1))

           vtens(i,j,2,k,ie)= spheremp(i,j)*&
            (ulatlon(i,j,1)*grade_n2(i,j,1)+&
             ulatlon(i,j,2)*grade_n2(i,j,2))
             !up(i,j,1)*(fcor(i,j) + zeta_n(i,j)) + grade_n3(i,j,2))




         ptens(i,j,k,ie) = spheremp(i,j)*(&
             dh_n(i,j)*(grade_n1(i,j,1) +grade_n2(i,j,2)))
        !    ( ulatlon(i,j,1)*grade_n3(i,j,1)+&
        !     ulatlon(i,j,2)*grade_n3(i,j,2))+ &


! calculate nonlinear residual 

    ! Backward Euler 1st order method

!       vtens(i,j,1,k,ie) = spheremp(i,j)* &
!         (up(i,j,1)-ulatlon(i,j,1))*dti - vtens(i,j,1,k,ie)
       vtens(i,j,1,k,ie) = spheremp(i,j)*up(i,j,1)*dti  + vtens(i,j,1,k,ie)

!       vtens(i,j,2,k,ie) = spheremp(i,j)* &
!         (up(i,j,2)-ulatlon(i,j,2))*dti - vtens(i,j,2,k,ie)
       vtens(i,j,2,k,ie) = spheremp(i,j)*up(i,j,2)*dti  + vtens(i,j,2,k,ie)

!       ptens(i,j,k,ie)   = spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,np1)- &
!          fptr%base(ie)%state%p(i,j,k,n0))*dti - ptens(i,j,k,ie)
       ptens(i,j,k,ie)   = spheremp(i,j)*dh(i,j)*dti 
       !+ ptens(i,j,k,ie)

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

              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

           end if
              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
!              if (n==3) fx(lx) = xs(lx)
            end do
          end do
       end do
     end do !ie
     end do

    call t_stopf('picard diag op')


  end subroutine picard_diag_op




  subroutine sw_picard_diag_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_diag')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : EdgeBuffer_t, edgevpack, edgerotate, edgevunpack
    use hybrid_mod, only : hybrid_t
    use derivative_mod, only : derivative_t, gradient_sphere, &
          divergence_sphere, vorticity_sphere, divergence_sphere_wk, &
          vorticity_sphere_diag
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : precon_type,derived_type
    use perf_mod, only : t_startf, t_stopf
    use physical_constants, only : g

    implicit none
!Aaron fix this after testing should be intent in
    real (c_double)        :: xs(nelemd)
    !real (c_double) ,intent(in)        :: xs(nelemd)
!end fix
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    !type(precon_type) ,pointer        :: fptr=>NULL()
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
!    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf ! basic schur complement
    real (kind=real_kind), dimension(np,np,2,nlev,nelem) :: diagf_n ! basic schur complement
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: schur,schur_n 
    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: gradef,gradef_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,Ef,E_n,Ef_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: zetaf,zetaf_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um, utmpx, utmpy,uco,upco
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv
    real (kind=real_kind), dimension(np,np,nlev,nelem)  :: dp
!    real (kind=real_kind), dimension(np,np)     :: metdet

    real (kind=real_kind) ::  v1, v2, v1p, v2p, v1m, v2m
    real (kind=real_kind) ::  vtens1, vtens2
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, usum1, usum2, psum, tot_nv

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx,lenp, lx

!call flush(6)
!write(6,*)'sw startf'
!call flush(6)

!    call t_startf('sw jacobian op')


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    lenx       = fptr%n 
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

vtens=0.0d0
ptens=0.0d0


    lx = 0
     do n=1,nvar
      do ie=ns,ne
        do k=1,nlev
          do j=1,np
            do i=1,np
             lx = lx+1
             if (n==1) fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
             !if (n==1) base_vector(lx)=fptr%base(ie)%state%v(i,j,1,k,n0) 
             if (n==2) fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx)
             !if (n==2) base_vector(lx)=fptr%base(ie)%state%v(i,j,2,k,n0) 
             if (n==3) fptr%base(ie)%state%p(i,j,k,np1)  = xs(lx)!p =g*h in notes
             !if (n==3) base_vector(lx)=fptr%base(ie)%state%p(i,j,k,n0) 
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
     end do !nvar

    do ie=ns,ne



      do j=1,np
       do i=1,np
        fcor(i,j)        = fptr%base(ie)%fcor(i,j)
        mv(i,j)          = fptr%base(ie)%mp(i,j)
        spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
        rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
        ps(i,j)          = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
!        metdet(i,j)      = fptr%base(ie)%metdet(i,j)
       end do !np
      end do !np



     do k=1,nlev

         ! ==============================================
         ! Compute kinetic energy term at each time level
         ! ==============================================

      do j=1,np
       do i=1,np

! set u
          ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
          ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
! set du
          up(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
          up(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,np1)   ! 



!upco=delu_i=D'Dv=g_{ij}delv
!          upco(i,j,1)=fptr%base(ie)%D(1,1,i,j)*up(i,j,1) + fptr%base(ie)%D(2,1,i,j)*up(i,j,2) ! latlon->covariant
!          upco(i,j,2)=fptr%base(ie)%D(1,2,i,j)*up(i,j,1) + fptr%base(ie)%D(2,2,i,j)*up(i,j,2) ! latlon->covariant

          delv1(i,j,1)=up(i,j,1)
          delv1(i,j,2)=0.0d0
          delv2(i,j,1)=0.0d0
          delv2(i,j,2)=up(i,j,2)



! set dh_n=p_notes=g(h+hmean)=phi
          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)


          pv1(i,j,1) = up(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
          pv1(i,j,2) = 0.d00

          pv2(i,j,1) = 0.d00
          pv2(i,j,2) = up(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )

          dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
          dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )

       end do !np
      end do !np

! create operators for full residual calc
     grade_n1 = gradient_sphere(up(:,:,1),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_1
     grade_n2 = gradient_sphere(up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector (grad)delu_2
     grade_n3 = gradient_sphere(dh,fptr%deriv,fptr%base(ie)%Dinv)  ! scalar -> latlon vector        (grad)del phi

     !newton_grade_n1 = gradient_sphere(ulatlon(:,:,1)*up(:,:,1) ,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_1
     !newton_grade_n2 = gradient_sphere(ulatlon(:,:,2)*up(:,:,2),fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector (grad)u_2
     !newton_grade_n3 = gradient_sphere(dh_n,fptr%deriv,fptr%base(ie))  ! scalar -> latlon vector        (grad)phi

     zeta_n = vorticity_sphere(ulatlon,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv2 = vorticity_sphere(delv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w
     newton_vortdelv1 = vorticity_sphere(delv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar  w

     divpv1 = divergence_sphere(pv1,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     divpv2 = divergence_sphere(pv2,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 
     !divdpv = divergence_sphere(dpv,fptr%deriv,fptr%base(ie)) ! latlon vector -> scalar 



          ! ==============================================
          ! Compute residual terms
          ! ==============================================

      do j=1,np
        do i=1,np

! Jacobian


          vtens(i,j,1,k,ie)= spheremp(i,j)*( &
             up(i,j,1)*dti - &                 !delv1/dt
             ulatlon(i,j,2)*newton_vortdelv1(i,j)& !u2*d(delv1)/dy & dx
!!            newton_grade_n1(i,j,1)+ &      !delv1*d(v1)/dx
!            -up(i,j,2)*(zeta_n(i,j)+fcor(i,j))+& !delv2*(w+f)
!            -ulatlon(i,j,2)*newton_vortdelv2(i,j)+& !v2*d(delv2)/dx &dy
!!            newton_grade_n2(i,j,1)+ &      !-delv2*d(v2)/dx
!            grade_n3(i,j,1)) ! + d(delp)/dx
              )



           vtens(i,j,2,k,ie)= spheremp(i,j)*( &
!             up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+& !delu1*(w+f)
!             ulatlon(i,j,1)*newton_vortdelv1(i,j)+&  !u1*d(delu1)/dx & dy
!!             newton_grade_n2(i,j,2)+&  !delu1*d(u2)/dy
             up(i,j,2)*dti +&                !delu2/dt
             ulatlon(i,j,1)*newton_vortdelv2(i,j)& !u1*d(delu2)/dy & dx
!!            newton_grade_n1(i,j,2)+ &      !delu1*d(u1)/dy
!             grade_n3(i,j,2))!d(delp)/dy
             ) 


         ptens(i,j,k,ie) = spheremp(i,j)*(&
             dh(i,j)*dti)                        !delp/dt
             !dh(i,j)*dti+&                         !delp/dt
             !divpv1(i,j)+divpv2(i,j))

            end do !np
          end do !np
       end do !nlev


      ! ===================================================
      ! Pack cube edges of tendencies, rotate velocities
      ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,fptr%base(ie)%desc)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,fptr%base(ie)%desc)
  end do !ie


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

              if(n==1) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
              if(n==2) fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

           end if
              if (n==3) fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
            end do
          end do
       end do
     end do !ie
     end do


  end subroutine sw_picard_diag_op


end module implicit_mod

