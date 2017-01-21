#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef TRILINOS

!#define DEBUG_PRINT_ON   

module prim_jacobian_operator_mod

  use iso_c_binding, only : c_int

  use kinds, only                 : real_kind
  use dimensions_mod, only        : np, nlev, nelemd
  use edgetype_mod, only          : EdgeBuffer_t
  use edge_mod, only              : edgevpack, edgevunpack
  use bndry_mod, only             : bndry_exchangev
  use prim_derived_type_mod, only : derived_type
  use parallel_mod, only          : abortmp, iam
  use perf_mod, only              : t_startf, t_stopf
  
  implicit none
  private
  save

  ! initialize and finalize routines
  public :: prim_jacobian_operator_init
  public :: prim_jacobian_operator_finalize
  public :: prim_jacobian_operator_update

  ! exchange buffer
  type(EdgeBuffer_t) :: edge3p1

  real(kind=real_kind), allocatable, dimension(:,:,:,:,:) :: vel_n
  real(kind=real_kind), allocatable, dimension(:,:,:,:)   :: temp_n
  real(kind=real_kind), allocatable, dimension(:,:,:)     :: ps_n

  ! real(kind=real_kind), allocatable, dimension(:,:,:,:,:) :: grad_p_n

  ! real(kind=real_kind), allocatable, dimension(:,:,:,:) :: vort_n
  ! real(kind=real_kind), allocatable, dimension(:,:,:,:) :: pressi_n
  ! real(kind=real_kind), allocatable, dimension(:,:,:,:) :: pressi2_n
  ! real(kind=real_kind), allocatable, dimension(:,:,:,:) :: dp_n
   
contains 

  subroutine prim_jacobian_operator_init(par, elem)

    use parallel_mod, only   : parallel_t
    use element_mod, only    : element_t
    use edge_mod, only       : initEdgeBuffer
    use dimensions_mod, only : nlev
    use control_mod, only    : qsplit, rsplit

    implicit none

    type (parallel_t), intent(in) :: par
    type (element_t),  intent(in), target :: elem(:)

    if (rsplit==0) then
       call initEdgeBuffer(par, edge3p1, elem, 3*nlev+1)
    else
       ! need extra buffer space for dp3d
       call initEdgeBuffer(par, edge3p1, elem, 4*nlev+1)
    end if

    allocate(vel_n(np,np,2,nlev,nelemd))
    allocate(temp_n(np,np,nlev,nelemd))
    allocate(ps_n(np,np,nelemd))

    ! allocate(grad_p_n(np,np,2,nlev,nelemd))
    
    ! allocate(vort_n(np,np,nlev,nelemd))
    ! allocate(pressi_n(np,np,nlev,nelemd))
    ! allocate(pressi2_n(np,np,nlev,nelemd))
    ! allocate(dp_n(np,np,nlev,nelemd))

  end subroutine prim_jacobian_operator_init
  
  
  
  subroutine prim_jacobian_operator_finalize()

    use edge_mod, only : FreeEdgeBuffer
    
    implicit none

    call FreeEdgeBuffer(edge3p1)

    deallocate(vel_n, temp_n, ps_n)

    ! deallocate(grad_p_n)

    ! deallocate(vort_n)
    ! deallocate(pressi_n)
    ! deallocate(pressi2_n)
    ! deallocate(dp_n)
    
  end subroutine prim_jacobian_operator_finalize



  subroutine prim_jacobian_operator_update(xstate, nxstate, c_ptr_to_object)
    ! ---------------------------------------------------------------------------
    ! update current state with most recent iteration value and compute the blocks
    ! of the Jacobian matrix
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    ! use derivative_mod, only : vorticity_sphere, gradient_sphere, divergence_sphere

    implicit none

    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    type(c_ptr)                       :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    ! real(kind=real_kind), dimension(np,np,2,nelemd) :: grad_ps_n

    integer :: nets, nete
    integer :: i, j, k, ie, lx
    ! ---------------------------------------------------------------------------

#ifdef DEBUG_PRINT_ON   
    write(*,*) 'Start: prim_jacobian_operator_update'
#endif

    call t_startf('prim_jacobian_operator_update')

    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr
    
    nets = fptr%nets
    nete = fptr%nete

    ! initialize 1d array index
    lx = 1

    ! copy V1
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                vel_n(i,j,1,k,ie) = xstate(lx)
                lx = lx+1
             end do
          end do
       end do
    end do

    ! copy V2
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                vel_n(i,j,2,k,ie) = xstate(lx)
                lx = lx+1
             end do
          end do
       end do
    end do

    ! copy Ps
    do ie=nets,nete
       do j=1,np
          do i=1,np
             ps_n(i,j,ie) = xstate(lx)
             lx = lx+1
          end do
       end do
    end do

    ! copy T
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                temp_n(i,j,k,ie) = xstate(lx)
                lx = lx+1
             end do
          end do
       end do
    end do

    ! vorticity
    ! do ie = nets,nete
    !   do k = 1,nlev
    !      vort_n(:,:,k,ie) = vorticity_sphere(vel_n(:,:,:,k,ie), &
    !           fptr%deriv, fptr%base(ie))
    !   end do
    !end do

    ! pressure
    ! do ie = nets,nete
    !    do k = 1,nlev
    !       pressi_n(:,:,k,ie) = 1.0d0 / (fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 &
    !                          + fptr%hvcoord%hybm(k) * ps_n(:,:,ie))
    !    end do
    ! end do

    ! pressi2_n = pressi_n**2

    ! gradient of ps
    ! do ie = nets,nete
    !    grad_ps_n(:,:,:,ie) = gradient_sphere(ps_n(:,:,ie), &
    !         fptr%deriv, fptr%base(ie)%Dinv)
    ! end do

    ! do k = 1,nlev
    !    grad_p_n(:,:,:,k,ie) = fptr%hvcoord%hybm(k) * grad_ps_n(:,:,:,ie)
    ! end do

    ! dpressure / deta
    ! do ie = nets,nete
    !    do k = 1,nlev
    !       dp_n(:,:,k,ie) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
    !                      + fptr%hvcoord%hybi(k+1) * ps_n(:,:,ie)) &
    !                      - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
    !                      + fptr%hvcoord%hybi(k) * ps_n(:,:,ie))
    !    end do
    ! end do

    ! v_grad_p(:,:,k) = vel_n(:,:,1,k,ie) * grad_p(:,:,1) &
    !      + vel_n(:,:,2,k,ie) * grad_p(:,:,2)

    ! ! divergence of vertical derivative of pressure * vel                     
    ! do k = 1,nlev
    !    dp = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
    !         + fptr%hvcoord%hybi(k+1) * ps_n(:,:,ie)) &
    !         - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
    !         + fptr%hvcoord%hybi(k) * ps_n(:,:,ie))
       
    !    ! dpvel(:,:,1) = dp * fptr%base(ie)%state%v(:,:,1,k,n0)
    !    ! dpvel(:,:,2) = dp * fptr%base(ie)%state%v(:,:,2,k,n0)
    !    dpvel(:,:,1) = dp * vel_n(:,:,1,k,ie)
    !    dpvel(:,:,2) = dp * vel_n(:,:,2,k,ie)
       
    !    div_dpvel(:,:,k) = divergence_sphere(dpvel, fptr%deriv, fptr%base(ie))
    ! end do
    
    ! ! pressure vertical velocity
    ! call preq_omega_ps(omega_p, fptr%hvcoord, press, v_grad_p, div_dpvel)
    
    ! ! essentially everything needed for omega

    ! ! virtual temperature

    call t_stopf('prim_jacobian_operator_update')

#ifdef DEBUG_PRINT_ON   
    write(*,*) 'End: prim_jacobian_operator_update'
#endif
    
  end subroutine prim_jacobian_operator_update



  subroutine prim_Ablock_op(xstate, nxstate, fxstate, c_ptr_to_object) &
       bind(C,name='prim_Ablock_op')
    
    use, intrinsic :: iso_c_binding 

    use derivative_mod, only : vorticity_sphere, gradient_sphere
    use control_mod, only    : tstep_type

    implicit none

    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(out)       :: fxstate(nxstate)
    type(c_ptr)                       :: c_ptr_to_object

    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,2,nlev,nelemd) :: vtens

    real(kind=real_kind), dimension(np,np,2) :: dvel
    !real(kind=real_kind), dimension(np,np,2) :: grad_ke

    !real(kind=real_kind), dimension(np,np) :: veldvel
    real(kind=real_kind), dimension(np,np) :: vort, dvort
    real(kind=real_kind), dimension(np,np) :: fcor
    real(kind=real_kind), dimension(np,np) :: spheremp

    real(kind=real_kind) :: dti ! 1/dt

    integer :: nets, nete  ! element start and end index
    integer :: np1         ! time levels
    integer :: i, j, k, ie ! loop counters
    integer :: lx          ! 

#ifdef DEBUG_PRINT_ON
    write(*,*) "End: prim_Ablock_op"
#endif

    call t_startf('prim_Ablock_op')
    
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
    
    np1  = fptr%tl%np1  ! next Newton iter
    nets = fptr%nets    ! stating element index
    nete = fptr%nete    ! ending element index

    dti = 1.0d0/fptr%dt 
        
    if (tstep_type == 12) then
       dti = (3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2
    end if

    ! unpack velocities from input array
    lx = 1
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%v(i,j,1,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%v(i,j,2,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    ! compute A block matvec with spectral element operators
    do ie = nets,nete

       ! DSS weights
       spheremp = fptr%base(ie)%spheremp(:,:)

       ! coriolis force
       fcor = fptr%base(ie)%fcor(:,:)
       
       do k = 1,nlev
          
          ! delta velocity
          dvel = fptr%base(ie)%state%v(:,:,:,k,np1)

          ! vorticity
          vort  = vorticity_sphere(vel_n(:,:,:,k,ie), fptr%deriv, fptr%base(ie))
          dvort = vorticity_sphere(dvel, fptr%deriv, fptr%base(ie))

          ! kinetic energy
          ! veldvel = vel_n(:,:,1,k,ie) * dvel(:,:,1) &
          !         + vel_n(:,:,2,k,ie) * dvel(:,:,2)
          ! grad_ke = gradient_sphere(veldvel, fptr%deriv, fptr%base(ie)%Dinv)

          ! vertical advection

          ! compute Ablock * dvel 
          ! NOTE: shallow water drops KE contribution
          ! NOTE: dropped vertical advection contribution 
          vtens(:,:,1,k,ie) = spheremp(:,:) * &
               (dvel(:,:,1) * dti &
               - vel_n(:,:,2,k,ie) * dvort(:,:) &
               - dvel(:,:,2) * (vort(:,:) + fcor(:,:)))! &
               !+ grad_ke(:,:,1))
               
          vtens(:,:,2,k,ie) = spheremp(:,:) * &
               (dvel(:,:,2) * dti &
               + vel_n(:,:,1,k,ie) * dvort(:,:) &
               + dvel(:,:,1) * (vort(:,:) + fcor(:,:)))! &
               !+ grad_ke(:,:,2))
       end do
    end do

    ! pack edge data
    do ie = nets,nete
       call edgeVpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, 0, ie)
    end do

    ! exchange edge data
    call bndry_exchangeV(fptr%hybrid, edge3p1)

    ! unpack edge data
    do ie = nets,nete
       call edgeVunpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, 0, ie)
    end do

    ! pack result into output array
    lx = 1
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                lx = lx + 1
             end do
          end do
       end do
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
                lx = lx + 1
             end do
          end do
       end do
    end do
    
    call t_stopf('prim_Ablock_op')

#ifdef DEBUG_PRINT_ON
    write(*,*) "End: prim_Ablock_op"
#endif

  end subroutine prim_Ablock_op



  subroutine prim_Bblock_op(xstate, nxstate, fxstate, nfxstate, c_ptr_to_object) &
       bind(C,name='prim_Bblock_op')
    
    use, intrinsic :: iso_c_binding 

    use derivative_mod, only     : gradient_sphere
    use physical_constants, only : Rgas
    use physics_mod, only        : virtual_temperature

    implicit none

    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(out)       :: fxstate(nfxstate)
    integer(c_int), intent(in), value :: nfxstate
    type(c_ptr)                       :: c_ptr_to_object

    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,2,nlev,nelemd) :: vtens

    real(kind=real_kind), dimension(np,np,nlev) :: Tv
    !real(kind=real_kind), dimension(np,np,nlev) :: dp, Qt

    real(kind=real_kind), dimension(np,np,2) :: grad_ps, grad_dps

    real(kind=real_kind), dimension(np,np) :: dps
    real(kind=real_kind), dimension(np,np) :: press, pressi, pressi2    
    real(kind=real_kind), dimension(np,np) :: spheremp

    real(kind=real_kind) :: beta

    integer :: nets, nete  ! element start and end index
    integer :: np1         ! time levels
    integer :: qn0         ! moisture
    integer :: i, j, k, ie ! loop counters
    integer :: lx          ! 

#ifdef DEBUG_PRINT_ON
    write(*,*) "Start: prim_Bblock_op"
#endif

    call t_startf('prim_Bblock_op')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    np1  = fptr%tl%np1
    nets = fptr%nets 
    nete = fptr%nete 
    qn0  = fptr%n_Q
    
    ! unpack surface pressure from input array
    lx = 1
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             fptr%base(ie)%state%ps_v(i,j,np1) = xstate(lx)
             lx = lx + 1
          end do
       end do
    end do

    ! compute B block matvec with spectral element operators
    do ie = nets,nete

       ! DSS weights
       spheremp = fptr%base(ie)%spheremp(:,:)

       ! compute T_v for timelevel n0
       if (qn0 == -1 ) then
          Tv(:,:,:) = temp_n(:,:,:,ie)
       else
          call abortmp('ERROR: in prim_Bblock_op, moist run')
          ! need to add scaling factor for derivative of Tv
          ! do k = 1,nlev
          !    dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
          !         + fptr%hvcoord%hybi(k+1) * ps_n(:,:,ie)) &
          !         - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
          !         + fptr%hvcoord%hybi(k) * ps_n(:,:,ie))
          ! end do

          ! Qt = fptr%base(ie)%state%Qdp(:,:,:,1,qn0)/dp(:,:,:)
          ! Tv = virtual_temperature(temp_n(:,:,:,ie), Qt)
       end if

       ! delta surface pressure
       dps = fptr%base(ie)%state%ps_v(:,:,np1)

       ! gradient of delta surface pressure
       grad_dps = gradient_sphere(dps, fptr%deriv, fptr%base(ie)%Dinv)
       
       grad_ps = gradient_sphere(ps_n(:,:,ie), fptr%deriv, fptr%base(ie)%Dinv)

       do k = 1,nlev

          ! surface pressure and pressure
          press = fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybm(k) * ps_n(:,:,ie)
          
          beta = fptr%hvcoord%hybm(k)

          pressi  = 1.0d0 / press
          pressi2 = 1.0d0 / press**2

          ! geopotential 
          ! call preq_hydrostatic(fptr%base(ie)%derived%phi, & 
          !      fptr%base(ie)%state%phis,T_v,p,dp)
          
          vtens(:,:,1,k,ie) = spheremp * Rgas * Tv(:,:,k) * beta * &
               (-1.0d0 * pressi2 * beta * grad_ps(:,:,1) * dps &
               + pressi * grad_dps(:,:,1))
          
          vtens(:,:,2,k,ie) = spheremp * Rgas * Tv(:,:,k) * beta * &
               (-1.0d0 * pressi2 * beta * grad_ps(:,:,2) * dps &
               + pressi * grad_dps(:,:,2))
       end do
    end do

    ! pack edge data
    do ie = nets,nete
       call edgeVpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, 0, ie)
    end do

    ! exchange edge data
    call bndry_exchangeV(fptr%hybrid, edge3p1)

    ! unpack edge data
    do ie = nets,nete
       call edgeVunpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, 0, ie)
    end do

    ! pack result into output array
    lx = 1
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                lx = lx + 1
             end do
          end do
       end do
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
                lx = lx + 1
             end do
          end do
       end do
    end do
    
    call t_stopf('prim_Bblock_op')

#ifdef DEBUG_PRINT_ON
    write(*,*) "End: prim_Bblock_op"
#endif

  end subroutine prim_Bblock_op
  
  
  
  subroutine prim_Cblock_op(xstate, nxstate, fxstate, nfxstate, c_ptr_to_object) &
       bind(C,name='prim_Cblock_op')
    
    use, intrinsic :: iso_c_binding 

    use derivative_mod, only     : gradient_sphere
    use physical_constants, only : Rgas
    use physics_mod, only        : virtual_temperature

    implicit none

    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(out)       :: fxstate(nfxstate)
    integer(c_int), intent(in), value :: nfxstate
    type(c_ptr)                       :: c_ptr_to_object

    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,2,nlev,nelemd) :: vtens

    !real(kind=real_kind), dimension(np,np,nlev) :: dp
    real(kind=real_kind), dimension(np,np,nlev) :: Qt
    real(kind=real_kind), dimension(np,np,nlev) :: Ttemp, Tvfactor

    real(kind=real_kind), dimension(np,np,2) :: grad_ps

    real(kind=real_kind), dimension(np,np) :: press, pressi
    real(kind=real_kind), dimension(np,np) :: spheremp

    real(kind=real_kind) :: beta

    integer :: nets, nete  ! element start and end index
    integer :: np1         ! time levels
    integer :: qn0         ! moisture
    integer :: i, j, k, ie ! loop counters
    integer :: lx          ! 

#ifdef DEBUG_PRINT_ON
    write(*,*) "Start: prim_Cblock_op"
#endif

    call t_startf('prim_Cblock_op')
    
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
    
    np1  = fptr%tl%np1
    nets = fptr%nets 
    nete = fptr%nete 
    qn0  = fptr%n_Q
   
    ! unpack temperature from input array
    lx = 1
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%T(i,j,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    Ttemp    = 1.0d0          
    Tvfactor = 1.0d0
    
    do ie = nets,nete

       ! DSS weights
       spheremp = fptr%base(ie)%spheremp(:,:)

       ! compute T_v for timelevel n0
       if ( qn0 /= -1 ) then
          call abortmp('ERROR: in prim_Cblock_op, moist run')
          ! double check scaling factor
          ! do k = 1,nlev
          !    dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
          !         + fptr%hvcoord%hybi(k+1) * ps_n(:,:,ie)) &
          !         - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
          !         + fptr%hvcoord%hybi(k) * ps_n(:,:,ie))
          ! end do

          ! Qt       = fptr%base(ie)%state%Qdp(:,:,:,1,qn0)/dp(:,:,:)
          ! Tvfactor = virtual_temperature(Ttemp, Qt)
       end if

       grad_ps = gradient_sphere(ps_n(:,:,ie), fptr%deriv, fptr%base(ie)%Dinv)

       do k = 1,nlev

          press = fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybm(k) * ps_n(:,:,ie)
          
          pressi = 1.0d0 / press

          beta = fptr%hvcoord%hybm(k)
          
          ! geopotential 
          !call preq_hydrostatic(fptr%base(ie)%derived%phi,fptr%base(ie)%state%phis,T_v,p,dp)
          
          vtens(:,:,1,k,ie) = spheremp * &
               (Rgas * pressi * beta * grad_ps(:,:,1) * &
               Tvfactor(:,:,k) * fptr%base(ie)%state%T(:,:,k,np1))
          
          vtens(:,:,2,k,ie) = spheremp * &
               (Rgas * pressi * beta * grad_ps(:,:,2) * &
               Tvfactor(:,:,k) * fptr%base(ie)%state%T(:,:,k,np1))
       end do
    end do

    ! pack edge data
    do ie = nets,nete
       call edgeVpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, 0, ie)
    end do

    ! exchange edge data
    call bndry_exchangeV(fptr%hybrid, edge3p1)

    ! unpack edge data
    do ie = nets,nete
       call edgeVunpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, 0, ie)
    end do

    ! pack result into output array
    lx = 1
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                lx = lx + 1
             end do
          end do
       end do
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
                lx = lx + 1
             end do
          end do
       end do
    end do
    
    call t_stopf('prim_Cblock_op')

#ifdef DEBUG_PRINT_ON
    write(*,*) "End: prim_Cblock_op"
#endif
    
  end subroutine prim_Cblock_op



  subroutine prim_Dblock_op(xstate, nxstate, fxstate, nfxstate, c_ptr_to_object) &
       bind(C,name='prim_Dblock_op')

    use, intrinsic :: iso_c_binding 

    use derivative_mod, only : divergence_sphere

    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(out)       :: fxstate(nfxstate)
    integer(c_int), intent(in), value :: nfxstate
    type(c_ptr)                       :: c_ptr_to_object

    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,nelemd) :: pstens  

    !real(kind=real_kind), dimension(np,np,2) :: dvel
    real(kind=real_kind), dimension(np,np,2) :: dpdvel

    real(kind=real_kind), dimension(np,np) :: div_dpdvel
    real(kind=real_kind), dimension(np,np) :: dp
    real(kind=real_kind), dimension(np,np) :: spheremp

    integer :: nets, nete  ! element start and end index
    integer :: np1         ! time levels
    integer :: i, j, k, ie ! loop counters
    integer :: lx          ! 

#ifdef DEBUG_PRINT_ON
    write(*,*) "Start: prim_Dblock_op"
#endif

    call t_startf('prim_Dblock_op')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    np1  = fptr%tl%np1
    nets = fptr%nets 
    nete = fptr%nete 

    ! unpack velocity from input array
    lx = 1
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np 
             do i = 1,np
                fptr%base(ie)%state%v(i,j,1,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np 
             do i = 1,np
                fptr%base(ie)%state%v(i,j,2,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    ! initialize matvec
    pstens = 0.0d0

    ! compute matvec
    do ie = nets,nete

       ! DSS weights
       spheremp = fptr%base(ie)%spheremp(:,:)

       do k = 1,nlev

          dp = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybi(k+1) * ps_n(:,:,ie)) &
               - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybi(k) * ps_n(:,:,ie))

          dpdvel(:,:,1) = dp * fptr%base(ie)%state%v(:,:,1,k,np1)
          dpdvel(:,:,2) = dp * fptr%base(ie)%state%v(:,:,2,k,np1)
          
          div_dpdvel = divergence_sphere(dpdvel, fptr%deriv, fptr%base(ie))

          pstens(:,:,ie) = pstens(:,:,ie) + div_dpdvel
       end do

       pstens(:,:,ie) = spheremp * pstens(:,:,ie)
       
    end do

    ! pack edge data
    do ie = nets,nete
       call edgeVpack(edge3p1, pstens(:,:,ie), 1, 0, ie)
    end do

    ! exchange edge data
    call bndry_exchangeV(fptr%hybrid, edge3p1)

    ! unpack edge data
    do ie = nets,nete
       call edgeVunpack(edge3p1, pstens(:,:,ie), 1, 0, ie)
    end do

    ! pack result into output array
    lx = 1
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             fxstate(lx) = fptr%base(ie)%rspheremp(i,j)*pstens(i,j,ie)
             lx = lx + 1
          end do
       end do
    end do
 
    call t_stopf('prim_Dblock_op')

#ifdef DEBUG_PRINT_ON
    write(*,*) "End: prim_Dblock_op"
#endif

  end subroutine prim_Dblock_op



  subroutine prim_Eblock_op(xstate, nxstate, fxstate, c_ptr_to_object) &
       bind(C,name='prim_Eblock_op')
    
    use, intrinsic :: iso_c_binding 
    
    use control_mod, only    : tstep_type
    use derivative_mod, only : divergence_sphere
    
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(out)       :: fxstate(nxstate)
    type(c_ptr)                       :: c_ptr_to_object
    
    type(derived_type), pointer :: fptr=>NULL()
    
    real(kind=real_kind), dimension(np,np,nelemd) :: pstens

    !real(kind=real_kind), dimension(np,np,2) :: vel
    real(kind=real_kind), dimension(np,np,2) :: dpsvel

    real(kind=real_kind), dimension(np,np) :: div_dpsvel
    real(kind=real_kind), dimension(np,np) :: dps
    real(kind=real_kind), dimension(np,np) :: spheremp

    real(kind=real_kind) :: dti

    integer :: nets, nete  ! element start and end index
    integer :: np1         ! time levels
    integer :: i, j, k, ie ! loop counters
    integer :: lx          ! 

#ifdef DEBUG_PRINT_ON
    write(*,*) "Start: prim_Eblock_op"
#endif

    call t_startf('prim_Eblock_op')
    
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
    
    np1  = fptr%tl%np1
    nets = fptr%nets 
    nete = fptr%nete 

    dti = 1.0d0/fptr%dt
        
    if (tstep_type == 12) then
       dti = (3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2
    end if

    ! unpack surface pressure from input array
    lx = 1
    do ie = nets,nete
       do j = 1,np 
          do i = 1,np
             fptr%base(ie)%state%ps_v(i,j,np1) = xstate(lx)
             lx = lx + 1
          end do
       end do
    end do

    do ie = nets,nete

       ! DSS weights
       spheremp = fptr%base(ie)%spheremp(:,:)

       dps = fptr%base(ie)%state%ps_v(:,:,np1)

       ! initialize with time contribution
       pstens(:,:,ie) = dti * dps
       
       do k = 1,nlev
          
          dpsvel(:,:,1) = (fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)) * &
               dps * vel_n(:,:,1,k,ie)

          dpsvel(:,:,2) = (fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)) * &
               dps * vel_n(:,:,2,k,ie)
          
          div_dpsvel = divergence_sphere(dpsvel, fptr%deriv, fptr%base(ie))

          pstens(:,:,ie) = pstens(:,:,ie) + div_dpsvel
       end do
       
       pstens(:,:,ie) = spheremp * pstens(:,:,ie)
    end do

    ! pack edge data
    do ie = nets,nete
       call edgeVpack(edge3p1, pstens(:,:,ie), 1, 0, ie)
    end do

    ! exchange edge data
    call bndry_exchangeV(fptr%hybrid, edge3p1)

    ! unpack edge data
    do ie = nets,nete
       call edgeVunpack(edge3p1, pstens(:,:,ie), 1, 0, ie)
    end do

    ! pack result into output array
    lx = 1
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             fxstate(lx) = fptr%base(ie)%rspheremp(i,j)*pstens(i,j,ie)
             lx = lx + 1
          end do
       end do
    end do
 
    call t_stopf('prim_Eblock_op')

#ifdef DEBUG_PRINT_ON
    write(*,*) "End: prim_Eblock_op"
#endif

  end subroutine prim_Eblock_op

  ! Add Fblock

  ! Add Gblock

  subroutine prim_Hblock_op(xstate, nxstate, fxstate, c_ptr_to_object) &
       bind(C,name='prim_Hblock_op')
    
    use, intrinsic :: iso_c_binding 
    
    use control_mod, only        : tstep_type   
    use derivative_mod, only     : gradient_sphere, divergence_sphere
    use prim_si_mod, only        : preq_omega_ps
    use physical_constants, only : kappa
    
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(out)       :: fxstate(nxstate)
    type(c_ptr)                       :: c_ptr_to_object
    
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,nlev,nelemd) :: ttens

    real(kind=real_kind), dimension(np,np,nlev) :: v_grad_p
    real(kind=real_kind), dimension(np,np,nlev) :: press
    real(kind=real_kind), dimension(np,np,nlev) :: omega_p
    real(kind=real_kind), dimension(np,np,nlev) :: div_dpvel
    real(kind=real_kind), dimension(np,np,nlev) :: kappa_star

    real(kind=real_kind), dimension(np,np,2) :: grad_ps, grad_p
    real(kind=real_kind), dimension(np,np,2) :: dpvel
    !real(kind=real_kind), dimension(np,np,2) :: vel
    real(kind=real_kind), dimension(np,np,2) :: grad_dtemp

    real(kind=real_kind), dimension(np,np) :: dp
    real(kind=real_kind), dimension(np,np) :: dtemp
    real(kind=real_kind), dimension(np,np) :: v_grad_dtemp
    real(kind=real_kind), dimension(np,np) :: spheremp

    real(kind=real_kind) :: dti ! 1/dt

    integer :: nets, nete  ! element start and end index
    integer :: np1         ! time levels
    integer :: qn0         ! moisture
    integer :: i, j, k, ie ! loop counters
    integer :: lx          ! 

#ifdef DEBUG_PRINT_ON
    write(*,*) "Start: prim_Hblock_op"
#endif
    
    call t_startf('prim_Hblock_op')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    np1  = fptr%tl%np1
    nets = fptr%nets 
    nete = fptr%nete 
    qn0  = fptr%n_Q

    dti = 1.0d0/fptr%dt 

    if (tstep_type == 12) then
       dti = (3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2
    end if

    ! unpack temperature from input array
    lx = 1
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np 
             do i = 1,np
                fptr%base(ie)%state%T(i,j,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    if ( qn0 == -1 ) then
       kappa_star = kappa
    else
       call abortmp('ERROR: in prim_Hblock_op, moist run')
       ! Qt = fptr%base(ie)%state%Qdp(i,j,k,1,qn0)/dp(i,j,k)
       ! if (use_cpstar==1) then
       !    kappa_star(i,j,k) =  Rgas/Virtual_Specific_Heat(Qt)
       ! else
       !    kappa_star(i,j,k) = kappa
       ! end if
    end if

    
    do ie = nets,nete
       
       ! DSS weights
       spheremp = fptr%base(ie)%spheremp(:,:)

       ! pressure
       do k = 1,nlev
          press(:,:,k) = fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybm(k) * ps_n(:,:,ie)
       end do

       ! vel dot gradient of pressure
       grad_ps = gradient_sphere(ps_n(:,:,ie), fptr%deriv, fptr%base(ie)%Dinv)

       do k = 1,nlev
          grad_p = fptr%hvcoord%hybm(k) * grad_ps

          v_grad_p(:,:,k) = vel_n(:,:,1,k,ie) * grad_p(:,:,1) &
               + vel_n(:,:,2,k,ie) * grad_p(:,:,2)
       end do

       ! divergence of vertical derivative of pressure * vel                     
       do k = 1,nlev
          dp = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybi(k+1) * ps_n(:,:,ie)) &
               - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybi(k) * ps_n(:,:,ie))
          
          dpvel(:,:,1) = dp * vel_n(:,:,1,k,ie)
          dpvel(:,:,2) = dp * vel_n(:,:,2,k,ie)

          div_dpvel(:,:,k) = divergence_sphere(dpvel, fptr%deriv, fptr%base(ie))
       end do

       ! pressure vertical velocity
       call preq_omega_ps(omega_p, fptr%hvcoord, press, v_grad_p, div_dpvel)
              
       do k = 1,nlev
          
          dtemp = fptr%base(ie)%state%T(:,:,k,np1)
          
          grad_dtemp = gradient_sphere(dtemp, fptr%deriv, fptr%base(ie)%Dinv)

          ! vel dot gradient of temperature
          v_grad_dtemp = vel_n(:,:,1,k,ie) * grad_dtemp(:,:,1) &
               + vel_n(:,:,2,k,ie) * grad_dtemp(:,:,2)

          ! vertical advection of temperature

          ttens(:,:,k,ie) = spheremp * &
               (dti * dtemp &
               + v_grad_dtemp &
               - kappa_star(:,:,k) * dtemp(:,:) * omega_p(:,:,k))

       end do
    end do

    ! pack edge data
    do ie = nets,nete
       call edgeVpack(edge3p1, ttens(:,:,:,ie), nlev, 0, ie)
    end do

    ! exchange edge data
    call bndry_exchangeV(fptr%hybrid, edge3p1)
    
    ! unpack edge data
    do ie = nets,nete
       call edgeVunpack(edge3p1, ttens(:,:,:,ie), nlev, 0, ie)
    end do
    
    ! pack result into output array
    lx = 1
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = fptr%base(ie)%rspheremp(i,j)*ttens(i,j,k,ie)
                lx = lx + 1
             end do
          end do
       end do
    end do

    call t_stopf('prim_Hblock_op')

#ifdef DEBUG_PRINT_ON
    write(*,*) "End: prim_Hblock_op"
#endif
    
  end subroutine prim_Hblock_op



  ! Add AhatInv * B operator
  subroutine prim_AhatInvB_op(xstate, nxstate, fxstate, nfxstate, c_ptr_to_object) &
       bind(C,name="prim_AhatInvB_op")

    use, intrinsic :: iso_c_binding 
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nelemd
    use edge_mod, only              : edgevpack, edgevunpack
    use bndry_mod, only             : bndry_exchangev
    use prim_derived_type_mod, only : derived_type
    
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(out)       :: fxstate(nfxstate)
    integer(c_int), intent(in), value :: nfxstate
    type(c_ptr)                       :: c_ptr_to_object
    
    type(derived_type), pointer :: fptr=>NULL()

#ifdef DEBUG_PRINT_ON
    write(*,*) "Start: prim_AhatInvB_op"
#endif
    
    call t_startf('prim_AhatInvB_op')
    
    call prim_Bblock_op(xstate, nxstate, fxstate, nfxstate, c_ptr_to_object)

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    fxstate = fptr%dt * fxstate

    call t_stopf('prim_AhatInvB_op')

#ifdef DEBUG_PRINT_ON
    write(*,*) "End: prim_AhatInvB_op"
#endif

  end subroutine prim_AhatInvB_op



  ! Add D * AhatInv operator
  subroutine prim_DAhatInv_op(xstate, nxstate, fxstate, nfxstate, c_ptr_to_object) &
       bind(C,name="prim_DAhatInv_op")
    
    use, intrinsic :: iso_c_binding 
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nelemd
    use edge_mod, only              : edgevpack, edgevunpack
    use bndry_mod, only             : bndry_exchangev
    use prim_derived_type_mod, only : derived_type
    
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(out)       :: fxstate(nfxstate)
    integer(c_int), intent(in), value :: nfxstate
    type(c_ptr)                       :: c_ptr_to_object
    
    type(derived_type), pointer :: fptr=>NULL()

    real(c_double) :: xstate_temp(nxstate)

    
#ifdef DEBUG_PRINT_ON
    write(*,*) "Start: prim_DAhatInv_op"
#endif
    
    call t_startf('prim_DAhatInv_op')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    xstate_temp = fptr%dt * xstate
    
    call prim_Dblock_op(xstate_temp, nxstate, fxstate, nfxstate, c_ptr_to_object)

    call t_stopf('prim_DAhatInv_op')

#ifdef DEBUG_PRINT_ON
    write(*,*) "End: prim_DAhatInv_op"
#endif

  end subroutine prim_DAhatInv_op



  ! Add Schur Complement operator S = E - DAB
  subroutine prim_SchurS_op(xstate, nxstate, fxstate, c_ptr_to_object) &
       bind(C,name="prim_SchurS_op")

    use, intrinsic :: iso_c_binding 
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nelemd
    use edge_mod, only              : edgevpack, edgevunpack
    use bndry_mod, only             : bndry_exchangev
    use prim_derived_type_mod, only : derived_type
        
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(out)       :: fxstate(nxstate)
    type(c_ptr)                       :: c_ptr_to_object
    
    type(derived_type), pointer :: fptr=>NULL()
    
    real(c_double) :: ftemp1(2*np*np*nlev*nelemd)
    real(c_double) :: ftemp2(np*np*nelemd)
    real(c_double) :: ftemp3(np*np*nelemd)

    real(kind=real_kind), dimension(np,np,nelemd) :: pstens
    
    integer :: nets, nete
    integer :: i, j, ie
    integer :: lx

#ifdef DEBUG_PRINT_ON
    write(*,*) "Start: prim_SchurS_op"
#endif

    call prim_AhatInvB_op(xstate, nxstate, ftemp1, 2*np*np*nlev*nelemd, &
         c_ptr_to_object)

    call prim_Dblock_op(ftemp1, 2*np*np*nlev*nelemd, ftemp2, np*np*nelemd, &
         c_ptr_to_object)

    call prim_Eblock_op(xstate, nxstate, ftemp3, c_ptr_to_object)

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    nets = fptr%nets 
    nete = fptr%nete 

    ! combine
    lx = 1
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             pstens(i,j,ie) = fptr%base(ie)%spheremp(i,j) * &
                  (ftemp3(lx) - ftemp2(lx))
             lx = lx + 1
          end do
       end do
    end do  

    ! pack edge data
    do ie = nets,nete
       call edgeVpack(edge3p1, pstens(:,:,ie), 1, 0, ie)
    end do
    
    ! exchange edge data
    call bndry_exchangeV(fptr%hybrid, edge3p1)
    
    ! unpack edge data
    do ie = nets,nete
       call edgeVunpack(edge3p1, pstens(:,:,ie), 1, 0, ie)
    end do
    
    ! output  
    lx = 1
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             fxstate(lx) = fptr%base(ie)%rspheremp(i,j) * pstens(i,j,ie)
             lx = lx + 1
          end do
       end do
    end do

#ifdef DEBUG_PRINT_ON
    write(*,*) "End: prim_SchurS_op"
#endif

  end subroutine prim_SchurS_op

  ! approximation of full Jacobian operator
  subroutine prim_approx_jacobian_op(xstate, nxstate, fxstate, c_ptr_to_object) &
       bind(C,name='prim_approx_jacobian_op')

    use, intrinsic :: iso_c_binding 
    
    use control_mod, only        : tstep_type
    use physical_constants, only : Rgas, kappa
    use prim_si_mod, only        : preq_omega_ps
    use derivative_mod, only     : vorticity_sphere, gradient_sphere, &
                                   divergence_sphere

    implicit none
    
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(out)       :: fxstate(nxstate)
    type(c_ptr)                       :: c_ptr_to_object
    
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,2,nlev,nelemd) :: vtens
    real(kind=real_kind), dimension(np,np,nelemd)        :: pstens
    real(kind=real_kind), dimension(np,np,nlev,nelemd)   :: ttens

    real(kind=real_kind), dimension(np,np,2) :: dvel
    real(kind=real_kind), dimension(np,np)   :: dps
    real(kind=real_kind), dimension(np,np)   :: dtemp

    real(kind=real_kind), dimension(np,np,nlev) :: Tv
    real(kind=real_kind), dimension(np,np,nlev) :: kappa_star
    real(kind=real_kind), dimension(np,np,nlev) :: press
    real(kind=real_kind), dimension(np,np,nlev) :: v_grad_p
    real(kind=real_kind), dimension(np,np,nlev) :: dp
    real(kind=real_kind), dimension(np,np,nlev) :: omega_p
    real(kind=real_kind), dimension(np,np,nlev) :: div_dpvel

    real(kind=real_kind), dimension(np,np,2) :: grad_dps
    real(kind=real_kind), dimension(np,np,2) :: grad_ps
    real(kind=real_kind), dimension(np,np,2) :: grad_p
    real(kind=real_kind), dimension(np,np,2) :: grad_dtemp
    real(kind=real_kind), dimension(np,np,2) :: dpvel
    real(kind=real_kind), dimension(np,np,2) :: dpdvel
    real(kind=real_kind), dimension(np,np,2) :: dpsvel
    real(kind=real_kind), dimension(np,np,2) :: grad_ke

    real(kind=real_kind), dimension(np,np) :: fcor
    real(kind=real_kind), dimension(np,np) :: vort, dvort
    real(kind=real_kind), dimension(np,np) :: pressi, pressi2
    real(kind=real_kind), dimension(np,np) :: div_dpdvel
    real(kind=real_kind), dimension(np,np) :: div_dpsvel
    real(kind=real_kind), dimension(np,np) :: v_grad_dtemp
    real(kind=real_kind), dimension(np,np) :: veldvel
    real(kind=real_kind), dimension(np,np) :: spheremp

    real(kind=real_kind), dimension(nlev) :: beta

    real(kind=real_kind) :: dti

    integer :: nets, nete  ! element start and end index
    integer :: np1         ! time levels
    integer :: i, j, k, ie ! loop counters
    integer :: lx          ! 
    integer :: qn0         ! moisture

#ifdef DEBUG_PRINT_ON
    write(*,*) 'Start: prim_approx_jacobian_op'
#endif

    call t_startf('prim_approx_jacobian_op')
    
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    np1  = fptr%tl%np1  ! Newton update
    nets = fptr%nets    ! stating element index
    nete = fptr%nete    ! ending element index
    qn0  = fptr%n_Q     ! moist run flag

    dti = 1.0d0/fptr%dt 

    if (tstep_type == 12) then
       dti = (3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2
    end if

    ! initialize 1d array index
    lx = 1

    ! copy V1
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%v(i,j,1,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    ! copy V2
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%v(i,j,2,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    ! copy Ps
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             fptr%base(ie)%state%ps_v(i,j,np1) = xstate(lx)
             lx = lx + 1
          end do
       end do
    end do

    ! copy V2
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%T(i,j,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    ! Ttemp    = 1.0d0          
    ! Tvfactor = 1.0d0

    do k = 1,nlev
       beta(k) = fptr%hvcoord%hybm(k)
    end do

    ! compute Jacobian matvec with spectral element operators
    do ie = nets,nete

       ! DSS weights
       spheremp = fptr%base(ie)%spheremp(:,:)
       
       ! coriolis force
       fcor = fptr%base(ie)%fcor(:,:)
       
       ! compute T_v for timelevel n0
       if ( qn0 == -1 ) then
          Tv(:,:,:)  = temp_n(:,:,:,ie)
          kappa_star = kappa
       else
          call abortmp('ERROR: in prim_Bblock_op, moist run')
          ! need to add scaling factor for derivative of Tv
          ! do k = 1,nlev
          !    dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
          !         + fptr%hvcoord%hybi(k+1) * ps_n(:,:,ie)) &
          !         - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
          !         + fptr%hvcoord%hybi(k) * ps_n(:,:,ie))
          ! end do
          
          ! Qt = fptr%base(ie)%state%Qdp(:,:,:,1,qn0)/dp(:,:,:)
          ! Tv = virtual_temperature(temp_n(:,:,:,ie), Qt)
       end if
              
       ! delta surface pressure
       dps = fptr%base(ie)%state%ps_v(:,:,np1)
              
       ! gradient of delta surface pressure
       grad_dps = gradient_sphere(dps, fptr%deriv, fptr%base(ie)%Dinv)
              
       grad_ps = gradient_sphere(ps_n(:,:,ie), fptr%deriv, fptr%base(ie)%Dinv)
              
       ! pressure
       do k = 1,nlev
          press(:,:,k) = fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybm(k) * ps_n(:,:,ie)
       end do
              
       ! vel dot gradient of pressure
       do k = 1,nlev
          grad_p = fptr%hvcoord%hybm(k) * grad_ps
          
          v_grad_p(:,:,k) = vel_n(:,:,1,k,ie) * grad_p(:,:,1) &
               + vel_n(:,:,2,k,ie) * grad_p(:,:,2)
       end do
              
       ! divergence of vertical derivative of pressure * vel                     
       do k = 1,nlev
          dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybi(k+1) * ps_n(:,:,ie)) &
               - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybi(k) * ps_n(:,:,ie))
          
          dpvel(:,:,1) = dp(:,:,k) * vel_n(:,:,1,k,ie)
          dpvel(:,:,2) = dp(:,:,k) * vel_n(:,:,2,k,ie)
          
          div_dpvel(:,:,k) = divergence_sphere(dpvel, fptr%deriv, fptr%base(ie))
       end do
              
       ! pressure vertical velocity
       call preq_omega_ps(omega_p, fptr%hvcoord, press, v_grad_p, div_dpvel)
              
       ! initialize with time contribution
       pstens(:,:,ie) = dti * dps
              
       do k = 1,nlev
          
          ! delta velocity and temperature
          dvel  = fptr%base(ie)%state%v(:,:,:,k,np1)
          dtemp = fptr%base(ie)%state%T(:,:,k,np1)
          
          ! ---------------------------------------------------------------------
          ! compute vtens
          ! ---------------------------------------------------------------------
         
          ! vorticity
          vort  = vorticity_sphere(vel_n(:,:,:,k,ie), fptr%deriv, fptr%base(ie))
          dvort = vorticity_sphere(dvel, fptr%deriv, fptr%base(ie))
          
          ! kinetic energy (ignore for now)
          veldvel = vel_n(:,:,1,k,ie) * dvel(:,:,1) &
                  + vel_n(:,:,2,k,ie) * dvel(:,:,2)
          grad_ke = gradient_sphere(veldvel, fptr%deriv, fptr%base(ie)%Dinv)
          
          ! vertical advection (ignore for now)
          
          ! surface pressure and pressure                    
          pressi  = 1.0d0 / press(:,:,k)
          pressi2 = 1.0d0 / press(:,:,k)**2

          ! geopotential (ignore for now)
          ! call preq_hydrostatic(fptr%base(ie)%derived%phi, & 
          !      fptr%base(ie)%state%phis,T_v,p,dp)
          
          vtens(:,:,1,k,ie) = spheremp(:,:) * &
               (dvel(:,:,1) * dti &
               - vel_n(:,:,2,k,ie) * dvort(:,:) &
               - dvel(:,:,2) * (vort(:,:) + fcor(:,:)) &
               + Rgas * Tv(:,:,k) * beta(k) * &
               (-1.0d0 * pressi2 * beta(k) * grad_ps(:,:,1) * dps &
               + pressi * grad_dps(:,:,1)) &
               + Rgas * pressi * beta(k) * grad_ps(:,:,1) * dtemp(:,:) &
               + grad_ke(:,:,1))
               
          vtens(:,:,2,k,ie) = spheremp(:,:) * &
               (dvel(:,:,2) * dti &
               + vel_n(:,:,1,k,ie) * dvort(:,:) &
               + dvel(:,:,1) * (vort(:,:) + fcor(:,:)) &
               + Rgas * Tv(:,:,k) * beta(k) * &
               (-1.0d0 * pressi2 * beta(k) * grad_ps(:,:,2) * dps &
               + pressi * grad_dps(:,:,2)) &
               + Rgas * pressi * beta(k) * grad_ps(:,:,2) * dtemp(:,:) &
               + grad_ke(:,:,2))

          ! ---------------------------------------------------------------------
          ! compute pstens
          ! ---------------------------------------------------------------------
          dpdvel(:,:,1) = dp(:,:,k) * dvel(:,:,1)
          dpdvel(:,:,2) = dp(:,:,k) * dvel(:,:,2)
          
          div_dpdvel = divergence_sphere(dpdvel, fptr%deriv, fptr%base(ie))
          
          dpsvel(:,:,1) = (fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)) * &
               dps * vel_n(:,:,1,k,ie)
          
          dpsvel(:,:,2) = (fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)) * &
               dps * vel_n(:,:,2,k,ie)
          
          div_dpsvel = divergence_sphere(dpsvel, fptr%deriv, fptr%base(ie))
          
          pstens(:,:,ie) = pstens(:,:,ie) + div_dpdvel + div_dpsvel
          
          ! ---------------------------------------------------------------------
          ! compute ttens
          ! ---------------------------------------------------------------------
          grad_dtemp = gradient_sphere(dtemp, fptr%deriv, fptr%base(ie)%Dinv)

          ! vel dot gradient of temperature
          v_grad_dtemp = vel_n(:,:,1,k,ie) * grad_dtemp(:,:,1) &
               + vel_n(:,:,2,k,ie) * grad_dtemp(:,:,2)

          ! vertical advection of temperature (ignoring for now)

          ttens(:,:,k,ie) = spheremp * &
               (dti * dtemp &
               + v_grad_dtemp &
               - kappa_star(:,:,k) * dtemp(:,:) * omega_p(:,:,k))
       end do

       pstens(:,:,ie) = spheremp(:,:) * pstens(:,:,ie)

    end do

    ! pack edge data
    do ie = nets,nete
       call edgeVpack(edge3p1, pstens(:,:,ie),         1,      0, ie)
       call edgeVpack(edge3p1, ttens(:,:,:,ie),     nlev,      1, ie)
       call edgeVpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, nlev+1, ie)
    end do
    
    ! exchange edge data
    call bndry_exchangeV(fptr%hybrid, edge3p1)
        
    ! unpack edge data
    do ie = nets,nete
       call edgeVunpack(edge3p1, pstens(:,:,ie),         1,      0, ie)
       call edgeVunpack(edge3p1, ttens(:,:,:,ie),     nlev,      1, ie)
       call edgeVunpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, nlev+1, ie)
    end do
        
    do ie = nets,nete
       pstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:) * pstens(:,:,ie)
       
       do k=1,nlev
          ttens(:,:,k,ie)   = fptr%base(ie)%rspheremp(:,:) * ttens(:,:,k,ie)
          vtens(:,:,1,k,ie) = fptr%base(ie)%rspheremp(:,:) * vtens(:,:,1,k,ie)
          vtens(:,:,2,k,ie) = fptr%base(ie)%rspheremp(:,:) * vtens(:,:,2,k,ie)
       end do
    end do
    
    ! initialize 1d array index
    lx = 1
    
    ! copy V1  
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = vtens(i,j,1,k,ie)
                lx = lx+1
             end do
          end do
       end do
    end do
    
    ! copy V2
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = vtens(i,j,2,k,ie)
                lx = lx+1
             end do
          end do
       end do
    end do
    
    ! copy Ps
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             fxstate(lx) = pstens(i,j,ie)
             lx = lx+1
          end do
       end do
    end do

    ! copy T
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = ttens(i,j,k,ie)
                lx = lx+1
             end do
          end do
       end do
    end do

    call t_stopf('prim_approx_jacobian_op')

#ifdef DEBUG_PRINT_ON
    write(*,*) 'End: prim_approx_jacobian_op'
#endif
    
  end subroutine prim_approx_jacobian_op



  ! approximation of full Jacobian operator
  subroutine prim_approx_jacobian_op_2(xstate, nxstate, fxstate, c_ptr_to_object) &
       bind(C,name='prim_approx_jacobian_op_2')
    
    use, intrinsic :: iso_c_binding 
    
    use control_mod, only        : tstep_type
    use physical_constants, only : Rgas, kappa
    use prim_si_mod, only        : preq_omega_ps
    use derivative_mod, only     : vorticity_sphere, gradient_sphere, &
                                   divergence_sphere

    implicit none
    
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(out)       :: fxstate(nxstate)
    type(c_ptr)                       :: c_ptr_to_object
    
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,2,nlev,nelemd) :: vtens
    real(kind=real_kind), dimension(np,np,nelemd)        :: pstens
    real(kind=real_kind), dimension(np,np,nlev,nelemd)   :: ttens

    real(kind=real_kind), dimension(np,np,2) :: dvel
    real(kind=real_kind), dimension(np,np)   :: dps
    real(kind=real_kind), dimension(np,np)   :: dtemp

    real(kind=real_kind), dimension(np,np,nlev) :: Tv
    real(kind=real_kind), dimension(np,np,nlev) :: kappa_star

    real(kind=real_kind), dimension(np,np,2) :: grad_dps
    real(kind=real_kind), dimension(np,np,2) :: grad_ps
    real(kind=real_kind), dimension(np,np,2) :: grad_temp
    real(kind=real_kind), dimension(np,np,2) :: grad_dtemp
    real(kind=real_kind), dimension(np,np,2) :: dpdvel
    real(kind=real_kind), dimension(np,np,2) :: dpsvel

    real(kind=real_kind), dimension(np,np) :: dp
    real(kind=real_kind), dimension(np,np) :: fcor
    real(kind=real_kind), dimension(np,np) :: pressi, pressi2
    real(kind=real_kind), dimension(np,np) :: div_dpdvel
    real(kind=real_kind), dimension(np,np) :: div_dpsvel
    real(kind=real_kind), dimension(np,np) :: dv_grad_temp
    real(kind=real_kind), dimension(np,np) :: v_grad_dtemp
    real(kind=real_kind), dimension(np,np) :: spheremp

    real(kind=real_kind), dimension(nlev) :: beta

    real(kind=real_kind) :: dti

    integer :: nets, nete  ! element start and end index
    integer :: np1         ! time levels
    integer :: i, j, k, ie ! loop counters
    integer :: lx          ! 
    integer :: qn0         ! moisture

#ifdef DEBUG_PRINT_ON
    write(*,*) 'Start: prim_approx_jacobian_op_2'
#endif

    call t_startf('prim_approx_jacobian_op_2')
    
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    np1  = fptr%tl%np1  ! Newton update
    nets = fptr%nets    ! stating element index
    nete = fptr%nete    ! ending element index
    qn0  = fptr%n_Q     ! moist run flag

    dti = 1.0d0/fptr%dt 

    if (tstep_type == 12) then
       dti = (3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2
    end if

    ! initialize 1d array index
    lx = 1

    ! copy V1
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%v(i,j,1,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    ! copy V2
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%v(i,j,2,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    ! copy Ps
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             fptr%base(ie)%state%ps_v(i,j,np1) = xstate(lx)
             lx = lx + 1
          end do
       end do
    end do

    ! copy V2
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%T(i,j,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    ! Ttemp    = 1.0d0          
    ! Tvfactor = 1.0d0

    do k = 1,nlev
       beta(k) = fptr%hvcoord%hybm(k)
    end do

    ! compute Jacobian matvec with spectral element operators
    do ie = nets,nete

       ! DSS weights
       spheremp = fptr%base(ie)%spheremp(:,:)
       
       ! coriolis force
       fcor = fptr%base(ie)%fcor(:,:)
       
       ! compute T_v for timelevel n0
       if ( qn0 == -1 ) then
          Tv(:,:,:)  = temp_n(:,:,:,ie)
          kappa_star = kappa
       else
          call abortmp('ERROR: in prim_Bblock_op, moist run')
          ! need to add scaling factor for derivative of Tv
       end if
              
       ! delta surface pressure
       dps = fptr%base(ie)%state%ps_v(:,:,np1)
              
       ! gradient of delta surface pressure
       grad_dps = gradient_sphere(dps, fptr%deriv, fptr%base(ie)%Dinv)
       grad_ps  = gradient_sphere(ps_n(:,:,ie), fptr%deriv, fptr%base(ie)%Dinv)
                                   
       ! initialize with time contribution
       pstens(:,:,ie) = dti * dps
              
       do k = 1,nlev
          
          ! delta velocity and temperature
          dvel  = fptr%base(ie)%state%v(:,:,:,k,np1)
          dtemp = fptr%base(ie)%state%T(:,:,k,np1)
          
          ! ---------------------------------------------------------------------
          ! compute vtens
          ! ---------------------------------------------------------------------
          
          pressi = 1.0d0 / (fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 &
                 + fptr%hvcoord%hybm(k) * ps_n(:,:,ie))
                   
          pressi2 = pressi**2
          
          vtens(:,:,1,k,ie) = spheremp(:,:) * &
               (dvel(:,:,1) * dti &
               - dvel(:,:,2) * fcor(:,:) &
               + Rgas * Tv(:,:,k) * beta(k) * (-1.0d0 * pressi2 * beta(k) * grad_ps(:,:,1) * dps + pressi * grad_dps(:,:,1)) &
               + Rgas * pressi * beta(k) * grad_ps(:,:,1) * dtemp(:,:))
               
          vtens(:,:,2,k,ie) = spheremp(:,:) * &
               (dvel(:,:,2) * dti &
               + dvel(:,:,1) * fcor(:,:) &
               + Rgas * Tv(:,:,k) * beta(k) * (-1.0d0 * pressi2 * beta(k) * grad_ps(:,:,2) * dps + pressi * grad_dps(:,:,2)) &
               + Rgas * pressi * beta(k) * grad_ps(:,:,2) * dtemp(:,:))

          ! ---------------------------------------------------------------------
          ! compute pstens
          ! ---------------------------------------------------------------------
          dp = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
             + fptr%hvcoord%hybi(k+1) * ps_n(:,:,ie)) &
             - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
             + fptr%hvcoord%hybi(k) * ps_n(:,:,ie))

          dpdvel(:,:,1) = dp * dvel(:,:,1)
          dpdvel(:,:,2) = dp * dvel(:,:,2)
          
          div_dpdvel = divergence_sphere(dpdvel, fptr%deriv, fptr%base(ie))
          
          dpsvel(:,:,1) = (fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)) * &
               dps * vel_n(:,:,1,k,ie)

          dpsvel(:,:,2) = (fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)) * &
               dps * vel_n(:,:,2,k,ie)
          
          div_dpsvel = divergence_sphere(dpsvel, fptr%deriv, fptr%base(ie))

          pstens(:,:,ie) = pstens(:,:,ie) + div_dpdvel + div_dpsvel

          ! ---------------------------------------------------------------------
          ! compute ttens
          ! ---------------------------------------------------------------------
          grad_temp  = gradient_sphere(temp_n(:,:,k,ie), fptr%deriv, fptr%base(ie)%Dinv)
          grad_dtemp = gradient_sphere(dtemp, fptr%deriv, fptr%base(ie)%Dinv)

          ! vel dot gradient of temperature
          v_grad_dtemp = vel_n(:,:,1,k,ie) * grad_dtemp(:,:,1) &
               + vel_n(:,:,2,k,ie) * grad_dtemp(:,:,2)

          dv_grad_temp = dvel(:,:,1) * grad_temp(:,:,1) &
               + dvel(:,:,2) * grad_temp(:,:,2)

          ! vertical advection of temperature (ignoring for now)

          ttens(:,:,k,ie) = spheremp * &
               (dti * dtemp + dv_grad_temp + v_grad_dtemp)
       end do

       pstens(:,:,ie) = spheremp(:,:) * pstens(:,:,ie)

    end do

    ! pack edge data
    do ie = nets,nete
       call edgeVpack(edge3p1, pstens(:,:,ie),         1,      0, ie)
       call edgeVpack(edge3p1, ttens(:,:,:,ie),     nlev,      1, ie)
       call edgeVpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, nlev+1, ie)
    end do
    
    ! exchange edge data
    call bndry_exchangeV(fptr%hybrid, edge3p1)
        
    ! unpack edge data
    do ie = nets,nete
       call edgeVunpack(edge3p1, pstens(:,:,ie),         1,      0, ie)
       call edgeVunpack(edge3p1, ttens(:,:,:,ie),     nlev,      1, ie)
       call edgeVunpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, nlev+1, ie)
    end do
        
    do ie = nets,nete
       pstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:) * pstens(:,:,ie)
       
       do k=1,nlev
          ttens(:,:,k,ie)   = fptr%base(ie)%rspheremp(:,:) * ttens(:,:,k,ie)
          vtens(:,:,1,k,ie) = fptr%base(ie)%rspheremp(:,:) * vtens(:,:,1,k,ie)
          vtens(:,:,2,k,ie) = fptr%base(ie)%rspheremp(:,:) * vtens(:,:,2,k,ie)
       end do
    end do
    
    ! initialize 1d array index
    lx = 1
    
    ! copy V1  
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = vtens(i,j,1,k,ie)
                lx = lx+1
             end do
          end do
       end do
    end do
    
    ! copy V2
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = vtens(i,j,2,k,ie)
                lx = lx+1
             end do
          end do
       end do
    end do
    
    ! copy Ps
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             fxstate(lx) = pstens(i,j,ie)
             lx = lx+1
          end do
       end do
    end do

    ! copy T
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = ttens(i,j,k,ie)
                lx = lx+1
             end do
          end do
       end do
    end do

    call t_stopf('prim_approx_jacobian_op_2')

#ifdef DEBUG_PRINT_ON
    write(*,*) 'End: prim_approx_jacobian_op_2'
#endif
    
  end subroutine prim_approx_jacobian_op_2



  ! approximation of full Jacobian operator
  subroutine prim_approx_jacobian_op_3(xstate, nxstate, fxstate, c_ptr_to_object) &
       bind(C,name='prim_approx_jacobian_op_3')
    
    use, intrinsic :: iso_c_binding 
    
    use control_mod, only        : tstep_type
    use physical_constants, only : Rgas, kappa
    use prim_si_mod, only        : preq_omega_ps
    use derivative_mod, only     : vorticity_sphere, gradient_sphere, &
                                   divergence_sphere

    implicit none
    
    real(c_double), intent(in)        :: xstate(nxstate)
    integer(c_int), intent(in), value :: nxstate
    real(c_double), intent(out)       :: fxstate(nxstate)
    type(c_ptr)                       :: c_ptr_to_object
    
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,2,nlev,nelemd) :: vtens
    real(kind=real_kind), dimension(np,np,nelemd)        :: pstens
    real(kind=real_kind), dimension(np,np,nlev,nelemd)   :: ttens

    real(kind=real_kind), dimension(np,np,2) :: dvel
    real(kind=real_kind), dimension(np,np)   :: dps
    real(kind=real_kind), dimension(np,np)   :: dtemp

    real(kind=real_kind), dimension(np,np,nlev) :: Tv
    real(kind=real_kind), dimension(np,np,nlev) :: kappa_star

    real(kind=real_kind), dimension(np,np,2) :: grad_dps
    real(kind=real_kind), dimension(np,np,2) :: grad_ps
    real(kind=real_kind), dimension(np,np,2) :: grad_temp
    real(kind=real_kind), dimension(np,np,2) :: grad_dtemp
    real(kind=real_kind), dimension(np,np,2) :: dpdvel
    real(kind=real_kind), dimension(np,np,2) :: dpsvel

    real(kind=real_kind), dimension(np,np) :: dp
    real(kind=real_kind), dimension(np,np) :: fcor
    real(kind=real_kind), dimension(np,np) :: pressi, pressi2
    real(kind=real_kind), dimension(np,np) :: div_dpdvel
    real(kind=real_kind), dimension(np,np) :: div_dpsvel
    real(kind=real_kind), dimension(np,np) :: dv_grad_temp
    real(kind=real_kind), dimension(np,np) :: v_grad_dtemp
    real(kind=real_kind), dimension(np,np) :: spheremp

    real(kind=real_kind), dimension(nlev) :: beta

    real(kind=real_kind) :: dti

    integer :: nets, nete  ! element start and end index
    integer :: np1         ! time levels
    integer :: i, j, k, ie ! loop counters
    integer :: lx          ! 
    integer :: qn0         ! moisture

#ifdef DEBUG_PRINT_ON
    write(*,*) 'Start: prim_approx_jacobian_op_3'
#endif

    call t_startf('prim_approx_jacobian_op_3')
    
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    np1  = fptr%tl%np1  ! Newton update
    nets = fptr%nets    ! stating element index
    nete = fptr%nete    ! ending element index
    qn0  = fptr%n_Q     ! moist run flag

    dti = 1.0d0/fptr%dt 

    if (tstep_type == 12) then
       dti = (3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2
    end if

    ! initialize 1d array index
    lx = 1

    ! copy V1
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%v(i,j,1,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    ! copy V2
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%v(i,j,2,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    ! copy Ps
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             fptr%base(ie)%state%ps_v(i,j,np1) = xstate(lx)
             lx = lx + 1
          end do
       end do
    end do

    ! copy V2
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fptr%base(ie)%state%T(i,j,k,np1) = xstate(lx)
                lx = lx + 1
             end do
          end do
       end do
    end do

    ! Ttemp    = 1.0d0          
    ! Tvfactor = 1.0d0

    do k = 1,nlev
       beta(k) = fptr%hvcoord%hybm(k)
    end do

    ! compute Jacobian matvec with spectral element operators
    do ie = nets,nete

       ! DSS weights
       spheremp = fptr%base(ie)%spheremp(:,:)
       
       ! coriolis force
       fcor = fptr%base(ie)%fcor(:,:)
       
       ! compute T_v for timelevel n0
       if ( qn0 == -1 ) then
          Tv(:,:,:)  = temp_n(:,:,:,ie)
          kappa_star = kappa
       else
          call abortmp('ERROR: in prim_Bblock_op, moist run')
          ! need to add scaling factor for derivative of Tv
       end if
              
       ! delta surface pressure
       dps = fptr%base(ie)%state%ps_v(:,:,np1)
              
       ! gradient of delta surface pressure
       grad_dps = gradient_sphere(dps, fptr%deriv, fptr%base(ie)%Dinv)
       grad_ps  = gradient_sphere(ps_n(:,:,ie), fptr%deriv, fptr%base(ie)%Dinv)
                                   
       ! initialize with time contribution
       pstens(:,:,ie) = dti * dps
              
       do k = 1,nlev
          
          ! delta velocity and temperature
          dvel  = fptr%base(ie)%state%v(:,:,:,k,np1)
          dtemp = fptr%base(ie)%state%T(:,:,k,np1)
          
          ! ---------------------------------------------------------------------
          ! compute vtens
          ! ---------------------------------------------------------------------
          
          pressi = 1.0d0 / (fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 &
                 + fptr%hvcoord%hybm(k) * ps_n(:,:,ie))
                   
          pressi2 = pressi**2
          
          vtens(:,:,1,k,ie) = spheremp(:,:) * &
               (dvel(:,:,1) * dti &
               - dvel(:,:,2) * fcor(:,:) &
               + Rgas * Tv(:,:,k) * beta(k) * (-1.0d0 * pressi2 * beta(k) * grad_ps(:,:,1) * dps + pressi * grad_dps(:,:,1)) &
               + Rgas * pressi * beta(k) * grad_ps(:,:,1) * dtemp(:,:))

               
          vtens(:,:,2,k,ie) = spheremp(:,:) * &
               (dvel(:,:,2) * dti &
               + dvel(:,:,1) * fcor(:,:) &
               + Rgas * Tv(:,:,k) * beta(k) * (-1.0d0 * pressi2 * beta(k) * grad_ps(:,:,2) * dps + pressi * grad_dps(:,:,2)) &
               + Rgas * pressi * beta(k) * grad_ps(:,:,2) * dtemp(:,:))

          ! ---------------------------------------------------------------------
          ! compute pstens
          ! ---------------------------------------------------------------------
          pstens(:,:,ie) = pstens(:,:,ie) + &
               beta(k) * grad_dps(:,:,1) * vel_n(:,:,1,k,ie) + &
               beta(k) * grad_dps(:,:,2) * vel_n(:,:,2,k,ie)
          
          ! ---------------------------------------------------------------------
          ! compute ttens
          ! ---------------------------------------------------------------------
          grad_temp  = gradient_sphere(temp_n(:,:,k,ie), fptr%deriv, fptr%base(ie)%Dinv)
          grad_dtemp = gradient_sphere(dtemp, fptr%deriv, fptr%base(ie)%Dinv)

          ! vel dot gradient of temperature
          v_grad_dtemp = vel_n(:,:,1,k,ie) * grad_dtemp(:,:,1) &
               + vel_n(:,:,2,k,ie) * grad_dtemp(:,:,2)

          dv_grad_temp = dvel(:,:,1) * grad_temp(:,:,1) &
               + dvel(:,:,2) * grad_temp(:,:,2)

          ! vertical advection of temperature (ignoring for now)

          ttens(:,:,k,ie) = spheremp * &
               (dti * dtemp + dv_grad_temp + v_grad_dtemp)
       end do

       pstens(:,:,ie) = spheremp(:,:) * pstens(:,:,ie)

    end do

    ! pack edge data
    do ie = nets,nete
       call edgeVpack(edge3p1, pstens(:,:,ie),         1,      0, ie)
       call edgeVpack(edge3p1, ttens(:,:,:,ie),     nlev,      1, ie)
       call edgeVpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, nlev+1, ie)
    end do
    
    ! exchange edge data
    call bndry_exchangeV(fptr%hybrid, edge3p1)
        
    ! unpack edge data
    do ie = nets,nete
       call edgeVunpack(edge3p1, pstens(:,:,ie),         1,      0, ie)
       call edgeVunpack(edge3p1, ttens(:,:,:,ie),     nlev,      1, ie)
       call edgeVunpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, nlev+1, ie)
    end do
        
    do ie = nets,nete
       pstens(:,:,ie) = fptr%base(ie)%rspheremp(:,:) * pstens(:,:,ie)
       
       do k=1,nlev
          ttens(:,:,k,ie)   = fptr%base(ie)%rspheremp(:,:) * ttens(:,:,k,ie)
          vtens(:,:,1,k,ie) = fptr%base(ie)%rspheremp(:,:) * vtens(:,:,1,k,ie)
          vtens(:,:,2,k,ie) = fptr%base(ie)%rspheremp(:,:) * vtens(:,:,2,k,ie)
       end do
    end do
    
    ! initialize 1d array index
    lx = 1
    
    ! copy V1  
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = vtens(i,j,1,k,ie)
                lx = lx+1
             end do
          end do
       end do
    end do
    
    ! copy V2
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = vtens(i,j,2,k,ie)
                lx = lx+1
             end do
          end do
       end do
    end do
    
    ! copy Ps
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             fxstate(lx) = pstens(i,j,ie)
             lx = lx+1
          end do
       end do
    end do

    ! copy T
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                fxstate(lx) = ttens(i,j,k,ie)
                lx = lx+1
             end do
          end do
       end do
    end do

    call t_stopf('prim_approx_jacobian_op_3')

#ifdef DEBUG_PRINT_ON
    write(*,*) 'End: prim_approx_jacobian_op_3'
#endif
    
  end subroutine prim_approx_jacobian_op_3
  
  
  
  ! subroutine prim_diag_prec(xstate, nxstate, fxstate, c_ptr_to_object) &
  !      bind(C,name='prim_diag_prec')
  !   ! ---------------------------------------------------------------------------
  !   ! scaled identity function
  !   ! ---------------------------------------------------------------------------
  !   use, intrinsic :: iso_c_binding 
    
  !   use kinds, only                 : real_kind
  !   use prim_derived_type_mod, only : derived_type
  !   use control_mod, only           : tstep_type
    
  !   implicit none 
    
  !   integer(c_int), intent(in), value :: nxstate
  !   real(c_double), intent(in)        :: xstate(nxstate)
  !   real(c_double), intent(out)       :: fxstate(nxstate)
  !   type(c_ptr)                       :: c_ptr_to_object
    
  !   type(derived_type), pointer :: fptr=>NULL()
    
  !   real(kind=real_kind) :: dti

  !   call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

  !   dti = 1.0d0/fptr%dt 
        
  !   if (tstep_type == 12) then
  !      dti = (3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2
  !   end if

  !   fxstate = dt * xstate
    
  ! end subroutine prim_dt_ident


end module prim_jacobian_operator_mod

#endif
! TRILINOS
