#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
  
#ifdef TRILINOS

!#define DEBUG_PRINT_ON

! ==============================================================================
! This module contains subroutines for computing the entries of the analytic 
! Jacobian of the nonlinear residual equation of the hydrostatic primitive 
! equations. 
! ==============================================================================
  
module prim_jacobian_sparse_mod
  
  use kinds, only          : real_kind, long_kind
  use dimensions_mod, only : np, nlev
  use edgetype_mod, only   : EdgeBuffer_t
  use parallel_mod, only   : haltmp, abortmp, iam
  use perf_mod, only       : t_startf, t_stopf
  
  implicit none
  private

#ifdef _MPI
#include <mpif.h>
#endif

  public :: prim_jacobian_sparse_init, prim_jacobian_sparse_finalize
  public :: print_global_IDs
  public :: uoffset, voffset, toffset, psoffset

  ! global index offsets from gdofp for 3D variables (long_kind = 8)
  integer(kind=long_kind), dimension(nlev), save :: Uoffset, Voffset, Toffset

  ! global index offsets from gdofp for 2D variables (long_kind = 8)
  integer(kind=long_kind), save :: Psoffset

contains

  subroutine prim_jacobian_sparse_init(par, elem)
    ! ---------------------------------------------------------------------------
    ! Computes a global ID value for each unknown (U, V, Ps, T) based on gdofp
    ! in order to construct the global sparse Jacobian matrix. Variables that are
    ! shared at element edges and corners are given the same ID. 
    !
    ! Since gdofp stores a unique saptial node ID for one level of the mesh this 
    ! subroutines computes the necessary offset values to create a unique ID for 
    ! all variables.
    ! ---------------------------------------------------------------------------
    use kinds, only          : long_kind   
    use dimensions_mod, only : np, nlev, nelemd
    use element_mod, only    : element_t
    use control_mod, only    : rsplit
    use parallel_mod, only   : parallel_t
 
    implicit none

    ! --------------------------------Arguments----------------------------------
    type(parallel_t), intent(in) :: par
    type(element_t),  intent(in), target :: elem(:)
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    integer(kind=long_kind), dimension(nelemd) :: LocalMaxVals

    integer(kind=long_kind) :: LocalMaxDof_2D, GlobalMaxDof_2D
    integer(kind=long_kind) :: VariableOffset_3D, LevelOffset

    integer :: ie, k, ierr
    ! --------------------------------------------------------------------------- 
    
    ! find the global max node ID, this is not the same as the number of unique 
    ! nodes in the mesh since some ID values are skipped in gdofp.
    do ie = 1,nelemd
       LocalMaxVals(ie) = maxval(elem(ie)%gdofp(:,:))
    end do
    LocalMaxDof_2D = maxval(LocalMaxVals)

    call MPI_Allreduce(LocalMaxDof_2D, GlobalMaxDof_2D, 1, &
         MPI_INTEGER8, MPI_MAX, par%comm, ierr)

    ! ID offset between different 3D variables
    VariableOffset_3D = nlev * GlobalMaxDof_2D

    do k = 1,nlev
       ! ID offset between the same variable on different levels
       LevelOffset = (k-1) * GlobalMaxDof_2D

       ! offset from gdofp to get unique ID for each variable
       Uoffset(k) = LevelOffset
       Voffset(k) = VariableOffset_3D + LevelOffset
       Toffset(k) = 2 * VariableOffset_3D + GlobalMaxDof_2D + LevelOffset
    end do

    Psoffset = 2 * VariableOffset_3D

  end subroutine prim_jacobian_sparse_init
 
  
  
  subroutine prim_jacobian_sparse_finalize()
    ! ---------------------------------------------------------------------------
    ! Dummy subrutine does nothing for sparse jacobian
    ! ---------------------------------------------------------------------------
    implicit none

  end subroutine prim_jacobian_sparse_finalize


  subroutine print_global_IDs(fptr)
    ! ---------------------------------------------------------------------------
    ! Outputs variable offsets to file
    ! ---------------------------------------------------------------------------
    use dimensions_mod, only        : np, nlev, nelemd
    use prim_derived_type_mod, only : derived_type

    implicit none
    
    ! --------------------------------Arguments----------------------------------
    type(derived_type), pointer, intent(in) :: fptr
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    integer :: i, j, k, ie ! equation loop
    integer :: idx         ! row and column index
    integer :: fid         ! file id for output
    
    character(len=1024) :: fname ! output file name
    ! ---------------------------------------------------------------------------

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*) "Writing global IDs to file"
#endif

    fid = 1000 + iam 
    write(fname,'(a,i2.2,a,i2.2,a)') 'U_globalIDs_proc',iam,'.txt'
    open(fid, file=fname, form='formatted')

    do ie = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

                idx = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
                write(fid,'(5i23)') i,j,k,ie,idx 

             end do
          end do
       end do
    end do

    close(fid)

    ! ---------------------------------------------------------------------------
    write(fname,'(a,i2.2,a,i2.2,a)') 'V_globalIDs_proc',iam,'.txt'
    open(fid, file=fname, form='formatted')

    do ie = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

                idx = fptr%base(ie)%gdofP(i,j) + Voffset(k)
                write(fid,'(5i23)') i,j,k,ie,idx 

             end do
          end do
       end do
    end do

    close(fid)

    ! ---------------------------------------------------------------------------
    write(fname,'(a,i2.2,a,i2.2,a)') 'Ps_globalIDs_proc',iam,'.txt'
    open(fid, file=fname, form='formatted')

    do ie = 1,nelemd
       do j = 1,np
          do i = 1,np
             
             idx = fptr%base(ie)%gdofP(i,j) + Psoffset
             write(fid,'(4i23)') i,j,ie,idx 
             
          end do
       end do
    end do

    close(fid)

    ! ---------------------------------------------------------------------------
    write(fname,'(a,i2.2,a,i2.2,a)') 'T_globalIDs_proc',iam,'.txt'
    open(fid, file=fname, form='formatted')

    do ie = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

                idx = fptr%base(ie)%gdofP(i,j) + Toffset(k)
                write(fid,'(5i23)') i,j,k,ie,idx 

             end do
          end do
       end do
    end do

    close(fid)

  end subroutine print_global_IDs



  ! =============================================================================
  ! Private subroutines for building Trilinos matrix maps 
  ! =============================================================================



  subroutine GlobalIDs_Vel(gidx, nidx, c_ptr_to_object) &
       bind(C,name='GlobalIDs_Vel')
    ! ---------------------------------------------------------------------------
    ! Global IDs for Velocity (U, V) variables on this process
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use prim_derived_type_mod, only : derived_type
    use dimensions_mod, only        : np, nlev

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: nidx ! number of global indices

    integer(c_int), intent(out) :: gidx(nidx) ! global indices

    type(c_ptr) :: c_ptr_to_object ! HOMME data
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    integer :: ie, nets, nete
    integer :: i, j, k
    integer :: idx
    !----------------------------------------------------------------------------

    call t_startf('GlobalIDs_Vel')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*) "Computing global velocity IDs"
#endif

    ! elements on this process
    nets = fptr%nets
    nete = fptr%nete

    ! Initialize arrays and storage index
    gidx = 0
    idx  = 0

    ! global indices for U unknowns
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                
                ! increment storage index, get global row/column index
                idx = idx + 1
                gidx(idx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
                
             end do
          end do
       end do
    end do

    ! global indices for V unknowns
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                
                ! increment storage index, get global row/column index
                idx = idx + 1
                gidx(idx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)
                
             end do
          end do
       end do
    end do

    call t_stopf('GlobalIDs_Vel')

  end subroutine GlobalIDs_Vel



  subroutine GlobalIDs_Ps(gidx, nidx, c_ptr_to_object) &
       bind(C,name='GlobalIDs_Ps')
    ! ---------------------------------------------------------------------------
    ! Global IDs for Surface Pressure (Ps) variables on this process
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use prim_derived_type_mod, only : derived_type
    use dimensions_mod, only        : np, nlev

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: nidx ! number of global indices

    integer(c_int), intent(out) :: gidx(nidx) ! global indices

    type(c_ptr) :: c_ptr_to_object ! HOMME data
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    integer :: ie, nets, nete
    integer :: i, j
    integer :: idx
    !----------------------------------------------------------------------------

    call t_startf('GlobalIDs_Ps')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*) "Computing global surface pressure IDs"
#endif

    ! elements on this process
    nets = fptr%nets
    nete = fptr%nete

    ! Initialize arrays and storage index
    gidx = 0
    idx  = 0

    ! global indices for Ps unknowns
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             
             ! increment storage index, get global row/column index
             idx = idx + 1
             gidx(idx) = fptr%base(ie)%gdofP(i,j) + Psoffset
             
          end do
       end do
    end do

    call t_stopf('GlobalIDs_Ps')

  end subroutine GlobalIDs_Ps



  subroutine GlobalIDs_T(gidx, nidx, c_ptr_to_object) &
       bind(C,name='GlobalIDs_T')
    ! ---------------------------------------------------------------------------
    ! Global IDs for Temperature (T) variables on this process
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use prim_derived_type_mod, only : derived_type
    use dimensions_mod, only        : np, nlev

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: nidx ! number of global indices

    integer(c_int), intent(out) :: gidx(nidx) ! global indices

    type(c_ptr) :: c_ptr_to_object ! HOMME data
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    integer :: ie, nets, nete
    integer :: i, j, k
    integer :: idx
    !----------------------------------------------------------------------------

    call t_startf('GlobalIDs_T')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*) "Computing global temperature IDs"
#endif

    ! elements on this process
    nets = fptr%nets
    nete = fptr%nete

    ! Initialize arrays and storage index
    gidx = 0
    idx  = 0

    ! global indices for T unknowns
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                
                ! increment storage index, get global row/column index
                idx = idx + 1
                gidx(idx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)
                
             end do
          end do
       end do
    end do

    call t_stopf('GlobalIDs_T')

  end subroutine GlobalIDs_T



  ! =============================================================================
  ! Private subroutines for setting node weights
  ! =============================================================================
  
  
  
  subroutine Weights_Vel(weights, nwgt, c_ptr_to_object) &
       bind(C,name='Weights_Vel')
    ! ---------------------------------------------------------------------------
    ! Weights for Velocity (U, V) variables on this process
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use prim_derived_type_mod, only : derived_type
    use dimensions_mod, only        : np, nlev

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: nwgt ! number of weights

    real(c_double), intent(out) :: weights(nwgt) ! weights

    type(c_ptr) :: c_ptr_to_object ! HOMME data
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    integer :: ie, nets, nete
    integer :: i, j, k
    integer :: idx
    !----------------------------------------------------------------------------

    call t_startf('Weights_Vel')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*) "Computing velocity weights"
#endif

    ! elements on this process
    nets = fptr%nets
    nete = fptr%nete

    ! initialize arrays and storage index
    weights = 1.0d0
    idx = 0

    ! weights for U unknowns
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                
                ! increment storage index, set weight
                idx = idx + 1
                weights(idx) = fptr%base(ie)%spheremp(i,j) &
                             * fptr%base(ie)%rspheremp(i,j)
                
             end do
          end do
       end do
    end do

    ! weights for V unknowns
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                
                ! increment storage index, set weight
                idx = idx + 1
                weights(idx) = fptr%base(ie)%spheremp(i,j) &
                             * fptr%base(ie)%rspheremp(i,j)
                
             end do
          end do
       end do
    end do

    call t_stopf('Weights_Vel')

  end subroutine Weights_Vel



  subroutine Weights_Ps(weights, nwgt, c_ptr_to_object) &
       bind(C,name='Weights_Ps')
    ! ---------------------------------------------------------------------------
    ! Weights for surface pressure variables on this process
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use prim_derived_type_mod, only : derived_type
    use dimensions_mod, only        : np, nlev

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: nwgt ! number of weights

    real(c_double), intent(out) :: weights(nwgt) ! weights

    type(c_ptr) :: c_ptr_to_object ! HOMME data
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    integer :: ie, nets, nete
    integer :: i, j
    integer :: idx
    !----------------------------------------------------------------------------

    call t_startf('Weights_Ps')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*) "Computing Surface Pressure weights"
#endif

    ! elements on this process
    nets = fptr%nets
    nete = fptr%nete

    ! initialize arrays and storage index
    weights = 1.0d0
    idx = 0

    ! weights for Ps unknowns
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             
             ! increment storage index, set weight
             idx = idx + 1
             weights(idx) = fptr%base(ie)%spheremp(i,j) &
                          * fptr%base(ie)%rspheremp(i,j)
             
          end do
       end do
    end do
    
    call t_stopf('Weights_Ps')

  end subroutine Weights_Ps



  subroutine Weights_T(weights, nwgt, c_ptr_to_object) &
       bind(C,name='Weights_T')
    ! ---------------------------------------------------------------------------
    ! Weights for temperature variables on this process
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use prim_derived_type_mod, only : derived_type
    use dimensions_mod, only        : np, nlev

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: nwgt ! number of weights

    real(c_double), intent(out) :: weights(nwgt) ! weights

    type(c_ptr) :: c_ptr_to_object ! HOMME data
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    integer :: ie, nets, nete
    integer :: i, j, k
    integer :: idx
    !----------------------------------------------------------------------------

    call t_startf('Weights_T')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*) "Computing temperature weights"
#endif

    ! elements on this process
    nets = fptr%nets
    nete = fptr%nete
    
    ! initialize arrays and storage index
    weights = 1.0d0
    idx = 0

    ! weights for T unknowns
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                
                ! increment storage index, set weight
                idx = idx + 1
                weights(idx) = fptr%base(ie)%spheremp(i,j) &
                             * fptr%base(ie)%rspheremp(i,j)
                
             end do
          end do
       end do
    end do

    call t_stopf('Weights_T')

  end subroutine Weights_T
  
  
  
  ! =============================================================================
  ! Private subroutines for computing Jacobian matrix entries
  !
  !      [ du/du   du/dv  | du/dps  | du/dT ] 
  !      [                |         |       ]   
  !      [ dv/du   dv/dv  | dv/dps  | dv/dT ]   [ A | B | C ]
  !  J = [----------------------------------] = [ D | E | - ]
  !      [ dps/du  dps/dv | dps/dps |   0   ]   [ F | G | H ]
  !      [----------------------------------]
  !      [ dT/du   dT/dv  | dT/dps  | dT/dT ]
  !
  ! =============================================================================

  ! =============================================================================
  ! A Block of the Jacobian
  ! =============================================================================
  subroutine Jacobian_Ablock(ie, nnz_max, row_nnz, Urows, Ucols, Uvals, &
       Vrows, Vcols, Vvals, c_ptr_to_object) bind(C,name="Jacobian_Ablock")
    !----------------------------------------------------------------------------
    ! Computes the A block of the Jacobian matrix
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type

    use control_mod, only           : tstep_type
    use derivative_mod, only        : divergence_sphere, vorticity_sphere
    use physical_constants, only    : rrearth
    
    implicit none

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    ! row indices 
    integer(c_int), dimension(np*np*nlev), intent(out) :: Urows
    integer(c_int), dimension(np*np*nlev), intent(out) :: Vrows

    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Ucols
    integer(c_int), dimension(nnz_max), intent(out) :: Vcols

    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Uvals   
    real(c_double), dimension(nnz_max), intent(out) :: Vvals      

    type(c_ptr) :: c_ptr_to_object ! HOMME data
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    ! u and v derivatives of etadot_dpdn
    real(kind=real_kind), dimension(np,np,nlev,np,np,nlev+1) :: du_etadot_dpdn
    real(kind=real_kind), dimension(np,np,nlev,np,np,nlev+1) :: dv_etadot_dpdn

    real(kind=real_kind), dimension(np,np,nlev+1) :: etadot_dpdn

    ! vertical derivative of pressure
    real(kind=real_kind), dimension(np,np,nlev) :: dpdn         
    
    ! 0.5 * 1/dpdn
    real(kind=real_kind), dimension(np,np,nlev) :: half_rdpdn   
    
    ! divergence(dpdn * velocity) 
    real(kind=real_kind), dimension(np,np,nlev) :: div_dpdn_vel 

    ! dpdn * velocity  
    real(kind=real_kind), dimension(np,np,2) :: dpdn_vel 

    real(kind=real_kind), dimension(np,np) :: sdot_sum
    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    real(kind=real_kind), dimension(np,np,nlev) :: zeta ! vorticity

    real(kind=real_kind) :: dti ! 1/dt
    real(kind=real_kind) :: temp1_Vort, temp2_Vort, temp3_Vort ! temp vort vars
    real(kind=real_kind) :: temp1_KE, temp2_KE                 ! temp KE vars

    real(kind=real_kind) :: half_rdpdn_diffu1, half_rdpdn_diffu2
    real(kind=real_kind) :: half_rdpdn_diffv1, half_rdpdn_diffv2
    real(kind=real_kind) :: temp1, temp2
    
    integer :: np1        ! current time level 
    integer :: i, j, k    ! local equation index
    integer :: a, b, c    ! local derivative index
    integer :: ridx, cidx ! index for output U and V output arrays
    !----------------------------------------------------------------------------
    
    call t_startf('Jacobian_ABlock')    

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*) "Computing Jacobian_Ablock"
#endif

    ! initialize arrays
    Urows = 0
    Ucols = 0
    Uvals = 0.0d0

    Vrows = 0
    Vcols = 0
    Vvals = 0.0d0

    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = 2*(np+np-1)*nlev

    ! time level where current iteration is stored
    np1 = fptr%tl%np1

    ! inverse of time step size
    dti = 1.0d0 / fptr%dt

    if (tstep_type==12) then
       dti = (3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2 on step size
    endif

    ! vorticity 
    do k = 1,nlev
       zeta(:,:,k) = vorticity_sphere(fptr%base(ie)%state%v(:,:,:,k,np1), fptr%deriv, fptr%base(ie))
    end do
       
    ! vertical derivative of pressure
    do k = 1,nlev
       dpdn(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
            - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))
    end do
    
    half_rdpdn = 0.5d0 / dpdn
    
    ! divergence(dpdn * velocity)
    do k = 1,nlev
       do j = 1,np
          do i = 1,np             
             dpdn_vel(i,j,1) = dpdn(i,j,k) * fptr%base(ie)%state%v(i,j,1,k,np1)
             dpdn_vel(i,j,2) = dpdn(i,j,k) * fptr%base(ie)%state%v(i,j,2,k,np1)
          end do
       end do
       
       div_dpdn_vel(:,:,k) = divergence_sphere(dpdn_vel, fptr%deriv, fptr%base(ie))
    end do
    
    ! compute etadot_dpdn
    sdot_sum = 0.0d0

    do k=1,nlev
       sdot_sum(:,:)        = sdot_sum(:,:) + div_dpdn_vel(:,:,k)
       etadot_dpdn(:,:,k+1) = sdot_sum(:,:)
    end do

    do k=1,nlev-1
       etadot_dpdn(:,:,k+1) = fptr%hvcoord%hybi(k+1)*sdot_sum(:,:) &
            - etadot_dpdn(:,:,k+1)
    end do

    etadot_dpdn(:,:,1)      = 0.0d0
    etadot_dpdn(:,:,nlev+1) = 0.0d0

    ! velocity derivatives of etadot_dpdn
    call DerivVel_etadot_dpdn(ie, fptr, du_etadot_dpdn, dv_etadot_dpdn)

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp

    ! ---------------------------------------------------------------------------
    ! compute Jacobian entires 
    ! ---------------------------------------------------------------------------

    ! Level 1
    ! ---------------------------------------------------------------------------
    k = 1
    do j = 1,np
       do i = 1,np
          
          ! get global row index for U and V equation
          ridx = ridx + 1            
          Urows(ridx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
          Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)
          
          half_rdpdn_diffu1 = half_rdpdn(i,j,k) * &
               (fptr%base(ie)%state%v(i,j,1,k+1,np1) &
               - fptr%base(ie)%state%v(i,j,1,k,np1))
          
          half_rdpdn_diffv1 = half_rdpdn(i,j,k) &
               * (fptr%base(ie)%state%v(i,j,2,k+1,np1) &
               - fptr%base(ie)%state%v(i,j,2,k,np1))
          
          temp1_Vort = rrearth * fptr%base(ie)%rmetdet(i,j)
          
          do c = 1,nlev

             ! a /= i, b = j, c = * (skips zeros in d*_etadot_dpdn)
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                if (c == k) then
                   temp2_Vort = temp1_Vort * fptr%base(ie)%D(a,j,1,2) &
                                * fptr%deriv%Dvv(a,i)  

                   temp3_Vort = temp1_Vort * fptr%base(ie)%D(a,j,2,2) &
                                * fptr%deriv%Dvv(a,i)  
                   
                   temp1_KE = rrearth * fptr%base(ie)%Dinv(i,j,1,1) &
                              * fptr%deriv%Dvv(a,i) 

                   temp2_KE = rrearth * fptr%base(ie)%Dinv(i,j,1,2) &
                              * fptr%deriv%Dvv(a,i)
                end if

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)

                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(a,j,c,i,j,k+1) &
                              * half_rdpdn_diffu1

                ! dV(i,j,k) / dU(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(a,j,c,i,j,k+1) &
                              * half_rdpdn_diffv1

                if (c == k) then
                   ! dU(i,j,k) / dU(a,b,c)
                   Uvals(cidx) = Uvals(cidx) + massmatrix(i,j) * &
                        (-1.0d0 * temp2_Vort * fptr%base(ie)%state%v(i,j,2,k,np1) &
                        + temp1_KE * fptr%base(ie)%state%v(a,j,1,k,np1))
                   
                   ! dV(i,j,k) / dU(a,b,c)                  
                   Vvals(cidx) = Vvals(cidx) + massmatrix(i,j) * &
                        (temp2_Vort * fptr%base(ie)%state%v(i,j,1,k,np1) &
                        + temp2_KE * fptr%base(ie)%state%v(a,j,1,k,np1))
                end if

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)

                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(a,j,c,i,j,k+1) &
                              * half_rdpdn_diffu1

                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(a,j,c,i,j,k+1) &
                              * half_rdpdn_diffv1

                if (c == k) then
                   ! dU(i,j,k) / dV(a,b,c)
                   Uvals(cidx) = Uvals(cidx) + massmatrix(i,j) * &
                        (-1.0d0 * temp3_Vort * fptr%base(ie)%state%v(i,j,2,k,np1) &
                        + temp1_KE * fptr%base(ie)%state%v(a,j,2,k,np1))
                   
                   ! dV(i,j,k) / dV(a,b,c)
                   Vvals(cidx) = Vvals(cidx) + massmatrix(i,j) * &
                        (temp3_Vort * fptr%base(ie)%state%v(i,j,1,k,np1) &
                        + temp2_KE * fptr%base(ie)%state%v(a,j,2,k,np1))
                end if               

             end do

             ! a = i, b /= j, c = * (skips zeros in d*_etadot_dpdn)
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                if (c == k) then
                   temp2_Vort = -1.0d0 * temp1_Vort * fptr%base(ie)%D(i,b,1,1) &
                                 * fptr%deriv%Dvv(b,j)

                   temp3_Vort = -1.0d0 * temp1_Vort * fptr%base(ie)%D(i,b,2,1) &
                                * fptr%deriv%Dvv(b,j)
                   
                   temp1_KE = rrearth * fptr%base(ie)%Dinv(i,j,2,1) &
                              * fptr%deriv%Dvv(b,j) 

                   temp2_KE = rrearth * fptr%base(ie)%Dinv(i,j,2,2) &
                              * fptr%deriv%Dvv(b,j) 
                end if

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)

                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,b,c,i,j,k+1) &
                              * half_rdpdn_diffu1
                
                ! dV(i,j,k) / dU(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,b,c,i,j,k+1) &
                              * half_rdpdn_diffv1

                if (c == k) then
                   ! dU(i,j,k) / dU(a,b,c)
                   Uvals(cidx) = Uvals(cidx) + massmatrix(i,j) * &
                        (-1.0d0 * temp2_Vort * fptr%base(ie)%state%v(i,j,2,k,np1) &
                        + temp1_KE * fptr%base(ie)%state%v(i,b,1,k,np1)) 
                   
                   ! dV(i,j,k) / dU(a,b,c)
                   Vvals(cidx) = Vvals(cidx) + massmatrix(i,j) * &
                        (temp2_Vort * fptr%base(ie)%state%v(i,j,1,k,np1) &
                        + temp2_KE * fptr%base(ie)%state%v(i,b,1,k,np1))
                end if

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)

                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,b,c,i,j,k+1) &
                              * half_rdpdn_diffu1

                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,b,c,i,j,k+1) &
                              * half_rdpdn_diffv1

                if (c == k) then
                   ! dU(i,j,k) / dV(a,b,c)
                   Uvals(cidx) = Uvals(cidx) + massmatrix(i,j) * &
                        (-1.0d0 * temp3_Vort * fptr%base(ie)%state%v(i,j,2,k,np1) &
                        + temp1_KE * fptr%base(ie)%state%v(i,b,2,k,np1))
                   
                   ! dV(i,j,k) / dV(a,b,c)
                   Vvals(cidx) = Vvals(cidx) + massmatrix(i,j) * & 
                        (temp3_Vort * fptr%base(ie)%state%v(i,j,1,k,np1) &
                        + temp2_KE * fptr%base(ie)%state%v(i,b,2,k,np1))
                end if

             end do

             ! a = i, b = j
             ! ------------------------------------------------------------------

             temp1 = half_rdpdn(i,j,k) * etadot_dpdn(i,j,k+1) ! etadot_dpdn k+1/2
             
             if (c == k) then

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)
                
                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * & 
                     (du_etadot_dpdn(i,j,k,i,j,k+1) * half_rdpdn_diffu1 - temp1)
                
                ! dV(i,j,k) / dU(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,j,c,i,j,k+1) &
                              * half_rdpdn_diffv1
                
                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)
                
                
                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,j,c,i,j,k+1) &
                     * half_rdpdn_diffu1

                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * & 
                     (dv_etadot_dpdn(i,j,k,i,j,k+1) * half_rdpdn_diffv1 - temp1)
                
             else if (c == k+1) then
                
                ! a = i, b = j, c = k+1
                ! ---------------------------------------------------------------
                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k+1)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k+1)
                
                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * & 
                     (du_etadot_dpdn(i,j,k+1,i,j,k+1) * half_rdpdn_diffu1 + temp1)

                ! dV(i,j,k) / dU(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,j,c,i,j,k+1) &
                              * half_rdpdn_diffv1
                
                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k+1)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k+1)

                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,j,c,i,j,k+1) &
                     * half_rdpdn_diffu1
                
                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * & 
                     (dv_etadot_dpdn(i,j,k+1,i,j,k+1) * half_rdpdn_diffv1 + temp1)

             else
                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)
                
                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,j,c,i,j,k+1) &
                              * half_rdpdn_diffu1

                ! dV(i,j,k) / dU(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,j,c,i,j,k+1) &
                              * half_rdpdn_diffv1

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)               

                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,j,c,i,j,k+1) &
                              * half_rdpdn_diffu1
                
                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,j,c,i,j,k+1) &
                              * half_rdpdn_diffv1
             end if

          end do
       end do
    end do

    ! Levels 2 to nelv-1
    ! ---------------------------------------------------------------------------
    do k = 2,nlev-1
       do j = 1,np
          do i = 1,np

             ridx = ridx + 1            
             Urows(ridx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
             Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)
             
             half_rdpdn_diffu1 = half_rdpdn(i,j,k) &
                  * (fptr%base(ie)%state%v(i,j,1,k+1,np1) &
                  - fptr%base(ie)%state%v(i,j,1,k,np1))

             half_rdpdn_diffv1 = half_rdpdn(i,j,k) &
                  * (fptr%base(ie)%state%v(i,j,2,k+1,np1) &
                  - fptr%base(ie)%state%v(i,j,2,k,np1))

             half_rdpdn_diffu2 = half_rdpdn(i,j,k) &
                  * (fptr%base(ie)%state%v(i,j,1,k,np1) &
                  - fptr%base(ie)%state%v(i,j,1,k-1,np1))

             half_rdpdn_diffv2 = half_rdpdn(i,j,k) &
                  * (fptr%base(ie)%state%v(i,j,2,k,np1) &
                  - fptr%base(ie)%state%v(i,j,2,k-1,np1))

             do c = 1,nlev

                ! a /= i, b = j, c = * (skips zeros in d*_etadot_dpdn)
                ! ---------------------------------------------------------------
                do a = 1,np
                   if (a == i) cycle

                   cidx = cidx + 1              
                   Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)
                   
                   ! dU(i,j,k) / dU(a,b,c)
                   Uvals(cidx) = massmatrix(i,j) * & 
                        (du_etadot_dpdn(a,j,c,i,j,k+1) * half_rdpdn_diffu1 &
                        + du_etadot_dpdn(a,j,c,i,j,k) * half_rdpdn_diffu2)

                   ! dV(i,j,k) / dU(a,b,c) 
                   Vvals(cidx) = massmatrix(i,j) * &
                        (du_etadot_dpdn(a,j,c,i,j,k+1) * half_rdpdn_diffv1 &
                        + du_etadot_dpdn(a,j,c,i,j,k) * half_rdpdn_diffv2)

                   cidx = cidx + 1              
                   Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)

                   ! dU(i,j,k) / dV(a,b,c)
                   Uvals(cidx) = massmatrix(i,j) * & 
                        (dv_etadot_dpdn(a,j,c,i,j,k+1) * half_rdpdn_diffu1 &
                        + dv_etadot_dpdn(a,j,c,i,j,k) * half_rdpdn_diffu2)

                   ! dV(i,j,k) / dV(a,b,c)
                   Vvals(cidx) = massmatrix(i,j) * & 
                        (dv_etadot_dpdn(a,j,c,i,j,k+1) * half_rdpdn_diffv1 &
                        + dv_etadot_dpdn(a,j,c,i,j,k) * half_rdpdn_diffv2)
                end do

                ! a = i, b /= j, c = * (skips zeros in d*_etadot_dpdn)
                ! ---------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle

                   cidx = cidx + 1              
                   Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)
                   
                   ! dU(i,j,k) / dU(a,b,c)
                   Uvals(cidx) = massmatrix(i,j) * &
                        (du_etadot_dpdn(i,b,c,i,j,k+1) * half_rdpdn_diffu1 &
                        + du_etadot_dpdn(i,b,c,i,j,k) * half_rdpdn_diffu2)
                   
                   ! dV(i,j,k) / dU(a,b,c)
                   Vvals(cidx) = massmatrix(i,j) * &
                        (du_etadot_dpdn(i,b,c,i,j,k+1) * half_rdpdn_diffv1 &
                        + du_etadot_dpdn(i,b,c,i,j,k) * half_rdpdn_diffv2)

                   cidx = cidx + 1              
                   Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)

                   ! dU(i,j,k) / dV(a,b,c)
                   Uvals(cidx) = massmatrix(i,j) * &
                        (dv_etadot_dpdn(i,b,c,i,j,k+1) * half_rdpdn_diffu1 &
                        + dv_etadot_dpdn(i,b,c,i,j,k) * half_rdpdn_diffu2)
                   
                   ! dV(i,j,k) / dV(a,b,c)
                   Vvals(cidx) = massmatrix(i,j) * &
                        (dv_etadot_dpdn(i,b,c,i,j,k+1) * half_rdpdn_diffv1 &
                        + dv_etadot_dpdn(i,b,c,i,j,k) * half_rdpdn_diffv2)
                end do

                ! a = i, b = j
                ! ---------------------------------------------------------------
                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)

                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * &
                     (dv_etadot_dpdn(i,j,c,i,j,k+1) * half_rdpdn_diffu1 &
                     + dv_etadot_dpdn(i,j,c,i,j,k) * half_rdpdn_diffu2)
                                
                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * &
                     (du_etadot_dpdn(i,j,c,i,j,k+1) * half_rdpdn_diffv1 &
                     + du_etadot_dpdn(i,j,c,i,j,k) * half_rdpdn_diffv2)
                
                ! check if this is a special case (product rule)
                if ( c == k-1 .or. c == k .or. c == k+1 ) then
                   cycle
                else
                   cidx = cidx + 1              
                   Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)
                   
                   ! dU(i,j,k) / dU(a,b,c)
                   Uvals(cidx) = massmatrix(i,j) * &
                        (du_etadot_dpdn(i,j,c,i,j,k+1) * half_rdpdn_diffu1 &
                        + du_etadot_dpdn(i,j,c,i,j,k) * half_rdpdn_diffu2)

                   ! dV(i,j,k) / dV(a,b,c)
                   Vvals(cidx) = massmatrix(i,j) * &
                        (dv_etadot_dpdn(i,j,c,i,j,k+1) * half_rdpdn_diffv1 &
                        + dv_etadot_dpdn(i,j,c,i,j,k) * half_rdpdn_diffv2)
                end if
             end do

             ! special cases
             temp1 = half_rdpdn(i,j,k) * etadot_dpdn(i,j,k+1) ! etadot_dpdn l+1/2
             temp2 = half_rdpdn(i,j,k) * etadot_dpdn(i,j,k)   ! etadot_dpdn l-1/2

             ! a = i, b = j, c = k-1
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k-1)          
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k-1)

             ! dU(i,j,k) / dU(a,b,c)
             Uvals(cidx) = massmatrix(i,j) * &
                  (du_etadot_dpdn(i,j,k-1,i,j,k+1) * half_rdpdn_diffu1 &
                  + du_etadot_dpdn(i,j,k-1,i,j,k) * half_rdpdn_diffu2 - temp2)
             
             ! dV(i,j,k) / dV(a,b,c)
             Vvals(cidx) = massmatrix(i,j) * &
                  (dv_etadot_dpdn(i,j,k-1,i,j,k+1) * half_rdpdn_diffv1 &
                  + dv_etadot_dpdn(i,j,k-1,i,j,k) * half_rdpdn_diffv2 - temp2)

             ! a = i, b = j, c = k
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)
             
             ! dU(i,j,k) / du(a,b,c)
             Uvals(cidx) = massmatrix(i,j) * &
                  (du_etadot_dpdn(i,j,k,i,j,k+1) * half_rdpdn_diffu1 &
                  + du_etadot_dpdn(i,j,k,i,j,k) * half_rdpdn_diffu2 &
                  - temp1 + temp2)

             ! dV(i,j,k) / dV(a,b,c)
             Vvals(cidx) = massmatrix(i,j) * &
                  (dv_etadot_dpdn(i,j,k,i,j,k+1) * half_rdpdn_diffv1 &
                  + dv_etadot_dpdn(i,j,k,i,j,k) * half_rdpdn_diffv2 &
                  - temp1 + temp2)

             ! a = i, b = j, c = k+1
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k+1)
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k+1)
             
             ! dU(i,j,k) / dU(a,b,c)
             Uvals(cidx) = massmatrix(i,j) * &
                  (du_etadot_dpdn(i,j,k+1,i,j,k+1) * half_rdpdn_diffu1 &
                  + du_etadot_dpdn(i,j,k+1,i,j,k) * half_rdpdn_diffu2 + temp1)
             
             ! dV(i,j,k) / dV(a,b,c)
             Vvals(cidx) = massmatrix(i,j) * &
                  (dv_etadot_dpdn(i,j,k+1,i,j,k+1) * half_rdpdn_diffv1 &
                  + dv_etadot_dpdn(i,j,k+1,i,j,k) * half_rdpdn_diffv2 + temp1)
          end do
       end do
    end do

    ! Level nlev
    ! ---------------------------------------------------------------------------
    k = nlev
    do j = 1,np
       do i = 1,np

          ridx = ridx + 1            
          Urows(ridx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
          Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)

          half_rdpdn_diffu2 = half_rdpdn(i,j,k) &
               * (fptr%base(ie)%state%v(i,j,1,k,np1) &
               - fptr%base(ie)%state%v(i,j,1,k-1,np1))

          half_rdpdn_diffv2 = half_rdpdn(i,j,k) &
               * (fptr%base(ie)%state%v(i,j,2,k,np1) &
               - fptr%base(ie)%state%v(i,j,2,k-1,np1))

          do c = 1,nlev

             ! a /= i, b = j, c = * (skip zeros in d*_etadot_dpdn)
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)
                
                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(a,j,c,i,j,k) &
                              * half_rdpdn_diffu2

                ! dV(i,j,k) / dU(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(a,j,c,i,j,k) &
                              * half_rdpdn_diffv2

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)

                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(a,j,c,i,j,k) &
                              * half_rdpdn_diffu2

                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(a,j,c,i,j,k) &
                              * half_rdpdn_diffv2
             end do

             ! a = i, b /= j, c = * (skips zeros in d*_etadot_dpdn)
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)

                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,b,c,i,j,k) &
                              * half_rdpdn_diffu2

                ! dV(i,j,k) / dU(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,b,c,i,j,k) &
                              * half_rdpdn_diffv2

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)

                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,b,c,i,j,k) &
                              * half_rdpdn_diffu2

                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,b,c,i,j,k) &
                              * half_rdpdn_diffv2
             end do

             ! a = i, b = j, c = *
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)

             ! dU(i,j,k) / dV(a,b,c)
             Uvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,j,c,i,j,k) &
                           * half_rdpdn_diffu2

             ! dV(i,j,k) / dU(a,b,c)
             Vvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,j,c,i,j,k) &
                           * half_rdpdn_diffv2

             ! check if this is a special case (product rule) 
             if ( c == k-1 .or. c == k) then
                cycle
             else
                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)

                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,j,c,i,j,k) &
                              * half_rdpdn_diffu2
                
                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,j,c,i,j,k) &
                              * half_rdpdn_diffv2
             end if
          end do

          ! special cases
          temp2 = half_rdpdn(i,j,k) * etadot_dpdn(i,j,k) ! etadot_dpdn l-1/2

          ! a = i, b = j, c = k-1
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k-1)
          Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k-1)
          
          ! dU(i,j,k) / dU(a,b,c)
          Uvals(cidx) = massmatrix(i,j) * &
               (du_etadot_dpdn(i,j,k-1,i,j,k) * half_rdpdn_diffu2 - temp2)

          ! dV(i,j,k) / dV(a,b,c)
          Vvals(cidx) = massmatrix(i,j) * &
               (dv_etadot_dpdn(i,j,k-1,i,j,k) * half_rdpdn_diffv2 - temp2)

          ! a = i, b = j, c = k
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
          Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)
          
          ! dU(i,j,k) / dU(a,b,c)
          Uvals(cidx) = massmatrix(i,j) * &
               (du_etadot_dpdn(i,j,k,i,j,k) * half_rdpdn_diffu2 + temp2)

          ! dV(i,j,k) / dV(a,b,c)
          Vvals(cidx) = massmatrix(i,j) * &
               (dv_etadot_dpdn(i,j,k,i,j,k) * half_rdpdn_diffv2 + temp2)
       end do
    end do
  
    call t_stopf('Jacobian_ABlock')    
    
  end subroutine Jacobian_Ablock


  subroutine Jac_Ablock(ie, nnz_max, row_nnz, Urows, Ucols, Uvals, &
       Vrows, Vcols, Vvals, c_ptr_to_object) bind(C,name="Jac_Ablock")
    !----------------------------------------------------------------------------
    ! Computes the A block of the Jacobian matrix contributions from:
    !     - Time Derivative
    !     - Vorticity
    !     - Coriolis
    !     - Kinetic Energy
    !
    ! Contributions from vertical advection are computed in a separate subroutine
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use control_mod, only           : tstep_type
    use dimensions_mod, only        : np, nlev
    use physical_constants, only    : rrearth
    use prim_derived_type_mod, only : derived_type
    use derivative_mod, only        : vorticity_sphere

    implicit none

    !---------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    ! row indices 
    integer(c_int), dimension(np*np*nlev), intent(out) :: Urows
    integer(c_int), dimension(np*np*nlev), intent(out) :: Vrows

    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Ucols
    integer(c_int), dimension(nnz_max), intent(out) :: Vcols

    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Uvals   
    real(c_double), dimension(nnz_max), intent(out) :: Vvals
    
    type(c_ptr) :: c_ptr_to_object
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()
 
    real(kind=real_kind) :: dti                                ! 1/dt
    real(kind=real_kind) :: temp1_Vort, temp2_Vort, temp3_Vort ! temp vort vars
    real(kind=real_kind) :: temp1_KE, temp2_KE                 ! temp KE vars

    real(kind=real_kind), dimension(np,np,nlev) :: zeta ! vorticity

    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    integer :: np1         ! time level
    integer :: a, b        ! derivative index (c not needed)
    integer :: i, j, k     ! equation index
    integer :: ridx, cidx  ! global row/column index 
    !----------------------------------------------------------------------------

    call t_startf('Jac_Ablock')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*) "Computing Ablock"
#endif

    ! initialize arrays
    Urows = 0
    Ucols = 0
    Uvals = 0.0d0

    Vrows = 0
    Vcols = 0
    Vvals = 0.0d0
    
    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = 2*(np+np-1)

    ! time level where current iteration is stored
    np1 = fptr%tl%np1

    ! inverse of time step size
    dti = 1.0d0 / fptr%dt

    if (tstep_type==12) then
       dti = (3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2 on step size
    endif

    ! vorticity 
    do k = 1,nlev
       zeta(:,:,k) = vorticity_sphere(fptr%base(ie)%state%v(:,:,:,k,np1), fptr%deriv, fptr%base(ie))
    end do

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp

    ! ---------------------------------------------------------------------------
    ! compute Jacobian entires 
    ! ---------------------------------------------------------------------------
    do k = 1,nlev
       do j = 1,np
          do i = 1,np

             ! set row index for U and V equations
             ridx = ridx + 1
             Urows(ridx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
             Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)

             temp1_Vort = rrearth * fptr%base(ie)%rmetdet(i,j)

             ! a /= i, b = j, c = k
             ! ---------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                temp2_Vort = temp1_Vort * fptr%base(ie)%D(a,j,1,2) * fptr%deriv%Dvv(a,i)  
                temp3_Vort = temp1_Vort * fptr%base(ie)%D(a,j,2,2) * fptr%deriv%Dvv(a,i)  

                temp1_KE = rrearth * fptr%base(ie)%Dinv(i,j,1,1) * fptr%deriv%Dvv(a,i) 
                temp2_KE = rrearth * fptr%base(ie)%Dinv(i,j,1,2) * fptr%deriv%Dvv(a,i) 

                cidx = cidx + 1
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(k)
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(k)

                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * &
                     (-1.0d0 * temp2_Vort * fptr%base(ie)%state%v(i,j,2,k,np1) &
                     + temp1_KE * fptr%base(ie)%state%v(a,j,1,k,np1))

                ! dV(i,j,k) / dU(a,b,c)                  
                Vvals(cidx) = massmatrix(i,j) * &
                     (temp2_Vort * fptr%base(ie)%state%v(i,j,1,k,np1) &
                     + temp2_KE * fptr%base(ie)%state%v(a,j,1,k,np1))

                cidx = cidx + 1
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(k)
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(k)

                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * &
                     (-1.0d0 * temp3_Vort * fptr%base(ie)%state%v(i,j,2,k,np1) &
                     + temp1_KE * fptr%base(ie)%state%v(a,j,2,k,np1))

                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * &
                     (temp3_Vort * fptr%base(ie)%state%v(i,j,1,k,np1) &
                     + temp2_KE * fptr%base(ie)%state%v(a,j,2,k,np1))
             end do


             ! a = i, b /= j, c = k
             ! ---------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                temp2_Vort = -1.0d0 * temp1_Vort * fptr%base(ie)%D(i,b,1,1) * fptr%deriv%Dvv(b,j)
                temp3_Vort = -1.0d0 * temp1_Vort * fptr%base(ie)%D(i,b,2,1) * fptr%deriv%Dvv(b,j)

                temp1_KE = rrearth * fptr%base(ie)%Dinv(i,j,2,1) * fptr%deriv%Dvv(b,j) 
                temp2_KE = rrearth * fptr%base(ie)%Dinv(i,j,2,2) * fptr%deriv%Dvv(b,j) 

                cidx = cidx + 1
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(k)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(k)

                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * &
                     (-1.0d0 * temp2_Vort * fptr%base(ie)%state%v(i,j,2,k,np1) &
                     + temp1_KE * fptr%base(ie)%state%v(i,b,1,k,np1)) 

                ! dV(i,j,k) / dU(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * &
                     (temp2_Vort * fptr%base(ie)%state%v(i,j,1,k,np1) &
                     + temp2_KE * fptr%base(ie)%state%v(i,b,1,k,np1))

                cidx = cidx + 1
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(k)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(k)

                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * &
                     (-1.0d0 * temp3_Vort * fptr%base(ie)%state%v(i,j,2,k,np1) &
                     + temp1_KE * fptr%base(ie)%state%v(i,b,2,k,np1))

                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * & 
                     (temp3_Vort * fptr%base(ie)%state%v(i,j,1,k,np1) &
                     + temp2_KE * fptr%base(ie)%state%v(i,b,2,k,np1))
             end do


             ! a = i, b = j, c = k
             ! ---------------------------------------------------------------
             temp1_Vort = rrearth * fptr%base(ie)%rmetdet(i,j) &
                  * (fptr%base(ie)%D(i,j,1,2) * fptr%deriv%Dvv(i,i) &
                  -  fptr%base(ie)%D(i,j,1,1) * fptr%deriv%Dvv(j,j))

             temp2_Vort = rrearth * fptr%base(ie)%rmetdet(i,j) & 
                  * (fptr%base(ie)%D(i,j,2,2) * fptr%deriv%Dvv(i,i) &
                  -  fptr%base(ie)%D(i,j,2,1) * fptr%deriv%Dvv(j,j))

             temp1_KE = rrearth &
                  * (fptr%base(ie)%Dinv(i,j,1,1) * fptr%deriv%Dvv(i,i) &
                  +  fptr%base(ie)%Dinv(i,j,2,1) * fptr%deriv%Dvv(j,j))

             temp2_KE = rrearth &
                  * (fptr%base(ie)%Dinv(i,j,1,2) * fptr%deriv%Dvv(i,i) &
                  +  fptr%base(ie)%Dinv(i,j,2,2) * fptr%deriv%Dvv(j,j))

             cidx = cidx + 1
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)

             ! dU(i,j,k) / dU(a,b,c)
             Uvals(cidx) = massmatrix(i,j) * & 
                  (dti - temp1_Vort * fptr%base(ie)%state%v(i,j,2,k,np1) &
                  + temp1_KE * fptr%base(ie)%state%v(i,j,1,k,np1)) 

             ! dV(i,j,k) / dU(a,b,c)
             Vvals(cidx) = massmatrix(i,j) * & 
                  (fptr%base(ie)%fcor(i,j) &
                  + temp1_Vort * fptr%base(ie)%state%v(i,j,1,k,np1) &
                  + zeta(i,j,k) &
                  + temp2_KE * fptr%base(ie)%state%v(i,j,1,k,np1))

             cidx = cidx + 1
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)

             ! dU(i,j,k) / dV(a,b,c)
             Uvals(cidx) = massmatrix(i,j) * & 
                  (-1.0d0 * fptr%base(ie)%fcor(i,j) &
                  - temp2_Vort * fptr%base(ie)%state%v(i,j,2,k,np1) &  
                  - zeta(i,j,k) &
                  + temp1_KE * fptr%base(ie)%state%v(i,j,2,k,np1))

             ! dV(i,j,k) / dV(a,b,c)
             Vvals(cidx) = massmatrix(i,j) * & 
                  (dti + temp2_Vort * fptr%base(ie)%state%v(i,j,1,k,np1) &
                  + temp2_KE * fptr%base(ie)%state%v(i,j,2,k,np1))
          end do
       end do
    end do

    call t_stopf('Jac_Ablock')

  end subroutine Jac_Ablock



  subroutine Jac_Ablock_VertAdv(ie, nnz_max, row_nnz, Urows, Ucols, Uvals, &
       Vrows, Vcols, Vvals, c_ptr_to_object) bind(C,name="Jac_Ablock_VertAdv")
    !----------------------------------------------------------------------------
    ! Computes the A block of the Jacobian matrix contributions from:
    !     - Vertical Advection
    !
    ! Contributions from Vorticity, Coriolis, Kinetic Energy
    ! are computed in a separate subroutine
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use derivative_mod, only        : divergence_sphere
    
    implicit none

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    ! row indices 
    integer(c_int), dimension(np*np*nlev), intent(out) :: Urows
    integer(c_int), dimension(np*np*nlev), intent(out) :: Vrows

    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Ucols
    integer(c_int), dimension(nnz_max), intent(out) :: Vcols

    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Uvals   
    real(c_double), dimension(nnz_max), intent(out) :: Vvals      

    type(c_ptr) :: c_ptr_to_object ! HOMME data
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,nlev,np,np,nlev+1) :: du_etadot_dpdn ! u derivative of etadot_dpdn
    real(kind=real_kind), dimension(np,np,nlev,np,np,nlev+1) :: dv_etadot_dpdn ! v derivative of etadot_dpdn

    real(kind=real_kind), dimension(np,np,nlev+1) :: etadot_dpdn

    real(kind=real_kind), dimension(np,np,nlev) :: dpdn         ! vertical derivative of pressure
    real(kind=real_kind), dimension(np,np,nlev) :: half_rdpdn   ! 0.5 * 1/dpdn
    real(kind=real_kind), dimension(np,np,nlev) :: div_dpdn_vel ! divergence(dpdn * velocity) 

    real(kind=real_kind), dimension(np,np,2) :: dpdn_vel ! dpdn * velocity  

    real(kind=real_kind), dimension(np,np) :: sdot_sum
    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    real(kind=real_kind) :: half_rdpdn_diffu1, half_rdpdn_diffu2
    real(kind=real_kind) :: half_rdpdn_diffv1, half_rdpdn_diffv2
    real(kind=real_kind) :: temp1, temp2
    
    integer :: np1        ! current time level 
    integer :: i, j, k    ! local equation index
    integer :: a, b, c    ! local derivative index
    integer :: ridx, cidx ! index for output U and V output arrays
    !----------------------------------------------------------------------------
    
    call t_startf('Jac_ABlock_VertAdv')    

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*) "Computing Ablock_VertAdv"
#endif

    ! initialize arrays
    Urows = 0
    Ucols = 0
    Uvals = 0.0d0

    Vrows = 0
    Vcols = 0
    Vvals = 0.0d0

    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = 2*(np+np-1)*nlev

    ! time level where current iteration is stored
    np1 = fptr%tl%np1
       
    ! vertical derivative of pressure
    do k = 1,nlev
       dpdn(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
            - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))
    end do
    
    half_rdpdn = 0.5d0 / dpdn
    
    ! divergence(dpdn * velocity)
    do k = 1,nlev
       do j = 1,np
          do i = 1,np             
             dpdn_vel(i,j,1) = dpdn(i,j,k) * fptr%base(ie)%state%v(i,j,1,k,np1)
             dpdn_vel(i,j,2) = dpdn(i,j,k) * fptr%base(ie)%state%v(i,j,2,k,np1)
          end do
       end do
       
       div_dpdn_vel(:,:,k) = divergence_sphere(dpdn_vel, fptr%deriv, fptr%base(ie))
    end do
    
    ! compute etadot_dpdn
    sdot_sum = 0.0d0

    do k=1,nlev
       sdot_sum(:,:)        = sdot_sum(:,:) + div_dpdn_vel(:,:,k)
       etadot_dpdn(:,:,k+1) = sdot_sum(:,:)
    end do

    do k=1,nlev-1
       etadot_dpdn(:,:,k+1) = fptr%hvcoord%hybi(k+1)*sdot_sum(:,:) &
            - etadot_dpdn(:,:,k+1)
    end do

    etadot_dpdn(:,:,1)      = 0.0d0
    etadot_dpdn(:,:,nlev+1) = 0.0d0

    ! velocity derivatives of etadot_dpdn
    call DerivVel_etadot_dpdn(ie, fptr, du_etadot_dpdn, dv_etadot_dpdn)

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp

    ! ---------------------------------------------------------------------------
    ! compute Jacobian entires 
    ! ---------------------------------------------------------------------------

    ! Level 1
    ! ---------------------------------------------------------------------------
    k = 1
    do j = 1,np
       do i = 1,np

          ! get global row index for U and V equation
          ridx = ridx + 1            
          Urows(ridx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
          Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)
         
          half_rdpdn_diffu1 = half_rdpdn(i,j,k) * &
               (fptr%base(ie)%state%v(i,j,1,k+1,np1) &
               - fptr%base(ie)%state%v(i,j,1,k,np1))

          half_rdpdn_diffv1 = half_rdpdn(i,j,k) &
               * (fptr%base(ie)%state%v(i,j,2,k+1,np1) &
               - fptr%base(ie)%state%v(i,j,2,k,np1))

          do c = 1,nlev

             ! a /= i, b = j, c = * (skips zeros in d*_etadot_dpdn)
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)

                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(a,j,c,i,j,k+1) &
                              * half_rdpdn_diffu1

                ! dV(i,j,k) / dU(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(a,j,c,i,j,k+1) &
                              * half_rdpdn_diffv1

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)

                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(a,j,c,i,j,k+1) &
                              * half_rdpdn_diffu1

                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(a,j,c,i,j,k+1) &
                              * half_rdpdn_diffv1

             end do

             ! a = i, b /= j, c = * (skips zeros in d*_etadot_dpdn)
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)

                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,b,c,i,j,k+1) &
                              * half_rdpdn_diffu1
                
                ! dV(i,j,k) / dU(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,b,c,i,j,k+1) &
                              * half_rdpdn_diffv1

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)

                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,b,c,i,j,k+1) &
                              * half_rdpdn_diffu1

                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,b,c,i,j,k+1) &
                              * half_rdpdn_diffv1

             end do

             ! a = i, b = j, c /= k or c /= k-1
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)

             ! dU(i,j,k) / dV(a,b,c)
             Uvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,j,c,i,j,k+1) &
                           * half_rdpdn_diffu1
                
             ! dV(i,j,k) / dU(a,b,c)
             Vvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,j,c,i,j,k+1) &
                           * half_rdpdn_diffv1

             ! check if this is a special case (product rule)
             if (c == k .or. c == k+1) then                                  
                cycle
             else
                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)
                
                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,j,c,i,j,k+1) &
                              * half_rdpdn_diffu1
                
                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,j,c,i,j,k+1) &
                              * half_rdpdn_diffv1
             end if

          end do

          ! special cases
          temp1 = half_rdpdn(i,j,k) * etadot_dpdn(i,j,k+1) ! etadot_dpdn k+1/2

          ! a = i, b = j, c = k
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
          Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)

          ! dU(i,j,k) / dU(a,b,c)
          Uvals(cidx) = massmatrix(i,j) * & 
               (du_etadot_dpdn(i,j,k,i,j,k+1) * half_rdpdn_diffu1 - temp1)
          
          ! dV(i,j,k) / dV(a,b,c)
          Vvals(cidx) = massmatrix(i,j) * & 
               (dv_etadot_dpdn(i,j,k,i,j,k+1) * half_rdpdn_diffv1 - temp1)

          ! a = i, b = j, c = k+1
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k+1)
          Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k+1)
          
          ! dU(i,j,k) / dU(a,b,c)
          Uvals(cidx) = massmatrix(i,j) * & 
               (du_etadot_dpdn(i,j,k+1,i,j,k+1) * half_rdpdn_diffu1 + temp1)

          ! dV(i,j,k) / dV(a,b,c)
          Vvals(cidx) = massmatrix(i,j) * & 
               (dv_etadot_dpdn(i,j,k+1,i,j,k+1) * half_rdpdn_diffv1 + temp1)

       end do
    end do

    ! Levels 2 to nelv-1
    ! ---------------------------------------------------------------------------
    do k = 2,nlev-1
       do j = 1,np
          do i = 1,np

             ridx = ridx + 1            
             Urows(ridx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
             Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)
             
             half_rdpdn_diffu1 = half_rdpdn(i,j,k) &
                  * (fptr%base(ie)%state%v(i,j,1,k+1,np1) &
                  - fptr%base(ie)%state%v(i,j,1,k,np1))

             half_rdpdn_diffv1 = half_rdpdn(i,j,k) &
                  * (fptr%base(ie)%state%v(i,j,2,k+1,np1) &
                  - fptr%base(ie)%state%v(i,j,2,k,np1))

             half_rdpdn_diffu2 = half_rdpdn(i,j,k) &
                  * (fptr%base(ie)%state%v(i,j,1,k,np1) &
                  - fptr%base(ie)%state%v(i,j,1,k-1,np1))

             half_rdpdn_diffv2 = half_rdpdn(i,j,k) &
                  * (fptr%base(ie)%state%v(i,j,2,k,np1) &
                  - fptr%base(ie)%state%v(i,j,2,k-1,np1))

             do c = 1,nlev

                ! a /= i, b = j, c = * (skips zeros in d*_etadot_dpdn)
                ! ---------------------------------------------------------------
                do a = 1,np
                   if (a == i) cycle

                   cidx = cidx + 1              
                   Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)
                   
                   ! dU(i,j,k) / dU(a,b,c)
                   Uvals(cidx) = massmatrix(i,j) * & 
                        (du_etadot_dpdn(a,j,c,i,j,k+1) * half_rdpdn_diffu1 &
                        + du_etadot_dpdn(a,j,c,i,j,k) * half_rdpdn_diffu2)

                   ! dV(i,j,k) / dU(a,b,c) 
                   Vvals(cidx) = massmatrix(i,j) * &
                        (du_etadot_dpdn(a,j,c,i,j,k+1) * half_rdpdn_diffv1 &
                        + du_etadot_dpdn(a,j,c,i,j,k) * half_rdpdn_diffv2)

                   cidx = cidx + 1              
                   Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)

                   ! dU(i,j,k) / dV(a,b,c)
                   Uvals(cidx) = massmatrix(i,j) * & 
                        (dv_etadot_dpdn(a,j,c,i,j,k+1) * half_rdpdn_diffu1 &
                        + dv_etadot_dpdn(a,j,c,i,j,k) * half_rdpdn_diffu2)

                   ! dV(i,j,k) / dV(a,b,c)
                   Vvals(cidx) = massmatrix(i,j) * & 
                        (dv_etadot_dpdn(a,j,c,i,j,k+1) * half_rdpdn_diffv1 &
                        + dv_etadot_dpdn(a,j,c,i,j,k) * half_rdpdn_diffv2)
                end do

                ! a = i, b /= j, c = * (skips zeros in d*_etadot_dpdn)
                ! ---------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle

                   cidx = cidx + 1              
                   Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)
                   
                   ! dU(i,j,k) / dU(a,b,c)
                   Uvals(cidx) = massmatrix(i,j) * &
                        (du_etadot_dpdn(i,b,c,i,j,k+1) * half_rdpdn_diffu1 &
                        + du_etadot_dpdn(i,b,c,i,j,k) * half_rdpdn_diffu2)
                   
                   ! dV(i,j,k) / dU(a,b,c)
                   Vvals(cidx) = massmatrix(i,j) * &
                        (du_etadot_dpdn(i,b,c,i,j,k+1) * half_rdpdn_diffv1 &
                        + du_etadot_dpdn(i,b,c,i,j,k) * half_rdpdn_diffv2)

                   cidx = cidx + 1              
                   Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)

                   ! dU(i,j,k) / dV(a,b,c)
                   Uvals(cidx) = massmatrix(i,j) * &
                        (dv_etadot_dpdn(i,b,c,i,j,k+1) * half_rdpdn_diffu1 &
                        + dv_etadot_dpdn(i,b,c,i,j,k) * half_rdpdn_diffu2)
                   
                   ! dV(i,j,k) / dV(a,b,c)
                   Vvals(cidx) = massmatrix(i,j) * &
                        (dv_etadot_dpdn(i,b,c,i,j,k+1) * half_rdpdn_diffv1 &
                        + dv_etadot_dpdn(i,b,c,i,j,k) * half_rdpdn_diffv2)
                end do

                ! a = i, b = j
                ! ---------------------------------------------------------------
                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)

                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * &
                     (dv_etadot_dpdn(i,j,c,i,j,k+1) * half_rdpdn_diffu1 &
                     + dv_etadot_dpdn(i,j,c,i,j,k) * half_rdpdn_diffu2)
                                
                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * &
                     (du_etadot_dpdn(i,j,c,i,j,k+1) * half_rdpdn_diffv1 &
                     + du_etadot_dpdn(i,j,c,i,j,k) * half_rdpdn_diffv2)
                
                ! check if this is a special case (product rule)
                if ( c == k-1 .or. c == k .or. c == k+1 ) then
                   cycle
                else
                   cidx = cidx + 1              
                   Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)
                   
                   ! dU(i,j,k) / dU(a,b,c)
                   Uvals(cidx) = massmatrix(i,j) * &
                        (du_etadot_dpdn(i,j,c,i,j,k+1) * half_rdpdn_diffu1 &
                        + du_etadot_dpdn(i,j,c,i,j,k) * half_rdpdn_diffu2)

                   ! dV(i,j,k) / dV(a,b,c)
                   Vvals(cidx) = massmatrix(i,j) * &
                        (dv_etadot_dpdn(i,j,c,i,j,k+1) * half_rdpdn_diffv1 &
                        + dv_etadot_dpdn(i,j,c,i,j,k) * half_rdpdn_diffv2)
                end if
             end do

             ! special cases
             temp1 = half_rdpdn(i,j,k) * etadot_dpdn(i,j,k+1) ! etadot_dpdn l+1/2
             temp2 = half_rdpdn(i,j,k) * etadot_dpdn(i,j,k)   ! etadot_dpdn l-1/2

             ! a = i, b = j, c = k-1
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k-1)          
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k-1)

             ! dU(i,j,k) / dU(a,b,c)
             Uvals(cidx) = massmatrix(i,j) * &
                  (du_etadot_dpdn(i,j,k-1,i,j,k+1) * half_rdpdn_diffu1 &
                  + du_etadot_dpdn(i,j,k-1,i,j,k) * half_rdpdn_diffu2 - temp2)
             
             ! dV(i,j,k) / dV(a,b,c)
             Vvals(cidx) = massmatrix(i,j) * &
                  (dv_etadot_dpdn(i,j,k-1,i,j,k+1) * half_rdpdn_diffv1 &
                  + dv_etadot_dpdn(i,j,k-1,i,j,k) * half_rdpdn_diffv2 - temp2)

             ! a = i, b = j, c = k
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)
             
             ! dU(i,j,k) / du(a,b,c)
             Uvals(cidx) = massmatrix(i,j) * &
                  (du_etadot_dpdn(i,j,k,i,j,k+1) * half_rdpdn_diffu1 &
                  + du_etadot_dpdn(i,j,k,i,j,k) * half_rdpdn_diffu2 &
                  - temp1 + temp2)

             ! dV(i,j,k) / dV(a,b,c)
             Vvals(cidx) = massmatrix(i,j) * &
                  (dv_etadot_dpdn(i,j,k,i,j,k+1) * half_rdpdn_diffv1 &
                  + dv_etadot_dpdn(i,j,k,i,j,k) * half_rdpdn_diffv2 &
                  - temp1 + temp2)

             ! a = i, b = j, c = k+1
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k+1)
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k+1)
             
             ! dU(i,j,k) / dU(a,b,c)
             Uvals(cidx) = massmatrix(i,j) * &
                  (du_etadot_dpdn(i,j,k+1,i,j,k+1) * half_rdpdn_diffu1 &
                  + du_etadot_dpdn(i,j,k+1,i,j,k) * half_rdpdn_diffu2 + temp1)
             
             ! dV(i,j,k) / dV(a,b,c)
             Vvals(cidx) = massmatrix(i,j) * &
                  (dv_etadot_dpdn(i,j,k+1,i,j,k+1) * half_rdpdn_diffv1 &
                  + dv_etadot_dpdn(i,j,k+1,i,j,k) * half_rdpdn_diffv2 + temp1)
          end do
       end do
    end do

    ! Level nlev
    ! ---------------------------------------------------------------------------
    k = nlev
    do j = 1,np
       do i = 1,np

          ridx = ridx + 1            
          Urows(ridx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
          Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)

          half_rdpdn_diffu2 = half_rdpdn(i,j,k) &
               * (fptr%base(ie)%state%v(i,j,1,k,np1) &
               - fptr%base(ie)%state%v(i,j,1,k-1,np1))

          half_rdpdn_diffv2 = half_rdpdn(i,j,k) &
               * (fptr%base(ie)%state%v(i,j,2,k,np1) &
               - fptr%base(ie)%state%v(i,j,2,k-1,np1))

          do c = 1,nlev

             ! a /= i, b = j, c = * (skip zeros in d*_etadot_dpdn)
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)
                
                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(a,j,c,i,j,k) &
                              * half_rdpdn_diffu2

                ! dV(i,j,k) / dU(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(a,j,c,i,j,k) &
                              * half_rdpdn_diffv2

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)

                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(a,j,c,i,j,k) &
                              * half_rdpdn_diffu2

                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(a,j,c,i,j,k) &
                              * half_rdpdn_diffv2
             end do

             ! a = i, b /= j, c = * (skips zeros in d*_etadot_dpdn)
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)

                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,b,c,i,j,k) &
                              * half_rdpdn_diffu2

                ! dV(i,j,k) / dU(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,b,c,i,j,k) &
                              * half_rdpdn_diffv2

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)

                ! dU(i,j,k) / dV(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,b,c,i,j,k) &
                              * half_rdpdn_diffu2

                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,b,c,i,j,k) &
                              * half_rdpdn_diffv2
             end do

             ! a = i, b = j, c = *
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)

             ! dU(i,j,k) / dV(a,b,c)
             Uvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,j,c,i,j,k) &
                           * half_rdpdn_diffu2

             ! dV(i,j,k) / dU(a,b,c)
             Vvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,j,c,i,j,k) &
                           * half_rdpdn_diffv2

             ! check if this is a special case (product rule) 
             if ( c == k-1 .or. c == k) then
                cycle
             else
                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)

                ! dU(i,j,k) / dU(a,b,c)
                Uvals(cidx) = massmatrix(i,j) * du_etadot_dpdn(i,j,c,i,j,k) &
                              * half_rdpdn_diffu2
                
                ! dV(i,j,k) / dV(a,b,c)
                Vvals(cidx) = massmatrix(i,j) * dv_etadot_dpdn(i,j,c,i,j,k) &
                              * half_rdpdn_diffv2
             end if
          end do

          ! special cases
          temp2 = half_rdpdn(i,j,k) * etadot_dpdn(i,j,k) ! etadot_dpdn l-1/2

          ! a = i, b = j, c = k-1
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k-1)
          Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k-1)
          
          ! dU(i,j,k) / dU(a,b,c)
          Uvals(cidx) = massmatrix(i,j) * &
               (du_etadot_dpdn(i,j,k-1,i,j,k) * half_rdpdn_diffu2 - temp2)

          ! dV(i,j,k) / dV(a,b,c)
          Vvals(cidx) = massmatrix(i,j) * &
               (dv_etadot_dpdn(i,j,k-1,i,j,k) * half_rdpdn_diffv2 - temp2)

          ! a = i, b = j, c = k
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
          Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)
          
          ! dU(i,j,k) / dU(a,b,c)
          Uvals(cidx) = massmatrix(i,j) * &
               (du_etadot_dpdn(i,j,k,i,j,k) * half_rdpdn_diffu2 + temp2)

          ! dV(i,j,k) / dV(a,b,c)
          Vvals(cidx) = massmatrix(i,j) * &
               (dv_etadot_dpdn(i,j,k,i,j,k) * half_rdpdn_diffv2 + temp2)
       end do
    end do
  
    call t_stopf('Jac_ABlock_VertAdv')    
    
  end subroutine Jac_Ablock_VertAdv



  ! =============================================================================
  ! B Block of the Jacobian
  ! =============================================================================



  subroutine Jac_Bblock_GP(ie, nnz_max, row_nnz, Urows, Ucols, Uvals, &
       Vrows, Vcols, Vvals, c_ptr_to_object) bind(C,name="Jac_Bblock_GP")
    !----------------------------------------------------------------------------
    ! Computes the B block of the Jacobian matrix contributions from:
    !     - Geopotential
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth, rgas
    use physics_mod, only           : Virtual_Temperature
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    ! row indices 
    integer(c_int), dimension(np*np*nlev), intent(out) :: Urows
    integer(c_int), dimension(np*np*nlev), intent(out) :: Vrows

    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Ucols
    integer(c_int), dimension(nnz_max), intent(out) :: Vcols

    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Uvals   
    real(c_double), dimension(nnz_max), intent(out) :: Vvals
    
    type(c_ptr) :: c_ptr_to_object
    ! ---------------------------------------------------------------------------
    
    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,nlev) :: T_v  ! virtual temperature

    real(kind=real_kind), dimension(np,np,nlev) :: dpdn ! vertical derivative of pressure
    real(kind=real_kind), dimension(np,np,nlev) :: p    ! pressure
    real(kind=real_kind), dimension(np,np,nlev) :: p2   ! pressure squared
    real(kind=real_kind), dimension(np,np,nlev) :: Qt   !

    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    real(kind=real_kind), dimension(nlev) :: mBm, dBi


    real(kind=real_kind) :: temp, temp11, temp12, temp21, temp22

    integer :: np1         ! time level
    integer :: qn0         ! dry or moist flag
    integer :: m           ! loop index for vertical integrals
    integer :: a, b        ! derivative index
    integer :: i, j, k     ! equation index
    integer :: ridx, cidx  ! global row/column index 
    !----------------------------------------------------------------------------

    call t_startf('Jac_Bblock_GP')    
    
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*) "Computing Bblock_GP"
#endif

    ! initialize arrays
    Urows = 0
    Ucols = 0
    Uvals = 0.0d0

    Vrows = 0
    Vcols = 0
    Vvals = 0.0d0
    
    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = np+np-1

    ! ---------------------------------------------------------------------------
    ! precompute quantities needed for the Jacobian entries
    ! ---------------------------------------------------------------------------

    ! time level for current state
    np1 = fptr%tl%np1
    qn0 = fptr%n_Q 

    do k = 1,nlev
             
       ! pressure
       p(:,:,k) = (fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 & 
            + fptr%hvcoord%hybm(k) * fptr%base(ie)%state%ps_v(:,:,np1))
       
       p2(:,:,k) = p(:,:,k)**2

       ! vertical derivative of pressure (dp/dn)
       dpdn(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
            - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))
       
       ! -1 * dp/dps
       mBm(k) = -1.0d0 * fptr%hvcoord%hybm(k)
       
       ! d(dp/dn)/dps
       dBi(k) = fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)
    end do
    
    ! virtual temperature
    if (qn0 == -1) then
       T_v = fptr%base(ie)%state%T(:,:,:,np1)
    else
       Qt  = fptr%base(ie)%state%Qdp(:,:,:,1,qn0) / dpdn
       T_v = Virtual_Temperature(fptr%base(ie)%state%T(:,:,:,np1), Qt)
    end if

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp
    
    ! ---------------------------------------------------------------------------
    ! compute Jacobian entires 
    ! ---------------------------------------------------------------------------
    do k = 1,nlev
       do j = 1,np
          do i = 1,np
             
             ! get global row index for U and V equation
             ridx = ridx + 1            
             Urows(ridx) = fptr%base(ie)%gdofP(i,j) + UOffSet(k)
             Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + VOffSet(k)
             
             temp11 = massmatrix(i,j) * rrearth * fptr%base(ie)%Dinv(i,j,1,1) * rgas
             temp12 = massmatrix(i,j) * rrearth * fptr%base(ie)%Dinv(i,j,1,2) * rgas

             temp21 = massmatrix(i,j) * rrearth * fptr%base(ie)%Dinv(i,j,2,1) * rgas
             temp22 = massmatrix(i,j) * rrearth * fptr%base(ie)%Dinv(i,j,2,2) * rgas
                         
             ! a /= i, b = j
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle
                
                temp = T_v(a,j,k) * 0.5d0 & 
                     * (mBm(k) * dpdn(a,j,k)/p2(a,j,k) + dBi(k)/p(a,j,k)) 
                
                do m = k+1,nlev
                   temp = temp + T_v(a,j,m) &
                        * (mBm(m) * dpdn(a,j,m)/p2(a,j,m) + dBi(m)/p(a,j,m))
                end do
                
                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffSet
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffSet
                
                ! dU(i,j,k) / dPs(a,b)
                Uvals(cidx) = temp11 * temp * fptr%deriv%Dvv(a,i)

                ! dV(i,j,k) / dPs(a,b)
                Vvals(cidx) = temp12 * temp * fptr%deriv%Dvv(a,i)
                
             end do
             
             ! a = i, b /= j
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle
                
                temp = T_v(i,b,k) * 0.5d0 &
                     * (mBm(k) * dpdn(i,b,k)/p2(i,b,k) + dBi(k)/p(i,b,k))
                
                do m = k+1,nlev
                   temp = temp + T_v(i,b,m) &
                        * (mBm(m) * dpdn(i,b,m)/p2(i,b,m) + dBi(m)/p(i,b,m))
                end do

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffSet
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffSet

                ! du(i,j,k) / dps(a,b)               
                Uvals(cidx) = temp21 * temp * fptr%deriv%Dvv(b,j)

                ! dv(i,j,k) / dps(a,b)
                Vvals(cidx) = temp22 * temp * fptr%deriv%Dvv(b,j)
                
             end do
             
             ! a = i, b = j, c = k (x and y-derivatives)
             ! ---------------------------------------------------------------
             temp = T_v(i,j,k) * 0.5d0 &
                  * (mBm(k) * dpdn(i,j,k)/p2(i,j,k) + dBi(k)/p(i,j,k))
             
             do m = k+1,nlev
                temp = temp + T_v(i,j,m) &
                     * (mBm(m) * dpdn(i,j,m)/p2(i,j,m) + dBi(m)/p(i,j,m))
             end do

             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffSet
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffSet
             
             ! dU(i,j,k) / dPs(a,b)
             Uvals(cidx) = massmatrix(i,j) * rrearth * rgas * temp * &
                  (fptr%base(ie)%Dinv(i,j,1,1) * fptr%deriv%Dvv(i,i) &
                  + fptr%base(ie)%Dinv(i,j,2,1) * fptr%deriv%Dvv(j,j))

             ! dV(i,j,k) / dPs(a,b)            
             Vvals(cidx) = massmatrix(i,j) * rrearth * rgas * temp * &
                  (fptr%base(ie)%Dinv(i,j,1,2) * fptr%deriv%Dvv(i,i) &
                  + fptr%base(ie)%Dinv(i,j,2,2) * fptr%deriv%Dvv(j,j))
             
          end do
       end do
    end do
   
    call t_stopf('Jac_Bblock_GP')    

  end subroutine Jac_Bblock_GP



  subroutine Jac_Bblock_VertAdv(ie, nnz_max, row_nnz, Urows, Ucols, Uvals, &
       Vrows, Vcols, Vvals, c_ptr_to_object) bind(C,name="Jac_Bblock_VertAdv")
    !----------------------------------------------------------------------------
    ! Computes contributions to the A block of the Jacobian from the vertical
    ! advection term of the residual equation
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use derivative_mod, only        : divergence_sphere
    
    implicit none

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    ! row indices 
    integer(c_int), dimension(np*np*nlev), intent(out) :: Urows
    integer(c_int), dimension(np*np*nlev), intent(out) :: Vrows

    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Ucols
    integer(c_int), dimension(nnz_max), intent(out) :: Vcols

    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Uvals   
    real(c_double), dimension(nnz_max), intent(out) :: Vvals
    
    type(c_ptr) :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,np,np,nlev+1) :: dps_etadot_dpdn

    real(kind=real_kind), dimension(np,np,nlev+1) :: etadot_dpdn

    real(kind=real_kind), dimension(np,np,nlev) :: dpdn
    real(kind=real_kind), dimension(np,np,nlev) :: rdpdn
    real(kind=real_kind), dimension(np,np,nlev) :: half_rdpdn
    real(kind=real_kind), dimension(np,np,nlev) :: div_dpdn_vel

    real(kind=real_kind), dimension(np,np,2) :: dpdn_vel   

    real(kind=real_kind), dimension(np,np) :: sdot_sum
    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    real(kind=real_kind) :: half_rdpdn_diffu1, half_rdpdn_diffu2
    real(kind=real_kind) :: half_rdpdn_diffv1, half_rdpdn_diffv2
    real(kind=real_kind) :: temp1, temp2, diffB
    
    integer :: np1        ! current time level
    integer :: i, j, k    ! local equation index
    integer :: a, b       ! local derivative index
    integer :: ridx, cidx ! index of U and V output arrays
    !----------------------------------------------------------------------------
    
    call t_startf('Jac_Bblock_VertAdv')    

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Urows = 0
    Ucols = 0
    Uvals = 0.0d0

    Vrows = 0
    Vcols = 0
    Vvals = 0.0d0

    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = np+np-1

    ! ---------------------------------------------------------------------------
    ! precompute quantities needed for the Jacobian entries
    ! ---------------------------------------------------------------------------

    ! time level where current iteration is stored
    np1 = fptr%tl%np1
       
    ! vertical derivative of pressure, dpdn
    do k = 1,nlev
       dpdn(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
            - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))
    end do

    rdpdn = 1.0d0 / dpdn
    
    half_rdpdn = 0.5d0 / dpdn
    
    ! divergence of dpdn * velocity, needed for etadot_dpdn
    do k = 1,nlev
       do j = 1,np
          do i = 1,np             
             dpdn_vel(i,j,1) = dpdn(i,j,k) * fptr%base(ie)%state%v(i,j,1,k,np1)
             dpdn_vel(i,j,2) = dpdn(i,j,k) * fptr%base(ie)%state%v(i,j,2,k,np1)
          end do
       end do
       
       div_dpdn_vel(:,:,k) = divergence_sphere(dpdn_vel, fptr%deriv, fptr%base(ie))
    end do
    
    ! compute etadot_dpdn
    sdot_sum = 0.0d0

    do k=1,nlev
       sdot_sum(:,:)        = sdot_sum(:,:) + div_dpdn_vel(:,:,k)
       etadot_dpdn(:,:,k+1) = sdot_sum(:,:)
    end do

    do k=1,nlev-1
       etadot_dpdn(:,:,k+1) = fptr%hvcoord%hybi(k+1)*sdot_sum(:,:) &
            - etadot_dpdn(:,:,k+1)
    end do

    etadot_dpdn(:,:,1)      = 0.0d0
    etadot_dpdn(:,:,nlev+1) = 0.0d0

    ! compute velocity derivatives of etadot_dpdn
    call DerivPs_etadot_dpdn(ie, fptr, dps_etadot_dpdn)

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp

    ! ---------------------------------------------------------------------------
    ! compute Jacobian entires 
    ! ---------------------------------------------------------------------------

    ! Level 1
    ! ---------------------------------------------------------------------------
    k = 1
    do j = 1,np
       do i = 1,np

          ! get global row index for U and V equation
          ridx = ridx + 1            
          Urows(ridx) = fptr%base(ie)%gdofP(i,j) + UOffSet(k)
          Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + VOffSet(k)

          half_rdpdn_diffu1 = half_rdpdn(i,j,k) * &
               (fptr%base(ie)%state%v(i,j,1,k+1,np1) &
               - fptr%base(ie)%state%v(i,j,1,k,np1))

          half_rdpdn_diffv1 = half_rdpdn(i,j,k) * &
               (fptr%base(ie)%state%v(i,j,2,k+1,np1) &
               - fptr%base(ie)%state%v(i,j,2,k,np1))

          ! a /= i, b = j (skips zeros in dps_etadot_dpdn)
          ! ------------------------------------------------------------------
          do a = 1,np
             if (a == i) cycle

             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffSet
             Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffSet

             ! dU(i,j,k) / dPs(a,b)
             Uvals(cidx) = massmatrix(i,j) * dps_etadot_dpdn(a,j,i,j,k+1) &
                           * half_rdpdn_diffu1

             ! dV(i,j,k) / dPs(a,b)
             Vvals(cidx) = massmatrix(i,j) * dps_etadot_dpdn(a,j,i,j,k+1) &
                           * half_rdpdn_diffv1
          end do

          ! a = i, b /= j (skips zeros in dps_etadot_dpdn)
          ! ------------------------------------------------------------------
          do b = 1,np
             if (b == j) cycle
             
             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffSet
             Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffSet

             ! dU(i,j,k) / dPs(a,b)
             Uvals(cidx) = massmatrix(i,j) *  dps_etadot_dpdn(i,b,i,j,k+1) &
                           * half_rdpdn_diffu1
             
             ! dV(i,j,k) / dPs(a,b)
             Vvals(cidx) = massmatrix(i,j) * dps_etadot_dpdn(i,b,i,j,k+1) &
                           * half_rdpdn_diffv1
          end do

          ! a = i, b = j
          ! ------------------------------------------------------------------
          diffB = fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)

          temp1 = -0.5d0 * diffB * rdpdn(i,j,k)**2 * etadot_dpdn(i,j,k+1)

          cidx = cidx + 1              
          Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffSet 
          Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffSet

          ! dU(i,j,k) / dPs(a,b)
          Uvals(cidx) = massmatrix(i,j) * & 
               (dps_etadot_dpdn(i,j,i,j,k+1) * half_rdpdn_diffu1 &
               + temp1 * (fptr%base(ie)%state%v(i,j,1,k+1,np1) &
               - fptr%base(ie)%state%v(i,j,1,k,np1)))

          ! dV(i,j,k) / dPs(a,b)
          Vvals(cidx) = massmatrix(i,j) * &
               (dps_etadot_dpdn(i,j,i,j,k+1) * half_rdpdn_diffv1 &
               + temp1 * (fptr%base(ie)%state%v(i,j,2,k+1,np1) &
               - fptr%base(ie)%state%v(i,j,2,k,np1)))
       end do
    end do

    ! Levels 2 to nelv-1
    ! ---------------------------------------------------------------------------
    do k = 2,nlev-1
       do j = 1,np
          do i = 1,np

             ridx = ridx + 1            
             Urows(ridx) = fptr%base(ie)%gdofP(i,j) + UOffSet(k)
             Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + VOffSet(k)

             half_rdpdn_diffu1 = half_rdpdn(i,j,k) &
                  * (fptr%base(ie)%state%v(i,j,1,k+1,np1) &
                  - fptr%base(ie)%state%v(i,j,1,k,np1))

             half_rdpdn_diffv1 = half_rdpdn(i,j,k) &
                  * (fptr%base(ie)%state%v(i,j,2,k+1,np1) &
                  - fptr%base(ie)%state%v(i,j,2,k,np1))

             half_rdpdn_diffu2 = half_rdpdn(i,j,k) &
                  * (fptr%base(ie)%state%v(i,j,1,k,np1) & 
                  - fptr%base(ie)%state%v(i,j,1,k-1,np1))

             half_rdpdn_diffv2 = half_rdpdn(i,j,k) &
                  * (fptr%base(ie)%state%v(i,j,2,k,np1) &
                  - fptr%base(ie)%state%v(i,j,2,k-1,np1))

             ! a /= i, b = j (skips zeros in dps_etadot_dpdn)
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffSet
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffSet

                ! dU(i,j,k) / dPs(a,b)
                Uvals(cidx) = massmatrix(i,j) * &
                     (dps_etadot_dpdn(a,j,i,j,k+1) * half_rdpdn_diffu1 &
                     + dps_etadot_dpdn(a,j,i,j,k) * half_rdpdn_diffu2)

                ! dV(i,j,k) / dPs(a,b)
                Vvals(cidx) = massmatrix(i,j) * &
                     (dps_etadot_dpdn(a,j,i,j,k+1) * half_rdpdn_diffv1 &
                     + dps_etadot_dpdn(a,j,i,j,k) * half_rdpdn_diffv2)
             end do

             ! a = i, b /= j (skips zeros in dps_etadot_dpdn)
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffSet
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffSet

                ! dU(i,j,k) / dPs(a,b)
                Uvals(cidx) = massmatrix(i,j) * &
                     (dps_etadot_dpdn(i,b,i,j,k+1) * half_rdpdn_diffu1 &
                     + dps_etadot_dpdn(i,b,i,j,k) * half_rdpdn_diffu2)

                ! dV(i,j,k) / dPs(a,b)
                Vvals(cidx) = massmatrix(i,j) * &
                     (dps_etadot_dpdn(i,b,i,j,k+1) * half_rdpdn_diffv1 &
                     + dps_etadot_dpdn(i,b,i,j,k) * half_rdpdn_diffv2)
             end do

             ! a = i and b = j
             ! ------------------------------------------------------------------
             diffB = fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)

             temp1 = -0.5d0 * diffB * rdpdn(i,j,k)**2 * etadot_dpdn(i,j,k+1) ! etadot_dpdn at interface k+1/2
             temp2 = -0.5d0 * diffB * rdpdn(i,j,k)**2 * etadot_dpdn(i,j,k)   ! etadot_dpdn at interface k-1/2

             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffSet
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffSet

             ! dU(i,j,k) / dPs(a,b)
             Uvals(cidx) = massmatrix(i,j) * &
                  (dps_etadot_dpdn(i,j,i,j,k+1) * half_rdpdn_diffu1 &
                  + dps_etadot_dpdn(i,j,i,j,k) * half_rdpdn_diffu2 &
                  + temp1 * (fptr%base(ie)%state%v(i,j,1,k+1,np1) - fptr%base(ie)%state%v(i,j,1,k,np1)) &
                  + temp2 * (fptr%base(ie)%state%v(i,j,1,k,np1)   - fptr%base(ie)%state%v(i,j,1,k-1,np1)))

             ! dV(i,j,k) / dPs(a,b)
             Vvals(cidx) = massmatrix(i,j) * &
                  (dps_etadot_dpdn(i,j,i,j,k+1) * half_rdpdn_diffv1 &
                  + dps_etadot_dpdn(i,j,i,j,k) * half_rdpdn_diffv2 &
                  + temp1 * (fptr%base(ie)%state%v(i,j,2,k+1,np1) - fptr%base(ie)%state%v(i,j,2,k,np1)) &
                  + temp2 * (fptr%base(ie)%state%v(i,j,2,k,np1)   - fptr%base(ie)%state%v(i,j,2,k-1,np1)))
          end do
       end do
    end do

    ! Level nlev
    ! ---------------------------------------------------------------------------
    k = nlev
    do j = 1,np
       do i = 1,np

          ridx = ridx + 1            
          Urows(ridx) = fptr%base(ie)%gdofP(i,j) + UOffSet(k)
          Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + VOffSet(k)
          
          half_rdpdn_diffu2 = half_rdpdn(i,j,k) &
               * (fptr%base(ie)%state%v(i,j,1,k,np1) &
               - fptr%base(ie)%state%v(i,j,1,k-1,np1))
          
          half_rdpdn_diffv2 = half_rdpdn(i,j,k) &
               * (fptr%base(ie)%state%v(i,j,2,k,np1) &
               - fptr%base(ie)%state%v(i,j,2,k-1,np1))

          ! a /= i, b = j (skip zeros in dps_etadot_dpdn)
          ! ---------------------------------------------------------------------
          do a = 1,np
             if (a == i) cycle

             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffSet
             Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffSet

             ! dU(i,j,k) / dps(a,b)
             Uvals(cidx) = massmatrix(i,j) * dps_etadot_dpdn(a,j,i,j,k) &
                           * half_rdpdn_diffu2

             ! dV(i,j,k) / dps(a,b)
             Vvals(cidx) = massmatrix(i,j) * dps_etadot_dpdn(a,j,i,j,k) &
                           * half_rdpdn_diffv2
          end do

          ! a = i, b /= j (to skip zeros in dps_etadot_dpdn)
          ! ---------------------------------------------------------------------
          do b = 1,np
             if (b == j) cycle

             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffSet
             Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffSet

             ! dU(i,j,k) / dPs(a,b)
             Uvals(cidx) = massmatrix(i,j) * dps_etadot_dpdn(i,b,i,j,k) &
                           * half_rdpdn_diffu2

             ! dV(i,j,k) / dPs(a,b)
             Vvals(cidx) = massmatrix(i,j) * dps_etadot_dpdn(i,b,i,j,k) &
                           * half_rdpdn_diffv2
          end do

          ! a = i, b = j
          ! ---------------------------------------------------------------------
          diffB = fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)

          temp2 = -0.5d0 * diffB * rdpdn(i,j,k)**2 * etadot_dpdn(i,j,k) ! etadot_dpdn at interface k-1/2

          cidx = cidx + 1              
          Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffSet
          Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffSet
          
          ! dU(i,j,k) / dPs(a,b)
          Uvals(cidx) = massmatrix(i,j) * &
               (dps_etadot_dpdn(i,j,i,j,k) * half_rdpdn_diffu2 &
               + temp2 * (fptr%base(ie)%state%v(i,j,1,k,np1) &
               - fptr%base(ie)%state%v(i,j,1,k-1,np1)))

          ! dV(i,j,k) / dPs(a,b)
          Vvals(cidx) = massmatrix(i,j) * &
               (dps_etadot_dpdn(i,j,i,j,k) * half_rdpdn_diffv2 &
               + temp2 * (fptr%base(ie)%state%v(i,j,2,k,np1) &
               - fptr%base(ie)%state%v(i,j,2,k-1,np1)))
       end do
    end do

    call t_stopf('Jac_Bblock_VertAdv')    

  end subroutine Jac_Bblock_VertAdv



  subroutine Jac_Bblock_gradP(ie, nnz_max, row_nnz, Urows, Ucols, Uvals, &
       Vrows, Vcols, Vvals, c_ptr_to_object) bind(C,name="Jac_Bblock_gradP")
    !----------------------------------------------------------------------------
    ! Computes contributions to the B block of the Jacobian from grad(p)
    ! term of the residual equation
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth, rgas
    use physics_mod, only           : Virtual_Temperature
    use derivative_mod, only        : gradient_sphere
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    ! row indices 
    integer(c_int), dimension(np*np*nlev), intent(out) :: Urows
    integer(c_int), dimension(np*np*nlev), intent(out) :: Vrows

    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Ucols
    integer(c_int), dimension(nnz_max), intent(out) :: Vcols

    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Uvals   
    real(c_double), dimension(nnz_max), intent(out) :: Vvals
    
    type(c_ptr) :: c_ptr_to_object
    ! ---------------------------------------------------------------------------
    
    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,2,nlev) :: grad_p ! gradient of pressure

    real(kind=real_kind), dimension(np,np,nlev)   :: T_v ! virtual temperature
    real(kind=real_kind), dimension(np,np,nlev)   :: p   ! pressure
    real(kind=real_kind), dimension(np,np,nlev)   :: p2  ! pressure squared
    real(kind=real_kind), dimension(np,np,nlev)   :: dp  ! vertical derivative of pressure
    real(kind=real_kind), dimension(np,np,nlev)   :: Qt

    real(kind=real_kind), dimension(np,np,2) :: grad_ps ! gradient of surface pressure

    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    real(kind=real_kind) :: temp1, temp2

    integer :: np1         ! time level
    integer :: qn0         ! dry or moist flag

    integer :: a, b        ! derivative index
    integer :: i, j, k     ! equation index
    integer :: ridx, cidx  ! global row/column index 
    !----------------------------------------------------------------------------

    call t_startf('Jac_Bblock_gradP')    
    
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Urows = 0
    Ucols = 0
    Uvals = 0.0d0

    Vrows = 0
    Vcols = 0
    Vvals = 0.0d0
    
    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = np+np-1

    ! ---------------------------------------------------------------------------
    ! precompute quantities needed for the Jacobian entries
    ! ---------------------------------------------------------------------------

    ! time level for current state
    np1 = fptr%tl%np1
    qn0 = fptr%n_Q

    ! virtual temperature
    if (qn0 == -1) then
       T_v(:,:,:)   = fptr%base(ie)%state%T(:,:,:,np1)   
    else

       do k = 1,nlev
          dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
               - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))
       end do
       
       Qt(:,:,:)  = fptr%base(ie)%state%Qdp(:,:,:,1,qn0)/dp(:,:,:)
       T_v(:,:,:) = Virtual_Temperature(fptr%base(ie)%state%T(:,:,:,np1),Qt) 
    end if
    
    ! gradient of surface pressure and pressure
    grad_ps = gradient_sphere(fptr%base(ie)%state%ps_v(:,:,np1), fptr%deriv, fptr%base(ie)%Dinv)

    do k = 1,nlev
       grad_p(:,:,:,k) = fptr%hvcoord%hybm(k) * grad_ps(:,:,:)
    end do
    
    ! pressure and pressure squared
    do k = 1,nlev
       p(:,:,k)  = (fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybm(k) * fptr%base(ie)%state%ps_v(:,:,np1))

       p2(:,:,k) = p(:,:,k)**2
    end do

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp

    ! ---------------------------------------------------------------------------
    ! compute Jacobian entires 
    ! ---------------------------------------------------------------------------
    do k = 1,nlev
       do j = 1,np
          do i = 1,np
             
             ! get global row index for U and V equation
             ridx = ridx + 1            
             Urows(ridx) = fptr%base(ie)%gdofP(i,j) + UOffSet(k)
             Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + VOffSet(k)
             
             temp1 = massmatrix(i,j) * fptr%hvcoord%hybm(k) * Rgas * T_v(i,j,k) &
                     * rrearth * fptr%base(ie)%Dinv(i,j,1,1) / p(i,j,k)

             temp2 = massmatrix(i,j) * fptr%hvcoord%hybm(k) * Rgas * T_v(i,j,k) &
                     * rrearth * fptr%base(ie)%Dinv(i,j,1,2) / p(i,j,k)
             
             ! a /= i, b = j (x-derivatives in gradient)
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                ! increment storage index, get global row/column index
                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffSet
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffSet

                ! dU(i,j,k) / dPs(a,b)
                Uvals(cidx) = temp1 * fptr%deriv%Dvv(a,i)

                ! dV(i,j,k) / dPs(a,b)
                Vvals(cidx) = temp2 * fptr%deriv%Dvv(a,i)

             end do

             temp1 = massmatrix(i,j) * fptr%hvcoord%hybm(k) * Rgas * T_v(i,j,k) &
                     * rrearth * fptr%base(ie)%Dinv(i,j,2,1) / p(i,j,k)

             temp2 = massmatrix(i,j) * fptr%hvcoord%hybm(k) * Rgas * T_v(i,j,k) &
                     * rrearth * fptr%base(ie)%Dinv(i,j,2,2) / p(i,j,k)

             ! a = i, b /= j (y-derivatives in gradient) 
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle  

                ! increment storage index, get global row/column index
                cidx = cidx + 1              
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffSet
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffSet

                ! dU(i,j,k) / dPs(a,b)
                Uvals(cidx) = temp1 * fptr%deriv%Dvv(b,j)

                ! dV(i,j,k) / dPs(a,b)
                Vvals(cidx) = temp2 * fptr%deriv%Dvv(b,j)

             end do

             ! a = i, b = j(x- and y-derivatives in gradient)
             ! ------------------------------------------------------------------
             temp1 = massmatrix(i,j) * (-1.0d0) * fptr%hvcoord%hybm(k) * Rgas &
                     * T_v(i,j,k) / p2(i,j,k)

             temp2 = massmatrix(i,j) * fptr%hvcoord%hybm(k) * Rgas &
                     * T_v(i,j,k) * rrearth / p(i,j,k)

             ! increment storage index, get global row/column index
             cidx = cidx + 1              
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffSet
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffSet

             ! dU(i,j,k) / dPs(a,b)
             Uvals(cidx) = temp1 * grad_p(i,j,1,k) + temp2 &
                  * (fptr%base(ie)%Dinv(i,j,1,1) * fptr%deriv%Dvv(i,i) &
                  + fptr%base(ie)%Dinv(i,j,2,1) * fptr%deriv%Dvv(j,j))

             ! dV(i,j,k) / dPs(a,b)
             Vvals(cidx) = temp1 * grad_p(i,j,2,k) + temp2 &
                  * (fptr%base(ie)%Dinv(i,j,1,2) * fptr%deriv%Dvv(i,i) &
                  + fptr%base(ie)%Dinv(i,j,2,2) * fptr%deriv%Dvv(j,j))

          end do
       end do
    end do

    call t_stopf('Jac_Bblock_gradP')    

  end subroutine Jac_Bblock_gradP



  ! =============================================================================
  ! C Block of the Jacobian
  ! =============================================================================



  subroutine Jac_Cblock_GP(ie, nnz_max, row_nnz, Urows, Ucols, Uvals, &
       Vrows, Vcols, Vvals, c_ptr_to_object) bind(C,name="Jac_Cblock_GP")
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the gradient of the
    ! geopotential term in the velocity equations
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth, rgas
    use physics_mod, only           : Virtual_Temperature

    implicit none

    !---------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), dimension(np*np*nlev), intent(out) :: row_nnz 

    ! row indices 
    integer(c_int), dimension(np*np*nlev), intent(out) :: Urows
    integer(c_int), dimension(np*np*nlev), intent(out) :: Vrows

    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Ucols
    integer(c_int), dimension(nnz_max), intent(out) :: Vcols

    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Uvals   
    real(c_double), dimension(nnz_max), intent(out) :: Vvals
    
    type(c_ptr) :: c_ptr_to_object
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,nlev) :: p    ! pressure
    real(kind=real_kind), dimension(np,np,nlev) :: dpdn ! vertical derivative of pressure
    real(kind=real_kind), dimension(np,np,nlev) :: Qt   !
    real(kind=real_kind), dimension(np,np,nlev) :: c1   ! rgas * scaling factor for Tv

    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    real(kind=real_kind) :: temp, temp11, temp12, temp21, temp22

    integer :: np1         ! time level
    integer :: a, b, c     ! derivative index
    integer :: i, j, k     ! equation index
    integer :: ridx, cidx  ! global row/column index 
    integer :: qn0         ! dry / moist flag
    !----------------------------------------------------------------------------

    call t_startf('Jac_Cblock_GP')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Urows = 0
    Ucols = 0
    Uvals = 0.0d0

    Vrows = 0
    Vcols = 0
    Vvals = 0.0d0
    
    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! initialize number of nonzeros per row (varies)
    row_nnz = 0

    ! ---------------------------------------------------------------------------
    ! precompute quantities needed for the Jacobian entries
    ! ---------------------------------------------------------------------------

    ! time level where current iteration is stored
    np1 = fptr%tl%np1

    qn0 = fptr%n_Q 

    do k = 1,nlev
       ! pressure
       p(:,:,k) = (fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 & 
            + fptr%hvcoord%hybm(k) * fptr%base(ie)%state%ps_v(:,:,np1))

       ! vertical derivative of pressure (dp/dn)
       dpdn(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
            - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))       
    end do

    ! virtual temperature
    if (qn0 == -1) then
       c1 = rgas
    else
       Qt = fptr%base(ie)%state%Qdp(:,:,:,1,qn0) / dpdn
       c1 = rgas
       c1 = Virtual_Temperature(c1, Qt)
    end if

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp

    ! temperature derivatives, c = k and c > k
    ! ---------------------------------------------------------------------------
    do k = 1,nlev
       do j = 1,np
          do i = 1,np

             ! set row index for U and V equations
             ridx = ridx + 1
             Urows(ridx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
             Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)

             temp11 = massmatrix(i,j) * rrearth * fptr%base(ie)%Dinv(i,j,1,1)
             temp12 = massmatrix(i,j) * rrearth * fptr%base(ie)%Dinv(i,j,1,2)

             temp21 = massmatrix(i,j) * rrearth * fptr%base(ie)%Dinv(i,j,2,1)
             temp22 = massmatrix(i,j) * rrearth * fptr%base(ie)%Dinv(i,j,2,2)

             ! a /= i, b = j, c = k (x-derivative)
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Toffset(k)
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Toffset(k)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1

                temp = c1(a,j,k) * 0.5d0 * dpdn(a,j,k) &
                     * fptr%deriv%Dvv(a,i) / p(a,j,k)

                ! dU(i,j,k) / dT(a,b,c)
                Uvals(cidx) = temp11 * temp

                ! dV(i,j,k) / dT(a,b,c)
                Vvals(cidx) = temp12 * temp
             end do

             ! a = i, b /= j, c = k (y-derivative)
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Toffset(k)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Toffset(k)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1

                temp = c1(i,b,k) * 0.5d0 * dpdn(i,b,k) &
                     * fptr%deriv%Dvv(b,j) / p(i,b,k)

                ! dU(i,j,k) / dT(a,b,c)
                Uvals(cidx) = temp21 * temp

                ! dV(i,j,k) / dT(a,b,c)
                Vvals(cidx) = temp22 * temp
             end do

             ! a = i, b = j, c = k (x and y-derivatives)
             ! ------------------------------------------------------------------
             cidx = cidx + 1
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

             ! update nnz per row
             row_nnz(ridx) = row_nnz(ridx) + 1

             temp = massmatrix(i,j) * rrearth * c1(i,j,k) * 0.5d0 * dpdn(i,j,k)/p(i,j,k)

             ! dU(i,j,k) / dT(a,b,c)
             Uvals(cidx) = temp * (fptr%base(ie)%Dinv(i,j,1,1) &
                  * fptr%deriv%Dvv(i,i) + fptr%base(ie)%Dinv(i,j,2,1) &
                  * fptr%deriv%Dvv(j,j))

             ! dV(i,j,k) / dT(a,b,c)
             Vvals(cidx) = temp * (fptr%base(ie)%Dinv(i,j,1,2) &
                  * fptr%deriv%Dvv(i,i) + fptr%base(ie)%Dinv(i,j,2,2) &
                  * fptr%deriv%Dvv(j,j))

             do c = k+1,nlev

                ! a /= i, b = j, c > k
                ! ---------------------------------------------------------------
                do a = 1,np
                   if (a == i) cycle

                   cidx = cidx + 1
                   Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Toffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Toffset(c)

                   ! update nnz per row
                   row_nnz(ridx) = row_nnz(ridx) + 1

                   temp = c1(a,j,c) * dpdn(a,j,c) * fptr%deriv%Dvv(a,i)/p(a,j,c) 

                   ! dU(i,j,k) / dT(a,b,c)
                   Uvals(cidx) = temp11 * temp

                   ! dV(i,j,k) / dT(a,b,c)
                   Vvals(cidx) = temp12 * temp
                end do

                ! a = i, b /= j, c > k
                !----------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle

                   cidx = cidx + 1
                   Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Toffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Toffset(c)

                   ! update nnz per row
                   row_nnz(ridx) = row_nnz(ridx) + 1

                   temp = c1(i,b,c) * dpdn(i,b,c) * fptr%deriv%Dvv(b,j)/p(i,b,c)

                   ! dU(i,j,k) / dT(a,b,c)
                   Uvals(cidx) = temp21 * temp

                   ! dV(i,j,k) / dT(a,b,c)
                   Vvals(cidx) = temp22 * temp   
                end do

                ! a = i, b = j, c > k
                ! ---------------------------------------------------------------
                cidx = cidx + 1
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(c)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1

                temp = massmatrix(i,j) * rrearth * c1(i,j,c) * dpdn(i,j,c)/p(i,j,c)

                ! dU(i,j,k) / dT(a,b,c)
                Uvals(cidx) = temp * (fptr%base(ie)%Dinv(i,j,1,1) &
                     * fptr%deriv%Dvv(i,i) + fptr%base(ie)%Dinv(i,j,2,1) &
                     * fptr%deriv%Dvv(j,j))

                ! dV(i,j,k) / dT(a,b,c)
                Vvals(cidx) = temp * (fptr%base(ie)%Dinv(i,j,1,2) &
                     * fptr%deriv%Dvv(i,i) + fptr%base(ie)%Dinv(i,j,2,2) &
                     * fptr%deriv%Dvv(j,j))
             end do

          end do
       end do
    end do

    call t_stopf('Jac_Cblock_GP')

  end subroutine Jac_Cblock_GP



  subroutine Jac_Cblock_gradP(ie, nnz_max, row_nnz, Urows, Ucols, Uvals, &
       Vrows, Vcols, Vvals, c_ptr_to_object) bind(C,name="Jac_Cblock_gradP")
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the pressure gradient
    ! term in the velocity equations
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use derivative_mod, only        : gradient_sphere
    use physical_constants, only    : rrearth, Rgas
    use physics_mod, only           : Virtual_Temperature

    implicit none

    !---------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    ! row indices 
    integer(c_int), dimension(np*np*nlev), intent(out) :: Urows
    integer(c_int), dimension(np*np*nlev), intent(out) :: Vrows

    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Ucols
    integer(c_int), dimension(nnz_max), intent(out) :: Vcols

    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Uvals   
    real(c_double), dimension(nnz_max), intent(out) :: Vvals
    
    type(c_ptr) :: c_ptr_to_object
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,2,nlev) :: grad_p ! gradient of pressure

    real(kind=real_kind), dimension(np,np,nlev) :: p  ! pressure
    real(kind=real_kind), dimension(np,np,nlev) :: dp ! vertical derivative of pressure
    real(kind=real_kind), dimension(np,np,nlev) :: Qt !
    real(kind=real_kind), dimension(np,np,nlev) :: c1 ! 

    real(kind=real_kind), dimension(np,np,2) :: grad_ps ! gradient of surface pressure

    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    integer :: np1         ! time level

    integer :: i, j, k     ! equation index
    integer :: ridx, cidx  ! global row/column index 
    integer :: qn0         ! dry / moist flag
    !----------------------------------------------------------------------------

    call t_startf('Jac_Cblock_gradP')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Urows = 0
    Ucols = 0
    Uvals = 0.0d0

    Vrows = 0
    Vcols = 0
    Vvals = 0.0d0
    
    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = 1

    ! ---------------------------------------------------------------------------
    ! precompute quantities needed for the Jacobian entries
    ! ---------------------------------------------------------------------------

    ! time level where current iteration is stored
    np1 = fptr%tl%np1

    qn0 = fptr%n_Q 

    if (qn0 == -1) then
       c1 = Rgas
    else
       do k = 1,nlev
          dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
               - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))
       end do
       
       Qt = fptr%base(ie)%state%Qdp(:,:,:,1,qn0)/dp(:,:,:)
       c1 = Rgas
       c1 = Virtual_temperature(c1, Qt)
    end if
    
    ! gradient of surface pressure
    grad_ps = gradient_sphere(fptr%base(ie)%state%ps_v(:,:,np1), fptr%deriv, &
         fptr%base(ie)%Dinv)
    
    do k = 1,nlev
       ! pressure
       p(:,:,k) = (fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybm(k) * fptr%base(ie)%state%ps_v(:,:,np1))

       ! gradient of pressure
       grad_p(:,:,:,k) = fptr%hvcoord%hybm(k) * grad_ps(:,:,:)
    end do

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp

    ! ---------------------------------------------------------------------------
    ! compute Jacobian entries (temperature derivative)
    ! ---------------------------------------------------------------------------
    do k = 1,nlev
       do j = 1,np
          do i = 1,np

             ridx = ridx + 1
             Urows(ridx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
             Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)

             ! a = i, b = j, c = k
             ! ------------------------------------------------------------------
             cidx = cidx + 1
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

             ! dU(i,j,k) / dT(a,b,c)
             Uvals(cidx) = massmatrix(i,j) * c1(i,j,k) * grad_p(i,j,1,k) / p(i,j,k)

             ! dV(i,j,k) / dT(a,b,c)
             Vvals(cidx) = massmatrix(i,j) * c1(i,j,k) * grad_p(i,j,2,k) / p(i,j,k)

          end do
       end do
    end do

    call t_stopf('Jac_Cblock_gradP')

  end subroutine Jac_Cblock_gradP



  ! =============================================================================
  ! D Block of the Jacobian
  ! =============================================================================



  subroutine Jac_Dblock(ie, nnz_max, row_nnz, Rows, Cols, Vals, c_ptr_to_object) &
       bind(C,name="Jac_Dblock")
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the integral of the 
    ! divergence of dp/dn * velocity in the surface pressure equation
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type

    implicit none

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows ! row indices 
    integer(c_int), dimension(nnz_max), intent(out)    :: Cols ! column indices
    real(c_double), dimension(nnz_max), intent(out)    :: Vals ! matrix entries
    
    type(c_ptr) :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------  
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    ! velocity derivatives of divergence(dpdn * velocity)
    real(kind=real_kind) :: du_div_dpdn_vel(np,np,nlev,np,np,nlev) ! a,b,c,i,j,k
    real(kind=real_kind) :: dv_div_dpdn_vel(np,np,nlev,np,np,nlev) ! a,b,c,i,j,k

    integer :: np1        ! current time level
    integer :: i, j       ! local equation index
    integer :: a, b, c    ! local derivative index
    integer :: ridx, cidx ! index of Ps output arrays
    !----------------------------------------------------------------------------

    call t_startf('Jac_Dblock') 

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0

    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = 2*(np+np-1)*nlev

    ! ---------------------------------------------------------------------------
    ! precompute quantities needed for the Jacobian entries
    ! ---------------------------------------------------------------------------

    ! time level where current iteration is stored
    np1 = fptr%tl%np1

    ! compute velocity derivatives of div(dpdn * velocity)
    call DerivVel_div_dpdn_vel(ie, fptr, du_div_dpdn_vel, dv_div_dpdn_vel)

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp

    ! ---------------------------------------------------------------------------
    ! compute Jacobian entires (velocity derivatives)
    ! ---------------------------------------------------------------------------
    do j = 1,np
       do i = 1,np              

          ! get global row index for Ps equation
          ridx = ridx + 1            
          Rows(ridx) = fptr%base(ie)%gdofP(i,j) + PsOffset

          ! U derivatives 
          ! ---------------------------------------------------------------------
          do c = 1,nlev
             
             ! a /= i, b = j, c = k
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)                

                Vals(cidx) = massmatrix(i,j) * du_div_dpdn_vel(a,j,c,i,j,c)
             end do
             
             ! a = i, b /= j, c = k
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)                

                Vals(cidx) = massmatrix(i,j) * du_div_dpdn_vel(i,b,c,i,j,c)
             end do
             
             ! a = i, b = j, c = k
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)

             Vals(cidx) = massmatrix(i,j) * du_div_dpdn_vel(i,j,c,i,j,c)
          end do

          ! V derivatives 
          ! ---------------------------------------------------------------------
          do c = 1,nlev
             
             ! a /= i, b = j, c = k
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)
                
                Vals(cidx) = massmatrix(i,j) * dv_div_dpdn_vel(a,j,c,i,j,c)
             end do
             
             ! a = i, b /= j, c = k
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)
                
                Vals(cidx) = massmatrix(i,j) * dv_div_dpdn_vel(i,b,c,i,j,c)
             end do
             
             ! a = i, b = j, c = k
             ! ------------------------------------------------------------------             
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)

             Vals(cidx) = massmatrix(i,j) * dv_div_dpdn_vel(i,j,c,i,j,c)            
          end do

       end do
    end do

    call t_stopf('Jac_Dblock') 
    
  end subroutine Jac_Dblock



  ! =============================================================================
  ! E Block of the Jacobian
  ! =============================================================================



  subroutine Jac_Eblock(ie, nnz_max, row_nnz, Rows, Cols, Vals, c_ptr_to_object) &
       bind(C,name="Jac_Eblock")
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the integral of the 
    ! time derivative and the divergence of dp/dn * velocity in the surface pressure equation
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use control_mod, only           : tstep_type
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type

    implicit none

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows ! row indices 
    integer(c_int), dimension(nnz_max), intent(out)    :: Cols ! column indices
    real(c_double), dimension(nnz_max), intent(out)    :: Vals ! matrix entries
    
    type(c_ptr) :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------  
    type(derived_type), pointer :: fptr=>NULL()

    ! surface pressure derivative of divergence(dpdn * velocity)
    real(kind=real_kind), dimension(np,np,np,np,nlev) :: dps_div_dpdn_vel ! a,b,-,h,k,l

    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    real(kind=real_kind) :: dti

    integer :: np1         ! current time level
    integer :: i, j        ! local equation index
    integer :: a, b        ! local derivative index
    integer :: ridx, cidx  ! index of Ps output arrays
    integer :: n           ! loop index for vertical integral
    !----------------------------------------------------------------------------

    call t_startf('Jac_Eblock') 

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0

    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = np+np-1

    ! ---------------------------------------------------------------------------
    ! precompute quantities needed for the Jacobian entries
    ! ---------------------------------------------------------------------------

    ! time level where current iteration is stored
    np1 = fptr%tl%np1

    dti = 1.0d0 / fptr%dt

    if (tstep_type==12) then
       dti = (3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2 on step size
    endif

    ! surface pressure derivatives of div(dpdn * velocity)
    call DerivPs_div_dpdn_vel(ie, fptr, dps_div_dpdn_vel)

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp

    ! ---------------------------------------------------------------------------
    ! compute Jacobian entires, surface pressure derivatives
    ! ---------------------------------------------------------------------------
    do j = 1,np
       do i = 1,np              
          
          ! get global row index for Ps equation
          ridx = ridx + 1            
          Rows(ridx) = fptr%base(ie)%gdofP(i,j) + PsOffset
          
          ! a /= i, b = j
          ! ---------------------------------------------------------------------
          do a = 1,np
             if (a == i) cycle                                      
             
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffset

             do n = 1,nlev
                Vals(cidx) = Vals(cidx) + dps_div_dpdn_vel(a,j,i,j,n)
             end do
             Vals(cidx) = massmatrix(i,j) * Vals(cidx)
          end do
          
          ! a = i, b /= j
          ! ---------------------------------------------------------------------
          do b = 1,np
             if (b == j) cycle                  

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffset

             do n = 1,nlev
                Vals(cidx) = Vals(cidx) + dps_div_dpdn_vel(i,b,i,j,n) 
             end do
             Vals(cidx) = massmatrix(i,j) * Vals(cidx)
          end do
          
          ! a = i, b = j
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Cols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffset

          Vals(cidx) = dti
          do n = 1,nlev
             Vals(cidx) = Vals(cidx) + dps_div_dpdn_vel(i,j,i,j,n) 
          end do
          Vals(cidx) = massmatrix(i,j) * Vals(cidx)
          
       end do
    end do
    
    call t_stopf('Jac_Eblock') 
    
  end subroutine Jac_Eblock



  ! =============================================================================
  ! F Block of the Jacobian
  ! =============================================================================



  subroutine Jac_Fblock_HorizAdv(ie, nnz_max, row_nnz, Rows, Cols, Vals, &
       c_ptr_to_object) bind(C,name="Jac_Fblock_HorizAdv")
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the horizontal advection 
    ! terms in the temperature equation
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth
    use derivative_mod, only        : gradient_sphere

    implicit none

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows ! row indices 
    integer(c_int), dimension(nnz_max), intent(out)    :: Cols ! column indices
    real(c_double), dimension(nnz_max), intent(out)    :: Vals ! matrix entries
    
    type(c_ptr) :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------  
    type(derived_type), pointer :: fptr=>NULL()
  
    real(kind=real_kind), dimension(np,np,2) :: grad_T ! gradient of temperature

    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    integer :: np1         ! current time level
    integer :: i, j, k     ! local equation index
    integer :: ridx, cidx  ! index of Ps output arrays
    !----------------------------------------------------------------------------

    call t_startf('Jac_Fblock_HorizAdv')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0

    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = 2

    ! time level where current iteration is stored
    np1 = fptr%tl%np1

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp
    
    ! ---------------------------------------------------------------------------
    ! compute Jacobian entires, velocity derivatives 
    ! ---------------------------------------------------------------------------
    do k = 1,nlev

       ! gradient of temperature
       grad_T = gradient_sphere(fptr%base(ie)%state%T(:,:,k,np1), &
            fptr%deriv, fptr%base(ie)%Dinv)
       
       do j = 1,np
          do i = 1,np

             ! get global row index for Ps equation
             ridx = ridx + 1            
             Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)
                          
             ! a = i, b = j, c = k
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)

             Vals(cidx) = massmatrix(i,j) * grad_T(i,j,1)

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)
             
             Vals(cidx) = massmatrix(i,j) * grad_T(i,j,2)
                          
          end do
       end do
    end do
    
    call t_stopf('Jac_Fblock_HorizAdv')

  end subroutine Jac_Fblock_HorizAdv



  subroutine Jac_Fblock_VertAdv(ie, nnz_max, row_nnz, Rows, Cols, Vals, &
       c_ptr_to_object) bind(C,name="Jac_Fblock_VertAdv")
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the vertical advection 
    ! terms in the velocity and temperature equations
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use derivative_mod, only        : divergence_sphere

    implicit none

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows ! row indices 
    integer(c_int), dimension(nnz_max), intent(out)    :: Cols ! column indices
    real(c_double), dimension(nnz_max), intent(out)    :: Vals ! matrix entries
    
    type(c_ptr) :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()
    
    ! velocity derivatives of etadot_dpdn
    real(kind=real_kind), dimension(np,np,nlev,np,np,nlev+1) :: du_eta_dot_dpdn
    real(kind=real_kind), dimension(np,np,nlev,np,np,nlev+1) :: dv_eta_dot_dpdn

    real(kind=real_kind), dimension(np,np,nlev+1) :: eta_dot_dpdn

    real(kind=real_kind), dimension(np,np,nlev) :: dpdn         ! vertical derivative of pressure

    real(kind=real_kind), dimension(np,np,nlev) :: half_rdpdn   ! 0.5 / dpdn
    real(kind=real_kind), dimension(np,np,nlev) :: div_dpdn_vel ! divergence(dpdn * velocity)

    real(kind=real_kind), dimension(np,np,2) :: dpdn_vel ! dpdn * velocity

    real(kind=real_kind), dimension(np,np) :: sdot_sum
    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    real(kind=real_kind) :: half_rdpdn_diffT1, half_rdpdn_diffT2


    integer :: np1        ! current time level
    integer :: i, j, k    ! local equation index
    integer :: a, b, c    ! local derivative index
    integer :: ridx, cidx ! index of Ps output arrays
    !----------------------------------------------------------------------------

    call t_startf('Jac_Fblock_VertAdv')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0

    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = 2*(np+np-1)*nlev

    ! ---------------------------------------------------------------------------
    ! precompute quantities needed for the Jacobian entries
    ! ---------------------------------------------------------------------------

    ! time level where current iteration is stored
    np1 = fptr%tl%np1

    ! vertical derivative of pressure, dpdn
    do k = 1,nlev
       dpdn(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
            - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))

       half_rdpdn(:,:,k) = 0.5d0 / dpdn(:,:,k)
    end do

    ! divergence of dpdn * velocity, needed for eta_dot_dpdn
    do k = 1,nlev
       do j=1,np
          do i=1,np             
             dpdn_vel(i,j,1) = dpdn(i,j,k) * fptr%base(ie)%state%v(i,j,1,k,np1)
             dpdn_vel(i,j,2) = dpdn(i,j,k) * fptr%base(ie)%state%v(i,j,2,k,np1)
          end do
       end do

       div_dpdn_vel(:,:,k) = divergence_sphere(dpdn_vel, fptr%deriv, fptr%base(ie))
    end do

    ! compute eta_dot_dpdn
    ! ------------------------------------------------------------------------
    sdot_sum = 0.0d0

    do k=1,nlev
       sdot_sum(:,:)         = sdot_sum(:,:) + div_dpdn_vel(:,:,k)
       eta_dot_dpdn(:,:,k+1) = sdot_sum(:,:)
    end do

    do k=1,nlev-1
       eta_dot_dpdn(:,:,k+1) = fptr%hvcoord%hybi(k+1)*sdot_sum(:,:) - eta_dot_dpdn(:,:,k+1)
    end do

    eta_dot_dpdn(:,:,1)      = 0.0d0
    eta_dot_dpdn(:,:,nlev+1) = 0.0d0

    ! compute the Jacobian of etadot_dpdn
    call DerivVel_etadot_dpdn(ie, fptr, du_eta_dot_dpdn, dv_eta_dot_dpdn) 

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp

    ! ------------------------------------------------------------------------
    ! Level 1 Jacobian entries, velocity derivatives
    ! ------------------------------------------------------------------------
    k = 1
    do j = 1,np
       do i = 1,np

          ! get global row index for T equation
          ridx = ridx + 1            
          Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

          half_rdpdn_diffT1 = massmatrix(i,j) * half_rdpdn(i,j,k) * &
               (fptr%base(ie)%state%T(i,j,k+1,np1) &
               - fptr%base(ie)%state%T(i,j,k,np1))

          do c = 1,nlev

             ! a /= i, b = j, c = * (skips zeros in dvel_eta_dot_dpdn)
             ! ------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)

                Vals(cidx) = du_eta_dot_dpdn(a,j,c,i,j,k+1) * half_rdpdn_diffT1

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)

                Vals(cidx) = dv_eta_dot_dpdn(a,j,c,i,j,k+1) * half_rdpdn_diffT1
             end do

             ! a = i, b /= j, c = * (skips zeros in dvel_eta_dot_dpdn)
             ! ------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)

                Vals(cidx) = du_eta_dot_dpdn(i,b,c,i,j,k+1) * half_rdpdn_diffT1

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)

                Vals(cidx) = dv_eta_dot_dpdn(i,b,c,i,j,k+1) * half_rdpdn_diffT1
             end do

             ! a = i, b = j
             ! ------------------------------------------------------------
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)

             Vals(cidx) = du_eta_dot_dpdn(i,j,c,i,j,k+1) * half_rdpdn_diffT1

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)

             Vals(cidx) = dv_eta_dot_dpdn(i,j,c,i,j,k+1) * half_rdpdn_diffT1
          end do

       end do
    end do

    ! ------------------------------------------------------------------------
    ! Level 2 to nlev-1 Jacobian Entries
    ! ------------------------------------------------------------------------
    do k = 2,nlev-1
       do j = 1,np
          do i = 1,np

             ! get global row index for Ps equation
             ridx = ridx + 1            
             Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

             half_rdpdn_diffT1 = massmatrix(i,j) * half_rdpdn(i,j,k) * &
                  (fptr%base(ie)%state%T(i,j,k+1,np1) &
                  - fptr%base(ie)%state%T(i,j,k,np1))

             half_rdpdn_diffT2 = massmatrix(i,j) * half_rdpdn(i,j,k) * &
                  (fptr%base(ie)%state%T(i,j,k,np1) &
                  - fptr%base(ie)%state%T(i,j,k-1,np1))

             do c = 1,nlev

                ! a /= i, b = j, c = * (skips zeros in dvel_eta_dot_dpdn)
                ! ------------------------------------------------------------
                do a = 1,np
                   if (a == i) cycle

                   cidx = cidx + 1              
                   Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)

                   Vals(cidx) = du_eta_dot_dpdn(a,j,c,i,j,k+1) * half_rdpdn_diffT1 &
                        + du_eta_dot_dpdn(a,j,c,i,j,k) * half_rdpdn_diffT2

                   cidx = cidx + 1              
                   Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)

                   Vals(cidx) = dv_eta_dot_dpdn(a,j,c,i,j,k+1) * half_rdpdn_diffT1 &
                        + dv_eta_dot_dpdn(a,j,c,i,j,k) * half_rdpdn_diffT2
                end do

                ! a = i, b /= j, c = * (skips zeros in dvel_eta_dot_dpdn)
                ! ------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle

                   cidx = cidx + 1              
                   Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)

                   Vals(cidx) = du_eta_dot_dpdn(i,b,c,i,j,k+1) * half_rdpdn_diffT1 &
                        + du_eta_dot_dpdn(i,b,c,i,j,k) * half_rdpdn_diffT2

                   cidx = cidx + 1              
                   Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)

                   Vals(cidx) = dv_eta_dot_dpdn(i,b,c,i,j,k+1) * half_rdpdn_diffT1 &
                        + dv_eta_dot_dpdn(i,b,c,i,j,k) * half_rdpdn_diffT2
                end do

                ! a = i, b = j
                ! ------------------------------------------------------------
                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)

                Vals(cidx) = du_eta_dot_dpdn(i,j,c,i,j,k+1) * half_rdpdn_diffT1 &
                     + du_eta_dot_dpdn(i,j,c,i,j,k) * half_rdpdn_diffT2

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)

                Vals(cidx) = dv_eta_dot_dpdn(i,j,c,i,j,k+1) * half_rdpdn_diffT1 &
                     + dv_eta_dot_dpdn(i,j,c,i,j,k) * half_rdpdn_diffT2
             end do

          end do
       end do
    end do

    ! ---------------------------------------------------------------------------
    ! Level nlev Jacobian Entries
    ! ---------------------------------------------------------------------------
    k = nlev
    do j=1,np
       do i=1,np

          ! get global row index for Ps equation
          ridx = ridx + 1            
          Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

          half_rdpdn_diffT2 = massmatrix(i,j) * half_rdpdn(i,j,k) * &
               (fptr%base(ie)%state%T(i,j,k,np1) &
               - fptr%base(ie)%state%T(i,j,k-1,np1))

          do c = 1,nlev

             ! a /= i, b = j (skip zeros in dvel_eta_dot_dpdn)
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)
                
                Vals(cidx) = du_eta_dot_dpdn(a,j,c,i,j,k) * half_rdpdn_diffT2

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)

                Vals(cidx) = dv_eta_dot_dpdn(a,j,c,i,j,k) * half_rdpdn_diffT2
             end do

             ! a = i, b /= j (skips zeros in dvel_eta_dot_dpdn)
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)

                Vals(cidx) = du_eta_dot_dpdn(i,b,c,i,j,k) * half_rdpdn_diffT2

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)

                Vals(cidx) = dv_eta_dot_dpdn(i,b,c,i,j,k) * half_rdpdn_diffT2
             end do

             ! a = i, b = j, c = *
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)

             Vals(cidx) = du_eta_dot_dpdn(i,j,c,i,j,k) * half_rdpdn_diffT2

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)

             Vals(cidx) = dv_eta_dot_dpdn(i,j,c,i,j,k) * half_rdpdn_diffT2
          end do

       end do
    end do

    call t_stopf('Jac_Fblock_VertAdv')

  end subroutine Jac_Fblock_VertAdv



  subroutine Jac_Fblock_PVertVel(ie, nnz_max, row_nnz, Rows, Cols, Vals, &
       c_ptr_to_object) bind(C,name="Jac_Fblock_PVertVel")
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the pressure vertical
    ! velocity term in the temperature equation
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth, kappa
    use prim_si_mod, only           : preq_omega_ps
    use derivative_mod, only        : gradient_sphere, divergence_sphere

    implicit none

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), dimension(np*np*nlev), intent(out) :: row_nnz 

    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows ! row indices 
    integer(c_int), dimension(nnz_max), intent(out)    :: Cols ! column indices
    real(c_double), dimension(nnz_max), intent(out)    :: Vals ! matrix entries
    
    type(c_ptr) :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    ! velocity derivatives of divergence(dpdn * velocity)
    real(kind=real_kind), dimension(np,np,nlev,np,np,nlev) :: du_div_dpdn_vel ! a,b,c,i,j,k
    real(kind=real_kind), dimension(np,np,nlev,np,np,nlev) :: dv_div_dpdn_vel ! a,b,c,i,j,k

    real(kind=real_kind), dimension(np,np,2,nlev) :: grad_p ! gradient of pressure

    real(kind=real_kind), dimension(np,np,nlev) :: p          ! pressure
    real(kind=real_kind), dimension(np,np,nlev) :: T_v        ! virtural temperature
    real(kind=real_kind), dimension(np,np,nlev) :: kappa_star ! variable of convenience

    real(kind=real_kind), dimension(np,np,2) :: grad_ps ! gradient of surface pressure


    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    real(kind=real_kind) :: temp1

    integer :: np1         ! current time level
    integer :: i, j, k     ! local equation index
    integer :: a, b, c     ! local derivative index
    integer :: ridx, cidx  ! index of Ps output arrays
    integer :: qn0         ! dry or moist flag
    !----------------------------------------------------------------------------

    call t_startf('Jac_Fblock_PVertVel')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0

    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! initialize number of nonzeros per row (varies)
    row_nnz = 0

    ! ---------------------------------------------------------------------------
    ! precompute quantities needed for the Jacobian entries
    ! ---------------------------------------------------------------------------

    ! time level where current iteration is stored
    np1 = fptr%tl%np1

    qn0 = fptr%n_Q

    ! compute omega = v grad p  - integral_etatop^eta[ divdp ]
    ! ------------------------------------------------------------------------

    ! gradient of surface pressure
    grad_ps = gradient_sphere(fptr%base(ie)%state%ps_v(:,:,np1), fptr%deriv, fptr%base(ie)%Dinv)
    
    do k = 1,nlev
       ! pressure
       p(:,:,k) = fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybm(k) * fptr%base(ie)%state%ps_v(:,:,np1)
       
       ! gradient of pressure
       grad_p(:,:,:,k) = fptr%hvcoord%hybm(k)*grad_ps(:,:,:)
    end do
        
    ! virtual temperature
    ! ------------------------------------------------------------------------
    if (qn0 == -1 ) then
       T_v(:,:,:) = fptr%base(ie)%state%T(:,:,:,np1)
       kappa_star = kappa
    else          
       ! Qt = fptr%base(ie)%state%Qdp(:,:,:,1,qn0)/dp(:,:,:)
       ! T_v(:,:,:) = Virtual_Temperature(fptr%base(ie)%state%T(:,:,:,np1),Qt)
       ! if (use_cpstar==1) then
       !    kappa_star = Rgas/Virtual_Specific_Heat(Qt)
       ! else
       !    kappa_star = kappa
       ! endif
       
       ! temp3(:,:,:) = Virtual_temperature(temp3, Qt)
       ! kappa_star = kappa_star * temp3
    end if
    
    ! velocity derivatives of div(dpdn v)
    call DerivVel_div_dpdn_vel(ie, fptr, du_div_dpdn_vel, dv_div_dpdn_vel)

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp
    
    ! ------------------------------------------------------------------------
    ! Level 1 
    ! ------------------------------------------------------------------------
    k = 1
    do j = 1,np
       do i = 1,np

          ! get global row index for Ps equation
          ridx = ridx + 1            
          Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)
          
          ! note the signs because pressure vertical velocity is subtracted
          ! in the resiudal equation

          temp1 = massmatrix(i,j) * kappa_star(i,j,k) * T_v(i,j,k) / p(i,j,k)

          ! a /= i, b = j, c = k
          ! ---------------------------------------------------------------------
          do a = 1,np
             if (a == i) cycle

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(k)

             ! update nnz per row
             row_nnz(ridx) = row_nnz(ridx) + 1

             Vals(cidx) = 0.5d0 * temp1 * du_div_dpdn_vel(a,j,k,i,j,k)

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(k)

             ! update nnz per row
             row_nnz(ridx) = row_nnz(ridx) + 1

             Vals(cidx) = 0.5d0 * temp1 * dv_div_dpdn_vel(a,j,k,i,j,k)
          end do

          ! a = i, b /= j, c = k
          ! ---------------------------------------------------------------------
          do b = 1,np
             if (b == j) cycle

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(k)

             ! update nnz per row
             row_nnz(ridx) = row_nnz(ridx) + 1

             Vals(cidx) = 0.5d0 * temp1 * du_div_dpdn_vel(i,b,k,i,j,k)

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(k)

             ! update nnz per row
             row_nnz(ridx) = row_nnz(ridx) + 1

             Vals(cidx) = 0.5d0 * temp1 * dv_div_dpdn_vel(i,b,k,i,j,k)
          end do

          ! a = i, b = j, c = k
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)

          ! update nnz per row
          row_nnz(ridx) = row_nnz(ridx) + 1

          Vals(cidx) = -1.0d0 * temp1 * &
               (grad_p(i,j,1,k) - 0.5d0 * du_div_dpdn_vel(i,j,k,i,j,k))

          cidx = cidx + 1              
          Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)

          ! update nnz per row
          row_nnz(ridx) = row_nnz(ridx) + 1

          Vals(cidx) = -1.0d0 * temp1 * &
               (grad_p(i,j,2,k) - 0.5d0 * dv_div_dpdn_vel(i,j,k,i,j,k))
       end do
    end do

    ! ------------------------------------------------------------------------
    ! Level 2 to nlev
    ! ------------------------------------------------------------------------
    do k = 2,nlev
       do j = 1,np
          do i = 1,np              

             ! get global row index for Ps equation
             ridx = ridx + 1            
             Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)
             
             temp1 = massmatrix(i,j) * kappa_star(i,j,k) * T_v(i,j,k) / p(i,j,k)

             ! case: c < k
             do c = 1,k-1

                ! a /= i, b = j, c < k
                ! ---------------------------------------------------------------
                do a = 1,np
                   if (a == i) cycle

                   cidx = cidx + 1              
                   Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)

                   ! update nnz per row
                   row_nnz(ridx) = row_nnz(ridx) + 1
                   
                   Vals(cidx) = temp1 * du_div_dpdn_vel(a,j,c,i,j,c) ! note k = c

                   cidx = cidx + 1              
                   Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)

                   ! update nnz per row
                   row_nnz(ridx) = row_nnz(ridx) + 1

                   Vals(cidx) = temp1 * dv_div_dpdn_vel(a,j,c,i,j,c) ! note k = c
                end do

                ! a = i, b /= j, c < k
                ! ---------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle

                   cidx = cidx + 1              
                   Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)

                   ! update nnz per row
                   row_nnz(ridx) = row_nnz(ridx) + 1

                   Vals(cidx) = temp1 * du_div_dpdn_vel(i,b,c,i,j,c) ! note k = c

                   cidx = cidx + 1              
                   Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)

                   ! update nnz per row
                   row_nnz(ridx) = row_nnz(ridx) + 1

                   Vals(cidx) = temp1 * dv_div_dpdn_vel(i,b,c,i,j,c) ! note k = c
                end do

                ! a = i, b = j, c < k
                ! ---------------------------------------------------------------
                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1

                Vals(cidx) = temp1 * du_div_dpdn_vel(i,j,c,i,j,c) ! note k = c

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1

                Vals(cidx) = temp1 * dv_div_dpdn_vel(i,j,c,i,j,c) ! note k = c

             end do

             ! a /= i, b = j, c = k
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(k)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1

                Vals(cidx) = 0.5d0 * temp1 * du_div_dpdn_vel(a,j,k,i,j,k)

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(k)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1
                
                Vals(cidx) = 0.5d0 * temp1 * dv_div_dpdn_vel(a,j,k,i,j,k)
             end do

             ! a = i, b /= j, c = k
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(k)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1

                Vals(cidx) = 0.5d0 * temp1 * du_div_dpdn_vel(i,b,k,i,j,k)

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(k)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1

                Vals(cidx) = 0.5d0 * temp1 * dv_div_dpdn_vel(i,b,k,i,j,k)
             end do

             ! a = i, b = j, c = k
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)

             ! update nnz per row
             row_nnz(ridx) = row_nnz(ridx) + 1

             Vals(cidx) = -1.0d0 * temp1 * &
                  (grad_p(i,j,1,k) - 0.5d0 * du_div_dpdn_vel(i,j,k,i,j,k))

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)

             ! update nnz per row
             row_nnz(ridx) = row_nnz(ridx) + 1

             Vals(cidx) = -1.0d0 * temp1 * &
                  (grad_p(i,j,2,k) - 0.5d0 * dv_div_dpdn_vel(i,j,k,i,j,k))
          end do
       end do
    end do

    call t_stopf('Jac_Fblock_PVertVel')

  end subroutine Jac_Fblock_PVertVel



  ! =============================================================================
  ! G Block of the Jacobian
  ! =============================================================================



  subroutine Jac_Gblock_VertAdv(ie, nnz_max, row_nnz, Rows, Cols, Vals, &
       c_ptr_to_object) bind(C,name="Jac_Gblock_VertAdv")
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the vertical advection 
    ! terms in the velocity and temperature equations
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use derivative_mod, only        : divergence_sphere

    implicit none

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows ! row indices 
    integer(c_int), dimension(nnz_max), intent(out)    :: Cols ! column indices
    real(c_double), dimension(nnz_max), intent(out)    :: Vals ! matrix entries
    
    type(c_ptr) :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! ------------------------------Local workspace------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,np,np,nlev+1) :: dps_eta_dot_dpdn

    real(kind=real_kind), dimension(np,np,nlev+1) :: eta_dot_dpdn

    real(kind=real_kind), dimension(np,np,nlev) :: dpdn
    real(kind=real_kind), dimension(np,np,nlev) :: rdpdn
    real(kind=real_kind), dimension(np,np,nlev) :: half_rdpdn
    real(kind=real_kind), dimension(np,np,nlev) :: div_dpdn_vel

    real(kind=real_kind), dimension(np,np,2) :: dpdn_vel   

    real(kind=real_kind), dimension(np,np) :: sdot_sum
    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    real(kind=real_kind) :: half_rdpdn_diffT1, half_rdpdn_diffT2
    real(kind=real_kind) :: temp1, temp2, diffB

    integer :: np1         ! current time level
    integer :: i, j, k     ! local equation index
    integer :: a, b        ! local derivative index
    integer :: ridx, cidx  ! index of Ps output arrays
    ! ---------------------------------------------------------------------------

    call t_startf('Jac_Gblock_VertAdv')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0

    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = np+np-1

    ! ---------------------------------------------------------------------------
    ! precompute quantities needed for the Jacobian entries
    ! ---------------------------------------------------------------------------

    ! time level where current iteration is stored
    np1 = fptr%tl%np1

    ! vertical derivative of pressure, dpdn
    do k = 1,nlev
       dpdn(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
            - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))

       rdpdn(:,:,k) = 1.0d0 / dpdn(:,:,k)

       half_rdpdn(:,:,k) = 0.5d0 / dpdn(:,:,k)
    end do

    ! divergence of dpdn * velocity, needed for eta_dot_dpdn
    do k = 1,nlev
       do j=1,np
          do i=1,np             
             dpdn_vel(i,j,1) = dpdn(i,j,k) * fptr%base(ie)%state%v(i,j,1,k,np1)
             dpdn_vel(i,j,2) = dpdn(i,j,k) * fptr%base(ie)%state%v(i,j,2,k,np1)
          end do
       end do

       div_dpdn_vel(:,:,k) = divergence_sphere(dpdn_vel, fptr%deriv, fptr%base(ie))
    end do

    ! compute eta_dot_dpdn
    ! ---------------------------------------------------------------------------
    sdot_sum = 0.0d0

    do k=1,nlev
       sdot_sum(:,:)         = sdot_sum(:,:) + div_dpdn_vel(:,:,k)
       eta_dot_dpdn(:,:,k+1) = sdot_sum(:,:)
    end do

    do k=1,nlev-1
       eta_dot_dpdn(:,:,k+1) = fptr%hvcoord%hybi(k+1)*sdot_sum(:,:) - eta_dot_dpdn(:,:,k+1)
    end do

    eta_dot_dpdn(:,:,1)      = 0.0d0
    eta_dot_dpdn(:,:,nlev+1) = 0.0d0

    ! compute the Jacobian of etadot_dpdn
    call DerivPs_etadot_dpdn(ie, fptr, dps_eta_dot_dpdn) 

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp

    ! ---------------------------------------------------------------------------
    ! Level 1 Jacobian entries, surface pressure derivatives
    ! ---------------------------------------------------------------------------
    k = 1
    do j = 1,np
       do i = 1,np

          ! get global row index for Ps equation
          ridx = ridx + 1            
          Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

          half_rdpdn_diffT1 = massmatrix(i,j) * half_rdpdn(i,j,k) * &
               (fptr%base(ie)%state%T(i,j,k+1,np1) &
               - fptr%base(ie)%state%T(i,j,k,np1))

          ! a /= i, b = j (skips zeros in dvel_eta_dot_dpdn)
          ! ---------------------------------------------------------------------
          do a = 1,np
             if (a == i) cycle

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffset
             
             Vals(cidx) = dps_eta_dot_dpdn(a,j,i,j,k+1) * half_rdpdn_diffT1
          end do

          ! a = i, b /= j (skips zeros in dvel_eta_dot_dpdn)
          ! ---------------------------------------------------------------------
          do b = 1,np
             if (b == j) cycle

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffset

             Vals(cidx) = dps_eta_dot_dpdn(i,b,i,j,k+1) * half_rdpdn_diffT1
          end do

          ! a = i, b = j
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Cols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffset

          diffB = fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)

          temp1 = massmatrix(i,j) * &
               (-0.5d0) * diffB * rdpdn(i,j,k)**2 * eta_dot_dpdn(i,j,k+1)

          Vals(cidx) = dps_eta_dot_dpdn(i,j,i,j,k+1) * half_rdpdn_diffT1 &
               + temp1 * (fptr%base(ie)%state%T(i,j,k+1,np1) &
               - fptr%base(ie)%state%T(i,j,k,np1))
       end do
    end do

    ! ------------------------------------------------------------------------
    ! Level 2 to nlev-1 Jacobian Entries, surface pressure derivatives
    ! ------------------------------------------------------------------------
    do k = 2,nlev-1
       do j = 1,np
          do i = 1,np

             ! get global row index for Ps equation
             ridx = ridx + 1            
             Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

             half_rdpdn_diffT1 = massmatrix(i,j) * half_rdpdn(i,j,k) * &
                  (fptr%base(ie)%state%T(i,j,k+1,np1) &
                  - fptr%base(ie)%state%T(i,j,k,np1))

             half_rdpdn_diffT2 = massmatrix(i,j) * half_rdpdn(i,j,k) * &
                  (fptr%base(ie)%state%T(i,j,k,np1) &
                  - fptr%base(ie)%state%T(i,j,k-1,np1))

             ! a /= i, b = j (skips zeros in dvel_eta_dot_dpdn)
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffset

                Vals(cidx) = dps_eta_dot_dpdn(a,j,i,j,k+1) * half_rdpdn_diffT1 &
                     + dps_eta_dot_dpdn(a,j,i,j,k) * half_rdpdn_diffT2
             end do

             ! a = i, b /= j (skips zeros in dvel_eta_dot_dpdn)
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffset

                Vals(cidx) = dps_eta_dot_dpdn(i,b,i,j,k+1) * half_rdpdn_diffT1 &
                     + dps_eta_dot_dpdn(i,b,i,j,k) * half_rdpdn_diffT2
             end do

             ! a = i and b = j
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffset

             diffB = fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)

             temp1 = massmatrix(i,j) * (-0.5d0) * diffB * rdpdn(i,j,k)**2 * eta_dot_dpdn(i,j,k+1) ! eta_dot_dpdn k+1/2
             temp2 = massmatrix(i,j) * (-0.5d0) * diffB * rdpdn(i,j,k)**2 * eta_dot_dpdn(i,j,k)   ! eta_dot_dpdn k-1/2

             Vals(cidx) = dps_eta_dot_dpdn(i,j,i,j,k+1) * half_rdpdn_diffT1 &
                  + dps_eta_dot_dpdn(i,j,i,j,k) * half_rdpdn_diffT2 &
                  + temp1 * (fptr%base(ie)%state%T(i,j,k+1,np1) - fptr%base(ie)%state%T(i,j,k,np1)) &
                  + temp2 * (fptr%base(ie)%state%T(i,j,k,np1)   - fptr%base(ie)%state%T(i,j,k-1,np1))
          end do
       end do
    end do

    ! ------------------------------------------------------------------------
    ! Level nlev Jacobian Entries, surface pressure derivatives
    ! ------------------------------------------------------------------------
    k = nlev
    do j=1,np
       do i=1,np

          ! get global row index for Ps equation
          ridx = ridx + 1            
          Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

          half_rdpdn_diffT2 = massmatrix(i,j) * half_rdpdn(i,j,k) * &
               (fptr%base(ie)%state%T(i,j,k,np1) &
               - fptr%base(ie)%state%T(i,j,k-1,np1))

          ! a \= i, b = j (skips zeros in dvel_eta_dot_dpdn)
          ! ---------------------------------------------------------------------
          do a = 1,np
             if (a == i) cycle

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffset

             Vals(cidx) = dps_eta_dot_dpdn(a,j,i,j,k) * half_rdpdn_diffT2
          end do

          ! a = i, b \= j (skips zeros in dvel_eta_dot_dpdn)
          ! ---------------------------------------------------------------------
          do b = 1,np
             if (b == j) cycle

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffset

             Vals(cidx) = dps_eta_dot_dpdn(i,b,i,j,k) * half_rdpdn_diffT2
          end do

          ! a = i, b = j 
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Cols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffset

          diffB = fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)

          temp2 = massmatrix(i,j) * (-0.5d0) * diffB * rdpdn(i,j,k)**2 * eta_dot_dpdn(i,j,k) ! eta_dot_dpdn k-1/2

          Vals(cidx) = dps_eta_dot_dpdn(i,j,i,j,k) * half_rdpdn_diffT2 &
               + temp2 * (fptr%base(ie)%state%T(i,j,k,np1) &
               - fptr%base(ie)%state%T(i,j,k-1,np1))
       end do
    end do

    call t_stopf('Jac_Gblock_VertAdv')

  end subroutine Jac_Gblock_VertAdv



  subroutine Jac_Gblock_PVertVel(ie, nnz_max, row_nnz, Rows, Cols, Vals, &
       c_ptr_to_object) bind(C,name="Jac_Gblock_PVertVel")
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the pressure vertical
    ! velocity term in the temperature equation
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth, kappa
    use prim_si_mod, only           : preq_omega_ps
    use derivative_mod, only        : gradient_sphere, divergence_sphere

    implicit none

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows ! row indices 
    integer(c_int), dimension(nnz_max), intent(out)    :: Cols ! column indices
    real(c_double), dimension(nnz_max), intent(out)    :: Vals ! matrix entries
    
    type(c_ptr) :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()
  
    real(kind=real_kind), dimension(np,np,np,np,nlev) :: dps_div_dpdn_vel ! a,b,-,i,j,k

    real(kind=real_kind), dimension(np,np,2,nlev) :: grad_p ! gradient of pressure

    real(kind=real_kind), dimension(np,np,nlev) :: p          ! pressure
    real(kind=real_kind), dimension(np,np,nlev) :: dp         ! vertical derivative of pressure
    real(kind=real_kind), dimension(np,np,nlev) :: vgrad_p    ! v dot gradient(p)
    real(kind=real_kind), dimension(np,np,nlev) :: divdp      ! divergence(v*dp)
    real(kind=real_kind), dimension(np,np,nlev) :: T_v        ! virtural temperature
    real(kind=real_kind), dimension(np,np,nlev) :: omega_p    ! pressure vertical velocity / pressure
    real(kind=real_kind), dimension(np,np,nlev) :: kappa_star ! variable of convenience

    real(kind=real_kind), dimension(np,np,2) :: grad_ps ! gradient of surface pressure
    real(kind=real_kind), dimension(np,np,2) :: vtemp   ! velocity * dp

    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    real(kind=real_kind) :: v1, v2, temp1, temp2, temp3, temp_sum

    integer :: np1         ! current time level
    integer :: i, j, k     ! local equation index
    integer :: a, b, c     ! local derivative index
    integer :: ridx, cidx  ! index of Ps output arrays
    integer :: qn0         ! dry or moist flag
    !----------------------------------------------------------------------------

    call t_startf('Jac_Gblock_PVertVel')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0

    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = np+np-1

    ! ---------------------------------------------------------------------------
    ! precompute quantities needed for the Jacobian entries
    ! ---------------------------------------------------------------------------

    ! time level where current iteration is stored
    np1 = fptr%tl%np1

    qn0 = fptr%n_Q

    ! compute omega = v grad p  - integral_etatop^eta[ divdp ]
    ! ------------------------------------------------------------------------

    ! gradient of surface pressure
    grad_ps = gradient_sphere(fptr%base(ie)%state%ps_v(:,:,np1), fptr%deriv, fptr%base(ie)%Dinv)
    
    do k = 1,nlev
       ! vertical derivative of pressure
       dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
            - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))
       
       ! pressure
       p(:,:,k) = fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybm(k) * fptr%base(ie)%state%ps_v(:,:,np1)
       
       ! gradient of pressure
       grad_p(:,:,:,k) = fptr%hvcoord%hybm(k)*grad_ps(:,:,:)
       
       do j = 1,np
          do i = 1,np
             v1 = fptr%base(ie)%state%v(i,j,1,k,np1)
             v2 = fptr%base(ie)%state%v(i,j,2,k,np1)
             
             vgrad_p(i,j,k) = (v1*grad_p(i,j,1,k) + v2*grad_p(i,j,2,k))
             vtemp(i,j,1) = v1*dp(i,j,k)
             vtemp(i,j,2) = v2*dp(i,j,k)
          end do
       end do
       
       divdp(:,:,k) = divergence_sphere(vtemp, fptr%deriv, fptr%base(ie))
       
    end do
    
    ! compute pressure vertical velocity / pressure (omega/p)
    call preq_omega_ps(omega_p, fptr%hvcoord, p, vgrad_p, divdp)
    
    ! virtual temperature
    if (qn0 == -1 ) then
       T_v(:,:,:) = fptr%base(ie)%state%T(:,:,:,np1)
       kappa_star = kappa
    else          
       ! Qt = fptr%base(ie)%state%Qdp(:,:,:,1,qn0)/dp(:,:,:)
       ! T_v(:,:,:) = Virtual_Temperature(fptr%base(ie)%state%T(:,:,:,np1),Qt)
       ! if (use_cpstar==1) then
       !    kappa_star = Rgas/Virtual_Specific_Heat(Qt)
       ! else
       !    kappa_star = kappa
       ! endif
       
       ! temp3(:,:,:) = Virtual_temperature(temp3, Qt)
       ! kappa_star = kappa_star * temp3
    end if
    
    ! compute derivatives div(dpdn v)
    call DerivPs_div_dpdn_vel(ie, fptr, dps_div_dpdn_vel)

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp
    
    ! ------------------------------------------------------------------------
    ! Level 1 Jacobian entries 
    ! ------------------------------------------------------------------------
    k = 1
    do j = 1,np
       do i = 1,np

          ! get global row index for Ps equation
          ridx = ridx + 1            
          Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

          ! note the signs because pressure vertical velocity is subtracted
          ! in the resiudal equation

          temp1 = massmatrix(i,j) * kappa_star(i,j,k) * T_v(i,j,k) / p(i,j,k)

          temp2 = rrearth * fptr%hvcoord%hybm(k) &
               * (fptr%base(ie)%state%v(i,j,1,k,np1) * fptr%base(ie)%Dinv(i,j,1,1) &
               + fptr%base(ie)%state%v(i,j,2,k,np1) * fptr%base(ie)%Dinv(i,j,1,2))
          
          temp3 = rrearth * fptr%hvcoord%hybm(k) &
               * (fptr%base(ie)%state%v(i,j,1,k,np1) * fptr%base(ie)%Dinv(i,j,2,1) &
               + fptr%base(ie)%state%v(i,j,2,k,np1) * fptr%base(ie)%Dinv(i,j,2,2))

          ! a /= i, b = j
          ! ---------------------------------------------------------------------
          do a = 1,np
             if (a == i) cycle

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffset

             Vals(cidx) = -1.0d0 * temp1 * (temp2 * fptr%deriv%Dvv(a,i) &
                  - 0.5d0 * dps_div_dpdn_vel(a,j,i,j,k))
          end do

          ! a = i, b /= j
          ! ---------------------------------------------------------------------
          do b = 1,np
             if (b == j) cycle

             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffset

             Vals(cidx) = -1.0d0 * temp1 * (temp3 * fptr%deriv%Dvv(b,j) &
                  - 0.5d0 * dps_div_dpdn_vel(i,b,i,j,k))
          end do

          ! a = i, b = j
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Cols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffset

          Vals(cidx) = temp1 * fptr%hvcoord%hybm(k) * omega_p(i,j,k) &
               - temp1 * (temp2 * fptr%deriv%Dvv(i,i) + temp3 * fptr%deriv%Dvv(j,j) &
               - 0.5d0 * dps_div_dpdn_vel(i,j,i,j,k))

       end do
    end do

    ! ------------------------------------------------------------------------
    ! Level 2 to nlev Jacobian entries
    ! ------------------------------------------------------------------------
    do k = 2,nlev
       do j = 1,np
          do i = 1,np              

             ! get global row index for Ps equation
             ridx = ridx + 1            
             Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)
             
             temp1 = massmatrix(i,j) * kappa_star(i,j,k) * T_v(i,j,k) / p(i,j,k)

             temp2 = rrearth * fptr%hvcoord%hybm(k) &
                  * (fptr%base(ie)%state%v(i,j,1,k,np1) * fptr%base(ie)%Dinv(i,j,1,1) &
                  + fptr%base(ie)%state%v(i,j,2,k,np1) * fptr%base(ie)%Dinv(i,j,1,2))
             
             temp3 = rrearth * fptr%hvcoord%hybm(k) &
                  * (fptr%base(ie)%state%v(i,j,1,k,np1) * fptr%base(ie)%Dinv(i,j,2,1) &
                  + fptr%base(ie)%state%v(i,j,2,k,np1) * fptr%base(ie)%Dinv(i,j,2,2))
             
             ! a /= i, b = j
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffset

                temp_sum = 0.0d0
                do c = 1,k-1
                   temp_sum = temp_sum + dps_div_dpdn_vel(a,j,i,j,c)
                end do

                Vals(cidx) = -1.0d0 * temp1 * (temp2 * fptr%deriv%Dvv(a,i) &
                     - temp_sum - 0.5d0 * dps_div_dpdn_vel(a,j,i,j,k))
             end do

             ! a = i, b /= j
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffset

                temp_sum = 0.0d0
                do c = 1,k-1
                   temp_sum = temp_sum + dps_div_dpdn_vel(i,b,i,j,c)
                end do

                Vals(cidx) = -1.0d0 * temp1 * (temp3 * fptr%deriv%Dvv(b,j) &
                     - temp_sum - 0.5d0 * dps_div_dpdn_vel(i,b,i,j,k))
             end do

             ! a = i, b = j
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffset

             temp_sum = 0.0d0
             do c = 1,k-1
                temp_sum = temp_sum + dps_div_dpdn_vel(i,j,i,j,c)
             end do

             Vals(cidx) = temp1 * fptr%hvcoord%hybm(k) * omega_p(i,j,k) &
                  - temp1 * (temp2 * fptr%deriv%Dvv(i,i) + temp3 * fptr%deriv%Dvv(j,j) &
                  - temp_sum - 0.5d0 * dps_div_dpdn_vel(i,j,i,j,k))
          end do
       end do
    end do

    call t_stopf('Jac_Gblock_PVertVel')

  end subroutine Jac_Gblock_PVertVel



  ! =============================================================================
  ! H Block of the Jacobian
  ! =============================================================================



  subroutine Jac_Hblock_Time(ie, nnz_max, row_nnz, Rows, Cols, Vals, &
       c_ptr_to_object) bind(C,name="Jac_Hblock_Time")
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the pressure vertical
    ! velocity term in the temperature equation
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use control_mod, only           : tstep_type
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type

    implicit none

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows ! row indices 
    integer(c_int), dimension(nnz_max), intent(out)    :: Cols ! column indices
    real(c_double), dimension(nnz_max), intent(out)    :: Vals ! matrix entries
    
    type(c_ptr) :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling
    
    real(kind=real_kind) :: dti ! 1 / time step


    integer :: i, j, k     ! local equation index
    integer :: ridx, cidx  ! index of Ps output arrays
    !----------------------------------------------------------------------------

    call t_startf('Jac_Hblock_Time')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0

    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = 1 

    ! inverse of time step size
    dti = 1.0d0 / fptr%dt

    if (tstep_type==12) then
       dti = (3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2 on step size
    endif

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp
        
    ! ------------------------------------------------------------------------
    ! compute Jacobian entries, temperature derivative
    ! ------------------------------------------------------------------------
    do k = 1,nlev
       do j = 1,np
          do i = 1,np              

             ! get global row index for temperature equation
             ridx = ridx + 1            
             Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)            

             ! a = i, b = j, c = k               
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

             Vals(cidx) = massmatrix(i,j) * dti

          end do
       end do
    end do

    call t_stopf('Jac_Hblock_Time')

  end subroutine Jac_Hblock_Time



  subroutine Jac_Hblock_HorizAdv(ie, nnz_max, row_nnz, Rows, Cols, Vals, &
       c_ptr_to_object) bind(C,name="Jac_Hblock_HorizAdv")
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the horizontal advection 
    ! terms in the temperature equation
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth
    use derivative_mod, only        : gradient_sphere

    implicit none

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows ! row indices 
    integer(c_int), dimension(nnz_max), intent(out)    :: Cols ! column indices
    real(c_double), dimension(nnz_max), intent(out)    :: Vals ! matrix entries
    
    type(c_ptr) :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------  
    type(derived_type), pointer :: fptr=>NULL()
  
    real(kind=real_kind), dimension(np,np,2) :: grad_T ! gradient of temperature

    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    real(kind=real_kind) :: temp1, temp2

    integer :: np1         ! current time level
    integer :: i, j, k     ! local equation index
    integer :: a, b        ! local derivative index
    integer :: ridx, cidx  ! index of Ps output arrays
    !----------------------------------------------------------------------------

    call t_startf('Jac_Hblock_HorizAdv')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0

    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = np+np-1

    ! time level where current iteration is stored
    np1 = fptr%tl%np1

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp

    ! ---------------------------------------------------------------------------
    ! compute Jacobian entires (temperature derivatives)
    ! ---------------------------------------------------------------------------
    do k = 1,nlev

       ! gradient of temperature
       grad_T = gradient_sphere(fptr%base(ie)%state%T(:,:,k,np1), fptr%deriv, &
            fptr%base(ie)%Dinv)
       
       do j = 1,np
          do i = 1,np

             ! get global row index for Ps equation
             ridx = ridx + 1            
             Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)
                          
             temp1 = massmatrix(i,j) * rrearth * &
                  (fptr%base(ie)%state%v(i,j,1,k,np1) * fptr%base(ie)%Dinv(i,j,1,1) &
                  + fptr%base(ie)%state%v(i,j,2,k,np1) * fptr%base(ie)%Dinv(i,j,1,2))
             
             temp2 = massmatrix(i,j) * rrearth * &
                  (fptr%base(ie)%state%v(i,j,1,k,np1) * fptr%base(ie)%Dinv(i,j,2,1) &
                  + fptr%base(ie)%state%v(i,j,2,k,np1) * fptr%base(ie)%Dinv(i,j,2,2))
             
             ! a /= i, b = j, c = k (x-derivative terms)
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Toffset(k)

                Vals(cidx) = temp1 * fptr%deriv%Dvv(a,i)
             end do
             
             ! a = i, b /= j, c = k (y-derivative terms)
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle                  

                cidx = cidx + 1              
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Toffset(k)

                Vals(cidx) = temp2 * fptr%deriv%Dvv(b,j)
             end do
             
             ! a = i, b = j, c = k (x and y-derivative terms)
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

             Vals(cidx) = temp1 * fptr%deriv%Dvv(i,i) + temp2 * fptr%deriv%Dvv(j,j)
             
          end do
       end do
    end do
    
    call t_stopf('Jac_Hblock_HorizAdv')

  end subroutine Jac_Hblock_HorizAdv



  subroutine Jac_Hblock_VertAdv(ie, nnz_max, row_nnz, Rows, Cols, Vals, &
       c_ptr_to_object) bind(C,name="Jac_Hblock_VertAdv")
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the vertical advection 
    ! terms in the velocity and temperature equations
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use derivative_mod, only        : divergence_sphere

    implicit none

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), dimension(np*np*nlev), intent(out) :: row_nnz 

    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows ! row indices 
    integer(c_int), dimension(nnz_max), intent(out)    :: Cols ! column indices
    real(c_double), dimension(nnz_max), intent(out)    :: Vals ! matrix entries
    
    type(c_ptr) :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,nlev+1) :: eta_dot_dpdn

    real(kind=real_kind), dimension(np,np,nlev) :: dpdn
    real(kind=real_kind), dimension(np,np,nlev) :: half_rdpdn
    real(kind=real_kind), dimension(np,np,nlev) :: div_dpdn_vel

    real(kind=real_kind), dimension(np,np,2) :: dpdn_vel   

    real(kind=real_kind), dimension(np,np) :: sdot_sum
    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling


    real(kind=real_kind) :: temp1, temp2

    integer :: np1         ! current time level
    integer :: i, j, k     ! local equation index

    integer :: ridx, cidx  ! index of Ps output arrays
    !----------------------------------------------------------------------------

    call t_startf('Jac_Hblock_VertAdv')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0

    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! initialize number of nonzeros per row (varies: 2,3,2)
    row_nnz = 0 

    ! ---------------------------------------------------------------------------
    ! precompute quantities needed for the Jacobian entries
    ! ---------------------------------------------------------------------------

    ! time level where current iteration is stored
    np1 = fptr%tl%np1

    ! vertical derivative of pressure, dpdn
    do k = 1,nlev
       dpdn(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
            - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))

       half_rdpdn(:,:,k) = 0.5d0 / dpdn(:,:,k)
    end do

    ! divergence of dpdn * velocity, needed for eta_dot_dpdn
    do k = 1,nlev
       do j=1,np
          do i=1,np             
             dpdn_vel(i,j,1) = dpdn(i,j,k) * fptr%base(ie)%state%v(i,j,1,k,np1)
             dpdn_vel(i,j,2) = dpdn(i,j,k) * fptr%base(ie)%state%v(i,j,2,k,np1)
          end do
       end do

       div_dpdn_vel(:,:,k) = divergence_sphere(dpdn_vel, fptr%deriv, fptr%base(ie))
    end do

    ! compute eta_dot_dpdn
    sdot_sum = 0.0d0

    do k=1,nlev
       sdot_sum(:,:)         = sdot_sum(:,:) + div_dpdn_vel(:,:,k)
       eta_dot_dpdn(:,:,k+1) = sdot_sum(:,:)
    end do

    do k=1,nlev-1
       eta_dot_dpdn(:,:,k+1) = fptr%hvcoord%hybi(k+1)*sdot_sum(:,:) - eta_dot_dpdn(:,:,k+1)
    end do

    eta_dot_dpdn(:,:,1)      = 0.0d0
    eta_dot_dpdn(:,:,nlev+1) = 0.0d0

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp

    ! ---------------------------------------------------------------------------
    ! Level 1 Jacobian entries (temperature derivative)
    ! ---------------------------------------------------------------------------
    k = 1
    do j = 1,np
       do i = 1,np

          ! get global row index for Ps equation
          ridx = ridx + 1            
          Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

          temp1 = massmatrix(i,j) * half_rdpdn(i,j,k) * eta_dot_dpdn(i,j,k+1) ! eta_dot_dpdn k+1/2

          ! a = i, b = j, c = k
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

          ! update nnz per row
          row_nnz(ridx) = row_nnz(ridx) + 1

          Vals(cidx) = -1.0d0 * temp1

          ! a = i, b = j, c = k+1
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k+1)

          ! update nnz per row
          row_nnz(ridx) = row_nnz(ridx) + 1

          Vals(cidx) = temp1

       end do
    end do

    ! ------------------------------------------------------------------------
    ! Level 2 to nlev-1 Jacobian Entries, temperature derivative
    ! ------------------------------------------------------------------------
    do k = 2,nlev-1
       do j = 1,np
          do i = 1,np

             ! get global row index for Ps equation
             ridx = ridx + 1            
             Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

             temp1 = massmatrix(i,j) * half_rdpdn(i,j,k) * eta_dot_dpdn(i,j,k+1) ! eta_dot_dpdn k+1/2
             temp2 = massmatrix(i,j) * half_rdpdn(i,j,k) * eta_dot_dpdn(i,j,k)   ! eta_dot_dpdn k-1/2              

             ! a = i, b = j, c = k-1
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k-1)

             ! update nnz per row
             row_nnz(ridx) = row_nnz(ridx) + 1
             
             Vals(cidx) = -1.0d0 * temp2

             ! a = i, b = j, c = k
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

             ! update nnz per row
             row_nnz(ridx) = row_nnz(ridx) + 1

             Vals(cidx) =  temp2 - temp1

             ! a = i, b = j, c = k+1
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k+1)

             ! update nnz per row
             row_nnz(ridx) = row_nnz(ridx) + 1

             Vals(cidx) = temp1

          end do
       end do
    end do

    ! ---------------------------------------------------------------------------
    ! Level nlev Jacobian Entries, temperature derivative
    ! ---------------------------------------------------------------------------
    k = nlev
    do j=1,np
       do i=1,np

          ! get global row index for Ps equation
          ridx = ridx + 1            
          Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

          temp2 = massmatrix(i,j) * half_rdpdn(i,j,k) * eta_dot_dpdn(i,j,k)   ! eta_dot_dpdn k-1/2

          ! a = i, b = j, c = k
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

          ! update nnz per row
          row_nnz(ridx) = row_nnz(ridx) + 1

          Vals(cidx) = temp2

          ! a = i, b = j, c = k-1
          ! ---------------------------------------------------------------------
          cidx = cidx + 1              
          Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k-1)

          ! update nnz per row
          row_nnz(ridx) = row_nnz(ridx) + 1

          Vals(cidx) = -1.0d0 * temp2

       end do
    end do

    call t_stopf('Jac_Hblock_VertAdv')

  end subroutine Jac_Hblock_VertAdv



  subroutine Jac_Hblock_PVertVel(ie, nnz_max, row_nnz, Rows, Cols, Vals, &
       c_ptr_to_object) bind(C,name="Jac_Hblock_PVertVel")
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the pressure vertical
    ! velocity term in the temperature equation
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth, kappa
    use prim_si_mod, only           : preq_omega_ps
    use derivative_mod, only        : gradient_sphere, divergence_sphere

    implicit none

    ! --------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per row

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows ! row indices 
    integer(c_int), dimension(nnz_max), intent(out)    :: Cols ! column indices
    real(c_double), dimension(nnz_max), intent(out)    :: Vals ! matrix entries
    
    type(c_ptr) :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()
    
    real(kind=real_kind), dimension(np,np,2,nlev) :: grad_p ! gradient of pressure

    real(kind=real_kind), dimension(np,np,nlev) :: p          ! pressure
    real(kind=real_kind), dimension(np,np,nlev) :: dp         ! vertical derivative of pressure
    real(kind=real_kind), dimension(np,np,nlev) :: vgrad_p    ! v dot gradient(p)
    real(kind=real_kind), dimension(np,np,nlev) :: divdp      ! divergence(v*dp)
    real(kind=real_kind), dimension(np,np,nlev) :: T_v        ! virtural temperature
    real(kind=real_kind), dimension(np,np,nlev) :: omega_p    ! pressure vertical velocity / pressure
    real(kind=real_kind), dimension(np,np,nlev) :: kappa_star ! variable of convenience

    real(kind=real_kind), dimension(np,np,2) :: grad_ps ! gradient of surface pressure
    real(kind=real_kind), dimension(np,np,2) :: vtemp   ! velocity * dp

    real(kind=real_kind), dimension(np,np) :: massmatrix ! DSS scaling

    real(kind=real_kind) :: v1, v2

    integer :: np1        ! current time level
    integer :: i, j, k    ! local equation index

    integer :: ridx, cidx ! index of Ps output arrays
    integer :: qn0        ! dry or moist flag
    !----------------------------------------------------------------------------

    call t_startf('Jac_Hblock_PVertVel')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0

    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = 1 

    ! ---------------------------------------------------------------------------
    ! precompute quantities needed for the Jacobian entries
    ! ---------------------------------------------------------------------------

    ! time level where current iteration is stored
    np1 = fptr%tl%np1

    qn0 = fptr%n_Q

    ! compute omega = v grad p  - integral_etatop^eta[ divdp ]
    ! ------------------------------------------------------------------------

    ! gradient of surface pressure
    grad_ps = gradient_sphere(fptr%base(ie)%state%ps_v(:,:,np1), fptr%deriv, fptr%base(ie)%Dinv)
    
    do k = 1,nlev
       ! vertical derivative of pressure
       dp(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
            - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))
       
       ! pressure
       p(:,:,k) = fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybm(k) * fptr%base(ie)%state%ps_v(:,:,np1)
       
       ! gradient of pressure
       grad_p(:,:,:,k) = fptr%hvcoord%hybm(k)*grad_ps(:,:,:)
       
       do j = 1,np
          do i = 1,np
             v1 = fptr%base(ie)%state%v(i,j,1,k,np1)
             v2 = fptr%base(ie)%state%v(i,j,2,k,np1)
             
             vgrad_p(i,j,k) = (v1*grad_p(i,j,1,k) + v2*grad_p(i,j,2,k))
             vtemp(i,j,1) = v1*dp(i,j,k)
             vtemp(i,j,2) = v2*dp(i,j,k)
          end do
       end do
       
       divdp(:,:,k) = divergence_sphere(vtemp, fptr%deriv, fptr%base(ie))
       
    end do
    
    ! compute omega/p
    call preq_omega_ps(omega_p, fptr%hvcoord, p, vgrad_p, divdp)
    
    ! virtual temperature
    if (qn0 == -1) then
       T_v(:,:,:) = fptr%base(ie)%state%T(:,:,:,np1)
       kappa_star = kappa
    else          
       ! Qt = fptr%base(ie)%state%Qdp(:,:,:,1,qn0)/dp(:,:,:)
       ! T_v(:,:,:) = Virtual_Temperature(fptr%base(ie)%state%T(:,:,:,np1),Qt)
       ! if (use_cpstar==1) then
       !    kappa_star = Rgas/Virtual_Specific_Heat(Qt)
       ! else
       !    kappa_star = kappa
       ! endif
       
       ! temp3(:,:,:) = Virtual_temperature(temp3, Qt)
       ! kappa_star = kappa_star * temp3
    end if

    ! mass matrix scaling
    massmatrix = fptr%base(ie)%spheremp * fptr%base(ie)%rspheremp
        
    ! ------------------------------------------------------------------------
    ! compute Jacobian entries, temperature derivative
    ! ------------------------------------------------------------------------
    do k = 1,nlev
       do j = 1,np
          do i = 1,np              

             ! get global row index for temperature equation
             ridx = ridx + 1            
             Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)            

             ! a = i, b = j, c = k               
             ! ------------------------------------------------------------------
             cidx = cidx + 1              
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

             Vals(cidx) = massmatrix(i,j) * (-1.0d0) * kappa_star(i,j,k) * omega_p(i,j,k)

          end do
       end do
    end do

    call t_stopf('Jac_Hblock_PVertVel')

  end subroutine Jac_Hblock_PVertVel



  ! =============================================================================
  ! Utility subroutines for computing Jacobian entries
  ! =============================================================================



  subroutine DerivVel_div_dpdn_vel(ie, fptr, du_div_dpdn_vel, dv_div_dpdn_vel)
    !----------------------------------------------------------------------------
    ! compute U and V derivatives of divergence(dpdn * vel)
    !----------------------------------------------------------------------------
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth
    
    implicit none

    !---------------------------------Arguments----------------------------------
    integer, intent(in) :: ie 

    type(derived_type), pointer, intent(in) :: fptr

    real(kind=real_kind), intent(out) :: du_div_dpdn_vel(np,np,nlev,np,np,nlev) ! a,b,c,i,j,k
    real(kind=real_kind), intent(out) :: dv_div_dpdn_vel(np,np,nlev,np,np,nlev) ! a,b,c,i,j,k
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    real(kind=real_kind), dimension(np,np,nlev) :: dpdn
    real(kind=real_kind), dimension(nlev)       :: diffBi

    real(kind=real_kind) :: temp1, temp2, temp3

    integer :: np1     ! time level
    integer :: i, j, k ! equation index
    integer :: a, b    ! derivative index
    !----------------------------------------------------------------------------

    call t_startf('DerivVel_div_dpdn_vel')

    ! get time level
    np1 = fptr%tl%np1

    ! initialize outputs
    du_div_dpdn_vel = 0.0d0
    dv_div_dpdn_vel = 0.0d0
    
    ! compute dp/deta and d(dp/deta)/dps
    do k = 1,nlev
       dpdn(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
            - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))

       diffBi(k) = fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)
    end do

    do k = 1,nlev       
       do j = 1,np
          do i = 1,np

             temp1 = rrearth * fptr%base(ie)%rmetdet(i,j)
             
             ! a /= i, b = j, c = k (x-derivative)
             do a = 1,np
                if (a == i) cycle
                
                temp2 = temp1 * fptr%base(ie)%metdet(a,j) * dpdn(a,j,k) * fptr%deriv%Dvv(a,i)
                
                du_div_dpdn_vel(a,j,k,i,j,k) = temp2 * fptr%base(ie)%Dinv(a,j,1,1)

                dv_div_dpdn_vel(a,j,k,i,j,k) = temp2 * fptr%base(ie)%Dinv(a,j,1,2)

             end do
             
             ! a = i, b /= j, c = k (y-derivative)
             do b = 1,np
                if (b == j) cycle

                temp2 = temp1 * fptr%base(ie)%metdet(i,b) * dpdn(i,b,k) * fptr%deriv%Dvv(b,j)
                
                du_div_dpdn_vel(i,b,k,i,j,k) = temp2 * fptr%base(ie)%Dinv(i,b,2,1)

                dv_div_dpdn_vel(i,b,k,i,j,k) = temp2 * fptr%base(ie)%Dinv(i,b,2,2)

             end do

             ! a = i, b = j, c = k (x and y-derivative)
             temp2 = rrearth * dpdn(i,j,k) * fptr%deriv%Dvv(i,i)
             
             temp3 = rrearth * dpdn(i,j,k) * fptr%deriv%Dvv(j,j)
             
             du_div_dpdn_vel(i,j,k,i,j,k) = temp2 * fptr%base(ie)%Dinv(i,j,1,1) &
                  + temp3 * fptr%base(ie)%Dinv(i,j,2,1)
             
             dv_div_dpdn_vel(i,j,k,i,j,k) = temp2 * fptr%base(ie)%Dinv(i,j,1,2) &
                  + temp3 * fptr%base(ie)%Dinv(i,j,2,2)
          
          end do
       end do
    end do
   
    call t_stopf('DerivVel_div_dpdn_vel')
    
  end subroutine DerivVel_div_dpdn_vel



  subroutine DerivPs_div_dpdn_vel(ie, fptr, dps_div_dpdn_vel)
    !----------------------------------------------------------------------------
    ! compute Ps derivatives of divergence(dpdn * vel)
    !----------------------------------------------------------------------------
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth
    
    implicit none

    !---------------------------------Arguments----------------------------------
    integer, intent(in) :: ie 

    type(derived_type), pointer, intent(in) :: fptr

    real(kind=real_kind), intent(out) :: dps_div_dpdn_vel(np,np,np,np,nlev) ! a,b,c,i,j,k
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    real(kind=real_kind), dimension(np,np,nlev) :: dpdn
    real(kind=real_kind), dimension(nlev)       :: diffBi

    real(kind=real_kind) :: temp1

    integer :: np1     ! time level
    integer :: i, j, k ! equation index
    integer :: a, b    ! derivative index
    !----------------------------------------------------------------------------

    call t_startf('DerivPs_div_dpdn_vel')

    ! get time level
    np1 = fptr%tl%np1

    ! initialize output
    dps_div_dpdn_vel = 0.0d0
    
    ! compute dp/deta and d(dp/deta)/dps
    do k = 1,nlev
       dpdn(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
            - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))

       diffBi(k) = fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)
    end do

    do k = 1,nlev       
       do j = 1,np
          do i = 1,np

             temp1 = rrearth * fptr%base(ie)%rmetdet(i,j)
             
             ! a /= i, b = j (x-derivative)
             do a = 1,np
                if (a == i) cycle
                
                dps_div_dpdn_vel(a,j,i,j,k)  = temp1 * fptr%base(ie)%metdet(a,j) * diffBi(k) * fptr%deriv%Dvv(a,i) &
                     * (fptr%base(ie)%Dinv(a,j,1,1) * fptr%base(ie)%state%v(a,j,1,k,np1) &
                     + fptr%base(ie)%Dinv(a,j,1,2) * fptr%base(ie)%state%v(a,j,2,k,np1))
             end do
             
             ! a = i, b /= j (y-derivative)
             do b = 1,np
                if (b == j) cycle

                dps_div_dpdn_vel(i,b,i,j,k)  = temp1 * fptr%base(ie)%metdet(i,b) * diffBi(k) * fptr%deriv%Dvv(b,j) &
                     * (fptr%base(ie)%Dinv(i,b,2,1) * fptr%base(ie)%state%v(i,b,1,k,np1) &
                     + fptr%base(ie)%Dinv(i,b,2,2) * fptr%base(ie)%state%v(i,b,2,k,np1))
             end do

             ! a = i, b = j (x and y-derivative)             
             dps_div_dpdn_vel(i,j,i,j,k) = rrearth * diffBi(k) * fptr%deriv%Dvv(i,i) &
                  * (fptr%base(ie)%Dinv(i,j,1,1) * fptr%base(ie)%state%v(i,j,1,k,np1) &
                  + fptr%base(ie)%Dinv(i,j,1,2) * fptr%base(ie)%state%v(i,j,2,k,np1)) &
                  + rrearth * diffBi(k) * fptr%deriv%Dvv(j,j) &
                  * (fptr%base(ie)%Dinv(i,j,2,1) * fptr%base(ie)%state%v(i,j,1,k,np1) &
                  + fptr%base(ie)%Dinv(i,j,2,2) * fptr%base(ie)%state%v(i,j,2,k,np1))
          end do
       end do
    end do
   
    call t_stopf('DerivPs_div_dpdn_vel')
    
  end subroutine DerivPs_div_dpdn_vel



  subroutine DerivVel_etadot_dpdn(ie, fptr, du_etadot_dpdn, dv_etadot_dpdn)
    !----------------------------------------------------------------------------
    ! compute Du and Dv of (etadot * dpdn)
    !----------------------------------------------------------------------------
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
   
    implicit none

    !---------------------------------Arguments----------------------------------
    integer, intent(in) :: ie

    type(derived_type), pointer, intent(in) :: fptr 

    ! derivative matrix d(etadot_dpdn(i,j,k)) / du(a,b,c) => (a,b,c,i,j,k)
    real(kind=real_kind), intent(out) :: du_etadot_dpdn(np,np,nlev,np,np,nlev+1) 

    ! derivative matrix d(etadot_dpdn(i,j,k)) / dv(a,b,c) => (a,b,c,i,j,k)
    real(kind=real_kind), intent(out) :: dv_etadot_dpdn(np,np,nlev,np,np,nlev+1)
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    ! derivative matrix d(div_dpdn_vel(i,j,k)) / du(a,b,c) => (a,b,c,i,j,k)
    real(kind=real_kind) :: du_div_dpdn_vel(np,np,nlev,np,np,nlev) 

    ! derivative matrix d(div_dpdn_vel(i,j,k)) / dv(a,b,c) => (a,b,c,i,j,k)
    real(kind=real_kind) :: dv_div_dpdn_vel(np,np,nlev,np,np,nlev)

    real(kind=real_kind) :: Bi, Bim1

    integer :: i, j, k ! equation (node) index
    integer :: a, b, c ! derivative index
    !----------------------------------------------------------------------------

    call t_startf('DerivVel_etadot_dpdn')
    
    du_etadot_dpdn = 0.0d0
    dv_etadot_dpdn = 0.0d0
    
    call DerivVel_div_dpdn_vel(ie, fptr, du_div_dpdn_vel, dv_div_dpdn_vel)

    ! loop over equations (nodes)
    do k = 1,nlev-1
       
       Bim1 = fptr%hvcoord%hybi(k+1) - 1.0d0
       Bi   = fptr%hvcoord%hybi(k+1)
       
       do j = 1,np
          do i = 1,np
                 
             ! c <= k (contribution from both summations in etadot_dpdn)
             do c = 1,k

                ! a /= i, b = j (skips zeros in div_dpdn_vel)
                do a = 1,np
                   du_etadot_dpdn(a,j,c,i,j,k+1) = Bim1 * du_div_dpdn_vel(a,j,c,i,j,c)
                   dv_etadot_dpdn(a,j,c,i,j,k+1) = Bim1 * dv_div_dpdn_vel(a,j,c,i,j,c)
                end do

                ! a = i, b /= j (skips zeros in div_dpdn_vel)
                do b = 1,np
                   du_etadot_dpdn(i,b,c,i,j,k+1) = Bim1 * du_div_dpdn_vel(i,b,c,i,j,c)
                   dv_etadot_dpdn(i,b,c,i,j,k+1) = Bim1 * dv_div_dpdn_vel(i,b,c,i,j,c)
                end do
             end do

             ! c > k (contribution from first summation in etadot_dpdn)
             do c = k+1,nlev

                ! a /= i, b = j (skips zeros in div_dpdn_vel)
                do a = 1,np
                   du_etadot_dpdn(a,j,c,i,j,k+1) = Bi * du_div_dpdn_vel(a,j,c,i,j,c)
                   dv_etadot_dpdn(a,j,c,i,j,k+1) = Bi * dv_div_dpdn_vel(a,j,c,i,j,c)
                end do

                ! a = i, b /= j (skips zeros in div_dpdn_vel)
                do b = 1,np
                   du_etadot_dpdn(i,b,c,i,j,k+1) = Bi * du_div_dpdn_vel(i,b,c,i,j,c)
                   dv_etadot_dpdn(i,b,c,i,j,k+1) = Bi * dv_div_dpdn_vel(i,b,c,i,j,c)                   
                end do
             end do

          end do
       end do
    end do
   
    call t_stopf('DerivVel_etadot_dpdn')
    
  end subroutine DerivVel_etadot_dpdn



  subroutine DerivPs_etadot_dpdn(ie, fptr, dps_etadot_dpdn)
    !----------------------------------------------------------------------------
    ! compute Jacobian of (etadot * dpdn) REVISED VERSION
    !----------------------------------------------------------------------------
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
   
    implicit none

    !---------------------------------Arguments----------------------------------
    integer, intent(in) :: ie

    type(derived_type), pointer, intent(in) :: fptr 

    ! derivative matrix d(etadot_dpdn(i,j,k)) / dps(a,b) => (a,b,i,j,k)
    real(kind=real_kind), intent(out) :: dps_etadot_dpdn(np,np,np,np,nlev+1)
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    ! derivative matrix d(div_dpdn_vel(i,j,k)) / dps(a,b) => (a,b,i,j,k)
    real(kind=real_kind) :: dps_div_dpdn_vel(np,np,np,np,nlev)

    real(kind=real_kind) :: Bi

    integer :: i, j, k ! node index
    integer :: a, b, c ! derivative index
    !----------------------------------------------------------------------------

    call t_startf('DerivPs_etadot_dpdn')
    
    dps_etadot_dpdn = 0.0d0

    call DerivPs_div_dpdn_vel(ie, fptr, dps_div_dpdn_vel)

    do k = 1,nlev-1
       
       Bi = fptr%hvcoord%hybi(k+1)
       
       do j = 1,np
          do i = 1,np

             ! ==================================================================
             ! integral of div_dpdn_vel from 1 to nlev
             ! ==================================================================
             do c = 1,nlev

                ! a /= i, b = j
                do a = 1,np                 
                   if (a == i) cycle
                   
                   dps_etadot_dpdn(a,j,i,j,k+1) = dps_etadot_dpdn(a,j,i,j,k+1) &
                        + dps_div_dpdn_vel(a,j,i,j,c)
                end do

                ! a = i, b /= j
                do b = 1,np
                   if (b == j) cycle
                   
                   dps_etadot_dpdn(i,b,i,j,k+1) = dps_etadot_dpdn(i,b,i,j,k+1) &
                        + dps_div_dpdn_vel(i,b,i,j,c)
                end do

                ! a = i, b = j
                dps_etadot_dpdn(i,j,i,j,k+1) = dps_etadot_dpdn(i,j,i,j,k+1) &
                     + dps_div_dpdn_vel(i,j,i,j,c)
                
             end do

             ! ==================================================================  
             ! multiply by Bi
             ! ==================================================================

             ! a /= i, b = j
             do a = 1,np                 
                if (a == i) cycle
                dps_etadot_dpdn(a,j,i,j,k+1) = Bi * dps_etadot_dpdn(a,j,i,j,k+1)
             end do
             
             ! a = i, b /= j
             do b = 1,np
                if (b == j) cycle
                dps_etadot_dpdn(i,b,i,j,k+1) = Bi * dps_etadot_dpdn(i,b,i,j,k+1)
             end do            

             ! a = i, b = j
             dps_etadot_dpdn(i,j,i,j,k+1) = Bi * dps_etadot_dpdn(i,j,i,j,k+1)

             ! ==================================================================
             ! integral of div_dpdn_vel from 1 to k
             ! ==================================================================
             do c = 1,k

                ! a /= i, b = j 
                do a = 1,np                 
                   if (a == i) cycle
                   
                   dps_etadot_dpdn(a,j,i,j,k+1) = dps_etadot_dpdn(a,j,i,j,k+1) &
                        - dps_div_dpdn_vel(a,j,i,j,c)
                end do
             
                ! a = i, b /= j
                do b = 1,np
                   if (b == j) cycle

                   dps_etadot_dpdn(i,b,i,j,k+1) = dps_etadot_dpdn(i,b,i,j,k+1) &
                        - dps_div_dpdn_vel(i,b,i,j,c)
                end do

                ! a = i, b = j
                dps_etadot_dpdn(i,j,i,j,k+1) = dps_etadot_dpdn(i,j,i,j,k+1) &
                     - dps_div_dpdn_vel(i,j,i,j,c)

             end do
             
          end do
       end do
    end do
   
    call t_stopf('DerivPs_etadot_dpdn')
    
  end subroutine DerivPs_etadot_dpdn



  ! =============================================================================
  ! The following subroutines are for zeroing out blocks of the Jacobian matrix
  ! =============================================================================



  subroutine Zero_Ablock(ie, nnz_max, row_nnz, Urows, Ucols, Uvals, &
       Vrows, Vcols, Vvals, c_ptr_to_object) bind(C,name="Zero_Ablock")
    !----------------------------------------------------------------------------
    ! Zeros out nonzero entries of the A block of the Jacobian
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type
    
    implicit none
    
    !---------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per element
    
    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 
    
    ! row indices 
    integer(c_int), dimension(np*np*nlev), intent(out) :: Urows
    integer(c_int), dimension(np*np*nlev), intent(out) :: Vrows
    
    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Ucols
    integer(c_int), dimension(nnz_max), intent(out) :: Vcols
    
    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Uvals   
    real(c_double), dimension(nnz_max), intent(out) :: Vvals
    
    type(c_ptr) :: c_ptr_to_object
    !----------------------------------------------------------------------------
    
    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr
    
    integer :: a, b, c     ! derivative index
    integer :: i, j, k     ! equation index
    integer :: ridx, cidx  ! global row/column index 
    !----------------------------------------------------------------------------

#ifdef DEBUG_PRINT_ON   
    write(*,*) "Start: Zero_Ablock",iam
#endif
    
    call t_startf('Zero_Ablock')
    
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
    
    ! initialize arrays
    Urows = 0
    Ucols = 0
    Uvals = 0.0d0
    
    Vrows = 0
    Vcols = 0
    Vvals = 0.0d0
    
    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0
    
    ! number of nonzeros per row
    row_nnz = 2*(np+np-1)*nlev
    
    do k = 1,nlev
       do j = 1,np
          do i = 1,np
             
             ! set row index for U and V equations
             ridx = ridx + 1
             Urows(ridx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
             Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)
             
             do c = 1,nlev 
                
                ! a /= i, b = j, c = *
                ! ------------------------------------------------------------
                do a = 1,np
                   if (a == i) cycle
                   
                   cidx = cidx + 1
                   Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)
                   
                   cidx = cidx + 1
                   Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)
                   
                end do
                
                ! a = i, b /= j, c = *
                ! ------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle
                   
                   cidx = cidx + 1
                   Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)
                   
                   cidx = cidx + 1
                   Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)
                   
                end do
                
                ! a = i, b = j, c = *
                ! ------------------------------------------------------------
                cidx = cidx + 1
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)
                
                cidx = cidx + 1
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)
                
             end do
             
          end do
       end do
    end do
        
    call t_stopf('Zero_Ablock')

#ifdef DEBUG_PRINT_ON   
    write(*,*) "End: Zero_Ablock",iam
#endif

  end subroutine Zero_Ablock



  subroutine Zero_Bblock(ie, nnz_max, row_nnz, Urows, Ucols, Uvals, &
       Vrows, Vcols, Vvals, c_ptr_to_object) bind(C,name="Zero_Bblock")
    !----------------------------------------------------------------------------
    ! Zeros out nonzero entries of the B block of the Jacobian
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type

    implicit none

    !---------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per process

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    ! row indices 
    integer(c_int), dimension(np*np*nlev), intent(out) :: Urows
    integer(c_int), dimension(np*np*nlev), intent(out) :: Vrows

    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Ucols
    integer(c_int), dimension(nnz_max), intent(out) :: Vcols

    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Uvals   
    real(c_double), dimension(nnz_max), intent(out) :: Vvals
    
    type(c_ptr) :: c_ptr_to_object
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr

    integer :: a, b        ! derivative index
    integer :: i, j, k     ! equation index
    integer :: ridx, cidx  ! global row/column index 
    !----------------------------------------------------------------------------
    
    call t_startf('Zero_Bblock')
    
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
        
    ! initialize arrays
    Urows = 0
    Ucols = 0
    Uvals = 0.0d0

    Vrows = 0
    Vcols = 0
    Vvals = 0.0d0
    
    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = np+np-1
    
    do k = 1,nlev
       do j = 1,np
          do i = 1,np

             ! set row index for U and V equations
             ridx = ridx + 1
             Urows(ridx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
             Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)

             ! a /= i, b = j
             ! ---------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Psoffset
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Psoffset

             end do

             ! a = i, b /= j
             ! ---------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Psoffset
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Psoffset

             end do

             ! a = i, b = j
             ! ---------------------------------------------------------------
             cidx = cidx + 1
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Psoffset
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Psoffset

          end do
       end do
    end do
    
    call t_stopf('Zero_Bblock')
    
  end subroutine Zero_Bblock



  subroutine Zero_Cblock(ie, nnz_max, row_nnz, Urows, Ucols, Uvals, &
       Vrows, Vcols, Vvals, c_ptr_to_object) bind(C,name="Zero_Cblock")
    !----------------------------------------------------------------------------
    ! Zeros out nonzero entries of the C block of the Jacobian
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type

    implicit none

    !---------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per process

    ! number of nonzero values per row
    integer(c_int), dimension(np*np*nlev), intent(out) :: row_nnz 

    ! row indices 
    integer(c_int), dimension(np*np*nlev), intent(out) :: Urows
    integer(c_int), dimension(np*np*nlev), intent(out) :: Vrows

    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Ucols
    integer(c_int), dimension(nnz_max), intent(out) :: Vcols

    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Uvals   
    real(c_double), dimension(nnz_max), intent(out) :: Vvals
    
    type(c_ptr) :: c_ptr_to_object
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr

    integer :: a, b, c     ! derivative index
    integer :: i, j, k     ! equation index
    integer :: ridx, cidx  ! global row/column index 
    !----------------------------------------------------------------------------

    call t_startf('Zero_Cblock')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
    
    ! initialize arrays
    Urows = 0
    Ucols = 0
    Uvals = 0.0d0

    Vrows = 0
    Vcols = 0
    Vvals = 0.0d0
    
    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! initialize number of nonzeros per row (varies)
    row_nnz = 0
       
    do k = 1,nlev
       do j = 1,np
          do i = 1,np

             ! set row index for U and V equations
             ridx = ridx + 1
             Urows(ridx) = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
             Vrows(ridx) = fptr%base(ie)%gdofP(i,j) + Voffset(k)

             ! a /= i, b = j, c = k
             ! ---------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1
                Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Toffset(k)
                Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Toffset(k)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1

             end do

             ! a = i, b /= j, c = k
             ! ---------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1
                Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Toffset(k)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Toffset(k)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1

             end do

             ! a = i, b = j, c = k
             ! ---------------------------------------------------------------
             cidx = cidx + 1
             Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)
             Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

             ! update nnz per row
             row_nnz(ridx) = row_nnz(ridx) + 1

             do c = k+1,nlev 

                ! a /= i, b = j, c = *
                ! ------------------------------------------------------------
                do a = 1,np
                   if (a == i) cycle

                   cidx = cidx + 1
                   Ucols(cidx) = fptr%base(ie)%gdofP(a,j) + Toffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(a,j) + Toffset(c)

                   ! update nnz per row
                   row_nnz(ridx) = row_nnz(ridx) + 1

                end do

                ! a = i, b /= j, c = *
                ! ------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle

                   cidx = cidx + 1
                   Ucols(cidx) = fptr%base(ie)%gdofP(i,b) + Toffset(c)
                   Vcols(cidx) = fptr%base(ie)%gdofP(i,b) + Toffset(c)

                   ! update nnz per row
                   row_nnz(ridx) = row_nnz(ridx) + 1

                end do

                ! a = i, b = j, c = *
                ! ------------------------------------------------------------
                cidx = cidx + 1
                Ucols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(c)
                Vcols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(c)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1

             end do

          end do
       end do
    end do

    call t_stopf('Zero_Cblock')

  end subroutine Zero_Cblock



  subroutine Zero_Dblock(ie, nnz_max, row_nnz, Rows, Cols, Vals, &
       c_ptr_to_object) bind(C,name="Zero_Dblock")
    !----------------------------------------------------------------------------
    ! Zeros out nonzero entries of the D block of the Jacobian
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type

    implicit none

    !---------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per process

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    ! row indices
    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows
    
    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Cols
    
    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Vals
    
    type(c_ptr) :: c_ptr_to_object
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr

    integer :: a, b, c     ! derivative index
    integer :: i, j        ! equation index
    integer :: ridx, cidx  ! global row/column index 
    !----------------------------------------------------------------------------

    call t_startf('Zero_Dblock')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0
    
    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = 2*(np+np-1)*nlev

    do j = 1,np
       do i = 1,np

          ! set row index for Ps equation
          ridx = ridx + 1
          Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Psoffset

          do c = 1,nlev 

             ! a /= i, b = j, c = *
             ! ---------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)

                cidx = cidx + 1
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)

             end do

             ! a = i, b /= j, c = *
             ! ---------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)

                cidx = cidx + 1
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)

             end do

             ! a = i, b = j, c = *
             ! ---------------------------------------------------------------
             cidx = cidx + 1
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)

             cidx = cidx + 1
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)

          end do

       end do
    end do
    
    call t_stopf('Zero_Dblock')

  end subroutine Zero_Dblock



  subroutine Zero_Eblock(ie, nnz_max, row_nnz, Rows, Cols, Vals, &
       c_ptr_to_object) bind(C,name="Zero_Eblock")
    !----------------------------------------------------------------------------
    ! Zeros out nonzero entries of the E block of the Jacobian
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type

    implicit none

    !---------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per process

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 

    ! row indices
    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows
    
    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Cols

    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Vals
    
    type(c_ptr) :: c_ptr_to_object
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr

    integer :: a, b        ! derivative index
    integer :: i, j        ! equation index
    integer :: ridx, cidx  ! global row/column index 
    !----------------------------------------------------------------------------

    call t_startf('Zero_Eblock')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0
    
    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = np+np-1

    do j = 1,np
       do i = 1,np

          ! set row index for Ps equation
          ridx = ridx + 1
          Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Psoffset

          ! a /= i, b = j
          ! ---------------------------------------------------------------
          do a = 1,np
             if (a == i) cycle

             cidx = cidx + 1
             Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Psoffset

          end do

          ! a = i, b /= j
          ! ---------------------------------------------------------------
          do b = 1,np
             if (b == j) cycle

             cidx = cidx + 1
             Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Psoffset

          end do

          ! a = i, b = j
          ! ---------------------------------------------------------------
          cidx = cidx + 1
          Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Psoffset

       end do
    end do
    
    call t_stopf('Zero_Eblock')

  end subroutine Zero_Eblock
  
  
  
  subroutine Zero_Fblock(ie, nnz_max, row_nnz, Rows, Cols, Vals, &
       c_ptr_to_object) bind(C,name="Zero_Fblock")
    !----------------------------------------------------------------------------
    ! Zeros out nonzero entries of the F block of the Jacobian
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type

    implicit none

    !---------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per process

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz 
    
    ! row indices
    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows 

    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Cols
    
    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Vals
    
    type(c_ptr) :: c_ptr_to_object
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr

    integer :: a, b, c    ! derivative index
    integer :: i, j, k    ! equation index
    integer :: ridx, cidx ! global row/column index 
    !----------------------------------------------------------------------------

    call t_startf('Zero_Fblock')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0
    
    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = 2*(np+np-1)*nlev

    do k = 1,nlev
       do j = 1,np
          do i = 1,np

             ! set row index for T equations
             ridx = ridx + 1
             Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

             do c = 1,nlev 

                ! a /= i, b = j, c = *
                ! ------------------------------------------------------------
                do a = 1,np
                   if (a == i) cycle

                   cidx = cidx + 1
                   Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Uoffset(c)

                   cidx = cidx + 1
                   Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Voffset(c)

                end do

                ! a = i, b /= j, c = *
                ! ------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle

                   cidx = cidx + 1
                   Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Uoffset(c)

                   cidx = cidx + 1
                   Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Voffset(c)

                end do

                ! a = i, b = j, c = *
                ! ------------------------------------------------------------
                cidx = cidx + 1
                Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Uoffset(c)

                cidx = cidx + 1
                Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Voffset(c)
             end do

          end do
       end do
    end do
    
    call t_stopf('Zero_Fblock')

  end subroutine Zero_Fblock



  subroutine Zero_Gblock(ie, nnz_max, row_nnz, Rows, Cols, Vals, &
       c_ptr_to_object) bind(C,name="Zero_Gblock")
    !----------------------------------------------------------------------------
    ! Zeros out nonzero entries of the G block of the Jacobian
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type

    implicit none

    !---------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per process

    ! number of nonzero values per row
    integer(c_int), intent(out) :: row_nnz

    ! row indices
    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows 
    
    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Cols

    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Vals
    
    type(c_ptr) :: c_ptr_to_object
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr

    integer :: a, b        ! derivative index
    integer :: i, j, k     ! equation index
    integer :: ridx, cidx  ! global row/column index 
    !----------------------------------------------------------------------------

    call t_startf('Zero_Gblock')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0
    
    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! number of nonzeros per row
    row_nnz = np+np-1
    
    do k = 1,nlev
       do j = 1,np
          do i = 1,np

             ! set row index for T equations
             ridx = ridx + 1
             Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

             ! a /= i, b = j, c = *
             ! ---------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + PsOffset
             end do

             ! a = i, b /= j, c = *
             ! ---------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + PsOffset
             end do

             ! a = i, b = j, c = *
             ! ---------------------------------------------------------------
             cidx = cidx + 1
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + PsOffset

          end do
       end do
    end do
    
    call t_stopf('Zero_Gblock')

  end subroutine Zero_Gblock



  subroutine Zero_Hblock(ie, nnz_max, row_nnz, Rows, Cols, Vals, &
       c_ptr_to_object) bind(C,name="Zero_Hblock")
    !----------------------------------------------------------------------------
    ! Zeros out nonzero entries of the H block of the Jacobian
    !----------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev
    use prim_derived_type_mod, only : derived_type

    implicit none

    !---------------------------------Arguments----------------------------------
    integer(c_int), intent(in), value :: ie      ! element index
    integer(c_int), intent(in), value :: nnz_max ! max nonzeros per process

    ! number of nonzero values per row
    integer(c_int), dimension(np*np*nlev), intent(out) :: row_nnz 

    ! row indices
    integer(c_int), dimension(np*np*nlev), intent(out) :: Rows

    ! column indices
    integer(c_int), dimension(nnz_max), intent(out) :: Cols

    ! matrix entries
    real(c_double), dimension(nnz_max), intent(out) :: Vals
    
    type(c_ptr) :: c_ptr_to_object
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr

    integer :: a, b        ! derivative index
    integer :: i, j, k     ! equation index
    integer :: ridx, cidx  ! global row/column index 
    !----------------------------------------------------------------------------

    call t_startf('Zero_Hblock')

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    ! initialize arrays
    Rows = 0
    Cols = 0
    Vals = 0.0d0
    
    ! initialize row/column storage indicies 
    ridx = 0
    cidx = 0

    ! initialize number of nonzeros per row (varies)
    row_nnz = 0

    do k = 1,nlev
       do j = 1,np
          do i = 1,np

             ! set row index for T equations
             ridx = ridx + 1
             Rows(ridx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

             ! a /= i, b = j, c = k
             ! ---------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle

                cidx = cidx + 1
                Cols(cidx) = fptr%base(ie)%gdofP(a,j) + Toffset(k)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1
             end do

             ! a = i, b /= j, c = k
             ! ---------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                cidx = cidx + 1
                Cols(cidx) = fptr%base(ie)%gdofP(i,b) + Toffset(k)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1
             end do

             ! a = i, b = j, c = k
             ! ---------------------------------------------------------------
             cidx = cidx + 1
             Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k)

             ! update nnz per row
             row_nnz(ridx) = row_nnz(ridx) + 1

             if (k > 1) then
                ! a = i, b = j, c = k-1
                ! ------------------------------------------------------------
                cidx = cidx + 1
                Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k-1)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1
             end if

             if (k < nlev) then
                ! a = i, b = j, c = k+1
                ! ------------------------------------------------------------
                cidx = cidx + 1
                Cols(cidx) = fptr%base(ie)%gdofP(i,j) + Toffset(k+1)

                ! update nnz per row
                row_nnz(ridx) = row_nnz(ridx) + 1
             end if

          end do
       end do
    end do
    
    call t_stopf('Zero_Hblock')

  end subroutine Zero_Hblock

  
end module prim_jacobian_sparse_mod

#endif 
! TRILINOS
