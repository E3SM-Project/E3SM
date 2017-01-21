#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
  
#ifdef TRILINOS

! ==============================================================================
! This module contains subroutines for computing the blocks of the analytic 
! Jacobian for each element as well as subroutines for computing matrix-vector 
! products with the Jacbian blocks for preconditioning.
! ==============================================================================

!#define DEBUG_PRINT_ON   

! output Jacobian enteries using local index index rather than global index
!#define LOCAL_IDX
  
module prim_jacobian_dense_mod
  
  use kinds, only        : real_kind
  use edgetype_mod, only : EdgeBuffer_t
  use parallel_mod, only : abortmp, iam
  use perf_mod, only     : t_startf, t_stopf
 
  implicit none
  private
  save

  public :: prim_jacobian_dense_init, prim_jacobian_dense_finalize
  public :: compute_jacobian

  ! make subroutines for building Jacobian blocks public for testing
  public :: print_analytic_jacobian
  public :: zero_jacobian
  public :: Jac_div_dpdn_vel, Jac_etadot_dpdn
  public :: Jac_TimeDeriv, Jac_Vorticity
  public :: Jac_Coriolis, Jac_Kinetic_Energy
  public :: Jac_Geopotential, Jac_press_grad
  public :: Jac_VertAdv, Jac_HorizAdv_Temp
  public :: Jac_pressure_vert_vel
  public :: Jac_int_div_dpdn_vel

  ! exchange buffer
  type(EdgeBuffer_t) :: edge3p1

  ! =============================================================================
  ! For each element (ie) the local Jacobian matrix blocks are stored as 
  ! dX(i,j,k)/Y(a,b,c) = dX_dY(a,b,c,i,j,k,ie)
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

  real(kind=real_kind), allocatable :: dU_dU(:,:,:,:,:,:,:) ! A = J(1,1) block
  real(kind=real_kind), allocatable :: dU_dV(:,:,:,:,:,:,:) ! A = J(1,1) block
  real(kind=real_kind), allocatable :: dV_dU(:,:,:,:,:,:,:) ! A = J(1,1) block
  real(kind=real_kind), allocatable :: dV_dV(:,:,:,:,:,:,:) ! A = J(1,1) block

  real(kind=real_kind), allocatable :: dU_dPs(:,:,:,:,:,:)  ! B = J(1,2) block
  real(kind=real_kind), allocatable :: dV_dPs(:,:,:,:,:,:)  ! B = J(1,2) block

  real(kind=real_kind), allocatable :: dU_dT(:,:,:,:,:,:,:) ! C = J(1,3) block
  real(kind=real_kind), allocatable :: dV_dT(:,:,:,:,:,:,:) ! C = J(1,3) block

  real(kind=real_kind), allocatable :: dPs_dU(:,:,:,:,:,:)  ! D = J(2,1) block
  real(kind=real_kind), allocatable :: dPs_dV(:,:,:,:,:,:)  ! D = J(2,1) block

  real(kind=real_kind), allocatable :: dPs_dPs(:,:,:,:,:)   ! E = J(2,2) block

  real(kind=real_kind), allocatable :: dT_dU(:,:,:,:,:,:,:) ! F = J(3,1) block
  real(kind=real_kind), allocatable :: dT_dV(:,:,:,:,:,:,:) ! F = J(3,1) block

  real(kind=real_kind), allocatable :: dT_dPs(:,:,:,:,:,:)  ! G = J(3,1) block

  real(kind=real_kind), allocatable :: dT_dT(:,:,:,:,:,:,:) ! H = J(3,3) block

contains

  subroutine prim_jacobian_dense_init(par, elem)
    ! ---------------------------------------------------------------------------
    ! Allocate Jacobian blocks data
    ! ---------------------------------------------------------------------------
    use parallel_mod, only   : parallel_t
    use element_mod, only    : element_t
    use edge_mod, only       : initEdgeBuffer
    use dimensions_mod, only : np, nlev, nelemd
    use control_mod, only    : rsplit

    implicit none

    type(parallel_t), intent(in) :: par
    type(element_t),  intent(in), target :: elem(:)

    ! allocate exchange buffe
    if (rsplit==0) then
       call initEdgeBuffer(par, edge3p1, elem, 3*nlev+1)
    else
       ! need extra buffer space for dp3d
       call initEdgeBuffer(par, edge3p1, elem, 4*nlev+1)
    endif

    ! allocate Jacobian blocks
    ! A = J(1,1) blocks 
    allocate(dU_dU(np,np,nlev,np,np,nlev,nelemd))
    allocate(dU_dV(np,np,nlev,np,np,nlev,nelemd))

    allocate(dV_dU(np,np,nlev,np,np,nlev,nelemd))
    allocate(dV_dV(np,np,nlev,np,np,nlev,nelemd))

    ! B = J(1,2) block
    allocate(dU_dPs(np,np,np,np,nlev,nelemd))
    allocate(dV_dPs(np,np,np,np,nlev,nelemd))

    ! C = J(1,3) block
    allocate(dU_dT(np,np,nlev,np,np,nlev,nelemd))
    allocate(dV_dT(np,np,nlev,np,np,nlev,nelemd))

    ! D = J(2,1) block
    allocate(dPs_dU(np,np,nlev,np,np,nelemd))
    allocate(dPs_dV(np,np,nlev,np,np,nelemd))

    ! E = J(2,2) block
    allocate(dPs_dPs(np,np,np,np,nelemd))

    ! F = J(3,1) block
    allocate(dT_dU(np,np,nlev,np,np,nlev,nelemd))
    allocate(dT_dV(np,np,nlev,np,np,nlev,nelemd))

    ! G = J(3,1) block
    allocate(dT_dPs(np,np,np,np,nlev,nelemd))

    ! H = J(3,3) block
    allocate(dT_dT(np,np,nlev,np,np,nlev,nelemd))

  end subroutine prim_jacobian_dense_init



  subroutine prim_jacobian_dense_finalize()
    ! ---------------------------------------------------------------------------
    ! Free exchange buffer and Jacobian blocks
    ! ---------------------------------------------------------------------------
    use edge_mod, only : FreeEdgeBuffer

    implicit none

    call FreeEdgeBuffer(edge3p1)

    deallocate(dU_dU, dU_dV, dV_dU, dV_dV)
    deallocate(dU_dPs, dV_dPs)
    deallocate(dU_dT, dV_dT)
    deallocate(dPs_dU, dPs_dV)
    deallocate(dPs_dPs)
    deallocate(dT_dU, dT_dV)
    deallocate(dT_dPs)
    deallocate(dT_dT)

  end subroutine prim_jacobian_dense_finalize
  
  
  subroutine zero_jacobian()
    
    ! set jacobian blcks to zero
    dU_dU   = 0.0d0
    dU_dV   = 0.0d0    
    dV_dU   = 0.0d0 
    dV_dV   = 0.0d0
    dU_dPs  = 0.0d0  
    dV_dPs  = 0.0d0 
    dU_dT   = 0.0d0
    dV_dT   = 0.0d0
    dPs_dU  = 0.0d0
    dPs_dV  = 0.0d0
    dPs_dPs = 0.0d0
    dT_dU   = 0.0d0
    dT_dV   = 0.0d0
    dT_dPs  = 0.0d0
    dT_dT   = 0.0d0

  end subroutine zero_jacobian

  
  subroutine compute_jacobian(fptr)
    ! ---------------------------------------------------------------------------
    ! Fill Jacobian matrix blocks
    ! ---------------------------------------------------------------------------

    use prim_derived_type_mod, only : derived_type

    implicit none

    ! --------------------------------Arguments----------------------------------
    type(derived_type), pointer, intent(in) :: fptr
    ! ---------------------------------------------------------------------------

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*) ">>> HOMME: compute_jacobian"
#endif

    call t_startf('Compute_Jac')

    call zero_jacobian()

    ! time derivatives in velocity, temperature, and surface pressure equ
    call Jac_TimeDeriv(fptr)

    ! vorticity in velocity equ
    call Jac_Vorticity(fptr)

    ! Coriolis force in velocity equ
    call Jac_Coriolis(fptr)

    ! gradient of kinetic energy in velocity equ
    call Jac_Kinetic_Energy(fptr)

    ! gradient of geopotential in velocity equ
    call Jac_Geopotential(fptr)

    ! pressure gradient in velocity equ
    call Jac_press_grad(fptr) 

    ! vertical advection in velocity and temperature equ
    call Jac_VertAdv(fptr)

    ! horizontal advection in temperature equ
    call Jac_HorizAdv_Temp(fptr)

    ! pressure vertical velocity in temperature equ
    call Jac_pressure_vert_vel(fptr)

    ! integral of divergence(dp/dn * velocity) in surface pressure equ
    call Jac_int_div_dpdn_vel(fptr)

    call t_stopf('Compute_Jac')

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*) "<<< HOMME: compute_jacobian"
#endif

  end subroutine compute_jacobian



  subroutine print_analytic_jacobian(fptr, Jname)
    ! ---------------------------------------------------------------------------
    ! Outputs analytic Jacobian to file
    ! ---------------------------------------------------------------------------
    use dimensions_mod, only           : np, nlev, nvar, nelemd
    use prim_derived_type_mod, only    : derived_type
    use prim_jacobian_sparse_mod, only : Uoffset, Voffset, Toffset, Psoffset
    
    implicit none
    
    ! --------------------------------Arguments----------------------------------
    type(derived_type), pointer, intent(in) :: fptr

    character(len=1024), intent(in) :: Jname ! name of Jacobian subroutine
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    integer :: i, j, k, ie ! equation loop
    integer :: a, b, c     ! derivative loop
    integer :: ridx, cidx  ! row and column index
    integer :: fid         ! file id for output
    
    character(len=1024) :: fname ! output file name
    ! ---------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------
    ! A Block 
    ! ---------------------------------------------------------------------------

    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'Ablock_',trim(Jname),'_proc',iam,'.txt'
    open(fid, file=fname, form='formatted')

    write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 

    ! dU_dU
    do ie = 1,nelemd
       
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                
#ifdef LOCAL_IDX
                ! local row index
                ridx = i + (j-1)*np + (k-1)*np*np
#else
                ! global row index
                ridx = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
#endif
                
                do c = 1,nlev
                   do b = 1,np
                      do a = 1,np
                         
#ifdef LOCAL_IDX
                         ! local column index
                         cidx = a + (b-1)*np + (c-1)*np*np
#else
                         ! global column index
                         cidx = fptr%base(ie)%gdofP(a,b) + Uoffset(c)
#endif
                         write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dU_dU(a,b,c,i,j,k,ie) 

                      end do
                   end do
                end do

             end do
          end do
       end do

    end do

    ! dU_dV
    do ie = 1,nelemd
       
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = i + (j-1)*np + (k-1)*np*np
#else
                ! global row index
                ridx = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
#endif
                
                do c = 1,nlev
                   do b = 1,np
                      do a = 1,np
                         
#ifdef LOCAL_IDX
                         ! local column index
                         cidx = np*np*nlev + a + (b-1)*np + (c-1)*np*np
#else
                         ! global column index                        
                         cidx = fptr%base(ie)%gdofP(a,b) + Voffset(c)
#endif
                         write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dU_dV(a,b,c,i,j,k,ie) 
                         
                      end do
                   end do
                end do

             end do
          end do
       end do

    end do

    ! dV_dU
    do ie = 1,nelemd

       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = np*np*nlev + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index
                ridx = fptr%base(ie)%gdofP(i,j) + Voffset(k)
#endif
                
                do c = 1,nlev
                   do b = 1,np
                      do a = 1,np

#ifdef LOCAL_IDX
                         ! local column index
                         cidx = a + (b-1)*np + (c-1)*np*np
#else
                         ! global column index
                         cidx = fptr%base(ie)%gdofP(a,b) + Uoffset(c)
#endif
                         write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dV_dU(a,b,c,i,j,k,ie) 
                         
                      end do
                   end do
                end do

             end do
          end do
       end do

    end do

    ! dV_dV
    do ie = 1,nelemd

       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = np*np*nlev + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index
                ridx = fptr%base(ie)%gdofP(i,j) + Voffset(k)
#endif
                
                do c = 1,nlev
                   do b = 1,np
                      do a = 1,np

#ifdef LOCAL_IDX
                         ! local column index
                         cidx = np*np*nlev + a + (b-1)*np + (c-1)*np*np
#else
                         ! global column index
                         cidx = fptr%base(ie)%gdofP(a,b) + Voffset(c)
#endif
                         write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dV_dV(a,b,c,i,j,k,ie) 
                         
                      end do
                   end do
                end do

             end do
          end do
       end do

    end do

    close(fid)

    ! ---------------------------------------------------------------------------
    ! B Block 
    ! ---------------------------------------------------------------------------

    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'Bblock_',trim(Jname),'_proc',iam,'.txt'
    open(fid, file=fname, form='formatted')

    write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 

    ! dU_dPs 
    do ie = 1,nelemd

       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = i + (j-1)*np + (k-1)*np*np
#else
                ! global row index
                ridx = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
#endif 
                
                do b = 1,np
                   do a = 1,np
#ifdef LOCAL_IDX
                      ! local column index
                      cidx = 2*np*np*nlev + a + (b-1)*np
#else                     
                      ! global column index
                      cidx = fptr%base(ie)%gdofP(a,b) + Psoffset
#endif
                      write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dU_dPs(a,b,i,j,k,ie) 
                      
                   end do
                end do

             end do
          end do
       end do
       
    end do

    ! dV_dPs
    do ie = 1,nelemd

       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = np*np*nlev + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index
                ridx = fptr%base(ie)%gdofP(i,j) + Voffset(k)
#endif
                
                do b = 1,np
                   do a = 1,np
                      
#ifdef LOCAL_IDX
                      ! local column index
                      cidx = 2*np*np*nlev + a + (b-1)*np
#else 
                      ! global column index
                      cidx = fptr%base(ie)%gdofP(a,b) + Psoffset
#endif
                      write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dV_dPs(a,b,i,j,k,ie) 
                      
                   end do
                end do
                
             end do
          end do
       end do
       
    end do
    
    close(fid)

    ! ---------------------------------------------------------------------------
    ! C Block 
    ! ---------------------------------------------------------------------------
    
    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'Cblock_',trim(Jname),'_proc',iam,'.txt'
    open(fid, file=fname, form='formatted')

    write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 

    ! dU_dT
    do ie = 1,nelemd

       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = i + (j-1)*np + (k-1)*np*np
#else
                ! global row index
                ridx = fptr%base(ie)%gdofP(i,j) + Uoffset(k)
#endif
                
                do c = 1,nlev
                   do b = 1,np
                      do a = 1,np

#ifdef LOCAL_IDX
                         ! local column index
                         cidx = 2*np*np*nlev + np*np + a + (b-1)*np + (c-1)*np*np
#else
                         ! global column index
                         cidx = fptr%base(ie)%gdofP(a,b) + Toffset(c)
#endif
                         write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dU_dT(a,b,c,i,j,k,ie) 
                         
                      end do
                   end do
                end do

             end do
          end do
       end do

    end do

    ! dV_dT
    do ie = 1,nelemd

       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = np*np*nlev + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index
                ridx = fptr%base(ie)%gdofP(i,j) + Voffset(k)
#endif
                
                do c = 1,nlev
                   do b = 1,np
                      do a = 1,np

#ifdef LOCAL_IDX
                         ! local column index
                         cidx = 2*np*np*nlev + np*np + a + (b-1)*np + (c-1)*np*np
#else
                         ! global column index
                         cidx = fptr%base(ie)%gdofP(a,b) + Toffset(c)
#endif
                         write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dV_dT(a,b,c,i,j,k,ie) 
                         
                      end do
                   end do
                end do

             end do
          end do
       end do

    end do

    close(fid)

    ! ---------------------------------------------------------------------------
    ! D Block 
    ! ---------------------------------------------------------------------------

    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'Dblock_',trim(Jname),'_proc',iam,'.txt'
    open(fid, file=fname, form='formatted')

    write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 

    ! dPs_dU
    do ie = 1,nelemd

       do j = 1,np
          do i = 1,np

#ifdef LOCAL_IDX
             ! local row index
             ridx = 2*np*np*nlev + i + (j-1)*np
#else
             ! global row index            
             ridx = fptr%base(ie)%gdofP(i,j) + Psoffset
#endif
             
             do c = 1,nlev
                do b = 1,np
                   do a = 1,np
      
#ifdef LOCAL_IDX
                      ! local column index
                      cidx = a + (b-1)*np + (c-1)*np*np
#else
                      ! global column index                      
                      cidx = fptr%base(ie)%gdofP(a,b) + Uoffset(c)
#endif 
                      write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dPs_dU(a,b,c,i,j,ie) 
                      
                   end do
                end do
             end do
             
          end do
       end do

    end do

    ! dPs_dV
    do ie = 1,nelemd

       do j = 1,np
          do i = 1,np

#ifdef LOCAL_IDX
             ! local row index
             ridx = 2*np*np*nlev + i + (j-1)*np
#else
             ! global row index            
             ridx = fptr%base(ie)%gdofP(i,j) + Psoffset
#endif
             
             do c = 1,nlev
                do b = 1,np
                   do a = 1,np

#ifdef LOCAL_IDX
                      ! local column index
                      cidx = np*np*nlev + a + (b-1)*np + (c-1)*np*np
#else
                      ! global column index                                           
                      cidx = fptr%base(ie)%gdofP(a,b) + Voffset(c)
#endif
                      write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dPs_dV(a,b,c,i,j,ie) 
                      
                   end do
                end do
             end do
             
          end do
       end do

    end do
    
    close(fid)

    ! ---------------------------------------------------------------------------
    ! E Block 
    ! ---------------------------------------------------------------------------

    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'Eblock_',trim(Jname),'_proc',iam,'.txt'
    open(fid, file=fname, form='formatted')

    write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 

    ! dPs_dPs
    do ie = 1,nelemd

       do j = 1,np
          do i = 1,np

#ifdef LOCAL_IDX
             ! local row index
             ridx = 2*np*np*nlev + i + (j-1)*np
#else
             ! global row index            
             ridx = fptr%base(ie)%gdofP(i,j) + Psoffset
#endif
             
             do b = 1,np
                do a = 1,np

#ifdef LOCAL_IDX
                   ! local column index
                   cidx = 2*np*np*nlev + a + (b-1)*np
#else
                   ! global column index
                   cidx = fptr%base(ie)%gdofP(a,b) + Psoffset
#endif
                   write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dPs_dPs(a,b,i,j,ie) 
                   
                end do
             end do
             
          end do
       end do

    end do

    close(fid)

    ! ---------------------------------------------------------------------------
    ! F Block 
    ! ---------------------------------------------------------------------------

    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'Fblock_',trim(Jname),'_proc',iam,'.txt'
    open(fid, file=fname, form='formatted')

    write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 

    ! dT_dU
    do ie = 1,nelemd

       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = 2*np*np*nlev + np*np + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index            
                ridx = fptr%base(ie)%gdofP(i,j) + Toffset(k)
#endif
                
                do c = 1,nlev
                   do b = 1,np
                      do a = 1,np

#ifdef LOCAL_IDX
                         ! local column index
                         cidx = a + (b-1)*np + (c-1)*np*np
#else
                         ! global column index            
                         cidx = fptr%base(ie)%gdofP(a,b) + Uoffset(c)
#endif
                         write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dT_dU(a,b,c,i,j,k,ie) 
                         
                      end do
                   end do
                end do

             end do
          end do
       end do

    end do

    ! dT_dV
    do ie = 1,nelemd

       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = 2*np*np*nlev + np*np + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index            
                ridx = fptr%base(ie)%gdofP(i,j) + Toffset(k)
#endif
                
                do c = 1,nlev
                   do b = 1,np
                      do a = 1,np

#ifdef LOCAL_IDX
                         ! local column index
                         cidx = np*np*nlev + a + (b-1)*np + (c-1)*np*np
#else
                         ! global column index            
                         cidx = fptr%base(ie)%gdofP(a,b) + Voffset(c)
#endif
                         write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dT_dV(a,b,c,i,j,k,ie) 
                         
                      end do
                   end do
                end do

             end do
          end do
       end do

    end do

    close(fid)

    ! ---------------------------------------------------------------------------
    ! G Block 
    ! ---------------------------------------------------------------------------

    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'Gblock_',trim(Jname),'_proc',iam,'.txt'
    open(fid, file=fname, form='formatted')

    write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 
    
    ! dT_dPs
    do ie = 1,nelemd

       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = 2*np*np*nlev + np*np + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index            
                ridx = fptr%base(ie)%gdofP(i,j) + Toffset(k)
#endif
                
                do b = 1,np
                   do a = 1,np
                      
#ifdef LOCAL_IDX
                      ! local column index
                      cidx = 2*np*np*nlev + a + (b-1)*np
#else
                      ! global column index            
                      cidx = fptr%base(ie)%gdofP(a,b) + Psoffset
#endif
                      write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dT_dPs(a,b,i,j,k,ie) 
                      
                   end do
                end do

             end do
          end do
       end do

    end do

    close(fid)

    ! ---------------------------------------------------------------------------
    ! H Block 
    ! ---------------------------------------------------------------------------

    fid = 1000 + iam
    write(fname,'(a,a,a,i2.2,a,i2.2,a)') 'Hblock_',trim(Jname),'_proc',iam,'.txt'
    open(fid, file=fname, form='formatted')

    write(fid,'(a23,4x,a23,4x,a23)') "Row Index","Col Index","Value" 

    ! dT_dT
    do ie = 1,nelemd

       do k = 1,nlev
          do j = 1,np
             do i = 1,np

#ifdef LOCAL_IDX
                ! local row index
                ridx = 2*np*np*nlev + np*np + i + (j-1)*np + (k-1)*np*np
#else
                ! global row index            
                ridx = fptr%base(ie)%gdofP(i,j) + Toffset(k)
#endif
                
                do c = 1,nlev
                   do b = 1,np
                      do a = 1,np
                         
#ifdef LOCAL_IDX
                         ! local column index
                         cidx = 2*np*np*nlev + np*np + a + (b-1)*np + (c-1)*np*np
#else
                         ! global column index            
                         cidx = fptr%base(ie)%gdofP(a,b) + Toffset(c)
#endif
                         write(fid,'(i23,4x,i23,4x,e23.16)') ridx, cidx, dT_dT(a,b,c,i,j,k,ie) 
                         
                      end do
                   end do
                end do

             end do
          end do
       end do

    end do

    close(fid)

  end subroutine print_analytic_jacobian
  
  
  
  ! =============================================================================
  ! The following subroutines compute entries for blocks of the Jacobian matrix
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



  subroutine Jac_TimeDeriv(fptr)
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the time derivative
    ! terms in the velocity, surface pressure, and temperature equations
    !----------------------------------------------------------------------------

    use kinds, only                 : real_kind
    use control_mod, only           : tstep_type
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use prim_derived_type_mod, only : derived_type

    implicit none

    !---------------------------------Arguments----------------------------------
    type(derived_type), pointer, intent(in) :: fptr
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------  
    real(kind=real_kind) :: dti  ! inverse of time step size

    integer :: i, j, k, ie ! equation index
    !----------------------------------------------------------------------------

    call t_startf('Jac_TimeDeriv')

    dti = 1.0d0 / fptr%dt

    if (tstep_type==12) then
       dti = (3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2 on step size
    endif

    do ie = 1,nelemd
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

                ! a = i, b = j, c = k
                dU_dU(i,j,k,i,j,k,ie) = dU_dU(i,j,k,i,j,k,ie) + dti
                dV_dV(i,j,k,i,j,k,ie) = dV_dV(i,j,k,i,j,k,ie) + dti

                ! a = i, b = j, c = k
                dT_dT(i,j,k,i,j,k,ie) = dT_dT(i,j,k,i,j,k,ie) + dti

             end do
          end do
       end do
    end do

    do ie = 1,nelemd
       do j = 1,np
          do i = 1,np

             ! a = i, b = j
             dPs_dPs(i,j,i,j,ie) = dPs_dPs(i,j,i,j,ie) + dti

          end do
       end do
    end do

    call t_stopf('Jac_TimeDeriv')

  end subroutine Jac_TimeDeriv



  subroutine Jac_Vorticity(fptr)
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the vorcitity x velocity
    ! terms in the velocity, surface pressure, and temperature equations
    !----------------------------------------------------------------------------
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth
    use derivative_mod, only        : vorticity_sphere

    implicit none

    !---------------------------------Arguments----------------------------------
    type(derived_type), pointer, intent(in) :: fptr
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    real(kind=real_kind), dimension(np,np,nlev) :: zeta

    real(kind=real_kind) :: temp1, temp2, temp3

    integer :: np1
    integer :: a,b
    integer :: i,j,k,ie    ! loop counters
    !----------------------------------------------------------------------------

    call t_startf('Jac_Vorticity')    

    np1 = fptr%tl%np1

    do ie = 1,nelemd

       ! vorticity 
       do k = 1,nlev
          zeta(:,:,k) = vorticity_sphere(fptr%base(ie)%state%v(:,:,:,k,np1), fptr%deriv, fptr%base(ie))
       end do

       ! ---------------------------------------------------------------
       ! Jacobian entries
       ! ---------------------------------------------------------------
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

                temp1 = rrearth * fptr%base(ie)%rmetdet(i,j)

                ! a /= i, b = j, c = k (contribution from x-derivative)
                ! ---------------------------------------------------------------
                do a = 1,np
                   if (a == i) cycle

                   temp2 = temp1 * fptr%base(ie)%D(a,j,1,2) * fptr%deriv%Dvv(a,i) ! for du 
                   temp3 = temp1 * fptr%base(ie)%D(a,j,2,2) * fptr%deriv%Dvv(a,i) ! for dv 

                   dU_dU(a,j,k,i,j,k,ie) = dU_dU(a,j,k,i,j,k,ie) - &
                        temp2 * fptr%base(ie)%state%v(i,j,2,k,np1) ! note - sign <<<

                   dU_dV(a,j,k,i,j,k,ie) = dU_dV(a,j,k,i,j,k,ie) - &
                        temp3 * fptr%base(ie)%state%v(i,j,2,k,np1) ! note - sign <<<

                   dV_dU(a,j,k,i,j,k,ie) = dV_dU(a,j,k,i,j,k,ie) + &
                        temp2 * fptr%base(ie)%state%v(i,j,1,k,np1)

                   dV_dV(a,j,k,i,j,k,ie) = dV_dV(a,j,k,i,j,k,ie) + &
                        temp3 * fptr%base(ie)%state%v(i,j,1,k,np1)
                end do

                ! a = i, b /= j, c = k (contribution from y-derivative)
                ! ---------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle

                   temp2 = -1.0d0 * temp1 * fptr%base(ie)%D(i,b,1,1) * fptr%deriv%Dvv(b,j)
                   temp3 = -1.0d0 * temp1 * fptr%base(ie)%D(i,b,2,1) * fptr%deriv%Dvv(b,j)

                   dU_dU(i,b,k,i,j,k,ie) = dU_dU(i,b,k,i,j,k,ie) - &
                        temp2 * fptr%base(ie)%state%v(i,j,2,k,np1) ! note - sign <<<

                   dU_dV(i,b,k,i,j,k,ie) = dU_dV(i,b,k,i,j,k,ie) - &
                        temp3 * fptr%base(ie)%state%v(i,j,2,k,np1) ! note - sign <<<

                   dV_dU(i,b,k,i,j,k,ie) = dV_dU(i,b,k,i,j,k,ie) + &
                        temp2 * fptr%base(ie)%state%v(i,j,1,k,np1)

                   dV_dV(i,b,k,i,j,k,ie) = dV_dV(i,b,k,i,j,k,ie) + &
                        temp3 * fptr%base(ie)%state%v(i,j,1,k,np1)
                end do

                ! a = i, b = j, c = k (contribution from x and y-derivatives)
                ! ---------------------------------------------------------------
                temp3 = rrearth * fptr%base(ie)%rmetdet(i,j) &
                     * (fptr%base(ie)%D(i,j,1,2) * fptr%deriv%Dvv(i,i) &
                     - fptr%base(ie)%D(i,j,1,1) * fptr%deriv%Dvv(j,j))

                dU_dU(i,j,k,i,j,k,ie) = dU_dU(i,j,k,i,j,k,ie) &
                     - temp3 * fptr%base(ie)%state%v(i,j,2,k,np1) ! note - sign <<<

                dV_dU(i,j,k,i,j,k,ie) = dV_dU(i,j,k,i,j,k,ie) &
                     + temp3 * fptr%base(ie)%state%v(i,j,1,k,np1) + zeta(i,j,k)

                temp3 = rrearth * fptr%base(ie)%rmetdet(i,j) & 
                     * (fptr%base(ie)%D(i,j,2,2) * fptr%deriv%Dvv(i,i) &
                     - fptr%base(ie)%D(i,j,2,1) * fptr%deriv%Dvv(j,j))

                dU_dV(i,j,k,i,j,k,ie) = dU_dV(i,j,k,i,j,k,ie) &  
                     - temp3 * fptr%base(ie)%state%v(i,j,2,k,np1) - zeta(i,j,k) ! note - signS <<<

                dV_dV(i,j,k,i,j,k,ie) = dV_dV(i,j,k,i,j,k,ie) &
                     + temp3 * fptr%base(ie)%state%v(i,j,1,k,np1)

             end do
          end do
       end do

    end do

    call t_stopf('Jac_Vorticity')

  end subroutine Jac_Vorticity



  subroutine Jac_Coriolis(fptr)
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the coriolis x velocity
    ! terms in the velocity, surface pressure, and temperature equations
    !----------------------------------------------------------------------------
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use prim_derived_type_mod, only : derived_type

    implicit none

    !---------------------------------Arguments----------------------------------
    type(derived_type), pointer, intent(in) :: fptr
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    integer :: i, j, k, ie ! equation index
    !----------------------------------------------------------------------------

    call t_startf('Jac_Coriolis')    

    do ie = 1,nelemd

       do k = 1,nlev
          do j = 1,np
             do i = 1,np

                ! a = i, b = j, c = k
                ! ---------------------------------------------------------------
                dU_dV(i,j,k,i,j,k,ie) = dU_dV(i,j,k,i,j,k,ie) - fptr%base(ie)%fcor(i,j) ! note - sign <<<
                dV_dU(i,j,k,i,j,k,ie) = dV_dU(i,j,k,i,j,k,ie) + fptr%base(ie)%fcor(i,j)

             end do
          end do
       end do

    end do

    call t_stopf('Jac_Coriolis')

  end subroutine Jac_Coriolis



  subroutine Jac_Kinetic_Energy(fptr)
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the kinetic energy
    ! terms in the velocity, surface pressure, and temperature equations
    !----------------------------------------------------------------------------
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth

    implicit none

    !---------------------------------Arguments----------------------------------
    type(derived_type), pointer, intent(in) :: fptr
    !----------------------------------------------------------------------------

    !------------------------------Local workspace------------------------------- 
    real(kind=real_kind) :: temp1, temp2

    integer :: np1
    integer :: a, b        ! derivative index
    integer :: i, j, k, ie ! equation index
    !----------------------------------------------------------------------------

    call t_startf('Jac_Kinetic_Energy')    

    np1 = fptr%tl%np1

    do ie = 1,nelemd

       do k = 1,nlev
          do j = 1,np
             do i = 1,np

                ! a /= i, b = j, c = k
                ! ---------------------------------------------------------------
                do a = 1,np
                   if (a == i) cycle

                   temp1 = rrearth * fptr%base(ie)%Dinv(i,j,1,1) * fptr%deriv%Dvv(a,i) ! for u equ
                   temp2 = rrearth * fptr%base(ie)%Dinv(i,j,1,2) * fptr%deriv%Dvv(a,i) ! for v equ

                   dU_dU(a,j,k,i,j,k,ie) = dU_dU(a,j,k,i,j,k,ie) + &
                        temp1 * fptr%base(ie)%state%v(a,j,1,k,np1) 

                   dU_dV(a,j,k,i,j,k,ie) = dU_dV(a,j,k,i,j,k,ie) + &
                        temp1 * fptr%base(ie)%state%v(a,j,2,k,np1)

                   dV_dU(a,j,k,i,j,k,ie) = dV_dU(a,j,k,i,j,k,ie) + &
                        temp2 * fptr%base(ie)%state%v(a,j,1,k,np1)

                   dV_dV(a,j,k,i,j,k,ie) = dV_dV(a,j,k,i,j,k,ie) + &
                        temp2 * fptr%base(ie)%state%v(a,j,2,k,np1)
                end do

                ! a = i, b /= j, c = k
                ! ---------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle

                   temp1 = rrearth * fptr%base(ie)%Dinv(i,j,2,1) * fptr%deriv%Dvv(b,j) ! for u equ
                   temp2 = rrearth * fptr%base(ie)%Dinv(i,j,2,2) * fptr%deriv%Dvv(b,j) ! for v equ

                   dU_dU(i,b,k,i,j,k,ie) = dU_dU(i,b,k,i,j,k,ie) + &
                        temp1 * fptr%base(ie)%state%v(i,b,1,k,np1) 

                   dU_dV(i,b,k,i,j,k,ie) = dU_dV(i,b,k,i,j,k,ie) + &
                        temp1 * fptr%base(ie)%state%v(i,b,2,k,np1)

                   dV_dU(i,b,k,i,j,k,ie) = dV_dU(i,b,k,i,j,k,ie) + &
                        temp2 * fptr%base(ie)%state%v(i,b,1,k,np1)

                   dV_dV(i,b,k,i,j,k,ie) = dV_dV(i,b,k,i,j,k,ie) + &
                        temp2 * fptr%base(ie)%state%v(i,b,2,k,np1)
                end do

                ! a = i, b = j, c = k
                ! ---------------------------------------------------------------
                temp1 = rrearth * (fptr%base(ie)%Dinv(i,j,1,1) * fptr%deriv%Dvv(i,i) &
                     + fptr%base(ie)%Dinv(i,j,2,1) * fptr%deriv%Dvv(j,j))

                dU_dU(i,j,k,i,j,k,ie) = dU_dU(i,j,k,i,j,k,ie) + &
                     temp1 * fptr%base(ie)%state%v(i,j,1,k,np1) 

                dU_dV(i,j,k,i,j,k,ie) = dU_dV(i,j,k,i,j,k,ie) + &
                     temp1 * fptr%base(ie)%state%v(i,j,2,k,np1) 

                temp2 = rrearth * (fptr%base(ie)%Dinv(i,j,1,2) * fptr%deriv%Dvv(i,i) &
                     + fptr%base(ie)%Dinv(i,j,2,2) * fptr%deriv%Dvv(j,j))

                dV_dU(i,j,k,i,j,k,ie) = dV_dU(i,j,k,i,j,k,ie) + &
                     temp2 * fptr%base(ie)%state%v(i,j,1,k,np1) 

                dV_dV(i,j,k,i,j,k,ie) = dV_dV(i,j,k,i,j,k,ie) + &
                     temp2 * fptr%base(ie)%state%v(i,j,2,k,np1) 

             end do
          end do
       end do

    end do

    call t_stopf('Jac_Kinetic_Energy')

  end subroutine Jac_Kinetic_Energy



  subroutine Jac_Geopotential(fptr)
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the gradient of the
    ! geopotential term in the velocity equations
    !----------------------------------------------------------------------------
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth, rgas
    use physics_mod, only           : Virtual_Temperature

    implicit none

    !---------------------------------Arguments----------------------------------
    type(derived_type), pointer, intent(in) :: fptr ! pointer with state data
    !----------------------------------------------------------------------------

    !------------------------------Local workspace------------------------------- 
    real(kind=real_kind), dimension(np,np,nlev) :: T_v, c1, dpdn, p, p2, Qt
    real(kind=real_kind), dimension(nlev)       :: mBm, dBi

    real(kind=real_kind) :: temp, temp11, temp12, temp21, temp22

    integer :: np1, qn0
    integer :: a, b, c ! derivative index
    integer :: i, j, k ! equation index
    integer :: ie      ! element index
    integer :: m       ! vertical sum
    !----------------------------------------------------------------------------

    call t_startf('Jac_Geopotential')

    np1 = fptr%tl%np1
    qn0 = fptr%n_Q 

    do ie = 1,nelemd

       do k = 1,nlev

          ! vertical derivative of pressure (dp/dn)
          dpdn(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
               - (fptr%hvcoord%hyai(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(k) * fptr%base(ie)%state%ps_v(:,:,np1))       

          ! pressure
          p(:,:,k)  = (fptr%hvcoord%hyam(k) * fptr%hvcoord%ps0 + fptr%hvcoord%hybm(k) * fptr%base(ie)%state%ps_v(:,:,np1))
          p2(:,:,k) = p(:,:,k)**2

          ! -1 * dp/dps
          mBm(k) = -1.0d0 * fptr%hvcoord%hybm(k)

          ! d(dp/dn)/dps
          dBi(k) = fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)
       end do

       ! virtual temperature
       if (qn0 == -1) then
          T_v = fptr%base(ie)%state%T(:,:,:,np1)
          c1  = rgas
       else
          Qt  = fptr%base(ie)%state%Qdp(:,:,:,1,qn0) / dpdn
          T_v = Virtual_Temperature(fptr%base(ie)%state%T(:,:,:,np1), Qt)

          c1 = rgas
          c1 = Virtual_Temperature(c1, Qt)
       end if

       ! Levels: 1 to nlev-1
       ! ------------------------------------------------------------------------
       do k = 1,nlev
          do j = 1,np
             do i = 1,np

                temp11 = rrearth * fptr%base(ie)%Dinv(i,j,1,1)
                temp12 = rrearth * fptr%base(ie)%Dinv(i,j,1,2)

                temp21 = rrearth * fptr%base(ie)%Dinv(i,j,2,1)
                temp22 = rrearth * fptr%base(ie)%Dinv(i,j,2,2)

                ! ===============================================================
                ! surface pressure derivatives 
                ! ===============================================================

                ! a /= i, b = j (x-derivative)
                ! ---------------------------------------------------------------
                do a = 1,np
                   if (a == i) cycle

                   temp = T_v(a,j,k) * 0.5d0 & 
                        * (mBm(k) * dpdn(a,j,k) + dBi(k) * p(a,j,k))/p2(a,j,k) 

                   do m = k+1,nlev
                      temp = temp + T_v(a,j,m) &
                           * (mBm(m) * dpdn(a,j,m) + dBi(m) * p(a,j,m))/p2(a,j,m)
                   end do

                   temp = temp * rgas * fptr%deriv%Dvv(a,i)

                   dU_dPs(a,j,i,j,k,ie) = dU_dPs(a,j,i,j,k,ie) + temp11 * temp
                   dV_dPs(a,j,i,j,k,ie) = dV_dPs(a,j,i,j,k,ie) + temp12 * temp
                end do

                ! a = i, b /= j (y-derivative)
                ! ---------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle

                   temp = T_v(i,b,k) * 0.5d0 &
                        * (mBm(k) * dpdn(i,b,k) + dBi(k) * p(i,b,k))/p2(i,b,k)

                   do m = k+1,nlev
                      temp = temp + T_v(i,b,m) &
                           * (mBm(m) * dpdn(i,b,m) + dBi(m) * p(i,b,m))/p2(i,b,m)
                   end do

                   temp = temp * rgas * fptr%deriv%Dvv(b,j)

                   dU_dPs(i,b,i,j,k,ie) = dU_dPs(i,b,i,j,k,ie) + temp21 * temp
                   dV_dPs(i,b,i,j,k,ie) = dV_dPs(i,b,i,j,k,ie) + temp22 * temp
                end do

                ! a = i, b = j, c = k (x and y-derivatives)
                ! ---------------------------------------------------------------
                temp = T_v(i,j,k) * 0.5d0 &
                     * (mBm(k) * dpdn(i,j,k) + dBi(k) * p(i,j,k))/p2(i,j,k)

                do m = k+1,nlev
                   temp = temp + T_v(i,j,m) &
                        * (mBm(m) * dpdn(i,j,m) + dBi(m) * p(i,j,m))/p2(i,j,m)
                end do

                dU_dPs(i,j,i,j,k,ie) = dU_dPs(i,j,i,j,k,ie) &
                     + rrearth * rgas * temp * (fptr%base(ie)%Dinv(i,j,1,1) * fptr%deriv%Dvv(i,i) &
                     + fptr%base(ie)%Dinv(i,j,2,1) * fptr%deriv%Dvv(j,j))

                dV_dPs(i,j,i,j,k,ie) = dV_dPs(i,j,i,j,k,ie) &
                     + rrearth * rgas * temp * (fptr%base(ie)%Dinv(i,j,1,2) * fptr%deriv%Dvv(i,i) &
                     + fptr%base(ie)%Dinv(i,j,2,2) * fptr%deriv%Dvv(j,j))

                ! ===============================================================
                ! Temperature Derivatives, c = k and c > k
                ! ===============================================================

                ! a /= i, b = j, c = k (x-derivative)
                ! ---------------------------------------------------------------
                do a = 1,np
                   if (a == i) cycle

                   temp = c1(a,j,k) * 0.5d0 * dpdn(a,j,k) * fptr%deriv%Dvv(a,i)/p(a,j,k)

                   dU_dT(a,j,k,i,j,k,ie) = dU_dT(a,j,k,i,j,k,ie) + temp11 * temp
                   dV_dT(a,j,k,i,j,k,ie) = dV_dT(a,j,k,i,j,k,ie) + temp12 * temp
                end do

                ! a = i, b /= j, c = k (y-derivative)
                ! ---------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle
                   
                   temp = c1(i,b,k) * 0.5d0 * dpdn(i,b,k) * fptr%deriv%Dvv(b,j)/p(i,b,k)

                   dU_dT(i,b,k,i,j,k,ie) = dU_dT(i,b,k,i,j,k,ie) + temp21 * temp
                   dV_dT(i,b,k,i,j,k,ie) = dV_dT(i,b,k,i,j,k,ie) + temp22 * temp
                end do

                ! a = i, b = j, c = k (x and y-derivatives)
                ! ---------------------------------------------------------------
                temp = rrearth * c1(i,j,k) * 0.5d0 * dpdn(i,j,k)/p(i,j,k)

                dU_dT(i,j,k,i,j,k,ie) = dU_dT(i,j,k,i,j,k,ie)          &
                     + temp * (fptr%base(ie)%Dinv(i,j,1,1) * fptr%deriv%Dvv(i,i) &
                     + fptr%base(ie)%Dinv(i,j,2,1) * fptr%deriv%Dvv(j,j))

                dV_dT(i,j,k,i,j,k,ie) = dV_dT(i,j,k,i,j,k,ie)          &
                     + temp * (fptr%base(ie)%Dinv(i,j,1,2) * fptr%deriv%Dvv(i,i) &
                     + fptr%base(ie)%Dinv(i,j,2,2) * fptr%deriv%Dvv(j,j))

                do c = k+1,nlev

                   ! a /= i, b = j, c > k
                   ! ------------------------------------------------------------
                   do a = 1,np
                      if (a == i) cycle

                      temp = c1(a,j,c) * dpdn(a,j,c) * fptr%deriv%Dvv(a,i)/p(a,j,c)

                      dU_dT(a,j,c,i,j,k,ie) = dU_dT(a,j,c,i,j,k,ie) + temp11 * temp
                      dV_dT(a,j,c,i,j,k,ie) = dV_dT(a,j,c,i,j,k,ie) + temp12 * temp
                   end do

                   ! a = i, b /= j, c > k
                   !-------------------------------------------------------------
                   do b = 1,np
                      if (b == j) cycle

                      temp = c1(i,b,c) * dpdn(i,b,c) * fptr%deriv%Dvv(b,j)/p(i,b,c)

                      dU_dT(i,b,c,i,j,k,ie) = dU_dT(i,b,c,i,j,k,ie) + temp21 * temp
                      dV_dT(i,b,c,i,j,k,ie) = dV_dT(i,b,c,i,j,k,ie) + temp22 * temp
                   end do

                   ! a = i, b = j, c > k
                   ! ------------------------------------------------------------
                   temp = rrearth * c1(i,j,c) * dpdn(i,j,c)/p(i,j,c)

                   dU_dT(i,j,c,i,j,k,ie) = dU_dT(i,j,c,i,j,k,ie)          &
                        + temp * (fptr%base(ie)%Dinv(i,j,1,1) * fptr%deriv%Dvv(i,i) &
                        + fptr%base(ie)%Dinv(i,j,2,1) * fptr%deriv%Dvv(j,j))

                   dV_dT(i,j,c,i,j,k,ie) = dV_dT(i,j,c,i,j,k,ie)          &
                        + temp * (fptr%base(ie)%Dinv(i,j,1,2) * fptr%deriv%Dvv(i,i) &
                        + fptr%base(ie)%Dinv(i,j,2,2) * fptr%deriv%Dvv(j,j))
                end do

             end do
          end do
       end do

       ! Level K
       !-------------------------------------------------------------------------
       ! k = nlev
       ! do j = 1,np
       !    do i = 1,np

       !       temp11 = rrearth * fptr%base(ie)%Dinv(i,j,1,1) ! u equation
       !       temp12 = rrearth * fptr%base(ie)%Dinv(i,j,1,2) ! v equation   

       !       temp21 = rrearth * fptr%base(ie)%Dinv(i,j,2,1) ! u equation 
       !       temp22 = rrearth * fptr%base(ie)%Dinv(i,j,2,2) ! v equation  

       !       ! ===============================================================
       !       ! surface pressure derivatives 
       !       ! ===============================================================

       !       ! a /= i, b = j
       !       ! ------------------------------------------------------------
       !       do a = 1,np
       !          if (a == i) cycle

       !          temp = 0.5d0 * rgas * T_v(a,j,k) * fptr%deriv%Dvv(a,i) &
       !               * (mBm(k) * dpdn(a,j,k)/p2(a,j,k) + dBi(k)/p(a,j,k)) 

       !          dU_dPs(a,j,i,j,k,ie) = dU_dPs(a,j,i,j,k,ie) + temp11 * temp
       !          dV_dPs(a,j,i,j,k,ie) = dV_dPs(a,j,i,j,k,ie) + temp12 * temp
       !       end do

       !       ! a = i, b /= j
       !       ! ------------------------------------------------------------
       !       do b = 1,np
       !          if (b == j) cycle

       !          temp = 0.5d0* rgas * T_v(i,b,k) * fptr%deriv%Dvv(b,j) &
       !               * (mBm(k) * dpdn(i,b,k)/p2(i,b,k) + dBi(k)/p(i,b,k)) 

       !          dU_dPs(i,b,i,j,k,ie) = dU_dPs(i,b,i,j,k,ie) + temp21 * temp
       !          dV_dPs(i,b,i,j,k,ie) = dV_dPs(i,b,i,j,k,ie) + temp22 * temp
       !       end do

       !       ! a = i, b = j
       !       ! ------------------------------------------------------------
       !       temp = rrearth * rgas * T_v(i,j,k) * 0.5d0 &
       !            * (mBm(k) * dpdn(i,j,k)/p2(i,j,k) + dBi(k)/p(i,j,k)) 

       !       dU_dPs(i,j,i,j,k,ie) = dU_dPs(i,j,i,j,k,ie)            & 
       !            + temp * (fptr%base(ie)%Dinv(i,j,1,1) * fptr%deriv%Dvv(i,i) & 
       !            + fptr%base(ie)%Dinv(i,j,2,1) * fptr%deriv%Dvv(j,j))

       !       dV_dPs(i,j,i,j,k,ie) = dV_dPs(i,j,i,j,k,ie)            & 
       !            + temp * (fptr%base(ie)%Dinv(i,j,1,2) * fptr%deriv%Dvv(i,i) &
       !            + fptr%base(ie)%Dinv(i,j,2,2) * fptr%deriv%Dvv(j,j))

       !       ! ===============================================================
       !       ! temperature derivatives 
       !       ! ===============================================================

       !       ! a /= i, b = j, c = k
       !       ! ------------------------------------------------------------
       !       do a = 1,np
       !          if (a == i) cycle

       !          temp = c1(a,j,k) * 0.5d0 * dpdn(a,j,k) * fptr%deriv%Dvv(a,i)/p(a,j,k)

       !          dU_dT(a,j,k,i,j,k,ie) = dU_dT(a,j,k,i,j,k,ie) + temp11 * temp
       !          dV_dT(a,j,k,i,j,k,ie) = dV_dT(a,j,k,i,j,k,ie) + temp12 * temp
       !       end do

       !       ! a = i, b /= j, c = k
       !       ! ------------------------------------------------------------
       !       do b = 1,np
       !          if (b == j) cycle
                
       !          temp = c1(i,b,k) * 0.5d0 * dpdn(i,b,k) * fptr%deriv%Dvv(b,j)/p(i,b,k)

       !          dU_dT(i,b,k,i,j,k,ie) = dU_dT(i,b,k,i,j,k,ie) + temp21 * temp
       !          dV_dT(i,b,k,i,j,k,ie) = dV_dT(i,b,k,i,j,k,ie) + temp22 * temp
       !       end do

       !       ! a = i, b = j, c = k
       !       ! ------------------------------------------------------------
       !       temp = rrearth * c1(i,j,k) * 0.5d0 * dpdn(i,j,k)/p(i,j,k)

       !       dU_dT(i,j,k,i,j,k,ie) = dU_dT(i,j,k,i,j,k,ie)          &
       !            + temp * (fptr%base(ie)%Dinv(i,j,1,1) * fptr%deriv%Dvv(i,i) &
       !            + fptr%base(ie)%Dinv(i,j,2,1) * fptr%deriv%Dvv(j,j))

       !       dV_dT(i,j,k,i,j,k,ie) = dV_dT(i,j,k,i,j,k,ie)          &
       !            + temp * (fptr%base(ie)%Dinv(i,j,1,2) * fptr%deriv%Dvv(i,i) &
       !            + fptr%base(ie)%Dinv(i,j,2,2) * fptr%deriv%Dvv(j,j))

       !    end do
       ! end do

    end do

    call t_stopf('Jac_Geopotential')

  end subroutine Jac_Geopotential



  subroutine Jac_press_grad(fptr)
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the pressure gradient
    ! term in the velocity equations
    !----------------------------------------------------------------------------

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use prim_derived_type_mod, only : derived_type
    use derivative_mod, only        : gradient_sphere
    use physical_constants, only    : rrearth, Rgas
    use physics_mod, only           : Virtual_Temperature

    implicit none

    !---------------------------------Arguments----------------------------------
    type(derived_type), pointer, intent(in) :: fptr ! pointer with state data
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    real(kind=real_kind), dimension(np,np,nlev)   :: T_v, p, p2, dp, Qt, temp3
    real(kind=real_kind), dimension(np,np,2,nlev) :: grad_p
    real(kind=real_kind), dimension(np,np,2)      :: grad_ps

    real(kind=real_kind) :: temp1, temp2

    integer :: a,b      ! derivative index
    integer :: h,k,l,ie ! equation index
    integer :: np1, qn0 ! time level 
    !----------------------------------------------------------------------------

    call t_startf('Jac_press_grad')

    np1  = fptr%tl%np1
    qn0 = fptr%n_Q

    do ie = 1,nelemd

       ! virtual temperature
       ! ------------------------------------------------------------------------
       if (qn0 == -1 ) then
          T_v(:,:,:)   = fptr%base(ie)%state%T(:,:,:,np1)
          temp3(:,:,:) = Rgas
       else
          do l = 1,nlev
             dp(:,:,l) = (fptr%hvcoord%hyai(l+1)*fptr%hvcoord%ps0 + fptr%hvcoord%hybi(l+1)*fptr%base(ie)%state%ps_v(:,:,np1)) &
                  - (fptr%hvcoord%hyai(l)*fptr%hvcoord%ps0 + fptr%hvcoord%hybi(l)*fptr%base(ie)%state%ps_v(:,:,np1))
          end do

          Qt(:,:,:)  = fptr%base(ie)%state%Qdp(:,:,:,1,qn0)/dp(:,:,:)
          T_v(:,:,:) = Virtual_Temperature(fptr%base(ie)%state%T(:,:,:,np1),Qt)

          temp3(:,:,:) = Rgas
          temp3(:,:,:) = Virtual_temperature(temp3, Qt)
       end if

       ! gradient of surface pressure, pressure, 1/pressure, and gradient of pressure
       ! ------------------------------------------------------------------------
       grad_ps = gradient_sphere(fptr%base(ie)%state%ps_v(:,:,np1), fptr%deriv, fptr%base(ie)%Dinv)

       do l = 1,nlev

          p(:,:,l)  = (fptr%hvcoord%hyam(l) * fptr%hvcoord%ps0 + fptr%hvcoord%hybm(l) * fptr%base(ie)%state%ps_v(:,:,np1))
          p2(:,:,l) = p(:,:,l)**2

          grad_p(:,:,:,l) = fptr%hvcoord%hybm(l) * grad_ps(:,:,:)
       end do

       ! ------------------------------------------------------------------------
       ! Jacobian entries
       ! ------------------------------------------------------------------------
       do l = 1,nlev
          do k = 1,np
             do h = 1,np

                ! ---------------------------------------------------------------
                ! surface pressure derivatives
                ! ---------------------------------------------------------------

                temp1 = fptr%hvcoord%hybm(l) * Rgas * T_v(h,k,l) * rrearth * fptr%base(ie)%Dinv(h,k,1,1) / p(h,k,l)
                temp2 = fptr%hvcoord%hybm(l) * Rgas * T_v(h,k,l) * rrearth * fptr%base(ie)%Dinv(h,k,1,2) / p(h,k,l)

                ! a /= h, b = k
                b = k
                do a = 1,np
                   if (a == h) cycle
                   dU_dPs(a,b,h,k,l,ie) = dU_dPs(a,b,h,k,l,ie) + temp1 * fptr%deriv%Dvv(a,h)
                   dV_dPs(a,b,h,k,l,ie) = dV_dPs(a,b,h,k,l,ie) + temp2 * fptr%deriv%Dvv(a,h)
                end do

                temp1 = fptr%hvcoord%hybm(l) * Rgas * T_v(h,k,l) * rrearth * fptr%base(ie)%Dinv(h,k,2,1) / p(h,k,l)
                temp2 = fptr%hvcoord%hybm(l) * Rgas * T_v(h,k,l) * rrearth * fptr%base(ie)%Dinv(h,k,2,2) / p(h,k,l)

                ! a = h, b /= k
                a = h 
                do b = 1,np
                   if (b == k) cycle  
                   dU_dPs(a,b,h,k,l,ie) = dU_dPs(a,b,h,k,l,ie) + temp1 * fptr%deriv%Dvv(b,k)
                   dV_dPs(a,b,h,k,l,ie) = dV_dPs(a,b,h,k,l,ie) + temp2 * fptr%deriv%Dvv(b,k)
                end do

                ! a = h, b = k
                a = h; b = k

                temp1 = -1.0d0 * fptr%hvcoord%hybm(l) * Rgas * T_v(h,k,l) / p2(h,k,l)
                temp2 = fptr%hvcoord%hybm(l) * Rgas * T_v(h,k,l) * rrearth / p(h,k,l)

                dU_dPs(a,b,h,k,l,ie) = dU_dPs(a,b,h,k,l,ie) + temp1 * grad_p(h,k,1,l) &
                     + temp2 * (fptr%base(ie)%Dinv(h,k,1,1) * fptr%deriv%Dvv(a,h) + fptr%base(ie)%Dinv(h,k,2,1) * fptr%deriv%Dvv(b,k))

                dV_dPs(a,b,h,k,l,ie) = dV_dPs(a,b,h,k,l,ie) + temp1 * grad_p(h,k,2,l) &
                     + temp2 * (fptr%base(ie)%Dinv(h,k,1,2) * fptr%deriv%Dvv(a,h) + fptr%base(ie)%Dinv(h,k,2,2) * fptr%deriv%Dvv(b,k))

                ! ---------------------------------------------------------------
                ! Temperature derivatives
                ! ---------------------------------------------------------------

                ! a = h, b = k, c = l
                dU_dT(h,k,l,h,k,l,ie) = dU_dT(h,k,l,h,k,l,ie) + temp3(h,k,l) * grad_p(h,k,1,l) / p(h,k,l)
                dV_dT(h,k,l,h,k,l,ie) = dV_dT(h,k,l,h,k,l,ie) + temp3(h,k,l) * grad_p(h,k,2,l) / p(h,k,l)

             end do
          end do
       end do

    end do

    call t_stopf('Jac_press_grad')

  end subroutine Jac_press_grad



  subroutine Jac_VertAdv(fptr)
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the vertical advection 
    ! terms in the velocity and temperature equations
    !----------------------------------------------------------------------------

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nelemd
    use prim_derived_type_mod, only : derived_type
    use derivative_mod, only        : divergence_sphere

    implicit none

    !---------------------------------Arguments----------------------------------
    type(derived_type), pointer, intent(in) :: fptr
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    real(kind=real_kind), dimension(np,np,nlev,np,np,nlev+1) :: du_eta_dot_dpdn, dv_eta_dot_dpdn
    real(kind=real_kind), dimension(np,np,np,np,nlev+1)      :: dps_eta_dot_dpdn
    real(kind=real_kind), dimension(np,np,nlev+1)            :: eta_dot_dpdn
    real(kind=real_kind), dimension(np,np,nlev)              :: dpdn, rdpdn, half_rdpdn, div_dpdn_vel
    real(kind=real_kind), dimension(np,np,2)                 :: dpdn_vel   
    real(kind=real_kind), dimension(np,np)                   :: sdot_sum(np,np)

    real(kind=real_kind) :: half_rdpdn_diffu1, half_rdpdn_diffu2
    real(kind=real_kind) :: half_rdpdn_diffv1, half_rdpdn_diffv2
    real(kind=real_kind) :: half_rdpdn_diffT1, half_rdpdn_diffT2
    real(kind=real_kind) :: temp1, temp2, diffB

    integer :: np1
    integer :: a,b,c
    integer :: h,k,l,ie
    !----------------------------------------------------------------------------

    call t_startf('Jac_VertAdv')

    np1 = fptr%tl%np1    

    do ie = 1,nelemd

       ! vertical derivative of pressure, dpdn
       do l = 1,nlev
          dpdn(:,:,l) = (fptr%hvcoord%hyai(l+1)*fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybi(l+1)*fptr%base(ie)%state%ps_v(:,:,np1)) &
               - (fptr%hvcoord%hyai(l)*fptr%hvcoord%ps0 &
               + fptr%hvcoord%hybi(l)*fptr%base(ie)%state%ps_v(:,:,np1))

          rdpdn(:,:,l) = 1.0d0 / dpdn(:,:,l)

          half_rdpdn(:,:,l) = 0.5d0 / dpdn(:,:,l)
       end do

       ! divergence of dpdn * velocity, needed for eta_dot_dpdn
       do l = 1,nlev
          do k=1,np
             do h=1,np             
                dpdn_vel(h,k,1) = dpdn(h,k,l) * fptr%base(ie)%state%v(h,k,1,l,np1)
                dpdn_vel(h,k,2) = dpdn(h,k,l) * fptr%base(ie)%state%v(h,k,2,l,np1)
             end do
          end do

          div_dpdn_vel(:,:,l) = divergence_sphere(dpdn_vel, fptr%deriv, fptr%base(ie))
       end do

       ! compute eta_dot_dpdn
       ! ------------------------------------------------------------------------
       sdot_sum = 0.0d0

       do l=1,nlev
          sdot_sum(:,:)         = sdot_sum(:,:) + div_dpdn_vel(:,:,l)
          eta_dot_dpdn(:,:,l+1) = sdot_sum(:,:)
       end do

       do l=1,nlev-1
          eta_dot_dpdn(:,:,l+1) = fptr%hvcoord%hybi(l+1)*sdot_sum(:,:) &
               - eta_dot_dpdn(:,:,l+1)
       end do

       eta_dot_dpdn(:,:,1)      = 0.0d0
       eta_dot_dpdn(:,:,nlev+1) = 0.0d0

       ! compute the Jacobian of etadot_dpdn
       call Jac_etadot_dpdn(du_eta_dot_dpdn, dv_eta_dot_dpdn, dps_eta_dot_dpdn, fptr, ie) 

       ! ------------------------------------------------------------------------
       ! Level 1 Jacobian entries
       ! ------------------------------------------------------------------------
       l = 1
       do k = 1,np
          do h = 1,np

             half_rdpdn_diffu1 = half_rdpdn(h,k,l) * &
                  (fptr%base(ie)%state%v(h,k,1,l+1,np1) &
                  - fptr%base(ie)%state%v(h,k,1,l,np1))

             half_rdpdn_diffv1 = half_rdpdn(h,k,l) * &
                  (fptr%base(ie)%state%v(h,k,2,l+1,np1) &
                  - fptr%base(ie)%state%v(h,k,2,l,np1))

             half_rdpdn_diffT1 = half_rdpdn(h,k,l) * &
                  (fptr%base(ie)%state%T(h,k,l+1,np1) &
                  - fptr%base(ie)%state%T(h,k,l,np1))

             ! ==================================================================
             ! velocity derivatives
             ! ==================================================================
             do c = 1,nlev

                ! a /= h, b = k, c = * (skips zeros in dvel_eta_dot_dpdn)
                ! ------------------------------------------------------------
                b = k
                do a = 1,np
                   if (a == h) cycle

                   dU_dU(a,b,c,h,k,l,ie) = dU_dU(a,b,c,h,k,l,ie) &
                        + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1

                   dU_dV(a,b,c,h,k,l,ie) = dU_dV(a,b,c,h,k,l,ie) &
                        + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1

                   dV_dU(a,b,c,h,k,l,ie) = dV_dU(a,b,c,h,k,l,ie) &
                        + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1

                   dV_dV(a,b,c,h,k,l,ie) = dV_dV(a,b,c,h,k,l,ie) &
                        + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1

                   dT_dU(a,b,c,h,k,l,ie)  = dT_dU(a,b,c,h,k,l,ie) &
                        + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffT1

                   dT_dV(a,b,c,h,k,l,ie)  = dT_dV(a,b,c,h,k,l,ie) &
                        + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffT1
                end do

                ! a = h, b /= k, c = * (skips zeros in dvel_eta_dot_dpdn)
                ! ------------------------------------------------------------
                a = h
                do b = 1,np
                   if (b == k) cycle

                   dU_dU(a,b,c,h,k,l,ie) = dU_dU(a,b,c,h,k,l,ie) &
                        + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1

                   dU_dV(a,b,c,h,k,l,ie) = dU_dV(a,b,c,h,k,l,ie) &
                        + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1

                   dV_dU(a,b,c,h,k,l,ie) = dV_dU(a,b,c,h,k,l,ie) &
                        + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1

                   dV_dV(a,b,c,h,k,l,ie) = dV_dV(a,b,c,h,k,l,ie) &
                        + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1

                   dT_dU(a,b,c,h,k,l,ie)  = dT_dU(a,b,c,h,k,l,ie)  &
                        + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffT1

                   dT_dV(a,b,c,h,k,l,ie)  = dT_dV(a,b,c,h,k,l,ie)  &
                        + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffT1
                end do

                ! a = h, b = k
                ! ------------------------------------------------------------
                a = h; b = k

                dU_dV(a,b,c,h,k,l,ie) = dU_dV(a,b,c,h,k,l,ie) &
                     + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1

                dV_dU(a,b,c,h,k,l,ie) = dV_dU(a,b,c,h,k,l,ie) &
                     + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1

                if (c /= l .and. c /= l+1) then               
                   dU_dU(a,b,c,h,k,l,ie) = dU_dU(a,b,c,h,k,l,ie) &
                        + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1

                   dV_dV(a,b,c,h,k,l,ie) = dV_dV(a,b,c,h,k,l,ie) &
                        + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1
                end if

                dT_dU(a,b,c,h,k,l,ie)  = dT_dU(a,b,c,h,k,l,ie)  &
                     + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffT1

                dT_dV(a,b,c,h,k,l,ie)  = dT_dV(a,b,c,h,k,l,ie)  &
                     + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffT1
             end do

             ! special cases
             temp1 = half_rdpdn(h,k,l) * eta_dot_dpdn(h,k,l+1) ! eta_dot_dpdn l+1/2

             ! a = h, b = k, c = l
             ! ------------------------------------------------------------
             a = h; b = k; c = l

             dU_dU(a,b,c,h,k,l,ie) = dU_dU(a,b,c,h,k,l,ie) &
                  + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1 - temp1

             dV_dV(a,b,c,h,k,l,ie) = dV_dV(a,b,c,h,k,l,ie) &
                  + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1 - temp1

             ! a = h, b = k, c = l+1
             ! ------------------------------------------------------------
             a = h; b = k; c = l+1

             dU_dU(a,b,c,h,k,l,ie) = dU_dU(a,b,c,h,k,l,ie) &
                  + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1 + temp1

             dV_dV(a,b,c,h,k,l,ie) = dV_dV(a,b,c,h,k,l,ie) &
                  + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1 + temp1

             ! ==================================================================
             ! surface pressure derivatives
             ! ==================================================================

             ! a /= h, b = k (skips zeros in dvel_eta_dot_dpdn)
             ! ------------------------------------------------------------
             b = k
             do a = 1,np
                if (a == h) cycle

                dU_dPs(a,b,h,k,l,ie) = dU_dPs(a,b,h,k,l,ie) &
                     + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffu1

                dV_dPs(a,b,h,k,l,ie) = dV_dPs(a,b,h,k,l,ie) &
                     + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffv1

                dT_dPs(a,b,h,k,l,ie)  = dT_dPs(a,b,h,k,l,ie)  &
                     + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffT1
             end do

             ! a = h, b /= k (skips zeros in dvel_eta_dot_dpdn)
             ! ------------------------------------------------------------
             a = h
             do b = 1,np
                if (b == k) cycle

                dU_dPs(a,b,h,k,l,ie) = dU_dPs(a,b,h,k,l,ie) &
                     + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffu1

                dV_dPs(a,b,h,k,l,ie) = dV_dPs(a,b,h,k,l,ie) &
                     + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffv1

                dT_dPs(a,b,h,k,l,ie)  = dT_dPs(a,b,h,k,l,ie)  &
                     + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffT1
             end do

             ! a = h, b = k
             ! ------------------------------------------------------------
             a = h; b = k 

             diffB = fptr%hvcoord%hybi(l+1) - fptr%hvcoord%hybi(l)

             temp1 = -0.5d0 * diffB * rdpdn(h,k,l)**2 * eta_dot_dpdn(h,k,l+1)

             dU_dPs(a,b,h,k,l,ie) = dU_dPs(a,b,h,k,l,ie)    &
                  + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffu1 &
                  + temp1 * (fptr%base(ie)%state%v(h,k,1,l+1,np1) &
                  - fptr%base(ie)%state%v(h,k,1,l,np1))

             dV_dPs(a,b,h,k,l,ie) = dV_dPs(a,b,h,k,l,ie)    &
                  + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffv1 &
                  + temp1 * (fptr%base(ie)%state%v(h,k,2,l+1,np1) &
                  - fptr%base(ie)%state%v(h,k,2,l,np1))

             dT_dPs(a,b,h,k,l,ie) = dT_dPs(a,b,h,k,l,ie)      &
                  + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffT1 &
                  + temp1 * (fptr%base(ie)%state%T(h,k,l+1,np1) &
                  - fptr%base(ie)%state%T(h,k,l,np1))

             ! ==================================================================
             ! temperature derivatives
             ! ==================================================================

             temp1 = half_rdpdn(h,k,l) * eta_dot_dpdn(h,k,l+1) ! eta_dot_dpdn l+1/2

             ! a = h, b = k, c = l
             dT_dT(h,k,l,h,k,l,ie) = dT_dT(h,k,l,h,k,l,ie) - temp1

             ! a = h, b = k, c = l+1
             dT_dT(h,k,l+1,h,k,l,ie) = dT_dT(h,k,l+1,h,k,l,ie) + temp1

          end do
       end do

       ! ------------------------------------------------------------------------
       ! Level 2 to nlev-1 Jacobian Entries
       ! ------------------------------------------------------------------------
       do l=2,nlev-1
          do k=1,np
             do h=1,np

                half_rdpdn_diffu1 = half_rdpdn(h,k,l) &
                     * (fptr%base(ie)%state%v(h,k,1,l+1,np1) &
                     - fptr%base(ie)%state%v(h,k,1,l,np1))

                half_rdpdn_diffu2 = half_rdpdn(h,k,l) &
                     * (fptr%base(ie)%state%v(h,k,1,l,np1) &
                     - fptr%base(ie)%state%v(h,k,1,l-1,np1))

                half_rdpdn_diffv1 = half_rdpdn(h,k,l) &
                     * (fptr%base(ie)%state%v(h,k,2,l+1,np1) &
                     - fptr%base(ie)%state%v(h,k,2,l,np1))

                half_rdpdn_diffv2 = half_rdpdn(h,k,l) &
                     * (fptr%base(ie)%state%v(h,k,2,l,np1) &
                     - fptr%base(ie)%state%v(h,k,2,l-1,np1))

                half_rdpdn_diffT1 = half_rdpdn(h,k,l) &
                     * (fptr%base(ie)%state%T(h,k,l+1,np1) &
                     - fptr%base(ie)%state%T(h,k,l,np1))

                half_rdpdn_diffT2 = half_rdpdn(h,k,l) &
                     * (fptr%base(ie)%state%T(h,k,l,np1) &
                     - fptr%base(ie)%state%T(h,k,l-1,np1))

                ! ===============================================================
                ! velocity derivatives
                ! ===============================================================
                do c = 1,nlev

                   ! a /= h, b = k, c = * (skips zeros in dvel_eta_dot_dpdn)
                   ! ------------------------------------------------------------
                   b = k
                   do a = 1,np
                      if (a == h) cycle

                      dU_dU(a,b,c,h,k,l,ie) = dU_dU(a,b,c,h,k,l,ie)   &
                           + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1 &
                           + du_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffu2

                      dU_dV(a,b,c,h,k,l,ie) = dU_dV(a,b,c,h,k,l,ie)   &
                           + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1 &
                           + dv_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffu2

                      dV_dU(a,b,c,h,k,l,ie) = dV_dU(a,b,c,h,k,l,ie)   &
                           + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1 &
                           + du_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffv2

                      dV_dV(a,b,c,h,k,l,ie) = dV_dV(a,b,c,h,k,l,ie)   &
                           + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1 &
                           + dv_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffv2

                      dT_dU(a,b,c,h,k,l,ie)  = dT_dU(a,b,c,h,k,l,ie)    &
                           + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffT1 &
                           + du_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffT2

                      dT_dV(a,b,c,h,k,l,ie)  = dT_dV(a,b,c,h,k,l,ie)    &
                           + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffT1 &
                           + dv_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffT2
                   end do

                   ! a = h, b /= k, c = * (skips zeros in dvel_eta_dot_dpdn)
                   ! ------------------------------------------------------------
                   a = h
                   do b = 1,np
                      if (b == k) cycle

                      dU_dU(a,b,c,h,k,l,ie) = dU_dU(a,b,c,h,k,l,ie)   &
                           + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1 &
                           + du_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffu2

                      dU_dV(a,b,c,h,k,l,ie) = dU_dV(a,b,c,h,k,l,ie)   &
                           + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1 &
                           + dv_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffu2

                      dV_dU(a,b,c,h,k,l,ie) = dV_dU(a,b,c,h,k,l,ie)   &
                           + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1 &
                           + du_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffv2

                      dV_dV(a,b,c,h,k,l,ie) = dV_dV(a,b,c,h,k,l,ie)   &
                           + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1 &
                           + dv_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffv2

                      dT_dU(a,b,c,h,k,l,ie)  = dT_dU(a,b,c,h,k,l,ie)    &
                           + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffT1 &
                           + du_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffT2

                      dT_dV(a,b,c,h,k,l,ie)  = dT_dV(a,b,c,h,k,l,ie)    &
                           + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffT1 &
                           + dv_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffT2
                   end do

                   ! a = h, b = k
                   ! ------------------------------------------------------------
                   a = h; b = k

                   dU_dV(a,b,c,h,k,l,ie) = dU_dV(a,b,c,h,k,l,ie)   &
                        + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1 &
                        + dv_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffu2

                   dV_dU(a,b,c,h,k,l,ie) = dV_dU(a,b,c,h,k,l,ie)   &
                        + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1 &
                        + du_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffv2

                   if ( c /= l-1 .and. c /= l .and. c /= l+1 ) then               
                      dU_dU(a,b,c,h,k,l,ie) = dU_dU(a,b,c,h,k,l,ie)   &
                           + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1 &
                           + du_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffu2

                      dV_dV(a,b,c,h,k,l,ie) = dV_dV(a,b,c,h,k,l,ie)   & 
                           + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1 &
                           + dv_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffv2
                   end if

                   dT_dU(a,b,c,h,k,l,ie) = dT_dU(a,b,c,h,k,l,ie)     &
                        + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffT1 &
                        + du_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffT2

                   dT_dV(a,b,c,h,k,l,ie) = dT_dV(a,b,c,h,k,l,ie)     &
                        + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffT1 &
                        + dv_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffT2
                end do

                ! special cases
                temp1 = half_rdpdn(h,k,l) * eta_dot_dpdn(h,k,l+1) ! eta_dot_dpdn l+1/2
                temp2 = half_rdpdn(h,k,l) * eta_dot_dpdn(h,k,l)   ! eta_dot_dpdn l-1/2

                ! a = h, b = k, c = l-1
                ! ------------------------------------------------------------
                a = h; b = k; c = l-1

                dU_dU(a,b,c,h,k,l,ie) = dU_dU(a,b,c,h,k,l,ie)   &
                     + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1 &
                     + du_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffu2 - temp2

                dV_dV(a,b,c,h,k,l,ie) = dV_dV(a,b,c,h,k,l,ie)   &
                     + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1 &
                     + dv_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffv2 - temp2

                ! a = h, b = k, c = l
                ! ------------------------------------------------------------
                a = h; b = k; c = l

                dU_dU(a,b,c,h,k,l,ie) = dU_dU(a,b,c,h,k,l,ie)   &
                     + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1 &
                     + du_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffu2   &
                     - temp1 + temp2

                dV_dV(a,b,c,h,k,l,ie) = dV_dV(a,b,c,h,k,l,ie)   &
                     + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1 &
                     + dv_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffv2   &
                     - temp1 + temp2

                ! a = h, b = k, c = l+1
                ! ------------------------------------------------------------
                a = h; b = k; c = l+1

                dU_dU(a,b,c,h,k,l,ie) = dU_dU(a,b,c,h,k,l,ie)   &
                     + du_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffu1 &
                     + du_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffu2 + temp1

                dV_dV(a,b,c,h,k,l,ie) = dV_dV(a,b,c,h,k,l,ie)   &
                     + dv_eta_dot_dpdn(a,b,c,h,k,l+1) * half_rdpdn_diffv1 &
                     + dv_eta_dot_dpdn(a,b,c,h,k,l)   * half_rdpdn_diffv2 + temp1

                ! ===============================================================
                ! surface pressure derivatives
                ! ===============================================================

                ! a /= h, b = k (skips zeros in dvel_eta_dot_dpdn)
                ! ------------------------------------------------------------
                b = k
                do a = 1,np
                   if (a == h) cycle

                   dU_dPs(a,b,h,k,l,ie) = dU_dPs(a,b,h,k,l,ie)    &
                        + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffu1 &
                        + dps_eta_dot_dpdn(a,b,h,k,l)   * half_rdpdn_diffu2

                   dV_dPs(a,b,h,k,l,ie) = dV_dPs(a,b,h,k,l,ie)    &
                        + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffv1 &
                        + dps_eta_dot_dpdn(a,b,h,k,l)   * half_rdpdn_diffv2

                   dT_dPs(a,b,h,k,l,ie) = dT_dPs(a,b,h,k,l,ie)      &
                        + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffT1 &
                        + dps_eta_dot_dpdn(a,b,h,k,l)   * half_rdpdn_difft2
                end do

                ! a = h, b /= k (skips zeros in dvel_eta_dot_dpdn)
                ! ------------------------------------------------------------
                a = h
                do b = 1,np
                   if (b == k) cycle

                   dU_dPs(a,b,h,k,l,ie) = dU_dPs(a,b,h,k,l,ie)    &
                        + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffu1 &
                        + dps_eta_dot_dpdn(a,b,h,k,l)   * half_rdpdn_diffu2

                   dV_dPs(a,b,h,k,l,ie) = dV_dPs(a,b,h,k,l,ie)    &
                        + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffv1 &
                        + dps_eta_dot_dpdn(a,b,h,k,l)   * half_rdpdn_diffv2

                   dT_dPs(a,b,h,k,l,ie) = dT_dPs(a,b,h,k,l,ie)      &
                        + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffT1 &
                        + dps_eta_dot_dpdn(a,b,h,k,l)   * half_rdpdn_difft2
                end do

                ! a = h and b = k
                ! ------------------------------------------------------------
                a = h; b = k

                diffB = fptr%hvcoord%hybi(l+1) - fptr%hvcoord%hybi(l)

                temp1 = -0.5d0 * diffB * rdpdn(h,k,l)**2 * eta_dot_dpdn(h,k,l+1) ! eta_dot_dpdn l+1/2
                temp2 = -0.5d0 * diffB * rdpdn(h,k,l)**2 * eta_dot_dpdn(h,k,l)   ! eta_dot_dpdn l-1/2

                dU_dPs(a,b,h,k,l,ie) = dU_dPs(a,b,h,k,l,ie)    &
                     + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffu1 &
                     + dps_eta_dot_dpdn(a,b,h,k,l)   * half_rdpdn_diffu2 &
                     + temp1 * (fptr%base(ie)%state%v(h,k,1,l+1,np1) - fptr%base(ie)%state%v(h,k,1,l,np1)) &
                     + temp2 * (fptr%base(ie)%state%v(h,k,1,l,np1)   - fptr%base(ie)%state%v(h,k,1,l-1,np1))

                dV_dPs(a,b,h,k,l,ie) = dV_dPs(a,b,h,k,l,ie)    &
                     + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffv1 &
                     + dps_eta_dot_dpdn(a,b,h,k,l)   * half_rdpdn_diffv2 &
                     + temp1 * (fptr%base(ie)%state%v(h,k,2,l+1,np1) - fptr%base(ie)%state%v(h,k,2,l,np1)) &
                     + temp2 * (fptr%base(ie)%state%v(h,k,2,l,np1)   - fptr%base(ie)%state%v(h,k,2,l-1,np1))

                dT_dPs(a,b,h,k,l,ie) = dT_dPs(a,b,h,k,l,ie)      &
                     + dps_eta_dot_dpdn(a,b,h,k,l+1) * half_rdpdn_diffT1 &
                     + dps_eta_dot_dpdn(a,b,h,k,l)   * half_rdpdn_difft2 &
                     + temp1 * (fptr%base(ie)%state%T(h,k,l+1,np1) - fptr%base(ie)%state%T(h,k,l,np1)) &
                     + temp2 * (fptr%base(ie)%state%T(h,k,l,np1)   - fptr%base(ie)%state%T(h,k,l-1,np1))

                ! ==================================================================
                ! temperature derivatives
                ! ==================================================================

                temp1 = half_rdpdn(h,k,l) * eta_dot_dpdn(h,k,l+1) ! eta_dot_dpdn l+1/2
                temp2 = half_rdpdn(h,k,l) * eta_dot_dpdn(h,k,l)   ! eta_dot_dpdn l-1/2              

                ! a = h, b = k, c = l-1
                dT_dT(h,k,l-1,h,k,l,ie) = dT_dT(h,k,l-1,h,k,l,ie) - temp2

                ! a = h, b = k, c = l
                dT_dT(h,k,l,h,k,l,ie) = dT_dT(h,k,l,h,k,l,ie) - temp1 + temp2

                ! a = h, b = k, c = l+1
                dT_dT(h,k,l+1,h,k,l,ie) = dT_dT(h,k,l+1,h,k,l,ie) + temp1

             end do
          end do
       end do

       ! ------------------------------------------------------------------------
       ! Level nlev Jacobian Entries
       ! ------------------------------------------------------------------------
       l = nlev
       do k=1,np
          do h=1,np

             half_rdpdn_diffu2 = half_rdpdn(h,k,l) * (fptr%base(ie)%state%v(h,k,1,l,np1) - fptr%base(ie)%state%v(h,k,1,l-1,np1))
             half_rdpdn_diffv2 = half_rdpdn(h,k,l) * (fptr%base(ie)%state%v(h,k,2,l,np1) - fptr%base(ie)%state%v(h,k,2,l-1,np1))

             half_rdpdn_diffT2 = half_rdpdn(h,k,l) * (fptr%base(ie)%state%T(h,k,l,np1)   - fptr%base(ie)%state%T(h,k,l-1,np1))

             ! ==================================================================
             ! velocity derivatives
             ! ==================================================================

             do c = 1,nlev

                ! a /= h, b = k (skip zeros in dvel_eta_dot_dpdn)
                ! ------------------------------------------------------------
                b = k
                do a = 1,np
                   if (a == h) cycle

                   dU_dU(a,b,c,h,k,l,ie) = dU_dU(a,b,c,h,k,l,ie) + du_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffu2
                   dU_dV(a,b,c,h,k,l,ie) = dU_dV(a,b,c,h,k,l,ie) + dv_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffu2

                   dV_dU(a,b,c,h,k,l,ie) = dV_dU(a,b,c,h,k,l,ie) + du_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffv2
                   dV_dV(a,b,c,h,k,l,ie) = dV_dV(a,b,c,h,k,l,ie) + dv_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffv2

                   dT_dU(a,b,c,h,k,l,ie)  = dT_dU(a,b,c,h,k,l,ie)  + du_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffT2
                   dT_dV(a,b,c,h,k,l,ie)  = dT_dV(a,b,c,h,k,l,ie)  + dv_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffT2
                end do

                ! a = h, b /= k (skips zeros in dvel_eta_dot_dpdn)
                ! ------------------------------------------------------------
                a = h
                do b = 1,np
                   if (b == k) cycle

                   dU_dU(a,b,c,h,k,l,ie) = dU_dU(a,b,c,h,k,l,ie) + du_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffu2
                   dU_dV(a,b,c,h,k,l,ie) = dU_dV(a,b,c,h,k,l,ie) + dv_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffu2

                   dV_dU(a,b,c,h,k,l,ie) = dV_dU(a,b,c,h,k,l,ie) + du_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffv2
                   dV_dV(a,b,c,h,k,l,ie) = dV_dV(a,b,c,h,k,l,ie) + dv_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffv2

                   dT_dU(a,b,c,h,k,l,ie)  = dT_dU(a,b,c,h,k,l,ie)  + du_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffT2
                   dT_dV(a,b,c,h,k,l,ie)  = dT_dV(a,b,c,h,k,l,ie)  + dv_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffT2
                end do

                ! a = h, b = k, c = *
                ! ------------------------------------------------------------
                a = h; b = k

                dU_dV(a,b,c,h,k,l,ie) = dU_dV(a,b,c,h,k,l,ie) + dv_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffu2
                dV_dU(a,b,c,h,k,l,ie) = dV_dU(a,b,c,h,k,l,ie) + du_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffv2

                if ( c /= l-1 .and. c /= l) then               
                   dU_dU(a,b,c,h,k,l,ie) = dU_dU(a,b,c,h,k,l,ie) + du_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffu2
                   dV_dV(a,b,c,h,k,l,ie) = dV_dV(a,b,c,h,k,l,ie) + dv_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffv2
                end if

                dT_dU(a,b,c,h,k,l,ie)  = dT_dU(a,b,c,h,k,l,ie)  + du_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffT2
                dT_dV(a,b,c,h,k,l,ie)  = dT_dV(a,b,c,h,k,l,ie)  + dv_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffT2
             end do

             ! special cases
             ! ------------------------------------------------------------
             temp2 = half_rdpdn(h,k,l) * eta_dot_dpdn(h,k,l) ! eta_dot_dpdn l-1/2

             ! a = h, b = k, c = l-1
             ! ------------------------------------------------------------
             a = h; b = k; c = l-1

             dU_dU(h,k,l-1,h,k,l,ie) = dU_dU(h,k,l-1,h,k,l,ie) &
                  + du_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffu2 - temp2

             dV_dV(h,k,l-1,h,k,l,ie) = dV_dV(h,k,l-1,h,k,l,ie) &
                  + dv_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffv2 - temp2 

             ! a = h, b = k, c = l
             ! ------------------------------------------------------------
             a = h; b = k; c = l

             dU_dU(h,k,l,h,k,l,ie) = dU_dU(h,k,l,h,k,l,ie) &
                  + du_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffu2 + temp2

             dV_dV(h,k,l,h,k,l,ie) = dV_dV(h,k,l,h,k,l,ie) &
                  + dv_eta_dot_dpdn(a,b,c,h,k,l) * half_rdpdn_diffv2 + temp2

             ! ==================================================================
             ! surface pressure derivatives
             ! ==================================================================

             ! b = k to skip zeros in dvel_eta_dot_dpdn
             b = k
             do a = 1,np
                if (a == h) cycle

                dU_dPs(a,b,h,k,l,ie) = dU_dPs(a,b,h,k,l,ie) &
                     + dps_eta_dot_dpdn(a,b,h,k,l) * half_rdpdn_diffu2

                dV_dPs(a,b,h,k,l,ie) = dV_dPs(a,b,h,k,l,ie) &
                     + dps_eta_dot_dpdn(a,b,h,k,l) * half_rdpdn_diffv2

                dT_dPs(a,b,h,k,l,ie)  = dT_dPs(a,b,h,k,l,ie)  &
                     + dps_eta_dot_dpdn(a,b,h,k,l) * half_rdpdn_diffT2
             end do

             ! a = h to skip zeros in dvel_eta_dot_dpdn
             a = h
             do b = 1,np
                if (b == k) cycle

                dU_dPs(a,b,h,k,l,ie) = dU_dPs(a,b,h,k,l,ie) &
                     + dps_eta_dot_dpdn(a,b,h,k,l) * half_rdpdn_diffu2

                dV_dPs(a,b,h,k,l,ie) = dV_dPs(a,b,h,k,l,ie) &
                     + dps_eta_dot_dpdn(a,b,h,k,l) * half_rdpdn_diffv2

                dT_dPs(a,b,h,k,l,ie)  = dT_dPs(a,b,h,k,l,ie)  &
                     + dps_eta_dot_dpdn(a,b,h,k,l) * half_rdpdn_diffT2
             end do

             ! a = h, b = k
             a = h; b = k

             diffB = fptr%hvcoord%hybi(l+1) - fptr%hvcoord%hybi(l)

             temp2 = -0.5d0 * diffB * rdpdn(h,k,l)**2 * eta_dot_dpdn(h,k,l)   ! eta_dot_dpdn l-1/2

             dU_dPs(h,k,h,k,l,ie) = dU_dPs(h,k,h,k,l,ie)  &
                  + dps_eta_dot_dpdn(h,k,h,k,l) * half_rdpdn_diffu2 &
                  + temp2 * (fptr%base(ie)%state%v(h,k,1,l,np1) &
                  - fptr%base(ie)%state%v(h,k,1,l-1,np1))

             dV_dPs(h,k,h,k,l,ie) = dV_dPs(h,k,h,k,l,ie)  &
                  + dps_eta_dot_dpdn(h,k,h,k,l) * half_rdpdn_diffv2 &
                  + temp2 * (fptr%base(ie)%state%v(h,k,2,l,np1) &
                  - fptr%base(ie)%state%v(h,k,2,l-1,np1))          

             dT_dPs(h,k,h,k,l,ie) = dT_dPs(h,k,h,k,l,ie)    &
                  + dps_eta_dot_dpdn(h,k,h,k,l) * half_rdpdn_diffT2 &
                  + temp2 * (fptr%base(ie)%state%T(h,k,l,np1) &
                  - fptr%base(ie)%state%T(h,k,l-1,np1))

             ! ==================================================================
             ! temperature derivatives
             ! ==================================================================

             temp2 = half_rdpdn(h,k,l) * eta_dot_dpdn(h,k,l)   ! eta_dot_dpdn l-1/2

             ! a = h, b = k, c = l
             dT_dT(h,k,l,h,k,l,ie) = dT_dT(h,k,l,h,k,l,ie) + temp2

             ! a = h, b = k, c = l-1
             dT_dT(h,k,l-1,h,k,l,ie) = dT_dT(h,k,l-1,h,k,l,ie) - temp2

          end do
       end do

    end do

    call t_stopf('Jac_VertAdv')

  end subroutine Jac_VertAdv



  subroutine Jac_HorizAdv_Temp(fptr)
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the horizontal advection 
    ! terms in the temperature equation
    !----------------------------------------------------------------------------

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth
    use derivative_mod, only        : gradient_sphere

    implicit none

    !---------------------------------Arguments----------------------------------
    type(derived_type), pointer, intent(in) :: fptr
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    real(kind=real_kind) :: grad_T(np,np,2)
    real(kind=real_kind) :: temp1, temp2

    integer :: np1         ! time level
    integer :: a, b        ! derivative index
    integer :: h, k, l, ie ! equation index
    !----------------------------------------------------------------------------

    call t_startf('Jac_HorizAdv_Temp')

    np1 = fptr%tl%np1

    do ie = 1,nelemd

       do l = 1,nlev

          ! compute gradient of Temperature
          grad_T = gradient_sphere(fptr%base(ie)%state%T(:,:,l,np1), fptr%deriv, fptr%base(ie)%Dinv)

          do k = 1,np
             do h = 1,np

                ! ---------------------------------------------------------------
                ! velocity derivatives
                ! ---------------------------------------------------------------

                ! non-zero when a = h, b = k, c = l
                dT_dU(h,k,l,h,k,l,ie) = dT_dU(h,k,l,h,k,l,ie) + grad_T(h,k,1)
                dT_dV(h,k,l,h,k,l,ie) = dT_dV(h,k,l,h,k,l,ie) + grad_T(h,k,2)

                ! ---------------------------------------------------------------
                ! temperature derivatives
                ! ---------------------------------------------------------------

                temp1 = rrearth * (fptr%base(ie)%state%v(h,k,1,l,np1) * fptr%base(ie)%Dinv(h,k,1,1) &
                     + fptr%base(ie)%state%v(h,k,2,l,np1) * fptr%base(ie)%Dinv(h,k,1,2))

                temp2 = rrearth * (fptr%base(ie)%state%v(h,k,1,l,np1) * fptr%base(ie)%Dinv(h,k,2,1) &
                     + fptr%base(ie)%state%v(h,k,2,l,np1) * fptr%base(ie)%Dinv(h,k,2,2))

                ! a /= h, b = k, c = l
                do a = 1,np
                   if (a == h) cycle
                   dT_dT(a,k,l,h,k,l,ie) = dT_dT(a,k,l,h,k,l,ie) + temp1 * fptr%deriv%Dvv(a,h)
                end do

                ! a = h, b /= k, c = l
                do b = 1,np
                   if (b == k) cycle 
                   dT_dT(h,b,l,h,k,l,ie) = dT_dT(h,b,l,h,k,l,ie) + temp2 * fptr%deriv%Dvv(b,k)
                end do

                ! a = h, b = k, c = l (x and y-derivative terms)
                dT_dT(h,k,l,h,k,l,ie) = dT_dT(h,k,l,h,k,l,ie) &
                     + temp1 * fptr%deriv%Dvv(h,h) + temp2 * fptr%deriv%Dvv(k,k)

             end do
          end do
       end do

    end do

    call t_stopf('Jac_HorizAdv_Temp')

  end subroutine Jac_HorizAdv_Temp




  subroutine Jac_pressure_vert_vel(fptr)
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the pressure vertical
    ! velocity term in the temperature equation
    !----------------------------------------------------------------------------

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth, kappa
    use prim_si_mod, only           : preq_omega_ps
    use derivative_mod, only        : gradient_sphere, divergence_sphere

    implicit none

    !---------------------------------Arguments----------------------------------
    type(derived_type), pointer, intent(in) :: fptr
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    real(kind=real_kind) :: du_vel_grad_p(np,np,nlev,np,np,nlev) ! a,b,c,h,k,l
    real(kind=real_kind) :: dv_vel_grad_p(np,np,nlev,np,np,nlev) ! a,b,c,h,k,l
    real(kind=real_kind) :: dps_vel_grad_p(np,np,np,np,nlev)     ! a,b,  h,k,l

    real(kind=real_kind) :: du_div_dpdn_vel(np,np,nlev,np,np,nlev) ! a,b,c,h,k,l
    real(kind=real_kind) :: dv_div_dpdn_vel(np,np,nlev,np,np,nlev) ! a,b,c,h,k,l
    real(kind=real_kind) :: dps_div_dpdn_vel(np,np,np,np,nlev)     ! a,b,  h,k,l

    real(kind=real_kind) :: vtemp(np,np,2), grad_ps(np,np,2), grad_p(np,np,2,nlev)

    real(kind=real_kind) :: dp(np,np,nlev), p(np,np,nlev)
    real(kind=real_kind) :: vgrad_p(np,np,nlev), divdp(np,np,nlev)
    real(kind=real_kind) :: kappa_star(np,np,nlev), omega_p(np,np,nlev)
    real(kind=real_kind) :: ones(np,np,nlev) = 1.0d0

    real(kind=real_kind) :: T_v(np,np,nlev)

    real(kind=real_kind) :: v1, v2, temp1, temp2

    integer :: a,b,c,h,k,l,ie    ! loop counters
    integer :: np1, qn0
    !----------------------------------------------------------------------------

    call t_startf('Jac_pressure_vert_vel')

    np1  = fptr%tl%np1
    qn0 = fptr%n_Q

    do ie = 1,nelemd

       ! compute omega = v grad p  - integral_etatop^eta[ divdp ]
       ! ------------------------------------------------------------------------
       grad_ps = gradient_sphere(fptr%base(ie)%state%ps_v(:,:,np1), fptr%deriv, fptr%base(ie)%Dinv)

       do l = 1,nlev

          dp(:,:,l) = (fptr%hvcoord%hyai(l+1) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(l+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
               - (fptr%hvcoord%hyai(l) * fptr%hvcoord%ps0 + fptr%hvcoord%hybi(l) * fptr%base(ie)%state%ps_v(:,:,np1))

          p(:,:,l) = fptr%hvcoord%hyam(l) * fptr%hvcoord%ps0 + fptr%hvcoord%hybm(l) * fptr%base(ie)%state%ps_v(:,:,np1)

          grad_p(:,:,:,l) = fptr%hvcoord%hybm(l)*grad_ps(:,:,:)

          do k = 1,np
             do h = 1,np
                v1 = fptr%base(ie)%state%v(h,k,1,l,np1)
                v2 = fptr%base(ie)%state%v(h,k,2,l,np1)

                vgrad_p(h,k,l) = (v1*grad_p(h,k,1,l) + v2*grad_p(h,k,2,l))
                vtemp(h,k,1) = v1*dp(h,k,l)
                vtemp(h,k,2) = v2*dp(h,k,l)
             end do
          end do

          divdp(:,:,l) = divergence_sphere(vtemp, fptr%deriv, fptr%base(ie))

       end do

       ! note the ones array is used get just omega not omega/p
       call preq_omega_ps(omega_p, fptr%hvcoord, ones, vgrad_p, divdp)

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

       ! compute derivatives of v . grad(p) and div(dpdn v)
       ! ------------------------------------------------------------------------
       call Jac_vel_grad_p(du_vel_grad_p, dv_vel_grad_p, dps_vel_grad_p, fptr, ie)

       call Jac_div_dpdn_vel(du_div_dpdn_vel, dv_div_dpdn_vel, dps_div_dpdn_vel, fptr, ie)

       ! ------------------------------------------------------------------------
       ! Level 1 
       ! ------------------------------------------------------------------------
       l = 1
       do k = 1,np
          do h = 1,np

             ! ------------------------------------------------------------------
             ! velocity derivatives
             ! ------------------------------------------------------------------

             ! note the signs because pressure vertical velocity is subtracted
             ! in the resiudal equation

             temp1 = kappa_star(h,k,l) * T_v(h,k,l) / p(h,k,l)

             ! b = k, c = l
             do a = 1,np
                if (a == h) cycle

                dT_dU(a,k,l,h,k,l,ie) = dT_dU(a,k,l,h,k,l,ie) &
                     + 0.5d0 * temp1 * du_div_dpdn_vel(a,k,l,h,k,l)

                dT_dV(a,k,l,h,k,l,ie) = dT_dV(a,k,l,h,k,l,ie) &
                     + 0.5d0 * temp1 * dv_div_dpdn_vel(a,k,l,h,k,l)
             end do

             ! a = h, c = l
             do b = 1,np
                if (b == k) cycle

                dT_dU(h,b,l,h,k,l,ie) = dT_dU(h,b,l,h,k,l,ie) &
                     + 0.5d0 * temp1 * du_div_dpdn_vel(h,b,l,h,k,l)

                dT_dV(h,b,l,h,k,l,ie) = dT_dV(h,b,l,h,k,l,ie) &
                     + 0.5d0 * temp1 * dv_div_dpdn_vel(h,b,l,h,k,l)
             end do

             ! a = h, b = k, c = l
             dT_dU(h,k,l,h,k,l,ie) = dT_dU(h,k,l,h,k,l,ie) &
                  - temp1 * (du_vel_grad_p(h,k,l,h,k,l) &
                  - 0.5d0 * du_div_dpdn_vel(h,k,l,h,k,l))

             dT_dV(h,k,l,h,k,l,ie) = dT_dV(h,k,l,h,k,l,ie) &
                  - temp1 * (dv_vel_grad_p(h,k,l,h,k,l) &
                  - 0.5d0 * dv_div_dpdn_vel(h,k,l,h,k,l))

             ! ------------------------------------------------------------------
             ! surface pressure derivatives
             ! ------------------------------------------------------------------

             ! b = k
             do a = 1,np
                if (a == h) cycle
                dT_dPs(a,k,h,k,l,ie) = dT_dPs(a,k,h,k,l,ie) &
                     - temp1 * (dps_vel_grad_p(a,k,h,k,l) &
                     - 0.5d0 * dps_div_dpdn_vel(a,k,h,k,l))
             end do

             ! a = h
             do b = 1,np
                if (b == k) cycle
                dT_dPs(h,b,h,k,l,ie) = dT_dPs(h,b,h,k,l,ie) &
                     - temp1 * (dps_vel_grad_p(h,b,h,k,l) &
                     - 0.5d0 * dps_div_dpdn_vel(h,b,h,k,l))
             end do

             ! a = h, b = k
             dT_dPs(h,k,h,k,l,ie) = dT_dPs(h,k,h,k,l,ie) &
                  + temp1 * fptr%hvcoord%hybm(l) * omega_p(h,k,l) / p(h,k,l) &
                  - temp1 * (dps_vel_grad_p(h,k,h,k,l) &
                  - 0.5d0 * dps_div_dpdn_vel(h,k,h,k,l))

             ! ------------------------------------------------------------------
             ! temperature derivatives
             ! ------------------------------------------------------------------

             ! a = h, b = k, c = l
             dT_dT(h,k,l,h,k,l,ie) = dT_dT(h,k,l,h,k,l,ie) &
                  - kappa_star(h,k,l) * omega_p(h,k,l) / p(h,k,l)

          end do
       end do

       ! ------------------------------------------------------------------------
       ! Level 2 to nlev
       ! ------------------------------------------------------------------------
       do l = 2,nlev
          do k = 1,np
             do h = 1,np              

                ! ---------------------------------------------------------------
                ! velocity derivatives
                ! ---------------------------------------------------------------

                temp1 = kappa_star(h,k,l) * T_v(h,k,l) / p(h,k,l)

                ! case: c < l
                do c = 1,l-1

                   ! b = k, c < l
                   b = k
                   do a = 1,np
                      if (a == h) cycle

                      dT_dU(a,b,c,h,k,l,ie) = dT_dU(a,b,c,h,k,l,ie) &
                           + temp1 * du_div_dpdn_vel(a,b,c,h,k,c) ! note l = c

                      dT_dV(a,b,c,h,k,l,ie) = dT_dV(a,b,c,h,k,l,ie) &
                           + temp1 * dv_div_dpdn_vel(a,b,c,h,k,c) ! note l = c
                   end do

                   ! a = h, c < l
                   a = h
                   do b = 1,np
                      if (b == k) cycle

                      dT_dU(a,b,c,h,k,l,ie) = dT_dU(a,b,c,h,k,l,ie) &
                           + temp1 * du_div_dpdn_vel(a,b,c,h,k,c) ! note l = c

                      dT_dV(a,b,c,h,k,l,ie) = dT_dV(a,b,c,h,k,l,ie) &
                           + temp1 * dv_div_dpdn_vel(a,b,c,h,k,c) ! note l = c
                   end do

                   ! a = h, b = k, c < l
                   a = h; b = k; 

                   dT_dU(a,b,c,h,k,l,ie) = dT_dU(a,b,c,h,k,l,ie) &
                        + temp1 * du_div_dpdn_vel(a,b,c,h,k,c) ! note l = c

                   dT_dV(a,b,c,h,k,l,ie) = dT_dV(a,b,c,h,k,l,ie) &
                        + temp1 * dv_div_dpdn_vel(a,b,c,h,k,c) ! note l = c

                end do

                ! a /= h, b = k, c = l
                b = k; c = l
                do a = 1,np
                   if (a == h) cycle

                   dT_dU(a,b,c,h,k,l,ie) = dT_dU(a,b,c,h,k,l,ie) &
                        + 0.5d0 * temp1 * du_div_dpdn_vel(a,b,c,h,k,l)

                   dT_dV(a,b,c,h,k,l,ie) = dT_dV(a,b,c,h,k,l,ie) &
                        + 0.5d0 * temp1 * dv_div_dpdn_vel(a,b,c,h,k,l)
                end do

                ! a = h, b /= k, c = l
                a = h; c = l
                do b = 1,np
                   if (b == k) cycle

                   dT_dU(a,b,c,h,k,l,ie) = dT_dU(a,b,c,h,k,l,ie) &
                        + 0.5d0 * temp1 * du_div_dpdn_vel(a,b,c,h,k,l)

                   dT_dV(a,b,c,h,k,l,ie) = dT_dV(a,b,c,h,k,l,ie) &
                        + 0.5d0 * temp1 * dv_div_dpdn_vel(a,b,c,h,k,l)
                end do

                ! a = h, b = k, c = l
                a = h; b = k; c = l

                dT_dU(a,b,c,h,k,l,ie) = dT_dU(a,b,c,h,k,l,ie) &
                     - temp1 * (du_vel_grad_p(a,b,c,h,k,l) &
                     - 0.5d0 * du_div_dpdn_vel(a,b,c,h,k,l))

                dT_dV(a,b,c,h,k,l,ie) = dT_dV(a,b,c,h,k,l,ie) &
                     - temp1 * (dv_vel_grad_p(a,b,c,h,k,l) &
                     - 0.5d0 * dv_div_dpdn_vel(a,b,c,h,k,l))

                ! ---------------------------------------------------------------
                ! surface pressure derivatives
                ! ---------------------------------------------------------------

                ! a /= h, b = k
                b = k
                do a = 1,np
                   if (a == h) cycle

                   temp2 = 0.0d0
                   do c = 1,l-1
                      temp2 = temp2 + dps_div_dpdn_vel(a,b,h,k,c)
                   end do

                   dT_dPs(a,b,h,k,l,ie) = dT_dPs(a,b,h,k,l,ie) &
                        - temp1 * (dps_vel_grad_p(a,b,h,k,l) - temp2 &
                        - 0.5d0 * dps_div_dpdn_vel(a,b,h,k,l))
                end do

                ! a = h, b /= k
                a = h
                do b = 1,np
                   if (b == k) cycle

                   temp2 = 0.0d0
                   do c = 1,l-1
                      temp2 = temp2 + dps_div_dpdn_vel(a,b,h,k,c)
                   end do

                   dT_dPs(a,b,h,k,l,ie) = dT_dPs(a,b,h,k,l,ie) &
                        - temp1 * (dps_vel_grad_p(a,b,h,k,l) - temp2 &
                        - 0.5d0 * dps_div_dpdn_vel(a,b,h,k,l))
                end do

                ! a = h, b = k
                a = h; b = k

                temp2 = 0.0d0
                do c = 1,l-1
                   temp2 = temp2 + dps_div_dpdn_vel(a,b,h,k,c)
                end do

                dT_dPs(a,b,h,k,l,ie) = dT_dPs(a,b,h,k,l,ie) &
                     + temp1 * fptr%hvcoord%hybm(l) * omega_p(h,k,l) / p(h,k,l) &
                     - temp1 * (dps_vel_grad_p(a,b,h,k,l) - temp2 &
                     - 0.5d0 * dps_div_dpdn_vel(a,b,h,k,l))

                ! ---------------------------------------------------------------
                ! temperature derivatives
                ! ---------------------------------------------------------------

                ! a = h, b = k, c = l               
                dT_dT(h,k,l,h,k,l,ie) = dT_dT(h,k,l,h,k,l,ie) &
                     - kappa_star(h,k,l) * omega_p(h,k,l) / p(h,k,l)

             end do
          end do
       end do

    end do  

    call t_stopf('Jac_pressure_vert_vel')

  end subroutine Jac_pressure_vert_vel



  subroutine Jac_int_div_dpdn_vel(fptr)
    !----------------------------------------------------------------------------
    ! Computes contributions to the Jacobian matrix from the integral of the 
    ! divergence of dp/dn * velocity in the surface pressure equation
    !----------------------------------------------------------------------------

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use prim_derived_type_mod, only : derived_type

    implicit none

    !---------------------------------Arguments----------------------------------
    type(derived_type), pointer, intent(in) :: fptr
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------  
    real(kind=real_kind) :: du_div_dpdn_vel(np,np,nlev,np,np,nlev) ! a,b,c,h,k,l
    real(kind=real_kind) :: dv_div_dpdn_vel(np,np,nlev,np,np,nlev) ! a,b,c,h,k,l
    real(kind=real_kind) :: dps_div_dpdn_vel(np,np,np,np,nlev)     ! a,b,  h,k,l

    integer :: a, b, c     ! derivative index
    integer :: i, j, k, ie ! equation index
    integer :: np1
    !----------------------------------------------------------------------------

    call t_startf('Jac_int_div_dpdn_vel') 

    np1 = fptr%tl%np1

    do ie = 1,nelemd

       call Jac_div_dpdn_vel(du_div_dpdn_vel, dv_div_dpdn_vel, dps_div_dpdn_vel, fptr, ie)

       ! ---------------------------------------------------------------
       ! velocity derivatives
       ! ---------------------------------------------------------------
       do j = 1,np
          do i = 1,np              

             do c = 1,nlev

                ! a /= i, b = j, c = k
                ! ---------------------------------------------------------------
                do a = 1,np
                   if (a == i) cycle

                   dPs_dU(a,j,c,i,j,ie) = dPs_dU(a,j,c,i,j,ie) + &
                        du_div_dpdn_vel(a,j,c,i,j,c)

                   dPs_dV(a,j,c,i,j,ie) = dPs_dV(a,j,c,i,j,ie) + &
                        dv_div_dpdn_vel(a,j,c,i,j,c)
                end do

                ! a = i, b /= j, c = k
                ! ---------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle

                   dPs_dU(i,b,c,i,j,ie) = dPs_dU(i,b,c,i,j,ie) + &
                        du_div_dpdn_vel(i,b,c,i,j,c)

                   dPs_dV(i,b,c,i,j,ie) = dPs_dV(i,b,c,i,j,ie) + &
                        dv_div_dpdn_vel(i,b,c,i,j,c)
                end do

                ! a = i, b = j, c = k
                ! ---------------------------------------------------------------
                dPs_dU(i,j,c,i,j,ie) = dPs_dU(i,j,c,i,j,ie) + &
                     du_div_dpdn_vel(i,j,c,i,j,c)

                dPs_dV(i,j,c,i,j,ie) = dPs_dV(i,j,c,i,j,ie) + &
                     dv_div_dpdn_vel(i,j,c,i,j,c)

             end do
          end do
       end do


       ! ------------------------------------------------------------------
       ! surface pressure derivatives
       ! ------------------------------------------------------------------
       do k = 1,nlev ! vertical integral 
          do j = 1,np
             do i = 1,np              
                
                ! a /= i, b = j
                ! ---------------------------------------------------------------
                do a = 1,np
                   if (a == i) cycle                                      
                   dPs_dPs(a,j,i,j,ie) = dPs_dPs(a,j,i,j,ie) + dps_div_dpdn_vel(a,j,i,j,k)
                end do
                                
                ! a = i, b /= j
                ! ---------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle                  
                   dPs_dPs(i,b,i,j,ie) = dPs_dPs(i,b,i,j,ie) + dps_div_dpdn_vel(i,b,i,j,k) 
                end do
                
                ! a = i, b = j
                ! ---------------------------------------------------------------
                dPs_dPs(i,j,i,j,ie) = dPs_dPs(i,j,i,j,ie) + dps_div_dpdn_vel(i,j,i,j,k) 
                
             end do
          end do   
       end do
              
    end do

    call t_stopf('Jac_int_div_dpdn_vel') 

  end subroutine Jac_int_div_dpdn_vel
  
  
  
  ! =============================================================================
  ! Utility subroutines necessary for computing Jacobian entries
  ! =============================================================================
  
  
  
  subroutine Jac_vel_grad_p(du_vel_grad_p, dv_vel_grad_p, dps_vel_grad_p, fptr, ie)
    !----------------------------------------------------------------------------
    ! compute Jacobian of (v dot grad p), used in pressure vert velocity
    !----------------------------------------------------------------------------
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth
    use derivative_mod, only        : gradient_sphere
    
    implicit none

    !---------------------------------Arguments----------------------------------
    real(kind=real_kind), intent(out) :: du_vel_grad_p(np,np,nlev,np,np,nlev) ! a,b,c,h,k,l
    real(kind=real_kind), intent(out) :: dv_vel_grad_p(np,np,nlev,np,np,nlev) ! a,b,c,h,k,l
    real(kind=real_kind), intent(out) :: dps_vel_grad_p(np,np,np,np,nlev)     ! a,b,  h,k,l

    type(derived_type), pointer, intent(in) :: fptr

    integer, intent(in) :: ie
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    ! real(kind=real_kind), dimension(np,np,nlev)   :: p
    real(kind=real_kind), dimension(np,np,2,nlev) :: grad_p
    real(kind=real_kind), dimension(np,np,2)      :: grad_ps

    real(kind=real_kind) :: temp1, temp2

    integer :: a,b,c ! derivative index
    integer :: h,k,l ! equation index
    integer :: np1   ! current iteration
    !----------------------------------------------------------------------------

    call t_startf('Jac_vel_grad_p')

    du_vel_grad_p  = 0.0d0
    dv_vel_grad_p  = 0.0d0
    dps_vel_grad_p = 0.0d0
    
    np1 = fptr%tl%np1

    grad_ps = gradient_sphere(fptr%base(ie)%state%ps_v(:,:,np1), fptr%deriv, fptr%base(ie)%Dinv)

    do l = 1,nlev

       ! p(:,:,l) = fptr%hvcoord%hyam(l) * fptr%hvcoord%ps0 &
       !      + fptr%hvcoord%hybm(l) * fptr%base(ie)%state%ps_v(:,:,np1)

       grad_p(:,:,:,l) = fptr%hvcoord%hybm(l)*grad_ps(:,:,:)
    end do
       
    do l = 1,nlev
       do k = 1,np
          do h = 1,np
                 
             ! ==================================================================
             ! velocity derivatives
             ! ==================================================================

             ! a = h, b = k, c = l
             a = h; b = k; c = l
             
             du_vel_grad_p(a,b,c,h,k,l) = grad_p(h,k,1,l)
             dv_vel_grad_p(a,b,c,h,k,l) = grad_p(h,k,2,l)
             
             ! ==================================================================
             ! surface pressure derivatives
             ! ==================================================================

             temp1 = rrearth * fptr%hvcoord%hybm(l) &
                  * (fptr%base(ie)%state%v(h,k,1,l,np1) * fptr%base(ie)%Dinv(h,k,1,1) &
                  + fptr%base(ie)%state%v(h,k,2,l,np1) * fptr%base(ie)%Dinv(h,k,1,2))

             temp2 = rrearth * fptr%hvcoord%hybm(l) &
                  * (fptr%base(ie)%state%v(h,k,1,l,np1) * fptr%base(ie)%Dinv(h,k,2,1) &
                  + fptr%base(ie)%state%v(h,k,2,l,np1) * fptr%base(ie)%Dinv(h,k,2,2))

             ! a /= h, b = k (x-derivative)
             b = k 
             do a = 1,np
                if (a == h) cycle
                dps_vel_grad_p(a,b,h,k,l) = temp1 * fptr%deriv%Dvv(a,h) 
             end do
             
             ! a = h, b /= k (y-derivative)
             a = h
             do b = 1,np
                if (b == k) cycle
                dps_vel_grad_p(a,b,h,k,l) = temp2 * fptr%deriv%Dvv(b,k)
             end do

             ! a = h, b = k (x and y-derivative)
             a = h; b = k
             dps_vel_grad_p(a,b,h,k,l) = temp1 * fptr%deriv%Dvv(a,h) + temp2 * fptr%deriv%Dvv(b,k)

          end do
       end do
    end do
          
    call t_stopf('Jac_vel_grad_p')
    
  end subroutine Jac_vel_grad_p
  
  
  
  subroutine Jac_div_dpdn_vel(du_vals, dv_vals, dps_vals, fptr, ie)
    !----------------------------------------------------------------------------
    ! compute Jacobian of divergence(dpdn * vel)
    !----------------------------------------------------------------------------

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar
    use prim_derived_type_mod, only : derived_type
    use physical_constants, only    : rrearth
    
    implicit none

    !---------------------------------Arguments----------------------------------
    real(kind=real_kind), intent(out) :: du_vals(np,np,nlev,np,np,nlev) ! a,b,c,h,k,l
    real(kind=real_kind), intent(out) :: dv_vals(np,np,nlev,np,np,nlev) ! a,b,c,h,k,l
    real(kind=real_kind), intent(out) :: dps_vals(np,np,np,np,nlev)     ! a,b,  h,k,l
    type(derived_type), pointer, intent(in) :: fptr

    integer, intent(in) :: ie 
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    real(kind=real_kind), dimension(np,np,nlev) :: dpdn
    real(kind=real_kind), dimension(nlev)       :: diffBi

    real(kind=real_kind) :: temp1, temp2, temp3

    integer :: a,b   ! derivative index
    integer :: i,j,k ! equation index
    integer :: np1   ! time level
    !----------------------------------------------------------------------------

    call t_startf('Jac_div_dpdn_vel')

    ! get time level
    np1 = fptr%tl%np1

    ! initialize outputs
    du_vals  = 0.0d0
    dv_vals  = 0.0d0
    dps_vals = 0.0d0
    
    ! compute dp/deta and d(dp/deta)/dps
    do k = 1,nlev
       dpdn(:,:,k) = (fptr%hvcoord%hyai(k+1) * fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k+1) * fptr%base(ie)%state%ps_v(:,:,np1)) &
            - (fptr%hvcoord%hyai(k)*fptr%hvcoord%ps0 &
            + fptr%hvcoord%hybi(k)*fptr%base(ie)%state%ps_v(:,:,np1))

       diffBi(k) = fptr%hvcoord%hybi(k+1) - fptr%hvcoord%hybi(k)
    end do

    do k = 1,nlev       
       do j = 1,np
          do i = 1,np

             temp1 = rrearth * fptr%base(ie)%rmetdet(i,j)
             
             ! a /= i, b = j, c = k (x-derivative)
             ! ------------------------------------------------------------------
             do a = 1,np
                if (a == i) cycle
                
                temp2 = temp1 * fptr%base(ie)%metdet(a,j) * dpdn(a,j,k) * fptr%deriv%Dvv(a,i)
                
                du_vals(a,j,k,i,j,k) = temp2 * fptr%base(ie)%Dinv(a,j,1,1)

                dv_vals(a,j,k,i,j,k) = temp2 * fptr%base(ie)%Dinv(a,j,1,2)

                dps_vals(a,j,i,j,k)  = temp1 * fptr%base(ie)%metdet(a,j) * diffBi(k) * fptr%deriv%Dvv(a,i) &
                     * (fptr%base(ie)%Dinv(a,j,1,1) * fptr%base(ie)%state%v(a,j,1,k,np1) &
                     + fptr%base(ie)%Dinv(a,j,1,2) * fptr%base(ie)%state%v(a,j,2,k,np1))
             end do
             
             ! a = i, b /= j, c = k (y-derivative)
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle

                temp2 = temp1 * fptr%base(ie)%metdet(i,b) * dpdn(i,b,k) * fptr%deriv%Dvv(b,j)
                
                du_vals(i,b,k,i,j,k) = temp2 * fptr%base(ie)%Dinv(i,b,2,1)

                dv_vals(i,b,k,i,j,k) = temp2 * fptr%base(ie)%Dinv(i,b,2,2)

                dps_vals(i,b,i,j,k)  = temp1 * fptr%base(ie)%metdet(i,b) * diffBi(k) * fptr%deriv%Dvv(b,j) &
                     * (fptr%base(ie)%Dinv(i,b,2,1) * fptr%base(ie)%state%v(i,b,1,k,np1) &
                     + fptr%base(ie)%Dinv(i,b,2,2) * fptr%base(ie)%state%v(i,b,2,k,np1))
             end do

             ! a = h, b = k, c = l (x and y-derivative)
             ! ------------------------------------------------------------------                            
             temp2 = rrearth * dpdn(i,j,k) * fptr%deriv%Dvv(i,i)

             temp3 = rrearth * dpdn(i,j,k) * fptr%deriv%Dvv(j,j)
             
             du_vals(i,j,k,i,j,k) = temp2 * fptr%base(ie)%Dinv(i,j,1,1) &
                  + temp3 * fptr%base(ie)%Dinv(i,j,2,1)
             
             dv_vals(i,j,k,i,j,k) = temp2 * fptr%base(ie)%Dinv(i,j,1,2) &
                  + temp3 * fptr%base(ie)%Dinv(i,j,2,2)

             dps_vals(i,j,i,j,k) = rrearth * diffBi(k) * fptr%deriv%Dvv(i,i) &
                  * (fptr%base(ie)%Dinv(i,j,1,1) * fptr%base(ie)%state%v(i,j,1,k,np1) &
                  + fptr%base(ie)%Dinv(i,j,1,2) * fptr%base(ie)%state%v(i,j,2,k,np1)) &
                  + rrearth * diffBi(k) * fptr%deriv%Dvv(j,j) &
                  * (fptr%base(ie)%Dinv(i,j,2,1) * fptr%base(ie)%state%v(i,j,1,k,np1) &
                  + fptr%base(ie)%Dinv(i,j,2,2) * fptr%base(ie)%state%v(i,j,2,k,np1))
          end do
       end do
    end do
   
    call t_stopf('Jac_div_dpdn_vel')
    
  end subroutine Jac_div_dpdn_vel



  subroutine Jac_etadot_dpdn(du_vals, dv_vals, dps_vals, fptr, ie)
    !----------------------------------------------------------------------------
    ! compute Jacobian of (eta_dot * dpdn) REVISED VERSION
    !----------------------------------------------------------------------------
    
    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar
    use prim_derived_type_mod, only : derived_type
   
    implicit none

    !---------------------------------Arguments----------------------------------
    real(kind=real_kind), intent(out) :: du_vals(np,np,nlev,np,np,nlev+1) ! a,b,c,h,k,l
    real(kind=real_kind), intent(out) :: dv_vals(np,np,nlev,np,np,nlev+1) ! a,b,c,h,k,l
    real(kind=real_kind), intent(out) :: dps_vals(np,np,np,np,nlev+1)     ! a,b,  h,k,l
    type(derived_type), pointer, intent(in) :: fptr 
    integer, intent(in) :: ie
    !----------------------------------------------------------------------------

    !------------------------------Local workspace-------------------------------
    real(kind=real_kind) :: du_div_dpdn_vel(np,np,nlev,np,np,nlev) ! a,b,c,h,k,l
    real(kind=real_kind) :: dv_div_dpdn_vel(np,np,nlev,np,np,nlev) ! a,b,c,h,k,l
    real(kind=real_kind) :: dps_div_dpdn_vel(np,np,np,np,nlev)     ! a,b,  h,k,l

    real(kind=real_kind) :: Bi, Bim1

    integer :: a,b,c
    integer :: i,j,k
    !----------------------------------------------------------------------------

    call t_startf('Jac_etadot_dpdn')
    
    du_vals  = 0.0d0
    dv_vals  = 0.0d0
    dps_vals = 0.0d0

    call Jac_div_dpdn_vel(du_div_dpdn_vel, dv_div_dpdn_vel, dps_div_dpdn_vel, fptr, ie)

    ! ===========================================================================
    ! velocity derivatives
    ! ===========================================================================
    do k = 1,nlev-1
       
       Bim1 = fptr%hvcoord%hybi(k+1) - 1.0d0
       Bi   = fptr%hvcoord%hybi(k+1)
       
       do j = 1,np
          do i = 1,np
                 
             ! c <= k (contribution from both summations in etadot_dpdn)
             ! ------------------------------------------------------------------
             do c = 1,k

                ! a /= i, b = j (skips zeros in div_dpdn_vel)
                ! ---------------------------------------------------------------
                do a = 1,np
                   du_vals(a,j,c,i,j,k+1) = Bim1 * du_div_dpdn_vel(a,j,c,i,j,c)
                   dv_vals(a,j,c,i,j,k+1) = Bim1 * dv_div_dpdn_vel(a,j,c,i,j,c)
                end do

                ! a = i, b /= j (skips zeros in div_dpdn_vel)
                ! ---------------------------------------------------------------
                do b = 1,np
                   du_vals(i,b,c,i,j,k+1) = Bim1 * du_div_dpdn_vel(i,b,c,i,j,c)
                   dv_vals(i,b,c,i,j,k+1) = Bim1 * dv_div_dpdn_vel(i,b,c,i,j,c)
                end do
             end do

             ! c > k (contribution from first summation in etadot_dpdn)
             ! ------------------------------------------------------------------
             do c = k+1,nlev

                ! a /= i, b = j (skips zeros in div_dpdn_vel)
                ! ---------------------------------------------------------------
                do a = 1,np
                   du_vals(a,j,c,i,j,k+1) = Bi * du_div_dpdn_vel(a,j,c,i,j,c)
                   dv_vals(a,j,c,i,j,k+1) = Bi * dv_div_dpdn_vel(a,j,c,i,j,c)
                end do

                ! a = i, b /= j (skips zeros in div_dpdn_vel)
                ! ---------------------------------------------------------------
                do b = 1,np
                   du_vals(i,b,c,i,j,k+1) = Bi * du_div_dpdn_vel(i,b,c,i,j,c)
                   dv_vals(i,b,c,i,j,k+1) = Bi * dv_div_dpdn_vel(i,b,c,i,j,c)

                end do
             end do

          end do
       end do
    end do

    ! ===========================================================================
    ! surface pressure derivatives
    ! ===========================================================================

    do k = 1,nlev-1
       
       Bi = fptr%hvcoord%hybi(k+1)
       
       do j = 1,np
          do i = 1,np

             ! ==================================================================
             ! integral of div_dpdn_vel from 1 to nlev
             ! ==================================================================
             do c = 1,nlev

                ! a /= i, b = j
                ! ---------------------------------------------------------------
                do a = 1,np                 
                   if (a == i) cycle
                   dps_vals(a,j,i,j,k+1) = dps_vals(a,j,i,j,k+1) &
                        + dps_div_dpdn_vel(a,j,i,j,c)
                end do

                ! a = i, b /= j
                ! ---------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle
                   dps_vals(i,b,i,j,k+1) = dps_vals(i,b,i,j,k+1) &
                        + dps_div_dpdn_vel(i,b,i,j,c)
                end do

                ! a = i, b = j
                ! ---------------------------------------------------------------
                dps_vals(i,j,i,j,k+1) = dps_vals(i,j,i,j,k+1) &
                     + dps_div_dpdn_vel(i,j,i,j,c)
                
             end do

             ! ==================================================================
             ! multiply by Bi
             ! ==================================================================

             ! a /= i, b = j
             ! ------------------------------------------------------------------
             do a = 1,np                 
                if (a == i) cycle
                dps_vals(a,j,i,j,k+1) = Bi * dps_vals(a,j,i,j,k+1)
             end do
             
             ! a = i, b /= j
             ! ------------------------------------------------------------------
             do b = 1,np
                if (b == j) cycle
                dps_vals(i,b,i,j,k+1) = Bi * dps_vals(i,b,i,j,k+1)
             end do            

             ! a = i, b = j
             ! ------------------------------------------------------------------
             dps_vals(i,j,i,j,k+1) = Bi * dps_vals(i,j,i,j,k+1)

             ! ==================================================================
             ! integral of div_dpdn_vel from 1 to l
             ! ==================================================================
             do c = 1,k

                ! a /= i, b = j 
                ! ---------------------------------------------------------------
                do a = 1,np                 
                   if (a == i) cycle
                   dps_vals(a,j,i,j,k+1) = dps_vals(a,j,i,j,k+1) &
                        - dps_div_dpdn_vel(a,j,i,j,c)
                end do
             
                ! a = i, b /= j
                ! ---------------------------------------------------------------
                do b = 1,np
                   if (b == j) cycle
                   dps_vals(i,b,i,j,k+1) = dps_vals(i,b,i,j,k+1) &
                        - dps_div_dpdn_vel(i,b,i,j,c)
                end do

                ! a = h, b = k
                ! ---------------------------------------------------------------
                dps_vals(i,j,i,j,k+1) = dps_vals(i,j,i,j,k+1) &
                     - dps_div_dpdn_vel(i,j,i,j,c)

             end do
             
          end do
       end do
    end do
   
    call t_stopf('Jac_etadot_dpdn')
    
  end subroutine Jac_etadot_dpdn


  ! =============================================================================
  ! Subroutines to compute preconditioner operations using Jacobian blocks
  ! =============================================================================
    
  
  subroutine Ablock_matvec(xvec, nxvec, yvec, c_ptr_to_object) &
       bind(C,name='Ablock_matvec')
    ! ---------------------------------------------------------------------------
    ! Operator for A = dv/dv block of the Jacobian matrix
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use prim_derived_type_mod, only : derived_type
    use edge_mod, only              : edgevpack, edgevunpack
    use bndry_mod, only             : bndry_exchangev
    
    implicit none

    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xvec(nxvec)    
    integer(c_int), intent(in), value :: nxvec
    real(c_double), intent(out)       :: yvec(nxvec)    
    type(c_ptr)                       :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr => NULL()

    real(kind=real_kind), dimension(np,np,nlev)          :: vel1, vel2
    real(kind=real_kind), dimension(np,np,nlev)          :: utend1, utend2
    real(kind=real_kind), dimension(np,np,nlev)          :: vtend1, vtend2
    real(kind=real_kind), dimension(np,np,2,nlev, nelemd) :: vtens
           
    integer :: nets, nete
    integer :: np1
    integer :: kptr
    integer :: i,j,k,ie,lx
    ! ---------------------------------------------------------------------------
    
    call t_startf('Ablock_matvec')
    
    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*)  ">>> HOMME: Ablock_matvec"
#endif
    
    np1   = fptr%tl%np1   
    nets  = fptr%nets     
    nete  = fptr%nete     

    ! unpack 1D input array
    lx = 0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,1,k,np1) = xvec(lx)
             end do
          end do
       end do
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,2,k,np1) = xvec(lx)
             end do
          end do
       end do
    end do

    
    do ie = nets,nete

       vel1 = fptr%base(ie)%state%v(:,:,1,:,np1)
       vel2 = fptr%base(ie)%state%v(:,:,2,:,np1)
       
       ! compute MatVecProd for vertical advection
       call MatVecProd3Dx3D(utend1, dU_dU, vel1, ie)
       call MatVecProd3Dx3D(utend2, dU_dV, vel2, ie)

       call MatVecProd3Dx3D(vtend1, dV_dU, vel1, ie)
       call MatVecProd3Dx3D(vtend2, dV_dV, vel2, ie)
       
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                ! u vel 
                vtens(i,j,1,k,ie) = fptr%base(ie)%spheremp(i,j)*(utend1(i,j,k) + utend2(i,j,k))
                
                ! v vel 
                vtens(i,j,2,k,ie) = fptr%base(ie)%spheremp(i,j)*(vtend1(i,j,k) + vtend2(i,j,k))   
             end do
          end do
       end do
       
       kptr = 0
       call edgeVpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
    end do !ie

    call bndry_exchangeV(fptr%hybrid,edge3p1)

    kptr = 0
    do ie = nets,nete
       call edgeVunpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
    end do

    ! pack 1D output array
    lx = 0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                yvec(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
             end do
          end do
       end do
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                yvec(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
             end do
          end do
       end do
    end do

    call t_stopf('Ablock_matvec')

#ifdef DEBUG_PRINT_ON
    if (fptr%hybrid%masterthread) write(*,*)  "<<< HOMME: Ablock_matvec"
#endif
   
  end subroutine Ablock_matvec


  
  subroutine Bblock_matvec(xvec, nxvec, yvec, nyvec, c_ptr_to_object) &
       bind(C,name='Bblock_matvec')
    ! ---------------------------------------------------------------------------
    ! Operator for B = dv/dps block of the Jacobian matrix
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use prim_derived_type_mod, only : derived_type
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use edge_mod, only              : edgevpack, edgevunpack
    use bndry_mod, only             : bndry_exchangev
    use control_mod, only           : tstep_type
    
    implicit none

    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xvec(nxvec)
    integer(c_int), intent(in), value :: nxvec
    real(c_double), intent(out)       :: yvec(nyvec)
    integer(c_int), intent(in), value :: nyvec
    type(c_ptr)                       :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()
    
    real(kind=real_kind), dimension(np,np)              :: ps
    real(kind=real_kind), dimension(np,np,nlev)         :: MatVecProd1, MatVecProd2
    real(kind=real_kind), dimension(np,np,2,nlev, nelemd) :: vtens
    
    integer :: nets, nete
    integer :: np1
    integer :: kptr
    integer :: i,j,k,ie,lx
    ! ---------------------------------------------------------------------------

    call t_startf('Bblock_matvec')
    
    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*)  ">>> HOMME: Bblock_matvec"
#endif
    
    np1  = fptr%tl%np1   
    nets = fptr%nets     
    nete = fptr%nete     
    
    ! unpack 1D input array
    lx = 0
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             lx = lx+1
             fptr%base(ie)%state%ps_v(i,j,np1) = xvec(lx)
          end do
       end do
    end do

    do ie = nets,nete

       ps = fptr%base(ie)%state%ps_v(:,:,np1)

       call MatVecProd3Dx2D(MatVecProd1, dU_dPs, ps, ie)
       call MatVecProd3Dx2D(MatVecProd2, dV_dPs, ps, ie)

       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                ! u vel
                vtens(i,j,1,k,ie) = fptr%base(ie)%spheremp(i,j) * MatVecProd1(i,j,k)

                ! v vel
                vtens(i,j,2,k,ie) = fptr%base(ie)%spheremp(i,j) * MatVecProd2(i,j,k)
             end do
          end do
       end do

       kptr = 0
       call edgeVpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
    end do

    call bndry_exchangeV(fptr%hybrid, edge3p1)

    kptr = 0
    do ie = nets,nete
       call edgeVunpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
    end do

    ! pack 1D output array
    lx = 0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                yvec(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
             end do
          end do
       end do
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                yvec(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
             end do
          end do
       end do
    end do

    call t_stopf('Bblock_matvec')

#ifdef DEBUG_PRINT_ON
    if (fptr%hybrid%masterthread) write(*,*)  "<<< HOMME: Bblock_matvec"
#endif
    
  end subroutine Bblock_matvec


  
  subroutine Cblock_matvec(xvec, nxvec, yvec, nyvec, c_ptr_to_object) &
       bind(C,name='Cblock_matvec')
    ! ---------------------------------------------------------------------------
    ! Operator for C = dv/dT block of the Jacobian matrix
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use prim_derived_type_mod, only : derived_type
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use edge_mod, only              : edgevpack, edgevunpack
    use bndry_mod, only             : bndry_exchangev
    use control_mod, only           : tstep_type

    implicit none

    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xvec(nxvec)
    integer(c_int), intent(in), value :: nxvec
    real(c_double), intent(out)       :: yvec(nyvec)
    integer(c_int), intent(in), value :: nyvec
    type(c_ptr)                       :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,nlev) :: T
    real(kind=real_kind), dimension(np,np,nlev) :: MatVecProd1, MatVecProd2
    real(kind=real_kind), dimension(np,np,2,nlev, nelemd) :: vtens
    
    integer :: nets, nete
    integer :: np1
    integer :: kptr
    integer :: i,j,k,ie,lx    
    ! ---------------------------------------------------------------------------

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*)  ">>> HOMME: Cblock_matvec"
#endif

    call t_startf('Cblock_matvec')
    
    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr
    
    np1  = fptr%tl%np1   
    nets = fptr%nets     
    nete = fptr%nete     
    
    ! unpack 1D input array
    lx = 0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                fptr%base(ie)%state%T(i,j,k,np1) = xvec(lx)
             end do
          end do
       end do
    end do

    do ie = nets,nete

       T = fptr%base(ie)%state%T(:,:,:,np1)

       call MatVecProd3Dx3D(MatVecProd1, dU_dT, T, ie)
       call MatVecProd3Dx3D(MatVecProd2, dV_dT, T, ie)

       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                ! u vel
                vtens(i,j,1,k,ie) = fptr%base(ie)%spheremp(i,j) * MatVecProd1(i,j,k)

                ! v vel
                vtens(i,j,2,k,ie) = fptr%base(ie)%spheremp(i,j) * MatVecProd2(i,j,k)
             end do
          end do
       end do

       kptr = 0
       call edgeVpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
    end do

    call bndry_exchangeV(fptr%hybrid, edge3p1)

    kptr = 0
    do ie = nets,nete
       call edgeVunpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
    end do

    ! pack 1D output array
    lx = 0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                yvec(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
             end do
          end do
       end do
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                yvec(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
             end do
          end do
       end do
    end do

    call t_stopf('Cblock_matvec')

#ifdef DEBUG_PRINT_ON
    if (fptr%hybrid%masterthread) write(*,*)  "<<< HOMME: Cblock_matvec"
#endif  

  end subroutine Cblock_matvec


  
  subroutine Dblock_matvec(xvec, nxvec, yvec, nyvec, c_ptr_to_object) &
       bind(C,name='Dblock_matvec')
    ! ---------------------------------------------------------------------------
    ! Operator for D = dps/dv block of the Jacobian matrix
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding 

    use kinds, only                 : real_kind
    use prim_derived_type_mod, only : derived_type
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use edge_mod, only              : edgevpack, edgevunpack
    use bndry_mod, only             : bndry_exchangev
    use control_mod, only           : tstep_type

    implicit none

    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xvec(nxvec)
    integer(c_int), intent(in), value :: nxvec
    real(c_double), intent(out)       :: yvec(nyvec)
    integer(c_int), intent(in), value :: nyvec
    type(c_ptr)                       :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np)       :: MatVecProd1, MatVecProd2   
    real(kind=real_kind), dimension(np,np,nlev)  :: vel1, vel2
    real(kind=real_kind), dimension(np,np, nelemd) :: ptens
    
    integer :: nets, nete
    integer :: np1                
    integer :: kptr               
    integer :: i,j,k,ie,lx 
    ! ---------------------------------------------------------------------------
    
    call t_startf('Dblock_matvec')
    
    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*)  ">>> HOMME: Dblock_matvec"
#endif
    
    np1   = fptr%tl%np1   
    nets  = fptr%nets     
    nete  = fptr%nete     
    
    ! unpack 1D input array
    lx = 0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,1,k,np1) = xvec(lx)
             end do
          end do
       end do
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,2,k,np1) = xvec(lx)
             end do
          end do
       end do
    end do

    do ie = nets,nete
       
       vel1 = fptr%base(ie)%state%v(:,:,1,:,np1)
       vel2 = fptr%base(ie)%state%v(:,:,2,:,np1)
       
       call MatVecProd2Dx3D(MatVecProd1, dPs_dU, vel1, ie)
       call MatVecProd2Dx3D(MatVecProd2, dPs_dV, vel2, ie)
       
       do j = 1,np
          do i = 1,np
             ptens(i,j,ie) = fptr%base(ie)%spheremp(i,j)*(MatVecProd1(i,j) + MatVecProd2(i,j))
          end do
       end do

       ! pack tendencies into exchange buffer for DSS
       kptr = 0
       call edgeVpack(edge3p1, ptens(:,:,ie), 1, kptr, ie)
    end do

    call bndry_exchangeV(fptr%hybrid, edge3p1)

    ! unpack tendencies from exchange buffer for DSS
    kptr = 0
    do ie = nets,nete
       call edgeVunpack(edge3p1, ptens(:,:,ie), 1, kptr, ie)
    end do !ie

    ! pack 1D output array
    lx = 0
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             lx = lx+1
             yvec(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,ie) 
          end do
       end do
    end do
    
    call t_stopf('Dblock_matvec')

#ifdef DEBUG_PRINT_ON
    if (fptr%hybrid%masterthread) write(*,*)  "<<< HOMME: Dblock_matvec"
#endif

  end subroutine Dblock_matvec


  
  subroutine Eblock_matvec(xvec, nxvec, yvec, c_ptr_to_object) &
       bind(C,name='Eblock_matvec')
    ! ---------------------------------------------------------------------------
    ! Operator for E = dv/dps block of the Jacobian matrix
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use prim_derived_type_mod, only : derived_type
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use edge_mod, only              : edgevpack, edgevunpack
    use bndry_mod, only             : bndry_exchangev
    use control_mod, only           : tstep_type
    
    implicit none   

    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xvec(nxvec)
    integer(c_int), intent(in), value :: nxvec
    real(c_double), intent(out)       :: yvec(nxvec)
    type(c_ptr)                       :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np)        :: ps, MatVecProd
    real(kind=real_kind), dimension(np,np, nelemd) :: ptens
    
    integer :: nets, nete
    integer :: np1                
    integer :: kptr               
    integer :: i,j,k,ie,lx 
    ! ---------------------------------------------------------------------------

    call t_startf('Eblock_matvec')
    
    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*)  ">>> HOMME: Eblock_matvec"
#endif
    
    np1   = fptr%tl%np1   
    nets  = fptr%nets     
    nete  = fptr%nete     
    
    ! unpack 1D input array
    lx = 0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                fptr%base(ie)%state%ps_v(i,j,np1) = xvec(lx)
             end do
          end do
       end do
    end do

    do ie = nets,nete

       ps = fptr%base(ie)%state%ps_v(:,:,np1)
       
       call MatVecProd2Dx2D(MatVecProd, dPs_dPs, ps, ie)
       
       do j = 1,np
          do i = 1,np
             ptens(i,j,ie) = fptr%base(ie)%spheremp(i,j)*MatVecProd(i,j)
          end do
       end do

       ! pack tendencies into exchange buffer for DSS
       kptr = 0
       call edgeVpack(edge3p1, ptens(:,:,ie), 1, kptr, ie)
    end do

    call bndry_exchangeV(fptr%hybrid, edge3p1)

    ! unpack tendencies from exchange buffer for DSS
    kptr = 0
    do ie = nets,nete
       call edgeVunpack(edge3p1, ptens(:,:,ie), 1, kptr, ie)
    end do !ie

    ! pack 1D output array
    lx = 0
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             lx = lx+1
             yvec(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,ie) 
          end do
       end do
    end do
    
    call t_stopf('Eblock_matvec')

#ifdef DEBUG_PRINT_ON
    if (fptr%hybrid%masterthread) write(*,*)  "<<< HOMME: Eblock_matvec"
#endif
    
  end subroutine Eblock_matvec

  

  subroutine Fblock_matvec(xvec, nxvec, yvec, nyvec, c_ptr_to_object) &
       bind(C,name='Fblock_matvec')
    ! ---------------------------------------------------------------------------
    ! Operator for F = dT/dv block of the Jacobian matrix
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use prim_derived_type_mod, only : derived_type
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use edge_mod, only              : edgevpack, edgevunpack
    use bndry_mod, only             : bndry_exchangev
    use control_mod, only           : tstep_type

    implicit none

    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xvec(nxvec)
    integer(c_int), intent(in), value :: nxvec
    real(c_double), intent(out)       :: yvec(nyvec)
    integer(c_int), intent(in), value :: nyvec
    type(c_ptr)                       :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()

    real(kind=real_kind), dimension(np,np,nlev)        :: MatVecProd1, MatVecProd2
    real(kind=real_kind), dimension(np,np,nlev)        :: vel1, vel2
    real(kind=real_kind), dimension(np,np,nlev, nelemd) :: ttens
    
    integer :: nets, nete
    integer :: np1                
    integer :: kptr               
    integer :: i,j,k,ie,lx 
    ! ---------------------------------------------------------------------------

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*)  ">>> HOMME: Fblock_matvec"
#endif
    
    call t_startf('Fblock_matvec')
    
    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr
    
    np1   = fptr%tl%np1   
    nets  = fptr%nets     
    nete  = fptr%nete     
    
    ! unpack 1D input array
    lx = 0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,1,k,np1) = xvec(lx)
             end do
          end do
       end do
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,2,k,np1) = xvec(lx)
             end do
          end do
       end do
    end do

    do ie = nets,nete

       vel1 = fptr%base(ie)%state%v(:,:,1,:,np1)
       vel2 = fptr%base(ie)%state%v(:,:,2,:,np1)
       
       call MatVecProd3Dx3D(MatVecProd1, dT_dU, vel1, ie)
       call MatVecProd3Dx3D(MatVecProd2, dT_dV, vel2, ie)
       
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                ttens(i,j,k,ie) = fptr%base(ie)%spheremp(i,j)*(MatVecProd1(i,j,k) + MatVecProd2(i,j,k))
             end do
          end do
       end do

       ! pack tendencies into exchange buffer for DSS
       kptr = 0
       call edgeVpack(edge3p1, ttens(:,:,:,ie), nlev, kptr, ie)
    end do

    call bndry_exchangeV(fptr%hybrid, edge3p1)

    ! unpack tendencies from exchange buffer for DSS
    kptr = 0
    do ie = nets,nete
       call edgeVunpack(edge3p1, ttens(:,:,:,ie), nlev, kptr, ie)
    end do !ie

    ! pack 1D output array
    lx = 0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                yvec(lx) = fptr%base(ie)%rspheremp(i,j)*ttens(i,j,k,ie) 
             end do
          end do
       end do
    end do
    
    call t_stopf('Fblock_matvec')

#ifdef DEBUG_PRINT_ON
    if (fptr%hybrid%masterthread) write(*,*)  "<<< HOMME: Fblock_matvec"
#endif
    
  end subroutine Fblock_matvec


  
  subroutine Gblock_matvec(xvec, nxvec, yvec, nyvec, c_ptr_to_object) &
       bind(C,name='Gblock_matvec')
    ! ---------------------------------------------------------------------------
    ! Operator for G = dT/dps block of the Jacobian matrix
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use prim_derived_type_mod, only : derived_type
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use edge_mod, only              : edgevpack, edgevunpack
    use bndry_mod, only             : bndry_exchangev
    use control_mod, only           : tstep_type

    implicit none

    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xvec(nxvec)
    integer(c_int), intent(in), value :: nxvec
    real(c_double), intent(out)       :: yvec(nyvec)
    integer(c_int), intent(in), value :: nyvec
    type(c_ptr)                       :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr=>NULL()
    
    real(kind=real_kind), dimension(np,np)             :: ps
    real(kind=real_kind), dimension(np,np,nlev)        :: MatVecProd
    real(kind=real_kind), dimension(np,np,nlev, nelemd) :: ttens
    
    integer :: nets, nete
    integer :: np1
    integer :: kptr
    integer :: i,j,k,ie,lx
    ! ---------------------------------------------------------------------------

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*)  ">>> HOMME: Gblock_matvec"
#endif

    call t_startf('Gblock_matvec')
    
    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr
    
    np1  = fptr%tl%np1   
    nets = fptr%nets     
    nete = fptr%nete     
    
    ! unpack 1D input array
    lx = 0
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             lx = lx+1
             fptr%base(ie)%state%ps_v(i,j,np1) = xvec(lx)
          end do
       end do
    end do

    do ie = nets,nete

       ps = fptr%base(ie)%state%ps_v(:,:,np1)

       call MatVecProd3Dx2D(MatVecProd, dT_dPs, ps, ie)

       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                ttens(i,j,k,ie) = fptr%base(ie)%spheremp(i,j) * MatVecProd(i,j,k)
             end do
          end do
       end do

       kptr = 0
       call edgeVpack(edge3p1, ttens(:,:,:,ie), nlev, kptr, ie)
    end do

    call bndry_exchangeV(fptr%hybrid, edge3p1)

    kptr = 0
    do ie = nets,nete
       call edgeVunpack(edge3p1, ttens(:,:,:,ie), nlev, kptr, ie)
    end do

    ! pack 1D output array
    lx = 0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                yvec(lx) = fptr%base(ie)%rspheremp(i,j)*ttens(i,j,k,ie)
             end do
          end do
       end do
    end do

    call t_stopf('Gblock_matvec')

#ifdef DEBUG_PRINT_ON
    if (fptr%hybrid%masterthread) write(*,*)  "<<< HOMME: Gblock_matvec"
#endif
    
  end subroutine Gblock_matvec

  

  subroutine Hblock_matvec(xvec, nxvec, yvec, c_ptr_to_object) &
       bind(C,name='Hblock_matvec')
    ! ---------------------------------------------------------------------------
    ! Operator for H = dT/dT block of the Jacobian matrix
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use prim_derived_type_mod, only : derived_type
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use edge_mod, only              : edgevpack, edgevunpack
    use bndry_mod, only             : bndry_exchangev
    use control_mod, only           : tstep_type
    
    implicit none

    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xvec(nxvec)    
    integer(c_int), intent(in), value :: nxvec
    real(c_double), intent(out)       :: yvec(nxvec)    
    type(c_ptr)                       :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr => NULL()

    real(kind=real_kind), dimension(np,np,nlev)        :: Temp
    real(kind=real_kind), dimension(np,np,nlev)        :: MatVecProd
    real(kind=real_kind), dimension(np,np,nlev, nelemd) :: ttens
    
    real(kind=real_kind) :: dti
    
    integer :: nets, nete
    integer :: i, j, k, ie, lx
    integer :: np1
    integer :: kptr
    ! ---------------------------------------------------------------------------
    
    call t_startf('Hblock_matvec')
    
    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*)  ">>> HOMME: Hblock_matvec"
#endif
    
    np1   = fptr%tl%np1   
    nets  = fptr%nets     
    nete  = fptr%nete     

    dti = 1.0d0/fptr%dt
    
    if (tstep_type==12) then
       dti=(3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2
    endif

    ! unpack 1D input array
    lx = 0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                fptr%base(ie)%state%T(i,j,k,np1) = xvec(lx)
             end do
          end do
       end do
    end do

    do ie = nets,nete

       Temp = fptr%base(ie)%state%T(:,:,:,np1)

       ! compute matvec using full dTdT block
       call MatVecProd3Dx3D(MatVecProd, dT_dT, Temp, ie)
      
       ! compute matvec using diagonal of dTdT block
       !call MatVecProd3Dx3D_Diag(MatVecProd, dT_dT, Temp, ie)
       
       ! comput matvec using just dti
       !MatVecProd = dti * Temp
       
       do k=1,nlev
          do j=1,np
             do i=1,np

                ttens(i,j,k,ie) = fptr%base(ie)%spheremp(i,j) * MatVecProd(i,j,k)

             end do
          end do
       end do
       
       kptr = 0
       call edgeVpack(edge3p1, ttens(:,:,:,ie), nlev, kptr, ie)
    end do !ie

    call bndry_exchangeV(fptr%hybrid,edge3p1)

    kptr = 0
    do ie = nets,nete
       call edgeVunpack(edge3p1, ttens(:,:,:,ie), nlev, kptr, ie)
    end do

    ! pack first ttens component
    lx = 0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                yvec(lx) = fptr%base(ie)%rspheremp(i,j) * ttens(i,j,k,ie)
             end do
          end do
       end do
    end do

    call t_stopf('Hblock_matvec')

#ifdef DEBUG_PRINT_ON
    if (fptr%hybrid%masterthread) write(*,*)  "<<< HOMME: Hblock_matvec"
#endif
   
  end subroutine Hblock_matvec
  


  subroutine AinvB_matvec(xvec, nxvec, yvec, nyvec, c_ptr_to_object) &
       bind(C,name='AinvB_matvec')
    ! ---------------------------------------------------------------------------
    ! Operator for Ahat^{-1} B
    ! --------------------------------------------------------------------------- 
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use prim_derived_type_mod, only : derived_type
    use edge_mod, only              : edgevpack, edgevunpack
    use bndry_mod, only             : bndry_exchangev
    
    implicit none

    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xvec(nxvec)    
    integer(c_int), intent(in), value :: nxvec
    real(c_double), intent(out)       :: yvec(nyvec)    
    integer(c_int), intent(in), value :: nyvec
    type(c_ptr)                       :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr => NULL()

    real(kind=real_kind), dimension(np,np,2,nlev, nelemd) :: vtens
    real(kind=real_kind), dimension(np,np,nlev)          :: vel1, vel2
    real(kind=real_kind), dimension(np,np,nlev)          :: utend1, vtend1
    real(kind=real_kind), dimension(np,np)               :: ps
           
    integer :: nets, nete
    integer :: np1
    integer :: kptr
    integer :: i,j,k,ie,lx
    ! ---------------------------------------------------------------------------
    
    call t_startf('AinvB_matvec')
    
    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(6,*)  ">>> HOMME: AinvB_matvec"
#endif
    
    np1  = fptr%tl%np1   
    nets = fptr%nets     
    nete = fptr%nete     
    
    ! unpack 1D input array
    lx = 0
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             lx = lx+1
             fptr%base(ie)%state%ps_v(i,j,np1) = xvec(lx)
          end do
       end do
    end do

    ! ---------------------------------------------------------------------------
    ! compute B * ps_v
    ! ---------------------------------------------------------------------------

    do ie = nets,nete

       ps = fptr%base(ie)%state%ps_v(:,:,np1)

       call MatVecProd3Dx2D(utend1, dU_dPs, ps, ie)
       call MatVecProd3Dx2D(vtend1, dV_dPs, ps, ie)

       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                ! u vel
                vtens(i,j,1,k,ie) = fptr%base(ie)%spheremp(i,j) * utend1(i,j,k)

                ! v vel
                vtens(i,j,2,k,ie) = fptr%base(ie)%spheremp(i,j) * vtend1(i,j,k)
             end do
          end do
       end do

       kptr = 0
       call edgeVpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
    end do

    call bndry_exchangeV(fptr%hybrid,edge3p1)

    kptr = 0
    do ie = nets,nete
       call edgeVunpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                vtens(i,j,1,k,ie) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
             end do
          end do
       end do
    end do

    ! ---------------------------------------------------------------------------
    ! compute Ahat^{-1} * B * ps_v
    ! ---------------------------------------------------------------------------

    do ie = nets,nete

       vel1 = vtens(:,:,1,:,ie)
       vel2 = vtens(:,:,2,:,ie)
     
       ! compute MatVecProd for vertical advection
       call MatVecProd3Dx3D_DiagInv(utend1, dU_dU, vel1, ie)
       call MatVecProd3Dx3D_DiagInv(vtend1, dV_dV, vel2, ie)
       
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                ! u vel 
                vtens(i,j,1,k,ie) = fptr%base(ie)%spheremp(i,j)*utend1(i,j,k)
                
                ! v vel 
                vtens(i,j,2,k,ie) = fptr%base(ie)%spheremp(i,j)*vtend1(i,j,k)
             end do
          end do
       end do
       
       kptr = 0
       call edgeVpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
    end do !ie

    call bndry_exchangeV(fptr%hybrid, edge3p1)
    
    kptr = 0
    do ie = nets,nete
       call edgeVunpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
    end do

    ! ---------------------------------------------------------------------------
    ! pack 1D output array
    ! ---------------------------------------------------------------------------

    lx = 0
    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                yvec(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
             end do
          end do
       end do
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                lx = lx+1
                yvec(lx) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
             end do
          end do
       end do
    end do

    call t_stopf('AinvB_matvec')

#ifdef DEBUG_PRINT_ON
    if (fptr%hybrid%masterthread) write(6,*)  "<<< HOMME: AinvB_matvec"
#endif
   
  end subroutine AinvB_matvec
 


  subroutine SchurS_matvec(xvec, nxvec, yvec, c_ptr_to_object) &
       bind(C,name='SchurS_matvec')
    ! ---------------------------------------------------------------------------
    ! Schur Block: S_hat = E - D A_hat^{-1} B 
    !
    ! This version DOES NOT call the individual subroutines for E, D A_hat^{-1}, 
    ! B rather it performs the individual MatVecProds and has fewer exchanges
    ! ---------------------------------------------------------------------------
    use, intrinsic :: iso_c_binding

    use kinds, only                 : real_kind
    use prim_derived_type_mod, only : derived_type
    use dimensions_mod, only        : np, nlev, nvar, nelemd
    use edge_mod, only              : edgevpack, edgevunpack
    use bndry_mod, only             : bndry_exchangev
    use control_mod, only           : tstep_type
    
    implicit none

    ! --------------------------------Arguments----------------------------------
    real(c_double), intent(in)        :: xvec(nxvec)    
    integer(c_int), intent(in), value :: nxvec
    real(c_double), intent(out)       :: yvec(nxvec)    
    type(c_ptr)                       :: c_ptr_to_object
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    type(derived_type), pointer :: fptr => NULL()
    
    real(kind=real_kind), dimension(np,np)        :: ps
    real(kind=real_kind), dimension(np,np,nelemd) :: ptens1, ptens2

    real(kind=real_kind), dimension(np,np,nlev)          :: vel1, vel2
    real(kind=real_kind), dimension(np,np,2,nlev,nelemd) :: vtens

    real(kind=real_kind), dimension(np,np,nlev) :: MatVecProd1, MatVecProd2
    real(kind=real_kind), dimension(np,np)      :: MatVecProd3, MatVecProd4
    
    real(kind=real_kind) :: dti
    
    integer :: nets, nete
    integer :: np1
    integer :: kptr
    integer :: i,j,k,ie,lx
    ! ---------------------------------------------------------------------------
    
    call t_startf('SchurS_matvec')
    
    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr

#ifdef DEBUG_PRINT_ON   
    if (fptr%hybrid%masterthread) write(*,*)  ">>> HOMME: SchurS_matvec"
#endif

    nets  = fptr%nets     
    nete  = fptr%nete     
    np1   = fptr%tl%np1   

    dti = 1.0d0/fptr%dt 

    if (tstep_type==12) then
       dti=(3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2 in the time step coefficient
    endif
    
    ! unpack 1D input array
    lx = 0
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             lx = lx+1
             fptr%base(ie)%state%ps_v(i,j,np1) = xvec(lx)
          end do
       end do
    end do

    ! ---------------------------------------------------------------------------
    ! vtens = B * xvec (dv/dps)
    ! ---------------------------------------------------------------------------
    
    do ie = nets,nete

       ps = fptr%base(ie)%state%ps_v(:,:,np1)

       call MatVecProd3Dx2D(MatVecProd1, dU_dPs, ps, ie)
       call MatVecProd3Dx2D(MatVecProd2, dV_dPs, ps, ie)

       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                ! u vel
                vtens(i,j,1,k,ie) = fptr%base(ie)%spheremp(i,j) * MatVecProd1(i,j,k)

                ! v vel
                vtens(i,j,2,k,ie) = fptr%base(ie)%spheremp(i,j) * MatVecProd2(i,j,k)
             end do
          end do
       end do

       kptr = 0
       call edgeVpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
    end do

    call bndry_exchangeV(fptr%hybrid, edge3p1)

    kptr = 0
    do ie = nets,nete
       call edgeVunpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                vtens(i,j,1,k,ie) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
             end do
          end do
       end do
    end do

    ! ---------------------------------------------------------------------------
    ! apply rspheremp and vtens = A_hat^{-1} * vtens 
    ! ---------------------------------------------------------------------------

    do ie = nets,nete

       vel1 = vtens(:,:,1,:,ie)
       vel2 = vtens(:,:,2,:,ie)
       
       ! compute MatVecProd for vertical advection
       call MatVecProd3Dx3D_DiagInv(MatVecProd1, dU_dU, vel1, ie)
       call MatVecProd3Dx3D_DiagInv(MatVecProd2, dV_dV, vel2, ie)
       
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                ! u vel 
                vtens(i,j,1,k,ie) = fptr%base(ie)%spheremp(i,j)*MatVecProd1(i,j,k)
                
                ! v vel 
                vtens(i,j,2,k,ie) = fptr%base(ie)%spheremp(i,j)*MatVecProd2(i,j,k)
             end do
          end do
       end do
       
       kptr = 0
       call edgeVpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
    end do !ie

    call bndry_exchangeV(fptr%hybrid,edge3p1)

    kptr = 0
    do ie = nets,nete
       call edgeVunpack(edge3p1, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
    end do

    do ie = nets,nete
       do k = 1,nlev
          do j = 1,np
             do i = 1,np
                vtens(i,j,1,k,ie) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
             end do
          end do
       end do
    end do

    ! ---------------------------------------------------------------------------
    ! ptens1 = D * vetns (dps/dv)
    ! ---------------------------------------------------------------------------

    do ie = nets,nete
      
       vel1 = vtens(:,:,1,:,ie)
       vel2 = vtens(:,:,2,:,ie)
       
       call MatVecProd2Dx3D(MatVecProd3, dPs_dU, vel1, ie)
       call MatVecProd2Dx3D(MatVecProd4, dPs_dV, vel2, ie)
       
       do j = 1,np
          do i = 1,np
             ptens1(i,j,ie) = MatVecProd3(i,j) + MatVecProd4(i,j)
          end do
       end do
    end do
    
    ! ---------------------------------------------------------------------------
    ! ptens2 = E * xvec (dps/dps)
    ! ---------------------------------------------------------------------------

    do ie = nets,nete
    
       ps = fptr%base(ie)%state%ps_v(:,:,np1)
       
       call MatVecProd2Dx2D(MatVecProd3, dPs_dPs, ps, ie)
       
       do j = 1,np
          do i = 1,np
             ptens2(i,j,ie) = MatVecProd3(i,j)
          end do
       end do

    end do
  
    ! ---------------------------------------------------------------------------
    ! yvec = ptens2 - ptens1 
    ! ---------------------------------------------------------------------------

    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             ptens2(i,j,ie) = fptr%base(ie)%spheremp(i,j)*(ptens2(i,j,ie) - ptens1(i,j,ie))
          end do
       end do
         
       kptr = 0
       call edgeVpack(edge3p1, ptens2(:,:,ie), 1, kptr, ie)
    end do !ie
    
    call bndry_exchangeV(fptr%hybrid,edge3p1)
    
    kptr = 0
    do ie = nets,nete
       call edgeVunpack(edge3p1, ptens2(:,:,ie), 1, kptr, ie)
    end do

    ! pack 1D output array
    lx = 0
    do ie = nets,nete
       do j = 1,np
          do i = 1,np
             lx = lx+1
             yvec(lx) = fptr%base(ie)%rspheremp(i,j) * ptens2(i,j,ie)
          end do
       end do
    end do

    call t_stopf('SchurS_matvec')

#ifdef DEBUG_PRINT_ON
    if (fptr%hybrid%masterthread) write(*,*)  "<<< HOMME: SchurS_matvec"
#endif
    
  end subroutine SchurS_matvec



  ! =============================================================================
  ! The following subroutines are for computing matrix-vector products with
  ! the analytic Jacobian 
  ! =============================================================================

  

  subroutine MatVecProd3Dx3D(MatVecProd, matrix, vector, ie)
    ! ---------------------------------------------------------------------------
    ! Compute matrix vector product for 3D matrix and 3D variables
    ! ---------------------------------------------------------------------------
    use kinds, only          : real_kind
    use dimensions_mod, only : np, nlev, nelemd

    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(kind=real_kind), intent(out) :: MatVecProd(np,np,nlev)
    real(kind=real_kind), intent(in)  :: matrix(np,np,nlev,np,np,nlev,nelemd)
    real(kind=real_kind), intent(in)  :: vector(np,np,nlev)

    integer, intent(in) :: ie
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    integer :: a, b, c
    integer :: i, j, k
    ! ---------------------------------------------------------------------------

    call t_startf('MatVecProd3Dx3D')

    MatVecProd = 0.0d0
    
    ! loop over rows (residual equations / node indices)
    do k = 1,nlev
       do j = 1,np
          do i = 1,np
             
             ! loop over columns (derivatives)
             do c = 1,nlev
                do b = 1,np
                   do a = 1,np
                      MatVecProd(i,j,k) = MatVecProd(i,j,k) + matrix(a,b,c,i,j,k,ie) * vector(a,b,c)
                   end do
                end do
             end do
             
          end do
       end do
    end do

    call t_stopf('MatVecProd3Dx3D')
    
  end subroutine MatVecProd3Dx3D



  subroutine MatVecProd3Dx2D(MatVecProd, matrix, vector, ie)
    !----------------------------------------------------------------------------
    ! compute matrix vector product for 3D matrix and 2D variables
    !----------------------------------------------------------------------------
    use kinds, only          : real_kind
    use dimensions_mod, only : np, nlev, nelemd

    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(kind=real_kind), intent(out) :: MatVecProd(np,np,nlev)
    real(kind=real_kind), intent(in)  :: matrix(np,np,np,np,nlev,nelemd)
    real(kind=real_kind), intent(in)  :: vector(np,np)

    integer, intent(in) :: ie
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    integer :: a, b
    integer :: i, j, k
    ! ---------------------------------------------------------------------------

    call t_startf('MatVecProd3Dx2D')
    
    MatVecProd = 0.0d0
    
    ! loop over rows (residual equations / node indices)
    do k = 1,nlev
       do j = 1,np
          do i = 1,np
             
             ! loop over columns (derivatives)
             do b = 1,np
                do a = 1,np
                   MatVecProd(i,j,k) = MatVecProd(i,j,k) + matrix(a,b,i,j,k,ie) * vector(a,b)
                end do
             end do
             
          end do
       end do
    end do

    call t_stopf('MatVecProd3Dx2D')
    
  end subroutine MatVecProd3Dx2D



  subroutine MatVecProd2Dx3D(MatVecProd, matrix, vector, ie)
    ! ---------------------------------------------------------------------------
    ! Compute matrix vector product of a 2Dx3D matrix with a 3D vector
    ! ---------------------------------------------------------------------------
    use kinds, only          : real_kind
    use dimensions_mod, only : np, nlev, nelemd

    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(kind=real_kind), intent(out) :: MatVecProd(np,np)
    real(kind=real_kind), intent(in)  :: matrix(np,np,nlev,np,np,nelemd)
    real(kind=real_kind), intent(in)  :: vector(np,np,nlev)

    integer, intent(in) :: ie
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    integer :: a, b, c
    integer :: i, j
    ! ---------------------------------------------------------------------------

    call t_startf('MatVecProd2Dx3D')
    
    MatVecProd = 0.0d0
    
    ! loop over rows (residual equations / node indices)
    do j = 1,np
       do i = 1,np
          
          ! loop over columns (derivatives)
          do c = 1,nlev
             do b = 1,np
                do a = 1,np
                   MatVecProd(i,j) = MatVecProd(i,j) + matrix(a,b,c,i,j,ie) * vector(a,b,c)
                end do
             end do
          end do
             
       end do
    end do

    call t_stopf('MatVecProd2Dx3D')
    
  end subroutine MatVecProd2Dx3D



  subroutine MatVecProd2Dx2D(MatVecProd, matrix, vector, ie)
    ! ---------------------------------------------------------------------------
    ! Compute matrix vector product for at 2Dx2D matrix with as 2D vector
    ! ---------------------------------------------------------------------------
    use kinds, only          : real_kind
    use dimensions_mod, only : np, nlev, nelemd

    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(kind=real_kind), intent(out) :: MatVecProd(np,np)
    real(kind=real_kind), intent(in)  :: matrix(np,np,np,np,nelemd)
    real(kind=real_kind), intent(in)  :: vector(np,np)

    integer, intent(in) :: ie
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    integer :: a, b
    integer :: i, j
    ! ---------------------------------------------------------------------------

    call t_startf('MatVecProd2Dx2D')
    
    MatVecProd = 0.0d0
    
    ! loop over rows (residual equations / node indices)
    do j = 1,np
       do i = 1,np
          
          ! loop over columns (derivatives)
          do b = 1,np
             do a = 1,np
                MatVecProd(i,j) = MatVecProd(i,j) + matrix(a,b,i,j,ie) * vector(a,b)
             end do
          end do
       
       end do
    end do

    call t_stopf('MatVecProd2Dx2D')
    
  end subroutine MatVecProd2Dx2D



  subroutine MatVecProd3Dx3D_Diag(MatVecProd, matrix, vector, ie)
    ! ---------------------------------------------------------------------------
    ! Compute matrix vector product for 3D matrix and 3D variables
    ! ---------------------------------------------------------------------------
    use kinds, only          : real_kind
    use dimensions_mod, only : np, nlev, nelemd

    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(kind=real_kind), intent(out) :: MatVecProd(np,np,nlev)
    real(kind=real_kind), intent(in)  :: matrix(np,np,nlev,np,np,nlev,nelemd)
    real(kind=real_kind), intent(in)  :: vector(np,np,nlev)

    integer, intent(in) :: ie
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    integer :: i, j, k
    ! ---------------------------------------------------------------------------

    call t_startf('MatVecProd3Dx3D_Diag')

    MatVecProd = 0.0d0
   
    do k = 1,nlev
       do j = 1,np
          do i = 1,np
             
             MatVecProd(i,j,k) = matrix(i,j,k,i,j,k,ie) * vector(i,j,k)
             
          end do
       end do
    end do

    call t_stopf('MatVecProd3Dx3D_Diag')
    
  end subroutine MatVecProd3Dx3D_Diag



  subroutine MatVecProd3Dx3D_DiagInv(MatVecProd, matrix, vector, ie)
    ! ---------------------------------------------------------------------------
    ! Compute matrix vector product for 3D matrix and 3D variables
    ! ---------------------------------------------------------------------------
    use kinds, only          : real_kind
    use dimensions_mod, only : np, nlev, nelemd

    implicit none
    
    ! --------------------------------Arguments----------------------------------
    real(kind=real_kind), intent(out) :: MatVecProd(np,np,nlev)
    real(kind=real_kind), intent(in)  :: matrix(np,np,nlev,np,np,nlev,nelemd)
    real(kind=real_kind), intent(in)  :: vector(np,np,nlev)

    integer, intent(in) :: ie
    ! ---------------------------------------------------------------------------

    ! -----------------------------Local workspace-------------------------------
    integer :: i, j, k
    ! ---------------------------------------------------------------------------

    call t_startf('MatVecProd3Dx3D_DiagInv')

    MatVecProd = 0.0d0
   
    do k = 1,nlev
       do j = 1,np
          do i = 1,np
             
             MatVecProd(i,j,k) = vector(i,j,k) / matrix(i,j,k,i,j,k,ie)
             
          end do
       end do
    end do

    call t_stopf('MatVecProd3Dx3D_DiagInv')
    
  end subroutine MatVecProd3Dx3D_DiagInv

end module prim_jacobian_dense_mod

#endif 
! TRILINOS
