!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_sparse_pcg.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module cism_sparse_pcg

  use glimmer_global, only: dp
  use glimmer_sparse_type, only: sparse_matrix_type
  use parallel

  implicit none

  private
  public :: pcg_solver_structured

!WHL - debug
    logical :: verbose_pcg

!WHL - debug
    integer, parameter :: &
        ndiagmax = 10      ! number of values to print out for debugging

!WHL - debug

   integer, parameter :: &
!       itest = 24, jtest = 17, ktest = 1
!       itest = 17, jtest = 10, ktest = 1
       itest = 10, jtest = 17, ktest = 1

   integer, parameter :: &
!       ntest = 2371   ! nodeID for (24,17,1)
!       ntest = 411    ! nodeID for (17,10,1)
       ntest = 1771   ! nodeID for (10,17,1)

    integer, parameter :: rtest = 0    ! rank for (17,10) with 4 procs

contains

!****************************************************************************

  subroutine pcg_solver_structured(nx,        ny,            &
                                   nz,        nhalo,         &
                                   indxA,     active_vertex, &
                                   Auu,       Auv,           &
                                   Avu,       Avv,           &
                                   bu,        bv,            &
                                   xu,        xv,            &
                                   err,       niters,        &
                                   verbose) 

    !---------------------------------------------------------------
    !  This subroutine uses a preconditioned conjugate-gradient solver
    !  to solve the equation $Ax=b$.
    !  Convergence is checked every {\em ncheck} steps.
    !
    !  It is based on the barotropic solver in the POP ocean model 
    !  (author Phil Jones, LANL).  Input and output arrays are located
    !  on a structured (i,j,k) grid as defined in the glissade_velo_higher
    !  module.  The global matrix is sparse, but its nonzero elements
    !  are stored in four dense matrices called Auu, Avv, Auv, and Avu.
    !  Each matrix has 3x3x3 = 27 potential nonzero elements per
    !  node (i,j,k).
    !
    !  The current preconditioning options are
    !  (0) no preconditioning
    !  (1) diagonal preconditioning
    !  (2) preconditioning using a physics-based SIA solver
    ! 
    !  For the dome test case with higher-order dynamics, option (2) is best. 
    !
    !  Here is a schematic of the method implemented below for solving Ax = b:
    !
    !  r0 = b - A*x0
    !  d0 = 0
    !  eta0 = 1
    !
    !  while (not converged)
    !     z = (Minv)r
    !     eta1 = (r,z)
    !     beta = eta1/eta0
    !     d = z + beta*d
    !     eta0 = eta1
    !     y = Ad
    !     eta2 = (d,y)
    !     alpha = eta1/eta2
    !     x = x + alpha*d
    !     r = r - alpha*y (or occasionally, r = b - Ax)
    !     Check for convergence: err = sqrt(r,r)/sqrt(b,b) < tolerance
    !  end while
    !
    !  where x = solution (initial value = x0)
    !        d = conjugate direction vector (initial value = d0)
    !        r = residual vector (initial value = r0)
    !     Minv = inverse of preconditioning matrix M
    !            (can be implemented by solving Mz = r without forming Minv)
    !    (r,z) = dot product of vectors r and z
    !            and similarly for (d,y)
    !       
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions (for scalars)
                                ! velocity grid has dimensions (nx-1,ny-1)
       nz,                   &  ! number of vertical levels where velocity is computed
       nhalo                    ! number of halo layers (for scalars)

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
       indxA                 ! maps relative (x,y,z) coordinates to an index between 1 and 27

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F
 
    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 1st dimension = 27 (node and its nearest neighbors in x, y and z direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::   &
       xu, xv             ! u and v components of solution (i.e., uvel and vvel)

    real(dp), intent(out) ::  &
       err                               ! error (L2 norm of residual) in final solution

    integer, intent(out) ::   &
       niters                            ! iterations needed to solution

    logical, intent(in), optional :: &
       verbose                           ! if true, print diagnostic output
    
    !---------------------------------------------------------------
    ! Local variables and parameters   
    !---------------------------------------------------------------

    integer ::  i, j, k      ! grid indices
    integer ::  iA, jA, kA   ! grid offsets ranging from -1 to 1
    integer ::  m            ! matrix element index
    integer ::  n            ! iteration counter

    real(dp) ::           &
       eta0, eta1, eta2,  &! scalar inner product results
       alpha,             &! eta1/eta2 = term in expression for new residual
       beta                ! eta1/eta0 = term in expression for new direction vector

    ! vectors (each of these is split into u and v components)
    real(dp), dimension(nz,nx-1,ny-1) ::  &
       Adiagu, Adiagv,    &! diagonal terms of matrices Auu and Avv
       ru, rv,            &! residual vector (b-Ax)
       du, dv,            &! conjugate direction vector
       yu, yv,            &! result of a matvec multiply
       zu, zv,            &! solution of Mz = r
       work0u, work0v      ! cg intermediate results

    real(dp) ::  &
       L2_resid,          &! L2 norm of residual vector Ax-b
       L2_rhs              ! L2 norm of rhs vector b
                           ! solver converges when L2_resid/L2_rhs < tolerance

    real(dp), dimension(-1:1,nz,nx-1,ny-1) ::  &
       Muu, Mvv            ! simplified SIA matrices for preconditioning

    !---------------------------------------------------------------
    ! Solver parameters
    ! TODO: Pass these in as arguments?
    !---------------------------------------------------------------

    real(dp), parameter ::   &
       tolerance = 1.d-11    ! tolerance for linear solver

    ! WHL: With SIA preconconditioning, up to 210 iterations are needed
    !       for convergence with ismip-hom test A, 20 km resolution
    !      More are probably needed at higher resolutions.
    ! TODO: Implement a better preconditioner for higher-order problems!
 
    integer, parameter ::    &
!       maxiters = 200        ! max number of linear iterations before quitting
       maxiters = 300        ! max number of linear iterations before quitting
                             
    integer, parameter :: &
       solv_ncheck = 1       ! check for convergence every solv_ncheck iterations
                             !TODO - See if performance improves for less frequent checks

    integer, parameter :: &
       solv_resid = 50       ! solve for residual as r = b - Ax every solv_resid iterations
                             ! else compute residual as r = r - alpha*(Ad)
    integer, parameter :: &
       precond = 2           ! = 0 for no preconditioning
                             ! = 1 for diagonal preconditioning (does not work well for SIA-dominated flow)
                             ! = 2 for preconditioning with SIA solver (works well for SIA-dominated flow)

    if (present(verbose)) then
       verbose_pcg = verbose
    else
       verbose_pcg = .true.   ! for debugging
    endif

    if (verbose_pcg .and. this_rank==rtest) then
       print*, ' '
       print*, 'In structured pcg solver'
       print*, 'tolerance, maxiters =', tolerance, maxiters
    endif

    ! Set up preconditioning

    if (precond == 1) then    ! form diagonal matrix for preconditioning

       m = indxA(0,0,0)
       Adiagu(:,:,:) = Auu(m,:,:,:)
       Adiagv(:,:,:) = Avv(m,:,:,:)

       if (verbose_pcg .and. this_rank==rtest) then
          print*, ' '
          print*, 'Using diagonal solver for preconditioning'
!          print*, ' '
!          print*, 'i, j, k, diagonal entries, initial guess, residual:'
!          m = 0
!          do j = 1, ny-1
!          do i = 1, nx-1
!          do k = 1, nz
!             if (Adiagu(k,i,j) > 0.d0 .and. m < ndiagmax) then
!                m = m + 2
!                print*, i, j, k, Adiagu(k,i,j), xu(k,i,j), ru(k,i,j)
!                print*, i, j, k, Adiagv(k,i,j), xv(k,i,j), rv(k,i,j)
!             else
!                exit
!             endif
!          enddo
!          enddo
!          enddo
       endif  ! verbose_pcg

    elseif (precond == 2) then  ! form SIA matrices Muu and Mvv with vertical coupling only

       Muu(:,:,:,:) = 0.d0
       Mvv(:,:,:,:) = 0.d0

       do j = 1, ny-1
       do i = 1, nx-1
       do k = 1, nz
             ! Remove horizontal coupling by summing over iA = -1:1, jA = -1:1 for each layer
             ! Note: In dome tests this is slightly better than setting Muu and Muv to
             !       the matrix elements corresponding to iA = jA = 0

           Muu(-1,k,i,j) = sum(Auu(1:9,k,i,j))
           Mvv(-1,k,i,j) = sum(Avv(1:9,k,i,j))

           Muu( 0,k,i,j) = sum(Auu(10:18,k,i,j))
           Mvv( 0,k,i,j) = sum(Avv(10:18,k,i,j))

           Muu( 1,k,i,j) = sum(Auu(19:27,k,i,j))
           Mvv( 1,k,i,j) = sum(Avv(19:27,k,i,j))

       enddo
       enddo
       enddo

!WHL - debug
       if (verbose_pcg .and. this_rank==rtest) then
          print*, ' '
          print*, 'Using easy SIA solver for preconditioning'
!          i = itest
!          j = jtest
!          print*, ' '
!          print*, 'i, j =', i, j
!          print*, ' '
!          print*, 'k, Muu(-1:1):'
!          do k = 1, nz
!             print*, k, Muu(-1:1,k,i,j)
!          enddo
!          print*, ' '
!          print*, 'k, Mvv(-1:1):'
!          do k = 1, nz
!             print*, k, Mvv(-1:1,k,i,j)
!          enddo
       endif

    endif      ! precond

    ! Compute initial residual and initialize the direction vector d
    ! Note: The matrix A must be complete for all rows corresponding to locally 
    !        owned vertices, and x must have the correct values in
    !        halo vertices bordering the locally owned vertices.
    !       Then y = Ax will be correct for locally owned vertices.

    !TODO - Could use active_vertex array to set up indirect addressing
    !       and avoid an 'if' in the matvec subroutine

    ! Halo update for x (initial guess for velocity solution)

    call staggered_parallel_halo(xu)
    call staggered_parallel_halo(xv)

    ! Compute y = Ax

    call matvec_multiply_structured(nx,        ny,            &
                                    nz,        nhalo,         &
                                    indxA,     active_vertex, &
                                    Auu,       Auv,           &
                                    Avu,       Avv,           &
                                    xu,        xv,            &
                                    yu,        yv)

    ! Compute the initial residual r(0) = b - Ax(0)
    ! This will be correct for locally owned vertices.

    ru(:,:,:) = bu(:,:,:) - yu(:,:,:)
    rv(:,:,:) = bv(:,:,:) - yv(:,:,:)

    ! Initialize scalars and vectors

    niters = maxiters 
    eta0 = 1.d0

    du(:,:,:) = 0.d0
    dv(:,:,:) = 0.d0

    zu(:,:,:) = 0.d0
    zv(:,:,:) = 0.d0

    ! Compute the L2 norm of the RHS vectors
    ! (Goal is to obtain L2_resid/L2_rhs < tolerance)

    work0u(:,:,:) = bu(:,:,:)*bu(:,:,:)    ! terms of dot product (b, b)
    work0v(:,:,:) = bv(:,:,:)*bv(:,:,:)

    ! find global sum of the squared L2 norm

    call global_sum_staggered(nx,     ny,     &
                              nz,     nhalo,  &
                              L2_rhs,         &
                              work0u, work0v)

    ! take square root

    L2_rhs = sqrt(L2_rhs)       ! L2 norm of RHS
    if (verbose_pcg .and. this_rank==rtest) print*, 'Global L2_rhs =', L2_rhs

    ! Iterate to solution

    iter_loop: do n = 1, maxiters

       ! Compute (PC)r = solution z of Mz = r

       if (precond == 0) then      ! no preconditioning

           zu(:,:,:) = ru(:,:,:)         ! PC(r) = r     
           zv(:,:,:) = rv(:,:,:)         ! PC(r) = r    

       elseif (precond == 1 ) then  ! diagonal preconditioning

          do j = 1, ny-1
          do i = 1, nx-1
          do k = 1, nz
             if (Adiagu(k,i,j) /= 0.d0) then
                zu(k,i,j) = ru(k,i,j) / Adiagu(k,i,j)   ! PC(r), where PC is formed from diagonal elements of A
             else                                        
                zu(k,i,j) = 0.d0
             endif
             if (Adiagv(k,i,j) /= 0.d0) then
                zv(k,i,j) = rv(k,i,j) / Adiagv(k,i,j)  
             else                                        
                zv(k,i,j) = 0.d0
             endif
          enddo    ! k
          enddo    ! i
          enddo    ! j

!WHL - debug
          if (verbose_pcg .and. this_rank==rtest) then
!             i = itest
!             j = jtest
!             print*, ' '
!             print*, 'Diagonal preconditioning, i, j =', i, j
!             print*, ' '
!             print*, 'k, ru, PC(ru):'
!             do k = 1, nz
!                print*, k, ru(k,i,j), zu(k,i,j)
!             enddo
!             print*, ' '
!             print*, 'k, rv, PC(rv):'
!             do k = 1, nz
!                print*, k, rv(k,i,j), zv(k,i,j)
!             enddo
           endif

       elseif (precond == 2) then   ! local vertical shallow-ice solver for preconditioning

          call easy_sia_solver(nx,   ny,   nz,        &
                               active_vertex,         &
                               Muu,  ru,   zu)      ! solve Muu*zu = ru 

          call easy_sia_solver(nx,   ny,   nz,        &
                               active_vertex,         &
                               Mvv,  rv,   zv)      ! solve Mvv*zv = rv 

!WHL - debug
          if (verbose_pcg .and. this_rank==rtest) then
!             i = itest
!             j = jtest
!             print*, ' '
!             print*, 'SIA preconditioning, i, j =', i, j
!             print*, ' '
!             print*, 'k, ru, PC(ru):'
!             do k = 1, nz
!                print*, k, ru(k,i,j), zu(k,i,j)
!             enddo
!             print*, ' '
!             print*, 'k, rv, PC(rv):'
!             do k = 1, nz
!                print*, k, rv(k,i,j), zv(k,i,j)
!             enddo
           endif

       endif    ! precond

       ! Compute the dot product eta1 = (r, PC(r))

       work0u(:,:,:) = ru(:,:,:)*zu(:,:,:)    ! terms of dot product (r, PC(r))
       work0v(:,:,:) = rv(:,:,:)*zv(:,:,:)    

       call global_sum_staggered(nx,     ny,     &
                                 nz,     nhalo,  &
                                 eta1,           &
                                 work0u, work0v)

       ! Update the conjugate direction vector d

       beta = eta1/eta0

       if (verbose_pcg .and. this_rank==rtest) then
!          print*, 'eta0, eta1, beta =', eta0, eta1, beta
       endif

       du(:,:,:) = zu(:,:,:) + beta*du(:,:,:)       ! d_(i+1) = PC(r_(i+1)) + beta_(i+1)*d_i
       dv(:,:,:) = zv(:,:,:) + beta*dv(:,:,:)       !
                                                    !                    (r_(i+1), PC(r_(i+1)))
                                                    ! where beta_(i+1) = --------------------  
                                                    !                        (r_i, PC(r_i)) 
                                                    ! Initially eta0 = 1  
                                                    ! For n >=2, eta0 = old eta1

       ! Halo update for d

       call staggered_parallel_halo(du)
       call staggered_parallel_halo(dv)
  
       ! Compute y = A*d
       ! This is the one matvec multiply required for each iteration

       call matvec_multiply_structured(nx,        ny,            &
                                       nz,        nhalo,         &
                                       indxA,     active_vertex, &
                                       Auu,       Auv,           &
                                       Avu,       Avv,           &
                                       du,        dv,            &
                                       yu,        yv)

       ! Copy old eta1 = (r, PC(r)) to eta0

       eta0 = eta1               ! (r_(i+1), PC(r_(i+1))) --> (r_i, PC(r_i)) 

       ! Compute the dot product eta2 = (d, A*d)

       work0u(:,:,:) = yu(:,:,:) * du(:,:,:)       ! terms of dot product (d, Ad)
       work0v(:,:,:) = yv(:,:,:) * dv(:,:,:)

       call global_sum_staggered(nx,     ny,     &
                                 nz,     nhalo,  &
                                 eta2,           &
                                 work0u, work0v)

       ! Compute alpha

                              !          (r, PC(r))
       alpha = eta1/eta2      ! alpha = ----------
                              !          (d, A*d)
       

       if (verbose_pcg .and. this_rank==rtest) then
!          print*, '(r, PC(r)) =', eta1
!          print*, '(d, Ad) =', eta2
!          print*, 'alpha =', alpha
       endif

       ! Compute the new solution

       xu(:,:,:) = xu(:,:,:) + alpha * du(:,:,:)    ! new solution, x_(i+1) = x_i + alpha*d
       xv(:,:,:) = xv(:,:,:) + alpha * dv(:,:,:)

       ! Compute the new residual

       ! The cheap way is to take r = r - alpha*(Ad), since Ad is already available.
       ! The expensive way is to take r = b - Ax, which requires a new computation of Ax.
       ! Most of the time we do things the cheap way, but occasionally we use the expensive way
       !  to avoid cumulative loss of accuracy.

       if (mod(n,solv_resid) == 0) then    ! r = b - Ax every solv_resid iterations

          ! Halo update for x

          call staggered_parallel_halo(xu)
          call staggered_parallel_halo(xv)

          ! Compute y = Ax
           
          call matvec_multiply_structured(nx,        ny,            &
                                          nz,        nhalo,         &
                                          indxA,     active_vertex, &
                                          Auu,       Auv,           &
                                          Avu,       Avv,           &
                                          xu,        xv,            &
                                          yu,        yv)

          ! Compute residual r = b - Ax

          ru(:,:,:) = bu(:,:,:) - yu(:,:,:)
          rv(:,:,:) = bv(:,:,:) - yv(:,:,:)

       else    ! r = r - alpha*(Ad)

          ru(:,:,:) = ru(:,:,:) - alpha * yu(:,:,:)    ! new residual, r_(i+1) = r_i - alpha*(Ad)
          rv(:,:,:) = rv(:,:,:) - alpha * yv(:,:,:)

       endif

       if (verbose_pcg .and. this_rank==rtest) then
!          i = itest
!          j = jtest
!          print*, ' '
!          print*, 'i, j, =', i, j
!          print*, 'k, xu, ru:'
!          do k = 1, nz
!             print*, k, xu(k,i,j), ru(k,i,j)
!          enddo
       endif

       ! Check for convergence

       if (mod(n,solv_ncheck) == 0) then

          ! Compute squared L2 norm of (r, r)

          work0u(:,:,:) = ru(:,:,:)*ru(:,:,:)   ! terms of dot product (r, r)
          work0v(:,:,:) = rv(:,:,:)*rv(:,:,:)

          call global_sum_staggered(nx,     ny,       &
                                    nz,     nhalo,    &
                                    L2_resid,         &
                                    work0u, work0v)

          ! take square root
          L2_resid = sqrt(L2_resid)       ! L2 norm of residual

          err = L2_resid/L2_rhs           ! normalized error

          if (verbose_pcg .and. main_task) then
!             print*, ' '
             print*, 'iter, L2_resid, error =', n, L2_resid, err
          endif

          if (err < tolerance) then
             niters = n
             exit iter_loop
          endif            

       endif    ! solv_ncheck

    enddo iter_loop

!WHL - Without good preconditioning, convergence is slow, but the solution after maxiters might be good enough.
 
    if (niters == maxiters) then
       if (main_task) then
          print*, 'Glissade PCG solver not converged'
          print*, 'niters, err, tolerance:', niters, err, tolerance
       endif
!!!     stop   ! TODO - Abort cleanly
    endif

  end subroutine pcg_solver_structured

!****************************************************************************

  subroutine global_sum_staggered(nx,     ny,      &
                                  nz,     nhalo,   &
                                  global_sum,      &
                                  work1,  work2)

     ! Sum one or two local arrays on the staggered grid, then take the global sum.

     integer, intent(in) :: &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
       nz,                 &  ! number of vertical layers at which velocity is computed
       nhalo                  ! number of halo layers (for scalars)

     real(dp), intent(out) :: global_sum   ! global sum

     real(dp), intent(in), dimension(nz,nx-1,ny-1) :: work1            ! local array
     real(dp), intent(in), dimension(nz,nx-1,ny-1), optional :: work2  ! local array

     integer :: i, j, k
     real(dp) :: local_sum

     local_sum = 0.d0

     ! sum over locally owned velocity points

     if (present(work2)) then
        do j = nhalo+1, ny-nhalo
        do i = nhalo+1, nx-nhalo
        do k = 1, nz
           local_sum = local_sum + work1(k,i,j) + work2(k,i,j)
        enddo
        enddo
        enddo
     else
        do j = nhalo+1, ny-nhalo
        do i = nhalo+1, nx-nhalo
        do k = 1, nz
           local_sum = local_sum + work1(k,i,j)    
        enddo
        enddo
        enddo
     endif

     ! take the global sum

     global_sum = parallel_reduce_sum(local_sum)


     if (verbose_pcg .and. this_rank==rtest) then
!        print*, 'this_rank, local sum =', this_rank, local_sum 
!        print*, 'global sum =', global_sum 
     endif

    end subroutine global_sum_staggered

!****************************************************************************

  subroutine matvec_multiply_structured(nx,        ny,            &
                                        nz,        nhalo,         &
                                        indxA,     active_vertex, &
                                        Auu,       Auv,           &
                                        Avu,       Avv,           &
                                        xu,        xv,            &
                                        yu,        yv)

    !---------------------------------------------------------------
    ! Compute the matrix-vector product $y = Ax$.
    !
    ! The A matrices should have complete matrix elements for all
    !  rows corresponding to locally owned vertices.
    ! The terms of x should be correct for all locally owned vertices
    !  and also for all halo vertices adjacent to locally owned vertices.
    ! The resulting y will then be correct for locally owned vertices.
    !
    ! TODO: It's likely that more time will be spent in this subroutine than
    !       in any other, so we should think about how to optimize it,
    !       especially for GPUs.
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
       nz,                 &  ! number of vertical layers at which velocity is computed
       nhalo                  ! number of halo layers (for scalars)

    integer, dimension(-1:1,-1:1,-1:1), intent(in) :: &
       indxA                 ! maps relative (x,y,z) coordinates to an index between 1 and 27
    
    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 1st dimension = 27 (node and its nearest neighbors in x, y and z direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::   &
       xu, xv             ! current guess for solution


    real(dp), dimension(nz,nx-1,ny-1), intent(out) ::  &
       yu, yv             ! y = Ax

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: i, j, k, m
    integer :: iA, jA, kA
    
!WHL= debug
    real(dp) :: maxyu, maxyv
    integer :: imaxu, jmaxu, kmaxu
    integer :: imaxv, jmaxv, kmaxv

    ! Initialize the result vector.

    yu(:,:,:) = 0.d0
    yv(:,:,:) = 0.d0

    ! Compute y = Ax

    ! Loop over locally owned vertices

    do j = nhalo+1, ny-nhalo
    do i = nhalo+1, nx-nhalo

       !TODO - Use indirect addressing to avoid an 'if' here?
       if (active_vertex(i,j)) then

          do k = 1, nz

             !TODO - Replace these three short loops with long multadds for better GPU efficiency?
             do kA = -1,1
             do jA = -1,1
             do iA = -1,1

             !TODO - Can we somehow get rid of this 'if' statement and still keep the xu/xv indices in bounds?
                if ( (k+kA >= 1 .and. k+kA <= nz)         &
                                .and.                     &
                     (i+iA >= 1 .and. i+iA <= nx-1)       &
                                .and.                     &
                     (j+jA >= 1 .and. j+jA <= ny-1) ) then

                   m = indxA(iA,jA,kA)

                   yu(k,i,j) = yu(k,i,j)   &
                             + Auu(m,k,i,j)*xu(k+kA,i+iA,j+jA)  &
                             + Auv(m,k,i,j)*xv(k+kA,i+iA,j+jA)

                   yv(k,i,j) = yv(k,i,j)   &
                             + Avu(m,k,i,j)*xu(k+kA,i+iA,j+jA)  &
                             + Avv(m,k,i,j)*xv(k+kA,i+iA,j+jA)
 
                endif   ! k+kA, i+iA, j+jA in bounds

             enddo   ! kA
             enddo   ! iA
             enddo   ! jA

          enddo   ! k

       endif   ! active_vertex

    enddo   ! i
    enddo   ! j

    if (verbose_pcg) then
!       print*, ' '
!       print*, 'Find max yu, yv'
!       maxyu = -9999.d0
!       maxyv = -9999.d0
!       do j = 1, ny-1
!       do i = 1, nx-1
!       do k = 1, nz
!          if (yu(k,i,j) > maxyu) then
!             maxyu = yu(k,i,j)
!             imaxu = i
!             jmaxu = j
!             kmaxu = k
!          endif
!          if (yv(k,i,j) > maxyv) then
!             maxyv = yv(k,i,j)
!             imaxv = i
!             jmaxv = j
!             kmaxv = k
!          endif
!       enddo
!       enddo
!       enddo

    endif  ! verbose_pcg

  end subroutine matvec_multiply_structured
 
!****************************************************************************

  subroutine easy_sia_solver(nx,   ny,   nz,    &
                             active_vertex,     & 
                             A,    b,    x)

    !---------------------------------------------------------------
    ! Solve the problem Ax = b where A is a local shallow-ice matrix,
    !  with coupling in the vertical but not the horizontal.
    ! We simply solve a tridiagonal matrix for each column.
    !---------------------------------------------------------------
   
    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
       nz                     ! number of vertical levels
           
    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(-1:1,nz,nx-1,ny-1), intent(in) ::   &
       A                      ! matrix with vertical coupling only
                              ! 1st dimension = node and its upper and lower neighbors

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::   &
       b                      ! right-hand side

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::   &
       x                      ! solution

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    real(dp), dimension(nz) ::  &
       sbdiag,         &      ! subdiagonal matrix entries
       diag,           &      ! diagonal matrix entries
       spdiag,         &      ! superdiagonal matrix entries
       rhs,            &      ! right-hand side
       soln                   ! tridiagonal solution

    integer :: i, j, k

    do j = 1, ny-1
    do i = 1, nx-1

       if (active_vertex(i,j)) then

          ! initialize rhs and solution 

          rhs(:)  = b(:,i,j)
          soln(:) = x(:,i,j)

          ! top layer

          k = 1
          sbdiag(k) = 0.d0
          diag(k)   = A(0,k,i,j)
          spdiag(k) = A(1,k,i,j)

          ! intermediate layers

          do k = 2, nz-1
             sbdiag(k) = A(-1,k,i,j)
             diag(k)   = A( 0,k,i,j)
             spdiag(k) = A( 1,k,i,j)
          enddo

          ! bottom layer

          k = nz 
          sbdiag(k) = A(-1,k,i,j)
          diag(k)   = A( 0,k,i,j)
          spdiag(k) = 0.d0           
         
          ! solve

          call tridiag_solver(nz,    sbdiag,   &
                              diag,  spdiag,   &
                              rhs,   soln)

          x(:,i,j) = soln(:)

       endif  ! active_vertex

    enddo     ! i
    enddo     ! j

  end subroutine easy_sia_solver
  
!****************************************************************************

  subroutine tridiag_solver(order,    sbdiag,   &
                            diag,     spdiag,   &
                            rhs,      soln)

    !---------------------------------------------------------------
    ! Solve a 1D tridiagonal matrix problem.
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
       order            ! matrix dimension

    real(dp), dimension(order), intent(in) :: &
       sbdiag,        & ! sub-diagonal matrix elements
       diag,          & ! diagonal matrix elements
       spdiag,        & ! super-diagonal matrix elements
       rhs              ! right-hand side

    real(dp), dimension(order), intent(inout) :: &
       soln             ! solution vector

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: &
       k               ! row counter

    real(dp) :: &
       beta            ! temporary matrix variable

    real(dp), dimension(order) :: &
       gamma           ! temporary matrix variable

    ! Solve

    beta = diag(1)
    soln(1) = rhs(1) / beta

    do k = 2, order
       gamma(k) = spdiag(k-1) / beta
       beta = diag(k) - sbdiag(k)*gamma(k)
       soln(k) = (rhs(k) - sbdiag(k)*soln(k-1)) / beta
    enddo

    do k = order-1, 1, -1
       soln(k) = soln(k) - gamma(k+1)*soln(k+1)
    enddo

  end subroutine tridiag_solver

!****************************************************************************

end module cism_sparse_pcg

!****************************************************************************
