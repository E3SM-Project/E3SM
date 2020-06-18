!----------------------------------------------------------------------------
! $Id: gmres_wrap.F90 7012 2014-07-07 14:18:31Z schemena@uwm.edu $
!==============================================================================

module gmres_wrap

#ifdef MKL

  ! Description:
  ! This module wraps the MKL version of GMRES, an iterative solver. Note that
  ! this will only work for the MKL-specific version of GMRES--any other GMRES
  ! implementations will require retooling of this code!
  !
  ! The primary subroutine, gmres_solve utilizes GMRES to solve a given matrix.
  !
  ! There is also a gmres_init, which initializes some of the internal data
  ! used for the wrapper.
  !
  ! This wrapper automatically keeps prior solutions to use the previous data
  ! to speed up the solves. For the purposes of allowing this solver to be used
  ! with more than one matrix type, the wrapper has a "solve index" variable.
  ! Pass in the proper solve index variable to associate your solve with
  ! previous solves of the same matrix.

  use gmres_cache, only: &
    maximum_gmres_idx ! Variable

  implicit none

  public :: gmres_solve, gmres_init

  private ! Default scope

  contains

  subroutine gmres_init(max_numeqns, max_elements) ! Intent(in)

    ! Description:
    ! Initialization subroutine for the GMRES iterative matrix equation solver
    !
    ! This subroutine initializes the previous memory handles for the GMRES
    ! routines, for the purpose of speeding up calculations.
    ! These handles are initialized to a size specified by the number of
    ! equations specified in this subroutine.
    !
    ! WARNING: Once initialized, only use the specified gmres_idx for that
    ! particular matrix! Failure to do so could result in greatly decreased
    ! performance, incorrect solutions, or both!
    !
    ! Once this is called, the proper prev_soln_<matrix> and prev_lu_<matrix>
    ! handles in the gmres_cache module can be used, and will need to be passed
    ! in to gmres_solve for that matrix.
    !
    ! References:
    !   None

    use gmres_cache, only: &
      gmres_cache_matrix_init ! Subroutines

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      max_numeqns, & ! Maximum number of equations for a matrix that will be
                     ! solved with GMRES
      max_elements   ! Maximum number of non-zero elements for a matrix that
                     ! will be solved with GMRES

    call gmres_cache_matrix_init( max_numeqns, max_elements, maximum_gmres_idx )

  end subroutine gmres_init

  subroutine gmres_solve(elements, numeqns, &                 !Intent(in)
                         csr_a, csr_ia, csr_ja, tempsize, &   !Intent(in)
                         prev_soln, prev_lu, rhs, temp, &     !Intent(in/out)
                         solution, err_code)                  !Intent(out)

    ! Description:
    ! Solves a matrix equation using GMRES. On the first timestep and every
    ! fifth timestep afterward, a preconditioner is computed for the matrix
    ! and stored. In addition, on the first timestep the matrix is solved using
    ! LAPACK, which is used as the estimate for GMRES for the first timestep.
    ! After this, the previous solution found is used as the estimate.
    !
    ! To use the proper cached preconditioner and solution, make sure you pass
    ! the proper gmres_idx corresponding to the matrix you're solving--using a
    ! value different than what has been used in the past will cause, at best,
    ! a slower solve, and at worst, an incorrect one.
    !
    ! References:
    !   None

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd

    implicit none

    include "mkl_rci.fi"

    ! Input variables
    integer, intent(in) :: &
      elements, &  ! Number of elements in the csr_a/csr_ja arrays
      numeqns      ! Number of equations in the matrix

    real( kind = core_rknd ), dimension(elements), intent(in) :: &
      csr_a        ! A-array description of the matrix in CSR format. This
                   ! will be converted to double precision for the purposes
                   ! of running GMRES.

    integer, dimension(numeqns + 1), intent(in) :: &
      csr_ia       ! IA-array portion of the matrix description in CSR format.
                   ! This describes the indices of the JA-array that start
                   ! new rows. For more details, check the documentation in
                   ! the csr_matrix_module.

    integer, dimension(elements), intent(in) :: &
      csr_ja       ! JA-array portion of the matrix description in CSR format.
                   ! This describes which columns of a are nonzero. For more
                   ! details, check the documentation in the csr_matrix_module.

    integer, intent(in) :: &
      tempsize     ! Denotes the size of the temporary array used for GMRES
                   ! calculations.

    ! Input/Output variables
    real( kind = core_rknd ), dimension(numeqns), intent(inout) :: &
      rhs          ! Right-hand-side vectors to solve the equation for.

    real( kind = dp ), dimension(numeqns), intent(inout) :: &
      prev_soln    ! Previous solution cache vector for the matrix to be solved
                   ! for--pass the proper handle from the gmres_cache module

    real( kind = dp ), dimension(elements), intent(inout) :: &
      prev_lu      ! Previous LU-decomposition a-array for the matrix to be
                   ! solved for--pass the proper handle from the gmres_cache
                   ! module

    real( kind = dp ), dimension(tempsize), intent(inout) :: &
      temp         ! Temporary array that stores working values while the GMRES
                   ! solver iterates

    ! Output variables
    real( kind = core_rknd ), dimension(numeqns), intent(out) :: &
      solution     ! Solution vector, output of solver routine

    integer, intent(out) :: &
      err_code     ! Error code, nonzero if errors occurred.

    ! Local variables
    logical :: l_gmres_run ! Variable denoting if we need to loop and run
                           ! a GMRES iteration again.

    integer :: &
      rci_req, &   ! RCI_Request for GMRES--allows us to take action based
                   ! on what the iterative solver requests to be done.
      iters        ! Total number of iterations GMRES has run.

    integer, dimension(128) :: &
      ipar         ! Parameter array for the GMRES iterative solver

    real( kind = dp ), dimension(128) :: &
      dpar         ! Parameter array for the GMRES iterative solver

    ! The following local variables are double-precision so we can use GMRES
    ! as there is only double-precision support for GMRES. 
    ! We will need to convert our single-precision numbers to double precision
    ! for the duration of the calculations.
    real( kind = dp ), dimension(elements) :: &
      csr_dbl_a    ! Double-precision version of the CSR-format A array

    real( kind = dp ), dimension(numeqns) :: &
      dbl_rhs, &   ! Double-precision version of the rhs vector
      dbl_soln, &  ! Double-precision version of the solution vector
      tempvec      ! Temporary vector for applying inverse LU-decomp matrix
      !tmp_rhs

    ! Variables used to solve the preconditioner the first time with PARDISO.
    !integer, parameter :: &
    !pardiso_size_arrays = 64, &
    !real_nonsymm = 11

    !integer(kind=8), dimension(pardiso_size_arrays) :: &
    !  pt ! PARDISO internal pointer array

    !integer(kind=4), dimension(pardiso_size_arrays) :: &
    !  iparm

    !integer(kind=4), dimension(numeqns) :: &
    !  perm

    ! integer :: i, j

    ! We want to be running, initially.
    l_gmres_run = .true.

    ! Set the default error code to 0 (no errors)
    ! This is to make the default explicit; Fortran initializes
    ! values to 0.
    err_code = 0

    ! Convert our A array and rhs vector to double precision...
    csr_dbl_a = real(csr_a, kind=dp)
    dbl_rhs = real(rhs, kind=dp)

    ! DEBUG: Set our a_array so it represents the identity matrix, and
    ! set the RHS so we can get a meaningful answer.
!    csr_dbl_a = 1_dp
!    csr_dbl_a(1) = 1D1
!    csr_dbl_a(5) = 1D1
!    csr_dbl_a(elements) = 1D1
!    csr_dbl_a(elements - 4) = 1D1
!    do i=10,elements - 9,5
!      csr_dbl_a(i) = 1D1
!    end do
!    do i=1,numeqns,1
!      dbl_rhs(i) = i * 1_dp
!    end do
!    dbl_rhs = 9D3
!    dbl_rhs = 1D1

    ! DEBUG: Make sure our a_array isn't wrong
!    do i=1,elements,1
!      print *, "csr_dbl_a idx",i,"=",csr_dbl_a(i)
!    end do

    ! Figure out the default value for ipar(15) and put it in our ipar_15 int.
    !ip_15 = min(150, numeqns)

    ! Figure out the size of the temp array.
    !tempsize = ((((2*numeqns + 1)*numeqns)+(numeqns*(numeqns+9))/2) + 1)
      ! This ugly equation was lifted from the Intel documentation of dfgmres:
      ! http://www.intel.com/software/products/mkl/docs/webhelp/ssr/functn_rci_dfgmres.html
      ! All of the ipar(15)s have been replaced with "numeqns", as the code
      ! examples seemed to use N (numeqns) in place of ipar(15).

    ! Allocate the temp array.
    !allocate(temp(1:tempsize))

    ! Generate our preconditioner matrix with the ILU0 subroutine.
    call dcsrilu0( numeqns, csr_dbl_a, csr_ia, csr_ja, &
                     prev_lu, ipar, dpar, err_code )

    ! On the first timestep we need to solve our preconditioner to give us
    ! our first solution estimate. After this, the previous solution will
    ! suffice as an estimate.
!    if (iteration_num(gmres_idx) == 0) then
      !solve with precond_a, csr_ia, csr_ja.
      !One thing to test, too: try just setting the solution vector to 1
      ! for the first timestep and see if it's not too unreasonably slow?
!      call pardisoinit( pt, real_nonsymm, iparm )
#ifdef _OPENMP
!      iparm(3) = omp_get_max_threads()
#else
!      iparm(3) = 1
#endif

!      call pardiso( pt, 1, 1, real_nonsymm, 13, numeqns,        & !Intent(in)
!                    prev_lu, csr_ia, csr_ja, perm, 1, iparm, 0, & !Intent(in)
!                    dbl_rhs,                                    & !Intent(inout)
!                    prev_soln, err_code )                         !Intent(out)
!    end if !iteration_num == 1

    !DEBUG: Set apporximate solution vector to 0.9 (?) for now
    !prev_soln(:) = 0.9_dp

    !do i=1,numeqns,1
    !  print *, "Current approximate solution idx",i,"=",prev_soln(i)
    !end do

    ! Initialize our solution vector to the previous solution passed in
    dbl_soln = prev_soln

    ! Set up the GMRES solver.
    call dfgmres_init( numeqns, dbl_soln, dbl_rhs, &
                       rci_req, ipar, dpar, temp ) 

    ! Set the parameters that tell GMRES to handle stopping tests
    ipar(9) = 1
    ipar(10) = 0
    ipar(12) = 1

    ! Set the parameter that tells GMRES to use a preconditioner
    ipar(11) = 1

    ! Check our GMRES settings.
    call dfgmres_check( numeqns, dbl_soln, dbl_rhs, &
                        rci_req, ipar, dpar, temp )

    ! Start the GMRES solver. We set up a while loop which will be broken when
    ! the GMRES solver indicates that a solution has been found.
    do while(l_gmres_run)
      !print *, "********************************************************"
      !print *, "BEGINNING ANOTHER ITERATION..."
      !print *, "========================================================"
      ! Run a GMRES iteration.
      call dfgmres( numeqns, dbl_soln, dbl_rhs, &
                    rci_req, ipar, dpar, temp )

      select case(rci_req)
        case (0)
          l_gmres_run = .false.
        case (1)
          ! Multiply our left-hand side by the vector placed in the temp array,
          ! at ipar(22), and place the result in the temp array at ipar(23).
          ! Display temp(ipar(22))
          ! print *, "------------------------------------------------"
          ! print *, "RCI_REQ=1: MULTIPLY VECTOR BY A MATRIX"
          ! do i=1,numeqns,1
          !   print *, "Tempvec before, idx",i,"=",temp(ipar(22)+i-1)
          ! end do
          call mkl_dcsrgemv( 'N', numeqns, csr_dbl_a, csr_ia, csr_ja, &
                             temp(ipar(22)), temp(ipar(23)) ) ! Known magic number
          ! do i=1,numeqns,1
          !   print *, "Tempvec after, idx",i,"=",temp(ipar(23)+i-1)
          ! end do
          ! print *, "------------------------------------------------"
        case (2)
          ! Ignore this for now, see if GMRES ever escapes.
        case (3)
          ! Apply the inverse of the preconditioner to the vector placed in the
          ! temp array at ipar(22), and place the result in the temp array at
          ! ipar(23).
          !print *, "------------------------------------------------"
          !print *, "RCI_REQ=3: APPLY PRECONDITION TO VECTOR"
          !do i=1,numeqns,1
          !  print *, "Tempvec before, idx",i,"=",temp(ipar(22)+i-1)
          !end do
          call mkl_dcsrtrsv( 'L', 'N', 'U', numeqns, &
                             prev_lu, csr_ia, csr_ja, &
                             temp(ipar(22)), tempvec ) ! Known magic number
          call mkl_dcsrtrsv( 'U', 'N', 'N', numeqns, &
                             prev_lu, csr_ia, csr_ja, &
                             tempvec, temp(ipar(23)) ) ! Known magic number
          !do i=1,numeqns,1
          !  print *, "Tempvec after, idx",i,"=",temp(ipar(23)+i-1)
          !end do
          !print *, "------------------------------------------------"

        case (4)
!          if (dpar(7) < GMRES_TOL) then
!            l_gmres_run = .false.
!          else
!            ! Keep running, we aren't there yet.
!            l_gmres_run = .true.
!          end if
        case default
          ! We got a response we weren't expecting. This is probably bad.
          ! (Then again, maybe it's just not something we accounted for?)
          ! Regardless, let's set an error code and break out of here.
          print *, "Unknown rci_request returned from GMRES:", rci_req
          l_gmres_run = .false.
          err_code = -1
      end select
      ! Report current iteration
!      call dfgmres_get( numeqns, dbl_soln, dbl_rhs, rci_req, &
!                        ipar, dpar, temp, iters )
!      print *, "========================================================"
!      print *, "END OF LOOP: REPORTING INFORMATION"
!      print *, "Current number of GMRES iterations: ", iters
!      do i=1,numeqns,1
!        print *, "double value of soln so far, idx",i,"=",dbl_soln(i)
!      end do
!      print *, "========================================================"
!      print *, "********************************************************"
    end do
    !if (err_code == 0) then

      ! Get the answer, convert it to single-precision
      call dfgmres_get( numeqns, dbl_soln, dbl_rhs, rci_req, &
                        ipar, dpar, temp, iters )

      !print *, "Total iterations for GMRES:",iters

      !do i=1,numeqns,1
      !  print *, "double value of soln, idx",i,"=",dbl_soln(i)
      !end do

      ! Store our solution as the previous solution for use in the next
      ! simulation timestep.
      prev_soln = dbl_soln

      solution = real(dbl_soln)
    !end if
    
  end subroutine gmres_solve

#endif /* MKL */

end module gmres_wrap
