!----------------------------------------------------------------------------
! $Id: gmres_cache.F90 8014 2016-03-12 00:54:18Z raut@uwm.edu $
!==============================================================================
module gmres_cache

#ifdef MKL

  use clubb_precision, only: &
    dp ! double precision

  ! Description:
  ! This module contains cache data structures for the GMRES wrapper class.
  !
  ! This is mostly to allow us to get around some...odd errors when it was
  ! integrated into the gmres_wrap module. The cache variables are public, as
  ! they will need to be passed in whenever gmres_solve is called.

  implicit none

  public :: gmres_cache_matrix_init, gmres_cache_soln, &
            gmres_cache_temp_init

  private ! Default scope

  real( kind = dp ), public, allocatable, dimension(:,:) :: &
    gmres_prev_soln, &    ! Stores the previous solution vectors from earlier
                          ! GMRES solve runs. The first dimension is for the
                          ! actual vector; the second dimension is to determine
                          ! which cache to access--this is done via the GMRES
                          ! indices for each of the different matrices.
    gmres_prev_precond_a  ! Stores the previous preconditioner matrix from
                          ! earlier GMRES solve runs. The first dimension is
                          ! for the a-array itself; the second dimension is to
                          ! determine which cached array to access--this is
                          ! done via the GMRES indices for each of the
                          ! different matrices.

!$omp threadprivate( gmres_prev_soln, gmres_prev_precond_a )

  real( kind = dp ), public, allocatable, dimension(:) :: &
    gmres_temp_intlc, &   ! Temporary array that stores GMRES internal values
                          ! for the interlaced matrices (2 x gr%nz grid
                          ! levels)
    gmres_temp_norm       ! Temporary array that stores GMRES internal values
                          ! for the non-interlaced matrices (gr%nz grid
                          ! levels)

!$omp threadprivate( gmres_tmp_intlc, gmres_temp_norm )

  integer, public :: &
    gmres_tempsize_norm, &     ! Size of the temporary array for
                               ! non-interlaced matrices
    gmres_tempsize_intlc       ! Size of the temporary array for
                               ! interlaced matrices

!$omp threadprivate( gmres_tempsize_norm, gmres_tempsize_intlc )

  integer, public, parameter :: &
    maximum_gmres_idx = 1 ! Maximum number of different types of solves the
                          ! wrapper can keep memory for. If new matrices are
                          ! added that GMRES is to be used for, increase this
                          ! number and add a public parameter corresponding to
                          ! the matrix below:

  integer, public, parameter :: &
    gmres_idx_wp2wp3 = 1    ! GMRES wrapper index for the wp2_wp3 matrices

  logical, public, dimension(maximum_gmres_idx) :: &
    l_gmres_soln_ok ! Stores if the current solution is "okay"--that is, if an
                    ! initial solution has been passed in for that particular
                    ! cache index. This defaults to false and is set to true
                    ! when a solution is updated.

!$omp threadprivate(l_gmres_soln_ok)

  contains

  subroutine gmres_cache_temp_init(numeqns) ! Intent(in)
    ! Description:
    ! Initialization subroutine for the temporary arrays for GMRES
    !
    ! This subroutine initializes the temporary arrays that are used to work
    ! the GMRES solver.
    !
    ! These temporary arrays are used for all GMRES solves.
    !
    ! References:
    !   None

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      numeqns    ! Number of equations for non-interlaced matrices (gr%nz)

    integer :: &
      numeqns_intlc     ! Number of equations for interlaced matrices

    numeqns_intlc = numeqns * 2

    ! Figure out the sizes of the temporary arrays
    ! The equations were lifted from the Intel documentation of dfgmres:
    ! http://www.intel.com/software/products/mkl/docs/webhelp/ssr/functn_rci_dfgmres.html
    ! All of the ipar(15)s have been replaced with "numeqns", as the code
    ! examples seemed to use N (numeqns) in place of ipar(15).
    gmres_tempsize_norm = ((((2*numeqns + 1)*numeqns) &
                          + (numeqns*(numeqns+9))/2) + 1) ! Known magic number

    gmres_tempsize_intlc = ((((2*numeqns_intlc + 1)*numeqns_intlc) &
                       + (numeqns_intlc*(numeqns_intlc+9))/2) + 1) ! Known magic number

    ! Allocate the temporary arrays
    allocate( gmres_temp_intlc(1:gmres_tempsize_intlc), &
              gmres_temp_norm(1:gmres_tempsize_norm) )

  end subroutine gmres_cache_temp_init

  subroutine gmres_cache_matrix_init(max_numeqns, max_elements, &  ! Intent(in)
                                     max_gmres_idx)                ! Intent(in)
    ! Description:
    ! Initialization subroutine for the caches for GMRES.
    !
    ! This initializes the cache that stores the previous solution and
    ! previous preconditioner values for all GMRES solves.
    !
    ! References:
    !   None

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      max_numeqns, &  ! Maximum number of equations for a matrix that will be
                      ! solved with GMRES
      max_elements, & ! Maximum number of non-zero elements for a matrix that
                      ! will be solved with GMRES
      max_gmres_idx   ! Maximum number of distinct matrices that will be solved
                      ! with GMRES

    allocate( gmres_prev_soln(1:max_numeqns,1:max_gmres_idx), &
              gmres_prev_precond_a(1:max_elements,1:max_gmres_idx) )

    l_gmres_soln_ok = .false.

  end subroutine gmres_cache_matrix_init

  subroutine gmres_cache_soln(numeqns, gmres_idx, solution) ! Intent(in)
    ! Description:
    ! Subroutine that caches a previous solution for a particular GMRES-solved
    ! matrix.
    !
    ! Stores the current solution in the cache so it can be referenced for
    ! the next GMRES solve. This subroutine will also set the solution_ok
    ! flag for that particular GMRES index.
    !
    ! References:
    !  None

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    integer, intent(in) :: &
      numeqns, &  ! The number of equations in the solution vector
      gmres_idx   ! The index for the particular matrix solved by GMRES

    real( kind = core_rknd ), dimension(numeqns), intent(in) :: &
      solution    ! The solution vector to be cached

    gmres_prev_soln(1:numeqns,gmres_idx) = solution

    l_gmres_soln_ok(gmres_idx) = .true.

  end subroutine gmres_cache_soln

#endif /* MKL */

end module gmres_cache
