!$Id: corr_matrix_module.F90 5997 2012-12-18 20:47:09Z raut@uwm.edu $
!---------------------------------------------------------------------------------------------------
module corr_matrix_module

  implicit none

  ! Latin hypercube indices / Correlation array indices
  integer, public :: &
    iiLH_s_mellor = -1, &
    iiLH_t_mellor = -1, &
    iiLH_w        = -1
!$omp threadprivate(iiLH_s_mellor, iiLH_t_mellor, iiLH_w)

  integer, public :: &
   iiLH_rrain    = -1, &
   iiLH_rsnow    = -1, &
   iiLH_rice     = -1, &
   iiLH_rgraupel = -1
!$omp threadprivate(iiLH_rrain, iiLH_rsnow, iiLH_rice, iiLH_rgraupel)

  integer, public :: &
   iiLH_Nr       = -1, &
   iiLH_Nsnow    = -1, &
   iiLH_Ni       = -1, &
   iiLH_Ngraupel = -1, &
   iiLH_Nc       = -1
!$omp threadprivate(iiLH_Nr, iiLH_Nsnow, iiLH_Ni, iiLH_Ngraupel, iiLH_Nc)

  public :: read_correlation_matrix

  private :: get_corr_var_index

  private

  contains

  !-----------------------------------------------------------------------------
  subroutine read_correlation_matrix( iunit, input_file, d_variables, &
                                      corr_array )

    ! Description:
    !   Reads a correlation variance array from a file and stores it in an array.
    !-----------------------------------------------------------------------------

    use input_reader, only: &
      one_dim_read_var, & ! Variable(s)
      read_one_dim_file, deallocate_one_dim_vars, count_columns ! Procedure(s)

    use matrix_operations, only: set_lower_triangular_matrix ! Procedure(s)

    use constants_clubb, only: fstderr ! Variable(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      iunit, &    ! File I/O unit
      d_variables ! number of variables in the array

    character(len=*), intent(in) :: input_file ! Path to the file

    ! Input/Output Variable(s)
    real( kind = core_rknd ), dimension(d_variables,d_variables), intent(inout) :: &
      corr_array ! Correlation variance array

    ! Local Variable(s)

    type(one_dim_read_var), allocatable, dimension(:) :: &
      retVars ! stores the variables read in from the corr_varnce.in file

    integer ::   &
      var_index1,    & ! variable index
      var_index2,    & ! variable index
      nCols,         & ! the number of columns in the file
      i, j         ! Loop index


    !--------------------------- BEGIN CODE -------------------------

    nCols = count_columns( iunit, input_file )

    ! Allocate all arrays based on d_variables
    allocate( retVars(1:nCols) )

    ! Initializing to zero means that correlations we don't have
    ! (e.g. Nc and any variable other than s_mellor ) are assumed to be 0.
    corr_array(:,:) = 0.0_core_rknd

    ! Set main diagonal to 1
    do i=1, d_variables
      corr_array(i,i) = 1.0_core_rknd
    end do

    ! Read the values from the specified file
    call read_one_dim_file( iunit, nCols, input_file, retVars )

    if( size( retVars(1)%values ) /= nCols ) then
      write(fstderr, *) "Correlation matrix must have an equal number of rows and cols in file ", &
            input_file
      stop "Bad data in correlation file."
    end if

    ! Start at 2 because the first index is always just 1.0 in the first row
    ! and the rest of the rows are ignored
    do i=2, nCols
      var_index1 = get_corr_var_index( retVars(i)%name )
      if( var_index1 > -1 ) then
        do j=1, (i-1)
          var_index2 = get_corr_var_index( retVars(j)%name )
          if( var_index2 > -1 ) then
            call set_lower_triangular_matrix &
                 ( d_variables, var_index1, var_index2, retVars(i)%values(j), &
                   corr_array )
          end if
        end do
      end if
    end do

    call deallocate_one_dim_vars( nCols, retVars )

    return
  end subroutine read_correlation_matrix

  !--------------------------------------------------------------------------
  function get_corr_var_index( var_name ) result( i )

    ! Definition:
    !   Returns the index for a variable based on its name.
    !--------------------------------------------------------------------------

    implicit none

    character(len=*), intent(in) :: var_name ! The name of the variable

    ! Output variable
    integer :: i

    !------------------ BEGIN CODE -----------------------------
    i = -1

    select case( trim(var_name) )

    case( "s" )
      i = iiLH_s_mellor

    case( "t" )
      i = iiLH_t_mellor

    case( "w" )
      i = iiLH_w

    case( "Nc" )
      i = iiLH_Nc

    case( "rrain" )
      i = iiLH_rrain

    case( "Nr" )
      i = iiLH_Nr

    case( "rice" )
      i = iiLH_rice

    case( "Ni" )
      i = iiLH_Ni

    case( "rsnow" )
      i = iiLH_rsnow

    case( "Nsnow" )
      i = iiLH_Nsnow

    end select

    return

  end function get_corr_var_index
end module corr_matrix_module
