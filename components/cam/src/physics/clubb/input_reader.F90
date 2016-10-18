!-----------------------------------------------------------------------
!$Id: input_reader.F90 8014 2016-03-12 00:54:18Z raut@uwm.edu $
!===============================================================================
module input_reader

! Description:
!  This module is respondsible for the procedures and structures necessary to
!  read in "SAM-Like" case specific files. Currently only the
!  <casename>_sounding.in file is formatted to be used by this module.
!
! References:
!   None
!---------------------------------------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  private

  public :: one_dim_read_var, &
            read_one_dim_file, &
            two_dim_read_var, &
            read_two_dim_file, &
            fill_blanks_one_dim_vars, &
            fill_blanks_two_dim_vars, &
            deallocate_one_dim_vars, &
            deallocate_two_dim_vars, &
            read_x_table, &
            read_x_profile, &
            get_target_index, &
            count_columns

  ! Derived type for representing a rank 1 variable that has been read in by one
  ! of the procedures.
  type one_dim_read_var

    character(len=30) :: name               ! Name of the variable

    character(len=30) :: dim_name           ! Name of the dimension that the
    !                                         variable varies along

    real( kind = core_rknd ), dimension(:), allocatable :: values   ! Values of that variable

  end type one_dim_read_var

  ! Derived type for representing a rank 2 variable that has been read in by one
  ! of the procedures.
  type two_dim_read_var

    character(len=30) :: name                ! Name of the variable

    character(len=30) :: dim1_name           ! Name of one of the dimensions
    !                                          that the variable varies along

    character(len=30) :: dim2_name           ! Name of the other variable that
    !                                          the variable varies along

    real( kind = core_rknd ), dimension(:,:), allocatable :: values  ! Values of that variable

  end type two_dim_read_var


  ! Constant Parameter(s)
  real( kind = core_rknd ), parameter, private :: &
    blank_value = -999.9_core_rknd ! Used to denote if a value is missing from the file

  contains

  !-------------------------------------------------------------------------------------------------
  subroutine read_two_dim_file( iunit, nCol, filename, read_vars, other_dim )
    !
    ! Description: This subroutine reads from a file containing data that varies
    !              in two dimensions. These are dimensions are typically height
    !              and time.
    !
    !-----------------------------------------------------------------------------------------------
    use constants_clubb, only: &
      fstderr ! Constant(s)

    use input_names, only: &
      time_name ! Constant(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: trim, index

    ! Input Variable(s)

    integer, intent(in) :: iunit ! File I/O unit

    integer, intent(in) :: nCol ! Number of columns expected in the data file


    character(len=*), intent(in) :: filename ! Name of the file being read from

    ! Output Variable(s)
    type (two_dim_read_var), dimension(nCol),intent(out) :: read_vars ! Structured information
    !                                                                   from the file

    type (one_dim_read_var), intent(out) :: other_dim  ! Structured information
    !                                                    on the dimesion not stored in read_vars

    ! Local Variables
    character(len=30),dimension(nCol) :: names ! Names of variables

    integer nRowI ! Inner row

    integer nRowO ! Outer row

    integer :: k, j, i

    logical :: isComment

    character(len=200) :: tmpline

    real( kind = core_rknd ), dimension(nCol) :: tmp

    integer :: input_status ! The status of a read statement

    ! ---- Begin Code ----

    ! First run through, take names and determine how large the data file is.
    open(unit=iunit, file=trim( filename ), status = 'old', action='read' )

    isComment = .true.

    ! Skip all the comments at the top of the file
    do while ( isComment )
      read(iunit,fmt='(A)') tmpline
      k = index( tmpline, "!" )
      isComment = .false.
      if ( k > 0 ) then
        isComment = .true.
      end if
    end do

    ! Go back to the line that wasn't a comment.
    backspace(iunit)

    read(iunit, fmt=*) names

    nRowO = 0
    do while(.true.)
      read(iunit, *, iostat=input_status) tmp(1), nRowI

      ! If input_status shows an end of data, then exit the loop
      if( input_status < 0 ) then
        exit
      else if ( input_status > 0 ) then
        write(fstderr,*) "Error reading data from file: " //trim( filename )
        stop "Fatal error input_reader"
      end if

      if( nRowI < 1 ) then
        stop "Number of elements must be an integer and greater than zero in two-dim  input file."
      end if

      do k =1, nRowI
        read(iunit, *) tmp
      end do
      nRowO = nRowO + 1
    end do

    do i=1, nRowO

      backspace(iunit)

      do j=1, nRowI

        backspace(iunit)

      end do

    end do

    backspace(iunit)

    ! Store the names into the structure and allocate accordingly
    do k =1, nCol
      read_vars(k)%name = names(k)
      read_vars(k)%dim1_name = time_name
      read_vars(k)%dim2_name = names(1)

      allocate( read_vars(k)%values(nRowI, nRowO) )
    end do

    other_dim%name = time_name
    other_dim%dim_name = time_name

    allocate( other_dim%values(nRowO) )

    ! Read in the data again to the newly allocated arrays
    do k=1, nRowO
      read(iunit,*) other_dim%values(k)
      do j=1, nRowI
        read(iunit,*) ( read_vars(i)%values(j,k), i=1, nCol)
      end do
    end do

    close(iunit)

    ! Eliminate a compiler warning
    if ( .false. ) print *, tmp

    return
  end subroutine read_two_dim_file

  !------------------------------------------------------------------------------------------------
  subroutine read_one_dim_file( iunit, nCol, filename, read_vars )
    !
    ! Description: 
    !   This subroutine reads from a file containing data that varies
    !   in one dimension. The dimension is typically time.
    ! 
    ! References:
    !   None
    !----------------------------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External

    intrinsic :: trim, index

    ! Input Variable(s)

    integer, intent(in) :: iunit ! I/O unit

    integer, intent(in) :: nCol ! Number of columns expected in the data file

    character(len=*), intent(in)  :: filename ! Name of the file being read from

    ! Output Variable(s)

    type (one_dim_read_var), dimension(nCol),intent(out) :: &
      read_vars ! Structured information from the file

    ! Local Variable(s)
    character(len=30),dimension(nCol) :: names

    character(len=200) :: tmpline

    integer nRow

    integer :: k, j

    real( kind = core_rknd ), dimension(nCol) :: tmp

    logical :: isComment

    integer :: input_status ! The status of a read statement

    ! ---- Begin Code ----

    isComment = .true.

    ! First run through, take names and determine how large the data file is.
    open(unit=iunit, file=trim( filename ), status = 'old' )

    ! Skip all the comments at the top of the file
    do while(isComment)
      read(iunit,fmt='(A)') tmpline
      k = index( tmpline, "!" )
      isComment = .false.
      if(k > 0) then
        isComment = .true.
      end if
    end do

    ! Go back to the line that wasn't a comment.
    backspace(iunit)

    read(iunit, fmt=*) names

    ! Count up the number of rows
    nRow = 0
    do while(.true.)
      read(iunit, *, iostat=input_status) tmp

      ! If input_status shows an end of file, exit the loop
      if( input_status < 0 ) then
        exit
      end if

      nRow = nRow+1
    end do

    ! Rewind that many rows
    do k = 0, nRow
      backspace(iunit)
    end do

    ! Store the names into the structure and allocate accordingly
    do k = 1, nCol
      read_vars(k)%name = names(k)
      read_vars(k)%dim_name = names(1)
      allocate( read_vars(k)%values(nRow) )
    end do

    ! Read in the data again to the newly allocated arrays
    do k=1, nRow
      read(iunit,*) ( read_vars(j)%values(k), j=1, nCol)
    end do

    close(iunit)

    ! Avoiding compiler warning
    if ( .false. ) print *, tmp

    return

  end subroutine read_one_dim_file

  !------------------------------------------------------------------------------------------------
  subroutine fill_blanks_one_dim_vars( num_vars, one_dim_vars )
    !
    ! Description: 
    !   This subroutine fills in the blank spots (signified by constant blank_value)
    !   with values linearly interpolated using the first element of the array as a
    !   guide.
    !
    ! References:
    !   None
    !----------------------------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: size

    ! Input Variable(s)
    integer, intent(in) :: num_vars ! Number of elements in one_dim_vars

    ! Input/Output Variable(s)
    type(one_dim_read_var), dimension(num_vars), intent(inout) :: &
      one_dim_vars ! Read data that may have gaps.

    ! Local variable(s)
    integer :: i

    ! ---- Begin Code ----

    do i=1, num_vars
      one_dim_vars(i)%values = linear_fill_blanks( size( one_dim_vars(i)%values ), &
                                                   one_dim_vars(1)%values, one_dim_vars(i)%values, &
                                                   0.0_core_rknd )
    end do

    return

  end subroutine fill_blanks_one_dim_vars

  !------------------------------------------------------------------------------------------------
  subroutine fill_blanks_two_dim_vars( num_vars, other_dim, two_dim_vars )
    !
    ! Description: 
    !   This subroutine fills in the blank spots (signified by the
    !   constant blank_value with values linearly interpolated using the first 
    !   element of the array and the values in the other_dim argument as a guide.
    !
    !   This is a two step process. First we assume that the other_dim values
    !   have no holes, but there are blanks for that variable across that
    !   dimension. Then we fill holes across the dimension whose values are first
    !   in the array of two_dim_vars.
    !
    !   Ex. Time is the 'other_dim' and Height in meters is the first element in
    !   two_dim_vars.
    !
    ! References:
    !   None
    !----------------------------------------------------------------------------------------------

    implicit none

    ! External
    intrinsic :: size

    ! Input Variable(s)
    integer, intent(in) :: num_vars ! Number of elements in one_dim_vars

    ! Input/Output Variable(s)
    type(one_dim_read_var), intent(in) :: other_dim ! Read data

    type(two_dim_read_var), dimension(num_vars), intent(inout) ::  &
      two_dim_vars ! Read data that may have gaps.

    ! Local variables
    integer :: i,j ! Loop iterators

    integer :: & 
      dim_size, &    ! 1st dimension size
      other_dim_size ! 2nd dimension size

    ! ---- Begin Code ----

    dim_size = size( two_dim_vars(1)%values, 1 )

    other_dim_size = size( other_dim%values )

    do i=2, num_vars
      ! Interpolate along main dim
      do j=1, other_dim_size
        two_dim_vars(i)%values(:,j) = linear_fill_blanks( dim_size, &
                                        two_dim_vars(1)%values(:,j), &
                                        two_dim_vars(i)%values(:,j), blank_value ) 
      end do ! j = 1 .. other_dim_size

      ! Interpolate along other dim
      do j=1, dim_size
        two_dim_vars(i)%values(j,:) = linear_fill_blanks( other_dim_size, &
                                        other_dim%values, &
                                        two_dim_vars(i)%values(j,:), blank_value ) 
      end do ! j = 1 .. dim_size

    end do ! i = 2 .. num_vars

    return

  end subroutine fill_blanks_two_dim_vars


  !------------------------------------------------------------------------------------------------
  function linear_fill_blanks( dim_grid, grid, var, default_value ) &
    result( var_out )
  !
  ! Description: 
  !   This function fills blanks in array var using the grid
  !   as a guide. Blank values in var are signified by being
  !   less than or equal to the constant blank_value.
  !
  ! References:
  !   None
  !-----------------------------------------------------------------------------------------------

    use interpolation, only: zlinterp_fnc

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: dim_grid ! Size of grid

    real( kind = core_rknd ), dimension(dim_grid), intent(in) :: &
      grid ! Array that var is being interpolated to.

    real( kind = core_rknd ), dimension(dim_grid), intent(in) :: &
      var ! Array that may contain gaps.

    real( kind = core_rknd ), intent(in) :: &
      default_value ! Default value if entire profile == blank_value

    ! Output Variable(s)
    real( kind = core_rknd ), dimension(dim_grid) :: &
      var_out ! Return variable

    ! Local Variables
    real( kind = core_rknd ), dimension(dim_grid) :: temp_grid
    real( kind = core_rknd ), dimension(dim_grid) :: temp_var

    integer :: i
    integer :: amt

    logical :: reversed

    ! ---- Begin Code ----

    reversed = .false.

    ! Essentially this code leverages the previously written zlinterp function.
    ! A smaller temporary grid and var variable are being created to pass to
    ! zlinterp. zlinterp then performs the work of taking the temporary var
    ! array and interpolating it to the actual grid array.

    amt = 0
    do i=1, dim_grid
      if ( var(i) > blank_value ) then
        amt = amt + 1
        temp_var(amt) = var(i)
        temp_grid(amt) = grid(i)
      end if
      if ( i > 1 ) then
        if ( grid(i) < grid(i-1) ) then
          reversed = .true.
        end if
      end if
    end do


    if ( amt == 0 ) then
      var_out = default_value
    else if (amt < dim_grid) then
      if ( reversed ) then
        var_out = zlinterp_fnc( dim_grid, amt, -grid, -temp_grid(1:amt), temp_var(1:amt) )
      else
        var_out = zlinterp_fnc( dim_grid, amt, grid, temp_grid(1:amt), temp_var(1:amt) )
      end if
    else
      var_out = var
    end if

    return
  end function linear_fill_blanks
  !----------------------------------------------------------------------------
  subroutine deallocate_one_dim_vars( num_vars, one_dim_vars )
    !
    !  Description: 
    !    This subroutine deallocates the pointer stored in
    !    one_dim_vars%value for the whole array.
    !
    !------------------------------------------------------------------------------
    implicit none

    ! External functions
    intrinsic :: allocated

    ! Input Variable(s)
    integer, intent(in) :: num_vars ! Number of elements in one_dim_vars

    type(one_dim_read_var), dimension(num_vars), intent(inout) :: &
      one_dim_vars ! Read data that may have gaps.

    ! Local Variable(s)
    integer :: i

    ! Begin Code

    do i=1, num_vars

      if ( allocated( one_dim_vars(i)%values ) ) then

        deallocate( one_dim_vars(i)%values )

      end if

    end do ! 1 .. num_vars

    return
  end subroutine deallocate_one_dim_vars

  !------------------------------------------------------------------------------------------------
  subroutine deallocate_two_dim_vars( num_vars, two_dim_vars, other_dim )
    !
    ! Description: 
    !   This subroutine deallocates the pointer stored in
    !   two_dim_vars%value for the whole array
    !
    ! References:
    !   None
    !----------------------------------------------------------------------------------------------
    implicit none

    ! External Functions
    intrinsic :: allocated

    ! Input Variable(s)
    integer, intent(in) :: num_vars ! Number of elements in one_dim_vars

    ! Input/Output Variables
    type(one_dim_read_var), intent(inout) :: other_dim

    type(two_dim_read_var), dimension(num_vars), intent(inout) :: &
       two_dim_vars ! Read data that may have gaps.

    ! Local Variable(s)
    integer :: i

    ! ---- Begin Code ----

    do i=1, num_vars

      if ( allocated( two_dim_vars(i)%values ) ) then

        deallocate(two_dim_vars(i)%values)

      end if

    end do

    if ( allocated( other_dim%values ) ) then

      deallocate(other_dim%values)

    end if

    return
  end subroutine deallocate_two_dim_vars
  !------------------------------------------------------------------------------------------------
  function read_x_table( nvar, xdim, ydim, target_name, retVars ) result( x )
    !
    ! Description: 
    !   Searches for the variable specified by target_name in the
    !   collection of retVars. If the function finds the variable then it returns
    !   it. If it does not the program using this function will exit gracefully
    !   with a warning message.
    !
    ! References:
    !   None
    !-----------------------------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use constants_clubb, only: &
      fstderr ! Constant(s)

    implicit none

    ! External Functions
    intrinsic :: trim

    ! Input Variable(s)
    integer, intent(in) :: nvar ! Number of variables in retVars

    integer, intent(in) :: xdim, ydim

    character(len=*), intent(in) :: &
      target_name ! Name of the variable that is being searched for

    type(two_dim_read_var), dimension(nvar), intent(in) :: &
      retVars ! Collection of data being searched through

    ! Output Variable(s)
    real( kind = core_rknd ), dimension( xdim, ydim ) :: x

    ! Local Variables
    integer :: i ! Loop iterator

    logical :: l_found

    ! ---- Begin Code ----

    l_found = .false.

    i = 1

    do while( i <= nvar .and. .not. l_found)

      if( retVars(i)%name == target_name ) then

        l_found = .true.

        x = retVars(i)%values

      end if

      i=i+1

    end do ! i <= nvar .and. not l_found

    if ( .not. l_found ) then

      write(fstderr,*) trim( target_name )//" could not be found."

      stop "Fatal error in function read_x_table"

    end if

    return

  end function read_x_table


  !------------------------------------------------------------------------------------------------
  function read_x_profile( nvar, dim_size, target_name, retVars, &
                           input_file ) result( x )
    !
    !  Description: 
    !    Searches for the variable specified by target_name in the
    !    collection of retVars. If the function finds the variable then it returns
    !    it. If it does not the program using this function will exit gracefully
    !    with a warning message.
    !
    ! Modified by Cavyn, June 2010
    !----------------------------------------------------------------------------------------------

    use constants_clubb, only: &
      fstderr  ! Variable for writing to error stream

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External Functions
    intrinsic :: present, size, trim

    ! Input Variable(s)
    integer, intent(in) :: & 
      nvar,  & ! Number of variables in retVars
      dim_size ! Size of the array returned

    character(len=*), intent(in) :: &
      target_name  ! Name of the variable that is being searched for

    type(one_dim_read_var), dimension(nvar), intent(in) :: &
      retVars ! Collection being searched

    character(len=*), optional, intent(in) :: &
      input_file ! Name of the input file containing the variables

    ! Output Variable(s)
    real( kind = core_rknd ), dimension(dim_size) :: x

    ! Local Variables
    integer :: i

    ! ---- Begin Code ----    

    i = get_target_index( nvar, target_name, retVars )
    
    if ( i > 0 ) then
        x(1:size(retVars(i)%values)) = retVars(i)%values
        
    else
      if( present( input_file ) ) then
        write(fstderr,*) trim( target_name ), ' could not be found. Check the file ', input_file
      else
        write(fstderr,*) trim( target_name ), ' could not be found. Check your sounding.in file.'
      end if ! present( input_file )
      stop "Fatal error in read_x_profile"
      
    end if ! target_exists_in_array

    return

  end function read_x_profile
  
  !------------------------------------------------------------------------------
  function get_target_index( nvar, target_name, retVars) result( i )
    !
    ! Description:
    !   Returns the index of the variable specified by target_name in the
    !   collection of retVars. Returns -1 if variable does not exist in retVars
    !
    ! References:
    !   None
    !
    ! Created by Cavyn, July 2010
    !----------------------------------------------------------------------------------------------

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: nvar                    ! Number of variables in retVars
    character(len=*), intent(in) :: target_name           ! Variable being searched for
    type(one_dim_read_var), dimension(nvar), intent(in) :: retVars ! Collection being searched

    ! Output Variable
    integer :: i

    ! Local Variable(s)
    logical :: l_found
    
    !----------------BEGIN CODE------------------
    
    l_found = .false.
    
    i = 0
    do while ( i < nvar .and. .not. l_found )
      i = i+1
      if( retVars(i)%name == target_name ) then
        l_found = .true.
      end if
    end do
    
    if( .not. l_found ) then
      i = -1
    end if

    return

  end function get_target_index

  !=============================================================================
  function count_columns( iunit, filename ) result( nCols )
  ! Description:
  !   This function counts the number of columns in a file, assuming that the
  !   first line of the file contains only column headers. (Comments are OK)

  ! References:
  !   None

  ! Created by Cavyn, July 2010
  !-----------------------------------------------------------------------------

    implicit none

    ! External
    intrinsic :: index, trim, size

    ! Input Variables
    integer, intent(in) :: iunit ! I/O unit
    character(len=*), intent(in) :: filename ! Name of the file being read from
    
    ! Output Variable
    integer :: nCols ! The number of data columns in the selected file
    
    ! Local Variables
    integer :: i, k                               ! Loop Counter
    character(len=200) :: tmp                     ! Temporary char buffer
    character(len=200), dimension(50) :: colArray ! Max of 50 columns
    logical :: isComment
    integer :: status_var ! IO status for read statement


    ! -------------------------BEGIN CODE-------------------------------------
    
    isComment = .true.

    open(unit=iunit, file=trim( filename ), status = 'old' )

    ! Skip all the comments at the top of the file
    do while(isComment)
      read(iunit,fmt='(A)') tmp
      k = index( tmp, "!" )
      isComment = .false.
      if(k > 0) then
        isComment = .true.
      end if
    end do

    ! Go back to the line that wasn't a comment.
    backspace(iunit)
    
    ! Count the number of columns
    nCols = 0
    colArray = ""
    read(iunit,fmt='(A)',iostat=status_var) tmp
    ! Only continue if there was no IO error or end of data
    if( status_var == 0 ) then
      ! Move all words into an array
      read(tmp,*,iostat=status_var) (colArray(i), i=1,size( colArray )) 

    else if ( status_var > 0 ) then
      ! Handle the case where we have an error before the EOF marker is found
      stop "Fatal error reading data in time_dependent_input function count_columns"

    end if
    
    do i=1,size( colArray )
      if( colArray(i) /= "" ) then ! Increment number of columns until array is blank
        nCols = nCols+1
      end if
    end do
    
    close(iunit)

  end function count_columns

!------------------------------------------------------------------------------
end module input_reader
