!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module file_functions

  implicit none

  public :: file_read_1d, file_read_2d

  private ! Default Scope

  contains

!===============================================================================
  subroutine file_read_1d( file_unit, path_and_filename,  & 
                           num_datapts, entries_per_line, &
                           variable )

!     Description:
!     This subroutine reads in values from a data file with a number of
!     rows and a declared number of columns (entries_per_line) of data.
!     It reads in the data in the form of:
!     1 ==> (row 1, column 1); 2 ==> (row 1, column 2); etc.
!
!     Example:  a diagram of a data file with 18 total data points
!               (DP1 to DP18), with 4 data points per row.
!
!              i = 1     i = 2     i = 3     i = 4
!           ---------------------------------------
!     k = 1 |   DP1       DP2       DP3       DP4
!           |
!     k = 2 |   DP5       DP6       DP7       DP8
!           |
!     k = 3 |   DP9       DP10      DP11      DP12
!           |
!     k = 4 |   DP13      DP14      DP15      DP16
!           |
!     k = 5 |   DP17      DP18
!
!     See Michael Falk's comments below for more information.
!-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: fstderr ! Constant(s)

    implicit none

    integer, intent(in) :: & 
     file_unit,          & ! Unit number of file being read.
     num_datapts,        & ! Total number of data points being read in.
     entries_per_line   ! Number of data points
    ! on one line of the file being read.

    character(*), intent(in) :: & 
     path_and_filename  ! Path to file and filename of file being read.

    real( kind = core_rknd ), dimension(num_datapts), intent(out) :: & 
     variable           ! Data values output into variable

    integer :: k        ! Data file row number.
    integer :: i        ! Data file column number.
    integer :: ierr

    ! ---- Begin Code ----
! A ThreadLock is necessary here because FORTRAN can only have each file open on
! one file_unit at a time. For example, suppose we are running CLUBB in parallel
! with OpenMP using two threads. Suppose the first thread opens the file with file_unit = 0
! (file_unit is assigned a value based on thread number).
! Then suppose, that before thread 1 exits, thread 2 opens the same file with file_unit = 1.
! This would cause FORTRAN to crash.
!$omp critical

    ! Open data file.
    open( unit=file_unit, file=path_and_filename, action='read', status='old', &
          iostat=ierr )
    if ( ierr /= 0 ) then
      write(fstderr,*) "CLUBB encountered an error trying to open "//path_and_filename
      error stop "Error opening forcings file"
    end if

    ! Michael Falk wrote this routine to read data files in a particular format for mpace_a.
    ! Each line has a specific number of values, until the last line in the file, which
    ! has the last few values and then ends.  This reads the correct number of values on
    ! each line.  24 September 2007

    ! Loop over each full line of the input file.
    do k = 1, (num_datapts/entries_per_line), 1
      read(file_unit,*) ( variable( ((k-1)*entries_per_line) + i ), & 
                                    i=1,entries_per_line )
    enddo
    ! Read any partial line remaining.
    if ( mod(num_datapts,entries_per_line) /= 0 ) then
      k = (num_datapts/entries_per_line)
      read(file_unit,*) ( variable( (k*entries_per_line) + i ), & 
                          i=1,(mod(num_datapts,entries_per_line)) )
    endif

    ! Close data file.
    close( file_unit )

!$omp end critical

    return

  end subroutine file_read_1d

!===============================================================================
  subroutine file_read_2d( device, file_path, file_dimension1, & 
                           file_dimension2, file_per_line, variable )

! Description:
!   Michael Falk wrote this routine to read data files in a particular format for mpace_a.
!   The 2d mpace_a files list the (file_dimension2) values on a given vertical level, then
!   moves to the next level to list its values.
!   Each line has a specific number of values, until the last line on a level, which
!   is short-- it has the last few values and then a line break.  The next line, beginning
!   the next level, is full-sized again.  24 September 2007
!
! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    integer, intent(in) :: & 
      device, & 
      file_dimension1, & 
      file_dimension2, & 
      file_per_line

    character(*), intent(in) :: & 
      file_path

    real( kind = core_rknd ), dimension(file_dimension1,file_dimension2), intent(out) :: & 
      variable

    integer i, j, k

    ! ---- Begin Code ----

    variable = -999._core_rknd ! Initialize to nonsense values

    open(device,file=file_path,action='read')

    do k=1,(file_dimension1)                ! For each level in the data file,
      do j=0,((file_dimension2/file_per_line)-1)
        read(device,*) (variable(k,(j*file_per_line)+i), & ! read file_per_line values in,
            i=1,file_per_line)
      end do
      read (device,*) (variable(k,(j*file_per_line)+i),           & ! then read the partial line
            i=1,(mod(file_dimension2,file_per_line)))
    end do                                              ! and then start over at the next level.

    close(device)

    return
  end subroutine file_read_2d

!===============================================================================

end module file_functions
