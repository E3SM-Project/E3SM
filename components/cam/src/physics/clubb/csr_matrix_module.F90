!-----------------------------------------------------------------------
! $Id: csr_matrix_module.F90 8014 2016-03-12 00:54:18Z raut@uwm.edu $
!===============================================================================
module csr_matrix_module

  ! Description:
  ! This module contains some of the matrix description arrays required by
  ! PARDISO, GMRES, and other sparse matrix solvers. The format is called CSR
  ! (compressed sparse row) format, and is currently leveraged through PARDISO
  ! and GMRES.
  ! These are all 1 dimensional arrays that describe a matrix that
  ! will be passed to the solver. The _ja arrays describe which
  ! columns in the matrix have nonzero values--for our purposes, all the
  ! elements on the appropriate diagonals have values. The _ia arrays describe
  ! which _ja array elements correspond to new rows.
  ! Further description of this format can be found in the PARDISO manual, or
  ! alternately, in Intel MKL's documentation.
  ! For our purposes, the _ia and _ja arrays will be fixed for the types
  ! of matrices we have, so we calculate these initially using
  ! initialize_csr_class and simply use the pointers, similar to how
  ! the grid pointers are initialized. This should save a fair amount of time,
  ! as we do not have to recalculate the arrays.
  !
  ! A description of the CSR matrix format:
  ! The CSR matrix format requires three arrays--an a array, 
  ! a ja array, and an ia array.
  !
  ! The a array stores, in sequential order, the actual values in the matrix.
  ! Essentially, just copy the matrix into a 1-dimensional array as you move
  ! from left to right, top down through the matrix. The a array changes
  ! frequently for our purposes in CLUBB, and is not useful to be initialized
  ! here.
  !
  ! The ja array stores, in sequential order, the columns of each element in
  ! the matrix that is nonzero. Essentially, you take the column of each
  ! element that is nonzero as you move from left to right, top down through
  ! the matrix.
  !
  ! An example follows to illustrate the point:
  ! [3.0 2.0 0.0 0.0 0.0 0.0
  !  2.5 1.7 3.6 0.0 0.0 0.0
  !  0.0 5.2 1.7 3.6 0.0 0.0
  !  0.0 0.0 4.7 2.9 0.6 0.0
  !  0.0 0.0 0.0 8.9 4.6 1.2
  !  0.0 0.0 0.0 0.0 5.8 3.7]
  !
  ! Our ja array would look like the following--a pipe denotes a new row:
  ! [1 2 | 1 2 3 | 2 3 4 | 3 4 5 | 4 5 6 | 5 6]
  !
  ! The ia array stores the indices of the ja array that correspond to new rows
  ! in the matrix, with a final entry just beyond the end of the ja matrix
  ! that signifies the end of the matrix.
  ! In our example, the ia array would look like this:
  !
  ! [1 3 6 9 12 15 17]
  !
  ! Similar principles can be applied to find the ia and ja matrices for all
  ! of the general cases CLUBB uses. In addition, because CLUBB typically
  ! uses similar matrices for its calculations, we can simply initialize
  ! the ia and ja matrices in this module rather than repeatedly initialize
  ! them. This should save on compute time and provide a centralized location
  ! to acquire ia and ja arrays.
  
  implicit none

  public :: csr_tridiag_ia, csr_tridiag_ja, &
            csr_banddiag5_135_ia, csr_banddiag5_135_ja, &
            csr_banddiag5_12345_ia, csr_banddiag5_12345_ja, &
            initialize_csr_matrix, &
            ia_size, tridiag_ja_size, band12345_ja_size, band135_ja_size, &
            csr_intlc_s3b_f5b_ia, csr_intlc_s3b_f5b_ja, &
            csr_intlc_trid_5b_ia, csr_intlc_trid_5b_ja, &
            csr_intlc_5b_5b_ia, csr_intlc_5b_5b_ja, &
            intlc_ia_size, intlc_s3d_5d_ja_size, intlc_5d_5d_ja_size, &
            intlc_td_5d_ja_size

  private ! Default scope

  integer, allocatable, dimension(:) :: &
    csr_tridiag_ia, & !_ia array description for a tridiagonal matrix
    csr_tridiag_ja, & !_ja array description for a tridiagonal matrix
    csr_banddiag5_135_ia, & !_ia array description for a 5-band matrix
    !                        with the first upper and lower bands as 0.
    csr_banddiag5_135_ja, & !_ja array description for a 5-band matrix
    !                        with the first upper and lower bands as 0.
    csr_banddiag5_12345_ia, & !_ia array description for a 5-band matrix
    csr_banddiag5_12345_ja, & !_ja array description for a 5-band matrix
    csr_intlc_s3b_f5b_ia,   & !_ia array description for interlaced 5-band
    !                          matrix ("spaced 3-band, full 5-band")
    csr_intlc_s3b_f5b_ja,   & !_ja array description for interlaced 5-band
    !                          matrix ("spaced 3-band, full 5-band")
    csr_intlc_trid_5b_ia,   & !_ia array description for interlaced tridiag
    !                          and 5-band matrix (tridiag, 5-band)
    csr_intlc_trid_5b_ja,   & !_ja array description for interlaced tridiag
    !                          and 5-band matrix (tridiag, 5-band)
    csr_intlc_5b_5b_ia,     & !_ia array description for "interlaced"
    !                          5-band and 5-band matrix (double-size 5-band)
    csr_intlc_5b_5b_ja        !_ja array description for "interlaced"
    !                          5-band and 5-band matrix (double-size 5-band)

  integer :: &
    ia_size, & ! Size of the _ia arrays.
    tridiag_ja_size, &   ! Size of the tridiagonal ja array.
    band12345_ja_size, & ! Size of the 5-band-with-first-bands-0 ja array.
    band135_ja_size, &   ! Size of the 5-band ja array.
    intlc_ia_size, &     ! Size of the interlaced _ia arrays.
    intlc_s3d_5d_ja_size, & ! Size of the interlaced spaced 
                            ! 3-diag+5-diag ja arrays.
    intlc_5d_5d_ja_size, & ! Size of the interlaced 5-diag+5-diag ja arrays.
    intlc_td_5d_ja_size  ! Size of the interlaced tridiag+5-diag ja arrays.

!$omp threadprivate (csr_tridiag_ia, csr_tridiag_ja)
!$omp threadprivate (csr_banddiag5_135_ia, csr_banddiag5_135_ja)
!$omp threadprivate (csr_banddiag5_12345_ia, csr_banddiag5_12345_ja)
!$omp threadprivate (ia_size, tridiag_ja_size, band12345_ja_size, band135_ja_size)
!$omp threadprivate (csr_intlc_s3b_f5b_ia, csr_intlc_s3b_f5b_ja)
!$omp threadprivate (csr_intlc_trid_5b_ia, csr_intlc_trid_5b_ja)
!$omp threadprivate (csr_intlc_5b_5b_ia, csr_intlc_5b_5b_ja)
!$omp threadprivate (intlc_ia_size, intlc_s3d_5d_ja_size, intlc_5d_5d_ja_size)
!$omp threadprivate (intlc_td_5d_ja_size)

  contains

  !============================================================================
  subroutine initialize_csr_matrix

    ! Description:
    ! PARDISO matrix array initialization
    !
    ! This subroutine creates the _ia and _ja arrays, and calculates their
    ! required values for the current gr%nz.
    !
    ! References:
    !   None
    !------------------------------------------------------------------------

    use constants_clubb, only: &
      fstderr ! Variable(s)

    use grid_class, only: &
      gr ! Variable(s)

    implicit none

    ! Local variables
    integer :: &
      i, j, &              ! Loop indices
      error, &             ! Status for allocation
      num_bands, &         ! Number of diagonals for allocation
      num_diags, &         ! Number of non-empty diagonals for allocation
      cur_row, &           ! Current row--used in initialization
      cur_diag, & ! Current diagonal--num_diags/2 + 1 is the main diagonal
                  ! Note: At the boundaries, less diagonals are in scope.
                  ! At the lower boundaries, the subdiagonals aren't in scope.
                  ! At the upper boundaries, the superdiagonals aren't in scope.
      counter     ! Counter used to initialize the interlaced matrices

    logical :: l_print_ia_ja ! Debug flag to print the ia and ja arrays after
                           ! initialization is complete.

    ! ---- Begin Code ----

    ! Define the array sizes
    ia_size = gr%nz + 1
    intlc_ia_size = (2 * gr%nz) + 1
    
    ! Tridiagonal case and 5-band with 2 empty diagonals have 3 full diagonals
    num_diags = 3
    tridiag_ja_size = (gr%nz * num_diags) - 2
    band135_ja_size = (gr%nz * num_diags) - 4

    ! 5-band with all diagonals has 5 full diagonals
    num_diags = 5
    band12345_ja_size = (gr%nz * num_diags) - 6

    ! Interlaced arrays are tricky--there is an average of 4 diagonals for
    ! the 3/5band, but we need to take into account the fact that the
    ! tridiagonal and spaced 3-band will have different boundary indices.
    num_diags = 4
    intlc_td_5d_ja_size = (gr%nz * 2 * num_diags) - 4
    intlc_s3d_5d_ja_size = (gr%nz * 2 * num_diags) - 5

    ! The double-sized "interlaced" 5-band is similar to the standard 5-band
    num_diags = 5
    intlc_5d_5d_ja_size = (gr%nz * 2 * num_diags) - 6

    ! Allocate the correct amount of space for the actual _ia and _ja arrays
    allocate( csr_tridiag_ia(1:ia_size), &
              csr_tridiag_ja(1:tridiag_ja_size), &
              csr_banddiag5_12345_ia(1:ia_size), &
              csr_banddiag5_12345_ja(1:band12345_ja_size), &
              csr_banddiag5_135_ia(1:ia_size), &
              csr_banddiag5_135_ja(1:band135_ja_size), &
              csr_intlc_s3b_f5b_ia(1:intlc_ia_size), &
              csr_intlc_s3b_f5b_ja(1:intlc_s3d_5d_ja_size), &
              csr_intlc_trid_5b_ia(1:intlc_ia_size), &
              csr_intlc_trid_5b_ja(1:intlc_td_5d_ja_size), &
              csr_intlc_5b_5b_ia(1:intlc_ia_size), &
              csr_intlc_5b_5b_ja(1:intlc_5d_5d_ja_size), &
              stat=error )

    if ( error /= 0 ) then
      write(fstderr,*) "Allocation of CSR matrix arrays failed."
      stop "Fatal error--allocation of CSR matrix arrays failed."
    end if

    ! Initialize the tridiagonal matrix arrays
    num_bands = 3
    do i = 2, (gr%nz - 1), 1
      cur_row = (i - 1) * num_bands
      do j = 1, num_bands, 1
        cur_diag = j - 1
        csr_tridiag_ja(cur_row + cur_diag) = i + j - 2
      end do
      csr_tridiag_ia(i) = cur_row
    end do ! i = 2...gr%nz-1

    ! Handle boundary conditions for the tridiagonal matrix arrays
    ! These conditions have been hand-calculated bearing in mind that the
    ! matrix in question is tridiagonal.
    
    ! Make sure we don't crash if someone sets up gr%nz as 1.
    if ( gr%nz > 1 ) then
      ! Lower boundaries
      csr_tridiag_ja(1) = 1
      csr_tridiag_ja(2) = 2
      csr_tridiag_ia(1) = 1

      ! Upper boundaries
      csr_tridiag_ja(tridiag_ja_size - 1) = gr%nz - 1
      csr_tridiag_ja(tridiag_ja_size) = gr%nz
      csr_tridiag_ia(ia_size - 1) = tridiag_ja_size - 1

      ! This final boundary is to signify the end of the matrix, and is
      ! intended to be beyond the bound of the ja array.
      csr_tridiag_ia(ia_size) = tridiag_ja_size + 1
    end if ! gr%nz > 1

    ! Initialize the 5-band matrix arrays
    num_bands = 5
    do i = 3, (gr%nz - 2), 1

      ! Full 5-band matrix has 5 diagonals to initialize
      num_diags = 5
      cur_row = num_diags * (i - 1)
      do j = 1, num_diags, 1
        cur_diag = j - 3
        csr_banddiag5_12345_ja(cur_row + cur_diag) = i + cur_diag
      end do

      csr_banddiag5_12345_ia(i) = cur_row - 2

      ! 5-band matrix with 2 zero bands has 3 diagonals to initialize
      num_diags = 3
      cur_row = num_diags * (i - 1)
      do j = 1, num_diags, 1
        cur_diag = j - 2
        ! The first upper and first lower bands are zero, so there needs to be
        ! special handling to account for this. The j * 2 takes into account
        ! the spaces between diagonals.
        csr_banddiag5_135_ja(cur_row + cur_diag) = i + ((j * 2) - 1) - num_diags
      end do

      csr_banddiag5_135_ia(i) = cur_row - 1

    end do ! i = 3...gr%nz-2

    ! Handle boundary conditions for the 5-band matrix arrays
    ! These values have been hand-calculated bearing in mind the two different
    ! types of 5-band matrices.

    ! Make sure we don't crash if someone sets up gr%nz as less than 3.
    if ( gr%nz > 2 ) then

      ! -------------- (full) 5-band matrix boundaries ---------------

      ! Lower boundaries for the (full) 5-band matrix.
      do i = 1, 3, 1
        csr_banddiag5_12345_ja(i) = i
      end do
      do i = 1, 4, 1
        csr_banddiag5_12345_ja(i + 3) = i
      end do
      csr_banddiag5_12345_ia(1) = 1
      csr_banddiag5_12345_ia(2) = 4

      ! Upper boundaries for the (full) 5-band matrix.
      ! 7 and 3 are the number of elements from the "end" of the matrix if we
      ! travel right to left, bottom up. Because the ja matrices correspond to
      ! the column the element is in, we go 3 or 4 elements from the end for the
      ! second to last row (both superdiagonals absent on last row),
      ! and 3 for the last row (both superdiagonals absent). The indices are
      ! similarly calculated, except that in the case of the second to last
      ! row, it is necessary to offset for the last row as well (hence,
      ! 7 = 4+3).
      do i = 1, 4, 1
        csr_banddiag5_12345_ja(band12345_ja_size - 7 + i) = gr%nz + i - 4
      end do
      do i = 1, 3, 1
        csr_banddiag5_12345_ja(band12345_ja_size - 3 + i) = gr%nz + i - 3
      end do
      csr_banddiag5_12345_ia(ia_size - 2) = band12345_ja_size - 6
      csr_banddiag5_12345_ia(ia_size - 1) = band12345_ja_size - 2

      ! This final boundary is to signify the end of the matrix, and is
      ! intended to be beyond the bound of the ja array.
      csr_banddiag5_12345_ia(ia_size) = band12345_ja_size + 1

      ! ------------ end (full) 5-band matrix boundaries ---------------

      ! --------- 5-band matrix w/ empty first bands boundaries ----------

      ! Lower boundaries for the 5-band w/ empty first bands matrix
      ! The 2 * i is present because of the space between the main diagonal
      ! and the superdiagonal that actually have nonzero values.
      do i = 1, 2, 1
        csr_banddiag5_135_ja(i) = (2 * i) - 1
        csr_banddiag5_135_ja(i + 2) = (2 * i)
        csr_banddiag5_135_ia(i) = (2 * i) - 1
      end do

      ! Upper boundaries for the 5-band w/ empty first bands matrix
      ! The values for the boundaries are tricky, as the indices and values
      ! are not equal. The indices are 2 and 4 away from the end, as there are
      ! only two nonzero values at the two final rows.
      ! The values, on the other hand, are different, because of the
      ! aforementioned space, this time between the main and subdiagonal.
      do i = 1, 2, 1
        csr_banddiag5_135_ja(band135_ja_size - 4 + i) = gr%nz + (i * 2) - 5
        csr_banddiag5_135_ja(band135_ja_size - 2 + i) = gr%nz + (i * 2) - 4
      end do
      csr_banddiag5_135_ia(ia_size - 2) = band135_ja_size - 3
      csr_banddiag5_135_ia(ia_size - 1) = band135_ja_size + 1

      ! This final boundary is to signify the end of the matrix, and is
      ! intended to be beyond the bound of the ja array.
     csr_banddiag5_135_ia(ia_size) = band135_ja_size + 1

      ! ------- end 5-band matrix w/ empty first bands boundaries --------

    end if ! gr%nz > 2

    ! Initialize the interlaced arrays--all of them are 5-band right now.
    num_bands = 5

    ! Our counter starts at 2--this is used for the 3/5 interlaced matrices.
    ! We start at 2 so when we enter the odd row and increment by 5,
    ! it becomes 7.
    counter = 2

    do i = 3, ((gr%nz * 2) - 2), 1
      if (mod( i,2 ) == 1) then
        ! Odd row, this is the potentially non 5-band row.
        ! Increment counter. Last row was an even row, so we'll need to add 5.
        counter = counter + 5

        ! For our tridiag and spaced 3-band arrays, this will be a
        ! 3-diagonal row.
        num_diags = 3
        cur_row = counter + 1
        do j = 1, num_diags, 1
          cur_diag = j - 2
          csr_intlc_s3b_f5b_ja(cur_row + cur_diag) &
            = i + ((j * 2) - 1) - num_diags
          csr_intlc_trid_5b_ja(cur_row + cur_diag) = i + cur_diag
        end do
        csr_intlc_s3b_f5b_ia(i) = counter
        csr_intlc_trid_5b_ia(i) = counter

        ! For our 5-band interlaced-size array, this will be a
        ! 5-diagonal row (obviously!).
        num_diags = 5
        cur_row = num_diags * (i - 1)
        do j = 1, num_diags, 1
          cur_diag = j - 3
          csr_intlc_5b_5b_ja(cur_row + cur_diag) = i + cur_diag
        end do

        csr_intlc_5b_5b_ia(i) = cur_row - 2

      else
        ! Even row, this is the "guaranteed" 5-band row.
        ! Increment counter. Last row was an odd row, so we'll need to add 3.
        counter = counter + 3

        ! For our tridiag and spaced 3-band arrays, this will be a
        ! 5-diagonal row.
        num_diags = 5
        cur_row = counter + 2
        do j = 1, num_diags, 1
          cur_diag = j - 3
          csr_intlc_s3b_f5b_ja(cur_row + cur_diag) = i + cur_diag
          csr_intlc_trid_5b_ja(cur_row + cur_diag) = i + cur_diag
        end do

        csr_intlc_s3b_f5b_ia(i) = counter
        csr_intlc_trid_5b_ia(i) = counter
 
        ! For our 5-band "interlaced" array, this will also be a
        ! 5-diagonal row. However, we need to change the cur_row to match
        ! what we're expecting for the 5-band.
        num_diags = 5
        cur_row = num_diags * (i - 1)
        do j = 1, num_diags, 1
          cur_diag = j - 3
          csr_intlc_5b_5b_ja(cur_row + cur_diag) = i + cur_diag
        end do

        csr_intlc_5b_5b_ia(i) = cur_row - 2

      end if ! mod(i,2) == 1
    end do ! i = 3...(gr%nz*2)-2

    ! Handle boundary conditions for the interlaced matrix arrays
    ! These conditions have been hand-calculated bearing in mind
    ! the structure of the interlaced matrices.

    ! Make sure we don't crash if someone sets up gr%nz as less than 3.
    if (gr%nz > 2) then
      ! Lower boundaries

      ! First row
      do i = 1, 2, 1
        csr_intlc_s3b_f5b_ja(i) = (i * 2) - 1
        csr_intlc_trid_5b_ja(i) = i
      end do
      do i = 1, 3, 1
        csr_intlc_5b_5b_ja(i) = i
      end do
      csr_intlc_s3b_f5b_ia(1) = 1
      csr_intlc_trid_5b_ia(1) = 1
      csr_intlc_5b_5b_ia(1) = 1

      ! Second row
      do i = 1, 4, 1
        csr_intlc_s3b_f5b_ja(i + 2) = i
        csr_intlc_trid_5b_ja(i + 2) = i
        csr_intlc_5b_5b_ja(i + 3) = i
      end do
      csr_intlc_s3b_f5b_ia(2) = 3
      csr_intlc_trid_5b_ia(2) = 3
      csr_intlc_5b_5b_ia(2) = 4

      ! Upper boundaries

      ! Last two rows
      ! Note that in comparison to the other upper boundaries, we have to use
      ! intlc_ia_size - 1 for our upper index limit as the matrix is
      ! double-sized.

      ! Second-to-last row
      do i = 1, 2, 1
        csr_intlc_s3b_f5b_ja(intlc_s3d_5d_ja_size - 5 + i) &
          = intlc_ia_size - 1 + (i * 2) - 5
      end do
      do i = 1, 3, 1
        csr_intlc_trid_5b_ja(intlc_td_5d_ja_size - 6 + i) &
          = intlc_ia_size - 1 + i - 3
      end do
      do i = 1, 4, 1
        csr_intlc_5b_5b_ja(intlc_5d_5d_ja_size - 7 + i) &
          = intlc_ia_size-1 + i - 4
      end do

      ! Last row
      do i = 1, 3, 1
        csr_intlc_s3b_f5b_ja(intlc_s3d_5d_ja_size - 3 + i) &
          = intlc_ia_size-1 + i - 3
        csr_intlc_trid_5b_ja(intlc_td_5d_ja_size - 3 + i) &
          = intlc_ia_size-1 + i - 3
        csr_intlc_5b_5b_ja(intlc_5d_5d_ja_size - 3 + i) &
          = intlc_ia_size-1 + i - 3
      end do

      ! Lastly, take care of the ia arrays.
      csr_intlc_s3b_f5b_ia(intlc_ia_size - 2) = intlc_s3d_5d_ja_size - 4
      csr_intlc_s3b_f5b_ia(intlc_ia_size - 1) = intlc_s3d_5d_ja_size - 2
      csr_intlc_s3b_f5b_ia(intlc_ia_size) = intlc_s3d_5d_ja_size + 1

      csr_intlc_trid_5b_ia(intlc_ia_size - 2) = intlc_td_5d_ja_size - 5
      csr_intlc_trid_5b_ia(intlc_ia_size - 1) = intlc_td_5d_ja_size - 2
      csr_intlc_trid_5b_ia(intlc_ia_size) = intlc_td_5d_ja_size + 1

      csr_intlc_5b_5b_ia(intlc_ia_size - 2) = intlc_5d_5d_ja_size - 6
      csr_intlc_5b_5b_ia(intlc_ia_size - 1) = intlc_5d_5d_ja_size - 2
      csr_intlc_5b_5b_ia(intlc_ia_size) = intlc_5d_5d_ja_size + 1

        
    end if ! gr%nz > 2

    ! Enable printing the ia/ja arrays for debug purposes
    l_print_ia_ja = .false.
    if (l_print_ia_ja) then
      do i = 1, ia_size, 1
        print *, "tridiag ia idx", i, "=", csr_tridiag_ia(i)
        print *, "banddiag12345 ia idx", i, "=", csr_banddiag5_12345_ia(i)
        print *, "banddiag135 ia idx", i, "=", csr_banddiag5_135_ia(i)
      end do
      do i = 1, intlc_ia_size, 1
        print *, "interlaced tridiag w/ 5-band ia idx", i, &
          "=", csr_intlc_trid_5b_ia(i)
        print *, "interlaced spaced-3-band+5-band ia idx", i, &
          "=", csr_intlc_s3b_f5b_ia(i)
        print *, "interlaced 5-band w/ 5-band ia idx", i, "=", &
          csr_intlc_5b_5b_ia(i)
      end do
      do i = 1, tridiag_ja_size, 1
        print *, "tridiag ja idx", i, "=", csr_tridiag_ja(i)
      end do
      do i = 1, band12345_ja_size, 1
        print *, "band12345 ja idx", i, "=", csr_banddiag5_12345_ja(i)
      end do
      do i = 1, band135_ja_size, 1
        print *, "band135 ja idx", i, "=", csr_banddiag5_135_ja(i)
      end do
      do i = 1, intlc_td_5d_ja_size, 1
        print *, "interlaced tridiag w/ 5-band ja idx", i, &
          "=", csr_intlc_trid_5b_ja(i)
      end do
      do i = 1, intlc_s3d_5d_ja_size, 1
        print *, "interlaced spaced-3-band+5-band ja idx", i, &
          "=", csr_intlc_s3b_f5b_ja(i)
      end do
      do i = 1, intlc_5d_5d_ja_size, 1
        print *, "interlaced 5-band w/ 5-band ja idx", i, "=", &
          csr_intlc_5b_5b_ja(i)
      end do
    end if ! l_print_ia_ja

  end subroutine initialize_csr_matrix

end module csr_matrix_module
