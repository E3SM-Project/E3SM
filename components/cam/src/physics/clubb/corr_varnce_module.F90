!-----------------------------------------------------------------------
!$Id: corr_varnce_module.F90 8001 2016-03-03 22:42:36Z raut@uwm.edu $
!-------------------------------------------------------------------------------
module corr_varnce_module

  use clubb_precision, only: &
      core_rknd

  implicit none

  type hmp2_ip_on_hmm2_ip_ratios_type

    ! In CLUBB standalone, these parameters can be set based on the value for a
    ! given case in the CASE_model.in file.

    ! Prescribed parameters for hydrometeor values of <hm|_ip'^2> / <hm|_ip>^2,
    ! where <hm|_ip> is the in-precip. mean of the hydrometeor and <hm|_ip'^2>
    ! is the in-precip. variance of the hydrometeor.
    ! They can be set based on values for a given case in the CASE_model.in file.
    real( kind = core_rknd ) :: &
      rrp2_ip_on_rrm2_ip = 1.0_core_rknd, & ! Ratio <rr|_ip'^2> / <rr|_ip>^2 [-]
      Nrp2_ip_on_Nrm2_ip = 1.0_core_rknd, & ! Ratio <Nr|_ip'^2> / <Nr|_ip>^2 [-]
      rip2_ip_on_rim2_ip = 1.0_core_rknd, & ! Ratio <ri|_ip'^2> / <ri|_ip>^2 [-]
      Nip2_ip_on_Nim2_ip = 1.0_core_rknd, & ! Ratio <Ni|_ip'^2> / <Ni|_ip>^2 [-]
      rsp2_ip_on_rsm2_ip = 1.0_core_rknd, & ! Ratio <rs|_ip'^2> / <rs|_ip>^2 [-]
      Nsp2_ip_on_Nsm2_ip = 1.0_core_rknd, & ! Ratio <Ns|_ip'^2> / <Ns|_ip>^2 [-]
      rgp2_ip_on_rgm2_ip = 1.0_core_rknd, & ! Ratio <rg|_ip'^2> / <rg|_ip>^2 [-]
      Ngp2_ip_on_Ngm2_ip = 1.0_core_rknd    ! Ratio <Ng|_ip'^2> / <Ng|_ip>^2 [-]

  end type hmp2_ip_on_hmm2_ip_ratios_type

  ! Prescribed parameter for <N_cn'^2> / <N_cn>^2.
  ! NOTE: In the case that l_const_Nc_in_cloud is true, Ncn is constant
  !       throughout the entire grid box, so the parameter below should be
  !       ignored.
  real( kind = core_rknd ), public :: &
    Ncnp2_on_Ncnm2 = 1.0_core_rknd   ! Prescribed ratio <N_cn'^2> / <N_cn>^2 [-]

!$omp threadprivate(Ncnp2_on_Ncnm2)

  ! Latin hypercube indices / Correlation array indices
  integer, public :: &
    iiPDF_chi = -1, &
    iiPDF_eta = -1, &
    iiPDF_w   = -1
!$omp threadprivate(iiPDF_chi, iiPDF_eta, iiPDF_w)

  integer, public :: &
   iiPDF_rr = -1, &
   iiPDF_rs = -1, &
   iiPDF_ri = -1, &
   iiPDF_rg = -1
!$omp threadprivate(iiPDF_rr, iiPDF_rs, iiPDF_ri, iiPDF_rg)

  integer, public :: &
   iiPDF_Nr  = -1, &
   iiPDF_Ns  = -1, &
   iiPDF_Ni  = -1, &
   iiPDF_Ng  = -1, &
   iiPDF_Ncn = -1
!$omp threadprivate(iiPDF_Nr, iiPDF_Ns, iiPDF_Ni, iiPDF_Ng, iiPDF_Ncn)

  integer, parameter, public :: &
    d_var_total = 12 ! Size of the default correlation arrays

  integer, public :: &
    d_variables
!$omp threadprivate(d_variables)

  real( kind = core_rknd ), dimension(:), allocatable, public :: &
     hmp2_ip_on_hmm2_ip

!$omp threadprivate(hmp2_ip_on_hmm2_ip)

  real( kind = core_rknd ), public, dimension(:,:), allocatable :: &
    corr_array_n_cloud, &
    corr_array_n_below
!$omp threadprivate(corr_array_n_cloud, corr_array_n_below)

  real( kind = core_rknd ), public, dimension(:,:), allocatable :: &
      corr_array_n_cloud_def, &
      corr_array_n_below_def
!$omp threadprivate( corr_array_n_cloud_def, corr_array_n_below_def )


  private

  public :: hmp2_ip_on_hmm2_ip_ratios_type, &
            read_correlation_matrix, setup_pdf_indices, &
            setup_corr_varnce_array, cleanup_corr_matrix_arrays, &
            assert_corr_symmetric, print_corr_matrix

  private :: get_corr_var_index, def_corr_idx


  contains

  !-----------------------------------------------------------------------------
  subroutine init_default_corr_arrays(  ) 

    ! Description:
    ! Initializes the default correlation arrays.
    !---------------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Constant(s)
        zero

    implicit none

    integer:: indx

    ! This "renaming" is used to shorten the matrix declarations below.
    integer, parameter :: c = core_rknd

    ! ---- Begin Code ----
 
    ! Allocate Arrays.
    allocate( corr_array_n_cloud_def(d_var_total,d_var_total) )
    allocate( corr_array_n_below_def(d_var_total,d_var_total) )

    ! Initialize all values to 0.
    corr_array_n_cloud_def = zero
    corr_array_n_below_def = zero

    ! Set the correlation of any variable with itself to 1.
    do indx = 1, d_var_total, 1
       corr_array_n_cloud_def(indx,indx) = one
       corr_array_n_below_def(indx,indx) = one
    enddo

    ! Set up default normal space correlation arrays.
    ! The default normal space correlation arrays used here are the normal space
    ! correlation arrays used for the ARM 97 case.  Any changes should be made
    ! concurrently here and in
    ! ../../input/case_setups/arm_97_corr_array_cloud.in (for "in-cloud") and
    ! in ../../input/case_setups/arm_97_corr_array_cloud.in (for "below-cloud").
    corr_array_n_cloud_def = reshape( &

(/1._c,-.6_c, .09_c , .09_c , .788_c, .675_c, .240_c, .222_c, .240_c, .222_c, .240_c, .222_c, &! chi
  0._c, 1._c, .027_c, .027_c, .114_c, .115_c,-.029_c, .093_c, .022_c, .013_c, 0._c  , 0._c  , &! eta
  0._c, 0._c, 1._c  , .34_c , .315_c, .270_c, .120_c, .167_c, 0._c  , 0._c  , 0._c  , 0._c  , &! w
  0._c, 0._c, 0._c  , 1._c  , 0._c  , 0._c  , .464_c, .320_c, .168_c, .232_c, 0._c  , 0._c  , &! Ncn
  0._c, 0._c, 0._c  , 0._c  , 1._c  , .821_c, 0._c  , 0._c  , .173_c, .164_c, .319_c, .308_c, &! rr
  0._c, 0._c, 0._c  , 0._c  , 0._c  , 1._c  , .152_c, .143_c, 0._c  , 0._c  , .285_c, .273_c, &! Nr
  0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .758_c, .585_c, .571_c, .379_c, .363_c, &! ri
  0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .571_c, .550_c, .363_c, .345_c, &! Ni
  0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .758_c, .485_c, .470_c, &! rs
  0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .470_c, .450_c, &! Ns
  0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .758_c, &! rg
  0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c/), &! Ng

    shape(corr_array_n_cloud_def) )
!  chi   eta    w      Ncn     rr      Nr      ri      Ni      rs      Ns      rg      Ng

    corr_array_n_cloud_def = transpose( corr_array_n_cloud_def )


    corr_array_n_below_def = reshape( &

(/1._c, .3_c, .09_c , .09_c , .788_c, .675_c, .240_c, .222_c, .240_c, .222_c, .240_c, .222_c, &! chi
  0._c, 1._c, .027_c, .027_c, .114_c, .115_c,-.029_c, .093_c, .022_c, .013_c, 0._c  , 0._c  , &! eta
  0._c, 0._c, 1._c  , .34_c , .315_c, .270_c, .120_c, .167_c, 0._c  , 0._c  , 0._c  , 0._c  , &! w
  0._c, 0._c, 0._c  , 1._c  , 0._c  , 0._c  , .464_c, .320_c, .168_c, .232_c, 0._c  , 0._c  , &! Ncn
  0._c, 0._c, 0._c  , 0._c  , 1._c  , .821_c, 0._c  , 0._c  , .173_c, .164_c, .319_c, .308_c, &! rr
  0._c, 0._c, 0._c  , 0._c  , 0._c  , 1._c  , .152_c, .143_c, 0._c  , 0._c  , .285_c, .273_c, &! Nr
  0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .758_c, .585_c, .571_c, .379_c, .363_c, &! ri
  0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .571_c, .550_c, .363_c, .345_c, &! Ni
  0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .758_c, .485_c, .470_c, &! rs
  0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .470_c, .450_c, &! Ns
  0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .758_c, &! rg
  0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c/), &! Ng

    shape(corr_array_n_below_def) )
!  chi   eta    w      Ncn     rr      Nr      ri      Ni      rs      Ns      rg      Ng

    corr_array_n_below_def = transpose( corr_array_n_below_def )


    return

  end subroutine init_default_corr_arrays

  !-----------------------------------------------------------------------------
  pure function def_corr_idx( iiPDF_x ) result(ii_def_corr)

    ! Description:
    !   Map from a iiPDF index to the corresponding index in the default 
    !   correlation arrays.
    !-----------------------------------------------------------------------------

    implicit none

    ! Constant Parameters

    ! Indices that represent the order in the default corr arrays
    ! (chi (old s), eta (old t), w, Ncn, rr, Nr, ri, Ni, rs, Ns, rg, Ng)
    integer, parameter :: &
    ii_chi = 1, &
    ii_eta = 2, &
    ii_w = 3,   &
    ii_Ncn = 4, &
    ii_rr = 5,  &
    ii_Nr = 6,  &
    ii_ri = 7,  &
    ii_Ni = 8,  &
    ii_rs = 9,  &
    ii_Ns = 10, &
    ii_rg = 11, &
    ii_Ng = 12

    ! Input Variables

    integer, intent(in) :: iiPDF_x

    ! Return Variable

    integer :: ii_def_corr

    ! ---- Begin Code ----

    ii_def_corr = -1

      if (iiPDF_x == iiPDF_chi) then
         ii_def_corr = ii_chi

      elseif (iiPDF_x == iiPDF_eta) then
        ii_def_corr = ii_eta

      elseif (iiPDF_x == iiPDF_w) then
        ii_def_corr = ii_w

      elseif (iiPDF_x == iiPDF_Ncn) then
        ii_def_corr = ii_Ncn

      elseif (iiPDF_x == iiPDF_rr) then
        ii_def_corr = ii_rr

      elseif (iiPDF_x == iiPDF_Nr) then
        ii_def_corr = ii_Nr

      elseif (iiPDF_x == iiPDF_ri) then
        ii_def_corr = ii_ri

      elseif (iiPDF_x == iiPDF_Ni) then
        ii_def_corr = ii_Ni

      elseif (iiPDF_x == iiPDF_rs) then
        ii_def_corr = ii_rs

      elseif (iiPDF_x == iiPDF_Ns) then
        ii_def_corr = ii_Ns

      elseif (iiPDF_x == iiPDF_rg) then
        ii_def_corr = ii_rg

      elseif (iiPDF_x == iiPDF_Ng) then
        ii_def_corr = ii_Ng

      endif
  end function def_corr_idx

  !-----------------------------------------------------------------------------
  subroutine set_corr_arrays_to_default(  ) 

    ! Description:
    !   If there are no corr_array.in files for the current case, default 
    !   correlations are used. 
    !-----------------------------------------------------------------------------
  
    use constants_clubb, only: &
        zero, &
        one

    implicit none

    ! Local Variables
    integer :: i, j ! Loop iterators


    ! ---- Begin Code ----

    corr_array_n_cloud = zero
    corr_array_n_below = zero

    do i = 1, d_variables
       corr_array_n_cloud(i,i) = one
       corr_array_n_below(i,i) = one
    enddo

    do i = 1, d_variables-1
       do j = i+1, d_variables
          if ( def_corr_idx(i) > def_corr_idx(j) ) then
             corr_array_n_cloud(j, i) = corr_array_n_cloud_def(def_corr_idx(j), def_corr_idx(i))
             corr_array_n_below(j, i) = corr_array_n_below_def(def_corr_idx(j), def_corr_idx(i))
          else
             corr_array_n_cloud(j, i) = corr_array_n_cloud_def(def_corr_idx(i), def_corr_idx(j))
             corr_array_n_below(j, i) = corr_array_n_below_def(def_corr_idx(i), def_corr_idx(j))
          endif
       enddo
    enddo

  end subroutine set_corr_arrays_to_default


  !-----------------------------------------------------------------------------
  subroutine read_correlation_matrix( iunit, input_file, d_variables, &
                                      corr_array_n )

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
      corr_array_n ! Normal space correlation array

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

    ! Initializing to zero means that correlations we don't have are assumed to be 0.
    corr_array_n(:,:) = 0.0_core_rknd

    ! Set main diagonal to 1
    do i=1, d_variables
      corr_array_n(i,i) = 1.0_core_rknd
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
                   corr_array_n )
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

    case( "chi" )
      i = iiPDF_chi

    case( "eta" )
      i = iiPDF_eta

    case( "w" )
      i = iiPDF_w

    case( "Ncn" )
      i = iiPDF_Ncn

    case( "rr" )
      i = iiPDF_rr

    case( "Nr" )
      i = iiPDF_Nr

    case( "ri" )
      i = iiPDF_ri

    case( "Ni" )
      i = iiPDF_Ni

    case( "rs" )
      i = iiPDF_rs

    case( "Ns" )
      i = iiPDF_Ns
        
    case( "rg" )
      i = iiPDF_rg

    case( "Ng" )
      i = iiPDF_Ng

    end select

    return

  end function get_corr_var_index

  !-----------------------------------------------------------------------
  subroutine setup_pdf_indices( hydromet_dim, iirrm, iiNrm, &
                                iirim, iiNim, iirsm, iiNsm, &
                                iirgm, iiNgm )

    ! Description:
    !
    ! Setup for the iiPDF indices. These indices are used to address chi(s), eta(t), w
    ! and the hydrometeors in the mean/stdev/corr arrays
    !
    ! References:
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      hydromet_dim    ! Total number of hydrometeor species.

    integer, intent(in) :: &
      iirrm, & ! Index of rain water mixing ratio
      iiNrm, & ! Index of rain drop concentration
      iirim, & ! Index of ice mixing ratio
      iiNim, & ! Index of ice crystal concentration
      iirsm, & ! Index of snow mixing ratio
      iiNsm, & ! Index of snow flake concentration
      iirgm, & ! Index of graupel mixing ratio
      iiNgm    ! Index of graupel concentration

    ! Local Variables
    integer :: &
      pdf_count, & ! Count number of PDF variables
      i            ! Hydrometeor loop index

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    iiPDF_chi = 1 ! Extended liquid water mixing ratio, chi
    iiPDF_eta = 2 ! 'eta' orthogonal to 'chi'
    iiPDF_w   = 3 ! vertical velocity
    iiPDF_Ncn = 4 ! Simplified cloud nuclei concentration or extended Nc.

    pdf_count = iiPDF_Ncn

    ! Loop over hydrometeors.
    ! Hydrometeor indices in the PDF arrays should be in the same order as
    ! found in the hydrometeor arrays.
    if ( hydromet_dim > 0 ) then

       do i = 1, hydromet_dim, 1

          if ( i == iirrm ) then
             pdf_count = pdf_count + 1
             iiPDF_rr = pdf_count
          endif

          if ( i == iiNrm ) then
             pdf_count = pdf_count + 1
             iiPDF_Nr = pdf_count
          endif

          if ( i == iirim ) then
             pdf_count = pdf_count + 1
             iiPDF_ri = pdf_count
          endif

          if ( i == iiNim ) then
             pdf_count = pdf_count + 1
             iiPDF_Ni = pdf_count
          endif

          if ( i == iirsm ) then
             pdf_count = pdf_count + 1
             iiPDF_rs = pdf_count
          endif

          if ( i == iiNsm ) then
             pdf_count = pdf_count + 1
             iiPDF_Ns = pdf_count
          endif

          if ( i == iirgm ) then
             pdf_count = pdf_count + 1
             iiPDF_rg = pdf_count
          endif
        
          if ( i == iiNgm ) then
             pdf_count = pdf_count + 1
             iiPDF_Ng = pdf_count
          endif   

       enddo ! i = 1, hydromet_dim, 1

    endif ! hydromet_dim > 0

    d_variables = pdf_count


    return

  end subroutine setup_pdf_indices
  !-----------------------------------------------------------------------

!===============================================================================
  subroutine setup_corr_varnce_array( input_file_cloud, input_file_below, &
                                      iunit )

! Description:
!   Setup an array with the x'^2/xm^2 variables on the diagonal and the other
!   elements to be correlations between various variables.

! References:
!   None.
!-------------------------------------------------------------------------------

    use model_flags, only: &
      l_fix_chi_eta_correlations    ! Variable(s)

    use matrix_operations, only: mirror_lower_triangular_matrix ! Procedure

    use constants_clubb, only: &
      fstderr, &  ! Constant(s)
      zero

    use error_code, only: &
      clubb_debug, & ! Procedure(s)
      clubb_at_least_debug_level

    implicit none

    ! External
    intrinsic :: max, epsilon, trim

    character(len=*), intent(in) :: &
      input_file_cloud, &    ! Path to the in cloud correlation file
      input_file_below       ! Path to the out of cloud correlation file

    ! Input Variables
    integer, intent(in) :: &
      iunit ! The file unit

    ! Local variables
    logical :: l_warning, l_corr_file_1_exist, l_corr_file_2_exist
    integer :: i

    ! ---- Begin Code ----

    allocate( corr_array_n_cloud(d_variables,d_variables) )
    allocate( corr_array_n_below(d_variables,d_variables) )

    inquire( file = input_file_cloud, exist = l_corr_file_1_exist )
    inquire( file = input_file_below, exist = l_corr_file_2_exist )

    if ( l_corr_file_1_exist .and. l_corr_file_2_exist ) then

       call read_correlation_matrix( iunit, trim( input_file_cloud ), d_variables, & ! In
                                     corr_array_n_cloud ) ! Out

       call read_correlation_matrix( iunit, trim( input_file_below ), d_variables, & ! In
                                     corr_array_n_below ) ! Out

    else ! Read in default correlation matrices

       call clubb_debug( 1, "Warning: "//trim( input_file_cloud )//" was not found! " // &
                        "The default correlation arrays will be used." )

       call init_default_corr_arrays( )

       call set_corr_arrays_to_default( )

    endif

    ! Mirror the correlation matrices
    call mirror_lower_triangular_matrix( d_variables, corr_array_n_cloud )
    call mirror_lower_triangular_matrix( d_variables, corr_array_n_below )

    ! Sanity check to avoid confusing non-convergence results.
    if ( clubb_at_least_debug_level( 2 ) ) then

      if ( .not. l_fix_chi_eta_correlations .and. iiPDF_Ncn > 0 ) then
        l_warning = .false.
        do i = 1, d_variables
          if ( ( corr_array_n_cloud(i,iiPDF_Ncn) /= zero .or.  &
                 corr_array_n_below(i,iiPDF_Ncn) /= zero ) .and. &
               i /= iiPDF_Ncn ) then
            l_warning = .true.
          end if
        end do ! 1..d_variables
        if ( l_warning ) then
          write(fstderr,*) "Warning: the specified correlations for chi" &
                           // " (old s) and Ncn are non-zero."
          write(fstderr,*) "The latin hypercube code will not converge to" &
                           // " the analytic solution using these settings."
        end if
       end if ! l_fix_chi_eta_correlations .and. iiPDF_Ncn > 0

    end if ! clubb_at_least_debug_level( 2 )


    return

  end subroutine setup_corr_varnce_array

  !-----------------------------------------------------------------------------
  subroutine cleanup_corr_matrix_arrays( )

    ! Description:
    !   De-allocate latin hypercube arrays
    ! References:
    !   None
    !---------------------------------------------------------------------------
    implicit none

    ! External
    intrinsic :: allocated

    ! ---- Begin Code ----

    if ( allocated( corr_array_n_cloud ) ) then
      deallocate( corr_array_n_cloud )
    end if

    if ( allocated( corr_array_n_below ) ) then
      deallocate( corr_array_n_below )
    end if

    if ( allocated( corr_array_n_cloud_def ) ) then
      deallocate( corr_array_n_cloud_def )
    end if

    if ( allocated( corr_array_n_below_def ) ) then
      deallocate( corr_array_n_below_def )
    end if


    return

  end subroutine cleanup_corr_matrix_arrays

  !-----------------------------------------------------------------------------
  subroutine assert_corr_symmetric( corr_array_n, & ! intent(in)
                                    d_variables ) ! intent(in)

    ! Description:
    !   Asserts that corr_matrix(i,j) == corr_matrix(j,i) for all indeces
    !   in the correlation array. If this is not the case, stops the program.
    ! References:
    !   None
    !---------------------------------------------------------------------------

    use constants_clubb, only: fstderr ! Constant(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables    ! Number of variables in the correlation array

    real( kind = core_rknd ), dimension(d_variables, d_variables), &
      intent(in) :: corr_array_n ! Normal space correlation array to be checked

    ! Local Variables

    ! tolerance used for real precision testing
    real( kind = core_rknd ), parameter :: tol = 1.0e-6_core_rknd

    integer :: n_row, n_col !indeces

    logical :: l_error !error found between the two arrays

    !----- Begin Code -----

    l_error = .false.

    !Do the check
    do n_col = 1, d_variables
      do n_row = 1, d_variables
        if (abs(corr_array_n(n_col, n_row) - corr_array_n(n_row, n_col)) > tol) then
          l_error = .true.
        end if
        if (n_col == n_row .and. corr_array_n(n_col, n_row) /= 1.0_core_rknd) then
          l_error = .true.
        end if
      end do
    end do

    !Report if any errors are found
    if (l_error) then
      write(fstderr,*) "Error: Correlation array is non symmetric or formatted incorrectly."
      write(fstderr,*) corr_array_n
      stop
    end if

  end subroutine assert_corr_symmetric

  !-----------------------------------------------------------------------------
  subroutine print_corr_matrix( d_variables, & ! intent(in)
                                corr_array_n ) ! intent(in)

    ! Description:
    !   Prints the correlation matrix to the console.
    ! References:
    !   None
    !---------------------------------------------------------------------------

    use clubb_precision, only: core_rknd

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables    ! Number of variables in the correlation array

    real( kind = core_rknd ), dimension(d_variables, d_variables), &
      intent(in) :: corr_array_n ! Normal space correlation array to be printed

    ! Local Variables
    integer :: n, & ! Loop indeces
               m, &
               current_character_index ! keeps track of the position in the string

    character(LEN=72) :: current_line ! The current line to be printed
    character(LEN=10) :: str_array_value

    !----- Begin Code -----

    current_character_index = 0

    do n = 1, d_variables
      do m = 1, d_variables
        write(str_array_value,'(F5.2)') corr_array_n(m,n)
        current_line = current_line(1:current_character_index)//str_array_value
        current_character_index = current_character_index + 6
      end do
      write(*, *) current_line
      current_line = ""
      current_character_index = 0
    end do

  end subroutine print_corr_matrix
  !-----------------------------------------------------------------------------

end module corr_varnce_module
