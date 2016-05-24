!-----------------------------------------------------------------------
!$Id: corr_varnce_module.F90 7130 2014-07-29 23:29:54Z raut@uwm.edu $
!-------------------------------------------------------------------------------
module corr_varnce_module

  use clubb_precision, only: &
      core_rknd

  implicit none

  type sigma2_on_mu2_ratios_type

    ! In CLUBB standalone, these parameters can be set based on the value for a
    ! given case in the CASE_model.in file.

    ! Prescribed parameters for hydrometeor in-precip values of
    ! sigma_hm_i^2 / mu_hm_i^2 at grid levels that have some cloud.
    ! They can be set based on values for a given case in the CASE_model.in file.
    real( kind = core_rknd ) :: &
      rr_sigma2_on_mu2_ip_cloud = 1.0_core_rknd, & ! sigma_rr_i^2/mu_rr_i^2  [-]
      Nr_sigma2_on_mu2_ip_cloud = 1.0_core_rknd    ! sigma_Nr_i^2/mu_Nr_i^2  [-]

    ! Prescribed parameters for hydrometeor in-precip values of
    ! sigma_hm_i^2 / mu_hm_i^2 at grid levels that are entirely clear.
    ! They can be set based on values for a given case in the CASE_model.in file.
    real( kind = core_rknd ) :: &
      rr_sigma2_on_mu2_ip_below = 1.0_core_rknd, & ! sigma_rr_i^2/mu_rr_i^2  [-]
      Nr_sigma2_on_mu2_ip_below = 1.0_core_rknd    ! sigma_Nr_i^2/mu_Nr_i^2  [-]

    ! Parameters added for ice microphysics and latin hypercube sampling
    real( kind = core_rknd ) :: &
      rs_sigma2_on_mu2_ip_cloud = 1.0_core_rknd, & ! sigma_rs_i^2/mu_rs_i^2  [-]
      Ns_sigma2_on_mu2_ip_cloud = 1.0_core_rknd, & ! sigma_Ns_i^2/mu_Ns_i^2  [-]
      ri_sigma2_on_mu2_ip_cloud = 1.0_core_rknd, & ! sigma_ri_i^2/mu_ri_i^2  [-]
      Ni_sigma2_on_mu2_ip_cloud = 1.0_core_rknd, & ! sigma_Ni_i^2/mu_Ni_i^2  [-]
      rg_sigma2_on_mu2_ip_cloud = 1.0_core_rknd, & ! sigma_rg_i^2/mu_rg_i^2  [-]
      Ng_sigma2_on_mu2_ip_cloud = 1.0_core_rknd    ! sigma_Ng_i^2/mu_Ng_i^2  [-]

    ! Parameters added for ice microphysics and latin hypercube sampling
    real( kind = core_rknd ) :: &
      rs_sigma2_on_mu2_ip_below = 1.0_core_rknd, & ! sigma_rs_i^2/mu_rs_i^2  [-]
      Ns_sigma2_on_mu2_ip_below = 1.0_core_rknd, & ! sigma_Ns_i^2/mu_Ns_i^2  [-]
      ri_sigma2_on_mu2_ip_below = 1.0_core_rknd, & ! sigma_ri_i^2/mu_ri_i^2  [-]
      Ni_sigma2_on_mu2_ip_below = 1.0_core_rknd, & ! sigma_Ni_i^2/mu_Ni_i^2  [-]
      rg_sigma2_on_mu2_ip_below = 1.0_core_rknd, & ! sigma_rg_i^2/mu_rg_i^2  [-]
      Ng_sigma2_on_mu2_ip_below = 1.0_core_rknd    ! sigma_Ng_i^2/mu_Ng_i^2  [-]

    ! Prescribed parameter for <N_cn'^2> / <N_cn>^2 at any grid level.
    ! NOTE: In the case that l_const_Nc_in_cloud is true, Ncn is constant
    !       throughout the entire grid box, so the parameter below should be
    !       ignored.
    real( kind = core_rknd ) :: &
      Ncnp2_on_Ncnm2 = 1.0_core_rknd   ! Prescribed ratio <N_cn'^2>/<N_cn>^2 [-]

  end type sigma2_on_mu2_ratios_type

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

  real( kind = core_rknd ), public, dimension(:), allocatable :: &
    sigma2_on_mu2_ip_array_cloud, &
    sigma2_on_mu2_ip_array_below

  real( kind = core_rknd ), public, dimension(:,:), allocatable :: &
    corr_array_cloud, &
    corr_array_below
!$omp threadprivate(sigma2_on_mu2_ip_array_cloud, sigma2_on_mu2_ip_array_below, &
!$omp   corr_array_cloud, corr_array_below)

  real( kind = core_rknd ), public, dimension(:,:), allocatable :: &
      corr_array_cloud_def, &
      corr_array_below_def
!$omp threadprivate( corr_array_cloud_def, corr_array_below_def )


  private

  public :: sigma2_on_mu2_ratios_type, read_correlation_matrix, setup_pdf_indices, &
            setup_corr_varnce_array, cleanup_corr_matrix_arrays, &
            assert_corr_symmetric, print_corr_matrix

  private :: get_corr_var_index, return_pdf_index, def_corr_idx


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
    allocate( corr_array_cloud_def(d_var_total,d_var_total) )
    allocate( corr_array_below_def(d_var_total,d_var_total) )

    ! Initialize all values to 0.
    corr_array_cloud_def = zero
    corr_array_below_def = zero

    ! Set the correlation of any variable with itself to 1.
    do indx = 1, d_var_total, 1
       corr_array_cloud_def(indx,indx) = one
       corr_array_below_def(indx,indx) = one
    enddo

    ! Set up default correlation arrays.
    ! The default correlation arrays used here are the correlation arrays used
    ! for the ARM 97 case.  Any changes should be made concurrently here and in
    ! ../../input/case_setups/arm_97_corr_array_cloud.in (for "in-cloud") and
    ! in ../../input/case_setups/arm_97_corr_array_cloud.in (for "below-cloud").
    corr_array_cloud_def = reshape( &

(/1._c, -.6_c, .09_c , .09_c , .5_c   , .5_c   , .2_c   , .2_c  , .2_c  , .2_c  , .2_c, .2_c, &! chi
  0._c, 1._c , .027_c, .027_c, .0726_c, .0855_c, -.024_c, .084_c, .018_c, .012_c, 0._c, 0._c, &! eta
  0._c, 0._c , 1._c  , .34_c , 0.2_c  , 0.2_c  ,  .1_c  , .15_c , 0._c  , 0._c  , 0._c, 0._c, &! w
  0._c, 0._c , 0._c  , 1._c  , 0._c   , 0._c   ,  .39_c , .29_c , .14_c , .21_c , 0._c, 0._c, &! Ncn
  0._c, 0._c , 0._c  , 0._c  , 1._c   , .7_c   ,  0._c  , 0._c  , .1_c  , .1_c  , .2_c, .2_c, &! rr
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 1._c   ,  .1_c  , .1_c  , 0._c  , 0._c  , .2_c, .2_c, &! Nr
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  1._c  , .7_c  , .5_c  , .5_c  , .3_c, .3_c, &! ri
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 1._c  , .5_c  , .5_c  , .3_c, .3_c, &! Ni
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 1._c  , .7_c  , .4_c, .4_c, &! rs
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 1._c  , .4_c, .4_c, &! Ns
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 0._c  , 1._c, .7_c, &! rg
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 0._c  , 0._c, 1._c/),&!Ng

    shape(corr_array_cloud_def) )
!  chi   eta     w      Ncn      rr       Nr        ri      Ni      rs      Ns      rg    Ng

    corr_array_cloud_def = transpose( corr_array_cloud_def )


    corr_array_below_def = reshape( &

(/1._c, .3_c , .09_c , .09_c , .5_c   , .5_c   , .2_c   , .2_c  , .2_c  , .2_c  , .2_c, .2_c, &! chi
  0._c, 1._c , .027_c, .027_c, .0726_c, .0855_c, -.024_c, .084_c, .018_c, .012_c, 0._c, 0._c, &! eta
  0._c, 0._c , 1._c  , .34_c , 0.2_c  , 0.2_c  ,  .1_c  , .15_c , 0._c  , 0._c  , 0._c, 0._c, &! w
  0._c, 0._c , 0._c  , 1._c  , 0._c   , 0._c   ,  .39_c , .29_c , .14_c , .21_c , 0._c, 0._c, &! Ncn
  0._c, 0._c , 0._c  , 0._c  , 1._c   , .7_c   ,  0._c  , 0._c  , .1_c  , .1_c  , .2_c, .2_c, &! rr
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 1._c   ,  .1_c  , .1_c  , 0._c  , 0._c  , .2_c, .2_c, &! Nr
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  1._c  , .7_c  , .5_c  , .5_c  , .3_c, .3_c, &! ri
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 1._c  , .5_c  , .5_c  , .3_c, .3_c, &! Ni
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 1._c  , .7_c  , .4_c, .4_c, &! rs
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 1._c  , .4_c, .4_c, &! Ns
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 0._c  , 1._c, .7_c, &! rg
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 0._c  , 0._c, 1._c/),&!Ng

    shape(corr_array_below_def) )
!  chi   eta     w      Ncn      rr       Nr        ri      Ni      rs      Ns      rg    Ng

    corr_array_below_def = transpose( corr_array_below_def )


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

    corr_array_cloud = zero
    corr_array_below = zero

    do i = 1, d_variables
       corr_array_cloud(i,i) = one
       corr_array_below(i,i) = one
    enddo

    do i = 1, d_variables-1
       do j = i+1, d_variables
          if ( def_corr_idx(i) > def_corr_idx(j) ) then
             corr_array_cloud(j, i) = corr_array_cloud_def(def_corr_idx(j), def_corr_idx(i))
             corr_array_below(j, i) = corr_array_below_def(def_corr_idx(j), def_corr_idx(i))
          else
             corr_array_cloud(j, i) = corr_array_cloud_def(def_corr_idx(i), def_corr_idx(j))
             corr_array_below(j, i) = corr_array_below_def(def_corr_idx(i), def_corr_idx(j))
          endif
       enddo
    enddo

  end subroutine set_corr_arrays_to_default


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

    ! Initializing to zero means that correlations we don't have are assumed to be 0.
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

  !-------------------------------------------------------------------------------
  subroutine return_pdf_index( hydromet_index, pdf_count, pdf_index )

  ! Description:
  !   Set the Latin hypercube variable index if the hydrometeor exists
  ! References:
  !   None
  !-------------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      hydromet_index

    ! Input/Output Variables
    integer, intent(inout) :: &
      pdf_count

    ! Output Variables
    integer, intent(out) :: &
      pdf_index

    ! ---- Begin Code ----

    if ( hydromet_index > 0 ) then
      pdf_count = pdf_count + 1
      pdf_index = pdf_count
    else
      pdf_index = -1
    end if

    return
  end subroutine return_pdf_index

!===============================================================================
  subroutine setup_corr_varnce_array( input_file_cloud, input_file_below, &
                                      iunit, sigma2_on_mu2_ratios )

! Description:
!   Setup an array with the x'^2/xm^2 variables on the diagonal and the other
!   elements to be correlations between various variables.

! References:
!   None.
!-------------------------------------------------------------------------------

    use model_flags, only: &
      l_fix_chi_eta_correlations, & ! Variable(s)
      l_const_Nc_in_cloud

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

    type(sigma2_on_mu2_ratios_type), intent(in) :: &
      sigma2_on_mu2_ratios   ! Prescribed sigma^2 / mu^2 terms

    ! Local variables
    logical :: l_warning, corr_file_exist
    integer :: i

    ! ---- Begin Code ----

    allocate( corr_array_cloud(d_variables,d_variables) )
    allocate( corr_array_below(d_variables,d_variables) )

    allocate( sigma2_on_mu2_ip_array_cloud(d_variables) )
    allocate( sigma2_on_mu2_ip_array_below(d_variables) )

    sigma2_on_mu2_ip_array_cloud(:) = zero
    sigma2_on_mu2_ip_array_below(:) = zero

    ! corr_file_exist is true if the *_corr_array_cloud.in file exists
    ! Note: It is assumed that if the *_corr_array_cloud.in file exists
    !       then *_corr_array_below.in also exists
    inquire( file = input_file_cloud, exist = corr_file_exist )

    if ( corr_file_exist ) then

       call read_correlation_matrix( iunit, trim( input_file_cloud ), d_variables, & ! In
                                     corr_array_cloud ) ! Out

       call read_correlation_matrix( iunit, trim( input_file_below ), d_variables, & ! In
                                     corr_array_below ) ! Out

    else ! Read in default correlation matrices

       call clubb_debug( 1, "Warning: "//trim( input_file_cloud )//" was not found! " // &
                        "The default correlation arrays will be used." )

       call init_default_corr_arrays( )

       call set_corr_arrays_to_default( )

    endif

    ! Mirror the correlation matrices
    call mirror_lower_triangular_matrix( d_variables, corr_array_cloud )
    call mirror_lower_triangular_matrix( d_variables, corr_array_below )

    ! Sanity check to avoid confusing non-convergence results.
    if ( clubb_at_least_debug_level( 2 ) ) then

      if ( .not. l_fix_chi_eta_correlations .and. iiPDF_Ncn > 0 ) then
        l_warning = .false.
        do i = 1, d_variables
          if ( ( corr_array_cloud(i,iiPDF_Ncn) /= zero .or.  &
                 corr_array_below(i,iiPDF_Ncn) /= zero ) .and. &
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

    if ( iiPDF_Ncn > 0 ) then

      if ( l_const_Nc_in_cloud ) then
        ! Ncn is constant throughout every grid box!
        sigma2_on_mu2_ip_array_cloud(iiPDF_Ncn) = zero
        sigma2_on_mu2_ip_array_below(iiPDF_Ncn) = zero
      else
        sigma2_on_mu2_ip_array_cloud(iiPDF_Ncn) = sigma2_on_mu2_ratios%Ncnp2_on_Ncnm2
        sigma2_on_mu2_ip_array_below(iiPDF_Ncn) = sigma2_on_mu2_ratios%Ncnp2_on_Ncnm2
      end if

    end if

    if ( iiPDF_rr > 0 ) then
      sigma2_on_mu2_ip_array_cloud(iiPDF_rr) = sigma2_on_mu2_ratios%rr_sigma2_on_mu2_ip_cloud
      if ( iiPDF_Nr > 0 ) then
        sigma2_on_mu2_ip_array_cloud(iiPDF_Nr) = sigma2_on_mu2_ratios%Nr_sigma2_on_mu2_ip_cloud
      end if ! iiPDF_Nr > 0
    end if ! iiPDF_rr > 0

    if ( iiPDF_rs > 0 ) then
      sigma2_on_mu2_ip_array_cloud(iiPDF_rs) = sigma2_on_mu2_ratios%rs_sigma2_on_mu2_ip_cloud


      if ( iiPDF_Ns > 0 ) then
        sigma2_on_mu2_ip_array_cloud(iiPDF_Ns) = sigma2_on_mu2_ratios%Ns_sigma2_on_mu2_ip_cloud


      end if ! iiPDF_Ns > 0
    end if ! iiPDF_rs > 0

    if ( iiPDF_ri > 0 ) then
      sigma2_on_mu2_ip_array_cloud(iiPDF_ri) = sigma2_on_mu2_ratios%ri_sigma2_on_mu2_ip_cloud


      if ( iiPDF_Ni > 0 ) then
        sigma2_on_mu2_ip_array_cloud(iiPDF_Ni) = sigma2_on_mu2_ratios%Ni_sigma2_on_mu2_ip_cloud

      end if ! iiPDF_Ni > 0
    end if ! iiPDF_ri > 0

    ! Sampling for graupel (disabled)
    if ( iiPDF_rg > 0 ) then
      sigma2_on_mu2_ip_array_cloud(iiPDF_rg) = sigma2_on_mu2_ratios%rg_sigma2_on_mu2_ip_cloud


      if ( iiPDF_Ng > 0 ) then
        sigma2_on_mu2_ip_array_cloud(iiPDF_Ng) = sigma2_on_mu2_ratios%Ng_sigma2_on_mu2_ip_cloud


      end if ! iiPDF_Ng > 0
    end if ! iiPDF_rg > 0

    if ( iiPDF_rr > 0 ) then
      sigma2_on_mu2_ip_array_below(iiPDF_rr) = sigma2_on_mu2_ratios%rr_sigma2_on_mu2_ip_below



      if ( iiPDF_Nr > 0 ) then
        sigma2_on_mu2_ip_array_below(iiPDF_Nr) = sigma2_on_mu2_ratios%Nr_sigma2_on_mu2_ip_below


      end if ! iiPDF_Nr > 0
    end if ! iiPDF_rr > 0

    if ( iiPDF_rs > 0 ) then
      sigma2_on_mu2_ip_array_below(iiPDF_rs) = sigma2_on_mu2_ratios%rs_sigma2_on_mu2_ip_below


      if ( iiPDF_Ns > 0 ) then
        sigma2_on_mu2_ip_array_below(iiPDF_Ns) = sigma2_on_mu2_ratios%Ns_sigma2_on_mu2_ip_below

      end if ! iiPDF_Ns > 0
    end if ! iiPDF_rs > 0

    if ( iiPDF_ri > 0 ) then
      sigma2_on_mu2_ip_array_below(iiPDF_ri) = sigma2_on_mu2_ratios%ri_sigma2_on_mu2_ip_below


      if ( iiPDF_Ni > 0 ) then
        sigma2_on_mu2_ip_array_below(iiPDF_Ni) =  sigma2_on_mu2_ratios%Ni_sigma2_on_mu2_ip_below
      end if ! iiPDF_Ni > 0

    end if ! iiPDF_ri > 0

    if ( iiPDF_rg > 0 ) then
      sigma2_on_mu2_ip_array_below(iiPDF_rg) = sigma2_on_mu2_ratios%rg_sigma2_on_mu2_ip_below


      if ( iiPDF_Ng > 0 ) then
        sigma2_on_mu2_ip_array_below(iiPDF_Ng) = sigma2_on_mu2_ratios%Ng_sigma2_on_mu2_ip_below


      end if ! iiPDF_Ng > 0
    end if ! iiPDF_rg > 0

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

    if ( allocated( corr_array_cloud ) ) then
      deallocate( corr_array_cloud )
    end if

    if ( allocated( corr_array_below ) ) then
      deallocate( corr_array_below )
    end if

    if ( allocated( sigma2_on_mu2_ip_array_cloud ) ) then
      deallocate( sigma2_on_mu2_ip_array_cloud )
    end if

    if ( allocated( sigma2_on_mu2_ip_array_below ) ) then
      deallocate( sigma2_on_mu2_ip_array_below )
    end if

    if ( allocated( corr_array_cloud_def ) ) then
      deallocate( corr_array_cloud_def )
    end if

    if ( allocated( corr_array_below_def ) ) then
      deallocate( corr_array_below_def )
    end if

    return
  end subroutine cleanup_corr_matrix_arrays

  !-----------------------------------------------------------------------------
  subroutine assert_corr_symmetric( corr_array, & ! intent(in)
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
      intent(in) :: corr_array ! Correlation array to be checked

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
        if (abs(corr_array(n_col, n_row) - corr_array(n_row, n_col)) > tol) then
          l_error = .true.
        end if
        if (n_col == n_row .and. corr_array(n_col, n_row) /= 1.0_core_rknd) then
          l_error = .true.
        end if
      end do
    end do

    !Report if any errors are found
    if (l_error) then
      write(fstderr,*) "Error: Correlation array is non symmetric or formatted incorrectly."
      write(fstderr,*) corr_array
      stop
    end if

  end subroutine assert_corr_symmetric

  !-----------------------------------------------------------------------------
  subroutine print_corr_matrix( d_variables, & ! intent(in)
                                corr_array ) ! intent(in)

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
      intent(in) :: corr_array ! Correlation array to be printed

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
        write(str_array_value,'(F5.2)') corr_array(m,n)
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
