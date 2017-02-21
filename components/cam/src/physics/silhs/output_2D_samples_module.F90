!-------------------------------------------------------------------------------
! $Id: output_2D_samples_module.F90 7994 2016-02-24 20:22:34Z raut@uwm.edu $
!===============================================================================
module output_2D_samples_module

  use stat_file_module, only : stat_file ! Type

  implicit none

  public :: open_2D_samples_file, close_2D_samples_file, &
    output_2D_uniform_dist_file, output_2D_lognormal_dist_file

  private ! Default scope

  type(stat_file), public :: &
    lognormal_sample_file, &
    uniform_sample_file

  !$omp threadprivate( lognormal_sample_file, uniform_sample_file )

  contains
!-------------------------------------------------------------------------------
  subroutine open_2D_samples_file( nz, num_samples, n_2D_variables, &
                                   fname_prefix, fdir, &
                                   time, dtwrite, zgrid, variable_names, &
                                   variable_descriptions, variable_units, &
                                   sample_file )
! Description:
!   Open a 2D sample file
! References:
!   None
!-------------------------------------------------------------------------------
#ifdef NETCDF
    use output_netcdf, only: open_netcdf_for_writing ! Procedure(s)
#endif

    use clubb_precision, only: time_precision, core_rknd ! Constant(s)

    implicit none

    ! Parameter Constants
    integer, parameter :: &
      day   = 1, & ! Made up times for GrADS
      month = 1, &
      year  = 1900

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of vertical levels
      num_samples, & ! Number of samples per variable
      n_2D_variables ! Number variables to output

    character(len=*), intent(in) :: &
      fdir,      & ! Output directory
      fname_prefix ! Prefix for the netCDF output

    character(len=*), intent(in), dimension(n_2D_variables) :: &
      variable_names,        & ! Names of the variables to be used in the 2D netCDF file
      variable_descriptions, & ! Description of the variables in the 2D file
      variable_units           ! Units on the variables

    real(kind=time_precision), intent(in) :: & 
      time      ! Start time                      [s] 
    
    real(kind=core_rknd), intent(in) :: &
      dtwrite   ! Interval for writing to disk    [s] 

    real( kind = core_rknd ), intent(in), dimension(nz) :: & 
      zgrid ! Vertical grid levels [m]

    ! Input/Output Variables
    type(stat_file), intent(inout) :: &
      sample_file ! File that is being opened

    ! Local Variables
    integer :: nlat, nlon ! Not actually latitudes and longitudes

    real( kind = core_rknd ), dimension(num_samples) :: rlat

    real( kind = core_rknd ), dimension(1) :: rlon

    character(len=100) :: fname
    integer :: i

    ! ---- Begin Code ----

    fname = trim( fname_prefix )//"_lh_sample_points_2D"
    i =1  ! This assignment prevents a g 95 compiler warning

    ! We need to set this like a latitude to trick GrADS and allow of viewing of
    ! the sample points with the GrADS application and sdfopen.
    nlat = num_samples
    nlon = 1

    allocate( sample_file%rlat(num_samples), sample_file%rlon(1) )
    allocate( sample_file%var(n_2D_variables) )
    allocate( sample_file%z(nz) )

    forall( i=1:num_samples )
      rlat(i) = real( i, kind = core_rknd ) ! Use made up arbitrary values for degrees north
    end forall

    rlon = 1.0_core_rknd ! Also made up

    do i=1, n_2D_variables
      sample_file%var(i)%name = trim( variable_names(i) )
      sample_file%var(i)%description = trim( variable_descriptions(i) )
      sample_file%var(i)%units = trim( variable_units(i) )
    end do

#ifdef NETCDF
    call open_netcdf_for_writing( nlat, nlon, fdir, fname, 1, nz, zgrid, &
                      day, month, year, rlat, rlon, &
                      time, dtwrite, n_2D_variables, sample_file )
#else
    stop "This version of CLUBB was not compiled for netCDF output"
#endif

    return
  end subroutine open_2D_samples_file

!-------------------------------------------------------------------------------
  subroutine output_2D_lognormal_dist_file &
             ( nz, num_samples, d_variables, X_nl_all_levs )
! Description:
!   Output a 2D snapshot of latin hypercube samples
! References:
!   None
!-------------------------------------------------------------------------------
#ifdef NETCDF
    use output_netcdf, only: write_netcdf ! Procedure(s)
#endif

    use clubb_precision, only: stat_rknd, core_rknd ! Constant(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of vertical levels
      num_samples, & ! Number of samples per variable
      d_variables    ! Number variates being sampled

    real(kind=stat_rknd), intent(in), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer :: sample, j

    ! ---- Begin Code ----

    do j = 1, d_variables
      allocate( lognormal_sample_file%var(j)%ptr(num_samples,1,nz) )
    end do

    do sample = 1, num_samples
      do j = 1, d_variables
        lognormal_sample_file%var(j)%ptr(sample,1,1:nz) = X_nl_all_levs(1:nz,sample,j)
      end do
    end do

#ifdef NETCDF
    call write_netcdf( lognormal_sample_file )
#else
    stop "This version of CLUBB was not compiled for netCDF output"
#endif

    do j = 1, d_variables
      deallocate( lognormal_sample_file%var(j)%ptr )
    end do

    return
  end subroutine output_2D_lognormal_dist_file

!-------------------------------------------------------------------------------
  subroutine output_2D_uniform_dist_file &
             ( nz, num_samples, dp2, X_u_all_levs, X_mixt_comp_all_levs, &
               lh_sample_point_weights )
! Description:
!   Output a 2D snapshot of latin hypercube uniform distribution, i.e. (0,1)
! References:
!   None
!-------------------------------------------------------------------------------
#ifdef NETCDF
    use output_netcdf, only: write_netcdf ! Procedure(s)
#endif

    use clubb_precision, only: &
      core_rknd, &          ! Precision(s)
      stat_rknd

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of vertical levels
      num_samples, & ! Number of samples per variable
      dp2            ! Number of variates being sampled + 2

    real(kind=core_rknd), intent(in), dimension(nz,num_samples,dp2) :: &
      X_u_all_levs ! Uniformly distributed numbers between (0,1)

    integer, intent(in), dimension(nz,num_samples) :: &
      X_mixt_comp_all_levs ! Either 1 or 2

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights ! Weight of each sample

    integer :: sample, j, k

    ! ---- Begin Code ----

    do j = 1, dp2+2
      allocate( uniform_sample_file%var(j)%ptr(num_samples,1,nz) )
    end do

    do sample = 1, num_samples
      do j = 1, dp2
        uniform_sample_file%var(j)%ptr(sample,1,1:nz) = &
          real( X_u_all_levs(1:nz,sample,j), kind = stat_rknd )
      end do
      uniform_sample_file%var(dp2+1)%ptr(sample,1,1:nz) = &
        real( X_mixt_comp_all_levs(1:nz,sample), kind=stat_rknd )
      do k = 1, nz 
        uniform_sample_file%var(dp2+2)%ptr(sample,1,k) = &
          real( lh_sample_point_weights(sample), kind=stat_rknd )
      end do
    end do

#ifdef NETCDF
    call write_netcdf( uniform_sample_file )
#else
    stop "This version of CLUBB was not compiled for netCDF output"
#endif

    do j = 1, dp2+2
      deallocate( uniform_sample_file%var(j)%ptr )
    end do

    return
  end subroutine output_2D_uniform_dist_file

!-------------------------------------------------------------------------------
  subroutine close_2D_samples_file( sample_file )
! Description:
!   Close a 2D sample file
! References:
!   None
!-------------------------------------------------------------------------------
#ifdef NETCDF
    use output_netcdf, only: close_netcdf ! Procedure
#endif

    implicit none

    type(stat_file), intent(inout) :: &
      sample_file ! File that is being opened

    ! ---- Begin Code ----

#ifdef NETCDF
    call close_netcdf( sample_file )
#else
    stop "This version of CLUBB was not compiled for netCDF output"
#endif

    deallocate( sample_file%rlat, sample_file%rlon )
    deallocate( sample_file%var)
    deallocate( sample_file%z)

    return
  end subroutine close_2D_samples_file

end module output_2D_samples_module
