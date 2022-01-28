!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module glc_indexing_info
  
  !BOP
  ! !MODULE: glc_indexing_info

  ! !DESCRIPTION:
  ! Contains information about the indexing of the points owned by each processor.
  !
  ! This includes local indices (translation between (i,j) and a scalar 1..n) and global
  ! indices (unique indices across all procs).

  implicit none
  private
  save

  ! !PUBLIC ROUTINES:
  public :: glc_indexing_init

  ! !PUBLIC MODULE VARIABLES:
  integer, public :: nx       ! number of columns owned by this proc
  integer, public :: ny       ! number of rows owned by this proc
  integer, public :: npts     ! total number of points owned by this proc
  integer, public :: nx_tot   ! total number of columns in full grid (all procs)
  integer, public :: ny_tot   ! total number of rows in full grid (all procs)
  integer, public :: npts_tot ! total number of points in full grid (all procs)

  integer, allocatable, public :: local_indices(:,:)  ! mapping from (i,j) to 1..npts
  integer, allocatable, public :: global_indices(:,:) ! unique indices across all procs (matches indexing on mapping files)
  
contains

  !-----------------------------------------------------------------------
  subroutine glc_indexing_init(params, instance_index)
    !
    ! !DESCRIPTION:
    ! Initialize indices stored here.
    !
    ! Note that the global indexing needs to match the indexing on the SCRIP grid file
    ! that is used to generate GLC mapping files for the coupler.
    !
    ! !USES:
    use glad_main, only : glad_params, glad_get_grid_size, glad_get_grid_indices
    !
    ! !ARGUMENTS:
    type(glad_params), intent(in) :: params
    integer, intent(in) :: instance_index  ! index of current ice sheet index
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'glc_indexing_init'
    !-----------------------------------------------------------------------

    call glad_get_grid_size(params, instance_index, &
         ewn = nx, nsn = ny, npts = npts, &
         ewn_tot = nx_tot, nsn_tot = ny_tot, npts_tot = npts_tot)
    
    allocate(local_indices(nx, ny))
    allocate(global_indices(nx, ny))
    
    call glad_get_grid_indices(params, instance_index, global_indices, local_indices)
    
  end subroutine glc_indexing_init

end module glc_indexing_info
