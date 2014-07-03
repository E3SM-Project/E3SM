module spmd_dyn

  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: SPMD implementation of CAM SE finite element dynamics.
  ! 
  !-----------------------------------------------------------------------

  use parallel_mod, only: initmp,  par
  use abortutils,   only: endrun
  use spmd_utils,   only: masterproc, npes

  implicit none
  private

  public spmd_readnl, spmdbuf

  ! These variables are not used locally, but are set and used in phys_grid.
  ! They probably should be moved there.
  logical, public :: local_dp_map=.true.  ! flag indicates that mapping between dynamics 
                                          ! and physics decompositions does not require 
                                          ! interprocess communication
  integer, public :: block_buf_nrecs      ! number of local grid points (lon,lat,lev)
                                          !  in dynamics decomposition (including level 0)
  integer, public :: chunk_buf_nrecs      ! number of local grid points (lon,lat,lev)
                                          !  in physics decomposition (including level 0)
                                          ! assigned in phys_grid.F90

!========================================================================
CONTAINS
!========================================================================

  subroutine spmd_readnl(nlfilename, dyn_npes)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use cam_logfile,     only: iulog
    use mpishorthand

    implicit none

    character(len=*), intent(in) :: nlfilename
    integer, intent(out) :: dyn_npes
    integer :: ierr           ! error code
    integer :: unitn          ! namelist unit number
    integer :: color, nproc_tmp
    character(len=*), parameter ::  subname = "spmd_readnl"

    logical :: dyn_equi_by_col
    integer :: dyn_alltoall
    integer :: dyn_npes_stride
    integer :: dyn_allgather



!   Note that only dyn_npes is currently used by the SE dycore
    namelist /spmd_dyn_inparm/ dyn_alltoall, &
             dyn_allgather,  &
             dyn_equi_by_col,&
             dyn_npes,       &
             dyn_npes_stride 



    dyn_npes = npes

    if (masterproc) then
       write(iulog,*) 'Read in spmd_dyn_inparm namelist from: ', trim(nlfilename)
       unitn = getunit()
       open( unitn, file=trim(nlfilename), status='old' )
       ! Look for spmd_dyn_inparm group name in the input file.  If found, leave the
       ! file positioned at that namelist group.
       call find_group_name(unitn, 'spmd_dyn_inparm', status=ierr)
       if (ierr == 0) then  ! found spmd_dyn_inparm
          read(unitn, spmd_dyn_inparm, iostat=ierr)  ! read the spmd_dyn_inparm namelist group
          if (ierr /= 0) then
             call endrun( subname//':: namelist read returns an'// &
                  ' error condition for spmd_dyn_inparm' )
          end if
       end if
       if (dyn_npes .lt. 1 .or. dyn_npes .gt. npes) then
          call endrun( subname//':: namelist read returns a'// &
               ' bad value for dyn_npes' )
       endif
       write (iulog,*) 'SE will use ', dyn_npes, '  tasks'
       close( unitn )
       call freeunit( unitn )
    endif

    call mpibcast (dyn_npes,1,mpiint,0,mpicom)

  end subroutine spmd_readnl

!================================================================================================

subroutine spmdbuf
end subroutine spmdbuf

!================================================================================================

end module spmd_dyn
