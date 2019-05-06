module restart_physics

  use shr_kind_mod,       only: r8 => shr_kind_r8
  use spmd_utils,         only: masterproc
  use constituents,       only: pcnst

  use cam_abortutils,     only: endrun
  use camsrfexch,         only: cam_in_t, cam_out_t
  use cam_logfile,        only: iulog
  use pio,                only: file_desc_t, io_desc_t, var_desc_t, &
                                pio_double, pio_int, pio_noerr, &
                                pio_seterrorhandling, pio_bcast_error, &
                                pio_inq_varid, &
                                pio_def_var, pio_def_dim, &
                                pio_put_var, pio_get_var

  implicit none
  private
  save
!
! Public interfaces
!
  public :: write_restart_physics    ! Write the physics restart info out
  public :: read_restart_physics     ! Read the physics restart info in
  public :: get_abs_restart_filepath ! Get the name of the restart filepath
  public :: init_restart_physics

!
! Private data
!
  character(len=256) :: pname  ! Full abs-ems restart filepath

  CONTAINS
    subroutine init_restart_physics ( File, pbuf2d)
      
    use physics_buffer,      only: pbuf_init_restart, physics_buffer_desc
    use cam_grid_support,    only: cam_grid_write_attr, cam_grid_id
    use cam_grid_support,    only: cam_grid_header_info_t

    type(file_desc_t), intent(inout) :: file
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer                      :: grid_id
    integer                      :: hdimcnt, ierr, i, vsize
    integer                      :: dimids(4)
    integer, allocatable         :: hdimids(:)
    integer                      :: ndims, pver_id, pverp_id
    integer                      :: kiss_seed_dim

    type(cam_grid_header_info_t) :: info

    call pio_seterrorhandling(File, PIO_BCAST_ERROR)
    ! Probably should have the grid write this out.
    grid_id = cam_grid_id('physgrid')
    call cam_grid_write_attr(File, grid_id, info)
    hdimcnt = info%num_hdims()

    do i = 1, hdimcnt
      dimids(i) = info%get_hdimid(i)
    end do
    allocate(hdimids(hdimcnt))
    hdimids(1:hdimcnt) = dimids(1:hdimcnt)

    call pbuf_init_restart(File, pbuf2d)
      
  end subroutine init_restart_physics

  subroutine write_restart_physics (File, cam_in, cam_out, pbuf2d)

      !-----------------------------------------------------------------------
      use physics_buffer,      only: physics_buffer_desc, pbuf_write_restart
      use phys_grid,           only: phys_decomp
      use ppgrid,              only: begchunk, endchunk, pcols, pverp
      use cam_grid_support,    only: cam_grid_write_var
      !
      ! Input arguments
      !
      type(file_desc_t), intent(inout) :: File
      type(cam_in_t),    intent(in)    :: cam_in(begchunk:endchunk)
      type(cam_out_t),   intent(in)    :: cam_out(begchunk:endchunk)
      type(physics_buffer_desc), pointer        :: pbuf2d(:,:)
      !-----------------------------------------------------------------------

      ! Write grid vars
      call cam_grid_write_var(File, phys_decomp)

      call pbuf_write_restart(File, pbuf2d)
      
    end subroutine write_restart_physics

!#######################################################################

    subroutine read_restart_physics(File, cam_in, cam_out, pbuf2d)

     !-----------------------------------------------------------------------
     use physics_buffer,      only: physics_buffer_desc, pbuf_read_restart
     
     use ppgrid,              only: begchunk, endchunk, pcols, pver, pverp
     use cam_grid_support,    only: cam_grid_read_dist_array, cam_grid_id
     use cam_grid_support,    only: cam_grid_get_decomp, cam_grid_dimensions
     use cam_history_support, only: fillvalue
     use pio,                 only: pio_read_darray
     !
     ! Arguments
     !
     type(file_desc_t),   intent(inout) :: File
     type(cam_in_t),            pointer :: cam_in(:)
     type(cam_out_t),           pointer :: cam_out(:)
     type(physics_buffer_desc), pointer :: pbuf2d(:,:)
     !-----------------------------------------------------------------------

     call pbuf_read_restart(File, pbuf2d)

   end subroutine read_restart_physics


   character(len=256) function get_abs_restart_filepath ( )
     !	
     ! Return the full filepath to the abs-ems restart file
     !	
     get_abs_restart_filepath = pname
   end function get_abs_restart_filepath

 end module restart_physics
