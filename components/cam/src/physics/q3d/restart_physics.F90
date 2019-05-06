module restart_physics

  use shr_kind_mod,       only: r8 => shr_kind_r8
  use spmd_utils,         only: masterproc
  use co2_cycle,          only: co2_transport
  use constituents,       only: pcnst

!  use radiation,          only: radiation_define_restart
!  use radiation,          only: radiation_write_restart
!  use radiation,          only: radiation_read_restart

  use ioFileMod
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
  public :: init_restart_physics
  public :: get_abs_restart_filepath

!
! Private data
!
    character(len=256) :: pname  ! Full abs-ems restart filepath
    type(var_desc_t) :: flwds_desc, &
         solld_desc, co2prog_desc, co2diag_desc, sols_desc, soll_desc, &
         solsd_desc

    type(var_desc_t) :: bcphidry_desc, bcphodry_desc, ocphidry_desc, ocphodry_desc, &
       dstdry1_desc, dstdry2_desc, dstdry3_desc, dstdry4_desc

    type(var_desc_t) :: cflx_desc(pcnst)

    type(var_desc_t) :: wsx_desc
    type(var_desc_t) :: wsy_desc
    type(var_desc_t) :: shf_desc

  CONTAINS
    subroutine init_restart_physics ( File, pbuf2d)

    use physics_buffer,      only: pbuf_init_restart, physics_buffer_desc
    use q3d_restart,         only: init_q3d_restart
    use cam_grid_support,    only: cam_grid_write_attr, cam_grid_id
    use cam_grid_support,    only: cam_grid_header_info_t
    use cam_pio_utils,       only: cam_pio_def_dim
    use subcol_utils,        only: is_subcol_on
    use subcol,              only: subcol_init_restart

    type(file_desc_t), intent(inout) :: file
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer                      :: grid_id
    integer                      :: hdimcnt, ierr, i
    integer                      :: dimids(4)
    integer, allocatable         :: hdimids(:)
    type(cam_grid_header_info_t) :: info
    character(len=4)   :: num

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

    call init_q3d_restart(File)

    ierr = pio_def_var(File, 'FLWDS', pio_double, hdimids, flwds_desc)
    ierr = pio_def_var(File, 'SOLS', pio_double, hdimids, sols_desc)
    ierr = pio_def_var(File, 'SOLL', pio_double, hdimids, soll_desc)
    ierr = pio_def_var(File, 'SOLSD', pio_double, hdimids, solsd_desc)
    ierr = pio_def_var(File, 'SOLLD', pio_double, hdimids, solld_desc)

    ierr = pio_def_var(File, 'BCPHIDRY', pio_double, hdimids, bcphidry_desc)
    ierr = pio_def_var(File, 'BCPHODRY', pio_double, hdimids, bcphodry_desc)
    ierr = pio_def_var(File, 'OCPHIDRY', pio_double, hdimids, ocphidry_desc)
    ierr = pio_def_var(File, 'OCPHODRY', pio_double, hdimids, ocphodry_desc)
    ierr = pio_def_var(File, 'DSTDRY1',  pio_double, hdimids, dstdry1_desc)
    ierr = pio_def_var(File, 'DSTDRY2',  pio_double, hdimids, dstdry2_desc)
    ierr = pio_def_var(File, 'DSTDRY3',  pio_double, hdimids, dstdry3_desc)
    ierr = pio_def_var(File, 'DSTDRY4',  pio_double, hdimids, dstdry4_desc)

    if(co2_transport()) then
      ierr = pio_def_var(File, 'CO2PROG', pio_double, hdimids, co2prog_desc)
      ierr = pio_def_var(File, 'CO2DIAG', pio_double, hdimids, co2diag_desc)
    end if

    ! cam_import variables -- write the constituent surface fluxes as individual 2D arrays
    ! rather than as a single variable with a pcnst dimension.  Note that the cflx components
    ! are only needed for those constituents that are not passed to the coupler.  The restart
    ! for constituents passed through the coupler are handled by the .rs. restart file.  But
    ! we don't currently have a mechanism to know whether the constituent is handled by the
    ! coupler or not, so we write all of cflx to the CAM restart file.
    do i = 1, pcnst
      write(num,'(i4.4)') i
      ierr = pio_def_var(File, 'CFLX'//num,  pio_double, hdimids, cflx_desc(i))
    end do

    ierr = pio_def_var(File, 'wsx',  pio_double, hdimids, wsx_desc)
    ierr = pio_def_var(File, 'wsy',  pio_double, hdimids, wsy_desc)
    ierr = pio_def_var(File, 'shf',  pio_double, hdimids, shf_desc)

!    call radiation_define_restart(file)

    if (is_subcol_on()) then
      call subcol_init_restart(file, hdimids)
    end if

  end subroutine init_restart_physics

  subroutine write_restart_physics (File, cam_in, cam_out, pbuf2d)

      !-----------------------------------------------------------------------
      use physics_buffer,      only: physics_buffer_desc, pbuf_write_restart
      use phys_grid,           only: phys_decomp

      use ppgrid,              only: begchunk, endchunk, pcols
      use q3d_restart,         only: write_q3d_restart

      use cam_history_support, only: fillvalue
      use spmd_utils,          only: iam
      use cam_grid_support,    only: cam_grid_write_dist_array, cam_grid_id
      use cam_grid_support,    only: cam_grid_get_decomp, cam_grid_dimensions
      use cam_grid_support,    only: cam_grid_write_var
      use pio,                 only: pio_write_darray
      use subcol_utils,        only: is_subcol_on
      use subcol,              only: subcol_write_restart
      !
      ! Input arguments
      !
      type(file_desc_t), intent(inout)   :: File
      type(cam_in_t),    intent(in)      :: cam_in(begchunk:endchunk)
      type(cam_out_t),   intent(in)      :: cam_out(begchunk:endchunk)
      type(physics_buffer_desc), pointer :: pbuf2d(:,:)
      !
      ! Local workspace
      !
      type(io_desc_t), pointer :: iodesc
      real(r8):: tmpfield(pcols, begchunk:endchunk)
      integer :: i, m          ! loop index
      integer :: ncol          ! number of vertical columns
      integer :: ierr
      integer :: physgrid
      integer :: dims(3), gdims(3)
      integer :: nhdims
      !-----------------------------------------------------------------------

      ! Write grid vars
      call cam_grid_write_var(File, phys_decomp)

      if (is_subcol_on()) then
         call subcol_write_restart(File)
      end if
      ! Physics buffer
      call pbuf_write_restart(File, pbuf2d)

      physgrid = cam_grid_id('physgrid')
      call cam_grid_dimensions(physgrid, gdims(1:2), nhdims)

      call write_q3d_restart(File)

      ! cam_in/out variables
      ! This is a group of surface variables so can reuse dims
      dims(1) = pcols
      dims(2) = endchunk - begchunk + 1
      call cam_grid_get_decomp(physgrid, dims(1:2), gdims(1:nhdims), &
           pio_double, iodesc)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%flwds(:ncol)
        ! Only have to do this once (cam_in/out vars all same shape)
        if (ncol < pcols) then
          tmpfield(ncol+1:, i) = fillvalue
        end if
      end do
      call pio_write_darray(File, flwds_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%sols(:ncol)
      end do
      call pio_write_darray(File, sols_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%soll(:ncol)
      end do
      call pio_write_darray(File, soll_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%solsd(:ncol)
      end do
      call pio_write_darray(File, solsd_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%solld(:ncol)
      end do
      call pio_write_darray(File, solld_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%bcphidry(:ncol)
      end do
      call pio_write_darray(File, bcphidry_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%bcphodry(:ncol)
      end do
      call pio_write_darray(File, bcphodry_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%ocphidry(:ncol)
      end do
      call pio_write_darray(File, ocphidry_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%ocphodry(:ncol)
      end do
      call pio_write_darray(File, ocphodry_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%dstdry1(:ncol)
      end do
      call pio_write_darray(File, dstdry1_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%dstdry2(:ncol)
      end do
      call pio_write_darray(File, dstdry2_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%dstdry3(:ncol)
      end do
      call pio_write_darray(File, dstdry3_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_out(i)%ncol
        tmpfield(:ncol, i) = cam_out(i)%dstdry4(:ncol)
      end do
      call pio_write_darray(File, dstdry4_desc, iodesc, tmpfield, ierr)

      if (co2_transport()) then
        do i = begchunk, endchunk
          ncol = cam_out(i)%ncol
          tmpfield(:ncol, i) = cam_out(i)%co2prog(:ncol)
        end do
        call pio_write_darray(File, co2prog_desc, iodesc, tmpfield, ierr)

        do i = begchunk, endchunk
          ncol = cam_out(i)%ncol
          tmpfield(:ncol, i) = cam_out(i)%co2diag(:ncol)
        end do
        call pio_write_darray(File, co2diag_desc, iodesc, tmpfield, ierr)
      end if

      ! cam_in components
      do m = 1, pcnst
        do i = begchunk, endchunk
          ncol = cam_in(i)%ncol
          tmpfield(:ncol, i) = cam_in(i)%cflx(:ncol, m)
        end do
        call pio_write_darray(File, cflx_desc(m), iodesc, tmpfield, ierr)
      end do

      do i = begchunk, endchunk
        ncol = cam_in(i)%ncol
        tmpfield(:ncol,i) = cam_in(i)%wsx(:ncol)
      end do
      call pio_write_darray(File, wsx_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_in(i)%ncol
        tmpfield(:ncol,i) = cam_in(i)%wsy(:ncol)
      end do
      call pio_write_darray(File, wsy_desc, iodesc, tmpfield, ierr)

      do i = begchunk, endchunk
        ncol = cam_in(i)%ncol
        tmpfield(:ncol,i) = cam_in(i)%shf(:ncol)
      end do
      call pio_write_darray(File, shf_desc, iodesc, tmpfield, ierr)

!      call radiation_write_restart(file)

    end subroutine write_restart_physics

!#######################################################################

    subroutine read_restart_physics(File, cam_in, cam_out, pbuf2d)

     !-----------------------------------------------------------------------
     use physics_buffer,      only: physics_buffer_desc, pbuf_read_restart

     use ppgrid,              only: begchunk, endchunk, pcols
     use cam_grid_support,    only: cam_grid_read_dist_array, cam_grid_id
     use cam_grid_support,    only: cam_grid_get_decomp, cam_grid_dimensions
     use cam_history_support, only: fillvalue

     use q3d_restart,         only: read_q3d_restart
     use subcol_utils,        only: is_subcol_on
     use subcol,              only: subcol_read_restart
     use pio,                 only: pio_read_darray
     !
     ! Arguments
     !
     type(file_desc_t), intent(inout)   :: File
     type(cam_in_t),            pointer :: cam_in(:)
     type(cam_out_t),           pointer :: cam_out(:)
     type(physics_buffer_desc), pointer :: pbuf2d(:,:)
     !
     ! Local workspace
     !
     real(r8), allocatable :: tmpfield2(:,:)
     integer :: i, c, m           ! loop index
     integer :: ierr             ! I/O status
     type(io_desc_t), pointer :: iodesc
     type(var_desc_t)         :: vardesc
     integer                  :: csize, vsize
     character(len=4)         :: num
     integer                  :: dims(3), gdims(3), nhdims
     integer                  :: err_handling
     integer                  :: physgrid
     !-----------------------------------------------------------------------

     ! subcol_read_restart must be called before pbuf_read_restart
     if (is_subcol_on()) then
        call subcol_read_restart(File)
     end if

     call pbuf_read_restart(File, pbuf2d)

     csize=endchunk-begchunk+1
     dims(1) = pcols
     dims(2) = csize

     physgrid = cam_grid_id('physgrid')

     call cam_grid_dimensions(physgrid, gdims(1:2))

     if (gdims(2) == 1) then
       nhdims = 1
     else
       nhdims = 2
     end if
     call cam_grid_get_decomp(physgrid, dims(1:2), gdims(1:nhdims), pio_double, &
          iodesc)

     call read_q3d_restart(File)

     allocate(tmpfield2(pcols, begchunk:endchunk))
     tmpfield2 = fillvalue

     ierr = pio_inq_varid(File, 'FLWDS', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%flwds(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'SOLS', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%sols(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'SOLL', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%soll(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'SOLSD', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%solsd(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'SOLLD', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%solld(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'BCPHIDRY', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%bcphidry(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'BCPHODRY', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%bcphodry(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'OCPHIDRY', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%ocphidry(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'OCPHODRY', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%ocphodry(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'DSTDRY1', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%dstdry1(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'DSTDRY2', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%dstdry2(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'DSTDRY3', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%dstdry3(i) = tmpfield2(i, c)
       end do
     end do

     ierr = pio_inq_varid(File, 'DSTDRY4', vardesc)
     call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
     do c=begchunk,endchunk
       do i=1,pcols
         cam_out(c)%dstdry4(i) = tmpfield2(i, c)
       end do
     end do

     if (co2_transport()) then
       ierr = pio_inq_varid(File, 'CO2PROG', vardesc)
       call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
       do c=begchunk,endchunk
         do i=1,pcols
           cam_out(c)%co2prog(i) = tmpfield2(i, c)
         end do
       end do

       ierr = pio_inq_varid(File, 'CO2DIAG', vardesc)
       call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
       do c=begchunk,endchunk
         do i=1,pcols
           cam_out(c)%co2diag(i) = tmpfield2(i, c)
         end do
       end do
     end if

     ! Reading the CFLX* components from the restart is optional for
     ! backwards compatibility.  These fields were not needed for an
     ! exact restart until the UNICON scheme was added.  More generally,
     ! these components are only needed if they are not handled by the
     ! coupling layer restart (the ".rs." file), and if the values are
     ! used in the tphysbc physics before the tphysac code has a chance
     ! to update the values that are coming from boundary datasets.
     do m = 1, pcnst

       write(num,'(i4.4)') m

       call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)
       ierr = pio_inq_varid(File, 'CFLX'//num, vardesc)
       call pio_seterrorhandling(File, err_handling)

       if (ierr == PIO_NOERR) then ! CFLX variable found on restart file
         call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
         do c= begchunk, endchunk
           do i = 1, pcols
             cam_in(c)%cflx(i,m) = tmpfield2(i, c)
           end do
         end do
       end if

     end do

     call pio_seterrorhandling(File, PIO_BCAST_ERROR, err_handling)
     ierr = pio_inq_varid(File, 'wsx', vardesc)
     if (ierr == PIO_NOERR) then ! variable found on restart file
       call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
       do c= begchunk, endchunk
         do i = 1, pcols
           cam_in(c)%wsx(i) = tmpfield2(i, c)
         end do
       end do
     end if
     ierr = pio_inq_varid(File, 'wsy', vardesc)
     if (ierr == PIO_NOERR) then ! variable found on restart file
       call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
       do c= begchunk, endchunk
         do i = 1, pcols
           cam_in(c)%wsy(i) = tmpfield2(i, c)
         end do
       end do
     end if
     ierr = pio_inq_varid(File, 'shf', vardesc)
     if (ierr == PIO_NOERR) then ! variable found on restart file
       call pio_read_darray(File, vardesc, iodesc, tmpfield2, ierr)
       do c= begchunk, endchunk
         do i = 1, pcols
           cam_in(c)%shf(i) = tmpfield2(i, c)
         end do
       end do
     endif
     call pio_seterrorhandling(File, err_handling)

     deallocate(tmpfield2)

!     call radiation_read_restart(file)

   end subroutine read_restart_physics

   character(len=256) function get_abs_restart_filepath ( )
     !  
     ! Return the full filepath to the abs-ems restart file
     !  
     get_abs_restart_filepath = pname
   end function get_abs_restart_filepath

 end module restart_physics
