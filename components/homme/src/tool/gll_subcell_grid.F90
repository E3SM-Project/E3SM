module gll_subcell_grid
  
  implicit none
  private
  public :: write_gll_subcell_grid

contains
  
  subroutine write_gll_subcell_grid(par, elem)
    use parallel_mod,     only: parallel_t, haltmp
    use element_mod, only : element_t
    use pio, only : file_desc_t, pio_def_dim, var_desc_t, pio_int, pio_def_var, &
         pio_enddef, pio_closefile, pio_initdecomp, io_desc_t, pio_write_darray, &
         pio_freedecomp, pio_setdebuglevel, pio_rearr_box, PIO_CLOBBER, iotype_netcdf, &
         pio_init, pio_createfile
    use dimensions_mod, only : ne, np, nelem, nelemd
    use dof_mod, only : createmetadata
    use common_io_mod,    only: infilenames, PIOFS,io_stride, num_io_procs, num_agg, output_dir
    
    type(parallel_t) :: par
    type(element_t) :: elem(:)
    type(file_desc_t) :: nc
    type(var_desc_t) :: vid
    type(io_desc_t) :: iodesc
    integer :: dim1, dim2, ierr, i, j, ie, cc, base, ii, jj
    integer, parameter :: npm12 = (np-1)*(np-1)
    integer :: subelement_corners(npm12*nelemd,4)
    integer :: dof(npm12*nelemd*4)
    character(len=12) :: nestr,npstr
    character(len=240) :: mappingfile

    ! Create a CS grid mapping file for postprocessing tools

    ! write meta data for physics on GLL nodes
    if(.not.par%dynproc) then
       ! The special case of npes_se < npes_cam is not worth dealing with here
       call haltmp('Native mapping code requires npes_se==npes_cam')
    end if
    call PIO_Init(par%rank, par%comm, num_io_procs, num_agg, &
         io_stride, pio_rearr_box, PIOFS)

    write(nestr,*) ne
    write(npstr,*) np
    mappingfile = trim(output_dir) // &
         'gll_subcell_grid_ne' // trim(adjustl(nestr)) // 'np' // trim(adjustl(npstr)) // ".nc"
    
    ierr = pio_createfile(PIOFS, nc, iotype_netcdf, mappingfile,PIO_CLOBBER)

    ierr = pio_def_dim(nc, 'ncenters', npm12*nelem, dim1)
    ierr = pio_def_dim(nc, 'ncorners', 4, dim2)
    ierr = pio_def_var(nc, 'element_corners', PIO_INT, (/dim1,dim2/),vid)

    ierr = pio_enddef(nc)
    if (par%dynproc) then
       call createmetadata(par, elem, subelement_corners)
    end if

    jj=0
    do cc=0,3
       do ie=1,nelemd
          base = ((elem(ie)%globalid-1)+cc*nelem)*npm12
          ii=0
          do j=1,np-1
             do i=1,np-1
                ii=ii+1
                jj=jj+1
                dof(jj) = base+ii
             end do
          end do
       end do
    end do

    call pio_initdecomp(nc%iosystem, pio_int, (/nelem*npm12,4/), dof, iodesc)

    call pio_write_darray(nc, vid, iodesc, reshape(subelement_corners,(/nelemd*npm12*4/)), ierr)

    call pio_freedecomp(nc, iodesc)

    call pio_closefile(nc)

  end subroutine
end module

