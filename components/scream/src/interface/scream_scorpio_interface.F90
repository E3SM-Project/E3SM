module scream_scorpio_interface
  
  use pio_mods, only: io_desc_t, iosystem_desc_t, file_desc_t, var_desc_t, pio_global
                     
  public :: eam_init_pio_subsystem, &  ! Get pio subsystem info from main code
            eam_h_define, &      ! Create a new NetCDF file for PIO writing
            eam_h_finalize

  
  integer               :: pio_iotype
  type(file_desc_t)     :: pioFile
  type(iosystem_desc_t), pointer, public :: pio_subsystem => null()

contains

!=====================================================================!
  subroutine eam_h_define()

   character(len=100) :: fname

   fname = "eam_pio_example.nc"

   ! Create the file
   call eam_pio_createfile(pioFile,trim(fname)) ! TODO set up File and fname inputs
   ! Create netCDF Header info (like caseid, title, etc.)
   call eam_pio_createHeader(pioFile)
   ! Define dimensions
!   call eam_pio_define_basic_dimensions(pioFile) ! TODO set up basic_vars definition
   ! TODO define vars from field manager
   ! TODO define optional dimensions for nonstandard variables.  Alternatively,
   ! define dimensions based on field manager (second option is probably better).
   ! Define Grid attribute
   ! TODO Define grid attribute routine ala cam_grid_write_attr


  end subroutine eam_h_define
!=====================================================================!
  subroutine eam_h_finalize()

    call eam_pio_finalize()

  end subroutine eam_h_finalize
!=====================================================================!
  subroutine eam_pio_createHeader(File)

    use pio_mods, only : PIO_put_att

    type(file_desc_t), intent(in) :: File             ! Pio file Handle
    integer :: retval

    ! TODO change options below to match specific simulation case
    retval=pio_put_att (File, PIO_GLOBAL, 'source', 'SCREAM')
    retval=pio_put_att (File, PIO_GLOBAL, 'case', 'TEST')
    retval=pio_put_att (File, PIO_GLOBAL, 'title', 'SCORPIO TEST')
    retval=pio_put_att (File, PIO_GLOBAL, 'logname','THE GIT LOG HASH')
    retval=pio_put_att (File, PIO_GLOBAL, 'host', 'THE HOST')
    retval=pio_put_att (File, PIO_GLOBAL, 'Version', &
           '0')
    retval=pio_put_att (File, PIO_GLOBAL, 'revision_Id', &
           'None')
    retval=pio_put_att (File, PIO_GLOBAL, 'initial_file', 'NONE FOR NOW')
    retval=pio_put_att (File, PIO_GLOBAL, 'topography_file', 'NONE FOR NOW')
    
  end subroutine eam_pio_createHeader
!=====================================================================!
  subroutine eam_init_pio_subsystem(atm_id)
    use shr_pio_mod,   only: shr_pio_getiosys, shr_pio_getiotype
    
    integer, intent(in) :: atm_id

    pio_subsystem => shr_pio_getiosys(atm_id)
    pio_iotype =  shr_pio_getiotype(atm_id)

  end subroutine eam_init_pio_subsystem
!=====================================================================!
  subroutine eam_pio_createfile(File,fname)
    use pio_mods, only:  pio_createfile, pio_clobber

    type(file_desc_t), intent(inout) :: File             ! Pio file Handle
    character(len=*),  intent(in)    :: fname
    !--
    integer :: retval
    integer                                   :: mode

    mode = pio_clobber ! Set to CLOBBER for now, TODO: fix to allow for optional mode type like in CAM
    retval = pio_createfile(pio_subsystem,File,pio_iotype,fname,mode) 

  end subroutine eam_pio_createfile
!=====================================================================!
  subroutine eam_pio_finalize()

    use pio_mods, only: pio_finalize

    integer :: ierr

    call PIO_finalize(pio_subsystem, ierr)

  end subroutine eam_pio_finalize
!=====================================================================!

end module scream_scorpio_interface
