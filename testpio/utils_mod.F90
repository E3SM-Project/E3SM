module utils_mod

  use pio   ! _EXTERNAL
  use kinds_mod
  implicit none
  private

  public :: WriteHeader

#if 0
  interface WriteMeta
     module procedure WriteMeta_single, &
		      WriteMeta_multiple
  end interface
#endif

contains

subroutine WriteHeader(File,nx,ny,nz,dimid_x,dimid_y,dimid_z)

       type (File_desc_t), intent(inout) :: File
       integer(i4), intent(in) :: nx,ny,nz
       integer(i4), intent(out) :: dimid_x,dimid_y,dimid_z

       integer(i4) :: itmp,iostat

       iostat = PIO_put_att(File,pio_global,'title','Test NetCDF file')
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error writing TITLE to netCDF file'
       endif

       iostat = PIO_put_att(File,pio_global,'ivalue', 4)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error writing iVALUE to netCDF file'
       endif

       iostat = PIO_def_dim(File,'X',nx,dimid_x)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining dimension X for netCDF file'
       endif

       iostat = PIO_def_dim(File,'Y',ny,dimid_y)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining dimension Y for netCDF file'
       endif

       iostat = PIO_def_dim(File,'Z',nz,dimid_z)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining dimension Z for netCDF file'
       endif

end subroutine WriteHeader

subroutine ReadHeader(File,nx,ny,nz,dimid_x,dimid_y,dimid_z)
!! NOTE: This subroutine is screwed up, it should be get, not put (tcraig 2/18/09)

       type (File_desc_t), intent(inout) :: File
       integer(i4), intent(in) :: nx,ny,nz
       integer(i4), intent(out) :: dimid_x,dimid_y,dimid_z

       integer(i4) :: itmp,iostat

       iostat = PIO_put_att(File,pio_global,'title','Test NetCDF file')
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error writing TITLE to netCDF file'
       endif

       itmp = 4
       iostat = PIO_put_att(File,pio_global,'ivalue',(/itmp/))
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error writing iVALUE to netCDF file'
       endif

       iostat = PIO_def_dim(File,'X',nx,dimid_x)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining dimension X for netCDF file'
       endif

       iostat = PIO_def_dim(File,'Y',ny,dimid_y)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining dimension Y for netCDF file'
       endif

       iostat = PIO_def_dim(File,'Z',nz,dimid_z)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining dimension Z for netCDF file'
       endif

end subroutine ReadHeader

#if 0
subroutine WriteMeta_single(File, IODesc, xtype, nx,ny)

       type (File_desc_t), intent(inout) :: File
       type (IO_desc_t), intent(inout) :: IODesc
       integer(i4), intent(in) :: xtype
       integer(i4), intent(in) :: nx,ny

       integer :: iostat,dimid_x,dimid_y,itmp

       call WriteHeader(File,nx,ny,dimid_x,dimid_y)

       iostat = PIO_def_var(File,'myTaskID',xtype,(/dimid_x,dimid_y/),IODesc%varID)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining variable Tvel for netCDF file'
       endif

       iostat = PIO_enddef(File)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error ending define-mode'
       endif

end subroutine WriteMeta_single


subroutine WriteMeta_multiple(File, IODescA, IODescB,IODescC, nx,ny)

       type (File_desc_t), intent(inout) :: File
       type (IO_desc_t), intent(inout) :: IODescA,IODescB,IODescC
       integer(i4), intent(in) :: nx,ny

       integer :: iostat,dimid_x,dimid_y,itmp

       call WriteHeader(File,nx,ny,dimid_x,dimid_y)

       iostat = PIO_def_var(File,'myTaskID_r8',PIO_double,(/dimid_x,dimid_y/),IODescA%varID)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining variable myTaskID_r8 for netCDF file'
       endif

       iostat = PIO_def_var(File,'myTaskID_r4',PIO_real,(/dimid_x,dimid_y/),IODescB%varID)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining variable myTaskID_r4 for netCDF file'
       endif

       iostat = PIO_def_var(File,'myTaskID_i4',PIO_int,(/dimid_x,dimid_y/),IODescC%varID)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error defining variable myTaskID_i4 for netCDF file'
       endif

       iostat = PIO_enddef(File)
       if(iostat /= pio_noerr) then
          write(*,*) 'testPIO:  Error ending define-mode'
       endif

end subroutine WriteMeta_multiple
#endif

end module utils_mod
