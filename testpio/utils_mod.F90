module utils_mod

  use pio   ! _EXTERNAL
  use kinds_mod
  implicit none
  private

  public :: WriteHeader

contains

!>
!! @private
!! @brief Writes netcdf header information for testpio. 
!! @param File @copydoc file_desc_t
!! @param nx
!! @param ny
!! @param nz
!! @param dimid_x
!! @param dimid_y
!! @param dimid_z
!<
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

end module utils_mod
