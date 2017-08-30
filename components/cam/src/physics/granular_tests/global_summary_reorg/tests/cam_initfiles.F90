module cam_initfiles
!----------------------------------------------------------------------- 
! 
! Open, close, and provide access to the initial conditions and topography files.
! 
!-----------------------------------------------------------------------

implicit none
private
save

! Public methods

public :: &
   cam_initfiles_open,   &! open initial and topo files
   initial_file_get_id,  &! returns filehandle for initial file
   cam_initfiles_close     ! close initial and topo files

   integer, save :: fh_ini, fh_topo

!======================================================================= 
contains
!======================================================================= 

integer function initial_file_get_id()
  fh_ini = 12
  initial_file_get_id = fh_ini
end function initial_file_get_id

!======================================================================= 

subroutine cam_initfiles_open()

   ! Open the initial conditions and topography files.

!!$   use cam_pio_utils,    only: cam_pio_openfile
!!$   use pio,              only: pio_nowrite

!!$   use readinitial,      only: read_initial

   character(len=256) :: ncdata_loc     ! filepath of initial file on local disk
   !----------------------------------------------------------------------- 
   
   ! Open initial, topography, and landfrac datasets
!!$   call getfil (ncdata, ncdata_loc)
   open(file="state_snapshot",unit=fh_ini)
!!$   call cam_pio_openfile(fh_ini, ncdata_loc, PIO_NOWRITE)
   ! Backward compatibility: look for topography data on initial file if topo file name not provided.
!!$   if (trim(bnd_topo) /= 'bnd_topo' .and. len_trim(bnd_topo) > 0) then
!!$      allocate(fh_topo)
!!$      call getfil(bnd_topo, bnd_topo_loc)
!!$      call cam_pio_openfile(fh_topo, bnd_topo_loc, PIO_NOWRITE)
!!$   else
!!$      fh_topo => fh_ini
!!$   end if

   ! Check for consistent settings on initial dataset -- this is dycore
   ! dependent -- should move to dycore interface
!!$   call read_initial (fh_ini)

end subroutine cam_initfiles_open

!======================================================================= 

subroutine cam_initfiles_close()

!!$  use pio,          only: pio_closefile

!!$  if(associated(fh_ini)) then
!!$     if(.not. associated(fh_ini, target=fh_topo)) then
!!$        call pio_closefile(fh_topo)
!!$        deallocate(fh_topo)
!!$     end if
!!$     
!!$     call pio_closefile(fh_ini)
!!$     deallocate(fh_ini)
!!$     nullify(fh_ini)
!!$     nullify(fh_topo)
!!$  end if
  close(fh_ini)
end subroutine cam_initfiles_close

!======================================================================= 

end module cam_initfiles
