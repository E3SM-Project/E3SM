module module_mosaic_support
  !Purpose: This module contains subroutines which have codes which depend upon 
  !         the host code (CAM, WRF etc.). #defines are used to seprate codes
  !         which depends on the host code

  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun

  implicit none
  private

  public:: mosaic_warn_mess
  public:: mosaic_err_mess
  
  
contains

  subroutine mosaic_warn_mess(message)
    !Purpose: Print out the warning messages from Mosaic code

    character(len=*), intent(in) :: message

    !Local variables
    character(len=16), parameter :: warn_str = 'MOSAIC WARNING: ' 

    ! write(iulog,*)warn_str,message!BALLI -comment out to avoid exxcessive warning messages.   

  end subroutine mosaic_warn_mess


  subroutine mosaic_err_mess(message)
    !Purpose
    character(len=*), intent(in) :: message

    !Local variables
    character(len=14), parameter :: err_str = 'MOSAIC ERROR: ' 
    character(len=500) :: str_to_prnt 

    write(str_to_prnt,*)err_str,message
    
    call endrun(trim(adjustl(str_to_prnt)))    

  end subroutine mosaic_err_mess

end module module_mosaic_support
