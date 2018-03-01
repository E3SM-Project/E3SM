module mpp_shr_log_mod

  implicit none

  character(len=*), parameter :: mod_filename = &
       __FILE__

  private

  ! !PUBLIC TYPES:

  ! no public types

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: shr_log_errMsg

contains

  function shr_log_errMsg(file, line)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(len=512)   :: shr_log_errMsg
    character(len=*), intent(in) :: file
    integer         , intent(in) :: line

    !EOP

    write(shr_log_errMsg, '(a, a, a, i0)') 'ERROR in ', trim(file), ' at line ', line

  end function shr_log_errMsg
end module mpp_shr_log_mod
