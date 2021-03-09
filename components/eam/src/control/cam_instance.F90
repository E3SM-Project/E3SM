module cam_instance

use seq_comm_mct, only: seq_comm_suffix, seq_comm_inst, seq_comm_name

implicit none
private
save

public :: cam_instance_init

integer,           public :: atm_id
integer,           public :: inst_index
character(len=16), public :: inst_name
character(len=16), public :: inst_suffix

!===============================================================================
CONTAINS
!===============================================================================

subroutine cam_instance_init(in_atm_id)

   integer, intent(in) :: in_atm_id

   atm_id      = in_atm_id
   inst_name   = seq_comm_name(atm_id)
   inst_index  = seq_comm_inst(atm_id)
   inst_suffix = seq_comm_suffix(atm_id)

end subroutine cam_instance_init

end module cam_instance
