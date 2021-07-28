module elm_instance

use seq_comm_mct, only: seq_comm_suffix, seq_comm_inst, seq_comm_name

implicit none
private
save

public :: elm_instance_init

integer,           public :: lnd_id      ! land component identifier for the current instance
integer,           public :: inst_index
character(len=16), public :: inst_name
character(len=16), public :: inst_suffix

!===============================================================================
CONTAINS
!===============================================================================

subroutine elm_instance_init(in_lnd_id)

   integer, intent(in) :: in_lnd_id

   lnd_id      = in_lnd_id
   inst_name   = seq_comm_name(lnd_id)
   inst_index  = seq_comm_inst(lnd_id)
   inst_suffix = seq_comm_suffix(lnd_id)

end subroutine elm_instance_init

end module elm_instance
