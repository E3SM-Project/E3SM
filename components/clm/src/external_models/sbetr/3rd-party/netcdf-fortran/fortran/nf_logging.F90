#ifdef LOGGING
! Function to turn on logging
!-------------------------------- nf_set_log_level ----------------------------
 Function nf_set_log_level(new_level) Result(status)

 USE ISO_C_BINDING, ONLY: C_INT

 Implicit NONE

 Integer, Intent(IN) :: new_level

 Integer             :: status

 Integer(C_INT) :: cnew_level, cstatus

 Interface  ! define binding here instead of nc_interfaces since its conditional
  Function nc_set_log_level(new_level) BIND(C)
   USE ISO_C_BINDING, ONLY: C_INT

   Integer(C_INT), VALUE :: new_level
   Integer(C_INT)        :: nc_set_log_level
 End Function nc_set_log_level
End Interface

 cnew_level = new_level
 cstatus = nc_set_log_level(cnew_level)

 status = cstatus

End Function nf_set_log_level

#endif
