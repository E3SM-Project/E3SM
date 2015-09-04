MODULE shr_log_mod

   use shr_kind_mod

   !----------------------------------------------------------------------------
   ! low-level shared variables for logging, these may not be parameters
   !----------------------------------------------------------------------------
   public

   integer(SHR_KIND_IN) :: shr_log_Level = 1
   integer(SHR_KIND_IN) :: shr_log_Unit  = 6

END MODULE shr_log_mod
