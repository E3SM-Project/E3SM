module ocn_comp_mct

   ! seq mods
   use seq_cdata_mod, only: seq_cdata

   ! toolkits mods
   use esmf, only: ESMF_Clock
   use mct_mod, only: mct_aVect

   implicit none
   save
   private

   !---------------------------------------------------------------------------
   ! Public interfaces
   !---------------------------------------------------------------------------
   public :: ocn_init_mct
   public :: ocn_run_mct
   public :: ocn_final_mct

contains
   subroutine ocn_init_mct(EClock, cdata, x2o, o2x, NLFilename)
      ! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Clock), intent(inout) :: EClock
      type(seq_cdata), intent(inout) :: cdata
      type(mct_aVect), intent(inout) :: x2o, o2x
      character(len=*), optional, intent(in) :: NLFilename
   end subroutine ocn_init_mct

   subroutine ocn_run_mct(EClock, cdata, x2o, o2x)
      ! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Clock), intent(inout) :: EClock
      type(seq_cdata), intent(inout) :: cdata
      type(mct_aVect), intent(inout) :: x2o, o2x
   end subroutine ocn_run_mct

   subroutine ocn_final_mct(EClock, cdata, x2o, o2x)
      ! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_Clock), intent(inout) :: EClock
      type(seq_cdata), intent(inout) :: cdata
      type(mct_aVect), intent(inout) :: x2o, o2x
   end subroutine ocn_final_mct
end module ocn_comp_mct
