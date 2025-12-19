module RDycoreRestFile

  use shr_kind_mod   , only : r8 => shr_kind_r8, i8=>shr_kind_i8, shr_kind_cl
  use shr_sys_mod    , only : shr_sys_flush, shr_sys_abort
  use rdycore_varctl , only : iulog, caseid,inst_suffix
  use rdycoreSpmdMod , only : masterproc

  implicit none
  private

  public :: RDycoreRestFileNameBase
  public :: RDycoreRestFileName

contains

  !------------------------------------------------------------------------
  character(len=1024) function RDycoreRestFileNameBase( rdate )
    use rdycore_varctl, only : inst_suffix

    implicit none
    character(len=*), intent(in) :: rdate   ! input date for restart file name

    RDycoreRestFileNameBase = "./"//trim(caseid)//".rdycore"//trim(inst_suffix)//".r."//trim(rdate)!//".h5"

  end function RDycoreRestFileNameBase

  !------------------------------------------------------------------------
  character(len=1024) function RDycoreRestFileName( rdate )
    use rdycore_varctl, only : inst_suffix

    implicit none
    character(len=*), intent(in) :: rdate   ! input date for restart file name

    RDycoreRestFileName = "./"//trim(caseid)//".rdycore"//trim(inst_suffix)//".r."//trim(rdate)//".h5"

  end function RDycoreRestFileName

end module RDycoreRestFile
