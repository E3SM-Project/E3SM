module prep_iac_mod

#include "shr_assert.h"
  use shr_kind_mod,     only: r8 => SHR_KIND_R8
  use shr_kind_mod,     only: cs => SHR_KIND_CS
  use shr_kind_mod,     only: cl => SHR_KIND_CL
  use shr_kind_mod,     only: cxx => SHR_KIND_CXX
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct,     only: num_inst_lnd, num_inst_iac, num_inst_frc
  use seq_comm_mct,     only: CPLID, ROFID, logunit
  use seq_comm_mct,     only: seq_comm_getData=>seq_comm_setptrs
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata
  use shr_log_mod     , only: errMsg => shr_log_errMsg
  use seq_map_type_mod
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: iac, lnd
  use prep_lnd_mod, only: prep_lnd_get_mapper_Fr2l

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_iac_init
  public :: prep_iac_mrg

  public :: prep_iac_accum
  public :: prep_iac_accum_avg

  public :: prep_iac_calc_l2x_zx

  public :: prep_iac_get_l2zacc_lx
  public :: prep_iac_get_l2zacc_lx_cnt
  public :: prep_iac_get_mapper_Fl2z

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_Fl2z

  ! attribute vectors
  type(mct_aVect), pointer :: l2x_zx(:)

  ! accumulation variables
  type(mct_aVect), pointer :: l2zacc_lx(:)   ! lnd export, lnd grid, cpl pes
  integer        , target  :: l2zacc_lx_cnt  ! l2racc_lx: number of time samples accumulated

  ! other module variables
  integer :: mpicom_CPLID                            ! MPI cpl communicator

  !================================================================================================

contains

  !================================================================================================

  subroutine prep_iac_init(infodata, lnd_c2_iac)

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and all other non-mapping
    ! module variables
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    logical                 , intent(in)    :: lnd_c2_iac ! .true.  => lnd to iac coupling on
    !
    ! Local Variables

  end subroutine prep_iac_init

  !================================================================================================

  subroutine prep_iac_accum(timer)

    !---------------------------------------------------------------
    ! Description
    ! Accumulate land input to iac
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables

  end subroutine prep_iac_accum

  !================================================================================================

  subroutine prep_iac_accum_avg(timer)

    !---------------------------------------------------------------
    ! Description
    ! Finalize accumulation of land input to river component
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables

  end subroutine prep_iac_accum_avg

  !================================================================================================

  subroutine prep_iac_mrg(infodata, fractions_zx, timer_mrg)

    !---------------------------------------------------------------
    ! Description
    ! Merge iac inputs
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    type(mct_aVect)         , intent(in)    :: fractions_zx(:)
    character(len=*)        , intent(in)    :: timer_mrg
    !
    ! Local Variables

  end subroutine prep_iac_mrg

  !================================================================================================

  !================================================================================================

  subroutine prep_iac_calc_l2x_zx(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create l2x_zx (note that l2x_zx is a local module variable)
    !
    ! Arguments
    ! Don't know if we need these fractions just yet
    ! type(mct_aVect) , intent(in) :: fractions_lx(:)
    character(len=*), intent(in) :: timer
    !
    ! Local Variables

  end subroutine prep_iac_calc_l2x_zx

  !================================================================================================

  function prep_iac_get_l2zacc_lx()
    type(mct_aVect), pointer :: prep_iac_get_l2zacc_lx(:)
    prep_iac_get_l2zacc_lx => l2zacc_lx(:)
  end function prep_iac_get_l2zacc_lx

  function prep_iac_get_l2zacc_lx_cnt()
    integer, pointer :: prep_iac_get_l2zacc_lx_cnt
    prep_iac_get_l2zacc_lx_cnt => l2zacc_lx_cnt
  end function prep_iac_get_l2zacc_lx_cnt

  function prep_iac_get_mapper_Fl2z()
    type(seq_map), pointer :: prep_iac_get_mapper_Fl2z
    prep_iac_get_mapper_Fl2z => mapper_Fl2z
  end function prep_iac_get_mapper_Fl2z

end module prep_iac_mod
