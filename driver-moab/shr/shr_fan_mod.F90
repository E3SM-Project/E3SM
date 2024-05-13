!================================================================================
module shr_fan_mod

  use shr_kind_mod,only : r8 => shr_kind_r8
  use shr_kind_mod,only : CL => SHR_KIND_CL, CX => SHR_KIND_CX, CS => SHR_KIND_CS
  use shr_sys_mod, only : shr_sys_abort
  use shr_log_mod, only : loglev  => shr_log_Level
  use shr_log_mod, only : logunit => shr_log_Unit
  use shr_file_mod,   only : shr_file_getUnit, shr_file_freeUnit
  use seq_comm_mct,   only : seq_comm_iamroot, seq_comm_setptrs
  use shr_mpi_mod,    only : shr_mpi_bcast

  implicit none
  private

  public shr_fan_readnl

  logical, save, public :: shr_fan_to_atm = .false.
  character(len=CS), save, public :: shr_fan_fields_token = ''

contains

  subroutine shr_fan_readnl(nlfilename, id, fan_fields, have_fields)
    use shr_mpi_mod, only : shr_mpi_bcast
    use shr_nl_mod   , only : shr_nl_find_group_name

    character(len=*), intent(in)  :: nlfilename
    integer, intent(in) :: id ! seq_comm ID
    character(len=*), intent(out) :: fan_fields
    logical, intent(out) :: have_fields

    integer :: mpicomm, iostat, fileunit
    logical :: exists, fan_nh3_to_atm
    character(*),parameter :: subname = '(shr_fan_reanl) '

    namelist /fan_inparm/ fan_nh3_to_atm

    call seq_comm_setptrs(id, mpicom=mpicomm)
    ! Need to initilize to default value in case drv_flds_in doesn't exist
    fan_nh3_to_atm = .false. 
    if (seq_comm_iamroot(id)) then
       inquire(file=trim(nlfilename), exist=exists)
       if (exists) then
          fileunit = shr_file_getUnit()
          open(fileunit, file=trim(nlfilename), status='old' )
          call shr_nl_find_group_name(fileunit, 'fan_inparm', iostat)
          if (iostat /= 0) then
             write(logunit, *) subname, 'FAN/CAM coupling not specified'
             fan_nh3_to_atm = .false.
             !call shr_sys_abort(subName//'Error reading namelist')
          else          
             read(fileunit, fan_inparm, iostat=iostat)
             if (iostat /= 0) then
                call shr_sys_abort(subName//'Error reading namelist')
             end if
          end if
          close(fileunit)
          call shr_file_freeunit(fileunit)
       end if
    end if ! root
    call shr_mpi_bcast(fan_nh3_to_atm, mpicomm)
    have_fields = fan_nh3_to_atm
    if (fan_nh3_to_atm) then
       fan_fields = 'Fall_FANNH3'
    else
       fan_fields = ''
    end if
    shr_fan_to_atm = have_fields
    shr_fan_fields_token = fan_fields

  end subroutine shr_fan_readnl

endmodule shr_fan_mod
