module shr_carma_mod

  !================================================================================
  ! This reads the carma_inparm namelist in drv_flds_in and makes the relavent
  ! information available to CAM, CLM, and driver. 
  !================================================================================

  use shr_kind_mod , only : r8 => shr_kind_r8, CX => SHR_KIND_CX
  use shr_sys_mod  , only : shr_sys_abort
  use shr_log_mod  , only : logunit => shr_log_Unit
  use shr_nl_mod   , only : shr_nl_find_group_name

  implicit none
  private

  public :: shr_carma_readnl ! reads carma_inparm namelist

!-------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------

  subroutine shr_carma_readnl( NLFileName, carma_fields)

    !-------------------------------------------------------------------------
    ! This reads the carma_emis_nl namelist group in drv_flds_in and parses the
    ! namelist information for the driver, CLM, and CAM.
    !-------------------------------------------------------------------------

    use ESMF, only : ESMF_VM, ESMF_VMGetCurrent, ESMF_VMGet, ESMF_VMBroadcast

    character(len=*) , intent(in)  :: NLFileName
    character(len=CX), intent(out) :: carma_fields

    type(ESMF_VM) :: vm
    integer :: localPet
    integer :: rc
    integer :: unitn            ! namelist unit number
    integer :: ierr             ! error code
    logical :: exists           ! if file exists or not
    integer :: i, tmp(1)
    character(*),parameter :: F00   = "('(shr_carma_readnl) ',2a)"

    namelist /carma_inparm/ carma_fields

    carma_fields = ' '
    call ESMF_VMGetCurrent(vm, rc=rc)
    call ESMF_VMGet(vm, localpet=localpet, rc=rc)
    tmp = 0
    if (localpet==0) then
       inquire( file=trim(NLFileName), exist=exists)
       if ( exists ) then
          open(newunit=unitn, file=trim(NLFilename), status='old' )
          write(logunit,F00) 'Read in carma_inparm namelist from: ', trim(NLFilename)
          call shr_nl_find_group_name(unitn, 'carma_inparm', status=ierr)
          if (ierr == 0) then
             read(unitn, carma_inparm, iostat=ierr)
             if (ierr > 0) then
                call shr_sys_abort( 'problem on read of carma_inparm namelist in shr_carma_readnl' )
             endif
          else
             write(logunit,*) 'shr_carma_readnl:  no carma_inparm namelist found in ',NLFilename
          end if
          close( unitn )
       else
          write(logunit,*) 'shr_carma_readnl:  no file ',NLFilename, ' found'
       end if
       if (len_trim(carma_fields) > 0) tmp(1)=1
    end if
    call ESMF_VMBroadcast(vm, tmp, 1, 0, rc=rc)
    if(tmp(1) == 1) then
       call ESMF_VMBroadcast(vm, carma_fields, CX, 0, rc=rc)
    endif

  end subroutine shr_carma_readnl

endmodule shr_carma_mod
