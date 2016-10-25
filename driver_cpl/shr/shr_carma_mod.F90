!================================================================================
! This reads the carma_inparm namelist in drv_flds_in and makes the relavent 
! information available to CAM, CLM, and driver. The driver sets up CLM to CAM 
! communication for the  VOC flux fields. CLM needs to know what specific VOC
! fluxes need to be passed to the coupler and how to assimble the fluxes.  
! CAM needs to know what specific VOC fluxes to expect from CLM.
!
! Mariana Vertenstein -- 24 Sep 2012
!================================================================================
module shr_carma_mod

  use shr_kind_mod,only : r8 => shr_kind_r8
  use shr_kind_mod,only : CL => SHR_KIND_CL, CX => SHR_KIND_CX, CS => SHR_KIND_CS
  use shr_sys_mod, only : shr_sys_abort
  use shr_log_mod, only : loglev  => shr_log_Level
  use shr_log_mod, only : logunit => shr_log_Unit

  implicit none
  save
  private

  public :: shr_carma_readnl           ! reads carma_inparm namelist

contains

  !-------------------------------------------------------------------------
  ! This reads the carma_emis_nl namelist group in drv_flds_in and parses the
  ! namelist information for the driver, CLM, and CAM.
  !-------------------------------------------------------------------------
  subroutine shr_carma_readnl( NLFileName, carma_fields )

    use shr_nl_mod,     only : shr_nl_find_group_name
    use shr_file_mod,   only : shr_file_getUnit, shr_file_freeUnit

    character(len=*),  intent(in)  :: NLFileName
    character(len=CX), intent(out) :: carma_fields	

    integer :: unitn            ! namelist unit number
    integer :: ierr             ! error code
    logical :: exists           ! if file exists or not
    character(*),parameter :: F00   = "('(shr_carma_readnl) ',2a)" 

    namelist /carma_inparm/ carma_fields

    carma_fields = ' '

    inquire( file=trim(NLFileName), exist=exists)
    if ( exists ) then
       unitn = shr_file_getUnit()
       open( unitn, file=trim(NLFilename), status='old' )
       if ( loglev > 0 ) write(logunit,F00) &
            'Read in carma_inparm namelist from: ', trim(NLFilename)

       call shr_nl_find_group_name(unitn, 'carma_inparm', status=ierr)
       ! If ierr /= 0, no namelist present.

       if (ierr == 0) then
          read(unitn, carma_inparm, iostat=ierr)
          if (ierr > 0) then
             call shr_sys_abort( 'problem on read of carma_inparm namelist in shr_carma_readnl' )
          endif
       end if

       close( unitn )
       call shr_file_freeUnit( unitn )
    end if

  end subroutine shr_carma_readnl

endmodule shr_carma_mod
