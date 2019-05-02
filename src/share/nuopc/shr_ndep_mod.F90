module shr_ndep_mod

  !========================================================================
  ! Module for handling nitrogen depostion of tracers.
  ! This module is shared by land and atmosphere models for the computations of
  ! dry deposition of tracers
  !========================================================================

  use ESMF         , only : ESMF_VMGetCurrent, ESMF_VM, ESMF_VMBroadcast, ESMF_VMGet
  use ESMF         , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
  use shr_sys_mod,   only : shr_sys_abort
  use shr_log_mod  , only : s_logunit => shr_log_Unit
  use shr_kind_mod,  only : r8 => shr_kind_r8
  use shr_file_mod , only : shr_file_getUnit, shr_file_freeUnit
  use shr_nl_mod   , only : shr_nl_find_group_name

  implicit none
  private

  ! !PUBLIC MEMBER FUNCTIONS
  public :: shr_ndep_readnl       ! Read namelist

  character(len=*), parameter :: &
       u_FILE_u=__FILE__

!====================================================================================
CONTAINS
!====================================================================================

  subroutine shr_ndep_readnl(NLFilename, ndep_nflds)

    !========================================================================
    ! reads ndep_inparm namelist and sets up driver list of fields for
    ! atmosphere -> land and atmosphere -> ocn communications.
    !========================================================================

    ! input/output variables
    character(len=*), intent(in)  :: NLFilename ! Namelist filename
    integer         , intent(out) :: ndep_nflds

    !----- local -----
    type(ESMF_VM)      :: vm
    integer            :: i                      ! Indices
    integer            :: unitn                  ! namelist unit number
    integer            :: ierr                   ! error code
    integer            :: tmp(1)
    logical            :: exists                 ! if file exists or not
    integer            :: rc
    integer, parameter :: maxspc = 100           ! Maximum number of species
    character(len=32)  :: ndep_list(maxspc) = '' ! List of ndep species
    integer            :: localpet
    character(*),parameter :: subName = '(shr_ndep_read) '
    character(*),parameter :: F00   = "('(shr_ndep_read) ',8a)"
    character(*),parameter :: FI1   = "('(shr_ndep_init) ',a,I2)"
    ! ------------------------------------------------------------------

    namelist /ndep_inparm/ ndep_list

    !-----------------------------------------------------------------------------
    ! Read namelist and figure out the ndep field list to pass
    ! First check if file exists and if not, n_ndep will be zero
    !-----------------------------------------------------------------------------

    !--- Open and read namelist ---
    if ( len_trim(NLFilename) == 0 ) then
       call shr_sys_abort( subName//'ERROR: nlfilename not set' )
    end if

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return 

    call ESMF_VMGet(vm, localpet=localpet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return 

    ndep_nflds=0
    if (localpet==0) then
       inquire( file=trim(NLFileName), exist=exists)
       if ( exists ) then
          unitn = shr_file_getUnit()
          open( unitn, file=trim(NLFilename), status='old' )
          write(s_logunit,F00) 'Read in ndep_inparm namelist from: ', trim(NLFilename)
          call shr_nl_find_group_name(unitn, 'ndep_inparm', ierr)
          if (ierr == 0) then
             ierr = 1
             do while ( ierr /= 0 )
                read(unitn, ndep_inparm, iostat=ierr)
                if (ierr < 0) then
                   call shr_sys_abort( subName//'ERROR: encountered end-of-file on namelist read' )
                endif
             end do
          else
             write(s_logunit,*) 'shr_ndep_readnl:  no ndep_inparm namelist found in ',NLFilename
          endif
          close( unitn )
          call shr_file_freeUnit( unitn )
          do i=1,maxspc
             if (len_trim(ndep_list(i)) > 0) then
                ndep_nflds = ndep_nflds+1
             endif
          enddo
       end if
    end if
    tmp = ndep_nflds
    call ESMF_VMBroadcast(vm, tmp, 1, 0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return 
    ndep_nflds=tmp(1)

  end subroutine shr_ndep_readnl

end module shr_ndep_mod
