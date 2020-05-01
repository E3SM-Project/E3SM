module shr_ndep_mod

  !========================================================================
  ! Module for handling nitrogen depostion of tracers.
  ! This module is shared by land and atmosphere models for the computations of
  ! dry deposition of tracers
  !========================================================================

  use ESMF         , only : ESMF_VMGetCurrent, ESMF_VM, ESMF_VMGet
  use ESMF         , only : ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU, ESMF_SUCCESS
  use shr_sys_mod  , only : shr_sys_abort
  use shr_log_mod  , only : s_logunit => shr_log_Unit
  use shr_kind_mod , only : r8 => shr_kind_r8
  use shr_nl_mod   , only : shr_nl_find_group_name
  use shr_mpi_mod  , only : shr_mpi_bcast

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
    logical            :: exists                 ! if file exists or not
    integer            :: rc
    integer, parameter :: maxspc = 100           ! Maximum number of species
    character(len=32)  :: ndep_list(maxspc) = '' ! List of ndep species
    integer            :: localpet
    integer            :: mpicom
    character(*),parameter :: F00   = "('(shr_ndep_read) ',8a)"
    character(*),parameter :: FI1   = "('(shr_ndep_init) ',a,I2)"
    character(*),parameter :: subName = '(shr_ndep_read) '
    ! ------------------------------------------------------------------

    namelist /ndep_inparm/ ndep_list

    !-----------------------------------------------------------------------------
    ! Read namelist and figure out the ndep field list to pass
    ! First check if file exists and if not, n_ndep will be zero
    !-----------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    !--- Open and read namelist ---
    if ( len_trim(NLFilename) == 0 ) then
       call shr_sys_abort( subName//'ERROR: nlfilename not set' )
    end if

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return 

    call ESMF_VMGet(vm, localPet=localPet, mpiCommunicator=mpicom, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return 

    ! Note the following still needs to be called on all processors since the mpi_bcast is a collective 
    ! call on all the pes of mpicom
    if (localpet==0) then
       inquire( file=trim(NLFileName), exist=exists)
       if ( exists ) then
          open(newunit=unitn, file=trim(NLFilename), status='old' )
          write(s_logunit,F00) 'Read in ndep_inparm namelist from: ', trim(NLFilename)
          call shr_nl_find_group_name(unitn, 'ndep_inparm', ierr)
          if (ierr == 0) then
             ! Note that ierr /= 0, no namelist is present.
             read(unitn, ndep_inparm, iostat=ierr)
             if (ierr > 0) then
                call shr_sys_abort(trim(subName) //'problem of read of ndep_inparm ')
             endif
          endif
          close( unitn )
       end if
    end if
    call shr_mpi_bcast( ndep_list, mpicom )

    ndep_nflds = 0
    do i=1,maxspc
       if (len_trim(ndep_list(i)) > 0) then
          ndep_nflds = ndep_nflds+1
       endif
    enddo

  end subroutine shr_ndep_readnl

end module shr_ndep_mod
