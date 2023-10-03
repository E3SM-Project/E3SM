!=====================================================================================
! Module for handling dust-related information shared by different components models.
!
! History:
! - Hui Wan, Jan 2023, following examples in driver-mct/shr/shr_*_mod.F90
!=====================================================================================
module shr_dust_mod

  use shr_sys_mod,   only : shr_sys_abort
  use shr_log_mod,   only : s_loglev  => shr_log_Level

  implicit none
  save
  private

  public :: shr_dust_readnl       ! subroutine that reads namelist shr_dust_nl
  public :: dust_emis_scheme      ! module parameter

  integer :: dust_emis_scheme = 1

CONTAINS

  !-------------------------------------------------------------------------
  ! Reads namelist shr_dust_nl
  !-------------------------------------------------------------------------
  subroutine shr_dust_readnl(NLFilename, ID)

    use shr_file_mod , only : shr_file_getUnit, shr_file_freeUnit
    use shr_log_mod  , only : s_logunit => shr_log_Unit
    use seq_comm_mct , only : seq_comm_iamroot, seq_comm_setptrs
    use shr_mpi_mod  , only : shr_mpi_bcast
    use shr_nl_mod   , only : shr_nl_find_group_name

    implicit none

    character(len=*), intent(in)  :: NLFilename ! Namelist filename
    integer         , intent(in)  :: ID         ! seq_comm ID

    !----- local -----
    integer :: unitn            ! namelist unit number
    integer :: ierr             ! error code
    logical :: exists           ! if file exists or not
    integer :: mpicom           ! MPI communicator

    !----- formats -----
    character(*),parameter :: subName = '(shr_dust_readnl) '
    character(*),parameter :: F00   = "('(shr_dust_readnl) ',8a)"

    namelist /shr_dust_nl/ dust_emis_scheme

    !--- Open and read namelist ---
    if ( len_trim(NLFilename) == 0 ) then
       call shr_sys_abort( subName//'ERROR: nlfilename not set' )
    end if
    call seq_comm_setptrs(ID,mpicom=mpicom)
    if (seq_comm_iamroot(ID)) then
       inquire( file=trim(NLFileName), exist=exists)
       if ( exists ) then
          unitn = shr_file_getUnit()
          open( unitn, file=trim(NLFilename), status='old' )
          if ( s_loglev > 0 ) then
             write(s_logunit,F00) 'Read in shr_dust_nl namelist from: ', trim(NLFilename)
          end if
          call shr_nl_find_group_name(unitn, 'shr_dust_nl', ierr)
          if (ierr == 0) then
             ierr = 1
             do while ( ierr /= 0 )
                read(unitn, shr_dust_nl, iostat=ierr)
                if (ierr < 0) then
                   call shr_sys_abort( subName//'ERROR: encountered end-of-file on namelist read' )
                endif
             end do
          else
             write(s_logunit,*) 'shr_dust_readnl:  no shr_dust_nl namelist found in ',NLFilename
          endif
          close( unitn )
          call shr_file_freeUnit( unitn )
       end if
    end if
    call shr_mpi_bcast( dust_emis_scheme, mpicom )

  end subroutine shr_dust_readnl

end module shr_dust_mod
