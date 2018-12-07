module seq_comm_mct

!---------------------------------------------------------------------
!
! Purpose: Set up necessary communications
!          Note that if no MPI, will call MCTs fake version
!          (including mpif.h) will be utilized
!
!---------------------------------------------------------------------


!!! NOTE: If all atmospheres are identical in number of processes,
!!! number of threads, and grid layout, we should check that the
!!! user-provided number of processes and threads are consistent
!!! (or else, only accept one entry for these quantities when reading
!!! the namelist).  ARE OTHER PROTECTIONS/CHECKS NEEDED???
  use ESMF, only : ESMF_LogKind_Flag
  implicit none

  private

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public seq_comm_setcomm
  public seq_comm_iamin
  public seq_comm_iamroot
  public seq_comm_mpicom
  public seq_comm_iam
  public seq_comm_gloroot
  public seq_comm_name
  public seq_comm_inst
  public seq_comm_suffix
  public seq_comm_petlist
  public seq_comm_setptrs
  public seq_comm_setnthreads
  public seq_comm_getnthreads
  public seq_comm_printcomms

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

  integer, public :: logunit  = 6     ! log unit number
  integer, public :: loglevel = 1     ! log level

  ! NOTE: NUM_COMP_INST_XXX are cpp variables set in buildlib.csm_share
  integer, parameter, public :: num_inst_atm = NUM_COMP_INST_ATM
  integer, parameter, public :: num_inst_lnd = NUM_COMP_INST_LND
  integer, parameter, public :: num_inst_ocn = NUM_COMP_INST_OCN
  integer, parameter, public :: num_inst_ice = NUM_COMP_INST_ICE
  integer, parameter, public :: num_inst_glc = NUM_COMP_INST_GLC
  integer, parameter, public :: num_inst_wav = NUM_COMP_INST_WAV
  integer, parameter, public :: num_inst_rof = NUM_COMP_INST_ROF
  integer, parameter, public :: num_inst_esp = NUM_COMP_INST_ESP

  integer, public :: num_inst_min, num_inst_max

  integer, parameter, public :: num_inst_total = &
       num_inst_atm + num_inst_lnd + num_inst_ocn + num_inst_ice + &
       num_inst_glc + num_inst_wav + num_inst_rof + num_inst_esp + 1

  integer, parameter :: ncouplers  = 1  ! number of couplers
  integer, parameter :: ncomps = (ncouplers + num_inst_total)

  integer, public :: GLOID
  integer, public :: CPLID
  integer, public :: ATMID(num_inst_atm)
  integer, public :: LNDID(num_inst_lnd)
  integer, public :: OCNID(num_inst_ocn)
  integer, public :: ICEID(num_inst_ice)
  integer, public :: GLCID(num_inst_glc)
  integer, public :: ROFID(num_inst_rof)
  integer, public :: WAVID(num_inst_wav)
  integer, public :: ESPID(num_inst_esp)

  type(ESMF_LogKind_Flag), public :: esmf_logfile_kind

  integer, parameter, public :: seq_comm_namelen=16

  type seq_comm_type
    character(len=seq_comm_namelen) :: name     ! my name
    character(len=seq_comm_namelen) :: suffix   ! recommended suffix
    integer :: inst            ! my inst index
    integer :: ID              ! my id number
    integer :: mpicom          ! mpicom
    integer :: mpigrp          ! mpigrp
    integer :: npes            ! number of mpi tasks in comm
    integer :: nthreads        ! number of omp threads per task
    integer :: iam             ! my task number in mpicom
    logical :: iamroot         ! am i the root task in mpicom
    integer :: gloroot         ! the global task number of each comps root on all pes
    integer :: pethreads       ! max number of threads on my task
    logical :: set             ! has this datatype been set
    integer, pointer    :: petlist(:)  ! esmf pet list
    logical :: petlist_allocated ! whether the petlist pointer variable was allocated
  end type seq_comm_type

  type(seq_comm_type) :: seq_comms(ncomps)

  character(*), parameter :: layout_concurrent = 'concurrent'
  character(*), parameter :: layout_sequential = 'sequential'

  character(*), parameter :: F11 = "(a,a,'(',i3,' ',a,')',a,   3i6,' (',a,i6,')',' (',a,i3,')')"
  character(*), parameter :: F12 = "(a,a,'(',i3,' ',a,')',a,2i6,6x,' (',a,i6,')',' (',a,i3,')','(',a,2i6,')')"
  character(*), parameter :: F13 = "(a,a,'(',i3,' ',a,')',a,2i6,6x,' (',a,i6,')',' (',a,i3,')')"
  character(*), parameter :: F14 = "(a,a,'(',i3,' ',a,')',a,    6x,' (',a,i6,')',' (',a,i3,')')"
  integer :: Global_Comm


  character(len=32), public :: &
       atm_layout, lnd_layout, ice_layout, glc_layout, rof_layout, &
       ocn_layout, wav_layout, esp_layout

  logical :: seq_comm_mct_initialized = .false.  ! whether this module has been initialized

!=======================================================================
contains
!======================================================================

  subroutine seq_comm_init(Comm_in, nmlfile)
    use shr_sys_mod, only : shr_sys_abort
    use shr_file_mod, only : shr_file_getUnit, shr_file_freeUnit
    use shr_mpi_mod, only : shr_mpi_max, shr_mpi_chkerr, shr_mpi_bcast
    use mpi, only : mpi_comm_rank, mpi_comm_size, mpi_comm_null, MPI_INTEGER, MPI_CHARACTER
    use mpi, only : mpi_comm_null, mpi_group_null
    use ESMF, only : ESMF_LOGKIND_NONE, ESMF_LOGKIND_MULTI, ESMF_LOGKIND_SINGLE
    use mct_mod, only: mct_die, mct_world_init
    !----------------------------------------------------------
    ! Arguments
    integer, intent(in) :: Comm_in
    character(len=*), intent(IN) :: nmlfile
    !
    ! Local variables
    !
    logical :: error_state
    integer :: ierr, n, count
    integer :: mpi_group_world   ! MPI_COMM_WORLD group
    integer :: mype,numpes,myncomps,max_threads,gloroot
    integer :: atm_inst_tasks, lnd_inst_tasks, ocn_inst_tasks, ice_inst_tasks
    integer :: glc_inst_tasks, rof_inst_tasks, wav_inst_tasks, esp_inst_tasks
    integer :: current_task_rootpe, droot
    integer :: amin(num_inst_atm), amax(num_inst_atm), astr(num_inst_atm)
    integer :: lmin(num_inst_lnd), lmax(num_inst_lnd), lstr(num_inst_lnd)
    integer :: imin(num_inst_ice), imax(num_inst_ice), istr(num_inst_ice)
    integer :: omin(num_inst_ocn), omax(num_inst_ocn), ostr(num_inst_ocn)
    integer :: gmin(num_inst_glc), gmax(num_inst_glc), gstr(num_inst_glc)
    integer :: wmin(num_inst_wav), wmax(num_inst_wav), wstr(num_inst_wav)
    integer :: rmin(num_inst_rof), rmax(num_inst_rof), rstr(num_inst_rof)
    integer :: emin(num_inst_esp), emax(num_inst_esp), estr(num_inst_esp)
    integer :: cmin,cmax,cstr
    integer :: pelist(3,1)       ! start, stop, stride for group
    integer :: nu, i
    character(len=24) :: esmf_logging
    integer, pointer :: comps(:) ! array with component ids
    integer, pointer :: comms(:) ! array with mpicoms
    character(*), parameter :: subName =   '(seq_comm_init) '

    integer :: &
         atm_ntasks, atm_rootpe, atm_pestride, atm_nthreads, &
         lnd_ntasks, lnd_rootpe, lnd_pestride, lnd_nthreads, &
         ice_ntasks, ice_rootpe, ice_pestride, ice_nthreads, &
         glc_ntasks, glc_rootpe, glc_pestride, glc_nthreads, &
         wav_ntasks, wav_rootpe, wav_pestride, wav_nthreads, &
         rof_ntasks, rof_rootpe, rof_pestride, rof_nthreads, &
         ocn_ntasks, ocn_rootpe, ocn_pestride, ocn_nthreads, &
         esp_ntasks, esp_rootpe, esp_pestride, esp_nthreads, &
         cpl_ntasks, cpl_rootpe, cpl_pestride, cpl_nthreads

    namelist /cime_pes/  &
         atm_ntasks, atm_rootpe, atm_pestride, atm_nthreads, atm_layout, &
         lnd_ntasks, lnd_rootpe, lnd_pestride, lnd_nthreads, lnd_layout, &
         ice_ntasks, ice_rootpe, ice_pestride, ice_nthreads, ice_layout, &
         glc_ntasks, glc_rootpe, glc_pestride, glc_nthreads, glc_layout, &
         wav_ntasks, wav_rootpe, wav_pestride, wav_nthreads, wav_layout, &
         rof_ntasks, rof_rootpe, rof_pestride, rof_nthreads, rof_layout, &
         ocn_ntasks, ocn_rootpe, ocn_pestride, ocn_nthreads, ocn_layout, &
         esp_ntasks, esp_rootpe, esp_pestride, esp_nthreads, esp_layout, &
         cpl_ntasks, cpl_rootpe, cpl_pestride, cpl_nthreads, esmf_logging
    !----------------------------------------------------------

    ! make sure this is first pass and set comms unset
    if (seq_comm_mct_initialized) then
       write(logunit,*) trim(subname),' ERROR seq_comm_init already called '
       call shr_sys_abort()
    endif
    seq_comm_mct_initialized = .true.
    Global_Comm = Comm_in

    !! Initialize seq_comms elements

    do n = 1,ncomps
       seq_comms(n)%name = 'unknown'
       seq_comms(n)%suffix = ' '
       seq_comms(n)%inst = 0
       seq_comms(n)%set = .false.
       seq_comms(n)%petlist_allocated = .false.
       seq_comms(n)%mpicom = MPI_COMM_NULL    ! do some initialization here
       seq_comms(n)%iam = -1
       seq_comms(n)%iamroot = .false.
       seq_comms(n)%npes = -1
       seq_comms(n)%nthreads = -1
       seq_comms(n)%gloroot = -1
       seq_comms(n)%pethreads = -1
    enddo

    ! Initialize MPI -  Note that if no MPI, will call MCTs fake version
    call mpi_comm_rank(GLOBAL_COMM, mype  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank comm_world')

    call mpi_comm_size(GLOBAL_COMM, numpes, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size comm_world')

    ! Set ntasks, rootpe, pestride, nthreads for all components
    if (mype == 0) then

       !! Set up default component process parameters
       atm_ntasks = numpes
       atm_rootpe = 0
       atm_pestride = 1
       atm_nthreads = 1
       atm_layout = trim(layout_concurrent)

       lnd_ntasks = numpes
       lnd_rootpe = 0
       lnd_pestride = 1
       lnd_nthreads = 1
       lnd_layout = trim(layout_concurrent)

       ocn_ntasks = numpes
       ocn_rootpe = 0
       ocn_pestride = 1
       ocn_nthreads = 1
       ocn_layout = trim(layout_concurrent)

       ice_ntasks = numpes
       ice_rootpe = 0
       ice_pestride = 1
       ice_nthreads = 1
       ice_layout = trim(layout_concurrent)

       glc_ntasks = numpes
       glc_rootpe = 0
       glc_pestride = 1
       glc_nthreads = 1
       glc_layout = trim(layout_concurrent)

       rof_ntasks = numpes
       rof_rootpe = 0
       rof_pestride = 1
       rof_nthreads = 1
       rof_layout = trim(layout_concurrent)

       wav_ntasks = numpes
       wav_rootpe = 0
       wav_pestride = 1
       wav_nthreads = 1
       wav_layout = trim(layout_concurrent)

       esp_ntasks = numpes
       esp_rootpe = 0
       esp_pestride = 1
       esp_nthreads = 1
       esp_layout = trim(layout_concurrent)

       cpl_ntasks = numpes
       cpl_rootpe = 0
       cpl_pestride = 1
       cpl_nthreads = 1

       esmf_logging = "ESMF_LOGKIND_NONE"

       ! Read namelist if it exists
       ! TODO: obtain this from attributes
       nu = shr_file_getUnit()
       open(nu, file=trim(nmlfile), status='old', iostat=ierr)

       if (ierr == 0) then
          ierr = 1
          do while( ierr > 0 )
             read(nu, nml=cime_pes, iostat=ierr)
          end do
          close(nu)
       end if
       call shr_file_freeUnit(nu)
    end if

    !--- compute num_inst_min, num_inst_max
    !--- instances must be either 1 or a constant across components
    !--- checks for prognostic/present consistency in the driver

    error_state = .false.
    num_inst_min = num_inst_atm
    num_inst_min = min(num_inst_min, num_inst_lnd)
    num_inst_min = min(num_inst_min, num_inst_ocn)
    num_inst_min = min(num_inst_min, num_inst_ice)
    num_inst_min = min(num_inst_min, num_inst_glc)
    num_inst_min = min(num_inst_min, num_inst_wav)
    num_inst_min = min(num_inst_min, num_inst_rof)
    ! ESP is currently limited to one instance, should not affect other comps
    !    num_inst_min = min(num_inst_min, num_inst_esp)
    num_inst_max = num_inst_atm
    num_inst_max = max(num_inst_max, num_inst_lnd)
    num_inst_max = max(num_inst_max, num_inst_ocn)
    num_inst_max = max(num_inst_max, num_inst_ice)
    num_inst_max = max(num_inst_max, num_inst_glc)
    num_inst_max = max(num_inst_max, num_inst_wav)
    num_inst_max = max(num_inst_max, num_inst_rof)
    num_inst_max = max(num_inst_max, num_inst_esp)

    if (num_inst_min /= num_inst_max .and. num_inst_min /= 1) error_state = .true.
    if (num_inst_atm /= num_inst_min .and. num_inst_atm /= num_inst_max) error_state = .true.
    if (num_inst_lnd /= num_inst_min .and. num_inst_lnd /= num_inst_max) error_state = .true.
    if (num_inst_ocn /= num_inst_min .and. num_inst_ocn /= num_inst_max) error_state = .true.
    if (num_inst_ice /= num_inst_min .and. num_inst_ice /= num_inst_max) error_state = .true.
    if (num_inst_glc /= num_inst_min .and. num_inst_glc /= num_inst_max) error_state = .true.
    if (num_inst_wav /= num_inst_min .and. num_inst_wav /= num_inst_max) error_state = .true.
    if (num_inst_rof /= num_inst_min .and. num_inst_rof /= num_inst_max) error_state = .true.
    if (num_inst_esp /= 1) then
       write(logunit,*) trim(subname),' ERROR: ESP restricted to one instance'
       error_state = .true.
    end if

    if (error_state) then
       write(logunit,*) trim(subname),' ERROR: num_inst inconsistent'
       call shr_sys_abort(trim(subname)//' ERROR: num_inst inconsistent')
    endif

    ! Initialize IDs

    count = 0
    count = count + 1
    GLOID = count
    count = count + 1
    CPLID = count
    do n = 1, num_inst_atm
       count = count + 1
       ATMID(n) = count
    end do
    do n = 1, num_inst_lnd
       count = count + 1
       LNDID(n) = count
    end do
    do n = 1, num_inst_ocn
       count = count + 1
       OCNID(n) = count
    end do
    do n = 1, num_inst_ice
       count = count + 1
       ICEID(n) = count
    end do
    do n = 1, num_inst_glc
       count = count + 1
       GLCID(n) = count
    end do
    do n = 1, num_inst_rof
       count = count + 1
       ROFID(n) = count
    end do
    do n = 1, num_inst_wav
       count = count + 1
       WAVID(n) = count
    end do
    do n = 1, num_inst_esp
       count = count + 1
       ESPID(n) = count
    end do
    if (count /= ncomps) then
       write(logunit,*) trim(subname),' ERROR in ID count ',count,ncomps
       call shr_sys_abort(trim(subname)//' ERROR in ID count')
    endif

    if (mype == 0) then
       !--- validation of inputs ---
       ! rootpes >= 0

       error_state = .false.

       if (atm_rootpe < 0) error_state = .true.
       if (lnd_rootpe < 0) error_state = .true.
       if (ice_rootpe < 0) error_state = .true.
       if (ocn_rootpe < 0) error_state = .true.
       if (glc_rootpe < 0) error_state = .true.
       if (wav_rootpe < 0) error_state = .true.
       if (rof_rootpe < 0) error_state = .true.
       if (esp_rootpe < 0) error_state = .true.
       if (cpl_rootpe < 0) error_state = .true.

       if (error_state) then
          write(logunit,*) trim(subname),' ERROR: rootpes must be >= 0'
          call shr_sys_abort(trim(subname)//' ERROR: rootpes >= 0')
       endif

       !! Determine the process layout
       !!
       !! We will assign atm_ntasks / num_inst_atm tasks to each atmosphere
       !! instance.  (This may lead to unallocated tasks if atm_ntasks is
       !! not an integer multiple of num_inst_atm.)

       if (trim(atm_layout) == trim(layout_concurrent)) then
          atm_inst_tasks = atm_ntasks / num_inst_atm
          droot = (atm_inst_tasks * atm_pestride)
       elseif (trim(atm_layout) == trim(layout_sequential)) then
          atm_inst_tasks = atm_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid atm_layout ')
       endif
       current_task_rootpe = atm_rootpe
       do n = 1, num_inst_atm
          amin(n) = current_task_rootpe
          amax(n) = current_task_rootpe &
                    + ((atm_inst_tasks - 1) * atm_pestride)
          astr(n) = atm_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! Land instance tasks

       if (trim(lnd_layout) == trim(layout_concurrent)) then
          lnd_inst_tasks = lnd_ntasks / num_inst_lnd
          droot = (lnd_inst_tasks * lnd_pestride)
       elseif (trim(lnd_layout) == trim(layout_sequential)) then
          lnd_inst_tasks = lnd_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid lnd_layout ')
       endif
       current_task_rootpe = lnd_rootpe
       do n = 1, num_inst_lnd
          lmin(n) = current_task_rootpe
          lmax(n) = current_task_rootpe &
                    + ((lnd_inst_tasks - 1) * lnd_pestride)
          lstr(n) = lnd_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! Ocean instance tasks

       if (trim(ocn_layout) == trim(layout_concurrent)) then
          ocn_inst_tasks = ocn_ntasks / num_inst_ocn
          droot = (ocn_inst_tasks * ocn_pestride)
       elseif (trim(ocn_layout) == trim(layout_sequential)) then
          ocn_inst_tasks = ocn_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid ocn_layout ')
       endif
       current_task_rootpe = ocn_rootpe
       do n = 1, num_inst_ocn
          omin(n) = current_task_rootpe
          omax(n) = current_task_rootpe &
                    + ((ocn_inst_tasks - 1) * ocn_pestride)
          ostr(n) = ocn_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! Sea ice instance tasks

       if (trim(ice_layout) == trim(layout_concurrent)) then
          ice_inst_tasks = ice_ntasks / num_inst_ice
          droot = (ice_inst_tasks * ice_pestride)
       elseif (trim(ice_layout) == trim(layout_sequential)) then
          ice_inst_tasks = ice_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid ice_layout ')
       endif
       current_task_rootpe = ice_rootpe
       do n = 1, num_inst_ice
          imin(n) = current_task_rootpe
          imax(n) = current_task_rootpe &
                    + ((ice_inst_tasks - 1) * ice_pestride)
          istr(n) = ice_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! Glacier instance tasks

       if (trim(glc_layout) == trim(layout_concurrent)) then
          glc_inst_tasks = glc_ntasks / num_inst_glc
          droot = (glc_inst_tasks * glc_pestride)
       elseif (trim(glc_layout) == trim(layout_sequential)) then
          glc_inst_tasks = glc_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid glc_layout ')
       endif
       current_task_rootpe = glc_rootpe
       do n = 1, num_inst_glc
          gmin(n) = current_task_rootpe
          gmax(n) = current_task_rootpe &
                    + ((glc_inst_tasks - 1) * glc_pestride)
          gstr(n) = glc_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! Runoff instance tasks

       if (trim(rof_layout) == trim(layout_concurrent)) then
          rof_inst_tasks = rof_ntasks / num_inst_rof
          droot = (rof_inst_tasks * rof_pestride)
       elseif (trim(rof_layout) == trim(layout_sequential)) then
          rof_inst_tasks = rof_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid rof_layout ')
       endif
       current_task_rootpe = rof_rootpe
       do n = 1, num_inst_rof
          rmin(n) = current_task_rootpe
          rmax(n) = current_task_rootpe &
                    + ((rof_inst_tasks - 1) * rof_pestride)
          rstr(n) = rof_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! Wave instance tasks

       if (trim(wav_layout) == trim(layout_concurrent)) then
          wav_inst_tasks = wav_ntasks / num_inst_wav
          droot = (wav_inst_tasks * wav_pestride)
       elseif (trim(wav_layout) == trim(layout_sequential)) then
          wav_inst_tasks = wav_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid wav_layout ')
       endif
       current_task_rootpe = wav_rootpe
       do n = 1, num_inst_wav
          wmin(n) = current_task_rootpe
          wmax(n) = current_task_rootpe &
                    + ((wav_inst_tasks - 1) * wav_pestride)
          wstr(n) = wav_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! External System Processing instance tasks

       if (trim(esp_layout) == trim(layout_concurrent)) then
          esp_inst_tasks = esp_ntasks / num_inst_esp
          droot = (esp_inst_tasks * esp_pestride)
       elseif (trim(esp_layout) == trim(layout_sequential)) then
          esp_inst_tasks = esp_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid esp_layout ')
       endif
       current_task_rootpe = esp_rootpe
       do n = 1, num_inst_esp
          emin(n) = current_task_rootpe
          emax(n) = current_task_rootpe &
                    + ((esp_inst_tasks - 1) * esp_pestride)
          estr(n) = esp_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do

       !! Coupler tasks

       cmin = cpl_rootpe
       cmax = cpl_rootpe + (cpl_ntasks-1)*cpl_pestride
       cstr = cpl_pestride
    end if

    call shr_mpi_bcast(atm_nthreads,GLOBAL_COMM,'atm_nthreads')
    call shr_mpi_bcast(lnd_nthreads,GLOBAL_COMM,'lnd_nthreads')
    call shr_mpi_bcast(ocn_nthreads,GLOBAL_COMM,'ocn_nthreads')
    call shr_mpi_bcast(ice_nthreads,GLOBAL_COMM,'ice_nthreads')
    call shr_mpi_bcast(glc_nthreads,GLOBAL_COMM,'glc_nthreads')
    call shr_mpi_bcast(wav_nthreads,GLOBAL_COMM,'wav_nthreads')
    call shr_mpi_bcast(rof_nthreads,GLOBAL_COMM,'rof_nthreads')
    call shr_mpi_bcast(esp_nthreads,GLOBAL_COMM,'esp_nthreads')
    call shr_mpi_bcast(cpl_nthreads,GLOBAL_COMM,'cpl_nthreads')

    ! Create MPI communicator groups

    if (mype == 0) then
       pelist(1,1) = 0
       pelist(2,1) = numpes-1
       pelist(3,1) = 1
    end if
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
    call seq_comm_setcomm(GLOID, pelist,iname='GLOBAL')

    if (mype == 0) then
       pelist(1,1) = cmin
       pelist(2,1) = cmax
       pelist(3,1) = cstr
    end if
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
    call seq_comm_setcomm(CPLID,pelist,cpl_nthreads,'CPL')

    do n = 1, num_inst_atm
       if (mype == 0) then
          pelist(1,1) = amin(n)
          pelist(2,1) = amax(n)
          pelist(3,1) = astr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(ATMID(n), pelist, atm_nthreads, 'ATM', n, num_inst_atm)
    end do

    do n = 1, num_inst_lnd
       if (mype == 0) then
          pelist(1,1) = lmin(n)
          pelist(2,1) = lmax(n)
          pelist(3,1) = lstr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(LNDID(n), pelist, lnd_nthreads, 'LND', n, num_inst_lnd)
    end do

    do n = 1, num_inst_ocn
       if (mype == 0) then
          pelist(1,1) = omin(n)
          pelist(2,1) = omax(n)
          pelist(3,1) = ostr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(OCNID(n), pelist, ocn_nthreads, 'OCN', n, num_inst_ocn)
    end do

    do n = 1, num_inst_ice
       if (mype == 0) then
          pelist(1,1) = imin(n)
          pelist(2,1) = imax(n)
          pelist(3,1) = istr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(ICEID(n), pelist, ice_nthreads, 'ICE', n, num_inst_ice)
    end do

    do n = 1, num_inst_glc
       if (mype == 0) then
          pelist(1,1) = gmin(n)
          pelist(2,1) = gmax(n)
          pelist(3,1) = gstr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(GLCID(n), pelist, glc_nthreads, 'GLC', n, num_inst_glc)
    end do

    do n = 1, num_inst_rof
       if (mype == 0) then
          pelist(1,1) = rmin(n)
          pelist(2,1) = rmax(n)
          pelist(3,1) = rstr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(ROFID(n), pelist, rof_nthreads, 'ROF', n, num_inst_rof)
    end do

    do n = 1, num_inst_wav
       if (mype == 0) then
          pelist(1,1) = wmin(n)
          pelist(2,1) = wmax(n)
          pelist(3,1) = wstr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(WAVID(n), pelist, wav_nthreads, 'WAV', n, num_inst_wav)
    end do

    do n = 1, num_inst_esp
       if (mype == 0) then
          pelist(1,1) = emin(n)
          pelist(2,1) = emax(n)
          pelist(3,1) = estr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(ESPID(n), pelist, esp_nthreads, 'ESP', n, num_inst_esp)
    end do

    !! Count the total number of threads

    max_threads = -1
    do n = 1,ncomps
       max_threads = max(max_threads,seq_comms(n)%nthreads)
    enddo
    do n = 1,ncomps
       seq_comms(n)%pethreads = max_threads
    enddo

    ! compute each components root pe global id and broadcast so all pes have info
    do n = 1,ncomps
       gloroot = -999
       call shr_mpi_max(gloroot,seq_comms(n)%gloroot,GLOBAL_COMM, trim(subname)//' gloroot',all=.true.)
    enddo

    !------------------------------------------
    ! Initialize MCT
    !------------------------------------------

    ! add up valid comps on local pe
    myncomps = 0
    do n = 1,ncomps
       if (seq_comms(n)%mpicom /= MPI_COMM_NULL) then
          myncomps = myncomps + 1
       endif
    enddo

    ! set comps and comms
    allocate(comps(myncomps),comms(myncomps),stat=ierr)
    if(ierr/=0) call mct_die(subName,'allocate comps comms',ierr)

    myncomps = 0
    do n = 1,ncomps
       if (seq_comms(n)%mpicom /= MPI_COMM_NULL) then
          myncomps = myncomps + 1
          if (myncomps > size(comps)) then
             write(logunit,*) trim(subname),' ERROR in myncomps ',myncomps,size(comps)
             call shr_sys_abort()
          endif
          comps(myncomps) = seq_comms(n)%ID
          comms(myncomps) = seq_comms(n)%mpicom
       endif
    enddo
    if (myncomps /= size(comps)) then
       write(logunit,*) trim(subname),' ERROR in myncomps ',myncomps,size(comps)
       call shr_sys_abort()
    endif

    call mct_world_init(ncomps, GLOBAL_COMM, comms, comps)

    deallocate(comps,comms)

    !------------------------------------------
    ! ESMF logging (only has effect if ESMF libraries are used)
    !------------------------------------------

    call mpi_bcast(esmf_logging, len(esmf_logging), MPI_CHARACTER, 0, GLOBAL_COMM, ierr)

    select case(esmf_logging)
    case ("ESMF_LOGKIND_SINGLE")
       esmf_logfile_kind = ESMF_LOGKIND_SINGLE
    case ("ESMF_LOGKIND_MULTI")
       esmf_logfile_kind = ESMF_LOGKIND_MULTI
    case ("ESMF_LOGKIND_NONE")
       esmf_logfile_kind = ESMF_LOGKIND_NONE
    case default
       if (mype == 0) then
          write(logunit,*) trim(subname),' ERROR: Invalid value for esmf_logging, ',esmf_logging
       endif
       call shr_sys_abort(trim(subname)//' ERROR: Invalid value for esmf_logging '//esmf_logging)
    end select

  end subroutine seq_comm_init

!---------------------------------------------------------
  subroutine seq_comm_setcomm(ID,pelist,nthreads,iname,inst,tinst, comm_in)
    use shr_sys_mod, only : shr_sys_abort
    use mpi, only : MPI_COMM_NULL, mpi_comm_group, mpi_comm_create, mpi_group_range_incl
    use mpi, only: mpi_comm_size, mpi_comm_rank
    use shr_mpi_mod, only : shr_mpi_chkerr
    implicit none
    integer,intent(IN) :: ID
    integer,intent(IN) :: pelist(:,:)
    integer,intent(IN),optional :: nthreads
    character(len=*),intent(IN),optional :: iname  ! name of component
    integer,intent(IN),optional :: inst  ! instance of component
    integer,intent(IN),optional :: tinst ! total number of instances for this component
    integer,intent(in),optional :: comm_in

    integer :: mpigrp_world
    integer :: mpigrp
    integer :: mpicom
    integer :: ntask,ntasks,cnt
    integer :: ierr
    character(len=seq_comm_namelen) :: cname
    logical :: set_suffix
    character(*),parameter :: subName =   '(seq_comm_setcomm) '

    if (ID < 1 .or. ID > ncomps) then
       write(logunit,*) subname,' ID out of range, abort ',ID
       call shr_sys_abort()
    endif
    if(present(comm_in)) then
       GLOBAL_COMM=comm_in
    endif

    call mpi_comm_group(GLOBAL_COMM, mpigrp_world, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_group mpigrp_world')
    call mpi_group_range_incl(mpigrp_world, 1, pelist, mpigrp,ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_group_range_incl mpigrp')
    call mpi_comm_create(GLOBAL_COMM, mpigrp, mpicom, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_create mpigrp')

    ntasks = ((pelist(2,1) - pelist(1,1)) / pelist(3,1)) + 1
    allocate(seq_comms(ID)%petlist(ntasks))
    seq_comms(ID)%petlist_allocated = .true.
    cnt = 0
    do ntask = pelist(1,1),pelist(2,1),pelist(3,1)
        cnt = cnt + 1
        if (cnt > ntasks) then
           write(logunit,*) subname,' ERROR in petlist init ',ntasks,pelist(1:3,1),ntask,cnt
           call shr_sys_abort(subname//' ERROR in petlist init')
        endif
        seq_comms(ID)%petlist(cnt) = ntask
    enddo

    seq_comms(ID)%set = .true.
    seq_comms(ID)%ID = ID

    if (present(inst)) then
       seq_comms(ID)%inst = inst
       set_suffix = .true.
    else
       seq_comms(ID)%inst = 1
       set_suffix = .false.
    endif

    if (present(tinst)) then
       if (tinst == 1) set_suffix = .false.
    endif

    if (present(iname)) then
       seq_comms(ID)%name = trim(iname)
       if (set_suffix) then
          call seq_comm_mkname(cname,iname,seq_comms(ID)%inst)
          seq_comms(ID)%name = trim(cname)
       endif
    endif

    if (set_suffix) then
       call seq_comm_mkname(cname,'_',seq_comms(ID)%inst)
       seq_comms(ID)%suffix = trim(cname)
    else
       seq_comms(ID)%suffix = ' '
    endif

    seq_comms(ID)%mpicom = mpicom
    seq_comms(ID)%mpigrp = mpigrp
    if (present(nthreads)) then
       seq_comms(ID)%nthreads = nthreads
    else
       seq_comms(ID)%nthreads = 1
    endif

    if (mpicom /= MPI_COMM_NULL) then
       call mpi_comm_size(mpicom,seq_comms(ID)%npes,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_size')
       call mpi_comm_rank(mpicom,seq_comms(ID)%iam,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank')
       if (seq_comms(ID)%iam == 0) then
          seq_comms(ID)%iamroot = .true.
       else
          seq_comms(ID)%iamroot = .false.
       endif
    else
       seq_comms(ID)%npes = -1
       seq_comms(ID)%iam = -1
       seq_comms(ID)%nthreads = 1
       seq_comms(ID)%iamroot = .false.
    endif

    if (seq_comms(ID)%iamroot) then
       write(logunit,F11) trim(subname),'  initialize ID ',ID,seq_comms(ID)%name, &
         ' pelist   =',pelist,' npes =',seq_comms(ID)%npes,' nthreads =',seq_comms(ID)%nthreads
    endif

  end subroutine seq_comm_setcomm

!---------------------------------------------------------
  subroutine seq_comm_printcomms()
    use shr_sys_mod, only : shr_sys_flush
    integer :: n
    character(*),parameter :: subName =   '(seq_comm_printcomms) '

    do n = 1,ncomps
       write(logunit,'(a,4i6,2x,3a)') trim(subName), n, &
            seq_comms(n)%gloroot, seq_comms(n)%npes, seq_comms(n)%nthreads, &
            trim(seq_comms(n)%name),':',trim(seq_comms(n)%suffix)
    enddo
    call shr_sys_flush(logunit)

  end subroutine seq_comm_printcomms

!---------------------------------------------------------

  subroutine seq_comm_setptrs(ID,mpicom,mpigrp,npes,nthreads,iam,iamroot,gloroot, pethreads, name)
    use mpi, only : mpi_comm_null, mpi_group_null
    implicit none
    integer,intent(in) :: ID
    integer,intent(out),optional :: mpicom
    integer,intent(out),optional :: mpigrp
    integer,intent(out),optional :: npes
    integer,intent(out),optional :: nthreads
    integer,intent(out),optional :: iam
    logical,intent(out),optional :: iamroot
    integer,intent(out),optional :: gloroot
    integer,intent(out),optional :: pethreads
    character(len=seq_comm_namelen)  , intent(out), optional :: name
    character(*),parameter :: subName =   '(seq_comm_setptrs) '

    ! Negative ID means there is no comm, return default or inactive values
    if ((ID == 0) .or. (ID > ncomps)) then
       write(logunit,*) subname,' ID out of range, return ',ID
       return
    endif

    if (present(mpicom)) then
       if (ID > 0) then
          mpicom = seq_comms(ID)%mpicom
       else
          mpicom = MPI_COMM_NULL
       end if
    endif

    if (present(mpigrp)) then
       if (ID > 0) then
          mpigrp = seq_comms(ID)%mpigrp
       else
          mpigrp = MPI_GROUP_NULL
       end if
    endif

    if (present(npes)) then
       if (ID > 0) then
          npes = seq_comms(ID)%npes
       else
          npes = 0
       end if
    endif

    if (present(nthreads)) then
       if (ID > 0) then
          nthreads = seq_comms(ID)%nthreads
       else
          nthreads = 1
       end if
    endif

    if (present(iam)) then
       if (ID > 0) then
          iam = seq_comms(ID)%iam
       else
          iam = -1
       end if
    endif

    if (present(iamroot)) then
       if (ID > 0) then
          iamroot = seq_comms(ID)%iamroot
       else
          iamroot = .false.
       end if
    endif

    if (present(gloroot)) then
       if (ID > 0) then
          gloroot = seq_comms(ID)%gloroot
       else
          gloroot = -1
       end if
    endif

    if (present(pethreads)) then
       if (ID > 0) then
          pethreads = seq_comms(ID)%pethreads
       else
          pethreads = 1
       end if
    endif

    if(present(name)) then
       if (ID > 0) then
          name = seq_comms(ID)%name
       else
          name = ''
       end if
    end if

  end subroutine seq_comm_setptrs
!---------------------------------------------------------
  subroutine seq_comm_setnthreads(nthreads)
    use shr_sys_mod, only : shr_sys_abort

    implicit none
    integer,intent(in) :: nthreads
    character(*),parameter :: subName =   '(seq_comm_setnthreads) '

#ifdef _OPENMP
    if (nthreads < 1) then
       call shr_sys_abort(subname//' ERROR: nthreads less than one')
    endif
    call omp_set_num_threads(nthreads)
#endif

  end subroutine seq_comm_setnthreads
!---------------------------------------------------------
  integer function seq_comm_getnthreads()

    implicit none
    integer :: omp_get_num_threads
    character(*),parameter :: subName =   '(seq_comm_getnthreads) '

    seq_comm_getnthreads = -1
#ifdef _OPENMP
!$OMP PARALLEL
    seq_comm_getnthreads = omp_get_num_threads()
!$OMP END PARALLEL
#endif

  end function seq_comm_getnthreads
!---------------------------------------------------------
  logical function seq_comm_iamin(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_iamin) '

    if ((ID < 1) .or. (ID > ncomps)) then
      seq_comm_iamin = .false.
    else if (seq_comms(ID)%iam >= 0) then
       seq_comm_iamin = .true.
    else
       seq_comm_iamin = .false.
    endif

  end function seq_comm_iamin
!---------------------------------------------------------
  logical function seq_comm_iamroot(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_iamroot) '

    if ((ID < 1) .or. (ID > ncomps)) then
       seq_comm_iamroot = .false.
    else
       seq_comm_iamroot = seq_comms(ID)%iamroot
    end if

  end function seq_comm_iamroot
!---------------------------------------------------------
  integer function seq_comm_mpicom(ID)
    use mpi, only : mpi_comm_null
    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_mpicom) '

    if ((ID < 1) .or. (ID > ncomps)) then
       seq_comm_mpicom = MPI_COMM_NULL
    else
       seq_comm_mpicom = seq_comms(ID)%mpicom
    end if

  end function seq_comm_mpicom
!---------------------------------------------------------
  integer function seq_comm_iam(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_iam) '

    if ((ID < 1) .or. (ID > ncomps)) then
       seq_comm_iam = -1
    else
       seq_comm_iam = seq_comms(ID)%iam
    end if

  end function seq_comm_iam

!---------------------------------------------------------
  integer function seq_comm_gloroot(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_gloroot) '

    if ((ID < 1) .or. (ID > ncomps)) then
       seq_comm_gloroot = -1
    else
       seq_comm_gloroot = seq_comms(ID)%gloroot
    end if

  end function seq_comm_gloroot

!---------------------------------------------------------
  character(len=seq_comm_namelen) function seq_comm_name(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_name) '

    if ((ID < 1) .or. (ID > ncomps)) then
       seq_comm_name = ''
    else
       seq_comm_name = trim(seq_comms(ID)%name)
    end if

  end function seq_comm_name
!---------------------------------------------------------
  character(len=seq_comm_namelen) function seq_comm_suffix(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_suffix) '

    if ((ID < 1) .or. (ID > ncomps)) then
       seq_comm_suffix = ''
    else
       seq_comm_suffix = trim(seq_comms(ID)%suffix)
    end if

  end function seq_comm_suffix
!---------------------------------------------------------
  subroutine seq_comm_petlist(ID,petlist)

    implicit none
    integer,intent(in) :: ID
    integer,pointer :: petlist(:)
    character(*),parameter :: subName =   '(seq_comm_petlist) '

    if ((ID < 1) .or. (ID > ncomps)) then
       nullify(petlist)
    else
       petlist => seq_comms(ID)%petlist
    end if

  end subroutine seq_comm_petlist
!---------------------------------------------------------
  integer function seq_comm_inst(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_inst) '

    if ((ID < 1) .or. (ID > ncomps)) then
      seq_comm_inst = 0
    else
      seq_comm_inst = seq_comms(ID)%inst
    end if

  end function seq_comm_inst
!---------------------------------------------------------
  subroutine seq_comm_mkname(oname,str1,num)
    use shr_sys_mod, only : shr_sys_abort
    implicit none
    character(len=*),intent(out) :: oname
    character(len=*),intent(in)  :: str1
    integer,intent(in)           :: num
    character(*),parameter :: subName =   '(seq_comm_mkname) '

    character(len=8) :: cnum

    write(cnum,'(i4.4)') num
    if (len_trim(str1) + len_trim(cnum) > len(oname)) then
       write(logunit,*) trim(subname),' ERROR in str lens ',len(oname),trim(str1),trim(cnum)
       call shr_sys_abort(trim(subname))
    endif
    oname = trim(str1)//trim(cnum)

  end subroutine seq_comm_mkname
!---------------------------------------------------------
end module seq_comm_mct
