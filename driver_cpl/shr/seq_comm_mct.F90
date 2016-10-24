module seq_comm_mct

!---------------------------------------------------------------------
!
! Purpose: MCT utitlity functions used in sequential CCSM.
!          Note that if no MPI, will call MCTs fake version
!          (including mpif.h) will be utilized
!
!---------------------------------------------------------------------


!!! NOTE: If all atmospheres are identical in number of processes,
!!! number of threads, and grid layout, we should check that the
!!! user-provided number of processes and threads are consistent
!!! (or else, only accept one entry for these quantities when reading
!!! the namelist).  ARE OTHER PROTECTIONS/CHECKS NEEDED???


  use mct_mod     , only : mct_world_init, mct_world_clean, mct_die
  use shr_sys_mod , only : shr_sys_abort, shr_sys_flush
  use shr_mpi_mod , only : shr_mpi_chkerr, shr_mpi_bcast, shr_mpi_max
  use shr_file_mod, only : shr_file_getUnit, shr_file_freeUnit
#ifdef USE_ESMF_LIB
  use esmf
#endif

  implicit none

  private
#include <mpif.h>  
  save

!--------------------------------------------------------------------------
! Public interfaces
!--------------------------------------------------------------------------

  public seq_comm_init
  public seq_comm_clean
  public seq_comm_iamin
  public seq_comm_iamroot
  public seq_comm_mpicom
  public seq_comm_iam
  public seq_comm_gloiam
  public seq_comm_gloroot
  public seq_comm_cplpe
  public seq_comm_cmppe
  public seq_comm_name
  public seq_comm_inst
  public seq_comm_suffix
  public seq_comm_petlist
  public seq_comm_setptrs
  public seq_comm_setnthreads
  public seq_comm_getnthreads
  public seq_comm_printcomms
  public seq_comm_get_ncomps
#ifdef USE_ESMF_LIB
  public seq_comm_getcompstates
  public seq_comm_setcompstates
#endif

!--------------------------------------------------------------------------
! Public data
!--------------------------------------------------------------------------

  integer, public, parameter :: default_logunit = 6
  integer, public :: logunit  = default_logunit     ! log unit number
  integer, public :: loglevel = 1     ! log level

  integer, public :: global_mype = -1  !! To be initialized

  !!! Note - NUM_COMP_INST_XXX are cpp variables set in buildlib.csm_share

  integer, parameter :: nphysmod = 7  ! number of physical models
  integer, parameter, public :: num_inst_atm = NUM_COMP_INST_ATM
  integer, parameter, public :: num_inst_lnd = NUM_COMP_INST_LND
  integer, parameter, public :: num_inst_ocn = NUM_COMP_INST_OCN
  integer, parameter, public :: num_inst_ice = NUM_COMP_INST_ICE
  integer, parameter, public :: num_inst_glc = NUM_COMP_INST_GLC
  integer, parameter, public :: num_inst_wav = NUM_COMP_INST_WAV
  integer, parameter, public :: num_inst_rof = NUM_COMP_INST_ROF

  integer, parameter, public :: num_inst_total= num_inst_atm + &
                                                num_inst_lnd + &
                                                num_inst_ocn + &
                                                num_inst_ice + &
                                                num_inst_glc + &
                                                num_inst_wav + &
                                                num_inst_rof + 1

  integer, public :: num_inst_min, num_inst_max
  integer, public :: num_inst_xao    ! for xao flux
  integer, public :: num_inst_frc    ! for fractions

  !!! Each component instance needs two communicators: one internal to the
  !!! instance, and one for communicating with the coupler.
  !!! Additionally, one communicator is needed for the coupler's
  !!! internal communications, and one is needed for the global space.

  integer, parameter, public :: num_inst_phys = num_inst_atm + num_inst_lnd + &
                                                num_inst_ocn + num_inst_ice + &
                                                num_inst_glc + num_inst_rof + &
                                                num_inst_wav
  integer, parameter :: ncomps = (2 + 2*nphysmod + (2 * num_inst_phys))

  integer, public :: GLOID
  integer, public :: CPLID

  integer, public :: ALLATMID
  integer, public :: ALLLNDID
  integer, public :: ALLOCNID
  integer, public :: ALLICEID
  integer, public :: ALLGLCID
  integer, public :: ALLROFID
  integer, public :: ALLWAVID

  integer, public :: CPLALLATMID
  integer, public :: CPLALLLNDID
  integer, public :: CPLALLOCNID
  integer, public :: CPLALLICEID
  integer, public :: CPLALLGLCID
  integer, public :: CPLALLROFID
  integer, public :: CPLALLWAVID

  integer, public :: ATMID(num_inst_atm)
  integer, public :: LNDID(num_inst_lnd)
  integer, public :: OCNID(num_inst_ocn)
  integer, public :: ICEID(num_inst_ice)
  integer, public :: GLCID(num_inst_glc)
  integer, public :: ROFID(num_inst_rof)
  integer, public :: WAVID(num_inst_wav)

  integer, public :: CPLATMID(num_inst_atm)
  integer, public :: CPLLNDID(num_inst_lnd)
  integer, public :: CPLOCNID(num_inst_ocn)
  integer, public :: CPLICEID(num_inst_ice)
  integer, public :: CPLGLCID(num_inst_glc)
  integer, public :: CPLROFID(num_inst_rof)
  integer, public :: CPLWAVID(num_inst_wav)

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
    integer :: gloiam          ! my task number in mpi_comm_world
    integer :: gloroot         ! the global task number of each comps root on all pes
    integer :: pethreads       ! max number of threads on my task
    integer :: cplpe           ! a common task in mpicom from the cpl group for join mpicoms
    integer :: cmppe           ! a common task in mpicom from the component group for join mpicoms
    logical :: set             ! has this datatype been set
    integer, pointer    :: petlist(:)  ! esmf pet list
    logical :: petlist_allocated ! whether the petlist pointer variable was allocated
#ifdef USE_ESMF_LIB
    type(ESMF_GridComp) :: esmf_comp   ! esmf gridded component
                                       ! The following state members are not needed in 520r.
    type(ESMF_State)    :: imp_state   ! esmf import state for the gridded component
    type(ESMF_State)    :: exp_state   ! esmf export state for the gridded component
#endif
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
       ocn_layout, wav_layout
  
  logical :: seq_comm_mct_initialized = .false.  ! whether this module has been initialized

!=======================================================================
contains
!======================================================================
  integer function seq_comm_get_ncomps() 
    seq_comm_get_ncomps = ncomps
  end function seq_comm_get_ncomps

  
  subroutine seq_comm_init(Comm_in, nmlfile)
      
    !----------------------------------------------------------
    !
    ! Arguments
    implicit none
    integer, intent(in) :: Comm_in
    character(len=*), intent(IN) :: nmlfile
    !
    ! Local variables
    !
    logical :: error_state
    integer :: ierr, n, count
    character(*), parameter :: subName =   '(seq_comm_init) '
    integer :: mpi_group_world   ! MPI_COMM_WORLD group
    integer :: mype,numpes,myncomps,max_threads,gloroot
    integer :: atm_inst_tasks, lnd_inst_tasks, ocn_inst_tasks, ice_inst_tasks, &
               glc_inst_tasks, rof_inst_tasks, wav_inst_tasks
    integer :: current_task_rootpe, droot
    integer :: amin(num_inst_atm), amax(num_inst_atm), astr(num_inst_atm)
    integer :: lmin(num_inst_lnd), lmax(num_inst_lnd), lstr(num_inst_lnd)
    integer :: imin(num_inst_ice), imax(num_inst_ice), istr(num_inst_ice)
    integer :: omin(num_inst_ocn), omax(num_inst_ocn), ostr(num_inst_ocn)
    integer :: gmin(num_inst_glc), gmax(num_inst_glc), gstr(num_inst_glc)
    integer :: wmin(num_inst_wav), wmax(num_inst_wav), wstr(num_inst_wav)
    integer :: rmin(num_inst_rof), rmax(num_inst_rof), rstr(num_inst_rof)
    integer :: cmin,cmax,cstr
    integer :: pelist(3,1)       ! start, stop, stride for group
    integer, pointer :: comps(:) ! array with component ids
    integer, pointer :: comms(:) ! array with mpicoms
    integer :: nu, i

    integer :: &
         atm_ntasks, atm_rootpe, atm_pestride, atm_nthreads, &
         lnd_ntasks, lnd_rootpe, lnd_pestride, lnd_nthreads, &
         ice_ntasks, ice_rootpe, ice_pestride, ice_nthreads, &
         glc_ntasks, glc_rootpe, glc_pestride, glc_nthreads, &
         wav_ntasks, wav_rootpe, wav_pestride, wav_nthreads, &
         rof_ntasks, rof_rootpe, rof_pestride, rof_nthreads, &
         ocn_ntasks, ocn_rootpe, ocn_pestride, ocn_nthreads, &
         cpl_ntasks, cpl_rootpe, cpl_pestride, cpl_nthreads
    namelist /ccsm_pes/  &
         atm_ntasks, atm_rootpe, atm_pestride, atm_nthreads, atm_layout, &
         lnd_ntasks, lnd_rootpe, lnd_pestride, lnd_nthreads, lnd_layout, &
         ice_ntasks, ice_rootpe, ice_pestride, ice_nthreads, ice_layout, &
         glc_ntasks, glc_rootpe, glc_pestride, glc_nthreads, glc_layout, &
         wav_ntasks, wav_rootpe, wav_pestride, wav_nthreads, wav_layout, &
         rof_ntasks, rof_rootpe, rof_pestride, rof_nthreads, rof_layout, &
         ocn_ntasks, ocn_rootpe, ocn_pestride, ocn_nthreads, ocn_layout, &
         cpl_ntasks, cpl_rootpe, cpl_pestride, cpl_nthreads 
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
       seq_comms(n)%gloiam = -1
       seq_comms(n)%gloroot = -1
       seq_comms(n)%pethreads = -1
       seq_comms(n)%cplpe = -1
       seq_comms(n)%cmppe = -1
    enddo


    ! Initialize MPI
    ! Note that if no MPI, will call MCTs fake version

    call mpi_comm_rank(GLOBAL_COMM, mype  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank comm_world')
    call mpi_comm_size(GLOBAL_COMM, numpes, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size comm_world')

    ! Initialize gloiam on all IDs

    global_mype = mype
 
    do n = 1,ncomps
       seq_comms(n)%gloiam = mype
    enddo

    ! Set ntasks, rootpe, pestride, nthreads for all components

    if (mype == 0) then

       !! Set up default atmosphere process parameters

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

       cpl_ntasks = numpes
       cpl_rootpe = 0
       cpl_pestride = 1
       cpl_nthreads = 1

       ! Read namelist if it exists

       nu = shr_file_getUnit()
       open(nu, file=trim(nmlfile), status='old', iostat=ierr)

       if (ierr == 0) then
          ierr = 1
          do while( ierr > 0 )
             read(nu, nml=ccsm_pes, iostat=ierr)
          end do
          close(nu)
       end if
       call shr_file_freeUnit(nu)

    end if

    !--- compute some other num_inst values
       
    num_inst_xao = max(num_inst_atm,num_inst_ocn)
    num_inst_frc = num_inst_ice

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
    num_inst_max = num_inst_atm
    num_inst_max = max(num_inst_max, num_inst_lnd)
    num_inst_max = max(num_inst_max, num_inst_ocn)
    num_inst_max = max(num_inst_max, num_inst_ice)
    num_inst_max = max(num_inst_max, num_inst_glc)
    num_inst_max = max(num_inst_max, num_inst_wav)
    num_inst_max = max(num_inst_max, num_inst_rof)

    if (num_inst_min /= num_inst_max .and. num_inst_min /= 1) error_state = .true.
    if (num_inst_atm /= num_inst_min .and. num_inst_atm /= num_inst_max) error_state = .true.
    if (num_inst_lnd /= num_inst_min .and. num_inst_lnd /= num_inst_max) error_state = .true.
    if (num_inst_ocn /= num_inst_min .and. num_inst_ocn /= num_inst_max) error_state = .true.
    if (num_inst_ice /= num_inst_min .and. num_inst_ice /= num_inst_max) error_state = .true.
    if (num_inst_glc /= num_inst_min .and. num_inst_glc /= num_inst_max) error_state = .true.
    if (num_inst_wav /= num_inst_min .and. num_inst_wav /= num_inst_max) error_state = .true.
    if (num_inst_rof /= num_inst_min .and. num_inst_rof /= num_inst_max) error_state = .true.

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

    count = count + 1
    ALLATMID = count
    count = count + 1
    ALLLNDID = count
    count = count + 1
    ALLOCNID = count
    count = count + 1
    ALLICEID = count
    count = count + 1
    ALLGLCID = count
    count = count + 1
    ALLROFID = count
    count = count + 1
    ALLWAVID = count

    count = count + 1
    CPLALLATMID = count
    count = count + 1
    CPLALLLNDID = count
    count = count + 1
    CPLALLOCNID = count
    count = count + 1
    CPLALLICEID = count
    count = count + 1
    CPLALLGLCID = count
    count = count + 1
    CPLALLROFID = count
    count = count + 1
    CPLALLWAVID = count

    do n = 1, num_inst_atm
       count = count + 1
       ATMID(n) = count
       count = count + 1
       CPLATMID(n) = count
    end do
         
    do n = 1, num_inst_lnd
       count = count + 1
       LNDID(n) = count
       count = count + 1
       CPLLNDID(n) = count
    end do
       
    do n = 1, num_inst_ocn
       count = count + 1
       OCNID(n) = count
       count = count + 1
       CPLOCNID(n) = count
    end do
       
    do n = 1, num_inst_ice
       count = count + 1
       ICEID(n) = count
       count = count + 1
       CPLICEID(n) = count
    end do

    do n = 1, num_inst_glc
       count = count + 1
       GLCID(n) = count
       count = count + 1
       CPLGLCID(n) = count
    end do

    do n = 1, num_inst_rof
       count = count + 1
       ROFID(n) = count
       count = count + 1
       CPLROFID(n) = count
    end do

    do n = 1, num_inst_wav
       count = count + 1
       WAVID(n) = count
       count = count + 1
       CPLWAVID(n) = count
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
       if (cpl_rootpe < 0) error_state = .true.
       
       if (error_state) then
          write(logunit,*) trim(subname),' ERROR: rootpes must be >= 0'
          call shr_sys_abort(trim(subname)//' ERROR: rootpes >= 0')
       endif

!       ! nthreads = 1, temporary
!       if (atm_nthreads /= 1 .or. lnd_nthreads /= 1 .or. ice_nthreads /= 1 .or. &
!           ocn_nthreads /= 1 .or. cpl_nthreads /= 1) then
!          write(logunit,*) trim(subname),' ERROR: nthreads must be 1'
!          call shr_sys_abort()
!       endif

!       ! nthreads should be 1 or something consistent, compute max nthreads
!       amax = max(atm_nthreads,lnd_nthreads)
!       amax = max(amax        ,ice_nthreads)
!       amax = max(amax        ,ocn_nthreads)
!       amax = max(amax        ,cpl_nthreads)

!       ! check that everything is either 1 or max nthreads
!       if ((atm_nthreads /= 1 .and. atm_nthreads /= amax) .or. &
!           (lnd_nthreads /= 1 .and. lnd_nthreads /= amax) .or. &
!           (ice_nthreads /= 1 .and. ice_nthreads /= amax) .or. &
!           (ocn_nthreads /= 1 .and. ocn_nthreads /= amax) .or. &
!           (cpl_nthreads /= 1 .and. cpl_nthreads /= amax)) then
!          write(logunit,*) trim(subname),' ERROR: nthreads must be consistent'
!          call shr_sys_abort()
!       endif

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

       !! Coupler tasks

       cmin = cpl_rootpe
       cmax = cpl_rootpe + (cpl_ntasks-1)*cpl_pestride
       cstr = cpl_pestride
    end if

#if (1 == 0)
    ! create petlist for ESMF components, doesn't work for ensembles 
    if(present(atm_petlist)) then
        call shr_mpi_bcast(atm_ntasks, GLOBAL_COMM, 'atm_ntasks')
        call shr_mpi_bcast(atm_rootpe, GLOBAL_COMM, 'atm_rootpe')
        call shr_mpi_bcast(atm_pestride, GLOBAL_COMM, 'atm_pestride')
        allocate(atm_petlist(atm_ntasks))
        do i = 1, atm_ntasks
            atm_petlist(i) = atm_rootpe + (i-1)*atm_pestride
        enddo
    endif

    if(present(lnd_petlist)) then
        call shr_mpi_bcast(lnd_ntasks, GLOBAL_COMM, 'lnd_ntasks')
        call shr_mpi_bcast(lnd_rootpe, GLOBAL_COMM, 'lnd_rootpe')
        call shr_mpi_bcast(lnd_pestride, GLOBAL_COMM, 'lnd_pestride')
        allocate(lnd_petlist(lnd_ntasks))
        do i = 1, lnd_ntasks
            lnd_petlist(i) = lnd_rootpe + (i-1)*lnd_pestride
        enddo
    endif

    if(present(ice_petlist)) then
        call shr_mpi_bcast(ice_ntasks, GLOBAL_COMM, 'ice_ntasks')
        call shr_mpi_bcast(ice_rootpe, GLOBAL_COMM, 'ice_rootpe')
        call shr_mpi_bcast(ice_pestride, GLOBAL_COMM, 'ice_pestride')
        allocate(ice_petlist(ice_ntasks))
        do i = 1, ice_ntasks
            ice_petlist(i) = ice_rootpe + (i-1)*ice_pestride
        enddo
    endif

    if(present(ocn_petlist)) then
        call shr_mpi_bcast(ocn_ntasks, GLOBAL_COMM, 'ocn_ntasks')
        call shr_mpi_bcast(ocn_rootpe, GLOBAL_COMM, 'ocn_rootpe')
        call shr_mpi_bcast(ocn_pestride, GLOBAL_COMM, 'ocn_pestride')
        allocate(ocn_petlist(ocn_ntasks))
        do i = 1, ocn_ntasks
            ocn_petlist(i) = ocn_rootpe + (i-1)*ocn_pestride
        enddo
    endif

    if(present(glc_petlist)) then
        call shr_mpi_bcast(glc_ntasks, GLOBAL_COMM, 'glc_ntasks')
        call shr_mpi_bcast(glc_rootpe, GLOBAL_COMM, 'glc_rootpe')
        call shr_mpi_bcast(glc_pestride, GLOBAL_COMM, 'glc_pestride')
        allocate(glc_petlist(glc_ntasks))
        do i = 1, glc_ntasks
            glc_petlist(i) = glc_rootpe + (i-1)*glc_pestride
        enddo
    endif

    if(present(rof_petlist)) then
        call shr_mpi_bcast(rof_ntasks, GLOBAL_COMM, 'rof_ntasks')
        call shr_mpi_bcast(rof_rootpe, GLOBAL_COMM, 'rof_rootpe')
        call shr_mpi_bcast(rof_pestride, GLOBAL_COMM, 'rof_pestride')
        allocate(rof_petlist(rof_ntasks))
        do i = 1, rof_ntasks
            rof_petlist(i) = rof_rootpe + (i-1)*rof_pestride
        enddo
    endif

    if(present(wav_petlist)) then
        call shr_mpi_bcast(wav_ntasks, GLOBAL_COMM, 'wav_ntasks')
        call shr_mpi_bcast(wav_rootpe, GLOBAL_COMM, 'wav_rootpe')
        call shr_mpi_bcast(wav_pestride, GLOBAL_COMM, 'wav_pestride')
        allocate(wav_petlist(wav_ntasks))
        do i = 1, wav_ntasks
            wav_petlist(i) = wav_rootpe + (i-1)*wav_pestride
        enddo
    endif
#endif
       
    call shr_mpi_bcast(atm_nthreads,GLOBAL_COMM,'atm_nthreads')
    call shr_mpi_bcast(lnd_nthreads,GLOBAL_COMM,'lnd_nthreads')
    call shr_mpi_bcast(ocn_nthreads,GLOBAL_COMM,'ocn_nthreads')
    call shr_mpi_bcast(ice_nthreads,GLOBAL_COMM,'ice_nthreads')
    call shr_mpi_bcast(glc_nthreads,GLOBAL_COMM,'glc_nthreads')
    call shr_mpi_bcast(wav_nthreads,GLOBAL_COMM,'wav_nthreads')
    call shr_mpi_bcast(rof_nthreads,GLOBAL_COMM,'rof_nthreads')
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
       call seq_comm_joincomm(CPLID, ATMID(n), CPLATMID(n), 'CPLATM', n, num_inst_atm)
    end do
    call seq_comm_jcommarr(ATMID,ALLATMID,'ALLATMID',1,1)
    call seq_comm_joincomm(CPLID,ALLATMID,CPLALLATMID,'CPLALLATMID',1,1)

    do n = 1, num_inst_lnd
       if (mype == 0) then
          pelist(1,1) = lmin(n)
          pelist(2,1) = lmax(n)
          pelist(3,1) = lstr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(LNDID(n), pelist, lnd_nthreads, 'LND', n, num_inst_lnd)
       call seq_comm_joincomm(CPLID, LNDID(n), CPLLNDID(n), 'CPLLND', n, num_inst_lnd)
    end do
    call seq_comm_jcommarr(LNDID,ALLLNDID,'ALLLNDID',1,1)
    call seq_comm_joincomm(CPLID,ALLLNDID,CPLALLLNDID,'CPLALLLNDID',1,1)

    do n = 1, num_inst_ocn
       if (mype == 0) then
          pelist(1,1) = omin(n)
          pelist(2,1) = omax(n)
          pelist(3,1) = ostr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(OCNID(n), pelist, ocn_nthreads, 'OCN', n, num_inst_ocn)
       call seq_comm_joincomm(CPLID, OCNID(n), CPLOCNID(n), 'CPLOCN', n, num_inst_ocn)
    end do
    call seq_comm_jcommarr(OCNID,ALLOCNID,'ALLOCNID',1,1)
    call seq_comm_joincomm(CPLID,ALLOCNID,CPLALLOCNID,'CPLALLOCNID',1,1)

    do n = 1, num_inst_ice
       if (mype == 0) then
          pelist(1,1) = imin(n)
          pelist(2,1) = imax(n)
          pelist(3,1) = istr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(ICEID(n), pelist, ice_nthreads, 'ICE', n, num_inst_ice)
       call seq_comm_joincomm(CPLID, ICEID(n), CPLICEID(n), 'CPLICE', n, num_inst_ice)
    end do
    call seq_comm_jcommarr(ICEID,ALLICEID,'ALLICEID',1,1)
    call seq_comm_joincomm(CPLID,ALLICEID,CPLALLICEID,'CPLALLICEID',1,1)

    do n = 1, num_inst_glc
       if (mype == 0) then
          pelist(1,1) = gmin(n)
          pelist(2,1) = gmax(n)
          pelist(3,1) = gstr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(GLCID(n), pelist, glc_nthreads, 'GLC', n, num_inst_glc)
       call seq_comm_joincomm(CPLID, GLCID(n), CPLGLCID(n), 'CPLGLC', n, num_inst_glc)
    end do
    call seq_comm_jcommarr(GLCID,ALLGLCID,'ALLGLCID',1,1)
    call seq_comm_joincomm(CPLID,ALLGLCID,CPLALLGLCID,'CPLALLGLCID',1,1)

    do n = 1, num_inst_rof
       if (mype == 0) then
          pelist(1,1) = rmin(n)
          pelist(2,1) = rmax(n)
          pelist(3,1) = rstr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(ROFID(n), pelist, rof_nthreads, 'ROF', n, num_inst_rof)
       call seq_comm_joincomm(CPLID, ROFID(n), CPLROFID(n), 'CPLROF', n, num_inst_rof)
    end do
    call seq_comm_jcommarr(ROFID,ALLROFID,'ALLROFID',1,1)
    call seq_comm_joincomm(CPLID,ALLROFID,CPLALLROFID,'CPLALLROFID',1,1)

    do n = 1, num_inst_wav
       if (mype == 0) then
          pelist(1,1) = wmin(n)
          pelist(2,1) = wmax(n)
          pelist(3,1) = wstr(n)
       end if
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, GLOBAL_COMM, ierr)
       call seq_comm_setcomm(WAVID(n), pelist, wav_nthreads, 'WAV', n, num_inst_wav)
       call seq_comm_joincomm(CPLID, WAVID(n), CPLWAVID(n), 'CPLWAV', n, num_inst_wav)
    end do
    call seq_comm_jcommarr(WAVID,ALLWAVID,'ALLWAVID',1,1)
    call seq_comm_joincomm(CPLID,ALLWAVID,CPLALLWAVID,'CPLALLWAVID',1,1)

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
       if (seq_comms(n)%iamroot) gloroot = seq_comms(n)%gloiam
       call shr_mpi_max(gloroot,seq_comms(n)%gloroot,GLOBAL_COMM, &
                        trim(subname)//' gloroot',all=.true.)
    enddo

    ! Initialize MCT

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

    call seq_comm_printcomms()

  end subroutine seq_comm_init

!---------------------------------------------------------
  subroutine seq_comm_clean()
    ! Resets this module - freeing memory, etc.
    !
    ! This potentially allows seq_comm_init can be called again, e.g., from unit tests.
    !
    ! Also calls mct_world_clean, to be symmetric with the mct_world_init call from
    ! seq_comm_init.

    integer :: id
    
    character(*), parameter :: subName =   '(seq_comm_clean) '
    !----------------------------------------------------------
    
    if (.not. seq_comm_mct_initialized) then
       write(logunit,*) trim(subname),' ERROR seq_comm_init has not been called '
       call shr_sys_abort()
    end if
    seq_comm_mct_initialized = .false.

    do id = 1, ncomps
       if (seq_comms(id)%petlist_allocated) then
          deallocate(seq_comms(id)%petlist)
       end if
    end do
    
    call mct_world_clean()

  end subroutine seq_comm_clean
  
!---------------------------------------------------------
  subroutine seq_comm_setcomm(ID,pelist,nthreads,iname,inst,tinst)

    implicit none
    integer,intent(IN) :: ID
    integer,intent(IN) :: pelist(:,:)
    integer,intent(IN),optional :: nthreads
    character(len=*),intent(IN),optional :: iname  ! name of component
    integer,intent(IN),optional :: inst  ! instance of component
    integer,intent(IN),optional :: tinst ! total number of instances for this component

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
  subroutine seq_comm_joincomm(ID1,ID2,ID,iname,inst,tinst)

    implicit none
    integer,intent(IN) :: ID1    ! src id
    integer,intent(IN) :: ID2    ! srd id
    integer,intent(IN) :: ID     ! computed id
    character(len=*),intent(IN),optional :: iname  ! comm name
    integer,intent(IN),optional :: inst
    integer,intent(IN),optional :: tinst

    integer :: mpigrp
    integer :: mpicom
    integer :: ierr
    integer :: n,nsize
    character(len=seq_comm_namelen) :: cname
    logical :: set_suffix
    integer,allocatable :: pe_t1(:),pe_t2(:)
    character(*),parameter :: subName =   '(seq_comm_joincomm) '

    ! check that IDs are in valid range, that ID1 and ID2 have
    ! been set, and that ID has not been set

    if (ID1 < 1 .or. ID1 > ncomps) then
       write(logunit,*) subname,' ID1 out of range, abort ',ID1
       call shr_sys_abort()
    endif 
    if (ID2 < 1 .or. ID2 > ncomps) then
       write(logunit,*) subname,' ID2 out of range, abort ',ID2
       call shr_sys_abort()
    endif 
    if (ID < 1 .or. ID > ncomps) then
       write(logunit,*) subname,' ID out of range, abort ',ID
       call shr_sys_abort()
    endif
    if (.not. seq_comms(ID1)%set .or. .not. seq_comms(ID2)%set) then
       write(logunit,*) subname,' ID1 or ID2 not set ',ID1,ID2
       call shr_sys_abort()
    endif
    if (seq_comms(ID)%set) then
       write(logunit,*) subname,' ID already set ',ID
       call shr_sys_abort()
    endif

    call mpi_group_union(seq_comms(ID1)%mpigrp,seq_comms(ID2)%mpigrp,mpigrp,ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_union mpigrp')
    call mpi_comm_create(GLOBAL_COMM, mpigrp, mpicom, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_create mpigrp')

    seq_comms(ID)%set = .true.
    seq_comms(ID)%ID = ID

    if (present(inst)) then
       seq_comms(ID)%inst = inst
    else
       seq_comms(ID)%inst = 1
    endif

    set_suffix = .true.
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
    seq_comms(ID)%nthreads = max(seq_comms(ID1)%nthreads,seq_comms(ID2)%nthreads)
    seq_comms(ID)%nthreads = max(seq_comms(ID)%nthreads,1)

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
       seq_comms(ID)%iamroot = .false.
    endif

! needs to be excluded until mpi_group_size is added to serial mpi in mct
#if (1 == 0)
    if (loglevel > 3) then
       ! some debug code to prove the join is working ok
       ! when joining mpicomms, the local rank may be quite different
       !   from either the global or local ranks of the joining comms
       call mpi_group_size(seq_comms(ID1)%mpigrp,nsize,ierr)
       allocate(pe_t1(nsize),pe_t2(nsize))
       do n = 1,nsize
          pe_t1(n) = n-1
          pe_t2(n) = -1
       enddo
       call mpi_group_translate_ranks(seq_comms(ID1)%mpigrp, nsize, pe_t1, mpigrp, pe_t2, ierr)
       write(logunit,*) 'ID1      ranks ',pe_t1
       write(logunit,*) 'ID1-JOIN ranks ',pe_t2
       deallocate(pe_t1,pe_t2)

       call mpi_group_size(seq_comms(ID2)%mpigrp,nsize,ierr)
       allocate(pe_t1(nsize),pe_t2(nsize))
       do n = 1,nsize
          pe_t1(n) = n-1
          pe_t2(n) = -1
       enddo
       call mpi_group_translate_ranks(seq_comms(ID2)%mpigrp, nsize, pe_t1, mpigrp, pe_t2, ierr)
       write(logunit,*) 'ID2      ranks ',pe_t1
       write(logunit,*) 'ID2-JOIN ranks ',pe_t2
       deallocate(pe_t1,pe_t2)
    endif
#endif

    allocate(pe_t1(1),pe_t2(1))
    pe_t1(1) = 0
    call mpi_group_translate_ranks(seq_comms(ID1)%mpigrp, 1, pe_t1, mpigrp, pe_t2, ierr)
    seq_comms(ID)%cplpe = pe_t2(1)
    pe_t1(1) = 0
    call mpi_group_translate_ranks(seq_comms(ID2)%mpigrp, 1, pe_t1, mpigrp, pe_t2, ierr)
    seq_comms(ID)%cmppe = pe_t2(1)
    deallocate(pe_t1,pe_t2)

    if (seq_comms(ID)%iamroot) then
       if (loglevel > 1) then
          write(logunit,F12) trim(subname),' initialize ID ',ID,seq_comms(ID)%name, &
          ' join IDs =',ID1,ID2,' npes =',seq_comms(ID)%npes, &
          ' nthreads =',seq_comms(ID)%nthreads, &
          ' cpl/cmp pes =',seq_comms(ID)%cplpe,seq_comms(ID)%cmppe
       else
          write(logunit,F13) trim(subname),' initialize ID ',ID,seq_comms(ID)%name, &
          ' join IDs =',ID1,ID2,' npes =',seq_comms(ID)%npes, &
          ' nthreads =',seq_comms(ID)%nthreads
       endif
    endif

  end subroutine seq_comm_joincomm

!---------------------------------------------------------
  subroutine seq_comm_jcommarr(IDs,ID,iname,inst,tinst)

    implicit none
    integer,intent(IN) :: IDs(:) ! src id
    integer,intent(IN) :: ID     ! computed id
    character(len=*),intent(IN),optional :: iname  ! comm name
    integer,intent(IN),optional :: inst
    integer,intent(IN),optional :: tinst

    integer :: mpigrp, mpigrpp
    integer :: mpicom, nids
    integer :: ierr
    integer :: n,nsize
    character(len=seq_comm_namelen) :: cname
    logical :: set_suffix
    integer,allocatable :: pe_t1(:),pe_t2(:)
    character(*),parameter :: subName =   '(seq_comm_jcommarr) '

    ! check that IDs are in valid range, that IDs have
    ! been set, and that ID has not been set

    nids = size(IDs)
    do n = 1,nids
       if (IDs(n) < 1 .or. IDs(n) > ncomps) then
          write(logunit,*) subname,' IDs out of range, abort ',n,IDs(n)
          call shr_sys_abort()
       endif 
       if (.not. seq_comms(IDs(n))%set) then
          write(logunit,*) subname,' IDs not set ',n,IDs(n)
          call shr_sys_abort()
       endif
    enddo

    if (ID < 1 .or. ID > ncomps) then
       write(logunit,*) subname,' ID out of range, abort ',ID
       call shr_sys_abort()
    endif
    if (seq_comms(ID)%set) then
       write(logunit,*) subname,' ID already set ',ID
       call shr_sys_abort()
    endif

    mpigrp = seq_comms(IDs(1))%mpigrp
    do n = 1,nids
       mpigrpp = mpigrp
       call mpi_group_union(mpigrpp,seq_comms(IDs(n))%mpigrp,mpigrp,ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_comm_union mpigrp')
    enddo
    call mpi_comm_create(GLOBAL_COMM, mpigrp, mpicom, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_create mpigrp')

    seq_comms(ID)%set = .true.
    seq_comms(ID)%ID = ID

    if (present(inst)) then
       seq_comms(ID)%inst = inst
    else
       seq_comms(ID)%inst = 1
    endif

    set_suffix = .true.
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
    
    seq_comms(ID)%nthreads = 1
    do n = 1,nids
       seq_comms(ID)%nthreads = max(seq_comms(ID)%nthreads,seq_comms(IDs(n))%nthreads)
    enddo

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
       seq_comms(ID)%iamroot = .false.
    endif

    seq_comms(ID)%cplpe = -1
    seq_comms(ID)%cmppe = -1

    if (seq_comms(ID)%iamroot) then
       if (loglevel > 1) then
          write(logunit,F14) trim(subname),' initialize ID ',ID,seq_comms(ID)%name, &
          ' join multiple comp IDs',' npes =',seq_comms(ID)%npes, &
          ' nthreads =',seq_comms(ID)%nthreads
       else
          write(logunit,F14) trim(subname),' initialize ID ',ID,seq_comms(ID)%name, &
          ' join multiple comp IDs',' npes =',seq_comms(ID)%npes, &
          ' nthreads =',seq_comms(ID)%nthreads
       endif
    endif

  end subroutine seq_comm_jcommarr

!---------------------------------------------------------
  subroutine seq_comm_printcomms()

    implicit none
    character(*),parameter :: subName =   '(seq_comm_printcomms) '
    integer :: m,n,mype,npes,ierr
    character(len=256) :: iamstring
    character(*),parameter :: F01 = "(4x,a4,4x   ,40(1x,a8))"
    character(*),parameter :: F02 = "(4x,i4,3x,a1,40(2x,i6,1x))"
    character(*),parameter :: F03 = "(4x,i4,3x,a1,a)"

    call mpi_comm_size(GLOBAL_COMM, npes  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size comm_world')
    call mpi_comm_rank(GLOBAL_COMM, mype  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank comm_world')

    call shr_sys_flush(logunit)
    call mpi_barrier(GLOBAL_COMM,ierr)
    if (mype == 0) then
       do n = 1,ncomps
          write(logunit,'(a,4i6,2x,3a)') trim(subName),n, &
             seq_comms(n)%gloroot,seq_comms(n)%npes,seq_comms(n)%nthreads, &
             trim(seq_comms(n)%name),':',trim(seq_comms(n)%suffix)
       enddo
!       write(logunit,*) ' '
!       write(logunit,*) trim(subName),' ID layout : global pes vs local pe for each ID'
!       write(logunit,F01) ' gpe',(seq_comms(n)%name,n=1,ncomps),'nthrds'
!       write(logunit,F01) ' ---',(' ------ '       ,n=1,ncomps),'------'
       call shr_sys_flush(logunit)
    endif
!    iamstring = ' '
!   do n = 1,ncomps
!      if (seq_comms(n)%iam >= 0) then
!         write(iamstring((n-1)*9+1:n*9),"(2x,i6,1x)") seq_comms(n)%iam
!      endif
!   enddo
!   n = ncomps + 1
!   write(iamstring((n-1)*9+1:n*9),"(2x,i6,1x)") seq_comms(GLOID)%pethreads

!    call shr_sys_flush(logunit)
!    call mpi_barrier(GLOBAL_COMM,ierr)
!   do m = 0,npes-1
!      if (mype == m) then
!!          write(logunit,F02) mype,':',(seq_comms(n)%iam,n=1,ncomps)
!         write(logunit,F03) mype,':',trim(iamstring)
!         if (m == npes-1) then
!            write(logunit,*) ' '
!         endif
!      endif
!      call shr_sys_flush(logunit)
!      call mpi_barrier(GLOBAL_COMM,ierr)
!   enddo

  end subroutine seq_comm_printcomms

!---------------------------------------------------------
  subroutine seq_comm_setptrs(ID,mpicom,mpigrp,npes,nthreads,iam,iamroot,gloiam,gloroot, &
                                 cplpe,cmppe,pethreads, name)

    implicit none
    integer,intent(in) :: ID
    integer,intent(out),optional :: mpicom
    integer,intent(out),optional :: mpigrp
    integer,intent(out),optional :: npes
    integer,intent(out),optional :: nthreads
    integer,intent(out),optional :: iam
    logical,intent(out),optional :: iamroot
    integer,intent(out),optional :: gloiam
    integer,intent(out),optional :: gloroot
    integer,intent(out),optional :: cplpe
    integer,intent(out),optional :: cmppe
    integer,intent(out),optional :: pethreads
    character(len=seq_comm_namelen)  , intent(out), optional :: name
    character(*),parameter :: subName =   '(seq_comm_setptrs) '

    if (ID < 1 .or. ID > ncomps) then
       write(logunit,*) subname,' ID out of range, return ',ID
       return
    endif 

    if (present(mpicom)) then
       mpicom = seq_comms(ID)%mpicom
    endif

    if (present(mpigrp)) then
       mpigrp = seq_comms(ID)%mpigrp
    endif

    if (present(npes)) then
       npes = seq_comms(ID)%npes
    endif

    if (present(nthreads)) then
       nthreads = seq_comms(ID)%nthreads
    endif

    if (present(iam)) then
       iam = seq_comms(ID)%iam
    endif

    if (present(iamroot)) then
       iamroot = seq_comms(ID)%iamroot
    endif

    if (present(gloiam)) then
       gloiam = seq_comms(ID)%gloiam
    endif

    if (present(gloroot)) then
       gloroot = seq_comms(ID)%gloroot
    endif

    if (present(cplpe)) then
       cplpe = seq_comms(ID)%cplpe
    endif

    if (present(cmppe)) then
       cmppe = seq_comms(ID)%cmppe
    endif

    if (present(pethreads)) then
       pethreads = seq_comms(ID)%pethreads
    endif

    if(present(name)) then
       name = seq_comms(ID)%name
    end if



  end subroutine seq_comm_setptrs
!---------------------------------------------------------
  subroutine seq_comm_setnthreads(nthreads)

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

    if (seq_comms(ID)%iam >= 0) then
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

    seq_comm_iamroot = seq_comms(ID)%iamroot

  end function seq_comm_iamroot
!---------------------------------------------------------
  integer function seq_comm_mpicom(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_mpicom) '

    seq_comm_mpicom = seq_comms(ID)%mpicom

  end function seq_comm_mpicom
!---------------------------------------------------------
  integer function seq_comm_iam(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_iam) '

    seq_comm_iam = seq_comms(ID)%iam

  end function seq_comm_iam
!---------------------------------------------------------
  integer function seq_comm_gloiam(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_gloiam) '

    seq_comm_gloiam = seq_comms(ID)%gloiam

  end function seq_comm_gloiam
!---------------------------------------------------------
  integer function seq_comm_gloroot(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_gloroot) '

    seq_comm_gloroot = seq_comms(ID)%gloroot

  end function seq_comm_gloroot
!---------------------------------------------------------
  integer function seq_comm_cplpe(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_cplpe) '

    seq_comm_cplpe = seq_comms(ID)%cplpe

  end function seq_comm_cplpe
!---------------------------------------------------------
  integer function seq_comm_cmppe(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_cmppe) '

    seq_comm_cmppe = seq_comms(ID)%cmppe

  end function seq_comm_cmppe
!---------------------------------------------------------
  character(len=seq_comm_namelen) function seq_comm_name(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_name) '

    seq_comm_name = trim(seq_comms(ID)%name)

  end function seq_comm_name
!---------------------------------------------------------
  character(len=seq_comm_namelen) function seq_comm_suffix(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_suffix) '

    seq_comm_suffix = trim(seq_comms(ID)%suffix)

  end function seq_comm_suffix
!---------------------------------------------------------
  subroutine seq_comm_petlist(ID,petlist)

    implicit none
    integer,intent(in) :: ID
    integer,pointer :: petlist(:)
    character(*),parameter :: subName =   '(seq_comm_petlist) '

    petlist => seq_comms(ID)%petlist

  end subroutine seq_comm_petlist
!---------------------------------------------------------
#ifdef USE_ESMF_LIB
  subroutine seq_comm_getcompstates(ID,comp,imp_state,exp_state)

    implicit none
    integer,                       intent(in)  :: ID
    type(ESMF_GridComp), optional, intent(out) :: comp
    type(ESMF_State),    optional, intent(out) :: imp_state, exp_state
    character(*),parameter :: subName =   '(seq_comm_getcompstates) '

    if(present(comp))      comp      = seq_comms(ID)%esmf_comp
    if(present(imp_state)) imp_state = seq_comms(ID)%imp_state
    if(present(imp_state)) exp_state = seq_comms(ID)%exp_state

  end subroutine seq_comm_getcompstates
!---------------------------------------------------------
  subroutine seq_comm_setcompstates(ID,comp,imp_state,exp_state)

    implicit none
    integer,                       intent(in) :: ID
    type(ESMF_GridComp), optional, intent(in) :: comp
    type(ESMF_State)   , optional, intent(in) :: imp_state, exp_state
    character(*),parameter :: subName =   '(seq_comm_setcompstates) '

    if(present(comp))      seq_comms(ID)%esmf_comp = comp
    if(present(imp_state)) seq_comms(ID)%imp_state = imp_state
    if(present(imp_state)) seq_comms(ID)%exp_state = exp_state

  end subroutine seq_comm_setcompstates
#endif
!---------------------------------------------------------
  integer function seq_comm_inst(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_inst) '

    seq_comm_inst = seq_comms(ID)%inst

  end function seq_comm_inst
!---------------------------------------------------------
  subroutine seq_comm_mkname(oname,str1,num)
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
