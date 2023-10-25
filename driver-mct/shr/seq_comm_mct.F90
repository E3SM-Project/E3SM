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

  use mct_mod        , only : mct_world_init, mct_world_clean, mct_die
  use shr_sys_mod    , only : shr_sys_abort, shr_sys_flush
  use shr_mpi_mod    , only : shr_mpi_chkerr, shr_mpi_bcast, shr_mpi_max
  use shr_file_mod   , only : shr_file_getUnit, shr_file_freeUnit
  ! gptl timing library is not built for unit tests but it is on
  ! by default for Makefile (model) builds.
#ifdef TIMING
  use shr_taskmap_mod, only : shr_taskmap_write
  use perf_mod       , only : t_startf, t_stopf
#endif
  use esmf           , only : ESMF_LogKind_Flag, ESMF_LOGKIND_NONE
  use esmf           , only : ESMF_LOGKIND_SINGLE, ESMF_LOGKIND_MULTI

  implicit none

  private
#include <mpif.h>

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
  public seq_comm_setptrs
  public seq_comm_setnthreads
  public seq_comm_getnthreads
  public seq_comm_printcomms
  public seq_comm_get_ncomps

  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------

  integer, public, parameter :: default_logunit = 6
  integer, public :: logunit  = default_logunit     ! log unit number
  integer, public :: loglevel = 1     ! log level

  integer, public :: global_mype = -1  !! To be initialized

 !!! Note - NUM_COMP_INST_XXX are cpp variables set in buildlib.csm_share

  integer, parameter :: ncomptypes = 9  ! total number of component types
  integer, parameter :: ncouplers  = 1  ! number of couplers
  integer, parameter, public :: num_inst_atm = NUM_COMP_INST_ATM
  integer, parameter, public :: num_inst_lnd = NUM_COMP_INST_LND
  integer, parameter, public :: num_inst_ocn = NUM_COMP_INST_OCN
  integer, parameter, public :: num_inst_ice = NUM_COMP_INST_ICE
  integer, parameter, public :: num_inst_glc = NUM_COMP_INST_GLC
  integer, parameter, public :: num_inst_wav = NUM_COMP_INST_WAV
  integer, parameter, public :: num_inst_rof = NUM_COMP_INST_ROF
  integer, parameter, public :: num_inst_iac = NUM_COMP_INST_IAC
  integer, parameter, public :: num_inst_esp = NUM_COMP_INST_ESP

  integer, parameter, public :: num_inst_total= num_inst_atm + &
       num_inst_lnd + &
       num_inst_ocn + &
       num_inst_ice + &
       num_inst_glc + &
       num_inst_wav + &
       num_inst_rof + &
       num_inst_iac + &
       num_inst_esp + 1

  integer, public :: num_inst_min, num_inst_max
  integer, public :: num_inst_xao    ! for xao flux
  integer, public :: num_inst_frc    ! for fractions
  integer, public :: num_inst_driver = 1

!!! Each component instance needs two communicators: one internal to the
!!! instance, and one for communicating with the coupler.
!!! Additionally, one communicator is needed for the coupler's
!!! internal communications, and one is needed for the global space.
!!! All instances of a component type also share a separate communicator
!!! All instances of a component type share a communicator with the coupler

  integer, parameter, public :: num_inst_phys = num_inst_atm + num_inst_lnd + &
       num_inst_ocn + num_inst_ice + &
       num_inst_glc + num_inst_rof + &
       num_inst_wav + num_inst_esp + &
       num_inst_iac
  integer, parameter, public :: num_cpl_phys  = num_inst_atm + num_inst_lnd + &
       num_inst_ocn + num_inst_ice + &
       num_inst_glc + num_inst_rof + &
       num_inst_wav + num_inst_esp + &
       num_inst_iac
  integer, parameter :: ncomps = (1 + ncouplers + 2*ncomptypes + num_inst_phys + num_cpl_phys)

  integer, public :: GLOID
  integer, public :: CPLID

  integer, public :: ALLATMID
  integer, public :: ALLLNDID
  integer, public :: ALLOCNID
  integer, public :: ALLICEID
  integer, public :: ALLGLCID
  integer, public :: ALLROFID
  integer, public :: ALLWAVID
  integer, public :: ALLIACID
  integer, public :: ALLESPID

  integer, public :: CPLALLATMID
  integer, public :: CPLALLLNDID
  integer, public :: CPLALLOCNID
  integer, public :: CPLALLICEID
  integer, public :: CPLALLGLCID
  integer, public :: CPLALLROFID
  integer, public :: CPLALLWAVID
  integer, public :: CPLALLIACID
  integer, public :: CPLALLESPID

  integer, public :: ATMID(num_inst_atm)
  integer, public :: LNDID(num_inst_lnd)
  integer, public :: OCNID(num_inst_ocn)
  integer, public :: ICEID(num_inst_ice)
  integer, public :: GLCID(num_inst_glc)
  integer, public :: ROFID(num_inst_rof)
  integer, public :: WAVID(num_inst_wav)
  integer, public :: IACID(num_inst_iac)
  integer, public :: ESPID(num_inst_esp)

  integer, public :: CPLATMID(num_inst_atm)
  integer, public :: CPLLNDID(num_inst_lnd)
  integer, public :: CPLOCNID(num_inst_ocn)
  integer, public :: CPLICEID(num_inst_ice)
  integer, public :: CPLGLCID(num_inst_glc)
  integer, public :: CPLROFID(num_inst_rof)
  integer, public :: CPLWAVID(num_inst_wav)
  integer, public :: CPLIACID(num_inst_iac)
  integer, public :: CPLESPID(num_inst_esp)

  integer, parameter, public :: seq_comm_namelen=16

  ! taskmap output level specifications for components
  ! (0:no output, 1:compact, 2:verbose)
  integer, public :: info_taskmap_model, info_taskmap_comp
  integer, public :: driver_nnodes
  integer, public, allocatable :: driver_task_node_map(:)
  integer, public :: info_mprof, info_mprof_dt

  ! suffix for log and timing files if multi coupler driver
  character(len=seq_comm_namelen), public  :: cpl_inst_tag

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

     integer :: gloiam          ! my task number in global_comm
     integer :: gloroot         ! the global task number of each comps root on all pes

     integer :: pethreads       ! max number of threads on my task
     integer :: cplpe           ! a common task in mpicom from the cpl group for join mpicoms
     ! cplpe is used to broadcast information from the coupler to the component
     integer :: cmppe           ! a common task in mpicom from the component group for join mpicoms
     ! cmppe is used to broadcast information from the component to the coupler
     logical :: set             ! has this datatype been set
     integer :: excl_group      ! mpi group of tasks owned exclusively by this component

  end type seq_comm_type

  type(seq_comm_type) :: seq_comms(ncomps)

  character(*), parameter :: layout_concurrent = 'concurrent'
  character(*), parameter :: layout_sequential = 'sequential'

  character(*), parameter :: F11 = "(a,a,'(',i3,' ',a,')',a,   3i6,' (',a,i6,')',' (',a,i3,')','(',a,a,')')"
  character(*), parameter :: F12 = "(a,a,'(',i3,' ',a,')',a,2i6,6x,' (',a,i6,')',' (',a,i3,')','(',a,2i6,')')"
  character(*), parameter :: F13 = "(a,a,'(',i3,' ',a,')',a,2i6,6x,' (',a,i6,')',' (',a,i3,')')"
  character(*), parameter :: F14 = "(a,a,'(',i3,' ',a,')',a,    6x,' (',a,i6,')',' (',a,i3,')')"

  ! Exposed for use in the esp component, please don't use this elsewhere
  integer, public :: Global_Comm
  integer :: driver_comm

  character(len=32), public :: &
       atm_layout, lnd_layout, ice_layout, glc_layout, rof_layout, &
       ocn_layout, wav_layout, esp_layout, iac_layout

  logical :: seq_comm_mct_initialized = .false.  ! whether this module has been initialized

  !=======================================================================
contains
  !======================================================================
  integer function seq_comm_get_ncomps()
    seq_comm_get_ncomps = ncomps
  end function seq_comm_get_ncomps

  subroutine seq_comm_init(global_comm_in, driver_comm_in, nmlfile, drv_comm_id)
    !----------------------------------------------------------
    !
    ! Arguments
    implicit none
    integer, intent(in) :: global_comm_in
    integer, intent(in) :: driver_comm_in
    character(len=*), intent(IN) :: nmlfile
    integer, intent(in), optional :: drv_comm_id
    !
    ! Local variables
    !
    logical :: error_state
    integer :: ierr, n, count, xcount
    character(*), parameter :: subName =   '(seq_comm_init) '
    integer :: mype,numpes,myncomps,max_threads,gloroot, global_numpes
    integer :: pelist(3,1)       ! start, stop, stride for group
    integer :: exlist(3), mpigrp_world, exgrp
    integer, pointer :: comps(:) ! array with component ids
    integer, pointer :: comms(:) ! array with mpicoms
    integer :: nu
    logical :: verbose_taskmap_output
    integer :: drv_inst
    character(len=8) :: c_drv_inst      ! driver instance number
    character(len=8) :: c_driver_numpes ! number of pes in driver
    character(len=16):: c_comm_name     ! comm. name
    character(len=seq_comm_namelen) :: valid_comps(ncomps)

    integer :: &
         atm_ntasks, atm_rootpe, atm_pestride, atm_excl_stride, atm_nthreads, &
         lnd_ntasks, lnd_rootpe, lnd_pestride, lnd_excl_stride, lnd_nthreads, &
         ice_ntasks, ice_rootpe, ice_pestride, ice_excl_stride, ice_nthreads, &
         glc_ntasks, glc_rootpe, glc_pestride, glc_excl_stride, glc_nthreads, &
         wav_ntasks, wav_rootpe, wav_pestride, wav_excl_stride, wav_nthreads, &
         rof_ntasks, rof_rootpe, rof_pestride, rof_excl_stride, rof_nthreads, &
         ocn_ntasks, ocn_rootpe, ocn_pestride, ocn_excl_stride, ocn_nthreads, &
         esp_ntasks, esp_rootpe, esp_pestride, esp_excl_stride, esp_nthreads, &
         iac_ntasks, iac_rootpe, iac_pestride, iac_excl_stride, iac_nthreads, &
         cpl_ntasks, cpl_rootpe, cpl_pestride, cpl_excl_stride, cpl_nthreads

    namelist /cime_pes/  &
         atm_ntasks, atm_rootpe, atm_pestride, atm_excl_stride, atm_nthreads, atm_layout, &
         lnd_ntasks, lnd_rootpe, lnd_pestride, lnd_excl_stride, lnd_nthreads, lnd_layout, &
         ice_ntasks, ice_rootpe, ice_pestride, ice_excl_stride, ice_nthreads, ice_layout, &
         glc_ntasks, glc_rootpe, glc_pestride, glc_excl_stride, glc_nthreads, glc_layout, &
         wav_ntasks, wav_rootpe, wav_pestride, wav_excl_stride, wav_nthreads, wav_layout, &
         rof_ntasks, rof_rootpe, rof_pestride, rof_excl_stride, rof_nthreads, rof_layout, &
         ocn_ntasks, ocn_rootpe, ocn_pestride, ocn_excl_stride, ocn_nthreads, ocn_layout, &
         esp_ntasks, esp_rootpe, esp_pestride, esp_excl_stride, esp_nthreads, esp_layout, &
         iac_ntasks, iac_rootpe, iac_pestride, iac_excl_stride, iac_nthreads, iac_layout, &
         cpl_ntasks, cpl_rootpe, cpl_pestride, cpl_excl_stride, cpl_nthreads,             &
         info_taskmap_model, info_taskmap_comp, info_mprof, info_mprof_dt
    !----------------------------------------------------------

    ! make sure this is first pass and set comms unset
    if (seq_comm_mct_initialized) then
       write(logunit,*) trim(subname),' ERROR seq_comm_init already called '
       call shr_sys_abort()
    endif
    seq_comm_mct_initialized = .true.
    global_comm = global_comm_in
    driver_comm = driver_comm_in

    !! Initialize seq_comms elements

    do n = 1,ncomps
       seq_comms(n)%name = 'unknown'
       seq_comms(n)%suffix = ' '
       seq_comms(n)%inst = 0
       seq_comms(n)%set = .false.
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
       seq_comms(n)%excl_group = -1
    enddo


    ! Initialize MPI
    ! Note that if no MPI, will call MCTs fake version

    call mpi_comm_size(GLOBAL_COMM_IN, global_numpes , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size comm_world')
    call mpi_comm_rank(DRIVER_COMM, mype  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank driver')
    call mpi_comm_size(DRIVER_COMM, numpes, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size driver')

    if (mod(global_numpes, numpes) .ne. 0) then
       write(logunit,*) trim(subname),' ERROR: numpes driver: ', numpes, ' should divide global_numpes: ',global_numpes
       call shr_sys_abort(trim(subname)//' ERROR decomposition error ')
    endif

    ! Initialize gloiam on all IDs

    global_mype = mype

    do n = 1,ncomps
       seq_comms(n)%gloiam = mype
    enddo

    ! Set ntasks, rootpe, pestride, nthreads for all components

    if (mype == 0) then
       !! Set up default component process parameters
       call comp_pelayout_init(numpes, atm_ntasks, atm_rootpe, atm_pestride, atm_nthreads, atm_layout)
       call comp_pelayout_init(numpes, lnd_ntasks, lnd_rootpe, lnd_pestride, lnd_nthreads, lnd_layout)
       call comp_pelayout_init(numpes, ice_ntasks, ice_rootpe, ice_pestride, ice_nthreads, ice_layout)
       call comp_pelayout_init(numpes, ocn_ntasks, ocn_rootpe, ocn_pestride, ocn_nthreads, ocn_layout)
       call comp_pelayout_init(numpes, rof_ntasks, rof_rootpe, rof_pestride, rof_nthreads, rof_layout)
       call comp_pelayout_init(numpes, wav_ntasks, wav_rootpe, wav_pestride, wav_nthreads, wav_layout)
       call comp_pelayout_init(numpes, glc_ntasks, glc_rootpe, glc_pestride, glc_nthreads, glc_layout)
       call comp_pelayout_init(numpes, esp_ntasks, esp_rootpe, esp_pestride, esp_nthreads, esp_layout)
       call comp_pelayout_init(numpes, iac_ntasks, iac_rootpe, iac_pestride, iac_nthreads, iac_layout)
       call comp_pelayout_init(numpes, cpl_ntasks, cpl_rootpe, cpl_pestride, cpl_nthreads)
       info_taskmap_model = 0
       info_taskmap_comp  = 0
       info_mprof         = 0
       info_mprof_dt      = 86400

       ! Read namelist if it exists

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

    call shr_mpi_bcast(atm_nthreads,DRIVER_COMM,'atm_nthreads')
    call shr_mpi_bcast(lnd_nthreads,DRIVER_COMM,'lnd_nthreads')
    call shr_mpi_bcast(ocn_nthreads,DRIVER_COMM,'ocn_nthreads')
    call shr_mpi_bcast(ice_nthreads,DRIVER_COMM,'ice_nthreads')
    call shr_mpi_bcast(glc_nthreads,DRIVER_COMM,'glc_nthreads')
    call shr_mpi_bcast(wav_nthreads,DRIVER_COMM,'wav_nthreads')
    call shr_mpi_bcast(rof_nthreads,DRIVER_COMM,'rof_nthreads')
    call shr_mpi_bcast(esp_nthreads,DRIVER_COMM,'esp_nthreads')
    call shr_mpi_bcast(iac_nthreads,DRIVER_COMM,'iac_nthreads')
    call shr_mpi_bcast(cpl_nthreads,DRIVER_COMM,'cpl_nthreads')

    call shr_mpi_bcast(atm_layout,DRIVER_COMM,'atm_layout')
    call shr_mpi_bcast(lnd_layout,DRIVER_COMM,'lnd_layout')
    call shr_mpi_bcast(ocn_layout,DRIVER_COMM,'ocn_layout')
    call shr_mpi_bcast(ice_layout,DRIVER_COMM,'ice_layout')
    call shr_mpi_bcast(glc_layout,DRIVER_COMM,'glc_layout')
    call shr_mpi_bcast(wav_layout,DRIVER_COMM,'wav_layout')
    call shr_mpi_bcast(rof_layout,DRIVER_COMM,'rof_layout')
    call shr_mpi_bcast(iac_layout,DRIVER_COMM,'iac_layout')
    call shr_mpi_bcast(esp_layout,DRIVER_COMM,'esp_layout')

    call shr_mpi_bcast(info_taskmap_model,DRIVER_COMM,'info_taskmap_model')
    call shr_mpi_bcast(info_taskmap_comp, DRIVER_COMM,'info_taskmap_comp' )
    call shr_mpi_bcast(info_mprof,   DRIVER_COMM,'info_mprof')
    call shr_mpi_bcast(info_mprof_dt,DRIVER_COMM,'info_mprof_dt')

#ifdef TIMING
    if (info_taskmap_model > 0) then
       ! output task-to-node mapping

       if (info_taskmap_model == 1) then
         verbose_taskmap_output = .false.
       else
         verbose_taskmap_output = .true.
       endif

       if (present(drv_comm_id)) then
          drv_inst = drv_comm_id
          write(c_drv_inst,'(i8)') drv_inst
       else
          drv_inst = 0
       endif

       if (mype == 0) then
          write(c_driver_numpes,'(i8)') numpes
          if (drv_inst == 0) then
             write(logunit,'(2A)') trim(adjustl(c_driver_numpes)), &
                                   ' pes participating in computation of coupled model'
          else
             write(logunit,'(3A)') trim(adjustl(c_driver_numpes)), &
                                   ' pes participating in computation of DRIVER instance #', &
                                   trim(adjustl(c_drv_inst))
          endif
          call shr_sys_flush(logunit)
       endif

       if (info_mprof > 2) then
          allocate( driver_task_node_map(0:global_numpes-1), stat=ierr)
          if (ierr /= 0) call shr_sys_abort(trim(subname)//' allocate driver_task_node_map failed ')
       endif

       call t_startf("shr_taskmap_write")
       if (drv_inst == 0) then
          c_comm_name = 'GLOBAL'
       else
          c_comm_name = 'DRIVER #'//trim(adjustl(c_drv_inst))
       endif
       if (info_mprof > 2) then
          call shr_taskmap_write(logunit, DRIVER_COMM, &
                                 c_comm_name, &
                                 verbose=verbose_taskmap_output, &
                                 save_nnodes=driver_nnodes, &
                                 save_task_node_map=driver_task_node_map)
       else
          call shr_taskmap_write(logunit, DRIVER_COMM, &
                                 c_comm_name, &
                                 verbose=verbose_taskmap_output)
       endif
       call t_stopf("shr_taskmap_write")
    endif
#endif

    !--- compute some other num_inst values

    num_inst_xao = max(num_inst_atm,num_inst_ocn)
    num_inst_frc = num_inst_ice

    !--- compute num_inst_min, num_inst_max
    !--- instances must be either 1 or a constant across components
    !--- checks for prognostic/present consistency in the driver

    error_state = .false.
    num_inst_min = min(num_inst_atm, num_inst_lnd, num_inst_ocn,&
         num_inst_ice, num_inst_glc, num_inst_wav, num_inst_rof,&
         num_inst_esp, num_inst_iac)
    num_inst_max = max(num_inst_atm, num_inst_lnd, num_inst_ocn,&
         num_inst_ice, num_inst_glc, num_inst_wav, num_inst_rof,&
         num_inst_esp, num_inst_iac)

    if (num_inst_min /= num_inst_max .and. num_inst_min /= 1) error_state = .true.
    if (num_inst_atm /= num_inst_min .and. num_inst_atm /= num_inst_max) error_state = .true.
    if (num_inst_lnd /= num_inst_min .and. num_inst_lnd /= num_inst_max) error_state = .true.
    if (num_inst_ocn /= num_inst_min .and. num_inst_ocn /= num_inst_max) error_state = .true.
    if (num_inst_ice /= num_inst_min .and. num_inst_ice /= num_inst_max) error_state = .true.
    if (num_inst_glc /= num_inst_min .and. num_inst_glc /= num_inst_max) error_state = .true.
    if (num_inst_wav /= num_inst_min .and. num_inst_wav /= num_inst_max) error_state = .true.
    if (num_inst_rof /= num_inst_min .and. num_inst_rof /= num_inst_max) error_state = .true.
    if (num_inst_iac /= num_inst_min .and. num_inst_iac /= num_inst_max) error_state = .true.
    if (num_inst_esp /= num_inst_min .and. num_inst_esp /= num_inst_max) error_state = .true.

    if (error_state) then
       write(logunit,*) trim(subname),' ERROR: num_inst inconsistent'
       write(logunit,*) num_inst_atm, num_inst_lnd, num_inst_ocn,&
            num_inst_ice, num_inst_glc, num_inst_wav, num_inst_rof,&
            num_inst_esp, num_inst_min, num_inst_max
       call shr_sys_abort(trim(subname)//' ERROR: num_inst inconsistent')
    endif

    ! Initialize IDs

    count = 0

    count = count + 1
    GLOID = count
    count = count + 1
    CPLID = count

    if (mype == 0) then
       pelist(1,1) = 0
       pelist(2,1) = numpes-1
       pelist(3,1) = 1
    end if
    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, DRIVER_COMM, ierr)
    call seq_comm_setcomm(GLOID, pelist,iname='GLOBAL')

    if (mype == 0) then
       pelist(1,1) = cpl_rootpe
       pelist(2,1) = cpl_rootpe + (cpl_ntasks -1) * cpl_pestride
       pelist(3,1) = cpl_pestride
       exlist(1) = pelist(1,1)
       exlist(2) = pelist(2,1)
       exlist(3) = cpl_excl_stride
    end if

    call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, DRIVER_COMM, ierr)
    call mpi_bcast(exlist, size(exlist), MPI_INTEGER, 0, DRIVER_COMM, ierr)
    call mpi_comm_group(DRIVER_COMM, mpigrp_world, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_group mpigrp_world')
    if (exlist(3) > 0) then
       call mpi_group_range_incl(mpigrp_world, 1, exlist, exgrp, ierr)
       call shr_mpi_chkerr(ierr,subname//' mpi_group_range_incl CPLID')
       seq_comms(CPLID)%excl_group = exgrp
    endif
    call seq_comm_setcomm(CPLID,pelist,nthreads=cpl_nthreads,iname='CPL')

    ! init excl-strides
    xcount = CPLID + 2*num_inst_atm + 1
    call comp_exstride_init(driver_comm, lnd_rootpe, lnd_ntasks, lnd_pestride, &
         lnd_excl_stride, num_inst_lnd, xcount, mpigrp_world)
    call comp_exstride_init(driver_comm, ice_rootpe, ice_ntasks, ice_pestride, &
         ice_excl_stride, num_inst_ice, xcount, mpigrp_world)
    call comp_exstride_init(driver_comm, ocn_rootpe, ocn_ntasks, ocn_pestride, &
         ocn_excl_stride, num_inst_ocn, xcount, mpigrp_world)
    call comp_exstride_init(driver_comm, rof_rootpe, rof_ntasks, rof_pestride, &
         rof_excl_stride, num_inst_rof, xcount, mpigrp_world)
    call comp_exstride_init(driver_comm, glc_rootpe, glc_ntasks, glc_pestride, &
         glc_excl_stride, num_inst_glc, xcount, mpigrp_world)
    call comp_exstride_init(driver_comm, wav_rootpe, wav_ntasks, wav_pestride, &
         wav_excl_stride, num_inst_wav, xcount, mpigrp_world)
    call comp_exstride_init(driver_comm, esp_rootpe, esp_ntasks, esp_pestride, &
         esp_excl_stride, num_inst_esp, xcount, mpigrp_world)
    call comp_exstride_init(driver_comm, iac_rootpe, iac_ntasks, iac_pestride, &
         iac_excl_stride, num_inst_iac, xcount, mpigrp_world)

    call comp_comm_init(driver_comm, atm_rootpe, atm_nthreads, atm_layout, &
         atm_ntasks, atm_pestride, atm_excl_stride, num_inst_atm, &
         CPLID, ATMID, CPLATMID, ALLATMID, CPLALLATMID, 'ATM', count, drv_comm_id)
    call comp_comm_init(driver_comm, lnd_rootpe, lnd_nthreads, lnd_layout, &
         lnd_ntasks, lnd_pestride, lnd_excl_stride, num_inst_lnd, &
         CPLID, LNDID, CPLLNDID, ALLLNDID, CPLALLLNDID, 'LND', count, drv_comm_id)
    call comp_comm_init(driver_comm, ice_rootpe, ice_nthreads, ice_layout, &
         ice_ntasks, ice_pestride, ice_excl_stride, num_inst_ice, &
         CPLID, ICEID, CPLICEID, ALLICEID, CPLALLICEID, 'ICE', count, drv_comm_id)
    call comp_comm_init(driver_comm, ocn_rootpe, ocn_nthreads, ocn_layout, &
         ocn_ntasks, ocn_pestride, ocn_excl_stride, num_inst_ocn, &
         CPLID, OCNID, CPLOCNID, ALLOCNID, CPLALLOCNID, 'OCN', count, drv_comm_id)
    call comp_comm_init(driver_comm, rof_rootpe, rof_nthreads, rof_layout, &
         rof_ntasks, rof_pestride, rof_excl_stride, num_inst_rof, &
         CPLID, ROFID, CPLROFID, ALLROFID, CPLALLROFID, 'ROF', count, drv_comm_id)
    call comp_comm_init(driver_comm, glc_rootpe, glc_nthreads, glc_layout, &
         glc_ntasks, glc_pestride, glc_excl_stride, num_inst_glc, &
         CPLID, GLCID, CPLGLCID, ALLGLCID, CPLALLGLCID, 'GLC', count, drv_comm_id)
    call comp_comm_init(driver_comm, wav_rootpe, wav_nthreads, wav_layout, &
         wav_ntasks, wav_pestride, wav_excl_stride, num_inst_wav, &
         CPLID, WAVID, CPLWAVID, ALLWAVID, CPLALLWAVID, 'WAV', count, drv_comm_id)
    call comp_comm_init(driver_comm, esp_rootpe, esp_nthreads, esp_layout, &
         esp_ntasks, esp_pestride, esp_excl_stride, num_inst_esp, &
         CPLID, ESPID, CPLESPID, ALLESPID, CPLALLESPID, 'ESP', count, drv_comm_id)
    call comp_comm_init(driver_comm, iac_rootpe, iac_nthreads, iac_layout, &
         iac_ntasks, iac_pestride, iac_excl_stride, num_inst_iac, &
         CPLID, IACID, CPLIACID, ALLIACID, CPLALLIACID, 'IAC', count, drv_comm_id)

    if (count /= ncomps) then
       write(logunit,*) trim(subname),' ERROR in ID count ',count,ncomps
       call shr_sys_abort(trim(subname)//' ERROR in ID count')
    endif
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
       call shr_mpi_max(gloroot,seq_comms(n)%gloroot,DRIVER_COMM, &
            trim(subname)//' gloroot',all=.true.)
    enddo

    ! Initialize MCT

    ! ensure that all driver_comm processes initialized their comms
    call mpi_barrier(DRIVER_COMM,ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_barrier driver pre-mct-init')

    ! add up valid comps on local pe

    valid_comps = '*'
    myncomps = 0
    do n = 1,ncomps
       if (seq_comms(n)%mpicom /= MPI_COMM_NULL) then
          myncomps = myncomps + 1
          valid_comps(n) = seq_comms(n)%name
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
       write(logunit,*) trim(subname),' ERROR in myncomps ',myncomps,size(comps),comps,valid_comps
       call shr_sys_abort()
    endif

    call mct_world_init(ncomps, DRIVER_COMM, comms, comps)

    deallocate(comps,comms)


    call seq_comm_printcomms()

  end subroutine seq_comm_init

  subroutine comp_exstride_init(driver_comm, comp_rootpe, comp_ntasks, comp_pestride, &
       comp_exstride, num_inst_comp, xcount, drv_grp)
    integer, intent(in) :: driver_comm
    integer, intent(in) :: comp_rootpe
    integer, intent(in) :: comp_ntasks
    integer, intent(in) :: comp_pestride
    integer, intent(in) :: comp_exstride
    integer, intent(in) :: num_inst_comp
    integer, intent(inout) :: xcount
    integer, intent(in) :: drv_grp

    integer :: exlist(3), exgrp, ierr

    xcount = xcount + 2*num_inst_comp + 2

    exlist(1) = comp_rootpe
    exlist(2) = comp_rootpe + (comp_ntasks - 1) * comp_pestride
    exlist(3) = comp_exstride
    call mpi_bcast(exlist, 3, MPI_INTEGER, 0, driver_comm, ierr)

    if (exlist(3) > 0) then
       call mpi_group_range_incl(drv_grp, 1, exlist, exgrp, ierr)
       seq_comms(xcount)%excl_group = exgrp
    endif

  end subroutine comp_exstride_init

  subroutine comp_comm_init(driver_comm, comp_rootpe, comp_nthreads, comp_layout, &
       comp_ntasks, comp_pestride, comp_exstride, num_inst_comp, &
       CPLID, COMPID, CPLCOMPID, ALLCOMPID, CPLALLCOMPID, name, count, drv_comm_id)
    integer, intent(in) :: driver_comm
    integer, intent(in) :: comp_rootpe
    integer, intent(in) :: comp_nthreads
    character(len=*), intent(in) :: comp_layout
    integer, intent(in) :: comp_ntasks
    integer, intent(in) :: comp_pestride
    integer, intent(in) :: comp_exstride
    integer, intent(in) :: num_inst_comp
    integer, intent(in) :: CPLID
    integer, intent(out) :: COMPID(num_inst_comp)
    integer, intent(out) :: CPLCOMPID(num_inst_comp)
    integer, intent(out) :: ALLCOMPID
    integer, intent(out) :: CPLALLCOMPID
    integer, intent(inout) :: count
    integer, intent(in), optional :: drv_comm_id
    character(len=*), intent(in) :: name

    character(len=*), parameter :: subname = "(comp_comm_init) "
    integer :: comp_inst_tasks
    integer :: droot
    integer :: current_task_rootpe
    integer :: cmin(num_inst_comp), cmax(num_inst_comp), cstr(num_inst_comp)
    integer :: n
    integer :: pelist(3,1), exlist(3), mpigrp_world, exgrp
    integer :: ierr
    integer :: mype

    call mpi_comm_rank(driver_comm, mype, ierr)

    count = count + 1
    ALLCOMPID = count
    count = count + 1
    CPLALLCOMPID = count
    do n = 1, num_inst_comp
       count = count + 1
       COMPID(n) = count
       count = count + 1
       CPLCOMPID(n) = count
    enddo

    if (mype == 0) then
       !--- validation of inputs ---
       ! rootpes >= 0
       !! Determine the process layout
       !!
       !! We will assign comp_ntasks / num_inst_comp tasks to each component
       !! instance.  (This may lead to unallocated tasks if comp_ntasks is
       !! not an integer multiple of num_inst_comp.)

       if (comp_rootpe < 0) then
          call shr_sys_abort(trim(subname)//' ERROR: rootpes must be >= 0 for component '//trim(name))
       endif

       if (trim(comp_layout) == trim(layout_concurrent)) then
          comp_inst_tasks = comp_ntasks / num_inst_comp
          droot = (comp_inst_tasks * comp_pestride)
       elseif (trim(comp_layout) == trim(layout_sequential)) then
          comp_inst_tasks = comp_ntasks
          droot = 0
       else
          call shr_sys_abort(subname//' ERROR invalid comp_layout for component '//trim(name))
       endif
       current_task_rootpe = comp_rootpe
       do n = 1, num_inst_comp
          cmin(n) = current_task_rootpe
          cmax(n) = current_task_rootpe &
               + ((comp_inst_tasks - 1) * comp_pestride)
          cstr(n) = comp_pestride
          current_task_rootpe = current_task_rootpe + droot
       end do
    endif

    do n = 1, num_inst_comp
       if (mype==0) then
          pelist(1,1) = cmin(n)
          pelist(2,1) = cmax(n)
          pelist(3,1) = cstr(n)
          exlist(1)   = cmin(n)
          exlist(2)   = cmax(n)
          exlist(3)   = comp_exstride
       endif
       call mpi_bcast(pelist, size(pelist), MPI_INTEGER, 0, DRIVER_COMM, ierr)
       call mpi_bcast(exlist, size(exlist), MPI_INTEGER, 0, DRIVER_COMM, ierr)
       if (exlist(3) > 0) then
          call mpi_comm_group(DRIVER_COMM, mpigrp_world, ierr)
          call shr_mpi_chkerr(ierr,subname//' mpi_comm_group mpigrp_world')
          call mpi_group_range_incl(mpigrp_world, 1, exlist, exgrp, ierr)
          call shr_mpi_chkerr(ierr,subname//' mpi_group_range_incl COMPID')
          seq_comms(COMPID(n))%excl_group = exgrp
       endif
       if (present(drv_comm_id)) then
          call seq_comm_setcomm(COMPID(n),pelist,nthreads=comp_nthreads,iname=name,inst=drv_comm_id)
       else
          call seq_comm_setcomm(COMPID(n),pelist,nthreads=comp_nthreads,iname=name,inst=n,tinst=num_inst_comp)
       endif
       call seq_comm_joincomm(CPLID, COMPID(n), CPLCOMPID(n), 'CPL'//name, n, num_inst_comp)
    enddo
    call seq_comm_jcommarr(COMPID, ALLCOMPID, 'ALL'//name//'ID', 1, 1)
    call seq_comm_joincomm(CPLID, ALLCOMPID, CPLALLCOMPID, 'CPLALL'//name//'ID', 1, 1)

  end subroutine comp_comm_init

  subroutine comp_pelayout_init(numpes, ntasks, rootpe, pestride, nthreads, layout)
    integer,intent(in) :: numpes
    integer,intent(out) :: ntasks, rootpe, pestride, nthreads
    character(len=*),optional :: layout

    ntasks = numpes
    rootpe = 0
    pestride = 1
    nthreads = 1
    if(present(layout)) then
       layout = trim(layout_concurrent)
    endif
  end subroutine comp_pelayout_init

  !---------------------------------------------------------
  subroutine seq_comm_clean()
    ! Resets this module - freeing memory, etc.
    !
    ! This potentially allows seq_comm_init can be called again, e.g., from unit tests.
    !
    ! Also calls mct_world_clean, to be symmetric with the mct_world_init call from
    ! seq_comm_init.

    character(*), parameter :: subName =   '(seq_comm_clean) '
    !----------------------------------------------------------

    if (.not. seq_comm_mct_initialized) then
       write(logunit,*) trim(subname),' ERROR seq_comm_init has not been called '
       call shr_sys_abort()
    end if
    seq_comm_mct_initialized = .false.

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
    integer :: mpigrp, newgrp
    integer :: mpicom
    integer :: ntasks
    integer :: ierr, n
    character(len=seq_comm_namelen) :: cname
    logical :: set_suffix
    character(*),parameter :: subName =   '(seq_comm_setcomm) '

    if (ID < 1 .or. ID > ncomps) then
       write(logunit,*) subname,' ID out of range, abort ',ID
       call shr_sys_abort()
    endif

    call mpi_comm_group(DRIVER_COMM, mpigrp_world, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_group mpigrp_world')
    call mpi_group_range_incl(mpigrp_world, 1, pelist, mpigrp,ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_group_range_incl mpigrp')

    ! exclude tasks dedicated to other components
    do n = 2, ncomps
       if (seq_comms(n)%excl_group /= -1) then
          if (n == ID) cycle ! don't exclude self
          call mpi_group_difference(mpigrp, seq_comms(n)%excl_group, newgrp, ierr)
          call shr_mpi_chkerr(ierr,subname//' mpi_group_difference excl_group')
          mpigrp = newgrp
       endif
    enddo

    call mpi_comm_create(DRIVER_COMM, mpigrp, mpicom, ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_create mpigrp')

    ntasks = ((pelist(2,1) - pelist(1,1)) / pelist(3,1)) + 1

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
       write(logunit,F11) trim(subname),'  init ID ',ID,seq_comms(ID)%name, &
            ' pelist   =',pelist,' npes =',seq_comms(ID)%npes,' nthreads =',seq_comms(ID)%nthreads,&
            ' suffix =',trim(seq_comms(ID)%suffix)
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
    call mpi_comm_create(DRIVER_COMM, mpigrp, mpicom, ierr)
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
          write(logunit,F12) trim(subname),' init ID ',ID,seq_comms(ID)%name, &
               ' join IDs =',ID1,ID2,' npes =',seq_comms(ID)%npes, &
               ' nthreads =',seq_comms(ID)%nthreads, &
               ' cpl/cmp pes =',seq_comms(ID)%cplpe,seq_comms(ID)%cmppe
       else
          write(logunit,F13) trim(subname),' init ID ',ID,seq_comms(ID)%name, &
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
    integer :: n
    character(len=seq_comm_namelen) :: cname
    logical :: set_suffix
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
    ! The allcompid is created across multiple drivers.
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
          write(logunit,F14) trim(subname),' init ID ',ID,seq_comms(ID)%name, &
               ' join multiple comp IDs',' npes =',seq_comms(ID)%npes, &
               ' nthreads =',seq_comms(ID)%nthreads
       else
          write(logunit,F14) trim(subname),' init ID ',ID,seq_comms(ID)%name, &
               ' join multiple comp IDs',' npes =',seq_comms(ID)%npes, &
               ' nthreads =',seq_comms(ID)%nthreads
       endif
    endif

  end subroutine seq_comm_jcommarr

  !---------------------------------------------------------
  subroutine seq_comm_printcomms()

    implicit none
    character(*),parameter :: subName =   '(seq_comm_printcomms) '
    integer :: n,mype,npes,ierr

    call mpi_comm_size(DRIVER_COMM, npes  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_size comm_world')
    call mpi_comm_rank(DRIVER_COMM, mype  , ierr)
    call shr_mpi_chkerr(ierr,subname//' mpi_comm_rank comm_world')

    call shr_sys_flush(logunit)
    call mpi_barrier(DRIVER_COMM,ierr)
    if (mype == 0) then
       do n = 1,ncomps
          write(logunit,'(a,4i6,2x,3a)') trim(subName),n, &
               seq_comms(n)%gloroot,seq_comms(n)%npes,seq_comms(n)%nthreads, &
               trim(seq_comms(n)%name),':',trim(seq_comms(n)%suffix)
       enddo
       call shr_sys_flush(logunit)
    endif

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

    if (present(gloiam)) then
       if (ID > 0) then
          gloiam = seq_comms(ID)%gloiam
       else
          gloiam = -1
       end if
    endif

    if (present(gloroot)) then
       if (ID > 0) then
          gloroot = seq_comms(ID)%gloroot
       else
          gloroot = -1
       end if
    endif

    if (present(cplpe)) then
       if (ID > 0) then
          cplpe = seq_comms(ID)%cplpe
       else
          cplpe = -1
       end if
    endif

    if (present(cmppe)) then
       if (ID > 0) then
          cmppe = seq_comms(ID)%cmppe
       else
          cmppe = -1
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
    character(*),parameter :: subName =   '(seq_comm_getnthreads) '
#ifdef _OPENMP
    integer :: omp_get_num_threads
    seq_comm_getnthreads = -1

    !$OMP PARALLEL
    seq_comm_getnthreads = omp_get_num_threads()
    !$OMP END PARALLEL
#else
    seq_comm_getnthreads = -1
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
  integer function seq_comm_gloiam(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_gloiam) '

    if ((ID < 1) .or. (ID > ncomps)) then
       seq_comm_gloiam = -1
    else
       seq_comm_gloiam = seq_comms(ID)%gloiam
    end if

  end function seq_comm_gloiam
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
  integer function seq_comm_cplpe(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_cplpe) '

    if ((ID < 1) .or. (ID > ncomps)) then
       seq_comm_cplpe = -1
    else
       seq_comm_cplpe = seq_comms(ID)%cplpe
    end if

  end function seq_comm_cplpe
  !---------------------------------------------------------
  integer function seq_comm_cmppe(ID)

    implicit none
    integer,intent(in) :: ID
    character(*),parameter :: subName =   '(seq_comm_cmppe) '

    if ((ID < 1) .or. (ID > ncomps)) then
       seq_comm_cmppe = -1
    else
       seq_comm_cmppe = seq_comms(ID)%cmppe
    end if

  end function seq_comm_cmppe
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
