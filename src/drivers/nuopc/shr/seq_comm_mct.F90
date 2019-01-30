module seq_comm_mct

  !---------------------------------------------------------------------
  ! Purpose: Set up necessary communications
  !---------------------------------------------------------------------
  !NOTE: If all atmospheres are identical in number of processes,
  !number of threads, and grid layout, we should check that the
  !user-provided number of processes and threads are consistent
  !(or else, only accept one entry for these quantities when reading
  !the namelist).  ARE OTHER PROTECTIONS/CHECKS NEEDED???
  
  use shr_sys_mod, only : shr_sys_abort

  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public seq_comm_setcomm
  public seq_comm_setptrs
  public seq_comm_name
  public seq_comm_suffix
  public seq_comm_inst

  private seq_comm_mkname

  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------

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
    integer, pointer    :: petlist(:)  ! esmf pet list
    logical :: petlist_allocated ! whether the petlist pointer variable was allocated
  end type seq_comm_type

  type(seq_comm_type) :: seq_comms(ncomps)
  integer             :: Global_Comm
  integer, public     :: logunit  = 6     ! log unit number

!=======================================================================
contains
!======================================================================

  subroutine seq_comm_setcomm(ID,pelist,nthreads,iname,inst,tinst, comm_in)

    use mpi         , only : MPI_COMM_NULL, mpi_comm_group, mpi_comm_create, mpi_group_range_incl
    use mpi         , only : mpi_comm_size, mpi_comm_rank
    use shr_mpi_mod , only : shr_mpi_chkerr

    integer          ,intent(in)          :: ID
    integer          ,intent(in)          :: pelist(:,:)
    integer          ,intent(in),optional :: nthreads
    character(len=*) ,intent(in),optional :: iname  ! name of component
    integer          ,intent(in),optional :: inst  ! instance of component
    integer          ,intent(in),optional :: tinst ! total number of instances for this component
    integer          ,intent(in),optional :: comm_in

    integer :: mpigrp_world
    integer :: mpigrp
    integer :: mpicom
    integer :: ntask,ntasks,cnt
    integer :: ierr
    logical :: set_suffix
    character(len=seq_comm_namelen) :: cname
    character(*), parameter :: F11 = "(a,a,'(',i3,' ',a,')',a,   3i6,' (',a,i6,')',' (',a,i3,')')"
    character(*), parameter :: subName =   '(seq_comm_setcomm) '

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
  subroutine seq_comm_setptrs(ID,mpicom,mpigrp,npes,nthreads,iam,iamroot,gloroot, pethreads, name)

    use mpi, only : mpi_comm_null, mpi_group_null

    integer                         ,intent(in)            :: ID
    integer                         ,intent(out), optional :: mpicom
    integer                         ,intent(out), optional :: mpigrp
    integer                         ,intent(out), optional :: npes
    integer                         ,intent(out), optional :: nthreads
    integer                         ,intent(out), optional :: iam
    logical                         ,intent(out), optional :: iamroot
    integer                         ,intent(out), optional :: gloroot
    integer                         ,intent(out), optional :: pethreads
    character(len=seq_comm_namelen) ,intent(out), optional :: name

    character(*),parameter :: subName =   '(seq_comm_setptrs) '

    ! Negative ID means there is no comm, return default or inactive values
    if ((ID == 0) .or. (ID > ncomps)) then
       write(logunit,*) 'seq_comm_setptrs: ID out of range, return ',ID
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
  character(len=seq_comm_namelen) function seq_comm_suffix(ID)
    integer,intent(in) :: ID
    if ((ID < 1) .or. (ID > ncomps)) then
       seq_comm_suffix = ''
    else
       seq_comm_suffix = trim(seq_comms(ID)%suffix)
    end if
  end function seq_comm_suffix

  !---------------------------------------------------------
  integer function seq_comm_inst(ID)
    integer,intent(in) :: ID
    if ((ID < 1) .or. (ID > ncomps)) then
      seq_comm_inst = 0
    else
      seq_comm_inst = seq_comms(ID)%inst
    end if
  end function seq_comm_inst

  !---------------------------------------------------------
  character(len=seq_comm_namelen) function seq_comm_name(ID)
    integer,intent(in) :: ID
    if ((ID < 1) .or. (ID > ncomps)) then
       seq_comm_name = ''
    else
       seq_comm_name = trim(seq_comms(ID)%name)
    end if
  end function seq_comm_name

  !---------------------------------------------------------
  subroutine seq_comm_mkname(oname,str1,num)

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

end module seq_comm_mct
