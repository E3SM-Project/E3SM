module dead_nuopc_mod

 use esmf               , only : esmf_clock
 use shr_kind_mod       , only : IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
 use shr_sys_mod        , only : shr_sys_abort, shr_sys_flush
 use shr_const_mod      , only : shr_const_pi
 use shr_string_mod     , only : shr_string_listGetIndexF
 use dead_data_mod      , only : dead_grid_lat, dead_grid_lon, dead_grid_area
 use dead_data_mod      , only : dead_grid_mask, dead_grid_frac, dead_grid_index
 use dead_mod           , only : dead_setnewgrid, dead_read_inparms
 use seq_timemgr_mod    , only : seq_timemgr_EClockGetData

  implicit none
  private
  save

  public  :: dead_init_nuopc, dead_run_nuopc, dead_final_nuopc

!===============================================================================
contains
!===============================================================================

  subroutine dead_init_nuopc(model, mpicom, my_task, master_task, &
         inst_index, inst_suffix, inst_name, logunit, lsize, gbuf, nxg, nyg)

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*) , intent(in)    :: model
    integer(IN)      , intent(in)    :: mpicom      ! mpi communicator
    integer(IN)      , intent(in)    :: my_task     ! my task in mpi communicator mpicom
    integer(IN)      , intent(in)    :: master_task ! task number of master task
    integer(IN)      , intent(in)    :: inst_index  ! instance number
    character(len=*) , intent(in)    :: inst_suffix ! char string associated with instance
    character(len=*) , intent(in)    :: inst_name   ! fullname of current instance (ie. "lnd_0001")
    integer(IN)      , intent(in)    :: logunit     ! logging unit number
    integer(IN)      , intent(out)   :: lsize       ! logging unit number
    real(r8)         , pointer       :: gbuf(:,:)   ! model grid
    integer          , pointer       :: gindex(:)   ! global index space
    integer(IN)      , intent(out)   :: nxg         ! global dim i-direction
    integer(IN)      , intent(out)   :: nyg         ! global dim j-direction

    !--- local variables ---
    integer(IN)              :: ierr          ! error code
    integer(IN)              :: local_comm    ! local communicator
    integer(IN)              :: mype          ! pe info
    integer(IN)              :: totpe         ! total number of pes
    integer(IN)              :: nproc_x
    integer(IN)              :: seg_len
    integer(IN)              :: decomp_type
    logical                  :: flood=.false. ! rof flood flag

    !--- formats ---
    character(*), parameter :: F00   = "('(',a,'_init_nuopc) ',8a)"
    character(*), parameter :: F01   = "('(',a,'_init_nuopc) ',a,4i8)"
    character(*), parameter :: F02   = "('(',a,'_init_nuopc) ',a,4es13.6)"
    character(*), parameter :: F90   = "('(',a,'_init_nuopc) ',73('='))"
    character(*), parameter :: F91   = "('(',a,'_init_nuopc) ',73('-'))"
    character(*), parameter :: subName = "(dead_init_nuopc) "
    !-------------------------------------------------------------------------------

    ! Determine communicator groups and sizes

    local_comm = mpicom
    call MPI_COMM_RANK(local_comm,mype ,ierr)
    call MPI_COMM_SIZE(local_comm,totpe,ierr)

    ! Read input parms

    call dead_read_inparms(model, mpicom, my_task, master_task, &
       inst_index, inst_suffix, inst_name, logunit, &
       nxg, nyg, decomp_type, nproc_x, seg_len, flood)

    ! Initialize grid

    call dead_setNewGrid(decomp_type, nxg, nyg, totpe, mype, logunit, &
                         lsize, gbuf, seg_len, nproc_x)

  end subroutine dead_init_nuopc

  !===============================================================================
  subroutine dead_run_nuopc(model, EClock, x2d, d2x, gbuf, flds_d2x, my_task, master_task, logunit)

    implicit none

    ! !DESCRIPTION: run method for dead model

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*) , intent(in)    :: model
    type(ESMF_Clock) , intent(inout) :: EClock
    real(r8)         , intent(inout) :: x2d(:,:)    ! driver -> dead
    real(r8)         , intent(inout) :: d2x(:,:)    ! dead   -> driver
    real(r8)         , pointer       :: gbuf(:,:)   ! model grid
    character(len=*) , intent(in)    :: flds_d2x    ! list of fields to dead -> driver
    integer(IN)      , intent(in)    :: my_task     ! my task in mpi communicator mpicom
    integer(IN)      , intent(in)    :: master_task ! task number of master task
    integer(IN)      , intent(in)    :: logunit     ! logging unit number

    !--- local ---
    integer(IN)             :: CurrentYMD        ! model date
    integer(IN)             :: CurrentTOD        ! model sec into model date
    integer(IN)             :: n                 ! index
    integer(IN)             :: nf                ! fields loop index
    integer(IN)             :: ki                ! index
    integer(IN)             :: lsize             ! size of AttrVect
    real(R8)                :: lat               ! latitude
    real(R8)                :: lon               ! longitude
    integer                 :: nflds_x2d
    integer                 :: nflds_d2x
    integer                 :: ncomp
    character(*), parameter :: F04   = "('(',a,'_run_nuopc) ',2a,2i8,'s')"
    character(*), parameter :: subName = "(dead_run_nuopc) "
    !-------------------------------------------------------------------------------

    nflds_x2d = size(x2d, dim=1)
    nflds_d2x = size(d2x, dim=1)

    ! PACK (currently no unpacking)

    selectcase(model)
    case('atm')
       ncomp = 1
    case('lnd')
       ncomp = 2
    case('ice')
       ncomp = 3
    case('ocn')
       ncomp = 4
    case('glc')
       ncomp = 5
    case('rof')
       ncomp = 6
    case('wav')
       ncomp = 7
    end select

    nflds_d2x = size(d2x, dim=1)
    nflds_x2d = size(x2d, dim=1)
    lsize = size(d2x, dim=2)

    if (model.eq.'rof') then

       do nf=1,nflds_d2x
          do n=1,lsize
             d2x(nf,n) = (nf+1) * 1.0_r8
          enddo
       enddo

    else if (model.eq.'glc') then

       do nf=1,nflds_d2x
          do n=1,lsize
             lon = gbuf(n,dead_grid_lon)
             lat = gbuf(n,dead_grid_lat)
             d2x(nf,n) = (nf*100)                          &
                  *  cos (SHR_CONST_PI*lat/180.0_R8)       &
                  *  cos (SHR_CONST_PI*lat/180.0_R8)       &
                  *  sin (SHR_CONST_PI*lon/180.0_R8)       &
                  *  sin (SHR_CONST_PI*lon/180.0_R8)       &
                  + (ncomp*10.0_R8)
          enddo
       enddo

    else

       do nf=1,nflds_d2x
          do n=1,lsize
             lon = gbuf(n,dead_grid_lon)
             lat = gbuf(n,dead_grid_lat)
             d2x(nf,n) = (nf*100)                          &
                  *  cos (SHR_CONST_PI*lat/180.0_R8)       &
                  *  sin((SHR_CONST_PI*lon/180.0_R8)       &
                  -      (ncomp-1)*(SHR_CONST_PI/3.0_R8) ) &
                  + (ncomp*10.0_R8)
          enddo
       enddo

    endif

    selectcase(model)
    case('ice')

       ki = shr_string_listGetIndexF(flds_d2x, "Si_ifrac")
       d2x(ki,:) = min(1.0_R8,max(0.0_R8,d2x(ki,:)))

       ki = shr_string_listGetIndexF(flds_d2x, "Si_imask")
       d2x(ki,:) = float(nint(min(1.0_R8,max(0.0_R8,d2x(ki,:)))))

    case('ocn')

       ki = shr_string_listGetIndexF(flds_d2x, "So_omask")
       d2x(ki,:) = float(nint(min(1.0_R8,max(0.0_R8,d2x(ki,:)))))

    case('lnd')

       ki = shr_string_listGetIndexF(flds_d2x, "Sl_lfrin")
       d2x(ki,:) = 1.0_R8

    case('glc')

       ki = shr_string_listGetIndexF(flds_d2x, "Sg_icemask")
       d2x(ki,:) = 1.0_R8

       ki = shr_string_listGetIndexF(flds_d2x, "Sg_icemask_coupled_fluxes")
       d2x(ki,:) = 1.0_R8

       ki = shr_string_listGetIndexF(flds_d2x, "Sg_ice_covered")
       d2x(ki,:) = 1.0_R8

    end select

    ! log output for model date

    if (my_task == master_task) then
       call seq_timemgr_EClockGetData(EClock,curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
       write(logunit,F04) model,trim(model),': model date ', CurrentYMD,CurrentTOD
       call shr_sys_flush(logunit)
    end if

  end subroutine dead_run_nuopc

  !===============================================================================
  subroutine dead_final_nuopc(model, my_task, master_task, logunit)

    ! !DESCRIPTION: finalize method for datm model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*) , intent(in) :: model
    integer(IN)      , intent(in) :: my_task     ! my task in mpi communicator mpicom
    integer(IN)      , intent(in) :: master_task ! task number of master task
    integer(IN)      , intent(in) :: logunit     ! logging unit number

    !--- formats ---
    character(*), parameter :: F00   = "('(dead_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(dead_comp_final) ',73('-'))"
    character(*), parameter :: subName = "(dead_comp_final) "
    !-------------------------------------------------------------------------------

    if (my_task == master_task) then
       write(logunit,F91)
       write(logunit,F00) trim(model),': end of main integration loop'
       write(logunit,F91)
    end if

  end subroutine dead_final_nuopc

end module dead_nuopc_mod
