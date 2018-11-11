module dead_mct_mod

  use esmf            , only : esmf_clock
  use shr_kind_mod    , only : IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_sys_mod     , only : shr_sys_abort, shr_sys_flush
  use shr_const_mod   , only : shr_const_pi
  use mct_mod         , only : mct_gsmap, mct_ggrid, mct_avect, mct_ggrid_init, mct_gsmap_lsize, mct_ggrid_lsize
  use mct_mod         , only : mct_avect_lsize, MCT_AVECT_NRATTR, mct_avect_indexra,mct_avect_zero
  use mct_mod         , only : mct_ggrid_importiattr, mct_ggrid_importrattr, mct_gsmap_init, mct_aVect_init
  use mct_mod         , only : mct_gsmap_orderedpoints
  use dead_data_mod   , only : dead_grid_lat, dead_grid_lon, dead_grid_area, dead_grid_mask, dead_grid_frac, dead_grid_index
  use dead_mod        , only : dead_setnewgrid, dead_read_inparms
  use seq_flds_mod    , only : seq_flds_dom_coord, seq_flds_dom_other
  use seq_timemgr_mod , only : seq_timemgr_EClockGetData

  implicit none
  private
  save

  public  :: dead_init_mct, dead_run_mct, dead_final_mct
  private :: dead_domain_mct

!===============================================================================
contains
!===============================================================================

  subroutine dead_init_mct(model, Eclock, x2d, d2x, &
         flds_x2d, flds_d2x, &
         gsmap, ggrid, gbuf, mpicom, compid, my_task, master_task, &
         inst_index, inst_suffix, inst_name, logunit, nxg, nyg)

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*) , intent(in)    :: model
    type(ESMF_Clock) , intent(inout) :: EClock
    type(mct_aVect)  , intent(inout) :: x2d         ! driver -> dead
    type(mct_aVect)  , intent(inout) :: d2x         ! dead   -> driver
    character(len=*) , intent(in)    :: flds_x2d
    character(len=*) , intent(in)    :: flds_d2x
    type(mct_gsMap)  , pointer       :: gsMap       ! model global sep map (output)
    type(mct_gGrid)  , pointer       :: ggrid       ! model ggrid (output)
    real(r8)         , pointer       :: gbuf(:,:)   ! model grid (output)
    integer(IN)      , intent(in)    :: mpicom      ! mpi communicator
    integer(IN)      , intent(in)    :: compid      ! mct comp id
    integer(IN)      , intent(in)    :: my_task     ! my task in mpi communicator mpicom
    integer(IN)      , intent(in)    :: master_task ! task number of master task
    integer(IN)      , intent(in)    :: inst_index  ! instance number
    character(len=*) , intent(in)    :: inst_suffix ! char string associated with instance
    character(len=*) , intent(in)    :: inst_name   ! fullname of current instance (ie. "lnd_0001")
    integer(IN)      , intent(in)    :: logunit     ! logging unit number
    integer(IN)      , intent(out)   :: nxg         ! global dim i-direction
    integer(IN)      , intent(out)   :: nyg         ! global dim j-direction

    !--- local variables ---
    integer(IN)              :: ierr          ! error code
    integer(IN)              :: local_comm    ! local communicator
    integer(IN)              :: mype          ! pe info
    integer(IN)              :: totpe         ! total number of pes
    integer(IN), allocatable :: gindex(:)     ! global index
    integer(IN)              :: lsize
    integer(IN)              :: nproc_x
    integer(IN)              :: seg_len
    integer(IN)              :: decomp_type
    logical                  :: flood=.false. ! rof flood flag

    !--- formats ---
    character(*), parameter :: F00   = "('(',a,'_init_mct) ',8a)"
    character(*), parameter :: F01   = "('(',a,'_init_mct) ',a,4i8)"
    character(*), parameter :: F02   = "('(',a,'_init_mct) ',a,4es13.6)"
    character(*), parameter :: F90   = "('(',a,'_init_mct) ',73('='))"
    character(*), parameter :: F91   = "('(',a,'_init_mct) ',73('-'))"
    character(*), parameter :: subName = "(dead_init_mct) "
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

    ! Initialize MCT global seg map

    allocate(gindex(lsize))
    gindex(:) = gbuf(:,dead_grid_index)
    call mct_gsMap_init( gsMap, gindex, mpicom, compid, lsize, nxg*nyg )
    deallocate(gindex)

    ! Initialize MCT domain

    call dead_domain_mct(mpicom, gbuf, gsMap, logunit, ggrid)

    ! Initialize MCT attribute vectors

    call mct_aVect_init(d2x, rList=flds_d2x, lsize=lsize)
    call mct_aVect_zero(d2x)

    call mct_aVect_init(x2d, rList=flds_x2d, lsize=lsize)
    call mct_aVect_zero(x2d)

  end subroutine dead_init_mct

  !===============================================================================
  subroutine dead_run_mct(model, EClock, x2d, d2x, &
       gsmap, ggrid, gbuf, mpicom, compid, my_task, master_task, logunit)

    implicit none

    ! !DESCRIPTION: run method for dead model

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*) , intent(in)    :: model
    type(ESMF_Clock) , intent(inout) :: EClock
    type(mct_aVect)  , intent(inout) :: x2d         ! driver -> dead
    type(mct_aVect)  , intent(inout) :: d2x         ! dead   -> driver
    type(mct_gsMap)  , pointer       :: gsMap       ! model global sep map (output)
    type(mct_gGrid)  , pointer       :: ggrid       ! model ggrid (output)
    real(r8)         , pointer       :: gbuf(:,:)   ! model grid
    integer(IN)      , intent(in)    :: mpicom      ! mpi communicator
    integer(IN)      , intent(in)    :: compid      ! mct comp id
    integer(IN)      , intent(in)    :: my_task     ! my task in mpi communicator mpicom
    integer(IN)      , intent(in)    :: master_task ! task number of master task
    integer(IN)      , intent(in)    :: logunit     ! logging unit number

    !--- local ---
    integer(IN)             :: CurrentYMD        ! model date
    integer(IN)             :: CurrentTOD        ! model sec into model date
    integer(IN)             :: n                 ! index
    integer(IN)             :: nf                ! fields loop index
    integer(IN)             :: ki                ! index of ifrac
    integer(IN)             :: lsize             ! size of AttrVect
    real(R8)                :: lat               ! latitude
    real(R8)                :: lon               ! longitude
    integer                 :: nflds_d2x, nflds_x2d
    integer                 :: ncomp
    character(*), parameter :: F04   = "('(',a,'_run_mct) ',2a,2i8,'s')"
    character(*), parameter :: subName = "(dead_run_mct) "
    !-------------------------------------------------------------------------------

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

    lsize = mct_avect_lsize(x2d)
    nflds_d2x = mct_avect_nRattr(d2x)
    nflds_x2d = mct_avect_nRattr(x2d)

    if (model.eq.'rof') then

       do nf=1,nflds_d2x
          do n=1,lsize
             d2x%rAttr(nf,n) = (nf+1) * 1.0_r8
          enddo
       enddo

    else if (model.eq.'glc') then

       do nf=1,nflds_d2x
          do n=1,lsize
             lon = gbuf(n,dead_grid_lon)
             lat = gbuf(n,dead_grid_lat)
             d2x%rAttr(nf,n) = (nf*100)                   &
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
             d2x%rAttr(nf,n) = (nf*100)                   &
                  *  cos (SHR_CONST_PI*lat/180.0_R8)       &
                  *  sin((SHR_CONST_PI*lon/180.0_R8)       &
                  -      (ncomp-1)*(SHR_CONST_PI/3.0_R8) ) &
                  + (ncomp*10.0_R8)
          enddo
       enddo

    endif

    selectcase(model)
    case('ice')

       ki = mct_aVect_indexRA(d2x,"Si_ifrac",perrWith=subname)
       d2x%rAttr(ki,:) = min(1.0_R8,max(0.0_R8,d2x%rAttr(ki,:)))

    case('glc')

       ki = mct_aVect_indexRA(d2x,"Sg_icemask",perrWith=subname)
       d2x%rAttr(ki,:) = 1.0_R8

       ki = mct_aVect_indexRA(d2x,"Sg_icemask_coupled_fluxes",perrWith=subname)
       d2x%rAttr(ki,:) = 1.0_R8

       ki = mct_aVect_indexRA(d2x,"Sg_ice_covered",perrWith=subname)
       d2x%rAttr(ki,:) = 1.0_R8

    end select

    ! log output for model date

    if (my_task == master_task) then
       call seq_timemgr_EClockGetData(EClock,curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
       write(logunit,F04) model,trim(model),': model date ', CurrentYMD,CurrentTOD
       call shr_sys_flush(logunit)
    end if

  end subroutine dead_run_mct

  !===============================================================================
  subroutine dead_final_mct(model, my_task, master_task, logunit)

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

  end subroutine dead_final_mct

  !===============================================================================
  subroutine dead_domain_mct( mpicom, gbuf, gsMap, logunit, domain )

    !-------------------------------------------------------------------
    !---arguments---
    integer(IN)    , intent(in)  :: mpicom
    real(R8)       , intent(in)  :: gbuf(:,:)
    type(mct_gsMap), intent(in)  :: gsMap
    integer(IN)    , intent(in)  :: logunit
    type(mct_ggrid), intent(out) :: domain

    !---local variables---
    integer(IN)          :: my_task     ! mpi task within communicator
    integer(IN)          :: lsize       ! temporary
    integer(IN)          :: n    ,j,i   ! indices
    integer(IN)          :: ier         ! error status
    real(R8), pointer    :: data(:)     ! temporary
    integer(IN), pointer :: idata(:)    ! temporary
    !-------------------------------------------------------------------
    !
    ! Initialize mct dead domain
    !
    call mct_gGrid_init( GGrid=domain, CoordChars=trim(seq_flds_dom_coord), &
         OtherChars=trim(seq_flds_dom_other), &
         lsize=mct_gsMap_lsize(gsMap, mpicom) )
    call mct_aVect_zero(domain%data)
    !
    ! Allocate memory
    !
    lsize = mct_gGrid_lsize(domain)
    if (size(gbuf,dim=1) /= lsize) then
       call shr_sys_abort('mct_dead_domain size error')
    endif
    allocate(data(lsize))
    allocate(idata(lsize))
    !
    ! Initialize attribute vector with special value
    !
    call mpi_comm_rank(mpicom, my_task, ier)
    call mct_gsMap_orderedPoints(gsMap, my_task, idata)
    call mct_gGrid_importIAttr(domain,'GlobGridNum',idata,lsize)
    !
    call mct_aVect_zero(domain%data)
    data(:) = -9999.0_R8  ! generic special value
    call mct_gGrid_importRAttr(domain,"lat" ,data,lsize)
    call mct_gGrid_importRAttr(domain,"lon" ,data,lsize)
    call mct_gGrid_importRAttr(domain,"area",data,lsize)
    call mct_gGrid_importRAttr(domain,"frac",data,lsize)

    data(:) = 0.0_R8  ! generic special value
    call mct_gGrid_importRAttr(domain,"mask" ,data,lsize)
    call mct_gGrid_importRAttr(domain,"aream",data,lsize)
    !
    ! Fill in correct values for domain components
    !
    do n = 1,lsize
       data(n) = gbuf(n,dead_grid_lat)
    enddo
    call mct_gGrid_importRAttr(domain,"lat",data,lsize)

    do n = 1,lsize
       data(n) = gbuf(n,dead_grid_lon)
    enddo
    call mct_gGrid_importRAttr(domain,"lon",data,lsize)

    do n = 1,lsize
       data(n) = gbuf(n,dead_grid_area)
    enddo
    call mct_gGrid_importRAttr(domain,"area",data,lsize)
    call mct_gGrid_importRAttr(domain,"aream",data,lsize)

    do n = 1,lsize
       data(n) = gbuf(n,dead_grid_mask)
    enddo
    call mct_gGrid_importRAttr(domain,"mask"   ,data,lsize)

    do n = 1,lsize
       data(n) = gbuf(n,dead_grid_frac)
    enddo
    call mct_gGrid_importRAttr(domain,"frac"   ,data,lsize)

    deallocate(data)
    deallocate(idata)

  end subroutine dead_domain_mct
  !===============================================================================

end module dead_mct_mod
