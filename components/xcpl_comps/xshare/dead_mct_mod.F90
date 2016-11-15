module dead_mct_mod

  use shr_kind_mod		, only : IN =>SHR_KIND_IN, R8=>SHR_KIND_R8
  use shr_sys_mod		, only : shr_sys_abort, shr_sys_flush
  use seq_flds_mod		, only : seq_flds_dom_coord, seq_flds_dom_other
  use shr_file_mod		, only : shr_file_getlogunit, shr_file_getunit, shr_file_setio, shr_file_getloglevel, &
                                         shr_file_setlogunit, shr_file_freeunit,shr_file_setloglevel
  use esmf                      , only : esmf_clock
  use seq_cdata_mod		, only : seq_cdata, seq_cdata_setptrs
  use seq_infodata_mod		, only : seq_infodata_type, seq_infodata_getData, seq_infodata_PutData
  use mct_mod			, only : mct_gsmap, mct_ggrid, mct_avect, mct_ggrid_init, mct_gsmap_lsize, mct_ggrid_lsize, &
                                         mct_avect_lsize, MCT_AVECT_NRATTR, mct_avect_indexra,mct_avect_zero, &
                                         mct_ggrid_importiattr, mct_ggrid_importrattr, mct_gsmap_init, mct_aVect_init, &
                                         mct_gsmap_orderedpoints
  use dead_data_mod		, only : dead_grid_lat, dead_grid_lon, dead_grid_area, dead_grid_mask, dead_grid_frac, dead_grid_index
  use dead_mod			, only : dead_setnewgrid
  use seq_comm_mct		, only : seq_comm_suffix, seq_comm_name,seq_comm_inst, seq_comm_get_ncomps
  use shr_const_mod		, only : shr_const_pi
  use seq_timemgr_mod		, only : seq_timemgr_EClockGetData
  use shr_mpi_mod		, only : shr_mpi_bcast
  implicit none
  private
  save 

  public		:: dead_domain_mct, dead_init_mct, dead_run_mct, dead_final_mct
  integer, parameter	:: master_task=0
  integer, pointer	:: logunits(:) => null()

  interface dead_domain_mct ; module procedure &
       dead_domain_mct_new
  end interface

  real(r8), pointer	:: gbuf_atm(:,:)
  real(r8), pointer	:: gbuf_lnd(:,:)
  real(r8), pointer	:: gbuf_ice(:,:)
  real(r8), pointer	:: gbuf_ocn(:,:)
  real(r8), pointer	:: gbuf_glc(:,:)
  real(r8), pointer	:: gbuf_rof(:,:)
  real(r8), pointer	:: gbuf_wav(:,:)

  !===============================================================================
contains
  !===============================================================================

  subroutine dead_domain_mct_new( mpicom, gbuf, gsMap, domain )

    !-------------------------------------------------------------------
    !---arguments--- 
    integer(IN)		, intent(in)		:: mpicom
    real(R8)		, intent(in)		:: gbuf(:,:)
    type(mct_gsMap)	, intent(in)		:: gsMap
    type(mct_ggrid)	, intent(out)	        :: domain  

    !---local variables---
    integer(IN)				:: my_task               ! mpi task within communicator
    integer(IN)				:: lsize                 ! temporary
    integer(IN)				:: n	,j,i                 ! indices	
    integer(IN)				:: ier                   ! error status
    real(R8), pointer		        :: data(:)     ! temporary
    integer(IN), pointer		:: idata(:)    ! temporary
    integer(IN)				:: logunit
    !-------------------------------------------------------------------

    call shr_file_getlogunit(logunit)

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

  end subroutine dead_domain_mct_new

  subroutine dead_init_mct( model, flds_d2x, flds_x2d, EClock, cdata, x2d, d2x, NLFilename )

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*), intent(in) :: model, flds_d2x, flds_x2d

    type(ESMF_Clock)            , intent(inout) :: EClock
    type(seq_cdata)             , intent(inout) :: cdata
    type(mct_aVect)             , intent(inout) :: x2d, d2x
    character(len=*), optional  , intent(in)    :: NLFilename ! Namelist filename

    !EOP

    !--- local variables ---
    integer(IN)   :: unitn       ! Unit for namelist file
    integer(IN)   :: ierr        ! error code
    integer(IN)   :: local_comm  ! local communicator
    integer(IN)   :: mype        ! pe info
    integer(IN)   :: totpe       ! total number of pes
    integer(IN), allocatable :: gindex(:)  ! global index

    integer(IN)                         :: COMPID
    integer(IN)                         :: mpicom
    integer(IN)                         :: lsize
    type(mct_gsMap)        , pointer	:: gsMap
    type(mct_gGrid)        , pointer	:: dom
    type(seq_infodata_type), pointer	:: infodata   ! Input init object
    integer(IN)				:: shrlogunit, shrloglev ! original log unit and level
    integer(IN)				:: nproc_x
    integer(IN)				:: seg_len
    integer(IN)				:: nxg           ! global dim i-direction
    integer(IN)				:: nyg           ! global dim j-direction
    integer(IN)				:: decomp_type
    integer(IN)				:: my_task
    character(len=16)			:: inst_suffix       ! char string associated with instance 
    integer(IN)				:: phase
    integer(IN)				:: logunit
    logical				:: flood=.false.     ! rof flood flag
    real(r8), pointer			:: gbuf(:,:)

    !--- formats ---
    character(*), parameter		:: F00   = "('(',a,'_init_mct) ',8a)"
    character(*), parameter		:: F01   = "('(',a,'_init_mct) ',a,4i8)"
    character(*), parameter		:: F02   = "('(',a,'_init_mct) ',a,4es13.6)"
    character(*), parameter		:: F03   = "('(',a,'_init_mct) ',a,i8,a)"
    character(*), parameter              :: F04   = "('(rof_init_mct) ',a,l4)"
    character(*), parameter		:: F90   = "('(',a,'_init_mct) ',73('='))"
    character(*), parameter		:: F91   = "('(',a,'_init_mct) ',73('-'))"
    character(*), parameter		:: subName = "(dead_init_mct) "
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    ! Set cdata pointers

    call seq_cdata_setptrs(cdata, ID=COMPID, mpicom=mpicom, &
         gsMap=gsMap, dom=dom, infodata=infodata)

    call mpi_comm_rank(mpicom, my_task, ierr)


    if(model .eq. 'atm') then
       call seq_infodata_getData(infodata,atm_phase=phase)
       if (phase > 1) return
    endif


    inst_suffix = seq_comm_suffix(COMPID)
    if(.not. associated(logunits)) then
       allocate(logunits(seq_comm_get_ncomps()))
    endif
    !--- open log file ---
    if (my_task == master_task) then
       logunits(compid) = shr_file_getUnit()
       call shr_file_setIO(model//'_modelio.nml'//trim(inst_suffix),logunits(compid))
    else
       logunits(compid) = 6
    endif

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logunits(compid))

    ! read the namelist input (used to configure model)

    nxg            =  -9999
    nyg            =  -9999
    nproc_x        =  -9999
    seg_len        =  -9999
    decomp_type    =  -9999


    if (my_task == master_task) then
       unitn = shr_file_getUnit()
       open( unitn, file='x'//model//'_in'//trim(inst_suffix), status='old' )
       read(unitn,*) nxg
       read(unitn,*) nyg
       read(unitn,*) decomp_type
       read(unitn,*) nproc_x
       read(unitn,*) seg_len
       if(model.eq.'rof') then
          read(unitn,*) flood
       end if
       close (unitn)
       call shr_file_freeunit(unitn)
    endif

    call shr_mpi_bcast(nxg        ,mpicom,'x'//model//' nxg')
    call shr_mpi_bcast(nyg        ,mpicom,'x'//model//' nyg')
    call shr_mpi_bcast(decomp_type,mpicom,'x'//model//' decomp_type')
    call shr_mpi_bcast(nproc_x    ,mpicom,'x'//model//' nproc_x')
    call shr_mpi_bcast(seg_len    ,mpicom,'x'//model//' seg_len')
    if(model.eq.'rof') then
       call shr_mpi_bcast(flood      ,mpicom,'xrof flood')
    end if

    if (my_task == master_task) then
       write(logunits(compid),*)   ' Read in X'//model//' input from file= x'//model//'_in'
       write(logunits(compid),F00) model
       write(logunits(compid),F00) model,'         Model  :  ',model
       write(logunits(compid),F01) model,'           NGX  :  ',nxg
       write(logunits(compid),F01) model,'           NGY  :  ',nyg
       write(logunits(compid),F01) model,' Decomposition  :  ',decomp_type
       write(logunits(compid),F03) model,' Num pes in X   :  ',nproc_x,'  (type 3 only)'
       write(logunits(compid),F03) model,' Segment Length :  ',seg_len,'  (type 11 only)'
       write(logunits(compid),F01) model,'    inst_index  :  ',   seq_comm_inst(COMPID)
       write(logunits(compid),F00) model,'    inst_name   :  ',trim(seq_comm_name(COMPID))
       write(logunits(compid),F00) model,'    inst_suffix :  ',trim(inst_suffix)
       if(model.eq.'rof') then
          write(logunits(compid),F04) ' Flood mode     :  ',flood
       endif
       write(logunits(compid),F00) model
       call shr_sys_flush(logunits(compid))
    end if

    ! Set flag to specify dead components
    selectcase(model)
    case('ice')
       if (nxg == 0 .and. nyg == 0) then
          call seq_infodata_PutData(infodata, dead_comps=.true., ice_present=.false., &
               ice_prognostic = .false., iceberg_prognostic=.false., ice_nx=nxg, ice_ny=nyg)
       else
          call seq_infodata_PutData(infodata, dead_comps=.true., ice_present=.true., &
               ice_prognostic = .true., iceberg_prognostic=.true., ice_nx=nxg, ice_ny=nyg)
       endif
    case('atm')
       if (nxg == 0 .and. nyg == 0) then
          call seq_infodata_PutData( infodata, dead_comps=.true., atm_present=.false., &
               atm_prognostic=.false., atm_nx=nxg, atm_ny=nyg)
       else
          call seq_infodata_PutData( infodata, dead_comps=.true., atm_present=.true., &
               atm_prognostic=.true., atm_nx=nxg, atm_ny=nyg)
       endif
    case('lnd')
       if (nxg == 0 .and. nyg == 0) then
          call seq_infodata_PutData(infodata, dead_comps=.true., lnd_present=.false., &
               lnd_prognostic=.false., lnd_nx=nxg, lnd_ny=nyg)
       else
          call seq_infodata_PutData(infodata, dead_comps=.true., lnd_present=.true., &
               lnd_prognostic=.true., lnd_nx=nxg, lnd_ny=nyg)
       endif
    case('wav')
       if (nxg == 0 .and. nyg == 0) then
          call seq_infodata_PutData( infodata, dead_comps=.true.,wav_present=.false., &
               wav_prognostic=.false., wav_nx=nxg, wav_ny=nyg)
       else
          call seq_infodata_PutData( infodata, dead_comps=.true.,wav_present=.true., &
               wav_prognostic=.true., wav_nx=nxg, wav_ny=nyg)
       endif
    case('glc')
       if (nxg == 0 .and. nyg == 0) then
          call seq_infodata_PutData(infodata, dead_comps=.true., glc_present=.false., &
               glclnd_present = .false., glcocn_present=.false., glcice_present=.false., &
               glc_prognostic = .false., glc_nx=nxg, glc_ny=nyg)
       else
          call seq_infodata_PutData(infodata, dead_comps=.true., glc_present=.true., &
               glclnd_present = .true., glcocn_present=.false., glcice_present=.false., &
               glc_prognostic = .true., glc_nx=nxg, glc_ny=nyg)
       endif
    case('rof')
       if (nxg == 0 .and. nyg == 0) then
          call seq_infodata_PutData(infodata, dead_comps=.true., rof_present=.false., &
               rof_prognostic = .false., rofice_present = .false., flood_present = .false., &
               rof_nx=nxg, rof_ny=nyg)
       else
          call seq_infodata_PutData(infodata, dead_comps=.true., rof_present=.true., &
               rof_prognostic = .true., rofice_present = .false., flood_present = flood, &
               rof_nx=nxg, rof_ny=nyg)
       endif
    case('ocn')
       if (nxg == 0 .and. nyg == 0) then
          call seq_infodata_PutData( infodata, dead_comps=.true.,ocn_present=.false., &
               ocn_prognostic=.false., ocnrof_prognostic=.false., ocn_nx=nxg, ocn_ny=nyg)
       else
          call seq_infodata_PutData( infodata, dead_comps=.true.,ocn_present=.true., &
               ocn_prognostic=.true., ocnrof_prognostic=.true., ocn_nx=nxg, ocn_ny=nyg)
       endif

    end select

    ! Determine communicator groups and sizes

    local_comm = mpicom
    call MPI_COMM_RANK(local_comm,mype ,ierr)
    call MPI_COMM_SIZE(local_comm,totpe,ierr)

    ! Initialize grid

    call dead_setNewGrid(decomp_type,nxg,nyg,totpe,mype,lsize,gbuf,seg_len,nproc_x)

    selectcase(model)
    case('ice')
       gbuf_ice => gbuf
    case('atm')
       gbuf_atm => gbuf
    case('lnd')
       gbuf_lnd => gbuf
    case('wav')
       gbuf_wav => gbuf
    case('glc')
       gbuf_glc => gbuf
    case('rof')
       gbuf_rof => gbuf
    case('ocn')
       gbuf_ocn => gbuf
    end select

    ! Initialize MCT global seg map

    allocate(gindex(lsize))
    gindex(:) = gbuf(:,dead_grid_index)
    call mct_gsMap_init( gsMap, gindex, mpicom, COMPID, lsize, nxg*nyg)
    deallocate(gindex)

    ! Initialize MCT domain

    call dead_domain_mct(mpicom, gbuf, gsMap, dom)

    ! Initialize MCT attribute vectors

    call mct_aVect_init(d2x, rList=flds_d2x, lsize=lsize)
    call mct_aVect_zero(d2x)

    call mct_aVect_init(x2d, rList=flds_x2d, lsize=lsize)
    call mct_aVect_zero(x2d)

    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunits(compid))

  end subroutine dead_init_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: dead_run_mct
  !
  ! !DESCRIPTION:
  !     run method for dead model
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine dead_run_mct( model, EClock, cdata, x2d, d2x)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*)            ,intent(in)    :: model
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d        ! driver -> dead
    type(mct_aVect)             ,intent(inout) :: d2x        ! dead   -> driver

    !EOP
    !--- local ---
    integer(IN)				:: CurrentYMD        ! model date
    integer(IN)				:: CurrentTOD        ! model sec into model date
    integer(IN)				:: n                 ! index
    integer(IN)				:: nf                ! fields loop index
    integer(IN)				:: ki                ! index of ifrac
    integer(IN)				:: lsize             ! size of AttrVect
    real(R8)				:: lat               ! latitude
    real(R8)				:: lon               ! longitude
    real(r8), pointer			:: gbuf(:,:)
    integer(IN)				:: shrlogunit, shrloglev ! original log unit and level
    integer				:: nflds_d2x, nflds_x2d
    integer				:: ncomp
    integer				:: compid
    integer				:: mpicom, ierr, my_task
    real(R8)				:: nextsw_cday ! calendar of next atm shortwave
    type(seq_infodata_type), pointer	:: infodata   ! Input init object
    character(*), parameter		:: F00   = "('(',a,'_run_mct) ',8a)"
    character(*), parameter		:: F04   = "('(',a,'_run_mct) ',2a,2i8,'s')"
    character(*), parameter		:: subName = "(dead_run_mct) "
    !-------------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    ! Reset shr logging to my log file
    !----------------------------------------------------------------------------
    call seq_cdata_setptrs(cdata, ID=COMPID,  infodata=infodata, mpicom=mpicom)
    call shr_file_getLogUnit (shrlogunit)
    call shr_file_getLogLevel(shrloglev)
    call shr_file_setLogUnit (logunits(compid))
    call mpi_comm_rank(mpicom, my_task, ierr)

    lsize = mct_avect_lsize(x2d)
    nflds_d2x = mct_avect_nRattr(d2x)
    nflds_x2d = mct_avect_nRattr(x2d)

    ! PACK

!    do nf=1,nflds_x2d
!       do n=1,lsize
!          ?? = x2d%rAttr(nf,n)
!       enddo
!    enddo

    ! UNPACK
    selectcase(model)
    case('atm')
       gbuf => gbuf_atm
       ncomp = 1
    case('lnd')
       gbuf => gbuf_lnd
       ncomp = 2
    case('ice')
       gbuf => gbuf_ice
       ncomp = 3
    case('ocn')
       gbuf => gbuf_ocn
       ncomp = 4
    case('glc')
       gbuf => gbuf_glc
       ncomp = 5
    case('rof')
       gbuf => gbuf_rof
       ncomp = 6
    case('wav')
       gbuf => gbuf_wav
       ncomp = 7
    end select
    if(model.eq.'rof') then
       do nf=1,nflds_d2x
          do n=1,lsize
             d2x%rAttr(nf,n) = (nf+1) * 1.0_r8
          enddo
       enddo
    else if(model.eq.'glc') then
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
    case('atm')
       ! Set time of next radiadtion computation
       call seq_timemgr_EClockGetData (EClock, next_cday=nextsw_cday)
       call seq_infodata_PutData(infodata, nextsw_cday=nextsw_cday)
    end select

    ! log output for model date

    if (my_task == master_task) then
       call seq_timemgr_EClockGetData(EClock,curr_ymd=CurrentYMD, curr_tod=CurrentTOD)
       write(logunits(compid),F04) model,trim(model),': model date ', CurrentYMD,CurrentTOD
       call shr_sys_flush(logunits(compid))
    end if


    !----------------------------------------------------------------------------
    ! Reset shr logging to original values
    !----------------------------------------------------------------------------
    call shr_file_setLogUnit (shrlogunit)
    call shr_file_setLogLevel(shrloglev)
    call shr_sys_flush(logunits(compid))

  end subroutine dead_run_mct

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: dead_final_mct
  !
  ! !DESCRIPTION:
  !     finalize method for dead model
  !
  ! !REVISION HISTORY:
  !
  ! !INTERFACE: ------------------------------------------------------------------
  !
  subroutine dead_final_mct( model, EClock, cdata, x2d, d2x)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*)            ,intent(in)    :: model
    type(ESMF_Clock)            ,intent(inout) :: EClock
    type(seq_cdata)             ,intent(inout) :: cdata
    type(mct_aVect)             ,intent(inout) :: x2d
    type(mct_aVect)             ,intent(inout) :: d2x

    !EOP

    !--- formats ---
    character(*), parameter :: F00   = "('(',a,'_final_mct) ',8a)"
    character(*), parameter :: F91   = "('(',a,'_final_mct) ',73('-'))"
    character(*), parameter :: subName = "(dead_final_mct) "
    integer(IN)             :: mpicom, my_task,ier, compid
    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------
    call seq_cdata_setptrs(cdata, id=compid, mpicom=mpicom)
    call mpi_comm_rank(mpicom, my_task, ier)

    if (my_task == master_task) then
       write(logunits(compid),F91) model
       write(logunits(compid),F00) model,': end of main integration loop'
       write(logunits(compid),F91) model
    end if

  end subroutine dead_final_mct
  !===============================================================================
end module dead_mct_mod
