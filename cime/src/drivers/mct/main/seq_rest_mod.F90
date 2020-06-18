!===============================================================================
!BOP ===========================================================================
!
! !MODULE: seq_rest_mod -- cpl7 restart reading/writing routines
!
! !DESCRIPTION:
!
!    Reads & writes cpl7 restart files
!
! !REMARKS:
!
!    aVect, domain, and fraction info accessed via seq_avdata_mod
!    to avoid excessively long routine arg lists.
!
! !REVISION HISTORY:
!     2009-Sep-25 - B. Kauffman - move from cpl7 main program into rest module
!     2007-mmm-dd - T. Craig - initial restart functionality
!
! !INTERFACE: ------------------------------------------------------------------

module seq_rest_mod

  ! !USES:

  use shr_kind_mod,      only: R8 => SHR_KIND_R8, IN => SHR_KIND_IN
  use shr_kind_mod,      only: CL => SHR_KIND_CL, CS => SHR_KIND_CS
  use shr_sys_mod,       only: shr_sys_abort, shr_sys_flush
  use shr_mpi_mod,       only: shr_mpi_bcast
  use shr_cal_mod,       only: shr_cal_date2ymd
  use shr_file_mod,      only: shr_file_getunit, shr_file_freeunit
  use mct_mod
  use ESMF
  use component_type_mod

  ! diagnostic routines
  use seq_diag_mct, only : budg_dataG, budg_ns

  ! Sets mpi communicators, logunit and loglevel
  use seq_comm_mct, only: seq_comm_getdata=>seq_comm_setptrs, seq_comm_setnthreads, &
       seq_comm_iamin, CPLID, GLOID, logunit, loglevel

  ! "infodata" gathers various control flags into one datatype
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getData

  ! clock & alarm routines
  use seq_timemgr_mod, only: seq_timemgr_type, seq_timemgr_EClockGetData

  ! diagnostic routines
  use seq_diag_mct, only: budg_datag

  ! lower level io routines
  use seq_io_mod, only: seq_io_read, seq_io_write, seq_io_enddef
  use seq_io_mod, only: seq_io_wopen, seq_io_close

  ! prep modules - coupler communication between different components
  use prep_ocn_mod,    only: prep_ocn_get_x2oacc_ox
  use prep_ocn_mod,    only: prep_ocn_get_x2oacc_ox_cnt
#ifdef SUMMITDEV_PGI
  use prep_ocn_mod,    only: dummy_pgibugfix
#endif
  use prep_rof_mod,    only: prep_rof_get_l2racc_lx
  use prep_rof_mod,    only: prep_rof_get_l2racc_lx_cnt
  use prep_glc_mod,    only: prep_glc_get_l2gacc_lx
  use prep_glc_mod,    only: prep_glc_get_l2gacc_lx_cnt
  use prep_glc_mod,    only: prep_glc_get_x2gacc_gx
  use prep_glc_mod,    only: prep_glc_get_x2gacc_gx_cnt

  use prep_aoflux_mod, only: prep_aoflux_get_xao_ox
  use prep_aoflux_mod, only: prep_aoflux_get_xao_ax

  implicit none

  private

  ! !PUBLIC TYPES:

  ! no public types

  ! !PUBLIC MEMBER FUNCTIONS

  public :: seq_rest_read   ! read  cpl7 restart data
  public :: seq_rest_write  ! write cpl7 restart data

  ! !PUBLIC DATA MEMBERS:

  ! no public data

  !EOP

  !----------------------------------------------------------------------------
  ! local data
  !----------------------------------------------------------------------------

  logical     :: iamin_CPLID            ! pe associated with CPLID
  integer(IN) :: mpicom_GLOID           ! MPI global communicator
  integer(IN) :: mpicom_CPLID           ! MPI cpl communicator

  integer(IN) :: nthreads_GLOID         ! OMP global number of threads
  integer(IN) :: nthreads_CPLID         ! OMP cpl number of threads
  logical     :: drv_threading          ! driver threading control

  logical     :: atm_present            ! .true.  => atm is present
  logical     :: lnd_present            ! .true.  => land is present
  logical     :: ice_present            ! .true.  => ice is present
  logical     :: ocn_present            ! .true.  => ocn is present
  logical     :: rof_present            ! .true.  => land runoff is present
  logical     :: rof_prognostic         ! .true.  => rof comp expects input
  logical     :: glc_present            ! .true.  => glc is present
  logical     :: wav_present            ! .true.  => wav is present
  logical     :: esp_present            ! .true.  => esp is present
  logical     :: iac_present            ! .true.  => iac is present

  logical     :: atm_prognostic         ! .true.  => atm comp expects input
  logical     :: lnd_prognostic         ! .true.  => lnd comp expects input
  logical     :: ice_prognostic         ! .true.  => ice comp expects input
  logical     :: ocn_prognostic         ! .true.  => ocn comp expects input
  logical     :: ocnrof_prognostic      ! .true.  => ocn comp expects runoff input
  logical     :: glc_prognostic         ! .true.  => glc comp expects input
  logical     :: wav_prognostic         ! .true.  => wav comp expects input
  logical     :: esp_prognostic         ! .true.  => esp comp expects input
  logical     :: iac_prognostic         ! .true.  => iac comp expects input

  logical     :: ocn_c2_glcshelf        ! .true.  => ocn to glcshelf coupling on

  !--- temporary pointers ---
  type(mct_gsMap), pointer :: gsmap
  type(mct_aVect), pointer :: x2oacc_ox(:)
  integer        , pointer :: x2oacc_ox_cnt
  type(mct_aVect), pointer :: l2racc_lx(:)
  integer        , pointer :: l2racc_lx_cnt
  type(mct_aVect), pointer :: l2gacc_lx(:)
  integer        , pointer :: l2gacc_lx_cnt
  type(mct_aVect), pointer :: x2gacc_gx(:)
  integer        , pointer :: x2gacc_gx_cnt
  type(mct_aVect), pointer :: xao_ox(:)
  type(mct_aVect), pointer :: xao_ax(:)

  !===============================================================================
contains
  !===============================================================================

  subroutine seq_rest_read(rest_file, infodata, &
       atm, lnd, ice, ocn, rof, glc, wav, esp, iac, &
       fractions_ax, fractions_lx, fractions_ix, fractions_ox, &
       fractions_rx, fractions_gx, fractions_wx, fractions_zx)

    implicit none

    character(*)           , intent(in) :: rest_file  ! restart file path/name
    type(seq_infodata_type), intent(in) :: infodata
    type (component_type) , intent(inout) :: atm(:)
    type (component_type) , intent(inout) :: lnd(:)
    type (component_type) , intent(inout) :: ice(:)
    type (component_type) , intent(inout) :: ocn(:)
    type (component_type) , intent(inout) :: rof(:)
    type (component_type) , intent(inout) :: glc(:)
    type (component_type) , intent(inout) :: wav(:)
    type (component_type) , intent(inout) :: esp(:)
    type (component_type) , intent(inout) :: iac(:)
    type(mct_aVect)  , intent(inout) :: fractions_ax(:)   ! Fractions on atm grid/decomp
    type(mct_aVect)  , intent(inout) :: fractions_lx(:)   ! Fractions on lnd grid/decomp
    type(mct_aVect)  , intent(inout) :: fractions_ix(:)   ! Fractions on ice grid/decomp
    type(mct_aVect)  , intent(inout) :: fractions_ox(:)   ! Fractions on ocn grid/decomp
    type(mct_aVect)  , intent(inout) :: fractions_rx(:)   ! Fractions on rof grid/decomp
    type(mct_aVect)  , intent(inout) :: fractions_gx(:)   ! Fractions on glc grid/decomp
    type(mct_aVect)  , intent(inout) :: fractions_wx(:)   ! Fractions on wav grid/decomp
    type(mct_aVect)  , intent(inout) :: fractions_zx(:)   ! Fractions on iac grid/decomp

    integer(IN)          :: n,n1,n2,n3
    real(r8),allocatable :: ds(:)         ! for reshaping diag data for restart file
    real(r8),allocatable :: ns(:)         ! for reshaping diag data for restart file
    character(len=*), parameter :: subname = "(seq_rest_read) "

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    ! get required infodata
    !----------------------------------------------------------------------------
    iamin_CPLID  = seq_comm_iamin(CPLID)

    call seq_comm_getdata(GLOID,&
         mpicom=mpicom_GLOID, nthreads=nthreads_GLOID)
    call seq_comm_getdata(CPLID, &
         mpicom=mpicom_CPLID, nthreads=nthreads_CPLID)

    call seq_infodata_getData(infodata,      &
         drv_threading=drv_threading,        &
         atm_present=atm_present,        &
         lnd_present=lnd_present,        &
         rof_present=rof_present,        &
         ice_present=ice_present,        &
         ocn_present=ocn_present,        &
         glc_present=glc_present,        &
         wav_present=wav_present,        &
         esp_present=esp_present,        &
         iac_present=iac_present,        &
         atm_prognostic=atm_prognostic,      &
         lnd_prognostic=lnd_prognostic,      &
         ice_prognostic=ice_prognostic,      &
         ocn_prognostic=ocn_prognostic,      &
         rof_prognostic=rof_prognostic,      &
         ocnrof_prognostic=ocnrof_prognostic,    &
         glc_prognostic=glc_prognostic,      &
         wav_prognostic=wav_prognostic,      &
         iac_prognostic=iac_prognostic,      &
         esp_prognostic=esp_prognostic,      &
         ocn_c2_glcshelf=ocn_c2_glcshelf)

    if (iamin_CPLID) then
       if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)
       if (atm_present) then
          gsmap => component_get_gsmap_cx(atm(1))
          xao_ax        => prep_aoflux_get_xao_ax()
          call seq_io_read(rest_file, gsmap, fractions_ax, 'fractions_ax')
          call seq_io_read(rest_file, atm, 'c2x', 'a2x_ax')
          call seq_io_read(rest_file, gsmap, xao_ax, 'xao_ax')
       endif
       if (lnd_present) then
          gsmap => component_get_gsmap_cx(lnd(1))
          call seq_io_read(rest_file, gsmap, fractions_lx, 'fractions_lx')
       endif
       if (lnd_present .and. rof_prognostic) then
          gsmap         => component_get_gsmap_cx(lnd(1))
          l2racc_lx     => prep_rof_get_l2racc_lx()
          l2racc_lx_cnt => prep_rof_get_l2racc_lx_cnt()
          call seq_io_read(rest_file, gsmap, l2racc_lx, 'l2racc_lx')
          call seq_io_read(rest_file, l2racc_lx_cnt ,'l2racc_lx_cnt')
       end if
       if (lnd_present .and. glc_prognostic) then
          gsmap         => component_get_gsmap_cx(lnd(1))
          l2gacc_lx     => prep_glc_get_l2gacc_lx()
          l2gacc_lx_cnt => prep_glc_get_l2gacc_lx_cnt()
          call seq_io_read(rest_file, gsmap, l2gacc_lx, 'l2gacc_lx')
          call seq_io_read(rest_file, l2gacc_lx_cnt ,'l2gacc_lx_cnt')
       end if

       if (ocn_c2_glcshelf) then
          gsmap         => component_get_gsmap_cx(glc(1))
          x2gacc_gx     => prep_glc_get_x2gacc_gx()
          x2gacc_gx_cnt => prep_glc_get_x2gacc_gx_cnt()
          call seq_io_read(rest_file, gsmap, x2gacc_gx, 'x2gacc_gx')
          call seq_io_read(rest_file, x2gacc_gx_cnt ,'x2gacc_gx_cnt')
       end if

       if (ocn_present) then
          gsmap         => component_get_gsmap_cx(ocn(1))
          x2oacc_ox     => prep_ocn_get_x2oacc_ox()
#ifdef SUMMITDEV_PGI
          dummy_pgibugfix = associated(x2oacc_ox)
#endif
          x2oacc_ox_cnt => prep_ocn_get_x2oacc_ox_cnt()
          xao_ox        => prep_aoflux_get_xao_ox()
          call seq_io_read(rest_file, gsmap, fractions_ox, 'fractions_ox')
          call seq_io_read(rest_file, ocn, 'c2x', 'o2x_ox')  ! get o2x_ox
          call seq_io_read(rest_file, gsmap, x2oacc_ox, 'x2oacc_ox')
          call seq_io_read(rest_file, x2oacc_ox_cnt, 'x2oacc_ox_cnt')
          call seq_io_read(rest_file, gsmap, xao_ox, 'xao_ox')
       endif
       if (ice_present) then
          gsmap => component_get_gsmap_cx(ice(1))
          call seq_io_read(rest_file, gsmap, fractions_ix, 'fractions_ix')
          call seq_io_read(rest_file, ice, 'c2x', 'i2x_ix')
       endif
       if (rof_present) then
          gsmap => component_get_gsmap_cx(rof(1))
          call seq_io_read(rest_file, gsmap, fractions_rx, 'fractions_rx')
          call seq_io_read(rest_file, rof, 'c2x', 'r2x_rx')
       endif
       if (glc_present) then
          gsmap => component_get_gsmap_cx(glc(1))
          call seq_io_read(rest_file, gsmap, fractions_gx, 'fractions_gx')
          call seq_io_read(rest_file, glc, 'c2x', 'g2x_gx')
       endif
       if (wav_present) then
          gsmap => component_get_gsmap_cx(wav(1))
          call seq_io_read(rest_file, gsmap, fractions_wx, 'fractions_wx')
          call seq_io_read(rest_file, wav, 'c2x', 'w2x_wx')
       endif
       if (iac_present) then
          gsmap => component_get_gsmap_cx(iac(1))
          call seq_io_read(rest_file, gsmap, fractions_zx, 'fractions_zx')
          call seq_io_read(rest_file, iac, 'c2x', 'z2x_zx')
       endif
       ! Add ESP restart read here

       n = size(budg_dataG)
       allocate(ds(n),ns(n))
       call seq_io_read(rest_file, ds, 'budg_dataG')
       call seq_io_read(rest_file, ns, 'budg_ns')

       n = 0
       do n1 = 1,size(budg_dataG,dim=1)
          do n2 = 1,size(budg_dataG,dim=2)
             do n3 = 1,size(budg_dataG,dim=3)
                n = n + 1
                budg_dataG(n1,n2,n3) = ds(n)
                budg_ns   (n1,n2,n3) = ns(n)
             enddo
          enddo
       enddo
       !     call shr_mpi_bcast(budg_dataG,cpl_io_root) ! not necessary, io lib does bcast

       deallocate(ds,ns)

       if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)

    endif

  end subroutine seq_rest_read

  !===============================================================================

  subroutine seq_rest_write(EClock_d, seq_SyncClock, infodata, &
       atm, lnd, ice, ocn, rof, glc, wav, esp, iac,            &
       fractions_ax, fractions_lx, fractions_ix, fractions_ox, &
       fractions_rx, fractions_gx, fractions_wx, fractions_zx, &
       tag, rest_file)

    implicit none

    type(ESMF_Clock)       , intent(in)    :: EClock_d      ! driver clock
    type(seq_timemgr_type) , intent(inout) :: seq_SyncClock ! contains ptr to driver clock
    type(seq_infodata_type), intent(in)    :: infodata
    type (component_type)       , intent(inout) :: atm(:)
    type (component_type)  , intent(inout) :: lnd(:)
    type (component_type)  , intent(inout) :: ice(:)
    type (component_type)  , intent(inout) :: ocn(:)
    type (component_type)  , intent(inout) :: rof(:)
    type (component_type)  , intent(inout) :: glc(:)
    type (component_type)  , intent(inout) :: wav(:)
    type (component_type)  , intent(inout) :: esp(:)
    type (component_type)  , intent(inout) :: iac(:)
    type(mct_aVect)        , intent(inout) :: fractions_ax(:)   ! Fractions on atm grid/decomp
    type(mct_aVect)        , intent(inout) :: fractions_lx(:)   ! Fractions on lnd grid/decomp
    type(mct_aVect)        , intent(inout) :: fractions_ix(:)   ! Fractions on ice grid/decomp
    type(mct_aVect)        , intent(inout) :: fractions_ox(:)   ! Fractions on ocn grid/decomp
    type(mct_aVect)        , intent(inout) :: fractions_rx(:)   ! Fractions on rof grid/decomp
    type(mct_aVect)        , intent(inout) :: fractions_gx(:)   ! Fractions on glc grid/decomp
    type(mct_aVect)        , intent(inout) :: fractions_wx(:)   ! Fractions on wav grid/decomp
    type(mct_aVect)        , intent(inout) :: fractions_zx(:)   ! Fractions on iac grid/decomp
    character(len=*)       , intent(in)    :: tag
    character(len=CL)      , intent(out)   :: rest_file         ! Restart filename

    integer(IN)   :: n,n1,n2,n3,fk
    integer(IN)   :: curr_ymd         ! Current date YYYYMMDD
    integer(IN)   :: curr_tod         ! Current time-of-day (s)
    integer(IN)   :: yy,mm,dd         ! year, month, day
    character(CL) :: case_name        ! case name
    character(CL) :: cvar             ! char variable
    integer(IN)   :: ivar             ! integer variable
    real(r8)      :: rvar             ! real variable
    logical       :: whead,wdata      ! flags header/data writing
    logical       :: cplroot          ! root pe on cpl id
    integer(IN)   :: iun              ! unit number
    type(mct_gsMap),pointer :: gsmap
    character(len=6) :: year_char

    real(r8),allocatable :: ds(:)     ! for reshaping diag data for restart file
    real(r8),allocatable :: ns(:)     ! for reshaping diag data for restart file
    character(CL) :: model_doi_url
    character(len=*),parameter :: subname = "(seq_rest_write) "

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    !----------------------------------------------------------------------------
    ! get required infodata
    !----------------------------------------------------------------------------
    iamin_CPLID  = seq_comm_iamin(CPLID)

    call seq_comm_getdata(GLOID,&
         mpicom=mpicom_GLOID, nthreads=nthreads_GLOID)

    call seq_comm_getdata(CPLID,&
         mpicom=mpicom_CPLID, nthreads=nthreads_CPLID, iamroot=cplroot)

    call seq_infodata_getData(infodata,      &
         drv_threading=drv_threading,        &
         atm_present=atm_present,        &
         lnd_present=lnd_present,        &
         rof_present=rof_present,        &
         ice_present=ice_present,        &
         ocn_present=ocn_present,        &
         glc_present=glc_present,        &
         wav_present=wav_present,        &
         esp_present=esp_present,        &
         iac_present=iac_present,        &
         atm_prognostic=atm_prognostic,      &
         lnd_prognostic=lnd_prognostic,      &
         ice_prognostic=ice_prognostic,      &
         rof_prognostic=rof_prognostic,      &
         ocn_prognostic=ocn_prognostic,      &
         ocnrof_prognostic=ocnrof_prognostic,    &
         glc_prognostic=glc_prognostic,      &
         wav_prognostic=wav_prognostic,      &
         esp_prognostic=esp_prognostic,      &
         iac_prognostic=iac_prognostic,      &
         ocn_c2_glcshelf=ocn_c2_glcshelf,    &
         case_name=case_name,                &
         model_doi_url=model_doi_url)

    ! Write out infodata and time manager data to restart file

    call seq_timemgr_EClockGetData( EClock_d, curr_ymd=curr_ymd, curr_tod=curr_tod)
    call shr_cal_date2ymd(curr_ymd,yy,mm,dd)
    write(year_char,'(i6.4)') yy
    write(rest_file,"(4a,i2.2,a,i2.2,a,i5.5,a)") &
         trim(case_name), '.cpl'//trim(tag)//'.r.',trim(adjustl(year_char)),'-',mm,'-',dd,'-',curr_tod,'.nc'

    ! Write driver data to restart file

    if (iamin_CPLID) then

       if (drv_threading) call seq_comm_setnthreads(nthreads_CPLID)

       ! copy budg_dataG into 1d array
       n = size(budg_dataG)
       allocate(ds(n),ns(n))
       call shr_mpi_bcast(budg_dataG,mpicom_CPLID) ! pio requires data on all pe's?

       n = 0
       do n1 = 1,size(budg_dataG,dim=1)
          do n2 = 1,size(budg_dataG,dim=2)
             do n3 = 1,size(budg_dataG,dim=3)
                n = n + 1
                ds(n) = budg_dataG(n1,n2,n3)
                ns(n) = budg_ns(n1,n2,n3)
             enddo
          enddo
       enddo

       if (cplroot) then
          iun = shr_file_getUnit()
          call seq_infodata_GetData(infodata,restart_pfile=cvar)
          if (loglevel > 0) write(logunit,"(3A)") subname," write rpointer file ", &
               trim(cvar)
          open(iun, file=cvar, form='FORMATTED')
          write(iun,'(a)') rest_file
          close(iun)
          call shr_file_freeUnit( iun )
       endif

       call shr_mpi_bcast(rest_file,mpicom_CPLID)
       call seq_io_wopen(rest_file,clobber=.true., model_doi_url=model_doi_url)

       ! loop twice (for perf), first time write header, second time write data
       do fk = 1,2
          if (fk == 1) then
             whead = .true.
             wdata = .false.
          elseif (fk == 2) then
             whead = .false.
             wdata = .true.
             call seq_io_enddef(rest_file)
          else
             call shr_sys_abort('driver_write_rstart fk illegal')
          end if
          call seq_infodata_GetData(infodata,nextsw_cday=rvar)
          call seq_io_write(rest_file,rvar,'seq_infodata_nextsw_cday',whead=whead,wdata=wdata)
          call seq_infodata_GetData(infodata,precip_fact=rvar)
          call seq_io_write(rest_file,rvar,'seq_infodata_precip_fact',whead=whead,wdata=wdata)
          call seq_infodata_GetData(infodata,case_name=cvar)
          call seq_io_write(rest_file,trim(cvar),'seq_infodata_case_name',whead=whead,wdata=wdata)

          call seq_timemgr_EClockGetData( EClock_d, start_ymd=ivar)
          call seq_io_write(rest_file,ivar,'seq_timemgr_start_ymd',whead=whead,wdata=wdata)
          call seq_timemgr_EClockGetData( EClock_d, start_tod=ivar)
          call seq_io_write(rest_file,ivar,'seq_timemgr_start_tod',whead=whead,wdata=wdata)
          call seq_timemgr_EClockGetData( EClock_d, ref_ymd=ivar)
          call seq_io_write(rest_file,ivar,'seq_timemgr_ref_ymd'  ,whead=whead,wdata=wdata)
          call seq_timemgr_EClockGetData( EClock_d, ref_tod=ivar)
          call seq_io_write(rest_file,ivar,'seq_timemgr_ref_tod'  ,whead=whead,wdata=wdata)
          call seq_timemgr_EClockGetData( EClock_d, curr_ymd=ivar)
          call seq_io_write(rest_file,ivar,'seq_timemgr_curr_ymd' ,whead=whead,wdata=wdata)
          call seq_timemgr_EClockGetData( EClock_d, curr_tod=ivar)
          call seq_io_write(rest_file,ivar,'seq_timemgr_curr_tod' ,whead=whead,wdata=wdata)

          call seq_io_write(rest_file,ds,'budg_dataG',whead=whead,wdata=wdata)
          call seq_io_write(rest_file,ns,'budg_ns',whead=whead,wdata=wdata)

          if (atm_present) then
             gsmap => component_get_gsmap_cx(atm(1))
             xao_ax        => prep_aoflux_get_xao_ax()
             call seq_io_write(rest_file, gsmap, fractions_ax, 'fractions_ax', &
                  whead=whead, wdata=wdata)
             call seq_io_write(rest_file, atm, 'c2x', 'a2x_ax', &
                  whead=whead, wdata=wdata)
             call seq_io_write(rest_file, gsmap, xao_ax, 'xao_ax', &
                  whead=whead, wdata=wdata)
          endif
          if (lnd_present) then
             gsmap => component_get_gsmap_cx(lnd(1))
             call seq_io_write(rest_file, gsmap, fractions_lx, 'fractions_lx', &
                  whead=whead, wdata=wdata)
          endif
          if (lnd_present .and. rof_prognostic) then
             gsmap         => component_get_gsmap_cx(lnd(1))
             l2racc_lx     => prep_rof_get_l2racc_lx()
             l2racc_lx_cnt =>  prep_rof_get_l2racc_lx_cnt()
             call seq_io_write(rest_file, gsmap, l2racc_lx, 'l2racc_lx', &
                  whead=whead, wdata=wdata)
             call seq_io_write(rest_file, l2racc_lx_cnt, 'l2racc_lx_cnt', &
                  whead=whead, wdata=wdata)
          end if
          if (lnd_present .and. glc_prognostic) then
             gsmap         => component_get_gsmap_cx(lnd(1))
             l2gacc_lx     => prep_glc_get_l2gacc_lx()
             l2gacc_lx_cnt => prep_glc_get_l2gacc_lx_cnt()
             call seq_io_write(rest_file, gsmap, l2gacc_lx, 'l2gacc_lx', &
                  whead=whead, wdata=wdata)
             call seq_io_write(rest_file, l2gacc_lx_cnt, 'l2gacc_lx_cnt', &
                  whead=whead, wdata=wdata)
          end if
          if (ocn_c2_glcshelf) then
             gsmap         => component_get_gsmap_cx(glc(1))
             x2gacc_gx => prep_glc_get_x2gacc_gx()
             x2gacc_gx_cnt => prep_glc_get_x2gacc_gx_cnt()
             call seq_io_write(rest_file, gsmap, x2gacc_gx , 'x2gacc_gx', &
                  whead=whead, wdata=wdata)
             call seq_io_write(rest_file, x2gacc_gx_cnt, 'x2gacc_gx_cnt', &
                  whead=whead, wdata=wdata)
          end if
          if (ocn_present) then
             gsmap         => component_get_gsmap_cx(ocn(1))
             x2oacc_ox     => prep_ocn_get_x2oacc_ox()
#ifdef SUMMITDEV_PGI
             dummy_pgibugfix = associated(x2oacc_ox)
#endif
             x2oacc_ox_cnt => prep_ocn_get_x2oacc_ox_cnt()
             xao_ox        => prep_aoflux_get_xao_ox()
             call seq_io_write(rest_file, gsmap, fractions_ox, 'fractions_ox', &
                  whead=whead, wdata=wdata)
             call seq_io_write(rest_file, ocn, 'c2x', 'o2x_ox', &
                  whead=whead, wdata=wdata)
             call seq_io_write(rest_file, gsmap, x2oacc_ox, 'x2oacc_ox', &
                  whead=whead, wdata=wdata)
             call seq_io_write(rest_file, x2oacc_ox_cnt, 'x2oacc_ox_cnt', &
                  whead=whead, wdata=wdata)
             call seq_io_write(rest_file, gsmap, xao_ox, 'xao_ox', &
                  whead=whead, wdata=wdata)
          endif
          if (ice_present) then
             gsmap  => component_get_gsmap_cx(ice(1))
             call seq_io_write(rest_file, gsmap, fractions_ix, 'fractions_ix', &
                  whead=whead, wdata=wdata)
             call seq_io_write(rest_file, ice, 'c2x', 'i2x_ix', &
                  whead=whead, wdata=wdata)
          endif
          if (rof_present) then
             gsmap  => component_get_gsmap_cx(rof(1))
             call seq_io_write(rest_file, gsmap, fractions_rx, 'fractions_rx', &
                  whead=whead, wdata=wdata)
             call seq_io_write(rest_file, rof, 'c2x', 'r2x_rx', &
                  whead=whead, wdata=wdata)
          endif
          if (glc_present) then
             gsmap  => component_get_gsmap_cx(glc(1))
             call seq_io_write(rest_file, gsmap, fractions_gx, 'fractions_gx', &
                  whead=whead, wdata=wdata)
             call seq_io_write(rest_file, glc, 'c2x', 'g2x_gx', &
                  whead=whead, wdata=wdata)
          endif
          if (wav_present) then
             gsmap  => component_get_gsmap_cx(wav(1))
             call seq_io_write(rest_file, gsmap, fractions_wx, 'fractions_wx', &
                  whead=whead, wdata=wdata)
             call seq_io_write(rest_file, wav, 'c2x', 'w2x_wx', &
                  whead=whead, wdata=wdata)
          endif
          if (iac_present) then
             gsmap  => component_get_gsmap_cx(iac(1))
             call seq_io_write(rest_file, gsmap, fractions_zx, 'fractions_zx', &
                  whead=whead, wdata=wdata)
             call seq_io_write(rest_file, iac, 'c2x', 'z2x_zx', &
                  whead=whead, wdata=wdata)
          endif
          ! Write ESP restart data here
       enddo

       call seq_io_close(rest_file)
       deallocate(ds,ns)

       if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
    endif
  end subroutine seq_rest_write

  !===============================================================================

end module seq_rest_mod
