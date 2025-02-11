#ifdef AIX
@PROCESS ALIAS_SIZE(805306368)
#endif
module docn_comp_mod

  ! !USES:

  use esmf
  use mct_mod
  use perf_mod
  use shr_pcdf_mod
  use shr_const_mod
  use shr_sys_mod
  use shr_kind_mod    , only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_file_mod    , only: shr_file_getunit, shr_file_freeunit
  use shr_mpi_mod     , only: shr_mpi_bcast
  use shr_frz_mod     , only: shr_frz_freezetemp
  use shr_strdata_mod , only: shr_strdata_type, shr_strdata_pioinit, shr_strdata_init
  use shr_strdata_mod , only: shr_strdata_print, shr_strdata_restRead
  use shr_strdata_mod , only: shr_strdata_advance, shr_strdata_restWrite
  use shr_dmodel_mod  , only: shr_dmodel_gsmapcreate, shr_dmodel_rearrGGrid
  use shr_dmodel_mod  , only: shr_dmodel_translate_list, shr_dmodel_translateAV_list, shr_dmodel_translateAV
  use shr_cal_mod     , only: shr_cal_datetod2string
  use seq_timemgr_mod , only: seq_timemgr_EClockGetData

  use docn_shr_mod   , only: datamode       ! namelist input
  use docn_shr_mod   , only: aquap_option   ! derived from datamode namelist input
  use docn_shr_mod   , only: decomp         ! namelist input
  use docn_shr_mod   , only: rest_file      ! namelist input
  use docn_shr_mod   , only: rest_file_strm ! namelist input
  use docn_shr_mod   , only: sst_constant_value ! namelist input
  use docn_shr_mod   , only: RSO_relax_tau      ! namelist input for relaxed slab ocean (RSO)
  use docn_shr_mod   , only: RSO_fixed_MLD      ! namelist input for relaxed slab ocean (RSO)
  use docn_shr_mod   , only: nullstr

#ifdef HAVE_MOAB
!   character(1024)         :: domain_file        ! file containing domain info (set my input)
  use seq_comm_mct,     only: mpoid  ! iMOAB pid for ocean mesh on component pes

  use iso_c_binding
#endif

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: docn_comp_init
  public :: docn_comp_run
  public :: docn_comp_final

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  character(CS) :: myModelName = 'ocn'   ! user defined model name
  logical       :: firstcall = .true.    ! first call logical

  character(len=*),parameter :: rpfile = 'rpointer.ocn'
  integer(IN)   :: dbug = 0              ! debug level (higher is more)

  real(R8),parameter :: cpsw    = shr_const_cpsw        ! specific heat of sea h2o ~ J/kg/K
  real(R8),parameter :: rhosw   = shr_const_rhosw       ! density of sea water ~ kg/m^3
  real(R8),parameter :: TkFrz   = shr_const_TkFrz       ! freezing point, fresh water (Kelvin)
  real(R8),parameter :: TkFrzSw = shr_const_TkFrzSw     ! freezing point, sea   water (Kelvin)
  real(R8),parameter :: latice  = shr_const_latice      ! latent heat of fusion
  real(R8),parameter :: ocnsalt = shr_const_ocn_ref_sal ! ocean reference salinity

  integer(IN)   :: kt,ks,ku,kv,kdhdx,kdhdy,kq,kswp  ! field indices
  integer(IN)   :: kswnet,klwup,klwdn,ksen,klat,kmelth,ksnow,krofi
  integer(IN)   :: kh,kqbot,kfraz
  integer(IN)   :: k10uu           ! index for u10
  integer(IN)   :: kRSO_bckgrd_sst ! index for background SST (relaxed slab ocean)
  integer(IN)   :: index_lat, index_lon
  integer(IN)   :: kmask, kfrac ! frac and mask field indices of docn domain
  integer(IN)   :: ksomask      ! So_omask field index

  type(mct_rearr)        :: rearr
  type(mct_avect)        :: avstrm       ! av of data from stream
  real(R8), pointer      :: somtp(:)
  real(R8), pointer      :: tfreeze(:)
  integer(IN), pointer   :: imask(:)
  real(R8), pointer      :: xc(:), yc(:) ! arryas of model latitudes and longitudes

#ifdef HAVE_MOAB
  integer ::    mdpoid ! data: ocean local component
#endif

  !--------------------------------------------------------------------------
  integer(IN)     , parameter :: ktrans = 8
  character(12)   , parameter :: avifld(1:ktrans) = &
       (/ "t           ","u           ","v           ","dhdx        ",&
          "dhdy        ","s           ","h           ","qbot        "/)
  character(12)   , parameter  :: avofld(1:ktrans) = &
       (/ "So_t        ","So_u        ","So_v        ","So_dhdx     ",&
          "So_dhdy     ","So_s        ","strm_h      ","strm_qbot   "/)
  character(len=*),parameter :: flds_strm = 'strm_h:strm_qbot:So_t'
  !--------------------------------------------------------------------------

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifdef HAVE_MOAB
  SUBROUTINE errorout(ierr, message)
  integer ierr
  character*(*) message
  if (ierr.ne.0) then
    print *, message
    call exit (1)
  end if
  return
  end subroutine
#endif

  !===============================================================================
  subroutine docn_comp_init(Eclock, x2o, o2x, &
       seq_flds_x2o_fields, seq_flds_o2x_fields, &
       SDOCN, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, inst_name, logunit, read_restart, &
       scmMode, scm_multcols, scmlat, scmlon, scm_nx, scm_ny)

    ! !DESCRIPTION: initialize docn model
    use pio        , only : iosystem_desc_t
    use shr_pio_mod, only : shr_pio_getiosys, shr_pio_getiotype
#ifdef HAVE_MOAB
#include "moab/MOABConfig.h"
    use iMOAB, only: iMOAB_DefineTagStorage, iMOAB_GetDoubleTagStorage, &
                     iMOAB_SetIntTagStorage, iMOAB_SetDoubleTagStorage, &
                     iMOAB_ResolveSharedEntities, iMOAB_CreateVertices, &
                     iMOAB_GetMeshInfo, iMOAB_UpdateMeshInfo
#endif
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)       , intent(in)    :: EClock
    type(mct_aVect)        , intent(inout) :: x2o, o2x            ! input/output attribute vectors
    character(len=*)       , intent(in)    :: seq_flds_x2o_fields ! fields from mediator
    character(len=*)       , intent(in)    :: seq_flds_o2x_fields ! fields to mediator
    type(shr_strdata_type) , intent(inout) :: SDOCN               ! model shr_strdata instance (output)
    type(mct_gsMap)        , pointer       :: gsMap               ! model global seg map (output)
    type(mct_gGrid)        , pointer       :: ggrid               ! model ggrid (output)
    integer(IN)            , intent(in)    :: mpicom              ! mpi communicator
    integer(IN)            , intent(in)    :: compid              ! mct comp id
    integer(IN)            , intent(in)    :: my_task             ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task         ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix         ! char string associated with instance
    character(len=*)       , intent(in)    :: inst_name           ! fullname of current instance (ie. "lnd_0001")
    integer(IN)            , intent(in)    :: logunit             ! logging unit number
    logical                , intent(in)    :: read_restart        ! start from restart
    logical                , intent(in)    :: scmMode             ! single column mode
    logical                , intent(in)    :: scm_multcols        ! single column functionality but
                                                                  ! extrapolated over multiple columns
    real(R8)               , intent(in)    :: scmLat              ! single column lat
    real(R8)               , intent(in)    :: scmLon              ! single column lon
    integer(IN)            , intent(in)    :: scm_nx              ! number of points for SCM
                                                                  ! functionality (x direction)
    integer(IN)            , intent(in)    :: scm_ny              ! same but for y direction

    !--- local variables ---
    integer(IN)   :: n,k      ! generic counters
    integer(IN)   :: ierr     ! error code
    integer(IN)   :: lsize    ! local size
    logical       :: exists, exists1   ! file existance
    integer(IN)   :: nu       ! unit number
    character(CL) :: calendar ! model calendar
    integer(IN)   :: currentYMD    ! model date
    integer(IN)   :: currentTOD    ! model sec into model date
    logical       :: write_restart=.false.
    type(iosystem_desc_t), pointer :: ocn_pio_subsystem

#ifdef HAVE_MOAB
    character*400  tagname
    real(R8) latv, lonv
    integer iv, tagindex
    real(R8), allocatable, target :: data(:)
    integer(IN), pointer :: idata(:)   ! temporary
    real(r8), dimension(:), allocatable :: moab_vert_coords  ! temporary
#endif

    !--- formats ---
    character(*), parameter :: F00   = "('(docn_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(docn_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(docn_comp_init) ',a,5i8)"
    character(*), parameter :: F02   = "('(docn_comp_init) ',a,4es13.6)"
    character(*), parameter :: F03   = "('(docn_comp_init) ',a,i8,a)"
    character(*), parameter :: F04   = "('(docn_comp_init) ',2a,2i8,'s')"
    character(*), parameter :: F05   = "('(docn_comp_init) ',a,2f10.4)"
    character(*), parameter :: F06   = "('(docn_comp_init) ',a,5i8)"
    character(*), parameter :: F90   = "('(docn_comp_init) ',73('='))"
    character(*), parameter :: F91   = "('(docn_comp_init) ',73('-'))"
    character(*), parameter :: subName = "(docn_comp_init) "
    !-------------------------------------------------------------------------------

    call t_startf('DOCN_INIT')

    !----------------------------------------------------------------------------
    ! Initialize pio
    !----------------------------------------------------------------------------

    call shr_strdata_pioinit(SDOCN, COMPID)

    !----------------------------------------------------------------------------
    ! Initialize SDOCN
    !----------------------------------------------------------------------------

    call t_startf('docn_strdata_init')

    call seq_timemgr_EClockGetData( EClock, calendar=calendar )

    ! NOTE: shr_strdata_init calls shr_dmodel_readgrid which reads the data model
    ! grid and from that computes SDOCN%gsmap and SDOCN%ggrid. DOCN%gsmap is created
    ! using the decomp '2d1d' (1d decomp of 2d grid)

    if (scmmode) then
       if (my_task == master_task) &
            write(logunit,F05) ' scm lon lat = ',scmlon,scmlat
       call shr_strdata_init(SDOCN,mpicom,compid,name='ocn', &
            scmmode=scmmode,scm_multcols=scm_multcols,scmlon=scmlon,scmlat=scmlat, &
            scm_nx=scm_nx,scm_ny=scm_ny, &
            calendar=calendar, reset_domain_mask=.true.)
    else
       if (datamode == 'SST_AQUAPANAL' .or. datamode == 'SST_AQUAPFILE' .or. &
           datamode == 'SOM_AQUAP' .or. datamode == 'SST_AQUAP_CONSTANT' ) then
          ! Special logic for either prescribed or som aquaplanet - overwrite and
          call shr_strdata_init(SDOCN,mpicom,compid,name='ocn', calendar=calendar, reset_domain_mask=.true.)
       else
          call shr_strdata_init(SDOCN,mpicom,compid,name='ocn', calendar=calendar)
       end if
    endif

    if (my_task == master_task) then
       call shr_strdata_print(SDOCN,'SDOCN data')
    endif

    call t_stopf('docn_strdata_init')

    !----------------------------------------------------------------------------
    ! Initialize data model MCT global seg map, 1d decomp
    !----------------------------------------------------------------------------

    call t_startf('docn_initgsmaps')
    if (my_task == master_task) write(logunit,F00) ' initialize gsmaps'
    call shr_sys_flush(logunit)

    ! create a data model global seqmap (gsmap) given the data model global grid sizes
    ! NOTE: gsmap is initialized using the decomp read in from the docn_in namelist
    ! (which by default is "1d")
    call shr_dmodel_gsmapcreate(gsmap, SDOCN%nxg*SDOCN%nyg, compid, mpicom, decomp)
    lsize = mct_gsmap_lsize(gsmap, mpicom)

    ! create a rearranger from the data model SDOCN%gsmap to gsmap
    call mct_rearr_init(SDOCN%gsmap, gsmap, mpicom, rearr)
    call t_stopf('docn_initgsmaps')

    !----------------------------------------------------------------------------
    ! Initialize data model MCT domain
    !----------------------------------------------------------------------------

    call t_startf('docn_initmctdom')
    if (my_task == master_task) write(logunit,F00) 'copy domains'
    call shr_sys_flush(logunit)

    call shr_dmodel_rearrGGrid(SDOCN%grid, ggrid, gsmap, rearr, mpicom)
    call t_stopf('docn_initmctdom')

    !----------------------------------------------------------------------------
    ! Initialize MCT attribute vectors
    !----------------------------------------------------------------------------

    call t_startf('docn_initmctavs')
    if (my_task == master_task) write(logunit,F00) 'allocate AVs'
    call shr_sys_flush(logunit)

    call mct_aVect_init(o2x, rList=seq_flds_o2x_fields, lsize=lsize)
    call mct_aVect_zero(o2x)

    kt    = mct_aVect_indexRA(o2x,'So_t')
    ks    = mct_aVect_indexRA(o2x,'So_s')
    ku    = mct_aVect_indexRA(o2x,'So_u')
    kv    = mct_aVect_indexRA(o2x,'So_v')
    kdhdx = mct_aVect_indexRA(o2x,'So_dhdx')
    kdhdy = mct_aVect_indexRA(o2x,'So_dhdy')
    kswp  = mct_aVect_indexRA(o2x,'So_fswpen', perrwith='quiet')
    kq    = mct_aVect_indexRA(o2x,'Fioo_q') ! ocn freezing melting potential
    kfraz = mct_aVect_indexRA(o2x,'Fioo_frazil') ! ocn frazil

    call mct_aVect_init(x2o, rList=seq_flds_x2o_fields, lsize=lsize)
    call mct_aVect_zero(x2o)

    kswnet = mct_aVect_indexRA(x2o,'Foxx_swnet')
    klwup  = mct_aVect_indexRA(x2o,'Foxx_lwup')
    ksen   = mct_aVect_indexRA(x2o,'Foxx_sen')
    klat   = mct_aVect_indexRA(x2o,'Foxx_lat')
    krofi  = mct_aVect_indexRA(x2o,'Foxx_rofi')
    klwdn  = mct_aVect_indexRA(x2o,'Faxa_lwdn')
    ksnow  = mct_aVect_indexRA(x2o,'Faxa_snow')
    kmelth = mct_aVect_indexRA(x2o,'Fioi_melth')
    k10uu  = mct_aVect_indexRA(x2o,'So_duu10n')

    call mct_aVect_init(avstrm, rList=flds_strm, lsize=lsize)
    call mct_aVect_zero(avstrm)

    kh    = mct_aVect_indexRA(avstrm,'strm_h')
    kqbot = mct_aVect_indexRA(avstrm,'strm_qbot')
    kRSO_bckgrd_sst = mct_aVect_indexRA(avstrm,'So_t')

    allocate(somtp(lsize))
    allocate(tfreeze(lsize))
    allocate(imask(lsize))
    allocate(xc(lsize))
    allocate(yc(lsize))

    kfrac = mct_aVect_indexRA(ggrid%data,'frac')

    ksomask = mct_aVect_indexRA(o2x,'So_omask', perrwith='quiet')
    if (ksomask /= 0) then
       o2x%rAttr(ksomask, :) = ggrid%data%rAttr(kfrac,:)
    end if

    kmask = mct_aVect_indexRA(ggrid%data,'mask')
    imask(:) = nint(ggrid%data%rAttr(kmask,:))

    index_lon = mct_aVect_indexRA(ggrid%data,'lon')
    xc(:) = ggrid%data%rAttr(index_lon,:)

    index_lat = mct_aVect_indexRA(ggrid%data,'lat')
    yc(:) = ggrid%data%rAttr(index_lat,:)

    call t_stopf('docn_initmctavs')

#ifdef HAVE_MOAB

   allocate(moab_vert_coords(lsize*3))
   do iv = 1, lsize
      lonv = xc(iv) * SHR_CONST_PI/180.
      latv = yc(iv) * SHR_CONST_PI/180.
      moab_vert_coords(3*iv-2)=COS(latv)*COS(lonv)
      moab_vert_coords(3*iv-1)=COS(latv)*SIN(lonv)
      moab_vert_coords(3*iv  )=SIN(latv)
   enddo

   ! create the vertices with coordinates from MCT domain
   ierr = iMOAB_CreateVertices(mpoid, lsize*3, 3, moab_vert_coords)
   if (ierr .ne. 0)  &
      call shr_sys_abort('Error: fail to create MOAB vertices in land model')

   tagname='GLOBAL_ID'//C_NULL_CHAR
   ierr = iMOAB_DefineTagStorage(mpoid, tagname, &
                                 0, & ! dense, integer
                                 1, & ! number of components
                                 tagindex )
   if (ierr .ne. 0)  &
      call shr_sys_abort('Error: fail to retrieve GLOBAL_ID tag ')

   ! get list of global IDs for Dofs
   call mct_gsMap_orderedPoints(gsMap, my_task, idata)

   ierr = iMOAB_SetIntTagStorage ( mpoid, tagname, lsize, &
                                    0, & ! vertex type
                                    idata)
   if (ierr .ne. 0)  &
      call shr_sys_abort('Error: fail to set GLOBAL_ID tag ')

   ierr = iMOAB_ResolveSharedEntities( mpoid, lsize, idata );
   if (ierr .ne. 0)  &
      call shr_sys_abort('Error: fail to resolve shared entities')

   deallocate(moab_vert_coords)
   deallocate(idata)

   ierr = iMOAB_UpdateMeshInfo( mpoid )
   if (ierr .ne. 0)  &
      call shr_sys_abort('Error: fail to update mesh info ')

   allocate(data(lsize))
   ierr = iMOAB_DefineTagStorage( mpoid, "area:aream:frac:mask"//C_NULL_CHAR, &
                                     1, & ! dense, double
                                     1, & ! number of components
                                     tagindex )
   if (ierr > 0 )  &
      call errorout(ierr, 'Error: fail to create tag: area:aream:frac:mask' )

   data(:) = ggrid%data%rAttr(mct_aVect_indexRA(ggrid%data,'area'),:)
   tagname='area'//C_NULL_CHAR
   ierr = iMOAB_SetDoubleTagStorage ( mpoid, tagname, lsize, &
                                      0, & ! set data on vertices
                                      data)
   if (ierr > 0 )  &
      call errorout(ierr, 'Error: fail to get area tag ')

   ! set the same data for aream (model area) as area
   ! data(:) = ggrid%data%rAttr(mct_aVect_indexRA(ggrid%data,'aream'),:)
   tagname='aream'//C_NULL_CHAR
   ierr = iMOAB_SetDoubleTagStorage ( mpoid, tagname, lsize, &
                                      0, & ! set data on vertices
                                      data)
   if (ierr > 0 )  &
      call errorout(ierr, 'Error: fail to set aream tag ')

   data(:) = ggrid%data%rAttr(kmask,:)
   tagname='mask'//C_NULL_CHAR
   ierr = iMOAB_SetDoubleTagStorage ( mpoid, tagname, lsize, &
                                      0, & ! set data on vertices
                                      data)
   if (ierr > 0 )  &
      call errorout(ierr, 'Error: fail to set mask tag ')

   data(:) = ggrid%data%rAttr(kfrac,:)
   tagname='frac'//C_NULL_CHAR
   ierr = iMOAB_SetDoubleTagStorage ( mpoid, tagname, lsize, &
                                      0, & ! set data on vertices
                                      data)
   if (ierr > 0 )  &
      call errorout(ierr, 'Error: fail to set frac tag ')

   deallocate(data)

   ! define tags
   ierr = iMOAB_DefineTagStorage( mpoid, trim(seq_flds_x2o_fields)//C_NULL_CHAR, &
                                     1, & ! dense, double
                                     1, & ! number of components
                                     tagindex )
   if (ierr > 0 )  &
      call errorout(ierr, 'Error: fail to create seq_flds_x2o_fields tags ')

   ierr = iMOAB_DefineTagStorage( mpoid, trim(seq_flds_o2x_fields)//C_NULL_CHAR, &
                                     1, & ! dense, double
                                     1, & ! number of components
                                     tagindex )
   if (ierr > 0 )  &
      call errorout(ierr, 'Error: fail to create seq_flds_o2x_fields tags ')

   ierr = iMOAB_DefineTagStorage( mpoid, trim(flds_strm)//C_NULL_CHAR, &
                                     1, & ! dense, double
                                     1, & ! number of components
                                     tagindex )
   if (ierr > 0 )  &
      call errorout(ierr, 'Error: fail to create flds_strm tags ')
#endif
    !----------------------------------------------------------------------------
    ! Read restart
    !----------------------------------------------------------------------------

    if (read_restart) then
       exists = .false.
       exists1 = .false.
       if (trim(rest_file)      == trim(nullstr) .and. &
           trim(rest_file_strm) == trim(nullstr)) then
          if (my_task == master_task) then
             write(logunit,F00) ' restart filenames from rpointer = ',trim(rpfile)
             call shr_sys_flush(logunit)
             inquire(file=trim(rpfile)//trim(inst_suffix),exist=exists)
             if (exists) then
                nu = shr_file_getUnit()
                open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
                read(nu,'(a)') rest_file
                read(nu,'(a)') rest_file_strm
                close(nu)
                call shr_file_freeUnit(nu)
                inquire(file=trim(rest_file_strm),exist=exists)
                inquire(file=trim(rest_file),exist=exists1)
             endif
          endif
          call shr_mpi_bcast(rest_file,mpicom,'rest_file')
          call shr_mpi_bcast(rest_file_strm,mpicom,'rest_file_strm')
       else
          ! use namelist already read
          if (my_task == master_task) then
             write(logunit,F00) ' restart filenames from namelist '
             call shr_sys_flush(logunit)
             inquire(file=trim(rest_file_strm),exist=exists)
          endif
       endif

       call shr_mpi_bcast(exists,mpicom,'exists')
       call shr_mpi_bcast(exists1,mpicom,'exists1')

       if (     trim(datamode) == 'SOM' &       ! Traditional slab ocean
           .or. trim(datamode) == 'RSO' &       ! Relaxed slab ocean
           .or. trim(datamode) == 'SOM_AQUAP' & ! Aquaplanet slab ocean
          ) then
          if (exists1) then
             if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file)
             call shr_pcdf_readwrite('read',SDOCN%pio_subsystem, SDOCN%io_type, &
                  trim(rest_file), mpicom, gsmap=gsmap, rf1=somtp, rf1n='somtp', io_format=SDOCN%io_format)
          else
             if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file)
          endif
       endif

       if (exists) then
          if (my_task == master_task) write(logunit,F00) ' reading ',trim(rest_file_strm)
          call shr_strdata_restRead(trim(rest_file_strm),SDOCN,mpicom)
       else
          if (my_task == master_task) write(logunit,F00) ' file not found, skipping ',trim(rest_file_strm)
       endif
       call shr_sys_flush(logunit)
    endif

    !----------------------------------------------------------------------------
    ! Set initial ocn state
    !----------------------------------------------------------------------------

    call t_adj_detailf(+2)

    call seq_timemgr_EClockGetData( EClock, curr_ymd=CurrentYMD, curr_tod=CurrentTOD)

    call docn_comp_run(EClock, x2o, o2x, &
         SDOCN, gsmap, ggrid, mpicom, compid, my_task, master_task, &
         inst_suffix, logunit, read_restart, write_restart, &
         currentYMD, currentTOD)

    if (my_task == master_task) write(logunit,F00) 'docn_comp_init done'
    call shr_sys_flush(logunit)

    call t_adj_detailf(-2)

    if (dbug > 0 .and. my_task == master_task) then
       do n = 1,lsize
          write(logunit,F06)'n,ofrac = ',mct_aVect_indexRA(ggrid%data,'frac')
       end do
    end if

    call t_stopf('DOCN_INIT')

  end subroutine docn_comp_init

  subroutine docn_comp_run(EClock, x2o, o2x, &
       SDOCN, gsmap, ggrid, mpicom, compid, my_task, master_task, &
       inst_suffix, logunit, read_restart, write_restart, &
       target_ymd, target_tod, case_name)

    ! !DESCRIPTION:  run method for docn model
#ifdef HAVE_MOAB
    use iMOAB, only: iMOAB_GetMeshInfo, &
                     iMOAB_SetDoubleTagStorage, &
                     iMOAB_WriteMesh
    use seq_flds_mod, only: moab_set_tag_from_av
#endif

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_Clock)       , intent(in)    :: EClock
    type(mct_aVect)        , intent(inout) :: x2o
    type(mct_aVect)        , intent(inout) :: o2x
    type(shr_strdata_type) , intent(inout) :: SDOCN
    type(mct_gsMap)        , pointer       :: gsMap
    type(mct_gGrid)        , pointer       :: ggrid
    integer(IN)            , intent(in)    :: mpicom        ! mpi communicator
    integer(IN)            , intent(in)    :: compid        ! mct comp id
    integer(IN)            , intent(in)    :: my_task       ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task   ! task number of master task
    character(len=*)       , intent(in)    :: inst_suffix   ! char string associated with instance
    integer(IN)            , intent(in)    :: logunit       ! logging unit number
    logical                , intent(in)    :: read_restart  ! start from restart
    logical                , intent(in)    :: write_restart ! restart alarm is on
    integer(IN)            , intent(in)    :: target_ymd    ! model date
    integer(IN)            , intent(in)    :: target_tod    ! model sec into model date
    character(len=*)  , intent(in),optional :: case_name ! case name

    !--- local ---
    integer(IN)   :: yy,mm,dd,tod          ! year month day time-of-day
    integer(IN)   :: n                     ! indices
    integer(IN)   :: nf                    ! fields loop index
    integer(IN)   :: nl                    ! ocn frac index
    integer(IN)   :: lsize                 ! size of attr vect
    integer(IN)   :: idt                   ! integer timestep
    real(R8)      :: dt                    ! timestep
    integer(IN)   :: nu                    ! unit number
    real(R8)      :: hn                    ! h field - mixed layer depth (MLD)
    ! relaxed slab ocean mode variables
    real(R8)      :: RSO_bckgrd_sst        ! background SST 
    real(R8)      :: RSO_X_cool            ! logistics function weight
    real(R8)      :: u10                   ! 10 m wind
    ! relaxed slab ocean fixed parameters
    integer,  parameter :: RSO_slab_option = 0                  ! Option for setting RSO_X_cool
    real(R8), parameter :: RSO_R_cool      = 11.75_r8/86400._r8 ! base cooling rate [K/s]
    real(R8), parameter :: RSO_Tdeep       = 271.0_r8           ! deep water temperature [K]
    real(R8), parameter :: RSO_dT_o        = 27.0_r8            ! scaling temperature gradient
    real(R8), parameter :: RSO_h_o         = 30.0_r8            ! scaling mixed layer depth

    character(len=18) :: date_str
    character(len=CL) :: local_case_name
    real(R8), parameter :: &
         swp = 0.67_R8*(exp((-1._R8*shr_const_zsrflyr) /1.0_R8)) + 0.33_R8*exp((-1._R8*shr_const_zsrflyr)/17.0_R8)

#ifdef HAVE_MOAB
    integer :: ierr     ! error code
    integer :: kgg
    character*100  tagname
    integer tagindex
    real(R8), allocatable, target :: data(:)
#ifdef MOABDEBUG
    integer  :: cur_docn_stepno
    character*100 outfile, wopts, lnum
#endif

#endif

    character(*), parameter :: F00   = "('(docn_comp_run) ',8a)"
    character(*), parameter :: F01   = "('(docn_comp_run) ',a, i7,2x,i5,2x,i5,2x,d21.14)"
    character(*), parameter :: F04   = "('(docn_comp_run) ',2a,2i8,'s')"
    character(*), parameter :: subName = "(docn_comp_run) "
    !-------------------------------------------------------------------------------

    call t_startf('DOCN_RUN')

    call t_startf('docn_run1')
    call seq_timemgr_EClockGetData( EClock, dtime=idt)
    dt = idt * 1.0_R8
    call t_stopf('docn_run1')
    if(present(case_name)) then
       local_case_name = case_name
    else
       local_case_name = ""
    endif

    !--------------------
    ! ADVANCE OCN
    !--------------------

    call t_barrierf('docn_BARRIER',mpicom)
    call t_startf('docn')

    !--- defaults, copy all fields from streams to o2x ---

    lsize = mct_avect_lsize(o2x)
    do n = 1,lsize
       if (ksomask /= 0) then
          o2x%rAttr(ksomask, n) = ggrid%data%rAttr(kfrac,n)
       end if
       o2x%rAttr(kt   ,n) = TkFrz
       o2x%rAttr(ks   ,n) = ocnsalt
       o2x%rAttr(ku   ,n) = 0.0_R8
       o2x%rAttr(kv   ,n) = 0.0_R8
       o2x%rAttr(kdhdx,n) = 0.0_R8
       o2x%rAttr(kdhdy,n) = 0.0_R8
       o2x%rAttr(kq   ,n) = 0.0_R8
! make sure frazil is 0. MPAS-seaice will still use it.
       o2x%rAttr(kfraz,n) = 0.0_R8
       if (kswp /= 0) then
          o2x%rAttr(kswp ,n) = swp
       end if
    enddo

    ! NOTE: for SST_AQUAPANAL, the docn buildnml sets the stream to "null"
    ! and thereby shr_strdata_advance does nothing

    call t_startf('docn_strdata_advance')
    call shr_strdata_advance(SDOCN, target_ymd, target_tod, mpicom, 'docn')
    call t_stopf('docn_strdata_advance')

    !--- copy streams to o2x ---
    call t_barrierf('docn_scatter_BARRIER', mpicom)
    call t_startf('docn_scatter')
    do n = 1, SDOCN%nstreams
       call shr_dmodel_translateAV(SDOCN%avs(n), o2x, avifld, avofld, rearr)
    enddo
    call t_stopf('docn_scatter')

    !-------------------------------------------------
    ! Determine data model behavior based on the mode
    !-------------------------------------------------

    if (my_task .EQ. master_task) &
         write(logunit,*) "DOCN datamode case = ", trim(datamode)
    call t_startf('docn_datamode')
    select case (trim(datamode))

    case('COPYALL')
       ! do nothing extra

    case('SSTDATA')
       lsize = mct_avect_lsize(o2x)
       do n = 1,lsize
          if (ksomask /= 0) then
             o2x%rAttr(ksomask, n) = ggrid%data%rAttr(kfrac,n)
          end if
          o2x%rAttr(kt   ,n) = o2x%rAttr(kt,n) + TkFrz
          o2x%rAttr(ks   ,n) = ocnsalt
          o2x%rAttr(ku   ,n) = 0.0_R8
          o2x%rAttr(kv   ,n) = 0.0_R8
          o2x%rAttr(kdhdx,n) = 0.0_R8
          o2x%rAttr(kdhdy,n) = 0.0_R8
          o2x%rAttr(kq   ,n) = 0.0_R8
          if (kswp /= 0) then
             o2x%rAttr(kswp ,n) = swp
          end if
       enddo

    case('SST_AQUAPANAL')
       lsize = mct_avect_lsize(o2x)
       ! Zero out the attribute vector before calling the prescribed_sst
       ! function - so this also zeroes out the So_omask if it is needed
       ! so need to re-introduce it
       do n = 1,lsize
          o2x%rAttr(:,n) = 0.0_r8
       end do
       call prescribed_sst(xc, yc, lsize, aquap_option, o2x%rAttr(kt,:))
       do n = 1,lsize
          o2x%rAttr(kt,n) = o2x%rAttr(kt,n) + TkFrz
          if (ksomask /= 0) then
             o2x%rAttr(ksomask, n) = ggrid%data%rAttr(kfrac,n)
          end if
       enddo

    case('SST_AQUAPFILE')
       lsize = mct_avect_lsize(o2x)
       do n = 1,lsize
          o2x%rAttr(kt   ,n) = o2x%rAttr(kt,n) + TkFrz
          o2x%rAttr(ks   ,n) = ocnsalt
          o2x%rAttr(ku   ,n) = 0.0_R8
          o2x%rAttr(kv   ,n) = 0.0_R8
          o2x%rAttr(kdhdx,n) = 0.0_R8
          o2x%rAttr(kdhdy,n) = 0.0_R8
          o2x%rAttr(kq   ,n) = 0.0_R8
          if (kswp /= 0) then
             o2x%rAttr(kswp ,n) = swp
          end if
       enddo

    case('SST_AQUAP_CONSTANT')
       lsize = mct_avect_lsize(o2x)
       ! Zero out the attribute vector except for temperature
       do n = 1,lsize
          o2x%rAttr(:,n) = 0.0_r8
       end do
       ! Set temperature and re-set omask
       do n = 1,lsize
          o2x%rAttr(kt,n) = sst_constant_value
          if (ksomask /= 0) then
             o2x%rAttr(ksomask, n) = ggrid%data%rAttr(kfrac,n)
          end if
       enddo

    case('IAF')
       lsize = mct_avect_lsize(o2x)
       do n = 1,lsize
          o2x%rAttr(kt   ,n) = o2x%rAttr(kt,n) + TkFrz
          o2x%rAttr(ks   ,n) = ocnsalt
          o2x%rAttr(ku   ,n) = 0.0_R8
          o2x%rAttr(kv   ,n) = 0.0_R8
          o2x%rAttr(kdhdx,n) = 0.0_R8
          o2x%rAttr(kdhdy,n) = 0.0_R8
          o2x%rAttr(kq   ,n) = 0.0_R8
          if (kswp /= 0) then
             o2x%rAttr(kswp ,n) = swp
          end if
       enddo

    case('SOM')
       lsize = mct_avect_lsize(o2x)
       do n = 1,SDOCN%nstreams
          call shr_dmodel_translateAV(SDOCN%avs(n),avstrm,avifld,avofld,rearr)
       enddo
       if (firstcall) then
          do n = 1,lsize
             if (.not. read_restart) then
                somtp(n) = o2x%rAttr(kt,n) + TkFrz
             endif
             o2x%rAttr(kt,n) = somtp(n)
             o2x%rAttr(kq,n) = 0.0_R8
          enddo
       else   ! firstcall
          tfreeze = shr_frz_freezetemp(o2x%rAttr(ks,:)) + TkFrz
          do n = 1,lsize
             if (imask(n) /= 0) then
                !--- pull out h from av for resuse below ---
                hn = avstrm%rAttr(kh,n)
                !--- compute new temp ---
                o2x%rAttr(kt,n) = somtp(n) + &
                     (x2o%rAttr(kswnet,n) + &  ! shortwave
                      x2o%rAttr(klwup ,n) + &  ! longwave
                      x2o%rAttr(klwdn ,n) + &  ! longwave
                      x2o%rAttr(ksen  ,n) + &  ! sensible
                      x2o%rAttr(klat  ,n) + &  ! latent
                      x2o%rAttr(kmelth,n) - &  ! ice melt
                      avstrm%rAttr(kqbot ,n) - &  ! flux at bottom
                     (x2o%rAttr(ksnow,n)+x2o%rAttr(krofi,n))*latice) * &  ! latent by prec and roff
                     dt/(cpsw*rhosw*hn)
                !--- compute ice formed or melt potential ---
                o2x%rAttr(kq,n) = (tfreeze(n) - o2x%rAttr(kt,n))*(cpsw*rhosw*hn)/dt  ! ice formed q>0
                o2x%rAttr(kt,n) = max(tfreeze(n),o2x%rAttr(kt,n))                    ! reset temp
                somtp(n) = o2x%rAttr(kt,n)                                           ! save temp
             endif
          enddo
       endif   ! firstcall

    ! Relaxed Slab Ocean based on Zarzycki(2016)
    ! Zarzycki, C. M., 2016: Tropical Cyclone Intensity Errors Associated with Lack of Two-Way Ocean Coupling in High-Resolution Global Simulations. J. Climate, 29, 8589–8610.
    ! https://journals.ametsoc.org/view/journals/clim/29/23/jcli-d-16-0273.1.xml
    case('RSO')
      lsize = mct_avect_lsize(o2x)
      do n = 1,SDOCN%nstreams
        call shr_dmodel_translateAV(SDOCN%avs(n),avstrm,avifld,avofld,rearr)
      enddo
      if (firstcall) then
        do n = 1,lsize
          if (.not. read_restart) then
            somtp(n) = o2x%rAttr(kt,n) + TkFrz
          endif
          o2x%rAttr(kt,n) = somtp(n)
          o2x%rAttr(kq,n) = 0.0_R8
        enddo
      else   ! firstcall
        tfreeze = shr_frz_freezetemp(o2x%rAttr(ks,:)) + TkFrz
        do n = 1,lsize
          if (imask(n) /= 0) then
            !*******************************************************************
            if (RSO_fixed_MLD>=0) then
              hn = RSO_fixed_MLD
            else
              hn = avstrm%rAttr(kh,n)
            endif
            ! Get "background" temperature for relaxation
            RSO_bckgrd_sst = avstrm%rAttr(kRSO_bckgrd_sst,n) + TkFrz
            u10 = SQRT(x2o%rAttr(k10uu,n))
            !*******************************************************************
            ! Calculate scaling function - see Eq 3 in Zarzycki (2016)
            if (RSO_slab_option==0) RSO_X_cool = 1._r8/(1._r8+EXP(-0.5_r8*(u10-30._r8)) )                      ! SLAB1
            if (RSO_slab_option==1) RSO_X_cool =(1._r8/(1._r8+EXP(-0.2_r8*(u10-30._r8)) ))*(u10*2.4_r8/80._r8) ! SLAB2
            if (RSO_slab_option==2) RSO_X_cool = 0.0_r8                                                        ! THERMO
            !*******************************************************************
            ! compute new ocean surface temperature
            o2x%rAttr(kt,n) = somtp(n) &
                              ! Thermodynamic terms
                              +( x2o%rAttr(kswnet,n)        & ! shortwave net
                                +x2o%rAttr(klwup ,n)        & ! longwave up
                                +x2o%rAttr(klwdn ,n)        & ! longwave down
                                +x2o%rAttr(ksen  ,n)        & ! sfc sensible heat flux
                                +x2o%rAttr(klat  ,n)        & ! sfc latent heat flux
                                -x2o%rAttr(ksnow ,n)*latice & ! latent heat from snow
                                -x2o%rAttr(krofi ,n)*latice & ! latent heat from runoff
                              ) * dt/(cpsw*rhosw*hn) & 
                              - RSO_X_cool*RSO_R_cool*((somtp(n)-RSO_Tdeep)/RSO_dT_o)*(RSO_h_o/hn)*dt & ! Turb mixing
                              + (1_r8/RSO_relax_tau)*(RSO_bckgrd_sst - somtp(n))*dt ! Newtonian Relaxation
            !*******************************************************************
            ! Ignore ice formed or melt potential
            o2x%rAttr(kq,n) = 0.0
            ! Cap SSTs to freezing
            o2x%rAttr(kt,n) = max( TkFrzSw, o2x%rAttr(kt,n) )
            ! Save temperature to send back to coupler          
            somtp(n) = o2x%rAttr(kt,n)
          endif ! imask /= 0
        enddo ! lsize
      endif ! firstcall

    case('SOM_AQUAP')
       lsize = mct_avect_lsize(o2x)
       do n = 1,SDOCN%nstreams
          call shr_dmodel_translateAV(SDOCN%avs(n),avstrm,avifld,avofld,rearr)
       enddo
       if (firstcall) then
          do n = 1,lsize
             if (.not. read_restart) then
                somtp(n) = o2x%rAttr(kt,n) + TkFrz
             endif
             o2x%rAttr(kt,n) = somtp(n)
             o2x%rAttr(kq,n) = 0.0_R8
          enddo
       else   ! firstcall
          tfreeze = shr_frz_freezetemp(o2x%rAttr(ks,:)) + TkFrz
          do n = 1,lsize
             !--- pull out h from av for resuse below ---
             hn = avstrm%rAttr(kh,n)
             !--- compute new temp ---
             o2x%rAttr(kt,n) = somtp(n) + &
                  (x2o%rAttr(kswnet,n) + &  ! shortwave
                  x2o%rAttr(klwup ,n) + &  ! longwave
                  x2o%rAttr(klwdn ,n) + &  ! longwave
                  x2o%rAttr(ksen  ,n) + &  ! sensible
                  x2o%rAttr(klat  ,n) + &  ! latent
                  x2o%rAttr(kmelth,n) - &  ! ice melt
                  avstrm%rAttr(kqbot ,n) - &  ! flux at bottom
                  (x2o%rAttr(ksnow,n)+x2o%rAttr(krofi,n))*latice) * &  ! latent by prec and roff
                  dt/(cpsw*rhosw*hn)
             !--- compute ice formed or melt potential ---
             o2x%rAttr(kq,n) = (tfreeze(n) - o2x%rAttr(kt,n))*(cpsw*rhosw*hn)/dt  ! ice formed q>0
             somtp(n) = o2x%rAttr(kt,n)                                        ! save temp
          enddo
       endif   ! firstcall

    end select

    call t_stopf('docn_datamode')

#ifdef HAVE_MOAB

   allocate(data(lsize))
   data(:) = 0.0

   ! set dense double tags on vertices of the temporary DOCN app
   ! first set o2x data
   call moab_set_tag_from_av('So_t'//C_NULL_CHAR, o2x, kt, mpoid, data, lsize) 

   call moab_set_tag_from_av('So_s'//C_NULL_CHAR, o2x, ks, mpoid, data, lsize)

   call moab_set_tag_from_av( 'So_u'//C_NULL_CHAR, o2x, ku, mpoid, data, lsize)

   call moab_set_tag_from_av( 'So_v'//C_NULL_CHAR, o2x, kv, mpoid, data, lsize)

   call moab_set_tag_from_av( 'So_dhdx'//C_NULL_CHAR, o2x, kdhdx, mpoid, data, lsize)

   call moab_set_tag_from_av( 'So_dhdy'//C_NULL_CHAR, o2x, kdhdy, mpoid, data, lsize)

   call moab_set_tag_from_av( 'Fioo_q'//C_NULL_CHAR, o2x, kq, mpoid, data, lsize)

   call moab_set_tag_from_av( 'Fioo_frazil'//C_NULL_CHAR, o2x, kq, mpoid, data, lsize)

   if (kswp /= 0) then
      call moab_set_tag_from_av( 'So_fswpen'//C_NULL_CHAR, o2x, kswp, mpoid, data, lsize)
   endif

   ! next set x2o data
   call moab_set_tag_from_av( 'Foxx_swnet'//C_NULL_CHAR, x2o, kswnet, mpoid, data, lsize)

   call moab_set_tag_from_av( 'Foxx_lwup'//C_NULL_CHAR, x2o, klwup, mpoid, data, lsize)

   call moab_set_tag_from_av( 'Foxx_sen'//C_NULL_CHAR, x2o, ksen, mpoid, data, lsize)

   call moab_set_tag_from_av( 'Foxx_lat'//C_NULL_CHAR, x2o, klat, mpoid, data, lsize)

   call moab_set_tag_from_av( 'Foxx_rofi'//C_NULL_CHAR, x2o, krofi, mpoid, data, lsize)

   call moab_set_tag_from_av( 'Faxa_lwdn'//C_NULL_CHAR, x2o, klwdn, mpoid, data, lsize)

   call moab_set_tag_from_av( 'Faxa_snow'//C_NULL_CHAR, x2o, ksnow, mpoid, data, lsize)

   call moab_set_tag_from_av( 'Fioi_melth'//C_NULL_CHAR, x2o, kmelth, mpoid, data, lsize)

   ! next set avstrm data
   call moab_set_tag_from_av( 'strm_h'//C_NULL_CHAR, avstrm, kh, mpoid, data, lsize)

   call moab_set_tag_from_av( 'strm_qbot'//C_NULL_CHAR, avstrm, kqbot, mpoid, data, lsize)


#ifdef MOABDEBUG
    call seq_timemgr_EClockGetData( EClock, stepno=cur_docn_stepno )
    write(lnum,"(I0.2)")cur_docn_stepno
    outfile = 'docn_comp_run_'//trim(lnum)//'.h5m'//C_NULL_CHAR
    wopts   = 'PARALLEL=WRITE_PART'//C_NULL_CHAR
    ierr = iMOAB_WriteMesh(mpoid, outfile, wopts)
    if (ierr > 0 )  then
       write(logunit,*) 'Failed to write ocean component state '
    endif
#endif

#endif

    !----------------------------------------------------------
    ! Debug output
    !----------------------------------------------------------
    if (dbug > 0 .and. my_task == master_task) then
       do n = 1,lsize
          write(logunit,F01)'import: ymd,tod,n,Foxx_swnet = ', target_ymd, target_tod, n, x2o%rattr(kswnet,n)
          write(logunit,F01)'import: ymd,tod,n,Foxx_lwup  = ', target_ymd, target_tod, n, x2o%rattr(klwup,n)
          write(logunit,F01)'import: ymd,tod,n,Foxx_lwdn  = ', target_ymd, target_tod, n, x2o%rattr(klwdn,n)
          write(logunit,F01)'import: ymd,tod,n,Fioi_melth = ', target_ymd, target_tod, n, x2o%rattr(kmelth,n)
          write(logunit,F01)'import: ymd,tod,n,Foxx_sen   = ', target_ymd, target_tod, n, x2o%rattr(ksen,n)
          write(logunit,F01)'import: ymd,tod,n,Foxx_lat   = ', target_ymd, target_tod, n, x2o%rattr(klat,n)
          write(logunit,F01)'import: ymd,tod,n,Foxx_rofi  = ', target_ymd, target_tod, n, x2o%rattr(krofi,n)

          write(logunit,F01)'export: ymd,tod,n,So_t       = ', target_ymd, target_tod, n, o2x%rattr(kt,n)
          write(logunit,F01)'export: ymd,tod,n,So_s       = ', target_ymd, target_tod, n, o2x%rattr(ks,n)
          write(logunit,F01)'export: ymd,tod,n,So_u       = ', target_ymd, target_tod, n, o2x%rattr(ku,n)
          write(logunit,F01)'export: ymd,tod,n,So_v       = ', target_ymd, target_tod, n, o2x%rattr(kv,n)
          write(logunit,F01)'export: ymd,tod,n,So_dhdx    = ', target_ymd, target_tod, n, o2x%rattr(kdhdx,n)
          write(logunit,F01)'export: ymd,tod,n,So_dhdy    = ', target_ymd, target_tod, n, o2x%rattr(kdhdy,n)
          write(logunit,F01)'export: ymd,tod,n,Fioo_q     = ', target_ymd, target_tod, n, o2x%rattr(kq,n)
       end do
    end if

    !--------------------
    ! Write restart
    !--------------------

    if (write_restart) then
       call t_startf('docn_restart')
       call shr_cal_datetod2string(date_str, target_ymd, target_tod)
       write(rest_file,"(6a)") &
            trim(local_case_name), '.docn',trim(inst_suffix),'.r.', &
            trim(date_str),'.nc'
       write(rest_file_strm,"(6a)") &
            trim(local_case_name), '.docn',trim(inst_suffix),'.rs1.', &
            trim(date_str),'.bin'
       if (my_task == master_task) then
          nu = shr_file_getUnit()
          open(nu,file=trim(rpfile)//trim(inst_suffix),form='formatted')
          write(nu,'(a)') rest_file
          write(nu,'(a)') rest_file_strm
          close(nu)
          call shr_file_freeUnit(nu)
       endif
       if (     trim(datamode) == 'SOM' &       ! Traditional slab ocean
           .or. trim(datamode) == 'RSO' &       ! Relaxed slab ocean
           .or. trim(datamode) == 'SOM_AQUAP' & ! Aquaplanet slab ocean
          ) then
          if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file),target_ymd,target_tod
          call shr_pcdf_readwrite('write', SDOCN%pio_subsystem, SDOCN%io_type,&
               trim(rest_file), mpicom, gsmap, clobber=.true., rf1=somtp,rf1n='somtp')
       endif
       if (my_task == master_task) write(logunit,F04) ' writing ',trim(rest_file_strm),target_ymd,target_tod
       call shr_strdata_restWrite(trim(rest_file_strm), SDOCN, mpicom, trim(local_case_name), 'SDOCN strdata')
       call shr_sys_flush(logunit)
       call t_stopf('docn_restart')
    endif

    call t_stopf('docn')

    !----------------------------------------------------------------------------
    ! Log output for model date
    !----------------------------------------------------------------------------

    call t_startf('docn_run2')
    if (my_task == master_task) then
       write(logunit,F04) trim(myModelName),': model date ', target_ymd,target_tod
       call shr_sys_flush(logunit)
    end if

    firstcall = .false.
    call t_stopf('docn_run2')

    call t_stopf('DOCN_RUN')

  end subroutine docn_comp_run

  !===============================================================================
  subroutine docn_comp_final(my_task, master_task, logunit)

    ! !DESCRIPTION:  finalize method for docn model
    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    integer(IN) , intent(in) :: my_task     ! my task in mpi communicator mpicom
    integer(IN) , intent(in) :: master_task ! task number of master task
    integer(IN) , intent(in) :: logunit     ! logging unit number

    !--- formats ---
    character(*), parameter :: F00   = "('(docn_comp_final) ',8a)"
    character(*), parameter :: F91   = "('(docn_comp_final) ',73('-'))"
    character(*), parameter :: subName = "(docn_comp_final) "
    !-------------------------------------------------------------------------------

    call t_startf('DOCN_FINAL')

    if (my_task == master_task) then
       write(logunit,F91)
       write(logunit,F00) trim(myModelName),': end of main integration loop'
       write(logunit,F91)
    end if

    call t_stopf('DOCN_FINAL')

  end subroutine docn_comp_final

  !===============================================================================
  subroutine prescribed_sst(xc, yc, lsize, sst_option, sst)

    real(R8)     , intent(in)    :: xc(:)  !degrees
    real(R8)     , intent(in)    :: yc(:)  !degrees
    integer(IN)  , intent(in)    :: lsize
    integer(IN)  , intent(in)    :: sst_option
    real(R8)     , intent(inout) :: sst(:)

    ! local
    integer  :: i
    real(r8) :: tmp, tmp1, pi
    real(r8) :: rlon(lsize), rlat(lsize)
    real(r8) :: mean_SST, delta_SST

    real(r8), parameter :: pio180 = SHR_CONST_PI/180._r8

    ! Parameters for zonally symmetric experiments
    real(r8), parameter ::   t0_max     = 27._r8
    real(r8), parameter ::   t0_min     = 0._r8
    real(r8), parameter ::   maxlat     = 60._r8*pio180
    real(r8), parameter ::   shift      = 5._r8*pio180
    real(r8), parameter ::   shift9     = 10._r8*pio180
    real(r8), parameter ::   shift10    = 15._r8*pio180

    ! Parameters for zonally asymmetric experiments
    real(r8), parameter ::   t0_max6    = 1._r8
    real(r8), parameter ::   t0_max7    = 3._r8
    real(r8), parameter ::   latcen     = 0._r8*pio180
    real(r8), parameter ::   loncen     = 0._r8*pio180
    real(r8), parameter ::   latrad6    = 15._r8*pio180
    real(r8), parameter ::   latrad8    = 30._r8*pio180
    real(r8), parameter ::   lonrad     = 30._r8*pio180

    pi = SHR_CONST_PI

    ! convert xc and yc from degrees to radians

    rlon(:) = xc(:) * pio180
    rlat(:) = yc(:) * pio180

    ! Control

    if (sst_option < 1 .or. sst_option > 15) then
       call shr_sys_abort ('prescribed_sst: ERROR: sst_option must be between 1 and 15')
    end if

    if (sst_option == 1 .or. sst_option == 6 .or. sst_option == 7 .or. sst_option == 8) then
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if

    ! Flat

    if (sst_option == 2) then
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
             tmp = 1._r8 - tmp*tmp*tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if

    ! Qobs

    if (sst_option == 3) then
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
             tmp = (2._r8 - tmp*tmp*tmp*tmp - tmp*tmp)*0.5_r8
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if

    ! Peaked

    if (sst_option == 4) then
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = (maxlat - abs(rlat(i)))/maxlat
             tmp1 = 1._r8 - tmp
             sst(i) = t0_max*tmp + t0_min*tmp1
          end if
       end do
    end if

    ! Control-5N

    if (sst_option == 5) then
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else if (rlat(i) > shift) then
             tmp = sin((rlat(i)-shift)*pi*0.5_r8/(maxlat-shift))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          else
             tmp = sin((rlat(i)-shift)*pi*0.5_r8/(maxlat+shift))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if

    ! 1KEQ

    if (sst_option == 6) then
       do i = 1,lsize
          if (abs(rlat(i)-latcen) <= latrad6) then
             tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad6)
             tmp1 = tmp1*tmp1
             tmp = abs(rlon(i)-loncen)
             tmp = min(tmp , 2._r8*pi-tmp)
             if(tmp <= lonrad) then
                tmp = cos(tmp*pi*0.5_r8/lonrad)
                tmp = tmp*tmp
                sst(i) = sst(i) + t0_max6*tmp*tmp1
             end if
          end if
       end do
    end if

    ! 3KEQ

    if (sst_option == 7) then
       do i = 1, lsize
          if (abs(rlat(i)-latcen) <= latrad6) then
             tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad6)
             tmp1 = tmp1*tmp1
             tmp = abs(rlon(i)-loncen)
             tmp = min(tmp , 2._r8*pi-tmp)
             if (tmp <= lonrad) then
                tmp = cos(tmp*pi*0.5_r8/lonrad)
                tmp = tmp*tmp
                sst(i) = sst(i) + t0_max7*tmp*tmp1
             end if
          end if
       end do
    end if

    ! 3KW1

    if (sst_option == 8) then
       do i = 1, lsize
          if (abs(rlat(i)-latcen) <= latrad8) then
             tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad8)
             tmp1 = tmp1*tmp1
             tmp = cos(rlon(i)-loncen)
             sst(i) = sst(i) + t0_max7*tmp*tmp1
          end if
       end do
    end if

    ! Control-10N

    if (sst_option == 9) then
       do i = 1, lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else if (rlat(i) > shift9) then
             tmp = sin((rlat(i)-shift9)*pi*0.5_r8/(maxlat-shift9))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          else
             tmp = sin((rlat(i)-shift9)*pi*0.5_r8/(maxlat+shift9))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if

    ! Control-15N

    if (sst_option == 10) then
       do i = 1, lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else if(rlat(i) > shift10) then
             tmp = sin((rlat(i)-shift10)*pi*0.5_r8/(maxlat-shift10))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          else
             tmp = sin((rlat(i)-shift10)*pi*0.5_r8/(maxlat+shift10))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if

    !-------------------------------------------------------------------------------
    ! RCEMIP phase 2 - Mock-Walker
    if (sst_option>=11 .and. sst_option<=15) then
      if (sst_option==11) then; mean_SST = 295 - TkFrz; delta_SST = 1.250; end if ! MW_295dT1p25
      if (sst_option==12) then; mean_SST = 300 - TkFrz; delta_SST = 0.625; end if ! MW_300dT0p625
      if (sst_option==13) then; mean_SST = 300 - TkFrz; delta_SST = 1.250; end if ! MW_300dT1p25
      if (sst_option==14) then; mean_SST = 300 - TkFrz; delta_SST = 2.500; end if ! MW_300dT2p5
      if (sst_option==15) then; mean_SST = 305 - TkFrz; delta_SST = 1.250; end if ! MW_305dT1p25
      do i = 1, lsize
        sst(i) = mean_SST + (delta_SST/2) * cos( rlat(i) * 360/54 )
      end do
    end if
    !-------------------------------------------------------------------------------

  end subroutine prescribed_sst

end module docn_comp_mod
