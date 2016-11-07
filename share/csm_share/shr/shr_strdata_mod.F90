!===============================================================================
! SVN: $Id: shr_strdata.F90 11584 2008-09-08 03:16:24Z mvertens $
! SVN: $HeadURL: https://svn-ccsm-models.cgd.ucar.edu/dlnd7/trunk_tags/dlnd7_090320/shr_strdata.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_strdata_mod -- holds data and methods to advance data models
!
! !DESCRIPTION:
!    holds data and methods to advance data models
!
! !REVISION HISTORY:
!    2009-Apr-15 - T. Craig initial version
!
! !INTERFACE: ------------------------------------------------------------------

module shr_strdata_mod

   use shr_const_mod, only: SHR_CONST_PI
   use shr_kind_mod, only: IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, &
                           CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, &
                           CXX=>SHR_KIND_CXX
   use shr_sys_mod, only : shr_sys_abort, shr_sys_flush
   use shr_mpi_mod, only : shr_mpi_bcast
   use shr_file_mod, only : shr_file_getunit, shr_file_freeunit
   use shr_log_mod, only: loglev  => shr_log_Level
   use shr_log_mod, only: logunit => shr_log_Unit
   use shr_stream_mod  ! stream data type and methods
   use shr_string_mod
   use shr_map_mod
   use shr_cal_mod, only: shr_cal_calendarname, shr_cal_timeSet, &
                          shr_cal_noleap, shr_cal_gregorian, &
                          shr_cal_date2ymd, shr_cal_ymd2date
   use shr_orb_mod, only: shr_orb_decl, shr_orb_cosz, shr_orb_undef_real
   use shr_nl_mod , only: shr_nl_find_group_name
   use shr_tinterp_mod

   use shr_dmodel_mod    ! shr data model stuff

   use shr_mct_mod
   use mct_mod         ! mct
   use perf_mod        ! timing
   use pio             ! pio
   use esmf

   implicit none

   private

! !PUBLIC TYPES:

   public shr_strdata_type

! !PUBLIC MEMBER FUNCTIONS:

   public shr_strdata_readnml
   public shr_strdata_bcastnml
   public shr_strdata_restRead
   public shr_strdata_restWrite
   public shr_strdata_setOrbs
   public shr_strdata_print
   public shr_strdata_init
   public shr_strdata_create
   public shr_strdata_advance
   public shr_strdata_clean
   public shr_strdata_setlogunit
   public shr_strdata_pioinit

! !PUBLIC DATA MEMBERS:


! !PRIVATE:

  integer(SHR_KIND_IN)  :: debug    = 0  ! local debug flag
  integer(IN),parameter :: nStrMax = 30
  integer(IN),parameter :: nVecMax = 30
  character(len=*),public,parameter :: shr_strdata_nullstr = 'null'
  character(len=*),parameter :: shr_strdata_unset = 'NOT_SET'
  real(R8),parameter,private :: dtlimit_default = 1.5_R8

  type shr_strdata_type
     ! --- set by input ---
    character(CL)                  :: dataMode          ! flags physics options wrt input data
    character(CL)                  :: domainFile        ! file   containing domain info
    character(CL)                  :: streams (nStrMax) ! stream description file names
    character(CL)                  :: taxMode (nStrMax) ! time axis cycling mode
    real(R8)                       :: dtlimit (nStrMax) ! dt max/min limit
    character(CL)                  :: vectors (nVecMax) ! define vectors to vector map
    character(CL)                  :: fillalgo(nStrMax) ! fill algorithm
    character(CL)                  :: fillmask(nStrMax) ! fill mask
    character(CL)                  :: fillread(nStrMax) ! fill mapping file to read
    character(CL)                  :: fillwrit(nStrMax) ! fill mapping file to write
    character(CL)                  :: mapalgo (nStrMax) ! scalar map algorithm
    character(CL)                  :: mapmask (nStrMax) ! scalar map mask
    character(CL)                  :: mapread (nStrMax) ! regrid mapping file to read
    character(CL)                  :: mapwrit (nStrMax) ! regrid mapping file to write
    character(CL)                  :: tintalgo(nStrMax) ! time interpolation algorithm
    integer(IN)                    :: io_type

    !--- data required by cosz t-interp method, set by user ---
    real(R8)                       :: eccen
    real(R8)                       :: mvelpp
    real(R8)                       :: lambm0
    real(R8)                       :: obliqr
    integer(IN)                    :: modeldt           ! model dt in seconds

    ! --- internal, public ---
    integer(IN)                    :: nxg               ! model grid lon size
    integer(IN)                    :: nyg               ! model grid lat size
    integer(IN)                    :: nzg               ! model grid vertical size
    integer(IN)                    :: lsize             ! model grid local size
    type(mct_gsmap)                :: gsmap             ! model grid global seg map
    type(mct_ggrid)                :: grid              ! model grid ggrid
    type(mct_avect)                :: avs(nStrMax)      ! model grid stream attribute vectors

    ! --- internal ---
    type(shr_stream_streamType)    :: stream(nStrMax)
    type(iosystem_desc_t), pointer :: pio_subsystem => null()
    type(io_desc_t)                :: pio_iodesc(nStrMax)
    integer(IN)                    :: nstreams          ! number of streams
    integer(IN)                    :: strnxg(nStrMax)
    integer(IN)                    :: strnyg(nStrMax)
    integer(IN)                    :: strnzg(nStrMax)
    logical                        :: dofill(nStrMax)
    logical                        :: domaps(nStrMax)
    integer(IN)                    :: lsizeR(nStrMax)
    type(mct_gsmap)                :: gsmapR(nStrMax)
    type(mct_rearr)                :: rearrR(nStrMax)
    type(mct_ggrid)                :: gridR(nStrMax)
    type(mct_avect)                :: avRLB(nStrMax)    ! Read attrvect
    type(mct_avect)                :: avRUB(nStrMax)    ! Read attrvect
    type(mct_avect)                :: avFUB(nStrMax)    ! Final attrvect
    type(mct_avect)                :: avFLB(nStrMax)    ! Final attrvect
    type(mct_avect)                :: avCoszen(nStrMax) ! data assocaited with coszen time interp
    type(mct_sMatP)                :: sMatPf(nStrMax)
    type(mct_sMatP)                :: sMatPs(nStrMax)
    integer(IN)                    :: ymdLB(nStrMax),todLB(nStrMax)
    integer(IN)                    :: ymdUB(nStrMax),todUB(nStrMax)
    real(R8)                       :: dtmin(nStrMax)
    real(R8)                       :: dtmax(nStrMax)
    integer(IN)                    :: ymd  ,tod
    character(CL)                  :: calendar          ! model calendar for ymd,tod
    integer(IN)                    :: nvectors          ! number of vectors
    integer(IN)                    :: ustrm (nVecMax)
    integer(IN)                    :: vstrm (nVecMax)
    character(CL)                  :: allocstring
  end type shr_strdata_type

  real(R8),parameter,private :: deg2rad = SHR_CONST_PI/180.0_R8
  character(len=*),parameter :: allocstring_value = 'strdata_allocated'

!===============================================================================
  contains
!===============================================================================

  subroutine shr_strdata_init(SDAT,mpicom,compid,name,scmmode,scmlon,scmlat, &
                              gsmap,ggrid,nxg,nyg,nzg,calendar)

  implicit none

  type(shr_strdata_type),intent(inout) :: SDAT
  integer(IN)           ,intent(in)    :: mpicom
  integer(IN)           ,intent(in)    :: compid
  character(len=*)      ,intent(in),optional :: name
  logical               ,intent(in),optional :: scmmode
  real(R8)              ,intent(in),optional :: scmlon
  real(R8)              ,intent(in),optional :: scmlat
  type(mct_gsmap)       ,intent(in),optional :: gsmap
  type(mct_ggrid)       ,intent(in),optional :: ggrid
  integer(IN)           ,intent(in),optional :: nxg
  integer(IN)           ,intent(in),optional :: nyg
  integer(IN)           ,intent(in),optional :: nzg
  character(len=*)      ,intent(in),optional :: calendar

  integer(IN) :: n,m,k         ! generic index
  integer(IN) :: nu,nv         ! u,v index
  integer(IN) :: my_task,npes  ! my task, total pes
  integer(IN),parameter :: master_task = 0
  character(CS) :: lname       ! local name
  character(CL) :: filePath    ! generic file path
  character(CL) :: fileName    ! generic file name
  character(CS) :: timeName    ! domain file: time variable name
  character(CS) ::  lonName    ! domain file: lon  variable name
  character(CS) ::  latName    ! domain file: lat  variable name
  character(CS) ::  hgtName    ! domain file: hgt  variable name
  character(CS) :: maskName    ! domain file: mask variable name
  character(CS) :: areaName    ! domain file: area variable name
  character(CS) :: uname       ! u vector field name
  character(CS) :: vname       ! v vector field name
  character(CXX):: fldList     ! list of fields
  integer(IN)   :: lsize
  integer(IN)   :: nfiles
  integer(IN)   :: ierr
  integer(IN)   :: method
  integer(IN), pointer :: dof(:)
  type(mct_sMat):: sMati
  logical       :: lscmmode

  character(len=*),parameter :: subname = "(shr_strdata_init) "
  character(*),parameter ::   F00 = "('(shr_strdata_init) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  call MPI_COMM_SIZE(mpicom,npes,ierr)
  call MPI_COMM_RANK(mpicom,my_task,ierr)
  !--- Count streams again in case user made changes ---
  if (my_task == master_task) then
     do n=1,nStrMax
        !--- check if a streams string is defined in strdata
        if (trim(SDAT%streams(n)) /= trim(shr_strdata_nullstr)) SDAT%nstreams = max(SDAT%nstreams,n)
        !--- check if a filename is defined in the stream
        call shr_stream_getNFiles(SDAT%stream(n),nfiles)
        if (nfiles > 0) SDAT%nstreams = max(SDAT%nstreams,n)

        if (trim(SDAT%taxMode(n)) == trim(shr_stream_taxis_extend)) SDAT%dtlimit(n) = 1.0e30
     end do
     SDAT%nvectors = 0
     do n=1,nVecMax
        if (trim(SDAT%vectors(n)) /= trim(shr_strdata_nullstr)) SDAT%nvectors = n
     end do
  endif

  call shr_mpi_bcast(SDAT%nstreams  ,mpicom,'nstreams')
  call shr_mpi_bcast(SDAT%nvectors  ,mpicom,'nvectors')
  call shr_mpi_bcast(SDAT%dtlimit   ,mpicom,'dtlimit')

  n = 0
  if (present(gsmap)) then
     n = n + 1
  endif
  if (present(ggrid)) then
     n = n + 1
  endif
  if (present(nxg)) then
     n = n + 1
  endif
  if (present(nyg)) then
     n = n + 1
  endif
  if (present(nzg)) then
     n = n + 1
  endif

  if ( n == 0 .or. n == 5) then
     ! either all set or none set, this is OK
  else
     write(logunit,*) subname,' ERROR: gsmap, ggrid, nxg, nyg and nzg  must be specified together'
     call shr_sys_abort(subname//' ERROR: gsmap, ggrid, nxg, nyg, nzg  set together')
  endif

  lscmmode = .false.
  if (present(scmmode)) then
     lscmmode = scmmode
     if (lscmmode) then
        if (.not.present(scmlon) .or. .not.present(scmlat)) then
           write(logunit,*) subname,' ERROR: scmmode requires scmlon and scmlat'
           call shr_sys_abort(subname//' ERROR: scmmode1 lon lat')
        endif
     endif
  endif

  lname = ""
  if (present(name)) then
     lname = "_"//trim(name)
  endif

  if (present(calendar)) then
     SDAT%calendar = trim(shr_cal_calendarName(trim(calendar)))
  endif

  ! --- initialize streams and stream domains ---

  do n = 1,SDAT%nstreams
     if (my_task == master_task) then
        call shr_stream_getDomainInfo(SDAT%stream(n),filePath,fileName,timeName,lonName, &
             latName,hgtName,maskName,areaName)
        call shr_stream_getFile(filePath,fileName)
        write(logunit,*) subname,' stream ',n
        write(logunit,*) subname,' filePath = ',n,trim(filePath)
        write(logunit,*) subname,' fileName = ',n,trim(fileName)
        write(logunit,*) subname,' timeName = ',n,trim(timeName)
        write(logunit,*) subname,'  lonName = ',n,trim(lonName)
        write(logunit,*) subname,'  latName = ',n,trim(latName)
        write(logunit,*) subname,'  hgtName = ',n,trim(hgtName)
        write(logunit,*) subname,' maskName = ',n,trim(maskName)
        write(logunit,*) subname,' areaName = ',n,trim(areaName)
     endif
     call shr_mpi_bcast(fileName,mpicom)
     call shr_mpi_bcast(lonName,mpicom)
     call shr_mpi_bcast(latName,mpicom)
     call shr_mpi_bcast(hgtName,mpicom)
     call shr_mpi_bcast(maskName,mpicom)
     call shr_mpi_bcast(areaName,mpicom)
     call shr_dmodel_readgrid(SDAT%gridR(n),SDAT%gsmapR(n),SDAT%strnxg(n),SDAT%strnyg(n),SDAT%strnzg(n), &
          fileName, compid, mpicom, '2d1d', lonName, latName, hgtName, maskName, areaName)
     SDAT%lsizeR(n) = mct_gsmap_lsize(SDAT%gsmapR(n),mpicom)
     call mct_gsmap_OrderedPoints(SDAT%gsmapR(n), my_task, dof)
     if (SDAT%strnzg(n) <= 0) then
       call pio_initdecomp(SDAT%pio_subsystem, pio_double, &
            (/SDAT%strnxg(n),SDAT%strnyg(n)/), dof, SDAT%pio_iodesc(n))
     else
       call pio_initdecomp(SDAT%pio_subsystem, pio_double, &
            (/SDAT%strnxg(n),SDAT%strnyg(n),SDAT%strnzg(n)/), dof, SDAT%pio_iodesc(n))
     endif
     deallocate(dof)

     call shr_mpi_bcast(SDAT%stream(n)%calendar,mpicom)
  enddo

  ! --- initialize model domain ---

  if (present(gsmap)) then
     SDAT%nxg = nxg
     SDAT%nyg = nyg
     SDAT%nzg = nzg
     lsize = mct_gsmap_lsize(gsmap,mpicom)
     call mct_gsmap_Copy(gsmap,SDAT%gsmap)
     call mct_ggrid_init(SDAT%grid, ggrid, lsize)
     call mct_aVect_copy(ggrid%data, SDAT%grid%data)
  else
     if (trim(SDAT%domainfile) == trim(shr_strdata_nullstr)) then
        if (SDAT%nstreams > 0) then
           if (my_task == master_task) then
              call shr_stream_getDomainInfo(SDAT%stream(1),filePath,fileName,timeName,lonName, &
                      latName,hgtName,maskName,areaName)
              call shr_stream_getFile(filePath,fileName)
           endif
           call shr_mpi_bcast(fileName,mpicom)
           call shr_mpi_bcast(lonName,mpicom)
           call shr_mpi_bcast(latName,mpicom)
           call shr_mpi_bcast(hgtName,mpicom)
           call shr_mpi_bcast(maskName,mpicom)
           call shr_mpi_bcast(areaName,mpicom)
           if (lscmmode) then
              call shr_dmodel_readgrid(SDAT%grid,SDAT%gsmap,SDAT%nxg,SDAT%nyg,SDAT%nzg, &
                   fileName, compid, mpicom, '2d1d', lonName, latName, hgtName, maskName, areaName, &
                   scmmode=lscmmode,scmlon=scmlon,scmlat=scmlat)
           else
              call shr_dmodel_readgrid(SDAT%grid,SDAT%gsmap,SDAT%nxg,SDAT%nyg,SDAT%nzg, &
                   fileName, compid, mpicom, '2d1d', lonName, latName, hgtName, maskName, areaName)
          endif
        endif
     else
         if (lscmmode) then
            call shr_dmodel_readgrid(SDAT%grid,SDAT%gsmap,SDAT%nxg,SDAT%nyg,SDAT%nzg, &
                 SDAT%domainfile, compid, mpicom, '2d1d', readfrac=.true., &
                 scmmode=lscmmode,scmlon=scmlon,scmlat=scmlat)
         else
            call shr_dmodel_readgrid(SDAT%grid,SDAT%gsmap,SDAT%nxg,SDAT%nyg,SDAT%nzg, &
                 SDAT%domainfile, compid, mpicom, '2d1d', readfrac=.true.)
         endif
     endif
  endif
  SDAT%lsize = mct_gsmap_lsize(SDAT%gsmap,mpicom)

  ! --- setup mapping ---

  do n = 1,SDAT%nstreams
     if (shr_dmodel_gGridCompare(SDAT%gridR(n),SDAT%gsmapR(n),SDAT%grid,SDAT%gsmap, &
         shr_dmodel_gGridCompareMaskSubset,mpicom) .or. trim(SDAT%fillalgo(n))=='none') then
        SDAT%dofill(n) = .false.
     else
        SDAT%dofill(n) = .true.
     endif
     if (trim(SDAT%mapmask(n)) == 'dstmask') then
        method = shr_dmodel_gGridCompareXYabsMask
     else
        method = shr_dmodel_gGridCompareXYabs
     endif
     if (shr_dmodel_gGridCompare(SDAT%gridR(n),SDAT%gsmapR(n),SDAT%grid,SDAT%gsmap, &
         method,mpicom,0.01_r8) .or. trim(SDAT%mapalgo(n))=='none') then
        SDAT%domaps(n) = .false.
     else
        SDAT%domaps(n) = .true.
     endif

     if (SDAT%dofill(n)) then
        if (SDAT%strnzg(n) > 1) then
           write(logunit,*) trim(subname),' do fill called with 3d data, not allowed'
           call shr_sys_abort(subname//': do fill called with 3d data, not allowed')
        endif
        if (trim(SDAT%fillread(n)) == trim(shr_strdata_unset)) then
           if (my_task == master_task) then
             write(logunit,F00) ' calling shr_dmodel_mapSet for fill'
             call shr_sys_flush(logunit)
           endif
           call shr_dmodel_mapSet(SDAT%sMatPf(n), &
             SDAT%gridR(n),SDAT%gsmapR(n),SDAT%strnxg(n),SDAT%strnyg(n), &
             SDAT%gridR(n),SDAT%gsmapR(n),SDAT%strnxg(n),SDAT%strnyg(n), &
             name='mapFill', type='cfill', &
             algo=trim(SDAT%fillalgo(n)),mask=trim(SDAT%fillmask(n)),vect='scalar', &
             compid=compid,mpicom=mpicom)
           if (trim(SDAT%fillwrit(n)) /= trim(shr_strdata_unset)) then
             if (my_task == master_task) then
               write(logunit,F00) ' writing ',trim(SDAT%fillwrit(n))
               call shr_sys_flush(logunit)
             endif
             call shr_mct_sMatWritednc(SDAT%sMatPf(n)%Matrix,SDAT%pio_subsystem,sdat%io_type,  SDAT%fillwrit(n),compid,mpicom)
           endif
        else
           if (my_task == master_task) then
             write(logunit,F00) ' reading ',trim(SDAT%fillread(n))
             call shr_sys_flush(logunit)
           endif
           call shr_mct_sMatReaddnc(sMati,SDAT%gsmapR(n),SDAT%gsmapR(n),'src', &
             filename=trim(SDAT%fillread(n)),mytask=my_task,mpicom=mpicom)
           call mct_sMatP_Init(SDAT%sMatPf(n),sMati,SDAT%gsMapR(n),SDAT%gsmapR(n),0, mpicom, compid)
           call mct_sMat_Clean(sMati)
        endif
     endif

     if (SDAT%domaps(n)) then
        if (SDAT%strnzg(n) > 1) then
           write(logunit,*) trim(subname),' do maps called with 3d data, not allowed'
           call shr_sys_abort(subname//': do maps called with 3d data, not allowed')
        endif
        if (trim(SDAT%mapread(n)) == trim(shr_strdata_unset)) then
           if (my_task == master_task) then
             write(logunit,F00) ' calling shr_dmodel_mapSet for remap'
             call shr_sys_flush(logunit)
           endif
           call shr_dmodel_mapSet(SDAT%sMatPs(n), &
             SDAT%gridR(n),SDAT%gsmapR(n),SDAT%strnxg(n),SDAT%strnyg(n), &
             SDAT%grid ,SDAT%gsmap ,SDAT%nxg      ,SDAT%nyg, &
             name='mapScalar', type='remap', &
             algo=trim(SDAT%mapalgo(n)),mask=trim(SDAT%mapmask(n)), vect='scalar', &
             compid=compid,mpicom=mpicom)
           if (trim(SDAT%mapwrit(n)) /= trim(shr_strdata_unset)) then
             if (my_task == master_task) then
               write(logunit,F00) ' writing ',trim(SDAT%mapwrit(n))
               call shr_sys_flush(logunit)
             endif
             call shr_mct_sMatWritednc(SDAT%sMatPs(n)%Matrix,sdat%pio_subsystem,sdat%io_type,SDAT%mapwrit(n),compid,mpicom)
           endif
        else
           if (my_task == master_task) then
             write(logunit,F00) ' reading ',trim(SDAT%mapread(n))
             call shr_sys_flush(logunit)
           endif
           call shr_mct_sMatReaddnc(sMati,SDAT%gsmapR(n),SDAT%gsmap,'src', &
             filename=trim(SDAT%mapread(n)),mytask=my_task,mpicom=mpicom)
           call mct_sMatP_Init(SDAT%sMatPs(n),sMati,SDAT%gsMapR(n),SDAT%gsmap,0, mpicom, compid)
           call mct_sMat_Clean(sMati)
        endif
     else
        call mct_rearr_init(SDAT%gsmapR(n), SDAT%gsmap, mpicom, SDAT%rearrR(n))
     endif
  enddo

  ! --- setup datatypes ---

  do n = 1,SDAT%nstreams
     if (my_task == master_task) then
        call shr_stream_getModelFieldList(SDAT%stream(n),fldList)
     endif
     call shr_mpi_bcast(fldList,mpicom)
     call mct_aVect_init(SDAT%avs(n)  ,rlist=fldList,lsize=SDAT%lsize)
     call mct_aVect_init(SDAT%avFLB(n),rlist=fldList,lsize=SDAT%lsize)
     call mct_aVect_init(SDAT%avFUB(n),rlist=fldList,lsize=SDAT%lsize)
     call mct_aVect_init(SDAT%avRLB(n),rlist=fldList,lsize=SDAT%lsizeR(n))
     call mct_aVect_init(SDAT%avRUB(n),rlist=fldList,lsize=SDAT%lsizeR(n))
     if (trim(SDAT%tintalgo(n)) == 'coszen') then
        call mct_aVect_init(SDAT%avCoszen(n),rlist="tavCosz",lsize=SDAT%lsize)
     endif
  enddo

  ! --- check vectors and compute ustrm,vstrm ---

  do m = 1,SDAT%nvectors
     if (.not. shr_string_listIsValid(SDAT%vectors(m))) then
        write(logunit,*) trim(subname),' vec fldlist invalid m=',m,trim(SDAT%vectors(m))
        call shr_sys_abort(subname//': vec fldlist invalid:'//trim(SDAT%vectors(m)))
     endif
     if (shr_string_listGetNum(SDAT%vectors(m)) /= 2) then
        write(logunit,*) trim(subname),' vec fldlist ne 2 m=',m,trim(SDAT%vectors(m))
        call shr_sys_abort(subname//': vec fldlist ne 2:'//trim(SDAT%vectors(m)))
     endif
     call shr_string_listGetName(SDAT%vectors(m),1,uname)
     call shr_string_listGetName(SDAT%vectors(m),2,vname)
     nu = 0
     nv = 0
     do n = 1,SDAT%nstreams
        k = mct_aVect_indexRA(SDAT%avRLB(n),trim(uname),perrWith='quiet')
        if (k > 0) nu = n
        k = mct_aVect_indexRA(SDAT%avRLB(n),trim(vname),perrWith='quiet')
        if (k > 0) nv = n
     enddo
     if (nu == 0  .or. nv == 0) then
        write(logunit,*) trim(subname),' vec flds not found  m=',m,trim(SDAT%vectors(m))
        call shr_sys_abort(subname//': vec flds not found:'//trim(SDAT%vectors(m)))
     endif
     if (nu /= nv) then
        if ((.not. shr_dmodel_gGridCompare(SDAT%gridR(nu),SDAT%gsmapR(nu), &
                                         SDAT%gridR(nv),SDAT%gsmapR(nv), &
                                         shr_dmodel_gGridCompareXYabs,mpicom,0.01_r8)) .or. &
            (.not. shr_dmodel_gGridCompare(SDAT%gridR(nu),SDAT%gsmapR(nu), &
                                         SDAT%gridR(nv),SDAT%gsmapR(nv), &
                                         shr_dmodel_gGridCompareMaskZeros,mpicom))) then
           write(logunit,*) trim(subname),' vec fld doms not same m=',m,trim(SDAT%vectors(m))
           call shr_sys_abort(subname//': vec fld doms not same:'//trim(SDAT%vectors(m)))
        endif
     endif
     SDAT%ustrm(m) = nu
     SDAT%vstrm(m) = nv
  enddo

  end subroutine shr_strdata_init

!===============================================================================

  subroutine shr_strdata_advance(SDAT,ymd,tod,mpicom,istr,timers)

  implicit none

  type(shr_strdata_type),intent(inout) :: SDAT
  integer(IN),intent(in) :: ymd    ! current model date
  integer(IN),intent(in) :: tod    ! current model date
  integer(IN),intent(in) :: mpicom
  character(len=*),intent(in),optional :: istr
  logical         ,intent(in),optional :: timers

  integer(IN) :: n,m,i,k,l,kf           ! generic index
  integer(IN) :: my_task,npes
  integer(IN),parameter :: master_task = 0
  logical :: mssrmlf
  logical,allocatable :: newData(:)
  integer(IN) :: ierr
  integer(IN) :: nu,nv
  integer(IN) :: lsize,lsizeR,lsizeF
  integer(IN),allocatable :: ymdmod(:)   ! modified model dates to handle Feb 29
  integer(IN) :: todmod                  ! modified model dates to handle Feb 29
  type(mct_avect) :: avRtmp
  type(mct_avect) :: avRV,avFV
  character(len=32) :: lstr
  logical  :: ltimers
  real(R8) :: flb,fub            ! factor for lb and ub

  !--- for cosz method ---
  real(R8) :: calday             ! julian day of year
  real(R8) :: declin             ! solar declination (radians)
  real(R8) :: eccf               ! earth sun distance factor
  real(R8),pointer :: lonr(:)    ! lon radians
  real(R8),pointer :: latr(:)    ! lat radians
  real(R8),pointer :: cosz(:)    ! cosz
  real(R8),pointer :: tavCosz(:) ! cosz, time avg over [LB,UB]
  real(R8),pointer :: xlon(:),ylon(:)
  real(R8),parameter :: solZenMin    = 0.001_R8 ! minimum solar zenith angle

  type(ESMF_Time) :: timeLB, timeUB  ! lb and ub times
  type(ESMF_TimeInterval) :: timeint ! delta time
  integer(IN)   :: dday          ! delta days
  real(R8)      :: dtime         ! delta time
  integer(IN)   :: uvar,vvar
  logical       :: someNewData ! newData test
  character(CS) :: uname       ! u vector field name
  character(CS) :: vname       ! v vector field name
  integer(IN)   :: year,month,day  ! date year month day
  character(len=*),parameter :: timname = "_strd_adv"
  integer(IN),parameter :: tadj = 2

   !----- formats -----
  character(*),parameter :: subname = "(shr_strdata_advance) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  if (SDAT%nstreams < 1) return

  lstr = ''
  if (present(istr)) then
    lstr = trim(istr)
  endif

  ltimers = .true.
  if (present(timers)) then
    ltimers = timers
  endif

  if (.not.ltimers) call t_adj_detailf(tadj)

  call t_barrierf(trim(lstr)//trim(timname)//'_total_BARRIER',mpicom)
  call t_startf(trim(lstr)//trim(timname)//'_total')

  call MPI_COMM_SIZE(mpicom,npes,ierr)
  call MPI_COMM_RANK(mpicom,my_task,ierr)

  mssrmlf = .false.

  SDAT%ymd = ymd
  SDAT%tod = tod

  if (SDAT%nstreams > 0) then
  allocate(newData(SDAT%nstreams))
  allocate(ymdmod(SDAT%nstreams))

  do n = 1,SDAT%nstreams
     ! ------------------------------------------------------- !
     ! tcraig, Oct 11 2010.  Mismatching calendars: 4 cases    !
     ! ------------------------------------------------------- !
     ! ymdmod and todmod are the ymd and tod to time           !
     ! interpolate to.  Generally, these are just the model    !
     ! date and time.  Also, always use the stream calendar    !
     ! for time interpolation for reasons described below.     !
     ! When there is a calendar mismatch, support Feb 29 in a  !
     ! special way as needed to get reasonable values.         !
     ! Note that when Feb 29 needs to be treated specially,    !
     ! a discontinuity will be introduced.  The size of that   !
     ! discontinuity will depend on the time series input data.!
     ! ------------------------------------------------------- !
     ! (0) The stream calendar and model calendar are          !
     ! identical.  Proceed in the standard way.                !
     ! ------------------------------------------------------- !
     ! (1) If the stream is a no leap calendar and the model   !
     ! is gregorian, then time interpolate on the noleap       !
     ! calendar.  Then if the model date is Feb 29, compute    !
     ! stream data for Feb 28 by setting ymdmod and todmod to  !
     ! Feb 28.  This results in duplicate stream data on       !
     ! Feb 28 and Feb 29 and a discontinuity at the start of   !
     ! Feb 29.                                                 !
     ! This could be potentially updated by using the gregorian!
     ! calendar for time interpolation when the input data is  !
     ! relatively infrequent (say greater than daily) with the !
     ! following concerns.
     !   - The forcing will not be reproduced identically on   !
     !     the same day with climatological inputs data        !
     !   - Input data with variable input frequency might      !
     !     behave funny
     !   - An arbitrary discontinuity will be introduced in    !
     !     the time interpolation method based upon the        !
     !     logic chosen to transition from reproducing Feb 28  !
     !     on Feb 29 and interpolating to Feb 29.              !
     !   - The time gradient of data will change by adding a   !
     !     day arbitrarily.
     ! ------------------------------------------------------- !
     ! (2) If the stream is a gregorian calendar and the model !
     ! is a noleap calendar, then just time interpolate on the !
     ! gregorian calendar.  The causes Feb 29 stream data      !
     ! to be skipped and will lead to a discontinuity at the   !
     ! start of March 1.                                       !
     ! ------------------------------------------------------- !
     ! (3) If the calendars mismatch and neither of the three  !
     ! cases above are recognized, then abort.                 !
     ! ------------------------------------------------------- !

     ! case(0)
     ymdmod(n) = ymd
     todmod    = tod
     if (trim(SDAT%calendar) /= trim(SDAT%stream(n)%calendar)) then
        if ((trim(SDAT%calendar) == trim(shr_cal_gregorian)) .and. &
            (trim(SDAT%stream(n)%calendar) == trim(shr_cal_noleap))) then
            ! case (1), set feb 29 = feb 28
            call shr_cal_date2ymd (ymd,year,month,day)
            if (month == 2 .and. day == 29) then
               call shr_cal_ymd2date(year,2,28,ymdmod(n))
            endif
        else if ((trim(SDAT%calendar) == trim(shr_cal_noleap)) .and. &
            (trim(SDAT%stream(n)%calendar) == trim(shr_cal_gregorian))) then
            ! case (2), feb 29 input data will be skipped automatically
        else
            ! case (3), abort
            write(logunit,*) trim(subname),' ERROR: mismatch calendar ', &
               trim(SDAT%calendar),':',trim(SDAT%stream(n)%calendar)
            call shr_sys_abort(trim(subname)//' ERROR: mismatch calendar ')
        endif
     endif

     call t_barrierf(trim(lstr)//trim(timname)//'_readLBUB_BARRIER',mpicom)
     call t_startf(trim(lstr)//trim(timname)//'_readLBUB')

     call shr_dmodel_readLBUB(SDAT%stream(n),SDAT%pio_subsystem,SDAT%io_type,SDAT%pio_iodesc(n), &
          ymdmod(n),todmod,mpicom,SDAT%gsmapR(n),&
          SDAT%avRLB(n),SDAT%ymdLB(n),SDAT%todLB(n), &
          SDAT%avRUB(n),SDAT%ymdUB(n),SDAT%todUB(n),newData(n), &
          istr=trim(lstr)//'_readLBUB')

     if (debug > 0) then
        write(logunit,*) trim(subname),' newData flag = ',n,newData(n)
        write(logunit,*) trim(subname),' LB ymd,tod = ',n,SDAT%ymdLB(n),SDAT%todLB(n)
        write(logunit,*) trim(subname),' UB ymd,tod = ',n,SDAT%ymdUB(n),SDAT%todUB(n)
     endif

     if (newData(n)) then
        if (debug > 0) then
           write(logunit,*) trim(subname),' newData RLB = ',n,minval(SDAT%avRLB(n)%rAttr), &
	                                 maxval(SDAT%avRLB(n)%rAttr),sum(SDAT%avRLB(n)%rAttr)
           write(logunit,*) trim(subname),' newData RUB = ',n,minval(SDAT%avRUB(n)%rAttr), &
	                                 maxval(SDAT%avRUB(n)%rAttr),sum(SDAT%avRUB(n)%rAttr)
        endif
        call shr_cal_date2ymd(SDAT%ymdLB(n),year,month,day)
        call shr_cal_timeSet(timeLB,SDAT%ymdLB(n),0,SDAT%stream(n)%calendar)
        call shr_cal_timeSet(timeUB,SDAT%ymdUB(n),0,SDAT%stream(n)%calendar)
        timeint = timeUB-timeLB
        call ESMF_TimeIntervalGet(timeint,StartTimeIn=timeLB,d=dday)
        dtime = abs(real(dday,R8) + real(SDAT%todUB(n)-SDAT%todLB(n),R8)/shr_const_cDay)

        SDAT%dtmin(n) = min(SDAT%dtmin(n),dtime)
        SDAT%dtmax(n) = max(SDAT%dtmax(n),dtime)
        if ((SDAT%dtmax(n)/SDAT%dtmin(n)) > SDAT%dtlimit(n)) then
           write(logunit,*) trim(subName),' ERROR: dt limit1 ',SDAT%dtmax(n),SDAT%dtmin(n),SDAT%dtlimit(n)
           write(logunit,*) trim(subName),' ERROR: dt limit2 ',SDAT%ymdLB(n),SDAT%todLB(n), &
                                                               SDAT%ymdUB(n),SDAT%todUB(n)
           call shr_sys_abort(trim(subName)//' ERROR dt limit')
        endif
     endif
     call t_stopf(trim(lstr)//trim(timname)//'_readLBUB')
  enddo

  do n = 1,SDAT%nstreams
     if (newData(n)) then

        if (SDAT%doFill(n)) then
           call t_startf(trim(lstr)//trim(timname)//'_fill')
           lsize = mct_aVect_lsize(SDAT%avRLB(n))
           call mct_aVect_init(avRtmp,SDAT%avRLB(n),lsize)
           call mct_aVect_copy(SDAT%avRLB(n),avRtmp)
           call mct_sMat_avMult(avRtmp,SDAT%sMatPf(n),SDAT%avRLB(n))
           call mct_aVect_copy(SDAT%avRUB(n),avRtmp)
           call mct_sMat_avMult(avRtmp,SDAT%sMatPf(n),SDAT%avRUB(n))
           call mct_aVect_clean(avRtmp)
           call t_stopf(trim(lstr)//trim(timname)//'_fill')
        endif

        if (SDAT%domaps(n)) then
           call t_startf(trim(lstr)//trim(timname)//'_map')
           call mct_sMat_avMult(SDAT%avRLB(n),SDAT%sMatPs(n),SDAT%avFLB(n))
           call mct_sMat_avMult(SDAT%avRUB(n),SDAT%sMatPs(n),SDAT%avFUB(n))
           call t_stopf(trim(lstr)//trim(timname)//'_map')
        else
           call t_startf(trim(lstr)//trim(timname)//'_rearr')
           call mct_rearr_rearrange(SDAT%avRLB(n),SDAT%avFLB(n),SDAT%rearrR(n))
           call mct_rearr_rearrange(SDAT%avRUB(n),SDAT%avFUB(n),SDAT%rearrR(n))
           call t_stopf(trim(lstr)//trim(timname)//'_rearr')
        endif

        if (debug > 0) then
           write(logunit,*) trim(subname),' newData FLB = ',n,minval(SDAT%avFLB(n)%rAttr), &
	                                maxval(SDAT%avFLB(n)%rAttr),sum(SDAT%avFLB(n)%rAttr)
           write(logunit,*) trim(subname),' newData FUB = ',n,minval(SDAT%avFUB(n)%rAttr), &
	                                maxval(SDAT%avFUB(n)%rAttr),sum(SDAT%avFUB(n)%rAttr)
        endif
     endif
  enddo

  do m = 1,SDAT%nvectors
     nu = SDAT%ustrm(m)
     nv = SDAT%vstrm(m)
     if ((SDAT%domaps(nu) .or. SDAT%domaps(nv)) .and. &
         (newdata(nu) .or. newdata(nv))) then

        call t_startf(trim(lstr)//trim(timname)//'_vect')
        call shr_string_listGetName(SDAT%vectors(m),1,uname)
        call shr_string_listGetName(SDAT%vectors(m),2,vname)
        lsizeR = mct_aVect_lsize(SDAT%avRLB(nu))
        lsizeF = mct_aVect_lsize(SDAT%avFLB(nu))
        call mct_aVect_init(avRV,rlist=SDAT%vectors(m),lsize=lsizeR)
        call mct_aVect_init(avFV,rlist=SDAT%vectors(m),lsize=lsizeF)
        allocate(xlon(lsizeR))
        allocate(ylon(lsizeF))
        call mct_aVect_exportRattr(SDAT%gridR(nu)%data,'lon',xlon)
        call mct_aVect_exportRattr(SDAT%grid     %data,'lon',ylon)
        xlon = xlon * deg2rad
        ylon = ylon * deg2rad

        !--- map LB ---

        uvar = mct_aVect_indexRA(SDAT%avRLB(nu),trim(uname))
        vvar = mct_aVect_indexRA(SDAT%avRLB(nv),trim(vname))
        do i = 1,lsizeR
           avRV%rAttr(1,i) =  SDAT%avRLB(nu)%rAttr(uvar,i) * cos(xlon(i))  &
                             -SDAT%avRLB(nv)%rAttr(vvar,i) * sin(xlon(i))
           avRV%rAttr(2,i) =  SDAT%avRLB(nu)%rAttr(uvar,i) * sin(xlon(i))  &
                             +SDAT%avRLB(nv)%rAttr(vvar,i) * cos(xlon(i))
        enddo
        call mct_sMat_avMult(avRV,SDAT%sMatPs(nu),avFV)
! ---   don't need to recompute uvar and vvar, should be the same
!       uvar = mct_aVect_indexRA(SDAT%avFLB(nu),trim(uname))
!       vvar = mct_aVect_indexRA(SDAT%avFLB(nv),trim(vname))
        do i = 1,lsizeF
           SDAT%avFLB(nu)%rAttr(uvar,i) =  avFV%rAttr(1,i) * cos(ylon(i))  &
                                          +avFV%rAttr(2,i) * sin(ylon(i))
           SDAT%avFLB(nv)%rAttr(vvar,i) = -avFV%rAttr(1,i) * sin(ylon(i))  &
                                          +avFV%rAttr(2,i) * cos(ylon(i))
        enddo

        !--- map UB ---

        uvar = mct_aVect_indexRA(SDAT%avRUB(nu),trim(uname))
        vvar = mct_aVect_indexRA(SDAT%avRUB(nv),trim(vname))
        do i = 1,lsizeR
           avRV%rAttr(1,i) =  SDAT%avRUB(nu)%rAttr(uvar,i) * cos(xlon(i))  &
                             -SDAT%avRUB(nv)%rAttr(vvar,i) * sin(xlon(i))
           avRV%rAttr(2,i) =  SDAT%avRUB(nu)%rAttr(uvar,i) * sin(xlon(i))  &
                             +SDAT%avRUB(nv)%rAttr(vvar,i) * cos(xlon(i))
        enddo
        call mct_sMat_avMult(avRV,SDAT%sMatPs(nu),avFV)
! ---   don't need to recompute uvar and vvar, should be the same
!       uvar = mct_aVect_indexRA(SDAT%avFUB(nu),trim(uname))
!       vvar = mct_aVect_indexRA(SDAT%avFUB(nv),trim(vname))
        do i = 1,lsizeF
           SDAT%avFUB(nu)%rAttr(uvar,i) =  avFV%rAttr(1,i) * cos(ylon(i))  &
                                          +avFV%rAttr(2,i) * sin(ylon(i))
           SDAT%avFUB(nv)%rAttr(vvar,i) = -avFV%rAttr(1,i) * sin(ylon(i))  &
                                          +avFV%rAttr(2,i) * cos(ylon(i))
        enddo

        call mct_aVect_clean(avRV)
        call mct_aVect_clean(avFV)
        deallocate(xlon,ylon)

        call t_stopf(trim(lstr)//trim(timname)//'_vect')
     endif
  enddo

  do n = 1,SDAT%nstreams

     !--- method: coszen -------------------------------------------------------
     if (trim(SDAT%tintalgo(n)) == 'coszen') then
        call t_startf(trim(lstr)//trim(timname)//'_coszen')

        !--- make sure orb info has been set ---
        if (SDAT%eccen == SHR_ORB_UNDEF_REAL) then
           call shr_sys_abort(subname//' ERROR in orb params for coszen tinterp')
        else if (SDAT%modeldt < 1) then
           call shr_sys_abort(subname//' ERROR: model dt < 1 for coszen tinterp')
        endif

        !--- allocate avg cosz array ---
        lsizeF = mct_aVect_lsize(SDAT%avFLB(n))
        allocate(tavCosz(lsizeF),cosz(lsizeF),lonr(lsizeF),latr(lsizeF))

        !--- get lat/lon data ---
        kf = mct_aVect_indexRA(SDAT%grid%data,'lat')
        latr(:) = SDAT%grid%data%rAttr(kf,:) * deg2rad
        kf = mct_aVect_indexRA(SDAT%grid%data,'lon')
        lonr(:) = SDAT%grid%data%rAttr(kf,:) * deg2rad

        call t_startf(trim(lstr)//trim(timname)//'_coszenC')
        cosz = 0.0_r8
        call shr_tInterp_getCosz(cosz,lonr,latr,ymdmod(n),todmod, &
             SDAT%eccen,SDAT%mvelpp,SDAT%lambm0,SDAT%obliqr,SDAT%stream(n)%calendar)
        call t_stopf(trim(lstr)//trim(timname)//'_coszenC')

        if (newdata(n)) then
           !--- compute a new avg cosz ---
           call t_startf(trim(lstr)//trim(timname)//'_coszenN')
           call shr_tInterp_getAvgCosz(tavCosz,lonr,latr, &
                SDAT%ymdLB(n),SDAT%todLB(n), SDAT%ymdUB(n),SDAT%todUB(n), &
                SDAT%eccen,SDAT%mvelpp,SDAT%lambm0,SDAT%obliqr,SDAT%modeldt,&
	        SDAT%stream(n)%calendar)
           call mct_avect_importRAttr(SDAT%avCoszen(n),'tavCosz',tavCosz,lsizeF)
           call t_stopf(trim(lstr)//trim(timname)//'_coszenN')
        else
           !--- reuse existing avg cosz ---
           call mct_avect_exportRAttr(SDAT%avCoszen(n),'tavCosz',tavCosz)
        endif

        !--- t-interp is LB data normalized with this factor: cosz/tavCosz ---
        do i = 1,lsizeF
           if (cosz(i) > solZenMin) then
              SDAT%avs(n)%rAttr(:,i) = SDAT%avFLB(n)%rAttr(:,i)*cosz(i)/tavCosz(i)
           else
              SDAT%avs(n)%rAttr(:,i) =  0._r8
           endif
        enddo
        deallocate(tavCosz,cosz,lonr,latr)
        call t_stopf(trim(lstr)//trim(timname)//'_coszen')

     !--- method: not coszen ---------------------------------------------------
     elseif (trim(SDAT%tintalgo(n)) /= trim(shr_strdata_nullstr)) then

        call t_startf(trim(lstr)//trim(timname)//'_tint')
        call shr_tInterp_getFactors(SDAT%ymdlb(n),SDAT%todlb(n),SDAT%ymdub(n),SDAT%todub(n), &
                                    ymdmod(n),todmod,flb,fub, &
                                    calendar=SDAT%stream(n)%calendar,algo=trim(SDAT%tintalgo(n)))
        if (debug > 0) then
           write(logunit,*) trim(subname),' interp = ',n,flb,fub
        endif
        SDAT%avs(n)%rAttr(:,:) = SDAT%avFLB(n)%rAttr(:,:)*flb + SDAT%avFUB(n)%rAttr(:,:)*fub
        call t_stopf(trim(lstr)//trim(timname)//'_tint')

     else
        call t_startf(trim(lstr)//trim(timname)//'_zero')
        call mct_avect_zero(SDAT%avs(n))
        call t_stopf(trim(lstr)//trim(timname)//'_zero')
     endif
     if (debug > 0) then
        write(logunit,*) trim(subname),' SDAT av = ',n,minval(SDAT%avs(n)%rAttr),maxval(SDAT%avs(n)%rAttr),sum(SDAT%avs(n)%rAttr)
     endif

  enddo

  deallocate(newData)
  deallocate(ymdmod)

  endif    ! nstreams > 0

  call t_stopf(trim(lstr)//trim(timname)//'_total')
  if (.not.ltimers) call t_adj_detailf(-tadj)

  end subroutine shr_strdata_advance

!===============================================================================
  subroutine shr_strdata_clean(SDAT)

  implicit none

  type(shr_strdata_type),intent(inout) :: SDAT

  integer(IN) :: n

   !----- formats -----
  character(len=*),parameter :: subname = "(shr_strdata_clean) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  if (SDAT%nxg * SDAT%nyg == 0) then
     return
  endif

  ! Free MCT and PIO data first, while we still know which objects were
  ! allocated for which streams.
  call mct_ggrid_clean(SDAT%grid)
  call mct_gsmap_clean(SDAT%gsmap)

  do n = 1, SDAT%nstreams
     call pio_freedecomp(SDAT%pio_subsystem, SDAT%pio_iodesc(n))
     call mct_avect_clean(SDAT%avs(n))
     call mct_avect_clean(SDAT%avRLB(n))
     call mct_avect_clean(SDAT%avRUB(n))
     call mct_avect_clean(SDAT%avFLB(n))
     call mct_avect_clean(SDAT%avFUB(n))
     call mct_ggrid_clean(SDAT%gridR(n))
     if (SDAT%dofill(n)) call mct_sMatP_clean(SDAT%sMatPf(n))
     if (SDAT%domaps(n)) call mct_sMatP_clean(SDAT%sMatPs(n))
     call mct_gsmap_clean(SDAT%gsmapR(n))
  end do

  ! Now that all sub-objects are freed, clear components of the strdata
  ! object itself.
  SDAT%nxg = 0
  SDAT%nyg = 0
  SDAT%nzg = 0
  SDAT%strnxg = 0
  SDAT%strnyg = 0
  SDAT%strnzg = 0

  SDAT%nstreams = 0
  SDAT%nvectors = 0
  SDAT%ustrm = 0
  SDAT%vstrm = 0

  SDAT%dofill = .false.
  SDAT%domaps = .false.

  end subroutine shr_strdata_clean

!===============================================================================
!===============================================================================
  subroutine shr_strdata_restWrite(filename,SDAT,mpicom,str1,str2)

  implicit none

  character(len=*)      ,intent(in)    :: filename
  type(shr_strdata_type),intent(inout) :: SDAT
  integer(IN)           ,intent(in)    :: mpicom
  character(len=*)      ,intent(in)    :: str1
  character(len=*)      ,intent(in)    :: str2

  !--- local ----
  type(shr_stream_streamtype),pointer :: streams(:)
  integer(IN) :: n,my_task,ier

   !----- formats -----
  character(len=*),parameter :: subname = "(shr_strdata_restWrite) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  call MPI_COMM_RANK(mpicom,my_task,ier)

  if (my_task == 0) then
     call shr_stream_restWrite(SDAT%stream,trim(filename),trim(str1),trim(str2),SDAT%nstreams)
  endif

  end subroutine shr_strdata_restWrite

!===============================================================================
!===============================================================================
  subroutine shr_strdata_restRead(filename,SDAT,mpicom)

  implicit none

  character(len=*)      ,intent(in)    :: filename
  type(shr_strdata_type),intent(inout) :: SDAT
  integer(IN)           ,intent(in)    :: mpicom

  !--- local ----
  type(shr_stream_streamtype),pointer :: streams(:)
  integer(IN) :: n,my_task,ier

   !----- formats -----
  character(len=*),parameter :: subname = "(shr_strdata_restRead) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  call MPI_COMM_RANK(mpicom,my_task,ier)

  if (my_task == 0) then
     call shr_stream_restRead(SDAT%stream,trim(filename),SDAT%nstreams)
  endif

  end subroutine shr_strdata_restRead

!===============================================================================
!===============================================================================
  subroutine shr_strdata_setOrbs(SDAT,eccen,mvelpp,lambm0,obliqr,modeldt)

  implicit none

  type(shr_strdata_type),intent(inout) :: SDAT
  real(R8),intent(in) :: eccen
  real(R8),intent(in) :: mvelpp
  real(R8),intent(in) :: lambm0
  real(R8),intent(in) :: obliqr
  integer(IN),intent(in) :: modeldt

  !----- formats -----
  character(len=*),parameter :: subname = "(shr_strdata_setOrbs) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

  SDAT%eccen   = eccen
  SDAT%mvelpp  = mvelpp
  SDAT%lambm0  = lambm0
  SDAT%obliqr  = obliqr
  SDAT%modeldt = modeldt

  end subroutine shr_strdata_setOrbs

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_strdata_readnml -- read control strdata
!
! !DESCRIPTION:
!     Reads strdata common to all data models
!
! !REVISION HISTORY:
!     2004-Dec-15 - J. Schramm - first version
!     2009-Apr-16 - T. Craig - add minimal parallel support
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_strdata_readnml(SDAT,file,rc,mpicom)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_strdata_type),intent(inout):: SDAT  ! strdata data data-type
   character(*),optional ,intent(in)   :: file ! file to read strdata from
   integer(IN),optional  ,intent(out)  :: rc   ! return code
   integer(IN),optional  ,intent(in)   :: mpicom ! mpi comm

!EOP

   integer(IN)    :: rCode         ! return code
   integer(IN)    :: nUnit         ! fortran i/o unit number
   integer(IN)    :: n             ! generic loop index
   integer(IN)    :: my_task       ! my task number, 0 is default
   integer(IN)    :: master_task   ! master task number, 0 is default
   integer(IN)    :: ntasks        ! total number of tasks

   !----- temporary/local namelist vars to read int -----
   character(CL) :: dataMode          ! flags physics options wrt input data
   character(CL) :: domainFile        ! file   containing domain info
   character(CL) :: streams(nStrMax)  ! stream description file names
   character(CL) :: taxMode(nStrMax)  ! time axis cycling mode
   real(R8)      :: dtlimit(nStrMax)  ! delta time limiter
   character(CL) :: vectors(nVecMax)  ! define vectors to vector map
   character(CL) :: fillalgo(nStrMax) ! fill algorithm
   character(CL) :: fillmask(nStrMax) ! fill mask
   character(CL) :: fillread(nStrMax) ! fill mapping file to read
   character(CL) :: fillwrite(nStrMax)! fill mapping file to write
   character(CL) :: mapalgo(nStrMax)  ! scalar map algorithm
   character(CL) :: mapmask(nStrMax)  ! scalar map mask
   character(CL) :: mapread(nStrMax)  ! regrid mapping file to read
   character(CL) :: mapwrite(nStrMax) ! regrid mapping file to write
   character(CL) :: tintalgo(nStrMax) ! time interpolation algorithm
   character(CL) :: io_type
   integer(IN)   :: num_iotasks
   integer(IN)   :: io_root
   integer(IN)   :: io_stride
   integer(IN)   :: num_agg
   character(CL) :: fileName    ! generic file name
   integer(IN)   :: yearFirst   ! first year to use in data stream
   integer(IN)   :: yearLast    ! last  year to use in data stream
   integer(IN)   :: yearAlign   ! data year that aligns with yearFirst

   !----- define namelist -----
   namelist / shr_strdata_nml / &
        dataMode        &
      , domainFile      &
      , streams         &
      , taxMode         &
      , dtlimit         &
      , vectors         &
      , fillalgo        &
      , fillmask        &
      , fillread        &
      , fillwrite       &
      , mapalgo         &
      , mapmask         &
      , mapread         &
      , mapwrite        &
      , tintalgo
   !----- formats -----
   character(*),parameter :: subName = "(shr_strdata_readnml) "
   character(*),parameter ::   F00 = "('(shr_strdata_readnml) ',8a)"
   character(*),parameter ::   F01 = "('(shr_strdata_readnml) ',a,i6,a)"
   character(*),parameter ::   F02 = "('(shr_strdata_readnml) ',a,es13.6)"
   character(*),parameter ::   F03 = "('(shr_strdata_readnml) ',a,l6)"
   character(*),parameter ::   F04 = "('(shr_strdata_readnml) ',a,i2,a,a)"
   character(*),parameter ::   F20 = "('(shr_strdata_readnml) ',a,i6,a)"
   character(*),parameter ::   F90 = "('(shr_strdata_readnml) ',58('-'))"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (present(rc)) rc = 0

   my_task = 0
   master_task = 0
   ntasks = 1
   if (present(mpicom)) then
      call MPI_COMM_RANK(mpicom,my_task,rCode)
      call MPI_COMM_SIZE(mpicom,ntasks,rCode)
   endif

!--master--task--
   if (my_task == master_task) then

   !----------------------------------------------------------------------------
   ! set default values for namelist vars
   !----------------------------------------------------------------------------
   dataMode    = 'NULL'
   domainFile  = trim(shr_strdata_nullstr)
   streams(:)  = trim(shr_strdata_nullstr)
   taxMode(:)  = trim(shr_stream_taxis_cycle)
   dtlimit(:)  = dtlimit_default
   vectors(:)  = trim(shr_strdata_nullstr)
   fillalgo(:) = 'nn'
   fillmask(:) = 'nomask'
   fillread(:) = trim(shr_strdata_unset)
   fillwrite(:)= trim(shr_strdata_unset)
   mapalgo(:)  = 'bilinear'
   mapmask(:)  = 'dstmask'
   mapread(:)  = trim(shr_strdata_unset)
   mapwrite(:) = trim(shr_strdata_unset)
   tintalgo(:) = 'linear'


   !----------------------------------------------------------------------------
   ! read input namelist
   !----------------------------------------------------------------------------
   if (present(file)) then
      write(logunit,F00) 'reading input namelist file: ',trim(file)
      call shr_sys_flush(logunit)
      nUnit = shr_file_getUnit() ! get unused fortran i/o unit number
      open (nUnit,file=trim(file),status="old",action="read")
      call shr_nl_find_group_name(nUnit, 'shr_strdata_nml', status=rCode)
      if (rCode == 0) then
         read (nUnit, nml=shr_strdata_nml, iostat=rCode)
         if (rCode /= 0) then
            write(logunit,F01) 'ERROR: reading input namelist shr_strdata_input from file, &
                 &'//trim(file)//' iostat=',rCode
            call shr_sys_abort(subName//": namelist read error "//trim(file))
         end if
      end if
      close(nUnit)
      call shr_file_freeUnit(nUnit)
   endif

   !----------------------------------------------------------------------------
   ! copy temporary/local namelist vars into data structure
   !----------------------------------------------------------------------------
   SDAT%nstreams    = 0
   do n=1,nStrMax
      call shr_stream_default(SDAT%stream(n))
   enddo
   SDAT%dataMode    = dataMode
   SDAT%domainFile  = domainFile
   SDAT%streams(:)  = streams(:)
   SDAT%taxMode(:)  = taxMode(:)
   SDAT%dtlimit(:)  = dtlimit(:)
   SDAT%vectors(:)  = vectors(:)
   SDAT%fillalgo(:) = fillalgo(:)
   SDAT%fillmask(:) = fillmask(:)
   SDAT%fillread(:) = fillread(:)
   SDAT%fillwrit(:) = fillwrite(:)
   SDAT%mapalgo(:)  = mapalgo(:)
   SDAT%mapmask(:)  = mapmask(:)
   SDAT%mapread(:)  = mapread(:)
   SDAT%mapwrit(:)  = mapwrite(:)
   SDAT%tintalgo(:) = tintalgo(:)
   do n=1,nStrMax
      if (trim(streams(n)) /= trim(shr_strdata_nullstr)) SDAT%nstreams = max(SDAT%nstreams,n)
      if (trim(SDAT%taxMode(n)) == trim(shr_stream_taxis_extend)) SDAT%dtlimit(n) = 1.0e30
   end do
   SDAT%nvectors = 0
   do n=1,nVecMax
      if (trim(vectors(n)) /= trim(shr_strdata_nullstr)) SDAT%nvectors = n
   end do

   do n = 1,SDAT%nstreams
      call shr_stream_parseInput(SDAT%streams(n),fileName,yearAlign,yearFirst,yearLast)
      call shr_stream_init(SDAT%stream(n),fileName,yearFirst,yearLast,yearAlign, &
                   trim(SDAT%taxMode(n)))
   enddo

!   call shr_strdata_print(SDAT,trim(file)//' NML_ONLY')

   endif   ! master_task
!--master--task--

   if (present(mpicom)) then
      call shr_strdata_bcastnml(SDAT,mpicom)
   endif

   SDAT%ymdLB = -1
   SDAT%todLB = -1
   SDAT%ymdUB = -1
   SDAT%todUB = -1
   SDAT%dtmin = 1.0e30
   SDAT%dtmax = 0.0
   SDAT%nxg   = 0
   SDAT%nyg   = 0
   SDAT%nzg   = 0
   SDAT%eccen  = SHR_ORB_UNDEF_REAL
   SDAT%mvelpp = SHR_ORB_UNDEF_REAL
   SDAT%lambm0 = SHR_ORB_UNDEF_REAL
   SDAT%obliqr = SHR_ORB_UNDEF_REAL
   SDAT%modeldt = 0
   SDAT%calendar = shr_cal_noleap

end subroutine shr_strdata_readnml

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_strdata_pioinit -- initialize pio layer
!
! !DESCRIPTION:
!     Initialize PIO for a component model
!
! !REVISION HISTORY:
!     2010-10-26    Jim Edwards
!
! !INTERFACE: ------------------------------------------------------------------
subroutine shr_strdata_pioinit(SDAT,io_subsystem, io_type )
  type(shr_strdata_type),intent(inout):: SDAT  ! strdata data data-type
  type(iosystem_desc_t), pointer :: io_subsystem
  integer, intent(in) :: io_type

  SDAT%pio_subsystem => io_subsystem
  SDAT%io_type=io_type

end subroutine shr_strdata_pioinit

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_strdata_create -- set strdata and stream info from interface
!
! !DESCRIPTION:
!     Set strdata and stream info from fortran interface.
!     Note: When this is called, previous settings are reset to defaults
!           and then the values passed are used.
!
! !REVISION HISTORY:
!     2004-Dec-15 - J. Schramm - first version
!     2009-Apr-16 - T. Craig - add minimal parallel support
!
! !INTERFACE: ------------------------------------------------------------------
subroutine shr_strdata_create(SDAT, name, mpicom, compid, gsmap, ggrid, nxg, nyg, &
!--- streams stuff required ---
           yearFirst, yearLast, yearAlign, offset,          &
           domFilePath, domFileName,                        &
           domTvarName, domXvarName, domYvarName, domAreaName, domMaskName, &
           filePath, filename, fldListFile, fldListModel,   &
           pio_subsystem, pio_iotype,                       &
!--- strdata optional ---
           nzg, domZvarName,                                &
           taxMode, dtlimit, tintalgo,                      &
           fillalgo, fillmask, fillread, fillwrite,         &
           mapalgo, mapmask, mapread, mapwrite,             &
           calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_strdata_type),intent(inout):: SDAT  ! strdata data data-type
   character(*)          ,intent(in)   :: name  ! name of strdata
   integer(IN)           ,intent(in)   :: mpicom ! mpi comm
   integer(IN)           ,intent(in)   :: compid
   type(mct_gsmap)       ,intent(in)   :: gsmap
   type(mct_ggrid)       ,intent(in)   :: ggrid
   integer(IN)           ,intent(in)   :: nxg
   integer(IN)           ,intent(in)   :: nyg

   integer(IN)           ,intent(in)   :: yearFirst ! first year to use
   integer(IN)           ,intent(in)   :: yearLast  ! last  year to use
   integer(IN)           ,intent(in)   :: yearAlign ! align yearFirst with this model year
   integer(IN)           ,intent(in)   :: offset    ! offset in seconds of stream data
   character(*)          ,intent(in)   :: domFilePath  ! domain file path
   character(*)          ,intent(in)   :: domFileName  ! domain file name
   character(*)          ,intent(in)   :: domTvarName  ! domain time dim name
   character(*)          ,intent(in)   :: domXvarName  ! domain x dim name
   character(*)          ,intent(in)   :: domYvarName  ! domain y dim name
   character(*)          ,intent(in)   :: domAreaName  ! domain area name
   character(*)          ,intent(in)   :: domMaskName  ! domain mask name
   character(*)          ,intent(in)   :: filePath     ! path to filenames
   character(*)          ,intent(in)   :: filename(:)  ! filename for index filenumber
   character(*)          ,intent(in)   :: fldListFile  ! file field names, colon delim list
   character(*)          ,intent(in)   :: fldListModel ! model field names, colon delim list
   type(iosystem_desc_t), pointer      :: pio_subsystem ! PIO subsystem pointer
   integer(IN)          , intent(in)   :: pio_iotype    ! PIO file type

   integer(IN) ,optional ,intent(in)   :: nzg
   character(*),optional ,intent(in)   :: domZvarName  ! domain z dim name
   character(*),optional ,intent(in)   :: taxMode
   real(R8)    ,optional ,intent(in)   :: dtlimit
   character(*),optional ,intent(in)   :: fillalgo  ! fill algorithm
   character(*),optional ,intent(in)   :: fillmask  ! fill mask
   character(*),optional ,intent(in)   :: fillread  ! fill mapping file to read
   character(*),optional ,intent(in)   :: fillwrite ! fill mapping file to write
   character(*),optional ,intent(in)   :: mapalgo   ! scalar map algorithm
   character(*),optional ,intent(in)   :: mapmask   ! scalar map mask
   character(*),optional ,intent(in)   :: mapread   ! regrid mapping file to read
   character(*),optional ,intent(in)   :: mapwrite  ! regrid mapping file to write
   character(*),optional ,intent(in)   :: tintalgo  ! time interpolation algorithm
   character(*),optional, intent(in)   :: calendar

!EOP

!  --- local ---
!  --- formats ---
   character(*),parameter :: subName = "(shr_strdata_create) "
   character(*),parameter ::   F00 = "('(shr_strdata_create) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call shr_strdata_readnml(SDAT,mpicom=mpicom)

   SDAT%nstreams = 1

   call shr_strdata_pioinit(sdat, pio_subsystem, pio_iotype)

   if (present(taxMode)) then
      SDAT%taxMode(1) = taxMode
      if (trim(SDAT%taxMode(1)) == trim(shr_stream_taxis_extend)) SDAT%dtlimit(1) = 1.0e30
   endif
   if (present(dtlimit)) then
      SDAT%dtlimit(1) = dtlimit
   endif
   if (present(fillalgo)) then
      SDAT%fillalgo(1) = fillalgo
   endif
   if (present(fillmask)) then
      SDAT%fillmask(1) = fillmask
   endif
   if (present(fillread)) then
      SDAT%fillread(1) = fillread
   endif
   if (present(fillwrite)) then
      SDAT%fillwrit(1) = fillwrite
   endif
   if (present(mapalgo)) then
      SDAT%mapalgo(1) = mapalgo
   endif
   if (present(mapmask)) then
      SDAT%mapmask(1) = mapmask
   endif
   if (present(mapread)) then
      SDAT%mapread(1) = mapread
   endif
   if (present(mapwrite)) then
      SDAT%mapwrit(1) = mapwrite
   endif
   if (present(tintalgo)) then
      SDAT%tintalgo(1) = tintalgo
   endif
   if (present(mapmask)) then
      SDAT%mapmask(1) = mapmask
   endif
   if (present(calendar)) then
     SDAT%calendar = trim(shr_cal_calendarName(trim(calendar)))
   endif

   !---- Backwards compatibility requires Z be optional

   if (present(domZvarName)) then
      call shr_stream_set(SDAT%stream(1),yearFirst,yearLast,yearAlign,offset,taxMode, &
                          fldListFile,fldListModel,domFilePath,domFileName, &
                          domTvarName,domXvarName,domYvarName,domZvarName, &
                          domAreaName,domMaskName, &
                          filePath,filename,trim(name))
   else
      call shr_stream_set(SDAT%stream(1),yearFirst,yearLast,yearAlign,offset,taxMode, &
                          fldListFile,fldListModel,domFilePath,domFileName, &
                          domTvarName,domXvarName,domYvarName,'undefined', &
                          domAreaName,domMaskName, &
                          filePath,filename,trim(name))
   endif

   if (present(nzg)) then
      call shr_strdata_init(SDAT, mpicom, compid, &
           gsmap=gsmap, ggrid=ggrid, nxg=nxg, nyg=nyg, nzg=nzg)
   else
      call shr_strdata_init(SDAT, mpicom, compid, &
           gsmap=gsmap, ggrid=ggrid, nxg=nxg, nyg=nyg, nzg=1)
   endif

end subroutine shr_strdata_create

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_strdata_print -- read control strdata
!
! !DESCRIPTION:
!     Reads strdata common to all data models
!
! !REVISION HISTORY:
!     2004-Dec-15 - J. Schramm - first version
!     2009-Apr-16 - T. Craig - add minimal parallel support
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_strdata_print(SDAT,name)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_strdata_type)  ,intent(in) :: SDAT  ! strdata data data-type
   character(len=*),optional,intent(in) :: name  ! just a name for tracking

!EOP

   integer(IN)   :: n
   character(CL) :: lname

   !----- formats -----
   character(*),parameter :: subName = "(shr_strdata_print) "
   character(*),parameter ::   F00 = "('(shr_strdata_print) ',8a)"
   character(*),parameter ::   F01 = "('(shr_strdata_print) ',a,i6,a)"
   character(*),parameter ::   F02 = "('(shr_strdata_print) ',a,es13.6)"
   character(*),parameter ::   F03 = "('(shr_strdata_print) ',a,l6)"
   character(*),parameter ::   F04 = "('(shr_strdata_print) ',a,i2,a,a)"
   character(*),parameter ::   F05 = "('(shr_strdata_print) ',a,i2,a,i6)"
   character(*),parameter ::   F06 = "('(shr_strdata_print) ',a,i2,a,l2)"
   character(*),parameter ::   F07 = "('(shr_strdata_print) ',a,i2,a,es13.6)"
   character(*),parameter ::   F20 = "('(shr_strdata_print) ',a,i6,a)"
   character(*),parameter ::   F90 = "('(shr_strdata_print) ',58('-'))"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lname = 'unknown'
   if (present(name)) then
      lname = trim(name)
   endif
   !----------------------------------------------------------------------------
   ! document datatype settings
   !----------------------------------------------------------------------------
   write(logunit,F90)
   write(logunit,F00) "name        = ",trim(lname)
   write(logunit,F00) "dataMode    = ",trim(SDAT%dataMode)
   write(logunit,F00) "domainFile  = ",trim(SDAT%domainFile)
   write(logunit,F01) "nxg         = ",SDAT%nxg
   write(logunit,F01) "nyg         = ",SDAT%nyg
   write(logunit,F01) "nzg         = ",SDAT%nzg
   write(logunit,F00) "calendar    = ",trim(SDAT%calendar)
   write(logunit,F01) "io_type     = ",SDAT%io_type
   write(logunit,F02) "eccen       = ",SDAT%eccen
   write(logunit,F02) "mvelpp      = ",SDAT%mvelpp
   write(logunit,F02) "lambm0      = ",SDAT%lambm0
   write(logunit,F02) "obliqr      = ",SDAT%obliqr
   write(logunit,F01) "nstreams    = ",SDAT%nstreams
   write(logunit,F01) "pio_iotype  = ",sdat%io_type

   do n=1, SDAT%nstreams
     write(logunit,F04) "  streams (",n,") = ",trim(SDAT%streams(n))
     write(logunit,F04) "  taxMode (",n,") = ",trim(SDAT%taxMode(n))
     write(logunit,F07) "  dtlimit (",n,") = ",SDAT%dtlimit(n)
     write(logunit,F05) "  strnxg  (",n,") = ",SDAT%strnxg(n)
     write(logunit,F05) "  strnyg  (",n,") = ",SDAT%strnyg(n)
     write(logunit,F05) "  strnzg  (",n,") = ",SDAT%strnzg(n)
     write(logunit,F06) "  dofill  (",n,") = ",SDAT%dofill(n)
     write(logunit,F04) "  fillalgo(",n,") = ",trim(SDAT%fillalgo(n))
     write(logunit,F04) "  fillmask(",n,") = ",trim(SDAT%fillmask(n))
     write(logunit,F04) "  fillread(",n,") = ",trim(SDAT%fillread(n))
     write(logunit,F04) "  fillwrit(",n,") = ",trim(SDAT%fillwrit(n))
     write(logunit,F06) "  domaps  (",n,") = ",SDAT%domaps(n)
     write(logunit,F04) "  mapalgo (",n,") = ",trim(SDAT%mapalgo(n))
     write(logunit,F04) "  mapmask (",n,") = ",trim(SDAT%mapmask(n))
     write(logunit,F04) "  mapread (",n,") = ",trim(SDAT%mapread(n))
     write(logunit,F04) "  mapwrit (",n,") = ",trim(SDAT%mapwrit(n))
     write(logunit,F04) "  tintalgo(",n,") = ",trim(SDAT%tintalgo(n))
     write(logunit,F01) " "
   end do
   write(logunit,F01) "nvectors    = ",SDAT%nvectors
   do n=1, SDAT%nvectors
     write(logunit,F04) "  vectors (",n,") = ",trim(SDAT%vectors(n))
   end do
   write(logunit,F90)
   call shr_sys_flush(logunit)

end subroutine shr_strdata_print

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_strdata_bcastnml -- broadcast control strdata
!
! !DESCRIPTION:
!     Broadcast strdata
!
! !REVISION HISTORY:
!     2009-Apr-16 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_strdata_bcastnml(SDAT,mpicom,rc)

   use shr_mpi_mod, only : shr_mpi_bcast

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_strdata_type),intent(inout) :: SDAT  ! strdata data data-type
   integer(IN)           ,intent(in)  :: mpicom ! mpi communicator
   integer(IN),optional  ,intent(out) :: rc   ! return code

!EOP

   !----- local -----
   integer(IN) :: lrc

   !----- formats -----
   character(*),parameter :: subName = "(shr_strdata_bcastnml) "

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   lrc = 0

   call shr_mpi_bcast(SDAT%dataMode  ,mpicom,'dataMode')
   call shr_mpi_bcast(SDAT%domainFile,mpicom,'domainFile')
   call shr_mpi_bcast(SDAT%calendar  ,mpicom,'calendar')
   call shr_mpi_bcast(SDAT%nstreams  ,mpicom,'nstreams')
   call shr_mpi_bcast(SDAT%nvectors  ,mpicom,'nvectors')
   call shr_mpi_bcast(SDAT%streams   ,mpicom,'streams')
   call shr_mpi_bcast(SDAT%taxMode   ,mpicom,'taxMode')
   call shr_mpi_bcast(SDAT%dtlimit   ,mpicom,'dtlimit')
   call shr_mpi_bcast(SDAT%vectors   ,mpicom,'vectors')
   call shr_mpi_bcast(SDAT%fillalgo  ,mpicom,'fillalgo')
   call shr_mpi_bcast(SDAT%fillmask  ,mpicom,'fillmask')
   call shr_mpi_bcast(SDAT%fillread  ,mpicom,'fillread')
   call shr_mpi_bcast(SDAT%fillwrit  ,mpicom,'fillwrit')
   call shr_mpi_bcast(SDAT%mapalgo   ,mpicom,'mapalgo')
   call shr_mpi_bcast(SDAT%mapmask   ,mpicom,'mapmask')
   call shr_mpi_bcast(SDAT%mapread   ,mpicom,'mapread')
   call shr_mpi_bcast(SDAT%mapwrit   ,mpicom,'mapwrit')
   call shr_mpi_bcast(SDAT%tintalgo  ,mpicom,'tintalgo')

   if (present(rc)) then
      rc = lrc
   endif

end subroutine shr_strdata_bcastnml

!===============================================================================

subroutine shr_strdata_setlogunit(nu)

  integer(IN),intent(in) :: nu
  character(len=*),parameter :: subname = "(shr_strdata_setlogunit) "

  ! tcx DOES NOTHING, REMOVE

end subroutine shr_strdata_setlogunit

!===============================================================================
!===============================================================================

end module shr_strdata_mod

