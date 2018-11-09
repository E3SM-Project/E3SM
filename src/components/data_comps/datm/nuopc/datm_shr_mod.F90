module datm_shr_mod

  ! !USES:
  use shr_kind_mod   , only : IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, I8=>SHR_KIND_I8
  use shr_kind_mod   , only : CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_const_mod  , only : SHR_CONST_CDAY,SHR_CONST_TKFRZ,SHR_CONST_SPVAL
  use shr_file_mod   , only : shr_file_getlogunit, shr_file_getunit, shr_file_freeunit
  use shr_sys_mod    , only : shr_sys_flush, shr_sys_abort
  use shr_strdata_mod, only : shr_strdata_readnml
  use shr_dmodel_mod , only : shr_dmodel_mapset
  use shr_cal_mod    , only : shr_cal_date2julian
  use shr_ncread_mod , only : shr_ncread_varExists, shr_ncread_varDimSizes, shr_ncread_field4dG
  use shr_strdata_mod, only : shr_strdata_type
  use mct_mod

  ! !PUBLIC TYPES:
  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  interface datm_shr_getNextRadCday  
     module procedure datm_shr_getNextRadCDay_i8
     module procedure datm_shr_getNextRadCDay_i4
  end interface datm_shr_getNextRadCday

  public :: datm_shr_getNextRadCDay
  public :: datm_shr_CORE2getFactors
  public :: datm_shr_TN460getFactors
  public :: datm_shr_eSat
  public :: datm_shr_read_namelists

  !--------------------------------------------------------------------------
  ! Public data
  !--------------------------------------------------------------------------

  ! Note that model decomp will now come from reading in the mesh directly

  ! input namelist variables
  character(CL) , public :: restfilm              ! model restart file namelist
  character(CL) , public :: restfils              ! stream restart file namelist
  character(CL) , public :: bias_correct          ! true => send bias correction fields to coupler
  character(CL) , public :: anomaly_forcing(8)    ! true => send anomaly forcing fields to coupler
  logical       , public :: force_prognostic_true ! if true set prognostic true
  logical       , public :: wiso_datm = .false.   ! expect isotopic forcing from file?
  integer(IN)   , public :: iradsw                ! radiation interval
  character(CL) , public :: factorFn              ! file containing correction factors
  logical       , public :: presaero              ! true => send valid prescribe aero fields to coupler

  ! variables obtained from namelist read
  character(CL) , public :: rest_file             ! restart filename
  character(CL) , public :: rest_file_strm        ! restart filename for streams
  character(CL) , public :: datamode              ! mode
  character(len=*), public, parameter :: nullstr = 'undefined'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine datm_shr_read_namelists(filename, mpicom, my_task, master_task, &
       logunit, SDATM, atm_prognostic)

    ! !INPUT/OUTPUT PARAMETERS:
    character(len=*)       , intent(in)    :: filename       ! input namelist filename
    integer(IN)            , intent(in)    :: mpicom         ! mpi communicator
    integer(IN)            , intent(in)    :: my_task        ! my task in mpi communicator mpicom
    integer(IN)            , intent(in)    :: master_task    ! task number of master task
    integer(IN)            , intent(in)    :: logunit        ! logging unit number
    type(shr_strdata_type) , intent(inout) :: SDATM
    logical                , intent(out)   :: atm_prognostic ! flag

    !--- local variables ---
    integer(IN)   :: nunit       ! unit number
    integer(IN)   :: ierr        ! error code
    character(CL) :: decomp      ! decomp strategy - not used for NUOPC - but still needed in namelist for now

    !--- formats ---
    character(*), parameter :: F00   = "('(datm_comp_init) ',8a)"
    character(*), parameter :: F0L   = "('(datm_comp_init) ',a, l2)"
    character(*), parameter :: F01   = "('(datm_comp_init) ',a,5i8)"
    character(*), parameter :: subName = "(shr_datm_read_namelists) "
    !-------------------------------------------------------------------------------

    !----- define namelist -----
    namelist / datm_nml / decomp, &
         iradsw, factorFn, restfilm, restfils, presaero, bias_correct, &
         anomaly_forcing, force_prognostic_true, wiso_datm

    !----------------------------------------------------------------------------
    ! Read datm_in
    !----------------------------------------------------------------------------

    iradsw = 0
    factorFn = 'null'
    restfilm = trim(nullstr)
    restfils = trim(nullstr)
    presaero = .false.
    force_prognostic_true = .false.

    if (my_task == master_task) then
       nunit = shr_file_getUnit() ! get unused unit number
       open (nunit,file=trim(filename),status="old",action="read")
       read (nunit,nml=datm_nml,iostat=ierr)
       close(nunit)

       call shr_file_freeUnit(nunit)
       if (ierr > 0) then
          write(logunit,F01) 'ERROR: reading input namelist, '//trim(filename)//' iostat=',ierr
          call shr_sys_abort(subName//': namelist read error '//trim(filename))
       end if
       write(logunit,F01)' iradsw   = ',iradsw
       write(logunit,F00)' factorFn = ',trim(factorFn)
       write(logunit,F00)' restfilm = ',trim(restfilm)
       write(logunit,F00)' restfils = ',trim(restfils)
       write(logunit,F0L)' presaero = ',presaero
       write(logunit,F0L)' force_prognostic_true = ',force_prognostic_true
       write(logunit,F0L)' wiso_datm   = ',wiso_datm
       call shr_sys_flush(logunit)
    endif
    call shr_mpi_bcast(iradsw                ,mpicom, 'iradsw')
    call shr_mpi_bcast(factorFn              ,mpicom, 'factorFn')
    call shr_mpi_bcast(restfilm              ,mpicom, 'restfilm')
    call shr_mpi_bcast(restfils              ,mpicom, 'restfils')
    call shr_mpi_bcast(presaero              ,mpicom, 'presaero')
    call shr_mpi_bcast(force_prognostic_true ,mpicom, 'force_prognostic_true')
    call shr_mpi_bcast(wiso_datm             ,mpicom, 'wiso_datm')

    rest_file = trim(restfilm)
    rest_file_strm = trim(restfils)

    !----------------------------------------------------------------------------
    ! Read dshr namelist
    !----------------------------------------------------------------------------

    call shr_strdata_readnml(SDATM, trim(filename), mpicom=mpicom)

    ! Validate mode

    datamode = trim(SDATM%dataMode)
    if (trim(datamode) == 'NULL'      .or. &
        trim(datamode) == 'CORE2_NYF' .or. &
        trim(datamode) == 'CORE2_IAF' .or. &
        trim(datamode) == 'CORE_IAF_JRA' .or. &
        trim(datamode) == 'CLMNCEP'   .or. &
        trim(datamode) == 'COPYALL'   ) then
       if (my_task == master_task) then
          write(logunit,F00) ' datm datamode = ',trim(datamode)
          call shr_sys_flush(logunit)
       end if
    else
       write(logunit,F00) ' ERROR illegal datm datamode = ',trim(datamode)
       call shr_sys_abort()
    endif

    !----------------------------------------------------------------------------
    ! Determine present and prognostic flag
    !----------------------------------------------------------------------------

    atm_prognostic = .false.
    if (force_prognostic_true) then
       atm_prognostic = .true.
    endif

  end subroutine datm_shr_read_namelists

  !===============================================================================
  real(R8) function datm_shr_getNextRadCDay_i8( ymd, tod, stepno, dtime, iradsw, calendar )

    !  Return the calendar day of the next radiation time-step.
    !  General Usage: nextswday = datm_shr_getNextRadCDay(curr_date)

    ! !INPUT/OUTPUT PARAMETERS:
    integer(IN), intent(in)    :: ymd
    integer(IN), intent(in)    :: tod
    integer(I8), intent(in)    :: stepno
    integer(IN), intent(in)    :: dtime
    integer(IN), intent(in)    :: iradsw
    character(*),intent(in)    :: calendar

    !----- local -----
    real(R8) :: nextsw_cday
    real(R8) :: julday
    integer  :: liradsw
    integer  :: yy,mm,dd
    character(*),parameter :: subName =  '(datm_shr_getNextRadCDay) '
    !-------------------------------------------------------------------------------

    liradsw = iradsw
    if (liradsw < 0) liradsw  = nint((-liradsw *3600._r8)/dtime)

    call shr_cal_date2julian(ymd,tod,julday,calendar)

    if (liradsw > 1) then
       if (mod(stepno+1,liradsw) == 0 .and. stepno > 0) then
          nextsw_cday = julday + 2*dtime/SHR_CONST_CDAY
       else
          nextsw_cday = -1._r8
       end if
    else
       nextsw_cday = julday + dtime/SHR_CONST_CDAY
    end if
    datm_shr_getNextRadCDay_i8 = nextsw_cday

  end function datm_shr_getNextRadCDay_i8

  real(R8) function datm_shr_getNextRadCDay_i4( ymd, tod, stepno, dtime, iradsw, calendar )

    !  Return the calendar day of the next radiation time-step.
    !  General Usage: nextswday = datm_shr_getNextRadCDay(curr_date)

    ! !INPUT/OUTPUT PARAMETERS:
    integer(IN), intent(in)    :: ymd
    integer(IN), intent(in)    :: tod
    integer    , intent(in)    :: stepno
    integer(IN), intent(in)    :: dtime
    integer(IN), intent(in)    :: iradsw
    character(*),intent(in)    :: calendar

    !----- local -----
    real(R8) :: nextsw_cday
    real(R8) :: julday
    integer  :: liradsw
    integer  :: yy,mm,dd
    character(*),parameter :: subName =  '(datm_shr_getNextRadCDay) '
    !-------------------------------------------------------------------------------

    liradsw = iradsw
    if (liradsw < 0) liradsw  = nint((-liradsw *3600._r8)/dtime)

    call shr_cal_date2julian(ymd,tod,julday,calendar)

    if (liradsw > 1) then
       if (mod(stepno+1,liradsw) == 0 .and. stepno > 0) then
          nextsw_cday = julday + 2*dtime/SHR_CONST_CDAY
       else
          nextsw_cday = -1._r8
       end if
    else
       nextsw_cday = julday + dtime/SHR_CONST_CDAY
    end if
    datm_shr_getNextRadCDay_i4 = nextsw_cday

  end function datm_shr_getNextRadCDay_i4

  !===============================================================================
  subroutine datm_shr_CORE2getFactors(fileName,windF,winddF,qsatF,mpicom,compid, &
       gsmap,ggrid,nxg,nyg)

    !--- arguments ---
    character(*)    ,intent(in)    :: fileName   ! file name string
    real(R8)        ,intent(inout) :: windF(:)   ! wind adjustment factor
    real(R8)        ,intent(inout) :: winddF(:)  ! wind adjustment factor
    real(R8)        ,intent(inout) :: qsatF(:)   ! rel humidty adjustment factors
    integer(IN)     ,intent(in)    :: mpicom     ! mpi comm
    integer(IN)     ,intent(in)    :: compid     ! mct compid
    type(mct_gsmap) ,intent(in)    :: gsmap      ! decomp of wind,windd,qsat
    type(mct_ggrid) ,intent(in)    :: ggrid      ! ggrid of grid info
    integer(IN)     ,intent(in)    :: nxg        ! size of input grid
    integer(IN)     ,intent(in)    :: nyg        ! size of input grid

    !--- local ---
    integer(IN) :: my_task,logunit,ier
    character(*),parameter :: subName =  '(datm_shr_CORE2getFactors) '
    character(*),parameter :: F00    = "('(datm_shr_CORE2getFactors) ',4a) "
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ier)
    call shr_file_getLogUnit(logunit)

    if (my_task == 0) then

       !--- verify necessary data is in input file ---
       if ( .not. shr_ncread_varExists(fileName,'lat'       )  &
            .or. .not. shr_ncread_varExists(fileName,'lon'       )  &
            .or. .not. shr_ncread_varExists(fileName,'mask'      )  &
            .or. .not. shr_ncread_varExists(fileName,'windFactor')  &
            .or. .not. shr_ncread_varExists(fileName,'winddFactor') &
            .or. .not. shr_ncread_varExists(fileName,'qsatFactor')  ) then
          write(logunit,F00) "ERROR: invalid correction factor data file"
          call shr_sys_abort(subName//"invalid correction factor data file")
       end if
    endif

    call datm_shr_getFactors(fileName,windF,winddF,qsatF,mpicom,compid, &
         gsmap,ggrid,nxg,nyg)

  end subroutine datm_shr_CORE2getFactors

  !===============================================================================

  subroutine datm_shr_TN460getFactors(fileName,windF,qsatF,mpicom,compid, &
       gsmap,ggrid,nxg,nyg)

    !--- arguments ---
    character(*)    ,intent(in)    :: fileName   ! file name string
    real(R8)        ,intent(inout) :: windF(:)   ! wind adjustment factor
    real(R8)        ,intent(inout) :: qsatF(:)   ! rel humidty adjustment factors
    integer(IN)     ,intent(in)    :: mpicom     ! mpi comm
    integer(IN)     ,intent(in)    :: compid     ! mct compid
    type(mct_gsmap) ,intent(in)    :: gsmap      ! decomp of wind,windd,qsat
    type(mct_ggrid) ,intent(in)    :: ggrid      ! ggrid of grid info
    integer(IN)     ,intent(in)    :: nxg        ! size of input grid
    integer(IN)     ,intent(in)    :: nyg        ! size of input grid

    !--- local ---
    integer(IN) :: my_task,logunit,ier
    real(R8),pointer :: winddF(:)  ! wind adjustment factor
    character(*),parameter :: subName =  '(datm_shr_TN460getFactors) '
    character(*),parameter :: F00    = "('(datm_shr_TN460getFactors) ',4a) "
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ier)
    call shr_file_getLogUnit(logunit)

    if (my_task == 0) then

       !--- verify necessary data is in input file ---
       if ( .not. shr_ncread_varExists(fileName,'lat'       )  &
            .or. .not. shr_ncread_varExists(fileName,'lon'       )  &
            .or. .not. shr_ncread_varExists(fileName,'mask'      )  &
            .or. .not. shr_ncread_varExists(fileName,'windFactor')  &
            .or. .not. shr_ncread_varExists(fileName,'qsatFactor')  ) then
          write(logunit,F00) "ERROR: invalid correction factor data file"
          call shr_sys_abort(subName//"invalid correction factor data file")
       end if
    endif

    call datm_shr_getFactors(fileName,windF,winddF,qsatF,mpicom,compid, &
         gsmap,ggrid,nxg,nyg)

  end subroutine datm_shr_TN460getFactors

  !===============================================================================

  subroutine datm_shr_getFactors(fileName,windF,winddF,qsatF,mpicom,compid, &
       gsmapo,ggrido,nxgo,nygo)

    use shr_map_mod

    !--- arguments ---
    character(*)    ,intent(in)    :: fileName   ! file name string
    real(R8)        ,intent(inout) :: windF(:)   ! wind adjustment factor
    real(R8)        ,intent(inout) :: winddF(:)  ! wind adjustment factor
    real(R8)        ,intent(inout) :: qsatF(:)   ! rel humidty adjustment factors
    integer(IN)     ,intent(in)    :: mpicom     ! mpi comm
    integer(IN)     ,intent(in)    :: compid     ! mct compid
    type(mct_gsmap) ,intent(in)    :: gsmapo     ! decomp of wind,windd,qsat
    type(mct_ggrid) ,intent(in)    :: ggrido     ! ggrid of grid info
    integer(IN)     ,intent(in)    :: nxgo       ! size of input grid
    integer(IN)     ,intent(in)    :: nygo       ! size of input grid

    !--- data that describes the local model domain ---
    integer(IN)      :: ni0,nj0       ! dimensions of global bundle0
    integer(IN)      :: ni1,nj1,nf1   ! dimensions of global bundle1
    integer(IN)      :: i,j,n         ! generic indicies
    integer(IN)      :: my_task       ! local pe number
    integer(IN)      :: ier           ! error code
    integer(IN)      :: logunit       ! logging unit
    type(mct_ggrid)  :: ggridi        ! input file grid
    type(mct_ggrid)  :: ggridoG       ! output grid gathered
    type(mct_gsmap)  :: gsmapi        ! input file gsmap
    type(mct_sMatp)  :: smatp         ! sparse matrix weights
    type(mct_avect)  :: avi           ! input attr vect
    type(mct_avect)  :: avo           ! output attr vect
    integer(IN)      :: lsizei        ! local size of input
    integer(IN)      :: lsizeo        ! local size of output
    integer(IN),pointer :: start(:)   ! start list
    integer(IN),pointer :: length(:)  ! length list
    integer(IN)      :: gsizei        ! input global size
    integer(IN)      :: numel         ! number of elements in start list
    real(R8)         :: dadd          ! lon correction
    logical          :: domap         ! map or not
    integer(IN)      :: klon,klat     ! lon lat fld index

    !--- temp arrays for data input ---
    real(R8)   ,allocatable :: tempR4D(:,:,:,:)   ! 4D data array
    real(R8)   ,pointer     :: tempR1D(:)         ! 1D data array
    integer(IN),allocatable :: tempI4D(:,:,:,:)   ! 4D data array
    character(*),parameter  :: subName =  '(datm_shr_getFactors) '
    character(*),parameter  :: F00    = "('(datm_shr_getFactors) ',4a) "
    character(*),parameter  :: F01    = "('(datm_shr_getFactors) ',a,2i5)"
    character(*),parameter  :: F02    = "('(datm_shr_getFactors) ',a,6e12.3)"

    !-------------------------------------------------------------------------------
    !   Note: gsmapi is all gridcells on root pe
    !-------------------------------------------------------------------------------

    call MPI_COMM_RANK(mpicom,my_task,ier)
    call shr_file_getLogUnit(logunit)

    ni0 = 0
    nj0 = 0
    allocate(start(1),length(1))
    start = 0
    length = 0
    numel = 0

    !----------------------------------------------------------------------------
    ! read in and map global correction factors
    !----------------------------------------------------------------------------
    if (my_task == 0) then

       !--- verify necessary data is in input file ---
       if ( .not. shr_ncread_varExists(fileName,'lat'       )  &
            .or. .not. shr_ncread_varExists(fileName,'lon'       )  &
            .or. .not. shr_ncread_varExists(fileName,'mask'      )  &
            .or. .not. shr_ncread_varExists(fileName,'windFactor')  ) then
          write(logunit,F00) "ERROR: invalid correction factor data file"
          call shr_sys_abort(subName//"invalid correction factor data file")
       end if
       call shr_ncread_varDimSizes(fileName,"windFactor",ni0,nj0)
       start = 1
       length = ni0*nj0
       numel = 1
    endif
    call shr_mpi_bcast(ni0,mpicom,subname//' ni0')
    call shr_mpi_bcast(nj0,mpicom,subname//' nj0')
    gsizei = ni0*nj0

    !--- allocate datatypes for input data ---
    call mct_gsmap_init(gsmapi,start,length,0,mpicom,compid,gsize=gsizei,numel=numel)
    deallocate(start,length)
    lsizei = mct_gsmap_lsize(gsmapi,mpicom)
    lsizeo = mct_gsmap_lsize(gsmapo,mpicom)
    call mct_gGrid_init(GGrid=gGridi, &
         CoordChars='lat:lon:hgt', OtherChars='area:aream:mask:frac', lsize=lsizei )
    call mct_aVect_init(avi,rList="wind:windd:qsat",lsize=lsizei)
    avi%rAttr = SHR_CONST_SPVAL

    !--- gather output grid for map logic ---
    call mct_ggrid_gather(ggrido, ggridoG, gsmapo, 0, mpicom)

    if (my_task == 0) then
       allocate(tempR1D(ni0*nj0))

       !--- read domain data: lon ---
       allocate(tempR4D(ni0,1,1,1))
       call shr_ncread_field4dG(fileName,'lon' ,rfld=tempR4D)
       !--- needs to be monotonically increasing, add 360 at wraparound+ ---
       dadd = 0.0_R8
       do i = 2,ni0
          if (tempR4D(i-1,1,1,1) > tempR4D(i,1,1,1)) dadd = 360.0_R8
          tempR4D(i,1,1,1) = tempR4D(i,1,1,1) + dadd
       enddo
       n = 0
       do j=1,nj0
          do i=1,ni0
             n = n + 1
             tempR1D(n) = tempR4D(i,1,1,1)
          end do
       end do
       deallocate(tempR4D)
       call mct_gGrid_importRattr(gGridi,'lon',tempR1D,lsizei)

       !--- read domain data: lat ---
       allocate(tempR4D(nj0,1,1,1))
       call shr_ncread_field4dG(fileName,'lat' ,rfld=tempR4D)
       n = 0
       do j=1,nj0
          do i=1,ni0
             n = n + 1
             tempR1D(n) = tempR4D(j,1,1,1)
          end do
       end do
       deallocate(tempR4D)
       call mct_gGrid_importRattr(gGridi,'lat',tempR1D,lsizei)

       !--- read domain mask---
       allocate(tempI4D(ni0,nj0,1,1))
       call shr_ncread_field4dG(fileName,'mask',ifld=tempI4D)
       n = 0
       do j=1,nj0
          do i=1,ni0
             n = n + 1
             tempR1D(n) = real(tempI4D(i,j,1,1),R8)
          end do
       end do
       deallocate(tempI4D)
       call mct_gGrid_importRattr(gGridi,'mask',tempR1D,lsizei)

       !--- read bundle data: wind factor ---
       allocate(tempR4D(ni0,nj0,1,1))
       if (shr_ncread_varExists(fileName,'windFactor')  ) then
          call shr_ncread_field4dG(fileName,'windFactor',rfld=tempR4D)
          n = 0
          do j=1,nj0
             do i=1,ni0
                n = n + 1
                tempR1D(n) = tempR4D(i,j,1,1)
             end do
          end do
          call mct_aVect_importRattr(avi,'wind',tempR1D,lsizei)
       endif

       !--- read bundle data: windd factor ---
       if (shr_ncread_varExists(fileName,'winddFactor')  ) then
          call shr_ncread_field4dG(fileName,'winddFactor',rfld=tempR4D)
          n = 0
          do j=1,nj0
             do i=1,ni0
                n = n + 1
                tempR1D(n) = tempR4D(i,j,1,1)
             end do
          end do
          call mct_aVect_importRattr(avi,'windd',tempR1D,lsizei)
       endif

       !--- read bundle data: qsat factor ---
       if (shr_ncread_varExists(fileName,'qsatFactor')  ) then
          call shr_ncread_field4dG(fileName,'qsatFactor',rfld=tempR4D)
          n = 0
          do j=1,nj0
             do i=1,ni0
                n = n + 1
                tempR1D(n) = tempR4D(i,j,1,1)
             end do
          end do
          call mct_aVect_importRattr(avi,'qsat',tempR1D,lsizei)
       endif

       deallocate(tempR4D)
       deallocate(tempR1D)

       domap = .false.
       if (ni0 /= nxgo .or. nj0 /= nygo) then
          domap = .true.
       else
          klon = mct_aVect_indexRA(ggridi%data,'lon')
          klat = mct_aVect_indexRA(ggrido%data,'lat')
          do n = 1,lsizei
             if (abs(ggridi%data%rAttr(klon,n)-ggridoG%data%rAttr(klon,n)) > 0.01_R8) domap=.true.
             if (abs(ggridi%data%rAttr(klat,n)-ggridoG%data%rAttr(klat,n)) > 0.01_R8) domap=.true.
          enddo
       endif

       call mct_gGrid_clean(ggridoG)

    endif

    call shr_mpi_bcast(domap,mpicom,subname//' domap')

    if (domap) then
       call shr_dmodel_mapSet(smatp,ggridi,gsmapi,ni0 ,nj0 , &
            ggrido,gsmapo,nxgo,nygo, &
            'datmfactor',shr_map_fs_remap,shr_map_fs_bilinear,shr_map_fs_srcmask, &
            shr_map_fs_scalar,compid,mpicom,'Xonly')
       call mct_aVect_init(avo,avi,lsizeo)
       call mct_sMat_avMult(avi,smatp,avo)
       call mct_sMatP_clean(smatp)
    else
       call mct_aVect_scatter(avi,avo,gsmapo,0,mpicom)
    endif

    !--- fill the interface arrays, only if they are the right size ---
    allocate(tempR1D(lsizeo))
    if (size(windF ) >= lsizeo) then
       call mct_aVect_exportRattr(avo,'wind' ,tempR1D,lsizeo)
       windF = tempR1D
    endif
    if (size(winddF) >= lsizeo) then
       call mct_aVect_exportRattr(avo,'windd',tempR1D,lsizeo)
       winddF = tempR1D
    endif
    if (size(qsatF ) >= lsizeo) then
       call mct_aVect_exportRattr(avo,'qsat' ,tempR1D,lsizeo)
       qsatF = tempR1D
    endif
    deallocate(tempR1D)

    call mct_aVect_clean(avi)
    call mct_aVect_clean(avo)
    call mct_gGrid_clean(ggridi)
    call mct_gsmap_clean(gsmapi)

  end subroutine datm_shr_getFactors

  !===============================================================================

  real(R8) function datm_shr_eSat(tK,tKbot)

    !--- arguments ---
    real(R8),intent(in) :: tK    ! temp used in polynomial calculation
    real(R8),intent(in) :: tKbot ! bottom atm temp

    !--- local ---
    real(R8)           :: t     ! tK converted to Celcius
    real(R8),parameter :: tkFrz = SHR_CONST_TKFRZ  ! freezing T of fresh water ~ K

    !--- coefficients for esat over water ---
    real(R8),parameter :: a0=6.107799961_R8
    real(R8),parameter :: a1=4.436518521e-01_R8
    real(R8),parameter :: a2=1.428945805e-02_R8
    real(R8),parameter :: a3=2.650648471e-04_R8
    real(R8),parameter :: a4=3.031240396e-06_R8
    real(R8),parameter :: a5=2.034080948e-08_R8
    real(R8),parameter :: a6=6.136820929e-11_R8

    !--- coefficients for esat over ice ---
    real(R8),parameter :: b0=6.109177956_R8
    real(R8),parameter :: b1=5.034698970e-01_R8
    real(R8),parameter :: b2=1.886013408e-02_R8
    real(R8),parameter :: b3=4.176223716e-04_R8
    real(R8),parameter :: b4=5.824720280e-06_R8
    real(R8),parameter :: b5=4.838803174e-08_R8
    real(R8),parameter :: b6=1.838826904e-10_R8

    !----------------------------------------------------------------------------
    ! use polynomials to calculate saturation vapor pressure and derivative with
    ! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
    ! required to convert relative humidity to specific humidity
    !----------------------------------------------------------------------------

    t = min( 50.0_R8, max(-50.0_R8,(tK-tKfrz)) )
    if ( tKbot < tKfrz) then
       datm_shr_eSat = 100.0_R8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
    else
       datm_shr_eSat = 100.0_R8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    end if

  end function datm_shr_eSat

end module datm_shr_mod
