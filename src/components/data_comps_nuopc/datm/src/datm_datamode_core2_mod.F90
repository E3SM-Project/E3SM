module datm_datamode_core2_mod

  use ESMF
  use NUOPC            , only : NUOPC_Advertise
  use shr_kind_mod     , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod      , only : shr_sys_abort
  use shr_cal_mod      , only : shr_cal_date2julian
  use shr_const_mod    , only : shr_const_tkfrz, shr_const_pi
  use dshr_strdata_mod , only : shr_strdata_get_stream_pointer, shr_strdata_type
  use dshr_methods_mod , only : dshr_state_getfldptr, dshr_fldbun_getfldptr, dshr_fldbun_regrid, chkerr
  use dshr_mod         , only : dshr_restart_read, dshr_restart_write
  use dshr_strdata_mod , only : shr_strdata_type
  use dshr_fldlist_mod , only : fldlist_type, dshr_fldlist_add
  use pio

  implicit none
  private ! except

  public  :: datm_datamode_core2_advertise
  public  :: datm_datamode_core2_init_pointers
  public  :: datm_datamode_core2_advance
  public  :: datm_datamode_core2_restart_write
  public  :: datm_datamode_core2_restart_read

  private :: datm_get_adjustment_factors

  ! export state pointers
  real(r8), pointer :: Sa_u(:)       => null()
  real(r8), pointer :: Sa_v(:)       => null()
  real(r8), pointer :: Sa_z(:)       => null()
  real(r8), pointer :: Sa_tbot(:)    => null()
  real(r8), pointer :: Sa_ptem(:)    => null()
  real(r8), pointer :: Sa_shum(:)    => null()
  real(r8), pointer :: Sa_pbot(:)    => null()
  real(r8), pointer :: Sa_pslv(:)    => null()
  real(r8), pointer :: Faxa_lwdn(:)  => null()
  real(r8), pointer :: Faxa_rainc(:) => null()
  real(r8), pointer :: Faxa_rainl(:) => null()
  real(r8), pointer :: Faxa_snowc(:) => null()
  real(r8), pointer :: Faxa_snowl(:) => null()
  real(r8), pointer :: Faxa_swndr(:) => null()
  real(r8), pointer :: Faxa_swndf(:) => null()
  real(r8), pointer :: Faxa_swvdr(:) => null()
  real(r8), pointer :: Faxa_swvdf(:) => null()
  real(r8), pointer :: Faxa_swnet(:) => null()
  real(r8), pointer :: Faxa_swdn(:)  => null()

  ! stream data
  real(r8), pointer :: strm_prec(:)      => null()
  real(r8), pointer :: strm_swdn(:)      => null()
  real(r8), pointer :: strm_tarcf(:)     => null()

  ! othe module arrays
  real(R8), pointer :: windFactor(:)
  real(R8), pointer :: winddFactor(:)
  real(R8), pointer :: qsatFactor(:)
  real(R8), pointer :: yc(:)                 ! array of model latitudes

  ! constants
  real(R8) , parameter :: tKFrz    = SHR_CONST_TKFRZ
  real(R8) , parameter :: degtorad = SHR_CONST_PI/180.0_R8
  real(R8) , parameter :: avg_c0   =  61.846_R8
  real(R8) , parameter :: avg_c1   =   1.107_R8
  real(R8) , parameter :: amp_c0   = -21.841_R8
  real(R8) , parameter :: amp_c1   =  -0.447_R8
  real(R8) , parameter :: phs_c0   =   0.298_R8
  real(R8) , parameter :: dLWarc   =  -5.000_R8

  real(R8) :: dTarc(12)
  data   dTarc      / 0.49_R8, 0.06_R8,-0.73_R8,  -0.89_R8,-0.77_R8,-1.02_R8, &
                     -1.99_R8,-0.91_R8, 1.72_R8,   2.30_R8, 1.81_R8, 1.06_R8/

  character(*), parameter :: nullstr = 'undefined'
  character(*), parameter :: rpfile  = 'rpointer.atm'
  character(*), parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine datm_datamode_core2_advertise(exportState, fldsexport, flds_scalar_name, &
       flds_co2, flds_wiso, flds_presaero, rc)

    ! input/output variables
    type(esmf_State)   , intent(inout) :: exportState
    type(fldlist_type) , pointer       :: fldsexport
    character(len=*)   , intent(in)    :: flds_scalar_name
    logical            , intent(in)    :: flds_co2
    logical            , intent(in)    :: flds_wiso
    logical            , intent(in)    :: flds_presaero
    integer            , intent(out)   :: rc

    ! local variables
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call dshr_fldList_add(fldsExport, trim(flds_scalar_name))
    call dshr_fldList_add(fldsExport, 'Sa_z'       )
    call dshr_fldList_add(fldsExport, 'Sa_u'       )
    call dshr_fldList_add(fldsExport, 'Sa_v'       )
    call dshr_fldList_add(fldsExport, 'Sa_ptem'    )
    call dshr_fldList_add(fldsExport, 'Sa_dens'    )
    call dshr_fldList_add(fldsExport, 'Sa_pslv'    )
    call dshr_fldList_add(fldsExport, 'Sa_tbot'    )
    call dshr_fldList_add(fldsExport, 'Sa_pbot'    )
    call dshr_fldList_add(fldsExport, 'Sa_shum'    )
    call dshr_fldList_add(fldsExport, 'Faxa_rainc' )
    call dshr_fldList_add(fldsExport, 'Faxa_rainl' )
    call dshr_fldList_add(fldsExport, 'Faxa_snowc' )
    call dshr_fldList_add(fldsExport, 'Faxa_snowl' )
    call dshr_fldList_add(fldsExport, 'Faxa_swndr' )
    call dshr_fldList_add(fldsExport, 'Faxa_swvdr' )
    call dshr_fldList_add(fldsExport, 'Faxa_swndf' )
    call dshr_fldList_add(fldsExport, 'Faxa_swvdf' )
    call dshr_fldList_add(fldsExport, 'Faxa_swnet' )
    call dshr_fldList_add(fldsExport, 'Faxa_lwdn'  )
    call dshr_fldList_add(fldsExport, 'Faxa_swdn'  )

    if (flds_co2) then
       call dshr_fldList_add(fldsExport, 'Sa_co2prog')
       call dshr_fldList_add(fldsExport, 'Sa_co2diag')
    end if
    if (flds_presaero) then
       call dshr_fldList_add(fldsExport, 'Faxa_bcph'   , ungridded_lbound=1, ungridded_ubound=3)
       call dshr_fldList_add(fldsExport, 'Faxa_ocph'   , ungridded_lbound=1, ungridded_ubound=3)
       call dshr_fldList_add(fldsExport, 'Faxa_dstwet' , ungridded_lbound=1, ungridded_ubound=4)
       call dshr_fldList_add(fldsExport, 'Faxa_dstdry' , ungridded_lbound=1, ungridded_ubound=4)
    end if
    if (flds_wiso) then
       call dshr_fldList_add(fldsExport, 'Faxa_rainc_wiso', ungridded_lbound=1, ungridded_ubound=3)
       call dshr_fldList_add(fldsExport, 'Faxa_rainl_wiso', ungridded_lbound=1, ungridded_ubound=3)
       call dshr_fldList_add(fldsExport, 'Faxa_snowc_wiso', ungridded_lbound=1, ungridded_ubound=3)
       call dshr_fldList_add(fldsExport, 'Faxa_snowl_wiso', ungridded_lbound=1, ungridded_ubound=3)
       call dshr_fldList_add(fldsExport, 'Faxa_shum_wiso' , ungridded_lbound=1, ungridded_ubound=3)
    end if

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(datm_comp_advertise): Fr_atm'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

  end subroutine datm_datamode_core2_advertise

  !===============================================================================
  subroutine datm_datamode_core2_init_pointers(exportState, sdat, datamode, factorfn_mesh, factorfn_data, rc)

    ! input/output variables
    type(ESMF_State)       , intent(inout) :: exportState
    type(shr_strdata_type) , intent(in)    :: sdat
    character(len=*)       , intent(in)    :: datamode
    character(len=*)       , intent(in)    :: factorfn_mesh
    character(len=*)       , intent(in)    :: factorfn_data
    integer                , intent(out)   :: rc

    ! local variables
    integer           :: n
    integer           :: lsize
    integer           :: spatialDim         ! number of dimension in mesh
    integer           :: numOwnedElements   ! size of mesh
    real(r8), pointer :: ownedElemCoords(:) ! mesh lat and lons
    character(len=*), parameter :: subname='(datm_init_pointers): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lsize = sdat%model_lsize

    ! allocate module arrays
    allocate(windFactor(lsize))
    allocate(winddFactor(lsize))
    allocate(qsatFactor(lsize))

    call ESMF_MeshGet(sdat%model_mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(ownedElemCoords(spatialDim*numOwnedElements))
    allocate(yc(numOwnedElements))
    call ESMF_MeshGet(sdat%model_mesh, ownedElemCoords=ownedElemCoords)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    do n = 1,numOwnedElements
       yc(n) = ownedElemCoords(2*n)
    end do

    ! get stream pointers
    call shr_strdata_get_stream_pointer( sdat, 'Faxa_prec'  , strm_prec  , rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_strdata_get_stream_pointer( sdat, 'Faxa_swdn'  , strm_swdn  , rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_strdata_get_stream_pointer( sdat, 'tarcf'      , strm_tarcf , rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! get export state pointers
    call dshr_state_getfldptr(exportState, 'Sa_z'       , fldptr1=Sa_z       , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Sa_u'       , fldptr1=Sa_u       , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Sa_v'       , fldptr1=Sa_v       , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Sa_tbot'    , fldptr1=Sa_tbot    , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Sa_pbot'    , fldptr1=Sa_pbot    , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Sa_pslv'    , fldptr1=Sa_pslv    , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Sa_ptem'    , fldptr1=Sa_ptem    , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Sa_shum'    , fldptr1=Sa_shum    , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Faxa_rainc' , fldptr1=Faxa_rainc , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Faxa_rainl' , fldptr1=Faxa_rainl , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Faxa_snowc' , fldptr1=Faxa_snowc , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Faxa_snowl' , fldptr1=Faxa_snowl , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Faxa_swvdr' , fldptr1=Faxa_swvdr , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Faxa_swvdf' , fldptr1=Faxa_swvdf , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Faxa_swndr' , fldptr1=Faxa_swndr , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Faxa_swndf' , fldptr1=Faxa_swndf , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Faxa_swnet' , fldptr1=Faxa_swnet , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Faxa_lwdn'  , fldptr1=Faxa_lwdn  , rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (.not. associated(strm_prec) .or. .not. associated(strm_swdn)) then
       call shr_sys_abort(trim(subname)//'ERROR: prec and swdn must be in streams for CORE2')
    endif

    if (trim(datamode) == 'CORE2_IAF' ) then
       if (.not. associated(strm_tarcf)) then
          call shr_sys_abort(trim(subname)//'tarcf must be in an input stream for CORE2_IAF')
       endif
    endif

    ! create adjustment factor arrays
    call datm_get_adjustment_factors(sdat, factorFn_mesh, factorFn_data, windFactor, winddFactor, qsatFactor, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine datm_datamode_core2_init_pointers

  !===============================================================================
  subroutine datm_datamode_core2_advance(exportstate, datamode, target_ymd, target_tod, target_mon, &
       model_calendar, factorfn_mesh, factorfn_data, rc)

    ! input/output variables
    type(ESMF_State)       , intent(inout) :: exportState
    character(len=*)       , intent(in)    :: datamode
    integer                , intent(in)    :: target_ymd
    integer                , intent(in)    :: target_tod
    integer                , intent(in)    :: target_mon
    character(len=*)       , intent(in)    :: model_calendar
    character(len=*)       , intent(in)    :: factorfn_mesh
    character(len=*)       , intent(in)    :: factorfn_data
    integer                , intent(out)   :: rc

    ! local variables
    integer  :: n
    integer  :: lsize
    real(R8) :: avg_alb            ! average albedo
    real(R8) :: rday               ! elapsed day
    real(R8) :: cosFactor          ! cosine factor
    real(R8) :: factor             ! generic/temporary correction factor
    real(R8) :: tMin               ! minimum temperature
    real(R8) :: uprime,vprime
    character(len=*), parameter :: subname='(datm_datamode_core2): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lsize = size(Sa_z)

    call shr_cal_date2julian(target_ymd, target_tod, rday, model_calendar)
    rday = mod((rday - 1.0_R8),365.0_R8)
    cosfactor = cos((2.0_R8*SHR_CONST_PI*rday)/365 - phs_c0)

    do n = 1,lsize
       Sa_z(n) = 10.0_R8

       !--- correction to NCEP winds based on QSCAT ---
       uprime = Sa_u(n)*windFactor(n)
       vprime = Sa_v(n)*windFactor(n)
       Sa_u(n) = uprime*cos(winddFactor(n)*degtorad) - vprime*sin(winddFactor(n)*degtorad)
       Sa_v(n) = uprime*sin(winddFactor(n)*degtorad) + vprime*cos(winddFactor(n)*degtorad)

       !--- density and pslv taken directly from input stream, set pbot ---
       Sa_pbot(n) = Sa_pslv(n)

       !--- correction to NCEP Arctic & Antarctic air T & potential T ---
       if      ( yc(n) < -60.0_R8 ) then
          tMin = (avg_c0 + avg_c1*yc(n)) + (amp_c0 + amp_c1*yc(n))*cosFactor + tKFrz
          Sa_tbot(n) = max(Sa_tbot(n), tMin)
       else if ( yc(n) > 60.0_R8 ) then
          factor = MIN(1.0_R8, 0.1_R8*(yc(n)-60.0_R8) )
          Sa_tbot(n) = Sa_tbot(n) + factor * dTarc(target_mon)
       endif
       Sa_ptem(n) = Sa_tbot(n)

       !---  correction to NCEP relative humidity for heat budget balance ---
       Sa_shum(n) = Sa_shum(n) + qsatFactor(n)

       !--- Dupont correction to NCEP Arctic air T  ---
       !--- don't correct during summer months (July-September)
       !--- ONLY correct when forcing year is 1997->2004
       if (trim(datamode) == 'CORE2_IAF' ) then
          Sa_tbot(n) = Sa_tbot(n) +  strm_tarcf(n)
          Sa_ptem(n) = Sa_tbot(n)
       end if

       ! PRECIPITATION DATA
       strm_prec(n) = strm_prec(n)/86400.0_R8        ! convert mm/day to kg/m^2/s
       !  only correct satellite products, do not correct Serreze Arctic data
       if ( yc(n) < 58. ) then
          strm_prec(n) = strm_prec(n)*1.14168_R8
       endif
       if ( yc(n) >= 58. .and. yc(n) < 68. ) then
          factor = MAX(0.0_R8, 1.0_R8 - 0.1_R8*(yc(n)-58.0_R8) )
          strm_prec(n) = strm_prec(n)*(factor*(1.14168_R8 - 1.0_R8) + 1.0_R8)
       endif
       Faxa_rainc(n) = 0.0_R8               ! default zero
       Faxa_snowc(n) = 0.0_R8
       if (Sa_tbot(n) < tKFrz ) then        ! assign precip to rain/snow components
          Faxa_rainl(n) = 0.0_R8
          Faxa_snowl(n) = strm_prec(n)
       else
          Faxa_rainl(n) = strm_prec(n)
          Faxa_snowl(n) = 0.0_R8
       endif

       ! RADIATION DATA
       !--- fabricate required swdn components from net swdn ---
       Faxa_swvdr(n) = strm_swdn(n)*(0.28_R8)
       Faxa_swndr(n) = strm_swdn(n)*(0.31_R8)
       Faxa_swvdf(n) = strm_swdn(n)*(0.24_R8)
       Faxa_swndf(n) = strm_swdn(n)*(0.17_R8)
       !--- compute net short-wave based on LY08 latitudinally-varying albedo ---
       avg_alb = ( 0.069 - 0.011*cos(2.0_R8*yc(n)*degtorad ) )
       Faxa_swnet(n) = strm_swdn(n)*(1.0_R8 - avg_alb)
       !--- corrections to GISS sswdn for heat budget balancing ---
       factor = 1.0_R8
       if      ( -60.0_R8 < yc(n) .and. yc(n) < -50.0_R8 ) then
          factor = 1.0_R8 - (yc(n) + 60.0_R8)*(0.05_R8/10.0_R8)
       else if ( -50.0_R8 < yc(n) .and. yc(n) <  30.0_R8 ) then
          factor = 0.95_R8
       else if (  30.0_R8 < yc(n) .and. yc(n) <  40._R8 ) then
          factor = 1.0_R8 - (40.0_R8 - yc(n))*(0.05_R8/10.0_R8)
       endif
       Faxa_swnet(n) = Faxa_swnet(n)*factor
       Faxa_swvdr(n) = Faxa_swvdr(n)*factor
       Faxa_swndr(n) = Faxa_swndr(n)*factor
       Faxa_swvdf(n) = Faxa_swvdf(n)*factor
       Faxa_swndf(n) = Faxa_swndf(n)*factor
       !--- correction to GISS lwdn in Arctic ---
       if ( yc(n) > 60._R8 ) then
          factor = MIN(1.0_R8, 0.1_R8*(yc(n)-60.0_R8) )
          Faxa_lwdn(n) = Faxa_lwdn(n) + factor * dLWarc
       endif

    enddo   ! lsize

  end subroutine datm_datamode_core2_advance

  !===============================================================================
  subroutine datm_datamode_core2_restart_write(case_name, inst_suffix, ymd, tod, &
       logunit, mpicom, my_task, sdat)
    
    ! input/output variables
    character(len=*)            , intent(in)    :: case_name
    character(len=*)            , intent(in)    :: inst_suffix
    integer                     , intent(in)    :: ymd       ! model date
    integer                     , intent(in)    :: tod       ! model sec into model date
    integer                     , intent(in)    :: logunit
    integer                     , intent(in)    :: my_task
    integer                     , intent(in)    :: mpicom
    type(shr_strdata_type)      , intent(inout) :: sdat
    !-------------------------------------------------------------------------------

    call dshr_restart_write(rpfile, case_name, 'datm', inst_suffix, ymd, tod, &
         logunit, mpicom, my_task, sdat)

  end subroutine datm_datamode_core2_restart_write

  !===============================================================================
  subroutine datm_datamode_core2_restart_read(rest_filem, inst_suffix, logunit, my_task, mpicom, sdat)

    ! input/output arguments
    character(len=*)            , intent(inout) :: rest_filem
    character(len=*)            , intent(in)    :: inst_suffix
    integer                     , intent(in)    :: logunit
    integer                     , intent(in)    :: my_task
    integer                     , intent(in)    :: mpicom
    type(shr_strdata_type)      , intent(inout) :: sdat
    !-------------------------------------------------------------------------------

    call dshr_restart_read(rest_filem, rpfile, inst_suffix, nullstr, logunit, my_task, mpicom, sdat)

  end subroutine datm_datamode_core2_restart_read

  !===============================================================================
  subroutine datm_get_adjustment_factors(sdat, fileName_mesh, fileName_data, windF, winddF, qsatF, rc)

    ! input/output variables
    type(shr_strdata_type) , intent(in)  :: sdat 
    character(*)           , intent(in)  :: fileName_mesh ! file name string
    character(*)           , intent(in)  :: fileName_data ! file name string
    real(R8)               , pointer     :: windF(:)      ! wind adjustment factor
    real(R8)               , pointer     :: winddF(:)     ! wind adjustment factor
    real(r8)               , pointer     :: qsatF(:)      ! rel humidty adjustment factor
    integer                , intent(out) :: rc

    ! local variables
    type(ESMF_Mesh)        :: mesh      ! mesh read in from fileName_mesh
    type(ESMF_DistGrid)    :: distgrid
    type(ESMF_FieldBundle) :: fldbun_src
    type(ESMF_FieldBundle) :: fldbun_dst
    type(ESMF_RouteHandle) :: route_handle
    type(ESMF_Field)       :: field_src
    type(ESMF_Field)       :: field_dst
    integer                :: lsize
    integer                :: gsize
    integer                :: n         ! generic indicies
    logical                :: domap     ! map or not
    integer, pointer       :: gindex(:) ! domain decomposition of data
    integer                :: ndims     ! number of dims
    integer, allocatable   :: dimid(:)
    type(var_desc_t)       :: varid
    type(file_desc_t)      :: pioid
    type(io_desc_t)        :: pio_iodesc
    integer                :: rcode
    integer                :: nxg, nyg
    real(r8), pointer      :: data(:)
    integer                :: srcTermProcessing_Value = 0
    character(*) ,parameter :: subName =  '(datm_get_adjustment_factors) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! read in the factors mesh
    mesh = ESMF_MeshCreate(trim(filename_mesh), fileformat=ESMF_FILEFORMAT_ESMFMESH, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create input and output field bundles
    fldbun_src = ESMF_FieldBundleCreate(rc=rc) ! input field bundle
    field_src = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name='windFactor', &
         meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(fldbun_src, (/field_src/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_src = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name='winddFactor', &
         meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(fldbun_src, (/field_src/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_src = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name='qsatFactor', &
         meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(fldbun_src, (/field_src/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    fldbun_dst = ESMF_FieldBundleCreate(rc=rc) ! output field bundle
    field_dst = ESMF_FieldCreate(sdat%model_mesh, ESMF_TYPEKIND_R8, name='windFactor', &
         meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(fldbun_dst, (/field_dst/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_dst = ESMF_FieldCreate(sdat%model_mesh, ESMF_TYPEKIND_R8, name='winddFactor', &
         meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(fldbun_dst, (/field_dst/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    field_dst = ESMF_FieldCreate(sdat%model_mesh, ESMF_TYPEKIND_R8, name='qsatFactor', &
         meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call ESMF_FieldBundleAdd(fldbun_dst, (/field_dst/), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get mesh info
    call ESMF_MeshGet(mesh, elementdistGrid=distGrid, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_DistGridGet(distGrid, localDe=0, elementCount=lsize, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    allocate(gindex(lsize))
    call ESMF_DistGridGet(distGrid, localDe=0, seqIndexList=gindex, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create_pio_iodesc
    rcode = pio_openfile(sdat%pio_subsystem, pioid, sdat%io_type, trim(filename_data), pio_nowrite)
    call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
    rcode = pio_inq_varid(pioid, 'windFactor', varid)
    rcode = pio_inq_varndims(pioid, varid, ndims)
    allocate(dimid(ndims))
    rcode = pio_inq_varid(pioid, 'windFactor', varid)
    rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
    rcode = pio_inq_dimlen(pioid, dimid(1), nxg)
    rcode = pio_inq_dimlen(pioid, dimid(2), nyg)
    call pio_initdecomp(sdat%pio_subsystem, pio_double, (/nxg,nyg/), gindex, pio_iodesc)
    deallocate(gindex)

    ! Read in the data into the appropriate field bundle pointers
    call dshr_fldbun_getFldPtr(fldbun_src, 'windFactor', data, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    rcode = pio_inq_varid(pioid, 'windFactor', varid)
    call pio_read_darray(pioid, varid, pio_iodesc, data, rcode)

    call dshr_fldbun_getFldPtr(fldbun_src, 'winddFactor', data, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    rcode = pio_inq_varid(pioid, 'winddFactor', varid)
    call pio_read_darray(pioid, varid, pio_iodesc, data, rcode)

    call dshr_fldbun_getFldPtr(fldbun_src, 'qsatFactor', data, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    rcode = pio_inq_varid(pioid, 'qsatFactor', varid)
    call pio_read_darray(pioid, varid, pio_iodesc, data, rcode)

    if (nxg*nyg /= sdat%model_gsize) then
       ! TODO: this needs a mask that needs to be read in to have the mapping be accurate
       ! create bilinear route handle -
       call ESMF_FieldRegridStore(field_src, field_dst, routehandle=route_handle, &
            regridmethod=ESMF_REGRIDMETHOD_BILINEAR, &
            polemethod=ESMF_POLEMETHOD_ALLAVG, &
            extrapMethod=ESMF_EXTRAPMETHOD_NEAREST_STOD, &
            srcTermProcessing=srcTermProcessing_Value, &
            ignoreDegenerate=.true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_FieldRedistStore(field_src, field_dst, routehandle=route_handle, &
            ignoreUnmatchedIndices=.true., rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! apply the route handle
    call dshr_fldbun_regrid(fldbun_src, fldbun_dst, route_handle, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! obtain output data
    call dshr_fldbun_getFldPtr(fldbun_dst, 'windFactor', windF, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_fldbun_getFldPtr(fldbun_dst, 'winddFactor', winddF, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_fldbun_getFldPtr(fldbun_dst, 'qsatFactor', qsatF, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call pio_closefile(pioid)
    call pio_freedecomp(sdat%pio_subsystem, pio_iodesc)
    call ESMF_RouteHandleDestroy(route_handle, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

  end subroutine datm_get_adjustment_factors

end module datm_datamode_core2_mod
