module datm_comp_mod

  use NUOPC                 , only : NUOPC_Advertise
  use ESMF                  , only : ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_LogWrite
  use ESMF                  , only : ESMF_State, ESMF_StateGet, ESMF_StateItem_Flag, ESMF_STATEITEM_NOTFOUND
  use ESMF                  , only : ESMF_Mesh, ESMF_MeshGet
  use ESMF                  , only : operator(/=), operator(==)
  use perf_mod              , only : t_startf, t_stopf, t_adj_detailf, t_barrierf
  use shr_kind_mod          , only : r8=>shr_kind_r8, cxx=>shr_kind_cxx, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod           , only : shr_sys_abort
  use shr_mpi_mod           , only : shr_mpi_bcast, shr_mpi_max
  use shr_cal_mod           , only : shr_cal_date2julian
  use shr_precip_mod        , only : shr_precip_partition_rain_snow_ramp
  use shr_const_mod         , only : shr_const_spval, shr_const_tkfrz, shr_const_pi
  use shr_const_mod         , only : shr_const_pstd, shr_const_stebol, shr_const_rdair
  use shr_strdata_mod       , only : shr_strdata_type, shr_strdata_advance, shr_strdata_setOrbs
  use dshr_methods_mod      , only : chkerr, state_getfldptr
  use dshr_dfield_mod       , only : dfield_type, dshr_dfield_add, dshr_dfield_copy
  use dshr_fldlist_mod      , only : fldlist_type, dshr_fldlist_add, dshr_fldlist_realize
  use dshr_nuopc_mod        , only : dshr_get_griddata, dshr_set_griddata, dshr_get_atm_adjustment_factors

  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: datm_comp_advertise
  public :: datm_comp_realize
  public :: datm_comp_run

  private :: datm_esat  ! determine saturation vapor pressure

  !--------------------------------------------------------------------------
  ! Module data
  !--------------------------------------------------------------------------

  ! linked lists
  type(fldList_type) , pointer :: fldsImport => null()
  type(fldList_type) , pointer :: fldsExport => null()
  type(dfield_type)  , pointer :: dfields    => null()

  logical                      :: get_importdata
  logical                      :: flds_co2a
  logical                      :: flds_co2b
  logical                      :: flds_co2c
  logical                      :: flds_wiso
  logical                      :: flds_presaero
  character(CL)                :: factorFn              ! file containing correction factors
  real(R8)                     :: tbotmax               ! units detector
  real(R8)                     :: tdewmax               ! units detector
  real(R8)                     :: anidrmax              ! existance detector
  real(R8), pointer            :: yc(:)                 ! array of model latitudes
  real(R8), pointer            :: windFactor(:)
  real(R8), pointer            :: winddFactor(:)
  real(R8), pointer            :: qsatFactor(:)

  ! constants
  real(R8),parameter           :: tKFrz    = SHR_CONST_TKFRZ
  real(R8),parameter           :: degtorad = SHR_CONST_PI/180.0_R8
  real(R8),parameter           :: pstd     = SHR_CONST_PSTD     ! standard pressure ~ Pa
  real(R8),parameter           :: stebol   = SHR_CONST_STEBOL   ! Stefan-Boltzmann constant ~ W/m^2/K^4
  real(R8),parameter           :: rdair    = SHR_CONST_RDAIR    ! dry air gas constant   ~ J/K/kg
  real(R8),parameter           :: avg_c0 =  61.846_R8
  real(R8),parameter           :: avg_c1 =   1.107_R8
  real(R8),parameter           :: amp_c0 = -21.841_R8
  real(R8),parameter           :: amp_c1 =  -0.447_R8
  real(R8),parameter           :: phs_c0 =   0.298_R8
  real(R8),parameter           :: dLWarc =  -5.000_R8

  real(R8) :: dTarc(12)
  data   dTarc      / 0.49_R8, 0.06_R8,-0.73_R8,  -0.89_R8,-0.77_R8,-1.02_R8, &
                     -1.99_R8,-0.91_R8, 1.72_R8,   2.30_R8, 1.81_R8, 1.06_R8/

  ! export state data
  real(r8), pointer :: Sa_topo(:)           => null()
  real(r8), pointer :: Sa_z(:)              => null()
  real(r8), pointer :: Sa_u(:)              => null()
  real(r8), pointer :: Sa_v(:)              => null()
  real(r8), pointer :: Sa_tbot(:)           => null()
  real(r8), pointer :: Sa_ptem(:)           => null()
  real(r8), pointer :: Sa_shum(:)           => null()
  real(r8), pointer :: Sa_dens(:)           => null()
  real(r8), pointer :: Sa_pbot(:)           => null()
  real(r8), pointer :: Sa_pslv(:)           => null()
  real(r8), pointer :: Faxa_lwdn(:)         => null()
  real(r8), pointer :: Faxa_rainc(:)        => null()
  real(r8), pointer :: Faxa_rainl(:)        => null()
  real(r8), pointer :: Faxa_snowc(:)        => null()
  real(r8), pointer :: Faxa_snowl(:)        => null()
  real(r8), pointer :: Faxa_swndr(:)        => null()
  real(r8), pointer :: Faxa_swndf(:)        => null()
  real(r8), pointer :: Faxa_swvdr(:)        => null()
  real(r8), pointer :: Faxa_swvdf(:)        => null()
  real(r8), pointer :: Faxa_swnet(:)        => null()
  real(r8), pointer :: Sa_co2prog(:)        => null() ! co2
  real(r8), pointer :: Sa_co2diag(:)        => null() ! co2
  real(r8), pointer :: Faxa_bcph(:,:)       => null() ! prescribed aerosols
  real(r8), pointer :: Faxa_ocph(:,:)       => null() ! prescribed aerosols
  real(r8), pointer :: Faxa_dstwet(:,:)     => null() ! prescribed aerosols
  real(r8), pointer :: Faxa_dstdry(:,:)     => null() ! prescribed aerosols
  real(r8), pointer :: Faxa_rainc_wiso(:,:) => null() ! water isotopes
  real(r8), pointer :: Faxa_rainl_wiso(:,:) => null() ! water isotopes
  real(r8), pointer :: Faxa_snowc_wiso(:,:) => null() ! water isotopes
  real(r8), pointer :: Faxa_snowl_wiso(:,:) => null() ! water isotopes
  real(r8), pointer :: Sa_shum_wiso(:,:)    => null() ! water isotopes
  real(r8), pointer :: So_t(:)              => null() ! ocean temperature

  ! import state data
  ! Note - do not need to allocate memory for fields that are advertised if
  ! they are not used - this will permit IAF compsets to work and be written
  ! out to coupler history files and sent to datm - but not actually used
  real(r8), pointer :: Sx_avsdr(:) => null()
  real(r8), pointer :: Sx_anidr(:) => null()
  real(r8), pointer :: Sx_avsdf(:) => null()
  real(r8), pointer :: Sx_anidf(:) => null()

  ! stream internal data
  real(r8), pointer :: strm_z(:)         => null()
  real(r8), pointer :: strm_wind(:)      => null()
  real(r8), pointer :: strm_tdew(:)      => null()
  real(r8), pointer :: strm_tbot(:)      => null()
  real(r8), pointer :: strm_pbot(:)      => null()
  real(r8), pointer :: strm_shum(:)      => null()
  real(r8), pointer :: strm_lwdn(:)      => null()
  real(r8), pointer :: strm_rh(:)        => null()
  real(r8), pointer :: strm_swdn(:)      => null()
  real(r8), pointer :: strm_swdndf(:)    => null()
  real(r8), pointer :: strm_swdndr(:)    => null()
  real(r8), pointer :: strm_prec(:)      => null()
  real(r8), pointer :: strm_precc(:)     => null()
  real(r8), pointer :: strm_precl(:)     => null()
  real(r8), pointer :: strm_precn(:)     => null()
  real(r8), pointer :: strm_swup(:)      => null()
  real(r8), pointer :: strm_tarcf(:)     => null()
  real(r8), pointer :: strm_precsf(:)    => null()
  real(r8), pointer :: strm_rh_16O(:)    => null() ! water isoptopes
  real(r8), pointer :: strm_rh_18O(:)    => null() ! water isoptopes
  real(r8), pointer :: strm_rh_HDO(:)    => null() ! water isoptopes
  real(r8), pointer :: strm_precn_16O(:) => null() ! water isoptopes
  real(r8), pointer :: strm_precn_18O(:) => null() ! water isoptopes
  real(r8), pointer :: strm_precn_HDO(:) => null() ! water isoptopes
  real(r8), pointer :: strm_u_af(:)      => null() ! anomoly forcing
  real(r8), pointer :: strm_v_af(:)      => null() ! anomoly forcing
  real(r8), pointer :: strm_prec_af(:)   => null() ! anomoly forcing
  real(r8), pointer :: strm_tbot_af(:)   => null() ! anomoly forcing
  real(r8), pointer :: strm_pbot_af(:)   => null() ! anomoly forcing
  real(r8), pointer :: strm_shum_af(:)   => null() ! anomoly forcing
  real(r8), pointer :: strm_swdn_af(:)   => null() ! anomoly forcing
  real(r8), pointer :: strm_lwdn_af(:)   => null() ! anomoly forcing

  character(*),parameter       :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine datm_comp_advertise(importState, exportState, flds_scalar_name, &
       get_importdata_in, flds_presaero_in, flds_wiso_in, &
       flds_co2a_in, flds_co2b_in, flds_co2c_in, factorfn_in, rc)

    ! --------------------------------------------------------------
    ! determine export and import fields to advertise to mediator
    ! --------------------------------------------------------------

    ! input/output arguments
    type(ESMF_State)     , intent(inout) :: importState
    type(ESMF_State)     , intent(inout) :: exportState
    character(len=*)     , intent(in)    :: flds_scalar_name
    logical              , intent(in)    :: get_importdata_in
    logical              , intent(in)    :: flds_presaero_in
    logical              , intent(in)    :: flds_wiso_in
    logical              , intent(in)    :: flds_co2a_in
    logical              , intent(in)    :: flds_co2b_in
    logical              , intent(in)    :: flds_co2c_in
    character(len=*)     , intent(in)    :: factorfn_in
    integer              , intent(out)   :: rc

    ! local variables
    integer :: n
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    get_importdata = get_importdata_in
    flds_presaero  = flds_presaero_in
    flds_wiso      = flds_wiso_in
    flds_co2a      = flds_co2a_in
    flds_co2b      = flds_co2b_in
    flds_co2c      = flds_co2c_in
    factorfn       = factorfn_in

    !-------------------
    ! Advertise export fields
    !-------------------

    call dshr_fldList_add(fldsExport, trim(flds_scalar_name))
    call dshr_fldList_add(fldsExport, 'Sa_topo'    )
    call dshr_fldList_add(fldsExport, 'Sa_z'       )
    call dshr_fldList_add(fldsExport, 'Sa_u'       )
    call dshr_fldList_add(fldsExport, 'Sa_v'       )
    call dshr_fldList_add(fldsExport, 'Sa_ptem'    )
    call dshr_fldList_add(fldsExport, 'Sa_dens'    )
    call dshr_fldList_add(fldsExport, 'Sa_pslv'    )
    call dshr_fldList_add(fldsExport, 'Faxa_rainc' )
    call dshr_fldList_add(fldsExport, 'Faxa_rainl' )
    call dshr_fldList_add(fldsExport, 'Faxa_snowc' )
    call dshr_fldList_add(fldsExport, 'Faxa_snowl' )
    call dshr_fldList_add(fldsExport, 'Faxa_swndr' )
    call dshr_fldList_add(fldsExport, 'Faxa_swvdr' )
    call dshr_fldList_add(fldsExport, 'Faxa_swndf' )
    call dshr_fldList_add(fldsExport, 'Faxa_swvdf' )
    call dshr_fldList_add(fldsExport, 'Faxa_swnet' )
    call dshr_fldList_add(fldsExport, 'Sa_tbot'    )
    call dshr_fldList_add(fldsExport, 'Sa_pbot'    )
    call dshr_fldList_add(fldsExport, 'Sa_shum'    )
    call dshr_fldList_add(fldsExport, 'Faxa_lwdn'  )
    if (flds_co2a .or. flds_co2b .or. flds_co2c) then
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

    !-------------------
    ! advertise import fields
    !-------------------

    if (get_importdata) then
       call dshr_fldList_add(fldsImport, trim(flds_scalar_name))
       call dshr_fldList_add(fldsImport, "Sx_tref"       )
       call dshr_fldList_add(fldsImport, "Sx_qref"       )
       call dshr_fldList_add(fldsImport, "Sx_t"          )
       call dshr_fldList_add(fldsImport, "So_t"          )
       call dshr_fldList_add(fldsImport, "Sl_snowh"      )
       call dshr_fldList_add(fldsImport, "Sl_lfrac"      )
       call dshr_fldList_add(fldsImport, "Si_ifrac"      )
       call dshr_fldList_add(fldsImport, "So_ofrac"      )
       call dshr_fldList_add(fldsImport, "Faxx_taux"     )
       call dshr_fldList_add(fldsImport, "Faxx_tauy"     )
       call dshr_fldList_add(fldsImport, "Faxx_lat"      )
       call dshr_fldList_add(fldsImport, "Faxx_sen"      )
       call dshr_fldList_add(fldsImport, "Faxx_lwup"     )
       call dshr_fldList_add(fldsImport, "Faxx_evap"     )
     ! call dshr_fldList_add(fldsImport, "Fall_fco2_lnd" )
     ! call dshr_fldList_add(fldsImport, "Faoo_fco2_ocn" )

       fldlist => fldsImport ! the head of the linked list
       do while (associated(fldlist))
          call NUOPC_Advertise(importState, standardName=fldlist%stdname, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_LogWrite('(datm_comp_advertise): Fr_atm'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
          fldList => fldList%next
       enddo
    end if

  end subroutine datm_comp_advertise

  !===============================================================================

  subroutine datm_comp_realize(sdat, importState, exportState, flds_scalar_name, flds_scalar_num, mesh, &
       logunit, masterproc, rc)

    ! -----------------------------
    ! Initialize dfields arrays
    ! -----------------------------

    ! input/output parameters
    type(shr_strdata_type) , intent(inout) :: sdat
    type(ESMF_State)       , intent(inout) :: importState
    type(ESMF_State)       , intent(inout) :: exportState
    character(len=*)       , intent(in)    :: flds_scalar_name
    integer                , intent(in)    :: flds_scalar_num
    type(ESMF_Mesh)        , intent(in)    :: mesh
    integer                , intent(in)    :: logunit
    logical                , intent(in)    :: masterproc
    integer                , intent(out)   :: rc

    ! local variables
    integer  :: n, lsize, kf
    character(CS), allocatable :: strm_flds(:)
    character(*), parameter   :: subName = "(datm_comp_realize) "
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! -------------------------------------
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    ! -------------------------------------

    call dshr_fldlist_realize( exportState, fldsExport, flds_scalar_name, flds_scalar_num, mesh, &
         subname//':datmExport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_fldlist_realize( importState, fldsImport, flds_scalar_name, flds_scalar_num, mesh, &
         subname//':datmImport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !-----------------------------
    ! initialize dfields for export fields that have a corresponding stream field
    !-----------------------------

    call dshr_dfield_add(dfields, sdat, state_fld='Sa_topo'    , strm_fld='topo'  , &
         state=exportState, state_ptr=Sa_topo    , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat, state_fld='Sa_z'       , strm_fld='z'     , &
         state=exportState, state_ptr=Sa_z       , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat, state_fld='Sa_u'       , strm_fld='u'     , &
         state=exportState, state_ptr=Sa_u       , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Sa_v'       , strm_fld='v'     , &
         state=exportState, state_ptr=Sa_v       , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Sa_ptem'    , strm_fld='ptem'  , &
         state=exportState, state_ptr=Sa_ptem    , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Sa_dens'    , strm_fld='dens'  , &
         state=exportState, state_ptr=Sa_dens    , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Sa_pslv'    , strm_fld='pslv'  , &
         state=exportState, state_ptr=Sa_pslv    , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Faxa_rainc' , strm_fld='rainc' , &
         state=exportState, state_ptr=Faxa_rainc , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Faxa_rainl' , strm_fld='rainl' , &
         state=exportState, state_ptr=Faxa_rainl , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Faxa_snowc' , strm_fld='snowc' , &
         state=exportState, state_ptr=Faxa_snowc , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Faxa_snowl' , strm_fld='snowl' , &
         state=exportState, state_ptr=Faxa_snowl , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Faxa_swndr' , strm_fld='swndr' , &
         state=exportState, state_ptr=Faxa_swndr , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Faxa_swvdr' , strm_fld='swvdr' , &
         state=exportState, state_ptr=Faxa_swvdr , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Faxa_swndf' , strm_fld='swndf' , &
         state=exportState, state_ptr=Faxa_swndf , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Faxa_swvdf' , strm_fld='swvdf' , &
         state=exportState, state_ptr=Faxa_swvdf , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Faxa_swnet' , strm_fld='swnet' , &
         state=exportState, state_ptr=Faxa_swnet , logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (flds_presaero) then
       allocate(strm_flds(3))
       strm_flds = (/'bcphidry', 'bcphodry', 'bcphiwet'/)
       call dshr_dfield_add( dfields, sdat, state_fld='Faxa_bcph', strm_flds=strm_flds, &
            state=exportState, state_ptr=Faxa_bcph, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       strm_flds = (/'ocphidry', 'ocphodry', 'ocphiwet'/)
       call dshr_dfield_add( dfields, sdat, state_fld='Faxa_ocph', strm_flds=strm_flds, &
            state=exportState, state_ptr=Faxa_ocph, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       deallocate(strm_flds)
       allocate(strm_flds(4))
       strm_flds = (/'dstwet1', 'dstwet2', 'dstwet3', 'dstwet4'/)
       call dshr_dfield_add( dfields, sdat, state_fld='Faxa_dstwet', strm_flds=strm_flds, &
            state=exportState, state_ptr=Faxa_dstwet, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       strm_flds = (/'dstdry1', 'dstdry2', 'dstdry3', 'dstdry4'/)
       call dshr_dfield_add( dfields, sdat, state_fld='Faxa_dstdry', strm_flds=strm_flds, &
            state=exportState, state_ptr=Faxa_dstdry, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       deallocate(strm_flds)
    end if

    if (flds_wiso) then ! isopic forcing
       allocate(strm_flds(3))
       strm_flds = (/'rainc_16O', 'rainc_18O', 'rainc_HDO'/)
       call dshr_dfield_add( dfields, sdat, state_fld='Faxa_rainc_wiso', strm_flds=strm_flds, &
            state=exportState, state_ptr=Faxa_rainc_wiso, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
       strm_flds = (/'rainl_16O', 'rainl_18O', 'rainl_HDO'/)
       call dshr_dfield_add( dfields, sdat, state_fld='Faxa_rainl_wiso', strm_flds=strm_flds, &
            state=exportState, state_ptr=Faxa_rainl_wiso, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
       strm_flds = (/'snowc_16O', 'snowc_18O', 'snowc_HDO'/)
       call dshr_dfield_add( dfields, sdat, state_fld='Faxa_snowc_wiso', strm_flds=strm_flds, &
            state=exportState, state_ptr=Faxa_snowc_wiso, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
       strm_flds = (/'snowl_16O', 'snowl_18O', 'snowl_HDO'/)
       call dshr_dfield_add( dfields, sdat, state_fld='Faxa_snowl_wiso', strm_flds=strm_flds, &
            state=exportState, state_ptr=Faxa_snowl_wiso, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
       strm_flds = (/'shum_16O', 'shum_18O', 'shum_HDO'/)
       call dshr_dfield_add( dfields, sdat, state_fld='Sa_shum_wiso', strm_flds=strm_flds, &
            state=exportState, state_ptr=Sa_shum_wiso, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !-----------------------------
    ! initialize dfields for export fields that have a corresponding stream field AND that have a corresponding module field
    !-----------------------------

    call dshr_dfield_add( dfields, sdat, state_fld='Sa_tbot'   , strm_fld='tbot', &
         state=exportState, state_ptr=Sa_tbot, strm_ptr=strm_tbot, logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Sa_pbot'   , strm_fld='pbot', &
         state=exportState, state_ptr=Sa_pbot, strm_ptr=strm_pbot, logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Sa_shum'   , strm_fld='shum', &
         state=exportState, state_ptr=Sa_shum, strm_ptr=strm_shum, logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add( dfields, sdat, state_fld='Faxa_lwdn' , strm_fld='lwdn', &
         state=exportState, state_ptr=Faxa_lwdn, strm_ptr=strm_lwdn, logunit=logunit, masterproc=masterproc, rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    if (flds_co2a .or. flds_co2b .or. flds_co2c) then
       call dshr_dfield_add( dfields, sdat, state_fld='Sa_co2prog', strm_fld='co2prog', &
            state=exportState, state_ptr=Sa_co2prog, logunit=logunit, masterproc=masterproc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call dshr_dfield_add( dfields, sdat, state_fld='Sa_co2diag', strm_fld='co2diag', &
            state=exportState, state_ptr=Sa_co2diag, logunit=logunit, masterproc=masterproc, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !-----------------------------
    ! initialize dfields for stream fields that have no corresponding import or export fields
    !-----------------------------

    call dshr_dfield_add(dfields, sdat, 'wind'  , strm_wind   , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'tdew'  , strm_tdew   , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'tbot'  , strm_tbot   , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'pbot'  , strm_pbot   , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'shum'  , strm_shum   , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'lwdn'  , strm_lwdn   , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'rh'    , strm_rh     , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'swdn'  , strm_swdn   , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'swdndf', strm_swdndf , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'swdndr', strm_swdndr , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'prec'  , strm_prec   , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'precc' , strm_precc  , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'precl' , strm_precl  , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'precn' , strm_precn  , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'swup'  , strm_swup   , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'tarcf' , strm_tarcf  , logunit=logunit, masterproc=masterproc)

    ! water isotopes
    call dshr_dfield_add(dfields, sdat, 'rh_16O'   , strm_rh_16O    , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'rh_18O'   , strm_rh_18O    , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'rh_HDO'   , strm_rh_HDO    , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'precn_16O', strm_precn_16O , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'precn_18O', strm_precn_18O , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'precn_HDO', strm_precn_HDO , logunit=logunit, masterproc=masterproc)

    ! values for optionalcorrection / anomaly forcing (add Sa_precsf for precip scale factor)
    call dshr_dfield_add(dfields, sdat, 'precsf'   , strm_precsf    , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'prec_af'  , strm_prec_af   , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'u_af'     , strm_u_af      , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'v_af'     , strm_v_af      , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'tbot_af'  , strm_tbot_af   , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'pbot_af'  , strm_pbot_af   , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'shum_af'  , strm_shum_af   , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'swdn_af'  , strm_swdn_af   , logunit=logunit, masterproc=masterproc)
    call dshr_dfield_add(dfields, sdat, 'lwdn_af'  , strm_lwdn_af   , logunit=logunit, masterproc=masterproc)

  end subroutine datm_comp_realize

  !===============================================================================

  subroutine datm_comp_run(mpicom, compid, masterproc, logunit, target_ymd, target_tod, sdat, &
       mesh, target_mon, orbEccen, orbMvelpp, orbLambm0, orbObliqr, idt, rc)

    ! ----------------------------------
    ! run method for datm model
    ! ----------------------------------

    ! input/output variables
    integer                , intent(in)    :: mpicom           ! mpi communicator
    integer                , intent(in)    :: compid           ! mct comp id
    logical                , intent(in)    :: masterproc       ! true => master task
    integer                , intent(in)    :: logunit          ! logging unit number
    integer                , intent(in)    :: target_ymd       ! model date
    integer                , intent(in)    :: target_tod       ! model sec into model date
    type(shr_strdata_type) , intent(inout) :: sdat
    type(ESMF_Mesh)        , intent(in)    :: mesh
    integer                , intent(in)    :: target_mon       ! model month
    integer                , intent(in)    :: idt              ! integer model time step
    real(R8)               , intent(in)    :: orbEccen         ! orb eccentricity (unit-less)
    real(R8)               , intent(in)    :: orbMvelpp        ! orb moving vernal eq (radians)
    real(R8)               , intent(in)    :: orbLambm0        ! orb mean long of perhelion (radians)
    real(R8)               , intent(in)    :: orbObliqr        ! orb obliquity (radians)
    integer                , intent(out)   :: rc

    ! local variables
    integer                 :: n,kf                ! indices
    integer                 :: lsize               ! size of attr vect
    integer                 :: eday                ! elapsed day
    real(R8)                :: rday                ! elapsed day
    real(R8)                :: cosFactor           ! cosine factor
    real(R8)                :: factor              ! generic/temporary correction factor
    real(R8)                :: avg_alb             ! average albedo
    real(R8)                :: tMin                ! minimum temperature
    integer                 :: spatialDim          ! number of dimension in mesh
    integer                 :: numOwnedElements    ! size of mesh
    real(r8), pointer       :: ownedElemCoords(:)  ! mesh lat and lons
    real(R8)                :: uprime,vprime
    real(r8)                :: swndr,swndf,swvdr,swvdf,ratio_rvrf
    real(R8)                :: tbot,pbot,rtmp,vp,ea,e,qsat,frac,qsatT
    logical                 :: first_time = .true. ! first call logical
    character(*), parameter :: subName = '(datm_comp_run) '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('DATM_RUN')

    !--------------------
    ! Advance datm streams
    !--------------------

    ! set data needed for cosz t-interp method
    call shr_strdata_setOrbs(sdat, orbEccen, orbMvelpp, orbLambm0, orbObliqr, idt)

    ! time and spatially interpolate to model time and grid
    call t_barrierf('datm_BARRIER',mpicom)
    call t_startf('datm_strdata_advance')
    call shr_strdata_advance(sdat, target_ymd, target_tod, mpicom, 'datm')
    call t_stopf('datm_strdata_advance')

    !--------------------
    ! copy all fields from streams to export state as default
    !--------------------

    ! This automatically will update the fields in the export state
    call t_barrierf('datm_comp_dfield_copy_BARRIER', mpicom)
    call t_startf('datm_dfield_copy')
    call dshr_dfield_copy(dfields,  sdat, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('datm_dfield_copy')

    !-------------------------------------------------
    ! Determine data model behavior based on the mode
    !-------------------------------------------------

    lsize = size(Sa_z)

    if (first_time) then
       ! overwrite mask and frac
       call dshr_set_griddata(sdat, 'mask', rvalue=1.0_r8)
       call dshr_set_griddata(sdat, 'frac', rvalue=1.0_r8)

       ! allocate module arrays
       allocate(windFactor(lsize))
       allocate(winddFactor(lsize))
       allocate(qsatFactor(lsize))

       call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       allocate(ownedElemCoords(spatialDim*numOwnedElements))
       allocate(yc(numOwnedElements))
       call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       do n = 1,numOwnedElements
          yc(n) = ownedElemCoords(2*n)
       end do
    end if

    call t_startf('datm_datamode')
    select case (trim(sdat%datamode))

    case('COPYALL')
       ! do nothing extra

    case('CORE2_NYF','CORE2_IAF')
       if (first_time) then
          if (.not. associated(strm_prec) .or. .not. associated(strm_swdn)) then
             call shr_sys_abort(trim(subname)//'ERROR: prec and swdn must be in streams for CORE2')
          endif
          if (trim(sdat%datamode) == 'CORE2_IAF' ) then
             if (.not. associated(strm_tarcf)) then
                call shr_sys_abort(trim(subname)//'tarcf must be in an input stream for CORE2_IAF')
             endif
          endif
          call dshr_get_atm_adjustment_factors(factorFn, windFactor, winddFactor, qsatFactor,  &
               mpicom, compid, masterproc, logunit, sdat)
       endif
       call shr_cal_date2julian(target_ymd, target_tod, rday, sdat%calendar)
       rday = mod((rday - 1.0_R8),365.0_R8)
       cosfactor = cos((2.0_R8*SHR_CONST_PI*rday)/365 - phs_c0)

       do n = 1,lsize
          Sa_z(n) = 10.0_R8

          !--- correction to NCEP winds based on QSCAT ---
          uprime = Sa_u(n)*windFactor(n)
          vprime = Sa_v(n)*windFactor(n)
          Sa_u(n) = uprime*cos(winddFactor(n)*degtorad) - vprime*sin(winddFactor(n)*degtorad)
          Sa_v(n) = uprime*sin(winddFactor(n)*degtorad) + vprime*cos(winddFactor(n)*degtorad)

          !--- density, tbot, & pslv taken directly from input stream, set pbot ---
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
          if (trim(sdat%datamode) == 'CORE2_IAF' ) then
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

    case('CORE_IAF_JRA')
       if (first_time) then
          if (.not. associated(strm_prec) .or. .not. associated(strm_swdn)) then
             call shr_sys_abort(trim(subname)//'ERROR: prec and swdn must be in streams for CORE_IAF_JRA')
          endif
          if (trim(sdat%datamode) == 'CORE_IAF_JRA' ) then
             if (.not. associated(strm_tarcf)) then
                call shr_sys_abort(trim(subname)//'ERROR: tarcf must be in an input stream for CORE_IAF_JRA')
             endif
          endif
          if (trim(factorFn) == 'null') then
            windFactor  = 1.0_R8
            winddFactor = 1.0_R8
            qsatFactor  = 1.0_R8
          else
             call dshr_get_atm_adjustment_factors(factorFn, windFactor, winddFactor, qsatFactor,  &
                  mpicom, compid, masterproc, logunit, sdat)
          endif
       endif
       call shr_cal_date2julian(target_ymd, target_tod, rday, sdat%calendar)
       rday = mod((rday - 1.0_R8),365.0_R8)
       cosfactor = cos((2.0_R8*SHR_CONST_PI*rday)/365 - phs_c0)

       do n = 1,lsize
          Sa_z(n) = 10.0_R8
          Sa_pbot(n) = Sa_pslv(n)
          Sa_ptem(n) = Sa_tbot(n)

          !--- density computation for JRA55 forcing ---
          Sa_dens(n) = Sa_pbot(n)/(rdair*Sa_tbot(n)*(1+0.608* Sa_shum(n)))

          ! PRECIPITATION DATA
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

       enddo   ! lsize

    case('CLMNCEP')
       if (first_time) then
          if (.not. associated(strm_wind) .or. .not. associated(strm_tbot)) then
             call shr_sys_abort(trim(subname)//' ERROR: wind and tbot must be in streams for CLMNCEP')
          endif
          rtmp = maxval(Sa_tbot(:))
          call shr_mpi_max(rtmp,tbotmax,mpicom,'datm_tbot',all=.true.)
          if (get_importdata) then
             rtmp = maxval(Sx_anidr(:))
             call shr_mpi_max(rtmp,anidrmax,mpicom,'datm_ani',all=.true.)
          else
             anidrmax = SHR_CONST_SPVAL ! see below for use
          end if
          if (associated(strm_tdew)) then
             rtmp = maxval(strm_tdew(:))
             call shr_mpi_max(rtmp,tdewmax,mpicom,'datm_tdew',all=.true.)
          endif
          if (masterproc) write(logunit,*) trim(subname),' max values = ',tbotmax,tdewmax,anidrmax
       endif
       do n = 1,lsize
          !--- bottom layer height ---
          if (.not. associated(strm_z)) Sa_z(n) = 30.0_R8

          !--- temperature ---
          if (tbotmax < 50.0_R8) Sa_tbot(n) = Sa_tbot(n) + tkFrz
          ! Limit very cold forcing to 180K
          Sa_tbot(n) = max(180._r8, Sa_tbot(n))
          Sa_ptem(n) = Sa_tbot(n)

          !--- pressure ---
          if (.not. associated(strm_pbot)) Sa_pbot(n) = pstd
          Sa_pslv(n) = Sa_pbot(n)

          !--- u, v wind velocity ---
          Sa_u(n) = strm_wind(n)/sqrt(2.0_R8)
          Sa_v(n) = Sa_u(n)

          !--- specific humidity ---
          tbot = Sa_tbot(n)
          pbot = Sa_pbot(n)
          if (associated(strm_shum)) then
             e = datm_esat(tbot,tbot)
             qsat = (0.622_R8 * e)/(pbot - 0.378_R8 * e)
             if (qsat < Sa_shum(n)) then
                Sa_shum(n) = qsat
             endif
          else if (associated(strm_rh)) then
             e = strm_rh(n) * 0.01_R8 * datm_esat(tbot,tbot)
             qsat = (0.622_R8 * e)/(pbot - 0.378_R8 * e)
             Sa_shum(n) = qsat
             if (flds_wiso) then
                ! for isotopic tracer specific humidity, expect a delta, so
                ! just keep the delta from the input file
                Sa_shum_wiso(1,n) = strm_rh_16O(n)
                Sa_shum_wiso(2,n) = strm_rh_18O(n)
                Sa_shum_wiso(3,n) = strm_rh_HDO(n)
             end if
          else if (associated(strm_tdew)) then
             if (tdewmax < 50.0_R8) strm_tdew(n) = strm_tdew(n) + tkFrz
             e = datm_esat(strm_tdew(n),tbot)
             qsat = (0.622_R8 * e)/(pbot - 0.378_R8 * e)
             Sa_shum(n) = qsat
          else
             call shr_sys_abort(subname//'ERROR: cannot compute shum')
          endif

          !--- density ---
          vp = (Sa_shum(n)*pbot) / (0.622_R8 + 0.378_R8 * Sa_shum(n))
          Sa_dens(n) = (pbot - 0.378_R8 * vp) / (tbot*rdair)

          !--- downward longwave ---
          if (.not. associated(strm_lwdn)) then
             e  = Sa_pslv(n) * Sa_shum(n) / (0.622_R8 + 0.378_R8 * Sa_shum(n))
             ea = 0.70_R8 + 5.95e-05_R8 * 0.01_R8 * e * exp(1500.0_R8/tbot)
             Faxa_lwdn(n) = ea * stebol * tbot**4
          endif

          !--- shortwave radiation ---
          if (associated(strm_swdndf) .and. associated(strm_swdndr)) then
             Faxa_swndr(n) = strm_swdndr(n) * 0.50_R8
             Faxa_swvdr(n) = strm_swdndr(n) * 0.50_R8
             Faxa_swndf(n) = strm_swdndf(n) * 0.50_R8
             Faxa_swvdf(n) = strm_swdndf(n) * 0.50_R8
          elseif (associated(strm_swdn)) then
             ! relationship between incoming NIR or VIS radiation and ratio of
             ! direct to diffuse radiation calculated based on one year's worth of
             ! hourly CAM output from CAM version cam3_5_55
             swndr = strm_swdn(n) * 0.50_R8
             ratio_rvrf =  min(0.99_R8,max(0.29548_R8 + 0.00504_R8*swndr  &
                  -1.4957e-05_R8*swndr**2 + 1.4881e-08_R8*swndr**3,0.01_R8))
             Faxa_swndr(n) = ratio_rvrf*swndr
             swndf = strm_swdn(n) * 0.50_R8
             Faxa_swndf(n) = (1._R8 - ratio_rvrf)*swndf

             swvdr = strm_swdn(n) * 0.50_R8
             ratio_rvrf =   min(0.99_R8,max(0.17639_R8 + 0.00380_R8*swvdr  &
                  -9.0039e-06_R8*swvdr**2 + 8.1351e-09_R8*swvdr**3,0.01_R8))
             Faxa_swvdr(n) = ratio_rvrf*swvdr
             swvdf = strm_swdn(n) * 0.50_R8
             Faxa_swvdf(n) = (1._R8 - ratio_rvrf)*swvdf
          else
             call shr_sys_abort(subName//'ERROR: cannot compute short-wave down')
          endif

          !--- swnet: a diagnostic quantity ---
          if (anidrmax < 1.0e-8 .or. anidrmax > SHR_CONST_SPVAL * 0.9_R8) then
             Faxa_swnet(n) = 0.0_R8
          else if ( associated(Sx_anidr) .and. associated(Sx_anidf) .and. &
                    associated(Sx_avsdr) .and. associated(Sx_avsdf)) then
             Faxa_swnet(n) = (1.0_R8-Sx_anidr(n))*Faxa_swndr(n) + &
                             (1.0_R8-Sx_avsdr(n))*Faxa_swvdr(n) + &
                             (1.0_R8-Sx_anidf(n))*Faxa_swndf(n) + &
                             (1.0_R8-Sx_avsdf(n))*Faxa_swvdf(n)
          else
             Faxa_swnet(n) = Faxa_swndr(n) + Faxa_swvdr(n) + Faxa_swndf(n) + Faxa_swvdf(n)
          endif

          !--- rain and snow ---
          if (associated(strm_precc) .and. associated(strm_precl)) then
             Faxa_rainc(n) = strm_precc(n)
             Faxa_rainl(n) = strm_precl(n)
          else if (associated(strm_precn)) then
             Faxa_rainc(n) = strm_precn(n)*0.1_R8
             Faxa_rainl(n) = strm_precn(n)*0.9_R8
          else
             call shr_sys_abort(subName//'ERROR: cannot compute rain and snow')
          endif

          !--- split precip between rain & snow ---
          call shr_precip_partition_rain_snow_ramp(tbot, frac)
          Faxa_snowc(n) = max(0.0_R8, Faxa_rainc(n)*(1.0_R8 - frac))
          Faxa_snowl(n) = max(0.0_R8, Faxa_rainl(n)*(1.0_R8 - frac))
          Faxa_rainc(n) = max(0.0_R8, Faxa_rainc(n)*(         frac))
          Faxa_rainl(n) = max(0.0_R8, Faxa_rainl(n)*(         frac))

       end do

    end select

    !----------------------------------------------------------
    ! bias correction / anomaly forcing ( start block )
    ! modify atmospheric input fields if streams exist
    !----------------------------------------------------------

    ! bias correct precipitation relative to observed
    ! (via bias_correct nameslist option)
    if (associated(strm_precsf)) then
       Faxa_snowc(:) = Faxa_snowc(:) * min(1.e2_r8,strm_precsf(:))
       Faxa_snowl(:) = Faxa_snowl(:) * min(1.e2_r8,strm_precsf(:))
       Faxa_rainc(:) = Faxa_rainc(:) * min(1.e2_r8,strm_precsf(:))
       Faxa_rainl(:) = Faxa_rainl(:) * min(1.e2_r8,strm_precsf(:))
    endif

    ! adjust atmospheric input fields if anomaly forcing streams exist
    ! (via anomaly_forcing namelist option)

    ! wind
    if (associated(strm_u_af) .and. associated(strm_v_af)) then
       Sa_u(:) = Sa_u(:) + strm_u_af(:)
       Sa_v(:) = Sa_v(:) + strm_v_af(:)
    endif

    ! specific humidity
    if (associated(strm_shum_af)) then
       Sa_shum(:) = Sa_shum(:) + strm_shum_af(:)
       ! avoid possible negative q values
       where (Sa_shum < 0._r8)
          Sa_shum = 1.e-6_r8
       end where
    endif

    ! pressure
    if (associated(strm_pbot_af)) then
       Sa_pbot(:) = Sa_pbot(:) + strm_pbot_af(:)
    endif

    ! temperature
    if (associated(strm_tbot_af)) then
       Sa_tbot(:) = Sa_tbot(:) + strm_tbot_af(:)
    endif

    ! longwave
    if (associated(strm_lwdn_af)) then
       Faxa_lwdn(:) = Faxa_lwdn(:) * strm_lwdn_af(:)
    endif

    ! precipitation
    if (associated(strm_prec_af)) then
       Faxa_snowc(:) = Faxa_snowc(:) * strm_prec_af(:)
       Faxa_snowl(:) = Faxa_snowl(:) * strm_prec_af(:)
       Faxa_rainc(:) = Faxa_rainc(:) * strm_prec_af(:)
       Faxa_rainl(:) = Faxa_rainl(:) * strm_prec_af(:)
    end if

    ! shortwave
    if (associated(strm_swdn_af)) then
       Faxa_swndr(:) = Faxa_swndr(:) * strm_swdn_af(:)
       Faxa_swvdr(:) = Faxa_swvdr(:) * strm_swdn_af(:)
       Faxa_swndf(:) = Faxa_swndf(:) * strm_swdn_af(:)
       Faxa_swvdf(:) = Faxa_swvdf(:) * strm_swdn_af(:)
    endif
    ! bias correction / anomaly forcing ( end block )


    ! Log output for model date
    if (masterproc) write(logunit,*) 'atm : model date ', target_ymd, target_tod
    first_time = .false.

    call t_stopf('datm_datamode')
    call t_stopf('DATM_RUN')

  end subroutine datm_comp_run

!===============================================================================

  real(R8) function datm_eSat(tK,tKbot)

    !----------------------------------------------------------------------------
    ! use polynomials to calculate saturation vapor pressure and derivative with
    ! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
    ! required to convert relative humidity to specific humidity
    !----------------------------------------------------------------------------

    ! input/output variables
    real(R8),intent(in) :: tK    ! temp used in polynomial calculation
    real(R8),intent(in) :: tKbot ! bottom atm temp

    ! local variables
    real(R8)           :: t     ! tK converted to Celcius
    real(R8),parameter :: tkFrz = shr_const_tkfrz  ! freezing T of fresh water ~ K

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

    t = min( 50.0_R8, max(-50.0_R8,(tK-tKfrz)) )
    if ( tKbot < tKfrz) then
       datm_eSat = 100.0_R8*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
    else
       datm_eSat = 100.0_R8*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
    end if

  end function datm_eSat

end module datm_comp_mod
