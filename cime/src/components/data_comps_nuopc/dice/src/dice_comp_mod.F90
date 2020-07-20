module dice_comp_mod

  use NUOPC                 , only : NUOPC_Advertise
  use ESMF                  , only : ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_LogWrite
  use ESMF                  , only : ESMF_State, ESMF_StateGet, ESMF_StateItem_Flag, ESMF_STATEITEM_NOTFOUND
  use ESMF                  , only : ESMF_Mesh, ESMF_MeshGet
  use ESMF                  , only : operator(/=), operator(==)
  use perf_mod              , only : t_startf, t_stopf, t_adj_detailf, t_barrierf
  use shr_kind_mod          , only : r8=>shr_kind_r8, cxx=>shr_kind_cxx, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod           , only : shr_sys_abort
  use shr_const_mod         , only : shr_const_pi, shr_const_spval, shr_const_tkfrz, shr_const_latice
  use shr_frz_mod           , only : shr_frz_freezetemp
  use shr_cal_mod           , only : shr_cal_ymd2julian
  use shr_strdata_mod       , only : shr_strdata_type, shr_strdata_advance
  use dshr_methods_mod      , only : chkerr, state_getfldptr
  use dshr_dfield_mod       , only : dfield_type, dshr_dfield_add, dshr_dfield_copy
  use dshr_fldlist_mod      , only : fldlist_type, dshr_fldlist_add, dshr_fldlist_realize
  use dshr_nuopc_mod        , only : dshr_get_griddata
  use dice_flux_atmice_mod  , only : dice_flux_atmice

  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public/Private interfaces
  !--------------------------------------------------------------------------

  public :: dice_comp_advertise
  public :: dice_comp_realize
  public :: dice_comp_run

  !--------------------------------------------------------------------------
  ! Module data
  !--------------------------------------------------------------------------

  ! linked lists
  type(fldList_type) , pointer :: fldsImport => null()
  type(fldList_type) , pointer :: fldsExport => null()
  type(dfield_type)  , pointer :: dfields    => null()

  real(r8),parameter         :: pi     = shr_const_pi      ! pi
  real(r8),parameter         :: spval  = shr_const_spval   ! flags invalid data
  real(r8),parameter         :: tFrz   = shr_const_tkfrz   ! temp of freezing
  real(r8),parameter         :: latice = shr_const_latice  ! latent heat of fusion
  real(r8),parameter         :: waterMax = 1000.0_r8       ! wrt iFrac comp & frazil ice (kg/m^2)

  ! surface albedo constants
  real(r8),parameter         :: snwfrac = 0.286_r8 ! snow cover fraction ~ [0,1]
  real(r8),parameter         :: as_nidf = 0.950_r8 ! albedo: snow,near-infr,diffuse
  real(r8),parameter         :: as_vsdf = 0.700_r8 ! albedo: snow,visible  ,diffuse
  real(r8),parameter         :: as_nidr = 0.960_r8 ! albedo: snow,near-infr,direct
  real(r8),parameter         :: as_vsdr = 0.800_r8 ! albedo: snow,visible  ,direct
  real(r8),parameter         :: ai_nidf = 0.700_r8 ! albedo: ice, near-infr,diffuse
  real(r8),parameter         :: ai_vsdf = 0.500_r8 ! albedo: ice, visible  ,diffuse
  real(r8),parameter         :: ai_nidr = 0.700_r8 ! albedo: ice, near-infr,direct
  real(r8),parameter         :: ai_vsdr = 0.500_r8 ! albedo: ice, visible  ,direct
  real(r8),parameter         :: ax_nidf = ai_nidf*(1.0_r8-snwfrac) + as_nidf*snwfrac
  real(r8),parameter         :: ax_vsdf = ai_vsdf*(1.0_r8-snwfrac) + as_vsdf*snwfrac
  real(r8),parameter         :: ax_nidr = ai_nidr*(1.0_r8-snwfrac) + as_nidr*snwfrac
  real(r8),parameter         :: ax_vsdr = ai_vsdr*(1.0_r8-snwfrac) + as_vsdr*snwfrac

  ! restart fields
  real(r8), pointer, public :: water(:) => null()

  ! internal fields
  real(r8), pointer :: yc(:)      => null() ! mesh lats (degrees)
  integer , pointer :: imask(:)   => null()
  real(r8), pointer :: tfreeze(:) => null()
  !real(r8), pointer:: ifrac0(:)  => null()

  ! export fields
  real(r8), pointer ::  Si_imask(:)      => null()
  real(r8), pointer ::  Si_ifrac(:)      => null()
  real(r8), pointer ::  Si_t(:)          => null()
  real(r8), pointer ::  Si_tref(:)       => null()
  real(r8), pointer ::  Si_qref(:)       => null()
  real(r8), pointer ::  Si_avsdr(:)      => null()
  real(r8), pointer ::  Si_anidr(:)      => null()
  real(r8), pointer ::  Si_avsdf(:)      => null()
  real(r8), pointer ::  Si_anidf(:)      => null()
  real(r8), pointer ::  Faii_swnet(:)    => null()
  real(r8), pointer ::  Faii_sen(:)      => null()
  real(r8), pointer ::  Faii_lat(:)      => null()
  real(r8), pointer ::  Faii_lwup(:)     => null()
  real(r8), pointer ::  Faii_evap(:)     => null()
  real(r8), pointer ::  Faii_taux(:)     => null()
  real(r8), pointer ::  Faii_tauy(:)     => null()
  real(r8), pointer ::  Fioi_melth(:)    => null()
  real(r8), pointer ::  Fioi_meltw(:)    => null()
  real(r8), pointer ::  Fioi_swpen(:)    => null()
  real(r8), pointer ::  Fioi_taux(:)     => null()
  real(r8), pointer ::  Fioi_tauy(:)     => null()
  real(r8), pointer ::  Fioi_salt(:)     => null()
  real(r8), pointer ::  Fioi_bcpho(:)    => null()
  real(r8), pointer ::  Fioi_bcphi(:)    => null()
  real(r8), pointer ::  Fioi_flxdst(:)   => null()
  real(r8), pointer ::  Si_ifrac_n(:,:)  => null()
  real(r8), pointer ::  Fioi_swpen_ifrac_n(:,:) => null()

  ! import fields
  real(r8), pointer :: Faxa_swvdr(:)    => null()
  real(r8), pointer :: Faxa_swvdf(:)    => null()
  real(r8), pointer :: Faxa_swndr(:)    => null()
  real(r8), pointer :: Faxa_swndf(:)    => null()
  real(r8), pointer :: Fioo_q(:)        => null()
  real(r8), pointer :: Sa_z(:)          => null()
  real(r8), pointer :: Sa_u(:)          => null()
  real(r8), pointer :: Sa_v(:)          => null()
  real(r8), pointer :: Sa_ptem(:)       => null()
  real(r8), pointer :: Sa_shum(:)       => null()
  real(r8), pointer :: Sa_dens(:)       => null()
  real(r8), pointer :: Sa_tbot(:)       => null()
  real(r8), pointer :: So_s(:)          => null()
  real(r8), pointer :: Faxa_bcph(:,:)   => null()
  real(r8), pointer :: Faxa_ocph(:,:)   => null()
  real(r8), pointer :: Faxa_dstdry(:,:) => null()
  real(r8), pointer :: Faxa_dstwet(:,:) => null()

  logical :: flds_i2o_per_cat

  character(*) , parameter    :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine dice_comp_advertise(importState, exportState, flds_scalar_name, flds_i2o_per_cat_in, rc)

    ! --------------------------------------------------------------
    ! determine export and import fields to advertise to mediator
    ! --------------------------------------------------------------

    ! input/output arguments
    type(ESMF_State)     , intent(inout) :: importState
    type(ESMF_State)     , intent(inout) :: exportState
    character(len=*)     , intent(in)    :: flds_scalar_name
    logical              , intent(in)    :: flds_i2o_per_cat_in
    integer              , intent(out)   :: rc

    ! local variables
    integer         :: n
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    ! Save as module variables
    flds_i2o_per_cat = flds_i2o_per_cat_in

    ! Advertise export fields

    call dshr_fldList_add(fldsExport , trim(flds_scalar_name))
    call dshr_fldList_add(fldsExport ,'Si_ifrac'    )
    call dshr_fldList_add(fldsExport ,'Si_imask'    )
    call dshr_fldList_add(fldsExport ,'Si_t'        )
    call dshr_fldList_add(fldsExport ,'Si_tref'     )
    call dshr_fldList_add(fldsExport ,'Si_qref'     )
    call dshr_fldList_add(fldsExport ,'Si_avsdr'    )
    call dshr_fldList_add(fldsExport ,'Si_anidr'    )
    call dshr_fldList_add(fldsExport ,'Si_avsdf'    )
    call dshr_fldList_add(fldsExport ,'Si_anidf'    )
    call dshr_fldList_add(fldsExport ,'Faii_swnet'  )
    call dshr_fldList_add(fldsExport ,'Faii_sen'    )
    call dshr_fldList_add(fldsExport ,'Faii_lat'    )
    call dshr_fldList_add(fldsExport ,'Faii_lwup'   )
    call dshr_fldList_add(fldsExport ,'Faii_evap'   )
    call dshr_fldList_add(fldsExport ,'Faii_taux'   )
    call dshr_fldList_add(fldsExport ,'Faii_tauy'   )
    call dshr_fldList_add(fldsExport ,'Fioi_melth'  )
    call dshr_fldList_add(fldsExport ,'Fioi_meltw'  )
    call dshr_fldList_add(fldsExport ,'Fioi_swpen'  )
    call dshr_fldList_add(fldsExport ,'Fioi_taux'   )
    call dshr_fldList_add(fldsExport ,'Fioi_tauy'   )
    call dshr_fldList_add(fldsExport ,'Fioi_salt'   )
    call dshr_fldList_add(fldsExport ,'Fioi_bcpho'  )
    call dshr_fldList_add(fldsExport ,'Fioi_bcphi'  )
    call dshr_fldList_add(fldsExport ,'Fioi_flxdst' )
    if (flds_i2o_per_cat) then
       call dshr_fldList_add(fldsExport, 'Si_ifrac_n'        , ungridded_lbound=1, ungridded_ubound=1)
       call dshr_fldList_add(fldsExport, 'Fioi_swpen_ifrac_n', ungridded_lbound=1, ungridded_ubound=1)
    end if

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(dice_comp_advertise): Fr_ice'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

    ! Advertise import fields

    call dshr_fldList_add(fldsImport , trim(flds_scalar_name))
    call dshr_fldList_add(fldsImport, 'Faxa_swvdr' )
    call dshr_fldList_add(fldsImport, 'Faxa_swvdf' )
    call dshr_fldList_add(fldsImport, 'Faxa_swndr' )
    call dshr_fldList_add(fldsImport, 'Faxa_swndf' )
    call dshr_fldList_add(fldsImport, 'Fioo_q'     )
    call dshr_fldList_add(fldsImport, 'Sa_z'       )
    call dshr_fldList_add(fldsImport, 'Sa_u'       )
    call dshr_fldList_add(fldsImport, 'Sa_v'       )
    call dshr_fldList_add(fldsImport, 'Sa_ptem'    )
    call dshr_fldList_add(fldsImport, 'Sa_shum'    )
    call dshr_fldList_add(fldsImport, 'Sa_dens'    )
    call dshr_fldList_add(fldsImport, 'Sa_tbot'    )
    call dshr_fldList_add(fldsImport, 'So_s'       )
    call dshr_fldList_add(fldsImport, 'Faxa_bcph'  , ungridded_lbound=1, ungridded_ubound=3)
    call dshr_fldList_add(fldsImport, 'Faxa_ocph'  , ungridded_lbound=1, ungridded_ubound=3)
    call dshr_fldList_add(fldsImport, 'Faxa_dstdry', ungridded_lbound=1, ungridded_ubound=4)
    call dshr_fldList_add(fldsImport, 'Faxa_dstwet', ungridded_lbound=1, ungridded_ubound=4)

    fldlist => fldsImport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(importState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(dice_comp_advertise): Fr_ice'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

  end subroutine dice_comp_advertise

!===============================================================================

  subroutine dice_comp_realize(sdat, importState, exportState, flds_scalar_name, flds_scalar_num, mesh, &
       logunit, masterproc, rc)

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
    character(*), parameter   :: subName = "(dice_comp_realize) "
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! -------------------------------------
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    ! -------------------------------------

    call dshr_fldlist_realize( exportState, fldsExport, flds_scalar_name, flds_scalar_num, mesh, &
         subname//':diceExport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_fldlist_realize( importState, fldsImport, flds_scalar_name, flds_scalar_num, mesh, &
         subname//':diceImport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -------------------------------------
    ! Set pointers to exportState fields
    ! -------------------------------------

    ! Initialize dfields with export state data that has corresponding stream field
    call dshr_dfield_add(dfields, sdat, state_fld='Si_ifrac', strm_fld='ifrac', &
         state=exportState, state_ptr=Si_ifrac, logunit=logunit, masterproc=masterproc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (flds_i2o_per_cat) then
       call state_getfldptr(exportState, fldname='Si_ifrac_n'        , fldptr2=Si_ifrac_n        , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(exportState, fldname='Fioi_swpen_ifrac_n', fldptr2=Fioi_swpen_ifrac_n, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    ! Set pointers to exportState fields that have no corresponding stream field
    call state_getfldptr(exportState, fldname='Si_imask'    , fldptr1=Si_imask    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_t'        , fldptr1=Si_t        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_tref'     , fldptr1=Si_tref     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_qref'     , fldptr1=Si_qref     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_avsdr'    , fldptr1=Si_avsdr    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_anidr'    , fldptr1=Si_anidr    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_avsdf'    , fldptr1=Si_avsdf    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Si_anidf'    , fldptr1=Si_anidf    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Faii_swnet'  , fldptr1=Faii_swnet  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Faii_sen'    , fldptr1=Faii_sen    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Faii_lat'    , fldptr1=Faii_lat    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Faii_lwup'   , fldptr1=Faii_lwup   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Faii_evap'   , fldptr1=Faii_evap   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Faii_taux'   , fldptr1=Faii_taux   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Faii_tauy'   , fldptr1=Faii_tauy   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_melth'  , fldptr1=Fioi_melth  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_meltw'  , fldptr1=Fioi_meltw  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_swpen'  , fldptr1=Fioi_swpen  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_taux'   , fldptr1=Fioi_taux   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_tauy'   , fldptr1=Fioi_tauy   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_salt'   , fldptr1=Fioi_salt   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_bcpho'  , fldptr1=Fioi_bcpho  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_bcphi'  , fldptr1=Fioi_bcphi  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(exportState, fldname='Fioi_flxdst' , fldptr1=Fioi_flxdst , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Set Si_imask
    call dshr_get_griddata(sdat, 'mask', Si_imask)

    ! -------------------------------------
    ! Set pointers to importState fields
    ! -------------------------------------

    call state_getfldptr(importState, fldname='Faxa_swvdr'  , fldptr1=Faxa_swvdr  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Faxa_swvdf'  , fldptr1=Faxa_swvdf  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Faxa_swndr'  , fldptr1=Faxa_swndr  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Faxa_swndf'  , fldptr1=Faxa_swndf  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Faxa_bcph'   , fldptr2=Faxa_bcph   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Faxa_ocph'   , fldptr2=Faxa_ocph   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Faxa_dstdry' , fldptr2=Faxa_dstdry , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Faxa_dstwet' , fldptr2=Faxa_dstwet , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Fioo_q'      , fldptr1=Fioo_q      , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Sa_z'        , fldptr1=Sa_z        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Sa_u'        , fldptr1=Sa_u        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Sa_v'        , fldptr1=Sa_v        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Sa_ptem'     , fldptr1=Sa_ptem     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Sa_tbot'     , fldptr1=Sa_tbot     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Sa_shum'     , fldptr1=Sa_shum     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='Sa_dens'     , fldptr1=Sa_dens     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call state_getfldptr(importState, fldname='So_s'        , fldptr1=So_s        , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! initialize import arrays
    ! On initial call, these are used for the first use in generating the export state
    ! and should have no impact on the solution!!

    Faxa_swvdr(:)    = 0._r8
    Faxa_swvdf(:)    = 0._r8
    Faxa_swndr(:)    = 0._r8
    Faxa_swndf(:)    = 0._r8
    Faxa_bcph(:,:)   = 0._r8
    Faxa_ocph(:,:)   = 0._r8
    Faxa_dstdry(:,:) = 0._r8
    Faxa_dstwet(:,:) = 0._r8
    Fioo_q(:)        = 0._r8
    Sa_z(:)          = 10.0_r8
    Sa_u(:)          = 5.0_r8
    Sa_v(:)          = 5.0_r8
    Sa_ptem(:)       = 260.0_r8
    Sa_tbot(:)       = 260.0_r8
    Sa_shum(:)       = 0.0014_r8
    Sa_dens(:)       = 1.3_r8
    So_s(:)          = 0._r8

    ! -------------------------------------
    ! allocate module arrays that are not part of import and export state
    ! -------------------------------------

    lsize = size(Si_ifrac)
    allocate(yc(lsize))
    allocate(water(lsize))
    allocate(tfreeze(lsize))
    allocate(imask(lsize))
    imask(:) = nint(Si_imask)

  end subroutine dice_comp_realize

!===============================================================================

  subroutine dice_comp_run(mpicom, my_task, master_task, logunit, target_ymd, target_tod, sdat, &
       mesh, cosarg, flux_swpf, flux_Qmin, flux_Qacc, flux_Qacc0, dt, read_restart, rc)

    ! --------------------------
    ! advance dice
    ! --------------------------

    ! input/output variables:
    integer                , intent(in)    :: mpicom     ! mpi communicator
    integer                , intent(in)    :: my_task
    integer                , intent(in)    :: master_task
    integer                , intent(in)    :: logunit
    integer                , intent(in)    :: target_ymd ! model date
    integer                , intent(in)    :: target_tod ! model sec into model date
    type(shr_strdata_type) , intent(inout) :: sdat
    type(ESMF_Mesh)        , intent(in)    :: mesh
    real(R8)               , intent(in)    :: cosarg     ! for setting ice temp pattern
    real(r8)               , intent(in)    :: flux_swpf  ! short-wave penatration factor
    real(r8)               , intent(in)    :: flux_Qmin  ! bound on melt rate
    logical                , intent(in)    :: flux_Qacc  ! activates water accumulation/melt wrt Q
    real(r8)               , intent(in)    :: flux_Qacc0 ! initial water accumulation value
    logical                , intent(in)    :: read_restart
    real(r8)               , intent(in)    :: dt
    integer                , intent(out)   :: rc

    ! local variables
    integer           :: n, lsize
    real(r8)          :: jday, jday0        ! elapsed day counters
    integer           :: spatialDim         ! number of dimension in mesh
    integer           :: numOwnedElements   ! size of mesh
    real(r8), pointer :: ownedElemCoords(:) ! mesh lat and lons
    real(r8)          :: qmeltall          ! q that would melt all accumulated water
    logical           :: first_time = .true.
    character(*), parameter :: subName = "(dice_comp_run) "
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('DICE_RUN')

    !--------------------
    ! advance dice streams
    !--------------------

    ! time and spatially interpolate to model time and grid
    call t_barrierf('dice_BARRIER',mpicom)
    call t_startf('dice_strdata_advance')
    call shr_strdata_advance(sdat, target_ymd, target_tod, mpicom, 'dice')
    call t_stopf('dice_strdata_advance')

    !--------------------
    ! copy all fields from streams to export state as default
    !--------------------

    ! This automatically will update the fields in the export state
    call t_barrierf('dice_comp_dfield_copy_BARRIER', mpicom)
    call t_startf('dice_dfield_copy')
    call dshr_dfield_copy(dfields, sdat, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('dice_dfield_copy')

    !-------------------------------------------------
    ! Determine data model behavior based on the mode
    !-------------------------------------------------

    lsize = size(Si_ifrac)

    call t_startf('dice_datamode')
    select case (trim(sdat%datamode))

    case('COPYALL')
       ! do nothing extra

    case('SSTDATA')
       if (first_time) then
          if (.not. read_restart) then
             do n = 1,lsize
                if (Si_ifrac(n) > 0.0_r8) then
                   water(n) = flux_Qacc0
                else
                   water(n) = 0.0_r8
                end if
             end do
             ! iFrac0 = iFrac  ! previous step's ice fraction
          endif

          call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(ownedElemCoords(spatialDim*numOwnedElements))
          call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do n = 1,numOwnedElements
             yc(n) = ownedElemCoords(2*n)
          end do
          deallocate(ownedElemCoords)
       end if

       tfreeze(:) = shr_frz_freezetemp(So_s(:)) + tFrz ! convert to Kelvin

       do n = 1,lsize

          !--- fix erroneous iFrac ---
          Si_ifrac(n) = min(1.0_r8,max(0.0_r8,Si_ifrac(n)))

          !--- fabricate ice surface T, fix erroneous iFrac ---
          if ( yc(n) > 0.0_r8) then
             Si_t(n) = 260.0_r8 + 10.0_r8*cos(cosArg)
          else
             Si_t(n) = 260.0_r8 - 10.0_r8*cos(cosArg)
          end if

          !--- set albedos (constant) ---
          Si_avsdr(n) = ax_vsdr
          Si_anidr(n) = ax_nidr
          Si_avsdf(n) = ax_vsdf
          Si_anidf(n) = ax_nidf

          !--- swnet is sent to cpl as a diagnostic quantity only ---
          !--- newly recv'd swdn goes with previously sent albedo ---
          !--- but albedos are (currently) time invariant         ---
          Faii_swnet(n)   = (1.0_r8 - Si_avsdr(n))*Faxa_swvdr(n) &
                          + (1.0_r8 - Si_anidr(n))*Faxa_swndr(n) &
                          + (1.0_r8 - Si_avsdf(n))*Faxa_swvdf(n) &
                          + (1.0_r8 - Si_anidf(n))*Faxa_swndf(n)

          !--- compute melt/freeze water balance, adjust iFrac  -------------
          if ( .not. flux_Qacc ) then ! Q accumulation option is OFF
             Fioi_melth(n) = min(Fioo_q(n),0.0_r8 )          ! q<0 => melt potential
             Fioi_melth(n) = max(Fioi_melth(n),Flux_Qmin   ) ! limit the melt rate
             Fioi_meltw(n) =    -Fioi_melth(n)/latice        ! corresponding water flux

          else                                 ! Q accumulation option is ON
             !--------------------------------------------------------------
             ! 1a) Q<0 & iFrac > 0  =>  infinite supply of water to melt
             ! 1b) Q<0 & iFrac = 0  =>  melt accumulated water only
             ! 2a) Q>0 & iFrac > 0  =>  zero-out accumulated water
             ! 2b) Q>0 & iFrac = 0  =>  accumulated water
             !--------------------------------------------------------------

             if ( Fioo_q(n) <  0.0_r8 ) then ! Q<0 => melt
                if (Si_ifrac(n) > 0.0_r8 ) then
                   Fioi_melth(n) = Si_ifrac(n)*max(Fioo_q(n),Flux_Qmin)
                   Fioi_meltw(n) =    -Fioi_melth(n)/latice
                   !  water(n) = < don't change this value >
                else
                   Qmeltall   = -water(n)*latice/dt
                   Fioi_melth(n) = max(Fioo_q(n), Qmeltall, Flux_Qmin )
                   Fioi_meltw(n) = -Fioi_melth(n)/latice
                   water(n) =  water(n) - Fioi_meltw(n)*dt
                end if
             else                       ! Q>0 => freeze
                if (Si_ifrac(n) > 0.0_r8 ) then
                   Fioi_melth(n) = 0.0_r8
                   Fioi_meltw(n) = 0.0_r8
                   water(n) = 0.0_r8
                else
                   Fioi_melth(n) = 0.0_r8
                   Fioi_meltw(n) = 0.0_r8
                   water(n) = water(n) + dt*Fioo_q(n)/latice
                end if
             end if

             if (water(n) < 1.0e-16_r8 ) water(n) = 0.0_r8

             !--- non-zero water => non-zero iFrac ---
             if (Si_ifrac(n) <= 0.0_r8  .and.  water(n) > 0.0_r8) then
                Si_ifrac(n) = min(1.0_r8,water(n)/waterMax)
                ! Si_t(n) = tfreeze(n)     ! T can be above freezing?!?
             end if

             !--- cpl multiplies Fioi_melth & Fioi_meltw by iFrac ---
             !--- divide by iFrac here => fixed quantity flux (not per area) ---
             if (Si_ifrac(n) > 0.0_r8) then
                Si_ifrac(n) = max( 0.01_r8, Si_ifrac(n)) ! min iFrac
                Fioi_melth(n) = Fioi_melth(n)/Si_ifrac(n)
                Fioi_meltw(n) = Fioi_meltw(n)/Si_ifrac(n)
             else
                Fioi_melth(n) = 0.0_r8
                Fioi_meltw(n) = 0.0_r8
             end if
          end if

          !--- modify T wrt iFrac: (iFrac -> 0) => (T -> tfreeze) ---
          Si_t(n) = tfreeze(n) + Si_ifrac(n)*(Si_t(n)-tfreeze(n))
       end do

       ! compute ice/ice surface fluxes
       call dice_flux_atmice( &
            imask     ,Sa_z      ,Sa_u      ,Sa_v     , &
            Sa_ptem   ,Sa_shum   ,Sa_dens   ,Sa_tbot  , &
            Si_t      ,Faii_sen  ,Faii_lat  ,Faii_lwup, &
            Faii_evap ,Faii_taux ,Faii_tauy ,Si_tref  , &
            Si_qref   ,logunit )

       ! compute ice/oce surface fluxes (except melth & meltw, see above)
       do n=1,lsize
          if (imask(n) == 0) then
             Fioi_swpen(n) = spval
             Fioi_melth(n) = spval
             Fioi_meltw(n) = spval
             Fioi_salt (n) = spval
             Fioi_taux(n)  = spval
             Fioi_tauy(n)  = spval
             Si_ifrac(n)   = 0.0_r8
          else
             !--- penetrating short wave ---
             Fioi_swpen(n) = max(0.0_r8, flux_swpf*Faii_swnet(n) ) ! must be non-negative

             !--- i/o surface stress ( = atm/ice stress) ---
             Fioi_taux(n) = Faii_taux(n)
             Fioi_tauy(n) = Faii_tauy(n)

             !--- salt flux ---
             Fioi_salt(n) = 0.0_r8
          end if
          ! !--- save ifrac for next timestep
          ! iFrac0(n) = Si_ifrac(n)
       end do

       ! Compute outgoing aerosol fluxes
       do n = 1,lsize
          Fioi_bcpho(n) = Faxa_bcph(2,n)
          Fioi_bcphi(n) = Faxa_bcph(1,n) + Faxa_bcph(3,n)
          Fioi_flxdst(n) =  Faxa_dstdry(1,n) + Faxa_dstwet(1,n) + Faxa_dstdry(2,n) + Faxa_dstwet(2,n) &
                          + Faxa_dstdry(3,n) + Faxa_dstwet(3,n) + Faxa_dstdry(4,n) + Faxa_dstwet(4,n)
       end do

    end select

    !-------------------------------------------------
    ! optional per thickness category fields
    !-------------------------------------------------

    if (flds_i2o_per_cat) then
       do n = 1,lsize
          Si_iFrac_n(1,n)         = Si_ifrac(n)
          Fioi_swpen_iFrac_n(1,n) = Fioi_swpen(n) * Si_ifrac(n)
       end do
    end if

    first_time = .false.

    call t_stopf('dice_datamode')
    call t_stopf('DICE_RUN')

  end subroutine dice_comp_run

end module dice_comp_mod
