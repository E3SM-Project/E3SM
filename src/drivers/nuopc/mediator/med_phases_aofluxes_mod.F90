module med_phases_aofluxes_mod

  use ESMF
  use NUOPC
  use shr_kind_mod          , only : R8=>SHR_KIND_R8, IN=>SHR_KIND_IN
  use shr_kind_mod          , only : CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
  use shr_sys_mod           , only : shr_sys_abort
  use esmFlds               , only : flds_scalar_name
  use esmFlds               , only : flds_scalar_num
  use esmFlds               , only : fldListFr, fldListTo
  use esmFlds               , only : compatm, compocn, compname
  use esmFlds               , only : fldListMed_aoflux_o
  use shr_nuopc_fldList_mod , only : shr_nuopc_fldlist_getfldnames
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_init
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFieldN
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_GetScalar
  use shr_flux_mod          , only : shr_flux_atmocn, shr_flux_atmocn_diurnal, shr_flux_adjust_constants
  use shr_const_mod         , only : shr_const_spval, shr_const_pi
  use seq_timemgr_mod       , only : seq_timemgr_EclockGetData
  use med_constants_mod     , only : med_constants_dbug_flag
  use med_constants_mod     , only : med_constants_czero
  use med_map_mod           , only : med_map_FB_Regrid_Norm
  use med_internalstate_mod , only : InternalState

  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public routines
  !--------------------------------------------------------------------------

  public  :: med_phases_aofluxes_init
  public  :: med_phases_aofluxes_run

  !--------------------------------------------------------------------------
  ! Private routines
  !--------------------------------------------------------------------------

  private :: med_aofluxes_init
  private :: med_aofluxes_run

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  real(ESMF_KIND_R8) , pointer :: lats        (:) ! latitudes  (degrees)
  real(ESMF_KIND_R8) , pointer :: lons        (:) ! longitudes (degrees)
  integer            , pointer :: mask        (:) ! ocn domain mask: 0 <=> inactive cell
  real(r8)           , pointer :: rmask       (:) ! ocn domain mask: 0 <=> inactive cell
  real(ESMF_KIND_R8) , pointer :: uocn        (:) ! ocn velocity, zonal
  real(ESMF_KIND_R8) , pointer :: vocn        (:) ! ocn velocity, meridional
  real(ESMF_KIND_R8) , pointer :: tocn        (:) ! ocean temperature
  real(ESMF_KIND_R8) , pointer :: zbot        (:) ! atm level height
  real(ESMF_KIND_R8) , pointer :: ubot        (:) ! atm velocity, zonal
  real(ESMF_KIND_R8) , pointer :: vbot        (:) ! atm velocity, meridional
  real(ESMF_KIND_R8) , pointer :: thbot       (:) ! atm potential T
  real(ESMF_KIND_R8) , pointer :: shum        (:) ! atm specific humidity
  real(ESMF_KIND_R8) , pointer :: shum_16O    (:) ! atm H2O tracer
  real(ESMF_KIND_R8) , pointer :: shum_HDO    (:) ! atm HDO tracer
  real(ESMF_KIND_R8) , pointer :: shum_18O    (:) ! atm H218O tracer
  real(ESMF_KIND_R8) , pointer :: roce_16O    (:) ! ocn H2O ratio
  real(ESMF_KIND_R8) , pointer :: roce_HDO    (:) ! ocn HDO ratio
  real(ESMF_KIND_R8) , pointer :: roce_18O    (:) ! ocn H218O ratio
  real(ESMF_KIND_R8) , pointer :: dens        (:) ! atm density
  real(ESMF_KIND_R8) , pointer :: tbot        (:) ! atm bottom surface T
  real(ESMF_KIND_R8) , pointer :: sen         (:) ! heat flux: sensible
  real(ESMF_KIND_R8) , pointer :: lat         (:) ! heat flux: latent
  real(ESMF_KIND_R8) , pointer :: lwup        (:) ! lwup over ocean
  real(ESMF_KIND_R8) , pointer :: evap        (:) ! water flux: evaporation
  real(ESMF_KIND_R8) , pointer :: evap_16O    (:) ! H2O flux: evaporation
  real(ESMF_KIND_R8) , pointer :: evap_HDO    (:) ! HDO flux: evaporation
  real(ESMF_KIND_R8) , pointer :: evap_18O    (:) ! H218O flux: evaporation
  real(ESMF_KIND_R8) , pointer :: taux        (:) ! wind stress, zonal
  real(ESMF_KIND_R8) , pointer :: tauy        (:) ! wind stress, meridional
  real(ESMF_KIND_R8) , pointer :: tref        (:) ! diagnostic:  2m ref T
  real(ESMF_KIND_R8) , pointer :: qref        (:) ! diagnostic:  2m ref Q
  real(ESMF_KIND_R8) , pointer :: u10         (:) ! diagnostic: 10m wind speed
  real(ESMF_KIND_R8) , pointer :: duu10n      (:) ! diagnostic: 10m wind speed squared
  real(ESMF_KIND_R8) , pointer :: fswpen      (:) ! fraction of sw penetrating ocn surface layer
  real(ESMF_KIND_R8) , pointer :: ocnsal      (:) ! ocean salinity
  real(ESMF_KIND_R8) , pointer :: lwdn        (:) ! long  wave, downward
  real(ESMF_KIND_R8) , pointer :: swdn        (:) ! short wave, downward
  real(ESMF_KIND_R8) , pointer :: swup        (:) ! short wave, upward
  real(ESMF_KIND_R8) , pointer :: rainc       (:) ! rainc
  real(ESMF_KIND_R8) , pointer :: rainl       (:) ! rainl
  real(ESMF_KIND_R8) , pointer :: snowc       (:) ! snowc
  real(ESMF_KIND_R8) , pointer :: snowl       (:) ! snowl
  real(ESMF_KIND_R8) , pointer :: tbulk       (:) ! diurnal diagnostic: ocn bulk T
  real(ESMF_KIND_R8) , pointer :: tskin       (:) ! diurnal diagnostic: ocn skin T
  real(ESMF_KIND_R8) , pointer :: tskin_night (:) ! diurnal diagnostic: ocn skin T
  real(ESMF_KIND_R8) , pointer :: tskin_day   (:) ! diurnal diagnostic: ocn skin T
  real(ESMF_KIND_R8) , pointer :: cSkin       (:) ! diurnal diagnostic: ocn cool skin
  real(ESMF_KIND_R8) , pointer :: cSkin_night (:) ! diurnal diagnostic: ocn cool skin
  real(ESMF_KIND_R8) , pointer :: warm        (:) ! diurnal diagnostic: ocn warming
  real(ESMF_KIND_R8) , pointer :: warmMax     (:) ! diurnal diagnostic: ocn warming, max daily value
  real(ESMF_KIND_R8) , pointer :: warmMaxInc  (:) ! diurnal diagnostic: ocn warming, max daily value, increment
  real(ESMF_KIND_R8) , pointer :: salt        (:) ! diurnal diagnostic: ocn salting
  real(ESMF_KIND_R8) , pointer :: speed       (:) ! diurnal diagnostic: ocn speed
  real(ESMF_KIND_R8) , pointer :: regime      (:) ! diurnal diagnostic: ocn regime
  real(ESMF_KIND_R8) , pointer :: windMax     (:) ! diurnal diagnostic: ocn wind   , max daily value
  real(ESMF_KIND_R8) , pointer :: windAvg     (:) ! diurnal diagnostic: ocn wind   , daily avg
  real(ESMF_KIND_R8) , pointer :: windMaxInc  (:) ! diurnal diagnostic: ocn wind   , max daily value, increment
  real(ESMF_KIND_R8) , pointer :: QsolAvg     (:) ! diurnal diagnostic: ocn Qsol   , daily avg
  real(ESMF_KIND_R8) , pointer :: qSolInc     (:) ! diurnal diagnostic: ocn Qsol   , daily avg, increment
  real(ESMF_KIND_R8) , pointer :: windInc     (:) ! diurnal diagnostic: ocn wind   , daily avg, increment
  real(ESMF_KIND_R8) , pointer :: nInc        (:) ! diurnal diagnostic: a/o flux   , increment
  real(ESMF_KIND_R8) , pointer :: ustar       (:) ! saved ustar
  real(ESMF_KIND_R8) , pointer :: re          (:) ! saved re
  real(ESMF_KIND_R8) , pointer :: ssq         (:) ! saved sq

  ! Fields that are not obtained via GetFldPtr
  real(ESMF_KIND_R8) , pointer :: uGust       (:) ! wind gust
  real(ESMF_KIND_R8) , pointer :: prec        (:) ! precip
  real(ESMF_KIND_R8) , pointer :: prec_gust   (:) ! atm precip for convective gustiness (kg/m^3)

  ! The following three variables are obtained as attributes from gcomp
  logical       :: flds_wiso  ! use case
  logical       :: do_flux_diurnal
  character(CL) :: aoflux_grid

  ! Conversion from degrees to radians
  integer                        :: dbug_flag      = med_constants_dbug_flag
  real(ESMF_KIND_R8) , parameter :: czero          = med_constants_czero
  real(ESMF_KIND_R8) , parameter :: const_deg2rad  = shr_const_pi/180.0_r8  ! deg to rads
  character(*)       , parameter :: u_FILE_u       = __FILE__
  integer                        :: dbrc
  character(len=1024)            :: tmpstr
  logical                        :: mastertask

  ! Mediator field bundles used for atm/ocn flux computation
  type(ESMF_FieldBundle), public  :: FBMed_aoflux_accum_o  ! Ocn/Atm flux accumulator on ocn grid
  type(ESMF_FieldBundle), public  :: FBMed_aoflux_diurnl_o ! Ocn/Atm flux fields only needed for history
  type(ESMF_FieldBundle), public  :: FBMed_aoflux_diurnl_a ! Ocn/Atm flux fields only needed for history

!================================================================================
contains
!================================================================================

  subroutine med_phases_aofluxes_init(gcomp, rc)

    ! Initialize ocn/atm flux calculations

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    character(CL)          :: cvalue
    character(CL)          :: aoflux_grid
    type(InternalState)    :: is_local
    integer                :: nflds
    integer                :: localPet
    type(ESMF_VM)          :: vm
    character(CL), pointer :: fldnames(:)
    character(len=*),parameter :: subname='(med_phases_aofluxes_init)'
    !---------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    mastertask = .false.
    if (localPet == 0) mastertask=.true.

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Create module field bundles
    
    call NUOPC_CompAttributeGet(gcomp, name='flux_diurnal', value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    read(cvalue,*) do_flux_diurnal

    nflds = size(fldListMed_aoflux_o%flds)
    allocate(fldnames(nflds))
    call shr_nuopc_fldList_getfldnames(fldListMed_aoflux_o%flds, fldnames)

    call shr_nuopc_methods_FB_init(FBMed_aoflux_accum_o, flds_scalar_name, &
         STgeom=is_local%wrap%NStateImp(compocn), fieldnamelist=fldnames, name='FBMed_aoflux_o_accum', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (do_flux_diurnal) then
       call shr_nuopc_methods_FB_init(FBMed_aoflux_diurnl_o, flds_scalar_name, &
            STgeom=is_local%wrap%NStateImp(compocn), fieldnamelist=fldnames, name='FBMed_diurnl_o_accum', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_init(FBMed_aoflux_diurnl_a, flds_scalar_name, &
            STgeom=is_local%wrap%NStateImp(compatm), fieldnamelist=fldnames, name='FBMed_diurnl_a_accum', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    deallocate(fldnames)

    ! Determine src and dst comps depending on the aoflux_grid setting

    call NUOPC_CompAttributeGet(gcomp, name='aoflux_grid', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) aoflux_grid

    if (trim(aoflux_grid) == 'ocn') then

       ! Create FBMed_aoflux_o (field bundle on the ocean grid)
       call med_aofluxes_init(gcomp, &
            FBAtm=is_local%wrap%FBImp(compatm,compocn), &
            FBOcn=is_local%wrap%FBImp(compocn,compocn), &
            FBFrac=is_local%wrap%FBfrac(compocn), &
            FBMed_ocnalb=is_local%wrap%FBMed_ocnalb_o, &
            FBMed_aoflux=is_local%wrap%FBMed_aoflux_o, &
            FBMed_aoflux_diurnl=FBMed_aoflux_diurnl_o, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    else if (trim(aoflux_grid) == 'atm') then

       ! Create FBMed_aoflux_a (field bundle on the atmosphere grid)
       call med_aofluxes_init(gcomp, &
            FBAtm=is_local%wrap%FBImp(compatm,compatm), &
            FBOcn=is_local%wrap%FBImp(compocn,compatm), &
            FBFrac=is_local%wrap%FBfrac(compatm), &
            FBMed_ocnalb=is_local%wrap%FBMed_ocnalb_a, &
            FBMed_aoflux=is_local%wrap%FBMed_aoflux_a, &
            FBMed_aoflux_diurnl=FBMed_aoflux_diurnl_a, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    else

       call ESMF_LogWrite(trim(subname)//' aoflux_grid = '//trim(aoflux_grid)//' not available', &
            ESMF_LOGMSG_INFO, rc=dbrc)
       return

    end if

  end subroutine med_phases_aofluxes_init

!================================================================================

  subroutine med_phases_aofluxes_run(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    type(ESMF_Clock)    :: clock
    character(CL)       :: cvalue
    logical             :: ocn_prognostic
    character(CL)       :: aoflux_grid
    character(len=*),parameter :: subname='(med_phases_aofluxes)'
    !---------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! Get the clock from the mediator Component
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from the mediator Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine if ocean is prognostic
    ocn_prognostic = NUOPC_IsConnected(is_local%wrap%NStateImp(compocn))

    ! Determine source and destination comps depending on the aoflux_grid setting
    call NUOPC_CompAttributeGet(gcomp, name='aoflux_grid', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) aoflux_grid

    if (trim(aoflux_grid) == 'ocn') then

       ! Regrid atm import field bundle from atm to ocn grid as input for ocn/atm flux calculation
       call med_map_FB_Regrid_Norm( &
            fldListFr(compatm)%flds, compatm, compocn, &
            is_local%wrap%FBImp(compatm,compatm), &
            is_local%wrap%FBImp(compatm,compocn), &
            is_local%wrap%FBFrac(compatm), &
            is_local%wrap%FBNormOne(compatm,compocn,:), &
            is_local%wrap%RH(compatm,compocn,:), &
            string=trim(compname(compatm))//'2'//trim(compname(compocn)), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Calculate atm/ocn fluxes on the destination grid
       call med_aofluxes_run(gcomp, clock, ocn_prognostic, rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBMed_aoflux_o, string=trim(subname) //' FBAMed_aoflux_o' , rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

    else if (trim(aoflux_grid) == 'atm') then

       call med_map_FB_Regrid_Norm( &
            fldListFr(compocn)%flds, compocn, compatm, &
            is_local%wrap%FBImp(compocn,compocn), &
            is_local%wrap%FBImp(compocn,compatm), &
            is_local%wrap%FBFrac(compocn), &
            is_local%wrap%FBNormOne(compocn,compatm,:), &
            is_local%wrap%RH(compocn,compatm,:), &
            string=trim(compname(compocn))//'2'//trim(compname(compatm)), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImp(compocn,compatm), &
               string=trim(subname) //' FBImp('//trim(compname(compocn))//','//trim(compname(compatm))//') ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! Calculate atm/ocn fluxes on the destination grid
       call med_aofluxes_run(gcomp, clock, ocn_prognostic, rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImp(compocn,compatm), &
               string=trim(subname) //' FBImp('//trim(compname(compocn))//','//trim(compname(compatm))//') ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

    else

       call ESMF_LogWrite(trim(subname)//' aoflux_grid = '//trim(aoflux_grid)//' not available', &
            ESMF_LOGMSG_INFO, rc=dbrc)
       return

    end if

  end subroutine med_phases_aofluxes_run

!================================================================================

  subroutine med_aofluxes_init(gcomp, FBAtm, FBOcn, FBFrac, FBMed_ocnalb, &
       FBMed_aoflux, FBMed_aoflux_diurnl, rc)

    !-----------------------------------------------------------------------
    ! Initialize pointers to the module variables 
    !-----------------------------------------------------------------------

    ! Arguments
    type(ESMF_GridComp)                    :: gcomp
    type(ESMF_FieldBundle) , intent(in)    :: FBAtm               ! Atm Import fields on aoflux grid
    type(ESMF_FieldBundle) , intent(in)    :: FBOcn               ! Ocn Import fields on aoflux grid
    type(ESMF_FieldBundle) , intent(in)    :: FBfrac              ! Fraction data for various components, on their grid
    type(ESMF_FieldBundle) , intent(in)    :: FBMed_ocnalb        ! Ocn albedos computed in mediator 
    type(ESMF_FieldBundle) , intent(inout) :: FBMed_aoflux        ! Ocn albedos computed in mediator 
    type(ESMF_FieldBundle) , intent(inout) :: FBMed_aoflux_diurnl ! Ocn albedos computed in mediator 
    integer                , intent(out)   :: rc
    !
    ! Local variables
    type(ESMF_VM)            :: vm
    integer                  :: iam
    type(ESMF_Field)         :: lfield
    type(ESMF_Grid)          :: lgrid
    type(ESMF_Mesh)          :: lmesh
    type(ESMF_GeomType_Flag) :: geomtype
    integer                  :: n
    integer                  :: lsize
    real(r8), pointer        :: ofrac(:)
    real(r8), pointer        :: ifrac(:)
    integer                  :: dimCount
    integer                  :: spatialDim
    integer                  :: numOwnedElements
    character(CL)            :: cvalue
    real(ESMF_KIND_R8), pointer :: ownedElemCoords(:)
    character(*),parameter   :: subName =   '(med_aofluxes_init) '
    !-----------------------------------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! The following is for debugging
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! get attributes that are set as module variables
    !----------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='flds_wiso', value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    read(cvalue,*) flds_wiso

    call NUOPC_CompAttributeGet(gcomp, name='aoflux_grid', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) aoflux_grid

    !----------------------------------
    ! fields calculated in flux_ocnalb
    !----------------------------------

    call shr_nuopc_methods_FB_GetFldPtr(FBMed_ocnalb, fldname='Faox_swdn', fldptr1=swdn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_ocnalb, fldname='Faox_swup', fldptr1=swup, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! atm/ocn fields
    !----------------------------------

    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_tref', fldptr1=tref, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_qref', fldptr1=qref, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_ustar', fldptr1=ustar, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_re', fldptr1=re, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_ssq', fldptr1=ssq, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_u10', fldptr1=u10, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_duu10n', fldptr1=duu10n, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_taux', fldptr1=taux, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_tauy', fldptr1=tauy, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_lat', fldptr1=lat, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_sen', fldptr1=sen, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_evap', fldptr1=evap, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    lsize = size(evap)
    if (flds_wiso) then
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_evap_16O', fldptr1=evap_16O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_evap_18O', fldptr1=evap_18O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_evap_HDO', fldptr1=evap_HDO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(evap_16O(lsize)); evap_16O(:) = 0._r8
       allocate(evap_18O(lsize)); evap_18O(:) = 0._r8
       allocate(evap_HDO(lsize)); evap_HDO(:) = 0._r8
    end if

    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_lwup', fldptr1=lwup, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_fswpen', fldptr1=fswpen, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (do_flux_diurnal) then
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_warm_diurn', fldptr1=warm, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_salt_diurn', fldptr1=salt, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_speed_diurn', fldptr1=speed, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_regime_diurn', fldptr1=regime, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_warmmax_diurn', fldptr1=warmMax, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_windmax_diurn', fldptr1=windMax, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_warmmaxinc_diurn', fldptr1=warmMaxInc, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_qsolinc_diurn', fldptr1=qSolInc, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_windinc_diurn', fldptr1=windInc, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_ninc_diurn', fldptr1=nInc, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_tbulk_diurn', fldptr1=tbulk, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_tskin_diurn', fldptr1=tskin, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_tskin_day_diurn', fldptr1=tskin_day, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_tskin_night_diurn', fldptr1=tskin_night, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_cskin_diurn', fldptr1=cskin, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_cskin_night_diurn', fldptr1=cskin_night, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_windavg_diurn', fldptr1=windAvg, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_qsolavg_diurn', fldptr1=qSolAvg, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux_diurnl, fldname='So_windmaxinc_diurn', fldptr1=windmaxinc, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    !----------------------------------
    ! Ocn import fields
    !----------------------------------

    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_omask', fldptr1=rmask, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_t', fldptr1=tocn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_u', fldptr1=uocn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_v', fldptr1=vocn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_s', fldptr1=ocnsal, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_fswpen', fldptr1=fswpen, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (flds_wiso) then
       call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_roce_16O', fldptr1=roce_16O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_roce_18O', fldptr1=roce_18O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_roce_HDO', fldptr1=roce_HDO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(roce_16O(lsize)); roce_16O(:) = 0._r8
       allocate(roce_18O(lsize)); roce_18O(:) = 0._r8
       allocate(roce_HDO(lsize)); roce_HDO(:) = 0._r8
    end if

    !----------------------------------
    ! Atm import fields
    !----------------------------------

    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_z', fldptr1=zbot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_u', fldptr1=ubot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_v', fldptr1=vbot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_tbot', fldptr1=tbot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_ptem', fldptr1=thbot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_shum', fldptr1=shum, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (flds_wiso) then
       call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_shum_16O', fldptr1=shum_16O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_shum_18O', fldptr1=shum_18O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_shum_HDO', fldptr1=shum_HDO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(shum_16O(lsize)); shum_16O(:) = 0._r8
       allocate(shum_18O(lsize)); shum_18O(:) = 0._r8
       allocate(shum_HDO(lsize)); shum_HDO(:) = 0._r8
    end if

    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_dens', fldptr1=dens, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_lwdn', fldptr1=lwdn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_rainc', fldptr1=rainc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_rainl', fldptr1=rainl, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_snowc', fldptr1=snowc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_snowl', fldptr1=snowl, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! !----------------------------------
    ! ! Get lat, lon, which are time-invariant
    ! !----------------------------------

    if (do_flux_diurnal) then
       ! Get the first field from the field bundle - assumes that all fields
       ! in FBMed_aoflux have the same grid - so only need to query field 1
       call shr_nuopc_methods_FB_getFieldN(FBMed_aoflux, fieldnum=1, field=lfield, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Determine if first field is on a grid or a mesh - default will be mesh
       call ESMF_FieldGet(lfield, geomtype=geomtype, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (geomtype == ESMF_GEOMTYPE_MESH) then
          if (dbug_flag > 5) then
             call ESMF_LogWrite(trim(subname)//" : FBAtm is on a mesh ", ESMF_LOGMSG_INFO, rc=rc)
          end if
          call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call ESMF_MeshGet(lmesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          if (numOwnedElements /= lsize) then
             call ESMF_LogWrite(trim(subname)//": ERROR numOwnedElements not equal to local size", ESMF_LOGMSG_INFO, rc=rc)
             rc = ESMF_FAILURE
             return
          end if
          allocate(ownedElemCoords(spatialDim*numOwnedElements))
          allocate(lons(numOwnedElements))
          allocate(lats(numOwnedElements))
          call ESMF_MeshGet(lmesh, ownedElemCoords=ownedElemCoords)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
             !TODO: the following is only needed for MOM6 until ESMF is updated - uncomment if you are using MOM6
             lons(:) = 0.0
             lats(:) = 0.0
          else
             do n = 1,lsize
                lons(n) = ownedElemCoords(2*n-1)
                lats(n) = ownedElemCoords(2*n)
             end do
          end if
       else
          call ESMF_LogWrite(trim(subname)//": ERROR FBATM must be either on a mesh", ESMF_LOGMSG_INFO, rc=rc)
          rc = ESMF_FAILURE
          return
       end if
    end if

    !----------------------------------
    ! Fields that are not obtained via GetFldPtr
    !----------------------------------
    allocate(uGust(lsize))     ; uGust(:) =  0.0_r8
    allocate(prec(lsize))      ; prec(:) =  0.0_r8
    allocate(prec_gust(lsize)) ; prec_gust(:) =  0.0_r8

    !----------------------------------
    ! setup the compute mask.
    !----------------------------------

    ! allocate grid mask fields
    ! default compute everywhere, then "turn off" gridcells
    allocate(mask(lsize))
    mask(:) = 1

    write(tmpstr,'(i12,g22.12,i12)') lsize,sum(rmask),sum(mask)
    call ESMF_LogWrite(trim(subname)//" : maskA= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    where (rmask(:) == 0._r8) mask(:) = 0   ! like nint

    write(tmpstr,'(i12,g22.12,i12)') lsize,sum(rmask),sum(mask)
    call ESMF_LogWrite(trim(subname)//" : maskB= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    ! TODO: need to check if this logic is correct
    ! then check ofrac + ifrac
    ! call shr_nuopc_methods_FB_getFldPtr(FBFrac , fldname='ofrac' , fldptr1=ofrac, rc=rc)
    ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call shr_nuopc_methods_FB_getFldPtr(FBFrac , fldname='ifrac' , fldptr1=ifrac, rc=rc)
    ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! where (ofrac(:) + ifrac(:) <= 0.0_r8) mask(:) = 0

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_aofluxes_init

!===============================================================================

  subroutine med_aofluxes_run(gcomp, clock, ocn_prognostic, rc)

    !-----------------------------------------------------------------------
    ! Determine atm/ocn fluxes eother on atm or on ocean grid
    ! The module arrays are set via pointers the the mediator internal states
    ! in med_ocnatm_init and are used below.
    ! gcomp (the mediator gridded component) is only needed to retreive the
    ! attributes
    !-----------------------------------------------------------------------

    ! Arguments
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_Clock)      :: clock
    logical , intent(in)  :: ocn_prognostic
    integer , intent(out) :: rc
    !
    ! Local variables
    character(CL) :: cvalue
    integer       :: tod, dt
    integer       :: n,i                     ! indices
    integer       :: lsize                   ! local size
    real(r8)      :: gust_fac = huge(1.0_r8) ! wind gust factor
    logical       :: cold_start              ! .true. to initialize internal fields in shr_flux diurnal
    logical       :: read_restart            ! .true. => continue run
    real(r8)      :: flux_convergence        ! convergence criteria for imlicit flux computation
    integer       :: flux_max_iteration      ! maximum number of iterations for convergence
    logical       :: coldair_outbreak_mod    ! cold air outbreak adjustment  (Mahrt & Sun 1995,MWR)
    logical,save  :: first_call = .true.
    character(*),parameter :: F01 = "('(med_aofluxes_run) ',a,i4,2x,d21.14)"
    character(*),parameter :: F02 = "('(med_aofluxes_run) ',a,i4,2x,i4)"
    character(*),parameter :: subName =   '(med_fluxes_run) '
    !-----------------------------------------------------------------------

    call seq_timemgr_EClockGetData( clock, curr_tod=tod, dtime=dt )

    write(tmpstr,'(2i12)') tod,dt
    call ESMF_LogWrite(trim(subname)//" : tod,dt= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    cold_start = .false.   ! use restart data or data from last timestep
    if (first_call) then
       call NUOPC_CompAttributeGet(gcomp, name='gust_fac', value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) gust_fac

       call NUOPC_CompAttributeGet(gcomp, name='coldair_outbreak_mod', value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) coldair_outbreak_mod

       call NUOPC_CompAttributeGet(gcomp, name='flux_max_iteration', value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flux_max_iteration

       call NUOPC_CompAttributeGet(gcomp, name='flux_convergence', value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flux_convergence

       call shr_flux_adjust_constants(&
            flux_convergence_tolerance=flux_convergence, &
            flux_convergence_max_iteration=flux_max_iteration, &
            coldair_outbreak_mod=coldair_outbreak_mod)

       call NUOPC_CompAttributeGet(gcomp, name='read_restart', value=cvalue, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) read_restart
       if (.not.read_restart) then
          cold_start = .true.
       end if

       first_call = .false.
    end if

    !----------------------------------
    ! Determine the compute mask
    !----------------------------------

    ! Prefer to compute just where ocean exists, so setup a mask here.
    ! this could be run with either the ocean or atm grid so need to be careful.
    ! really want the ocean mask on ocean grid or ocean mask mapped to atm grid,
    ! but do not have access to the ocean mask mapped to the atm grid.
    ! the dom mask is a good place to start, on ocean grid, it should be what we want,
    ! on the atm grid, it's just all 1's so not very useful.
    ! next look at ofrac+ifrac in fractions.  want to compute on all non-land points.
    ! using ofrac alone will exclude points that are currently all sea ice but that later
    ! could be less that 100% covered in ice.

    lsize = size(mask)

    write(tmpstr,'(3L4,i12)') cold_start,first_call,lsize
    call ESMF_LogWrite(trim(subname)//" : coldstart,firstcall,lsize= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    write(tmpstr,'(i12,g22.12,i12)') lsize,sum(rmask),sum(mask)
    call ESMF_LogWrite(trim(subname)//" : maskA= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    mask(:) = 1
    where (rmask(:) == 0._r8) mask(:) = 0   ! like nint

    write(tmpstr,'(i12,g22.12,i12)') lsize,sum(rmask),sum(mask)
    call ESMF_LogWrite(trim(subname)//" : maskB= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    !----------------------------------
    ! Update ocean surface fluxes
    !----------------------------------

    do n = 1,lsize
       if (mask(n) /= 0) then
          !--- mask missing atm or ocn data if found
          if (dens(n) < 1.0e-12 .or. tocn(n) < 1.0) then
             mask(n) = 0
          endif

          !!uGust(n) = 1.5_r8*sqrt(uocn(n)**2 + vocn(n)**2) ! there is no wind gust data from ocn
          uGust(n) = 0.0_r8
          prec(n) = rainc(n) + rainl(n) + snowc(n) + snowl(n)
          prec_gust(n) = rainc(n)
          ! Note: swdn and swup are set in flux_ocnalb using data from previous timestep
       end if
       sen(n)  = shr_const_spval
       lat(n)  = shr_const_spval
       lwup(n) = shr_const_spval
       evap(n) = shr_const_spval
       taux(n) = shr_const_spval
       tauy(n) = shr_const_spval
       tref(n) = shr_const_spval
       qref(n) = shr_const_spval
    end do

    if (do_flux_diurnal) then
       do n = 1,lsize
          nInc(n) = 0._r8 ! needed for minval/maxval calculation
       end do
       call shr_flux_atmocn_diurnal (lsize , zbot , ubot, vbot, thbot, &
            shum , shum_16O , shum_HDO, shum_18O, dens , tbot, uocn, vocn , &
            tocn , mask, sen , lat , lwup , &
            roce_16O, roce_HDO, roce_18O,    &
            evap , evap_16O, evap_HDO, evap_18O, taux , tauy, tref, qref , &
            uGust, lwdn , swdn , swup, prec, &
            fswpen, ocnsal, ocn_prognostic, do_flux_diurnal,    &
            lats, lons , warm , salt , speed, regime,       &
            warmMax, windMax, qSolAvg, windAvg,             &
            warmMaxInc, windMaxInc, qSolInc, windInc, nInc, &
            tbulk, tskin, tskin_day, tskin_night, &
            cskin, cskin_night, tod, dt,          &
            duu10n,ustar, re  , ssq, &
            !missval should not be needed if flux calc
            !consistent with mrgx2a fraction
            !duu10n,ustar, re  , ssq, missval = 0.0_r8, &
            cold_start=cold_start)
    else
       write(tmpstr,'(3i12)') lsize,size(mask),sum(mask)
       call ESMF_LogWrite(trim(subname)//" : mask= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

       call shr_flux_atmocn (lsize , zbot , ubot, vbot, thbot, prec_gust, gust_fac, &
            shum , shum_16O , shum_HDO, shum_18O, dens , &
            tbot, uocn, vocn , &
            tocn , mask, sen , lat , lwup , &
            roce_16O, roce_HDO, roce_18O, &
            evap , evap_16O, evap_HDO, evap_18O, taux , tauy, tref, qref , &
            duu10n, ustar, re, ssq)

       !        do n = 1,lsize
       !           write(6,100)'import: n,zbot      = ',n,zbot(n)
       !           write(6,100)'import: n,ubot      = ',n,ubot(n)
       !           write(6,100)'import: n,vbot      = ',n,vbot(n)
       !           write(6,100)'import: n,thbot     = ',n,thbot(n)
       !           write(6,100)'import: n,prec_gust = ',n,prec_gust(n)
       !           write(6,100)'import: n,tocn      = ',n,tocn(n)
       !           write(6,100)'import: n,uocn      = ',n,uocn(n)
       !           write(6,100)'import: n,vocn      = ',n,vocn(n)
       !           write(6,100)'export: n,latent    = ',n,lat(n)
       !           write(6,100)'export: n,sensible  = ',n,lat(n)
       !           write(6,100)'export: n,taux      = ',n,taux(n)
       !           write(6,100)'export: n,tauy      = ',n,tauy(n)
       !        end do
       ! 100    format ('(atmocn_flux) ',a,i8,d21.14)
    endif

    do n = 1,lsize
       if (mask(n) /= 0) then
          u10(n) = sqrt(duu10n(n))
       end if
    enddo

  end subroutine med_aofluxes_run

end module med_phases_aofluxes_mod
