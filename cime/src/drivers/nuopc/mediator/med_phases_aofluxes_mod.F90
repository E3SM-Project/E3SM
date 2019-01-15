module med_phases_aofluxes_mod

  use med_constants_mod , only : R8
  use med_constants_mod , only : dbug_flag => med_constants_dbug_flag
  use shr_nuopc_utils_mod, only : shr_nuopc_memcheck
  use med_internalstate_mod, only : mastertask

  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public routines
  !--------------------------------------------------------------------------

  public  :: med_phases_aofluxes_run

  !--------------------------------------------------------------------------
  ! Private routines
  !--------------------------------------------------------------------------

  private :: med_phases_aofluxes_init
  private :: med_aofluxes_init
  private :: med_aofluxes_run

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  type aoflux_type
     integer  , pointer :: mask        (:) ! ocn domain mask: 0 <=> inactive cell
     real(R8) , pointer :: rmask       (:) ! ocn domain mask: 0 <=> inactive cell
     real(R8) , pointer :: lats        (:) ! latitudes  (degrees)
     real(R8) , pointer :: lons        (:) ! longitudes (degrees)
     real(R8) , pointer :: uocn        (:) ! ocn velocity, zonal
     real(R8) , pointer :: vocn        (:) ! ocn velocity, meridional
     real(R8) , pointer :: tocn        (:) ! ocean temperature
     real(R8) , pointer :: zbot        (:) ! atm level height
     real(R8) , pointer :: ubot        (:) ! atm velocity, zonal
     real(R8) , pointer :: vbot        (:) ! atm velocity, meridional
     real(R8) , pointer :: thbot       (:) ! atm potential T
     real(R8) , pointer :: shum        (:) ! atm specific humidity
     real(R8) , pointer :: shum_16O    (:) ! atm H2O tracer
     real(R8) , pointer :: shum_HDO    (:) ! atm HDO tracer
     real(R8) , pointer :: shum_18O    (:) ! atm H218O tracer
     real(R8) , pointer :: roce_16O    (:) ! ocn H2O ratio
     real(R8) , pointer :: roce_HDO    (:) ! ocn HDO ratio
     real(R8) , pointer :: roce_18O    (:) ! ocn H218O ratio
     real(R8) , pointer :: dens        (:) ! atm density
     real(R8) , pointer :: tbot        (:) ! atm bottom surface T
     real(R8) , pointer :: sen         (:) ! heat flux: sensible
     real(R8) , pointer :: lat         (:) ! heat flux: latent
     real(R8) , pointer :: lwup        (:) ! lwup over ocean
     real(R8) , pointer :: evap        (:) ! water flux: evaporation
     real(R8) , pointer :: evap_16O    (:) ! H2O flux: evaporation
     real(R8) , pointer :: evap_HDO    (:) ! HDO flux: evaporation
     real(R8) , pointer :: evap_18O    (:) ! H218O flux: evaporation
     real(R8) , pointer :: taux        (:) ! wind stress, zonal
     real(R8) , pointer :: tauy        (:) ! wind stress, meridional
     real(R8) , pointer :: tref        (:) ! diagnostic:  2m ref T
     real(R8) , pointer :: qref        (:) ! diagnostic:  2m ref Q
     real(R8) , pointer :: u10         (:) ! diagnostic: 10m wind speed
     real(R8) , pointer :: duu10n      (:) ! diagnostic: 10m wind speed squared
     real(R8) , pointer :: lwdn        (:) ! long  wave, downward
     real(R8) , pointer :: swdn        (:) ! short wave, downward
     real(R8) , pointer :: swup        (:) ! short wave, upward
     real(R8) , pointer :: rainc       (:) ! rainc
     real(R8) , pointer :: rainl       (:) ! rainl
     real(R8) , pointer :: snowc       (:) ! snowc
     real(R8) , pointer :: snowl       (:) ! snowl
     real(R8) , pointer :: ustar       (:) ! saved ustar
     real(R8) , pointer :: re          (:) ! saved re
     real(R8) , pointer :: ssq         (:) ! saved sq

     ! Fields that are not obtained via GetFldPtr
     real(R8) , pointer :: uGust       (:) ! wind gust
     real(R8) , pointer :: prec        (:) ! precip
     real(R8) , pointer :: prec_gust   (:) ! atm precip for convective gustiness (kg/m^3)
  end type aoflux_type

  ! The following three variables are obtained as attributes from gcomp
  logical       :: flds_wiso  ! use case
  character(3)  :: aoflux_grid
  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine med_phases_aofluxes_init(gcomp, aoflux, rc)

    ! Initialize ocn/atm flux calculations

    use ESMF                  , only : ESMF_GridComp, ESMF_VM, ESMF_VMGet, ESMF_GridCompGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LOGERR_PASSTHRU
    use ESMF                  , only : ESMF_SUCCESS, ESMF_LogFoundError
    use NUOPC                 , only : NUOPC_CompAttributeGet
    use esmFlds               , only : fldListMed_aoflux_o
    use esmFlds               , only : compatm, compocn
    use med_internalstate_mod , only : InternalState, mastertask
    use shr_nuopc_scalars_mod , only : flds_scalar_name
    use shr_nuopc_scalars_mod , only : flds_scalar_num
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_init
    use shr_nuopc_fldList_mod , only : shr_nuopc_fldlist_getfldnames
    use perf_mod              , only : t_startf, t_stopf
    ! input/output variables
    type(ESMF_GridComp)               :: gcomp
    type(aoflux_type) , intent(inout) :: aoflux
    integer           , intent(out)   :: rc

    ! local variables
    character(3)        :: aoflux_grid
    character(len=256)  :: cvalue
    type(InternalState) :: is_local
    integer             :: nflds
    integer             :: localPet
    type(ESMF_VM)       :: vm
    integer             :: dbrc
    character(len=*),parameter :: subname='(med_phases_aofluxes_init)'
    !---------------------------------------
    call t_startf('MED:'//subname)

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

    ! Determine src and dst comps depending on the aoflux_grid setting

    call NUOPC_CompAttributeGet(gcomp, name='aoflux_grid', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) aoflux_grid

    if (trim(aoflux_grid) == 'ocn') then

       ! Create FBMed_aoflux_o (field bundle on the ocean grid)
       call med_aofluxes_init(gcomp, aoflux, &
            FBAtm=is_local%wrap%FBImp(compatm,compocn), &
            FBOcn=is_local%wrap%FBImp(compocn,compocn), &
            FBFrac=is_local%wrap%FBfrac(compocn), &
            FBMed_ocnalb=is_local%wrap%FBMed_ocnalb_o, &
            FBMed_aoflux=is_local%wrap%FBMed_aoflux_o, &
            rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    else if (trim(aoflux_grid) == 'atm') then

       ! Create FBMed_aoflux_a (field bundle on the atmosphere grid)
       call med_aofluxes_init(gcomp, aoflux, &
            FBAtm=is_local%wrap%FBImp(compatm,compatm), &
            FBOcn=is_local%wrap%FBImp(compocn,compatm), &
            FBFrac=is_local%wrap%FBfrac(compatm), &
            FBMed_ocnalb=is_local%wrap%FBMed_ocnalb_a, &
            FBMed_aoflux=is_local%wrap%FBMed_aoflux_a, &
            rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    else

       call ESMF_LogWrite(trim(subname)//' aoflux_grid = '//trim(aoflux_grid)//' not available', &
            ESMF_LOGMSG_INFO, rc=dbrc)
       return

    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_aofluxes_init

!================================================================================

  subroutine med_phases_aofluxes_run(gcomp, rc)
    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_GridCompGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use NUOPC                 , only : NUOPC_IsConnected, NUOPC_CompAttributeGet
    use med_internalstate_mod , only : InternalState
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use esmFlds               , only : fldListFr
    use esmFlds               , only : compatm, compocn, compname
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_init
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use med_constants_mod     , only : CL
    use perf_mod              , only : t_startf, t_stopf
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)     :: is_local
    type(ESMF_Clock)        :: clock
    character(CL)           :: cvalue
    character(CL)           :: aoflux_grid
    type(aoflux_type), save :: aoflux
    logical, save           :: first_call = .true.
    integer                 :: dbrc
    character(len=*),parameter :: subname='(med_phases_aofluxes)'
    !---------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS
    call shr_nuopc_memcheck(subname, 5, mastertask)
    ! Get the clock from the mediator Component
    call ESMF_GridCompGet(gcomp, clock=clock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from the mediator Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Initialize aoflux instance
    if (first_call) then
       call med_phases_aofluxes_init(gcomp, aoflux, rc)
       first_call = .false.
    end if

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
       call med_aofluxes_run(gcomp, aoflux, rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBMed_aoflux_o, &
               string=trim(subname) //' FBAMed_aoflux_o' , rc=rc)
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
       call med_aofluxes_run(gcomp, aoflux, rc)
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
    call t_stopf('MED:'//subname)

  end subroutine med_phases_aofluxes_run

!================================================================================

  subroutine med_aofluxes_init(gcomp, aoflux, FBAtm, FBOcn, FBFrac, FBMed_ocnalb, FBMed_aoflux, rc)
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LogFoundError
    use ESMF                  , only : ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU, ESMF_FAILURE
    use ESMF                  , only : ESMF_GridComp, ESMF_FieldBundle, ESMF_VM, ESMF_Field
    use ESMF                  , only : ESMF_FieldGet, ESMF_VMGet, ESMF_GridCompGet
    use ESMF                  , only : ESMF_Grid, ESMF_Mesh, ESMF_MeshGet, ESMF_GEOMTYPE_FLAG
    use ESMF                  , only : operator(==), ESMF_GEOMTYPE_MESH
    use NUOPC                 , only : NUOPC_CompAttributeGet
    use med_constants_mod     , only : CL, CX
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use perf_mod              , only : t_startf, t_stopf
    !-----------------------------------------------------------------------
    ! Initialize pointers to the module variables
    !-----------------------------------------------------------------------

    ! Arguments
    type(ESMF_GridComp)                    :: gcomp
    type(aoflux_type)      , intent(inout) :: aoflux
    type(ESMF_FieldBundle) , intent(in)    :: FBAtm               ! Atm Import fields on aoflux grid
    type(ESMF_FieldBundle) , intent(in)    :: FBOcn               ! Ocn Import fields on aoflux grid
    type(ESMF_FieldBundle) , intent(in)    :: FBfrac              ! Fraction data for various components, on their grid
    type(ESMF_FieldBundle) , intent(in)    :: FBMed_ocnalb        ! Ocn albedos computed in mediator
    type(ESMF_FieldBundle) , intent(inout) :: FBMed_aoflux        ! Ocn albedos computed in mediator
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
    real(R8), pointer        :: ofrac(:)
    real(R8), pointer        :: ifrac(:)
    integer                  :: dimCount
    integer                  :: spatialDim
    integer                  :: numOwnedElements
    character(CL)            :: cvalue
    real(R8), pointer        :: ownedElemCoords(:)
    character(*),parameter   :: subName =   '(med_aofluxes_init) '
    logical                  :: flds_wiso  ! use case
    integer                  :: dbrc
    character(len=CX)        :: tmpstr

    !-----------------------------------------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS
    call shr_nuopc_memcheck(subname, 5, mastertask)
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
    ! atm/ocn fields
    !----------------------------------

    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_tref', fldptr1=aoflux%tref, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_qref', fldptr1=aoflux%qref, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_ustar', fldptr1=aoflux%ustar, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_re', fldptr1=aoflux%re, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_ssq', fldptr1=aoflux%ssq, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_u10', fldptr1=aoflux%u10, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='So_duu10n', fldptr1=aoflux%duu10n, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_taux', fldptr1=aoflux%taux, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_tauy', fldptr1=aoflux%tauy, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_lat', fldptr1=aoflux%lat, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_sen', fldptr1=aoflux%sen, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_evap', fldptr1=aoflux%evap, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    lsize = size(aoflux%evap)
    if (flds_wiso) then
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_evap_16O', fldptr1=aoflux%evap_16O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_evap_18O', fldptr1=aoflux%evap_18O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_evap_HDO', fldptr1=aoflux%evap_HDO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%evap_16O(lsize)); aoflux%evap_16O(:) = 0._R8
       allocate(aoflux%evap_18O(lsize)); aoflux%evap_18O(:) = 0._R8
       allocate(aoflux%evap_HDO(lsize)); aoflux%evap_HDO(:) = 0._R8
    end if

    call shr_nuopc_methods_FB_GetFldPtr(FBMed_aoflux, fldname='Faox_lwup', fldptr1=aoflux%lwup, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Ocn import fields
    !----------------------------------

    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_omask', fldptr1=aoflux%rmask, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_t', fldptr1=aoflux%tocn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_u', fldptr1=aoflux%uocn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_v', fldptr1=aoflux%vocn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    if (flds_wiso) then
       call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_roce_16O', fldptr1=aoflux%roce_16O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_roce_18O', fldptr1=aoflux%roce_18O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBOcn, fldname='So_roce_HDO', fldptr1=aoflux%roce_HDO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%roce_16O(lsize)); aoflux%roce_16O(:) = 0._R8
       allocate(aoflux%roce_18O(lsize)); aoflux%roce_18O(:) = 0._R8
       allocate(aoflux%roce_HDO(lsize)); aoflux%roce_HDO(:) = 0._R8
    end if

    !----------------------------------
    ! Atm import fields
    !----------------------------------

    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_z', fldptr1=aoflux%zbot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_u', fldptr1=aoflux%ubot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_v', fldptr1=aoflux%vbot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_tbot', fldptr1=aoflux%tbot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_ptem', fldptr1=aoflux%thbot, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_shum', fldptr1=aoflux%shum, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (flds_wiso) then
       call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_shum_16O', fldptr1=aoflux%shum_16O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_shum_18O', fldptr1=aoflux%shum_18O, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_shum_HDO', fldptr1=aoflux%shum_HDO, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%shum_16O(lsize)); aoflux%shum_16O(:) = 0._R8
       allocate(aoflux%shum_18O(lsize)); aoflux%shum_18O(:) = 0._R8
       allocate(aoflux%shum_HDO(lsize)); aoflux%shum_HDO(:) = 0._R8
    end if

    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Sa_dens', fldptr1=aoflux%dens, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_lwdn', fldptr1=aoflux%lwdn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_rainc', fldptr1=aoflux%rainc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_rainl', fldptr1=aoflux%rainl, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_snowc', fldptr1=aoflux%snowc, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(FBAtm, fldname='Faxa_snowl', fldptr1=aoflux%snowl, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Fields that are not obtained via GetFldPtr
    !----------------------------------
    allocate(aoflux%uGust(lsize))     ; aoflux%uGust(:)     =  0.0_R8
    allocate(aoflux%prec(lsize))      ; aoflux%prec(:)      =  0.0_R8
    allocate(aoflux%prec_gust(lsize)) ; aoflux%prec_gust(:) =  0.0_R8

    !----------------------------------
    ! setup the compute mask.
    !----------------------------------

    ! allocate grid mask fields
    ! default compute everywhere, then "turn off" gridcells
    allocate(aoflux%mask(lsize))
    aoflux%mask(:) = 1

    write(tmpstr,'(i12,g22.12,i12)') lsize,sum(aoflux%rmask),sum(aoflux%mask)
    call ESMF_LogWrite(trim(subname)//" : maskA= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    where (aoflux%rmask(:) == 0._R8) aoflux%mask(:) = 0   ! like nint

    write(tmpstr,'(i12,g22.12,i12)') lsize,sum(aoflux%rmask),sum(aoflux%mask)
    call ESMF_LogWrite(trim(subname)//" : maskB= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    ! TODO: need to check if this logic is correct
    ! then check ofrac + ifrac
    ! call shr_nuopc_methods_FB_getFldPtr(FBFrac , fldname='ofrac' , fldptr1=ofrac, rc=rc)
    ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! call shr_nuopc_methods_FB_getFldPtr(FBFrac , fldname='ifrac' , fldptr1=ifrac, rc=rc)
    ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    ! where (ofrac(:) + ifrac(:) <= 0.0_R8) mask(:) = 0

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_aofluxes_init

!===============================================================================

  subroutine med_aofluxes_run(gcomp, aoflux, rc)
    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time, ESMF_TimeInterval
    use ESMF                  , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_TimeGet, ESMF_TimeIntervalGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LogMsg_Info
    use NUOPC                 , only : NUOPC_CompAttributeGet
    use med_constants_mod     , only : CX, CL
    use shr_const_mod         , only : shr_const_spval
    use shr_flux_mod          , only : shr_flux_atmocn, shr_flux_atmocn_diurnal, shr_flux_adjust_constants
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use perf_mod              , only : t_startf, t_stopf
    !-----------------------------------------------------------------------
    ! Determine atm/ocn fluxes eother on atm or on ocean grid
    ! The module arrays are set via pointers the the mediator internal states
    ! in med_ocnatm_init and are used below.
    ! gcomp (the mediator gridded component) is only needed to retreive the
    ! attributes
    !-----------------------------------------------------------------------

    ! Arguments
    type(ESMF_GridComp)               :: gcomp
    type(aoflux_type) , intent(inout) :: aoflux
    integer           , intent(out)   :: rc
    !
    ! Local variables
    type(ESMF_clock)        :: EClock
    type(ESMF_Time)         :: ETime
    type(ESMF_TimeInterval) :: timeStep
    integer                 :: tod, dt
    character(CL)           :: cvalue
    integer                 :: n,i                     ! indices
    integer                 :: lsize                   ! local size
    real(R8)                :: gust_fac = huge(1.0_R8) ! wind gust factor
    real(R8)                :: flux_convergence        ! convergence criteria for imlicit flux computation
    integer                 :: flux_max_iteration      ! maximum number of iterations for convergence
    logical                 :: coldair_outbreak_mod    ! cold air outbreak adjustment  (Mahrt & Sun 1995,MWR)
    character(len=CX)       :: tmpstr
    logical,save            :: first_call = .true.
    character(*),parameter  :: F01 = "('(med_aofluxes_run) ',a,i4,2x,d21.14)"
    character(*),parameter  :: F02 = "('(med_aofluxes_run) ',a,i4,2x,i4)"
    character(*),parameter  :: subName = '(med_fluxes_run) '
    !-----------------------------------------------------------------------
    call t_startf('MED:'//subname)

    call ESMF_GridCompGet(gcomp, clock=Eclock, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get tod and dt
    call ESMF_ClockGet( Eclock, currTime=ETime, timeStep=timeStep, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeGet( ETime, s=tod, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call ESMF_TimeIntervalGet( timeStep, s=dt, rc=rc )
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    write(tmpstr,'(2i12)') tod,dt
    call ESMF_LogWrite(trim(subname)//" : tod,dt= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

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

    lsize = size(aoflux%mask)

    write(tmpstr,'(i12,g22.12,i12)') lsize,sum(aoflux%rmask),sum(aoflux%mask)
    call ESMF_LogWrite(trim(subname)//" : maskA= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    aoflux%mask(:) = 1
    where (aoflux%rmask(:) == 0._R8) aoflux%mask(:) = 0   ! like nint

    write(tmpstr,'(i12,g22.12,i12)') lsize,sum(aoflux%rmask),sum(aoflux%mask)
    call ESMF_LogWrite(trim(subname)//" : maskB= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    write(tmpstr,'(g22.12,g22.12)') sum(aoflux%dens),sum(aoflux%tocn)
    call ESMF_LogWrite(trim(subname)//" : dens,tocn= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)
    !----------------------------------
    ! Update atmosphere/ocean surface fluxes
    !----------------------------------

    do n = 1,lsize
       if (aoflux%mask(n) /= 0) then
          !--- mask missing atm or ocn data if found
          if (aoflux%dens(n) < 1.0e-12 .or. aoflux%tocn(n) < 1.0) then
             aoflux%mask(n) = 0
          endif
          !!uGust(n) = 1.5_R8*sqrt(uocn(n)**2 + vocn(n)**2) ! there is no wind gust data from ocn
          aoflux%uGust(n) = 0.0_R8
          aoflux%prec(n) = aoflux%rainc(n) + aoflux%rainl(n) + aoflux%snowc(n) + aoflux%snowl(n)
          aoflux%prec_gust(n) = aoflux%rainc(n)
          ! Note: swdn and swup are set in flux_ocnalb using data from previous timestep
       end if
       aoflux%sen(n)  = shr_const_spval
       aoflux%lat(n)  = shr_const_spval
       aoflux%lwup(n) = shr_const_spval
       aoflux%evap(n) = shr_const_spval
       aoflux%taux(n) = shr_const_spval
       aoflux%tauy(n) = shr_const_spval
       aoflux%tref(n) = shr_const_spval
       aoflux%qref(n) = shr_const_spval
    end do

    write(tmpstr,'(3i12)') lsize,size(aoflux%mask),sum(aoflux%mask)
    call ESMF_LogWrite(trim(subname)//" : mask= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    call shr_flux_atmocn (&
         lsize, aoflux%zbot, aoflux%ubot, aoflux%vbot, aoflux%thbot, aoflux%prec_gust, gust_fac, &
         aoflux%shum, aoflux%shum_16O, aoflux%shum_HDO, aoflux%shum_18O, aoflux%dens , &
         aoflux%tbot, aoflux%uocn, aoflux%vocn, &
         aoflux%tocn, aoflux%mask, aoflux%sen, aoflux%lat, aoflux%lwup, &
         aoflux%roce_16O, aoflux%roce_HDO, aoflux%roce_18O, &
         aoflux%evap, aoflux%evap_16O, aoflux%evap_HDO, aoflux%evap_18O, &
         aoflux%taux, aoflux%tauy, aoflux%tref, aoflux%qref, &
         aoflux%duu10n, ustar_sv=aoflux%ustar, re_sv=aoflux%re, ssq_sv=aoflux%ssq)

    do n = 1,lsize
       if (aoflux%mask(n) /= 0) then
          aoflux%u10(n) = sqrt(aoflux%duu10n(n))
       end if
    enddo
    call t_stopf('MED:'//subname)

  end subroutine med_aofluxes_run

end module med_phases_aofluxes_mod
