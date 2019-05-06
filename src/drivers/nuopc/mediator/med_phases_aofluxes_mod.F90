module med_phases_aofluxes_mod

  use med_constants_mod     , only : R8, CL, CX
  use med_constants_mod     , only : shr_const_zvir
  use med_constants_mod     , only : shr_const_zvir
  use med_constants_mod     , only : shr_const_cpdair
  use med_constants_mod     , only : shr_const_cpvir
  use med_constants_mod     , only : shr_const_karman
  use med_constants_mod     , only : shr_const_g
  use med_constants_mod     , only : shr_const_latvap
  use med_constants_mod     , only : shr_const_latice
  use med_constants_mod     , only : shr_const_stebol
  use med_constants_mod     , only : shr_const_spval
  use med_internalstate_mod , only : mastertask, logunit
  use med_constants_mod     , only : dbug_flag    => med_constants_dbug_flag
  use shr_nuopc_utils_mod   , only : memcheck     => shr_nuopc_memcheck
  use shr_nuopc_utils_mod   , only : chkerr       => shr_nuopc_utils_chkerr
  use shr_nuopc_methods_mod , only : FB_fldchk    => shr_nuopc_methods_FB_FldChk
  use shr_nuopc_methods_mod , only : FB_GetFldPtr => shr_nuopc_methods_FB_GetFldPtr
  use shr_nuopc_methods_mod , only : FB_diagnose  => shr_nuopc_methods_FB_diagnose
  use shr_nuopc_methods_mod , only : FB_init      => shr_nuopc_methods_FB_init
  use water_isotopes        , only : wiso_flxoce ! calculate water isotope fluxes.

  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public routines
  !--------------------------------------------------------------------------

  public  :: med_phases_aofluxes_run

  !--------------------------------------------------------------------------
  ! Private routines
  !--------------------------------------------------------------------------

  private :: med_aofluxes_init
  private :: med_aofluxes_run
  private :: med_aoflux_adjust_constants ! adjust constant values used in flux calculations.
  private :: med_aoflux_compute          ! computes atm/ocn fluxes

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
     real(R8) , pointer :: pbot        (:) ! atm bottom pressure
     real(R8) , pointer :: qbot        (:) ! atm bottom specific humidity
     real(R8) , pointer :: dens        (:) ! atm bottom density
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
     real(R8) , pointer :: ustar       (:) ! saved ustar
     real(R8) , pointer :: re          (:) ! saved re
     real(R8) , pointer :: ssq         (:) ! saved sq

     ! Fields that are not obtained via GetFldPtr
     logical            :: created         ! has this data type been created
  end type aoflux_type

  ! The follow variables are not declared as parameters so that they can be
  ! adjusted to support aquaplanet and potentially other simple model modes.
  ! The shr_flux_adjust_constants subroutine is called to set the desired
  ! values.  The default values are from shr_const_mod.  Currently they are
  ! only used by the shr_flux_atmocn and shr_flux_atmice routines.

  real(R8) :: loc_zvir   = shr_const_zvir
  real(R8) :: loc_cpdair = shr_const_cpdair
  real(R8) :: loc_cpvir  = shr_const_cpvir
  real(R8) :: loc_karman = shr_const_karman
  real(R8) :: loc_g      = shr_const_g
  real(R8) :: loc_latvap = shr_const_latvap
  real(R8) :: loc_latice = shr_const_latice
  real(R8) :: loc_stebol = shr_const_stebol

  ! These control convergence of the iterative flux calculation
  real(r8) :: flux_con_tol = 0.0_R8
  integer  :: flux_con_max_iter = 2

  ! cold air outbreak parameters  (Mahrt & Sun 1995,MWR)
  logical             :: use_coldair_outbreak_mod = .false.
  real(R8),parameter  :: alpha = 1.4_R8
  real(R8),parameter  :: maxscl =2._R8  ! maximum wind scaling for flux
  real(R8),parameter  :: td0 = -10._R8  ! start t-ts for scaling

  ! The following three variables are obtained as attributes from gcomp
  logical       :: flds_wiso  ! use case
  logical       :: compute_atm_dens
  logical       :: compute_atm_thbot
  character(*), parameter :: u_FILE_u = &
       __FILE__

!================================================================================
contains
!================================================================================

  subroutine med_phases_aofluxes_run(gcomp, rc)

    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_GridCompGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_FieldBundleIsCreated
    use NUOPC                 , only : NUOPC_IsConnected, NUOPC_CompAttributeGet
    use med_internalstate_mod , only : InternalState
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use esmFlds               , only : shr_nuopc_fldList_GetNumFlds, shr_nuopc_fldList_GetFldNames
    use esmFlds               , only : fldListFr, fldListMed_aoflux, compatm, compocn, compname
    use perf_mod              , only : t_startf, t_stopf

    !-----------------------------------------------------------------------
    ! Compute atm/ocn fluxes
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)     :: is_local
    type(aoflux_type), save :: aoflux
    logical, save           :: first_call = .true.
    character(len=*),parameter :: subname='(med_phases_aofluxes)'
    !---------------------------------------

    rc = ESMF_SUCCESS

    ! Get the internal state from the mediator Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (first_call) then
       ! If field bundles have been created for the ocean/atmosphere flux computation
       if ( ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_a, rc=rc) .and. &
            ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_o, rc=rc)) then

          ! Allocate memoroy for the aoflux module data type (mediator atm/ocn field bundle on the ocean grid)
          call med_aofluxes_init(gcomp, aoflux, &
               FBAtm=is_local%wrap%FBImp(compatm,compocn), &
               FBOcn=is_local%wrap%FBImp(compocn,compocn), &
               FBFrac=is_local%wrap%FBfrac(compocn), &
               FBMed_aoflux=is_local%wrap%FBMed_aoflux_o, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
          aoflux%created = .true.
       else
          aoflux%created = .false.
       end if

       ! Now set first_call to .false.
       first_call = .false.
    end if

    ! Return if there is no aoflux has not been created
    if (.not. aoflux%created) then
       RETURN
    end if

    ! Start time timer
    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif

    call memcheck(subname, 5, mastertask)

    ! TODO(mvertens, 2019-01-12): ONLY regrid atm import fields that are needed for the atm/ocn flux calculation

    ! Regrid atm import field bundle from atm to ocn grid as input for ocn/atm flux calculation
    call med_map_FB_Regrid_Norm( &
         fldListFr(compatm)%flds, compatm, compocn, &
         is_local%wrap%FBImp(compatm,compatm), &
         is_local%wrap%FBImp(compatm,compocn), &
         is_local%wrap%FBFrac(compatm), &
         is_local%wrap%FBFrac(compocn), &
         is_local%wrap%FBNormOne(compatm,compocn,:), &
         is_local%wrap%RH(compatm,compocn,:), &
         string=trim(compname(compatm))//'2'//trim(compname(compocn)), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Calculate atm/ocn fluxes on the destination grid
    call med_aofluxes_run(gcomp, aoflux, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 1) then
       call FB_diagnose(is_local%wrap%FBMed_aoflux_o, &
            string=trim(subname) //' FBAMed_aoflux_o' , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_phases_aofluxes_run

!================================================================================

  subroutine med_aofluxes_init(gcomp, aoflux, FBAtm, FBOcn, FBFrac, FBMed_aoflux, rc)

    use ESMF     , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_LogFoundError
    use ESMF     , only : ESMF_SUCCESS, ESMF_LOGERR_PASSTHRU
    use ESMF     , only : ESMF_GridComp, ESMF_GridCompGet, ESMF_VM
    use ESMF     , only : ESMF_Field, ESMF_FieldGet, ESMF_FieldBundle, ESMF_VMGet
    use NUOPC    , only : NUOPC_CompAttributeGet
    use perf_mod , only : t_startf, t_stopf

    !-----------------------------------------------------------------------
    ! Initialize pointers to the module variables
    !-----------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp)                    :: gcomp
    type(aoflux_type)      , intent(inout) :: aoflux
    type(ESMF_FieldBundle) , intent(in)    :: FBAtm               ! Atm Import fields on aoflux grid
    type(ESMF_FieldBundle) , intent(in)    :: FBOcn               ! Ocn Import fields on aoflux grid
    type(ESMF_FieldBundle) , intent(in)    :: FBfrac              ! Fraction data for various components, on their grid
    type(ESMF_FieldBundle) , intent(inout) :: FBMed_aoflux        ! Ocn albedos computed in mediator
    integer                , intent(out)   :: rc

    ! local variables
    integer                  :: iam
    integer                  :: n
    integer                  :: lsize
    real(R8), pointer        :: ofrac(:)
    real(R8), pointer        :: ifrac(:)
    character(CL)            :: cvalue
    logical                  :: flds_wiso  ! use case
    character(len=CX)        :: tmpstr
    character(*),parameter   :: subName =   '(med_aofluxes_init) '
    !-----------------------------------------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO)
    endif
    rc = ESMF_SUCCESS
    call memcheck(subname, 5, mastertask)

    !----------------------------------
    ! get attributes that are set as module variables
    !----------------------------------

    call NUOPC_CompAttributeGet(gcomp, name='flds_wiso', value=cvalue, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    read(cvalue,*) flds_wiso

    !----------------------------------
    ! atm/ocn fields
    !----------------------------------

    call FB_GetFldPtr(FBMed_aoflux, fldname='So_tref', fldptr1=aoflux%tref, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBMed_aoflux, fldname='So_qref', fldptr1=aoflux%qref, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBMed_aoflux, fldname='So_ustar', fldptr1=aoflux%ustar, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBMed_aoflux, fldname='So_re', fldptr1=aoflux%re, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBMed_aoflux, fldname='So_ssq', fldptr1=aoflux%ssq, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBMed_aoflux, fldname='So_u10', fldptr1=aoflux%u10, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBMed_aoflux, fldname='So_duu10n', fldptr1=aoflux%duu10n, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call FB_GetFldPtr(FBMed_aoflux, fldname='Faox_taux', fldptr1=aoflux%taux, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBMed_aoflux, fldname='Faox_tauy', fldptr1=aoflux%tauy, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBMed_aoflux, fldname='Faox_lat', fldptr1=aoflux%lat, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBMed_aoflux, fldname='Faox_sen', fldptr1=aoflux%sen, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBMed_aoflux, fldname='Faox_evap', fldptr1=aoflux%evap, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    lsize = size(aoflux%evap)
    if (flds_wiso) then
       call FB_GetFldPtr(FBMed_aoflux, fldname='Faox_evap_16O', fldptr1=aoflux%evap_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBMed_aoflux, fldname='Faox_evap_18O', fldptr1=aoflux%evap_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBMed_aoflux, fldname='Faox_evap_HDO', fldptr1=aoflux%evap_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%evap_16O(lsize)); aoflux%evap_16O(:) = 0._R8
       allocate(aoflux%evap_18O(lsize)); aoflux%evap_18O(:) = 0._R8
       allocate(aoflux%evap_HDO(lsize)); aoflux%evap_HDO(:) = 0._R8
    end if

    call FB_GetFldPtr(FBMed_aoflux, fldname='Faox_lwup', fldptr1=aoflux%lwup, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Ocn import fields
    !----------------------------------

    call FB_GetFldPtr(FBOcn, fldname='So_omask', fldptr1=aoflux%rmask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBOcn, fldname='So_t', fldptr1=aoflux%tocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBOcn, fldname='So_u', fldptr1=aoflux%uocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBOcn, fldname='So_v', fldptr1=aoflux%vocn, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (flds_wiso) then
       call FB_GetFldPtr(FBOcn, fldname='So_roce_16O', fldptr1=aoflux%roce_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBOcn, fldname='So_roce_18O', fldptr1=aoflux%roce_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBOcn, fldname='So_roce_HDO', fldptr1=aoflux%roce_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%roce_16O(lsize)); aoflux%roce_16O(:) = 0._R8
       allocate(aoflux%roce_18O(lsize)); aoflux%roce_18O(:) = 0._R8
       allocate(aoflux%roce_HDO(lsize)); aoflux%roce_HDO(:) = 0._R8
    end if

    !----------------------------------
    ! Atm import fields
    !----------------------------------

    call FB_GetFldPtr(FBAtm, fldname='Sa_z', fldptr1=aoflux%zbot, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBAtm, fldname='Sa_u', fldptr1=aoflux%ubot, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBAtm, fldname='Sa_v', fldptr1=aoflux%vbot, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call FB_GetFldPtr(FBAtm, fldname='Sa_tbot', fldptr1=aoflux%tbot, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! bottom level potential temperature will need to be computed if not received from the atm
    if (FB_fldchk(FBAtm, 'Sa_ptem', rc=rc)) then
       call FB_GetFldPtr(FBAtm, fldname='Sa_ptem', fldptr1=aoflux%thbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       compute_atm_thbot = .false.
    else
       allocate(aoflux%thbot(lsize))
       compute_atm_thbot = .true.
    end if

    ! bottom level density will need to be computed if not received from the atm
    if (FB_fldchk(FBAtm, 'Sa_dens', rc=rc)) then
       call FB_GetFldPtr(FBAtm, fldname='Sa_dens', fldptr1=aoflux%dens, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       compute_atm_dens = .false.
    else
       compute_atm_dens = .true.
       allocate(aoflux%dens(lsize))
    end if

    ! if either density or potential temperature are computed, will need bottom level pressure
    if (compute_atm_dens .or. compute_atm_thbot) then
       call FB_GetFldPtr(FBAtm, fldname='Sa_pbot', fldptr1=aoflux%pbot, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    call FB_GetFldPtr(FBAtm, fldname='Sa_shum', fldptr1=aoflux%shum, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (flds_wiso) then
       call FB_GetFldPtr(FBAtm, fldname='Sa_shum_16O', fldptr1=aoflux%shum_16O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBAtm, fldname='Sa_shum_18O', fldptr1=aoflux%shum_18O, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(FBAtm, fldname='Sa_shum_HDO', fldptr1=aoflux%shum_HDO, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    else
       allocate(aoflux%shum_16O(lsize)); aoflux%shum_16O(:) = 0._R8
       allocate(aoflux%shum_18O(lsize)); aoflux%shum_18O(:) = 0._R8
       allocate(aoflux%shum_HDO(lsize)); aoflux%shum_HDO(:) = 0._R8
    end if

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
    ! call FB_getFldPtr(FBFrac , fldname='ofrac' , fldptr1=ofrac, rc=rc)
    ! if (chkerr(rc,__LINE__,u_FILE_u)) return
    ! call FB_getFldPtr(FBFrac , fldname='ifrac' , fldptr1=ifrac, rc=rc)
    ! if (chkerr(rc,__LINE__,u_FILE_u)) return
    ! where (ofrac(:) + ifrac(:) <= 0.0_R8) mask(:) = 0

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_aofluxes_init

!===============================================================================

  subroutine med_aofluxes_run(gcomp, aoflux, rc)

    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time, ESMF_TimeInterval
    use ESMF                  , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_TimeGet, ESMF_TimeIntervalGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LogMsg_Info
    use NUOPC                 , only : NUOPC_CompAttributeGet
    use perf_mod              , only : t_startf, t_stopf

    !-----------------------------------------------------------------------
    ! Determine atm/ocn fluxes eother on atm or on ocean grid
    ! The module arrays are set via pointers the the mediator internal states
    ! in med_ocnatm_init and are used below.
    !-----------------------------------------------------------------------

    ! Arguments
    type(ESMF_GridComp)               :: gcomp
    type(aoflux_type) , intent(inout) :: aoflux
    integer           , intent(out)   :: rc
    !
    ! Local variables
    character(CL)           :: cvalue
    integer                 :: n,i                     ! indices
    integer                 :: lsize                   ! local size
    real(R8)                :: gust_fac = huge(1.0_R8) ! wind gust factor
    real(R8)                :: flux_convergence        ! convergence criteria for imlicit flux computation
    integer                 :: flux_max_iteration      ! maximum number of iterations for convergence
    logical                 :: coldair_outbreak_mod    ! cold air outbreak adjustment  (Mahrt & Sun 1995,MWR)
    character(len=CX)       :: tmpstr
    logical,save            :: first_call = .true.
    character(*),parameter  :: subName = '(med_aofluxes_run) '
    !-----------------------------------------------------------------------

    call t_startf('MED:'//subname)

    !----------------------------------
    ! Get config variables on first call
    !----------------------------------

    if (first_call) then
       call NUOPC_CompAttributeGet(gcomp, name='gust_fac', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) gust_fac

       call NUOPC_CompAttributeGet(gcomp, name='coldair_outbreak_mod', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) coldair_outbreak_mod

       call NUOPC_CompAttributeGet(gcomp, name='flux_max_iteration', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flux_max_iteration

       call NUOPC_CompAttributeGet(gcomp, name='flux_convergence', value=cvalue, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       read(cvalue,*) flux_convergence

       call med_aoflux_adjust_constants(&
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

    write(tmpstr,'(3i12)') lsize,size(aoflux%mask),sum(aoflux%mask)
    call ESMF_LogWrite(trim(subname)//" : mask= "//trim(tmpstr), ESMF_LOGMSG_INFO, rc=rc)

    !----------------------------------
    ! Update atmosphere/ocean surface fluxes
    !----------------------------------

    if (compute_atm_thbot) then
       do n = 1,lsize
          if (aoflux%mask(n) /= 0._r8) then
             aoflux%thbot(n) = aoflux%tbot(n)*((100000._R8/aoflux%pbot(n))**0.286_R8)
          end if
       end do
    end if
    if (compute_atm_dens) then
       do n = 1,lsize
          if (aoflux%mask(n) /= 0._r8) then
             aoflux%dens(n) = aoflux%pbot(n)/(287.058_R8*(1._R8 + 0.608_R8*aoflux%shum(n))*aoflux%tbot(n))
          end if
       end do
    end if

    call med_aoflux_compute (&
         lsize, aoflux%zbot, aoflux%ubot, aoflux%vbot, aoflux%thbot, &
         aoflux%shum, aoflux%shum_16O, aoflux%shum_HDO, aoflux%shum_18O, aoflux%dens , &
         aoflux%tbot, aoflux%uocn, aoflux%vocn, &
         aoflux%tocn, aoflux%mask, aoflux%sen, aoflux%lat, aoflux%lwup, &
         aoflux%roce_16O, aoflux%roce_HDO, aoflux%roce_18O, &
         aoflux%evap, aoflux%evap_16O, aoflux%evap_HDO, aoflux%evap_18O, &
         aoflux%taux, aoflux%tauy, aoflux%tref, aoflux%qref, &
         aoflux%duu10n, ustar_sv=aoflux%ustar, re_sv=aoflux%re, ssq_sv=aoflux%ssq, &
         missval=0.0_r8, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do n = 1,lsize
       if (aoflux%mask(n) /= 0) then
          aoflux%u10(n) = sqrt(aoflux%duu10n(n))
       end if
    enddo
    call t_stopf('MED:'//subname)

  end subroutine med_aofluxes_run

  !===============================================================================

  subroutine med_aoflux_adjust_constants( zvir, cpair, cpvir, karman, gravit, &
       latvap, latice, stebol, flux_convergence_tolerance, &
       flux_convergence_max_iteration, coldair_outbreak_mod)

    ! Adjust local constants.  Used to support simple models.

    real(R8)    , optional, intent(in) :: zvir
    real(R8)    , optional, intent(in) :: cpair
    real(R8)    , optional, intent(in) :: cpvir
    real(R8)    , optional, intent(in) :: karman
    real(R8)    , optional, intent(in) :: gravit
    real(R8)    , optional, intent(in) :: latvap
    real(R8)    , optional, intent(in) :: latice
    real(R8)    , optional, intent(in) :: stebol
    real(r8)    , optional, intent(in) :: flux_convergence_tolerance
    integer     , optional, intent(in) :: flux_convergence_max_iteration
    logical     , optional, intent(in) :: coldair_outbreak_mod
    !----------------------------------------------------------------------------

    if (present(zvir))   loc_zvir   = zvir
    if (present(cpair))  loc_cpdair = cpair
    if (present(cpvir))  loc_cpvir  = cpvir
    if (present(karman)) loc_karman = karman
    if (present(gravit)) loc_g      = gravit
    if (present(latvap)) loc_latvap = latvap
    if (present(latice)) loc_latice = latice
    if (present(stebol)) loc_stebol = stebol
    if (present(flux_convergence_tolerance     )) flux_con_tol = flux_convergence_tolerance
    if (present(flux_convergence_max_iteration )) flux_con_max_iter = flux_convergence_max_iteration
    if (present(coldair_outbreak_mod           )) use_coldair_outbreak_mod = coldair_outbreak_mod

  end subroutine med_aoflux_adjust_constants

  !===============================================================================

  subroutine med_aoflux_compute(nMax  ,zbot  ,ubot  ,vbot  ,thbot ,  &
                                qbot  ,s16O  ,sHDO  ,s18O  ,rbot  ,  &
                                tbot  ,us    ,vs    ,                &
                                ts    ,mask  ,sen   ,lat   ,lwup  ,  &
                                r16O, rhdo, r18O,                    &
                                evap  ,evap_16O, evap_HDO, evap_18O, &
                                taux  ,tauy  ,tref  ,qref  ,         &
                                duu10n,  ustar_sv   ,re_sv ,ssq_sv,  &
                                missval, rc)

    !------------------------------------
    ! Internal atm/ocn flux calculation
    !------------------------------------

    use ESMF, only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LogWrite

    ! input/output variables
    integer    ,intent(in) ::       nMax       ! data vector length
    integer    ,intent(in) :: mask (nMax)      ! ocn domain mask        0 <=> out of domain
    real(R8)   ,intent(in) :: zbot (nMax)      ! atm level height       (m)
    real(R8)   ,intent(in) :: ubot (nMax)      ! atm u wind             (m/s)
    real(R8)   ,intent(in) :: vbot (nMax)      ! atm v wind             (m/s)
    real(R8)   ,intent(in) :: thbot(nMax)      ! atm potential T        (K)
    real(R8)   ,intent(in) :: qbot (nMax)      ! atm specific humidity  (kg/kg)
    real(R8)   ,intent(in) :: s16O (nMax)      ! atm H216O tracer conc. (kg/kg)
    real(R8)   ,intent(in) :: sHDO (nMax)      ! atm HDO tracer conc.   (kg/kg)
    real(R8)   ,intent(in) :: s18O (nMax)      ! atm H218O tracer conc. (kg/kg)
    real(R8)   ,intent(in) :: r16O (nMax)      ! ocn H216O tracer ratio/Rstd
    real(R8)   ,intent(in) :: rHDO (nMax)      ! ocn HDO tracer ratio/Rstd
    real(R8)   ,intent(in) :: r18O (nMax)      ! ocn H218O tracer ratio/Rstd
    real(R8)   ,intent(in) :: rbot (nMax)      ! atm air density        (kg/m^3)
    real(R8)   ,intent(in) :: tbot (nMax)      ! atm T                  (K)
    real(R8)   ,intent(in) :: us   (nMax)      ! ocn u-velocity         (m/s)
    real(R8)   ,intent(in) :: vs   (nMax)      ! ocn v-velocity         (m/s)
    real(R8)   ,intent(in) :: ts   (nMax)      ! ocn temperature        (K)

    ! output variables
    real(R8),intent(out)  ::  sen  (nMax)     ! heat flux: sensible        (W/m^2)
    real(R8),intent(out)  ::  lat  (nMax)     ! heat flux: latent          (W/m^2)
    real(R8),intent(out)  ::  lwup (nMax)     ! heat flux: lw upward       (W/m^2)
    real(R8),intent(out)  ::  evap (nMax)     ! water flux: evap           ((kg/s)/m^2)
    real(R8),intent(out)  ::  evap_16O (nMax) ! water flux: evap           ((kg/s/m^2)
    real(R8),intent(out)  ::  evap_HDO (nMax) ! water flux: evap           ((kg/s)/m^2)
    real(R8),intent(out)  ::  evap_18O (nMax) ! water flux: evap           ((kg/s/m^2)
    real(R8),intent(out)  ::  taux (nMax)     ! surface stress, zonal      (N)
    real(R8),intent(out)  ::  tauy (nMax)     ! surface stress, maridional (N)
    real(R8),intent(out)  ::  tref (nMax)     ! diag:  2m ref height T     (K)
    real(R8),intent(out)  ::  qref (nMax)     ! diag:  2m ref humidity     (kg/kg)
    real(R8),intent(out)  :: duu10n(nMax)     ! diag: 10m wind speed squared (m/s)^2
    integer ,intent(out)  :: rc

    real(R8),intent(out),optional :: ustar_sv(nMax) ! diag: ustar
    real(R8),intent(out),optional :: re_sv   (nMax) ! diag: sqrt of exchange coefficient (water)
    real(R8),intent(out),optional :: ssq_sv  (nMax) ! diag: sea surface humidity  (kg/kg)
    real(R8),intent(in) ,optional :: missval        ! masked value

    !--- local constants --------------------------------
    real(R8),parameter :: umin  =  0.5_R8 ! minimum wind speed       (m/s)
    real(R8),parameter :: zref  = 10.0_R8 ! reference height           (m)
    real(R8),parameter :: ztref =  2.0_R8 ! reference height for air T (m)

    !--- local variables --------------------------------
    integer     :: n      ! vector loop index
    integer     :: iter
    real(R8)    :: vmag   ! surface wind magnitude   (m/s)
    real(R8)    :: ssq    ! sea surface humidity     (kg/kg)
    real(R8)    :: delt   ! potential T difference   (K)
    real(R8)    :: delq   ! humidity difference      (kg/kg)
    real(R8)    :: stable ! stability factor
    real(R8)    :: rdn    ! sqrt of neutral exchange coeff (momentum)
    real(R8)    :: rhn    ! sqrt of neutral exchange coeff (heat)
    real(R8)    :: ren    ! sqrt of neutral exchange coeff (water)
    real(R8)    :: rd     ! sqrt of exchange coefficient (momentum)
    real(R8)    :: rh     ! sqrt of exchange coefficient (heat)
    real(R8)    :: re     ! sqrt of exchange coefficient (water)
    real(R8)    :: ustar  ! ustar
    real(r8)    :: ustar_prev
    real(R8)    :: qstar  ! qstar
    real(R8)    :: tstar  ! tstar
    real(R8)    :: hol    ! H (at zbot) over L
    real(R8)    :: xsq    ! ?
    real(R8)    :: xqq    ! ?
    real(R8)    :: psimh  ! stability function at zbot (momentum)
    real(R8)    :: psixh  ! stability function at zbot (heat and water)
    real(R8)    :: psix2  ! stability function at ztref reference height
    real(R8)    :: alz    ! ln(zbot/zref)
    real(R8)    :: al2    ! ln(zref/ztref)
    real(R8)    :: u10n   ! 10m neutral wind
    real(R8)    :: tau    ! stress at zbot
    real(R8)    :: cp     ! specific heat of moist air
    real(R8)    :: fac    ! vertical interpolation factor
    real(R8)    :: spval  ! local missing value

    !--- local functions --------------------------------
    real(R8)    :: qsat   ! function: the saturation humididty of air (kg/m^3)
    real(R8)    :: cdn    ! function: neutral drag coeff at 10m
    real(R8)    :: psimhu ! function: unstable part of psimh
    real(R8)    :: psixhu ! function: unstable part of psimx
    real(R8)    :: Umps   ! dummy arg ~ wind velocity (m/s)
    real(R8)    :: Tk     ! dummy arg ~ temperature (K)
    real(R8)    :: xd     ! dummy arg ~ ?
    real(R8)    :: gprec  ! dummy arg ~ ?
    !--- for cold air outbreak calc --------------------------------
    real(R8)    :: tdiff(nMax)               ! tbot - ts
    real(R8)    :: vscl

    qsat(Tk)   = 640380.0_R8 / exp(5107.4_R8/Tk)
    cdn(Umps)  =   0.0027_R8 / Umps + 0.000142_R8 + 0.0000764_R8 * Umps
    psimhu(xd) = log((1.0_R8+xd*(2.0_R8+xd))*(1.0_R8+xd*xd)/8.0_R8) - 2.0_R8*atan(xd) + 1.571_R8
    psixhu(xd) = 2.0_R8 * log((1.0_R8 + xd*xd)/2.0_R8)

    !--- formats ----------------------------------------
    character(*),parameter :: subName = '(shr_flux_atmOcn) '
    character(*),parameter ::   F00 = "('(shr_flux_atmOcn) ',4a)"
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    if (present(missval)) then
       spval = missval
    else
       spval = shr_const_spval
    endif
    u10n = spval
    rh = spval
    psixh = spval
    hol=spval

    !--- for cold air outbreak calc --------------------------------
    tdiff= tbot - ts

    al2 = log(zref/ztref)
    DO n=1,nMax
       if (mask(n) /= 0) then

          ! compute some needed quantities
          vmag   = max(umin, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) )

          ! Cold Air Outbreak Modification: increase windspeed for negative tbot-ts
          ! based on Mahrt & Sun 1995,MWR
          if (use_coldair_outbreak_mod) then
             if (tdiff(n).lt.td0) then
                vscl=min((1._R8+alpha*(abs(tdiff(n)-td0)**0.5_R8/abs(vmag))),maxscl)
                vmag=vmag*vscl
             endif
          endif

          ssq    = 0.98_R8 * qsat(ts(n)) / rbot(n)   ! sea surf hum (kg/kg)
          delt   = thbot(n) - ts(n)                  ! pot temp diff (K)
          delq   = qbot(n) - ssq                     ! spec hum dif (kg/kg)
          alz    = log(zbot(n)/zref)
          cp     = loc_cpdair*(1.0_R8 + loc_cpvir*ssq)

          !------------------------------------------------------------
          ! first estimate of Z/L and ustar, tstar and qstar
          !------------------------------------------------------------
          !--- neutral coefficients, z/L = 0.0 ---
          stable = 0.5_R8 + sign(0.5_R8 , delt)
          rdn    = sqrt(cdn(vmag))
          rhn    = (1.0_R8-stable) * 0.0327_R8 + stable * 0.018_R8
          ren    = 0.0346_R8

          !--- ustar, tstar, qstar ---
          ustar = rdn * vmag
          tstar = rhn * delt
          qstar = ren * delq
          ustar_prev = ustar*2.0_R8
          iter = 0
          do while( abs((ustar - ustar_prev)/ustar) > flux_con_tol .and. iter < flux_con_max_iter)
             iter = iter + 1
             ustar_prev = ustar
             !--- compute stability & evaluate all stability functions ---
             hol  = loc_karman*loc_g*zbot(n)*  &
                  (tstar/thbot(n)+qstar/(1.0_R8/loc_zvir+qbot(n)))/ustar**2
             hol  = sign( min(abs(hol),10.0_R8), hol )
             stable = 0.5_R8 + sign(0.5_R8 , hol)
             xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
             xqq    = sqrt(xsq)
             psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
             psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)

             !--- shift wind speed using old coefficient ---
             rd   = rdn / (1.0_R8 + rdn/loc_karman*(alz-psimh))
             u10n = vmag * rd / rdn

             !--- update transfer coeffs at 10m and neutral stability ---
             rdn = sqrt(cdn(u10n))
             ren = 0.0346_R8
             rhn = (1.0_R8-stable)*0.0327_R8 + stable * 0.018_R8

             !--- shift all coeffs to measurement height and stability ---
             rd = rdn / (1.0_R8 + rdn/loc_karman*(alz-psimh))
             rh = rhn / (1.0_R8 + rhn/loc_karman*(alz-psixh))
             re = ren / (1.0_R8 + ren/loc_karman*(alz-psixh))

             !--- update ustar, tstar, qstar using updated, shifted coeffs --
             ustar = rd * vmag
             tstar = rh * delt
             qstar = re * delq
          enddo
          if (iter < 1) then
             write(logunit,*) ustar,ustar_prev,flux_con_tol,flux_con_max_iter
             call ESMF_LogWrite('No iterations performed - ERROR in med_calc_aofluxe')
             rc=ESMF_Failure
             return
          end if

          !------------------------------------------------------------
          ! compute the fluxes
          !------------------------------------------------------------

          tau = rbot(n) * ustar * ustar

          !--- momentum flux ---
          taux(n) = tau * (ubot(n)-us(n)) / vmag
          tauy(n) = tau * (vbot(n)-vs(n)) / vmag

          !--- heat flux ---
          sen (n) =          cp * tau * tstar / ustar
          lat (n) =  loc_latvap * tau * qstar / ustar
          lwup(n) = -loc_stebol * ts(n)**4

          !--- water flux ---
          evap(n) = lat(n)/loc_latvap

          !---water isotope flux ---
          call wiso_flxoce(2,rbot(n),zbot(n),s16O(n),ts(n),r16O(n),ustar,re,ssq,evap_16O(n), &
               qbot(n),evap(n))
          call wiso_flxoce(3,rbot(n),zbot(n),sHDO(n),ts(n),rHDO(n),ustar,re,ssq, evap_HDO(n),&
               qbot(n),evap(n))
          call wiso_flxoce(4,rbot(n),zbot(n),s18O(n),ts(n),r18O(n),ustar,re,ssq, evap_18O(n), &
               qbot(n),evap(n))

          !------------------------------------------------------------
          ! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
          !------------------------------------------------------------
          hol = hol*ztref/zbot(n)
          xsq = max( 1.0_R8, sqrt(abs(1.0_R8-16.0_R8*hol)) )
          xqq = sqrt(xsq)
          psix2   = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
          fac     = (rh/loc_karman) * (alz + al2 - psixh + psix2 )
          tref(n) = thbot(n) - delt*fac
          tref(n) = tref(n) - 0.01_R8*ztref   ! pot temp to temp correction
          fac     = (re/loc_karman) * (alz + al2 - psixh + psix2 )
          qref(n) =  qbot(n) - delq*fac

          duu10n(n) = u10n*u10n ! 10m wind speed squared

          !------------------------------------------------------------
          ! optional diagnostics, needed for water tracer fluxes (dcn)
          !------------------------------------------------------------
          if (present(ustar_sv)) ustar_sv(n) = ustar
          if (present(re_sv   )) re_sv(n)    = re
          if (present(ssq_sv  )) ssq_sv(n)   = ssq

       else
          !------------------------------------------------------------
          ! no valid data here -- out of domain
          !------------------------------------------------------------
          sen   (n) = spval  ! sensible         heat flux  (W/m^2)
          lat   (n) = spval  ! latent           heat flux  (W/m^2)
          lwup  (n) = spval  ! long-wave upward heat flux  (W/m^2)
          evap  (n) = spval  ! evaporative water flux ((kg/s)/m^2)
          evap_16O (n) = spval !water tracer flux (kg/s)/m^2)
          evap_HDO (n) = spval !HDO tracer flux  (kg/s)/m^2)
          evap_18O (n) = spval !H218O tracer flux (kg/s)/m^2)
          taux  (n) = spval  ! x surface stress (N)
          tauy  (n) = spval  ! y surface stress (N)
          tref  (n) = spval  !  2m reference height temperature (K)
          qref  (n) = spval  !  2m reference height humidity (kg/kg)
          duu10n(n) = spval  ! 10m wind speed squared (m/s)^2

          if (present(ustar_sv)) ustar_sv(n) = spval
          if (present(re_sv   )) re_sv   (n) = spval
          if (present(ssq_sv  )) ssq_sv  (n) = spval
       endif
    end DO

  end subroutine med_aoflux_compute

end module med_phases_aofluxes_mod
