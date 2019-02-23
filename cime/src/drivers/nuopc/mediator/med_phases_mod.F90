module med_phases_mod

  !-----------------------------------------------------------------------------
  ! Mediator Phases
  !-----------------------------------------------------------------------------

  implicit none
  private

  public  :: med_phases_init

!-----------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------

  subroutine med_phases_init(gcomp, llogunit, rc)
    use ESMF                    ,only : ESMF_GridCompGet, ESMF_VMGet, ESMF_LogWrite, ESMF_LogFlush
    use ESMF                    ,only : ESMF_GRIDCOMP, ESMF_VM, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use med_constants_mod            ,only : CL, R8
    use med_constants_mod            ,only : dbug_flag => med_constants_dbug_flag
    use esmFlds                 , only : compatm, complnd, compocn
    use esmFlds                 , only : compice, comprof, compglc
    use esmFlds                 , only : ncomps, compname
    use esmFlds                 , only : flds_scalar_name
    use esmFlds                 , only : fldListFr, fldListTo
    use esmFlds                 , only : fldListMed_aoflux_a
    use esmFlds                 , only : fldListMed_aoflux_o
    use esmFlds                 , only : fldListMed_ocnalb_o
    use shr_nuopc_fldList_mod   , only : shr_nuopc_fldList_GetFldNames
    use shr_nuopc_methods_mod   , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_init
    use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_diagnose
    use shr_nuopc_methods_mod   , only : shr_nuopc_methods_FB_FldChk
    use med_fraction_mod        , only : med_fraction_init
    use med_constants_mod       , only : med_constants_dbug_flag
    use med_constants_mod       , only : med_constants_czero
    use med_merge_mod           , only : med_merge_auto
    use med_map_mod             , only : med_map_FB_Regrid_Norm
    use med_internalstate_mod   , only : InternalState
    use perf_mod                , only : t_startf, t_stopf
    !----------------------------------------------------------
    ! Initialize field bundles, etc. that are needed as part of
    ! the med_phases routines
    !----------------------------------------------------------

    type(ESMF_GridComp)  :: gcomp
    integer, intent(in)  :: llogunit
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)    :: is_local
    type(ESMF_VM)          :: vm
    integer                :: localPet
    integer                :: n, n1, n2, ncomp, nflds
    character(CL), pointer :: fldnames(:)
    logical                       :: mastertask
    character(*)      , parameter :: u_FILE_u  = __FILE__
    character(len=*)  , parameter :: subname="med_phases_init"
    !-----------------------------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 1) then
       call ESMF_LogWrite("Starting to initialize mediator phases", ESMF_LOGMSG_INFO)
       call ESMF_LogFlush()
    endif

    rc = ESMF_SUCCESS

    ! Determine mastertask
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    call ESMF_VMGet(vm, localPet=localPet, rc=rc)
    mastertask = .false.
    if (localPet == 0) mastertask=.true.

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------------------------------
    ! Create FBfrac field bundles and initialize fractions
    !----------------------------------------------------------

    call med_fraction_init(gcomp,rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Initialize field bundles needed for ocn albedo and ocn/atm flux calculations
    !---------------------------------------

    if (is_local%wrap%med_coupling_active(compocn,compatm) .and. &
        is_local%wrap%med_coupling_active(compatm,compocn)) then

       ! NOTE: the NStateImp(compocn) or NStateImp(compatm) used below
       ! rather than NStateExp(n2), since the export state might only
       ! contain control data and no grid information if if the target
       ! component (n2) is not prognostic only receives control data back

       ! Create field bundles for ocean albedo computation

       nflds = size(fldListMed_ocnalb_o%flds)
       allocate(fldnames(nflds))
       call shr_nuopc_fldList_getfldnames(fldListMed_ocnalb_o%flds, fldnames)

       call shr_nuopc_methods_FB_init(is_local%wrap%FBMed_ocnalb_a, flds_scalar_name, &
            STgeom=is_local%wrap%NStateImp(compatm), fieldnamelist=fldnames, name='FBMed_ocnalb_a', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_init(is_local%wrap%FBMed_ocnalb_o, flds_scalar_name, &
            STgeom=is_local%wrap%NStateImp(compocn), fieldnamelist=fldnames, name='FBMed_ocnalb_o', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       deallocate(fldnames)

       ! Create field bundles for ocean/atmosphere flux computation

       nflds = size(fldListMed_aoflux_o%flds)
       allocate(fldnames(nflds))
       call shr_nuopc_fldList_getfldnames(fldListMed_aoflux_a%flds, fldnames)

       call shr_nuopc_methods_FB_init(is_local%wrap%FBMed_aoflux_a, flds_scalar_name, &
            STgeom=is_local%wrap%NStateImp(compatm), fieldnamelist=fldnames, name='FBMed_aoflux_a', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_init(is_local%wrap%FBMed_aoflux_o, flds_scalar_name, &
            STgeom=is_local%wrap%NStateImp(compocn), fieldnamelist=fldnames, name='FBMed_aoflux_o', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       deallocate(fldnames)

    end if

    !----------------------------------------------------------
    ! Create mediator specific field bundles needed in phases routines
    ! TODO: this needs to be filled in
    !----------------------------------------------------------

    ! FBs for lnd <-> glc accumulation and elevation class downscaling
    if (is_local%wrap%comp_present(complnd) .and. is_local%wrap%comp_present(compglc)) then
       ! call shr_nuopc_methods_FB_init(is_local%wrap%FBMed_l2x_to_glc_accum, &
       !      STgeom=is_local%wrap%NStateImp(complnd), fieldnamelist=flds_l2x_to_glc, name='FBMed_l2g_l_accum', rc=rc)
       ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! call shr_nuopc_methods_FB_init(is_local%wrap%FBMed_g2x_to_lnd, &
       !      STgeom=is_local%wrap%NStateImp(complnd), fieldnamelist=flds_g2x_to_lnd, name='FBMed_g2x_to_lnd', rc=rc)
       ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_init

end module med_phases_mod
