module med_phases_prep_ocn_mod

  use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag
  use shr_nuopc_utils_mod   , only : shr_nuopc_memcheck
  use med_internalstate_mod , only : mastertask

  !-----------------------------------------------------------------------------
  ! Carry out fast accumulation for the ocean
  !-----------------------------------------------------------------------------

  implicit none
  private

  public :: med_phases_prep_ocn_map
  public :: med_phases_prep_ocn_merge
  public :: med_phases_prep_ocn_accum_fast
  public :: med_phases_prep_ocn_accum_avg

  character(*), parameter :: u_FILE_u  = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_prep_ocn_map(gcomp, rc)

    !---------------------------------------
    ! Map all fields in from relevant source components to the ocean grid
    !---------------------------------------

    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO,ESMF_SUCCESS
    use ESMF                  , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_TimeGet, ESMF_ClockPrint
    use ESMF                  , only : ESMF_FieldBundleGet, ESMF_FieldBundleIsCreated
    use med_internalstate_mod , only : InternalState
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getNumFlds
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use esmFlds               , only : fldListFr
    use esmFlds               , only : compocn, ncomps, compname
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)         :: is_local
    integer                     :: n1, ncnt
    integer                     :: dbrc
    character(len=*), parameter :: subname='(med_phases_prep_ocn_map)'
    !-------------------------------------------------------------------------------

    call t_startf('MED:'//subname)
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)
    rc = ESMF_SUCCESS
    call shr_nuopc_memcheck(subname, 5, mastertask)

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------
    call shr_nuopc_methods_FB_getNumFlds(is_local%wrap%FBExp(compocn), trim(subname)//"FBexp(compocn)", ncnt, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       !---------------------------------------
       !--- map all fields in FBImp that have active ocean coupling to the ocean grid
       !---------------------------------------

       do n1 = 1,ncomps
          if (is_local%wrap%med_coupling_active(n1,compocn)) then
             call med_map_FB_Regrid_Norm( &
                  fldListFr(n1)%flds, n1, compocn, &
                  is_local%wrap%FBImp(n1,n1), &
                  is_local%wrap%FBImp(n1,compocn), &
                  is_local%wrap%FBFrac(n1), &
                  is_local%wrap%FBNormOne(n1,compocn,:), &
                  is_local%wrap%RH(n1,compocn,:), &
                  string=trim(compname(n1))//'2'//trim(compname(compocn)), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       enddo
    endif

    call t_stopf('MED:'//subname)
    call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine med_phases_prep_ocn_map

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_ocn_merge(gcomp, rc)

    use ESMF                  , only : ESMF_GridComp, ESMF_FieldBundleGet, ESMF_FieldBundleIsCreated
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_FAILURE,  ESMF_LOGMSG_ERROR
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : fldchk => shr_nuopc_methods_FB_FldChk
    use shr_nuopc_methods_mod , only : FB_GetFldPtr => shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getNumFlds
    use med_constants_mod     , only : R8
    use med_internalstate_mod , only : InternalState, mastertask, logunit
    use med_merge_mod         , only : med_merge_auto, med_merge_field
    use esmFlds               , only : fldListTo
    use esmFlds               , only : compocn, compname, compatm, compice
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n, ncnt
    real(R8)            :: c1,c2,c3,c4
    real(R8), pointer   :: dataptr1(:)
    real(R8), pointer   :: ifrac(:), ofrac(:)
    real(R8), pointer   :: ifracr(:), ofracr(:)
    real(R8), pointer   :: avsdr(:), avsdf(:)
    real(R8), pointer   :: anidr(:), anidf(:)
    real(R8), pointer   :: Faxa_swvdf(:), Faxa_swndf(:)
    real(R8), pointer   :: Faxa_swvdr(:), Faxa_swndr(:)
    real(R8), pointer   :: Foxx_swnet(:)
    real(R8), pointer   :: Foxx_swnet_vdr(:), Foxx_swnet_vdf(:)
    real(R8), pointer   :: Foxx_swnet_idr(:), Foxx_swnet_idf(:)
    real(R8), pointer   :: Fioi_swpen_vdr(:), Fioi_swpen_vdf(:)
    real(R8), pointer   :: Fioi_swpen_idr(:), Fioi_swpen_idf(:)
    real(R8), pointer   :: Fioi_swpen(:)
    real(R8), pointer   :: Foxx_evap(:)
    real(R8), pointer   :: Foxx_lwnet(:)
    real(R8), pointer   :: Faox_lwup(:)
    real(R8), pointer   :: Faxa_lwdn(:)
    real(R8)            :: ifrac_scaled, ofrac_scaled
    real(R8)            :: ifracr_scaled, ofracr_scaled
    real(R8)            :: frac_sum
    real(R8)            :: albvis_dir, albvis_dif
    real(R8)            :: albnir_dir, albnir_dif
    real(R8)            :: fswabsv, fswabsi
    real(R8)            :: flux_epbalfact
    logical             :: compute_ocnalb_in_med
    logical             :: compute_aoflux_in_med
    logical             :: export_swnet_by_bands
    logical             :: import_swpen_by_bands
    logical             :: nems_orig
    logical             :: first_call = .true.
    integer             :: lsize
    integer             :: dbrc
    ! NEMS-orig
    real(R8), pointer   :: atmwgt(:)
    real(R8), pointer   :: icewgt(:)
    real(R8), pointer   :: customwgt(:)
    real(R8), pointer   :: atmwgt1(:)
    real(R8), pointer   :: icewgt1(:)
    real(R8), pointer   :: wgtp01(:)
    real(R8), pointer   :: wgtm01(:)
    !
    real(R8)        , parameter :: const_lhvap = 2.501e6_R8  ! latent heat of evaporation ~ J/kg
    real(R8)        , parameter :: albdif = 0.06_r8          ! 60 deg reference albedo, diffuse
    character(len=*), parameter :: subname='(med_phases_prep_ocn_merge)'
    !---------------------------------------

    call t_startf('MED:'//subname)
    call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)
    rc = ESMF_SUCCESS
    call shr_nuopc_memcheck(subname, 5, mastertask)

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    call shr_nuopc_methods_FB_getNumFlds(is_local%wrap%FBExp(compocn), trim(subname)//"FBexp(compocn)", ncnt, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt >= 0) then

       !---------------------------------------
       !--- auto merges to ocn
       !---------------------------------------

       compute_aoflux_in_med = (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_o, rc=rc))
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (compute_aoflux_in_med) then
          call med_merge_auto(trim(compname(compocn)), &
               is_local%wrap%FBExp(compocn), is_local%wrap%FBFrac(compocn), &
               is_local%wrap%FBImp(:,compocn), fldListTo(compocn), &
               FBMed1=is_local%wrap%FBMed_aoflux_o, &
               document=first_call, string='(merge_to_ocn)', mastertask=mastertask, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          call med_merge_auto(trim(compname(compocn)), &
               is_local%wrap%FBExp(compocn), is_local%wrap%FBFrac(compocn), &
               is_local%wrap%FBImp(:,compocn), fldListTo(compocn), &
               document=first_call, string='(merge_to_ocn)', mastertask=mastertask, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       !--- custom calculations
       !---------------------------------------

       !-------------
       ! scale precipitation to ocean
       !-------------

       ! TODO (mvertens, 2018-12-16): the calculation needs to be set at run time based on receiving it from the ocean
       flux_epbalfact = 1.0_r8

       if (fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rain', rc=rc)) then
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rain' , dataptr1, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr1(:) = dataptr1(:) * flux_epbalfact
          if (first_call .and. mastertask) then
             write(logunit,'(a)')'(merge_to_ocn): Scaling Foxx_rain by flux_epbalfact '
          end if
       end if
       if (fldchk(is_local%wrap%FBExp(compocn), 'Foxx_snow', rc=rc)) then
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_snow' , dataptr1, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr1(:) = dataptr1(:) * flux_epbalfact
          if (first_call .and. mastertask) then
             write(logunit,'(a)')'(merge_to_ocn): Scaling Foxx_snow by flux_epbalfact '
          end if
       end if
       if (fldchk(is_local%wrap%FBExp(compocn), 'Foxx_prec', rc=rc)) then
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_prec' , dataptr1, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr1(:) = dataptr1(:) * flux_epbalfact
          if (first_call .and. mastertask) then
             write(logunit,'(a)')'(merge_to_ocn): Scaling Foxx_prec by flux_epbalfact '
          end if
       end if
       if (fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofl', rc=rc)) then
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rofl' , dataptr1, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr1(:) = dataptr1(:) * flux_epbalfact
          if (first_call .and. mastertask) then
             write(logunit,'(a)')'(merge_to_ocn): Scaling Foxx_rofl by flux_epbalfact '
          end if
       end if
       if (fldchk(is_local%wrap%FBExp(compocn), 'Foxx_rofi', rc=rc)) then
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rofi' , dataptr1, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr1(:) = dataptr1(:) * flux_epbalfact
          if (first_call .and. mastertask) then
             write(logunit,'(a)')'(merge_to_ocn): Scaling Foxx_rofi by flux_epbalfact '
          end if
       end if

       !-------------
       ! Compute netsw for ocean
       !-------------

       ! netsw_for_ocn = downsw_from_atm * (1-ocn_albedo) * (1-ice_fraction) + pensw_from_ice * (ice_fraction)

       ! Input from atm
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swvdr', Faxa_swvdr, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swndr', Faxa_swndr, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swvdf', Faxa_swvdf, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swndf', Faxa_swndf, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       lsize = size(Faxa_swvdr)

       ! Input from mediator, ice-covered ocean and open ocean fractions
       call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrac' , ifrac, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac' , ofrac, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Input from mediator, ocean albedos
       compute_ocnalb_in_med = ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_o, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (compute_ocnalb_in_med) then
          call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_avsdr' , avsdr, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_anidr' , anidr, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_avsdf' , avsdf, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_anidf' , anidf, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrad' , ifracr, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ofrad' , ofracr, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! Input from ice
       if (is_local%wrap%comp_present(compice)) then
          call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen', Fioi_swpen, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          if (fldchk(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_vdr', rc=rc)) then
             import_swpen_by_bands = .true.
             call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_vdr', Fioi_swpen_vdr, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_vdf', Fioi_swpen_vdf, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_idr', Fioi_swpen_idr, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_idf', Fioi_swpen_idf, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             import_swpen_by_bands = .false.
         end if
       end if

       ! Output to ocean
       if (fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet', rc=rc)) then
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet',  Foxx_swnet, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          write(6,*)'DEBUG: i am here'
       else
          lsize = size(Faxa_swvdr)
          allocate(Foxx_swnet(lsize))
       end if

       if (fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr', rc=rc)) then
          export_swnet_by_bands = .true.
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr', Foxx_swnet_vdr, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdf', Foxx_swnet_vdf, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idr', Foxx_swnet_idr, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idf', Foxx_swnet_idf, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          export_swnet_by_bands = .false.
       end if

       do n = 1,lsize
          ! Determine ocean albedos
          if (compute_ocnalb_in_med) then
             albvis_dir = avsdr(n)
             albvis_dif = avsdf(n)
             albnir_dir = anidr(n)
             albnir_dif = anidf(n)
          else
             albvis_dir = albdif
             albvis_dif = albdif
             albnir_dir = albdif
             albnir_dif = albdif
          end if

          ! Compute total swnet to ocean independent of swpen from sea-ice
          fswabsv  = Faxa_swvdr(n) * (1.0_R8 - albvis_dir) + Faxa_swvdf(n) * (1.0_R8 - albvis_dif)
          fswabsi  = Faxa_swndr(n) * (1.0_R8 - albnir_dir) + Faxa_swndf(n) * (1.0_R8 - albnir_dif)
          Foxx_swnet(n) = fswabsv + fswabsi

          ! Add swpen from sea ice if sea ice is present
          if (is_local%wrap%comp_present(compice)) then
             if (compute_ocnalb_in_med) then
                ifrac_scaled = ifrac(n)
                ofrac_scaled = ofrac(n)
                frac_sum = ifrac(n) + ofrac(n)
                if (frac_sum /= 0._R8) then
                   ifrac_scaled = ifrac(n) / (frac_sum)
                   ofrac_scaled = ofrac(n) / (frac_sum)
                endif
                ifracr_scaled = ifracr(n)
                ofracr_scaled = ofracr(n)
                frac_sum = ifracr(n) + ofracr(n)
                if (frac_sum /= 0._R8) then
                   ifracr_scaled = ifracr(n) / (frac_sum)
                   ofracr_scaled = ofracr(n) / (frac_sum)
                endif
             else
                ofracr_scaled = ofrac(n)
                ifrac_scaled  = ifrac(n)
             end if
             Foxx_swnet(n) = ofracr_scaled*Foxx_swnet(n) + ifrac_scaled*Fioi_swpen(n)

             if (export_swnet_by_bands) then
                if (import_swpen_by_bands) then
                   ! use each individual band for swpen coming from the sea-ice
                   Foxx_swnet_vdr(n) = Faxa_swvdr(n)*(1.0_R8-albvis_dir)*ofracr_scaled + Fioi_swpen_vdr(n)*ifrac_scaled
                   Foxx_swnet_vdf(n) = Faxa_swvdf(n)*(1.0_R8-albvis_dif)*ofracr_scaled + Fioi_swpen_vdf(n)*ifrac_scaled
                   Foxx_swnet_idr(n) = Faxa_swndr(n)*(1.0_R8-albnir_dir)*ofracr_scaled + Fioi_swpen_idr(n)*ifrac_scaled
                   Foxx_swnet_idf(n) = Faxa_swndf(n)*(1.0_R8-albnir_dif)*ofracr_scaled + Fioi_swpen_idf(n)*ifrac_scaled
                else
                   ! scale total Foxx_swnet to get contributions from each band
                   c1 = 0.285
                   c2 = 0.285
                   c3 = 0.215
                   c4 = 0.215
                   Foxx_swnet_vdr(n) = c1 * Foxx_swnet(n)
                   Foxx_swnet_vdf(n) = c2 * Foxx_swnet(n)
                   Foxx_swnet_idr(n) = c3 * Foxx_swnet(n)
                   Foxx_swnet_idf(n) = c4 * Foxx_swnet(n)
                end if
             end if
             ! TODO (mvertens, 2018-12-16): fill in the following
             ! if (i2o_per_cat) then
             !   Sf_ofrac(n)  = ofrac(n)
             !   Sf_ofracr(n) = ofracr(n)
             !   Foxx_swnet_ofracr(n) = (fswabsv + fswabsi) * ofracr_scaled
             ! end if
          end if  ! if sea-ice is present
       end do

       !-------------
       ! determine if mediator was NEMS-orig coupling
       !-------------

       if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Foxx_sen'  , rc=rc) .and. &
            fldchk(is_local%wrap%FBMed_aoflux_o         , 'Faox_sen'  , rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_melth', rc=rc) .and. &
            fldchk(is_local%wrap%FBImp(compatm, compocn), 'Faxa_sen'  , rc=rc)) then
          nems_orig = .true.
       else
          nems_orig = .false.
       end if

       !-------------
       ! determine evaporation to send to ocean if not computed in mediator and not NEMS-orig
       !-------------
       ! Note - don't need to scale the calculated evap by ofrac - since the merged latent heat
       ! to the ocean has already had this scaling done

       if (.not. nems_orig) then
          if ( fldchk(is_local%wrap%FBImp(compatm,compatm), 'Faxa_lat' , rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_evap', rc=rc)) then
             allocate(customwgt(lsize))

             ! NEMS-frac
             ! TODO (mvertens, 2018-12-16): is this the right sign? Minus here is based on nems mediator

             customwgt(:) = - 1._r8 / const_lhvap
             call med_merge_field(is_local%wrap%FBExp(compocn), 'Foxx_evap', &
                  FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_lat',  wgtA=customwgt, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             deallocate(customwgt)
          end if
       end if

       !-------------
       ! Field merges for NEMS-orig coupling
       !-------------

       if (nems_orig) then

          allocate(atmwgt(lsize))
          allocate(atmwgt1(lsize))
          allocate(icewgt1(lsize))
          allocate(wgtp01(lsize))
          allocate(wgtm01(lsize))
          allocate(customwgt(lsize))

          do n = 1,lsize
             atmwgt(n)  = ofrac(n)
             atmwgt1(n) = ofrac(n)
             icewgt1(n) = ifrac(n)
             wgtp01(n)  = 0.0_R8
             wgtm01(n)  = 0.0_R8
             if (ifrac(n) <= 0._R8) then
                atmwgt1(n) =  0.0_R8
                icewgt1(n) =  0.0_R8
                wgtp01(n)  =  1.0_R8
                wgtm01(n)  = -1.0_R8
             end if

             ! check wgts do add to 1 as expected
             if (abs(atmwgt(n)  + icewgt(n) - 1.0_R8) > 1.0e-12 .or. &
                 abs(atmwgt1(n) + icewgt1(n) + wgtp01(n) - 1.0_R8) > 1.0e-12 .or. &
                 abs(atmwgt1(n) + icewgt1(n) - wgtm01(n) - 1.0_R8) > 1.0e-12) then

                call ESMF_LogWrite(trim(subname)//": ERROR atm + ice fracs inconsistent", &
                     ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
                rc = ESMF_FAILURE
                return
             endif
          end do

          ! determine evap to to ocean NEMS-orig
          if ( fldchk(is_local%wrap%FBMed_aoflux_o        , 'Faox_evap', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm,compocn), 'Faxa_lat' , rc=rc) .and. &
               fldchk(is_local%wrap%FBExp(compocn)        , 'Foxx_evap', rc=rc)) then

             customwgt(:) = wgtm01(:) / const_lhvap
             call med_merge_field(is_local%wrap%FBExp(compocn), 'Foxx_evap', &
                  FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_evap', wgtA=atmwgt1, &
                  FBinB=is_local%wrap%FBImp(compatm,compocn), fnameB='Faxa_lat' , wgtB=customwgt, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! determine sensible heat flux to ocean - NEMS orig
          if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Foxx_sen'  , rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o         , 'Faox_sen'  , rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_melth', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm, compocn), 'Faxa_sen'  , rc=rc)) then

             call med_merge_field(is_local%wrap%FBExp(compocn), 'Foxx_sen',    &
                  FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_sen '  , wgtA=atmwgt1, &
                  FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_melth' , wgtB=icewgt1, &
                  FBinC=is_local%wrap%FBImp(compice,compocn), fnameC='Faxa_sen' ,   wgtc=wgtm01, rc=rc)
          end if

          ! determine zonal stress to ocean - NEMS orig
          if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Foxx_taux', rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o         , 'Faox_taux', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_taux', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm, compocn), 'Faxa_taux', rc=rc)) then

             call med_merge_field(is_local%wrap%FBExp(compocn), 'Foxx_taux',  &
                  FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_taux ', wgtA=atmwgt1, &
                  FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_taux' , wgtB=icewgt1, &
                  FBinC=is_local%wrap%FBImp(compatm,compocn), fnameC='Faxa_taux' , wgtc=wgtm01, rc=rc)
          end if

          ! determine meridional stress to ocean - NEMS orig
          if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Foxx_tauy', rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o         , 'Faox_tauy', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compice, compice), 'Fioi_tauy', rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm, compocn), 'Faxa_tauy', rc=rc)) then

             call med_merge_field(is_local%wrap%FBExp(compocn), 'Foxx_tauy',  &
                  FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_tauy ', wgtA=atmwgt1, &
                  FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_tauy' , wgtB=icewgt1, &
                  FBinC=is_local%wrap%FBImp(compatm,compocn), fnameC='Faxa_tauy' , wgtc=wgtm01, rc=rc)
          end if

          ! determine net longwave to ocean - NEMS orig
          if ( fldchk(is_local%wrap%FBexp(compocn)         , 'Foxx_lwnet', rc=rc) .and. &
               fldchk(is_local%wrap%FBMed_aoflux_o         , 'Faox_lwup' , rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm, compocn), 'Faxa_lwup' , rc=rc) .and. &
               fldchk(is_local%wrap%FBImp(compatm, compocn), 'Faxa_lwdn' , rc=rc)) then

             call med_merge_field(is_local%wrap%FBExp(compocn), 'Foxx_lwnet', &
                  FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_lwup ', wgtA=atmwgt1, &
                  FBinB=is_local%wrap%FBImp(compatm,compocn), fnameB='Faxa_lwdn' , wgtB=atmwgt1, &
                  FBinC=is_local%wrap%FBImp(compatm,compocn), fnameC='Faxa_lwup' , wgtc=wgtp01, rc=rc)
          end if

          deallocate(atmwgt)
          deallocate(icewgt)
          deallocate(atmwgt1)
          deallocate(icewgt1)
          deallocate(wgtp01)
          deallocate(wgtm01)
          deallocate(customwgt)

       end if  ! end of NEMS-orig ocn prep phase

       !---------------------------------------
       !--- diagnose output
       !---------------------------------------

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(compocn), string=trim(subname)//' FBexp(compocn) ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

       ! TODO (mvertens, 2018-12-16): document above custom calculation

       !---------------------------------------
       !--- clean up
       !---------------------------------------

       first_call = .false.
    endif

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ocn_merge

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_ocn_accum_fast(gcomp, rc)

    ! Carry out fast accumulation for the ocean

    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_TimeGet, ESMF_ClockPrint
    use ESMF                  , only : ESMF_FieldBundleGet, ESMF_FieldBundleIsCreated
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_accum
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getNumFlds
    use med_internalstate_mod , only : InternalState, mastertask
    use esmFlds               , only : compocn
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_Time)             :: time
    character(len=64)           :: timestr
    type(InternalState)         :: is_local
    integer                     :: i,j,n,ncnt
    integer                     :: dbrc
    character(len=*), parameter :: subname='(med_phases_accum_fast)'
    !---------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------
    call shr_nuopc_methods_FB_getNumFlds(is_local%wrap%FBExp(compocn), trim(subname)//"FBexp(compocn)", ncnt, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       !---------------------------------------
       !--- Get the current time from the clock
       !---------------------------------------

       call ESMF_GridCompGet(gcomp, clock=clock)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockGet(clock,currtime=time,rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeGet(time,timestring=timestr)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (dbug_flag > 1) then
          call ESMF_LogWrite(trim(subname)//": time = "//trim(timestr), ESMF_LOGMSG_INFO, rc=dbrc)
       endif
#if DEBUG
       if (mastertask) then
          call ESMF_ClockPrint(clock, options="currTime", preString="-------->"//trim(subname)//" mediating for: ", rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
#endif
       !---------------------------------------
       !--- ocean accumulator
       !---------------------------------------

       call shr_nuopc_methods_FB_accum(is_local%wrap%FBExpAccum(compocn), is_local%wrap%FBExp(compocn), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       is_local%wrap%FBExpAccumCnt(compocn) = is_local%wrap%FBExpAccumCnt(compocn) + 1

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExpAccum(compocn), &
               string=trim(subname)//' FBaccOcn_AFaccum ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

       !---------------------------------------
       !--- clean up
       !---------------------------------------
    endif

    call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ocn_accum_fast

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_ocn_accum_avg(gcomp, rc)

    ! Prepare the OCN import Fields.

    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_FieldBundleGet
    use med_constants_mod     , only : czero=>med_constants_czero
    use med_internalstate_mod , only : InternalState
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_average
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_copy
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getNumFlds
    use esmFlds               , only : compocn
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)           :: clock
    type(ESMF_Time)            :: time
    character(len=64)          :: timestr
    type(InternalState)        :: is_local
    integer                    :: i,j,n,ncnt
    integer                    :: dbrc
    character(len=*),parameter :: subname='(med_phases_prep_ocn_accum_avg)'
    !---------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------
    call shr_nuopc_methods_FB_getNumFlds(is_local%wrap%FBExpAccum(compocn), trim(subname)//"FBExpAccum(compocn)", ncnt, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       !---------------------------------------
       !--- average ocn accumulator
       !---------------------------------------

       if (dbug_flag > 5) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExpAccum(compocn), &
               string=trim(subname)//' FBExpAccum(compocn) before avg ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

       call shr_nuopc_methods_FB_average(is_local%wrap%FBExpAccum(compocn), &
            is_local%wrap%FBExpAccumCnt(compocn), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 5) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExpAccum(compocn), &
               string=trim(subname)//' FBaccO_avg ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

       !---------------------------------------
       !--- copy to FBExp(compocn)
       !---------------------------------------

       call shr_nuopc_methods_FB_copy(is_local%wrap%FBExp(compocn), is_local%wrap%FBExpAccum(compocn), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       !--- zero accumulator
       !---------------------------------------

       is_local%wrap%FBExpAccumFlag(compocn) = .true.
       is_local%wrap%FBExpAccumCnt(compocn) = 0
       call shr_nuopc_methods_FB_reset(is_local%wrap%FBExpAccum(compocn), value=czero, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 5) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExpAccum(compocn), &
               string=trim(subname)//' FBExpAccum(compocn) after avg ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

    end if

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ocn_accum_avg

end module med_phases_prep_ocn_mod
