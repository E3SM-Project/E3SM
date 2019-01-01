module med_phases_prep_ocn_mod

  use med_constants_mod, only : dbug_flag=>med_constants_dbug_flag
  use shr_nuopc_utils_mod, only : shr_nuopc_memcheck
  use med_internalstate_mod, only : mastertask
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

    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_TimeGet, ESMF_ClockPrint
    use ESMF                  , only : ESMF_FieldBundleGet
    use med_internalstate_mod , only : InternalState
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
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
    !--- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    ! Note - the scalar field has been removed from all mediator field bundles - so this is why we check if the
    ! fieldCount is 0 and not 1 here

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compocn), fieldCount=ncnt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt == 0) then
       call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBexp(compocn), returning", &
            ESMF_LOGMSG_INFO, rc=dbrc)
    else

       !---------------------------------------
       !--- map all fields in FBImp that have
       !--- active ocean coupling to the ocean grid
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
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FldChk
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use med_constants_mod     , only : R8
    use med_internalstate_mod , only : InternalState, mastertask, logunit
    use med_merge_mod         , only : med_merge_auto
    use esmFlds               , only : fldListMed_ocnalb_o
    use esmFlds               , only : fldListTo
    use esmFlds               , only : compocn, compname, compatm, compice
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)         :: is_local
    integer                     :: n, ncnt
    real(R8)                    :: c1,c2,c3,c4
    real(R8), pointer           :: dataptr1(:)
    real(R8), pointer           :: ifrac(:), ofrac(:)
    real(R8), pointer           :: ifracr(:), ofracr(:)
    real(R8), pointer           :: avsdr(:), avsdf(:)
    real(R8), pointer           :: anidr(:), anidf(:)
    real(R8), pointer           :: swvdf(:), swndf(:)
    real(R8), pointer           :: swvdr(:), swndr(:)
    real(R8), pointer           :: Foxx_swnet(:)
    real(R8), pointer           :: Foxx_swnet_vdr(:)
    real(R8), pointer           :: Foxx_swnet_vdf(:)
    real(R8), pointer           :: Foxx_swnet_idr(:)
    real(R8), pointer           :: Foxx_swnet_idf(:)
    real(R8), pointer           :: Fioi_swpen(:)
    real(R8), pointer           :: Fioi_swpen_vdr(:)
    real(R8), pointer           :: Fioi_swpen_vdf(:)
    real(R8), pointer           :: Fioi_swpen_idr(:)
    real(R8), pointer           :: Fioi_swpen_idf(:)
    real(R8), pointer           :: latent(:)
    real(R8), pointer           :: evap(:)
    real(R8)                    :: ifrac_scaled, ofrac_scaled
    real(R8)                    :: ifracr_scaled, ofracr_scaled
    real(R8)                    :: frac_sum
    real(R8)                    :: fswabsv, fswabsi
    real(R8)                    :: flux_epbalfact
    logical                     :: compute_ocnalb_in_med
    logical                     :: compute_aoflux_in_med
    logical                     :: compute_evap_in_med
    logical                     :: export_swnet_by_bands
    logical                     :: import_swpen_by_bands
    logical                     :: first_call = .true.
    integer                     :: lsize
    integer                     :: dbrc
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
    !--- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    ! Note - the scalar field has been removed from all mediator field bundles - so this is why we check if the
    ! fieldCount is 0 and not 1 here

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compocn), fieldCount=ncnt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt == 0) then
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBexp(compocn), returning", &
               ESMF_LOGMSG_INFO, rc=dbrc)
       endif
    else
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

       if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compocn), 'Foxx_rain', rc=rc)) then
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rain' , dataptr1, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr1(:) = dataptr1(:) * flux_epbalfact
          if (first_call .and. mastertask) then
             write(logunit,'(a)')'(merge_to_ocn): Scaling Foxx_rain by flux_epbalfact '
          end if
       end if
       if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compocn), 'Foxx_snow', rc=rc)) then
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_snow' , dataptr1, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr1(:) = dataptr1(:) * flux_epbalfact
          if (first_call .and. mastertask) then
             write(logunit,'(a)')'(merge_to_ocn): Scaling Foxx_snow by flux_epbalfact '
          end if
       end if
       if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compocn), 'Foxx_prec', rc=rc)) then
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_prec' , dataptr1, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr1(:) = dataptr1(:) * flux_epbalfact
          if (first_call .and. mastertask) then
             write(logunit,'(a)')'(merge_to_ocn): Scaling Foxx_prec by flux_epbalfact '
          end if
       end if
       if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compocn), 'Foxx_rofl', rc=rc)) then
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rofl' , dataptr1, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr1(:) = dataptr1(:) * flux_epbalfact
          if (first_call .and. mastertask) then
             write(logunit,'(a)')'(merge_to_ocn): Scaling Foxx_rofl by flux_epbalfact '
          end if
       end if
       if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compocn), 'Foxx_rofi', rc=rc)) then
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_rofi' , dataptr1, rc=rc)
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

       ! ----------------
       ! Input from atm
       ! ----------------

       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swvdr', swvdr, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swndr', swndr, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swvdf', swvdf, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swndf', swndf, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       lsize = size(swvdr)

       ! ----------------
       ! Input from mediator
       ! ----------------

       ! get ice-covered ocean and open ocean fractions
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrac' , ifrac, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac' , ofrac, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! determine if ocean albedos are computed in mediator
       compute_ocnalb_in_med = ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_ocnalb_o, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       if (compute_ocnalb_in_med) then
          ! ocean albedos are computed in mediator
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_avsdr' , avsdr, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_anidr' , anidr, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_avsdf' , avsdf, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_anidf' , anidf, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrad' , ifracr, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ofrad' , ofracr, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! ----------------
       ! Input from ice
       ! ----------------

       if (is_local%wrap%comp_present(compice)) then
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen', Fioi_swpen, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_vdr', rc=rc)) then
             import_swpen_by_bands = .true.
             call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_vdr', Fioi_swpen_vdr, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_vdf', Fioi_swpen_vdf, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_idr', Fioi_swpen_idr, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_idf', Fioi_swpen_idf, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             import_swpen_by_bands = .false.
         end if
       end if

       ! ----------------
       ! Output to ocean
       ! ----------------

       if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compocn), 'Foxx_swnet', rc=rc)) then
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet',  Foxx_swnet, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          lsize = size(swvdr)
          allocate(Foxx_swnet(lsize))
       end if

       if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr', rc=rc)) then
          export_swnet_by_bands = .true.
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr', Foxx_swnet_vdr, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdf', Foxx_swnet_vdf, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idr', Foxx_swnet_idr, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idf', Foxx_swnet_idf, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          export_swnet_by_bands = .false.
       end if

       do n = 1,lsize
          ! Compute total swnet to ocean independent of swpen from sea-ice
          if (compute_ocnalb_in_med) then
             fswabsv  = swvdr(n) * (1.0_R8 - avsdr(n)) + swvdf(n) * (1.0_R8 - avsdf(n))
             fswabsi  = swndr(n) * (1.0_R8 - anidr(n)) + swndf(n) * (1.0_R8 - anidf(n))
          else
             fswabsv  = swvdr(n) * (1.0_R8 - albdif) + swvdf(n) * (1.0_R8 - albdif)
             fswabsi  = swndr(n) * (1.0_R8 - albdif) + swndf(n) * (1.0_R8 - albdif)
          end if
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
                   Foxx_swnet_vdr(n) = swvdr(n)*(1.0_R8-avsdr(n))*ofracr_scaled + Fioi_swpen_vdr(n)*ifrac_scaled
                   Foxx_swnet_vdf(n) = swvdf(n)*(1.0_R8-avsdf(n))*ofracr_scaled + Fioi_swpen_vdf(n)*ifrac_scaled
                   Foxx_swnet_idr(n) = swndr(n)*(1.0_R8-avsdr(n))*ofracr_scaled + Fioi_swpen_idr(n)*ifrac_scaled
                   Foxx_swnet_idf(n) = swndf(n)*(1.0_R8-avsdf(n))*ofracr_scaled + Fioi_swpen_idf(n)*ifrac_scaled
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
       ! determine evaporation to send to ocean
       !-------------

       ! In both cesm and nems, evaporation is not sent by atm component
       ! In cesm, its computed in the mediator in the atm/ocn flux calculation
       ! In nems (i.e. fv3 is the atm), it will be computed here using the merged latent heat flux
       ! that is sent to the ocean
       ! Note - don't need to scale the calculated evap by ofrac - since the merged latent heat
       ! to the ocean has already had this scaling done

       ! determine if evaporation need to be computed outside of mediator aoflux computation
       compute_evap_in_med  = (ESMF_FieldBundleIsCreated(is_local%wrap%FBMed_aoflux_o, rc=rc))
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       if (.not. compute_evap_in_med) then
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_lat', latent, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_evap', evap, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          ! TODO (mvertens, 2018-12-16): is this the right sign? Minus here is based on nems mediator
          do n = 1, size(evap)
             evap(n) = - latent(n)/const_lhvap
          end do
       end if

       ! TODO (mvertens, 2018-12-16): document above custom calculation

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(compocn), string=trim(subname)//' FBexp(compocn) ', rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

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

    use ESMF                  , only: ESMF_GridComp, ESMF_Clock, ESMF_Time
    use ESMF                  , only: ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only: ESMF_GridCompGet, ESMF_ClockGet, ESMF_TimeGet, ESMF_ClockPrint
    use ESMF                  , only: ESMF_FieldBundleGet
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_accum
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
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
    !--- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    ! Note - the scalar field has been removed from all mediator field bundles - so this is why we check if the
    ! fieldCount is 0 and not 1 here

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compocn), fieldCount=ncnt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt == 0) then
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBexp(compocn), returning", &
               ESMF_LOGMSG_INFO, rc=dbrc)
       endif
    else

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
    !--- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    ! Note - the scalar field has been removed from all mediator field bundles - so this is why we check if the
    ! fieldCount is 0 and not 1 here

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(compocn), fieldCount=ncnt, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt == 0) then
       if (dbug_flag > 5) then
          call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBexp(compocn), returning", &
               ESMF_LOGMSG_INFO, rc=dbrc)
       endif
       RETURN
    end if

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

    call shr_nuopc_methods_FB_copy(is_local%wrap%FBExp(compocn), &
         is_local%wrap%FBExpAccum(compocn), rc=rc)
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
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ocn_accum_avg

end module med_phases_prep_ocn_mod
