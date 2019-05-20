module med_phases_prep_ocn_mod

  !-----------------------------------------------------------------------------
  ! Mediator phases for preparing ocn export from mediator
  !-----------------------------------------------------------------------------

  use med_internalstate_mod , only : mastertask
  use med_constants_mod     , only : dbug_flag     => med_constants_dbug_flag
  use shr_nuopc_utils_mod   , only : memcheck      => shr_nuopc_memcheck
  use shr_nuopc_utils_mod   , only : chkerr        => shr_nuopc_utils_ChkErr
  use shr_nuopc_methods_mod , only : FB_diagnose   => shr_nuopc_methods_FB_diagnose
  use shr_nuopc_methods_mod , only : FB_getNumFlds => shr_nuopc_methods_FB_getNumFlds
  use shr_nuopc_methods_mod , only : FB_fldchk     => shr_nuopc_methods_FB_FldChk
  use shr_nuopc_methods_mod , only : FB_GetFldPtr  => shr_nuopc_methods_FB_GetFldPtr
  use shr_nuopc_methods_mod , only : FB_accum      => shr_nuopc_methods_FB_accum
  use shr_nuopc_methods_mod , only : FB_average    => shr_nuopc_methods_FB_average
  use shr_nuopc_methods_mod , only : FB_copy       => shr_nuopc_methods_FB_copy
  use shr_nuopc_methods_mod , only : FB_reset      => shr_nuopc_methods_FB_reset

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
    use ESMF                  , only : ESMF_FieldBundleGet
    use med_internalstate_mod , only : InternalState
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
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    rc = ESMF_SUCCESS
    call memcheck(subname, 5, mastertask)

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------
    call FB_getNumFlds(is_local%wrap%FBExp(compocn), trim(subname)//"FBexp(compocn)", ncnt, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

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
                  is_local%wrap%FBFrac(compocn), &
                  is_local%wrap%FBNormOne(n1,compocn,:), &
                  is_local%wrap%RH(n1,compocn,:), &
                  string=trim(compname(n1))//'2'//trim(compname(compocn)), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       enddo
    endif

    call t_stopf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' done', ESMF_LOGMSG_INFO, rc=dbrc)
    end if

  end subroutine med_phases_prep_ocn_map

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_ocn_merge(gcomp, rc)

    use ESMF                  , only : ESMF_GridComp, ESMF_FieldBundleGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_FAILURE,  ESMF_LOGMSG_ERROR
    use med_constants_mod     , only : R8, CS
    use med_internalstate_mod , only : InternalState, mastertask, logunit
    use med_merge_mod         , only : med_merge_auto, med_merge_field
    use esmFlds               , only : fldListTo
    use esmFlds               , only : compocn, compname, compatm, compice
    use esmFlds               , only : coupling_mode
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: n, ncnt
    real(R8)            :: c1,c2,c3,c4
    real(R8), pointer   :: dataptr(:)
    real(R8), pointer   :: ifrac(:), ofrac(:)
    real(R8), pointer   :: ifracr(:), ofracr(:)
    real(R8), pointer   :: avsdr(:), avsdf(:)
    real(R8), pointer   :: anidr(:), anidf(:)
    real(R8), pointer   :: Faxa_swvdf(:), Faxa_swndf(:)
    real(R8), pointer   :: Faxa_swvdr(:), Faxa_swndr(:)
    real(R8), pointer   :: Foxx_swnet(:)
    real(R8), pointer   :: Foxx_swnet_afracr(:)
    real(R8), pointer   :: Foxx_swnet_vdr(:), Foxx_swnet_vdf(:)
    real(R8), pointer   :: Foxx_swnet_idr(:), Foxx_swnet_idf(:)
    real(R8), pointer   :: Fioi_swpen_vdr(:), Fioi_swpen_vdf(:)
    real(R8), pointer   :: Fioi_swpen_idr(:), Fioi_swpen_idf(:)
    real(R8), pointer   :: Fioi_swpen(:)
    real(R8), pointer   :: Foxx_evap(:)
    real(R8), pointer   :: Foxx_lwnet(:)
    real(R8), pointer   :: Faox_lwup(:)
    real(R8), pointer   :: Faxa_lwdn(:)
    real(R8), pointer   :: dataptr_i(:), dataptr_o(:)
    real(R8), pointer   :: dataptr2d_i(:,:), dataptr2d_o(:,:)
    real(R8)            :: ifrac_scaled, ofrac_scaled
    real(R8)            :: ifracr_scaled, ofracr_scaled
    real(R8)            :: frac_sum
    real(R8)            :: albvis_dir, albvis_dif
    real(R8)            :: albnir_dir, albnir_dif
    real(R8)            :: fswabsv, fswabsi
    logical             :: export_swnet_by_bands
    logical             :: import_swpen_by_bands
    logical             :: export_swnet_afracr
    logical             :: first_precip_fact_call = .true.
    real(R8)            :: precip_fact
    integer             :: lsize
    integer             :: dbrc
    character(CS)       :: cvalue
    ! NEMS-orig
    real(R8), pointer   :: ocnwgt1(:)
    real(R8), pointer   :: icewgt1(:)
    real(R8), pointer   :: wgtp01(:)
    real(R8), pointer   :: wgtm01(:)
    real(R8), pointer   :: customwgt(:)
    !
    character(len=64), allocatable :: fldnames(:)
    real(R8)        , parameter    :: const_lhvap = 2.501e6_R8  ! latent heat of evaporation ~ J/kg
    real(R8)        , parameter    :: albdif = 0.06_r8          ! 60 deg reference albedo, diffuse
    character(len=*), parameter    :: subname='(med_phases_prep_ocn_merge)'
    !---------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(subname//' called', ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    rc = ESMF_SUCCESS
    call memcheck(subname, 5, mastertask)

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    call FB_getNumFlds(is_local%wrap%FBExp(compocn), trim(subname)//"FBexp(compocn)", ncnt, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       !---------------------------------------
       !--- auto merges to ocn
       !---------------------------------------

       if (trim(coupling_mode) == 'cesm' .or. trim(coupling_mode) == 'nems_orig') then
          call med_merge_auto(trim(compname(compocn)), &
               is_local%wrap%FBExp(compocn), is_local%wrap%FBFrac(compocn), &
               is_local%wrap%FBImp(:,compocn), fldListTo(compocn), &
               FBMed1=is_local%wrap%FBMed_aoflux_o, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else if (trim(coupling_mode) == 'nems_frac') then
          call med_merge_auto(trim(compname(compocn)), &
               is_local%wrap%FBExp(compocn), is_local%wrap%FBFrac(compocn), &
               is_local%wrap%FBImp(:,compocn), fldListTo(compocn), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       !--- custom calculations
       !---------------------------------------

       !-------------
       ! Compute netsw for ocean
       !-------------

       ! netsw_for_ocn = downsw_from_atm * (1-ocn_albedo) * (1-ice_fraction) + pensw_from_ice * (ice_fraction)

       ! Input from atm
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swvdr', Faxa_swvdr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swndr', Faxa_swndr, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swvdf', Faxa_swvdf, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), 'Faxa_swndf', Faxa_swndf, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       lsize = size(Faxa_swvdr)

       ! Input from mediator, ice-covered ocean and open ocean fractions
       call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrac' , ifrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac' , ofrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Input from mediator, ocean albedos
       if (trim(coupling_mode) == 'cesm') then
          call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_avsdr' , avsdr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_anidr' , anidr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_avsdf' , avsdf, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, 'So_anidf' , anidf, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ifrad' , ifracr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBfrac(compocn), 'ofrad' , ofracr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! Input from ice
       if (is_local%wrap%comp_present(compice)) then
          call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen', Fioi_swpen, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          if (FB_fldchk(is_local%wrap%FBImp(compice,compice), 'Fioi_swpen_vdr', rc=rc)) then
             import_swpen_by_bands = .true.
             call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_vdr', Fioi_swpen_vdr, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_vdf', Fioi_swpen_vdf, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_idr', Fioi_swpen_idr, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call FB_GetFldPtr(is_local%wrap%FBImp(compice,compocn), 'Fioi_swpen_idf', Fioi_swpen_idf, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          else
             import_swpen_by_bands = .false.
         end if
       end if

       ! Output to ocean swnet 
       if (FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet', rc=rc)) then
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet',  Foxx_swnet, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          lsize = size(Faxa_swvdr)
          allocate(Foxx_swnet(lsize))
       end if

       ! Output to ocean swnet by radiation bands
       if (FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr', rc=rc)) then
          export_swnet_by_bands = .true.
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdr', Foxx_swnet_vdr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_vdf', Foxx_swnet_vdf, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idr', Foxx_swnet_idr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_idf', Foxx_swnet_idf, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       else
          export_swnet_by_bands = .false.
       end if

       ! Swnet without swpen from sea-ice
       if ( FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet_afracr',rc=rc)) then
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Foxx_swnet_afracr', Foxx_swnet_afracr, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          export_swnet_afracr = .true.
       else
          export_swnet_afracr = .false.
       end if

       do n = 1,lsize

          ! Determine ocean albedos
          if (trim(coupling_mode) == 'cesm') then
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

             if (trim(coupling_mode) == 'cesm') then
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

             Foxx_swnet(n) = ofracr_scaled*(fswabsv + fswabsi) + ifrac_scaled*Fioi_swpen(n)

             if (export_swnet_afracr) then
                Foxx_swnet_afracr(n) = ofracr_scaled*(fswabsv + fswabsi) 
             end if

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
          end if  ! if sea-ice is present
       end do

       ! Deallocate Foxx_swnet if it was allocated in this subroutine
       if (.not. FB_fldchk(is_local%wrap%FBExp(compocn), 'Foxx_swnet', rc=rc)) then
          deallocate(Foxx_swnet)
       end if

       ! Output to ocean per ice thickness fraction and sw penetrating into ocean
       if ( FB_fldchk(is_local%wrap%FBExp(compocn), 'Sf_afrac', rc=rc)) then
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Sf_afrac', fldptr1=dataptr_o, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr_o(:) = ofrac(:)
       end if
       if ( FB_fldchk(is_local%wrap%FBExp(compocn), 'Sf_afracr', rc=rc)) then
          call FB_GetFldPtr(is_local%wrap%FBExp(compocn), 'Sf_afracr', fldptr1=dataptr_o, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          dataptr_o(:) = ofracr(:)
       end if

       !-------------
       ! application of precipitation factor from ocean
       !-------------
       precip_fact = 1.0_R8
       if (precip_fact /= 1.0_R8) then
          if (first_precip_fact_call .and. mastertask) then
             write(logunit,'(a)')'(merge_to_ocn): Scaling rain, snow, liquid and ice runoff by precip_fact '
             first_precip_fact_call = .false.
          end if
          write(cvalue,*) precip_fact
          call ESMF_LogWrite(trim(subname)//" precip_fact is "//trim(cvalue), ESMF_LOGMSG_INFO, rc=dbrc)

          allocate(fldnames(4))
          fldnames = (/'Faxa_rain','Faxa_snow', 'Foxx_rofl', 'Foxx_rofi'/)
          do n = 1,size(fldnames)
             if (FB_fldchk(is_local%wrap%FBExp(compocn), trim(fldnames(n)), rc=rc)) then
                call FB_GetFldPtr(is_local%wrap%FBExp(compocn), trim(fldnames(n)) , dataptr, rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
                dataptr(:) = dataptr(:) * precip_fact
             end if
          end do
          deallocate(fldnames)
       end if

       !-------------
       ! custom calculation for nems_frac coupling
       !-------------
       if (trim(coupling_mode) == 'nems_frac') then

          ! determine evaporation to send to ocean 
          ! Note - don't need to scale the calculated evap by ofrac - since the merged latent heat
          ! to the ocean has already had this scaling done
          ! TODO (mvertens, 2018-12-16): is this the right sign below? Minus here is based on nems mediator

          allocate(customwgt(lsize))
          customwgt(:) = - 1._r8 / const_lhvap
          call med_merge_field(is_local%wrap%FBExp(compocn), 'Foxx_evap', &
               FBinA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_lat',  wgtA=customwgt, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          deallocate(customwgt)
       end if

       !-------------
       ! Custom calculation for nems_orig coupling
       !-------------
       if (trim(coupling_mode) == 'nems_orig') then

          ! open ocean (i.e. atm)  and ice fraction
          ! ocnwgt and icewgt are the "normal" fractions
          ! ocnwgt1, icewgt1, and wgtp01 are the fractions that switch between atm and mediator fluxes
          ! wgtp01 and wgtm01 are the same just one is +1 and the other is -1 to change sign depending on the ice fraction.  
          !   ocnwgt1+icewgt1+wgtp01 = 1.0 always 
          !   wgtp01 = 1 and wgtm01 = -1 when ice fraction = 0  
          !   wgtp01 = 0 and wgtm01 =  0 when ice fraction > 0
       
          allocate(ocnwgt1(lsize))
          allocate(icewgt1(lsize))
          allocate(wgtp01(lsize))
          allocate(wgtm01(lsize))
          allocate(customwgt(lsize))

          do n = 1,lsize
             if (ifrac(n) <= 0._R8) then
                ! ice fraction is 0
                ocnwgt1(n) =  0.0_R8
                icewgt1(n) =  0.0_R8
                wgtp01(n)  =  1.0_R8
                wgtm01(n)  = -1.0_R8
             else
                ! ice fraction is > 0
                ocnwgt1(n) = ofrac(n)
                icewgt1(n) = ifrac(n)
                wgtp01(n)  = 0.0_R8
                wgtm01(n)  = 0.0_R8
             end if

             ! check wgts do add to 1 as expected
             if ( abs( ofrac(n) + ifrac(n) - 1.0_R8) > 1.0e-12 .or. &
                  abs( ocnwgt1(n) + icewgt1(n) + wgtp01(n) - 1.0_R8) > 1.0e-12 .or. &
                  abs( ocnwgt1(n) + icewgt1(n) - wgtm01(n) - 1.0_R8) > 1.0e-12) then

                write(6,100)trim(subname)//'ERROR: n, ofrac, ifrac, sum',&
                     n,ofrac(n),ifrac(n),ofrac(n)+ifrac(n)
                write(6,101)trim(subname)//'ERROR: n, ocnwgt1, icewgt1, wgtp01, sum ', &
                     n,ocnwgt1(n),icewgt1(n),wgtp01(n),ocnwgt1(n)+icewgt1(n)+wgtp01(n)  
                write(6,101)trim(subname)//'ERROR: n, ocnwgt1, icewgt1, -wgtm01, sum ', &
                     n,ocnwgt1(n),icewgt1(n),-wgtp01(n),ocnwgt1(n)+icewgt1(n)-wgtm01(n)  
100             format(a,i8,2x,3(d20.13,2x))
101             format(a,i8,2x,4(d20.13,2x))

                call ESMF_LogWrite(trim(subname)//": ERROR atm + ice fracs inconsistent", &
                     ESMF_LOGMSG_ERROR, line=__LINE__, file=__FILE__, rc=dbrc)
                rc = ESMF_FAILURE
                return
             endif
          end do

          customwgt(:) = wgtm01(:) / const_lhvap
          call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_evap', &
               FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_evap', wgtA=ocnwgt1, &
               FBinB=is_local%wrap%FBImp(compatm,compocn), fnameB='Faxa_lat' , wgtB=customwgt, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_sen',    &
               FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_sen '  , wgtA=ocnwgt1, &
               FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_melth' , wgtB=icewgt1, &
               FBinC=is_local%wrap%FBImp(compatm,compocn), fnameC='Faxa_sen'   , wgtc=wgtm01, rc=rc)

          call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_taux',  &
               FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_taux ', wgtA=ocnwgt1, &
               FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_taux' , wgtB=icewgt1, &
               FBinC=is_local%wrap%FBImp(compatm,compocn), fnameC='Faxa_taux' , wgtc=wgtm01, rc=rc)

          call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_tauy',  &
               FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_tauy ', wgtA=ocnwgt1, &
               FBinB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_tauy' , wgtB=icewgt1, &
               FBinC=is_local%wrap%FBImp(compatm,compocn), fnameC='Faxa_tauy' , wgtc=wgtm01, rc=rc)

          ! If there is no ice on the ocn gridcell (ocnwgt1=0) - sum Faxa_lwdn and Faxa_lwup
          ! If there is ice on the ocn gridcell -  merge Faox_lwup and Faxa_lwdn and ignore Faxa_lwup
          call med_merge_field(is_local%wrap%FBExp(compocn),      'Foxx_lwnet', &
               FBinA=is_local%wrap%FBMed_aoflux_o        , fnameA='Faox_lwup ', wgtA=ocnwgt1, &
               FBinB=is_local%wrap%FBImp(compatm,compocn), fnameB='Faxa_lwdn' , wgtB=ocnwgt1, &
               FBinC=is_local%wrap%FBImp(compatm,compocn), fnameC='Faxa_lwnet' , wgtc=wgtp01, rc=rc)

          call med_merge_field(is_local%wrap%FBExp(compocn),      'Faxa_rain' , &
               FBInA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_rain' , wgtA=ofrac, &
               FBInB=is_local%wrap%FBImp(compice,compocn), fnameB='Fioi_meltw', wgtB=ifrac, rc=rc)

          call med_merge_field(is_local%wrap%FBExp(compocn),      'Faxa_snow' , &
               FBInA=is_local%wrap%FBImp(compatm,compocn), fnameA='Faxa_snow' , wgtA=ofrac, rc=rc)

          deallocate(ocnwgt1)
          deallocate(icewgt1)
          deallocate(wgtp01)
          deallocate(wgtm01)
          deallocate(customwgt)

       end if  ! end of NEMS-orig ocn prep phase

       !---------------------------------------
       !--- diagnose output
       !---------------------------------------

       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExp(compocn), &
               string=trim(subname)//' FBexp(compocn) ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       ! TODO (mvertens, 2018-12-16): document above custom calculation

       !---------------------------------------
       !--- clean up
       !---------------------------------------

    endif

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ocn_merge

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_ocn_accum_fast(gcomp, rc)

    ! Carry out fast accumulation for the ocean

    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_TimeGet, ESMF_ClockPrint
    use ESMF                  , only : ESMF_FieldBundleGet
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

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------
    call FB_getNumFlds(is_local%wrap%FBExp(compocn), trim(subname)//"FBexp(compocn)", ncnt, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       !---------------------------------------
       !--- ocean accumulator
       !---------------------------------------

       call FB_accum(is_local%wrap%FBExpAccum(compocn), is_local%wrap%FBExp(compocn), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       is_local%wrap%FBExpAccumCnt(compocn) = is_local%wrap%FBExpAccumCnt(compocn) + 1

       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExpAccum(compocn), &
               string=trim(subname)//' FBExpAccum accumulation ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       !--- clean up
       !---------------------------------------
    endif

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
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

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! --- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------
    call FB_getNumFlds(is_local%wrap%FBExpAccum(compocn), trim(subname)//"FBExpAccum(compocn)", ncnt, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (ncnt > 0) then

       !---------------------------------------
       !--- average ocn accumulator
       !---------------------------------------

       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExpAccum(compocn), &
               string=trim(subname)//' FBExpAccum(compocn) before avg ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       call FB_average(is_local%wrap%FBExpAccum(compocn), &
            is_local%wrap%FBExpAccumCnt(compocn), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call FB_diagnose(is_local%wrap%FBExp(compocn), &
               string=trim(subname)//' FBExpAccum(compocn) after avg ', rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       !--- copy to FBExp(compocn)
       !---------------------------------------

       call FB_copy(is_local%wrap%FBExp(compocn), is_local%wrap%FBExpAccum(compocn), rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       !--- zero accumulator
       !---------------------------------------

       is_local%wrap%FBExpAccumFlag(compocn) = .true.
       is_local%wrap%FBExpAccumCnt(compocn) = 0
       call FB_reset(is_local%wrap%FBExpAccum(compocn), value=czero, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_ocn_accum_avg

end module med_phases_prep_ocn_mod
