module med_phases_prep_rof_mod

  !-----------------------------------------------------------------------------
  ! Create rof export fields
  ! - accumulate import lnd fields on the land grid that are sent to rof 
  !   this will be done in med_phases_prep_rof_accum_fast
  ! - time avergage accumulated import lnd fields when necessary
  !   map the time averaged accumulated lnd fields to the rof grid
  !   merge the mapped lnd fields to create FBExp(comprof)
  !   this will be done in med_phases_prep_rof_avg
  !-----------------------------------------------------------------------------

  use ESMF                  , only : ESMF_FieldBundle
  use esmFlds               , only : ncomps, complnd, comprof, compname, mapconsf
  use med_constants_mod     , only : R8, CS
  use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag
  use shr_nuopc_methods_mod , only : chkerr => shr_nuopc_methods_chkerr
  use perf_mod              , only : t_startf, t_stopf

  implicit none
  private

  public  :: med_phases_prep_rof_accum_fast
  public  :: med_phases_prep_rof_avg

  private :: med_phases_prep_rof_irrig       

  type(ESMF_FieldBundle)      :: FBlndVolr          ! needed for lnd2rof irrigation
  type(ESMF_FieldBundle)      :: FBrofVolr          ! needed for lnd2rof irrigation
  type(ESMF_FieldBundle)      :: FBlndIrrig         ! needed for lnd2rof irrigation
  type(ESMF_FieldBundle)      :: FBrofIrrig         ! needed for lnd2rof irrigation

  character(len=*), parameter :: volr_field             = 'Flrr_volrmch'
  character(len=*), parameter :: irrig_flux_field       = 'Flrl_irrig'
  character(len=*), parameter :: irrig_normalized_field = 'Flrl_irrig_normalized'
  character(len=*), parameter :: irrig_volr0_field      = 'Flrl_irrig_volr0'

  character(*)    , parameter :: u_FILE_u = &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_phases_prep_rof_accum_fast(gcomp, rc)

    ! Carry out fast accumulation for the river (rof) component
    ! Accumulation and averaging is done on the land input to the river component on the land grid
    ! Mapping from the land to the rof grid is then done with the time averaged fields

    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_FieldBundleGet, ESMF_StateIsCreated, ESMF_StateGet
    use ESMF                  , only : ESMF_FieldBundleIsCreated
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_accum
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use med_internalstate_mod , only : InternalState

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState) :: is_local
    integer             :: i,j,n,ncnt
    integer             :: dbrc
    character(len=*), parameter :: subname='(med_phases_prep_rof_mod: med_phases_prep_rof_accum_fast)'
    !---------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(complnd,complnd))) then
       ncnt = 0
    else 
       ! The scalar field has been removed from all mediator field bundles - so check if the fieldCount is
       ! 0 and not 1 here
       call ESMF_FieldBundleGet(is_local%wrap%FBImp(complnd,complnd), fieldCount=ncnt, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

    if (ncnt == 0) then
       call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBimp(complnd), returning", &
            ESMF_LOGMSG_INFO, rc=dbrc)
    else

       !---------------------------------------
       ! Accumulate lnd input on lnd grid to send to rof
       !---------------------------------------
       call shr_nuopc_methods_FB_accum(is_local%wrap%FBImpAccum(complnd,complnd), &
            is_local%wrap%FBImp(complnd,complnd), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       is_local%wrap%FBImpAccumCnt(complnd) = is_local%wrap%FBImpAccumCnt(complnd) + 1

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImpAccum(complnd,complnd), &
               string=trim(subname)//' FBImpAccum(complnd,complnd) ', rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if

    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_rof_accum_fast

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_rof_avg(gcomp, rc)

    ! Prepare the ROF export Fields from the mediator

    use NUOPC                 , only : NUOPC_IsConnected
    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet 
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_FieldBundleGet
    use esmFlds               , only : fldListTo, fldListFr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_average
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
    use med_merge_mod         , only : med_merge_auto
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use med_internalstate_mod , only : InternalState, mastertask
    use med_constants_mod     , only : czero => med_constants_czero

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)         :: is_local
    integer                     :: i,j,n,n1,ncnt
    integer                     :: dbrc
    logical                     :: connected
    real(r8), pointer           :: dataptr(:)
    character(len=*),parameter  :: subname='(med_phases_prep_rof_mod: med_phases_prep_rof_avg)'
    !---------------------------------------

    call t_startf('MED:'//subname)
    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- Count the number of fields outside of scalar data, if zero, then return
    !---------------------------------------

    ! Note - the scalar field has been removed from all mediator field bundles - so this is why we check if the
    ! fieldCount is 0 and not 1 here

    call ESMF_FieldBundleGet(is_local%wrap%FBExp(comprof), fieldCount=ncnt, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (ncnt == 0) then

       call ESMF_LogWrite(trim(subname)//": only scalar data is present in FBexp(comprof), returning", &
            ESMF_LOGMSG_INFO, rc=dbrc)
    else

       !---------------------------------------
       !--- average import from land accumuled FB
       !---------------------------------------

       call shr_nuopc_methods_FB_average(is_local%wrap%FBImpAccum(complnd,complnd), &
                                         is_local%wrap%FBImpAccumCnt(complnd), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImpAccum(complnd,complnd), &
               string=trim(subname)//' FBImpAccum(complnd,complnd) after avg ', rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       !--- map to create FBImpAccum(complnd,comprof)
       !---------------------------------------

       ! The following assumes that only land import fields are needed to create the
       ! export fields for the river component and that ALL mappings are done with mapconsf

       if (is_local%wrap%med_coupling_active(complnd,comprof)) then

          call med_map_FB_Regrid_Norm( &
               fldListFr(complnd)%flds, complnd, comprof, &
               is_local%wrap%FBImpAccum(complnd,complnd), &
               is_local%wrap%FBImpAccum(complnd,comprof), &
               is_local%wrap%FBFrac(complnd), &
               is_local%wrap%FBNormOne(complnd,comprof,:), &
               is_local%wrap%RH(complnd,comprof,:), &
               string=trim(compname(complnd))//'2'//trim(compname(comprof)), rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          if (dbug_flag > 1) then
             call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBImpAccum(complnd,comprof), &
                  string=trim(subname)//' FBImpAccum(complnd,comprof) after avg ', rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Reset the irrig_flux_field with the map_lnd2rof_irrig calculation below if appropriate
          if ( NUOPC_IsConnected(is_local%wrap%NStateImp(complnd), fieldname=trim(irrig_flux_field))) then
             call med_phases_prep_rof_irrig( gcomp, rc=rc )
             if (chkerr(rc,__LINE__,u_FILE_u)) return
          else
             ! This will ensure that no irrig is sent from the land
             call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImpAccum(complnd,comprof), &
                  trim(irrig_flux_field), dataptr, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return
             dataptr(:) = 0._r8
          end if
       endif

       !---------------------------------------
       !--- auto merges to create FBExp(comprof)
       !---------------------------------------

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBFrac(comprof), &
               string=trim(subname)//' FBFrac(comprof) before merge ', rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if

       call med_merge_auto(trim(compname(comprof)), &
            is_local%wrap%FBExp(comprof), &
            is_local%wrap%FBFrac(comprof), &
            is_local%wrap%FBImpAccum(:,comprof), &
            fldListTo(comprof), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       if (dbug_flag > 1) then
          call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBExp(comprof), &
               string=trim(subname)//' FBexp(comprof) ', rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return
       end if

       !---------------------------------------
       !--- zero accumulator
       !---------------------------------------

       is_local%wrap%FBImpAccumCnt(complnd) = 0

       call shr_nuopc_methods_FB_reset(is_local%wrap%FBImpAccum(complnd,complnd), value=czero, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       !---------------------------------------
       !--- custom calculations
       !---------------------------------------

       !---------------------------------------
       !--- clean up
       !---------------------------------------

    endif

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_rof_avg

  !-----------------------------------------------------------------------------

  subroutine med_phases_prep_rof_irrig(gcomp, rc)

    !---------------------------------------------------------------
    ! Description
    ! Do custom mapping for the irrigation flux, from land -> rof.
    !
    ! The basic idea is that we want to pull irrigation out of ROF cells proportionally to
    ! the river volume (volr) in each cell. This is important in cases where the various
    ! ROF cells overlapping a CTSM cell have very different volr: If we didn't do this
    ! volr-normalized remapping, we'd try to extract the same amount of water from each
    ! of the ROF cells, which would be more likely to have withdrawals exceeding
    ! available volr.
    !
    ! (Both RTM and MOSART have code to handle excess withdrawals by pulling the excess
    ! directly out of the ocean. We'd like to avoid resorting to this if possible.
    !
    ! This mapping works by:
    ! (1) Normalizing the land's irrigation flux by volr
    ! (2) Mapping this volr-normalized flux to the rof grid
    ! (3) Converting the mapped, volr-normalized flux back to a normal
    !     (non-volr-normalized) flux on the rof grid.
    !---------------------------------------------------------------

    use ESMF                  , only : ESMF_GridComp, ESMF_Field, ESMF_FieldRegrid
    use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleGet, ESMF_FieldBundleIsCreated
    use ESMF                  , only : ESMF_SUCCESS, ESMF_FAILURE, ESMF_RouteHandleIsCreated
    use ESMF                  , only : ESMF_LOGMSG_INFO, ESMF_LogWrite
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_init
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_clean
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FieldRegrid
    use shr_nuopc_scalars_mod , only : flds_scalar_name
    use med_internalstate_mod , only : InternalState, mastertask
    use med_map_mod           , only : med_map_FB_Regrid_norm

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    integer                     :: r,l
    integer                     :: dbrc
    type(InternalState)         :: is_local
    real(r8), pointer           :: volr_l(:)
    real(r8), pointer           :: volr_r(:), volr_r_import(:)
    real(r8), pointer           :: irrig_normalized_l(:)
    real(r8), pointer           :: irrig_normalized_r(:)
    real(r8), pointer           :: irrig_volr0_l(:)
    real(r8), pointer           :: irrig_volr0_r(:)
    real(r8), pointer           :: irrig_flux_l(:)
    real(r8), pointer           :: irrig_flux_r(:)
    character(len=*), parameter :: subname='(med_phases_prep_rof_mod: med_phases_prep_rof_irrig)'
    !---------------------------------------------------------------

    call t_startf('MED:'//subname)

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    rc = ESMF_SUCCESS

    !---------------------------------------
    ! --- Get the internal state
    !---------------------------------------

    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(complnd,comprof,mapconsf), rc=rc)) then
       call ESMF_LogWrite(trim(subname)//": ERROR conservativing route handle not created for lnd->rof mapping", &
            ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    end if
    if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(comprof,complnd,mapconsf), rc=rc)) then
       call ESMF_LogWrite(trim(subname)//": ERROR conservativing route handle not created for rof->lnd mapping", &
            ESMF_LOGMSG_INFO, rc=rc)
       rc = ESMF_FAILURE
       return
    end if

    ! ------------------------------------------------------------------------
    ! Initialize module field bundles if not already initialized
    ! ------------------------------------------------------------------------

    if (.not. ESMF_FieldBundleIsCreated(FBlndVolr)  .and. &
        .not. ESMF_FieldBundleIsCreated(FBrofVolr)  .and. &
        .not. ESMF_FieldBundleIsCreated(FBlndIrrig) .and. &
        .not. ESMF_FieldBundleIsCreated(FBrofIrrig)) then

       call shr_nuopc_methods_FB_init(FBout=FBlndVolr, flds_scalar_name=flds_scalar_name, &
            FBgeom=is_local%wrap%FBImp(complnd,complnd), &
            fieldNameList=(/trim(volr_field)/), rc=rc)
       if (chkerr(rc,__line__,u_file_u)) return

       call shr_nuopc_methods_FB_init(FBout=FBrofVolr, flds_scalar_name=flds_scalar_name, &
            FBgeom=is_local%wrap%FBImp(comprof,comprof), &
            fieldNameList=(/trim(volr_field)/), rc=rc)
       if (chkerr(rc,__line__,u_file_u)) return

       call shr_nuopc_methods_FB_init(FBout=FBlndIrrig, flds_scalar_name=flds_scalar_name, &
            FBgeom=is_local%wrap%FBImp(complnd,complnd), &
            fieldNameList=(/trim(irrig_normalized_field), trim(irrig_volr0_field)/), rc=rc)
       if (chkerr(rc,__line__,u_file_u)) return

       call shr_nuopc_methods_FB_init(FBout=FBrofIrrig, flds_scalar_name=flds_scalar_name, &
            FBgeom=is_local%wrap%FBImp(comprof,comprof), &
            fieldNameList=(/trim(irrig_normalized_field), trim(irrig_volr0_field)/), rc=rc)
       if (chkerr(rc,__line__,u_file_u)) return
    end if

    ! ------------------------------------------------------------------------
    ! 1) Create volr_l: Adjust volr_r, and map it to the land grid
    ! ------------------------------------------------------------------------

    ! Treat any rof point with volr < 0 as if it had volr = 0. Negative volr values can
    ! arise in RTM. This fix is needed to avoid mapping negative irrigation to those
    ! cells: while conservative, this would be unphysical (it would mean that irrigation
    ! actually adds water to those cells).

    call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(comprof,comprof), &
         trim(volr_field), volr_r_import, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_FB_getFldPtr(FBrofVolr, trim(volr_field), volr_r, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do r = 1, size(volr_r)
       if (volr_r_import(r) < 0._r8) then
          volr_r(r) = 0._r8
       else
          volr_r(r) = volr_r_import(r)
       end if
    end do

    ! Map volr_r to volr_l (rof->lnd) using conservative mapping without any fractional weighting
    call shr_nuopc_methods_FB_FieldRegrid(FBrofVolr, trim(volr_field), FBlndVolr, trim(volr_field), &
         is_local%wrap%RH(comprof, complnd, mapconsf), rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get volr_l
    call shr_nuopc_methods_FB_getFldPtr(FBlndVolr, trim(volr_field), volr_l, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! ------------------------------------------------------------------------
    ! (2) Determine irrigation from land on land grid normalized by volr_l
    ! ------------------------------------------------------------------------

    ! In order to avoid possible divide by 0, as well as to handle non-sensical negative
    ! volr on the land grid, we divide the land's irrigation flux into two separate flux
    ! components:
    ! - a component where we have positive volr on the land grid (put in
    !   irrig_normalized_l, which is mapped using volr-normalization)
    ! - a component where we have zero or negative volr on the land
    !   grid (put in irrig_volr0_l, which is mapped as a standard flux).
    ! We then remap both of these components to the rof grid, and then
    ! finally add the two components to determine the total irrigation
    ! flux on the rof grid.

    ! First extract accumulated irrigation flux from land
    call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImpAccum(complnd,complnd), &
         trim(irrig_flux_field), irrig_flux_l, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Fill in values for irrig_normalized_l and irrig_volr0_l in temporary FBlndIrrig field bundle
    call shr_nuopc_methods_FB_getFldPtr(FBlndIrrig, trim(irrig_normalized_field), irrig_normalized_l, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_FB_getFldPtr(FBlndIrrig, trim(irrig_volr0_field), irrig_volr0_l, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do l = 1, size(volr_l)
       if (volr_l(l) > 0._r8) then
          irrig_normalized_l(l) = irrig_flux_l(l) / volr_l(l)
          irrig_volr0_l(l)      = 0._r8
       else
          irrig_normalized_l(l) = 0._r8
          irrig_volr0_l(l)      = irrig_flux_l(l)
       end if
    end do

    ! ------------------------------------------------------------------------
    ! (3) Map normalized irrigation from land to rof grid and
    !     convert to a total irrigation flux on the ROF grid
    ! ------------------------------------------------------------------------

    call med_map_FB_Regrid_Norm((/trim(irrig_normalized_field), trim(irrig_volr0_field)/), &
         FBlndIrrig, FBrofIrrig, &
         is_local%wrap%FBFrac(complnd), 'lfrin', &
         is_local%wrap%RH(complnd, comprof, mapconsf), &
         string='mapping normalized irrig from lnd to to rof', rc=rc)

    call shr_nuopc_methods_FB_getFldPtr(FBrofIrrig, trim(irrig_normalized_field), irrig_normalized_r, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_FB_getFldPtr(FBrofIrrig, trim(irrig_volr0_field), irrig_volr0_r, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Convert to a total irrigation flux on the ROF grid, and put this in the pre-merge FBImpAccum(complnd,comprof)
    call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImpAccum(complnd,comprof), &
         trim(irrig_flux_field), irrig_flux_r, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    do r = 1, size(irrig_flux_r)
       irrig_flux_r(r) = (irrig_normalized_r(r) * volr_r(r)) + irrig_volr0_r(r)
    end do

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_prep_rof_irrig

end module med_phases_prep_rof_mod
