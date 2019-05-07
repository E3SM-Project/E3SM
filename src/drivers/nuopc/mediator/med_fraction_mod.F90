module med_fraction_mod

  !-----------------------------------------------------------------------------
  ! Mediator Component.
  ! Sets fracations on all component grids
  !  the fractions fields are now afrac, ifrac, ofrac, lfrac, and lfrin.
  !    afrac = fraction of atm on a grid
  !    lfrac = fraction of lnd on a grid
  !    ifrac = fraction of ice on a grid
  !    ofrac = fraction of ocn on a grid
  !    lfrin = land fraction defined by the land model
  !    ifrad = fraction of ocn on a grid at last radiation time
  !    ofrad = fraction of ice on a grid at last radiation time
  !
  !    afrac, lfrac, ifrac, and ofrac:
  !       are the self-consistent values in the system
  !    lfrin:
  !       is the fraction on the land grid and is allowed to
  !       vary from the self-consistent value as descibed below.
  !    ifrad and ofrad:
  !       are needed for the swnet calculation.
  !
  !  the fractions fields are defined for each grid in the fraction bundles as
  !    needed as follows.
  !    character(*),parameter :: fraclist_a = 'afrac:ifrac:ofrac:lfrac:lfrin'
  !    character(*),parameter :: fraclist_o = 'afrac:ifrac:ofrac:ifrad:ofrad'
  !    character(*),parameter :: fraclist_i = 'afrac:ifrac:ofrac'
  !    character(*),parameter :: fraclist_l = 'afrac:lfrac:lfrin'
  !    character(*),parameter :: fraclist_g = 'gfrac:lfrac'
  !    character(*),parameter :: fraclist_r = 'lfrac:rfrac'
  !
  !  we assume ocean and ice are on the same grids, same masks
  !  we assume ocn2atm and ice2atm are masked maps
  !  we assume lnd2atm is a global map
  !  we assume that the ice fraction evolves in time but that
  !    the land model fraction does not.  the ocean fraction then
  !    is just the complement of the ice fraction over the region
  !    of the ocean/ice mask.
  !  we assume that component fractions sent at runtime
  !    are always the relative fraction covered.
  !    for example, if an ice cell can be up to 50% covered in
  !    ice and 50% land, then the ice domain should have a fraction
  !    value of 0.5 at that grid cell.  at run time though, the ice
  !    fraction will be between 0.0 and 1.0 meaning that grid cells
  !    is covered with between 0.0 and 0.5 by ice.  the "relative" fractions
  !    sent at run-time are corrected by the model to be total fractions
  !    such that in general, on every grid,
  !       fractions_*(afrac) = 1.0
  !       fractions_*(ifrac) + fractions_*(ofrac) + fractions_*(lfrac) = 1.0
  !  where fractions_* are a bundle of fractions on a particular grid and
  !    *frac (ie afrac) is the fraction of a particular component in the bundle.
  !
  !  the fractions are computed fundamentally as follows (although the
  !    detailed implementation might be slightly different)
  !
  !  initialization:
  !    afrac is set on all grids
  !      fractions_a(afrac) = 1.0
  !      fractions_o(afrac) = mapa2o(fractions_a(afrac))
  !      fractions_i(afrac) = mapa2i(fractions_a(afrac))
  !      fractions_l(afrac) = mapa2l(fractions_a(afrac))
  !    initially assume ifrac on all grids is zero
  !      fractions_*(ifrac) = 0.0
  !    fractions/masks provided by surface components
  !      fractions_o(ofrac) = ocean "mask" provided by ocean
  !      fractions_l(lfrin) = land "fraction  provided by land
  !    then mapped to the atm model
  !      fractions_a(ofrac) = mapo2a(fractions_o(ofrac))
  !      fractions_a(lfrin) = mapl2a(fractions_l(lfrin))
  !    and a few things are then derived
  !      fractions_a(lfrac) = 1.0 - fractions_a(ofrac)
  !           this is truncated to zero for very small values (< 0.001)
  !           to attempt to preserve non-land gridcells.
  !      fractions_l(lfrac) = mapa2l(fractions_a(lfrac))
  !      fractions_r(lfrac) = mapl2r(fractions_l(lfrac))
  !      fractions_g(lfrac) = mapl2g(fractions_l(lfrac))
  !
  !  run-time (frac_set):
  !    update fractions on ice grid
  !      fractions_i(ifrac) = i2x_i(Si_ifrac)  ! ice frac from ice model
  !      fractions_i(ofrac) = 1.0 - fractions_i(ifrac)
  !        note: the relative fractions are corrected to total fractions
  !      fractions_o(ifrac) = mapi2o(fractions_i(ifrac))
  !      fractions_o(ofrac) = mapi2o(fractions_i(ofrac))
  !      fractions_a(ifrac) = mapi2a(fractions_i(ifrac))
  !      fractions_a(ofrac) = mapi2a(fractions_i(ofrac))
  !
  !  fractions used in merging are as follows
  !  merge to atm   uses fractions_a(lfrac,ofrac,ifrac)
  !  merge to ocean uses fractions_o(ofrac,ifrac) normalized to one
  !
  !  fraction corrections in mapping are as follows
  !    mapo2a uses *fractions_o(ofrac) and /fractions_a(ofrac)
  !    mapi2a uses *fractions_i(ifrac) and /fractions_a(ifrac)
  !    mapl2a uses *fractions_l(lfrin) and /fractions_a(lfrin)
  !    mapl2g weights by fractions_l(lfrac) with normalization and multiplies by fractions_g(lfrac)
  !    mapa2* should use *fractions_a(afrac) and /fractions_*(afrac) but this
  !      has been defered since the ratio always close to 1.0
  !
  !  run time:
  !      fractions_a(lfrac) + fractions_a(ofrac) + fractions_a(ifrac) ~ 1.0
  !      0.0-eps < fractions_*(*) < 1.0+eps
  !
  ! Note that the following FBImp field names are current hard-wired below
  ! TODO: this needs to be generalized - these names should be set dynamically at run time in the
  ! source component
  !    is_local%wrap%FBImp(compglc,compglc) => 'frac'
  !    is_local%wrap%FBImp(complnd,complnd) => 'Sl_lfrin'
  !    is_local%wrap%FBImp(compice,compice) => 'Si_imask'
  !    is_local%wrap%FBImp(compocn,compocn) => 'So_omask'
  !    is_local%wrap%FBImp(compice,compice) => 'Si_ifrac' (runtime)
  !
  !-----------------------------------------------------------------------------

  use esmFlds               , only : ncomps
  use med_constants_mod     , only : R8
  use med_constants_mod     , only : dbug_flag      => med_constants_dbug_flag
  use med_constants_mod     , only : czero          => med_constants_czero
  use shr_nuopc_utils_mod   , only : chkErr         => shr_nuopc_utils_ChkErr
  use shr_nuopc_methods_mod , only : FB_init        => shr_nuopc_methods_FB_init
  use shr_nuopc_methods_mod , only : FB_reset       => shr_nuopc_methods_FB_reset
  use shr_nuopc_methods_mod , only : FB_getFldPtr   => shr_nuopc_methods_FB_getFldPtr
  use shr_nuopc_methods_mod , only : FB_FieldRegrid => shr_nuopc_methods_FB_FieldRegrid
  use shr_nuopc_methods_mod , only : FB_diagnose    => shr_nuopc_methods_FB_diagnose
  use shr_nuopc_methods_mod , only : FB_fldChk      => shr_nuopc_methods_FB_fldChk

  implicit none
  private

  ! Note - everything is private in this module other than these routines
  public med_fraction_init
  public med_fraction_set

  integer, parameter                      :: nfracs = 5
  character(len=5)                        :: fraclist(nfracs,ncomps)
  character(len=5),parameter,dimension(5) :: fraclist_a = (/'afrac','ifrac','ofrac','lfrac','lfrin'/)
  character(len=5),parameter,dimension(5) :: fraclist_o = (/'afrac','ifrac','ofrac','ifrad','ofrad'/)
  character(len=5),parameter,dimension(3) :: fraclist_i = (/'afrac','ifrac','ofrac'/)
  character(len=5),parameter,dimension(3) :: fraclist_l = (/'afrac','lfrac','lfrin'/)
  character(len=5),parameter,dimension(2) :: fraclist_g = (/'gfrac','lfrac'/)
  character(len=5),parameter,dimension(2) :: fraclist_r = (/'rfrac','lfrac'/)
  character(len=5),parameter,dimension(1) :: fraclist_w = (/'wfrac'/)

  !--- standard ---
  real(R8),parameter :: eps_fracsum = 1.0e-02      ! allowed error in sum of fracs
  real(R8),parameter :: eps_fracval = 1.0e-02      ! allowed error in any frac +- 0,1
  real(R8),parameter :: eps_fraclim = 1.0e-03      ! truncation limit in fractions_a(lfrac)
  logical ,parameter :: atm_frac_correct = .false. ! turn on frac correction on atm grid

  !--- standard plus atm fraction consistency ---
  !  real(R8),parameter :: eps_fracsum = 1.0e-12   ! allowed error in sum of fracs
  !  real(R8),parameter :: eps_fracval = 1.0e-02   ! allowed error in any frac +- 0,1
  !  real(R8),parameter :: eps_fraclim = 1.0e-03   ! truncation limit in fractions_a(lfrac)
  !  logical ,parameter :: atm_frac_correct = .true. ! turn on frac correction on atm grid

  !--- unconstrained and area conserving? ---
  !  real(R8),parameter :: eps_fracsum = 1.0e-12   ! allowed error in sum of fracs
  !  real(R8),parameter :: eps_fracval = 1.0e-02   ! allowed error in any frac +- 0,1
  !  real(R8),parameter :: eps_fraclim = 1.0e-20   ! truncation limit in fractions_a(lfrac)
  !  logical ,parameter :: atm_frac_correct = .true. ! turn on frac correction on atm grid

  character(*), parameter :: u_FILE_u =  &
       __FILE__

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

  subroutine med_fraction_init(gcomp, rc)

    ! Initialize FBFrac(:) field bundles

    use ESMF                  , only : ESMF_GridComp, ESMF_Field
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_GridCompGet, ESMF_StateIsCreated, ESMF_RouteHandleIsCreated
    use ESMF                  , only : ESMF_FieldBundle, ESMF_FieldBundleIsCreated, ESMF_FieldBundleDestroy
    use esmFlds               , only : compatm, compocn, compice, complnd
    use esmFlds               , only : comprof, compglc, compwav, compname
    use esmFlds               , only : mapconsf, mapfcopy
    use med_map_mod           , only : med_map_Fractions_init
    use med_internalstate_mod , only : InternalState
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    type(ESMF_FieldBundle)     :: FBtemp
    real(R8), pointer          :: frac(:)
    real(R8), pointer          :: ofrac(:)
    real(R8), pointer          :: lfrac(:)
    real(R8), pointer          :: ifrac(:)
    real(R8), pointer          :: afrac(:)
    real(R8), pointer          :: gfrac(:)
    real(R8), pointer          :: lfrin(:)
    real(R8), pointer          :: rfrac(:)
    real(R8), pointer          :: wfrac(:)
    real(R8), pointer          :: Sl_lfrin(:)
    real(R8), pointer          :: Si_imask(:)
    real(R8), pointer          :: So_omask(:)
    integer                    :: i,j,n,n1
    integer                    :: maptype
    integer                    :: dbrc
    logical, save              :: first_call = .true.
    character(len=*),parameter :: subname='(med_fraction_init)'
    !---------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    rc = ESMF_SUCCESS

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    if (first_call) then

       !---------------------------------------
       ! Initialize the fraclist arrays
       !---------------------------------------

       fraclist(:,:) = ' '
       fraclist(1:size(fraclist_a),compatm) = fraclist_a
       fraclist(1:size(fraclist_o),compocn) = fraclist_o
       fraclist(1:size(fraclist_i),compice) = fraclist_i
       fraclist(1:size(fraclist_l),complnd) = fraclist_l
       fraclist(1:size(fraclist_r),comprof) = fraclist_r
       fraclist(1:size(fraclist_w),compwav) = fraclist_w
       fraclist(1:size(fraclist_g),compglc) = fraclist_g

       !---------------------------------------
       ! Initialize FBFrac(:) to zero
       !---------------------------------------

       ! Note - must use import state here - since export state might not
       ! contain anything other than scalar data if the component is not prognostic
       do n1 = 1,ncomps
          if (is_local%wrap%comp_present(n1) .and. ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc)) then

             call FB_init(is_local%wrap%FBfrac(n1), is_local%wrap%flds_scalar_name, &
                  STgeom=is_local%wrap%NStateImp(n1), fieldNameList=fraclist(:,n1), &
                  name='FBfrac'//trim(compname(n1)), rc=rc)

             ! zero out FBfracs
             call FB_reset(is_local%wrap%FBfrac(n1), value=czero, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end do
       first_call = .false.
    endif

    !---------------------------------------
    ! Set 'afrac' for FBFrac(compatm), FBFrac(compice), FBFrac(compocn), FBFrac(complnd)
    !---------------------------------------

    if (is_local%wrap%comp_present(compatm)) then

       ! Set 'afrac' for FBFrac(compatm) to 1
       call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'afrac', afrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       afrac(:) = 1.0_R8

       ! Set 'afrac' for FBFrac(compice), FBFrac(compocn) and FBFrac(complnd)
       do n = 1,ncomps
          if (n == compice .or. n == compocn .or. n == complnd) then
             if (is_local%wrap%med_coupling_active(compatm,n)) then
                if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,n,mapfcopy), rc=rc)) then
                   maptype = mapfcopy
                else
                   maptype = mapconsf
                   if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,n,mapconsf), rc=rc)) then
                      call med_map_Fractions_init( gcomp, compatm, n, &
                           FBSrc=is_local%wrap%FBImp(compatm,compatm), &
                           FBDst=is_local%wrap%FBImp(compatm,n), &
                           RouteHandle=is_local%wrap%RH(compatm,n,mapconsf), rc=rc)
                      if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   end if
                end if
                call FB_FieldRegrid(&
                     is_local%wrap%FBfrac(compatm), 'afrac', &
                     is_local%wrap%FBfrac(n), 'afrac', &
                     is_local%wrap%RH(compatm,n,maptype), rc=rc)
                if(ChkErr(rc,__LINE__,u_FILE_u)) return
             endif
          end if
       end do
    end if

    !---------------------------------------
    ! Set 'lfrin' for FBFrac(complnd) and FBFrac(compatm)
    !---------------------------------------

    ! The following is just an initial "guess", updated later

    if (is_local%wrap%comp_present(complnd)) then

       ! Set 'lfrin' for FBFrac(complnd)
       call FB_getFldPtr(is_local%wrap%FBImp(complnd,complnd) , 'Sl_lfrin' , Sl_lfrin, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBfrac(complnd), 'lfrin', lfrin, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       lfrin(:) = Sl_lfrin(:)

       ! Set 'lfrin for FBFrac(compatm)
       if (is_local%wrap%comp_present(compatm) .and. (is_local%wrap%med_coupling_active(compatm,complnd))) then
          ! Note - need to do the following if compatm->complnd is active, even if complnd->compatm is not active

          ! Create a temporary field bundle if one does not exists
          if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(complnd,compatm))) then
             call FB_init(FBout=FBtemp, &
                  flds_scalar_name=is_local%wrap%flds_scalar_name, &
                  FBgeom=is_local%wrap%FBImp(compatm,compatm), &
                  fieldNameList=(/'Fldtemp'/), name='FBtemp', rc=rc)
             if (chkerr(rc,__line__,u_file_u)) return
          end if

          ! Determine map type
          if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,complnd,mapfcopy), rc=rc)) then
             maptype = mapfcopy
          else
             maptype = mapconsf
          end if

          ! Create route handle from lnd->atm if necessary
          if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(complnd,compatm,maptype), rc=rc)) then
             if (ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(complnd,compatm))) then
                call med_map_Fractions_init( gcomp, complnd, compatm, &
                  FBSrc=is_local%wrap%FBImp(complnd,complnd), &
                  FBDst=is_local%wrap%FBImp(complnd,compatm), &
                  RouteHandle=is_local%wrap%RH(complnd,compatm,maptype), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             else
                call med_map_Fractions_init( gcomp, complnd, compatm, &
                     FBSrc=is_local%wrap%FBImp(complnd,complnd), &
                     FBDst=FBtemp, &
                     RouteHandle=is_local%wrap%RH(complnd,compatm,maptype), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end if

          ! Regrid 'lfrin' from FBFrac(complnd) -> FBFrac(compatm)
          call FB_FieldRegrid(&
               is_local%wrap%FBfrac(complnd), 'lfrin', &
               is_local%wrap%FBfrac(compatm), 'lfrin', &
               is_local%wrap%RH(complnd,compatm,maptype), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Destroy temporary field bundle if created
          if (ESMF_FieldBundleIsCreated(FBTemp)) then
             call ESMF_FieldBundleDestroy(FBtemp, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end if
    end if

    !---------------------------------------
    ! Set 'ifrac' in FBFrac(compice) and BFrac(compatm)
    !---------------------------------------

    if (is_local%wrap%comp_present(compice)) then

       ! Set 'ifrac' FBFrac(compice)
       call FB_getFldPtr(is_local%wrap%FBImp(compice,compice) , 'Si_imask' , Si_imask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ifrac', ifrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ifrac(:) = Si_imask(:)

       ! Set 'ifrac' in  FBFrac(compatm)
       if (is_local%wrap%comp_present(compatm)) then
          if (is_local%wrap%med_coupling_active(compice,compatm)) then
             if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compice,compatm,mapfcopy), rc=rc)) then
                maptype = mapfcopy
             else
                maptype = mapconsf
             end if
             if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compice,compatm,maptype), rc=rc)) then
                call med_map_Fractions_init( gcomp, compice, compatm, &
                     FBSrc=is_local%wrap%FBImp(compice,compice), &
                     FBDst=is_local%wrap%FBImp(compice,compatm), &
                     RouteHandle=is_local%wrap%RH(compice,compatm,maptype), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
             call FB_FieldRegrid(&
                  is_local%wrap%FBfrac(compice), 'ifrac', &
                  is_local%wrap%FBfrac(compatm), 'ifrac', &
                  is_local%wrap%RH(compice,compatm,maptype), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       endif
    endif

    !---------------------------------------
    ! Set 'ofrac' in FBFrac(compocn) and FBFrac(compatm)
    !---------------------------------------

    if (is_local%wrap%comp_present(compocn)) then

       call FB_getFldPtr(is_local%wrap%FBImp(compocn,compocn) , 'So_omask', So_omask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac', ofrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ofrac(:) = So_omask(:)

       if (is_local%wrap%med_coupling_active(compocn,compatm)) then
          if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compocn,compatm,mapconsf), rc=rc)) then
             call med_map_Fractions_init( gcomp, compocn, compatm, &
                  FBSrc=is_local%wrap%FBImp(compocn,compocn), &
                  FBDst=is_local%wrap%FBImp(compocn,compatm), &
                  RouteHandle=is_local%wrap%RH(compocn,compatm,mapconsf), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call FB_FieldRegrid(&
               is_local%wrap%FBfrac(compocn), 'ofrac', &
               is_local%wrap%FBfrac(compatm), 'ofrac', &
               is_local%wrap%RH(compocn,compatm,mapconsf), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if


    !---------------------------------------
    ! Set 'lfrac' in FBFrac(compatm) and correct 'ofrac' in FBFrac(compatm)
    ! ---------------------------------------

    ! These should actually be mapo2a of ofrac and lfrac but we can't
    ! map lfrac from o2a due to masked mapping weights.  So we have to
    ! settle for a residual calculation that is truncated to zero to
    ! try to preserve "all ocean" cells.

    if (is_local%wrap%comp_present(compatm)) then

       if (is_local%wrap%comp_present(compocn) .or. is_local%wrap%comp_present(compice)) then
          call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrac', lfrac, rc=rc)
          call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ofrac', ofrac, rc=rc)

          if (.not. is_local%wrap%comp_present(complnd)) then
             lfrac(:) = 0.0_R8
             if (atm_frac_correct) then
                ofrac(:) = 1.0_R8
             end if
          else
             do n = 1,size(lfrac)
                lfrac(n) = 1.0_R8 - ofrac(n)
                if (abs(lfrac(n)) < eps_fraclim) then
                   lfrac(n) = 0.0_R8
                   if (atm_frac_correct) then
                      ofrac(n) = 1.0_R8
                   end if
                end if
             end do
          end if

       else if (is_local%wrap%comp_present(complnd)) then

          ! If the ocean or ice are absent, then simply set 'lfrac' to 'lfrin' for FBFrac(compatm)
          call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrin', lfrin, rc=rc)
          call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrac', lfrac, rc=rc)
          call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ofrac', ofrac, rc=rc)
          do n = 1,size(lfrac)
             lfrac(n) = lfrin(n)
             ofrac(n) = 1.0_R8 - lfrac(n)
             if (abs(ofrac(n)) < eps_fraclim) then
                ofrac(n) = 0.0_R8
                if (atm_frac_correct) then
                   lfrac(n) = 1.0_R8
                endif
             end if
          end do

       end if
    end if

    !---------------------------------------
    ! Set 'lfrac' in FBFrac(complnd)
    !---------------------------------------

    if (is_local%wrap%comp_present(complnd)) then

       ! Set 'lfrac' in FBFrac(complnd)
       if (is_local%wrap%comp_present(compatm)) then
          ! If atm -> lnd coupling is active - map 'lfrac' from FBFrac(compatm) to FBFrac(complnd)
          if (is_local%wrap%med_coupling_active(compatm,complnd)) then
             if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,complnd,mapconsf), rc=rc)) then
                call med_map_Fractions_init( gcomp, compatm, complnd, &
                     FBSrc=is_local%wrap%FBImp(compatm,compatm), &
                     FBDst=is_local%wrap%FBImp(compatm,complnd), &
                     RouteHandle=is_local%wrap%RH(compatm,complnd,mapconsf), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
             call FB_FieldRegrid(&
                  is_local%wrap%FBfrac(compatm), 'lfrac', &
                  is_local%wrap%FBfrac(complnd), 'lfrac', &
                  is_local%wrap%RH(compatm,complnd,mapconsf), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       else
          ! If the atm ->lnd coupling is not active - simply set 'lfrac' to 'lfrin' 
          call FB_getFldPtr(is_local%wrap%FBfrac(complnd), 'lfrin', lfrin, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call FB_getFldPtr(is_local%wrap%FBfrac(complnd), 'lfrac', lfrac, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          lfrac(:) = lfrin(:)
       endif

    endif

    !---------------------------------------
    ! Set 'rfrac' and 'lfrac' for FBFrac(comprof)
    !---------------------------------------

    if (is_local%wrap%comp_present(comprof)) then

       ! Set 'rfrac' in FBFrac(comprof)
       if ( FB_FldChk(is_local%wrap%FBfrac(comprof)        , 'rfrac', rc=rc) .and. &
            FB_FldChk(is_local%wrap%FBImp(comprof, comprof), 'frac' , rc=rc)) then
          call FB_getFldPtr(is_local%wrap%FBfrac(comprof)       , 'rfrac', rfrac, rc=rc)
          call FB_getFldPtr(is_local%wrap%FBImp(comprof,comprof), 'frac' , frac, rc=rc)
          rfrac(:) = frac(:)
       else
          ! Set 'rfrac' in FBfrac(comprof) to 1.
          call FB_getFldPtr(is_local%wrap%FBfrac(comprof), 'rfrac', rfrac, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          rfrac(:) = 1.0_R8
       endif

       ! Set 'lfrac' in FBFrac(comprof)
       if (is_local%wrap%comp_present(complnd)) then
          if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(complnd,comprof,mapconsf), rc=rc)) then
             call med_map_Fractions_init( gcomp, complnd, comprof, &
                  FBSrc=is_local%wrap%FBImp(complnd,complnd), &
                  FBDst=is_local%wrap%FBImp(complnd,comprof), &
                  RouteHandle=is_local%wrap%RH(complnd,comprof,mapconsf), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call FB_FieldRegrid(&
               is_local%wrap%FBfrac(complnd), 'lfrac', &
               is_local%wrap%FBfrac(comprof), 'lfrac', &
               is_local%wrap%RH(complnd,comprof,mapconsf), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
    endif

    !---------------------------------------
    ! Set 'gfrac' and 'lfrac' for FBFrac(compglc)
    !---------------------------------------

    if (is_local%wrap%comp_present(compglc)) then
       ! Set 'gfrac' in FBFrac(compglc)
       if ( FB_FldChk(is_local%wrap%FBfrac(compglc)        , 'gfrac', rc=rc) .and. &
            FB_FldChk(is_local%wrap%FBImp(compglc, compglc), 'frac' , rc=rc)) then
          call FB_getFldPtr(is_local%wrap%FBfrac(compglc)       , 'gfrac', gfrac, rc=rc)
          call FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), 'frac' , frac, rc=rc)
          gfrac(:) = frac(:)
       else
          ! Set 'gfrac' in FBfrac(compglc) to 1.
          call FB_getFldPtr(is_local%wrap%FBfrac(compglc), 'gfrac', gfrac, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          gfrac(:) = 1.0_R8
       endif

       ! Set 'lfrac' in FBFrac(compglc)
       if (is_local%wrap%comp_present(complnd)) then
          if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(complnd,compglc,mapconsf), rc=rc)) then
             call med_map_Fractions_init( gcomp, complnd, compglc, &
                  FBSrc=is_local%wrap%FBImp(complnd,complnd), &
                  FBDst=is_local%wrap%FBImp(complnd,compglc), &
                  RouteHandle=is_local%wrap%RH(complnd,compglc,mapconsf), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call FB_FieldRegrid(&
               is_local%wrap%FBfrac(complnd), 'lfrac', &
               is_local%wrap%FBfrac(compglc), 'lfrac', &
               is_local%wrap%RH(complnd,compglc,mapconsf), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif
    endif

    !---------------------------------------
    ! Set 'wfrac' for FBFrac(compwav)
    !---------------------------------------

    if (is_local%wrap%comp_present(compwav)) then
       ! Set 'wfrac' in FBfrac(compwav) to 1.
       call FB_getFldPtr(is_local%wrap%FBfrac(compwav), 'wfrac', wfrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       wfrac(:) = 1.0_R8
    endif

    !---------------------------------------
    ! Diagnostic output
    !---------------------------------------

    if (dbug_flag > 1) then
       do n = 1,ncomps
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBfrac(n),rc=rc)) then
             call FB_diagnose(is_local%wrap%FBfrac(n), &
                  trim(subname) // trim(compname(n)), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end do
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_fraction_init

  !-----------------------------------------------------------------------------

  subroutine med_fraction_set(gcomp, rc)

    ! Update time varying fractions

    use ESMF                  , only : ESMF_GridComp, ESMF_GridCompGet
    use ESMF                  , only : ESMF_RouteHandleIsCreated, ESMF_FieldBundleIsCreated
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_REGION_TOTAL, ESMF_REGION_SELECT
    use esmFlds               , only : compatm, compocn, compice, compname
    use esmFlds               , only : mapconsf, mapnstod, mapfcopy
    use esmFlds               , only : coupling_mode
    use med_internalstate_mod , only : InternalState
    use med_map_mod           , only : med_map_Fractions_init
    use perf_mod              , only : t_startf, t_stopf

    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(InternalState)        :: is_local
    real(r8), pointer          :: lfrac(:)
    real(r8), pointer          :: ifrac(:)
    real(r8), pointer          :: ofrac(:)
    real(r8), pointer          :: Si_ifrac(:)
    real(r8), pointer          :: Si_imask(:)
    integer                    :: n
    integer                    :: dbrc
    integer                    :: maptype
    character(len=*),parameter :: subname='(med_fraction_set)'
    !---------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    end if
    rc = ESMF_SUCCESS

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    ! Update FBFrac(compice), FBFrac(compocn) and FBFrac(compatm) field bundles
    !---------------------------------------

    if (is_local%wrap%comp_present(compice) .and. is_local%wrap%comp_present(compocn)) then
       if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compice,compocn,mapfcopy), rc=rc)) then
          if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(compice,compocn))) then
             call FB_init(is_local%wrap%FBImp(compice,compocn), is_local%wrap%flds_scalar_name, &
                  STgeom=is_local%wrap%NStateImp(compocn), &
                  STflds=is_local%wrap%NStateImp(compice), &
                  name='FBImp'//trim(compname(compice))//'_'//trim(compname(compocn)), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call med_map_Fractions_init( gcomp, compice, compocn, &
               FBSrc=is_local%wrap%FBImp(compice,compice), &
               FBDst=is_local%wrap%FBImp(compice,compocn), &
               RouteHandle=is_local%wrap%RH(compice,compocn,mapfcopy), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compocn,compice,mapfcopy), rc=rc)) then
          if (.not. ESMF_FieldBundleIsCreated(is_local%wrap%FBImp(compocn,compice))) then
             call FB_init(is_local%wrap%FBImp(compocn,compice), is_local%wrap%flds_scalar_name, &
                  STgeom=is_local%wrap%NStateImp(compice), &
                  STflds=is_local%wrap%NStateImp(compocn), &
                  name='FBImp'//trim(compname(compocn))//'_'//trim(compname(compice)), rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
          call med_map_Fractions_init( gcomp, compocn, compice, &
               FBSrc=is_local%wrap%FBImp(compocn,compocn), &
               FBDst=is_local%wrap%FBImp(compocn,compice), &
               RouteHandle=is_local%wrap%RH(compocn,compice,mapfcopy), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    if (is_local%wrap%comp_present(compice)) then

       ! -------------------------------------------
       ! Set FBfrac(compice)
       ! -------------------------------------------

       ! Si_imask is the ice domain mask which is constant over time
       ! Si_ifrac is the time evolving ice fraction on the ice grid

       call FB_getFldPtr(is_local%wrap%FBImp(compice,compice) , 'Si_ifrac', Si_ifrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBImp(compice,compice) , 'Si_imask' , Si_imask, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       call FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ifrac', ifrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ofrac', ofrac, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return

       ! set ifrac = Si_ifrac * Si_imask
       ifrac(:) = Si_ifrac(:) * Si_imask(:)

       if (trim(coupling_mode) == 'nems_orig') then
          ofrac(:) = 1._r8 - ifrac(:)
       else
          ! set ofrac = Si_imask - ifrac
          ofrac(:) = Si_imask(:) - ifrac(:)
       end if

       ! -------------------------------------------
       ! Set FBfrac(compocn)
       ! -------------------------------------------

       ! The following is just a redistribution from FBFrac(compice)

       if (is_local%wrap%comp_present(compocn)) then
          ! Map 'ifrac' from FBfrac(compice) to FBfrac(compocn)
          call FB_FieldRegrid(&
               is_local%wrap%FBfrac(compice), 'ifrac', &
               is_local%wrap%FBfrac(compocn), 'ifrac', &
               is_local%wrap%RH(compice,compocn,mapfcopy), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return

          ! Map 'ofrac' from FBfrac(compice) to FBfrac(compocn)
          call FB_FieldRegrid(&
               is_local%wrap%FBfrac(compice), 'ofrac', &
               is_local%wrap%FBfrac(compocn), 'ofrac', &
               is_local%wrap%RH(compice,compocn,mapfcopy), rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

       ! -------------------------------------------
       ! Set FBfrac(compatm)
       ! -------------------------------------------
       if (is_local%wrap%comp_present(compatm)) then

          if (trim(coupling_mode) == 'nems_orig') then

             ! Map 'ifrac' from FBfrac(compice) to FBfrac(compatm)
             call FB_FieldRegrid(&
                  is_local%wrap%FBfrac(compice), 'ifrac', &
                  is_local%wrap%FBfrac(compatm), 'ifrac', &
                  is_local%wrap%RH(compice,compatm,mapnstod), &
                  zeroregion=ESMF_REGION_TOTAL, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             call FB_FieldRegrid(&
                  is_local%wrap%FBfrac(compice), 'ifrac', &
                  is_local%wrap%FBfrac(compatm), 'ifrac', &
                  is_local%wrap%RH(compice,compatm,mapconsf), &
                  zeroregion=ESMF_REGION_SELECT, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             ! Now set ofrac=1-ifrac and lfrac=0 on the atm grid
             call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ifrac', ifrac, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ofrac', ofrac, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
             call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrac', lfrac, rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return

             ofrac = 1.0_R8 - ifrac
             lfrac = 0.0_R8

          else

             if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compice,compatm,mapfcopy), rc=rc)) then
                maptype = mapfcopy
             else
                maptype = mapconsf
             end if

             ! Map 'ifrac' from FBfrac(compice) to FBfrac(compatm)
             if (is_local%wrap%med_coupling_active(compice,compatm)) then
                call FB_FieldRegrid(&
                     is_local%wrap%FBfrac(compice), 'ifrac', &
                     is_local%wrap%FBfrac(compatm), 'ifrac', &
                     is_local%wrap%RH(compice,compatm,maptype), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if

             ! Map 'ofrac' from FBfrac(compice) to FBfrac(compatm)
             if (is_local%wrap%med_coupling_active(compocn,compatm)) then
                call FB_FieldRegrid(&
                     is_local%wrap%FBfrac(compice), 'ofrac', &
                     is_local%wrap%FBfrac(compatm), 'ofrac', &
                     is_local%wrap%RH(compice,compatm,maptype), rc=rc)
                if (ChkErr(rc,__LINE__,u_FILE_u)) return
             end if

             ! Note: 'lfrac' from FBFrac(compatm) is just going to be in the init
             if ( is_local%wrap%med_coupling_active(compice,compatm) .and. &
                  is_local%wrap%med_coupling_active(compocn,compatm) ) then

                if (atm_frac_correct) then
                   call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ifrac', ifrac, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return

                   call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ofrac', ofrac, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return

                   call FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrac', lfrac, rc=rc)
                   if (ChkErr(rc,__LINE__,u_FILE_u)) return
                   where (ifrac + ofrac > 0.0_R8)
                      ifrac = ifrac * ((1.0_R8 - lfrac)/(ofrac+ifrac))
                      ofrac = ofrac * ((1.0_R8 - lfrac)/(ofrac+ifrac))
                   elsewhere
                      ifrac = 0.0_R8
                      ofrac = 0.0_R8
                   end where
                endif
             endif

          end if
       end if
    end if

    !---------------------------------------
    ! Diagnostic output
    !---------------------------------------

    if (dbug_flag > 1) then
       do n = 1,ncomps
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBfrac(n),rc=rc)) then
             call FB_diagnose(is_local%wrap%FBfrac(n), &
                  trim(subname) // trim(compname(n))//' frac', rc=rc)
             if (ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       enddo
    end if

    if (dbug_flag > 20) then
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    end if

    call t_stopf('MED:'//subname)

  end subroutine med_fraction_set

end module med_fraction_mod
