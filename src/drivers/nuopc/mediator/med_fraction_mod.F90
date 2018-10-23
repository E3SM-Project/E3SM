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

  use med_constants_mod      , only : R8
  use esmFlds, only : ncomps

  implicit none
  private

  character(*)      , parameter :: u_FILE_u =  __FILE__

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
  character(len=5),parameter,dimension(2) :: fraclist_r = (/'lfrac','rfrac'/)
  character(len=5),parameter,dimension(1) :: fraclist_w = (/'wfrac'/)

  !--- standard ---
  real(R8),parameter :: eps_fracsum = 1.0e-02      ! allowed error in sum of fracs
  real(R8),parameter :: eps_fracval = 1.0e-02      ! allowed error in any frac +- 0,1
  real(R8),parameter :: eps_fraclim = 1.0e-03      ! truncation limit in fractions_a(lfrac)
  logical           ,parameter :: atm_frac_correct = .false. ! turn on frac correction on atm grid

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

!-----------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------

  subroutine med_fraction_init(gcomp, rc)

    ! Initialize the fractions

    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time, ESMF_State, ESMF_Field
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_GridCompGet, ESMF_StateIsCreated, ESMF_RouteHandleIsCreated
    use med_constants_mod     , only : czero=>med_constants_czero
    use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag
    use esmFlds               , only : compatm, compocn, compice, complnd
    use esmFlds               , only : comprof, compglc, compwav, compname
    use shr_nuopc_scalars_mod , only : flds_scalar_name
    use shr_nuopc_fldList_mod , only : mapconsf, mapfcopy
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_init
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_reset
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FieldRegrid
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_fldChk
    use med_map_mod           , only : med_map_Fractions_init
    use med_internalstate_mod , only : InternalState
    use perf_mod              , only : t_startf, t_stopf
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)           :: clock
    type(ESMF_Time)            :: time
    character(len=64)          :: timestr
    type(ESMF_State)           :: importState, exportState
    type(ESMF_Field)           :: field
    type(InternalState)        :: is_local
    real(R8), pointer          :: dataPtr(:)
    real(R8), pointer          :: dataPtr1(:),dataPtr2(:),dataPtr3(:),dataPtr4(:)
    real(R8), pointer          :: dataPtr_lfrac(:)
    real(R8), pointer          :: dataPtr_lfrin(:)
    real(R8), pointer          :: dataPtr_ofrac(:)
    integer                    :: i,j,n,n1
    logical, save              :: first_call = .true.
    integer                    :: maptype
    integer                    :: dbrc
    character(len=*),parameter :: subname='(med_fraction_init)'
    !---------------------------------------
    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

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
       !--- Initialize mediator FBfrac entries
       !---------------------------------------

       ! Note - must use import state here - since export state might not
       ! contain anything other than scalar data if the component is not prognostic
       do n1 = 1,ncomps
          if (is_local%wrap%comp_present(n1) .and. &
               ESMF_StateIsCreated(is_local%wrap%NStateImp(n1),rc=rc) .and. &
               ESMF_StateIsCreated(is_local%wrap%NStateExp(n1),rc=rc)) then

             call shr_nuopc_methods_FB_init(is_local%wrap%FBfrac(n1), flds_scalar_name, &
                  STgeom=is_local%wrap%NStateImp(n1), fieldNameList=fraclist(:,n1), &
                  name='FBfrac'//trim(compname(n1)), rc=rc)

             ! zero out FBfracs
             call shr_nuopc_methods_FB_reset(is_local%wrap%FBfrac(n1), value=czero, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       end do
       first_call = .false.
    endif

    !---------------------------------------
    !--- Initialize fractions on atm grid/decomp
    !---------------------------------------

    if (is_local%wrap%comp_present(compatm)) then

      ! Set atm 'afrac' to 1.
      call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'afrac', dataPtr, rc=rc)
      if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
      dataPtr(:) = 1.0_R8

      ! map atm 'afrac' to ocn 'afrac' conservatively or redist
      if (is_local%wrap%med_coupling_active(compatm,compocn)) then
          if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,compocn,mapfcopy), rc=rc)) then
             maptype = mapfcopy
          else
             maptype = mapconsf
             if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,compocn,mapconsf), rc=rc)) then
                call med_map_Fractions_init( gcomp, compatm, compocn, &
                     FBSrc=is_local%wrap%FBImp(compatm,compatm), &
                     FBDst=is_local%wrap%FBImp(compatm,compocn), &
                     RouteHandle=is_local%wrap%RH(compatm,compocn,mapconsf), rc=rc)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end if
          call shr_nuopc_methods_FB_FieldRegrid(&
               is_local%wrap%FBfrac(compatm), 'afrac', &
               is_local%wrap%FBfrac(compocn), 'afrac', &
               is_local%wrap%RH(compatm,compocn,maptype), rc=rc)
          if(shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

       ! map atm 'afrac' to ice 'afrac' conservatively or redist
       if (is_local%wrap%med_coupling_active(compatm,compice)) then
          if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,compice,mapfcopy), rc=rc)) then
             maptype = mapfcopy
          else
             maptype = mapconsf
             if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,compice,mapconsf), rc=rc)) then
                call med_map_Fractions_init( gcomp, compatm, compocn, &
                     FBSrc=is_local%wrap%FBImp(compatm,compatm), &
                     FBDst=is_local%wrap%FBImp(compatm,compice), &
                     RouteHandle=is_local%wrap%RH(compatm,compice,mapconsf), rc=rc)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
          end if
          call shr_nuopc_methods_FB_FieldRegrid(&
               is_local%wrap%FBfrac(compatm), 'afrac', &
               is_local%wrap%FBfrac(compice), 'afrac', &
               is_local%wrap%RH(compatm,compice,maptype), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

    endif

    !---------------------------------------
    !--- Initialize fractions on glc grid decomp
    !---------------------------------------

    if (is_local%wrap%comp_present(compglc)) then
       ! If 'gfrac' and 'frac' exists, then copy 'frac' to 'gfrac'
       ! TODO: implement a more general scheme that hard-wiring the name 'frac'
       if (shr_nuopc_methods_FB_FldChk(is_local%wrap%FBfrac(compglc), 'gfrac', rc=rc) .and. &
            shr_nuopc_methods_FB_FldChk(is_local%wrap%FBImp(compglc, compglc), 'frac', rc=rc)) then
          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compglc), 'gfrac', dataPtr1, rc=rc)
          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(compglc,compglc), 'frac' , dataPtr2, rc=rc)
          dataPtr1 = dataPtr2
       endif
    endif

    !---------------------------------------
    !--- Initialize fractions on land grid decomp, just an initial "guess", updated later
    !---------------------------------------

    if (is_local%wrap%comp_present(complnd)) then

       ! Set 'lfrin' (copy FBImp 'Sl_lfrin' to FBFrac 'lfrin')
       ! TODO: implement a more general scheme that hard-wiring the name 'Sl_lfrin'
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(complnd,complnd) , 'Sl_lfrin' , dataPtr2, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(complnd), 'lfrin', dataPtr1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       dataPtr1(:) = dataPtr2(:)

       if (is_local%wrap%comp_present(compatm)) then

          ! map atm 'afrac' to lnd 'afrac' conservatively or redist
          if (is_local%wrap%med_coupling_active(compatm,complnd)) then
             if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,complnd,mapfcopy), rc=rc)) then
                maptype = mapfcopy
             else
                maptype = mapconsf
                if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,complnd,mapconsf), rc=rc)) then
                   call med_map_Fractions_init( gcomp, compatm, complnd, &
                        FBSrc=is_local%wrap%FBImp(compatm,compatm), &
                        FBDst=is_local%wrap%FBImp(compatm,complnd), &
                        RouteHandle=is_local%wrap%RH(compatm,complnd,mapconsf), rc=rc)
                   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                end if
             end if
             call shr_nuopc_methods_FB_FieldRegrid(&
                  is_local%wrap%FBfrac(compatm), 'afrac', &
                  is_local%wrap%FBfrac(complnd), 'afrac', &
                  is_local%wrap%RH(compatm,complnd,maptype), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! map lnd 'lfrin' to atm 'lfrin' conservatively or redist
          if (is_local%wrap%med_coupling_active(complnd,compatm)) then
             if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,complnd,mapfcopy), rc=rc)) then
                maptype = mapfcopy
             else
                maptype = mapconsf
                if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(complnd,compatm,maptype), rc=rc)) then
                   call med_map_Fractions_init( gcomp, complnd, compatm, &
                        FBSrc=is_local%wrap%FBImp(complnd,complnd), &
                        FBDst=is_local%wrap%FBImp(complnd,compatm), &
                        RouteHandle=is_local%wrap%RH(complnd,compatm,maptype), rc=rc)
                   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                end if
             end if
             call shr_nuopc_methods_FB_FieldRegrid(&
                  is_local%wrap%FBfrac(complnd), 'lfrin', &
                  is_local%wrap%FBfrac(compatm), 'lfrin', &
                  is_local%wrap%RH(complnd,compatm,maptype), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

       end if
    end if

    !---------------------------------------
    !--- Initialize fractions on rof grid/decomp
    !---------------------------------------

    if (is_local%wrap%comp_present(comprof)) then

       ! Set 'frac' in FBfrac(comprof) to 1.
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(comprof), 'rfrac', dataPtr1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       dataPtr1(:) = 1.0_R8

       ! TODO: should this be uncommented?
       ! call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(comprof,comprof) , 'frac' , dataPtr2, rc=rc)
       ! if (.not. shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) then
       !   dataPtr1 = dataPtr2
       ! endif

    endif

    !---------------------------------------
    !--- Initialize fractions on wav grid decomp - FBFrac(compwav)
    !---------------------------------------

    if (is_local%wrap%comp_present(compwav)) then

       ! Set 'wfrac' in FBfrac(compwav) to 1.
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compwav), 'wfrac', dataPtr, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       dataPtr(:) = 1.0_R8

    endif

    !---------------------------------------
    !--- Initialize fractions on ice grid/decomp - FBFrac(compice)
    !--- Reset FBFrac(compatm) 'ofrac'
    !--- Reset FBFrac(compice) 'afrac'
    !---------------------------------------

    if (is_local%wrap%comp_present(compice)) then

       ! copy ice FBImp 'Si_imask' to FBFrac 'ofrac'
       ! set ofrac = Si_imask in FBFrac(compice)
       ! TODO: implement a more general scheme that hard-wiring the name 'frac'
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(compice,compice) , 'Si_imask' , dataPtr2, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ofrac', dataPtr1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       dataPtr1(:) = dataPtr2(:)

       if (is_local%wrap%comp_present(compatm)) then

          ! Reset FBFrac(compatm) 'ofrac' by mapping ice 'ofrac' to atm 'ofrac' conservatively or redist
          if (is_local%wrap%med_coupling_active(compice,compatm)) then
             if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compice,compatm,mapfcopy), rc=rc)) then
                maptype = mapfcopy
             else
                maptype = mapconsf
                if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compice,compatm,mapconsf), rc=rc)) then
                   call med_map_Fractions_init( gcomp, compice, compatm, &
                        FBSrc=is_local%wrap%FBImp(compice,compice), &
                        FBDst=is_local%wrap%FBImp(compice,compatm), &
                        RouteHandle=is_local%wrap%RH(compice,compatm,mapconsf), rc=rc)
                   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                end if
             end if
             call shr_nuopc_methods_FB_FieldRegrid(&
                  is_local%wrap%FBfrac(compice), 'ofrac', &
                  is_local%wrap%FBfrac(compatm), 'ofrac', &
                  is_local%wrap%RH(compice,compatm,maptype), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Reset FBFrac(compice) 'afrac' by mapping atm 'afrac' to ice 'afrac' conservatively or redist
          if (is_local%wrap%med_coupling_active(compatm,compice)) then
             if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compice,compatm,mapfcopy), rc=rc)) then
                maptype = mapfcopy
             else
                maptype = mapconsf
                if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,compice,mapconsf), rc=rc)) then
                   call med_map_Fractions_init( gcomp, compatm, compice, &
                        FBSrc=is_local%wrap%FBImp(compatm,compatm), &
                        FBDst=is_local%wrap%FBImp(compatm,compice), &
                        RouteHandle=is_local%wrap%RH(compatm,compice,mapconsf), rc=rc)
                   if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
                end if
             end if
             call shr_nuopc_methods_FB_FieldRegrid(&
                  is_local%wrap%FBfrac(compatm), 'afrac', &
                  is_local%wrap%FBfrac(compice), 'afrac', &
                  is_local%wrap%RH(compatm,compice,maptype), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

       endif
    endif

    !---------------------------------------
    !--- Initialize fractions on ocean grid/decomp
    !--- These are initialized the same as for ice
    !---------------------------------------

    if (is_local%wrap%med_coupling_active(compice,compocn)) then
       if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compice,compocn,mapfcopy), rc=rc)) then
          call med_map_Fractions_init( gcomp, compice, compocn, &
               FBSrc=is_local%wrap%FBImp(compice,compice), &
               FBDst=is_local%wrap%FBImp(compice,compocn), &
               RouteHandle=is_local%wrap%RH(compice,compocn,mapfcopy), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compocn,compice,mapfcopy), rc=rc)) then
          call med_map_Fractions_init( gcomp, compocn, compice, &
               FBSrc=is_local%wrap%FBImp(compocn,compocn), &
               FBDst=is_local%wrap%FBImp(compocn,compice), &
               RouteHandle=is_local%wrap%RH(compocn,compice,mapfcopy), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    ! map ocn 'ofrac' to atm 'ofrac' conservatively
    if (is_local%wrap%med_coupling_active(compocn,compatm)) then
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(compocn,compocn) , 'So_omask' , dataPtr2, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compocn), 'ofrac', dataPtr1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       ! Copy 'So_omask' to 'ofrac'
       dataPtr1(:) = dataPtr2(:)

       if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compocn,compatm,mapconsf), rc=rc)) then
          call med_map_Fractions_init( gcomp, compocn, compatm, &
               FBSrc=is_local%wrap%FBImp(compocn,compocn), &
               FBDst=is_local%wrap%FBImp(compocn,compatm), &
               RouteHandle=is_local%wrap%RH(compocn,compatm,mapconsf), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       call shr_nuopc_methods_FB_FieldRegrid(&
            is_local%wrap%FBfrac(compocn), 'ofrac', &
            is_local%wrap%FBfrac(compatm), 'ofrac', &
            is_local%wrap%RH(compocn,compatm,mapconsf), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

    ! map atm 'afrac' to ocn 'afrac' conservatively
    if (is_local%wrap%med_coupling_active(compatm,compocn)) then
       if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,compocn,mapconsf), rc=rc)) then
          call med_map_Fractions_init( gcomp, compatm, compocn, &
               FBSrc=is_local%wrap%FBImp(compatm,compatm), &
               FBDst=is_local%wrap%FBImp(compatm,compocn), &
               RouteHandle=is_local%wrap%RH(compatm,compocn,mapconsf), rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
       call shr_nuopc_methods_FB_FieldRegrid(&
            is_local%wrap%FBfrac(compatm), 'afrac', &
            is_local%wrap%FBfrac(compocn), 'afrac', &
            is_local%wrap%RH(compatm,compocn,mapconsf), rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    endif

    !---------------------------------------
    !--- Set ofrac and lfrac on atm grid.  These should actually be mapo2a of
    !--- ofrac and lfrac but we can't map lfrac from o2a due to masked mapping
    !--- weights.  So we have to settle for a residual calculation that is
    !--- truncated to zero to try to preserve "all ocean" cells.
    !---------------------------------------

    if (is_local%wrap%comp_present(compatm)) then

       if (is_local%wrap%comp_present(compocn) .or. is_local%wrap%comp_present(compice)) then
          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrac', dataPtr_lfrac, rc=rc)
          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ofrac', dataPtr_ofrac, rc=rc)
          do n = 1,size(dataPtr_lfrac)
             dataPtr_lfrac(n) = 1.0_R8 - dataPtr_ofrac(n)
             if (abs(dataPtr_lfrac(n)) < eps_fraclim) then
                dataPtr_lfrac(n) = 0.0_R8
                if (atm_frac_correct) then
                   dataPtr_ofrac(n) = 1.0_R8
                end if
             end if
          end do
       else if (is_local%wrap%comp_present(complnd)) then
          ! If the atmosphere is absent, then simply set lfrac=lfrin on atm grid
          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrac', dataPtr_lfrac, rc=rc)
          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrin', dataPtr_lfrin, rc=rc)
          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ofrac', dataPtr_ofrac, rc=rc)
          do n = 1,size(dataPtr_lfrac)
             dataPtr_lfrac(n) = dataPtr_lfrin(n)
             dataPtr_ofrac(n) = 1.0_R8 - dataPtr_lfrac(n)
             if (abs(dataPtr_ofrac(n)) < eps_fraclim) then
                dataPtr_ofrac(n) = 0.0_R8
                if (atm_frac_correct) then
                   dataPtr_lfrac(n) = 1.0_R8
                endif
             end if
          end do
       end if

    end if

    !---------------------------------------
    !--- finally:
    !--- set fractions_l(lfrac) from fractions_a(lfrac)
    !--- set fractions_r(lfrac) from fractions_l(lfrac)
    !--- set fractions_g(lfrac) from fractions_l(lfrac)
    !---------------------------------------

    if (is_local%wrap%comp_present(complnd)) then

       ! set fractions_l(lfrac) from fractions_a(lfrac)
       if (is_local%wrap%comp_present(compatm)) then
          if (is_local%wrap%med_coupling_active(compatm,complnd)) then
             if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,complnd,mapconsf), rc=rc)) then
                call med_map_Fractions_init( gcomp, compatm, complnd, &
                     FBSrc=is_local%wrap%FBImp(compatm,compatm), &
                     FBDst=is_local%wrap%FBImp(compatm,complnd), &
                     RouteHandle=is_local%wrap%RH(compatm,complnd,mapconsf), rc=rc)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
             call shr_nuopc_methods_FB_FieldRegrid(&
                  is_local%wrap%FBfrac(compatm), 'lfrac', &
                  is_local%wrap%FBfrac(complnd), 'lfrac', &
                  is_local%wrap%RH(compatm,complnd,mapconsf), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end if
       else
          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(complnd), 'lfrin', dataPtr_lfrin, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(complnd), 'lfrac', dataPtr_lfrac, rc=rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          dataPtr_lfrac(:) = dataPtr_lfrin(:)
       endif

       ! set fractions_r(lfrac) from fractions_l(lfrac)
       if (is_local%wrap%comp_present(comprof)) then
          if (is_local%wrap%med_coupling_active(complnd,comprof)) then
             if (.not. ESMF_RouteHandleIsCreated(is_local%wrap%RH(complnd,comprof,mapconsf), rc=rc)) then
                call med_map_Fractions_init( gcomp, complnd, comprof, &
                     FBSrc=is_local%wrap%FBImp(complnd,complnd), &
                     FBDst=is_local%wrap%FBImp(complnd,comprof), &
                     RouteHandle=is_local%wrap%RH(complnd,comprof,mapconsf), rc=rc)
                if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             end if
             call shr_nuopc_methods_FB_FieldRegrid(&
                  is_local%wrap%FBfrac(complnd), 'lfrac', &
                  is_local%wrap%FBfrac(comprof), 'lfrac', &
                  is_local%wrap%RH(complnd,comprof,mapconsf), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       end if

       ! set fractions_g(lfrac) from fractions_l(lfrac)
       if (is_local%wrap%comp_present(compglc)) then
          ! TODO: l2g_consf does not exist yet
          ! call shr_nuopc_methods_FB_FieldRegrid(&
          !      is_local%wrap%FBfrac(complnd), 'lfrac', &
          !      is_local%wrap%FBfrac(compglc), 'lfrac', &
          !      is_local%wrap%RH(complnd,compglc,mapconsf), rc=rc)
          ! if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       endif

    endif

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_fraction_init

  !-----------------------------------------------------------------------------

  subroutine med_fraction_set(gcomp, rc)

    ! Update time varying fractions

    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time, ESMF_State, ESMF_Field
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use ESMF                  , only : ESMF_GridCompGet, ESMF_FieldBundleIsCreated
    use esmFlds               , only : compatm, compocn, compice, complnd
    use esmFlds               , only : comprof, compglc, compwav, compname
    use med_constants_mod     , only : dbug_flag=>med_constants_dbug_flag
    use med_internalstate_mod , only : InternalState
    use shr_nuopc_fldList_mod , only : mapconsf, mapfcopy
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FieldRegrid
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use perf_mod              , only : t_startf, t_stopf
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)           :: clock
    type(ESMF_Time)            :: time
    character(len=64)          :: timestr
    type(ESMF_State)           :: importState, exportState
    type(ESMF_Field)           :: field
    type(InternalState)        :: is_local
    real(R8), pointer          :: dataPtr1(:),dataPtr2(:),dataPtr3(:),dataPtr4(:)
    integer                    :: i,j,n,n1
    integer                    :: dbrc
    character(len=*),parameter :: subname='(med_fraction_set)'
    !---------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
         exportState=exportState, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !---------------------------------------
    !--- update ice fraction
    !---------------------------------------

    if (is_local%wrap%comp_present(compice)) then

       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(compice,compice), 'Si_ifrac', dataPtr1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ifrac', dataPtr1, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compice), 'ofrac', dataPtr4, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! Note Si_imask is the ice domain real fraction which is a constant over time
       ! and  Si_ifrac is the time evolving ice fraction on the ice grid
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(compice,compice) , 'Si_ifrac', dataPtr3, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBImp(compice,compice) , 'Si_imask' , dataPtr2, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       ! for FBfrac(compice): set ifrac = Si_ifrac * Si_imask
       dataPtr1(:) = dataptr3(:) * dataPtr2(:)

       ! for FBfrac(compice): set ofrac = Si_imask - ifrac
       dataPtr4(:) = dataPtr2(:) - dataPtr1(:)

       ! Set ocean grid fractions
       if (is_local%wrap%comp_present(compocn)) then
          ! Map 'ifrac' from FBfrac(compice) to FBfrac(compocn)
          if (is_local%wrap%med_coupling_active(compice,compocn)) then
             call shr_nuopc_methods_FB_FieldRegrid(&
                  is_local%wrap%FBfrac(compice), 'ifrac', &
                  is_local%wrap%FBfrac(compocn), 'ifrac', &
                  is_local%wrap%RH(compice,compocn,mapfcopy), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Map 'ofrac' from FBfrac(compice) to FBfrac(comp)
          if (is_local%wrap%med_coupling_active(compice,compocn)) then
             call shr_nuopc_methods_FB_FieldRegrid(&
                  is_local%wrap%FBfrac(compice), 'ofrac', &
                  is_local%wrap%FBfrac(compocn), 'ofrac', &
                  is_local%wrap%RH(compice,compocn,mapfcopy), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       end if

       ! Set atm grid fractions for ice and ocean
       if (is_local%wrap%comp_present(compatm)) then

          ! Map 'ifrac' from FBfrac(compice) to FBfrac(compatm)
          if (is_local%wrap%med_coupling_active(compice,compatm)) then
             call shr_nuopc_methods_FB_FieldRegrid(&
                  is_local%wrap%FBfrac(compice), 'ifrac', &
                  is_local%wrap%FBfrac(compatm), 'ifrac', &
                  is_local%wrap%RH(compice,compatm,mapconsf), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Map 'ofrac' from FBfrac(compice) to FBfrac(compatm)
          if (is_local%wrap%med_coupling_active(compocn,compatm)) then
             call shr_nuopc_methods_FB_FieldRegrid(&
                  is_local%wrap%FBfrac(compice), 'ofrac', &
                  is_local%wrap%FBfrac(compatm), 'ofrac', &
                  is_local%wrap%RH(compice,compatm,mapconsf), rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          end if

          ! Note: 'lfrac' from FBFrac(compatm) is just going to be in the init
          if (atm_frac_correct) then
             call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ifrac', dataPtr1, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'ofrac', dataPtr2, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

             call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compatm), 'lfrac', dataPtr3, rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
             where (dataPtr1 + dataPtr2 > 0.0_R8)
                dataPtr1 = dataPtr1 * ((1.0_R8 - dataPtr3)/(dataPtr2+dataPtr1))
                dataPtr2 = dataPtr2 * ((1.0_R8 - dataPtr3)/(dataPtr2+dataPtr1))
             elsewhere
                dataPtr1 = 0.0_R8
                dataPtr2 = 0.0_R8
             end where
          endif
       endif

    end if

    !---------------------------------------
    !--- clean up
    !---------------------------------------

    if (dbug_flag > 5) then
       do n1 = 1,ncomps
          if (ESMF_FieldBundleIsCreated(is_local%wrap%FBfrac(n1),rc=rc)) then
             call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBfrac(n1), subname // trim(compname(n1))//' frac', rc=rc)
             if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
          endif
       enddo
       call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_fraction_set

end module med_fraction_mod
