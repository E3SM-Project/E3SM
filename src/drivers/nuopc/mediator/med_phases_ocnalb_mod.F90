module med_phases_ocnalb_mod

  use ESMF
  use NUOPC
  use shr_kind_mod          , only : r8=>shr_kind_r8, in=>shr_kind_in
  use shr_kind_mod          , only : cs=>shr_kind_cs, cl=>shr_kind_cl
  use shr_sys_mod           , only : shr_sys_abort
  use shr_orb_mod           , only : shr_orb_cosz, shr_orb_decl
  use shr_const_mod         , only : shr_const_pi, shr_const_spval
  use esmFlds               , only : fldListMed_ocnalb_o
  use esmFlds               , only : flds_scalar_name
  use esmFlds               , only : flds_scalar_num
  use esmFlds               , only : compatm, compocn
  use esmFlds               , only : flds_scalar_index_nextsw_cday
  use shr_nuopc_fldList_mod , only : mapconsf, mapnames
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_init
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFieldN
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FieldRegrid
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_GetScalar
  use shr_nuopc_methods_mod , only : shr_nuopc_methods_ChkErr
  use seq_timemgr_mod       , only : seq_timemgr_EclockGetData
  use med_map_mod           , only : med_map_FB_Regrid_Norm 
  use med_constants_mod     , only : med_constants_dbug_flag
  use med_constants_mod     , only : med_constants_czero
  use med_internalstate_mod , only : InternalState 

  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public med_phases_ocnalb_init
  public med_phases_ocnalb_run
  public med_phases_ocnalb_mapo2a

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  real(r8) , pointer :: lats  (:) ! latitudes  (degrees)
  real(r8) , pointer :: lons  (:) ! longitudes (degrees)
  integer  , pointer :: mask  (:) ! ocn domain mask: 0 <=> inactive cell
  real(r8) , pointer :: anidr (:) ! albedo: near infrared, direct
  real(r8) , pointer :: avsdr (:) ! albedo: visible      , direct
  real(r8) , pointer :: anidf (:) ! albedo: near infrared, diffuse
  real(r8) , pointer :: avsdf (:) ! albedo: visible      , diffuse
  real(r8) , pointer :: swndr (:) ! direct near-infrared  incident solar radiation
  real(r8) , pointer :: swndf (:) ! diffuse near-infrared incident solar radiation
  real(r8) , pointer :: swvdr (:) ! direct visible  incident solar radiation
  real(r8) , pointer :: swvdf (:) ! diffuse visible incident solar radiation
  real(r8) , pointer :: swdn  (:) ! short wave, downward (only used for diurnal calc)
  real(r8) , pointer :: swup  (:) ! short wave, upward (only used for diurnal calc)

  ! Conversion from degrees to radians
  integer                :: dbug_flag = med_constants_dbug_flag
  integer                :: dbrc
  character(len=1024)    :: tmpstr
  real(ESMF_KIND_R8)    , parameter :: czero = med_constants_czero
  real(ESMF_KIND_R8)    , parameter :: const_deg2rad = shr_const_pi/180.0_ESMF_KIND_R8  ! deg to rads
  character(*),parameter :: u_FILE_u = __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine med_phases_ocnalb_init(gcomp, rc)

    !-----------------------------------------------------------------------
    ! Initialize pointers to the module variables and then use the module
    ! variables in the med_ocnalb phase
    ! All input field bundles are ASSUMED to be on the ocean grid
    !-----------------------------------------------------------------------

    ! Arguments
    type(ESMF_GridComp)    :: gcomp
    integer, intent(out)   :: rc
    !
    ! Local variables
    type(ESMF_VM)               :: vm
    integer(in)                 :: iam
    type(ESMF_Field)            :: lfield
    type(ESMF_Grid)             :: lgrid
    type(ESMF_Mesh)             :: lmesh
    type(ESMF_GeomType_Flag)    :: geomtype
    integer                     :: n
    integer                     :: lsize
    real(ESMF_KIND_R8), pointer :: rmask(:)  ! ocn domain mask
    integer                     :: dimCount
    integer                     :: spatialDim
    integer                     :: numOwnedElements
    type(InternalState)         :: is_local
    real(ESMF_KIND_R8), pointer :: ownedElemCoords(:)
    character(len=CL)           :: tempc1,tempc2
    character(*), parameter     :: subname = '(med_phases_ocnalb_init) '
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

    ! Get the internal state from gcomp
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Set pointers to fields needed for albedo calculations 
    !----------------------------------

    ! These must must be on the ocean grid since the ocean albedo computation is on the ocean grid
    ! The following sets pointers to the module arrays

    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Faxa_swndr', fldptr1=swndr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Faxa_swndf', fldptr1=swndf, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Faxa_swvdr', fldptr1=swvdr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBImp(compatm,compocn), fldname='Faxa_swvdf', fldptr1=swvdf, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, fldname='So_avsdr', fldptr1=avsdr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, fldname='So_avsdf', fldptr1=avsdf, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, fldname='So_anidr', fldptr1=anidr, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, fldname='So_anidf', fldptr1=anidf, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, fldname='Faox_swdn', fldptr1=swdn, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, fldname='Faox_swup', fldptr1=swup, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Get lat, lon, which are time-invariant
    !----------------------------------

    ! The following assumes that all fields in FBMed_ocnalb_o have the same grid - so
    ! only need to query field 1
    call shr_nuopc_methods_FB_getFieldN(is_local%wrap%FBMed_ocnalb_o, fieldnum=1, field=lfield, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Determine if first field is on a grid or a mesh - default will be mesh
    call ESMF_FieldGet(lfield, geomtype=geomtype, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (geomtype == ESMF_GEOMTYPE_MESH) then
       call ESMF_LogWrite(trim(subname)//" : FBAtm is on a mesh ", ESMF_LOGMSG_INFO, rc=rc)
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_MeshGet(lmesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       lsize = size(swndr)
       if (numOwnedElements /= lsize) then
          write(tempc1,'(i10)') numOwnedElements
          write(tempc2,'(i10)') lsize
          call ESMF_LogWrite(trim(subname)//": ERROR numOwnedElements "// trim(tempc1) // &
               " not equal to local size "// trim(tempc2), ESMF_LOGMSG_INFO, rc=rc)
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
      call ESMF_LogWrite(trim(subname)//": ERROR FBATM must be either on a grid or a mesh", ESMF_LOGMSG_INFO, rc=rc)
      rc = ESMF_FAILURE
      return
    end if

    ! Compute ocean albedoes
    ! This will update the relevant module arrays 

    call med_phases_ocnalb_run(gcomp, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif

  end subroutine med_phases_ocnalb_init

  !===============================================================================
  
  subroutine med_phases_ocnalb_run(gcomp, rc)

    ! Compute ocean albedos (on the ocean grid)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    logical                       :: first_time = .true.
    logical                       :: update_alb
    type(InternalState)           :: is_local
    type(ESMF_Clock)              :: clock
    type(ESMF_Time)               :: currTime
    character(CL)                 :: cvalue
    character(CS)                 :: starttype        ! config start type
    character(CL)                 :: runtype          ! initial, continue, hybrid, branch
    character(CL)                 :: aoflux_grid
    logical                       :: flux_albav       ! flux avg option
    real(ESMF_KIND_R8)            :: nextsw_cday      ! calendar day of next atm shortwave
    real(ESMF_KIND_R8), pointer   :: ofrac(:)
    real(ESMF_KIND_R8), pointer   :: ofrad(:)
    real(ESMF_KIND_R8), pointer   :: ifrac(:)
    real(ESMF_KIND_R8), pointer   :: ifrad(:)
    integer                       :: lsize            ! local size
    integer                       :: n,i              ! indices
    real(ESMF_KIND_R8)            :: rlat             ! gridcell latitude in radians
    real(ESMF_KIND_R8)            :: rlon             ! gridcell longitude in radians
    real(ESMF_KIND_R8)            :: cosz             ! Cosine of solar zenith angle
    real(ESMF_KIND_R8)            :: eccen            ! Earth orbit eccentricity
    real(ESMF_KIND_R8)            :: mvelpp           ! Earth orbit
    real(ESMF_KIND_R8)            :: lambm0           ! Earth orbit
    real(ESMF_KIND_R8)            :: obliqr           ! Earth orbit
    real(ESMF_KIND_R8)            :: delta            ! Solar declination angle  in radians
    real(ESMF_KIND_R8)            :: eccf             ! Earth orbit eccentricity factor
    real(ESMF_KIND_R8), parameter :: albdif = 0.06_r8 ! 60 deg reference albedo, diffuse
    real(ESMF_KIND_R8), parameter :: albdir = 0.07_r8 ! 60 deg reference albedo, direct
    character(len=*)  , parameter :: subname='(med_phases_ocnalb_run)'
    !---------------------------------------

    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Note that in the mct version the atm was initialized first so
    ! that nextsw_cday could be passed to the other components - this
    ! assumed that atmosphere component was ALWAYS initialized first.
    ! In the nuopc version it will be easier to assume that on startup
    ! - nextsw_cday is just what cam was setting it as the current
    ! calendar day

    ! TODO: need to determine how to handle restart and branch runs -
    ! for now will just assume that nextsw_cday is not computed on
    ! initialization and is read from the restart file.

    if (first_time) then

       call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       read(cvalue,*) starttype

       if (trim(starttype) == trim('startup')) then
          runtype = "initial"
       else if (trim(starttype) == trim('continue') ) then
          runtype = "continue"
       else if (trim(starttype) == trim('branch')) then
          runtype = "continue"
       else
          call shr_sys_abort( subname//' ERROR: unknown starttype' )
       end if

       if (trim(runtype) /= 'initial') then
          nextsw_cday = -1.0_ESMF_KIND_R8
       else
          call ESMF_GridCompGet(gcomp, clock=clock)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          call ESMF_ClockGet( clock, currTime=currTime, rc=rc )
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
          call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       end if
       first_time = .false.

    else

       ! Note that shr_nuopc_methods_State_GetScalar includes a broadcast to all other pets im mpicom
       call shr_nuopc_methods_State_GetScalar(is_local%wrap%NstateImp(compatm), &
            flds_scalar_name=flds_scalar_name, flds_scalar_num=flds_scalar_num, &
            scalar_id=flds_scalar_index_nextsw_cday, value=nextsw_cday, mpicom=is_local%wrap%mpicom, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    end if

    call NUOPC_CompAttributeGet(gcomp, name='flux_albav', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flux_albav

    call NUOPC_CompAttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) eccen

    call NUOPC_CompAttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) obliqr

    call NUOPC_CompAttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) lambm0

    call NUOPC_CompAttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) mvelpp

    ! Map relevant atm fluxes to the ocean grid
    ! NOTE: there are already pointers in place as module variables in med_ocnalb_mod.F90
    ! to the arrays in FBMed_aoflux_o below - so do not need to pass them as arguments to med_ocnalb_run

    if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compatm,compocn,mapconsf), rc=rc)) then
       call shr_nuopc_methods_FB_FieldRegrid(&
            is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdf', &
            is_local%wrap%FBImp(compatm,compocn), 'Faxa_swvdf', &
            is_local%wrap%RH(compatm,compocn,mapconsf), rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_FieldRegrid(&
            is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndf', &
            is_local%wrap%FBImp(compatm,compocn), 'Faxa_swndf', &
            is_local%wrap%RH(compatm,compocn,mapconsf), rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_FieldRegrid(&
            is_local%wrap%FBImp(compatm,compatm), 'Faxa_swvdr', &
            is_local%wrap%FBImp(compatm,compocn), 'Faxa_swvdr', &
            is_local%wrap%RH(compatm,compocn,mapconsf), rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

       call shr_nuopc_methods_FB_FieldRegrid(&
            is_local%wrap%FBImp(compatm,compatm), 'Faxa_swndr', &
            is_local%wrap%FBImp(compatm,compocn), 'Faxa_swndr', &
            is_local%wrap%RH(compatm,compocn,mapconsf), rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    else
       call ESMF_LogWrite(trim(subname)//": ERROR RH not available for "//trim(mapnames(mapconsf)), &
            ESMF_LOGMSG_ERROR, line=__LINE__, file=u_FILE_u, rc=dbrc)
       rc = ESMF_FAILURE
    end if

    ! Calculate ocean albedos on the ocean grid

    update_alb = .false.
    lsize = size(anidr)

    if (flux_albav) then
       do n = 1,lsize
          anidr(n) = albdir
          avsdr(n) = albdir
          anidf(n) = albdif
          avsdf(n) = albdif

          ! Albedo is now function of latitude (will be new implementation)
          !rlat = const_deg2rad * lats(n)
          !anidr = 0.069_r8 - 0.011_r8 * cos(2._r8 * rlat)
          !avsdr = anidr
          !anidf = anidr
          !avsdf = anidr
       end do
       update_alb = .true.
    else
       ! need swdn & swup = swdn*(-albedo)
       ! swdn & albedos are time-aligned  BEFORE albedos get updated below ---

       do n=1,lsize
          if ( anidr(n) == 1.0_r8 ) then ! dark side of earth
             swup(n) = 0.0_r8
             swdn(n) = 0.0_r8
          else
             swup(n) = swndr(n) * (-anidr(n)) + swndf(n) * (-anidf(n)) + &
                       swvdr(n) * (-avsdr(n)) + swvdf(n) * (-avsdf(n))
             swdn(n) = swndr(n) + swndf(n) + swvdr(n) + swvdf(n)
          end if
       end do

       ! Solar declination
       ! Will only do albedo calculation if nextsw_cday is not -1.
       if (nextsw_cday >= -0.5_r8) then

          call shr_orb_decl(nextsw_cday, eccen, mvelpp,lambm0, obliqr, delta, eccf)

          ! Compute albedos
          do n = 1,lsize
             rlat = const_deg2rad * lats(n)
             rlon = const_deg2rad * lons(n)
             cosz = shr_orb_cosz( nextsw_cday, rlat, rlon, delta )
             if (cosz  >  0.0_r8) then !--- sun hit --
                anidr(n) = (.026_r8/(cosz**1.7_r8 + 0.065_r8)) +   &
                           (.150_r8*(cosz         - 0.100_r8 ) *   &
                           (cosz - 0.500_r8 ) * (cosz - 1.000_r8 )  )
                avsdr(n) = anidr(n)
                anidf(n) = albdif
                avsdf(n) = albdif
             else !--- dark side of earth ---
                anidr(n) = 1.0_r8
                avsdr(n) = 1.0_r8
                anidf(n) = 1.0_r8
                avsdf(n) = 1.0_r8
             end if
          end do
          update_alb = .true.
       endif    ! nextsw_cday
    end if   ! flux_albav

    ! If the aoflux grid is the atm - then need to map swdn and swup from the ocean grid (where
    ! they were calculated to the atmosphere grid)
    call NUOPC_CompAttributeGet(gcomp, name='aoflux_grid', value=cvalue, rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) aoflux_grid

    if (trim(aoflux_grid) == 'atm') then
       if (ESMF_RouteHandleIsCreated(is_local%wrap%RH(compocn,compatm,mapconsf), rc=rc)) then
          call shr_nuopc_methods_FB_FieldRegrid(&
               is_local%wrap%FBMed_ocnalb_o, 'Faox_swdn', &
               is_local%wrap%FBMed_ocnalb_a, 'Faox_swdn', &
               is_local%wrap%RH(compocn,compatm,mapconsf), rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

          call shr_nuopc_methods_FB_FieldRegrid(&
               is_local%wrap%FBMed_ocnalb_o, 'Faox_swup', &
               is_local%wrap%FBMed_ocnalb_a, 'Faox_swup', &
               is_local%wrap%RH(compocn,compatm,mapconsf), rc)
          if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       end if
    end if

    ! Update current ifrad/ofrad values if albedo was updated in field bundle
    if (update_alb) then
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compocn), fldname='ifrac', fldptr1=ifrac, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compocn), fldname='ifrad', fldptr1=ifrad, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compocn), fldname='ofrac', fldptr1=ofrac, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compocn), fldname='ofrad', fldptr1=ofrad, rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
       ifrad(:) = ifrac(:)
       ofrad(:) = ofrac(:)
    endif

    if (dbug_flag > 1) then
       call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBMed_ocnalb_o, string=trim(subname)//' FBMed_ocnalb_o', rc=rc)
       if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    end if

  end subroutine med_phases_ocnalb_run

  !===============================================================================

  subroutine med_phases_ocnalb_mapo2a(gcomp, rc)

    !----------------------------------------------------------
    ! Map ocean albedos from ocn to atm grid
    ! the med_phases routines
    !----------------------------------------------------------

    ! Arguments
    type(ESMF_GridComp)    :: gcomp
    integer, intent(out)   :: rc

    ! Local variables
    type(InternalState) :: is_local
    character(*), parameter :: subName =   '(med_ocnalb_mapo2a) '
    !-----------------------------------------------------------------------

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! Get the internal state from gcomp
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return

    ! Map the field bundle from the ocean to the atm grid
    call med_map_FB_Regrid_Norm( &
         fldListMed_ocnalb_o%flds, compocn, compatm, &
         is_local%wrap%FBMed_ocnalb_o, &
         is_local%wrap%FBMed_ocnalb_a, &
         is_local%wrap%FBFrac(compocn), &
         is_local%wrap%FBNormOne(compocn,compatm,:), &
         is_local%wrap%RH(compocn,compatm,:), &
         string='FBMed_ocnalb_o_To_FBMed_ocnalb_a', rc=rc)
    if (shr_nuopc_methods_ChkErr(rc,__LINE__,u_FILE_u)) return
    
  end subroutine med_phases_ocnalb_mapo2a

end module med_phases_ocnalb_mod
