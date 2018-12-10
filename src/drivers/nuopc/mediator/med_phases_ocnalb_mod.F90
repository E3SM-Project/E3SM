module med_phases_ocnalb_mod

  use med_constants_mod, only : R8

  implicit none
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public med_phases_ocnalb_run
  public med_phases_ocnalb_mapo2a

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private med_phases_ocnalb_init

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  type ocnalb_type
     real(r8) , pointer :: lats  (:) ! latitudes  (degrees)
     real(r8) , pointer :: lons  (:) ! longitudes (degrees)
     integer  , pointer :: mask  (:) ! ocn domain mask: 0 <=> inactive cell
     real(r8) , pointer :: anidr (:) ! albedo: near infrared, direct
     real(r8) , pointer :: avsdr (:) ! albedo: visible      , direct
     real(r8) , pointer :: anidf (:) ! albedo: near infrared, diffuse
     real(r8) , pointer :: avsdf (:) ! albedo: visible      , diffuse
  end type ocnalb_type

  ! Conversion from degrees to radians
  character(*),parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine med_phases_ocnalb_init(gcomp, ocnalb, rc)

    !-----------------------------------------------------------------------
    ! Initialize pointers to the module variables and then use the module
    ! variables in the med_ocnalb phase
    ! All input field bundles are ASSUMED to be on the ocean grid
    !-----------------------------------------------------------------------

    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_FAILURE
    use ESMF                  , only : ESMF_GridComp, ESMF_VM, ESMF_Field, ESMF_Grid, ESMF_Mesh, ESMF_GeomType_Flag
    use ESMF                  , only : ESMF_GridCompGet, ESMF_VMGet, ESMF_FieldGet, ESMF_GEOMTYPE_MESH
    use ESMF                  , only : ESMF_MeshGet
    use ESMF                  , only : operator(==)
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_getFieldN
    use shr_nuopc_utils_mod , only : shr_nuopc_utils_chkerr
    use med_internalstate_mod , only : InternalState
    use med_constants_mod     , only : CL, R8
    use med_constants_mod     , only : dbug_flag =>med_constants_dbug_flag
    use esmFlds               , only : compatm, compocn
    use perf_mod              , only : t_startf, t_stopf
    ! Arguments
    type(ESMF_GridComp)               :: gcomp
    type(ocnalb_type) , intent(inout) :: ocnalb
    integer           , intent(out)   :: rc
    !
    ! Local variables
    type(ESMF_VM)            :: vm
    integer                  :: iam
    type(ESMF_Field)         :: lfield
    type(ESMF_Mesh)          :: lmesh
    type(ESMF_GeomType_Flag) :: geomtype
    integer                  :: n
    integer                  :: lsize
    integer                  :: dimCount
    integer                  :: spatialDim
    integer                  :: numOwnedElements
    type(InternalState)      :: is_local
    real(R8), pointer        :: ownedElemCoords(:)
    character(len=CL)        :: tempc1,tempc2
    integer                  :: dbrc
    character(*), parameter  :: subname = '(med_phases_ocnalb_init) '
    !-----------------------------------------------------------------------
    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! The following is for debugging
    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    call ESMF_VMGet(vm, localPet=iam, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    ! Get the internal state from gcomp
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Set pointers to fields needed for albedo calculations
    !----------------------------------

    ! These must must be on the ocean grid since the ocean albedo computation is on the ocean grid
    ! The following sets pointers to the module arrays

    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, fldname='So_avsdr', fldptr1=ocnalb%avsdr, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, fldname='So_avsdf', fldptr1=ocnalb%avsdf, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, fldname='So_anidr', fldptr1=ocnalb%anidr, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    call shr_nuopc_methods_FB_GetFldPtr(is_local%wrap%FBMed_ocnalb_o, fldname='So_anidf', fldptr1=ocnalb%anidf, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    !----------------------------------
    ! Get lat, lon, which are time-invariant
    !----------------------------------

    ! The following assumes that all fields in FBMed_ocnalb_o have the same grid - so
    ! only need to query field 1
    call shr_nuopc_methods_FB_getFieldN(is_local%wrap%FBMed_ocnalb_o, fieldnum=1, field=lfield, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    ! Determine if first field is on a grid or a mesh - default will be mesh
    call ESMF_FieldGet(lfield, geomtype=geomtype, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    if (geomtype == ESMF_GEOMTYPE_MESH) then
       call ESMF_LogWrite(trim(subname)//" : FBAtm is on a mesh ", ESMF_LOGMSG_INFO, rc=rc)
       call ESMF_FieldGet(lfield, mesh=lmesh, rc=rc)
       if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_MeshGet(lmesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
       if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
       lsize = size(ocnalb%anidr)
       if (numOwnedElements /= lsize) then
          write(tempc1,'(i10)') numOwnedElements
          write(tempc2,'(i10)') lsize
          call ESMF_LogWrite(trim(subname)//": ERROR numOwnedElements "// trim(tempc1) // &
               " not equal to local size "// trim(tempc2), ESMF_LOGMSG_INFO, rc=rc)
          rc = ESMF_FAILURE
          return
       end if
       allocate(ownedElemCoords(spatialDim*numOwnedElements))
       allocate(ocnalb%lons(numOwnedElements))
       allocate(ocnalb%lats(numOwnedElements))
       call ESMF_MeshGet(lmesh, ownedElemCoords=ownedElemCoords)
       if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
       do n = 1,lsize
          ocnalb%lons(n) = ownedElemCoords(2*n-1)
          ocnalb%lats(n) = ownedElemCoords(2*n)
       end do
    else
      call ESMF_LogWrite(trim(subname)//": ERROR field bundle must be either on mesh", ESMF_LOGMSG_INFO, rc=rc)
      rc = ESMF_FAILURE
      return
    end if

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": done", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    call t_stopf('MED:'//subname)

  end subroutine med_phases_ocnalb_init

  !===============================================================================

  subroutine med_phases_ocnalb_run(gcomp, rc)

    !-----------------------------------------------------------------------
    ! Compute ocean albedos (on the ocean grid)
    !-----------------------------------------------------------------------

    use ESMF                  , only : ESMF_GridComp, ESMF_Clock, ESMF_Time, ESMF_TimeInterval
    use ESMF                  , only : ESMF_GridCompGet, ESMF_ClockGet, ESMF_TimeGet
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS, ESMF_LogFoundError
    use ESMF                  , only : ESMF_RouteHandleIsCreated
    use ESMF                  , only : operator(+)
    use NUOPC                 , only : NUOPC_CompAttributeGet
    use shr_const_mod         , only : shr_const_pi
    use shr_sys_mod           , only : shr_sys_abort
    use shr_orb_mod           , only : shr_orb_cosz, shr_orb_decl
    use shr_nuopc_fldList_mod , only : mapconsf, mapnames
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_GetFldPtr
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_diagnose
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_State_GetScalar
    use shr_nuopc_methods_mod , only : shr_nuopc_methods_FB_FieldRegrid
    use shr_nuopc_utils_mod , only : shr_nuopc_utils_chkerr
    use med_constants_mod     , only : CS, CL, R8
    use med_constants_mod     , only : dbug_flag =>med_constants_dbug_flag
    use med_internalstate_mod , only : InternalState, logunit
    use shr_nuopc_scalars_mod , only : flds_scalar_name
    use shr_nuopc_scalars_mod , only : flds_scalar_num
    use shr_nuopc_scalars_mod , only : flds_scalar_index_nextsw_cday
    use esmFlds               , only : compatm, compocn
    use perf_mod              , only : t_startf, t_stopf
    ! input/output variables
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    ! local variables
    type(ocnalb_type), save :: ocnalb
    logical                 :: update_alb
    type(InternalState)     :: is_local
    type(ESMF_Clock)        :: clock
    type(ESMF_Time)         :: currTime
    type(ESMF_Time)         :: nextTime
    type(ESMF_TimeInterval) :: timeStep
    character(CL)           :: cvalue
    character(CS)           :: starttype        ! config start type
    character(CL)           :: runtype          ! initial, continue, hybrid, branch
    character(CL)           :: aoflux_grid
    logical                 :: flux_albav       ! flux avg option
    real(R8)                :: nextsw_cday      ! calendar day of next atm shortwave
    real(R8), pointer       :: ofrac(:)
    real(R8), pointer       :: ofrad(:)
    real(R8), pointer       :: ifrac(:)
    real(R8), pointer       :: ifrad(:)
    integer                 :: lsize            ! local size
    integer                 :: n,i              ! indices
    real(R8)                :: rlat             ! gridcell latitude in radians
    real(R8)                :: rlon             ! gridcell longitude in radians
    real(R8)                :: cosz             ! Cosine of solar zenith angle
    real(R8)                :: eccen            ! Earth orbit eccentricity
    real(R8)                :: mvelpp           ! Earth orbit
    real(R8)                :: lambm0           ! Earth orbit
    real(R8)                :: obliqr           ! Earth orbit
    real(R8)                :: delta            ! Solar declination angle  in radians
    real(R8)                :: eccf             ! Earth orbit eccentricity factor
    real(R8), parameter     :: albdif = 0.06_r8 ! 60 deg reference albedo, diffuse
    real(R8), parameter     :: albdir = 0.07_r8 ! 60 deg reference albedo, direct
    real(R8), parameter     :: const_deg2rad = shr_const_pi/180.0_R8  ! deg to rads
    integer                 :: dbrc
    logical                 :: first_call = .true.
    character(len=*)  , parameter :: subname='(med_phases_ocnalb_run)'
    !---------------------------------------
    call t_startf('MED:'//subname)
    if (dbug_flag > 5) then
       call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! Get the internal state from Component.
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    ! Note that in the mct version the atm was initialized first so
    ! that nextsw_cday could be passed to the other components - this
    ! assumed that atmosphere component was ALWAYS initialized first.
    ! In the nuopc version it will be easier to assume that on startup
    ! - nextsw_cday is just what cam was setting it as the current calendar day

    if (first_call) then

       ! Initialize ocean albedo calculation
       call med_phases_ocnalb_init(gcomp, ocnalb, rc)
       if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

       call NUOPC_CompAttributeGet(gcomp, name='start_type', value=cvalue, rc=rc)
       if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
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

       call ESMF_GridCompGet(gcomp, clock=clock)
       if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_ClockGet( clock, currTime=currTime, timeStep=timeStep, rc=rc)
       if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

       if (trim(runtype) == 'initial') then
          call ESMF_TimeGet( currTime, dayOfYear_r8=nextsw_cday, rc=rc )
          if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
       else
          call shr_nuopc_methods_State_GetScalar(is_local%wrap%NstateImp(compatm), &
               flds_scalar_name=flds_scalar_name, flds_scalar_num=flds_scalar_num, &
               scalar_id=flds_scalar_index_nextsw_cday, value=nextsw_cday, rc=rc)
          if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
       end if

       first_call = .false.

    else

       ! Note that shr_nuopc_methods_State_GetScalar includes a broadcast to all other pets
       call shr_nuopc_methods_State_GetScalar(is_local%wrap%NstateImp(compatm), &
            flds_scalar_name=flds_scalar_name, flds_scalar_num=flds_scalar_num, &
            scalar_id=flds_scalar_index_nextsw_cday, value=nextsw_cday, rc=rc)
       if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    end if

    call NUOPC_CompAttributeGet(gcomp, name='flux_albav', value=cvalue, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) flux_albav

    call NUOPC_CompAttributeGet(gcomp, name='orb_eccen', value=cvalue, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) eccen
    call NUOPC_CompAttributeGet(gcomp, name='orb_obliqr', value=cvalue, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) obliqr
    call NUOPC_CompAttributeGet(gcomp, name='orb_lambm0', value=cvalue, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) lambm0
    call NUOPC_CompAttributeGet(gcomp, name='orb_mvelpp', value=cvalue, rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    read(cvalue,*) mvelpp

    ! Calculate ocean albedos on the ocean grid

    update_alb = .false.
    lsize = size(ocnalb%anidr)

    if (flux_albav) then
       do n = 1,lsize
          ocnalb%anidr(n) = albdir
          ocnalb%avsdr(n) = albdir
          ocnalb%anidf(n) = albdif
          ocnalb%avsdf(n) = albdif
       end do
       update_alb = .true.
    else
       ! Solar declination
       ! Will only do albedo calculation if nextsw_cday is not -1.
       if (nextsw_cday >= -0.5_r8) then

          call shr_orb_decl(nextsw_cday, eccen, mvelpp,lambm0, obliqr, delta, eccf)

          ! Compute albedos
          do n = 1,lsize
             rlat = const_deg2rad * ocnalb%lats(n)
             rlon = const_deg2rad * ocnalb%lons(n)
             cosz = shr_orb_cosz( nextsw_cday, rlat, rlon, delta )
             if (cosz  >  0.0_r8) then !--- sun hit --
                ocnalb%anidr(n) = (.026_r8/(cosz**1.7_r8 + 0.065_r8)) +   &
                                  (.150_r8*(cosz         - 0.100_r8 ) *   &
                                  (cosz - 0.500_r8 ) * (cosz - 1.000_r8 )  )
                ocnalb%avsdr(n) = ocnalb%anidr(n)
                ocnalb%anidf(n) = albdif
                ocnalb%avsdf(n) = albdif
             else !--- dark side of earth ---
                ocnalb%anidr(n) = 1.0_r8
                ocnalb%avsdr(n) = 1.0_r8
                ocnalb%anidf(n) = 1.0_r8
                ocnalb%avsdf(n) = 1.0_r8
             end if
          end do
          update_alb = .true.

       endif    ! nextsw_cday
    end if   ! flux_albav

    ! Update current ifrad/ofrad values if albedo was updated in field bundle
    if (update_alb) then
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compocn), fldname='ifrac', fldptr1=ifrac, rc=rc)
       if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compocn), fldname='ifrad', fldptr1=ifrad, rc=rc)
       if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compocn), fldname='ofrac', fldptr1=ofrac, rc=rc)
       if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
       call shr_nuopc_methods_FB_getFldPtr(is_local%wrap%FBfrac(compocn), fldname='ofrad', fldptr1=ofrad, rc=rc)
       if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
       ifrad(:) = ifrac(:)
       ofrad(:) = ofrac(:)
    endif

    if (dbug_flag > 1) then
       call shr_nuopc_methods_FB_diagnose(is_local%wrap%FBMed_ocnalb_o, string=trim(subname)//' FBMed_ocnalb_o', rc=rc)
       if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    end if
    call t_stopf('MED:'//subname)

  end subroutine med_phases_ocnalb_run

  !===============================================================================

  subroutine med_phases_ocnalb_mapo2a(gcomp, rc)

    !----------------------------------------------------------
    ! Map ocean albedos from ocn to atm grid
    !----------------------------------------------------------

    use ESMF                  , only : ESMF_GridComp
    use ESMF                  , only : ESMF_LogWrite, ESMF_LOGMSG_INFO, ESMF_SUCCESS
    use shr_nuopc_utils_mod , only : shr_nuopc_utils_chkerr
    use med_map_mod           , only : med_map_FB_Regrid_Norm
    use med_internalstate_mod , only : InternalState
    use med_constants_mod     , only : R8
    use med_constants_mod     , only : dbug_flag =>med_constants_dbug_flag
    use esmFlds               , only : fldListMed_ocnalb_o
    use esmFlds               , only : compatm, compocn
    use perf_mod              , only : t_startf, t_stopf
    ! Arguments
    type(ESMF_GridComp)    :: gcomp
    integer, intent(out)   :: rc

    ! Local variables
    type(InternalState) :: is_local
    integer             :: dbrc
    character(*), parameter :: subName =   '(med_ocnalb_mapo2a) '
    !-----------------------------------------------------------------------
    call t_startf('MED:'//subname)

    if (dbug_flag > 5) then
      call ESMF_LogWrite(trim(subname)//": called", ESMF_LOGMSG_INFO, rc=dbrc)
    endif
    rc = ESMF_SUCCESS

    ! Get the internal state from gcomp
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return

    ! Map the field bundle from the ocean to the atm grid
    call med_map_FB_Regrid_Norm( &
         fldListMed_ocnalb_o%flds, compocn, compatm, &
         is_local%wrap%FBMed_ocnalb_o, &
         is_local%wrap%FBMed_ocnalb_a, &
         is_local%wrap%FBFrac(compocn), &
         is_local%wrap%FBNormOne(compocn,compatm,:), &
         is_local%wrap%RH(compocn,compatm,:), &
         string='FBMed_ocnalb_o_To_FBMed_ocnalb_a', rc=rc)
    if (shr_nuopc_utils_chkerr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('MED:'//subname)

  end subroutine med_phases_ocnalb_mapo2a

end module med_phases_ocnalb_mod
