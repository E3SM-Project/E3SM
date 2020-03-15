module docn_comp_mod

  use NUOPC                 , only : NUOPC_Advertise
  use ESMF                  , only : ESMF_SUCCESS, ESMF_LOGMSG_INFO, ESMF_LogWrite
  use ESMF                  , only : ESMF_State, ESMF_StateGet, ESMF_StateItem_Flag, ESMF_STATEITEM_NOTFOUND
  use ESMF                  , only : ESMF_Mesh, ESMF_MeshGet
  use ESMF                  , only : operator(/=), operator(==)
  use perf_mod              , only : t_startf, t_stopf, t_adj_detailf, t_barrierf
  use shr_kind_mod          , only : r8=>shr_kind_r8, cxx=>shr_kind_cxx, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod           , only : shr_sys_abort
  use shr_const_mod         , only : shr_const_cpsw, shr_const_rhosw, shr_const_TkFrz
  use shr_const_mod         , only : shr_const_TkFrzSw, shr_const_latice, shr_const_ocn_ref_sal
  use shr_const_mod         , only : shr_const_zsrflyr, shr_const_pi
  use shr_frz_mod           , only : shr_frz_freezetemp
  use shr_strdata_mod       , only : shr_strdata_type, shr_strdata_advance
  use dshr_methods_mod      , only : chkerr, state_getfldptr
  use dshr_dfield_mod       , only : dfield_type, dshr_dfield_add, dshr_dfield_copy
  use dshr_fldlist_mod      , only : fldlist_type, dshr_fldlist_add, dshr_fldlist_realize
  use dshr_nuopc_mod        , only : dshr_get_griddata

  implicit none
  private ! except

  !--------------------------------------------------------------------------
  ! Public/Private interfaces
  !--------------------------------------------------------------------------

  public :: docn_comp_advertise
  public :: docn_comp_realize
  public :: docn_comp_run

  !--------------------------------------------------------------------------
  ! Module data
  !--------------------------------------------------------------------------

  ! linked lists
  type(fldList_type) , pointer :: fldsImport => null()
  type(fldList_type) , pointer :: fldsExport => null()
  type(dfield_type)  , pointer :: dfields    => null()

  real(r8)     , parameter :: cpsw    = shr_const_cpsw        ! specific heat of sea h2o ~ j/kg/k
  real(r8)     , parameter :: rhosw   = shr_const_rhosw       ! density of sea water ~ kg/m^3
  real(r8)     , parameter :: tkfrz   = shr_const_tkfrz       ! freezing point, fresh water (kelvin)
  real(r8)     , parameter :: tkfrzsw = shr_const_tkfrzsw     ! freezing point, sea   water (kelvin)
  real(r8)     , parameter :: latice  = shr_const_latice      ! latent heat of fusion
  real(r8)     , parameter :: ocnsalt = shr_const_ocn_ref_sal ! ocean reference salinity

  ! prognostic flag set in docn_comp_advertise
  logical :: ocn_prognostic

  ! internal fields
  real(r8),         pointer :: xc(:), yc(:) ! mesh lats and lons - needed for aquaplanet analytical
  real(R8), public, pointer :: somtp(:)     ! SOM ocean temperature needed for restart

  ! export fields
  integer , pointer :: imask(:)     => null()    ! integer ocean mask
  real(r8), pointer :: So_omask(:)  => null()
  real(r8), pointer :: So_t(:)      => null()
  real(r8), pointer :: So_s(:)      => null()
  real(r8), pointer :: So_u(:)      => null()
  real(r8), pointer :: So_v(:)      => null()
  real(r8), pointer :: So_dhdx(:)   => null()
  real(r8), pointer :: So_dhdy(:)   => null()
  real(r8), pointer :: So_fswpen(:) => null()
  real(r8), pointer :: Fioo_q(:)    => null()

  ! import  fields
  real(r8), pointer :: Foxx_swnet(:) => null()
  real(r8), pointer :: Foxx_lwup(:)  => null()
  real(r8), pointer :: Foxx_sen(:)   => null()
  real(r8), pointer :: Foxx_lat(:)   => null()
  real(r8), pointer :: Faxa_lwdn(:)  => null()
  real(r8), pointer :: Faxa_snow(:)  => null()
  real(r8), pointer :: Fioi_melth(:) => null()
  real(r8), pointer :: Foxx_rofi(:)  => null()

  ! internal stream type
  real(r8), pointer :: strm_h(:)    => null()
  real(r8), pointer :: strm_qbot(:) => null()

  character(*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine docn_comp_advertise(importState, exportState, flds_scalar_name, ocn_prognostic_in, rc)

    ! --------------------------------------------------------------
    ! determine export and import fields to advertise to mediator
    ! --------------------------------------------------------------

    ! input/output arguments
    type(ESMF_State)     , intent(inout) :: importState
    type(ESMF_State)     , intent(inout) :: exportState
    character(len=*)     , intent(in)    :: flds_scalar_name
    logical              , intent(in)    :: ocn_prognostic_in
    integer              , intent(out)   :: rc

    ! local variables
    integer         :: n
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    ! Save as module variables for use in debugging
    ocn_prognostic = ocn_prognostic_in

    ! Advertise export fields

    call dshr_fldList_add(fldsExport, trim(flds_scalar_name))
    call dshr_fldList_add(fldsExport, 'So_omask'            )
    call dshr_fldList_add(fldsExport, 'So_t'                )
    call dshr_fldList_add(fldsExport, 'So_s'                )
    call dshr_fldList_add(fldsExport, 'So_u'                )
    call dshr_fldList_add(fldsExport, 'So_v'                )
    call dshr_fldList_add(fldsExport, 'So_dhdx'             )
    call dshr_fldList_add(fldsExport, 'So_dhdy'             )
    call dshr_fldList_add(fldsExport, 'Fioo_q'              )
    call dshr_fldList_add(fldsExport, 'So_fswpen'           )

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(docn_comp_advertise): Fr_ocn'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

    ! Advertise import fields

    if (ocn_prognostic) then
       call dshr_fldList_add(fldsImport, trim(flds_scalar_name))
       call dshr_fldList_add(fldsImport, 'Foxx_swnet'          )
       call dshr_fldList_add(fldsImport, 'Foxx_lwup'           )
       call dshr_fldList_add(fldsImport, 'Foxx_sen'            )
       call dshr_fldList_add(fldsImport, 'Foxx_lat'            )
       call dshr_fldList_add(fldsImport, 'Faxa_lwdn'           )
       call dshr_fldList_add(fldsImport, 'Faxa_snow'           )
       call dshr_fldList_add(fldsImport, 'Fioi_melth'          )
       call dshr_fldList_add(fldsImport, 'Foxx_rofi'           )

       fldlist => fldsImport ! the head of the linked list
       do while (associated(fldlist))
          call NUOPC_Advertise(importState, standardName=fldlist%stdname, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          call ESMF_LogWrite('(docn_comp_advertise): Fr_ocn'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
          fldList => fldList%next
       enddo
    end if

  end subroutine docn_comp_advertise

!===============================================================================

  subroutine docn_comp_realize(sdat, importState, exportState, flds_scalar_name, flds_scalar_num, mesh, &
       logunit, masterproc, rc)

    ! input/output parameters
    type(shr_strdata_type) , intent(inout) :: sdat
    type(ESMF_State)       , intent(inout) :: importState
    type(ESMF_State)       , intent(inout) :: exportState
    character(len=*)       , intent(in)    :: flds_scalar_name
    integer                , intent(in)    :: flds_scalar_num
    type(ESMF_Mesh)        , intent(in)    :: mesh
    integer                , intent(in)    :: logunit
    logical                , intent(in)    :: masterproc
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_StateItem_Flag) :: itemFlag
    integer                   :: n, lsize, kf
    character(*), parameter   :: subName = "(docn_comp_realize) "
    real(R8)    , parameter   :: &
         swp = 0.67_R8*(exp((-1._R8*shr_const_zsrflyr) /1.0_R8)) + 0.33_R8*exp((-1._R8*shr_const_zsrflyr)/17.0_R8)
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! -------------------------------------
    ! NUOPC_Realize "realizes" a previously advertised field in the importState and exportState
    ! by replacing the advertised fields with the newly created fields of the same name.
    ! -------------------------------------

    call dshr_fldlist_realize( exportState, fldsExport, flds_scalar_name, flds_scalar_num, mesh, &
         subname//':docnExport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call dshr_fldlist_realize( importState, fldsImport, flds_scalar_name, flds_scalar_num, mesh, &
         subname//':docnImport', rc=rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    ! -------------------------------------
    ! Set pointers to exportState fields
    ! -------------------------------------

    ! Set pointers to exportState fields that have no corresponding stream field
    call state_getfldptr(exportState, fldname='So_omask', fldptr1=So_omask, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_get_griddata(sdat, 'frac', So_omask)

    call state_getfldptr(exportState, fldname='Fioo_q', fldptr1=Fioo_q, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    Fioo_q(:) = 0._r8

    ! Initialize export state data that has corresponding stream field
    call dshr_dfield_add(dfields, sdat, state_fld='So_t', strm_fld='t', &
         state=exportState, state_ptr=So_t, logunit=logunit, masterproc=masterproc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat, state_fld='So_s', strm_fld='s', &
         state=exportState, state_ptr=So_s, logunit=logunit, masterproc=masterproc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat,  state_fld='So_u', strm_fld='u', &
         state=exportState, state_ptr=So_u, logunit=logunit, masterproc=masterproc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat,  state_fld='So_v', strm_fld='v', &
         state=exportState, state_ptr=So_v, logunit=logunit, masterproc=masterproc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat, state_fld='So_dhdx', strm_fld='dhdx', &
         state=exportState, state_ptr=So_dhdx, logunit=logunit, masterproc=masterproc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_dfield_add(dfields, sdat, state_fld='So_dhdy', strm_fld='dhdy', &
         state=exportState, state_ptr=So_dhdy, logunit=logunit, masterproc=masterproc, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! Initialize dfields stream fields that have no corresponding export fields
    call dshr_dfield_add(dfields, sdat,  strm_fld='qbot', strm_ptr=strm_qbot)
    call dshr_dfield_add(dfields, sdat,  strm_fld='h'   , strm_ptr=strm_h)

    ! For So_fswpen is only needed for diurnal cycle calculation of atm/ocn fluxes - and 
    ! currently this is not implemented in cmeps
    call ESMF_StateGet(exportState, 'So_fswpen', itemFlag, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
       call state_getfldptr(exportState, 'So_fswpen', fldptr1=So_fswpen, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       So_fswpen(:) = swp
    end if

    ! Initialize export state pointers to non-zero
    So_t(:) = TkFrz
    So_s(:) = ocnsalt

    ! determine module mask array (imask)
    lsize = size(So_omask)
    allocate(imask(lsize))
    imask = nint(So_omask)
    allocate(somtp(lsize))

    ! -------------------------------------
    ! Set pointers to importState fields
    ! -------------------------------------

    if (ocn_prognostic) then
       call state_getfldptr(importState, 'Foxx_swnet' , fldptr1=Foxx_swnet , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Foxx_lwup'  , fldptr1=Foxx_lwup  , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Foxx_lwup'  , fldptr1=Foxx_lwup  , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Foxx_sen'   , fldptr1=Foxx_sen   , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Foxx_lat'   , fldptr1=Foxx_lat   , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Faxa_lwdn'  , fldptr1=Faxa_lwdn  , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Faxa_snow'  , fldptr1=Faxa_snow  , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Fioi_melth' , fldptr1=Fioi_melth , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call state_getfldptr(importState, 'Foxx_rofi'  , fldptr1=Foxx_rofi  , rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
    end if

  end subroutine docn_comp_realize

!===============================================================================

  subroutine docn_comp_run(mpicom, my_task, master_task, logunit, target_ymd, target_tod, sdat, &
       mesh, sst_constant_value, dt, aquap_option, read_restart, rc)

    ! --------------------------
    ! advance docn
    ! --------------------------

    ! input/output variables:
    integer                , intent(in)    :: mpicom           ! mpi communicator
    integer                , intent(in)    :: my_task
    integer                , intent(in)    :: master_task
    integer                , intent(in)    :: logunit
    integer                , intent(in)    :: target_ymd       ! model date
    integer                , intent(in)    :: target_tod       ! model sec into model date
    type(shr_strdata_type) , intent(inout) :: sdat
    type(ESMF_Mesh)        , intent(in)    :: mesh
    real(r8)               , intent(in)    :: sst_constant_value
    real(r8)               , intent(in)    :: dt
    integer                , intent(in)    :: aquap_option
    logical                , intent(in)    :: read_restart
    integer                , intent(out)   :: rc

    ! local variables
    integer               :: n, lsize
    logical               :: first_time = .true.
    integer               :: spatialDim         ! number of dimension in mesh
    integer               :: numOwnedElements   ! size of mesh
    real(r8), pointer     :: ownedElemCoords(:) ! mesh lat and lons
    real(r8), allocatable :: tfreeze(:)         ! SOM ocean freezing temperature
    character(*), parameter :: subName = "(docn_comp_run) "
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call t_startf('DOCN_RUN')

    !--------------------
    ! advance docn streams
    !--------------------

    ! time and spatially interpolate to model time and grid
    call t_barrierf('docn_BARRIER',mpicom)
    call t_startf('docn_strdata_advance')
    call shr_strdata_advance(sdat, target_ymd, target_tod, mpicom, 'docn')
    call t_stopf('docn_strdata_advance')

    !--------------------
    ! copy all fields from streams to export state as default
    !--------------------

    ! This automatically will update the fields in the export state
    call t_barrierf('docn_dfield_copy_BARRIER', mpicom)
    call t_startf('docn_dfield_copy')
    call dshr_dfield_copy(dfields, sdat, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return
    call t_stopf('docn_dfield_copy')

    !-------------------------------------------------
    ! Determine additional data model behavior based on the mode
    !-------------------------------------------------

    lsize = size(So_t)

    call t_startf('docn_datamode')
    select case (trim(sdat%datamode))

    case('COPYALL')
       ! do nothing extra

    case('SSTDATA')
       So_t(:) = So_t(:) + TkFrz

    case('SST_AQUAPANAL')
       So_s(:)      = 0.0_r8
       if (associated(So_fswpen)) then
          So_fswpen(:) = 0.0_r8
       end if
       if (first_time) then
          call ESMF_MeshGet(mesh, spatialDim=spatialDim, numOwnedElements=numOwnedElements, rc=rc)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          allocate(ownedElemCoords(spatialDim*numOwnedElements))
          allocate(xc(numOwnedElements), yc(numOwnedElements))
          call ESMF_MeshGet(mesh, ownedElemCoords=ownedElemCoords)
          if (ChkErr(rc,__LINE__,u_FILE_u)) return
          do n = 1,numOwnedElements
             xc(n) = ownedElemCoords(2*n-1)
             yc(n) = ownedElemCoords(2*n)
          end do
       end if
       call prescribed_sst(xc, yc, lsize, aquap_option, So_t)
       So_t(:) = So_t(:) + TkFrz

    case('SST_AQUAPFILE')
       So_s(:)      = 0.0_r8
       if (associated(So_fswpen)) then
          So_fswpen(:) = 0.0_r8
       end if
       So_t(:) = So_t(:) + TkFrz

    case('SST_AQUAP_CONSTANT')
       So_s(:)      = 0.0_r8
       if (associated(So_fswpen)) then
          So_fswpen(:) = 0.0_r8
       end if
       So_t(:) = sst_constant_value

    case('IAF')
       So_t(:) = So_t(:) + TkFrz

    case('SOM')
       if (first_time) then
          do n = 1,lsize
             if (.not. read_restart) then
                somtp(n) = So_t(n) + TkFrz
             endif
             So_t(n) = somtp(n)
             Fioo_q(n) = 0.0_R8
          enddo
       else
          allocate(tfreeze(lsize))
          tfreeze(:) = shr_frz_freezetemp(So_s(:)) + TkFrz
          do n = 1,lsize
             if (imask(n) /= 0) then
                ! compute new temp (last term is latent by prec and roff)
                So_t(n) = somtp(n) +  &
                     ( Foxx_swnet(n) + Foxx_lwup(n) + Faxa_lwdn(n) + Foxx_sen(n) + Foxx_lat(n) + &
                       Fioi_melth(n) - strm_qbot(n) - (Faxa_snow(n)+Foxx_rofi(n))*latice) * dt/(cpsw*rhosw* strm_h(n))

                ! compute ice formed or melt potential
                Fioo_q(n) = (tfreeze(n) - So_t(n))*(cpsw*rhosw*strm_h(n))/dt ! ice formed q>0

                ! reset temp
                So_t(n)  = max(tfreeze(n),So_t(n))

                ! save somtp to restart file
                somtp(n) = So_t(n)
             endif
          end do
          deallocate(tfreeze)
       endif   ! first_time

    case('SOM_AQUAP')
       if (first_time) then
          do n = 1,lsize
             if (.not. read_restart) then
                somtp(n) = So_t(n) + TkFrz
             endif
             So_t(n) = somtp(n)
             Fioo_q(n) = 0.0_R8
          enddo
       else
          allocate(tfreeze(lsize))
          tfreeze(:) = shr_frz_freezetemp(So_s(:)) + TkFrz
          do n = 1,lsize
             ! compute new temp (last term is latent by prec and roff)
             So_t(n) = somtp(n) + &
                  ( Foxx_swnet(n) + Foxx_lwup(n) + Faxa_lwdn(n) + Foxx_sen(n) + Foxx_lat(n) + &
                    Fioi_melth(n) - strm_qbot(n) - (Faxa_snow(n)+Foxx_rofi(n))*latice) * dt/(cpsw*rhosw*strm_h(n))

             ! compute ice formed or melt potential
             Fioo_q(n) = (tfreeze(n) - So_t(n))*(cpsw*rhosw*strm_h(n))/dt  ! ice formed q>0

             ! save somtp on restart file
             somtp(n) = So_t(n)
          enddo
          deallocate(tfreeze)
       endif   ! first_time

    end select

    first_time= .false.

    call t_stopf('docn_datamode')
    call t_stopf('DOCN_RUN')

  end subroutine docn_comp_run

!===============================================================================

  subroutine prescribed_sst(xc, yc, lsize, sst_option, sst)

    ! input/output variables
    real(R8)     , intent(in)    :: xc(:)  !degrees
    real(R8)     , intent(in)    :: yc(:)  !degrees
    integer      , intent(in)    :: lsize
    integer      , intent(in)    :: sst_option
    real(R8)     , intent(inout) :: sst(:)

    ! local variables
    integer  :: i
    real(r8) :: tmp, tmp1, pi
    real(r8) :: rlon(lsize), rlat(lsize)

    real(r8), parameter :: pio180 = SHR_CONST_PI/180._r8

    ! Parameters for zonally symmetric experiments
    real(r8), parameter ::   t0_max     = 27._r8
    real(r8), parameter ::   t0_min     = 0._r8
    real(r8), parameter ::   maxlat     = 60._r8*pio180
    real(r8), parameter ::   shift      = 5._r8*pio180
    real(r8), parameter ::   shift9     = 10._r8*pio180
    real(r8), parameter ::   shift10    = 15._r8*pio180

    ! Parameters for zonally asymmetric experiments
    real(r8), parameter ::   t0_max6    = 1._r8
    real(r8), parameter ::   t0_max7    = 3._r8
    real(r8), parameter ::   latcen     = 0._r8*pio180
    real(r8), parameter ::   loncen     = 0._r8*pio180
    real(r8), parameter ::   latrad6    = 15._r8*pio180
    real(r8), parameter ::   latrad8    = 30._r8*pio180
    real(r8), parameter ::   lonrad     = 30._r8*pio180
    !-------------------------------------------------------------------------------

    pi = SHR_CONST_PI

    ! convert xc and yc from degrees to radians

    rlon(:) = xc(:) * pio180
    rlat(:) = yc(:) * pio180

    ! Control
    if (sst_option < 1 .or. sst_option > 10) then
       call shr_sys_abort ('prescribed_sst: ERROR: sst_option must be between 1 and 10')
    end if

    if (sst_option == 1 .or. sst_option == 6 .or. sst_option == 7 .or. sst_option == 8) then
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if
    if (sst_option == 2) then ! Flat
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
             tmp = 1._r8 - tmp*tmp*tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if
    if (sst_option == 3) then ! Qobs
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = sin(rlat(i)*pi*0.5_r8/maxlat)
             tmp = (2._r8 - tmp*tmp*tmp*tmp - tmp*tmp)*0.5_r8
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if
    if (sst_option == 4) then ! Peaked
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else
             tmp = (maxlat - abs(rlat(i)))/maxlat
             tmp1 = 1._r8 - tmp
             sst(i) = t0_max*tmp + t0_min*tmp1
          end if
       end do
    end if
    if (sst_option == 5) then ! Control-5N
       do i = 1,lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else if (rlat(i) > shift) then
             tmp = sin((rlat(i)-shift)*pi*0.5_r8/(maxlat-shift))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          else
             tmp = sin((rlat(i)-shift)*pi*0.5_r8/(maxlat+shift))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if
    if (sst_option == 6) then ! 1KEQ
       do i = 1,lsize
          if (abs(rlat(i)-latcen) <= latrad6) then
             tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad6)
             tmp1 = tmp1*tmp1
             tmp = abs(rlon(i)-loncen)
             tmp = min(tmp , 2._r8*pi-tmp)
             if(tmp <= lonrad) then
                tmp = cos(tmp*pi*0.5_r8/lonrad)
                tmp = tmp*tmp
                sst(i) = sst(i) + t0_max6*tmp*tmp1
             end if
          end if
       end do
    end if
    if (sst_option == 7) then ! 3KEQ
       do i = 1, lsize
          if (abs(rlat(i)-latcen) <= latrad6) then
             tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad6)
             tmp1 = tmp1*tmp1
             tmp = abs(rlon(i)-loncen)
             tmp = min(tmp , 2._r8*pi-tmp)
             if (tmp <= lonrad) then
                tmp = cos(tmp*pi*0.5_r8/lonrad)
                tmp = tmp*tmp
                sst(i) = sst(i) + t0_max7*tmp*tmp1
             end if
          end if
       end do
    end if
    if (sst_option == 8) then ! 3KW1
       do i = 1, lsize
          if (abs(rlat(i)-latcen) <= latrad8) then
             tmp1 = cos((rlat(i)-latcen)*pi*0.5_r8/latrad8)
             tmp1 = tmp1*tmp1
             tmp = cos(rlon(i)-loncen)
             sst(i) = sst(i) + t0_max7*tmp*tmp1
          end if
       end do
    end if
    if (sst_option == 9) then ! Control-10N
       do i = 1, lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else if (rlat(i) > shift9) then
             tmp = sin((rlat(i)-shift9)*pi*0.5_r8/(maxlat-shift9))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          else
             tmp = sin((rlat(i)-shift9)*pi*0.5_r8/(maxlat+shift9))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if
    if (sst_option == 10) then ! Control-15N
       do i = 1, lsize
          if (abs(rlat(i)) > maxlat) then
             sst(i) = t0_min
          else if(rlat(i) > shift10) then
             tmp = sin((rlat(i)-shift10)*pi*0.5_r8/(maxlat-shift10))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          else
             tmp = sin((rlat(i)-shift10)*pi*0.5_r8/(maxlat+shift10))
             tmp = 1._r8 - tmp*tmp
             sst(i) = tmp*(t0_max - t0_min) + t0_min
          end if
       end do
    end if

  end subroutine prescribed_sst

end module docn_comp_mod
