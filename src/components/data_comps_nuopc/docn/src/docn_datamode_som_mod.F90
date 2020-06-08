module docn_datamode_som_mod

  use ESMF
  use NUOPC            , only : NUOPC_Advertise
  use shr_kind_mod     , only : r8=>shr_kind_r8, i8=>shr_kind_i8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_sys_mod      , only : shr_sys_abort
  use shr_cal_mod      , only : shr_cal_date2julian
  use shr_const_mod    , only : shr_const_cpsw, shr_const_rhosw, shr_const_TkFrz
  use shr_const_mod    , only : shr_const_TkFrzSw, shr_const_latice, shr_const_ocn_ref_sal
  use shr_const_mod    , only : shr_const_zsrflyr, shr_const_pi
  use shr_frz_mod      , only : shr_frz_freezetemp
  use dshr_strdata_mod , only : shr_strdata_get_stream_pointer, shr_strdata_type
  use dshr_methods_mod , only : dshr_state_getfldptr, dshr_fldbun_getfldptr, chkerr
  use dshr_strdata_mod , only : shr_strdata_type
  use dshr_mod         , only : dshr_restart_read, dshr_restart_write
  use dshr_fldlist_mod , only : fldlist_type, dshr_fldlist_add
  use pio

  implicit none
  private ! except

  public :: docn_datamode_som_advertise
  public :: docn_datamode_som_init_pointers
  public :: docn_datamode_som_advance
  public :: docn_datamode_som_restart_read
  public :: docn_datamode_som_restart_write

  ! export fields
  real(r8), pointer :: So_omask(:)  => null()    ! real ocean fraction sent to mediator
  real(r8), pointer :: So_t(:)      => null()
  real(r8), pointer :: So_s(:)      => null()
  real(r8), pointer :: So_u(:)      => null()
  real(r8), pointer :: So_v(:)      => null()
  real(r8), pointer :: So_dhdx(:)   => null()
  real(r8), pointer :: So_dhdy(:)   => null()
  real(r8), pointer :: Fioo_q(:)    => null()
  real(r8), pointer :: So_fswpen(:) => null()

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

  ! restart fields
  real(R8), public, pointer :: somtp(:)     ! SOM ocean temperature needed for restart

  real(R8) :: dt                            ! real model timestep

  ! constants
  real(r8) , parameter :: cpsw    = shr_const_cpsw        ! specific heat of sea h2o ~ j/kg/k
  real(r8) , parameter :: rhosw   = shr_const_rhosw       ! density of sea water ~ kg/m^3
  real(r8) , parameter :: tkfrz   = shr_const_tkfrz       ! freezing point, fresh water (kelvin)
  real(r8) , parameter :: tkfrzsw = shr_const_tkfrzsw     ! freezing point, sea   water (kelvin)
  real(r8) , parameter :: latice  = shr_const_latice      ! latent heat of fusion
  real(r8) , parameter :: ocnsalt = shr_const_ocn_ref_sal ! ocean reference salinity

  character(*) , parameter :: nullstr = 'undefined'
  character(*) , parameter :: rpfile  = 'rpointer.ocn'
  character(*) , parameter :: u_FILE_u = &
       __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine docn_datamode_som_advertise(importState, exportState, fldsimport, fldsexport, flds_scalar_name, rc)

    ! input/output variables
    type(esmf_State)   , intent(inout) :: importState
    type(esmf_State)   , intent(inout) :: exportState
    type(fldlist_type) , pointer       :: fldsimport
    type(fldlist_type) , pointer       :: fldsexport
    character(len=*)   , intent(in)    :: flds_scalar_name
    integer            , intent(out)   :: rc

    ! local variables
    type(fldlist_type), pointer :: fldList
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

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

    ! Advertise import fields
    call dshr_fldList_add(fldsImport, trim(flds_scalar_name))
    call dshr_fldList_add(fldsImport, 'Foxx_swnet'          )
    call dshr_fldList_add(fldsImport, 'Foxx_lwup'           )
    call dshr_fldList_add(fldsImport, 'Foxx_sen'            )
    call dshr_fldList_add(fldsImport, 'Foxx_lat'            )
    call dshr_fldList_add(fldsImport, 'Faxa_lwdn'           )
    call dshr_fldList_add(fldsImport, 'Faxa_snow'           )
    call dshr_fldList_add(fldsImport, 'Fioi_melth'          )
    call dshr_fldList_add(fldsImport, 'Foxx_rofi'           )

    fldlist => fldsExport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(exportState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(docn_comp_advertise): Fr_ocn'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

    ! Advertise import fields
    fldlist => fldsImport ! the head of the linked list
    do while (associated(fldlist))
       call NUOPC_Advertise(importState, standardName=fldlist%stdname, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_LogWrite('(docn_comp_advertise): Fr_ocn'//trim(fldList%stdname), ESMF_LOGMSG_INFO)
       fldList => fldList%next
    enddo

  end subroutine docn_datamode_som_advertise

  !===============================================================================
  subroutine docn_datamode_som_init_pointers(importState, exportState, sdat, ocn_fraction, rc)

    ! input/output variables
    type(ESMF_State)       , intent(inout) :: exportState
    type(ESMF_State)       , intent(inout) :: importState
    type(shr_strdata_type) , intent(in)    :: sdat
    real(r8)               , intent(in)    :: ocn_fraction(:)
    integer                , intent(out)   :: rc

    ! local variables
    type(ESMF_StateItem_Flag)   :: itemFlag
    character(len=*), parameter :: subname='(docn_init_pointers): '
    real(R8)        , parameter   :: &
         swp = 0.67_R8*(exp((-1._R8*shr_const_zsrflyr) /1.0_R8)) + 0.33_R8*exp((-1._R8*shr_const_zsrflyr)/17.0_R8)
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! initialize pointers to stream fields
    call shr_strdata_get_stream_pointer( sdat, 'So_qbot', strm_qbot, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call shr_strdata_get_stream_pointer( sdat, 'So_h'   , strm_h   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! initialize pointers to import fields
    call dshr_state_getfldptr(importState, 'Foxx_swnet' , fldptr1=Foxx_swnet , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, 'Foxx_lwup'  , fldptr1=Foxx_lwup  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, 'Foxx_lwup'  , fldptr1=Foxx_lwup  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, 'Foxx_sen'   , fldptr1=Foxx_sen   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, 'Foxx_lat'   , fldptr1=Foxx_lat   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, 'Faxa_lwdn'  , fldptr1=Faxa_lwdn  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, 'Faxa_snow'  , fldptr1=Faxa_snow  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, 'Fioi_melth' , fldptr1=Fioi_melth , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(importState, 'Foxx_rofi'  , fldptr1=Foxx_rofi  , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! initialize pointers to export fields
    call dshr_state_getfldptr(exportState, 'So_omask'   , fldptr1=So_omask   , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'So_t'       , fldptr1=So_t       , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'So_s'       , fldptr1=So_s       , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'So_u'       , fldptr1=So_u       , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'So_v'       , fldptr1=So_v       , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'So_dhdx'    , fldptr1=So_dhdx    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'So_dhdy'    , fldptr1=So_dhdy    , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    call dshr_state_getfldptr(exportState, 'Fioo_q'     , fldptr1=Fioo_q     , rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ! For So_fswpen is only needed for diurnal cycle calculation of atm/ocn fluxes 
    ! Currently this is not implemented in cmeps
    call ESMF_StateGet(exportState, 'So_fswpen', itemFlag, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    if (itemFlag /= ESMF_STATEITEM_NOTFOUND) then
       call dshr_state_getfldptr(exportState, 'So_fswpen', fldptr1=So_fswpen, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       So_fswpen(:) = swp
    end if

    ! Set export state ocean fraction (So_omask)
    So_omask(:) = ocn_fraction(:)

    ! Allocate memory for somtp
    allocate(somtp(sdat%model_lsize))

    ! Initialize export state pointers to non-zero
    So_t(:) = TkFrz
    So_s(:) = ocnsalt

  end subroutine docn_datamode_som_init_pointers

  !===============================================================================
  subroutine docn_datamode_som_advance(importState, exportState, clock, restart_read, datamode, rc)

    ! input/output variables
    type(ESMF_State)       , intent(inout) :: importState
    type(ESMF_State)       , intent(inout) :: exportState
    type(ESMF_Clock)       , intent(in)    :: clock 
    logical                , intent(in)    :: restart_read
    character(len=*)       , intent(in)    :: datamode
    integer                , intent(out)   :: rc

    ! local variables
    logical                 :: first_time = .true.
    type(ESMF_TimeInterval) :: timeStep
    integer                 :: idt          ! integer model timestep
    real(r8), allocatable   :: tfreeze(:)         ! SOM ocean freezing temperature
    integer                 :: lsize
    integer                 :: n
    logical                 :: reset_temp
    character(len=*), parameter :: subname='(docn_datamode_som): '
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    lsize = size(So_t)

    if (first_time) then

       ! set model timestep
       call ESMF_ClockGet( clock, timeStep=timeStep, rc=rc)
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       call ESMF_TimeIntervalGet( timeStep, s=idt, rc=rc )
       if (ChkErr(rc,__LINE__,u_FILE_u)) return
       dt = idt * 1.0_r8

       do n = 1,lsize
          if (.not. restart_read) then
             somtp(n) = So_t(n) + TkFrz
          endif
          So_t(n) = somtp(n)
          Fioo_q(n) = 0.0_R8
       enddo

    else

       ! determine if will reset temp
       if (datamode == 'SOM_AQUAP') then
          reset_temp = .false.
       else
          reset_temp = .true.
       end if

       allocate(tfreeze(lsize))
       tfreeze(:) = shr_frz_freezetemp(So_s(:)) + TkFrz
       do n = 1,lsize
          if (So_omask(n) /= 0._r8) then
             ! compute new temp (last term is latent by prec and roff)
             So_t(n) = somtp(n) +  &
                  ( Foxx_swnet(n) + Foxx_lwup(n) + Faxa_lwdn(n) + Foxx_sen(n) + Foxx_lat(n) + &
                  Fioi_melth(n) - strm_qbot(n) - (Faxa_snow(n)+Foxx_rofi(n))*latice) * dt/(cpsw*rhosw* strm_h(n))

             ! compute ice formed or melt potential
             Fioo_q(n) = (tfreeze(n) - So_t(n))*(cpsw*rhosw*strm_h(n))/dt ! ice formed q>0
             
             ! reset temp
             if (reset_temp) then
                So_t(n)  = max(tfreeze(n),So_t(n))
             end if

             ! save somtp to restart file
             somtp(n) = So_t(n)
          endif
       end do
       deallocate(tfreeze)

    endif   ! first_time

    first_time = .false.

  end subroutine docn_datamode_som_advance

  !===============================================================================
  subroutine docn_datamode_som_restart_write(case_name, inst_suffix, ymd, tod, &
       logunit, mpicom, my_task, sdat)
    
    ! write restart file

    ! input/output variables
    character(len=*)            , intent(in)    :: case_name
    character(len=*)            , intent(in)    :: inst_suffix
    integer                     , intent(in)    :: ymd       ! model date
    integer                     , intent(in)    :: tod       ! model sec into model date
    integer                     , intent(in)    :: logunit
    integer                     , intent(in)    :: my_task
    integer                     , intent(in)    :: mpicom
    type(shr_strdata_type)      , intent(inout) :: sdat
    !-------------------------------------------------------------------------------

    call dshr_restart_write(rpfile, case_name, 'docn', inst_suffix, ymd, tod, &
         logunit, mpicom, my_task, sdat, fld=somtp, fldname='somtp')

  end subroutine docn_datamode_som_restart_write

  !===============================================================================
  subroutine docn_datamode_som_restart_read(rest_filem, inst_suffix, logunit, my_task, mpicom, sdat)

    ! read restart file

    ! input/output arguments
    character(len=*)            , intent(inout) :: rest_filem
    character(len=*)            , intent(in)    :: inst_suffix
    integer                     , intent(in)    :: logunit
    integer                     , intent(in)    :: my_task
    integer                     , intent(in)    :: mpicom
    type(shr_strdata_type)      , intent(inout) :: sdat
    !-------------------------------------------------------------------------------

    ! allocate module memory for restart fields that are read in
    allocate(somtp(sdat%model_lsize))

    ! read restart
    call dshr_restart_read(rest_filem, rpfile, inst_suffix, nullstr, logunit, my_task, mpicom, sdat, &
         fld=somtp, fldname='somtp')

  end subroutine docn_datamode_som_restart_read

end module docn_datamode_som_mod
