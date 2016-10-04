module component_mod

  !----------------------------------------------------------------------------
  ! share code & libs
  !----------------------------------------------------------------------------
  use shr_kind_mod,     only: r8 => SHR_KIND_R8 
  use shr_kind_mod,     only: cs => SHR_KIND_CS
  use shr_kind_mod,     only: cl => SHR_KIND_CL
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  use shr_const_mod,    only: shr_const_cday
  use shr_file_mod,     only: shr_file_setLogLevel, shr_file_setLogUnit
  use shr_file_mod,     only: shr_file_setIO, shr_file_getUnit
  use shr_scam_mod,     only: shr_scam_checkSurface
  use shr_mpi_mod,      only: shr_mpi_min, shr_mpi_max
  use shr_mem_mod,      only: shr_mem_init, shr_mem_getusage
  use shr_cal_mod,      only: shr_cal_date2ymd
  use shr_orb_mod,      only: shr_orb_params
  use shr_reprosum_mod, only: shr_reprosum_setopts
  use seq_comm_mct,     only: GLOID, CPLID, logunit
  use seq_comm_mct,     only: seq_comm_iamin, seq_comm_namelen, num_inst_frc
  use seq_comm_mct,     only: seq_comm_suffix, seq_comm_name, seq_comm_setnthreads
  use seq_comm_mct,     only: seq_comm_getinfo => seq_comm_setptrs
  use seq_comm_mct,     only: seq_comm_petlist 
  use seq_infodata_mod, only: seq_infodata_putData, seq_infodata_GetData
  use seq_infodata_mod, only: seq_infodata_exchange, seq_infodata_type
  use seq_diag_mct,     only: seq_diag_avect_mct 
  use seq_map_type_mod  
  use seq_map_mod
  use t_drv_timers_mod
  use component_type_mod
  use seq_cdata_mod,    only : seq_cdata
  use mct_mod   ! mct_ wrappers for mct lib
  use perf_mod
  use ESMF
#ifdef ESMF_INTERFACE
  use esmfshr_mod
#endif 

  implicit none

#include <mpif.h>

  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: component_init_pre
  public :: component_init_cc            ! mct and esmf versions
  public :: component_init_cx
  public :: component_init_aream
  public :: component_init_areacor
#ifdef ESMF_INTERFACE
  public :: component_init_update_petlist !esmf only
#endif
  public :: component_run                 ! mct and esmf versions
  public :: component_final               ! mct and esmf versions 
  public :: component_exch
  public :: component_diag


  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

   logical  :: iamroot_GLOID, iamroot_CPLID         ! GLOID, CPLID masterproc
   logical  :: iamin_CPLID                          ! true => pe associated with CPLID
   integer  :: mpicom_GLOID, mpicom_CPLID           ! GLOID, CPLID mpi communicator 
   integer  :: nthreads_GLOID, nthreads_CPLID
   logical  :: drv_threading

  !===============================================================================

contains

  !===============================================================================

  subroutine component_init_pre(comp, compid, cplcompid, cplallcompid, &
       infodata, ntype)

    !---------------------------------------------------------------
    ! Initialize driver rearrangers and AVs on driver
    ! Initialize cdata_*x data
    ! Zero out x2*_** in case it never gets used then it'll produce zeros in diags
    ! For ensembles, create only a single dom_*x for the coupler based on the
    !   first ensemble member.  otherwise, just extend the dom_** and dom_*x to
    !   other ensemble members.
    !
    ! Arguments
    type(component_type)     , intent(inout)         :: comp(:)
    integer                  , intent(in)            :: compid(:)
    integer                  , intent(in)            :: cplcompid(:)
    integer                  , intent(in)            :: cplallcompid
    type (seq_infodata_type) , intent(inout), target :: infodata
    character(len=3)         , intent(in)            :: ntype
    !
    ! Local Variables
    integer  :: eci       ! index
    character(*), parameter :: subname = '(component_init_pre)'
    !---------------------------------------------------------------

    ! initialize module variables (this is repetitive here- but does not require a different routine)

    call seq_infodata_getdata(infodata, drv_threading=drv_threading)
    call seq_comm_getinfo(GLOID, mpicom=mpicom_GLOID, iamroot=iamroot_GLOID, nthreads=nthreads_GLOID)
    call seq_comm_getinfo(CPLID, mpicom=mpicom_CPLID, iamroot=iamroot_CPLID, nthreads=nthreads_CPLID)
    iamin_CPLID = seq_comm_iamin(CPLID)
   
    ! Initialize component type variables
    do eci = 1,size(comp)

       comp(eci)%compid       = compid(eci)
       comp(eci)%cplcompid    = cplcompid(eci)
       comp(eci)%cplallcompid = cplallcompid

       call seq_comm_getinfo(comp(eci)%cplallcompid, mpicom=comp(eci)%mpicom_cplallcompid)
       call seq_comm_getinfo(comp(eci)%cplcompid   , mpicom=comp(eci)%mpicom_cplcompid)
       call seq_comm_getinfo(comp(eci)%compid      , mpicom=comp(eci)%mpicom_compid)
       call seq_comm_getinfo(comp(eci)%compid      , iamroot=comp(eci)%iamroot_compid)
       call seq_comm_getinfo(comp(eci)%compid      , nthreads=comp(eci)%nthreads_compid)

       comp(eci)%iamin_compid       =  seq_comm_iamin (comp(eci)%compid)
       comp(eci)%iamin_cplcompid    =  seq_comm_iamin (comp(eci)%cplcompid)
       comp(eci)%iamin_cplallcompid =  seq_comm_iamin (comp(eci)%cplallcompid)
       comp(eci)%suffix             =  seq_comm_suffix(comp(eci)%compid)
       comp(eci)%name               =  seq_comm_name  (comp(eci)%compid)
       comp(eci)%ntype              =  ntype(1:3)
       comp(eci)%oneletterid        =  ntype(1:1)

       if (eci == 1) then
          allocate(comp(1)%dom_cx)
          allocate(comp(1)%gsmap_cx)
       else
          comp(eci)%dom_cx   => comp(1)%dom_cx 
          comp(eci)%gsmap_cx => comp(1)%gsmap_cx 
       end if

       ! Set cdata_cc - unique for each instance
       allocate(comp(eci)%dom_cc)
       allocate(comp(eci)%gsmap_cc)
       allocate(comp(eci)%cdata_cc)
       comp(eci)%cdata_cc%name     = 'cdata_'//ntype(1:1)//ntype(1:1)
       comp(eci)%cdata_cc%ID       =  comp(eci)%compid
       comp(eci)%cdata_cc%mpicom   =  comp(eci)%mpicom_compid
       comp(eci)%cdata_cc%dom      => comp(eci)%dom_cc 
       comp(eci)%cdata_cc%gsmap    => comp(eci)%gsmap_cc
       comp(eci)%cdata_cc%infodata => infodata

       ! Determine initial value of comp_present in infodata - to do - add this to component 

       if (comp(1)%oneletterid == 'a') call seq_infodata_getData(infodata, atm_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'l') call seq_infodata_getData(infodata, lnd_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'i') call seq_infodata_getData(infodata, ice_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'o') call seq_infodata_getData(infodata, ocn_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'r') call seq_infodata_getData(infodata, rof_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'g') call seq_infodata_getData(infodata, glc_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'w') call seq_infodata_getData(infodata, wav_present=comp(eci)%present)

    end do

  end subroutine component_init_pre

  !===============================================================================

#ifndef ESMF_INTERFACE
  subroutine component_init_cc(Eclock, comp, comp_init, infodata, NLFilename, &
       seq_flds_x2c_fluxes, seq_flds_c2x_fluxes)

    !---------------------------------------------------------------
    !
    ! Arguments
    type(ESMF_Clock)         , intent(inout) :: EClock
    type(component_type)     , intent(inout) :: comp(:)
    interface 
       subroutine comp_init( Eclock, cdata, x2c, c2x, nlfilename)
         use ESMF         , only: ESMF_Clock
         use seq_cdata_mod, only: seq_cdata
         use mct_mod      , only: mct_avect
         implicit none
         type(ESMF_Clock), intent(inout) :: EClock
         type(seq_cdata) , intent(inout) :: cdata
         type(mct_aVect) , intent(inout) :: x2c
         type(mct_aVect) , intent(inout) :: c2x   
         character(len=*), optional, intent(IN) :: NLFilename ! Namelist filename
       end subroutine comp_init
    end interface 
    type (seq_infodata_type) , intent(inout)        :: infodata
    character(len=*)         , intent(in)           :: NLFilename 
    character(len=*)         , intent(in), optional :: seq_flds_x2c_fluxes
    character(len=*)         , intent(in), optional :: seq_flds_c2x_fluxes
    !
    ! Local Variables
    integer :: k1, k2
    integer :: eci
    character(*), parameter :: subname = '(component_init_cc:mct)'
    character(*), parameter :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    ! **** Initialize component - this initializes  x2c_cc and c2x_cc ***
    ! the following will call the appropriate comp_init_mct routine

    if (comp(1)%iamin_cplallcompid) then
       call seq_infodata_exchange(infodata, comp(1)%cplallcompid, &
            'cpl2'//comp(1)%ntype(1:3)//'_init')
    end if

    ! The following initializes the component instance cdata_cc (gsmap and dom), 
    ! x2c_cc and c2x_cc

    do eci = 1,size(comp)
       if (iamroot_CPLID .and. comp(eci)%present) then
          write(logunit,F00) 'Initialize component '//trim(comp(eci)%ntype)
          call shr_sys_flush(logunit)
       endif

       if (.not. associated(comp(eci)%x2c_cc)) allocate(comp(eci)%x2c_cc)
       if (.not. associated(comp(eci)%c2x_cc)) allocate(comp(eci)%c2x_cc)

       if (comp(eci)%iamin_compid .and. comp(eci)%present) then
          if (drv_threading) call seq_comm_setnthreads(comp(eci)%nthreads_compid)
          call shr_sys_flush(logunit)
          
          if (present(seq_flds_x2c_fluxes)) then
             call mct_avect_vecmult(comp(eci)%x2c_cc, comp(eci)%drv2mdl, seq_flds_x2c_fluxes, mask_spval=.true.)
          end if

          call t_set_prefixf(comp(1)%oneletterid//"_i:")
          call comp_init( EClock, comp(eci)%cdata_cc, comp(eci)%x2c_cc, comp(eci)%c2x_cc, &
               NLFilename=NLFilename )
          call t_unset_prefixf()
          
          if (present(seq_flds_c2x_fluxes)) then
             call mct_avect_vecmult(comp(eci)%c2x_cc, comp(eci)%mdl2drv, seq_flds_c2x_fluxes, mask_spval=.true.)
          end if
          
          if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
       end if
    end do

    if (comp(1)%iamin_cplcompid) then
       call seq_infodata_exchange(infodata, comp(1)%cplcompid, &
            comp(1)%ntype(1:3)//'2cpl_init')
    endif

    ! Determine final value of comp_present in infodata (after component initialization)

    do eci = 1,size(comp) 
       if (comp(1)%oneletterid == 'a') call seq_infodata_getData(infodata, atm_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'l') call seq_infodata_getData(infodata, lnd_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'i') call seq_infodata_getData(infodata, ice_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'o') call seq_infodata_getData(infodata, ocn_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'r') call seq_infodata_getData(infodata, rof_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'g') call seq_infodata_getData(infodata, glc_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'w') call seq_infodata_getData(infodata, wav_present=comp(eci)%present)
    end do


    ! Initialize aream, set it to area for now until maps are read
    !   in some cases, maps are not read at all !!
    ! Entire domain must have reasonable values before calling xxx2xxx init

    do eci = 1,size(comp)
       if (comp(eci)%iamin_compid .and. comp(eci)%present) then
          if (drv_threading) call seq_comm_setnthreads(comp(eci)%nthreads_compid)
          k1 = mct_aVect_indexRa(comp(eci)%cdata_cc%dom%data, "area"  ,perrWith='aa area ')
          k2 = mct_aVect_indexRa(comp(eci)%cdata_cc%dom%data, "aream" ,perrWith='aa aream')

          comp(eci)%cdata_cc%dom%data%rAttr(k2,:) = comp(eci)%cdata_cc%dom%data%rAttr(k1,:)

          if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
       endif
    end do

  end subroutine component_init_cc
#endif

  !===============================================================================

#ifdef ESMF_INTERFACE
  subroutine component_init_cc(Eclock, drvcomp, comp, gridcomp_register,  &
       infodata, NlFilename, seq_flds_x2c_fields, seq_flds_c2x_fields, &
       seq_flds_x2c_fluxes, seq_flds_c2x_fluxes) 

    !---------------------------------------------------------------
    ! Uses
    use esmf2mct_mod,         only: esmf2mct_init, esmf2mct_copy
    use mct2esmf_mod,         only: mct2esmf_init
    use seq_flds_mod,         only: seq_flds_dom_coord, seq_flds_dom_other, seq_flds_dom_fields
    !
    ! Arguments
    type(ESMF_Clock)     , intent(inout) :: EClock
    type(ESMF_CplComp)   , intent(inout) :: drvComp
    type(component_type) , intent(inout) :: comp(:)
    interface 
       subroutine gridcomp_register(gridcomp, rc)
         use ESMF
         implicit none
         type(ESMF_GridComp)  :: gridcomp
         integer, intent(out) :: rc
       end subroutine gridcomp_register
    end interface 
    type (seq_infodata_type) , intent(inout)        :: infodata
    character(len=*)         , intent(in), optional :: NLFilename 
    character(len=*)         , intent(in), optional :: seq_flds_x2c_fields
    character(len=*)         , intent(in), optional :: seq_flds_c2x_fields
    character(len=*)         , intent(in), optional :: seq_flds_x2c_fluxes
    character(len=*)         , intent(in), optional :: seq_flds_c2x_fluxes
    !
    ! Local Variables
    integer             :: k1, k2
    integer             :: nx, ny
    integer             :: eci
    integer             :: rc, urc
    integer             :: init_phase
    integer             :: lsize
    type(ESMF_DistGrid) :: distgrid_cc
    type(ESMF_DistGrid) :: distgrid_cx
    type(ESMF_Array)    :: dom_cc_array
    type(ESMF_Array)    :: dom_cx_array
    type(ESMF_Array)    :: x2c_cc_array
    type(ESMF_Array)    :: x2c_cx_array
    type(ESMF_Array)    :: c2x_cc_array
    type(ESMF_Array)    :: c2x_cx_array
    integer , pointer   :: petlist(:)
    real(R8), pointer   :: fptr(:,:)            ! pointer into    array data
    character(len=1)    :: cid
    character(len=8196) :: mct_names_x2c, mct_names_c2x, mct_names_dom
    character(*), parameter :: subname = '(component_init_cc:esmf)'
    character(*), parameter :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    if (present(seq_flds_x2c_fluxes) .and. present(seq_flds_c2x_fluxes)) then
       init_phase = 2
    else
       init_phase = 1
    end if

    do eci = 1,size(comp)

       if (init_phase == 1) then

          ! Create gridcomp for this instance

          call seq_comm_petlist(comp(eci)%compid, petlist) 

          comp(eci)%gridcomp_cc = ESMF_GridCompCreate(name=trim(comp(eci)%name), petList=petlist, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          call ESMF_GridCompSetServices(comp(eci)%gridcomp_cc, userRoutine=gridcomp_register, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          ! create import and export states state_x2c_cc and state_c2x_cc

          comp(eci)%x2c_cc_state = ESMF_StateCreate(name=trim(comp(eci)%ntype)//" x2c_cc", &
               stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          comp(eci)%c2x_cc_state = ESMF_StateCreate(name=trim(comp(eci)%ntype)//" c2x_cc", &
               stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          ! link attributes

          call ESMF_AttributeLink(drvcomp, comp(eci)%gridcomp_cc, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

       end if

       !--------------------------------------------------
       ! COMPONENT-CPL PES
       !--------------------------------------------------
       ! Transfer infodata from coupler pes -> component pes

       if (comp(1)%iamin_cplallcompid) then
          call seq_infodata_exchange(infodata, comp(1)%cplallcompid, &
               'cpl2'//comp(1)%ntype(1:3)//'_init')
       end if

       !--------------------------------------------------
       ! COMPONENT PES
       !--------------------------------------------------

       ! The following initializes the component instance values of x2c_cc and c2x_cc
       
       if (iamroot_CPLID .and. comp(eci)%present) then
          write(logunit,F00) 'Initialize component '//trim(comp(eci)%ntype)
          call shr_sys_flush(logunit)
       end if

       if (comp(eci)%iamin_compid) then

          if (comp(eci)%present) then

             if (drv_threading) call seq_comm_setnthreads(comp(eci)%nthreads_compid)
             call shr_sys_flush(logunit)

             if (init_phase == 1) then 
                call ESMF_AttributeSet(comp(eci)%c2x_cc_state, name="ID", &
                     value=comp(eci)%compid, rc=rc)
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
             end if

             ! Set init_phase attribute value
             call ESMF_AttributeSet(comp(eci)%c2x_cc_state, name=trim(comp(eci)%ntype)//"_phase", &
                  value=init_phase, rc=rc)
             if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

             ! Convert infodata information into appropriate export state attributes
             call esmfshr_infodata_infodata2state(infodata, comp(eci)%c2x_cc_state, &
                  id=comp(eci)%compid, rc=rc)
             if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

             ! Rescale the attribute vector before sending for phase 2 of atm initialization
             ! Note that x2c_cc attribute vector and x2c_cc_array share the same memory 
             ! (set in init_phase=1)

             if (init_phase == 2) then ! phase 2 (only for atm for now)
                call ESMF_StateGet(comp(eci)%x2c_cc_state, itemName="x2d", array=x2c_cc_array, rc=rc)
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call mct_avect_vecmult(comp(eci)%x2c_cc, comp(eci)%drv2mdl, seq_flds_x2c_fluxes, mask_spval=.true.)
             end if

             !-----------------------------------
             ! *** call into ESMF init method ***
             !-----------------------------------
             call t_set_prefixf(comp(1)%oneletterid//"_i:")
             call ESMF_GridCompInitialize(comp(eci)%gridcomp_cc, &
                  importState=comp(eci)%x2c_cc_state, exportState=comp(eci)%c2x_cc_state, &
                  clock=EClock, userRc=urc, rc=rc)
             if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
             if (rc  /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc , endflag=ESMF_END_ABORT)
             call t_unset_prefixf()
             !-----------------------------------

             if (init_phase == 2) then
                ! Rescale the attribute vector after receiving for phase 2 of atm initialization
                call ESMF_StateGet(comp(eci)%c2x_cc_state, itemName="d2x", array=c2x_cc_array, rc=rc)
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call mct_avect_vecmult(comp(eci)%c2x_cc, comp(eci)%mdl2drv, seq_flds_c2x_fluxes, mask_spval=.true.)
             end if

             ! Convert appropriate export state attributes back to infodata, 
             ! the new nextsw_cday is updated in infodata
             call esmfshr_infodata_state2infodata(comp(eci)%c2x_cc_state, infodata, rc=rc)
             if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

             ! The following is only from infodata updates on the component pes

             cid = comp(1)%oneletterid
             if (cid=='a') call seq_infodata_getData(infodata, atm_present=comp(eci)%present, atm_nx=nx, atm_ny=ny)
             if (cid=='l') call seq_infodata_getData(infodata, lnd_present=comp(eci)%present, lnd_nx=nx, lnd_ny=ny)
             if (cid=='i') call seq_infodata_getData(infodata, ice_present=comp(eci)%present, ice_nx=nx, ice_ny=ny)
             if (cid=='o') call seq_infodata_getData(infodata, ocn_present=comp(eci)%present, ocn_nx=nx, ocn_ny=ny)
             if (cid=='r') call seq_infodata_getData(infodata, rof_present=comp(eci)%present, rof_nx=nx, rof_ny=ny)
             if (cid=='g') call seq_infodata_getData(infodata, glc_present=comp(eci)%present, glc_nx=nx, glc_ny=ny)
             if (cid=='w') call seq_infodata_getData(infodata, wav_present=comp(eci)%present, wav_nx=nx, wav_ny=ny)

             if (init_phase == 1 .and. comp(eci)%present) then

                ! allocate memory and initialize dom_cc%data, x2c_cc and c2x_cc attribute vectors

                if (.not. associated(comp(eci)%x2c_cc)) allocate(comp(eci)%x2c_cc)
                if (.not. associated(comp(eci)%c2x_cc)) allocate(comp(eci)%c2x_cc)
                if (.not. associated(comp(eci)%dom_cc)) allocate(comp(eci)%dom_cc)

                call ESMF_StateGet(comp(eci)%x2c_cc_state, itemName="x2d", array=x2c_cc_array, rc=rc)
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call ESMF_StateGet(comp(eci)%c2x_cc_state, itemName="d2x", array=c2x_cc_array, rc=rc)
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call ESMF_StateGet(comp(eci)%c2x_cc_state, itemName="domain", array=dom_cc_array, rc=rc)
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                ! initialize MCT gsmap_cc global seg map from ESMF distgrid_cc

                call ESMF_ArrayGet(c2x_cc_array, distgrid=distgrid_cc, rc=rc) 
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call esmf2mct_init(distgrid_cc, comp(eci)%compid, comp(eci)%gsmap_cc, &
                     comp(eci)%mpicom_compid, nx*ny, rc)
                if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                ! initialize MCT x2c_cc and c2x_cc attribute vectors from ESMF arrays

                call esmf2mct_init(x2c_cc_array, comp(eci)%x2c_cc, rc)
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
                call esmf2mct_copy(x2c_cc_array, comp(eci)%x2c_cc, rc=rc)
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call esmf2mct_init(c2x_cc_array, comp(eci)%c2x_cc, rc)
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
                call esmf2mct_copy(c2x_cc_array, comp(eci)%c2x_cc, rc=rc)
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                ! initialize MCT ggrid dom_cc from ESMF array

                call ESMF_DistGridGet(distgrid_cc, localDe=0, elementCount=lsize, rc=rc)
                if(rc /= 0) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call mct_gGrid_init(comp(eci)%dom_cc, coordchars=seq_flds_dom_coord, &
                     otherchars=seq_flds_dom_other, lsize=lsize )

                call esmf2mct_copy(dom_cc_array, comp(eci)%dom_cc%data, rc=rc)
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                ! destroy original ESMF arrays, x2c_cc_array, c2x_cc_array and c2x_dom_array 
                ! in preparation for creating new ones that share memory with MCT attribute vecs
                ! *** But this will remove any attributes that were originally in the ESMF array
                ! so need to extract this info out first ****

                call ESMF_AttributeGet(x2c_cc_array, name='mct_names', value=mct_names_x2c, rc=rc)
                if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call ESMF_ArrayDestroy(x2c_cc_array, rc=rc) ! destroy the Array
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call ESMF_AttributeGet(c2x_cc_array, name='mct_names', value=mct_names_c2x, rc=rc)
                if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call ESMF_ArrayDestroy(c2x_cc_array, rc=rc) ! destroy the Array
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call ESMF_AttributeGet(dom_cc_array, name='mct_names', value=mct_names_dom, rc=rc)
                if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call ESMF_ArrayDestroy(dom_cc_array, rc=rc) ! destroy the Array
                if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                ! create new ESMF arrays
                ! - x2c_cc_array  shares memory with comp(eci)%x2c_cc attribute vector 
                ! - c2x_cc_array  shares memory with comp(eci)%c2x_cc attribute vector 
                ! - c2x_dom_array shares memory with comp(eci)%dom_cc%data attribute vector 

                x2c_cc_array = ESMF_ArrayCreate(distgrid=distgrid_cc, farrayPtr=comp(eci)%x2c_cc%rattr, &
                     distgridToArrayMap=(/2/), name="x2d", rc=rc)
                if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                c2x_cc_array = ESMF_ArrayCreate(distgrid=distgrid_cc, farrayPtr=comp(eci)%c2x_cc%rattr, &
                     distgridToArrayMap=(/2/), name="d2x", rc=rc)
                if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                dom_cc_array = ESMF_ArrayCreate(distgrid=distgrid_cc, farrayPtr=comp(eci)%dom_cc%data%rattr, &
                     distgridToArrayMap=(/2/), name="domain", rc=rc)
                if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call ESMF_AttributeSet(x2c_cc_array, name="mct_names", value=trim(mct_names_x2c), rc=rc)
                if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call ESMF_AttributeSet(c2x_cc_array, name="mct_names", value=trim(mct_names_c2x), rc=rc)
                if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call ESMF_AttributeSet(dom_cc_array, name="mct_names", value=trim(mct_names_dom), rc=rc)
                if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

                call ESMF_StateReplace(comp(eci)%x2c_cc_state, (/x2c_cc_array/), rc=rc)
                call ESMF_StateReplace(comp(eci)%c2x_cc_state, (/c2x_cc_array/), rc=rc)
                call ESMF_StateReplace(comp(eci)%c2x_cc_state, (/dom_cc_array/), rc=rc)

             end if

             ! Convert appropriate export state attributes back to infodata, 
             ! the new nextsw_cday is updated in infodata
             call esmfshr_infodata_state2infodata(comp(eci)%c2x_cc_state, infodata, rc=rc)
             if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

             if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)

          end if   ! end of comp(eci)%present
       end if   ! end of comp(eci)%iamin_compid

       ! allocate memory for attribute vectors that are in cpl id - if compid and cplid
       ! are not the smae
       if (comp(eci)%iamin_cplcompid) then
          if (init_phase == 1 .and. comp(eci)%present) then
             if (.not. associated(comp(eci)%x2c_cc)) allocate(comp(eci)%x2c_cc)
             if (.not. associated(comp(eci)%c2x_cc)) allocate(comp(eci)%c2x_cc)
             if (.not. associated(comp(eci)%dom_cc)) allocate(comp(eci)%dom_cc)
          end if
       end if

    end do   ! end of loop over instances

    !--------------------------------------------------
    ! COMPONENT -CPL PES
    !--------------------------------------------------
    ! Transfer infodata between component pes -> coupler pes

    if (comp(1)%iamin_cplcompid) then
       call seq_infodata_exchange(infodata, comp(1)%cplcompid, &
            comp(1)%ntype(1:3)//'2cpl_init')
    endif

    ! Determine final value of comp_present in infodata (after component initialization)

    do eci = 1,size(comp) 
       if (comp(1)%oneletterid == 'a') call seq_infodata_getData(infodata, atm_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'l') call seq_infodata_getData(infodata, lnd_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'i') call seq_infodata_getData(infodata, ice_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'o') call seq_infodata_getData(infodata, ocn_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'r') call seq_infodata_getData(infodata, rof_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'g') call seq_infodata_getData(infodata, glc_present=comp(eci)%present)
       if (comp(1)%oneletterid == 'w') call seq_infodata_getData(infodata, wav_present=comp(eci)%present)
    end do

    !--------------------------------------------------
    ! COMPONENT PES
    !--------------------------------------------------
    ! Initialize aream, set it to area for now until maps are read
    !   in some cases, maps are not read at all !!
    ! Entire domain must have reasonable values before calling xxx2xxx init

    do eci = 1,size(comp)
       if (comp(eci)%iamin_compid .and. comp(eci)%present) then
          if (drv_threading) call seq_comm_setnthreads(comp(eci)%nthreads_compid)
          k1 = mct_aVect_indexRa(comp(eci)%cdata_cc%dom%data, "area"  ,perrWith='aa area ')
          k2 = mct_aVect_indexRa(comp(eci)%cdata_cc%dom%data, "aream" ,perrWith='aa aream')

          comp(eci)%cdata_cc%dom%data%rAttr(k2,:) = comp(eci)%cdata_cc%dom%data%rAttr(k1,:)

          if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
       endif
    end do

  end subroutine component_init_cc
#endif

  !===============================================================================

#ifdef ESMF_INTERFACE
  subroutine component_init_update_petlist(comp, vm)

    !---------------------------------------------------------------
    ! Arguments
    type(component_type), intent(inout) :: comp(:) 
    type(ESMF_VM) :: vm
    !
    ! Local variables
    integer :: eci
    integer :: rc
    integer, pointer :: petlist(:)
    character(*), parameter :: subname = '(component_init_update_petlist)'
    !---------------------------------------------------------------

   ! Update petlist attribute
    do eci = 1, size(comp)
!BUG       write(6,*)'DEBUG: comp is ',trim(comp(eci)%name)
       call seq_comm_petlist(comp(eci)%compid, petlist)
!BUG       write(6,*)'DEBUG: petlist is ',petlist
!BUG       call ESMF_AttributeUpdate(comp(eci)%gridcomp_cc, vm, rootList=petlist, rc=rc)
!BUG       if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
!BUG       write(6,*)'DEBUG: successfully updated attributes'
    enddo

  end subroutine component_init_update_petlist
#endif

  !===============================================================================
    
  subroutine component_init_cx(comp, infodata) 

    !---------------------------------------------------------------
    ! Uses
    use cplcomp_exchange_mod, only: seq_mctext_gsmapinit, seq_mctext_avInit
    use cplcomp_exchange_mod, only: seq_mctext_avExtend, seq_mctext_gGridInit
    use cplcomp_exchange_mod, only: seq_map_init_exchange, seq_map_map_exchange
    use seq_domain_mct,       only: seq_domain_compare
    use mct_mod,              only: mct_ggrid_clean
    !
    ! Arguments
    type(component_type)     , intent(inout) :: comp(:)
    type (seq_infodata_type) , intent(inout) :: infodata
    !
    ! Local Variables
    integer         :: eci
    integer         :: rc        ! return code
    type(mct_gGrid) :: dom_tmp   ! temporary
    character(*), parameter :: subname = '(component_init_cx)'
    character(*), parameter :: F0I = "('"//subname//" : ', A, 2i8 )"
    !---------------------------------------------------------------

    ! Initialize driver rearrangers and AVs on driver
    ! Initialize cdata_*x data
    ! Zero out x2*_** in case it never gets used then it'll produce zeros in diags
    ! For ensembles, create only a single dom_*x for the coupler based on the
    !   first ensemble member.  otherwise, just extend the dom_** and dom_*x to
    !   other ensemble members.

    do eci = 1,size(comp)
       if (comp(eci)%present) then

          if (iamroot_CPLID) then
             write(logunit,*) ' '
             call shr_sys_flush(logunit)
          end if

          if (comp(eci)%iamin_cplcompid) then

             ! Create gsmap_cx (note that comp(eci)%gsmap_cx all point to comp(1)%gsmap_cx
             ! This will only be valid on the coupler pes
             if (eci == 1) then
                if (iamroot_CPLID) then
                   write(logunit,F0I) 'creating gsmap_cx for '//comp(eci)%ntype(1:3)
                   call shr_sys_flush(logunit)
                end if
                call seq_mctext_gsmapInit(comp(1))
             endif

             ! Create mapper_Cc2x and mapper_Cx2c
             allocate(comp(eci)%mapper_Cc2x, comp(eci)%mapper_Cx2c)  
             if (iamroot_CPLID) then
                write(logunit,F0I) 'Initializing mapper_C'//comp(eci)%ntype(1:1)//'2x',eci
                call shr_sys_flush(logunit)
             end if
             call seq_map_init_exchange(comp(eci), flow='c2x', mapper=comp(eci)%mapper_Cc2x)    
             if (iamroot_CPLID) then
                write(logunit,F0I) 'Initializing mapper_Cx2'//comp(eci)%ntype(1:1),eci
                call shr_sys_flush(logunit)
             end if
             call seq_map_init_exchange(comp(eci), flow='x2c', mapper=comp(eci)%mapper_Cx2c)  

             ! Create x2c_cx and c2x_cx
             allocate(comp(eci)%x2c_cx, comp(eci)%c2x_cx)
             call seq_mctext_avinit(comp(eci), flow='x2c')
             call seq_mctext_avinit(comp(eci), flow='c2x')

             ! Create dom_cx (note that  comp(eci)%dom_cx all point to  comp(1)%dom_cx
             ! Then verify other ensembles have same domain by comparing to dom_cx
             if (eci == 1) then  ! create dom_cx
                if (iamroot_CPLID) then
                   write(logunit,F0I) 'creating dom_cx'
                   call shr_sys_flush(logunit)
                end if
                call seq_mctext_gGridInit(comp(1))
                call seq_map_map_exchange(comp(1), flow='c2x', dom_flag=.true., msgtag=comp(1)%cplcompid*100+1*10+1)
             else if (eci > 1) then  
                if (iamroot_CPLID) then
                   write(logunit,F0I) 'comparing comp domain ensemble number ',eci
                   call shr_sys_flush(logunit)
                end if
                call seq_mctext_avExtend(comp(eci)%dom_cx%data, cplid, comp(eci)%cplcompid)
                call seq_mctext_gGridInit(comp(eci), dom_tmp)
                call seq_map_map_exchange(comp(eci), flow='c2x', dom_flag=.true., dom_tmp=dom_tmp)
                if (iamin_CPLID) then
                   call seq_domain_compare(comp(eci)%dom_cx, dom_tmp, mpicom_CPLID)
                end if
                call mct_ggrid_clean(dom_tmp,rc)
             endif

             call mct_avect_zero(comp(eci)%x2c_cc)
             call mct_avect_zero(comp(eci)%x2c_cx)

          end if ! if comp(eci)%iamin_cplcompid
       end if  ! if comp(eci)%present
    end do  ! end of eci loop

  end subroutine component_init_cx

  !===============================================================================

  subroutine component_init_aream(infodata, rof_c2_ocn, samegrid_ao, samegrid_al, &
       samegrid_ro, samegrid_lg)

    !---------------------------------------------------------------
    ! Description
    ! Update (read) aream in domains where appropriate - ON cpl pes
    !
    ! Uses
    use prep_ocn_mod,       only : prep_ocn_get_mapper_Fa2o
    use prep_lnd_mod,       only : prep_lnd_get_mapper_Sa2l
    use prep_ice_mod,       only : prep_ice_get_mapper_SFo2i
    use prep_glc_mod,       only : prep_glc_get_mapper_Sl2g
    use component_type_mod, only : atm, lnd, ice, ocn, rof, glc
    !
    ! Arguments
    type (seq_infodata_type) , intent(inout) :: infodata
    logical                  , intent(in)    :: rof_c2_ocn
    logical                  , intent(in)    :: samegrid_ao
    logical                  , intent(in)    :: samegrid_al
    logical                  , intent(in)    :: samegrid_ro
    logical                  , intent(in)    :: samegrid_lg  ! lnd & glc on same grid
    !
    ! Local variables
    type(mct_gsmap), pointer :: gsmap_s, gsmap_d
    type(mct_ggrid), pointer :: dom_s, dom_d
    type(seq_map)  , pointer :: mapper_Fa2o
    type(seq_map)  , pointer :: mapper_Sa2l
    type(seq_map)  , pointer :: mapper_SFo2i
    type(seq_map)  , pointer :: mapper_Sl2g
    logical                  :: atm_present ! atm present flag
    logical                  :: lnd_present ! lnd present flag
    logical                  :: ocn_present ! ocn present flag
    logical                  :: ice_present ! ice present flag
    logical                  :: glc_present ! glc present flag
    integer                  :: ka,km
    character(*), parameter :: subname = '(component_init_aream)'
    !---------------------------------------------------------------

    ! Note that the following is assumed to hold - all gsmaps_cx for a given 
    ! instance of a component (e.g. atm(i)) are identical on the coupler processes

    mapper_Fa2o  => prep_ocn_get_mapper_Fa2o()
    mapper_Sa2l  => prep_lnd_get_mapper_Sa2l()
    mapper_SFo2i => prep_ice_get_mapper_SFo2i()
    mapper_Sl2g  => prep_glc_get_mapper_Sl2g()

    call seq_infodata_GetData( infodata, &
         atm_present=atm_present,        &
         ocn_present=ocn_present,        &
         ice_present=ice_present,        &
         lnd_present=lnd_present,        &
         glc_present=glc_present)

    if (atm_present .and. ocn_present) then
       if (samegrid_ao) then
          dom_s  => component_get_dom_cx(atm(1))   !dom_ax
          dom_d  => component_get_dom_cx(ocn(1))   !dom_ox
          ka = mct_aVect_indexRa(dom_s%data, "area" )
          km = mct_aVect_indexRa(dom_s%data, "aream" )
          dom_s%data%rAttr(km,:) = dom_s%data%rAttr(ka,:)

          call seq_map_map(mapper_Fa2o, av_s=dom_s%data, av_d=dom_d%data, fldlist='aream')
       else
          gsmap_s => component_get_gsmap_cx(ocn(1)) ! gsmap_ox
          gsmap_d => component_get_gsmap_cx(atm(1)) ! gsmap_ax
          dom_s   => component_get_dom_cx(ocn(1))   ! dom_ox
          dom_d   => component_get_dom_cx(atm(1))   ! dom_ax

          call seq_map_readdata('seq_maps.rc','ocn2atm_fmapname:', mpicom_CPLID, CPLID, &
               gsmap_s=gsmap_s, av_s=dom_s%data, avfld_s='aream', filefld_s='area_a', &
               gsmap_d=gsmap_d, av_d=dom_d%data, avfld_d='aream', filefld_d='area_b', &
               string='ocn2atm aream initialization')
       endif
    end if

    if (ice_present .and. ocn_present) then
       dom_s  => component_get_dom_cx(ocn(1))   !dom_ox
       dom_d  => component_get_dom_cx(ice(1))   !dom_ix

       call seq_map_map(mapper_SFo2i, av_s=dom_s%data, av_d=dom_d%data, fldlist='aream') 
    endif

    if (rof_c2_ocn) then
       if (.not.samegrid_ro) then
          gsmap_s => component_get_gsmap_cx(rof(1)) ! gsmap_rx
          dom_s   => component_get_dom_cx(rof(1))   ! dom_rx

          call seq_map_readdata('seq_maps.rc', 'rof2ocn_rmapname:',mpicom_CPLID, CPLID, &
               gsmap_s=gsmap_s, av_s=dom_s%data, avfld_s='aream', filefld_s='area_a', &
               string='rof2ocn aream initialization')
       endif
    end if

    if (lnd_present .and. atm_present) then
       if (samegrid_al) then
          dom_s  => component_get_dom_cx(atm(1))   !dom_ax
          dom_d  => component_get_dom_cx(lnd(1))   !dom_lx

          call seq_map_map(mapper_Sa2l, av_s=dom_s%data, av_d=dom_d%data, fldlist='aream')
       else
          gsmap_d => component_get_gsmap_cx(lnd(1)) ! gsmap_lx
          dom_d   => component_get_dom_cx(lnd(1))   ! dom_lx

          call seq_map_readdata('seq_maps.rc','atm2lnd_fmapname:',mpicom_CPLID, CPLID, &
               gsmap_d=gsmap_d, av_d=dom_d%data, avfld_d='aream', filefld_d='area_b', &
               string='atm2lnd aream initialization')
       endif
    end if

    if (lnd_present .and. glc_present) then
       if (samegrid_lg) then
          dom_s  => component_get_dom_cx(lnd(1))   !dom_lx
          dom_d  => component_get_dom_cx(glc(1))   !dom_gx

          call seq_map_map(mapper_Sl2g, av_s=dom_s%data, av_d=dom_d%data, fldlist='aream')
       else
          gsmap_d => component_get_gsmap_cx(glc(1)) ! gsmap_gx
          dom_d   => component_get_dom_cx(glc(1))   ! dom_gx

          call seq_map_readdata('seq_maps.rc','lnd2glc_fmapname:',mpicom_CPLID, CPLID, &
               gsmap_d=gsmap_d, av_d=dom_d%data, avfld_d='aream', filefld_d='area_b', &
               string='lnd2glc aream initialization')
       endif
    endif

  end subroutine component_init_aream

  !===============================================================================

  subroutine component_init_areacor(comp, samegrid, seq_flds_c2x_fluxes)
    !---------------------------------------------------------------
    ! COMPONENT PES and CPL/COMPONENT (for exchange only)
    !
    ! Uses
    use seq_domain_mct, only : seq_domain_areafactinit
    !
    ! Arguments
    type(component_type) , intent(inout) :: comp(:)
    logical              , intent(in)    :: samegrid
    character(len=*)     , intent(in)    :: seq_flds_c2x_fluxes
    !
    ! Local Variables
    integer :: eci, num_inst
    character(*), parameter :: subname = '(component_init_areacor)'
    !---------------------------------------------------------------

    num_inst = size(comp) 
    do eci = 1,num_inst

       ! For joint cpl-component pes
       if (comp(eci)%iamin_cplcompid) then

          ! Map component domain from coupler to component processes
          call seq_map_map(comp(eci)%mapper_Cx2c, comp(eci)%dom_cx%data, &
               comp(eci)%dom_cc%data, msgtag=comp(eci)%cplcompid*100+eci*10+5)

          ! For only component pes
          if (comp(eci)%iamin_compid) then

             ! Allocate and initialize area correction factors on component processes
             ! Note that the following call allocates comp(eci)%mld2drv(:) and comp(eci)%drv2mdl(:)
             call seq_domain_areafactinit(comp(eci)%dom_cc,           &
                  comp(eci)%mdl2drv, comp(eci)%drv2mdl, samegrid, &
                  comp(eci)%mpicom_compid, comp(eci)%iamroot_compid,  &
                  'areafact_'//comp(eci)%oneletterid//'_'//trim(comp(eci)%name))

             ! Area correct component initialization output fields
             call mct_avect_vecmult(comp(eci)%c2x_cc, comp(eci)%mdl2drv, seq_flds_c2x_fluxes, mask_spval=.true.)

          endif

          ! Map corrected initial component AVs from component to coupler pes
          call seq_map_map(comp(eci)%mapper_cc2x, comp(eci)%c2x_cc, &
               comp(eci)%c2x_cx, msgtag=comp(eci)%cplcompid*100+eci*10+7)

       endif
    enddo

  end subroutine component_init_areacor

  !===============================================================================

#ifndef ESMF_INTERFACE
  subroutine component_run(Eclock, comp, comp_run, infodata,  &
       seq_flds_x2c_fluxes, seq_flds_c2x_fluxes, &
       comp_prognostic, comp_num, timer_barrier, timer_comp_run, &
       run_barriers, ymd, tod, comp_layout)

    !---------------------------------------------------------------
    ! Description
    ! Run component model
    !
    ! Arguments
    type(ESMF_Clock)     , intent(inout)   :: EClock
    type(component_type) , intent(inout)   :: comp(:)
    interface 
       subroutine comp_run( Eclock, cdata, x2c, c2x)
         use ESMF,          only : ESMF_Clock
         use seq_cdata_mod, only : seq_cdata
         use mct_mod,       only : mct_avect
         implicit none
         type(ESMF_Clock), intent(inout) :: EClock
         type(seq_cdata) , intent(inout) :: cdata
         type(mct_aVect) , intent(inout) :: x2c
         type(mct_aVect) , intent(inout) :: c2x   
       end subroutine comp_run
    end interface 
    type (seq_infodata_type) , intent(inout)        :: infodata
    character(len=*)         , intent(in)           :: seq_flds_x2c_fluxes
    character(len=*)         , intent(in)           :: seq_flds_c2x_fluxes
    logical                  , intent(in)           :: comp_prognostic
    integer                  , intent(in), optional :: comp_num
    character(len=*)         , intent(in), optional :: timer_barrier   
    character(len=*)         , intent(in), optional :: timer_comp_run
    logical                  , intent(in), optional :: run_barriers
    integer                  , intent(in), optional :: ymd  ! Current date (YYYYMMDD)
    integer                  , intent(in), optional :: tod  ! Current time of day (seconds) 
    character(len=*)         , intent(in), optional :: comp_layout
    !
    ! Local Variables
    integer  :: eci
    integer  :: ierr
    integer  :: num_inst
    real(r8) :: time_brun         ! Start time
    real(r8) :: time_erun         ! Ending time
    real(r8) :: cktime            ! delta time
    real(r8) :: cktime_acc(10)    ! cktime accumulator array 1 = all, 2 = atm, etc
    integer  :: cktime_cnt(10)    ! cktime counter array
    logical  :: seq_multi_inst    ! a special case of running multiinstances on the same pes. 
    integer  :: phase, phasemin, phasemax  ! phase support
    logical  :: firstloop         ! first time around phase loop
    character(*), parameter :: subname = '(component_run:mct)'
    !---------------------------------------------------------------

    num_inst = size(comp)
    seq_multi_inst = .false.
    phasemin = 1
    phasemax = 1

    if(present(comp_layout)) then
       if(comp_layout .eq. "sequential" .and. num_inst > 1) then
          seq_multi_inst=.true.
          phasemin = 0
       endif
    endif

    do phase = phasemin,phasemax
       if (phase == phasemin) then
          firstloop = .true.
       else
          firstloop = .false.
       endif
       if (comp(1)%oneletterid == 'a') call seq_infodata_putData(infodata, atm_phase=phase)
       if (comp(1)%oneletterid == 'l') call seq_infodata_putData(infodata, lnd_phase=phase)
       if (comp(1)%oneletterid == 'i') call seq_infodata_putData(infodata, ice_phase=phase)
       if (comp(1)%oneletterid == 'o') call seq_infodata_putData(infodata, ocn_phase=phase)
       if (comp(1)%oneletterid == 'r') call seq_infodata_putData(infodata, rof_phase=phase)
       if (comp(1)%oneletterid == 'g') call seq_infodata_putData(infodata, glc_phase=phase)
       if (comp(1)%oneletterid == 'w') call seq_infodata_putData(infodata, wav_phase=phase)

       do eci = 1,num_inst
          if (comp(eci)%iamin_compid) then

             if (present(timer_barrier))  then
                if (present(run_barriers)) then
                   if (run_barriers) then
                      call t_drvstartf (trim(timer_barrier))
                      call mpi_barrier(comp(eci)%mpicom_compid, ierr) 
                      call t_drvstopf (trim(timer_barrier))
                      time_brun = mpi_wtime()
                   endif
                end if
             end if

             if (present(timer_comp_run)) then
                call t_drvstartf (trim(timer_comp_run), barrier=comp(eci)%mpicom_compid)
             end if
             if (drv_threading) call seq_comm_setnthreads(comp(1)%nthreads_compid) 

             if (comp_prognostic .and. firstloop) then
                call mct_avect_vecmult(comp(eci)%x2c_cc, comp(eci)%drv2mdl, seq_flds_x2c_fluxes, mask_spval=.true.)
             end if

             call t_set_prefixf(comp(1)%oneletterid//":")
             call comp_run(EClock, comp(eci)%cdata_cc, comp(eci)%x2c_cc, comp(eci)%c2x_cc)
             call t_unset_prefixf()

             if (phase == 1) then
                call mct_avect_vecmult(comp(eci)%c2x_cc, comp(eci)%mdl2drv, seq_flds_c2x_fluxes, mask_spval=.true.)
             endif

             if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)

             if (present(timer_comp_run)) then
                call t_drvstopf (trim(timer_comp_run))
             end if

             if (present(comp_num)) then
                if (present(run_barriers)) then
                   if (run_barriers) then
                      time_erun = mpi_wtime()
                      cktime = time_erun - time_brun
                      cktime_acc(comp_num) = cktime_acc(comp_num) + cktime
                      cktime_cnt(comp_num) = cktime_cnt(comp_num) + 1
                      if (present(ymd) .and. present(tod)) then
                         write(logunit,107) ' rstamp ',trim(comp(eci)%name),          &
                              '_run_time: model date = ',ymd,tod,                     &
                              ' avg dt = ',cktime_acc(comp_num)/cktime_cnt(comp_num), &
                              ' dt = ',cktime, ' phase = ',phase
                      end if
                   endif
                end if
             end if

          endif
       enddo   ! eci

    enddo   ! phase

107 format( 3A, 2i8, A, f12.4, A, f12.4 )

  end subroutine component_run
#endif

  !===============================================================================

#ifdef ESMF_INTERFACE
  subroutine component_run(Eclock, comp, infodata, &
       seq_flds_x2c_fluxes, seq_flds_c2x_fluxes, &
       comp_prognostic, comp_num, timer_barrier, timer_comp_run, &
       run_barriers, ymd, tod, comp_layout)

    !---------------------------------------------------------------
    !
    ! Arguments
    type(ESMF_Clock)         , intent(inout)        :: EClock
    type(component_type)     , intent(inout)        :: comp(:)
    type (seq_infodata_type) , intent(inout)        :: infodata
    character(len=*)         , intent(in)           :: seq_flds_x2c_fluxes
    character(len=*)         , intent(in)           :: seq_flds_c2x_fluxes
    logical                  , intent(in)           :: comp_prognostic
    integer                  , intent(in), optional :: comp_num
    character(len=*)         , intent(in), optional :: timer_barrier   
    character(len=*)         , intent(in), optional :: timer_comp_run
    logical                  , intent(in), optional :: run_barriers
    integer                  , intent(in), optional :: ymd  ! Current date (YYYYMMDD)
    integer                  , intent(in), optional :: tod  ! Current time of day (seconds) 
    character(len=*)         , intent(in), optional :: comp_layout
    !
    ! Local Variables
    type(ESMF_Array)         :: x2d_array
    type(ESMF_Array)         :: d2x_array
    integer                  :: rc, urc
    integer                  :: eci
    integer                  :: ierr
    integer                  :: num_inst
    real(r8)                 :: time_brun         ! Start time
    real(r8)                 :: time_erun         ! Ending time
    real(r8)                 :: cktime            ! delta time
    real(r8)                 :: cktime_acc(10)    ! cktime accumulator array 1 = all, 2 = atm, etc
    integer                  :: cktime_cnt(10)    ! cktime counter array
    character(*), parameter :: subname = '(component_run:esmf)'
    !---------------------------------------------------------------

    num_inst = size(comp)
    do eci = 1,num_inst
       if (comp(eci)%iamin_compid) then

          if (present(timer_barrier))  then
             if (present(run_barriers)) then
                if (run_barriers) then
                   call t_drvstartf (trim(timer_barrier))
                   call mpi_barrier(comp(eci)%mpicom_compid, ierr) 
                   call t_drvstopf (trim(timer_barrier))
                   time_brun = mpi_wtime()
                endif
             end if
          end if

          if (present(timer_comp_run)) then
             call t_drvstartf (trim(timer_comp_run), barrier=comp(eci)%mpicom_compid)
          end if
          if (drv_threading) call seq_comm_setnthreads(comp(1)%nthreads_compid) 

          ! Put infodata information into export state (NOTE - not into import state)
          call esmfshr_infodata_infodata2state(infodata, comp(eci)%c2x_cc_state, &
               id=comp(eci)%compid, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          ! Determine import state into component
          ! Remember that import state array and import attribute vector share memory now
          call ESMF_StateGet(comp(eci)%x2c_cc_state, itemName="x2d", array=x2d_array, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          ! Apply area correction factor from x2c on mct attribute vector
          if (comp_prognostic) then
             call mct_avect_vecmult(comp(eci)%x2c_cc, comp(eci)%drv2mdl, seq_flds_x2c_fluxes, mask_spval=.true.)
          end if

          ! Convert mct attribute vector to esmf array
          call mct2esmf_copy(comp(eci)%x2c_cc, x2d_array, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          !----------------------------------------------
          ! *** Run the component on component pes***
          !----------------------------------------------
          call t_set_prefixf(comp(1)%oneletterid//":")
          call ESMF_GridCompRun(comp(eci)%gridcomp_cc, &
               importState=comp(eci)%x2c_cc_state, exportState=comp(eci)%c2x_cc_state, &
               clock=EClock, userRc=urc, rc=rc)
          if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
          if (rc  /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc , endflag=ESMF_END_ABORT)
          call t_unset_prefixf()
          !----------------------------------------------

          call ESMF_AttributeSet(comp(eci)%c2x_cc_state, name="ID", value=comp(eci)%compid, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          ! Convert export state back to infodata, the new nextsw_cday is updated in infodata
          call esmfshr_infodata_state2infodata(comp(eci)%c2x_cc_state, infodata, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          ! Determine export state and obtain output esmf array 
          call ESMF_StateGet(comp(eci)%c2x_cc_state, itemName="d2x", array=d2x_array, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

          ! Convert output esmf array to mct attribute vector
          call esmf2mct_copy(d2x_array, comp(eci)%c2x_cc, rc=rc)
          if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
       
          ! Apply area correction for c2x on mct attribute vector
          call mct_avect_vecmult(comp(eci)%c2x_cc, comp(eci)%mdl2drv, seq_flds_c2x_fluxes, mask_spval=.true.)

          if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
          
          if (present(timer_comp_run)) then
             call t_drvstopf (trim(timer_comp_run))
          end if
       end if
    end do

  end subroutine component_run
#endif

  !===============================================================================

#ifndef ESMF_INTERFACE
  subroutine component_final(Eclock, comp, comp_final)

    !---------------------------------------------------------------
    ! Description
    ! Run component model
    !
    ! Arguments
    type(ESMF_Clock)     , intent(inout) :: EClock
    type(component_type) , intent(inout) :: comp(:)
    interface
       subroutine comp_final( Eclock, cdata, x2c, c2x)
         use ESMF,          only : ESMF_Clock
         use seq_cdata_mod, only : seq_cdata
         use mct_mod,       only : mct_avect
         implicit none
         type(ESMF_Clock), intent(inout) :: EClock
         type(seq_cdata) , intent(inout) :: cdata
         type(mct_aVect) , intent(inout) :: x2c
         type(mct_aVect) , intent(inout) :: c2x   
       end subroutine comp_final
    end interface 
    !
    ! Local Variables
    integer :: eci
    integer :: num_inst
    character(*), parameter :: subname = '(component_final:mct)'
    !---------------------------------------------------------------

    num_inst = size(comp)
    do eci = 1,num_inst
       if (comp(eci)%iamin_compid) then
          if (drv_threading) call seq_comm_setnthreads(comp(1)%nthreads_compid) 
          call t_set_prefixf(comp(1)%oneletterid//"_f:")
          call comp_final(EClock, comp(eci)%cdata_cc, comp(eci)%x2c_cc, comp(eci)%c2x_cc)
          call t_unset_prefixf()
          if (drv_threading) call seq_comm_setnthreads(nthreads_GLOID)
       end if
    end do
    
  end subroutine component_final
#endif

  !===============================================================================

#ifdef ESMF_INTERFACE
  subroutine component_final(Eclock, comp)

    !---------------------------------------------------------------
    !
    ! Arguments
    type(ESMF_Clock)     , intent(inout) :: EClock
    type(component_type) , intent(inout) :: comp(:)
    !
    ! Local Variables
    integer             :: eci
    integer             :: rc, urc
    integer             :: num_inst
    character(*), parameter :: subname = '(component_final:esmf)'
    !---------------------------------------------------------------
    call t_set_prefixf(comp(1)%oneletterid//"_f:")

    num_inst = size(comp)
    do eci = 1,num_inst

       ! This calls xxx_final_esmf
       call ESMF_GridCompFinalize(comp(eci)%gridcomp_cc, &
            importState=comp(eci)%x2c_cc_state, exportState=comp(eci)%c2x_cc_state, &
            userRc=urc, rc=rc)
       if (urc /= ESMF_SUCCESS) call ESMF_Finalize(rc=urc, endflag=ESMF_END_ABORT)
       if (rc  /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc , endflag=ESMF_END_ABORT)

    end do

    call t_unset_prefixf()
  end subroutine component_final
#endif

  !===============================================================================

  subroutine component_exch(comp, flow, infodata, infodata_string, &
       mpicom_barrier, run_barriers, &
       timer_barrier, timer_comp_exch, timer_map_exch, timer_infodata_exch) 

    !---------------------------------------------------------------
    ! Description
    ! Map x2m_mx to x2m_mm (component input av from 
    ! coupler processes to component model processes)
    !
    ! Arguments
    implicit none
    type(component_type)    , intent(inout)        :: comp(:)
    character(len=3)        , intent(in)           :: flow
    type(seq_infodata_type) , intent(inout)        :: infodata
    character(len=*)        , intent(in)           :: infodata_string
    integer                 , intent(in), optional :: mpicom_barrier      ! mpicom for barrier call
    logical                 , intent(in), optional :: run_barriers
    character(len=*)        , intent(in), optional :: timer_barrier       ! timer
    character(len=*)        , intent(in), optional :: timer_comp_exch
    character(len=*)        , intent(in), optional :: timer_map_exch
    character(len=*)        , intent(in), optional :: timer_infodata_exch
    !
    ! Local Variables
    integer :: eci
    integer :: ierr
    character(*), parameter :: subname = '(component_exch)'
    !---------------------------------------------------------------

    if (present(timer_barrier))  then
       if (run_barriers) then
          call t_drvstartf (trim(timer_barrier))
          call mpi_barrier(comp(1)%mpicom_cplallcompid,ierr)
          call t_drvstopf (trim(timer_barrier))
       endif
    end if

    if (present(timer_comp_exch)) then
       if (present(mpicom_barrier)) then
          call t_drvstartf (trim(timer_comp_exch), cplcom=.true., barrier=mpicom_barrier)
       end if
    end if

    do eci = 1,size(comp)
       if (comp(eci)%iamin_cplcompid) then
          if (present(timer_map_exch)) then
             call t_drvstartf (trim(timer_map_exch), barrier=comp(eci)%mpicom_cplcompid)
          end if

          if (flow == 'x2c') then ! coupler to component
             call seq_map_map(comp(eci)%mapper_Cx2c, comp(eci)%x2c_cx, comp(eci)%x2c_cc, &
                  msgtag=comp(eci)%cplcompid*100+eci*10+2)
          else if (flow == 'c2x') then ! component to coupler
             call seq_map_map(comp(eci)%mapper_Cc2x, comp(eci)%c2x_cc, comp(eci)%c2x_cx, &
                  msgtag=comp(eci)%cplcompid*100+eci*10+4)
          end if

          if (present(timer_map_exch)) then
             call t_drvstopf (trim(timer_map_exch))
          end if
       endif
    enddo

    if (present(timer_infodata_exch)) then
       call t_drvstartf (trim(timer_infodata_exch), barrier=mpicom_barrier)
    end if
    if (flow == 'c2x') then
       if (comp(1)%iamin_cplcompid) then
          call seq_infodata_exchange(infodata, comp(1)%cplcompid, trim(infodata_string))
       end if
    else if (flow == 'x2c') then
       if (comp(1)%iamin_cplallcompid) then
          call seq_infodata_exchange(infodata, comp(1)%cplallcompid, trim(infodata_string))
       end if
    endif
    if (present(timer_infodata_exch)) then
       call t_drvstopf (trim(timer_infodata_exch))
    end if

    if (present(timer_comp_exch)) then
       if (present(mpicom_barrier)) then
          call t_drvstopf (trim(timer_comp_exch), cplcom=.true.) 
       end if
    end if

  end subroutine component_exch

  !===============================================================================

  subroutine component_diag(infodata, comp, flow, comment, info_debug, timer_diag ) 

    !---------------------------------------------------------------
    ! Description
    ! Component diagnostics for send/recv to coupler
    !
    ! Arguments
    type (seq_infodata_type) , intent(inout)        :: infodata
    type(component_type)     , intent(in)           :: comp(:)
    character(len=3)         , intent(in)           :: flow 
    character(len=*)         , intent(in)           :: comment
    integer                  , intent(in)           :: info_debug
    character(len=*)         , intent(in), optional :: timer_diag
    !
    ! Local Variables
    integer :: eci
    character(*), parameter :: subname = '(component_diag)'
    !---------------------------------------------------------------

    if (info_debug > 1) then
       if (present(timer_diag)) then
          call t_drvstartf (trim(timer_diag), barrier=mpicom_CPLID)
       end if

       do eci = 1,size(comp)
          if (flow == 'x2c') then  ! coupler to component
             call seq_diag_avect_mct(infodata, CPLID, comp(eci)%x2c_cx, &
                  comp(eci)%dom_cx, comp(eci)%gsmap_cx, trim(comment)//comp(eci)%suffix)
          end if
          if (flow == 'c2x') then  ! component to coupler
             call seq_diag_avect_mct(infodata, CPLID, comp(eci)%c2x_cx, &
                  comp(eci)%dom_cx, comp(eci)%gsmap_cx, trim(comment)//comp(eci)%suffix)
          end if
       enddo

       if (present(timer_diag)) then
          call t_drvstopf (trim(timer_diag))
       end if
    endif

  end subroutine component_diag

end module component_mod


