module cpl_comp_esmf

#ifdef ESMF_INTERFACE

  ! !USES:
  use esmf
  use esmfshr_mod
  use component_type_mod, only: atm, lnd, ice, ocn, wav, rof, glc
  use mct_mod
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private ! except

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: cpl_esmf_init
  public :: cpl_esmf_run
  public :: cpl_esmf_final
  public :: cpl_esmf_register

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------
  private :: cpl_esmf_init_states

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CONTAINS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  subroutine cpl_esmf_register(comp, rc)

    implicit none

    type(ESMF_GridComp), intent(inout)  :: comp
    integer,             intent(out)    :: rc

    rc = ESMF_SUCCESS

    print *, "In seq map register routine"
    ! Register the callback routines.

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_INITIALIZE, &
         cpl_esmf_init, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_RUN, &
         cpl_esmf_run, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_GridCompSetEntryPoint(comp, ESMF_METHOD_FINALIZE, &
         cpl_esmf_final, phase=1, rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine cpl_esmf_register

  !===============================================================================

  subroutine cpl_esmf_init(cplgc, c2x_cx_state, x2c_cx_state, EClock, rc)

    !---------------------------------------------------------------
    use seq_flds_mod
    !
    ! Arguments
    type(ESMF_GridComp)  :: cplgc
    type(ESMF_State)     :: c2x_cx_state
    type(ESMF_State)     :: x2c_cx_state
    type(ESMF_Clock)     :: EClock
    integer, intent(out) :: rc
    !
    ! Local Variables
    integer :: eci
    logical, save :: first_time = .true.
    !---------------------------------------------------------------

    if (first_time) then
       do eci = 1,size(atm)
          if (atm(eci)%present) then
             call cpl_esmf_init_states(atm(eci)%gsmap_cx, &
                  atm(eci)%x2c_cx, atm(eci)%c2x_cx, atm(eci)%dom_cx, &
                  seq_flds_a2x_fields, seq_flds_x2a_fields, &
                  atm(eci)%c2x_cx_array, atm(eci)%x2c_cx_array, atm(eci)%dom_cx_array, &
                  c2x_cx_state, x2c_cx_state)
          end if
       end do

       do eci = 1,size(lnd)
          if (lnd(eci)%present) then
             call cpl_esmf_init_states(lnd(eci)%gsmap_cx, &
                  lnd(eci)%x2c_cx, lnd(eci)%c2x_cx, lnd(eci)%dom_cx, &
                  seq_flds_l2x_fields, seq_flds_x2l_fields, &
                  lnd(eci)%c2x_cx_array, lnd(eci)%x2c_cx_array, lnd(eci)%dom_cx_array, &
                  c2x_cx_state, x2c_cx_state)
          end if
       end do

       do eci = 1,size(ice)
          if (ice(eci)%present) then
             call cpl_esmf_init_states(ice(eci)%gsmap_cx, &
                  ice(eci)%x2c_cx, ice(eci)%c2x_cx, ice(eci)%dom_cx, &
                  seq_flds_i2x_fields, seq_flds_x2i_fields, &
                  ice(eci)%c2x_cx_array, ice(eci)%x2c_cx_array, ice(eci)%dom_cx_array, &
                  c2x_cx_state, x2c_cx_state)
          end if
       end do

       do eci = 1,size(ocn)
          if (ocn(eci)%present) then
             call cpl_esmf_init_states(ocn(eci)%gsmap_cx, &
                  ocn(eci)%x2c_cx, ocn(eci)%c2x_cx, ocn(eci)%dom_cx, &
                  seq_flds_o2x_fields, seq_flds_x2o_fields, &
                  ocn(eci)%c2x_cx_array, ocn(eci)%x2c_cx_array, ocn(eci)%dom_cx_array, &
                  c2x_cx_state, x2c_cx_state)
          end if
       end do

       do eci = 1,size(rof)
          if (rof(eci)%present) then
             call cpl_esmf_init_states(rof(eci)%gsmap_cx, &
                  rof(eci)%x2c_cx, rof(eci)%c2x_cx, rof(eci)%dom_cx, &
                  seq_flds_r2x_fields, seq_flds_x2r_fields, &
                  rof(eci)%c2x_cx_array, rof(eci)%x2c_cx_array, rof(eci)%dom_cx_array, &
                  c2x_cx_state, x2c_cx_state)
          end if
       end do

       do eci = 1,size(glc)
          if (glc(eci)%present) then
             call cpl_esmf_init_states(glc(eci)%gsmap_cx, &
                  glc(eci)%x2c_cx, glc(eci)%c2x_cx, glc(eci)%dom_cx, &
                  seq_flds_g2x_fields, seq_flds_x2g_fields, &
                  glc(eci)%c2x_cx_array, glc(eci)%x2c_cx_array, glc(eci)%dom_cx_array, &
                  c2x_cx_state, x2c_cx_state)
          end if
       end do

       do eci = 1,size(wav)
          if (wav(eci)%present) then
             call cpl_esmf_init_states(wav(eci)%gsmap_cx, &
                  wav(eci)%x2c_cx, wav(eci)%c2x_cx, wav(eci)%dom_cx, &
                  seq_flds_w2x_fields, seq_flds_x2w_fields, &
                  wav(eci)%c2x_cx_array, wav(eci)%x2c_cx_array, wav(eci)%dom_cx_array, &
                  c2x_cx_state, x2c_cx_state)
          end if
       end do
    else
       !TODO - call prep functions
    end if

    rc = ESMF_SUCCESS

  end subroutine cpl_esmf_init

  !===============================================================================

  subroutine cpl_esmf_run(comp, import_state, export_state, EClock, rc)

    !---------------------------------------------------------------
    !
    ! Arguments
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc
    !
    ! Local Variables
    !-------------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! For now this does nothing

  end subroutine cpl_esmf_run

  !===============================================================================

  subroutine cpl_esmf_final(comp, import_state, export_state, EClock, rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:
    type(ESMF_GridComp)          :: comp
    type(ESMF_State)             :: import_state
    type(ESMF_State)             :: export_state
    type(ESMF_Clock)             :: EClock
    integer, intent(out)         :: rc
    !EOP

  end subroutine cpl_esmf_final

  !===============================================================================

  subroutine cpl_esmf_init_states(gsmap_cx, &
       x2c_cx, c2x_cx, dom_cx, &
       seq_flds_c2x_fields, seq_flds_x2c_fields, &
       c2x_cx_array, x2c_cx_array, dom_cx_array, c2x_cx_state, x2c_cx_state)

    !---------------------------------------------------------------
    use seq_flds_mod, only: seq_flds_dom_fields
    !
    ! Arguments
    type(mct_gsmap)  , intent(inout) :: gsmap_cx
    type(mct_avect)  , intent(inout) :: c2x_cx 
    type(mct_avect)  , intent(inout) :: x2c_cx 
    type(mct_ggrid)  , intent(inout) :: dom_cx 
    character(len=*) , intent(in)    :: seq_flds_x2c_fields
    character(len=*) , intent(in)    :: seq_flds_c2x_fields
    type(ESMF_Array) , intent(inout) :: c2x_cx_array 
    type(ESMF_Array) , intent(inout) :: x2c_cx_array 
    type(ESMF_Array) , intent(inout) :: dom_cx_array 
    type(ESMF_State) , intent(inout) :: c2x_cx_state
    type(ESMF_State) , intent(inout) :: x2c_cx_state
    !
    ! Local Variables
    integer             :: lsize
    integer             :: lpet
    integer             :: vm_comm_id, mct_comm_id
    type(ESMF_VM)       :: vm
    integer             :: rc
    integer, pointer    :: gindex(:)
    type(ESMF_DistGrid) :: distgrid_cx
    !---------------------------------------------------------------

    call ESMF_VMGetCurrent(vm, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_VMGet(vm, localPet=lpet, mpiCommunicator=vm_comm_id, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call MPI_Comm_Dup(vm_comm_id, mct_comm_id, rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    lsize = mct_gsMap_lsize(gsmap_cx, comm=mct_comm_id)
    allocate(gindex(lsize), stat=rc)
    if(rc /= 0) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call mct_gsMap_OrderedPoints(gsmap_cx, peno=lpet, points=gindex)

    distgrid_cx = ESMF_DistGridCreate(arbSeqIndexList=gindex, rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    deallocate(gindex)

    c2x_cx_array = ESMF_ArrayCreate(distgrid=distgrid_cx, farrayPtr=c2x_cx%rattr,      &
         distgridToArrayMap=(/2/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    x2c_cx_array = ESMF_ArrayCreate(distgrid=distgrid_cx, farrayPtr=x2c_cx%rattr,      &
         distgridToArrayMap=(/2/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
 
    dom_cx_array = ESMF_ArrayCreate(distgrid=distgrid_cx, farrayPtr=dom_cx%data%rattr, &
         distgridToArrayMap=(/2/), rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(c2x_cx_array, name="mct_names", value=trim(seq_flds_c2x_fields), rc=rc)  
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(x2c_cx_array, name="mct_names", value=trim(seq_flds_x2c_fields), rc=rc)  
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_AttributeSet(dom_cx_array, name="mct_names", value=trim(seq_flds_dom_fields), rc=rc)  
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateAdd(c2x_cx_state, (/c2x_cx_array/), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateAdd(x2c_cx_state, (/x2c_cx_array/), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_StateAdd(c2x_cx_state, (/dom_cx_array/), rc=rc)
    if (rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

  end subroutine cpl_esmf_init_states

#endif

end module cpl_comp_esmf
