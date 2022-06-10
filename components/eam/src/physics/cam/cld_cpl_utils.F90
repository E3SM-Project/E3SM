module cld_cpl_utils

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun

  implicit none
  public

contains

   subroutine cld_cpl_register( cld_cpl_opt )

     use physics_buffer, only: pbuf_add_field, dtype_r8
     use ppgrid,         only: pcols, pver

     integer, intent(in) :: cld_cpl_opt

     integer :: idxtmp   ! pbuf component index
 
     if (cld_cpl_opt > 1) then
        call pbuf_add_field(    'T_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field(    'Q_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field(   'QL_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field(   'QI_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field(   'NL_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
        call pbuf_add_field(   'NI_AFT_MACMIC', 'global', dtype_r8, (/pcols,pver/), idxtmp)
     end if

   end subroutine cld_cpl_register

   !---------------------------------------------------------------------------------------------------------
   subroutine set_state_and_ptend( state, pbuf, ztodt, cpair, ptend_dribble )

   use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
   use physics_buffer,   only: physics_buffer_desc, pbuf_get_field
   use physics_buffer,   only: pbuf_get_index
   use constituents,     only: pcnst, cnst_get_ind
   use ppgrid,           only: pcols, pver

   type(physics_state), intent(inout) :: state
   type(physics_buffer_desc), pointer :: pbuf(:)
   real(r8), intent(in)               :: ztodt, cpair
   type(physics_ptend),intent(out)    :: ptend_dribble

   ! local variables

   integer :: ifld
   integer :: ncol
   integer :: ixcldliq, ixcldice, ixnumliq, ixnumice, ixq
   logical :: lq(pcnst)

    real(r8), pointer, dimension(:,:) :: t_after_macmic
    real(r8), pointer, dimension(:,:) :: q_after_macmic
    real(r8), pointer, dimension(:,:) :: ql_after_macmic
    real(r8), pointer, dimension(:,:) :: qi_after_macmic
    real(r8), pointer, dimension(:,:) :: nl_after_macmic
    real(r8), pointer, dimension(:,:) :: ni_after_macmic

    !-----------------------------------------------
    ncol = state%ncol

    ! Look up tracer indices

    call cnst_get_ind('Q',      ixq)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('NUMLIQ', ixnumliq)
    call cnst_get_ind('NUMICE', ixnumice)

    ! Look up indices of pbuf variables

    ifld = pbuf_get_index( 'T_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld,  t_after_macmic )
    ifld = pbuf_get_index( 'Q_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld,  q_after_macmic )
    ifld = pbuf_get_index('QL_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, ql_after_macmic )
    ifld = pbuf_get_index('QI_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, qi_after_macmic )
    ifld = pbuf_get_index('NL_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, nl_after_macmic )
    ifld = pbuf_get_index('NI_AFT_MACMIC'); call pbuf_get_field(pbuf, ifld, ni_after_macmic )

    ! Calculate tendencies to be dribbled, and save in ptend_dribble

    lq(:)        = .FALSE.
    lq(ixq)      = .TRUE.
    lq(ixcldliq) = .TRUE.
    lq(ixcldice) = .TRUE.
    lq(ixnumliq) = .TRUE.
    lq(ixnumice) = .TRUE.

    call physics_ptend_init(ptend_dribble, state%psetcols, 'macmic_dribble_tend', ls=.true., lq=lq)

    ptend_dribble%s(:ncol,:pver)          = (state%t(:ncol,:pver)          -   t_after_macmic(:ncol,:pver))  / ztodt *cpair
    ptend_dribble%q(:ncol,:pver,ixq)      = (state%q(:ncol,:pver,ixq)      -   q_after_macmic(:ncol,:pver))  / ztodt
    ptend_dribble%q(:ncol,:pver,ixcldliq) = (state%q(:ncol,:pver,ixcldliq) -  ql_after_macmic(:ncol,:pver))  / ztodt
    ptend_dribble%q(:ncol,:pver,ixcldice) = (state%q(:ncol,:pver,ixcldice) -  qi_after_macmic(:ncol,:pver))  / ztodt
    ptend_dribble%q(:ncol,:pver,ixnumliq) = (state%q(:ncol,:pver,ixnumliq) -  nl_after_macmic(:ncol,:pver))  / ztodt
    ptend_dribble%q(:ncol,:pver,ixnumice) = (state%q(:ncol,:pver,ixnumice) -  ni_after_macmic(:ncol,:pver))  / ztodt

    ! Reset state back to an old snapshot

    state%t(:ncol,:pver)          =  t_after_macmic(:ncol,:pver)
    state%q(:ncol,:pver,ixq)      =  q_after_macmic(:ncol,:pver)
    state%q(:ncol,:pver,ixcldliq) = ql_after_macmic(:ncol,:pver)
    state%q(:ncol,:pver,ixcldice) = qi_after_macmic(:ncol,:pver)
    state%q(:ncol,:pver,ixnumliq) = nl_after_macmic(:ncol,:pver)
    state%q(:ncol,:pver,ixnumice) = ni_after_macmic(:ncol,:pver)

   end subroutine set_state_and_ptend

   !---------------------------------------------------------------------------------------------------
   subroutine save_state_snapshot_to_pbuf( state, pbuf )

     use physics_types,  only: physics_state
     use physics_buffer, only: physics_buffer_desc, pbuf_get_field, pbuf_get_index
     use ppgrid,         only: pcols, pver
     use constituents,   only: cnst_get_ind

     ! Arguments

     type(physics_state), intent(in)    :: state
     type(physics_buffer_desc), pointer :: pbuf(:)

     ! Local variables

     real(r8), pointer, dimension(:,:) :: ptr2d
     character(len=20) :: varname
     integer :: ncol, ix

     !---------------------
     ncol = state%ncol

     varname = 'T_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
     ptr2d(:ncol,:pver) = state%t(:ncol,:pver)

     varname = 'Q_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
     call cnst_get_ind('Q',      ix); ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ix)

     varname = 'QL_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
     call cnst_get_ind('CLDLIQ', ix); ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ix)

     varname = 'QI_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
     call cnst_get_ind('CLDICE', ix); ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ix)

     varname = 'NL_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
     call cnst_get_ind('NUMLIQ', ix); ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ix)

     varname = 'NI_AFT_MACMIC'; call pbuf_get_field(pbuf, pbuf_get_index(trim(varname)), ptr2d)
     call cnst_get_ind('NUMICE', ix); ptr2d(:ncol,:pver) = state%q(:ncol,:pver,ix)

   end subroutine save_state_snapshot_to_pbuf
   !---------------------------------------------------------------------------------------------------------

end module cld_cpl_utils
