module cflx

  !-----------------------------------------------------------------------------------------------------
  ! History:
  !  Separated from module clubb_intr, subroutine clubb_surface by Hui Wan (PNNL), 2022
  !-----------------------------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8

  implicit none
  public

contains

    subroutine cflx_tend (state, cam_in, ztodt, ptend, skip_co2, co2_only)

    use physics_types,          only: physics_state, physics_ptend, &
                                      physics_ptend_init, &
                                      set_dry_to_wet, set_wet_to_dry
    use physconst,              only: gravit
    use ppgrid,                 only: pver, pcols
    use constituents,           only: pcnst, cnst_get_ind, cnst_type
    use co2_cycle,              only: co2_cycle_set_cnst_type, co2_transport, c_i
    use camsrfexch,             only: cam_in_t

    implicit none

    ! Input Auguments

    type(physics_state), intent(inout)  :: state                ! Physics state variables
    type(cam_in_t),      intent(in)     :: cam_in               ! contains surface fluxes of constituents
    real(r8),            intent(in)     :: ztodt                ! 2 delta-t        [ s ]
    logical,             intent(in), optional :: skip_co2        ! if .true., skip CO2 tracers (apply all others)
    logical,             intent(in), optional :: co2_only        ! if .true., apply CO2 tracers only

    ! Output Auguments

    type(physics_ptend), intent(out)    :: ptend                ! Individual parameterization tendencies

    ! Local Variables

    integer :: i                                                ! indicees
    integer :: ncol                                             ! # of atmospheric columns

    real(r8) :: tmp1(pcols)
    real(r8) :: rztodt                                          ! 1./ztodt
    integer  :: m, k

    logical  :: lq(pcnst)
    logical  :: l_skip_co2, l_co2_only
    logical  :: co2_mask(pcnst)                                 ! .true. for CO2 tracer indices

    character(len=3), dimension(pcnst) :: cnst_type_loc         ! local override option for constituents cnst_type


    ncol = state%ncol

    ! Process optional arguments
    l_skip_co2 = .false.
    if (present(skip_co2)) l_skip_co2 = skip_co2
    l_co2_only = .false.
    if (present(co2_only)) l_co2_only = co2_only

    ! Build a mask identifying CO2 tracer indices
    co2_mask(:) = .false.
    if (co2_transport()) then
       do k = 1, size(c_i)
          if (c_i(k) >= 1 .and. c_i(k) <= pcnst) co2_mask(c_i(k)) = .true.
       end do
    end if

    !-------------------------------------------------------
    ! Assume 'wet' mixing ratios in surface diffusion code.
    ! don't convert co2 tracers to wet mixing ratios

    cnst_type_loc(:) = cnst_type(:)
    call co2_cycle_set_cnst_type(cnst_type_loc, 'wet')
    call set_dry_to_wet(state, cnst_type_loc)

    !-------------------------------------------------------
    ! Initialize ptend with appropriate tracer mask
    ! skip_co2=.true.: apply all tracers except CO2 (used by tphysbc when cflx_cpl_opt==2)
    ! co2_only=.true.: apply CO2 tracers only (used by tphysac when cflx_cpl_opt==2)
    ! default (both .false.): apply all tracers

    if (l_co2_only) then
       lq(:) = co2_mask(:)
    else
       lq(:) = .TRUE.
       if (l_skip_co2) lq(:) = lq(:) .and. (.not. co2_mask(:))
    end if
    call physics_ptend_init(ptend, state%psetcols, 'clubb_srf', lq=lq)

    !-------------------------------------------------------
    ! Calculate tracer mixing ratio tendencies from cflx

    rztodt                 = 1._r8/ztodt
    ptend%q(:ncol,:pver,:) = state%q(:ncol,:pver,:)
    tmp1(:ncol)            = ztodt * gravit * state%rpdel(:ncol,pver)

    do m = 2, pcnst
      if (.not. lq(m)) cycle
      ptend%q(:ncol,pver,m) = ptend%q(:ncol,pver,m) + tmp1(:ncol) * cam_in%cflx(:ncol,m)
    enddo

    ptend%q(:ncol,:pver,:) = (ptend%q(:ncol,:pver,:) - state%q(:ncol,:pver,:)) * rztodt

    ! Convert tendencies of dry constituents to dry basis.
    do m = 1,pcnst
       if (cnst_type(m).eq.'dry') then
          ptend%q(:ncol,:pver,m) = ptend%q(:ncol,:pver,m)*state%pdel(:ncol,:pver)/state%pdeldry(:ncol,:pver)
       endif
    end do

    !-------------------------------------------------------
    ! convert wet mmr back to dry before conservation check
    ! avoid converting co2 tracers again

    cnst_type_loc(:) = cnst_type(:)
    call co2_cycle_set_cnst_type(cnst_type_loc, 'wet')
    call set_wet_to_dry(state, cnst_type_loc)

    return

  end subroutine cflx_tend

end module
