module homme_inputs_mod

  implicit none
  private

  public :: prim_set_test_initial_conditions_f90

contains

  subroutine prim_set_test_initial_conditions_f90 () bind(c)
    use test_mod,          only: set_test_initial_conditions
    use dimensions_mod,    only: nelemd, qsize
    use element_state,     only: FM=>elem_derived_FM, FQ=>elem_derived_FQ, FT=>elem_derived_FT,          &
                                 FPHI=>elem_derived_FPHI, FVTheta=>elem_derived_FVTheta
    use theta_f2c_mod,     only: push_forcing_to_c
    use homme_context_mod, only: elem, hybrid, hvcoord, deriv, tl
    !
    ! Local(s)
    !
    integer :: ie, q

    ! Compute initial condition
    call set_test_initial_conditions(elem,deriv,hybrid,hvcoord,tl,1,nelemd)

    ! Convert state variable Q to (Qdp) because initial conditon reads in Q, not Qdp
    do ie=1,nelemd
      do q=1,qsize
        elem(ie)%state%Qdp(:,:,:,q,1)=elem(ie)%state%Q(:,:,:,q)*elem(ie)%state%dp3d(:,:,:,tl%n0)
        elem(ie)%state%Qdp(:,:,:,q,2)=elem(ie)%state%Q(:,:,:,q)*elem(ie)%state%dp3d(:,:,:,tl%n0)
      enddo
    enddo

    ! Copy all forcing to cxx views
    call push_forcing_to_c(FM, FVTheta, FT, FPHI, FQ)

  end subroutine prim_set_test_initial_conditions_f90

end module homme_inputs_mod
