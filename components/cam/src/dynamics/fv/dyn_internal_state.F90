module dyn_internal_state
  use dynamics_vars, only : T_FVDYCORE_STATE, T_FVDYCORE_GRID, &
                            T_FVDYCORE_CONSTANTS, T_FVDYCORE_VARS

  type (T_FVDYCORE_STATE), save, target :: dyn_state 
  private
  public:: get_dyn_state, get_dyn_state_grid, get_dyn_state_vars, get_dyn_state_constants

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function get_dyn_state() result(dynstate)
    type(T_FVDYCORE_STATE), pointer :: dynstate
    dynstate => dyn_state
  end function get_dyn_state

  function get_dyn_state_grid() result(grid)
    type(T_FVDYCORE_GRID), pointer :: grid
    grid => dyn_state%grid
  end function get_dyn_state_grid

  function get_dyn_state_vars() result(vars)
    type(T_FVDYCORE_VARS), pointer :: vars
    vars => dyn_state%vars
  end function get_dyn_state_vars

  function get_dyn_state_constants() result(constants)
    type(T_FVDYCORE_CONSTANTS), pointer :: constants
    constants => dyn_state%constants
  end function get_dyn_state_constants

  subroutine get_dyn_state_grid_bounds(im, jm, km, nq, jfirst, kfirst, jlast, klast, klastp, &
       ifirstxy, ilastxy, jfirstxy, jlastxy)
  integer, intent(out), optional :: im, jm, km, nq, jfirst, kfirst, jlast, klast, klastp, &
       ifirstxy, ilastxy, jfirstxy, jlastxy

  type(T_FVDYCORE_GRID), pointer :: grid

  grid => get_dyn_state_grid()

  if(present(im)) im=grid%im
  if(present(jm)) jm=grid%jm
  if(present(km)) km=grid%km
  if(present(nq)) nq=grid%nq


  if(present(jfirst)) jfirst=grid%jfirst
  if(present(kfirst)) kfirst=grid%kfirst
  if(present(jlast)) jlast=grid%jlast
  if(present(klast)) klast=grid%klast
  if(present(klastp)) klastp=grid%klastp

  if(present(ifirstxy)) ifirstxy=grid%ifirstxy
  if(present(jfirstxy)) jfirstxy=grid%jfirstxy
  if(present(ilastxy)) ilastxy=grid%ilastxy
  if(present(jlastxy)) jlastxy=grid%jlastxy
end subroutine get_dyn_state_grid_bounds



end module dyn_internal_state
