module ocn_types

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols
  implicit none
  public

  type ocn_out_t
     real(r8) :: ts(pcols)        ! surface temperature
  end type ocn_out_t

end module ocn_types
