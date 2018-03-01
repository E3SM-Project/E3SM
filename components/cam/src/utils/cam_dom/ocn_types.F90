module ocn_types

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, begchunk, endchunk
  implicit none
  public

  type frac_t
     real(r8) :: land(pcols)
  end type frac_t

  type ocn_out_t
     real(r8) :: ts(pcols)        ! surface temperature
  end type ocn_out_t

contains

  subroutine ocn_types_alloc(ocn_out)

    type(ocn_out_t), pointer :: ocn_out(:)

    integer :: c

    allocate (ocn_out(begchunk:endchunk))
    do c = begchunk,endchunk
       ocn_out(c)%ts   (:)  = 0.0_r8
    end do

  end subroutine ocn_types_alloc

end module ocn_types
