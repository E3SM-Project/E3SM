module dynEDMod

  !-----------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use landunit_varcon, only : istsoil
  use PatchType      , only : patch
  use ColumnType     , only : col
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  !
  public :: dyn_ED     ! transfers weights calculated internally by ED into wtcol. 
  !------------------------------------------------------------------------
 
contains

  !------------------------------------------------------------------------
  subroutine dyn_ED( bounds )
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! bounds
    
    ! !LOCAL VARIABLES:
    integer  ::  p,c           ! indices
    !------------------------------------------------------------------------
    
    do p = bounds%begp,bounds%endp
       c = patch%column(p)
       if (col%itype(c) == istsoil) then 
          if (patch%is_veg(p) .or. patch%is_bareground(p)) then
             patch%wtcol(p) = patch%wt_ed(p)
          else
             patch%wtcol(p)  = 0.0_r8 
          end if
       end if
    end do

  end subroutine dyn_ED

end module dynEDMod
