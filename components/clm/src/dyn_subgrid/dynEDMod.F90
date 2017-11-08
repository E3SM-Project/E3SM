module dynEDMod

  !-----------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use landunit_varcon, only : istsoil
  use VegetationType      , only : veg_pp
  use ColumnType     , only : col_pp
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save

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
       c = veg_pp%column(p)
       if (col_pp%itype(c) == istsoil) then 
          if ( veg_pp%is_veg(p) .or. veg_pp%is_bareground(p)) then
             veg_pp%wtcol(p) = veg_pp%wt_ed(p)
          else
             veg_pp%wtcol(p) = 0.0_r8 
          end if
       end if
    end do

  end subroutine dyn_ED

end module dynEDMod
