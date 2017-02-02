module dynEDMod

  !-----------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use landunit_varcon, only : istsoil
  use PatchType      , only : pft_pp
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
       c = pft_pp%column(p)
       if (col_pp%itype(c) == istsoil) then 
          if ( pft_pp%is_veg(p) .or. pft_pp%is_bareground(p)) then
             pft_pp%wtcol(p) = pft_pp%wt_ed(p)
          else
             pft_pp%wtcol(p) = 0.0_r8 
          end if
       end if
    end do

  end subroutine dyn_ED

end module dynEDMod
