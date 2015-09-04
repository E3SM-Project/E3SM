module dynEDMod

  !-----------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type
  use landunit_varcon, only : istsoil
  use PatchType      , only : pft
  use ColumnType     , only : col
  use EDVecPatchType , only : EDpft
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
       c = pft%column(p)
       if (col%itype(c) == istsoil) then 
          if ((EDpft%ED_patch(p) == 1 ) .or. (EDpft%ED_bareground(p) == 1)) then
             pft%wtcol(p) = EDpft%wtED(p)
          else
             pft%wtcol(p)  = 0.0_r8 
          end if
       end if
    end do

  end subroutine dyn_ED

end module dynEDMod
