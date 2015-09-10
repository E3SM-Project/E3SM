module column_types_mod

  use hybvcoord_mod,   only : hvcoord_t 
  use hybrid_mod,      only : hybrid_t
  use time_mod,        only : TimeLevel_t
  use dimensions_mod,  only : nlev, nlevp, np
  use kinds,           only : real_kind, int_kind

  implicit none

  public 


  type, public :: HeldSuarezForcing_t
     logical                    :: INIT = .FALSE.
  end type HeldSuarezForcing_t

  type, public :: ColumnModel_t
     type (HeldSuarezForcing_t)        :: cm_hs
     type (hvcoord_t)                  :: hvcoord
     type (hybrid_t), pointer          :: hybrid
     type (TimeLevel_t), pointer       :: tl
     integer                           :: nets,nete
  end type ColumnModel_t
  

end module column_types_mod
