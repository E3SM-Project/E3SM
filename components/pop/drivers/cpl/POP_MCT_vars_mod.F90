module POP_MCT_vars_mod

   use mct_mod
   use kinds_mod

   implicit none
   save
   public

   integer(int_kind) :: POP_MCT_OCNID
   type(mct_gsMap), pointer :: POP_MCT_gsMap_o
   type(mct_gGrid), pointer :: POP_MCT_dom_o

end module POP_MCT_vars_mod
