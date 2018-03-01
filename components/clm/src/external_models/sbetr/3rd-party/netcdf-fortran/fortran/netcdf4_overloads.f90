  ! Overload fill value functions

  interface nf90_def_var_fill
    module procedure nf90_def_var_fill_OneByteInt,   &
                     nf90_def_var_fill_TwoByteInt,   &
                     nf90_def_var_fill_FourByteInt,  &
                     nf90_def_var_fill_EightByteInt, &
                     nf90_def_var_fill_FourByteReal, &
                     nf90_def_var_fill_EightByteReal
   end interface

  interface nf90_inq_var_fill
    module procedure nf90_inq_var_fill_OneByteInt,   &
                     nf90_inq_var_fill_TwoByteInt,   &
                     nf90_inq_var_fill_FourByteInt,  &
                     nf90_inq_var_fill_EightByteInt, &
                     nf90_inq_var_fill_FourByteReal, &
                     nf90_inq_var_fill_EightByteReal
   end interface
  
