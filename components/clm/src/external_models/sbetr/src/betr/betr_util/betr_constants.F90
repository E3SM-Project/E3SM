module betr_constants

  implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  integer, parameter :: stdout = 6
  integer, parameter :: betr_var_name_length = 36
  integer, parameter :: betr_string_length = 128
  integer, parameter :: betr_string_length_long = 256
  integer, parameter :: betr_filename_length = 256
  integer, parameter :: betr_namelist_buffer_size = 4096
  integer, parameter :: betr_namelist_buffer_size_ext = 12288
  integer, parameter :: betr_errmsg_len = 4096
end module betr_constants
