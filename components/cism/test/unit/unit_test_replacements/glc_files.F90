! Trimmed-down version of glc_files including just what is needed for cism unit tests,
! in order to avoid dependencies

module glc_files

  implicit none
  public 
  save

  character(len=*), parameter :: nml_filename = 'cism_in'
end module glc_files
