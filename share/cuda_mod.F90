!This is where all of the PGI CUDA FORTRAN code will go, and these routines will be called from prim_advection_mod.
!This is compiled regardless, but PGI-specific calls are always wrapped in the _ACCEL ifdefs that are automagically
!activated when -Mcuda is specified during compilation with a PGI compiler. Thus, it will be ignored unless explicitly
!activated by the user

! Please pay attention to this all caps passive aggresive banner.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!                     !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  STATUS INCOMPLETE  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  DO NOT USE YET     !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  UNTIL THIS BANNER  !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!  IS REMOVED         !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!                     !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module cuda_mod
  implicit none
#ifdef _ACCEL
  !Put everything CUDA-specific in here so it doesn't get compiled without -Mcuda enabled on a PGI compiler
  use cudafor





#endif
end module cuda_mod


