!This is where all of the PGI CUDA FORTRAN code will go, and these routines will be called from prim_advection_mod.
!This is compiled regardless, but PGI-specific calls are always wrapped in the _ACCEL ifdefs that are automagically
!activated when -Mcuda is specified during compilation with a PGI compiler. Thus, it will be ignored unless explicitly
!activated by the user
!
!As a general rule, all of the routines in here will be called within a threaded context (assuming ELEMENT_OPENMP is not
!deifned), and therefore, we enforce BARRIERS, MASTERS, and SINGLES from within these routines rather than outside them.
!This is to minimize the visible code impacts on the existing CPU code.

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
#ifdef _ACCEL
!Put everything CUDA-specific in here so it doesn't get compiled without -Mcuda enabled on a PGI compiler
  use cudafor
  use dimensions_mod, only: np,nlevp,nlev,qsize,qsize_d,max_corner_elem,max_neigh_edges,nelemd
  use element_mod, only: timelevels
  implicit none
  private

  !First listed are all externally accescible routines
  public :: cuda_mod_init



contains



  !The point of this is to initialize any data required in other routines of this module as well
  !as to run one initial CUDA kernel just to get those overheads out of the way so that subsequent
  !timing routines are accurage.
  subroutine cuda_mod_init()
    use dimensions_mod, only : np, nlev, qsize
    implicit none
    integer     :: ierr, ie
    type (dim3) :: griddim,blockdim
#ifndef ELEMENT_OPENMP
!$OMP BARRIER
!$OMP MASTER
#endif
    write(*,*) "cuda_mod_init"
    blockdim = dim3(1,1,1)
    griddim  = dim3(1,1,1)
    call warmup <<< griddim , blockdim >>> ( ie )
    ierr = cudaThreadSynchronize()
    write(*,*)"done cuda_mod_init"
#ifndef ELEMENT_OPENMP
!$OMP END MASTER
!$OMP BARRIER
#endif
  end subroutine cuda_mod_init

  !Meaningless kernel just to get initial kernel overheads out of the way.
  attributes(global) subroutine warmup(a)
    integer,value :: a
    a = 2.0 * a
  end subroutine warmup



#endif
end module cuda_mod


