module gbarrier_mod
  use gbarriertype_mod, only: gbarrier_t
  implicit none

  integer, parameter :: LOG2MAX = 6 
  integer, parameter :: MAXTHREADS = 64

  public :: gbarrier_init
  public :: gbarrier_info
  public :: gbarrier

  contains

    subroutine gbarrier_init(barrier, nthreads)
      type (gbarrier_t), intent(out) :: barrier
      integer, intent(in) :: nthreads

      interface 
        subroutine gbarrier_initialize(c_barrier, nthreads) bind(C)
          use, intrinsic :: ISO_C_Binding, only: C_ptr, C_int
          implicit none

          type (C_ptr), intent(out) :: c_barrier
          integer (C_int), intent(in), value :: nthreads
        end subroutine gbarrier_initialize
      end interface

      call gbarrier_initialize(barrier%c_barrier, nthreads)
    end subroutine gbarrier_init

    subroutine gbarrier_delete(barrier)
      type (gbarrier_t), intent(in) :: barrier

      interface
        subroutine gbarrier_free(c_barrier) bind(C)
          use, intrinsic :: ISO_C_Binding, only: C_ptr
          implicit none

          type (C_ptr), intent(in) :: c_barrier
        end subroutine gbarrier_free
      end interface

      call gbarrier_free(barrier%c_barrier)
    end subroutine gbarrier_delete

    subroutine gbarrier_info(barrier)
      type (gbarrier_t), intent(in) :: barrier

      interface
        subroutine gbarrier_print(c_barrier) bind(C)
          use, intrinsic :: ISO_C_Binding, only: C_ptr
          implicit none
          type (C_ptr), value :: c_barrier
        end subroutine gbarrier_print
      end interface

      call gbarrier_print(barrier%c_barrier)
    end subroutine gbarrier_info


    subroutine gbarrier(barrier, threadID)
      type (gbarrier_t), intent(in) :: barrier
      integer, intent(in) :: threadID

      interface 
        subroutine gbarrier_synchronize(c_barrier, thread) bind(C)
          use, intrinsic :: ISO_C_Binding, only: C_ptr, C_int
          implicit none

          type (C_ptr), intent(in), value :: c_barrier
          integer (C_int), intent(in), value :: thread
        end subroutine gbarrier_synchronize
      end interface

      call gbarrier_synchronize(barrier%c_barrier, threadID)
    end subroutine gbarrier

end module gbarrier_mod

