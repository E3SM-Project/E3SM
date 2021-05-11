!-----------------------------------------------------------------------
! MPI Context:
!-----------------------------------------------------------------------
!
!  Now only contains cpp defines.  This file should really be called
!     mpi_defines.h
!

#define  MAX_TRF      4
#define  MAX_SMP      4
#define  MAX_PAX      4

#define CPP_MPI_INTEGER MPI_INTEGER
#define CPP_MPI_REAL8 MPI_DOUBLE_PRECISION
#define CPP_MPI_REAL4 MPI_REAL

#if defined(STAND_ALONE)
#  define CPP_REAL4 selected_real_kind(5)
#  define CPP_REAL8 selected_real_kind(12)
#  define CPP_INTEGER4 selected_int_kind(6)
#  define CPP_INTEGER8 selected_int_kind(13)
#else
#  define CPP_REAL4     r4
#  define CPP_REAL8     r8
#  define CPP_INTEGER4  i4
#  define CPP_INTEGER8  i8
#endif
