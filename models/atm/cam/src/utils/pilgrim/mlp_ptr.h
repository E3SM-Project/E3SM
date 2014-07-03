#if defined ( LAHEY )
#if defined ( PTR_INT )
      integer ptrg_t1
      integer ptrg_t2
      integer ptrg_t3
      integer ptrg_t4
      integer ptrg_r8
      integer ptrg_r4
      integer ptrg_i4
#else

!
! Main vars:
!
#define ga_t1_r ga_t1
#define ga_t1_s ga_t1
      pointer (wing_t1, ga_t1)
      real    :: ga_t1(PLON*PLAT*(PLEV+1)*max_nq*max_call)
!
#define ga_t2_r ga_t2
#define ga_t2_s ga_t2
      pointer (wing_t2, ga_t2)
      real    :: ga_t2(PLON*PLAT*(PLEV+1)*max_nq*max_call)
!
#define ga_t3_r ga_t3
#define ga_t3_s ga_t3
      pointer (wing_t3, ga_t3)
      real    :: ga_t3(PLON*PLAT*(PLEV+1)*max_nq*max_call)
!
#define ga_t4_r ga_t4
#define ga_t4_s ga_t4
      pointer (wing_t4, ga_t4)
      real    :: ga_t4(PLON*PLAT*(PLEV+1)*max_nq*max_call)
!
#define ga_r8_r ga_r8
#define ga_r8_s ga_r8
      pointer (wing_r8, ga_r8)
      real    :: ga_r8(PLON*PLAT*(PLEV+1)*max_nq)
!
#define ga_r4_r ga_r4
#define ga_r4_s ga_r4
      pointer (wing_r4, ga_r4)
      real*4  :: ga_r4(PLON*PLAT*(PLEV+1)*max_nq)
!
#define ga_i4_r ga_i4
#define ga_i4_s ga_i4
      pointer (wing_i4, ga_i4)
      integer :: ga_i4(PLON*PLAT*PLEV)
#endif

#if ( ! defined NOT_ASSIGNED )
      wing_t1=ptrg_t1
      wing_t2=ptrg_t2
      wing_t3=ptrg_t3
      wing_t4=ptrg_t4
      wing_r8=ptrg_r8
      wing_r4=ptrg_r4
      wing_i4=ptrg_i4
#endif
#endif
