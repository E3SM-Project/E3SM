#ifdef sysIRIX64
#define Const_MPItype_
#else
#ifdef sysAIX
#define Const_MPItype_
#else
#ifdef sysIRIX
#define Const_MPItype_
#else
#ifdef sysLinux
#define Const_MPItype_
#else
#ifdef sysOSF1
#define Const_MPItype_
#else
#ifdef sysSunOS
#define Const_MPItype_
#else
#ifdef sysSUPER_UX
#define Const_MPItype_
#else
#ifdef sysFujitsu
#define Const_MPItype_
#else
#ifdef sysUNICOS
#define Const_MPItype_
#else
#ifdef Const_MPItype_
#undef Const_MPItype_
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
