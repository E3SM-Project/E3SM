#ifdef sysIRIX64
#define Const_MPItype_
#else
#ifdef sysAIX
#define Const_MPItype_
#else
#ifdef sysIRIX
#define Const_MPItype_
#else
#ifdef Const_MPItype_
#undef Const_MPItype_
#endif
#endif
#endif
#endif
