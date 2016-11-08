#ifdef _OPENMP

#if _OPENMP >= 201307

#define OMP_SIMD !$omp simd

#else

#define OMP_SIMD !dir$ simd
#endif


#else

#define OMP_SIMD

#endif

