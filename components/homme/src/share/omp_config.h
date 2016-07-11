#ifdef _OPENMP

#if _OPENMP >= 201307
/* OpenMP 4.0 version */
#define OMP_SIMD !$omp simd

#else

#define OMP_SIMD !dir$ simd
#endif


#else
/* in pure-MPI mode, simd macro is a blank */
#define OMP_SIMD 


#endif

