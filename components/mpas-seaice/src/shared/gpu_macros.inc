#ifdef MPAS_OPENACC

#define GPU  acc
#define GPUC acc               /* continuation line macro */
#define GPUF acc )             /* final line of GPU directive */
#define ENTER_DATA  enter data
#define EXIT_DATA   exit data
#define DATA        data
#define DATA_END    end data
#define COPY_IN_LP  copyin(    /* for multi-line variable lists */
#define COPY_OUT_LP copyout(   /* for multi-line variable lists */
#define COPY_DEL_LP delete(
#define UPDATE_D(v) update device v
#define UPDATE_D_LP update device(
#define UPDATE_H(v) update host v   /* !$GPU UPDATE_H((stressDivergenceU, stressDivergenceV)) -> !$acc update host (stressDivergenceU, stressDivergenceV) */
#define UPDATE_H_LP update host(    /* for multi-line variable lists */

#else

#define GPU         omp
#define GPUC        omp&
#define GPUF        omp& )
#define ENTER_DATA  target enter data
#define EXIT_DATA   target exit data
#define DATA        target data
#define DATA_END    end target data
#define COPY_IN_LP  map(to:
#define COPY_OUT_LP map(from:
#define COPY_DEL_LP map(delete:
#define UPDATE_D(v) target update to v
#define UPDATE_D_LP target update to(
#define UPDATE_H(v) target update from v
#define UPDATE_H_LP target update from(

#endif

