#ifdef MPAS_OPENACC

#define GPU 1
#define GPUC 1
#define GPUF 1
#define ENTER_DATA 1
#define EXIT_DATA 1
#define DATA 1
#define DATA_END 1
#define COPY_IN_LP 1
#define COPY_OUT_LP 1
#define COPY_DEL_LP 1
#define UPDATE_D(v) v
#define UPDATE_D_LP 1
#define UPDATE_H(v) v
#define UPDATE_H_LP 1

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

