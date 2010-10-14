
#ifndef _MPI_H_
#define _MPI_H_


typedef int MPI_Comm;
typedef int MPI_Request;


#define MPI_COMM_WORLD (1)
#define MPI_COMM_NULL (0)      /* handle 0 maps to NULL */


typedef int MPI_Group;

/* MPI_GROUP_EMPTY and MPI_GROUP_NULL must not conflict with MPI_GROUP_ONE */
#define MPI_GROUP_EMPTY (-1)  
#define MPI_GROUP_NULL  (0)


/*
 * Return codes
 *   On error, mpi-serial aborts so the values don't really matter
 *   as long as they are different than MPI_SUCCESS
 *
 */

#define MPI_SUCCESS        (0)
#define MPI_ERR_BUFFER     (-1)
#define MPI_ERR_COUNT      (-1)
#define MPI_ERR_TYPE       (-1)
#define MPI_ERR_TAG        (-1)
#define MPI_ERR_COMM       (-1)
#define MPI_ERR_RANK       (-1)
#define MPI_ERR_REQUEST    (-1)
#define MPI_ERR_ROOT       (-1)
#define MPI_ERR_GROUP      (-1)
#define MPI_ERR_OP         (-1)
#define MPI_ERR_TOPOLOGY   (-1)
#define MPI_ERR_DIMS       (-1)
#define MPI_ERR_ARG        (-1)
#define MPI_ERR_UNKNOWN    (-1)
#define MPI_ERR_TRUNCATE   (-1)
#define MPI_ERR_OTHER      (-1)
#define MPI_ERR_INTERN     (-1)
#define MPI_PENDING        (-1)
#define MPI_ERR_IN_STATUS  (-1)
#define MPI_ERR_LASTCODE   (-1)


/*
 * MPI_UNDEFINED
 *
 * Uses:
 *   value for "color" in e.g. comm_split
 *   value for rank in Group_translate_ranks
 *
 */


#define MPI_UNDEFINED (-1)


/*
 * Data types etc.
 */

typedef unsigned long int MPI_Aint;
#define MPI_BOTTOM (0)
typedef int MPI_Datatype;


/* The type's value is its size in bytes */

#define MPI_BYTE           (sizeof(char))
#define MPI_CHAR           (sizeof(char))
#define MPI_UNSIGNED_CHAR  (sizeof(unsigned char))
#define MPI_SHORT          (sizeof(short))
#define MPI_UNSIGNED_SHORT (sizeof(unsigned short))
#define MPI_INT            (sizeof(int))
#define MPI_UNSIGNED       (sizeof(unsigned))
#define MPI_LONG           (sizeof(long))
#define MPI_UNSIGNED_LONG  (sizeof(unsigned long))
#define MPI_FLOAT          (sizeof(float))
#define MPI_DOUBLE         (sizeof(double))
#define MPI_LONG_DOUBLE    (sizeof(long double))
#define MPI_PACKED         (sizeof(char))

/* types for MINLOC and MAXLOC */

#define MPI_FLOAT_INT        (sizeof(struct{float a; int b;}))
#define MPI_DOUBLE_INT       (sizeof(struct{double a; int b;}))
#define MPI_LONG_INT         (sizeof(struct{long a; int b;}))
#define MPI_2INT             (sizeof(struct{int a; int b;}))
#define MPI_SHORT_INT        (sizeof (struct{short a; int b;}))
#define MPI_LONG_DOUBLE_INT  (sizeof (struct{long double a; int b;}))

/* size-specific types */

#define MPI_INTEGER1       (1)
#define MPI_INTEGER2       (2)
#define MPI_INTEGER4       (4)
#define MPI_INTEGER8       (8)
#define MPI_INTEGER16      (16)

#define MPI_REAL4          (4)
#define MPI_REAL8          (8)
#define MPI_REAL16         (16)



/*
 * Fortran int size
 *
 */

typedef int MPI_Fint;



#define MPI_ANY_TAG (-1)
#define MPI_ANY_SOURCE (-1)

#define MPI_REQUEST_NULL (0)

#define MPI_MAX_ERROR_STRING (128)
#define MPI_MAX_PROCESSOR_NAME (128)


/*
 * MPI_Status
 *
 * definition must be compatible with the mpif.h values for
 * MPI_STATUS_SIZE, MPI_SOURCE, MPI_TAG, and MPI_ERROR.
 *
 * Note: The type used for MPI_Status_int must be chosen to match
 * Fortran INTEGER.
 *
 */

typedef int MPI_Status_int;

typedef struct                  /* Fortran: INTEGER status(MPI_STATUS_SIZE) */
{
  MPI_Status_int MPI_SOURCE;    /* Fortran: status(MPI_SOURCE) */
  MPI_Status_int MPI_TAG;       /* Fortran: status(MPI_TAG) */
  MPI_Status_int MPI_ERROR;     /* Fortran: status(MPI_ERROR) */

} MPI_Status;


/*
 * Collective operations
 */


typedef int MPI_Op;

typedef void MPI_User_function( void *invec, void *inoutvec, int *len,
				MPI_Datatype *datatype); 

#define MPI_OP_NULL (0)

#define MPI_MAX     (0)
#define MPI_MIN     (0)
#define MPI_SUM     (0)
#define MPI_PROD    (0)
#define MPI_LAND    (0)
#define MPI_BAND    (0)
#define MPI_LOR     (0)
#define MPI_BOR     (0)
#define MPI_LXOR    (0)
#define MPI_BXOR    (0)
#define MPI_MAXLOC  (0)
#define MPI_MINLOC  (0)



/*
 * These are provided for Fortran...
 */


#define MPI_INTEGER           MPI_INT
#define MPI_REAL              MPI_FLOAT
#define MPI_DOUBLE_PRECISION  MPI_DOUBLE

#define MPI_STATUS_SIZE       (sizeof(MPI_Status) / sizeof(int))


/**********************************************************
 *
 * Note: if you need to regenerate the prototypes below,
 * you can use 'protify.awk' and paste the output here.
 *
 */                                                      


extern int MPI_Barrier(MPI_Comm comm );
extern int MPI_Bcast(void* buffer, int count, MPI_Datatype datatype,
                     int root, MPI_Comm comm );
extern int MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                      void* recvbuf, int recvcount, MPI_Datatype recvtype,
                      int root, MPI_Comm comm);
extern int MPI_Gatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, 
                       void* recvbuf, int *recvcounts, int *displs,
                       MPI_Datatype recvtype, int root, MPI_Comm comm);
extern int MPI_Allgather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                         void* recvbuf, int recvcount, MPI_Datatype recvtype,
                         MPI_Comm comm);
extern int MPI_Allgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                          void* recvbuf, int *recvcounts, int *displs,
                          MPI_Datatype recvtype, MPI_Comm comm);
extern int MPI_Scatter( void* sendbuf, int sendcount, MPI_Datatype sendtype, 
                        void* recvbuf, int recvcount, MPI_Datatype recvtype,
                        int root, MPI_Comm comm);
extern int MPI_Scatterv(void* sendbuf, int *sendcounts, int *displs, 
                        MPI_Datatype sendtype, void* recvbuf, int recvcount, 
                        MPI_Datatype recvtype, int root, MPI_Comm comm);
extern int MPI_Reduce(void* sendbuf, void* recvbuf, int count, 
                      MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
extern int MPI_Allreduce(void* sendbuf, void* recvbuf, int count, 
                         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
extern int MPI_Scan( void* sendbuf, void* recvbuf, int count, 
                     MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
extern int MPI_Alltoall(void *sendbuf, int sendcount, MPI_Datatype sendtype,
                        void *recvbuf, int recvcount, MPI_Datatype recvtype,
                        MPI_Comm comm);
extern int MPI_Alltoallv(void *sendbuf, int *sendcounts,
                         int *sdispls, MPI_Datatype sendtype,
                         void *recvbuf, int *recvcounts,
                         int *rdispls, MPI_Datatype recvtype,
                         MPI_Comm comm) ;
extern int MPI_Op_create(MPI_User_function *function, int commute, 
                         MPI_Op *op);
extern MPI_Op MPI_Op_f2c(MPI_Fint op);
extern MPI_Fint MPI_Op_c2f(MPI_Op op);
extern MPI_Comm mpi_comm_new(void);
extern int MPI_Comm_free(MPI_Comm *comm);
extern int MPI_Comm_size(MPI_Comm comm, int *size);
extern int MPI_Comm_rank(MPI_Comm comm, int *rank);
extern int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);
extern int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm);
extern int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm);
extern int MPI_Comm_group(MPI_Comm comm, MPI_Group *group);
extern MPI_Comm MPI_Comm_f2c(MPI_Fint comm);
extern MPI_Fint MPI_Comm_c2f(MPI_Comm comm);
extern int MPI_Group_incl(MPI_Group group, int n, int *ranks, MPI_Group *newgroup);
extern int MPI_Group_range_incl(MPI_Group group, int n, int ranges[][3],
                                MPI_Group *newgroup);
extern int MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup);
extern int MPI_Group_intersection(MPI_Group group1, MPI_Group group2,
                                  MPI_Group *newgroup);
extern int MPI_Group_difference(MPI_Group group1, MPI_Group group2,
                                MPI_Group *newgroup);
extern int MPI_Group_free(MPI_Group *group);
extern int MPI_Group_translate_ranks(MPI_Group group1, int n, int *ranks1,
                                     MPI_Group group2, int *ranks2);
extern MPI_Group MPI_Group_f2c(MPI_Fint group);
extern MPI_Fint MPI_Group_c2f(MPI_Group group);

extern int MPI_Init(int *argc, char **argv[]) ;
extern int MPI_Finalize(void);
extern int MPI_Abort(MPI_Comm comm, int errorcode);
extern int MPI_Error_string(int errorcode, char *string, int *resultlen);
extern int MPI_Get_processor_name(char *name, int *resultlen);
extern int MPI_Initialized(int *flag);
extern int MPI_Pack( void *inbuf, int incount, MPI_Datatype datatype,
                     void *outbuf, int outsize, int *position, MPI_Comm comm);
extern int MPI_Unpack( void *inbuf, int insize, int *position,
                       void *outbuf, int outcount, MPI_Datatype datatype,
                       MPI_Comm comm );
extern int MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
                     int source, int tag, MPI_Comm comm, MPI_Request *request);
extern int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source,
                    int tag, MPI_Comm comm, MPI_Status *status);

extern int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status);
extern int MPI_Wait(MPI_Request *request, MPI_Status *status);
extern int MPI_Waitany(int count, MPI_Request *array_of_requests,
                       int *index, MPI_Status *status);
extern int MPI_Waitall(int count, MPI_Request *array_of_requests,
                       MPI_Status *array_of_statuses);
extern MPI_Request MPI_Request_f2c(MPI_Fint request);
extern MPI_Fint MPI_Request_c2f(MPI_Request request);
extern int MPI_Isend(void *buf, int count, MPI_Datatype datatype,
                     int dest, int tag, MPI_Comm comm, MPI_Request *request) ;
extern int MPI_Send(void* buf, int count, MPI_Datatype datatype,
                    int dest, int tag, MPI_Comm comm);
extern int MPI_Ssend(void* buf, int count, MPI_Datatype datatype,
                     int dest, int tag, MPI_Comm comm);
extern int MPI_Rsend(void* buf, int count, MPI_Datatype datatype,
                     int dest, int tag, MPI_Comm comm);
extern int MPI_Irsend(void *buf, int count, MPI_Datatype datatype,
                     int dest, int tag, MPI_Comm comm, MPI_Request *request) ;
extern double MPI_Wtime(void);;
extern double MPI_Wtime(void);

#endif
