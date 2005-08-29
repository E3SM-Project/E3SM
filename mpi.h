
typedef int MPI_Comm;
typedef int MPI_Request;

typedef int MPI_Datatype;
typedef int MPI_Op;


#define MPI_COMM_WORLD (1)
#define MPI_COMM_NULL (0)      /* handle 0 maps to NULL */

#define MPI_SUCCESS   (0)


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

/* types for MINLOC and MAXLOC */

#define MPI_FLOAT_INT        (sizeof(struct{float a; int b;}))
#define MPI_DOUBLE_INT       (sizeof(struct{double a; int b;}))
#define MPI_LONG_INT         (sizeof(struct{long a; int b;}))
#define MPI_2INT             (sizeof(struct{int a; int b;}))
#define MPI_SHORT_INT        (sizeof (struct{short a; int b;}))
#define MPI_LONG_DOUBLE_INT  (sizeof (struct{long double a; int b;}))


#define MPI_ANY_TAG (-1)
#define MPI_ANY_SOURCE (-1)

#define MPI_REQUEST_NULL (0)

#define MPI_MAX_ERROR_STRING (128)


typedef struct        /* Order and size MUST match mpif.h */
{
  int MPI_SOURCE;
  int MPI_TAG;
  int MPI_ERROR;

} MPI_Status;



/* These are provided for Fortran... */

#define MPI_INTEGER           MPI_INT
#define MPI_REAL              MPI_FLOAT
#define MPI_DOUBLE_PRECISION  MPI_DOUBLE

#define MPI_STATUS_SIZE       (sizeof(MPI_Status) / sizeof(int))


/******************************************************  
 * WARNING: The below is automatically generated.     *  
 ******************************************************  
 */                                                      



extern void *mpi_malloc(int size);
extern void mpi_free(void *ptr);

extern int MPI_Comm_free(MPI_Comm *comm);
extern int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm);
extern int MPI_Comm_size(MPI_Comm comm, int *size);
extern int MPI_Comm_rank(MPI_Comm comm, int *rank);
extern int MPI_Init(int *argc, char **argv[]) ;
extern int MPI_Finalize(void);
extern int MPI_Abort(MPI_Comm comm, int errorcode);
extern int MPI_Error_string(int errorcode, char *string, int *resultlen);
extern int MPI_Initialized(int *flag);

extern int MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
                     int source, int tag, MPI_Comm comm, MPI_Request *request);
extern int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source,
                    int tag, MPI_Comm comm, MPI_Status *status);
extern int MPI_Isend(void *buf, int count, MPI_Datatype datatype,
                     int dest, int tag, MPI_Comm comm, MPI_Request *request) ;
extern int MPI_Send(void* buf, int count, MPI_Datatype datatype,
                    int dest, int tag, MPI_Comm comm);

extern int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status);
extern int MPI_Wait(MPI_Request *request, MPI_Status *status);
extern int MPI_Waitany(int count, MPI_Request *array_of_requests,
                       int *index, MPI_Status *status);
extern int MPI_Waitall(int count, MPI_Request *array_of_requests,
                       MPI_Status *array_of_statuses);

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
extern int MPI_Scatterv(void* sendbuf, int *sendcounts, int *displs, 
                        MPI_Datatype sendtype, void* recvbuf, int recvcount, 
                        MPI_Datatype recvtype, int root, MPI_Comm comm);
extern int MPI_Reduce(void* sendbuf, void* recvbuf, int count, 
                      MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);
extern int MPI_Allreduce(void* sendbuf, void* recvbuf, int count, 
                         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);


extern void mpi_destroy_handles(void);
extern void mpi_alloc_handle(int *handle, void **data);
extern void *mpi_handle_to_ptr(int handle);
extern void mpi_free_handle(int handle);

