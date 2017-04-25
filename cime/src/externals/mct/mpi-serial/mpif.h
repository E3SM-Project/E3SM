
!
! MPI_COMM_WORLD
!

INTEGER MPI_COMM_WORLD
parameter (mpi_comm_world=1)

!
!
!

integer MPI_BOTTOM
parameter (MPI_BOTTOM=0)


!
! source,tag
  !

  integer MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_TAG_UB
  parameter (mpi_any_source=-1, mpi_any_tag= -1, mpi_tag_ub=1681915906)

  integer MPI_PROC_NULL, MPI_ROOT
  parameter (MPI_PROC_NULL=-2, MPI_ROOT=-3)

  integer MPI_COMM_NULL, MPI_REQUEST_NULL
  parameter (MPI_COMM_NULL=0, MPI_REQUEST_NULL=0)

  integer MPI_GROUP_NULL, MPI_GROUP_EMPTY
  parameter (MPI_GROUP_NULL=0, MPI_GROUP_EMPTY= -1)

  integer MPI_MAX_ERROR_STRING
  parameter (MPI_MAX_ERROR_STRING=128)

  integer MPI_MAX_PROCESSOR_NAME
  parameter (MPI_MAX_PROCESSOR_NAME=128)

  !
  ! Return codes
  !

  integer MPI_SUCCESS
  parameter (MPI_SUCCESS=0)

  integer MPI_ERR_BUFFER
  parameter (MPI_ERR_BUFFER= -1)

  integer MPI_ERR_COUNT
  parameter (MPI_ERR_COUNT= -1)

  integer MPI_ERR_TYPE
  parameter (MPI_ERR_TYPE= -1)

  integer MPI_ERR_TAG
  parameter (MPI_ERR_TAG= -1)

  integer MPI_ERR_COMM
  parameter (MPI_ERR_COMM= -1)

  integer MPI_ERR_RANK
  parameter (MPI_ERR_RANK= -1)

  integer MPI_ERR_REQUEST
  parameter (MPI_ERR_REQUEST= -1)

  integer MPI_ERR_ROOT
  parameter (MPI_ERR_ROOT= -1)

  integer MPI_ERR_GROUP
  parameter (MPI_ERR_GROUP= -1)

  integer MPI_ERR_OP
  parameter (MPI_ERR_OP= -1)

  integer MPI_ERR_TOPOLOGY
  parameter (MPI_ERR_TOPOLOGY= -1)

  integer MPI_ERR_DIMS
  parameter (MPI_ERR_DIMS= -1)

  integer MPI_ERR_ARG
  parameter (MPI_ERR_ARG= -1)

  integer MPI_ERR_UNKNOWN
  parameter (MPI_ERR_UNKNOWN= -1)

  integer MPI_ERR_TRUNCATE
  parameter (MPI_ERR_TRUNCATE= -1)

  integer MPI_ERR_OTHER
  parameter (MPI_ERR_OTHER= -1)

  integer MPI_ERR_INTERN
  parameter (MPI_ERR_INTERN= -1)

  integer MPI_PENDING
  parameter (MPI_PENDING= -1)

  integer MPI_ERR_IN_STATUS
  parameter (MPI_ERR_IN_STATUS= -1)

  integer MPI_ERR_LASTCODE
  parameter (MPI_ERR_LASTCODE= -1)

  integer MPI_ERRORS_RETURN
  parameter (MPI_ERRORS_RETURN= -1)

  !
  !


  integer MPI_UNDEFINED
  parameter (MPI_UNDEFINED= -1)


  !
  ! MPI_Status
  !
  ! The values in this section MUST match the struct definition
  ! in mpi.h
  !


  INTEGER MPI_STATUS_SIZE
  PARAMETER (MPI_STATUS_SIZE=4)

  INTEGER MPI_SOURCE, MPI_TAG, MPI_ERROR
  PARAMETER(MPI_SOURCE=1, MPI_TAG=2, MPI_ERROR=3)
  ! There is a 4th value only used internally

  INTEGER MPI_STATUS_IGNORE(MPI_STATUS_SIZE)
  INTEGER MPI_STATUSES_IGNORE(MPI_STATUS_SIZE,1)
  COMMON /MPISERIAL/ MPI_STATUS_IGNORE
  COMMON /MPISERIAL/ MPI_STATUSES_IGNORE

  !
  ! MPI_IN_PLACE
  !

  INTEGER MPI_IN_PLACE
  COMMON /MPISERIAL/ MPI_IN_PLACE

  SAVE   /MPISERIAL/   ! Technically needed in case goes out of scope


  !
  ! MPI_Datatype values
  !
  !  New datatype values
  !  Type constants represent integer handles, matching up to the index of the
  !  type array equal to the absolute value of the constant plus one.  For
  !  example, MPI_BYTE=-12, corresponding to type index 11.
  !  (Array in type_const.c)
  !


  INTEGER MPI_DATATYPE_NULL
  PARAMETER (MPI_DATATYPE_NULL=0)

  INTEGER MPI_BYTE
  PARAMETER (MPI_BYTE=-12)

  INTEGER MPI_PACKED
  PARAMETER (MPI_PACKED=-13)

  INTEGER MPI_LB
  PARAMETER (MPI_LB=-14)

  INTEGER MPI_UB
  PARAMETER (MPI_UB=-15)

  INTEGER MPI_INTEGER
  PARAMETER (MPI_INTEGER=-16)

  INTEGER MPI_REAL
  PARAMETER (MPI_REAL=-17)

  INTEGER MPI_DOUBLE_PRECISION
  PARAMETER (MPI_DOUBLE_PRECISION=-18)

  INTEGER MPI_COMPLEX
  PARAMETER (MPI_COMPLEX=-19)

  INTEGER MPI_DOUBLE_COMPLEX
  PARAMETER (MPI_DOUBLE_COMPLEX=-20)

  INTEGER MPI_LOGICAL
  PARAMETER (MPI_LOGICAL=-21)

  INTEGER MPI_CHARACTER
  PARAMETER (MPI_CHARACTER=-22)

  integer MPI_2REAL
  parameter (MPI_2REAL= -23)

  integer MPI_2DOUBLE_PRECISION
  parameter (MPI_2DOUBLE_PRECISION= -24)

  integer MPI_2INTEGER
  parameter (MPI_2INTEGER= -25)


  !
  ! Size-specific types
  !

  INTEGER MPI_INTEGER1
  PARAMETER (MPI_INTEGER1= -32 )

  INTEGER MPI_INTEGER2
  PARAMETER (MPI_INTEGER2= -33 )

  INTEGER MPI_INTEGER4
  PARAMETER (MPI_INTEGER4= -34 )

  INTEGER MPI_INTEGER8
  PARAMETER (MPI_INTEGER8= -35 )

  INTEGER MPI_INTEGER16
  PARAMETER (MPI_INTEGER16= -36 )


  INTEGER MPI_REAL4
  PARAMETER (MPI_REAL4= -37 )

  INTEGER MPI_REAL8
  PARAMETER (MPI_REAL8= -38 )

  INTEGER MPI_REAL16
  PARAMETER (MPI_REAL16= -39 )


  integer MPI_COMPLEX8
  parameter (MPI_COMPLEX8= -40 )

  integer MPI_COMPLEX16
  parameter (MPI_COMPLEX16= -41 )

  integer MPI_COMPLEX32
  parameter (MPI_COMPLEX32= -42 )

  integer MPI_LONG_LONG_INT
  parameter (MPI_LONG_LONG_INT= -43)

  integer MPI_LONG_LONG
  parameter (MPI_LONG_LONG= MPI_LONG_LONG_INT)

  integer MPI_UNSIGNED_LONG_LONG
  parameter (MPI_UNSIGNED_LONG_LONG= -44)

  integer MPI_OFFSET
  parameter (MPI_OFFSET= -45)

  !
  ! MPI_Op values
  !
  ! (All are handled as no-op so no value is necessary; but provide one
     ! anyway just in case.)
  !

  INTEGER MPI_SUM
  PARAMETER (MPI_SUM=0)
  INTEGER MPI_MAX
  PARAMETER (MPI_MAX=0)
  INTEGER MPI_MIN
  PARAMETER (MPI_MIN=0)
  INTEGER MPI_PROD
  PARAMETER (MPI_PROD=0)
  INTEGER MPI_LAND
  PARAMETER (MPI_LAND=0)
  INTEGER MPI_BAND
  PARAMETER (MPI_BAND=0)
  INTEGER MPI_LOR
  PARAMETER (MPI_LOR=0)
  INTEGER MPI_BOR
  PARAMETER (MPI_BOR=0)
  INTEGER MPI_LXOR
  PARAMETER (MPI_LXOR=0)
  INTEGER MPI_BXOR
  PARAMETER (MPI_BXOR=0)
  INTEGER MPI_MINLOC
  PARAMETER (MPI_MINLOC=0)
  INTEGER MPI_MAXLOC
  PARAMETER (MPI_MAXLOC=0)
  INTEGER MPI_OP_NULL
  PARAMETER (MPI_OP_NULL=0)

  !
  ! MPI_Wtime
  !

  DOUBLE PRECISION MPI_WTIME
  EXTERNAL MPI_WTIME


  !
  ! Kinds
  !

  INTEGER MPI_OFFSET_KIND
  PARAMETER (MPI_OFFSET_KIND=selected_int_kind(13))

  INTEGER MPI_MODE_RDONLY
  PARAMETER (MPI_MODE_RDONLY=0)

  INTEGER MPI_MODE_CREATE
  PARAMETER (MPI_MODE_CREATE=1)

  INTEGER MPI_MODE_RDWR
  PARAMETER (MPI_MODE_RDWR=2)


  !
  ! Info
  !

  INTEGER MPI_INFO_NULL
  PARAMETER (MPI_INFO_NULL=0)


  !
  ! Library version string (must match C value)
  !

  INTEGER MPI_MAX_LIBRARY_VERSION_STRING
  PARAMETER (MPI_MAX_LIBRARY_VERSION_STRING=80)
