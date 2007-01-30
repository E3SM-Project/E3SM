
!!!
!!!   NOTE: The files mpif.realXdoubleY.h are generated from
!!!   mpif.master.h using make-mpif and later copied to mpif.h
!!!   during the library make.  All modifications should be
!!!   made to mpif.master.h
!!!


!
! MPI_COMM_WORLD
!

	INTEGER MPI_COMM_WORLD
        parameter (mpi_comm_world=1)


!
! source,tag
!

	integer MPI_ANY_SOURCE, MPI_ANY_TAG
        parameter (mpi_any_source=-1, mpi_any_tag= -1)


        integer MPI_COMM_NULL, MPI_REQUEST_NULL
        parameter (MPI_COMM_NULL=0, MPI_REQUEST_NULL=0)

        integer MPI_GROUP_NULL, MPI_GROUP_EMPTY
        parameter (MPI_GROUP_NULL=0, MPI_GROUP_EMPTY= -1)

        integer MPI_MAX_ERROR_STRING
        parameter (MPI_MAX_ERROR_STRING=128)

        integer MPI_MAX_PROCESSOR_NAME
        parameter (MPI_MAX_PROCESSOR_NAME=128)


        integer MPI_SUCCESS
        parameter (MPI_SUCCESS=0)

        integer MPI_UNDEFINED
        parameter (MPI_UNDEFINED= -1)


!
! MPI_Status
!
! The values in this section MUST match the struct definition
! in mpi.h
!


        INTEGER MPI_STATUS_SIZE
        PARAMETER (MPI_STATUS_SIZE=3)

        INTEGER MPI_SOURCE, MPI_TAG, MPI_ERROR
        PARAMETER(MPI_SOURCE=1, MPI_TAG=2, MPI_ERROR=3)



!
! MPI_Datatype values
!  
! The value is the size of the datatype in bytes.
! Change if necessary for the machine in question.
! (The mpi.h file uses sizeof(), so it should be more
! portable).
! 
!


	INTEGER MPI_BYTE
	PARAMETER (MPI_BYTE=1)

	INTEGER MPI_CHARACTER
	PARAMETER (MPI_CHARACTER=1)

	INTEGER MPI_REAL4
	PARAMETER (MPI_REAL4=4)

	INTEGER MPI_REAL8
	PARAMETER (MPI_REAL8=8)

	INTEGER MPI_INTEGER
	PARAMETER (MPI_INTEGER=4)

	INTEGER MPI_LOGICAL
	PARAMETER (MPI_LOGICAL=4)

!!!!!!!
	INTEGER MPI_REAL
	PARAMETER (MPI_REAL=4)

	INTEGER MPI_DOUBLE_PRECISION
	PARAMETER (MPI_DOUBLE_PRECISION=8)
!!!!!!!

	integer MPI_COMPLEX
	parameter (MPI_COMPLEX=2*MPI_REAL)

        integer MPI_2REAL
        parameter (MPI_2REAL=2*MPI_REAL)

        integer MPI_2DOUBLE_PRECISION
        parameter (MPI_2DOUBLE_PRECISION=2*MPI_DOUBLE_PRECISION)

        integer MPI_2INTEGER
        parameter (MPI_2INTEGER=2*MPI_INTEGER)

!
! MPI_Op values
!
! (All are handled as no-op so no value is necessary)
!

        INTEGER MPI_SUM, MPI_MAX, MPI_MIN, MPI_PROD, MPI_LAND, MPI_BAND
        INTEGER MPI_LOR, MPI_BOR, MPI_LXOR, MPI_BXOR, MPI_MINLOC
        INTEGER MPI_MAXLOC
        INTEGER MPI_OP_NULL

!
! MPI_Wtime
!

        DOUBLE PRECISION MPI_WTIME
        EXTERNAL MPI_WTIME
