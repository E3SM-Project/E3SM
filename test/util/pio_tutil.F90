! PIO Testing framework utilities module
MODULE pio_tutil
  USE pio
  IMPLICIT NONE
  ! Error/Return values
  INTEGER ::  pio_tf_nerrs_total_
  INTEGER ::  pio_tf_retval_utest_
  INTEGER, PARAMETER :: PIO_TF_ERR = 1

  ! IO Processing stuff
  INTEGER, PARAMETER :: MAX_STDIN_ARG_LEN=100
  ENUM, BIND(C)
    ENUMERATOR  ::  IARG_STRIDE_SIDX = 1
    ENUMERATOR  ::  IARG_NUM_IO_TASKS_SIDX
    ENUMERATOR  ::  IARG_NUM_AGGREGATORS_SIDX
    ! Unfortunately since fortran starts with index 1 we need the
    ! hack below. Don't forget to update when adding more argvs
    ENUMERATOR  ::  NUM_IARGS = IARG_NUM_AGGREGATORS_SIDX
    ENUMERATOR  ::  IARG_MAX_SIDX = NUM_IARGS
  ENDENUM

  ! PIO specific info
  INTEGER :: pio_tf_stride_, pio_tf_num_io_tasks_, pio_tf_num_aggregators_
  TYPE(iosystem_desc_t), save :: pio_tf_iosystem_

  ! MPI info
  INTEGER ::  pio_tf_world_rank_, pio_tf_world_sz_
  INTEGER :: pio_tf_comm_

  ! REAL types
  INTEGER, PARAMETER, PUBLIC :: fc_real   = selected_real_kind(6)
  INTEGER, PARAMETER, PUBLIC :: fc_double = selected_real_kind(13)

  ! Misc constants
  INTEGER, PARAMETER :: PIO_TF_MAX_STR_LEN=100

  ! Logging
  INTEGER :: pio_tf_log_level_
  ! PIO_TF_TEST_RES_FMT is used for formatted test result output
  ! -- Useful for writes like
  ! HEADER_STRING, TEST_DESC, FOOTER_STRING, TEST_STATUS
  ! "PIO_TF: Test no: 12", "Test name, desc etc", "-----", "PASSED" 
  CHARACTER(LEN=*), PARAMETER :: PIO_TF_TEST_RES_FMT = "(A20,T22,A40,T64,A6,T72,A6)"
  ! -- Useful for writes like
  ! HEADER_STRING, NUMBER_OF_TESTS, FOOTER_STRING, TEST_STATUS
  ! "PIO_TF: [", 3, "] -----", "FAILED" 
  CHARACTER(LEN=*), PARAMETER :: PIO_TF_TEST_RES_FMT2 = "(A20,T22,I5,T28,A10,T62,A16)"
  CHARACTER(LEN=PIO_TF_MAX_STR_LEN) :: pio_tf_tmp_log_str_

  ! Access modifiers
  ! Public variables
  PUBLIC  :: pio_tf_nerrs_total_, pio_tf_retval_utest_
  PUBLIC  :: pio_tf_iosystem_
  PUBLIC  :: pio_tf_world_rank_, pio_tf_world_sz_
  PUBLIC  :: pio_tf_log_level_
  PUBLIC  :: PIO_TF_TEST_RES_FMT
  PUBLIC  :: pio_tf_tmp_log_str_
  PUBLIC  :: PIO_TF_MAX_STR_LEN
  ! Public functions
  PUBLIC  :: PIO_TF_Init_, PIO_TF_Finalize_, PIO_TF_Passert_
  PUBLIC  :: PIO_TF_Is_netcdf
  PUBLIC  :: PIO_TF_Get_nc_iotypes, PIO_TF_Get_undef_nc_iotypes
  PUBLIC  :: PIO_TF_Get_iotypes, PIO_TF_Get_undef_iotypes
  PUBLIC  :: PIO_TF_Get_data_types
  PUBLIC  :: PIO_TF_Check_val_
  ! Private functions
  PRIVATE :: PIO_TF_Check_arr_int_val, PIO_TF_Check_arr_real_val
  PRIVATE :: PIO_TF_Check_arr_double_val, PIO_TF_Check_char_str_str

  INTERFACE PIO_TF_Check_val_
    MODULE PROCEDURE                  &
        PIO_TF_Check_arr_int_val,     &
        PIO_TF_Check_arr_real_val,    &
        PIO_TF_Check_arr_double_val,  &
        PIO_TF_Check_char_str_str
  END INTERFACE

CONTAINS
  ! Initialize Testing framework - Internal (Not directly used by unit tests)
  SUBROUTINE  PIO_TF_Init_
#ifndef NO_MPIMOD
    use mpi
#else
    include 'mpif.h'
#endif
    INTEGER ierr

    CALL MPI_COMM_DUP(MPI_COMM_WORLD, pio_tf_comm_, ierr);
    CALL MPI_COMM_RANK(pio_tf_comm_, pio_tf_world_rank_, ierr)
    CALL MPI_COMM_SIZE(pio_tf_comm_, pio_tf_world_sz_, ierr)

    pio_tf_log_level_ = 0
    pio_tf_num_aggregators_ = 0
    pio_tf_num_io_tasks_ = 0
    pio_tf_stride_ = 1
    ! Now read input args from rank 0 and bcast it
    ! Args supported are --num-io-tasks, --num-aggregators,
    !   --stride

    CALL Read_input()
    IF (pio_tf_num_io_tasks_ == 0) THEN
      pio_tf_num_io_tasks_ = pio_tf_world_sz_ / pio_tf_stride_;
    END IF
    !IF (pio_tf_world_rank_ == 0) THEN
    !  PRINT *, "PIO_TF: stride=", pio_tf_stride_, ", io_tasks=",&
    !    pio_tf_num_io_tasks_, ", no of aggregators=", &
    !    pio_tf_num_aggregators_
    !END IF

    ! FIXME: Do we need to test with different types of aggregators?
    ! Initialize PIO
    CALL PIO_init(pio_tf_world_rank_, &
          pio_tf_comm_,               &
          pio_tf_num_io_tasks_,       &
          pio_tf_num_aggregators_,    &
          pio_tf_stride_,             &
          PIO_rearr_box,              &
          pio_tf_iosystem_,           &
          base=1)

    ! Set PIO to bcast error
    CALL PIO_seterrorhandling(pio_tf_iosystem_, PIO_BCAST_ERROR)
  END SUBROUTINE PIO_TF_Init_

  ! Finalize Testing framework - Internal (Not directly used by unit tests)
  SUBROUTINE  PIO_TF_Finalize_
#ifndef NO_MPIMOD
    use mpi
#else
    include 'mpif.h'
#endif
    INTEGER ierr
    IF (pio_tf_world_rank_ == 0) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE, pio_tf_nerrs_total_, 1, MPI_INTEGER, MPI_MAX, 0, pio_tf_comm_, ierr)
    ELSE
      CALL MPI_REDUCE(pio_tf_nerrs_total_, pio_tf_nerrs_total_, 1, MPI_INTEGER, MPI_MAX, 0, pio_tf_comm_, ierr)
    END IF

    ! Finalize PIO
    CALL PIO_finalize(pio_tf_iosystem_, ierr);
  END SUBROUTINE PIO_TF_Finalize_

  ! Collective assert - Internal (Not directly used by unit tests)
  ! Each processes passes in its local assert condition and the function
  ! returns the global assert condition
  LOGICAL FUNCTION PIO_TF_Passert_(local_result)
#ifndef NO_MPIMOD
    use mpi
#else
    include 'mpif.h'
#endif
    LOGICAL, INTENT(IN) ::  local_result
    LOGICAL :: global_result
    LOGICAL :: failed, all_failed
    INTEGER :: rank, ierr
    LOGICAL, DIMENSION(:), ALLOCATABLE :: failed_ranks

    CALL MPI_ALLREDUCE(local_result, global_result, 1, MPI_LOGICAL, MPI_LAND, pio_tf_comm_, ierr)
    IF (.NOT. global_result) THEN
      failed = .NOT. local_result
      IF (pio_tf_world_rank_ == 0) THEN
        ALLOCATE(failed_ranks(pio_tf_world_sz_))
      END IF
      ! Gather the ranks where assertion failed
      CALL MPI_GATHER(failed, 1, MPI_LOGICAL, failed_ranks, 1, MPI_LOGICAL, 0, pio_tf_comm_, ierr)

      ! Display the ranks where the assertion failed
      IF (pio_tf_world_rank_ == 0) THEN
        all_failed = .TRUE.
        DO rank=1,pio_tf_world_sz_
          IF (.NOT. failed_ranks(rank)) THEN
            all_failed = .FALSE.
            ! Thank you - f90
            EXIT
          END IF
        END DO
        IF (all_failed) THEN
          PRINT *, "PIO_TF: Fatal Error: Assertion failed on ALL processes"
        ELSE
          PRINT *, "PIO_TF: Fatal Error: Assertion failed on following processes: "
          WRITE(*,"(A8)",ADVANCE="NO") "PIO_TF: "
          DO rank=1,pio_tf_world_sz_
            IF (failed_ranks(rank)) THEN
              WRITE(*,"(I5,A)",ADVANCE="NO") rank-1, ","
            END IF
          END DO
          PRINT *, ""
        END IF
        DEALLOCATE(failed_ranks)
      END IF
    END IF

    PIO_TF_Passert_ = global_result
  END FUNCTION PIO_TF_Passert_

  ! Returns true if iotype is a netcdf type, false otherwise
  LOGICAL FUNCTION PIO_TF_Is_netcdf(iotype)
    INTEGER,  INTENT(IN)  :: iotype
    PIO_TF_Is_netcdf = (iotype == PIO_iotype_netcdf) .or. &
                        (iotype == PIO_iotype_netcdf4p) .or. &
                        (iotype == PIO_iotype_netcdf4c) .or. &
                        (iotype == PIO_iotype_pnetcdf)

  END FUNCTION PIO_TF_Is_netcdf

  ! Returns a list of defined netcdf iotypes
  ! iotypes : After the routine returns contains a list of defined
  !             netcdf types
  ! iotype_descs : After the routine returns contains description of
  !                 the netcdf types returned in iotypes
  ! num_iotypes : After the routine returns contains the number of
  !                 of defined netcdf types, i.e., size of iotypes and
  !                 iotype_descs arrays
  SUBROUTINE PIO_TF_Get_nc_iotypes(iotypes, iotype_descs, num_iotypes)
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: iotypes
    CHARACTER(LEN=*), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: iotype_descs
    INTEGER, INTENT(OUT) :: num_iotypes
    INTEGER :: i

    num_iotypes = 0
    ! First find the number of io types
#ifdef _NETCDF4
      ! netcdf, netcdf4p, netcdf4c
      num_iotypes = num_iotypes + 3
#elif _NETCDF
      ! netcdf
      num_iotypes = num_iotypes + 1
#endif
#ifdef _PNETCDF
      ! pnetcdf
      num_iotypes = num_iotypes + 1
#endif

    ! ALLOCATE with 0 elements ok?
    ALLOCATE(iotypes(num_iotypes))
    ALLOCATE(iotype_descs(num_iotypes))

    i = 1
#ifdef _NETCDF4
      ! netcdf, netcdf4p, netcdf4c
      iotypes(i) = PIO_iotype_netcdf
      iotype_descs(i) = "NETCDF"
      i = i + 1
      iotypes(i) = PIO_iotype_netcdf4p
      iotype_descs(i) = "NETCDF4P"
      i = i + 1
      iotypes(i) = PIO_iotype_netcdf4c
      iotype_descs(i) = "NETCDF4C"
      i = i + 1
#elif _NETCDF
      ! netcdf
      iotypes(i) = PIO_iotype_netcdf
      iotype_descs(i) = "NETCDF"
      i = i + 1
#endif
#ifdef _PNETCDF
      ! pnetcdf
      iotypes(i) = PIO_iotype_pnetcdf
      iotype_descs(i) = "PNETCDF"
      i = i + 1
#endif
  END SUBROUTINE

  ! Returns a list of undefined netcdf iotypes
  ! e.g. This list could be used by a test to make sure that PIO
  !       fails gracefully for undefined types
  ! iotypes : After the routine returns contains a list of undefined
  !             netcdf types
  ! iotype_descs : After the routine returns contains description of
  !                 the netcdf types returned in iotypes
  ! num_iotypes : After the routine returns contains the number of
  !                 of undefined netcdf types, i.e., size of iotypes and
  !                 iotype_descs arrays
  SUBROUTINE PIO_TF_Get_undef_nc_iotypes(iotypes, iotype_descs, num_iotypes)
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: iotypes
    CHARACTER(LEN=*), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: iotype_descs
    INTEGER, INTENT(OUT) :: num_iotypes
    INTEGER :: i

    num_iotypes = 0
    ! First find the number of io types
#ifndef _NETCDF
      ! netcdf
      num_iotypes = num_iotypes + 1
#ifndef _NETCDF4
        ! netcdf4p, netcdf4c
        num_iotypes = num_iotypes + 2
#endif
#endif
#ifndef _PNETCDF
      ! pnetcdf
      num_iotypes = num_iotypes + 1
#endif

    ! ALLOCATE with 0 elements ok?
    ALLOCATE(iotypes(num_iotypes))
    ALLOCATE(iotype_descs(num_iotypes))

    i = 1
#ifndef _NETCDF
      ! netcdf
      iotypes(i) = PIO_iotype_netcdf
      iotype_descs(i) = "NETCDF"
      i = i + 1
#ifndef _NETCDF4
        ! netcdf4p, netcdf4c
        iotypes(i) = PIO_iotype_netcdf4p
        iotype_descs(i) = "NETCDF4P"
        i = i + 1
        iotypes(i) = PIO_iotype_netcdf4c
        iotype_descs(i) = "NETCDF4C"
        i = i + 1
#endif
#endif
#ifndef _PNETCDF
      ! pnetcdf
      iotypes(i) = PIO_iotype_pnetcdf
      iotype_descs(i) = "PNETCDF"
      i = i + 1
#endif
  END SUBROUTINE

  ! Returns a list of defined iotypes
  ! iotypes : After the routine returns contains a list of defined
  !             types
  ! iotype_descs : After the routine returns contains description of
  !                 the types returned in iotypes
  ! num_iotypes : After the routine returns contains the number of
  !                 of defined types, i.e., size of iotypes and
  !                 iotype_descs arrays
  SUBROUTINE PIO_TF_Get_iotypes(iotypes, iotype_descs, num_iotypes)
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: iotypes
    CHARACTER(LEN=*), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: iotype_descs
    INTEGER, INTENT(OUT) :: num_iotypes
    INTEGER :: i

    ! First find the number of io types
    ! binary
    num_iotypes = 1
#ifdef _USEMPIIO
      ! pbinary, direct_pbinary
      num_iotypes = num_iotypes + 2
#endif
#ifdef _NETCDF4
      ! netcdf, netcdf4p, netcdf4c
      num_iotypes = num_iotypes + 3
#elif _NETCDF
      ! netcdf
      num_iotypes = num_iotypes + 1
#endif
#ifdef _PNETCDF
      ! pnetcdf
      num_iotypes = num_iotypes + 1
#endif

    ! ALLOCATE with 0 elements ok?
    ALLOCATE(iotypes(num_iotypes))
    ALLOCATE(iotype_descs(num_iotypes))

    i = 1
    iotypes(i) = PIO_iotype_binary
    iotype_descs(i) = "BINARY"
    i = i + 1

#ifdef _USEMPIIO
      ! pbinary, direct_pbinary
      iotypes(i) = PIO_iotype_pbinary
      iotype_descs(i) = "PBINARY"
      i = i + 1
      iotypes(i) = PIO_iotype_direct_pbinary
      iotype_descs(i) = "DIRECT_PBINARY"
      i = i + 1
#endif
#ifdef _NETCDF4
      ! netcdf, netcdf4p, netcdf4c
      iotypes(i) = PIO_iotype_netcdf
      iotype_descs(i) = "NETCDF"
      i = i + 1
      iotypes(i) = PIO_iotype_netcdf4p
      iotype_descs(i) = "NETCDF4P"
      i = i + 1
      iotypes(i) = PIO_iotype_netcdf4c
      iotype_descs(i) = "NETCDF4C"
      i = i + 1
#elif _NETCDF
      ! netcdf
      iotypes(i) = PIO_iotype_netcdf
      iotype_descs(i) = "NETCDF"
      i = i + 1
#endif
#ifdef _PNETCDF
      ! pnetcdf
      iotypes(i) = PIO_iotype_pnetcdf
      iotype_descs(i) = "PNETCDF"
      i = i + 1
#endif
  END SUBROUTINE

  ! Returns a list of undefined iotypes
  ! e.g. This list could be used by a test to make sure that PIO
  !       fails gracefully for undefined types
  ! iotypes : After the routine returns contains a list of undefined
  !             types
  ! iotype_descs : After the routine returns contains description of
  !                 the types returned in iotypes
  ! num_iotypes : After the routine returns contains the number of
  !                 of undefined types, i.e., size of iotypes and
  !                 iotype_descs arrays
  SUBROUTINE PIO_TF_Get_undef_iotypes(iotypes, iotype_descs, num_iotypes)
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: iotypes
    CHARACTER(LEN=*), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: iotype_descs
    INTEGER, INTENT(OUT) :: num_iotypes
    INTEGER :: i

    ! First find the number of io types
    ! binary is always defined
    num_iotypes = 0
#ifndef _USEMPIIO
      ! pbinary, direct_pbinary
      num_iotypes = num_iotypes + 2
#endif
#ifndef _NETCDF
      ! netcdf
      num_iotypes = num_iotypes + 1
#ifndef _NETCDF4
      ! netcdf4p, netcdf4c
      num_iotypes = num_iotypes + 2
#endif
#endif

#ifndef _PNETCDF
      ! pnetcdf
      num_iotypes = num_iotypes + 1
#endif

    ! ALLOCATE with 0 elements ok?
    ALLOCATE(iotypes(num_iotypes))
    ALLOCATE(iotype_descs(num_iotypes))

    i = 1
#ifndef _USEMPIIO
      ! pbinary, direct_pbinary
      iotypes(i) = PIO_iotype_pbinary
      iotype_descs(i) = "PBINARY"
      i = i + 1
      iotypes(i) = PIO_iotype_direct_pbinary
      iotype_descs(i) = "DIRECT_PBINARY"
      i = i + 1
#endif
#ifndef _NETCDF
      ! netcdf
      iotypes(i) = PIO_iotype_netcdf
      iotype_descs(i) = "NETCDF"
      i = i + 1
#ifndef _NETCDF4
      ! netcdf4p, netcdf4c
      iotypes(i) = PIO_iotype_netcdf4p
      iotype_descs(i) = "NETCDF4P"
      i = i + 1
      iotypes(i) = PIO_iotype_netcdf4c
      iotype_descs(i) = "NETCDF4C"
      i = i + 1
#endif
#endif

#ifndef _PNETCDF
      ! pnetcdf
      iotypes(i) = PIO_iotype_pnetcdf
      iotype_descs(i) = "PNETCDF"
      i = i + 1
#endif
  END SUBROUTINE

  ! Returns a list of PIO base types
  SUBROUTINE PIO_TF_Get_data_types(data_types, data_type_descs, num_data_types)
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: data_types
    CHARACTER(LEN=*), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: data_type_descs
    INTEGER, INTENT(OUT) :: num_data_types

    !num_data_types = 4
    num_data_types = 3
    ALLOCATE(data_types(4))
    ALLOCATE(data_type_descs(4))

    data_types(1) = PIO_double
    data_type_descs(1) = "PIO_double"
    data_types(2) = PIO_real
    data_type_descs(2) = "PIO_real"
    data_types(3) = PIO_int
    data_type_descs(3) = "PIO_int"
    ! FIXME: Check why some PIO functions don't support the char type
    !data_types(4) = PIO_char
    !data_type_descs(4) = "PIO_char"

  END SUBROUTINE PIO_TF_Get_data_types

  LOGICAL FUNCTION PIO_TF_Check_arr_int_val(arr, val)
#ifndef NO_MPIMOD
    USE mpi
#else
    include 'mpif.h'
#endif
    INTEGER, DIMENSION(:), INTENT(IN) :: arr
    INTEGER, INTENT(IN) :: val
    INTEGER :: arr_sz, i, ierr
    ! Not equal at id = nequal_idx
    INTEGER :: nequal_idx
    ! Local and global equal bools
    LOGICAL :: lequal, gequal
    LOGICAL :: failed
    TYPE failed_info
      SEQUENCE
      INTEGER :: idx
      INTEGER :: val
      INTEGER :: exp_val
    END TYPE failed_info
    TYPE (failed_info) :: lfail_info
    TYPE (failed_info), DIMENSION(:), ALLOCATABLE :: gfail_info

    arr_sz = SIZE(arr)
    lequal = .TRUE.;
    gequal = .TRUE.;
    nequal_idx = -1;
    DO i=1, arr_sz
      IF (arr(i) /= val) THEN
        lequal = .FALSE.
        nequal_idx = i
      END IF
    END DO
    CALL MPI_ALLREDUCE(lequal, gequal, 1, MPI_LOGICAL, MPI_LAND, pio_tf_comm_, ierr)
    IF (.NOT. gequal) THEN
      lfail_info % idx = nequal_idx
      IF (nequal_idx /= -1) THEN
        lfail_info % val     = arr(nequal_idx)
        lfail_info % exp_val = val
      END IF
      IF (pio_tf_world_rank_ == 0) THEN
        ALLOCATE(gfail_info(pio_tf_world_sz_))
      END IF
      ! Gather the ranks where assertion failed
      CALL MPI_GATHER(lfail_info, 3, MPI_INTEGER, gfail_info, 3, MPI_INTEGER, 0, pio_tf_comm_, ierr)
      DO i=1,pio_tf_world_sz_
        IF(gfail_info(i) % idx /= -1) THEN
          PRINT *, "PIO_TF: Fatal Error: rank =", i, ", Val[", gfail_info(i) % idx, "]=", gfail_info(i) % val, ", Expected = ", gfail_info(i) % exp_val 
        END IF
      END DO
    END IF

    PIO_TF_Check_arr_int_val = gequal
  END FUNCTION

  LOGICAL FUNCTION PIO_TF_Check_arr_real_val(arr, val)
#ifndef NO_MPIMOD
    USE mpi
#else
    include 'mpif.h'
#endif
    REAL(KIND=fc_real), DIMENSION(:), INTENT(IN) :: arr
    REAL(KIND=fc_real), INTENT(IN) :: val
    INTEGER :: arr_sz, i, ierr
    ! Not equal at id = nequal_idx
    REAL(KIND=fc_real) :: nequal_idx
    ! Local and global equal bools
    LOGICAL :: lequal, gequal
    LOGICAL :: failed
    TYPE failed_info
      SEQUENCE
      REAL(KIND=fc_real) :: idx
      REAL(KIND=fc_real) :: val
      REAL(KIND=fc_real) :: exp_val
    END TYPE failed_info
    TYPE (failed_info) :: lfail_info
    TYPE (failed_info), DIMENSION(:), ALLOCATABLE :: gfail_info

    arr_sz = SIZE(arr)
    lequal = .TRUE.;
    gequal = .TRUE.;
    nequal_idx = -1;
    DO i=1, arr_sz
      IF (arr(i) /= val) THEN
        lequal = .FALSE.
        nequal_idx = i
      END IF
    END DO
    CALL MPI_ALLREDUCE(lequal, gequal, 1, MPI_LOGICAL, MPI_LAND, pio_tf_comm_, ierr)
    IF (.NOT. gequal) THEN
      lfail_info % idx = nequal_idx
      IF (INT(nequal_idx) /= -1) THEN
        lfail_info % val     = arr(INT(nequal_idx))
        lfail_info % exp_val = val
      END IF
      IF (pio_tf_world_rank_ == 0) THEN
        ALLOCATE(gfail_info(pio_tf_world_sz_))
      END IF
      ! Gather the ranks where assertion failed
      CALL MPI_GATHER(lfail_info, 3, MPI_REAL, gfail_info, 3, MPI_REAL, 0, pio_tf_comm_, ierr)
      DO i=1,pio_tf_world_sz_
        IF(INT(gfail_info(i) % idx) /= -1) THEN
          PRINT *, "PIO_TF: Fatal Error: rank =", i, ", Val[", INT(gfail_info(i) % idx), "]=", gfail_info(i) % val, ", Expected = ", gfail_info(i) % exp_val 
        END IF
      END DO
    END IF

    PIO_TF_Check_arr_real_val = gequal
  END FUNCTION

  LOGICAL FUNCTION PIO_TF_Check_arr_double_val(arr, val)
#ifndef NO_MPIMOD
    USE mpi
#else
    include 'mpif.h'
#endif
    REAL(KIND=fc_double), DIMENSION(:), INTENT(IN) :: arr
    REAL(KIND=fc_double), INTENT(IN) :: val
    INTEGER :: arr_sz, i, ierr
    ! Not equal at id = nequal_idx
    REAL(KIND=fc_double) :: nequal_idx
    ! Local and global equal bools
    LOGICAL :: lequal, gequal
    LOGICAL :: failed
    TYPE failed_info
      SEQUENCE
      REAL(KIND=fc_double) :: idx
      REAL(KIND=fc_double) :: val
      REAL(KIND=fc_double) :: exp_val
    END TYPE failed_info
    TYPE (failed_info) :: lfail_info
    TYPE (failed_info), DIMENSION(:), ALLOCATABLE :: gfail_info

    arr_sz = SIZE(arr)
    lequal = .TRUE.;
    gequal = .TRUE.;
    nequal_idx = -1;
    DO i=1, arr_sz
      IF (arr(i) /= val) THEN
        lequal = .FALSE.
        nequal_idx = i
      END IF
    END DO
    CALL MPI_ALLREDUCE(lequal, gequal, 1, MPI_LOGICAL, MPI_LAND, pio_tf_comm_, ierr)
    IF (.NOT. gequal) THEN
      lfail_info % idx = nequal_idx
      IF (INT(nequal_idx) /= -1) THEN
        lfail_info % val     = arr(INT(nequal_idx))
        lfail_info % exp_val = val
      END IF
      IF (pio_tf_world_rank_ == 0) THEN
        ALLOCATE(gfail_info(pio_tf_world_sz_))
      END IF
      ! Gather the ranks where assertion failed
      CALL MPI_GATHER(lfail_info, 3, MPI_DOUBLE, gfail_info, 3, MPI_DOUBLE, 0, pio_tf_comm_, ierr)
      DO i=1,pio_tf_world_sz_
        IF(INT(gfail_info(i) % idx) /= -1) THEN
          PRINT *, "PIO_TF: Fatal Error: rank =", i, ", Val[", INT(gfail_info(i) % idx), "]=", gfail_info(i) % val, ", Expected = ", gfail_info(i) % exp_val 
        END IF
      END DO
    END IF

    PIO_TF_Check_arr_double_val = gequal
  END FUNCTION

  LOGICAL FUNCTION PIO_TF_Check_char_str_str(str1, str2)
    CHARACTER(LEN=*), INTENT(IN)  :: str1
    CHARACTER(LEN=*), INTENT(IN)  :: str2

    PIO_TF_Check_char_str_str = .TRUE.
  END FUNCTION

  ! Parse and process input arguments like "--pio-tf-stride=2" passed
  ! to the unit tests - PRIVATE function
  ! FIXME: Do we need to move input argument processing to a new module?
  SUBROUTINE Parse_and_process_input(argv)
    CHARACTER(LEN=*), INTENT(IN)  :: argv
    INTEGER :: pos

    ! All input arguments are of the form <INPUT_ARG_NAME>=<INPUT_ARG>
    PRINT *, argv
    pos = INDEX(argv, "=")
    IF (pos == 0) THEN
      ! Ignore unrecognized args
      RETURN
    ELSE
      ! Check if it an input to PIO testing framework
      IF (argv(:pos) == "--pio-tf-num-io-tasks=") THEN
        READ(argv(pos+1:), *) pio_tf_num_io_tasks_
      ELSE IF (argv(:pos) == "--pio-tf-num-aggregators=") THEN
        READ(argv(pos+1:), *) pio_tf_num_aggregators_
      ELSE IF (argv(:pos) == "--pio-tf-stride=") THEN
        READ(argv(pos+1:), *) pio_tf_stride_
      ELSE IF (argv(:pos) == "--pio-tf-input-file=") THEN
        PRINT *, "This option is not implemented yet"
      END IF
    END IF

  END SUBROUTINE Parse_and_process_input

  ! Read input arguments - command line, namelist files - common
  ! to all unit test cases and make sure all MPI processes have
  ! access to it - PRIVATE function
  SUBROUTINE Read_input()
#ifndef NO_MPIMOD
    use mpi
#else
    include 'mpif.h'
#endif
    INTEGER :: i, nargs, ierr
    CHARACTER(LEN=MAX_STDIN_ARG_LEN) :: argv
    ! Need to send pio_tf_stride_, pio_tf_num_io_tasks_, pio_tf_num_aggregators_
    INTEGER :: send_buf(NUM_IARGS)

    IF (pio_tf_world_rank_ == 0) THEN
      nargs = COMMAND_ARGUMENT_COUNT()
      DO i=1, nargs
        CALL GET_COMMAND_ARGUMENT(i, argv)
        CALL Parse_and_process_input(argv)
      END DO
      send_buf(IARG_STRIDE_SIDX) = pio_tf_stride_
      send_buf(IARG_NUM_IO_TASKS_SIDX) = pio_tf_num_io_tasks_
      send_buf(IARG_NUM_AGGREGATORS_SIDX) = pio_tf_num_aggregators_
    END IF
    ! Make sure all processes get the input args
    CALL MPI_BCAST(send_buf, NUM_IARGS, MPI_INTEGER, 0, pio_tf_comm_, ierr)
    pio_tf_stride_ = send_buf(IARG_STRIDE_SIDX)
    pio_tf_num_io_tasks_ = send_buf(IARG_NUM_IO_TASKS_SIDX)
    pio_tf_num_aggregators_ = send_buf(IARG_NUM_AGGREGATORS_SIDX)

  END SUBROUTINE Read_input
END MODULE pio_tutil
