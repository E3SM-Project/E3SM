! nc write 1d array with fillvalues (the holes are explicitly specified)

SUBROUTINE nc_read_1d_implicit_fval_PIO_double_real_kind_fc_double___
USE pio_tutil
   ! pio_decomp_fillval.F90.in:274
  implicit none   ! pio_decomp_fillval.F90.in:275
  integer, parameter :: VEC_LOCAL_SZ = 7   ! pio_decomp_fillval.F90.in:276
  type(var_desc_t)  :: pio_var   ! pio_decomp_fillval.F90.in:277
  type(file_desc_t) :: pio_file   ! pio_decomp_fillval.F90.in:278
  character(len=PIO_TF_MAX_STR_LEN) :: filename   ! pio_decomp_fillval.F90.in:279
  type(io_desc_t) :: wiodesc, riodesc   ! pio_decomp_fillval.F90.in:280
  integer, dimension(:), allocatable :: rcompdof   ! pio_decomp_fillval.F90.in:281
  integer :: rcompdof_sz   ! pio_decomp_fillval.F90.in:282
  integer, dimension(VEC_LOCAL_SZ) :: wcompdof, compdof_rel_disps   ! pio_decomp_fillval.F90.in:283
  real(kind=fc_double), dimension(:), allocatable :: rbuf, exp_val   ! pio_decomp_fillval.F90.in:284
  real(kind=fc_double), dimension(VEC_LOCAL_SZ) :: wbuf   ! pio_decomp_fillval.F90.in:285
  integer, parameter :: BUF_FILLVAL = -2   ! pio_decomp_fillval.F90.in:286
  integer, dimension(1) :: dims   ! pio_decomp_fillval.F90.in:287
  integer :: pio_dim   ! pio_decomp_fillval.F90.in:288
  integer :: i, ierr, lsz   ! pio_decomp_fillval.F90.in:289
  ! iotypes = valid io types   ! pio_decomp_fillval.F90.in:290
  integer, dimension(:), allocatable :: iotypes   ! pio_decomp_fillval.F90.in:291
  character(len=PIO_TF_MAX_STR_LEN), dimension(:), allocatable :: iotype_descs   ! pio_decomp_fillval.F90.in:292
  integer :: num_iotypes   ! pio_decomp_fillval.F90.in:293
   ! pio_decomp_fillval.F90.in:294
  ! compdof is only specified for valid data values, the data holes are   ! pio_decomp_fillval.F90.in:295
  ! implicitly stated (by not specifying them rather than filling it with 0s)   ! pio_decomp_fillval.F90.in:296
  rcompdof_sz = min(pio_tf_world_rank_+1, VEC_LOCAL_SZ)    ! pio_decomp_fillval.F90.in:297
  allocate(rcompdof(rcompdof_sz))   ! pio_decomp_fillval.F90.in:298
  allocate(rbuf(rcompdof_sz))   ! pio_decomp_fillval.F90.in:299
  allocate(exp_val(rcompdof_sz))   ! pio_decomp_fillval.F90.in:300
   ! pio_decomp_fillval.F90.in:301
  do i=1,VEC_LOCAL_SZ   ! pio_decomp_fillval.F90.in:302
    compdof_rel_disps(i) = i   ! pio_decomp_fillval.F90.in:303
  end do   ! pio_decomp_fillval.F90.in:304
  dims(1) = VEC_LOCAL_SZ * pio_tf_world_sz_   ! pio_decomp_fillval.F90.in:305
  ! We only read 1 value on rank0, 2 values on rank1, ...   ! pio_decomp_fillval.F90.in:306
  do i=1,rcompdof_sz   ! pio_decomp_fillval.F90.in:307
    rcompdof(i) = VEC_LOCAL_SZ * pio_tf_world_rank_ + compdof_rel_disps(i)   ! pio_decomp_fillval.F90.in:308
  end do   ! pio_decomp_fillval.F90.in:309
  ! Write everything - we only read some of these values   ! pio_decomp_fillval.F90.in:310
  wcompdof = VEC_LOCAL_SZ * pio_tf_world_rank_ + compdof_rel_disps   ! pio_decomp_fillval.F90.in:311
  wbuf = wcompdof   ! pio_decomp_fillval.F90.in:312
   ! pio_decomp_fillval.F90.in:313
  rbuf = BUF_FILLVAL   ! pio_decomp_fillval.F90.in:314
  do i=1,rcompdof_sz   ! pio_decomp_fillval.F90.in:315
    exp_val(i) = wbuf(i)   ! pio_decomp_fillval.F90.in:316
  end do   ! pio_decomp_fillval.F90.in:317
   ! pio_decomp_fillval.F90.in:318
  call PIO_initdecomp(pio_tf_iosystem_, PIO_double, dims, wcompdof, wiodesc)   ! pio_decomp_fillval.F90.in:319
  call PIO_initdecomp(pio_tf_iosystem_, PIO_double, dims, rcompdof, riodesc)   ! pio_decomp_fillval.F90.in:320
   ! pio_decomp_fillval.F90.in:321
  num_iotypes = 0   ! pio_decomp_fillval.F90.in:322
  call PIO_TF_Get_nc_iotypes(iotypes, iotype_descs, num_iotypes)   ! pio_decomp_fillval.F90.in:323
  filename = "test_pio_decomp_fillval_tests.testfile"   ! pio_decomp_fillval.F90.in:324
  do i=1,1  !num_iotypes   ! pio_decomp_fillval.F90.in:325

    IF (pio_tf_world_rank_ == 0) THEN
      IF (pio_tf_log_level_ >= 0) THEN
        WRITE(*,"(A)",ADVANCE="NO") "PIO_TF: "
        WRITE(*,*)  "Testing : PIO_double : ", iotype_descs(i)
      END IF
    END IF   ! pio_decomp_fillval.F90.in:326
    ierr = PIO_createfile(pio_tf_iosystem_, pio_file, iotypes(i), filename, PIO_CLOBBER)    ! pio_decomp_fillval.F90.in:327

    IF (.NOT. (PIO_TF_Passert_((ierr) == PIO_NOERR))) THEN
      pio_tf_retval_utest_ = -1
      IF (pio_tf_world_rank_ == 0) THEN
        PRINT *, "PIO_TF: PIO Function failed:",&
           "Could not create file " // trim(filename),&
          ":", __FILE__, ":", __LINE__,&
          "(pio_decomp_fillval.F90.in:328)"
      END IF
      RETURN
    END IF   ! pio_decomp_fillval.F90.in:328
   ! pio_decomp_fillval.F90.in:329
    ierr = PIO_def_dim(pio_file, 'PIO_TF_test_dim', dims(1), pio_dim)   ! pio_decomp_fillval.F90.in:330

    IF (.NOT. (PIO_TF_Passert_((ierr) == PIO_NOERR))) THEN
      pio_tf_retval_utest_ = -1
      IF (pio_tf_world_rank_ == 0) THEN
        PRINT *, "PIO_TF: PIO Function failed:",&
           "Failed to define a dim : " // trim(filename),&
          ":", __FILE__, ":", __LINE__,&
          "(pio_decomp_fillval.F90.in:331)"
      END IF
      RETURN
    END IF   ! pio_decomp_fillval.F90.in:331
   ! pio_decomp_fillval.F90.in:332
    ierr = PIO_def_var(pio_file, 'PIO_TF_test_var', PIO_double, (/pio_dim/), pio_var)   ! pio_decomp_fillval.F90.in:333

    IF (.NOT. (PIO_TF_Passert_((ierr) == PIO_NOERR))) THEN
      pio_tf_retval_utest_ = -1
      IF (pio_tf_world_rank_ == 0) THEN
        PRINT *, "PIO_TF: PIO Function failed:",&
           "Failed to define a var : " // trim(filename),&
          ":", __FILE__, ":", __LINE__,&
          "(pio_decomp_fillval.F90.in:334)"
      END IF
      RETURN
    END IF   ! pio_decomp_fillval.F90.in:334
   ! pio_decomp_fillval.F90.in:335
    ierr = PIO_enddef(pio_file)   ! pio_decomp_fillval.F90.in:336

    IF (.NOT. (PIO_TF_Passert_((ierr) == PIO_NOERR))) THEN
      pio_tf_retval_utest_ = -1
      IF (pio_tf_world_rank_ == 0) THEN
        PRINT *, "PIO_TF: PIO Function failed:",&
           "Failed to end redef mode : " // trim(filename),&
          ":", __FILE__, ":", __LINE__,&
          "(pio_decomp_fillval.F90.in:337)"
      END IF
      RETURN
    END IF   ! pio_decomp_fillval.F90.in:337
   ! pio_decomp_fillval.F90.in:338
    call PIO_write_darray(pio_file, pio_var, wiodesc, wbuf, ierr)   ! pio_decomp_fillval.F90.in:339

    IF (.NOT. (PIO_TF_Passert_((ierr) == PIO_NOERR))) THEN
      pio_tf_retval_utest_ = -1
      IF (pio_tf_world_rank_ == 0) THEN
        PRINT *, "PIO_TF: PIO Function failed:",&
           "Failed to write darray : " // trim(filename),&
          ":", __FILE__, ":", __LINE__,&
          "(pio_decomp_fillval.F90.in:340)"
      END IF
      RETURN
    END IF   ! pio_decomp_fillval.F90.in:340
   ! pio_decomp_fillval.F90.in:341
    call PIO_syncfile(pio_file)   ! pio_decomp_fillval.F90.in:342
   ! pio_decomp_fillval.F90.in:343
    ! Read only part of the written data   ! pio_decomp_fillval.F90.in:344
    call PIO_read_darray(pio_file, pio_var, riodesc, rbuf, ierr)   ! pio_decomp_fillval.F90.in:345

    IF (.NOT. (PIO_TF_Passert_((ierr) == PIO_NOERR))) THEN
      pio_tf_retval_utest_ = -1
      IF (pio_tf_world_rank_ == 0) THEN
        PRINT *, "PIO_TF: PIO Function failed:",&
           "Failed to read darray : " // trim(filename),&
          ":", __FILE__, ":", __LINE__,&
          "(pio_decomp_fillval.F90.in:346)"
      END IF
      RETURN
    END IF   ! pio_decomp_fillval.F90.in:346
   ! pio_decomp_fillval.F90.in:347

    IF (.NOT. PIO_TF_Check_val_(rbuf, exp_val)) THEN
      pio_tf_retval_utest_ = -1
      IF (pio_tf_world_rank_ == 0) THEN
        PRINT *, "PIO_TF: PIO Check failed:",&
           "Got wrong val",&
          ":", __FILE__, ":", __LINE__,&
          "(pio_decomp_fillval.F90.in:348)"
      END IF
      RETURN
    END IF   ! pio_decomp_fillval.F90.in:348
   ! pio_decomp_fillval.F90.in:349
    call PIO_closefile(pio_file)   ! pio_decomp_fillval.F90.in:350
    call PIO_deletefile(pio_tf_iosystem_, filename)   ! pio_decomp_fillval.F90.in:351
  end do   ! pio_decomp_fillval.F90.in:352
  if(allocated(iotypes)) then   ! pio_decomp_fillval.F90.in:353
    deallocate(iotypes)   ! pio_decomp_fillval.F90.in:354
    deallocate(iotype_descs)   ! pio_decomp_fillval.F90.in:355
  end if   ! pio_decomp_fillval.F90.in:356
   ! pio_decomp_fillval.F90.in:357
  call PIO_freedecomp(pio_tf_iosystem_, riodesc)   ! pio_decomp_fillval.F90.in:358
  call PIO_freedecomp(pio_tf_iosystem_, wiodesc)   ! pio_decomp_fillval.F90.in:359
  deallocate(exp_val)   ! pio_decomp_fillval.F90.in:360
  deallocate(rbuf)   ! pio_decomp_fillval.F90.in:361
  deallocate(rcompdof)   ! pio_decomp_fillval.F90.in:362
END SUBROUTINE nc_read_1d_implicit_fval_PIO_double_real_kind_fc_double___   ! pio_decomp_fillval.F90.in:363
   ! pio_decomp_fillval.F90.in:363


  SUBROUTINE PIO_TF_Test_driver_
    USE pio_tutil
    IMPLICIT NONE
    pio_tf_retval_utest_ = 0
    IF (pio_tf_world_rank_ == 0) THEN
      PRINT *, "PIO_TF: Starting nc_read_1d_implicit_fval_PIO_double_real_kind_fc_double___"
    END IF
    CALL nc_read_1d_implicit_fval_PIO_double_real_kind_fc_double___()
    IF (pio_tf_retval_utest_ /= 0) THEN
      pio_tf_nerrs_total_ = pio_tf_nerrs_total_ + 1
    END IF
    IF (pio_tf_world_rank_ == 0) THEN
      IF (pio_tf_retval_utest_ == 0) THEN
        WRITE(*,PIO_TF_TEST_RES_FMT) "PIO_TF:Test 30:",&
          "nc_read_1d_implicit_fval_PIO_double_real_kind_fc_double___","-----------", "PASSED"
      ELSE
        WRITE(*,PIO_TF_TEST_RES_FMT) "PIO_TF:Test 30:",&
          "nc_read_1d_implicit_fval_PIO_double_real_kind_fc_double___","-----------", "FAILED"
      END IF
    END IF
  END SUBROUTINE PIO_TF_Test_driver_


  PROGRAM PIO_TF_Test_main_
    USE pio_tutil
    IMPLICIT NONE
    INTEGER ierr

    pio_tf_nerrs_total_=0
    pio_tf_retval_utest_=0
    CALL MPI_Init(ierr)
    CALL PIO_TF_Init_()
    CALL PIO_TF_Test_driver_()
    CALL PIO_TF_Finalize_()
    IF (pio_tf_world_rank_ == 0) THEN
      IF (pio_tf_nerrs_total_ == 0) THEN
        IF (pio_tf_retval_utest_ == 0) THEN
          WRITE(*,PIO_TF_TEST_RES_FMT) "PIO_TF: ",&
           "All tests", "---------", "PASSED"
        ELSE
          WRITE(*,PIO_TF_TEST_RES_FMT) "PIO_TF: ",&
           "Test driver", "---------", "FAILED"
        END IF
      ELSE
        WRITE(*,PIO_TF_TEST_RES_FMT2) "PIO_TF:[",&
          pio_tf_nerrs_total_,"] Tests",&
          "----- FAILED"
      END IF
    END IF
    CALL MPI_Finalize(ierr)
    IF (pio_tf_nerrs_total_ /= 0) THEN
      STOP 99
    END IF
  END PROGRAM
