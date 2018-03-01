      Module tests
!/*********************************************************************
! *   Copyright 1996, UCAR/Unidata
! *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
! *   $Id: tests.inc,v 1.15 2007/01/19 16:08:47 ed Exp $
! *********************************************************************/

! Implementation of test.inc in a module
! Modified by: Richard Weed,Ph.D
!              Center for Advanced Vehicular Systems
!              Mississippi State University
!              rweed@cavs.msstate.edu

!!!!
! Do not tabify this unless you like hitting the 72 char limit !!!
!!!
#ifndef UD_TESTS_INC
#define UD_TESTS_INC


!/* The following prevents non-FORTRAN code from appearing in the output. */
#if defined(__osf__)
#   undef _POSIX_SOURCE
#   define _POSIX_SOURCE
#endif

#if defined(NO_NETCDF_2)
#  undef NO_NETCDF_2
#  define NO_NETCDF_2 1
#else
#define NO_NETCDF_2 1
#endif

#include "nfconfig.inc"
#ifdef USE_NETCDF4
       USE netcdf4_f03
#else
       USE netcdf_f03
#endif
 
       Implicit NONE

       SAVE
!/* Parameters of test data */

#ifdef  NF_INT1_T
#   define NF_B 1
#else
#   define NF_B 0
#endif
#ifdef  NF_INT2_T
#   define NF_S 1
#else
#   define NF_S 0
#endif
! Total number of FORTRAN types:
       Integer, Parameter :: NUM_FORTRAN_TYPES = (3 + NF_S + NF_B)
#undef NF_B
#undef NF_S

      Integer, Parameter :: NTYPES=6
      Integer, Parameter :: NDIMS=5
      Integer, Parameter :: NVARS=136
      Integer, Parameter :: NRECS=2
      Integer, Parameter :: NGATTS=NTYPES
      Integer, Parameter :: RECDIM=1
      Integer, Parameter :: MAX_RANK=3
      Integer, Parameter :: MAX_NELS=64
      Integer, Parameter :: MAX_DIM_LEN=4
      Integer, Parameter :: MAX_NATTS=3


!/*
! * Limits of external types (based on those in ncx.h):
! */
      Integer, Parameter :: X_CHAR_MIN=0
      Integer, Parameter :: X_CHAR_MAX=127
      Integer, Parameter :: X_INT1_MIN=(-128)
      Integer, Parameter :: X_INT1_MAX=127
      Integer, Parameter :: X_INT2_MIN=(-32768)
      Integer, Parameter :: X_INT2_MAX=32767
      Integer, Parameter :: X_INT_MIN=(-2147483647-1)
      Integer, Parameter :: X_INT_MAX=2147483647
#if 0
      Real(RK4), Parameter :: X_REAL_MAX=3.4028234663852886e+38
#else
      Real(RK4), Parameter :: X_REAL_MAX=3.4028234663852886e+37
#endif
      Real(RK4), Parameter :: X_REAL_MIN=(-X_REAL_MAX)
#if 0
      Real(RK8), Parameter :: X_DOUBLE_MAX=1.7976931348623157E+308
#else
      Real(RK8), Parameter :: X_DOUBLE_MAX=1.7976931348623157D+200
#endif
      Real(RK8), Parameter :: X_DOUBLE_MIN=(-X_DOUBLE_MAX)

      Integer, Parameter :: X_BYTE_MIN=X_INT1_MIN
      Integer, Parameter :: X_BYTE_MAX=X_INT1_MAX
      Integer, Parameter :: X_SHORT_MIN=X_INT2_MIN
      Integer, Parameter :: X_SHORT_MAX=X_INT2_MAX
      Real(RK4), Parameter :: X_FLOAT_MIN=X_REAL_MIN
      Real(RK4), Parameter :: X_FLOAT_MAX=X_REAL_MAX


!/*
! * Examples of invalid argument values:
! */
      Integer, Parameter :: BAD_ID=-1
      Integer, Parameter :: BAD_DIMID=-1
      Integer, Parameter :: BAD_VARID=-2
      Integer, Parameter :: BAD_ATTNUM=-1
      Integer, Parameter :: BAD_TYPE=0
      Integer, Parameter :: BAD_FILLMODE=-1
      Character*3, Parameter :: BAD_NAME='a/b'


!/*
! * Internal data types:
! */
      Integer, Parameter :: NFT_UNSPECIFIED=0
      Integer, Parameter :: NFT_TEXT=16
      Integer, Parameter :: NFT_CHAR=NFT_TEXT
      Integer, Parameter :: NFT_INT1=17
      Integer, Parameter :: NFT_INT2=18
      Integer, Parameter :: NFT_INT=20
      Integer, Parameter :: NFT_REAL=36
      Integer, Parameter :: NFT_DOUBLE=40


!/*
! * Define a macro for trimming trailing blanks from character variables.
! */
#define TRIM(string) string(1:len_trim(string))


!
! FORTRAN GETARG() subroutine:
!
#ifdef __hpux
#   define      getarg  getarg_
#endif


#endif /* UD_TESTS_INC */


! #include "../fortran/netcdf.inc"


!    /* Global variables - filenames */

      CHARACTER(LEN=80) ::   testfile   !/* netCDF read-only test data */
      CHARACTER(LEN=80) ::   scratch    !/* netCDF test file for writing */

!    /* Global variables - command-line arguments */

      LOGICAL :: CREATE_FILE
      LOGICAL :: READONLY
      LOGICAL :: VERBOSE
      INTEGER :: NFAILS
      INTEGER :: MAX_NMPT        !/* max num messages per test */

!    /* Global variables - test data */

      CHARACTER(LEN=2)            :: DIM_NAME(NDIMS)
      INTEGER                     :: DIM_LEN(NDIMS)
      CHARACTER(LEN=(2+MAX_RANK)) :: VAR_NAME(NVARS)
      INTEGER                     :: VAR_TYPE(NVARS)
      INTEGER                     :: VAR_RANK(NVARS)
      INTEGER                     :: VAR_DIMID(MAX_RANK,NVARS)
      INTEGER                     :: VAR_SHAPE(MAX_RANK,NVARS)
      INTEGER                     :: VAR_NELS(NVARS)
      INTEGER                     :: VAR_NATTS(NVARS)
      CHARACTER(LEN=2)            :: ATTNAME(MAX_NATTS,NVARS)
      CHARACTER(LEN=2)            :: GATT_NAME(NGATTS)
      INTEGER                     :: ATTTYPE(NGATTS,NVARS)
      INTEGER                     :: GATT_TYPE(NGATTS)
      INTEGER                     :: ATTLEN(MAX_NATTS,NVARS)
      INTEGER                     :: GATT_LEN(NGATTS)

!    /* Miscellaneous global variables: */
      CHARACTER(LEN=80) :: PROGNAME        !/* name of the program */
      INTEGER           :: NFAILSTOTAL

!    /* Common blocks for global variables: */

!     COMMON    /LOGCOM/        CREATE_FILE,                            &
!    &                          READONLY,                               & 
!    &                          VERBOSE        

!     COMMON    /TXTCOM/        TESTFILE,                               &
!    &                          SCRATCH,                                &
!    &                          DIM_NAME,                               &
!    &                          VAR_NAME,                               &
!    &                          ATTNAME,                                &
!    &                          GATT_NAME,                              &
!    &                          PROGNAME

!     COMMON    /INTCOM/        NFAILS,                                 &  
!    &                          DIM_LEN,                                &
!    &                          VAR_TYPE,                               &
!    &                          VAR_RANK,                               &
!    &                          VAR_DIMID,                              &
!    &                          VAR_SHAPE,                              &
!    &                          VAR_NELS,                               &
!    &                          VAR_NATTS,                              &
!    &                          ATTTYPE,                                &
!    &                          GATT_TYPE,                              &
!    &                          ATTLEN,                                 &
!    &                          GATT_LEN,                               &
!    &                          MAX_NMPT,                               &
!    &                          NFAILSTOTAL


!    /* Functions for accessing attribute test data */
!    /* varid is -1 for NC_GLOBAL so can do global atts in same loop */

!/*      EXTERNAL       ATT_NAME */

      INTEGER, EXTERNAL           :: VARID
      INTEGER, EXTERNAL           :: NATTS
      CHARACTER(LEN=2), EXTERNAL  :: ATT_NAME
      INTEGER, EXTERNAL           :: ATT_TYPE
      INTEGER, EXTERNAL           :: ATT_LEN


      LOGICAL, EXTERNAL   :: INRANGE
      LOGICAL, EXTERNAL   :: INRANGE_UCHAR
      LOGICAL, EXTERNAL   :: INRANGE_FLOAT
      LOGICAL, EXTERNAL   :: INRANGE3
      LOGICAL, EXTERNAL   :: IN_INTERNAL_RANGE
      LOGICAL, EXTERNAL   :: EQUAL
      LOGICAL, EXTERNAL   :: INT_VEC_EQ

      INTEGER, EXTERNAL   :: ROLL
      INTEGER, EXTERNAL   :: INDEX2INDEXES
      INTEGER, EXTERNAL   :: INDEX2NCINDEXES
      INTEGER, EXTERNAL   :: INDEXES2INDEX
      INTEGER, EXTERNAL   :: NC2DBL
      INTEGER, EXTERNAL   :: DBL2NC
!      INTEGER, EXTERNAL   :: LEN_TRIM
      REAL(RK8), EXTERNAL :: HASH
      REAL(RK8), EXTERNAL :: HASH4
      REAL(RK8), EXTERNAL :: HASH_TEXT
      REAL(RK8), EXTERNAL :: HASH_INT1
      REAL(RK8), EXTERNAL :: HASH_INT2
      REAL(RK8), EXTERNAL :: HASH_INT
      REAL(RK8), EXTERNAL :: HASH_REAL
      REAL(RK8), EXTERNAL :: HASH_DOUBLE
      REAL(RK8), EXTERNAL :: INTERNAL_MIN
      REAL(RK8), EXTERNAL :: INTERNAL_MAX
      REAL(RK8), EXTERNAL :: EXTERNAL_MIN
      REAL(RK8), EXTERNAL :: EXTERNAL_MAX

      End Module tests
