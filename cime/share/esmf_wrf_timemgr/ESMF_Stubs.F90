! Various dummy type definitions and routines for the sole purpose of 
! mimicking newer ESMF interface features without necessarily implementing 
! them.  

MODULE ESMF_Stubs

   IMPLICIT NONE

   PRIVATE

! Bogus typedefs
   TYPE ESMF_Grid
      INTEGER :: dummy
   END TYPE

   TYPE ESMF_GridComp
      INTEGER :: dummy
   END TYPE

   TYPE ESMF_State
      INTEGER :: dummy
   END TYPE

   TYPE ESMF_VM
      INTEGER :: dummy
   END TYPE

   TYPE ESMF_END_FLAG
      INTEGER :: dummy
   END TYPE
   TYPE(ESMF_END_FLAG), PARAMETER ::  &
      ESMF_END_ABORT   = ESMF_END_FLAG(1), &
      ESMF_END_NORMAL  = ESMF_END_FLAG(2), &
      ESMF_END_KEEPMPI = ESMF_END_FLAG(3)

   TYPE ESMF_MsgType
      INTEGER :: mtype
   END TYPE
   TYPE(ESMF_MsgType), PARAMETER  ::      &
      ESMF_LOG_INFO  =   ESMF_MsgType(1), &
      ESMF_LOG_WARNING = ESMF_MsgType(2), &
      ESMF_LOG_ERROR =   ESMF_MsgType(3)

   TYPE ESMF_LOG
      INTEGER :: dummy
   END TYPE

   LOGICAL, private, save :: initialized = .false.

   PUBLIC ESMF_Grid, ESMF_GridComp, ESMF_State, ESMF_VM
   PUBLIC ESMF_Initialize, ESMF_Finalize, ESMF_IsInitialized
   PUBLIC ESMF_LogWrite, ESMF_LOG, ESMF_MsgType, ESMF_END_FLAG
   PUBLIC ESMF_LOG_INFO, ESMF_LOG_WARNING, ESMF_LOG_ERROR
   PUBLIC ESMF_END_ABORT, ESMF_END_NORMAL, ESMF_END_KEEPMPI

CONTAINS


! NOOP
   SUBROUTINE ESMF_Initialize( vm, defaultCalendar, rc )
      USE ESMF_BaseMod
      USE ESMF_CalendarMod
!     USE ESMF_TimeMod,     only: defaultCal
      TYPE(ESMF_VM),           INTENT(IN   ), OPTIONAL :: vm
      TYPE(ESMF_CalKind_Flag), INTENT(IN   ), OPTIONAL :: defaultCalendar
      INTEGER,                 INTENT(  OUT), OPTIONAL :: rc

      TYPE(ESMF_CalKind_Flag) :: defaultCalType
      INTEGER :: status

      IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
      ! Initialize the default time manager calendar
      IF ( PRESENT(defaultCalendar) )THEN
         defaultCalType = defaultCalendar
      ELSE
         defaultCalType = ESMF_CALKIND_NOLEAP
      END IF
      allocate( defaultCal )
!      write(6,*) 'tcx1 ESMF_Stubs defcal ',defaultcaltype%caltype
!      call flush(6)
      defaultCal = ESMF_CalendarCreate( calkindflag=defaultCalType, &
                        rc=status)
!      write(6,*) 'tcx2 ESMF_Stubs defcal ',defaultcal%type%caltype
!      call flush(6)
      allocate( gregorianCal )
!      write(6,*) 'tcx1 ESMF_Stubs grcal ',esmf_calkind_gregorian%caltype
!      call flush(6)
      gregorianCal = ESMF_CalendarCreate( calkindflag=ESMF_CALKIND_GREGORIAN, &
                        rc=status)
!      write(6,*) 'tcx2 ESMF_Stubs grcal ',gregoriancal%type%caltype
!      call flush(6)
      allocate( noleapCal )
!      write(6,*) 'tcx1 ESMF_Stubs nlcal ',esmf_calkind_noleap%caltype
!      call flush(6)
      noleapCal = ESMF_CalendarCreate( calkindflag=ESMF_CALKIND_NOLEAP, &
                        rc=status)
!      write(6,*) 'tcx2 ESMF_Stubs nlcal ',noleapcal%type%caltype
!      call flush(6)

      ! initialize tables in time manager
      CALL initdaym

      IF (status .ne. ESMF_SUCCESS) THEN
          PRINT *, "Error initializing the default time manager calendar"
          RETURN
      END IF
      initialized = .true.

      IF ( PRESENT( rc ) ) rc = ESMF_SUCCESS
   END SUBROUTINE ESMF_Initialize


   FUNCTION ESMF_IsInitialized()
      LOGICAL ESMF_IsInitialized
      ESMF_IsInitialized = initialized
   END FUNCTION ESMF_IsInitialized


! NOOP
   SUBROUTINE ESMF_Finalize( endflag, rc )
      USE ESMF_BaseMod
      type(ESMF_END_FLAG), intent(in), optional  :: endflag
      INTEGER, INTENT(  OUT), OPTIONAL :: rc
#ifndef HIDE_MPI
#include <mpif.h>
#endif
      INTEGER :: ier

      IF ( PRESENT( rc ) ) rc = ESMF_SUCCESS
#ifndef HIDE_MPI
      CALL MPI_Finalize( ier ) 
      IF ( ier .ne. mpi_success )THEN
        IF ( PRESENT( rc ) ) rc = ESMF_FAILURE
      END IF
#endif
   END SUBROUTINE ESMF_Finalize

! NOOP
   SUBROUTINE ESMF_LogWrite( msg, MsgType, line, file, method, log, rc )
      USE ESMF_BaseMod
      CHARACTER(LEN=*), INTENT(IN) :: msg
      TYPE(ESMF_MsgType), INTENT(IN) :: msgtype
      INTEGER, INTENT(IN), OPTIONAL :: line
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: file
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: method
      TYPE(ESMF_LOG),TARGET,OPTIONAL :: log
      INTEGER, INTENT(OUT),OPTIONAL :: rc
      IF ( PRESENT( rc ) ) rc = ESMF_SUCCESS
   END SUBROUTINE ESMF_LogWrite


END MODULE ESMF_Stubs


