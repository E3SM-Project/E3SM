!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !INCLUDE: Debug.h --- Debugging Macros
!
! !DESCRIPTION: 
!
!      The following CPP macros define the calls to the F90ass, Enter
!      and Leave in DebugUtilitiesModule.F if the -DDEBUG\_ON option
!      is set on the compile line.  The line \#include <Debug.h> makes
!      use of these facilities.  In production compilations where
!      DEBUG\_ON is not set, Debug.h defines blank lines and thus
!      does not affect code performance.
!
!      Note that, unlike other include statements, "Debug.h" must be
!      included {\it before} the IMPLICIT NONE statement (since
!      "Debug.h" contains a USE DebugUtilitiesModule statement).
!
!      Compile options used:  {\tt DEBUG\_ON}
!
! !SEE ALSO: 
!
!  DebugUtilitiesModule.F   -   Debugging module (implementation)
!
! !REVISION HISTORY: 
!
!  97.09.30  Sawyer   Creation
!  98.03.09  Sawyer   Renamed Debug.h, prepared for Walkthrough
!  98.09.03  Sawyer   DEBUG_ON does not work with f90 -cpp on SGI
!
!EOP
!-------------------------------------------------------------------------
!BOC
#if defined( DEBUG_ON ) && !defined( IRIX64 )
!
! This include file is unusual in that it MUST BE PLACED BEFORE 
! IMPLICIT NONE, not it the usual location for include files
!
      USE DebugUtilitiesModule, ONLY: DumAssert, DumEnter, DumLeave
#endif

#if defined( DEBUG_ON ) && !defined( IRIX64 )
#define CPP_ASSERT_INFO(cond_) PRINT *, "Assert: ", cond_," l:",  __LINE__
#define CPP_ASSERT_F90(cond_) CALL DumAssert(cond_,__FILE__,__LINE__)
#else
#define CPP_ASSERT_INFO(cond_) ! PRINT *, "Assert: ", Excised Test," l:",  __LINE__
#define CPP_ASSERT_F90(cond_) ! CALL DumAssert(Excised Test,__FILE__,__LINE__)
#endif

#if defined( DEBUG_ON ) && !defined( IRIX64 )
#define CPP_ENTER_PROCEDURE(ROUTINE_NAME) CALL DumEnter( ROUTINE_NAME )
#define CPP_LEAVE_PROCEDURE(ROUTINE_NAME) CALL DumLeave( ROUTINE_NAME )
#else
#define CPP_ENTER_PROCEDURE(ROUTINE_NAME) ! CALL DumEnter( ROUTINE_NAME )
#define CPP_LEAVE_PROCEDURE(ROUTINE_NAME) ! CALL DumLeave( ROUTINE_NAME )
#endif
!EOC
!-------------------------------------------------------------------------
