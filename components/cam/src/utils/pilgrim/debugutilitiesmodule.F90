!-------------------------------------------------------------------------
!         Nasa/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS
!-------------------------------------------------------------------------
      MODULE debugutilitiesmodule
!BOP
!
! !MODULE: debugutilitiesmodule
!
! !USES:
#if defined( STAND_ALONE )
# define iulog 6
#else
      use cam_logfile, only: iulog
#endif
      IMPLICIT NONE

#define  MAX_STACK_LEVEL 20
#define  MAX_STRING_LEN  40
!
! !PUBLIC MEMBER FUNCTIONS:
      PUBLIC     DumAssert, DumEnter, DumLeave

!
! !DESCRIPTION:
!
!      This module provides the basic utilities to support debugging
!
!      \begin{tabular}{|l|l|} \hline \hline
!        DumAssert         & Make an assertion \\ \hline
!        DumEnter          & Tracing: enter a subroutine \\ \hline
!        DumLeave          & Tracing: leave a subroutine   \\ \hline
!      \end{tabular}
!
!      The DumAssert makes an assertion (i.e., claims that a boolean
!      argument is true) for a given line of code in a given source
!      file.  DumEnter and DumLeave to be used as a pair and placed at the
!      beginning and end of routines to be traced.
!
!      It is not intended for the user to make use of these routines
!      directly but rather in conjunction with the CPP macros defined
!      in the "Debug.h" file in the INCLUDE directory.  The CPP 
!      macros define the calls to the three above-mention routines if
!      the -DDEBUG\_ON option is set on the compile line.  The line
!      \#include "Debug.h" statement in any routine which makes use 
!      of these facilities.  In production compilations where DEBUG\_ON
!      is not set, "Debug.h" defines blank lines and thus does not
!      affect code performance.  The CPP definition of DEBUG\_LEVEL in
!      the compile line, e.g., -DDEBUG\_LEVEL=2, denotes the level
!      of debugging performed.  A higher level performs all the 
!      debugging at the lower levels and then some.
!
!      Note that, unlike other include statements, "Debug.h" must be 
!      included {\it before} the IMPLICIT NONE statement (since "Debug.h"
!      contains a USE DebugModule statement.
!
!      Compile options used:  {\tt MPI\_VER}, {\tt DEBUG\_LEVEL}
!
! !LOCAL VARIABLES:
      CHARACTER(len=MAX_STRING_LEN) :: TraceStack( MAX_STACK_LEVEL )
      INTEGER    :: StackLevel = 0
!
! !REVISION HISTORY:
!   97.09.30   Sawyer     Creation
!   98.03.09   Sawyer     Added documentation for walkthrough
!   01.02.12   Sawyer     Converted to free format
!
! !BUGS:
!
!EOP
      CONTAINS
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: DumAssert --- Raise Assertion
!
! !INTERFACE: 
      SUBROUTINE DumAssert ( Condition, FileName, Linenumber )
!
! !USES:
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      LOGICAL, INTENT(IN)       :: Condition      ! Condition asserted
      CHARACTER(*), INTENT(IN)  :: FileName       ! Source file
      INTEGER, INTENT( IN )     :: LineNumber     ! Source line

! !DESCRIPTION:
!     Condition is claimed by the calling Routine in Filename at
!     Linenumber to be true.  If it is, do nothing.  If not, print
!     as much information as possible.  
!
!      \begin{tabular}{|c|l|} \hline \hline
!        {\bf Debug Level} & {\bf Action} \\ \hline \hline
!        0                 & Return immediately \\ \hline
!        1                 & Print assertion failed \\ \hline
!        2                 & Print assertion failed and trace stack \\ \hline
!      \end{tabular}
!
! !LOCAL VARIABLES:
      INTEGER I, MyID, Ierror
!
! !SYSTEM ROUTINES:
!
! !REVISION HISTORY:
!   97.09.30   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC

#if !defined(DEBUG_LEVEL)
#define DEBUG_LEVEL 1
#endif

#if ( DEBUG_LEVEL > 0 )
      IF (.NOT. Condition) THEN
        write(iulog,*) 'Assertion failed:',                                    &
                  ' source file: ', FileName,                            &
                  ' source line: ', LineNumber
!
!  Check if trace available
!
#if ( DEBUG_LEVEL > 1 )
        PRINT *, "Printing Trace: "
        DO I = 1, StackLevel
          PRINT *, "Level ", StackLevel,                                 &
                   " Called ", TraceStack( StackLevel )
        ENDDO
#endif
      ENDIF
#endif
      RETURN
!EOC
      END SUBROUTINE DumAssert
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: DumEnter --- Tracing: Enter a Subroutine
!
! !INTERFACE: 
      SUBROUTINE DumEnter ( RoutineName )
!
! !USES:
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      CHARACTER(*), INTENT(IN)     :: RoutineName    ! Source file

! !DESCRIPTION:
!      This routine marks the beginning of a region to be traced,
!      usually a subroutine.  
!      
!      \begin{tabular}{|c|l|} \hline \hline
!        {\bf Debug Level} & {\bf Action} \\ \hline \hline
!        0                 & Return immediately \\ \hline
!        1                 & Perform bookkeeping \\ \hline
!        2                 & Perform bookkeeping, print trace \\ \hline
!      \end{tabular}
!
! !LOCAL VARIABLES:
      INTEGER MyID, Ierror
! !REVISION HISTORY:
!   97.09.30   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
#if !defined(DEBUG_LEVEL)
#define DEBUG_LEVEL 1
#endif

#if ( DEBUG_LEVEL > 0 ) 
      StackLevel = StackLevel + 1
      IF ( StackLevel .GT. MAX_STACK_LEVEL ) THEN
        PRINT *, "StackLevel overflow: ", StackLevel, " Stopping"
        STOP
      ENDIF
      TraceStack( StackLevel ) = RoutineName      
#if ( DEBUG_LEVEL > 1 )
      PRINT *, "Level ", StackLevel, " Entering ", RoutineName
#endif
#endif
      RETURN
!EOC
      END SUBROUTINE DumEnter
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: DumLeave --- Tracing: Leave a Subroutine
!
! !INTERFACE: 
      SUBROUTINE DumLeave ( RoutineName )
!
! !USES:
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
      CHARACTER(*), INTENT(IN)     :: RoutineName    ! Source file

! !DESCRIPTION:
!     Tracing facility: leave a subroutine, remove the history trail.
!     Depending on the debugging level, do nothing (0), update the
!     stack only (1), or update stack and print trace message (2) to
!     stdout.  The CALL to Leave should be placed just before every
!     egress of the subroutine (hopefully the exit point is unique).
!
!      \begin{tabular}{|c|l|} \hline \hline
!        {\bf Debug Level} & {\bf Action} \\ \hline \hline
!        0                 & Return immediately \\ \hline
!        1                 & Perform bookkeeping, consistency check \\ \hline
!        2                 & Bookkeeping, consistency, print trace \\ \hline
!      \end{tabular}
!
! !LOCAL VARIABLES:
      INTEGER MyID, Ierror
! !REVISION HISTORY:
!   97.09.30   Sawyer     Creation
!
!EOP
!-----------------------------------------------------------------------
!BOC
#if !defined(DEBUG_LEVEL)
#define DEBUG_LEVEL 1
#endif

#if ( DEBUG_LEVEL > 0 ) 
!
!  Make sure that the Enter and Leave correspond
!
      IF ( TraceStack(StackLevel) .NE. RoutineName ) THEN
        PRINT *, "Expected: ", TraceStack(StackLevel),                   &
     &           "Got: ", RoutineName, " STOPPING "
        STOP 
      ENDIF
#if ( DEBUG_LEVEL > 1 )
      PRINT *, "Level ", StackLevel, " Leaving ", RoutineName
#endif
      IF ( StackLevel .LE. 0 ) THEN
        PRINT *, "StackLevel underflow: ", StackLevel, " Stopping"
        STOP
      ENDIF
      TraceStack( StackLevel ) = ""
      StackLevel = StackLevel - 1
#endif
      RETURN
!EOC
      END SUBROUTINE DumLeave
!-----------------------------------------------------------------------

      END MODULE debugutilitiesmodule
