!/ ------------------------------------------------------------------- /
!/    Preprocessing macros
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |        T. J. Campbell, NRL        |
!/                  |                               CPP |
!/                  | Last update :         26-Oct-2015 |
!/                  +-----------------------------------+
!/
!/    10-Dec-2014 : Origination.                     ( version 5.04 )
!/    26-Oct-2015 : Replace C style comments with Fortran
!/                  style comments.                  ( version 5.09 )
!/
!/ 1. Purpose :
!/
!/    Define preprocessor macros for WW3 ftn source code.
!/
!/ 2. Method :
!/
!/ 3. Parameters :
!/
!/ 4. Subroutines used :
!/
!/ 5. Called by :
!/
!/ 6. Error messages :
!/
!/ 7. Remarks :
!/
!/    This file uses Fortran style comments, and hence, can only be
!/    included in the Fortran (ftn) source files.  The Fortran style
!/    comments are used because not all Fortran pre-processors recognize
!/    the C style comments.
!/
!/    The __FILE__ and __LINE__ macros are defined by CPP.
!/
!/ 8. Structure :
!/
!/    See source code.
!/
!/ 9. Source code :
!/
!/ ------------------------------------------------------------------- /

!/
!/ Macros to wrap checking allocate/deallocate status
!/
#define CHECK_ALLOC_STATUS( STAT ) \
   IF ( STAT .NE. 0 ) \
   CALL EXTCDE ( 99, MSG="ALLOCATE FAILED", FILE=__FILE__, LINE=__LINE__ )
#define CHECK_DEALLOC_STATUS( STAT ) \
   IF ( STAT .NE. 0 ) \
   CALL EXTCDE ( 99, MSG="DEALLOCATE FAILED", FILE=__FILE__, LINE=__LINE__ )

!/
!/ End of w3macros.h ------------------------------------------------- /
!/
