! FORTRAN interfaces to the C utilties defined in fortlib.c
!
! Written by: Richard Weed, Ph.D
!             Center for Advanced Vehicular Systems
!             Mississippi State University
!             rweed@cavs.msstate.edu




! License (and other Lawyer Language)
 
! This software is released under the Apache 2.0 Open Source License. The
! full text of the License can be viewed at :
!
!   http:www.apache.org/licenses/LICENSE-2.0.html
!
! The author grants to the University Corporation for Atmospheric Research
! (UCAR), Boulder, CO, USA the right to revise and extend the software
! without restriction. However, the author retains all copyrights and
! intellectual property rights explicitly stated in or implied by the
! Apache license

! Version 1.  June 2006
!             Unchanged for netCDF 4.1.1

!-------------------------------  udexit --------------------------------------
 Subroutine udexit(status)
!
 USE ISO_C_BINDING, ONLY: C_INT
 Implicit NONE

 
 Integer, Intent(IN) :: status
 Integer(KIND=C_INT) :: cstatus
 Interface
  Subroutine exit(status) BIND(C)
  USE ISO_C_BINDING, ONLY: C_INT
  Integer(KIND=C_INT), VALUE :: status  
  End Subroutine exit
 End Interface

 cstatus = status
 Call exit(cstatus)

 End Subroutine udexit
!-------------------------------  udabort --------------------------------------
 Subroutine udabort()
 USE ISO_C_BINDING
 Implicit NONE
 Interface
  Subroutine abort() BIND(C)
  End Subroutine abort
 End Interface

 Call abort()

 End Subroutine udabort
!-------------------------------  udrand -------------------------------------
 Function udrand(iflag) RESULT(rannum)

 USE ISO_C_BINDING

 Implicit NONE
 Integer, Parameter :: RK8=SELECTED_REAL_KIND(P=13, R=307)  ! double
 Integer, Intent(IN) :: iflag
 Real(RK8) :: rannum
 Integer(KIND=C_INT) :: ciflag
 Real(KIND=C_DOUBLE) :: crannum
 Interface
  Function myrand(iflag) BIND(C)
   USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
   Integer(KIND=C_INT), VALUE :: iflag
   Real(KIND=C_DOUBLE) :: myrand
  End Function myrand
 End Interface

 ciflag = iflag
  
 crannum = myrand(ciflag)
 rannum = crannum

 End Function udrand
!-------------------------------  udshift -------------------------------------
 Function udshift(ivalue, amount) RESULT(shiftval)

 USE ISO_C_BINDING

 Implicit NONE

 Integer, Intent(IN) :: ivalue, amount
 Integer(KIND=C_INT) :: cvalue, camount 
 Integer :: shiftval 
 Integer(KIND=C_INT) :: cshiftval
 Interface
  Function myshift(cvalue, camount) BIND(C)
   USE ISO_C_BINDING, ONLY: C_INT
   Integer(KIND=C_INT), VALUE :: cvalue, camount 
   Integer(KIND=C_INT) :: myshift
  End Function myshift
 End Interface

 cvalue = ivalue
 camount = amount
  
 cshiftval = myshift(cvalue, camount)
 shiftval = cshiftval

 End Function udshift
!-------------------------------  ignorefpe  ----------------------------------
 Subroutine ignorefpe(idoit)

 USE ISO_C_BINDING

 Implicit NONE

 Integer, Intent(IN) :: idoit 
 Integer(KIND=C_INT) :: cdoit 
 Interface
  Subroutine nc_ignorefpe(cdoit) BIND(C)
   USE ISO_C_BINDING, ONLY: C_INT
   Integer(KIND=C_INT), VALUE :: cdoit 
  End Subroutine nc_ignorefpe 
 End Interface

 cdoit = idoit
 Call nc_ignorefpe(cdoit)

 End Subroutine ignorefpe
!-------------------------------  max_uchar  ----------------------------------
 Function max_uchar() RESULT(cmax)

 USE ISO_C_BINDING

 Implicit NONE
 Integer, Parameter :: RK8=SELECTED_REAL_KIND(P=13, R=307)  ! double
 Real(RK8) :: cmax
 Real(KIND=C_DOUBLE) :: ccmax
 Interface
  Function cmax_uchar() BIND(C)
   USE ISO_C_BINDING, ONLY: C_DOUBLE
   Real(KIND=C_DOUBLE) :: cmax_uchar
  End Function cmax_uchar
 End Interface

 ccmax = cmax_uchar()
 cmax = ccmax

 End Function max_uchar
!-------------------------------  min_schar  ----------------------------------
 Function min_schar() RESULT(cmin)

 USE ISO_C_BINDING

 Implicit NONE
 Integer, Parameter :: RK8=SELECTED_REAL_KIND(P=13, R=307)  ! double
 Real(RK8) :: cmin
 Real(KIND=C_DOUBLE) :: ccmin
 Interface
  Function cmin_schar() BIND(C)
   USE ISO_C_BINDING, ONLY: C_DOUBLE
   Real(KIND=C_DOUBLE) :: cmin_schar
  End Function cmin_schar
 End Interface

 ccmin = cmin_schar()
 cmin = ccmin

 End Function min_schar
!-------------------------------  max_schar  ----------------------------------
 Function max_schar() RESULT(cmax)

 USE ISO_C_BINDING

 Implicit NONE
 Integer, Parameter :: RK8=SELECTED_REAL_KIND(P=13, R=307)  ! double
 Real(RK8) :: cmax
 Real(KIND=C_DOUBLE) :: ccmax
 Interface
  Function cmax_schar() BIND(C)
   USE ISO_C_BINDING, ONLY: C_DOUBLE
   Real(KIND=C_DOUBLE) :: cmax_schar
  End Function cmax_schar
 End Interface

 ccmax = cmax_schar()
 cmax = ccmax

 End Function max_schar
!-------------------------------  min_short  ----------------------------------
 Function min_short() RESULT(cmin)

 USE ISO_C_BINDING

 Implicit NONE
 Integer, Parameter :: RK8=SELECTED_REAL_KIND(P=13, R=307)  ! double
 Real(RK8) :: cmin
 Real(KIND=C_DOUBLE) :: ccmin
 Interface
  Function cmin_short() BIND(C)
   USE ISO_C_BINDING, ONLY: C_DOUBLE
   Real(KIND=C_DOUBLE) :: cmin_short
  End Function cmin_short
 End Interface

 ccmin = cmin_short()
 cmin = ccmin

 End Function min_short
!-------------------------------  max_short  ----------------------------------
 Function max_short() RESULT(cmax)

 USE ISO_C_BINDING

 Implicit NONE
 Integer, Parameter :: RK8=SELECTED_REAL_KIND(P=13, R=307)  ! double
 Real(RK8) :: cmax
 Real(KIND=C_DOUBLE) :: ccmax
 Interface
  Function cmax_short() BIND(C)
   USE ISO_C_BINDING, ONLY: C_DOUBLE
   Real(KIND=C_DOUBLE) :: cmax_short
  End Function cmax_short
 End Interface

 ccmax = cmax_short()
 cmax = ccmax

 End Function max_short
!-------------------------------  min_int  ----------------------------------
 Function min_int() RESULT(cmin)

 USE ISO_C_BINDING

 Implicit NONE
 Integer, Parameter :: RK8=SELECTED_REAL_KIND(P=13, R=307)  ! double
 Real(RK8) :: cmin
 Real(KIND=C_DOUBLE) :: ccmin
 Interface
  Function cmin_int() BIND(C)
   USE ISO_C_BINDING, ONLY: C_DOUBLE
   Real(KIND=C_DOUBLE) :: cmin_int
  End Function cmin_int
 End Interface

 ccmin = cmin_int()
 cmin = ccmin

 End Function min_int
!-------------------------------  max_int  ----------------------------------
 Function max_int() RESULT(cmax)

 USE ISO_C_BINDING

 Implicit NONE
 Integer, Parameter :: RK8=SELECTED_REAL_KIND(P=13, R=307)  ! double
 Real(RK8) :: cmax
 Real(KIND=C_DOUBLE) :: ccmax
 Interface
  Function cmax_int() BIND(C)
   USE ISO_C_BINDING, ONLY: C_DOUBLE
   Real(KIND=C_DOUBLE) :: cmax_int
  End Function cmax_int
 End Interface

 ccmax = cmax_int()
 cmax = ccmax

 End Function max_int
!-------------------------------  min_long  ----------------------------------
 Function min_long() RESULT(cmin)

 USE ISO_C_BINDING

 Implicit NONE
 Integer, Parameter :: RK8=SELECTED_REAL_KIND(P=13, R=307)  ! double
 Real(RK8) :: cmin
 Real(KIND=C_DOUBLE) :: ccmin
 Interface
  Function cmin_long() BIND(C)
   USE ISO_C_BINDING, ONLY: C_DOUBLE
   Real(KIND=C_DOUBLE) :: cmin_long
  End Function cmin_long
 End Interface

 ccmin = cmin_long()
 cmin = ccmin

 End Function min_long
!-------------------------------  max_long  ----------------------------------
 Function max_long() RESULT(cmax)

 USE ISO_C_BINDING

 Implicit NONE
 Integer, Parameter :: RK8=SELECTED_REAL_KIND(P=13, R=307)  ! double
 Real(RK8) :: cmax
 Real(KIND=C_DOUBLE) :: ccmax
 Interface
  Function cmax_long() BIND(C)
   USE ISO_C_BINDING, ONLY: C_DOUBLE
   Real(KIND=C_DOUBLE) :: cmax_long
  End Function cmax_long
 End Interface

 ccmax = cmax_long()
 cmax = ccmax

 End Function max_long
!-------------------------------  max_float  ----------------------------------
 Function max_float() RESULT(cmax)

 USE ISO_C_BINDING

 Implicit NONE
 Integer, Parameter :: RK8=SELECTED_REAL_KIND(P=13, R=307)  ! double
 Real(RK8) :: cmax
 Real(KIND=C_DOUBLE) :: ccmax
 Interface
  Function cmax_float() BIND(C)
   USE ISO_C_BINDING, ONLY: C_DOUBLE
   Real(KIND=C_DOUBLE) :: cmax_float
  End Function cmax_float
 End Interface

 ccmax = cmax_float()
 cmax = ccmax

 End Function max_float
!-------------------------------  max_double  ----------------------------------
 Function max_double() RESULT(cmax)

 USE ISO_C_BINDING

 Implicit NONE
 Integer, Parameter :: RK8=SELECTED_REAL_KIND(P=13, R=307)  ! double
 Real(RK8) :: cmax
 Real(KIND=C_DOUBLE) :: ccmax
 Interface
  Function cmax_double() BIND(C)
   USE ISO_C_BINDING, ONLY: C_DOUBLE
   Real(KIND=C_DOUBLE) :: cmax_double
  End Function cmax_double
 End Interface

 ccmax = cmax_double()
 cmax = ccmax

 End Function max_double
