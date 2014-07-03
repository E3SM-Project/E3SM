!BOP
! !MODULE: ice_spacecurve

module ice_spacecurve

! !DESCRIPTION:
!  This module contains routines necessary to 
!  create space-filling curves.   
!
! !REVISION HISTORY:
!  SVN:$Id: $
!
! author: John Dennis, NCAR

! !USES:
   use ice_kinds_mod
   use ice_communicate, only: my_task, master_task
   use ice_exit, only: abort_ice
   use ice_fileunits

   implicit none

! !PUBLIC TYPES: 

   type, public :: factor_t
        integer(int_kind)        :: numfact ! The # of factors for a value
        integer(int_kind), dimension(:),pointer :: factors ! The factors
        integer(int_kind), dimension(:), pointer :: used
   end type

! !PUBLIC MEMBER FUNCTIONS: 

   public :: GenSpaceCurve,     &
	     IsLoadBalanced

   public :: Factor,            &
             IsFactorable,      &
             PrintFactor,       &
             ProdFactor,        &
             MatchFactor

! !PRIVATE MEMBER FUNCTIONS:

   private :: map,    		&
	      PeanoM, 		&
	      Hilbert, 		&
	      Cinco,  		&
              GenCurve

   private :: FirstFactor,      &
              FindandMark

   integer(int_kind), dimension(:,:), allocatable ::  &
	dir,      &! direction to move along each level
        ordered    ! the ordering 
   integer(int_kind), dimension(:), allocatable ::  &
	pos        ! position along each of the axes
   
   integer(int_kind) ::  &
	maxdim,	  &! dimensionality of entire space
	vcnt       ! visitation count

   logical           :: verbose=.FALSE. 
   
   type (factor_t),  public,save :: fact  ! stores the factorization

!EOP
!EOC
!***********************************************************************

contains 

!***********************************************************************
!BOP
! !IROUTINE: Cinco
! !INTERFACE:

   recursive function Cinco(l,type,ma,md,ja,jd) result(ierr)

! !DESCRIPTION:
!  This subroutine implements a Cinco space-filling curve.
!  Cinco curves connect a Nb x Nb block of points where 
!  		
!        Nb = 5^p 
!
! !REVISION HISTORY:
!  same as module
!


! !INPUT PARAMETERS 

   integer(int_kind), intent(in) ::  &
	l, 	& ! level of the space-filling curve 
        type,   & ! type of SFC curve
	ma,     & ! Major axis [0,1]
	md,  	& ! direction of major axis [-1,1]
	ja,	& ! joiner axis [0,1]
	jd	  ! direction of joiner axis [-1,1]

! !OUTPUT PARAMETERS

   integer(int_kind) :: ierr  ! error return code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer(int_kind) :: &
	lma,		&! local major axis (next level)
	lmd,		&! local major direction (next level)
	lja,		&! local joiner axis (next level)
	ljd,		&! local joiner direction (next level)
	ltype,          &! type of SFC on next level 
        ll		 ! next level down 

   logical     :: debug = .FALSE.

!-----------------------------------------------------------------------
     ll = l
     if(ll .gt. 1) ltype = fact%factors(ll-1) ! Set the next type of space curve

     !--------------------------------------------------------------
     !  Position [0,0]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,21) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'Cinco: After Position [0,0] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [1,0]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,22) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [1,0] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [2,0]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,23) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [2,0] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [2,1]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,24) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [2,1] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [2,2]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = ma
     ljd       = -md

     if(ll .gt. 1) then
        if(debug) write(*,25) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [2,2] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [1,2]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = -md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,26) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [1,2] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [1,1]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = -md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,27) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [1,1] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [0,1]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = -md
     lja       = MOD(ma+1,maxdim)
     ljd       = md

     if(ll .gt. 1) then
        if(debug) write(*,28) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [0,1] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [0,2]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,29) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [0,2] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [0,3]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,30) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [0,3] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [0,4]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,31) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [0,4] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [1,4]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = MOD(ma+1,maxdim)
     ljd       = -md

     if(ll .gt. 1) then
        if(debug) write(*,32) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [1,4] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [1,3]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = -md
     lja       = ma
     ljd       = md

     if(ll .gt. 1) then
        if(debug) write(*,33) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [1,3] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [2,3]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,34) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [2,3] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [2,4]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,35) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [2,4] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [3,4]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,36) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [3,4] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [4,4]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = MOD(ma+1,maxdim)
     ljd       = -md

     if(ll .gt. 1) then
        if(debug) write(*,37) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [4,4] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [4,3]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = -md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,38) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [4,3] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [3,3]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = -md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,39) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [3,3] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [3,2]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = -md
     lja       = ma
     ljd       = md

     if(ll .gt. 1) then
        if(debug) write(*,40) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [3,2] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [4,2]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = MOD(ma+1,maxdim)
     ljd       = -md

     if(ll .gt. 1) then
        if(debug) write(*,41) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [4,2] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [4,1]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = -md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,42) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [4,1] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [3,1]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = -md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,43) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [3,1] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [3,0]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = -md
     lja       = ma
     ljd       = md

     if(ll .gt. 1) then
        if(debug) write(*,44) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [3,0] ',pos
     endif

     !--------------------------------------------------------------
     !  Position [4,0]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = ja
     ljd       = jd

     if(ll .gt. 1) then
        if(debug) write(*,45) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'After Position [4,0] ',pos
     endif

 21   format('Call Cinco Pos [0,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
 22   format('Call Cinco Pos [1,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
 23   format('Call Cinco Pos [2,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
 24   format('Call Cinco Pos [2,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 25   format('Call Cinco Pos [2,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 26   format('Call Cinco Pos [1,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 27   format('Call Cinco Pos [1,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 28   format('Call Cinco Pos [0,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 29   format('Call Cinco Pos [0,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 30   format('Call Cinco Pos [0,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
 31   format('Call Cinco Pos [0,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
 32   format('Call Cinco Pos [1,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
 33   format('Call Cinco Pos [1,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
 34   format('Call Cinco Pos [2,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
 35   format('Call Cinco Pos [2,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
 36   format('Call Cinco Pos [3,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
 37   format('Call Cinco Pos [4,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
 38   format('Call Cinco Pos [4,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
 39   format('Call Cinco Pos [3,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
 40   format('Call Cinco Pos [3,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 41   format('Call Cinco Pos [4,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 42   format('Call Cinco Pos [4,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 43   format('Call Cinco Pos [3,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 44   format('Call Cinco Pos [3,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
 45   format('Call Cinco Pos [4,0] Level ',i1,' at (',i2,',',i2,')',4(i3))

!EOC
!-----------------------------------------------------------------------

   end function Cinco

!***********************************************************************
!BOP
! !IROUTINE: PeanoM
! !INTERFACE:

   recursive function PeanoM(l,type,ma,md,ja,jd) result(ierr)

! !DESCRIPTION:
!  This function implements a meandering Peano 
!  space-filling curve. A meandering Peano curve 
!  connects a Nb x Nb block of points where
!
!        Nb = 3^p
!
! !REVISION HISTORY:
!  same as module
!

! !INPUT PARAMETERS

   integer(int_kind), intent(in) ::  &
        l,      & ! level of the space-filling curve
        type,   & ! type of SFC curve
        ma,     & ! Major axis [0,1]
        md,     & ! direction of major axis [-1,1]
        ja,     & ! joiner axis [0,1]
        jd        ! direction of joiner axis [-1,1]

! !OUTPUT PARAMETERS

   integer(int_kind) :: ierr  ! error return code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


   integer(int_kind) :: &
        lma,            &! local major axis (next level)
        lmd,            &! local major direction (next level)
        lja,            &! local joiner axis (next level)
        ljd,            &! local joiner direction (next level)
        ltype,          &! type of SFC on next level
        ll               ! next level down

   logical     :: debug = .FALSE.

!-----------------------------------------------------------------------

     ll = l
     if(ll .gt. 1) ltype = fact%factors(ll-1) ! Set the next type of space curve
     !--------------------------------------------------------------
     !  Position [0,0]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,21) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'PeanoM: After Position [0,0] ',pos
     endif


     !--------------------------------------------------------------
     ! Position [0,1]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd
     if(ll .gt. 1) then
        if(debug) write(*,22) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'PeanoM: After Position [0,1] ',pos
     endif

     !--------------------------------------------------------------
     ! Position [0,2]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd
     if(ll .gt. 1) then
        if(debug) write(*,23) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'PeanoM: After Position [0,2] ',pos
     endif

     !--------------------------------------------------------------
     ! Position [1,2]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd
     if(ll .gt. 1) then
        if(debug) write(*,24) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'PeanoM: After Position [1,2] ',pos
     endif


     !--------------------------------------------------------------
     ! Position [2,2]
     !--------------------------------------------------------------
     lma        = ma
     lmd        = md
     lja        = MOD(lma+1,maxdim)
     ljd        = -lmd

     if(ll .gt. 1) then
        if(debug) write(*,25) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'PeanoM: After Position [2,2] ',pos
     endif

     !--------------------------------------------------------------
     ! Position [2,1]
     !--------------------------------------------------------------
     lma        = ma
     lmd        = -md
     lja        = lma
     ljd        = lmd

     if(ll .gt. 1) then
        if(debug) write(*,26) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'PeanoM: After Position [2,1] ',pos
     endif

     !--------------------------------------------------------------
     ! Position [1,1]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = -md
     lja        = lma
     ljd        = lmd

     if(ll .gt. 1) then
        if(debug) write(*,27) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'PeanoM: After Position [1,1] ',pos
     endif


     !--------------------------------------------------------------
     ! Position [1,0]
     !--------------------------------------------------------------
     lma        = MOD(ma+1,maxdim)
     lmd        = -md
     lja        = MOD(lma+1,maxdim)
     ljd        = -lmd

     if(ll .gt. 1) then
        if(debug) write(*,28) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'PeanoM: After Position [1,0] ',pos
     endif

     !--------------------------------------------------------------
     ! Position [2,0]
     !--------------------------------------------------------------
     lma        = ma
     lmd        = md
     lja        = ja
     ljd        = jd

     if(ll .gt. 1) then
        if(debug) write(*,29) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'PeanoM: After Position [2,0] ',pos
     endif

 21   format('Call PeanoM Pos [0,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
 22   format('Call PeanoM Pos [0,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 23   format('Call PeanoM Pos [0,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 24   format('Call PeanoM Pos [1,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 25   format('Call PeanoM Pos [2,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
 26   format('Call PeanoM Pos [2,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 27   format('Call PeanoM Pos [1,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 28   format('Call PeanoM Pos [1,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
 29   format('Call PeanoM Pos [2,0] Level ',i1,' at (',i2,',',i2,')',4(i3))

!EOC
!-----------------------------------------------------------------------

   end function PeanoM

!***********************************************************************
!BOP
! !IROUTINE: Hilbert
! !INTERFACE:

   recursive function Hilbert(l,type,ma,md,ja,jd) result(ierr)

! !DESCRIPTION:
!  This function implements a Hilbert space-filling curve.
!  A Hilbert curve connect a Nb x Nb block of points where
!
!        Nb = 2^p
!
! !REVISION HISTORY:
!  same as module
!


! !INPUT PARAMETERS

   integer(int_kind), intent(in) ::  &
        l,      & ! level of the space-filling curve
        type,   & ! type of SFC curve
        ma,     & ! Major axis [0,1]
        md,     & ! direction of major axis [-1,1]
        ja,     & ! joiner axis [0,1]
        jd        ! direction of joiner axis [-1,1]

! !OUTPUT PARAMETERS

   integer(int_kind) :: ierr  ! error return code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------


   integer(int_kind) :: &
        lma,            &! local major axis (next level)
        lmd,            &! local major direction (next level)
        lja,            &! local joiner axis (next level)
        ljd,            &! local joiner direction (next level)
        ltype,          &! type of SFC on next level
        ll               ! next level down

   logical     :: debug = .FALSE.

!-----------------------------------------------------------------------
     ll = l
     if(ll .gt. 1) ltype = fact%factors(ll-1) ! Set the next type of space curve
     !--------------------------------------------------------------
     !  Position [0,0]
     !--------------------------------------------------------------
     lma       = MOD(ma+1,maxdim)
     lmd       = md
     lja       = lma
     ljd       = lmd

     if(ll .gt. 1) then
        if(debug) write(*,21) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'Hilbert: After Position [0,0] ',pos
     endif


     !--------------------------------------------------------------
     ! Position [0,1]
     !--------------------------------------------------------------
     lma       = ma
     lmd       = md
     lja       = lma
     ljd       = lmd
     if(ll .gt. 1) then
        if(debug) write(*,22) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'Hilbert: After Position [0,1] ',pos
     endif


     !--------------------------------------------------------------
     ! Position [1,1]
     !--------------------------------------------------------------
     lma        = ma
     lmd        = md
     lja        = MOD(ma+1,maxdim)
     ljd        = -md

     if(ll .gt. 1) then
        if(debug) write(*,23) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'Hilbert: After Position [1,1] ',pos
     endif

     !--------------------------------------------------------------
     ! Position [1,0]
     !--------------------------------------------------------------
     lma        = MOD(ma+1,maxdim)
     lmd        = -md
     lja        = ja
     ljd        = jd

     if(ll .gt. 1) then
        if(debug) write(*,24) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
        ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
        if(debug) call PrintCurve(ordered)
     else
        ierr = IncrementCurve(lja,ljd)
        if(debug) print *,'Hilbert: After Position [1,0] ',pos
     endif

 21   format('Call Hilbert Pos [0,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
 22   format('Call Hilbert Pos [0,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 23   format('Call Hilbert Pos [1,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
 24   format('Call Hilbert Pos [1,0] Level ',i1,' at (',i2,',',i2,')',4(i3))

!EOC
!-----------------------------------------------------------------------

   end function hilbert

!***********************************************************************
!BOP
! !IROUTINE: IncrementCurve
! !INTERFACE:

   function IncrementCurve(ja,jd) result(ierr)

! !DESCRIPTION:
!   This function creates the curve which is stored in the 
!   the ordered array.  The curve is implemented by 
!   incrementing the curve in the direction [jd] of axis [ja].
!
! !REVISION HISTORY:
!  same as module
!

! !INPUT PARAMETERS:
     integer(int_kind)  :: &
	ja, 	&! axis to increment
	jd	 ! direction along axis

! !OUTPUT PARAMETERS:
     integer(int_kind) :: ierr ! error return code

     !-----------------------------
     ! mark the newly visited point
     !-----------------------------
     ordered(pos(0)+1,pos(1)+1) = vcnt
	
     !------------------------------------
     ! increment curve and update position
     !------------------------------------
     vcnt  = vcnt + 1
     pos(ja) = pos(ja) + jd

     ierr = 0
!EOC
!-----------------------------------------------------------------------

   end function IncrementCurve

!***********************************************************************
!BOP
! !IROUTINE: log2
! !INTERFACE:

   function log2( n)

! !DESCRIPTION:
!  This function calculates the log2 of its integer 
!  input.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer(int_kind), intent(in) :: n  ! integer value to find the log2
   
! !OUTPUT PARAMETERS: 

   integer(int_kind) :: log2

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer(int_kind) ::  tmp

   !-------------------------------
   !  Find the log2 of input value
   !  Abort if n < 1
   !-------------------------------

   if (n < 1) then
      call abort_ice ('ice: spacecurve log2 error')

   elseif (n == 1) then
      log2 = 0

   else  ! n > 1
      log2 = 1
      tmp =n
      do while (tmp > 1 .and. tmp/2 .ne. 1) 
         tmp=tmp/2
         log2=log2+1
      enddo 
   endif

!EOP
!-----------------------------------------------------------------------

   end function log2

!***********************************************************************
!BOP
! !IROUTINE: IsLoadBalanced
! !INTERFACE:

   function  IsLoadBalanced(nelem,npart)
   
! !DESCRIPTION:
!  This function determines if we can create 
!  a perfectly load-balanced partitioning.
!
! !REVISION HISTORY:
!  same as module

! !INTPUT PARAMETERS:

   integer(int_kind), intent(in) ::  &
	nelem,		&  ! number of blocks/elements to partition
	npart              ! size of partition

! !OUTPUT PARAMETERS:
   logical        :: IsLoadBalanced   ! .TRUE. if a perfectly load balanced 
				      ! partition is possible
!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
	
   integer(int_kind)   :: tmp1 ! temporary int

!-----------------------------------------------------------------------
   tmp1 = nelem/npart

   if(npart*tmp1 == nelem ) then 
	IsLoadBalanced=.TRUE.
   else
        IsLoadBalanced=.FALSE.
   endif

!EOP
!-----------------------------------------------------------------------

   end function IsLoadBalanced

!***********************************************************************
!BOP
! !IROUTINE: GenCurve
! !INTERFACE:

   function GenCurve(l,type,ma,md,ja,jd) result(ierr)

! !DESCRIPTION:
!  This subroutine generates the next level down
!  space-filling curve
!
! !REVISION HISTORY:
!  same as module
!

! !INPUT PARAMETERS

   integer(int_kind), intent(in) ::  &
        l,      & ! level of the space-filling curve
        type,   & ! type of SFC curve
        ma,     & ! Major axis [0,1]
        md,     & ! direction of major axis [-1,1]
        ja,     & ! joiner axis [0,1]
        jd        ! direction of joiner axis [-1,1]

! !OUTPUT PARAMETERS

   integer(int_kind) :: ierr  ! error return code

!EOP
!BOC
!-----------------------------------------------------------------------

   !-------------------------------------------------
   ! create the space-filling curve on the next level  
   !-------------------------------------------------

   if(type == 2) then
      ierr = Hilbert(l,type,ma,md,ja,jd)
   elseif ( type == 3) then
      ierr = PeanoM(l,type,ma,md,ja,jd)
   elseif ( type == 5) then 
      ierr = Cinco(l,type,ma,md,ja,jd)
   endif

!EOP
!-----------------------------------------------------------------------

   end function GenCurve


    function FirstFactor(fac) result(res)
       type (factor_t) :: fac
       integer :: res
       logical :: found
       integer (int_kind) :: i

       found = .false.
       i=1
       do while (i<=fac%numfact .and. (.not. found))
          if(fac%used(i) == 0) then
                res = fac%factors(i)
                found = .true.
          endif
          i=i+1
        enddo

    end function FirstFactor

    function FindandMark(fac,val,f2) result(found)
       type (factor_t) :: fac
       integer :: val
       logical :: found
       logical :: f2
       integer (int_kind) :: i

       found = .false.
       i=1
       do while (i<=fac%numfact .and. (.not. found))
          if(fac%used(i) == 0) then
                if(fac%factors(i) .eq. val) then
                   if(f2)  then
                      fac%used(i) = 1
                      found = .true.
                   else if( .not. f2) then
                      fac%used(i) = -1
                      found = .false.
                   endif
                endif
          endif
          i=i+1
        enddo

    end function FindandMark


   subroutine MatchFactor(fac1,fac2,val,found)
      type (factor_t) :: fac1
      type (factor_t) :: fac2
      integer :: val
      integer :: val1
      logical :: found
      logical :: tmp

      found = .false.

      val1 = FirstFactor(fac1)
!JMD      print *,'Matchfactor: found value: ',val1
      found = FindandMark(fac2,val1,.true.)
      tmp = FindandMark(fac1,val1,found)
      if (found) then
        val = val1
      else
        val = 1
      endif

   end subroutine MatchFactor

   function ProdFactor(fac) result(res)

   type (factor_t) :: fac
   integer :: res
   integer (int_kind) :: i

     res = 1
     do i=1,fac%numfact
        if(fac%used(i) <= 0) then
          res = res * fac%factors(i)
        endif
     enddo

   end function ProdFactor

   subroutine PrintFactor(msg,fac)


      character(len=*) :: msg
      type (factor_t) :: fac
      integer (int_kind) :: i

      write(*,*) ' '
      write(*,*) 'PrintFactor: ',msg
      write(*,*) (fac%factors(i),i=1,fac%numfact)
      write(*,*) (fac%used(i),i=1,fac%numfact)


   end subroutine PrintFactor

!***********************************************************************
!BOP
! !IROUTINE: Factor
! !INTERFACE:

   function Factor(num) result(res)

! !DESCRIPTION:
!  This function factors the input value num into a 
!  product of 2,3, and 5.
!
! !REVISION HISTORY:
!  same as module
!

! !INPUT PARAMETERS:

   integer(int_kind), intent(in)  :: num  ! number to factor

! !OUTPUT PARAMETERS:

   type (factor_t)     :: res

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer(int_kind)   ::  &
	tmp,tmp2,tmp3,tmp5   ! tempories for the factorization algorithm
   integer(int_kind)   :: i,n    ! loop tempories
   logical             :: found  ! logical temporary

   ! --------------------------------------
   ! Allocate allocate for max # of factors
   ! --------------------------------------
   tmp = num
   tmp2 = max(log2(num),1)
   allocate(res%factors(tmp2))
   allocate(res%used(tmp2))

   res%used = 0
   n=0


   !-----------------------
   !  Look for factors of 2
   !-----------------------
   found=.TRUE.
   do while (found)
      found = .FALSE.
      tmp2 = tmp/2
      if( tmp2*2 == tmp ) then
        n = n + 1
        res%factors(n) = 2
        found = .TRUE.
        tmp = tmp2
      endif
   enddo

   !-----------------------
   !  Look for factors of 3
   !-----------------------
   found=.TRUE.
   do while (found)
      found = .FALSE.
      tmp3 = tmp/3
      if( tmp3*3 == tmp ) then
        n = n + 1
        res%factors(n) = 3
        found = .TRUE.
        tmp = tmp3
      endif
   enddo

   !-----------------------
   !  Look for factors of 5
   !-----------------------
   found=.TRUE.
   do while (found)
      found = .FALSE.
      tmp5 = tmp/5
      if( tmp5*5 == tmp ) then
        n = n + 1
        res%factors(n) = 5
        found = .TRUE.
        tmp = tmp5
      endif
   enddo

   !------------------------------------
   ! make sure that the input value 
   ! only contains factors of 2,3,and 5  
   !------------------------------------
   tmp=1
   do i=1,n
     tmp = tmp * res%factors(i)
   enddo
   if(tmp == num) then
     res%numfact = n
   else
     res%numfact = -1
   endif

!EOP
!---------------------------------------------------------
   end function Factor

!***********************************************************************
!BOP
! !IROUTINE: IsFactorable
! !INTERFACE:

   function IsFactorable(n)
   
! !DESCRIPTION:
!  This function determines if we can factor
!   n into 2,3,and 5.  
!
! !REVISION HISTORY:
!  same as module


! !INTPUT PARAMETERS:

   integer(int_kind), intent(in)  :: n  ! number to factor

! !OUTPUT PARAMETERS:
   logical  :: IsFactorable  ! .TRUE. if it is factorable

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   type (factor_t)     :: fact  ! data structure to store factor information

   fact = Factor(n)
   if(fact%numfact .ne. -1) then
     IsFactorable = .TRUE.
   else
     IsFactorable = .FALSE.
   endif

!EOP
!-----------------------------------------------------------------------

   end function IsFactorable

!***********************************************************************
!BOP
! !IROUTINE: map
! !INTERFACE:

   subroutine map(l)

! !DESCRIPTION:
!   Interface routine between internal subroutines and public 
!   subroutines.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:
   integer(int_kind)  :: l   ! level of space-filling curve


!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer(int_kind)  :: &
	d, 		 & ! dimension of curve only 2D is supported
	type,		 & ! type of space-filling curve to start off
        ierr   		   ! error return code

   d = SIZE(pos)

   pos=0
   maxdim=d
   vcnt=0

   type = fact%factors(l)
   ierr = GenCurve(l,type,0,1,0,1)


!EOP
!-----------------------------------------------------------------------

   end subroutine map

!***********************************************************************
!BOP
! !IROUTINE: PrintCurve
! !INTERFACE:

   subroutine PrintCurve(Mesh)


! !DESCRIPTION:
!  This subroutine prints the several low order 
!  space-filling curves in an easy to read format
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:

     integer(int_kind), intent(in), target ::  Mesh(:,:) ! SFC to be printed

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
     integer(int_kind) ::  &
        gridsize,	   &! order of space-filling curve
        i		    ! loop temporary

!-----------------------------------------------------------------------

     gridsize = SIZE(Mesh,dim=1)

     if(gridsize == 2) then
        write (*,*) "A Level 1 Hilbert Curve:"
        write (*,*) "------------------------"
        do i=1,gridsize
           write(*,2) Mesh(1,i),Mesh(2,i)
        enddo
     else if(gridsize == 3) then
        write (*,*) "A Level 1 Peano Meandering Curve:"
        write (*,*) "---------------------------------"
        do i=1,gridsize
           write(*,3) Mesh(1,i),Mesh(2,i),Mesh(3,i)
        enddo
     else if(gridsize == 4) then
        write (*,*) "A Level 2 Hilbert Curve:"
        write (*,*) "------------------------"
        do i=1,gridsize
           write(*,4) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i)
        enddo
     else if(gridsize == 5) then
        write (*,*) "A Level 1 Cinco Curve:"
        write (*,*) "------------------------"
        do i=1,gridsize
           write(*,5) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i),Mesh(5,i)
        enddo
     else if(gridsize == 6) then
        write (*,*) "A Level 1 Hilbert and Level 1 Peano Curve:"
        write (*,*) "------------------------------------------"
        do i=1,gridsize
           write(*,6) Mesh(1,i),Mesh(2,i),Mesh(3,i), &
	    	      Mesh(4,i),Mesh(5,i),Mesh(6,i)
        enddo
     else if(gridsize == 8) then
        write (*,*) "A Level 3 Hilbert Curve:"
        write (*,*) "------------------------"
        do i=1,gridsize
           write(*,8) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                      Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i)
         enddo
     else if(gridsize == 9) then
        write (*,*) "A Level 2 Peano Meandering Curve:"
        write (*,*) "---------------------------------"
        do i=1,gridsize
           write(*,9) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                      Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i), &
                      Mesh(9,i)
         enddo
     else if(gridsize == 10) then
        write (*,*) "A Level 1 Hilbert and Level 1 Cinco Curve:"
        write (*,*) "---------------------------------"
        do i=1,gridsize
           write(*,10) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                      Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i), &
                      Mesh(9,i),Mesh(10,i)
         enddo
     else if(gridsize == 12) then
        write (*,*) "A Level 2 Hilbert and Level 1 Peano Curve:"
        write (*,*) "------------------------------------------"
        do i=1,gridsize
           write(*,12) Mesh(1,i),Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                       Mesh(5,i),Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                       Mesh(9,i),Mesh(10,i),Mesh(11,i),Mesh(12,i)
        enddo
     else if(gridsize == 15) then
        write (*,*) "A Level 1 Peano and Level 1 Cinco Curve:"
        write (*,*) "------------------------"
        do i=1,gridsize
           write(*,15) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                       Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i), &
                       Mesh(9,i),Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i)
        enddo
     else if(gridsize == 16) then
        write (*,*) "A Level 4 Hilbert Curve:"
        write (*,*) "------------------------"
        do i=1,gridsize
           write(*,16) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                       Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i), &
                       Mesh(9,i),Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i)
        enddo
     else if(gridsize == 18) then
        write (*,*) "A Level 1 Hilbert and Level 2 Peano Curve:"
        write (*,*) "------------------------------------------"
        do i=1,gridsize
           write(*,18) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                       Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                       Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                       Mesh(17,i),Mesh(18,i)
        enddo
     else if(gridsize == 20) then
        write (*,*) "A Level 2 Hilbert and Level 1 Cinco Curve:"
        write (*,*) "------------------------------------------"
        do i=1,gridsize
           write(*,20) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                       Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                       Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                       Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i)
        enddo
     else if(gridsize == 24) then
        write (*,*) "A Level 3 Hilbert and Level 1 Peano Curve:"
        write (*,*) "------------------------------------------"
        do i=1,gridsize
           write(*,24) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                       Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                       Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                       Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i), &
                       Mesh(21,i),Mesh(22,i),Mesh(23,i),Mesh(24,i)
        enddo
     else if(gridsize == 25) then
        write (*,*) "A Level 2 Cinco Curve:"
        write (*,*) "------------------------------------------"
        do i=1,gridsize
           write(*,25) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                       Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                       Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                       Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i), &
                       Mesh(21,i),Mesh(22,i),Mesh(23,i),Mesh(24,i), &
		       Mesh(25,i)
        enddo
     else if(gridsize == 27) then
        write (*,*) "A Level 3 Peano Meandering Curve:"
        write (*,*) "---------------------------------"
        do i=1,gridsize
           write(*,27) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                       Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                       Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                       Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i), &
                       Mesh(21,i),Mesh(22,i),Mesh(23,i),Mesh(24,i), &
                       Mesh(25,i),Mesh(26,i),Mesh(27,i)
        enddo
     else if(gridsize == 32) then
        write (*,*) "A Level 5 Hilbert Curve:"
        write (*,*) "------------------------"
        do i=1,gridsize
           write(*,32) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i),  &
                       Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i),  &
                       Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                       Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                       Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i), &
                       Mesh(21,i),Mesh(22,i),Mesh(23,i),Mesh(24,i), &
                       Mesh(25,i),Mesh(26,i),Mesh(27,i),Mesh(28,i), &
                       Mesh(29,i),Mesh(30,i),Mesh(31,i),Mesh(32,i)
        enddo
     endif
 2 format('|',2(i2,'|'))
 3 format('|',3(i2,'|'))
 4 format('|',4(i2,'|'))
 5 format('|',5(i2,'|'))
 6 format('|',6(i2,'|'))
 8 format('|',8(i2,'|'))
 9 format('|',9(i2,'|'))
10 format('|',10(i2,'|'))
12 format('|',12(i3,'|'))
15 format('|',15(i3,'|'))
16 format('|',16(i3,'|'))
18 format('|',18(i3,'|'))
20 format('|',20(i3,'|'))
24 format('|',24(i3,'|'))
25 format('|',25(i3,'|'))
27 format('|',27(i3,'|'))
32 format('|',32(i4,'|'))

!EOC
!-----------------------------------------------------------------------

   end subroutine PrintCurve

!***********************************************************************
!BOP
! !IROUTINE: GenSpaceCurve
! !INTERFACE:

  subroutine  GenSpaceCurve(Mesh)

! !DESCRIPTION:
!  This subroutine is the public interface into the 
!  space-filling curve functionality
!
! !REVISION HISTORY:
!  same as module
!

! !INPUT/OUTPUT PARAMETERS:
   integer(int_kind), target,intent(inout) :: &
	Mesh(:,:)		! The SFC ordering in 2D array

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer(int_kind) ::  &
	level,   &! Level of space-filling curve		
	dim       ! dimension of SFC... currently limited to 2D

   integer(int_kind) :: gridsize   ! number of points on a side
   
!-----------------------------------------------------------------------

   !-----------------------------------------
   !  Setup the size of the grid to traverse
   !-----------------------------------------
   dim = 2
   gridsize = SIZE(Mesh,dim=1)
   fact     = factor(gridsize)
   level    = fact%numfact

   if(verbose) print *,'GenSpacecurve: level is ',level
   allocate(ordered(gridsize,gridsize))

   !--------------------------------------------
   ! Setup the working arrays for the traversal
   !--------------------------------------------
   allocate(pos(0:dim-1))
   
   !-----------------------------------------------------
   !  The array ordered will contain the visitation order
   !-----------------------------------------------------
   ordered(:,:) = 0

   call map(level) 

   Mesh(:,:) = ordered(:,:)

   deallocate(pos,ordered)

!EOP
!-----------------------------------------------------------------------

  end subroutine GenSpaceCurve 

  recursive subroutine qsort(a)
   
    integer, intent(inout) :: a(:)
    integer :: split
   
    if(SIZE(a) > 1) then 
      call partition(a,split)
      call qsort(a(:split-1))
      call qsort(a(split:))
    endif 

  end subroutine qsort

  subroutine partition(a,marker)

     INTEGER, INTENT(IN OUT) :: a(:)
     INTEGER, INTENT(OUT) :: marker
     INTEGER :: left, right, pivot, temp
 
     pivot = (a(1) + a(size(a))) / 2  ! Average of first and last elements to prevent quadratic 
     left = 0                         ! behavior with sorted or reverse sorted data
     right = size(a) + 1
 
     DO WHILE (left < right)
        right = right - 1
        DO WHILE (a(right) > pivot)
           right = right-1
        END DO
        left = left + 1
        DO WHILE (a(left) < pivot)
           left = left + 1
        END DO
        IF (left < right) THEN 
           temp = a(left)
           a(left) = a(right)
           a(right) = temp
        END IF
     END DO
 
     IF (left == right) THEN
        marker = left + 1
     ELSE
        marker = left
     END IF

  end subroutine partition
     


end module ice_spacecurve

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
