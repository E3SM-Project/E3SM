#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module spacecurve_mod
  use kinds, only : iulog
  implicit none
  private

  type, public :: factor_t
     integer                       :: numfact
     integer, dimension(:),pointer :: factors
  end type factor_t


  integer,public, dimension(:,:), allocatable :: ordered
  integer,public, dimension(:,:), allocatable :: dir     ! direction to move along each level
  integer,public, dimension(:)  , allocatable :: pos     ! position along each of the axes

  integer,public                              :: maxdim  ! dimensionality of entire space
  integer,public                              :: vcnt   ! visitation count
  logical,private                             :: verbose=.FALSE. 

  type (factor_t),  public                      :: fact

  !JMD new addition
  SAVE:: fact
  public :: map 
  public :: hilbert_old
  public :: PeanoM,hilbert, Cinco
  public :: GenCurve
  public :: GenSpaceCurve
  public :: log2,Factor
  public :: PrintCurve
  public :: IsFactorable,IsLoadBalanced
  public :: genspacepart
contains 
  !---------------------------------------------------------
  recursive function Cinco(l,type,ma,md,ja,jd) result(ierr)

    implicit none
    integer,intent(in) ::   l,type,ma,md,ja,jd

    integer ::  lma,lmd,lja,ljd,ltype
    integer :: ll
    integer :: ierr
    logical     :: debug = .FALSE.

    ll = l
    if(ll .gt. 1) ltype = fact%factors(ll-1)  ! Set the next type of space curve

    !--------------------------------------------------------------
    !  Position [0,0]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,21) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'Cinco: After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [1,0]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,22) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [1,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [2,0]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,23) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [2,1]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,24) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [2,2]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = md
    lja       = ma
    ljd       = -md

    if(ll .gt. 1) then
       if(debug) write(iulog,25) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif


    !--------------------------------------------------------------
    !  Position [1,2]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = -md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,26) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [1,1]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = -md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,27) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif
    !--------------------------------------------------------------
    !  Position [0,1]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = -md
    lja       = MOD(ma+1,maxdim)
    ljd       = md

    if(ll .gt. 1) then
       if(debug) write(iulog,28) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [0,2]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,29) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [0,3]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,30) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [0,4]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,31) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [1,4]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = md
    lja       = MOD(ma+1,maxdim)
    ljd       = -md

    if(ll .gt. 1) then
       if(debug) write(iulog,32) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [1,3]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = -md
    lja       = ma
    ljd       = md

    if(ll .gt. 1) then
       if(debug) write(iulog,33) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [2,3]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,34) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [2,4]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,35) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [3,4]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,36) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [4,4]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = md
    lja       = MOD(ma+1,maxdim)
    ljd       = -md

    if(ll .gt. 1) then
       if(debug) write(iulog,37) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [4,3]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = -md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,38) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [3,3]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = -md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,39) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [3,2]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = -md
    lja       = ma
    ljd       = md

    if(ll .gt. 1) then
       if(debug) write(iulog,40) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [4,2]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = md
    lja       = MOD(ma+1,maxdim)
    ljd       = -md

    if(ll .gt. 1) then
       if(debug) write(iulog,41) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [4,1]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = -md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,42) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [3,1]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = -md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,43) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [3,0]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = -md
    lja       = ma
    ljd       = md

    if(ll .gt. 1) then
       if(debug) write(iulog,44) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

    !--------------------------------------------------------------
    !  Position [4,0]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = md
    lja       = ja
    ljd       = jd

    if(ll .gt. 1) then
       if(debug) write(iulog,45) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif

21  format('Call Cinco Pos [0,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
22  format('Call Cinco Pos [1,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
23  format('Call Cinco Pos [2,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
24  format('Call Cinco Pos [2,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
25  format('Call Cinco Pos [2,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
26  format('Call Cinco Pos [1,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
27  format('Call Cinco Pos [1,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
28  format('Call Cinco Pos [0,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
29  format('Call Cinco Pos [0,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
30  format('Call Cinco Pos [0,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
31  format('Call Cinco Pos [0,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
32  format('Call Cinco Pos [1,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
33  format('Call Cinco Pos [1,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
34  format('Call Cinco Pos [2,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
35  format('Call Cinco Pos [2,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
36  format('Call Cinco Pos [3,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
37  format('Call Cinco Pos [4,4] Level ',i1,' at (',i2,',',i2,')',4(i3))
38  format('Call Cinco Pos [4,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
39  format('Call Cinco Pos [3,3] Level ',i1,' at (',i2,',',i2,')',4(i3))
40  format('Call Cinco Pos [3,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
41  format('Call Cinco Pos [4,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
42  format('Call Cinco Pos [4,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
43  format('Call Cinco Pos [3,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
44  format('Call Cinco Pos [3,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
45  format('Call Cinco Pos [4,0] Level ',i1,' at (',i2,',',i2,')',4(i3))

  end function Cinco

  !---------------------------------------------------------
  recursive function PeanoM(l,type,ma,md,ja,jd) result(ierr)

    implicit none
    integer,intent(in) ::   l,type,ma,md,ja,jd

    integer ::  lma,lmd,lja,ljd,ltype
    integer :: ll
    integer :: ierr
    logical     :: debug = .FALSE.

    ll = l
    if(ll .gt. 1) ltype = fact%factors(ll-1)  ! Set the next type of space curve
    !--------------------------------------------------------------
    !  Position [0,0]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,21) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif


    !--------------------------------------------------------------
    ! Position [0,1]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = md
    lja       = lma
    ljd       = lmd
    if(ll .gt. 1) then
       if(debug) write(iulog,22) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,1] ',pos
    endif

    !--------------------------------------------------------------
    ! Position [0,2]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = md
    lja       = lma
    ljd       = lmd
    if(ll .gt. 1) then
       if(debug) write(iulog,23) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,2] ',pos
    endif

    !--------------------------------------------------------------
    ! Position [1,2]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = md
    lja       = lma
    ljd       = lmd
    if(ll .gt. 1) then
       if(debug) write(iulog,24) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [1,2] ',pos
    endif


    !--------------------------------------------------------------
    ! Position [2,2]
    !--------------------------------------------------------------
    lma        = ma
    lmd        = md
    lja        = MOD(lma+1,maxdim)
    ljd        = -lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,25) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [2,2] ',pos
    endif

    !--------------------------------------------------------------
    ! Position [2,1]
    !--------------------------------------------------------------
    lma        = ma
    lmd        = -md
    lja        = lma
    ljd        = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,26) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [2,1] ',pos
    endif

    !--------------------------------------------------------------
    ! Position [1,1]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = -md
    lja        = lma
    ljd        = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,27) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [1,1] ',pos
    endif


    !--------------------------------------------------------------
    ! Position [1,0]
    !--------------------------------------------------------------
    lma        = MOD(ma+1,maxdim)
    lmd        = -md
    lja        = MOD(lma+1,maxdim)
    ljd        = -lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,28) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [1,0] ',pos
    endif

    !--------------------------------------------------------------
    ! Position [2,0]
    !--------------------------------------------------------------
    lma        = ma
    lmd        = md
    lja        = ja
    ljd        = jd

    if(ll .gt. 1) then
       if(debug) write(iulog,29) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [2,0] ',pos
    endif

21  format('Call PeanoM Pos [0,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
22  format('Call PeanoM Pos [0,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
23  format('Call PeanoM Pos [0,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
24  format('Call PeanoM Pos [1,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
25  format('Call PeanoM Pos [2,2] Level ',i1,' at (',i2,',',i2,')',4(i3))
26  format('Call PeanoM Pos [2,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
27  format('Call PeanoM Pos [1,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
28  format('Call PeanoM Pos [1,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
29  format('Call PeanoM Pos [2,0] Level ',i1,' at (',i2,',',i2,')',4(i3))

  end function PeanoM
  !---------------------------------------------------------
  recursive function hilbert(l,type,ma,md,ja,jd) result(ierr)

    implicit none
    integer,intent(in) ::   l,type,ma,md,ja,jd

    integer ::  lma,lmd,lja,ljd,ltype
    integer :: ll
    integer :: ierr
    logical     :: debug = .FALSE.

    ll = l
    if(ll .gt. 1) ltype = fact%factors(ll-1)  ! Set the next type of space curve
    !--------------------------------------------------------------
    !  Position [0,0]
    !--------------------------------------------------------------
    lma       = MOD(ma+1,maxdim)
    lmd       = md
    lja       = lma
    ljd       = lmd

    if(ll .gt. 1) then
       if(debug) write(iulog,21) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,0] ',pos
    endif


    !--------------------------------------------------------------
    ! Position [0,1]
    !--------------------------------------------------------------
    lma       = ma
    lmd       = md
    lja       = lma
    ljd       = lmd
    if(ll .gt. 1) then
       if(debug) write(iulog,22) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [0,1] ',pos
    endif


    !--------------------------------------------------------------
    ! Position [1,1]
    !--------------------------------------------------------------
    lma        = ma
    lmd        = md
    lja        = MOD(ma+1,maxdim)
    ljd        = -md

    if(ll .gt. 1) then
       if(debug) write(iulog,23) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [1,1] ',pos
    endif

    !--------------------------------------------------------------
    ! Position [1,0]
    !--------------------------------------------------------------
    lma        = MOD(ma+1,maxdim)
    lmd        = -md
    lja        = ja
    ljd        = jd

    if(ll .gt. 1) then
       if(debug) write(iulog,24) ll-1,pos(0),pos(1),lma,lmd,lja,ljd
       ierr  = GenCurve(ll-1,ltype,lma,lmd,lja,ljd)
       if(debug) call PrintCurve(ordered)
    else
       ierr = IncrementCurve(lja,ljd)
       if(debug) write(iulog,*)'After Position [1,0] ',pos
    endif

21  format('Call Hilbert Pos [0,0] Level ',i1,' at (',i2,',',i2,')',4(i3))
22  format('Call Hilbert Pos [0,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
23  format('Call Hilbert Pos [1,1] Level ',i1,' at (',i2,',',i2,')',4(i3))
24  format('Call Hilbert Pos [1,0] Level ',i1,' at (',i2,',',i2,')',4(i3))

  end function hilbert
  !---------------------------------------------------------
  function IncrementCurve(ja,jd) result(ierr)

    implicit none

    integer    :: ja,jd
    integer    :: ierr

    ordered(pos(0)+1,pos(1)+1) = vcnt
    vcnt  = vcnt + 1
    pos(ja) = pos(ja) + jd

    ierr = 0
  end function IncrementCurve
  !---------------------------------------------------------
  recursive function hilbert_old(l,d,ma,md,ja,jd) result(ierr)

    integer  :: l,d             ! log base 2 of levels and dimensions left
    integer  :: ma,md           ! main axis and direction 
    integer  :: ja,jd           ! joiner axis and direction 

    integer  :: ierr
    integer  :: axis
    integer  :: ll

    if(verbose) write(iulog,10) l,d,ma,md,ja,jd,pos(0),pos(1)
    ll = l     ! Copy this to a temporary variable 
    if(d == 0) then 
       ll=ll-1
       if(ll == 0) then 
          return 
       endif
       axis = ja
       if(dir(ll,axis) /= jd) then     ! do not move away from joiner plane 
          axis = MOD(axis+1,maxdim)    ! next axis
       endif
       if(verbose) write(iulog,*)'hilbert_old: call hilbert_old(l,d) #1:'
       ierr = hilbert_old(ll,maxdim,axis,dir(ll,axis),ja,jd) 
       dir(ll,ja) = -dir(ll,ja)
       return 
    endif
    axis = MOD(ma+1,maxdim)
    if(verbose) write(iulog,*)'hilbert_old: before call hilbert_old(l,d) #2:'
    ierr = hilbert_old(ll,d-1,axis,dir(ll,axis),ma,md)
    if(verbose) write(iulog,*)'hilbert_old:  after call hilbert_old(l,d) #2:'
    if(verbose) write(iulog,30) l,d,ma,md,ja,jd,pos(0),pos(1)


    pos(ma) = pos(ma) + md
    dir(ll,ma) = - dir(ll,ma)

    !----------------------------------
    !  Mark this node as visited 
    !----------------------------------
    if(verbose) write(iulog,20) l,d,ma,md,ja,jd,pos(0),pos(1)
    vcnt=vcnt+1
    if(verbose) write(iulog,15) pos(0)+1,pos(1)+1,vcnt 
    if(verbose) write(iulog,*)'  '
    if(verbose) write(iulog,*)'  '
    ordered(pos(0)+1,pos(1)+1)=vcnt

    if(verbose) write(iulog,*)'hilbert_old: before call hilbert_old(l,d) #3:'
    ierr = hilbert_old(ll,d-1,axis,dir(ll,axis),ja,jd)
    if(verbose) write(iulog,*)'hilbert_old:  after call hilbert_old(l,d) #3:'

10  format('hilbert_old: Entering hilbert_old (l,d,ma,md,ja,jd) are: ', &
         2(i4),'  [',2(i3),'][',2(i3),']',2(i3))
15  format('hilbert_old: mark element {x,y,ordered}:',3(i4))
20  format('hilbert_old: Before visit code (l,d,ma,md,ja,jd) are:', &
         2(i4),'  [',2(i3),'][',2(i3),']',2(i3))

30  format('hilbert_old:  after call hilbert_old(l,d) #2: (l,d,ma,md,ja,jd are:', &
         2(i4),'  [',2(i3),'][',2(i3),']',2(i3))

  end function hilbert_old
  !---------------------------------------------------------
  function log2( n)

    implicit none

    integer  :: n

    integer  :: log2,tmp
    ! 
    !  Find the log2 of input value
    !
    log2 = 1
    tmp =n
    do while (tmp/2 .ne. 1) 
       tmp=tmp/2
       log2=log2+1
    enddo

  end function log2
  !---------------------------------------------------------
  function  IsLoadBalanced(nelem,npart)

    implicit none 

    integer        :: nelem,npart

    logical        :: IsLoadBalanced

    integer        :: tmp1

    tmp1 = nelem/npart

    if(npart*tmp1 == nelem ) then 
       IsLoadBalanced=.TRUE.
    else
       IsLoadBalanced=.FALSE.
    endif

  end function IsLoadBalanced
  !---------------------------------------------------------
  recursive function GenCurve(l,type,ma,md,ja,jd) result(ierr)

    implicit none
    integer,intent(in)               :: l,type,ma,md,ja,jd
    integer                          :: ierr

    if(type == 2) then
       ierr = hilbert(l,type,ma,md,ja,jd)
    elseif ( type == 3) then
       ierr = PeanoM(l,type,ma,md,ja,jd)
    elseif ( type == 5) then
       ierr = Cinco(l,type,ma,md,ja,jd)
    endif

  end function GenCurve
  !---------------------------------------------------------
  function Factor(num) result(res)

    implicit none
    integer,intent(in)  :: num

    type (factor_t)     :: res
    integer             :: tmp,tmp2,tmp3,tmp5
    integer             :: i,n
    logical             :: found

    ! --------------------------------------
    ! Allocate for max # of factors
    ! --------------------------------------
    tmp = num
    tmp2 = log2(num)
    allocate(res%factors(tmp2))

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

    tmp=1
    do i=1,n
       tmp = tmp * res%factors(i)
    enddo
    if(tmp == num) then
       res%numfact = n
    else
       res%numfact = -1
    endif

  end function Factor
  !---------------------------------------------------------
  function IsFactorable(n)

    implicit none 

    integer,intent(in)  :: n
    type (factor_t)     :: fact

    logical  :: IsFactorable

    fact = Factor(n)
    if(fact%numfact .ne. -1) then
       IsFactorable = .TRUE.
    else
       IsFactorable = .FALSE.
    endif

  end function IsFactorable
  !------------------------------------------------
  subroutine map(l)

    implicit none
    integer :: l,d
    integer :: type, ierr

    d = SIZE(pos)

    pos=0
    maxdim=d
    vcnt=0

    type = fact%factors(l)
       ierr = GenCurve(l,type,0,1,0,1)

     end subroutine map
     !---------------------------------------------------------
     subroutine GenSpaceCurve(Mesh) 

       implicit none

       integer,target,intent(inout) :: Mesh(:,:)
       integer :: level,dim 

       integer :: gridsize

       !  Setup the size of the grid to traverse

       dim = 2
       gridsize = SIZE(Mesh,dim=1)
       fact     = factor(gridsize)
       level    = fact%numfact

       if(verbose) write(iulog,*)'GenSpacecurve: level is ',level
       allocate(ordered(gridsize,gridsize))

       ! Setup the working arrays for the traversal
       allocate(pos(0:dim-1))

       !  The array ordered will contain the visitation order
       ordered(:,:) = 0

       call map(level) 

       Mesh(:,:) = ordered(:,:)

     end subroutine GenSpaceCurve
     !------------------------------------------------------------------------------------------------------- 
     subroutine PrintCurve(Mesh)
       implicit none
       integer,target ::  Mesh(:,:)
       integer :: gridsize,i

       gridsize = SIZE(Mesh,dim=1)

       if(gridsize == 2) then
          write (iulog,*) "A Level 1 Hilbert Curve:"
          write (iulog,*) "------------------------"
          do i=1,gridsize
             write(iulog,2) Mesh(1,i),Mesh(2,i)
          enddo
       else if(gridsize == 3) then
          write (iulog,*) "A Level 1 Peano Meandering Curve:"
          write (iulog,*) "---------------------------------"
          do i=1,gridsize
             write(iulog,3) Mesh(1,i),Mesh(2,i),Mesh(3,i)
          enddo
       else if(gridsize == 4) then
          write (iulog,*) "A Level 2 Hilbert Curve:"
          write (iulog,*) "------------------------"
          do i=1,gridsize
             write(iulog,4) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i)
          enddo
       else if(gridsize == 5) then
          write (iulog,*) "A Level 1 Cinco Curve:"
          write (iulog,*) "------------------------"
          do i=1,gridsize
             write(iulog,5) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i),Mesh(5,i)
          enddo
       else if(gridsize == 6) then
          write (iulog,*) "A Level 1 Hilbert and Level 1 Peano Curve:"
          write (iulog,*) "------------------------------------------"
          do i=1,gridsize
             write(iulog,6) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i),Mesh(5,i),Mesh(6,i)
          enddo
       else if(gridsize == 8) then
          write (iulog,*) "A Level 3 Hilbert Curve:"
          write (iulog,*) "------------------------"
          do i=1,gridsize
             write(iulog,8) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                  Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i)
          enddo
       else if(gridsize == 9) then
          write (iulog,*) "A Level 2 Peano Meandering Curve:"
          write (iulog,*) "---------------------------------"
          do i=1,gridsize
             write(iulog,9) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                  Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i), &
                  Mesh(9,i)
          enddo
       else if(gridsize == 10) then
          write (iulog,*) "A Level 1 Hilbert and Level 1 Cinco Curve:"
          write (iulog,*) "---------------------------------"
          do i=1,gridsize
             write(iulog,10) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                  Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i), &
                  Mesh(9,i),Mesh(10,i)
          enddo
       else if(gridsize == 12) then
          write (iulog,*) "A Level 2 Hilbert and Level 1 Peano Curve:"
          write (iulog,*) "------------------------------------------"
          do i=1,gridsize
             write(iulog,12) Mesh(1,i),Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                  Mesh(5,i),Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                  Mesh(9,i),Mesh(10,i),Mesh(11,i),Mesh(12,i)
          enddo
       else if(gridsize == 15) then
          write (iulog,*) "A Level 1 Peano and Level 1 Cinco Curve:"
          write (iulog,*) "------------------------"
          do i=1,gridsize
             write(iulog,15) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                  Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i), &
                  Mesh(9,i),Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                  Mesh(13,i),Mesh(14,i),Mesh(15,i)
          enddo
       else if(gridsize == 16) then
          write (iulog,*) "A Level 4 Hilbert Curve:"
          write (iulog,*) "------------------------"
          do i=1,gridsize
             write(iulog,16) Mesh(1,i),Mesh(2,i),Mesh(3,i),Mesh(4,i), &
                  Mesh(5,i),Mesh(6,i),Mesh(7,i),Mesh(8,i), &
                  Mesh(9,i),Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                  Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i)
          enddo
       else if(gridsize == 18) then
          write (iulog,*) "A Level 1 Hilbert and Level 2 Peano Curve:"
          write (iulog,*) "------------------------------------------"
          do i=1,gridsize
             write(iulog,18) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                  Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                  Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                  Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                  Mesh(17,i),Mesh(18,i)
          enddo
       else if(gridsize == 20) then
          write (iulog,*) "A Level 2 Hilbert and Level 1 Cinco Curve:"
          write (iulog,*) "------------------------------------------"
          do i=1,gridsize
             write(iulog,20) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                  Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                  Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                  Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                  Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i)
          enddo
       else if(gridsize == 24) then
          write (iulog,*) "A Level 3 Hilbert and Level 1 Peano Curve:"
          write (iulog,*) "------------------------------------------"
          do i=1,gridsize
             write(iulog,24) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                  Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                  Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                  Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                  Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i), &
                  Mesh(21,i),Mesh(22,i),Mesh(23,i),Mesh(24,i)
          enddo
       else if(gridsize == 25) then
          write (iulog,*) "A Level 2 Cinco Curve:"
          write (iulog,*) "------------------------------------------"
          do i=1,gridsize
             write(iulog,25) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                  Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                  Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                  Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                  Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i), &
                  Mesh(21,i),Mesh(22,i),Mesh(23,i),Mesh(24,i), &
                  Mesh(25,i)
          enddo
       else if(gridsize == 27) then
          write (iulog,*) "A Level 3 Peano Meandering Curve:"
          write (iulog,*) "---------------------------------"
          do i=1,gridsize
             write(iulog,27) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i), &
                  Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i), &
                  Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                  Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                  Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i), &
                  Mesh(21,i),Mesh(22,i),Mesh(23,i),Mesh(24,i), &
                  Mesh(25,i),Mesh(26,i),Mesh(27,i)
          enddo
       else if(gridsize == 32) then
          write (iulog,*) "A Level 5 Hilbert Curve:"
          write (iulog,*) "------------------------"
          do i=1,gridsize
             write(iulog,32) Mesh(1,i), Mesh(2,i), Mesh(3,i), Mesh(4,i),  &
                  Mesh(5,i), Mesh(6,i), Mesh(7,i), Mesh(8,i),  &
                  Mesh(9,i), Mesh(10,i),Mesh(11,i),Mesh(12,i), &
                  Mesh(13,i),Mesh(14,i),Mesh(15,i),Mesh(16,i), &
                  Mesh(17,i),Mesh(18,i),Mesh(19,i),Mesh(20,i), &
                  Mesh(21,i),Mesh(22,i),Mesh(23,i),Mesh(24,i), &
                  Mesh(25,i),Mesh(26,i),Mesh(27,i),Mesh(28,i), &
                  Mesh(29,i),Mesh(30,i),Mesh(31,i),Mesh(32,i)
          enddo
       endif
2      format('|',2(i2,'|'))
3      format('|',3(i2,'|'))
4      format('|',4(i2,'|'))
5      format('|',5(i2,'|'))
6      format('|',6(i2,'|'))
8      format('|',8(i2,'|'))
9      format('|',9(i2,'|'))
10     format('|',10(i2,'|'))
12     format('|',12(i3,'|'))
15     format('|',15(i3,'|'))
16     format('|',16(i3,'|'))
18     format('|',18(i3,'|'))
20     format('|',20(i3,'|'))
24     format('|',24(i3,'|'))
25     format('|',25(i3,'|'))
27     format('|',27(i3,'|'))
32     format('|',32(i4,'|'))

     end subroutine PrintCurve

     !------------------------------------------------------------------------------------------------------- 
     subroutine genspacepart(GridEdge,GridVertex)
       use dimensions_mod, only : npart
       use gridgraph_mod, only : gridedge_t, gridvertex_t


       implicit none 

       type (GridVertex_t), intent(inout) :: GridVertex(:)
       type (GridEdge_t),   intent(inout) :: GridEdge(:)

       integer               :: nelem,nelem_edge,nelemd
       integer               :: head_part,tail_part
       integer               :: j,k,tmp1,id,s1,extra

       nelem      = SIZE(GridVertex(:))
       nelem_edge = SIZE(GridEdge(:))

       nelemd = nelem/npart
       ! every cpu gets nelemd elements, but the first 'extra' get nelemd+1
       extra = mod(nelem,npart)
       s1 = extra*(nelemd+1)

       ! split curve into two curves:
       ! 1 ... s1  s2 ... nelem
       !
       !  s1 = extra*(nelemd+1)         (count be 0)
       !  s2 = s1+1 
       !
       ! First region gets nelemd+1 elements per Processor
       ! Second region gets nelemd elements per Processor

       ! ===========================================
       !  Add the partitioning information into the
       !    Grid Vertex and Grid Edge structures
       ! ===========================================

       do k=1,nelem
          id=GridVertex(k)%SpaceCurve
          if (id<=s1) then
             tmp1 = id/(nelemd+1)
             GridVertex(k)%processor_number = tmp1+1
          else
             id=id-s1
             tmp1 = id/nelemd
             GridVertex(k)%processor_number = extra + tmp1+1
          endif
       enddo
#if 0
       write(iulog,*)'Space-Filling Curve Parititioning: '
       do k=1,nelem
          write(iulog,*) k,GridVertex(k)%processor_number
       enddo
       stop 'halting: at the end of genspacepart'
#endif

     end subroutine genspacepart

   end module spacecurve_mod
