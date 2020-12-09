#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module spacecurve_mod
!
!
! Revisions:
! John Dennis, NCAR: Initial
! Mark Taylor: 2018/3  Remove memory leaks, make most routines private
! Mark Taylor: 2018/10 add more deallocates
! AMB: 2018/10  Add sfcmap_* (i,j) <-> SFC index routines
!
  use kinds, only : iulog
  implicit none
  private

  type, public :: factor_t
     integer                       :: numfact
     integer, dimension(:),pointer :: factors
  end type factor_t


  integer, dimension(:,:), allocatable :: ordered
  integer, dimension(:)  , allocatable :: pos     ! position along each of the axes

  integer,public                              :: maxdim  ! dimensionality of entire space
  integer,public                              :: vcnt   ! visitation count
  logical,private                             :: verbose=.FALSE.

  type (factor_t), private                  :: fact

  SAVE:: fact
  public :: GenSpaceCurve
  public :: PrintCurve
  public :: IsFactorable,IsLoadBalanced
  public :: genspacepart
  public :: GilbertCurve

  ! Map (i,j) <-> SFC index in O(log ne) time. Unlike the above routines,
  ! nothing like a mesh(ne,ne) is allocated; these routines use O(log ne) memory
  ! rather than O(ne^2).
  type, public :: sfcmap_t
     type (factor_t) :: fact
     integer :: n, index_os, index, pos_os(2), pos(2)
     logical :: pos2i
     ! For use when working on a range of SFC indices:
     integer :: idxs, idxe
     integer, pointer :: positions(:,:)
  end type sfcmap_t
  ! Call sfcmap_init before calling the map routines sfcmap_i2pos/pos2i and
  ! sfcmap_indexrange2pos. Call sfcmap_finalize at the end. sfcmap_test is a
  ! non-MPI correctness check, so it makes sense to call it only if
  ! par%masterproc.
  !   In input and outuput, SFC index space starts at 0, and position starts at
  ! (1,1). (Internally, pos starts at (0,0), but this is not apparent to the
  ! caller.)
  public :: &
       ! Set up and tear down the internal data structure.
       sfcmap_init, sfcmap_finalize, &
       ! Map one SFC index to (i,j), and one (i,j) to SFC index.
       sfcmap_i2pos, sfcmap_pos2i, &
       ! Map a range of SFC indices, [idxs, idxe], to a list of positions
       ! pos(idxe-idxs+1,2). This is faster than calling sfcmap_i2pos on each
       ! index individually.
       sfcmap_indexrange2pos, &
       ! Run unit tests. This is run alone, i.e., w/o initialization or
       ! finalization.
       sfcmap_test

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
  ! note: this function is allocating memory in the 'fact' struct
  ! to avoid memory leaks, the calling program needs to deallocate
  ! this array.  poor design choice.

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
    deallocate(fact%factors)

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

       deallocate(fact%factors)
       deallocate(ordered)
       deallocate(pos)

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

  !-----------------------------------------------------------------------------
  ! O(log ne) (i,j) <-> SFC index maps.
  !
  ! index space starts at 0
  ! pos internally starts at (0,0); externally it starts at (1,1)

  function sfcmap_init(ne, sfcmap) result(ierr)
    integer, intent(in) :: ne
    type (sfcmap_t), intent(out) :: sfcmap
    integer :: i, ierr, f

    sfcmap%n = ne
    sfcmap%fact = factor(ne)
    ierr = 0
    if (sfcmap%fact%numfact == -1) ierr = -1
  end function sfcmap_init

  subroutine sfcmap_finalize(sfcmap)
    type (sfcmap_t), intent(inout) :: sfcmap
    deallocate(sfcmap%fact%factors)
  end subroutine sfcmap_finalize

  function sfcmap_pos2i(s, pos, index) result(status)
    type (sfcmap_t), intent(inout) :: s
    integer, intent(in) :: pos(2)
    integer, intent(out) :: index
    integer :: status

    s%pos2i = .true.
    s%idxs = -1
    s%pos_os = (/ 0, 0 /)
    s%pos = pos - 1
    s%index = 0
    status = sfcmap_impl(s, s%pos_os, s%fact%numfact, s%n, 0, 1, 0, 1)
    index = s%index
  end function sfcmap_pos2i

  function sfcmap_i2pos(s, index, pos) result(status)
    type (sfcmap_t), intent(inout) :: s
    integer, intent(in) :: index
    integer, intent(out) :: pos(2)
    integer :: status

    s%pos2i = .false.
    s%idxs = -1
    s%index_os = 0
    s%index = index
    s%pos = (/ 0, 0 /)
    status = sfcmap_impl(s, s%pos, s%fact%numfact, s%n, 0, 1, 0, 1)
    pos = s%pos + 1
  end function sfcmap_i2pos

  function sfcmap_indexrange2pos(s, idxs, idxe, pos) result(status)
    type (sfcmap_t), intent(inout) :: s
    integer, intent(in) :: idxs, idxe
    integer, target, intent(inout) :: pos(:,:)
    integer :: index, status

    s%pos2i = .false.
    s%idxs = idxs
    s%idxe = idxe
    s%index_os = 0
    s%index = idxs
    s%pos = (/ 0, 0 /)
    s%positions => pos
    status = sfcmap_impl(s, s%pos, s%fact%numfact, s%n, 0, 1, 0, 1)
  end function sfcmap_indexrange2pos

  ! Top level implementation routine for the caller-level functions above.
  recursive function sfcmap_impl(s, pos, k, k_n, ma, md, ja, jd) result(o)
    type (sfcmap_t), intent(inout) :: s
    integer, intent(in) :: k, k_n, ma, md, ja, jd
    integer, intent(inout) :: pos(2)
    integer :: km1, km1_n, km1_n2, region, next_region, o

    if (k == 0) then
       ! In this level, each region is a point, so we're done.
       o = 0
       return
    end if

    km1 = k - 1                      ! child level
    km1_n = k_n / s%fact%factors(k)  ! size of child level
    km1_n2 = km1_n*km1_n

    ! Calculate the region (of the factor^2 possible) in which the input
    ! position or index lies.
    if (s%pos2i) then
       region = sfcmap_pos2region(s%fact%factors(k), km1_n, ma, s%pos, s%pos_os)
       s%index = s%index + km1_n2*region
    else
       region = (s%index - s%index_os) / km1_n2
       s%index_os = s%index_os + km1_n2*region
    end if

    ! Find either the only index/pos requested, or find the first in the range.
    o = sfcmap_impl_find_first(s, pos, k, km1, km1_n, region, ma, md, ja, jd)

    ! For indexrange2pos, now that we found the first point, follow the curve to
    ! the last point. Doing this for a range is faster than finding one position
    ! at a time.
    if (o == 0 .and. s%idxs >= 0) then
       if (k > 1) then
          next_region = region + 1
       else
          next_region = region ! repeat the first index b/c of new base case
       end if
       o = sfcmap_impl_follow(s, k, k_n, next_region, ma, md, ja, jd)
    end if
  end function sfcmap_impl

  recursive function sfcmap_impl_find_first(s, pos, k, km1, km1_n, region, ma,  md, ja, jd) &
       result(o)
    type (sfcmap_t), intent(inout) :: s
    integer, intent(in) :: k, km1, km1_n, region, ma, md, ja, jd
    integer, intent(inout) :: pos(2)
    integer :: ima, km1_n_md, o

    ! The pos increment line moves the cursor to the starting point of the child
    ! region containing the input position or index. Thus, if you're trying to
    ! understand the code, draw a 2-level curve so that the starting point is
    ! relevant. (In a 1-level curve, the child region is just a point, and thus
    ! the starting point is not revealed.) Then the function recurses on that
    ! region.
    ima = modulo(ma+1, 2) ! other indexer; one is for x, the other for y
    km1_n_md = km1_n*md
    select case (s%fact%factors(k))
    case (2)
       !  _
       ! | |
       ! 0 3
       select case(region)
       case (0)
          o = sfcmap_impl(s, pos, km1, km1_n, ima,  md, ima,  md)
       case (1)
          pos(ima+1) = pos(ima+1) +   km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md,  ma,  md)
       case (2)
          pos( ma+1) = pos( ma+1) +   km1_n_md
          pos(ima+1) = pos(ima+1) +   km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md, ima, -md)
       case (3)
          pos( ma+1) = pos( ma+1) + 2*km1_n_md - md
          pos(ima+1) = pos(ima+1) +   km1_n_md - md
          o = sfcmap_impl(s, pos, km1, km1_n, ima, -md,  ja,  jd)
       end select
    case (3)
       !  _ _
       ! |  _|
       ! | |_
       ! 0   8
       select case (region)
       case (0)
          o = sfcmap_impl(s, pos, km1, km1_n, ima,  md, ima,  md)
       case (1)
          pos(ima+1) = pos(ima+1) +   km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n, ima,  md, ima,  md)
       case (2)
          pos(ima+1) = pos(ima+1) + 2*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md,  ma,  md)
       case (3)
          pos( ma+1) = pos( ma+1) +   km1_n_md
          pos(ima+1) = pos(ima+1) + 2*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md,  ma,  md)
       case (4)
          pos( ma+1) = pos( ma+1) + 2*km1_n_md
          pos(ima+1) = pos(ima+1) + 2*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md, ima, -md)
       case (5)
          pos( ma+1) = pos( ma+1) + 3*km1_n_md - md
          pos(ima+1) = pos(ima+1) + 2*km1_n_md - md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma, -md,  ma, -md)
       case (6)
          pos( ma+1) = pos( ma+1) + 2*km1_n_md - md
          pos(ima+1) = pos(ima+1) + 2*km1_n_md - md
          o = sfcmap_impl(s, pos, km1, km1_n, ima, -md, ima, -md)
       case (7)
          pos( ma+1) = pos( ma+1) + 2*km1_n_md - md
          pos(ima+1) = pos(ima+1) +   km1_n_md - md
          o = sfcmap_impl(s, pos, km1, km1_n, ima, -md,  ma,  md)
       case (8)
          pos( ma+1) = pos( ma+1) + 2*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md,  ja,  jd)
       end select
    case (5)
       !  _   _ _
       ! | |_|  _|
       ! |  _  |_
       ! |_| |  _|
       !  _ _| |_
       ! 0       24
       select case (region)
       case (0)
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md,  ma,  md)
       case (1)
          pos( ma+1) = pos( ma+1) +   km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md,  ma,  md)
       case (2)
          pos( ma+1) = pos( ma+1) + 2*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n, ima,  md, ima,  md)
       case (3)
          pos( ma+1) = pos( ma+1) + 2*km1_n_md
          pos(ima+1) = pos(ima+1) +   km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n, ima,  md, ima,  md)
       case (4)
          pos( ma+1) = pos( ma+1) + 2*km1_n_md
          pos(ima+1) = pos(ima+1) + 2*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n, ima,  md,  ma, -md)
       case (5)
          pos( ma+1) = pos( ma+1) + 2*km1_n_md - md
          pos(ima+1) = pos(ima+1) + 3*km1_n_md - md
          o = sfcmap_impl(s, pos, km1, km1_n, ima, -md, ima, -md)
       case (6)
          pos( ma+1) = pos( ma+1) + 2*km1_n_md - md
          pos(ima+1) = pos(ima+1) + 2*km1_n_md - md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma, -md,  ma, -md)
       case (7)
          pos( ma+1) = pos( ma+1) +   km1_n_md - md
          pos(ima+1) = pos(ima+1) + 2*km1_n_md - md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma, -md, ima,  md)
       case (8)
          pos( ma+1) = pos( ma+1)
          pos(ima+1) = pos(ima+1) + 2*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n, ima,  md, ima,  md)
       case (9)
          pos( ma+1) = pos( ma+1)
          pos(ima+1) = pos(ima+1) + 3*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n, ima,  md, ima,  md)
       case (10)
          pos( ma+1) = pos( ma+1)
          pos(ima+1) = pos(ima+1) + 4*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md,  ma,  md)
       case (11)
          pos( ma+1) = pos( ma+1) +   km1_n_md
          pos(ima+1) = pos(ima+1) + 4*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md, ima, -md)
       case (12)
          pos( ma+1) = pos( ma+1) + 2*km1_n_md - md
          pos(ima+1) = pos(ima+1) + 4*km1_n_md - md
          o = sfcmap_impl(s, pos, km1, km1_n, ima, -md,  ma,  md)
       case (13)
          pos( ma+1) = pos( ma+1) + 2*km1_n_md
          pos(ima+1) = pos(ima+1) + 3*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n, ima,  md, ima,  md)
       case (14)
          pos( ma+1) = pos( ma+1) + 2*km1_n_md
          pos(ima+1) = pos(ima+1) + 4*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md,  ma,  md)
       case (15)
          pos( ma+1) = pos( ma+1) + 3*km1_n_md
          pos(ima+1) = pos(ima+1) + 4*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md,  ma,  md)
       case (16)
          pos( ma+1) = pos( ma+1) + 4*km1_n_md
          pos(ima+1) = pos(ima+1) + 4*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md, ima, -md)
       case (17)
          pos( ma+1) = pos( ma+1) + 5*km1_n_md - md
          pos(ima+1) = pos(ima+1) + 4*km1_n_md - md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma, -md,  ma, -md)
       case (18)
          pos( ma+1) = pos( ma+1) + 4*km1_n_md - md
          pos(ima+1) = pos(ima+1) + 4*km1_n_md - md
          o = sfcmap_impl(s, pos, km1, km1_n, ima, -md, ima, -md)
       case (19)
          pos( ma+1) = pos( ma+1) + 4*km1_n_md - md
          pos(ima+1) = pos(ima+1) + 3*km1_n_md - md
          o = sfcmap_impl(s, pos, km1, km1_n, ima, -md,  ma,  md)
       case (20)
          pos( ma+1) = pos( ma+1) + 4*km1_n_md
          pos(ima+1) = pos(ima+1) + 2*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md, ima, -md)
       case (21)
          pos( ma+1) = pos( ma+1) + 5*km1_n_md - md
          pos(ima+1) = pos(ima+1) + 2*km1_n_md - md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma, -md,  ma, -md)
       case (22)
          pos( ma+1) = pos( ma+1) + 4*km1_n_md - md
          pos(ima+1) = pos(ima+1) + 2*km1_n_md - md
          o = sfcmap_impl(s, pos, km1, km1_n, ima, -md, ima, -md)
       case (23)
          pos( ma+1) = pos( ma+1) + 4*km1_n_md - md
          pos(ima+1) = pos(ima+1) +   km1_n_md - md
          o = sfcmap_impl(s, pos, km1, km1_n, ima, -md,  ma,  md)
       case (24)
          pos( ma+1) = pos( ma+1) + 4*km1_n_md
          o = sfcmap_impl(s, pos, km1, km1_n,  ma,  md,  ja,  jd)
       end select
    case default
       print *, 'error: factors =', s%fact%factors
       o = -1
    end select
  end function sfcmap_impl_find_first

  function sfcmap_pos2region(factor, km1_n, ma, pos, pos_os) result(reg)
    integer, intent(in) :: factor, km1_n, ma, pos(2), pos_os(2)
    integer :: reg, x, y, tmp

    x = abs(pos(1) - pos_os(1)) / km1_n
    y = abs(pos(2) - pos_os(2)) / km1_n
    if (ma == 1) then
       tmp = x
       x = y
       y = tmp
    end if
    select case (factor)
    case (2)
       if (x == 0) then
          if (y == 0) then
             reg = 0
          else
             reg = 1
          end if
       else
          if (y == 0) then
             reg = 3
          else
             reg = 2
          end if
       end if
    case (3)
       select case (x)
       case (0)
          reg = y
       case (1)
          select case (y)
          case (0); reg =  7
          case (1); reg =  6
          case (2); reg =  3
          end select
       case (2)
          select case (y)
          case (0); reg =  8
          case (1); reg =  5
          case (2); reg =  4
          end select
       end select
    case (5)
       select case (x)
       case (0)
          select case (y)
          case (0); reg =  0
          case (1); reg =  7
          case (2); reg =  8
          case (3); reg =  9
          case (4); reg = 10
          end select
       case (1)
          select case (y)
          case (0); reg =  1
          case (1); reg =  6
          case (2); reg =  5
          case (3); reg = 12
          case (4); reg = 11
          end select
       case (2)
          select case (y)
          case (0); reg =  2
          case (1); reg =  3
          case (2); reg =  4
          case (3); reg = 13
          case (4); reg = 14
          end select
       case (3)
          select case (y)
          case (0); reg = 23
          case (1); reg = 22
          case (2); reg = 19
          case (3); reg = 18
          case (4); reg = 15
          end select
       case (4)
          select case (y)
          case (0); reg = 24
          case (1); reg = 21
          case (2); reg = 20
          case (3); reg = 17
          case (4); reg = 16
          end select
       end select
    end select
  end function sfcmap_pos2region

  recursive function sfcmap_impl_follow(s, k, k_n, region_min, &
       ma, md, ja, jd) result(o)
    type (sfcmap_t), intent(inout) :: s
    integer, intent(in) :: k, k_n, region_min, ma, md, ja, jd
    integer :: km1, km1_n, km1_n2, ima, region, o

    if (k == 0) then ! base case
       s%index = s%index + 1
       s%positions(s%index-s%idxs,1) = s%pos(1) + 1
       s%positions(s%index-s%idxs,2) = s%pos(2) + 1
       s%pos(ja+1) = s%pos(ja+1) + jd
       o = 0
       if (s%index == s%idxe+1) o = 1
       return
    end if

    km1 = k - 1
    km1_n = k_n / s%fact%factors(k)
    ima = modulo(ma+1, 2)

    select case (s%fact%factors(k))
    case (2)
       do region = region_min, 3
          select case(region)
          case (0); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima,  md, ima,  md)
          case (1); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md,  ma,  md)
          case (2); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md, ima, -md)
          case (3); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima, -md,  ja,  jd)
          end select
          if (o == 1) return
       end do
    case (3)
       do region = region_min, 8
          select case (region)
          case (0); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima,  md, ima,  md)
          case (1); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima,  md, ima,  md)
          case (2); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md,  ma,  md)
          case (3); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md,  ma,  md)
          case (4); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md, ima, -md)
          case (5); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma, -md,  ma, -md)
          case (6); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima, -md, ima, -md)
          case (7); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima, -md,  ma,  md)
          case (8); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md,  ja,  jd)
          end select
          if (o == 1) return
       end do
    case (5)
       do region = region_min, 24
          select case (region)
          case ( 0); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md,  ma,  md)
          case ( 1); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md,  ma,  md)
          case ( 2); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima,  md, ima,  md)
          case ( 3); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima,  md, ima,  md)
          case ( 4); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima,  md,  ma, -md)
          case ( 5); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima, -md, ima, -md)
          case ( 6); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma, -md,  ma, -md)
          case ( 7); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma, -md, ima,  md)
          case ( 8); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima,  md, ima,  md)
          case ( 9); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima,  md, ima,  md)
          case (10); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md,  ma,  md)
          case (11); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md, ima, -md)
          case (12); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima, -md,  ma,  md)
          case (13); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima,  md, ima,  md)
          case (14); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md,  ma,  md)
          case (15); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md,  ma,  md)
          case (16); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md, ima, -md)
          case (17); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma, -md,  ma, -md)
          case (18); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima, -md, ima, -md)
          case (19); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima, -md,  ma,  md)
          case (20); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md, ima, -md)
          case (21); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma, -md,  ma, -md)
          case (22); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima, -md, ima, -md)
          case (23); o = sfcmap_impl_follow(s, km1, km1_n, 0, ima, -md,  ma,  md)
          case (24); o = sfcmap_impl_follow(s, km1, km1_n, 0,  ma,  md,  ja,  jd)
          end select
          if (o == 1) return
       end do
    end select
    o = 0
  end function sfcmap_impl_follow

  subroutine sfcmap_test(verbose)
    logical, intent(in) :: verbose
    integer, allocatable :: mesh(:,:)
    integer :: ne, ine, index, pos(2), ierr, i, j, len, idxs, idxe
    real :: start, finish
    type (sfcmap_t) :: sfcmap
    integer, allocatable :: positions(:,:)

    integer, parameter :: nes(8) = (/ 11, 8, 27, 75, 30, 36, 750, 1080 /)

    do ine = 1, size(nes)
       ne = nes(ine)
       ierr = sfcmap_init(ne, sfcmap)
       call cpu_time(start)
       if (ne /= 11 .and. ierr /= 0) then
          print *, 'SGI> not impled for ne',ne,'with factors',sfcmap%fact%factors
       end if
       if (ierr /= 0) then
          call sfcmap_finalize(sfcmap)
          cycle
       end if
       allocate(mesh(ne,ne))
       call cpu_time(start)
       call genspacecurve(mesh)
       call cpu_time(finish)
       if (verbose) print *, 'SGI>  gsc  et', ne, finish-start
       do index = 0, ne**2 - 1
          ierr = sfcmap_i2pos(sfcmap, index, pos)
          if (ierr /= 0) then
             print *, 'not impled for ne',ne,'with factors',sfcmap%fact%factors
             exit
          end if
          if (mesh(pos(1),pos(2)) /= index) then
             print '(a,i4,i6,i6,i3,i3)', 'SGI>', ne, index, mesh(pos(1),pos(2)), pos(1), pos(2)
          end if
       end do
       call cpu_time(finish)
       if (verbose) print *, 'SGI>  p2i  et', ne, finish-start
       index = -1
       call cpu_time(start)
       do j = 1, ne
          do i = 1, ne
             pos = (/ i, j /)
             ierr = sfcmap_pos2i(sfcmap, pos, index)
             if (mesh(i,j) /= index) then
                print '(a,i4,i6,i6,i3,i3)', 'SGI>', ne, mesh(i,j), index, i, j
             end if
          end do
       end do
       call cpu_time(finish)
       if (verbose) print *, 'SGI>  i2p  et', ne, finish-start
       if (ne >= 128) then
          idxs = ne*(ne/5);
          idxe = min(ne**2 - 1, idxs + (ne*ne)/32)
       else
          idxs = 0
          idxe = ne**2 - 1
       end if
       allocate(positions(idxe-idxs+1, 2))
       call cpu_time(start)
       ierr = sfcmap_indexrange2pos(sfcmap, idxs, idxe, positions)
       call cpu_time(finish)
       if (verbose) print *, 'SGI> is2ps et', real(size(positions,1))/(ne**2), finish-start
       do i = 1, idxe-idxs+1
          ierr = sfcmap_pos2i(sfcmap, positions(i,:), index)
          if (mesh(positions(i,1), positions(i,2)) /= index) then
             print '(a,i4,i6,i6,i3,i3)', 'SGI>', ne, &
                  mesh(positions(i,1), positions(i,2)), &
                  index, positions(i,1), positions(i,2)
          end if
       end do
       positions = 0
       len = 7
       if (ne**2 > len) then
          call cpu_time(start)
          do j = 0, ne**2 - len, len
             ierr = sfcmap_indexrange2pos(sfcmap, j, j+len-1, positions)
             do i = 1,len
                pos = positions(i,:)
                ierr = sfcmap_pos2i(sfcmap, pos, index)
                if (mesh(pos(1), pos(2)) /= index) then
                   print '(a,i4,i4,i4,i6,i6,i3,i3)', 'SGI>', &
                        ne, j, i, mesh(pos(1), pos(2)), index, pos(1), pos(2)
                end if
             end do
          end do
          call cpu_time(finish)
          if (verbose) print *, 'SGI> is2ps*et', ne, finish-start
       end if
       deallocate(positions)
       deallocate(mesh)
       call sfcmap_finalize(sfcmap)
    end do
  end subroutine sfcmap_test


  ! Generalized Hilbert ('gilbert') space-filling curve for arbitrary-sized
  ! 2D rectangular grids.
  ! Follows algorithm from https://github.com/jakubcerveny/gilbert

  subroutine GilbertCurve(Mesh)

    use dimensions_mod, only : ne_x, ne_y

    implicit none

    integer,target,intent(inout) :: Mesh(:,:)
    integer :: global_index_ctr

    global_index_ctr = 0
    if (ne_x >= ne_y) then
        call gilbert(Mesh, 0, 0, ne_x, 0, 0, ne_y, global_index_ctr)
    else
        call gilbert(Mesh, 0, 0, 0, ne_y, ne_x, 0, global_index_ctr)
    end if

  end subroutine GilbertCurve

  recursive subroutine gilbert(Mesh, ix, iy, iax, iay, ibx, iby, global_index_ctr)

    use dimensions_mod, only : ne_x, ne_y

    implicit none

    integer, intent(in) :: ix,iy,iax,iay,ibx,iby
    integer, intent(inout) :: global_index_ctr
    integer,target,intent(inout) :: Mesh(:,:)

    integer :: x, y, ax, ay, bx, by
    integer :: dax, day, dbx, dby, width, height, width2, height2, ax2, bx2, ay2, by2
    integer :: ii

    x=ix; y=iy; ax=iax; ay=iay; bx=ibx; by=iby ! we want the values here, since we are recursing


    width = abs(ax + ay)
    height = abs(bx + by)

    ! ax/ay are major direction, bx/by are minor direction
    ! merge gives 0 if ax==0 and sign(1,ax) if ax /=0
    ! end result is 1 if ax>0, -1 if ax<0 and 0 if ax=0
    dax = merge(0,sign(1,ax),ax==0)
    day = merge(0,sign(1,ay),ay==0)
    dbx = merge(0,sign(1,bx),bx==0)
    dby = merge(0,sign(1,by),by==0)

    !print *,x, y, ax, ay, bx, by, dax,day,dbx,dby

    !trivial row fill
    if (height == 1) then
      do ii=0,width-1,1
        ! print *, x,y,global_index_ctr
          Mesh(x+1,y+1) = global_index_ctr ! Fortran arrays are indexed starting with 1
          global_index_ctr = global_index_ctr + 1
          x = x + dax
          y = y + day
      end do
      return
    end if

    ! trivial column fill
    if (width == 1) then
      do ii=0,height-1,1
        ! print *, x,y,global_index_ctr
          Mesh(x+1,y+1) = global_index_ctr ! Fortran arrays are indexed starting with 1
          global_index_ctr = global_index_ctr + 1
          x = x + dbx
          y = y + dby
      end do
      return
    end if

! This is required to get "floor" division ie always rounds DOWN to nearest int, even for negative numbers
    ax2 = floor(real(ax)/2.0D0)
    ay2 = floor(real(ay)/2.0D0)
    bx2 = floor(real(bx)/2.0D0)
    by2 = floor(real(by)/2.0D0)

    width2 = abs(ax2 + ay2)
    height2 = abs(bx2 + by2)

    if (2*width > 3*height) then
      ! if width2 is odd, shift since even steps are prefered
        if ((MOD(width2,2) /= 0) .and. (width > 2)) then
             ax2 = ax2 + dax
             ay2 = ay2 + day
        end if

        !long case: split in two parts only
        !print *,"call 1"
        call gilbert(Mesh, x, y, ax2, ay2, bx, by, global_index_ctr)
        !print *,"call 2"
        call gilbert(Mesh, x+ax2, y+ay2, ax-ax2, ay-ay2, bx, by, global_index_ctr)
    else
       ! if height2 is odd, shift since even steps are prefered
        if ((MOD(height2,2) /= 0) .and. (height > 2)) then
            bx2 = bx2 + dbx
            by2 = by2 + dby
        end if
        ! standard case: one step up, one long horizontal, one step down
        !print *,"call 3"
        call gilbert(Mesh, x, y, bx2, by2, ax2, ay2, global_index_ctr)
        !print *,"call 4"
        call gilbert(Mesh, x+bx2, y+by2, ax, ay, bx-bx2, by-by2, global_index_ctr)
        !print *,"call 5"
        call gilbert(Mesh, x+(ax-dax)+(bx2-dbx), y+(ay-day)+(by2-dby),-bx2, -by2, -(ax-ax2), -(ay-ay2), global_index_ctr)
    end if

  end subroutine gilbert

end module spacecurve_mod
