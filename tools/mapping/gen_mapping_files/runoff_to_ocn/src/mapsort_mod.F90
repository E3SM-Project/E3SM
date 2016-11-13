!===============================================================================
! SVN $Id: mapsort_mod.F90 35698 2012-03-22 23:59:57Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/tools/mapping/trunk_tags/mapping_141106/gen_mapping_files/runoff_to_ocn/src/mapsort_mod.F90 $
!===============================================================================

MODULE mapsort_mod

   use shr_sys_mod
   use shr_timer_mod
   use kind_mod
   use map_mod

   implicit none

   private ! default private

   public :: mapsort_checkSort  ! checks if map matrix is properly sorted
   public :: mapsort_sort       ! sort the map

   interface mapsort_sort; module procedure mapsort_byrow; end interface

   integer(IN),parameter :: debug = 1

!===============================================================================
CONTAINS
!===============================================================================

SUBROUTINE mapsort_checkSort(map,rowOK,colOK)

   implicit none

   !--- arguments ---
   type(sMatrix), intent(inout) :: map
   logical,optional,intent(out) :: rowOK,colOK

   !--- local ---
   integer(IN) :: n
   logical     :: sorted

   !--- formats ---
   character(*),parameter :: subName = "(mapsort_checkSort) "
   character(*),parameter :: F00 =   "('(mapsort_checkSort) ',8a)"

!-------------------------------------------------------------------------------
! PURPOSE:
! o checks whether elements of an sMatrix are properly sorted
!-------------------------------------------------------------------------------

   !--- are rows monotonically increasing? ---
   sorted = .true.
   do n=1,map%n_s-1
     if (map%row(n) > map%row(n+1) ) sorted = .false.
   enddo
   if (      sorted ) write(6,F00) "VERIFIED: rows are monotonically increasing"
   if (.not. sorted ) write(6,F00) "WARNING : rows are not monotonically increasing"
   if (present(rowOK)) rowOK = sorted

   !--- for a given row, are cols monotonically increasing? ---
   sorted = .true.
   do n=1,map%n_s-1
      if ((map%row(n) == map%row(n+1)) .and. (map%col(n) > map%col(n+1))) sorted = .false.
   enddo
   if (      sorted ) write(6,F00) "VERIFIED: cols are monotonically increasing"
   if (.not. sorted ) write(6,F00) "WARNING : cols are not monotonically increasing"
   if (present(colOK)) colOK = sorted

END SUBROUTINE mapsort_checkSort

!===============================================================================
!===============================================================================

SUBROUTINE mapsort_byrow(map)

   implicit none

   !--- arguments ---
   type(sMatrix), intent(inout) :: map

   !--- local ---
   integer(IN)             :: i,j          ! index into 2d matrix S(i,j)
   integer(IN)             :: n            ! index into 1d vector S(n)
   integer(IN)             :: k            ! index into tempS(i,k) = S(i,j) where j = tempCol(n,k)
   integer(IN)             :: na,nb        ! size of uncompressed matrix nb*na
   integer(IN)             :: ns           ! size of compressed matrix

   real   (R8),allocatable :: tempS  (:,:) ! matrix elements in temp non-standard format
   integer(IN),allocatable :: tempCol(:,:) ! matrix columns  in temp non-standard format
   integer(IN),allocatable :: nCols  (  :) ! number of active cols in each row of tempS(i,k)
   integer(IN)             :: kMax         ! size of k dim allocated/available in tempS(i,k)
   integer(IN)             :: count        ! counts occurances of things for debugging

   logical                 :: noSwap       ! flags col swapping when col sorting
   integer(IN)             :: nMin ,nMax   !
   integer(IN)             :: nMin2,nMax2  !
   integer(IN)             :: nMin0        !
   integer(IN),parameter   :: nMod=100000  ! for debug info
   integer(IN)             :: iTmp         ! temp var for col sorting
   real   (R8)             :: rTmp         ! temp var for col sorting
   logical                 :: done         ! col sorting is complete

   integer(IN)             :: nNew         ! new S(n) index wrt combining duplicate i/j

   character( 8)           :: dstr         ! wall clock date
   character(10)           :: tstr         ! wall clock time

#ifdef  _OPENMP
   integer(IN)              :: omp_get_max_threads !  $OMP function call
#endif

   !--- formats ---
   character(*),parameter :: subName = '(mapsort_byrow) '
   character(*),parameter :: F00   = "('(mapsort_byrow) ',8a)"
   character(*),parameter :: F01   = "('(mapsort_byrow) ',a,4i11)"
   character(*),parameter :: F12   = "('(mapsort_byrow) date & time:',1x,a4,2('-',a2),2x,a2,2(':',a2))"

!-------------------------------------------------------------------------------
! PURPOSE:
!   sort all elements of S so that row(n) <= row(n+1)
!   and if  row(n) == row(n+1), then col(n) <= col(n+1)
!-------------------------------------------------------------------------------

   ns = map%n_s
   na = map%n_a
   nb = map%n_b

   write(6,F01) "na, nb, ns = ",map%n_a,map%n_b,map%n_s
   write(6,F01) "min/max row = ",minval(map%row),maxval(map%row)
   write(6,F01) "min/max col = ",minval(map%col),maxval(map%col)

   !----------------------------------------------------------------------------
   write(6,F00) "fill a temp matrix (non-standard format) with row-sorted data"
   !----------------------------------------------------------------------------

   call date_and_time(dstr,tstr)
   write(6,F12) dstr(1:4),dstr(5:6),dstr(7:8) ,tstr(1:2),tstr(3:4),tstr(5:6)
   call shr_sys_flush(6)

   !--- determine how big temp matrix must be, allocate temp matrix ---
   allocate(nCols(nb))

   nCols = 0
   kMax = 0
   do n=1,ns
      i = map%row(n)
      nCols(i) = nCols(i)+1
      kMax = max(kMax,nCols(i))
   end do
   write(6,F01) "max number of cols per row = ",kMax

   allocate(tempS  (nb,kMax))
   allocate(tempCol(nb,kMax))

   !--- fill the temp matrix arrays with data ---
   nCols = 0
   do n=1,ns

      if (debug > 0  .and. mod(n,nMod) == 0) then
         write(6,F01) "... working, max elements, current element = ",ns,n
         call shr_sys_flush(6)
      endif

      i = map%row(n)  ! s(n) belongs in row i
      k = nCols(i)+1  ! s(n) is k-th non-zero element in row i

      if (nCols(i) > kMax) then
         write(6,F01) "ERROR: temp S,col arrays not large enough, kMax = ",kMax
         call shr_sys_abort()
      endif

      nCols  (i)   = k           ! number of elements in row i
      tempCol(i,k) = map%col(n)  ! associated col
      tempS  (i,k) = map%s  (n)  ! matrix element
   end do

   !----------------------------------------------------------------------------
   write(6,F00) "fill a input/output matrix with row-sorted data"
   !----------------------------------------------------------------------------

   call date_and_time(dstr,tstr)
   write(6,F12) dstr(1:4),dstr(5:6),dstr(7:8) ,tstr(1:2),tstr(3:4),tstr(5:6)
   call shr_sys_flush(6)

   n = 0
   do i=1,nb
   do k=1,nCols(i)
      n = n+1
      map%row(n) = i
      map%col(n) = tempCol(i,k)
      map%s  (n) = tempS  (i,k)
   end do
   end do

   deallocate(nCols,tempS,tempCol)

   call mapsort_checkSort(map) ! confirm rows are sorted

   !----------------------------------------------------------------------------
   write(6,F00) "sort columns within each row (assumes rows are sorted)"
   !----------------------------------------------------------------------------

   call date_and_time(dstr,tstr)
   write(6,F12) dstr(1:4),dstr(5:6),dstr(7:8) ,tstr(1:2),tstr(3:4),tstr(5:6)
   call shr_sys_flush(6)

   call map_Sn1(map) ! determine where rows begin & end in s(n)

#ifdef  _OPENMP
   n = omp_get_max_threads()
   write(6,F01) 'FYI: this routine is threaded, omp_get_max_threads() = ',n
#else
   write(6,F00) 'FYI: this routine is not threaded'
#endif

!$OMP PARALLEL DO PRIVATE(i,done,n,nMin,nMax,nMin0,nMin2,nMax2,iTmp,rTmp,noSwap,dstr,tstr)
   do i=1,nb
      if (debug > 0  .and. mod(i,nMod) == 0) then
         write(6,F01) "... working on row ",i
         call date_and_time(dstr,tstr)
         write(6,F12) dstr(1:4),dstr(5:6),dstr(7:8) ,tstr(1:2),tstr(3:4),tstr(5:6)
         call shr_sys_flush(6)
      endif

      nMin0 = map%sn2(i) + 1          ! 1st element in row i
      nMin  = nMin0                   ! 1st  potentially out-of-order element in row i
      nMax  = map%sn2(i) + map%sn1(i) ! last potentially out-of-order element in row i
      done = .false.
      do while ( .not. done )
         noSwap = .true.
         do n=nMin,nMax-1 ! all elements should be in the same row
         !  if (map%row(n) == map%row(n+1) .and. map%col(n) > map%col(n+1)) then
            if (map%col(n) > map%col(n+1)) then

               !--- swap adjacent out-of-order columns (within same row) ---
               iTmp         = map%col(n+1)
               map%col(n+1) = map%col(n)
               map%col(n)   = iTmp
               rTmp         = map%s(n+1)
               map%s(n+1)   = map%s(n)
               map%s(n)     = rTmp

               !--- note subset of out-of-order elements ---
               if (noSwap) then ! 1st swap
                  nMin2 = max(n-1,nMin0)
                  noSwap = .false.
               end if
               nMax2=n          ! last swap
            endif
         enddo
         if (noSwap) then
            done = .true.
         else ! narrow the search
            nMin=nMin2  ! 1st  potentially out-of-order element
            nMax=nMax2  ! last potentially out-of-order element
         endif
      enddo
   end do

   call mapsort_checkSort(map) ! confirm cols are sorted

   !----------------------------------------------------------------------------
   write(6,F00) "combine elements with same row/col (assumes matrix is sorted)"
   !----------------------------------------------------------------------------

   call date_and_time(dstr,tstr)
   write(6,F12) dstr(1:4),dstr(5:6),dstr(7:8) ,tstr(1:2),tstr(3:4),tstr(5:6)
   call shr_sys_flush(6)

   write(6,F01) "... orig size of S = ",map%n_s
   i = 0
   j = 0
   nNew = 0
   count = 0
   do n=1,map%n_s
      if (map%row(n) == i  .and.  map%col(n) == j) then
         !--- same i,j as previous element => combine them ---
         count = count + 1
         map%s(nNew) = map%s(nNew) + map%s(n)
      else
         !--- different i,j from previous element => don't combine ---
         nNew = nNew + 1
         map%row(nNew) = map%row(n)
         map%col(nNew) = map%col(n)
         map%s  (nNew) = map%s  (n)
         i = map%row(n)
         j = map%col(n)
      endif
   enddo
   map%n_s = nNew
   write(6,F01) "... new  size of S = ",map%n_s
   write(6,F01) "... number of duplicate row/col = ",count

   call map_Sn1(map) ! determine where rows begin & end in s(n)

   write(6,F00) "Done."
   call date_and_time(dstr,tstr)
   write(6,F12) dstr(1:4),dstr(5:6),dstr(7:8) ,tstr(1:2),tstr(3:4),tstr(5:6)
   call shr_sys_flush(6)

END SUBROUTINE mapsort_byrow

!===============================================================================
!===============================================================================

SUBROUTINE mapsort_sortMatRC(map)

   implicit none

   !--- arguments ---
   type(sMatrix), intent(inout) :: map

   !--- local ---
   integer(IN) :: i,m,n
   logical     :: rowOK,colOK
   integer(IN) :: t01 = 0         ! system-clock-timer number

#ifdef  _OPENMP
   integer(IN) :: omp_get_max_threads !  $OMP function call
#endif
   character(8)  :: dstr          ! F90 wall clock date string yyyymmdd
   character(10) :: tstr          ! F90 wall clock time string hhmmss.sss

   !--- formats ---
   character(*),parameter :: subName = "(mapsort_sortMatRC) "
   character(*),parameter :: F00   = "('(mapsort_sortMatRC) ',8a)"
   character(*),parameter :: F01   = "('(mapsort_sortMatRC) ',a,4i11)"
   character(*),parameter :: F10   = "('(mapsort_sortMatRC) ',a,i7, &
                                   &   '  date & time:',1x, &
                                   &   a4,2('-',a2),2x,a2,2(':',a2))"

!-------------------------------------------------------------------------------
! PURPOSE:
! o sorts the elements of an sMatrix to achieve uniformly increasing columns
!   within each given row
! o assumes the elements of an sMatrix already have uniformly increasing rows
!-------------------------------------------------------------------------------

   write(6,F01) "sort matrix, size = ",map%n_s

   if (t01 == 0) call shr_timer_get(t01,subName//"sort matrix")

   !----------------------------------------------------------------------------
   !check for correct row/column sorting
   !----------------------------------------------------------------------------
   call mapsort_checkSort(map)

   !----------------------------------------------------------------------------
   ! sort entire matrix
   !----------------------------------------------------------------------------
   write(6,F00) "sort entire matrix"
   call shr_timer_zero (t01)
   call shr_timer_start(t01)

   call date_and_time(dstr,tstr) ! F90 intrinsic routine
   write(6,F10) 'i=',0,dstr(1:4),dstr(5:6),dstr(7:8),tstr(1:2),tstr(3:4),tstr(5:6)
   call shr_sys_flush(6)

   n = 1
   m = map%n_s
   call mapsort_qsortRC(map%s(n:m),map%row(n:m),map%col(n:m))

   call shr_timer_stop(t01)
   call shr_timer_print(t01)

   !----------------------------------------------------------------------------
   ! verify matrix has been properly sorted
   !----------------------------------------------------------------------------
   call mapsort_checkSort(map,rowOK=rowOK,colOK=colOK)
   if (.not. rowOK ) write(6,F00) "ERROR: rows are NOT monotonically increasing"
   if (.not. colOK ) write(6,F00) "ERROR: cols are NOT monotonically increasing"
   if (.not. colOK .or. .not. rowOK) stop

END SUBROUTINE mapsort_sortMatRC

!===============================================================================
!===============================================================================

recursive SUBROUTINE mapsort_qsortRC(S,row,col)

   implicit none

   !--- arguments ---
   real   (R8) ::   S(:)
   integer(IN) :: row(:),col(:)

   !--- local ---
   integer(IN) :: ns  ! size of data array to sort
   integer(IN) :: m   ! marker that divides data array into two

   !--- formats ---
   character(*),parameter :: F00 = "('(mapsort_qsortRC) ',8a)"
   character(*),parameter :: F01 = "('(mapsort_qsortRC) ',a,5i6)"

!-------------------------------------------------------------------------------
! PURPOSE:
!   sort all elements of S so that col(n) is strickly increasing
!   assumes that row(n) is a constant (ie. row(n) = row(m) for every n,m)
!-------------------------------------------------------------------------------

   ns = size(S)
   if ( ns > 1) then
      call mapsort_partitianRC(S,row,col,m)
      call mapsort_qsortRC   (S(1:m-1),row(1:m-1),col(1:m-1))
      call mapsort_qsortRC   (S(m:ns ),row(m:ns ),col(m:ns ))
   end if

END SUBROUTINE mapsort_qsortRC

!===============================================================================
!===============================================================================

SUBROUTINE mapsort_partitianRC(S,row,col,marker)

   implicit none

   !--- arguments ---
   real   (R8) :: S(:)     ! must rearrange along with col values
   integer(IN) :: row(:)   ! unused: assumed constant
   integer(IN) :: col(:)   ! sort wrt column value
   integer(IN) :: marker   ! partitian wrt this index

   !--- local ---
   integer(IN) :: i,j      ! generic index
   integer(IN) :: ns       ! size of S(:),row(:),col(:)
   integer(IN) :: pivotR   ! partitian wrt this row value...
   integer(IN) :: pivotC   ! ... and   wrt this col value
   real   (R8) :: rTemp    ! temp var for swapping
   integer(IN) :: iTemp    ! temp var for swapping

   !--- formats ---
   character(*),parameter :: F00 = "('(mapsort_partitianRC) ',8a)"
   character(*),parameter :: F01 = "('(mapsort_partitianRC) ',a,2i11)"

!-------------------------------------------------------------------------------
! PURPOSE:
! o DOES NOT assume rows are sorted, assumes a completely unsorted matrix
! o rearrange elements in S & col and compute marker such that
!   if  n = marker  and  i <= n <= j   we have col(i) <= col(n) <= col(j)
!-------------------------------------------------------------------------------

   ns = size(S)
   if (ns < 2) then
      write(6,F01) "ERROR: ns < 2, ns = ",ns
      stop
   end if

   pivotC = col(ns/2) ! note: fortran truncation
   pivotR = row(ns/2) ! note: fortran truncation
!  pivotC = col(1)
!  pivotR = row(1)

   i = 0
   j = ns+1
   do
      !--- move j downward ---
      j = j-1
      do
        if  (row(j) <  pivotR)                        exit
        if  (row(j) == pivotR .and. col(j) <= pivotC) exit
        j = j-1
      end do
      !--- move i upward ---
      i = i+1
      do
        if  (row(j) >  pivotR)                         exit
        if  (row(j) == pivotR .and. col(j) >= pivotC ) exit
        i = i+1
      end do
      !--- if (i<j) then swap & continue, else return ---
      if ( i < j ) then
         !--- swap matrix elements ---
         rTemp  = S  (i)
         S  (i) = S  (j)
         S  (j) = rTemp
         !--- swap column index ---
         iTemp  = col(i)
         col(i) = col(j)
         col(j) = iTemp
         !--- swap row    index ---
         iTemp  = row(i)
         row(i) = row(j)
         row(j) = iTemp
      else if (i == j) then
         marker = i+1
         return
      else
         marker = i
         return
      end if
   end do

END SUBROUTINE mapsort_partitianRC

!===============================================================================
!===============================================================================
!===============================================================================

SUBROUTINE mapsort_sortMatC(map)

   implicit none

   !--- arguments ---
   type(sMatrix), intent(inout) :: map

   !--- local ---
   integer(IN) :: i,m,n
   logical     :: rowOK,colOK
   integer(IN) :: t01 = 0         ! system-clock-timer number

#ifdef  _OPENMP
   integer(IN) :: omp_get_max_threads !  $OMP function call
#endif
   character(8)  :: dstr          ! F90 wall clock date string yyyymmdd
   character(10) :: tstr          ! F90 wall clock time string hhmmss.sss

   !--- formats ---
   character(*),parameter :: subName = "(mapsort_sortMatC) "
   character(*),parameter :: F00   = "('(mapsort_sortMatC) ',8a)"
   character(*),parameter :: F01   = "('(mapsort_sortMatC) ',a,4i11)"
   character(*),parameter :: F10   = "('(mapsort_sortMatC) ',a,i7, &
                                   &   '  date & time:',1x, &
                                   &   a4,2('-',a2),2x,a2,2(':',a2))"

!-------------------------------------------------------------------------------
! PURPOSE:
! o sorts the elements of an sMatrix to achieve uniformly increasing columns
!   within each given row
! o assumes the elements of an sMatrix already have uniformly increasing rows
!-------------------------------------------------------------------------------

   write(6,F01) "sort matrix, size = ",map%n_s

   if (t01 == 0) then
      call shr_timer_get(t01,subName//"sort matrix")
#ifdef  _OPENMP
      n = omp_get_max_threads()
      write(6,F01) 'FYI: this routine is threaded, omp_get_max_threads() = ',n
#else
      write(6,F01) 'FYI: this routine is not threaded'
#endif
   end if

   !----------------------------------------------------------------------------
   !check for correct row/column sorting
   !----------------------------------------------------------------------------

   call mapsort_checkSort(map,rowOK=rowOK)
   if (.not. rowOK ) then
      write(6,F00) "ERROR: rows are NOT monotonically increasing"
      stop
   end if

   !----------------------------------------------------------------------------
   ! sort the columns within each row
   !----------------------------------------------------------------------------
   write(6,F00) "sort the columns within each row of a matrix"
   write(6,F01) "min/max Sn1 = ",minval(map%Sn1),maxval(map%Sn1)
   write(6,F01) "min/max Sn2 = ",minval(map%Sn2),maxval(map%Sn2)
   call shr_timer_zero (t01)
   call shr_timer_start(t01)

   write(6,F01) "loop: i=1,nRows    where nRows = ",map%n_b
   call date_and_time(dstr,tstr) ! F90 intrinsic routine
   write(6,F10) 'i=',0,dstr(1:4),dstr(5:6),dstr(7:8),tstr(1:2),tstr(3:4),tstr(5:6)
   call shr_sys_flush(6)

!$OMP PARALLEL DO PRIVATE(i,n,m)
   do i=1,map%n_b
      !--- sort the columns in row i ------------------------
      if (mod(i,1000)==0) then
         call date_and_time(dstr,tstr) ! F90 intrinsic routine
         write(6,F10) 'i=',i,dstr(1:4),dstr(5:6),dstr(7:8),tstr(1:2),tstr(3:4),tstr(5:6)
         call shr_sys_flush(6)
      end if
      !------------------------------------------------------
      if ( map%sn1(i) > 1) then
         n = map%sn2(i) + 1          ! S(n) is 1st  element in row i
         m = map%sn2(i) + map%sn1(i) ! S(m) is last element in row i
         call mapsort_QsortC(map%s(n:m),map%row(n:m),map%col(n:m))
      end if
      !------------------------------------------------------
   end do
   call shr_timer_stop(t01)
   call shr_timer_print(t01)

   !----------------------------------------------------------------------------
   ! verify matrix has been properly sorted
   !----------------------------------------------------------------------------
   call mapsort_checkSort(map,rowOK=rowOK,colOK=colOK)
   if (.not. rowOK ) write(6,F00) "ERROR: rows are NOT monotonically increasing"
   if (.not. colOK ) write(6,F00) "ERROR: cols are NOT monotonically increasing"
   if (.not. colOK .or. .not. rowOK) stop

END SUBROUTINE mapsort_sortMatC

!===============================================================================
!===============================================================================

recursive SUBROUTINE mapsort_QsortC(S,row,col)

   implicit none

   !--- arguments ---
   real   (R8) ::   S(:)
   integer(IN) :: row(:),col(:)

   !--- local ---
   integer(IN) :: ns  ! size of data array to sort
   integer(IN) :: m   ! marker that divides data array into two

   !--- formats ---
   character(*),parameter :: F00 = "('(mapsort_QsortC) ',8a)"
   character(*),parameter :: F01 = "('(mapsort_QsortC) ',a,5i6)"

!-------------------------------------------------------------------------------
! PURPOSE:
!   sort all elements of S so that col(n) is strickly increasing
!   assumes that row(n) is a constant (ie. row(n) = row(m) for every n,m)
!-------------------------------------------------------------------------------

   ns = size(S)
   if ( ns > 1) then
     call mapsort_partitianC(S,row,col,m)
     call mapsort_QsortC   (S(1:m-1),row(1:m-1),col(1:m-1))
     call mapsort_QsortC   (S(m:ns ),row(m:ns ),col(m:ns ))
   end if

END SUBROUTINE mapsort_QsortC

!===============================================================================
!===============================================================================

SUBROUTINE mapsort_partitianC(S,row,col,marker)

   implicit none

   !--- arguments ---
   real   (R8) :: S(:)     ! must rearrange along with col values
   integer(IN) :: row(:)   ! unused: assumed constant
   integer(IN) :: col(:)   ! sort wrt column value
   integer(IN) :: marker   ! partitian wrt this index

   !--- local ---
   integer(IN) :: i,j      ! generic index
   integer(IN) :: ns       ! size of S(:),row(:),col(:)
   integer(IN) :: pivot    ! partitian wrt this value
   real   (R8) :: rTemp    ! temp var for swapping
   integer(IN) :: iTemp    ! temp var for swapping

   !--- formats ---
   character(*),parameter :: F00 = "('(mapsort_partitianC) ',8a)"
   character(*),parameter :: F01 = "('(mapsort_partitianC) ',a,2i11)"

!-------------------------------------------------------------------------------
! PURPOSE:
! o assumes rows are already sorted
! o rearrange elements in S & col and compute marker such that
!   if  n = marker  and  i <= n <= j   we have col(i) <= col(n) <= col(j)
!-------------------------------------------------------------------------------

   ns = size(S)
   if (ns < 2) then
      write(6,F01) "ERROR: ns < 2, ns = ",ns
      stop
   end if
   if (minval(row) /= maxval(row)) then
      write(6,F01) "ERROR: all elements not in same row"
      write(6,F01) "ERROR: row min/max = ",minval(row),maxval(row)
      stop
   end if

!  pivot = col(ns/2) ! note: fortran truncation
   pivot = col(1)

   i = 0
   j = ns+1
   do
      !--- move j downward ---
      j = j-1
      do
        if (col(j) <= pivot) exit
        j = j-1
      end do
      !--- move i upward ---
      i = i+1
      do
        if (col(i) >= pivot) exit
        i = i+1
      end do
      !--- if (i<j) then swap & continue, else return ---
      if ( i < j ) then
         rTemp  = S  (i)
         S  (i) = S  (j)
         S  (j) = rTemp
         iTemp  = col(i)
         col(i) = col(j)
         col(j) = iTemp
      else if (i == j) then
         marker = i+1
         return
      else
         marker = i
         return
      end if
   end do

END SUBROUTINE mapsort_partitianC

!===============================================================================
!===============================================================================

SUBROUTINE mapsort_sortH(map)

   implicit none

   !--- arguments ---
   type(sMatrix), intent(inout) :: map

   !--- local ---
   integer n,thisrow,thiscol,n0,n1
   integer imin,imax,cnt,itmp,imin2,imax2
   double precision rtmp

!-------------------------------------------------------------------------------
! PURPOSE:
! o  subroutinized version of /fs/cgd/home0/hecht/csm/runoff/sort_map.F
! o  sorts the links in sMatrix "map" to achieve a uniformly increasing
!       map%row array.
! o  sorts the links in sMatrix "map" to achieve uniformly increasing
!       map%col for a given map%row value.
! o  assures that the array of links contains unique row-col pairs.
!-------------------------------------------------------------------------------

   imin=1
   imax=map%n_s
   cnt=map%n_s

   print *,'(mapsort_sortH): starting sort...'
   do while (cnt.gt.0)
   cnt=0
   do n=imin,imax-1
     if (map%row(n).gt.map%row(n+1)) then
       itmp=map%row(n+1)
       map%row(n+1)=map%row(n)
       map%row(n)=itmp
       itmp=map%col(n+1)
       map%col(n+1)=map%col(n)
       map%col(n)=itmp
       rtmp=map%s(n+1)
       map%s(n+1)=map%s(n)
       map%s(n)=rtmp
       cnt=cnt+1
       if (cnt.eq.1) imin2=n-1
       imax2=n
     endif
   enddo
   imin=max(imin2,1)
   imax=imax2
   enddo

   print *,'(mapsort_sortH): halfway done...'

   imin=1
   imax=map%n_s
   cnt=map%n_s

   do while (cnt.gt.0)
   cnt=0
   do n=imin,imax-1
     if (map%col(n).gt.map%col(n+1).and.map%row(n).eq.map%row(n+1)) then
       itmp=map%col(n+1)
       map%col(n+1)=map%col(n)
       map%col(n)=itmp
       rtmp=map%s(n+1)
       map%s(n+1)=map%s(n)
       map%s(n)=rtmp
       cnt=cnt+1
       if (cnt.eq.1) imin2=n-1
       imax2=n
     endif
   enddo
   imin=max(imin2,1)
   imax=imax2
   enddo
   print *,'(mapsort_sortH): done sorting for both row AND col'

   thisrow = 0
   thiscol = 0
   n0 = 0
   n1 = 0
   do n=1,map%n_s
     if (map%row(n).eq.thisrow.and.map%col(n).eq.thiscol) then
       n1 = n1 + 1
       map%s(n0) = map%s(n0)+map%s(n)
     else
       n0 = n0 + 1
       thisrow = map%row(n)
       thiscol = map%col(n)
       map%row(n0) = map%row(n)
       map%col(n0) = map%col(n)
       map%s(n0) = map%s(n)
     endif
   enddo
   map%n_s = n0
   print *,'(mapsort_sortH): # of repeat col->row links found ',n1
   print *,'(mapsort_sortH): now n_s = ',n0

END SUBROUTINE mapsort_sortH

!===============================================================================

SUBROUTINE mapsort_sortY(map)

   implicit none

   !--- arguments ---
   type(sMatrix), intent(inout) :: map

   !--- local ---
   integer i,j,n,thisrow,thiscol,n0,n1,ncnt,nskip,nrow
   integer imin,imax,cnt,itmp,imin2,imax2
   integer, allocatable :: rowtmp(:),coltmp(:),stmp(:)
   logical rowok,colok
   double precision rtmp

!-------------------------------------------------------------------------------
! PURPOSE:
! o  subroutinized version of /fs/cgd/home0/hecht/csm/runoff/sort_map.F
! o  sorts the links in sMatrix "map" to achieve a uniformly increasing
!       map%row array.
! o  sorts the links in sMatrix "map" to achieve uniformly increasing
!       map%col for a given map%row value.
! o  assures that the array of links contains unique row-col pairs.
!-------------------------------------------------------------------------------

   allocate(stmp(map%n_s))
   allocate(rowtmp(map%n_s))
   allocate(coltmp(map%n_s))
   stmp = map%s
   rowtmp = map%row
   coltmp = map%col

!  check for correct row/column sorting  (copied from map_check)
   rowok = .true.
   colok = .true.
!  is row monotonically increasing?
   do n=2,map%n_s
     if (map%row(n-1) > map%row(n) ) rowok = .false.
   enddo
!  is col monotonically increasing for a given row?
   do i=1,map%n_s-1
      if (map%col(i).gt.map%col(i+1).and.map%row(i+1).eq.map%row(i)) then
         colok = .false.
      endif
   enddo

!  if (.not. rowok .and. .not. colok) then
!    imin=1
!    imax=map%n_s
!    cnt=map%n_s
!    print *,'(mapsort_sortY): slow-sorting rows...'
!    do while (cnt.gt.0)
!    cnt=0
!    do n=imin,imax-1
!      if (map%row(n).gt.map%row(n+1)) then
!        itmp=map%row(n+1)
!        map%row(n+1)=map%row(n)
!        map%row(n)=itmp
!        itmp=map%col(n+1)
!        map%col(n+1)=map%col(n)
!        map%col(n)=itmp
!        rtmp=map%s(n+1)
!        map%s(n+1)=map%s(n)
!        map%s(n)=rtmp
!        cnt=cnt+1
!        if (cnt.eq.1) imin2=n-1
!        imax2=n
!      endif
!    enddo
!    imin=max(imin2,1)
!    imax=imax2
!    enddo
!  endif

!  if (.not. rowok .and. colok) then
!    cnt=1
!    print *,'(mapsort_sortY): fast-sorting rows...'
!    do j=1,map%n_b
!      n0 = 0
!      n1 = 0
!      do n=cnt,map%n_s
!        if (map%row(n).eq.j .and. n0.eq.0) n0=n
!        if (map%row(n).ne.j .and. n0.ne.0 .and. n1.eq.0) n1=n-1
!      enddo
!      if (n0.ne.0) then
!        if (n1.eq.0) n1=n0
!        nrow = n1-n0+1
!        nskip = n0-cnt
!        map%row(cnt:cnt+nrow-1)=rowtmp(n0:n1)
!        map%col(cnt:cnt+nrow-1)=coltmp(n0:n1)
!        map%s(cnt:cnt+nrow-1)=stmp(n0:n1)
!        if (nskip.gt.0) then
!          map%row(cnt+nrow:cnt+nrow+nskip-1) = rowtmp(cnt:n0-1)
!          map%col(cnt+nrow:cnt+nrow+nskip-1) = coltmp(cnt:n0-1)
!          map%s(cnt+nrow:cnt+nrow+nskip-1) = stmp(cnt:n0-1)
!        endif
!        cnt = cnt+nrow
!        if (n1.lt.map%n_s) then
!          if (cnt+nskip .ne. n1+1) write(6,*) 'ERROR cnt+nskip .ne. n1+1'
!          map%row(n1+1:map%n_s) = rowtmp(n1+1:map%n_s)
!          map%col(n1+1:map%n_s) = coltmp(n1+1:map%n_s)
!          map%s(n1+1:map%n_s) = stmp(n1+1:map%n_s)
!        endif
!        rowtmp = map%row
!        coltmp = map%col
!        stmp = map%s
!      endif
!    enddo
!  endif

   if (.not. rowok) then
     call mapsort_byrow(map)
   endif

   if (.not. colok) then
     imin=1
     imax=map%n_s
     cnt=map%n_s
     print *,'(mapsort_sortY): sorting columns...'
     do while (cnt.gt.0)
     cnt=0
     do n=imin,imax-1
       if (map%col(n).gt.map%col(n+1).and.map%row(n).eq.map%row(n+1)) then
         itmp=map%col(n+1)
         map%col(n+1)=map%col(n)
         map%col(n)=itmp
         rtmp=map%s(n+1)
         map%s(n+1)=map%s(n)
         map%s(n)=rtmp
         cnt=cnt+1
         if (cnt.eq.1) imin2=n-1
         imax2=n
       endif
     enddo
     imin=max(imin2,1)
     imax=imax2
     enddo
   endif

   print *,'(mapsort_sortY): done sorting for both row AND col'

   thisrow = 0
   thiscol = 0
   n0 = 0
   n1 = 0
   do n=1,map%n_s
     if (map%row(n).eq.thisrow.and.map%col(n).eq.thiscol) then
       n1 = n1 + 1
       map%s(n0) = map%s(n0)+map%s(n)
     else
       n0 = n0 + 1
       thisrow = map%row(n)
       thiscol = map%col(n)
       map%row(n0) = map%row(n)
       map%col(n0) = map%col(n)
       map%s(n0) = map%s(n)
     endif
   enddo
   map%n_s = n0
   print *,'(mapsort_sortY): # of repeat col->row links found ',n1
   print *,'(mapsort_sortY): now n_s = ',n0

END SUBROUTINE mapsort_sortY

!===============================================================================
!===============================================================================

SUBROUTINE mapsort_bubble(map)

   implicit none

   !--- arguments ---
   type(sMatrix), intent(inout) :: map

   !--- local ---
   integer(IN) :: k,m,n          !
   logical     :: rowOK,colOK    !
   integer(IN) :: t01 = 0        ! system-clock-timer number
   real   (R8) :: rTemp          ! temp var for swapping
   integer(IN) :: iTemp          ! temp var for swapping
   real   (R8),pointer :: S  (:) !
   integer(IN),pointer :: row(:) !
   integer(IN),pointer :: col(:) !

   character(8)  :: dstr          ! F90 wall clock date string yyyymmdd
   character(10) :: tstr          ! F90 wall clock time string hhmmss.sss

   !--- formats ---
   character(*),parameter :: subName = "(mapsort_bubble) "
   character(*),parameter :: F00   = "('(mapsort_bubble) ',8a)"
   character(*),parameter :: F01   = "('(mapsort_bubble) ',a,4i11)"
   character(*),parameter :: F10   = "('(mapsort_bubble) ',a,i7, &
                                   &   '  date & time:',1x, &
                                   &   a4,2('-',a2),2x,a2,2(':',a2))"

!-------------------------------------------------------------------------------
! PURPOSE:
!-------------------------------------------------------------------------------

   write(6,F01) "sort entire matrix, size = ",map%n_s

   if (t01 == 0) call shr_timer_get(t01,subName//"sort matrix")

   !----------------------------------------------------------------------------
   !check for correct row/column sorting
   !----------------------------------------------------------------------------
   call mapsort_checkSort(map)

   !----------------------------------------------------------------------------
   ! sort entire matrix
   !----------------------------------------------------------------------------
   call shr_timer_zero (t01)
   call shr_timer_start(t01)

   S   => map%S
   row => map%row
   col => map%col

   do n = 1, map%n_s - 1
      !--- progress report ----------------------------------
      if (mod(n-1,10000)==0) then
         call date_and_time(dstr,tstr) ! F90 intrinsic routine
         write(6,F10) 'n=',n,dstr(1:4),dstr(5:6),dstr(7:8),tstr(1:2),tstr(3:4),tstr(5:6)
         call shr_sys_flush(6)
      end if
      !--- find desired element ---------------------------
      k = n
      do m = n+1, map%n_s
         if      (row(m) <  row(k)) then
            k = m
         else if (row(m) == row(k) .and. col(m) < col(k)) then
            k = m
         end if
      end do
      !--- put desired element into place -----------------
      if ( n /= k ) then
         !--- swap matrix elements ---
         rTemp  = S  (n)
         S  (n) = S  (k)
         S  (k) = rTemp
         !--- swap column index ---
         iTemp  = col(n)
         col(n) = col(k)
         col(k) = iTemp
         !--- swap row    index ---
         iTemp  = row(n)
         row(n) = row(k)
         row(k) = iTemp
      end if
   end do

   call shr_timer_stop (t01)
   call shr_timer_print(t01)

   !----------------------------------------------------------------------------
   ! verify matrix has been properly sorted
   !----------------------------------------------------------------------------
   call mapsort_checkSort(map,rowOK=rowOK,colOK=colOK)
   if (.not. rowOK ) write(6,F00) "ERROR: rows are NOT monotonically increasing"
   if (.not. colOK ) write(6,F00) "ERROR: cols are NOT monotonically increasing"
   if (.not. colOK .or. .not. rowOK) stop

END SUBROUTINE mapsort_bubble

!===============================================================================
!===============================================================================

END MODULE mapsort_mod
