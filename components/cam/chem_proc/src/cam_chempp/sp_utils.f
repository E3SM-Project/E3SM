
      recursive subroutine STCO( vertex )
!----------------------------------------------------------------
!	... Permutate sparse matrix in order to find the
!           the strongly connected blocks in a lower triangular
!           block form (LTBF)
!----------------------------------------------------------------

      use SP_MODS

      implicit none

!----------------------------------------------------------------
!	... Dummy args
!----------------------------------------------------------------
      integer, intent(in) :: vertex

!----------------------------------------------------------------
!	... Local variables
!----------------------------------------------------------------
      integer :: i, j, vn
      integer :: adj_vertex

      nb = nb + 1
      lowlink(vertex) = nb
      number(vertex) = nb
      sp = sp + 1
      vstack(sp) = vertex
      do vn = rp(vertex),rp(vertex+1)-1
	 adj_vertex = ci(vn)
	 if( number(adj_vertex) == 0 ) then
	    call STCO( adj_vertex )
	    lowlink(vertex) = MIN( lowlink(vertex),lowlink(adj_vertex) )
	 else if( number(adj_vertex) < number(vertex) ) then
            do j = 1,sp
	       if( vstack(j) == adj_vertex ) then
		  lowlink(vertex) = MIN( lowlink(vertex),number(adj_vertex) )
		  exit
	       end if
	    end do
	 end if
      end do

      if( lowlink(vertex) == number(vertex) ) then
	 blkcnt = blkcnt + 1
	 stcoblk(blkcnt) = pp + 1
	 do
	    if( sp == 0 ) then
	       exit
	    end if
	    vn = vstack(sp)
	    if( number(vn) >= number(vertex) ) then
	       sp = sp - 1
	       pp = pp + 1
	       perm(pp) = vn
	       blkmemcnt(blkcnt) = blkmemcnt(blkcnt) + 1
	    else
	       exit
	    end if
	 end do
      end if

      end subroutine STCO

      recursive subroutine STCO_PAT( vertex, order )
!----------------------------------------------------------------
!	... Permutate sparse matrix in order to find the
!           the strongly connected blocks in a lower triangular
!           block form (LTBF)
!----------------------------------------------------------------

      use SP_MODS

      implicit none

!----------------------------------------------------------------
!	... Dummy args
!----------------------------------------------------------------
      integer, intent(in) :: vertex, order

!----------------------------------------------------------------
!	... Local variables
!----------------------------------------------------------------
      integer :: i, j, vn
      integer :: adj_vertex

      nb = nb + 1
      lowlink(vertex) = nb
      number(vertex) = nb
      sp = sp + 1
      vstack(sp) = vertex
      do vn = 1,order
	 if( matrix(vertex,vn) ) then
	    adj_vertex = vn
	    if( number(adj_vertex) == 0 ) then
	       call STCO_PAT( adj_vertex, order )
	       lowlink(vertex) = MIN( lowlink(vertex),lowlink(adj_vertex) )
	    else if( number(adj_vertex) < number(vertex) ) then
               do j = 1,sp
	          if( vstack(j) == adj_vertex ) then
		     lowlink(vertex) = MIN( lowlink(vertex),number(adj_vertex) )
		     exit
	          end if
	       end do
	    end if
	 end if
      end do

      if( lowlink(vertex) == number(vertex) ) then
	 blkcnt = blkcnt + 1
	 stcoblk(blkcnt) = pp + 1
	 do
	    if( sp == 0 ) then
	       exit
	    end if
	    vn = vstack(sp)
	    if( number(vn) >= number(vertex) ) then
	       sp = sp - 1
	       pp = pp + 1
	       perm(pp) = vn
	       blkmemcnt(blkcnt) = blkmemcnt(blkcnt) + 1
	    else
	       exit
	    end if
	 end do
      end if

      end subroutine STCO_PAT

      subroutine DIAG_MARK( order, matrix, perm )
!---------------------------------------------------------------------------
!	... Find permuatation via diagonal markowitz to produce
!           near optimal LU fillin
!---------------------------------------------------------------------------

      implicit none

!---------------------------------------------------------------------------
!	... Dummy args
!---------------------------------------------------------------------------
      integer, intent(in)  :: order
      integer, intent(out) :: perm(order)
      logical, intent(in)  :: matrix(order,order)

!---------------------------------------------------------------------------
!	... Local variables
!---------------------------------------------------------------------------
      integer :: row, rowp1, col, maxrow
      integer :: beta, cnt
      integer :: i, j
      logical :: holder(order)
      logical :: pattern(order,order)

      cnt = COUNT( matrix )
      perm(:order) = (/ (row,row=1,order) /)        ! no permutations
      pattern = matrix
      do row = 1,order-1
	 rowp1 = row + 1
	 beta = (order - 1)**2
	 maxrow = row
	 do col = row,order
	    cnt = (COUNT(pattern(col,row:order)) - 1) * &
	          (COUNT(pattern(row:order,col)) - 1)
	    if( cnt < beta ) then
	       beta = cnt
	       maxrow = col
	    end if
	 end do
	 if( maxrow /= row ) then
!----------------------------------------------------------------
!	... Row and column permuatation
!----------------------------------------------------------------
	    holder(row:order) = pattern(row,row:order)
	    pattern(row,row:order) = pattern(maxrow,row:order)
	    pattern(maxrow,row:order) = holder(row:order)
	    holder(row:order) = pattern(row:order,row)
	    pattern(row:order,row) = pattern(row:order,maxrow)
	    pattern(row:order,maxrow) = holder(row:order)
	    beta = perm(row)
	    perm(row) = perm(maxrow)
	    perm(maxrow) = beta
	 end if
!----------------------------------------------------------------
!	... Now do "symbolic" decomposition on sub-matrix
!----------------------------------------------------------------
	 do col = rowp1,order
	    if( pattern(row,col) ) then
	       pattern(rowp1:order,col) = pattern(rowp1:order,row) .or. &
					  pattern(rowp1:order,col)
	    end if
	 end do
      end do

      end subroutine DIAG_MARK

      subroutine SYM_FAC( order, matrix, adds, mults )
!----------------------------------------------------------------
!	... Do a symbolic LU decomposition on sparse matrix
!----------------------------------------------------------------

      implicit none

!----------------------------------------------------------------
!	... Dummy args
!----------------------------------------------------------------
      integer, intent(in) :: order                     ! matrix order
      logical, intent(inout) :: matrix(order,order)    ! sparse matrix
      integer, intent(out) :: adds(2), mults(2)        ! operation counts

!----------------------------------------------------------------
!	... Local variables
!----------------------------------------------------------------
      integer :: i, im1, j, rcnt, ccnt

      adds = 0; mults = 0;
      do i = 2,order
	 im1 = i - 1
	 ccnt = COUNT( matrix(i:order,im1) )
	 rcnt = COUNT( matrix(im1,i:order) )
	 mults(1) = mults(1) + ccnt * (1 + rcnt)
	 mults(2) = mults(2) + ccnt
	 do j = i,order
	    if( matrix(im1,j) ) then
	       adds(1) = adds(1) + COUNT( matrix(i:order,j) .and. matrix(i:order,im1) )
	       matrix(i:order,j) = matrix(i:order,j) .or. matrix(i:order,im1)
	    end if
	 end do
      end do

      do i = order,2,-1
	 ccnt = COUNT( matrix(1:i-1,i) )
	 mults(2) = mults(2) + ccnt
      end do
      adds(2) = mults(2)

      end subroutine SYM_FAC

      subroutine GEN_PAT( order, pattern, rp, ci )
!----------------------------------------------------------------
!	... Generate sparsity pattern from row storage
!----------------------------------------------------------------

      implicit none
!----------------------------------------------------------------
!	... Dummy args
!----------------------------------------------------------------
      integer, intent(in) :: order
      integer, intent(in) :: rp(order+1)
      integer, intent(in) :: ci(*)
      logical, intent(out) :: pattern(order,order)

!----------------------------------------------------------------
!	... Local variables
!----------------------------------------------------------------
      integer :: i, row

      pattern = .false.
      do i = 1,order
	 do row = rp(i),rp(i+1)-1
	    pattern(i,ci(row)) = .true.
	 end do
      end do

      end subroutine GEN_PAT

      subroutine PERM_MAT( order, matrix, perm)
!----------------------------------------------------------------
!	... Matrix row and column permutation
!----------------------------------------------------------------

      implicit none

!----------------------------------------------------------------
!	... Dummy args
!----------------------------------------------------------------
      integer, intent(in) :: order
      integer, intent(in) :: perm(order)
      logical, intent(inout) :: matrix(order,order)

!----------------------------------------------------------------
!	... Local variables
!----------------------------------------------------------------
      integer :: i
      logical :: copy(order,order)

      copy = matrix
!----------------------------------------------------------------
!	... Row permutation
!----------------------------------------------------------------
      do i = 1,order
	 if( perm(i) /= i ) then
	    copy(i,:order) = matrix(perm(i),:order)
	 end if
      end do
      matrix = copy
!----------------------------------------------------------------
!	... Col permutation
!----------------------------------------------------------------
      do i = 1,order
	 if( perm(i) /= i ) then
	    matrix(:order,i) = copy(:order,perm(i))
	 end if
      end do

      end subroutine PERM_MAT

      subroutine DRAW_MAT( order, matrix )
!----------------------------------------------------------------
!	... Draw matrix sparsity
!----------------------------------------------------------------

      implicit none

!----------------------------------------------------------------
!	... Dummy args
!----------------------------------------------------------------
      integer, intent(in) :: order
      logical, intent(in) :: matrix(order,order)

!----------------------------------------------------------------
!	... Local variables
!----------------------------------------------------------------
      integer :: i, row
      character(len=1)    :: matrow(order)
      character(len=32)   :: frmt
      character(len=132)  :: line

      write(*,*) ' '
      line = ' '
      do i = 10,order,10
	 write(line(5+i*2:5+i*2),'(i1)') i/10
      end do
      write(*,'(a)') line(:LEN_TRIM(line))
      line = ' '
      do i = 1,order
	 write(line(5+i*2:5+i*2),'(i1)') MOD( i,10 )
      end do
      write(*,'(a)') line(:LEN_TRIM(line))
      write(*,*) ' '
      frmt = '(1x,i2,2x,'
      i = LEN_TRIM( frmt ) + 1
      if( order < 10 ) then
	 write(frmt(i:i),'(i1)') order
      else
	 write(frmt(i:i+1),'(i2)') order
      end if
      frmt(LEN_TRIM(frmt)+1:) = '(1x,a))'
      do row = 1,order
	 matrow = ' '
	 where( matrix(row,:order) )
	    matrow(:order) = 'X'
	 endwhere
	 write(*,frmt) row, matrow
      end do

      end subroutine DRAW_MAT
