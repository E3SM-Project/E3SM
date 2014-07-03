!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_sparse_type.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glimmer_sparse_type
 
  use glimmer_global, only:dp
  implicit none
  
  !*FD sparse matrix type
  type sparse_matrix_type
     integer :: nonzeros  !*FD number of nonzero elements currently stored
     integer :: order     !*FD order of the matrix (e.g. number of rows)
     logical :: symmetric !*FD True only if triangle of the symmetric matrix stored
     integer, dimension(:), pointer :: col => NULL()        !*FD column index
     integer, dimension(:), pointer :: row => NULL()        !*FD row index
     real(kind=dp), dimension(:), pointer :: val => NULL()  !*FD values
  end type sparse_matrix_type

  type sparse_solver_options_base
        real(kind=dp) :: tolerance !*FD Error tolerance
        integer :: maxiters        !*FD Max iterations before giving up
        integer :: method
  end type

  ! size of sparse matrix 
  integer, parameter, private :: chunksize=1000

  !MAKE_RESTART
!EIB!#ifdef RESTARTS
!EIB!#define RST_GLIMMER_SPARSE
!EIB!#include "glimmer_rst_head.inc"
!EIB!#undef RST_GLIMMER_SPARSE
!EIB!#endif

contains

!EIB!#ifdef RESTARTS
!EIB!#define RST_GLIMMER_SPARSE
!EIB!#include "glimmer_rst_body.inc"
!EIB!#undef RST_GLIMMER_SPARSE
!EIB!#endif

  subroutine new_sparse_matrix(order,n,mat)
    !*FD create a new sparse matrix
    implicit none
    integer, intent(in) :: n          !*FD initial size of matrix (non-zeros)
    type(sparse_matrix_type) :: mat   !*FD matrix
    integer, intent(in) :: order      !*FD Order (number of rows and columns) of the matrix
    
    if (.not.associated(mat%col)) then
       allocate(mat%row(n))
       !SLAP's sparse column scheme looks past the assumed bounds of col to see
       !what sparse storage format we're in.  To avoid array bounds problems, we
       !add 2 to the column size.  See mailing list discussion at:
       !http://forge.nesc.ac.uk/pipermail/glimmer-discuss/2005-February/000078.html
       allocate(mat%col(n+2))
       allocate(mat%val(n))
    else
       if (size(mat%row).lt.n) then
          call del_sparse_matrix(mat)
          allocate(mat%row(n))
          allocate(mat%col(n+2))
          allocate(mat%val(n))
       end if
    end if
    mat%nonzeros = 0
    mat%order = order
    mat%symmetric = .false.
  end subroutine new_sparse_matrix

  subroutine copy_sparse_matrix(inmat,outmat)
    !*FD copy a sparse matrix.
    !*FD Slap workspace allocation on the new
    !*FD matrix is *not* done.
    implicit none
    type(sparse_matrix_type) :: inmat  !*FD matrix to be copied
    type(sparse_matrix_type) :: outmat !*FD result matrix

    call new_sparse_matrix(inmat%order,inmat%nonzeros,outmat)
    outmat%row(:) = inmat%row(:)
    outmat%col(:) = inmat%col(:)
    outmat%val(:) = inmat%val(:)
    outmat%nonzeros = inmat%nonzeros
    outmat%symmetric = inmat%symmetric
  end subroutine copy_sparse_matrix

  subroutine grow_sparse_matrix(matrix)
    !*FD grow sparse matrix
    implicit none
    type(sparse_matrix_type) :: matrix !*FD matrix

    integer, dimension(:), pointer :: newrow,newcol
    real(kind=dp), dimension(:), pointer :: newval
    integer oldsize

    oldsize = size(matrix%val)
    
    allocate(newrow(chunksize+oldsize))
    allocate(newcol(chunksize+oldsize))
    allocate(newval(chunksize+oldsize))
    write(*,*)size(matrix%col), size(matrix%row), size(matrix%val), size(newcol), size(newrow), size(newval)
    newcol(1:oldsize) = matrix%col(:)
    newrow(1:oldsize) = matrix%row(:)
    newval(1:oldsize) = matrix%val(:)

    deallocate(matrix%col)
    deallocate(matrix%row)
    deallocate(matrix%val)

    matrix%col => newcol
    matrix%row => newrow
    matrix%val => newval

  end subroutine grow_sparse_matrix

  subroutine del_sparse_matrix(matrix)
    !*FD delete sparse matrix
    implicit none
    type(sparse_matrix_type) :: matrix !*FD matrix

    if (associated(matrix%col)) then
       deallocate(matrix%col)
       deallocate(matrix%row)
       deallocate(matrix%val)
    end if

  end subroutine del_sparse_matrix

  subroutine print_sparse(matrix, unit)
    !*FD print sparse matrix
    implicit none
    type(sparse_matrix_type) :: matrix !*FD matrix
    integer, intent(in) :: unit        !*FD unit to be printed to

    integer i
    do i = 1, matrix%nonzeros
       write(unit,*) matrix%col(i), matrix%row(i), matrix%val(i)
    end do
  end subroutine print_sparse

  subroutine sparse_matrix_vec_prod(matrix, vec, res)
    !*FD sparse matrix vector product
    implicit none
    type(sparse_matrix_type) :: matrix                !*FD matrix
    real(kind=dp), intent(in), dimension(:) :: vec    !*FD input vector
    real(kind=dp), intent(out), dimension(:) :: res   !*FD result vector

    integer i

    res = 0.
    do i=1,matrix%nonzeros
       res(matrix%col(i)) = res(matrix%col(i)) + vec(matrix%row(i))*matrix%val(i)
    end do
  end subroutine sparse_matrix_vec_prod

  subroutine sparse_insert_val(matrix, i, j, val)
    !*FD insert value into sparse matrix.  This is safe to call even if val=0
    implicit none
    type(sparse_matrix_type) :: matrix !*FD matrix
    integer, intent(in) :: i,j         !*FD column and row
    real(kind=dp), intent(in) :: val   !*FD value
    if (val /= 0.0 .and. i > 0 .and. j > 0 .and. i <= matrix%order .and. j <= matrix%order) then
        matrix%nonzeros =  matrix%nonzeros + 1
        matrix%row(matrix%nonzeros) = i
        matrix%col(matrix%nonzeros) = j
        matrix%val(matrix%nonzeros) = val

        if (matrix%nonzeros .eq. size(matrix%val)) then
            call grow_sparse_matrix(matrix)
        end if
    end if
  end subroutine sparse_insert_val

  subroutine sparse_clear(matrix)
    !*FD Clears the sparse matrix, without deallocating any of the
    !*FD previously used memory
    type(sparse_matrix_type) :: matrix
    
    matrix%nonzeros = 0
    !Clearing these shouldn't be strictly necessary, but SLAP barfs if we don't
    matrix%row = 0
    matrix%col = 0
    matrix%val = 0
  end subroutine

  function is_triad_format(matrix)
    type(sparse_matrix_type) :: matrix
    logical :: is_triad_format

    is_triad_format = .not. is_column_format(matrix) .and. .not.  is_row_format(matrix)
  end function

  function is_row_format(matrix)
    type(sparse_matrix_type) :: matrix
    logical :: is_row_format

    is_row_format = matrix%row(matrix%order + 1) == matrix%nonzeros + 1
  end function
!----------------------------------------------------------------------- 
      subroutine coicsr (n,nnz,job,a,ja,ia,iwk)
          use glimmer_global, only : dp
          implicit none
          integer,intent(in) :: n,nnz,job
          real(dp),dimension(:),intent(inout) :: a
          integer, dimension(:),intent(inout) :: ja,ia
          integer, dimension(:),intent(inout) :: iwk

          !Local
          real(kind=dp) :: t,tnext
          logical :: values
          integer :: i,j,k,init,ipos,inext,jnext

!------------------------------------------------------------------------
! IN-PLACE coo-csr conversion routine.
!------------------------------------------------------------------------
! this subroutine converts a matrix stored in coordinate format into 
! the csr format. The conversion is done in place in that the arrays 
! a,ja,ia of the result are overwritten onto the original arrays.
!------------------------------------------------------------------------
! on entry:
!--------- 
! n	= integer. row dimension of A.
! nnz	= integer. number of nonzero elements in A.
! job   = integer. Job indicator. when job=1, the real values in a are
!         filled. Otherwise a is not touched and the structure of the
!         array only (i.e. ja, ia)  is obtained.
! a	= real array of size nnz (number of nonzero elements in A)
!         containing the nonzero elements 
! ja	= integer array of length nnz containing the column positions
! 	  of the corresponding elements in a.
! ia	= integer array of length nnz containing the row positions
! 	  of the corresponding elements in a.
! iwk	= integer work array of length n+1 
! on return:
!----------
! a
! ja 
! ia	= contains the compressed sparse row data structure for the 
!         resulting matrix.
! Note: 
!-------
!         the entries of the output matrix are not sorted (the column
!         indices in each are not in increasing order) use coocsr
!         if you want them sorted.
!----------------------------------------------------------------------c
!  Coded by Y. Saad, Sep. 26 1989                                      c
!  Released under the LGPL   
!       
! Converted to F90 by JVJ -- 11/3/09
!----------------------------------------------------------------------c
!----------------------------------------------------------------------- 
      values = (job .eq. 1) 
! find pointer array for resulting matrix. 
      do i=1,n+1
         iwk(i) = 0
      end do
      do k=1,nnz
         i = ia(k)
         iwk(i+1) = iwk(i+1)+1
      end do 
!------------------------------------------------------------------------
      iwk(1) = 1 
      do i=2,n
         iwk(i) = iwk(i-1) + iwk(i)
      end do 
!
!     loop for a cycle in chasing process. 
!
      init = 1
      k = 0
 5    if (values) t = a(init)
      i = ia(init)
      j = ja(init)
      ia(init) = -1
!------------------------------------------------------------------------
 6    k = k+1
!     current row number is i.  determine  where to go. 
      ipos = iwk(i)
!     save the chased element. 
      if (values) tnext = a(ipos)
      inext = ia(ipos)
      jnext = ja(ipos)
!     then occupy its location.
      if (values) a(ipos)  = t
      ja(ipos) = j
!     update pointer information for next element to come in row i. 
      iwk(i) = ipos+1
!     determine  next element to be chased,
      if (ia(ipos) .lt. 0) goto 65
      t = tnext
      i = inext
      j = jnext 
      ia(ipos) = -1
      if (k .lt. nnz) goto 6
      goto 70
 65   init = init+1
      if (init .gt. nnz) goto 70
      if (ia(init) .lt. 0) goto 65
!     restart chasing --	
      goto 5
 70   do i=1,n 
         ia(i+1) = iwk(i)
      end do 
      ia(1) = 1
      return
  end subroutine
!----------------- end of coicsr ----------------------------------------


!----------------------------------------------------------------------- 
      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
          use glimmer_global, only : dp
          implicit none
          integer, intent(in) :: nrow,nnz
          real(dp),dimension(:),intent(in) :: a
          integer,dimension(:),intent(in) :: ir
          integer,dimension(:),intent(in) :: jc
          real(dp),dimension(:),intent(out) :: ao
          integer, dimension(:),intent(out) :: jao
          integer, dimension(:),intent(out) :: iao

          ! Local
          real(dp) :: x
          integer :: i,k,j,k0,iad
!----------------------j-------------------------------------------------
!  Coordinate to Compressed Sparse Row 
!  Written by Yousef Saad as part of SparseKit2
!  Released under the LGPL
! 
! Converted to F90 by JVJ -- 10/21/09
!----------------------------------------------------------------------- 
! converts a matrix that is stored in coordinate format
!  a, ir, jc into a row general sparse ao, jao, iao format.
!
! on entry:
!--------- 
! nrow	= dimension of the matrix 
! nnz	= number of nonzero elements in matrix
! a,
! ir, 
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
! 	  the elements, ir(k) = its row number and jc(k) = its column 
!	  number. The order of the elements is arbitrary. 
!
! on return:
!----------- 
! ir 	is destroyed
!
! ao, jao, iao = matrix in general sparse matrix format with ao 
! 	continung the real values, jao containing the column indices, 
!	and iao being the pointer to the beginning of the row, 
!	in arrays ao, jao.
!------------------------------------------------------------------------
      iao = 0
! determine row-lengths.
      do k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
      end do
! starting position of each row..
      k = 1
      do j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
      end do
! go through the structure  once more. Fill in output matrix.
      do k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
      end do
! shift back iao
      do j=nrow,1,-1
         iao(j+1) = iao(j)
      end do
      iao(1) = 1
      return
      end subroutine
!------------- end of coocsr ------------------------------------------- 
  function is_column_format(matrix)
    type(sparse_matrix_type) :: matrix
    logical :: is_column_format

    is_column_format = matrix%col(matrix%order + 1) == matrix%nonzeros + 1
  end function

  subroutine to_column_format(matrix)
    type(sparse_matrix_type) :: matrix
     
    if(is_triad_format(matrix)) then
        call ds2y(matrix%order, matrix%nonzeros, matrix%row, matrix%col, matrix%val, 0)
    end if
  end subroutine

  subroutine sort_column_format(matrix)
    !*FD Takes a column format matrix and sorts the row indices within each column
    !*FD This is not strictly needed in some compressed-column matrices
    !*FD (e.g. those used in SLAP), but it *is* necessary in some other libraries
    !*FD (e.g. UMFPACK).  For this reason, it is not done automatically in
    !*FD to_column_format.
    implicit none
    type(sparse_matrix_type) :: matrix
    integer :: i

    do i=1,matrix%order !Loop through each column index
      call sort(matrix%val, matrix%row, matrix%col(i), matrix%col(i+1)-1)
    end do
  end subroutine

  subroutine sort_row_format(matrix)
    !*FD Takes a row format matrix and sorts the column indices within each row
    !*FD This is not strictly needed in some compressed-row matrices
    !*FD (e.g. those used in SLAP), but it *is* necessary in some other libraries
    !*FD (e.g. PARDISO).  
    implicit none
    type(sparse_matrix_type),intent(inout) :: matrix
    integer :: i

    do i=1,matrix%order !Loop through each column index
      call sort(matrix%val, matrix%col, matrix%row(i), matrix%row(i+1)-1)
    end do
  end subroutine


  subroutine sort(values, indices, startindex, endindex)
    implicit none
    real(dp),dimension(:),intent(inout) :: values
    integer,dimension(:),intent(inout) :: indices
    integer, intent(in) :: startindex
    integer, intent(in) :: endindex
    integer :: currentindex
    real(dp) :: currentvalue
    integer :: i,j
    
    !Insertion Sort 
    do i=startindex+1,endindex
        currentindex = indices(i)
        currentvalue = values(i)

        j = i-1
        do while (j >= startindex .and. indices(j) > currentindex)
            indices(j+1) = indices(j)
            values(j+1) = values(j)
            j = j - 1
        end do
        indices(j+1) = currentindex
        values(j+1) = currentvalue
    end do
  end subroutine

end module glimmer_sparse_type
