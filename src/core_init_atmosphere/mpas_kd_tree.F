module mpas_kd_tree

   !***********************************************************************
   !
   !  module mpas_kd_tree
   !
   !> \brief   MPAS KD-Tree module
   !> \author  Miles A. Curry
   !> \date    01/28/20
   !> A KD-Tree implementation to create and search perfectly balanced
   !> KD-Trees.
   !>
   !> Use `mpas_kd_type` dervied type to construct points for mpas_kd_construct:
   !>
   !> real (kind=RKIND), dimension(:,:), allocatable :: array
   !> type (mpas_kd_type), pointer :: tree => null()
   !> type (mpas_kd_type), dimension(:), pointer :: points => null()
   !>
   !> allocate(array(k,n)) ! K dims and n points
   !> allocate(points(n))
   !> array(:,:) = (/.../)  ! Fill array with values
   !>
   !> do i = 1, n
   !>    allocate(points(i) % point(k))    ! Allocate point with k dimensions
   !>    points(i) % point(:) = array(:,i)
   !>    points(i) % id = i                ! Or a value of your choice
   !> enddo
   !>
   !> tree => mpas_kd_construct(points, k)
   !>
   !> call mpas_kd_free(tree)
   !> deallocate(points)
   !> deallocate(array)
   !>
   !
   !-----------------------------------------------------------------------
   use mpas_kind_types, only : RKIND

   implicit none

   private

   public :: mpas_kd_type

   ! Public Subroutines
   public :: mpas_kd_construct
   public :: mpas_kd_search
   public :: mpas_kd_free

   type mpas_kd_type
      type (mpas_kd_type), pointer :: left => null()
      type (mpas_kd_type), pointer :: right => null()

      integer :: split_dim
      real (kind=RKIND), dimension(:), pointer :: point => null()

      integer :: id
   end type mpas_kd_type

   contains

   !***********************************************************************
   !
   !  recursive routine mpas_kd_construct_internal
   !
   !> \brief   Create a KD-Tree from a set of k-Dimensional points
   !> \author  Miles A. Curry
   !> \date    01/28/20
   !> \details
   !> Private, recursive function to construct a KD-Tree from an array
   !> of mpas_kd_type, points, and return the root of the tree.
   !>
   !> ndims should be the dimensioned of each individual point found
   !> in points and npoints should be the number of points. dim represents
   !> the current split dimensioned and is used internally. Upon calling
   !> this function, dim should always be set to 0.
   !
   !-----------------------------------------------------------------------
   recursive function mpas_kd_construct_internal(points, ndims, npoints, dim) result(tree)

      implicit none

      ! Input Variables
      type (mpas_kd_type), dimension(:), target :: points
      integer, intent(in) :: ndims
      integer, value :: npoints
      integer, value :: dim

      ! Return Value
      type (mpas_kd_type), pointer :: tree

      ! Local Variables
      integer :: median

      if (npoints < 1) then
         tree => null()
         return
      endif

      ! Sort the points at the split dimension
      dim = mod(dim, ndims) + 1
      call quickSort(points, dim, 1, npoints, ndims)

      median = (1 + npoints) / 2

      points(median) % split_dim = dim
      tree => points(median)

      ! Build the right and left sub-trees but do not include the median
      ! point (the root of the current tree)
      if (npoints /= 1) then
          points(median) % left => mpas_kd_construct_internal(points(1:median-1), ndims, median - 1, points(median) % split_dim)
          points(median) % right => mpas_kd_construct_internal(points(median+1:npoints), ndims, npoints - median, &
                                                                                             points(median) % split_dim)
      endif

   end function mpas_kd_construct_internal


   !***********************************************************************
   !
   !  routine mpas_kd_construct
   !
   !> \brief   Construct a balanced KD-Tree
   !> \author  Miles A. Curry
   !> \date    01/28/20
   !> \details
   !> Create and return a perfectly balanced KD-Tree from an array of
   !> mpas_kd_type, points. The point member of every element of the points
   !> array should be allocated and set to the points desired to be in the
   !> KD-Tree and ndims should be the dimensions of the points.
   !>
   !> Upon error, the returned tree will be unassociated.
   !
   !-----------------------------------------------------------------------
   function mpas_kd_construct(points, ndims) result(tree)

      implicit none

      ! Input Varaibles
      type (mpas_kd_type), dimension(:) :: points
      integer, intent(in) :: ndims

      ! Return Value
      type (mpas_kd_type), pointer :: tree

      ! Local Varaibles
      integer :: npoints

      npoints = size(points)

      if (npoints < 1) then
         tree => null()
         return
      endif

      tree => mpas_kd_construct_internal(points(:), ndims, npoints, 0)

   end function mpas_kd_construct

   !***********************************************************************
   !
   !  routine break_tie
   !
   !> \brief   Break a tie for two n-dim points
   !> \author  Miles A. Curry
   !> \date    07/07/20
   !> \details
   !> Compare 1..n dimensions of p1 and p2 and return -1 if p1(i) is less than
   !> p2(i) and return 1 if p1(i) is greater than p2(i). If p1(i) and p2(i) are
   !> equal, then the same comparison will be done on p1(i+1) and p2(i+1) until
   !> p1(n) and p2(n). If p1(:) and p2(:) are equal across all n, then 0 will
   !> be returned.
   !
   !-----------------------------------------------------------------------
   function break_tie(p1, p2) result(tie)

      implicit none

      ! Input Variables
      type (mpas_kd_type), intent(in) :: p1
      type (mpas_kd_type), intent(in) :: p2
      integer :: tie

      integer :: i

      tie = 0
      do i = 1, size(p1 % point(:))
          if (p1 % point(i) < p2 % point(i)) then
              tie = -1
              return
          else if (p1 % point(i) > p2 % point(i)) then
              tie = 1
              return
          endif
      enddo

   end function break_tie


   !***********************************************************************
   !
   !  recursive routine mpas_kd_search_internal
   !
   !> \brief   Recursively search the KD-Tree for query
   !> \author  Miles A. Curry
   !> \date    01/28/20
   !> \details
   !> Private, recursive function to search kdtree for query. Upon succes
   !> res will point to the nearest neighbor to query and distance will hold
   !> the squared distance between query and res.
   !>
   !> Distance is calculated and compared as squared distance to increase
   !> efficiency.
   !
   !-----------------------------------------------------------------------
   recursive subroutine mpas_kd_search_internal(kdtree, query, res, distance)

      implicit none

      ! Input Variables
      type (mpas_kd_type), pointer, intent(in) :: kdtree
      real (kind=RKIND), dimension(:), intent(in) :: query
      type (mpas_kd_type), pointer, intent(inout) :: res
      real (kind=RKIND), intent(inout) :: distance

      ! Local Values
      real (kind=RKIND) :: current_distance

      current_distance = sum((kdtree % point(:) - query(:))**2)
      if (current_distance < distance) then
         distance = current_distance
         res => kdtree
      else if (current_distance == distance) then
          !
          ! Consistently break a tie if a query is equidistant from two points
          !
          if (associated(res)) then
              if (break_tie(res, kdtree) == 1) then
                 res => kdtree
              endif
          endif
      endif

      !
      ! To find the nearest neighbor, first serach the tree in a similar manner
      ! as a single dimensioned BST, by comparing points on the current split
      ! dimension.
      !
      ! If the distance between the current node and the query is less then the
      ! minimum distance found within the subtree we just searched, then the nearest
      ! neighbor might be in the opposite subtree, so search it.
      !

      if (query(kdtree % split_dim) > kdtree % point(kdtree % split_dim)) then
         if (associated(kdtree % right)) then ! Search right
            call mpas_kd_search_internal(kdtree % right, query, res, distance)
         endif
         if ((kdtree % point(kdtree % split_dim) - query(kdtree % split_dim))**2 <= distance .and. associated(kdtree % left)) then
            call mpas_kd_search_internal(kdtree % left, query, res, distance) ! Check the other subtree
         endif
      else if (query(kdtree % split_dim) < kdtree % point(kdtree % split_dim)) then
         if (associated(kdtree % left)) then ! Search left
            call mpas_kd_search_internal(kdtree % left, query, res, distance)
         endif
         if ((kdtree % point(kdtree % split_dim) - query(kdtree % split_dim))**2 <= distance .and. associated(kdtree % right)) then
            call mpas_kd_search_internal(kdtree % right, query, res, distance) ! Check the other subtree
         endif
      else ! Nearest point could be in either left or right subtree, so search both
         if (associated(kdtree % right)) then
            call mpas_kd_search_internal(kdtree % right, query, res, distance)
         endif
         if (associated(kdtree % left)) then
            call mpas_kd_search_internal(kdtree % left, query, res, distance)
         endif
      endif

   end subroutine mpas_kd_search_internal

   !***********************************************************************
   !
   !  routine mpas_kd_search
   !
   !> \brief   Find the nearest point in a KD-Tree to a query
   !> \author  Miles A. Curry
   !> \date    01/28/20
   !> \details
   !> Search kdtree and returned the nearest point to query into the
   !> res argument. Optionally, if distance is present, returned the
   !> squared distance between query and res.
   !>
   !> If the dimension of query does not match the dimensions of points
   !> within kdtree, then res will be returned as unassociated. Likewise,
   !> if kdtree is empty/unassociated, res will be returned as unassociated.
   !
   !-----------------------------------------------------------------------
   subroutine mpas_kd_search(kdtree, query, res, distance)

      implicit none
      type (mpas_kd_type), pointer, intent(in) :: kdtree
      real (kind=RKIND), dimension(:), intent(in) :: query
      type (mpas_kd_type), pointer, intent(inout) :: res
      real (kind=RKIND), intent(out), optional :: distance

      real (kind=RKIND) :: dis

      if (.not. associated(kdtree)) then
         res => null()
         return
      end if

      if (size(kdtree % point) /= size(query)) then
         res => null()
         return
      endif

      dis = huge(dis)
      call mpas_kd_search_internal(kdtree, query, res, dis)

      if (present(distance)) then
         distance = dis
      endif

   end subroutine mpas_kd_search

   !***********************************************************************
   !
   !  routine mpas_kd_free
   !
   !> \brief   Free all nodes within a tree.
   !> \author  Miles A. Curry
   !> \date    01/28/20
   !> \details
   !> Deallocate and nullify all point nodes of kdtree and nullify the
   !> left and right pointers.
   !>
   !> After calling this function, the array of mpas_kd_type that was used
   !> to construct kdtree will still be allocated and will need to be
   !> deallocated separate from this routine.
   !
   !-----------------------------------------------------------------------
   recursive subroutine mpas_kd_free(kdtree)

      implicit none
      type (mpas_kd_type), pointer :: kdtree

      if (.not. associated(kdtree)) then
         return
      endif

      if (associated(kdtree % left)) then
         call mpas_kd_free(kdtree % left)
      endif

      if (associated(kdtree % right)) then
         call mpas_kd_free(kdtree % right)
      endif

      deallocate(kdtree % point)
      nullify(kdtree % left)
      nullify(kdtree % right)
      nullify(kdtree)

   end subroutine mpas_kd_free


   !***********************************************************************
   !
   !  routine mpas_kd_quicksort
   !
   !> \brief   Sort an array along a dimension
   !> \author  Miles A. Curry
   !> \date    01/28/20
   !> \details
   !> Sort points starting from arrayStart, to arrayEnd along the given dimension
   !> `dim`. If two points are swapped, the entire K-Coordinate point are swapped.
   !
   !-----------------------------------------------------------------------
   recursive subroutine quickSort(array, dim, arrayStart, arrayEnd, ndims)

      implicit none

      ! Input Variables
      type (mpas_kd_type), dimension(:) :: array
      integer, intent(in), value :: dim
      integer, intent(in), value :: arrayStart, arrayEnd
      integer, intent(in) :: ndims

      ! Local Variables
      type (mpas_kd_type) :: temp
      real (kind=RKIND), dimension(ndims) :: pivot_value

      integer :: l, r, pivot, s

      if ((arrayEnd - arrayStart) < 1) then
         return
      endif

      ! Create the left, right, and start pointers
      l = arrayStart
      r = arrayEnd - 1
      s = l

      pivot = (l+r)/2
      pivot_value = array(pivot) % point

      ! Move the pivot to the far right
      temp = array(pivot)
      array(pivot) = array(arrayEnd)
      array(arrayEnd) = temp

      do while (.true.)
         ! Advance the left pointer until it is a value less then our pivot_value(dim)
         do while (.true.)
            if (array(l) % point(dim) < pivot_value(dim)) then
               l = l + 1
            else
               exit
            endif
         enddo

         ! Advance the right pointer until it is a value more then our pivot_value(dim)
         do while (.true.)
            if (r <= 0) then
               exit
            endif

            if(array(r) % point(dim) >= pivot_value(dim)) then
               r = r - 1
            else
               exit
            endif
         enddo

         if (l >= r) then
            exit
         else ! Swap elements about the pivot
            temp = array(l)
            array(l) = array(r)
            array(r) = temp
         endif
      enddo

      ! Move the pivot to l ended up
      temp = array(l)
      array(l) = array(arrayEnd)
      array(arrayEnd) = temp

      ! Quick Sort on the lower partition
      call quickSort(array(:), dim, s, l-1, ndims)

      ! Quick sort on the upper partition
      call quickSort(array(:), dim, l+1, arrayEnd, ndims)

   end subroutine quicksort

end module mpas_kd_tree
