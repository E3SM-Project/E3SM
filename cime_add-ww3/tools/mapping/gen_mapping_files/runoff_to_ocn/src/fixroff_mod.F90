MODULE fixroff_mod

   use shr_sys_mod
   use kind_mod
   use map_mod

   implicit none

   private ! default private

   public  :: fixroff_modifyMap   ! modify map wrt relocated points
   private :: fixroff_findActive  ! determine if & where to relocate runoff

   logical :: debug = .false. ! trigger debug statements

!===============================================================================
CONTAINS
!===============================================================================

SUBROUTINE fixroff_modifyMap(map, map2)

   implicit none

   !--- arguments ---
   type(sMatrix), intent(inout) :: map        ! map: original
   type(sMatrix), intent(inout) :: map2       ! map: modified wrt newpt

   !--- local ---
   integer(IN),allocatable :: newpt(:,:)          ! cell relocation info

   integer(IN)             :: i,j,i2,j2           ! 2d row,col indecies
   integer(IN)             :: n                   ! 2d vector index
   integer(IN),allocatable :: badpts (:)          ! 0 <=> no relocation necessary
   integer(IN),allocatable :: badptsf(:)          ! 0 <=> relocation necessary but not done
   integer(IN)             :: newptr(map%n_b)     !
   integer(IN)             :: cnt1,cnt2,cnt3,cnt4 ! counts various cases
   real(R8)   ,allocatable :: wcol(:)             ! cnt2 weight corrections
   real(R8)                :: dist,distMin        ! distance, min distance
   real(R8)                :: element             ! a matrix element

   !--- formats ---
   character(*),parameter :: subName = "(fixroff_modifyMap) "
   character(*),parameter :: F00 =   "('(fixroff_modifyMap) ',4a)"
   character(*),parameter :: F01 =   "('(fixroff_modifyMap) ',a,i8)"

!-------------------------------------------------------------------------------
! PURPOSE:
! o Returns a new version of "map", "map2", which corrects three kinds
!   of scrip mapping "errors" (errors wrt runoff mapping):
!
!   (case 1) cnt1 points: mask_a is non-zero, but frac_a = 0.0
!      issue:
!         the source cell is active, but 0% participates in the mapping;
!         these represent active source cells which
!         do not intersect any destination cells.
!      correction:
!         make frac_a = 1.0; add links for these points to
!         closest destination grid cell with s weight based
!         on simple area ratio. (NOTE:  may result in S>1.0)
!
!   (case 2) cnt2 points: mask_a is non-zero  and  0 < frac_a < 1
!      issue:
!         the source cell is active but only partially intersects
!         the active destination domain.
!      correction:
!         make frac_a = 1.0 and apply correction factor to s weight of 1/old_frac_a.
!         NOTE:  may result in S>1.0
!
!   (case 3) cnt3 points: mask_a is zero  and  0 < frac_a 1
!      issue:
!         the source cell is inactive, but the scrip matrix will map the data anyway
!         This indicates a script error that deserves attention.
!      correction:
!         no correction. (?but should set matrix element to 0?)
!
!   (case 4) cnt4 points: mask_b is zero, any runoff to this cell needs to be relocated
!      issue:
!         since scrip creates links on the basis of mask_a values,
!         (ie. assumes mask_b /= 0 everywhere), a link to this cell exists.
!         Since this is for runoff mapping, we cannot accept
!         such links to land cells on the destination grid.
!      correction:
!         Shift the flux to the nearest mask_b /= 0 point, as determined previously by
!         fixroff_findActive.  NOTE: this may result in a "doubling-up" of links
!         (non-unique links) between a given source and destination cell.
!-------------------------------------------------------------------------------

   !----------------------------------------------------------------------------
   ! preliminary data: determine where to relocate runoff, if necessary
   !----------------------------------------------------------------------------
   allocate(newpt(map%ni_b,map%nj_b))
   call fixroff_findActive(map, newpt)

   !----------------------------------------------------------------------------
   ! determine how many, and what types of, corrections are necessary
   !----------------------------------------------------------------------------
   cnt1 = 0 ! counts cells wrt case (1)
   cnt2 = 0 ! counts cells wrt case (2)
   cnt3 = 0 ! counts cells wrt case (3)
   cnt4 = 0 ! counts cells wrt case (4)

   allocate(badpts (map%n_b))
   allocate(badptsf(map%n_b))
   allocate(   wcol(map%n_a))

   wcol = 1.0_R8
   do n=1,map%n_a
      if (map%frac_a(n) == 0.0_R8 .and. map%mask_a(n) /= 0) then
         cnt1=cnt1+1
      endif
      if (map%frac_a(n) /= 0.0_R8 .and. abs(map%frac_a(n)-1.0_R8) > 1.e-20_R8 .and. &
      &   map%mask_a(n) /= 0) then
         wcol(n)       = 1.0_r8/map%frac_a(n) ! apply this correction to s(n)
         map%frac_a(n) = 1.0_R8
         cnt2=cnt2+1
      endif
      if (map%frac_a(n) /= 0.0_R8 .and. map%frac_a(n) /= 1.0_R8 .and. &
      &   map%mask_a(n) == 0) then
         cnt3=cnt3+1
      endif
      if (map%mask_b(n) == 0) then
         cnt4 = cnt4+1
      endif
   enddo
   write(6,F01) 'case(1): mask_a /= 0 AND  frac_a == 0.0     : ',cnt1
   write(6,F01) 'case(2): mask_a /= 0 AND (0 < frac_a < 1)   : ',cnt2
   write(6,F01) 'case(3): mask_a == 0 AND (0 < frac_a < 1)   : ',cnt3
   write(6,F01) 'case(4): mask_b == 0 AND is being mapped to : ',cnt4

   !----------------------------------------------------------------------------
   ! initialize up new map that is sufficiently larger (will add cnt1 links)
   !----------------------------------------------------------------------------

   !--- duplicate ---
   call map_dup(map,map2)

   !--- resize ---
   write(6,F01) "Original matrix size = ",map2%n_s
   deallocate( map2%s   )
   deallocate( map2%row )
   deallocate( map2%col )
   map2%n_s  = map%n_s + cnt1
   allocate(map2%s  (map2%n_s) )
   allocate(map2%row(map2%n_s) )
   allocate(map2%col(map2%n_s) )
   write(6,F01) "Modified matrix size = ",map2%n_s

   !--- copy initial data ---
   map2%s  (1:map%n_s) = map%s
   map2%row(1:map%n_s) = map%row
   map2%col(1:map%n_s) = map%col

   !----------------------------------------------------------------------------
   ! case (1) add new links for active src cells that are currently unmapped
   !----------------------------------------------------------------------------
   cnt1 = 0
   n    = map%n_s
   do j=1,map%n_a ! for each column
      if (map%frac_a(j) == 0.0_R8 .and. map%mask_a(j) /= 0) then
         !--- none of this src grid cell is mapped to the dest grid ---
         n    = n    + 1
         cnt1 = cnt1 + 1

         !--- find closest grid_b(i) cell to this grid_a(j) cell ---
         i = 1
         distMin = 1.e30_R8
         do i2 = 1,map%n_b
            dist = (map%yc_b(i2)-map%yc_a(j))**2 + (map%xc_b(i2)-map%xc_a(j))**2
            if (dist < distMin) then
               i = i2
               distMin = dist
            endif
         enddo

         !--- compute new matrix element for this i,j pair ---
         if (map%normal == 'fracarea') then
            element = (map%area_a(j)*map%frac_a(j))/(map%area_b(i)*map%frac_b(i))
         elseif (map%normal == 'destarea') then
            element =                map%area_a(j) / map%area_b(i)
         endif

         !--- add new element to matrix ---
         map2%frac_a(j) = 1.0_R8
         map2%row(n)    = i
         map2%col(n)    = j
         map2%S  (n)    = element

         if (debug) write(6,*) subName,'added element: ',n, &
         &  map2%xc_a(j),map2%yc_a(j),map2%xc_b(i),map2%yc_b(i),map2%s(n)
       endif
   enddo
   write(6,F01) 'case (1): new links added: ',cnt1

   !----------------------------------------------------------------------------
   ! case (2) adjust matrix weights for non-unity src fractions
   ! case (4) relocate links that map to mask_b=0 points.
   !    Do this by changing the element's row to the row specified by newpt index.
   !    Also must correct the matrix element value wrt src & dest cell areas
   !----------------------------------------------------------------------------
   badpts (:) = 0 ! 1 <=> cell needs relocation
   badptsf(:) = 0 ! 1 <=> cell relocation error
   newptr = reshape(newpt,(/map%n_b/))
   do n=1,map2%n_s

      i = map2%row(n) ! original/unmodified row
      j = map2%col(n) ! original col (this won't change)

      !--- case (2) adjust for non-unity src fractions ---
      map2%s(n)=map2%s(n)*wcol(j)

      !--- case (4) relocate mapping to masked out destination cell ---
      if (map2%mask_b(i) == 0) then
         badpts   (i) = 1 ! 1 <=> cell needs relocation
         badptsf  (i) = 1 ! 1 <=> cell relocation error
         i2=newptr(i)     ! new/modified row
         if (map2%mask_b(i2) /= 0) then ! OK to relocate row to here
            badptsf(i)=0 ! row was successfully relocated
            if ( i  < 1 .or. map2%n_b < i ) write(6,*) subName,'ERROR1 n,i1,i2:',n,i,i2
            if ( i2 < 1 .or. map2%n_b < i2) write(6,*) subName,'ERROR2:n,i1,i2:',n,i,i2
            if (map2%normal == 'fracarea') then
               if (map2%frac_b(i) == 0.0_R8 .or. map2%frac_b(i2) == 0.0_R8) then
                  write(6,*) subName,'ij3 ',n,i2,map2%frac_b(i),map2%frac_b(i2)
                  element = map2%s(n)*map2%area_b(i)/map2%area_b(i2)
               else
                  element = map2%s(n)*(map2%area_b(i )*map2%frac_b(i )) / &
                  &                   (map2%area_b(i2)*map2%frac_b(i2))
               endif
            elseif (map2%normal == 'destarea') then
                element = map2%s(n)*map2%area_b(i)/map2%area_b(i2)
            endif
            map2%s  (n) = element  ! corrected matrix element
            map2%row(n) = i2       ! new row
        endif
      endif
   enddo

   open(20,file="fixroff_modifyMap_info.txt",form="formatted")
   do n=1,map2%n_b
      if (badptsf(n) /= 0) write(6,F01) 'CAUTION: badpt /= 0 ',n
      write(20,*) map2%mask_b(n),badpts(n),badptsf(n),newptr(n)
   enddo
   close(20)

END SUBROUTINE fixroff_modifyMap

!===============================================================================
!===============================================================================

SUBROUTINE fixroff_findActive(map, newpt)

   implicit none

   !--- arguments ---
   type(sMatrix),intent(in)    :: map        ! map in question
   integer(IN)  ,intent(inout) :: newpt(:,:) ! indicates where to relocate runoff

   !--- local ---
   integer(IN),allocatable :: newpt_old(:,:) ! temporary work array

   integer(IN) :: ni,nj       ! size of grid_b
   integer(IN) :: i,j         ! loop index into 2D array
   integer(IN) :: n           ! loop index into 1D vector
   integer(IN) :: i2,j2       ! index of neighbor cells

   integer(IN) :: nFound      ! # cells whose nearest active neighbor was found
   integer(IN) :: nNotFound   ! # cells whose nearest active neighbor is unknown
   logical     :: found       ! an active cell has been found
   logical     :: done        ! T => active cells have been found for all cells

   !--- formats ---
   character(*),parameter :: subName = "(fixroff_findActive) "

!-------------------------------------------------------------------------------
! PURPOSE:
! o  For each destination grid cell, find the nearest active (unmasked) grid cell.
!    For destination grid points which are  active  (ie. mask_b(n) /= 0)
!    the nearest active (unmasked) grid cell is itself.
!    For destination grid points which are inactive (ie. mask_b(n) == 0)
!    the nearest active cell must be searched for.
!
! NOTES:
! o  The sMatrix "map", which maps from grid_a -> grid_b,
!    where the size of grid_b, shaped as a 2d array,  is (ni_b,nj_b), and
!    where the size of grid_b, shaped as a 1d vector, is (n_b)...
!
! o  The "newPt(:,:)" data array is a 2d array dimensioned newPt(ni_b,nj_b)
!    assign a 1d index n, between 1->n_b, for each point of grid_b.
!    If mask_b(n) /= 0, the index is n
!    If mask_b(n) == 0, the index is that of the nearest neigbor m
!                       for which mask_b(m)==1
! HISTORY:
!    2007-01-05 - B. Kauffman, code cleanup
!    2002-06-12 - M. Hecht, original version
!-------------------------------------------------------------------------------

   ni = map%ni_b
   nj = map%nj_b

   !--- nearest neighbors are known for all active cells ---
   do j=1,nj
   do i=1,ni
      n = (j-1)*ni + i ! corresponding 1d vector index
      if (map%mask_b(n) /= 0) then
         newpt(i,j) = n  ! nearest neighbor is itself
      else
         newpt(i,j) = 0  ! nearest neighbor is unknown
      endif
   enddo
   enddo

   !--- search for nearest neighbors for all inactive cells ---
   allocate(newpt_old(ni,nj))
   done = .false.
   do while ( .not. done )
      newpt_old(:,:) = newpt(:,:)

      nNotFound = 0
      nFound    = 0

      do j=1,nj
      do i=1,ni
         if (newpt(i,j) == 0) then ! nearest active neighbor not found
            found = .false.
 loop:      do n=1,4 ! check immediate neighbors
               i2=i
               j2=j
               if (n == 1) i2=mod(i     ,ni)+1 ! right neighbor
               if (n == 2) i2=mod(i-2+ni,ni)+1 ! left  neighbor
               if (n == 3) j2=min(j+1   ,nj)   ! neighbor above
               if (n == 4) j2=max(j-1   , 1)   ! neighbor below
               if (newpt_old(i2,j2) /= 0) then
                  found = .true.
                  newpt(i,j)=newpt_old(i2,j2)
                  exit loop
               endif
            enddo loop
            if (      found) nFound    = nFound    + 1
            if (.not. found) nNotFound = nNotFound + 1
         endif
      enddo
      enddo
      if (nNotFound == 0) done = .true.
      write(6,'(2a,2i8)') subName,"num of cells found, unknown = ",nFound,nNotFound
   enddo
   deallocate(newpt_old)

!  !--- write debugging info --- DELETE THIS VERSION
!  open(21,form="formatted")
!  do j=1,nj
!  do i=1,ni
!     n = (j-1)*map%ni_b+i
!     write(21,*)   map%mask_b(n),newpt(i,j)
!  enddo
!  enddo
!  close(21)
!  !---------------------------- DELETE THIS VERSION

   !--- write debugging info ---
   open(10,file="fixroff_findActive_info.txt",form="formatted")
   write(10,*) "cell index, mask value, nearest active neighbor index"
   do j=1,nj
   do i=1,ni
      n = (j-1)*ni+i
      write(10,*) n,map%mask_b(n),newpt(i,j)
   enddo
   enddo
   close(10)

END SUBROUTINE fixroff_findActive

!===============================================================================
!===============================================================================

END MODULE fixroff_mod
